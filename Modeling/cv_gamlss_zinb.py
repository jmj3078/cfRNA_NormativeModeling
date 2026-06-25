#!/usr/bin/env python
"""Stratified 5-fold CV benchmark using NB-ZINB (Zero-Inflated NBI) via rpy2.

ZINB is a zero-inflated mixture model:
    P(Y = 0)   = nu(X) + (1-nu(X)) * f_NBI(0 | mu(X), sigma(X))
    P(Y = k>0) = (1-nu(X)) * f_NBI(k | mu(X), sigma(X))

All three sub-parameters (mu, sigma, nu) are modelled as linear functions of
the 10 batch-bias covariates.

Three complementary z-scores are computed per gene per sample:

    z_binary     : Bernoulli RQR using the combined P(Y=0) from ZINB
                   (all N_hc samples)
    z_count_cond : Conditional NBI RQR given Y > 0
                   (only detected samples; NA for zeros)
    z_full       : Full ZINB CDF residual (all samples)

Usage:
    python cv_zanbi.py --limit 3
    python cv_zanbi.py --binary-z-threshold 0.10
    python cv_zanbi.py --nbi-stats CV_Results/cv_gamlss_stats.csv
"""

import argparse
import csv
import pickle
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import norm, skew, kurtosis
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects.numpy2ri as rpyn

warnings.filterwarnings("ignore", category=UserWarning, module="rpy2")
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

# ── Paths / Constants (root config.py 단일 소스) ───────────────────
import sys
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
import config

BASE_DIR  = config.MODELING_DIR
DATA_DIR  = config.DATA_DIR
H5AD_PATH = config.H5AD_PATH
SAVE_DIR  = config.CV_RESULTS_DIR
R_HELPER  = config.R_HELPER

BIAS_COLUMNS   = config.BIAS_COLUMNS
STRATIFY_COL   = config.MODELING_PARAMS["stratify_col"]
DET_RATE_MIN   = config.MODELING_PARAMS["det_rate_min"]   # 0.01로 통일
MEAN_COUNT_MIN = config.MODELING_PARAMS["mean_count_min"]
N_SPLITS       = config.MODELING_PARAMS["n_splits"]

META_FIELDS = [
    "gene", "n_hc", "det_rate_hc", "mean_count_hc",
    # binary z-score (detectability, all samples)
    "w1_binary",
    "mean_binary_z", "std_binary_z", "skew_binary_z", "kurt_binary_z", "n_binary",
    # count conditional z-score (detected samples only)
    "w1_count_cond",
    "mean_count_cond_z", "std_count_cond_z",
    "skew_count_cond_z", "kurt_count_cond_z", "n_count_cond",
    # full ZINB z-score (all samples)
    "w1_full",
    "mean_full_z", "std_full_z", "skew_full_z", "kurt_full_z", "n_full",
    # model parameters
    "mu_mean", "sigma_mean", "nu_mean",
    # outlier refinement
    "n_removed",
    # primary z selection
    "primary_z_type",
    "fold_success_rate",
    "time_s",
    # ── Failure / anomaly flags ──────────────────────────────────
    "flag_fit_failure",     # fold_success_rate < 1.0
    "flag_nu_explosion",    # nu_mean > 0.95 or nu_mean < 0.01
    "flag_sigma_explosion", # sigma_mean > 50
    "flag_mu_explosion",    # mu_mean > 1e5 or mu_mean < 1e-3
    "flag_z_bias",          # |mean_z| > 0.3 on binary or full z-score
    "flag_z_unstable",      # std_full_z > 1.5 or std_full_z < 0.5
    "flag_high_removal",    # n_removed > 5% of total training data
    "flag_nonnormal",       # w1_full > 0.25
    "any_flag",             # OR of all above flags
]

# Failure thresholds (all adjustable via CLI where noted)
_THR = dict(
    fold_success_rate = 1.0,   # < this → fit_failure
    nu_hi             = 0.95,  # nu_mean > this → nu_explosion
    nu_lo             = 0.01,  # nu_mean < this → nu_explosion
    sigma_hi          = 50.0,  # sigma_mean > this → sigma_explosion
    mu_hi             = 1e5,   # mu_mean > this → mu_explosion
    mu_lo             = 1e-3,  # mu_mean < this → mu_explosion
    z_bias            = 0.3,   # |mean_z| > this → z_bias
    std_hi            = 1.5,   # std_z > this → z_unstable
    std_lo            = 0.5,   # std_z < this → z_unstable
    removal_frac      = 0.05,  # n_removed/(n_hc*0.8) > this → high_removal
    w1                = 0.25,  # w1_full > this → nonnormal
)


# ── Failure detection ─────────────────────────────────────────────

def detect_flags(row: dict, n_hc: int, n_folds: int) -> dict:
    """Compute anomaly flags for one gene's result row.

    Returns a dict of {flag_name: int (0/1)} plus "any_flag".
    """
    T = _THR
    flags = {}

    flags["flag_fit_failure"]     = int(row["fold_success_rate"] < T["fold_success_rate"])
    nu = row["nu_mean"]
    flags["flag_nu_explosion"]    = int(
        (np.isfinite(nu) and (nu > T["nu_hi"] or nu < T["nu_lo"]))
        or not np.isfinite(nu)
    )
    sigma = row["sigma_mean"]
    flags["flag_sigma_explosion"] = int(
        (np.isfinite(sigma) and sigma > T["sigma_hi"])
        or not np.isfinite(sigma)
    )
    mu = row["mu_mean"]
    flags["flag_mu_explosion"]    = int(
        (np.isfinite(mu) and (mu > T["mu_hi"] or mu < T["mu_lo"]))
        or not np.isfinite(mu)
    )

    mean_bin  = row["mean_binary_z"]
    mean_full = row["mean_full_z"]
    flags["flag_z_bias"] = int(
        (np.isfinite(mean_bin)  and abs(mean_bin)  > T["z_bias"])
        or (np.isfinite(mean_full) and abs(mean_full) > T["z_bias"])
    )

    std_full = row["std_full_z"]
    flags["flag_z_unstable"] = int(
        (np.isfinite(std_full) and (std_full > T["std_hi"] or std_full < T["std_lo"]))
        or not np.isfinite(std_full)
    )

    # Expected training samples per fold ≈ n_hc * (1 - 1/n_folds)
    n_tr_expected = n_hc * (1 - 1 / n_folds)
    flags["flag_high_removal"] = int(
        row["n_removed"] > n_tr_expected * T["removal_frac"]
    )

    w1 = row["w1_full"]
    flags["flag_nonnormal"] = int(
        (np.isfinite(w1) and w1 > T["w1"])
        or not np.isfinite(w1)
    )

    flags["any_flag"] = int(any(v for k, v in flags.items() if k != "any_flag"))
    return flags


# ── R initialisation ───────────────────────────────────────────────

def init_r():
    ro.r(f'source("{R_HELPER}")')
    return ro.globalenv["fit_zinb_gene"]


# ── Data loading ───────────────────────────────────────────────────

def load_hc_data():
    adata = sc.read_h5ad(H5AD_PATH)
    adata = adata[adata.obs["QC_Passed"] == True]
    adata = adata[adata.obs["Phenotype_Processed"].notna()]
    adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
    adata = adata[adata.obs["broad_protocol_category"] != "Exome-based (EB)"]  # WTS only
    is_hc = (adata.obs["Phenotype_Processed"].astype(str) == "Healthy Control").values

    batch_raw  = adata.obs[STRATIFY_COL].astype(object)
    strata_all = batch_raw.fillna("Unknown").astype(str).values

    X_raw = adata.obs[BIAS_COLUMNS].values.astype(np.float64)
    scaler = StandardScaler()
    X_hc_scaled = scaler.fit_transform(X_raw[is_hc])

    Y_raw = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
    Y_raw = np.round(Y_raw).astype(np.float64)

    is_pc         = (adata.var["GeneType"] == "protein_coding").values
    pc_gene_names = adata.var_names[is_pc].tolist()
    pc_indices    = np.where(is_pc)[0]

    strata_hc = strata_all[is_hc]
    return X_hc_scaled, Y_raw, is_hc, strata_hc, pc_gene_names, pc_indices


def select_candidate_genes(Y_raw, is_hc, pc_gene_names, pc_indices,
                            det_rate_min):
    Y_hc  = Y_raw[is_hc][:, pc_indices]
    det_r = (Y_hc > 0).mean(axis=0)
    cand  = det_r >= det_rate_min
    return np.array(pc_gene_names)[cand].tolist(), pc_indices[cand]


# ── CV folds ───────────────────────────────────────────────────────

def make_stratified_folds(strata, n_splits=N_SPLITS):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    return list(skf.split(np.zeros(len(strata)), strata))


# ── Evaluation helpers ─────────────────────────────────────────────

def wasserstein1_normal(z_arr: np.ndarray) -> float:
    """1st Wasserstein distance between empirical z distribution and N(0,1).
    N-robust: converges without inflating with n. Threshold: > 0.25 = poor calibration.
    """
    z_valid = z_arr[np.isfinite(z_arr)]
    n = len(z_valid)
    if n < 8:
        return np.nan
    z_ref = norm.ppf(np.linspace(1 / (2 * n), 1 - 1 / (2 * n), n))
    return float(np.mean(np.abs(np.sort(z_valid) - z_ref)))


def zscore_stats(z_arr):
    z_valid = z_arr[np.isfinite(z_arr)]
    n = len(z_valid)
    if n < 2:
        return np.nan, np.nan, np.nan, np.nan, n
    return (float(z_valid.mean()), float(z_valid.std()),
            float(skew(z_valid)), float(kurtosis(z_valid)), n)


# ── rpy2 data conversion helpers ───────────────────────────────────

def _to_r_matrix(arr: np.ndarray, col_names: list):
    with localconverter(ro.default_converter + rpyn.converter):
        r_mat = ro.conversion.py2rpy(np.ascontiguousarray(arr, dtype=np.float64))
    return ro.r["matrix"](
        r_mat,
        nrow=arr.shape[0], ncol=arr.shape[1],
        dimnames=ro.r["list"](ro.NULL, ro.StrVector(col_names)),
    )


def _to_r_vec(arr: np.ndarray):
    with localconverter(ro.default_converter + rpyn.converter):
        return ro.conversion.py2rpy(np.ascontiguousarray(arr, dtype=np.float64))


# ── Primary z-score decision ───────────────────────────────────────

def _build_primary_lookup(nbi_stats_path):
    """Returns dict {gene: "binary"} for genes that should use Binary_Z
    based on prior NBI modeling quality."""
    if nbi_stats_path is None:
        return {}
    df = pd.read_csv(nbi_stats_path)
    bad = df[(df["fold_success_rate"] < 1.0) | (df["ad_pass"] == 0)]
    return {g: "binary" for g in bad["gene"]}


def decide_primary(det_rate, gene, binary_threshold, nbi_lookup):
    """Returns 'binary', 'count_cond', or 'full'."""
    use_binary = False
    if binary_threshold is not None and det_rate < binary_threshold:
        use_binary = True
    if gene in nbi_lookup:
        use_binary = True
    if use_binary:
        return "binary"
    if binary_threshold is not None or nbi_lookup:
        return "count_cond"
    return "full"


# ── Per-gene CV ────────────────────────────────────────────────────

def eval_gene_cv(y_hc, X_hc_scaled, folds, r_fit_fn, col_names, args):
    n = len(y_hc)
    z_binary     = np.full(n, np.nan)
    z_count_cond = np.full(n, np.nan)
    z_full       = np.full(n, np.nan)
    mu_all       = np.full(n, np.nan)
    sigma_all    = np.full(n, np.nan)
    nu_all       = np.full(n, np.nan)
    n_success    = 0
    n_removed    = 0

    for fold_idx, (tr_idx, te_idx) in enumerate(folds):
        seed_r          = ro.IntVector([42 + fold_idx])
        n_cyc_r         = ro.IntVector([50])
        outlier_z       = ro.FloatVector([args.outlier_z])
        max_iter        = ro.IntVector([args.max_iter])
        max_remove_frac = ro.FloatVector([args.max_remove_frac])
        nu_type_r       = ro.StrVector([args.nu_formula])
        lambda_sigma    = ro.FloatVector([args.lambda_sigma])

        res = r_fit_fn(
            _to_r_vec(y_hc[tr_idx]),
            _to_r_vec(y_hc[te_idx]),
            _to_r_matrix(X_hc_scaled[tr_idx], col_names),
            _to_r_matrix(X_hc_scaled[te_idx], col_names),
            seed_r, n_cyc_r, outlier_z, max_iter, max_remove_frac, nu_type_r, lambda_sigma,
        )

        if res.rx2("success")[0]:
            z_binary[te_idx]     = np.array(res.rx2("z_binary"))
            z_count_cond[te_idx] = np.array(res.rx2("z_count_cond"))
            z_full[te_idx]       = np.array(res.rx2("z_full"))
            mu_all[te_idx]       = np.array(res.rx2("mu_test"))
            sigma_all[te_idx]    = np.array(res.rx2("sigma_test"))
            nu_all[te_idx]       = np.array(res.rx2("nu_test"))
            n_removed           += int(res.rx2("n_removed")[0])
            n_success           += 1

    return (z_binary, z_count_cond, z_full,
            mu_all, sigma_all, nu_all,
            n_success / len(folds), n_removed)


# ── Resume helper ──────────────────────────────────────────────────

def _load_done_genes(meta_path):
    done = set()
    if meta_path.exists():
        with open(meta_path, newline="") as f:
            for row in csv.DictReader(f):
                done.add(row["gene"])
    return done


# ── Main ───────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--det-rate-min",        type=float, default=DET_RATE_MIN)
    parser.add_argument("--n-folds",              type=int,   default=N_SPLITS)
    parser.add_argument("--binary-z-threshold",   type=float, default=None,
                        metavar="T",
                        help="genes with det_rate < T use Binary_Z as primary (e.g. 0.10)")
    parser.add_argument("--nbi-stats",            type=str,   default=None,
                        metavar="PATH",
                        help="path to cv_gamlss_stats.csv; genes with poor NBI fit "
                             "get Binary_Z as primary")
    parser.add_argument("--nu-formula",       choices=["intercept", "full"],
                        default="intercept",
                        help="nu sub-model formula: 'intercept' (nu~1, default, prevents "
                             "overfitting) or 'full' (nu~all covariates)")
    parser.add_argument("--outlier-z",        type=float, default=5.0,
                        help="remove training samples with |z_train| > this (default 5.0)")
    parser.add_argument("--max-iter",         type=int,   default=2,
                        help="max outlier-removal iterations per fold (default 2)")
    parser.add_argument("--max-remove-frac",   type=float, default=0.10,
                        help="max fraction of training samples removable per iteration (default 0.10)")
    parser.add_argument("--lambda-sigma",   type=float, default=0.05,
                        help="L2 ridge penalty on sigma submodel coefficients (default 0.05; 0=disabled)")
    parser.add_argument("--no-resume",  action="store_true")
    parser.add_argument("--limit",      type=int,   default=None,
                        help="process only first N genes (for testing)")
    args = parser.parse_args()

    SAVE_DIR.mkdir(parents=True, exist_ok=True)
    meta_path    = SAVE_DIR / "cv_zinb_stats.csv"
    zscores_path = SAVE_DIR / "cv_zinb_zscores.pkl"

    print("Initialising R / ZINB...")
    r_fit_fn = init_r()

    nbi_lookup = _build_primary_lookup(args.nbi_stats)
    if args.nbi_stats:
        print(f"NBI stats loaded: {len(nbi_lookup)} genes flagged for Binary_Z")

    print("Loading HC data...")
    X_hc_scaled, Y_raw, is_hc, strata_hc, pc_gene_names, pc_indices = load_hc_data()
    gene_names, gene_indices = select_candidate_genes(
        Y_raw, is_hc, pc_gene_names, pc_indices,
        args.det_rate_min,
    )
    Y_hc = Y_raw[is_hc]

    print(f"HC samples      : {is_hc.sum()}")
    print(f"Candidate genes : {len(gene_names)}"
          f"  (det>={args.det_rate_min})")
    print(f"nu formula      : {args.nu_formula}")
    if args.binary_z_threshold is not None:
        print(f"Binary_Z threshold (det_rate) : < {args.binary_z_threshold:.0%}")

    if args.limit is not None:
        gene_names   = gene_names[:args.limit]
        gene_indices = gene_indices[:args.limit]
        print(f"Limited to {len(gene_names)} gene(s) for testing")

    folds = make_stratified_folds(strata_hc, n_splits=args.n_folds)

    ppc_path  = SAVE_DIR / "cv_zinb_ppc.pkl"   # per-sample mu, sigma, nu for PPC

    done_genes   = set() if args.no_resume else _load_done_genes(meta_path)
    zscores_dict = {}
    if not args.no_resume and zscores_path.exists():
        try:
            with open(zscores_path, "rb") as f:
                zscores_dict = pickle.load(f)
        except (EOFError, pickle.UnpicklingError):
            print(f"[Warning] {zscores_path} is corrupted (truncated write). Starting fresh.")
            zscores_path.unlink()

    ppc_dict = {}
    if not args.no_resume and ppc_path.exists():
        try:
            with open(ppc_path, "rb") as f:
                ppc_dict = pickle.load(f)
        except (EOFError, pickle.UnpicklingError):
            ppc_path.unlink()

    write_header = args.no_resume or not meta_path.exists()
    meta_file    = open(meta_path, "w" if args.no_resume else "a", newline="")
    meta_writer  = csv.DictWriter(meta_file, fieldnames=META_FIELDS)
    if write_header:
        meta_writer.writeheader()
        meta_file.flush()

    n_done = n_skipped = 0
    t_start = time.perf_counter()

    for i, (g_name, g_idx) in enumerate(zip(gene_names, gene_indices)):
        if g_name in done_genes:
            n_skipped += 1
            continue

        y_hc_gene  = Y_hc[:, g_idx]
        det_rate   = float((y_hc_gene > 0).mean())
        mean_count = float(y_hc_gene.mean())

        t0 = time.perf_counter()
        (z_bin, z_cond, z_full,
         mu_all, sigma_all, nu_all, fold_ok, n_removed) = eval_gene_cv(
            y_hc_gene, X_hc_scaled, folds, r_fit_fn, BIAS_COLUMNS, args
        )
        elapsed = time.perf_counter() - t0

        # ── W1 calibration metrics ────────────────────────────────
        w1_bin  = wasserstein1_normal(z_bin)
        w1_cond = wasserstein1_normal(z_cond)
        w1_full = wasserstein1_normal(z_full)

        # ── Distribution stats ─────────────────────────────────────
        m_b, s_b, sk_b, ku_b, nv_b = zscore_stats(z_bin)
        m_c, s_c, sk_c, ku_c, nv_c = zscore_stats(z_cond)
        m_f, s_f, sk_f, ku_f, nv_f = zscore_stats(z_full)

        # ── Primary z selection ────────────────────────────────────
        ptype     = decide_primary(det_rate, g_name, args.binary_z_threshold, nbi_lookup)
        z_primary = {"binary": z_bin, "count_cond": z_cond, "full": z_full}[ptype]

        row = {
            "gene": g_name, "n_hc": int(is_hc.sum()),
            "det_rate_hc": det_rate, "mean_count_hc": mean_count,
            "w1_binary":   w1_bin,
            "mean_binary_z":  m_b,  "std_binary_z":  s_b,
            "skew_binary_z":  sk_b, "kurt_binary_z":  ku_b, "n_binary": nv_b,
            "w1_count_cond": w1_cond,
            "mean_count_cond_z": m_c, "std_count_cond_z": s_c,
            "skew_count_cond_z": sk_c, "kurt_count_cond_z": ku_c, "n_count_cond": nv_c,
            "w1_full":     w1_full,
            "mean_full_z":  m_f,  "std_full_z":  s_f,
            "skew_full_z":  sk_f, "kurt_full_z":  ku_f, "n_full": nv_f,
            "mu_mean":    float(np.nanmean(mu_all)),
            "sigma_mean": float(np.nanmean(sigma_all)),
            "nu_mean":    float(np.nanmean(nu_all)),
            "n_removed":  n_removed,
            "primary_z_type": ptype,
            "fold_success_rate": fold_ok,
            "time_s": elapsed,
        }
        flags = detect_flags(row, int(is_hc.sum()), args.n_folds)
        row.update(flags)

        meta_writer.writerow(row)
        meta_file.flush()

        zscores_dict[g_name] = {
            "binary":     z_bin,
            "count_cond": z_cond,
            "full":       z_full,
            "primary":    z_primary,
        }
        with open(zscores_path, "wb") as f:
            pickle.dump(zscores_dict, f)

        ppc_dict[g_name] = {
            'mu':    mu_all.astype(np.float32),
            'sigma': sigma_all.astype(np.float32),
            'nu':    nu_all.astype(np.float32),
        }
        with open(ppc_path, "wb") as f:
            pickle.dump(ppc_dict, f)

        n_done += 1
        flag_str = " [FLAG]" if flags["any_flag"] else ""
        print(
            f"[{i+1:4d}/{len(gene_names)}] {g_name:<22s} "
            f"det={det_rate:.2f} rmvd={n_removed} | "
            f"W1 bin={w1_bin:.3f} cnd={w1_cond:.3f} full={w1_full:.3f} | "
            f"primary={ptype:<10s} folds={fold_ok:.0%} {elapsed:.1f}s{flag_str}"
        )

    meta_file.close()

    failures_path = SAVE_DIR / "cv_zinb_failures.csv"
    try:
        import pandas as pd
        df_all = pd.read_csv(meta_path)
        df_fail = df_all[df_all["any_flag"] == 1].copy()

        # Summarise which flags were triggered per gene
        flag_cols = [c for c in df_fail.columns if c.startswith("flag_") and c != "any_flag"]
        df_fail["failure_types"] = df_fail[flag_cols].apply(
            lambda r: "|".join(c.replace("flag_", "") for c in flag_cols if r[c] == 1),
            axis=1,
        )
        df_fail.to_csv(failures_path, index=False)
        n_fail = len(df_fail)
        print(f"\nFailures ({n_fail}/{n_done} genes) → {failures_path}")
        if n_fail > 0:
            top_flags = df_fail[flag_cols].sum().sort_values(ascending=False)
            print("  Flag counts:")
            for col, cnt in top_flags[top_flags > 0].items():
                print(f"    {col:<30s}: {int(cnt)}")
    except Exception as exc:
        print(f"[Warning] Failure summary failed: {exc}")

    print(f"\nDone. {n_done} genes, {n_skipped} skipped — "
          f"{time.perf_counter()-t_start:.1f}s total")
    print(f"Stats    → {meta_path}")
    print(f"Scores   → {zscores_path}")
    print(f"Failures → {failures_path}")


if __name__ == "__main__":
    main()
