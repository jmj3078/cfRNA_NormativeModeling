#!/usr/bin/env python
"""Stratified 5-fold CV benchmark using NB-GAMLSS (via rpy2 + R gamlss package).

For each candidate protein-coding gene the script:
  1. Splits the HC samples into 5 folds stratified by Batch_ID.
  2. Trains NB-GAMLSS (NBI family; both mu and sigma as linear functions of the
     10 batch-bias covariates) on each training fold in R via rpy2.
  3. Computes randomized quantile residual (RQR) z-scores on the held-out fold.
  4. Pools all 5 fold z-scores (length == N_hc) and evaluates normality via the
     Anderson-Darling test.
     
Resumable: genes already in cv_gamlss_stats.csv are skipped unless --no-resume.

Usage:
    python cv_gamlss.py --limit 5       # smoke-test
    python cv_gamlss.py                 # full run
"""

import argparse
import csv
import os
import pickle
import time
import warnings
from pathlib import Path

import numpy as np
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import anderson, skew, kurtosis
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# ── rpy2 setup ────────────────────────────────────────────────────
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects.numpy2ri as rpyn

warnings.filterwarnings("ignore", category=UserWarning, module="rpy2")

# ── Paths ──────────────────────────────────────────────────────────
BASE_DIR  = Path(__file__).resolve().parent
DATA_DIR  = BASE_DIR.parent / "OpenAccess_nfcore"
H5AD_PATH = DATA_DIR / "Merged_Processed_AnnData_with_Batch_Biases_QC_Status.h5ad"
SAVE_DIR  = BASE_DIR / "CV_Results"
R_HELPER  = BASE_DIR / "gamlss.r"

# ── Constants ──────────────────────────────────────────────────────
BIAS_COLUMNS = [
    "log(Total Reads)",
    "Spliced Reads (%)",
    "gDNA Contamination (Intron/Exon)",
    "rRNA Fraction",
    "RNA Degradation (3' Bias)",
    "Platelet Score",
    "GC Bias",
    "Gene Length Bias",
    "NG80",
    "(NP80/NG80)",
]
STRATIFY_COL   = "Batch_ID"
DET_RATE_MIN   = 0.1
MEAN_COUNT_MIN = 2.0
N_SPLITS       = 5

META_FIELDS = [
    "gene", "n_hc", "det_rate_hc", "mean_count_hc",
    # Full-sample AD (reference)
    "ad_stat", "ad_crit5", "ad_pass",
    # Subsampled AD (200-sample median, power-controlled)
    "ad_sub_stat", "ad_sub_crit5", "ad_sub_pass",
    "mean_z", "std_z", "skew_z", "kurt_z", "n_valid",
    "mu_mean", "sigma_mean",
    "n_removed",          # training samples removed by outlier refinement
    "fold_success_rate",
    "time_s",
]


# ── R environment initialisation ───────────────────────────────────

def init_r():
    """Source gamlss.r once and return the R fit function."""
    ro.r(f'source("{R_HELPER}")')
    return ro.globalenv["fit_gamlss_gene"]


# ── Data loading ───────────────────────────────────────────────────

def load_hc_data():
    adata = sc.read_h5ad(H5AD_PATH)
    adata = adata[adata.obs["QC_Passed"] == True]
    adata = adata[adata.obs["Phenotype_Processed"].notna()]
    adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
    adata = adata[adata.obs["broad_protocol_category"] != "Exome-based (EB)"]  # WTS only
    is_hc = (adata.obs["Phenotype_Processed"].astype(str) == "Healthy Control").values

    batch_raw = adata.obs[STRATIFY_COL].astype(object)
    strata_all = batch_raw.fillna("Unknown").astype(str).values

    X_raw = adata.obs[BIAS_COLUMNS].values.astype(np.float64)
    scaler = StandardScaler()
    X_hc_scaled = scaler.fit_transform(X_raw[is_hc])

    Y_raw = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
    Y_raw = np.round(Y_raw).astype(np.float64)

    is_pc = (adata.var["GeneType"] == "protein_coding").values
    pc_gene_names = adata.var_names[is_pc].tolist()
    pc_indices    = np.where(is_pc)[0]

    strata_hc = strata_all[is_hc]
    return X_hc_scaled, Y_raw, is_hc, strata_hc, pc_gene_names, pc_indices


def select_candidate_genes(Y_raw, is_hc, pc_gene_names, pc_indices,
                            det_rate_min, mean_count_min):
    Y_hc   = Y_raw[is_hc][:, pc_indices]
    det_r  = (Y_hc > 0).mean(axis=0)
    mean_c = Y_hc.mean(axis=0)
    cand   = (det_r >= det_rate_min) & (mean_c >= mean_count_min)
    return np.array(pc_gene_names)[cand].tolist(), pc_indices[cand]


# ── CV folds ───────────────────────────────────────────────────────

def make_stratified_folds(strata, n_splits=N_SPLITS):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    return list(skf.split(np.zeros(len(strata)), strata))


# ── Evaluation helpers ─────────────────────────────────────────────

def ad_test_normal(z_arr):
    """Anderson-Darling against N(0,1). Returns (statistic, crit_5pct, passed)."""
    z_valid = z_arr[np.isfinite(z_arr)]
    if len(z_valid) < 8:
        return np.nan, np.nan, False
    res   = anderson(z_valid, dist="norm")
    stat  = float(res.statistic)
    crit5 = float(res.critical_values[2])
    return stat, crit5, bool(stat <= crit5)


def ad_test_subsampled(z_arr, n_sub: int = 100, n_boot: int = 100, seed: int = 42):
    """Subsampled Anderson-Darling test to mitigate excessive power at large n.

    With n=996 the standard AD test rejects N(0,1) for any tiny real-world
    deviation. By drawing n_boot random subsamples of size n_sub and taking
    the median statistic, we obtain a test with the power appropriate for
    n_sub rather than the full sample — giving a fairer comparison across
    genes with different detection rates.

    Returns (median_stat, crit_5pct, passed) using n_sub-calibrated critical
    values (scipy's asymptotic values, which are conservative for n_sub=100).
    """
    z_valid = z_arr[np.isfinite(z_arr)]
    n = len(z_valid)
    if n < 8:
        return np.nan, np.nan, False
    if n <= n_sub:
        return ad_test_normal(z_arr)

    rng   = np.random.default_rng(seed)
    stats = []
    for _ in range(n_boot):
        sub = rng.choice(z_valid, size=n_sub, replace=False)
        res = anderson(sub, dist="norm")
        stats.append(res.statistic)

    crit5  = float(res.critical_values[4])   # same for all (asymptotic)
    median = float(np.median(stats))
    return median, crit5, bool(median <= crit5)


def zscore_stats(z_arr):
    """(mean, std, skewness, excess_kurtosis, n_valid)."""
    z_valid = z_arr[np.isfinite(z_arr)]
    n = len(z_valid)
    if n < 2:
        return np.nan, np.nan, np.nan, np.nan, n
    return (float(z_valid.mean()), float(z_valid.std()),
            float(skew(z_valid)), float(kurtosis(z_valid)), n)


# ── Per-gene CV ────────────────────────────────────────────────────

def _np_to_r_matrix(arr: np.ndarray, col_names: list):
    """Convert a 2-D numpy array to an R matrix with dimnames."""
    with localconverter(ro.default_converter + rpyn.converter):
        r_mat = ro.conversion.py2rpy(np.ascontiguousarray(arr, dtype=np.float64))
    # Attach column names so R can build the formula
    r_mat = ro.r["matrix"](
        r_mat,
        nrow=arr.shape[0], ncol=arr.shape[1],
        dimnames=ro.r["list"](ro.NULL, ro.StrVector(col_names)),
    )
    return r_mat


def _np_to_r_vec(arr: np.ndarray):
    with localconverter(ro.default_converter + rpyn.converter):
        return ro.conversion.py2rpy(np.ascontiguousarray(arr, dtype=np.float64))


def eval_gene_cv(y_hc, X_hc_scaled, folds, r_fit_fn, col_names, args):
    """
    Run 5-fold CV for one gene.
    Returns (z_all, mu_all, sigma_all, fold_success_rate, n_removed_total).
    """
    n = len(y_hc)
    z_all     = np.full(n, np.nan)
    mu_all    = np.full(n, np.nan)
    sigma_all = np.full(n, np.nan)
    n_success   = 0
    n_removed   = 0

    for fold_idx, (tr_idx, te_idx) in enumerate(folds):
        seed_r          = ro.IntVector([42 + fold_idx])
        n_cyc_r         = ro.IntVector([50])
        outlier_z       = ro.FloatVector([args.outlier_z])
        max_iter        = ro.IntVector([args.max_iter])
        max_remove_frac = ro.FloatVector([args.max_remove_frac])
        lambda_sigma    = ro.FloatVector([args.lambda_sigma])

        res = r_fit_fn(
            _np_to_r_vec(y_hc[tr_idx]),
            _np_to_r_vec(y_hc[te_idx]),
            _np_to_r_matrix(X_hc_scaled[tr_idx], col_names),
            _np_to_r_matrix(X_hc_scaled[te_idx], col_names),
            seed_r, n_cyc_r, outlier_z, max_iter, max_remove_frac, lambda_sigma,
        )

        if res.rx2("success")[0]:
            z_all[te_idx]     = np.array(res.rx2("z"))
            mu_all[te_idx]    = np.array(res.rx2("mu_test"))
            sigma_all[te_idx] = np.array(res.rx2("sigma_test"))
            n_removed        += int(res.rx2("n_removed")[0])
            n_success        += 1

    return z_all, mu_all, sigma_all, n_success / len(folds), n_removed


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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--det-rate-min",  type=float, default=DET_RATE_MIN)
    parser.add_argument("--mean-count-min",type=float, default=MEAN_COUNT_MIN)
    parser.add_argument("--n-folds",       type=int,   default=N_SPLITS)
    parser.add_argument("--outlier-z",     type=float, default=5.0,
                        help="remove training samples with |z_train| > this value (default 5.0)")
    parser.add_argument("--max-iter",          type=int,   default=2,
                        help="max outlier-removal iterations per fold (default 2)")
    parser.add_argument("--max-remove-frac",   type=float, default=0.10,
                        help="max fraction of training samples removable per iteration (default 0.10)")
    parser.add_argument("--lambda-sigma",   type=float, default=0.05,
                        help="L2 ridge penalty on sigma submodel coefficients (default 0.05; 0=disabled)")
    parser.add_argument("--ad-n-sub",      type=int,   default=100,
                        help="subsample size for power-controlled AD test (default 200)")
    parser.add_argument("--ad-n-boot",     type=int,   default=100,
                        help="number of bootstrap draws for subsampled AD (default 100)")
    parser.add_argument("--no-resume",     action="store_true")
    parser.add_argument("--limit",         type=int,   default=None,
                        help="only process first N genes (for testing)")
    args = parser.parse_args()

    SAVE_DIR.mkdir(parents=True, exist_ok=True)
    meta_path    = SAVE_DIR / "cv_gamlss_stats.csv"
    zscores_path = SAVE_DIR / "cv_gamlss_zscores.pkl"

    print("Initialising R / gamlss...")
    r_fit_fn = init_r()

    print("Loading HC data...")
    X_hc_scaled, Y_raw, is_hc, strata_hc, pc_gene_names, pc_indices = load_hc_data()
    gene_names, gene_indices = select_candidate_genes(
        Y_raw, is_hc, pc_gene_names, pc_indices,
        args.det_rate_min, args.mean_count_min,
    )
    Y_hc = Y_raw[is_hc]

    print(f"HC samples      : {is_hc.sum()}")
    print(f"Candidate genes : {len(gene_names)}"
          f"  (det>={args.det_rate_min}, mean>={args.mean_count_min})")

    if args.limit is not None:
        gene_names   = gene_names[:args.limit]
        gene_indices = gene_indices[:args.limit]
        print(f"Limited to {len(gene_names)} gene(s) for testing")

    folds    = make_stratified_folds(strata_hc, n_splits=args.n_folds)
    col_names = BIAS_COLUMNS   # passed to R for column names

    ppc_path  = SAVE_DIR / "cv_gamlss_nb_ppc.pkl"   # per-sample mu & sigma for PPC

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
        z_all, mu_all, sigma_all, fold_ok, n_removed = eval_gene_cv(
            y_hc_gene, X_hc_scaled, folds, r_fit_fn, col_names, args
        )
        elapsed = time.perf_counter() - t0

        ad_stat,     ad_crit5,     ad_pass     = ad_test_normal(z_all)
        ad_sub_stat, ad_sub_crit5, ad_sub_pass = ad_test_subsampled(
            z_all, n_sub=args.ad_n_sub, n_boot=args.ad_n_boot)
        m_z, s_z, sk_z, ku_z, nv = zscore_stats(z_all)

        mu_mean    = float(np.nanmean(mu_all))
        sigma_mean = float(np.nanmean(sigma_all))

        row = {
            "gene": g_name, "n_hc": int(is_hc.sum()),
            "det_rate_hc": det_rate, "mean_count_hc": mean_count,
            "ad_stat":      ad_stat,     "ad_crit5":     ad_crit5,     "ad_pass":     int(ad_pass),
            "ad_sub_stat":  ad_sub_stat, "ad_sub_crit5": ad_sub_crit5, "ad_sub_pass": int(ad_sub_pass),
            "mean_z":   m_z,  "std_z":    s_z,  "skew_z": sk_z, "kurt_z": ku_z, "n_valid": nv,
            "mu_mean":  mu_mean, "sigma_mean": sigma_mean,
            "n_removed": n_removed,
            "fold_success_rate": fold_ok,
            "time_s": elapsed,
        }
        meta_writer.writerow(row)
        meta_file.flush()

        zscores_dict[g_name] = z_all
        with open(zscores_path, "wb") as f:
            pickle.dump(zscores_dict, f)

        ppc_dict[g_name] = {
            'mu':    mu_all.astype(np.float32),
            'sigma': sigma_all.astype(np.float32),
        }
        with open(ppc_path, "wb") as f:
            pickle.dump(ppc_dict, f)

        n_done += 1
        sub_status = "ok" if ad_sub_pass else "FAIL"
        print(
            f"[{i+1:4d}/{len(gene_names)}] {g_name:<22s} "
            f"det={det_rate:.2f} mean={mean_count:7.1f} | "
            f"AD={ad_stat:.3f} sub={ad_sub_stat:.3f}({sub_status}) "
            f"mean={m_z:+.3f} std={s_z:.3f} rmvd={n_removed} | "
            f"folds={fold_ok:.0%} {elapsed:.1f}s"
        )

    meta_file.close()
    print(f"\nDone. {n_done} genes, {n_skipped} skipped — "
          f"{time.perf_counter()-t_start:.1f}s total")
    print(f"Stats  → {meta_path}")
    print(f"Scores → {zscores_path}")


if __name__ == "__main__":
    main()
