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
import pickle
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as rpyn
import scanpy as sc
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter
from scipy.sparse import issparse
from scipy.stats import kurtosis, norm, skew
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=UserWarning, module="rpy2")

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
import config

BIAS_COLUMNS = config.BIAS_COLUMNS
STRATIFY_COL = config.MODELING_PARAMS["stratify_col"]
DET_RATE_MIN = config.MODELING_PARAMS["det_rate_min"]
N_SPLITS = config.MODELING_PARAMS["n_splits"]

META_FIELDS = [
    "gene", "n_hc", "det_rate_hc", "mean_count_hc",
    "w1",
    "mean_z", "std_z", "skew_z", "kurt_z", "n_valid",
    "mu_mean", "sigma_mean",
    "n_removed",
    "fold_success_rate",
    "time_s",
]


def init_r():
    """Source gamlss.r once and return the R fit function."""
    ro.r(f'source("{config.R_HELPER}")')
    return ro.globalenv["fit_gamlss_gene"]


def load_hc_data():
    adata = sc.read_h5ad(config.H5AD_PATH)
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
    pc_indices = np.where(is_pc)[0]

    strata_hc = strata_all[is_hc]
    return X_hc_scaled, Y_raw, is_hc, strata_hc, pc_gene_names, pc_indices


def select_candidate_genes(Y_raw, is_hc, pc_gene_names, pc_indices, det_rate_min, det_rate_max=None):
    Y_hc = Y_raw[is_hc][:, pc_indices]
    det_r = (Y_hc > 0).mean(axis=0)
    cand = det_r >= det_rate_min
    if det_rate_max is not None:
        cand &= det_r < det_rate_max
    return np.array(pc_gene_names)[cand].tolist(), pc_indices[cand]


def make_stratified_folds(strata, n_splits=N_SPLITS):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    return list(skf.split(np.zeros(len(strata)), strata))


def wasserstein1_normal(z_arr):
    """1st Wasserstein distance between empirical z distribution and N(0,1).

    W1 = mean(|z_(i:n) - Phi^{-1}((i-0.5)/n)|)

    N-robust: converges to the true population W1 without inflating with n.
    Threshold guideline: > 0.25 indicates poor calibration.
    """
    z_valid = z_arr[np.isfinite(z_arr)]
    n = len(z_valid)
    if n < 8:
        return np.nan
    z_ref = norm.ppf(np.linspace(1 / (2 * n), 1 - 1 / (2 * n), n))
    return float(np.mean(np.abs(np.sort(z_valid) - z_ref)))


def zscore_stats(z_arr):
    """(mean, std, skewness, excess_kurtosis, n_valid)."""
    z_valid = z_arr[np.isfinite(z_arr)]
    n = len(z_valid)
    if n < 2:
        return np.nan, np.nan, np.nan, np.nan, n
    return (float(z_valid.mean()), float(z_valid.std()),
            float(skew(z_valid)), float(kurtosis(z_valid)), n)


def _np_to_r_matrix(arr, col_names):
    """Convert a 2-D numpy array to an R matrix with dimnames."""
    with localconverter(ro.default_converter + rpyn.converter):
        r_mat = ro.conversion.py2rpy(np.ascontiguousarray(arr, dtype=np.float64))
    r_mat = ro.r["matrix"](
        r_mat,
        nrow=arr.shape[0], ncol=arr.shape[1],
        dimnames=ro.r["list"](ro.NULL, ro.StrVector(col_names)),
    )
    return r_mat


def _np_to_r_vec(arr):
    with localconverter(ro.default_converter + rpyn.converter):
        return ro.conversion.py2rpy(np.ascontiguousarray(arr, dtype=np.float64))


def eval_gene_cv(y_hc, X_hc_scaled, folds, r_fit_fn, col_names, args):
    """
    Run 5-fold CV for one gene.
    Returns (z_all, mu_all, sigma_all, fold_success_rate, n_removed_total).
    """
    n = len(y_hc)
    z_all = np.full(n, np.nan)
    mu_all = np.full(n, np.nan)
    sigma_all = np.full(n, np.nan)
    n_success = 0
    n_removed = 0

    for fold_idx, (tr_idx, te_idx) in enumerate(folds):
        seed_r = ro.IntVector([42 + fold_idx])
        n_cyc_r = ro.IntVector([50])
        outlier_z = ro.FloatVector([args.outlier_z])
        max_iter = ro.IntVector([args.max_iter])
        max_remove_frac = ro.FloatVector([args.max_remove_frac])
        lambda_sigma = ro.FloatVector([args.lambda_sigma])

        res = r_fit_fn(
            _np_to_r_vec(y_hc[tr_idx]),
            _np_to_r_vec(y_hc[te_idx]),
            _np_to_r_matrix(X_hc_scaled[tr_idx], col_names),
            _np_to_r_matrix(X_hc_scaled[te_idx], col_names),
            seed_r, n_cyc_r, outlier_z, max_iter, max_remove_frac, lambda_sigma,
        )

        if res.rx2("success")[0]:
            z_all[te_idx] = np.array(res.rx2("z"))
            mu_all[te_idx] = np.array(res.rx2("mu_test"))
            sigma_all[te_idx] = np.array(res.rx2("sigma_test"))
            n_removed += int(res.rx2("n_removed")[0])
            n_success += 1

    return z_all, mu_all, sigma_all, n_success / len(folds), n_removed


def _load_done_genes(meta_path):
    done = set()
    if meta_path.exists():
        with open(meta_path, newline="") as f:
            for row in csv.DictReader(f):
                done.add(row["gene"])
    return done


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--det-rate-min", type=float, default=DET_RATE_MIN)
    parser.add_argument("--det-rate-max", type=float, default=None,
                        help="upper bound on HC detection rate (exclusive); restricts to a band")
    parser.add_argument("--out-dir", type=str, default=None,
                        help="output directory (default config.CV_RESULTS_DIR); use to keep canonical results separate")
    parser.add_argument("--n-folds", type=int, default=N_SPLITS)
    parser.add_argument("--outlier-z", type=float, default=5.0,
                        help="remove training samples with |z_train| > this value (default 5.0)")
    parser.add_argument("--max-iter", type=int, default=2,
                        help="max outlier-removal iterations per fold (default 2)")
    parser.add_argument("--max-remove-frac", type=float, default=0.10,
                        help="max fraction of training samples removable per iteration (default 0.10)")
    parser.add_argument("--lambda-sigma", type=float, default=0.05,
                        help="L2 ridge penalty on sigma submodel coefficients (default 0.05; 0=disabled)")
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--limit", type=int, default=None,
                        help="only process first N genes (for testing)")
    args = parser.parse_args()

    out_dir = Path(args.out_dir) if args.out_dir else config.CV_RESULTS_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    meta_path = out_dir / "cv_gamlss_stats.csv"
    zscores_path = out_dir / "cv_gamlss_zscores.pkl"

    print("Initialising R / gamlss...")
    r_fit_fn = init_r()

    print("Loading HC data...")
    X_hc_scaled, Y_raw, is_hc, strata_hc, pc_gene_names, pc_indices = load_hc_data()
    gene_names, gene_indices = select_candidate_genes(
        Y_raw, is_hc, pc_gene_names, pc_indices,
        args.det_rate_min, args.det_rate_max,
    )
    Y_hc = Y_raw[is_hc]

    band = f">={args.det_rate_min}" if args.det_rate_max is None \
        else f"[{args.det_rate_min}, {args.det_rate_max})"
    print(f"HC samples      : {is_hc.sum()}")
    print(f"Candidate genes : {len(gene_names)}  (det {band})")

    if args.limit is not None:
        gene_names = gene_names[:args.limit]
        gene_indices = gene_indices[:args.limit]
        print(f"Limited to {len(gene_names)} gene(s) for testing")

    folds = make_stratified_folds(strata_hc, n_splits=args.n_folds)
    col_names = BIAS_COLUMNS

    ppc_path = out_dir / "cv_gamlss_nb_ppc.pkl"

    done_genes = set() if args.no_resume else _load_done_genes(meta_path)
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
    meta_file = open(meta_path, "w" if args.no_resume else "a", newline="")
    meta_writer = csv.DictWriter(meta_file, fieldnames=META_FIELDS)
    if write_header:
        meta_writer.writeheader()
        meta_file.flush()

    n_done = n_skipped = 0
    t_start = time.perf_counter()

    for i, (g_name, g_idx) in enumerate(zip(gene_names, gene_indices)):
        if g_name in done_genes:
            n_skipped += 1
            continue

        y_hc_gene = Y_hc[:, g_idx]
        det_rate = float((y_hc_gene > 0).mean())
        mean_count = float(y_hc_gene.mean())

        t0 = time.perf_counter()
        z_all, mu_all, sigma_all, fold_ok, n_removed = eval_gene_cv(
            y_hc_gene, X_hc_scaled, folds, r_fit_fn, col_names, args
        )
        elapsed = time.perf_counter() - t0

        w1 = wasserstein1_normal(z_all)
        m_z, s_z, sk_z, ku_z, nv = zscore_stats(z_all)

        mu_mean = float(np.nanmean(mu_all))
        sigma_mean = float(np.nanmean(sigma_all))

        row = {
            "gene": g_name, "n_hc": int(is_hc.sum()),
            "det_rate_hc": det_rate, "mean_count_hc": mean_count,
            "w1": w1,
            "mean_z": m_z, "std_z": s_z, "skew_z": sk_z, "kurt_z": ku_z, "n_valid": nv,
            "mu_mean": mu_mean, "sigma_mean": sigma_mean,
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
            'mu': mu_all.astype(np.float32),
            'sigma': sigma_all.astype(np.float32),
        }
        with open(ppc_path, "wb") as f:
            pickle.dump(ppc_dict, f)

        n_done += 1
        print(
            f"[{i+1:4d}/{len(gene_names)}] {g_name:<22s} "
            f"det={det_rate:.2f} mean={mean_count:7.1f} | "
            f"W1={w1:.3f} mean={m_z:+.3f} std={s_z:.3f} rmvd={n_removed} | "
            f"folds={fold_ok:.0%} {elapsed:.1f}s"
        )

    meta_file.close()
    print(f"\nDone. {n_done} genes, {n_skipped} skipped — "
          f"{time.perf_counter()-t_start:.1f}s total")
    print(f"Stats  → {meta_path}")
    print(f"Scores → {zscores_path}")


if __name__ == "__main__":
    main()
