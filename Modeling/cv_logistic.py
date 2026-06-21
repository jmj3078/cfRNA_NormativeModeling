#!/usr/bin/env python
"""Stratified 5-fold CV benchmark for low-detectability genes using L2 Logistic Regression.

Applies to protein-coding genes where det_rate < det_rate_max (default 0.10).
For these genes, count-based models (NBI, ZINBI) have too little information;
instead we model whether the gene is detected at all (binary outcome) using
logistic regression, and score each sample with a Bernoulli randomized quantile
residual (RQR) z-score.

Z-score:
    p = sigmoid(X @ beta)          P(gene detected | covariates)
    y = 0 → u ~ Uniform(0,   1-p)
    y > 0 → u ~ Uniform(1-p, 1  )
    z = Φ⁻¹(u)

Under the correct logistic model, z ~ N(0,1).

Output (CV_Results/):
    cv_logistic_stats.csv   — per-gene AD test, z-score moments, failure flags
    cv_logistic_zscores.pkl — dict {gene: z_array (length N_hc)}

Resumable via --no-resume.

Usage:
    python cv_logistic.py --limit 10        # smoke-test
    python cv_logistic.py                   # full run
"""

import argparse
import csv
import pickle
import time
import warnings
from pathlib import Path

from joblib import Parallel, delayed

import numpy as np
from scipy.stats import anderson, norm, skew, kurtosis
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning
import scanpy as sc
from scipy.sparse import issparse

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

# ── Paths ──────────────────────────────────────────────────────────
BASE_DIR  = Path(__file__).resolve().parent
DATA_DIR  = BASE_DIR.parent / "OpenAccess_nfcore"
H5AD_PATH = DATA_DIR / "Merged_Processed_AnnData_with_Batch_Biases_QC_Status.h5ad"
SAVE_DIR  = BASE_DIR / "CV_Results"

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
DET_RATE_MAX   = 0.10   # upper bound — genes BELOW this threshold
DET_RATE_MIN   = 0.01   # lower bound — genes with <1% detection rarely add signal
N_SPLITS       = 5
LR_C           = 1.0    # inverse regularization strength (L2)
LR_MAX_ITER    = 1000

META_FIELDS = [
    "gene", "n_hc", "det_rate_hc", "n_detected",
    # Full-sample AD
    "ad_stat", "ad_crit5", "ad_pass",
    # Subsampled AD (power-controlled)
    "ad_sub_stat", "ad_sub_crit5", "ad_sub_pass",
    "mean_z", "std_z", "skew_z", "kurt_z", "n_valid",
    "p_detect_mean",     # mean predicted P(detected) across HC samples
    "fold_success_rate",
    # Flags
    "flag_fit_failure",  # fold_success_rate < 1.0
    "flag_p_extreme",    # p_detect_mean < 0.01 or > 0.99
    "flag_z_bias",       # |mean_z| > 0.3
    "flag_z_unstable",   # std_z > 1.5 or std_z < 0.5
    "flag_nonnormal",    # ads > 3.0
    "any_flag",
    "time_s",
]


# ── Data loading ───────────────────────────────────────────────────

def load_hc_data():
    adata = sc.read_h5ad(H5AD_PATH)
    adata = adata[adata.obs["QC_Passed"] == True]
    adata = adata[adata.obs["Phenotype_Processed"].notna()]
    adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
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


def select_low_det_genes(Y_raw, is_hc, pc_gene_names, pc_indices,
                          det_rate_min, det_rate_max):
    Y_hc  = Y_raw[is_hc][:, pc_indices]
    det_r = (Y_hc > 0).mean(axis=0)
    cand  = (det_r >= det_rate_min) & (det_r < det_rate_max)
    return np.array(pc_gene_names)[cand].tolist(), pc_indices[cand]


# ── Logistic model ─────────────────────────────────────────────────

def train_logistic(X_train, y_train, C=LR_C):
    lr = LogisticRegression(
        penalty="l2", C=C, solver="lbfgs",
        max_iter=LR_MAX_ITER, random_state=42,
    )
    y_bin = (y_train > 0).astype(int)
    if y_bin.sum() == 0 or y_bin.sum() == len(y_bin):
        return None   # no variation → can't fit
    lr.fit(X_train, y_bin)
    return lr


def bernoulli_rqr(p_detect: np.ndarray, y: np.ndarray, seed=None) -> np.ndarray:
    """Bernoulli RQR z-score.

    P(Y=0) = 1 - p_detect
    y=0  → u ~ Uniform(0,   1-p)
    y>0  → u ~ Uniform(1-p, 1  )
    """
    detected = (y > 0)
    a = np.where(detected, 1 - p_detect, 0.0)
    b = np.where(detected, 1.0,          1 - p_detect)
    lo = np.clip(np.minimum(a, b), 1e-8, 1 - 1e-8)
    hi = np.clip(np.maximum(a, b), 1e-8, 1 - 1e-8)
    rng = np.random.default_rng(seed)
    u = rng.uniform(lo, hi)
    return norm.ppf(u).astype(np.float32)


def score_logistic(lr, X_test, y_test, seed=None):
    if lr is None:
        return np.full(len(y_test), np.nan, dtype=np.float32)
    p_detect = lr.predict_proba(X_test)[:, 1]
    return bernoulli_rqr(p_detect, y_test, seed=seed)


# ── CV folds ───────────────────────────────────────────────────────

def make_stratified_folds(strata, n_splits=N_SPLITS):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    return list(skf.split(np.zeros(len(strata)), strata))


# ── Evaluation helpers ─────────────────────────────────────────────

def ad_test_normal(z_arr):
    z_valid = z_arr[np.isfinite(z_arr)]
    if len(z_valid) < 8:
        return np.nan, np.nan, False
    res   = anderson(z_valid, dist="norm")
    stat  = float(res.statistic)
    crit5 = float(res.critical_values[4])
    return stat, crit5, bool(stat <= crit5)


def ad_test_subsampled(z_arr, n_sub=200, n_boot=100, seed=42):
    z_valid = z_arr[np.isfinite(z_arr)]
    n = len(z_valid)
    if n < 8:
        return np.nan, np.nan, False
    if n <= n_sub:
        return ad_test_normal(z_arr)
    rng   = np.random.default_rng(seed)
    stats = [anderson(rng.choice(z_valid, size=n_sub, replace=False),
                      dist="norm").statistic for _ in range(n_boot)]
    crit5 = float(anderson(rng.standard_normal(n_sub), dist="norm").critical_values[4])
    median = float(np.median(stats))
    return median, crit5, bool(median <= crit5)


def zscore_stats(z_arr):
    z_valid = z_arr[np.isfinite(z_arr)]
    n = len(z_valid)
    if n < 2:
        return np.nan, np.nan, np.nan, np.nan, n
    return (float(z_valid.mean()), float(z_valid.std()),
            float(skew(z_valid)), float(kurtosis(z_valid)), n)


# ── Failure flags ─────────────────────────────────────────────────

def detect_flags(row):
    f = {}
    f["flag_fit_failure"] = int(row["fold_success_rate"] < 1.0)
    p = row["p_detect_mean"]
    f["flag_p_extreme"]   = int(not np.isfinite(p) or p < 0.01 or p > 0.99)
    f["flag_z_bias"]      = int(np.isfinite(row["mean_z"]) and abs(row["mean_z"]) > 0.3)
    std = row["std_z"]
    f["flag_z_unstable"]  = int(not np.isfinite(std) or std > 1.5 or std < 0.5)
    ads = row["ad_sub_stat"]
    f["flag_nonnormal"]   = int(not np.isfinite(ads) or ads > 3.0)
    f["any_flag"] = int(any(v for k, v in f.items() if k != "any_flag"))
    return f


# ── Per-gene CV ────────────────────────────────────────────────────

def eval_gene_cv(y_hc, X_hc_scaled, folds, C=LR_C):
    n         = len(y_hc)
    z_all     = np.full(n, np.nan, dtype=np.float32)
    p_all     = np.full(n, np.nan, dtype=np.float64)
    n_success = 0

    for fold_idx, (tr_idx, te_idx) in enumerate(folds):
        seed = 42 + fold_idx
        lr   = train_logistic(X_hc_scaled[tr_idx], y_hc[tr_idx], C=C)
        if lr is None:
            continue
        p_te         = lr.predict_proba(X_hc_scaled[te_idx])[:, 1]
        z_all[te_idx] = bernoulli_rqr(p_te, y_hc[te_idx], seed=seed)
        p_all[te_idx] = p_te
        n_success    += 1

    return z_all, p_all, n_success / len(folds)


# ── Parallel worker (top-level required for pickling) ─────────────

def _eval_one_gene(g_name, y_hc, X_hc_scaled, folds, n_hc, ad_n_sub, ad_n_boot, lr_C):
    """Single-gene evaluation. Called by joblib workers."""
    import time, numpy as np
    det_rate   = float((y_hc > 0).mean())
    n_detected = int((y_hc > 0).sum())
    t0 = time.perf_counter()
    z_all, p_all, fold_ok = eval_gene_cv(y_hc, X_hc_scaled, folds, C=lr_C)
    elapsed = time.perf_counter() - t0

    ad_stat,     ad_crit5,     ad_pass     = ad_test_normal(z_all)
    ad_sub_stat, ad_sub_crit5, ad_sub_pass = ad_test_subsampled(
        z_all, n_sub=ad_n_sub, n_boot=ad_n_boot)
    m_z, s_z, sk_z, ku_z, nv = zscore_stats(z_all)
    p_mean = float(np.nanmean(p_all))

    row = {
        "gene": g_name, "n_hc": n_hc,
        "det_rate_hc": det_rate, "n_detected": n_detected,
        "ad_stat":      ad_stat,     "ad_crit5":     ad_crit5,     "ad_pass":     int(ad_pass),
        "ad_sub_stat":  ad_sub_stat, "ad_sub_crit5": ad_sub_crit5, "ad_sub_pass": int(ad_sub_pass),
        "mean_z": m_z, "std_z": s_z, "skew_z": sk_z, "kurt_z": ku_z, "n_valid": nv,
        "p_detect_mean": p_mean,
        "fold_success_rate": fold_ok,
        "time_s": elapsed,
    }
    flags = detect_flags(row)
    row.update(flags)
    return g_name, z_all, row


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
    parser.add_argument("--det-rate-min",  type=float, default=DET_RATE_MIN,
                        help=f"minimum detection rate (default {DET_RATE_MIN})")
    parser.add_argument("--det-rate-max",  type=float, default=DET_RATE_MAX,
                        help=f"upper detection rate threshold (default {DET_RATE_MAX})")
    parser.add_argument("--n-folds",       type=int,   default=N_SPLITS)
    parser.add_argument("--lr-C",          type=float, default=LR_C,
                        help="logistic regression L2 regularization inverse strength (default 1.0)")
    parser.add_argument("--ad-n-sub",      type=int,   default=200)
    parser.add_argument("--ad-n-boot",     type=int,   default=100)
    parser.add_argument("--n-jobs",         type=int,   default=1,
                        help="parallel workers: -1 = all cores, 1 = sequential (default)")
    parser.add_argument("--no-resume",     action="store_true")
    parser.add_argument("--limit",         type=int,   default=None)
    args = parser.parse_args()

    SAVE_DIR.mkdir(parents=True, exist_ok=True)
    meta_path    = SAVE_DIR / "cv_logistic_stats.csv"
    zscores_path = SAVE_DIR / "cv_logistic_zscores.pkl"
    failures_path= SAVE_DIR / "cv_logistic_failures.csv"

    print("Loading HC data...")
    X_hc_scaled, Y_raw, is_hc, strata_hc, pc_gene_names, pc_indices = load_hc_data()
    gene_names, gene_indices = select_low_det_genes(
        Y_raw, is_hc, pc_gene_names, pc_indices,
        args.det_rate_min, args.det_rate_max,
    )
    Y_hc = Y_raw[is_hc]

    print(f"HC samples      : {is_hc.sum()}")
    print(f"Low-det genes   : {len(gene_names)}"
          f"  ({args.det_rate_min:.0%} ≤ det_rate < {args.det_rate_max:.0%})")

    if args.limit is not None:
        gene_names   = gene_names[:args.limit]
        gene_indices = gene_indices[:args.limit]
        print(f"Limited to {len(gene_names)} gene(s)")

    folds = make_stratified_folds(strata_hc, n_splits=args.n_folds)

    done_genes   = set() if args.no_resume else _load_done_genes(meta_path)
    zscores_dict = {}
    if not args.no_resume and zscores_path.exists():
        try:
            with open(zscores_path, "rb") as f:
                zscores_dict = pickle.load(f)
        except (EOFError, pickle.UnpicklingError):
            print(f"[Warning] {zscores_path} is corrupted (truncated write). Starting fresh.")
            zscores_path.unlink()

    write_header = args.no_resume or not meta_path.exists()
    meta_file    = open(meta_path, "w" if args.no_resume else "a", newline="")
    meta_writer  = csv.DictWriter(meta_file, fieldnames=META_FIELDS)
    if write_header:
        meta_writer.writeheader()
        meta_file.flush()

    n_hc = int(is_hc.sum())
    # Pre-filter: skip already-done genes
    todo = [(g, int(gi)) for g, gi in zip(gene_names, gene_indices)
            if g not in done_genes]
    n_skipped = len(gene_names) - len(todo)
    print(f"To evaluate: {len(todo)} genes  (skipped {n_skipped} already done)"
          f"  n_jobs={args.n_jobs}")

    t_start = time.perf_counter()

    # ── Parallel execution ─────────────────────────────────────────
    results = Parallel(
        n_jobs=args.n_jobs,
        backend="loky",
        verbose=0,
    )(
        delayed(_eval_one_gene)(
            g_name, Y_hc[:, g_idx], X_hc_scaled, folds,
            n_hc, args.ad_n_sub, args.ad_n_boot, args.lr_C,
        )
        for g_name, g_idx in todo
    )

    # ── Write results (sequential after parallel completes) ────────
    n_done = 0
    for g_name, z_all, row in results:
        meta_writer.writerow(row)
        zscores_dict[g_name] = z_all
        n_done += 1
        ads   = row["ad_sub_stat"]
        sub_ok = "ok" if row["ad_sub_pass"] else "F"
        flag_str = " [FLAG]" if row["any_flag"] else ""
        print(
            f"[{n_done:4d}/{len(todo)}] {g_name:<22s} "
            f"det={row['det_rate_hc']:.3f}({row['n_detected']:3d}) | "
            f"ADsub={ads:.3f}({sub_ok}) "
            f"mean={row['mean_z']:+.3f} std={row['std_z']:.3f} | "
            f"folds={row['fold_success_rate']:.0%} {row['time_s']:.2f}s{flag_str}"
        )

    meta_file.flush()
    meta_file.close()

    # Persist z-scores
    with open(zscores_path, "wb") as f:
        pickle.dump(zscores_dict, f)

    # ── Post-run failure summary ───────────────────────────────────
    try:
        import pandas as pd
        df_all  = pd.read_csv(meta_path)
        df_fail = df_all[df_all["any_flag"] == 1].copy()
        flag_cols = [c for c in df_fail.columns if c.startswith("flag_") and c != "any_flag"]
        df_fail["failure_types"] = df_fail[flag_cols].apply(
            lambda r: "|".join(c.replace("flag_", "") for c in flag_cols if r[c] == 1), axis=1)
        df_fail.to_csv(failures_path, index=False)
        print(f"\nFailures ({len(df_fail)}/{n_done}) → {failures_path}")
        if len(df_fail) > 0:
            for col, cnt in df_fail[flag_cols].sum().sort_values(ascending=False).items():
                if cnt > 0:
                    print(f"  {col:<30s}: {int(cnt)}")
    except Exception as e:
        print(f"[Warning] Failure summary failed: {e}")

    print(f"\nDone. {n_done} genes, {n_skipped} skipped — "
          f"{time.perf_counter()-t_start:.1f}s total")
    print(f"Stats    → {meta_path}")
    print(f"Scores   → {zscores_path}")


if __name__ == "__main__":
    main()
