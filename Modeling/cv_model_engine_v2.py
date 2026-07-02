#!/usr/bin/env python
"""Stratified 5-fold CV calibration check for the trained NormativeModelEngineV2.

The route (pool/model) AND the specific stage within route "model" (nbi,
nb_fixed with its GAIC full-vs-intercept choice, or the final intercept stage --
see stage/mean_model_chosen in training_summary.csv) each gene ended up on are
taken as FIXED from the full-data training -- route/stage selection happens once
on all HC data, and CV only re-evaluates held-out calibration for that same
model, never re-decides the chain. This mirrors cv_gamlss_nb.py's
W1/mean/std/skew/kurt diagnostics, computed per gene from pooled held-out z.

Usage:
    python cv_model_engine_v2.py               # full run
    python cv_model_engine_v2.py --limit 300   # smoke test
"""

import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.robjects as ro
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import kurtosis, norm, skew
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config
from model_engine_v2 import (_nb_rqr, _nbi_rqr_from_coeffs, _poisson_rqr,
                             _to_r_matrix, _to_r_vec, fit_intercept_only_gene, fit_route_b_gene)
from dispersion_trend import load_trend

MP2 = config.MODELING_PARAMS_V2
STRATIFY_COL = MP2["stratify_col"]
N_SPLITS = MP2["n_splits"]


def load_hc():
    adata = sc.read_h5ad(config.H5AD_PATH)
    m = ((adata.obs["QC_Passed"] == True) & (adata.obs["Phenotype_Processed"].notna()) &
         (adata.obs["Phenotype_Processed"] != "Unknown") &
         (adata.obs["broad_protocol_category"] != "Exome-based (EB)"))
    a = adata[m]
    is_hc = (a.obs["Phenotype_Processed"].astype(str) == "Healthy Control").values
    X_raw = a.obs[config.BIAS_COLUMNS].values.astype(np.float64)[is_hc]
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X_raw)
    Y = a.X.toarray() if issparse(a.X) else np.asarray(a.X)
    Y = np.round(Y[is_hc]).astype(np.float64)
    strata = a.obs[STRATIFY_COL].astype(object).fillna("NA").astype(str).values[is_hc]
    var_names = np.array(a.var_names.tolist())
    return Xs, Y, var_names, strata


def w1_normal(z):
    v = z[np.isfinite(z)]
    n = len(v)
    if n < 8:
        return np.nan
    ref = norm.ppf(np.linspace(1 / (2 * n), 1 - 1 / (2 * n), n))
    return float(np.mean(np.abs(np.sort(v) - ref)))


def z_stats(z):
    v = z[np.isfinite(z)]
    n = len(v)
    if n < 2:
        return dict(w1=np.nan, mean_z=np.nan, std_z=np.nan, skew_z=np.nan, kurt_z=np.nan, n_valid=n)
    return dict(w1=w1_normal(v), mean_z=float(v.mean()), std_z=float(v.std()),
               skew_z=float(skew(v)), kurt_z=float(kurtosis(v)), n_valid=n)


def cv_pool(y, Xs, folds, mean_hc_full, rare_glm_full, seed):
    """Held-out z for a single route "pool" (rare-pooled) gene, refitting the
    shared pooled beta per fold is out of scope here -- reuse the full-data
    pooled GLM (consistent with 'route decided from full data' but scored on
    held-out y)."""
    n = len(y)
    z = np.full(n, np.nan)
    eps = rare_glm_full["eps"]
    for fi, (tr, te) in enumerate(folds):
        Xc = np.column_stack([np.ones(len(te)), Xs[te]])
        mu = np.clip((mean_hc_full + eps) * np.exp(Xc @ rare_glm_full["beta"]), 1e-12, 1e8)
        if rare_glm_full["family"] == "poisson":
            z[te] = _poisson_rqr(y[te], mu, seed + fi)
        else:
            z[te] = _nb_rqr(y[te], mu, rare_glm_full["alpha"], seed + fi)
    return z


def cv_nb_fixed(y, Xs, folds, alpha_fn, outlier_z, max_iter, max_remove_frac, seed,
                beta_explode_thr=None, gaic_k=None):
    """Re-fits stage nb_fixed (full-vs-intercept GAIC comparison) per fold,
    falling through to fit_intercept_only_gene if the full IRLS diverges."""
    n = len(y)
    z = np.full(n, np.nan)
    for fi, (tr, te) in enumerate(folds):
        res = fit_route_b_gene(y[tr], Xs[tr], alpha_fn, outlier_z, max_iter, max_remove_frac,
                               beta_explode_thr=beta_explode_thr, gaic_k=gaic_k)
        Xa_te = np.column_stack([np.ones(len(te)), Xs[te]])
        if res["success"]:
            if res["chosen"] == "full":
                mu = np.clip(np.exp(Xa_te @ res["beta"]), 1e-6, 1e8)
            else:
                mu = np.full(len(te), np.exp(res["beta_null"][0])).clip(1e-6, 1e8)
            z[te] = _nb_rqr(y[te], mu, res["alpha"], seed + fi)
        else:
            fb = fit_intercept_only_gene(y[tr], alpha_fn)
            if fb["success"]:
                mu = np.full(len(te), np.exp(fb["beta"][0])).clip(1e-6, 1e8)
                z[te] = _nb_rqr(y[te], mu, fb["alpha"], seed + fi)
    return z


def cv_intercept(y, alpha_fn, folds, seed):
    """Re-evaluates genes whose final stage IS "intercept"
    (stage == 'intercept' in training_summary.csv): closed-form,
    per fold, mirroring fit_intercept_only_gene exactly."""
    n = len(y)
    z = np.full(n, np.nan)
    for fi, (tr, te) in enumerate(folds):
        res = fit_intercept_only_gene(y[tr], alpha_fn)
        if not res["success"]:
            continue
        mu = np.full(len(te), np.exp(res["beta"][0])).clip(1e-6, 1e8)
        z[te] = _nb_rqr(y[te], mu, res["alpha"], seed + fi)
    return z


def cv_nbi(y, Xs, folds, r_fit_fn, col_names, outlier_z, max_iter, max_remove_frac,
          lambda_sigma, seed):
    n = len(y)
    z = np.full(n, np.nan)
    for fi, (tr, te) in enumerate(folds):
        try:
            res = r_fit_fn(
                _to_r_vec(y[tr]), _to_r_vec(y[te]),
                _to_r_matrix(Xs[tr], col_names), _to_r_matrix(Xs[te], col_names),
                ro.IntVector([seed + fi]), ro.IntVector([50]),
                ro.FloatVector([outlier_z]), ro.IntVector([max_iter]),
                ro.FloatVector([max_remove_frac]), ro.FloatVector([lambda_sigma]),
            )
            if res.rx2("success")[0]:
                z[te] = np.array(res.rx2("z"))
        except Exception:
            pass
    return z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    out_dir = config.CV_RESULTS_DIR_V2
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading trained engine summary + config...")
    summary_path = config.ENGINE_DIR_V2 / "training_summary.csv"
    summary = pd.read_csv(summary_path, index_col="gene")
    summary = summary[summary["attempted"] & (summary["route"] != "excluded")]
    if args.limit:
        summary = summary.iloc[:args.limit]

    # Use the parameters the engine actually trained with (config.pkl), NOT the
    # current config.MODELING_PARAMS_V2 -- config.py may have changed since training,
    # and CV must re-evaluate under the same settings that produced the routes.
    with open(config.ENGINE_DIR_V2 / "config.pkl", "rb") as f:
        engine_cfg = pickle.load(f)
    print(f"  engine config: {engine_cfg}")

    with open(config.ENGINE_DIR_V2 / "rare_glm.pkl", "rb") as f:
        rare_glm_full = pickle.load(f)
    alpha_fn = load_trend()

    print("Initialising R / gamlss...")
    ro.r(f'source("{config.R_HELPER}")')
    r_fit_fn = ro.globalenv["fit_gamlss_gene"]

    print("Loading HC data...")
    Xs, Y, var_names, strata = load_hc()
    name2col = {g: i for i, g in enumerate(var_names)}
    n_hc = Xs.shape[0]
    folds = list(StratifiedKFold(N_SPLITS, shuffle=True, random_state=args.seed)
                .split(np.zeros(n_hc), strata))

    stage_key = summary["stage"].where(summary["stage"] != "", summary["route"])
    print(f"Starting CV on {len(summary)} genes  "
         f"(stage counts: {stage_key.value_counts().to_dict()})")

    rows = []
    zdict = {}
    route_counts = {}
    t0 = time.perf_counter()
    for i, (gene, row) in enumerate(summary.iterrows()):
        j = name2col.get(gene)
        if j is None:
            continue
        y = Y[:, j]
        route = row["route"]
        stage = row.get("stage", "")
        if route == "pool":
            z = cv_pool(y, Xs, folds, y.mean(), rare_glm_full, args.seed)
        elif route == "model" and stage == "intercept":
            z = cv_intercept(y, alpha_fn, folds, args.seed)
        elif route == "model" and stage == "nb_fixed":
            z = cv_nb_fixed(y, Xs, folds, alpha_fn, engine_cfg["outlier_z"],
                            engine_cfg["max_outlier_iter"], engine_cfg["max_remove_frac"], args.seed,
                            beta_explode_thr=engine_cfg["beta_explode_thr"], gaic_k=engine_cfg["gaic_k"])
        elif route == "model" and stage == "nbi":
            z = cv_nbi(y, Xs, folds, r_fit_fn, config.BIAS_COLUMNS, engine_cfg["outlier_z"],
                      engine_cfg["max_outlier_iter"], engine_cfg["max_remove_frac"],
                      engine_cfg["ridge_lambda_sigma"], args.seed)
        else:
            continue
        zdict[gene] = z.astype(np.float32)
        st = z_stats(z)
        st.update(gene=gene, route=route, stage=stage, nz=int(row["nz"]))
        rows.append(st)
        route_counts[stage or route] = route_counts.get(stage or route, 0) + 1

        if (i + 1) % 50 == 0 or (i + 1) == len(summary):
            elapsed = time.perf_counter() - t0
            rate = (i + 1) / elapsed
            eta = (len(summary) - (i + 1)) / rate if rate > 0 else float("nan")
            print(f"  [{i+1}/{len(summary)}]  elapsed={elapsed:.1f}s  "
                 f"rate={rate:.2f} genes/s  eta={eta/60:.1f}min  routes_so_far={route_counts}",
                 flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "cv_stats.csv", index=False)
    with open(out_dir / "cv_zscores.pkl", "wb") as f:
        pickle.dump(zdict, f)

    print("\nCalibration by stage (median):")
    group_key = df["stage"].where(df["stage"] != "", df["route"])
    print(df.groupby(group_key)[["w1", "mean_z", "std_z", "skew_z", "kurt_z"]].median().to_string())
    print(f"\nSaved -> {out_dir}/cv_stats.csv, cv_zscores.pkl")


if __name__ == "__main__":
    main()
