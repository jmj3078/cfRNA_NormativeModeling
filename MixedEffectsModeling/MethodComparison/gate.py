"""Phase 0 fairness gate.

Re-runs the existing 3-arm CV comparison with only three changes -- shared evaluation
universe, shared perturbed-gene draw (by gene name), and BH restricted to that universe
-- to test whether the engine's margin over OUTRIDER survives them. Everything else
(oracle-mu scoring for the OUTRIDER arms, 500 perturbed genes, fold-shared draw) is left
exactly as it was so the diff isolates the three fixes.

The engine arm uses CV_Results_mixed/fold_params/, i.e. a genuine per-fold refit on the
same StratifiedKFold(5, seed 42) split the OUTRIDER folds were built from.

Run:  python gate.py
"""
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.dispersion_trend import load_trend
from MixedEffectsModeling.core.marginal_rqr import marginal_nb_rqr
from MixedEffectsModeling.core.shash import shash_transform_to_z
from MixedEffectsModeling.validation.cv_engine import _mu_alpha

from common import (CACHE, LOG2FCS, count_rejections, draw_perturbed, eval_universe,
                    injection_reference, nb_pvals, perturb, summarize, z_pvals)

N_PERTURB = 500
FOLD_PARAMS = config.CV_MIXED_DIR / "fold_params"


def load_counts(universe):
    # keyed on the universe -- a stale cache from a different gene set would silently
    # misalign every downstream column
    key = hashlib.sha1("\n".join(universe).encode()).hexdigest()[:12]
    path = CACHE / f"hc_counts_{key}.npz"
    if path.exists():
        d = np.load(path, allow_pickle=True)
        return d["samples"].tolist(), d["y"]
    df = pd.read_csv(config.OUTRIDER_COMPARISON_DIR / "hc_counts.csv", index_col=0)
    df = df.loc[universe]
    samples = df.columns.tolist()
    y = df.to_numpy(dtype=np.float64).T
    np.savez_compressed(path, samples=np.array(samples), y=y)
    return samples, y


def load_covariates(samples):
    with h5py.File(config.H5AD_PATH, "r") as f:
        names = [n.decode() if isinstance(n, bytes) else n for n in f["obs/_index"][()]]
        X = np.column_stack([f[f"obs/{c}"][()] for c in config.BIAS_COLUMNS]).astype(np.float64)
    pos = {n: i for i, n in enumerate(names)}
    return X[[pos[s] for s in samples]]


def engine_fold(fi, universe, tr_rows, te_rows, X_raw):
    """Per-fold refit coefficients + the fold's own train-fit SHASH, hoisted out of the
    log2fc loop. mu/alpha depend only on covariates, so they are computed once and reused
    across every injection level -- which also gives common random numbers across log2fc."""
    fits = pd.read_csv(FOLD_PARAMS / f"model_fold{fi}.csv").set_index("gene")
    shash = pd.read_csv(FOLD_PARAMS / f"shash_model_fold{fi}.csv").set_index("gene")
    alpha_fn = load_trend(config.DISPERSION_TREND_PATH)
    scaler = StandardScaler().fit(X_raw[tr_rows])
    Xa = np.column_stack([np.ones(len(te_rows)), scaler.transform(X_raw[te_rows])])

    mu_cols = [c for c in fits.columns if c.startswith("mu_coef_")]
    disp_cols = [c for c in fits.columns if c.startswith("disp_coef_")]
    params = {}
    for j, g in enumerate(universe):
        if g not in fits.index:
            continue
        row = fits.loc[g]
        if not bool(row["ok"]):
            continue
        mu, alpha = _mu_alpha(Xa, row[mu_cols].to_numpy(dtype=float),
                              row[disp_cols].to_numpy(dtype=float), row, alpha_fn)
        s = shash.loc[g] if g in shash.index and bool(shash.loc[g, "ok"]) else None
        params[j] = (mu, alpha, float(row["tau2"]), s)
    return params


def engine_z(fi, params, n_samples, n_genes, y_pert):
    z = np.full((n_samples, n_genes), np.nan)
    for j, (mu, alpha, tau2, s) in params.items():
        zg = marginal_nb_rqr(y_pert[:, j], mu, alpha, tau2, seed=1000 * fi + j)
        if s is not None:
            zg = shash_transform_to_z(zg, s.xi, s.eta, s.eps, s.delta)
        z[:, j] = np.clip(zg, -50, 50)
    return z


def engine_pvals(fi, params, n_samples, n_genes, y_pert):
    return z_pvals(engine_z(fi, params, n_samples, n_genes, y_pert))


def outrider_fold(arm_dir, fi, universe, test_names):
    mu = pd.read_csv(arm_dir / f"cv_fold{fi}_mu.csv", index_col=0).loc[test_names, universe]
    theta = pd.read_csv(arm_dir / f"cv_fold{fi}_theta.csv").set_index("gene")["theta"].loc[universe]
    return mu.to_numpy(dtype=np.float64), theta.to_numpy(dtype=np.float64)


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    universe = eval_universe()
    print(f"EVAL_UNIVERSE: {len(universe)} genes", flush=True)

    samples, Y = load_counts(universe)
    X_raw = load_covariates(samples)
    row_of = {s: i for i, s in enumerate(samples)}
    folds = json.load(open(config.OUTRIDER_COMPARISON_DIR / "cv_folds.json"))

    rows = []
    for fi_key in sorted(folds, key=int):
        fi = int(fi_key)
        tr_names, te_names = folds[fi_key]["train"], folds[fi_key]["test"]
        tr_rows = np.array([row_of[s] for s in tr_names])
        te_rows = np.array([row_of[s] for s in te_names])

        m_ref = injection_reference(Y[tr_rows], Y[te_rows])
        pert_idx = draw_perturbed(universe, N_PERTURB, fi)
        labels = np.zeros(len(universe), dtype=bool)
        labels[pert_idx] = True

        outr = {name: outrider_fold(d, fi, universe, te_names) for name, d in
                [("outrider_insample", config.OUTRIDER_COMPARISON_DIR),
                 ("outrider_heldout", config.OUTRIDER_HELD_OUT_DIR)]}
        params = engine_fold(fi, universe, tr_rows, te_rows, X_raw)

        for log2fc in LOG2FCS:
            y_pert = perturb(Y[te_rows], m_ref, pert_idx, log2fc)
            arms = {"engine_cv": engine_pvals(fi, params, len(te_rows), len(universe), y_pert)}
            for name, (mu, theta) in outr.items():
                arms[name] = nb_pvals(y_pert, mu, theta[None, :])
            for arm, p in arms.items():
                for c in count_rejections(p, labels):
                    rows.append(dict(arm=arm, fold=fi, log2fc=log2fc,
                                     n_samples=len(te_rows), n_universe=len(universe),
                                     n_perturbed=len(pert_idx), **c))
            print(f"  fold{fi} log2fc={log2fc} done", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(CACHE / "gate_counts.csv", index=False)
    summary = summarize(df)
    summary.to_csv(CACHE / "gate_summary.csv")
    print(summary.to_string())


if __name__ == "__main__":
    main()
