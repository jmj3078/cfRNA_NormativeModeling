"""Normative-mode sweep: every method is fit on train HC only, frozen, then projects
held-out samples whose counts carry an injected signal. Records performance per model
across the (log2fc x FDR q) grid, plus the full injected-count and Z matrices.

Two splits, both genuinely inductive:
  cv   -- StratifiedKFold(5, seed 42), the engine's own CV split. Held-out samples are
          unseen, but every batch is represented in train.
  lobo -- leave-one-batch-out on the five batches 3_detection_limit.ipynb evaluates the
          engine on. Train has never seen the held-out batch at all: the normative claim's
          actual setting.

Normative mode means full-pipeline: the perturbed counts go through each method's own
scoring path, so size factors and any expression-derived latent reference are recomputed
from the perturbed sample -- what happens when a frozen model meets a real outlier. The
engine's covariates come from h5ad obs and cannot move, so for it this is identical to
scoring unperturbed covariates; that asymmetry is a property of the designs, not of this
harness.

Every arm receives bit-identical injected counts, shares one evaluation universe, and has
BH applied within that universe. Each grid cell is cached, so a rerun only does what is
missing.

Run:  python sweep.py [cv|lobo|both]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.dispersion_trend import load_trend
from MixedEffectsModeling.core.marginal_rqr import marginal_nb_rqr
from MixedEffectsModeling.core.shash import shash_transform_to_z
from MixedEffectsModeling.validation.cv_engine import _mu_alpha

import outsingle
import peer_arm
from common import (CACHE, LOG2FCS, Q_LEVELS, count_rejections, draw_perturbed,
                    injection_reference, perturb, z_pvals)
from frozen_latent import FITS, load_fit, score_frozen
from gate import load_counts, load_covariates

N_PERTURB = 500
LATENT_ARMS = ["autoencoder", "pca"]
MATRIX_DIR = CACHE / "matrices"
CELL_DIR = CACHE / "cells"
SPLITS = ("cv", "lobo")


def fold_spec(split):
    """[(key, tag, train_names, test_names, engine_source)] for the split."""
    if split == "cv":
        folds = json.load(open(config.OUTRIDER_COMPARISON_DIR / "cv_folds.json"))
        return [(k, f"fold{k}", folds[k]["train"], folds[k]["test"],
                 config.CV_MIXED_DIR / "fold_params")
                for k in sorted(folds, key=int)]
    folds = json.load(open(CACHE / "lobo_folds.json"))
    return [(k, f"lobo_{folds[k]['batch_dir']}", folds[k]["train"], folds[k]["test"],
             config.LOBO_MIXED_DIR / folds[k]["batch_dir"])
            for k in sorted(folds)]


def eval_universe_nm():
    """Engine model-route genes intersected with every latent fit's surviving gene set,
    across BOTH splits, so cv and lobo numbers sit on the same gene set."""
    path = CACHE / "eval_universe_nm.txt"
    if path.exists():
        return path.read_text().split()
    summary = pd.read_csv(config.ENGINE_MIXED_DIR / "training_summary.csv")
    genes = set(summary.loc[summary["ok"].fillna(False) & (summary["route"] == "model"), "gene"])
    for split in SPLITS:
        for _, tag, _, _, _ in fold_spec(split):
            for impl in LATENT_ARMS:
                f = FITS / f"{impl}_{tag}_genes.csv"
                if not f.exists():
                    raise FileNotFoundError(f"missing latent fit {f} -- run fit_latent.R first")
                genes &= set(pd.read_csv(f)["gene"])
    universe = sorted(genes)
    path.write_text("\n".join(universe))
    return universe


def engine_params(source, tag, universe, tr_rows, te_rows, X_raw):
    """Per-fold frozen GLMM coefficients + that fold's train-fit SHASH. mu/alpha depend
    only on covariates, so they are computed once and reused across every injection level,
    which also gives common random numbers across log2fc."""
    if (source / "model_fits.csv").exists():        # LOBO layout
        fits = pd.read_csv(source / "model_fits.csv").set_index("gene")
        shash = pd.read_csv(source / "shash_params.csv").set_index("gene")
    else:                                            # CV fold_params layout
        fi = tag.replace("fold", "")
        fits = pd.read_csv(source / f"model_fold{fi}.csv").set_index("gene")
        shash = pd.read_csv(source / f"shash_model_fold{fi}.csv").set_index("gene")

    alpha_fn = load_trend(config.DISPERSION_TREND_PATH)
    scaler = StandardScaler().fit(X_raw[tr_rows])
    Xa = np.column_stack([np.ones(len(te_rows)), scaler.transform(X_raw[te_rows])])
    mu_cols = [c for c in fits.columns if c.startswith("mu_coef_")]
    disp_cols = [c for c in fits.columns if c.startswith("disp_coef_")]

    params = {}
    for j, g in enumerate(universe):
        if g not in fits.index or not bool(fits.loc[g, "ok"]):
            continue
        row = fits.loc[g]
        mu, alpha = _mu_alpha(Xa, row[mu_cols].to_numpy(dtype=float),
                              row[disp_cols].to_numpy(dtype=float), row, alpha_fn)
        s = shash.loc[g] if g in shash.index and bool(shash.loc[g, "ok"]) else None
        params[j] = (mu, alpha, float(row["tau2"]), s)
    return params


def engine_z(seed_base, params, n_samples, n_genes, y_pert):
    z = np.full((n_samples, n_genes), np.nan)
    for j, (mu, alpha, tau2, s) in params.items():
        zg = marginal_nb_rqr(y_pert[:, j], mu, alpha, tau2, seed=seed_base + j)
        if s is not None:
            zg = shash_transform_to_z(zg, s.xi, s.eta, s.eps, s.delta)
        z[:, j] = np.clip(zg, -50, 50)
    return z


def save_cell(split, tag, log2fc, samples, universe, y_pert, labels, z, p):
    """Persist the injected counts and every arm's Z and p, so later questions never need
    the sweep re-run."""
    MATRIX_DIR.mkdir(parents=True, exist_ok=True)
    payload = dict(samples=np.array(samples), genes=np.array(universe),
                   y_pert=y_pert.astype(np.int32), labels=labels,
                   log2fc=log2fc, split=split, tag=tag)
    for arm in z:
        payload[f"z__{arm}"] = z[arm].astype(np.float32)
        payload[f"p__{arm}"] = p[arm].astype(np.float32)
    np.savez_compressed(MATRIX_DIR / f"{split}_{tag}_log2fc{log2fc:g}.npz", **payload)


def run_split(split, universe, all_counts, samples_u, Y_u, X_raw):
    CELL_DIR.mkdir(parents=True, exist_ok=True)
    row_of = {s: i for i, s in enumerate(samples_u)}
    out = []

    for seed_i, (key, tag, tr_names, te_names, src) in enumerate(fold_spec(split)):
        tr_rows = np.array([row_of[s] for s in tr_names])
        te_rows = np.array([row_of[s] for s in te_names])

        pert_idx = draw_perturbed(universe, N_PERTURB, seed_i)
        pert_genes = [universe[j] for j in pert_idx]
        labels = np.zeros(len(universe), dtype=bool)
        labels[pert_idx] = True
        m_ref = injection_reference(Y_u[tr_rows], Y_u[te_rows])

        pending = [f for f in LOG2FCS
                   if not (CELL_DIR / f"{split}_{tag}_log2fc{f:g}.csv").exists()]
        if not pending:
            print(f"  {split}/{tag}: all cells cached", flush=True)
            out += [pd.read_csv(CELL_DIR / f"{split}_{tag}_log2fc{f:g}.csv") for f in LOG2FCS]
            continue

        params = engine_params(src, tag, universe, tr_rows, te_rows, X_raw)
        latent = {}
        for impl in LATENT_ARMS:
            fit = load_fit(impl, tag)
            pos = {g: i for i, g in enumerate(fit["genes"])}
            latent[impl] = dict(
                fit=fit,
                y_te=all_counts.loc[fit["genes"], te_names].to_numpy(dtype=np.float64).T,
                y_tr=all_counts.loc[fit["genes"], tr_names].to_numpy(dtype=np.float64).T,
                gi=[pos[g] for g in universe], pi=[pos[g] for g in pert_genes])

        # OutSingle and PEER have no expression filter of their own, so they get the same
        # input matrix the autoencoder was fit on: the three then differ only in how the
        # latent basis is learned.
        ae = latent["autoencoder"]
        residual_fits = {"outsingle": outsingle.fit_train(ae["y_tr"], f"{split}_{tag}"),
                         "peer": peer_arm.fit_train(ae["y_tr"], f"{split}_{tag}")}

        for log2fc in LOG2FCS:
            cell_csv = CELL_DIR / f"{split}_{tag}_log2fc{log2fc:g}.csv"
            if cell_csv.exists():
                out.append(pd.read_csv(cell_csv))
                continue

            y_pert = perturb(Y_u[te_rows], m_ref, pert_idx, log2fc)
            z, p = {}, {}
            z["engine"] = engine_z(1000 * (seed_i + 1), params, len(te_rows), len(universe), y_pert)
            p["engine"] = z_pvals(z["engine"])

            for impl, L in latent.items():
                y_arm = L["y_te"].copy()
                y_arm[:, L["pi"]] = y_pert[:, pert_idx]   # bit-identical injected counts
                o = score_frozen(L["fit"], y_arm)
                z[impl], p[impl] = o["z"][:, L["gi"]], o["pval"][:, L["gi"]]

            y_res = ae["y_te"].copy()
            y_res[:, ae["pi"]] = y_pert[:, pert_idx]
            for arm, mod in (("outsingle", outsingle), ("peer", peer_arm)):
                o = mod.score_frozen(residual_fits[arm], y_res)
                z[arm], p[arm] = o["z"][:, ae["gi"]], o["pval"][:, ae["gi"]]

            save_cell(split, tag, log2fc, te_names, universe, y_pert, labels, z, p)
            cell_rows = [dict(split=split, arm=arm, fold=key, log2fc=log2fc,
                              n_samples=len(te_rows), n_universe=len(universe),
                              n_perturbed=len(pert_idx), **c)
                         for arm, pv in p.items() for c in count_rejections(pv, labels, Q_LEVELS)]
            df = pd.DataFrame(cell_rows)
            df.to_csv(cell_csv, index=False)
            out.append(df)
            print(f"  {split}/{tag} log2fc={log2fc} done", flush=True)
    return out


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    splits = SPLITS if which == "both" else (which,)

    universe = eval_universe_nm()
    print(f"EVAL_UNIVERSE (NM, shared across splits): {len(universe)} genes", flush=True)

    all_counts = pd.read_csv(config.OUTRIDER_COMPARISON_DIR / "hc_counts.csv", index_col=0)
    samples_u, Y_u = load_counts(universe)
    X_raw = load_covariates(samples_u)

    parts = []
    for split in splits:
        parts += run_split(split, universe, all_counts, samples_u, Y_u, X_raw)

    df = pd.concat(parts, ignore_index=True)
    df.to_csv(CACHE / "sweep_nm_counts.csv", index=False)
    summary = df.groupby(["split", "arm", "q", "log2fc"])[["tp", "fp", "fn", "tn"]].sum()
    tp, fp, fn, tn = (summary[c] for c in ("tp", "fp", "fn", "tn"))
    summary = pd.DataFrame(dict(
        TP=tp, FP=fp, FN=fn, TN=tn,
        Precision=tp / (tp + fp), Recall=tp / (tp + fn), FDR=fp / (tp + fp),
        F1=2 * tp / (2 * tp + fp + fn), RejRate=(tp + fp) / (tp + fp + fn + tn),
        Prevalence=(tp + fn) / (tp + fp + fn + tn))).round(4)
    summary.to_csv(CACHE / "sweep_nm_summary.csv")
    print(summary.to_string())


if __name__ == "__main__":
    main()
