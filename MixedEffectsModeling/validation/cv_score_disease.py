"""Score disease samples with the cached CV fold models (fold_params/), so disease Z come from the same
80%-HC fold models as the out-of-fold HC Z in cv_zscores_shash.pkl. No refit: reuses model_fold*/pool_fold*
fits and their per-fold SHASH params. Each disease sample is scored by one fold (KFold, seed 42).

Run from the repo root: NM_PIPELINE=p0 python -m MixedEffectsModeling.validation.cv_score_disease
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.covariates import Covariates
from MixedEffectsModeling.core.dispersion_trend import load_trend
from MixedEffectsModeling.core.marginal_rqr import _poisson_rqr, marginal_nb_rqr
from MixedEffectsModeling.core.model_engine_mixed import NormativeModelEngineMixed
from MixedEffectsModeling.core.shash import shash_transform_to_z
from MixedEffectsModeling.validation.cv_engine import MP, POISSON_ALPHA_EPS, _mu_alpha
from MixedEffectsModeling.validation.lobo_engine import load_full_data


def shash(z, srow):
    return shash_transform_to_z(z, srow["xi"], srow["eta"], srow["eps"], srow["delta"]) if bool(srow["ok"]) else z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine-dir", type=Path, default=config.ENGINE_MIXED_DIR)
    ap.add_argument("--cv-dir", type=Path, default=config.CV_MIXED_DIR)
    args = ap.parse_args()
    fp = args.cv_dir / "fold_params"
    out = args.cv_dir / "Z_disease_cv_shash.npy"

    e2 = NormativeModelEngineMixed()
    e2.load_hc_data()
    e2.alpha_fn = load_trend(args.engine_dir / "dispersion_trend.json")
    folds = list(StratifiedKFold(MP["n_splits"], shuffle=True, random_state=42).split(np.zeros(len(e2.names)), e2.batch))
    saved = json.loads((config.OUTRIDER_COMPARISON_DIR / "cv_folds.json").read_text())
    assert all(e2.names[te].tolist() == saved[str(i)]["test"] for i, (_, te) in enumerate(folds))

    data = load_full_data()
    dis = np.where(~data["is_hc"])[0]
    meta = pd.read_csv(config.ZSCORES_MIXED_DIR / "sample_meta.csv")
    assert (meta["sample"].values == data["names"][dis]).all()
    genes = pd.read_pickle(config.ZSCORES_MIXED_DIR / "gene_names.pkl")
    col = {g: j for j, g in enumerate(genes)}
    obs, Y, gcol = data["obs"][dis], data["Y"][dis], data["gene_col"]

    dfold = np.zeros(len(dis), dtype=int)
    for fi, (_, te) in enumerate(KFold(MP["n_splits"], shuffle=True, random_state=42).split(dis)):
        dfold[te] = fi
    Z = np.full((len(dis), len(genes)), np.nan, dtype=np.float32)

    for fi, (tr, _) in enumerate(folds):
        t0 = time.time()
        rows = np.where(dfold == fi)[0]
        cov = Covariates(e2._gene_col).fit(e2.Y_hc[tr])
        scaler = StandardScaler().fit(cov.transform(e2.obs_hc, e2.Y_hc)[tr])
        Xd = scaler.transform(cov.transform(obs[rows], Y[rows]))
        Xa = np.column_stack([np.ones(len(rows)), Xd])

        fits = pd.read_csv(fp / f"model_fold{fi}.csv").set_index("gene")
        sh = pd.read_csv(fp / f"shash_model_fold{fi}.csv").set_index("gene")
        mu_cols = [c for c in fits.columns if c.startswith("mu_coef_")]
        dp_cols = [c for c in fits.columns if c.startswith("disp_coef_")]
        for g, row in fits[fits["ok"].astype(bool)].iterrows():
            if g not in col or g not in sh.index:
                continue
            mu, alpha = _mu_alpha(Xa, row[mu_cols].values.astype(float), row[dp_cols].values.astype(float), row, e2.alpha_fn)
            z = marginal_nb_rqr(Y[rows, gcol[g]], mu, alpha, float(row["tau2"]), seed=2042 + fi)
            Z[rows, col[g]] = shash(z, sh.loc[g])

        pool = json.loads((fp / f"pool_fold{fi}.json").read_text())
        if pool["ok"]:
            shp = pd.read_csv(fp / f"shash_pool_fold{fi}.csv").set_index("gene")
            beta = np.asarray(pool["beta"])
            tau2 = float(pool["tau2"]) if pool.get("tau2") is not None else 0.0
            alpha_eff = pool["alpha"] if pool["family"] == "negbin" else POISSON_ALPHA_EPS
            mult = np.exp(Xd @ beta[1:])
            if pool.get("mult_lo") is not None:
                mult = np.clip(mult, pool["mult_lo"], pool["mult_hi"])
            for g, m in zip(pool["gene"], pool["mean_hc"]):
                if g not in col or g not in shp.index:
                    continue
                mu = np.clip((m + pool["eps"]) * np.exp(beta[0]) * mult, 1e-6, 1e8)
                y = Y[rows, gcol[g]].astype(np.float64)
                z = _poisson_rqr(y, mu, 2042 + fi) if pool["family"] == "poisson" else marginal_nb_rqr(y, mu, alpha_eff, tau2, seed=2042 + fi)
                Z[rows, col[g]] = shash(z, shp.loc[g])
        print(f"fold {fi}: {len(rows)} disease samples, {time.time() - t0:.0f}s", flush=True)

    np.save(out, Z)
    pd.DataFrame({"sample": meta["sample"], "cv_fold": dfold}).to_csv(args.cv_dir / "disease_cv_fold.csv", index=False)
    print(f"Saved -> {out} {Z.shape}, finite frac {np.isfinite(Z).mean():.4f}")


if __name__ == "__main__":
    main()
