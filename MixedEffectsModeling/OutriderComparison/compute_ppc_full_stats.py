import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from MixedEffectsModeling.validation.ppc_simulate import ppc_moment_pvalues

DIR = Path(__file__).parent
N_REPS = 200


def calib_stats(o, p, logscale):
    o, p = np.asarray(o, float), np.asarray(p, float)
    m = np.isfinite(o) & np.isfinite(p)
    o, p = o[m], p[m]
    if logscale:
        o, p = np.log10(o + 1), np.log10(p + 1)
    r, _ = pearsonr(o, p)
    rho, _ = spearmanr(o, p)
    rmse = np.sqrt(np.mean((p - o) ** 2))
    mae = np.mean(np.abs(p - o))
    return dict(pearson_r=r, r2=r ** 2, spearman_rho=rho, rmse=rmse, mae=mae, n=len(o))


def gene_level_calib_table(cal_csv):
    cal = pd.read_csv(cal_csv)
    panels = [("obs_mean", "pred_mean", "mean", True), ("obs_var", "pred_var", "var", True),
             ("obs_zero", "pred_zero", "zero_frac", False)]
    return {lab: calib_stats(cal[oc], cal[pc], logsc) for oc, pc, lab, logsc in panels}


def pvalue_summary(y_mu_alpha_tau2_iter, n_genes, seed0=2000):
    rows = []
    for i, (g, y, mu, alpha, tau2) in enumerate(y_mu_alpha_tau2_iter):
        if y.size == 0 or not (np.isfinite(mu).all() and np.isfinite(alpha).all() and np.isfinite(tau2).all()):
            continue
        pv = ppc_moment_pvalues(y, mu, alpha, tau2, n_reps=N_REPS, seed=seed0 + i)
        pv["gene"] = g
        rows.append(pv)
        if (i + 1) % 2000 == 0:
            print(f"  {i+1}/{n_genes}")
    return pd.DataFrame(rows)


def our_engine_iter():
    d = pickle.load(open(DIR.parent.parent / "CV_Results_mixed" / "cv_ppc.pkl", "rb"))
    for g, v in d.items():
        yield g, np.asarray(v["y"], float), np.asarray(v["mu"], float), np.asarray(v["alpha"], float), np.asarray(v["tau2"], float)


def outrider_iter():
    mus, ys, thetas = [], [], []
    for fi in range(5):
        mus.append(pd.read_csv(DIR / f"cv_fold{fi}_mu.csv", index_col=0))
        ys.append(pd.read_csv(DIR / f"cv_fold{fi}_y.csv", index_col=0))
        thetas.append(pd.read_csv(DIR / f"cv_fold{fi}_theta.csv").set_index("gene")["theta"])
    common_genes = set(mus[0].columns)
    for d in mus[1:]:
        common_genes &= set(d.columns)
    common_genes = sorted(common_genes)
    Y = pd.concat([d[common_genes] for d in ys], axis=0)
    MU = pd.concat([d[common_genes] for d in mus], axis=0)
    theta_rows = [pd.DataFrame(np.tile(thetas[fi][common_genes].values, (len(ys[fi]), 1)), columns=common_genes)
                 for fi in range(5)]
    TH = pd.concat(theta_rows, axis=0)
    TH.index = Y.index
    for g in common_genes:
        y = Y[g].values.astype(float)
        mu = MU[g].values.astype(float)
        alpha = 1.0 / np.maximum(TH[g].values.astype(float), 1e-8)
        tau2 = np.zeros_like(mu)
        yield g, y, mu, alpha, tau2


if __name__ == "__main__":
    out_eng = DIR / "our_engine_ppc_pvalues.csv"
    out_outr = DIR / "outrider_ppc_pvalues.csv"

    if not out_eng.exists():
        print("computing our-engine PPC p-values...")
        df = pvalue_summary(our_engine_iter(), 19858)
        df.to_csv(out_eng, index=False)
        print(f"saved -> {out_eng}")
    else:
        print(f"already cached -> {out_eng}")

    if not out_outr.exists():
        print("computing OUTRIDER PPC p-values...")
        df = pvalue_summary(outrider_iter(), 12305)
        df.to_csv(out_outr, index=False)
        print(f"saved -> {out_outr}")
    else:
        print(f"already cached -> {out_outr}")

    eng_calib = gene_level_calib_table(DIR.parent.parent / "CV_Results_mixed" / "cv_calibration_moments.csv")
    outr_calib = gene_level_calib_table(DIR / "outrider_cv_calibration_moments.csv")
    print("\nour engine calib_stats (full 19858g):", eng_calib)
    print("\nOUTRIDER calib_stats (12305g):", outr_calib)

    rows = []
    for engine_name, table in [("our_engine_full19858", eng_calib), ("outrider_12305", outr_calib)]:
        for panel, stats in table.items():
            rows.append(dict(engine=engine_name, panel=panel, **stats))
    pd.DataFrame(rows).to_csv(DIR / "ppc_calib_stats.csv", index=False)
    print(f"\nsaved -> {DIR / 'ppc_calib_stats.csv'}")

    # our-engine calib_stats restricted to the same 12305-gene subset OUTRIDER covers, for a
    # like-for-like comparison alongside the full-19858 number.
    eng_cal_full = pd.read_csv(DIR.parent.parent / "CV_Results_mixed" / "cv_calibration_moments.csv")
    outr_genes = pd.read_csv(DIR / "outrider_cv_calibration_moments.csv")["gene"]
    eng_cal_sub = eng_cal_full[eng_cal_full["gene"].isin(outr_genes)]
    sub_path = DIR / "our_engine_cv_calibration_moments_12305subset.csv"
    eng_cal_sub.to_csv(sub_path, index=False)
    panels = [("obs_mean", "pred_mean", "mean", True), ("obs_var", "pred_var", "var", True),
             ("obs_zero", "pred_zero", "zero_frac", False)]
    eng_calib_sub = {lab: calib_stats(eng_cal_sub[oc], eng_cal_sub[pc], logsc) for oc, pc, lab, logsc in panels}
    rows_sub = [dict(engine="our_engine_12305subset", panel=panel, **stats) for panel, stats in eng_calib_sub.items()]
    pd.concat([pd.read_csv(DIR / "ppc_calib_stats.csv"), pd.DataFrame(rows_sub)], ignore_index=True).to_csv(
        DIR / "ppc_calib_stats.csv", index=False)
    print("our engine calib_stats (12305g subset):", eng_calib_sub)
