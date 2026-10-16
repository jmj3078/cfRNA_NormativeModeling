"""Always-on diagnostic report for the two calibration-derived hyperparameters:
the covariate-adjusted dispersion trend and the per-covariate EB slope prior sd.

Written by prepare_hyperparams on every run that executes a calibration fit, so
a trend is never deployed without its calibration record. The residual table it saves is the
acceptance criterion: log(alpha_fit) - log(alpha_trend) must be centred at 0 in
every expression bin, which the covariate-free MoM trend fails by up to 2.8 log
units at high expression.
"""

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from viz_style import apply_style

COL_FIT = "#2a9d8f"
COL_REF = "#d1495b"
COL_PT = "#9aa5b1"


def curve_of(trend):
    lm = np.asarray(trend["lowess_logmu"], dtype=np.float64)
    ls = np.asarray(trend["lowess_logsigma"], dtype=np.float64)
    return lm, ls


def eval_trend(trend, mean_hc):
    lm, ls = curve_of(trend)
    x = np.log(np.maximum(np.asarray(mean_hc, dtype=np.float64), 1e-8))
    return np.clip(np.exp(np.interp(x, lm, ls, left=ls[0], right=ls[-1])),
                   trend["alpha_floor"], trend["alpha_cap"])


def residual_table(mean_hc, alpha_fit, trend, n_bins=12):
    d = pd.DataFrame({"mean": mean_hc, "alpha_fit": alpha_fit})
    d = d[(d["mean"] > 0) & np.isfinite(d.alpha_fit) & (d.alpha_fit > 0)]
    d["alpha_trend"] = eval_trend(trend, d["mean"].values)
    d["resid"] = np.log(d.alpha_fit) - np.log(d.alpha_trend)
    d["bin"] = pd.qcut(np.log(d["mean"]), n_bins, duplicates="drop")
    t = d.groupby("bin", observed=True).agg(
        n=("mean", "size"), mu=("mean", "median"),
        alpha_fit=("alpha_fit", "median"), alpha_trend=("alpha_trend", "median"),
        resid_med=("resid", "median"),
        resid_q25=("resid", lambda s: s.quantile(0.25)),
        resid_q75=("resid", lambda s: s.quantile(0.75)),
        frac_pos=("resid", lambda s: (s > 0).mean())).reset_index(drop=True)
    return d, t


def trend_report(mean_hc, alpha_fit, trend, disp_prior, out_dir, mom_trend=None, tau_d2=None):
    """Saves Figures/dispersion_trend.png + trend_residuals.csv, returns summary stats."""
    apply_style()
    out_dir = Path(out_dir)
    fig_dir = out_dir / "Figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    d, t = residual_table(mean_hc, alpha_fit, trend)
    t.to_csv(out_dir / "trend_residuals.csv", index=False)
    lm, ls = curve_of(trend)
    rv = (1.4826 * np.median(np.abs(d.resid - np.median(d.resid)))) ** 2

    fig, ax = plt.subplots(1, 4, figsize=(20, 4.3))

    ax[0].scatter(d["mean"], d.alpha_fit, s=2, alpha=0.15, color=COL_PT, rasterized=True,
                  label=f"calib fits (n={len(d)})")
    if mom_trend is not None:
        mlm, mls = curve_of(mom_trend)
        ax[0].plot(np.exp(mlm), np.exp(mls), color=COL_REF, lw=2, ls="--",
                   label="covariate-free MoM (reference)")
    ax[0].plot(np.exp(lm), np.exp(ls), color=COL_FIT, lw=2.2,
               label=f"trend: lowess frac={trend.get('lowess_frac')} it={trend.get('lowess_it')}")
    ax[0].set(xscale="log", yscale="log", xlabel="HC mean count", ylabel="dispersion alpha",
              title="A  Conditional dispersion trend")
    ax[0].legend(frameon=False, fontsize=7.5, loc="lower left")

    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].fill_between(t.mu, t.resid_q25, t.resid_q75, color=COL_FIT, alpha=0.2, lw=0)
    ax[1].plot(t.mu, t.resid_med, "s-", color=COL_FIT, label="bin median")
    if mom_trend is not None:
        mr = np.log(d.alpha_fit) - np.log(eval_trend(mom_trend, d["mean"].values))
        mt = mr.groupby(d["bin"], observed=True).median()
        ax[1].plot(t.mu, mt.values, "o-", color=COL_REF, label="covariate-free MoM")
    ax[1].set(xscale="log", xlabel="HC mean count", ylabel="log(alpha_fit) - log(alpha_trend)",
              title=f"B  Bias vs expression (max |med| = {t.resid_med.abs().max():.3f})")
    ax[1].legend(frameon=False, fontsize=7.5, loc="lower left")

    ax[2].hist(d.resid, bins=80, color=COL_FIT, alpha=0.75)
    ax[2].axvline(0, color="k", lw=0.8)
    ax[2].axvline(d.resid.median(), color=COL_REF, lw=1.5,
                  label=f"median={d.resid.median():+.3f}")
    td2 = rv if tau_d2 is None else tau_d2
    ax[2].set(xlabel="log(alpha_fit) - log(alpha_trend)", ylabel="genes",
              title=f"C  Residual spread (tau_d={np.sqrt(td2):.3f})")
    ax[2].legend(frameon=False, fontsize=7.5)

    tau = np.asarray(disp_prior["tau_slope"], dtype=np.float64)
    y = np.arange(len(tau))
    ax[3].barh(y, tau, color=COL_FIT, alpha=0.85)
    ax[3].axvline(0.05, color=COL_REF, lw=1.5, ls="--", label="v2 blanket prior sd")
    ax[3].set_yticks(y)
    ax[3].set_yticklabels(disp_prior["covariates"], fontsize=7)
    ax[3].invert_yaxis()
    ax[3].set(xlabel="EB prior sd tau_k", title="D  Dispersion slope prior")
    ax[3].legend(frameon=False, fontsize=7.5, loc="lower right")

    fig.tight_layout()
    fig.savefig(fig_dir / "dispersion_trend.png", dpi=170)
    plt.close(fig)
    return {"n_genes": int(len(d)), "resid_median": float(d.resid.median()),
            "resid_robust_sd": float(np.sqrt(rv)),
            "max_abs_bin_median": float(t.resid_med.abs().max()),
            "figure": str(fig_dir / "dispersion_trend.png")}
