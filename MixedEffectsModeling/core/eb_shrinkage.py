"""Empirical-Bayes shrinkage of the NB2 dispersion submodel (limma/edgeR style).

Both targets use the same moment decomposition: a per-gene MLE phi_hat_g carries
its own estimation error, so Var(phi_hat) = tau^2 + mean(SE^2) and the prior
scale is recovered by subtracting the error component. MAD/median replace
variance/mean because a handful of near-divergent genes would otherwise inflate
tau and silently disable the shrinkage.

  - estimate_slope_prior: tau_k for the dispersion SLOPES, from a no-prior calibration
    run. Fed back into glmmTMB as normal(0, tau_k) per covariate.
  - squeeze_log_theta: precision-weighted posterior mean of the dispersion
    INTERCEPT toward the Phase-0 lowess trend. SE=NaN (unusable sdreport) maps to
    SE^2=inf, i.e. exactly the trend value -- v2's hard-fixed nb_fixed stage is
    the limiting case of this rule rather than a separate stage.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

import MixedEffectsModeling.config as config

EB = config.EB_PARAMS


def robust_var(x):
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if len(x) < 8:
        return np.nan
    return (1.4826 * np.median(np.abs(x - np.median(x)))) ** 2


def _tau2(est_spread_var, se, floor):
    err = np.nanmedian(np.asarray(se, dtype=np.float64) ** 2)
    err = 0.0 if not np.isfinite(err) else err
    base = 0.0 if not np.isfinite(est_spread_var) else est_spread_var
    return max(base - err, floor ** 2)


def estimate_slope_prior(calib_csv, n_cov=None, floor=None):
    """tau_k per dispersion slope from a --mode calib run (no dispersion prior)."""
    n_cov = len(config.COVARIATES) if n_cov is None else n_cov
    floor = EB["tau_floor"] if floor is None else floor
    df = pd.read_csv(calib_csv)
    df = df[df["ok"].astype(bool) & (df["stage"] == "nbi_full_eb")]
    tau = []
    for k in range(1, n_cov + 1):
        est = df[f"disp_coef_{k}"].to_numpy(dtype=np.float64)
        se = df[f"disp_se_{k}"].to_numpy(dtype=np.float64)
        tau.append(float(np.sqrt(_tau2(robust_var(est), se, floor))))
    return {"n_calib_genes": int(len(df)), "covariates": list(config.COVARIATES), "tau_slope": tau}


def save_disp_prior(prior, path=None):
    path = Path(path or config.DISP_PRIOR_PATH)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(prior, indent=2))


def load_disp_prior(path=None):
    path = Path(path or config.DISP_PRIOR_PATH)
    return json.loads(path.read_text()) if path.exists() else None


def squeeze_log_theta(log_theta_hat, se, log_theta_trend, floor=None):
    """Precision-weighted posterior mean of log(theta) toward the lowess trend.

    Returns (log_theta_post, tau_d2). tau_d2 is estimated once from the pooled
    residual spread of all supplied genes, so nbi_full_eb and nbi_intercept_eb must be passed
    together (nbi_intercept_eb alone is too small for a stable estimate).
    """
    floor = EB["tau_floor"] if floor is None else floor
    hat = np.asarray(log_theta_hat, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    trend = np.asarray(log_theta_trend, dtype=np.float64)
    tau_d2 = _tau2(robust_var(hat - trend), se, floor)
    se2 = np.where(np.isfinite(se) & (se > 0), se ** 2, np.inf)
    with np.errstate(divide="ignore", invalid="ignore"):
        prec = 1.0 / se2
        post = (np.nan_to_num(hat * prec) + trend / tau_d2) / (prec + 1.0 / tau_d2)
    return np.where(np.isfinite(hat), post, trend), float(tau_d2)
