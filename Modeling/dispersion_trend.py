"""Covariate-free mean-dispersion trend for the Route B (Tier 2) engine.

Ignore covariates entirely. For every HC gene compute the closed-form NB2
method-of-moments dispersion from the raw counts:
    sigma_j = max(0, (var_j - mean_j) / mean_j^2)      (Var = mu + sigma*mu^2)
Genes with too few nonzero HC samples (nz < trend_min_nz) give a structurally
noisy MoM estimate (a single outlying observation can dominate the sample
variance) and are excluded from trend fitting -- not trimmed as outliers, since
the noise is not an outlier problem but an information-poverty problem.

The reliable genes are summarized by nonzero-weighted median in log(mu) bins,
then smoothed with lowess in log-log space (log(sigma) ~ log(mu)); a single
log-log line underfits both the low-mu plateau and the high-mu asymptote, so
lowess is the canonical trend, matching edgeR/DESeq2 trended-dispersion shape.

Route B modeling fixes each gene's dispersion at alpha_of(mean_train), so the
covariates spend all their degrees of freedom on the mean.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

import config

MP2 = config.MODELING_PARAMS_V2


def weighted_median(x, w):
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = np.cumsum(w)
    return float(x[np.searchsorted(c, 0.5 * c[-1])])


def build_trend(Y_hc, min_nz=None, n_bins=25, min_bin=20, lowess_frac=0.5):
    """Y_hc: (n_hc, n_genes) raw HC counts. Returns dict with lowess log-log curve."""
    min_nz = MP2["trend_min_nz"] if min_nz is None else min_nz
    nz = (Y_hc > 0).sum(0).astype(int)
    mean_c = Y_hc.mean(0)
    var_c = Y_hc.var(0)
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_mom = np.where(mean_c > 0, (var_c - mean_c) / mean_c ** 2, np.nan)
    sigma_mom = np.clip(sigma_mom, 0, None)

    reliable = (nz >= min_nz) & (mean_c > 0) & np.isfinite(sigma_mom)
    mu = mean_c[reliable]
    sig = sigma_mom[reliable]
    w = nz[reliable].astype(float)

    edges = np.geomspace(mu.min(), mu.max(), n_bins + 1)
    binid = np.clip(np.digitize(mu, edges) - 1, 0, n_bins - 1)
    rows = []
    for b in range(n_bins):
        m = binid == b
        if m.sum() < min_bin:
            continue
        rows.append({"mu_bin": weighted_median(mu[m], w[m]),
                     "sigma_wmed": weighted_median(sig[m], w[m]), "n": int(m.sum())})
    bins = pd.DataFrame(rows).sort_values("mu_bin")

    blx, bly = np.log(bins["mu_bin"].values), np.log(bins["sigma_wmed"].values)
    sm_curve = lowess(bly, blx, frac=lowess_frac, return_sorted=True)

    return {
        "a0": None, "a1": None,  # legacy parametric slot, unused (lowess is canonical)
        "alpha_floor": MP2["alpha_floor"], "alpha_cap": MP2["alpha_cap"],
        "min_nz": min_nz, "n_reliable": int(reliable.sum()), "n_bins_used": len(bins),
        "lowess_logmu": sm_curve[:, 0].tolist(), "lowess_logsigma": sm_curve[:, 1].tolist(),
    }


def save_trend(trend, path=None):
    path = Path(path or config.DISPERSION_TREND_PATH)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(trend, indent=2))


def load_trend(path=None):
    """Returns alpha_of(mean) -> fixed NB2 dispersion for Route B scoring/training."""
    path = Path(path or config.DISPERSION_TREND_PATH)
    cf = json.loads(path.read_text())
    logmu = np.asarray(cf["lowess_logmu"])
    logsig = np.asarray(cf["lowess_logsigma"])
    floor, cap = cf["alpha_floor"], cf["alpha_cap"]

    def alpha_of(mean):
        lm = np.log(max(float(mean), 1e-8))
        s = float(np.exp(np.interp(lm, logmu, logsig, left=logsig[0], right=logsig[-1])))
        return float(np.clip(s, floor, cap))

    return alpha_of
