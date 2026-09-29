import sys
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

_TAU2_EPS = 1e-12


def simulate_marginal_nb(mu, alpha, tau2, n_reps, seed):
    # Posterior-predictive replicates matching marginal_rqr's mixture assumption:
    # b ~ N(0,tau2) per (rep, sample), mu_b = mu*exp(b), y_rep ~ NB(mu_b, alpha).
    # Draws the batch intercept fresh each rep so rep_var/rep_zero carry the same
    # across-sample mu spread the observed data has (a fixed-mu draw would not).
    mu = np.asarray(mu, dtype=np.float64)
    alpha = np.broadcast_to(np.asarray(alpha, dtype=np.float64), mu.shape)
    tau2 = np.broadcast_to(np.asarray(tau2, dtype=np.float64), mu.shape)
    rng = np.random.default_rng(seed)
    n_obs = len(mu)
    b = np.where(tau2 < _TAU2_EPS, 0.0, rng.normal(0.0, np.sqrt(np.maximum(tau2, 0.0)), size=(n_reps, n_obs)))
    mu_b = np.clip(mu[None, :] * np.exp(b), 1e-8, 1e10)
    n_ = np.broadcast_to(1.0 / np.maximum(alpha, 1e-8), (n_reps, n_obs))
    p_ = np.clip(n_ / (n_ + mu_b), 1e-10, 1 - 1e-10)
    return rng.negative_binomial(n_, p_)


def ppc_moment_pvalues(y, mu, alpha, tau2, n_reps=500, seed=0):
    # Bayesian posterior-predictive p = P(rep_stat >= obs_stat) for mean/var/zero/max.
    # p ~ U(0,1) under a correct model; p near 0/1 flags a moment the model cannot reproduce.
    y = np.asarray(y, dtype=np.float64)
    y_rep = simulate_marginal_nb(mu, alpha, tau2, n_reps, seed)
    return {
        "p_mean": float(np.mean(y_rep.mean(1) >= y.mean())),
        "p_var": float(np.mean(y_rep.var(1) >= y.var())),
        "p_zero": float(np.mean((y_rep == 0).mean(1) >= (y == 0).mean())),
        "p_max": float(np.mean(y_rep.max(1) >= y.max())),
    }


def calib_stats(o, p, logscale):
    o, p = np.asarray(o, float), np.asarray(p, float)
    m = np.isfinite(o) & np.isfinite(p)
    o, p = o[m], p[m]
    if logscale:
        o, p = np.log10(o + 1), np.log10(p + 1)
    r, _ = pearsonr(o, p)
    rho, _ = spearmanr(o, p)
    return dict(pearson_r=r, r2=r ** 2, spearman_rho=rho,
                rmse=float(np.sqrt(np.mean((p - o) ** 2))), mae=float(np.mean(np.abs(p - o))), n=len(o))


def binned_medians(o, p, nb=20):
    o, p = np.asarray(o, float), np.asarray(p, float)
    m = np.isfinite(o) & np.isfinite(p) & (o > 0)
    o, p = o[m], p[m]
    q = np.unique(np.quantile(o, np.linspace(0, 1, nb + 1)))
    idx = np.clip(np.digitize(o, q[1:-1]), 0, len(q) - 2)
    return (np.array([np.median(o[idx == k]) for k in range(len(q) - 1)]),
            np.array([np.median(p[idx == k]) for k in range(len(q) - 1)]))
