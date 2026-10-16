import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.stats import nbinom, norm, poisson

RQR_EPS = 1e-8
RQR_EPS_DIRECT = 1e-15  # tau2-collapsed genes skip GH marginalization -> nbinom.cdf tail is precise enough for a tighter bound


def _poisson_rqr(y, mu, seed=None):
    y = np.asarray(y)
    lo = np.where(y > 0, poisson.cdf(y - 1, mu), 0.0)
    hi = poisson.cdf(y, mu)
    lo = np.clip(lo, RQR_EPS, 1 - RQR_EPS); hi = np.clip(hi, RQR_EPS, 1 - RQR_EPS)
    rng = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(np.minimum(lo, hi), np.maximum(lo, hi))).astype(np.float32)


def _nb_cdf(y, mu, alpha):
    n = 1.0 / alpha
    p = np.clip(n / (n + mu), RQR_EPS, 1 - RQR_EPS)
    return nbinom.cdf(y, n, p)


def _nb_rqr(y, mu, alpha, seed=None):
    y = np.asarray(y)
    n = 1.0 / alpha
    p = np.clip(n / (n + mu), RQR_EPS, 1 - RQR_EPS)
    lo = np.where(y > 0, nbinom.cdf(y - 1, n, p), 0.0)
    hi = nbinom.cdf(y, n, p)
    lo = np.clip(lo, RQR_EPS_DIRECT, 1 - RQR_EPS_DIRECT); hi = np.clip(hi, RQR_EPS_DIRECT, 1 - RQR_EPS_DIRECT)
    rng = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(np.minimum(lo, hi), np.maximum(lo, hi))).astype(np.float32)


def _nb_logpmf(y, mu, alpha):
    n = 1.0 / alpha
    p = np.clip(n / (n + mu), RQR_EPS, 1 - RQR_EPS)
    return nbinom.logpmf(y, n, p)


def marginal_nb_loglik(y, mu, alpha, tau2, n_nodes=7):
    y = np.asarray(y)
    tau2 = np.asarray(tau2, dtype=np.float64)
    if np.all(tau2 < 1e-6):
        return _nb_logpmf(y, mu, alpha)
    nodes, weights = hermegauss(n_nodes)
    weights = weights / weights.sum()
    sd = np.sqrt(np.maximum(tau2, 0.0))
    logpmf_k = np.stack([
        _nb_logpmf(y, mu * np.exp(sd * node), alpha) + np.log(w)
        for node, w in zip(nodes, weights)
    ])
    m = logpmf_k.max(axis=0)
    return m + np.log(np.exp(logpmf_k - m).sum(axis=0))


def nb_marginal_mean_var(mu, alpha, tau2):
    tau2 = np.asarray(tau2, dtype=np.float64)
    if np.all(tau2 < 1e-6):
        return mu, mu + alpha * mu ** 2
    mean = mu * np.exp(tau2 / 2)
    ey2 = mu * np.exp(tau2 / 2) + (1 + alpha) * mu ** 2 * np.exp(2 * tau2)
    return mean, ey2 - mean ** 2


def nb_marginal_pmf0(mu, alpha, tau2, n_nodes=7):
    tau2 = np.asarray(tau2, dtype=np.float64)
    n_ = 1.0 / np.maximum(alpha, 1e-8)
    if np.all(tau2 < 1e-6):
        p_ = np.clip(n_ / (n_ + mu), RQR_EPS, 1 - RQR_EPS)
        return nbinom.pmf(0, n_, p_)
    nodes, weights = hermegauss(n_nodes)
    weights = weights / weights.sum()
    sd = np.sqrt(np.maximum(tau2, 0.0))
    total = np.zeros_like(mu)
    for node, w in zip(nodes, weights):
        mu_b = mu * np.exp(sd * node)
        p_ = np.clip(n_ / (n_ + mu_b), RQR_EPS, 1 - RQR_EPS)
        total += w * nbinom.pmf(0, n_, p_)
    return total


def marginal_nb_rqr(y, mu, alpha, tau2, seed, n_nodes=7):
    y = np.asarray(y)
    tau2 = np.asarray(tau2, dtype=np.float64)
    if np.all(tau2 < 1e-6):
        return _nb_rqr(y, mu, alpha, seed)

    nodes, weights = hermegauss(n_nodes)
    weights = weights / weights.sum()
    sd = np.sqrt(np.maximum(tau2, 0.0))
    lo = np.zeros_like(y, dtype=np.float64)
    hi = np.zeros_like(y, dtype=np.float64)
    for node, w in zip(nodes, weights):
        mu_b = mu * np.exp(sd * node)
        lo += w * np.where(y > 0, _nb_cdf(y - 1, mu_b, alpha), 0.0)
        hi += w * _nb_cdf(y, mu_b, alpha)
    lo = np.clip(lo, RQR_EPS, 1 - RQR_EPS)
    hi = np.clip(hi, RQR_EPS, 1 - RQR_EPS)
    rng = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(np.minimum(lo, hi), np.maximum(lo, hi))).astype(np.float32)
