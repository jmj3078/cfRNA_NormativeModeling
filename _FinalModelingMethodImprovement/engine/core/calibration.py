import numpy as np
from scipy.stats import kurtosis, norm, skew

from _FinalModelingMethodImprovement.engine.core.shash import fit_and_correct, shash_quantile


def bh_fdr_reject(pvals, q=0.05):
    p = np.asarray(pvals, dtype=np.float64)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    thresh = q * (np.arange(1, n + 1) / n)
    passed = ranked <= thresh
    reject = np.zeros(n, dtype=bool)
    if passed.any():
        k_max = np.nonzero(passed)[0].max()
        reject[order[:k_max + 1]] = True
    return reject


# naive_fdr_reject_rate/corr_fdr_reject_rate directly measure false-positive
# inflation under the naive N(0,1) assumption vs after SHASH warping (Fraza et
# al. 2021 NeuroImage; Efron 2007 Annals of Statistics). z_eval must be a true
# null (held-out HC, never disease) for these rates to mean anything.
def gene_shash_fit_eval(z_fit, z_eval):
    """Fit SHASH on z_fit (train/in-sample), evaluate calibration on the
    INDEPENDENT z_eval (held-out) -- params must never be fit on the same Z
    they are graded against, or the calibration numbers are circular."""
    z_eval = np.asarray(z_eval, dtype=np.float64)
    z_eval = z_eval[np.isfinite(z_eval)]
    params, z_corr = fit_and_correct(z_fit, z_eval)
    xi, eta, eps, delta, ok = params["xi"], params["eta"], params["eps"], params["delta"], params["ok"]
    z_lo, z_hi = shash_quantile(np.array([0.025, 0.975]), xi, eta, eps, delta)
    naive_exceed = float(np.mean(np.abs(z_eval) > 1.96))
    shash_exceed = float(np.mean((z_eval < z_lo) | (z_eval > z_hi)))
    p_naive = 2 * norm.sf(np.abs(z_eval))
    p_corr = 2 * norm.sf(np.abs(z_corr))
    naive_fdr = float(bh_fdr_reject(p_naive).mean())
    corr_fdr = float(bh_fdr_reject(p_corr).mean())
    return dict(
        shash_ok=ok, shash_xi=float(xi), shash_eta=float(eta), shash_eps=float(eps), shash_delta=float(delta),
        z_lo=float(z_lo), z_hi=float(z_hi),
        raw_skew=float(skew(z_eval)), raw_kurtosis=float(kurtosis(z_eval)),
        corrected_skew=float(skew(z_corr)), corrected_kurtosis=float(kurtosis(z_corr)),
        naive_exceed=naive_exceed, shash_exceed=shash_exceed,
        naive_fdr_reject_rate=naive_fdr, corr_fdr_reject_rate=corr_fdr,
    )


def calibration_metrics(z_raw, z_corr):
    """Naive-vs-SHASH-corrected calibration summary given an ALREADY pooled
    (raw, corrected) held-out pair -- for CV/LOBO, which fit SHASH per fold/
    batch on that fold's train Z and pool the resulting per-fold corrected
    held-out Z before reporting one number per gene."""
    z_raw = np.asarray(z_raw, dtype=np.float64)
    z_corr = np.asarray(z_corr, dtype=np.float64)
    p_naive = 2 * norm.sf(np.abs(z_raw))
    p_corr = 2 * norm.sf(np.abs(z_corr))
    return dict(
        raw_skew=float(skew(z_raw)), raw_kurtosis=float(kurtosis(z_raw)),
        corrected_skew=float(skew(z_corr)), corrected_kurtosis=float(kurtosis(z_corr)),
        naive_exceed=float(np.mean(np.abs(z_raw) > 1.96)), shash_exceed=float(np.mean(np.abs(z_corr) > 1.96)),
        naive_fdr_reject_rate=float(bh_fdr_reject(p_naive).mean()),
        corr_fdr_reject_rate=float(bh_fdr_reject(p_corr).mean()),
    )


def gene_shash_calibration(z):
    """In-sample fit+self-evaluation -- for the production engine's own SHASH
    fit, where there is no independent held-out split (CV/LOBO use
    gene_shash_fit_eval with a genuine train/held-out pair instead)."""
    return gene_shash_fit_eval(z, z)
