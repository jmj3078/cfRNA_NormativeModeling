"""OutSingle (Salkovic et al., Bioinformatics 2023, doi 10.1093/bioinformatics/btad142)
reimplemented in the frozen/inductive form this benchmark needs.

Log-normal count model + SVD, with the latent dimension chosen by the Gavish-Donoho
optimal hard threshold rather than a tuned hyperparameter. The published implementation
fits and scores one cohort at a time; here the train basis is frozen so held-out samples
project through it without refitting, matching frozen_outrider.R::fit_train/score_frozen.

Correctness is checked two ways: gd_rank against the closed-form threshold on a matrix of
known rank, and the frozen path against the in-sample path on the training samples
themselves (the same contract validate_frozen.R enforces for OUTRIDER).
"""
import numpy as np
from scipy.stats import norm

from common import CACHE, size_factors


def gd_rank(sv, n_rows, n_cols):
    """Gavish & Donoho (2014) optimal hard threshold for singular values, unknown noise
    level: tau = omega(beta) * median(sv), with the published polynomial approximation."""
    beta = min(n_rows, n_cols) / max(n_rows, n_cols)
    omega = 0.56 * beta ** 3 - 0.95 * beta ** 2 + 1.82 * beta + 1.43
    return max(int((sv > omega * np.median(sv)).sum()), 1)


def _standardize(counts, loggeomeans, center, scale):
    sf = size_factors(counts, loggeomeans)
    lg = np.log2(counts / sf[:, None] + 1)
    return (lg - center) / scale, sf


def fit_train(y_train, tag=None):
    """y_train: samples x genes of raw counts. Returns the frozen reference.
    With a tag, the basis is cached so reruns and later analyses reuse the same one."""
    path = (CACHE / "outsingle_fits" / f"{tag}.npz") if tag else None
    if path is not None and path.exists():
        return dict(np.load(path))

    with np.errstate(divide="ignore"):
        loggeomeans = np.log(y_train).mean(axis=0)
    sf = size_factors(y_train, loggeomeans)
    lg = np.log2(y_train / sf[:, None] + 1)
    center = lg.mean(axis=0)
    scale = np.where(lg.std(axis=0) > 1e-8, lg.std(axis=0), 1.0)
    x = (lg - center) / scale

    u, sv, vt = np.linalg.svd(x, full_matrices=False)
    k = gd_rank(sv, *x.shape)
    vk = vt[:k]
    resid = x - (x @ vk.T) @ vk
    resid_sd = np.where(resid.std(axis=0) > 1e-8, resid.std(axis=0), 1.0)
    fit = dict(loggeomeans=loggeomeans, center=center, scale=scale, vk=vk, k=k,
               resid_sd=resid_sd)
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, **fit)
    return fit


def score_frozen(fit, y_test):
    """Project unseen samples through the frozen basis. No test data enters any
    fitted quantity -- basis, centering, scaling, and residual sd all come from train."""
    x, sf = _standardize(y_test, fit["loggeomeans"], fit["center"], fit["scale"])
    resid = x - (x @ fit["vk"].T) @ fit["vk"]
    z = resid / fit["resid_sd"]
    return dict(z=z, pval=2 * norm.sf(np.abs(z)), sf=sf)


def _demo():
    rng = np.random.default_rng(0)
    n_train, n_test, n_genes, true_k = 200, 40, 400, 4

    # counts with a shared low-rank confounder structure plus per-gene depth
    loadings = rng.normal(size=(true_k, n_genes)) * 1.5
    base = rng.uniform(3, 8, size=n_genes)

    def make(n):
        scores = rng.normal(size=(n, true_k))
        lg = base + scores @ loadings + rng.normal(scale=0.3, size=(n, n_genes))
        return np.round(np.maximum(2 ** lg - 1, 0))

    y_tr, y_te = make(n_train), make(n_test)

    # 1. gd_rank recovers a planted rank on a clean low-rank-plus-noise matrix
    signal = rng.normal(size=(n_train, true_k)) @ rng.normal(size=(true_k, n_genes)) * 20
    noisy = signal + rng.normal(size=(n_train, n_genes))
    assert gd_rank(np.linalg.svd(noisy, compute_uv=False), *noisy.shape) == true_k

    fit = fit_train(y_tr)
    assert 1 <= fit["k"] <= min(y_tr.shape), fit["k"]

    # 2. frozen scoring of the TRAIN samples reproduces the in-sample residuals exactly
    x_in = (np.log2(y_tr / size_factors(y_tr, fit["loggeomeans"])[:, None] + 1)
            - fit["center"]) / fit["scale"]
    z_in = (x_in - (x_in @ fit["vk"].T) @ fit["vk"]) / fit["resid_sd"]
    assert np.abs(score_frozen(fit, y_tr)["z"] - z_in).max() < 1e-10

    # 3. held-out z is standardized and an injected spike is detected
    out = score_frozen(fit, y_te)
    assert out["z"].shape == (n_test, n_genes)
    assert abs(out["z"].std() - 1) < 0.35, out["z"].std()

    y_spiked = y_te.copy()
    y_spiked[0, :20] = np.round(y_spiked[0, :20] * 32 + 50)
    spiked = score_frozen(fit, y_spiked)["z"][0, :20]
    assert np.abs(spiked).mean() > np.abs(out["z"][0, :20]).mean(), "spike not detected"

    print(f"outsingle OK | gd_rank={fit['k']} | held-out z sd={out['z'].std():.3f} "
          f"| spiked |z| mean={np.abs(spiked).mean():.2f}")


if __name__ == "__main__":
    _demo()
