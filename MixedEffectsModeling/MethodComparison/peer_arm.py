"""PEER arm in frozen/normative form, following the GTEx expression-outlier construction
(Li et al. 2017, Nature, doi 10.1038/nature24267): correct the expression matrix for
inferred hidden factors, then score each individual as a Z against the reference.

PEER itself only fits a cohort -- PEER_getX returns factor values for the samples it was
given and there is no API for a new sample. The inductive step is the posterior mean of a
held-out sample's factors under the FROZEN weights, which for y = Wx + e, e ~ N(0, diag(s2)),
x ~ N(0, I) is

    x = (I + W' L W)^-1 W' L y,     L = diag(1/s2)

i.e. a ridge projection onto the frozen basis. Residual = y - Wx, standardized by the
train residual sd. Structurally identical to the OutSingle arm except the basis comes from
ARD factor analysis rather than SVD + optimal hard threshold.

The PEER fit runs in peer_env via fit_peer.R (R 3.4.1); everything else is numpy here.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from common import CACHE, size_factors

HERE = Path(__file__).resolve().parent
PEER_DIR = CACHE / "peer_fits"
CONDA_SH = "/home/mjcho/miniconda3/etc/profile.d/conda.sh"
# Matched to the autoencoder and PCA arms' q, so the latent dimension is not a free
# variable between arms. ARD prunes below this anyway. (K=60, the GTEx rule-of-thumb for
# n~540, was measured at >75 min/fold single-threaded and would also have made PEER the
# only arm with a different latent dimension.)
K_FACTORS = 20
MAX_ITER = 100


def _standardize_train(y_train):
    with np.errstate(divide="ignore"):
        loggeomeans = np.log(y_train).mean(axis=0)
    sf = size_factors(y_train, loggeomeans)
    lg = np.log2(y_train / sf[:, None] + 1)
    center = lg.mean(axis=0)
    scale = np.where(lg.std(axis=0) > 1e-8, lg.std(axis=0), 1.0)
    return (lg - center) / scale, dict(loggeomeans=loggeomeans, center=center, scale=scale)


def fit_train(y_train, tag, k=K_FACTORS, max_iter=MAX_ITER):
    """Returns the frozen PEER reference. Caches the R fit, so a rerun is free."""
    PEER_DIR.mkdir(parents=True, exist_ok=True)
    prefix = PEER_DIR / tag
    x_train, tr = _standardize_train(y_train)

    if not (PEER_DIR / f"{tag}_W.csv").exists():
        train_csv = PEER_DIR / f"{tag}_train.csv"
        pd.DataFrame(x_train).to_csv(train_csv)
        # stream to a log rather than capturing: PEER prints per-iteration progress and
        # these fits run for tens of minutes, so silent capture hides whether it is moving
        log = PEER_DIR / f"{tag}_fit.out"
        cmd = (f"source {CONDA_SH} && conda activate peer_env && "
               f"Rscript {HERE / 'fit_peer.R'} {train_csv} {prefix} {k} {max_iter}")
        with open(log, "w") as fh:
            subprocess.run(["bash", "-lc", cmd], stdout=fh, stderr=subprocess.STDOUT, text=True)
        tail = log.read_text()
        if "PEER_FIT_DONE" not in tail:
            raise RuntimeError(f"PEER fit failed for {tag}, see {log}:\n{tail[-2000:]}")
        print([l for l in tail.splitlines() if l.startswith("W:")][-1], flush=True)
        train_csv.unlink()

    w = pd.read_csv(f"{prefix}_W.csv").to_numpy(dtype=np.float64)
    eps = np.clip(pd.read_csv(f"{prefix}_eps.csv")["eps"].to_numpy(dtype=np.float64), 1e-8, None)

    fit = dict(w=w, eps=eps, **tr)
    resid = x_train - _project(x_train, fit)
    fit["resid_sd"] = np.where(resid.std(axis=0) > 1e-8, resid.std(axis=0), 1.0)
    return fit


def _project(x, fit):
    """Posterior-mean reconstruction under the frozen basis. eps is PEER's per-gene noise
    precision, so it enters directly as the weighting Lambda."""
    w, lam = fit["w"], fit["eps"]
    a = np.eye(w.shape[1]) + (w.T * lam) @ w
    factors = np.linalg.solve(a, (w.T * lam) @ x.T).T
    return factors @ w.T


def score_frozen(fit, y_test):
    sf = size_factors(y_test, fit["loggeomeans"])
    x = (np.log2(y_test / sf[:, None] + 1) - fit["center"]) / fit["scale"]
    z = (x - _project(x, fit)) / fit["resid_sd"]
    return dict(z=z, pval=2 * norm.sf(np.abs(z)), sf=sf)


def _demo():
    """Same contract the other arms are held to: the frozen path applied to the training
    samples must reproduce the in-sample residuals it was calibrated on."""
    rng = np.random.default_rng(0)
    n_tr, n_te, n_g, true_k = 120, 30, 300, 3
    loadings = rng.normal(size=(true_k, n_g)) * 1.2
    base = rng.uniform(3, 8, size=n_g)
    make = lambda n: np.round(np.maximum(
        2 ** (base + rng.normal(size=(n, true_k)) @ loadings + rng.normal(scale=0.3, size=(n, n_g))) - 1, 0))
    y_tr, y_te = make(n_tr), make(n_te)

    fit = fit_train(y_tr, "demo", k=10, max_iter=50)
    x_tr, _ = _standardize_train(y_tr)
    z_in = (x_tr - _project(x_tr, fit)) / fit["resid_sd"]
    assert np.abs(score_frozen(fit, y_tr)["z"] - z_in).max() < 1e-10

    out = score_frozen(fit, y_te)
    y_spiked = y_te.copy()
    y_spiked[0, :20] = np.round(y_spiked[0, :20] * 32 + 50)
    spiked = score_frozen(fit, y_spiked)["z"][0, :20]
    assert np.abs(spiked).mean() > np.abs(out["z"][0, :20]).mean(), "spike not detected"
    print(f"peer_arm OK | W={fit['w'].shape} | held-out z sd={out['z'].std():.3f} "
          f"| spiked |z| mean={np.abs(spiked).mean():.2f}")


if __name__ == "__main__":
    _demo()
