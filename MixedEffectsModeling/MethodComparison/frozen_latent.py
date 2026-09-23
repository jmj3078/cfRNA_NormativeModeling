"""NumPy port of frozen_outrider.R::score_frozen, for the autoencoder and PCA arms.

The decoder is the only compiled piece in the R path; OUTRIDER:::predictMatC was verified
to be exactly  normF = exp((x @ E) @ D.T + b) * sf  (2e-13 relative on random inputs), so
the whole scoring path ports to numpy. That matters because normative-mode scoring has to
re-run for every injection level: in R that is a round trip per level, here it is a matmul.

validate() checks the port end-to-end against the mu matrices the R pipeline already wrote
(held_out_comparison/cv_fold{i}_mu.csv), i.e. on real data, not a toy.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import nbinom

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from common import CACHE, size_factors

FITS = CACHE / "latent_fits"


def load_fit(implementation, fold_tag):
    """fold_tag is the split-aware fold key written by fit_latent.R: 'fold0' for CV,
    'lobo_<safe batch dir>' for LOBO."""
    tag = f"{implementation}_{fold_tag}"
    g = pd.read_csv(FITS / f"{tag}_genes.csv")
    fit = {c: g[c].to_numpy(dtype=np.float64) for c in
           ("b", "gene_mu", "theta", "loggeomeans", "x_center", "l2fc_mean", "l2fc_sd")}
    fit["genes"] = g["gene"].tolist()
    fit["E"] = pd.read_csv(FITS / f"{tag}_E.csv").to_numpy(dtype=np.float64)
    fit["D"] = pd.read_csv(FITS / f"{tag}_D.csv").to_numpy(dtype=np.float64)
    return fit


def score_frozen(fit, y):
    """y: samples x genes, aligned to fit['genes']. Returns mu_nb, p, z, sf.

    Normative mode: the sample's own counts drive sf and the encoder input, so an injected
    outlier propagates into its own expected value. That is the behaviour under test."""
    sf = size_factors(y, fit["loggeomeans"])
    x = np.log((1.0 + y) / sf[:, None]) - fit["x_center"]
    normF = np.exp((x @ fit["E"]) @ fit["D"].T + fit["b"]) * sf[:, None]
    mu_nb = fit["gene_mu"] * normF

    theta = np.clip(fit["theta"], 1e-8, None)
    p_nb = theta / (theta + mu_nb)
    pless = nbinom.cdf(y, theta, p_nb)
    dval = nbinom.pmf(y, theta, p_nb)
    pval = 2 * np.minimum(np.minimum(0.5, pless), 1 - pless + dval)

    l2fc = np.log2(y + 1) - np.log2(normF + 1)
    z = (l2fc - fit["l2fc_mean"]) / fit["l2fc_sd"]
    return dict(mu_nb=mu_nb, pval=np.clip(pval, 0, 1), z=z, sf=sf, normF=normF)


def validate(folds=(0,)):
    """Reproduce the R pipeline's own mu on unperturbed held-out counts."""
    counts = pd.read_csv(config.OUTRIDER_COMPARISON_DIR / "hc_counts.csv", index_col=0)
    worst = 0.0
    for fi in folds:
        ref = pd.read_csv(config.OUTRIDER_HELD_OUT_DIR / f"cv_fold{fi}_mu.csv", index_col=0)
        fit = load_fit("autoencoder", f"fold{fi}")
        genes = [g for g in fit["genes"] if g in ref.columns]
        gi = [fit["genes"].index(g) for g in genes]
        y = counts.loc[fit["genes"], ref.index].to_numpy(dtype=np.float64).T

        mu = score_frozen(fit, y)["mu_nb"][:, gi]
        rel = np.abs(mu - ref[genes].to_numpy(dtype=np.float64)) / np.abs(ref[genes].to_numpy(dtype=np.float64))
        worst = max(worst, float(np.nanmax(rel)))
        print(f"  fold{fi}: {len(genes)} genes, max relative mu diff = {np.nanmax(rel):.3e}")
    print(f"{'PASS' if worst < 1e-6 else 'FAIL'}: worst relative diff {worst:.3e}")
    return worst


if __name__ == "__main__":
    validate(folds=range(5))
