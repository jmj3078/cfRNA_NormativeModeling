"""Adaptation LOBO: every method may use half of the held-out batch's HC to recalibrate.

Generic, refit-free recalibration applied identically to every method: p -> signed z
(sign of the method's own z), per-gene mean m of the adaptation half, z' = z - m,
p' = 2 sf(|z'|). Mean only: on 2026-10-16 mean+SD / median / median+MAD were also tried
on the in-silico LOBO cells; scaling raised the engine's null FP (42 -> 72 per sample,
MAD -> 221) and location-only was the only variant that did not hurt any method's
precision much, so it is the fixed choice.

Cross-fitted 2-fold: each HC half is shifted by the other half's mean, so no sample
recalibrates itself. The adaptation half is read from the unperturbed (log2fc=0) cell, so
injected signal never enters m. BLUP pipelines (rqr="blup") adapt inside the model and are
not passed through this.

Run: python eval/adaptation.py insilico
"""
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _FinalModelingMethodImprovement.config as config
from _FinalModelingMethodImprovement.eval.insilico_sweep import (
    LOG2FCS, PROD_MC, Q_LEVELS, count_rejections, fold_spec)
from _FinalModelingMethodImprovement.eval.results import append

ADAPT_SEED = 20261016


def halves(n, key):
    """Fixed 2-fold split of a batch's HC, keyed by batch name."""
    perm = np.random.default_rng([ADAPT_SEED, zlib.crc32(key.encode())]).permutation(n)
    return np.sort(perm[: n // 2]), np.sort(perm[n // 2:])


def signed_z(z, p):
    return np.sign(z) * norm.isf(np.clip(p.astype(np.float64) / 2, 1e-300, 0.5))


def loc(z_ref):
    m = np.nanmean(z_ref, axis=0)
    return np.where(np.isfinite(m), m, 0.0)


def cross_fit(z_eval, z_ref, idx):
    """z_eval/z_ref: samples x genes over the same samples. Each half is shifted by the
    other half's per-gene mean of z_ref."""
    out = np.full_like(z_eval, np.nan, dtype=np.float64)
    a, b = idx
    for tgt, src in ((a, b), (b, a)):
        out[tgt] = z_eval[tgt] - loc(z_ref[src])
    return out


def method_cells(method, tag):
    """{log2fc: signed z} for one LOBO batch: outlier arms from the production matrices,
    engine pipelines from their own zero-shot in-silico Z."""
    out = {}
    for f in LOG2FCS:
        if method in config.OUTLIER_ARMS:
            c = np.load(PROD_MC / "matrices" / f"lobo_{tag}_log2fc{f:g}.npz")
            out[f] = signed_z(c[f"z__{method}"], c[f"p__{method}"])
        else:
            z = np.load(config.FMMI / "runs" / method / "insilico" / f"lobo_{tag}_log2fc{f:g}.npz")["z"]
            out[f] = signed_z(z, 2 * norm.sf(np.abs(z.astype(np.float64))))
    return out


def insilico():
    methods = list(config.OUTLIER_ARMS) + [
        p for p, spec in config.PIPELINES.items()
        if spec["rqr"] == "marginal" and (config.FMMI / "runs" / p / "insilico").exists()]
    rows = []
    for key, tag, _, te in fold_spec("lobo"):
        cells = {f: np.load(config.INSILICO_DIR / "sweep" / f"lobo_{tag}_log2fc{f:g}.npz") for f in LOG2FCS}
        idx = halves(len(te), key)
        for method in methods:
            z = method_cells(method, tag)
            for f, c in cells.items():
                p = 2 * norm.sf(np.abs(cross_fit(z[f], z[0.0], idx)))
                rows += [dict(split="lobo", method=method, scenario="adaptation", fold=key, log2fc=f,
                              n_samples=len(te), n_universe=p.shape[1], n_perturbed=int(c["labels"].sum()), **r)
                         for r in count_rejections(p, c["labels"], Q_LEVELS)]
        print(f"adapted {tag}: {methods}", flush=True)
    append(config.INSILICO_CSV, pd.DataFrame(rows))


if __name__ == "__main__":
    {"insilico": insilico}[sys.argv[1]]()
