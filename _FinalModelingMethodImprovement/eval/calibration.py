"""HC calibration of the engine pipelines, 1_cv_analysis.ipynb (CV) and 2_lobo_calibration.ipynb
(LOBO) metrics, written as long tables for 2_calibration.ipynb. Engine pipelines only.

splits   cv        every HC of the 5-fold CV (nb1)
         cv_paired CV Z of exactly the LOBO HC kept by the OOD filter (nb2 "CV (seen batch)")
         lobo      LOBO Z of the same HC (nb2 "LOBO (unseen batch)")
scenario zero_shot as scored; adaptation (lobo only) = shifted by the other half's per-gene
         mean within each batch (eval/adaptation.py).

Outputs (config.RESULTS_DIR): calibration_gene.csv, calibration_sample.csv,
calibration_exceedance.csv, calibration_lobo_batch.csv, calibration_ppc.csv.
Run: FMMI_PIPELINE=<pipe> python eval/calibration.py
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kurtosis, norm, skew

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _FinalModelingMethodImprovement.config as config
from _FinalModelingMethodImprovement.engine.core.calibration import bh_fdr_reject
from _FinalModelingMethodImprovement.engine.validation.lobo_engine import safe_dir
from _FinalModelingMethodImprovement.engine.validation.ppc_simulate import simulate_marginal_nb
from _FinalModelingMethodImprovement.eval.adaptation import halves, loc
from _FinalModelingMethodImprovement.eval.results import append

Q_GENE = (0.01, 0.005)
Q_SAMPLE = 0.05
ALPHAS = np.logspace(-4, np.log10(0.5), 25)
PPC_N_REPS = 200
OUT = {k: config.RESULTS_DIR / f"calibration_{k}.csv" for k in ("gene", "sample", "exceedance", "lobo_batch", "ppc")}


def bh_count(v, q):
    v = v[np.isfinite(v)]
    return int(bh_fdr_reject(2 * norm.sf(np.abs(v)), q=q).sum()) if len(v) else 0


def gene_table(z, genes, min_n=1):
    """nb1 cell 5 / nb2 cell 5: per-gene moments over finite samples, exceedance, BH outliers."""
    rows = []
    for j, g in enumerate(genes):
        v = z[:, j][np.isfinite(z[:, j])].astype(np.float64)
        if len(v) < min_n:
            continue
        rows.append(dict(gene=g, n=len(v), mean=v.mean(), std=v.std(), skew=skew(v), kurt=kurtosis(v),
                         exceed_196=(np.abs(v) > 1.96).mean(),
                         **{f"nout_q{int(round(q * 1000)):03d}": bh_count(v, q) for q in Q_GENE}))
    return pd.DataFrame(rows)


def load_cv(run):
    folds = json.loads(config.CV_FOLDS_PATH.read_text())
    order = np.concatenate([folds[k]["test"] for k in sorted(folds, key=int)])
    raw = pickle.load(open(run / "cv" / "cv_zscores.pkl", "rb"))
    sh = pickle.load(open(run / "cv" / "cv_zscores_shash.pkl", "rb"))
    return order, raw, sh


def load_lobo(run, genes_cv):
    """nb2 cell 3: OOD-kept HC of each frozen LOBO batch, genes shared by CV and every batch."""
    batches = sorted(json.loads(config.LOBO_FOLDS_PATH.read_text()))
    metas = [json.loads((run / "lobo" / safe_dir(b) / "meta.json").read_text()) for b in batches]
    gl = [pickle.load(open(run / "lobo" / safe_dir(b) / "gene_names.pkl", "rb")) for b in batches]
    genes = sorted(set.intersection(*map(set, gl)) & set(genes_cv))
    parts = []
    for b, m, g in zip(batches, metas, gl):
        bdir = run / "lobo" / safe_dir(b)
        keep = np.array(m["test_is_hc"]) & np.load(bdir / "ood_mask.npy")
        gi = pd.Index(g).get_indexer(genes)
        parts.append(dict(batch=b, names=np.array(m["test_names"])[keep],
                          raw=np.load(bdir / "Z_test.npy")[keep][:, gi].astype(float),
                          shash=np.load(bdir / "Z_test_shash.npy")[keep][:, gi].astype(float)))
    return genes, parts


def adapt(parts, kind):
    """Cross-fitted per-batch mean shift on the OOD-kept HC (split keyed by batch)."""
    out = []
    for p in parts:
        z = p[kind]
        a, b = halves(len(z), p["batch"])
        za = np.empty_like(z)
        za[a], za[b] = z[a] - loc(z[b]), z[b] - loc(z[a])
        out.append(za)
    return np.vstack(out)


def ppc(run):
    """nb1 cell 4: per-gene posterior-predictive mean / SD / zero fraction on the CV held-out HC."""
    rows = []
    for i, (g, d) in enumerate(pickle.load(open(run / "cv" / "cv_ppc.pkl", "rb")).items()):
        y, mu, alpha, tau2 = (np.asarray(d[k], float) for k in ("y", "mu", "alpha", "tau2"))
        if y.size == 0 or not (np.isfinite(mu).all() and np.isfinite(alpha).all() and np.isfinite(tau2).all()):
            continue
        y_rep = simulate_marginal_nb(mu, alpha, tau2, PPC_N_REPS, seed=1000 + i)
        rows.append(dict(gene=g, obs_mean=y.mean(), pred_mean=y_rep.mean(1).mean(), obs_sd=y.std(),
                         pred_sd=y_rep.std(1).mean(), obs_zero=(y == 0).mean(), pred_zero=(y_rep == 0).mean(1).mean()))
    return pd.DataFrame(rows)


def main():
    pipe = config.PIPELINE
    run = config.RUN_DIR
    order, cv_raw, cv_sh = load_cv(run)
    full = [g for g, z in cv_raw.items() if len(z) == len(order)]
    genes, parts = load_lobo(run, full)
    names = np.concatenate([p["names"] for p in parts])
    batch = np.concatenate([np.repeat(p["batch"], len(p["names"])) for p in parts])
    pos = pd.Series(np.arange(len(order)), index=order).loc[names].values

    Z = {("cv", "zero_shot", "raw"): np.column_stack([cv_raw[g] for g in full]).astype(float),
         ("cv", "zero_shot", "shash"): np.column_stack([cv_sh[g] for g in full]).astype(float),
         ("cv_paired", "zero_shot", "raw"): np.column_stack([cv_raw[g][pos] for g in genes]).astype(float),
         ("cv_paired", "zero_shot", "shash"): np.column_stack([cv_sh[g][pos] for g in genes]).astype(float)}
    for kind in ("raw", "shash"):
        Z[("lobo", "zero_shot", kind)] = np.vstack([p[kind] for p in parts])
        Z[("lobo", "adaptation", kind)] = adapt(parts, kind)
    for z in Z.values():
        z[~np.isfinite(z)] = np.nan

    tag = dict(method=pipe)
    gene_rows, sample_rows, exc_rows, batch_rows = [], [], [], []
    for (split, scen, kind), z in Z.items():
        g = full if split == "cv" else genes
        gene_rows.append(gene_table(z, g, min_n=20 if split == "cv" else 1).assign(split=split, scenario=scen, kind=kind, **tag))
        if split == "cv" or kind != "shash":
            continue
        nsig = np.array([bh_count(r, Q_SAMPLE) for r in z])
        sample_rows.append(pd.DataFrame(dict(sample=names, batch=batch, n_sig=nsig, n_genes=len(genes),
                                             split=split, scenario=scen, **tag)))
        for b in dict.fromkeys(batch):
            zb = z[batch == b]
            p = 2 * norm.sf(np.abs(zb[np.isfinite(zb)]))
            exc_rows += [dict(batch=b, alpha=a, exceed=(p < a).mean(), split=split, scenario=scen, **tag) for a in ALPHAS]
            nb = nsig[batch == b]
            batch_rows.append(dict(batch=b, split=split, scenario=scen, n_hc=int((batch == b).sum()),
                                   median_gene_mean=np.nanmedian(np.nanmean(zb, 0)),
                                   median_gene_sd=np.nanmedian(np.nanstd(zb, 0)),
                                   exceed_196=np.nanmean(np.abs(zb) > 1.96), exceed_329=np.nanmean(np.abs(zb) > 3.29),
                                   nsig_median=np.median(nb), nsig_mean=nb.mean(),
                                   frac_samples_nsig0=(nb == 0).mean(), **tag))
    append(OUT["gene"], pd.concat(gene_rows, ignore_index=True))
    append(OUT["sample"], pd.concat(sample_rows, ignore_index=True))
    append(OUT["exceedance"], pd.DataFrame(exc_rows))
    append(OUT["lobo_batch"], pd.DataFrame(batch_rows))
    ppc_cache = run / "cv" / "cv_ppc_moments.csv"
    if not ppc_cache.exists():
        ppc(run).to_csv(ppc_cache, index=False)
    append(OUT["ppc"], pd.read_csv(ppc_cache).assign(split="cv", scenario="zero_shot", **tag))
    print(f"calibration tables written for {pipe}: {list(OUT.values())}")


if __name__ == "__main__":
    main()
