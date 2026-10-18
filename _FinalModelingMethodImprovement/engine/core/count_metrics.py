"""Covariates that are functions of a sample's own counts, recomputed exactly as the merged h5ad
defines them (EDA/1_overview.ipynb with EDA/analysis/bias_metrics.py; reproduced on the h5ad:
NG80 / ratio / GC / length exact, Platelet Score 3e-6, log(Total Reads) 2e-3 (float32 storage)).
Used to move these covariates with an in-silico injection. Read-level metrics (Spliced Reads,
gDNA, rRNA Fraction, 3' bias) are not functions of the count matrix and stay fixed.

CPM = protein-coding counts / their sum * 1e6, log1p (the CPM_log1p layer). Platelet Score =
sc.tl.score_genes on that layer over the whole h5ad; its control genes (random_state 0) are frozen
in splits/platelet_score_control_genes.txt, so a new sample is scored against the same reference.
"""
import h5py
import numpy as np
from scipy.stats import median_abs_deviation
from statsmodels.nonparametric.smoothers_lowess import lowess

import _FinalModelingMethodImprovement.config as config

N_BINS, LOESS_FRAC, OUTLIER_PCT, MIN_EXPRESSED = 20, 0.7, 99, 50  # root config PARAMS
COUNT_METRICS = ("log(Total Reads)", "Platelet Score", "GC Bias", "Gene Length Bias", "NG80", "(NP80/NG80)")
_META = {}


def _meta():
    if not _META:
        from _FinalModelingMethodImprovement.engine.validation.lobo_engine import _decode_cat
        with h5py.File(config.H5AD_PATH, "r") as f:
            names = f["var/_index"][()].astype(str)
            pc = _decode_cat(f["var/GeneType"]) == "protein_coding"
            plt = f["var/is_platelet"][()].astype(bool)
            gc, ln = f["var/GC_Percent"][()], f["var/log10_Length"][()]
        pos = {g: i for i, g in enumerate(names[pc])}
        ctrl = (config.SPLITS_DIR / "platelet_score_control_genes.txt").read_text().split()
        _META.update(pc=pc, gc=gc[pc].astype(float), ln=ln[pc].astype(float), plt=np.where(plt[pc])[0],
                     ctrl=np.array([pos[g] for g in ctrl if g in pos]), n_ctrl=len(ctrl))
    return _META


def n80(Y):
    out = np.zeros(len(Y))
    for i, row in enumerate(Y):
        cs = np.cumsum(np.sort(row[row > 0])[::-1])
        if len(cs) and cs[-1] > 0:
            out[i] = np.argmax(cs >= cs[-1] * 0.8) + 1
    return out


def curve_score(L, feat):
    """bias_metrics._compute_score_sparse on dense rows."""
    import pandas as pd
    out = []
    for row in L:
        pos = row > 0
        expr, f = row[pos], feat[pos]
        if len(expr) < MIN_EXPRESSED:
            out.append(0.0)
            continue
        keep = expr <= np.percentile(expr, OUTLIER_PCT)
        df = pd.DataFrame({"expr": expr[keep], "feat": f[keep]})
        try:
            df["bin"] = pd.qcut(df["feat"], q=N_BINS, duplicates="drop")
        except ValueError:
            out.append(0.0)
            continue
        b = df.groupby("bin", observed=True).agg({"expr": "median", "feat": "mean"}).dropna()
        if len(b) < 2:
            out.append(0.0)
            continue
        curve = median_abs_deviation(lowess(b["expr"], b["feat"], frac=LOESS_FRAC, it=0)[:, 1], scale="normal")
        total = median_abs_deviation(df["expr"], scale="normal")
        out.append(curve / total if total > 0 else 0.0)
    return np.array(out)


def count_metrics(Y, names):
    """{name: per-row value} for the requested COUNT_METRICS; Y = raw counts over all h5ad genes."""
    m = _meta()
    Yp = Y[:, m["pc"]]
    L = np.log1p(Yp / Yp.sum(1, keepdims=True) * 1e6)
    out = {}
    if "log(Total Reads)" in names:
        out["log(Total Reads)"] = np.log(Y.sum(1))
    if "Platelet Score" in names:
        out["Platelet Score"] = L[:, m["plt"]].mean(1) - L[:, m["ctrl"]].sum(1) / m["n_ctrl"]
    if "GC Bias" in names:
        out["GC Bias"] = curve_score(L, m["gc"])
    if "Gene Length Bias" in names:
        out["Gene Length Bias"] = curve_score(L, m["ln"])
    if "NG80" in names or "(NP80/NG80)" in names:
        ng = n80(Y)
        out["NG80"] = ng
        with np.errstate(divide="ignore", invalid="ignore"):
            out["(NP80/NG80)"] = np.where(ng == 0, 0, n80(Yp) / ng)
    return out
