"""Covariates that are functions of a sample's own counts, recomputed exactly as the merged h5ad
defines them (EDA/1_overview.ipynb with EDA/analysis/bias_metrics.py). Used to move these
covariates with an in-silico injection, the way they would move in real data. Read-level metrics
(Spliced Reads, gDNA, rRNA Fraction, 3' bias) are not functions of the count matrix and stay fixed.

CPM = protein-coding counts / their sum * 1e6, log1p (the CPM_log1p layer). Platelet Score =
sc.tl.score_genes on that layer over the whole h5ad; its control genes (random_state 0) are frozen
in splits/platelet_score_control_genes.txt and checked against the stored scores on
first use. After an h5ad rebuild: freeze_platelet_controls().
"""
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import median_abs_deviation
from statsmodels.nonparametric.smoothers_lowess import lowess

import _FinalModelingMethodImprovement.config as config

N_BINS, LOESS_FRAC, OUTLIER_PCT, MIN_EXPRESSED = 20, 0.7, 99, 50  # root config PARAMS
COUNT_METRICS = ("log(Total Reads)", "Platelet Score", "GC Bias", "Gene Length Bias", "NG80", "(NP80/NG80)")
CONTROL_LIST = config.SPLITS_DIR / "platelet_score_control_genes.txt"
_META = {}


def _decode(grp):
    from _FinalModelingMethodImprovement.engine.validation.lobo_engine import _decode_cat
    return _decode_cat(grp)


def freeze_platelet_controls():
    """Control genes of sc.tl.score_genes(CPM_log1p, is_platelet genes, random_state=0) on the h5ad."""
    import scanpy as sc
    from scanpy.tools._score_genes import _check_score_genes_args, _score_genes_bins
    ad = sc.read_h5ad(config.H5AD_PATH)
    tmp = sc.AnnData(X=ad.layers["CPM_log1p"], var=ad.var)
    np.random.seed(0)
    genes, pool, get_subset = _check_score_genes_args(tmp, ad.var_names[ad.var["is_platelet"]].tolist(), None,
                                                      use_raw=False, layer=None)
    ctrl = set()
    for r in _score_genes_bins(genes, pool, ctrl_as_ref=True, ctrl_size=50, n_bins=25, get_subset=get_subset):
        ctrl |= set(r)
    CONTROL_LIST.write_text("\n".join(sorted(ctrl)))
    _META.clear()


def _meta():
    if not _META:
        with h5py.File(config.H5AD_PATH, "r") as f:
            names = f["var/_index"][()].astype(str)
            pc = _decode(f["var/GeneType"]) == "protein_coding"
            plt = f["var/is_platelet"][()].astype(bool)
            gc, ln = f["var/GC_Percent"][()], f["var/log10_Length"][()]
            ptr = f["X/indptr"][()]
            rows = np.linspace(0, len(ptr) - 2, 20).astype(int)
            X = np.zeros((len(rows), len(names)))
            for i, r in enumerate(rows):
                X[i, f["X/indices"][ptr[r]:ptr[r + 1]]] = f["X/data"][ptr[r]:ptr[r + 1]]
            stored = f["obs/Platelet Score"][()][rows]
        pos = {g: i for i, g in enumerate(names[pc])}
        ctrl = CONTROL_LIST.read_text().split()
        _META.update(pc=pc, gc=gc[pc].astype(float), ln=ln[pc].astype(float), plt=np.where(plt[pc])[0],
                     ctrl=np.array([pos[g] for g in ctrl if g in pos]), n_ctrl=len(ctrl))
        err = np.abs(count_metrics(X, ["Platelet Score"])["Platelet Score"] - stored).max()
        assert err < 1e-4, f"Platelet Score controls do not reproduce the h5ad ({err:.2g}); run freeze_platelet_controls()"
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
    m = _META if _META else _meta()
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


def shift_covariates(X, columns, Y0, Y1):
    """X (rows x columns) with every count-derived column moved by f(Y1) - f(Y0): the covariates of
    the same samples after their counts went from Y0 to Y1. Unchanged rows stay bit-identical."""
    X = np.array(X, dtype=np.float64)
    cols = [c for c in columns if c in COUNT_METRICS]
    changed = np.where((Y0 != Y1).any(1))[0]
    if cols and len(changed):
        new, old = count_metrics(Y1[changed], cols), count_metrics(Y0[changed], cols)
        for c in cols:
            X[changed, list(columns).index(c)] += new[c] - old[c]
    return X
