"""What the group-level pancreatic-neoplasm signal looks like inside one controlled contrast.

Two scopes side by side: CEDAR alone (Moore Batch_1) and CEDAR+BCC (the older Batch_1+2+3 run),
each PancreaticNeoplasm vs its own controls under four DESeq2 designs. No normative model here.
"""
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config
from analysis.plot_utils import _save

DESIGNS = ["no_covariate", "ruvg_k1", "ruvg_k2", "ruvg_k3"]
DESIGN_SHORT = {"no_covariate": "none", "ruvg_k1": "k=1", "ruvg_k2": "k=2", "ruvg_k3": "k=3"}
COHORTS = {"cedar": config.PANCREATIC_DEG_DIR, "pooled": config.PANCREATIC_DEG_POOLED_DIR}
COHORT_LABEL = {"cedar": "CEDAR", "pooled": "CEDAR+BCC"}
CASE = "PancreaticNeoplasm"
CASE_PHENOTYPE = "Pancreatic Cancer"
HC_PHENOTYPE = "Healthy Control"
SUBSET = f"moore_b1__{CASE}"
DESEQ2_PADJ = 0.05
C_PAT, C_HC = "#B2182B", "#4D4D4D"
SEED = 0


def residualize(Y, W):
    return Y - W @ np.linalg.lstsq(W, Y, rcond=None)[0]


def load(cohort="cedar"):
    deg_dir = COHORTS[cohort]
    cache = pickle.load(open(deg_dir / "pool_cache.pkl", "rb"))
    W = pd.read_csv(deg_dir / "ruvg_W" / f"{SUBSET}.csv", index_col=0)
    obs = cache["obs"]
    keep = np.isin(obs["sample"].values, W.index.values)
    pheno = obs.loc[keep, "phenotype"].values
    tmm = cache["tmm"][keep]
    W_m = W.loc[obs.loc[keep, "sample"].values].values
    genes = np.array(cache["genes"])
    counts = cache["counts"].values[keep]

    res = {d: pd.read_csv(deg_dir / SUBSET / f"{d}.csv.gz", index_col=0) for d in DESIGNS}
    return dict(cohort=cohort, genes=genes, pos={g: i for i, g in enumerate(genes)},
                sym=pd.Series(cache["gene_names"], index=genes), counts=counts,
                cpm=counts / counts.sum(axis=1, keepdims=True) * 1e6,
                X={d: (tmm if d == "no_covariate" else residualize(tmm, W_m[:, :int(d[-1])]))
                   for d in DESIGNS},
                is_case=pheno == CASE_PHENOTYPE, is_hc=pheno == HC_PHENOTYPE, res=res)


def de_genes(D, design):
    r = D["res"][design]
    return sorted(r.index[(r["padj"].fillna(1.0) < DESEQ2_PADJ) & r.index.isin(D["genes"])])


CENTER_LABEL = {"z": "HC-standardized", "median": "TMM-log2 - HC median", "none": "TMM-log2"}


def standardized(D, design, gene_list, center="z"):
    idx = np.array([D["pos"][g] for g in gene_list])
    X = D["X"][design]
    Xh = X[D["is_hc"]][:, idx]
    if center == "none":
        m, s = 0.0, 1.0
    elif center == "median":
        m, s = np.median(Xh, axis=0), 1.0
    else:
        m, s = Xh.mean(axis=0), Xh.std(axis=0, ddof=1)
    return (X[D["is_case"]][:, idx] - m) / s, (Xh - m) / s


def deg_summary(Ds):
    rows = []
    for D in Ds:
        for d in DESIGNS:
            G = de_genes(D, d)
            lfc = D["res"][d]["log2FoldChange"].reindex(G)
            rows.append(dict(scope=COHORT_LABEL[D["cohort"]], design=DESIGN_SHORT[d],
                             n_case=int(D["is_case"].sum()), n_hc=int(D["is_hc"].sum()),
                             deg=len(G), up=int((lfc > 0).sum()), down=int((lfc < 0).sum())))
    return pd.DataFrame(rows)


def plot_top_degs(D, n_show=20, thr=1.96, save_path=None, design="ruvg_k2", rank="lfc", center="z"):
    G = de_genes(D, design)
    if not G:
        return None
    r = D["res"][design].reindex(G)
    order = r["padj"].sort_values()
    if rank == "lfc":
        order = order.reindex(r["log2FoldChange"].sort_values(ascending=False).index)
    if rank == "outliers":
        Za, _ = standardized(D, design, list(order.index), center)
        n_out = (Za > thr).sum(axis=0)
        order = order.iloc[np.argsort(-n_out, kind="stable")]
    top = list(order.index[:n_show])

    sym = [D["sym"].get(g, g.split(".")[0]) for g in top]
    sym = [x if len(x) <= 10 else f"{x[:4]}..{x[-4:]}" for x in sym]
    Zp, Zh = standardized(D, design, top, center)

    ext = Zp > thr
    ext_h = Zh > thr
    lead = np.abs(Zp).argmax(axis=0)
    carriers = sorted(set(np.where(ext)[0].tolist()))
    shuffled = np.random.default_rng(SEED).permutation(len(carriers))
    pat_col = {p: plt.get_cmap("Spectral")(shuffled[i] / max(len(carriers) - 1, 1))
               for i, p in enumerate(carriers)}
    ylo = float(np.floor(min(Zp.min(), Zh.min())))

    fig, (top_ax, ax) = plt.subplots(2, 1, figsize=(max(6.0, 0.56 * len(top) + 2.0), 6.2),
                                     sharex=True, constrained_layout=True,
                                     gridspec_kw=dict(height_ratios=[1.35, 2.0], hspace=0.06))

    zmax = max(float(Zp.max()), thr * 1.2)
    lab_lo, lab_hi = np.log10(thr * 1.15), np.log10(zmax * 1.9)
    for k in range(len(top)):
        oh = np.where(ext_h[:, k])[0]
        top_ax.plot(k + 0.19 + np.linspace(-0.05, 0.05, len(oh)), Zh[oh, k], "o", ms=3.2,
                    color=C_HC, mec="none", alpha=0.75, zorder=2)
        out = np.where(ext[:, k])[0]
        out = out[np.argsort(-Zp[out, k])]
        jit = np.linspace(-0.08, 0.08, len(out)) - 0.19 if len(out) > 1 else np.full(len(out), -0.19)
        slots = (np.linspace(lab_hi, lab_lo, len(out)) if len(out) > 1
                 else np.array([np.log10(Zp[out[0], k]) if len(out) else lab_lo]))
        for j, pi in enumerate(out):
            x, y = k + jit[j], Zp[pi, k]
            top_ax.plot([x, k - 0.36], [y, 10 ** slots[j]], "-", lw=0.5, color="0.7", zorder=1)
            top_ax.text(k - 0.39, 10 ** slots[j], f"P{pi}", va="center", ha="right", fontsize=7.5,
                        color="0.15", zorder=4)
            top_ax.plot(x, y, "o", ms=7 if pi == lead[k] else 5, color=pat_col[pi], mec="k",
                        mew=0.8 if pi == lead[k] else 0.35, zorder=3)

    top_ax.set_yscale("log")
    top_ax.set_ylim(thr, zmax * 2.6)
    top_ax.set_ylabel("outliers")
    top_ax.spines["bottom"].set_visible(False)
    top_ax.tick_params(axis="x", length=0)
    top_ax.set_title(f"{COHORT_LABEL[D['cohort']]}, RUVg {DESIGN_SHORT[design]}: "
                     f"{int(ext.sum())} patient and {int(ext_h.sum())} control values above "
                     f"{thr:g} {'SD' if center == 'z' else 'log2'}; {len(carriers)} of {Zp.shape[0]} "
                     f"patients carry one", fontsize=10, loc="left")

    for k in range(len(top)):
        for arr, off, col in [(Zp[:, k], -0.19, C_PAT), (Zh[:, k], 0.19, C_HC)]:
            inside = arr <= thr
            if np.ptp(arr[inside]) > 0:
                v = ax.violinplot(arr[inside], positions=[k + off], widths=0.34, showextrema=False)
                v["bodies"][0].set_facecolor(col)
                v["bodies"][0].set_alpha(0.32)
            jit = np.random.default_rng(SEED + k).normal(0, 0.03, len(arr))
            ax.plot(k + off + jit[inside], arr[inside], ".", ms=2.4, color=col, alpha=0.6)
    ax.axhline(0, color="k", lw=0.8)
    if thr > 2.1 and center == "z":
        ax.axhline(1.96, color="k", lw=0.7, ls=":")
        ax.text(len(top) - 0.4, 1.96, " 1.96 SD", va="center", ha="left", fontsize=10, color="0.4")
    ax.set_ylim(ylo, thr)
    ax.set_xlim(-0.9, len(top) - 0.4)
    ax.spines["top"].set_visible(False)
    ax.set_xticks(range(len(top)), sym, rotation=90)
    by = {"outliers": "most patient outliers", "lfc": "highest log$_2$FC"}.get(rank, "$p_{adj}$")
    ax.set_xlabel(f"{len(top)} of {len(G)} DE genes ({by}), patients (red) vs controls (grey)")
    ax.set_ylabel(CENTER_LABEL[center])

    for a, y in [(top_ax, 0), (ax, 1)]:
        a.plot([0, 1], [y, y], transform=a.transAxes, ls="none", marker=[(-1, -0.5), (1, 0.5)],
               ms=8, mec="k", mew=1, clip_on=False)
    _save(fig, save_path)
    return fig
