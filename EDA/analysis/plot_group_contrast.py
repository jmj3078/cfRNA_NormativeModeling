"""Group-contrast limits on the DESeq2 scale alone: no normative model anywhere in here.

Expression is TMM-log2 residualized on each run's own RUVg W, standardized on the
batch-matched healthy controls, so every statistic below is one DESeq2 could have computed
on its own input.
"""
import pickle

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config
from analysis.plot_utils import _save

DESIGNS = ["no_covariate", "ruvg_k1", "ruvg_k2", "ruvg_k3"]
DESIGN_LABEL = {"no_covariate": "DESeq2 Default", "ruvg_k1": "DESeq2 + RUVg (k=1)",
                "ruvg_k2": "DESeq2 + RUVg (k=2)", "ruvg_k3": "DESeq2 + RUVg (k=3)"}
DESIGN_C = {"no_covariate": "#A2A2A2", "ruvg_k1": "#489ACA", "ruvg_k2": "#009E73",
            "ruvg_k3": "#CC79A7"}
PICK = "ruvg_k2"
CASE_PHENOTYPE = "Pancreatic Cancer"
HC_PHENOTYPE = "Healthy Control"
SUBSET = "moore_b1__PancreaticNeoplasm"
DESEQ2_PADJ = 0.05
C_PAT, C_HC = "#B2182B", "#4D4D4D"
SEED = 0


def residualize(Y, W):
    return Y - W @ np.linalg.lstsq(W, Y, rcond=None)[0]


def load():
    cache = pickle.load(open(config.CTRL_COMP_DIR / "moore_b1_cache.pkl", "rb"))
    W = pd.read_csv(config.PANCREATIC_DEG_DIR / "ruvg_W" / f"{SUBSET}.csv", index_col=0)
    obs = cache["obs"]
    keep = np.isin(obs["sample"].values, W.index.values)
    pheno = obs.loc[keep, "phenotype"].values
    tmm = cache["layers"]["TMM_log2"][keep]
    W_m = W.loc[obs.loc[keep, "sample"].values].values
    genes = np.array(cache["genes"])

    deg, padj = {}, {}
    for f in sorted(config.PANCREATIC_DEG_DIR.glob("moore_b1__*/*.csv.gz")):
        r = pd.read_csv(f, index_col=0)
        key = (f.parent.name.split("__")[1], f.name.split(".")[0])
        deg[key] = set(r.index[r["padj"].fillna(1.0) < DESEQ2_PADJ]) & set(genes)
        padj[key] = r["padj"]

    return dict(genes=genes, pos={g: i for i, g in enumerate(genes)},
                X={d: (tmm if d == "no_covariate" else residualize(tmm, W_m[:, :int(d[-1])]))
                   for d in DESIGNS},
                is_case=pheno == CASE_PHENOTYPE, is_hc=pheno == HC_PHENOTYPE,
                deg=deg, padj=padj)


def de_genes(D, design):
    return sorted(D["deg"][("PancreaticNeoplasm", design)] | D["deg"][("PDAC", design)])


def standardized(D, design, gene_list):
    idx = np.array([D["pos"][g] for g in gene_list])
    X = D["X"][design]
    Xh = X[D["is_hc"]][:, idx]
    m, s = Xh.mean(axis=0), Xh.std(axis=0, ddof=1)
    return (X[D["is_case"]][:, idx] - m) / s, (Xh - m) / s


def representativeness(D, save_path=None):
    """|mean|/SD and sign concordance of the group mean over the run's own DE genes."""
    rows = []
    for design in DESIGNS:
        G = de_genes(D, design)
        Z, _ = standardized(D, design, G)
        mz = Z.mean(axis=0)
        ams = np.abs(mz) / Z.std(axis=0, ddof=1)
        sc = (np.sign(Z) == np.sign(mz)[None, :]).mean(axis=0)
        rows.append(dict(design=design, n_genes=len(G), n_case=int(D["is_case"].sum()),
                         abs_mean_over_sd_med=np.median(ams), frac_ams_lt_05=(ams < 0.5).mean(),
                         sign_concordance_med=np.median(sc), frac_sc_lt_05=(sc < 0.5).mean(),
                         frac_sc_lt_06=(sc < 0.6).mean()))
    out = pd.DataFrame(rows)
    if save_path:
        out.to_csv(save_path, index=False)
    return out


def plot_direction(D, save_path=None):
    G = de_genes(D, PICK)
    Z, _ = standardized(D, PICK, G)
    mz = Z.mean(axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4),
                             gridspec_kw=dict(width_ratios=[1, 1.0], wspace=0.35))
    ax = axes[0]
    im = ax.imshow(Z[:, np.argsort(mz)], cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto",
                   interpolation="nearest")
    ax.set_xlabel(f"{len(G)} group DE genes (sorted by group mean)")
    ax.set_ylabel(f"{Z.shape[0]} patients")
    ax.set_xticks([])
    ax.set_yticks([])
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("HC-standardized TMM-log2")
    for s in ax.spines.values():
        s.set_visible(False)

    ax = axes[1]
    groups, labels, cols = [], [], []
    for design in DESIGNS:
        Gd = de_genes(D, design)
        if not Gd:
            continue
        Zd, _ = standardized(D, design, Gd)
        md = Zd.mean(axis=0)
        groups.append((np.sign(Zd) == np.sign(md)[None, :]).mean(axis=0))
        labels.append(f"{DESIGN_LABEL[design]}\n{len(Gd)} genes")
        cols.append(DESIGN_C[design])

    bp = ax.boxplot(groups, positions=range(len(groups)), widths=0.55, patch_artist=True,
                    flierprops=dict(marker="o", ms=1.5, mfc="0.45", mec="none"),
                    medianprops=dict(color="k", lw=2.5))
    for box, c in zip(bp["boxes"], cols):
        box.set(facecolor=c, alpha=0.7, edgecolor=c, lw=1.2)
    for i, (v, c) in enumerate(zip(groups, cols)):
        ax.text(i, 1.03, f"{np.median(v):.2f}", ha="center", color=c)
    ax.axhline(0.5, color="k", lw=0.8, ls=":")
    ax.axhline(1.0, color="#333333", lw=0.8, ls="--")
    ax.set_xticks(range(len(groups)), labels)
    ax.set_ylim(min(0.0, min(v.min() for v in groups) - 0.05), 1.1)
    ax.set_ylabel("patients on the group side")
    for a, l in zip(axes, "ab"):
        a.text(-0.08, 1.10, l, transform=a.transAxes, fontweight="bold", va="top")
    _save(fig, save_path)
    return fig, pd.DataFrame(dict(design=DESIGNS, median=[np.median(v) for v in groups],
                                  frac_below_06=[(v < 0.6).mean() for v in groups]))


def plot_top_degs(D, n_show=20, save_path=None):
    G = de_genes(D, PICK)
    pa = D["padj"][("PancreaticNeoplasm", PICK)].combine_first(D["padj"][("PDAC", PICK)])
    top = list(pa.reindex(G).sort_values().index[:n_show])

    with h5py.File(config.H5AD_PATH, "r") as h:
        dec = lambda a: np.array([x.decode("utf-8", "replace") for x in a])
        gn = h["var/GeneName"]
        sym_of = pd.Series(pd.Categorical.from_codes(gn["codes"][()],
                                                     dec(gn["categories"][()])).astype(str),
                           index=dec(h["var/_index"][()]))
    sym = [sym_of.get(g, g.split(".")[0]) for g in top]
    Zp, Zh = standardized(D, PICK, top)

    ylo, yhi = np.percentile(np.concatenate([Zp.ravel(), Zh.ravel()]), [0.5, 99.5])
    ylo, yhi = float(np.floor(ylo - 1)), float(np.ceil(yhi + 1))
    n_off = int((Zp > yhi).sum() + (Zp < ylo).sum() + (Zh > yhi).sum() + (Zh < ylo).sum())

    fig, ax = plt.subplots(figsize=(14.0, 5.0))
    for k in range(len(top)):
        for arr, off, col in [(Zp[:, k], -0.2, C_PAT), (Zh[:, k], 0.2, C_HC)]:
            v = ax.violinplot(np.clip(arr, ylo, yhi), positions=[k + off], widths=0.36,
                              showextrema=False)
            v["bodies"][0].set_facecolor(col)
            v["bodies"][0].set_alpha(0.35)
            jit = k + off + np.random.default_rng(SEED + k).normal(0, 0.035, len(arr))
            inside = (arr >= ylo) & (arr <= yhi)
            ax.plot(jit[inside], arr[inside], ".", ms=2.2, color=col, alpha=0.7)
            for m, sel in [("^", arr > yhi), ("v", arr < ylo)]:
                ax.plot(jit[sel], np.full(sel.sum(), yhi if m == "^" else ylo), m, ms=3.5,
                        color=col, alpha=0.9, clip_on=False)
    ax.set_ylim(ylo, yhi)
    ax.text(0.0, 1.02, f"{n_off} points outside the axis, drawn as triangles at the edge",
            transform=ax.transAxes, fontsize=12, color="0.35")
    ax.axhline(0, color="k", lw=0.8)
    for y in (-1.96, 1.96):
        ax.axhline(y, color="k", lw=0.7, ls=":")
    ax.text(len(top) - 0.45, 2.1, "1.96 SD of healthy", va="bottom", ha="right", fontsize=12)
    ax.set_xticks(range(len(top)), sym, rotation=45, ha="right")
    ax.set_xlabel(f"top {len(top)} group DE genes (sorted by $p_{{adj}}$)")
    ax.set_ylabel("HC-standardized TMM-log2")
    _save(fig, save_path)
    return fig


def centering_contrast(D, design=PICK, save_path=None):
    """Why the centering is not a free choice: without it both statistics change question."""
    G = de_genes(D, design)
    idx = np.array([D["pos"][g] for g in G])
    X = D["X"][design]
    Xp, Xh = X[D["is_case"]][:, idx], X[D["is_hc"]][:, idx]
    variants = {"none (raw TMM-log2)": Xp,
                "HC mean": Xp - Xh.mean(axis=0),
                "HC mean and SD": (Xp - Xh.mean(axis=0)) / Xh.std(axis=0, ddof=1),
                "all-sample mean": Xp - X[:, idx].mean(axis=0)}
    rows = []
    for lab, Z in variants.items():
        m = Z.mean(axis=0)
        rows.append(dict(centering=lab, abs_mean_over_sd_med=np.median(np.abs(m) / Z.std(axis=0, ddof=1)),
                         sign_concordance_med=np.median((np.sign(Z) == np.sign(m)[None, :]).mean(axis=0))))
    out = pd.DataFrame(rows)
    if save_path:
        out.to_csv(save_path, index=False)
    return out


def outlier_decomposition(D, design=PICK, thr=6.0, save_path=None):
    """Is the spread a few outlier patients, or everyone? Leave-one-patient-out on the mean."""
    G = de_genes(D, design)
    Z, _ = standardized(D, design, G)
    ext = np.abs(Z) > thr
    m = Z.mean(axis=0)
    drop = np.array([np.abs(np.delete(Z, j, 0).mean(axis=0) - m).mean() for j in range(Z.shape[0])])
    out = pd.DataFrame([dict(
        design=design, thr=thr, n_extreme=int(ext.sum()), frac_extreme=float(ext.mean()),
        patients_with_extreme=int((ext.sum(axis=1) > 0).sum()), n_patients=Z.shape[0],
        genes_with_extreme=int((ext.sum(axis=0) > 0).sum()), n_genes=len(G),
        loo_mean_shift_med=float(np.median(drop)), loo_mean_shift_max=float(drop.max()))])
    if save_path:
        out.to_csv(save_path, index=False)
    return out
