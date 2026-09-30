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
DESIGN_SHORT = {"no_covariate": "none", "ruvg_k1": "k=1", "ruvg_k2": "k=2", "ruvg_k3": "k=3"}
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

    curves, labels, cols = [], [], []
    for design in DESIGNS:
        Gd = de_genes(D, design)
        if not Gd:
            continue
        Zd, _ = standardized(D, design, Gd)
        md = Zd.mean(axis=0)
        curves.append((np.sign(Zd) == np.sign(md)[None, :]).mean(axis=0))
        labels.append(design)
        cols.append(DESIGN_C[design])

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.8), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[1.3, 1.0]))

    ax = axes[0]
    im = ax.imshow(Z[:, np.argsort(mz)], cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto",
                   interpolation="nearest")
    ax.set_xlabel(f"{len(G)} group DE genes (sorted by group mean)")
    ax.set_ylabel(f"{Z.shape[0]} patients")
    ax.set_xticks([])
    ax.set_yticks([])
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("HC-standardized")
    for s in ax.spines.values():
        s.set_visible(False)

    ax = axes[1]
    parts = ax.violinplot(curves, positions=range(len(curves)), widths=0.85, showextrema=False)
    for b, c in zip(parts["bodies"], cols):
        b.set_facecolor(c)
        b.set_alpha(0.35)
    for i, (v, c) in enumerate(zip(curves, cols)):
        jit = np.random.default_rng(SEED + i).normal(0, 0.055, len(v))
        ax.plot(i + jit, v, ".", ms=2.6, color=c, alpha=0.65)
        ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color="k", lw=2.2, zorder=5)
        ax.text(i, 1.06, f"{np.median(v):.2f}", ha="center", fontsize=11, color=c)
        ax.text(i, -0.07, f"{(v < 0.5).mean():.0%} < 0.5", ha="center", fontsize=10,
                color="0.35")
    ax.axhline(0.5, color="k", lw=0.9, ls=":")
    ax.set_xticks(range(len(curves)),
                  [f"{DESIGN_SHORT[l]}\n{len(v)}" for l, v in zip(labels, curves)])
    ax.set_xlabel("RUVg factors, and the genes each design calls")
    ax.set_ylim(-0.12, 1.14)
    ax.set_ylabel("patients on the group side")
    ax.text(-0.45, 0.5, "0.5 = no shared direction", va="bottom", ha="left", fontsize=10,
            color="0.35")

    for a, l in zip(axes, "ab"):
        a.text(-0.05, 1.12, l, transform=a.transAxes, fontweight="bold", va="top")
    _save(fig, save_path)
    return fig, pd.DataFrame(dict(design=labels, median=[np.median(v) for v in curves],
                                  frac_below_05=[(v < 0.5).mean() for v in curves],
                                  frac_below_06=[(v < 0.6).mean() for v in curves]))


def plot_top_degs(D, n_show=20, thr=3.0, save_path=None, design=PICK):
    """Broken axis: the bulk keeps a linear panel, the outliers get a log panel of their own,
    each coloured by the patient carrying it."""
    G = de_genes(D, design)
    pa = D["padj"][("PancreaticNeoplasm", design)].combine_first(D["padj"][("PDAC", design)])
    top = list(pa.reindex(G).sort_values().index[:n_show])

    with h5py.File(config.H5AD_PATH, "r") as h:
        dec = lambda a: np.array([x.decode("utf-8", "replace") for x in a])
        gn = h["var/GeneName"]
        sym_of = pd.Series(pd.Categorical.from_codes(gn["codes"][()],
                                                     dec(gn["categories"][()])).astype(str),
                           index=dec(h["var/_index"][()]))
    sym = [sym_of.get(g, g.split(".")[0]) for g in top]
    sym = [x if len(x) <= 10 else f"{x[:4]}..{x[-4:]}" for x in sym]
    Zp, Zh = standardized(D, design, top)

    ext = Zp > thr
    lead = np.abs(Zp).argmax(axis=0)
    carriers = sorted(set(np.where(ext)[0].tolist()))
    shuffled = np.random.default_rng(SEED).permutation(len(carriers))
    pat_col = {p: plt.get_cmap("Spectral")(shuffled[i] / max(len(carriers) - 1, 1))
               for i, p in enumerate(carriers)}
    ylo = float(np.floor(np.percentile(np.concatenate([Zp.ravel(), Zh.ravel()]), 0.2)))

    fig, (top_ax, ax) = plt.subplots(2, 1, figsize=(14.0, 6.2), sharex=True,
                                     constrained_layout=True,
                                     gridspec_kw=dict(height_ratios=[1.35, 2.0], hspace=0.06))

    lab_lo, lab_hi = np.log10(thr * 1.15), np.log10(float(Zp.max()) * 1.9)
    for k in range(len(top)):
        out = np.where(ext[:, k])[0]
        out = out[np.argsort(-Zp[out, k])]
        jit = np.linspace(-0.11, 0.11, len(out)) if len(out) > 1 else np.zeros(len(out))
        slots = (np.linspace(lab_hi, lab_lo, len(out)) if len(out) > 1
                 else np.array([np.log10(Zp[out[0], k]) if len(out) else lab_lo]))
        for j, pi in enumerate(out):
            x, y = k + jit[j], Zp[pi, k]
            top_ax.plot([x, k + 0.30], [y, 10 ** slots[j]], "-", lw=0.5, color="0.7", zorder=1)
            top_ax.plot(x, y, "o", ms=7 if pi == lead[k] else 5, color=pat_col[pi], mec="k",
                        mew=0.8 if pi == lead[k] else 0.35, zorder=3)
            top_ax.text(k + 0.33, 10 ** slots[j], f"P{pi}", va="center", ha="left", fontsize=7.5,
                        color="0.15", zorder=4,
                        fontweight="bold" if pi == lead[k] else "normal")

    top_ax.set_yscale("log")
    top_ax.set_ylim(thr, float(Zp.max()) * 2.6)
    top_ax.set_ylabel("outliers (>3 SD)")
    top_ax.spines["bottom"].set_visible(False)
    top_ax.tick_params(axis="x", length=0)
    top_ax.text(0.0, 1.12, f"{int(ext.sum())} patient values above {thr:g} SD, one colour per "
                f"patient, no colour repeated ({len(carriers)} of {Zp.shape[0]} patients). The largest belongs to a "
                f"different patient (bold) in {len(set(lead.tolist()))} of the {len(top)} genes.",
                transform=top_ax.transAxes, fontsize=11, color="0.25")

    for k in range(len(top)):
        for arr, off, col in [(Zp[:, k], -0.19, C_PAT), (Zh[:, k], 0.19, C_HC)]:
            v = ax.violinplot(np.clip(arr, ylo, thr), positions=[k + off], widths=0.34,
                              showextrema=False)
            v["bodies"][0].set_facecolor(col)
            v["bodies"][0].set_alpha(0.32)
            inside = arr <= thr
            jit = np.random.default_rng(SEED + k).normal(0, 0.03, len(arr))
            ax.plot(k + off + jit[inside], arr[inside], ".", ms=2.4, color=col, alpha=0.6)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(1.96, color="k", lw=0.7, ls=":")
    ax.text(len(top) - 0.4, 1.96, " 1.96 SD", va="center", ha="left", fontsize=10, color="0.4")
    ax.set_ylim(ylo, thr)
    ax.set_xlim(-0.6, len(top) - 0.4)
    ax.spines["top"].set_visible(False)
    ax.set_xticks(range(len(top)), sym, rotation=90)
    ax.set_xlabel(f"top {len(top)} group DE genes (sorted by $p_{{adj}}$), "
                  f"patients (red) vs healthy (grey), RUVg {DESIGN_SHORT[design]}")
    ax.set_ylabel("HC-standardized TMM-log2")

    for a, y in [(top_ax, 0), (ax, 1)]:
        a.plot([0, 1], [y, y], transform=a.transAxes, ls="none", marker=[(-1, -0.5), (1, 0.5)],
               ms=8, mec="k", mew=1, clip_on=False)
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
