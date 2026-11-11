import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import norm

from analysis.deg_fragility import DESIGNS, K_BINS, LABEL, MAX_STEPS, NM_CLASS, SCOPES, SHORT, run_genes
from viz_style import FS_ANNOT_SM

DESIGN_C = {"no_covariate": "#A2A2A2", "ruvg_k1": "#489ACA", "ruvg_k2": "#009E73", "ruvg_k3": "#CC79A7"}
K_C = {"1": "#D55E00", "2": "#E69F00", "3": "#009E73", "4-10": "#0072B2"}
C_PAT, C_HC = "#B2182B", "#4D4D4D"
KS = np.arange(0, 7)
LEGEND_OUT = dict(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)


def _grid(ax, axis="both"):
    ax.grid(axis=axis, alpha=0.3, lw=0.5)
    ax.set_axisbelow(True)


def _label(sym, gene):
    return sym if isinstance(sym, str) else gene.split(".")[0]


def plot_volcano(R, key, ax, shrunk_lim=3.0, y_max=50.0):
    """x = DESeq2 log2FC (symlog), DE points coloured by the apeGLM-shrunken log2FC: a large fold change
    that stays pale is one the data do not support beyond a few samples. -log10 p above y_max is drawn
    as a triangle on the top edge."""
    F = R["fits"][key]
    res = F["res"].dropna(subset=["pvalue"])
    x = res["log2FoldChange"]
    p = res["pvalue"]
    y = -np.log10(p.clip(lower=p[p > 0].min()))
    over = y > y_max
    y = y.clip(upper=y_max)
    sig = res.index.isin(F["sig"])
    ax.scatter(x[~sig], y[~sig], s=3, color="0.85", rasterized=True)
    shr = R["shrunk"][key].reindex(res.index)
    kw = dict(c=None, cmap="RdBu_r", vmin=-shrunk_lim, vmax=shrunk_lim, edgecolor="k", linewidth=0.3, zorder=3)
    for m, mk in [(sig & ~over, "o"), (sig & over, "^")]:
        sc = ax.scatter(x[m], y[m], s=18, marker=mk, **{**kw, "c": shr[m]})
    ax.axhline(-np.log10(F["pstar"]), color="k", ls=":", lw=0.8)
    for v in (-1, 1):
        ax.axvline(v, color="k", ls="--", lw=0.6)
    ax.set_xscale("symlog", linthresh=1)
    ax.set_ylim(0, y_max * 1.04)
    _grid(ax)
    ax.set_xlabel("log2 fold change (DESeq2; symlog)")
    ax.set_ylabel(f"-log10 p (Wald; capped at {y_max:g})")
    cb = plt.colorbar(sc, ax=ax, extend="both", fraction=0.05, pad=0.02)
    cb.set_label("shrunken log2FC (apeGLM)")


def plot_trajectory(R, key, ax):
    """One line per DE gene. Every point that is back above the BH cut after the line first dropped
    below it is marked x (gene named once): removing further top expressers of the same group can flip
    the contrast and make the gene significant again, usually in the opposite direction."""
    F = R["fits"][key]
    t = run_genes(R, key)
    lstar = -np.log10(F["pstar"])
    rng = np.random.default_rng(0)
    n_again = 0
    for g, zs in enumerate(R["traj"][key]):
        lp = -np.log10(2 * norm.sf(np.abs(zs))) - lstar
        xs = np.arange(len(zs)) + rng.uniform(-0.12, 0.12)
        ax.plot(xs, lp, "-o", ms=2.5, lw=0.8, color=K_C[t.k_bin.iloc[g]], alpha=0.55)
        below = np.where(lp < 0)[0]
        if not len(below):
            continue
        again = np.where((np.arange(len(lp)) > below[0]) & (lp >= 0))[0]
        if not len(again):
            continue
        n_again += 1
        ax.plot(xs[again], lp[again], "x", ms=6, mew=1.4, color="k", zorder=5)
        i = again[0]
        flip = "reversed" if np.sign(zs[i]) != np.sign(zs[0]) else "same sign"
        ax.annotate(f"{_label(t.symbol.iloc[g], t.gene.iloc[g])} ({flip})", (xs[i], lp[i]), xytext=(4, 2),
                    textcoords="offset points", fontsize=FS_ANNOT_SM, zorder=6)
    ax.axhline(0, color="k", lw=1)
    _grid(ax)
    ax.text(0.99, 0, "BH 5% cut ", fontsize=FS_ANNOT_SM, va="bottom", ha="right", transform=ax.get_yaxis_transform())
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xticks(range(MAX_STEPS + 1))
    ax.set_xlim(-0.4, max(len(z) for z in R["traj"][key]) - 0.6)
    ax.set_xlabel("highest-expressing samples removed (DE-direction group)")
    ax.set_ylabel("-log10 p above the BH cut\n(dispersion re-estimated, symlog)")
    for b in K_BINS:
        ax.plot([], [], color=K_C[b], label=f"lost after {b}")
    ax.plot([], [], "x", color="k", label=f"significant again after loss\n({n_again} genes)")
    ax.legend(**LEGEND_OUT)


def plot_cpm(R, key, per_row=30):
    """Every DE gene of one run (by padj, per_row per row); coloured points = samples whose removal
    ends significance, one colour per sample across the whole figure."""
    F = R["fits"][key]
    t = run_genes(R, key).sort_values("padj")
    Y = np.log2(R["cpm"][key[0]][list(t.gene)].values + 1)
    pat = F["is_case"]
    removed = []
    for _, r in t.iterrows():
        o = R["frag"][key][F["sig"].index(r.gene)]["order"]
        removed.append(o[:int(r.k_reestimated_disp)] if np.isfinite(r.k_reestimated_disp) else o)
    carriers = sorted(set().union(*map(set, removed)))
    shuffled = np.random.default_rng(0).permutation(len(carriers))
    col_of = {c: plt.get_cmap("Spectral")(shuffled[i] / max(len(carriers) - 1, 1)) for i, c in enumerate(carriers)}
    pat_rank = {i: n for n, i in enumerate(np.where(pat)[0])}
    hc_rank = {i: n for n, i in enumerate(np.where(~pat)[0])}
    name = lambda i: f"P{pat_rank[i]}" if pat[i] else f"C{hc_rank[i]}"
    labels = [f'{_label(a, g)} ({int(b) if np.isfinite(b) else ">10"})'
              for a, g, b in zip(t.symbol, t.gene, t.k_reestimated_disp)]
    n_row = int(np.ceil(len(t) / per_row))
    n_col = min(len(t), per_row)
    fig, axes = plt.subplots(n_row, 1, figsize=(max(5, 0.42 * n_col + 1.5), 4.2 * n_row), sharey=True,
                             constrained_layout=True, squeeze=False)
    rng = np.random.default_rng(0)
    for row, ax in enumerate(axes[:, 0]):
        ks = range(row * per_row, min((row + 1) * per_row, len(t)))
        for x, k in enumerate(ks):
            for off, m, col in [(-0.2, pat, C_PAT), (0.2, ~pat, C_HC)]:
                v = Y[m, k]
                if np.ptp(v) > 0:
                    vp = ax.violinplot(v, positions=[x + off], widths=0.36, showextrema=False)
                    vp["bodies"][0].set_facecolor(col)
                    vp["bodies"][0].set_alpha(0.25)
                ax.plot(x + off + rng.normal(0, 0.03, m.sum()), v, ".", ms=2, color=col, alpha=0.5)
            for i in removed[k]:
                off = -0.2 if pat[i] else 0.2
                ax.plot(x + off, Y[i, k], "o", ms=7, color=col_of[i], mec="k", mew=0.8, zorder=5)
                ax.text(x + off + 0.08, Y[i, k], name(i), fontsize=FS_ANNOT_SM, va="center", zorder=6)
        ax.set_xticks(range(len(ks)), [labels[k] for k in ks], rotation=90)
        ax.set_xlim(-0.6, n_col - 0.4)
        ax.set_ylabel("log2(CPM + 1)")
        _grid(ax, "y")
    axes[-1, 0].set_xlabel("DE genes by padj (samples removed until significance is lost)")
    table = pd.DataFrame(dict(gene=t.gene.values, symbol=t.symbol.values,
                              removed=[", ".join(name(i) + ":" + F["samples"][i] for i in r) for r in removed]))
    return fig, table


def plot_all_conditions(R):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4), constrained_layout=True)
    count_rows = []
    for s, ls in zip(SCOPES, ["-o", "--s"]):
        for d in DESIGNS:
            if (s, d) not in R["frag"]:
                continue
            t = run_genes(R, (s, d))
            rd = [int((t.k_reestimated_disp.fillna(np.inf) > k).sum()) for k in KS]
            fx = [int((t.k_fixed_disp.fillna(np.inf) > k).sum()) for k in KS]
            axes[0].plot(KS, rd, ls, ms=3, color=DESIGN_C[d], label=f"{LABEL[s]}, RUVg {SHORT[d]}")
            count_rows += [dict(scope=LABEL[s], design=SHORT[d], removed=k, deg_fixed_disp=a, deg_reestimated_disp=b)
                           for k, a, b in zip(KS, fx, rd)]
    axes[0].set_xticks(KS)
    _grid(axes[0], "y")
    axes[0].set_xlabel("highest-expressing samples removed per gene\n(dispersion re-estimated)")
    axes[0].set_ylabel("DE genes still significant")
    axes[0].legend(frameon=False, loc="upper right", borderaxespad=0)

    g = R["genes"]
    tab = g.groupby(["scope", "design", "k_bin"], sort=False).size().unstack(fill_value=0).reindex(
        columns=K_BINS, fill_value=0)
    tab = tab.loc[[(LABEL[s], SHORT[d]) for s in SCOPES for d in DESIGNS if (LABEL[s], SHORT[d]) in tab.index]]
    y = np.arange(len(tab))[::-1]
    left = np.zeros(len(tab))
    ax = axes[1]
    for b in K_BINS:
        ax.barh(y, tab[b], left=left, color=K_C[b], edgecolor="k", linewidth=0.3, height=0.7, label=f"{b} sample(s)")
        for yi, (v, lf, tot) in enumerate(zip(tab[b], left, tab.sum(axis=1))):
            if v / tot >= 0.15 and v >= 20:
                ax.text(lf + v / 2, y[yi], f"{v} ({v / tot:.0%})", ha="center", va="center", fontsize=10)
        left += tab[b].values
    ax.set_yticks(y, [f"{s}, {d}" for s, d in tab.index])
    ax.set_xlabel("DE genes")
    _grid(ax, "x")
    ax.legend(frameon=False, title="significance lost after", loc="upper right", borderaxespad=0)
    return fig, pd.DataFrame(count_rows), tab.assign(total=tab.sum(axis=1))


def plot_driver_bias_nm(drivers, nm_meta, assoc, classed, cut=3):
    fig, axes = plt.subplots(4, 1, figsize=(7, 15), height_ratios=[1, 0.45, 1, 1], constrained_layout=True)
    pat = drivers[drivers.first_removed_group == "patient"]
    dis = nm_meta[nm_meta.Z == "disease"]
    ax = axes[0]
    d = pat.drop_duplicates("first_removed")
    ax.hist(d.bias_z_pool, bins=np.arange(0, d.bias_z_pool.max() + 0.5, 0.25), color=C_PAT, alpha=0.7, edgecolor="k",
            linewidth=0.3)
    ax.axvline(cut, color="k", ls=":", lw=0.8)
    n_ext = int((d.bias_z_pool > cut).sum())
    ax.text(0.98, 0.98, f"BIAS-extreme (|z| > {cut}): {n_ext}/{len(d)} = {n_ext / len(d):.0%}\n"
            f"all disease samples: {(dis.bias_z_pool > cut).mean():.0%}", transform=ax.transAxes, ha="right", va="top",
            fontsize=FS_ANNOT_SM)
    ax.set(xlabel="BIAS max |robust z| vs the NM HC pool (n=592)", ylabel="Patients driving a DEG",
           title="A. Most patients driving DEGs are not BIAS-extreme")
    _grid(ax, "y")

    ax = axes[1]
    labels = ["|log2FC| <= 10", "|log2FC| > 10"]
    for i, b in enumerate([False, True]):
        v = pat.loc[pat.boundary_lfc == b, "nm_z_de_dir"].dropna()
        n_kept, n_lost = int((v > 1.96).sum()), int((v <= 1.96).sum())
        ax.barh(i, n_kept, color=C_PAT, alpha=0.7, edgecolor="k", linewidth=0.3)
        ax.barh(i, n_lost, left=n_kept, color="white", edgecolor="k", linewidth=0.3, hatch="//")
        ax.text(n_kept + n_lost + 3, i, f"{n_kept} kept / {n_lost} normal in NM", va="center", fontsize=FS_ANNOT_SM)
    ax.legend(handles=[Patch(facecolor=C_PAT, alpha=0.7, edgecolor="k", label="still a deviation in NM (Z > 1.96)"),
                       Patch(facecolor="white", edgecolor="k", hatch="//", label="normal in NM (Z <= 1.96)")],
              frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=2, fontsize=FS_ANNOT_SM)
    ax.set_yticks([0, 1], labels)
    ax.invert_yaxis()
    ax.set_xlim(right=ax.get_xlim()[1] * 1.35)
    ax.set(xlabel="DEGs (judged at the driving patient)", title="B. Part of the DEG signal is normal under NM")

    ax = axes[2]
    edges = [0, 1, 2, 3, 5, np.inf]
    names = ["<1", "1-2", "2-3", "3-5", ">5"]
    y = 100 * (dis.frac_z196 - dis.groupby("batch").frac_z196.transform("median"))
    b = pd.cut(dis.bias_z_pool, edges, labels=names, right=False)
    v = [y[b == n].dropna() for n in names]
    bp = ax.boxplot(v, positions=range(len(names)), widths=0.5, showfliers=False, patch_artist=True,
                    medianprops=dict(color="k"))
    for box in bp["boxes"]:
        box.set(facecolor=C_PAT, alpha=0.5)
    for x, u in enumerate(v):
        ax.text(x, 0.02, f"n={len(u)}", transform=ax.get_xaxis_transform(), ha="center", fontsize=FS_ANNOT_SM)
    rho = assoc.loc[(assoc.Z == "disease") & (assoc.outcome == "frac_z196"), "rho_within_batch"].iloc[0]
    ax.text(0.02, 0.98, f"within-batch Spearman rho = {rho:.2f}", transform=ax.transAxes, va="top", fontsize=FS_ANNOT_SM)
    ax.axhline(0, color="k", ls=":", lw=0.8)
    ax.set_xticks(range(len(names)), names)
    ax.set(xlabel="BIAS max |robust z| vs the NM HC pool (n=592)",
           ylabel="NM hit rate minus batch median (%p)\n(% of genes with |Z| > 1.96)",
           title="C. NM hit rate does not rise with BIAS (disease samples)")
    _grid(ax, "y")

    ax = axes[3]
    rng = np.random.default_rng(0)
    for i, c in enumerate(NM_CLASS):
        v = classed[classed.nm_class == c]
        ax.boxplot(v.bias_z_pool, positions=[i], widths=0.5, showfliers=False, medianprops=dict(color="k"))
        ax.scatter(i + rng.normal(0, 0.08, len(v)), v.bias_z_pool, s=10, alpha=0.5, color=C_PAT)
    ax.axhline(cut, color="k", ls=":", lw=0.8)
    ax.set_xticks(range(len(NM_CLASS)), [f"{c.replace(' (', chr(10) + '(')}\n{(classed.nm_class == c).sum()} DEGs\n"
                                         f"{classed.loc[classed.nm_class == c, 'first_removed'].nunique()} patients"
                                         for c in NM_CLASS])
    ax.set(ylabel="BIAS max |robust z| vs the NM HC pool (n=592)",
           title="D. BIAS of the driving patient by NM call", xlabel="NM Z of the driving patient at its DE gene")
    _grid(ax, "y")
    return fig
