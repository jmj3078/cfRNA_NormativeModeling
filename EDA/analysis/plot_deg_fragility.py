import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

from analysis.deg_fragility import DESIGNS, K_BINS, LABEL, MAX_STEPS, SCOPES, SHORT, run_genes
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
