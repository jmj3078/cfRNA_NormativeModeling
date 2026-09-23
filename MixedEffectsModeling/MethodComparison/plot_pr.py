"""Figures for the in-silico method comparison.

Reads the aggregated cache written by sweep.py and draws three things:

  null_falsecalls  -- how many genes each method calls abnormal in a healthy sample
  null_panel       -- null p-value calibration and post-BH rejection rate vs level
  plot             -- precision-recall grid, rows = injected log2FC, columns = split

All aggregation is cached in cache/plot_data.pkl, so redrawing never re-reads the
several-GB matrix cache. The companion notebook 4_in_silico_method_comparison.ipynb
loads the same cache and only draws.

Run:  python plot_pr.py [--refresh]
"""
import pickle
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from viz_style import apply_style

from common import CACHE, Q_LEVELS
from sweep import MATRIX_DIR, SPLITS, fold_spec

apply_style()

ARMS = ["engine", "autoencoder", "pca", "outsingle", "peer"]
LABELS = {"engine": "Normative engine", "autoencoder": "OUTRIDER (AE)", "pca": "PCA",
          "outsingle": "OutSingle", "peer": "PEER"}
COLORS = {"engine": "#000000", "autoencoder": "#1B9E77", "pca": "#7570B3",
          "outsingle": "#D95F02", "peer": "#E7298A"}
Q_MARKERS = {0.05: "o", 0.10: "s", 0.20: "^", 0.25: "D"}
# only the normative engine is solid: the baselines read as one dashed family against it
STYLES = {"engine": "-", "autoencoder": (0, (5, 2)), "pca": (0, (3, 1.5)),
          "outsingle": (0, (1.5, 1.5)), "peer": (0, (6, 1.5, 1.5, 1.5))}
SPLIT_TITLE = {"cv": "5-fold CV", "lobo": "leave-one-batch-out"}
LOG2FCS_SHOWN = [2.0, 3.0, 4.0, 5.0]
RECALL_GRID = np.linspace(0.0, 1.0, 201)
# BH levels swept to trace the realized operating curve; the reported q values are
# marked on it, so curve and markers are the same procedure.
Q_GRID = np.concatenate([[0.0], np.logspace(-6, 0, 300)])


def bh_counts_vs_q(p, labels, q_grid):
    """tp/fp/fn after per-sample BH at each q in q_grid -- the same procedure every number
    in the sweep table uses, so the plotted curve and the marked operating points are one
    thing rather than two.

    Per sample, BH rejects the k largest-ranked p with k = max{i : p_(i) <= q i / m}.
    Writing c_i = p_(i) m / i and g_i = min_{j>=i} c_j (non-decreasing), that is
    k(q) = searchsorted(g, q), so one sort per sample answers the whole q grid."""
    tp = np.zeros(len(q_grid)); fp = np.zeros(len(q_grid))
    n_pos = 0
    for row in p:
        ok = np.isfinite(row)
        pv, lab = row[ok], labels[ok]
        m = len(pv)
        if m == 0:
            continue
        order = np.argsort(pv, kind="stable")
        ps, ls = pv[order], lab[order]
        c = ps * m / np.arange(1, m + 1)
        g = np.minimum.accumulate(c[::-1])[::-1]
        k = np.searchsorted(g, q_grid, side="right")
        cum_tp = np.concatenate([[0], np.cumsum(ls)])
        tp += cum_tp[k]
        fp += k - cum_tp[k]
        n_pos += int(ls.sum())
    return tp, fp, n_pos


def collect(split):
    """Per-sample BH operating curves, per fold and pooled.

    Both are kept because they answer different questions and, in LOBO, disagree. Pooling
    sums tp/fp across folds, so it is a ratio-of-sums dominated by the largest held-out
    batch (n = 116 vs 27 here); the per-fold median weights every unseen batch equally.
    A ratio of sums need not lie between the per-fold ratios, so the pooled curve can sit
    outside the fold IQR entirely -- a Simpson-type effect, not an error."""
    acc, per_fold = {}, {}
    for _, tag, _, _, _ in fold_spec(split):
        for f in LOG2FCS_SHOWN:
            d = np.load(MATRIX_DIR / f"{split}_{tag}_log2fc{f:g}.npz", allow_pickle=True)
            labels = d["labels"]
            for arm in ARMS:
                tp, fp, n_pos = bh_counts_vs_q(d[f"p__{arm}"], labels, Q_GRID)
                a = acc.setdefault((f, arm), [np.zeros(len(Q_GRID)), np.zeros(len(Q_GRID)), 0])
                a[0] += tp; a[1] += fp; a[2] += n_pos
                rej = tp + fp
                ok = rej > 0                      # q so small nothing is rejected: precision undefined
                rec, prec = tp[ok] / max(n_pos, 1), tp[ok] / rej[ok]
                curve = np.interp(RECALL_GRID, rec, prec)
                curve[(RECALL_GRID < rec[0]) | (RECALL_GRID > rec[-1])] = np.nan
                per_fold.setdefault((f, arm), []).append(curve)
            d.close()
        print(f"  {split}/{tag} curves done", flush=True)

    curves = {}
    for k, v in acc.items():
        ok = (v[0] + v[1]) > 0
        curves[k] = (v[0][ok] / max(v[2], 1), v[0][ok] / (v[0] + v[1])[ok])
    bands = {k: np.nanpercentile(np.vstack(v), [25, 50, 75], axis=0)
             for k, v in per_fold.items()}
    return curves, bands


def operating_points(split):
    """(recall, precision) of per-sample BH at each reported q, both pooled over folds and
    as the per-fold median -- read from the same counts the summary table is built from,
    so figure and table can be checked against each other."""
    df = pd.read_csv(CACHE / "sweep_nm_counts.csv")
    df = df[df.split == split]
    pooled = df.groupby(["arm", "q", "log2fc"])[["tp", "fp", "fn"]].sum()
    rp = lambda d: (d["tp"] / max(d["tp"] + d["fn"], 1), d["tp"] / max(d["tp"] + d["fp"], 1))
    out = {}
    for a in ARMS:
        for q in Q_LEVELS:
            for f in LOG2FCS_SHOWN:
                out[(q, f, a)] = rp(pooled.loc[(a, q, f)])
                sub = df[(df.arm == a) & (df.q == q) & (df.log2fc == f)]
                per = np.array([rp(r) for _, r in sub.iterrows()])
                out[("med", q, f, a)] = tuple(np.median(per, axis=0))
    return out


def plot(d):
    """Rows = injected effect size, columns = split. The BH operating curve does not depend
    on q -- only the operating point on it does -- so q is a marker shape rather than a row,
    which buys 4x the panel area the old q-by-log2FC grid wasted on redrawing one curve.

    The line is the per-fold MEDIAN, so it lies inside its own ribbon by construction and
    every unseen batch counts once; each marker is placed at the median recall reached at
    that BH level and read off the curve, so it lies on the line rather than beside it.
    One line per method, never two: the engine is the only solid curve and the four
    baselines are a dashed family, so the comparison reads at a glance. The pooled
    (count-weighted) curve is left out of the figure -- it is what sweep_nm_summary.csv
    reports and is kept in the cache as d[split]["curves"] for checking."""
    fig, axes = plt.subplots(len(LOG2FCS_SHOWN), len(SPLITS),
                             figsize=(2.6 * len(SPLITS), 2.0 * len(LOG2FCS_SHOWN)),
                             sharex=True, sharey=True, squeeze=False)
    for i, f in enumerate(LOG2FCS_SHOWN):
        for j, split in enumerate(SPLITS):
            ax = axes[i, j]
            bands, ops = d[split]["bands"], d[split]["ops"]
            for arm in ARMS:
                lo, med, hi = bands[(f, arm)]
                ax.fill_between(RECALL_GRID, lo, hi, color=COLORS[arm], alpha=0.15, lw=0)
                ax.plot(RECALL_GRID, med, color=COLORS[arm], ls=STYLES[arm],
                        lw=2.0 if arm == "engine" else 1.4,
                        label=LABELS[arm] if (i, j) == (0, 0) else None,
                        zorder=4 if arm == "engine" else 3)
                for q in Q_LEVELS:
                    r = ops[("med", q, f, arm)][0]
                    ax.plot(r, np.interp(r, RECALL_GRID, med), marker=Q_MARKERS[q], ms=4,
                            color=COLORS[arm], mec="white", mew=0.6, ls="none", zorder=5)
            ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
            ax.grid(alpha=0.2, lw=0.5)
            if i == 0:
                ax.set_title(SPLIT_TITLE[split])
            if j == len(SPLITS) - 1:
                ax.text(1.06, 0.5, f"injected log2FC = {f:g}", rotation=270, va="center",
                        ha="left", fontsize=9, transform=ax.transAxes,
                        bbox=dict(facecolor="0.92", edgecolor="none", pad=3))
            if j == 0:
                ax.set_ylabel("Precision")
            if i == len(LOG2FCS_SHOWN) - 1:
                ax.set_xlabel("Recall")

    h, l = axes[0, 0].get_legend_handles_labels()
    leg = fig.legend(h, l, loc="upper left", bbox_to_anchor=(1.02, 0.95), frameon=False,
                     fontsize=7, title="Method")
    leg.get_title().set_fontsize(7)
    qh = [plt.Line2D([], [], marker=m, ls="none", color="0.3", ms=4.5, mec="white", mew=0.7)
          for m in Q_MARKERS.values()]
    leg2 = fig.legend(qh, [f"{q}" for q in Q_MARKERS], loc="upper left",
                      bbox_to_anchor=(1.02, 0.55), frameon=False, fontsize=7, title="BH q")
    leg2.get_title().set_fontsize(7)
    fig.suptitle("Normative-mode detection after per-sample BH"
                 "\nline = per-fold median, ribbon = fold IQR, markers = BH level",
                 fontsize=8)
    fig.tight_layout()
    out = CACHE.parent / "Figures" / "pr_grid.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}", flush=True)


def collect_pcal():
    """Empirical CDF of the null p-values on Q_GRID, pooled over folds. Under a valid test
    every gene at log2FC = 0 is a true null, so P(p <= a) must equal a: the identity line is
    the target, and sitting above it means the p-values are anti-conservative -- any FDR
    claim built on them is void before BH ever runs."""
    out = {}
    for split in SPLITS:
        for arm in ARMS:
            hit, n = np.zeros(len(Q_GRID)), 0
            for _, tag, _, _, _ in fold_spec(split):
                d = np.load(MATRIX_DIR / f"{split}_{tag}_log2fc0.npz", allow_pickle=True)
                pv = d[f"p__{arm}"].ravel()
                pv = np.sort(pv[np.isfinite(pv)])
                hit += np.searchsorted(pv, Q_GRID, side="right"); n += pv.size
                d.close()
            out[(split, arm)] = hit / n
    return out


def collect_null():
    """Pooled null rejection rate vs Q_GRID for every (split, arm), from the log2FC = 0 cells."""
    out = {}
    for split in SPLITS:
        for arm in ARMS:
            rej, tot = np.zeros(len(Q_GRID)), 0
            for _, tag, _, _, _ in fold_spec(split):
                d = np.load(MATRIX_DIR / f"{split}_{tag}_log2fc0.npz", allow_pickle=True)
                pv = d[f"p__{arm}"]
                tp, fp, _ = bh_counts_vs_q(pv, d["labels"], Q_GRID)
                rej += tp + fp; tot += int(np.isfinite(pv).sum())
                d.close()
            out[(split, arm)] = rej / tot
    return out


def null_headline(q=0.05):
    """The false-positive story as one number per method: in a healthy sample where nothing
    is wrong, how many genes does each method call abnormal at BH q?

    Counts, not rates, because the count is the thing a reader can judge -- and the exact
    reference is free: with no multiple-testing correction at all you would expect q * G
    false calls, so any bar past that line is rejecting more AFTER BH than an uncorrected
    test would."""
    d = plot_data()
    n_genes = len((CACHE / "eval_universe_nm.txt").read_text().split())
    lg = np.log10(np.maximum(Q_GRID, 1e-12))
    rate = lambda sp, a: np.interp(np.log10(q), lg, d["null"][(sp, a)])

    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    h = 0.38
    xmax = max(rate(sp, a) for sp in ("cv", "lobo") for a in ARMS) * n_genes * 1.12
    for k, split in enumerate(SPLITS):
        y = np.arange(len(ARMS)) + (k - 0.5) * h
        val = np.array([rate(split, a) * n_genes for a in ARMS])
        ax.barh(y, val, height=h, color=["0.15" if split == "cv" else "0.55"][0],
                edgecolor="none", label=SPLIT_TITLE[split])
        for yi, v in zip(y, val):
            ax.text(v + xmax * 0.01, yi, f"{v:,.0f}", va="center", fontsize=8)
    ax.axvline(q * n_genes, color="#D95F02", lw=1.2, ls=(0, (3, 2)))
    ax.annotate(f"no multiple-testing\ncorrection at all\n({q:g} x {n_genes:,} = {q * n_genes:,.0f})",
                xy=(q * n_genes, -0.62), color="#D95F02", fontsize=7,
                ha="center", va="bottom", annotation_clip=False)
    ax.set_xlim(0, xmax); ax.set_ylim(len(ARMS) - 0.5, -1.35)
    ax.set_yticks(range(len(ARMS))); ax.set_yticklabels([LABELS[a] for a in ARMS])
    ax.set_xlabel(f"false gene calls per healthy sample  (BH q = {q:g})")
    ax.grid(axis="x", alpha=0.25, lw=0.5)
    ax.legend(frameon=False, fontsize=7, loc="upper right", ncol=1, bbox_to_anchor=(0.995, 1.0))
    ax.set_title("What each method does when nothing is wrong\n"
                 f"log2FC = 0, no signal injected; {n_genes:,} genes scored per sample", fontsize=8)
    fig.tight_layout()
    out = CACHE.parent / "Figures" / "null_falsecalls.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}", flush=True)


def null_panel():
    """Two rows, both from the log2FC = 0 cells where nothing was injected.

    Top: null p-value calibration. Every gene is a true null, so the observed fraction of
    p <= a must equal a -- the identity line is what a correct test follows. Above it means
    anti-conservative p-values.

    Bottom: what per-sample BH then does at level q. Under a complete null BH controls the
    family-wise rate, so it should reject essentially nothing; y = q is drawn as the
    no-correction reference, and a curve above it is rejecting more after BH than an
    uncorrected test would."""
    d = plot_data()
    pcal, null = d["pcal"], d["null"]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.4), sharex=True, sharey="row")
    for j, split in enumerate(SPLITS):
        for arm in ARMS:
            axes[0, j].plot(Q_GRID, np.maximum(pcal[(split, arm)], 1e-6), color=COLORS[arm],
                            lw=1.5, label=LABELS[arm] if j == 0 else None)
            axes[1, j].plot(Q_GRID, np.maximum(null[(split, arm)], 1e-6), color=COLORS[arm], lw=1.5)
        for i in (0, 1):
            ax = axes[i, j]
            ax.plot(Q_GRID, Q_GRID, color="0.4", lw=1.0, ls=(0, (2, 2)))
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlim(1e-4, 1); ax.set_ylim(1e-5, 1.5)
            ax.grid(alpha=0.2, lw=0.5)
        axes[0, j].set_title(SPLIT_TITLE[split])
        axes[1, j].set_xlabel("nominal level (alpha, or BH q)")
    axes[0, 0].set_ylabel("observed P(p <= alpha)\n[raw p-value calibration]")
    axes[1, 0].set_ylabel("rejection rate after\nper-sample BH at q")
    axes[0, 1].text(0.50, 0.42, "target: y = alpha", color="0.35", fontsize=6.5,
                    rotation=39, transform=axes[0, 1].transAxes)
    axes[1, 1].text(0.50, 0.42, "no correction (y = q)", color="0.35", fontsize=6.5,
                    rotation=39, transform=axes[1, 1].transAxes)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, frameon=False, fontsize=7,
               bbox_to_anchor=(0.5, -0.03))
    fig.suptitle("Calibration under a pure null (log2FC = 0, no signal injected)"
                 "\ndashed = the line a correctly calibrated test follows; above it = too many"
                 " false positives", fontsize=8)
    fig.tight_layout()
    out = CACHE.parent / "Figures" / "null_rejection.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}", flush=True)


PLOT_DATA = CACHE / "plot_data.pkl"


def plot_data(refresh=False):
    """Everything the figures need, already aggregated. Cached so redrawing (notebook,
    restyling) never re-reads the several-GB matrix cache."""
    if PLOT_DATA.exists() and not refresh:
        return pickle.loads(PLOT_DATA.read_bytes())
    d = {"null": collect_null(), "pcal": collect_pcal(), "q_grid": Q_GRID, "recall_grid": RECALL_GRID}
    for split in SPLITS:
        print(f"== {split} ==", flush=True)
        curves, bands = collect(split)
        d[split] = dict(curves=curves, bands=bands, ops=operating_points(split))
    PLOT_DATA.write_bytes(pickle.dumps(d))
    return d


def main():
    d = plot_data(refresh="--refresh" in sys.argv)
    null_headline()
    null_panel()
    plot(d)


if __name__ == "__main__":
    main()
