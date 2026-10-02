import json
import sys
from pathlib import Path

from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import gaussian_kde, mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config

LOBO_DIR = config.LOBO_MIXED_DIR

def _median_heuristic_gamma(pooled):
    d2 = squareform(pdist(pooled, metric="sqeuclidean"))
    med = np.median(d2[d2 > 0])
    return (1.0 / med if med > 0 else 1.0), d2


def _ref_centroid_direction(Z, is_hc):
    """Kernel-embedding distance to the held-out-HC mean embedding, in the SAME
    RKHS (same RBF kernel, same per-batch median-heuristic gamma) as the MMD
    test above -- not a raw Euclidean distance in gene space, which collapses
    19,804 genes with equal, uncorrelated weight and is only a Euclidean cousin
    of the mean|Z|/chi2 statistics already rejected for exactly that reason.
    The reference is built ONLY from this batch's own LOBO-scored held-out HC
    (a model that never trained on this batch), never from CV (CV's
    StratifiedKFold splits within every batch rather than excluding one, so it
    would leak this batch's own samples back into the reference). Each HC
    sample compares to the mean EMBEDDING of the OTHER held-out HC (leave-one-out,
    or it would be pulled toward its own position and its distance spuriously
    shrunk); disease samples never contribute to the reference, so they compare
    against the full HC mean embedding. Squared distance to a kernel mean
    embedding has a closed form in terms of the Gram matrix alone --
    ||phi(x) - mean_i phi(h_i)||^2 = k(x,x) - (2/n)sum_i k(x,h_i) + (1/n^2)sum_ij k(h_i,h_j)."""
    hc_Z, dis_Z = Z[is_hc], Z[~is_hc]
    n_hc = len(hc_Z)
    gamma, _ = _median_heuristic_gamma(np.vstack([hc_Z, dis_Z]))
    Khh = np.exp(-gamma * squareform(pdist(hc_Z, metric="sqeuclidean")))
    np.fill_diagonal(Khh, 1.0)
    Kdh = np.exp(-gamma * cdist(dis_Z, hc_Z, metric="sqeuclidean"))

    hc_row_sum = Khh.sum(axis=1) - 1.0  # sum_{k!=i} k(h_i,h_k), k(x,x)=1 excluded
    hc_total_full = Khh.sum()           # sum_{k,l} k(h_k,h_l) over ALL HC, diagonal included
    n = n_hc - 1
    # sum_{k!=i,l!=i} k(h_k,h_l): full sum minus row i, col i (each hc_row_sum+diag), plus the
    # doubly-removed (i,i)=1 term back once -- see derivation in the review that caught this.
    rest_total = hc_total_full - 2 * hc_row_sum - 1.0
    d2_hc = np.clip(1.0 - (2.0 / n) * hc_row_sum + rest_total / n**2, 0, None)
    d2_dis = np.clip(1.0 - (2.0 / n_hc) * Kdh.sum(axis=1) + hc_total_full / n_hc**2, 0, None)

    d_hc, d_dis = np.sqrt(d2_hc), np.sqrt(d2_dis)
    _, p_direction = mannwhitneyu(d_dis, d_hc, alternative="greater")
    return d_hc, d_dis, float(p_direction)


def p_to_asterisk(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def short_batch_label(batch):
    name, _, num = batch.replace(" et al.", "").partition("_Batch_")
    return f"{name}_{num}" if num else name


def plot_mmd_direction(df, raw, out_path, p_col="p_direction", sort_col="mmd2", row_h=1.4):
    """Per-batch kernel-embedding distance from the HC reference: held-out HC vs disease.
    Shared by 0_outrider_comparison and 2_lobo_validation."""
    d = df.sort_values(sort_col, ascending=False).reset_index(drop=True)
    fig, axes = plt.subplots(len(d), 1, figsize=(7, row_h * len(d)), sharex=True)
    if len(d) == 1:
        axes = [axes]

    color_hc, color_sig, color_ns = "#A4AFB8", "#00C78B", "#D64545"
    for ax, (_, row) in zip(axes, d.iterrows()):
        r = raw[row["batch"]]
        d_hc, d_dis = np.asarray(r["d_hc"]), np.asarray(r["d_dis"])

        x_min, x_max = min(d_hc.min(), d_dis.min()), max(d_hc.max(), d_dis.max())
        x_margin = (x_max - x_min) * 0.2
        x_grid = np.linspace(x_min - x_margin, x_max + x_margin, 300)
        kde_hc, kde_dis = gaussian_kde(d_hc)(x_grid), gaussian_kde(d_dis)(x_grid)

        mean_hc, mean_dis = d_hc.mean(), d_dis.mean()
        p_val = row.get(p_col, 1.0)
        color_dis = color_sig if p_val < 0.05 else color_ns

        ax.plot(x_grid, kde_hc, color=color_hc, lw=1.5)
        ax.fill_between(x_grid, kde_hc, color=color_hc, alpha=0.35)
        ax.plot(x_grid, kde_dis, color=color_dis, lw=1.5)
        ax.fill_between(x_grid, kde_dis, color=color_dis, alpha=0.35)

        max_y = max(kde_hc.max(), kde_dis.max())
        ax.axvline(mean_hc, color=color_hc, linestyle="--", lw=1.5, zorder=3)
        ax.axvline(mean_dis, color=color_dis, linestyle="--", lw=1.5, zorder=3)

        y_bar = max_y * 1.15
        ax.annotate("", xy=(mean_hc, y_bar), xytext=(mean_dis, y_bar),
                    arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))

        delta = mean_dis - mean_hc
        ax.text((mean_hc + mean_dis) / 2, y_bar + max_y * 0.08,
                f"delta={delta:+.3f} ({p_to_asterisk(p_val) or 'n.s.'})",
                ha="center", va="bottom", fontweight="bold",
                color=color_ns if delta < 0 else "black")

        ax.set_ylabel(short_batch_label(row["batch"]), rotation=0, ha="right", va="center")
        ax.set_ylim(0, max_y * 1.55)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="x", linestyle=":", alpha=0.4)

    axes[-1].set_xlabel("Kernel embedding distance from HC reference")
    axes[0].legend(handles=[
        Patch(facecolor=color_hc, edgecolor=color_hc, alpha=0.5, label="Held-out HC"),
        Patch(facecolor=color_sig, edgecolor=color_sig, alpha=0.5, label="Disease (sig.)"),
        Patch(facecolor=color_ns, edgecolor=color_ns, alpha=0.5, label="Disease (n.s.)"),
    ], loc="upper right", frameon=False)
    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    return fig


# ---------------------------------------------------------------------------
# Shared LOBO evaluation (2_lobo_validation and 0_outrider_comparison use the
# SAME criterion through these functions -- the tables used to be duplicated in
# each notebook, which is how they drifted apart).
# ---------------------------------------------------------------------------
FDR_Q, Z_CAP = 0.05, 10.0
MIN_N_HC_EVAL, MIN_N_DIS_EVAL = 25, 10
DISEASE_SUBSET = {"Moore et al._Batch_1": ("PDAC",)}
_COND = None


def condition_of(names):
    """Stage/Condition for the given sample names."""
    global _COND
    if _COND is None:
        import scanpy as sc
        _COND = sc.read_h5ad(config.H5AD_PATH, backed="r").obs["Stage/Condition"].astype(str)
    return _COND.loc[list(names)].values


def subset_keep(batch, names, is_hc):
    """Drop disease samples outside DISEASE_SUBSET[batch]; HC is never dropped.
    Moore_1 pools Pancreatitis/PDAC/IPMN/Islet/Other, and pooling them is what
    broke the two-sample tests there -- see memory project_lobo_moore1_pdac_only."""
    want = DISEASE_SUBSET.get(batch)
    if want is None:
        return np.ones(len(names), bool)
    return np.asarray(is_hc) | np.isin(condition_of(names), want)


def bh_q(p):
    p = np.asarray(p, dtype=float)
    m = len(p)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.clip(q, 0, 1)
    return out


def load_engine_arm(shash=True, ood_filter=True, lobo_dir=None):
    """{batch: dict(names, is_hc, Z)} from the LOBO fits, OOD-filtered and
    DISEASE_SUBSET-restricted. This is the sample set every arm must score on."""
    lobo_dir = lobo_dir or LOBO_DIR
    out = {}
    for path in sorted(lobo_dir.glob("*/meta.json")):
        meta = json.load(open(path))
        b = meta["batch_id"]
        Z = np.load(path.parent / ("Z_test_shash.npy" if shash else "Z_test.npy")).astype(np.float32)
        names, is_hc = np.array(meta["test_names"]), np.array(meta["test_is_hc"])
        if ood_filter:
            keep = np.load(path.parent / "ood_mask.npy")
            Z, names, is_hc = Z[keep], names[keep], is_hc[keep]
        sub = subset_keep(b, names, is_hc)
        out[b] = dict(names=names[sub], is_hc=is_hc[sub],
                      Z=np.nan_to_num(Z[sub], nan=0.0, posinf=Z_CAP, neginf=-Z_CAP))
    return out


def load_csv_arm(z_dir, ref, pattern="z_test_{safe}.csv"):
    """Same sample set as `ref` (an arm dict), Z read from per-batch CSVs --
    OUTRIDER's held-out Z. Reindexing to ref's names is what makes the arms
    comparable; the gene sets differ by construction and are reported."""
    out = {}
    for b, r in ref.items():
        f = z_dir / pattern.format(safe=b.replace(" ", "_"))
        if not f.exists():
            continue
        df = pd.read_csv(f, index_col=0).loc[list(r["names"])]
        out[b] = dict(names=r["names"], is_hc=r["is_hc"],
                      Z=np.nan_to_num(df.values, nan=0.0, posinf=Z_CAP, neginf=-Z_CAP))
    return out


def _scoreable(arm):
    for b, a in sorted(arm.items()):
        n_hc, n_dis = int(a["is_hc"].sum()), int((~a["is_hc"]).sum())
        if n_hc >= MIN_N_HC_EVAL and n_dis >= MIN_N_DIS_EVAL:
            yield b, a, n_hc, n_dis


def rkhs_stats(arm):
    """RKHS kernel-embedding distance to the held-out-HC mean embedding,
    disease vs HC by one-sided Mann-Whitney, BH across batches."""
    rows, raw = [], {}
    for b, a, n_hc, n_dis in _scoreable(arm):
        d_hc, d_dis, _ = _ref_centroid_direction(a["Z"], a["is_hc"])
        u, pu = mannwhitneyu(d_dis, d_hc, alternative="greater")
        rows.append(dict(batch=b, n_gene=a["Z"].shape[1], n_hc=n_hc, n_dis=n_dis,
                         auc=u / (n_dis * n_hc), mw_p=pu))
        raw[b] = dict(d_hc=d_hc, d_dis=d_dis)
    stats = pd.DataFrame(rows).sort_values("auc", ascending=False)
    stats["bh_q"] = bh_q(stats["mw_p"].values)
    return stats, raw


def nsig_table(arm):
    """Per-sample count of BH-significant genes."""
    from MixedEffectsModeling.core.calibration import bh_fdr_reject
    from scipy.stats import norm
    rows = []
    for b, a in sorted(arm.items()):
        P = 2 * norm.sf(np.abs(np.clip(a["Z"], -Z_CAP, Z_CAP)))
        ns = [int(bh_fdr_reject(r, q=FDR_Q).sum()) for r in P]
        rows += [(b, nm, bool(hc), v) for nm, hc, v in zip(a["names"], a["is_hc"], ns)]
    return pd.DataFrame(rows, columns=["batch", "sample", "is_hc", "n_sig"])


def nsig_stats(nsig):
    rows = []
    for b, g in nsig.groupby("batch"):
        h = g.loc[g["is_hc"], "n_sig"].values
        d = g.loc[~g["is_hc"], "n_sig"].values
        if len(h) < MIN_N_HC_EVAL or len(d) < MIN_N_DIS_EVAL:
            continue
        u, pu = mannwhitneyu(d, h, alternative="greater")
        rows.append(dict(batch=b, n_hc=len(h), n_dis=len(d), hc_med=np.median(h),
                         dis_med=np.median(d), auc=u / (len(d) * len(h)), mw_p=pu))
    stats = pd.DataFrame(rows).sort_values("auc", ascending=False)
    stats["bh_q"] = bh_q(stats["mw_p"].values)
    return stats


def plot_nsig_box(nsig, stats, title, out_path):
    GROUPS = ['Disease', 'HC (held-out)']
    order = stats.sort_values('dis_med', ascending=False)['batch'].tolist()
    plot = nsig[nsig['batch'].isin(order)].copy()
    plot['group'] = np.where(plot['is_hc'], GROUPS[1], GROUPS[0])
    plot['x'] = np.log10(plot['n_sig'] + 1)
    q_of = dict(zip(stats['batch'], stats['bh_q']))

    fig, ax = plt.subplots(figsize=(7.5, 0.45 * len(order) + 1.8))
    sns.boxplot(data=plot, y='batch', x='x', hue='group', order=order, hue_order=GROUPS,
                   palette={GROUPS[0]: "#be3b3b", GROUPS[1]: '#c8cdd2'}, orient='h', linecolor='black',
                   linewidth=0.8, width=0.65, ax=ax, showfliers=False)

    TICKS = [0, 1, 3, 10, 30, 100, 300, 1000, 10000]
    ax.set_xticks(np.log10(np.array(TICKS) + 1.0))
    ax.set_xticklabels(TICKS)
    x_max = plot['x'].max()
    ax.set_xlim(-0.12, x_max + 0.75)

    for i, b in enumerate(order):
        x_br = plot['x'].max() + 0.16
        ax.plot([x_br, x_br + 0.07, x_br + 0.07, x_br], [i - 0.21, i - 0.21, i + 0.21, i + 0.21],
                color='black', lw=0.9, clip_on=False)
        ax.text(x_br + 0.12, i, p_to_asterisk(q_of[b]) or 'ns', ha='left', va='center',
                 clip_on=False)

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([short_batch_label(b) for b in order])
    ax.set_xlabel('BH-significant genes per sample (q=0.05)')
    ax.set_ylabel('Batch (held-out)')
    ax.set_title(title, loc='left')
    h, l = ax.get_legend_handles_labels()
    ax.get_legend().remove()

    plt.tight_layout()
    fig.legend(h[:2], l[:2], frameon=False,
               loc='lower center', bbox_to_anchor=(0.5, -0.04), ncol=2)
    fig.savefig(out_path, bbox_inches='tight')
    return fig
