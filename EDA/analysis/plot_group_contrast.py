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


def private_calls(D, max_detect=0.2, min_reads=3):
    """A sample's private expression: gene detected (>=1 read) in fewer than max_detect of the
    reference controls, and the sample has >= min_reads with CPM above every reference control.
    Patients use all controls; each control uses the other controls (leave-one-out)."""
    C, P = D["counts"], D["cpm"]
    hc, pat = D["is_hc"], D["is_case"]
    Ch, Ph, n = C[hc], P[hc], int(hc.sum())
    nz = (Ch > 0).sum(axis=0)
    top2 = -np.sort(-Ph, axis=0)[:2]
    Mp = ((nz / n) < max_detect) & (P[pat] > top2[0]) & (C[pat] >= min_reads)
    loo_max = np.where(Ph >= top2[0], top2[1], top2[0])
    Mh = (((nz - (Ch > 0)) / (n - 1)) < max_detect) & (Ph > loo_max) & (Ch >= min_reads)
    return Mp, Mh


def _loo_median(A):
    n = A.shape[0]
    S = np.sort(A, axis=0)
    r = np.argsort(np.argsort(A, axis=0, kind="stable"), axis=0, kind="stable")
    pick = lambda k: np.where(k < r, S[k][None, :], S[k + 1][None, :])
    m = n - 1
    return pick(m // 2) if m % 2 else (pick(m // 2 - 1) + pick(m // 2)) / 2


def quiet_calls(D, spread=1.0, lift=2.0, min_reads=3):
    """Conservative private expression, restricted to genes that are quiet in the reference
    controls: no control lies more than `spread` log2(CPM + 1) above the control median (stable,
    no control outlier). A sample is called when it sits >= `lift` above that median with
    >= min_reads. Patients use all controls; each control uses the other controls."""
    Y, C = np.log2(D["cpm"] + 1), D["counts"]
    hc, pat = D["is_hc"], D["is_case"]
    Yh = Y[hc]
    med, top2 = np.median(Yh, axis=0), -np.sort(-Yh, axis=0)[:2]
    Mp = ((top2[0] - med) <= spread) & (Y[pat] - med >= lift) & (C[pat] >= min_reads)
    med_l = _loo_median(Yh)
    max_l = np.where(Yh >= top2[0], top2[1], top2[0])
    Mh = ((max_l - med_l) <= spread) & (Yh - med_l >= lift) & (C[hc] >= min_reads)
    return Mp, Mh


def _call_row(D, Mp, Mh, de_any, **rule):
    g = Mp.any(axis=0)
    return dict(scope=COHORT_LABEL[D["cohort"]], **rule, private_genes=int(g.sum()),
                single_patient=int((Mp.sum(0) == 1).sum()), also_DE_any_design=int((g & de_any).sum()),
                per_patient_med=np.median(Mp.sum(1)), per_control_med=np.median(Mh.sum(1)),
                patients_with_any=int(Mp.any(1).sum()), controls_with_any=int(Mh.any(1).sum()),
                n_patients=len(Mp), n_controls=len(Mh))


def private_sensitivity(Ds, detects=(0.05, 0.1, 0.2, 0.3), reads=(3, 5, 10)):
    rows = []
    for D in Ds:
        de_any = np.isin(D["genes"], list(set().union(*[de_genes(D, d) for d in DESIGNS])))
        for md in detects:
            for mr in reads:
                rows.append(_call_row(D, *private_calls(D, md, mr), de_any, rule="detect",
                                      max_hc_detect=md, min_reads=mr))
    return pd.DataFrame(rows)


def quiet_sensitivity(Ds, spreads=(0.5, 1.0, 1.5), lifts=(1.0, 2.0, 3.0), min_reads=3):
    rows = []
    for D in Ds:
        de_any = np.isin(D["genes"], list(set().union(*[de_genes(D, d) for d in DESIGNS])))
        for sp in spreads:
            for lf in lifts:
                if lf <= sp:
                    continue
                rows.append(_call_row(D, *quiet_calls(D, sp, lf, min_reads), de_any, rule="quiet",
                                      spread=sp, lift=lf, min_reads=min_reads))
    return pd.DataFrame(rows)


def plot_deg_vs_private(D, Mp, rule, design="no_covariate", save_path=None):
    Y = np.log2(D["cpm"] + 1)
    Y = Y - np.median(Y[D["is_hc"]], axis=0)
    hc_i, pat_i = np.where(D["is_hc"])[0], np.where(D["is_case"])[0]
    rows = np.concatenate([hc_i, pat_i])
    G = de_genes(D, design)
    deg = np.array([D["pos"][g] for g in G], dtype=int)
    deg = deg[np.argsort(Y[pat_i][:, deg].mean(axis=0))]
    n_pat = Mp.sum(axis=0)
    pv = np.where((n_pat > 0) & ~np.isin(D["genes"], G))[0]
    first = Mp[:, pv].argmax(axis=0)
    pv = pv[np.lexsort((first, -n_pat[pv]))]
    blocks = [(deg, f"{len(deg)} group DE genes (DESeq2, RUVg {DESIGN_SHORT[design]})"),
              (pv, f"{len(pv)} patient-level signals undetected by the group test "
                   f"({int((n_pat[pv] == 1).sum())} in a single patient)")]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True, constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[1, 2]))
    for ax, (cols, lab) in zip(axes, blocks):
        if len(cols):
            im = ax.imshow(Y[np.ix_(rows, cols)], cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto",
                           interpolation="nearest")
        ax.axhline(len(hc_i) - 0.5, color="k", lw=1.2)
        ax.set_xlabel(lab)
        ax.set_xticks([])
    axes[0].set_yticks([len(hc_i) / 2, len(hc_i) + len(pat_i) / 2],
                       [f"controls ({len(hc_i)})", f"patients ({len(pat_i)})"])
    fig.colorbar(im, ax=axes, fraction=0.02, pad=0.01).set_label("log2(CPM + 1) - control median")
    fig.suptitle(f"{COHORT_LABEL[D['cohort']]}: group DE genes vs patient-private expression ({rule})",
                 fontsize=11)
    _save(fig, save_path)
    return fig
