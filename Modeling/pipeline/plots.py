import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import config
from pipeline import signatures as sig
from viz_style import apply_style

apply_style()

MP = config.MODELING_PARAMS

UP, DN = '#d62728', '#1f77b4'
PALETTE = {'NBI': '#E41A1C', 'ZINBI': '#377EB8', 'Logistic': '#4DAF4A'}
BRANCH_COLOR = PALETTE
THEME_TINTS = ['#355475', '#AB5D10', '#3A7134', '#7C5471', '#9F3C3C', '#4F807C', '#A68D29', '#6D5141']
LEGEND_W = 0.15


# ── GSEA dotplots (gene_enrichment) ─────────────────────────────────────────
def _prep_plot_df(df_subset):
    if df_subset.empty:
        return df_subset
    d = df_subset.copy()
    d['Term_clean'] = d['Term'].str.split('__').str[1].fillna(d['Term'])
    d['tag_n'] = d['Tag %'].str.split('/').str[0].astype(float)
    d['FDR q-val'] = pd.to_numeric(d['FDR q-val'], errors='coerce').fillna(1.0)
    d['neg_log_q'] = -np.log10(d['FDR q-val'].clip(lower=1e-3))
    return d


def plot_gsea_dotplots(gsea_results, fdr_thr=None, top_n=None, fig_dir=None, sample_sizes=None):
    """Save per-phenotype bar/up/dn dotplots (gsea_bar/up/dn_{ph}.png)."""
    fdr_thr = MP['gsea_fdr_thr'] if fdr_thr is None else fdr_thr
    top_n = MP['gsea_top_n'] if top_n is None else top_n
    fig_dir = fig_dir or (config.GSEA_DIR / 'Figures' / 'gsea_dotplot')
    fig_dir.mkdir(parents=True, exist_ok=True)
    sample_sizes = sample_sizes or {}

    for ph, df in gsea_results.items():
        df = df.copy()
        df['NES'] = pd.to_numeric(df['NES'], errors='coerce')
        n_sig_pos = len(df[df['NES'] > 0])
        n_sig_neg = len(df[df['NES'] < 0])
        fname = ph.replace(' ', '_').replace('/', '_')
        n_samp = sample_sizes.get(ph, None)
        samp_str = f" (n={n_samp})" if n_samp else ""

        fig_bar, ax_bar = plt.subplots(figsize=(6, 2.5))
        bars = ax_bar.barh([0, 1], [n_sig_neg, n_sig_pos], color=[DN, UP], height=0.5, alpha=0.85)
        ax_bar.set_yticks([0, 1])
        ax_bar.set_yticklabels(['NES < 0 (Down)', 'NES > 0 (Up)'])
        ax_bar.set_xlabel('Total Significant Terms Count')
        ax_bar.set_title(f'GSEA Summary — {ph}{samp_str}\n(FDR < {fdr_thr})')
        max_count = max(n_sig_pos, n_sig_neg)
        for bar in bars:
            width = bar.get_width()
            if width > 0:
                ax_bar.text(width + (max_count * 0.02), bar.get_y() + bar.get_height() / 2,
                            f' {int(width)}', va='center', ha='left')
        ax_bar.spines['top'].set_visible(False)
        ax_bar.spines['right'].set_visible(False)
        if max_count > 0:
            ax_bar.set_xlim(0, max_count * 1.2)
        plt.tight_layout()
        plt.savefig(fig_dir / f'gsea_bar_{fname}.png', bbox_inches='tight')
        plt.close(fig_bar)

        pos_df = _prep_plot_df(df[df['NES'] > 0].nlargest(top_n, 'NES').copy())
        neg_df = _prep_plot_df(df[df['NES'] < 0].nsmallest(top_n, 'NES').copy())
        if not pos_df.empty:
            pos_df = pos_df.sort_values('NES', ascending=True)
        if not neg_df.empty:
            neg_df = neg_df.sort_values('NES', ascending=False)

        for sub, cmap, sign_lbl, color, xlim_fn, fn in [
            (pos_df, 'Reds', 'Upregulated Pathways (NES > 0)', UP,
             lambda d: (-0.2, d['NES'].max() + 0.5), f'gsea_up_{fname}.png'),
            (neg_df, 'Blues', 'Downregulated Pathways (NES < 0)', DN,
             lambda d: (d['NES'].min() - 0.5, 0.2), f'gsea_dn_{fname}.png'),
        ]:
            if sub.empty:
                continue
            fig, ax = plt.subplots(figsize=(10, max(6, len(sub) * 0.45)))
            ax.set_position([0.45, 0.18, 0.35, 0.75])
            scat = ax.scatter(sub['NES'], range(len(sub)), s=sub['tag_n'] * 15,
                              c=sub['neg_log_q'], cmap=cmap, vmin=0, alpha=0.85,
                              edgecolors='black', linewidths=1., zorder=3)
            ax.axvline(0, color='black', lw=1, ls='--', alpha=0.5)
            ax.set_yticks(range(len(sub)))
            ax.set_yticklabels(sub['Term_clean'])
            ax.set_xlim(*xlim_fn(sub))
            ax.set_title(f'{sign_lbl} — {ph}{samp_str}', color=color, pad=20)
            ax.set_xlabel('NES')
            ax.grid(alpha=0.3, axis='y')
            cax = fig.add_axes([0.83, 0.55, 0.02, 0.25])
            cb = fig.colorbar(scat, cax=cax, orientation='vertical')
            cb.ax.set_title('-log10\n(FDR)', pad=10, loc='left')
            tag_min, tag_max = int(sub['tag_n'].min()), int(sub['tag_n'].max())
            if tag_min == tag_max:
                size_vals = [tag_min]
            else:
                step = (tag_max - tag_min) / 3
                size_vals = sorted(set(int(tag_min + i * step) for i in range(4)))
            size_ex = [Line2D([0], [0], marker='o', color='w', markerfacecolor='#333333',
                              alpha=0.8, markersize=np.sqrt(n * 15), label=str(n)) for n in size_vals]
            ax.legend(handles=size_ex, title='Count', loc='upper left',
                      bbox_to_anchor=(1.05, 0.48), frameon=False, labelspacing=1.8)
            plt.savefig(fig_dir / fn, bbox_inches='tight')
            plt.close(fig)
        print(f'[{ph}] GSEA Plots generated successfully.')


# ── Heuristic signature figure (gsea_heuristic_signatures) ──────────────────
def plot_signature(ph, ctx, gsea_dir=None, fig_dir=None, themes=None, save=True):
    """Plot theme-shaded lollipop (left) and lead-gene specificity strip (right)."""
    fig_dir = fig_dir or (config.GSEA_DIR / 'Figures' / 'Heuristic_Signatures')
    fig_dir.mkdir(parents=True, exist_ok=True)
    rows, bands, lead_pool, yc = sig.theme_rows(ph, ctx, themes=themes, gsea_dir=gsea_dir)
    sym_to_idx, meanZ, phenos_u, samp_n = ctx.sym_to_idx, ctx.meanZ, ctx.phenos_u, ctx.samp_n

    rdf = pd.DataFrame(rows)
    H = max(8, yc * 0.35)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(20, H), gridspec_kw={'width_ratios': [1.55, 1]})

    for b in bands:
        axL.axhspan(b['y0'], b['y1'], color=THEME_TINTS[b['ti'] % len(THEME_TINTS)], alpha=0.13, zorder=0)
    cols = [UP if v > 0 else DN for v in rdf['nes']]
    axL.hlines(rdf['y'], 0, rdf['nes'], color=cols, lw=1.2, alpha=0.55, zorder=2)
    axL.scatter(rdf['nes'], rdf['y'], s=np.clip(rdf['tag'] * 7, 25, 300),
                c=cols, edgecolors='black', linewidths=0.7, zorder=3)
    axL.axvline(0, color='black', lw=1, ls='--', alpha=0.5)
    axL.set_yticks(rdf['y'])
    axL.set_yticklabels(rdf['term'])
    axL.tick_params(axis='y', length=0)
    axL.set_ylim(-1, yc - 0.4)
    axL.invert_yaxis()
    for tk, th in zip(axL.get_yticklabels(), rdf['theme']):
        tk.set_color(THEME_TINTS[int(th) % len(THEME_TINTS)])
    axL.set_xlabel('NES')
    axL.set_title(f'Signature pathways by theme — {ph} (n={samp_n[ph]})', pad=10)
    axL.grid(axis='x', alpha=0.3)
    axL.margins(x=0.12)

    seen = []
    for gg in lead_pool:
        if gg in sym_to_idx and gg not in seen:
            seen.append(gg)
    tgt0 = np.array([meanZ[ph][sym_to_idx[gg]] for gg in seen])
    genes = [seen[i] for i in np.argsort(-np.abs(tgt0))[:16]]
    allmat = np.array([[meanZ[p][sym_to_idx[gg]] for p in phenos_u] for gg in genes])
    tgt = np.array([meanZ[ph][sym_to_idx[gg]] for gg in genes])
    order = np.argsort(tgt)
    genes = [genes[i] for i in order]
    allmat = allmat[order]
    tgt = tgt[order]
    yy = np.arange(len(genes))
    for i in range(len(genes)):
        axR.scatter(allmat[i], [yy[i]] * allmat.shape[1], s=12, c='#bbbbbb', alpha=0.6, zorder=1)
    axR.scatter(tgt, yy, s=80, c=[UP if v > 0 else DN for v in tgt],
                edgecolors='black', linewidths=0.8, zorder=3)
    axR.axvline(0, color='black', lw=1, ls='--', alpha=0.5)
    axR.set_yticks(yy)
    axR.set_yticklabels(genes)
    axR.set_ylim(-1, len(genes))
    axR.set_xlabel('mean Z-score')
    axR.set_title('Lead-gene specificity')
    axR.grid(axis='x', alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1 - LEGEND_W, 1])
    handles = [Patch(facecolor=THEME_TINTS[b['ti'] % len(THEME_TINTS)], alpha=0.5, label=b['label'])
               for b in bands]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1 - LEGEND_W + 0.005, 0.5),
               frameon=False, handlelength=1.2, title='Theme  (n = FDR<0.05 Pathways)',
               alignment='left')
    if save:
        plt.savefig(fig_dir / f"signature_{ph.replace(' ', '_').replace('/', '_')}.png",
                    bbox_inches='tight')
    plt.show()


# ── gene_selection plots ────────────────────────────────────────────────────
def plot_zscore_outlier_hist(Z_dis, dis_pheno, thresh=None, fig_dir=None, save=True):
    thresh = MP['z_flag'] if thresh is None else thresh
    fig_dir = fig_dir or config.CV_FIG_DIR
    fig_dir.mkdir(parents=True, exist_ok=True)
    pheno_list = np.unique(dis_pheno)
    counts = {ph: (np.abs(Z_dis[dis_pheno == ph]) > thresh).sum(axis=1) for ph in pheno_list}
    order = sorted(pheno_list, key=lambda ph: -np.median(counts[ph]))
    n_pheno = len(order); n_cols = 4; n_rows = (n_pheno + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows))
    axes = np.array(axes).flatten()
    for ax, ph in zip(axes, order):
        c = counts[ph]
        ax.hist(c, bins=25, color='gray', edgecolor='white', linewidth=0.5)
        ax.axvline(np.median(c), color='tomato', linestyle='--', linewidth=1.0,
                   label=f'med={np.median(c):.0f}')
        ax.set_title(f'{ph} (n={len(c)})')
        ax.set_xlabel(f'# genes |Z|>{thresh}'); ax.set_ylabel('samples')
        ax.legend(frameon=False)
    for ax in axes[n_pheno:]:
        ax.set_visible(False)
    fig.suptitle(f'Per-phenotype distribution of |Z| > {thresh} gene count', y=1.01)
    plt.tight_layout()
    if save:
        plt.savefig(fig_dir / 'zscore_outlier_gene_dist.png', bbox_inches='tight')
    plt.show()


def _roc_shadow(ax, roc_d, color, label, auc, lw=1.5, alpha=0.18):
    if not roc_d or not roc_d.get('fprs'):
        ax.text(0.5, 0.5, 'N/A', ha='center', va='center', transform=ax.transAxes, color='grey')
        return
    tprs = [np.interp(np.linspace(0, 1, 101), f, t) for f, t in zip(roc_d['fprs'], roc_d['tprs'])]
    mean_t = np.mean(tprs, axis=0); std_t = np.std(tprs, axis=0)
    auc_str = f'{auc:.2f}' if not np.isnan(auc) else '—'
    base = np.linspace(0, 1, 101)
    ax.plot(base, mean_t, color=color, lw=lw, label=f'{label} AUC={auc_str}')
    ax.fill_between(base, np.clip(mean_t - std_t, 0, 1), np.clip(mean_t + std_t, 0, 1),
                    alpha=alpha, color=color)


def _pr_shadow(ax, pr_d, color, label, ap, prevalence=None, lw=1.5, alpha=0.18):
    if not pr_d or not pr_d.get('recs'):
        ax.text(0.5, 0.5, 'N/A', ha='center', va='center', transform=ax.transAxes, color='grey')
        return
    base = np.linspace(0, 1, 101)
    precs = []
    for rec, prec in zip(pr_d['recs'], pr_d['precs']):
        order = np.argsort(rec)
        precs.append(np.interp(base, rec[order], prec[order]))
    mean_p = np.mean(precs, axis=0); std_p = np.std(precs, axis=0)
    ap_str = f'{ap:.2f}' if not np.isnan(ap) else '—'
    ax.plot(base, mean_p, color=color, lw=lw, label=f'{label} AP={ap_str}')
    ax.fill_between(base, np.clip(mean_p - std_p, 0, 1), np.clip(mean_p + std_p, 0, 1),
                    alpha=alpha, color=color)
    if prevalence is not None:
        ax.axhline(prevalence, color='k', lw=0.7, ls='--', alpha=0.35)


def _curve_grid(curves, value_of, kind, color, label, suptitle, fname, fig_dir, save):
    """Familiar per-phenotype shadow-curve grid (ROC or PR), mirroring plot_roc_curves."""
    phenos = sorted(curves, key=lambda p: -(curves[p]['n'] if curves[p] else 0))
    ncols = 5; nrows = (len(phenos) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 3.5))
    axf = axes.flatten()
    for ax, ph in zip(axf, phenos):
        cv = curves[ph]
        if cv is None:
            ax.text(0.5, 0.5, '(excluded)', ha='center', va='center',
                    transform=ax.transAxes, color='grey')
            ax.set_title(f'{ph}'); ax.axis('off'); continue
        val = value_of.get(ph, np.nan)
        if kind == 'roc':
            _roc_shadow(ax, cv['roc'], color, label, val)
            ax.plot([0, 1], [0, 1], 'k--', lw=0.7, alpha=0.35)
            ax.set_xlabel('FPR'); ax.set_ylabel('TPR')
        else:
            _pr_shadow(ax, cv['pr'], color, label, val, prevalence=cv['prevalence'])
            ax.set_xlabel('Recall'); ax.set_ylabel('Precision')
        ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
        ax.set_title(f'{ph}  (n={cv["n"]})'); ax.tick_params(labelsize=6)
        ax.legend(frameon=False, loc='lower right' if kind == 'roc' else 'upper right',
                  fontsize=7)
    for ax in axf[len(phenos):]:
        ax.axis('off')
    fig.suptitle(suptitle, y=1.005)
    plt.tight_layout()
    if save:
        plt.savefig(fig_dir / fname, bbox_inches='tight', dpi=150)
    plt.show()


NEG_MODE_LABEL = {'real_hc': 'real HC cohort', 'null_fixed': 'healthy null N(0,1)',
                  'null_matched': 'healthy null N(0,1), count-matched'}


def plot_validation_curves(method, curves, binary_df, multi_df, controls,
                           neg_modes=None, fig_dir=None, save=True):
    """Nested-CV validation in the familiar shadow-curve idiom: binary ROC + PR grids for
    each negative mode (real_hc / null_fixed / null_matched), multiclass disease-vs-disease
    ROC + PR, and the real-HC-vs-null calibration panel (hugs diagonal/baseline if
    calibrated)."""
    fig_dir = fig_dir or config.CV_FIG_DIR
    fig_dir.mkdir(parents=True, exist_ok=True)
    cv = curves[method]
    neg_modes = neg_modes or list(cv['binary'].keys())
    mdf = multi_df[multi_df['method'] == method].set_index('phenotype')

    for nm in neg_modes:
        bdf = binary_df[(binary_df['method'] == method) & (binary_df['neg_mode'] == nm)] \
            .set_index('phenotype')
        lab = NEG_MODE_LABEL.get(nm, nm)
        _curve_grid(cv['binary'][nm], bdf['auc'].to_dict(), 'roc', '#377eb8', 'LogReg',
                    f'Disease vs {lab} — ROC (nested) — {method}  ({MP["n_splits"]}-fold, shadow=±1 SD)',
                    f'val_roc_binary_{method}_{nm}.png', fig_dir, save)
        _curve_grid(cv['binary'][nm], bdf['auprc'].to_dict(), 'pr', '#377eb8', 'LogReg',
                    f'Disease vs {lab} — PR (nested) — {method}  (dashed=prevalence baseline)',
                    f'val_pr_binary_{method}_{nm}.png', fig_dir, save)
    _curve_grid(cv['multi'], mdf['auc_mc'].to_dict(), 'roc', '#4daf4a', 'OVR',
                f'Disease vs Disease ROC (nested) — {method}',
                f'val_roc_multiclass_{method}.png', fig_dir, save)
    _curve_grid(cv['multi'], mdf['auprc_mc'].to_dict(), 'pr', '#4daf4a', 'OVR',
                f'Disease vs Disease PR (nested) — {method}  (dashed=prevalence baseline)',
                f'val_pr_multiclass_{method}.png', fig_dir, save)

    ctrl_auc, _, hc = controls[method]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    _roc_shadow(ax1, hc['roc'], '#999999', 'HC vs null', ctrl_auc)
    ax1.plot([0, 1], [0, 1], 'k--', lw=0.7, alpha=0.35)
    ax1.set_xlim(-0.05, 1.05); ax1.set_ylim(-0.05, 1.05)
    ax1.set_xlabel('FPR'); ax1.set_ylabel('TPR'); ax1.set_title('ROC')
    ax1.legend(frameon=False, loc='lower right')
    _pr_shadow(ax2, hc['pr'], '#999999', 'HC vs null', np.nan, prevalence=hc['prevalence'])
    ax2.set_xlim(-0.05, 1.05); ax2.set_ylim(-0.05, 1.05)
    ax2.set_xlabel('Recall'); ax2.set_ylabel('Precision'); ax2.set_title('PR')
    ax2.legend(frameon=False, loc='upper right')
    fig.suptitle(f'Calibration control: real HC vs healthy null — {method}  '
                 f'(~0.5 / baseline = calibrated)', y=1.02)
    plt.tight_layout()
    if save:
        plt.savefig(fig_dir / f'val_negcontrol_hc_{method}.png', bbox_inches='tight', dpi=150)
    plt.show()


def plot_selection_overview(method_name, all_results, Z_dis, dis_pheno, dis_names,
                            gene_names, fig_dir=None, save=True):
    """Save clustermap and UMAP for the selected gene set (heatmap_{m}.png, umap_{m}.png)."""
    import seaborn as sns
    from umap import UMAP
    fig_dir = fig_dir or config.CV_FIG_DIR
    fig_dir.mkdir(parents=True, exist_ok=True)
    res = all_results[method_name]
    g2i = {g: i for i, g in enumerate(gene_names)}
    idx = [g2i[g] for g in res['genes'] if g in g2i]
    X = Z_dis[:, idx]
    unique_phenos = sorted(np.unique(dis_pheno))
    p2c = dict(zip(unique_phenos, sns.color_palette('tab20', len(unique_phenos))))
    rc = [p2c[p] for p in dis_pheno]

    pivot = pd.DataFrame(X, index=dis_names, columns=res['genes'])
    g = sns.clustermap(
        pivot, method='ward', metric='euclidean',
        row_colors=rc, col_cluster=True, row_cluster=True,
        cmap='RdBu_r', vmin=-6, vmax=6, center=0,
        xticklabels=False, yticklabels=False,
        figsize=(20, 13), linewidths=0, rasterized=True,
        cbar_kws=dict(label='Z-score', shrink=0.35, ticks=[-6, -3, 0, 3, 6]))
    patches = [Patch(color=p2c[p], label=p) for p in unique_phenos]
    g.ax_heatmap.legend(handles=patches, loc='upper left',
                        bbox_to_anchor=(1.15, 0.7), frameon=False, title='Phenotype')
    g.fig.suptitle(f'Heatmap — {method_name}  ({len(res["genes"])} genes)', y=1.01)
    if save:
        plt.savefig(fig_dir / f'heatmap_{method_name}.png', bbox_inches='tight', dpi=300)
    plt.show()

    emb = UMAP(n_components=2, metric='euclidean', n_neighbors=15,
               min_dist=0.2, random_state=42).fit_transform(X)
    fig, ax = plt.subplots(figsize=(9, 7))
    for ph in unique_phenos:
        m = dis_pheno == ph
        ax.scatter(emb[m, 0], emb[m, 1], c=[p2c[ph]], s=12, alpha=0.75, label=ph, edgecolors='none')
    ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc='upper left', title='Phenotype')
    ax.set(title=f'UMAP — {method_name}  ({len(res["genes"])} genes)', xlabel='UMAP 1', ylabel='UMAP 2')
    plt.tight_layout()
    if save:
        plt.savefig(fig_dir / f'umap_{method_name}.png', bbox_inches='tight', dpi=300)
    plt.show()


# ── disease_scoring per-sample Manhattan (disease_scoring) ──────────────────
_SCORE_BRANCH_COLOR = {'count': PALETTE['NBI'], 'logistic': PALETTE['Logistic'], 'rare': '#984EA3'}


def plot_sample(df, sample_id, phenotype='', top_n=20, z_flag=None):
    z_flag = MP['z_flag'] if z_flag is None else z_flag
    df = df.dropna(subset=['score']).copy()
    df['abs_score'] = df['score'].abs()
    df_sorted = df.sort_values('score', ascending=False).reset_index(drop=True)
    flagged = df[df['abs_score'] >= z_flag].sort_values('abs_score', ascending=False)
    fig, axes = plt.subplots(1, 2, figsize=(18, 5), gridspec_kw={'width_ratios': [3, 1]})

    ax = axes[0]
    for branch, sub in df_sorted.groupby('branch'):
        ax.scatter(sub.index, sub['score'], s=0.2, alpha=0.25,
                   color=_SCORE_BRANCH_COLOR.get(branch, 'grey'), label=branch, rasterized=True)
    flag_sub = df_sorted[df_sorted['abs_score'] >= z_flag]
    ax.scatter(flag_sub.index, flag_sub['score'], s=5, color='black', zorder=5, alpha=0.8)
    for _, row in flag_sub.head(5).iterrows():
        ax.annotate(row['gene'], (row.name, row['score']), xytext=(5, 3),
                    textcoords='offset points', fontsize=7, alpha=0.85)
    ax.axhline(z_flag, color='red', lw=1, ls='--', alpha=0.6, label=f'|z|={z_flag}')
    ax.axhline(-z_flag, color='red', lw=1, ls='--', alpha=0.6)
    ax.axhline(0, color='grey', lw=0.8, ls='-', alpha=0.4)
    ax.set_xlabel('Genes (sorted by score)')
    ax.set_ylabel('Anomaly Score (z / rare_score)')
    ax.set_title(f'{sample_id}\n{phenotype}', fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')

    ax2 = axes[1]
    top = flagged.head(top_n)
    if len(top) == 0:
        ax2.text(0.5, 0.5, 'No flagged genes', ha='center', va='center',
                 transform=ax2.transAxes, color='grey')
        ax2.axis('off')
    else:
        colors = [_SCORE_BRANCH_COLOR.get(b, 'grey') for b in top['branch']]
        ax2.barh(range(len(top)), top['abs_score'].values, color=colors, alpha=0.8, edgecolor='white')
        ax2.set_yticks(range(len(top)))
        ax2.set_yticklabels([f"{g}  ({c:.0f}ct)" for g, c in
                             zip(top['gene'], top['raw_count'].fillna(0))], fontsize=7)
        ax2.invert_yaxis()
        ax2.axvline(z_flag, color='red', lw=1, ls='--', alpha=0.6)
        ax2.set_xlabel('|Score|')
        ax2.set_title(f'Top {len(top)} Flagged Genes', fontweight='bold')
        legend_els = [Patch(facecolor=c, label=b) for b, c in _SCORE_BRANCH_COLOR.items()]
        ax2.legend(handles=legend_els, fontsize=7, loc='lower right')
    plt.tight_layout()
    return fig


def _strip_panel(ax, detail, value_col, ensg_to_y, n_genes, rng, xlabel, xscale=None,
                 vlines=()):
    """One column of a per-gene strip plot: every disease sample as a point, jittered on
    y, flagged samples in red on top of grey others -- so multiple flagged samples per
    gene all stay visible instead of collapsing into a single summary bar."""
    for is_flag, color, size, alpha, z in [(False, 'grey', 14, 0.5, 2), (True, 'tomato', 45, 0.9, 3)]:
        sub = detail[detail['flagged'] == is_flag]
        y = sub['ensg'].map(ensg_to_y) + rng.uniform(-0.18, 0.18, len(sub))
        ax.scatter(sub[value_col], y, s=size, color=color, alpha=alpha, edgecolors='none',
                  zorder=z, label='flagged sample' if is_flag else 'other disease samples')
    if xscale:
        ax.set_xscale(xscale)
    for v in vlines:
        ax.axvline(v, color='black', lw=1, ls='--', alpha=0.5)
    ax.set_ylim(n_genes - 0.5, -0.5)
    ax.set_xlabel(xlabel)
    ax.legend(fontsize=7, loc='lower right', frameon=False)


def plot_rescued_genes(merged, rescued, phenotype, dd=None, z_flag=None, top_n=20,
                       fig_dir=None, save=True):
    """DESeq2-excluded genes vs normative signal: overview scatter + per-gene raw-count
    and Z-score strip plots (one point per disease sample, flagged ones highlighted)."""
    from pipeline.benchmark import gene_sample_detail
    z_flag = MP['z_flag'] if z_flag is None else z_flag
    fig_dir = fig_dir or config.BENCHMARK_DIR / 'Figures'
    fig_dir.mkdir(parents=True, exist_ok=True)
    excl = merged[merged['excluded']]
    kept = merged[~merged['excluded']]

    fig, axes = plt.subplots(1, 3, figsize=(19, 5))
    ax = axes[0]
    ax.scatter(kept['baseMean'] + 1, kept['norm_max_abs_z'], s=8, alpha=0.25,
               color='grey', label=f'DESeq2-tested (n={len(kept)})')
    ax.scatter(excl['baseMean'] + 1, excl['norm_max_abs_z'], s=14, alpha=0.7,
               color='tomato', label=f'DESeq2-excluded (n={len(excl)})')
    ax.axhline(z_flag, color='black', lw=1, ls='--', alpha=0.5)
    ax.set_xscale('log')
    ax.set_xlabel('DESeq2 baseMean + 1')
    ax.set_ylabel('Normative max |Z|\n(flagged calls only; 0 = never flagged)')
    ax.legend(fontsize=8, frameon=False)

    top = rescued.head(top_n)
    for ax in axes[1:]:
        if len(top) == 0:
            ax.text(0.5, 0.5, 'No rescued genes above threshold', ha='center', va='center',
                    transform=ax.transAxes, color='grey')
            ax.axis('off')
        elif dd is None:
            ax.text(0.5, 0.5, 'pass dd= for per-sample strip plot', ha='center', va='center',
                    transform=ax.transAxes, color='grey')
            ax.axis('off')

    if len(top) > 0 and dd is not None:
        detail = gene_sample_detail(dd, phenotype, top[['ensg', 'branch']])
        rng = np.random.default_rng(0)
        labels = [f"{s} ({b}, n_flag={n:.0f})" for s, b, n in
                  zip(top['symbol'], top['branch'].fillna('?'), top['norm_n_flagged'])]
        ensg_to_y = {g: i for i, g in enumerate(top['ensg'])}

        _strip_panel(axes[1], detail.assign(raw_count=detail['raw_count'] + 1), 'raw_count',
                    ensg_to_y, len(top), rng, 'Raw count + 1 (one point per disease sample)',
                    xscale='log')
        axes[1].set_yticks(range(len(top)))
        axes[1].set_yticklabels(labels, fontsize=7)
        axes[1].invert_yaxis()
        axes[1].set_title(f'Top rescued genes (n={len(rescued)} total)')

        _strip_panel(axes[2], detail, 'z', ensg_to_y, len(top), rng,
                    'Z-score (one point per disease sample)', vlines=(z_flag, -z_flag))
        axes[2].set_yticks(range(len(top)))
        axes[2].set_yticklabels([])
        axes[2].invert_yaxis()
        axes[2].set_title('Z-score distribution within disease group')

    plt.tight_layout()
    if save:
        fname = phenotype.replace(' ', '_').replace('/', '_')
        plt.savefig(fig_dir / f'rescued_genes_{fname}.png', bbox_inches='tight')
    plt.show()
    return fig


DB_METHOD_STYLE = {'deseq2': ('#DD8452', 'DESeq2'),
                   'deseq2_cov': ('#C44E52', 'DESeq2 + covariates'),
                   'no_filter': ('#4C72B0', 'Normative (no_filter)'),
                   'with_rare': ('#55A868', 'Normative (with_rare)')}


def plot_db_hit_rates(rates, summary, fig_dir=None, save=True):
    """Symmetric DB-support comparison, counts first. Left: per-phenotype DB-supported term
    counts (n_db) for each method. Right: pooled DB-supported (filled) inside all-significant
    (outline) with the pooled DB-hit rate annotated, so absolute coverage and precision are
    both visible. Only phenotypes with an Open Targets reference are shown."""
    fig_dir = fig_dir or (config.BENCHMARK_DIR / 'Figures')
    fig_dir.mkdir(parents=True, exist_ok=True)
    sub = rates[rates['has_ot_ref']].copy()
    methods = [m for m in DB_METHOD_STYLE if m in set(sub['method'])]
    order = sub.groupby('phenotype')['n_db'].sum().sort_values().index.tolist()
    piv = sub.pivot(index='phenotype', columns='method', values='n_db').reindex(order)
    y = np.arange(len(order))
    h = 0.8 / len(methods)
    fig, axes = plt.subplots(1, 2, figsize=(13, 8), gridspec_kw={'width_ratios': [2.4, 1]})
    ax = axes[0]
    for i, m in enumerate(methods):
        c, lab = DB_METHOD_STYLE[m]
        ax.barh(y + (i - (len(methods) - 1) / 2) * h, piv[m].values, height=h, color=c, label=lab)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=8)
    ax.set_xlabel('DB-supported significant terms (count)')
    ax.set_title('Per-phenotype DB-hit counts')
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    ax2 = axes[1]
    sm = summary.set_index('method').reindex(methods)
    xb = np.arange(len(methods))
    cols = [DB_METHOD_STYLE[m][0] for m in methods]
    ax2.bar(xb, sm['total_sig'].values, color='none', edgecolor='grey', lw=1.0)
    ax2.bar(xb, sm['total_db'].values, color=cols)
    for k, m in enumerate(methods):
        ax2.text(k, sm.loc[m, 'total_db'],
                 f"{int(sm.loc[m, 'total_db'])}\n({sm.loc[m, 'pooled_db_hit_rate']:.2f})",
                 ha='center', va='bottom', fontsize=8)
    ax2.set_xticks(xb)
    ax2.set_xticklabels([DB_METHOD_STYLE[m][1] for m in methods], rotation=30, ha='right', fontsize=8)
    ax2.set_ylabel('pooled term count')
    ax2.set_title('DB-supported (fill) vs all significant (outline)\nlabel: n_db (pooled DB-hit rate)')
    fig.tight_layout()
    if save:
        plt.savefig(fig_dir / 'db_hit_rates.png', bbox_inches='tight', dpi=200)
    plt.show()
    return fig
