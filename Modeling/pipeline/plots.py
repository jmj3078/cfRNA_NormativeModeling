import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import config
from pipeline import signatures as sig
from viz_style import apply_style

apply_style()

MP = config.MODELING_PARAMS

# 공용 시각화 상수
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
    """phenotype별 bar/up/dn dotplot 저장 (gsea_bar/up/dn_{ph}.png)."""
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
    """좌: theme 음영 lollipop / 우: lead-gene 특이성 strip."""
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


def plot_roc_curves(method_name, all_results, dis_pheno, cv=None, fig_dir=None, save=True):
    cv = MP['n_splits'] if cv is None else cv
    fig_dir = fig_dir or config.CV_FIG_DIR
    fig_dir.mkdir(parents=True, exist_ok=True)
    res = all_results[method_name]; pp = res['per_pheno']
    phenos = sorted(np.unique(dis_pheno))
    ncols = 5; nrows = (len(phenos) + ncols - 1) // ncols

    # (1) Binary
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 3.5))
    axf = axes.flatten()
    for ax, ph in zip(axf, phenos):
        rb = res['roc_bin'].get(ph); n = (dis_pheno == ph).sum()
        auc_lr = pp.loc[ph, 'auc_logreg'] if ph in pp.index else np.nan
        auc_rf = pp.loc[ph, 'auc_rf'] if ph in pp.index else np.nan
        if rb is None:
            ax.text(0.5, 0.5, f'n={n} < {cv}\n(excluded)', ha='center', va='center',
                    transform=ax.transAxes, color='grey')
            ax.set_title(f'{ph}  (n={n})'); ax.axis('off'); continue
        _roc_shadow(ax, rb.get('lr'), '#377eb8', 'LogReg', auc_lr)
        _roc_shadow(ax, rb.get('rf'), '#e41a1c', 'RF', auc_rf)
        ax.plot([0, 1], [0, 1], 'k--', lw=0.7, alpha=0.35)
        ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel('FPR'); ax.set_ylabel('TPR'); ax.set_title(f'{ph}  (n={n})')
        ax.legend(frameon=False, loc='lower right')
    for ax in axf[len(phenos):]:
        ax.axis('off')
    fig.suptitle(f'HC vs Disease ROC — {method_name}  ({cv}-fold, shadow = ±1 SD)', y=1.01)
    plt.tight_layout()
    if save:
        plt.savefig(fig_dir / f'roc_binary_{method_name}.png', bbox_inches='tight', dpi=150)
    plt.show()

    # (2) Multiclass OVR
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 3.5))
    axf = axes.flatten()
    for ax, ph in zip(axf, phenos):
        rm = res['roc_mc'].get(ph, {}); n = (dis_pheno == ph).sum()
        auc = pp.loc[ph, 'auc_multiclass'] if ph in pp.index else np.nan
        _roc_shadow(ax, rm, '#4daf4a', 'OVR LogReg', auc)
        ax.plot([0, 1], [0, 1], 'k--', lw=0.7, alpha=0.35)
        ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
        ax.set_title(f'{ph}  (n={n})'); ax.set_xlabel('FPR'); ax.set_ylabel('TPR')
        ax.tick_params(labelsize=5); ax.legend(frameon=False, loc='lower right')
    for ax in axf[len(phenos):]:
        ax.axis('off')
    fig.suptitle(f'Multiclass OVR ROC — {method_name}  ({cv}-fold, shadow = ±1 SD)', y=1.01)
    plt.tight_layout()
    if save:
        plt.savefig(fig_dir / f'roc_multiclass_{method_name}.png', bbox_inches='tight', dpi=150)
    plt.show()


def plot_selection_overview(method_name, all_results, Z_dis, dis_pheno, dis_names,
                            gene_names, fig_dir=None, save=True):
    """선택 gene set의 clustermap + UMAP (heatmap_{m}.png, umap_{m}.png)."""
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
