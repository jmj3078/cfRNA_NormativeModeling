"""
Out-of-Distribution (OOD) sample filtering based on Mahalanobis distance.

Filters disease samples whose covariate (X) profile lies outside the HC
distribution, using HC-empirical percentile thresholds.

Usage
-----
    from sample_filter import MahalanobisFilter

    filt = MahalanobisFilter(percentile=95)
    filt.fit(X_hc)

    keep         = filt.mask(X_dis)          # boolean (n_dis,)
    summary_df   = filt.summary(X_dis, dis_pheno)
    figs         = filt.plot(X_dis, dis_pheno, keep, FIG_DIR)
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import mahalanobis
from numpy.linalg import inv


class MahalanobisFilter:
    def __init__(self, percentile: float = 95, reg: float = 1e-8):
        """
        Parameters
        ----------
        percentile : HC 분포에서 threshold로 사용할 백분위수 (default 95)
        reg        : 공분산 행렬 정규화 항
        """
        self.percentile = percentile
        self.reg        = reg
        self._fitted    = False

    # ── Fit ───────────────────────────────────────────────────────

    def fit(self, X_hc: np.ndarray) -> "MahalanobisFilter":
        """HC 샘플로 Mahalanobis 파라미터 추정."""
        X_hc = np.asarray(X_hc, dtype=float)
        p    = X_hc.shape[1]
        self.mu_      = X_hc.mean(axis=0)
        self.cov_inv_ = inv(np.cov(X_hc.T) + np.eye(p) * self.reg)
        self.hc_dist_ = self._distances(X_hc)
        self.thr_p50_ = float(np.percentile(self.hc_dist_, 50))
        self.thr_p95_ = float(np.percentile(self.hc_dist_, 95))
        self.thr_p99_ = float(np.percentile(self.hc_dist_, 99))
        self.threshold_ = float(np.percentile(self.hc_dist_, self.percentile))
        self._fitted = True
        return self

    # ── Core ──────────────────────────────────────────────────────

    def _distances(self, X: np.ndarray) -> np.ndarray:
        return np.array([mahalanobis(x, self.mu_, self.cov_inv_) for x in X])

    def distances(self, X: np.ndarray) -> np.ndarray:
        """Mahalanobis distance of each sample from HC distribution."""
        self._check_fit()
        return self._distances(np.asarray(X, dtype=float))

    def mask(self, X: np.ndarray) -> np.ndarray:
        """Boolean keep-mask: True = within threshold (inlier)."""
        return self.distances(X) <= self.threshold_

    # ── Summary ───────────────────────────────────────────────────

    def summary(self, X_dis: np.ndarray, dis_pheno: np.ndarray) -> pd.DataFrame:
        """Per-phenotype filtering summary DataFrame."""
        self._check_fit()
        d    = self.distances(X_dis)
        keep = d <= self.threshold_
        rows = []
        for ph in np.unique(dis_pheno):
            m = dis_pheno == ph
            rows.append({
                'phenotype':    ph,
                'n_before':     int(m.sum()),
                'n_after':      int((m & keep).sum()),
                'n_removed':    int((m & ~keep).sum()),
                'pct_removed':  round((m & ~keep).sum() / m.sum() * 100, 1),
                'mahal_mean':   round(d[m].mean(), 2),
                'mahal_max':    round(d[m].max(),  2),
                'pct_>p95_HC':  round((d[m] > self.thr_p95_).mean() * 100, 1),
                'pct_>p99_HC':  round((d[m] > self.thr_p99_).mean() * 100, 1),
            })
        return (pd.DataFrame(rows)
                .set_index('phenotype')
                .sort_values('mahal_mean', ascending=False))

    def print_summary(self, X_dis: np.ndarray, dis_pheno: np.ndarray) -> None:
        d    = self.distances(X_dis)
        keep = d <= self.threshold_
        print(f"Mahalanobis filter  P50={self.thr_p50_:.2f}  "
              f"P95={self.thr_p95_:.2f}  P99={self.thr_p99_:.2f}")
        print(f"Threshold (P{self.percentile}={self.threshold_:.2f}):  "
              f"retained {keep.sum()} / {len(keep)}  "
              f"removed {(~keep).sum()} ({(~keep).mean()*100:.1f}%)")

    # ── Visualization ─────────────────────────────────────────────

    def plot(self, X_dis: np.ndarray, dis_pheno: np.ndarray,
             save_dir=None, method_label: str = '') -> dict:
        """
        Figure 1: Violin + OOR heatmap
        Figure 2: UMAP (standardized) with outlier highlights

        Returns dict of matplotlib Figure objects.
        """
        import matplotlib.pyplot as plt
        import seaborn as sns
        from matplotlib.lines import Line2D
        from sklearn.preprocessing import StandardScaler

        self._check_fit()
        d    = self.distances(X_dis)
        keep = d <= self.threshold_
        dis_pheno = np.asarray(dis_pheno)
        unique_ph = np.unique(dis_pheno)
        p2c       = dict(zip(sorted(unique_ph),
                             sns.color_palette('tab20', len(unique_ph))))
        figs = {}

        # ── Figure 1: Violin + OOR heatmap ────────────────────────
        hc_q01  = np.percentile(self._X_hc_ref, 1,  axis=0)
        hc_q99  = np.percentile(self._X_hc_ref, 99, axis=0)
        oor_hc  = (self._X_hc_ref < hc_q01) | (self._X_hc_ref > hc_q99)
        oor_dis = (X_dis < hc_q01) | (X_dis > hc_q99)

        ph_order = (pd.DataFrame({'ph': dis_pheno, 'd': d})
                    .groupby('ph')['d'].mean()
                    .sort_values(ascending=False).index.tolist())

        short_cols = self._col_labels if hasattr(self, '_col_labels') else \
            [f'X{i}' for i in range(X_dis.shape[1])]

        oor_pct = pd.DataFrame(
            np.vstack([oor_dis[dis_pheno == ph].mean(axis=0) * 100 for ph in ph_order]),
            index=ph_order, columns=short_cols
        )
        hc_base = oor_hc.mean(axis=0) * 100

        fig1, (ax_l, ax_r) = plt.subplots(
            1, 2, figsize=(17, 8), gridspec_kw={'width_ratios': [1, 1.6]})
        df_vio = pd.DataFrame({'phenotype': dis_pheno, 'mahal': d})
        sns.violinplot(data=df_vio, y='phenotype', x='mahal', order=ph_order,
                       color='#aec6e8', orient='h', inner=None, cut=0,
                       linewidth=0.8, alpha=0.6, ax=ax_l)
        sns.stripplot(data=df_vio, y='phenotype', x='mahal', order=ph_order,
                      color='#2a6099', orient='h', size=2.5, alpha=0.5,
                      jitter=True, ax=ax_l)
        for thr, col, lbl in [(self.thr_p50_, 'steelblue', f'HC P50 ({self.thr_p50_:.1f})'),
                               (self.thr_p95_, 'orange',    f'HC P95 ({self.thr_p95_:.1f})'),
                               (self.thr_p99_, 'red',       f'HC P99 ({self.thr_p99_:.1f})')]:
            ax_l.axvline(thr, color=col, lw=1.2, ls='--', alpha=0.85, label=lbl)
        ax_l.legend(frameon=False, fontsize=8, loc='lower right')
        ax_l.set_xlabel('Mahalanobis Distance from HC')
        ax_l.set_ylabel('')
        ax_l.set_title('Covariate Extrapolation\n(sorted by mean Mahalanobis)')

        hm_data = pd.concat([
            oor_pct,
            pd.DataFrame([hc_base], index=['── HC baseline ──'], columns=short_cols)
        ])
        sns.heatmap(hm_data, ax=ax_r, cmap='Reds', annot=True, fmt='.0f',
                    vmin=0, vmax=40,
                    cbar_kws={'label': '% samples outside HC [P1–P99]', 'shrink': 0.7},
                    linewidths=0.3, linecolor='white')
        ax_r.set_title('% Disease Samples Outside HC P1–P99\n(per covariate; HC baseline at bottom)')
        fig1.suptitle('Covariate Distribution Overlap: Disease vs HC', y=1.01)
        plt.tight_layout()
        if save_dir:
            fig1.savefig(save_dir / 'covariate_overlap.png', bbox_inches='tight')
        figs['overlap'] = fig1

        # ── Figure 2: UMAP (standardized) ─────────────────────────
        try:
            from umap import UMAP
        except ImportError:
            print("umap-learn not installed — skipping UMAP figure")
            return figs

        scaler   = StandardScaler().fit(self._X_hc_ref)
        X_hc_sc  = scaler.transform(self._X_hc_ref)
        X_dis_sc = scaler.transform(X_dis)
        print('Computing UMAP on standardized bias covariates ...')
        emb     = UMAP(n_components=2, n_neighbors=30, min_dist=0.3,
                       metric='euclidean', random_state=42).fit_transform(
                           np.vstack([X_hc_sc, X_dis_sc]))
        emb_hc, emb_dis = emb[:len(X_hc_sc)], emb[len(X_hc_sc):]
        is_out = ~keep

        fig2, ax = plt.subplots(figsize=(10, 7))
        ax.scatter(emb_hc[:,0], emb_hc[:,1], c='lightgrey', s=8,
                   alpha=0.35, label='HC', zorder=1, edgecolors='none')
        for ph in sorted(unique_ph):
            m   = dis_pheno == ph
            sel = m & keep
            if sel.any():
                ax.scatter(emb_dis[sel,0], emb_dis[sel,1], c=[p2c[ph]], s=15,
                           alpha=0.75, label=ph, edgecolors='none', zorder=2)
        for ph in sorted(unique_ph):
            m   = dis_pheno == ph
            sel = m & is_out
            if sel.any():
                ax.scatter(emb_dis[sel,0], emb_dis[sel,1], c=[p2c[ph]], s=70,
                           alpha=0.95, marker='X', edgecolors='black',
                           linewidths=0.8, zorder=4)
        leg_h, leg_l = ax.get_legend_handles_labels()
        leg_h.append(Line2D([0],[0], marker='X', color='w',
                            markerfacecolor='grey', markeredgecolor='black',
                            markersize=9))
        leg_l.append(f'Outlier (Mahal > HC P{self.percentile})')
        ax.legend(handles=leg_h, labels=leg_l, fontsize=7, frameon=False,
                  bbox_to_anchor=(1.01, 1), loc='upper left', title='Phenotype')
        ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
        ax.set_title(f'UMAP of Bias Covariates (standardized)\n'
                     f'HC (grey) · Disease · Outlier X (>HC P{self.percentile})')
        plt.tight_layout()
        if save_dir:
            fig2.savefig(save_dir / 'covariate_umap.png', bbox_inches='tight')
        figs['umap'] = fig2
        print(f'Outliers shown: {is_out.sum()} samples')
        return figs

    def set_col_labels(self, labels: list) -> "MahalanobisFilter":
        """Optional: set short column labels for OOR heatmap."""
        self._col_labels = labels
        return self

    def _check_fit(self):
        if not self._fitted:
            raise RuntimeError("Call fit(X_hc) first.")

    def fit(self, X_hc: np.ndarray) -> "MahalanobisFilter":
        X_hc = np.asarray(X_hc, dtype=float)
        self._X_hc_ref = X_hc.copy()
        p    = X_hc.shape[1]
        self.mu_       = X_hc.mean(axis=0)
        self.cov_inv_  = inv(np.cov(X_hc.T) + np.eye(p) * self.reg)
        self.hc_dist_  = self._distances(X_hc)
        self.thr_p50_  = float(np.percentile(self.hc_dist_, 50))
        self.thr_p95_  = float(np.percentile(self.hc_dist_, 95))
        self.thr_p99_  = float(np.percentile(self.hc_dist_, 99))
        self.threshold_= float(np.percentile(self.hc_dist_, self.percentile))
        self._fitted   = True
        return self
