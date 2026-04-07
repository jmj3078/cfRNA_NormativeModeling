import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf

from scipy.stats import median_abs_deviation
from statsmodels.nonparametric.smoothers_lowess import lowess

import math
import scipy.sparse as sp
import gc

def calculate_diversity_ratio(adata, 
                              gene_type_col='GeneType', 
                              coding_label='protein_coding'):
    
    def _get_n80_sparse(matrix):
        if not sp.isspmatrix_csr(matrix):
            matrix = sp.csr_matrix(matrix)
            
        n_samples = matrix.shape[0]
        n80 = np.zeros(n_samples)
        
        for i in range(n_samples):
            start_idx = matrix.indptr[i]
            end_idx = matrix.indptr[i+1]
            row_data = matrix.data[start_idx:end_idx]
            
            if len(row_data) == 0:
                continue
                
            row_data = row_data[row_data > 0]
            if len(row_data) == 0:
                continue
                
            sorted_expr = np.sort(row_data)[::-1]
            cumsum = np.cumsum(sorted_expr)
            total_counts = cumsum[-1]
            
            if total_counts <= 0:
                continue
                
            threshold = total_counts * 0.8
            n80[i] = np.argmax(cumsum >= threshold) + 1
            
        return n80
    X = adata.X  
    print("Calculating NG80 (Total Library, Sparse Mode)...")
    ng80 = _get_n80_sparse(X)
    
    print(f"Calculating NP80 (Subset: {coding_label})...")
    if gene_type_col in adata.var.columns:
        is_coding = (adata.var[gene_type_col] == coding_label).values
        if np.sum(is_coding) == 0:
            print(f"  [Warning] No genes found with type '{coding_label}'. NP80 set to 0.")
            np80 = np.zeros(adata.n_obs)
        else:
            # CSR 슬라이싱은 여전히 희소 행렬을 유지함
            X_coding = X[:, is_coding] 
            np80 = _get_n80_sparse(X_coding)
    else:
        print(f"  [Warning] Column '{gene_type_col}' not found. Cannot calculate NP80.")
        np80 = np.full(adata.n_obs, np.nan)
        
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np80 / ng80
        ratio[ng80 == 0] = 0
    
    return pd.DataFrame({
        'NG80': ng80,
        'NP80': np80,
        'NP80_NG80_ratio': ratio
    }, index=adata.obs_names)
    
def calculate_bias_metrics(adata, layer=None, 
                           gene_type_col='GeneType', 
                           target_type='protein_coding',
                           gc_col='GC_Percent', 
                           len_col='log10_Length',
                           platelet_col='is_platelet',
                           n_bins=20):

    print(f"--- Calculating Bias Metrics (Target Layer: {layer if layer else 'X'}, Sparse Mode) ---")
    
    # 1. Setup Data Matrix
    X_data = adata.layers[layer] if layer and layer in adata.layers else adata.X
        
    # [수정됨] Dense 변환 코드 삭제. 대신 무조건 CSR로 보장.
    if not sp.issparse(X_data):
        print("  [Note] Input is dense. Converting to CSR sparse matrix internally.")
        X_data = sp.csr_matrix(X_data)
    elif not sp.isspmatrix_csr(X_data):
        X_data = X_data.tocsr()
    
    metrics_df = pd.DataFrame(index=adata.obs_names)

    def _compute_score_sparse(X_csr, feat_vals):
        scores = []
        feat_vals = np.array(feat_vals, dtype=float)

        n_samples = X_csr.shape[0]
        for i in range(n_samples):
            start_idx = X_csr.indptr[i]
            end_idx = X_csr.indptr[i+1]
            row_data = X_csr.data[start_idx:end_idx]
            row_indices = X_csr.indices[start_idx:end_idx]
            pos_mask = row_data > 0
            valid_expr = row_data[pos_mask]
            valid_indices = row_indices[pos_mask]
            
            if len(valid_expr) < 50:
                scores.append(0.0)
                continue
            
            # Dense 배열에서의 퍼센타일과 달리, 0이 배제된 상태이므로 수학적 왜곡 방지
            q99_thresh = np.percentile(valid_expr, 99)
            final_mask = valid_expr <= q99_thresh
            final_expr = valid_expr[final_mask]
            final_feat = feat_vals[valid_indices[final_mask]]
            df_tmp = pd.DataFrame({'expr': final_expr, 'feat': final_feat})
            
            try:
                df_tmp['bin'] = pd.qcut(df_tmp['feat'], q=n_bins, duplicates='drop')
            except ValueError:
                scores.append(0.0)
                continue
                
            bin_stats = df_tmp.groupby('bin', observed=True).agg({
                                'expr': 'median', 
                                'feat': 'mean'
                            }).dropna()
            
            if len(bin_stats) < 2:
                scores.append(0.0)
                continue
            
            smoothed = lowess(bin_stats['expr'], bin_stats['feat'], frac=0.7, it=0)
            curve_dispersion = median_abs_deviation(smoothed[:, 1], scale='normal')
            total_dispersion = median_abs_deviation(final_expr, scale='normal')
            
            if total_dispersion > 0:
                bias_score = curve_dispersion / total_dispersion
            else:
                bias_score = 0.0
            scores.append(bias_score)
            
        return np.array(scores)

    if gene_type_col in adata.var.columns:
        coding_mask = (adata.var[gene_type_col] == target_type).values
        if not np.any(coding_mask):
            print(f"  [Warning] No genes found for type '{target_type}'. Using all genes.")
            coding_mask = np.ones(adata.n_vars, dtype=bool)
    else:
        coding_mask = np.ones(adata.n_vars, dtype=bool)
        
    # CSR 행렬에서의 열 슬라이싱
    subset_X = X_data[:, coding_mask]
    subset_var = adata.var.iloc[coding_mask]

    if gc_col in subset_var.columns:
        print("  > Computing GC bias score (LOESS)...")
        metrics_df['gc_bias_score'] = _compute_score_sparse(
            subset_X, subset_var[gc_col]
        )
    else:
        print(f"  [Skip] GC column '{gc_col}' not found.")

    if len_col in subset_var.columns:
        print("  > Computing Length bias score (LOESS)...")
        metrics_df['len_bias_score'] = _compute_score_sparse(
            subset_X, subset_var[len_col]
        )
    else:
        print(f"  [Skip] Length column '{len_col}' not found.")

    if platelet_col in adata.var.columns:
        platelet_genes = adata.var_names[adata.var[platelet_col]].tolist()
        if platelet_genes:
            print(f"  > Computing Platelet scores ({len(platelet_genes)} genes)...")
            temp_adata = sc.AnnData(X=X_data, obs=adata.obs, var=adata.var)
            sc.tl.score_genes(temp_adata, gene_list=platelet_genes, score_name='platelet_score')
            metrics_df['platelet_score'] = temp_adata.obs['platelet_score'].values
        else:
             print("  [Skip] No platelet genes identified in var.")

    print("Calculation Done.\n")
    return metrics_df


class DataAnalysisPipeline:
    def __init__(self, adata, bias_metrics_df=None, phenotype_col='Type', 
                 batch_col='Batch_Granular', analysis_metrics=None):
        """
        7GB 이상의 대용량 Sparse Matrix 환경에 최적화된 분석 파이프라인.
        메모리 피크를 관리하며 PCA 진단 및 공변량 영향력을 평가합니다.
        """
        self.adata = adata.copy()
        
        # 외부 메트릭 병합 (중복 방지)
        if bias_metrics_df is not None:
            cols_to_use = bias_metrics_df.columns.difference(self.adata.obs.columns)
            self.adata.obs = self.adata.obs.join(bias_metrics_df[cols_to_use])
            
        self.cols = {
            'phenotype': phenotype_col,
            'batch': batch_col
        }
        
        default_metrics = ['gc_bias_score', 'len_bias_score', 'platelet_score', 'log1p_total_counts']
        self.target_metrics = analysis_metrics if analysis_metrics else default_metrics
        
        self._prepare_metadata()
        
    def _prepare_metadata(self):
        print("--- [Init] Preparing Metadata ---")
        
        # Phenotype 결측치 샘플 제거
        if self.cols['phenotype'] in self.adata.obs.columns:
            n_orig = self.adata.n_obs
            self.adata = self.adata[~self.adata.obs[self.cols['phenotype']].isna()].copy()
            if self.adata.n_obs < n_orig:
                print(f"   Dropped {n_orig - self.adata.n_obs} samples with missing phenotype.")

        # Batch 결측치 처리 (Unknown 할당)
        if self.cols['batch'] in self.adata.obs.columns:
            if self.adata.obs[self.cols['batch']].isna().any():
                if self.adata.obs[self.cols['batch']].dtype.name == 'category':
                    self.adata.obs[self.cols['batch']] = self.adata.obs[self.cols['batch']].cat.add_categories(['Unknown'])
                self.adata.obs[self.cols['batch']] = self.adata.obs[self.cols['batch']].fillna('Unknown')
        
        print(f"   Active Metrics for Analysis: {self.get_active_metrics()}")

    def get_active_metrics(self):
        return [m for m in self.target_metrics if m in self.adata.obs.columns]

    def set_layer(self, layer_name=None):
        """특정 Layer로 X를 교체하고 Sparse 구조를 유지하며 유전자 필터링 수행"""
        if layer_name is None:
            print("--- Set Layer: Using current X (default) ---")
        elif layer_name in self.adata.layers:
            print(f"--- Set Layer: Switching to '{layer_name}' ---")
            # Sparse Matrix 포맷 유지 보장
            self.adata.X = self.adata.layers[layer_name].copy()
        else:
            raise ValueError(f"Layer '{layer_name}' not found.")

        # 발현이 전혀 없는 유전자 제거 (PCA 연산 에러 방지)
        n_vars_before = self.adata.n_vars
        sc.pp.filter_genes(self.adata, min_cells=1)
        
        if n_vars_before > self.adata.n_vars:
            print(f"   [Filtering] Dropped {n_vars_before - self.adata.n_vars} genes.")
        
        # 기존 PCA 결과 초기화
        for key in ['X_pca', 'pca', 'PCs']:
            if key in self.adata.obsm: del self.adata.obsm[key]
            if key in self.adata.uns: del self.adata.uns[key]
            if key in self.adata.varm: del self.adata.varm[key]
        
        gc.collect() # 메모리 즉시 정리
        print(f"   PCA reset. Active layer ready.\n")

    def run_pca_diagnostics(self):
        """Sparse 최적화 솔버를 사용한 PCA 진단 및 메트릭 Heatmap 출력"""
        # Sparse Matrix에 가장 효율적인 arpack 솔버 사용
        if 'X_pca' not in self.adata.obsm: 
            sc.tl.pca(self.adata, svd_solver='arpack')
        
        print("--- [Analysis] PCA Diagnostics ---")
        fig_top, ax_top = plt.subplots(1, 2, figsize=(20, 5), gridspec_kw={'width_ratios': [3, 7]})
        
        # Scree Plot
        var_ratios = self.adata.uns['pca']['variance_ratio']
        ax_top[0].plot(range(1, len(var_ratios)+1), var_ratios, 'o-k', alpha=0.7)
        ax_top[0].set_title("Scree Plot")
        
        # Bias Metrics Heatmap
        active_metrics = self.get_active_metrics()
        if active_metrics:
            df_plot = self.adata.obs[active_metrics].copy()
            df_scaled = (df_plot - df_plot.mean()) / df_plot.std()
            sns.heatmap(df_scaled.T, cmap='RdBu_r', center=0, ax=ax_top[1], cbar_kws={'label': 'Z-score'})
            ax_top[1].set_title("Bias Metrics per Sample")
        plt.tight_layout(); plt.show()

        # 배치 레이블 단순화 (Legend 겹침 방지)
        batch_col = self.cols['batch']
        plot_batch_key = batch_col
        if batch_col in self.adata.obs.columns:
            unique_batches = self.adata.obs[batch_col].cat.categories if hasattr(self.adata.obs[batch_col], 'cat') else self.adata.obs[batch_col].unique()
            batch_map = {orig: f"batch_{i+1}" for i, orig in enumerate(unique_batches)}
            simple_batch_key = f"{batch_col}_simple"
            self.adata.obs[simple_batch_key] = self.adata.obs[batch_col].map(batch_map).astype('category')
            plot_batch_key = simple_batch_key

        # PCA Scatter Plots
        print("\n--- [Analysis] PCA Scatter Plots ---")
        pc1_var, pc2_var = var_ratios[0], var_ratios[1]
        plot_keys = active_metrics + [self.cols['phenotype'], plot_batch_key]
        plot_keys = list(dict.fromkeys([k for k in plot_keys if k in self.adata.obs.columns]))
        
        n_cols = 3
        n_rows = math.ceil(len(plot_keys) / n_cols)
        fig_pca, axes_pca = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows))
        axes_flat = axes_pca.flatten()
        
        for i, key in enumerate(plot_keys):
            is_numeric = pd.api.types.is_numeric_dtype(self.adata.obs[key])
            sc.pl.pca(self.adata, color=key, ax=axes_flat[i], show=False, 
                      cmap='RdBu_r' if is_numeric else None, size=80, 
                      legend_loc='right margin', title=f"Batch (Simple)" if key == plot_batch_key else key)
            axes_flat[i].set_xlabel(f"PC1 ({pc1_var:.1%})")
            axes_flat[i].set_ylabel(f"PC2 ({pc2_var:.1%})")
            
        for j in range(i+1, len(axes_flat)): axes_flat[j].axis('off')
        plt.tight_layout(); plt.show()  

    def analyze_pc_impact(self, n_pcs=50):
        """각 변수가 PCA 주성분에 미치는 영향력을 가중 R^2로 계산"""
        if 'X_pca' not in self.adata.obsm:
            sc.tl.pca(self.adata, n_comps=n_pcs, svd_solver='arpack')
        
        n_pcs_actual = min(n_pcs, self.adata.obsm['X_pca'].shape[1])
        pc_coords = self.adata.obsm['X_pca'][:, :n_pcs_actual]
        var_ratios = self.adata.uns['pca']['variance_ratio'][:n_pcs_actual]
        
        target_vars = self.get_active_metrics() + [self.cols['phenotype'], self.cols['batch']]
        if 'Author' in self.adata.obs.columns: target_vars.append('Author')
        target_vars = list(dict.fromkeys([v for v in target_vars if v in self.adata.obs.columns]))
        
        impact_results = []
        for var in target_vars:
            total_impact = 0
            valid = self.adata.obs[var].notna()
            if valid.sum() < 2: continue
            
            curr_obs = self.adata.obs[valid][var]
            curr_pcs = pc_coords[valid]
            is_cat = not pd.api.types.is_numeric_dtype(curr_obs) or (var in [self.cols['phenotype'], self.cols['batch'], 'Author'])
            
            for i in range(n_pcs_actual):
                y = curr_pcs[:, i]
                try:
                    df_tmp = pd.DataFrame({'y': y, 'x': curr_obs})
                    formula = 'y ~ C(x)' if is_cat else 'y ~ x'
                    res = smf.ols(formula, data=df_tmp).fit()
                    total_impact += res.rsquared * var_ratios[i]
                except: continue
                    
            impact_results.append({'Variable': var, 'Impact_Score': total_impact})

        df_plot = pd.DataFrame(impact_results).sort_values('Impact_Score', ascending=False)
        
        # 시각화 (Barplot)
        fig, ax = plt.subplots(figsize=(10, len(df_plot) * 0.45 + 2))
        sns.barplot(data=df_plot, x='Impact_Score', y='Variable', palette='magma', ax=ax, edgecolor='0.3')
        ax.set_title("Global Impact Score (Weighted R^2)")
        plt.tight_layout(); plt.show()
        return df_plot

    def run_cascade_analysis(self, n_pcs=50, batch_method='combat'):
        """보정 단계별 공변량 영향력 변화 추적 (Sparse-safe Regression/ComBat)"""
        print(f"--- [Analysis] Global Cascade Analysis ---")
        adata_temp = self.adata.copy()
        
        initial_metrics = [m for m in self.get_active_metrics() if m != self.cols['batch']]
        steps = [('Raw', [])]
        for m in initial_metrics: steps.append((f"Reg__{m}", [m]))
        if self.cols['batch'] in adata_temp.obs.columns: steps.append(("Corr_Batch", [self.cols['batch']]))
            
        eval_vars = list(dict.fromkeys(initial_metrics + [self.cols['phenotype'], self.cols['batch']]))
        eval_vars = [v for v in eval_vars if v in adata_temp.obs.columns]
        
        cascade_results = []
        failed_steps = set()

        for step_name, drop_vars in steps:
            print(f"   > Processing: {step_name}")
            if drop_vars:
                try:
                    if step_name == "Corr_Batch" and batch_method == 'combat':
                        # Combat은 Sparse를 미지원하므로 임시 Dense 변환 (메모리 주의)
                        if sp.issparse(adata_temp.X):
                            print("     [Note] Converting to dense for ComBat...")
                            adata_temp.X = adata_temp.X.toarray()
                        sc.pp.combat(adata_temp, key=drop_vars[0])
                        adata_temp.X = sp.csr_matrix(adata_temp.X) # 즉시 Sparse 복구
                    else:
                        # regress_out은 Sparse를 지원하지만 연산량이 많음
                        sc.pp.regress_out(adata_temp, drop_vars)
                    
                    # Sparse .data 직접 접근을 통한 NaN 고속 처리
                    if sp.issparse(adata_temp.X):
                        adata_temp.X.data = np.nan_to_num(adata_temp.X.data, nan=0.0)
                    else:
                        adata_temp.X = np.nan_to_num(adata_temp.X, nan=0.0)
                except Exception as e:
                    print(f"     [Skip] Correction failed: {e}")
                    failed_steps.add(step_name)

            # 분산이 없는 레이어에서 PCA 에러 방지 (Sparse 전용 연산)
            if sp.issparse(adata_temp.X):
                mean_sq = adata_temp.X.power(2).mean(axis=0)
                mean = adata_temp.X.mean(axis=0)
                var = np.array(mean_sq - np.square(mean)).flatten()
            else:
                var = np.var(adata_temp.X, axis=0)

            if np.sum(var > 1e-12) < 2:
                failed_steps.add(step_name)
                cascade_results.append({'Step': step_name, **{v: 0 for v in eval_vars}})
                continue

            try:
                sc.tl.pca(adata_temp, n_comps=n_pcs, svd_solver='arpack')
                pcs, v_ratios = adata_temp.obsm['X_pca'], adata_temp.uns['pca']['variance_ratio']
                step_impacts = {'Step': step_name}
                for var in eval_vars:
                    total_r2 = 0
                    is_cat = not pd.api.types.is_numeric_dtype(adata_temp.obs[var]) or (var in [self.cols['phenotype'], self.cols['batch']])
                    for i in range(min(n_pcs, pcs.shape[1])):
                        df_tmp = pd.DataFrame({'y': pcs[:, i], 'x': adata_temp.obs[var]})
                        formula = 'y ~ C(x)' if is_cat else 'y ~ x'
                        total_r2 += smf.ols(formula, data=df_tmp).fit().rsquared * v_ratios[i]
                    step_impacts[var] = total_r2
                cascade_results.append(step_impacts)
            except:
                failed_steps.add(step_name)
                cascade_results.append({'Step': step_name, **{v: 0 for v in eval_vars}})
            
            gc.collect() # 각 단계 종료 후 메모리 해제

        df_cascade = pd.DataFrame(cascade_results).set_index('Step')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 8))
        
        # 1) Heatmap
        sns.heatmap(df_cascade.T, annot=True, fmt='.3f', cmap='YlGnBu', ax=ax1)
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, fontsize=10) # 90도 회전
        ax1.set_title("Global Impact Score Heatmap (Weighted R²)", fontsize=14, fontweight='bold')
        
        # 2) Line Trace Plot
        df_plot_trace = df_cascade.copy()
        df_plot_trace.plot(marker='o', ax=ax2, linewidth=2)
        
        ax2.set_title("Impacts Across Variable Correction Steps", fontsize=14, fontweight='bold')
        ax2.set_ylabel("Weighted R² Score", fontsize=12)
        ax2.set_xlabel("Correction Steps", fontsize=12)
        
        ax2.set_xticks(range(len(df_cascade.index)))
        ax2.set_xticklabels(df_cascade.index, rotation=90, fontsize=10) 
        
        for i, step in enumerate(df_cascade.index):
            if step in failed_steps:
                ax2.text(i, -0.02, "FAILED", color='red', ha='center', 
                         fontsize=9, fontweight='bold', rotation=90)
                
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Variables", fontsize=10)
        ax2.grid(axis='both', linestyle='--', alpha=0.5) 
        
        plt.tight_layout()
        plt.show()
        
        return df_cascade