import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import spearmanr, median_abs_deviation
import math


def calculate_diversity_ratio(adata, 
                              gene_type_col='GeneType', 
                              coding_label='protein_coding'):
    
    # Helper function: Calculate N80 for any given matrix
    def _get_n80(matrix):
        if matrix.shape[1] == 0: return np.zeros(matrix.shape[0])
        sorted_expr = np.sort(matrix, axis=1)[:, ::-1]
        cumsum = np.cumsum(sorted_expr, axis=1)
        total_counts = cumsum[:, -1]
        thresholds = total_counts * 0.8
        n80 = np.argmax(cumsum >= thresholds[:, None], axis=1) + 1
        n80[total_counts == 0] = 0
        return n80

    # ---------------------------------------------------------
    X = adata.X  
    print("Calculating NG80 (Total Library)...")
    ng80 = _get_n80(X)
    
    print(f"Calculating NP80 (Subset: {coding_label})...")
    if gene_type_col in adata.var.columns:
        is_coding = adata.var[gene_type_col] == coding_label
        if np.sum(is_coding) == 0:
            print(f" [Warning] No genes found with type '{coding_label}'. NP80 set to 0.")
            np80 = np.zeros(adata.n_obs)
        else:
            X_coding = X[:, is_coding.values] 
            np80 = _get_n80(X_coding)
    else:
        print(f" [Warning] Column '{gene_type_col}' not found. Cannot calculate NP80.")
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

    print(f"--- Calculating Bias Metrics (Target Layer: {layer if layer else 'X'}) ---")
    
    # 1. Setup Data Matrix (Dense array needed for calculations)
    if layer and layer in adata.layers:
        X_data = adata.layers[layer]
    else:
        X_data = adata.X
        
    # Ensure dense format if sparse
    if hasattr(X_data, "toarray"):
        X_data = X_data.toarray()
    
    # DataFrame to store results
    metrics_df = pd.DataFrame(index=adata.obs_names)

    # Helper: Compute correlation/bias score
    def _compute_score(X, feat_vals):
        scores = []
        # Ensure feature values are float array
        feat_vals = np.array(feat_vals, dtype=float)

        for i in range(X.shape[0]):
            sample_expr = np.ravel(X[i, :])
            expressed_mask = sample_expr > 0
            if np.sum(expressed_mask) < 50:
                scores.append(0)
                continue
            
            q99_thresh = np.percentile(sample_expr[expressed_mask], 99)
            final_mask = expressed_mask & (sample_expr <= q99_thresh)
                
            valid_expr = sample_expr[final_mask]
            valid_feat = feat_vals[final_mask]

            df_tmp = pd.DataFrame({'expr': valid_expr, 'feat': valid_feat})
            df_tmp['bin'] = pd.qcut(df_tmp['feat'], q=n_bins, duplicates='drop')
            bin_stats = df_tmp.groupby('bin', observed=True).agg({
                                'expr': 'median', 
                                'feat': 'mean'
                            }).dropna()
            if len(bin_stats) < 2:
                scores.append(0)
                continue
            
            smoothed = lowess(bin_stats['expr'], bin_stats['feat'], frac=0.7, it=0)
            curve_dispersion = median_abs_deviation(smoothed[:, 1], scale='normal')
            total_dispersion = median_abs_deviation(valid_expr, scale='normal')
            if total_dispersion > 0:
                bias_score = curve_dispersion / total_dispersion
            else:
                bias_score = 0 
            scores.append(bias_score)
        return np.array(scores)

    if gene_type_col in adata.var.columns:
        coding_mask = (adata.var[gene_type_col] == target_type).values
        if not np.any(coding_mask):
            print(f"  [Warning] No genes found for type '{target_type}'. Using all genes.")
            coding_mask = np.ones(adata.n_vars, dtype=bool)
    else:
        coding_mask = np.ones(adata.n_vars, dtype=bool)
        
    subset_X = X_data[:, coding_mask]
    subset_var = adata.var.iloc[coding_mask]

    if gc_col in subset_var.columns:
        print("  > Computing GC bias score (LOESS)...")
        metrics_df['gc_bias_score'] = _compute_score(
            subset_X, subset_var[gc_col]
        )
    else:
        print(f"  [Skip] GC column '{gc_col}' not found.")

    # 4. Calculate Length Bias
    if len_col in subset_var.columns:
        print("  > Computing Length bias score (LOESS)...")
        metrics_df['len_bias_score'] = _compute_score(
            subset_X, subset_var[len_col]
        )
    else:
        print(f"  [Skip] Length column '{len_col}' not found.")

    # 5. Calculate Platelet Score (Using scanpy score_genes)
    if platelet_col in adata.var.columns:
        platelet_genes = adata.var_names[adata.var[platelet_col]].tolist()
        if platelet_genes:
            print(f"  > Computing Platelet scores ({len(platelet_genes)} genes)...")
            # score_genes requires the full object, so we create a temporary one
            temp_adata = sc.AnnData(X=X_data, obs=adata.obs, var=adata.var)
            sc.tl.score_genes(temp_adata, gene_list=platelet_genes, score_name='platelet_score')
            metrics_df['platelet_score'] = temp_adata.obs['platelet_score'].values
        else:
             print("  [Skip] No platelet genes identified in var.")

    print("Calculation Done.\n")
    return metrics_df


class DataAnalysisPipeline:
    def __init__(self, adata, 
                 bias_metrics_df=None,
                 phenotype_col='Type', 
                 batch_col='Batch_Granular',
                 analysis_metrics=None):
        """
        Args:
            adata: AnnData object.
            bias_metrics_df: DataFrame containing calculated bias metrics (index must match adata.obs).
                             If None, assumes metrics are already in adata.obs.
            phenotype_col: Column name for biological condition of interest.
            batch_col: Column name for batch info.
            analysis_metrics: List of column names in obs to use for QC analysis. 
                              If None, defaults to ['gc_bias_score', 'len_bias_score', 'platelet_score', 'log1p_total_counts'].
        """
        self.adata = adata.copy()
        
        # Merge external metrics if provided
        if bias_metrics_df is not None:
            # Join carefully to avoid duplicates
            cols_to_use = bias_metrics_df.columns.difference(self.adata.obs.columns)
            self.adata.obs = self.adata.obs.join(bias_metrics_df[cols_to_use])
            
        self.cols = {
            'phenotype': phenotype_col,
            'batch': batch_col
        }
        
        # Define which metrics to track
        default_metrics = ['gc_bias_score', 'len_bias_score', 'platelet_score', 'log1p_total_counts']
        self.target_metrics = analysis_metrics if analysis_metrics else default_metrics
        
        # Validate metrics existence
        missing = [m for m in self.target_metrics if m not in self.adata.obs.columns]
        if missing:
            print(f"[Warning] The following metrics are missing from obs: {missing}")
        
        self._prepare_metadata()
        
    def _prepare_metadata(self):
        """Prepare categorical data for analysis (Unknown handling, etc.)"""
        print("--- [Init] Preparing Metadata ---")
        
        # Filter NaNs in phenotype
        if self.cols['phenotype'] in self.adata.obs.columns:
            n_orig = self.adata.n_obs
            self.adata = self.adata[~self.adata.obs[self.cols['phenotype']].isna()].copy()
            if self.adata.n_obs < n_orig:
                print(f"  Dropped {n_orig - self.adata.n_obs} samples with missing phenotype.")

        # Handle Batch NaNs
        if self.cols['batch'] in self.adata.obs.columns:
            if self.adata.obs[self.cols['batch']].isna().any():
                if self.adata.obs[self.cols['batch']].dtype.name == 'category':
                    self.adata.obs[self.cols['batch']] = self.adata.obs[self.cols['batch']].cat.add_categories(['Unknown'])
                self.adata.obs[self.cols['batch']] = self.adata.obs[self.cols['batch']].fillna('Unknown')
        
        print(f"  Active Metrics for Analysis: {[m for m in self.get_active_metrics()]}")

    def get_active_metrics(self):
        """Return list of metrics that actually exist in the current adata.obs"""
        return [m for m in self.target_metrics if m in self.adata.obs.columns]

    def set_layer(self, layer_name=None):
        if layer_name is None:
            print("--- Set Layer: Using current X (default) ---")
        elif layer_name in self.adata.layers:
            print(f"--- Set Layer: Switching to '{layer_name}' ---")
            self.adata.X = self.adata.layers[layer_name].copy()
        else:
            raise ValueError(f"Layer '{layer_name}' not found.")

        n_vars_before = self.adata.n_vars
        sc.pp.filter_genes(self.adata, min_cells=1)
        n_vars_after = self.adata.n_vars
        
        if n_vars_before > n_vars_after:
            print(f"   [Filtering] Dropped {n_vars_before - n_vars_after} genes with 0 expression in this layer.")
            print(f"   Remaining genes: {n_vars_after}")
        else:
            print("   [Filtering] No zero-expression genes found.")

        for key in ['X_pca', 'pca', 'PCs']:
            if key in self.adata.obsm: del self.adata.obsm[key]
            if key in self.adata.uns: del self.adata.uns[key]
            if key in self.adata.varm: del self.adata.varm[key]
        
        print(f"   PCA reset. Active layer ready.\n")
        
        
    def run_pca_diagnostics(self):
        if 'X_pca' not in self.adata.obsm: 
            sc.tl.pca(self.adata)
        
        print("--- [Analysis] PCA Diagnostics ---")
        fig_top, ax_top = plt.subplots(1, 2, figsize=(20, 5), gridspec_kw={'width_ratios': [3, 7]})
        
        # 1. Scree Plot
        var_ratios = self.adata.uns['pca']['variance_ratio']
        ax_top[0].plot(range(1, len(var_ratios)+1), var_ratios, 'o-k', alpha=0.7)
        ax_top[0].set_title("Scree Plot")
        ax_top[0].set_xlabel("PC"); ax_top[0].set_ylabel("Variance Ratio")
        
        # 2. Heatmap of Metrics
        active_metrics = self.get_active_metrics()
        if active_metrics:
            df_plot = self.adata.obs[active_metrics].copy()
            df_scaled = (df_plot - df_plot.mean()) / df_plot.std()
            sns.heatmap(df_scaled.T, cmap='RdBu_r', center=0, ax=ax_top[1], cbar_kws={'label': 'Z-score'})
            ax_top[1].set_title("Bias Metrics per Sample")
        
        plt.tight_layout(); plt.show()

        # ==============================================================================
        # [수정] 배치 레이블 단순화 (Legend 찌그러짐 방지)
        # ==============================================================================
        batch_col = self.cols['batch']
        plot_batch_key = batch_col  # 기본은 원본 사용

        if batch_col in self.adata.obs.columns:
            # (데이터 타입이 category면 .categories, 아니면 .unique() 사용)
            if hasattr(self.adata.obs[batch_col], 'cat'):
                unique_batches = self.adata.obs[batch_col].cat.categories
            else:
                unique_batches = self.adata.obs[batch_col].unique()
            
            # 2. 매핑 딕셔너리 생성 (Original Name -> batch_1, batch_2 ...)
            batch_map = {orig: f"batch_{i+1}" for i, orig in enumerate(unique_batches)}
            
            simple_batch_key = f"{batch_col}_simple"
            self.adata.obs[simple_batch_key] = self.adata.obs[batch_col].map(batch_map).astype('category')
            
            print(f"\n[Info] Simplified batch labels for plotting (Original -> Simple):")
            for orig, simple in batch_map.items():
                print(f"  * {simple} : {orig}")
            
            # 플롯 키를 단순화된 컬럼으로 교체
            plot_batch_key = simple_batch_key
        # ==============================================================================

        # 3. Scatter Plots
        print("\n--- [Analysis] PCA Scatter Plots ---")
        pc1_var = var_ratios[0]
        pc2_var = var_ratios[1]
        
        # plot_keys 구성 시 원본 batch 컬럼 대신, 단순화된 plot_batch_key 사용
        plot_keys = active_metrics + [self.cols['phenotype']]
        if plot_batch_key not in plot_keys: # 중복 방지
            plot_keys.append(plot_batch_key)
        
        n_cols = 3
        n_rows = math.ceil(len(plot_keys) / n_cols)
        fig_pca, axes_pca = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows))
        axes_pca = axes_pca.flatten()
        
        for i, key in enumerate(plot_keys):
            if key not in self.adata.obs.columns: continue
            
            is_numeric = pd.api.types.is_numeric_dtype(self.adata.obs[key])
            
            # 단순화된 배치 컬럼일 경우 제목에 알림
            title = key
            if key == plot_batch_key and key != batch_col:
                title = f"Batch (Simplified)"

            sc.pl.pca(
                self.adata, color=key, ax=axes_pca[i], show=False, 
                cmap='RdBu_r' if is_numeric else None,
                size=80, legend_loc='right margin', 
                title=title # 제목 설정
            )
            axes_pca[i].set_xlabel(f"PC1 ({pc1_var:.1%})")
            axes_pca[i].set_ylabel(f"PC2 ({pc2_var:.1%})")
            
        for j in range(i+1, len(axes_pca)): axes_pca[j].axis('off')
        plt.tight_layout(); plt.show()  

    def analyze_pc_impact(self, n_pcs=50):
        if 'X_pca' not in self.adata.obsm:
            sc.tl.pca(self.adata, n_comps=n_pcs)
        
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
            # Determine variable type for OLS formula
            is_cat = not pd.api.types.is_numeric_dtype(curr_obs) or (var in [self.cols['phenotype'], self.cols['batch'], 'Author'])
            
            for i in range(n_pcs_actual):
                y = curr_pcs[:, i]
                try:
                    # Unified R^2 estimation via OLS
                    df_tmp = pd.DataFrame({'y': y, 'x': curr_obs})
                    formula = 'y ~ C(x)' if is_cat else 'y ~ x'
                    res = smf.ols(formula, data=df_tmp).fit()
                    score = res.rsquared # This is Pearson r^2 for continuous variables
                except:
                    score = 0
                    
                # Weight by Eigenvalue (Explained Variance Ratio of the PC)
                total_impact += score * var_ratios[i]
                
            impact_results.append({'Variable': var, 'Impact_Score': total_impact})

        df_plot = pd.DataFrame(impact_results).sort_values('Impact_Score', ascending=False)
        fig, ax = plt.subplots(figsize=(10, len(df_plot) * 0.45 + 2))
        sns.barplot(
            data=df_plot, 
            x='Impact_Score', 
            y='Variable', 
            palette='magma', 
            ax=ax,
            edgecolor='0.3' # 막대 테두리를 살짝 주어 구분감 향상
        )

        max_val = df_plot['Impact_Score'].max()
        ax.set_xlim(0, max_val * 1.15)
        for i, val in enumerate(df_plot['Impact_Score']):
            # 0.0000인 경우 가독성을 위해 0으로 표시하거나 소수점 조정
            display_val = f"{val:.4f}" if val > 0 else "0.0000"
            ax.text(
                val + (max_val * 0.01), # 약간의 여백
                i, 
                display_val, 
                va='center', 
                ha='left', 
                fontsize=10, 
                fontweight='bold',
                color='#333333'
            )
            
        ax.set_title(r"Global Impact Score ($\sum VarRatio \times R^2$) - Top 50 PCs", 
                    fontsize=15, fontweight='bold', pad=20)
        ax.set_xlabel("Weighted Explained Variance ($R^2$)", fontsize=12, fontweight='bold')
        ax.set_ylabel("Variables", fontsize=12, fontweight='bold')
        
        # 그리드 설정 (X축 방향만 점선으로)
        ax.xaxis.grid(True, linestyle='--', alpha=0.7)
        ax.yaxis.grid(False) # Y축 그리드는 지워서 가로 막대 강조
        
        sns.despine(left=True, bottom=True) # 불필요한 테두리 제거
        plt.tight_layout()
        plt.show()
        
        return df_plot

    def run_cascade_analysis(self, n_pcs=50, batch_method='combat'):
        print(f"--- [Analysis] Global Cascade Analysis (Unified R² Tracking) ---")
        adata_temp = self.adata.copy()
        
        initial_metrics = [m for m in self.get_active_metrics() if m != self.cols['batch']]
        steps = [('Raw', [])]
        for m in initial_metrics: 
            steps.append((f"Reg__{m}", [m]))
        if self.cols['batch'] in adata_temp.obs.columns: 
            steps.append(("Corr_Batch", [self.cols['batch']]))
            
        eval_vars = list(dict.fromkeys(initial_metrics + [self.cols['phenotype'], self.cols['batch']]))
        eval_vars = [v for v in eval_vars if v in adata_temp.obs.columns]
        
        cascade_results = []
        failed_steps = set() # PCA가 불가능한 단계를 기록

        for step_name, drop_vars in steps:
            print(f"  > Processing: {step_name}")
            
            # 1. 보정 수행
            if drop_vars:
                try:
                    if step_name == "Corr_Batch" and batch_method == 'combat':
                        sc.pp.combat(adata_temp, key=drop_vars[0])
                    else:
                        sc.pp.regress_out(adata_temp, drop_vars)
                    adata_temp.X = np.nan_to_num(adata_temp.X, 0.0)
                except Exception as e:
                    print(f"    [Skip] Correction failed: {e}")
                    failed_steps.add(step_name)

            # 2. PCA 가능 여부 체크 (분산 확인)
            # 모든 유전자의 분산을 계산하여 0이 아닌 유전자가 최소 2개는 있어야 PCA 가능
            gene_vars = np.var(adata_temp.X, axis=0)
            if np.sum(gene_vars > 1e-12) < 2:
                print(f"    [Skip] No meaningful variance left in Step: {step_name}")
                failed_steps.add(step_name)
                
                # 모든 메트릭의 영향력을 0으로 기록하고 다음 단계로
                step_impacts = {'Step': step_name}
                for var in eval_vars: step_impacts[var] = 0
                cascade_results.append(step_impacts)
                continue

            # 3. PCA 및 R^2 계산
            try:
                sc.tl.pca(adata_temp, n_comps=n_pcs)
                pcs = adata_temp.obsm['X_pca']
                v_ratios = adata_temp.uns['pca']['variance_ratio']
                
                step_impacts = {'Step': step_name}
                for var in eval_vars:
                    total_r2 = 0
                    is_cat = not pd.api.types.is_numeric_dtype(adata_temp.obs[var]) or (var in [self.cols['phenotype'], self.cols['batch']])
                    
                    for i in range(min(n_pcs, pcs.shape[1])):
                        df_tmp = pd.DataFrame({'y': pcs[:, i], 'x': adata_temp.obs[var]})
                        formula = 'y ~ C(x)' if is_cat else 'y ~ x'
                        r2 = smf.ols(formula, data=df_tmp).fit().rsquared
                        total_r2 += r2 * v_ratios[i]
                    step_impacts[var] = total_r2
                
                cascade_results.append(step_impacts)
                
            except Exception as e:
                print(f"    [Error] PCA failed at {step_name}: {e}")
                failed_steps.add(step_name)
                step_impacts = {'Step': step_name}
                for var in eval_vars: step_impacts[var] = 0
                cascade_results.append(step_impacts)

        # 4. 시각화 (데이터가 있는 경우에만 그림)
        df_cascade = pd.DataFrame(cascade_results).set_index('Step')
        
        # Heatmap과 Line plot 구성
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 8))
        
        sns.heatmap(df_cascade.T, annot=True, fmt='.3f', cmap='YlGnBu', ax=ax1)
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
        ax1.set_title("Global Impact Score Heatmap (Weighted R²)")
        
        # Line plot에서는 실패한 단계를 제외하거나 끊어서 표시
        df_plot_trace = df_cascade.copy()
        # PCA가 실패한 지점은 선 그래프에서 시각적 왜곡을 줄이기 위해 NaN 처리하거나 그대로 0으로 둠
        df_plot_trace.plot(marker='o', ax=ax2, linewidth=2)
        ax2.set_title("Evolution of Variable Impacts (Weighted R² Trace)")
        ax2.set_ylabel("Weighted R² Score")
        ax2.set_xticks(range(len(df_cascade.index)))
        ax2.set_xticklabels(df_cascade.index, rotation=90)
        
        # 실패한 단계에 대한 텍스트 안내 추가
        for i, step in enumerate(df_cascade.index):
            if step in failed_steps:
                ax2.text(i, -0.05, "FAILED", color='red', ha='center', fontsize=8, fontweight='bold')
                
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Variables")
        plt.tight_layout(); plt.show()
        
        return df_cascade