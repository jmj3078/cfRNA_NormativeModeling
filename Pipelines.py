import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import spearmanr
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
    def _compute_score(X, feat_vals, mode='spearman'):
        scores = []
        # Ensure feature values are float array
        feat_vals = np.array(feat_vals, dtype=float)

        for i in range(X.shape[0]):
            sample_expr = np.ravel(X[i, :])
            mask = sample_expr > 0 # Filter for detected genes only
            
            # Skip if too few genes detected
            if np.sum(mask) < 50:
                scores.append(0)
                continue
            
            valid_expr = sample_expr[mask]
            valid_feat = feat_vals[mask]

            if mode == 'lowess':
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
                bias_val = np.std(smoothed[:, 1]) 
                scores.append(bias_val)
            else:
                corr, _ = spearmanr(valid_expr, valid_feat)
                scores.append(corr if not np.isnan(corr) else 0)
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
            subset_X, subset_var[gc_col], mode='lowess'
        )
    else:
        print(f"  [Skip] GC column '{gc_col}' not found.")

    # 4. Calculate Length Bias
    if len_col in subset_var.columns:
        print("  > Computing Length bias score (Spearman)...")
        metrics_df['len_bias_score'] = _compute_score(
            subset_X, subset_var[len_col], mode='spearman'
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

    def analyze_pc_associations(self, n_pcs=5):
        if 'X_pca' not in self.adata.obsm: sc.tl.pca(self.adata)
        pc_df = pd.DataFrame(self.adata.obsm['X_pca'][:, :n_pcs], 
                             columns=[f'PC{i+1}' for i in range(n_pcs)], index=self.adata.obs_names)
        
        cont_vars = self.get_active_metrics()
        cat_vars = [self.cols['phenotype'], self.cols['batch']]
        assoc_matrix = pd.DataFrame(index=cont_vars + cat_vars, columns=pc_df.columns)

        # Continuous: Spearman
        for col in cont_vars:
            for pc in pc_df.columns:
                corr, _ = spearmanr(self.adata.obs[col], pc_df[pc])
                assoc_matrix.loc[col, pc] = abs(corr) if not np.isnan(corr) else 0
        
        # Categorical: ANOVA
        for col in cat_vars:
            if col not in self.adata.obs.columns or self.adata.obs[col].nunique() < 2:
                assoc_matrix.loc[col, :] = 0
                continue
            for pc in pc_df.columns:
                temp = pd.concat([self.adata.obs[col], pc_df[pc]], axis=1).dropna()
                temp.columns = ['G', 'V']
                try:
                    model = smf.ols('V ~ C(G)', data=temp).fit()
                    anova = sm.stats.anova_lm(model, typ=2)
                    assoc_matrix.loc[col, pc] = anova.loc['C(G)', 'sum_sq'] / anova['sum_sq'].sum()
                except:
                    assoc_matrix.loc[col, pc] = 0

        plt.figure(figsize=(8, len(assoc_matrix)*0.5 + 2))
        sns.heatmap(assoc_matrix.astype(float), annot=True, cmap='Reds', vmin=0, vmax=1, fmt='.2f')
        plt.title(f"PC - Variable Association")
        plt.show()

    def analyze_partial_correlation(self, n_pcs=5):
        print(f"--- [Analysis] Partial Correlation (Control Confounders) ---")
        target = self.cols['phenotype']
        confounders = self.get_active_metrics() + [self.cols['batch']]
        
        if 'X_pca' not in self.adata.obsm: sc.tl.pca(self.adata)
        pc_df = pd.DataFrame(self.adata.obsm['X_pca'][:, :n_pcs], columns=[f'PC{i+1}' for i in range(n_pcs)], index=self.adata.obs_names)
        data = pd.concat([pc_df, self.adata.obs[[target] + [c for c in confounders if c in self.adata.obs.columns]]], axis=1)
        
        # Numerize target
        if data[target].dtype == 'object' or data[target].dtype.name == 'category':
            data['target_num'] = pd.factorize(data[target])[0]
        else: data['target_num'] = data[target]

        results = {}
        for pc in pc_df.columns:
            # Raw Correlation
            raw = spearmanr(data[pc], data['target_num'])[0]
            
            # Partial Correlation
            conf_terms = [f"C({c})" if data[c].dtype.name in ['category', 'object'] else c for c in confounders if c in data.columns]
            if not conf_terms: 
                partial = raw
            else:
                try:
                    res_pc = smf.ols(f"{pc} ~ {' + '.join(conf_terms)}", data=data).fit().resid
                    res_tg = smf.ols(f"target_num ~ {' + '.join(conf_terms)}", data=data).fit().resid
                    partial = spearmanr(res_pc, res_tg)[0]
                except: partial = 0
            results[pc] = {'Raw': abs(raw), 'Partial': abs(partial)}
            
        pd.DataFrame(results).T.plot(kind='bar', figsize=(8, 4), colormap='tab20')
        plt.grid(alpha=0.1); plt.title(f"Correlation with {target}: Raw vs Controlled"); plt.show()

    def run_cascade_analysis(self, pcs_to_check=['PC1'], batch_method='combat', corr_threshold=0.95):


        print(f"--- [Analysis] Sequential Cascade Correction (Step-by-Step) ---")
        
        adata_temp = self.adata.copy()
        batch_col = self.cols['batch']
        
        # 1. 초기 메트릭 목록 확보 및 다중공선성 필터링
        initial_metrics = [m for m in self.get_active_metrics() if m != batch_col]
        metric_vars = []
        
        if len(initial_metrics) > 1:
            # 수치 안정성을 위해 결측치 임시 대치 후 상관계수 계산
            obs_temp = adata_temp.obs[initial_metrics].copy()
            for col in initial_metrics:
                obs_temp[col] = pd.to_numeric(obs_temp[col], errors='coerce').fillna(obs_temp[col].median())
            
            corr_matrix = obs_temp.corr(method='spearman').abs()
            upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
            to_drop = [column for column in upper.columns if any(upper[column] > corr_threshold)]
            metric_vars = [m for m in initial_metrics if m not in to_drop]
            
            if to_drop:
                print(f" [Warning] Dropped collinear metrics to prevent matrix singularity: {to_drop}")
        else:
            metric_vars = initial_metrics

        # 2. 단계(Steps) 정의
        steps = [('Raw Data', [], False)]
        for metric in metric_vars:
            steps.append((f"-{metric}", [metric], False))
        if batch_col in adata_temp.obs.columns:
            steps.append((f"-Batch ({batch_col})", [batch_col], True))
        
        # 3. 평가 대상 변수 (Trend Plot에 표시될 모든 변수)
        # 보정에서 제외된 다중공선성 변수들도 추세 확인을 위해 포함
        eval_vars = list(set(initial_metrics + [self.cols['phenotype']]))
        if batch_col in adata_temp.obs.columns:
            eval_vars.append(batch_col)
        eval_vars = [v for v in eval_vars if v in adata_temp.obs.columns]

        pc_results = {pc: {'val': [], 'var': []} for pc in pcs_to_check}
        
        # PCA Grid Layout 설정
        n_steps = len(steps)
        n_cols = 4
        n_rows = math.ceil(n_steps / n_cols)
        fig_pca, axes_pca = plt.subplots(n_rows, min(n_steps, n_cols), 
                                         figsize=(5 * min(n_steps, n_cols), 4.5 * n_rows), 
                                         constrained_layout=True)
        axes_flat = axes_pca.flatten() if n_steps > 1 else [axes_pca]
        
        failed_cascade = False 
        
        # --- Main Loop ---
        for idx, (step_name, drop_vars, is_batch_step) in enumerate(steps):
            ax = axes_flat[idx]
            error_msg = None
            
            # A. Correction & PC Recalculate
            if not failed_cascade:
                try:
                    if drop_vars:
                        if is_batch_step and batch_method == 'combat':
                            sc.pp.combat(adata_temp, key=drop_vars[0])
                        else:
                            sc.pp.regress_out(adata_temp, drop_vars)
                        adata_temp.X = np.nan_to_num(adata_temp.X, 0.0)

                    # 분산 체크: 최소 2개 이상의 유전자가 유의미한 분산을 가져야 PCA 가능
                    gene_vars = np.var(adata_temp.X, axis=0)
                    if np.sum(gene_vars > 1e-12) < 2:
                        error_msg = "Error: No Residual Remained"
                        failed_cascade = True
                    else:
                        sc.tl.pca(adata_temp, n_comps=5)
                        sc.pl.pca(adata_temp, color=self.cols['phenotype'], ax=ax, show=False, 
                                  size=150, alpha=0.6, title=f"Step {idx}: {step_name}", legend_loc='none')
                except Exception as e:
                    error_msg = f"Computation Failed:\n{str(e)[:30]}"
                    failed_cascade = True
            
            # B. 에러 발생 시 시각화 처리
            if failed_cascade or error_msg:
                display_msg = error_msg if error_msg else "Cascade Terminated"
                ax.text(0.5, 0.5, display_msg, ha='center', va='center', color='darkred', 
                        fontweight='bold', fontsize=10, transform=ax.transAxes)
                ax.set_title(f"Step {idx}: {step_name}", color='gray')
                ax.set_xticks([]); ax.set_yticks([])

            # C. Association Analysis (결과 저장 로직 - 에러 시 0으로 채움)
            for pc in pcs_to_check:
                step_vals = {'Step': step_name}
                current_var_ratio = 0
                
                if not failed_cascade and 'X_pca' in adata_temp.obsm:
                    pc_idx = int(pc.replace('PC','')) - 1
                    current_var_ratio = adata_temp.uns['pca']['variance_ratio'][pc_idx]
                    
                    pc_mat = pd.DataFrame(adata_temp.obsm['X_pca'], index=adata_temp.obs_names, 
                                          columns=[f'PC{i+1}' for i in range(5)])
                    data_step = pd.concat([pc_mat, adata_temp.obs[eval_vars]], axis=1)

                    for v in eval_vars:
                        val = 0
                        try:
                            if data_step[v].dtype.name in ['category', 'object']:
                                model = smf.ols(f"{pc} ~ C({v})", data=data_step).fit()
                                anova = sm.stats.anova_lm(model, typ=2)
                                val = anova.loc[f'C({v})', 'sum_sq'] / anova['sum_sq'].sum()
                            else:
                                val = abs(spearmanr(data_step[pc], data_step[v])[0])
                        except: val = 0
                        step_vals[v] = val if not np.isnan(val) else 0
                else:
                    # 에러 상태인 경우 모든 metric과의 연관성을 0으로 기록
                    for v in eval_vars: step_vals[v] = 0
                
                pc_results[pc]['val'].append(step_vals)
                pc_results[pc]['var'].append({'Step': step_name, 'Var_Ratio': current_var_ratio})

        # 빈 서브플롯 숨기기
        for j in range(n_steps, len(axes_flat)): axes_flat[j].axis('off')
        plt.show()

        # --- Trend Plot (Robust Implementation) ---
        print("Displaying Metric Evolution...")
        fig, axes = plt.subplots(len(pcs_to_check), 2, figsize=(18, 6 * len(pcs_to_check)))
        if len(pcs_to_check) == 1: axes = axes.reshape(1, -1)
        
        for i, pc in enumerate(pcs_to_check):
            df_vals = pd.DataFrame(pc_results[pc]['val']).set_index('Step').fillna(0)
            df_var = pd.DataFrame(pc_results[pc]['var']).set_index('Step').fillna(0)
            
            # Heatmap
            sns.heatmap(df_vals.T, cmap='Reds', annot=True, fmt='.2f', vmin=0, vmax=1, ax=axes[i,0])
            axes[i,0].set_title(f"{pc} Association Change Matrix", fontweight='bold')
            
            # Trend Line + PC Variance Bar
            ax_twin = axes[i,1].twinx()
            df_var.plot(kind='bar', y='Var_Ratio', ax=axes[i,1], color='lightgray', alpha=0.3, legend=False)
            axes[i,1].set_ylabel("PC Variance Ratio (Gray Bar)")
            
            colors = sns.color_palette("tab10", len(eval_vars))
            for j, v in enumerate(eval_vars):
                if v in df_vals.columns:
                    ax_twin.plot(range(len(df_vals)), df_vals[v], label=v, color=colors[j % 10], marker='o', lw=2, alpha=0.8)
            
            ax_twin.set_ylim(-0.05, 1.1)
            ax_twin.set_ylabel("Correlation / Eta-sq (Lines)")
            ax_twin.legend(bbox_to_anchor=(1.15, 1), loc='upper left', title="Metrics", fontsize='small')
            axes[i,1].set_xticklabels(df_vals.index, rotation=45, ha='right')
            axes[i,1].grid(axis='y', linestyle='--', alpha=0.7)
            
        plt.tight_layout(); plt.show()