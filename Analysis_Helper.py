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


import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import scipy.sparse as sp
import statsmodels.formula.api as smf
import gc

class DataAnalysisPipeline:
    def __init__(self, adata, bias_metrics_df=None, phenotype_col='Type', 
                 batch_col='Batch_Granular', analysis_metrics=None, min_samples_per_study=10):
        self.adata = adata.copy()
        self.min_samples = min_samples_per_study
        
        if bias_metrics_df is not None:
            cols_to_use = bias_metrics_df.columns.difference(self.adata.obs.columns)
            self.adata.obs = self.adata.obs.join(bias_metrics_df[cols_to_use])
            
        self.cols = {'phenotype': phenotype_col, 'batch': batch_col}
        default_metrics = ['gc_bias_score', 'len_bias_score', 'platelet_score', 'log1p_total_counts']
        self.target_metrics = analysis_metrics if analysis_metrics else default_metrics
        
        self.valid_studies = [] # 분석 가능한 Study 목록
        self._prepare_metadata()
        
    def _prepare_metadata(self):
        print("--- [Init] Preparing Metadata ---")
        if self.cols['phenotype'] in self.adata.obs.columns:
            self.adata = self.adata[~self.adata.obs[self.cols['phenotype']].isna()].copy()
        if self.cols['batch'] in self.adata.obs.columns:
            if self.adata.obs[self.cols['batch']].isna().any():
                if hasattr(self.adata.obs[self.cols['batch']], 'cat'):
                    self.adata.obs[self.cols['batch']] = self.adata.obs[self.cols['batch']].cat.add_categories(['Unknown'])
                self.adata.obs[self.cols['batch']] = self.adata.obs[self.cols['batch']].fillna('Unknown')
        
        if 'Author' not in self.adata.obs.columns:
            raise ValueError("'Author' column is missing in adata.obs.")
        study_counts = self.adata.obs['Author'].value_counts()
        self.valid_studies = study_counts[study_counts >= self.min_samples].index.tolist()
        dropped_studies = study_counts[study_counts < self.min_samples].index.tolist()
        print(f"   Valid Studies (>= {self.min_samples} samples): {len(self.valid_studies)}")
        if dropped_studies:
            print(f"   [Skip] Dropped Studies (too few samples): {dropped_studies}")
        self.adata = self.adata[self.adata.obs['Author'].isin(self.valid_studies)].copy()

    def get_active_metrics(self):
        return [m for m in self.target_metrics if m in self.adata.obs.columns]

    def set_layer(self, layer_name=None):
        if layer_name is not None and layer_name in self.adata.layers:
            print(f"--- Set Layer: Switching to '{layer_name}' ---")
            self.adata.X = self.adata.layers[layer_name].copy()
        else:
            print("--- Set Layer: Using current X (default) ---")
            
        sc.pp.filter_genes(self.adata, min_cells=1)
        gc.collect()

    def _run_pca_for_study(self, adata_subset, use_hvg=False, n_top_genes=2000, n_pcs=50):
        tmp_adata = adata_subset.copy()
        
        if use_hvg:
            sc.pp.highly_variable_genes(tmp_adata, n_top_genes=n_top_genes, flavor='seurat')
            tmp_adata = tmp_adata[:, tmp_adata.var.highly_variable].copy()
            
        if sp.issparse(tmp_adata.X):
            tmp_adata.X.data = np.nan_to_num(tmp_adata.X.data, nan=0.0)
        else:
            tmp_adata.X = np.nan_to_num(tmp_adata.X, nan=0.0)
            
        sc.tl.pca(tmp_adata, n_comps=n_pcs, svd_solver='arpack')
        return tmp_adata
    
    def run_study_pca_diagnostics(self, use_hvg=False):
        print(f"\n--- [Analysis] PCA Diagnostics per Study (HVG: {use_hvg}) ---")
        
        active_metrics = self.get_active_metrics()
        
        for study in self.valid_studies:
            print(f"\n{'='*60}")
            print(f" > Analyzing Study: {study}")
            print(f"{'='*60}")
            
            adata_study = self.adata[self.adata.obs['Author'] == study].copy()
            adata_pca = self._run_pca_for_study(adata_study, use_hvg=use_hvg)
            var_ratios = adata_pca.uns['pca']['variance_ratio']
            fig_top, ax_top = plt.subplots(1, 2, figsize=(20, 5), gridspec_kw={'width_ratios': [3, 7]})
            
            ax_top[0].plot(range(1, len(var_ratios)+1), var_ratios, 'o-k', alpha=0.7)
            ax_top[0].set_title(f"Scree Plot: {study}", fontweight='bold')
            ax_top[0].set_xlabel("PC"); ax_top[0].set_ylabel("Variance Ratio")
            ax_top[0].grid(True, linestyle='--', alpha=0.5)
            
            if active_metrics:
                df_plot = adata_pca.obs[active_metrics].copy()
                df_plot = df_plot.fillna(df_plot.mean())
                df_scaled = (df_plot - df_plot.mean()) / (df_plot.std() + 1e-9) 
                
                sns.heatmap(df_scaled.T, cmap='RdBu_r', center=0, ax=ax_top[1], cbar_kws={'label': 'Z-score'})
                ax_top[1].set_title(f"Bias Metrics per Sample: {study}", fontweight='bold')
            else:
                ax_top[1].axis('off')
                
            plt.tight_layout(); plt.show()

            batch_col = self.cols['batch']
            plot_batch_key = batch_col
            
            if batch_col in adata_pca.obs.columns:
                if hasattr(adata_pca.obs[batch_col], 'cat'):
                    unique_batches = adata_pca.obs[batch_col].cat.categories
                else:
                    unique_batches = adata_pca.obs[batch_col].unique()
                    
                batch_map = {orig: f"batch_{i+1}" for i, orig in enumerate(unique_batches)}
                simple_batch_key = f"{batch_col}_simple"
                adata_pca.obs[simple_batch_key] = adata_pca.obs[batch_col].map(batch_map).astype('category')
                plot_batch_key = simple_batch_key

            pc1_var, pc2_var = var_ratios[0], var_ratios[1]
            plot_keys = active_metrics + [self.cols['phenotype'], plot_batch_key]
            plot_keys = list(dict.fromkeys([k for k in plot_keys if k in adata_pca.obs.columns]))
            
            if not plot_keys:
                del adata_study, adata_pca; gc.collect()
                continue

            n_cols = 3
            n_rows = math.ceil(len(plot_keys) / n_cols)
            fig_pca, axes_pca = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows))
            
            if n_rows * n_cols == 1:
                axes_flat = [axes_pca]
            else:
                axes_flat = axes_pca.flatten()
            
            for i, key in enumerate(plot_keys):
                is_numeric = pd.api.types.is_numeric_dtype(adata_pca.obs[key])
                title = f"Batch (Simplified)" if key == plot_batch_key and key != batch_col else key
                
                sc.pl.pca(
                    adata_pca, color=key, ax=axes_flat[i], show=False, 
                    cmap='RdBu_r' if is_numeric else None, size=80, 
                    legend_loc='right margin', title=title
                )
                axes_flat[i].set_xlabel(f"PC1 ({pc1_var:.1%})")
                axes_flat[i].set_ylabel(f"PC2 ({pc2_var:.1%})")
                
            for j in range(i+1, len(axes_flat)): 
                axes_flat[j].axis('off')
                
            plt.tight_layout(); plt.show()  
            
            del adata_study, adata_pca; gc.collect()

    def analyze_pc_impact_per_study(self, n_pcs=50, use_hvg=False, customized_x_order=None):
        print(f"\n--- [Analysis] Covariate Impact (Weighted R²) per Study ---")
        target_vars = self.get_active_metrics() + [self.cols['phenotype'], self.cols['batch']]
        target_vars = list(dict.fromkeys([v for v in target_vars if v in self.adata.obs.columns]))
        
        results = []
        for study in self.valid_studies:
            adata_study = self.adata[self.adata.obs['Author'] == study].copy()
            try:
                adata_pca = self._run_pca_for_study(adata_study, use_hvg=use_hvg, n_pcs=n_pcs)
                pc_coords = adata_pca.obsm['X_pca']
                var_ratios = adata_pca.uns['pca']['variance_ratio']
                
                for var in target_vars:
                    valid = adata_pca.obs[var].notna()
                    if valid.sum() < 2: continue
                    
                    curr_obs = adata_pca.obs[valid][var]
                    is_cat = not pd.api.types.is_numeric_dtype(curr_obs) or (var in [self.cols['phenotype'], self.cols['batch']])
                    
                    total_impact = 0
                    for i in range(pc_coords.shape[1]):
                        try:
                            df_tmp = pd.DataFrame({'y': pc_coords[valid, i], 'x': curr_obs})
                            formula = 'y ~ C(x)' if is_cat else 'y ~ x'
                            r2 = smf.ols(formula, data=df_tmp).fit().rsquared
                            total_impact += r2 * var_ratios[i]
                        except: pass
                    results.append({'Study': study, 'Variable': var, 'Impact': total_impact})
            except Exception as e:
                print(f"   [Error] Impact analysis failed for {study}: {e}")
            
            del adata_study; gc.collect()
            
        df_res = pd.DataFrame(results).pivot(index='Study', columns='Variable', values='Impact').fillna(0)
        
        priority_cols = []
        for p_col in [self.cols['batch'], self.cols['phenotype']]:
            if p_col in df_res.columns:
                priority_cols.append(p_col)
                
        other_cols = [c for c in df_res.columns if c not in priority_cols]
        final_columns = priority_cols + other_cols
        
        df_res = df_res[final_columns]
        if customized_x_order:
            df_res = df_res.reindex(columns=customized_x_order)
        fig_width = max(10, len(final_columns) * 1.2)
        fig_height = max(6, len(self.valid_studies) * 0.6 + 2)
        
        plt.figure(figsize=(fig_width, fig_height))
        
        sns.heatmap(df_res, annot=True, fmt='.3f', cmap='Reds', linewidths=0.5)
        plt.title(f"Covariate Impact (Weighted R²) by Study\n(HVG: {use_hvg})", fontweight='bold', pad=15)
        plt.ylabel("Studies (Authors)", fontweight='bold')
        plt.xlabel("Covariates", fontweight='bold')
        plt.xticks(rotation=90, ha='right')
        
        plt.tight_layout()
        plt.show()
        
        return df_res

    def run_cascade_analysis_per_study(self, n_pcs=50, use_hvg=False):
        print(f"\n--- [Analysis] Cascade Correction Tracking per Study ---")
        
        initial_metrics = [m for m in self.get_active_metrics() if m != self.cols['batch']]
        steps = [('Raw', [])]
        for m in initial_metrics: steps.append((f"Reg__{m}", [m]))
        steps.append(("Corr_Batch", [self.cols['batch']]))
            
        eval_vars = list(dict.fromkeys(initial_metrics + [self.cols['phenotype'], self.cols['batch']]))
        all_study_results = {}

        for study in self.valid_studies:
            print(f" > Processing Cascade for Study: {study}")
            adata_study = self.adata[self.adata.obs['Author'] == study].copy()
            cascade_results = []
            
            for step_name, drop_vars in steps:
                if drop_vars:
                    if drop_vars[0] not in adata_study.obs.columns or adata_study.obs[drop_vars[0]].nunique() <= 1:
                        continue
                    try:
                        sc.pp.regress_out(adata_study, drop_vars)
                        if sp.issparse(adata_study.X):
                            adata_study.X.data = np.nan_to_num(adata_study.X.data, nan=0.0)
                        else:
                            adata_study.X = np.nan_to_num(adata_study.X, nan=0.0)
                    except Exception as e:
                        print(f"   [Skip] {step_name} failed: {e}")
                        continue
                
                try:
                    adata_pca = self._run_pca_for_study(adata_study, use_hvg=use_hvg, n_pcs=n_pcs)
                    pcs, v_ratios = adata_pca.obsm['X_pca'], adata_pca.uns['pca']['variance_ratio']
                    
                    step_impacts = {'Step': step_name}
                    for var in eval_vars:
                        if var not in adata_pca.obs.columns: continue
                        total_r2 = 0
                        is_cat = not pd.api.types.is_numeric_dtype(adata_pca.obs[var]) or (var in [self.cols['phenotype'], self.cols['batch']])
                        for i in range(pcs.shape[1]):
                            try:
                                df_tmp = pd.DataFrame({'y': pcs[:, i], 'x': adata_pca.obs[var]})
                                formula = 'y ~ C(x)' if is_cat else 'y ~ x'
                                total_r2 += smf.ols(formula, data=df_tmp).fit().rsquared * v_ratios[i]
                            except: pass
                        step_impacts[var] = total_r2
                    cascade_results.append(step_impacts)
                except Exception as e:
                    print(f"   [Error] PCA failed at {step_name}: {e}")
            
            df_cascade = pd.DataFrame(cascade_results).set_index('Step')
            all_study_results[study] = df_cascade
            del adata_study; gc.collect()
            custom_palette = [
                '#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', 
                '#EDC948', '#B07AA1', '#FF9DA7', '#9C755F', '#BAB0AC',
                '#0072B2', '#D55E00', '#009E73', "#FF0000", '#1F1F1F'
            ]
            plt.figure(figsize=(12, 8))
            df_cascade.plot(marker='o', linewidth=2, ax=plt.gca(), color=custom_palette)
            plt.title(f"Cascade Correction: {study}", fontweight='bold')
            plt.ylabel("Weighted R²"); plt.xlabel("Steps")
            plt.xticks(range(len(df_cascade.index)), df_cascade.index, rotation=90)
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout(); plt.show()
            
        return all_study_results