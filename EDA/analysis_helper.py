import gc
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import median_abs_deviation
from scipy.stats import rankdata, t as t_dist, ks_2samp
from statsmodels.nonparametric.smoothers_lowess import lowess

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from config import PARAMS
from analysis_plot import (
    plot_pca_scree_and_bias,
    plot_pca_scatter_grid,
    plot_rda_unique_heatmap,
    plot_rda_variance_partition,
    plot_cascade_spaghetti,
    plot_normalization_partition,
    plot_hc_rda_results,
)

DEFAULT_METRICS = [
    "gc_bias_score",
    "len_bias_score",
    "platelet_score",
    "log1p_total_counts",
]

# ---------------------------------------------------------------------------
# Module-level sparse helpers
# ---------------------------------------------------------------------------

def _get_n80_sparse(matrix):
    if not sp.isspmatrix_csr(matrix):
        matrix = sp.csr_matrix(matrix)

    n80 = np.zeros(matrix.shape[0])
    for i in range(matrix.shape[0]):
        row_data = matrix.data[matrix.indptr[i]:matrix.indptr[i + 1]]
        row_data = row_data[row_data > 0]
        if len(row_data) == 0:
            continue
        cumsum = np.cumsum(np.sort(row_data)[::-1])
        if cumsum[-1] <= 0:
            continue
        n80[i] = np.argmax(cumsum >= cumsum[-1] * 0.8) + 1
    return n80


def _compute_score_sparse(X_csr, feat_vals, n_bins=None, loess_frac=None, outlier_pct=None):
    n_bins = n_bins or PARAMS["n_bins"]
    loess_frac = loess_frac or PARAMS["loess_frac"]
    outlier_pct = outlier_pct or PARAMS["outlier_pct"]

    feat_vals = np.array(feat_vals, dtype=float)
    scores = []

    for i in range(X_csr.shape[0]):
        row_data = X_csr.data[X_csr.indptr[i]:X_csr.indptr[i + 1]]
        row_indices = X_csr.indices[X_csr.indptr[i]:X_csr.indptr[i + 1]]
        pos_mask = row_data > 0
        valid_expr = row_data[pos_mask]
        valid_indices = row_indices[pos_mask]

        if len(valid_expr) < PARAMS["min_expressed"]:
            scores.append(0.0)
            continue

        q_thresh = np.percentile(valid_expr, outlier_pct)
        keep = valid_expr <= q_thresh
        final_expr = valid_expr[keep]
        final_feat = feat_vals[valid_indices[keep]]

        df_tmp = pd.DataFrame({"expr": final_expr, "feat": final_feat})
        try:
            df_tmp["bin"] = pd.qcut(df_tmp["feat"], q=n_bins, duplicates="drop")
        except ValueError:
            scores.append(0.0)
            continue

        bin_stats = (
            df_tmp.groupby("bin", observed=True)
            .agg({"expr": "median", "feat": "mean"})
            .dropna()
        )
        if len(bin_stats) < 2:
            scores.append(0.0)
            continue

        smoothed = lowess(bin_stats["expr"], bin_stats["feat"], frac=loess_frac, it=0)
        curve_disp = median_abs_deviation(smoothed[:, 1], scale="normal")
        total_disp = median_abs_deviation(final_expr, scale="normal")
        scores.append(curve_disp / total_disp if total_disp > 0 else 0.0)

    return np.array(scores)


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

def calculate_diversity_ratio(adata, gene_type_col="GeneType", coding_label="protein_coding"):
    X = adata.X
    print("Calculating NG80 (sparse)...")
    ng80 = _get_n80_sparse(X)

    print(f"Calculating NP80 (subset: {coding_label})...")
    if gene_type_col in adata.var.columns:
        is_coding = (adata.var[gene_type_col] == coding_label).values
        if np.sum(is_coding) == 0:
            print(f"  [Warning] No genes with type '{coding_label}'. NP80 set to 0.")
            np80 = np.zeros(adata.n_obs)
        else:
            np80 = _get_n80_sparse(X[:, is_coding])
    else:
        print(f"  [Warning] Column '{gene_type_col}' not found. NP80 set to NaN.")
        np80 = np.full(adata.n_obs, np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np80 / ng80
        ratio[ng80 == 0] = 0

    return pd.DataFrame(
        {"NG80": ng80, "NP80": np80, "NP80_NG80_ratio": ratio},
        index=adata.obs_names,
    )


def calculate_bias_metrics(
    adata,
    layer=None,
    gene_type_col="GeneType",
    target_type="protein_coding",
    gc_col="GC_Percent",
    len_col="log10_Length",
    platelet_col="is_platelet",
):
    print(f"--- Calculating Bfias Metrics (layer: {layer or 'X'}) ---")

    X_data = adata.layers[layer] if (layer and layer in adata.layers) else adata.X
    if not sp.issparse(X_data):
        X_data = sp.csr_matrix(X_data)
    elif not sp.isspmatrix_csr(X_data):
        X_data = X_data.tocsr()

    metrics_df = pd.DataFrame(index=adata.obs_names)

    if gene_type_col in adata.var.columns:
        coding_mask = (adata.var[gene_type_col] == target_type).values
        if not np.any(coding_mask):
            print(f"  [Warning] No genes for type '{target_type}'. Using all genes.")
            coding_mask = np.ones(adata.n_vars, dtype=bool)
    else:
        coding_mask = np.ones(adata.n_vars, dtype=bool)

    subset_X = X_data[:, coding_mask]
    subset_var = adata.var.iloc[coding_mask]

    if gc_col in subset_var.columns:
        print("  > GC bias (LOESS)...")
        metrics_df["gc_bias_score"] = _compute_score_sparse(subset_X, subset_var[gc_col])
    else:
        print(f"  [Skip] GC column '{gc_col}' not found.")

    if len_col in subset_var.columns:
        print("  > Length bias (LOESS)...")
        metrics_df["len_bias_score"] = _compute_score_sparse(subset_X, subset_var[len_col])
    else:
        print(f"  [Skip] Length column '{len_col}' not found.")

    if platelet_col in adata.var.columns:
        platelet_genes = adata.var_names[adata.var[platelet_col]].tolist()
        if platelet_genes:
            print(f"  > Platelet score ({len(platelet_genes)} genes)...")
            tmp = sc.AnnData(X=X_data, obs=adata.obs, var=adata.var)
            sc.tl.score_genes(tmp, gene_list=platelet_genes, score_name="platelet_score")
            metrics_df["platelet_score"] = tmp.obs["platelet_score"].values

    print("Done.\n")
    return metrics_df


# ---------------------------------------------------------------------------
# DataAnalysisPipeline
# ---------------------------------------------------------------------------

class DataAnalysisPipeline:

    DEFAULT_METRICS = DEFAULT_METRICS

    def __init__(
        self,
        adata,
        bias_metrics_df=None,
        phenotype_col="Type",
        batch_col="Batch_Granular",
        analysis_metrics=None,
        min_samples_per_study=None,
    ):
        self.min_samples = min_samples_per_study or PARAMS["min_study_samples"]
        self.cols = {"phenotype": phenotype_col, "batch": batch_col}
        self.target_metrics = analysis_metrics if analysis_metrics else self.DEFAULT_METRICS
        self.valid_studies = []

        keep = self._compute_keep_mask(adata)
        self.adata = adata[keep].copy()

        if bias_metrics_df is not None:
            new_cols = bias_metrics_df.columns.difference(self.adata.obs.columns)
            self.adata.obs = self.adata.obs.join(bias_metrics_df[new_cols])

        self._fix_batch_na()

    def _compute_keep_mask(self, adata):
        print("--- [Init] Preparing Metadata ---")
        pheno = self.cols["phenotype"]

        if "Author" not in adata.obs.columns:
            raise ValueError("'Author' column is missing in adata.obs.")

        keep = np.ones(adata.n_obs, dtype=bool)
        if pheno in adata.obs.columns:
            keep &= ~adata.obs[pheno].isna().values

        counts = adata.obs.loc[keep, "Author"].value_counts()
        self.valid_studies = counts[counts >= self.min_samples].index.tolist()
        dropped = counts[counts < self.min_samples].index.tolist()
        print(f"   Valid studies (>= {self.min_samples}): {len(self.valid_studies)}")
        if dropped:
            print(f"   [Skip] Too few samples: {dropped}")

        keep &= adata.obs["Author"].isin(self.valid_studies).values
        return keep

    def _fix_batch_na(self):
        batch = self.cols["batch"]
        if batch in self.adata.obs.columns and self.adata.obs[batch].isna().any():
            obs = self.adata.obs[batch]
            if hasattr(obs, "cat"):
                self.adata.obs[batch] = obs.cat.add_categories(["Unknown"])
            self.adata.obs[batch] = self.adata.obs[batch].fillna("Unknown")

    def _prepare_metadata(self):
        pass

    def get_active_metrics(self):
        return [m for m in self.target_metrics if m in self.adata.obs.columns]

    def set_layer(self, layer_name=None):
        if layer_name and layer_name in self.adata.layers:
            print(f"--- Layer: switching to '{layer_name}' ---")
            self.adata.X = self.adata.layers[layer_name]
        else:
            print("--- Layer: using current X ---")
        sc.pp.filter_genes(self.adata, min_cells=1)
        gc.collect()

    def _run_pca_for_study(self, adata_subset, use_hvg=False, n_top_genes=None, n_pcs=None):
        n_top_genes = n_top_genes or PARAMS["n_top_genes"]
        n_pcs = n_pcs or PARAMS["n_pcs"]

        tmp = adata_subset.copy()
        if use_hvg:
            sc.pp.highly_variable_genes(tmp, n_top_genes=n_top_genes, flavor="seurat")
            tmp = tmp[:, tmp.var.highly_variable].copy()

        if sp.issparse(tmp.X):
            tmp.X.data = np.nan_to_num(tmp.X.data, nan=0.0)
        else:
            tmp.X = np.nan_to_num(tmp.X, nan=0.0)

        sc.tl.pca(tmp, n_comps=n_pcs, svd_solver="arpack")
        return tmp

    def run_study_pca_diagnostics(self, use_hvg=False, save_dir=None):
        print(f"\n--- [Analysis] PCA Diagnostics per Study (HVG: {use_hvg}) ---")
        active_metrics = self.get_active_metrics()
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)

        for study in self.valid_studies:
            print(f"\n{'=' * 60}\n > {study}\n{'=' * 60}")
            adata_study = self.adata[self.adata.obs["Author"] == study].copy()
            adata_pca = self._run_pca_for_study(adata_study, use_hvg=use_hvg)
            var_ratios = adata_pca.uns["pca"]["variance_ratio"]

            slug = re.sub(r"[^\w.-]", "_", study)
            scree_path = os.path.join(save_dir, f"{slug}_scree_bias.png") if save_dir else None
            plot_pca_scree_and_bias(var_ratios, study, adata_pca.obs, active_metrics,
                                    save_path=scree_path)

            batch_col = self.cols["batch"]
            plot_batch_key = batch_col
            if batch_col in adata_pca.obs.columns:
                cats = (
                    adata_pca.obs[batch_col].cat.categories
                    if hasattr(adata_pca.obs[batch_col], "cat")
                    else adata_pca.obs[batch_col].unique()
                )
                simple_key = f"{batch_col}_simple"
                adata_pca.obs[simple_key] = (
                    adata_pca.obs[batch_col]
                    .map({orig: f"batch_{i + 1}" for i, orig in enumerate(cats)})
                    .astype("category")
                )
                plot_batch_key = simple_key

            plot_keys = list(dict.fromkeys(
                [k for k in active_metrics + [self.cols["phenotype"], plot_batch_key]
                 if k in adata_pca.obs.columns]
            ))
            if not plot_keys:
                del adata_study, adata_pca
                gc.collect()
                continue

            key_title_map = {plot_batch_key: "Batch (simplified)"} if plot_batch_key != batch_col else {}
            grid_path = os.path.join(save_dir, f"{slug}_pca_grid.png") if save_dir else None
            plot_pca_scatter_grid(adata_pca, plot_keys, var_ratios, key_title_map=key_title_map,
                                  save_path=grid_path)

            del adata_study, adata_pca
            gc.collect()

    @staticmethod
    def _encode_design_matrix(obs_df, target_vars, categorical_vars):
        col_parts = []
        col_groups = {}
        idx = 1

        for var in target_vars:
            if var in categorical_vars:
                dummies = pd.get_dummies(obs_df[var], drop_first=True, dtype=float)
                if dummies.shape[1] == 0:
                    print(f"   [Warning] '{var}' has only 1 unique level after encoding — skipped.")
                    continue
                col_parts.append(dummies.values)
                col_groups[var] = list(range(idx, idx + dummies.shape[1]))
                idx += dummies.shape[1]
            else:
                vals = obs_df[var].values.astype(float)
                std = vals.std()
                if std > 1e-12:
                    vals = (vals - vals.mean()) / std
                col_parts.append(vals.reshape(-1, 1))
                col_groups[var] = [idx]
                idx += 1

        n = len(obs_df)
        if col_parts:
            X = np.column_stack([np.ones(n)] + col_parts)
        else:
            X = np.ones((n, 1))
        return X, col_groups

    @staticmethod
    def _compute_matrix_r2(Y_expr, X):
        """Multivariate R² = SS_reg / SS_tot via QR. Returns (r2, sst)."""
        Y_c = Y_expr - Y_expr.mean(axis=0)
        ss_tot = float(np.einsum('ij,ij->', Y_c, Y_c))
        if X.shape[1] <= 1 or ss_tot < 1e-12:
            return 0.0, ss_tot
        Q, _ = np.linalg.qr(X, mode='reduced')
        QTYc = Q.T @ Y_c
        ss_reg = float(np.einsum('ij,ij->', QTYc, QTYc))
        return max(0.0, ss_reg / ss_tot), ss_tot

    @staticmethod
    def _adj_r2(Y_expr, X):
        """Adjusted R² to penalize DoF inflation. Returns (r2_adj, sst)."""
        r2, sst = DataAnalysisPipeline._compute_matrix_r2(Y_expr, X)
        n = Y_expr.shape[0]
        p = X.shape[1] - 1
        if n <= p + 1 or sst < 1e-12:
            return 0.0, sst
        return max(0.0, 1.0 - (1.0 - r2) * (n - 1) / (n - p - 1)), sst

    @staticmethod
    def _prepare_gene_matrix(adata_subset, use_hvg=False, n_top_genes=2000):
        """Dense float32 gene expression matrix. HVG selection optional. NaN → 0."""
        tmp = adata_subset.copy()
        if use_hvg:
            sc.pp.highly_variable_genes(tmp, n_top_genes=n_top_genes, flavor="seurat")
            tmp = tmp[:, tmp.var.highly_variable]
        X = tmp.X
        if sp.issparse(X):
            X = X.toarray()
        return np.nan_to_num(X.astype(np.float32), nan=0.0)

    @staticmethod
    def _run_partial_rda_core(Y_expr, obs_df, target_vars, categorical_vars,
                               phenotype_var, conf_vars):
        """
        Multivariate Partial RDA on full gene expression matrix Y_expr (n × p).
        Returns (unique_dict, partition_dict, r2_summary_dict).
        r2_summary_dict includes 'sst' for denominator tracking.
        """
        _nan_unique = {v: np.nan for v in target_vars}
        _nan_part = {"pheno_unique": np.nan, "conf_unique": np.nan, "shared": np.nan, "unexplained": np.nan}
        _nan_r2 = {"r2_all": np.nan, "r2_pheno_only": np.nan, "r2_conf_only": np.nan, "sst": np.nan}

        valid_mask = obs_df[target_vars].notna().all(axis=1)
        if valid_mask.sum() < len(target_vars) + 2:
            return _nan_unique, _nan_part, _nan_r2

        obs_sub = obs_df[valid_mask]
        Y_sub = Y_expr[valid_mask.values]
        pheno_in_vars = phenotype_var in target_vars

        enc = DataAnalysisPipeline._encode_design_matrix
        mr2 = DataAnalysisPipeline._adj_r2

        X_all, col_groups = enc(obs_sub, target_vars, categorical_vars)
        r2_all, sst = mr2(Y_sub, X_all)

        if pheno_in_vars and conf_vars:
            r2_pheno_only, _ = mr2(Y_sub, enc(obs_sub, [phenotype_var], categorical_vars)[0])
            r2_conf_only, _ = mr2(Y_sub, enc(obs_sub, conf_vars, categorical_vars)[0])
            partition = {
                "pheno_unique": max(0.0, r2_all - r2_conf_only),
                "conf_unique": max(0.0, r2_all - r2_pheno_only),
                "shared": max(0.0, r2_pheno_only + r2_conf_only - r2_all),
                "unexplained": max(0.0, 1.0 - r2_all),
            }
        elif pheno_in_vars:
            r2_pheno_only, _ = mr2(Y_sub, enc(obs_sub, [phenotype_var], categorical_vars)[0])
            r2_conf_only = np.nan
            partition = {"pheno_unique": np.nan, "conf_unique": np.nan,
                         "shared": np.nan, "unexplained": max(0.0, 1.0 - r2_all)}
        else:
            r2_pheno_only = np.nan
            r2_conf_only, _ = mr2(Y_sub, enc(obs_sub, conf_vars, categorical_vars)[0]) if conf_vars else (np.nan, sst)
            partition = {"pheno_unique": np.nan, "conf_unique": np.nan,
                         "shared": np.nan, "unexplained": max(0.0, 1.0 - r2_all)}

        unique_dict = {}
        for var in target_vars:
            if var not in col_groups:
                unique_dict[var] = np.nan
                continue
            remaining = [v for v in target_vars if v != var and v in col_groups]
            r2_minus_v, _ = mr2(Y_sub, enc(obs_sub, remaining, categorical_vars)[0]) if remaining else (0.0, sst)
            unique_dict[var] = max(0.0, r2_all - r2_minus_v)

        return unique_dict, partition, {
            "r2_all": r2_all, "r2_pheno_only": r2_pheno_only,
            "r2_conf_only": r2_conf_only, "sst": sst,
        }

    def analyze_partial_rda_per_study(self, vars, use_hvg=False,
                                       phenotype_var=None, batch_var=None,
                                       customized_x_order=None, save_dir=None):
        phenotype_var = phenotype_var or self.cols["phenotype"]
        batch_var = batch_var or self.cols["batch"]
        print("\n--- [Analysis] Partial RDA Variance Decomposition per Study ---")

        target_vars = [v for v in vars if v in self.adata.obs.columns]
        missing = [v for v in vars if v not in self.adata.obs.columns]
        if missing:
            print(f"   [Warning] Variables not found in obs and skipped: {missing}")
        if not target_vars:
            print("   [Error] No valid variables provided.")
            return None

        categorical_vars = {
            var for var in target_vars
            if not pd.api.types.is_numeric_dtype(self.adata.obs[var])
            or var in [phenotype_var, batch_var]
        }
        pheno_in_vars = phenotype_var in target_vars
        conf_vars = [v for v in target_vars if v != phenotype_var]
        unique_rows, partition_rows, r2_summary_rows = [], [], []

        for study in self.valid_studies:
            adata_study = self.adata[self.adata.obs["Author"] == study].copy()
            try:
                Y_expr = self._prepare_gene_matrix(adata_study, use_hvg=use_hvg)
                unique_d, part_d, r2_d = self._run_partial_rda_core(
                    Y_expr, adata_study.obs, target_vars, categorical_vars, phenotype_var, conf_vars,
                )
                if np.isnan(r2_d["r2_all"]):
                    n_valid = adata_study.obs[target_vars].notna().all(axis=1).sum()
                    print(f"   [Warning] {study}: too few complete-case samples ({n_valid}). Skipping.")
                for var in target_vars:
                    unique_rows.append({"Study": study, "Variable": var, "UniqueR2": unique_d[var]})
                partition_rows.append({"Study": study, **part_d})
                r2_summary_rows.append({"Study": study, **r2_d})
            except Exception as e:
                print(f"   [Error] {study}: {e}")
                for var in target_vars:
                    unique_rows.append({"Study": study, "Variable": var, "UniqueR2": np.nan})
                partition_rows.append({"Study": study, "pheno_unique": np.nan,
                                       "conf_unique": np.nan, "shared": np.nan, "unexplained": np.nan})
                r2_summary_rows.append({"Study": study, "r2_all": np.nan,
                                        "r2_pheno_only": np.nan, "r2_conf_only": np.nan})
            del adata_study
            gc.collect()

        df_unique = pd.DataFrame(unique_rows).pivot(
            index="Study", columns="Variable", values="UniqueR2"
        )
        priority = [c for c in [phenotype_var, batch_var] if c in df_unique.columns]
        other_cols = [c for c in df_unique.columns if c not in priority]
        df_unique = df_unique[priority + other_cols]
        if customized_x_order:
            df_unique = df_unique.reindex(columns=customized_x_order)

        df_partition = pd.DataFrame(partition_rows).set_index("Study")
        df_r2 = pd.DataFrame(r2_summary_rows).set_index("Study")

        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
        hvg_tag = "hvg" if use_hvg else "all"
        plot_rda_unique_heatmap(
            df_unique, use_hvg,
            save_path=os.path.join(save_dir, f"rda_unique_{hvg_tag}.png") if save_dir else None,
        )

        if pheno_in_vars and conf_vars:
            plot_rda_variance_partition(
                df_partition,
                title="Variance Partition: Phenotype vs Confounders per Study",
                save_path=os.path.join(save_dir, f"rda_partition_{hvg_tag}.png") if save_dir else None,
            )

        return {"unique_contributions": df_unique, "variance_partition": df_partition, "r2_summary": df_r2}

    def run_cascade_rda_per_study(self, vars, use_hvg=False,
                                   phenotype_var=None, batch_var=None, show_plots=False,
                                   save_path=None):
        phenotype_var = phenotype_var or self.cols["phenotype"]
        batch_var = batch_var or self.cols["batch"]
        print("\n--- [Analysis] Sequential Partial RDA (One-by-One Cascade) ---")

        target_vars = [v for v in vars if v in self.adata.obs.columns]
        if phenotype_var not in target_vars:
            print(f"   [Error] phenotype_var '{phenotype_var}' must be in vars.")
            return None

        conf_vars = [v for v in target_vars if v != phenotype_var]

        categorical_vars = {
            var for var in target_vars
            if not pd.api.types.is_numeric_dtype(self.adata.obs[var])
            or var in [phenotype_var, batch_var]
        }

        enc = self._encode_design_matrix
        all_metrics = {}

        for study in self.valid_studies:
            adata_work = self.adata[self.adata.obs["Author"] == study]
            obs = adata_work.obs

            valid_mask = obs[target_vars].notna().all(axis=1)
            if valid_mask.sum() < len(target_vars) + 2:
                print(f"   [Skip] {study}: too few complete-case samples.")
                continue

            obs_sub = obs[valid_mask]
            Y_sub = self._prepare_gene_matrix(adata_work[valid_mask].copy(), use_hvg=use_hvg)

            study_metrics = {}

            X_pheno, _ = enc(obs_sub, [phenotype_var], categorical_vars)
            pheno_r2_adj, _ = self._adj_r2(Y_sub, X_pheno)
            study_metrics["Raw"] = pheno_r2_adj

            cumulative_conf_vars = []
            for cv in conf_vars:
                cumulative_conf_vars.append(cv)

                X_conf, _ = enc(obs_sub, cumulative_conf_vars, categorical_vars)
                r2_adj_conf, _ = self._adj_r2(Y_sub, X_conf)

                current_all_vars = [phenotype_var] + cumulative_conf_vars
                X_all, _ = enc(obs_sub, current_all_vars, categorical_vars)
                r2_adj_all, _ = self._adj_r2(Y_sub, X_all)

                partial_pheno_r2 = max(0.0, r2_adj_all - r2_adj_conf)
                study_metrics[f"+ {cv}"] = partial_pheno_r2

            all_metrics[study] = study_metrics

        df_metrics = pd.DataFrame(all_metrics).T

        if show_plots or save_path:
            plot_cascade_spaghetti(df_metrics, save_path=save_path)

        return df_metrics

    def analyze_normalization_rda(self, vars, layers=None, use_hvg=False,
                                   phenotype_var=None, batch_var=None, show_plots=False,
                                   save_path=None):
        phenotype_var = phenotype_var or self.cols["phenotype"]
        batch_var = batch_var or self.cols["batch"]
        print("\n--- [Analysis] Normalization Comparison (Partial RDA) per Study ---")
        target_vars = [v for v in vars if v in self.adata.obs.columns]
        missing = [v for v in vars if v not in self.adata.obs.columns]
        if missing:
            print(f"   [Warning] Variables not found in obs and skipped: {missing}")
        if not target_vars:
            print("   [Error] No valid variables provided.")
            return None

        available_layers = list(self.adata.layers.keys())
        layers = [l for l in (layers or available_layers) if l in available_layers]
        if not layers:
            print("   [Error] No valid layers found.")
            return None

        categorical_vars = {
            var for var in target_vars
            if not pd.api.types.is_numeric_dtype(self.adata.obs[var])
            or var in [phenotype_var, batch_var]
        }
        conf_vars = [v for v in target_vars if v != phenotype_var]

        unique_rows, partition_rows, r2_rows = [], [], []

        for study in self.valid_studies:
            print(f" > {study}")
            adata_study = self.adata[self.adata.obs["Author"] == study].copy()
            sc.pp.filter_genes(adata_study, min_cells=1)
            for layer in layers:
                try:
                    adata_study.X = adata_study.layers[layer].copy()

                    Y_expr = self._prepare_gene_matrix(adata_study, use_hvg=use_hvg)

                    unique_d, part_d, r2_d = self._run_partial_rda_core(
                        Y_expr, adata_study.obs, target_vars, categorical_vars, phenotype_var, conf_vars,
                    )

                    for var in target_vars:
                        unique_rows.append({"Study": study, "Layer": layer,
                                            "Variable": var, "UniqueR2": unique_d[var]})
                    partition_rows.append({"Study": study, "Layer": layer, **part_d})
                    r2_rows.append({"Study": study, "Layer": layer, **r2_d})
                except Exception as e:
                    print(f"   [Error] {study} / {layer}: {e}")
            del adata_study
            gc.collect()

        df_unique_all = (pd.DataFrame(unique_rows).pivot_table(index=["Layer", "Study"], columns="Variable", values="UniqueR2", aggfunc="first"))
        df_partition_all = pd.DataFrame(partition_rows).set_index(["Layer", "Study"])
        df_r2_all = pd.DataFrame(r2_rows).set_index(["Layer", "Study"])

        layer_mapping = {
            "Raw": "Raw",
            "tpm": "TPM",
            "TPM_log2": "log2(TPM)",
            "CPM_log1p": "log1p(CPM)",
            "FPKM_log2": "log2(FPKM)",
            "TMM_log2": "log2(TMM)",
            "length_scaled": "Salmon(Length Scaled)",
            "scaled": "Salmon(Scaled)",
            "EDA_Full_All": "EDAseq (GC + Length)",
            "RUVg_Platelet_k1": "RUVg Platelet (k=1)",
            "RUVg_Platelet_k2": "RUVg Platelet (k=2)",
            "RUVg_Platelet_k3": "RUVg Platelet (k=3)",
            "Proposed_Full_k1": "EDAseq + RUVg (k=1)",
            "Proposed_Full_k2": "EDAseq + RUVg (k=2)",
            "Proposed_Full_k3": "EDAseq + RUVg (k=3)"
        }

        df_unique_all.rename(index=layer_mapping, level="Layer", inplace=True)
        df_partition_all.rename(index=layer_mapping, level="Layer", inplace=True)
        df_r2_all.rename(index=layer_mapping, level="Layer", inplace=True)

        logical_order = [
            "Raw", "TPM", "log2(TPM)", "log1p(CPM)", "log2(FPKM)",
            "log2(TMM)", "Salmon(Length Scaled)", "Salmon(Scaled)",
            "EDAseq (GC + Length)",
            "RUVg Platelet (k=1)", "RUVg Platelet (k=2)", "RUVg Platelet (k=3)",
            "EDAseq + RUVg (k=1)", "EDAseq + RUVg (k=2)", "EDAseq + RUVg (k=3)"
        ]

        existing_layers = [l for l in logical_order if l in df_partition_all.index.get_level_values("Layer")]
        df_unique_all = df_unique_all.reindex(existing_layers, level="Layer")
        df_partition_all = df_partition_all.reindex(existing_layers, level="Layer")
        df_r2_all = df_r2_all.reindex(existing_layers, level="Layer")

        if (show_plots or save_path) and not df_partition_all.empty and phenotype_var in target_vars:
            studies = df_partition_all.index.get_level_values("Study").unique()
            plot_normalization_partition(df_partition_all, studies, save_path=save_path)

        return {"unique_contributions": df_unique_all,
                "variance_partition": df_partition_all,
                "r2_summary": df_r2_all}


    def analyze_hc_partial_rda(self, vars, batch_col="Batch_ID",
                                hc_label="Healthy Control", phenotype_col="Phenotype_Processed",
                                use_hvg=False, layer=None, save_path=None):
        print("\n--- [Analysis] HC Cohort Partial RDA ---")
        adata_hc = self.adata[self.adata.obs[phenotype_col] == hc_label].copy()
        target_vars = [v for v in vars if v in adata_hc.obs.columns and v != phenotype_col]

        if adata_hc.n_obs < len(target_vars) + 2:
            print(f"   [Error] Too few HC samples ({adata_hc.n_obs}).")
            return None
        if layer and layer in adata_hc.layers:
            adata_hc.X = adata_hc.layers[layer].copy()

        sc.pp.filter_genes(adata_hc, min_cells=1)
        missing = [v for v in vars if v not in adata_hc.obs.columns]
        if missing:
            print(f"   [Warning] Variables not found and skipped: {missing}")
        if not target_vars:
            print("   [Error] No valid variables after filtering.")
            return None

        categorical_vars = {
            var for var in target_vars
            if not pd.api.types.is_numeric_dtype(adata_hc.obs[var]) or var == batch_col
        }
        Y_expr = self._prepare_gene_matrix(adata_hc, use_hvg=use_hvg)
        enc = DataAnalysisPipeline._encode_design_matrix
        adj_r2 = DataAnalysisPipeline._adj_r2
        valid_mask = adata_hc.obs[target_vars].notna().all(axis=1)
        obs_sub = adata_hc.obs[valid_mask]
        Y_sub = Y_expr[valid_mask.values]

        X_all, col_groups = enc(obs_sub, target_vars, categorical_vars)
        r2_all, sst = adj_r2(Y_sub, X_all)
        unique_dict = {}
        for var in target_vars:
            if var not in col_groups:
                unique_dict[var] = np.nan
                continue

            remaining = [v for v in target_vars if v != var and v in col_groups]
            if remaining:
                X_minus_v, _ = enc(obs_sub, remaining, categorical_vars)
                r2_minus_v, _ = adj_r2(Y_sub, X_minus_v)
            else:
                r2_minus_v = 0.0

            unique_dict[var] = max(0.0, r2_all - r2_minus_v)
        sr_unique = pd.Series(unique_dict).sort_values(ascending=True).dropna()
        plot_hc_rda_results(sr_unique, r2_all, batch_col, layer, unique_dict, save_path=save_path)

        return {"unique_contributions": sr_unique, "r2_all": r2_all}


def compute_gene_wise_bias_rda(
    adata,
    bias_metrics,
    layer="CPM_log1p",
    phenotype_col="Phenotype_Processed",
    target_labels="Healthy Control",
    group_name="HC",
    min_expressed_frac=0.1,
):
    print(f"\n--- [{group_name}] Vectorized Gene-wise Partial RDA ---")

    if target_labels is None:
        sample_mask = np.ones(adata.n_obs, dtype=bool)
    elif isinstance(target_labels, str):
        sample_mask = adata.obs[phenotype_col] == target_labels
    else:
        sample_mask = adata.obs[phenotype_col].isin(target_labels)

    adata_sub = adata[sample_mask].copy()
    obs_cols = [m for m in bias_metrics if m in adata_sub.obs.columns]
    valid_obs_mask = adata_sub.obs[obs_cols].notna().all(axis=1)
    adata_sub = adata_sub[valid_obs_mask].copy()

    n_samples = adata_sub.n_obs
    print(f"[{group_name}] Valid samples : {n_samples:,}")

    X_expr = adata_sub.layers[layer]
    if sp.issparse(X_expr):
        X_expr = X_expr.toarray()
    X_expr = X_expr.astype(np.float32)

    expressed_frac = (X_expr > 0).mean(axis=0)
    keep = expressed_frac >= min_expressed_frac
    Y = X_expr[:, keep]
    gene_names = adata_sub.var_names[keep].tolist()
    n_genes = len(gene_names)
    print(f"[{group_name}] Genes analyzed (≥{min_expressed_frac*100:.0f}% expressed) : {n_genes:,}")

    Y_c = Y - Y.mean(axis=0)
    sst = np.sum(Y_c ** 2, axis=0)
    valid_gene_mask = sst > 1e-12

    categorical_vars = {v for v in obs_cols if not pd.api.types.is_numeric_dtype(adata_sub.obs[v])}

    def _get_design_matrix(vars_list):
        return DataAnalysisPipeline._encode_design_matrix(adata_sub.obs, vars_list, categorical_vars)[0]

    def _vectorized_adj_r2(design_X):
        n, p_plus_1 = design_X.shape
        p = p_plus_1 - 1
        if n <= p + 1:
            return np.zeros(n_genes)

        Q, _ = np.linalg.qr(design_X, mode='reduced')
        ss_reg = np.sum((Q.T @ Y_c) ** 2, axis=0)

        r2 = np.zeros(n_genes)
        r2[valid_gene_mask] = ss_reg[valid_gene_mask] / sst[valid_gene_mask]

        r2_adj = 1.0 - (1.0 - r2) * (n - 1) / (n - p - 1)
        return np.clip(r2_adj, 0.0, 1.0)

    print(f"[{group_name}] Computing multivariate R² via orthogonal projection...")

    X_all = _get_design_matrix(obs_cols)
    r2_all = _vectorized_adj_r2(X_all)

    gene_records = {"Gene": gene_names, "Joint_R2_All_Biases": r2_all}
    sum_unique = np.zeros(n_genes)

    for bias in obs_cols:
        remaining = [v for v in obs_cols if v != bias]
        if remaining:
            X_minus = _get_design_matrix(remaining)
            r2_minus = _vectorized_adj_r2(X_minus)
        else:
            r2_minus = np.zeros(n_genes)

        unique_r2 = np.clip(r2_all - r2_minus, 0.0, 1.0)
        gene_records[f"Unique_{bias}"] = unique_r2
        sum_unique += unique_r2

    gene_records["Shared_Biases"] = np.clip(r2_all - sum_unique, 0.0, 1.0)
    gene_records["Unexplained"] = np.clip(1.0 - r2_all, 0.0, 1.0)

    df_detail = pd.DataFrame(gene_records)
    summary_data = []

    n_joint_contaminated = np.sum(r2_all > 0.10)
    summary_data.append({
        "Variance_Component": "ALL_BIASES_COMBINED (Joint R²)",
        "Max_R2": round(np.max(r2_all), 4),
        "Mean_R2": round(np.mean(r2_all), 4),
        "Genes_Highly_Biased": int(n_joint_contaminated),
        "Threshold": "> 10%"
    })

    for bias in obs_cols:
        unique_col = f"Unique_{bias}"
        vals = df_detail[unique_col].values
        n_contaminated = np.sum(vals > 0.05)

        summary_data.append({
            "Variance_Component": f"Unique: {bias}",
            "Max_R2": round(np.max(vals), 4),
            "Mean_R2": round(np.mean(vals), 4),
            "Genes_Highly_Biased": int(n_contaminated),
            "Threshold": "> 5%"
        })

    df_summary = pd.DataFrame(summary_data)

    print("=" * 75)
    print(" [Summary] Gene-level Contamination by Confounders")
    print("=" * 75)
    print(df_summary.to_string(index=False))
    print("=" * 75)

    return df_detail, df_summary
