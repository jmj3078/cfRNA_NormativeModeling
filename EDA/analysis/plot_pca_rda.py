import math
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from analysis.plot_utils import PALETTE, _save


def plot_hc_variance_and_heatmap(adata, bias_metrics, hc_label="Healthy Control",
                                  phenotype_col="Phenotype_Processed", n_pcs=50, save_path=None):
    adata_hc = adata[adata.obs[phenotype_col] == hc_label]
    v_ratios = adata_hc.uns["pca"]["variance_ratio"][:n_pcs]

    fig, ax = plt.subplots(1, 2, figsize=(20, 5), gridspec_kw={"width_ratios": [3, 7]})
    ax[0].plot(range(1, len(v_ratios) + 1), v_ratios, "o-k", alpha=0.8, markersize=5)
    ax[0].set_title("Scree Plot")
    ax[0].set_ylabel("Explained Variance Ratio")
    ax[0].set_xlabel("Principal Component")
    ax[0].grid(True, linestyle="--", alpha=0.4)

    num_metrics = [m for m in bias_metrics if m in adata_hc.obs.columns
                   and pd.api.types.is_numeric_dtype(adata_hc.obs[m])]
    df_m = adata_hc.obs[num_metrics].copy()
    df_scaled = (df_m - df_m.mean()) / (df_m.std() + 1e-9)
    sns.heatmap(df_scaled.T, cmap="RdBu_r", center=0, ax=ax[1],
                cbar_kws={"label": "Z-score"}, xticklabels=False)
    ax[1].set_xlabel(f"Samples (n={len(adata_hc)})")
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()


def plot_hc_pca_grid(adata, plot_vars, hc_label="Healthy Control",
                      phenotype_col="Phenotype_Processed", save_path=None):
    adata_hc = adata[adata.obs[phenotype_col] == hc_label].copy()
    v_ratios = adata_hc.uns["pca"]["variance_ratio"]
    pc1_v, pc2_v = v_ratios[0], v_ratios[1]

    n_cols = 4
    n_rows = math.ceil(len(plot_vars) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5.5 * n_rows))
    axes_flat = axes.flatten() if n_rows * n_cols > 1 else [axes]

    for i, var in enumerate(plot_vars):
        if var not in adata_hc.obs.columns:
            continue
        is_num = pd.api.types.is_numeric_dtype(adata_hc.obs[var])
        sc.pl.pca(adata_hc, color=var, ax=axes_flat[i], show=False,
                  cmap="RdBu_r" if is_num else None, size=70, alpha=0.7,
                  palette="Spectral" if not is_num else None,
                  legend_loc="right margin", wspace=0.8, hspace=0.8)
        axes_flat[i].set_title(var)
        axes_flat[i].set_xlabel(f"PC1 ({pc1_v:.1%})")
        axes_flat[i].set_ylabel(f"PC2 ({pc2_v:.1%})")

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis("off")
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()


def plot_pca_scree_and_bias(var_ratios, study, obs_df, active_metrics, save_path=None):
    """Scree plot + bias metrics heatmap for a single study."""
    fig, ax = plt.subplots(1, 2, figsize=(20, 5), gridspec_kw={"width_ratios": [3, 7]})
    ax[0].plot(range(1, len(var_ratios) + 1), var_ratios, "o-k", alpha=0.7)
    ax[0].set_title(f"Scree: {study}")
    ax[0].set_xlabel("PC")
    ax[0].set_ylabel("Variance Ratio")
    ax[0].grid(True, linestyle="--", alpha=0.5)

    if active_metrics:
        df_plot = obs_df[active_metrics].copy().fillna(0)
        df_scaled = (df_plot - df_plot.mean()) / (df_plot.std() + 1e-9)
        sns.heatmap(df_scaled.T, cmap="RdBu_r", center=0, ax=ax[1],
                    cbar_kws={"label": "Z-score"})
        ax[1].set_title(f"Bias Metrics: {study}")
    else:
        ax[1].axis("off")
    ax[1].set_xlabel(f"Samples (n={len(obs_df)})")
    ax[1].set_xticks([])
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()


def plot_pca_scatter_grid(adata_pca, plot_keys, var_ratios, key_title_map=None, save_path=None):
    """PCA scatter grid colored by each metadata key."""
    pc1_var, pc2_var = var_ratios[0], var_ratios[1]
    n_cols = 3
    n_rows = math.ceil(len(plot_keys) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows))
    axes_flat = [axes] if n_rows * n_cols == 1 else axes.flatten()

    for i, key in enumerate(plot_keys):
        is_numeric = pd.api.types.is_numeric_dtype(adata_pca.obs[key])
        title = (key_title_map or {}).get(key, key)
        sc.pl.pca(adata_pca, color=key, ax=axes_flat[i], show=False,
                  cmap="RdBu_r" if is_numeric else None, size=80,
                  legend_loc="right margin", title=title)
        axes_flat[i].set_xlabel(f"PC1 ({pc1_var:.1%})")
        axes_flat[i].set_ylabel(f"PC2 ({pc2_var:.1%})")

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis("off")

    plt.tight_layout()
    _save(fig, save_path)
    plt.show()


def plot_rda_unique_heatmap(df_unique, use_hvg, save_path=None):
    """Heatmap of per-variable unique R² contributions across studies."""
    fig_w = max(10, len(df_unique.columns) * .8)
    fig_h = max(6, len(df_unique) * 0.4 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(df_unique, annot=True, fmt=".3f", cmap="Reds", annot_kws={"size": 11},
                linewidths=0.5, mask=df_unique.isna(), ax=ax)
    sns.heatmap(df_unique.isna(), cmap=["white", "lightgrey"],
                cbar=False, ax=ax, mask=~df_unique.isna(), linewidths=0.5)
    ax.set_title(f"Per-Variable Unique Contribution (Partial RDA)\n(HVG: {use_hvg})", pad=15)
    ax.set_ylabel("Studies (Authors)")
    ax.set_xlabel("Covariates")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()


def plot_rda_variance_partition(df_partition, title, x_label="Studies (Authors)",
                                 figsize=None, save_path=None):
    """Stacked horizontal bar chart of variance partition."""
    part_cols = ["pheno_unique", "conf_unique", "shared", "unexplained"]
    part_labels = ["Phenotype Unique", "Covariate Unique", "Shared", "Unexplained"]
    part_colors = ["steelblue", "coral", "mediumseagreen", "lightgrey"]
    df_plot = df_partition[part_cols].fillna(0)
    fig_h = figsize[1] if figsize else max(5, len(df_plot) * 0.4 + 2)
    fig_w = figsize[0] if figsize else 10
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    bottom = np.zeros(len(df_plot))
    for col, label, color in zip(part_cols, part_labels, part_colors):
        vals = df_plot[col].values
        ax.barh(df_plot.index, vals, left=bottom, color=color, label=label, edgecolor="black")
        bottom += vals
    ax.set_xlabel("Proportion of Total Variance")
    ax.set_ylabel(x_label)
    ax.set_title(title, pad=15)
    ax.set_xlim(0, 1)
    ax.legend(loc="center left", frameon=False, bbox_to_anchor=(1.01, 0.5))
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()


def plot_cascade_spaghetti(df_metrics, palette=None, save_path=None):
    """Spaghetti plot of phenotype R² trajectory across sequential confounders."""
    if df_metrics.empty:
        return
    palette = palette or PALETTE
    all_steps = df_metrics.columns.tolist()
    x = np.arange(len(all_steps))

    fig, ax = plt.subplots(figsize=(max(10, len(all_steps) * 0.4), 6))
    for i, study in enumerate(df_metrics.index):
        color = palette[i % len(palette)]
        ax.plot(x, df_metrics.loc[study], marker="o", markersize=6,
                linewidth=2, alpha=0.85, color=color, label=study)

    ax.set_xticks(x)
    ax.set_xticklabels(all_steps, rotation=90, ha="right")
    ax.set_ylabel("Adjusted Partial R²")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", title="Studies", frameon=False)
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()


def plot_normalization_partition(df_partition_all, studies, save_path=None):
    """Per-study variance partition bars across normalization layers."""
    n_cols = min(5, len(studies))
    n_rows = math.ceil(len(studies) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4.5, n_rows * 6.5), sharey=True)
    axes_flat = [axes] if n_rows * n_cols == 1 else axes.flatten()

    part_cols = ["pheno_unique", "conf_unique", "shared", "unexplained"]
    part_labels = ["Phenotype Unique", "Covariate Unique", "Shared", "Unexplained"]
    part_colors = ["steelblue", "coral", "mediumseagreen", "lightgrey"]
    for i, study in enumerate(studies):
        ax = axes_flat[i]
        df_study = df_partition_all.xs(study, level="Study")[part_cols].fillna(0)
        layers_plot = df_study.index.tolist()
        x = np.arange(len(layers_plot))
        bottom = np.zeros(len(layers_plot))

        for col, color in zip(part_cols, part_colors):
            vals = df_study[col].values
            ax.bar(x, vals, bottom=bottom, color=color, width=0.6, edgecolor="black", linewidth=0.5)
            bottom += vals
        ax.set_xticks(x)
        ax.set_xticklabels(layers_plot, rotation=90, ha="right")
        if i % n_cols == 0:
            ax.set_ylabel("Proportion of Variance")
        ax.set_ylim(0, 1)
        ax.grid(axis="y", linestyle="--", alpha=0.4)

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    patches = [mpatches.Patch(color=c, label=l) for c, l in zip(part_colors, part_labels)]
    fig.legend(handles=patches, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.05), frameon=False)
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()


def plot_normalization_unique_r2_summary(r2_results, phenotype_var="Phenotype_Processed",
                                          rank_on="EDAseq|RUVg", save_path=None):
    """Phenotype vs technical-covariate unique R2 across normalization layers.

    Reads the partial-RDA unique-contribution table (layer x study x variable).
    Left: composition of total unique R2. Right: technical / phenotype fold.
    Covariates are ordered by their contribution in the corrected layers
    (`rank_on`), so the ordering is not dominated by depth in the raw layers.
    """
    if isinstance(r2_results, str):
        r2_results = pd.read_csv(r2_results, sep="\t")
    df = r2_results.groupby("Layer", sort=False).mean(numeric_only=True)

    tech = [c for c in df.columns if c != phenotype_var]
    rank_rows = df.loc[df.index.str.contains(rank_on)] if rank_on else df
    order = df.loc[rank_rows.index, tech].mean().sort_values(ascending=False).index.tolist()

    comp = df[[phenotype_var] + order].fillna(0.0)
    total = comp.sum(axis=1)
    frac = comp.div(total, axis=0)
    fold = (total - comp[phenotype_var]) / comp[phenotype_var]

    cols = ["#4E79A7"] + [PALETTE[1:][i % (len(PALETTE) - 1)] for i in range(len(order))]
    y = np.arange(len(df))[::-1]
    corrected = df.index.str.contains(rank_on) if rank_on else np.ones(len(df), bool)

    fs = 19
    lw = 1.8
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 0.62 * len(df) * 2 + 5.5))

    left = np.zeros(len(df))
    for col, c in zip(comp.columns, cols):
        ax1.barh(y, frac[col].values, left=left, color=c, edgecolor="white",
                 linewidth=0.6, height=0.72, label=col.replace(phenotype_var, "Phenotype"))
        left += frac[col].values
    for yi, p in enumerate(frac[phenotype_var].values):
        ax1.text(p + 0.015, y[yi], f"{p * 100:.1f}%", va="center", fontsize=fs - 2,
                 fontweight="bold", color="#4E79A7",
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    ax1.set_yticks(y)
    ax1.set_yticklabels(df.index, fontsize=fs)
    for tick, c in zip(ax1.get_yticklabels(), corrected):
        tick.set_color("#1F1F1F" if c else "#7A7A7A")
    ax1.set_xlim(0, 1)
    ax1.set_xlabel("Share of total unique R² (partial RDA)", fontsize=fs)
    ax1.set_title("Phenotype vs technical covariates", fontsize=fs + 1)
    ax1.tick_params(axis="x", labelsize=fs - 1)
    ax1.grid(axis="x", linestyle="--", alpha=0.4)
    for spine in ax1.spines.values():
        spine.set_linewidth(lw)

    ax2.barh(y, fold.values, color=["#4E79A7" if c else "#BAB0AC" for c in corrected],
             edgecolor="black", linewidth=0.6, height=0.72)
    for yi, v in enumerate(fold.values):
        ax2.text(v * 1.06, y[yi], f"{v:.0f}×", va="center", fontsize=fs - 2)
    ax2.set_xscale("log")
    ax2.set_xlim(1, fold.max() * 3)
    ax2.axvline(1, color="black", linewidth=lw)
    ax2.set_yticks(y)
    ax2.set_yticklabels(df.index, fontsize=fs)
    for tick, c in zip(ax2.get_yticklabels(), corrected):
        tick.set_color("#1F1F1F" if c else "#7A7A7A")
    ax2.set_xlabel("Technical / phenotype unique R²  (log scale)", fontsize=fs)
    ax2.set_title("Technical variance dominance", fontsize=fs + 1)
    ax2.tick_params(axis="x", labelsize=fs - 1)
    ax2.grid(axis="x", linestyle="--", alpha=0.4)
    for spine in ax2.spines.values():
        spine.set_linewidth(lw)

    fig.legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.08), frameon=False, fontsize=fs - 1)
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()
    return pd.DataFrame({"pheno_unique_R2": comp[phenotype_var], "total_unique_R2": total,
                         "pheno_share": frac[phenotype_var], "fold_technical": fold})


def plot_hc_rda_results(sr_unique, r2_all, batch_col, layer, unique_dict, save_path=None):
    """Two-panel HC RDA result: per-variable bar + joint variance partition."""
    stem, ext = (os.path.splitext(save_path) if save_path else (None, ".png"))
    ext = ext or ".png"

    fig1, ax1 = plt.subplots(figsize=(7, max(4, len(sr_unique) * 0.45 + 1.5)))
    ax1.barh(sr_unique.index, sr_unique.values, color="steelblue", edgecolor="black")
    ax1.set_xlabel("Unique R² (Partial RDA)")
    ax1.axvline(0, color="black", linewidth=0.8)
    ax1.grid(axis="x", linestyle="--", alpha=0.4)
    plt.tight_layout()
    _save(fig1, f"{stem}_bar{ext}" if stem else None)
    plt.show()

    conf_unique_sum = sum(v for v in unique_dict.values() if not np.isnan(v))
    shared = max(0.0, r2_all - conf_unique_sum)
    unexplained = max(0.0, 1.0 - r2_all)
    part_vals = [conf_unique_sum, shared, unexplained]
    part_labels = ["Covariate Unique (sum)", "Shared", "Unexplained"]
    part_colors = ["coral", "mediumseagreen", "lightgrey"]
    fig2, ax2 = plt.subplots(figsize=(6.5, 1.0))
    bottom = 0.0
    for val, label, color in zip(part_vals, part_labels, part_colors):
        ax2.barh(["HC"], [val], left=[bottom], color=color, label=label, height=0.3, edgecolor="black")
        bottom += val
    ax2.set_xlabel("Proportion of Total Variance")
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.6), ncol=3, frameon=False)
    plt.tight_layout()
    _save(fig2, f"{stem}_partition{ext}" if stem else None)
    plt.show()
