import math
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
import statsmodels.formula.api as smf
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import make_pipeline

try:
    from skbio.stats.distance import permanova
    from skbio import DistanceMatrix
    _HAS_SKBIO = True
except ImportError:
    _HAS_SKBIO = False


PALETTE = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#0072B2", "#D55E00", "#009E73", "#FF0000", "#1F1F1F",
]


def _save(fig_or_plt, save_path, dpi=300):
    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        obj = fig_or_plt if hasattr(fig_or_plt, "savefig") else plt
        obj.savefig(save_path, dpi=dpi, bbox_inches="tight")


# ---------------------------------------------------------------------------
# QC & Batch Visualization
# ---------------------------------------------------------------------------

def plot_knee(series, title="Knee Plot", threshold=None, save_path=None):
    sorted_vals = np.sort(series.dropna())
    plt.figure(figsize=(8, 5))
    plt.plot(sorted_vals, color="black", linewidth=2)
    plt.yscale("log")
    plt.ylabel(series.name or "Value")
    plt.xlabel("Samples (sorted)")
    plt.title(title)
    if threshold is not None:
        plt.axhline(threshold, color="red", linestyle="--", linewidth=1.8,
                    label=f"Threshold = {threshold}")
        plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    _save(plt, save_path)
    plt.show()


def plot_batch_metric_violins(df_obs, metrics, batch_col="Batch_ID", author_col="Author",
                               author_change_indices=None, batch_type_dict=None,
                               save_path=None):
    for col in metrics:
        if col not in df_obs.columns:
            continue
        plt.figure(figsize=(22, 5))
        ax = sns.violinplot(data=df_obs, x=batch_col, y=col,
                             hue=author_col, dodge=False, palette=PALETTE, inner=None, alpha=0.7)
        sns.boxplot(data=df_obs, x=batch_col, y=col,
                    width=0.2, color="white", ax=ax, showfliers=False, zorder=3)

        if author_change_indices is not None:
            for x in author_change_indices:
                ax.axvline(x=x, color="black", linestyle="--", alpha=0.3)

        if batch_type_dict:
            ax.set_xticklabels(
                [batch_type_dict.get(t.get_text(), t.get_text()) for t in ax.get_xticklabels()],
                rotation=45, ha="right",
            )
        ax.legend(title="Study", bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
        plt.title(f"Distribution of {col} (grouped by study)")
        plt.grid(axis="y", alpha=0.2)
        plt.tight_layout()
        if save_path:
            stem, ext = os.path.splitext(save_path)
            ext = ext or ".png"
            safe_col = col.replace(" ", "_").replace("/", "_")
            _save(plt, f"{stem}_{safe_col}{ext}")
        plt.show()


def analyze_batch_statistics(df_obs, bias_cols, batch_col):
    print("\n### Batch Effect Statistics ###")
    kw_res = []
    for col in bias_cols:
        if col not in df_obs.columns:
            continue
        groups = [g[col].dropna().values for _, g in df_obs.groupby(batch_col)]
        h_stat, p_val = stats.kruskal(*groups)
        eta_sq = max(0.0, (h_stat - len(groups) + 1) / (len(df_obs) - len(groups)))
        kw_res.append({"Metric": col, "H-stat": h_stat, "p-value": p_val, "Eta_Sq": eta_sq})

    df_kw = pd.DataFrame(kw_res).sort_values("Eta_Sq", ascending=False)
    print(df_kw.round(4))

    permanova_res = None
    if _HAS_SKBIO:
        data_z = df_obs[bias_cols].apply(stats.zscore).dropna()
        dm = DistanceMatrix(squareform(pdist(data_z, metric="euclidean")), ids=data_z.index)
        permanova_res = permanova(dm, grouping=df_obs.loc[data_z.index, batch_col], permutations=999)
        print(f"\n[PERMANOVA] Pseudo-F: {permanova_res['test statistic']:.4f}, "
              f"p={permanova_res['p-value']:.4f}")
    else:
        print("\n[Skip] scikit-bio not installed; PERMANOVA skipped.")

    return df_kw, permanova_res


def plot_bias_correlations(df_obs, bias_cols, group_col, save_path=None):
    n = df_obs[group_col].nunique()
    palette = sns.color_palette("turbo", n_colors=n)
    dark_rc = {
        "axes.facecolor": "#0D0D0D", "figure.facecolor": "#0D0D0D",
        "text.color": "white", "axes.labelcolor": "white",
        "xtick.color": "#B0B0B0", "ytick.color": "#B0B0B0",
        "axes.edgecolor": "#333333", "grid.color": "#1A1A1A",
        "legend.facecolor": "#151515", "legend.edgecolor": "#333333",
    }
    with plt.rc_context(dark_rc):
        g = sns.pairplot(
            df_obs[bias_cols + [group_col]],
            hue=group_col, diag_kind="kde",
            plot_kws={"alpha": 0.5, "s": 15, "edgecolor": "none"},
            palette=palette,
        )
        g.fig.suptitle(f"Correlations by {group_col}", y=1.02, color="white", fontsize=16)
        if g._legend:
            plt.setp(g._legend.get_texts(), color="white", fontsize=8)
            plt.setp(g._legend.get_title(), color="white", fontsize=10)
        _save(g.figure, save_path, dpi=150)
        plt.show()


# ---------------------------------------------------------------------------
# HC PCA & Variance
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Classification / Batch Power
# ---------------------------------------------------------------------------

def _build_classifiers():
    return {
        "LogReg": make_pipeline(StandardScaler(),
                                LogisticRegression(max_iter=1000, random_state=42)),
        "SVM":    make_pipeline(StandardScaler(),
                                SVC(probability=True, kernel="rbf", random_state=42)),
        "RF":     RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42),
        "GBM":    GradientBoostingClassifier(n_estimators=50, max_depth=3, random_state=42),
    }


def check_batch_classification_power(adata, bias_list, target_col="Batch_ID",
                                      hc_label="Healthy Control",
                                      phenotype_col="Phenotype_Processed",
                                      min_class_samples=10, n_repeats=20):
    df_hc = adata.obs[adata.obs[phenotype_col] == hc_label].copy()
    counts = df_hc[target_col].value_counts()
    keep = counts[counts >= min_class_samples].index.tolist()
    df_hc = df_hc[df_hc[target_col].isin(keep)]

    df_sub = df_hc[bias_list + [target_col]].dropna()
    if df_sub[target_col].nunique() < 2:
        print(f"[Skip] '{target_col}' has < 2 classes in HC data.")
        return None

    X = df_sub[bias_list].values
    le = LabelEncoder()
    y = le.fit_transform(df_sub[target_col])
    n_classes = len(le.classes_)
    models = _build_classifiers()
    sss = StratifiedShuffleSplit(n_splits=n_repeats, test_size=0.3, random_state=42)
    records = []

    for name, model in models.items():
        for i, (tr, te) in enumerate(sss.split(X, y)):
            if len(np.unique(y[te])) < 2:
                continue
            try:
                model.fit(X[tr], y[tr])
                probs = model.predict_proba(X[te])
                score = (roc_auc_score(y[te], probs[:, 1]) if n_classes == 2
                         else roc_auc_score(y[te], probs, multi_class="ovr", average="macro"))
                records.append({"Model": name, "AUC": score, "Iteration": i})
            except Exception:
                pass

    return pd.DataFrame(records)


def check_bias_power_dist(df_meta, bias_list, condition_col,
                           control_label="Healthy Control", n_repeats=20):
    df_sub = df_meta[bias_list + [condition_col]].dropna().copy()
    if len(df_sub) < 15:
        print("[Skip] Too few samples (N < 15).")
        return None

    X = df_sub[bias_list].values
    y_raw = df_sub[condition_col].values
    le = LabelEncoder()
    y = le.fit_transform(y_raw)
    n_classes = len(le.classes_)

    if np.unique(y, return_counts=True)[1].min() < 2:
        print("[Skip] Some class has < 2 samples.")
        return None

    models = _build_classifiers()
    sss = StratifiedShuffleSplit(n_splits=n_repeats, test_size=0.3, random_state=42)
    records = []

    for name, model in models.items():
        for i, (tr, te) in enumerate(sss.split(X, y)):
            if n_classes > 2 and len(np.unique(y[te])) < n_classes:
                continue
            if len(np.unique(y[te])) < 2:
                continue
            try:
                model.fit(X[tr], y[tr])
                probs = model.predict_proba(X[te])
                score = (roc_auc_score(y[te], probs[:, 1]) if n_classes == 2
                         else roc_auc_score(y[te], probs, multi_class="ovr", average="macro"))
                records.append({"Model": name, "AUC": score, "Iteration": i})
            except Exception:
                pass

    return pd.DataFrame(records)


def plot_model_auc_grid(final_df, author_col="Author", n_cols=6, save_path=None):
    authors = final_df[author_col].unique()
    n_rows = math.ceil(len(authors) / n_cols)
    my_pal = {"LogReg": "#A8D8EA", "SVM": "#AA96DA", "RF": "#FCBAD3", "GBM": "#FFFFD2"}

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3 * n_cols, 4 * n_rows), sharey=True)
    axes_flat = axes.flatten()

    for i, author in enumerate(authors):
        ax = axes_flat[i]
        sub = final_df[final_df[author_col] == author]
        sns.boxplot(data=sub, x="Model", y="AUC", hue="Model", palette=my_pal,
                    ax=ax, showfliers=False, width=0.6, legend=False)
        sns.stripplot(data=sub, x="Model", y="AUC", color="black",
                      alpha=0.3, jitter=True, size=3, ax=ax)
        ax.axhline(0.5, color="gray", linestyle="--", alpha=0.6)
        ax.axhline(0.7, color="red", linestyle=":", linewidth=1.5)
        ax.set_title(author)
        ax.set_ylim(0, 1.1)
        ax.set_xlabel("")
        ax.grid(alpha=0.2)
        if i % n_cols != 0:
            ax.set_ylabel("")

    for j in range(i + 1, len(axes_flat)):
        fig.delaxes(axes_flat[j])

    plt.suptitle("Bias ~ Phenotype AUC Scores (Binary)", fontsize=20)
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()


# ---------------------------------------------------------------------------
# Per-Study PCA Diagnostics (called from DataAnalysisPipeline)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Partial RDA Visualization (called from DataAnalysisPipeline)
# ---------------------------------------------------------------------------

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
    part_cols   = ["pheno_unique", "conf_unique", "shared", "unexplained"]
    part_labels = ["Phenotype Unique", "Confounder Unique", "Shared", "Unexplained"]
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

    part_cols   = ["pheno_unique", "conf_unique", "shared", "unexplained"]
    part_labels = ["Phenotype Unique", "Confounder Unique", "Shared", "Unexplained"]
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
    part_labels = ["Confounder Unique (sum)", "Shared", "Unexplained"]
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


# ---------------------------------------------------------------------------
# Gene-wise Bias RDA Visualization
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Discrete Covariate Batch-Effect Assessment
# ---------------------------------------------------------------------------

def check_discrete_covariate_batch_effects(
    adata,
    bias_list,
    covariates,
    hc_label="Healthy Control",
    phenotype_col="Phenotype_Processed",
    min_class_samples=10,
    n_repeats=20,
):
    """
    HC 샘플에서 bias_list를 feature로 사용해 각 discrete covariate를 분류.
    Returns (df_auc, df_meta).
      df_auc : Covariate / Model / AUC / Iteration
      df_meta: covariate / n_classes / n_samples / classes
    """
    obs_hc     = adata.obs[adata.obs[phenotype_col] == hc_label].copy()
    avail_bias = [b for b in bias_list if b in obs_hc.columns]
    all_records, meta_records = [], []
    models = _build_classifiers()

    for covar in covariates:
        if covar not in obs_hc.columns:
            print(f"[Skip] {covar}: column not found")
            continue
        df_sub = obs_hc.dropna(subset=avail_bias + [covar]).copy()
        df_sub[covar] = df_sub[covar].astype(str)
        counts = df_sub[covar].value_counts()
        keep   = counts[counts >= min_class_samples].index.tolist()
        if len(keep) < 2:
            print(f"[Skip] {covar}: < 2 classes with ≥{min_class_samples} samples "
                  f"({dict(counts.head(5))})")
            continue
        df_sub = df_sub[df_sub[covar].isin(keep)]
        X  = df_sub[avail_bias].values
        le = LabelEncoder()
        y  = le.fit_transform(df_sub[covar])
        n_cls = len(le.classes_)
        meta_records.append({
            "covariate": covar, "n_classes": n_cls,
            "n_samples": len(df_sub),
            "classes": ", ".join(str(c) for c in le.classes_),
        })
        sss = StratifiedShuffleSplit(n_splits=n_repeats, test_size=0.3, random_state=42)
        for name, model in models.items():
            for i, (tr, te) in enumerate(sss.split(X, y)):
                if len(np.unique(y[te])) < 2:
                    continue
                try:
                    model.fit(X[tr], y[tr])
                    probs = model.predict_proba(X[te])
                    score = (
                        roc_auc_score(y[te], probs[:, 1]) if n_cls == 2
                        else roc_auc_score(y[te], probs, multi_class="ovr", average="macro")
                    )
                    all_records.append({"Covariate": covar, "Model": name,
                                        "AUC": score, "Iteration": i})
                except Exception:
                    pass
        mean_auc = np.mean([r["AUC"] for r in all_records if r["Covariate"] == covar])
        print(f"  {covar:45s}: n_classes={n_cls:2d}  n={len(df_sub):4d}  mean_AUC={mean_auc:.3f}")

    return pd.DataFrame(all_records), pd.DataFrame(meta_records)


def plot_covariate_auc_heatmap(df_auc, df_meta=None, save_path=None):
    """Mean AUC heatmap (discrete covariates × classifiers)."""
    label_map = {
        "instrument": "Sequencer Model",
        "rna_extraction_kit_short_name": "RNA Extraction Kit",
        "plasma_tubes_short_name": "Plasma Tube Type",
        "library_prep_kit_short_name": "Library Prep Kit",
        "Centrifuge_Protocol": "Centrifuge Protocol",
        "broad_protocol_category": "Broad Protocol Category",
        "cdna_library_type": "cDNA Library Type",
        "dnase": "DNase Treatment",
        "UMI": "UMI",
        "librarylayout": "Library Layout",
        "library_selection": "Library Selection",
    }
    pivot = (df_auc.groupby(["Covariate", "Model"])["AUC"]
             .mean().unstack().fillna(np.nan))
    pivot.index = [label_map.get(c, c) for c in pivot.index]
    pivot["_best"] = pivot.max(axis=1)
    pivot = pivot.sort_values("_best", ascending=False).drop(columns="_best")

    fig, ax = plt.subplots(figsize=(8, max(4, len(pivot) * 0.55 + 2)))
    sns.heatmap(
        pivot, annot=True, fmt=".3f",
        cmap="RdYlGn", vmin=0.5, vmax=1.0, center=0.75,
        linewidths=0.5, ax=ax,
        cbar_kws={"label": "Mean AUC (macro-OVR)"},
    )
    ax.set_title(
        "Discrete Covariate Separability via Bias Metrics  (HC only)\n"
        "AUC > 0.7 indicates meaningful batch signal",
        fontweight="bold", pad=14,
    )
    ax.set_xlabel("Classifier")
    ax.set_ylabel("")
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()
    return pivot


def plot_covariate_auc_violins(df_auc, n_cols=4, save_path=None):
    """Box+strip grid: one subplot per discrete covariate, models on x-axis."""
    label_map = {
        "instrument": "Sequencer Model",
        "rna_extraction_kit_short_name": "RNA Extraction Kit",
        "plasma_tubes_short_name": "Plasma Tube Type",
        "library_prep_kit_short_name": "Library Prep Kit",
        "Centrifuge_Protocol": "Centrifuge Protocol",
        "broad_protocol_category": "Broad Protocol Category",
        "cdna_library_type": "cDNA Library Type",
        "dnase": "DNase Treatment",
        "UMI": "UMI",
        "librarylayout": "Library Layout",
        "library_selection": "Library Selection",
    }
    covariates = df_auc["Covariate"].unique()
    n_rows = math.ceil(len(covariates) / n_cols)
    my_pal = {"LogReg": "#A8D8EA", "SVM": "#AA96DA", "RF": "#FCBAD3", "GBM": "#FFFFD2"}

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4 * n_cols, 4 * n_rows), sharey=True)
    axes_flat = axes.flatten() if n_rows * n_cols > 1 else [axes]

    for i, covar in enumerate(covariates):
        ax   = axes_flat[i]
        sub  = df_auc[df_auc["Covariate"] == covar]
        sns.boxplot(data=sub, x="Model", y="AUC", hue="Model",
                    palette=my_pal, ax=ax, showfliers=False, width=0.6, legend=False)
        sns.stripplot(data=sub, x="Model", y="AUC",
                      color="black", alpha=0.3, jitter=True, size=3, ax=ax)
        ax.axhline(0.5, color="gray", linestyle="--", alpha=0.6, linewidth=1)
        ax.axhline(0.7, color="red",  linestyle=":", linewidth=1.5)
        ax.set_title(label_map.get(covar, covar), fontweight="bold", fontsize=10)
        ax.set_ylim(0, 1.1)
        ax.set_xlabel("")
        ax.grid(alpha=0.2)
        if i % n_cols != 0:
            ax.set_ylabel("")

    for j in range(len(covariates), len(axes_flat)):
        fig.delaxes(axes_flat[j])

    plt.suptitle(
        "Discrete Covariates ~ Bias Metrics AUC  (HC only, 20-repeat shuffle CV)\n"
        "Red dotted line = AUC 0.7 concern threshold",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    _save(fig, save_path)
    plt.show()


def plot_gene_wise_bias_summary(df_detail, group_name, save_path=None):
    """Joint R² histogram + contamination severity bar chart for one group."""
    joint_r2 = df_detail["Joint_R2_All_Biases"]
    bins_edges = [-np.inf, 0.05, 0.10, 0.30, np.inf]
    sev_labels = ["Minimal (< 5%)", "Moderate (5–10%)", "High (10–30%)", "Severe (> 30%)"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)

    sns.histplot(joint_r2, bins=75, color="steelblue", ax=ax1, edgecolor="black")
    ax1.axvline(0.10, color="red", linestyle="--", linewidth=1.5, alpha=0.8)
    ax1.set_xlim(0, 1.0)
    ax1.set_xlabel("Total Variance Explained by All Biases (Joint R²)")
    ax1.set_ylabel("Number of Genes")
    ax1.set_title(group_name)
    ax1.grid(axis="y", linestyle="--", alpha=0.3)

    cats = pd.cut(joint_r2, bins=bins_edges, labels=sev_labels)
    counts = cats.value_counts(sort=False)
    n_total = len(joint_r2)
    pcts = counts / n_total * 100
    df_sev = pd.DataFrame({"Severity": sev_labels, "Pct": pcts.values, "N": counts.values})
    df_sev = df_sev.iloc[::-1].reset_index(drop=True)

    sns.barplot(data=df_sev, x="Pct", y="Severity", ax=ax2,
                palette="Reds_r", hue="Severity", legend=False, edgecolor="black")
    max_pct = df_sev["Pct"].max()
    for idx, row in df_sev.iterrows():
        ax2.text(row["Pct"] + max_pct * 0.02, idx,
                 f"{row['Pct']:.1f}% ({int(row['N']):,})", va="center", fontsize=11, color="#333333")
    ax2.set_xlim(0, max_pct * 1.35)
    ax2.set_xlabel("Proportion of Total Genes (%)")
    ax2.set_ylabel("Contamination Severity (Joint R²)")
    ax2.grid(axis="x", linestyle="--", alpha=0.3)

    _save(fig, save_path)
    plt.show()

