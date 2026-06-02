import math

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
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
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
            plt.savefig(save_path, dpi=200, bbox_inches="tight")
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
        g.fig.suptitle(f"Correlations by {group_col}", y=1.02, color="white",
                       fontsize=16)
        if g._legend:
            plt.setp(g._legend.get_texts(), color="white", fontsize=8)
            plt.setp(g._legend.get_title(), color="white", fontsize=10)
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.show()


def plot_hc_variance_and_heatmap(adata, bias_metrics, hc_label="Healthy Control",
                                  phenotype_col="Phenotype_Processed", n_pcs=50):
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
    plt.show()


def plot_hc_pca_grid(adata, plot_vars, hc_label="Healthy Control",
                      phenotype_col="Phenotype_Processed"):
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
    plt.show()


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


def plot_model_auc_grid(final_df, author_col="Author", n_cols=6):
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
    plt.show()


def plot_pca_scree_and_bias(var_ratios, study, obs_df, active_metrics):
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
    plt.show()


def plot_pca_scatter_grid(adata_pca, plot_keys, var_ratios, key_title_map=None):
    """PCA scatter grid colored by each metadata key."""
    pc1_var, pc2_var = var_ratios[0], var_ratios[1]
    n_cols = 3
    n_rows = math.ceil(len(plot_keys) / n_cols)
    fig2, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows))
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
    plt.show()


def plot_rda_unique_heatmap(df_unique, use_hvg):
    """Heatmap of per-variable unique R² contributions across studies."""
    fig_w = max(10, len(df_unique.columns) * .8)
    fig_h = max(6, len(df_unique) * 0.4 + 2)
    plt.figure(figsize=(fig_w, fig_h))
    ax = sns.heatmap(df_unique, annot=True, fmt=".3f", cmap="Reds", annot_kws={"size": 11},
                     linewidths=0.5, mask=df_unique.isna())
    sns.heatmap(df_unique.isna(), cmap=["white", "lightgrey"],
                cbar=False, ax=ax, mask=~df_unique.isna(), linewidths=0.5)
    plt.title(f"Per-Variable Unique Contribution (Partial RDA)\n(HVG: {use_hvg})", pad=15)
    plt.ylabel("Studies (Authors)")
    plt.xlabel("Covariates")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()


def plot_rda_variance_partition(df_partition, title, x_label="Studies (Authors)", figsize=None):
    """Stacked horizontal bar chart of variance partition."""
    part_cols   = ["pheno_unique", "conf_unique", "shared", "unexplained"]
    part_labels = ["Phenotype Unique", "Confounder Unique", "Shared", "Unexplained"]
    part_colors = ["steelblue", "coral", "mediumseagreen", "lightgrey"]
    df_plot = df_partition[part_cols].fillna(0)
    fig_h = figsize[1] if figsize else max(5, len(df_plot) * 0.4 + 2)
    fig_w = figsize[0] if figsize else 10
    plt.figure(figsize=(fig_w, fig_h))
    bottom = np.zeros(len(df_plot))
    for col, label, color in zip(part_cols, part_labels, part_colors):
        vals = df_plot[col].values
        plt.barh(df_plot.index, vals, left=bottom, color=color, label=label, edgecolor="black")
        bottom += vals
    plt.xlabel("Proportion of Total Variance")
    plt.ylabel(x_label)
    plt.title(title, pad=15)
    plt.xlim(0, 1)

    plt.legend(loc="center left", frameon=False, bbox_to_anchor=(1.01, 0.5))
    plt.tight_layout()
    plt.show()


def plot_cascade_spaghetti(df_metrics, palette=None):
    """Spaghetti plot of phenotype R² trajectory across sequential confounders."""
    if df_metrics.empty:
        return
    palette = palette or PALETTE
    all_steps = df_metrics.columns.tolist()
    x = np.arange(len(all_steps))

    plt.figure(figsize=(max(10, len(all_steps) * 0.4), 6))
    for i, study in enumerate(df_metrics.index):
        color = palette[i % len(palette)]
        plt.plot(x, df_metrics.loc[study], marker="o", markersize=6,
                 linewidth=2, alpha=0.85, color=color, label=study)

    plt.xticks(x, all_steps, rotation=90, ha="right")
    plt.ylabel("Adjusted Partial R²")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left", title='Studies', frameon=False)
    plt.tight_layout()
    plt.show()


def plot_normalization_partition(df_partition_all, studies):
    """Per-study variance partition bars across normalization layers."""
    n_cols = min(5, len(studies))
    n_rows = math.ceil(len(studies) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4.5, n_rows * 4.5), sharey=True)
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
            ax.bar(x, vals, bottom=bottom, color=color, width=0.6, edgecolor="white", linewidth=0.5)
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
    plt.show()


def plot_hc_rda_results(sr_unique, r2_all, batch_col, layer, unique_dict):
    plt.figure(figsize=(7, max(4, len(sr_unique) * 0.45 + 1.5)))
    plt.barh(sr_unique.index, sr_unique.values, color="steelblue", edgecolor="black")
    plt.xlabel("Unique R² (Partial RDA)")
    plt.axvline(0, color="black", linewidth=0.8)
    plt.grid(axis="x", linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.show()

    conf_unique_sum = sum(v for v in unique_dict.values() if not np.isnan(v))
    shared = max(0.0, r2_all - conf_unique_sum)
    unexplained = max(0.0, 1.0 - r2_all)
    part_vals = [conf_unique_sum, shared, unexplained]
    part_labels = ["Confounder Unique (sum)", "Shared", "Unexplained"]
    part_colors = ["coral", "mediumseagreen", "lightgrey"]
    plt.figure(figsize=(6.5, 1.))
    bottom = 0.0
    for val, label, color in zip(part_vals, part_labels, part_colors):
        plt.barh(["HC"], [val], left=[bottom], color=color, label=label, height=0.3, edgecolor="black")
        bottom += val
    plt.xlabel("Proportion of Total Variance")
    plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.6), ncol=3, frameon=False)   
    plt.show()