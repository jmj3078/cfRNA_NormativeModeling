import math
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupShuffleSplit, StratifiedShuffleSplit
from sklearn.preprocessing import LabelEncoder

from analysis.plot_utils import PALETTE, _build_classifiers, _save

COVARIATE_LABELS = {
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


try:
    from skbio import DistanceMatrix
    from skbio.stats.distance import permanova
    _HAS_SKBIO = True
except ImportError:
    _HAS_SKBIO = False


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


def check_discrete_covariate_batch_effects(
    adata,
    bias_list,
    covariates,
    hc_label="Healthy Control",
    phenotype_col="Phenotype_Processed",
    min_class_samples=10,
    n_repeats=20,
    group_col="Batch_ID",
    min_batches_per_class=2,
):
    # Protocol labels are constant within a batch, so a random split lets the same
    # batch appear in train and test: the classifier can memorise batch clusters.
    # Splits are grouped by batch, and classes confined to a single batch are dropped
    # because for them the covariate and the batch are the same variable.
    obs_hc = adata.obs[adata.obs[phenotype_col] == hc_label].copy()
    avail_bias = [b for b in bias_list if b in obs_hc.columns]
    all_records, meta_records = [], []
    models = _build_classifiers()

    for covar in covariates:
        if covar not in obs_hc.columns:
            print(f"[Skip] {covar}: column not found")
            continue
        df_sub = obs_hc.dropna(subset=avail_bias + [covar, group_col]).copy()
        df_sub[covar] = df_sub[covar].astype(str)
        counts = df_sub[covar].value_counts()
        n_batches = df_sub.groupby(covar, observed=True)[group_col].nunique()
        keep = [c for c in counts.index
                if counts[c] >= min_class_samples and n_batches[c] >= min_batches_per_class]
        dropped = [c for c in counts.index if c not in keep]
        if len(keep) < 2:
            print(f"[Skip] {covar}: < 2 classes with ≥{min_class_samples} samples spanning "
                  f"≥{min_batches_per_class} batches (batches/class: {dict(n_batches)})")
            continue
        if dropped:
            print(f"  [Drop] {covar}: classes confined to <{min_batches_per_class} batches "
                  f"or too few samples: {dropped}")
        df_sub = df_sub[df_sub[covar].isin(keep)]
        X = df_sub[avail_bias].values
        le = LabelEncoder()
        y = le.fit_transform(df_sub[covar])
        groups = df_sub[group_col].values
        n_cls = len(le.classes_)

        gss = GroupShuffleSplit(n_splits=n_repeats * 10, test_size=0.3, random_state=42)
        folds = [(tr, te) for tr, te in gss.split(X, y, groups)
                 if np.unique(y[tr]).size == n_cls and np.unique(y[te]).size == n_cls][:n_repeats]
        meta_records.append({
            "covariate": covar, "n_classes": n_cls,
            "n_samples": len(df_sub),
            "n_batches": df_sub[group_col].nunique(),
            "n_folds": len(folds),
            "classes": ", ".join(str(c) for c in le.classes_),
        })
        if not folds:
            print(f"  [Skip] {covar}: no batch-disjoint split keeps all {n_cls} classes "
                  f"on both sides")
            continue

        for name, model in models.items():
            for i, (tr, te) in enumerate(folds):
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
        print(f"  {covar:45s}: n_classes={n_cls:2d}  n={len(df_sub):4d}  "
              f"n_batches={df_sub[group_col].nunique():2d}  folds={len(folds):2d}  "
              f"mean_AUC={mean_auc:.3f}")

    return pd.DataFrame(all_records), pd.DataFrame(meta_records)


def plot_covariate_auc_heatmap(df_auc, save_path=None):
    """Mean AUC heatmap (discrete covariates × classifiers)."""
    pivot = (df_auc.groupby(["Covariate", "Model"])["AUC"]
             .mean().unstack().fillna(np.nan))
    pivot.index = [COVARIATE_LABELS.get(c, c) for c in pivot.index]
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
    covariates = df_auc["Covariate"].unique()
    n_rows = math.ceil(len(covariates) / n_cols)
    my_pal = {"LogReg": "#A8D8EA", "SVM": "#AA96DA", "RF": "#FCBAD3", "GBM": "#FFFFD2"}

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4 * n_cols, 4 * n_rows), sharey=True)
    axes_flat = axes.flatten() if n_rows * n_cols > 1 else [axes]

    for i, covar in enumerate(covariates):
        ax = axes_flat[i]
        sub = df_auc[df_auc["Covariate"] == covar]
        sns.boxplot(data=sub, x="Model", y="AUC", hue="Model",
                    palette=my_pal, ax=ax, showfliers=False, width=0.6, legend=False)
        sns.stripplot(data=sub, x="Model", y="AUC",
                      color="black", alpha=0.3, jitter=True, size=3, ax=ax)
        ax.axhline(0.5, color="gray", linestyle="--", alpha=0.6, linewidth=1)
        ax.axhline(0.7, color="red", linestyle=":", linewidth=1.5)
        ax.set_title(COVARIATE_LABELS.get(covar, covar), fontweight="bold", fontsize=10)
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
