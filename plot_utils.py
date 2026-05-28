import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import scanpy as sc
import statsmodels.formula.api as smf
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import make_pipeline
from statsmodels.stats.multitest import multipletests
from tqdm.auto import tqdm

try:
    from skbio.stats.distance import permanova
    from skbio import DistanceMatrix
    _HAS_SKBIO = True
except ImportError:
    _HAS_SKBIO = False


# ---------------------------------------------------------------------------
# Batch / QC visualization
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
                             hue=author_col, dodge=False, palette="Set2", inner=None)
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
        ax.legend(title="Study", bbox_to_anchor=(1.01, 1), loc="upper left")
        plt.title(f"Distribution of {col} (grouped by study)", fontweight="bold")
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
                       fontsize=16, fontweight="bold")
        if g._legend:
            plt.setp(g._legend.get_texts(), color="white", fontsize=8)
            plt.setp(g._legend.get_title(), color="white", fontsize=10, fontweight="bold")
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.show()


# ---------------------------------------------------------------------------
# HC-level PCA and impact analysis
# ---------------------------------------------------------------------------

def plot_hc_impact_bar(adata, bias_metrics, hc_label="Healthy Control",
                        phenotype_col="Phenotype_Processed", n_pcs=20):
    adata_hc = adata[adata.obs[phenotype_col] == hc_label]
    pcs = adata_hc.obsm["X_pca"][:, :n_pcs]
    v_ratios = adata_hc.uns["pca"]["variance_ratio"][:n_pcs]

    results = []
    for var in bias_metrics:
        if var not in adata_hc.obs.columns:
            continue
        valid = adata_hc.obs[var].notna()
        if valid.sum() < 2:
            continue
        curr_obs = adata_hc.obs[valid][var]
        is_cat = not pd.api.types.is_numeric_dtype(curr_obs)
        impact = 0.0
        for i in range(pcs.shape[1]):
            try:
                df_tmp = pd.DataFrame({"y": pcs[valid, i], "x": curr_obs})
                formula = "y ~ C(x)" if is_cat else "y ~ x"
                impact += smf.ols(formula, data=df_tmp).fit().rsquared * v_ratios[i]
            except Exception:
                pass
        results.append({"Metric": var, "Impact": impact})

    df_res = pd.DataFrame(results).sort_values("Impact", ascending=False)
    plt.figure(figsize=(8, 7))
    sns.barplot(data=df_res, x="Metric", y="Impact", palette="Reds", edgecolor="0.3")
    plt.title("Global Covariate Impact on HC Data (Σ R² × VarRatio)", fontweight="bold")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Weighted R² Impact Score")
    for i, v in enumerate(df_res["Impact"]):
        plt.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontweight="bold")
    plt.tight_layout()
    plt.show()
    return df_res


def plot_hc_variance_and_heatmap(adata, bias_metrics, hc_label="Healthy Control",
                                  phenotype_col="Phenotype_Processed", n_pcs=30):
    adata_hc = adata[adata.obs[phenotype_col] == hc_label]
    v_ratios = adata_hc.uns["pca"]["variance_ratio"][:n_pcs]

    fig, ax = plt.subplots(1, 2, figsize=(30, 6), gridspec_kw={"width_ratios": [1, 2.5]})
    ax[0].plot(range(1, len(v_ratios) + 1), v_ratios, "o-k", alpha=0.8, markersize=5)
    ax[0].set_title("Scree Plot", fontweight="bold")
    ax[0].set_ylabel("Explained Variance Ratio")
    ax[0].set_xlabel("Principal Component")
    ax[0].grid(True, linestyle="--", alpha=0.4)

    num_metrics = [m for m in bias_metrics if m in adata_hc.obs.columns
                   and pd.api.types.is_numeric_dtype(adata_hc.obs[m])]
    df_m = adata_hc.obs[num_metrics].copy()
    df_scaled = (df_m - df_m.mean()) / (df_m.std() + 1e-9)
    sns.heatmap(df_scaled.T, cmap="RdBu_r", center=0, ax=ax[1],
                cbar_kws={"label": "Z-score"}, xticklabels=False)
    ax[1].set_title("Metric Distribution across HC Samples (Z-score)", fontweight="bold")
    ax[1].set_ylabel("Metrics")
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
        axes_flat[i].set_title(var, fontweight="bold")
        axes_flat[i].set_xlabel(f"PC1 ({pc1_v:.1%})")
        axes_flat[i].set_ylabel(f"PC2 ({pc2_v:.1%})")

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis("off")
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Normalization comparison
# ---------------------------------------------------------------------------

def analyze_normalization_efficiency(adata, fixed_metrics, bias_metrics, n_pcs=50):
    metrics = ([m for m in fixed_metrics if m in adata.obs.columns]
               + [m for m in bias_metrics if m in adata.obs.columns])
    results = []

    for author in adata.obs["Author"].unique():
        print(f"  Processing: {author}...")
        adata_sub = adata[adata.obs["Author"] == author].copy()

        for layer in adata_sub.layers.keys():
            adata_sub.X = adata_sub.layers[layer].copy()
            sc.pp.filter_genes(adata_sub, min_cells=1)
            if adata_sub.shape[1] < n_pcs:
                continue
            try:
                sc.tl.pca(adata_sub, n_comps=n_pcs)
                pcs = adata_sub.obsm["X_pca"]
                v_ratios = adata_sub.uns["pca"]["variance_ratio"]

                for metric in metrics:
                    if metric not in adata_sub.obs.columns:
                        continue
                    is_cat = not pd.api.types.is_numeric_dtype(adata_sub.obs[metric])
                    total_r2 = 0.0
                    for i in range(min(n_pcs, pcs.shape[1])):
                        try:
                            df_tmp = pd.DataFrame({"y": pcs[:, i], "x": adata_sub.obs[metric]})
                            formula = "y ~ C(x)" if is_cat else "y ~ x"
                            total_r2 += smf.ols(formula, data=df_tmp).fit().rsquared * v_ratios[i]
                        except Exception:
                            pass
                    results.append({"Study": author, "Method": layer,
                                    "Metric": metric, "Impact": total_r2})
            except Exception as e:
                print(f"    [Skip] {author} - {layer}: {e}")

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Classification utilities
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
        ax.set_title(author, fontweight="bold")
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


# ---------------------------------------------------------------------------
# GP model evaluation / marker discovery
# ---------------------------------------------------------------------------

def evaluate_hc_reconstruction(pipeline, X_hc, Y_mask_hc, Y_target_hc, valid_genes,
                                 gp_model_cls, device):
    import gpytorch
    import torch

    X_hc = pipeline._ensure_array(X_hc, dtype=np.float32)
    Y_mask_hc = pipeline._ensure_array(Y_mask_hc, dtype=bool)
    Y_target_hc = pipeline._ensure_array(Y_target_hc, dtype=np.float32)

    valid_set = set(valid_genes)
    n_valid = len(valid_genes)
    g2idx = {g: i for i, g in enumerate(valid_genes)}
    valid_col_ids = [np.where(pipeline.gene_names == g)[0][0] for g in valid_genes]

    Y_mask_v = Y_mask_hc[:, valid_col_ids]
    Y_tgt_v = Y_target_hc[:, valid_col_ids]

    emp_mean = np.array([Y_tgt_v[Y_mask_v[:, i], i].mean() if Y_mask_v[:, i].any() else 0
                         for i in range(n_valid)])
    emp_std  = np.array([Y_tgt_v[Y_mask_v[:, i], i].std()  if Y_mask_v[:, i].any() else 0
                         for i in range(n_valid)])
    auc_scores = np.full(n_valid, np.nan)

    pipeline.lr_model.to(device)
    pipeline.lr_model.eval()
    X_t = torch.tensor(X_hc, dtype=torch.float32, device=device)
    with torch.no_grad():
        probs = torch.sigmoid(pipeline.lr_model(X_t)).cpu().numpy()
    pipeline.lr_model.cpu()

    for i, g in enumerate(valid_genes):
        if g in pipeline.lr_col_map:
            lr_idx = pipeline.lr_col_map[g]
            y_true = Y_mask_v[:, i].astype(int)
            if len(np.unique(y_true)) > 1:
                auc_scores[i] = roc_auc_score(y_true, probs[:, lr_idx])

    pred_mean = np.zeros(n_valid)
    pred_std  = np.zeros(n_valid)

    for n_pos, chunks in tqdm(pipeline.batched_gp_models.items(), desc="Evaluating GP"):
        for chunk in chunks:
            target = [g for g in chunk["gene_names"] if g in valid_set]
            if not target:
                continue
            b_size = len(chunk["gene_names"])
            b_shape = torch.Size([b_size])
            train_x = chunk["train_x"].to(device)
            train_y = chunk["train_y"].to(device)
            likelihood = gpytorch.likelihoods.GaussianLikelihood(batch_shape=b_shape).to(device)
            model = gp_model_cls(train_x, train_y, likelihood,
                                  pipeline.num_features, b_shape).to(device)
            model.load_state_dict(chunk["model_state"])
            likelihood.load_state_dict(chunk["likelihood_state"])
            model.eval(); likelihood.eval()

            with torch.no_grad(), gpytorch.settings.fast_pred_var():
                X_exp = X_t.unsqueeze(0).expand(b_size, -1, -1)
                mu_sc = likelihood(model(X_exp)).mean.cpu().numpy()

            y_m, y_s = chunk["y_mean"], chunk["y_std"]
            for local_idx, g in enumerate(chunk["gene_names"]):
                if g in g2idx:
                    vi = g2idx[g]
                    gc = np.where(pipeline.gene_names == g)[0][0]
                    pos = Y_mask_hc[:, gc]
                    if pos.any():
                        mu_r = (mu_sc[local_idx, pos] * y_s[local_idx]) + y_m[local_idx]
                        pred_mean[vi] = mu_r.mean()
                        pred_std[vi] = np.sqrt(
                            chunk["noise_var"][local_idx] * y_s[local_idx] ** 2 + mu_r.std() ** 2
                        )
            del model, likelihood, train_x, train_y
            torch.cuda.empty_cache()

    df_eval = pd.DataFrame({
        "Gene": valid_genes,
        "Dropout_AUC": auc_scores,
        "Empirical_Mean": emp_mean, "Predicted_Mean": pred_mean,
        "Empirical_Std":  emp_std,  "Predicted_Std":  pred_std,
    })

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    sns.histplot(df_eval["Dropout_AUC"].dropna(), bins=40, color="rebeccapurple",
                 ax=axes[0], kde=True)
    axes[0].axvline(0.5, color="red", linestyle="--", lw=2)
    med = df_eval["Dropout_AUC"].dropna().median()
    axes[0].set_title(f"Track A: Dropout AUC\nMedian={med:.3f}")

    for col_e, col_p, color, ax in [
        ("Empirical_Mean", "Predicted_Mean", "teal",  axes[1]),
        ("Empirical_Std",  "Predicted_Std",  "coral", axes[2]),
    ]:
        plot_df = df_eval[(df_eval[col_e] != 0) & (df_eval[col_p] != 0)].dropna()
        sns.scatterplot(data=plot_df, x=col_e, y=col_p, ax=ax,
                        alpha=0.5, color=color, s=15, edgecolor="none")
        if not plot_df.empty:
            lo = min(plot_df[col_e].min(), plot_df[col_p].min())
            hi = max(plot_df[col_e].max(), plot_df[col_p].max())
            ax.plot([lo, hi], [lo, hi], "k--", lw=1.5, alpha=0.8)
            r2 = np.corrcoef(plot_df[col_e], plot_df[col_p])[0, 1] ** 2
            ax.set_title(f"Track B: {col_e.split('_')[1]}\n$R^2$={r2:.3f}")
        ax.set_xlabel(f"Empirical {col_e.split('_')[1]}")
        ax.set_ylabel(f"Predicted {col_p.split('_')[1]}")

    plt.tight_layout()
    plt.show()
    return df_eval


def plot_normative_markers(z_results, adata_disease, phenotype_col="Phenotype_Processed",
                            fdr_th=0.01, z_abs_min=2.5, top_k=20):
    from scipy import stats as scipy_stats

    df_z = z_results["total"]
    z_vals = df_z.values
    valid = ~np.isnan(z_vals)
    raw_p = np.full(z_vals.shape, np.nan)
    raw_p[valid] = 2 * scipy_stats.norm.sf(np.abs(z_vals[valid]))

    fdr = np.full(z_vals.shape, np.nan)
    if valid.any():
        _, corr_p, _, _ = multipletests(raw_p[valid], method="fdr_bh")
        fdr[valid] = corr_p

    df_fdr = pd.DataFrame(fdr, index=df_z.index, columns=df_z.columns)
    df_hit = (df_fdr < fdr_th) & (df_z.abs() > z_abs_min)

    phenotypes = adata_disease.obs[phenotype_col].astype(str).values
    groups = [g for g in np.unique(phenotypes) if g not in ("nan", "Unknown", "None")]
    global_prev = df_hit.mean(axis=0)

    final_markers = []
    for grp in groups:
        mask = phenotypes == grp
        if mask.sum() < 2:
            continue
        grp_prev = df_hit.loc[mask].mean(axis=0)
        spec = grp_prev / (global_prev + 0.05)
        cands = grp_prev[grp_prev >= 0.2]
        if not cands.empty:
            final_markers.extend(spec[cands.index].nlargest(top_k).index)

    selected = list(dict.fromkeys(final_markers))
    if not selected:
        print("[Result] No markers found.")
        return None

    palette = sns.color_palette("tab20", len(groups))
    lut = dict(zip(groups, palette))
    col_colors = pd.Series(phenotypes, index=df_z.index).map(lut).fillna("lightgray").values

    plot_data = df_z[selected].replace([np.inf, -np.inf], np.nan).fillna(0).T
    g = sns.clustermap(
        plot_data,
        cmap="vlag", center=0, vmin=-5, vmax=5,
        col_colors=col_colors, method="ward", metric="euclidean",
        figsize=(16, 14), xticklabels=False,
        yticklabels=len(selected) < 100,
        dendrogram_ratio=(0.1, 0.03),
        cbar_pos=(1.05, 0.2, 0.2, 0.015),
        cbar_kws={"label": "Normative Z-score", "orientation": "horizontal", "ticks": [-5, 0, 5]},
    )
    legend_elements = [mlines.Line2D([0], [0], marker="s", color="w", label=lbl,
                                     markerfacecolor=c, markersize=10)
                       for lbl, c in lut.items()]
    g.ax_heatmap.legend(handles=legend_elements, title="Disease Type",
                         bbox_to_anchor=(1.1, 1.0), loc="upper left")
    plt.suptitle(f"Normative Markers (N={len(selected)})", y=1.02, fontsize=16)
    plt.show()
    return selected


def evaluate_lasso_classification(df_z_hc, df_z_dis, adata_disease, valid_genes,
                                    phenotype_col="Phenotype_Processed", min_n=30):
    X_hc = df_z_hc[valid_genes].replace([np.inf, -np.inf], np.nan).fillna(0).values
    X_dis_full = df_z_dis[valid_genes].replace([np.inf, -np.inf], np.nan).fillna(0).values
    phenotypes = adata_disease.obs[phenotype_col].astype(str).values
    diseases = [d for d in np.unique(phenotypes) if d not in ("nan", "Unknown", "None")]

    results = []
    plt.figure(figsize=(12, 8))
    sns.set_palette("tab20")

    for disease in diseases:
        mask = phenotypes == disease
        X_dis = X_dis_full[mask]
        if len(X_dis) < min_n:
            continue
        idx_hc = np.random.choice(len(X_hc), len(X_dis), replace=False)
        X_comb = np.vstack([X_hc[idx_hc], X_dis])
        y_comb = np.concatenate([np.zeros(len(X_dis)), np.ones(len(X_dis))])

        clf = LogisticRegressionCV(cv=3, penalty="l1", solver="liblinear",
                                    scoring="roc_auc", max_iter=2000, random_state=42)
        cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
        probs = np.zeros(len(y_comb))

        for tr, te in cv.split(X_comb, y_comb):
            clf.fit(X_comb[tr], y_comb[tr])
            probs[te] = clf.predict_proba(X_comb[te])[:, 1]

        fpr, tpr, _ = roc_curve(y_comb, probs)
        roc_auc_val = auc(fpr, tpr)
        results.append({"Disease": disease, "N_Disease": len(X_dis), "Lasso_AUC": roc_auc_val})
        plt.plot(fpr, tpr, lw=2, label=f"{disease:15s} (AUC={roc_auc_val:.3f}, N={len(X_dis)})")

    plt.plot([0, 1], [0, 1], color="black", lw=1, linestyle="--", alpha=0.5)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Diagnostic Performance (Lasso L1 on Z-scores)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.show()
    return pd.DataFrame(results).sort_values("Lasso_AUC", ascending=False)


def plot_lasso_marker_heatmap(df_z_hc, df_z_dis, adata_disease, valid_genes,
                               phenotype_col="Phenotype_Processed", top_k=20, min_n=30):
    phenotypes = adata_disease.obs[phenotype_col].astype(str).values
    diseases = [d for d in np.unique(phenotypes) if d not in ("nan", "Unknown", "None")]
    all_genes = set()
    lut = dict(zip(diseases, sns.color_palette("tab20", len(diseases))))

    for disease in diseases:
        mask = phenotypes == disease
        X_dis = df_z_dis.loc[mask, valid_genes].replace([np.inf, -np.inf], np.nan).fillna(0).values
        if len(X_dis) < min_n:
            continue
        idx_hc = np.random.choice(len(df_z_hc), len(X_dis), replace=False)
        X_hc_s = df_z_hc.iloc[idx_hc][valid_genes].replace([np.inf, -np.inf], np.nan).fillna(0).values
        X_comb = np.vstack([X_hc_s, X_dis])
        y_comb = np.concatenate([np.zeros(len(X_dis)), np.ones(len(X_dis))])

        clf = LogisticRegression(penalty="l1", solver="liblinear", C=0.5, random_state=42)
        clf.fit(X_comb, y_comb)
        coef = clf.coef_[0]
        top_idx = np.argsort(np.abs(coef))[-top_k:]
        all_genes.update([valid_genes[i] for i in top_idx if coef[i] != 0])

    selected = list(all_genes)
    col_colors = pd.Series(phenotypes).map(lut).values

    plot_data = df_z_dis[selected].replace([np.inf, -np.inf], np.nan).fillna(0).T
    g = sns.clustermap(
        plot_data,
        cmap="vlag", center=0, vmin=-4, vmax=4,
        col_colors=col_colors, figsize=(14, 12),
        method="ward", metric="euclidean",
        yticklabels=False, xticklabels=False,
    )
    legend_elements = [mlines.Line2D([0], [0], marker="s", color="w", label=lbl,
                                     markerfacecolor=c, markersize=10)
                       for lbl, c in lut.items()]
    g.ax_heatmap.legend(handles=legend_elements, title="Disease Groups",
                         bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.suptitle(f"Top-{top_k} Lasso Markers across Diseases (HC excluded)",
                 y=1.02, fontsize=15)
    plt.show()
    return selected
