import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from scipy import stats as scipy_stats
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from statsmodels.stats.multitest import multipletests
from tqdm.auto import tqdm


# ---------------------------------------------------------------------------
# GP Model Evaluation
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


# ---------------------------------------------------------------------------
# Normative Marker Discovery
# ---------------------------------------------------------------------------

def plot_normative_markers(z_results, adata_disease, phenotype_col="Phenotype_Processed",
                            fdr_th=0.01, z_abs_min=2.5, top_k=20):
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


# ---------------------------------------------------------------------------
# Lasso Classification
# ---------------------------------------------------------------------------

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
