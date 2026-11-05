import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import LabelEncoder

from analysis.plot_utils import _build_classifiers, _save


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
