import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_samples, roc_auc_score, roc_curve
from sklearn.neighbors import NearestNeighbors
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder

import config
from gene_selectors import GeneSelector

MP = config.MODELING_PARAMS
CV = MP['n_splits']
K = 15
BASE_FPR = np.linspace(0, 1, 101)


def gene_index(gene_names):
    g2i = {g: i for i, g in enumerate(gene_names)}
    return g2i, (lambda genes: [g2i[g] for g in genes if g in g2i])


def eval_unsupervised(idx, Z_dis, dis_pheno):
    X = Z_dis[:, idx]
    sil = silhouette_samples(X, dis_pheno, metric='euclidean')
    nbrs = NearestNeighbors(n_neighbors=K + 1, metric='euclidean', n_jobs=-1).fit(X)
    nn = nbrs.kneighbors(X)[1][:, 1:]
    pur = np.array([np.mean(dis_pheno[nn[i]] == dis_pheno[i]) for i in range(len(dis_pheno))])
    rows = [{'phenotype': ph,
             'silhouette': float(sil[dis_pheno == ph].mean()),
             'knn_purity': float(pur[dis_pheno == ph].mean())}
            for ph in np.unique(dis_pheno)]
    return pd.DataFrame(rows).set_index('phenotype'), float(sil.mean())


def eval_binary(idx, Z_dis, Z_hc, dis_pheno):
    """HC vs Disease per phenotype — 5-fold CV, LogReg + RF."""
    Xhc = Z_hc[:, idx]
    skf = StratifiedKFold(n_splits=CV, shuffle=True, random_state=42)
    lr = LogisticRegression(max_iter=500, C=0.1)
    rf = RandomForestClassifier(n_estimators=30, max_depth=6, random_state=42, n_jobs=-1)
    rows, roc_dict = [], {}
    for ph in np.unique(dis_pheno):
        m = dis_pheno == ph
        if m.sum() < CV:
            rows.append({'phenotype': ph, 'auc_logreg': np.nan, 'auc_rf': np.nan})
            roc_dict[ph] = None
            continue
        X = np.vstack([Z_dis[m][:, idx], Xhc])
        y = np.array([1] * m.sum() + [0] * len(Xhc))
        lr_a, rf_a = [], []
        lr_roc = {'fprs': [], 'tprs': []}
        rf_roc = {'fprs': [], 'tprs': []}
        for tr, te in skf.split(X, y):
            try:
                lr.fit(X[tr], y[tr]); p = lr.predict_proba(X[te])[:, 1]
                fpr, tpr, _ = roc_curve(y[te], p)
                lr_roc['fprs'].append(fpr); lr_roc['tprs'].append(tpr)
                lr_a.append(roc_auc_score(y[te], p))
                rf.fit(X[tr], y[tr]); p = rf.predict_proba(X[te])[:, 1]
                fpr, tpr, _ = roc_curve(y[te], p)
                rf_roc['fprs'].append(fpr); rf_roc['tprs'].append(tpr)
                rf_a.append(roc_auc_score(y[te], p))
            except Exception:
                pass
        roc_dict[ph] = {'lr': lr_roc, 'rf': rf_roc}
        rows.append({'phenotype': ph,
                     'auc_logreg': float(np.mean(lr_a)) if lr_a else np.nan,
                     'auc_rf': float(np.mean(rf_a)) if rf_a else np.nan})
    return pd.DataFrame(rows).set_index('phenotype'), roc_dict


def eval_multiclass(idx, Z_dis, dis_pheno):
    """All-disease OVR multiclass — 5-fold CV LogReg."""
    le = LabelEncoder()
    y = le.fit_transform(dis_pheno)
    X = Z_dis[:, idx]
    skf = StratifiedKFold(n_splits=CV, shuffle=True, random_state=42)
    clf = LogisticRegression(max_iter=500, C=0.1, multi_class='ovr')
    fold_yt, fold_yp = [], []
    for tr, te in skf.split(X, y):
        clf.fit(X[tr], y[tr])
        fold_yt.append(y[te]); fold_yp.append(clf.predict_proba(X[te]))
    rows, roc_dict = [], {}
    for i, ph in enumerate(le.classes_):
        class_roc = {'fprs': [], 'tprs': []}
        aucs = []
        for yt, yp in zip(fold_yt, fold_yp):
            try:
                y_bin = (yt == i).astype(int)
                fpr, tpr, _ = roc_curve(y_bin, yp[:, i])
                class_roc['fprs'].append(fpr); class_roc['tprs'].append(tpr)
                aucs.append(roc_auc_score(y_bin, yp[:, i]))
            except Exception:
                pass
        roc_dict[ph] = class_roc
        rows.append({'phenotype': ph, 'auc_multiclass': float(np.mean(aucs)) if aucs else np.nan})
    df = pd.DataFrame(rows).set_index('phenotype')
    return df, float(df['auc_multiclass'].mean()), roc_dict


def run_selection(Z_dis, Z_hc, dis_pheno, gene_names, n_per_pheno=30):
    """SELECTORS별 unsup/binary/multiclass 평가 → all_results dict."""
    import time
    gs = GeneSelector(Z_dis, dis_pheno, gene_names)
    selectors = gs.get_selectors(n_per_pheno=n_per_pheno)
    _, idx_of = gene_index(gene_names)
    all_results = {}
    for name, selector in selectors.items():
        t0 = time.perf_counter()
        genes = selector()
        idx = idx_of(genes)
        df_unsup, macro_sil = eval_unsupervised(idx, Z_dis, dis_pheno)
        df_bin, roc_bin = eval_binary(idx, Z_dis, Z_hc, dis_pheno)
        df_mc, macro_mc, roc_mc = eval_multiclass(idx, Z_dis, dis_pheno)
        per_pheno = df_unsup.join(df_bin, how='left').join(df_mc, how='left')
        per_pheno.insert(0, 'n', pd.Series(dis_pheno).value_counts())
        all_results[name] = dict(
            genes=genes, per_pheno=per_pheno, macro_sil=macro_sil, macro_mc=macro_mc,
            roc_bin=roc_bin, roc_mc=roc_mc, n_genes=len(genes), t=time.perf_counter() - t0)
        print(f'{name}: sil={macro_sil:.3f}  mc_auc={macro_mc:.3f}  ({all_results[name]["t"]:.0f}s)')
    return all_results, selectors
