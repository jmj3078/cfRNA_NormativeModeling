import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             precision_recall_curve, roc_auc_score, roc_curve,
                             silhouette_samples)
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import NearestNeighbors
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


def _select_idx(Z, pheno, gene_names, method, n_per_pheno):
    """Run a GeneSelector method on an arbitrary row subset, return column indices.

    Single source of truth for selection logic so the nested-CV path cannot drift from
    the full-data selectors used elsewhere.
    """
    gs = GeneSelector(Z, pheno, gene_names)
    genes = {'proportion': gs.proportion,
             'effect_size': gs.effect_size,
             'svd': gs.svd_signature,
             'effect_size_specific': gs.effect_size_specific,
             'l1': gs.l1_logistic}[method](n_per_pheno=n_per_pheno)
    g2i = {g: i for i, g in enumerate(gene_names)}
    return [g2i[g] for g in genes if g in g2i]


CONTRAST_METHODS = ('effect_size_specific', 'l1')


def _fold_select(method, X, tr, yy, sub_ph, gene_names, n_per_pheno):
    """Leakage-free per-fold gene selection for the binary / calibration tasks.

    Magnitude selectors (proportion / effect_size / svd) rank within the positive rows
    only -- a within-disease ranking that matches their production use. The contrast /
    supervised selectors need >= 2 classes, so they select on the full train fold using
    the binary labels (disease vs negative); selection still touches train rows only, so
    no test-fold label leaks in.
    """
    if method in CONTRAST_METHODS:
        lab = np.where(yy[tr] == 1, '_d', '_hc')
        return _select_idx(X[tr], lab, gene_names, method, n_per_pheno)
    pos = tr[yy[tr] == 1]
    if len(pos) < 2:
        return None
    return _select_idx(X[pos], sub_ph[pos], gene_names, method, n_per_pheno)


def _new_curves(n, prevalence):
    """Per-fold ROC/PR curve accumulator for one binary task (familiar shadow-curve plots)."""
    return {'n': n, 'prevalence': prevalence,
            'roc': {'fprs': [], 'tprs': []}, 'pr': {'recs': [], 'precs': []}}


def _stash_curves(cv, y_true, score):
    """Append one fold's ROC and PR traces to a curve accumulator."""
    fpr, tpr, _ = roc_curve(y_true, score)
    cv['roc']['fprs'].append(fpr); cv['roc']['tprs'].append(tpr)
    prec, rec, _ = precision_recall_curve(y_true, score)
    cv['pr']['recs'].append(rec); cv['pr']['precs'].append(prec)


def _healthy_null(n, n_genes, rng):
    """Covariate-matched healthy reference for the normative Z. The RQR maps the
    healthy-at-covariates-X distribution to N(0,1) for EVERY X, so the matched reference
    needs no count sampling -- it is i.i.d. N(0,1) (genes are modelled independently).
    Using this as the negative class isolates disease deviation from the covariate /
    cohort differences that a real-HC negative would leak in."""
    return rng.standard_normal((n, n_genes)).astype(np.float32)


# negative-class modes for the binary task (user compares all three to pick one):
#   real_hc      -- disease vs the real HC cohort Z (the pre-healthy-null implementation)
#   null_fixed   -- disease vs N(0,1) healthy null, |neg| = |HC| (prevalence comparable)
#   null_matched -- disease vs N(0,1) healthy null, |neg| = n_disease (balanced case-control)
NEG_MODES = ('real_hc', 'null_fixed', 'null_matched')


def eval_binary_nested(Z_dis, dis_pheno, gene_names, method, neg_mode, Z_hc,
                       n_per_pheno=30, n_perm=20, seed=42):
    """Per-phenotype binary discrimination with gene selection nested INSIDE each CV fold,
    against one of NEG_MODES negatives.

    real_hc reproduces the earlier 'disease vs real HC cohort' result. null_fixed /
    null_matched use the covariate-matched healthy null: the Z is a randomized quantile
    residual, so 'healthy at covariates X' is N(0,1) for every X (sampling counts from the
    model and re-scoring just regenerates N(0,1)); using it as the negative isolates disease
    deviation from the covariate/cohort mismatch a real-HC negative leaks in. CAVEAT:
    covariate matching does NOT remove unmodelled study/batch -- see eval_multiclass_nested.
    Reports leakage-controlled AUC + AUPRC + balanced accuracy + label-permutation null.
    """
    rng = np.random.default_rng(seed)
    fixed_neg = Z_hc if neg_mode == 'real_hc' else (
        _healthy_null(len(Z_hc), Z_dis.shape[1], rng) if neg_mode == 'null_fixed' else None)
    rows, curves = [], {}
    for ph in np.unique(dis_pheno):
        m = dis_pheno == ph
        n = int(m.sum())
        if m.sum() < CV:
            rows.append({'phenotype': ph, 'auc': np.nan, 'auprc': np.nan,
                         'bal_acc': np.nan, 'perm_p': np.nan})
            curves[ph] = None
            continue
        neg = _healthy_null(n, Z_dis.shape[1], rng) if neg_mode == 'null_matched' else fixed_neg
        n_neg = len(neg)
        X = np.vstack([Z_dis[m], neg])
        y = np.r_[np.ones(m.sum()), np.zeros(n_neg)]
        sub_ph = np.array(['_d'] * m.sum() + ['_hc'] * n_neg)

        def _cv_auc(yy, cv=None):
            skf = StratifiedKFold(CV, shuffle=True, random_state=seed)
            a, ap, ba = [], [], []
            for tr, te in skf.split(X, yy):
                idx = _fold_select(method, X, tr, yy, sub_ph, gene_names, n_per_pheno)
                if not idx:
                    continue
                lr = LogisticRegression(max_iter=500, C=0.1, n_jobs=-1).fit(X[tr][:, idx], yy[tr])
                p = lr.predict_proba(X[te][:, idx])[:, 1]
                a.append(roc_auc_score(yy[te], p))
                ap.append(average_precision_score(yy[te], p))
                ba.append(balanced_accuracy_score(yy[te], (p >= 0.5).astype(int)))
                if cv is not None:
                    _stash_curves(cv, yy[te], p)
            return (np.mean(a), np.mean(ap), np.mean(ba)) if a else (np.nan,) * 3

        cv = _new_curves(n, float(y.mean()))
        auc, auprc, bal = _cv_auc(y, cv)
        null = np.array([_cv_auc(rng.permutation(y))[0] for _ in range(n_perm)])
        perm_p = (np.sum(null >= auc) + 1) / (n_perm + 1)
        rows.append({'phenotype': ph, 'auc': auc, 'auprc': auprc,
                     'bal_acc': bal, 'perm_p': perm_p})
        curves[ph] = cv
    return pd.DataFrame(rows).set_index('phenotype'), curves


def eval_multiclass_nested(Z_dis, dis_pheno, gene_names, method,
                           n_per_pheno=30, n_perm=20, seed=42):
    """Disease-vs-disease OVR multiclass with selection nested inside each CV fold.

    This is the non-circular discrimination test: separating one disease from the others
    is not trivially answered by 'anomalous vs HC'. Per-class AUC/AUPRC plus a
    label-permutation null. CAVEAT: in this cohort phenotype is near-perfectly confounded
    with study/batch (~1 study per disease), so high discrimination may reflect residual
    technical structure, not biology -- interpret with that limitation stated.
    """
    rng = np.random.default_rng(seed)
    le = LabelEncoder()
    y = le.fit_transform(dis_pheno)
    classes = le.classes_

    def _cv_perclass(yy, cv=None):
        skf = StratifiedKFold(CV, shuffle=True, random_state=seed)
        auc = {c: [] for c in range(len(classes))}
        ap = {c: [] for c in range(len(classes))}
        for tr, te in skf.split(Z_dis, yy):
            ph_tr = np.array([str(v) for v in yy[tr]])
            idx = _select_idx(Z_dis[tr], ph_tr, gene_names, method, n_per_pheno)
            clf = LogisticRegression(max_iter=500, C=0.1, multi_class='ovr', n_jobs=-1)
            clf.fit(Z_dis[tr][:, idx], yy[tr])
            proba = clf.predict_proba(Z_dis[te][:, idx])
            present = {c: j for j, c in enumerate(clf.classes_)}
            for c in range(len(classes)):
                yb = (yy[te] == c).astype(int)
                if yb.sum() == 0 or c not in present:
                    continue
                auc[c].append(roc_auc_score(yb, proba[:, present[c]]))
                ap[c].append(average_precision_score(yb, proba[:, present[c]]))
                if cv is not None:
                    _stash_curves(cv[c], yb, proba[:, present[c]])
        return ({c: np.mean(v) if v else np.nan for c, v in auc.items()},
                {c: np.mean(v) if v else np.nan for c, v in ap.items()})

    counts = np.bincount(y, minlength=len(classes))
    cv = {c: _new_curves(int(counts[c]), float(counts[c] / len(y))) for c in range(len(classes))}
    auc, ap = _cv_perclass(y, cv)
    null = {c: [] for c in range(len(classes))}
    for _ in range(n_perm):
        na, _ = _cv_perclass(rng.permutation(y))
        for c in range(len(classes)):
            null[c].append(na[c])
    rows, curves = [], {}
    for c, ph in enumerate(classes):
        nd = np.array([v for v in null[c] if np.isfinite(v)])
        p = (np.sum(nd >= auc[c]) + 1) / (len(nd) + 1) if len(nd) else np.nan
        rows.append({'phenotype': ph, 'auc_mc': auc[c], 'auprc_mc': ap[c],
                     'perm_p_mc': p})
        curves[ph] = cv[c]
    return pd.DataFrame(rows).set_index('phenotype'), curves


def calibration_control_hc(Z_hc, gene_names, method, n_per_pheno=30, seed=42):
    """Calibration sanity check: real HC Z vs the covariate-matched healthy null N(0,1),
    selection nested in folds. If the normative model is well calibrated the real HC are
    indistinguishable from their own healthy null -> AUC ~ 0.5. AUC well above 0.5 means
    the model is miscalibrated on HC (e.g. residual over-dispersion makes some genes
    consistently deviate), which would also inflate the disease contrast."""
    rng = np.random.default_rng(seed)
    n = len(Z_hc)
    neg = _healthy_null(n, Z_hc.shape[1], rng)
    X = np.vstack([Z_hc, neg])
    y = np.r_[np.ones(n), np.zeros(n)]
    sub_ph = np.array(['_d'] * n + ['_hc'] * n)
    skf = StratifiedKFold(CV, shuffle=True, random_state=seed)
    cv = _new_curves(n, 0.5)
    a = []
    for tr, te in skf.split(X, y):
        idx = _fold_select(method, X, tr, y, sub_ph, gene_names, n_per_pheno)
        if not idx:
            continue
        lr = LogisticRegression(max_iter=500, C=0.1, n_jobs=-1).fit(X[tr][:, idx], y[tr])
        p = lr.predict_proba(X[te][:, idx])[:, 1]
        a.append(roc_auc_score(y[te], p))
        _stash_curves(cv, y[te], p)
    return float(np.mean(a)), float(np.std(a)), cv


def run_validation(Z_dis, Z_hc, dis_pheno, gene_names, n_per_pheno=30, n_perm=20,
                   methods=('proportion', 'effect_size', 'svd',
                            'effect_size_specific', 'l1')):
    """Leakage-controlled validation across selectors. The binary task is computed under
    all three NEG_MODES (real_hc / null_fixed / null_matched) so they can be compared.

    Returns (binary_df, multiclass_df, controls, curves). binary_df has a 'neg_mode' column.
    controls[method] = (calib_auc_mean, std, hc_curve), calib = real-HC-vs-null (~0.5 =
    calibrated). curves[method] = {'binary': {neg_mode: {ph: curve}}, 'multi': {ph: curve},
    'hc': hc_curve} for the familiar shadow ROC/PR grids (plots.plot_validation_curves)."""
    dis_pheno = np.array(dis_pheno)
    binary, multi, controls, curves = [], {}, {}, {}
    for method in methods:
        ctrl_auc, ctrl_sd, hc_cv = calibration_control_hc(Z_hc, gene_names, method, n_per_pheno)
        controls[method] = (ctrl_auc, ctrl_sd, hc_cv)
        mc_df, mc_cv = eval_multiclass_nested(Z_dis, dis_pheno, gene_names, method,
                                              n_per_pheno, n_perm)
        multi[method] = mc_df
        bin_curves = {}
        med = {}
        for nm in NEG_MODES:
            bin_df, bin_cv = eval_binary_nested(Z_dis, dis_pheno, gene_names, method, nm,
                                                Z_hc, n_per_pheno, n_perm)
            bin_df = bin_df.reset_index()
            bin_df.insert(0, 'neg_mode', nm)
            bin_df.insert(0, 'method', method)
            binary.append(bin_df)
            bin_curves[nm] = bin_cv
            med[nm] = bin_df['auc'].median()
        curves[method] = {'binary': bin_curves, 'multi': mc_cv, 'hc': hc_cv}
        print(f'{method}: calib(HC-vs-null)={ctrl_auc:.3f}  median binary AUC '
              f'[real_hc={med["real_hc"]:.3f} null_fixed={med["null_fixed"]:.3f} '
              f'null_matched={med["null_matched"]:.3f}]  median mc AUC={mc_df["auc_mc"].median():.3f}')
    binary_df = pd.concat(binary, ignore_index=True)
    multi_df = pd.concat(multi, names=['method']).reset_index()
    return binary_df, multi_df, controls, curves


def run_selection(Z_dis, Z_hc, dis_pheno, gene_names, n_per_pheno=30):
    """Return (all_results dict, selectors) with the selected gene set and leakage-free
    unsupervised structure (silhouette / kNN purity) per selector.

    Supervised classifier performance is NOT computed here -- the old eval_binary /
    eval_multiclass selected genes on the full labelled set and scored them by CV
    (selection-bias leakage). Use run_validation() for the nested-CV, leakage-controlled
    AUC / AUPRC instead. Z_hc is kept in the signature for backward compatibility.
    """
    import time
    gs = GeneSelector(Z_dis, dis_pheno, gene_names)
    selectors = gs.get_selectors(n_per_pheno=n_per_pheno)
    _, idx_of = gene_index(gene_names)
    all_results = {}
    for name, selector in selectors.items():
        t0 = time.perf_counter()
        genes = selector()
        idx = idx_of(genes)
        per_pheno, macro_sil = eval_unsupervised(idx, Z_dis, dis_pheno)
        per_pheno.insert(0, 'n', pd.Series(dis_pheno).value_counts())
        all_results[name] = dict(
            genes=genes, per_pheno=per_pheno, macro_sil=macro_sil,
            n_genes=len(genes), t=time.perf_counter() - t0)
        print(f'{name}: sil={macro_sil:.3f}  n_genes={len(genes)}  ({all_results[name]["t"]:.0f}s)')
    return all_results, selectors
