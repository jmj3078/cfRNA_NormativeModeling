#!/usr/bin/env python
"""K-fold calibration benchmark for the pooled rare-event covariate GLM.

Mirrors cv_logistic.py's output schema so pipeline/comparison.py can include 'Rare'
alongside NBI/ZINBI/Logistic. One pooled GLM (offset=log(mean_hc_j+eps), shared beta
across all rare genes) is fit per train fold; held-out HC samples are scored with the
project-standard RQR. Poisson first, NB only if pooled deviance/df exceeds the lenient
rare_overdisp_thr (per-fold). The per-gene mean_hc offset is recomputed from the train
fold so it never leaks into the held-out score.

Output (CV_Results/):
    cv_rare_stats.csv    -- per-gene det_rate_hc, mean_count_hc, z moments, w1, any_flag
    cv_rare_zscores.pkl  -- dict {gene: held-out z array (length N_hc)}
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import kurtosis, nbinom, norm, poisson, skew
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from statsmodels.discrete.discrete_model import NegativeBinomial

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
import config
from pipeline import data_prep

MP = config.MODELING_PARAMS
N_SPLITS = MP['n_splits']
Z_FLAG = MP['z_flag']
RARE_DET_MAX = MP['rare_det_max']
OVERDISP_THR = MP['rare_overdisp_thr']
SEED = 42


def _rqr(y, mu, family, alpha, seed):
    y = np.asarray(y)
    if family == 'poisson':
        lo = np.where(y > 0, poisson.cdf(y - 1, mu), 0.0)
        hi = poisson.cdf(y, mu)
    else:
        n = 1.0 / alpha
        p = np.clip(n / (n + mu), 1e-12, 1 - 1e-12)
        lo = np.where(y > 0, nbinom.cdf(y - 1, n, p), 0.0)
        hi = nbinom.cdf(y, n, p)
    lo = np.clip(lo, 1e-12, 1 - 1e-12)
    hi = np.clip(hi, 1e-12, 1 - 1e-12)
    rng = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(np.minimum(lo, hi), np.maximum(lo, hi)))


def _w1_to_normal(z):
    v = z[np.isfinite(z)]
    n = len(v)
    if n < 8:
        return np.nan
    ref = norm.ppf(np.linspace(1 / (2 * n), 1 - 1 / (2 * n), n))
    return float(np.mean(np.abs(np.sort(v) - ref)))


def _fit_pooled(X, Y, mean_hc, eps):
    n, g = Y.shape
    si = np.repeat(np.arange(n), g)
    gi = np.tile(np.arange(g), n)
    Xc = np.column_stack([np.ones(n * g), X[si]])
    y = Y[si, gi]
    off = np.log(mean_hc[gi] + eps)
    pois = sm.GLM(y, Xc, family=sm.families.Poisson(), offset=off).fit()
    ratio = float(pois.deviance / pois.df_resid)
    if ratio <= OVERDISP_THR:
        return 'poisson', np.asarray(pois.params), None, ratio
    nb = NegativeBinomial(y, Xc, offset=off).fit(disp=False)
    return 'negbin', np.asarray(nb.params[:-1]), float(nb.params[-1]), ratio


def main():
    adata = data_prep.load_adata()
    is_hc, _, _ = data_prep.make_phenotypes(adata)
    X = data_prep.bias_matrix(adata)[is_hc]
    Y_all = data_prep.count_matrix(adata)[is_hc]
    is_pc = (adata.var['GeneType'] == 'protein_coding').values
    det = (Y_all > 0).mean(axis=0)
    mean_cnt = Y_all.mean(axis=0)
    rare_mask = (det < RARE_DET_MAX) & is_pc
    genes = np.array(adata.var_names)[rare_mask]
    Y = Y_all[:, rare_mask]
    det_g = det[rare_mask]
    mean_g = mean_cnt[rare_mask]

    n_hc = X.shape[0]
    eps = 1.0 / (2 * n_hc)
    Xs = StandardScaler().fit_transform(X)
    print(f'HC={n_hc}  rare genes={len(genes)}  eps={eps:.6f}')

    z_acc = np.full(Y.shape, np.nan, dtype=np.float64)
    kf = KFold(N_SPLITS, shuffle=True, random_state=SEED)
    for fold, (tr, te) in enumerate(kf.split(Xs)):
        mean_tr = Y[tr].mean(axis=0)
        family, beta, alpha, ratio = _fit_pooled(Xs[tr], Y[tr], mean_tr, eps)
        eta = np.column_stack([np.ones(len(te)), Xs[te]]) @ beta
        mu = np.clip((mean_tr[None, :] + eps) * np.exp(eta[:, None]), 1e-12, 1e8)
        z = _rqr(Y[te].ravel(), mu.ravel(), family, alpha, SEED + fold).reshape(mu.shape)
        z_acc[te] = z
        print(f'  fold {fold}: family={family}  deviance/df={ratio:.3f}')

    rows, zdict = [], {}
    for j, g in enumerate(genes):
        zj_full = z_acc[:, j].astype(np.float32)
        zdict[g] = zj_full
        zj = zj_full[np.isfinite(zj_full)]
        rows.append({
            'gene': g, 'det_rate_hc': float(det_g[j]), 'mean_count_hc': float(mean_g[j]),
            'mean_z': float(zj.mean()), 'std_z': float(zj.std()),
            'skew_z': float(skew(zj)), 'kurt_z': float(kurtosis(zj)),
            'w1': _w1_to_normal(zj), 'fold_success_rate': 1.0,
            'any_flag': bool((np.abs(zj) >= Z_FLAG).any()), 'n_removed': 0,
        })

    config.CV_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(config.CV_RESULTS_DIR / 'cv_rare_stats.csv', index=False)
    with open(config.CV_RESULTS_DIR / 'cv_rare_zscores.pkl', 'wb') as f:
        pickle.dump(zdict, f)
    print(f'Saved -> {config.CV_RESULTS_DIR}/cv_rare_stats.csv  + cv_rare_zscores.pkl')


if __name__ == '__main__':
    main()
