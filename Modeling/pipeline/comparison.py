import pickle

import numpy as np
import pandas as pd
from scipy.stats import norm

import config

PALETTE = {'NBI': '#E41A1C', 'ZINBI': '#377EB8', 'Logistic': '#4DAF4A'}
MODEL_ORDER = ['NBI', 'ZINBI', 'Logistic']
ETP_THR = 1.96
W1_EXPECTED = 0.048
METRIC_DEFS = [
    ('w1',     'W1 ↓',           (0.0,  0.40),  W1_EXPECTED),
    ('etp',    'ETP ↓ (→0.05)', (0.0,  0.20),  0.05),
    ('skew_z', 'Skewness → 0',  (-2.0, 2.0),   0.0),
    ('kurt_z', 'Kurtosis → 0',  (-1.0, 6.0),   0.0),
    ('mean_z', 'Z-mean → 0',    (-0.4, 0.4),   0.0),
    ('std_z',  'Z-std → 1',     (0.6,  1.4),   1.0),
]


def load_csv(path):
    return pd.read_csv(path) if path.exists() else None


def load_pkl(path):
    if not path.exists():
        return None
    with open(path, 'rb') as f:
        return pickle.load(f)


def load_cv_outputs(cv_dir=None):
    """cv_*.py 산출 (stats csv + zscores pkl) 로드 → (nb, zinb, logr, nb_z, zinb_z, logr_z)."""
    cv_dir = cv_dir or config.CV_RESULTS_DIR
    nb = load_csv(cv_dir / 'cv_gamlss_stats.csv')
    zinb = load_csv(cv_dir / 'cv_zinb_stats.csv')
    logr = load_csv(cv_dir / 'cv_logistic_stats.csv')
    nb_z = load_pkl(cv_dir / 'cv_gamlss_zscores.pkl')
    zinb_z = load_pkl(cv_dir / 'cv_zinb_zscores.pkl')
    logr_z = load_pkl(cv_dir / 'cv_logistic_zscores.pkl')
    return nb, zinb, logr, nb_z, zinb_z, logr_z


def _get_zarr(z_dict, gene, z_type):
    if z_dict is None:
        return None
    entry = z_dict.get(gene)
    if entry is None:
        return None
    return entry.get(z_type) if isinstance(entry, dict) else entry


def w1(arr):
    v = arr[np.isfinite(arr)]
    n = len(v)
    if n < 8:
        return np.nan
    ref = norm.ppf(np.linspace(1 / (2 * n), 1 - 1 / (2 * n), n))
    return float(np.mean(np.abs(np.sort(v) - ref)))


def etp(arr, thr=ETP_THR):
    v = arr[np.isfinite(arr)]
    return float(np.mean(np.abs(v) > thr)) if len(v) >= 8 else np.nan


def extract_metrics(df, label, mz_col, sz_col, sk_col, ku_col, success_col,
                    z_dict, z_type='full', w1_col=None, flag_col=None):
    if df is None:
        return None
    genes = df['gene'].values
    n = len(genes)
    w1_arr = df[w1_col].values.astype(float) if (w1_col and w1_col in df.columns) else np.full(n, np.nan)
    etp_arr = np.full(n, np.nan)
    for i, gene in enumerate(genes):
        arr = _get_zarr(z_dict, gene, z_type)
        if arr is None:
            continue
        if np.isnan(w1_arr[i]):
            w1_arr[i] = w1(arr)
        etp_arr[i] = etp(arr)
    return pd.DataFrame({
        'model': label, 'gene': genes,
        'det_rate': df['det_rate_hc'].values.astype(float),
        'mean_count': df['mean_count_hc'].values.astype(float) if 'mean_count_hc' in df.columns else np.full(n, np.nan),
        'w1': w1_arr, 'etp': etp_arr,
        'mean_z': df[mz_col].values.astype(float), 'std_z': df[sz_col].values.astype(float),
        'skew_z': df[sk_col].values.astype(float), 'kurt_z': df[ku_col].values.astype(float),
        'fold_ok': df[success_col].values.astype(float),
        'any_flag': df[flag_col].astype(bool).values if (flag_col and flag_col in df.columns) else np.zeros(n, bool),
        'n_removed': df['n_removed'].values.astype(float) if 'n_removed' in df.columns else np.zeros(n),
    })


def build_all_df(nb, zinb, logr, nb_z, zinb_z, logr_z):
    frames = []
    if nb is not None:
        frames.append(extract_metrics(nb, 'NBI', 'mean_z', 'std_z', 'skew_z', 'kurt_z',
                                       'fold_success_rate', nb_z, 'full', w1_col='w1'))
    if zinb is not None:
        frames.append(extract_metrics(zinb, 'ZINBI', 'mean_full_z', 'std_full_z', 'skew_full_z',
                                       'kurt_full_z', 'fold_success_rate', zinb_z, 'full',
                                       w1_col='w1_full', flag_col='any_flag'))
    if logr is not None:
        frames.append(extract_metrics(logr, 'Logistic', 'mean_z', 'std_z', 'skew_z', 'kurt_z',
                                       'fold_success_rate', logr_z, 'full', flag_col='any_flag'))
    all_df = pd.concat(frames, ignore_index=True)
    models_present = [m for m in MODEL_ORDER if m in all_df['model'].unique()]
    return all_df, models_present


# ── PPC simulation (Posterior Predictive Check) ─────────────────────────────
def sim_nbi(mu, sigma, seed=42):
    rng = np.random.default_rng(seed)
    theta = 1.0 / np.maximum(sigma, 1e-8)
    p = np.clip(theta / (theta + mu), 1e-8, 1 - 1e-8)
    return rng.negative_binomial(theta, p).astype(np.float32)


def sim_zinbi(mu, sigma, nu, seed=42):
    rng = np.random.default_rng(seed)
    theta = 1.0 / np.maximum(sigma, 1e-8)
    p = np.clip(theta / (theta + mu), 1e-8, 1 - 1e-8)
    y_nbi = rng.negative_binomial(theta, p)
    return np.where(rng.uniform(size=len(mu)) < np.clip(nu, 0, 1), 0, y_nbi).astype(np.float32)


SIM_FN = {
    'NBI':   lambda params, seed: sim_nbi(params['mu'], params['sigma'], seed),
    'ZINBI': lambda params, seed: sim_zinbi(params['mu'], params['sigma'], params['nu'], seed),
}
