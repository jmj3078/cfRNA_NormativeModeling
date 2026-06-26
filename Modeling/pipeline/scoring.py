import time

import numpy as np
import pandas as pd
from scipy.stats import norm as _norm, poisson as _poisson

import config
from model_engine import NormativeModelEngine
from rare_event_scorer import RareEventScorer

MP = config.MODELING_PARAMS


def load_engine(h5ad_path=None, engine_dir=None, rare_ref=None):
    """Load saved engine/rare_scorer, or train/build if absent."""
    h5ad_path = h5ad_path or config.H5AD_PATH
    engine_dir = engine_dir or config.ENGINE_DIR
    rare_ref = rare_ref or config.RARE_REF
    if (engine_dir / 'genes.pkl').exists():
        engine = NormativeModelEngine.load(engine_dir)
    else:
        engine = NormativeModelEngine(
            count_model='nbi', low_det_thr=MP['low_det_thr'], det_rate_min=MP['det_rate_min'],
            nbi_outlier_z=5.0, nbi_max_iter=2, nbi_max_remove_frac=0.10, lambda_sigma=0.05)
        engine.load_hc_data(h5ad_path)
        engine.assign_branches()
        engine.train(verbose=True)
        engine.save(engine_dir)
    rare_scorer = RareEventScorer.load(rare_ref) if rare_ref.exists() \
        else RareEventScorer.from_h5ad(h5ad_path)
    if not rare_ref.exists():
        rare_scorer.save(rare_ref)
    return engine, rare_scorer


def _rare_scores(rare_scorer, gene_names, Y_dis):
    """Return vectorized rare-event score matrix and aligned (genes, category) arrays."""
    ref = rare_scorer.ref
    r_genes = ref['gene'].values
    r_mean = ref['mean_hc'].values.astype(float)
    r_det = ref['det_rate_hc'].values.astype(float)
    r_cat = ref['category'].values
    gn_map = {g: i for i, g in enumerate(gene_names)}
    r_cols = np.array([gn_map.get(g, -1) for g in r_genes])
    valid = r_cols >= 0
    r_cols, r_genes, r_mean, r_det, r_cat = (
        r_cols[valid], r_genes[valid], r_mean[valid], r_det[valid], r_cat[valid])
    Y_rare = Y_dis[:, r_cols].astype(float)
    nonzero = Y_rare > 0
    silent = (r_det == 0)[np.newaxis, :]
    sc_rare = np.zeros_like(Y_rare, dtype=np.float32)
    sc_rare[silent & nonzero] = 10.0
    near = (~silent) & nonzero
    if near.any():
        k_n = Y_rare[near].astype(int)
        mu_n = np.broadcast_to(r_mean[np.newaxis, :], Y_rare.shape)[near]
        pv = _poisson.sf(k_n - 1, mu_n)
        sc_rare[near] = np.where(pv <= 0, 10.0, np.minimum(_norm.ppf(1.0 - pv), 10.0)).astype(np.float32)
    return sc_rare, Y_rare, r_genes, r_cat


def _engine_long(engine, gene_names, Z, Y_dis, sa_arr, ph_arr, min_abs_score):
    g_arr = np.array(gene_names)
    branch_of = np.array(['logistic' if (r := engine.genes.get(g)) and r.branch == 'logistic'
                          else 'count' for g in g_arr])
    stype_of = np.where(branch_of == 'logistic', 'logistic_z', 'nbi_z')
    mask = np.abs(Z) >= min_abs_score if min_abs_score > 0 else np.isfinite(Z)
    row_s, row_g = np.nonzero(mask)
    return pd.DataFrame({
        'sample': sa_arr[row_s], 'phenotype': ph_arr[row_s], 'gene': g_arr[row_g],
        'score': Z[row_s, row_g].astype(float), 'score_type': stype_of[row_g],
        'raw_count': Y_dis[row_s, row_g].astype(float), 'branch': branch_of[row_g]})


def score_all(engine, rare_scorer, gene_names, X_dis, Y_dis, sample_names, pheno_names,
              min_abs_score=3.0):
    """Return integrated long-format anomaly scores from engine and rare-event scorer."""
    sa_arr, ph_arr = np.array(sample_names), np.array(pheno_names)
    result = engine.score(X_dis, Y_dis, gene_names=gene_names, seed=42)
    Z = result['combined']
    df = _engine_long(engine, list(result['gene_names']), Z, Y_dis, sa_arr, ph_arr, min_abs_score)
    sc_rare, Y_rare, r_genes, r_cat = _rare_scores(rare_scorer, gene_names, Y_dis)
    rs, rg = np.nonzero(sc_rare >= max(min_abs_score, 0.01))
    if len(rs):
        df_rare = pd.DataFrame({
            'sample': sa_arr[rs], 'phenotype': ph_arr[rs], 'gene': r_genes[rg],
            'score': sc_rare[rs, rg].astype(float),
            'score_type': np.where(r_cat[rg] == 'silent', 'rare_fixed', 'rare_poisson'),
            'raw_count': Y_rare[rs, rg], 'branch': 'rare'})
        df = pd.concat([df, df_rare], ignore_index=True)
    return df


def score_full(engine, rare_scorer, gene_names, X_dis, Y_dis, dis_names, dis_pheno,
               thr=None, save=True):
    """Score all disease samples; save Z matrix and long-format parquet to Z_SCORES_DIR."""
    thr = MP['z_flag'] if thr is None else thr
    t0 = time.perf_counter()
    result_full = engine.score(X_dis, Y_dis, gene_names=gene_names, seed=42)
    Z_full = result_full['combined']
    print(f'Z matrix: {Z_full.shape}  ({time.perf_counter()-t0:.1f}s)')
    if save:
        config.Z_SCORES_DIR.mkdir(parents=True, exist_ok=True)
        np.save(config.Z_DISEASE, Z_full)
        np.save(config.Z_SAMPLE_NAMES, np.array(dis_names))
        np.save(config.Z_GENE_NAMES, np.array(gene_names))
    df = _engine_long(engine, list(result_full['gene_names']), Z_full, Y_dis,
                      np.array(dis_names), np.array(dis_pheno), thr)
    sc_rare, Y_rare, r_genes, r_cat = _rare_scores(rare_scorer, gene_names, Y_dis)
    rs, rg = np.nonzero(sc_rare >= thr)
    if len(rs):
        df_r = pd.DataFrame({
            'sample': np.array(dis_names)[rs], 'phenotype': np.array(dis_pheno)[rs],
            'gene': r_genes[rg], 'score': sc_rare[rs, rg].astype(float),
            'score_type': np.where(r_cat[rg] == 'silent', 'rare_fixed', 'rare_poisson'),
            'raw_count': Y_rare[rs, rg], 'branch': 'rare'})
        df = pd.concat([df, df_r], ignore_index=True)
    if save:
        config.Z_SCORES_DIR.mkdir(parents=True, exist_ok=True)
        df.to_parquet(config.Z_SCORES_DIR / 'disease_scores_flagged.parquet', index=False)
    return Z_full, df


def score_hc(engine, X_hc, Y_hc, gene_names, hc_names, save=True):
    """Score HC samples and optionally save Z_hc arrays to Z_SCORES_DIR."""
    res = engine.score(X_hc, Y_hc, gene_names=gene_names, seed=42)
    Z_hc = res['combined']
    if save:
        config.Z_SCORES_DIR.mkdir(parents=True, exist_ok=True)
        np.save(config.Z_HC, Z_hc)
        np.save(config.Z_HC_NAMES, np.array(hc_names))
    return Z_hc
