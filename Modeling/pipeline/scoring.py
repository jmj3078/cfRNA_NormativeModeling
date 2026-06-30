import time

import numpy as np
import pandas as pd

import config
from model_engine import NormativeModelEngine

MP = config.MODELING_PARAMS


def load_engine(h5ad_path=None, engine_dir=None):
    """Load saved engine (rare branch folded in), or train/build if absent."""
    h5ad_path = h5ad_path or config.H5AD_PATH
    engine_dir = engine_dir or config.ENGINE_DIR
    if (engine_dir / 'genes.pkl').exists():
        return NormativeModelEngine.load(engine_dir)
    engine = NormativeModelEngine(
        count_model='nbi', low_det_thr=MP['low_det_thr'], det_rate_min=MP['det_rate_min'],
        nbi_outlier_z=5.0, nbi_max_iter=2, nbi_max_remove_frac=0.10, lambda_sigma=0.05)
    engine.load_hc_data(h5ad_path)
    engine.assign_branches()
    engine.train(verbose=True)
    engine.save(engine_dir)
    return engine


def _scores_long(engine, gene_names, Z_all, Y_dis, sa_arr, ph_arr, min_abs_score):
    """Unified long-format flagged table across all branches (logistic/count/rare).

    Built from the full Z (combined_all, rare included). NaN cells (unfitted genes) are
    never flagged since NaN >= thr is False.
    """
    g_arr = np.array(gene_names)
    raw_branch = np.array([engine.genes[g].branch if g in engine.genes else 'none'
                           for g in g_arr])
    # historical parquet labels: logistic | count | rare (NBI/ZINBI collapse to 'count')
    branch_of = np.where(raw_branch == 'logistic', 'logistic',
                         np.where(raw_branch == 'rare', 'rare', 'count'))
    stype_of = np.where(raw_branch == 'logistic', 'logistic_z',
                        np.where(raw_branch == 'rare', 'rare_glm', 'nbi_z'))
    mask = np.abs(Z_all) >= min_abs_score if min_abs_score > 0 else np.isfinite(Z_all)
    row_s, row_g = np.nonzero(mask)
    return pd.DataFrame({
        'sample': sa_arr[row_s], 'phenotype': ph_arr[row_s], 'gene': g_arr[row_g],
        'score': Z_all[row_s, row_g].astype(float), 'score_type': stype_of[row_g],
        'raw_count': Y_dis[row_s, row_g].astype(float), 'branch': branch_of[row_g]})


def _save_rare(result, gene_names, z_rare_path, gene_names_path):
    """Persist the rare submatrix (aligned to gene_names order) as a unified artifact."""
    gidx = {g: i for i, g in enumerate(gene_names)}
    rare_genes = result['rare_gene_names']
    Z_rare = np.full((result['combined_all'].shape[0], len(rare_genes)), 0.0, np.float32)
    src = {g: j for j, g in enumerate(gene_names)}
    for j, g in enumerate(rare_genes):
        Z_rare[:, j] = result['combined_all'][:, src[g]]
    np.save(z_rare_path, Z_rare)
    np.save(gene_names_path, np.array(rare_genes))
    return Z_rare, rare_genes


def score_all(engine, gene_names, X_dis, Y_dis, sample_names, pheno_names, min_abs_score=3.0):
    """Long-format anomaly scores across all branches (logistic/count/rare) for a few
    samples. In-memory only -- used by the single-sample inspection notebook."""
    result = engine.score(X_dis, Y_dis, gene_names=gene_names, seed=42)
    return _scores_long(engine, list(result['gene_names']), result['combined_all'], Y_dis,
                        np.array(sample_names), np.array(pheno_names), min_abs_score)


def score_full(engine, gene_names, X_dis, Y_dis, dis_names, dis_pheno, thr=None, save=True):
    """Score all disease samples.

    Saves the canonical engine-only Z matrix (Z_disease.npy, rare columns zeroed -- the
    historical placeholder contract, preserved so existing downstream is untouched), the
    unified rare artifact (Z_rare_disease.npy), and the long-format flagged parquet (all
    branches, rare scored by the pooled covariate GLM). Returns (Z_combined, Z_rare, df).
    """
    thr = MP['z_flag'] if thr is None else thr
    t0 = time.perf_counter()
    result = engine.score(X_dis, Y_dis, gene_names=gene_names, seed=42)
    Z_combined = result['combined']
    print(f'Z matrix: {Z_combined.shape}  ({time.perf_counter()-t0:.1f}s)')

    df = _scores_long(engine, list(result['gene_names']), result['combined_all'], Y_dis,
                      np.array(dis_names), np.array(dis_pheno), thr)
    if save:
        config.Z_SCORES_DIR.mkdir(parents=True, exist_ok=True)
        np.save(config.Z_DISEASE, Z_combined)
        np.save(config.Z_SAMPLE_NAMES, np.array(dis_names))
        np.save(config.Z_GENE_NAMES, np.array(gene_names))
        Z_rare, _ = _save_rare(result, list(result['gene_names']),
                               config.Z_RARE_DISEASE, config.Z_RARE_GENE_NAMES)
        df.to_parquet(config.Z_SCORES_DIR / 'disease_scores_flagged.parquet', index=False)
    else:
        Z_rare = result['rare']
    return Z_combined, Z_rare, df


def score_hc(engine, X_hc, Y_hc, gene_names, hc_names, save=True):
    """Score HC samples; save engine-only Z_hc and the unified rare HC artifact."""
    result = engine.score(X_hc, Y_hc, gene_names=gene_names, seed=42)
    Z_hc = result['combined']
    if save:
        config.Z_SCORES_DIR.mkdir(parents=True, exist_ok=True)
        np.save(config.Z_HC, Z_hc)
        np.save(config.Z_HC_NAMES, np.array(hc_names))
        _save_rare(result, list(result['gene_names']), config.Z_RARE_HC,
                   config.Z_RARE_GENE_NAMES)
    return Z_hc


def load_z(with_rare=False):
    """Load the canonical disease Z (engine-only). When with_rare=True, overlay the saved
    rare covariate scores onto their gene columns. Returns (Z, dis_names, gene_names)."""
    from pipeline import data_prep
    Z, dis_names, gene_names = data_prep.load_z_disease()
    if with_rare and config.Z_RARE_DISEASE.exists():
        Z = _overlay_rare(Z, gene_names, dis_names)
    return Z, dis_names, gene_names


def _overlay_rare(Z, gene_names, row_names):
    """Overlay the saved rare disease scores onto a copy of Z, aligned by sample name."""
    Z = Z.copy()
    all_names = np.load(config.Z_SAMPLE_NAMES, allow_pickle=True).tolist()
    rare_genes = np.load(config.Z_RARE_GENE_NAMES, allow_pickle=True).tolist()
    Z_rare = np.load(config.Z_RARE_DISEASE)
    row_of = {n: i for i, n in enumerate(all_names)}
    rows = [row_of[n] for n in row_names]
    Z_rare = Z_rare[rows]
    gidx = {g: i for i, g in enumerate(gene_names)}
    for j, g in enumerate(rare_genes):
        if g in gidx:
            Z[:, gidx[g]] = Z_rare[:, j]
    return Z


def score_disease_with_rare(dd, engine=None):
    """Disease Z matrix with the saved rare covariate scores overlaid, aligned to dd's
    OOD/min-sample-filtered sample order. In-memory only -- never writes to disk. For an
    explicit, opt-in 'with_rare' GSEA / gene-selection run compared side-by-side against
    the canonical (engine-only) dd.Z_dis."""
    return _overlay_rare(dd.Z_dis, dd.gene_names, dd.dis_names)
