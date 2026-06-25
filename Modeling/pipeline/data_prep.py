from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

import config
from sample_filter import MahalanobisFilter

MP = config.MODELING_PARAMS


def first_author(s):
    return str(s).replace('\xa0', ' ').split()[0]


def load_adata(h5ad_path=None, protein_coding=True):
    adata = sc.read_h5ad(h5ad_path or config.H5AD_PATH)
    obs_mask = (
        (adata.obs['QC_Passed'] == True) &
        (adata.obs['Phenotype_Processed'].notna()) &
        (adata.obs['Phenotype_Processed'] != 'Unknown') &
        (adata.obs['broad_protocol_category'] != 'Exome-based (EB)')
    )
    var_mask = (adata.var['GeneType'] == 'protein_coding').values if protein_coding else slice(None)
    return adata[obs_mask, var_mask].copy()


def make_phenotypes(adata):
    """HC mask + study-split phenotype labels (multi-author -> '(FirstAuthor)' suffix)."""
    is_hc = (adata.obs['Phenotype_Processed'].astype(str) == 'Healthy Control').values
    phenos = adata.obs['Phenotype_Processed'].astype(str).values
    authors = adata.obs['Author'].astype(str).values
    pmap = defaultdict(set)
    for ph, au, hc in zip(phenos, authors, is_hc):
        if not hc:
            pmap[ph.strip()].add(first_author(au))
    multi = {ph for ph, aus in pmap.items() if len(aus) > 1}
    phenos = np.array([
        ph if ph == 'Healthy Control' else
        (f"{ph.strip()} ({first_author(au)})" if ph.strip() in multi else ph.strip())
        for ph, au in zip(phenos, authors)
    ])
    return is_hc, phenos, sorted(multi)


def bias_matrix(adata):
    return adata.obs[config.BIAS_COLUMNS].values.astype(np.float64)


def count_matrix(adata):
    X = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
    return np.round(X).astype(np.float32)


def clean_z(Z):
    n_nan = int(np.isnan(Z).sum() + np.isinf(Z).sum())
    return np.nan_to_num(Z.astype(np.float32), nan=0.0, posinf=10.0, neginf=-10.0), n_nan


def load_z_disease():
    Z = np.load(config.Z_DISEASE)
    dis_names = np.load(config.Z_SAMPLE_NAMES, allow_pickle=True).tolist()
    gene_names = np.load(config.Z_GENE_NAMES, allow_pickle=True).tolist()
    Z, _ = clean_z(Z)
    return Z, dis_names, gene_names


def load_z_hc():
    """HC Z-score 캐시 로드 (disease_scoring/scoring.score_hc 산출)."""
    if not config.Z_HC.exists():
        raise FileNotFoundError(
            f'{config.Z_HC} 없음. disease_scoring 노트북 또는 scoring.score_hc() 먼저 실행.')
    Z = np.load(config.Z_HC)
    names = np.load(config.Z_HC_NAMES, allow_pickle=True).tolist()
    Z, _ = clean_z(Z)
    return Z, names


def ood_min_samples_filter(Z_dis, dis_pheno, dis_names, X_hc, X_dis,
                           percentile=None, min_samples=None, ood=None):
    """OOD(Mahalanobis, HC-fit) + MIN_SAMPLES 제외. (filtered Z, pheno, names, ood, keep, excluded) 반환."""
    percentile = MP['ood_percentile'] if percentile is None else percentile
    min_samples = MP['min_samples'] if min_samples is None else min_samples
    if ood is None:
        ood = MahalanobisFilter(percentile=percentile)
        ood.fit(X_hc)
    keep = ood.mask(X_dis)
    Z_dis = Z_dis[keep]
    dis_pheno = dis_pheno[keep]
    dis_names = [dis_names[i] for i, k in enumerate(keep) if k]
    counts = pd.Series(dis_pheno).value_counts()
    excluded = counts[counts < min_samples].index.tolist()
    if excluded:
        keep2 = np.array([p not in excluded for p in dis_pheno])
        Z_dis = Z_dis[keep2]
        dis_pheno = dis_pheno[keep2]
        dis_names = [dis_names[i] for i, k in enumerate(keep2) if k]
    return Z_dis, dis_pheno, dis_names, ood, keep, excluded


@dataclass
class DiseaseData:
    Z_dis: np.ndarray
    dis_pheno: np.ndarray
    dis_names: list
    gene_names: list
    gene_syms: list
    adata: object
    is_hc: np.ndarray
    X_raw: np.ndarray


def load_disease_filtered(with_symbols=True, adata=None):
    """selection/enrichment/signatures 공용: Z 로드 + OOD/MIN_SAMPLES 필터까지."""
    if adata is None:
        adata = load_adata()
    is_hc, phenos, _ = make_phenotypes(adata)
    X_raw = bias_matrix(adata)
    Z_dis, dis_names, gene_names = load_z_disease()
    dis_pheno = phenos[~is_hc]
    Z_dis, dis_pheno, dis_names, *_ = ood_min_samples_filter(
        Z_dis, dis_pheno, dis_names, X_raw[is_hc], X_raw[~is_hc])
    gene_syms = adata.var['GeneName'][gene_names].tolist() if with_symbols else None
    return DiseaseData(Z_dis, dis_pheno, dis_names, gene_names, gene_syms, adata, is_hc, X_raw)
