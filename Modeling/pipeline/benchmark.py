import pickle
import re

import pandas as pd

import config


def deseq2_path(phenotype):
    stem = re.sub(r'[^A-Za-z0-9]+', '_', phenotype).strip('_')
    return config.DESEQ2_RESULTS_DIR / f'deseq2_{stem}.csv'


def load_deseq2(phenotype):
    """DESeq2 result for one phenotype. excluded = independent-filtering dropout (pvalue NA)."""
    df = pd.read_csv(deseq2_path(phenotype))
    df['excluded'] = df['pvalue'].isna()
    return df


def load_training_summary():
    return pd.read_csv(config.ENGINE_DIR / 'training_summary.csv')


def load_rare_ref():
    with open(config.RARE_REF, 'rb') as f:
        return pickle.load(f)


def gene_branch_table():
    """Per-gene branch + HC detection rate across all three normative sub-models"""
    train = load_training_summary().rename(columns={'gene': 'ensg'})[['ensg', 'branch', 'det_rate']]
    rare = load_rare_ref().rename(columns={'gene': 'ensg', 'det_rate_hc': 'det_rate'})
    rare['branch'] = 'rare_' + rare['category']
    return pd.concat([train, rare[['ensg', 'branch', 'det_rate']]], ignore_index=True)


def load_flagged_calls(dd, phenotype):
    """Per-gene aggregate of already-flagged (|z|>=z_flag) sample calls for one phenotype.

    Sourced from disease_scores_flagged.parquet, not Z_disease.npy: the rare-event branch
    (det_rate < low_det_thr) is scored separately by RareEventScorer and only ever written
    to this parquet -- Z_disease.npy holds a zero placeholder for those genes.
    """
    names = set(n for n, p in zip(dd.dis_names, dd.dis_pheno) if p == phenotype)
    flagged = pd.read_parquet(config.Z_SCORES_DIR / 'disease_scores_flagged.parquet')
    sub = flagged[flagged['sample'].isin(names)]
    agg = (sub.groupby('gene')['score']
              .agg(norm_n_flagged='size', norm_max_abs_z=lambda s: s.abs().max())
              .reset_index().rename(columns={'gene': 'ensg'}))
    agg['norm_flagged_prop'] = agg['norm_n_flagged'] / len(names)
    return agg, len(names)


def gene_level_compare(dd, phenotype):
    """Merge DESeq2 per-gene results with normative flagged-call stats for one phenotype's disease group."""
    deseq = load_deseq2(phenotype)
    branch = gene_branch_table()
    flagged, n_samples = load_flagged_calls(dd, phenotype)
    merged = (
        deseq[['ensg', 'symbol', 'baseMean', 'log2FoldChange', 'pvalue', 'padj', 'excluded']]
        .merge(branch, on='ensg', how='left')
        .merge(flagged, on='ensg', how='left')
    )
    fill_cols = ['norm_n_flagged', 'norm_max_abs_z', 'norm_flagged_prop']
    merged[fill_cols] = merged[fill_cols].fillna(0.0)
    merged.attrs['n_disease_samples'] = n_samples
    return merged


def rescued_genes(merged, min_n_flagged=1):
    """DESeq2-excluded (independent filtering) genes with at least one flagged disease sample."""
    sub = merged[merged['excluded'] & (merged['norm_n_flagged'] >= min_n_flagged)].copy()
    return sub.sort_values('norm_max_abs_z', ascending=False)
