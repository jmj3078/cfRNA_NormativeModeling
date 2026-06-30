import pickle
import re

import numpy as np
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


def gene_sample_detail(dd, phenotype, genes_df):
    """Per-sample raw count + Z-score for genes_df['ensg']/['branch'] within one
    phenotype's disease group, with a 'flagged' column from disease_scores_flagged.parquet.

    NBI/ZINBI/logistic genes pull Z straight from dd.Z_dis. Rare-branch genes are not in
    dd.Z_dis (RareEventScorer scores them separately -- see load_flagged_calls), so their
    Z is computed on the fly, restricted to just these genes for speed.
    """
    from pipeline import data_prep
    from pipeline.scoring import _rare_scores, load_engine

    names = [n for n, p in zip(dd.dis_names, dd.dis_pheno) if p == phenotype]
    row_of = {n: i for i, n in enumerate(dd.adata.obs_names)}
    rows = [row_of[n] for n in names]
    dis_row_of = {n: i for i, n in enumerate(dd.dis_names)}
    dis_rows = [dis_row_of[n] for n in names]
    col_of = {g: i for i, g in enumerate(dd.gene_names)}

    ensg_list = genes_df['ensg'].tolist()
    is_rare = genes_df['branch'].astype(str).str.startswith('rare')
    engine_genes = genes_df.loc[~is_rare, 'ensg'].tolist()
    rare_genes = genes_df.loc[is_rare, 'ensg'].tolist()

    Y = data_prep.count_matrix(dd.adata)[rows][:, [col_of[g] for g in ensg_list]]
    counts = pd.DataFrame(Y, index=names, columns=ensg_list)

    z_parts = []
    if engine_genes:
        cols = [col_of[g] for g in engine_genes]
        Z = dd.Z_dis[np.ix_(dis_rows, cols)]
        z_parts.append(pd.DataFrame(Z, index=names, columns=engine_genes))
    if rare_genes:
        _, rare_scorer = load_engine()
        Y_rare = counts[rare_genes].values
        sc_rare, _, r_genes, _ = _rare_scores(rare_scorer, rare_genes, Y_rare)
        z_parts.append(pd.DataFrame(sc_rare, index=names, columns=list(r_genes)))
    z = pd.concat(z_parts, axis=1)[ensg_list]

    long = (counts.rename_axis('sample').reset_index()
                  .melt(id_vars='sample', var_name='ensg', value_name='raw_count'))
    long = long.merge(
        z.rename_axis('sample').reset_index().melt(id_vars='sample', var_name='ensg', value_name='z'),
        on=['sample', 'ensg'])

    flagged = pd.read_parquet(config.Z_SCORES_DIR / 'disease_scores_flagged.parquet')
    flagged = flagged[flagged['sample'].isin(names) & flagged['gene'].isin(ensg_list)]
    flag_keys = set(zip(flagged['sample'], flagged['gene']))
    long['flagged'] = [(s, g) in flag_keys for s, g in zip(long['sample'], long['ensg'])]
    return long
