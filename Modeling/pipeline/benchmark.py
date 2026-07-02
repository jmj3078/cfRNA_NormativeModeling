import numpy as np
import pandas as pd
import config


def _pheno_stem(phenotype):
    return phenotype.replace('/', '_').strip()


def deseq2_path(phenotype, cov=False):
    d = config.DESEQ2_COV_RESULTS_DIR if cov else config.DESEQ2_RESULTS_DIR
    return d / f'deseq2_{_pheno_stem(phenotype)}.csv'


def sig_gene_set(phenotype, cov=False, padj_thr=0.05):
    """Set of significant ENSG IDs from DESeq2 (no-cov or w/cov)."""
    p = deseq2_path(phenotype, cov=cov)
    if not p.exists():
        return set()
    df = pd.read_csv(p)
    return set(df.loc[df['padj'].fillna(1) < padj_thr, 'ensg'].dropna())


def norm_sig_gene_set(dd, phenotype):
    """Set of ENSG IDs flagged (|z|>=z_flag) in any disease sample of this phenotype."""
    names = {n for n, p in zip(dd.dis_names, dd.dis_pheno) if p == phenotype}
    flagged = pd.read_parquet(config.Z_SCORES_DIR / 'disease_scores_flagged.parquet')
    sub = flagged[flagged['sample'].isin(names)]
    return set(sub['gene'].unique())


def build_venn_sets(dd, phenotype, padj_thr=0.05):
    """Dict of gene sets (ENSG) for 3-way Venn: DESeq2 no-cov, DESeq2 w/cov, Normative."""
    return {
        'DESeq2\n(no-cov)': sig_gene_set(phenotype, cov=False, padj_thr=padj_thr),
        'DESeq2\n(w/cov)': sig_gene_set(phenotype, cov=True, padj_thr=padj_thr),
        'Normative\nModel': norm_sig_gene_set(dd, phenotype),
    }
