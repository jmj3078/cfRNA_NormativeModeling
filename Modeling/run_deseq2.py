"""
DESeq2 within-study, no covariate adjustment (Author-only design).

Design: ~Author_clean + condition  (single-author: ~condition)

Outputs:
  config.DESEQ2_RESULTS_DIR / deseq2_{phenotype}.csv     -- per-gene stats
  config.DESEQ2_GSEA_DIR    / gsea_result_{phenotype}.csv -- FDR-filtered GSEA terms
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
from scipy.sparse import issparse

import gseapy as gp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

import config
from pipeline import data_prep

MP = config.MODELING_PARAMS

config.DESEQ2_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
config.DESEQ2_GSEA_DIR.mkdir(parents=True, exist_ok=True)

dd = data_prep.load_disease_filtered()
_adata = dd.adata

if _adata.raw is not None:
    try:
        _raw_X = _adata.raw[:, _adata.var_names].X
        raw_counts_full = _raw_X.toarray() if issparse(_raw_X) else np.asarray(_raw_X)
    except Exception:
        raw_counts_full = _adata.X.toarray() if issparse(_adata.X) else np.asarray(_adata.X)
else:
    raw_counts_full = _adata.X.toarray() if issparse(_adata.X) else np.asarray(_adata.X)
raw_counts_full = np.round(raw_counts_full).astype(np.int64)

obs2idx = {n: i for i, n in enumerate(_adata.obs_names)}
gene_ids = _adata.var_names.tolist()
id2sym = dict(zip(_adata.var_names, _adata.var['GeneName']))
name2pheno = dict(zip(dd.dis_names, dd.dis_pheno))
unique_phenos = sorted(np.unique(dd.dis_pheno))


def _safe_col(s):
    return re.sub(r'[^A-Za-z0-9]+', '_', s).strip('_')


def run_deseq2(phenotype, metadata, counts_df, design):
    inference = DefaultInference(n_cpus=4)
    dds = DeseqDataSet(counts=counts_df, metadata=metadata, design=design,
                       refit_cooks=True, inference=inference)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=['condition', 'Disease', 'HC'], inference=inference)
    ds.summary()
    res = ds.results_df.copy()
    res.index.name = 'ensg'
    res = res.reset_index()
    res['symbol'] = res['ensg'].map(id2sym).fillna(res['ensg'])
    return res


deseq2_computed = {}

for phenotype in unique_phenos:
    out_path = config.DESEQ2_RESULTS_DIR / f'deseq2_{phenotype.replace("/", "_").strip()}.csv'
    if out_path.exists():
        print(f'[{phenotype}] loading existing')
        deseq2_computed[phenotype] = pd.read_csv(out_path)
        continue

    dis_names = [n for n, ph in name2pheno.items() if ph == phenotype]
    if len(dis_names) < 2:
        print(f'[{phenotype}] skip: n_dis={len(dis_names)} < 2')
        continue

    dis_obs = _adata.obs.loc[dis_names]
    authors = dis_obs['Author'].unique().tolist()
    hc_mask = dd.is_hc & _adata.obs['Author'].isin(authors)
    hc_names = _adata.obs_names[hc_mask].tolist()
    if len(hc_names) < 2:
        print(f'[{phenotype}] skip: n_hc={len(hc_names)} < 2')
        continue

    sample_names = dis_names + hc_names
    row_idx = [obs2idx[n] for n in sample_names]
    counts_arr = raw_counts_full[row_idx, :]
    counts_df = pd.DataFrame(counts_arr, index=sample_names, columns=gene_ids)
    counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= 10]

    condition = ['Disease'] * len(dis_names) + ['HC'] * len(hc_names)
    author_clean = [_safe_col(a) for a in
                    (dis_obs['Author'].tolist() + _adata.obs.loc[hc_names, 'Author'].tolist())]
    n_authors = len(set(author_clean))

    metadata = pd.DataFrame({
        'condition': pd.Categorical(condition, categories=['HC', 'Disease']),
        'Author_clean': author_clean,
    }, index=sample_names)

    design = f'~Author_clean + condition' if n_authors > 1 else '~condition'

    try:
        res = run_deseq2(phenotype, metadata, counts_df, design)
    except Exception as e:
        print(f'[{phenotype}] FAILED: {e}')
        continue

    res.to_csv(out_path, index=False)
    deseq2_computed[phenotype] = res
    n_sig = int((res['padj'].dropna() < 0.05).sum())
    print(f'[{phenotype}] n_dis={len(dis_names)} n_hc={len(hc_names)} '
          f'design={design!r}  sig={n_sig}')

print(f'\nDone: {len(deseq2_computed)} / {len(unique_phenos)} phenotypes computed')

# GSEA
print('\n--- GSEA ---')
deseq2_results = {}
for phenotype, res_df in deseq2_computed.items():
    out_csv = config.DESEQ2_GSEA_DIR / f'gsea_result_{phenotype.replace("/", "_").strip()}.csv'
    if out_csv.exists():
        print(f'[{phenotype}] loading existing GSEA')
        deseq2_results[phenotype] = pd.read_csv(out_csv)
        continue

    res_df = res_df.copy()
    res_df['stat'] = res_df['stat'].fillna(0.0)
    res_df = (res_df.assign(_a=res_df['stat'].abs())
              .sort_values('_a', ascending=False)
              .drop_duplicates('symbol')
              .drop(columns='_a'))
    rnk = (res_df[['symbol', 'stat']]
           .rename(columns={'symbol': 'gene', 'stat': 'score'})
           .sort_values('score', ascending=False)
           .reset_index(drop=True))
    try:
        res = gp.prerank(rnk=rnk, gene_sets=MP['gsea_gene_sets'], outdir=None,
                         min_size=10, max_size=500,
                         permutation_num=MP['gsea_perm'], seed=MP['gsea_seed'], verbose=False)
        sig = res.res2d[res.res2d['FDR q-val'] < MP['gsea_fdr_thr']].copy()
        sig.to_csv(out_csv, index=False)
        deseq2_results[phenotype] = sig
        sig['NES'] = pd.to_numeric(sig['NES'], errors='coerce')
        n_up = int((sig['NES'] > 0).sum())
        n_dn = int((sig['NES'] < 0).sum())
        print(f'[{phenotype}] {len(sig):4d} sig  (NES>0: {n_up:3d}  NES<0: {n_dn:3d})')
    except Exception as e:
        print(f'[{phenotype}] GSEA ERROR -- {e}')

print(f'Done: {len(deseq2_results)} phenotypes with GSEA')
