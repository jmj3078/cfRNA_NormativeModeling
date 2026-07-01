"""
DESeq2 within-study with covariate adjustment (BIAS_COLUMNS x10, z-standardized).

Design: ~Author_clean + cov1 + ... + cov10 + condition
When a phenotype has too few samples to fit the full design (singular matrix or n < rank),
falls back to Author-only (~Author_clean + condition) and records the fallback.

Outputs:
  config.DESEQ2_COV_RESULTS_DIR / deseq2_{phenotype}.csv  -- per-gene stats
  config.DESEQ2_COV_GSEA_DIR    / gsea_result_{phenotype}.csv  -- FDR-filtered GSEA terms
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler

import gseapy as gp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

import config
from pipeline import data_prep

MP = config.MODELING_PARAMS
BIAS = config.BIAS_COLUMNS
COV_COLS = [re.sub(r'[^A-Za-z0-9]+', '_', b).strip('_') for b in BIAS]

config.DESEQ2_COV_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
config.DESEQ2_COV_GSEA_DIR.mkdir(parents=True, exist_ok=True)

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

# fit a single scaler on ALL samples for the 10 covariates (same as normative model convention)
scaler = StandardScaler()
scaler.fit(_adata.obs[BIAS].values)


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
fallback_log = {}

for phenotype in unique_phenos:
    out_path = config.DESEQ2_COV_RESULTS_DIR / f'deseq2_{phenotype.replace("/", "_").strip()}.csv'
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

    raw_cov = _adata.obs.loc[sample_names, BIAS].values.astype(float)
    scaled_cov = scaler.transform(raw_cov)
    cov_df = pd.DataFrame(scaled_cov, index=sample_names, columns=COV_COLS)

    condition = ['Disease'] * len(dis_names) + ['HC'] * len(hc_names)
    author_clean = [_safe_col(a) for a in
                    (dis_obs['Author'].tolist() + _adata.obs.loc[hc_names, 'Author'].tolist())]
    n_authors = len(set(author_clean))

    metadata = pd.DataFrame({
        'condition': pd.Categorical(condition, categories=['HC', 'Disease']),
        'Author_clean': author_clean,
    }, index=sample_names)
    metadata = pd.concat([metadata, cov_df], axis=1)

    author_term = 'Author_clean + ' if n_authors > 1 else ''
    full_design = f'~{author_term}' + ' + '.join(COV_COLS) + ' + condition'

    n_total = len(sample_names)
    n_author_levels = len(set(author_clean)) - 1 if n_authors > 1 else 0
    n_terms = n_author_levels + len(COV_COLS) + 1  # +1 condition

    used_design = full_design
    try:
        if n_total <= n_terms + 2:
            raise ValueError(f'n={n_total} <= rank={n_terms+1}: insufficient df')
        res = run_deseq2(phenotype, metadata, counts_df, full_design)
        fallback_log[phenotype] = 'full'
    except Exception as e_full:
        # fallback: Author-only (mirrors original deseq2_benchmark behaviour)
        author_only = f'~{author_term}condition' if n_authors > 1 else '~condition'
        print(f'[{phenotype}] full covariate design failed ({e_full}); '
              f'falling back to {author_only!r}')
        try:
            meta_fallback = metadata[['condition', 'Author_clean']]
            res = run_deseq2(phenotype, meta_fallback, counts_df, author_only)
            used_design = author_only
            fallback_log[phenotype] = 'author_only'
        except Exception as e_fb:
            print(f'[{phenotype}] fallback also failed: {e_fb}')
            continue

    res.to_csv(out_path, index=False)
    deseq2_computed[phenotype] = res
    n_sig = int((res['padj'].dropna() < 0.05).sum())
    print(f'[{phenotype}] n_dis={len(dis_names)} n_hc={len(hc_names)} '
          f'design={used_design!r}  sig={n_sig}')

print(f'\nDone: {len(deseq2_computed)} / {len(unique_phenos)} phenotypes computed')
if fallback_log:
    for ph, status in fallback_log.items():
        if status != 'full':
            print(f'  fallback -> {ph}: {status}')

# GSEA
print('\n--- GSEA ---')
deseq2_results = {}
for phenotype, res_df in deseq2_computed.items():
    out_csv = config.DESEQ2_COV_GSEA_DIR / f'gsea_result_{phenotype.replace("/", "_").strip()}.csv'
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
