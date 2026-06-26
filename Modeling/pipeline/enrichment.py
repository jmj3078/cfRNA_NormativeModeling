import numpy as np
import pandas as pd
import gseapy as gp

import config
from gene_selectors import GeneSelector

MP = config.MODELING_PARAMS


def gsea_run_dir(ubiquity_thr=None, ubiquity_abs_z=0.5):
    """Return the GSEA output directory for a given filter condition.

    ubiquity_thr=None  → config.GSEA_DIR / 'no_filter'
    ubiquity_thr=0.7, ubiquity_abs_z=0.7 → config.GSEA_DIR / 'ubiq70_abs7'

    Use this function in both the enrichment notebook (to decide where to save)
    and in downstream notebooks (to decide where to load from), so the path is
    always derived from the same parameters without hardcoding.
    """
    if ubiquity_thr is None:
        return config.GSEA_DIR / 'no_filter'
    t = f"ubiq{int(ubiquity_thr * 100)}_abs{int(ubiquity_abs_z * 10)}"
    return config.GSEA_DIR / t


def run_gsea_prerank(Z_dis, dis_pheno, gene_syms, outdir=None,
                     min_size=10, max_size=500, save=True,
                     ubiquity_thr=None, ubiquity_abs_z=0.5):
    """Run per-phenotype mean-Z GSEA prerank; return {ph: df} and save gsea_result_{ph}.csv.

    ubiquity_thr  : passed to GeneSelector.mean_z_ranking. When set, genes
                    with cross-disease ubiquity score >= ubiquity_thr are zeroed
                    out in every phenotype's ranking before GSEA. None = disabled.
    ubiquity_abs_z: |mean_Z| threshold per phenotype for ubiquity computation.
    """
    outdir = outdir or gsea_run_dir(ubiquity_thr, ubiquity_abs_z)
    outdir.mkdir(parents=True, exist_ok=True)
    gs = GeneSelector(Z_dis, dis_pheno, gene_syms)
    if ubiquity_thr is not None:
        ubiq = gs.compute_ubiquity(abs_z_thr=ubiquity_abs_z)
        n_excluded = int((ubiq >= ubiquity_thr).sum())
        print(f"ubiquity filter: zeroing {n_excluded} genes "
              f"(|mean_Z|>{ubiquity_abs_z} in >={ubiquity_thr*100:.0f}% of phenotypes)")
    results = {}
    for ph in np.unique(dis_pheno):
        ranking = gs.mean_z_ranking(ph, ubiquity_thr=ubiquity_thr,
                                    ubiquity_abs_z=ubiquity_abs_z)
        rnk_df = (pd.DataFrame(list(ranking.items()), columns=['gene', 'score'])
                    .sort_values('score', ascending=False)
                    .reset_index(drop=True))
        try:
            res = gp.prerank(
                rnk=rnk_df, gene_sets=MP['gsea_gene_sets'], outdir=None,
                min_size=min_size, max_size=max_size,
                permutation_num=MP['gsea_perm'], seed=MP['gsea_seed'], verbose=False)
            df = res.res2d[res.res2d['FDR q-val'] < MP['gsea_fdr_thr']].copy()
            results[ph] = df
            if save:
                df.to_csv(outdir / f"gsea_result_{ph.replace('/', '_')}.csv", index=False)
            print(f"{ph:25s}: {len(df):4d} sig  "
                  f"(NES>0: {(df['NES'] > 0).sum():3d}  NES<0: {(df['NES'] < 0).sum():3d})")
        except Exception as e:
            print(f"{ph}: ERROR — {e}")
    return results
