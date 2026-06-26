import numpy as np
import pandas as pd
import gseapy as gp

import config
from gene_selectors import GeneSelector

MP = config.MODELING_PARAMS


def run_gsea_prerank(Z_dis, dis_pheno, gene_syms, outdir=None,
                     min_size=10, max_size=500, save=True):
    """Run per-phenotype mean-Z GSEA prerank; return {ph: df} and save gsea_result_{ph}.csv."""
    outdir = outdir or config.GSEA_DIR
    outdir.mkdir(parents=True, exist_ok=True)
    gs = GeneSelector(Z_dis, dis_pheno, gene_syms)
    results = {}
    for ph in np.unique(dis_pheno):
        ranking = gs.mean_z_ranking(ph)
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
