import pickle
import re

import gseapy as gp
import numpy as np
import pandas as pd
import scanpy as sc

import MixedEffectsModeling.config as config

# Per-sample pathway ORA (run_phenotype/run_phenotype_directional/run_reoccurrence_detail and the
# Jaccard reoccurrence statistics) was retired 2026-09-23 -- see
# _legacy/PerSamplePathwayAnalysis_ORA/README.md. What remains here is the gene-level substrate
# (symbol vocabulary, pathway library, preranked GSEA) that the surviving analyses still use.

PP = config.PATHWAY_CONV_PARAMS
PCDIR = config.PATHWAY_CONV_DIR


# ENSG -> gene-symbol vocabulary is phenotype-independent (fixed by the H5AD var table), cached once
# and reused across phenotypes.
def load_symbol_vocab(gene_names):
    path = PCDIR / "symbol_vocab.pkl"
    if path.exists():
        return pickle.load(open(path, "rb"))
    sym_of = sc.read_h5ad(config.H5AD_PATH, backed="r").var["GeneName"].reindex(gene_names)
    syms_all = sym_of.values
    universe_syms = pd.unique(syms_all[sym_of.notna().values])
    sym2idx = {s: i for i, s in enumerate(universe_syms)}
    col2sym = np.array([sym2idx.get(s, -1) if pd.notna(s) else -1 for s in syms_all])
    pickle.dump((universe_syms, sym2idx, col2sym), open(path, "wb"))
    return universe_syms, sym2idx, col2sym


# pathway gene-set membership (KEGG+Reactome, housekeeping-excluded) is also phenotype-independent
def load_pathway_library():
    lib_path = PCDIR / "gene_set_matrix.pkl"
    if lib_path.exists():
        terms_raw, M_raw = pickle.load(open(lib_path, "rb"))
    else:
        libs = {}
        for lib in PP["gene_sets"]:
            libs.update(gp.get_library(lib, organism="human"))
        universe_syms, sym2idx, _ = load_symbol_vocab(None)
        N = len(universe_syms)
        terms_raw = list(libs.keys())
        M_raw = np.zeros((len(terms_raw), N), dtype=bool)
        for ti, t in enumerate(terms_raw):
            for g in libs[t]:
                j = sym2idx.get(g)
                if j is not None:
                    M_raw[ti, j] = True
        keep = M_raw.sum(axis=1) >= PP["min_pathway_size"]
        terms_raw, M_raw = [t for t, k in zip(terms_raw, keep) if k], M_raw[keep]
        pickle.dump((terms_raw, M_raw), open(lib_path, "wb"))

    ribo_idx = np.where(M_raw[terms_raw.index(PP["ribo_reference_term"])])[0]
    frac_ribo = M_raw[:, ribo_idx].sum(axis=1) / M_raw.sum(axis=1)
    kw_rx = re.compile("|".join(PP["exclude_keywords"]), re.IGNORECASE)
    keep_hk = (frac_ribo <= PP["ribo_frac_max"]) & np.array([not kw_rx.search(t) for t in terms_raw])
    terms, M = [t for t, k in zip(terms_raw, keep_hk) if k], M_raw[keep_hk]
    return terms, M


def ensg_to_symbol():
    """The backed read still opens a 7.4GB h5ad (~6s a call). Callers that reindex must .copy()
    first -- the Series is shared."""
    return sc.read_h5ad(config.H5AD_PATH, backed="r").var["GeneName"]


def gsea_prerank(rnk, terms, M, universe_syms, n_perm=1000, seed=42, min_size=5, max_size=1000,
                 threads=8):
    """Preranked GSEA on a gene-symbol -> score ranking, using the same housekeeping-filtered
    KEGG+Reactome library as load_pathway_library. Returns gseapy's res2d."""
    gene_sets = {t: [universe_syms[j] for j in np.where(M[ti])[0]] for ti, t in enumerate(terms)}
    rnk_s = pd.Series(rnk, index=universe_syms).dropna().sort_values(ascending=False)
    res = gp.prerank(rnk=rnk_s, gene_sets=gene_sets, min_size=min_size, max_size=max_size,
                     permutation_num=n_perm, seed=seed, outdir=None, no_plot=True, threads=threads)
    return res.res2d


