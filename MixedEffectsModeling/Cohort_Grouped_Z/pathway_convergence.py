import pickle
import re

import gseapy as gp
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import norm

import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.calibration import bh_fdr_reject

# Per-sample pathway ORA (run_phenotype/run_phenotype_directional/run_reoccurrence_detail and the
# Jaccard reoccurrence statistics) was retired 2026-09-23 -- see
# _legacy/PerSamplePathwayAnalysis_ORA/README.md. What remains here is the gene-level substrate
# (symbol collapse, pathway library, per-sample BH) that the surviving analyses still use.

PP = config.PATHWAY_CONV_PARAMS
PCDIR = config.PATHWAY_CONV_DIR


def slugify(phenotype):
    return phenotype.strip().replace(" ", "_").replace("/", "-")


def gene_z_path(label):
    return PCDIR / f"{slugify(label)}_gene_z.pkl"


# ENSG -> gene-symbol vocabulary is phenotype-independent (fixed by the H5AD var table), cached once
# and reused across phenotypes. Zu/Fm (the collapsed Z matrix) is NOT -- it depends on which patients
# are in the cohort, so it's computed fresh per phenotype in collapse_to_symbols below.
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


def collapse_to_symbols(Zc, col2sym, N):
    n_pat = Zc.shape[0]
    keep_cols = col2sym >= 0
    Zc_v, sym_idx_v = Zc[:, keep_cols], col2sym[keep_cols]
    sum_mat, finite_cnt = np.zeros((n_pat, N)), np.zeros((n_pat, N))
    for i in range(n_pat):
        np.add.at(sum_mat[i], sym_idx_v, np.nan_to_num(Zc_v[i], nan=0.0))
        np.add.at(finite_cnt[i], sym_idx_v, np.isfinite(Zc_v[i]).astype(float))
    Zu_raw = np.divide(sum_mat, finite_cnt, out=np.full_like(sum_mat, np.nan), where=finite_cnt > 0)
    Zu = np.nan_to_num(Zu_raw, nan=0.0)
    Fm = np.isfinite(Zu_raw).astype(float)
    return Zu, Fm


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


def gene_sig_at_q(Zu, Fm, q):
    p_all = 2 * norm.sf(np.abs(Zu))
    sig = np.zeros_like(Fm, dtype=bool)
    for i in range(Zu.shape[0]):
        present = Fm[i] > 0
        if present.any():
            sig[i, present] = bh_fdr_reject(p_all[i, present], q=q)
    return sig


def q_tag(q):
    return "" if q == PP["fdr_q"] else f"_q{q:g}".replace(".", "")
