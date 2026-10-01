"""Pancreas cell-type marker panel from PanglaoDB, shared by 4_disease_scoring and
6_pancreatic_cancer_cohort_deepdive.

Lives here rather than in either notebook because 6 builds marker_panel.csv from Z, which
4 produces -- reading that file from 4 would be circular. The panel itself depends only on
the PanglaoDB library, so it is cheap to rebuild on both sides.
"""
import pickle
from collections import Counter

import gseapy as gp
import numpy as np
import pandas as pd

import MixedEffectsModeling.config as config

PANGLAO_LIB = "PanglaoDB_Augmented_2021"
PANGLAO_GROUPS = {"acinar": ["Acinar Cells"],
                  "islet": ["Alpha Cells", "Beta Cells", "Delta Cells", "Gamma (PP) Cells"],
                  "ductal": ["Ductal Cells"]}
PANGLAO_PANCREAS_OTHER = ["Pancreatic Progenitor Cells", "Pancreatic Stellate Cells",
                          "Peri-islet Schwann Cells", "Neuroendocrine Cells",
                          "Enteroendocrine Cells"]
MAX_OTHER_SETS = 1


def load_library():
    path = config.GROUP_VS_INDIV_DIR / f"{PANGLAO_LIB}.pkl"
    if path.exists():
        return pickle.load(open(path, "rb"))
    lib = gp.get_library(PANGLAO_LIB, organism="human")
    path.parent.mkdir(parents=True, exist_ok=True)
    pickle.dump(lib, open(path, "wb"))
    return lib


def marker_sets(lib=None):
    """{set: [symbol]} keeping only symbols in <= MAX_OTHER_SETS non-pancreatic cell types,
    and dropping symbols that land in more than one pancreatic set."""
    lib = load_library() if lib is None else lib
    panc = set(sum(PANGLAO_GROUPS.values(), [])) | set(PANGLAO_PANCREAS_OTHER)
    other = Counter()
    for term, genes in lib.items():
        if term not in panc:
            other.update(set(genes))
    grp_of = {}
    for lab, terms in PANGLAO_GROUPS.items():
        for sym in set().union(*[set(lib[t]) for t in terms]):
            if other[sym] <= MAX_OTHER_SETS:
                grp_of.setdefault(sym, set()).add(lab)
    return {lab: sorted(s for s, g in grp_of.items() if g == {lab}) for lab in PANGLAO_GROUPS}


def panel_table(gene_names, sym_of):
    """One row per engine-modeled gene on the panel.
    gene_names: ENSG ids in Z column order. sym_of: Series of symbols aligned to them."""
    sets = marker_sets()
    want = {sym: lab for lab, syms in sets.items() for sym in syms}
    rows = [dict(ens=e, col=j, set=want[s], symbol=s)
            for j, (e, s) in enumerate(zip(gene_names, np.asarray(sym_of))) if s in want]
    return pd.DataFrame(rows)
