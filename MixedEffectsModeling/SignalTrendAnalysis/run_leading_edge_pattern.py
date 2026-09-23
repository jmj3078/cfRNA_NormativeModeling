import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import MixedEffectsModeling.config as config
from MixedEffectsModeling.PerSamplePathwayAnalysis.pathway_convergence import gene_z_path
from MixedEffectsModeling.SignalTrendAnalysis.sankey_helpers import match_pathway_index

PCDIR = config.PATHWAY_CONV_DIR
CURDIR = config.SIGNAL_TREND_CUR_DIR
GSEA_CACHE = config.ROOT / "MixedEffectsModeling" / "Benchmark" / "gsea_cache"

# curation-file basename -> PerSamplePathwayAnalysis/<slug>/ folder name (differ for multi-study phenotypes)
SLUG_MAP = {
    "Tuberculosis": "Tuberculosis",
    "Pancreatitis": "Pancreatitis",
    "Pancreatic_Cancer": "Pancreatic_Cancer",
    "Pre-eclampsia": "Pre-eclampsia",
    "Colorectal_Cancer": "Colorectal_Cancer",
    "Lung_Cancer": "Lung_Cancer",
    "Esophagus_Cancer": "Esophagus_Cancer",
    "Stomach_Cancer": "Stomach_Cancer",
    "Liver_Cancer_Roskams-Hieter": "Liver_Cancer_(Roskams-Hieter_B_et_al.)",
    "Liver_Cancer_Chen": "Liver_Cancer_(Chen_et_al.)",
}

Z_THRESH = 1.96


def parse_curation(md_path):
    text = open(md_path).read()
    src_m = re.search(r"`(normative|covariate|no_covariate)__[^`]+\.csv`", text)
    gsea_file = src_m.group(0).strip("`") if src_m else None
    sel = text.split("## Selected Pathways")[1].split("## Dropped candidates")[0]
    terms = re.findall(r"^### (.+)$", sel, re.M)
    return gsea_file, [t.strip() for t in terms]


def lead_genes_for(term, gsea_file):
    df = pd.read_csv(GSEA_CACHE / gsea_file)
    terms = df["Term"].tolist()
    j = match_pathway_index(term, terms)
    if j is None:
        return []
    return df.iloc[j]["Lead_genes"].split(";")


def gini(x):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    cum = np.cumsum(x)
    return (n + 1 - 2 * (cum / cum[-1]).sum()) / n * -1 + 1  # standard formula, positive


if __name__ == "__main__":
    rows = []
    heat_data = {}
    for md_path in sorted(CURDIR.glob("*.md")):
        if md_path.stem == "SUMMARY":
            continue
        gsea_file, terms = parse_curation(md_path)
        slug = SLUG_MAP[md_path.stem]
        pdir = PCDIR / slug
        universe_syms, sym2idx, _ = None, None, None
        d = pickle.load(open(pdir / "sig.pkl", "rb"))
        universe_syms = d["universe_syms"]
        sym2idx = {s: i for i, s in enumerate(universe_syms)}
        Zu, Fm = pickle.load(open(gene_z_path(slug), "rb"))
        n_pat = Zu.shape[0]

        for term in terms:
            genes = lead_genes_for(term, gsea_file)
            idx = [sym2idx[g] for g in genes if g in sym2idx]
            if len(idx) < 3:
                continue
            Zsub = Zu[:, idx]
            hit = np.abs(Zsub) > Z_THRESH
            gene_hit_rate = hit.mean(axis=0)
            pat_hit_count = hit.sum(axis=1)
            total_hits = hit.sum()
            top_gene_share = gene_hit_rate.max() / gene_hit_rate.sum() if gene_hit_rate.sum() > 0 else np.nan
            rows.append(dict(
                phenotype=md_path.stem, term=term, n_pat=n_pat, n_genes_matched=len(idx),
                total_hits=int(total_hits), pat_hit_median=float(np.median(pat_hit_count)),
                pat_hit_frac_any=float((pat_hit_count > 0).mean()),
                top_gene_share=float(top_gene_share), gini_gene_hits=float(gini(hit.sum(axis=0))),
            ))
            heat_data[(md_path.stem, term)] = (Zsub, [genes[i] for i, g in enumerate(genes) if g in sym2idx])

    summary = pd.DataFrame(rows).sort_values(["phenotype", "top_gene_share"], ascending=[True, False])
    summary.to_csv(CURDIR / "leading_edge_pattern_summary.csv", index=False)
    pickle.dump(heat_data, open(CURDIR / "leading_edge_pattern_data.pkl", "wb"))
    print(summary.to_string(index=False))
