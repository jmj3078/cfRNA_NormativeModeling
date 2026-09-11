"""Symmetric group-vs-normative DB-hit comparison (gene-level + pathway-level).

Every method is scored by the same rule -- own significant genes/pathways intersected
with an Open Targets reference -- so counts (recall) and rates (precision) are directly
comparable. See MixedEffectsModeling/CLAUDE.md and PerSamplePathwayAnalysis for the underlying
normative sample-level results this reuses.
"""
import functools
import pickle
import re

import numpy as np
import pandas as pd
import scanpy as sc

import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.calibration import bh_fdr_reject
from MixedEffectsModeling.PerSamplePathwayAnalysis.pathway_convergence import (
    collapse_to_symbols, load_pathway_library, load_symbol_vocab, slugify,
)

HERE = config.ROOT / "MixedEffectsModeling" / "Benchmark"
DESEQ_DIR = HERE / "DESeq2"
REF_DIR = HERE / "disease_reference"
PC_DIR = config.PATHWAY_CONV_DIR
NULL_DIR = HERE / "DESeq2_pathway_null"
GSEA_DIR = HERE / "gsea_cache"
ZDIR = config.ROOT / "MixedEffectsModeling" / "Z_scores_mixed"
DESIGN_LABELS = {"no_covariate": "deseq2_no_cov", "covariate": "deseq2_cov",
                 "ruvg_k1": "deseq2_ruvg_k1", "ruvg_k2": "deseq2_ruvg_k2", "ruvg_k3": "deseq2_ruvg_k3"}
# same reference definition as 4_disease_scoring.ipynb cell 10 (marker panel): top-N by OT score,
# floored -- keeps reference size comparable across phenotypes instead of scaling with literature
# volume (an unfloored/unranked reference balloons to ~80% of the genome for well-studied cancers)
MARKER_TOPN = 300
MARKER_SCORE_FLOOR = 0.05


def load_reference(topn=MARKER_TOPN, score_floor=MARKER_SCORE_FLOOR):
    """{phenotype: set(symbols)} from Open Targets, top-N by score above a floor."""
    out = {}
    for f in REF_DIR.glob("*.json"):
        import json
        r = json.load(open(f))
        if r["phenotype"] == "MGUS":  # OT association pool too thin (22 genes) for any panel size
            continue
        ranked = sorted(r["genes"], key=lambda gs: -gs[1])
        floored = [g for g, s in ranked if s >= score_floor]
        out[r["phenotype"]] = set(floored[:topn]) | set(r.get("supplement", {}).get("genes", []))
    return out


@functools.lru_cache(maxsize=1)
def ensg_to_symbol():
    """Memoized -- the backed read still opens a 7.4GB h5ad (~6s a call) and this is hit from a
    dozen call sites. Callers that reindex must .copy() first (the Series is shared)."""
    return sc.read_h5ad(config.H5AD_PATH, backed="r").var["GeneName"]


def deseq2_tag(study, pheno):
    return f"{study.replace(' ', '_')}__{pheno.replace('/', '-').replace(' ', '_')}"


@functools.lru_cache(maxsize=None)
def deseq2_study_results(design):
    """{(study, phenotype): results_df} for one design, indexed by base ENSG (no version).
    Memoized -- 18 x ~20k-row CSVs per design, and the q-sweep in group_level_z_test asks for the
    same design once per q. Callers must treat the returned frames as read-only (shared)."""
    summary = pd.read_csv(DESEQ_DIR / "summary.csv")
    summary = summary[summary.design == design]
    out = {}
    for _, row in summary.iterrows():
        path = DESEQ_DIR / design / f"{deseq2_tag(row.study, row.phenotype)}.csv"
        if not path.exists():
            continue
        res = pd.read_csv(path, index_col=0)
        res.index = res.index.str.split(".").str[0]
        out[(row.study, row.phenotype)] = res
    return out


def deseq2_gene_sets(design, sym_of, alpha=0.05):
    """{phenotype: set(symbols)} pooling padj<alpha genes across studies (union)."""
    sym_of = sym_of.copy()
    sym_of.index = sym_of.index.str.split(".").str[0]
    out = {}
    for (study, pheno), res in deseq2_study_results(design).items():
        sig = res.index[res["padj"] < alpha]
        syms = set(sym_of.reindex(sig).dropna())
        out.setdefault(pheno, set()).update(syms)
    return out


def pc_dirs_by_phenotype(sample_meta):
    """{phenotype: [pdir, ...]} -- resolves PerSamplePathwayAnalysis subdirs (some phenotypes are
    split per-study) back to the canonical phenotype via the samples actually inside sig.pkl,
    not by parsing directory names."""
    ph_of = sample_meta.set_index("sample")["phenotype"]
    out = {}
    for pdir in sorted(d for d in PC_DIR.iterdir() if d.is_dir() and (d / "sig.pkl").exists()):
        d = pickle.load(open(pdir / "sig.pkl", "rb"))
        phenos = ph_of.reindex(d["names_c"]).dropna().unique()
        if len(phenos) != 1:
            raise ValueError(f"{pdir} mixes phenotypes: {phenos}")
        out.setdefault(phenos[0], []).append(pdir)
    return out


def model_route_mask(gene_names):
    """Bool array over gene_names, True for individually-fitted ('model' route) genes --
    False excludes the ~2060 pooled-GLM ('rare') genes, for the no-rare-pooling comparison."""
    ts = pd.read_csv(config.ENGINE_MIXED_DIR / "training_summary.csv").set_index("gene")["route"]
    return (ts.reindex(gene_names) == "model").values


def normative_gene_hits(sample_meta, gene_names, sym_of, q=0.05, Z=None, exclude_pool=False):
    """{phenotype: [per-patient significant-symbol-set, ...]}, computed directly from
    Z_disease_shash.npy (SHASH-calibrated, same array PerSamplePathwayAnalysis/4_disease_scoring.ipynb
    score significance from) for every phenotype with OOD-kept samples -- NOT sourced from
    PerSamplePathwayAnalysis's sig.pkl, which only covers the ~9 phenotypes that batch script scoped in
    (missing HIV, ME/CFS, MM, MGUS, Liver Cirrhosis, HIV+TB, CAD_HF+/-, Other Cancer, ICI-*).
    Raw per-patient sets -- callers build union (K=1) or recurrence (K>=k) gene sets from these
    via gene_sets_from_hits without recomputing p-values per threshold.
    exclude_pool=True restricts to individually-fitted genes (drops the pooled-GLM rare route)."""
    from scipy.stats import norm
    if Z is None:
        Z = np.load(ZDIR / "Z_disease_shash.npy")
    p_all = 2 * norm.sf(np.abs(Z))
    sym_arr = sym_of.reindex(gene_names).values
    col_mask = model_route_mask(gene_names) if exclude_pool else None
    out = {}
    for pheno, sub in sample_meta[sample_meta.ood_keep].groupby("phenotype"):
        per_patient = []
        for i in sub.index:
            row = p_all[i]
            finite = np.isfinite(row)
            reject = np.zeros(len(row), dtype=bool)
            reject[finite] = bh_fdr_reject(row[finite], q=q)
            if col_mask is not None:
                reject = reject & col_mask
            per_patient.append({s for s in sym_arr[reject] if pd.notna(s)})
        out[pheno] = per_patient
    return out


def gene_sets_from_hits(per_patient_hits, K=1):
    """{phenotype: symbol_set} -- a symbol counts if significant in >=K patients. K=1 is the plain
    union; K>=3 is the recurrence variant (same convention as normative_pathway_recurrence)."""
    out = {}
    for pheno, per_patient in per_patient_hits.items():
        counts = {}
        for syms in per_patient:
            for s in syms:
                counts[s] = counts.get(s, 0) + 1
        out[pheno] = {s for s, c in counts.items() if c >= K}
    return out


def gene_venn_sets(design_b="covariate"):
    """{phenotype: (deseq2_no_cov, deseq2_<design_b>, normative_union)} symbol sets, phenotypes
    present in all three methods only. design_b='ruvg_k1' swaps the second circle to
    RUVg-corrected DESeq2 instead of explicit-covariate DESeq2."""
    sym_of = ensg_to_symbol()
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    nocov = deseq2_gene_sets("no_covariate", sym_of)
    db = deseq2_gene_sets(design_b, sym_of)
    norm_union = gene_sets_from_hits(normative_gene_hits(sm, gene_names, sym_of), K=1)
    return {ph: (nocov[ph], db[ph], norm_union[ph])
            for ph in sorted(set(nocov) & set(db) & set(norm_union))}


def gene_venn_sets_4way(design_c="ruvg_k2"):
    """{phenotype: {deseq2_no_cov, deseq2_covariate, deseq2_<design_c>, normative}} symbol sets --
    same 4 methods as 6_group_level_comparison.ipynb's DESEQ2_DESIGNS + normative_union, for
    checking whether normative overlaps DESeq2(no_cov) more than the covariate-adjusted designs
    (a sign covariate adjustment isn't actually being absorbed in the normative Z)."""
    sym_of = ensg_to_symbol()
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    nocov = deseq2_gene_sets("no_covariate", sym_of)
    cov = deseq2_gene_sets("covariate", sym_of)
    ruvg = deseq2_gene_sets(design_c, sym_of)
    norm_union = gene_sets_from_hits(normative_gene_hits(sm, gene_names, sym_of), K=1)
    common = sorted(set(nocov) & set(cov) & set(ruvg) & set(norm_union))
    return {ph: {"deseq2_no_cov": nocov[ph], "deseq2_covariate": cov[ph],
                 DESIGN_LABELS[design_c]: ruvg[ph], "normative": norm_union[ph]} for ph in common}


def novel_gene_table():
    """One row per (phenotype, gene) for genes normative flags (union, K=1) that neither DESeq2
    design (no_covariate/covariate) reaches, restricted to genes with an Open Targets association
    score -- i.e. disease-relevant per DB, missed by the group-wise comparison. Also checked
    against Benchmark/literature_biomarkers.json (the original source paper's own named
    biomarker panel, where available) via `in_literature_panel` -- rows where this is False
    (and `lit_panel_available` is True) are missed by BOTH DESeq2 and the source paper's own
    biomarker list, not just DESeq2."""
    import json
    sym_of = ensg_to_symbol()
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    nocov = deseq2_gene_sets("no_covariate", sym_of)
    cov = deseq2_gene_sets("covariate", sym_of)
    hits = normative_gene_hits(sm, gene_names, sym_of)
    norm_union = gene_sets_from_hits(hits, K=1)

    counts_by_pheno = {}
    for pheno, per_patient in hits.items():
        c = {}
        for syms in per_patient:
            for s in syms:
                c[s] = c.get(s, 0) + 1
        counts_by_pheno[pheno] = c

    ref = load_reference()  # same top-300/floor-0.05 marker panel as db_hit_row elsewhere
    ref_scores = {}
    for f in REF_DIR.glob("*.json"):
        r = json.load(open(f))
        ref_scores[r["phenotype"]] = dict(r["genes"])
    lit = json.load(open(HERE / "literature_biomarkers.json"))

    rows = []
    for pheno in sorted(set(nocov) & set(cov) & set(norm_union)):
        novel = norm_union[pheno] - nocov[pheno] - cov[pheno]
        scores = ref_scores.get(pheno, {})
        n_total = len(hits[pheno])
        lit_entry = lit.get(pheno, {})
        lit_available = len(lit_entry.get("genes", [])) > 0
        for g in novel & ref.get(pheno, set()):
            if g in scores:
                rows.append(dict(phenotype=pheno, gene=g, ot_score=scores[g],
                                  n_patients=counts_by_pheno[pheno].get(g, 0), n_total=n_total,
                                  lit_panel_available=lit_available,
                                  in_literature_panel=g in lit_entry.get("genes", [])))
    df = pd.DataFrame(rows)
    df["recur_pct"] = (100 * df.n_patients / df.n_total).round(1)
    return df.sort_values(["phenotype", "ot_score"], ascending=[True, False]).reset_index(drop=True)


def literature_overlap_table():
    """One row per (phenotype, literature-reported biomarker gene): whether normative and/or
    DESeq2 (either design) recover it. Source: Benchmark/literature_biomarkers.json, curated
    per-study from the original cfRNA papers -- see that file's 'confidence'/'notes' fields,
    several panels are low-confidence (automated PDF extraction) or explicitly unavailable."""
    import json
    lit = json.load(open(HERE / "literature_biomarkers.json"))
    sym_of = ensg_to_symbol()
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    nocov = deseq2_gene_sets("no_covariate", sym_of)
    cov = deseq2_gene_sets("covariate", sym_of)
    norm_union = gene_sets_from_hits(normative_gene_hits(sm, gene_names, sym_of), K=1)

    rows = []
    for pheno, entry in lit.items():
        genes = entry["genes"]
        for g in genes:
            rows.append(dict(
                phenotype=pheno, gene=g, confidence=entry["confidence"],
                normative_hit=g in norm_union.get(pheno, set()),
                deseq2_hit=g in nocov.get(pheno, set()) or g in cov.get(pheno, set()),
            ))
        if not genes:
            rows.append(dict(phenotype=pheno, gene=None, confidence=entry["confidence"],
                              normative_hit=None, deseq2_hit=None))
    return pd.DataFrame(rows)


def discovery_summary():
    """Per-phenotype: of all normative-significant genes (union, K=1, regardless of DESeq2
    status), how many are already named in the source paper's own biomarker panel
    (literature_biomarkers.json) vs newly discovered, and of the newly-discovered set how many
    carry Open Targets disease-association support (DB-validated). Answers 'how many
    DB-validated signals are new, not just literature reproduced' directly, without requiring
    the gene to also be missed by DESeq2 (see novel_gene_table for that stricter cut)."""
    import json
    sym_of = ensg_to_symbol()
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    norm_union = gene_sets_from_hits(normative_gene_hits(sm, gene_names, sym_of), K=1)
    ref = load_reference()
    lit = json.load(open(HERE / "literature_biomarkers.json"))

    rows = []
    for pheno, sig in norm_union.items():
        lit_entry = lit.get(pheno, {})
        lit_genes = set(lit_entry.get("genes", []))
        known = sig & lit_genes
        new = sig - lit_genes
        new_validated = new & ref.get(pheno, set())
        rows.append(dict(
            phenotype=pheno, n_sig=len(sig), lit_panel_available=len(lit_genes) > 0,
            n_literature_known=len(known), n_newly_discovered=len(new),
            n_newly_discovered_db_validated=len(new_validated),
            pct_new_db_validated=round(100 * len(new_validated) / len(new), 1) if new else np.nan,
        ))
    return pd.DataFrame(rows).sort_values("phenotype").reset_index(drop=True)


def discovery_gene_table():
    """Gene-level detail behind discovery_summary: one row per (phenotype, gene) for normative
    -significant genes not in the source paper's own biomarker panel, restricted to genes with
    Open Targets disease-association support -- the newly-discovered, DB-validated signal set."""
    import json
    sym_of = ensg_to_symbol()
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    nocov = deseq2_gene_sets("no_covariate", sym_of)
    cov = deseq2_gene_sets("covariate", sym_of)
    hits = normative_gene_hits(sm, gene_names, sym_of)
    norm_union = gene_sets_from_hits(hits, K=1)
    ref = load_reference()
    ref_scores = {}
    for f in REF_DIR.glob("*.json"):
        r = json.load(open(f))
        ref_scores[r["phenotype"]] = dict(r["genes"])
    lit = json.load(open(HERE / "literature_biomarkers.json"))

    counts_by_pheno = {}
    for pheno, per_patient in hits.items():
        c = {}
        for syms in per_patient:
            for s in syms:
                c[s] = c.get(s, 0) + 1
        counts_by_pheno[pheno] = c

    rows = []
    for pheno, sig in norm_union.items():
        lit_genes = set(lit.get(pheno, {}).get("genes", []))
        scores = ref_scores.get(pheno, {})
        n_total = len(hits[pheno])
        for g in (sig - lit_genes) & ref.get(pheno, set()):
            rows.append(dict(
                phenotype=pheno, gene=g, ot_score=scores.get(g),
                n_patients=counts_by_pheno[pheno].get(g, 0), n_total=n_total,
                deseq2_hit=g in nocov.get(pheno, set()) or g in cov.get(pheno, set()),
            ))
    df = pd.DataFrame(rows)
    df["recur_pct"] = (100 * df.n_patients / df.n_total).round(1)
    return df.sort_values(["phenotype", "ot_score"], ascending=[True, False]).reset_index(drop=True)


def final_discovery_panel(deseq2_designs=("no_covariate", "covariate")):
    """One row per phenotype: total normative-significant genes (union, K=1), how many are
    newly discovered (not in the source paper's own literature panel), how many of those carry
    Open Targets disease-association support (DB-validated), and how many of THOSE are also
    missed by every design in `deseq2_designs` (union'd -- a gene counts as a DESeq2 hit if ANY
    listed design reaches it) -- with the gene list. Default reproduces the original
    no_covariate/covariate cut; pass deseq2_designs=("no_covariate",) for the pre-RUVg
    'missing' criterion or deseq2_designs=("ruvg_k1",) for post-RUVg. Phenotypes with no Open
    Targets reference at all are dropped (nothing to validate against), not zero-filled."""
    import json
    sym_of = ensg_to_symbol()
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    deseq2_sets = [deseq2_gene_sets(d, sym_of) for d in deseq2_designs]
    norm_union = gene_sets_from_hits(normative_gene_hits(sm, gene_names, sym_of), K=1)
    ref = load_reference()
    lit = json.load(open(HERE / "literature_biomarkers.json"))

    rows = []
    for pheno, sig in norm_union.items():
        dref = ref.get(pheno, set())
        if not dref:
            continue
        lit_genes = set(lit.get(pheno, {}).get("genes", []))
        new = sig - lit_genes
        new_validated = new & dref
        deseq2_sig = set()
        for ds in deseq2_sets:
            deseq2_sig |= ds.get(pheno, set())
        new_validated_deseq2_missed = sorted(new_validated - deseq2_sig)
        rows.append(dict(
            phenotype=pheno,
            n_sig=len(sig),
            n_newly_discovered=len(new),
            n_new_db_validated=len(new_validated),
            n_new_db_validated_deseq2_missed=len(new_validated_deseq2_missed),
            genes_new_db_validated_deseq2_missed=", ".join(new_validated_deseq2_missed),
        ))
    return pd.DataFrame(rows).sort_values("phenotype").reset_index(drop=True)


def matched_threshold_sweep(qs=(0.05, 0.10, 0.15, 0.20), deseq2_design="no_covariate"):
    """Symmetric sweep: normative_union and DESeq2 (same design) at IDENTICAL q, both scored
    against the same Open Targets reference by hypergeometric enrichment + BH-FDR across the
    whole (method, q, phenotype) table -- tests whether the normative-vs-DESeq2 enrichment gap
    is a q=0.05-specific artifact or holds at every threshold."""
    import json
    from scipy.stats import hypergeom
    from statsmodels.stats.multitest import multipletests
    sym_of = ensg_to_symbol()
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    sym_arr = sym_of.reindex(gene_names).values
    ref = load_reference()
    N_genes = len(gene_names)
    Z = np.load(ZDIR / "Z_disease_shash.npy")
    from scipy.stats import norm
    p_all = 2 * norm.sf(np.abs(Z))

    rows = []
    for q in qs:
        norm_sets = {}
        for pheno, sub in sm[sm.ood_keep].groupby("phenotype"):
            sig_syms = set()
            for i in sub.index:
                row = p_all[i]
                finite = np.isfinite(row)
                reject = np.zeros(len(row), dtype=bool)
                reject[finite] = bh_fdr_reject(row[finite], q=q)
                sig_syms |= {s for s in sym_arr[reject] if pd.notna(s)}
            norm_sets[pheno] = sig_syms
        for pheno, sig in norm_sets.items():
            dref = ref.get(pheno, set())
            if not dref:
                continue
            K, n, x = len(dref), len(sig), len(sig & dref)
            pval = hypergeom.sf(x - 1, N_genes, K, n) if n > 0 else np.nan
            rows.append(dict(method="normative_union", q=q, phenotype=pheno, n_sig=n, overlap=x, ref_size=K, pval=pval))

        ds = deseq2_gene_sets(deseq2_design, sym_of, alpha=q)
        for pheno, sig in ds.items():
            dref = ref.get(pheno, set())
            if not dref:
                continue
            K, n, x = len(dref), len(sig), len(sig & dref)
            pval = hypergeom.sf(x - 1, N_genes, K, n) if n > 0 else np.nan
            rows.append(dict(method=f"deseq2_{deseq2_design}", q=q, phenotype=pheno, n_sig=n, overlap=x, ref_size=K, pval=pval))

    df = pd.DataFrame(rows)
    _, padj, _, _ = multipletests(df.pval.fillna(1), method="fdr_bh")
    df["padj"] = padj
    df["sig_fdr05"] = padj < 0.05
    return df


def cross_study_replication(recur_K=1):
    """DB-independent reproducibility check: for phenotypes split across >=2 independent
    PerSamplePathwayAnalysis study dirs (currently only Liver Cancer: Chen et al. vs Roskams-Hieter B
    et al.), hypergeometric-test each pair's significant gene/pathway sets against EACH OTHER
    (not against the Open Targets reference) -- tests whether two independent cohorts of the same
    phenotype converge on the same signal, at gene and pathway level."""
    from itertools import combinations

    from scipy.stats import hypergeom
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    pc_map = pc_dirs_by_phenotype(sm)
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    N_genes = len(gene_names)

    rows = []
    for pheno, pdirs in pc_map.items():
        if len(pdirs) < 2:
            continue
        reps = []
        for pdir in pdirs:
            d = pickle.load(open(pdir / "sig.pkl", "rb"))
            gene_idx = set(np.where(d["gene_sig"].sum(axis=0) >= recur_K)[0].tolist())
            path_idx = set(np.where(d["path_sig"].sum(axis=0) >= recur_K)[0].tolist())
            reps.append(dict(study=pdir.name, n_terms=len(d["terms"]), n_pat=d["gene_sig"].shape[0],
                              gene_idx=gene_idx, path_idx=path_idx))
        for a, b in combinations(reps, 2):
            for level, N, key in (("gene", N_genes, "gene_idx"), ("pathway", a["n_terms"], "path_idx")):
                sa, sb = a[key], b[key]
                K, n, x = len(sa), len(sb), len(sa & sb)
                pval = hypergeom.sf(x - 1, N, K, n) if (K and n) else np.nan
                rows.append(dict(phenotype=pheno, level=level, study_a=a["study"], study_b=b["study"],
                                  n_pat_a=a["n_pat"], n_pat_b=b["n_pat"], n_a=K, n_b=n, overlap=x,
                                  N_universe=N, pval=pval))
    return pd.DataFrame(rows)


def stouffer_group_z(sm=None, Z=None, gene_names=None, min_n=3):
    """{(phenotype, study): (n_genes,) array} -- Stouffer's Z (sum(Z)/sqrt(n_finite)) per gene,
    combining per-patient calibrated Z-scores into ONE group-level statistic per gene, analogous to
    DESeq2's per-gene Wald stat. Valid because SHASH calibration targets Z~N(0,1) under the null for
    every patient independently (Stouffer 1949 meta-analytic combination of independent Z-scores).
    Grouped by (phenotype, study) -- NOT phenotype alone -- because Liver Cancer pools 3 technically
    distinct studies (Chen/Roskams-Hieter B/Block, same confound PerSamplePathwayAnalysis splits out via
    per-study sig.pkl dirs, see run_pathway_convergence_batch.py); pooling them would mix each
    study's own technical variance into a single ranking statistic. min_n=3 drops Block et al.
    (n=2), too small for a group statistic (same threshold as PerSamplePathwayAnalysis)."""
    if Z is None:
        Z = np.load(ZDIR / "Z_disease_shash.npy")
    if sm is None:
        sm = pd.read_csv(ZDIR / "sample_meta.csv")
    sm = sm.copy()
    sm["study"] = sm["batch"].str.replace(r"_Batch_\d+$", "", regex=True)
    out = {}
    for (pheno, study), sub in sm[sm.ood_keep].groupby(["phenotype", "study"]):
        if len(sub) < min_n:
            continue
        Zc = Z[sub.index.values]
        finite = np.isfinite(Zc)
        n = finite.sum(axis=0)
        s = np.where(finite, Zc, 0.0).sum(axis=0)
        gz = np.divide(s, np.sqrt(n), out=np.full(n.shape, np.nan), where=n > 0)
        out[(pheno, study)] = gz
    return out


def extreme_mwu_group_z(sm=None, Z=None, gene_names=None, min_n=3, z_thresh=2.0):
    """{(phenotype, study): (n_genes,) array} -- alternative to stouffer_group_z, Rutherford 2023
    eLife design: per-gene extreme-deviation indicator (|Z|>z_thresh) compared between the group's
    patients and the HC group via Mann-Whitney U, converted to a signed Z-equivalent
    (norm.isf(p/2), signed by whether the group has a HIGHER extreme-deviation rate than HC) so it
    is a drop-in replacement for stouffer_group_z's output in group_level_z_test /
    group_level_pathway_gsea (same downstream BH-FDR / hypergeometric / GSEA pipeline).
    CAVEAT: Z_hc_shash.npy is HC scored in-sample against its own training fit (see
    NormativeModelEngineMixed.fit_shash docstring, 'deliberately in-sample') -- not a held-out/LOBO
    Z, so the HC extreme-rate baseline here is optimistically low. Fine for a relative-ranking
    comparison against stouffer_group_z; treat absolute significance with that caveat in mind."""
    from scipy.stats import mannwhitneyu, norm
    if Z is None:
        Z = np.load(ZDIR / "Z_disease_shash.npy")
    if sm is None:
        sm = pd.read_csv(ZDIR / "sample_meta.csv")
    sm = sm.copy()
    sm["study"] = sm["batch"].str.replace(r"_Batch_\d+$", "", regex=True)
    Z_hc = np.load(ZDIR / "Z_hc_shash.npy")
    ind_all = np.where(np.isfinite(Z), (np.abs(Z) > z_thresh).astype(float), np.nan)
    ind_hc = np.where(np.isfinite(Z_hc), (np.abs(Z_hc) > z_thresh).astype(float), np.nan)
    hc_rate = np.nanmean(ind_hc, axis=0)

    out = {}
    for (pheno, study), sub in sm[sm.ood_keep].groupby(["phenotype", "study"]):
        if len(sub) < min_n:
            continue
        ind_group = ind_all[sub.index.values]
        res = mannwhitneyu(ind_group, ind_hc, axis=0, nan_policy="omit", alternative="two-sided")
        sign = np.sign(np.nanmean(ind_group, axis=0) - hc_rate)
        out[(pheno, study)] = sign * norm.isf(res.pvalue / 2)
    return out


def wolfers_chi2_group_z(sm=None, Z=None, gene_names=None, min_n=3, z_thresh=2.6):
    """{(phenotype, study): (n_genes,) array} -- Wolfers et al. 2018 JAMA Psychiatry design
    (doi:10.1001/jamapsychiatry.2018.2467): binarize deviations at |Z|>2.6 (their P<.005 cut),
    then test diagnosis-vs-extreme association with a chi-square test on the per-gene 2x2
    (group x extreme) table. chi-square with 1 df is Z^2, so returning sign*sqrt(chi2) reproduces
    exactly the same two-sided p downstream and makes this a drop-in for stouffer_group_z.
    Sign = whether the group's extreme rate exceeds HC's.
    Two deliberate deviations from the paper, both for cross-method comparability:
      - Wolfers corrected with Bonferroni-Holm; this pipeline applies its own shared BH-FDR
        downstream (group_level_z_test) so every method is thresholded identically.
      - Wolfers' unit was a voxel-proportion per subject; the unit here is one gene, matching
        DESeq2's per-gene statistic.
    Shares extreme_mwu_group_z's caveat: Z_hc_shash.npy is in-sample HC (fit_shash docstring)."""
    if Z is None:
        Z = np.load(ZDIR / "Z_disease_shash.npy")
    if sm is None:
        sm = pd.read_csv(ZDIR / "sample_meta.csv")
    sm = sm.copy()
    sm["study"] = sm["batch"].str.replace(r"_Batch_\d+$", "", regex=True)
    Z_hc = np.load(ZDIR / "Z_hc_shash.npy")
    ind_all = np.where(np.isfinite(Z), (np.abs(Z) > z_thresh).astype(float), np.nan)
    ind_hc = np.where(np.isfinite(Z_hc), (np.abs(Z_hc) > z_thresh).astype(float), np.nan)
    c = np.nansum(ind_hc, axis=0)                    # HC extreme count per gene
    n2 = np.isfinite(ind_hc).sum(axis=0)             # HC samples scored per gene
    d = n2 - c

    out = {}
    for (pheno, study), sub in sm[sm.ood_keep].groupby(["phenotype", "study"]):
        if len(sub) < min_n:
            continue
        ind_group = ind_all[sub.index.values]
        a = np.nansum(ind_group, axis=0)             # disease extreme count
        n1 = np.isfinite(ind_group).sum(axis=0)
        b = n1 - a
        N = n1 + n2
        den = n1 * n2 * (a + c) * (b + d)            # 0 when a gene is all-extreme or all-normal
        chi2 = np.divide(N * (a * d - b * c) ** 2, den, out=np.zeros(len(c)), where=den > 0)
        rate_g = np.divide(a, n1, out=np.zeros(len(c)), where=n1 > 0)
        rate_h = np.divide(c, n2, out=np.zeros(len(c)), where=n2 > 0)
        out[(pheno, study)] = np.sign(rate_g - rate_h) * np.sqrt(chi2)
    return out


def wolfers_prop_group_z(sm=None, Z=None, gene_names=None, min_n=3, z_thresh=2.6):
    """{(phenotype, study): (n_genes,) array} -- the variant already trialled in
    7_insilico_perturbation.ipynb (cell 2193089c): same |Z|>2.6 binarization, but each gene's
    extreme rate is tested with a one-sample proportion z against the GROUP'S OWN pooled
    background rate p0 (mean extreme rate over all genes x samples in that group), not against a
    separate HC cohort:  z = (phat - p0) / sqrt(p0(1-p0)/n).

    Differs from wolfers_chi2_group_z in what the null is, and the two answer different questions:
      - wolfers_chi2: 'is this gene extreme more often in disease than in HC?' (two-sample,
        diagnosis x extreme -- matches the paper's own 'associations between diagnosis and those
        scores' chi-square, and is the apples-to-apples analog of DESeq2's two-group test)
      - wolfers_prop : 'is this gene extreme more often than the average gene in this cohort?'
        (one-sample, competitive across genes -- needs no HC array, so it is immune to the
        in-sample Z_hc_shash caveat, and cancels batch-wide drift, but it is a gene-relative null)
    Kept available rather than merged so 6_group_level_comparison and 7_insilico_perturbation can
    be reconciled deliberately -- add 'wolfers_prop' to GROUP_METHODS to score it."""
    if Z is None:
        Z = np.load(ZDIR / "Z_disease_shash.npy")
    if sm is None:
        sm = pd.read_csv(ZDIR / "sample_meta.csv")
    sm = sm.copy()
    sm["study"] = sm["batch"].str.replace(r"_Batch_\d+$", "", regex=True)

    out = {}
    for (pheno, study), sub in sm[sm.ood_keep].groupby(["phenotype", "study"]):
        if len(sub) < min_n:
            continue
        Zc = Z[sub.index.values]
        ext = np.where(np.isfinite(Zc), (np.abs(Zc) > z_thresh).astype(float), np.nan)
        n = np.isfinite(ext).sum(axis=0)
        phat = np.divide(np.nansum(ext, axis=0), n, out=np.full(ext.shape[1], np.nan), where=n > 0)
        p0 = np.nanmean(ext)
        out[(pheno, study)] = (phat - p0) / np.sqrt(p0 * (1 - p0) / np.maximum(n, 1) + 1e-12)
    return out


GROUP_STAT_METHODS = {"stouffer": stouffer_group_z, "extreme_mwu": extreme_mwu_group_z,
                      "wolfers_chi2": wolfers_chi2_group_z, "wolfers_prop": wolfers_prop_group_z}
# explicit, not derived by string surgery -- a third method broke the old .replace() derivation
GROUP_STAT_LABELS = {"stouffer": "normative_group_z", "extreme_mwu": "normative_group_mwu",
                     "wolfers_chi2": "normative_group_chi2", "wolfers_prop": "normative_group_prop"}
GROUP_STAT_GSEA_LABELS = {"stouffer": "normative_group_gsea", "extreme_mwu": "normative_group_gsea_mwu",
                          "wolfers_chi2": "normative_group_gsea_chi2",
                          "wolfers_prop": "normative_group_gsea_prop"}
_GROUP_STAT_CACHE = {}


def group_stat_cached(method):
    """Memoized per-process group statistic. The normative side is design-independent, so the
    3-DESeq2-design loops in group_level_z_test / group_level_pathway_gsea would otherwise
    recompute it 3x -- extreme_mwu's per-gene Mann-Whitney is ~23s a call."""
    if method not in _GROUP_STAT_CACHE:
        _GROUP_STAT_CACHE[method] = GROUP_STAT_METHODS[method]()
    return _GROUP_STAT_CACHE[method]


def gsea_cache_key(method, pheno, study):
    """gsea_cache/ filename stem. 'stouffer' keeps the legacy 'normative__' prefix so the 21
    already-cached preranks stay valid -- a renamed key silently orphans them and costs a ~30s
    GSEA rerun each."""
    tag = "normative" if method == "stouffer" else f"normative_{method}"
    return f"{tag}__{pheno}__{study}"


def group_level_z_test(qs=(0.05, 0.10, 0.15, 0.20), deseq2_design="no_covariate", method="stouffer"):
    """Apples-to-apples group-level comparison (Rutherford 2023 eLife design): the SAME per-gene
    univariate-test -> BH-FDR -> hypergeometric-vs-OT pipeline DESeq2 uses, fed normative Z-scores
    (combined across patients via `method` -- 'stouffer' (sum(Z)/sqrt(n)) or 'extreme_mwu'
    (|Z|>2 indicator vs HC, Mann-Whitney U) -- see GROUP_STAT_METHODS) instead of raw counts --
    isolates whether deviation scores carry more disease signal than raw counts at the SAME
    statistical unit (one p-value per gene, group-level), avoiding the K=1-union across-patient
    min-p problem of normative_union."""
    from scipy.stats import hypergeom, norm
    from statsmodels.stats.multitest import multipletests
    sym_of = ensg_to_symbol()
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    sym_arr = sym_of.reindex(gene_names).values
    ref = load_reference()
    N_genes = len(gene_names)
    Z = np.load(ZDIR / "Z_disease_shash.npy")
    group_z = group_stat_cached(method)

    rows = []
    for q in qs:
        norm_sets = {}
        for (pheno, study), gz in group_z.items():
            finite = np.isfinite(gz)
            p = np.full(len(gz), np.nan)
            p[finite] = 2 * norm.sf(np.abs(gz[finite]))
            reject = np.zeros(len(gz), dtype=bool)
            reject[finite] = bh_fdr_reject(p[finite], q=q)
            sig = {s for s in sym_arr[reject] if pd.notna(s)}
            norm_sets.setdefault(pheno, set()).update(sig)  # union across studies, same convention as deseq2_gene_sets

        for pheno, sig in norm_sets.items():
            dref = ref.get(pheno, set())
            if not dref:
                continue
            K, n, x = len(dref), len(sig), len(sig & dref)
            pval = hypergeom.sf(x - 1, N_genes, K, n) if n > 0 else np.nan
            rows.append(dict(method=GROUP_STAT_LABELS[method], q=q, phenotype=pheno, n_sig=n, overlap=x, ref_size=K, pval=pval))

        ds = deseq2_gene_sets(deseq2_design, sym_of, alpha=q)
        for pheno, sig in ds.items():
            dref = ref.get(pheno, set())
            if not dref:
                continue
            K, n, x = len(dref), len(sig), len(sig & dref)
            pval = hypergeom.sf(x - 1, N_genes, K, n) if n > 0 else np.nan
            rows.append(dict(method=f"deseq2_{deseq2_design}", q=q, phenotype=pheno, n_sig=n, overlap=x, ref_size=K, pval=pval))

    df = pd.DataFrame(rows)
    _, padj, _, _ = multipletests(df.pval.fillna(1), method="fdr_bh")
    df["padj"] = padj
    df["sig_fdr05"] = padj < 0.05
    return df


def db_hit_row(phenotype, method, sig_set, ref):
    dref = ref.get(phenotype, set())
    n_sig, n_db = len(sig_set), len(sig_set & dref)
    return dict(phenotype=phenotype, method=method, n_sig=n_sig, n_db=n_db,
                db_hit_rate=round(n_db / n_sig, 3) if n_sig else np.nan,
                has_ot_ref=len(dref) > 0)


def gene_level_db_hits(save=True, designs=("no_covariate", "covariate", "ruvg_k1", "ruvg_k2", "ruvg_k3"),
                       recur_ks=(1, 3), q=0.05):
    """Symmetric gene-level DB-hit table: every DESeq2 design in `designs` + normative gene sets at
    each recurrence threshold in `recur_ks` (K=1 -> normative_union, K=3 -> normative_recur3, same
    convention as pathway-level) + per-patient normative rate distribution. `q` is the per-patient
    BH-FDR threshold for normative significance (DESeq2 sets are unaffected -- those come from
    their own pre-computed padj<0.05 calls). q!=0.05 writes to q-tagged filenames instead of the
    default cache, so relaxed-threshold sweeps don't clobber the q=0.05 result."""
    ref = load_reference()
    sym_of = ensg_to_symbol()
    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))

    design_sets = {d: deseq2_gene_sets(d, sym_of) for d in designs}
    Z = np.load(ZDIR / "Z_disease_shash.npy")
    hits = normative_gene_hits(sm, gene_names, sym_of, Z=Z, q=q)
    norm_by_k = {K: gene_sets_from_hits(hits, K=K) for K in recur_ks}

    all_phenos = set(hits)
    for ds in design_sets.values():
        all_phenos |= set(ds)

    rows = []
    for pheno in sorted(all_phenos):
        for design, ds in design_sets.items():
            if pheno in ds:
                rows.append(db_hit_row(pheno, DESIGN_LABELS[design], ds[pheno], ref))
        for K, norm_sets in norm_by_k.items():
            if pheno in norm_sets:
                label = "normative_union" if K == 1 else f'normative_recur{K}'
                rows.append(db_hit_row(pheno, label, norm_sets[pheno], ref))

    rates = pd.DataFrame(rows)

    pp_rows = []
    for pheno, patients in hits.items():
        for syms in patients:
            r = db_hit_row(pheno, "normative_persample", syms, ref)
            pp_rows.append(r)
    persample = pd.DataFrame(pp_rows)
    persample_summary = (persample[persample.has_ot_ref]
                          .groupby("phenotype")["db_hit_rate"]
                          .agg(["median", "count"]).reset_index())

    if save:
        tag = "" if q == 0.05 else f"_q{q:g}".replace(".", "")
        rates.to_csv(HERE / f"gene_db_hit_rates{tag}.csv", index=False)
        persample.to_csv(HERE / f"gene_db_hit_rates_persample{tag}.csv", index=False)
    return rates, persample_summary


# --- pathway-level: same permutation-null scorer as pathway_convergence.run_phenotype,
# fed the DESeq2 Wald stat as a single "sample" instead of a per-patient Z row.

def deseq2_pathway_sig(stat_series, gene_names_base, universe_syms, sym2idx, col2sym, terms, M,
                       cache_key=None, n_perm=800, seed=42):
    # col2sym is indexed in gene_names order (the modeled 19858 genes), not the full H5AD var table
    full = pd.Series(np.nan, index=gene_names_base)
    common = stat_series.index.intersection(full.index)
    full.loc[common] = stat_series.loc[common]
    Z = full.values[None, :]

    N = len(universe_syms)
    Zu, Fm = collapse_to_symbols(Z, col2sym, N)
    Mf = M.astype(float)
    num, den = (Zu * Fm) @ Mf.T, Fm @ Mf.T
    T = np.divide(num, den, out=np.zeros_like(num), where=den > 0)

    # the permutation null re-shuffles this exact stat vector's own values -> deterministic given
    # (stat_series, seed), so it's cacheable per (design, study, phenotype) like pathway_convergence's null.npz
    null_path = (NULL_DIR / f"{cache_key}.npz") if cache_key else None
    if null_path and null_path.exists():
        d = np.load(null_path)
        null_mean, null_sd = d["null_mean"], d["null_sd"]
    else:
        rng = np.random.default_rng(seed)
        null_sum, null_sumsq = np.zeros_like(T), np.zeros_like(T)
        for _ in range(n_perm):
            perm = rng.permutation(N)
            Zp, Fp = Zu[:, perm], Fm[:, perm]
            nu, de = (Zp * Fp) @ Mf.T, Fp @ Mf.T
            Tn = np.divide(nu, de, out=np.zeros_like(nu), where=de > 0)
            null_sum += Tn
            null_sumsq += Tn ** 2
        null_mean = null_sum / n_perm
        null_sd = np.sqrt(np.clip(null_sumsq / n_perm - null_mean ** 2, 1e-12, None))
        if null_path:
            NULL_DIR.mkdir(parents=True, exist_ok=True)
            np.savez(null_path, null_mean=null_mean, null_sd=null_sd)

    Tz = (T - null_mean) / null_sd
    from scipy.stats import norm
    p = 2 * norm.sf(np.abs(Tz))
    reject = bh_fdr_reject(p[0], q=config.PATHWAY_CONV_PARAMS["fdr_q"])
    return {t for t, r in zip(terms, reject) if r}


def reference_pathways(ref_syms, sym2idx, terms, M, q=0.05):
    """Pathways where reference genes are significantly OVER-represented (hypergeometric ORA,
    BH-FDR), not just any-overlap -- 'any member gene' saturates at 80-90% of the library once
    a phenotype has >500 reference genes, since large gene sets share members with most pathways."""
    from scipy.stats import hypergeom
    idx = [sym2idx[s] for s in ref_syms if s in sym2idx]
    if not idx:
        return set()
    N = M.shape[1]
    K = len(idx)
    n = M.sum(axis=1)
    x = M[:, idx].sum(axis=1)
    p = hypergeom.sf(x - 1, N, K, n)
    reject = bh_fdr_reject(p, q=q)
    return {t for t, r in zip(terms, reject) if r}


def normative_pathway_recurrence(pdirs):
    """{term: n_patients_hit} pooled across pdirs (a phenotype may span multiple per-study dirs)."""
    counts = {}
    for pdir in pdirs:
        d = pickle.load(open(pdir / "sig.pkl", "rb"))
        hits = d["path_sig"].sum(axis=0)
        for t, c in zip(d["terms"], hits):
            counts[t] = counts.get(t, 0) + int(c)
    return counts


def pathway_level_db_hits(save=True, designs=("no_covariate", "covariate", "ruvg_k1", "ruvg_k2", "ruvg_k3"),
                          recur_ks=(1, 3)):
    """recur_ks: normative pathway is 'detected' if significant in >=K patients (K=1 is the plain
    union -- 'at least one patient', which saturates toward the whole library as n_pat grows since
    each patient is an independent 5%-FDR test; K=3 requires independent reproduction across
    patients, closer to what DESeq2's single pooled test is actually being compared against)."""
    ref = load_reference()
    universe_syms, sym2idx, col2sym = load_symbol_vocab(None)
    terms, M = load_pathway_library()
    ref_path = {ph: reference_pathways(syms, sym2idx, terms, M) for ph, syms in ref.items()}

    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    gene_names_base = pd.Index(gene_names).str.split(".").str[0]

    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    pc_map = pc_dirs_by_phenotype(sm)
    summary = pd.read_csv(DESEQ_DIR / "summary.csv")

    rows = []
    for pheno, pdirs in pc_map.items():
        recur = normative_pathway_recurrence(pdirs)
        for K in recur_ks:
            label = "normative_union" if K == 1 else f'normative_recur{K}'
            detected = {t for t, c in recur.items() if c >= K}
            rows.append(db_hit_row(pheno, label, detected, ref_path))

        for design in designs:
            label = DESIGN_LABELS[design]
            studies = summary[(summary.design == design) & (summary.phenotype == pheno)]
            union = set()
            for _, r in studies.iterrows():
                path = DESEQ_DIR / design / f"{deseq2_tag(r.study, pheno)}.csv"
                if not path.exists():
                    continue
                res = pd.read_csv(path, index_col=0)
                res.index = res.index.str.split(".").str[0]
                cache_key = f"{design}__{deseq2_tag(r.study, pheno)}"
                sig = deseq2_pathway_sig(res["stat"].dropna(), gene_names_base, universe_syms,
                                         sym2idx, col2sym, terms, M, cache_key=cache_key)
                union |= sig
            if studies.shape[0]:
                rows.append(db_hit_row(pheno, label, union, ref_path))

    rates = pd.DataFrame(rows)
    if save:
        rates.to_csv(HERE / "pathway_db_hit_rates.csv", index=False)
    return rates


# --- group-level pathway comparison via real preranked GSEA (Subramanian 2005), replacing the
# ad hoc mean-Z + gene-label-permutation scorer above with the standard weighted running-sum
# algorithm -- possible now that both methods produce ONE group-level ranking statistic per gene
# (Stouffer Z for normative, Wald stat for DESeq2), the same apples-to-apples design as
# group_level_z_test at gene level.

def gsea_prerank(rnk, terms, M, universe_syms, n_perm=1000, seed=42, min_size=5, max_size=1000,
                 threads=8):
    """Standard preranked GSEA on a gene-symbol -> score ranking, using the SAME KEGG+Reactome
    housekeeping-filtered library as PerSamplePathwayAnalysis (`load_pathway_library`). Returns gseapy's
    res2d (Term, NES, FDR q-val, ...)."""
    import gseapy as gp
    gene_sets = {t: [universe_syms[j] for j in np.where(M[ti])[0]] for ti, t in enumerate(terms)}
    rnk_s = pd.Series(rnk, index=universe_syms).dropna().sort_values(ascending=False)
    res = gp.prerank(rnk=rnk_s, gene_sets=gene_sets, min_size=min_size, max_size=max_size,
                     permutation_num=n_perm, seed=seed, outdir=None, no_plot=True, threads=threads)
    return res.res2d


def _cached_gsea(cache_key, rnk, terms, M, universe_syms, **kw):
    GSEA_DIR.mkdir(parents=True, exist_ok=True)
    path = GSEA_DIR / f"{slugify(cache_key)}.csv"
    if path.exists():
        return pd.read_csv(path)
    res2d = gsea_prerank(rnk, terms, M, universe_syms, **kw)
    res2d.to_csv(path, index=False)
    return res2d


def group_level_pathway_gsea(save=True, deseq2_design="no_covariate", gsea_qs=(0.05, 0.25), method="stouffer"):
    """Pathway-level apples-to-apples: preranked GSEA on group normative Z (`method` -- 'stouffer' or
    'extreme_mwu', see GROUP_STAT_METHODS) vs DESeq2 Wald stat, same library, same algorithm.
    Reports hit-rate vs the Open Targets reference pathway set (`reference_pathways`, hypergeometric
    ORA on the DB gene list itself) at GSEA's own q<0.05 confirmatory / q<0.25 discovery tiers
    (Subramanian 2005 convention, see 논문화/중간결과/Benchmark/FDR_THRESHOLD_RATIONALE.md)."""
    universe_syms, sym2idx, col2sym = load_symbol_vocab(None)
    terms, M = load_pathway_library()
    N = len(universe_syms)
    ref = load_reference()
    ref_path = {ph: reference_pathways(syms, sym2idx, terms, M) for ph, syms in ref.items()}

    sm = pd.read_csv(ZDIR / "sample_meta.csv")
    gene_names = pickle.load(open(ZDIR / "gene_names.pkl", "rb"))
    Z = np.load(ZDIR / "Z_disease_shash.npy")
    group_z = group_stat_cached(method)

    rows = []
    norm_by_q = {}
    for (pheno, study), gz in group_z.items():
        Zu, Fm = collapse_to_symbols(gz[None, :], col2sym, N)
        rnk = np.where(Fm[0] > 0, Zu[0], np.nan)
        res2d = _cached_gsea(gsea_cache_key(method, pheno, study), rnk, terms, M, universe_syms)
        for q in gsea_qs:
            sig = set(res2d.loc[res2d["FDR q-val"] < q, "Term"])
            norm_by_q.setdefault((pheno, q), set()).update(sig)  # union across studies

    label_base = GROUP_STAT_GSEA_LABELS[method]
    for (pheno, q), sig in norm_by_q.items():
        label = label_base if q == 0.05 else f"{label_base}_q{q:g}".replace(".", "")
        rows.append(db_hit_row(pheno, label, sig, ref_path))

    sym_of = ensg_to_symbol().copy()  # .copy() -- ensg_to_symbol is lru_cached, never reindex in place
    sym_of.index = sym_of.index.str.split(".").str[0]
    summary = pd.read_csv(DESEQ_DIR / "summary.csv")
    summary = summary[summary.design == deseq2_design]
    for pheno, sub in summary.groupby("phenotype"):
        by_q = {q: set() for q in gsea_qs}
        for _, r in sub.iterrows():
            path = DESEQ_DIR / deseq2_design / f"{deseq2_tag(r.study, pheno)}.csv"
            if not path.exists():
                continue
            res = pd.read_csv(path, index_col=0)
            res.index = res.index.str.split(".").str[0]
            stat = res["stat"].dropna()
            syms = sym_of.reindex(stat.index)
            rnk_df = pd.Series(stat.values, index=syms.values)
            rnk_df = rnk_df[pd.notna(rnk_df.index)].groupby(level=0).mean().reindex(universe_syms)
            res2d = _cached_gsea(f"{deseq2_design}__{deseq2_tag(r.study, pheno)}", rnk_df.values,
                                 terms, M, universe_syms)
            for q in gsea_qs:
                by_q[q] |= set(res2d.loc[res2d["FDR q-val"] < q, "Term"])
        if sub.shape[0]:
            for q in gsea_qs:
                label = f"{DESIGN_LABELS[deseq2_design]}_gsea" if q == 0.05 else \
                    f"{DESIGN_LABELS[deseq2_design]}_gsea_q{q:g}".replace(".", "")
                rows.append(db_hit_row(pheno, label, by_q[q], ref_path))

    rates = pd.DataFrame(rows)
    if save:
        tag = "" if method == "stouffer" else f"_{method}"
        rates.to_csv(HERE / f"pathway_gsea_db_hit_rates{tag}.csv", index=False)
    return rates
