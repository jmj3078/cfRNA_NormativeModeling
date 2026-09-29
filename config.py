import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "OpenAccess_nfcore"

PATHS = {
    "merged_raw":    DATA_DIR / "Merged_Processed_AnnData.h5ad",
    "merged_biases": DATA_DIR / "Merged_Processed_AnnData_with_Batch_Biases.h5ad",
    "merged_qc":     DATA_DIR / "Merged_Processed_AnnData_with_Batch_Biases_QC_Status.h5ad",
}

PARAMS = {
    "min_study_samples": 10,
    "n_top_genes":       2000,
    "n_pcs":             50,
    "loess_frac":        0.7,
    "n_bins":            20,
    "outlier_pct":       99,
    "min_expressed":     50,
}

H5AD_PATH = PATHS["merged_qc"]   # normative modeling

EDA_RESULTS_DIR   = ROOT / "EDA" / "Analysis_Results"
EDA_OVERVIEW_DIR  = EDA_RESULTS_DIR / "Overview"
EDA_QC_DIR        = EDA_RESULTS_DIR / "QC"
EDA_PCA_NOHVG_DIR = EDA_RESULTS_DIR / "PCA" / "NoHVG"
EDA_PCA_HVG_DIR   = EDA_RESULTS_DIR / "PCA" / "HVG"
EDA_RDA_DIR       = EDA_RESULTS_DIR / "RDA"
EDA_RDA_NOHVG_DIR = EDA_RDA_DIR / "NoHVG"
EDA_RDA_HVG_DIR   = EDA_RDA_DIR / "HVG"
EDA_RDA_HC_DIR    = EDA_RDA_DIR / "HC"
EDA_BIAS_BATCH_DIR = EDA_RESULTS_DIR / "BiasBatch"
EDA_BIAS_PHENOTYPE_DIR = EDA_RESULTS_DIR / "BiasPhenotype"
EDA_GENE_BIAS_HC_DIR = EDA_RESULTS_DIR / "GeneBiasHC"
EDA_VIF_DIR       = EDA_RESULTS_DIR / "VIF"

CTRL_COMP_DIR     = EDA_RESULTS_DIR / "Control_Composition"
CTRL_COMP_W_DIR   = CTRL_COMP_DIR / "ruvg_W"
CTRL_COMP_EXPR_DIR = CTRL_COMP_DIR / "expr"
CTRL_COMP_STAT_DIR = CTRL_COMP_DIR / "tstats"
CTRL_COMP_FIG_DIR = CTRL_COMP_DIR / "Figures"
CTRL_COMP_DESEQ2_DIR = CTRL_COMP_DIR / "deseq2_stats"

EDA_GROUP_CONTRAST_DIR = EDA_RESULTS_DIR / "GroupContrastLimits"
EDA_GROUP_CONTRAST_FIG_DIR = EDA_GROUP_CONTRAST_DIR / "Figures"

BIAS_COLUMNS = [
    "log(Total Reads)",
    "Spliced Reads (%)",
    "gDNA Contamination (Intron/Exon)",
    "rRNA Fraction",
    "RNA Degradation (3' Bias)",
    "Platelet Score",
    "GC Bias",
    "Gene Length Bias",
    "NG80",
    "(NP80/NG80)",
]

EDA_BIAS_METRICS = [
    "log(Total Reads)", "Spliced Reads (%)", "gDNA Contamination (Intron/Exon)",
    "rRNA Fraction", "Platelet Score", "GC Bias", "Gene Length Bias",
    "RNA Degradation (3' Bias)", "NG80", "NP80", "(NP80/NG80)",
]

EDA_COMBINED_METRICS = [
    "log(Total Reads)", "Spliced Reads (%)", "gDNA Contamination (Intron/Exon)",
    "rRNA Fraction", "RNA Degradation (3' Bias)", "Platelet Score",
    "GC Bias", "Gene Length Bias", "NG80", "NP80", "(NP80/NG80)",
    "Sample Volume (mL)", "Total Centrifugation Force (g)",
]

# Used by EDA/control_composition/run_control_composition.py's MahalanobisFilter.
MODELING_PARAMS = {"ood_percentile": 95}

_MEM = ROOT / "MixedEffectsModeling"


# subprocess.run(["Rscript", ...]) fails if launched without the conda env's
# bin/ on PATH (e.g. nohup without `conda activate`); fall back to the
# running interpreter's own env, which ships Rscript alongside python.
RSCRIPT = shutil.which("Rscript") or str(Path(sys.executable).resolve().parent / "Rscript")


ENGINE_MIXED_DIR        = _MEM / "engine_state_mixed"
CV_MIXED_DIR            = _MEM / "CV_Results_mixed"
LOBO_MIXED_DIR          = _MEM / "LOBO_Results_mixed"
ZSCORES_MIXED_DIR       = _MEM / "Z_scores_mixed"
THRESHOLD_SWEEP_DIR     = _MEM / "Threshold_Sweep"
PCIS_CAL_DIR            = _MEM / "PCIS_Calibration"
PCIS_CAL_FIG_DIR        = PCIS_CAL_DIR / "Figures"
PATHWAY_CONV_DIR        = _MEM / "Cohort_Grouped_Z"
DISEASE_SCORING_FIG_DIR = _MEM / "DiseaseScoring" / "Figures"
DISEASE_REF_DIR         = _MEM / "Benchmark" / "disease_reference"
GROUP_VS_INDIV_DIR      = _MEM / "GroupVsIndividual"
GROUP_VS_INDIV_FIG_DIR  = GROUP_VS_INDIV_DIR / "Figures"
PANCREATIC_DEG_DIR      = CTRL_COMP_DIR / "pancreatic_deg"
OUTRIDER_COMPARISON_DIR = _MEM / "OutriderComparison" / "insample_comparison"
OUTRIDER_HELD_OUT_DIR   = _MEM / "OutriderComparison" / "held_out_comparison"
DETECTION_LIMIT_DIR     = _MEM / "DetectionLimitResults"
METHOD_COMP_DIR         = _MEM / "MethodComparison"
METHOD_COMP_CACHE_DIR   = METHOD_COMP_DIR / "cache"
METHOD_COMP_FIG_DIR     = METHOD_COMP_DIR / "Figures"
GLMM_FIT_R     = _MEM / "core" / "glmm_fit.R"
GLMM_FIT_POOL_R = _MEM / "core" / "glmm_fit_pool.R"
DISPERSION_TREND_PATH = ENGINE_MIXED_DIR / "dispersion_trend.json"
DISP_PRIOR_PATH = ENGINE_MIXED_DIR / "disp_prior.json"


STRATIFY_COL = "Batch_ID"

# Pooling cutoff. Set from the nz_a_max=0 run (every gene through the individual cascade).
# Raised 25 -> 31 (2026-07-29): CV fold-level convergence (all 5 folds) drops to ~0.6 in the
# nz 25-30 bin, below the bar the nz_a_max choice is meant to hold.
NZ_A_MAX = 31
MIN_HC_BATCH_SIZE = 5

SPIKE_PARAMS = {
    "beta_explode_thr": 3.0,
    "seed": 42,
    "rare_overdisp_thr": 2.0,
    "alpha_floor": 1e-2,
    "alpha_cap": 50.0,
    "n_splits": 5,
    "trend_min_nz": 30,
}

# Empirical-Bayes dispersion shrinkage + PCIS (Prior-Conditioned Impact Score)
# outlier removal
EB_PARAMS = {
    "calib_n_genes": 2000,
    "calib_n_strata": 10,
    "tau_floor": 1e-3,
}

FIT_PARAMS = {
    "beta_explode_thr": SPIKE_PARAMS["beta_explode_thr"],
    "tau2_max": 3.0,
    "disp_intercept_max": 10.0,
    "pcis_cut": 2.25,
    "max_outlier_frac": 0.05,
    "chunk_size": 200,
    "cores": 12,
}

# Gene- vs pathway-level deviation convergence (4_gene_enrichment.ipynb): patient-level BH-sig
# genes are heterogeneous, but does the same deviation converge onto shared pathways? Mirrors
# Wolfers 2018 JAMA Psych / Segal 2023 Nat Neurosci deviation-overlap design (see
# EDA/normative_modeling_literature.md).
PATHWAY_CONV_PARAMS = {
    # GO_Biological_Process tried and dropped: checked its top-scoring (recur*eff) terms for
    # Tuberculosis and 24/25 were near-duplicates (Jaccard>=0.3) of an existing KEGG/Reactome term
    # or too generic (mRNA splicing, transcription regulation, glycolysis, mitosis) to be a disease
    # story -- GO's fine-grained hierarchy mostly re-slices signal KEGG/Reactome already carry.
    "gene_sets": ["KEGG_2021_Human", "Reactome_2022"],
    "min_pathway_size": 5,
    # per-sample gene-level cutoff feeding the pathway hypergeometric/Fisher ORA test. Method
    # comparison (_scratch_pathway_methods/, 2026-08) benchmarked HC-population-null mean-Z,
    # CAMERA-style PAGE, singscore, and this |Z|-threshold + Fisher ORA against a negative
    # control (held-out HC samples scored as if they were patients, true null): singscore was
    # badly anti-conservative (up to 14.5% of pathways "significant" in healthy controls),
    # HC-population-null was badly batch-confounded (r=0.67 between hit count and |global
    # sample-mean Z|, driven by 2 specific batches), CAMERA was underpowered even in real disease
    # samples. Fisher ORA was the only one clean on the negative control (median 0 across all
    # thresholds tested) while still detecting signal in disease samples -- adopted as the
    # pipeline default. z_thresh=1.96 (nominal two-sided p<0.05) empirically beat looser (1.64,
    # dilutes the enrichment ratio with background noise genes) and stricter (2.33/2.58, too few
    # genes left for hypergeometric power) alternatives in a 4-point sweep on Tuberculosis.
    "z_thresh": 1.96,
    # kept at the nominal 0.05 default for the pipeline's own path_sig/path_sig_up/path_sig_down --
    # p_path/p_up/p_down (pre-BH hypergeometric p-values) are cached in sig.pkl/sig_directional.pkl
    # regardless of q, so a q-sweep for reoccurrence analysis is done by re-thresholding those cached
    # p-values in the notebook (5_gene_pathway_reoccurence.ipynb sec. 1), not by rerunning the engine.
    "fdr_q": 0.05,
    "seed": 42,
    # Blood/cfRNA transcriptomics has a literature-recognized confound here, not just an in-house
    # observation: Chaussabel et al. 2008 Immunity (PMID 18631455) modular blood-transcriptomics
    # framework identifies a coordinately-expressed "protein synthesis / ribosomal protein" module
    # that dominates variance in whole-blood/PBMC data and reflects generic translational activity or
    # cell-composition shift, not disease-specific biology -- reused for the same purpose in
    # Rinchai/Chaussabel 2020 (PMID 32736569), Vegh/Chaussabel 2019 (PMID 31253760). Goeman & Buhlmann
    # 2007 (PMID 17303618) gives the general mechanism: gene sets sharing a highly co-regulated block
    # are vulnerable to spurious enrichment regardless of the set's nominal biology. Name-based keyword
    # match alone misses pathways that carry this module by gene COMPOSITION but not by NAME (Influenza
    # Infection, SLIT/ROBO signaling, Cellular Response To Starvation all came out >45% ribosomal-protein
    # genes empirically here) -- so exclusion is composition-based: any pathway sharing > ribo_frac_max
    # of its genes with the reference KEGG "Ribosome" set (as an operational proxy for the Chaussabel
    # module) is dropped. The KEGG-Ribosome proxy and the 0.15 cutoff are our own operational choices,
    # not literature-derived -- Chaussabel's framework flags the module qualitatively, no numeric cutoff.
    # Keyword list stays as a fast belt-and-suspenders for OXPHOS/neurodegeneration, which the
    # ribosome-composition check does not catch (feedback_gsea_interpretation).
    "ribo_reference_term": "Ribosome",
    "ribo_frac_max": 0.15,
    "exclude_keywords": [
        "oxidative phosphorylation", "electron transport", "respiratory chain",
        "alzheimer", "parkinson", "huntington", "prion disease", "amyotrophic lateral sclerosis",
    ],
}
