import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# subprocess.run(["Rscript", ...]) fails if launched without the conda env's
# bin/ on PATH (e.g. nohup without `conda activate`); fall back to the
# running interpreter's own env, which ships Rscript alongside python.
RSCRIPT = shutil.which("Rscript") or str(Path(sys.executable).resolve().parent / "Rscript")
H5AD_PATH = ROOT / "OpenAccess_nfcore" / "Merged_Processed_AnnData_with_Batch_Biases_QC_Status.h5ad"
GAMLSS_R_HELPER = ROOT / "Modeling" / "gamlss.r"

_HERE = Path(__file__).resolve().parent

ENGINE_MIXED_DIR        = _HERE / "engine_state_mixed"
CV_MIXED_DIR            = _HERE / "CV_Results_mixed"
CV_MIXED_FIG_DIR        = CV_MIXED_DIR / "Figures"
LOBO_MIXED_DIR          = _HERE / "LOBO_Results_mixed"
ZSCORES_MIXED_DIR       = _HERE / "Z_scores_mixed"
THRESHOLD_SWEEP_DIR     = _HERE / "Threshold_Sweep"
THRESHOLD_SWEEP_FIG_DIR = THRESHOLD_SWEEP_DIR / "Figures"
PCIS_CAL_DIR            = _HERE / "PCIS_Calibration"
PCIS_CAL_FIG_DIR        = PCIS_CAL_DIR / "Figures"
PATHWAY_CONV_DIR        = _HERE / "PerSamplePathwayAnalysis"
PATHWAY_CONV_FIG_DIR    = PATHWAY_CONV_DIR / "Figures"
SIGNAL_TREND_DIR        = _HERE / "SignalTrendAnalysis"
SIGNAL_TREND_CUR_DIR    = SIGNAL_TREND_DIR / "PathwayCuration"
SIGNAL_TREND_FIG_DIR    = SIGNAL_TREND_DIR / "Figures"
BENCHMARK_DIR           = _HERE / "Benchmark"
DISEASE_SCORING_DIR     = _HERE / "DiseaseScoring"
DISEASE_SCORING_FIG_DIR = DISEASE_SCORING_DIR / "Figures"
OUTRIDER_COMPARISON_DIR = _HERE / "OutriderComparison" / "insample_comparison"
OUTRIDER_HELD_OUT_DIR   = _HERE / "OutriderComparison" / "held_out_comparison"
PERTURBATION_DIR        = _HERE / "PerturbationResults"
GLMM_HELPERS_R = _HERE / "core" / "glmm_helpers.R"
PCIS_NULL_R    = _HERE / "core" / "pcis_null.R"
GLMM_FIT_R     = _HERE / "core" / "glmm_fit.R"
GLMM_FIT_POOL_R = _HERE / "core" / "glmm_fit_pool.R"
DISPERSION_TREND_PATH = ENGINE_MIXED_DIR / "dispersion_trend.json"
DISP_PRIOR_PATH = ENGINE_MIXED_DIR / "disp_prior.json"

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
    "n_null_perm": 800,  # used only by the legacy permutation-null path in 5_gene_pathway_reoccurence.ipynb
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
    "max_pathway_size_select": 300,
    "top_k_pathways": 6,
    "redundancy_jaccard_max": 0.5,
    "exclude_keywords": [
        "oxidative phosphorylation", "electron transport", "respiratory chain",
        "alzheimer", "parkinson", "huntington", "prion disease", "amyotrophic lateral sclerosis",
    ],
}
