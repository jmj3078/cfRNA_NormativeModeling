import os
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
CTRL_COMP_STAT_DIR = CTRL_COMP_DIR / "tstats"
CTRL_COMP_FIG_DIR = CTRL_COMP_DIR / "Figures"
CTRL_COMP_DESEQ2_DIR = CTRL_COMP_DIR / "deseq2_stats"

# Sample scopes for the control-composition experiment. "b1" is the single pre-analytical
# condition (CEDAR, 2x150); "b123" pools it with the two BCC batches, which is what the old
# Batch_ID definition called "Moore et al._Batch_1" before 2026-09-30. Each scope writes to
# its own subdirectory and CSV, so both can be held side by side.
CTRL_COMP_SCOPES = {
    "b1": ["Moore et al._Batch_1"],
    "b123": ["Moore et al._Batch_1", "Moore et al._Batch_2", "Moore et al._Batch_3"],
}
CTRL_COMP_DEFAULT_SCOPE = "b1"

# Disjoint control strata per random split, per scope. Strata never share a sample, so an
# overlapping subsample design is not on the table: drawing 30 of 48 twice forces >=12 shared
# controls and inflates agreement. b1 uses 2 (HC 48 -> 24 each, the largest disjoint arms
# available); b123 stays at 3 because its results predate the change and were not refitted.
CTRL_COMP_N_SPLITS = {"b1": 2, "b123": 3}


def ctrl_comp_n_splits(scope=None):
    return CTRL_COMP_N_SPLITS[scope or CTRL_COMP_DEFAULT_SCOPE]

EDA_GROUP_CONTRAST_DIR = EDA_RESULTS_DIR / "GroupContrastLimits"
EDA_GROUP_CONTRAST_FIG_DIR = EDA_GROUP_CONTRAST_DIR / "Figures"
# Two DEG runs of the same Moore pancreatic cohort, kept side by side rather than one
# superseding the other: "pancreatic_deg" is CEDAR alone (Batch_1, 153 samples, what Moore
# et al. fitted on), "pancreatic_deg_pooled" is the older CEDAR+BCC pool (Batch_1+2+3) from
# before the 2026-09-30 Batch_ID redefinition.
PANCREATIC_DEG_DIR = EDA_GROUP_CONTRAST_DIR / "pancreatic_deg"
PANCREATIC_DEG_POOLED_DIR = EDA_GROUP_CONTRAST_DIR / "pancreatic_deg_pooled"

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


MOORE2025_REPO_DIR = ROOT / "Data" / "2024_ngolab_cfRNA" / "Analysis_Scripts"
MOORE2025_DE_R = _MEM / "notebookscripts" / "GroupVsIndividual" / "moore2025_de.R"
GLMM_FIT_R     = _MEM / "core" / "glmm_fit.R"
GLMM_FIT_POOL_R = _MEM / "core" / "glmm_fit_pool.R"
HPA_SINGLE_CELL_ZIP = ROOT / "Data" / "HPA" / "rna_single_cell_type.tsv.zip"  # HPA 25.1, Data/HPA/VERSION.txt
HPA_ATLAS_ZIP = ROOT / "Data" / "HPA" / "proteinatlas.tsv.zip"  # HPA 25.1 gene-level specificity categories
GTEX_WB_GCT = ROOT / "Data" / "GTEx" / "gene_reads_v10_whole_blood.gct"
PANGLAO_TSV = ROOT / "Data" / "PalangoDB_CellTypeMarkers.tsv"
RUVSEQ_RSCRIPT = Path.home() / "miniconda3" / "envs" / "ruvseq_env" / "bin" / "Rscript"


STRATIFY_COL = "Batch_ID"

# Pooling cutoff. Set from the nz_a_max=0 run (every gene through the individual cascade).
# Raised 25 -> 31 (2026-07-29): CV fold-level convergence (all 5 folds) drops to ~0.6 in the
# nz 25-30 bin, below the bar the nz_a_max choice is meant to hold.
NZ_A_MAX = 31
MIN_HC_BATCH_SIZE = 5

# Held-out batch selection for LOBO-based evaluation (build_lobo_folds.py,
# 3_detection_limit.ipynb). Derived from the fitted HC pool instead of a hardcoded
# batch list, so a change to the Batch_ID definition propagates automatically.
LOBO_MIN_TEST_HC = 25       # held-out HC needed for a batch to be scoreable
LOBO_REQUIRE_DISEASE = True # tier A: batch must also carry disease samples

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
    "pcis_cut": 3.83,
    "max_outlier_frac": 0.05,
    "chunk_size": 200,
    "cores": 12,
}


# Covariate pipelines. "W_<i>" = RUVg factor i of control_genes (core/covariates.py); any other name is an obs column.
# Every output and cache of a pipeline lives under MixedEffectsModeling/Notebook_<pipeline>/, next to its notebooks.
# The pipeline is chosen by the NM_PIPELINE environment variable only; without it no output path is defined,
# so nothing can read or write the other pipeline's results.
PIPELINES = {
    "p0": dict(covariates=BIAS_COLUMNS, pcis_cut=3.83),  # frozen: from a null run before per-gene seeding, not regenerable
    "p3": dict(covariates=["log(Total Reads)"] + [f"W_{i}" for i in range(1, 7)] + ["Gene Length Bias", "GC Bias"],
               control_genes=_MEM / "core" / "control_genes_blood.txt",
               pcis_cut=3.641121),  # frozen: FMMI p3 fit, null run before per-gene seeding (seeded NREP 5 rerun: 3.6488)
}
PIPELINE = os.environ.get("NM_PIPELINE")
if PIPELINE is not None:
    assert PIPELINE in PIPELINES, f"unknown NM_PIPELINE {PIPELINE!r}"
    COVARIATES = PIPELINES[PIPELINE]["covariates"]
    CONTROL_GENES = PIPELINES[PIPELINE]["control_genes"].read_text().split() if "control_genes" in PIPELINES[PIPELINE] else []
    FIT_PARAMS["pcis_cut"] = PIPELINES[PIPELINE]["pcis_cut"]
    PIPELINE_DIR = _MEM / f"Notebook_{PIPELINE}"
    TMP_ROOT = Path(os.environ.get("NM_TMP", f"/tmp/nm_{PIPELINE}"))
    ENGINE_DIR = PIPELINE_DIR / "EngineState"
    CV_DIR = PIPELINE_DIR / "1_CVResults"
    LOBO_DIR = PIPELINE_DIR / "2_LOBOResults"
    ZSCORES_DIR = PIPELINE_DIR / "ZScores"
    PCIS_CAL_DIR = PIPELINE_DIR / "PCISCalibration"
    PCIS_CAL_FIG_DIR = PCIS_CAL_DIR / "Figures"
    DISEASE_SCORING_FIG_DIR = PIPELINE_DIR / "4_DiseaseScoring" / "Figures"
    DISEASE_REF_DIR = PIPELINE_DIR / "Benchmark" / "DiseaseReference"
    GROUP_VS_INDIV_DIR = PIPELINE_DIR / "7_PancreaticCancerCohortsAnalysis"
    GROUP_VS_INDIV_FIG_DIR = GROUP_VS_INDIV_DIR / "Figures"
    MOORE2025_DE_DIR = GROUP_VS_INDIV_DIR / "Moore2025DE"
    OUTRIDER_COMPARISON_DIR = PIPELINE_DIR / "0_OutriderComparison" / "InsampleComparison"
    OUTRIDER_HELD_OUT_DIR = PIPELINE_DIR / "0_OutriderComparison" / "HeldOutComparison"
    DETECTION_LIMIT_DIR = PIPELINE_DIR / "3_DetectionLimitResults"
    METHOD_COMP_DIR = PIPELINE_DIR / "5_MethodComparison"
    METHOD_COMP_CACHE_DIR = METHOD_COMP_DIR / "cache"
    METHOD_COMP_FIG_DIR = METHOD_COMP_DIR / "Figures"
    DISEASE_TISSUE_DIR = PIPELINE_DIR / "6_DiseaseTissue"
    CONTROL_GENE_DIR = PIPELINE_DIR / "ControlGenes"
    DISPERSION_TREND_PATH = ENGINE_DIR / "dispersion_trend.json"
    DISP_PRIOR_PATH = ENGINE_DIR / "disp_prior.json"
