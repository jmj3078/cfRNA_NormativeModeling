import os
import shutil
import sys
from pathlib import Path

FMMI = Path(__file__).resolve().parent
PROJECT_ROOT = FMMI.parent

H5AD_PATH = PROJECT_ROOT / "OpenAccess_nfcore" / "Merged_Processed_AnnData_with_Batch_Biases_QC_Status.h5ad"
RSCRIPT = shutil.which("Rscript") or str(Path(sys.executable).resolve().parent / "Rscript")

PIPELINES = {
    "p1_baseline": dict(covariates="bias10", rqr="marginal"),
    "p2_blup": dict(covariates="bias10", rqr="blup"),
}
PIPELINE = os.environ.get("FMMI_PIPELINE", "p1_baseline")
assert PIPELINE in PIPELINES, f"unknown pipeline {PIPELINE!r}"

SPLITS_DIR = FMMI / "splits"
CV_FOLDS_PATH = SPLITS_DIR / "cv_folds.json"
LOBO_FOLDS_PATH = SPLITS_DIR / "lobo_folds.json"
INSILICO_DIR = FMMI / "insilico"
RESULTS_DIR = FMMI / "results"

RUN_DIR = FMMI / "runs" / PIPELINE
ENGINE_MIXED_DIR = RUN_DIR / "full"
CV_MIXED_DIR = RUN_DIR / "cv"
LOBO_MIXED_DIR = RUN_DIR / "lobo"
TMP_DIR = RUN_DIR / "tmp"
DISPERSION_TREND_PATH = ENGINE_MIXED_DIR / "dispersion_trend.json"
DISP_PRIOR_PATH = ENGINE_MIXED_DIR / "disp_prior.json"
PCIS_CAL_DIR = RUN_DIR / "pcis_calibration"
PCIS_CAL_FIG_DIR = PCIS_CAL_DIR / "Figures"

GLMM_FIT_R = FMMI / "engine" / "core" / "glmm_fit.R"
GLMM_FIT_POOL_R = FMMI / "engine" / "core" / "glmm_fit_pool.R"

# Read-only: outputs of the production pipeline, consumed by baseline_import.py only.
BASELINE_MEM_DIR = PROJECT_ROOT / "MixedEffectsModeling"

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
NZ_A_MAX = 31
MIN_HC_BATCH_SIZE = 5
LOBO_MIN_TEST_HC = 25
LOBO_REQUIRE_DISEASE = True

SPIKE_PARAMS = {
    "beta_explode_thr": 3.0,
    "seed": 42,
    "rare_overdisp_thr": 2.0,
    "alpha_floor": 1e-2,
    "alpha_cap": 50.0,
    "n_splits": 5,
    "trend_min_nz": 30,
}

EB_PARAMS = {
    "calib_n_genes": 2000,
    "calib_n_strata": 10,
    "tau_floor": 1e-3,
}

FIT_PARAMS = {
    "beta_explode_thr": SPIKE_PARAMS["beta_explode_thr"],
    "tau2_max": 3.0,
    "disp_intercept_max": 10.0,
    "pcis_cut": 3.62,
    "max_outlier_frac": 0.05,
    "chunk_size": 200,
    "cores": 12,
}
