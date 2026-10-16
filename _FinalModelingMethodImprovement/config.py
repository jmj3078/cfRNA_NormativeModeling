import os
import shutil
import sys
from pathlib import Path

FMMI = Path(__file__).resolve().parent


def _main_worktree(root):
    """Data and production outputs are untracked, so they live only in the main worktree."""
    git = root / ".git"
    return Path(git.read_text().split("gitdir:")[1].strip().split("/.git/")[0]) if git.is_file() else root


PROJECT_ROOT = _main_worktree(FMMI.parent)

H5AD_PATH = PROJECT_ROOT / "OpenAccess_nfcore" / "Merged_Processed_AnnData_with_Batch_Biases_QC_Status.h5ad"
RSCRIPT = shutil.which("Rscript") or str(Path(sys.executable).resolve().parent / "Rscript")

SPLITS_DIR = FMMI / "splits"
CV_FOLDS_PATH = SPLITS_DIR / "cv_folds.json"
LOBO_FOLDS_PATH = SPLITS_DIR / "lobo_folds.json"
INSILICO_DIR = FMMI / "insilico"
RESULTS_DIR = FMMI / "results"
INSILICO_CSV = RESULTS_DIR / "insilico_sweep.csv"
DETECTION_LIMIT_CSV = RESULTS_DIR / "detection_limit.csv"
OUTLIER_ARMS = {"autoencoder": "OUTRIDER (AE)", "pca": "PCA", "outsingle": "OutSingle", "peer": "PEER"}
COVARIATE_DIR = RESULTS_DIR / "covariates"
PANGLAO_TSV = PROJECT_ROOT / "Data" / "PalangoDB_CellTypeMarkers.tsv"
RUVSEQ_RSCRIPT = Path.home() / "miniconda3" / "envs" / "ruvseq_env" / "bin" / "Rscript"

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

# covariates: obs columns as they are, "W_<i>" = RUVg factor i of control_genes (engine/core/covariates.py).
# A new covariate set is a new entry here, with its own runs/<pipeline>/ and its own pcis_cut
# (the PCIS null depends on the design; eval/fit_pipeline.sh <pipe> calib prints it).
PIPELINES = {
    "p0": dict(covariates=BIAS_COLUMNS),  # current production pipeline
    "p2": dict(covariates=["log(Total Reads)"] + [f"W_{i}" for i in range(1, 7)],
               control_genes=SPLITS_DIR / "control_genes_N500.txt",  # 0_build_covariates.ipynb, N=500 + sex genes
               pcis_cut=2.561541),
}
PIPELINE = os.environ.get("FMMI_PIPELINE", "p0")
assert PIPELINE in PIPELINES, f"unknown pipeline {PIPELINE!r}"
COVARIATES = PIPELINES[PIPELINE]["covariates"]
CONTROL_GENES = PIPELINES[PIPELINE]["control_genes"].read_text().split() if "control_genes" in PIPELINES[PIPELINE] else []

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
    "pcis_cut": PIPELINES[PIPELINE].get("pcis_cut", 3.62),  # 1e-4 per-observation null rate, eval/fit_pipeline.sh calib
    "max_outlier_frac": 0.05,
    "chunk_size": 200,
    "cores": 12,
}
