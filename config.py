from pathlib import Path

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "OpenAccess_nfcore"
PIPELINE_DIR = ROOT / "Saved_Pipeline"

PATHS = {
    "merged_raw":    DATA_DIR / "Merged_Processed_AnnData.h5ad",
    "merged_biases": DATA_DIR / "Merged_Processed_AnnData_with_Batch_Biases.h5ad",
    "merged_qc":     DATA_DIR / "Merged_Processed_AnnData_with_Batch_Biases_QC_Status.h5ad",
    "pipeline_meta": PIPELINE_DIR / "LogisticGP_meta.pkl",
    "pipeline_lr":   PIPELINE_DIR / "LogisticGP_lr_model.pt",
    "pipeline_gp":   PIPELINE_DIR / "LogisticGP_gp_states.pt",
    "z_total":       PIPELINE_DIR / "Z_Total_Matrix.pkl",
    "z_state":       PIPELINE_DIR / "Z_State_Matrix.pkl",
    "z_quant":       PIPELINE_DIR / "Z_Quant_Matrix.pkl",
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
