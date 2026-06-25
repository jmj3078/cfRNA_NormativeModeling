from pathlib import Path

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "OpenAccess_nfcore"
PIPELINE_DIR = ROOT / "Saved_Pipeline"

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

# ---------------------------------------------------------------------------
# Modeling pipeline 
# ---------------------------------------------------------------------------
MODELING_DIR   = ROOT / "Modeling"
ENGINE_DIR     = MODELING_DIR / "engine_state"
CV_RESULTS_DIR = MODELING_DIR / "CV_Results"
CV_FIG_DIR     = CV_RESULTS_DIR / "Figures"
Z_SCORES_DIR   = MODELING_DIR / "Z_scores"          # normative model Z-score 행렬 전용
GSEA_DIR       = MODELING_DIR / "GSEA"
GSEA_FIG_DIR   = GSEA_DIR / "Figures"
RARE_REF       = Z_SCORES_DIR / "rare_event_ref.pkl"
R_HELPER       = MODELING_DIR / "gamlss.r"

H5AD_PATH = PATHS["merged_qc"]   # normative modeling

# scoring 산출 Z-score 행렬 (disease_scoring/scoring.py) → Z_scores/ 디렉토리
Z_DISEASE      = Z_SCORES_DIR / "Z_disease.npy"
Z_SAMPLE_NAMES = Z_SCORES_DIR / "Z_sample_names.npy"
Z_GENE_NAMES   = Z_SCORES_DIR / "Z_gene_names.npy"
Z_HC           = Z_SCORES_DIR / "Z_hc.npy"
Z_HC_NAMES     = Z_SCORES_DIR / "Z_hc_names.npy"

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

MODELING_PARAMS = {
    "ood_percentile":  95,
    "min_samples":     5,
    "z_flag":          3.0,
    "stratify_col":    "Batch_ID",
    "n_splits":        5,
    "det_rate_min":    0.10,   
    "low_det_thr":     0.01,
    "mean_count_min":  2.0,
    "lr_c":            1.0,
    "lr_max_iter":     1000,
    "gsea_gene_sets":  ["KEGG_2021_Human", "GO_Biological_Process_2023", "Reactome_2022"],
    "gsea_fdr_thr":    0.05,
    "gsea_top_n":      30,
    "gsea_perm":       100,
    "gsea_seed":       42,
    "sig_cap_per_theme": 8,
    "emap_sim_thr":    0.50,
}
