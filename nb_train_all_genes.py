#!/usr/bin/env python
"""Train NB-GLM / Bayesian NB-GLM / Laplace NB-GP on Healthy Control samples
for every candidate protein-coding gene, and save the fitted models to
NB_GP/saved_models.

Candidate genes are those passing the same filter used for stratified gene
selection in the notebook: det_rate_hc >= det_rate_min AND mean_count_hc >= mean_count_min.

Resumable: a gene is skipped if its 3 model files already exist in SAVE_DIR
(unless --no-resume is given).
"""

import argparse
import csv
import os
import pickle
import time
from pathlib import Path

import numpy as np
import scanpy as sc
import torch
from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler

from nb_models import NBGLM, train_bayesian_nbglm, train_laplace_nbgp

torch.manual_seed(42)
np.random.seed(42)

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "OpenAccess_nfcore"
H5AD_PATH = DATA_DIR / "Merged_Processed_AnnData_with_Batch_Biases_QC_Status.h5ad"
SAVE_DIR = BASE_DIR / "NB_GP" / "saved_models"
META_PATH = BASE_DIR / "NB_GP" / "training_meta.csv"

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

DET_RATE_MIN = 0.1
MEAN_COUNT_MIN = 2.0
BGLM_MAX_EPOCHS = 150
LAP_MAX_EPOCHS = 150

META_FIELDS = [
    "gene", "det_rate_hc", "hc_mean_count",
    "time_glm_s", "theta_glm",
    "time_bglm_s", "theta_bglm",
    "time_lap_s", "theta_lap",
]

def load_data():
    adata = sc.read_h5ad(H5AD_PATH)
    adata = adata[adata.obs["QC_Passed"] == True]
    adata = adata[adata.obs["Phenotype_Processed"].notna()]
    adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
    phenotypes = adata.obs["Phenotype_Processed"].astype(str).values
    is_hc = phenotypes == "Healthy Control"

    X_raw = adata.obs[BIAS_COLUMNS].values.astype(np.float32)
    scaler_X = StandardScaler()
    X_hc_scaled = scaler_X.fit_transform(X_raw[is_hc])

    Y_raw = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
    Y_raw = np.round(Y_raw).astype(np.float32)

    is_pc = (adata.var["GeneType"] == "protein_coding").values
    pc_gene_names = adata.var_names[is_pc].tolist()
    pc_indices = np.where(is_pc)[0]

    return adata, is_hc, X_hc_scaled, Y_raw, pc_gene_names, pc_indices


def select_candidate_genes(Y_raw, is_hc, pc_gene_names, pc_indices, det_rate_min, mean_count_min):
    Y_hc = Y_raw[is_hc][:, pc_indices]
    det_r = (Y_hc > 0).mean(axis=0)
    mean_c = Y_hc.mean(axis=0)
    cand = (det_r >= det_rate_min) & (mean_c >= mean_count_min)
    return np.array(pc_gene_names)[cand].tolist(), pc_indices[cand]


def already_done(g_name):
    return all(
        (SAVE_DIR / fname).exists()
        for fname in (f"{g_name}_glm.pkl", f"{g_name}_bglm.pt", f"{g_name}_laplace.pt")
    )


def fmt_theta(v):
    return f"{v:.2f}" if v is not None else "NA"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--det-rate-min", type=float, default=DET_RATE_MIN)
    parser.add_argument("--mean-count-min", type=float, default=MEAN_COUNT_MIN)
    parser.add_argument("--bglm-epochs", type=int, default=BGLM_MAX_EPOCHS)
    parser.add_argument("--lap-epochs", type=int, default=LAP_MAX_EPOCHS)
    parser.add_argument("--device", default=None, help="'cuda' or 'cpu' (default: auto)")
    parser.add_argument("--no-resume", action="store_true",
                         help="retrain genes even if saved model files already exist")
    parser.add_argument("--limit", type=int, default=None,
                         help="only process the first N candidate genes (for testing)")
    args = parser.parse_args()

    device = torch.device(args.device) if args.device else torch.device("cuda" if torch.cuda.is_available() else "cpu")
    SAVE_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Device   : {device}")
    print(f"H5AD     : {H5AD_PATH}")
    print(f"Save dir : {SAVE_DIR}")

    adata, is_hc, X_hc_scaled, Y_raw, pc_gene_names, pc_indices = load_data()
    gene_names, gene_indices = select_candidate_genes(
        Y_raw, is_hc, pc_gene_names, pc_indices, args.det_rate_min, args.mean_count_min
    )

    print(f"HC samples       : {is_hc.sum()} / {len(adata)}")
    print(f"Protein-coding   : {len(pc_gene_names)}")
    print(f"Candidate genes  : {len(gene_names)} "
          f"(det_rate>={args.det_rate_min}, mean_count>={args.mean_count_min})")

    if args.limit is not None:
        gene_names, gene_indices = gene_names[:args.limit], gene_indices[:args.limit]
        print(f"Limited to first {len(gene_names)} candidate gene(s) for this run")

    # BGLM (small p, many tiny ops) is faster on CPU; Laplace GP (NxN Cholesky)
    # is faster on GPU. Keep both tensor copies and dispatch per step.
    X_cpu_t = torch.tensor(X_hc_scaled, dtype=torch.float32)
    X_tr_t = X_cpu_t.to(device)

    # Append metadata one gene at a time so progress survives interruption.
    write_header = args.no_resume or not META_PATH.exists()
    meta_file = open(META_PATH, "w" if args.no_resume else "a", newline="")
    meta_writer = csv.DictWriter(meta_file, fieldnames=META_FIELDS)
    if write_header:
        meta_writer.writeheader()
        meta_file.flush()

    n_done = 0
    n_skipped = 0
    t_start = time.perf_counter()

    for i, (g_name, g_idx) in enumerate(zip(gene_names, gene_indices)):
        if not args.no_resume and already_done(g_name):
            n_skipped += 1
            continue

        y_hc = Y_raw[is_hc, g_idx]
        y_cpu_t = torch.tensor(y_hc, dtype=torch.float32)
        y_tr_t = y_cpu_t.to(device)

        det_rate_hc = float((y_hc > 0).mean())
        hc_mean_count = float(y_hc.mean())
        meta = {"gene": g_name, "det_rate_hc": det_rate_hc, "hc_mean_count": hc_mean_count}

        # 1. GLM (statsmodels)
        t0 = time.perf_counter()
        glm = NBGLM().fit(X_hc_scaled, y_hc)
        meta["time_glm_s"] = time.perf_counter() - t0
        meta["theta_glm"] = glm.theta_

        with open(SAVE_DIR / f"{g_name}_glm.pkl", "wb") as f:
            pickle.dump({"model": glm, "theta": glm.theta_}, f)
        del glm

        # 2. Bayesian NB-GLM (small matrices -> CPU is faster than GPU)
        t0 = time.perf_counter()
        bglm, beta_map, neg_H = train_bayesian_nbglm(X_cpu_t, y_cpu_t, max_epochs=args.bglm_epochs)
        meta["time_bglm_s"] = time.perf_counter() - t0
        meta["theta_bglm"] = bglm.theta.item()

        torch.save({
            "state_dict": bglm.state_dict(),
            "beta_map": beta_map.cpu(),
            "neg_H": neg_H.cpu(),
        }, SAVE_DIR / f"{g_name}_bglm.pt")
        del bglm, beta_map, neg_H

        # 3. Laplace NB-GP
        t0 = time.perf_counter()
        lap_model, f_map = train_laplace_nbgp(X_tr_t, y_tr_t, max_epochs=args.lap_epochs)
        meta["time_lap_s"] = time.perf_counter() - t0
        meta["theta_lap"] = lap_model.theta.item()

        torch.save({
            "state_dict": lap_model.state_dict(),
            "f_map": f_map.cpu(),
        }, SAVE_DIR / f"{g_name}_laplace.pt")
        del lap_model, f_map

        meta_writer.writerow(meta)
        meta_file.flush()
        n_done += 1

        if device.type == "cuda":
            torch.cuda.empty_cache()

        print(
            f"[{i + 1:4d}/{len(gene_names)}] {g_name:<22s} "
            f"det={det_rate_hc:.2f} mean={hc_mean_count:7.1f}  "
            f"theta: glm={fmt_theta(meta['theta_glm'])} bglm={meta['theta_bglm']:.2f} lap={meta['theta_lap']:.2f}  |  "
            f"t: {meta['time_glm_s']:.1f}+{meta['time_bglm_s']:.1f}+{meta['time_lap_s']:.1f}s"
        )

    meta_file.close()

    elapsed = time.perf_counter() - t_start
    print(f"\nDone. Trained {n_done} genes, skipped {n_skipped} already-saved gene(s) "
          f"in {elapsed:.1f}s -> {SAVE_DIR}")
    print(f"Saved training metadata -> {META_PATH}")


if __name__ == "__main__":
    main()
