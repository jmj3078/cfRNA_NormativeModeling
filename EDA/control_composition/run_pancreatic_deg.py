"""Reference DEG runs for the Moore et al. Batch_1 pancreatic cohort.

Five case definitions (PDAC, IPMN, Islet Cell Tumor, Pancreatitis, and the three neoplasms
pooled as one "Pancreatic Cancer" group) against two HC definitions -- Moore Batch_1's own
71 controls, and the full non-Exome HC pool the normative engine trains on -- under four
DESeq2 designs (plain, plus RUVg with k=1,2,3 in the GLM design).

The two HC scopes are the point of running both: the difference between them is itself a
control-composition axis, and it is the axis the normative model handles by construction
(its reference is the full HC pool conditioned on covariates, never a matched subset).

RUVg W is fitted once per (hc_scope, case) sample set with k=3 and sliced for k=1,2, so the
three RUVg designs are nested rather than independently refitted.

Run:  python EDA/control_composition/run_pancreatic_deg.py [--scope both] [--force]
"""
import argparse
import pickle
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from scipy.sparse import issparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
import MixedEffectsModeling.config as mconfig
from MixedEffectsModeling.core.ood_filter import MahalanobisFilter, RangeFilter
from run_control_composition import RUVG_BATCH_R, RSCRIPT, platelet_mask

BATCH = "Moore et al._Batch_1"
CASES = {
    "PDAC": ["PDAC"],
    "IPMN": ["IPMN"],
    "IsletCellTumor": ["Islet Cell Tumor"],
    "Pancreatitis": ["Pancreatitis"],
    "PancreaticNeoplasm": ["PDAC", "IPMN", "Islet Cell Tumor"],
}
HC_SCOPES = ["moore_b1", "full_hc"]
DESIGNS = {"no_covariate": "~condition",
           "ruvg_k1": "~W_1+condition",
           "ruvg_k2": "~W_1+W_2+condition",
           "ruvg_k3": "~W_1+W_2+W_3+condition"}
RUVG_K_MAX = 3
MIN_COUNT_SUM = 10
N_CPUS = 12

OUT_DIR = config.CTRL_COMP_DIR / "pancreatic_deg"
CACHE = OUT_DIR / "pool_cache.pkl"
W_DIR = OUT_DIR / "ruvg_W"
SUMMARY_CSV = OUT_DIR / "summary.csv"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def build_pool():
    if CACHE.exists():
        with open(CACHE, "rb") as f:
            return pickle.load(f)
    log("building pancreatic pool cache from h5ad")
    adata = sc.read_h5ad(mconfig.H5AD_PATH)
    adata = adata[adata.obs["QC_Passed"] == True]
    adata = adata[adata.obs["Phenotype_Processed"].notna()]
    adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
    adata = adata[adata.obs["broad_protocol_category"] != "Exome-based (EB)"]

    pheno = adata.obs["Phenotype_Processed"].astype(str).values
    batch = adata.obs[mconfig.STRATIFY_COL].astype(str).values
    stage = adata.obs["Stage/Condition"].astype(str).values
    is_hc = pheno == "Healthy Control"
    bsize = pd.Series(batch[is_hc]).value_counts()
    small = set(bsize.loc[lambda v: v < mconfig.MIN_HC_BATCH_SIZE].index)
    train_hc = is_hc & ~np.isin(batch, list(small))

    X_all = adata.obs[config.BIAS_COLUMNS].values.astype(float)
    ood = MahalanobisFilter(percentile=config.MODELING_PARAMS["ood_percentile"]).fit(X_all[train_hc])
    rng_f = RangeFilter(n_out_thr=2).fit(X_all[train_hc])
    inlier = ood.mask(X_all) & rng_f.mask(X_all)

    case_stages = sorted({s for v in CASES.values() for s in v})
    is_case = (batch == BATCH) & np.isin(stage, case_stages)
    keep = inlier & (is_case | (is_hc & ~np.isin(batch, list(small))))

    sub = adata[keep]
    raw = sub.layers["Raw"]
    raw = raw.toarray() if issparse(raw) else np.asarray(raw)
    gene_ok = (sub.var["GeneType"] == "protein_coding").values & (raw.sum(axis=0) >= MIN_COUNT_SUM)

    tmm = sub.layers["TMM_log2"]
    tmm = tmm.toarray() if issparse(tmm) else np.asarray(tmm)
    gene_names = sub.var["GeneName"].astype(str).values[gene_ok]
    obs = pd.DataFrame({"sample": sub.obs_names.astype(str).values,
                        "phenotype": pheno[keep], "stage": stage[keep], "batch": batch[keep]})
    data = dict(obs=obs,
                counts=pd.DataFrame(np.round(raw[:, gene_ok]).astype(int),
                                    index=obs["sample"].values,
                                    columns=sub.var_names[gene_ok].astype(str).values),
                tmm=np.asarray(tmm[:, gene_ok], dtype=np.float64),
                genes=sub.var_names[gene_ok].astype(str).values,
                gene_names=gene_names, is_platelet=platelet_mask(gene_names))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(CACHE, "wb") as f:
        pickle.dump(data, f)
    log(f"pool: {len(obs)} samples x {gene_ok.sum()} genes, "
        f"{data['is_platelet'].sum()} platelet control genes")
    log(obs.groupby(['phenotype', 'stage']).size().to_string())
    return data


def subset_index(data, case, scope):
    obs = data["obs"]
    case_mask = obs["stage"].isin(CASES[case]).values & (obs["batch"] == BATCH).values
    hc_mask = (obs["phenotype"] == "Healthy Control").values
    if scope == "moore_b1":
        hc_mask &= (obs["batch"] == BATCH).values
    return np.where(case_mask)[0], np.where(hc_mask)[0]


def fit_ruvg(sid, samples, data):
    out_path = W_DIR / f"{sid}.csv"
    if out_path.exists():
        return
    W_DIR.mkdir(parents=True, exist_ok=True)
    idx = pd.Index(data["obs"]["sample"].values).get_indexer(samples)
    with tempfile.TemporaryDirectory() as tmp:
        tmm_path, ctrl_path, sub_path = f"{tmp}/tmm.csv.gz", f"{tmp}/controls.txt", f"{tmp}/subsets.csv"
        pd.DataFrame(data["tmm"][idx].T, index=data["genes"], columns=samples).to_csv(tmm_path)
        Path(ctrl_path).write_text("\n".join(data["genes"][data["is_platelet"]]))
        pd.DataFrame({"subset_id": sid, "sample": samples}).to_csv(sub_path, index=False)
        subprocess.run([str(RSCRIPT), str(RUVG_BATCH_R), tmm_path, ctrl_path, sub_path,
                        str(RUVG_K_MAX), str(W_DIR)], check=True)


def load_W(sid, samples):
    w = pd.read_csv(W_DIR / f"{sid}.csv").set_index("sample")
    return w.loc[samples, [f"W_{i}" for i in range(1, RUVG_K_MAX + 1)]].values


def fit_deseq2(counts, cond_df, design):
    keep = counts.columns[counts.sum(axis=0) >= MIN_COUNT_SUM]
    inference = DefaultInference(n_cpus=N_CPUS)
    dds = DeseqDataSet(counts=counts[keep], metadata=cond_df, design=design,
                       inference=inference, quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "case", "HC"], inference=inference, quiet=True)
    stat.summary()
    return stat.results_df


def run_one(data, scope, case, force=False):
    case_idx, hc_idx = subset_index(data, case, scope)
    sid = f"{scope}__{case}"
    samples = np.concatenate([data["obs"]["sample"].values[case_idx],
                              data["obs"]["sample"].values[hc_idx]])
    condition = np.array(["case"] * len(case_idx) + ["HC"] * len(hc_idx))
    counts = data["counts"].loc[samples]

    fit_ruvg(sid, samples, data)
    W = load_W(sid, samples)

    rows = []
    stat_dir = OUT_DIR / sid
    stat_dir.mkdir(parents=True, exist_ok=True)
    for name, design in DESIGNS.items():
        out_path = stat_dir / f"{name}.csv.gz"
        if out_path.exists() and not force:
            res = pd.read_csv(out_path, index_col=0)
        else:
            cond_df = pd.DataFrame({"condition": condition}, index=samples)
            if name.startswith("ruvg_k"):
                for i in range(int(name.rsplit("_k", 1)[1])):
                    cond_df[f"W_{i + 1}"] = W[:, i]
            t0 = time.time()
            res = fit_deseq2(counts, cond_df, design)
            res.to_csv(out_path)
            log(f"{sid}/{name}: {time.time() - t0:.0f}s")
        padj = res["padj"].fillna(1.0)
        rows.append(dict(hc_scope=scope, case=case, design=name,
                         n_case=len(case_idx), n_hc=len(hc_idx),
                         n_sig_q05=int((padj < 0.05).sum()), n_sig_q10=int((padj < 0.10).sum()),
                         n_up=int(((padj < 0.05) & (res["log2FoldChange"] > 0)).sum()),
                         n_down=int(((padj < 0.05) & (res["log2FoldChange"] < 0)).sum())))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scope", default="both", choices=HC_SCOPES + ["both"])
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    scopes = HC_SCOPES if args.scope == "both" else [args.scope]
    data = build_pool()

    rows = []
    for scope in scopes:
        for case in CASES:
            rows += run_one(data, scope, case, force=args.force)
            pd.DataFrame(rows).to_csv(SUMMARY_CSV, index=False)
    log(f"done -> {SUMMARY_CSV}")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
