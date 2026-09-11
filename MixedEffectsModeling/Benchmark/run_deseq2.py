"""Study-wise DESeq2 benchmark matching 4_disease_scoring.ipynb's sample set.

For every study (Author) that has training HC, each disease phenotype in that study is
tested against that study's HC, twice: without covariates and with the same bias
covariates the normative engine was trained on.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from scipy.sparse import issparse

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import MixedEffectsModeling.config as config

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "DESeq2"
W_DIR = OUT_DIR / "ruvg_W"
RUVG_R = HERE / "ruvg_w.R"
RSCRIPT = Path.home() / "miniconda3/envs/ruvseq_env/bin/Rscript"
MARKER_TSV = config.ROOT / "Data" / "PalangoDB_CellTypeMarkers.tsv"
RUVG_K = [1, 2, 3]
MIN_GROUP = 5
MIN_COUNT_SUM = 10
N_CPUS = 12


def load_cohort():
    adata = sc.read_h5ad(config.H5AD_PATH)
    adata = adata[adata.obs["QC_Passed"] == True]
    adata = adata[adata.obs["Phenotype_Processed"].notna()]
    adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
    adata = adata[adata.obs["broad_protocol_category"] != "Exome-based (EB)"]

    summary = pd.read_csv(config.ENGINE_MIXED_DIR / "training_summary.csv")
    genes = [g for g in summary.loc[summary["ok"], "gene"] if g in adata.var_names]

    obs = adata.obs
    is_hc = obs["Phenotype_Processed"].astype(str) == "Healthy Control"
    bsize = obs.loc[is_hc, config.STRATIFY_COL].astype(str).value_counts()
    small = set(bsize.loc[lambda v: v < config.MIN_HC_BATCH_SIZE].index)
    # HC in small batches are excluded from engine training, so they are excluded here too
    keep = ~(is_hc & obs[config.STRATIFY_COL].astype(str).isin(small))

    adata = adata[keep.values, genes]
    Y = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
    counts = pd.DataFrame(np.round(Y).astype(int), index=adata.obs_names.astype(str), columns=genes)
    meta = adata.obs[["Author", config.STRATIFY_COL, "Phenotype_Processed"] + config.BIAS_COLUMNS].copy()
    meta.columns = ["study", "batch", "phenotype"] + [f"cov{i}" for i in range(len(config.BIAS_COLUMNS))]
    meta.index = counts.index

    # disease samples must match the normative pipeline's OOD-filtered cohort
    # (Z_scores_mixed/sample_meta.csv); HC has no ood_keep entry so it defaults to kept
    ood = pd.read_csv(config.ROOT / "MixedEffectsModeling" / "Z_scores_mixed" / "sample_meta.csv",
                      index_col="sample")["ood_keep"]
    meta["ood_keep"] = meta.index.map(ood).fillna(True)

    return counts, meta, adata.var["GeneName"].astype(str)


def platelet_controls(sym_of):
    """Same control set as EDA/VariousNormalizationMethods_OpenAccess.R: PalangoDB platelet markers."""
    markers = pd.read_csv(MARKER_TSV, sep="\t")
    syms = set(markers.loc[markers["cell type"] == "Platelets", "official gene symbol"])
    return sym_of.index[sym_of.isin(syms)].tolist()


def ruvg_w(sub_counts, controls, k, path):
    if path.exists():
        return pd.read_csv(path, index_col=0)
    with tempfile.TemporaryDirectory() as tmp:
        c_path, g_path = f"{tmp}/counts.csv", f"{tmp}/controls.txt"
        sub_counts.T.to_csv(c_path)
        Path(g_path).write_text("\n".join(controls))
        subprocess.run([str(RSCRIPT), str(RUVG_R), c_path, g_path, str(k), str(path)], check=True)
    return pd.read_csv(path, index_col=0)


def run_one(counts, meta, design, contrast, inference):
    dds = DeseqDataSet(counts=counts, metadata=meta, design=design, inference=inference, quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=contrast, inference=inference, quiet=True)
    stat.summary()
    return stat.results_df.sort_values("padj")


def main():
    counts, meta, sym_of = load_cohort()
    controls = platelet_controls(sym_of)
    cov_terms = " + ".join(f"cov{i}" for i in range(len(config.BIAS_COLUMNS)))
    designs = {"no_covariate": "~condition", "covariate": f"~{cov_terms} + condition"}
    designs.update({f"ruvg_k{k}": "~" + " + ".join([f"W_{i + 1}" for i in range(k)] + ["condition"]) for k in RUVG_K})
    inference = DefaultInference(n_cpus=N_CPUS)
    W_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for study, m in meta.groupby("study", observed=True):
        hc = m.index[m["phenotype"] == "Healthy Control"]
        if len(hc) < MIN_GROUP:
            continue
        for pheno, mp in m[m["phenotype"] != "Healthy Control"].groupby("phenotype", observed=True):
            mp = mp[mp["ood_keep"]]
            if len(mp) < MIN_GROUP:
                continue
            idx = hc.append(mp.index)
            sub_meta = meta.loc[idx].copy()
            sub_meta["condition"] = np.where(sub_meta["phenotype"] == "Healthy Control", "HC", "disease")
            sub_counts = counts.loc[idx]
            sub_counts = sub_counts.loc[:, sub_counts.sum(axis=0) >= MIN_COUNT_SUM]
            # covariates are standardized within the comparison so the fit is conditioned
            # on the same scale the engine's StandardScaler used at training time
            cov_cols = [c for c in sub_meta.columns if c.startswith("cov")]
            sub_meta[cov_cols] = (sub_meta[cov_cols] - sub_meta[cov_cols].mean()) / sub_meta[cov_cols].std()

            tag = f"{study.replace(' ', '_')}__{pheno.replace('/', '-').replace(' ', '_')}"
            # RUVg runs inside the comparison cohort (HC + this disease group only), so W
            # captures the unwanted variation of exactly the contrast it will condition
            W = ruvg_w(sub_counts, controls, max(RUVG_K), W_DIR / f"{tag}_k{max(RUVG_K)}.csv")
            sub_meta[W.columns] = W.loc[sub_meta.index]

            for label, design in designs.items():
                out = OUT_DIR / label
                out.mkdir(parents=True, exist_ok=True)
                path = out / f"{tag}.csv"
                if path.exists():
                    res = pd.read_csv(path, index_col=0)
                else:
                    res = run_one(sub_counts, sub_meta, design, ["condition", "disease", "HC"], inference)
                    res.to_csv(path)
                rows.append(dict(study=study, phenotype=pheno, design=label, n_hc=len(hc), n_disease=len(mp),
                                 n_genes=len(res), shared_batch=bool(set(m.loc[hc, "batch"]) & set(mp["batch"])),
                                 n_sig=int((res["padj"] < 0.05).sum()),
                                 n_sig_up=int(((res["padj"] < 0.05) & (res["log2FoldChange"] > 0)).sum()),
                                 n_sig_down=int(((res["padj"] < 0.05) & (res["log2FoldChange"] < 0)).sum())))
                print(rows[-1], flush=True)

    pd.DataFrame(rows).to_csv(OUT_DIR / "summary.csv", index=False)


if __name__ == "__main__":
    main()
