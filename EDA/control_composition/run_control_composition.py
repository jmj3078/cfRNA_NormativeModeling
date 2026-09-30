"""Control-composition sensitivity of group-wise biomarker selection (Moore et al. Batch_1).

Case group fixed; the HC pool split three ways at random, each stratum giving its own marker
list. Agreement between those lists is the estimand: how much of a "biomarker" set is decided
by which healthy controls happened to be drawn.

Bias-axis (tertile) stratified splits were removed on 2026-09-30 -- within Moore's old
Batch_1 the bias axes tracked the hidden pre-analytical sub-batches (Cramer's V 0.35-0.58 on
8 of 10 axes), so a tertile split was partly a batch split. Random splits only.

RUVg factors are re-estimated inside every single comparison (case + that one stratum),
because the stored RUVg_Platelet_* / Proposed_Full_* layers were fitted on the full cohort
and carry information no single-comparison pipeline could have had. Per-sample scalings
(CPM/TPM/TMM) are unaffected and are read from the stored layers.

Run:  python EDA/control_composition/run_control_composition.py [--n-null 200] [--force]
Every comparison writes its own W / t-statistic file as it finishes,
so the run is resumable and inspectable while still going.
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
from scipy.sparse import issparse
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import config
import MixedEffectsModeling.config as mconfig
from MixedEffectsModeling.core.ood_filter import MahalanobisFilter, RangeFilter

HERE = Path(__file__).resolve().parent
RUVG_BATCH_R = HERE / "ruvg_batch.R"
RSCRIPT = Path.home() / "miniconda3/envs/ruvseq_env/bin/Rscript"
MARKER_TSV = ROOT / "Data" / "PalangoDB_CellTypeMarkers.tsv"

BATCH = "Moore et al._Batch_1"
DISEASES = ["Pancreatic Cancer", "Pancreatitis"]
STATIC_LAYERS = ["CPM_log1p", "TPM_log2", "TMM_log2"]
RUVG_K = [1, 2, 3]
DYNAMIC_LAYERS = [f"RUVg_Platelet_k{k}" for k in RUVG_K] + [f"Proposed_Full_k{k}" for k in RUVG_K]
BASE_LAYERS = STATIC_LAYERS + ["EDA_Full_All"]
K_GRID = [25, 50, 100, 200, 500]
MIN_COUNT_SUM = 10
SEED = 42

CACHE = config.CTRL_COMP_DIR / "moore_b1_cache.pkl"
SUBSET_CSV = config.CTRL_COMP_DIR / "subsets.csv"
NULL_CSV = config.CTRL_COMP_DIR / "null_distribution.csv"
LOG_PATH = config.CTRL_COMP_DIR / "progress.log"


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a") as f:
        f.write(line + "\n")


# --------------------------------------------------------------------------- data
def build_cache():
    if CACHE.exists():
        with open(CACHE, "rb") as f:
            return pickle.load(f)
    log("building Moore Batch_1 cache from h5ad")
    adata = sc.read_h5ad(mconfig.H5AD_PATH)
    adata = adata[adata.obs["QC_Passed"] == True]
    adata = adata[adata.obs["Phenotype_Processed"].notna()]
    adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
    adata = adata[adata.obs["broad_protocol_category"] != "Exome-based (EB)"]

    pheno = adata.obs["Phenotype_Processed"].astype(str).values
    batch = adata.obs[mconfig.STRATIFY_COL].astype(str).values
    is_hc = pheno == "Healthy Control"
    bsize = pd.Series(batch[is_hc]).value_counts()
    small = set(bsize.loc[lambda v: v < mconfig.MIN_HC_BATCH_SIZE].index)
    train_hc = is_hc & ~np.isin(batch, list(small))

    X_all = adata.obs[config.BIAS_COLUMNS].values.astype(float)
    ood = MahalanobisFilter(percentile=config.MODELING_PARAMS["ood_percentile"]).fit(X_all[train_hc])
    rng_f = RangeFilter(n_out_thr=2).fit(X_all[train_hc])
    inlier = ood.mask(X_all) & rng_f.mask(X_all)

    keep = (batch == BATCH) & inlier & np.isin(pheno, ["Healthy Control"] + DISEASES)
    sub = adata[keep]
    raw = sub.layers["Raw"]
    raw = raw.toarray() if issparse(raw) else np.asarray(raw)
    gene_ok = (sub.var["GeneType"] == "protein_coding").values & (raw.sum(axis=0) >= MIN_COUNT_SUM)

    layers = {}
    for name in BASE_LAYERS:
        L = sub.layers[name]
        L = L.toarray() if issparse(L) else np.asarray(L)
        layers[name] = np.asarray(L[:, gene_ok], dtype=np.float64)

    gene_names = sub.var["GeneName"].astype(str).values[gene_ok]
    obs = pd.DataFrame({"sample": sub.obs_names.astype(str).values, "phenotype": pheno[keep]})
    obs[config.BIAS_COLUMNS] = X_all[keep]
    data = dict(obs=obs, layers=layers, genes=sub.var_names[gene_ok].astype(str).values,
                gene_names=gene_names, is_platelet=platelet_mask(gene_names))
    config.CTRL_COMP_DIR.mkdir(parents=True, exist_ok=True)
    with open(CACHE, "wb") as f:
        pickle.dump(data, f)
    log(f"cache: {len(obs)} samples x {gene_ok.sum()} genes, "
        f"{data['is_platelet'].sum()} platelet control genes")
    log(obs["phenotype"].value_counts().to_string())
    return data


def platelet_mask(gene_names):
    """Same control set as EDA/control_composition/VariousNormalizationMethods_OpenAccess.R."""
    markers = pd.read_csv(MARKER_TSV, sep="\t")
    syms = set(markers.loc[markers["cell type"] == "Platelets", "official gene symbol"])
    return np.isin(gene_names, list(syms))


# --------------------------------------------------------------------------- splits
def slug(s):
    return "".join(c if c.isalnum() else "_" for c in s).strip("_")


def enumerate_groups(data, n_null, seed=SEED):
    """Every comparison group this analysis runs. One group = one case cohort split three
    ways at random; the three (case vs stratum) comparisons inside it are what get compared.

    Bias-axis (tertile) stratified splits were removed: within Moore's old Batch_1 the bias
    axes tracked the hidden pre-analytical sub-batches, so a tertile split was partly a batch
    split rather than a control-composition contrast. Random splits only from here on."""
    obs = data["obs"]
    hc_idx = np.where(obs["phenotype"].values == "Healthy Control")[0]
    n = len(hc_idx)
    cuts = [0, n // 3, 2 * n // 3, n]
    rng = np.random.default_rng(seed)
    groups = []
    for disease in DISEASES:
        case_idx = np.where(obs["phenotype"].values == disease)[0]
        for b in range(n_null):
            perm = rng.permutation(hc_idx)
            groups.append(dict(disease=disease, split="random", draw=b,
                               tag=f"{slug(disease)}__null_{b:04d}",
                               case=case_idx, strata=[perm[cuts[t]:cuts[t + 1]] for t in range(3)]))
    return groups


def subset_table(data, groups):
    samples = data["obs"]["sample"].values
    rows = []
    for g in groups:
        for t, ctrl in enumerate(g["strata"]):
            sid = f"{g['tag']}__T{t}"
            for s in samples[np.concatenate([g["case"], ctrl])]:
                rows.append(dict(subset_id=sid, sample=s))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- RUVg
def run_ruvg(data, tbl, k=max(RUVG_K)):
    """RUVSeq::RUVg per comparison. All subsets go through one R process -- a per-subset
    Rscript launch does not scale to the null -- but R writes each W file as it finishes,
    so the stage is resumable. TMM is a per-sample scaling and is not refitted; only the
    RUV factors are."""
    config.CTRL_COMP_W_DIR.mkdir(parents=True, exist_ok=True)
    done = {p.stem for p in config.CTRL_COMP_W_DIR.glob("*.csv")}
    todo = tbl[~tbl["subset_id"].isin(done)]
    if todo.empty:
        log(f"ruvg: all {tbl['subset_id'].nunique()} subsets already present")
        return
    log(f"ruvg: {todo['subset_id'].nunique()} subsets to fit "
        f"({len(done)} cached), k={k}")
    with tempfile.TemporaryDirectory() as tmp:
        tmm_path, ctrl_path, sub_path = f"{tmp}/tmm.csv.gz", f"{tmp}/controls.txt", f"{tmp}/subsets.csv"
        pd.DataFrame(data["layers"]["TMM_log2"].T, index=data["genes"],
                     columns=data["obs"]["sample"].values).to_csv(tmm_path)
        Path(ctrl_path).write_text("\n".join(data["genes"][data["is_platelet"]]))
        todo.to_csv(sub_path, index=False)
        subprocess.run([str(RSCRIPT), str(RUVG_BATCH_R), tmm_path, ctrl_path, sub_path,
                        str(k), str(config.CTRL_COMP_W_DIR)], check=True)


def load_W(subset_id, samples):
    w = pd.read_csv(config.CTRL_COMP_W_DIR / f"{subset_id}.csv").set_index("sample")
    return w.loc[samples, [f"W_{i}" for i in RUVG_K]].values


def residualize(Y, W, intercept):
    """RUVSeq normalizedCounts (intercept=False) and limma::removeBatchEffect
    (intercept=True) are both OLS removal of the W columns."""
    D = np.column_stack([np.ones(len(W)), W]) if intercept else W
    beta = np.linalg.lstsq(D, Y, rcond=None)[0]
    return Y - (D[:, 1:] @ beta[1:] if intercept else D @ beta)


def verify_residualize(data):
    """Self-check: the python residualization must reproduce RUVSeq's normalizedCounts on
    the subset ruvg_batch.R dumped. Fails loudly rather than silently diverging."""
    check = config.CTRL_COMP_W_DIR / "_check_normalizedCounts.csv.gz"
    if not check.exists():
        log("verify: no R reference dump found, skipped")
        return
    ref = pd.read_csv(check, index_col=0)
    samples = data["obs"]["sample"].values
    sid = sorted(p.stem for p in config.CTRL_COMP_W_DIR.glob("*.csv") if not p.stem.startswith("_"))
    for s in sid:
        w = pd.read_csv(config.CTRL_COMP_W_DIR / f"{s}.csv")
        if set(w["sample"]) == set(ref.columns):
            idx = np.array([np.where(samples == c)[0][0] for c in ref.columns])
            W = w.set_index("sample").loc[ref.columns, [f"W_{i}" for i in RUVG_K]].values
            got = residualize(data["layers"]["TMM_log2"][idx], W, intercept=False)
            dev = float(np.abs(got.T - ref.values).max())
            assert dev < 1e-6, f"python residualization deviates from RUVSeq by {dev:.3g}"
            log(f"verify: python residualization matches RUVSeq (max dev {dev:.2e})")
            return
    log("verify: could not match the R reference dump to a subset, skipped")


# --------------------------------------------------------------------------- metrics
def welch_t(A, B):
    ma, mb = A.mean(0), B.mean(0)
    se = np.sqrt(A.var(0, ddof=1) / A.shape[0] + B.var(0, ddof=1) / B.shape[0])
    return np.divide(ma - mb, se, out=np.zeros_like(ma), where=se > 0)


def comparison_layers(data, case_idx, ctrl_idx, W):
    idx = np.concatenate([case_idx, ctrl_idx])
    out = {name: data["layers"][name][idx] for name in STATIC_LAYERS}
    tmm, eda = data["layers"]["TMM_log2"][idx], data["layers"]["EDA_Full_All"][idx]
    for k in RUVG_K:
        out[f"RUVg_Platelet_k{k}"] = residualize(tmm, W[:, :k], intercept=False)
        out[f"Proposed_Full_k{k}"] = residualize(eda, W[:, :k], intercept=True)
    return out


def pair_metrics(ti, tj, k_grid=K_GRID):
    out = {}
    oi, oj = np.argsort(-np.abs(ti)), np.argsort(-np.abs(tj))
    for k in k_grid:
        a, b = set(oi[:k]), set(oj[:k])
        inter = np.array(sorted(a & b), dtype=int)
        out[f"jaccard_k{k}"] = len(inter) / len(a | b)
        out[f"signflip_k{k}"] = float((np.sign(ti[inter]) != np.sign(tj[inter])).mean()) if len(inter) else np.nan
    out["spearman"] = float(np.corrcoef(rankdata(ti), rankdata(tj))[0, 1])
    return out


def append_row(path, row):
    df = pd.DataFrame([row])
    df.to_csv(path, mode="a", header=not path.exists(), index=False)


# --------------------------------------------------------------------------- driver
def process_group(data, g):
    """One comparison group -> per-stratum t-statistics on disk, plus one metrics row
    averaged over the three pairwise list comparisons."""
    samples = data["obs"]["sample"].values
    stats = []
    for t, ctrl in enumerate(g["strata"]):
        sid = f"{g['tag']}__T{t}"
        idx = np.concatenate([g["case"], ctrl])
        W = load_W(sid, samples[idx])
        mats = comparison_layers(data, g["case"], ctrl, W)
        n_case = len(g["case"])
        tt = {name: welch_t(M[:n_case], M[n_case:]) for name, M in mats.items()}
        stats.append(tt)

        stat_dir = config.CTRL_COMP_STAT_DIR / g["tag"]
        stat_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(tt, index=data["genes"]).assign(GeneName=data["gene_names"]).to_csv(
            stat_dir / f"T{t}_welch_t.csv.gz")

    for layer in STATIC_LAYERS + DYNAMIC_LAYERS:
        acc = [pair_metrics(stats[i][layer], stats[j][layer])
               for i in range(3) for j in range(i + 1, 3)]
        row = pd.DataFrame(acc).mean().to_dict()
        row.update(disease=g["disease"], split=g["split"], draw=g["draw"],
                   tag=g["tag"], layer=layer, n_case=len(g["case"]),
                   n_ctrl=float(np.mean([len(s) for s in g["strata"]])))
        append_row(NULL_CSV, row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-null", type=int, default=200)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    config.CTRL_COMP_DIR.mkdir(parents=True, exist_ok=True)
    if args.force:
        NULL_CSV.unlink(missing_ok=True)

    data = build_cache()
    groups = enumerate_groups(data, n_null=args.n_null)
    tbl = subset_table(data, groups)
    tbl.to_csv(SUBSET_CSV, index=False)
    log(f"{len(groups)} comparison groups, {tbl['subset_id'].nunique()} comparisons")

    run_ruvg(data, tbl)
    verify_residualize(data)

    done = set()
    if NULL_CSV.exists():
        done = set(pd.read_csv(NULL_CSV, usecols=["tag"])["tag"].unique())
    todo = [g for g in groups if g["tag"] not in done]
    log(f"metrics: {len(todo)} groups to run ({len(done)} already done)")

    t0 = time.time()
    for i, g in enumerate(todo, 1):
        process_group(data, g)
        if i % 20 == 0 or i == len(todo):
            rate = (time.time() - t0) / i
            log(f"metrics: {i}/{len(todo)} groups ({rate:.1f}s/group, "
                f"eta {(len(todo) - i) * rate / 60:.0f} min)")
    log("done")


if __name__ == "__main__":
    main()
