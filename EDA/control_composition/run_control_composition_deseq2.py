"""DESeq2 arm of the control-composition sensitivity experiment (Moore et al. Batch_1).

Same subsets, same tags, same seed as run_control_composition.py -- reuses the RUVg W
factors that script already fitted and cached (no refit here). Two designs per subset:
plain (~condition) and RUVg-covariate (~W_1+W_2+condition), the latter putting W into the
GLM design the way DESeq2 actually consumes unwanted variation (residualizing raw counts
is not meaningful for a count model).

Null draw count defaults lower than the Welch-t arm (30 vs 200) -- DESeq2 is ~10s/fit vs
~instant for Welch-t, so the same 1266-subset null is a multi-hour run. First N draws of a
fixed seed are a prefix of any larger run, so this reuses cached W without recomputing.

Run:  python EDA/control_composition/run_control_composition_deseq2.py [--n-null 30] [--force]
"""
import argparse
import sys
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
from run_control_composition import (RUVG_K, SEED, build_cache, enumerate_groups, load_W, log,
                                     pair_metrics, cohens_d, bias_axes, append_row)

DESIGNS = {"no_covariate": "~condition",
           "ruvg_k1": "~W_1+condition",
           "ruvg_k2": "~W_1+W_2+condition",
           "ruvg_k3": "~W_1+W_2+W_3+condition"}
MIN_COUNT_SUM = 10
N_CPUS = 16

RAW_CACHE = config.CTRL_COMP_DIR / "moore_b1_raw_counts.pkl"
RESULT_CSV = config.CTRL_COMP_DIR / "deseq2_jaccard_results.csv"
NULL_CSV = config.CTRL_COMP_DIR / "deseq2_null_distribution.csv"


def load_raw_counts(data):
    """Raw counts aligned to the same sample/gene order as the shared cache."""
    if RAW_CACHE.exists():
        return pd.read_pickle(RAW_CACHE)
    log("loading raw counts for DESeq2 (aligned to existing sample/gene cache)")
    adata = sc.read_h5ad(mconfig.H5AD_PATH)
    sub = adata[data["obs"]["sample"].values, data["genes"]]
    raw = sub.layers["Raw"]
    raw = raw.toarray() if issparse(raw) else np.asarray(raw)
    counts = pd.DataFrame(np.round(raw).astype(int), index=sub.obs_names.astype(str), columns=data["genes"])
    counts.to_pickle(RAW_CACHE)
    return counts


def fit_one(counts, cond_df, design):
    keep = counts.columns[counts.sum(axis=0) >= MIN_COUNT_SUM]
    inference = DefaultInference(n_cpus=N_CPUS)
    dds = DeseqDataSet(counts=counts[keep], metadata=cond_df, design=design, inference=inference, quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "disease", "HC"], inference=inference, quiet=True)
    stat.summary()
    out = pd.Series(0.0, index=counts.columns)
    out.loc[stat.results_df.index] = stat.results_df["stat"].fillna(0.0).values
    return out.values


def process_group(data, counts, g, layers=None):
    layers = list(DESIGNS) if layers is None else list(layers)
    samples = data["obs"]["sample"].values
    axes = bias_axes(data["obs"])
    n_case = len(g["case"])
    stats = {name: [] for name in layers}
    ds = []
    for t, ctrl in enumerate(g["strata"]):
        sid = f"{g['tag']}__T{t}"
        idx = np.concatenate([g["case"], ctrl])
        W = load_W(sid, samples[idx])
        sub_counts = counts.iloc[idx]
        condition = np.array(["disease"] * n_case + ["HC"] * len(ctrl))

        stat_dir = config.CTRL_COMP_DESEQ2_DIR / g["tag"]
        stat_dir.mkdir(parents=True, exist_ok=True)
        for name in layers:
            out_path = stat_dir / f"T{t}_{name}.csv.gz"
            if out_path.exists():
                s = pd.read_csv(out_path, index_col=0)["stat"].values
            else:
                cond_df = pd.DataFrame({"condition": condition}, index=sub_counts.index)
                if name.startswith("ruvg_k"):
                    for i in range(int(name.rsplit("_k", 1)[1])):
                        cond_df[f"W_{i + 1}"] = W[:, i]
                s = fit_one(sub_counts, cond_df, DESIGNS[name])
                pd.DataFrame({"stat": s}, index=data["genes"]).to_csv(out_path)
            stats[name].append(s)
        if g["split"] == "tertile":
            v = axes[g["axis"]]
            ds.append(cohens_d(v[g["case"]], v[ctrl]))

    out_path = RESULT_CSV if g["split"] == "tertile" else NULL_CSV
    for name in layers:
        acc = [pair_metrics(np.asarray(stats[name][i]), np.asarray(stats[name][j]))
               for i in range(3) for j in range(i + 1, 3)]
        row = pd.DataFrame(acc).mean().to_dict()
        row.update(disease=g["disease"], split=g["split"], axis=g["axis"], draw=g["draw"],
                   tag=g["tag"], layer=name, n_case=n_case,
                   n_ctrl=float(np.mean([len(s) for s in g["strata"]])),
                   delta_d=float(max(ds) - min(ds)) if ds else np.nan)
        append_row(out_path, row)


def _fit(counts, cond_df, design):
    keep = counts.columns[counts.sum(axis=0) >= MIN_COUNT_SUM]
    inference = DefaultInference(n_cpus=N_CPUS)
    dds = DeseqDataSet(counts=counts[keep], metadata=cond_df, design=design, inference=inference, quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "disease", "HC"], inference=inference, quiet=True)
    stat.summary()
    out = pd.Series(0.0, index=counts.columns)
    out.loc[stat.results_df.index] = stat.results_df["stat"].fillna(0.0).values
    return out.values


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-null", type=int, default=30)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    config.CTRL_COMP_DIR.mkdir(parents=True, exist_ok=True)
    if args.force:
        for p in (RESULT_CSV, NULL_CSV):
            p.unlink(missing_ok=True)

    data = build_cache()
    counts = load_raw_counts(data)
    # enumerate_groups shares one rng across diseases: disease 2's null draws start where
    # disease 1's left off, so its tag<->sample mapping depends on the FULL n_null used by
    # whichever run built the rng sequence first (here, run_control_composition.py's 200).
    # Generate that same full sequence and only pick the first --n-null draws per disease,
    # so tags always resolve to the W files already cached under those tags.
    groups_full = enumerate_groups(data, n_null=200, seed=SEED)
    groups = [g for g in groups_full if g["split"] == "tertile" or g["draw"] < args.n_null]
    log(f"deseq2: {len(groups)} comparison groups ({len(DESIGNS)} designs each)")

    done = set()
    for path in (RESULT_CSV, NULL_CSV):
        if path.exists():
            d = pd.read_csv(path, usecols=["tag", "layer"])
            done |= set(zip(d["tag"], d["layer"]))
    todo = [(g, [n for n in DESIGNS if (g["tag"], n) not in done]) for g in groups]
    todo = [(g, miss) for g, miss in todo if miss]
    log(f"deseq2: {len(todo)} groups to run, "
        f"{sum(len(m) for _, m in todo)} (group, design) fits missing")

    t0 = time.time()
    for i, (g, miss) in enumerate(todo, 1):
        process_group(data, counts, g, miss)
        if i % 5 == 0 or i == len(todo):
            rate = (time.time() - t0) / i
            log(f"deseq2: {i}/{len(todo)} groups ({rate:.1f}s/group, "
                f"eta {(len(todo) - i) * rate / 60:.0f} min)")
    log("deseq2 done")


if __name__ == "__main__":
    main()
