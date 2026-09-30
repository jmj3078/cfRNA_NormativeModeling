"""Pathway-level reproducibility of random (null) control splits.

For each random HC split of Moore et al. Batch_1, runs prerank GSEA on the cached DESeq2
stat files (T0/T1/T2_{design}.csv.gz from run_control_composition_deseq2.py) and reports how
much the FDR<0.05 pathway sets agree across the three strata. No refitting, GSEA only.

Same prerank machinery (gsea_prerank, housekeeping-excluded KEGG+Reactome library) as
MixedEffectsModeling's own group_level_pathway_gsea, for direct comparability.

The bias-axis (tertile) stratified counterpart was dropped on 2026-09-30: the bias axes
tracked hidden pre-analytical sub-batches inside the old Batch_1, so a tertile split was
partly a batch split.

Run: python EDA/control_composition/run_gsea_split_jaccard_null.py [--n-null 10]
"""
import argparse
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from pathway_lib import ensg_to_symbol, gsea_prerank, load_pathway_library, load_symbol_vocab

DESIGNS = ["no_covariate", "ruvg_k1", "ruvg_k2", "ruvg_k3"]
DISEASE = "Pancreatic Cancer"
def out_csv(scope):
    return config.CTRL_COMP_DIR / "gsea_split_jaccard" / f"summary_null_{scope}.csv"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def jaccard(a, b):
    a, b = set(a), set(b)
    u = a | b
    return len(a & b) / len(u) if u else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-null", type=int, default=10)
    ap.add_argument("--draws", default=None, help="comma-separated draw indices (shard); default 0..n_null-1")
    ap.add_argument("--designs", default=None, help="comma-separated designs; default all four")
    ap.add_argument("--out", default=None, help="output CSV (shard); merged into the scope summary later")
    ap.add_argument("--scope", default=config.CTRL_COMP_DEFAULT_SCOPE,
                    choices=sorted(config.CTRL_COMP_SCOPES))
    args = ap.parse_args()
    draws = [int(d) for d in args.draws.split(",")] if args.draws else list(range(args.n_null))
    designs = args.designs.split(",") if args.designs else DESIGNS
    summary = out_csv(args.scope)
    out_path = Path(args.out) if args.out else summary

    universe_syms, sym2idx, col2sym = load_symbol_vocab(None)
    terms, M = load_pathway_library()
    log(f"library: {len(terms)} terms (housekeeping-excluded, matches pathway_lib)")

    sym_of = ensg_to_symbol()
    sym_of.index = sym_of.index.str.split(".").str[0]

    rows = []
    out_path.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if summary.exists():
        prev = pd.read_csv(summary)
        done = set(prev["tag"] + "/" + prev["method"])

    for draw in draws:
        tag = f"{DISEASE.replace(' ', '_')}__null_{draw:04d}"
        stat_dir = config.CTRL_COMP_DESEQ2_DIR / args.scope / tag
        if not stat_dir.exists():
            log(f"skip {tag}: no cached DESeq2 stats")
            continue
        for design in designs:
            key = f"{tag}/deseq2__{design}"
            if key in done:
                continue
            t0 = time.time()
            n_sig, sig_terms = [], []
            for t in range(config.ctrl_comp_n_splits(args.scope)):
                stat = pd.read_csv(stat_dir / f"T{t}_{design}.csv.gz", index_col=0)["stat"]
                stat.index = stat.index.str.split(".").str[0]
                syms = sym_of.reindex(stat.index)
                rnk = pd.Series(stat.values, index=syms.values)
                rnk = rnk[pd.notna(rnk.index)].groupby(level=0).mean().reindex(universe_syms).values
                res2d = gsea_prerank(rnk, terms, M, universe_syms, n_perm=1000, seed=42 + t)
                sig = set(res2d.loc[res2d["FDR q-val"] < 0.05, "Term"])
                n_sig.append(len(sig))
                sig_terms.append(sig)
            pairs = list(combinations(range(len(sig_terms)), 2))
            jaccs = [jaccard(sig_terms[i], sig_terms[j]) for i, j in pairs]
            row = dict(tag=tag, axis="null", method=f"deseq2__{design}",
                       **{f"n_sig_T{t}": n for t, n in enumerate(n_sig)},
                       **{f"jacc_T{i}T{j}": v for (i, j), v in zip(pairs, jaccs)},
                       jacc_mean=float(np.mean(jaccs)))
            rows.append(row)
            pd.DataFrame(rows).to_csv(out_path, index=False,
                                      mode="a" if out_path.exists() else "w",
                                      header=not out_path.exists())
            rows = []
            log(f"{tag}/deseq2__{design}: n_sig={','.join(map(str, n_sig))} "
                f"jacc_mean={row['jacc_mean']:.3f}  ({time.time()-t0:.0f}s)")
    log("done")


if __name__ == "__main__":
    main()
