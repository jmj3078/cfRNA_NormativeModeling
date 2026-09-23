"""Null-split GSEA reproducibility arm -- missing counterpart to gsea_split_jaccard's tertile
run (that script/log is uncommitted; only its output summary.csv survived). Answers a direct
question raised during review: does the pathway-level Jaccard bar used to call the tertile splits
"reproducible" ALSO get cleared by a genuinely random (null) HC split, which carries no real
control-composition signal? If so, the tertile "reproducibility" claim is not distinguishing real
effect from a generically loose threshold.

Reuses cached DESeq2 stat files (T0/T1/T2_{design}.csv.gz) from run_control_composition_deseq2.py's
null draws (Pancreatic_Cancer__null_0000..0029) -- no refitting, GSEA only. Same prerank machinery
(gsea_prerank, housekeeping-excluded KEGG+Reactome library) as MixedEffectsModeling's own
group_level_pathway_gsea, for direct comparability.

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
from MixedEffectsModeling.Cohort_Grouped_Z.pathway_convergence import (
    ensg_to_symbol, gsea_prerank, load_pathway_library, load_symbol_vocab)

DESIGNS = ["no_covariate", "ruvg_k2"]
DISEASE = "Pancreatic Cancer"
OUT_CSV = config.CTRL_COMP_DIR / "gsea_split_jaccard" / "summary_null.csv"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def jaccard(a, b):
    a, b = set(a), set(b)
    u = a | b
    return len(a & b) / len(u) if u else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-null", type=int, default=10)
    args = ap.parse_args()

    universe_syms, sym2idx, col2sym = load_symbol_vocab(None)
    terms, M = load_pathway_library()
    log(f"library: {len(terms)} terms (housekeeping-excluded, matches pathway_convergence)")

    sym_of = ensg_to_symbol()
    sym_of.index = sym_of.index.str.split(".").str[0]

    rows = []
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if OUT_CSV.exists():
        prev = pd.read_csv(OUT_CSV)
        done = set(prev["tag"] + "/" + prev["method"])

    for draw in range(args.n_null):
        tag = f"{DISEASE.replace(' ', '_')}__null_{draw:04d}"
        stat_dir = config.CTRL_COMP_DESEQ2_DIR / tag
        if not stat_dir.exists():
            log(f"skip {tag}: no cached DESeq2 stats")
            continue
        for design in DESIGNS:
            key = f"{tag}/deseq2__{design}"
            if key in done:
                continue
            t0 = time.time()
            n_sig, sig_terms = [], []
            for t in range(3):
                stat = pd.read_csv(stat_dir / f"T{t}_{design}.csv.gz", index_col=0)["stat"]
                stat.index = stat.index.str.split(".").str[0]
                syms = sym_of.reindex(stat.index)
                rnk = pd.Series(stat.values, index=syms.values)
                rnk = rnk[pd.notna(rnk.index)].groupby(level=0).mean().reindex(universe_syms).values
                res2d = gsea_prerank(rnk, terms, M, universe_syms, n_perm=1000, seed=42 + t)
                sig = set(res2d.loc[res2d["FDR q-val"] < 0.05, "Term"])
                n_sig.append(len(sig))
                sig_terms.append(sig)
            jaccs = [jaccard(sig_terms[i], sig_terms[j]) for i, j in combinations(range(3), 2)]
            row = dict(tag=tag, axis="null", method=f"deseq2__{design}",
                      n_sig_T0=n_sig[0], n_sig_T1=n_sig[1], n_sig_T2=n_sig[2],
                      jacc_T0T1=jaccs[0], jacc_T0T2=jaccs[1], jacc_T1T2=jaccs[2],
                      jacc_mean=float(np.mean(jaccs)))
            rows.append(row)
            pd.DataFrame(rows).to_csv(OUT_CSV, index=False, mode="a" if OUT_CSV.exists() else "w",
                                      header=not OUT_CSV.exists())
            rows = []
            log(f"{tag}/deseq2__{design}: n_sig={n_sig[0]},{n_sig[1]},{n_sig[2]} "
                f"jacc_mean={row['jacc_mean']:.3f}  ({time.time()-t0:.0f}s)")
    log("done")


if __name__ == "__main__":
    main()
