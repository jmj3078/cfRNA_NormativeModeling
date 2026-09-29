"""Ideal-reproducibility ceiling for the rank-stability check (9_control_split_variation.ipynb
sec. 11): two INDEPENDENT bootstrap draws of controls (same size as a tertile subset, drawn
WITH replacement so they can overlap -- no forced non-overlap, no compositional bias), case
group fixed. This isolates pure sampling noise from the forced-disjoint-partition and
bias-axis-composition effects the rest of the notebook measures, giving an upper bound on
what rho this sample size could realistically produce under the best case.

Run: python EDA/control_composition/run_ideal_ceiling.py [--n-pairs 4]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from run_control_composition import DISEASES, build_cache
from run_control_composition_deseq2 import fit_one, load_raw_counts

OUT_CSV = config.CTRL_COMP_DIR / "ideal_ceiling_bootstrap.csv"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-pairs", type=int, default=4)
    ap.add_argument("--seed", type=int, default=123)
    args = ap.parse_args()

    data = build_cache()
    counts = load_raw_counts(data)
    obs = data["obs"]
    samples = obs["sample"].values
    rng = np.random.default_rng(args.seed)

    rows = []
    for disease in DISEASES:
        case_idx = np.where(obs["phenotype"].values == disease)[0]
        hc_idx = np.where(obs["phenotype"].values == "Healthy Control")[0]
        n_ctrl = len(hc_idx) // 3  # match tertile subset size
        for p in range(args.n_pairs):
            t0 = time.time()
            stats = []
            for rep in range(2):
                boot_ctrl = rng.choice(hc_idx, size=n_ctrl, replace=True)
                idx = np.concatenate([case_idx, boot_ctrl])
                sub_samples = samples[idx]
                cnt = counts.loc[sub_samples].copy()
                cnt.index = [f"{s}__{i}" for i, s in enumerate(sub_samples)]
                cond = pd.DataFrame({"condition": ["disease"] * len(case_idx) + ["HC"] * n_ctrl}, index=cnt.index)
                stats.append(fit_one(cnt, cond, "~condition"))
            rho = spearmanr(stats[0], stats[1]).correlation
            rows.append(dict(disease=disease, pair=p, rho=rho))
            log(f"{disease} pair {p}: rho={rho:.3f} ({time.time()-t0:.0f}s)")

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    log("done")
    print(out.groupby("disease")["rho"].median())


if __name__ == "__main__":
    main()
