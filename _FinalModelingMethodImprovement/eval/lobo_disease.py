"""Unseen-batch direction check: in each held-out LOBO batch, are the disease samples farther from the
normative model than the same batch's HC? Two axes per pipeline (eval/lobo_mmd.py): BH-significant genes
per sample (primary) and the RKHS distance to the held-out-HC mean embedding (secondary), one-sided
Mann-Whitney disease > HC, BH across batches. OOD-kept samples, Moore_1 disease = PDAC only,
batches with >= 25 HC and >= 10 disease.

Output: config.RESULTS_DIR / lobo_disease.csv. Run: python eval/lobo_disease.py p0 p2 p3
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _FinalModelingMethodImprovement.config as config
from _FinalModelingMethodImprovement.eval.lobo_mmd import load_engine_arm, nsig_stats, nsig_table, rkhs_stats
from _FinalModelingMethodImprovement.eval.results import append

OUT = config.RESULTS_DIR / "lobo_disease.csv"


def main(pipes):
    for p in pipes:
        arm = load_engine_arm(lobo_dir=config.FMMI / "runs" / p / "lobo")
        ns = nsig_stats(nsig_table(arm)).assign(axis="n_sig")
        rk = rkhs_stats(arm)[0].assign(axis="rkhs_distance")
        append(OUT, pd.concat([ns, rk], ignore_index=True).assign(method=p, scenario="zero_shot"))
        print(f"{p} done", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
