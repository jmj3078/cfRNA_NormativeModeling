#!/usr/bin/env python
"""Train NormativeModelEngineV2 on ALL protein-coding genes and report
demotion-chain statistics: nbi -> nb_fixed (GAIC full-vs-intercept) -> intercept
-> excluded. Route pool (rare pooling) only ever comes from nz < nz_a_max.

  python run_model_engine_v2.py                # full run
  python run_model_engine_v2.py --limit 200     # smoke test
"""

import argparse
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config
from model_engine_v2 import NormativeModelEngineV2
from viz_style import apply_style

SAVE_DIR = config.ENGINE_DIR_V2

parser = argparse.ArgumentParser()
parser.add_argument("--limit", type=int, default=None, help="limit to first N genes (smoke test)")
parser.add_argument("--nz-a-max", type=int, default=None)
args = parser.parse_args()

t0 = time.perf_counter()

engine = NormativeModelEngineV2(nz_a_max=args.nz_a_max)
engine.load_hc_data()
engine.build_dispersion_trend()
engine.assign_routes()
engine.train(verbose=True, limit=args.limit)
engine.save(SAVE_DIR)

print(f"\nTotal time: {(time.perf_counter() - t0) / 3600:.2f}h")

# ---- Summary statistics ---------------------------------------------------
df = engine.training_summary()
print("\n" + "=" * 70)
print("DEMOTION-CHAIN SUMMARY")
print("=" * 70)

print("\nInitial gating (Phase 1):")
print(df["initial_route"].value_counts().to_string())

print("\nFinal route (after demotion chain):")
print(df["route"].value_counts().to_string())

print("\nDemotion cross-tab (initial_route x final route):")
print(pd.crosstab(df["initial_route"], df["route"]).to_string())

n_model_candidates = int((df["initial_route"] == "model").sum())
n_nbi_fitted = int((df["stage"] == "nbi").sum())
n_nb_fixed = int((df["stage"] == "nb_fixed").sum())
n_intercept = int((df["stage"] == "intercept").sum())
n_excluded = int((df["route"] == "excluded").sum())
n_total = len(df)

print(f"\nmodel-candidates                          : {n_model_candidates}")
print(f"  fitted at stage nbi (full NBI)           : {n_nbi_fitted}")
print(f"  demoted to stage nb_fixed (GAIC full/intercept): {n_nb_fixed}")
print(f"  demoted further to stage intercept        : {n_intercept}")
print(f"  EXCLUDED (stage intercept failed)         : {n_excluded}  "
     f"({n_excluded/n_total:.2%} of all genes)")

mean_model_by_stage = df.loc[df["route"] == "model"].groupby(["stage", "mean_model_chosen"]).size()
print(f"\nModel outcome by stage (mean_model_chosen only applies to nb_fixed/intercept):")
print(mean_model_by_stage.to_string())

# ---- Figure: demotion funnel -----------------------------------------------
apply_style()

funnel_labels = ["pool (nz<%d)" % engine.nz_a_max, "model: nbi", "model: nb_fixed (full)",
                "model: nb_fixed (intercept)", "model: intercept", "excluded"]
funnel_counts = [
    int((df["route"] == "pool").sum()),
    n_nbi_fitted,
    int(((df["stage"] == "nb_fixed") & (df["mean_model_chosen"] == "full")).sum()),
    int(((df["stage"] == "nb_fixed") & (df["mean_model_chosen"] == "intercept")).sum()),
    n_intercept,
    n_excluded,
]

fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(funnel_labels, funnel_counts, color=["#1b9e77", "#7570b3", "#d95f02", "#e7298a", "#66a61e", "#666666"])
ax.set(ylabel="gene count", title="Final route outcome (demotion chain: C -> B -> intercept fallback)")
ax.tick_params(axis="x", rotation=30)
for i, v in enumerate(funnel_counts):
    ax.text(i, v, str(v), ha="center", va="bottom")
fig.tight_layout()
fig_path = SAVE_DIR / "route_demotion_summary.png"
fig.savefig(fig_path, dpi=150)
print(f"\nFigure saved -> {fig_path}")
print(f"Engine saved  -> {SAVE_DIR}")
