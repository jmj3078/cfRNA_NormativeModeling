#!/usr/bin/env python
"""Train NormativeModelEngineV2 on ALL protein-coding genes and report route-level
demotion/exclusion statistics.

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
parser.add_argument("--nz-bc-split", type=int, default=None)
args = parser.parse_args()

t0 = time.perf_counter()

engine = NormativeModelEngineV2(nz_a_max=args.nz_a_max, nz_bc_split=args.nz_bc_split)
engine.load_hc_data()
engine.build_dispersion_trend()
engine.assign_routes()
engine.train(verbose=True, limit=args.limit)
engine.save(SAVE_DIR)

print(f"\nTotal time: {(time.perf_counter() - t0) / 3600:.2f}h")

# ---- Summary statistics ---------------------------------------------------
df = engine.training_summary()
print("\n" + "=" * 70)
print("ROUTE-LEVEL SUMMARY")
print("=" * 70)

print("\nInitial gating (Phase 1):")
print(df["initial_route"].value_counts().to_string())

print("\nFinal route (after demotion chain):")
print(df["route"].value_counts().to_string())

print("\nDemotion cross-tab (initial_route x final route):")
print(pd.crosstab(df["initial_route"], df["route"]).to_string())

n_c_to_b = int(((df["initial_route"] == "C") & (df["route"] == "B")).sum())
n_b_to_a = int(((df["initial_route"] == "B") & (df["route"] == "A")).sum())
n_excluded = int((df["route"] == "excluded").sum())
n_total = len(df)

print(f"\nC -> B demotions : {n_c_to_b}")
print(f"B -> A demotions : {n_b_to_a}")
print(f"EXCLUDED (C->B->fail): {n_excluded}  ({n_excluded/n_total:.2%} of all genes)")

nbi_chosen = df.loc[df["route"] == "C", "nbi_chosen"].value_counts()
b_chosen = df.loc[df["route"] == "B", "b_chosen"].value_counts()
print(f"\nRoute C GAIC choice: {nbi_chosen.to_dict()}")
print(f"Route B GAIC choice: {b_chosen.to_dict()}")

# ---- Figure: gating vs final route stacked bar -----------------------------
apply_style()

fig, ax = plt.subplots(figsize=(7, 5))
ct = pd.crosstab(df["initial_route"], df["route"])
ct = ct.reindex(index=["A", "B", "C"], columns=[c for c in ["A", "B", "C", "excluded"] if c in ct.columns])
ct.plot(kind="bar", stacked=True, ax=ax, colormap="tab10")
ax.set(xlabel="initial route (Phase 1 gating)", ylabel="gene count",
      title="Route demotion outcome by initial gating")
ax.legend(title="final route")
fig.tight_layout()
fig_path = SAVE_DIR / "route_demotion_summary.png"
fig.savefig(fig_path, dpi=150)
print(f"\nFigure saved -> {fig_path}")
print(f"Engine saved  -> {SAVE_DIR}")
