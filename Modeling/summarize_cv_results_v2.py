#!/usr/bin/env python
"""Summarize cv_model_engine_v2.py output: per-route calibration stats + figures.

Reads CV_Results_v2/cv_stats.csv (written by cv_model_engine_v2.py) and produces
a route-level summary table plus diagnostic figures, mirroring the style of the
existing cv_gamlss_nb.py Figures/ outputs.

Usage:
    python summarize_cv_results_v2.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config
from viz_style import apply_style

apply_style()

CV_DIR = config.CV_RESULTS_DIR_V2
FIG_DIR = config.CV_FIG_DIR_V2


def main():
    stats_path = CV_DIR / "cv_stats.csv"
    if not stats_path.exists():
        raise FileNotFoundError(f"{stats_path} not found -- run cv_model_engine_v2.py first")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(stats_path)
    print(f"Loaded {len(df)} genes from {stats_path}")

    # ---- Route-level summary table ----
    summary = df.groupby("route").agg(
        n_genes=("gene", "size"),
        median_w1=("w1", "median"), mean_w1=("w1", "mean"),
        median_std_z=("std_z", "median"),
        median_skew_z=("skew_z", "median"),
        median_kurt_z=("kurt_z", "median"),
        frac_w1_lt_025=("w1", lambda s: float((s < 0.25).mean())),
    ).reset_index()
    summary.to_csv(CV_DIR / "cv_summary_by_route.csv", index=False)
    print("\nRoute-level calibration summary:")
    print(summary.round(3).to_string(index=False))

    # ---- NZ-binned summary (within each route, is calibration stable across NZ?) ----
    edges = sorted(set([0, config.MODELING_PARAMS_V2["nz_a_max"], 40, 70, 100, 200, 400, 700]))
    df["nz_bin"] = pd.cut(df["nz"], bins=edges, right=False)
    nz_summary = df.groupby(["route", "nz_bin"], observed=True).agg(
        n=("gene", "size"), median_w1=("w1", "median"), median_std_z=("std_z", "median"),
    ).reset_index()
    nz_summary.to_csv(CV_DIR / "cv_summary_by_route_nz.csv", index=False)

    # ---- Figures ----
    routes = sorted(df["route"].dropna().unique())
    colors = {"A": "#1b9e77", "B": "#d95f02", "C": "#7570b3"}

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for r in routes:
        sub = df[df["route"] == r]
        c = colors.get(r, "gray")
        axes[0].hist(sub["w1"].clip(0, 1), bins=40, alpha=0.5, color=c, label=f"Route {r} (n={len(sub)})")
        axes[1].hist(sub["std_z"].clip(0, 3), bins=40, alpha=0.5, color=c, label=f"Route {r}")
        axes[2].hist(sub["kurt_z"].clip(-2, 15), bins=40, alpha=0.5, color=c, label=f"Route {r}")
    axes[0].set(xlabel="W1 to N(0,1)", ylabel="genes", title="Calibration (W1) by route")
    axes[0].legend()
    axes[1].axvline(1, color="k", ls="--", lw=1)
    axes[1].set(xlabel="std(z)", title="Spread of held-out z")
    axes[1].legend()
    axes[2].axvline(0, color="k", ls="--", lw=1)
    axes[2].set(xlabel="excess kurtosis", title="Tail behavior")
    axes[2].legend()
    fig.tight_layout(); fig.savefig(FIG_DIR / "cv_calibration_by_route.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for r in routes:
        sub = df[df["route"] == r].sort_values("nz")
        ax.scatter(sub["nz"], sub["w1"].clip(0, 1), s=4, alpha=0.3, color=colors.get(r, "gray"), label=f"Route {r}")
    ax.axvline(config.MODELING_PARAMS_V2["nz_a_max"], color="k", ls="--", lw=1)
    ax.set(xlabel="NZ", ylabel="W1 to N(0,1)", xscale="log", title="Calibration vs NZ across routes")
    ax.legend()
    fig.tight_layout(); fig.savefig(FIG_DIR / "cv_w1_vs_nz.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved -> {CV_DIR}/cv_summary_by_route.csv, cv_summary_by_route_nz.csv")
    print(f"Figures -> {FIG_DIR}/cv_calibration_by_route.png, cv_w1_vs_nz.png")


if __name__ == "__main__":
    main()
