import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

parent_dir = str(Path(__file__).resolve().parents[3])
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)
from viz_style import apply_style

from MixedEffectsModeling.validation.ppc_simulate import binned_medians, calib_stats

apply_style()

DIR = Path(__file__).parent


if __name__ == "__main__":
    eng_sub_cal = pd.read_csv(DIR / "our_engine_cv_calibration_moments_12305subset.csv")
    outr_cal = pd.read_csv(DIR / "outrider_cv_calibration_moments.csv")

    panels = [("obs_mean", "pred_mean", "Mean", True), ("obs_var", "pred_var", "Variance", True),
             ("obs_zero", "pred_zero", "Zero-Fraction", False)]
    engines = [("Our engine (12138g subset)", eng_sub_cal, "#2b5c8f"), ("OUTRIDER (held-out)", outr_cal, "#c0392b")]

    fig, axes = plt.subplots(3, 2, figsize=(10, 14))

    for row, (oc, pc, title_str, logsc) in enumerate(panels):
        for col, (eng_label, cal, color) in enumerate(engines):
            ax = axes[row, col]
            ax.scatter(cal[oc], cal[pc], s=6, alpha=0.15, color=color)
            bx, by = binned_medians(cal[oc], cal[pc])
            ax.plot(bx, by, '-', color='black', linewidth=2)

            if logsc:
                lo = max(min(cal[oc].min(), cal[pc].min()), 1e-3)
                hi = max(cal[oc].max(), cal[pc].max())
                ax.plot([lo, hi], [lo, hi], 'k--', linewidth=1, alpha=0.7)
                ax.set_xscale('log'); ax.set_yscale('log')
            else:
                ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.7)

            cs = calib_stats(cal[oc], cal[pc], logsc)
            ax.set_title(f"{title_str} -- {eng_label}", fontweight='bold', pad=8)
            ax.set_xlabel(f'Observed {title_str.lower()}')
            ax.set_ylabel(f'PPC replicate {title_str.lower()}')
            stats_text = ("r = {pearson_r:.3f}\nR2 = {r2:.3f}\nRMSE = {rmse:.3f}\nn = {n}".format(**cs))
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.85, edgecolor='#cccccc'))
            ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
            ax.grid(True, linestyle=':', alpha=0.4)
            if row == 1:
                for spine in ax.spines.values():
                    spine.set_visible(True)
                    spine.set_edgecolor('#c0392b' if col == 1 else '#2b5c8f')
                    spine.set_linewidth(2)

    plt.tight_layout()
    fig.savefig(DIR / 'Figures' / 'ppc_scatter_comparison_held_out.png', dpi=200, bbox_inches='tight')
    print(f"saved -> {DIR / 'Figures' / 'ppc_scatter_comparison_held_out.png'}")
