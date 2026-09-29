import pickle
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.validation.lobo_mmd import plot_mmd_direction
from viz_style import apply_style

apply_style()

DIR = Path(__file__).parent
(DIR / "Figures").mkdir(exist_ok=True)


if __name__ == "__main__":
    outr_mmd = pd.read_csv(DIR / "outrider_mmd_summary.csv")
    with open(DIR / "outrider_mmd_raw.pkl", "rb") as f:
        raw_outrider = pickle.load(f)
    plot_mmd_direction(outr_mmd, raw_outrider, DIR / "Figures" / "mmd_direction_outrider_held_out.png")

    eng_mmd_shash = pd.read_csv(config.LOBO_MIXED_DIR / "mmd_summary_shash.csv")
    with open(config.LOBO_MIXED_DIR / "mmd_raw_shash.pkl", "rb") as f:
        raw_ours_shash = pickle.load(f)
    plot_mmd_direction(eng_mmd_shash, raw_ours_shash, DIR / "Figures" / "mmd_direction_our_engine_shash.png")

    merged = eng_mmd_shash[["batch", "n_hc", "n_dis", "perm_p", "p_direction", "disease_farther"]].merge(
        outr_mmd[["batch", "perm_p", "p_direction", "disease_farther"]],
        on="batch", suffixes=("_ours", "_outrider_held_out"))
    merged = merged.sort_values("n_dis", ascending=False)
    merged.to_csv(DIR / "mmd_comparison_held_out.csv", index=False)
    print(f"direction-significant: ours={int((merged.p_direction_ours < 0.05).sum())}/6  "
          f"outrider_held_out={int((merged.p_direction_outrider_held_out < 0.05).sum())}/6")
    print(merged.to_string())
