"""Write lobo_folds.json using exactly the engine's own LOBO definition.

From validation/lobo_engine.py and 3_detection_limit.ipynb:
  train = all HC except the held-out batch (and except batches with < MIN_HC_BATCH_SIZE
          HC -- already excluded from hc_meta.csv, whose smallest batch is 5)
  test  = that batch's HC only (disease samples are never scored here: they carry real
          outliers that would be labelled negative)

Batches are the five 3_detection_limit.ipynb evaluates the engine on, i.e. the largest by
held-out HC count. Note these are NOT lobo_test_meta.json's six: that set is tier-A
(n_dis > 0) for the MMD analysis, which drops Moufarrej_Batch_2 (93 held-out HC, no
disease) and adds Moufarrej_Batch_1/_4. Injection needs HC, not disease, so the
detection-limit set is the right one to match.

One deliberate difference from the notebook: it scored 5 random half-subsamples of each
batch's held-out HC ("boot"), which only discards data and has no baseline counterpart.
Here every held-out HC sample is scored once.
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config

from common import CACHE

BATCHES = ["Ward Z et al._Batch_1", "Moufarrej et al._Batch_2", "Moore et al._Batch_1",
           "Chen et al._Batch_2", "Roskams-Hieter B et al._Batch_2"]


def safe_dir(batch):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", batch)


def main():
    hc = pd.read_csv(config.OUTRIDER_COMPARISON_DIR / "hc_meta.csv")
    folds = {}
    for b in BATCHES:
        bdir = config.LOBO_MIXED_DIR / safe_dir(b)
        if not (bdir / "model_fits.csv").exists():
            raise FileNotFoundError(f"no engine LOBO fit for {b} at {bdir}")
        test = hc.loc[hc["batch"] == b, "sample"].tolist()
        train = hc.loc[hc["batch"] != b, "sample"].tolist()
        folds[b] = dict(train=train, test=test, batch_dir=safe_dir(b))
        print(f"{b:38} train={len(train):4d}  test(HC)={len(test):4d}")

    path = CACHE / "lobo_folds.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(folds, indent=1))
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
