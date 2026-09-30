"""Write lobo_folds.json using exactly the engine's own LOBO definition.

From validation/lobo_engine.py and 3_detection_limit.ipynb:
  train = all HC except the held-out batch (and except batches with < MIN_HC_BATCH_SIZE
          HC -- already excluded from hc_meta.csv)
  test  = that batch's HC only (disease samples are never scored here: they carry real
          outliers that would be labelled negative)

Held-out batches come from lobo_engine.selected_batches(), i.e. from the fits on disk
plus config.LOBO_MIN_TEST_HC / LOBO_REQUIRE_DISEASE -- never a hardcoded list. The
previous hardcoded list survived the 2026-09-30 batch redefinition naming a batch that
no longer exists.

One deliberate difference from the notebook: it scored 5 random half-subsamples of each
batch's held-out HC ("boot"), which only discards data and has no baseline counterpart.
Here every held-out HC sample is scored once.
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.validation.lobo_engine import selected_batches

from common import CACHE


def main():
    hc = pd.read_csv(config.OUTRIDER_COMPARISON_DIR / "hc_meta.csv")
    metas = selected_batches()
    if not metas:
        raise SystemExit(f"no LOBO fit in {config.LOBO_MIXED_DIR} qualifies "
                         f"(n_test_hc >= {config.LOBO_MIN_TEST_HC}"
                         f"{', n_dis > 0' if config.LOBO_REQUIRE_DISEASE else ''})")
    print(f"held-out batches: {len(metas)} of {hc['batch'].nunique()} in the HC pool")

    folds = {}
    for m in metas:
        b = m["batch_id"]
        test = hc.loc[hc["batch"] == b, "sample"].tolist()
        train = hc.loc[hc["batch"] != b, "sample"].tolist()
        if not test:
            raise SystemExit(f"{b} has an engine LOBO fit but no HC in hc_meta.csv -- "
                             f"the fit predates the current HC pool, refit it")
        folds[b] = dict(train=train, test=test, batch_dir=m["batch_dir"])
        print(f"{b:38} train={len(train):4d}  test(HC)={len(test):4d}  n_dis={m['n_dis']:4d}")

    path = CACHE / "lobo_folds.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(folds, indent=1))
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
