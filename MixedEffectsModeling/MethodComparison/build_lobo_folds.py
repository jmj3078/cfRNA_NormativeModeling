"""Write lobo_folds.json using exactly the engine's own LOBO definition.

From validation/lobo_engine.py and 3_detection_limit.ipynb:
  train = all HC except the held-out batch (and except batches with < MIN_HC_BATCH_SIZE
          HC -- already excluded from hc_meta.csv)
  test  = that batch's HC only (disease samples are never scored here: they carry real
          outliers that would be labelled negative)

Held-out batches are derived from the fitted HC pool, not hardcoded: a batch qualifies
when it holds >= config.LOBO_MIN_TEST_HC HC, carries disease samples if
config.LOBO_REQUIRE_DISEASE, and has an engine LOBO fit on disk. Deriving them keeps the
set correct when the Batch_ID definition changes -- the previous hardcoded list survived
the 2026-09-30 batch redefinition naming a batch that no longer exists.

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


def safe_dir(batch):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", batch)


def select_batches(hc, disease):
    n_hc = hc["batch"].value_counts()
    n_dis = disease["batch"].value_counts()
    cand = n_hc[n_hc >= config.LOBO_MIN_TEST_HC]
    if config.LOBO_REQUIRE_DISEASE:
        cand = cand[[n_dis.get(b, 0) > 0 for b in cand.index]]
    missing = [b for b in cand.index
               if not (config.LOBO_MIXED_DIR / safe_dir(b) / "model_fits.csv").exists()]
    if missing:
        raise FileNotFoundError(
            f"no engine LOBO fit for {len(missing)} qualifying batch(es): {missing}\n"
            f"run the LOBO engine for these before building folds")
    return cand.index.tolist(), n_dis


def main():
    hc = pd.read_csv(config.OUTRIDER_COMPARISON_DIR / "hc_meta.csv")
    disease = pd.read_csv(config.ZSCORES_MIXED_DIR / "sample_meta.csv")
    batches, n_dis = select_batches(hc, disease)
    print(f"held-out batches: n_hc >= {config.LOBO_MIN_TEST_HC}"
          f"{' and n_dis > 0' if config.LOBO_REQUIRE_DISEASE else ''} "
          f"-> {len(batches)} of {hc['batch'].nunique()} in the HC pool")

    folds = {}
    for b in batches:
        test = hc.loc[hc["batch"] == b, "sample"].tolist()
        train = hc.loc[hc["batch"] != b, "sample"].tolist()
        folds[b] = dict(train=train, test=test, batch_dir=safe_dir(b))
        print(f"{b:38} train={len(train):4d}  test(HC)={len(test):4d}  n_dis={n_dis.get(b, 0):4d}")

    path = CACHE / "lobo_folds.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(folds, indent=1))
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
