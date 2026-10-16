"""Gate S1-a: with the production fits copied in, the forked cv_engine / lobo_engine must
rescore CV and LOBO Z bit for bit. Builds runs/p0 from the pre-Ward production state; refuses to overwrite an existing run."""
import os
import pickle
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

os.environ["FMMI_PIPELINE"] = "p0"
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import _FinalModelingMethodImprovement.config as config

PROD = config.BASELINE_MEM_DIR
LOBO_CACHE = ["model_fits.csv", "pool_fit.json", "shash_params.csv"]


def same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and bool(((a == b) | (np.isnan(a) & np.isnan(b))).all())


def stage():
    assert not config.RUN_DIR.exists(), f"{config.RUN_DIR} exists; move it away first"
    shutil.copytree(PROD / "engine_state_mixed", config.ENGINE_MIXED_DIR)
    shutil.copytree(PROD / "CV_Results_mixed" / "fold_params", config.CV_MIXED_DIR / "fold_params")
    for d in (PROD / "LOBO_Results_mixed").glob("*/meta.json"):
        (config.LOBO_MIXED_DIR / d.parent.name).mkdir(parents=True)
        for f in LOBO_CACHE:
            shutil.copy2(d.parent / f, config.LOBO_MIXED_DIR / d.parent.name / f)


def run(script):
    subprocess.run([sys.executable, str(config.FMMI / "engine" / "validation" / script)],
                   check=True, env=dict(os.environ, FMMI_PIPELINE="p0"))


def compare():
    ok = True
    for f in ("cv_zscores.pkl", "cv_zscores_shash.pkl"):
        a, b = (pickle.load(open(d / f, "rb")) for d in (config.CV_MIXED_DIR, PROD / "CV_Results_mixed"))
        r = a.keys() == b.keys() and all(same(a[g], b[g]) for g in b)
        print(f"cv {f}: {r}")
        ok &= r
    a, b = (pd.read_csv(d / "cv_stats.csv") for d in (config.CV_MIXED_DIR, PROD / "CV_Results_mixed"))
    r = a.equals(b)
    print(f"cv cv_stats.csv: {r}")
    ok &= r
    for d in sorted((PROD / "LOBO_Results_mixed").glob("*/meta.json")):
        for f in ("Z_test.npy", "Z_test_shash.npy", "ood_mask.npy"):
            r = same(np.load(config.LOBO_MIXED_DIR / d.parent.name / f), np.load(d.parent / f))
            print(f"lobo {d.parent.name} {f}: {r}")
            ok &= r
    assert ok


if __name__ == "__main__":
    if "--compare-only" not in sys.argv:
        stage()
        run("cv_engine.py")
        run("lobo_engine.py")
    compare()
