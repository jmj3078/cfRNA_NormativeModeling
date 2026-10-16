"""Freeze the CV and LOBO splits once. Every pipeline reads these files; nothing rewrites them.

cv_folds.json   = StratifiedKFold(5, shuffle, rs=42) on the HC pool by Batch_ID, as cv_engine.
lobo_folds.json = held-out batches with >= LOBO_MIN_TEST_HC HC and >= 1 disease sample;
                  train = every other HC of the pool, test = that batch's HC (MethodComparison
                  build_lobo_folds.py). Disease samples of the batch are scored by lobo_engine.
Both are checked against the production split files before being written.
"""
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _FinalModelingMethodImprovement.config as config
from _FinalModelingMethodImprovement.engine.core.model_engine_mixed import NormativeModelEngineMixed
from _FinalModelingMethodImprovement.engine.validation.lobo_engine import load_full_data, usable_batches

PROD_CV = config.BASELINE_MEM_DIR / "OutriderComparison" / "insample_comparison" / "cv_folds.json"
PROD_LOBO = config.BASELINE_MEM_DIR / "MethodComparison" / "cache" / "lobo_folds.json"


def cv_folds():
    e = NormativeModelEngineMixed()
    e.load_hc_data()
    folds = StratifiedKFold(config.SPIKE_PARAMS["n_splits"], shuffle=True, random_state=42).split(
        np.zeros(len(e.names)), e.batch)
    return {str(i): dict(train=e.names[tr].tolist(), test=e.names[te].tolist()) for i, (tr, te) in enumerate(folds)}


def lobo_folds():
    data = load_full_data()
    pool = data["is_hc"] & ~np.isin(data["batch"], list(data["small_hc_batches"]))
    names, batch = data["names"][pool], data["batch"][pool]
    ub = usable_batches(data)
    out = {}
    for b in ub.loc[ub["n_dis"] > 0, "batch"] if config.LOBO_REQUIRE_DISEASE else ub["batch"]:
        test = batch == b
        if test.sum() >= config.LOBO_MIN_TEST_HC:
            out[b] = dict(train=names[~test].tolist(), test=names[test].tolist())
    return out


def check(folds):
    for k, f in folds.items():
        assert not set(f["train"]) & set(f["test"]), f"train/test overlap in {k}"


def main():
    config.SPLITS_DIR.mkdir(parents=True, exist_ok=True)
    cv, lobo = cv_folds(), lobo_folds()
    check(cv), check(lobo)
    assert sorted(sum((f["test"] for f in cv.values()), [])) == sorted(cv["0"]["train"] + cv["0"]["test"])
    assert cv == json.loads(PROD_CV.read_text()), "CV folds differ from production"
    prod = {b: dict(train=f["train"], test=f["test"]) for b, f in json.loads(PROD_LOBO.read_text()).items()}
    assert lobo == prod, "LOBO folds differ from production"
    for path, folds in ((config.CV_FOLDS_PATH, cv), (config.LOBO_FOLDS_PATH, lobo)):
        if path.exists():
            assert json.loads(path.read_text()) == folds, f"{path} exists with different content"
        else:
            path.write_text(json.dumps(folds, indent=1))
    print("cv:", {k: len(f["test"]) for k, f in cv.items()})
    print("lobo:", {k: len(f["test"]) for k, f in lobo.items()})


if __name__ == "__main__":
    main()
