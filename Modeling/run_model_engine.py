#!/usr/bin/env python
"""Run NormativeModelEngine: assign branches → train → save."""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config
from model_engine import NormativeModelEngine

SAVE_DIR = config.ENGINE_DIR
_MP = config.MODELING_PARAMS

engine = NormativeModelEngine(
    count_model         = "nbi",
    low_det_thr         = _MP["low_det_thr"],
    det_rate_min        = _MP["det_rate_min"],
    mean_count_min      = _MP["mean_count_min"],
    lr_C                = _MP["lr_c"],
    nbi_outlier_z       = 5.0,
    nbi_max_iter        = 2,
    nbi_max_remove_frac = 0.05,
    lambda_sigma        = 0.05,
)

t0 = time.perf_counter()

engine.load_hc_data()
engine.assign_branches()
engine.train(verbose=True)
engine.save(SAVE_DIR)

print(f"\nTotal time: {(time.perf_counter() - t0) / 3600:.2f}h")
print(f"Engine saved → {SAVE_DIR}")
