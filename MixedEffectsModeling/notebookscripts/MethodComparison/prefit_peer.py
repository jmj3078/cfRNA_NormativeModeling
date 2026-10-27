"""Pre-fit and cache the PEER basis for whichever folds already have an autoencoder fit.

PEER is the one arm whose runtime at this scale (12k genes x ~540 samples, K=60) was
unknown, and it runs in its own env, so fitting it ahead of the sweep both de-risks and
parallelizes. Cached fits are reused by sweep.py untouched.

Usage: python prefit_peer.py [cv|lobo|both]
"""
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

import peer_arm
from frozen_latent import FITS, load_fit
from sweep import SPLITS, fold_spec

import MixedEffectsModeling.config as config

# PEER is single-threaded, so folds run concurrently; the work happens in R subprocesses
# so threads are enough. Capped to leave headroom -- each fit held ~1 GB.
MAX_PARALLEL = 5


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    splits = SPLITS if which == "both" else (which,)
    counts = pd.read_csv(config.OUTRIDER_COMPARISON_DIR / "hc_counts.csv", index_col=0)

    jobs = []
    for split in splits:
        for _, tag, tr_names, _, _ in fold_spec(split):
            if not (FITS / f"autoencoder_{tag}_genes.csv").exists():
                print(f"skip {split}/{tag}: autoencoder fit not ready", flush=True)
                continue
            genes = load_fit("autoencoder", tag)["genes"]
            jobs.append((f"{split}_{tag}",
                         counts.loc[genes, tr_names].to_numpy(dtype=float).T))

    def run(job):
        tag, y_tr = job
        t0 = time.time()
        fit = peer_arm.fit_train(y_tr, tag)
        print(f"{tag}: W={fit['w'].shape} in {time.time() - t0:.0f}s", flush=True)

    with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as ex:
        list(ex.map(run, jobs))


if __name__ == "__main__":
    main()
