"""Write the count/metadata inputs every non-engine arm reads.

These files had no producer in the repo -- they were built once by hand, which is how
the 2026-09-30 Batch_ID redefinition could silently invalidate them. Regenerating them
is now part of the refit chain:

  run_engine.py -> export_inputs.py -> run_outrider_cv.R / MethodComparison
  lobo_engine.py -> export_inputs.py --lobo -> run_outrider_lobo.R

Gene set = the engine's own genes minus route "excluded", in training_summary order, so
every arm scores the same genes. Sample set = NormativeModelEngineMixed.load_hc_data's
pool for HC, every other QC-passing sample for disease.

cv_folds.json is NOT written here -- cv_engine.py writes it from the split it actually
used.
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.validation.lobo_engine import load_full_data, selected_batches

OUT = config.OUTRIDER_COMPARISON_DIR


def export_lobo_test_meta():
    metas = selected_batches()
    if not metas:
        raise SystemExit(f"no qualifying LOBO fit in {config.LOBO_MIXED_DIR}")
    payload = {m["batch_id"]: dict(test_names=m["test_names"], test_is_hc=m["test_is_hc"])
               for m in metas}
    (OUT / "lobo_test_meta.json").write_text(json.dumps(payload, indent=1))
    for b, v in payload.items():
        print(f"{b:38} n_test={len(v['test_names']):4d}  n_test_hc={sum(v['test_is_hc']):4d}")
    print(f"wrote {OUT / 'lobo_test_meta.json'} ({len(payload)} batches)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lobo", action="store_true",
                    help="only rewrite lobo_test_meta.json (needs lobo_engine.py to have run)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    if args.lobo:
        export_lobo_test_meta()
        return

    summary = pd.read_csv(config.ENGINE_MIXED_DIR / "training_summary.csv")
    genes = summary.loc[summary["route"] != "excluded", "gene"].tolist()

    data = load_full_data()
    is_hc = data["is_hc"] & ~np.isin(data["batch"], list(data["small_hc_batches"]))
    cols = [data["gene_col"][g] for g in genes]
    Y = data["Y"][:, cols].astype(np.int64)

    for tag, mask in (("hc", is_hc), ("disease", ~data["is_hc"])):
        df = pd.DataFrame(Y[mask].T, index=genes, columns=data["names"][mask])
        df.to_csv(OUT / f"{tag}_counts.csv")
        print(f"{tag}_counts.csv: {df.shape[0]} genes x {df.shape[1]} samples")

    pd.DataFrame(dict(sample=data["names"][is_hc], batch=data["batch"][is_hc])).to_csv(
        OUT / "hc_meta.csv", index=False)
    print(f"hc_meta.csv: {int(is_hc.sum())} HC in {pd.Series(data['batch'][is_hc]).nunique()} batches")

    with h5py.File(config.H5AD_PATH, "r") as f:
        length = pd.Series(f["var/Length"][()],
                           index=[g.decode() if isinstance(g, bytes) else g
                                  for g in f["var/_index"][()]])
    length.loc[genes].rename("Length").to_csv(OUT / "gene_lengths.csv")
    print(f"gene_lengths.csv: {len(genes)} genes")


if __name__ == "__main__":
    main()
