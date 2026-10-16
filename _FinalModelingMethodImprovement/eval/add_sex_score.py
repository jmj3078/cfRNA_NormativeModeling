"""Add obs['Sex Score'] to the merged h5ad in place: mean log2(CPM + 1) of the chrY genes of the Staedtler et al. 2013
sex-marker panel (Biomarkers 18:436, doi:10.3109/1354750X.2013.811538; XIST, the fifth gene, is not used so the score
needs no X-chromosome normalisation). CPM over all genes of X (raw counts). Defined for every sample, so it covers the
HC whose sex metadata is missing. Rerun after any h5ad rebuild (e.g. the Ward lane-merge refit); an existing column is replaced.

Run: python eval/add_sex_score.py
"""
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _FinalModelingMethodImprovement.config as config
from _FinalModelingMethodImprovement.engine.validation.lobo_engine import _decode_cat

SEX_GENES = ["RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D"]
COLUMN = "Sex Score"


def sex_score(f):
    sym = _decode_cat(f["var/GeneName"])
    chrom = _decode_cat(f["var/Chromosome"])
    idx = np.where(np.isin(sym, SEX_GENES) & np.isin(chrom, ["chrY", "Y"]))[0]
    assert sorted(sym[idx]) == sorted(SEX_GENES), sym[idx]
    X = csr_matrix((f["X/data"][()], f["X/indices"][()], f["X/indptr"][()]), shape=tuple(f["X"].attrs["shape"]))
    lib = np.asarray(X.sum(1)).ravel()
    return np.log2(X[:, idx].toarray() / lib[:, None] * 1e6 + 1).mean(1)


def main():
    with h5py.File(config.H5AD_PATH, "r+") as f:
        score = sex_score(f)
        obs = f["obs"]
        if COLUMN in obs:
            del obs[COLUMN]
        ds = obs.create_dataset(COLUMN, data=score.astype(np.float64))
        ds.attrs["encoding-type"], ds.attrs["encoding-version"] = "array", "0.2.0"
        order = [c for c in obs.attrs["column-order"] if c != COLUMN] + [COLUMN]
        obs.attrs["column-order"] = np.array(order, dtype=h5py.string_dtype())
        sex = _decode_cat(obs["sex"])
    for s in ("male", "female"):
        print(f"{s}: n={int((sex == s).sum())}, median {COLUMN} {np.median(score[sex == s]):.3f}")
    print(f"wrote obs['{COLUMN}'] for {len(score)} samples to {config.H5AD_PATH}")


if __name__ == "__main__":
    main()
