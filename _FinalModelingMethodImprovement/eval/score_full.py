"""In-sample scores of the full engine (4_disease_scoring.ipynb cells 1, 4-6):
HC pool and disease Z, raw and SHASH-corrected, written to runs/<pipe>/full/scores/."""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _FinalModelingMethodImprovement.config as config
from _FinalModelingMethodImprovement.engine.core.model_engine_mixed import NormativeModelEngineMixed
from _FinalModelingMethodImprovement.engine.core.shash import load_shash_params, shash_correct_matrix
from _FinalModelingMethodImprovement.engine.validation.lobo_engine import load_full_data

Z_CAP = 10.0
SCORES_DIR = config.ENGINE_MIXED_DIR / "scores"


def shash_z(Z, gene_names, shash_params):
    Zs = np.clip(shash_correct_matrix(np.nan_to_num(Z, nan=0.0, posinf=10.0, neginf=-10.0), gene_names, shash_params),
                 -Z_CAP, Z_CAP)
    Zs[~np.isfinite(Z)] = np.nan
    return Zs


def main():
    SCORES_DIR.mkdir(parents=True, exist_ok=True)
    engine = NormativeModelEngineMixed.load(config.ENGINE_MIXED_DIR)
    gene_names = [g for g in engine.genes if engine.genes[g].ok]
    shash_params = load_shash_params(config.ENGINE_MIXED_DIR / "training_summary.csv")
    data = load_full_data()
    cols = [data["gene_col"][g] for g in gene_names]
    hc = data["is_hc"] & ~np.isin(data["batch"], list(data["small_hc_batches"]))
    for tag, rows in (("hc", np.where(hc)[0]), ("disease", np.where(~data["is_hc"])[0])):
        Z = np.clip(engine.score(data["X_raw"][rows], data["Y"][rows][:, cols], gene_names=gene_names), -Z_CAP, Z_CAP)
        np.save(SCORES_DIR / f"Z_{tag}.npy", Z)
        np.save(SCORES_DIR / f"Z_{tag}_shash.npy", shash_z(Z, gene_names, shash_params))
        pd.DataFrame({"sample": data["names"][rows], "batch": data["batch"][rows]}).to_csv(
            SCORES_DIR / f"{tag}_meta.csv", index=False)
    pickle.dump(gene_names, open(SCORES_DIR / "gene_names.pkl", "wb"))


if __name__ == "__main__":
    main()
