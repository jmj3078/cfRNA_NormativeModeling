"""Gate: the forked engine, loaded from the production engine state, must reproduce the
production in-sample HC Z bit for bit."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import _FinalModelingMethodImprovement.config as config
from _FinalModelingMethodImprovement.engine.core.model_engine_mixed import NormativeModelEngineMixed
from _FinalModelingMethodImprovement.engine.validation.lobo_engine import load_full_data

Z_CAP = 10.0
ref_dir = config.BASELINE_MEM_DIR / "Z_scores_mixed"

engine = NormativeModelEngineMixed.load(config.BASELINE_MEM_DIR / "engine_state_mixed")
gene_names = [g for g in engine.genes if engine.genes[g].ok]
data = load_full_data()
hc = data["is_hc"] & ~np.isin(data["batch"], list(data["small_hc_batches"]))
Y_hc = data["Y"][hc][:, [data["gene_col"][g] for g in gene_names]]
Z = np.clip(engine.score(data["X_raw"][hc], Y_hc, gene_names=gene_names), -Z_CAP, Z_CAP)

ref = np.load(ref_dir / "Z_hc.npy")
assert Z.shape == ref.shape, (Z.shape, ref.shape)
same = (Z == ref) | (np.isnan(Z) & np.isnan(ref))
print(f"bit-identical: {same.all()}  mismatched cells: {int((~same).sum())}/{same.size}")
assert same.all()
