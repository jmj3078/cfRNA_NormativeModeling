import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from MixedEffectsModeling.validation.ppc_simulate import simulate_marginal_nb

DIR = Path(__file__).parent
OUT = DIR / "outrider_cv_calibration_moments.csv"
PPC_N_REPS = 200

if OUT.exists():
    print(f"already cached -> {OUT}")
    raise SystemExit

mus, ys, thetas = [], [], []
for fi in range(5):
    mus.append(pd.read_csv(DIR / f"cv_fold{fi}_mu.csv", index_col=0))
    ys.append(pd.read_csv(DIR / f"cv_fold{fi}_y.csv", index_col=0))
    thetas.append(pd.read_csv(DIR / f"cv_fold{fi}_theta.csv").set_index("gene")["theta"])

common_genes = set(mus[0].columns)
for d in mus[1:]:
    common_genes &= set(d.columns)
common_genes = sorted(common_genes)

Y = pd.concat([d[common_genes] for d in ys], axis=0)
MU = pd.concat([d[common_genes] for d in mus], axis=0)
theta_rows = [pd.DataFrame(np.tile(thetas[fi][common_genes].values, (len(ys[fi]), 1)), columns=common_genes)
             for fi in range(5)]
TH = pd.concat(theta_rows, axis=0)
TH.index = Y.index

rows = []
for i, g in enumerate(common_genes):
    y = Y[g].values.astype(np.float64)
    mu = MU[g].values.astype(np.float64)
    alpha = 1.0 / np.maximum(TH[g].values.astype(np.float64), 1e-8)
    tau2 = np.zeros_like(mu)  # OUTRIDER has no mixed-effect batch term
    y_rep = simulate_marginal_nb(mu, alpha, tau2, PPC_N_REPS, seed=1000 + i)
    rows.append(dict(gene=g, obs_mean=y.mean(), pred_mean=y_rep.mean(1).mean(),
                     obs_var=y.std(), pred_var=y_rep.std(1).mean(),
                     obs_zero=(y == 0).mean(), pred_zero=(y_rep == 0).mean(1).mean()))
    if (i + 1) % 2000 == 0:
        print(i + 1, "/", len(common_genes))

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f"saved -> {OUT}")
