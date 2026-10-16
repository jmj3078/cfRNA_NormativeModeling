"""In-silico method comparison, MethodComparison/sweep.py (nb5) design, split in two steps.

freeze          write every injected-count cell once (universe genes, all held-out HC of each
                CV fold / LOBO batch). The arm-independent injection makes this possible.
score           score the FMMI_PIPELINE engine on the frozen cells (zero_shot), append its
                rows to config.INSILICO_CSV and keep its Z under runs/<pipe>/insilico/.
import          copy the production outlier-method arms (zero_shot) into the same csv.
Adaptation rows are added by eval/adaptation.py.

Run: FMMI_PIPELINE=<pipe> python eval/insilico_sweep.py freeze|score|import
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _FinalModelingMethodImprovement.config as config
from _FinalModelingMethodImprovement.engine.core.calibration import bh_fdr_reject
from _FinalModelingMethodImprovement.engine.core.dispersion_trend import load_trend
from _FinalModelingMethodImprovement.engine.core.marginal_rqr import marginal_nb_rqr
from _FinalModelingMethodImprovement.engine.core.shash import shash_transform_to_z
from _FinalModelingMethodImprovement.engine.validation.cv_engine import _mu_alpha
from _FinalModelingMethodImprovement.engine.validation.lobo_engine import load_full_data, safe_dir
from _FinalModelingMethodImprovement.eval.results import append

PROD_MC = config.BASELINE_MEM_DIR / "MethodComparison" / "cache"
SWEEP_DIR = config.INSILICO_DIR / "sweep"
ROOT_SEED = 20260914
N_PERTURB = 500
LOG2FCS = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
Q_LEVELS = [0.05, 0.10, 0.20, 0.25]
SPLITS = ("cv", "lobo")


def fold_spec(split):
    """[(key, tag, train_names, test_names)] in the production enumeration order."""
    if split == "cv":
        folds = json.loads(config.CV_FOLDS_PATH.read_text())
        return [(k, f"fold{k}", folds[k]["train"], folds[k]["test"]) for k in sorted(folds, key=int)]
    folds = json.loads(config.LOBO_FOLDS_PATH.read_text())
    return [(k, f"lobo_{safe_dir(k)}", folds[k]["train"], folds[k]["test"]) for k in sorted(folds)]


def draw_perturbed(universe, n_perturb, fold):
    rng = np.random.default_rng(np.random.SeedSequence(entropy=ROOT_SEED, spawn_key=(fold, n_perturb)))
    return sorted(rng.choice(len(universe), min(n_perturb, len(universe)), replace=False).tolist())


def size_factors(counts, loggeomeans):
    finite = np.isfinite(loggeomeans)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.array([np.exp(np.median((np.log(c) - loggeomeans)[finite & (c > 0)])) for c in counts])


def injection_reference(y_train, y_test):
    with np.errstate(divide="ignore"):
        loggeomeans = np.log(y_train).mean(axis=0)
    base = (y_train / size_factors(y_train, loggeomeans)[:, None]).mean(axis=0)
    return base[None, :] * size_factors(y_test, loggeomeans)[:, None]


def perturb(y, m_ref, pert_idx, log2fc):
    y_pert = y.copy()
    y_pert[:, pert_idx] = np.round(np.maximum(y[:, pert_idx] + m_ref[:, pert_idx] * (2 ** log2fc - 1), 0))
    return y_pert


def count_rejections(p, labels, qs=Q_LEVELS):
    out = []
    for q in qs:
        tp = fp = fn = tn = 0
        for row in p:
            finite = np.isfinite(row)
            reject = np.zeros(len(row), dtype=bool)
            reject[finite] = bh_fdr_reject(row[finite], q=q)
            tp += int((reject & labels).sum())
            fp += int((reject & ~labels).sum())
            fn += int((~reject & labels).sum())
            tn += int((~reject & ~labels).sum())
        out.append(dict(q=q, tp=tp, fp=fp, fn=fn, tn=tn))
    return out


def cell_path(split, tag, log2fc):
    return SWEEP_DIR / f"{split}_{tag}_log2fc{log2fc:g}.npz"


def universe():
    return (SWEEP_DIR / "universe.txt").read_text().split()


def freeze():
    """Cells are regenerated from raw counts and must equal the production matrices bit for bit."""
    SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    if not (SWEEP_DIR / "universe.txt").exists():
        shutil.copy2(PROD_MC / "eval_universe_nm.txt", SWEEP_DIR / "universe.txt")
    U = universe()
    data = load_full_data()
    row_of = {s: i for i, s in enumerate(data["names"])}
    Y_u = data["Y"][:, [data["gene_col"][g] for g in U]].astype(np.float64)
    for split in SPLITS:
        for seed_i, (key, tag, tr, te) in enumerate(fold_spec(split)):
            tr_rows, te_rows = [row_of[s] for s in tr], [row_of[s] for s in te]
            pert_idx = draw_perturbed(U, N_PERTURB, seed_i)
            labels = np.zeros(len(U), dtype=bool)
            labels[pert_idx] = True
            m_ref = injection_reference(Y_u[tr_rows], Y_u[te_rows])
            for f in LOG2FCS:
                y_pert = perturb(Y_u[te_rows], m_ref, pert_idx, f).astype(np.int32)
                prod = np.load(PROD_MC / "matrices" / f"{split}_{tag}_log2fc{f:g}.npz")
                assert prod["samples"].tolist() == te and prod["genes"].tolist() == U
                assert np.array_equal(prod["y_pert"], y_pert) and np.array_equal(prod["labels"], labels), (split, tag, f)
                np.savez_compressed(cell_path(split, tag, f), samples=np.array(te), labels=labels,
                                    y_pert=y_pert, pert_idx=np.array(pert_idx), key=key)
            print(f"froze {split}/{tag}", flush=True)


def engine_params(split, key, tag, U, tr_rows, te_rows, X_raw):
    if split == "cv":
        src = config.CV_MIXED_DIR / "fold_params"
        fits = pd.read_csv(src / f"model_fold{key}.csv").set_index("gene")
        shash = pd.read_csv(src / f"shash_model_fold{key}.csv").set_index("gene")
    else:
        src = config.LOBO_MIXED_DIR / safe_dir(key)
        fits = pd.read_csv(src / "model_fits.csv").set_index("gene")
        shash = pd.read_csv(src / "shash_params.csv").set_index("gene")
    alpha_fn = load_trend(config.DISPERSION_TREND_PATH)
    scaler = StandardScaler().fit(X_raw[tr_rows])
    Xa = np.column_stack([np.ones(len(te_rows)), scaler.transform(X_raw[te_rows])])
    mu_cols = [c for c in fits.columns if c.startswith("mu_coef_")]
    disp_cols = [c for c in fits.columns if c.startswith("disp_coef_")]
    params = {}
    for j, g in enumerate(U):
        if g not in fits.index or not bool(fits.loc[g, "ok"]):
            continue
        row = fits.loc[g]
        mu, alpha = _mu_alpha(Xa, row[mu_cols].to_numpy(dtype=float), row[disp_cols].to_numpy(dtype=float), row, alpha_fn)
        s = shash.loc[g] if g in shash.index and bool(shash.loc[g, "ok"]) else None
        params[j] = (mu, alpha, float(row["tau2"]), s)
    return params


def engine_z(seed_base, params, n_samples, n_genes, y_pert):
    z = np.full((n_samples, n_genes), np.nan)
    for j, (mu, alpha, tau2, s) in params.items():
        zg = marginal_nb_rqr(y_pert[:, j], mu, alpha, tau2, seed=seed_base + j)
        if s is not None:
            zg = shash_transform_to_z(zg, s.xi, s.eta, s.eps, s.delta)
        z[:, j] = np.clip(zg, -50, 50)
    return z


def score():
    U = universe()
    data = load_full_data()
    row_of = {s: i for i, s in enumerate(data["names"])}
    out_dir = config.RUN_DIR / "insilico"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for split in SPLITS:
        for seed_i, (key, tag, tr, te) in enumerate(fold_spec(split)):
            tr_rows, te_rows = [row_of[s] for s in tr], [row_of[s] for s in te]
            params = engine_params(split, key, tag, U, tr_rows, te_rows, data["X_raw"])
            for f in LOG2FCS:
                c = np.load(cell_path(split, tag, f))
                z = engine_z(1000 * (seed_i + 1), params, len(te_rows), len(U), c["y_pert"].astype(np.float64))
                p = 2 * norm.sf(np.abs(z))
                np.savez_compressed(out_dir / f"{split}_{tag}_log2fc{f:g}.npz", z=z.astype(np.float32))
                rows += [dict(split=split, method=config.PIPELINE, scenario="zero_shot", fold=key, log2fc=f,
                              n_samples=len(te_rows),
                              n_universe=len(U), n_perturbed=int(c["labels"].sum()), **r)
                         for r in count_rejections(p, c["labels"])]
            print(f"scored {split}/{tag}", flush=True)
    append(config.INSILICO_CSV, pd.DataFrame(rows))


def import_baseline():
    prod = pd.read_csv(PROD_MC / "sweep_nm_counts.csv")
    prod = prod[prod["arm"].isin(config.OUTLIER_ARMS)].rename(columns={"arm": "method"})
    append(config.INSILICO_CSV, prod.assign(scenario="zero_shot"))


if __name__ == "__main__":
    {"freeze": freeze, "score": score, "import": import_baseline}[sys.argv[1]]()
