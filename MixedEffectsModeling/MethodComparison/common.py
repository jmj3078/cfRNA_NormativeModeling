import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import nbinom, norm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.calibration import bh_fdr_reject

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
FIG_DIR = HERE / "Figures"

ROOT_SEED = 20260914
LOG2FCS = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
Q_LEVELS = [0.05, 0.10, 0.20, 0.25]


def eval_universe():
    """Genes every arm can score: engine model-route, intersected with each OUTRIDER
    directory's per-fold surviving gene sets. Frozen to disk so every arm reads bytes
    rather than re-deriving its own universe."""
    path = CACHE / "eval_universe.txt"
    if path.exists():
        return path.read_text().split()

    summary = pd.read_csv(config.ENGINE_MIXED_DIR / "training_summary.csv")
    genes = set(summary.loc[summary["ok"].fillna(False) & (summary["route"] == "model"), "gene"])
    for d in (config.OUTRIDER_COMPARISON_DIR, config.OUTRIDER_HELD_OUT_DIR):
        for f in sorted(glob.glob(str(d / "cv_fold*_theta.csv"))):
            genes &= set(pd.read_csv(f)["gene"])
    universe = sorted(genes)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(universe))
    return universe


def draw_perturbed(universe, n_perturb, fold):
    """Gene names, not indices -- arms hold the universe in different orders."""
    rng = np.random.default_rng(np.random.SeedSequence(entropy=ROOT_SEED, spawn_key=(fold, n_perturb)))
    return sorted(rng.choice(len(universe), min(n_perturb, len(universe)), replace=False).tolist())


def size_factors(counts, loggeomeans):
    """DESeq2 median-of-ratios against a frozen train reference, so a test sample's depth
    never depends on the other test samples. counts: samples x genes."""
    finite = np.isfinite(loggeomeans)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.array([
            np.exp(np.median((np.log(c) - loggeomeans)[finite & (c > 0)]))
            for c in counts
        ])


def injection_reference(y_train, y_test):
    """m_ref[s, g] = train mean of depth-normalized counts, rescaled to each test sample's
    depth. Arm-independent: every arm injects the same counts at the same log2fc.
    y_train/y_test: samples x genes."""
    with np.errstate(divide="ignore"):
        loggeomeans = np.log(y_train).mean(axis=0)
    sf_tr = size_factors(y_train, loggeomeans)
    sf_te = size_factors(y_test, loggeomeans)
    base = (y_train / sf_tr[:, None]).mean(axis=0)
    return base[None, :] * sf_te[:, None]


def perturb(y, m_ref, pert_idx, log2fc):
    y_pert = y.copy()
    y_pert[:, pert_idx] = np.round(np.maximum(
        y[:, pert_idx] + m_ref[:, pert_idx] * (2 ** log2fc - 1), 0))
    return y_pert


def nb_pvals(y, mu, theta):
    """OUTRIDER's own two-sided tail test (OUTRIDER:::pVal). Var = mu + mu^2/theta."""
    n = np.clip(theta, 1e-8, None)
    p = n / (n + mu)
    cdf = nbinom.cdf(y, n, p)
    sf = nbinom.sf(y - 1, n, p)
    return np.clip(2 * np.minimum(cdf, sf), 0, 1)


def z_pvals(z):
    return 2 * norm.sf(np.abs(z))


def count_rejections(p, labels, qs=Q_LEVELS):
    """Per-sample BH across the shared evaluation universe. labels: bool over genes.
    Returns one dict per q with pooled tp/fp/fn/tn over samples."""
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


def summarize(df):
    g = df.groupby(["arm", "q", "log2fc"])[["tp", "fp", "fn", "tn"]].sum()
    tp, fp, fn, tn = g["tp"], g["fp"], g["fn"], g["tn"]
    return pd.DataFrame(dict(
        TP=tp, FP=fp, FN=fn, TN=tn,
        Precision=tp / (tp + fp), Recall=tp / (tp + fn),
        FDR=fp / (tp + fp), F1=2 * tp / (2 * tp + fp + fn),
        RejRate=(tp + fp) / (tp + fp + fn + tn),
        Prevalence=(tp + fn) / (tp + fp + fn + tn),
    )).round(4)
