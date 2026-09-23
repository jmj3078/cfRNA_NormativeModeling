"""Gene-level Jaccard over the FDR<0.05 significant SET (not a fixed top-k).

The significant set is what a paper actually reports, so this is the estimand the
control-composition claim is about. Set sizes differ across strata by design -- that
instability is part of the result, not a nuisance to be normalized away -- so the
expected-Jaccard-under-independence baseline is reported alongside.

Cache-first: reads the CSV if present, recomputes only with force=True.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import config

STAT_DIR = config.CTRL_COMP_DIR / "deseq2_stats"
OUT_CSV = config.CTRL_COMP_DIR / "gene_fdr_set_jaccard.csv"
DESIGNS = ["no_covariate", "ruvg_k1", "ruvg_k2", "ruvg_k3"]
N_GENES = 18892


def bh_reject(p, q=0.05):
    p = np.asarray(p)
    n = len(p)
    o = np.argsort(p)
    ok = p[o] <= (np.arange(1, n + 1) / n) * q
    r = np.zeros(n, dtype=bool)
    if ok.any():
        r[o[:np.max(np.where(ok)[0]) + 1]] = True
    return r


def sig_sets(tag, design, q=0.05):
    out = []
    for t in range(3):
        f = STAT_DIR / tag / f"T{t}_{design}.csv.gz"
        if not f.exists():
            return None
        s = pd.read_csv(f, index_col=0)["stat"].dropna()
        out.append(set(s.index[bh_reject(2 * norm.sf(np.abs(s.values)), q)]))
    return out


def pair_stats(sets):
    """Observed Jaccard and its independence expectation, averaged over the 3 pairs."""
    obs, exp = [], []
    for i, a in enumerate(sets):
        for b in sets[i + 1:]:
            if not (a | b):
                continue
            obs.append(len(a & b) / len(a | b))
            e = len(a) * len(b) / N_GENES
            exp.append(e / (len(a) + len(b) - e))
    if not obs:
        return np.nan, np.nan
    return float(np.mean(obs)), float(np.mean(exp))


def compute():
    rows = []
    for tag in sorted(p.name for p in STAT_DIR.iterdir() if p.is_dir()):
        disease = "Pancreatic Cancer" if tag.startswith("Pancreatic_Cancer") else "Pancreatitis"
        kind = "random" if "__null_" in tag else "tertile"
        axis = tag.split("__", 1)[1]
        for design in DESIGNS:
            s = sig_sets(tag, design)
            if s is None:
                continue
            j_obs, j_exp = pair_stats(s)
            sizes = list(map(len, s))
            rows.append(dict(tag=tag, disease=disease, kind=kind, axis=axis, design=design,
                             jaccard_fdr=j_obs, jaccard_exp=j_exp,
                             enrichment=j_obs / j_exp if j_exp and j_exp > 0 else np.nan,
                             n_sig_min=min(sizes), n_sig_med=int(np.median(sizes)),
                             n_sig_max=max(sizes)))
    return pd.DataFrame(rows)


def load(force=False):
    if OUT_CSV.exists() and not force:
        return pd.read_csv(OUT_CSV)
    df = compute()
    df.to_csv(OUT_CSV, index=False)
    return df


if __name__ == "__main__":
    df = load(force="--force" in sys.argv)
    pd.set_option("display.width", 200)
    print(df.groupby(["disease", "design", "kind"]).agg(
        n=("jaccard_fdr", "size"), J_obs=("jaccard_fdr", "mean"),
        J_exp=("jaccard_exp", "mean"), enrichment=("enrichment", "median"),
        nsig_min=("n_sig_min", "mean"), nsig_med=("n_sig_med", "mean"),
        nsig_max=("n_sig_max", "mean")).round(4).to_string())
