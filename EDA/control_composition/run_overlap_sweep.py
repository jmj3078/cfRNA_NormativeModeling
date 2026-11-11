"""Control-overlap sweep (DESeq2 only).

Case group fixed (all Pancreatic Cancer of the scope); control subsets are half of the HC pool.
Per random permutation of the HC pool, N_WIN circular windows of size m are placed at offsets
whose pairwise circular distances cover every target overlap level (0-90% of m), so each
permutation contributes at least one subset pair per level and permutations are independent.
100% overlap is the identical subset (Jaccard 1 by construction) and is not fitted.

Run:  python EDA/control_composition/run_overlap_sweep.py --scope b1 [--n-perm 48]
Resumable: RUVg W and every DESeq2 stat file are written per subset and skipped if present.
"""
import argparse
import os
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
import config
import run_control_composition as rcc
import run_control_composition_deseq2 as rcd
from compute_fdr_set_jaccard import bh_reject

DISEASE = "Pancreatic Cancer"
N_WIN = 6
LEVELS = [i / 10 for i in range(10)]
SEED = 42


def out_dir(scope):
    return config.CTRL_COMP_OVERLAP_DIR / scope


def stat_path(scope, sid, design):
    return out_dir(scope) / "deseq2" / sid / f"{design}.csv.gz"


def ruler(P, m):
    tgt = {m - round(f * m) for f in LEVELS}
    for c in combinations(range(1, P), N_WIN - 1):
        s = (0,) + c
        if tgt <= {min(b - a, P - b + a) for a, b in combinations(s, 2)}:
            return s
    raise ValueError(f"no {N_WIN}-window ruler covers all levels for P={P}, m={m}")


def make_subsets(data, n_perm):
    pheno = data["obs"]["phenotype"].values
    hc, case = np.where(pheno == "Healthy Control")[0], np.where(pheno == DISEASE)[0]
    P, m = len(hc), len(hc) // 2
    starts = ruler(P, m)
    tgt = {round(f * m) for f in LEVELS}
    rng = np.random.default_rng(SEED)
    groups = []
    for r in range(n_perm):
        perm = rng.permutation(hc)
        wins = [perm[(s + np.arange(m)) % P] for s in starts]
        got = {len(np.intersect1d(a, b)) for a, b in combinations(wins, 2)}
        assert tgt <= got, f"perm {r}: overlaps {sorted(got)} miss {sorted(tgt - got)}"
        groups += [dict(sid=f"ovl_r{r:03d}_w{w}", perm=r, win=w, case=case, ctrl=c)
                   for w, c in enumerate(wins)]
    return groups


def subset_frame(data, groups):
    samples = data["obs"]["sample"].values
    return pd.DataFrame([dict(subset_id=g["sid"], perm=g["perm"], win=g["win"], sample=samples[i], role=role)
                         for g in groups for role, idx in (("case", g["case"]), ("hc", g["ctrl"]))
                         for i in idx])


def fit_subset(scope, data, counts, g):
    samples = data["obs"]["sample"].values
    idx = np.concatenate([g["case"], g["ctrl"]])
    W = rcc.load_W(g["sid"], samples[idx])
    sub_counts = counts.iloc[idx]
    condition = np.array(["disease"] * len(g["case"]) + ["HC"] * len(g["ctrl"]))
    for name, formula in rcd.DESIGNS.items():
        path = stat_path(scope, g["sid"], name)
        if path.exists():
            continue
        cond_df = pd.DataFrame({"condition": condition}, index=sub_counts.index)
        if name.startswith("ruvg_k"):
            for i in range(int(name.rsplit("_k", 1)[1])):
                cond_df[f"W_{i + 1}"] = W[:, i]
        s = rcd.fit_one(sub_counts, cond_df, formula)
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + ".part.gz")
        pd.DataFrame({"stat": s}, index=data["genes"]).to_csv(part)
        with open(part, "rb") as f:
            os.fsync(f.fileno())
        part.replace(path)


def sig_set(path, q=0.05):
    s = pd.read_csv(path, index_col=0)["stat"].dropna()
    return set(s.index[bh_reject(2 * norm.sf(np.abs(s.values)), q)])


def pair_table(scope, q=0.05):
    """One row per (permutation, subset pair, design). Only permutations whose N_WIN subsets
    have every design on disk are used; the cache is rebuilt when that set changes."""
    sub = pd.read_csv(out_dir(scope) / "subsets.csv")
    hc = sub[sub.role == "hc"].groupby("subset_id")["sample"].agg(frozenset)
    meta = sub.drop_duplicates("subset_id").set_index("subset_id")[["perm", "win"]]
    ok = meta[[all(stat_path(scope, s, d).exists() for d in rcd.DESIGNS) for s in meta.index]]
    perms = sorted(p for p, n in ok.groupby("perm").size().items() if n == N_WIN)
    cache = out_dir(scope) / "pair_jaccard.csv"
    if cache.exists():
        df = pd.read_csv(cache)
        if sorted(df.perm.unique()) == perms:
            return df
    m = len(hc.iloc[0])
    rows = []
    for p in perms:
        sids = sorted(meta.index[meta.perm == p])
        sig = {(s, d): sig_set(stat_path(scope, s, d), q) for s in sids for d in rcd.DESIGNS}
        for a, b in combinations(sids, 2):
            k = len(hc[a] & hc[b])
            for d in rcd.DESIGNS:
                ga, gb = sig[(a, d)], sig[(b, d)]
                u = ga | gb
                rows.append(dict(perm=p, sid_a=a, sid_b=b, design=d, shared_hc=k, n_hc=m,
                                 overlap=100 * k / m, n_sig_a=len(ga), n_sig_b=len(gb),
                                 jacc=len(ga & gb) / len(u) if u else np.nan))
    df = pd.DataFrame(rows, columns=["perm", "sid_a", "sid_b", "design", "shared_hc", "n_hc", "overlap",
                                     "n_sig_a", "n_sig_b", "jacc"])
    df.to_csv(cache, index=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scope", default=config.CTRL_COMP_DEFAULT_SCOPE, choices=sorted(config.CTRL_COMP_SCOPES))
    ap.add_argument("--n-perm", type=int, default=None)
    args = ap.parse_args()
    scope = args.scope
    n_perm = args.n_perm or config.CTRL_COMP_OVERLAP_N_PERM[scope]
    rcc.set_scope(scope)
    rcc.w_dir = lambda: out_dir(scope) / "ruvg_W"

    data = rcc.build_cache()
    counts = rcd.load_raw_counts(data)
    groups = make_subsets(data, n_perm)
    out_dir(scope).mkdir(parents=True, exist_ok=True)
    tbl = subset_frame(data, groups)
    tbl.to_csv(out_dir(scope) / "subsets.csv", index=False)
    rcc.log(f"overlap[{scope}]: {n_perm} perms x {N_WIN} windows, "
            f"m={tbl[tbl.role == 'hc'].groupby('subset_id').size().iloc[0]} HC, "
            f"{(tbl.role == 'case').sum() // len(groups)} cases")

    rcc.run_ruvg(data, tbl[["subset_id", "sample"]])
    todo = [g for g in groups if not all(stat_path(scope, g["sid"], d).exists() for d in rcd.DESIGNS)]
    rcc.log(f"overlap[{scope}]: {len(todo)}/{len(groups)} subsets to fit")
    t0 = time.time()
    for i, g in enumerate(todo, 1):
        fit_subset(scope, data, counts, g)
        if i % N_WIN == 0 or i == len(todo):
            rate = (time.time() - t0) / i
            rcc.log(f"overlap[{scope}]: {i}/{len(todo)} subsets ({rate:.0f}s/subset, "
                    f"eta {(len(todo) - i) * rate / 60:.0f} min)")
    rcc.log(f"overlap[{scope}]: done")


if __name__ == "__main__":
    main()
