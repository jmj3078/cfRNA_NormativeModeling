"""How many samples does each group DE call rest on?

PancreaticNeoplasm vs its own controls in two scopes (CEDAR alone, CEDAR+BCC), refit with pydeseq2
under four designs. For every DE gene, the highest-expressing samples of the DE-direction group
(patients for up, controls for down) are removed one at a time and the Wald test is recomputed with
the same IRLS solver until p exceeds the run's BH-5% raw-p cut. Dispersion is either held at
DESeq2's value or re-estimated as MAP around the fitted trend.
"""
import pickle

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import fit_alpha_mle, irls_solver
from scipy.stats import f as f_dist
from scipy.stats import norm

import config

SCOPES = {"cedar": config.PANCREATIC_DEG_DIR, "pooled": config.PANCREATIC_DEG_POOLED_DIR}
LABEL = {"cedar": "CEDAR", "pooled": "CEDAR+BCC"}
DESIGNS = {"no_covariate": "~condition", "ruvg_k1": "~W_1+condition", "ruvg_k2": "~W_1+W_2+condition",
           "ruvg_k3": "~W_1+W_2+W_3+condition"}
SHORT = {"no_covariate": "none", "ruvg_k1": "k=1", "ruvg_k2": "k=2", "ruvg_k3": "k=3"}
SUBSET = "moore_b1__PancreaticNeoplasm"
MIN_COUNT_SUM, PADJ, MAX_STEPS, N_JOBS = 10, 0.05, 10, 12
K_BINS = ["1", "2", "3", "4-10"]


def load_scope(scope):
    d = SCOPES[scope]
    cache = pickle.load(open(d / "pool_cache.pkl", "rb"))
    W = pd.read_csv(d / "ruvg_W" / f"{SUBSET}.csv", index_col=0)
    obs = cache["obs"].set_index("sample").loc[W.index]
    cond = np.where(obs["phenotype"].values == "Pancreatic Cancer", "case", "HC")
    return cache["counts"].loc[W.index], cond, W, pd.Series(cache["gene_names"], index=cache["genes"])


def run_deseq2(scope, design):
    counts, cond, W, sym = load_scope(scope)
    meta = pd.DataFrame({"condition": cond}, index=counts.index)
    for i in range(int(design[-1]) if design.startswith("ruvg") else 0):
        meta[f"W_{i + 1}"] = W[f"W_{i + 1}"].values
    keep = counts.columns[counts.sum(axis=0) >= MIN_COUNT_SUM]
    inf = DefaultInference(n_cpus=N_JOBS)
    dds = DeseqDataSet(counts=counts[keep], metadata=meta, design=DESIGNS[design], inference=inf, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "case", "HC"], inference=inf, quiet=True)
    st.summary()
    return dds, st, counts, cond, sym, keep


def fit(scope, design, cache_dir):
    f = cache_dir / f"{scope}__{design}__fit.pkl"
    if f.is_file():
        return pickle.load(open(f, "rb"))
    dds, st, counts, cond, sym, keep = run_deseq2(scope, design)
    res = st.results_df
    sig = list(res.index[res["padj"].fillna(1.0) < PADJ])
    Y = pd.DataFrame(np.asarray(dds.X, dtype=float), index=counts.index, columns=keep)
    rf = getattr(dds, "counts_to_refit", None)
    if rf is not None:
        Y.loc[:, rf.var_names] = np.asarray(rf.X, dtype=float)
    X = dds.obsm["design_matrix"]
    out = dict(scope=scope, design=design, res=res, sig=sig, samples=list(counts.index),
               Y=Y[sig].values, X=X.values.astype(float),
               j=[k for k, c in enumerate(X.columns) if c.startswith("condition")][0],
               sf=dds.obs["size_factors"].values, disp=dds.var["dispersions"].reindex(sig).values,
               trend=dds.var["fitted_dispersions"].reindex(sig).values,
               prior_var=float(dds.uns["prior_disp_var"]),
               cooks=pd.DataFrame(dds.layers["cooks"], index=counts.index, columns=keep)[sig].values,
               is_case=cond == "case", sym=sym.reindex(sig).values,
               pstar=float(res.loc[sig, "pvalue"].max()) if sig else np.nan)
    pickle.dump(out, open(f, "wb"))
    return out


def shrunk_lfc(scope, design, cache_dir):
    f = cache_dir / f"{scope}__{design}__shrunk_lfc.csv"
    if f.is_file():
        return pd.read_csv(f, index_col=0)["log2FoldChange"]
    dds, st, *_ = run_deseq2(scope, design)
    st.lfc_shrink(coeff=[c for c in dds.varm["LFC"].columns if c.startswith("condition")][0])
    st.results_df[["log2FoldChange"]].to_csv(f)
    return st.results_df["log2FoldChange"]


def wald_z(y, X, sf, a, j):
    beta, mu, _, _ = irls_solver(y, sf, X, a)
    w = mu / (1 + mu * a)
    M = (X.T * w) @ X
    H = np.linalg.inv(M + np.eye(X.shape[1]) * 1e-6)
    return beta[j] / np.sqrt(H[j] @ M @ H[:, j])


def wald_p(y, X, sf, a, j):
    return 2 * norm.sf(abs(wald_z(y, X, sf, a, j)))


def wald_z_redisp(y, X, sf, trend, prior_var, j):
    _, mu, _, _ = irls_solver(y, sf, X, trend)
    a, _ = fit_alpha_mle(y, X, mu, trend, 1e-8, max(10.0, len(y)), prior_disp_var=prior_var, prior_reg=True)
    return wald_z(y, X, sf, a, j)


def wald_p_redisp(y, X, sf, trend, prior_var, j):
    return 2 * norm.sf(abs(wald_z_redisp(y, X, sf, trend, prior_var, j)))


def fragility(y, X, sf, a, j, pstar, trend, prior_var, up, is_case):
    cand = np.where(is_case if up else ~is_case)[0]
    order = [int(i) for i in cand[np.argsort(-(y / sf)[cand], kind="stable")][:MAX_STEPS]]
    k_fixed = k_redisp = np.nan
    m = np.ones(len(y), bool)
    for k, i in enumerate(order, 1):
        m[i] = False
        if np.isnan(k_fixed) and wald_p(y[m], X[m], sf[m], a, j) >= pstar:
            k_fixed = k
        if np.isnan(k_redisp) and wald_p_redisp(y[m], X[m], sf[m], trend, prior_var, j) >= pstar:
            k_redisp = k
        if np.isfinite(k_fixed) and np.isfinite(k_redisp):
            break
    return dict(p0=wald_p(y, X, sf, a, j), k_fixed=k_fixed, k_redisp=k_redisp, order=order)


def run_fragility(F, cache_dir):
    f = cache_dir / f"{F['scope']}__{F['design']}__fragility_toporder.pkl"
    if f.is_file():
        return pickle.load(open(f, "rb"))
    up = F["res"].loc[F["sig"], "log2FoldChange"].values > 0
    out = Parallel(n_jobs=N_JOBS)(delayed(fragility)(F["Y"][:, g], F["X"], F["sf"], F["disp"][g], F["j"],
                                                     F["pstar"], F["trend"][g], F["prior_var"], up[g], F["is_case"])
                                  for g in range(len(F["sig"])))
    pickle.dump(out, open(f, "wb"))
    return out


def k_bin(k):
    return "1" if k <= 1 else "2" if k <= 2 else "3" if k <= 3 else "4-10"


def trajectory(F, g, order):
    y, X, sf, j = F["Y"][:, g], F["X"], F["sf"], F["j"]
    args = (F["trend"][g], F["prior_var"], j)
    m = np.ones(len(y), bool)
    zs = [wald_z_redisp(y, X, sf, *args)]
    for i in order:
        m[i] = False
        zs.append(wald_z_redisp(y[m], X[m], sf[m], *args))
    return np.array(zs)


def gene_table(fits, frag):
    rows = []
    for key, R in frag.items():
        F = fits[key]
        n, p = F["X"].shape
        cut = f_dist.ppf(0.99, p, n - p)
        res = F["res"].loc[F["sig"]]
        for g, r in enumerate(R):
            first = r["order"][0]
            k = r["k_fixed"]
            removed = r["order"][:int(k)] if np.isfinite(k) else r["order"]
            rows.append(dict(scope=LABEL[F["scope"]], design=SHORT[F["design"]], gene=F["sig"][g],
                             symbol=F["sym"][g], log2FC=res["log2FoldChange"].iloc[g], padj=res["padj"].iloc[g],
                             p_deseq2=res["pvalue"].iloc[g], p_recomputed=r["p0"],
                             k_fixed_disp=k, k_reestimated_disp=r["k_redisp"],
                             removed_patients=int(F["is_case"][removed].sum()),
                             removed_controls=int((~F["is_case"][removed]).sum()),
                             first_removed=F["samples"][first],
                             first_removed_group="patient" if F["is_case"][first] else "control",
                             first_removed_cooks=F["cooks"][first, g], cooks_cutoff=cut,
                             first_removed_over_cutoff=bool(F["cooks"][first, g] > cut)))
    genes = pd.DataFrame(rows)
    genes["k_bin"] = genes.k_reestimated_disp.map(k_bin)
    return genes


def summary_table(genes):
    return genes.groupby(["scope", "design"], sort=False).apply(lambda t: pd.Series(dict(
        n_deg=len(t),
        lost_by_1=(t.k_fixed_disp <= 1).mean(), lost_by_3=(t.k_fixed_disp <= 3).mean(),
        lost_by_5=(t.k_fixed_disp <= 5).mean(), not_lost_by_max=t.k_fixed_disp.isna().mean(),
        lost_by_1_redisp=(t.k_reestimated_disp <= 1).mean(), lost_by_3_redisp=(t.k_reestimated_disp <= 3).mean(),
        first_removed_patient=(t.first_removed_group == "patient").mean(),
        first_removed_over_cooks=t.first_removed_over_cutoff.mean())), include_groups=False).reset_index()


def load_all(cache_dir):
    cache_dir.mkdir(parents=True, exist_ok=True)
    keys = [(s, d) for s in SCOPES for d in DESIGNS]
    fits = {k: fit(*k, cache_dir) for k in keys}
    shrunk = {k: shrunk_lfc(*k, cache_dir) for k in keys}
    frag = {k: run_fragility(F, cache_dir) for k, F in fits.items() if F["sig"]}
    genes = gene_table(fits, frag)
    traj = {k: [trajectory(fits[k], g, frag[k][g]["order"]) for g in range(len(fits[k]["sig"]))] for k in frag}
    cpm = {}
    for s in SCOPES:
        counts = load_scope(s)[0]
        cpm[s] = counts.div(counts.sum(axis=1), axis=0) * 1e6
    return dict(fits=fits, shrunk=shrunk, frag=frag, genes=genes, traj=traj, cpm=cpm)


def run_genes(R, key):
    g = R["genes"]
    return g[(g.scope == LABEL[key[0]]) & (g.design == SHORT[key[1]])]


def refit_check(R):
    rows = []
    for (s, d), F in R["fits"].items():
        stored = pd.read_csv(SCOPES[s] / SUBSET / f"{d}.csv.gz", index_col=0)
        rows.append(dict(scope=LABEL[s], design=SHORT[d], deg_refit=len(F["sig"]),
                         deg_stored=int((stored["padj"].fillna(1) < PADJ).sum()), pstar=F["pstar"]))
    g = R["genes"]
    max_dlogp = float(np.abs(np.log10(g.p_recomputed) - np.log10(g.p_deseq2)).max())
    return pd.DataFrame(rows), max_dlogp
