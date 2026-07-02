"""NZ-gated 4-Phase normative model engine (normative-v2).

Gating (Phase 1, on HC nonzero sample count NZ):
  Route A (NZ < nz_a_max)                -> rare pooling (Phase 4)
  Route B (nz_a_max <= NZ < nz_bc_split)  -> penalized mean-only NB, dispersion
                                             borrowed from the covariate-free
                                             mean-dispersion trend (Phase 3)
  Route C (NZ >= nz_bc_split)             -> full NBI GAMLSS (mu and sigma both
                                             regressed on covariates) (Phase 2)

Demotion chain (information-loss priority: NBI -> mean-only -> pooled):
  Route C candidate fails (full AND intercept-only NBI) -> try Route B.
    Route B also fails                                  -> EXCLUDE (no rare fallback:
                                                            NZ>=nz_bc_split means this
                                                            gene was never rare-pooling
                                                            material).
  Route B candidate fails (full AND intercept-only mean-only NB) -> Route A (rare
                                                                     pooling, always
                                                                     succeeds).
  Route A candidate -> rare pooling directly.

Discrete covariates (batch, study, ...) are NEVER used as model covariates, only
the 10 continuous BIAS_COLUMNS. Batch is only used for CV fold stratification
elsewhere.
"""

import pickle
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as rpyn
import scanpy as sc
import statsmodels.api as sm
from rpy2.robjects.conversion import localconverter
from scipy.sparse import issparse
from scipy.stats import nbinom, norm, poisson
from sklearn.preprocessing import StandardScaler
from statsmodels.discrete.discrete_model import NegativeBinomial

warnings.filterwarnings("ignore", category=UserWarning, module="rpy2")

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
import config
from dispersion_trend import build_trend, load_trend, save_trend

MP2 = config.MODELING_PARAMS_V2


# ---- Pure-Python RQR helpers (mirrors model_engine.py) ------------------

def _poisson_rqr(y, mu, seed=None):
    y = np.asarray(y)
    lo = np.where(y > 0, poisson.cdf(y - 1, mu), 0.0)
    hi = poisson.cdf(y, mu)
    lo = np.clip(lo, 1e-12, 1 - 1e-12); hi = np.clip(hi, 1e-12, 1 - 1e-12)
    rng = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(np.minimum(lo, hi), np.maximum(lo, hi))).astype(np.float32)


def _nb_rqr(y, mu, alpha, seed=None):
    y = np.asarray(y)
    n = 1.0 / alpha
    p = np.clip(n / (n + mu), 1e-12, 1 - 1e-12)
    lo = np.where(y > 0, nbinom.cdf(y - 1, n, p), 0.0)
    hi = nbinom.cdf(y, n, p)
    lo = np.clip(lo, 1e-12, 1 - 1e-12); hi = np.clip(hi, 1e-12, 1 - 1e-12)
    rng = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(np.minimum(lo, hi), np.maximum(lo, hi))).astype(np.float32)


def _nbi_rqr_from_coeffs(mu_coef, sigma_coef, X_test, y_test, seed=None):
    n = len(y_test)
    Xa = np.column_stack([np.ones(n), X_test])
    mu = np.exp(Xa @ mu_coef).clip(1e-4, 1e6)
    sigma = np.exp(Xa @ sigma_coef).clip(1e-8, 1e3)
    theta = (1.0 / sigma).clip(1e-4, 1e4)
    p_nb = np.clip(theta / (theta + mu), 1e-8, 1 - 1e-8)
    yi = np.asarray(y_test, dtype=int)
    a = np.where(yi > 0, nbinom.cdf(yi - 1, n=theta, p=p_nb), 0.0)
    b = nbinom.cdf(yi, n=theta, p=p_nb)
    lo = np.clip(np.minimum(a, b), 1e-8, 1 - 1e-8)
    hi = np.clip(np.maximum(a, b), 1e-8, 1 - 1e-8)
    rng = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(lo, hi)).astype(np.float32)


def _nbi_rqr_intercept(mu_coef, sigma_coef, y_test, seed=None):
    n = len(y_test)
    mu = np.full(n, np.exp(mu_coef[0])).clip(1e-4, 1e6)
    sigma = np.full(n, np.exp(sigma_coef[0])).clip(1e-8, 1e3)
    theta = (1.0 / sigma).clip(1e-4, 1e4)
    p_nb = np.clip(theta / (theta + mu), 1e-8, 1 - 1e-8)
    yi = np.asarray(y_test, dtype=int)
    a = np.where(yi > 0, nbinom.cdf(yi - 1, n=theta, p=p_nb), 0.0)
    b = nbinom.cdf(yi, n=theta, p=p_nb)
    lo = np.clip(np.minimum(a, b), 1e-8, 1 - 1e-8)
    hi = np.clip(np.maximum(a, b), 1e-8, 1 - 1e-8)
    rng = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(lo, hi)).astype(np.float32)


# ---- rpy2 helpers ---------------------------------------------------------

def _to_r_matrix(arr, col_names):
    with localconverter(ro.default_converter + rpyn.converter):
        r_mat = ro.conversion.py2rpy(np.ascontiguousarray(arr, dtype=np.float64))
    return ro.r["matrix"](r_mat, nrow=arr.shape[0], ncol=arr.shape[1],
                          dimnames=ro.r["list"](ro.NULL, ro.StrVector(col_names)))


def _to_r_vec(arr):
    with localconverter(ro.default_converter + rpyn.converter):
        return ro.conversion.py2rpy(np.ascontiguousarray(arr, dtype=np.float64))


# ---- Route B: pure-Python penalized mean-only NB (fixed dispersion) ------

def _ridge_nb_irls(y, X, alpha, lam, max_iter=100, tol=1e-8):
    """L2-penalized NB2(alpha fixed) mean-model IRLS. X includes intercept col 0
    (unpenalized); columns 1: are penalized by lam. Returns (beta, converged)."""
    n, p = X.shape
    beta = np.zeros(p)
    beta[0] = np.log(max(y.mean(), 1e-3))
    pen = np.eye(p) * lam
    pen[0, 0] = 0.0
    for _ in range(max_iter):
        eta = X @ beta
        mu = np.clip(np.exp(eta), 1e-6, 1e8)
        w = mu / (1.0 + alpha * mu)
        z = eta + (y - mu) / np.clip(mu, 1e-6, None)
        WX = X * w[:, None]
        XtWX = X.T @ WX + pen
        XtWz = X.T @ (w * z)
        try:
            beta_new = np.linalg.solve(XtWX, XtWz)
        except np.linalg.LinAlgError:
            return beta, False
        if not np.all(np.isfinite(beta_new)):
            return beta, False
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            return beta, True
        beta = beta_new
    return beta, True


def _nb_deviance(y, mu, alpha):
    mu = np.clip(mu, 1e-8, None)
    y_safe = np.where(y > 0, y, 1e-8)
    term = y * np.log(y_safe / mu) - (y + 1.0 / alpha) * np.log((1 + alpha * y) / (1 + alpha * mu))
    return float(2 * term.sum())


def fit_route_b_gene(y_train, X_train, alpha_fn, lam, outlier_z, max_iter, max_remove_frac,
                     beta_explode_thr=None, gaic_k=None):
    """Penalized mean-only NB with fixed dispersion from the trend. Fits full
    (covariate) and intercept-only models, with |z|>5 outlier removal (<=5%),
    returns dict(success, beta, beta_null, alpha, gaic_full, gaic_null, chosen,
    n_removed, fail_reason)."""
    beta_explode_thr = MP2["beta_explode_thr"] if beta_explode_thr is None else beta_explode_thr
    gaic_k = MP2["gaic_k"] if gaic_k is None else gaic_k
    n = len(y_train)
    Xa = np.column_stack([np.ones(n), X_train])
    keep = np.ones(n, dtype=bool)
    n_removed = 0
    beta = beta_null = None
    alpha_g = alpha_fn(float(y_train.mean()))

    for _ in range(max_iter):
        y_k, X_k = y_train[keep], Xa[keep]
        alpha_g = alpha_fn(float(y_k.mean()))
        beta, ok = _ridge_nb_irls(y_k, X_k, alpha_g, lam)
        if not ok or not np.all(np.isfinite(beta)):
            return dict(success=False, fail_reason="irls_diverged", n_removed=n_removed)
        mu_k = np.clip(np.exp(X_k @ beta), 1e-6, 1e8)
        z_k = _nb_rqr(y_k, mu_k, alpha_g, seed=0)
        outlier = np.isfinite(z_k) & (np.abs(z_k) > outlier_z)
        if not outlier.any():
            break
        if outlier.sum() / n > max_remove_frac:
            break
        idx_keep = np.where(keep)[0]
        keep[idx_keep[outlier]] = False
        n_removed += int(outlier.sum())

    if beta is None or not np.all(np.isfinite(beta)):
        return dict(success=False, fail_reason="fit_failed", n_removed=n_removed)

    y_k, X_k = y_train[keep], Xa[keep]
    mu_full = np.clip(np.exp(X_k @ beta), 1e-6, 1e8)
    dev_full = _nb_deviance(y_k, mu_full, alpha_g)
    edf_full = X_k.shape[1]
    gaic_full = dev_full + gaic_k * edf_full

    mu0 = np.full(1, np.log(max(y_k.mean(), 1e-3)))
    beta_null = np.concatenate([mu0, np.zeros(Xa.shape[1] - 1)])
    mu_null = np.full(len(y_k), np.exp(mu0[0])).clip(1e-6, 1e8)
    dev_null = _nb_deviance(y_k, mu_null, alpha_g)
    gaic_null = dev_null + gaic_k * 1

    beta_max = float(np.abs(beta[1:]).max())
    if not np.isfinite(gaic_full) or beta_max > beta_explode_thr:
        chosen = "intercept"
    else:
        chosen = "full" if gaic_full < gaic_null else "intercept"

    return dict(success=True, beta=beta, beta_null=beta_null, alpha=alpha_g,
               gaic_full=gaic_full, gaic_null=gaic_null, chosen=chosen,
               beta_max=beta_max, n_removed=n_removed, fail_reason="")


# ---- Gene record ----------------------------------------------------------

@dataclass
class GeneRecordV2:
    name: str
    initial_route: str      # "A" | "B" | "C"  (Phase 1 gating)
    route: str = ""          # final route actually used: "A" | "B" | "C" | "excluded"
    nz: int = 0
    attempted: bool = False  # True once train() has processed this gene (vs. skipped by --limit)

    # Route C (NBI)
    mu_coef: np.ndarray = None
    sigma_coef: np.ndarray = None
    nbi_chosen: str = ""     # "full" | "intercept"
    nbi_explode: str = ""    # "" | "mu" | "sigma" | "mu+sigma"  (which submodel triggered demotion)

    # Route B (mean-only NB, fixed dispersion)
    beta: np.ndarray = None
    alpha: float = None
    b_chosen: str = ""       # "full" | "intercept"

    # Route A (rare pooling)
    mean_hc: float = None

    fit_ok: bool = False
    n_removed: int = 0
    fail_reason: str = ""


# ---- Engine ---------------------------------------------------------------

class NormativeModelEngineV2:
    def __init__(self, nz_a_max=None, nz_bc_split=None, trend_min_nz=None,
                ridge_lambda_mu=None, outlier_z=None, max_outlier_iter=None,
                max_remove_frac=None, beta_explode_thr=None, gaic_k=None,
                rare_overdisp_thr=None, rare_z_cap=None):
        self.nz_a_max = nz_a_max or MP2["nz_a_max"]
        self.nz_bc_split = nz_bc_split or MP2["nz_bc_split"]
        self.trend_min_nz = trend_min_nz or MP2["trend_min_nz"]
        self.ridge_lambda_mu = ridge_lambda_mu or MP2["ridge_lambda_mu"]
        self.gaic_k = gaic_k or MP2["gaic_k"]
        self.outlier_z = outlier_z or MP2["outlier_z"]
        self.max_outlier_iter = max_outlier_iter or MP2["max_outlier_iter"]
        self.max_remove_frac = max_remove_frac or MP2["max_remove_frac"]
        self.beta_explode_thr = beta_explode_thr or MP2["beta_explode_thr"]
        self.rare_overdisp_thr = rare_overdisp_thr or MP2["rare_overdisp_thr"]
        self.rare_z_cap = rare_z_cap or MP2["rare_z_cap"]

        self.X_hc_scaled = None
        self.Y_hc = None
        self.scaler = None
        self.is_hc = None
        self.pc_gene_names = []
        self.pc_indices = None

        self.genes = {}
        self.alpha_fn = None
        self.rare_glm = None

        self._r_nbi_fn = None
        self._r_nbi_null_fn = None

    # ---- Data loading -----------------------------------------------------

    def load_hc_data(self, h5ad_path=config.H5AD_PATH):
        print("Loading HC data...")
        adata = sc.read_h5ad(h5ad_path)
        adata = adata[adata.obs["QC_Passed"] == True]
        adata = adata[adata.obs["Phenotype_Processed"].notna()]
        adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
        adata = adata[adata.obs["broad_protocol_category"] != "Exome-based (EB)"]
        self.is_hc = (adata.obs["Phenotype_Processed"].astype(str) == "Healthy Control").values

        X_raw = adata.obs[config.BIAS_COLUMNS].values.astype(np.float64)
        self.scaler = StandardScaler()
        self.X_hc_scaled = self.scaler.fit_transform(X_raw[self.is_hc])

        Y_raw = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
        self.Y_hc = np.round(Y_raw[self.is_hc]).astype(np.float64)

        is_pc = (adata.var["GeneType"] == "protein_coding").values
        self.pc_gene_names = adata.var_names[is_pc].tolist()
        self.pc_indices = np.where(is_pc)[0]
        print(f"  HC={self.is_hc.sum()}  protein-coding={len(self.pc_gene_names)}")

    # ---- Phase 0: dispersion trend -----------------------------------------

    def build_dispersion_trend(self):
        Y_pc = self.Y_hc[:, self.pc_indices]
        trend = build_trend(Y_pc, min_nz=self.trend_min_nz)
        save_trend(trend)
        self.alpha_fn = load_trend()
        print(f"Dispersion trend built: n_reliable={trend['n_reliable']} "
              f"n_bins={trend['n_bins_used']}")

    # ---- Phase 1: gating ----------------------------------------------------

    def assign_routes(self):
        assert self.Y_hc is not None, "Call load_hc_data() first."
        Y_pc = self.Y_hc[:, self.pc_indices]
        nz = (Y_pc > 0).sum(axis=0)
        self.genes = {}
        for i, g in enumerate(self.pc_gene_names):
            n = int(nz[i])
            if n < self.nz_a_max:
                route = "A"
            elif n < self.nz_bc_split:
                route = "B"
            else:
                route = "C"
            self.genes[g] = GeneRecordV2(name=g, initial_route=route, nz=n)
        counts = pd.Series([r.initial_route for r in self.genes.values()]).value_counts()
        print(f"Phase 1 gating: A={counts.get('A',0)} B={counts.get('B',0)} "
              f"C={counts.get('C',0)}  (nz_a_max={self.nz_a_max}, nz_bc_split={self.nz_bc_split})")
        return counts

    # ---- R init -------------------------------------------------------------

    def _init_r(self):
        if self._r_nbi_fn is None:
            ro.r(f'source("{config.R_HELPER}")')
            self._r_nbi_fn = ro.globalenv["train_nbi_coeffs"]
            self._r_nbi_null_fn = ro.globalenv["train_nbi_coeffs_null"]

    def _gene_y(self, g):
        idx = self.pc_gene_names.index(g)
        return self.Y_hc[:, self.pc_indices[idx]]

    # ---- Phase 2: Route C (NBI, full vs intercept, GAIC) ---------------------

    def _train_route_c(self, rec):
        self._init_r()
        y = self._gene_y(rec.name)
        try:
            res_full = self._r_nbi_fn(
                _to_r_vec(y), _to_r_matrix(self.X_hc_scaled, config.BIAS_COLUMNS),
                ro.IntVector([50]), ro.FloatVector([self.outlier_z]),
                ro.IntVector([self.max_outlier_iter]), ro.FloatVector([self.max_remove_frac]),
                ro.FloatVector([self.ridge_lambda_mu]),
            )
        except Exception as exc:
            rec.fail_reason = f"nbi_full_error:{exc}"
            return False

        full_ok = bool(res_full.rx2("success")[0])
        beta_full = np.array(res_full.rx2("mu_coef")) if full_ok else None
        sigma_full = np.array(res_full.rx2("sigma_coef")) if full_ok else None
        mu_explode = full_ok and float(np.abs(beta_full[1:]).max()) > self.beta_explode_thr
        sigma_explode = full_ok and float(np.abs(sigma_full[1:]).max()) > self.beta_explode_thr
        full_explode = mu_explode or sigma_explode

        try:
            res_null = self._r_nbi_null_fn(_to_r_vec(y), ro.IntVector([50]))
        except Exception as exc:
            rec.fail_reason = f"nbi_null_error:{exc}"
            res_null = None
        null_ok = res_null is not None and bool(res_null.rx2("success")[0])

        if full_ok and full_explode:
            rec.nbi_explode = "+".join(t for t, e in [("mu", mu_explode), ("sigma", sigma_explode)] if e)

        if not null_ok and (not full_ok or full_explode):
            rec.fail_reason = rec.fail_reason or "nbi_full_and_null_failed"
            return False

        if full_ok and not full_explode and null_ok:
            gaic_full = float(res_full.rx2("gaic")[0])
            gaic_null = float(res_null.rx2("gaic")[0])
            chosen = "full" if gaic_full < gaic_null else "intercept"
        elif full_ok and not full_explode:
            chosen = "full"
        else:
            chosen = "intercept"

        if chosen == "full":
            rec.mu_coef = np.array(res_full.rx2("mu_coef"))
            rec.sigma_coef = np.array(res_full.rx2("sigma_coef"))
            rec.n_removed = int(res_full.rx2("n_removed")[0])
        else:
            rec.mu_coef = np.array(res_null.rx2("mu_coef"))
            rec.sigma_coef = np.array(res_null.rx2("sigma_coef"))
        rec.nbi_chosen = chosen
        rec.route = "C"
        rec.fit_ok = True
        return True

    # ---- Phase 3: Route B (penalized mean-only NB) ---------------------------

    def _train_route_b(self, rec):
        y = self._gene_y(rec.name)
        res = fit_route_b_gene(y, self.X_hc_scaled, self.alpha_fn, self.ridge_lambda_mu,
                               self.outlier_z, self.max_outlier_iter, self.max_remove_frac,
                               beta_explode_thr=self.beta_explode_thr, gaic_k=self.gaic_k)
        if not res["success"]:
            rec.fail_reason = res["fail_reason"]
            return False
        rec.beta = res["beta"] if res["chosen"] == "full" else res["beta_null"]
        rec.alpha = res["alpha"]
        rec.b_chosen = res["chosen"]
        rec.n_removed = res["n_removed"]
        rec.route = "B"
        rec.fit_ok = True
        return True

    # ---- Phase 4: Route A (rare pooling, pooled GLM, always succeeds) --------

    def train_rare(self, gene_list):
        if not gene_list:
            return
        n_hc = self.X_hc_scaled.shape[0]
        eps = 1.0 / (2 * n_hc)
        cols = [self.pc_indices[self.pc_gene_names.index(g)] for g in gene_list]
        Y_rare = self.Y_hc[:, cols]
        mean_hc = Y_rare.mean(axis=0)
        for g, m in zip(gene_list, mean_hc):
            self.genes[g].mean_hc = float(m)
            self.genes[g].route = "A"
            self.genes[g].fit_ok = True
            self.genes[g].attempted = True

        n_rare = len(gene_list)
        sample_idx = np.repeat(np.arange(n_hc), n_rare)
        gene_idx = np.tile(np.arange(n_rare), n_hc)
        Xc = np.column_stack([np.ones(n_hc * n_rare), self.X_hc_scaled[sample_idx]])
        y = Y_rare[sample_idx, gene_idx]
        offset = np.log(mean_hc[gene_idx] + eps)
        pois = sm.GLM(y, Xc, family=sm.families.Poisson(), offset=offset).fit()
        ratio = float(pois.deviance / pois.df_resid)
        if ratio <= self.rare_overdisp_thr:
            family, beta, alpha = "poisson", np.asarray(pois.params), None
        else:
            nb = NegativeBinomial(y, Xc, offset=offset).fit(disp=False)
            family, beta, alpha = "negbin", np.asarray(nb.params[:-1]), float(nb.params[-1])
        self.rare_glm = {"family": family, "beta": beta, "alpha": alpha,
                         "eps": eps, "overdisp_ratio": ratio}
        print(f"Route A (rare pooling): {n_rare} genes pooled, family={family}, "
              f"deviance/df={ratio:.3f}")

    # ---- Bulk training with demotion chain -----------------------------------

    def train(self, verbose=True, limit=None):
        assert self.genes, "Call assign_routes() first."
        if self.alpha_fn is None:
            self.build_dispersion_trend()

        all_genes = list(self.genes.keys())[:limit]
        route_c = [g for g in all_genes if self.genes[g].initial_route == "C"]
        route_b = [g for g in all_genes if self.genes[g].initial_route == "B"]
        route_a = [g for g in all_genes if self.genes[g].initial_route == "A"]
        print(f"Training: C-candidates={len(route_c)}  B-candidates={len(route_b)}  "
              f"A-candidates={len(route_a)}")

        # Phase 2: Route C candidates
        demoted_to_b = []
        for i, g in enumerate(route_c):
            rec = self.genes[g]
            rec.attempted = True
            try:
                ok = self._train_route_c(rec)
            except Exception as exc:
                ok = False
                rec.fail_reason = str(exc)
            if not ok:
                demoted_to_b.append(g)
            if verbose and (i + 1) % 500 == 0:
                print(f"  [Route C {i+1:5d}/{len(route_c)}] demoted_so_far={len(demoted_to_b)}")
        print(f"Phase 2 (Route C): {len(route_c)-len(demoted_to_b)} fitted, "
              f"{len(demoted_to_b)} demoted to Route B")

        # Phase 3: Route B candidates (original + demoted from C)
        b_pool = route_b + demoted_to_b
        demoted_to_a = []
        excluded = []
        for i, g in enumerate(b_pool):
            rec = self.genes[g]
            rec.attempted = True
            try:
                ok = self._train_route_b(rec)
            except Exception as exc:
                ok = False
                rec.fail_reason = str(exc)
            if not ok:
                if rec.initial_route == "C":
                    excluded.append(g)
                    rec.route = "excluded"
                else:
                    demoted_to_a.append(g)
            if verbose and (i + 1) % 500 == 0:
                print(f"  [Route B {i+1:5d}/{len(b_pool)}] "
                      f"demoted_to_A={len(demoted_to_a)} excluded={len(excluded)}")
        print(f"Phase 3 (Route B): {len(b_pool)-len(demoted_to_a)-len(excluded)} fitted, "
              f"{len(demoted_to_a)} demoted to Route A, {len(excluded)} EXCLUDED "
              f"(were Route C originally)")

        # Phase 4: Route A candidates (original + demoted from B)
        a_pool = route_a + demoted_to_a
        self.train_rare(a_pool)

        n_fitted = sum(1 for r in self.genes.values() if r.fit_ok)
        n_excluded = sum(1 for r in self.genes.values() if r.route == "excluded")
        print(f"Training complete: {n_fitted} fitted, {n_excluded} excluded, "
              f"total={len(self.genes)}")

    # ---- Scoring --------------------------------------------------------------

    def _rare_z(self, rec, X_test, y_col, seed):
        g = self.rare_glm
        Xc = np.column_stack([np.ones(len(X_test)), X_test])
        mu = (rec.mean_hc + g["eps"]) * np.exp(Xc @ g["beta"])
        mu = np.clip(mu, 1e-12, 1e8)
        if g["family"] == "poisson":
            z = _poisson_rqr(y_col, mu, seed)
        else:
            z = _nb_rqr(y_col, mu, g["alpha"], seed)
        return np.clip(z, -self.rare_z_cap, self.rare_z_cap).astype(np.float32)

    def score(self, X_test_raw, Y_test, gene_names=None, seed=42):
        gene_names = gene_names or [g for g in self.genes if self.genes[g].fit_ok]
        X_test = self.scaler.transform(X_test_raw.astype(np.float64))
        n_test, n_gene = len(X_test), len(gene_names)
        if Y_test.shape[1] != n_gene:
            raise ValueError(f"Y_test has {Y_test.shape[1]} columns, expected {n_gene}")

        Z = np.full((n_test, n_gene), np.nan, dtype=np.float32)
        Xa = np.column_stack([np.ones(n_test), X_test])
        for j, g in enumerate(gene_names):
            rec = self.genes.get(g)
            if rec is None or not rec.fit_ok:
                continue
            y_col = Y_test[:, j].astype(np.float64)
            try:
                if rec.route == "A":
                    Z[:, j] = self._rare_z(rec, X_test, y_col, seed + j)
                elif rec.route == "B":
                    if rec.b_chosen == "full":
                        mu = np.clip(np.exp(Xa @ rec.beta), 1e-6, 1e8)
                    else:
                        mu = np.full(n_test, np.exp(rec.beta[0])).clip(1e-6, 1e8)
                    Z[:, j] = _nb_rqr(y_col, mu, rec.alpha, seed + j)
                elif rec.route == "C":
                    if rec.nbi_chosen == "full":
                        Z[:, j] = _nbi_rqr_from_coeffs(rec.mu_coef, rec.sigma_coef, X_test, y_col, seed + j)
                    else:
                        Z[:, j] = _nbi_rqr_intercept(rec.mu_coef, rec.sigma_coef, y_col, seed + j)
            except Exception:
                pass
        return Z

    # ---- Diagnostics & persistence --------------------------------------------

    def training_summary(self, attempted_only=True):
        """attempted_only=True (default) reports only genes train() actually processed
        -- relevant when train(limit=N) was used for a smoke test, since the remaining
        genes were gated (initial_route set) but never attempted."""
        recs = [r for r in self.genes.values() if r.attempted] if attempted_only else self.genes.values()
        rows = [{"gene": r.name, "initial_route": r.initial_route, "route": r.route,
                 "nz": r.nz, "fit_ok": r.fit_ok, "attempted": r.attempted,
                 "nbi_chosen": r.nbi_chosen, "nbi_explode": r.nbi_explode, "b_chosen": r.b_chosen,
                 "n_removed": r.n_removed, "fail_reason": r.fail_reason}
                for r in recs]
        return pd.DataFrame(rows).set_index("gene")

    def save(self, directory):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        with open(directory / "genes.pkl", "wb") as f: pickle.dump(self.genes, f)
        with open(directory / "scaler.pkl", "wb") as f: pickle.dump(self.scaler, f)
        if self.rare_glm is not None:
            with open(directory / "rare_glm.pkl", "wb") as f: pickle.dump(self.rare_glm, f)
        _SKIP = {"genes", "scaler", "X_hc_scaled", "Y_hc", "is_hc", "rare_glm",
                 "pc_gene_names", "pc_indices", "alpha_fn"}
        cfg = {k: v for k, v in vars(self).items()
               if not k.startswith("_") and k not in _SKIP}
        with open(directory / "config.pkl", "wb") as f: pickle.dump(cfg, f)
        df = self.training_summary(attempted_only=False)
        df.to_csv(directory / "training_summary.csv")
        print(f"Engine saved to {directory}/")

    @classmethod
    def load(cls, directory):
        directory = Path(directory)
        with open(directory / "config.pkl", "rb") as f: cfg = pickle.load(f)
        engine = cls(**{k: v for k, v in cfg.items()
                        if k in cls.__init__.__code__.co_varnames})
        with open(directory / "genes.pkl", "rb") as f: engine.genes = pickle.load(f)
        with open(directory / "scaler.pkl", "rb") as f: engine.scaler = pickle.load(f)
        rare_glm_path = directory / "rare_glm.pkl"
        if rare_glm_path.exists():
            with open(rare_glm_path, "rb") as f: engine.rare_glm = pickle.load(f)
        engine.alpha_fn = load_trend()
        n_ok = sum(1 for r in engine.genes.values() if r.fit_ok)
        print(f"Engine loaded from {directory}/  ({n_ok} fitted genes)")
        return engine
