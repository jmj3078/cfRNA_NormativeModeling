"""Gene-level normative model branching engine.

Branch assignment
-----------------
  det_rate < low_det_thr (default 0.10)
      -> "logistic"  : L2 Logistic Regression
         z-score: Bernoulli RQR from P(detected | X)

  det_rate >= low_det_thr  AND  mean_count >= mean_count_min (default 2.0)
      -> count_model ("nbi" or "zinbi")  : NBI / ZINBI GAMLSS
         z-score: NBI / ZINBI full quantile residual

Outputs
-------
  score() returns a dict:
    "logistic" : (n_test, n_logistic)  -- low-det gene anomaly matrix  [SEPARATE]
    "count"    : (n_test, n_count)     -- high-det gene anomaly matrix
    "combined" : (n_test, n_all)       -- concatenated (logistic | count)
    "logistic_gene_names" / "count_gene_names" / "gene_names"

  save(directory) writes:
    genes.pkl, scaler.pkl, config.pkl
    training_summary.csv    -- per-gene branch, fit_ok, n_removed
    training_failures.csv   -- genes where fit_ok == False (for root-cause analysis)

Typical workflow
----------------
  engine = NormativeModelEngine(count_model="nbi")   # or "zinbi"
  engine.load_hc_data()
  engine.assign_branches()
  engine.train(verbose=True)
  engine.save("engine_state/")

  result = engine.score(X_disease_raw, Y_disease, gene_names=order)
  Z_low  = result["logistic"]   # low-det matrix (separate)
  Z_high = result["count"]      # high-det matrix
  Z_all  = result["combined"]
"""

from __future__ import annotations

import pickle
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import nbinom, norm
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning

import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
import rpy2.robjects.numpy2ri as rpyn

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="rpy2")
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

# ---- Constants -----------------------------------------------------------

BASE_DIR  = Path(__file__).resolve().parent
DATA_DIR  = BASE_DIR.parent / "OpenAccess_nfcore"
H5AD_PATH = DATA_DIR / "Merged_Processed_AnnData_with_Batch_Biases_QC_Status.h5ad"
R_HELPER  = BASE_DIR / "gamlss.r"

BIAS_COLUMNS: List[str] = [
    "log(Total Reads)",
    "Spliced Reads (%)",
    "gDNA Contamination (Intron/Exon)",
    "rRNA Fraction",
    "RNA Degradation (3' Bias)",
    "Platelet Score",
    "GC Bias",
    "Gene Length Bias",
    "NG80",
    "(NP80/NG80)",
]

LOW_DET_THR    = 0.10
DET_RATE_MIN   = 0.01
MEAN_COUNT_MIN = 2.0
LR_C           = 1.0
LR_MAX_ITER    = 1000


# ---- Pure-Python scoring helpers ----------------------------------------

def _bernoulli_rqr(p_detect: np.ndarray, y: np.ndarray,
                   seed: Optional[int] = None) -> np.ndarray:
    """Bernoulli RQR: z ~ N(0,1) when P(Y>0) = p_detect is correct."""
    detected = y > 0
    a  = np.where(detected, 1.0 - p_detect, 0.0)
    b  = np.where(detected, 1.0,             1.0 - p_detect)
    lo = np.clip(np.minimum(a, b), 1e-8, 1 - 1e-8)
    hi = np.clip(np.maximum(a, b), 1e-8, 1 - 1e-8)
    rng = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(lo, hi)).astype(np.float32)


def _nbi_rqr_from_coeffs(mu_coef: np.ndarray, sigma_coef: np.ndarray,
                           X_test: np.ndarray, y_test: np.ndarray,
                           seed: Optional[int] = None) -> np.ndarray:
    """gamlss NBI: log(mu)=X@mu_coef, log(sigma)=X@sigma_coef, theta=1/sigma."""
    n     = len(y_test)
    Xa    = np.column_stack([np.ones(n), X_test])
    mu    = np.exp(Xa @ mu_coef).clip(1e-4, 1e6)
    sigma = np.exp(Xa @ sigma_coef).clip(1e-8, 1e3)
    theta = (1.0 / sigma).clip(1e-4, 1e4)
    p_nb  = np.clip(theta / (theta + mu), 1e-8, 1 - 1e-8)
    yi    = np.asarray(y_test, dtype=int)
    a     = np.where(yi > 0, nbinom.cdf(yi - 1, n=theta, p=p_nb), 0.0)
    b     = nbinom.cdf(yi, n=theta, p=p_nb)
    lo    = np.clip(np.minimum(a, b), 1e-8, 1 - 1e-8)
    hi    = np.clip(np.maximum(a, b), 1e-8, 1 - 1e-8)
    rng   = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(lo, hi)).astype(np.float32)


def _zinbi_rqr_from_coeffs(mu_coef: np.ndarray, sigma_coef: np.ndarray,
                             nu_coef: np.ndarray,
                             X_test: np.ndarray, y_test: np.ndarray,
                             seed: Optional[int] = None) -> np.ndarray:
    """ZINBI full z-score: F_ZINBI(k) = nu + (1-nu)*F_NBI(k)."""
    n     = len(y_test)
    Xa    = np.column_stack([np.ones(n), X_test])
    mu    = np.exp(Xa @ mu_coef).clip(1e-4, 1e6)
    sigma = np.exp(Xa @ sigma_coef).clip(1e-8, 1e3)
    theta = (1.0 / sigma).clip(1e-4, 1e4)
    p_nb  = np.clip(theta / (theta + mu), 1e-8, 1 - 1e-8)
    # nu: intercept-only (1 element) or per-sample
    nu    = np.full(n, 1.0 / (1.0 + np.exp(-nu_coef[0]))) if len(nu_coef) == 1 \
            else 1.0 / (1.0 + np.exp(-Xa @ nu_coef))
    nu    = np.clip(nu, 1e-8, 1 - 1e-8)
    yi    = np.asarray(y_test, dtype=int)
    fn1   = np.where(yi > 0, nbinom.cdf(yi - 1, n=theta, p=p_nb), 0.0)
    fn    = nbinom.cdf(yi, n=theta, p=p_nb)
    a     = np.where(yi > 0, nu + (1 - nu) * fn1, 0.0)
    b     = nu + (1 - nu) * fn
    lo    = np.clip(np.minimum(a, b), 1e-8, 1 - 1e-8)
    hi    = np.clip(np.maximum(a, b), 1e-8, 1 - 1e-8)
    rng   = np.random.default_rng(seed)
    return norm.ppf(rng.uniform(lo, hi)).astype(np.float32)


# ---- rpy2 helpers -------------------------------------------------------

def _to_r_matrix(arr: np.ndarray, col_names: List[str]):
    with localconverter(ro.default_converter + rpyn.converter):
        r_mat = ro.conversion.py2rpy(np.ascontiguousarray(arr, dtype=np.float64))
    return ro.r["matrix"](r_mat, nrow=arr.shape[0], ncol=arr.shape[1],
                           dimnames=ro.r["list"](ro.NULL, ro.StrVector(col_names)))


def _to_r_vec(arr: np.ndarray):
    with localconverter(ro.default_converter + rpyn.converter):
        return ro.conversion.py2rpy(np.ascontiguousarray(arr, dtype=np.float64))


# ---- Gene record --------------------------------------------------------

@dataclass
class GeneRecord:
    name:     str
    branch:   str       # "logistic" | "nbi" | "zinbi"
    det_rate: float

    logistic_model: Optional[LogisticRegression] = None
    mu_coef:    Optional[np.ndarray] = None
    sigma_coef: Optional[np.ndarray] = None
    nu_coef:    Optional[np.ndarray] = None   # ZINBI only

    fit_ok:      bool = False
    n_removed:   int  = 0
    fail_reason: str  = ""


# ---- Engine -------------------------------------------------------------

class NormativeModelEngine:
    """
    Parameters
    ----------
    count_model : {"nbi", "zinbi"}
        GAMLSS family for high-detectability genes.
    zinbi_nu_formula : {"intercept", "full"}
        nu sub-model formula when count_model="zinbi".
        "intercept" regularizes nu as a global constant.
    low_det_thr : float
        Genes below this detection rate go to the logistic branch.
    nbi_outlier_z, nbi_max_iter, nbi_max_remove_frac
        Iterative outlier-removal settings for the count model.
    """

    def __init__(
        self,
        count_model:          str   = "nbi",
        zinbi_nu_formula:     str   = "intercept",
        low_det_thr:          float = LOW_DET_THR,
        det_rate_min:         float = DET_RATE_MIN,
        mean_count_min:       float = MEAN_COUNT_MIN,
        lr_C:                 float = LR_C,
        nbi_outlier_z:        float = 5.0,
        nbi_max_iter:         int   = 2,
        nbi_max_remove_frac:  float = 0.05,
    ):
        if count_model not in ("nbi", "zinbi"):
            raise ValueError("count_model must be 'nbi' or 'zinbi'")
        self.count_model         = count_model
        self.zinbi_nu_formula    = zinbi_nu_formula
        self.low_det_thr         = low_det_thr
        self.det_rate_min        = det_rate_min
        self.mean_count_min      = mean_count_min
        self.lr_C                = lr_C
        self.nbi_outlier_z       = nbi_outlier_z
        self.nbi_max_iter        = nbi_max_iter
        self.nbi_max_remove_frac = nbi_max_remove_frac

        self.X_hc_scaled:   Optional[np.ndarray]  = None
        self.Y_hc:          Optional[np.ndarray]  = None
        self.scaler:        Optional[StandardScaler] = None
        self.is_hc:         Optional[np.ndarray]  = None
        self.pc_gene_names: List[str]              = []
        self.pc_indices:    Optional[np.ndarray]  = None

        self.genes:          Dict[str, GeneRecord] = {}
        self.logistic_genes: List[str]             = []
        self.count_genes:    List[str]             = []

        self._r_nbi_fn   = None
        self._r_zinbi_fn = None

    # ---- Data loading ---------------------------------------------------

    def load_hc_data(self, h5ad_path: Path = H5AD_PATH):
        print("Loading HC data...")
        adata = sc.read_h5ad(h5ad_path)
        adata = adata[adata.obs["QC_Passed"] == True]
        adata = adata[adata.obs["Phenotype_Processed"].notna()]
        adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
        self.is_hc = (adata.obs["Phenotype_Processed"].astype(str) == "Healthy Control").values

        X_raw = adata.obs[BIAS_COLUMNS].values.astype(np.float64)
        self.scaler      = StandardScaler()
        self.X_hc_scaled = self.scaler.fit_transform(X_raw[self.is_hc])

        Y_raw = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
        self.Y_hc = np.round(Y_raw[self.is_hc]).astype(np.float64)

        is_pc = (adata.var["GeneType"] == "protein_coding").values
        self.pc_gene_names = adata.var_names[is_pc].tolist()
        self.pc_indices    = np.where(is_pc)[0]
        print(f"  HC={self.is_hc.sum()}  protein-coding={len(self.pc_gene_names)}")

    # ---- Branch assignment ---------------------------------------------

    def assign_branches(self) -> Dict[str, List[str]]:
        assert self.Y_hc is not None, "Call load_hc_data() first."
        Y_pc   = self.Y_hc[:, self.pc_indices]
        det_r  = (Y_pc > 0).mean(axis=0)
        mean_c = Y_pc.mean(axis=0)
        self.genes = {}

        for i, g in enumerate(self.pc_gene_names):
            dr, mc = float(det_r[i]), float(mean_c[i])
            if dr < self.det_rate_min:
                continue
            if dr < self.low_det_thr or mc < self.mean_count_min:
                branch = "logistic"
            else:
                branch = self.count_model
            self.genes[g] = GeneRecord(name=g, branch=branch, det_rate=dr)

        self.logistic_genes = [g for g, r in self.genes.items() if r.branch == "logistic"]
        self.count_genes    = [g for g, r in self.genes.items() if r.branch != "logistic"]
        print(f"Branches: logistic={len(self.logistic_genes)}"
              f"  {self.count_model}={len(self.count_genes)}"
              f"  (threshold={self.low_det_thr:.0%})")
        return {"logistic": self.logistic_genes, self.count_model: self.count_genes}

    # ---- R init --------------------------------------------------------

    def _init_r(self):
        if self._r_nbi_fn is None:
            ro.r(f'source("{R_HELPER}")')
            self._r_nbi_fn   = ro.globalenv["train_nbi_coeffs"]
            self._r_zinbi_fn = ro.globalenv["train_zinbi_coeffs"]

    # ---- Per-gene training ---------------------------------------------

    def _gene_y(self, g: str) -> np.ndarray:
        idx = self.pc_gene_names.index(g)
        return self.Y_hc[:, self.pc_indices[idx]]

    def _train_logistic(self, rec: GeneRecord):
        y_bin = (self._gene_y(rec.name) > 0).astype(int)
        if y_bin.sum() == 0 or y_bin.sum() == len(y_bin):
            rec.fail_reason = "no_variation"
            return
        lr = LogisticRegression(penalty="l2", C=self.lr_C, solver="lbfgs",
                                 max_iter=LR_MAX_ITER, random_state=42)
        lr.fit(self.X_hc_scaled, y_bin)
        rec.logistic_model = lr
        rec.fit_ok = True

    def _train_nbi(self, rec: GeneRecord):
        self._init_r()
        res = self._r_nbi_fn(
            _to_r_vec(self._gene_y(rec.name)),
            _to_r_matrix(self.X_hc_scaled, BIAS_COLUMNS),
            ro.IntVector([50]),
            ro.FloatVector([self.nbi_outlier_z]),
            ro.IntVector([self.nbi_max_iter]),
            ro.FloatVector([self.nbi_max_remove_frac]),
        )
        if res.rx2("success")[0]:
            rec.mu_coef    = np.array(res.rx2("mu_coef"))
            rec.sigma_coef = np.array(res.rx2("sigma_coef"))
            rec.n_removed  = int(res.rx2("n_removed")[0])
            rec.fit_ok     = True
        else:
            rec.fail_reason = str(res.rx2("msg")[0])

    def _train_zinbi(self, rec: GeneRecord):
        self._init_r()
        res = self._r_zinbi_fn(
            _to_r_vec(self._gene_y(rec.name)),
            _to_r_matrix(self.X_hc_scaled, BIAS_COLUMNS),
            ro.IntVector([50]),
            ro.FloatVector([self.nbi_outlier_z]),
            ro.IntVector([self.nbi_max_iter]),
            ro.FloatVector([self.nbi_max_remove_frac]),
            ro.StrVector([self.zinbi_nu_formula]),
        )
        if res.rx2("success")[0]:
            rec.mu_coef    = np.array(res.rx2("mu_coef"))
            rec.sigma_coef = np.array(res.rx2("sigma_coef"))
            rec.nu_coef    = np.array(res.rx2("nu_coef"))
            rec.n_removed  = int(res.rx2("n_removed")[0])
            rec.fit_ok     = True
        else:
            rec.fail_reason = str(res.rx2("msg")[0])

    # ---- Bulk training -------------------------------------------------

    def train(self, verbose: bool = True, limit: Optional[int] = None):
        assert self.genes, "Call assign_branches() first."
        all_genes = list(self.genes.keys())[:limit]
        n_log = sum(1 for g in all_genes if self.genes[g].branch == "logistic")
        print(f"Training {len(all_genes)} genes  "
              f"(logistic={n_log}, {self.count_model}={len(all_genes)-n_log})")

        for i, g in enumerate(all_genes):
            rec = self.genes[g]
            try:
                if   rec.branch == "logistic": self._train_logistic(rec)
                elif rec.branch == "zinbi":    self._train_zinbi(rec)
                else:                          self._train_nbi(rec)
            except Exception as exc:
                rec.fail_reason = str(exc)
                if verbose:
                    print(f"  [ERR] {g} ({rec.branch}): {exc}")
            if verbose and (i + 1) % 500 == 0:
                ok = sum(1 for g2 in all_genes if self.genes[g2].fit_ok)
                print(f"  [{i+1:5d}/{len(all_genes)}] fitted={ok}")

        ok = sum(1 for g in all_genes if self.genes[g].fit_ok)
        print(f"Training complete: {ok}/{len(all_genes)} succeeded.")

    # ---- Scoring -------------------------------------------------------

    def score(
        self,
        X_test_raw: np.ndarray,
        Y_test:     np.ndarray,
        gene_names: Optional[List[str]] = None,
        seed:       int = 42,
    ) -> Dict:
        """Score new samples. Returns dict with separate logistic/count matrices.

        Parameters
        ----------
        X_test_raw : (n_test, 10)   raw bias-metric covariates (unscaled)
        Y_test     : (n_test, n_genes)  raw count matrix
        gene_names : column order in Y_test (default: all fitted genes)

        Returns
        -------
        dict:
          "logistic"             : (n_test, n_logistic)  low-det matrix [SEPARATE]
          "count"                : (n_test, n_count)     high-det matrix
          "combined"             : (n_test, n_all)       logistic | count
          "logistic_gene_names"  : list
          "count_gene_names"     : list
          "gene_names"           : list (= combined order)
        """
        gene_names = gene_names or [g for g in self.genes if self.genes[g].fit_ok]
        X_test     = self.scaler.transform(X_test_raw.astype(np.float64))
        n_test, n_gene = len(X_test), len(gene_names)
        if Y_test.shape[1] != n_gene:
            raise ValueError(f"Y_test has {Y_test.shape[1]} columns, expected {n_gene}")

        Z = np.full((n_test, n_gene), np.nan, dtype=np.float32)
        for j, g in enumerate(gene_names):
            rec = self.genes.get(g)
            if rec is None or not rec.fit_ok:
                continue
            y_col = Y_test[:, j].astype(np.float64)
            try:
                if rec.branch == "logistic":
                    p = rec.logistic_model.predict_proba(X_test)[:, 1]
                    Z[:, j] = _bernoulli_rqr(p, y_col, seed + j)
                elif rec.branch == "zinbi":
                    Z[:, j] = _zinbi_rqr_from_coeffs(
                        rec.mu_coef, rec.sigma_coef, rec.nu_coef,
                        X_test, y_col, seed + j)
                else:
                    Z[:, j] = _nbi_rqr_from_coeffs(
                        rec.mu_coef, rec.sigma_coef, X_test, y_col, seed + j)
            except Exception:
                pass

        log_idx = [j for j, g in enumerate(gene_names)
                   if self.genes.get(g) and self.genes[g].branch == "logistic"]
        cnt_idx = [j for j, g in enumerate(gene_names)
                   if self.genes.get(g) and self.genes[g].branch != "logistic"]

        return {
            "logistic":            Z[:, log_idx] if log_idx else np.zeros((n_test, 0), np.float32),
            "count":               Z[:, cnt_idx] if cnt_idx else np.zeros((n_test, 0), np.float32),
            "combined":            Z,
            "gene_names":          gene_names,
            "logistic_gene_names": [gene_names[j] for j in log_idx],
            "count_gene_names":    [gene_names[j] for j in cnt_idx],
        }

    def to_dataframe(self, result: Dict, sample_ids=None, which="combined") -> pd.DataFrame:
        Z          = result[which]
        col_key    = "gene_names" if which == "combined" else f"{which}_gene_names"
        gene_names = result[col_key]
        sample_ids = sample_ids or [f"s{i}" for i in range(Z.shape[0])]
        return pd.DataFrame(Z, index=sample_ids, columns=gene_names)

    # ---- Diagnostics & persistence ------------------------------------

    def branch_summary(self) -> pd.DataFrame:
        rows = [{"gene": r.name, "branch": r.branch,
                 "det_rate": round(r.det_rate, 4),
                 "fit_ok": r.fit_ok, "n_removed": r.n_removed,
                 "fail_reason": r.fail_reason}
                for r in self.genes.values()]
        df = pd.DataFrame(rows).set_index("gene")
        agg = df.groupby("branch")[["fit_ok", "n_removed"]].agg(
            n_genes=("fit_ok", "count"), n_fitted=("fit_ok", "sum"),
            total_removed=("n_removed", "sum"))
        print(agg)
        return df

    def save(self, directory: str | Path):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)

        with open(directory / "genes.pkl",  "wb") as f: pickle.dump(self.genes,  f)
        with open(directory / "scaler.pkl", "wb") as f: pickle.dump(self.scaler, f)

        _SKIP = {"genes", "scaler", "X_hc_scaled", "Y_hc", "is_hc",
                 "pc_gene_names", "pc_indices", "logistic_genes", "count_genes"}
        cfg = {k: v for k, v in vars(self).items()
               if not k.startswith("_") and k not in _SKIP}
        with open(directory / "config.pkl", "wb") as f: pickle.dump(cfg, f)

        df = self.branch_summary()
        df.to_csv(directory / "training_summary.csv")

        df_fail = df[~df["fit_ok"]].copy()
        df_fail.to_csv(directory / "training_failures.csv")

        print(f"Engine saved to {directory}/")
        if len(df_fail):
            print(f"  Failures: {len(df_fail)} -> training_failures.csv")
            for branch, cnt in df_fail.groupby("branch").size().items():
                print(f"    {branch}: {cnt}")

    def load(self, directory: str | Path):
        directory = Path(directory)
        with open(directory / "genes.pkl",  "rb") as f: self.genes  = pickle.load(f)
        with open(directory / "scaler.pkl", "rb") as f: self.scaler = pickle.load(f)
        with open(directory / "config.pkl", "rb") as f: cfg = pickle.load(f)
        for k, v in cfg.items(): setattr(self, k, v)
        self.logistic_genes = [g for g, r in self.genes.items() if r.branch == "logistic"]
        self.count_genes    = [g for g, r in self.genes.items() if r.branch != "logistic"]
        n_ok = sum(1 for r in self.genes.values() if r.fit_ok)
        print(f"Engine loaded from {directory}/  ({n_ok} fitted genes)")
