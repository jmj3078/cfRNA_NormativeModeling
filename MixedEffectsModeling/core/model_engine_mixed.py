import csv
import json
import pickle
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.calibration import gene_shash_calibration
from MixedEffectsModeling.core.dispersion_trend import (
    build_trend, build_trend_from_fits, load_trend, save_trend,
)
from MixedEffectsModeling.core.eb_shrinkage import (
    estimate_slope_prior, load_disp_prior, save_disp_prior, squeeze_log_theta,
)
from MixedEffectsModeling.core.marginal_rqr import _poisson_rqr, marginal_nb_rqr
from MixedEffectsModeling.core.trend_report import trend_report

MP = config.SPIKE_PARAMS
EB = config.EB_PARAMS


@dataclass
class GeneRecordMixed:
    name: str
    route: str = ""
    stage: str = ""
    nz: int = 0
    ok: bool = False
    singular: bool = False
    tau2: float = 0.0
    mu_coef: np.ndarray = None
    disp_coef: np.ndarray = None
    disp_se: np.ndarray = None
    fail_reason: str = ""
    nbi_full_eb_reject_reason: str = ""
    nbi_intercept_eb_reject_reason: str = ""
    n_outliers: int = 0
    outlier_refit_failed: bool = False
    mean_hc: float = None
    trend_alpha: float = None
    log_theta_raw: float = None
    log_theta_eb: float = None
    cv_shash_ok: bool = None
    cv_shash_xi: float = None
    cv_shash_eta: float = None
    cv_shash_eps: float = None
    cv_shash_delta: float = None
    cv_shash_z_lo: float = None
    cv_shash_z_hi: float = None
    cv_raw_skew: float = None
    cv_raw_kurtosis: float = None
    cv_corrected_skew: float = None
    cv_corrected_kurtosis: float = None
    cv_naive_exceed: float = None
    cv_shash_exceed: float = None
    cv_naive_fdr_reject_rate: float = None
    cv_corr_fdr_reject_rate: float = None
    cv_obs_zero_frac: float = None
    cv_pred_zero_frac: float = None
    cv_zero_diff: float = None
    cv_pearson_chi2: float = None
    cv_obs_mean: float = None
    cv_pred_mean: float = None
    cv_mean_rel_err: float = None
    cv_obs_var: float = None
    cv_pred_var: float = None
    cv_var_rel_err: float = None


class NormativeModelEngineMixed:
    def __init__(self):
        self.X_hc_raw = None
        self.X_hc_scaled = None
        self.Y_hc = None
        self.scaler = None
        self.batch = None
        self.names = None
        self.pc_gene_names = []
        self._gene_col = {}
        self.genes = {}
        self.alpha_fn = None
        self.nz_a_max = None
        self.rare_glm = None
        self.disp_prior = None
        self.disp_tau_d2 = None
        self.trend_report = None
        self.trend_path = config.DISPERSION_TREND_PATH

    def load_hc_data(self, h5ad_path=config.H5AD_PATH):
        adata = sc.read_h5ad(h5ad_path)
        adata = adata[adata.obs["QC_Passed"] == True]
        adata = adata[adata.obs["Phenotype_Processed"].notna()]
        adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
        adata = adata[adata.obs["broad_protocol_category"] != "Exome-based (EB)"]
        is_hc = (adata.obs["Phenotype_Processed"].astype(str) == "Healthy Control").values
        batch_hc = adata.obs["Batch_ID"].astype(str).values[is_hc]

        # Small batch size samples can cause tau^2 explosion 
        bsize = pd.Series(batch_hc).value_counts()
        small = set(bsize.loc[lambda v: v < config.MIN_HC_BATCH_SIZE].index)
        keep = np.array([b not in small for b in batch_hc])

        X_raw = adata.obs[config.BIAS_COLUMNS].values.astype(np.float64)[is_hc][keep]
        self.X_hc_raw = X_raw
        self.scaler = StandardScaler()
        self.X_hc_scaled = self.scaler.fit_transform(X_raw)
        self.batch = batch_hc[keep]
        self.names = adata.obs_names.values[is_hc][keep]
        Y_raw = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
        self.Y_hc = np.round(Y_raw[is_hc][keep]).astype(np.float64)
        is_pc = (adata.var["GeneType"] == "protein_coding").values
        self.pc_gene_names = adata.var_names[is_pc].tolist()
        pc_indices = np.where(is_pc)[0]
        self._gene_col = {g: pc_indices[i] for i, g in enumerate(self.pc_gene_names)}

    def assign_routes(self):
        self.nz_a_max = config.NZ_A_MAX
        nz = (self.Y_hc[:, list(self._gene_col.values())] > 0).sum(axis=0)
        for i, g in enumerate(self.pc_gene_names):
            n = int(nz[i])
            route = "pool" if n < self.nz_a_max else "model"
            self.genes[g] = GeneRecordMixed(name=g, route=route, nz=n)

    def _write_r_inputs(self, tmp_dir):
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
        pd.DataFrame(self.X_hc_scaled, columns=config.BIAS_COLUMNS).to_csv(f"{tmp_dir}/X.csv.gz")
        pd.DataFrame({"Batch_ID": self.batch}).to_csv(f"{tmp_dir}/batch.csv.gz")

    def _write_gene_block(self, genes, tmp_dir, tag):
        Y = self.Y_hc[:, [self._gene_col[g] for g in genes]]
        pd.DataFrame(Y, columns=genes).to_csv(f"{tmp_dir}/Y_{tag}.csv.gz")
        pd.DataFrame({"gene": genes}).to_csv(f"{tmp_dir}/genes_{tag}.csv", index=False)

    def _run_glmm_fit(self, tmp_dir, tag, mode, out_csv, disp_prior_path=None):
        fit_params = Path(tmp_dir) / "fit_params.json"
        fit_params.write_text(json.dumps(config.FIT_PARAMS))
        cmd = [config.RSCRIPT, str(config.GLMM_FIT_R), "--x", f"{tmp_dir}/X.csv.gz",
               "--y", f"{tmp_dir}/Y_{tag}.csv.gz", "--batch", f"{tmp_dir}/batch.csv.gz",
               "--genes", f"{tmp_dir}/genes_{tag}.csv", "--trend", str(self.trend_path),
               "--mode", mode, "--fit-params", str(fit_params), "--out", out_csv]
        if disp_prior_path is not None:
            cmd += ["--disp-prior", str(disp_prior_path)]
        subprocess.run(cmd, check=True, cwd=str(config.GLMM_FIT_R.parent))

    def calib_genes(self, n=None, n_strata=None, seed=42):
        """Mean-expression-stratified subsample of model-route genes, so the EB
        prior scale is not dominated by the low-expression bulk."""
        n = EB["calib_n_genes"] if n is None else n
        n_strata = EB["calib_n_strata"] if n_strata is None else n_strata
        genes = [g for g, r in self.genes.items() if r.route == "model"]
        mean_hc = self.Y_hc[:, [self._gene_col[g] for g in genes]].mean(axis=0)
        rng = np.random.default_rng(seed)
        per = max(1, n // n_strata)
        picked = []
        for chunk in np.array_split(np.argsort(mean_hc), n_strata):
            k = min(per, len(chunk))
            picked.extend(chunk[rng.choice(len(chunk), k, replace=False)])
        return [genes[i] for i in sorted(picked)]

    def prepare_hyperparams(self, trend_path=None, disp_prior_path=None,
                            tmp_dir="/tmp/glmm_train", n_genes=None):
        trend_path = Path(trend_path or config.DISPERSION_TREND_PATH)
        disp_prior_path = Path(disp_prior_path or config.DISP_PRIOR_PATH)
        calib_path = trend_path.parent / "calib_fits.csv"
        if trend_path.exists() and disp_prior_path.exists():
            self.trend_path = trend_path
            self.alpha_fn = load_trend(trend_path)
            self.disp_prior = load_disp_prior(disp_prior_path)
            return False

        if calib_path.exists():
            calib = pd.read_csv(calib_path)
        else:
            genes = self.calib_genes(n=n_genes)
            self._write_r_inputs(tmp_dir)
            self._write_gene_block(genes, tmp_dir, "calib")
            out = f"{tmp_dir}/results_calib.csv"
            if not Path(out).exists():
                self._run_glmm_fit(tmp_dir, "calib", "calib", out)
            calib = pd.read_csv(out)
            calib_path.parent.mkdir(parents=True, exist_ok=True)
            calib.to_csv(calib_path, index=False)

        mean_hc = np.array([self.Y_hc[:, self._gene_col[g]].mean() for g in calib["gene"]])
        alpha_fit = np.exp(-calib["disp_coef_0"].to_numpy(dtype=float))
        ok = calib["ok"].to_numpy(dtype=bool)
        trend = build_trend_from_fits(mean_hc, alpha_fit, ok=ok)
        save_trend(trend, trend_path)
        self.trend_path = trend_path
        self.alpha_fn = load_trend(trend_path)

        # The calib fits predate the trend, so glmm_fit.R leaves trend_alpha NA. Fill it from
        # the trend just built: pcis_null.R reads this file as its RES and drops every gene
        # whose trend_alpha is not finite.
        calib["trend_alpha"] = [self.alpha_fn(m) for m in mean_hc]
        calib.to_csv(calib_path, index=False)

        self.disp_prior = estimate_slope_prior(calib_path)
        save_disp_prior(self.disp_prior, disp_prior_path)

        mom = build_trend(self.Y_hc[:, list(self._gene_col.values())], min_nz=MP["trend_min_nz"])
        self.trend_report = trend_report(mean_hc[ok], alpha_fit[ok], trend, self.disp_prior,
                                         trend_path.parent, mom_trend=mom)
        return True

    def apply_disp_squeeze(self):
        """EB squeeze of the dispersion intercept toward the lowess trend, pooled
        over both stages. Written back into disp_coef[0] so score() needs no
        stage branch: alpha = exp(-Xa @ disp_coef) already gives the squeezed
        constant for nbi_intercept_eb (NaN slopes -> 0) and the squeezed intercept plus real
        slopes for nbi_full_eb."""
        recs = [r for r in self.genes.values()
                if r.ok and r.stage in ("nbi_full_eb", "nbi_intercept_eb") and r.disp_coef is not None
                and r.trend_alpha is not None]
        if not recs:
            return
        hat = np.array([r.disp_coef[0] for r in recs], dtype=np.float64)
        se = np.array([r.disp_se[0] if r.disp_se is not None else np.nan for r in recs], dtype=np.float64)
        trend = np.array([-np.log(r.trend_alpha) for r in recs], dtype=np.float64)
        post, tau_d2 = squeeze_log_theta(hat, se, trend)
        self.disp_tau_d2 = tau_d2
        for r, h, p in zip(recs, hat, post):
            r.log_theta_raw = float(h) if np.isfinite(h) else None
            r.log_theta_eb = float(p)
            r.disp_coef[0] = p

    def train(self, limit=None, tmp_dir="/tmp/glmm_train", disp_prior_path=None):
        model_genes = [g for g, r in self.genes.items() if r.route == "model"][:limit]
        self._write_r_inputs(tmp_dir)
        self._write_gene_block(model_genes, tmp_dir, "model")
        self._run_glmm_fit(tmp_dir, "model", "cascade", f"{tmp_dir}/results.csv",
                           disp_prior_path=disp_prior_path or config.DISP_PRIOR_PATH)

        results = pd.read_csv(f"{tmp_dir}/results.csv").set_index("gene")
        mu_cols = [c for c in results.columns if c.startswith("mu_coef_")]
        disp_cols = [c for c in results.columns if c.startswith("disp_coef_")]
        se_cols = [c for c in results.columns if c.startswith("disp_se_")]
        txt = lambda row, k: row[k] if k in row.index and not pd.isna(row[k]) else ""
        for g, row in results.iterrows():
            rec = self.genes[g]
            rec.stage, rec.ok, rec.tau2 = row["stage"], bool(row["ok"]), float(row["tau2"])
            rec.singular = bool(row["singular"]) if not pd.isna(row["singular"]) else False
            rec.trend_alpha = float(row["trend_alpha"]) if not pd.isna(row["trend_alpha"]) else None
            rec.n_outliers = int(row["n_outliers"]) if not pd.isna(row["n_outliers"]) else 0
            rec.outlier_refit_failed = bool(row["outlier_refit_failed"])
            rec.mu_coef = row[mu_cols].values.astype(float)
            rec.disp_coef = row[disp_cols].values.astype(float)
            rec.disp_se = row[se_cols].values.astype(float)
            rec.fail_reason = txt(row, "fail_reason")
            rec.nbi_full_eb_reject_reason = txt(row, "nbi_full_eb_reject_reason")
            rec.nbi_intercept_eb_reject_reason = txt(row, "nbi_intercept_eb_reject_reason")
            if not rec.ok:
                rec.route = "excluded"

        self.apply_disp_squeeze()
        self.train_pool(tmp_dir=tmp_dir)
        self.fit_shash()

    def train_pool(self, tmp_dir="/tmp/glmm_train"):
        """Route "pool": one shared-beta pooled GLM (+ batch random intercept)
        fit jointly across all pool-route genes. Unused this round (nz_a_max
        defaults to 0, so no gene ever routes to "pool"); kept so re-enabling
        pooling later needs no changes here."""
        pool_genes = [g for g, r in self.genes.items() if r.route == "pool"]
        if not pool_genes:
            return
        Path(tmp_dir).mkdir(exist_ok=True)
        Y_pool = self.Y_hc[:, [self._gene_col[g] for g in pool_genes]]
        pd.DataFrame(Y_pool, columns=pool_genes).to_csv(f"{tmp_dir}/Y_pool.csv.gz")
        pd.DataFrame({"gene": pool_genes}).to_csv(f"{tmp_dir}/genes_pool.csv", index=False)

        subprocess.run([
            config.RSCRIPT, str(config.GLMM_FIT_POOL_R), "--x", f"{tmp_dir}/X.csv.gz", "--y", f"{tmp_dir}/Y_pool.csv.gz",
            "--batch", f"{tmp_dir}/batch.csv.gz", "--genes", f"{tmp_dir}/genes_pool.csv",
            "--rare-overdisp-thr", str(MP["rare_overdisp_thr"]), "--out", f"{tmp_dir}/results_pool.json",
        ], check=True, cwd=str(config.GLMM_FIT_POOL_R.parent))

        with open(f"{tmp_dir}/results_pool.json") as f:
            fit = json.load(f)

        if not fit["ok"]:
            for g in pool_genes:
                self.genes[g].route = "excluded"
                self.genes[g].fail_reason = "fit_pooled_glmm failed"
            return

        n_hc = self.X_hc_scaled.shape[0]
        self.rare_glm = {"family": fit["family"], "beta": np.asarray(fit["beta"]),
                         "alpha": fit["alpha"], "eps": 1.0 / (2 * n_hc),
                         "tau2": float(fit["tau2"]) if fit.get("tau2") is not None else 0.0,
                         "mult_lo": fit["mult_lo"], "mult_hi": fit["mult_hi"]}
        for g, m in zip(fit["gene"], fit["mean_hc"]):
            rec = self.genes[g]
            rec.mean_hc, rec.ok, rec.stage = float(m), True, "pool"

    SHASH_PROGRESS_FIELDS = [
        "gene", "cv_shash_ok", "cv_shash_xi", "cv_shash_eta", "cv_shash_eps", "cv_shash_delta",
        "cv_shash_z_lo", "cv_shash_z_hi", "cv_raw_skew", "cv_raw_kurtosis",
        "cv_corrected_skew", "cv_corrected_kurtosis", "cv_naive_exceed", "cv_shash_exceed",
        "cv_naive_fdr_reject_rate", "cv_corr_fdr_reject_rate",
    ]

    def _set_shash_fields(self, rec, calib):
        rec.cv_shash_ok, rec.cv_shash_xi = calib["cv_shash_ok"], calib["cv_shash_xi"]
        rec.cv_shash_eta, rec.cv_shash_eps = calib["cv_shash_eta"], calib["cv_shash_eps"]
        rec.cv_shash_delta = calib["cv_shash_delta"]
        rec.cv_shash_z_lo, rec.cv_shash_z_hi = calib["cv_shash_z_lo"], calib["cv_shash_z_hi"]
        rec.cv_raw_skew, rec.cv_raw_kurtosis = calib["cv_raw_skew"], calib["cv_raw_kurtosis"]
        rec.cv_corrected_skew, rec.cv_corrected_kurtosis = calib["cv_corrected_skew"], calib["cv_corrected_kurtosis"]
        rec.cv_naive_exceed, rec.cv_shash_exceed = calib["cv_naive_exceed"], calib["cv_shash_exceed"]
        rec.cv_naive_fdr_reject_rate = calib["cv_naive_fdr_reject_rate"]
        rec.cv_corr_fdr_reject_rate = calib["cv_corr_fdr_reject_rate"]

    def fit_shash(self, seed=42, progress_path=None, log_every=200):
        """Per-gene SHASH fit on THIS engine's own in-sample HC Z (scored under
        the params just trained above) -- the production null that disease Z
        is corrected against downstream. Deliberately in-sample: CV/LOBO exist
        to validate that this train-fit transform generalizes to genuinely
        held-out HC, not to supply the production params themselves (fitting
        on held-out Z here would mean the params never face a held-out check
        at all -- see validation/cv_engine.py, validation/lobo_engine.py).

        progress_path, if given, appends each gene's fit to a CSV the moment
        it's done (flushed immediately -- watchable with `tail -f`) and, on a
        rerun, skips genes already in that file instead of re-fitting them --
        this loop is ~19k independent Nelder-Mead calls with no other cheap
        checkpoint, so an interrupted run would otherwise lose all progress."""
        gene_names = [g for g, r in self.genes.items() if r.ok]
        if not gene_names:
            return
        Y_sub = self.Y_hc[:, [self._gene_col[g] for g in gene_names]]
        Z = self.score(self.X_hc_raw, Y_sub, gene_names=gene_names, seed=seed)

        done = set()
        if progress_path and Path(progress_path).exists():
            prev = pd.read_csv(progress_path)
            for _, prow in prev.iterrows():
                g = prow["gene"]
                done.add(g)
                rec = self.genes.get(g)
                if rec is not None:
                    self._set_shash_fields(rec, prow.to_dict())
            print(f"fit_shash: resuming, {len(done)}/{len(gene_names)} genes already in {progress_path}", flush=True)

        f = writer = None
        if progress_path:
            progress_path = Path(progress_path)
            progress_path.parent.mkdir(parents=True, exist_ok=True)
            write_header = not (progress_path.exists() and progress_path.stat().st_size > 0)
            f = open(progress_path, "a", newline="")
            writer = csv.DictWriter(f, fieldnames=self.SHASH_PROGRESS_FIELDS)
            if write_header:
                writer.writeheader()

        t0 = time.time()
        n_done, n_todo = len(done), len(gene_names) - len(done)
        try:
            for j, g in enumerate(gene_names):
                if g in done:
                    continue
                z = Z[:, j]
                z = z[np.isfinite(z)]
                if len(z) < 8:
                    continue
                fit = gene_shash_calibration(z)
                row = dict(gene=g, cv_shash_ok=fit["shash_ok"], cv_shash_xi=fit["shash_xi"],
                          cv_shash_eta=fit["shash_eta"], cv_shash_eps=fit["shash_eps"],
                          cv_shash_delta=fit["shash_delta"], cv_shash_z_lo=fit["z_lo"], cv_shash_z_hi=fit["z_hi"],
                          cv_raw_skew=fit["raw_skew"], cv_raw_kurtosis=fit["raw_kurtosis"],
                          cv_corrected_skew=fit["corrected_skew"], cv_corrected_kurtosis=fit["corrected_kurtosis"],
                          cv_naive_exceed=fit["naive_exceed"], cv_shash_exceed=fit["shash_exceed"],
                          cv_naive_fdr_reject_rate=fit["naive_fdr_reject_rate"],
                          cv_corr_fdr_reject_rate=fit["corr_fdr_reject_rate"])
                self._set_shash_fields(self.genes[g], row)
                if writer:
                    writer.writerow(row)
                    f.flush()
                n_done += 1
                if log_every and (n_done % log_every == 0 or n_done == len(gene_names)):
                    elapsed = time.time() - t0
                    rate = (n_done - len(done)) / max(elapsed, 1e-9)
                    eta_s = (len(gene_names) - n_done) / rate if rate > 0 else float("nan")
                    print(f"fit_shash: {n_done}/{len(gene_names)} ({100 * n_done / len(gene_names):.1f}%) "
                         f"elapsed={elapsed:.0f}s rate={rate:.1f} genes/s eta={eta_s:.0f}s", flush=True)
        finally:
            if f:
                f.close()

    def training_summary(self):
        rows = [dict(gene=r.name, route=r.route, stage=r.stage, nz=r.nz, ok=r.ok,
                    singular=r.singular, tau2=r.tau2, tau2_collapsed=bool(r.tau2 < 1e-4),
                    trend_alpha=r.trend_alpha, log_theta_raw=r.log_theta_raw, log_theta_eb=r.log_theta_eb,
                    n_outliers=r.n_outliers, outlier_refit_failed=r.outlier_refit_failed,
                    fail_reason=r.fail_reason,
                    nbi_full_eb_reject_reason=r.nbi_full_eb_reject_reason,
                    nbi_intercept_eb_reject_reason=r.nbi_intercept_eb_reject_reason,
                    cv_shash_ok=r.cv_shash_ok, cv_shash_xi=r.cv_shash_xi, cv_shash_eta=r.cv_shash_eta,
                    cv_shash_eps=r.cv_shash_eps, cv_shash_delta=r.cv_shash_delta,
                    cv_naive_exceed=r.cv_naive_exceed, cv_shash_exceed=r.cv_shash_exceed,
                    cv_naive_fdr_reject_rate=r.cv_naive_fdr_reject_rate, cv_corr_fdr_reject_rate=r.cv_corr_fdr_reject_rate)
               for r in self.genes.values()]
        return pd.DataFrame(rows).set_index("gene")

    def score(self, X_test_raw, Y_test, gene_names=None, seed=42, as_dict=False):
        gene_names = gene_names or [g for g in self.genes if self.genes[g].ok]
        X_test = self.scaler.transform(X_test_raw.astype(np.float64))
        Xa = np.column_stack([np.ones(len(X_test)), X_test])
        Z = np.full((len(X_test), len(gene_names)), np.nan, dtype=np.float32)
        for j, g in enumerate(gene_names):
            rec = self.genes.get(g)
            if rec is None or not rec.ok:
                continue
            if rec.route == "pool":
                mult = np.exp(X_test @ self.rare_glm["beta"][1:])
                if "mult_lo" in self.rare_glm and self.rare_glm["mult_lo"] is not None:
                    mult = np.clip(mult, self.rare_glm["mult_lo"], self.rare_glm["mult_hi"])
                mu = np.clip((rec.mean_hc + self.rare_glm["eps"]) * np.exp(self.rare_glm["beta"][0]) * mult, 1e-6, 1e8)
                if self.rare_glm["family"] == "poisson":
                    Z[:, j] = _poisson_rqr(Y_test[:, j].astype(np.float64), mu, seed + j)
                else:
                    Z[:, j] = marginal_nb_rqr(Y_test[:, j].astype(np.float64), mu, self.rare_glm["alpha"],
                                              self.rare_glm.get("tau2", 0.0), seed + j)
                continue
            mu = np.clip(np.exp(Xa @ np.nan_to_num(rec.mu_coef, nan=0.0)), 1e-6, 1e8)
            if not np.all(np.isnan(rec.disp_coef)):
                alpha = np.exp(-Xa @ np.nan_to_num(rec.disp_coef, nan=0.0))
            elif rec.trend_alpha is not None:
                alpha = np.full(len(X_test), rec.trend_alpha)
            else:
                alpha = np.full(len(X_test), self.alpha_fn(float(mu.mean())))
            Z[:, j] = marginal_nb_rqr(Y_test[:, j].astype(np.float64), mu, alpha, rec.tau2, seed + j)
        return Z if not as_dict else {"combined": Z, "gene_names": list(gene_names)}

    def save(self, directory):
        directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
        with open(directory / "genes.pkl", "wb") as f: pickle.dump(self.genes, f)
        with open(directory / "scaler.pkl", "wb") as f: pickle.dump(self.scaler, f)
        if self.rare_glm is not None:
            with open(directory / "rare_glm.pkl", "wb") as f: pickle.dump(self.rare_glm, f)
        if self.disp_prior is not None:
            save_disp_prior(self.disp_prior, directory / "disp_prior.json")
        (directory / "eb_meta.json").write_text(json.dumps(
            {"disp_tau_d2": self.disp_tau_d2, "disp_tau_d": None if self.disp_tau_d2 is None else self.disp_tau_d2 ** 0.5,
             "tau_slope": None if self.disp_prior is None else self.disp_prior["tau_slope"]}, indent=2))
        self.training_summary().to_csv(directory / "training_summary.csv")

    @classmethod
    def load(cls, directory):
        directory = Path(directory)
        engine = cls()
        with open(directory / "genes.pkl", "rb") as f: engine.genes = pickle.load(f)
        with open(directory / "scaler.pkl", "rb") as f: engine.scaler = pickle.load(f)
        rare_glm_path = directory / "rare_glm.pkl"
        if rare_glm_path.exists():
            with open(rare_glm_path, "rb") as f: engine.rare_glm = pickle.load(f)
        engine.disp_prior = load_disp_prior(directory / "disp_prior.json")
        eb_meta_path = directory / "eb_meta.json"
        if eb_meta_path.exists():
            engine.disp_tau_d2 = json.loads(eb_meta_path.read_text()).get("disp_tau_d2")
        trend_path = directory / "dispersion_trend.json"
        engine.trend_path = trend_path if trend_path.exists() else config.DISPERSION_TREND_PATH
        engine.alpha_fn = load_trend(engine.trend_path)
        return engine
