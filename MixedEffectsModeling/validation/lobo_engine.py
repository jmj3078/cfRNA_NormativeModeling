import argparse
import json
import pickle
import re
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.dispersion_trend import load_trend
from MixedEffectsModeling.core.eb_shrinkage import squeeze_log_theta
from MixedEffectsModeling.core.marginal_rqr import _poisson_rqr, marginal_nb_rqr
from MixedEffectsModeling.core.ood_filter import MahalanobisFilter, RangeFilter
from MixedEffectsModeling.core.shash import fit_and_correct, shash_transform_to_z
from MixedEffectsModeling.validation.cv_engine import squeeze_fold

MP = config.SPIKE_PARAMS


def _decode_cat(grp):
    codes = grp["codes"][()]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in grp["categories"][()]])
    out = np.full(len(codes), "", dtype=object)
    valid = codes >= 0
    out[valid] = cats[codes[valid]]
    return out


def load_full_data(h5ad_path=config.H5AD_PATH):
    """Same QC filters as NormativeModelEngineMixed.load_hc_data, but keeps every
    phenotype (not just HC) so held-out disease samples can be scored too.

    Reads X + needed obs/var columns directly via h5py instead of sc.read_h5ad --
    even backed='r' only defers X, not `layers` (~7.8GB of unused RUVg/TMM/FPKM/tpm/
    scaled variants on this h5ad), which was OOM-ing this pipeline."""
    with h5py.File(h5ad_path, "r") as f:
        qc_passed = f["obs/QC_Passed"][()]
        phenotype_full = _decode_cat(f["obs/Phenotype_Processed"])
        protocol = _decode_cat(f["obs/broad_protocol_category"])
        keep = qc_passed & (phenotype_full != "") & (phenotype_full != "Unknown") & (protocol != "Exome-based (EB)")

        X_raw = np.column_stack([f[f"obs/{c}"][()] for c in config.BIAS_COLUMNS]).astype(np.float64)[keep]
        names = f["obs/_index"][()].astype(str)[keep]
        batch = _decode_cat(f[f"obs/{config.STRATIFY_COL}"])[keep]
        phenotype = phenotype_full[keep]

        shape = tuple(f["X"].attrs["shape"])
        Xsp = csr_matrix((f["X/data"][()], f["X/indices"][()], f["X/indptr"][()]), shape=shape)
        Y = np.round(Xsp[keep].toarray()).astype(np.float64)

        gene_names_all = f["var/_index"][()].astype(str)
        is_pc = _decode_cat(f["var/GeneType"]) == "protein_coding"
        pc_gene_names = gene_names_all[is_pc].tolist()
        pc_idx = np.where(is_pc)[0]
        gene_col = {g: pc_idx[i] for i, g in enumerate(pc_gene_names)}

    is_hc = phenotype == "Healthy Control"
    bsize = pd.Series(batch[is_hc]).value_counts()
    small_hc_batches = set(bsize.loc[lambda v: v < config.MIN_HC_BATCH_SIZE].index)
    return dict(X_raw=X_raw, Y=Y, names=names, batch=batch, is_hc=is_hc, phenotype=phenotype,
               small_hc_batches=small_hc_batches, gene_col=gene_col)


def usable_batches(data):
    """Batch_IDs with >=1 HC sample surviving the small-HC-batch drop. tier A =
    disease samples also present in that batch, tier B = HC-only (noise-floor only)."""
    is_hc, batch = data["is_hc"], data["batch"]
    hc_batches = sorted(set(batch[is_hc]) - data["small_hc_batches"])
    dis_counts = pd.Series(batch[~is_hc]).value_counts()
    rows = [dict(batch=b, n_dis=int(dis_counts.get(b, 0)),
                tier="A" if dis_counts.get(b, 0) > 0 else "B") for b in hc_batches]
    return pd.DataFrame(rows)


def _refit_model_genes(tr_idx, data, genes, stage_of, scaler, tmp, disp_prior_path, cache_path=None):
    """cache_path persists the post-squeeze GLMM fit for this batch's
    leave-one-batch-out train set -- a rerun with the cache present skips
    Rscript entirely, same cache-first contract as validation/cv_engine.py."""
    if cache_path and cache_path.exists():
        return pd.read_csv(cache_path).set_index("gene")
    Xs_tr = scaler.transform(data["X_raw"][tr_idx])
    Y_tr = data["Y"][tr_idx][:, [data["gene_col"][g] for g in genes]]
    pd.DataFrame(Xs_tr, columns=config.BIAS_COLUMNS).to_csv(f"{tmp}/X.csv.gz")
    pd.DataFrame(Y_tr, columns=genes).to_csv(f"{tmp}/Y.csv.gz")
    pd.DataFrame({"Batch_ID": data["batch"][tr_idx]}).to_csv(f"{tmp}/batch.csv.gz")
    pd.DataFrame({"gene": genes, "stage": [stage_of[g] for g in genes]}).to_csv(f"{tmp}/genes.csv", index=False)
    fit_params = Path(tmp) / "fit_params.json"
    fit_params.write_text(json.dumps(config.FIT_PARAMS))
    cmd = [config.RSCRIPT, str(config.GLMM_FIT_R), "--x", f"{tmp}/X.csv.gz", "--y", f"{tmp}/Y.csv.gz",
          "--batch", f"{tmp}/batch.csv.gz", "--genes", f"{tmp}/genes.csv",
          "--trend", str(config.DISPERSION_TREND_PATH), "--fit-params", str(fit_params),
          "--mode", "fixed_stage", "--out", f"{tmp}/res.csv"]
    if disp_prior_path is not None:
        cmd += ["--disp-prior", str(disp_prior_path)]
    subprocess.run(cmd, check=True, cwd=str(config.GLMM_FIT_R.parent))
    fits = squeeze_fold(pd.read_csv(f"{tmp}/res.csv").set_index("gene"))
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        fits.to_csv(cache_path)
    return fits


def _score_model_genes(fits, genes, idx, data, scaler, alpha_fn, seed=42, report=True):
    """Score `genes` at `idx` (either the held-out batch or, for the SHASH
    train-fit below, the same train_idx the GLMM above was just fit on).
    report=False skips the per-gene fold_info rows (only needed once, for the
    held-out call)."""
    Xa = np.column_stack([np.ones(len(idx)), scaler.transform(data["X_raw"][idx])])
    zdict, rows = {}, []
    for g in genes:
        if g not in fits.index:
            if report:
                rows.append(dict(gene=g, route="model", ok=False, fail_reason="fold_output_missing"))
            continue
        row = fits.loc[g]
        ok = bool(row["ok"])
        if report:
            rows.append(dict(gene=g, route="model", stage=row["stage"], ok=ok,
                             singular=bool(row["singular"]) if not pd.isna(row["singular"]) else None,
                             tau2=float(row["tau2"]) if not pd.isna(row["tau2"]) else np.nan,
                             fail_reason=row["fail_reason"] if not pd.isna(row["fail_reason"]) else ""))
        if not ok:
            continue
        mu_coef = row[[c for c in fits.columns if c.startswith("mu_coef_")]].values.astype(float)
        disp_coef = row[[c for c in fits.columns if c.startswith("disp_coef_")]].values.astype(float)
        mu = np.clip(np.exp(Xa @ np.nan_to_num(mu_coef, nan=0.0)), 1e-6, 1e8)
        if not np.all(np.isnan(disp_coef)):
            alpha = np.exp(-Xa @ np.nan_to_num(disp_coef, nan=0.0))
        elif "trend_alpha" in row.index and not pd.isna(row["trend_alpha"]):
            alpha = np.full(len(idx), float(row["trend_alpha"]))
        else:
            alpha = np.full(len(idx), alpha_fn(float(mu.mean())))
        y = data["Y"][idx, data["gene_col"][g]]
        z = marginal_nb_rqr(y, mu, alpha, float(row["tau2"]), seed=seed)
        zdict[g] = z.astype(np.float32)
    return zdict, rows


def _refit_and_score_pool(tr_idx, te_idx, data, genes, scaler, tmp, cache_path=None):
    if cache_path and cache_path.exists():
        with open(cache_path) as f:
            fit = json.load(f)
    else:
        Xs_tr = scaler.transform(data["X_raw"][tr_idx])
        Y_tr = data["Y"][tr_idx][:, [data["gene_col"][g] for g in genes]]
        pd.DataFrame(Xs_tr, columns=config.BIAS_COLUMNS).to_csv(f"{tmp}/Xp.csv.gz")
        pd.DataFrame(Y_tr, columns=genes).to_csv(f"{tmp}/Yp.csv.gz")
        pd.DataFrame({"Batch_ID": data["batch"][tr_idx]}).to_csv(f"{tmp}/batchp.csv.gz")
        pd.DataFrame({"gene": genes}).to_csv(f"{tmp}/genesp.csv", index=False)
        subprocess.run([
            config.RSCRIPT, str(config.GLMM_FIT_POOL_R), "--x", f"{tmp}/Xp.csv.gz", "--y", f"{tmp}/Yp.csv.gz",
            "--batch", f"{tmp}/batchp.csv.gz", "--genes", f"{tmp}/genesp.csv",
            "--rare-overdisp-thr", str(MP["rare_overdisp_thr"]), "--out", f"{tmp}/pool_res.json",
        ], check=True, cwd=str(config.GLMM_FIT_POOL_R.parent))

        with open(f"{tmp}/pool_res.json") as f:
            fit = json.load(f)
        if cache_path:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "w") as f:
                json.dump(fit, f)
    if not fit["ok"]:
        return {}, {}, [dict(gene=g, route="pool", ok=False, fail_reason="pooled_glmm_fit_error") for g in genes]

    beta = np.asarray(fit["beta"])
    tau2 = float(fit["tau2"]) if fit.get("tau2") is not None else 0.0
    alpha_eff = fit["alpha"] if fit["family"] == "negbin" else 1e-8
    mean_hc = dict(zip(fit["gene"], fit["mean_hc"]))
    eps = fit["eps"]
    Xs_te = scaler.transform(data["X_raw"][te_idx])
    Xs_tr = scaler.transform(data["X_raw"][tr_idx])
    mult_te = np.exp(Xs_te @ beta[1:])
    mult_tr = np.exp(Xs_tr @ beta[1:])
    if fit.get("mult_lo") is not None:
        mult_te = np.clip(mult_te, fit["mult_lo"], fit["mult_hi"])
        mult_tr = np.clip(mult_tr, fit["mult_lo"], fit["mult_hi"])

    rqr = _poisson_rqr if fit["family"] == "poisson" else lambda y, mu, seed: marginal_nb_rqr(y, mu, alpha_eff, tau2, seed=seed)
    zdict, zdict_tr, rows = {}, {}, []
    for g in genes:
        mu_te = np.clip((mean_hc[g] + eps) * np.exp(beta[0]) * mult_te, 1e-6, 1e8)
        mu_tr = np.clip((mean_hc[g] + eps) * np.exp(beta[0]) * mult_tr, 1e-6, 1e8)
        y_te = data["Y"][te_idx, data["gene_col"][g]]
        y_tr = data["Y"][tr_idx, data["gene_col"][g]]
        zdict[g] = rqr(y_te, mu_te, 42).astype(np.float32)
        zdict_tr[g] = rqr(y_tr, mu_tr, 1042).astype(np.float32)
        rows.append(dict(gene=g, route="pool", ok=True, family=fit["family"]))
    return zdict, zdict_tr, rows


def run_one_batch(batch_id, data, summary, alpha_fn, disp_prior_path, tmp, limit_genes=None, cache_dir=None):
    """cache_dir (normally the batch's own output directory) persists this
    batch's leave-one-batch-out GLMM fit (model_fits.csv / pool_fit.json) and
    train-fit SHASH params (shash_params.csv) -- a rerun with these cached
    needs neither Rscript nor a SHASH refit, only the cheap RQR rescoring."""
    is_hc, batch = data["is_hc"], data["batch"]
    held = batch == batch_id
    small = np.array([b in data["small_hc_batches"] for b in batch])
    tr_idx = np.where(is_hc & ~held & ~small)[0]
    te_idx = np.where(held)[0]

    genes = summary if limit_genes is None else summary.iloc[:limit_genes]
    model_genes = genes.index[genes["route"] == "model"].tolist()
    pool_genes = genes.index[genes["route"] == "pool"].tolist()
    stage_of = genes["stage"].to_dict()
    scaler = StandardScaler().fit(data["X_raw"][tr_idx])

    t0 = time.perf_counter()
    zdict, zdict_tr, fold_rows = {}, {}, []
    if model_genes:
        model_cache = cache_dir / "model_fits.csv" if cache_dir else None
        fits = _refit_model_genes(tr_idx, data, model_genes, stage_of, scaler, tmp, disp_prior_path, cache_path=model_cache)
        z_m, rows_m = _score_model_genes(fits, model_genes, te_idx, data, scaler, alpha_fn)
        z_m_tr, _ = _score_model_genes(fits, model_genes, tr_idx, data, scaler, alpha_fn, seed=1042, report=False)
        zdict.update(z_m)
        zdict_tr.update(z_m_tr)
        fold_rows += rows_m
    if pool_genes:
        pool_cache = cache_dir / "pool_fit.json" if cache_dir else None
        z_p, z_p_tr, rows_p = _refit_and_score_pool(tr_idx, te_idx, data, pool_genes, scaler, tmp, cache_path=pool_cache)
        zdict.update(z_p)
        zdict_tr.update(z_p_tr)
        fold_rows += rows_p

    # Per-gene SHASH fit on THIS batch's own train-fold (HC minus the held-out
    # batch) in-sample Z, applied to the held-out batch Z -- the LOBO analogue
    # of NormativeModelEngineMixed.fit_shash, never fit on the held-out Z itself.
    shash_cache = cache_dir / "shash_params.csv" if cache_dir else None
    shash_cached = pd.read_csv(shash_cache).set_index("gene") if shash_cache and shash_cache.exists() else None
    genes_scored = [g for g in zdict if g in zdict_tr]
    Z_shash, shash_rows = {}, []
    for g in genes_scored:
        if shash_cached is not None and g in shash_cached.index:
            srow = shash_cached.loc[g]
            sok, sxi, seta, seps, sdelta = bool(srow["ok"]), srow["xi"], srow["eta"], srow["eps"], srow["delta"]
            Z_shash[g] = (shash_transform_to_z(zdict[g], sxi, seta, seps, sdelta) if sok else zdict[g]).astype(np.float32)
        else:
            params, z_corr = fit_and_correct(zdict_tr[g], zdict[g])
            Z_shash[g] = z_corr.astype(np.float32)
            shash_rows.append(dict(gene=g, ok=params["ok"], xi=params["xi"], eta=params["eta"],
                                   eps=params["eps"], delta=params["delta"]))
    if shash_cache and shash_rows:
        shash_cache.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(shash_rows).to_csv(shash_cache, index=False)

    gene_names = list(zdict.keys())
    Z = np.column_stack([zdict[g] for g in gene_names]) if gene_names else np.empty((len(te_idx), 0))
    Z_shash_arr = (np.column_stack([Z_shash.get(g, zdict[g]) for g in gene_names])
                  if gene_names else np.empty((len(te_idx), 0)))
    return dict(
        batch_id=batch_id, n_hc_train=len(tr_idx), n_test=len(te_idx),
        test_names=data["names"][te_idx].tolist(), test_is_hc=is_hc[te_idx].tolist(),
        gene_names=gene_names, Z=Z, Z_shash=Z_shash_arr, fold_info=pd.DataFrame(fold_rows),
        elapsed_s=time.perf_counter() - t0,
    )


def save_batch_result(res, tier, n_dis, out_dir):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", res["batch_id"])
    bdir = out_dir / safe
    bdir.mkdir(parents=True, exist_ok=True)
    np.save(bdir / "Z_test.npy", res["Z"])
    np.save(bdir / "Z_test_shash.npy", res["Z_shash"])
    with open(bdir / "gene_names.pkl", "wb") as f:
        pickle.dump(res["gene_names"], f)
    meta = dict(batch_id=res["batch_id"], tier=tier, n_dis=n_dis,
               n_hc_train=res["n_hc_train"], n_test=res["n_test"],
               test_names=res["test_names"], test_is_hc=res["test_is_hc"],
               n_genes_scored=len(res["gene_names"]), elapsed_s=res["elapsed_s"])
    with open(bdir / "meta.json", "w") as f:
        json.dump(meta, f, indent=2, default=str)
    res["fold_info"].to_csv(bdir / "fold_info.csv", index=False)
    return bdir


def compute_ood(data=None, out_dir=config.LOBO_MIXED_DIR, percentile=95, n_out_thr=2):
    """Inference-only OOD flag for each LOBO_Results_mixed/<batch>/ held-out
    sample -- run_one_batch already refits per batch, so this only needs a
    filter fit per batch (on that batch's OWN train-fold HC, excluding the
    held-out batch, never the held-out batch itself) applied to the
    already-scored test rows. Combines MahalanobisFilter with RangeFilter
    (>=n_out_thr covariates individually outside HC's p1-p99) -- a quadratic
    distance alone dilutes several mildly-extreme covariates across the
    normal ones and misses samples a per-axis count catches (see
    core/ood_filter.py:RangeFilter, tested 2026-07-31 against 4_disease_scoring.ipynb).
    Writes ood_mask.npy (True=inlier) + ood_distance.npy per batch, same
    convention as the v1 engine's compute_lobo_ood.py."""
    data = data or load_full_data()
    name2row = {n: i for i, n in enumerate(data["names"])}
    small = np.array([b in data["small_hc_batches"] for b in data["batch"]])

    for bdir in sorted(Path(out_dir).iterdir()):
        meta_path = bdir / "meta.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        batch_id = meta["batch_id"]
        tr_idx = np.where(data["is_hc"] & (data["batch"] != batch_id) & ~small)[0]
        X_tr = data["X_raw"][tr_idx]
        mahal = MahalanobisFilter(percentile=percentile).fit(X_tr)
        rng_filter = RangeFilter(n_out_thr=n_out_thr).fit(X_tr)

        te_rows = np.array([name2row[n] for n in meta["test_names"]])
        X_te = data["X_raw"][te_rows]
        d = mahal.distances(X_te)
        keep = (d <= mahal.threshold_) & rng_filter.mask(X_te)

        np.save(bdir / "ood_mask.npy", keep)
        np.save(bdir / "ood_distance.npy", d)
        meta["ood_percentile"] = percentile
        meta["ood_threshold"] = mahal.threshold_
        meta["ood_n_out_thr"] = n_out_thr
        meta["n_ood_removed"] = int((~keep).sum())
        meta_path.write_text(json.dumps(meta, indent=2, default=str))
        print(f"{batch_id:40s} n_test={len(keep):4d}  removed_OOD={int((~keep).sum()):3d} ({(~keep).mean()*100:.1f}%)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=str, default=None, help="run a single Batch_ID only (default: all usable batches)")
    ap.add_argument("--limit-genes", type=int, default=None, help="smoke test")
    ap.add_argument("--engine-dir", type=Path, default=config.ENGINE_MIXED_DIR)
    ap.add_argument("--out-dir", type=Path, default=config.LOBO_MIXED_DIR)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    disp_prior_path = args.engine_dir / "disp_prior.json"
    if not disp_prior_path.exists():
        disp_prior_path = None

    print("Loading data...")
    data = load_full_data()
    summary = pd.read_csv(args.engine_dir / "training_summary.csv", index_col="gene")
    summary = summary[summary["ok"] & summary["route"].isin(["model", "pool"])]
    alpha_fn = load_trend(args.engine_dir / "dispersion_trend.json")

    batches_df = usable_batches(data)
    batches_df.to_csv(args.out_dir / "batch_tier_assignment.csv", index=False)
    print(batches_df.to_string())

    todo = batches_df if args.batch is None else batches_df[batches_df["batch"] == args.batch]
    if todo.empty:
        raise SystemExit(f"batch {args.batch!r} not found in usable_batches()")

    tmp = "/tmp/lobo_glmm"
    Path(tmp).mkdir(exist_ok=True)
    for _, brow in todo.iterrows():
        b, tier, n_dis = brow["batch"], brow["tier"], int(brow["n_dis"])
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", b)
        if (args.out_dir / safe / "meta.json").exists() and (args.out_dir / safe / "Z_test_shash.npy").exists():
            print(f"[skip, already done] {b}")
            continue
        print(f"\n=== LOBO batch: {b}  (tier={tier}, n_dis={n_dis}) ===", flush=True)
        res = run_one_batch(b, data, summary, alpha_fn, disp_prior_path, tmp, limit_genes=args.limit_genes,
                            cache_dir=args.out_dir / safe)
        out_dir = save_batch_result(res, tier, n_dis, args.out_dir)
        print(f"  saved -> {out_dir}  ({res['elapsed_s']:.0f}s, "
             f"{len(res['gene_names'])} genes, n_test={res['n_test']})", flush=True)


if __name__ == "__main__":
    main()
