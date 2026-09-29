import argparse
import logging
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.calibration import gene_shash_calibration
from MixedEffectsModeling.core.dispersion_trend import load_trend
from MixedEffectsModeling.core.model_engine_mixed import NormativeModelEngineMixed
from MixedEffectsModeling.validation.cv_engine import cv_model_route, cv_pool_route, ppc_metrics

MP = config.SPIKE_PARAMS
THRESHOLDS = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
SHASH_MAX_N = 3000


def naive_calib(z):
    v = np.asarray(z)
    v = v[np.isfinite(v)]
    if len(v) < 8:
        return dict(naive_exceed=np.nan, naive_fdr_reject_rate=np.nan, raw_skew=np.nan, raw_kurtosis=np.nan)
    v = v if len(v) <= SHASH_MAX_N else np.random.default_rng(42).choice(v, SHASH_MAX_N, replace=False)
    calib = gene_shash_calibration(v)
    return dict(naive_exceed=calib["naive_exceed"], naive_fdr_reject_rate=calib["naive_fdr_reject_rate"],
               raw_skew=calib["raw_skew"], raw_kurtosis=calib["raw_kurtosis"])


EMPTY_PPC = dict(y=np.array([]), mu=np.array([]), alpha=np.array([]), tau2=np.array([]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit-genes", type=int, default=None,
                    help="smoke test: only use the first N genes of the nz<max(THRESHOLDS) superset")
    ap.add_argument("--thresholds", type=int, nargs="+", default=None,
                    help="smoke test: override THRESHOLDS with a shorter list")
    args = ap.parse_args()
    thresholds = args.thresholds or THRESHOLDS

    out_dir = config.THRESHOLD_SWEEP_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    for h in list(logging.getLogger().handlers):
        logging.getLogger().removeHandler(h)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(message)s",
        handlers=[logging.FileHandler(out_dir / "sweep.log"), logging.StreamHandler()],
    )
    log = logging.getLogger(__name__)

    summary = pd.read_csv(config.ENGINE_MIXED_DIR / "training_summary.csv", index_col="gene")
    summary = summary[summary["ok"]]
    superset = summary.index[summary["nz"] < max(thresholds)].tolist()[:args.limit_genes]
    nz_of = summary["nz"].to_dict()
    stage_of = summary["stage"].to_dict()
    log.info(f"superset (nz < {max(thresholds)}): {len(superset)} genes")

    e2 = NormativeModelEngineMixed()
    e2.load_hc_data()
    e2.alpha_fn = load_trend()
    n_hc = e2.X_hc_scaled.shape[0]
    folds = list(StratifiedKFold(MP["n_splits"], shuffle=True, random_state=42).split(np.zeros(n_hc), e2.batch))

    tmp = "/tmp/pool_threshold_sweep"
    Path(tmp).mkdir(exist_ok=True)
    ind_cache_path = out_dir / "individual_cv_cache.pkl"
    zdict_ind, ppc_dict_ind, ind_fold_info = {}, {}, {}
    if ind_cache_path.exists():
        with open(ind_cache_path, "rb") as f:
            cached_superset, zdict_ind, ppc_dict_ind, ind_fold_info = pickle.load(f)
        if set(superset) <= set(cached_superset):
            log.info(f"individual (model-route) CV baseline: loaded cache -> {ind_cache_path}")
        else:
            log.info("individual (model-route) CV cache stale for current superset -- recomputing")
            zdict_ind, ppc_dict_ind, ind_fold_info = {}, {}, {}

    if not (set(superset) <= set(ind_fold_info)):
        log.info(f"individual (model-route) CV baseline on {len(superset)} genes...")
        zdict_ind, _, ppc_dict_ind, ind_fold_stat_rows = cv_model_route(e2, superset, stage_of, folds, tmp)
        ind_fold_info = {g: dict(n_folds_ok=0, fail_reasons=[]) for g in superset}
        for r in ind_fold_stat_rows:
            info = ind_fold_info[r["gene"]]
            if r["ok"]:
                info["n_folds_ok"] += 1
            else:
                info["fail_reasons"].append(f"fold{r['fold']}:{r['fail_reason'] or 'unknown'}")
        with open(ind_cache_path, "wb") as f:
            pickle.dump((superset, zdict_ind, ppc_dict_ind, ind_fold_info), f)

    # pooled CV is threshold-specific and cached per threshold, so re-running the sweep
    # (e.g. after adding a new threshold) only computes the missing ones.
    per_t_dfs = []
    for t in thresholds:
        t_path = out_dir / f"pool_vs_individual_naive_z_t{t}.csv"
        if t_path.exists():
            log.info(f"threshold={t}: cached -> {t_path}")
            per_t_dfs.append(pd.read_csv(t_path))
            continue

        genes_t = [g for g in superset if nz_of[g] < t]
        if not genes_t:
            continue
        log.info(f"threshold={t}: pooled CV on {len(genes_t)} genes...")
        zdict_pool, _, ppc_dict_pool, pool_fold_stats = cv_pool_route(e2, genes_t, folds, tmp)
        pool_n_folds_ok = sum(1 for r in pool_fold_stats if r["ok"])
        pool_fail_reasons = ";".join(f"fold{r['fold']}:{r['fail_reason']}" for r in pool_fold_stats if not r["ok"])

        t_rows = []
        for g in genes_t:
            m_ind = naive_calib(zdict_ind.get(g, np.array([])))
            m_pool = naive_calib(zdict_pool.get(g, np.array([])))
            ppc_ind = ppc_metrics(ppc_dict_ind.get(g, EMPTY_PPC))
            ppc_pool = ppc_metrics(ppc_dict_pool.get(g, EMPTY_PPC))
            t_rows.append(dict(
                threshold=t, gene=g, nz=nz_of[g], ind_stage=stage_of[g],
                ind_n_folds_ok=ind_fold_info[g]["n_folds_ok"],
                ind_fail_reasons=";".join(ind_fold_info[g]["fail_reasons"]),
                **{f"ind_{k}": v for k, v in m_ind.items()},
                **{f"ind_{k}": v for k, v in ppc_ind.items()},
                pool_n_folds_ok=pool_n_folds_ok, pool_fail_reasons=pool_fail_reasons,
                **{f"pool_{k}": v for k, v in m_pool.items()},
                **{f"pool_{k}": v for k, v in ppc_pool.items()},
            ))
        t_df = pd.DataFrame(t_rows)
        t_df.to_csv(t_path, index=False)
        log.info(f"Saved -> {t_path} ({len(t_df)} rows)")
        per_t_dfs.append(t_df)

    df = pd.concat(per_t_dfs, ignore_index=True)
    df.to_csv(out_dir / "pool_vs_individual_naive_z.csv", index=False)
    log.info(f"Saved combined -> {out_dir}/pool_vs_individual_naive_z.csv ({len(df)} rows)")


if __name__ == "__main__":
    main()
