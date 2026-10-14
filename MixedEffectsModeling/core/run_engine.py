import argparse
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.model_engine_mixed import NormativeModelEngineMixed

TMP_DIR = "/tmp/glmm_train"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--calib-genes", type=int, default=None,
                    help="smoke test: shrink the EB-prior calibration subsample from its default gene count")
    ap.add_argument("--resume", action="store_true",
                    help="resume from an interrupted run's chunked results.csv in "
                         f"{TMP_DIR} instead of starting fresh (only safe against the same code version)")
    ap.add_argument("--out-dir", type=Path, default=config.ENGINE_MIXED_DIR,
                    help="directory to save the trained engine into (defaults to the production engine_state_mixed)")
    ap.add_argument("--calib-only", action="store_true",
                    help="stop after prepare_hyperparams (calib_fits.csv/trend/disp_prior) -- skip the "
                         "full cascade. Use to get a cheap RES for pcis_recalibrate.py before committing "
                         "to the full (hours-long) retrain.")
    ap.add_argument("--pcis-cut", type=float, default=None,
                    help="override config.FIT_PARAMS['pcis_cut'] for this run (e.g. from a prior "
                         "pcis_recalibrate.py pass) instead of the committed default")
    ap.add_argument("--nz-a-max", type=int, default=None,
                    help="override config.NZ_A_MAX for this run (e.g. 0 to route every gene to "
                         "the cascade, disabling the pooled-GLM route)")
    args = ap.parse_args()

    if args.pcis_cut is not None:
        config.FIT_PARAMS["pcis_cut"] = args.pcis_cut
        print(f"pcis_cut override: {args.pcis_cut}")

    if args.nz_a_max is not None:
        config.NZ_A_MAX = args.nz_a_max
        print(f"nz_a_max override: {args.nz_a_max}")

    if not args.resume:
        shutil.rmtree(TMP_DIR, ignore_errors=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    trend_path = args.out_dir / "dispersion_trend.json"
    disp_prior_path = args.out_dir / "disp_prior.json"

    engine = NormativeModelEngineMixed()
    engine.resume = args.resume
    engine.load_hc_data()
    engine.assign_routes()
    n_pool = sum(1 for r in engine.genes.values() if r.route == "pool")
    n_model = sum(1 for r in engine.genes.values() if r.route == "model")
    print(f"HC={engine.X_hc_scaled.shape[0]} genes={len(engine.genes)} nz_a_max={engine.nz_a_max} pool_route={n_pool} model_route={n_model}")

    # Always (re)write X/batch to TMP_DIR, independent of the trend/disp_prior cache check below --
    # prepare_hyperparams short-circuits (skips _write_r_inputs) on a cache hit, and --calib-only
    # stops before train() (which normally writes them). Without this, a cached-hyperparams +
    # --calib-only run leaves TMP_DIR without X.csv.gz/batch.csv.gz, breaking pcis_recalibrate.py's
    # WORK dir requirement.
    engine._write_r_inputs(TMP_DIR)

    ran = engine.prepare_hyperparams(trend_path, disp_prior_path, TMP_DIR, args.calib_genes)
    prior = engine.disp_prior
    print(f"hyperparams {'from calib' if ran else 'cached'}: trend -> {trend_path}, prior -> {disp_prior_path}")
    print(f"disp slope EB prior sd (n_calib={prior['n_calib_genes']})")
    for cov, tau in zip(prior["covariates"], prior["tau_slope"]):
        print(f"  tau_slope[{cov}] = {tau:.4f}")

    if args.calib_only:
        print(f"--calib-only: stopping before the full cascade. RES for pcis_recalibrate.py -> "
              f"{args.out_dir / 'calib_fits.csv'}, WORK -> {TMP_DIR}")
        return

    engine.train(limit=args.limit, tmp_dir=TMP_DIR, disp_prior_path=disp_prior_path)
    engine.save(args.out_dir)

    summary = engine.training_summary()
    print(summary.groupby(["route", "stage"]).size().to_string())
    print(f"ok={int(summary['ok'].sum())}/{len(summary)}")
    ok = summary[summary["ok"]]
    print(f"disp intercept EB: tau_d={engine.disp_tau_d2 ** 0.5:.4f} "
          f"shrink_frac_median={np.median(_shrink_frac(ok)):.3f}")
    print(f"Cook outliers: genes_with_any={int((ok['n_outliers'] > 0).sum())}/{len(ok)} "
          f"mean_removed={ok['n_outliers'].mean():.2f} max={int(ok['n_outliers'].max())} "
          f"refit_failed={int(ok['outlier_refit_failed'].sum())}")
    print(f"Saved -> {args.out_dir}")


def _shrink_frac(ok):
    """Fraction of the way from the raw MLE to the trend that the EB squeeze moved
    log(theta). 0 = kept the gene's own estimate, 1 = pinned to the trend."""
    raw = ok["log_theta_raw"].to_numpy(float)
    eb = ok["log_theta_eb"].to_numpy(float)
    trend = -np.log(ok["trend_alpha"].to_numpy(float))
    gap = raw - trend
    frac = np.where(np.abs(gap) > 1e-9, (eb - raw) / np.where(np.abs(gap) > 1e-9, -gap, 1.0), 0.0)
    return frac[np.isfinite(frac)]


if __name__ == "__main__":
    main()
