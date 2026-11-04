"""Convergence log of the individual GLMM cascade for low-nz genes.

The deployed engine routes genes with nz < NZ_A_MAX to the pooled GLM, so their individual-fit
behaviour is never recorded. This refits every gene with nz in [NZ_LO, NZ_HI) through the same
cascade on all HC, then refits the chosen stage on the 5 CV training folds (same folds as
cv_engine.py), and logs ok / stage / fail reason per gene. Nothing deployed is touched.

Run:  NM_PIPELINE=p0 nohup python MixedEffectsModeling/notebookscripts/EngineDiagnostics/low_nz_convergence.py > MixedEffectsModeling/Notebook_p0/0_EngineDiagnostics/logs/low_nz_convergence.txt 2>&1 &
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.dispersion_trend import load_trend
from MixedEffectsModeling.core.model_engine_mixed import NormativeModelEngineMixed

NZ_LO, NZ_HI = 1, 61


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=config.ENGINE_DIAG_DIR / "low_nz_convergence.csv")
    args = ap.parse_args()
    if args.out.exists():
        print(f"exists: {args.out}")
        return
    trend_path = config.ENGINE_DIR / "dispersion_trend.json"
    prior_path = config.ENGINE_DIR / "disp_prior.json"
    tmp = config.TMP_ROOT / "low_nz_convergence"
    tmp.mkdir(parents=True, exist_ok=True)
    fit_params = tmp / "fit_params.json"
    fit_params.write_text(json.dumps(config.FIT_PARAMS))

    e = NormativeModelEngineMixed()
    e.load_hc_data()
    e.assign_routes()
    e.trend_path = trend_path
    e.alpha_fn = load_trend(trend_path)
    genes = [g for g, r in e.genes.items() if NZ_LO <= r.nz < NZ_HI]
    print(f"HC={len(e.batch)} genes with nz in [{NZ_LO}, {NZ_HI}): {len(genes)}")

    full_out = tmp / "full.csv"
    if not full_out.exists():
        e._write_r_inputs(str(tmp))
        e._write_gene_block(genes, str(tmp), "full")
        e._run_glmm_fit(str(tmp), "full", "cascade", str(full_out), disp_prior_path=prior_path)
    full = pd.read_csv(full_out).set_index("gene")

    folds = list(StratifiedKFold(config.SPIKE_PARAMS["n_splits"], shuffle=True, random_state=42)
                 .split(np.zeros(len(e.batch)), e.batch))
    fold_ok = pd.DataFrame(index=full.index)
    for fi, (tr, _) in enumerate(folds):
        out = tmp / f"fold{fi}.csv"
        if not out.exists():
            X = e.fold_X(tr)
            pd.DataFrame(StandardScaler().fit(X[tr]).transform(X[tr]), columns=config.COVARIATES).to_csv(tmp / f"X_{fi}.csv.gz")
            pd.DataFrame(e.Y_hc[tr][:, [e._gene_col[g] for g in full.index]], columns=full.index).to_csv(tmp / f"Y_{fi}.csv.gz")
            pd.DataFrame({"Batch_ID": e.batch[tr]}).to_csv(tmp / f"batch_{fi}.csv.gz")
            pd.DataFrame({"gene": full.index, "stage": full["stage"].values}).to_csv(tmp / f"genes_{fi}.csv", index=False)
            subprocess.run([config.RSCRIPT, str(config.GLMM_FIT_R), "--x", str(tmp / f"X_{fi}.csv.gz"),
                            "--y", str(tmp / f"Y_{fi}.csv.gz"), "--batch", str(tmp / f"batch_{fi}.csv.gz"),
                            "--genes", str(tmp / f"genes_{fi}.csv"), "--trend", str(trend_path),
                            "--disp-prior", str(prior_path), "--fit-params", str(fit_params),
                            "--mode", "fixed_stage", "--out", str(out)],
                           check=True, cwd=str(config.GLMM_FIT_R.parent))
        fold_ok[f"fold{fi}_ok"] = pd.read_csv(out).set_index("gene")["ok"].astype(bool).reindex(full.index).fillna(False)

    res = pd.DataFrame({
        "nz": [e.genes[g].nz for g in full.index],
        "deployed_route": [e.genes[g].route for g in full.index],
        "full_ok": full["ok"].astype(bool), "full_stage": full["stage"],
        "full_tau2": full["tau2"], "full_n_outliers": full["n_outliers"],
        "nbi_full_eb_reject_reason": full["nbi_full_eb_reject_reason"],
        "nbi_intercept_eb_reject_reason": full["nbi_intercept_eb_reject_reason"],
    }, index=full.index).join(fold_ok)
    res["n_fold_ok"] = fold_ok.sum(axis=1).astype(int)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.out)
    print(f"saved -> {args.out}")
    print(res.groupby(pd.cut(res.nz, [0, 10, 20, 30, 40, 50, 60]), observed=True)
          .agg(n=("full_ok", "size"), full_ok=("full_ok", "mean"),
               all5=("n_fold_ok", lambda s: (s == len(folds)).mean())).round(3).to_string())


if __name__ == "__main__":
    main()
