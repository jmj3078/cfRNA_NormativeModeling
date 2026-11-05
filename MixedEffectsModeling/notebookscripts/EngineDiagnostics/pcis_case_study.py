"""PCIS case study with the deployed stage, prior, trend and cut: re-fit every deployed gene that had PCIS removals, with and without the removal,
   scoring the same observations by PCIS (trend dispersion), PCIS with the gene's own dispersion,
   and classical Cook's distance (own dispersion, fixed-effect leverage).

Run:  NM_PIPELINE=p0 nohup python MixedEffectsModeling/notebookscripts/EngineDiagnostics/pcis_case_study.py > MixedEffectsModeling/Notebook_p0/0_EngineDiagnostics/logs/pcis_case_study.txt 2>&1 &
"""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.model_engine_mixed import NormativeModelEngineMixed

R_SCRIPT = Path(__file__).resolve().parent / "pcis_case_study.R"
OUT_GENE = config.ENGINE_DIAG_DIR / "pcis_cases_gene.csv"
OUT_OBS = config.ENGINE_DIAG_DIR / "pcis_cases_obs.csv.gz"


def main():
    if OUT_GENE.exists():
        print(f"exists: {OUT_GENE}")
        return
    summ = pd.read_csv(config.ENGINE_DIR / "training_summary.csv", index_col="gene")
    cases = summ[(summ.route == "model") & summ.ok.astype(bool) & (summ.n_outliers > 0)]
    e = NormativeModelEngineMixed()
    e.load_hc_data()
    tmp = config.TMP_ROOT / "pcis_case_study"
    e._write_r_inputs(str(tmp))
    fp = tmp / "fit_params.json"
    fp.write_text(json.dumps(config.FIT_PARAMS))
    prior = str(config.ENGINE_DIR / "disp_prior.json")
    OUT_GENE.parent.mkdir(parents=True, exist_ok=True)
    if not OUT_GENE.exists():
        case_study(e, cases, tmp, fp, prior)


def case_study(e, cases, tmp, fp, prior):
    genes = cases.index.tolist()
    pd.DataFrame(e.Y_hc[:, [e._gene_col[g] for g in genes]], columns=genes).to_csv(tmp / "Y_cases.csv.gz")
    sym = sc.read_h5ad(config.H5AD_PATH, backed="r").var["GeneName"]
    pd.DataFrame({"gene": genes, "symbol": sym.reindex(genes).values, "stage": cases.stage.values,
                  "trend_alpha": cases.trend_alpha.values, "n_outliers_deployed": cases.n_outliers.values,
                  "full_obs": (cases.n_outliers >= 2).values}).to_csv(tmp / "genes_cases.csv", index=False)
    pd.DataFrame({"sample": e.names, "batch": e.batch}).to_csv(tmp / "samples.csv", index=False)
    print(f"HC={len(e.batch)} genes with removals: {len(genes)} (>=2 removals: {int((cases.n_outliers >= 2).sum())})")
    subprocess.run([config.RSCRIPT, str(R_SCRIPT), str(tmp), prior,
                    str(fp), str(OUT_GENE), str(OUT_OBS)], check=True)
    g = pd.read_csv(OUT_GENE).merge(pd.read_csv(tmp / "genes_cases.csv"), on=["gene", "stage", "trend_alpha"])
    g.to_csv(OUT_GENE, index=False)
    pd.read_csv(tmp / "samples.csv").to_csv(config.ENGINE_DIAG_DIR / "pcis_cases_samples.csv", index=False)
    print(f"saved -> {OUT_GENE}, {OUT_OBS}")
    print(f"removal count reproduced: {(g.n_drop == g.n_outliers_deployed).mean():.3f}")


if __name__ == "__main__":
    main()
