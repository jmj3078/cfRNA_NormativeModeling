import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import _FinalModelingMethodImprovement.config as config
from viz_style import apply_style

ENGINE_MIXED_DIR = config.ENGINE_MIXED_DIR
PCIS_CAL_DIR = config.PCIS_CAL_DIR
PCIS_CAL_FIG_DIR = config.PCIS_CAL_FIG_DIR
TARGET_RATES = [1e-3, 5e-4, 3e-4, 2e-4, 1.5e-4, 1e-4, 5e-5, 1e-5]
TOPK_PER_GENE = 50


def load_null(path):
    d = pd.read_csv(path)
    return d[np.isfinite(d.pcis)].copy()


def rate_table(d, real_summary=None):
    n = int(d.n_obs.iloc[0])
    n_genes = int(d.gene.nunique())
    n_total = n_genes * n
    p = np.sort(d.pcis.values)[::-1]
    rows = []
    for rate in TARGET_RATES:
        k = int(round(rate * n_total))
        if not (1 <= k <= len(p)):
            continue
        c = float(p[k - 1])
        ex = d.pcis > c
        per_gene = ex.groupby(d.gene).sum().reindex(d.gene.unique(), fill_value=0)
        rows.append({
            "target_per_obs_rate": rate,
            "population_percentile": 1 - rate,
            "cut": c,
            "null_removed_per_gene": float(per_gene.mean()),
            "max_removed_per_gene": int(per_gene.max()),
            "topk_saturation_frac": per_gene.max() / TOPK_PER_GENE,
        })
    t = pd.DataFrame(rows)
    if real_summary is not None:
        r = real_summary[real_summary.ok.astype(str).str.upper() == "TRUE"]
        t["real_removed_per_gene"] = float(r.n_outliers.mean())
        t["null_share_of_real_removals"] = t.null_removed_per_gene / t["real_removed_per_gene"]
    return t


def figure(d, outdir):
    apply_style()
    os.makedirs(outdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(np.log10(d.pcis), bins=150, color="0.6")
    for c in [0.5, 1.0, 2.0]:
        ax.axvline(np.log10(c), ls="--", label=f"cut={c:g}")
    ax.set_xlabel("PCIS (null simulation)"); ax.set_ylabel("count"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "null_pcis_distribution.png"), dpi=200); plt.close(fig)


def run_all(null_path, outdir=None, fig_dir=None, real_summary_path=None):
    outdir = str(PCIS_CAL_DIR if outdir is None else outdir)
    fig_dir = str(PCIS_CAL_FIG_DIR if fig_dir is None else fig_dir)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    d = load_null(null_path)
    real = None
    rsp = str(ENGINE_MIXED_DIR / "training_summary.csv") if real_summary_path is None else str(real_summary_path)
    if os.path.isfile(rsp):
        real = pd.read_csv(rsp)

    t = rate_table(d, real)
    t.to_csv(os.path.join(outdir, "pcis_rate_table.csv"), index=False)
    figure(d, fig_dir)
    return t


if __name__ == "__main__":
    t = run_all(sys.argv[1] if len(sys.argv) > 1 else "/tmp/null_pcis.csv")
    print(t.to_string(index=False))
