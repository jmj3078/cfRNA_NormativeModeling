"""Per-sample spike-in detection limit, 3_detection_limit.ipynb design (engine mu as the
injection reference, LOBO half-boot and in-sample KFold arms). Needs runs/<pipe>/{full,lobo}
and full/scores from score_full.py. Writes runs/<pipe>/detection_limit/ and appends
config.DETECTION_LIMIT_CSV (scenario zero_shot, or adaptation with --adapt)."""
import json
import pickle
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _FinalModelingMethodImprovement.config as config
from _FinalModelingMethodImprovement.engine.core.calibration import bh_fdr_reject
from _FinalModelingMethodImprovement.engine.core.marginal_rqr import marginal_nb_rqr
from _FinalModelingMethodImprovement.engine.core.shash import shash_transform_to_z
from _FinalModelingMethodImprovement.engine.validation.lobo_engine import load_full_data, safe_dir
from _FinalModelingMethodImprovement.eval.adaptation import loc
from _FinalModelingMethodImprovement.eval.results import append

N_PERTURB_GENES = 500
LOG2FCS = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
N_BOOTSTRAP = 5
N_FOLDS = 10
Q_LEVELS = [0.05, 0.10, 0.20, 0.25]
SPLITS = ["lobo", "insample_cv"]
OUT_DIR = config.RUN_DIR / "detection_limit"
SCORES_DIR = config.ENGINE_MIXED_DIR / "scores"

data = load_full_data()
name2row = {n: i for i, n in enumerate(data["names"])}


def mu_alpha_tau2(b, genes, sample_rows):
    Xs = b["scaler"].transform(data["X_raw"][sample_rows])
    Xa = np.column_stack([np.ones(len(sample_rows)), Xs])
    mu, alpha, tau2 = {}, {}, {}
    fits = b["fits"]
    mu_cols = [c for c in fits.columns if c.startswith("mu_coef_")]
    disp_cols = [c for c in fits.columns if c.startswith("disp_coef_")]
    for g in genes:
        row = fits.loc[g]
        mu_coef = row[mu_cols].values.astype(float)
        disp_coef = row[disp_cols].values.astype(float)
        mu[g] = np.clip(np.exp(Xa @ np.nan_to_num(mu_coef, nan=0.0)), 1e-6, 1e8)
        alpha[g] = (np.exp(-Xa @ np.nan_to_num(disp_coef, nan=0.0)) if not np.all(np.isnan(disp_coef))
                    else np.full(len(sample_rows), float(row["trend_alpha"])))
        tau2[g] = float(row["tau2"])
    return mu, alpha, tau2


def perturb_genes_fold(b, genes, sample_rows, log2fc, seed):
    mu, alpha, tau2 = mu_alpha_tau2(b, genes, sample_rows)
    z_of = {}
    for g in genes:
        y_real = data["Y"][sample_rows, data["gene_col"][g]]
        y_pert = np.round(np.maximum(y_real + mu[g] * (2 ** log2fc - 1), 0))
        z_raw = marginal_nb_rqr(y_pert, mu[g], alpha[g], tau2[g], seed=seed + zlib.crc32(g.encode()) % 9973)
        srow = b["shash_p"].loc[g]
        z = shash_transform_to_z(z_raw, srow.xi, srow.eta, srow.eps, srow.delta) if srow.ok else z_raw
        z_of[g] = np.clip(z, -50, 50)
    return z_of


def load_batch(batch_id):
    bdir = config.LOBO_MIXED_DIR / safe_dir(batch_id)
    meta = json.loads((bdir / "meta.json").read_text())
    gene_names = pickle.load(open(bdir / "gene_names.pkl", "rb"))
    Z_test = np.load(bdir / "Z_test_shash.npy")
    test_is_hc = np.array(meta["test_is_hc"])
    fits = pd.read_csv(bdir / "model_fits.csv").set_index("gene")
    fits = fits[fits["ok"]]
    shash_p = pd.read_csv(bdir / "shash_params.csv").set_index("gene")
    universe = [g for g in fits.index if g in shash_p.index and g in gene_names]
    is_hc, batch, small = data["is_hc"], data["batch"], data["small_hc_batches"]
    tr_idx = np.where(is_hc & (batch != batch_id) & ~np.isin(batch, list(small)))[0]
    assert sorted(data["names"][tr_idx]) == sorted(json.loads(config.LOBO_FOLDS_PATH.read_text())[batch_id]["train"])
    hc_test_rows = np.array([name2row[n] for n in np.array(meta["test_names"])[test_is_hc]])
    return dict(gene_names=gene_names, gene_pos={g: j for j, g in enumerate(gene_names)},
                Z_hc=Z_test[test_is_hc], universe=universe, fits=fits, shash_p=shash_p,
                scaler=StandardScaler().fit(data["X_raw"][tr_idx]), hc_test_rows=hc_test_rows)


def load_insample():
    hc_meta = pd.read_csv(SCORES_DIR / "hc_meta.csv")
    Z_hc = np.load(SCORES_DIR / "Z_hc_shash.npy")
    gene_names = pickle.load(open(SCORES_DIR / "gene_names.pkl", "rb"))
    genes = pickle.load(open(config.ENGINE_MIXED_DIR / "genes.pkl", "rb"))
    summary = pd.read_csv(config.ENGINE_MIXED_DIR / "training_summary.csv").set_index("gene")
    scaler = pickle.load(open(config.ENGINE_MIXED_DIR / "scaler.pkl", "rb"))
    universe = [g for g in summary.index[summary["ok"].fillna(False)] if g in genes and genes[g].mu_coef is not None]
    recs = [genes[g] for g in universe]
    mu_coef_mat = np.nan_to_num(np.stack([r.mu_coef for r in recs]).astype(float), nan=0.0)
    disp_coef_mat = np.stack([r.disp_coef for r in recs]).astype(float)
    has_disp = ~np.all(np.isnan(disp_coef_mat), axis=1)
    disp_coef_mat = np.nan_to_num(disp_coef_mat, nan=0.0)
    su = summary.loc[universe]
    ok_arr = su["cv_shash_ok"].fillna(False).infer_objects(copy=False).values.astype(bool)
    fits = pd.DataFrame(mu_coef_mat, columns=[f"mu_coef_{i}" for i in range(mu_coef_mat.shape[1])], index=universe)
    for i in range(disp_coef_mat.shape[1]):
        fits[f"disp_coef_{i}"] = np.where(has_disp, disp_coef_mat[:, i], np.nan)
    fits["trend_alpha"] = [r.trend_alpha for r in recs]
    fits["tau2"] = [r.tau2 for r in recs]
    shash_p = pd.DataFrame(dict(
        xi=np.where(ok_arr, su["cv_shash_xi"].fillna(0.0).values, 0.0),
        eta=np.where(ok_arr, su["cv_shash_eta"].fillna(1.0).values, 1.0),
        eps=np.where(ok_arr, su["cv_shash_eps"].fillna(0.0).values, 0.0),
        delta=np.where(ok_arr, su["cv_shash_delta"].fillna(1.0).values, 1.0), ok=ok_arr), index=universe)
    return dict(gene_names=gene_names, gene_pos={g: j for j, g in enumerate(gene_names)}, Z_hc=Z_hc,
                universe=universe, fits=fits, shash_p=shash_p, scaler=scaler,
                hc_test_rows=np.array([name2row[s] for s in hc_meta["sample"]]))


def iter_lobo():
    for batch_id in sorted(json.loads(config.LOBO_FOLDS_PATH.read_text())):
        b = load_batch(batch_id)
        n_hc = len(b["hc_test_rows"])
        for boot in range(N_BOOTSTRAP):
            perm = np.random.default_rng(1000 * boot + zlib.crc32(batch_id.encode()) % 9973).permutation(n_hc)
            idx_A = perm[: n_hc // 2]
            yield batch_id, b, b["hc_test_rows"][idx_A], idx_A, boot, perm[n_hc // 2:]


def iter_insample():
    b = load_insample()
    frac = (np.abs(b["Z_hc"]) > 1.96).mean(axis=1)
    q1, q3 = np.percentile(frac, [25, 75])
    keep = np.where(~(frac > q3 + 1.5 * (q3 - q1)))[0]
    for fold, (_, idx_k) in enumerate(KFold(N_FOLDS, shuffle=True, random_state=0).split(np.arange(len(keep)))):
        idx_A = keep[idx_k]
        yield f"insample_fold{fold}", b, b["hc_test_rows"][idx_A], idx_A, fold, None


def run(split, adapt=False):
    rows = []
    for group_id, b, rows_A, idx_A, boot, idx_B in {"lobo": iter_lobo, "insample_cv": iter_insample}[split]():
        for log2fc in LOG2FCS:
            rng = np.random.default_rng(int(1e6 * log2fc) + boot + zlib.crc32(group_id.encode()) % 9973)
            pert_genes = list(rng.choice(b["universe"], min(N_PERTURB_GENES, len(b["universe"])), replace=False))
            z_of = perturb_genes_fold(b, pert_genes, rows_A, log2fc, seed=int(1e6 * log2fc) + boot)
            Z_A = b["Z_hc"][idx_A].copy()
            labels = np.zeros(len(b["gene_names"]), dtype=bool)
            for g in pert_genes:
                Z_A[:, b["gene_pos"][g]] = z_of[g]
                labels[b["gene_pos"][g]] = True
            if adapt:
                Z_A = Z_A - loc(b["Z_hc"][idx_B])
            finite = np.isfinite(Z_A)
            p = np.where(finite, 2 * norm.sf(np.abs(np.where(finite, Z_A, 0.0))), np.nan)
            for q in Q_LEVELS:
                tp = fp = fn = tn = 0
                for i in range(Z_A.shape[0]):
                    rej = np.zeros(Z_A.shape[1], dtype=bool)
                    rej[finite[i]] = bh_fdr_reject(p[i, finite[i]], q=q)
                    tp += int((rej & labels).sum())
                    fp += int((rej & ~labels).sum())
                    fn += int((~rej & labels).sum())
                    tn += int((~rej & ~labels).sum())
                rows.append(dict(batch=group_id, boot=boot, log2fc=log2fc, q=q, n_samples=int(Z_A.shape[0]),
                                 n_universe=int(Z_A.shape[1]), n_perturbed=len(pert_genes), tp=tp, fp=fp, fn=fn, tn=tn))
        print(split, group_id, boot, "done", flush=True)
    return pd.DataFrame(rows)


def main():
    """--adapt: LOBO only; the perturbed half is shifted by the other (unperturbed) half's
    per-gene mean Z (eval/adaptation.py), i.e. the half-boot split doubles as the
    adaptation split."""
    adapt = "--adapt" in sys.argv
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    parts = []
    for split in (["lobo"] if adapt else SPLITS):
        df = run(split, adapt)
        df.to_csv(OUT_DIR / f"sweep_persample_fold__{split}{'_adapt' if adapt else ''}.csv", index=False)
        parts.append(df.assign(split=split, method=config.PIPELINE,
                               scenario="adaptation" if adapt else "zero_shot"))
    append(config.DETECTION_LIMIT_CSV, pd.concat(parts, ignore_index=True))


if __name__ == "__main__":
    main()
