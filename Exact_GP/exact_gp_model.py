import os
import sys
from pathlib import Path
from collections import defaultdict

# config.py lives one level up (project root)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import joblib
import numpy as np
import torch
import torch.nn as nn
import gpytorch
from scipy.sparse import issparse
from scipy.stats import norm
from tqdm.auto import tqdm

from config import PATHS


# ---------------------------------------------------------------------------
# GP model
# ---------------------------------------------------------------------------

class BatchedExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, num_features, batch_shape):
        super().__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean(batch_shape=batch_shape)
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.RBFKernel(ard_num_dims=num_features, batch_shape=batch_shape),
            batch_shape=batch_shape,
        )

    def forward(self, x):
        return gpytorch.distributions.MultivariateNormal(
            self.mean_module(x), self.covar_module(x)
        )


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------

class PearsonResidualScaler:
    def __init__(self, clip_to_sqrt=True):
        self.clip_to_sqrt = clip_to_sqrt
        self.gene_fractions = None
        self.theta_g = None
        self.clip_max = None
        self.clip_min = None

    def fit(self, X_counts):
        Y = X_counts.toarray() if issparse(X_counts) else np.asarray(X_counts)
        n_samples = Y.shape[0]
        lib = Y.sum(axis=1)
        self.gene_fractions = Y.sum(axis=0) / lib.sum()

        mu = lib[:, None] * self.gene_fractions[None, :]
        num = np.sum(mu ** 2, axis=0)
        den = np.sum((Y - mu) ** 2 - mu, axis=0)

        theta = np.where(den > 0, num / den, 100.0)
        self.theta_g = np.clip(theta, 1e-3, 100.0)

        if self.clip_to_sqrt:
            self.clip_max = np.sqrt(n_samples)
            self.clip_min = -np.sqrt(n_samples)
        return self

    def transform(self, X_counts):
        if self.gene_fractions is None:
            raise ValueError("Scaler has not been fitted.")
        Y = X_counts.toarray() if issparse(X_counts) else np.asarray(X_counts)
        lib = Y.sum(axis=1)
        mu = lib[:, None] * self.gene_fractions[None, :]
        var = mu + (mu ** 2 / self.theta_g[None, :])
        res = (Y - mu) / np.sqrt(var + 1e-8)
        if self.clip_to_sqrt:
            res = np.clip(res, self.clip_min, self.clip_max)
        return res

    def fit_transform(self, X_counts):
        return self.fit(X_counts).transform(X_counts)


# ---------------------------------------------------------------------------
# Hurdle pipeline
# ---------------------------------------------------------------------------

class IntegratedHurdlePipeline:

    def __init__(self, gene_names, num_features=10, device=None):
        self.gene_names = np.array(gene_names)
        self.num_features = num_features
        self.device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.gene_registry = {g: {} for g in self.gene_names}
        self.batched_gp_models = {}
        self.lr_model = None
        self.lr_col_map = {}

    @staticmethod
    def _ensure_array(arr, dtype=None):
        if hasattr(arr, "toarray"):
            arr = arr.toarray()
        arr = np.asarray(arr)
        return arr.astype(dtype) if dtype is not None else arr

    def fit_logistic_regression(self, X_train_hc, Y_mask_hc, epochs=500, lr=0.1):
        Y_mask_hc = self._ensure_array(Y_mask_hc, dtype=bool)
        n_samples = len(X_train_hc)
        valid_indices, valid_genes = [], []

        for i, g in enumerate(self.gene_names):
            n_pos = int(Y_mask_hc[:, i].sum())
            info = {"n_pos": n_pos, "total": n_samples}
            if n_pos == n_samples:
                info["status"] = "GP_Only"
            elif n_pos == 0:
                info["status"] = "Never_Expressed"
            else:
                info["status"] = "LR_Only" if n_pos < 50 else "Track_A_and_B"
                valid_indices.append(i)
                valid_genes.append(g)
            self.gene_registry[g].update(info)

        if not valid_indices:
            return

        self.lr_col_map = {g: idx for idx, g in enumerate(valid_genes)}
        X_t = torch.tensor(X_train_hc, dtype=torch.float32, device=self.device)
        Y_bin = torch.tensor(
            Y_mask_hc[:, valid_indices].astype(float), dtype=torch.float32, device=self.device
        )
        self.lr_model = nn.Linear(self.num_features, len(valid_indices)).to(self.device)
        opt = torch.optim.Adam(self.lr_model.parameters(), lr=lr)
        crit = nn.BCEWithLogitsLoss()

        for _ in range(epochs):
            opt.zero_grad()
            loss = crit(self.lr_model(X_t), Y_bin)
            loss.backward()
            opt.step()

        self.lr_model.cpu()
        del X_t, Y_bin
        torch.cuda.empty_cache()

    def fit_GP(self, X_train_hc, Y_mask_hc, Y_target_hc, epochs=50, lr=0.1):
        Y_mask_hc = self._ensure_array(Y_mask_hc, dtype=bool)
        Y_target_hc = self._ensure_array(Y_target_hc, dtype=np.float32)

        buckets = defaultdict(list)
        for g_idx, g in enumerate(self.gene_names):
            if self.gene_registry[g].get("status") in ("Track_A_and_B", "GP_Only"):
                buckets[self.gene_registry[g]["n_pos"]].append((g_idx, g))

        bucket_items = sorted(buckets.items())
        n_buckets = len(bucket_items)
        safe_elems = 150_000_000
        max_cap = 256

        for bkt_idx, (n_pos, gene_tuples) in enumerate(bucket_items):
            max_bs = max(1, min(safe_elems // max(1, n_pos ** 2), max_cap))
            self.batched_gp_models[n_pos] = []

            for c_idx, i in enumerate(range(0, len(gene_tuples), max_bs)):
                chunk = gene_tuples[i : i + max_bs]
                g_indices = [g[0] for g in chunk]
                g_names = [g[1] for g in chunk]
                batch_size = len(g_indices)

                for g_name in g_names:
                    self.gene_registry[g_name]["chunk_idx"] = c_idx

                X_batch = [X_train_hc[Y_mask_hc[:, g_idx].flatten()] for g_idx in g_indices]
                Y_vals = np.stack([
                    Y_target_hc[Y_mask_hc[:, g_idx].flatten(), g_idx] for g_idx in g_indices
                ])
                y_mean = Y_vals.mean(axis=1)
                y_std = Y_vals.std(axis=1) + 1e-6
                Y_scaled = (Y_vals - y_mean[:, None]) / y_std[:, None]

                X_t = torch.tensor(np.stack(X_batch), dtype=torch.float32, device=self.device)
                Y_t = torch.tensor(Y_scaled, dtype=torch.float32, device=self.device)
                b_shape = torch.Size([batch_size])

                likelihood = gpytorch.likelihoods.GaussianLikelihood(batch_shape=b_shape).to(self.device)
                model = BatchedExactGPModel(X_t, Y_t, likelihood, self.num_features, b_shape).to(self.device)
                model.train(); likelihood.train()

                opt = torch.optim.Adam(model.parameters(), lr=lr)
                mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)
                i_loss = f_loss = 0.0

                for epoch in range(epochs):
                    opt.zero_grad()
                    loss = -mll(model(X_t), Y_t).sum()
                    loss.backward()
                    opt.step()
                    if epoch == 0:
                        i_loss = loss.item()
                    if epoch == epochs - 1:
                        f_loss = loss.item()

                model.eval(); likelihood.eval()
                with torch.no_grad(), gpytorch.settings.fast_pred_var():
                    pred = likelihood(model(X_t))
                    Z = (Y_t - pred.mean) / torch.sqrt(pred.variance)
                    z_mean = Z.mean(dim=1).cpu().numpy()
                    z_std = Z.std(dim=1).cpu().numpy()
                    noise_var = np.atleast_1d(likelihood.noise.detach().cpu().numpy().squeeze())

                self.batched_gp_models[n_pos].append({
                    "model_state":      {k: v.cpu().clone() for k, v in model.state_dict().items()},
                    "likelihood_state": {k: v.cpu().clone() for k, v in likelihood.state_dict().items()},
                    "gene_names": g_names,
                    "train_x":    X_t.cpu(),
                    "train_y":    Y_t.cpu(),
                    "y_mean": y_mean, "y_std": y_std,
                    "z_mean": z_mean, "z_std": z_std, "noise_var": noise_var,
                })

                tqdm.write(
                    f"   Bucket [{bkt_idx+1:03d}/{n_buckets:03d}] "
                    f"N_pos:{n_pos:4d} chunk:{c_idx+1}(size:{batch_size:4d}) "
                    f"loss:{i_loss:.1f}→{f_loss:.1f}"
                )
                del X_t, Y_t, loss, pred, model, likelihood, opt, mll
                torch.cuda.empty_cache()

    def predict_disease_chunk(self, n_pos, c_idx, X_disease, Y_mask_chunk, Y_target_chunk):
        Y_mask_chunk = self._ensure_array(Y_mask_chunk, dtype=bool)
        Y_target_chunk = self._ensure_array(Y_target_chunk, dtype=np.float32)

        chunk = self.batched_gp_models[n_pos][c_idx]
        g_names = chunk["gene_names"]
        b_size = len(g_names)
        n_samples = X_disease.shape[0]

        z_state = np.zeros((n_samples, b_size))
        z_quant = np.full((n_samples, b_size), np.nan)

        # Track A: logistic dropout
        self.lr_model.to(self.device)
        self.lr_model.eval()
        X_t = torch.tensor(X_disease, dtype=torch.float32, device=self.device)
        with torch.no_grad():
            logits = self.lr_model(X_t)
            for local_idx, g in enumerate(g_names):
                status = self.gene_registry[g]["status"]
                if status == "GP_Only":
                    continue
                col_idx = self.lr_col_map[g]
                p_det = torch.sigmoid(logits[:, col_idx]).cpu().numpy()
                is_det = Y_mask_chunk[:, local_idx].astype(float)
                denom = np.sqrt(p_det * (1 - p_det) + 1e-6)
                z_state[:, local_idx] = (is_det - p_det) / denom
        self.lr_model.cpu()

        # Track B: GP quantity
        b_shape = torch.Size([b_size])
        train_x = chunk["train_x"].to(self.device)
        train_y = chunk["train_y"].to(self.device)
        likelihood = gpytorch.likelihoods.GaussianLikelihood(batch_shape=b_shape).to(self.device)
        model = BatchedExactGPModel(train_x, train_y, likelihood, self.num_features, b_shape).to(self.device)
        model.load_state_dict(chunk["model_state"])
        likelihood.load_state_dict(chunk["likelihood_state"])
        model.eval(); likelihood.eval()

        mu_list, std_list = [], []
        batch = 1000
        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            for i in range(0, n_samples, batch):
                X_b = torch.tensor(
                    X_disease[i : i + batch], dtype=torch.float32, device=self.device
                ).unsqueeze(0)
                pred = likelihood(model(X_b.expand(b_size, -1, -1)))
                mu_list.append(pred.mean.cpu().numpy())
                std_list.append(torch.sqrt(pred.variance).cpu().numpy())

        mu = np.concatenate(mu_list, axis=1)
        std = np.concatenate(std_list, axis=1)
        y_m, y_s = chunk["y_mean"], chunk["y_std"]

        for local_idx in range(b_size):
            mask = Y_mask_chunk[:, local_idx]
            if mask.any():
                y_scaled = (Y_target_chunk[mask, local_idx] - y_m[local_idx]) / y_s[local_idx]
                safe_std = np.clip(std[local_idx, mask], a_min=1e-6, a_max=None)
                z_quant[mask, local_idx] = (y_scaled - mu[local_idx, mask]) / safe_std

        # Combined Z
        z_total = np.copy(z_state)
        for local_idx in range(b_size):
            mask = Y_mask_chunk[:, local_idx]
            if mask.any():
                zs = z_state[mask, local_idx]
                zq = z_quant[mask, local_idx]
                p_s = 2 * (1 - norm.cdf(np.abs(zs)))
                p_q = 2 * (1 - norm.cdf(np.abs(zq)))
                rho = np.corrcoef(zs, zq)[0, 1]
                rho = np.clip(rho if np.isfinite(rho) else 0.0, -1 + 1e-6, 1 - 1e-6)
                z_total[mask, local_idx] = (zs + zq) / np.sqrt(2 + 2 * rho)

        del model, likelihood, train_x, train_y
        torch.cuda.empty_cache()
        return z_total, z_state, z_quant, g_names


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def save_pipeline(pipeline, dir_path=None, prefix="LogisticGP"):
    dir_path = dir_path or str(PATHS["pipeline_meta"].parent)
    os.makedirs(dir_path, exist_ok=True)

    joblib.dump(
        {"registry": pipeline.gene_registry, "lr_col_map": pipeline.lr_col_map},
        os.path.join(dir_path, f"{prefix}_meta.pkl"),
    )
    if pipeline.lr_model is not None:
        torch.save(pipeline.lr_model.state_dict(), os.path.join(dir_path, f"{prefix}_lr_model.pt"))
    torch.save(pipeline.batched_gp_models, os.path.join(dir_path, f"{prefix}_gp_states.pt"))
    print(f"Pipeline saved to: {dir_path}")


def load_pipeline(gene_names, num_features, device, dir_path=None, prefix="LogisticGP"):
    dir_path = dir_path or str(PATHS["pipeline_meta"].parent)
    pipeline = IntegratedHurdlePipeline(gene_names, num_features, device)

    meta = joblib.load(os.path.join(dir_path, f"{prefix}_meta.pkl"))
    pipeline.gene_registry = meta["registry"]
    pipeline.lr_col_map = meta["lr_col_map"]

    lr_path = os.path.join(dir_path, f"{prefix}_lr_model.pt")
    if os.path.exists(lr_path):
        pipeline.lr_model = nn.Linear(num_features, len(pipeline.lr_col_map)).to(device)
        pipeline.lr_model.load_state_dict(
            torch.load(lr_path, map_location=device, weights_only=True)
        )
        pipeline.lr_model.eval().cpu()

    pipeline.batched_gp_models = torch.load(
        os.path.join(dir_path, f"{prefix}_gp_states.pt"),
        map_location=device, weights_only=False,
    )
    print(f"Pipeline loaded from: {dir_path}")
    return pipeline


# ---------------------------------------------------------------------------
# Inference utilities
# ---------------------------------------------------------------------------

def run_disease_inference(pipeline, valid_genes, X_disease, Y_mask_disease, Y_target_disease):
    import pandas as pd
    n_samples = X_disease.shape[0]
    valid_set = set(valid_genes)

    df_total = pd.DataFrame(index=range(n_samples), columns=valid_genes, dtype=np.float32)
    df_state = pd.DataFrame(index=range(n_samples), columns=valid_genes, dtype=np.float32)
    df_quant = pd.DataFrame(index=range(n_samples), columns=valid_genes, dtype=np.float32)

    for n_pos, chunks in tqdm(pipeline.batched_gp_models.items(), desc="Processing Buckets"):
        for c_idx, chunk in enumerate(chunks):
            g_names = chunk["gene_names"]
            if not any(g in valid_set for g in g_names):
                continue
            col_indices = [np.where(pipeline.gene_names == g)[0][0] for g in g_names]
            z_total, z_state, z_quant, _ = pipeline.predict_disease_chunk(
                n_pos, c_idx, X_disease,
                Y_mask_disease[:, col_indices],
                Y_target_disease[:, col_indices],
            )
            for local_idx, g in enumerate(g_names):
                if g in valid_set:
                    df_total[g] = z_total[:, local_idx]
                    df_state[g] = z_state[:, local_idx]
                    df_quant[g] = z_quant[:, local_idx]

    print("Inference complete.")
    return {"total": df_total, "state": df_state, "quant": df_quant}


def get_gp_qc_metrics(pipeline):
    import pandas as pd
    records = []
    for n_pos, chunks in pipeline.batched_gp_models.items():
        for chunk in chunks:
            for i, g in enumerate(chunk["gene_names"]):
                records.append({
                    "Gene": g, "N_pos": n_pos,
                    "Z_Mean": chunk["z_mean"][i],
                    "Z_Std":  chunk["z_std"][i],
                    "Noise_Var": chunk["noise_var"][i],
                })
    return pd.DataFrame(records)


def filter_valid_genes(df_metrics, z_mean_th=0.2, z_std_min=0.4, z_std_max=1.4, noise_th=0.6):
    m_mean  = df_metrics["Z_Mean"].abs() <= z_mean_th
    m_std_l = df_metrics["Z_Std"] < z_std_min
    m_std_h = df_metrics["Z_Std"] > z_std_max
    m_noise = df_metrics["Noise_Var"] <= noise_th

    valid_mask = m_mean & ~m_std_l & ~m_std_h & m_noise
    df_valid   = df_metrics[valid_mask].copy()
    df_invalid = df_metrics[~valid_mask].copy()

    df_invalid["Failure_Reason"] = "Multiple"
    df_invalid.loc[~m_mean,  "Failure_Reason"] = "Mean_Shift"
    df_invalid.loc[m_std_l,  "Failure_Reason"] = "Low_Std(Overfit)"
    df_invalid.loc[m_std_h,  "Failure_Reason"] = "High_Std"
    df_invalid.loc[~m_noise, "Failure_Reason"] = "High_Noise"

    print(f"Total: {len(df_metrics)} | Valid: {len(df_valid)} | Invalid: {len(df_invalid)}")
    print(df_invalid["Failure_Reason"].value_counts())
    return df_valid, df_invalid
