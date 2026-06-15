"""Model definitions and Z-score scoring utilities for NB-GP normative modeling.

Extracted from ``nb_gaussian_processing_gpytorch_application.ipynb`` so that
standalone training/scoring scripts can reuse the exact same model code
without depending on the notebook.
"""

import warnings

import numpy as np
import torch
import torch.nn as nn
import statsmodels.api as sm
from scipy.stats import nbinom, norm


# =====================================================================
# Shared helpers
# =====================================================================
def _estimate_theta_mom(y_train):
    mean_y = y_train.float().mean()
    var_y = y_train.float().var(unbiased=True)
    if var_y <= mean_y + 1e-4:
        return 100.0

    theta_mom = (mean_y ** 2) / (var_y - mean_y)
    return torch.clamp(theta_mom, min=1e-4, max=1e3).item()


# =====================================================================
# GLM: statsmodels NegativeBinomialP(p=1)
# =====================================================================
class NBGLM:
    """NB GLM wrapper. Falls back to intercept-only model on convergence failure."""

    def __init__(self):
        self._result = None
        self._mu_fb = None  # fallback mean
        self.theta_ = None

    def fit(self, X_train: np.ndarray, y_train: np.ndarray):
        X_c = sm.add_constant(X_train, has_constant="add")
        try:
            with warnings.catch_warnings(), np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                warnings.simplefilter("ignore")
                res = sm.NegativeBinomial(y_train, X_c).fit(disp=False, maxiter=200)
            self._result = res

            alpha = float(res.params[-1])
            if alpha > 1e-6:
                self.theta_ = 1.0 / alpha
            else:
                self.theta_ = 1e4

        except Exception:
            self._result = None
            self._mu_fb = float(y_train.mean()) if len(y_train) > 0 else 1.0
            self.theta_ = None
        return self

    def predict(self, X_test: np.ndarray):
        """Returns (mu_pred: ndarray, theta: float)."""
        X_c = sm.add_constant(X_test, has_constant="add")
        if self._result is not None:
            mu = self._result.predict(X_c).astype(np.float32)
        else:
            mu = np.full(len(X_test), self._mu_fb, dtype=np.float32)
        return np.clip(mu, 1e-4, 1e6), self.theta_


# =====================================================================
# Bayesian NB-GLM (Laplace approximation over linear coefficients)
# =====================================================================
class BayesianNBGLM(nn.Module):
    def __init__(self, n_features, theta_init=2.0, log_sigma_init=0.0):
        super().__init__()
        self.raw_theta = nn.Parameter(torch.tensor(float(np.log(np.exp(theta_init) - 1.0 + 1e-6))))
        self.log_sigma = nn.Parameter(torch.tensor(float(log_sigma_init)))

    @property
    def theta(self):
        return nn.functional.softplus(self.raw_theta) + 0.1

    @property
    def sigma(self):
        return self.log_sigma.exp()

    @staticmethod
    def _aug(X):
        """Prepend intercept column."""
        return torch.cat([torch.ones(len(X), 1, device=X.device, dtype=X.dtype), X], dim=1)

    def find_map(self, X, y, beta0=None, max_iter=100, tol=1e-7):
        """Newton MAP for beta; also returns neg_H = -∇²log p(beta|y) (positive definite)."""
        Xa = self._aug(X)
        n, p = Xa.shape
        dev = Xa.device
        theta = float(self.theta.detach())
        sigma2 = float(self.sigma.detach().pow(2))

        beta = beta0.detach().clone() if beta0 is not None else torch.zeros(p, device=dev)

        with torch.no_grad():
            for _ in range(max_iter):
                mu = (Xa @ beta).exp().clamp(1e-4, 1e6)
                grad_ll = y - mu * (y + theta) / (mu + theta)
                W = (mu * theta * (y + theta) / (mu + theta).pow(2)).clamp(1e-8)
                grad_pos = Xa.T @ grad_ll - beta / sigma2
                neg_H = (Xa.T * W) @ Xa + torch.eye(p, device=dev) / sigma2
                try:
                    L = torch.linalg.cholesky(neg_H)
                    d = torch.cholesky_solve(grad_pos[:, None], L).squeeze()
                except Exception:
                    break
                beta_new = beta + d
                if (beta_new - beta).norm().item() < tol:
                    beta = beta_new
                    break
                beta = beta_new

            # recompute neg_H at final beta for accurate posterior covariance
            mu = (Xa @ beta).exp().clamp(1e-4, 1e6)
            W = (mu * theta * (y + theta) / (mu + theta).pow(2)).clamp(1e-8)
            neg_H = (Xa.T * W) @ Xa + torch.eye(p, device=dev) / sigma2

        return beta, neg_H

    def log_marginal_lik(self, X, y, beta_map, neg_H):
        """Laplace log marginal likelihood for optimising theta and sigma."""
        Xa = self._aug(X)
        p = Xa.shape[1]
        dev = Xa.device
        theta = self.theta
        sigma2 = self.sigma.pow(2)

        b = beta_map.detach()
        mu = (Xa @ b).exp().clamp(1e-4, 1e6)

        log_lik = (
            torch.lgamma(y + theta) - torch.lgamma(theta) - torch.lgamma(y + 1)
            + theta * (theta.log() - (theta + mu).log())
            + y * (mu.log() - (theta + mu).log())
        ).sum()

        log_prior = -0.5 * (b.pow(2).sum() / sigma2 + p * sigma2.log())

        try:
            L = torch.linalg.cholesky(neg_H.detach() + 1e-6 * torch.eye(p, device=dev))
            log_det = L.diagonal().log().sum()
        except Exception:
            log_det = torch.tensor(0.0, device=dev)

        return log_lik + log_prior - log_det

    def posterior_predictive(self, beta_map, neg_H, X_test):
        """Returns (f_mean, f_var) under q(beta) = N(beta_MAP, neg_H^{-1})."""
        Xa = self._aug(X_test)
        p = beta_map.shape[0]
        dev = Xa.device
        f_mean = Xa @ beta_map.detach()
        try:
            L = torch.linalg.cholesky(neg_H.detach() + 1e-6 * torch.eye(p, device=dev))
            V = torch.linalg.solve_triangular(L, Xa.T, upper=False)  # (p, N*)
            f_var = V.pow(2).sum(0).clamp(0)
        except Exception:
            f_var = torch.zeros(len(X_test), device=dev)
        return f_mean, f_var


def train_bayesian_nbglm(
    X_train, y_train,
    max_epochs=500, lr=0.01,
    patience=15, abs_tol=1e-4,
):
    dev = X_train.device

    # Align theta initialization with Laplace-GP using MoM exclusively on training data
    theta_init_val = _estimate_theta_mom(y_train)
    model = BayesianNBGLM(X_train.shape[-1], theta_init=theta_init_val).to(dev)

    opt = torch.optim.Adam(model.parameters(), lr=lr)
    beta_prev = None

    best_loss, no_improve = float("inf"), 0
    for _ in range(max_epochs):
        model.train()
        opt.zero_grad()

        with torch.no_grad():
            beta_map, neg_H = model.find_map(X_train, y_train, beta0=beta_prev)
        beta_prev = beta_map.detach()

        try:
            loss = -model.log_marginal_lik(X_train, y_train, beta_map, neg_H)
        except Exception:
            continue

        if torch.isnan(loss):
            continue

        loss.backward()
        opt.step()

        v = loss.item()
        if best_loss == float("inf"):
            best_loss = v
            continue

        # Absolute improvement
        if (best_loss - v) > abs_tol:
            best_loss, no_improve = v, 0
        else:
            no_improve += 1
            if no_improve >= patience:
                break

    model.eval()
    with torch.no_grad():
        beta_map, neg_H = model.find_map(X_train, y_train, beta0=beta_prev)

    return model, beta_map.detach(), neg_H.detach()


# =====================================================================
# Laplace Approximation NB-GP  (IFT gradient)
# =====================================================================
def _nb_grad_hess(f, y, theta):
    """NB log-likelihood gradient and diagonal Hessian magnitude W at f = log mu."""
    mu = f.exp().clamp(1e-4, 1e6)
    grad = y - mu * (y + theta) / (mu + theta)
    W = (mu * theta * (y + theta) / (mu + theta).pow(2)).clamp(1e-8)
    return grad, W


class ImplicitLaplaceMAP(torch.autograd.Function):
    @staticmethod
    def forward(ctx, K, y, theta, f0, max_iter=50, tol=1e-5):
        """
        Executes Newton-Raphson decoupled from the computation graph.
        Guarantees O(1) space complexity over iterations and exact MAP convergence.
        """
        with torch.no_grad():
            n = len(y)
            dev = K.device
            f = (f0.clone() if f0 is not None else torch.zeros(n, device=dev))
            for _ in range(max_iter):
                grad_ll, W = _nb_grad_hess(f, y, theta)
                sqrt_W = W.sqrt()

                B = torch.eye(n, device=dev) + sqrt_W[:, None] * K * sqrt_W[None, :]
                try:
                    L = torch.linalg.cholesky(B)
                except RuntimeError:
                    # Fallback to current state if unstable during internal search
                    break

                b = W * f + grad_ll
                inner = torch.cholesky_solve((sqrt_W * (K @ b))[:, None], L).squeeze()
                f_new_calc = K @ (b - sqrt_W * inner)
                # Damping factor of 0.5 for stability; can be tuned or made adaptive if desired
                f_new = f + 0.5 * (f_new_calc - f)

                if torch.max(torch.abs(f_new - f)) < tol:
                    f = f_new
                    break
                f = f_new
        ctx.save_for_backward(K, y, theta, f)
        return f

    @staticmethod
    def backward(ctx, grad_f_map):
        """
        Analytically computes gradients w.r.t K and theta using the fixed-point
        constraints and Woodbury matrix identity.
        """
        K, y, theta, f = ctx.saved_tensors

        # Recompute local gradient and Hessian at MAP with localized autograd
        with torch.enable_grad():
            theta_b = theta.detach().requires_grad_(True)
            mu = f.exp().clamp(1e-4, 1e6)
            grad_ll = y - mu * (y + theta_b) / (mu + theta_b)
            W = (mu * theta_b.detach() * (y + theta_b.detach()) / (mu + theta_b.detach()).pow(2)).clamp(1e-8)

        sqrt_W = W.sqrt()
        B = torch.eye(len(y), device=K.device) + sqrt_W[:, None] * K * sqrt_W[None, :]

        try:
            L = torch.linalg.cholesky(B)
        except RuntimeError:
            jitter = 1e-5 * B.diagonal().mean().clamp(min=1e-6)
            L = torch.linalg.cholesky(B + jitter * torch.eye(len(y), device=K.device))

        v = grad_f_map
        Kv = K @ v
        inner_v = torch.cholesky_solve((sqrt_W * Kv)[:, None], L).squeeze()

        # w = K^{-1} u  (where u = H^{-1} v, H = K^{-1} + W)
        w = v - sqrt_W * inner_v
        u = K @ w

        # Analytical gradients via surrogate objective, ensuring symmetry for K
        grad_K = 0.5 * (torch.outer(w, grad_ll.detach()) + torch.outer(grad_ll.detach(), w))
        grad_theta, = torch.autograd.grad(grad_ll, theta_b, grad_outputs=u)

        return grad_K, None, grad_theta, None, None, None


class LaplaceNBGP(nn.Module):
    def __init__(self, n_features, theta_init=2.0):
        super().__init__()
        self.log_output_scale = nn.Parameter(torch.zeros(1))
        self.log_length_scale = nn.Parameter(torch.zeros(n_features))
        self.raw_theta = nn.Parameter(
            torch.tensor(float(np.log(np.exp(theta_init) - 1.0 + 1e-6)))
        )

    @property
    def theta(self):
        return nn.functional.softplus(self.raw_theta) + 0.1

    def kernel(self, X1, X2=None):
        if X2 is None:
            X2 = X1
        ls = self.log_length_scale.exp() + 1e-8
        os = self.log_output_scale.exp() + 1e-8
        diff = (X1.unsqueeze(1) - X2.unsqueeze(0)) / ls
        return os * (-0.5 * diff.pow(2).sum(-1)).exp()

    def _safe_cholesky(self, A, label=None):
        """Cholesky with relative jitter escalation.

        Jitter is scaled by A's mean diagonal magnitude (rather than fixed
        absolute values) so it adapts to the kernel's output-scale. If
        escalation beyond rel_jit=1e-3 is needed, prints a one-line warning
        (when `label` is given) so silently-inflated posterior variances
        are visible.
        """
        dev = A.device
        n = len(A)
        scale = A.diagonal().abs().mean().clamp(min=1e-6)
        for rel_jit in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
            jitter = rel_jit * scale
            try:
                L = torch.linalg.cholesky(A + jitter * torch.eye(n, device=dev))
                if rel_jit >= 1e-3 and label is not None:
                    print(f"  [jitter] {label}: rel_jitter={rel_jit:.0e} (scale={scale.item():.3g})")
                return L
            except Exception:
                continue
        raise RuntimeError("Cholesky failed after jitter escalation")

    def find_map(self, K, y, f0=None, max_iter=15):
        return ImplicitLaplaceMAP.apply(K, y, self.theta, f0, max_iter)

    def log_marginal_lik(self, K, f_map, y):
        """Laplace log marginal likelihood. f_map is kept live — gradient flows via IFT."""
        theta = self.theta
        f = f_map
        dev = K.device

        mu = f.exp().clamp(1e-4, 1e6)
        W = (mu * theta * (y + theta) / (mu + theta).pow(2)).clamp(1e-8)
        sqrt_W = W.sqrt()

        B = torch.eye(len(y), device=dev) + sqrt_W[:, None] * K * sqrt_W[None, :]
        try:
            L = self._safe_cholesky(B, label="log_marginal_lik:B")
        except RuntimeError:
            return K.sum() * 0.0 - 1e8

        log_lik = (
            torch.lgamma(y + theta) - torch.lgamma(theta) - torch.lgamma(y + 1)
            + theta * (theta.log() - (theta + mu).log())
            + y * (mu.log() - (theta + mu).log())
        ).sum()

        K_chol = self._safe_cholesky(K, label="log_marginal_lik:K")
        alpha = torch.cholesky_solve(f[:, None], K_chol).squeeze()
        log_prior = -0.5 * (f * alpha).sum()
        log_det = -L.diagonal().log().sum()
        return log_lik + log_prior + log_det

    def posterior_predictive(self, K_train, f_map, y, K_cross, k_test_diag):
        f = f_map.detach()
        dev = K_train.device
        theta = self.theta.detach()

        _, W = _nb_grad_hess(f, y, theta)
        sqrt_W = W.sqrt()
        B = torch.eye(len(y), device=dev) + sqrt_W[:, None] * K_train * sqrt_W[None, :]
        try:
            L = self._safe_cholesky(B, label="posterior_predictive:B")
        except RuntimeError:
            return torch.zeros(k_test_diag.shape, device=dev), k_test_diag

        K_chol = self._safe_cholesky(K_train, label="posterior_predictive:K")
        alpha = torch.cholesky_solve(f[:, None], K_chol).squeeze()
        f_mean = K_cross.T @ alpha

        v = torch.linalg.solve_triangular(L, sqrt_W[:, None] * K_cross, upper=False)
        f_var = (k_test_diag - v.pow(2).sum(0)).clamp(0)
        return f_mean, f_var


def train_laplace_nbgp(
    X_train, y_train,
    max_epochs=500, lr=0.01,
    patience=15, abs_tol=1e-4,
):
    """
    Train LaplaceNBGP with exact IFT gradients.
    Note: Always perform strict hold-out validation to check data leakage
    and verify length-scale/theta non-identifiability.
    """
    dev = X_train.device
    theta_init_val = _estimate_theta_mom(y_train)
    model = LaplaceNBGP(X_train.shape[-1], theta_init=theta_init_val).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    f_prev = torch.zeros(len(y_train), device=dev)

    best_loss, no_improve = float("inf"), 0
    for _ in range(max_epochs):
        model.train()
        opt.zero_grad()
        K = model.kernel(X_train)
        f_map = model.find_map(K, y_train, f0=f_prev)
        try:
            loss = -model.log_marginal_lik(K, f_map, y_train)
        except Exception:
            f_prev = f_map.detach()
            continue

        if torch.isnan(loss):
            f_prev = f_map.detach()
            continue

        loss.backward()
        opt.step()
        f_prev = f_map.detach()

        v = loss.item()
        if best_loss == float("inf"):
            best_loss = v
            continue

        if (best_loss - v) > abs_tol:
            best_loss, no_improve = v, 0
        else:
            no_improve += 1
            if no_improve >= patience:
                break

    model.eval()
    with torch.no_grad():
        K = model.kernel(X_train)
        f_map = model.find_map(K, y_train, f0=f_prev)

    return model, f_map.detach()


# =====================================================================
# Z-score functions
# =====================================================================
def pearson_zscore(y: np.ndarray, mu: np.ndarray, theta) -> np.ndarray:
    """Pearson residual under NB(mu, theta). Asymptotically N(0,1). Reproducible.

    `theta` may be a scalar or a per-observation array (e.g. an
    effective theta from `_effective_theta` that already encodes extra
    predictive variance from latent-function uncertainty).
    """
    var = mu + mu ** 2 / np.maximum(theta, 1e-6)
    return (y - mu) / np.sqrt(var + 1e-8)


def quantile_zscore(
    y: np.ndarray, mu: np.ndarray, theta, seed: int | None = None
) -> np.ndarray:
    """Randomized quantile residual (Dunn & Smyth 1996).
    Exactly N(0,1) under the correct NB model.
    Stochastic — pass seed for reproducibility within a run.

    scipy nbinom(n, p) parameterization: p = theta / (theta + mu)

    `theta` may be a scalar or a per-observation array (see `pearson_zscore`).
    """
    theta = np.maximum(theta, 1e-4)
    p = np.clip(theta / (theta + mu + 1e-8), 1e-8, 1 - 1e-8)
    y_int = y.astype(int)

    # CDF at y-1: for y=0, CDF(-1) = 0 by convention
    a = np.where(y_int > 0, nbinom.cdf(y_int - 1, n=theta, p=p), 0.0)
    b = nbinom.cdf(y_int, n=theta, p=p)

    rng = np.random.default_rng(seed)
    u = rng.uniform(
        np.clip(a, 1e-8, 1 - 1e-8),
        np.clip(b, 1e-8, 1 - 1e-8),
    )
    return norm.ppf(u)


def effective_theta(mu: np.ndarray, mu2: np.ndarray, theta: float) -> np.ndarray:
    """Method-of-moments NB theta that absorbs latent-function uncertainty.

    Given posterior samples of mu_f = exp(f) with mean `mu = E[mu_f]` and
    second moment `mu2 = E[mu_f^2]`, the total predictive variance is

        Var[y] = E[Var(y|f)] + Var[E(y|f)]
               = (mu + mu2/theta) + (mu2 - mu^2)

    `theta_eff` is the NB dispersion that reproduces this total variance
    under NB(mu, theta_eff), i.e. solves mu + mu^2/theta_eff = Var[y].
    When the posterior is a point mass (mu2 == mu^2), theta_eff == theta.
    """
    latent_var = np.clip(mu2 - mu ** 2, 0.0, None)
    denom = mu2 / max(theta, 1e-6) + latent_var
    theta_eff = np.where(denom > 1e-8, mu ** 2 / np.maximum(denom, 1e-8), theta)
    return theta_eff


def compute_zscore(
    y: np.ndarray,
    mu: np.ndarray,
    theta,
    method: str = "pearson",
    seed: int | None = None,
) -> np.ndarray:
    """method: 'pearson' | 'quantile'

    `theta` may be a scalar (no latent uncertainty) or a per-observation
    array of effective thetas (see `effective_theta`).
    """
    y, mu = np.asarray(y, dtype=np.float32), np.asarray(mu, dtype=np.float32)
    if theta is None:
        # GLM fell back to intercept-only / non-convergent fit — no usable theta
        return np.full_like(y, np.nan, dtype=np.float32)
    if method == "pearson":
        return pearson_zscore(y, mu, theta)
    if method == "quantile":
        return quantile_zscore(y, mu, theta, seed=seed)
    raise ValueError(f"Unknown method: {method!r}. Choose 'pearson' or 'quantile'.")


# ── Per-model scoring helpers ─────────────────────────────────────
def score_glm(
    glm: NBGLM, X_all: np.ndarray, y_all: np.ndarray,
    method: str = "pearson", seed: int | None = None,
):
    """Returns (z, mu_pred, theta)."""
    mu, theta = glm.predict(X_all)
    return compute_zscore(y_all, mu, theta, method=method, seed=seed), mu, theta


def score_bayes_glm(
    bglm: BayesianNBGLM, beta_map: torch.Tensor, neg_H: torch.Tensor,
    X_all_t: torch.Tensor, y_all: np.ndarray,
    method: str = "pearson", n_samples: int = 300, seed: int | None = None,
):
    """MC scoring by sampling beta ~ q(beta) = N(beta_MAP, neg_H^{-1}).

    Samples the full coefficient vector (dim p=11), then projects to f* = X*_aug @ beta.
    Avoids the N*×N* covariance matrix; memory cost is only p×S.
    """
    if seed is not None:
        torch.manual_seed(seed)
    with torch.no_grad():
        Xa = bglm._aug(X_all_t)
        p = beta_map.shape[0]
        dev = beta_map.device

        L = torch.linalg.cholesky(neg_H + 1e-6 * torch.eye(p, device=dev))
        eps = torch.randn(p, n_samples, device=dev)
        # Sigma_q = L^{-1} L^{-T}  →  sample: beta = beta_MAP + L^{-T} eps
        d_beta = torch.linalg.solve_triangular(L.T, eps, upper=True)  # (p, S)
        betas = beta_map[:, None] + d_beta  # (p, S)
        f_samp = Xa @ betas  # (N_all, S)

    mu_f = f_samp.exp()  # (N_all, S)
    mu = mu_f.mean(1).cpu().numpy()
    mu2 = mu_f.pow(2).mean(1).cpu().numpy()
    theta = bglm.theta.item()
    theta_eff = effective_theta(mu, mu2, theta)
    return compute_zscore(y_all, mu, theta_eff, method=method, seed=seed), mu, theta


def score_laplace(
    lap_model: LaplaceNBGP, f_map: torch.Tensor,
    X_tr_t: torch.Tensor, y_tr_t: torch.Tensor,
    X_all_t: torch.Tensor, y_all: np.ndarray,
    method: str = "pearson", n_samples: int = 300, seed: int | None = None,
):
    """Returns (z, mu_pred, theta). mu_pred estimated via Laplace posterior MC."""
    with torch.no_grad():
        K_tr = lap_model.kernel(X_tr_t)
        K_cross = lap_model.kernel(X_tr_t, X_all_t)
        k_diag = lap_model.kernel(X_all_t).diagonal()
        f_mean, f_var = lap_model.posterior_predictive(K_tr, f_map, y_tr_t, K_cross, k_diag)
    f_samples = torch.distributions.Normal(
        f_mean, f_var.sqrt().clamp(1e-6)
    ).rsample(torch.Size([n_samples]))
    mu_f = f_samples.exp()  # (S, N_all)
    mu = mu_f.mean(0).cpu().numpy()
    mu2 = mu_f.pow(2).mean(0).cpu().numpy()
    theta = lap_model.theta.item()
    theta_eff = effective_theta(mu, mu2, theta)
    return compute_zscore(y_all, mu, theta_eff, method=method, seed=seed), mu, theta
