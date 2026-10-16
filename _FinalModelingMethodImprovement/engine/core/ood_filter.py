import numpy as np
from numpy.linalg import inv
from scipy.spatial.distance import mahalanobis


class MahalanobisFilter:
    """Flags samples whose covariate (X) profile falls outside the HC
    distribution the model was trained on -- a sample far from HC in Mahalanobis
    distance can show a large |Z| purely from covariate extrapolation, not
    disease biology. Threshold is an HC-empirical percentile (ported from
    _legacy/GLM_Modeling_Cascade_v1/Modeling/sample_filter.py, plotting dropped)."""

    def __init__(self, percentile=95, reg=1e-8):
        self.percentile = percentile
        self.reg = reg
        self._fitted = False

    def fit(self, X_hc):
        X_hc = np.asarray(X_hc, dtype=float)
        p = X_hc.shape[1]
        self.mu_ = X_hc.mean(axis=0)
        self.cov_inv_ = inv(np.cov(X_hc.T) + np.eye(p) * self.reg)
        hc_dist = self._distances(X_hc)
        self.threshold_ = float(np.percentile(hc_dist, self.percentile))
        self._fitted = True
        return self

    def _distances(self, X):
        return np.array([mahalanobis(x, self.mu_, self.cov_inv_) for x in X])

    def distances(self, X):
        if not self._fitted:
            raise RuntimeError("Call fit(X_hc) first.")
        return self._distances(np.asarray(X, dtype=float))

    def mask(self, X):
        """Boolean keep-mask: True = inlier (within threshold)."""
        return self.distances(X) <= self.threshold_


class RangeFilter:
    """Flags samples with >= n_out_thr covariates individually outside the HC
    [lo_pct, hi_pct] range. A quadratic-form distance (Mahalanobis or a
    regression-coefficient-weighted variant) dilutes a few mildly-extreme axes
    across many normal ones; a per-axis count catches exactly that pattern --
    tested empirically (2026-07-31) to track actual disease-sample |mean Z|
    shift far better than Mahalanobis distance (r=0.43 vs 0.16)."""

    def __init__(self, n_out_thr=2, lo_pct=1, hi_pct=99):
        self.n_out_thr = n_out_thr
        self.lo_pct = lo_pct
        self.hi_pct = hi_pct
        self._fitted = False

    def fit(self, X_hc):
        X_hc = np.asarray(X_hc, dtype=float)
        self.lo_ = np.percentile(X_hc, self.lo_pct, axis=0)
        self.hi_ = np.percentile(X_hc, self.hi_pct, axis=0)
        self._fitted = True
        return self

    def n_out(self, X):
        if not self._fitted:
            raise RuntimeError("Call fit(X_hc) first.")
        X = np.asarray(X, dtype=float)
        return ((X < self.lo_) | (X > self.hi_)).sum(axis=1)

    def mask(self, X):
        """Boolean keep-mask: True = inlier (fewer than n_out_thr axes out of range)."""
        return self.n_out(X) < self.n_out_thr
