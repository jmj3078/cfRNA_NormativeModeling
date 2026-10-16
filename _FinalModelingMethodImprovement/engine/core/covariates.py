"""Covariate matrix of the active pipeline: the columns of config.COVARIATES, in that order.

"W_<i>": RUVg factor i of the config.CONTROL_GENES (RUVSeq::RUVg(isLog=TRUE, center=TRUE),
0_build_covariates.ipynb section 3), i.e. the left singular vectors of the training-HC log2 CPM
control matrix centred by its gene means. fit() sees the training HC only; transform() projects
any sample from its own counts, W = (log2 CPM - HC mean) V D^-1, so injected counts move W too.
Any other name: the obs column as it is.
"""
import numpy as np

import _FinalModelingMethodImprovement.config as config

OBS_COLUMNS = [c for c in config.COVARIATES if not c.startswith("W_")]
W_INDEX = [int(c[2:]) - 1 for c in config.COVARIATES if c.startswith("W_")]


def control_logcpm(Y, idx):
    """Y: raw counts over all genes (library size = row sum)."""
    return np.log2(Y[:, idx] / Y.sum(1)[:, None] * 1e6 + 1)


class Covariates:
    def __init__(self, gene_col):
        self.idx = [gene_col[g] for g in config.CONTROL_GENES] if W_INDEX else []

    def fit(self, Y_hc):
        if W_INDEX:
            x = control_logcpm(Y_hc, self.idx)
            self.mean = x.mean(0)
            _, d, vt = np.linalg.svd(x - self.mean, full_matrices=False)
            k = max(W_INDEX) + 1
            self.v, self.d = vt[:k].T, d[:k]
        return self

    def transform(self, obs, Y):
        """obs: rows x OBS_COLUMNS; Y: the same rows' counts over all genes."""
        cols = dict(zip(OBS_COLUMNS, np.asarray(obs, dtype=np.float64).T))
        if W_INDEX:
            W = (control_logcpm(Y, self.idx) - self.mean) @ self.v / self.d
            cols.update({f"W_{i + 1}": W[:, i] for i in W_INDEX})
        return np.column_stack([cols[c] for c in config.COVARIATES]).astype(np.float64)
