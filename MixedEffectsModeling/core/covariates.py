"""Covariate matrix of the active pipeline: the columns of config.COVARIATES, in that order.

"W_<i>": RUVg factor i of the config.CONTROL_GENES (RUVSeq::RUVg(isLog=TRUE, center=TRUE),
Notebook_p3/control_genes.ipynb), i.e. the left singular vectors of the training-HC log2 CPM
control matrix centred by its gene means. fit() sees the training HC only; transform() projects
any sample from its own counts, W = (log2 CPM - HC mean) V D^-1, so injected counts move W too.
Any other name: the obs column as it is; with Y0 (the same rows' counts before an in-silico injection)
the count-derived ones (count_metrics.COUNT_METRICS) are shifted by f(Y) - f(Y0), so they move with
the injection exactly as recomputed from counts and are unchanged where Y == Y0.
"""
import numpy as np

import MixedEffectsModeling.config as config
from MixedEffectsModeling.core.count_metrics import COUNT_METRICS, count_metrics

OBS_COLUMNS = [c for c in config.COVARIATES if not c.startswith("W_")]
W_INDEX = [int(c[2:]) - 1 for c in config.COVARIATES if c.startswith("W_")]
COUNT_COLUMNS = [c for c in OBS_COLUMNS if c in COUNT_METRICS]
DEPENDS_ON_COUNTS = bool(W_INDEX or COUNT_COLUMNS)


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

    def transform(self, obs, Y, Y0=None):
        """obs: rows x OBS_COLUMNS; Y: the same rows' counts over all genes; Y0: their counts before injection."""
        cols = dict(zip(OBS_COLUMNS, np.asarray(obs, dtype=np.float64).copy().T))
        if Y0 is not None and COUNT_COLUMNS:
            changed = np.where((Y != Y0).any(1))[0]
            if len(changed):
                new, old = count_metrics(Y[changed], COUNT_COLUMNS), count_metrics(Y0[changed], COUNT_COLUMNS)
                for c in COUNT_COLUMNS:
                    cols[c][changed] += new[c] - old[c]
        if W_INDEX:
            W = (control_logcpm(Y, self.idx) - self.mean) @ self.v / self.d
            cols.update({f"W_{i + 1}": W[:, i] for i in W_INDEX})
        return np.column_stack([cols[c] for c in config.COVARIATES]).astype(np.float64)
