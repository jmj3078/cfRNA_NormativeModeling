"""
Gene selection methods for cfRNA normative modeling.

Usage
-----
    from gene_selectors import GeneSelector

    gs = GeneSelector(Z_dis, dis_pheno, gene_names)

    genes = gs.proportion(n_per_pheno=30)

    gs.set_id2sym(adata.var['GeneName'])
    genes_by_pheno = gs.proportion_per_pheno(n_per_pheno=30)

    SELECTORS = gs.get_selectors(n_per_pheno=30)
"""

import numpy as np
from scipy.linalg import svd
from sklearn.linear_model import LogisticRegression


class GeneSelector:
    def __init__(self, Z_dis, dis_pheno, gene_names):
        self.Z_dis = Z_dis
        self.pheno = np.array(dis_pheno)
        self.gn = np.array(gene_names)
        self.phenos = np.unique(self.pheno)
        self._id2sym = None

    def set_id2sym(self, mapping):
        if hasattr(mapping, 'to_dict'):
            self._id2sym = mapping.to_dict()
        else:
            self._id2sym = dict(mapping)
        return self

    def _to_sym(self, ids):
        if self._id2sym is None:
            return ids
        return np.array([self._id2sym.get(g, g) for g in ids])

    def proportion(self, n_per_pheno=30, thr=3.0):
        """Top-N per phenotype by proportion of samples with |z| >= thr."""
        flagged = np.abs(self.Z_dis) >= thr
        selected = set()
        for ph in self.phenos:
            m = self.pheno == ph
            prop = flagged[m].mean(axis=0)
            selected.update(self._to_sym(self.gn[np.argsort(-prop)[:n_per_pheno]]))
        return sorted(selected)

    def effect_size(self, n_per_pheno=30):
        """Top-N per phenotype by mean |z| (effect size vs HC baseline=0)."""
        selected = set()
        for ph in self.phenos:
            m = self.pheno == ph
            selected.update(self._to_sym(
                self.gn[np.argsort(-np.abs(self.Z_dis[m]).mean(axis=0))[:n_per_pheno]]
            ))
        return sorted(selected)

    def svd_signature(self, n_per_pheno=30):
        """Top-N per phenotype via SVD of phenotype centroids."""
        centroids = np.array([self.Z_dis[self.pheno == ph].mean(axis=0)
                              for ph in self.phenos])
        U, S, Vh = svd(centroids, full_matrices=False)
        selected = set()
        for i in range(len(self.phenos)):
            direction = U[i, :] @ np.diag(S) @ Vh
            selected.update(self._to_sym(
                self.gn[np.argsort(-np.abs(direction))[:n_per_pheno]]
            ))
        return sorted(selected)

    def _centroids(self):
        return np.array([self.Z_dis[self.pheno == ph].mean(axis=0) for ph in self.phenos])

    def _specific_scores(self, ubiquity_penalty=True, ubiquity_abs_z=0.5):
        """Per-phenotype specificity score: |centroid_ph - mean(other-disease centroids)|.

        The contrast is vs the other diseases (not vs HC=0), so a gene perturbed in this
        phenotype only ranks high while one perturbed across all diseases cancels out. When
        ubiquity_penalty, the score is further damped by (1 - ubiquity) so a gene that is
        still anomalous in most phenotypes is suppressed (integrated approach 1 then 2).
        Returns {phenotype: score_vec}.
        """
        cents = self._centroids()
        ubiq = self.compute_ubiquity(abs_z_thr=ubiquity_abs_z) if ubiquity_penalty else None
        out = {}
        for i, ph in enumerate(self.phenos):
            others = np.delete(cents, i, axis=0).mean(axis=0)
            score = np.abs(cents[i] - others)
            if ubiquity_penalty:
                score = score * (1.0 - ubiq)
            out[ph] = score
        return out

    def effect_size_specific(self, n_per_pheno=30, ubiquity_penalty=True, ubiquity_abs_z=0.5):
        """Top-N per phenotype by disease-specific effect size (contrast vs other diseases,
        then ubiquity damping). Cuts the broadly-anomalous genes the vs-HC effect_size keeps."""
        scores = self._specific_scores(ubiquity_penalty, ubiquity_abs_z)
        selected = set()
        for ph in self.phenos:
            selected.update(self._to_sym(self.gn[np.argsort(-scores[ph])[:n_per_pheno]]))
        return sorted(selected)

    def _l1_scores(self, C=0.1):
        out = {}
        for ph in self.phenos:
            y = (self.pheno == ph).astype(int)
            if y.sum() < 2 or (1 - y).sum() < 2:
                out[ph] = None
                continue
            lr = LogisticRegression(penalty='l1', solver='liblinear', C=C, max_iter=500)
            lr.fit(self.Z_dis, y)
            out[ph] = np.abs(lr.coef_[0])
        return out

    def l1_logistic(self, n_per_pheno=30, C=0.1):
        """Top-N per phenotype by |coef| from an OVR L1 logistic separating each disease
        from the OTHER diseases (HC excluded). Discriminative -> favors phenotype-specific
        genes; a gene perturbed across all diseases gets ~0 weight. Supervised: in
        validation it MUST be refit inside each CV fold (selection._select_idx builds the
        selector on the train rows only) to avoid selection-bias leakage."""
        scores = self._l1_scores(C)
        selected = set()
        for ph in self.phenos:
            if scores[ph] is None:
                continue
            selected.update(self._to_sym(self.gn[np.argsort(-scores[ph])[:n_per_pheno]]))
        return sorted(selected)

    def proportion_per_pheno(self, n_per_pheno=30, thr=3.0):
        """Returns {phenotype: [genes]} sorted by flagging proportion."""
        flagged = np.abs(self.Z_dis) >= thr
        result = {}
        for ph in self.phenos:
            m = self.pheno == ph
            prop = flagged[m].mean(axis=0)
            result[ph] = self._to_sym(self.gn[np.argsort(-prop)[:n_per_pheno]]).tolist()
        return result

    def effect_size_per_pheno(self, n_per_pheno=30):
        """Returns {phenotype: [genes]} sorted by mean |z|."""
        result = {}
        for ph in self.phenos:
            m = self.pheno == ph
            idx = np.argsort(-np.abs(self.Z_dis[m]).mean(axis=0))[:n_per_pheno]
            result[ph] = self._to_sym(self.gn[idx]).tolist()
        return result

    def svd_signature_per_pheno(self, n_per_pheno=30):
        """Returns {phenotype: [genes]} via SVD phenotype direction."""
        centroids = np.array([self.Z_dis[self.pheno == ph].mean(axis=0)
                              for ph in self.phenos])
        U, S, Vh = svd(centroids, full_matrices=False)
        result = {}
        for i, ph in enumerate(self.phenos):
            direction = U[i, :] @ np.diag(S) @ Vh
            result[ph] = self._to_sym(
                self.gn[np.argsort(-np.abs(direction))[:n_per_pheno]]
            ).tolist()
        return result

    def effect_size_specific_per_pheno(self, n_per_pheno=30, ubiquity_penalty=True,
                                       ubiquity_abs_z=0.5):
        """Returns {phenotype: [genes]} by disease-specific effect size."""
        scores = self._specific_scores(ubiquity_penalty, ubiquity_abs_z)
        return {ph: self._to_sym(self.gn[np.argsort(-scores[ph])[:n_per_pheno]]).tolist()
                for ph in self.phenos}

    def l1_logistic_per_pheno(self, n_per_pheno=30, C=0.1):
        """Returns {phenotype: [genes]} by OVR L1 logistic |coef| (disease vs other diseases)."""
        scores = self._l1_scores(C)
        return {ph: (self._to_sym(self.gn[np.argsort(-scores[ph])[:n_per_pheno]]).tolist()
                     if scores[ph] is not None else [])
                for ph in self.phenos}

    def top_ranked_per_pheno(self, n_per_pheno=30, method='proportion', **kwargs):
        """Unified interface for per-phenotype selection.

        method : 'proportion' | 'effect_size' | 'svd' | 'effect_size_specific' | 'l1'
        """
        dispatch = {
            'proportion': self.proportion_per_pheno,
            'effect_size': self.effect_size_per_pheno,
            'svd': self.svd_signature_per_pheno,
            'effect_size_specific': self.effect_size_specific_per_pheno,
            'l1': self.l1_logistic_per_pheno,
        }
        if method not in dispatch:
            raise ValueError(f"method must be one of {list(dispatch)}")
        return dispatch[method](n_per_pheno=n_per_pheno, **kwargs)

    def get_selectors(self, n_per_pheno=30):
        """Returns SELECTORS dict for the eval loop in gene_selection.ipynb."""
        n = n_per_pheno
        return {
            f'proportion_top{n}': lambda: self.proportion(n_per_pheno=n),
            f'effect_size_top{n}': lambda: self.effect_size(n_per_pheno=n),
            f'svd_top{n}': lambda: self.svd_signature(n_per_pheno=n),
            f'effect_specific_top{n}': lambda: self.effect_size_specific(n_per_pheno=n),
            f'l1_logistic_top{n}': lambda: self.l1_logistic(n_per_pheno=n),
        }

    def compute_ubiquity(self, abs_z_thr=0.5):
        """Per-gene fraction of phenotypes where |mean_Z| > abs_z_thr.

        Returns array of shape (n_genes,) with values in [0, 1].
        Genes close to 1.0 are anomalous across nearly all disease phenotypes
        regardless of direction — non-specific signal candidates.
        """
        counts = np.zeros(self.Z_dis.shape[1], dtype=np.float32)
        for ph in self.phenos:
            m = self.pheno == ph
            counts += np.abs(self.Z_dis[m].mean(axis=0)) > abs_z_thr
        return counts / len(self.phenos)

    def mean_z_ranking(self, phenotype, jitter=1e-7, seed=42,
                       ubiquity_thr=None, ubiquity_abs_z=0.5):
        """Returns {symbol: mean_z} for a given phenotype — GSEA prerank input.

        jitter        : tiny Gaussian noise to break ties among 0-filled genes.
        ubiquity_thr  : float in (0, 1] or None. When set, genes whose
                        cross-disease ubiquity score >= ubiquity_thr are zeroed
                        out before ranking. This suppresses non-specific signal
                        (genes anomalous in most diseases regardless of direction)
                        from dominating the GSEA leading edge.
                        Recommended starting value: 0.5 (anomalous in >50% of
                        phenotypes). Set to None (default) to disable.
        ubiquity_abs_z: |mean_Z| threshold per phenotype used inside
                        compute_ubiquity(). Default 1.0.
        """
        m = self.pheno == phenotype
        mean_z = self.Z_dis[m].mean(axis=0).copy()
        if ubiquity_thr is not None:
            ubiq = self.compute_ubiquity(abs_z_thr=ubiquity_abs_z)
            mean_z[ubiq >= ubiquity_thr] = 0.0
        if jitter > 0:
            rng = np.random.default_rng(seed)
            mean_z = mean_z + rng.normal(0, jitter, len(mean_z))
        syms = self._to_sym(self.gn).tolist()
        return dict(zip(syms, mean_z.tolist()))
