"""
Gene selection methods for cfRNA normative modeling.

Usage
-----
    from gene_selectors import GeneSelector

    gs = GeneSelector(Z_dis, dis_pheno, gene_names)

    # flat list (eval loop용) — Ensembl ID 반환
    genes = gs.proportion(n_per_pheno=30)

    # per-phenotype dict (GE 분석용) — symbol 반환 (id2sym 제공 시)
    gs.set_id2sym(adata.var['GeneName'])       # adata.var index = ENSG ID
    genes_by_pheno = gs.proportion_per_pheno(n_per_pheno=30)

    # gene_selection.ipynb SELECTORS dict 호환
    SELECTORS = gs.get_selectors(n_per_pheno=30)
"""

import numpy as np
from scipy.linalg import svd


class GeneSelector:
    def __init__(self, Z_dis: np.ndarray, dis_pheno: np.ndarray, gene_names: list):
        self.Z_dis   = Z_dis
        self.pheno   = np.array(dis_pheno)
        self.gn      = np.array(gene_names)
        self.phenos  = np.unique(self.pheno)
        self._id2sym = None   # ENSG → symbol 매핑 (set_id2sym으로 주입)

    def set_id2sym(self, mapping) -> "GeneSelector":
        if hasattr(mapping, 'to_dict'):
            self._id2sym = mapping.to_dict()
        else:
            self._id2sym = dict(mapping)
        return self

    def _to_sym(self, ids: np.ndarray) -> np.ndarray:
        if self._id2sym is None:
            return ids
        return np.array([self._id2sym.get(g, g) for g in ids])

    # ── flat list selectors ────────────────────────────────────────

    def proportion(self, n_per_pheno: int = 30, thr: float = 3.0) -> list:
        """Top-N per phenotype by proportion of samples with |z| >= thr."""
        flagged  = np.abs(self.Z_dis) >= thr
        selected = set()
        for ph in self.phenos:
            m    = self.pheno == ph
            prop = flagged[m].mean(axis=0)
            selected.update(self._to_sym(self.gn[np.argsort(-prop)[:n_per_pheno]]))
        return sorted(selected)

    def effect_size(self, n_per_pheno: int = 30) -> list:
        """Top-N per phenotype by mean |z| (effect size vs HC baseline=0)."""
        selected = set()
        for ph in self.phenos:
            m = self.pheno == ph
            selected.update(self._to_sym(
                self.gn[np.argsort(-np.abs(self.Z_dis[m]).mean(axis=0))[:n_per_pheno]]
            ))
        return sorted(selected)

    def svd_signature(self, n_per_pheno: int = 30) -> list:
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

    # ── per-phenotype dict selectors (GE 분석용) ──────────────────

    def proportion_per_pheno(self, n_per_pheno: int = 30, thr: float = 3.0) -> dict:
        """Returns {phenotype: [genes]} sorted by flagging proportion."""
        flagged = np.abs(self.Z_dis) >= thr
        result  = {}
        for ph in self.phenos:
            m    = self.pheno == ph
            prop = flagged[m].mean(axis=0)
            result[ph] = self._to_sym(self.gn[np.argsort(-prop)[:n_per_pheno]]).tolist()
        return result

    def effect_size_per_pheno(self, n_per_pheno: int = 30) -> dict:
        """Returns {phenotype: [genes]} sorted by mean |z|."""
        result = {}
        for ph in self.phenos:
            m   = self.pheno == ph
            idx = np.argsort(-np.abs(self.Z_dis[m]).mean(axis=0))[:n_per_pheno]
            result[ph] = self._to_sym(self.gn[idx]).tolist()
        return result

    def svd_signature_per_pheno(self, n_per_pheno: int = 30) -> dict:
        """Returns {phenotype: [genes]} via SVD phenotype direction."""
        centroids = np.array([self.Z_dis[self.pheno == ph].mean(axis=0)
                              for ph in self.phenos])
        U, S, Vh  = svd(centroids, full_matrices=False)
        result    = {}
        for i, ph in enumerate(self.phenos):
            direction = U[i, :] @ np.diag(S) @ Vh
            result[ph] = self._to_sym(
                self.gn[np.argsort(-np.abs(direction))[:n_per_pheno]]
            ).tolist()
        return result

    def top_ranked_per_pheno(self, n_per_pheno: int = 30,
                              method: str = 'proportion', **kwargs) -> dict:
        """Unified interface for per-phenotype selection.

        method : 'proportion' | 'effect_size' | 'svd'
        """
        dispatch = {
            'proportion':  self.proportion_per_pheno,
            'effect_size': self.effect_size_per_pheno,
            'svd':         self.svd_signature_per_pheno,
        }
        if method not in dispatch:
            raise ValueError(f"method must be one of {list(dispatch)}")
        return dispatch[method](n_per_pheno=n_per_pheno, **kwargs)

    # ── gene_selection.ipynb SELECTORS dict 호환 ──────────────────

    def get_selectors(self, n_per_pheno: int = 30) -> dict:
        """Returns SELECTORS dict for the eval loop in gene_selection.ipynb."""
        n = n_per_pheno
        return {
            f'proportion_top{n}':  lambda: self.proportion(n_per_pheno=n),
            f'effect_size_top{n}': lambda: self.effect_size(n_per_pheno=n),
            f'svd_top{n}':         lambda: self.svd_signature(n_per_pheno=n),
        }

    # ── ranking (GSEA prerank용) ───────────────────────────────────

    def mean_z_ranking(self, phenotype: str, jitter: float = 1e-7,
                       seed: int = 42) -> dict:
        """Returns {symbol: mean_z} for a given phenotype — GSEA prerank input.

        jitter: tiny Gaussian noise added to break ties (0.0-filled NaN genes).
                Default 1e-7 << typical Z-score differences (>0.01).
        """
        m      = self.pheno == phenotype
        mean_z = self.Z_dis[m].mean(axis=0)
        if jitter > 0:
            rng    = np.random.default_rng(seed)
            mean_z = mean_z + rng.normal(0, jitter, len(mean_z))
        syms   = self._to_sym(self.gn).tolist()
        return dict(zip(syms, mean_z.tolist()))
