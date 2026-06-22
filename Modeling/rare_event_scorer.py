"""Rare Event Scorer for genes with HC detection rate < 1%.

HC에서 거의 또는 전혀 발현되지 않는 유전자에 대해
disease 샘플의 발현량을 통계적 점수로 변환합니다.

Scoring rules
─────────────
  det_rate_hc == 0%      → z = Z_MAX  (완전 침묵 유전자, 발현 자체가 RARE EVENT)
  0% < det_rate_hc < 1%  → z = clip(Φ⁻¹(1 - Poisson.sf(k-1, λ)), Z_MAX)
                            λ = mean_count_hc (공변량 미보정)
  disease count == 0     → z = 0.0

공변량 미보정 이유:
  HC positive 샘플이 0~6개로 회귀계수 추정 불가.
  disease overexpression 규모(수백~수만 배)가 배치 효과를 압도함.

Usage
─────
    scorer = RareEventScorer.from_h5ad(H5AD_PATH)
    scorer.save(SAVE_DIR / "rare_event_ref.pkl")

    # 이후 재사용
    scorer = RareEventScorer.load(SAVE_DIR / "rare_event_ref.pkl")
    results = scorer.score(Y_disease, sample_names, gene_names_all)
"""

from __future__ import annotations

import pickle
import warnings
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import norm, poisson

warnings.filterwarnings("ignore")

Z_MAX          = 10.0   
DET_RATE_MAX   = 0.01   
BIAS_COLUMNS = [
    "log(Total Reads)", "Spliced Reads (%)",
    "gDNA Contamination (Intron/Exon)", "rRNA Fraction",
    "RNA Degradation (3' Bias)", "Platelet Score",
    "GC Bias", "Gene Length Bias", "NG80", "(NP80/NG80)",
]

def _rare_score(count: float, mean_hc: float, det_rate_hc: float) -> float:
    """단일 (유전자, 샘플) 쌍의 rare event z-score 계산."""
    if count <= 0:
        return 0.0
    if det_rate_hc == 0.0 or mean_hc == 0.0:
        return Z_MAX
    p_val = poisson.sf(int(count) - 1, mean_hc)
    if p_val <= 0.0:
        return Z_MAX
    return float(min(norm.ppf(1.0 - p_val), Z_MAX))

class RareEventScorer:
    def __init__(self, ref: pd.DataFrame):
        self.ref = ref.copy()

    @classmethod
    def from_h5ad(
        cls,
        h5ad_path: Path,
        det_rate_max: float = DET_RATE_MAX,
    ) -> "RareEventScorer":
        print(f"Loading h5ad: {h5ad_path.name}")
        adata = sc.read_h5ad(h5ad_path)
        adata = adata[adata.obs["QC_Passed"] == True]
        adata = adata[adata.obs["Phenotype_Processed"].notna()]
        adata = adata[adata.obs["Phenotype_Processed"] != "Unknown"]
        adata = adata[adata.obs["broad_protocol_category"] != "Exome-based (EB)"]

        is_hc = (adata.obs["Phenotype_Processed"].astype(str) == "Healthy Control").values
        print(f"  HC samples: {is_hc.sum()}")

        Y = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
        Y = np.round(Y).astype(np.float32)
        Y_hc = Y[is_hc]

        is_pc = (adata.var["GeneType"] == "protein_coding").values
        gene_names = adata.var_names[is_pc].tolist()
        Y_hc_pc = Y_hc[:, is_pc]

        det_rate = (Y_hc_pc > 0).mean(axis=0)
        mean_cnt = Y_hc_pc.mean(axis=0)

        rare_mask = det_rate < det_rate_max
        n_rare    = rare_mask.sum()
        print(f"  Rare genes (det < {det_rate_max:.0%}): {n_rare:,}")
        print(f"    det=0% (silent):      {(det_rate[rare_mask] == 0).sum():,}")
        print(f"    det 0-1% (near-silent): {(det_rate[rare_mask] > 0).sum():,}")

        ref = pd.DataFrame({
            "gene":        np.array(gene_names)[rare_mask],
            "det_rate_hc": det_rate[rare_mask].astype(np.float32),
            "mean_hc":     mean_cnt[rare_mask].astype(np.float32),
        })
        ref["category"] = np.where(
            ref["det_rate_hc"] == 0, "silent", "near_silent"
        )
        return cls(ref)

    # ── Persistence ─────────────────────────────────────────────
    def save(self, path: Path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump(self.ref, f)
        print(f"Saved reference → {path}")

    @classmethod
    def load(cls, path: Path) -> "RareEventScorer":
        with open(path, "rb") as f:
            ref = pickle.load(f)
        print(f"Loaded reference: {len(ref):,} rare genes")
        return cls(ref)
    # ── Scoring ─────────────────────────────────────────────────
    def score(
        self,
        Y_disease:    np.ndarray,
        sample_names: List[str],
        gene_names_all: List[str],
        only_nonzero: bool = True,
    ) -> pd.DataFrame:
        gene_to_col = {g: i for i, g in enumerate(gene_names_all)}

        records = []
        for _, row in self.ref.iterrows():
            col = gene_to_col.get(row["gene"])
            if col is None:
                continue
            counts = Y_disease[:, col]

            for si, cnt in enumerate(counts):
                if only_nonzero and cnt <= 0:
                    continue
                z = _rare_score(float(cnt), float(row["mean_hc"]),
                                float(row["det_rate_hc"]))
                score_type = (
                    "rare_fixed"   if row["category"] == "silent"     else
                    "rare_poisson" if cnt > 0                          else
                    "rare_zero"
                )
                records.append({
                    "sample":      sample_names[si],
                    "gene":        row["gene"],
                    "raw_count":   float(cnt),
                    "score":       round(z, 4),
                    "score_type":  score_type,
                    "det_rate_hc": float(row["det_rate_hc"]),
                    "mean_hc":     float(row["mean_hc"]),
                    "category":    row["category"],
                })

        df = pd.DataFrame(records)
        if df.empty:
            return df
        return df.sort_values(["sample", "score"], ascending=[True, False]).reset_index(drop=True)

    def score_vector(
        self,
        counts: np.ndarray,
        gene_names_all: List[str],
        sample_name: str = "sample",
    ) -> pd.DataFrame:
        return self.score(
            counts[np.newaxis, :],
            [sample_name],
            gene_names_all,
        )

    # ── Summary ─────────────────────────────────────────────────

    def summary(self) -> pd.DataFrame:
        return pd.DataFrame({
            "category":     self.ref["category"].value_counts(),
            "det_rate_mean": self.ref.groupby("category")["det_rate_hc"].mean().round(4),
            "mean_hc_mean":  self.ref.groupby("category")["mean_hc"].mean().round(4),
        })

if __name__ == "__main__":
    import argparse

    BASE_DIR  = Path(__file__).resolve().parent
    DATA_DIR  = BASE_DIR.parent / "OpenAccess_nfcore"
    H5AD_PATH = DATA_DIR / "Merged_Processed_AnnData_with_Batch_Biases_QC_Status.h5ad"
    SAVE_DIR  = BASE_DIR / "CV_Results"

    parser = argparse.ArgumentParser(description="Build HC rare-event reference stats.")
    parser.add_argument("--h5ad",          type=Path, default=H5AD_PATH)
    parser.add_argument("--save-dir",      type=Path, default=SAVE_DIR)
    parser.add_argument("--det-rate-max",  type=float, default=DET_RATE_MAX)
    args = parser.parse_args()

    scorer = RareEventScorer.from_h5ad(args.h5ad, args.det_rate_max)
    scorer.save(args.save_dir / "rare_event_ref.pkl")

    print("\nReference summary:")
    print(scorer.summary().to_string())
    print(f"\nSample scores for a zero vector (expect all 0.0):")
    dummy = np.zeros(20097)
    dummy_genes = scorer.ref["gene"].tolist() + ["DUMMY"] * (20097 - len(scorer.ref))
    res = scorer.score_vector(dummy, dummy_genes)
    print(f"  Non-zero entries: {len(res)}")
