# cfRNA Normative Modeling — Project Guidelines

## Visualization Style (MANDATORY)

모든 시각화 코드에서 반드시 아래 스타일을 적용할 것.

### 노트북/스크립트에서 적용 방법

```python
from viz_style import apply_style
apply_style()
```

`viz_style.py`는 `Modeling/` 디렉토리에 위치. 모든 시각화 셀 실행 전 import 및 호출 필수.
특별한 언급이 없다면 시각화 코드에 fontsize, fontweight 선언 금지
### 스타일 내용

```python
matplotlib.rcParams["font.family"] = "Arial"
plt.rcParams.update({
    'font.family':        'Arial',
    'figure.dpi':         100,
    'savefig.dpi':        300,
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'lines.linewidth':    1.0,
    'axes.linewidth':     1.0,
    'grid.linewidth':     1.0,
    'axes.titlesize':     15,
    'axes.labelsize':     13,
    'xtick.labelsize':    13,
    'ytick.labelsize':    13,
    'legend.fontsize':    13,
})
```

### 적용 시점

- 노트북: imports 셀 하단 또는 별도 style 셀
- Python 스크립트: import 블록 직후
- 새로운 시각화 코드를 작성할 때마다 이 스타일이 적용되어 있는지 확인

---

## Project Overview

cfRNA normative modeling pipeline for disease sample anomaly detection.

- **HC reference**: WTS-only, 693 samples (Exome-based excluded)
- **Disease samples**: 913 raw → 854 (OOD P95 + MIN_SAMPLES≥5 후), 20 phenotypes
- **Branching**: det≥10% → NBI GAMLSS | 1-10% → Logistic | <1% → Rare Event Scorer
- **Working directory**: `Modeling/`
- **Key outputs**: `CV_Results/Z_disease.npy`, `CV_Results/Z_hc.npy`, `engine_state/`

---

## Analysis Pipeline (실행 순서 및 의존성)

```
disease_scoring.ipynb              # GAMLSS 엔진으로 disease/HC scoring
  └→ CV_Results/Z_disease.npy, Z_sample_names.npy, Z_gene_names.npy
  └→ CV_Results/Z_hc.npy, Z_hc_names.npy   (HC 캐시; gene_selection이 재사용)
       │
       ├→ gene_selection.ipynb     # Z 기반 gene set 선택 + 3종 CV 평가
       │     · OOD 필터(MahalanobisFilter P95, 10 BIAS_COLUMNS) → MIN_SAMPLES=5 제외
       │     · SELECTORS: proportion_top30 / effect_size_top30 / svd_top30
       │     · eval_unsupervised / eval_binary / eval_multiclass (명시적 파라미터 전달)
       │     · |Z|>3 gene count 히스토그램 (CV_Results/Figures/zscore_outlier_gene_dist.png)
       │
       └→ gene_enrichment.ipynb    # GSEA prerank + theme 기반 signature 시각화
             · GSEApy prerank (KEGG_2021 · GO-BP_2023 · Reactome_2022, perm=100, FDR<0.05)
             · 출력 CSV: Modeling/GSEA/gsea_result_{phenotype}.csv  ← (GSEA_DIR = RESULTS_DIR)
             · THEMES dict: phenotype→[(label,[substrings],is_novel)]; 같은 축 반복=신호 견고
             · fig_phenotype(ph): 좌=theme 음영 lollipop / 우=lead-gene 특이성 strip
             · 출력: Modeling/GSEA/Figures/signature/signature_{ph}.png
```

### 공통 전처리 규칙 (3 노트북 동일)
- **obs_mask**: QC_Passed==True & Phenotype_Processed notna & !='Unknown' & broad_protocol_category!='Exome-based (EB)', var: GeneType=='protein_coding'
- **OOD 필터**: `MahalanobisFilter(percentile=95)`를 HC covariate에 fit, disease에 적용 (51개/913 제거)
- **MIN_SAMPLES=5**: OOD 후 n<5 phenotype 제외 (Liver Cirrhosis (RH), Liver Cancer (Block), Pancreatic Cancer (Reggiardo))
- **Study-split**: 동일 phenotype 내 multi-author → `(FirstAuthorLastName)` suffix. `_first_author()`는 `\xa0` 처리 필수
- **NaN/Inf 처리**: `np.nan_to_num(Z, nan=0, posinf=10, neginf=-10)`

### NES 해석 규약
- NES>0 = HC 대비 상향 이탈, NES<0 = 하향 이탈
- 거의 모든 만성질환에서 번역·리보솜·OXPHOS의 거대한 NES<0 배경이 공통 → **재현 판정은 NES>0만 인정**(Rework 보고서 원칙)

---

## GSEA 해석 산출물 (Modeling/GSEA/)

| 파일 | 목적 |
|---|---|
| `GSEA_Master_Report.md` | 21개 표현형별 상향/하향 signature 상세 + Novel 발견 + 요약 테이블 |
| `GSEA_Master_Report_Rework.md` | **검증된 signature가 cfRNA에서 재현되는가** 평가 (GRADE식 STRONG/PARTIAL/DIVERGENT/UNVERIFIABLE). Novelty = "비-cfRNA 연구서 검증, cfRNA서 신규 포착" |

### Signature theme 묶음 — 두 가지 독립 방식
1. **수동 큐레이션** (`gene_enrichment.ipynb`의 `THEMES` dict): 경로명 키워드 substring 매칭. **PPT/발표용** (해석 전달에 최적, 그대로 유지).
2. **데이터 기반** (`pathway_functional_clustering.ipynb`): FDR<0.05 경로를 **leading-edge 유전자 공유**(EnrichmentMap combined coefficient = 0.5·Jaccard+0.5·Overlap) 그래프로 만들고 **Leiden** 커뮤니티 검출 → 군집 자동 라벨. 수동 theme의 방법론적 교차검증용 (예: CAD_HF+ 상향이 데이터만으로 FGFR1·coagulation·collagen 모듈로 분리 → 수동 theme과 일치 확인). 출력: `GSEA/Figures/functional_clusters/emap_*.png`, `clusters_*.csv`.

- **재현 평가 도구**: Open Targets GraphQL(암 driver), PubMed/Crossref(`/paper-lookup`), 문헌 패널(Berry 2010 TB, Sweeney3 2016, Shaughnessy 2007 MM GEP70, Maynard 2003 PE sFLT1)
- **핵심 원칙**: cfRNA=abundance(≠mutation) → KRAS/TP53/APC 같은 변이 driver 부재는 정상. RTK/분비단백 driver는 포착됨.

---
