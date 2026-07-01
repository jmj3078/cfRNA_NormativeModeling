# DESeq2 GSEA vs Normative-model GSEA — 비교 리포트

**대상:** DESeq2 rank 기반 GSEA(within-study, PyDESeq2) vs Normative-model GSEA(mean-Z rank). 18개 표현형(ICI-m/ICI-treated 제외; DESeq2 within-study 대조군 부재).
**질문:** 두 방법이 동일한 질병 관련 시그널을 검출하는가, 얼마나 검출하지 못하는가(정성) + 겹침의 정량 요약.
**유의 기준:** 양측 FDR q-val < 0.05. **NES 부호:** normative=HC 대비 이탈 방향, DESeq2=질병 vs 대조 log2FC 방향.

> **해석상 주의(객관화):** 이 비교는 우열 판정이 아니라 *일치도(concordance)* 측정이다. 두 방법은 서로 다른 대상을 랭킹한다 — normative는 샘플별 HC 규범 이탈(mean-Z), DESeq2는 집단 평균 차이(Wald stat). 랭킹 지표·유전자 우주·정규화·independent filtering이 다르므로 불일치의 일부는 방법론적 차이에서 기인한다. "DB 지지"는 Open Targets association 상위 300 유전자와 lead genes가 겹치는지를 뜻하며 질병관련성의 *대리지표*일 뿐 정답(ground truth)이 아니다. 아래 수치는 이 정의 하에서의 관측이다.

---

## 정량 요약 (표현형별)

normative(with_rare)를 기준으로, 그 중 DB-지지 유의 경로를 DESeq2가 얼마나 함께 검출(capture)하는지, 그리고 공유 경로에서 NES 부호가 얼마나 불일치하는지.

| Phenotype | DESeq2 유의 | Normative 유의 | Norm DB-지지 | DESeq2 포착 | 포착률 | 공유 term | 부호불일치율 |
|---|--:|--:|--:|--:|--:|--:|--:|
| Pre-eclampsia | 54 | 282 | 41 | 0 | 0.000 | 13 | 0.000 |
| Liver Cancer (Chen) | 107 | 656 | 320 | 5 | 0.016 | 27 | 0.481 |
| Pancreatic Cancer (Moore) | 130 | 604 | 297 | 9 | 0.030 | 57 | 0.702 |
| ME_CFS | 98 | 226 | 142 | 6 | 0.042 | 9 | 0.000 |
| Esophagus Cancer (Chen) | 127 | 746 | 290 | 16 | 0.055 | 48 | 0.062 |
| Colorectal Cancer | 160 | 441 | 186 | 12 | 0.065 | 59 | 0.000 |
| Pancreatitis | 50 | 310 | 119 | 8 | 0.067 | 22 | 0.955 |
| MGUS | 188 | 127 | 13 | 1 | 0.077 | 16 | 0.062 |
| HIV | 70 | 478 | 254 | 37 | 0.146 | 42 | 0.000 |
| Liver Cancer (Roskams-Hieter) | 287 | 431 | 314 | 52 | 0.166 | 67 | 0.015 |
| CAD_HF+ | 596 | 1085 | 291 | 61 | 0.210 | 340 | 0.021 |
| Other Cancer | 72 | 323 | 178 | 38 | 0.213 | 52 | 1.000 |
| CAD_HF- | 596 | 984 | 264 | 59 | 0.223 | 316 | 0.032 |
| Lung Cancer | 319 | 548 | 214 | 55 | 0.257 | 186 | 0.000 |
| MM | 439 | 405 | 207 | 85 | 0.411 | 195 | 0.072 |
| Tuberculosis | 436 | 672 | 446 | 184 | 0.413 | 241 | 0.000 |
| Stomach Cancer | 437 | 563 | 266 | 135 | 0.508 | 258 | 0.000 |
| HIV + Tuberculosis | 437 | 434 | 266 | 158 | 0.594 | 215 | 0.005 |

**요약 통계:** DESeq2의 DB-지지 경로 포착률 중앙값 = **0.156**(with_rare 기준; no_filter 기준 0.181). 평균 검출 term 수 DESeq2 256 vs Normative 518. 즉 DESeq2 GSEA는 normative가 찾는 DB-지지 질병관련 경로의 다수를 유의 수준에서 재현하지 않는다(관측).

> 위 표는 DB-지지를 *normative 항목에서만* 정의한 뒤 DESeq2 포착을 보는 **비대칭 관점**이다. 아래 (b)는 세 방법을 동일 규칙으로 채점한 대칭 관점을 추가한다.

## 정량 요약 (b) — 대칭 DB-hit 비교

각 방법을 **동일 규칙**으로 채점한다: 자기 자신의 FDR-유의 term 중 lead genes가 Open Targets 참조(질병별 상위 300 유전자)와 겹치는 term = DB-hit. 이렇게 하면 DESeq2·Normative의 DB-hit **개수(절대 커버리지)**와 **비율(정밀도)**을 대칭적으로 비교할 수 있다. 개수가 비율보다 중요할 수 있어 둘 다 제시한다. (참조 DB가 있는 표현형만 집계; `db_hit_rates.csv`·`db_hit_rates_summary.csv`, 그림 `Figures/db_hit_rates.png`.)

**pooled 요약 (참조 있는 표현형 전체 합산):**

| Method | 유의 term 총수 | DB-hit term수 | pooled DB-hit율 | 표현형평균율 | 표현형수 |
|---|--:|--:|--:|--:|--:|
| DESeq2 | 4,603 | 2,236 | 0.486 | 0.454 | 18 |
| DESeq2 + covariates | 1,963 | 1,143 | **0.582** | 0.437 | 18 |
| Normative (no_filter) | 8,764 | 3,843 | 0.438 | 0.445 | 20 |
| Normative (with_rare) | 9,837 | 4,375 | 0.445 | 0.450 | 19 |

**관측:** 비율(정밀도)은 DESeq2가 근소하게 높으나(0.486 vs 0.438~0.445) 세 방법이 사실상 동등(표현형평균율 0.45 안팎)하다. 반면 **DB-지지 term의 절대 개수는 Normative가 ~1.7–2배 많다**(3,843~4,375 vs 2,236) — DESeq2가 유의 term을 더 보수적으로 적게(4,603 vs 8,764~9,837) 부르기 때문. 즉 "비슷한 정밀도로 훨씬 넓은 DB-지지 커버리지"로 관측된다.

DESeq2 + covariates는 pooled DB-hit율이 0.582로 가장 높으나 유의 term 총수가 1,963으로 가장 적다 — 공변량 보정이 정밀도를 높이는 대신 커버리지를 대폭 낮추는 것으로 관측된다. 아래 섹션 (c) 참조.

**표현형별 DB-hit 개수(n_db) | DB-hit율:**

| Phenotype | DESeq2 수 | no_filter 수 | with_rare 수 | DESeq2 율 | no_filter 율 | with_rare 율 |
|---|--:|--:|--:|--:|--:|--:|
| Tuberculosis | 279 | 407 | 446 | 0.640 | 0.672 | 0.664 |
| HIV + Tuberculosis | 322 | 233 | 266 | 0.737 | 0.587 | 0.613 |
| Stomach Cancer | 268 | 248 | 266 | 0.613 | 0.481 | 0.472 |
| CAD_HF+ | 183 | 252 | 291 | 0.307 | 0.262 | 0.268 |
| CAD_HF- | 183 | 247 | 264 | 0.307 | 0.268 | 0.268 |
| Liver Cancer (Chen) | 52 | 306 | 320 | 0.486 | 0.490 | 0.488 |
| Liver Cancer (Roskams-Hieter) | 194 | 157 | 314 | 0.676 | 0.724 | 0.729 |
| Esophagus Cancer (Chen) | 77 | 274 | 290 | 0.606 | 0.401 | 0.389 |
| MM | 255 | 157 | 207 | 0.581 | 0.471 | 0.511 |
| Pancreatic Cancer (Moore) | 20 | 280 | 297 | 0.154 | 0.500 | 0.492 |
| Lung Cancer | 156 | 203 | 214 | 0.489 | 0.390 | 0.391 |
| HIV | 51 | 263 | 254 | 0.729 | 0.546 | 0.531 |
| ICI-m | - | 247 | 267 | - | 0.516 | 0.511 |
| Colorectal Cancer | 95 | 157 | 186 | 0.594 | 0.408 | 0.422 |
| Other Cancer | 22 | 149 | 178 | 0.306 | 0.560 | 0.551 |
| ME_CFS | 31 | 95 | 142 | 0.316 | 0.562 | 0.628 |
| Pancreatitis | 16 | 105 | 119 | 0.320 | 0.399 | 0.384 |
| Pre-eclampsia | 10 | 36 | 41 | 0.185 | 0.157 | 0.145 |
| MGUS | 22 | 9 | 13 | 0.117 | 0.084 | 0.102 |
| Liver Cirrhosis | - | 18 | - | - | 0.419 | - |

표현형별로는 방향이 갈린다: 소화기 암종(대장/식도/위)은 DESeq2 정밀도가 높으나 개수는 Normative가 대체로 우세, ME/CFS·Other Cancer·Pancreatic Cancer는 개수·비율 모두 Normative가 뚜렷이 우세하다(이들은 앞 (a) 표에서 DESeq2 포착률이 낮았던 소표본군과 일치). DESeq2가 없는 표현형(ICI-m·Liver Cirrhosis)은 within-study 대조군 부재로 빈칸.

## 정량 요약 (c) — DESeq2 공변량 보정 전후 비교

`run_deseq2_covariate.py`: BIAS_COLUMNS 10개를 전체 샘플 StandardScaler로 표준화 후 design `~cov1+…+cov10+condition`(Author_clean 제외; 여기서 within-study HCs가 자연 대조). 결과는 `Benchmark/gsea_compare/deseq2_cov_vs_nocov_overlap.csv` + `deseq2_cov_db_hits.csv`.

**핵심 관측: 공변량 보정이 DESeq2 GSEA 결과를 극적으로 바꾼다.**

- term-level Jaccard 중앙값 = **0.052** — 보정 전후 유의 term의 대부분이 교체됨
- 공유 term 내 sign_agree = 1.0 — 겹치는 term은 방향 일치
- 10개 표현형에서 유의 term이 0~2개로 급감(MGUS/Other Cancer/Pancreatitis/ME_CFS 등)
- CAD_HF+·CAD_HF-는 596개 → 11/2개로 급감 → 공변량 보정이 CAD 표현형의 DESeq2 시그널을 거의 제거

**보정 후 유의 term을 유지하는 표현형(n_sig ≥ 50):**

| Phenotype | DESeq2 유의 | DESeq2+cov 유의 | Jaccard | DB-hit율(cov) |
|---|--:|--:|--:|--:|
| Stomach Cancer | 437 | 541 | 0.362 | 0.577 |
| Lung Cancer | 319 | 388 | 0.354 | 0.564 |
| MM | 439 | 228 | 0.372 | 0.531 |
| HIV + Tuberculosis | 437 | 289 | 0.443 | 0.616 |
| Tuberculosis | 436 | 209 | 0.375 | 0.842 |
| HIV | 70 | 111 | 0.090 | 0.486 |

위 6개 표현형은 공변량 보정 후에도 시그널을 유지하며(Jaccard 0.09~0.44), DB-hit율이 높다(특히 Tuberculosis 0.842). 반면 나머지 12개 표현형은 보정 후 유의 term이 소멸하거나 극소화 — DESeq2의 원 시그널이 공변량 교란에 크게 의존했음을 시사한다.

이는 Normative Model이 공변량을 **각 샘플별로 사전에 모델링**하여 교란을 개별 제거하는 것과 대비된다: within-study HC 없이도 교란-보정 Z-score를 모든 표현형에서 안정적으로 산출한다.

---

## 정성 관찰 (표현형 그룹별)

포착률과 부호 불일치율을 기준으로 세 유형으로 나뉜다. 아래 서술은 관측 기술이며 인과·우열 판정을 포함하지 않는다.

### A. 상대적으로 일치가 높은 표현형 (포착률 ≥ 0.4, 부호 불일치 낮음)

HIV+TB(0.594), Stomach Cancer(0.508), Tuberculosis(0.413), MM(0.411). 이들은 DESeq2 유의 term 수가 많고(437~439) 공유 term에서 부호 불일치가 ≤0.07이다. 즉 두 방법이 **핵심 축의 방향과 상당 부분을 공유**한다 — 결핵·공동감염의 자연면역/보체(TB×complement 문헌 2274건), 위암 상피/접착, MM 형질세포 축. 이 표현형들은 상대적으로 샘플 수·효과크기가 커서 집단 평균 검정(DESeq2)의 검정력이 확보된 경우로 관측된다.

### B. DESeq2 검출이 희소한 표현형 (포착률 < 0.07)

Pre-eclampsia(0.000), Liver Cancer Chen(0.016), ME/CFS(0.042), Colorectal(0.065), Esophagus(0.055), Pancreatic(0.030), Pancreatitis(0.067). 이들은 DESeq2 유의 term이 50~160개로 적고, normative의 DB-지지 경로 대부분을 유의 수준에서 재현하지 않는다. 공통적으로 **표본 수가 작거나 within-study 대조군 구성이 제한**되어 DESeq2 집단 검정의 검정력이 낮은 것과 관측상 연관된다(예: Pre-eclampsia DESeq2 유의 54, normative DB-지지 41 중 0 포착). normative의 샘플별 규범 이탈 정량은 이 조건에서도 신호를 유지하는 것으로 관측된다.

### C. 공유 term의 방향이 반대인 표현형 (부호 불일치율 ≥ 0.5)

Other Cancer(1.000, n=52), Pancreatitis(0.955, n=22), Pancreatic Cancer(0.702, n=57). 이 경우 두 방법이 함께 유의하게 부른 경로조차 NES 부호가 반대다. 불일치 term을 조사하면 **번역·리보솜·유비퀴틴 매개 분해 등 housekeeping 축**에 집중된다:

| Phenotype | 공유 term (예시) | Normative NES | DESeq2 NES |
|---|---|--:|--:|
| Other Cancer | Nonsense Mediated Decay (EJC-independent) | −3.04 | +3.38 |
| Other Cancer | Selenocysteine Synthesis | −2.98 | +3.62 |
| Pancreatitis | Ubiquitin-mediated Degradation Of Cdc25A | −2.24 | +2.41 |
| Pancreatic Cancer | NIK/Noncanonical NF-κB Signaling | −2.49 | +1.92 |

이 축들은 라이브러리 크기 정규화·RNA 조성(composition)에 민감한 것으로 알려져 있다(DESeq2 GSEA에서 리보솜 단백 경로가 최상위 음성 NES로 지배적으로 나타나는 것과 동일 맥락). 두 방법의 정규화·랭킹 기준 차이(normative=HC 규범 대비 샘플별 이탈, DESeq2=집단 log2FC)가 **조성 민감 housekeeping 축에서 반대 부호**로 귀결되는 것으로 관측된다. 이는 핵심 질병 특이 경로의 생물학적 상충이라기보다 **기술적 축에서의 방법론적 불일치**로 해석하는 것이 보수적이다. 다만 어느 방향이 실제 생물학인지는 본 데이터만으로 판별할 수 없으며 외부 검증이 필요하다.

---

## 종합 및 한계

1. **정량:** (a) 비대칭 관점 — DESeq2 GSEA는 normative가 검출하는 DB-지지 경로의 중앙값 ~16%만 유의 수준에서 함께 검출하며, 총 검출 term 수도 절반 수준(256 vs 518)이다. 겹침(Jaccard)은 표현형별 0.03~0.36으로 낮다(overlap_stats.csv). (b) 대칭 관점 — 세 방법을 동일 규칙으로 채점하면 DB-hit **비율(정밀도)**은 사실상 동등(pooled 0.44~0.49, 표현형평균 0.45 안팎)하나 DB-지지 term **개수**는 Normative가 ~1.7–2배 많다(2,236 vs 3,843~4,375). 정밀도 손실 없이 커버리지가 넓은 것으로 관측된다(db_hit_rates.csv). (c) DESeq2 공변량 보정(+cov)은 pooled DB-hit율을 0.486→0.582로 높이나 유의 term을 4,603→1,963으로 대폭 줄인다; 20표현형 중 12개에서 시그널이 소멸/급감 — DESeq2 원 시그널의 많은 부분이 공변량 교란에 의존했음을 시사한다. 공변량 보정 전후 Jaccard 중앙값 = 0.052로 결과가 극적으로 달라진다.
2. **정성:** 일치도는 표본 규모·효과크기와 관측상 연관된다 — 대형 코호트(TB, HIV+TB, 위암, MM)에서 높고, 소표본(PE, ME/CFS, 췌장 계열)에서 낮다.
3. **방향 불일치**는 소수 표현형에서 크며, 그 대상이 번역/분해 housekeeping 축에 집중되어 조성 정규화 차이에 기인할 가능성이 관측된다.
4. **한계(객관):** (a) "DB 지지"는 Open Targets 상위 300 유전자에 대한 lead-gene 겹침으로 정의한 대리지표이며 정답이 아니다. (b) 두 GSEA는 랭킹 지표·유전자 우주·정규화가 달라 직접 비교에 방법론적 차이가 내재한다. (c) DESeq2는 within-study 대조군이 있는 18개 표현형에 한정되며, 표본 수 차이가 검정력에 영향을 준다. (d) 본 비교는 concordance 측정이며 어느 방법이 정답에 가까운지에 대한 판정을 포함하지 않는다.
