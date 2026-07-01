# 분석 출처·재현 기록 (Provenance) — rare/DESeq2 GSEA 비교

이 문서는 `GSEA/with_rare/GSEA_Master_Report.md`(Task1)와 `Benchmark/DESeq2_vs_Normative_Report.md`(Task2)에 사용된 **모든 외부 DB·문헌 출처와 정확한 쿼리**를 추후 검증할 수 있도록 기록한다. 수집일: **2026-07-01**.

---

## 1. 데이터 출처 (내부)

| 항목 | 경로 | 비고 |
|---|---|---|
| Normative GSEA (rare 미포함) | `GSEA/no_filter/gsea_result_*.csv` | mean-Z prerank, FDR q<0.05 사전필터 |
| Normative GSEA (rare 포함) | `GSEA/with_rare/gsea_result_*.csv` | `scoring.load_z(with_rare=True)` 기반 |
| DESeq2 GSEA | `Benchmark/deseq2_gsea/gsea_result_*.csv` | PyDESeq2 within-study rank prerank |
| GSEA gene sets | `config.MODELING_PARAMS['gsea_gene_sets']` | KEGG_2021_Human · GO_Biological_Process_2023 · Reactome_2022 |
| 유의 기준 | `config.MODELING_PARAMS['gsea_fdr_thr']` = 0.05 | Term 집합 = FDR q-val < 0.05 |
| Ensembl→HGNC 심볼 맵 | `config.PATHS['merged_qc']` adata.var['GeneName'] | rare 559 유전자 심볼(GSEA lead와 동일 소스); DESeq2 csv 유전자우주는 극저발현 rare 대부분 제외하므로 사용 금지 |
| rare-분기 유전자 목록 | `engine_state/training_summary.csv` (branch=='rare', 559개) | |

---

## 2. Open Targets Platform (질병별 참조 유전자 DB)

- **엔드포인트:** `https://api.platform.opentargets.org/api/v4/graphql` (GraphQL v4, POST, 인증 불필요)
- **수집 스크립트(재현):** `Modeling/build_disease_reference.py` → `Benchmark/disease_reference/{phenotype}.json`
- **파라미터:** 질병당 `associatedTargets` **association score 상위 300** target(approvedSymbol + score)
- **질병 ID resolve 쿼리:**
  ```graphql
  { search(queryString:"<disease name>", entityNames:["disease"], page:{index:0,size:1}){ hits{ id name } } }
  ```
- **연관 유전자 쿼리:**
  ```graphql
  query($id:String!,$p:Pagination!){ disease(efoId:$id){ name
    associatedTargets(page:$p, orderByScore:"score"){ count rows{ target{approvedSymbol} score } } } }
  ```

### 표현형 → resolve된 Open Targets 질병 ID (검증용)

| Phenotype | 검색어(query) | resolve ID | OT disease name | 전체 assoc | 저장 top-N |
|---|---|---|---|--:|--:|
| CAD_HF+ | coronary artery disease | MONDO_0005010 | coronary artery disorder | 7238 | 300 |
| CAD_HF- | coronary artery disease | MONDO_0005010 | coronary artery disorder | 7238 | 300 |
| Colorectal Cancer | colorectal carcinoma | MONDO_0024331 | colorectal carcinoma | 13882 | 300 |
| Esophagus Cancer (Chen) | esophageal carcinoma | MONDO_0019086 | carcinoma of esophagus | 6967 | 300 |
| HIV | HIV infection | MONDO_0005109 | HIV infectious disease | 3738 | 300 |
| HIV + Tuberculosis | tuberculosis | MONDO_0018076 | tuberculosis | 4700 | 300 |
| Tuberculosis | tuberculosis | MONDO_0018076 | tuberculosis | 4700 | 300 |
| Liver Cancer (Chen) | hepatocellular carcinoma | MONDO_0007256 | hepatocellular carcinoma | 15470 | 300 |
| Liver Cancer (Roskams-Hieter) | hepatocellular carcinoma | MONDO_0007256 | hepatocellular carcinoma | 15470 | 300 |
| Liver Cirrhosis | liver cirrhosis | MONDO_0005155 | cirrhosis of liver | 3436 | 300 |
| Lung Cancer | lung carcinoma | MONDO_0005138 | lung carcinoma | 14357 | 300 |
| ME_CFS | chronic fatigue syndrome | MONDO_0005404 | myalgic encephalomyelitis/CFS | 1775 | 300 |
| MGUS | monoclonal gammopathy of undetermined significance | EFO_1000836 | benign monoclonal gammopathy | 45 | 45 |
| MM | multiple myeloma | MONDO_0009693 | plasma cell myeloma | 5305 | 300 |
| Pancreatic Cancer (Moore) | pancreatic carcinoma | MONDO_0005192 | exocrine pancreatic carcinoma | 10899 | 300 |
| Pancreatitis | pancreatitis | MONDO_0004982 | pancreatitis | 3058 | 300 |
| Pre-eclampsia | pre-eclampsia | MONDO_0005081 | preeclampsia | 4277 | 300 |
| Stomach Cancer | gastric carcinoma | MONDO_0004950 | gastric carcinoma | 6700 | 300 |
| Other Cancer | cancer | MONDO_0004992 | cancer | 22581 | 300 |
| ICI-m | myocarditis | MONDO_0004496 | myocarditis | 1521 | 300 |
| ICI-treated Cancer | — (해당 없음) | — | — | 0 | 0 |

> **주의:** MGUS는 전체 association 45개뿐(참조 유전자 45개). ICI-treated Cancer는 이질 혼합 코호트라 단일 질병으로 resolve 불가 → OT 참조 없음, 문헌전용.

---

## 3. 문헌 (PubMed / NCBI E-utilities)

- **엔드포인트:** `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term=<query>`
- **용도:** 각 rare-led/headline 신호의 gene×disease 선행 문헌 유무 정량(히트 수). "기지"=다수 존재, "후보(novel)"=≤2건.
- 아래는 리포트 서술의 근거가 된 **정확한 쿼리와 히트 수**(2026-07-01 기준; 시점에 따라 소폭 변동 가능).

| # | PubMed 쿼리(term) | 히트 | 리포트 판정 |
|--:|---|--:|---|
| 1 | (chronic fatigue syndrome) AND (insulin OR INS) | 116 | ME/CFS INS 기지 |
| 2 | (chronic fatigue syndrome) AND (lipid metabolism) | 80 | ME/CFS 지질 기지 |
| 3 | (chronic fatigue syndrome) AND (PLA2G10 OR phospholipase A2) | 0 | **후보** |
| 4 | (stomach neoplasms) AND (CLDN18 OR claudin 18) | 324 | 위암 CLDN18 기지(임상표적) |
| 5 | (stomach neoplasms) AND (CLDN17 OR claudin 17) | 27 | 위암 CLDN17 상대적 신규 |
| 6 | (pancreatic neoplasms) AND (PNLIP OR pancreatic lipase) | 44 | PDAC PNLIP 기지 |
| 7 | (pancreatic neoplasms) AND (myostatin OR MSTN) | 11 | PDAC MSTN 기지(악액질) |
| 8 | (pancreatic neoplasms) AND (retinoid metabolism) | 152 | PDAC retinoid 기지 |
| 9 | (multiple myeloma) AND (neutrophil extracellular trap) | 18 | MM NET 신흥 |
| 10 | (monoclonal gammopathy) AND (histone OR nucleosome) | 719 | MGUS 히스톤 기지 |
| 11 | (myocarditis) AND (metallothionein OR MT1) | 2 | **후보** |
| 12 | (myocarditis) AND (iron OR ferroptosis) | 138 | 심근염 iron 기지 |
| 13 | (pre-eclampsia) AND (ribosome OR RPL) | 98 | PE 리보솜 기지 |
| 14 | (pre-eclampsia) AND (IFNL3 OR interferon lambda) | 1 | **후보** |
| 15 | (heart failure) AND (nitric oxide OR NOS3) | 4300 | HF NO/eNOS 기지 |
| 16 | (coronary artery disease) AND (extracellular matrix OR collagen) | 2484 | CAD ECM 기지 |
| 17 | (carcinoma hepatocellular) AND (GPC3 OR glypican) | 1204 | HCC GPC3 기지 |
| 18 | (carcinoma hepatocellular) AND (DEFB114 OR beta defensin) | 3 | **후보** |
| 19 | (tuberculosis) AND (DKK4 OR Wnt signaling) | 68 | TB WNT 신흥 |
| 20 | (colorectal neoplasms) AND (interleukin 7 OR IL7) | 64 | CRC IL-7 기지 |
| 21 | (lung neoplasms) AND (ERBB OR EGFR) | 23568 | 폐암 ERBB/EGFR 기지 |
| 22 | (HIV) AND (Vpr AND nuclear import) | 112 | HIV Vpr 기지 |
| 23 | (tuberculosis) AND (complement OR MBL2) | 2274 | TB 보체 기지 |

---

## 4. 재현 절차

1. **참조 DB 재생성:** `python Modeling/build_disease_reference.py` → `Benchmark/disease_reference/*.json` (Open Targets 재조회, 수집일 갱신)
2. **정량·교차검증 재생성:** `pipeline/gsea_compare.py` (thin-runner 노트북 `_gsea_rare_deseq2_comparison.ipynb`에서 호출)
   - `compare_all()` → `Benchmark/gsea_compare/overlap_stats.csv` + 질병별 diff CSV
   - `validate_rare_novel(diffs)` → `rare_novel_validated.csv` · `rare_novel_summary.csv`
   - `deseq2_coverage()` → `deseq2_coverage.csv`
3. **문헌 검증 재현:** 위 표의 PubMed term을 E-utilities에 재조회(히트 수는 시점 의존).
4. **정성 리포트:** 위 정량 산출물 + DB/문헌 위에 수기 작성(`GSEA/with_rare/GSEA_Master_Report.md`, `Benchmark/DESeq2_vs_Normative_Report.md`).

## 5. 검증상 유의점

- Open Targets "association score"는 유전학·문헌·발현·경로 등 다중 evidence의 가중합으로, 질병관련성의 **대리지표**다(정답 아님). `datasourceScores`로 세부 근거 확인 가능.
- PubMed 히트 수는 검색어 표현·시점에 따라 변동한다. 판정 임계(≤2건=후보)는 보수적 컷오프일 뿐이다.
- 모든 GSEA term 집합은 FDR q<0.05 사전필터 결과이며, GSEA permutation 무작위성으로 경계 근처 term은 재실행 시 소폭 변동할 수 있다.
