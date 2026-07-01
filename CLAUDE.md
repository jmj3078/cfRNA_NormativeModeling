# CLAUDE.md
This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
항상 CLAUDE.md에 존재할 가치가 있는 줄만 유지하도록 할 것.

# cfRNA Normative Modeling — 공통 프로젝트 요구사항
### 분석 핵심 가정 및 목적
- 기존 cfRNA 전사체 분석은 주로 정상군과 질병군 간의 집단 수준 비교(Group-wise comparison)에 의존해 왔다. 그러나 생물학적 및 기술적 공변량(Covariates)에 의한 분산이 질병 고유의 신호를 압도하는 경우가 많다. 이로 인해 집단 단위의 일괄적인 공변량 보정은 질병 신호의 소실이나 교란 요인의 잔존을 초래하는 근본적인 한계가 존재한다. 이를 극복하기 위해, 본 연구는 전장 전사체(Whole-Transcriptome) 기반의 대규모 정상군(Healthy Control) 데이터를 활용한 규범적 모델링(Normative Modeling)을 도입한다. 개별 샘플의 공변량을 반영하여 정상 상태의 예상 분포를 추정하고, 이를 통계적 편차(Z-score)로 산출함으로써 교란 요인의 영향 없이 질병 특이적 신호를 정밀하게 정량화.

### 코드 작성 요령
- **간결성**: 최소한의 코드를 지향. 불필요한 추상화·방어 로직 금지.
- **타입 힌트 금지**: 함수/메서드의 입력 인자 dtype, 반환값 dtype 모두 표기하지 않는다.
- **정렬용 공백 금지**: 줄·등호를 맞추기 위한 인위적 띄어쓰기 금지 (`a = 1` O / `a   = 1` X).
- **이모지 금지**: 코드·주석·출력·문서 어디에도 사용하지 않는다.
- **주석**: 영어로만 작성. 단, 사용자가 명시적으로 요청하기 전에는 주석을 달지 않는다 (작업 최종 완료 후 일괄 추가 예정).
- **import 순서**: 알파벳 순.
- **경로·재사용 변수**: 절대 재선언하지 말고 루트 `config.py`에서 전역 import (구조는 아래 디렉토리 트리 참조).
- **시각화**: 모든 그림은 `apply_style()`로 공통 테마 적용. 노트북/스크립트 위치 기준으로 아래 패턴 사용.
  ```python
  if parent_dir not in sys.path:
      sys.path.insert(0, parent_dir)
  from viz_style import apply_style
  apply_style()
  ```

### 프로젝트 디렉토리 구조 (경로 = 간단설명 단일 소스)
유지 규칙: 모든 경로/구조 정보는 아래 트리에 **경로 : 한 줄 설명** 형식으로만 기록한다. 별도 상세 섹션을 만들어 중복시키지 말 것. config.py나 디렉토리가 바뀌면 본 트리만 갱신.
작업 환경 : conda env "scRNA" 에서 작업
./config.py : 경로·파라미터 단일 소스. 모든 코드가 재선언 없이 import. ROOT/DATA_DIR/PIPELINE_DIR · MODELING_DIR/ENGINE_DIR/CV_RESULTS_DIR/Z_SCORES_DIR/GSEA_DIR/BENCHMARK_DIR(DESEQ2_RESULTS_DIR/DESEQ2_GSEA_DIR) · H5AD_PATH · Z_DISEASE/Z_HC/Z_RARE_DISEASE/Z_RARE_HC/Z_RARE_GENE_NAMES 등 Z 경로 · RARE_GLM · PATHS{merged_raw/biases/qc} · BIAS_COLUMNS · MODELING_PARAMS{ood_percentile=95/min_samples=5/z_flag=3.0/det_rate_min=0.01/low_det_thr=0.10/rare_det_max=0.01/rare_overdisp_thr=2.0/rare_z_cap=10.0/mean_count_min=2.0/lr_c=1.0/gsea_*/emap_sim_thr=0.50} · PARAMS(EDA용)
./viz_style.py : apply_style() 공통 matplotlib 테마 (모든 시각화 필수)
./CLAUDE.md · ./README.md : 공용 문서
./_legacy/ : 폐기 파일·임시파일. 참고·수정·읽기 금지
./Data/ : 분석 핵심 데이터 (권한 없이 수정 금지)
./OpenAccess_nfcore/ : 주 분석 데이터(adata) 원본·전처리본 (권한 없이 수정 금지). config.PATHS.merged_qc = 주 h5ad
./RPM_nfcore/ : Validation 실험실 데이터 (추가 예정)
./Saved_Pipeline/ : config.PIPELINE_DIR (LogisticGP·Z matrix, 생성 예정)
./EDA/ : 코호트 QC·batch/bias 교란분석. cwd=EDA 가정, 루트 config/viz_style를 Path.cwd().parent로 sys.path 등록해 import
./EDA/analysis_cfrna_cohorts.ipynb : QC→PCA→RDA 교란분석 (helper+plot 호출)
./EDA/analysis_helper.py : QC·bias 정량화 + RDA 분산분해 엔진 (자체적으로 root를 sys.path 등록)
./EDA/analysis_plot.py : analysis_helper 전용 시각화
./EDA/VariousNormalizationMethods_OpenAccess.R : 정규화 레이어 생성 (R, PROJECT_ROOT 절대경로)
./EDA/Analysis_Results/ : 위 노트북 출력 (노트북에서 ./Analysis_Results/ 상대경로)
./Modeling/ : cfRNA Normative Modeling 본체
./Modeling/*.py : model_engine(NBI/ZINBI·logistic·rare 분기) · run_model_engine(→engine_state/, --rare-only로 rare GLM만 부착) · sample_filter(MahalanobisFilter OOD) · gene_selectors(proportion/effect_size/svd + effect_size_specific[방안1 질병간대조+방안2 ubiquity damping] + l1_logistic[OVR L1 판별]) · cv_gamlss_nb/zinb · cv_logistic · cv_rare(rare pooled GLM CV) · build_disease_reference(Open Targets 질병별 참조 유전자 JSON 재생성). 경로/BIAS_COLUMNS/임계값 전부 config import. (rare_event_scorer는 폐기→_legacy/)
./Modeling/pipeline/ : 분석·시각화 모듈 패키지. data_prep(공통 전처리:load_adata/study-split/OOD·MIN_SAMPLES/Z 로드)·scoring·selection·enrichment·signatures(THEMES+heuristic/emap 군집)·comparison·benchmark(DESeq2 vs Normative 유전자 단위 비교)·gsea_compare(GSEA term-level 비교: with_rare↔no_filter↔DESeq2 겹침 통계+diff+Open Targets DB 교차검증. load_sets/compare_all/validate_rare_novel/deseq2_coverage[비대칭:normative DB-지지 경로의 DESeq2 포착률]/db_hit_rates[대칭:세 방법 동일규칙 DB-hit 개수·비율 비교,개수 강조]) + plots(시각화 전용, plot_db_hit_rates 포함). __init__이 root를 sys.path 등록. 노트북은 thin runner(import+호출)로 동일 산출물 재현
./Modeling/engine_state/ : run_model_engine.py 산출 (학습된 engine) = config.ENGINE_DIR. genes/scaler/config.pkl · rare_glm.pkl(pooled rare GLM 계수) · training_summary.csv(rare 포함)
./Modeling/Z_scores/ : normative model Z-score 산출물 전용 (Z_disease/sample/gene/hc/hc_names.npy = engine-only canonical · Z_rare_disease/hc/gene_names.npy = rare 공변량 GLM 별도 아티팩트 · disease_scores_flagged.parquet = 전 분기 통합 long표) = config.Z_SCORES_DIR
./Modeling/CV_Results/ : CV·평가 출력 (cv_*_stats.csv · cv_*_zscores.pkl · Figures/). Z 행렬은 여기 아님(→Z_scores/)
./Modeling/GSEA/ : GSEA 산출 (조건별 하위폴더 no_filter/with_rare/ubiqNN_absN, 각 gsea_result_*.csv · Clusters/ · Figures/) + 해석 리포트. no_filter/GSEA_Master_Report.md(rare 미포함) · with_rare/GSEA_Master_Report.md(rare 포함, DB+20질병 문헌검증: rare-led 신규신호 카탈로그+미보고 후보) · Analysis_Provenance.md(rare/DESeq2 비교에 쓴 모든 DB 엔드포인트·쿼리·resolve된 질병ID·PubMed 쿼리/히트 기록, 검증·재현용)
./Modeling/Benchmark/ : Normative Modeling vs DESeq2 정성/정량 비교 전용 = config.BENCHMARK_DIR. deseq2_results/(deseq2_*.csv·gsea_result_*.csv per phenotype, PyDESeq2 within-study 공변량미보정=config.DESEQ2_RESULTS_DIR) · deseq2_gsea/(공변량미보정 GSEA=config.DESEQ2_GSEA_DIR) · deseq2_covariate_results/(공변량보정 DESeq2 결과=config.DESEQ2_COV_RESULTS_DIR, run_deseq2_covariate.py 산출) · deseq2_covariate_gsea/(공변량보정 GSEA=config.DESEQ2_COV_GSEA_DIR) · disease_reference/(질병별 Open Targets association 상위 300 유전자 JSON = DB 교차검증 참조) · gsea_compare/(gsea_compare.py 산출: overlap_stats.csv · rare_novel_validated/summary.csv · deseq2_coverage.csv · deseq2_cov_vs_nocov_overlap.csv · deseq2_cov_db_hits.csv · db_hit_rates*.csv · {comparison}__{which}__{pheno}.csv diff 리스트) · rescued_genes_*.csv(분석1 산출) · DESeq2_vs_Normative_Report.md(term 커버리지/방향불일치 비교) · Figures/
./Modeling/노트북 : disease_scoring → gene_selection → gene_enrichment → gsea_heuristic_signatures(수동 theme,PPT용) · gsea_enrichmentmap_signatures(데이터기반 군집) · model_comparison(CV 평가) · deseq2_benchmark(분석1: DESeq2 independent-filtering 제외 유전자 중 normative Z-score로 살아나는 시그널 탐색, phenotype 단위로 점진 확장 예정) · gsea_rare_deseq2_comparison(rare 포함/미포함/DESeq2 3자 GSEA term 비교 재현 thin-runner: gsea_compare 호출→gsea_compare/ 산출물 재생성)

### pipeline/ 주요 진입점 (노트북에서 실제 호출되는 함수)
- `data_prep.load_disease_filtered()` → `DiseaseData` dataclass (Z_dis · dis_pheno · dis_names · gene_names · gene_syms · adata · is_hc · X_raw) 반환. 분석 노트북의 공통 진입점.
- `data_prep.ood_min_samples_filter()` → OOD + min_samples 양방향 필터 적용 후 (Z, pheno, names, ood, keep, excluded) 반환.
- `scoring.load_engine()` → engine_state/가 있으면 NormativeModelEngine.load(), 없으면 학습 후 저장.
- `disease_scores_flagged.parquet` (Z_scores/) → 샘플×유전자 Z-score를 z_flag(3.0) 기준으로 이진화한 플래그 표.

### 핵심 아키텍처 (여러 파일을 읽어야 파악되는 큰 그림)
- **Normative Model = 분기 모델**. `NormativeModelEngine`(model_engine.py)이 protein-coding 유전자를 HC detection rate 기준으로 분기 학습 (경계: rare<`rare_det_max`=1% · logistic 1~`low_det_thr`=10% · NBI/ZINBI≥10%):
  - 극저발현(`rare_det_max` 미만, 559개) → **rare** 분기. 모든 rare 유전자를 풀링한 단일 공변량 GLM(offset=log(mean_hc_j+eps), shared beta). Poisson 우선, pooled deviance/df가 `rare_overdisp_thr`(2.0, 관대) 초과 시에만 NB. Poisson/NB RQR로 z-score (RARE_Z_CAP로 클립). 학습/스코어 모두 순수 파이썬(statsmodels), R 불필요.
  - 저발현(`low_det_thr` 미만) → **logistic** 분기. L2 LogisticRegression으로 P(detected|covariate) 추정, Bernoulli RQR로 z-score.
  - 고발현 → **NBI/ZINBI GAMLSS** 분기. R `gamlss`(gamlss.r)를 rpy2로 호출해 mu/sigma(/nu) 회귀계수 학습, full quantile residual로 z-score.
- **Z-score = randomized quantile residual (RQR)**. HC 규범 분포가 맞으면 z ~ N(0,1). 분기별 RQR 계산은 model_engine.py의 순수 파이썬 함수(`_bernoulli_rqr`/`_poisson_rqr`/`_nb_rqr`/`_nbi_rqr_from_coeffs`/`_zinbi_rqr_from_coeffs`)에 있어 scoring 시 R 불필요(NBI/ZINBI 학습 시에만 R 필요).
- **공변량(X) = BIAS_COLUMNS 10개**(config). HC로 fit한 StandardScaler로 표준화 후 모델 입력. disease 샘플은 학습된 scaler/계수로 score만.
- **데이터 흐름**: h5ad(config.H5AD_PATH=merged_qc) → 엔진 학습(engine_state/) → disease scoring(Z_scores/ 의 Z_disease.npy=engine-only canonical + Z_rare_disease.npy=rare 별도) → OOD(Mahalanobis, HC-fit) + MIN_SAMPLES 필터(data_prep.py) → gene_selection/enrichment/GSEA/comparison.
- **rare 저장/사용은 통합 + opt-in**. canonical Z_disease.npy는 engine-only(rare 컬럼 0 placeholder) 유지로 기존 downstream 불변. rare 공변량 GLM z는 Z_rare_*.npy로 따로 저장되고, disease_scores_flagged.parquet에는 전 분기 통합 long으로 들어감. GSEA/signature 등에서 rare를 합쳐 쓸지는 `scoring.score_disease_with_rare(dd)` / `scoring.load_z(with_rare=True)`로 분석 단계에서 선택(기본은 미포함).
- **pipeline/ 패키지가 분석 로직, 노트북은 thin runner**. 같은 분석을 재현하려면 노트북이 아니라 pipeline 모듈을 수정. `pipeline/__init__.py`와 각 엔트리 스크립트가 자체적으로 root를 sys.path에 등록하므로 config/모듈 import는 재선언 없이 동작.

### 실행 명령어
- **R 의존성**: 엔진 학습(run_model_engine.py)과 `cv_gamlss_nb/zinb.py`는 R + `gamlss` 패키지가 설치되고 rpy2가 인식 가능해야 함(gamlss.r를 source). logistic 분기·모든 scoring(RQR)은 순수 파이썬으로 R 불필요.
- 엔진 학습 → engine_state/: `python Modeling/run_model_engine.py`
- CV 벤치마크 → CV_Results/(cv_*_stats.csv · cv_*_zscores.pkl): `python Modeling/cv_gamlss_nb.py` · `python Modeling/cv_gamlss_zinb.py` · `python Modeling/cv_logistic.py` · `python Modeling/cv_rare.py`(rare pooled GLM, 순수 파이썬)
- rare 분기만 재학습(전체 NBI 재학습 회피) → `python Modeling/run_model_engine.py --rare-only`
- 분석 노트북 실행 순서: disease_scoring → gene_selection → gene_enrichment → gsea_heuristic_signatures / gsea_enrichmentmap_signatures, 그리고 model_comparison(CV 평가). EDA는 cwd=EDA 가정.
- 테스트 스위트·린터·빌드 시스템 없음(연구 코드). 검증은 노트북 재실행/스크립트 산출물 확인으로 수행.

### 신규 분석 노트북 추가 시 체크리스트
1. `pipeline/` 모듈에 로직 구현 → 노트북은 import+호출만 (thin runner 원칙)
2. `config.py`의 경로/파라미터 import, 재선언 금지
3. `apply_style()` 호출 확인

### 데이터베이스 참조/논문참조
skill 중 /paper-lookup, /database-lookup, /scientific-critical-thinking 을 효율적이게 활용하여, 사용자가 결과의 해석을 요청한 경우 반드시 fetching과 skill을 적절히 활용하여 기존 연구결과의 엄격한 검증을 통해 해석을 수행할 것. 반드시 과학적인 근거가 있는 내용만을 보수적으로 제공할 것.
