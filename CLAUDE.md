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

./config.py : 경로·파라미터 단일 소스. 모든 코드가 재선언 없이 import. ROOT/DATA_DIR/PIPELINE_DIR · MODELING_DIR/ENGINE_DIR/CV_RESULTS_DIR/Z_SCORES_DIR/GSEA_DIR · H5AD_PATH · Z_DISEASE/Z_HC 등 Z 경로 · RARE_REF · PATHS{merged_raw/biases/qc} · BIAS_COLUMNS · MODELING_PARAMS{ood_percentile/min_samples/z_flag/det_rate_min=0.01/lr_c=1.0/gsea_*/emap_sim_thr ...} · PARAMS(EDA용)
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
./Modeling/*.py : model_engine · run_model_engine(→engine_state/) · sample_filter(MahalanobisFilter OOD) · gene_selectors · cv_gamlss_nb/zinb · cv_logistic · rare_event_scorer. 경로/BIAS_COLUMNS/임계값 전부 config import
./Modeling/pipeline/ : 분석·시각화 모듈 패키지. data_prep(공통 전처리:load_adata/study-split/OOD·MIN_SAMPLES/Z 로드)·scoring·selection·enrichment·signatures(THEMES+heuristic/emap 군집)·comparison + plots(시각화 전용). __init__이 root를 sys.path 등록. 노트북은 thin runner(import+호출)로 동일 산출물 재현
./Modeling/engine_state/ : run_model_engine.py 산출 (학습된 engine) = config.ENGINE_DIR
./Modeling/Z_scores/ : normative model Z-score 산출물 전용 (Z_disease/sample/gene/hc/hc_names.npy · disease_scores_flagged.parquet · rare_event_ref.pkl) = config.Z_SCORES_DIR
./Modeling/CV_Results/ : CV·평가 출력 (cv_*_stats.csv · cv_*_zscores.pkl · Figures/). Z 행렬은 여기 아님(→Z_scores/)
./Modeling/GSEA/ : GSEA 산출 (gsea_result_*.csv · Clusters/ · Figures/) + 해석 리포트 GSEA_Master_Report.md · GSEA_Master_Report_Rework.md
./Modeling/노트북 : disease_scoring → gene_selection → gene_enrichment → gsea_heuristic_signatures(수동 theme,PPT용) · gsea_enrichmentmap_signatures(데이터기반 군집) · model_comparison(CV 평가)



## 데이터베이스 참조/논문참조
skill 중 /paper-lookup, /database-lookup, /scientific-critical-thinking 을 효율적이게 활용하여, 사용자가 결과의 해석을 요청한 경우 반드시 fetching과 skill을 적절히 활용하여 기존 연구결과의 엄격한 검증을 통해 해석을 수행할 것. 반드시 과학적인 근거가 있는 내용만을 보수적으로 제공할 것.

