# cfRNA Normative Modeling — 공통 프로젝트 요구사항
### 분석 핵심 가정 및 목적
기존 cfRNA 전사체 분석은 주로 정상군과 질병군 간의 집단 수준 비교(Group-wise comparison)에 의존해 왔다. 그러나 생물학적 및 기술적 공변량(Covariates)에 의한 분산이 질병 고유의 신호를 압도하는 경우가 많다. 이로 인해 집단 단위의 일괄적인 공변량 보정은 질병 신호의 소실이나 교란 요인의 잔존을 초래하는 근본적인 한계가 존재한다. 이를 극복하기 위해, 본 연구는 전장 전사체(Whole-Transcriptome) 기반의 대규모 정상군(Healthy Control) 데이터를 활용한 규범적 모델링(Normative Modeling)을 도입한다. 개별 샘플의 공변량을 반영하여 정상 상태의 예상 분포를 추정하고, 이를 통계적 편차(Z-score)로 산출함으로써 교란 요인의 영향 없이 질병 특이적 신호를 정밀하게 정량화.

### 코드 작성 요령
간결하고, 최소한의 코드를 작성하는 것을 지향
함수 선언 : 함수 및 메서드 선언 이후 입력 변수별 dtype설정 금지, 출력 값에 대한 dtype 제시 금지
띄어쓰기 : 줄 맞춤을 위한 인위적인 띄어쓰기 금지
이모지 : 이모지 절대 사용 금지
주석 : 영어로 작성하되, 사용자가 요청하기 이전에 주석 작성 금지 (최종 수정 후 일괄 추가 예정)
import : 알파벳 순으로 수행
데이터 시각화시 ./viz_style.py import 하여 일괄 테마 적용, 작업중인 notebook이나 스크립트의 경로를 바탕으로 다음과같이 import하여 적용할 것
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)
from viz_style import apply_style
apply_style() 
경로 및 재사용 가능성이 높은 관련 변수는 항상 ./config에서 전역 참조할 수 있도록 설정할것.

### 프로젝트 디렉토리 구조
./_legacy/ : 참고, 수정, 읽기 금지, 현재 폐기된 프로젝트 파일 및 임시파일 저장소
./Data/ : 권한 없이 수정 금지, 분석 핵심 데이터 저장공간
./OpenAccess_nfcore/ : 권한 없이 수정 금지, 주요 분석이 이뤄지는 데이터(adata)의 원본 및 전처리 이전 버전의 데이터 존재.
./RPM_nfcore/ : Validation 을 위한 추가 실험실 데이터 추가 예정
./EDA/ : 코호트 QC·batch/bias 교란 분석 파이프라인 + 결과물(Analysis_Results). 루트의 config/viz_style를 import (노트북 실행 cwd=EDA 기준)
./EDA/Analysis_Results/ : analysis_cfrna_cohorts.ipynb 출력 (노트북에서 ./Analysis_Results/ 상대경로)
./Saved_Pipeline/ : config.PIPELINE_DIR. LogisticGP 파이프라인·Z matrix 저장(생성 예정)
./Modeling/ : cfRNA Normative Modeling 전반에 대한 정보 포함
./Modeling/engine_state/ : run_model_engine.py 실행 결과 저장
./Modeling/CV_Results/ : cv_gamlss_nb/zinb/logistic, gene_selection 등의 CV·평가 출력(cv_*_stats.csv·cv_*_zscores.pkl·Figures). **Z-score 행렬은 여기 아님 → Z_scores/**
./Modeling/Z_scores/ : run_model_engine 기반 normative model Z-score 산출물 전용 (Z_disease/sample/gene/hc/hc_names.npy · disease_scores_flagged.parquet · rare_event_ref.pkl). config.Z_SCORES_DIR
./Modeling/GSEA/ : cfRNA Normative Modeling 결과를 바탕으로 Z-score 산출 -> GSEA 분석 수행결과 저장
./Modeling/pipeline : 엔진/스크립트: model_engine.py · run_model_engine.py(→engine_state/) · sample_filter.py(MahalanobisFilter OOD) · gene_selectors.py · cv_gamlss_nb/zinb.py · cv_logistic.py · rare_event_scorer.py (→CV_Results/). **경로/BIAS_COLUMNS/임계값은 모두 root config.py에서 import** (DET_RATE_MIN 0.01·LR_C 1.0 통일). 
- **`Modeling/pipeline/` 패키지 (분석/시각화 모듈화)**: `data_prep`(공통 전처리: load_adata/study-split/OOD·MIN_SAMPLES/Z 로드) · `scoring`(disease_scoring) · `selection`(gene_selection) · `enrichment`(GSEA prerank) · `signatures`(THEMES+heuristic / emap 군집) · `comparison`(CV 평가) · **`plots`(모든 시각화 전용)**. `__init__`이 root를 sys.path 등록 → config/viz_style import.
- 노트북(제자리, **thin runner**: 모듈 import+호출만): disease_scoring → gene_selection → gene_enrichment → gsea_heuristic_signatures(수동 theme, PPT용) · gsea_enrichmentmap_signatures(데이터기반 군집) · model_comparison. 실행 시 동일 산출물 재현.
- 경로/파라미터 변경은 **config.py 한 곳**(MODELING_DIR/ENGINE_DIR/CV_RESULTS_DIR/**Z_SCORES_DIR**/GSEA_DIR, Z_DISEASE/Z_HC 등 Z-score 경로, RARE_REF, BIAS_COLUMNS, MODELING_PARAMS)에서. Z-score 행렬·rare_event_ref는 Z_scores/.
- GSEA 해석 산출물: GSEA/GSEA_Master_Report.md · GSEA_Master_Report_Rework.md

### EDA/ (코호트 QC·교란 분석 파이프라인)
analysis_cfrna_cohorts.ipynb : 코호트 QC·batch/bias 교란 분석 (helper+plot 사용)
analysis_helper.py : QC·bias 정량화 + 분산분해(RDA) 엔진
analysis_plot.py : analysis_helper 전용 시각화 함수 모음
VariousNormalizationMethods_OpenAccess.R : 정규화 레이어 생성(R, PROJECT_ROOT 절대경로 기준)
- 경로 규약: 노트북은 cwd=EDA 가정. 루트 모듈(config/viz_style)은 `Path.cwd().parent`를 sys.path에 추가해 import. analysis_helper.py는 자체적으로 루트를 sys.path에 등록(config import용).

### 루트 레벨 (모든 프로젝트 공용 참조만 유지)
config.py : 경로/파라미터 단일 관리 (아래 별도 항목)
viz_style.py : apply_style() 공통 matplotlib 테마
CLAUDE.md · README.md

## 경로/파라미터 관리 (config.py)
프로젝트 내 중복 경로 선언을 통합 관리하는 단일 소스. **신규 스크립트/노트북은 경로·상수를 재선언하지 말고 config에서 import 할 것.**
- `from config import PATHS, PARAMS`
- `ROOT`(프로젝트 루트) · `DATA_DIR`(OpenAccess_nfcore) · `PIPELINE_DIR`(Saved_Pipeline)
- `PATHS` : merged_raw / merged_biases / merged_qc(주 h5ad) · pipeline_* · z_total/state/quant
- `PARAMS` : min_study_samples, n_top_genes, n_pcs, loess_frac, n_bins, outlier_pct, min_expressed
- 현재 analysis_helper.py·analysis_cfrna_cohorts.ipynb만 사용. Modeling/ 노트북은 자체 BASE_DIR 선언 중 → 점진적으로 config 통합 권장.
- config는 수정이 이뤄질 때 마다 CLAUDE.md 수정



## 데이터베이스 참조/논문참조
skill 중 /paper-lookup, /database-lookup, /scientific-critical-thinking 을 효율적이게 활용하여, 사용자가 결과의 해석을 요청한 경우 반드시 fetching과 skill을 적절히 활용하여 기존 연구결과의 엄격한 검증을 통해 해석을 수행할 것. 반드시 과학적인 근거가 있는 내용만을 보수적으로 제공할 것.

