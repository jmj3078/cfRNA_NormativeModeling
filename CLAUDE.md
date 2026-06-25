# cfRNA Normative Modeling — 공통 프로젝트 요구사항

## 코드 작성 요령
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

## 프로젝트 디렉토리 구조
./_legacy/ : 참고, 수정, 읽기 금지, 현재 폐기된 프로젝트 파일 및 임시파일 저장소
./Data/ : 권한 없이 수정 금지, 분석 핵심 데이터 저장공간
./OpenAccess_nfcore/ : 권한 없이 수정 금지, 주요 분석이 이뤄지는 데이터(adata)의 원본 및 전처리 이전 버전의 데이터 존재.
./RPM_nfcore/ : Validation 을 위한 추가 실험실 데이터 추가 예정
./EDA/ : 코호트 QC·batch/bias 교란 분석 파이프라인 + 결과물(Analysis_Results). 루트의 config/viz_style를 import (노트북 실행 cwd=EDA 기준)
./EDA/Analysis_Results/ : analysis_cfrna_cohorts.ipynb 출력 (노트북에서 ./Analysis_Results/ 상대경로)
./Saved_Pipeline/ : config.PIPELINE_DIR. LogisticGP 파이프라인·Z matrix 저장(생성 예정)
./Modeling/ : cfRNA Normative Modeling 전반에 대한 정보 포함
./Modeling/engine_state/ : run_model_engine.py 실행 결과 저장
./Modeling/CV_Results/ : cv_gamlss_nb, cv_gamlss_zinb, cv_logistic, disease_scoring, gene_selection, rare_event_scorer 등의 스크립트 및 노트북의 출력 저장 경로
./Modeling/GSEA/ : cfRNA Normative Modeling 결과를 바탕으로 Z-score 산출 -> GSEA 분석 수행결과 저장

### 루트 레벨 (모든 프로젝트 공용 참조만 유지)
config.py : 경로/파라미터 단일 관리 (아래 별도 항목)
viz_style.py : apply_style() 공통 matplotlib 테마
CLAUDE.md · README.md

### EDA/ (코호트 QC·교란 분석 파이프라인)
analysis_cfrna_cohorts.ipynb : 코호트 QC·batch/bias 교란 분석 (helper+plot 사용)
analysis_helper.py : QC·bias 정량화 + 분산분해(RDA) 엔진
analysis_plot.py : analysis_helper 전용 시각화 함수 모음
VariousNormalizationMethods_OpenAccess.R : 정규화 레이어 생성(R, PROJECT_ROOT 절대경로 기준)
- 경로 규약: 노트북은 cwd=EDA 가정. 루트 모듈(config/viz_style)은 `Path.cwd().parent`를 sys.path에 추가해 import. analysis_helper.py는 자체적으로 루트를 sys.path에 등록(config import용).

## 경로/파라미터 관리 (config.py)
프로젝트 내 중복 경로 선언을 통합 관리하는 단일 소스. **신규 스크립트/노트북은 경로·상수를 재선언하지 말고 config에서 import 할 것.**
- `from config import PATHS, PARAMS`
- `ROOT`(프로젝트 루트) · `DATA_DIR`(OpenAccess_nfcore) · `PIPELINE_DIR`(Saved_Pipeline)
- `PATHS` : merged_raw / merged_biases / merged_qc(주 h5ad) · pipeline_* · z_total/state/quant
- `PARAMS` : min_study_samples, n_top_genes, n_pcs, loess_frac, n_bins, outlier_pct, min_expressed
- 현재 analysis_helper.py·analysis_cfrna_cohorts.ipynb만 사용. Modeling/ 노트북은 자체 BASE_DIR 선언 중 → 점진적으로 config 통합 권장.

### Modeling/ (Normative modeling 본체)
- 엔진: model_engine.py · run_model_engine.py(→engine_state/) · sample_filter.py(MahalanobisFilter OOD) · gene_selectors.py
- CV 스크립트: cv_gamlss_nb.py · cv_gamlss_zinb.py · cv_logistic.py · rare_event_scorer.py (→CV_Results/)
- 노트북 흐름: disease_scoring(Z 산출) → gene_selection(gene set+CV) → gene_enrichment(GSEA prerank→GSEA/*.csv) → gsea_heuristic_signatures(수동 theme 시각화, PPT용) · gsea_enrichmentmap_signatures(leading-edge 데이터기반 군집 검증) · model_comparison(CV 평가 요약)
- GSEA 해석 산출물: GSEA/GSEA_Master_Report.md · GSEA_Master_Report_Rework.md

## 데이터베이스 참조/논문참조
skill 중 /paper-lookup, /database-lookup, /scientific-critical-thinking 을 효율적이게 활용하여, 사용자가 결과의 해석을 요청한 경우 반드시 fetching과 skill을 적절히 활용하여 기존 연구결과의 엄격한 검증을 통해 해석을 수행할 것. 반드시 과학적인 근거가 있는 내용만을 보수적으로 제공할 것.

