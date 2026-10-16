# Final Modeling Method Improvement -- 계획 (초안, 2026-10-16)

목표: 같은 분할, 같은 in-silico 데이터, 같은 지표 코드로 후보 파이프라인을 비교한다.
우선순위: LOBO > signal detection > FDR 방어.

## 후보

| id | 공변량 | RQR | 비고 |
|---|---|---|---|
| `p1_baseline` | BIAS_COLUMNS 10개 | random effect 주변화 | 현행 그대로 |
| `p2_blup` | BIAS_COLUMNS 10개 | 배치 BLUP 조건부 | Z 클리핑 범위는 p1과 동일. held-out 배치의 BLUP은 그 배치 HC로 추정(누수 감수) |
| `p3_<tag>` | RUV 요인 W(k), 단독 또는 기술 공변량과 조합 | 주변화 | background gene 선정에 질병 샘플 사용(누수 감수), W 적합은 학습 HC만 |

## 디렉토리

```
_FinalModelingMethodImprovement/
  config.py           자체 경로·파라미터·PIPELINES 정의. 루트 config / MixedEffectsModeling import 금지
  engine/             MixedEffectsModeling/core + validation/{cv,lobo}_engine 포크
                      바뀌는 곳은 두 군데: covariate provider, rqr mode(marginal|blup)
  covariates/         provider: fit(Y_hc_train, obs_train) -> state, transform(Y, obs) -> X
  eval/
    splits.py         CV fold, LOBO fold를 한 번 고정해 splits/*.json에 저장(다시 쓰지 않음)
    insilico.py       sweep.py(nb5)·nb3 설계 그대로, 주입 세트를 한 번 생성해 insilico/에 저장
    metrics_cv.py     moments, PPC(mean/sd/zero: Pearson/Spearman/RMSE/MAE), exceedance, outlier
    metrics_lobo.py   calibration + MMD(distance/AUC/n_sig, _legacy/_lobo_mmd 이식)
    metrics_insilico.py  TP/FP/FN, PR, F1, detection limit(FDR 방어)
  baseline_import.py  기존 산출물을 읽어 baseline 행 생성 + 재현 검증
  runs/<pipeline>/{cv,lobo,full}/   fit, Z, 실행별 tmp(공유 /tmp 사용 금지)
  results/metrics_{cv,lobo,insilico}.csv   long format: pipeline, split, unit, metric, value
  0_covariate_eda.ipynb   사용자가 background gene 집합과 k, 변수 선택을 결정
  1_compare.ipynb         results/*.csv만 읽어 시각화
```

## 단계와 게이트

| 단계 | 내용 | 통과 조건 |
|---|---|---|
| S0 | 분할 고정, in-silico 세트 생성 (Ward 데이터 확정 후) | 모든 split에서 train∩test=∅, held-out 배치가 train 배치 목록에 없음 |
| S1 | baseline 행 생성 | (a) 새 eval 코드가 기존 csv(nb2 `lobo_hc_calibration_summary`, `sweep_nm_counts`, nb1 `cv_shash_calibration`)를 정확히 재현, (b) 포크 엔진이 기존 fit으로 Z를 bit 단위로 재생산 |
| S2 | p2 구현: R fit에서 배치 BLUP 출력, held-out 배치 BLUP 추정 | gene 50개 smoke run, BLUP=0으로 두면 p1 Z와 일치 |
| S3 | 공변량 EDA 노트북 (사용자 결정) | -- |
| S4 | p3 provider 구현 -> smoke -> 전체 CV/LOBO/full 적합 -> 행 추가 | provider가 학습 HC 외 샘플을 보지 않음(assert) |

누수 방지(코드 assert): 분할 파일 해시를 모든 run 메타에 기록, 캐시 키 = 설정 해시(존재 여부 아님), run마다 tmp 디렉토리 분리.

## 기존 파이프라인에서 발견된 문제 (S1은 그대로 재현하고, 수정은 별도 변형으로)

- CV/sweep가 fold별 trend가 아니라 배포 엔진 trend(`DISPERSION_TREND_PATH`)를 씀
- `draw_perturbed`가 CV fold i와 i번째 LOBO 배치에 같은 유전자를 주입
- nb2만 OOD 필터 적용, nb3/sweep는 미적용
- nb3과 sweep.py의 주입 설계(참조 mu, 테스트 집합, BH universe)가 달라 서로 비교 불가
- nb1의 `obs_var`/`pred_var` 열은 실제로 SD

## 결정 사항 (2026-10-16)

1. in-silico: 기존 설계를 그대로 쓴다. 방법 비교 = nb5/`MethodComparison/sweep.py`, detection limit = nb3. 두 설계는 서로 섞지 않는다
2. LOBO OOD 필터: 지표마다 현재 방식 그대로 (nb2 적용, nb3/sweep 미적용)
3. 디렉토리: `_FinalModelingMethodImprovement`
4. 브랜치: `method-improvement` (`ward-lane-merge`에서 분기)
5. MMD: `_legacy/_lobo_mmd/lobo_mmd.py`를 가져와 쓴다

## 진행 상황

- [x] 엔진 포크 (`engine/`), 자체 `config.py` (값은 루트 config와 동일함을 확인), 공유 /tmp 제거, CV fold는 `splits/cv_folds.json`과 일치해야만 진행
- [x] `check_fork_repro.py`: 포크 엔진이 production engine state로 in-sample HC Z를 bit 단위 재생산 (679x19,851 전부 일치, 2026-10-16 pre-Ward state)
- [x] `eval/lobo_mmd.py`: `_legacy/_lobo_mmd`에서 이식 (import만 변경)
- [ ] S0 `eval/splits.py`, in-silico 세트
