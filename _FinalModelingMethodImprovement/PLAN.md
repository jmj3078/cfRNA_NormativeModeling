# Final Modeling Method Improvement -- 계획 (2026-10-16 갱신)

목표: 같은 분할, 같은 in-silico 데이터, 같은 지표 코드로 후보 파이프라인을 비교한다.
우선순위: LOBO > signal detection > FDR 방어.

## 후보 x 시나리오 (2026-10-16 확정)

| id | 내용 | zero_shot | adaptation |
|---|---|---|---|
| `p0` | 현행 파이프라인 (BIAS 10개, 주변화 RQR) | O | O (평균 shift 재보정) |
| `p2` | background gene RUV 공변량 | O | O (평균 shift 재보정) |
| 이상치 방법 4개 | OUTRIDER(AE), PCA, OutSingle, PEER, 재적합 없음 | O | O (평균 shift 재보정) |

비교 조합: p0_zeroshot, p0_adaptation, p2_zeroshot, p2_adaptation (+ 이상치 방법 x 2).

adaptation = held-out 배치 HC를 2-fold 교차적합으로 사용, 자기 자신은 보정에 쓰지 않음.
재보정 방식은 평균만 (평균+SD, 중앙값, 중앙값+MAD는 2026-10-16 비교 후 기각, `eval/adaptation.py` docstring).
CV는 adaptation 없음.

p1(배치 BLUP 조건부 RQR)은 2026-10-16 제외: "새 배치 HC로 위치만 보정하면 나아지는가"는 adaptation이
재적합 없이 이미 답하며(평균 shift로 엔진 null FP 42 -> 37, 이상치 방법도 거의 불변), BLUP은 SHASH/pool/rare
경로와 검증 게이트를 새로 요구한다. 한계: Z 척도 shift는 count 척도 배치 절편과 동일하지 않음(저발현 RQR 비선형).
리뷰 지적 시 대응안 = oracle BLUP 사전검정(held-out 배치 비교란 HC 전체로 b-hat 추정, plug-in, in-silico LOBO).

## 디렉토리

```
_FinalModelingMethodImprovement/
  config.py             자체 경로·파라미터·PIPELINES, FMMI_PIPELINE로 run 선택. 루트 config / MixedEffectsModeling import 금지
  engine/               MixedEffectsModeling/core + validation/{cv,lobo}_engine 포크 (바뀌는 곳: covariate provider)
  check_fork_repro.py   포크 엔진 in-sample Z bit 재현 게이트
  repro_scoring.py      production fit을 runs/p0에 staging 후 CV/LOBO Z bit 재현 게이트
  eval/
    splits.py           CV fold, LOBO fold 고정 (splits/*.json)
    insilico_sweep.py   nb5 sweep 설계: freeze(주입 세트 1회 생성) / score / import(이상치 방법)
    adaptation.py       교차적합 평균 shift (in-silico 전 방법)
    detection_limit.py  nb3 설계, --adapt
    score_full.py       full fit in-sample HC/disease Z
    calibration.py      CV/LOBO HC calibration 표 (엔진 파이프라인만)
    lobo_mmd.py         _legacy/_lobo_mmd 이식 (아직 미통합)
    build_results.sh <pipe>   위 단계를 순서대로 실행, 로그는 runs/<pipe>/log_*.txt
  runs/<pipeline>/{full,cv,lobo,insilico,detection_limit,tmp}/
  results/*.csv         method x scenario 행 추가 방식 (eval/results.append)
  0_covariate_eda.ipynb   사용자가 background gene 집합과 k, 변수 선택을 결정
  1_compare.ipynb         in-silico / detection limit: 엔진 + 이상치 방법
  2_calibration.ipynb     CV/LOBO calibration: 엔진 파이프라인만 (KDE, violin)
```

## 단계와 게이트

| 단계 | 내용 | 통과 조건 |
|---|---|---|
| S0 | 분할 고정, in-silico 세트 생성 | production 분할·주입 행렬과 정확히 일치 |
| S1 | p0 baseline 행 생성 | 포크 엔진·eval 코드가 production Z와 csv를 재현 |
| S2 | 공변량 EDA 노트북 (사용자 결정) | -- |
| S3 | p2 provider 구현 -> smoke -> 전체 CV/LOBO/full 적합 -> 행 추가 | provider가 학습 HC 외 샘플을 보지 않음(assert) |
| S4 | LOBO MMD 통합 | p0가 production MMD 결과 재현 |
| S5 | Ward lane-merge 재적합 후 전체 재실행 | check_fork_repro 재실행 |

누수 방지(코드 assert): 고정 분할과 일치해야만 진행, run마다 tmp 디렉토리 분리, staging은 기존 run 덮어쓰기 거부.

## 기존 파이프라인에서 발견된 문제 (S1은 그대로 재현하고, 수정은 별도 변형으로)

- CV/sweep가 fold별 trend가 아니라 배포 엔진 trend(`DISPERSION_TREND_PATH`)를 씀
- `draw_perturbed`가 CV fold i와 i번째 LOBO 배치에 같은 유전자를 주입
- nb2만 OOD 필터 적용, nb3/sweep는 미적용
- nb3과 sweep.py의 주입 설계(참조 mu, 테스트 집합, BH universe)가 달라 서로 비교 불가
- nb1의 `obs_var`/`pred_var` 열은 실제로 SD
- cache csv를 pandas 기본 float parser로 읽으면 ulp 차이 발생 -> 포크 엔진은 `float_precision="round_trip"`

## 결정 사항 (2026-10-16)

1. in-silico: 기존 설계를 그대로 쓴다. 방법 비교 = nb5/`MethodComparison/sweep.py`, detection limit = nb3. 두 설계는 서로 섞지 않는다
2. LOBO OOD 필터: 지표마다 현재 방식 그대로 (nb2 적용, nb3/sweep 미적용)
3. 디렉토리: `_FinalModelingMethodImprovement`, 브랜치: `method-improvement` (`ward-lane-merge`에서 분기)
4. MMD: `_legacy/_lobo_mmd/lobo_mmd.py`를 가져와 쓴다
5. adaptation LOBO가 새 기준, zero-shot 결과는 보존. 이상치 방법은 재적합 없이 adaptation만
6. adaptation 재보정은 평균만
7. p1(BLUP) 제외
8. calibration 비교는 엔진 파이프라인끼리만

## 진행 상황 (pre-Ward production state 기준)

- [x] 엔진 포크, 자체 `config.py`, 공유 /tmp 제거
- [x] `check_fork_repro.py`: in-sample HC Z bit 재현 (679 x 19,851)
- [x] S0: 분할, in-silico 60셀 bit 일치
- [x] S1: CV Z/cv_stats, LOBO Z/SHASH Z/OOD mask(5배치), in-silico 엔진 Z·counts, detection limit(lobo 600행, insample_cv 240행), calibration(nb2 배치 요약, nb1 유전자 moment/outlier, PPC) 재현
- [x] p0 zero_shot / adaptation 행: in-silico, detection limit, calibration
- [x] `1_compare.ipynb`, `2_calibration.ipynb`
- [ ] S2 공변량 EDA (사용자)
- [ ] S3 p2
- [ ] S4 MMD 통합
- [ ] S5 Ward 재실행
