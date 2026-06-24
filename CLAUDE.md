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
- **Disease samples**: 913 samples, 20 phenotypes
- **Branching**: det≥10% → NBI GAMLSS | 1-10% → Logistic | <1% → Rare Event Scorer
- **Working directory**: `Modeling/`
- **Key outputs**: `CV_Results/Z_disease.npy`, `CV_Results/Z_hc.npy`, `engine_state/`
