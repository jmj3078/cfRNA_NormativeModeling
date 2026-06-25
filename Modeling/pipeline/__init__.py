import sys
from pathlib import Path

# project root(config, viz_style) + Modeling(model_engine, sample_filter, gene_selectors)을
# import 경로에 등록 → 노트북 cwd와 무관하게 동작
_MODELING = Path(__file__).resolve().parent.parent
_ROOT = _MODELING.parent
for _p in (str(_ROOT), str(_MODELING)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
