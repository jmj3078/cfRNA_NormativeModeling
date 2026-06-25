import sys
from pathlib import Path

_MODELING = Path(__file__).resolve().parent.parent
_ROOT = _MODELING.parent
for _p in (str(_ROOT), str(_MODELING)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
