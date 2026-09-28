"""Re-export of the project-root config.py.

Loaded by path rather than `import config` so a script run with MixedEffectsModeling/ as its own
sys.path[0] cannot resolve `config` back to this shim.
"""
import importlib.util
from pathlib import Path

_path = Path(__file__).resolve().parent.parent / "config.py"
_spec = importlib.util.spec_from_file_location("_cfrna_root_config", _path)
_root = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_root)
globals().update({k: v for k, v in vars(_root).items() if not k.startswith("_")})
