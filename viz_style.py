import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Arial has no Hangul. Register the system Noto CJK collection once and let
# matplotlib's per-glyph fallback pick it up, so Latin stays Arial and Korean
# labels render instead of tofu. Silently skipped if the font is not installed.
_CJK_TTC = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
_CJK_NAME = None
try:
    font_manager.fontManager.addfont(_CJK_TTC)
    _CJK_NAME = next(f.name for f in font_manager.fontManager.ttflist if "CJK" in f.name)
except Exception:
    pass

# Named sizes for manual ax.text()/set_title(fontsize=...) annotations that
# rcParams can't reach. Two tiers only: pick the closer one, don't add a third.
FS_ANNOT_LG = 12
FS_ANNOT_SM = 8

def apply_style():
    family = ["Arial"] + ([_CJK_NAME] if _CJK_NAME else [])
    matplotlib.rcParams["font.family"] = family
    plt.rcParams.update({
        'font.family':        family,
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
