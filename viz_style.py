import matplotlib
import matplotlib.pyplot as plt

def apply_style():
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
