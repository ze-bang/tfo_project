# Publication style: true Computer Modern via usetex, thin spines, tight ticks.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def apply(fontsize=9):
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "text.latex.preamble": r"\usepackage{amsmath}",
        "font.size": fontsize,
        "axes.titlesize": fontsize,
        "axes.labelsize": fontsize,
        "xtick.labelsize": fontsize - 1,
        "ytick.labelsize": fontsize - 1,
        "legend.fontsize": fontsize - 1,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.edgecolor": "0.15",
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "figure.dpi": 200,
    })

MODES = [(r"$q_{\rm FM}$", 0.38), (r"$E_{12}$", 0.50), (r"$E_{23}$", 0.70),
         (r"$q_{\rm AFM}$", 0.90), (r"$E_{13}$", 1.20)]

def mode_guides(ax, xmax=1.85, ymax=1.4, labels_top=True, color="w", alpha=0.30):
    for nm, v in MODES:
        ax.axvline(v, color=color, ls=(0, (4, 3)), lw=0.5, alpha=alpha)
        ax.axhline(v, color=color, ls=(0, (4, 3)), lw=0.45, alpha=alpha * 0.85)
        ax.axhline(-v, color=color, ls=(0, (4, 3)), lw=0.45, alpha=alpha * 0.85)
        if labels_top:
            ax.text(v, ymax * 0.985, nm, ha="center", va="top",
                    color=color, alpha=0.75, fontsize=6.2)

def panel_tag(ax, s, color="k"):
    ax.text(0.012, 0.985, s, transform=ax.transAxes, ha="left", va="top",
            fontsize=10, fontweight="bold", color=color)
