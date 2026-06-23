"""
Shared plotting style for the perturbation-evaluation reproducibility notebooks.

Usage in any notebook:
    from perturb_style import *
    # K562_COLOR, RPE1_COLOR, FIG_MAIN, etc. are then available.
    # apply_style() is called automatically on import.

For publication-ready figures, save as PDF:
    fig.savefig('figs/my_figure.pdf', bbox_inches='tight')

The colour and marker scheme matches the manuscript's Figure 1:
    K562 = blue, RPE1 = orange
    random = circle, stratified = square
"""

import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ── Colour palette ──────────────────────────────────────────────────────────
# Cell-type axis (matches Figure 1 convention)
K562_COLOR    = '#1f77b4'    # matplotlib C0 blue
RPE1_COLOR    = '#ff7f0e'    # matplotlib C1 orange
NEUTRAL_COLOR = '#888888'    # grey for reference/background

CELL_TYPE_PALETTE = {
    'K562': K562_COLOR,
    'RPE1': RPE1_COLOR,
}

# Split-strategy axis (varied by marker shape rather than colour)
SPLIT_MARKER = {
    'random':     'o',
    'stratified': 's',
}

# Four condition combinations (cell type × split). Colour by cell type;
# distinguish split by marker shape, edge style, or alpha as needed per panel.
CONDITION_PALETTE = {
    'K562_random':     K562_COLOR,
    'K562_stratified': K562_COLOR,
    'RPE1_random':     RPE1_COLOR,
    'RPE1_stratified': RPE1_COLOR,
}

CONDITION_ORDER = ['K562_random', 'K562_stratified', 'RPE1_random', 'RPE1_stratified']

# ── Standard figure sizes (inches) ──────────────────────────────────────────
FIG_SINGLE = (3.5, 3.0)     # single panel
FIG_DOUBLE = (6.5, 3.0)     # double panel
FIG_MAIN   = (6.5, 2.5)     # Figure 1 (three panels in one row, PMLR single-column)
FIG_FULL   = (6.5, 5.0)     # full-page multi-panel
FIG_WIDE   = (10.0, 4.0)    # wide for presentations

# ── Output directories ──────────────────────────────────────────────────────
FIGS_DIR        = Path('figs')
DATA_DIR        = Path('data')
PRECOMPUTED_DIR = Path('precomputed')

# ── Style application ──────────────────────────────────────────────────────
def apply_style():
    """Apply the shared matplotlib/seaborn style. Called on import."""
    sns.set_style('whitegrid')
    plt.rcParams['figure.dpi']       = 100
    plt.rcParams['savefig.dpi']      = 300
    plt.rcParams['savefig.bbox']     = 'tight'
    plt.rcParams['font.size']        = 9
    plt.rcParams['axes.labelsize']   = 9
    plt.rcParams['axes.titlesize']   = 10
    plt.rcParams['xtick.labelsize']  = 8
    plt.rcParams['ytick.labelsize']  = 8
    plt.rcParams['legend.fontsize']  = 8
    plt.rcParams['figure.titlesize'] = 11

apply_style()
