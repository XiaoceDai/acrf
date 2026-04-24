"""
Cell / Nature Medicine inspired publication style — single source of truth.

Design rules (locked):

  * One semantic axis = one color.  Maximum 4 colour classes per figure:
      LA / RA chamber  -> red / blue  (the only chromatic axis)
      structure        -> dark / mid / light grey (text, lines, fills)
      highlight        -> single deep orange, used sparingly (≤ 2 per figure)
      support          -> single green, used only for stability / agreement

  * Panel letters: lowercase bold "a, b, c…" (Nature convention).
  * Typography hierarchy is fixed (see FS_*).  6 pt is the absolute floor;
    nothing smaller is allowed.
  * One colour bar per figure-level concept; never one per panel.
  * No decorative borders, no shadow, no 3D, no gradients except CMAP_*.

Usage:
    from plot_style import *
    set_nature_style()
    panel_letter(ax, "a")
    panel_title(ax, "Healthy LA / RA discriminates donors not chambers")
"""

from __future__ import annotations

from typing import Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ----------------------------------------------------------------------
# Page widths (Nature Medicine column geometry)
# ----------------------------------------------------------------------
MM = 1.0 / 25.4
COL_SINGLE  = 89  * MM   # 3.5"
COL_ONEHALF = 120 * MM   # 4.7"
COL_DOUBLE  = 183 * MM   # 7.2"
PAGE_HEIGHT = 240 * MM   # max usable single-page height


# ----------------------------------------------------------------------
# Final composite cell sizes (locked).
# Each panel is rendered at exactly this size (mm), saved with
# `bbox_inches=None`, and composed into Figure 1 with NO further scaling.
# ----------------------------------------------------------------------
def _in(mm: float) -> float:
    return mm * MM

# (width_mm, height_mm) for every Figure 1 panel
PANEL_A_MM = ( 60, 48)   # PCA
PANEL_B_MM = (115, 48)   # Volcano
PANEL_C_MM = (175, 62)   # NCP selection (full width)
PANEL_D_MM = ( 87, 62)   # aCAI donor summary (half)
PANEL_E_MM = ( 87, 62)   # axis families (half)
PANEL_F_MM = (175, 70)   # compartment attribution (full width)

PANEL_A_SIZE = (_in(PANEL_A_MM[0]), _in(PANEL_A_MM[1]))
PANEL_B_SIZE = (_in(PANEL_B_MM[0]), _in(PANEL_B_MM[1]))
PANEL_C_SIZE = (_in(PANEL_C_MM[0]), _in(PANEL_C_MM[1]))
PANEL_D_SIZE = (_in(PANEL_D_MM[0]), _in(PANEL_D_MM[1]))
PANEL_E_SIZE = (_in(PANEL_E_MM[0]), _in(PANEL_E_MM[1]))
PANEL_F_SIZE = (_in(PANEL_F_MM[0]), _in(PANEL_F_MM[1]))

# ----------------------------------------------------------------------
# Color tokens — locked semantic palette
# ----------------------------------------------------------------------

# Chromatic axis (the ONLY use of red / blue in any panel)
LA_COLOR = "#C97A74"   # muted warm red, LA
RA_COLOR = "#7D99B8"   # muted cool blue, RA

# Structural greys — for text, lines, background shading
FG_COLOR    = "#1A1A1A"   # primary text
MUTED_COLOR = "#555555"   # secondary text / annotations
FAINT_COLOR = "#888888"   # axis ticks, light annotations
LINE_COLOR  = "#CCCCCC"   # gridlines, structure separators
BG_COLOR    = "#F5F5F5"   # neutral panel fill

# Single-purpose accents — used at most twice per figure
HIGHLIGHT_COLOR = "#C8A076"   # muted warm orange — main emphasis
SUPPORT_COLOR   = "#6E9B72"   # muted green — stability / agreement / "8/8" callouts

# ----------------------------------------------------------------------
# Diverging colormaps — each tied to a single semantic axis
# ----------------------------------------------------------------------

# (1) chamber-asymmetry (aCAI) — positive = LA = red, negative = RA = blue
CHAMBER_CMAP = LinearSegmentedColormap.from_list(
    "chamber_red_blue",
    ["#6F89A8", "#7D99B8", "#D8E2EC",
     "#F7F7F7",
     "#F0DFDA", "#C97A74", "#B87A73"],
)

# (2) one-sided intensity — for things like compartment z-score where
# direction has no chamber meaning.  White → mid → dark grey, with a
# subtle warm tint at high values to draw the eye.
INTENSITY_CMAP = LinearSegmentedColormap.from_list(
    "intensity_grey_warm",
    ["#FFFFFF", "#EDEDED", "#C8C8C8", "#8E8E8E",
     "#4D4D4D", "#1A1A1A"],
)

# Optional: bidirectional grey diverging for "depleted vs enriched" without
# adding a new colour family.  Negative end picks up SUPPORT green only as a
# single tonal hint, positive end is HIGHLIGHT orange.  Both kept low-saturation
# so they read as 'structure', not 'category'.
DIVERGE_GREY = LinearSegmentedColormap.from_list(
    "diverge_grey",
    ["#3F6512", "#A6CEE3" if False else "#9CB48F",
     "#E5E5E5", "#F7F7F7", "#F0E0CC", "#D8A56B", "#A65A1B"],
)

# Backwards-compatibility re-exports (older fig scripts import these names)
LA_FILL = LA_COLOR
RA_FILL = RA_COLOR
CHAMBER_CMAP_R = CHAMBER_CMAP   # legacy alias
REWIRING_CMAP = ["#542788", "#998EC3", "#D8DAEB",
                 "#F7F7F7",
                 "#FEE0B6", "#F1A340", "#B35806"]
REWIRING_CLASS_COLORS = {
    "Amplified":        "#B48A66",
    "Amplified-weak":   "#D5B08A",
    "Preserved":        "#F1E3CF",
    "Attenuated-weak":  "#DDD9E7",
    "Attenuated":       "#AAA2C0",
    "Reversed-weak":    "#AAA2C0",
    "Reversed":         "#7D7399",
    "Low-support":      "#BDBDBD",
    "NA":               "#FFFFFF",
}
AFCTRL_CMAP = [
    "#5F8F8A", "#78A5A0", "#A9C7C3", "#D9E9E6",
    "#F5F5F5",
    "#EEE4D6", "#DFC8AC", "#C2A080", "#9E7A5A",
]
COMPARTMENT_COLORS = {
    "atrialCM":                "#C9A978",
    "fibroblast":              "#8FAF7A",
    "endothelial_or_vascular": "#8D86B4",
    "immune":                  "#C7926F",
    "unresolved":              "#BDBDBD",
    "Unresolved":              "#BDBDBD",
}
DURATION_CMAP = [
    "#8F78A9", "#A896BC", "#C8BAD7", "#E5DFF0",
    "#F7F7F7",
    "#DEECDC", "#BDD6B9", "#92B892", "#6E9B72",
]

# ----------------------------------------------------------------------
# Typography (locked).  6 pt floor.  Nothing under 6 pt is allowed.
# ----------------------------------------------------------------------
FS_PANEL_LETTER = 9    # bold lowercase a, b, c…
FS_PANEL_TITLE  = 8    # bold panel title (one line ideally)
FS_AXIS_LABEL   = 7
FS_TICK         = 6.5
FS_ANNOT        = 6.5  # in-figure annotations
FS_LEGEND       = 6.5
FS_FOOTNOTE     = 6.5  # italic muted footnote, figure-level

# Line weights (locked)
LW_AXIS = 0.5
LW_THIN = 0.4
LW_BOLD = 0.8


# ----------------------------------------------------------------------
# Style entry point
# ----------------------------------------------------------------------

def set_nature_style() -> None:
    """Apply the unified Cell / Nature Medicine style to matplotlib."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [
            "Arial", "Helvetica", "Liberation Sans", "Nimbus Sans",
            "DejaVu Sans",
        ],
        "font.size":         FS_TICK,
        "axes.titlesize":    FS_PANEL_TITLE,
        "axes.labelsize":    FS_AXIS_LABEL,
        "xtick.labelsize":   FS_TICK,
        "ytick.labelsize":   FS_TICK,
        "legend.fontsize":   FS_LEGEND,
        "legend.title_fontsize": FS_LEGEND,
        "figure.titlesize":  FS_PANEL_TITLE + 1,

        "axes.linewidth":    LW_AXIS,
        "xtick.major.width": LW_AXIS,
        "ytick.major.width": LW_AXIS,
        "xtick.minor.width": LW_THIN,
        "ytick.minor.width": LW_THIN,
        "xtick.major.size":  2.4,
        "ytick.major.size":  2.4,
        "xtick.minor.size":  1.4,
        "ytick.minor.size":  1.4,
        "lines.linewidth":   LW_BOLD,
        "patch.linewidth":   LW_THIN,

        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.grid":          False,

        "axes.edgecolor":     FG_COLOR,
        "axes.labelcolor":    FG_COLOR,
        "xtick.color":        FG_COLOR,
        "ytick.color":        FG_COLOR,
        "text.color":         FG_COLOR,

        "figure.dpi":         150,
        "savefig.dpi":        300,
        "savefig.bbox":       "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype":       42,
        "ps.fonttype":        42,
        "svg.fonttype":       "none",

        "legend.frameon":     False,
        "legend.borderpad":   0.3,
        "legend.labelspacing":0.25,
        "legend.handletextpad":0.4,
        "legend.handlelength": 1.1,
        "legend.columnspacing":0.8,

        "axes.titlepad":  3,
        "axes.labelpad":  2.0,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
    })


# ----------------------------------------------------------------------
# Helper utilities used across all panels
# ----------------------------------------------------------------------

def panel_letter(ax: plt.Axes, text: str,
                 x: float = -0.10, y: float = 1.06,
                 size: float = FS_PANEL_LETTER, **_kw) -> None:
    """Place a lowercase bold panel letter (a, b, c, …) in axes coords."""
    ax.text(x, y, text, transform=ax.transAxes,
            fontsize=size, fontweight="bold",
            ha="left", va="bottom", color=FG_COLOR)


# Backwards-compatibility alias for older scripts that imported `panel_label`.
panel_label = panel_letter


def panel_title(ax: plt.Axes, text: str, *, pad: float = 4) -> None:
    """Bold left-aligned panel title using the locked typography."""
    ax.set_title(text, loc="left", pad=pad,
                 fontsize=FS_PANEL_TITLE, fontweight="bold",
                 color=FG_COLOR)


def figure_footnote(fig: plt.Figure, text: str,
                    y: float = 0.012,
                    color: str = MUTED_COLOR) -> None:
    """Italic muted footnote at the very bottom of the figure."""
    fig.text(0.012, y, text, fontsize=FS_FOOTNOTE,
             ha="left", va="bottom", color=color, style="italic")


def strip_spines(ax: plt.Axes,
                 keep: Iterable[str] = ("left", "bottom")) -> None:
    """Keep only the specified spines visible (default: left + bottom)."""
    keep_set = set(keep)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep_set)


def shared_chamber_colorbar(fig: plt.Figure, im, *,
                             rect: Sequence[float] = (0.34, 0.045, 0.32, 0.012),
                             label: str = "aCAI = log\u2082(LA) − log\u2082(RA)"
                             ) -> None:
    """Standardised horizontal LA/RA colour bar, one per figure."""
    cax = fig.add_axes(list(rect))
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label(label, fontsize=FS_AXIS_LABEL, color=FG_COLOR, labelpad=2)
    cb.outline.set_linewidth(LW_THIN)
    cb.ax.tick_params(width=LW_THIN, length=2, labelsize=FS_TICK,
                      color=FG_COLOR, labelcolor=FG_COLOR)
    cb.ax.text(0.0, 0.5, "RA ←", transform=cb.ax.transAxes,
               ha="right", va="center", fontsize=FS_LEGEND,
               color=RA_COLOR, fontweight="bold")
    cb.ax.text(1.0, 0.5, "→ LA", transform=cb.ax.transAxes,
               ha="left", va="center", fontsize=FS_LEGEND,
               color=LA_COLOR, fontweight="bold")
