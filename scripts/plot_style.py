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

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Page widths (Nature Medicine column geometry)
MM = 1.0 / 25.4
COL_SINGLE  = 89  * MM   # 3.5"
COL_ONEHALF = 120 * MM   # 4.7"
COL_DOUBLE  = 183 * MM   # 7.2"
PAGE_HEIGHT = 240 * MM   # max usable single-page height


# Final composite cell sizes (locked).
# Each panel is rendered at exactly this size (mm), saved with
# `bbox_inches=None`, and composed into Figure 1 with NO further scaling.
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

# Color tokens — locked semantic palette

# Chromatic axis (the ONLY use of red / blue in any panel)
LA_COLOR = "#C97A74"   # muted warm red, LA
RA_COLOR = "#7D99B8"   # muted cool blue, RA

# Deeper endpoints for heatmaps / emphasis (same hue family as LA/RA)
LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"

# Rhythm axis (AF vs SR) — distinct semantic axis from LA/RA chamber colours.
# Use for stacked bars / boxplots of rhythm; keep consistent across Fig 2 / S1.
AF_RHYTHM_COLOR = "#B2182B"
SR_RHYTHM_COLOR = "#2166AC"

# Healthy reference donors (e.g. GSE112339) when shown beside AF/SR stacks
HEALTHY_REF_COLOR = "#6C7B8A"

# Sex composition (supplementary demographics)
SEX_MALE_COLOR = "#5AAE61"
SEX_FEMALE_COLOR = "#9970AB"
SEX_NA_COLOR = "#BDBDBD"

# Qualitative stacks (disease background, many categories) — colorblind-friendly
QUALITATIVE_CATEGORY_COLORS = [
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666",
    "#999999", "#BDBDBD",
]

# Structural greys — for text, lines, background shading
FG_COLOR    = "#1A1A1A"   # primary text
MUTED_COLOR = "#555555"   # secondary text / annotations
FAINT_COLOR = "#888888"   # axis ticks, light annotations
LINE_COLOR  = "#CCCCCC"   # gridlines, structure separators
BG_COLOR    = "#F5F5F5"   # neutral panel fill

# Single-purpose accents — used at most twice per figure
HIGHLIGHT_COLOR = "#C8A076"   # muted warm orange — main emphasis
SUPPORT_COLOR   = "#6E9B72"   # muted green — stability / agreement / "8/8" callouts

# Diverging colormaps — each tied to a single semantic axis

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

# Typography (locked).  6.5 pt absolute floor; 7 pt default for main-text readability.
FS_PANEL_LETTER = 10     # bold lowercase a, b, c…
FS_PANEL_TITLE  = 9      # bold panel title (one line ideally)
FS_AXIS_LABEL   = 8
FS_TICK         = 7      # ticks / legend / in-figure annotations
FS_ANNOT        = 7      # in-figure annotations
FS_LEGEND       = 7      # legend text
FS_FOOTNOTE     = 7      # italic muted footnote, figure-level

# Submission / main-text minimum visible font size (legacy exploratory scripts may use 6.5)
FONT_FLOOR_MAIN = 7.0
FONT_FLOOR_HARD = 6.5

# Line weights (locked)
LW_AXIS = 0.5
LW_THIN = 0.4
LW_BOLD = 0.8


# Style entry point

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
        "xtick.direction":   "out",
        "ytick.direction":   "out",
        "xtick.major.pad":   2.5,
        "ytick.major.pad":   2.5,
        "lines.linewidth":   LW_BOLD,
        "patch.linewidth":   LW_THIN,

        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.grid":          False,
        "axes.axisbelow":     True,
        "axes.facecolor":     "white",
        "figure.facecolor":   "white",

        "axes.edgecolor":     FG_COLOR,
        "axes.labelcolor":    FG_COLOR,
        "xtick.color":        FG_COLOR,
        "ytick.color":        FG_COLOR,
        "text.color":         FG_COLOR,
        "axes.unicode_minus": False,

        "figure.dpi":         150,
        "savefig.dpi":        300,
        "savefig.facecolor":  "white",
        "savefig.edgecolor":  "none",
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


# Helper utilities used across all panels

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


def cohort_axis_labels(ax: plt.Axes, cohorts: Sequence[str], *, rotation: float = 30.0) -> None:
    """Standard GEO-style cohort tick labels (rotation + padding)."""
    ax.set_xticks(np.arange(len(cohorts)))
    ax.set_xticklabels(cohorts, rotation=rotation, ha="right")
    ax.tick_params(axis="x", pad=2)


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


# Layout helpers for readability (added 2026-07-09)

def enforce_font_floor(size: float, *, main_text: bool = False) -> float:
    """Return the largest of the requested size and the locked floor.

    main_text=True enforces the 7 pt main-text minimum; otherwise 6.5 pt.
    """
    floor = FONT_FLOOR_MAIN if main_text else FONT_FLOOR_HARD
    return max(size, floor)


def bump_figure_size(figsize: tuple[float, float], *, dh_mm: float = 12, dw_mm: float = 0) -> tuple[float, float]:
    """Increase figure height (and optionally width) when larger fonts cause overflow.

    Keep the locked mm-to-inch conversion so canvas sizes stay in physical units.
    """
    w, h = figsize
    return (w + dw_mm * MM, h + dh_mm * MM)


def reserve_header_strip(figsize: tuple[float, float], n_lines: int = 2,
                         *, line_height_mm: float = 4.2, margin_mm: float = 2.5) -> tuple[float, float, float]:
    """Return (new_figsize, new_gridspec_top, header_y0) that reserves space for n_lines of title.

    Example: figsize, top, y0 = reserve_header_strip(figsize, n_lines=2)
             gs = fig.add_gridspec(..., top=top)
             fig.text(0.03, y0, title, ...)
             fig.text(0.03, y0 - 0.04, subtitle, ...)
    """
    strip_h_mm = n_lines * line_height_mm + margin_mm
    new_figsize = bump_figure_size(figsize, dh_mm=strip_h_mm)
    # Convert reserved strip to figure fraction (approximate, assuming ~100% height used)
    strip_frac = strip_h_mm * MM / new_figsize[1]
    top = 1.0 - strip_frac - 0.02
    y0 = 1.0 - margin_mm * MM / new_figsize[1] - 0.005
    return new_figsize, top, y0


def figure_header(fig: plt.Figure, panel_letter: str, title: str, *,
                  subtitle: str | None = None,
                  x: float = 0.034, top_reserve: float = 0.08) -> float:
    """Place a single unified header strip at the top of the figure.

    Returns the bottom y coordinate of the subtitle so callers can position content below.
    """
    y0 = 1.0 - top_reserve + 0.01
    fig.text(x, y0, panel_letter, fontsize=FS_PANEL_LETTER, fontweight="bold",
             ha="left", va="baseline", color=FG_COLOR)
    title_x = x + 0.045
    fig.text(title_x, y0, title, fontsize=FS_PANEL_TITLE + 1, fontweight="bold",
             ha="left", va="baseline", color=FG_COLOR)
    if subtitle:
        y1 = y0 - 0.035
        fig.text(title_x, y1, subtitle, fontsize=FS_TICK, ha="left", va="baseline",
                 color=MUTED_COLOR)
        return y1 - 0.015
    return y0 - 0.015


def figure_footnote(fig: plt.Figure, text: str,
                    y: float = 0.014,
                    color: str = MUTED_COLOR) -> None:
    """Italic muted footnote at the very bottom of the figure (now at 7 pt)."""
    fig.text(0.012, y, text, fontsize=FS_FOOTNOTE,
             ha="left", va="bottom", color=color, style="italic")
