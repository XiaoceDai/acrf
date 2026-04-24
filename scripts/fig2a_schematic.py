#!/usr/bin/env python3
"""
Figure 2A - Analysis schematic.

Three panels laid out horizontally, each given its own independent Axes and
a wide gutter in between so text never collides with graphics:

  panel 1  Healthy chamber prior  (10 LA/RA tiles + aCAI_bulk note)
  panel 2  Within-donor asymmetry (per-donor CAI, aCAI, aCAI_rel formulae)
  panel 3  AF rewiring classification (4 class tiles with mini SR/AF cartoon)

Arrow glyphs are drawn in fig coordinates in the gutters.

Output: vector PDF, no rasters.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

from plot_style import (
    COL_DOUBLE, LA_COLOR, RA_COLOR,
    REWIRING_CLASS_COLORS, set_nature_style,
)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "fig2a_schematic.pdf"

C_AMP = REWIRING_CLASS_COLORS["Amplified"]
C_PRE = REWIRING_CLASS_COLORS["Preserved"]
C_ATT = REWIRING_CLASS_COLORS["Attenuated"]
C_REV = REWIRING_CLASS_COLORS["Reversed"]


# ---------------------------------------------------------------------------
# Panel 1 - Healthy chamber prior
# ---------------------------------------------------------------------------

def draw_panel_1(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
    ax.text(0.02, 0.96, "1", fontsize=9, fontweight="bold",
            ha="left", va="top")
    ax.text(0.08, 0.96, "Healthy chamber prior",
            fontsize=8, fontweight="bold", ha="left", va="top")
    ax.text(0.02, 0.88,
            "10 stable LA/RA programs\n(Figure 1C, signed bulk NES)",
            fontsize=6.5, ha="left", va="top", color="#555555")

    # LA tiles
    tile_w, tile_h = 0.18, 0.042
    y0 = 0.74
    for i in range(5):
        y = y0 - i * (tile_h + 0.010)
        ax.add_patch(FancyBboxPatch(
            (0.04, y), tile_w, tile_h,
            boxstyle="round,pad=0,rounding_size=0.010",
            facecolor=LA_COLOR, edgecolor="white", linewidth=0.4,
        ))
    ax.text(0.04 + tile_w + 0.015, y0 + tile_h / 2, "LA-high",
            fontsize=7, ha="left", va="center",
            color=LA_COLOR, fontweight="bold")

    # RA tiles
    y0 = 0.44
    for i in range(5):
        y = y0 - i * (tile_h + 0.010)
        ax.add_patch(FancyBboxPatch(
            (0.04, y), tile_w, tile_h,
            boxstyle="round,pad=0,rounding_size=0.010",
            facecolor=RA_COLOR, edgecolor="white", linewidth=0.4,
        ))
    ax.text(0.04 + tile_w + 0.015, y0 + tile_h / 2, "RA-high",
            fontsize=7, ha="left", va="center",
            color=RA_COLOR, fontweight="bold")

    # Bulk reference note
    ax.text(0.02, 0.08,
            r"$|\mathrm{aCAI}_{bulk}(p)|$" + " from GSE112339\n"
            "(same singscore on Fig 1\nsource matrix)",
            fontsize=6.5, ha="left", va="bottom", color="#333333")


# ---------------------------------------------------------------------------
# Panel 2 - Within-donor asymmetry
# ---------------------------------------------------------------------------

def draw_panel_2(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
    ax.text(0.02, 0.96, "2", fontsize=9, fontweight="bold",
            ha="left", va="top")
    ax.text(0.08, 0.96, "Within-donor asymmetry",
            fontsize=8, fontweight="bold", ha="left", va="top")
    ax.text(0.02, 0.88, "per AF / SR donor, per program",
            fontsize=6.5, ha="left", va="top", color="#555555")

    # Donor cartoon: a pair of tissue disks.
    disk_r = 0.050
    la_x, la_y = 0.30, 0.74
    ra_x, ra_y = 0.70, 0.74
    ax.add_patch(plt.Circle((la_x, la_y), disk_r,
                            facecolor=LA_COLOR, edgecolor="none"))
    ax.add_patch(plt.Circle((ra_x, ra_y), disk_r,
                            facecolor=RA_COLOR, edgecolor="none"))
    ax.text(la_x, la_y - disk_r - 0.035, "LA", fontsize=7,
            ha="center", va="top", color=LA_COLOR, fontweight="bold")
    ax.text(ra_x, ra_y - disk_r - 0.035, "RA", fontsize=7,
            ha="center", va="top", color=RA_COLOR, fontweight="bold")
    ax.text(la_x, la_y + disk_r + 0.020, "score(LA)", fontsize=6.5,
            ha="center", va="bottom", color="#333333")
    ax.text(ra_x, ra_y + disk_r + 0.020, "score(RA)", fontsize=6.5,
            ha="center", va="bottom", color="#333333")

    # Minus sign between disks.
    ax.text((la_x + ra_x) / 2, la_y, "-", fontsize=12,
            ha="center", va="center", color="#333333", fontweight="bold")

    # Three formulas stacked.
    for y, expr in [
        (0.50, r"$\mathrm{CAI} = \mathrm{score}(LA) - \mathrm{score}(RA)$"),
        (0.40, r"$\mathrm{aCAI} = s_p \cdot \mathrm{CAI}$"),
        (0.28, r"$\mathrm{aCAI}_{rel} = \mathrm{aCAI}\ /\ |\mathrm{aCAI}_{bulk}(p)|$"),
    ]:
        ax.text(0.5, y, expr, fontsize=7.5, ha="center", va="center",
                color="#222222")

    # Interpretation key at the bottom.  Use a 2 x 2 grid so the swatch +
    # text fits comfortably inside each cell.
    ax.text(0.03, 0.22,
            r"$\mathrm{aCAI}_{rel}$ interpretation:",
            fontsize=6.5, ha="left", va="center", fontweight="bold",
            color="#333333")
    keys = [
        (">1",          "amplified",  C_AMP),
        (r"$\approx$1", "preserved",  C_PRE),
        ("0-1",         "attenuated", C_ATT),
        ("<0",          "reversed",   C_REV),
    ]
    sw_w, sw_h = 0.030, 0.035
    xs = [0.04, 0.52]
    ys = [0.12, 0.04]
    for k, (val, label, col) in enumerate(keys):
        cx = xs[k % 2]
        cy = ys[k // 2]
        ax.add_patch(Rectangle((cx, cy), sw_w, sw_h,
                               facecolor=col, edgecolor="#555555", linewidth=0.3))
        ax.text(cx + sw_w + 0.010, cy + sw_h / 2, f"{val}  {label}",
                fontsize=6, ha="left", va="center", color="#222222")


# ---------------------------------------------------------------------------
# Panel 3 - AF rewiring classification
# ---------------------------------------------------------------------------

def draw_panel_3(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
    ax.text(0.02, 0.96, "3", fontsize=9, fontweight="bold",
            ha="left", va="top")
    ax.text(0.08, 0.96, "AF rewiring classification",
            fontsize=8, fontweight="bold", ha="left", va="top")
    ax.text(0.02, 0.88,
            "SR vs AF on " + r"$\Delta\mathrm{aCAI}_{rel}$" +
            ", permutation",
            fontsize=6.5, ha="left", va="top", color="#555555")

    # Four class "cards" stacked vertically.  Deep-fill cards (Amplified,
    # Reversed) use white text; pale cards use dark text.
    classes = [
        ("Amplified",  C_AMP, "AF > SR",              "white"),
        ("Preserved",  C_PRE, r"AF $\approx$ SR",     "#222222"),
        ("Attenuated", C_ATT, "AF < SR (same sign)",  "#222222"),
        ("Reversed",   C_REV, "AF < 0 (sign flipped)", "white"),
    ]
    card_w, card_h = 0.94, 0.100
    x0 = 0.03
    y0 = 0.64
    for k, (name, col, desc, tcol) in enumerate(classes):
        y = y0 - k * (card_h + 0.022)
        ax.add_patch(FancyBboxPatch(
            (x0, y), card_w, card_h,
            boxstyle="round,pad=0,rounding_size=0.012",
            facecolor=col, edgecolor="#444444", linewidth=0.5,
        ))
        ax.text(x0 + 0.025, y + card_h / 2 + 0.013, name,
                fontsize=7.5, fontweight="bold", ha="left", va="center",
                color=tcol)
        desc_color = "white" if tcol == "white" else "#555555"
        ax.text(x0 + 0.025, y + card_h / 2 - 0.020, desc,
                fontsize=6, ha="left", va="center", color=desc_color)

    # Footer note about Low-support and trend-class.
    ax.text(0.03, 0.05,
            "Low-support: SR fails baseline support (dropped).\n"
            "Trend-class (exploratory): same rules, permutation gate ignored.",
            fontsize=5.5, ha="left", va="bottom", color="#666666")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    set_nature_style()

    fig_w = COL_DOUBLE
    fig_h = 2.9
    fig = plt.figure(figsize=(fig_w, fig_h))

    # Each panel occupies a column. Leave generous gutters for arrows.
    ax1 = fig.add_axes([0.02, 0.10, 0.24, 0.82])
    ax2 = fig.add_axes([0.32, 0.10, 0.32, 0.82])
    ax3 = fig.add_axes([0.71, 0.10, 0.28, 0.82])
    draw_panel_1(ax1)
    draw_panel_2(ax2)
    draw_panel_3(ax3)

    # Arrows in figure coords (horizontal, centred vertically).
    arr_y = 0.52
    for x0, x1 in [(0.265, 0.315), (0.65, 0.705)]:
        fig.patches.append(FancyArrowPatch(
            (x0, arr_y), (x1, arr_y),
            transform=fig.transFigure,
            arrowstyle="-|>", mutation_scale=12,
            linewidth=1.0, color="#555555",
            shrinkA=0, shrinkB=0,
        ))

    # Methods footnote (full figure width).
    fig.text(0.5, 0.03,
             "Because rhythm is a donor-level attribute (every donor is entirely SR or "
             "entirely AF), the design is rank-deficient for sample-level "
             r"donor $\times$ chamber $\times$ rhythm interaction; "
             "rewiring is therefore quantified by donor-level aCAI with "
             "rhythm-label permutation.",
             fontsize=5.5, ha="center", va="bottom", color="#555555",
             wrap=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
