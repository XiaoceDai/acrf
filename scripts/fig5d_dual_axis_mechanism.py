#!/usr/bin/env python3
"""
Figure 5D — Dual-axis mechanism cartoon (clean, NC-style).

Three columns: Healthy chamber → AF rewiring → Distributed substrate.

Design discipline (post-polish):
  - No BioRender-style icons; just typographic boxes + colored chips.
  - Strong vertical grid (3 conceptual rows: identity-LA, identity-RA, message).
  - Single uniform arrow style.
  - Colour codes consistent with Fig 1–6 (LA blue, RA red, axis green/brown).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from plot_style import COL_DOUBLE, COMPARTMENT_COLORS, LA_COLOR, RA_COLOR, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "fig5d_dual_axis_mechanism.pdf"

C_FIB_FILL = COMPARTMENT_COLORS["fibroblast"]
C_FIB_LINE = "#3F6512"
C_CM_FILL  = COMPARTMENT_COLORS["atrialCM"]
C_CM_LINE  = "#7A5800"
C_METAB    = "#8C510A"


def _round(ax, x, y, w, h, fc, ec="#BBBBBB", lw=0.5, z=1, alpha=1.0,
           rounding=0.012):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha, zorder=z,
    ))


def _arr(ax, p1, p2, color="#666666", lw=0.7, mut=8, curve=0.0, alpha=0.95):
    ax.add_patch(FancyArrowPatch(
        p1, p2,
        arrowstyle="-|>,head_length=2.6,head_width=1.6",
        connectionstyle=f"arc3,rad={curve}",
        mutation_scale=mut, lw=lw, color=color, alpha=alpha,
        capstyle="round", joinstyle="round",
        shrinkA=2, shrinkB=2, zorder=10,
    ))


def _chip(ax, x, y, w, h, fc, label, label_color="white",
          fontsize=5.6, weight="bold"):
    """A small colored chip with centered label (cell-type / program tag)."""
    _round(ax, x, y, w, h, fc=fc, ec="white", lw=0.4, z=4, rounding=0.008)
    ax.text(x + w / 2, y + h / 2, label,
            ha="center", va="center", fontsize=fontsize,
            color=label_color, fontweight=weight, zorder=5)


def main() -> None:
    set_nature_style()
    fig_w = COL_DOUBLE
    fig_h = 3.10
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # ============ LEFT COLUMN — healthy chamber identity ============
    _round(ax, 0.04, 0.18, 0.22, 0.74, "white", "#CCCCCC", lw=0.5)
    ax.text(0.15, 0.875, "Healthy chamber identity",
            ha="center", va="center", fontsize=6.6, fontweight="bold",
            color="#222222")

    # LA half
    _round(ax, 0.055, 0.555, 0.19, 0.275, "#EEF5FB", LA_COLOR,
           lw=0.45, alpha=0.45)
    ax.text(0.15, 0.78, "LA  ·  metabolic / translational",
            ha="center", va="center", fontsize=5.6, color=LA_COLOR,
            fontweight="bold")
    # tag chips for the 3 LA-biased families
    chip_w, chip_h = 0.052, 0.034
    for i, lab in enumerate(["OXPHOS", "Resp ETC", "Translation"]):
        _chip(ax, 0.069 + i * 0.058, 0.685, chip_w, chip_h, LA_COLOR, lab,
              fontsize=4.6)
    ax.text(0.15, 0.605,
            "donor-stable LA-biased program family",
            ha="center", va="center", fontsize=4.7, color="#444444",
            style="italic")

    # RA half
    _round(ax, 0.055, 0.225, 0.19, 0.275, "#FBEEEE", RA_COLOR,
           lw=0.45, alpha=0.45)
    ax.text(0.15, 0.45, "RA  ·  ECM / glycan",
            ha="center", va="center", fontsize=5.6, color=RA_COLOR,
            fontweight="bold")
    for i, lab in enumerate(["ECM-PG", "HS-GAG", "BMP / RUNX2"]):
        _chip(ax, 0.069 + i * 0.058, 0.355, chip_w, chip_h, RA_COLOR, lab,
              fontsize=4.6)
    ax.text(0.15, 0.275,
            "donor-stable RA-biased program family",
            ha="center", va="center", fontsize=4.7, color="#444444",
            style="italic")

    # ============ MIDDLE COLUMN — AF rewiring ============
    ax.text(0.475, 0.875, "AF chamber rewiring",
            ha="center", va="center", fontsize=6.6, fontweight="bold",
            color="#A55F00")

    # central locked sentence
    _round(ax, 0.295, 0.745, 0.36, 0.085, "#FFFAF0", "#C58A2D",
           lw=0.4, alpha=0.45)
    ax.text(0.475, 0.788,
            '"AF selectively rewires chamber identity,\nit does not replace it."',
            ha="center", va="center", fontsize=5.8, color="#5A3F00",
            style="italic", linespacing=1.4)

    # two transition rows
    _round(ax, 0.295, 0.495, 0.36, 0.215, "#F1F8EE", C_FIB_LINE,
           lw=0.4, alpha=0.30)
    ax.text(0.475, 0.665,
            "→ proteoglycan / matrix remodeling",
            ha="center", va="center", fontsize=5.7, color=C_FIB_LINE,
            fontweight="bold")
    # 3 small fibroblast tags
    for i, lab in enumerate(["fibroblast", "ECM-PG", "HS-GAG"]):
        _chip(ax, 0.330 + i * 0.105, 0.595, 0.092, 0.028,
              C_FIB_FILL, lab, fontsize=4.8)
    ax.text(0.475, 0.535,
            "PRELP · OMD · FMOD · CHST3",
            ha="center", va="center", fontsize=4.9, color=C_FIB_LINE,
            style="italic")

    _round(ax, 0.295, 0.225, 0.36, 0.215, "#FFF4E0", "#A55F00",
           lw=0.4, alpha=0.30)
    ax.text(0.475, 0.395,
            "→ metabolic attrition",
            ha="center", va="center", fontsize=5.7, color=C_METAB,
            fontweight="bold")
    for i, lab in enumerate(["atrialCM", "OXPHOS", "Translation"]):
        _chip(ax, 0.330 + i * 0.105, 0.325, 0.092, 0.028,
              C_CM_FILL, lab, label_color="#3F2A00", fontsize=4.8)
    ax.text(0.475, 0.265,
            "OXPHOS↓  ·  LAA-weighted  ·  duration-linked",
            ha="center", va="center", fontsize=4.9, color=C_METAB,
            style="italic")

    # arrows from healthy → middle
    _arr(ax, (0.246, 0.69), (0.295, 0.65),
         color=LA_COLOR, lw=0.85, alpha=0.55)
    _arr(ax, (0.246, 0.355), (0.295, 0.36),
         color=RA_COLOR, lw=0.85, alpha=0.45)

    # ============ RIGHT COLUMN — chronic substrate ============
    _round(ax, 0.69, 0.555, 0.275, 0.275, "#F1F8EE", C_FIB_LINE,
           lw=0.55, alpha=0.45)
    ax.text(0.827, 0.79,
            "Fibroblast proteoglycan /\nmatrix remodeling axis",
            ha="center", va="top", fontsize=6.2, fontweight="bold",
            color=C_FIB_LINE, linespacing=1.25)
    ax.text(0.827, 0.665,
            "ECM proteoglycans ↑\nfibroblast-localized · spatially anchored",
            ha="center", va="center", fontsize=5.0, color="#333333",
            linespacing=1.4, style="italic")
    ax.text(0.827, 0.575, "Fig 2 · 3 · 5C",
            ha="center", va="center", fontsize=4.6, color="#777777",
            style="italic")

    _round(ax, 0.69, 0.225, 0.275, 0.275, "#FFF9E6", "#BF812D",
           lw=0.55, alpha=0.45)
    ax.text(0.827, 0.46,
            "Appendage-weighted CM\nmetabolic attrition axis",
            ha="center", va="top", fontsize=6.2, fontweight="bold",
            color=C_METAB, linespacing=1.25)
    ax.text(0.827, 0.335,
            "OXPHOS ↓  ·  duration-associated decline\nstronger in LAA than LAPV",
            ha="center", va="center", fontsize=5.0, color="#333333",
            linespacing=1.4, style="italic")
    ax.text(0.827, 0.250, "Fig 2 · 4",
            ha="center", va="center", fontsize=4.6, color="#777777",
            style="italic")

    # AF rewiring → right axes
    _arr(ax, (0.65, 0.660), (0.690, 0.690),
         color=C_FIB_LINE, lw=0.95, mut=10)
    _arr(ax, (0.65, 0.340), (0.690, 0.360),
         color=C_METAB, lw=0.95, mut=10)

    # substrate pill
    _round(ax, 0.69, 0.075, 0.275, 0.085, "#EFEFEF", "#444444", lw=0.7)
    ax.text(0.827, 0.118,
            "Distributed atrial remodeling substrate",
            ha="center", va="center", fontsize=6.2, fontweight="bold",
            color="#1A1A1A")
    _arr(ax, (0.760, 0.555), (0.795, 0.160),
         color=C_FIB_LINE, lw=0.7, curve=+0.10, alpha=0.75)
    _arr(ax, (0.895, 0.225), (0.860, 0.160),
         color=C_METAB, lw=0.7, curve=-0.10, alpha=0.75)

    # bottom caveat
    ax.text(0.5, 0.025,
            "No single FDR-robust dominant carrier; remodeling is "
            "distributed but structured (Fig 3 Unresolved carriers).",
            ha="center", va="bottom", fontsize=5.0, color="#666666",
            style="italic")

    # Title + label
    fig.suptitle(
        "Distributed dual-axis mechanism",
        fontsize=7.4, y=0.965, fontweight="bold",
    )
    ax.text(0.012, 0.97, "d",
            fontsize=10, fontweight="bold", ha="left", va="top",
            transform=fig.transFigure)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
