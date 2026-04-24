#!/usr/bin/env python3
"""
Graphical abstract — distributed dual-axis AF substrate framework.

Top band:    environment (exposome) as dashed upstream modifier.
Centre band: Healthy chamber → AF rewiring → Dual-axis substrate.
Bottom band: clinical translation tagline.

Colour discipline (consistent with Figs 1–6):
  - LA blue / RA red    : healthy identity
  - fibroblast green / cardiomyocyte mustard : axis identity
  - PuOr semantic       : matrix axis = orange, metabolic axis = mustard/brown
  - Dashed arrows       : hypothesised / upstream
  - Solid arrows        : established within this study
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import (
    Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon,
)
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

from plot_style import COL_DOUBLE, COMPARTMENT_COLORS, LA_COLOR, RA_COLOR, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "graphical_abstract.pdf"

C_FIB = COMPARTMENT_COLORS["fibroblast"]
C_FIB_DARK = "#3F6512"
C_CM = COMPARTMENT_COLORS["atrialCM"]
C_CM_DARK = "#7A5800"


def _round(ax, x, y, w, h, fc, ec="none", lw=0.5, z=1, rounding=0.014):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    ))


def _arrow(ax, p1, p2, color="#666666", lw=0.85, mut=8, dashed=False,
           curve=0.0, alpha=0.95):
    style = (0, (3, 2)) if dashed else "-"
    ax.add_patch(FancyArrowPatch(
        p1, p2,
        arrowstyle="-|>,head_length=2.8,head_width=1.7",
        connectionstyle=f"arc3,rad={curve}",
        mutation_scale=mut, lw=lw, color=color, linestyle=style,
        capstyle="round", joinstyle="round",
        shrinkA=2, shrinkB=2, alpha=alpha, zorder=10,
    ))


def _atrial_silhouette(ax, cx, cy, w=0.080, h=0.060,
                       la_color=LA_COLOR, ra_color=RA_COLOR):
    """Pair of stylised atria: LA (left, blue) and RA (right, red)."""
    # LA bulb (left)
    ax.add_patch(Ellipse((cx - w * 0.42, cy + h * 0.15),
                         w * 0.7, h * 0.85,
                         facecolor=la_color, edgecolor="white",
                         linewidth=0.6, alpha=0.85, zorder=4))
    # RA bulb (right)
    ax.add_patch(Ellipse((cx + w * 0.42, cy + h * 0.15),
                         w * 0.7, h * 0.85,
                         facecolor=ra_color, edgecolor="white",
                         linewidth=0.6, alpha=0.85, zorder=4))
    ax.text(cx - w * 0.42, cy + h * 0.15, "LA",
            ha="center", va="center", fontsize=5.6,
            color="white", fontweight="bold", zorder=5)
    ax.text(cx + w * 0.42, cy + h * 0.15, "RA",
            ha="center", va="center", fontsize=5.6,
            color="white", fontweight="bold", zorder=5)


def _icon_cm(ax, cx, cy, w=0.060, h=0.013, dim=False):
    """Slim cardiomyocyte rod with off-centre nucleus."""
    color = C_CM
    nuc_color = "#1E1308"
    body = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.0,rounding_size=0.005",
        facecolor=color, edgecolor=C_CM_DARK, linewidth=0.4,
        alpha=0.95, zorder=4)
    ax.add_patch(body)
    ax.add_patch(Ellipse((cx + 0.1 * w, cy), w * 0.12, h * 0.55,
                         facecolor=nuc_color, edgecolor="none",
                         alpha=0.9, zorder=5))
    # Optional dim mitochondria dots
    if dim:
        ax.add_patch(Ellipse((cx - 0.20 * w, cy), w * 0.08, h * 0.30,
                             facecolor="#D9B58A", edgecolor="none",
                             alpha=0.9, zorder=5))


def _icon_fibroblast(ax, cx, cy, w=0.034, h=0.020):
    """Compact 5-point star fibroblast."""
    import numpy as np
    n = 5
    pts = []
    for k in range(n * 2):
        theta = (k / (n * 2)) * 2 * np.pi - np.pi / 2
        r_w = (w / 2) if (k % 2 == 0) else (w * 0.18)
        r_h = (h / 2) if (k % 2 == 0) else (h * 0.22)
        pts.append((cx + np.cos(theta) * r_w,
                    cy + np.sin(theta) * r_h))
    ax.add_patch(Polygon(pts, closed=True,
                         facecolor=C_FIB, edgecolor=C_FIB_DARK,
                         linewidth=0.4, alpha=0.95, zorder=4))
    ax.add_patch(Ellipse((cx, cy), w * 0.22, h * 0.30,
                         facecolor="#1E1308", edgecolor="none",
                         alpha=0.95, zorder=5))


def _ecm_lines(ax, cx, cy, w, h, n=4, color="#A86F2C", lw=0.45, alpha=0.7):
    import numpy as np
    rng = np.random.default_rng(int(abs(cx * 71 + cy * 31) * 1e3) % 2**31)
    for k in range(n):
        y = cy + ((k - (n - 1) / 2) / max(n - 1, 1)) * h * 0.85
        x0 = cx - w / 2
        x1 = cx + w / 2
        midx = (x0 + x1) / 2
        midy = y + (1 if k % 2 == 0 else -1) * h * 0.10
        verts = [(x0, y), (midx, midy), (x1, y)]
        codes = [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3]
        ax.add_patch(PathPatch(MplPath(verts, codes),
                               facecolor="none", edgecolor=color,
                               linewidth=lw, alpha=alpha, zorder=3))


def main() -> None:
    set_nature_style()
    fig_w = COL_DOUBLE
    fig_h = 3.40
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # ============ TOP BAND — environment / exposome ============
    _round(ax, 0.04, 0.86, 0.92, 0.08, "#F4F4F4", "#BBBBBB", lw=0.45,
           rounding=0.014, z=1)
    fig.text(0.06, 0.905, "Environmental modifiers",
             fontsize=6.3, fontweight="bold", color="#222222",
             ha="left", va="center")
    chips = ["air pollution", "noise", "heat / weather", "built · social"]
    for k, lab in enumerate(chips):
        x = 0.30 + k * 0.16
        _round(ax, x, 0.875, 0.14, 0.05, "white", "#888888", lw=0.4,
               rounding=0.010, z=2)
        fig.text(x + 0.07, 0.900, lab, fontsize=5.4, color="#444444",
                 ha="center", va="center")

    # ============ MIDDLE — three boxes ============
    # Left: healthy
    _round(ax, 0.04, 0.30, 0.24, 0.50, "white", "#BBBBBB", lw=0.6,
           rounding=0.014, z=1)
    fig.text(0.16, 0.755, "Healthy chamber identity",
             ha="center", va="center", fontsize=6.4, fontweight="bold",
             color="#222222")
    _atrial_silhouette(ax, 0.16, 0.62, w=0.16, h=0.13)
    fig.text(0.16, 0.49,
             "stable LA / RA program architecture\n"
             "(metabolic-LA, structural-RA)",
             ha="center", va="center", fontsize=5.4, color="#444444",
             style="italic", linespacing=1.45)
    fig.text(0.16, 0.36, "Fig 1",
             ha="center", va="center", fontsize=5.0,
             color="#888888", fontweight="bold")

    # Centre: AF rewiring
    _round(ax, 0.32, 0.30, 0.30, 0.50, "#FFFAF0", "#C58A2D", lw=0.6,
           rounding=0.014, z=1)
    fig.text(0.47, 0.755, "AF chamber rewiring",
             ha="center", va="center", fontsize=6.4, fontweight="bold",
             color="#A55F00")
    fig.text(0.47, 0.685,
             "selective, distributed,\nbut structured",
             ha="center", va="center", fontsize=5.6,
             color="#5A3F00", linespacing=1.45, style="italic")

    # Two pathway labels with icons in middle
    fig.text(0.345, 0.585, "→ matrix\nremodeling",
             ha="left", va="center", fontsize=5.4, fontweight="bold",
             color=C_FIB_DARK, linespacing=1.2)
    _ecm_lines(ax, 0.535, 0.585, 0.10, 0.04, n=4, alpha=0.85)
    _icon_fibroblast(ax, 0.490, 0.585, w=0.034, h=0.022)
    _icon_fibroblast(ax, 0.578, 0.585, w=0.034, h=0.022)

    fig.text(0.345, 0.45, "→ metabolic\nattrition",
             ha="left", va="center", fontsize=5.4, fontweight="bold",
             color=C_CM_DARK, linespacing=1.2)
    _icon_cm(ax, 0.485, 0.45, w=0.060, h=0.014, dim=True)
    _icon_cm(ax, 0.565, 0.45, w=0.060, h=0.014, dim=True)

    fig.text(0.47, 0.36,
             "Fig 2 (chamber)  ·  Fig 3 (cells)  ·  Fig 4 (chronicity)",
             ha="center", va="center", fontsize=4.8,
             color="#888888")

    # Right: dual-axis substrate (split into 2 sub-cards)
    _round(ax, 0.66, 0.55, 0.30, 0.25, "#F1F8EE", C_FIB_DARK, lw=0.6,
           rounding=0.014, z=1)
    fig.text(0.81, 0.755,
             "Fibroblast proteoglycan /\nmatrix remodeling axis",
             ha="center", va="top", fontsize=5.8, fontweight="bold",
             color=C_FIB_DARK, linespacing=1.25)
    _ecm_lines(ax, 0.81, 0.625, 0.20, 0.04, n=5, alpha=0.95)
    _icon_fibroblast(ax, 0.74, 0.625, w=0.034, h=0.022)
    _icon_fibroblast(ax, 0.88, 0.625, w=0.034, h=0.022)
    fig.text(0.81, 0.575,
             "Fig 2–3–5C  ·  Fig 6A meta",
             ha="center", va="center", fontsize=4.6, color="#888888",
             style="italic")

    _round(ax, 0.66, 0.30, 0.30, 0.22, "#FFF9E6", "#BF812D", lw=0.6,
           rounding=0.014, z=1)
    fig.text(0.81, 0.50,
             "Appendage-weighted\nCM metabolic attrition axis",
             ha="center", va="top", fontsize=5.8, fontweight="bold",
             color=C_CM_DARK, linespacing=1.25)
    _icon_cm(ax, 0.74, 0.395, w=0.06, h=0.014, dim=True)
    _icon_cm(ax, 0.81, 0.395, w=0.06, h=0.014, dim=True)
    _icon_cm(ax, 0.88, 0.395, w=0.06, h=0.014, dim=True)
    fig.text(0.81, 0.345,
             "Fig 2–4  ·  Fig 6A meta",
             ha="center", va="center", fontsize=4.6, color="#888888",
             style="italic")

    # ============ Bottom merge box ============
    _round(ax, 0.34, 0.13, 0.32, 0.10, "#EFEFEF", "#444444", lw=0.7,
           rounding=0.014, z=2)
    fig.text(0.50, 0.18, "Distributed AF substrate",
             ha="center", va="center", fontsize=6.4, fontweight="bold",
             color="#1A1A1A")
    fig.text(0.50, 0.143,
             "(structured, no single FDR-robust dominant carrier)",
             ha="center", va="center", fontsize=4.6, color="#555555",
             style="italic")

    # ============ Arrows ============
    # environment → axes (dashed)
    for x_src in [0.37, 0.46, 0.54, 0.61]:
        _arrow(ax, (x_src, 0.875), (x_src + 0.10, 0.795),
               color="#888888", dashed=True, lw=0.55, mut=6, alpha=0.75)

    # healthy → AF rewiring
    _arrow(ax, (0.28, 0.55), (0.318, 0.55), color="#777777", lw=0.85)

    # AF rewiring → dual axes
    _arrow(ax, (0.62, 0.62), (0.66, 0.65),
           color=C_FIB_DARK, lw=0.95, mut=9, curve=-0.08)
    _arrow(ax, (0.62, 0.45), (0.66, 0.38),
           color=C_CM_DARK, lw=0.95, mut=9, curve=+0.08)

    # axes → distributed substrate
    _arrow(ax, (0.74, 0.55), (0.46, 0.225),
           color=C_FIB_DARK, lw=0.85, mut=8, curve=+0.10)
    _arrow(ax, (0.74, 0.30), (0.55, 0.225),
           color=C_CM_DARK, lw=0.85, mut=8, curve=-0.10)

    # ============ Bottom tagline ============
    fig.text(0.5, 0.06,
             "Substrate phenotyping beyond a binary AF / SR state",
             ha="center", va="center", fontsize=6.4, color="#1A1A1A",
             fontweight="bold")
    fig.text(0.5, 0.025,
             "matrix–metabolic two-axis framework  ·  cohort-stable structural arm  ·  "
             "appendage-weighted metabolic decline",
             ha="center", va="bottom", fontsize=4.8, color="#666666",
             style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
