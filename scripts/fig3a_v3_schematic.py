#!/usr/bin/env python3
"""
Fig 3a v3 - schematic + data declaration.

A pure schematic (no data ingest) explaining how Fig 3 is computed:

    GSE112339 LA/RA bulk DE rank
        |
        |  fgsea LE extraction (Subramanian 2005 / fgsea default)
        v
    NCP_v3 leading-edge gene sets per program
        |
        |  rank singscore on GSE255612 LA snRNA pseudobulk
        v
    donor x compartment x program score
        |
        |  AF (n=18) vs Ctrl (n=16) per (program, ct), Monte-Carlo permutation
        v
    family x compartment dScore + p_perm  =  Fig 3b heatmap
        |
        |  classify against Fig 1F HCA baseline LA carrier
        v
    Fidelity / Relay / Extension / Below-detection (Fig 3c)
        |
        |  case-deep-dives: F15 (relay), F16 (extension), F4 BCAA (sub-family)
        v
    Fig 3d / 3e / 3f
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT_PDF = ROOT / "figures" / "fig3_v3_a_schematic.pdf"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
ACCENT = "#C8A076"
GREY = "#666666"
LIGHT = "#F2F2F2"


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def stage_box(ax, x, y, w, h, label_top, label_body, color=LIGHT, edgecolor="#888",
              body_align="center"):
    ax.add_patch(mpatches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.005,rounding_size=0.012",
        facecolor=color, edgecolor=edgecolor, linewidth=0.6,
    ))
    ax.text(x + w / 2.0, y + h - 0.025, label_top,
            ha="center", va="top", fontsize=7.2, fontweight="bold", color="#1A1A1A")
    if body_align == "left":
        ax.text(x + 0.012, y + h - 0.060, label_body,
                ha="left", va="top", fontsize=6.0, color="#333333")
    else:
        ax.text(x + w / 2.0, y + h / 2.0 - 0.005, label_body,
                ha="center", va="center", fontsize=6.0, color="#333333")


def arrow(ax, x1, y1, x2, y2, label=None):
    ax.annotate(
        "", xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle="->", color="#1A1A1A", lw=0.7),
    )
    if label:
        ax.text((x1 + x2) / 2.0, (y1 + y2) / 2.0 + 0.012, label,
                fontsize=5.6, ha="center", va="bottom", color="#444444", style="italic")


def main():
    set_style()

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 95 * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # ---- title ----
    ax.text(0.04, 0.95, "Fig 3a  AF cellular substrate validation pipeline (LA snRNA)",
            fontsize=9.0, fontweight="bold", color="#1A1A1A")

    # ---- LA-arm only data box (top right) ----
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.71, 0.78), 0.27, 0.14,
        boxstyle="round,pad=0.005,rounding_size=0.010",
        facecolor="#FFF4E6", edgecolor=ACCENT, linewidth=0.7,
    ))
    ax.text(0.845, 0.89, "GSE255612 (Hill et al. 2024)", fontsize=7.0, fontweight="bold",
            color="#1A1A1A", ha="center")
    ax.text(0.845, 0.85, "snRNA-seq, LEFT ATRIUM only", fontsize=6.0, color=ACCENT,
            ha="center", fontweight="bold")
    ax.text(0.845, 0.815, "18 AF + 16 Ctrl donors\n9 reliable cell compartments",
            fontsize=5.7, color="#444444", ha="center")

    # ---- Stage 1: bulk DE rank ----
    stage_box(ax, 0.03, 0.55, 0.22, 0.18,
              "1. Healthy bulk LA/RA DE rank",
              "GSE112339\nlimma signed t-stat\n(LA > 0,  RA < 0)",
              color="#EAF1F7", edgecolor=RA_DEEP)

    # ---- Stage 2: LE extraction ----
    stage_box(ax, 0.30, 0.55, 0.22, 0.18,
              "2. fgsea leading-edge",
              "extract LE per program\n(NCP_v3 panel:\n30 program / 16 family)",
              color=LIGHT, edgecolor="#888")

    # ---- Stage 3: singscore ----
    stage_box(ax, 0.57, 0.55, 0.22, 0.18,
              "3. rank singscore",
              "GSE255612 LA snRNA\npseudobulk (donor x ct)\nFig3 v3 panel LE",
              color="#FFF4E6", edgecolor=ACCENT)

    # ---- Stage 4: stats ----
    stage_box(ax, 0.03, 0.30, 0.22, 0.18,
              "4. AF vs Ctrl test",
              "Monte-Carlo permutation\n(n = 10,000, seed = 20260418)\n"
              "family- and program-level\n"
              "\u2192 \u0394Score (AF \u2212 Ctrl)\n   per (family, compartment)",
              color="#F2EAEA", edgecolor=LA_DEEP, body_align="left")

    # ---- Stage 5: heatmap (Fig 3b) ----
    stage_box(ax, 0.30, 0.30, 0.22, 0.18,
              "5. Family x compartment\nheatmap  ->  Fig 3b",
              "16 row x 9 col\n"
              "color = \u0394Score\n"
              "scale \u00b1 0.04\n"
              "p<0.05 boxed",
              color=LIGHT, edgecolor="#888", body_align="left")

    # ---- Stage 6: classify against Fig 1F carrier (vertical 4-class list) ----
    stage_box(ax, 0.57, 0.30, 0.22, 0.18,
              "6. Classify vs Fig 1F\nbaseline LA carrier",
              "\u251c Fidelity\n"
              "\u251c Relay\n"
              "\u251c Extension\n"
              "\u2514 Below-detection\n"
              "\u2192 Fig 3c",
              color="#E8F0E8", edgecolor="#6E9B72", body_align="left")

    # ---- Stage 7: case studies ----
    stage_box(ax, 0.84, 0.30, 0.14, 0.18,
              "7. Case studies",
              "Fig 3d  F15 NTRK2\n  (relay)\n"
              "Fig 3e  F16 Keratan\n  (extension)\n"
              "Fig 3f  F4 BCAA\n  (sub-family)",
              color=LIGHT, edgecolor="#888", body_align="left")

    # ---- arrows ----
    arrow(ax, 0.25, 0.64, 0.30, 0.64)
    arrow(ax, 0.52, 0.64, 0.57, 0.64)
    # Stage 3 -> Stage 4: route through a small Ctrl-vs-AF inset that
    # demonstrates the dScore primitive visually.
    arrow(ax, 0.68, 0.55, 0.50, 0.495)
    arrow(ax, 0.42, 0.495, 0.18, 0.48, label="\u0394Score = AF \u2212 Ctrl")
    arrow(ax, 0.25, 0.39, 0.30, 0.39, label="dScore matrix")
    arrow(ax, 0.52, 0.39, 0.57, 0.39)
    arrow(ax, 0.79, 0.39, 0.84, 0.39)

    # ---- Inset: Ctrl vs AF mini dot plot (between stage 3 and stage 4) ----
    ix, iy, iw, ih = 0.42, 0.49, 0.08, 0.06
    ax.add_patch(mpatches.FancyBboxPatch(
        (ix, iy), iw, ih,
        boxstyle="round,pad=0.003,rounding_size=0.006",
        facecolor="white", edgecolor="#888", linewidth=0.4,
    ))
    # 3 Ctrl dots (left) at one y-level, 3 AF dots (right) at slightly lower
    rng_demo = np.random.default_rng(7)
    ctrl_y = 0.62 + rng_demo.normal(0, 0.05, size=3) * ih
    af_y = 0.42 + rng_demo.normal(0, 0.05, size=3) * ih
    ctrl_x = np.linspace(ix + 0.012, ix + 0.030, 3)
    af_x = np.linspace(ix + 0.045, ix + 0.063, 3)
    ax.scatter(ctrl_x, iy + ih * np.array(ctrl_y), s=8, color=RA_DEEP, edgecolor="white", linewidth=0.4, zorder=5)
    ax.scatter(af_x, iy + ih * np.array(af_y), s=8, color=LA_DEEP, edgecolor="white", linewidth=0.4, zorder=5)
    # mean lines
    ax.plot([ix + 0.008, ix + 0.034], [iy + ih * 0.62, iy + ih * 0.62],
            color=RA_DEEP, lw=0.5)
    ax.plot([ix + 0.041, ix + 0.067], [iy + ih * 0.42, iy + ih * 0.42],
            color=LA_DEEP, lw=0.5)
    # Down arrow indicating delta
    ax.annotate("", xy=(ix + iw / 2.0, iy + ih * 0.42),
                xytext=(ix + iw / 2.0, iy + ih * 0.62),
                arrowprops=dict(arrowstyle="->", color="#1A1A1A", lw=0.5))
    ax.text(ix + iw / 2.0 + 0.005, iy + ih * 0.52, "\u0394",
            fontsize=6.4, color="#1A1A1A", fontweight="bold", va="center")
    ax.text(ix + iw / 2.0, iy + ih + 0.004, "Ctrl   AF",
            fontsize=5.0, color="#666666", ha="center", va="bottom")

    # ---- limitation banner ----
    ax.add_patch(mpatches.Rectangle((0.03, 0.04), 0.95, 0.18,
                                     facecolor="#FFF8F0", edgecolor=ACCENT, linewidth=0.6))
    ax.text(0.05, 0.20, "Methodological limitations (acknowledged in Discussion)",
            fontsize=7.0, fontweight="bold", color="#1A1A1A", va="top")
    ax.text(0.05, 0.165,
            "(i)  GSE255612 contains LEFT ATRIUM samples only; RA-side cellular activity changes are not observable in this dataset.",
            fontsize=5.9, color="#333333", va="top")
    ax.text(0.05, 0.13,
            "(ii) Compartment-level pseudobulk dilutes within-cell-type rare-population signals; bulk rewiring may reside below per-compartment detection.",
            fontsize=5.9, color="#333333", va="top")
    ax.text(0.05, 0.095,
            "(iii) Cellular relay vs sub-cell-type composition shift cannot be fully separated at the pseudobulk level (single-nucleus subtype analysis is reserved as future work).",
            fontsize=5.9, color="#333333", va="top")
    ax.text(0.05, 0.060,
            "(iv) Fig 3 baseline carrier reference comes from the HCA paired LA + RA atlas (Fig 1F); the AF-time activity-shift cell type is determined within GSE255612 LA arm only.",
            fontsize=5.9, color="#333333", va="top")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
