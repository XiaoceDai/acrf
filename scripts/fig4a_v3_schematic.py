#!/usr/bin/env python3
"""
Fig 4a v3 - schematic + power statement banner.

Mirrors fig3a_v3_schematic.py 7-stage flowchart style.

Pipeline:
    1. GSE41177 LAA + LAPV bulk RNA-seq
    2. Per-sample rank singscore on Fig 3 v3 LE genes
    3. Family-mean across programs in family
    4. Spearman rho vs duration_months (AF-only n = 12)
    5. Bootstrap 95% CI + permutation p
    6. 3-class effect-size partition (Strong / Weak / Stable)
    7. Case-study scatter (Fig 4d/e/f)

Limitations banner explicitly states n=12 power constraint.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT_PDF = ROOT / "figures" / "fig4_v3_a_schematic.pdf"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
ACCENT = "#C8A076"
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
    fig = plt.figure(figsize=(190 * MM, 100 * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.04, 0.95, "Fig 4a  AF-duration chronicity pipeline (LAA + LAPV bulk)",
            fontsize=9.0, fontweight="bold", color="#1A1A1A")

    # GSE41177 dataset banner (top right)
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.66, 0.78), 0.32, 0.14,
        boxstyle="round,pad=0.005,rounding_size=0.010",
        facecolor="#FFF4E6", edgecolor=ACCENT, linewidth=0.7,
    ))
    ax.text(0.82, 0.89, "GSE41177 (Yeh et al. 2013)", fontsize=7.0, fontweight="bold",
            color="#1A1A1A", ha="center")
    ax.text(0.82, 0.85, "Affy HG-U133 Plus 2.0,  paired LAA + LAPV", fontsize=6.0,
            color=ACCENT, ha="center", fontweight="bold")
    ax.text(0.82, 0.815, "12 valvular AF donors with numeric AF duration\n(months);  3 SR controls (visual anchor)",
            fontsize=5.7, color="#444444", ha="center")

    # Stage 1: bulk
    stage_box(ax, 0.03, 0.55, 0.20, 0.18,
              "1. GSE41177 LAA + LAPV bulk",
              "12 AF (numeric mo) + 3 SR\n+ 1 AF (>1M only,\nused in gt1m_eq1 sensitivity)\nAffy HG-U133 Plus 2.0",
              color="#EAF1F7", edgecolor=RA_DEEP, body_align="left")

    # Stage 2: rank singscore
    stage_box(ax, 0.27, 0.55, 0.20, 0.18,
              "2. Rank singscore",
              "Per-sample rank-based singscore\n(LE genes per program;\nrank/N centred to [-0.5, +0.5];\nFig 3 v3 LE; Foroutan 2018)",
              color=LIGHT, edgecolor="#888", body_align="left")

    # Stage 3: family mean
    stage_box(ax, 0.51, 0.55, 0.20, 0.18,
              "3. Family aggregation",
              "16 NCP_v3 families\nmean of program scores\nper donor",
              color="#FFF4E6", edgecolor=ACCENT)

    # Stage 4: Spearman + CI
    stage_box(ax, 0.03, 0.30, 0.20, 0.18,
              "4. Spearman \u03c1 vs duration",
              "AF-only fit, n = 12\nSR (n = 3) excluded from\nduration correlation;\nused as visual reference\n(month = 0) in 4d/e/f scatters.\nperm p (10,000) + bootstrap CI",
              color="#F2EAEA", edgecolor=LA_DEEP, body_align="left")

    # Stage 5: 3-class
    stage_box(ax, 0.27, 0.30, 0.20, 0.18,
              "5. Classify",
              "Strong  |\u03c1| \u2265 0.4\nWeak    0.2 \u2264 |\u03c1| < 0.4\nStable  |\u03c1| < 0.2\n\u2192 Fig 4c",
              color="#E8F0E8", edgecolor="#6E9B72", body_align="left")

    # Stage 6: Sensitivity matrix
    stage_box(ax, 0.51, 0.30, 0.20, 0.18,
              "6. Sensitivity matrix",
              "4 codings:\n  primary, gt1m_eq1, ordinal, LAPV\nsame-sign concordance n / 4\n\u2192 Fig 4b",
              color=LIGHT, edgecolor="#888", body_align="left")

    # Stage 7: case studies
    stage_box(ax, 0.75, 0.30, 0.20, 0.18,
              "7. Case studies",
              "Selection informed by Stage 6\nconcordance:\nFig 4d  F9  HS-GAG (4/4 sens)\nFig 4e  F2  Mito Ca\u00b2\u207a (anchor)\nFig 4f  F15 NTRK2 (cross-Fig 3)",
              color="#E8F0E8", edgecolor="#6E9B72", body_align="left")

    # arrows row 1 (Stage 1 -> 2 -> 3, then drop down to Stage 4)
    arrow(ax, 0.23, 0.64, 0.27, 0.64,
          label="all 16 donors\nincl. SR")
    arrow(ax, 0.47, 0.64, 0.51, 0.64)
    arrow(ax, 0.61, 0.55, 0.13, 0.48)
    # arrows row 2
    arrow(ax, 0.23, 0.39, 0.27, 0.39, label="3-class")
    arrow(ax, 0.47, 0.39, 0.51, 0.39, label="4 codings")
    # Stage 6 -> Fig 4b (heatmap output, downward into figure caption)
    arrow(ax, 0.61, 0.30, 0.61, 0.245, label="Fig 4b heatmap")
    # Stage 6 -> Stage 7 (case-study selection)
    arrow(ax, 0.71, 0.39, 0.75, 0.39,
          label="case selection\n(F9 4/4 etc.)")

    # Limitations banner
    ax.add_patch(mpatches.Rectangle((0.03, 0.04), 0.95, 0.20,
                                     facecolor="#FFF8F0", edgecolor=ACCENT, linewidth=0.6))
    ax.text(0.05, 0.215, "Methodological limitations (acknowledged in Discussion)",
            fontsize=7.0, fontweight="bold", color="#1A1A1A", va="top")
    ax.text(0.05, 0.180,
            "(i)  GSE41177 LAA AF n = 12 (donors with numeric duration); statistical power to detect |\u03c1| \u2265 0.4 (Spearman, \u03b1 = 0.05) is well below 80%.",
            fontsize=5.9, color="#333333", va="top")
    ax.text(0.05, 0.145,
            "     Fig 4 is therefore an effect-size-prioritised exploratory analysis, not a high-powered association screen.",
            fontsize=5.9, color="#333333", va="top")
    ax.text(0.05, 0.110,
            "(ii) No family reaches the strict criterion |\u03c1| \u2265 0.4 + p < 0.05; classification is by effect-size partition only (3-class above).",
            fontsize=5.9, color="#333333", va="top")
    ax.text(0.05, 0.075,
            "(iii) GSE41177 is single-cohort valvular AF (no cross-cohort meta-analysis at this duration resolution); LAPV is reported as a region sensitivity arm.",
            fontsize=5.9, color="#333333", va="top")
    ax.text(0.05, 0.040,
            "(iv) Fig 4 adopts a 3-class framework parallel to Fig 3 (cellular attribution) on the temporal axis.  The two are formally analogous, not biologically equivalent.",
            fontsize=5.9, color="#333333", va="top")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
