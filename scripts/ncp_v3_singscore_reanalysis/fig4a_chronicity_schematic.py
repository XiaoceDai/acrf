#!/usr/bin/env python3
"""Corrected Fig 4A — cross-sectional duration-analysis design schematic."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
OUT_BASE = FIGURES / "fig4a_chronicity_schematic"
WIDTH_MM, HEIGHT_MM = 178.0, 45.0

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
ACCENT = "#C8A076"
LIGHT = "#F2F2F2"


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def stage_box(ax, x, y, w, h, label_top, label_body, *, color=LIGHT, edgecolor="#888"):
    ax.add_patch(mpatches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.005,rounding_size=0.012",
        facecolor=color, edgecolor=edgecolor, linewidth=0.6,
    ))
    ax.text(x + w / 2.0, y + h - 0.03, label_top,
            ha="center", va="top", fontsize=6.2, fontweight="bold", color="#1A1A1A")
    ax.text(x + 0.01, y + h - 0.07, label_body,
            ha="left", va="top", fontsize=5.4, color="#333333")


def arrow(ax, x1, y1, x2, y2) -> None:
    ax.annotate(
        "", xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle="->", color="#1A1A1A", lw=0.6),
    )


def save_outputs(fig: plt.Figure) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BASE.with_suffix(".pdf"), bbox_inches=None, pad_inches=0.02)
    fig.savefig(OUT_BASE.with_suffix(".svg"), bbox_inches=None, pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    set_style()
    mm = 1.0 / 25.4
    fig = plt.figure(figsize=(WIDTH_MM * mm, HEIGHT_MM * mm))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.02, 0.94, "Fig 4a  Cross-sectional duration-analysis design (GSE41177)",
            fontsize=7.5, fontweight="bold", color="#1A1A1A")

    stage_box(ax, 0.02, 0.52, 0.15, 0.34,
              "Frozen NCP",
              "16 families\n30 programs\nfrozen LE sets",
              color="#EAF1F7", edgecolor=RA_DEEP)
    stage_box(ax, 0.19, 0.52, 0.15, 0.34,
              "GSE41177 bulk",
              "paired LAA + LAPV\nAffy HG-U133 Plus 2.0",
              color="#FFF4E6", edgecolor=ACCENT)
    stage_box(ax, 0.36, 0.52, 0.17, 0.34,
              "Exact-reference score",
              "python_exact_reference\nprogram-first family\nactivity score",
              color=LIGHT, edgecolor="#888")
    stage_box(ax, 0.55, 0.52, 0.17, 0.34,
              "Primary analysis",
              "12 AF-LAA donors\nfinite duration (mo)\nSR excluded from stats",
              color="#F2EAEA", edgecolor=LA_DEEP)
    stage_box(ax, 0.74, 0.52, 0.13, 0.34,
              "LAPV sensitivity",
              "regional ρ only\nno class assignment",
              color=LIGHT, edgecolor="#888")
    stage_box(ax, 0.89, 0.52, 0.09, 0.34,
              "Classify",
              "Strong / Weak / Stable\n|ρ_LAA| only\n→ Fig 4c",
              color="#E8F0E8", edgecolor="#6E9B72")

    for x1, x2 in [(0.17, 0.19), (0.34, 0.36), (0.53, 0.55), (0.72, 0.74), (0.87, 0.89)]:
        arrow(ax, x1, 0.69, x2, 0.69)
    arrow(ax, 0.81, 0.52, 0.81, 0.46)

    stage_box(ax, 0.55, 0.08, 0.43, 0.34,
              "Exemplar donors (Fig 4d)",
              "F9 HS-GAG / bone / GPCR\nF2 Mito Ca2+ transport\nF15 NTRK2 / TrkB\n(AF-LAA only; no SR anchors)",
              color="#E8F0E8", edgecolor="#6E9B72")
    arrow(ax, 0.935, 0.52, 0.76, 0.42)

    ax.add_patch(mpatches.Rectangle((0.02, 0.06), 0.50, 0.36,
                                     facecolor="#FFF8F0", edgecolor=ACCENT, linewidth=0.5))
    ax.text(0.03, 0.38, "Contract (reader-facing)", fontsize=5.8, fontweight="bold", color="#1A1A1A", va="top")
    notes = [
        "Family activity score (not aCAIrel).",
        "Descriptive chronicity proxy — not recurrence, progression, or treatment response.",
        "SR excluded from all duration statistics.",
        "LAPV ρ reported for regional sensitivity; does not determine chronicity class.",
        "n = 12 primary AF-LAA; exploratory effect-size partition.",
    ]
    for i, line in enumerate(notes):
        ax.text(0.03, 0.33 - i * 0.055, f"• {line}", fontsize=5.2, color="#333333", va="top")

    save_outputs(fig)
    print(f"wrote {OUT_BASE.with_suffix('.pdf').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
