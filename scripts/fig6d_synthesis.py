#!/usr/bin/env python3
"""
Figure 6D — Robustness-to-interpretation synthesis (closing panel).

Three-input → one-output synthesis of Figure 6:
  6A  cross-cohort axis-level meta  →  axis-level robustness  ✓
  6B  per-program pooled summary    →  ECM/HS-GAG and OXPHOS anchor the axes
  6C  pre-registered MR / coloc     →  single-gene mediation NOT supported
            ↓ ↓ ↓
       Distributed program-level AF remodeling
       (not master-regulator; consistent with Fig 3 carrier pattern)

This panel reuses the colour and shape discipline of 6A / 6B / 6C; no new
data is plotted.  The previous "exposome upstream framing" is moved to the
Graphical Abstract / Discussion schematic (narrative companion to meta + MR panels).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from plot_style import COL_DOUBLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "fig6d_synthesis.pdf"

# Card colours (echo 6A header tints + 6C verdict palette)
COL_OK   = "#1B7837"    # supports / robust
COL_WARN = "#9E9E9E"    # qualifier / no claim
COL_NO   = "#B2182B"    # not supported

CARD_BG_A = "#ECF4E8"   # robustness (axis meta)
CARD_BG_B = "#F4ECEC"   # pooled programs
CARD_BG_C = "#F4ECEC"   # MR (negative)
HEAD_A = "#D6E8CC"
HEAD_B = "#F2DFE0"
HEAD_C = "#F2DFE0"
ICON_A = COL_OK
ICON_B = COL_WARN
ICON_C = COL_NO


def _round(ax, x, y, w, h, fc, ec="none", lw=0.5, z=1, rounding=0.014):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    ))


def _arrow(ax, p1, p2, color="#444444", lw=1.0, mut=8, curve=0.0,
           dashed=False, alpha=0.95):
    style = (0, (3, 2)) if dashed else "-"
    ax.add_patch(FancyArrowPatch(
        p1, p2,
        arrowstyle="-|>,head_length=3.0,head_width=1.8",
        connectionstyle=f"arc3,rad={curve}",
        mutation_scale=mut, lw=lw, color=color, linestyle=style,
        capstyle="round", joinstyle="round",
        shrinkA=2, shrinkB=2, alpha=alpha, zorder=10,
    ))


def _input_card(ax, x, w, head_color, body_color, head_label,
                head_color_text, takeaway_text, takeaway_color,
                icon_color, fig_tag, status_label):
    """One column input card with a coloured head band, central takeaway,
    a small status chip + Fig reference at bottom corners."""
    # Card body
    _round(ax, x, 0.34, w, 0.42, fc=body_color, ec="#BBBBBB",
           lw=0.5, z=1, rounding=0.014)
    # Head band
    _round(ax, x, 0.65, w, 0.11, fc=head_color, ec="none",
           z=2, rounding=0.014)
    ax.text(x + w / 2, 0.705, head_label,
            ha="center", va="center", fontsize=6.4, fontweight="bold",
            color=head_color_text)

    # Takeaway (multi-line, italic, in card body)
    ax.text(x + w / 2, 0.555, takeaway_text,
            ha="center", va="center", fontsize=5.7,
            color=takeaway_color, style="italic", linespacing=1.45)

    # Status chip (top-right): small coloured pill instead of glyph
    chip_w, chip_h = 0.075, 0.030
    _round(ax, x + w - chip_w - 0.012, 0.715,
           chip_w, chip_h,
           fc=icon_color, ec="white", lw=0.4,
           z=4, rounding=0.010)
    ax.text(x + w - chip_w / 2 - 0.012, 0.730, status_label,
            ha="center", va="center", fontsize=4.8,
            color="white", fontweight="bold", zorder=5)

    # Fig reference (bottom-left of card)
    ax.text(x + 0.020, 0.355, fig_tag,
            ha="left", va="bottom", fontsize=4.7, color="#888888",
            style="italic", zorder=4)


def main() -> None:
    set_nature_style()
    fig_w = COL_DOUBLE
    fig_h = 2.85
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # -------- Three input cards --------
    cards = [
        dict(x=0.04, w=0.29,
             head_color=HEAD_A, body_color=CARD_BG_A,
             head_color_text="#1F4A12",
             head_label="Cross-cohort axis meta",
             takeaway_text=(
                 "Axis-level Hedges' g pools negatively\n"
                 "across 3 paired-chamber cohorts;\n"
                 "structural axis I² = 0%, metabolic\n"
                 "axis directional but heterogeneous."),
             takeaway_color="#1A1A1A",
             icon_color=ICON_A, fig_tag="Fig 6A",
             status_label="ROBUST"),
        dict(x=0.355, w=0.29,
             head_color=HEAD_B, body_color="#FAF6F1",
             head_color_text="#5A3F00",
             head_label="Per-program pooled summary",
             takeaway_text=(
                 "All 10 frozen programs pool g ≤ 0;\n"
                 "ECM/HS-GAG and OXPHOS-ETC clusters\n"
                 "carry the lowest I² and the largest\n"
                 "effect — they anchor the two axes."),
             takeaway_color="#1A1A1A",
             icon_color="#B35806", fig_tag="Fig 6B",
             status_label="ANCHOR"),
        dict(x=0.67, w=0.29,
             head_color=HEAD_C, body_color=CARD_BG_C,
             head_color_text="#5C0F0F",
             head_label="Pre-registered cis-eQTL × AF GWAS",
             takeaway_text=(
                 "All 4 anchors fail the cis-eQTL gate\n"
                 "in tissue-matched GTEx; coloc PP4 < 0.05.\n"
                 "Single-gene mediation is NOT supported;\n"
                 "anchors not replaced (spec §6.6)."),
             takeaway_color="#1A1A1A",
             icon_color=ICON_C, fig_tag="Fig 6C",
             status_label="NOT SUPP."),
    ]
    for c in cards:
        _input_card(ax, c["x"], c["w"],
                    c["head_color"], c["body_color"],
                    c["head_label"], c["head_color_text"],
                    c["takeaway_text"], c["takeaway_color"],
                    c["icon_color"], c["fig_tag"], c["status_label"])

    # -------- Three converging arrows --------
    # Endpoints aimed at the top edge of the synthesis box
    targets_y = 0.215
    _arrow(ax, (0.185, 0.34), (0.350, targets_y),
           color=ICON_A, lw=1.0, mut=9, curve=+0.06)
    _arrow(ax, (0.500, 0.34), (0.500, targets_y),
           color="#B35806", lw=1.0, mut=9, curve=0.0)
    _arrow(ax, (0.815, 0.34), (0.650, targets_y),
           color=ICON_C, lw=1.0, mut=9, curve=-0.06)

    # -------- Synthesis output card --------
    _round(ax, 0.20, 0.05, 0.60, 0.16, fc="#EFEFEF", ec="#444444",
           lw=0.85, z=1, rounding=0.018)
    ax.text(0.50, 0.16,
            "Distributed program-level AF remodeling",
            ha="center", va="center", fontsize=7.2, fontweight="bold",
            color="#1A1A1A")
    ax.text(0.50, 0.115,
            "axis-level robust  ·  no master cis-eQTL driver",
            ha="center", va="center", fontsize=5.6,
            color="#444444")
    ax.text(0.50, 0.078,
            "A chamber-aware, program-level substrate framework — not a master-gene model — "
            "best explains the observed AF remodeling.",
            ha="center", va="center", fontsize=4.7,
            color="#666666", style="italic", linespacing=1.4)

    # -------- Title + label --------
    fig.suptitle(
        "Distributed program-level remodeling rather than a master-gene model",
        fontsize=7.4, y=0.965, fontweight="bold",
    )
    fig.text(0.012, 0.97, "d",
             fontsize=10, fontweight="bold", ha="left", va="top",
             transform=fig.transFigure)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
