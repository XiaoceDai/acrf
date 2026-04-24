#!/usr/bin/env python3
"""
Figure S7 — Leave-one-cohort-out (LOCO) sensitivity for the axis-level
random-effects meta of Figure 6A.

Reads `results/figure6_meta_axis_loco.tsv` (written by h2_meta_axis.R) and
renders, per axis, the LOCO pooled estimate when each of the 3 cohorts is
excluded.  Confirms direction stability of the structural axis.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle

from plot_style import COL_ONEHALF, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "results" / "figure6_meta_axis_loco.tsv"
OUT = ROOT / "figures" / "figS7_loco_sensitivity.pdf"

TINT_STRUCT = "#ECF4E8"
TINT_METAB = "#FCF1DD"
HEAD_STRUCT = "#D6E8CC"
HEAD_METAB = "#F5DDB0"
COL_STRUCT_LINE = "#3F6512"
COL_METAB_LINE = "#A55F00"
COL_NEG = "#542788"
COL_POS = "#B35806"
COL_POOL = "#1A1A1A"


def _round(ax, x, y, w, h, fc, ec="none", lw=0.4, z=1, rounding=0.012):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    ))


def _draw_axis(ax, sub, axis_label, axis_tint, head_color, axis_line_color,
               x_min, x_max, left_lab_x, right_lab_x):
    n_rows = len(sub)
    head_h = 0.55
    _round(ax, x_min, -0.5, x_max - x_min, n_rows + head_h,
           fc=axis_tint, ec="#BBBBBB", lw=0.45, z=0, rounding=0.018)
    _round(ax, x_min, n_rows - 0.5 - 0.04,
           x_max - x_min, head_h + 0.04,
           fc=head_color, ec="none", z=1, rounding=0.018)
    ax.text(left_lab_x, n_rows - 0.5 + head_h / 2, axis_label,
            fontsize=6.4, fontweight="bold", color=axis_line_color,
            ha="left", va="center")
    ax.text(right_lab_x, n_rows - 0.5 + head_h / 2,
            "g (AF − SR)   95% CI",
            fontsize=5.4, color="#444444",
            ha="right", va="center", family="monospace")

    ax.axvline(0, color="#888888", lw=0.55, ls=(0, (1.5, 1.5)), zorder=2)

    for i, r in sub.iterrows():
        y = i
        if i % 2 == 1:
            ax.add_patch(Rectangle((x_min + 0.02, y - 0.42),
                                   x_max - x_min - 0.04, 0.84,
                                   facecolor="white", edgecolor="none",
                                   alpha=0.45, zorder=1))
        color = COL_NEG if r["g"] < 0 else COL_POS
        ax.plot([r["ci_lo"], r["ci_hi"]], [y, y],
                color=color, lw=1.0, alpha=0.95, zorder=3,
                solid_capstyle="round")
        # diamond marker (LOCO is a pooled estimate too)
        poly_x = [r["ci_lo"], r["g"], r["ci_hi"], r["g"]]
        poly_y = [y, y - 0.20, y, y + 0.20]
        ax.add_patch(Polygon(list(zip(poly_x, poly_y)),
                             facecolor=color, edgecolor=color,
                             linewidth=0.5, zorder=4))
        ax.text(left_lab_x, y, f"−   excl.  {r['excluded_cohort']}",
                fontsize=5.7, ha="left", va="center", color="#1A1A1A",
                fontweight="bold")
        ax.text(right_lab_x, y,
                f"{r['g']:+.2f}  [{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}]",
                fontsize=5.4, ha="right", va="center",
                color="#1A1A1A", family="monospace")
        ax.text(right_lab_x, y + 0.32,
                f"I² = {r['I2']:.0f}%   τ² = {r['tau2']:.2f}",
                fontsize=4.6, ha="right", va="center", color="#444444")

    ax.set_ylim(n_rows + head_h - 0.5, -0.5)
    ax.set_xlim(x_min, x_max)
    ax.set_yticks([])
    for s in ("top", "left", "right", "bottom"):
        ax.spines[s].set_visible(False)


def main() -> None:
    set_nature_style()
    loco = pd.read_csv(IN, sep="\t")

    all_lo = loco["ci_lo"].astype(float).min()
    all_hi = loco["ci_hi"].astype(float).max()
    pad = 0.5
    plot_lo = float(min(all_lo, -1.0)) - pad
    plot_hi = float(max(all_hi, +1.0)) + pad
    LEFT_PAD = 1.7; RIGHT_PAD = 1.7
    x_min = plot_lo - LEFT_PAD; x_max = plot_hi + RIGHT_PAD
    left_lab_x = plot_lo - LEFT_PAD + 0.06
    right_lab_x = plot_hi + RIGHT_PAD - 0.06

    fig_w = COL_ONEHALF * 1.4
    fig_h = 2.85
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.025, right=0.975, top=0.85, bottom=0.13,
                        hspace=0.55)

    ax_top = fig.add_subplot(2, 1, 1)
    ax_bot = fig.add_subplot(2, 1, 2, sharex=ax_top)

    s = loco[loco["axis"] == "STRUCTURAL"].sort_values("excluded_cohort").reset_index(drop=True)
    m = loco[loco["axis"] == "METABOLIC"].sort_values("excluded_cohort").reset_index(drop=True)

    _draw_axis(ax_top, s, "STRUCTURAL axis  (LOCO)",
               TINT_STRUCT, HEAD_STRUCT, COL_STRUCT_LINE,
               x_min, x_max, left_lab_x, right_lab_x)
    _draw_axis(ax_bot, m, "METABOLIC axis  (LOCO)",
               TINT_METAB, HEAD_METAB, COL_METAB_LINE,
               x_min, x_max, left_lab_x, right_lab_x)

    ax_bot.tick_params(axis="x", labelsize=5.6, pad=2, length=2,
                       labelbottom=True)
    ax_top.tick_params(axis="x", labelbottom=False, length=0)
    ax_bot.spines["bottom"].set_visible(True)
    ax_bot.spines["bottom"].set_linewidth(0.5)
    ax_bot.set_xlabel("Hedges' g  (AF − SR), LOCO pooled estimate",
                      fontsize=6.4, labelpad=4)

    fig.suptitle(
        "Leave-one-cohort-out sensitivity for axis-level meta-analysis",
        fontsize=7.4, y=0.965, fontweight="bold",
    )
    fig.text(0.5, 0.012,
             "Supplementary to Fig 6A.  Each row excludes the named cohort "
             "and re-runs the REML/Hartung–Knapp pool on the remaining two.",
             ha="center", va="bottom", fontsize=4.8,
             color="#666666", style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
