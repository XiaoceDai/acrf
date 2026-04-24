#!/usr/bin/env python3
"""
Figure 6A — Axis-level random-effects forest (2 axes × 3 cohorts).

Two stacked sub-panels (STRUCTURAL above METABOLIC); each sub-panel:
  - 3 cohort rows: cohort name + n_AF/n_SR + g + 95% CI line + ●
                   (point size ∝ inverse-variance weight)
  - 1 RE pooled row (HK): black diamond + numeric annotation
  - Right-hand text: g [CI]
  - Side annotation: pooled g, 95% CI, I², τ², p_HK

Sign convention (locked):
  Hedges' g (AF − SR), donor-level axis composite.
  g < 0  ⇒  AF lower than SR  ⇒  AF attenuates the healthy axis.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle

from plot_style import COL_DOUBLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "results" / "figure6_meta_axis.tsv"
OUT = ROOT / "figures" / "fig6a_axis_forest.pdf"

# Background tints (Fig 2 PuOr semantics, very pale)
TINT_STRUCT = "#ECF4E8"
TINT_METAB = "#FCF1DD"
HEAD_STRUCT = "#D6E8CC"
HEAD_METAB = "#F5DDB0"
COL_STRUCT_LINE = "#3F6512"
COL_METAB_LINE = "#A55F00"
COL_NEG = "#542788"   # PuOr negative (= attenuation / reversal)
COL_POS = "#B35806"   # PuOr positive (= preservation / amplification)
COL_POOL = "#1A1A1A"


def _round(ax, x, y, w, h, fc, ec="none", lw=0.4, z=1, rounding=0.012):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    ))


def _draw_subpanel(ax, sub_meta, axis_label, axis_tint, head_color,
                   axis_line_color, x_min, x_max,
                   left_lab_x, right_lab_x):
    rows = sub_meta.reset_index(drop=True)
    n_rows = len(rows)
    head_h = 0.55

    # Background tint for the whole sub-panel (rows + small head band)
    _round(ax, x_min, -0.5, x_max - x_min, n_rows + head_h,
           fc=axis_tint, ec="#BBBBBB", lw=0.45, z=0, rounding=0.018)
    # Header band
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

    # 0 line
    ax.axvline(0, color="#888888", lw=0.55, ls=(0, (1.5, 1.5)), zorder=2)

    max_w = float(np.nanmax(1.0 / rows.loc[rows["type"] == "cohort", "var_g"]))
    for i, r in rows.iterrows():
        y = i
        # row stripe
        if i % 2 == 1 and r["type"] == "cohort":
            ax.add_patch(Rectangle((x_min + 0.02, y - 0.42),
                                   x_max - x_min - 0.04, 0.84,
                                   facecolor="white", edgecolor="none",
                                   alpha=0.45, zorder=1))

        if r["type"] == "cohort":
            w = (1.0 / r["var_g"]) / max_w
            ms = 28 + 70 * w
            color = COL_NEG if r["g"] < 0 else COL_POS
            ax.plot([r["ci_lo"], r["ci_hi"]], [y, y],
                    color=color, lw=1.0, alpha=0.95, zorder=3,
                    solid_capstyle="round")
            ax.scatter([r["g"]], [y], s=ms,
                       facecolor=color, edgecolor="white", linewidth=0.55,
                       zorder=4)
            ax.text(left_lab_x, y, f"{r['cohort']}",
                    fontsize=5.9, ha="left", va="center", color="#1A1A1A",
                    fontweight="bold")
            ax.text(left_lab_x, y + 0.30,
                    f"AF n={int(r['n_AF'])}  ·  SR n={int(r['n_SR'])}",
                    fontsize=4.7, ha="left", va="bottom", color="#666666")
            ax.text(right_lab_x, y,
                    f"{r['g']:+.2f}  [{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}]",
                    fontsize=5.4, ha="right", va="center",
                    color="#1A1A1A", family="monospace")
        else:
            # divider above
            ax.plot([x_min + 0.04, x_max - 0.04], [y - 0.5, y - 0.5],
                    color="#BBBBBB", lw=0.5, zorder=3)
            poly_x = [r["ci_lo"], r["g"], r["ci_hi"], r["g"]]
            poly_y = [y, y - 0.24, y, y + 0.24]
            ax.add_patch(Polygon(list(zip(poly_x, poly_y)),
                                 facecolor=COL_POOL, edgecolor=COL_POOL,
                                 linewidth=0.5, zorder=4))
            ax.text(left_lab_x, y, "RE pooled (Hartung–Knapp)",
                    fontsize=5.9, ha="left", va="center", color=COL_POOL,
                    fontweight="bold")
            ax.text(right_lab_x, y,
                    f"{r['g']:+.2f}  [{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}]",
                    fontsize=5.4, ha="right", va="center",
                    color=COL_POOL, family="monospace", fontweight="bold")
            ax.text(right_lab_x, y + 0.40,
                    f"I² = {r['I2']:.0f}%   τ² = {r['tau2']:.2f}"
                    f"   p$_\\mathrm{{HK}}$ = {r['p_HK']:.2f}",
                    fontsize=4.8, ha="right", va="center", color="#444444")

    ax.set_ylim(n_rows + head_h - 0.5, -0.5)
    ax.set_xlim(x_min, x_max)
    ax.set_yticks([])
    for s in ("top", "left", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", labelsize=5.6, pad=2, length=2)


def main() -> None:
    set_nature_style()
    meta = pd.read_csv(IN, sep="\t")

    all_lo = meta["ci_lo"].astype(float).min()
    all_hi = meta["ci_hi"].astype(float).max()
    pad = 0.5
    plot_lo = float(min(all_lo, -1.0)) - pad
    plot_hi = float(max(all_hi, +1.0)) + pad
    # Reserve space on either side for left labels and right text
    LEFT_LAB_PAD  = 1.7
    RIGHT_LAB_PAD = 1.7
    x_min = plot_lo - LEFT_LAB_PAD
    x_max = plot_hi + RIGHT_LAB_PAD
    left_lab_x  = plot_lo - LEFT_LAB_PAD + 0.06
    right_lab_x = plot_hi + RIGHT_LAB_PAD - 0.06

    fig_w = COL_DOUBLE * 0.95
    fig_h = 3.30
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.025, right=0.975, top=0.85, bottom=0.13,
                        hspace=0.55)

    ax_top = fig.add_subplot(2, 1, 1)
    ax_bot = fig.add_subplot(2, 1, 2, sharex=ax_top)

    s = meta[meta["axis"] == "STRUCTURAL"].copy()
    s = s.sort_values(by=["type", "cohort"], key=lambda c: c.map(
        {"cohort": 0, "pooled": 1}) if c.name == "type" else c)
    m = meta[meta["axis"] == "METABOLIC"].copy()
    m = m.sort_values(by=["type", "cohort"], key=lambda c: c.map(
        {"cohort": 0, "pooled": 1}) if c.name == "type" else c)

    _draw_subpanel(ax_top, s, "STRUCTURAL axis  (matrix / proteoglycan)",
                   TINT_STRUCT, HEAD_STRUCT, COL_STRUCT_LINE,
                   x_min, x_max, left_lab_x, right_lab_x)
    _draw_subpanel(ax_bot, m, "METABOLIC axis  (OXPHOS / translation)",
                   TINT_METAB, HEAD_METAB, COL_METAB_LINE,
                   x_min, x_max, left_lab_x, right_lab_x)

    # x-axis only on bottom sub-panel
    ax_bot.tick_params(axis="x", labelsize=5.6, pad=2, length=2,
                       labelbottom=True)
    ax_top.tick_params(axis="x", labelbottom=False, length=0)
    ax_bot.spines["bottom"].set_visible(True)
    ax_bot.spines["bottom"].set_linewidth(0.5)

    ax_bot.set_xlabel(
        "Hedges' g  (AF − SR), donor-level axis composite",
        fontsize=6.4, labelpad=4,
    )

    ax_bot.text(0.0, 1.07,
                "← AF attenuates healthy axis        AF preserves / amplifies →",
                transform=ax_bot.transAxes,
                fontsize=4.9, color="#666666", ha="left", va="bottom",
                style="italic")

    fig.suptitle(
        "Cross-cohort consistency of dual-axis chamber rewiring",
        fontsize=7.6, y=0.965, fontweight="bold",
    )
    fig.text(0.012, 0.97, "a",
             fontsize=10, fontweight="bold", ha="left", va="top",
             transform=fig.transFigure)

    fig.text(0.5, 0.012,
             "k = 3 paired-chamber cohorts; pooled estimates use the "
             "Hartung–Knapp variance adjustment.  Read as a consistency "
             "check, not a high-powered pooled discovery.",
             ha="center", va="bottom", fontsize=4.8,
             color="#666666", style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
