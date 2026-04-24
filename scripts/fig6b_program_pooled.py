#!/usr/bin/env python3
"""
Figure 6B — Program-level pooled summary.

Compact "dot + CI" summary of the 10 frozen panel programs, rows ordered by
`program_order` (Fig 1C / Fig 2 / Fig 5A order).
Each row carries:
  - axis-membership color (left side strip + dot color):
        STRUCTURAL → #B35806 (PuOr +)   METABOLIC → #542788 (PuOr -)
        panel-only → #BDBDBD
  - horizontal CI line + dot (point size scales with 1/SE)
  - right-hand g [CI] + I²

Reads `results/figure6_meta_program.tsv` (pooled rows only).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

from plot_style import COL_DOUBLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "results" / "figure6_meta_program.tsv"
OUT = ROOT / "figures" / "fig6b_program_pooled.pdf"

AXIS_COL = {
    "STRUCTURAL": "#B35806",
    "METABOLIC":  "#542788",
    "panel-only": "#BDBDBD",
}


def _round(ax, x, y, w, h, fc, ec="none", lw=0.4, z=1, rounding=0.012):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    ))


def main() -> None:
    set_nature_style()
    df = pd.read_csv(IN, sep="\t")
    df = df[df["type"] == "pooled"].sort_values("program_order").reset_index(drop=True)
    n = len(df)

    # ---- layout / geometry ----------------------------------------------
    fig_w = COL_DOUBLE * 0.95
    row_h = 0.36
    fig_h = row_h * n + 1.55
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.275, right=0.965, top=0.84, bottom=0.16)

    x_min = min(df["ci_lo"].min(), -1.0) - 0.40
    x_max = max(df["ci_hi"].max(), +1.0) + 0.85   # extra room for right text
    plot_max_for_text = x_max - 0.05

    # ---- background card ------------------------------------------------
    _round(ax, x_min, -0.6, x_max - x_min, n + 0.2,
           fc="#FAFAFA", ec="#BBBBBB", lw=0.45, z=0, rounding=0.014)

    # ---- alternating stripes + axis-membership side strip ---------------
    strip_w = (x_max - x_min) * 0.012
    for i, r in df.iterrows():
        # alternating zebra
        if i % 2 == 1:
            ax.add_patch(Rectangle((x_min + 0.02, i - 0.42),
                                   x_max - x_min - 0.04, 0.84,
                                   facecolor="white", edgecolor="none",
                                   zorder=1))
        # left side strip with rounded ends
        _round(ax, x_min + 0.04, i - 0.34,
               strip_w * 1.6, 0.68,
               fc=AXIS_COL[r["axis"]], ec="white", lw=0.3,
               z=2, rounding=0.012)

    # ---- LA-high / RA-high separator (between row 4 and 5) -------------
    ax.plot([x_min + 0.04, x_max - 0.04], [4.5, 4.5],
            color="#CCCCCC", lw=0.45, zorder=2)

    # ---- 0 line ---------------------------------------------------------
    ax.axvline(0, color="#888888", lw=0.55, ls=(0, (1.5, 1.5)), zorder=2)

    # ---- marker scaling -------------------------------------------------
    se = df["se_g"].astype(float).to_numpy()
    inv_se = 1.0 / np.maximum(se, 1e-6)
    inv_se_norm = (inv_se - inv_se.min()) / max((inv_se.max() - inv_se.min()), 1e-9)

    # ---- forest dots + right-hand annotations --------------------------
    for i, r in df.iterrows():
        y = i
        col = AXIS_COL[r["axis"]]
        ax.plot([r["ci_lo"], r["ci_hi"]], [y, y],
                color=col, lw=1.0, alpha=0.95, zorder=3,
                solid_capstyle="round")
        ax.scatter([r["g"]], [y], s=24 + 60 * inv_se_norm[i],
                   facecolor=col, edgecolor="white", linewidth=0.55,
                   zorder=4)
        # g [CI]   I²
        ax.text(plot_max_for_text, y,
                f"{r['g']:+.2f}  [{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}]   "
                f"I²={r['I2']:.0f}%",
                fontsize=5.0, ha="right", va="center",
                color="#1A1A1A", family="monospace")

    # ---- y axis labels (program names) ---------------------------------
    ax.set_yticks(range(n))
    ax.set_yticklabels(df["program_label"].tolist(), fontsize=6.0)
    for tick, r in zip(ax.get_yticklabels(), df.itertuples(index=False)):
        if r.axis in ("STRUCTURAL", "METABOLIC"):
            tick.set_fontweight("bold")
            tick.set_color("#1A1A1A")
        else:
            tick.set_color("#555555")

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(n - 0.5, -0.6)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.tick_params(axis="x", labelsize=5.6, pad=2, length=2)
    ax.set_xlabel("Hedges' g (AF − SR), pooled (REML + Hartung–Knapp)",
                  fontsize=6.4, labelpad=4)

    # ---- top-right legend (axis membership) ----------------------------
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="s", linestyle="None", markersize=5.5,
               markerfacecolor=AXIS_COL["STRUCTURAL"],
               markeredgecolor="white", markeredgewidth=0.4,
               label="STRUCTURAL support"),
        Line2D([], [], marker="s", linestyle="None", markersize=5.5,
               markerfacecolor=AXIS_COL["METABOLIC"],
               markeredgecolor="white", markeredgewidth=0.4,
               label="METABOLIC support"),
        Line2D([], [], marker="s", linestyle="None", markersize=5.5,
               markerfacecolor=AXIS_COL["panel-only"],
               markeredgecolor="white", markeredgewidth=0.4,
               label="panel-only"),
    ]
    fig.legend(handles=handles, loc="upper right",
               bbox_to_anchor=(0.965, 0.985),
               fontsize=5.2, frameon=False, ncol=3,
               handlelength=1.0, handleheight=0.55, columnspacing=1.0,
               labelspacing=0.25)

    fig.suptitle(
        "Structural programs show the most consistent pooled attenuation",
        fontsize=7.4, y=0.95, fontweight="bold",
    )
    fig.text(0.012, 0.965, "b",
             fontsize=10, fontweight="bold", ha="left", va="top",
             transform=fig.transFigure)

    fig.text(0.5, 0.020,
             "Marker size scales with 1/SE; bars = 95% CI (Hartung–Knapp).  "
             "Pooled estimates k = 3 cohorts, read as consistency check.",
             ha="center", va="bottom", fontsize=4.8,
             color="#666666", style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
