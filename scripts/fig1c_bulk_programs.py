#!/usr/bin/env python3
"""
Figure 1C - Bulk chamber programs overview.

One-column heatmap:
    rows = 10 stable bulk chamber programs (normal_chamber_programs.tsv)
    col  = "Bulk chamber bias"
    val  = signed NES (positive if LA-biased, negative if RA-biased)
    cmap = diverging RdBu    (blue = LA-high, red = RA-high)
            -- matches the LA=blue / RA=red convention used in 1A / 1B.

Purpose: lift 1B's gene-level LA/RA differences up to the program level.

Output: figures/fig1c_bulk_programs.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

from plot_style import COL_SINGLE, LA_COLOR, RA_COLOR, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "figures" / "fig1c_bulk_programs.pdf"


def load_bulk_programs() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "normal_chamber_programs.tsv", sep="\t")
    for c in ("program", "program_id", "direction"):
        df[c] = df[c].astype(str).str.strip('"')
    # Bulk NES in this file is already signed (LA-biased positive, RA-biased negative).
    df["signed_NES"] = df["NES"]
    df = df.sort_values("signed_NES", ascending=False).reset_index(drop=True)
    return df


def main() -> None:
    set_nature_style()
    df = load_bulk_programs()

    n = len(df)
    fig_w = COL_SINGLE * 0.85
    fig_h = 0.28 * n + 0.9
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    # Lock the plot area vertically so rows align with 1D panel.
    fig.subplots_adjust(left=0.62, right=0.92, top=0.92, bottom=0.24)

    values = df["signed_NES"].to_numpy()
    vneg = float(values[values < 0].min()) if (values < 0).any() else -1e-6
    vpos = float(values[values > 0].max()) if (values > 0).any() else +1e-6
    norm = TwoSlopeNorm(vcenter=0.0, vmin=vneg, vmax=vpos)
    # Convention across Figure 1: LA = blue, RA = red.
    # LA-biased pathways have positive NES -> blue end; flipped from RdBu_r.
    cmap = plt.get_cmap("RdBu")

    # Draw each cell as a VECTOR Rectangle (not an imshow raster image).
    # Illustrator handles rectangles with plain RGB fills reliably; embedded
    # imshow rasters occasionally render with a wrong colorspace.
    for i, v in enumerate(values):
        rgba = cmap(norm(v))
        ax.add_patch(Rectangle(
            (-0.5, i - 0.5), 1.0, 1.0,
            facecolor=rgba, edgecolor="white", linewidth=0.6,
        ))

    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_aspect("auto")

    # Annotate each cell with the signed bulk NES.
    for i, nes in enumerate(df["NES"].to_numpy()):
        v = values[i]
        side_max = vpos if v >= 0 else abs(vneg)
        txt_color = "white" if abs(v) > 0.55 * side_max else "black"
        ax.text(0, i, f"{nes:.2f}",
                ha="center", va="center", fontsize=6,
                color=txt_color, fontweight="bold")

    ax.set_xticks([0])
    ax.set_xticklabels(["Bulk chamber bias"])
    ax.set_yticks(range(n))
    ax.set_yticklabels(df["program"].tolist())
    ax.tick_params(axis="both", length=0)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)

    # Color row labels by bulk direction, same key as 1D row labels.
    for tick_label, direction in zip(ax.get_yticklabels(), df["direction"]):
        tick_label.set_color(LA_COLOR if direction.startswith("LA") else RA_COLOR)

    ax.set_title(
        "Normal Chamber Prior: LA-biased metabolic and RA-biased matrix/glycan programs",
        fontsize=7.4, pad=4)

    # Colorbar placed at the same vertical band as 1D.  Shrink the bar itself
    # to leave room for "RA-high" / "LA-high" labels on its left / right.
    cbar_y, cbar_h = 0.09, 0.025
    cbar_x, cbar_w = 0.68, 0.20
    cax = fig.add_axes([cbar_x, cbar_y, cbar_w, cbar_h])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(width=0.4, length=2, pad=1, labelsize=6)

    # Flank labels: RA-high on the LEFT end (negative NES -> red),
    # LA-high on the RIGHT end (positive NES -> blue).
    label_y = cbar_y + cbar_h / 2
    fig.text(cbar_x - 0.015, label_y, "RA-high",
             color=RA_COLOR, fontsize=6, fontweight="bold",
             ha="right", va="center")
    fig.text(cbar_x + cbar_w + 0.015, label_y, "LA-high",
             color=LA_COLOR, fontsize=6, fontweight="bold",
             ha="left", va="center")
    fig.text(cbar_x + cbar_w / 2, cbar_y - 0.04, "signed bulk NES",
             fontsize=6, ha="center", va="top")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
