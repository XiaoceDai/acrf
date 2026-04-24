#!/usr/bin/env python3
"""
Figure 2C - program x cohort rewiring class matrix + consensus + 115574 replication.

Output:

  figures/fig2c_rewiring_classes_with_trend.pdf
    strict primary fill (§3.5), trend overlay as corner triangles on
    the two discovery columns only (§3.5b).

Columns: 128188, 79768, Consensus, 115574
Rows:    10 programs, Fig 1C order.
Cells:
  * fill colour = strict rewiring_class (or consensus_class for Consensus col,
                  or 115574's strict class for 115574 col)
  * if variant "with_trend": bottom-right triangle coloured by trend_class
    for the discovery cohorts whenever trend != strict
  * text markers: *  p_perm < 0.05
                  .  0.05 <= p_perm < 0.10
                  (applied to 128188 / 79768 / 115574; no marker in Consensus)
  * Reversed-weak gets a diagonal-line hatch on top of the Attenuated light-purple fill
  * Consensus "NA" cells are shown as white with a faint grey border
  * Low-support cells are grey (legend on the right)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon, Rectangle

from plot_style import (
    COL_ONEHALF, LA_COLOR, RA_COLOR, REWIRING_CLASS_COLORS,
    set_nature_style,
)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "fig2c_rewiring_classes_with_trend.pdf"

DISCOVERY = ["GSE128188", "GSE79768"]
COHORT_SHORT = {"GSE128188": "128188", "GSE79768": "79768", "GSE115574": "115574"}

# Order in which classes stack in the legend.
LEGEND_ORDER = [
    "Amplified", "Amplified-weak",
    "Preserved",
    "Attenuated-weak", "Attenuated",
    "Reversed-weak",   "Reversed",
    "Low-support",
    "NA",
]


def load_data() -> pd.DataFrame:
    return pd.read_csv(ROOT / "results" / "af_chamber_rewiring_programs.tsv",
                       sep="\t")


def pick_class(row, cohort: str, which: str) -> str:
    """Return strict/trend/consensus class label for a given column."""
    if cohort == "Consensus":
        v = row.get("consensus_class")
        if pd.isna(v):
            return "NA"
        return v
    col = f"{which}_class_{cohort}"
    return row[col]


def perm_marker(p) -> str:
    if pd.isna(p):
        return ""
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "."
    return ""


def draw_cell(ax, j, i, strict_cls, trend_cls=None, text=""):
    """Draw one cell: fill by strict class, optional trend corner, text marker."""
    fill = REWIRING_CLASS_COLORS.get(strict_cls, "#FFFFFF")
    if strict_cls == "NA":
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                               facecolor="white", edgecolor="#BBBBBB", linewidth=0.5))
    else:
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                               facecolor=fill, edgecolor="white", linewidth=0.5))
    # Reversed-weak gets a diagonal hatch to distinguish from Attenuated (same light purple).
    if strict_cls == "Reversed-weak":
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                               facecolor="none", edgecolor="#542788",
                               linewidth=0.0, hatch="////"))

    # Trend overlay: bottom-right triangle, only where trend differs and both valid.
    if trend_cls is not None and trend_cls not in (None, "", strict_cls):
        if trend_cls == "NA":
            pass
        else:
            tri = np.array([[j + 0.5, i + 0.5],
                            [j + 0.5, i - 0.05],
                            [j - 0.05, i + 0.5]])
            tri_color = REWIRING_CLASS_COLORS.get(trend_cls, "#FFFFFF")
            ax.add_patch(Polygon(tri, closed=True,
                                 facecolor=tri_color, edgecolor="#333333",
                                 linewidth=0.3))
    # centre text marker (* / .)
    if text:
        ax.text(j, i, text, ha="center", va="center",
                fontsize=8, fontweight="bold", color="#222222")


def render(tbl: pd.DataFrame, out_path: Path, show_trend: bool) -> None:
    set_nature_style()
    tbl = tbl.sort_values("program_order").reset_index(drop=True)

    col_labels = ["128188", "79768", "Consensus", "115574"]
    col_keys   = ["GSE128188", "GSE79768", "Consensus", "GSE115574"]
    n_row, n_col = len(tbl), len(col_keys)

    # Figure size: keep row pitch same as 2B (0.30 in per row) for stackability.
    fig_w = COL_ONEHALF
    fig_h = 0.30 * n_row + 1.35
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.40, right=0.74, top=0.84, bottom=0.20)

    # Identity strip on the left (same position as 2B for visual rhyme).
    STRIP_W = 0.18
    STRIP_X = -0.5 - STRIP_W - 0.1

    # Cells
    for i, row in tbl.iterrows():
        # identity strip
        col = LA_COLOR if row["baseline_direction"] == 1 else RA_COLOR
        ax.add_patch(Rectangle((STRIP_X, i - 0.5), STRIP_W, 1,
                               facecolor=col, edgecolor="white", linewidth=0.5))

        for j, ck in enumerate(col_keys):
            strict = pick_class(row, ck, "rewiring")
            trend  = None
            marker = ""
            if ck == "Consensus":
                # no marker for consensus (derived from lookup, not a test)
                pass
            else:
                p = row.get(f"p_perm_{ck}")
                marker = perm_marker(p)
                if show_trend and ck in DISCOVERY:
                    trend = pick_class(row, ck, "trend")
            draw_cell(ax, j, i, strict_cls=strict, trend_cls=trend, text=marker)

    ax.set_xlim(STRIP_X - 0.1, n_col - 0.5 + 0.1)
    ax.set_ylim(n_row - 0.5, -0.5)
    ax.set_aspect("auto")

    # x ticks / labels (slight rotation so "Consensus" doesn't touch "115574")
    ax.set_xticks(range(n_col))
    ax.set_xticklabels(col_labels, fontsize=6.5, rotation=30,
                       ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0, pad=2)
    # y ticks / labels
    ax.set_yticks(range(n_row))
    ax.set_yticklabels(tbl["program_label"].tolist())
    ax.tick_params(axis="y", length=0, pad=2)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    for lbl, dirn in zip(ax.get_yticklabels(), tbl["baseline_direction"]):
        lbl.set_color(LA_COLOR if dirn == 1 else RA_COLOR)

    # "dir" header for the identity strip (no labels above the cells;
    # the x-tick labels below carry the column identity).
    ax.text(STRIP_X + STRIP_W / 2, -0.80, "dir",
            fontsize=6, ha="center", va="bottom", color="#333333",
            fontweight="bold", clip_on=False)

    fig.suptitle("AF rewiring classification  -  strict fill, trend overlay",
                 fontsize=7.5, y=0.96)

    # --- legend on the right: class swatches ---
    leg_x = 0.76
    leg_y_top = 0.84
    leg_dy = 0.035
    fig.text(leg_x, leg_y_top, "class", fontsize=6, fontweight="bold",
             ha="left", va="bottom", color="#333333")
    for k, cls in enumerate(LEGEND_ORDER):
        y = leg_y_top - (k + 1) * leg_dy
        fig.patches.append(Rectangle(
            (leg_x, y), 0.028, 0.025,
            transform=fig.transFigure,
            facecolor=REWIRING_CLASS_COLORS[cls],
            edgecolor="#888888" if cls == "NA" else "white",
            linewidth=0.5,
        ))
        if cls == "Reversed-weak":
            fig.patches.append(Rectangle(
                (leg_x, y), 0.028, 0.025,
                transform=fig.transFigure,
                facecolor="none", edgecolor="#542788",
                linewidth=0.0, hatch="////",
            ))
        fig.text(leg_x + 0.033, y + 0.012, cls,
                 fontsize=5.5, ha="left", va="center", color="#333333")

    # significance legend
    sig_y = leg_y_top - (len(LEGEND_ORDER) + 2) * leg_dy
    fig.text(leg_x, sig_y + 0.025, "perm p", fontsize=6, fontweight="bold",
             ha="left", va="bottom", color="#333333")
    fig.text(leg_x + 0.000, sig_y, "*  < 0.05",
             fontsize=5.5, ha="left", va="center", color="#333333")
    fig.text(leg_x + 0.000, sig_y - leg_dy, ".  < 0.10",
             fontsize=5.5, ha="left", va="center", color="#333333")

    # identity legend
    id_y = sig_y - 3 * leg_dy
    fig.text(leg_x, id_y + 0.030, "healthy bias", fontsize=6, fontweight="bold",
             ha="left", va="bottom", color="#333333")
    fig.patches.append(Rectangle((leg_x, id_y + 0.005), 0.028, 0.017,
                                 transform=fig.transFigure, facecolor=LA_COLOR,
                                 edgecolor="white", linewidth=0.5))
    fig.text(leg_x + 0.033, id_y + 0.013, "LA-high",
             fontsize=5.5, ha="left", va="center", color=LA_COLOR)
    fig.patches.append(Rectangle((leg_x, id_y - 0.020), 0.028, 0.017,
                                 transform=fig.transFigure, facecolor=RA_COLOR,
                                 edgecolor="white", linewidth=0.5))
    fig.text(leg_x + 0.033, id_y - 0.012, "RA-high",
             fontsize=5.5, ha="left", va="center", color=RA_COLOR)

    # trend-overlay legend
    if show_trend:
        tr_y = id_y - 5 * leg_dy
        fig.text(leg_x, tr_y + 0.030, "trend overlay",
                 fontsize=6, fontweight="bold",
                 ha="left", va="bottom", color="#333333")
        fig.text(leg_x, tr_y - 0.002,
                 "corner triangle\nshows trend-class\nwhen different from\nstrict (discovery).",
                 fontsize=5.5, ha="left", va="top", color="#555555")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {out_path.relative_to(ROOT)}")


def main() -> None:
    tbl = load_data()
    render(tbl, OUT, show_trend=True)


if __name__ == "__main__":
    main()
