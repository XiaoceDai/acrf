#!/usr/bin/env python3
"""
Figure 2E - GSE115574 replication lollipop.

Rows: consensus programs (rows where consensus_class is not NA), in Fig 1C order.
Row annotations (left -> right):
    1. identity strip (LA=blue, RA=red)
    2. consensus class strip (orange/purple/grey discrete colours)
    3. program short label

x-axis: dCAI_rel in GSE115574 (AF - SR) on the healthy-fraction scale.
Glyphs per row:
    *  solid dot  = 115574 observed dCAI_rel, coloured by sign-concordance
         - orange  : sign matches consensus expectation
         - purple  : sign opposite to expectation
         - grey    : 115574 baseline not supported OR Preserved consensus
    *  hollow circle = discovery mean dCAI_rel
         (average of 128188 and 79768 dCAI_rel; whichever is baseline-supported)
    *  lollipop stick from x=0 to the 115574 dot.
    *  text marker to the right of the dot:  * p<0.05, . p<0.10 (on 115574 p_perm)

Vertical reference line at x=0.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from plot_style import (
    COL_ONEHALF, LA_COLOR, RA_COLOR,
    REWIRING_CLASS_COLORS,
    set_nature_style,
)

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "figures" / "fig2e_replication_lollipop.pdf"

# Sign expected for each consensus class:  + means "dCAI should be positive
# (AF stronger than SR along healthy direction)", - means "dCAI should be
# negative (AF weaker or reversed)", 0 means "near zero".
CLASS_EXPECTED_SIGN = {
    "Amplified":       +1,
    "Amplified-weak":  +1,
    "Preserved":        0,
    "Attenuated":      -1,
    "Attenuated-weak": -1,
    "Reversed":        -1,
    "Reversed-weak":   -1,
}

DELTA_REL = 0.30        # same as spec §3.5; used for Preserved sign check


def discovery_mean_dCAI(row) -> float:
    """Mean dCAI_rel across whichever discovery cohorts have baseline support.
    Falls back to unweighted mean of both if neither supported.
    """
    vals = []
    for c in ("GSE128188", "GSE79768"):
        if row[f"sr_baseline_supported_{c}"]:
            v = row[f"dCAI_rel_{c}"]
            if not math.isnan(v):
                vals.append(v)
    if not vals:
        for c in ("GSE128188", "GSE79768"):
            v = row[f"dCAI_rel_{c}"]
            if not math.isnan(v):
                vals.append(v)
    if not vals:
        return float("nan")
    return float(np.mean(vals))


def replication_colour(row, d115: float) -> str:
    """Colour for the 115574 lollipop dot."""
    cc = row["consensus_class"]
    if not row["sr_baseline_supported_GSE115574"]:
        return REWIRING_CLASS_COLORS["Low-support"]
    exp = CLASS_EXPECTED_SIGN.get(cc)
    if math.isnan(d115) or exp is None:
        return REWIRING_CLASS_COLORS["Low-support"]
    if exp == 0:
        # Preserved expectation: colour orange if within +/- delta (stays flat),
        # colour purple if deviates substantially.
        return (REWIRING_CLASS_COLORS["Preserved"] if abs(d115) <= DELTA_REL
                else REWIRING_CLASS_COLORS["Reversed-weak"])
    return (REWIRING_CLASS_COLORS["Amplified"] if np.sign(d115) == exp
            else REWIRING_CLASS_COLORS["Reversed"])


def perm_marker(p: float) -> str:
    if math.isnan(p):
        return ""
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "."
    return ""


def main() -> None:
    set_nature_style()

    progs = pd.read_csv(ROOT / "results" / "af_chamber_rewiring_programs.tsv",
                        sep="\t")
    tbl = progs[progs["consensus_class"].notna()].copy()
    tbl = tbl.sort_values("program_order").reset_index(drop=True)

    tbl["discovery_mean_dCAI"] = tbl.apply(discovery_mean_dCAI, axis=1)
    tbl["rep_color"] = tbl.apply(
        lambda r: replication_colour(r, r["dCAI_rel_GSE115574"]),
        axis=1,
    )

    n = len(tbl)
    if n == 0:
        print("no consensus programs; nothing to plot")
        return

    # Figure geometry: a wider figure than 2B/2C so the row labels, strips,
    # main lollipop axis, and legend each have their own breathing room.
    fig_w = 6.0    # inches
    fig_h = 0.55 * n + 1.60
    fig = plt.figure(figsize=(fig_w, fig_h))

    # ----- three horizontal zones -----
    # [0.04, 0.28)  row labels + strips (a blank axes we only use for text)
    # [0.30, 0.74]  main lollipop axes
    # [0.76, 0.99)  legend
    strip_ax = fig.add_axes([0.02, 0.18, 0.27, 0.68])
    ax       = fig.add_axes([0.30, 0.18, 0.44, 0.68])
    leg_ax   = fig.add_axes([0.75, 0.18, 0.24, 0.68])
    strip_ax.axis("off")
    leg_ax.axis("off")

    # x range symmetric around 0, covering both discovery and 115574.
    all_x = np.concatenate([
        tbl["discovery_mean_dCAI"].to_numpy(),
        tbl["dCAI_rel_GSE115574"].to_numpy(),
    ])
    max_abs = float(np.nanmax(np.abs(all_x))) if len(all_x) else 0.5
    xlim = max(0.7, max_abs * 1.15)
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(n - 0.5, -0.5)

    # Reference vertical line at 0 and pale band for |d| <= delta.
    ax.axvline(0.0, color="#888888", lw=0.5, ls="--", zorder=1)
    ax.axvspan(-DELTA_REL, DELTA_REL, color="#F7F7F7", alpha=0.6, zorder=0, lw=0)

    # ----- row annotation strips on strip_ax (axes-fraction coords) -----
    strip_ax.set_xlim(0, 1)
    strip_ax.set_ylim(n - 0.5, -0.5)
    CC_X, ID_X = 0.78, 0.90           # class strip, then dir strip
    STRIP_W = 0.08
    # Row labels, then class strip, then dir strip.
    for i, row in tbl.iterrows():
        id_col = LA_COLOR if row["baseline_direction"] == 1 else RA_COLOR
        cc_col = REWIRING_CLASS_COLORS.get(row["consensus_class"], "#FFFFFF")
        strip_ax.add_patch(Rectangle((CC_X, i - 0.45), STRIP_W, 0.9,
                                     facecolor=cc_col, edgecolor="#DDDDDD",
                                     linewidth=0.3))
        strip_ax.add_patch(Rectangle((ID_X, i - 0.45), STRIP_W, 0.9,
                                     facecolor=id_col, edgecolor="white",
                                     linewidth=0.3))
        # Row label aligned to the right, just left of the class strip.
        strip_ax.text(CC_X - 0.02, i, row["program_label"],
                      fontsize=6.5, ha="right", va="center",
                      color=id_col)

    # Column headers on the strips.
    strip_ax.text(CC_X + STRIP_W / 2, -0.85, "class",
                  fontsize=6, ha="center", va="bottom", color="#333333",
                  fontweight="bold")
    strip_ax.text(ID_X + STRIP_W / 2, -0.85, "dir",
                  fontsize=6, ha="center", va="bottom", color="#333333",
                  fontweight="bold")

    # --- lollipops ---
    for i, row in tbl.iterrows():
        d115 = row["dCAI_rel_GSE115574"]
        ddisc = row["discovery_mean_dCAI"]
        fill = row["rep_color"]
        # stick from 0 to 115574 dot
        if not math.isnan(d115):
            ax.plot([0, d115], [i, i], color=fill, lw=1.1, zorder=2,
                    solid_capstyle="round")
            ax.scatter(d115, i, s=44, marker="o",
                       facecolor=fill, edgecolor="#333333", linewidth=0.5,
                       zorder=3)
            mark = perm_marker(row["p_perm_GSE115574"])
            if mark:
                ax.text(d115 + 0.03 * xlim * np.sign(d115) + 0.03 * xlim,
                        i, mark, fontsize=9, fontweight="bold", color="#222222",
                        ha="left" if d115 >= 0 else "right", va="center",
                        clip_on=False)
        # discovery reference: hollow circle
        if not math.isnan(ddisc):
            ax.scatter(ddisc, i, s=36, marker="o",
                       facecolor="white", edgecolor="#444444", linewidth=0.7,
                       zorder=4)

    # main axes cosmetics
    ax.set_xlabel(r"$\Delta\mathrm{aCAI}_{rel}\ =\ \mathrm{mean}_{AF} - \mathrm{mean}_{SR}$",
                  fontsize=6.5)
    ax.set_yticks([])                   # labels live on strip_ax
    ax.tick_params(axis="x", labelsize=6)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.6)

    # Figure title and subtitle at the very top.
    fig.suptitle("Replication in GSE115574 (AF + mitral regurgitation)",
                 fontsize=8, y=0.96)
    fig.text(0.30, 0.91,
             "solid dot: 115574 observed ; hollow dot: discovery mean",
             fontsize=5.5, ha="left", va="bottom", color="#555555")

    # ----- legend on leg_ax -----
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", linestyle="None", markersize=5,
               markerfacecolor=REWIRING_CLASS_COLORS["Amplified"],
               markeredgecolor="#333333", markeredgewidth=0.5,
               label="sign matches consensus"),
        Line2D([], [], marker="o", linestyle="None", markersize=5,
               markerfacecolor=REWIRING_CLASS_COLORS["Reversed"],
               markeredgecolor="#333333", markeredgewidth=0.5,
               label="sign opposite to consensus"),
        Line2D([], [], marker="o", linestyle="None", markersize=5,
               markerfacecolor=REWIRING_CLASS_COLORS["Preserved"],
               markeredgecolor="#333333", markeredgewidth=0.5,
               label=r"Preserved + $|d| \leq 0.3$"),
        Line2D([], [], marker="o", linestyle="None", markersize=5,
               markerfacecolor=REWIRING_CLASS_COLORS["Low-support"],
               markeredgecolor="#333333", markeredgewidth=0.5,
               label="115574 low-support"),
        Line2D([], [], marker="o", linestyle="None", markersize=5,
               markerfacecolor="white", markeredgecolor="#444444", markeredgewidth=0.7,
               label="discovery mean"),
    ]
    leg_ax.legend(handles=handles, loc="center left",
                  fontsize=5.5, frameon=False,
                  handletextpad=0.4, labelspacing=0.6,
                  handlelength=1.0, borderaxespad=0)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}  ({n} consensus programs)")


if __name__ == "__main__":
    main()
