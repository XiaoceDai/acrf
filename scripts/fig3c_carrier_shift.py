#!/usr/bin/env python3
"""
Figure 3C - Carrier shift table (three orthogonal columns).

Column 1  Fig 1D dominant_cell_type    (signed attribution argmax)
Column 2  Fig 1 top-expressor          (AF-naive mean singscore argmax)
Column 3  AF major carrier             (Fig 3 Δscore argmax under q_perm < 0.10)

Each cell is filled by its compartment colour; if the value is
"unresolved" / "Unresolved" it is rendered grey.  The three columns are
NOT comparable one-for-one (different quantities), so this panel deliberately
uses a side-by-side reading rather than a shift arrow; differences between
columns are meant to be read as "these three metrics disagree, which is the
biology".

Hard rule: no LA/RA colours and no rewiring colours; only the discrete
compartment palette from plot_style.COMPARTMENT_COLORS.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

from plot_style import COL_ONEHALF, COMPARTMENT_COLORS, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "figures" / "fig3c_carrier_shift.pdf"

COL_LABELS = [
    "healthy\ndominant\n(Fig 1D)",
    "healthy\ntop expr.\n(Fig 1)",
    "AF major\ncarrier\n(Fig 3)",
]
COL_KEYS = [
    "healthy_dominant_CT",
    "healthy_top_expressor_CT",
    "AF_major_carrier_CT",
]


def main() -> None:
    set_nature_style()
    tbl = pd.read_csv(ROOT / "results" / "figure3_carrier_shift.tsv", sep="\t")
    tbl = tbl.sort_values("program_order").reset_index(drop=True)

    n_row = len(tbl)
    n_col = len(COL_KEYS)

    fig_w = COL_ONEHALF
    fig_h = 0.30 * n_row + 1.65
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.40, right=0.74, top=0.78, bottom=0.14)

    # cell-type to short label
    SHORT = {
        "atrialCM":                "atrialCM",
        "fibroblast":              "fibroblast",
        "endothelial_or_vascular": "endo/vasc",
        "immune":                  "immune",
        "unresolved":              "unresolved",
        "Unresolved":              "",        # empty tile in grey (explained in legend)
    }
    DARK = {"atrialCM", "fibroblast", "endothelial_or_vascular", "immune"}

    for i, row in tbl.iterrows():
        for j, key in enumerate(COL_KEYS):
            val = row[key]
            if pd.isna(val):
                fill, label_text, text_color = "#FFFFFF", "", "#222222"
            else:
                fill = COMPARTMENT_COLORS.get(val, "#BDBDBD")
                label_text = SHORT.get(val, str(val))
                text_color = "white" if val in DARK else "#333333"

            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0,
                                   facecolor=fill, edgecolor="white",
                                   linewidth=0.5))
            if label_text:
                ax.text(j, i, label_text,
                        ha="center", va="center", fontsize=6,
                        color=text_color, fontweight="bold")

    ax.set_xlim(-0.6, n_col - 0.5 + 0.1)
    ax.set_ylim(n_row - 0.5, -0.5)
    ax.set_aspect("auto")

    # Move column headers to the top so they don't collide with each other.
    ax.set_xticks([])
    for j, lbl in enumerate(COL_LABELS):
        ax.text(j, -0.85, lbl,
                fontsize=6, ha="center", va="bottom",
                linespacing=1.15, clip_on=False, color="#333333",
                fontweight="bold")
    ax.set_yticks(range(n_row))
    ax.set_yticklabels(tbl["program_label"].tolist())
    ax.tick_params(axis="y", length=0, pad=2)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)

    # legend on the right: compartment swatches
    leg_x = 0.77
    leg_y_top = 0.80
    leg_dy = 0.055
    fig.text(leg_x, leg_y_top, "cell type",
             fontsize=6, fontweight="bold", ha="left", va="bottom",
             color="#333333")
    legend_entries = [
        ("atrialCM",                 COMPARTMENT_COLORS["atrialCM"]),
        ("fibroblast",               COMPARTMENT_COLORS["fibroblast"]),
        ("endo / vasc",              COMPARTMENT_COLORS["endothelial_or_vascular"]),
        ("immune",                   COMPARTMENT_COLORS["immune"]),
        ("unresolved",               COMPARTMENT_COLORS["unresolved"]),
    ]
    for k, (name, col) in enumerate(legend_entries):
        y = leg_y_top - (k + 1) * leg_dy
        fig.patches.append(Rectangle((leg_x, y), 0.028, 0.028,
                                     transform=fig.transFigure,
                                     facecolor=col, edgecolor="white",
                                     linewidth=0.4))
        fig.text(leg_x + 0.034, y + 0.014, name,
                 fontsize=5.5, ha="left", va="center", color="#222222")

    # read-the-table hint
    hint_y = leg_y_top - (len(legend_entries) + 2) * leg_dy
    fig.text(leg_x, hint_y + 0.030, "read as",
             fontsize=6, fontweight="bold", ha="left", va="bottom",
             color="#333333")
    fig.text(leg_x, hint_y,
             "col 1  healthy concordance\n          argmax (Fig 1D)\n\n"
             "col 2  healthy absolute\n          expression argmax\n          (Fig 1 pseudobulk)\n\n"
             "col 3  AF Δ argmax\n          under q<0.10\n"
             "          grey = unresolved,\n          i.e. no cell type\n          passes the gate",
             fontsize=5.2, ha="left", va="top", color="#555555",
             linespacing=1.35)

    # title
    fig.suptitle("Healthy carriers  vs  AF major carrier  per program",
                 fontsize=8, y=0.97)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
