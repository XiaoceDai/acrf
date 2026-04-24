#!/usr/bin/env python3
"""
Figure 3B - program x cell_type Δsinglescore heatmap.

Rows: the 10-program panel, Fig 1C order.
Cols: 4 main compartments (atrialCM / fibroblast / endothelial_or_vascular / immune).
Cells: raw Δsinglescore (AF - Ctrl).  Diverging BrBG palette:
    + (brown)  = AF higher
    0 (white)  = no change
    − (teal)   = AF lower
`*` / `.` markers from R4's sig_marker column (p_perm < 0.05 / 0.10).

Hard rule: this script must NOT import LA_COLOR, RA_COLOR, REWIRING_CMAP,
or REWIRING_CLASS_COLORS.  The BrBG axis is the only fill axis.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

from plot_style import AFCTRL_CMAP, COL_ONEHALF, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "figures" / "fig3b_delta_heatmap.pdf"

COMPARTMENTS = ["atrialCM", "fibroblast", "endothelial_or_vascular", "immune"]
COMPARTMENT_LABELS = {
    "atrialCM": "atrialCM",
    "fibroblast": "fibroblast",
    "endothelial_or_vascular": "endo / vasc",
    "immune": "immune",
}


def main() -> None:
    set_nature_style()
    tbl = pd.read_csv(ROOT / "results" / "af_la_program_celltype_delta.tsv", sep="\t")
    tbl = tbl.sort_values(["program_order", "cell_type"]).reset_index(drop=True)

    programs = (tbl[["program_id", "program_label", "program_order",
                     "baseline_direction"]]
                .drop_duplicates()
                .sort_values("program_order")
                .reset_index(drop=True))
    n_row = len(programs)
    n_col = len(COMPARTMENTS)

    # Δscore matrix
    mat = np.full((n_row, n_col), np.nan)
    sig = np.full((n_row, n_col), "", dtype=object)
    for i, p in programs.iterrows():
        for j, ct in enumerate(COMPARTMENTS):
            r = tbl[(tbl["program_id"] == p["program_id"]) &
                    (tbl["cell_type"] == ct)]
            if not r.empty:
                mat[i, j] = r["dScore"].iloc[0]
                mk = r["sig_marker"].iloc[0]
                sig[i, j] = "" if pd.isna(mk) else mk

    vmax = max(0.02, float(np.nanmax(np.abs(mat))) * 1.1)
    vmin = -vmax
    norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list("afctrl_brbg", AFCTRL_CMAP, N=256)

    # geometry
    fig_w = COL_ONEHALF
    fig_h = 0.30 * n_row + 1.60
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.40, right=0.86, top=0.84, bottom=0.22)

    # Fig 3B axis is Δ singscore (AF - Ctrl) only; no LA/RA identity
    # overlay.  Cross-figure row alignment is carried by `program_order`
    # (shared with Fig 1C / Fig 2), which is enforced by sorting above.
    for i in range(n_row):
        for j in range(n_col):
            v = mat[i, j]
            rgba = (1, 1, 1, 0) if np.isnan(v) else cmap(norm(v))
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0,
                                   facecolor=rgba, edgecolor="white",
                                   linewidth=0.5))
            if not np.isnan(v):
                # cell label: numeric value
                txt_color = "white" if abs(norm(v) - 0.5) > 0.30 else "#222222"
                ax.text(j, i - 0.10, f"{v:+.3f}",
                        ha="center", va="center", fontsize=5.8,
                        color=txt_color, fontweight="bold")
                if sig[i, j]:
                    ax.text(j, i + 0.25, sig[i, j],
                            ha="center", va="center", fontsize=9,
                            fontweight="bold", color=txt_color)

    ax.set_xlim(-0.6, n_col - 0.5 + 0.1)
    ax.set_ylim(n_row - 0.5, -0.5)
    ax.set_aspect("auto")

    ax.set_xticks(range(n_col))
    ax.set_xticklabels([COMPARTMENT_LABELS[c] for c in COMPARTMENTS],
                       fontsize=6.5, rotation=20,
                       ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0, pad=2)

    ax.set_yticks(range(n_row))
    ax.set_yticklabels(programs["program_label"].tolist())
    ax.tick_params(axis="y", length=0, pad=2)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)

    # colour bar
    cbar_ax = fig.add_axes([0.42, 0.08, 0.38, 0.023])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cbar_ax, orientation="horizontal")
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(width=0.4, length=2, labelsize=5.5, pad=1)
    cbar.set_ticks([round(vmin, 3), round(vmin / 2, 3), 0,
                    round(vmax / 2, 3), round(vmax, 3)])
    cbar.ax.set_title(r"$\Delta$ singscore (AF $-$ Ctrl)",
                      fontsize=6, pad=3)
    fig.text(0.42, 0.14, "AF lower", fontsize=5.5, ha="left", va="bottom",
             color="#01665E")
    fig.text(0.80, 0.14, "AF higher", fontsize=5.5, ha="right", va="bottom",
             color="#8C510A")

    # significance key on the right gutter (only axis encoding in Fig 3B
    # besides the BrBG Δ scale)
    gutter_x = n_col - 0.5 + 0.45
    ax.text(gutter_x, -0.5, "perm p",
            fontsize=6, ha="left", va="top", fontweight="bold", clip_on=False)
    ax.text(gutter_x, 0.15,
            r"$\star$ < 0.05" + "\n" + ".  < 0.10",
            fontsize=5.5, ha="left", va="top", clip_on=False)

    # title
    fig.suptitle(
        "AF-associated LA-arm shifts are modest and distributed across compartments",
        fontsize=7.4, y=0.96, fontweight="bold",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
