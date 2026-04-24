#!/usr/bin/env python3
"""
Figure 2B - program x cohort_rhythm heatmap of mean(aCAI_rel).

Rows: 10 Fig 1 chamber programs, ordered by signed bulk NES (Fig 1C order).
Cols: six (cohort, rhythm) groups:
        128188_SR, 128188_AF | 79768_SR, 79768_AF | 115574_SR, 115574_AF

Cells: mean aCAI_rel across donors. Diverging PuOr palette, centred at 0.
       +1 = reproduces Fig 1 healthy bias exactly.
       0  = chamber bias absent.
       -1 = reversed at one healthy magnitude.

This panel encodes magnitude only; classification & significance live in 2C.
Output: vector PDF.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

from plot_style import (
    COL_ONEHALF, LA_COLOR, RA_COLOR, REWIRING_CMAP,
    set_nature_style,
)

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "figures" / "fig2b_rewiring_heatmap.pdf"

COHORTS = ["GSE128188", "GSE79768", "GSE115574"]
COHORT_SHORT = {"GSE128188": "128188", "GSE79768": "79768", "GSE115574": "115574"}


def load_table() -> pd.DataFrame:
    """Build a long (program x cohort_rhythm x value) frame from P4 outputs."""
    src = pd.read_csv(ROOT / "results" / "figure2_source_data.tsv", sep="\t")
    agg = (src.groupby(["cohort", "rhythm", "program_id", "program_label",
                         "program_order", "baseline_direction"], as_index=False)
              ["aCAI_rel"].mean()
              .rename(columns={"aCAI_rel": "mean_aCAI_rel"}))
    return agg


def cell_grid(agg: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame, list[tuple[str,str]]]:
    programs = (agg[["program_id", "program_label", "program_order", "baseline_direction"]]
                .drop_duplicates().sort_values("program_order").reset_index(drop=True))
    col_keys = []
    for c in COHORTS:
        col_keys.append((c, "SR"))
        col_keys.append((c, "AF"))

    mat = np.full((len(programs), len(col_keys)), np.nan)
    for i, p in programs.iterrows():
        for j, (c, r) in enumerate(col_keys):
            sub = agg[(agg["program_id"] == p["program_id"]) &
                      (agg["cohort"] == c) & (agg["rhythm"] == r)]
            if not sub.empty:
                mat[i, j] = sub["mean_aCAI_rel"].iloc[0]
    return mat, programs, col_keys


def main() -> None:
    set_nature_style()
    agg = load_table()
    mat, programs, col_keys = cell_grid(agg)
    n_row, n_col = mat.shape

    # Diverging palette: purple (reversed) - white - orange (preserved/amplified).
    cmap = LinearSegmentedColormap.from_list("rewiring_puor", REWIRING_CMAP, N=256)
    # Data-driven symmetric upper bound; pin vmin to -1 (= full reversal) and
    # vmax to max(data, 1.0) so the "healthy baseline = +1" lives at the
    # orange-mid shade rather than the palette extreme.
    vmax = max(1.0, float(np.nanmax(mat)))
    vmin = -1.0
    norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)

    # Figure geometry.  Keep the row pitch identical to Fig 1C so rows align
    # when the two panels are tiled vertically in a composition.
    fig_w = COL_ONEHALF
    fig_h = 0.30 * n_row + 1.35
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.40, right=0.86, top=0.84, bottom=0.20)

    # Identity strip lives just left of the heatmap; cell width matched to
    # heatmap cell width.
    STRIP_W = 0.18         # in data-coordinate units (cell width = 1.0)
    STRIP_X = -0.5 - STRIP_W - 0.1

    # --- cells (vector rectangles, not imshow) ---
    for i in range(n_row):
        for j in range(n_col):
            v = mat[i, j]
            rgba = (1, 1, 1, 0) if np.isnan(v) else cmap(norm(v))
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0,
                                   facecolor=rgba, edgecolor="white", linewidth=0.5))
            if not np.isnan(v):
                txt_color = "white" if abs(norm(v) - 0.5) > 0.30 else "#222222"
                ax.text(j, i, f"{v:+.2f}",
                        ha="center", va="center", fontsize=5.8,
                        color=txt_color, fontweight="bold")

    # --- identity strip (LA=blue, RA=red) ---
    for i, dirn in enumerate(programs["baseline_direction"]):
        col = LA_COLOR if dirn == 1 else RA_COLOR
        ax.add_patch(Rectangle((STRIP_X, i - 0.5), STRIP_W, 1.0,
                               facecolor=col, edgecolor="white", linewidth=0.5))

    # --- column axis ---
    ax.set_xlim(STRIP_X - 0.1, n_col - 0.5 + 0.1)
    ax.set_ylim(n_row - 0.5, -0.5)
    ax.set_aspect("auto")

    ax.set_xticks(range(n_col))
    ax.set_xticklabels([r for (_, r) in col_keys], fontsize=6)
    # Secondary axis: cohort brackets above.
    for idx, c in enumerate(COHORTS):
        x0 = 2 * idx - 0.5
        x1 = 2 * idx + 1.5
        y  = -0.95
        ax.plot([x0, x1], [y, y], color="#444444", lw=0.6, clip_on=False)
        ax.text((x0 + x1) / 2, y - 0.25, COHORT_SHORT[c],
                fontsize=7, ha="center", va="top", fontweight="bold", clip_on=False)

    # Separators between cohorts.
    for idx in (1, 2):
        ax.axvline(2 * idx - 0.5, color="white", lw=1.8)

    # --- y ticks / row labels ---
    ax.set_yticks(range(n_row))
    ax.set_yticklabels(programs["program_label"].tolist())
    ax.tick_params(axis="both", length=0, pad=2)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    # Row label colours follow the identity strip.
    for lbl, dirn in zip(ax.get_yticklabels(), programs["baseline_direction"]):
        lbl.set_color(LA_COLOR if dirn == 1 else RA_COLOR)

    # Title and identity-strip caption.
    ax.text(STRIP_X + STRIP_W / 2, -0.95, "dir",
            fontsize=6, ha="center", va="top", color="#333333", clip_on=False)

    # --- colour bar ---
    cbar_ax = fig.add_axes([0.42, 0.07, 0.32, 0.022])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cbar_ax, orientation="horizontal")
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(width=0.4, length=2, labelsize=5.5, pad=1)
    # Keep ticks sparse and always include the data max as the last tick so
    # it does not collide with "1.0"; drop the 1.0 tick if max is close.
    if vmax > 1.12:
        cbar.set_ticks([-1.0, -0.5, 0.0, 0.5, round(vmax, 1)])
    else:
        cbar.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    cbar.ax.set_title(r"mean $\mathrm{aCAI}_{rel}$", fontsize=6, pad=3)
    # Flank labels (both far enough apart from the ticks).
    fig.text(0.42, 0.13, "reversed",
             fontsize=5.5, ha="left", va="bottom", color="#542788")
    fig.text(0.74, 0.13, "healthy direction",
             fontsize=5.5, ha="right", va="bottom", color="#B35806")

    # Identity legend (top-right).
    handles = [Rectangle((0, 0), 1, 1, facecolor=LA_COLOR),
               Rectangle((0, 0), 1, 1, facecolor=RA_COLOR)]
    ax.legend(handles, ["LA-high", "RA-high"],
              loc="lower left", bbox_to_anchor=(1.02, 0.0),
              frameon=False, fontsize=6, handlelength=0.8, handleheight=0.6,
              title="Healthy\nbias", title_fontsize=6)

    fig.suptitle(
        "AF chamber rewiring landscape across paired atrial cohorts",
        fontsize=7.4, y=0.965, fontweight="bold",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
