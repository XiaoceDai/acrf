#!/usr/bin/env python3
"""
Figure 3D  Program × cell-type Δsinglescore heatmap (full 30 × 11).

Provides the raw transparency layer underlying Fig 3B:  every NCP primary
program's Δ singscore in every LA-arm cell type, blocked by axis
family.  Significance markers as in Fig 3B (* p < 0.05, . p < 0.10).

Inputs:
    results/fig3_primary_program_celltype_delta.tsv
    results/fig2_primary_axis_families.tsv
    results/fig2_primary_acai_long.tsv

Output:
    figures/fig3d_program_celltype_heatmap.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label

DELTA_TSV = ROOT / "results" / "fig3_primary_program_celltype_delta.tsv"
FAM_TSV   = ROOT / "results" / "fig2_primary_axis_families.tsv"
ACAI_TSV  = ROOT / "results" / "fig2_primary_acai_long.tsv"
FIG_OUT   = ROOT / "figures" / "fig3d_program_celltype_heatmap.pdf"

LA_COLOR  = "#D6604D"
RA_COLOR  = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID  = "#9C9C9C"
BRBG_BROWN = "#8C510A"
BRBG_TEAL  = "#01665E"

CT_LABELS = {
    "atrialCM":                "atrialCM",
    "fibroblast":              "fibroblast",
    "endothelial":             "endothelial",
    "lymphatic_endothelial":   "lymph endo",
    "pericyte":                "pericyte",
    "vascular_smooth_muscle":  "VSMC",
    "adipocyte":               "adipocyte",
    "macrophage":              "macrophage",
    "lymphocyte":              "lymphocyte",
    "mast_cell":               "mast cell",
    "neuronal":                "neuronal",
}
CT_ORDER = list(CT_LABELS.keys())

FAMILY_PALETTE = [
    "#4F86C6", "#6FA8DC", "#9EC4E0",
    "#82B5A1", "#A2C8B0", "#B7D2BC",
    "#D4B58E", "#C8C26C",
    "#E89C8B", "#D86A6A", "#B94B4B",
]


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "axes.titlesize":    9,
        "axes.labelsize":    8,
        "xtick.labelsize":   7,
        "ytick.labelsize":   6.5,
        "axes.linewidth":    0.6,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def main() -> None:
    set_publication_style()

    df = pd.read_csv(DELTA_TSV, sep="\t")

    # Order programs by family, then by program_order within family.
    ordered = (df.drop_duplicates("program_id")
                  [["program_id", "program_label", "program_order",
                    "family_id", "baseline_direction"]]
                  .sort_values(["family_id", "program_order"])
                  .reset_index(drop=True))
    program_order = ordered["program_id"].tolist()
    K = len(program_order)
    C = len(CT_ORDER)

    # Pivot Δ, p_perm
    M = df.pivot_table(index="program_id", columns="cell_type",
                          values="dScore").reindex(
                              index=program_order, columns=CT_ORDER)
    Mp = df.pivot_table(index="program_id", columns="cell_type",
                            values="p_perm").reindex(
                                index=program_order, columns=CT_ORDER)
    Mincl = df.pivot_table(index="program_id", columns="cell_type",
                                values="cell_type_included").reindex(
                                    index=program_order, columns=CT_ORDER)

    # Pretty short labels for programs
    prog_short = {pid: short_pathway_label(pid, 26) for pid in program_order}

    # ---- figure layout ------------------------------------------------
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(180 * MM, 180 * MM))
    grid = fig.add_gridspec(
        nrows=2, ncols=4,
        height_ratios=[2.0, 30],
        width_ratios=[1.0, 18, 0.4, 0.6],
        hspace=0.04, wspace=0.05,
        left=0.27, right=0.94, top=0.95, bottom=0.07)

    ax_top  = fig.add_subplot(grid[0, 1])
    ax_main = fig.add_subplot(grid[1, 1])
    ax_fam  = fig.add_subplot(grid[1, 0])
    ax_cbar = fig.add_subplot(grid[1, 3])

    # colour map (BrBG diverging)
    finite_vals = M.values[np.isfinite(M.values)]
    vmax = float(np.nanpercentile(np.abs(finite_vals), 95))
    vmax = max(0.005, round(vmax, 3))
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "BrBG_div", [(0.0, BRBG_TEAL), (0.5, "#FFFFFF"),
                       (1.0, BRBG_BROWN)])
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    # ---- heatmap cells ----
    for i, pid in enumerate(program_order):
        for j, ct in enumerate(CT_ORDER):
            v = M.iloc[i, j]
            p = Mp.iloc[i, j]
            included = bool(Mincl.iloc[i, j]) if not pd.isna(Mincl.iloc[i, j]) else False
            if pd.isna(v) or not included:
                facec = "#F2F2F2"
            else:
                facec = cmap(norm(v))
            ax_main.add_patch(mpatches.Rectangle(
                (j, i), 1.0, 1.0, facecolor=facec,
                edgecolor="white", linewidth=0.4))
            if pd.isna(v) or not included:
                continue
            if not pd.isna(p) and p < 0.10:
                if p < 0.05:
                    badge_face, mark = "#000000", "*"
                    text_color = "white"
                    fs = 7.0
                else:
                    badge_face, mark = "#FFFFFF", "."
                    text_color = "black"
                    fs = 8.0
                ax_main.add_patch(mpatches.Circle(
                    (j + 0.85, i + 0.18), 0.13,
                    facecolor=badge_face,
                    edgecolor="black", linewidth=0.4, zorder=4))
                ax_main.text(j + 0.85, i + 0.20, mark,
                                ha="center", va="center",
                                fontsize=fs, fontweight="bold",
                                color=text_color, zorder=5)

    ax_main.set_xlim(0, C); ax_main.set_ylim(K, 0)
    ax_main.set_xticks([]); ax_main.set_yticks([])
    for s in ax_main.spines.values():
        s.set_visible(True); s.set_linewidth(0.4); s.set_color("#999")

    # ---- top labels (cell types) — anchored 45° rotation ----
    for j, ct in enumerate(CT_ORDER):
        ax_top.plot([j + 0.5, j + 0.5], [0.0, 0.10],
                      color=GREY_MID, lw=0.7, zorder=2)
        ax_top.text(j + 0.5, 0.16, CT_LABELS[ct],
                      ha="left", va="bottom",
                      rotation=45, rotation_mode="anchor",
                      fontsize=7.0, fontweight="bold",
                      color="black")
    ax_top.set_xlim(0, C); ax_top.set_ylim(0, 1); ax_top.axis("off")

    # ---- left family colour bar ----
    fam_order = sorted(set(ordered["family_id"]))
    fam_to_color = {fid: FAMILY_PALETTE[k % len(FAMILY_PALETTE)]
                       for k, fid in enumerate(fam_order)}
    fam_blocks: dict[int, list[int]] = {}
    for i, pid in enumerate(program_order):
        fid = int(ordered.loc[ordered["program_id"] == pid,
                                  "family_id"].iloc[0])
        fam_blocks.setdefault(fid, []).append(i)
        ax_fam.add_patch(mpatches.Rectangle(
            (0, i), 1.0, 1.0, facecolor=fam_to_color[fid],
            edgecolor="none"))
    for fid, idxs in fam_blocks.items():
        yc = (idxs[0] + idxs[-1] + 1) / 2
        ax_fam.text(0.5, yc, f"F{fid}",
                      ha="center", va="center",
                      fontsize=6.5, fontweight="bold",
                      color="white" if fid in (1, 9, 10, 11, 4) else "black")
    ax_fam.set_xlim(0, 1); ax_fam.set_ylim(K, 0); ax_fam.axis("off")

    # ---- left program labels (figure-level) ----
    for i, pid in enumerate(program_order):
        fid = int(ordered.loc[ordered["program_id"] == pid,
                                  "family_id"].iloc[0])
        d   = int(ordered.loc[ordered["program_id"] == pid,
                                  "baseline_direction"].iloc[0])
        col = LA_COLOR if d > 0 else RA_COLOR
        badge = "LA" if d > 0 else "RA"
        ax_fam.text(-0.35, i + 0.5,
                      f"{badge}  {prog_short[pid]}",
                      ha="right", va="center",
                      fontsize=6.2, color=col,
                      transform=ax_fam.transData,
                      clip_on=False)

    # ---- colorbar ----
    cb = mpl.colorbar.ColorbarBase(
        ax_cbar, cmap=cmap, norm=norm, orientation="vertical")
    cb.set_label("Δ singscore  (AF − Ctrl)", fontsize=6.8)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)

    # ---- title ----
    fig.suptitle(
        f"Fig 3D  Program × cell-type Δ singscore (LA snRNA arm, GSE255612) "
        f"— {K} NCP_v2 programs × {C} cell types",
        fontsize=8.5, fontweight="bold", y=0.99, x=0.02, ha="left")

    # ---- caption ----
    fig.text(0.02, 0.022,
              "Cells: Δ singscore (AF − Ctrl) with 10 000 MC perm test.   "
              "Markers: black-filled  *  p_perm < 0.05;  white  .  p_perm < 0.10.   "
              "Grey cells: < 4 donors per arm (excluded).\n"
              "Rows: 30 NCP_v2 programs grouped by axis family (F1–F11; "
              "left colour bar). LA-biased programs (blue) on top, "
              "RA-biased (red) at bottom of each family block.",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
