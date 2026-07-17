#!/usr/bin/env python3
"""Figure 3B  Family × cell-type Δsinglescore heatmap."""
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

LA_COLOR = "#D6604D"
RA_COLOR = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
BRBG_BROWN = "#8C510A"
BRBG_TEAL  = "#01665E"

DELTA_TSV = ROOT / "results" / "fig3_primary_family_celltype_delta.tsv"
FAM_TSV   = ROOT / "results" / "fig2_primary_axis_families.tsv"
PROG_TSV  = ROOT / "results" / "fig3_primary_program_celltype_delta.tsv"
ACAI_TSV  = ROOT / "results" / "fig2_primary_acai_long.tsv"
FIG_OUT   = ROOT / "figures" / "fig3b_family_celltype_heatmap.pdf"

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
# logical grouping: structural cores → vasculature → mural → metabolic →
# immune → other.  This ordering is purely cosmetic; perm tests are
# independent across cells.
CT_ORDER = [
    "atrialCM", "fibroblast",
    "endothelial", "lymphatic_endothelial",
    "pericyte", "vascular_smooth_muscle",
    "adipocyte",
    "macrophage", "lymphocyte", "mast_cell",
    "neuronal",
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
        "ytick.labelsize":   7,
        "axes.linewidth":    0.6,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def main() -> None:
    set_publication_style()

    df = pd.read_csv(DELTA_TSV, sep="\t")
    df = df.dropna(subset=["dScore"])  # keep included cells

    long = pd.read_csv(ACAI_TSV, sep="\t")
    fam = pd.read_csv(FAM_TSV, sep="\t")[["pathway", "family_id"]].rename(
        columns={"pathway": "program_id"})
    long = long.merge(fam, on="program_id", how="left")
    fam_rep = (long.assign(absA=long["aCAI"].abs())
                    .sort_values("absA", ascending=False)
                    .groupby("family_id")
                    .agg(rep=("program_id", "first"))
                    .reset_index())
    sys.path.insert(0, str(ROOT / "scripts"))
    from _pathway_labels import short_pathway_label
    fam_rep["short"] = fam_rep["rep"].apply(
        lambda s: short_pathway_label(s, 22))

    df = df.merge(fam_rep[["family_id", "short"]], on="family_id",
                    how="left")

    fam_meta = (df.drop_duplicates("family_id")
                    [["family_id", "direction_label", "n_programs", "short"]]
                    .copy())
    fam_meta["family_id"] = fam_meta["family_id"].astype(int)
    la = fam_meta[fam_meta["direction_label"] == "LA"] \
              .sort_values("family_id")
    ra = fam_meta[fam_meta["direction_label"] == "RA"] \
              .sort_values("family_id")
    fam_order = la["family_id"].tolist() + ra["family_id"].tolist()
    n_la, n_ra = len(la), len(ra)

    M = df.pivot_table(index="family_id", columns="cell_type",
                          values="dScore").reindex(
                              index=fam_order, columns=CT_ORDER)
    Mp = df.pivot_table(index="family_id", columns="cell_type",
                            values="p_perm").reindex(
                                index=fam_order, columns=CT_ORDER)
    Mq = df.pivot_table(index="family_id", columns="cell_type",
                            values="q_perm").reindex(
                                index=fam_order, columns=CT_ORDER)

    K = len(fam_order); C = len(CT_ORDER)
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(185 * MM, 115 * MM))
    grid = fig.add_gridspec(
        nrows=2, ncols=4,
        height_ratios=[2.6, 18],
        width_ratios=[1.0, 18, 0.4, 0.6],
        hspace=0.04, wspace=0.05,
        left=0.27, right=0.94, top=0.95, bottom=0.10)

    ax_top  = fig.add_subplot(grid[0, 1])    # top labels (cell types)
    ax_main = fig.add_subplot(grid[1, 1])    # heatmap
    ax_dir  = fig.add_subplot(grid[1, 0])    # LA / RA direction badge
    ax_cbar = fig.add_subplot(grid[1, 3])    # colorbar

    # colour map
    vmax = float(np.nanmax(np.abs(M.values)))
    vmax = max(0.005, round(vmax, 3))
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "BrBG_div", [(0.0, BRBG_TEAL), (0.5, "#FFFFFF"),
                       (1.0, BRBG_BROWN)])
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    for i, fid in enumerate(fam_order):
        for j, ct in enumerate(CT_ORDER):
            v = M.iloc[i, j]
            p = Mp.iloc[i, j]
            q = Mq.iloc[i, j]
            facec = "#F2F2F2" if pd.isna(v) else cmap(norm(v))
            rect = mpatches.Rectangle(
                (j, i), 1.0, 1.0, facecolor=facec,
                edgecolor="white", linewidth=0.6)
            ax_main.add_patch(rect)
            if pd.isna(v):
                continue
            # text contrast: pick black or white based on cell value
            shade = abs(v) / vmax
            txt_color = "white" if shade > 0.65 else "black"
            # numeric Δ centered (3 decimals; small font)
            ax_main.text(j + 0.5, i + 0.5, f"{v:+.3f}",
                            ha="center", va="center",
                            fontsize=6.0, color=txt_color, zorder=3)
            # significance marker: small filled badge at top-right
            # of the cell so it sits clearly INSIDE the cell border.
            if not pd.isna(p) and p < 0.10:
                if p < 0.05:
                    badge_face, badge_text = "#000000", "*"
                    text_color = "white"
                    fs = 7.5; fw = "bold"
                else:
                    badge_face, badge_text = "#FFFFFF", "."
                    text_color = "black"
                    fs = 8.5; fw = "bold"
                ax_main.add_patch(mpatches.Circle(
                    (j + 0.85, i + 0.18), 0.10,
                    facecolor=badge_face,
                    edgecolor="black", linewidth=0.4, zorder=4))
                ax_main.text(j + 0.85, i + 0.20, badge_text,
                                ha="center", va="center",
                                fontsize=fs, fontweight=fw,
                                color=text_color, zorder=5)

    ax_main.set_xlim(0, C); ax_main.set_ylim(K, 0)
    ax_main.set_xticks([]); ax_main.set_yticks([])
    for s in ax_main.spines.values():
        s.set_visible(True); s.set_linewidth(0.4); s.set_color("#999")

    # Two-row design: a thin coloured tick at column centre at the
    # bottom, and the rotated label anchored above so the text sits in
    # clean white space without crowding the heatmap top edge.
    for j, ct in enumerate(CT_ORDER):
        ax_top.plot([j + 0.5, j + 0.5], [0.0, 0.10],
                      color=GREY_MID, lw=0.7, zorder=2)
        ax_top.text(j + 0.5, 0.16, CT_LABELS[ct],
                      ha="left", va="bottom",
                      rotation=45, rotation_mode="anchor",
                      fontsize=7.0, fontweight="bold",
                      color="black")
    ax_top.set_xlim(0, C); ax_top.set_ylim(0, 1); ax_top.axis("off")

    for i, fid in enumerate(fam_order):
        d = "LA" if i < n_la else "RA"
        col = LA_COLOR if d == "LA" else RA_COLOR
        ax_dir.add_patch(mpatches.Rectangle(
            (0, i), 1.0, 1.0, facecolor=col, alpha=0.18,
            edgecolor="none"))
        ax_dir.text(0.5, i + 0.5, d,
                      ha="center", va="center",
                      fontsize=6.5, fontweight="bold", color=col)
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(K, 0); ax_dir.axis("off")

    for i, fid in enumerate(fam_order):
        meta = fam_meta[fam_meta["family_id"] == fid].iloc[0]
        lbl = f"F{int(fid)}  {meta['short']}  (n={int(meta['n_programs'])})"
        d = "LA" if i < n_la else "RA"
        col = LA_COLOR if d == "LA" else RA_COLOR
        ax_dir.text(-0.5, i + 0.5, lbl,
                      ha="right", va="center",
                      fontsize=6.5, color=col,
                      fontweight="normal",
                      transform=ax_dir.transData,
                      clip_on=False)

    cb = mpl.colorbar.ColorbarBase(
        ax_cbar, cmap=cmap, norm=norm, orientation="vertical")
    cb.set_label("Δ singscore  (AF − Ctrl)", fontsize=6.8)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)

    fig.suptitle(
        f"Fig 3B  Family × cell-type AF rewiring shift "
        f"(LA snRNA arm, GSE255612)",
        fontsize=8.5, fontweight="bold", y=0.97, x=0.02, ha="left")

    fig.text(0.02, 0.025,
              "Cells: family-mean Δ singscore (AF − Ctrl) over LA donors.   "
              "Markers: * p_perm < 0.05,  . p_perm < 0.10  (10,000 MC perms)\n"
              "Direction band: LA-biased families (blue, top); "
              "RA-biased families (red, bottom).",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
