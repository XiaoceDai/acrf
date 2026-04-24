#!/usr/bin/env python3
"""
Figure 6C  Cellular carriers of chamber-identity reversal genes.

24 reversal genes × 11 LA-snRNA cell types (Project_P3/GSE255612/
pseudobulk_full).  Cell value = sign-aligned Δ expression
(AF − Ctrl in the gene's expected reversal direction);
positive = consistent with chamber-reversal direction in that cell type.

Carrier classification per gene:
    dominant cell type = argmax(delta_aligned) over included compartments
    (only call if delta_aligned > 0).

Inputs:
    results/n2b_primary_reversal_per_celltype_long.tsv      (gene × ct rows)
    results/n2b_primary_reversal_per_celltype_summary.tsv   (per-ct summary)
    results/n2_chamber_reversal_genes.tsv                    (gene metadata)

Output:
    figures/fig6c_reversal_carrier.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

LONG_TSV    = ROOT / "results" / "n2b_primary_reversal_per_celltype_long.tsv"
SUMMARY_TSV = ROOT / "results" / "n2b_primary_reversal_per_celltype_summary.tsv"
GENES_TSV   = ROOT / "results" / "n2_chamber_reversal_genes.tsv"
FIG_OUT     = ROOT / "figures" / "fig6c_reversal_carrier.pdf"

LA_COLOR  = "#D6604D"
RA_COLOR  = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"
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

    long  = pd.read_csv(LONG_TSV, sep="\t")
    summ  = pd.read_csv(SUMMARY_TSV, sep="\t")
    genes = pd.read_csv(GENES_TSV, sep="\t")

    # restrict cell types to the 11 included compartments (matches Fig 3)
    long = long[long["compartment"].isin(CT_ORDER)]
    long = long[long["compartment_included"]]

    # row order: 21 RA→LA (sorted by aCAI_healthy ascending = most negative first)
    # then 3 LA→RA (sorted by aCAI_healthy descending)
    ra = genes[genes["healthy_direction"] == "RA-biased"] \
              .sort_values("aCAI_healthy")
    la = genes[genes["healthy_direction"] == "LA-biased"] \
              .sort_values("aCAI_healthy", ascending=False)
    gene_order = ra["gene"].tolist() + la["gene"].tolist()
    n_ra = len(ra)
    K = len(gene_order)
    C = len(CT_ORDER)

    # pivot
    M = long.pivot_table(index="gene", columns="compartment",
                            values="delta_aligned").reindex(
                                index=gene_order, columns=CT_ORDER)

    # ---- per-gene dominant carrier ----
    carriers = []
    for g in gene_order:
        row = M.loc[g]
        if row.dropna().empty:
            carriers.append(("Unresolved", float("nan"))); continue
        # only positive aligned deltas count as carrier
        pos = row[row > 0]
        if pos.empty:
            carriers.append(("None", float(row.max()))); continue
        ct = pos.idxmax(); val = float(pos.max())
        carriers.append((ct, val))

    # ---- figure layout ------------------------------------------------
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(195 * MM, 175 * MM))
    grid = fig.add_gridspec(
        nrows=2, ncols=5,
        height_ratios=[2.5, 26],
        width_ratios=[1.0, 14, 0.4, 0.6, 4.5],
        hspace=0.04, wspace=0.05,
        left=0.13, right=0.97, top=0.92, bottom=0.07)

    ax_top  = fig.add_subplot(grid[0, 1])
    ax_main = fig.add_subplot(grid[1, 1])
    ax_dir  = fig.add_subplot(grid[1, 0])
    ax_cbar = fig.add_subplot(grid[1, 3])
    ax_carr = fig.add_subplot(grid[1, 4])    # right-side dominant-carrier strip

    # colour map BrBG: positive = aligned with reversal direction
    finite = M.values[np.isfinite(M.values)]
    vmax = float(np.nanpercentile(np.abs(finite), 95))
    vmax = max(0.05, round(vmax + 0.005, 3))
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "BrBG_align",
        [(0.0, BRBG_TEAL), (0.5, "#FFFFFF"), (1.0, BRBG_BROWN)])
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    # ---- main heatmap ----
    for i, g in enumerate(gene_order):
        for j, ct in enumerate(CT_ORDER):
            v = M.iloc[i, j]
            if pd.isna(v):
                facec = "#F2F2F2"
            else:
                facec = cmap(norm(v))
            ax_main.add_patch(mpatches.Rectangle(
                (j, i), 1.0, 1.0, facecolor=facec,
                edgecolor="white", linewidth=0.4))
            if pd.isna(v):
                continue
            shade = abs(v) / vmax
            txt_color = "white" if shade > 0.7 else "black"
            ax_main.text(j + 0.5, i + 0.5, f"{v:+.2f}",
                            ha="center", va="center",
                            fontsize=5.4, color=txt_color)

    # split between RA→LA and LA→RA gene blocks
    ax_main.axhline(n_ra, color="white", lw=2.5, zorder=4)
    ax_main.axhline(n_ra, color=GREY_DARK, lw=0.5, zorder=4.5, ls=":")

    ax_main.set_xlim(0, C); ax_main.set_ylim(K, 0)
    ax_main.set_xticks([]); ax_main.set_yticks([])
    for s in ax_main.spines.values():
        s.set_visible(True); s.set_linewidth(0.4); s.set_color("#999")

    # ---- top labels (cell types) ----
    for j, ct in enumerate(CT_ORDER):
        ax_top.plot([j + 0.5, j + 0.5], [0.0, 0.10],
                      color=GREY_MID, lw=0.7)
        ax_top.text(j + 0.5, 0.16, CT_LABELS[ct],
                      ha="left", va="bottom",
                      rotation=45, rotation_mode="anchor",
                      fontsize=6.8, fontweight="bold")
    ax_top.set_xlim(0, C); ax_top.set_ylim(0, 1); ax_top.axis("off")

    # ---- left direction band ----
    for i, g in enumerate(gene_order):
        d = "RA→LA" if i < n_ra else "LA→RA"
        col = RA_COLOR if i < n_ra else LA_COLOR
        ax_dir.add_patch(mpatches.Rectangle(
            (0, i), 1.0, 1.0, facecolor=col, alpha=0.18,
            edgecolor="none"))
        ax_dir.text(0.5, i + 0.5, d,
                      ha="center", va="center",
                      fontsize=5.5, fontweight="bold", color=col)
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(K, 0); ax_dir.axis("off")

    # left gene labels (figure-level via ax_dir transData)
    for i, g in enumerate(gene_order):
        col = RA_COLOR if i < n_ra else LA_COLOR
        ax_dir.text(-0.25, i + 0.5, g,
                      ha="right", va="center",
                      fontsize=6.0, color=col, fontweight="bold",
                      transform=ax_dir.transData, clip_on=False)

    # ---- right side: dominant-carrier strip ----
    CT_COLOR = {
        "atrialCM": "#E6AB02", "fibroblast": "#66A61E",
        "endothelial": "#7570B3", "lymphatic_endothelial": "#A6761D",
        "pericyte": "#E7298A", "vascular_smooth_muscle": "#1F78B4",
        "adipocyte": "#B2DF8A", "macrophage": "#D95F02",
        "lymphocyte": "#7B3294", "mast_cell": "#FB9A99",
        "neuronal": "#1B9E77", "Unresolved": "#BBBBBB",
        "None": "#DDDDDD",
    }
    for i, (ct, val) in enumerate(carriers):
        col = CT_COLOR.get(ct, "#999999")
        ax_carr.add_patch(mpatches.Rectangle(
            (0, i), 1.0, 1.0, facecolor=col, edgecolor="white",
            linewidth=0.4))
        if ct in CT_LABELS:
            lbl = CT_LABELS[ct]
        else:
            lbl = ct
        ax_carr.text(1.05, i + 0.5,
                        f"{lbl}  Δ={val:+.3f}" if not pd.isna(val) else lbl,
                        ha="left", va="center",
                        fontsize=5.8, color=GREY_DARK,
                        transform=ax_carr.transData, clip_on=False)
    ax_carr.set_xlim(0, 1); ax_carr.set_ylim(K, 0); ax_carr.axis("off")
    # title for carrier strip
    ax_carr.text(0.5, -0.4, "dominant carrier",
                    ha="center", va="bottom",
                    fontsize=6.5, fontweight="bold", color=GREY_DARK,
                    transform=ax_carr.transData, clip_on=False)

    # ---- colorbar ----
    cb = mpl.colorbar.ColorbarBase(
        ax_cbar, cmap=cmap, norm=norm, orientation="vertical")
    cb.set_label("Δ expression  (sign-aligned)", fontsize=6.5)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)

    # ---- title + bottom caption ----
    fig.suptitle(
        "Fig 6C  Cellular carriers of chamber-identity reversal genes  "
        f"(n = {K} genes × {C} LA-snRNA cell types)",
        fontsize=8.8, fontweight="bold", y=0.965, x=0.02, ha="left")

    # per-cell-type bottom summary bar
    ax_bot = fig.add_axes([0.13, 0.025, 0.50, 0.030])
    summ_in = summ[summ["compartment_included"]] \
                  .set_index("compartment").reindex(CT_ORDER)
    bar_h = summ_in["mean_delta_aligned"].to_numpy(dtype=float)
    cols2 = [BRBG_BROWN if v > 0 else BRBG_TEAL for v in bar_h]
    ax_bot.bar(np.arange(C), bar_h, color=cols2, edgecolor="black",
                  linewidth=0.3, width=0.85)
    ax_bot.axhline(0, color=GREY_DARK, lw=0.5)
    ax_bot.set_xticks([])
    ax_bot.set_yticks([])
    ax_bot.set_xlim(-0.5, C - 0.5)
    ax_bot.text(-0.6, 0.5, "mean\nΔaligned",
                  ha="right", va="center",
                  fontsize=5.5, color=GREY_DARK,
                  transform=ax_bot.transData, clip_on=False)
    for s in ax_bot.spines.values():
        s.set_visible(False)

    fig.text(0.02, 0.005,
              "Δ aligned > 0  →  expression in that cell type moves "
              "in the gene's expected reversal direction in AF.   "
              "Bottom bar: per-cell-type mean of Δ aligned across all 24 genes.",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")

    # ---- console summary ----
    print("\n=== gene → dominant carrier ===")
    for g, (ct, val) in zip(gene_order, carriers):
        print(f"  {g:8s}  →  {ct:24s}  Δ_aligned = "
                f"{val:+.4f}" if not pd.isna(val) else f"  {g:8s}  →  {ct}")
    print("\n=== per-cell-type mean Δ_aligned across 24 genes ===")
    print(summ_in[["mean_delta_aligned",
                       "frac_genes_aligned_pos"]].to_string())


if __name__ == "__main__":
    main()
