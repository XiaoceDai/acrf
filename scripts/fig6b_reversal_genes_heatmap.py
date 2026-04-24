#!/usr/bin/env python3
"""
Figure 6B  Chamber-identity reversal genes — per-gene per-cohort
aCAI heatmap (healthy / SR / AF).

For each of the 24 framework-discovered reversal genes, show:
    Healthy bulk aCAI      (single column; reference)
    SR aCAI per cohort     (3 columns)
    AF aCAI per cohort     (3 columns)

Sign convention (locked):  + = LA-biased,  - = RA-biased.

A reversal gene is one whose AF aCAI flips sign relative to healthy
(in ≥ 2 cohorts), while SR donors of the same cohorts preserve the
healthy direction.  The flip is exactly what the heatmap shows
visually:  the AF block colour reverses against the healthy block.

Inputs:
    results/n2_chamber_reversal_genes.tsv   (24-gene shortlist + summary)
    results/n2_chamber_reversal_long.tsv    (per-gene per-cohort aCAI)

Output:
    figures/fig6b_reversal_genes_heatmap.pdf
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

GENES_TSV = ROOT / "results" / "n2_chamber_reversal_genes.tsv"
LONG_TSV  = ROOT / "results" / "n2_chamber_reversal_long.tsv"
FIG_OUT   = ROOT / "figures" / "fig6b_reversal_genes_heatmap.pdf"

LA_COLOR  = "#C97A74"
RA_COLOR  = "#7D99B8"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"

COHORT_ORDER = ["GSE128188", "GSE79768", "GSE115574"]


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
    genes  = pd.read_csv(GENES_TSV, sep="\t")
    long   = pd.read_csv(LONG_TSV, sep="\t")

    # restrict long to selected genes
    long = long[long["gene"].isin(genes["gene"])]

    # row order: RA-biased first (largest |aCAI_healthy|), then LA-biased
    ra = genes[genes["healthy_direction"] == "RA-biased"] \
              .sort_values("aCAI_healthy")           # most negative on top
    la = genes[genes["healthy_direction"] == "LA-biased"] \
              .sort_values("aCAI_healthy", ascending=False)
    gene_order = ra["gene"].tolist() + la["gene"].tolist()
    n_ra, n_la = len(ra), len(la)
    K = len(gene_order)

    # build wide matrix: rows = genes, cols = healthy + 3 SR + 3 AF
    cols = ["healthy"] + [f"SR_{c}" for c in COHORT_ORDER] + \
              [f"AF_{c}" for c in COHORT_ORDER]
    M = pd.DataFrame(index=gene_order, columns=cols, dtype=float)

    healthy = genes.set_index("gene")["aCAI_healthy"]
    for g in gene_order:
        M.loc[g, "healthy"] = healthy.get(g, np.nan)

    for _, r in long.iterrows():
        g = r["gene"]; c = r["cohort"]
        if g not in M.index: continue
        if pd.notna(r["aCAI_SR"]):
            M.loc[g, f"SR_{c}"] = float(r["aCAI_SR"])
        if pd.notna(r["aCAI_AF"]):
            M.loc[g, f"AF_{c}"] = float(r["aCAI_AF"])

    C = len(cols)

    # ---- figure layout ------------------------------------------------
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(170 * MM, 175 * MM))
    grid = fig.add_gridspec(
        nrows=3, ncols=4,
        height_ratios=[1.4, 1.4, 24],
        width_ratios=[0.8, 18, 0.4, 0.6],
        hspace=0.04, wspace=0.05,
        left=0.18, right=0.92, top=0.92, bottom=0.10)

    ax_top1 = fig.add_subplot(grid[0, 1])   # block band (Healthy / SR / AF)
    ax_top2 = fig.add_subplot(grid[1, 1])   # cohort band
    ax_main = fig.add_subplot(grid[2, 1])   # heatmap
    ax_dir  = fig.add_subplot(grid[2, 0])   # LA/RA direction
    ax_cbar = fig.add_subplot(grid[2, 3])

    # diverging palette for aCAI
    finite = M.values[np.isfinite(M.values.astype(float))]
    vmax = float(np.nanpercentile(np.abs(finite), 95))
    vmax = max(0.5, round(vmax + 0.05, 1))
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "LA_RA", [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)])
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # ---- main heatmap ----
    for i, g in enumerate(gene_order):
        for j, c in enumerate(cols):
            v = M.iloc[i, j]
            if pd.isna(v):
                facec = "#F2F2F2"
            else:
                facec = cmap(norm(v))
            ax_main.add_patch(mpatches.Rectangle(
                (j, i), 1.0, 1.0,
                facecolor=facec,
                edgecolor="white", linewidth=0.5))
            if pd.isna(v):
                continue
            shade = abs(v) / vmax
            txt_color = "white" if shade > 0.7 else "black"
            ax_main.text(j + 0.5, i + 0.5, f"{v:+.2f}",
                            ha="center", va="center",
                            fontsize=5.5, color=txt_color)

    # block-separator vertical lines
    ax_main.axvline(1.0, color=GREY_DARK, lw=0.6, zorder=4)
    ax_main.axvline(4.0, color=GREY_DARK, lw=0.6, zorder=4)

    # split between RA-biased and LA-biased gene blocks
    ax_main.axhline(n_ra, color="white", lw=2.5, zorder=4)
    ax_main.axhline(n_ra, color=GREY_DARK, lw=0.5, zorder=4.5, ls=":")

    ax_main.set_xlim(0, C); ax_main.set_ylim(K, 0)
    ax_main.set_xticks([]); ax_main.set_yticks([])
    for s in ax_main.spines.values():
        s.set_visible(True); s.set_linewidth(0.4); s.set_color("#999")

    # ---- top block band (Healthy / SR / AF) ----
    block_specs = [
        ("Healthy",  0, 1, "#FFE9B0"),
        ("SR",       1, 4, "#E1ECF4"),
        ("AF",       4, 7, "#F4D9D9"),
    ]
    for lbl, x0, x1, col in block_specs:
        ax_top1.add_patch(mpatches.Rectangle(
            (x0, 0), x1 - x0, 1.0, facecolor=col, edgecolor="white"))
        ax_top1.text((x0 + x1) / 2, 0.5,
                        f"{lbl}  aCAI",
                        ha="center", va="center",
                        fontsize=7.5, fontweight="bold")
    ax_top1.set_xlim(0, C); ax_top1.set_ylim(0, 1); ax_top1.axis("off")

    # ---- top cohort band ----
    cohort_band = ["—"] + COHORT_ORDER + COHORT_ORDER
    for j, c in enumerate(cohort_band):
        ax_top2.add_patch(mpatches.Rectangle(
            (j, 0), 1.0, 1.0,
            facecolor=GREY_LIGHT, edgecolor="white"))
        ax_top2.text(j + 0.5, 0.5, c if c != "—" else "bulk",
                        ha="center", va="center",
                        fontsize=6.0, color=GREY_DARK,
                        fontstyle="italic")
    ax_top2.set_xlim(0, C); ax_top2.set_ylim(0, 1); ax_top2.axis("off")

    # ---- left direction band ----
    for i, g in enumerate(gene_order):
        d = "RA" if i < n_ra else "LA"
        col = RA_COLOR if d == "RA" else LA_COLOR
        ax_dir.add_patch(mpatches.Rectangle(
            (0, i), 1.0, 1.0, facecolor=col, alpha=0.18,
            edgecolor="none"))
        ax_dir.text(0.5, i + 0.5, d,
                      ha="center", va="center",
                      fontsize=5.8, fontweight="bold", color=col)
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(K, 0); ax_dir.axis("off")

    # ---- left gene labels (figure-level via ax_dir) ----
    for i, g in enumerate(gene_order):
        d = "RA" if i < n_ra else "LA"
        col = RA_COLOR if d == "RA" else LA_COLOR
        ax_dir.text(-0.3, i + 0.5, g,
                      ha="right", va="center",
                      fontsize=6.5, color=col,
                      fontweight="bold",
                      transform=ax_dir.transData, clip_on=False)

    # ---- colorbar ----
    cb = mpl.colorbar.ColorbarBase(
        ax_cbar, cmap=cmap, norm=norm, orientation="vertical")
    cb.set_label("aCAI  (LA − RA)", fontsize=6.8)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)

    # ---- title + caption ----
    fig.suptitle(
        f"Fig 6B  Chamber-identity reversal genes  (n = {K}: "
        f"{n_ra} RA→LA  +  {n_la} LA→RA)",
        fontsize=8.8, fontweight="bold", y=0.965, x=0.02, ha="left")
    fig.text(0.02, 0.058,
              "Each cell shows the gene's aCAI = score_LA − score_RA "
              "in that group  (+ = LA-biased, − = RA-biased).   "
              "AF block colour flips against the healthy block;  SR block "
              "preserves the healthy direction → reversal is AF-specific.",
              fontsize=5.8, color="#444", ha="left", va="bottom")
    fig.text(0.02, 0.020,
              "Discovery primitive (N2):  per-gene aCAI computed at donor "
              "level in each of 3 paired cohorts;  shortlisted as 'reversal' "
              "if AF-direction flips against healthy AND SR-direction is "
              "preserved in ≥ 2 of 3 cohorts.",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
