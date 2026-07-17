#!/usr/bin/env python3
"""N2c — Reversal-gene panel sensitivity analysis."""
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

LONG_TSV = ROOT / "results" / "n2_chamber_reversal_long.tsv"
OUT_TSV  = ROOT / "results" / "n2c_panel_sensitivity_grid.tsv"
FIG_OUT  = ROOT / "figures" / "figS_n2c_panel_sensitivity.pdf"

# Default panel anchor for highlighting in heatmap
DEFAULT_MAG_HEALTHY = 0.30
DEFAULT_MAG_AF      = 0.15
DEFAULT_MIN_COHORTS = 2

# Sweep grid
MAG_HEALTHY_GRID = [0.20, 0.25, 0.30, 0.35, 0.40, 0.50]
MAG_AF_GRID      = [0.10, 0.12, 0.15, 0.20, 0.25, 0.30]
MIN_COHORTS_GRID = [1, 2, 3]


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

    long = pd.read_csv(LONG_TSV, sep="\t")
    print(f"loaded long: {len(long)} rows, "
            f"{long['gene'].nunique()} genes × {long['cohort'].nunique()} cohorts")
    cohorts = sorted(long["cohort"].unique().tolist())
    print(f"cohorts: {cohorts}")

    # absolute aCAI values per row
    long["abs_h"] = long["aCAI_healthy"].abs()
    long["abs_af"] = long["aCAI_AF"].abs()

    # Pre-compute af_specific (sign-only, threshold-invariant)
    long["af_specific_recheck"] = long["af_specific"].astype(bool)

    rows = []
    for mh in MAG_HEALTHY_GRID:
        for ma in MAG_AF_GRID:
            # per (gene, cohort): af_specific AND above_magnitude
            long["above"] = ((long["abs_h"] >= mh)
                              & (long["abs_af"] >= ma))
            long["call"] = long["af_specific_recheck"] & long["above"]
            # per gene: number of cohorts called
            ncohorts = long.groupby("gene")["call"].sum()
            for mc in MIN_COHORTS_GRID:
                n_panel = int((ncohorts >= mc).sum())
                # decompose by direction
                kept_genes = ncohorts[ncohorts >= mc].index
                if len(kept_genes) > 0:
                    sub = long[long["gene"].isin(kept_genes)]
                    dir_per_gene = sub.groupby("gene")["healthy_dir"].first()
                    n_la = int((dir_per_gene == 1).sum())
                    n_ra = int((dir_per_gene == -1).sum())
                else:
                    n_la = n_ra = 0
                rows.append({
                    "MAG_HEALTHY": mh,
                    "MAG_AF":      ma,
                    "MIN_COHORTS": mc,
                    "n_panel":     n_panel,
                    "n_RA_to_LA":  n_ra,    # RA-biased gene → LA in AF
                    "n_LA_to_RA":  n_la,
                })
    grid = pd.DataFrame(rows)
    grid.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(grid)} cells)")

    print("\n=== panel size at default thresholds ===")
    default = grid[(grid["MAG_HEALTHY"] == DEFAULT_MAG_HEALTHY)
                       & (grid["MAG_AF"] == DEFAULT_MAG_AF)
                       & (grid["MIN_COHORTS"] == DEFAULT_MIN_COHORTS)]
    print(default.to_string(index=False))

    print("\n=== panel size sweep at default MAG_HEALTHY = 0.30, MAG_AF = 0.15 ===")
    sweep = grid[(grid["MAG_HEALTHY"] == DEFAULT_MAG_HEALTHY)
                   & (grid["MAG_AF"] == DEFAULT_MAG_AF)]
    print(sweep[["MIN_COHORTS", "n_panel", "n_RA_to_LA", "n_LA_to_RA"]] \
              .to_string(index=False))

    fig = plt.figure(figsize=(190 / 25.4, 70 / 25.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.05],
                            wspace=0.30,
                            left=0.07, right=0.94, top=0.85, bottom=0.18)

    fig.suptitle(
        "Fig S  Chamber-reversal panel size sensitivity grid  "
        "(MAG_HEALTHY × MAG_AF × MIN_COHORTS)",
        fontsize=8.8, fontweight="bold", y=0.985, x=0.02, ha="left")

    # shared color scale across all 3 panels
    vmax = float(grid["n_panel"].max())
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "sequential",
        ["#FFFFFF", "#FDE5D6", "#FBB089", "#E36F45", "#A03029"])

    for k, mc in enumerate(MIN_COHORTS_GRID):
        ax = fig.add_subplot(gs[k])
        sub = grid[grid["MIN_COHORTS"] == mc]
        M = sub.pivot(index="MAG_AF", columns="MAG_HEALTHY",
                          values="n_panel")
        # ensure ordering
        M = M.reindex(index=MAG_AF_GRID, columns=MAG_HEALTHY_GRID)

        im = ax.imshow(M.values, cmap=cmap,
                          vmin=0, vmax=vmax,
                          aspect="auto", origin="lower",
                          interpolation="nearest")
        # numeric overlay
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = int(M.values[i, j])
                shade = v / vmax if vmax > 0 else 0
                col = "white" if shade > 0.55 else "black"
                ax.text(j, i, str(v), ha="center", va="center",
                          fontsize=6.0, color=col, fontweight="bold")
        # highlight the default cell
        if mc == DEFAULT_MIN_COHORTS:
            j_def = MAG_HEALTHY_GRID.index(DEFAULT_MAG_HEALTHY)
            i_def = MAG_AF_GRID.index(DEFAULT_MAG_AF)
            ax.add_patch(mpatches.Rectangle(
                (j_def - 0.5, i_def - 0.5), 1, 1,
                facecolor="none", edgecolor="#1B7837", linewidth=1.6,
                zorder=5))
            ax.text(j_def, i_def + 0.40, "★ default",
                      ha="center", va="bottom",
                      fontsize=6.0, color="#1B7837",
                      fontweight="bold", zorder=6)

        ax.set_xticks(range(len(MAG_HEALTHY_GRID)))
        ax.set_xticklabels([f"{x:.2f}" for x in MAG_HEALTHY_GRID],
                              fontsize=6.5, rotation=0)
        ax.set_yticks(range(len(MAG_AF_GRID)))
        ax.set_yticklabels([f"{x:.2f}" for x in MAG_AF_GRID],
                              fontsize=6.5)
        ax.set_xlabel("MAG_HEALTHY", fontsize=7)
        if k == 0:
            ax.set_ylabel("MAG_AF", fontsize=7)
        ax.set_title(f"MIN_COHORTS = {mc}",
                       loc="left", fontsize=7.5, fontweight="bold", pad=4)

    # shared colorbar
    cax = fig.add_subplot(gs[3])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                       norm=mcolors.Normalize(vmin=0, vmax=vmax),
                                       orientation="vertical")
    cb.set_label("n  panel genes", fontsize=6.5)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)

    fig.text(0.02, 0.025,
              "Cells:  number of reversal-panel genes at each "
              "(MAG_HEALTHY, MAG_AF, MIN_COHORTS) combination.   "
              "Default panel (★) = MAG_HEALTHY = 0.30, MAG_AF = 0.15, "
              f"MIN_COHORTS = 2  →  n = "
              f"{int(default['n_panel'].iloc[0])} genes.   "
              "MIN_COHORTS = 3 collapses the panel  →  no gene satisfies "
              "the AF-specific flip across all 3 paired cohorts.",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
