#!/usr/bin/env python3
"""
Figure 4B  Family × duration-trend ρ heatmap.

Bridges Fig 2D (family-level forest of bulk aCAI rewiring) and Fig 4
(duration-dependent intensification on the LA-arm bulk side, GSE41177).
Each cell of the heatmap is the Spearman ρ (or Kendall τ) between
family-mean singscore and AF duration in the named analysis.

Inputs:
    Project_P4/GSE41177/results_primary/duration_trend_FAMILY_laa_primary.tsv
    Project_P4/GSE41177/results_primary/duration_trend_FAMILY_lapv_sensitivity.tsv
    results/fig2_primary_axis_families.tsv
    results/fig2_primary_acai_long.tsv

Output:
    figures/fig4b_family_duration_heatmap.pdf
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

LAA_TSV  = ROOT / "Project_P4" / "GSE41177" / "results_primary" / \
            "duration_trend_FAMILY_laa_primary.tsv"
LAPV_TSV = ROOT / "Project_P4" / "GSE41177" / "results_primary" / \
            "duration_trend_FAMILY_lapv_sensitivity.tsv"
FAM_TSV  = ROOT / "results" / "fig2_primary_axis_families.tsv"
ACAI_TSV = ROOT / "results" / "fig2_primary_acai_long.tsv"
FIG_OUT  = ROOT / "figures" / "fig4b_family_duration_heatmap.pdf"

LA_COLOR = "#D6604D"
RA_COLOR = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
PRGN_PURPLE = "#762A83"   # diverging palette for ρ (PRGn)
PRGN_GREEN  = "#1B7837"

ANALYSIS_LABELS = {
    "laa_primary":  "LAA  ρ vs duration\n(primary, n = 12)",
    "lapv_primary": "LAPV  ρ vs duration\n(sensitivity, n = 12)",
}
ANALYSIS_ORDER = ["laa_primary", "lapv_primary"]


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

    laa  = pd.read_csv(LAA_TSV, sep="\t").assign(analysis="laa_primary")
    lapv = pd.read_csv(LAPV_TSV, sep="\t").assign(analysis="lapv_primary")
    df = pd.concat([laa, lapv], ignore_index=True)

    # representative pathway per family
    long = pd.read_csv(ACAI_TSV, sep="\t")
    fam = pd.read_csv(FAM_TSV, sep="\t")[["pathway", "family_id"]].rename(
        columns={"pathway": "program_id"})
    long = long.merge(fam, on="program_id", how="left")
    fam_rep = (long.assign(absA=long["aCAI"].abs())
                    .sort_values("absA", ascending=False)
                    .groupby("family_id")
                    .agg(rep=("program_id", "first"))
                    .reset_index())
    fam_rep["short"] = fam_rep["rep"].apply(
        lambda s: short_pathway_label(s, 22))
    df = df.merge(fam_rep[["family_id", "short"]], on="family_id",
                    how="left")

    # row order: LA-biased first, then RA-biased
    fam_meta = df.drop_duplicates("family_id") \
                  [["family_id", "direction_label", "n_programs", "short"]] \
                  .sort_values("family_id")
    la = fam_meta[fam_meta["direction_label"] == "LA-biased"]
    ra = fam_meta[fam_meta["direction_label"] == "RA-biased"]
    fam_order = la["family_id"].tolist() + ra["family_id"].tolist()
    n_la, n_ra = len(la), len(ra)

    # pivots
    Mr = df.pivot_table(index="family_id", columns="analysis",
                            values="statistic").reindex(
                                index=fam_order, columns=ANALYSIS_ORDER)
    Mp = df.pivot_table(index="family_id", columns="analysis",
                            values="p_perm").reindex(
                                index=fam_order, columns=ANALYSIS_ORDER)

    K = len(fam_order); C = len(ANALYSIS_ORDER)
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(140 * MM, 110 * MM))
    grid = fig.add_gridspec(
        nrows=2, ncols=4,
        height_ratios=[2.5, 18],
        width_ratios=[1.0, 8, 0.4, 0.6],
        hspace=0.04, wspace=0.05,
        left=0.34, right=0.92, top=0.93, bottom=0.13)

    ax_top  = fig.add_subplot(grid[0, 1])
    ax_main = fig.add_subplot(grid[1, 1])
    ax_dir  = fig.add_subplot(grid[1, 0])
    ax_cbar = fig.add_subplot(grid[1, 3])

    vmax = float(np.nanmax(np.abs(Mr.values)))
    vmax = max(0.2, round(vmax + 0.05, 1))
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "PRGn", [(0.0, PRGN_PURPLE), (0.5, "#FFFFFF"),
                   (1.0, PRGN_GREEN)])
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    for i, fid in enumerate(fam_order):
        for j, a in enumerate(ANALYSIS_ORDER):
            v = Mr.iloc[i, j]
            p = Mp.iloc[i, j]
            facec = "#F2F2F2" if pd.isna(v) else cmap(norm(v))
            ax_main.add_patch(mpatches.Rectangle(
                (j, i), 1.0, 1.0, facecolor=facec,
                edgecolor="white", linewidth=0.6))
            if pd.isna(v):
                continue
            shade = abs(v) / vmax
            txt_color = "white" if shade > 0.65 else "black"
            ax_main.text(j + 0.5, i + 0.5, f"{v:+.2f}",
                            ha="center", va="center",
                            fontsize=8.0, color=txt_color, zorder=3,
                            fontweight="bold")
            if not pd.isna(p) and p < 0.10:
                if p < 0.05:
                    badge_face, mark, fs = "#000000", "*", 8.0
                    text_color = "white"
                else:
                    badge_face, mark, fs = "#FFFFFF", ".", 9.0
                    text_color = "black"
                ax_main.add_patch(mpatches.Circle(
                    (j + 0.85, i + 0.18), 0.10,
                    facecolor=badge_face, edgecolor="black",
                    linewidth=0.4, zorder=4))
                ax_main.text(j + 0.85, i + 0.20, mark,
                                ha="center", va="center",
                                fontsize=fs, fontweight="bold",
                                color=text_color, zorder=5)

    ax_main.set_xlim(0, C); ax_main.set_ylim(K, 0)
    ax_main.set_xticks([]); ax_main.set_yticks([])
    for s in ax_main.spines.values():
        s.set_visible(True); s.set_linewidth(0.4); s.set_color("#999")

    # top labels
    for j, a in enumerate(ANALYSIS_ORDER):
        ax_top.plot([j + 0.5, j + 0.5], [0.0, 0.10],
                      color=GREY_MID, lw=0.7)
        ax_top.text(j + 0.5, 0.18, ANALYSIS_LABELS[a],
                      ha="center", va="bottom",
                      fontsize=7.0, fontweight="bold")
    ax_top.set_xlim(0, C); ax_top.set_ylim(0, 1); ax_top.axis("off")

    # left direction badge + family labels (on figure-level outside dir)
    for i, fid in enumerate(fam_order):
        d = "LA" if i < n_la else "RA"
        col = LA_COLOR if d == "LA" else RA_COLOR
        ax_dir.add_patch(mpatches.Rectangle(
            (0, i), 1.0, 1.0, facecolor=col, alpha=0.18,
            edgecolor="none"))
        ax_dir.text(0.5, i + 0.5, d,
                      ha="center", va="center",
                      fontsize=6.5, fontweight="bold", color=col)
        meta = fam_meta[fam_meta["family_id"] == fid].iloc[0]
        lbl = f"F{int(fid)}  {meta['short']}  (n={int(meta['n_programs'])})"
        ax_dir.text(-0.5, i + 0.5, lbl,
                      ha="right", va="center",
                      fontsize=6.5, color=col,
                      transform=ax_dir.transData, clip_on=False)
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(K, 0); ax_dir.axis("off")

    # colorbar
    cb = mpl.colorbar.ColorbarBase(
        ax_cbar, cmap=cmap, norm=norm, orientation="vertical")
    cb.set_label("ρ  (or τ)  vs AF duration", fontsize=6.8)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)

    fig.suptitle(
        "Fig 4B  Family-level AF-duration trend  "
        "(GSE41177; LAA primary  &  LAPV sensitivity)",
        fontsize=8.5, fontweight="bold", y=0.985, x=0.02, ha="left")

    fig.text(0.02, 0.025,
              "Cells: Spearman ρ between family-mean singscore and AF duration "
              "in months  (Kendall τ used only for ordinal sensitivity).   "
              "Markers: * p_perm < 0.05,  . p_perm < 0.10  (10,000 MC perms).\n"
              "Direction band: LA-biased families top (blue), "
              "RA-biased bottom (red).",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
