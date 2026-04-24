#!/usr/bin/env python3
"""
Figure 4D  Program × duration-trend ρ heatmap (full 30 × 4 analyses).

Shows the per-program AF-duration trend across all four GSE41177
analyses:
    laa_primary     LAA AF samples with numeric duration   (n = 12, ρ)
    laa_sens_gt1m   LAA AF samples treating ">1M" = 1 mo   (n = 16, ρ)
    laa_sens_ord    LAA all samples on 3-tier SR/sht/lng   (n = 19, τ)
    lapv_primary    LAPV AF samples with numeric duration  (n = 12, ρ)

Output:  figures/fig4d_program_duration_heatmap.pdf
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

D4_DIR = ROOT / "Project_P4" / "GSE41177" / "results_primary"
FIG_OUT = ROOT / "figures" / "fig4d_program_duration_heatmap.pdf"

LA_COLOR = "#D6604D"
RA_COLOR = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
PRGN_PURPLE = "#762A83"
PRGN_GREEN  = "#1B7837"

ANALYSIS_FILES = {
    "laa_primary":      ("LAA primary\nρ  n = 12",       "duration_trend_laa_primary.tsv"),
    "laa_sens_gt1m":    ("LAA sens A\nρ  n = 16 (>1M ≡ 1)", "duration_trend_laa_sens_gt1m_eq1.tsv"),
    "laa_sens_ord":     ("LAA sens B\nτ  n = 19 (3 tiers)", "duration_trend_laa_sens_ordinal.tsv"),
    "lapv_primary":     ("LAPV primary\nρ  n = 12",       "duration_trend_lapv_sensitivity.tsv"),
}
ANALYSIS_ORDER = list(ANALYSIS_FILES.keys())

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

    parts = []
    for a, (_, fn) in ANALYSIS_FILES.items():
        d = pd.read_csv(D4_DIR / fn, sep="\t")
        d["analysis"] = a
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)

    # ---- order rows by family then program_order
    ordered = (df.drop_duplicates("program_id")
                  [["program_id", "program_label", "program_order",
                    "family_id", "baseline_direction"]]
                  .sort_values(["family_id", "program_order"])
                  .reset_index(drop=True))
    program_order = ordered["program_id"].tolist()
    K = len(program_order); C = len(ANALYSIS_ORDER)

    Mr = df.pivot_table(index="program_id", columns="analysis",
                            values="statistic").reindex(
                                index=program_order, columns=ANALYSIS_ORDER)
    Mp = df.pivot_table(index="program_id", columns="analysis",
                            values="p_perm").reindex(
                                index=program_order, columns=ANALYSIS_ORDER)

    prog_short = {pid: short_pathway_label(pid, 26) for pid in program_order}

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(160 * MM, 175 * MM))
    grid = fig.add_gridspec(
        nrows=2, ncols=4,
        height_ratios=[2.5, 30],
        width_ratios=[1.0, 12, 0.4, 0.6],
        hspace=0.04, wspace=0.05,
        left=0.27, right=0.92, top=0.94, bottom=0.07)
    ax_top  = fig.add_subplot(grid[0, 1])
    ax_main = fig.add_subplot(grid[1, 1])
    ax_fam  = fig.add_subplot(grid[1, 0])
    ax_cbar = fig.add_subplot(grid[1, 3])

    finite_vals = Mr.values[np.isfinite(Mr.values)]
    vmax = float(np.nanmax(np.abs(finite_vals)))
    vmax = max(0.3, round(vmax + 0.05, 1))
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "PRGn", [(0.0, PRGN_PURPLE), (0.5, "#FFFFFF"),
                   (1.0, PRGN_GREEN)])
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    for i, pid in enumerate(program_order):
        for j, a in enumerate(ANALYSIS_ORDER):
            v = Mr.iloc[i, j]
            p = Mp.iloc[i, j]
            facec = "#F2F2F2" if pd.isna(v) else cmap(norm(v))
            ax_main.add_patch(mpatches.Rectangle(
                (j, i), 1.0, 1.0, facecolor=facec,
                edgecolor="white", linewidth=0.4))
            if pd.isna(v):
                continue
            shade = abs(v) / vmax
            txt_color = "white" if shade > 0.65 else "black"
            ax_main.text(j + 0.5, i + 0.5, f"{v:+.2f}",
                            ha="center", va="center",
                            fontsize=6.5, color=txt_color, zorder=3)
            if not pd.isna(p) and p < 0.10:
                if p < 0.05:
                    badge_face, mark, fs = "#000000", "*", 7.0
                    text_color = "white"
                else:
                    badge_face, mark, fs = "#FFFFFF", ".", 8.0
                    text_color = "black"
                ax_main.add_patch(mpatches.Circle(
                    (j + 0.85, i + 0.18), 0.13,
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
        ax_top.text(j + 0.5, 0.18, ANALYSIS_FILES[a][0],
                      ha="center", va="bottom",
                      fontsize=7.0, fontweight="bold")
    ax_top.set_xlim(0, C); ax_top.set_ylim(0, 1); ax_top.axis("off")

    # left family colour bar with F# tags
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

    # left program labels (figure-level via ax_fam transData)
    for i, pid in enumerate(program_order):
        d   = int(ordered.loc[ordered["program_id"] == pid,
                                  "baseline_direction"].iloc[0])
        col = LA_COLOR if d > 0 else RA_COLOR
        badge = "LA" if d > 0 else "RA"
        ax_fam.text(-0.35, i + 0.5,
                      f"{badge}  {prog_short[pid]}",
                      ha="right", va="center",
                      fontsize=6.2, color=col,
                      transform=ax_fam.transData, clip_on=False)

    # colorbar
    cb = mpl.colorbar.ColorbarBase(
        ax_cbar, cmap=cmap, norm=norm, orientation="vertical")
    cb.set_label("ρ  (or τ)  vs AF duration", fontsize=6.8)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)

    fig.suptitle(
        f"Fig 4D  Program × duration-trend statistic across 4 GSE41177 analyses "
        f"({K} NCP_v2 programs)",
        fontsize=8.5, fontweight="bold", y=0.99, x=0.02, ha="left")

    fig.text(0.02, 0.022,
              "Cells: per-program Spearman ρ (or Kendall τ) vs AF duration; "
              "10,000 MC perm test, BH q within each analysis.   "
              "Markers: * p_perm < 0.05, . p_perm < 0.10.   "
              "Grey cells: stat undefined (n < 3).\n"
              "Rows: 30 NCP_v2 programs grouped by axis family (F1–F11; "
              "left colour bar). LA-biased programs in blue, RA-biased in red.",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
