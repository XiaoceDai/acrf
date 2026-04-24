#!/usr/bin/env python3
"""
Figure 2B — Full donor × NCP primary program aCAI heatmap.

Visualises the raw aCAI data underlying Figs 2C / 2D in a single
publication-grade plot:
    rows : 30 NCP primary programs, ordered by axis family (F1 → F16),
            with a left-side family colour bar and a left-side
            LA/RA direction badge
    cols : 49 donors, blocked by cohort (GSE128188 → GSE79768 →
            GSE115574), and within cohort by rhythm (SR | AF) with
            a thin gap between rhythm blocks
    cell : aCAI = score_LA − score_RA  (RA-blue → 0 → LA-red diverging)

Inputs:
    results/fig2_primary_acai_long.tsv
    results/fig2_primary_axis_families.tsv
Output:
    figures/fig2b_donor_acai_heatmap.pdf
"""
from __future__ import annotations

import argparse
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

DEFAULT_ACAI_TSV = ROOT / "results" / "fig2_primary_acai_long.tsv"
DEFAULT_FAM_TSV  = ROOT / "results" / "fig2_primary_axis_families.tsv"
DEFAULT_FIG_OUT  = ROOT / "figures" / "fig2b_donor_acai_heatmap.pdf"

# Style ---------------------------------------------------------------
LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
GREY     = "#9C9C9C"
COHORT_COLORS = {
    "GSE128188": "#A6A6A6",
    "GSE79768":  "#9DBAA9",
    "GSE115574": "#A09BC9",
}
RHY_COLORS = {"SR": "#FFFFFF", "AF": "#333333"}

# Family colours (consistent palette for 11 families)
FAMILY_PALETTE = [
    "#4F86C6", "#6FA8DC", "#9EC4E0",   # LA-leaning blues  (F1 F2 F4)
    "#82B5A1", "#A2C8B0", "#B7D2BC",   # LA-leaning teals  (F3 F5 F6)
    "#D4B58E", "#C8C26C",              # mixed yellows     (F7 F8)
    "#E89C8B", "#D86A6A", "#B94B4B",   # RA-leaning reds   (F9 F10 F11)
]


def sanitize_label(text: str) -> str:
    return (
        text.replace("²", "2")
        .replace("³", "3")
        .replace("⁺", "+")
        .replace("⁻", "-")
        .replace("β", "beta")
        .replace("α", "alpha")
    )


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "axes.titlesize":    9,
        "axes.labelsize":    8,
        "xtick.labelsize":   6,
        "ytick.labelsize":   7,
        "legend.fontsize":   6.5,
        "axes.linewidth":    0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Draw donor x program aCAI heatmap.")
    p.add_argument("--acai-tsv", default=str(DEFAULT_ACAI_TSV),
                   help="Input donor-level aCAI table.")
    p.add_argument("--fam-tsv", default=str(DEFAULT_FAM_TSV),
                   help="Input pathway-family mapping table.")
    p.add_argument("--out-pdf", default=str(DEFAULT_FIG_OUT),
                   help="Output PDF path.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    acai_tsv = Path(args.acai_tsv)
    fam_tsv = Path(args.fam_tsv)
    fig_out = Path(args.out_pdf)
    if not acai_tsv.is_absolute():
        acai_tsv = ROOT / acai_tsv
    if not fam_tsv.is_absolute():
        fam_tsv = ROOT / fam_tsv
    if not fig_out.is_absolute():
        fig_out = ROOT / fig_out

    set_publication_style()

    long_df = pd.read_csv(acai_tsv, sep="\t")
    fam = pd.read_csv(fam_tsv, sep="\t")[["pathway", "family_id"]].rename(
        columns={"pathway": "program_id"})
    long_df = long_df.merge(fam, on="program_id", how="left")
    long_df = long_df.dropna(subset=["family_id"]).copy()
    long_df["family_id"] = long_df["family_id"].astype(int)

    # --- row order : programs sorted by family, then by mean baseline aCAI
    sr_baseline = (long_df[long_df["rhythm"] == "SR"]
                    .groupby(["family_id", "program_id", "direction"])
                    ["aCAI"].mean().reset_index()
                    .rename(columns={"aCAI": "baseline"}))
    sr_baseline = sr_baseline.sort_values(
        ["family_id", "baseline"], ascending=[True, True])
    program_order = sr_baseline["program_id"].tolist()
    program_to_family = dict(
        zip(sr_baseline["program_id"], sr_baseline["family_id"]))
    program_to_direction = dict(
        zip(sr_baseline["program_id"], sr_baseline["direction"]))

    # --- column order : cohort, then SR before AF, then donor_id
    rhythm_order = {"SR": 0, "AF": 1}
    long_df["__rhy"] = long_df["rhythm"].map(rhythm_order)
    cohort_order = ["GSE128188", "GSE79768", "GSE115574"]
    long_df["__coh"] = long_df["cohort"].map(
        {c: i for i, c in enumerate(cohort_order)})

    # Build a unique sample key to handle donor IDs reused across cohorts
    long_df["sample_key"] = long_df["cohort"] + "::" + long_df["donor_id"].astype(str)

    donors = (long_df[["__coh", "cohort", "__rhy", "rhythm",
                          "donor_id", "sample_key"]]
                  .drop_duplicates()
                  .sort_values(["__coh", "__rhy", "donor_id"])
                  .reset_index(drop=True))
    donor_ids = donors["sample_key"].tolist()
    donor_meta = donors.copy()

    pivot = long_df.pivot_table(
        index="program_id", columns="sample_key",
        values="aCAI", aggfunc="mean")
    pivot = pivot.reindex(index=program_order, columns=donor_ids)

    # --- column gaps : insert NaN spacer between rhythm blocks within
    # the same cohort and between cohorts.  We render this via x-tick
    # offsets rather than re-shaping, by tracking column positions.
    col_x = []
    x = 0.0
    prev_coh = None
    prev_rhy = None
    for _, row in donor_meta.iterrows():
        if prev_coh is not None and row["cohort"] != prev_coh:
            x += 1.5  # bigger cohort gap
        elif prev_rhy is not None and row["rhythm"] != prev_rhy:
            x += 0.7  # smaller SR/AF gap
        col_x.append(x)
        x += 1.0
        prev_coh = row["cohort"]
        prev_rhy = row["rhythm"]
    col_x = np.array(col_x)

    # --- figure layout: program labels on the LEFT (outside grid),
    # colour bar on the RIGHT (outside grid), so neither overlaps the
    # heatmap.
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(200 * MM, 150 * MM))
    grid = fig.add_gridspec(
        nrows=3, ncols=5,
        height_ratios=[0.7, 0.5, 30],
        width_ratios=[1.0, 1.2, 30, 0.6, 0.6],
        hspace=0.04, wspace=0.04,
        left=0.27, right=0.92, top=0.92, bottom=0.07)

    ax_top_coh  = fig.add_subplot(grid[0, 2])   # cohort band
    ax_top_rhy  = fig.add_subplot(grid[1, 2])   # rhythm band
    ax_main     = fig.add_subplot(grid[2, 2])   # heatmap
    ax_left_dir = fig.add_subplot(grid[2, 0])   # LA/RA direction badge
    ax_left_fam = fig.add_subplot(grid[2, 1])   # family colour
    ax_cbar     = fig.add_subplot(grid[2, 4])   # colour bar (gap col 3)

    # --- main heatmap (use scatter-like pcolormesh per row)
    vmax = float(np.nanpercentile(np.abs(pivot.values), 95))
    vmax = max(0.3, round(vmax, 1))
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "LA_RA", [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)])
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    K = len(program_order)
    N = len(donor_ids)
    # use imshow with non-uniform x: easier is per-cell rectangles
    for j, did in enumerate(donor_ids):
        col_vals = pivot[did].values
        for i, val in enumerate(col_vals):
            if pd.isna(val):
                continue
            rect = mpatches.Rectangle(
                (col_x[j], i), 1.0, 1.0,
                facecolor=cmap(norm(val)),
                edgecolor="none")
            ax_main.add_patch(rect)
    ax_main.set_xlim(0, col_x[-1] + 1.0)
    ax_main.set_ylim(K, 0)  # invert so first program at top
    ax_main.set_xticks([])
    ax_main.set_yticks([])
    for s in ax_main.spines.values():
        s.set_visible(True); s.set_linewidth(0.4); s.set_color("#999")

    # --- top cohort band
    for cohort in cohort_order:
        idxs = donor_meta.index[donor_meta["cohort"] == cohort].tolist()
        if not idxs:
            continue
        x0 = col_x[idxs[0]]; x1 = col_x[idxs[-1]] + 1.0
        ax_top_coh.add_patch(mpatches.Rectangle(
            (x0, 0), x1 - x0, 1,
            facecolor=COHORT_COLORS[cohort], edgecolor="none"))
        ax_top_coh.text((x0 + x1) / 2, 0.5,
                          f"{cohort} (n={len(idxs)})",
                          ha="center", va="center",
                          fontsize=7, fontweight="bold", color="black")
    ax_top_coh.set_xlim(0, col_x[-1] + 1.0)
    ax_top_coh.set_ylim(0, 1); ax_top_coh.axis("off")

    # --- top rhythm band  (SR white, AF dark)
    for j, row in donor_meta.iterrows():
        ax_top_rhy.add_patch(mpatches.Rectangle(
            (col_x[j], 0), 1.0, 1,
            facecolor=RHY_COLORS[row["rhythm"]], edgecolor="#999",
            linewidth=0.2))
    # draw SR/AF labels per cohort/rhythm group
    for cohort in cohort_order:
        for rhythm in ["SR", "AF"]:
            idxs = donor_meta.index[(donor_meta["cohort"] == cohort)
                                    & (donor_meta["rhythm"] == rhythm)].tolist()
            if not idxs:
                continue
            xc = (col_x[idxs[0]] + col_x[idxs[-1]] + 1.0) / 2
            text_color = "white" if rhythm == "AF" else "black"
            ax_top_rhy.text(xc, 0.5, rhythm,
                              ha="center", va="center",
                              fontsize=6.0, fontweight="bold",
                              color=text_color)
    ax_top_rhy.set_xlim(0, col_x[-1] + 1.0)
    ax_top_rhy.set_ylim(0, 1); ax_top_rhy.axis("off")

    # --- left family colour band  (one cell per program row)
    fam_order = sorted(set(program_to_family.values()))
    fam_to_color = {fid: FAMILY_PALETTE[k % len(FAMILY_PALETTE)]
                     for k, fid in enumerate(fam_order)}
    for i, pid in enumerate(program_order):
        fid = program_to_family[pid]
        ax_left_fam.add_patch(mpatches.Rectangle(
            (0, i), 1.0, 1.0,
            facecolor=fam_to_color[fid], edgecolor="none"))
    # family labels at the centre of each contiguous family block
    fam_blocks = {}
    for i, pid in enumerate(program_order):
        fid = program_to_family[pid]
        fam_blocks.setdefault(fid, []).append(i)
    for fid, idxs in fam_blocks.items():
        yc = (idxs[0] + idxs[-1] + 1) / 2
        ax_left_fam.text(0.5, yc, f"F{fid}",
                          ha="center", va="center",
                          fontsize=6.5, fontweight="bold",
                          color="white" if fid in {1, 11, 9, 10} else "black")
    ax_left_fam.set_xlim(0, 1)
    ax_left_fam.set_ylim(K, 0); ax_left_fam.axis("off")

    # --- left direction band (LA/RA badge)
    for i, pid in enumerate(program_order):
        d = program_to_direction[pid]
        col = LA_COLOR if d == "LA-biased" else RA_COLOR
        label = "LA" if d == "LA-biased" else "RA"
        ax_left_dir.add_patch(mpatches.Rectangle(
            (0, i), 1.0, 1.0,
            facecolor=col, alpha=0.18, edgecolor="none"))
        ax_left_dir.text(0.5, i + 0.5, label,
                          ha="center", va="center",
                          fontsize=5.5, fontweight="bold", color=col)
    ax_left_dir.set_xlim(0, 1)
    ax_left_dir.set_ylim(K, 0); ax_left_dir.axis("off")

    # --- right colour bar
    cb = mpl.colorbar.ColorbarBase(
        ax_cbar, cmap=cmap, norm=norm, orientation="vertical")
    cb.set_label("aCAI  (LA − RA)", fontsize=7)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)

    # --- program name labels on the LEFT, drawn on the direction
    # band axis at negative x in axes coords so they sit outside.
    short_labels = [sanitize_label(short_pathway_label(p, 30)) for p in program_order]
    for i, lbl in enumerate(short_labels):
        ax_left_dir.text(-0.5, i + 0.5, lbl,
                            ha="right", va="center",
                            fontsize=6.0, color="black",
                            transform=ax_left_dir.transData,
                            clip_on=False)
    ax_main.tick_params(left=False, labelleft=False,
                          right=False, labelright=False,
                          length=0)

    fig.suptitle(
        f"Fig 2B  Locked v3 panel donor-level aCAI for {K} programs across "
        f"three paired AF cohorts (n = {N} donors)",
        fontsize=8.5, fontweight="bold", y=0.985, x=0.02, ha="left")

    fig.text(0.02, 0.025,
              "Rows: programs ordered by axis family (F1–F11) and "
              "baseline aCAI.   Columns: donors blocked by cohort, then "
              "by rhythm (SR | AF; AF = dark band).   "
              "Cells: aCAI = score_LA − score_RA  "
              "(blue = RA bias, red = LA bias).",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    fig_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_out)
    plt.close(fig)
    print(f"wrote {fig_out.relative_to(ROOT)}")
    print(f"     dim: {K} programs × {N} donors    vmax={vmax:.2f}")


if __name__ == "__main__":
    main()
