#!/usr/bin/env python3
"""
Figure 5B  Family-level triangulation heatmap across Fig 2 / 3 / 4.

Three evidence layers per axis family (rows = 11 families):

    L1  Bulk aCAI rewiring  (Fig 2D pooled Hedges' g over 3 cohorts)
    L2  LA snRNA carrier    (Fig 3B signed max-|Δ| over 11 cell types,
                              annotated with the winning cell type)
    L3  Duration trend      (Fig 4B LAA primary Spearman ρ)

Each column carries its own diverging palette range (PRGn-like) since
the three statistics are not directly commensurable; cells are labelled
with the raw value so the reader can calibrate across layers.  Rows are
ordered LA-biased on top, RA-biased on bottom.

Inputs:
    results/fig2d_family_forest.tsv
    results/fig3_primary_family_celltype_delta.tsv
    Project_P4/GSE41177/results_primary/duration_trend_FAMILY_laa_primary.tsv
    results/fig2_primary_axis_families.tsv
    results/fig2_primary_acai_long.tsv

Output:
    figures/fig5b_triangulation_heatmap.pdf
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

FIG2D_TSV = ROOT / "results" / "fig2d_family_forest.tsv"
FIG3B_TSV = ROOT / "results" / "fig3_primary_family_celltype_delta.tsv"
FIG4B_TSV = ROOT / "Project_P4" / "GSE41177" / "results_primary" / \
             "duration_trend_FAMILY_laa_primary.tsv"
FAM_TSV   = ROOT / "results" / "fig2_primary_axis_families.tsv"
ACAI_TSV  = ROOT / "results" / "fig2_primary_acai_long.tsv"
OUT_TSV   = ROOT / "results" / "fig5b_triangulation.tsv"
FIG_OUT   = ROOT / "figures" / "fig5b_triangulation_heatmap.pdf"

LA_COLOR = "#D6604D"
RA_COLOR = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"

PRGN_PURPLE = "#762A83"   # negative
PRGN_GREEN  = "#1B7837"   # positive
BRBG_BROWN  = "#8C510A"
BRBG_TEAL   = "#01665E"

LAYER_LABELS = {
    "bulk_g":  "L1  Fig 2D\nBulk aCAI\npooled  g",
    "snrna_d": "L2  Fig 3B\nLA snRNA\nmax  |Δ|  (signed)",
    "dur_rho": "L3  Fig 4B\nLAA primary\nρ  vs duration",
}
LAYER_ORDER = ["bulk_g", "snrna_d", "dur_rho"]


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

    f2 = pd.read_csv(FIG2D_TSV, sep="\t")
    f2 = f2[f2["cohort"] == "POOLED"][["family_id", "direction", "g",
                                            "ci_lo", "ci_hi"]] \
              .rename(columns={"g": "bulk_g"})

    f3 = pd.read_csv(FIG3B_TSV, sep="\t")
    f3 = f3[f3["cell_type_included"]]
    # Pick each family's most-significant cell type (smallest p_perm)
    # as the representative carrier.  Using max|Δ| is vulnerable to
    # small-sample high-variance compartments (e.g. mast_cell), which
    # produce noise-dominated extremes; p_perm is already donor-count
    # aware via permutation.
    f3_winner = (f3.sort_values(["family_id", "p_perm"],
                                    ascending=[True, True])
                    .groupby("family_id")
                    .head(1)[["family_id", "cell_type", "dScore",
                                "p_perm", "q_perm"]]
                    .rename(columns={"dScore": "snrna_d",
                                      "cell_type": "snrna_top_ct",
                                      "p_perm":    "snrna_p",
                                      "q_perm":    "snrna_q"}))

    f4 = pd.read_csv(FIG4B_TSV, sep="\t")
    f4 = f4[["family_id", "statistic", "p_perm", "q_perm"]] \
              .rename(columns={"statistic": "dur_rho",
                                "p_perm":    "dur_p",
                                "q_perm":    "dur_q"})

    tri = f2.merge(f3_winner, on="family_id", how="outer") \
             .merge(f4,         on="family_id", how="outer")

    # representative pathway labels
    long = pd.read_csv(ACAI_TSV, sep="\t")
    fam_map = pd.read_csv(FAM_TSV, sep="\t")[["pathway", "family_id"]].rename(
        columns={"pathway": "program_id"})
    long = long.merge(fam_map, on="program_id", how="left")
    fam_rep = (long.assign(absA=long["aCAI"].abs())
                    .sort_values("absA", ascending=False)
                    .groupby("family_id")
                    .agg(rep=("program_id", "first"),
                           n_programs=("program_id", "nunique"))
                    .reset_index())
    fam_rep["short"] = fam_rep["rep"].apply(
        lambda s: short_pathway_label(s, 22))
    tri = tri.merge(fam_rep[["family_id", "short", "n_programs"]],
                       on="family_id", how="left")
    tri["direction_label"] = tri["direction"].apply(
        lambda d: "LA" if d == "LA-biased" else "RA")

    # row ordering
    la = tri[tri["direction_label"] == "LA"].sort_values("family_id")
    ra = tri[tri["direction_label"] == "RA"].sort_values("family_id")
    fam_order = la["family_id"].tolist() + ra["family_id"].tolist()
    n_la, n_ra = len(la), len(ra)
    tri = tri.set_index("family_id").loc[fam_order].reset_index()

    tri.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)}  ({len(tri)} families)")

    K = len(fam_order); C = len(LAYER_ORDER)
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(160 * MM, 105 * MM))
    grid = fig.add_gridspec(
        nrows=2, ncols=5,
        height_ratios=[3.0, 18],
        width_ratios=[1.0, 6, 0.6, 6, 0.6],
        hspace=0.04, wspace=0.12,
        left=0.32, right=0.96, top=0.87, bottom=0.14)

    ax_top  = fig.add_subplot(grid[0, 1])
    ax_top2 = fig.add_subplot(grid[0, 3])
    ax_main = fig.add_subplot(grid[1, 1])
    ax_right = fig.add_subplot(grid[1, 3])
    ax_dir  = fig.add_subplot(grid[1, 0])
    ax_cbar1 = fig.add_subplot(grid[1, 2])
    ax_cbar2 = fig.add_subplot(grid[1, 4])

    # ---- left main block: L1 bulk_g + L2 snrna_d  (BrBG shared palette)
    cmap_brbg = mcolors.LinearSegmentedColormap.from_list(
        "BrBG", [(0.0, BRBG_TEAL), (0.5, "#FFFFFF"), (1.0, BRBG_BROWN)])
    # shared vmax
    v12 = max(
        float(np.nanmax(np.abs(tri["bulk_g"]))),
        float(np.nanmax(np.abs(tri["snrna_d"])) * 10),
    ) / 1.0
    # but L1 (g) and L2 (Δ) are on very different scales.  Normalise
    # each layer independently, then display — each column uses its
    # own norm so the colour range makes sense per statistic.
    vmax_g = max(0.5, float(np.nanmax(np.abs(tri["bulk_g"]))) + 0.05)
    vmax_d = max(0.02, float(np.nanmax(np.abs(tri["snrna_d"]))) + 0.005)
    vmax_r = max(0.3, float(np.nanmax(np.abs(tri["dur_rho"]))) + 0.05)

    norms = {
        "bulk_g":  mcolors.TwoSlopeNorm(vmin=-vmax_g, vcenter=0.0, vmax=vmax_g),
        "snrna_d": mcolors.TwoSlopeNorm(vmin=-vmax_d, vcenter=0.0, vmax=vmax_d),
        "dur_rho": mcolors.TwoSlopeNorm(vmin=-vmax_r, vcenter=0.0, vmax=vmax_r),
    }
    cmap_prgn = mcolors.LinearSegmentedColormap.from_list(
        "PRGn", [(0.0, PRGN_PURPLE), (0.5, "#FFFFFF"), (1.0, PRGN_GREEN)])
    cmaps = {
        "bulk_g":  cmap_brbg,
        "snrna_d": cmap_brbg,
        "dur_rho": cmap_prgn,
    }

    layer_to_sig_cols = {
        "bulk_g":  ("ci_lo", "ci_hi", None),    # CI excludes 0 → significant
        "snrna_d": (None, None, "snrna_p"),     # p_perm
        "dur_rho": (None, None, "dur_p"),
    }

    # paint main block L1+L2 on ax_main (side-by-side), and L3 on ax_right
    # We do this by mapping columns 0..1 → ax_main, col 2 → ax_right
    def paint_cell(ax, i, j_local, layer, v, p_or_ci):
        if pd.isna(v):
            facec = "#F2F2F2"
        else:
            facec = cmaps[layer](norms[layer](v))
        ax.add_patch(mpatches.Rectangle(
            (j_local, i), 1.0, 1.0,
            facecolor=facec, edgecolor="white", linewidth=0.6))
        if pd.isna(v):
            return
        shade = abs(v) / (norms[layer].vmax)
        txt_color = "white" if shade > 0.65 else "black"
        # numeric value
        if layer == "bulk_g":
            txt = f"{v:+.2f}"
        elif layer == "snrna_d":
            txt = f"{v:+.3f}"
        else:
            txt = f"{v:+.2f}"
        ax.text(j_local + 0.5, i + 0.5, txt,
                 ha="center", va="center",
                 fontsize=6.5, color=txt_color, zorder=3,
                 fontweight="bold")

        # significance marker
        sig = False
        if layer == "bulk_g":
            ci_lo = tri.loc[tri["family_id"] == fam_order[i], "ci_lo"].iloc[0]
            ci_hi = tri.loc[tri["family_id"] == fam_order[i], "ci_hi"].iloc[0]
            if pd.notna(ci_lo) and pd.notna(ci_hi) and (ci_lo > 0 or ci_hi < 0):
                sig = True
                mark, fs = "*", 7.5
        elif layer == "snrna_d":
            p = tri.loc[tri["family_id"] == fam_order[i], "snrna_p"].iloc[0]
            if pd.notna(p):
                if p < 0.05: sig, mark, fs = True, "*", 7.5
                elif p < 0.10: sig, mark, fs = True, ".", 8.0
        elif layer == "dur_rho":
            p = tri.loc[tri["family_id"] == fam_order[i], "dur_p"].iloc[0]
            if pd.notna(p):
                if p < 0.05: sig, mark, fs = True, "*", 7.5
                elif p < 0.10: sig, mark, fs = True, ".", 8.0

        if sig:
            face = "#000000" if mark == "*" else "#FFFFFF"
            tcol = "white"   if mark == "*" else "black"
            ax.add_patch(mpatches.Circle(
                (j_local + 0.85, i + 0.18), 0.10,
                facecolor=face, edgecolor="black", linewidth=0.4,
                zorder=4))
            ax.text(j_local + 0.85, i + 0.20, mark,
                     ha="center", va="center",
                     fontsize=fs, fontweight="bold",
                     color=tcol, zorder=5)

    # left block (L1 + L2)
    for i, fid in enumerate(fam_order):
        r = tri[tri["family_id"] == fid].iloc[0]
        paint_cell(ax_main, i, 0, "bulk_g",  r["bulk_g"],  None)
        paint_cell(ax_main, i, 1, "snrna_d", r["snrna_d"], None)
    ax_main.set_xlim(0, 2); ax_main.set_ylim(K, 0)
    ax_main.set_xticks([]); ax_main.set_yticks([])
    for s in ax_main.spines.values():
        s.set_visible(True); s.set_linewidth(0.4); s.set_color("#999")

    # right block (L3)
    for i, fid in enumerate(fam_order):
        r = tri[tri["family_id"] == fid].iloc[0]
        paint_cell(ax_right, i, 0, "dur_rho", r["dur_rho"], None)
    ax_right.set_xlim(0, 1); ax_right.set_ylim(K, 0)
    ax_right.set_xticks([]); ax_right.set_yticks([])
    for s in ax_right.spines.values():
        s.set_visible(True); s.set_linewidth(0.4); s.set_color("#999")

    for j, lyr in enumerate(["bulk_g", "snrna_d"]):
        ax_top.plot([j + 0.5, j + 0.5], [0.0, 0.10],
                      color=GREY_MID, lw=0.7)
        ax_top.text(j + 0.5, 0.14, LAYER_LABELS[lyr],
                      ha="center", va="bottom",
                      fontsize=6.8, fontweight="bold")
    ax_top.set_xlim(0, 2); ax_top.set_ylim(0, 1); ax_top.axis("off")

    ax_top2.plot([0.5, 0.5], [0.0, 0.10], color=GREY_MID, lw=0.7)
    ax_top2.text(0.5, 0.14, LAYER_LABELS["dur_rho"],
                   ha="center", va="bottom",
                   fontsize=6.8, fontweight="bold")
    ax_top2.set_xlim(0, 1); ax_top2.set_ylim(0, 1); ax_top2.axis("off")

    for i, fid in enumerate(fam_order):
        r = tri[tri["family_id"] == fid].iloc[0]
        d = r["direction_label"]
        col = LA_COLOR if d == "LA" else RA_COLOR
        ax_dir.add_patch(mpatches.Rectangle(
            (0, i), 1.0, 1.0, facecolor=col, alpha=0.18,
            edgecolor="none"))
        ax_dir.text(0.5, i + 0.5, d,
                      ha="center", va="center",
                      fontsize=6.5, fontweight="bold", color=col)
        lbl = (f"F{int(fid)}  {r['short']}  "
                 f"(n={int(r['n_programs'])})")
        ax_dir.text(-0.5, i + 0.5, lbl,
                      ha="right", va="center",
                      fontsize=6.2, color=col,
                      transform=ax_dir.transData, clip_on=False)
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(K, 0); ax_dir.axis("off")

    # Winning cell type labels to the right of L2 column (inside main
    # block col=1 right edge)
    for i, fid in enumerate(fam_order):
        r = tri[tri["family_id"] == fid].iloc[0]
        if pd.isna(r["snrna_top_ct"]):
            continue
        ct_short = {
            "atrialCM":"CM", "fibroblast":"fib",
            "endothelial":"endo", "lymphatic_endothelial":"lymEC",
            "pericyte":"peri", "vascular_smooth_muscle":"VSMC",
            "adipocyte":"adi", "macrophage":"mac",
            "lymphocyte":"lym", "mast_cell":"mast",
            "neuronal":"neu",
        }.get(r["snrna_top_ct"], r["snrna_top_ct"][:4])
        ax_main.text(2.03, i + 0.5, ct_short,
                        ha="left", va="center",
                        fontsize=5.5, color=GREY_DARK,
                        fontstyle="italic",
                        transform=ax_main.transData, clip_on=False)

    cb1 = mpl.colorbar.ColorbarBase(
        ax_cbar1, cmap=cmap_brbg, norm=norms["bulk_g"],
        orientation="vertical")
    cb1.set_label("L1 / L2  (BrBG)", fontsize=6.5)
    cb1.outline.set_linewidth(0.4)
    cb1.ax.tick_params(labelsize=6, length=2, width=0.4)

    cb2 = mpl.colorbar.ColorbarBase(
        ax_cbar2, cmap=cmap_prgn, norm=norms["dur_rho"],
        orientation="vertical")
    cb2.set_label("L3  (PRGn)", fontsize=6.5)
    cb2.outline.set_linewidth(0.4)
    cb2.ax.tick_params(labelsize=6, length=2, width=0.4)

    fig.suptitle(
        "Fig 5B  Family-level triangulation across three evidence layers",
        fontsize=8.8, fontweight="bold", y=0.965, x=0.02, ha="left")

    fig.text(0.02, 0.025,
              "L1  pooled Hedges' g of bulk aCAI (AF − SR, 3 cohorts; Fig 2D).   "
              "L2  signed max |Δ| of family-mean singscore across 11 LA snRNA cell types "
              "(Fig 3B; winning cell italic right of cell).   "
              "L3  LAA-primary Spearman ρ vs AF duration (Fig 4B).\n"
              "Markers: * p_perm < 0.05 (L2, L3) or 95 % CI excludes 0 (L1);  "
              ".  p_perm < 0.10.   Direction band: LA blue / RA red.",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
