#!/usr/bin/env python3
"""
Fig S3 v3 - supplementary full program x compartment matrix (LA snRNA).

Two heatmaps stacked:
    Top    family  x compartment  (16 x 9, same as Fig 3b but with
                                   adaptive vmin/vmax, no censoring)
    Bottom program x compartment  (30 x 9, fine-grained view)

Shows trending hits (p < 0.10 dot) and significant hits (p < 0.05 box +
asterisk) for full transparency; not constrained to the +/- 0.04 range
used in main Fig 3b.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label  # noqa: E402

PROG_TSV = ROOT / "results" / "fig3_v3_program_celltype_delta.tsv"
FAM_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
ATLAS_TSV = ROOT / "results" / "fig_integrated_af_rewiring_atlas.tsv"
OUT_PDF = ROOT / "figures" / "figS3_supp_fullmatrix.pdf"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"

COMPARTMENTS = [
    "atrialCM", "fibroblast", "endothelial", "macrophage", "lymphocyte",
    "pericyte", "vascular_smooth_muscle", "adipocyte", "neuronal",
]
COMP_SHORT = {
    "atrialCM": "aCM", "fibroblast": "Fibro", "endothelial": "Endo",
    "macrophage": "Macro", "lymphocyte": "Lymph", "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC", "adipocyte": "Adipo", "neuronal": "Neuro",
}


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def draw_heatmap(ax, mat, pmat, vmax, row_labels, col_labels, ylabel):
    cmap = LinearSegmentedColormap.from_list(
        "AF_d_supp", [(0.0, RA_DEEP), (0.5, "#FFFFFF"), (1.0, LA_DEEP)]
    )
    cmap.set_bad("#E0E0E0")
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    masked = np.ma.masked_invalid(mat)
    n_rows, n_cols = mat.shape
    im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto",
                   extent=[-0.5, n_cols - 0.5, n_rows - 0.5, -0.5])
    for i in range(n_rows):
        for j in range(n_cols):
            p = pmat[i, j]
            if not np.isfinite(p):
                continue
            if p < 0.05:
                ax.add_patch(mpatches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    facecolor="none", edgecolor="#1A1A1A", linewidth=1.0,
                ))
                ax.text(j, i + 0.10, "*", ha="center", va="center",
                        fontsize=9, fontweight="bold", color="#1A1A1A")
            elif p < 0.10:
                ax.text(j, i + 0.10, ".", ha="center", va="center",
                        fontsize=9, fontweight="bold", color="#444444")
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(col_labels, rotation=30, ha="right", fontsize=6.6)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels, fontsize=5.4)
    ax.set_ylabel(ylabel, fontsize=7.4)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def main():
    set_style()
    fam = pd.read_csv(FAM_TSV, sep="\t")
    prog = pd.read_csv(PROG_TSV, sep="\t")
    hca = pd.read_csv(HCA_TSV, sep="\t")
    canonical = pd.read_csv(ROOT / "results" / "ncp_v3_canonical_display_order.tsv", sep="\t")
    labels = pd.read_csv(ROOT / "results" / "ncp_v3_family_labels.tsv", sep="\t")

    fam_order = (
        canonical[["family_display_order", "family_id"]]
        .drop_duplicates("family_id")
        .sort_values("family_display_order")["family_id"]
        .astype(int)
        .tolist()
    )

    # FAMILY matrix
    fam_pivot = fam.pivot(index="family_id", columns="compartment", values="dScore")
    fam_pivot = fam_pivot.reindex(index=fam_order, columns=COMPARTMENTS)
    fam_pmat = fam.pivot(index="family_id", columns="compartment", values="p_perm")
    fam_pmat = fam_pmat.reindex(index=fam_order, columns=COMPARTMENTS)

    fam_short_lookup = {
        int(r.family_id): short_pathway_label(str(r.family_label), 36)
        for r in labels.itertuples(index=False)
    }
    fam_labels = [f"F{fid}  {fam_short_lookup.get(fid, '')}" for fid in fam_order]

    # PROGRAM matrix - sort by family order, then program_order within family
    prog["family_id"] = prog["family_id"].astype(int)
    prog["program_order"] = prog["program_order"].astype(int)
    prog["fam_rank"] = prog["family_id"].map({fid: i for i, fid in enumerate(fam_order)}).fillna(99).astype(int)
    prog = prog.sort_values(["fam_rank", "program_order"])
    prog["program_short"] = prog["program_id"].map(lambda s: short_pathway_label(str(s), 36))
    prog_order = (
        prog.drop_duplicates("program_id")
        .sort_values(["fam_rank", "program_order"])["program_id"]
        .tolist()
    )
    prog_pivot = prog.pivot(index="program_id", columns="compartment", values="dScore").reindex(
        index=prog_order, columns=COMPARTMENTS
    )
    prog_pmat = prog.pivot(index="program_id", columns="compartment", values="p_perm").reindex(
        index=prog_order, columns=COMPARTMENTS
    )
    prog_label_lookup = dict(zip(prog["program_id"], prog["program_short"]))
    prog_fam_lookup = dict(zip(prog["program_id"], prog["family_id"]))
    prog_labels = [
        f"F{prog_fam_lookup[p]}  {prog_label_lookup[p]}" for p in prog_order
    ]

    # vmax
    fam_arr = fam_pivot.to_numpy(dtype=float)
    prog_arr = prog_pivot.to_numpy(dtype=float)
    vmax_fam = max(0.04, float(np.nanmax(np.abs(fam_arr))) * 1.05)
    vmax_prog = max(0.04, float(np.nanmax(np.abs(prog_arr))) * 1.05)

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 280 * MM))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[10, 0.6], height_ratios=[2.6, 5.0],
        left=0.30, right=0.97, top=0.94, bottom=0.05, wspace=0.04, hspace=0.20,
    )
    ax_fam = fig.add_subplot(gs[0, 0])
    ax_fam_cb = fig.add_subplot(gs[0, 1])
    ax_prog = fig.add_subplot(gs[1, 0])
    ax_prog_cb = fig.add_subplot(gs[1, 1])

    fam_pmat_arr = fam_pmat.to_numpy(dtype=float)
    prog_pmat_arr = prog_pmat.to_numpy(dtype=float)

    im_fam = draw_heatmap(
        ax_fam, fam_arr, fam_pmat_arr, vmax_fam,
        fam_labels, [COMP_SHORT[c] for c in COMPARTMENTS],
        "16 NCP axis families  (atlas-priority order)",
    )
    cbar_fam = fig.colorbar(im_fam, cax=ax_fam_cb, orientation="vertical")
    cbar_fam.set_label("\u0394singscore  (AF \u2212 Ctrl)", fontsize=6.4)
    cbar_fam.ax.tick_params(labelsize=5.6)

    im_prog = draw_heatmap(
        ax_prog, prog_arr, prog_pmat_arr, vmax_prog,
        prog_labels, [COMP_SHORT[c] for c in COMPARTMENTS],
        "30 NCP pathway programs  (within-family order)",
    )
    cbar_prog = fig.colorbar(im_prog, cax=ax_prog_cb, orientation="vertical")
    cbar_prog.set_label("\u0394singscore  (AF \u2212 Ctrl)", fontsize=6.4)
    cbar_prog.ax.tick_params(labelsize=5.6)

    fig.suptitle(
        "Fig S3  Supplementary  -  full family- and program-level AF \u0394singscore matrices "
        "(LA snRNA, GSE255612)",
        x=0.30, y=0.97, ha="left", fontsize=8.6, fontweight="bold",
    )
    fig.text(
        0.30, 0.025,
        "Adaptive color scale (no +/- 0.04 censoring as in main Fig 3b).  "
        "Boxed cells with * mark p < 0.05; dots mark p < 0.10 (donor-permutation, n_perm = 10000).  "
        "BH-q within each compartment is reported in the source TSV but not visually censored here, "
        "given the limited single-cell power (n = 16 Ctrl + 18 AF).  Use the corresponding "
        "results/fig3_v3_*_celltype_delta.tsv tables for downstream filtering.",
        fontsize=5.7, color="#444444",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
