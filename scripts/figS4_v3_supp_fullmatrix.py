#!/usr/bin/env python3
"""
Fig S4 v3 - Supplementary full Spearman rho matrix vs AF duration.

Two stacked heatmaps:
    Top    16 family x 4 sensitivity codings   (mirror of Fig 4b but adaptive)
    Bottom 30 program x 4 sensitivity codings  (fine-grained)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label  # noqa: E402

PROG_PRIM = ROOT / "results" / "fig4_v3_PROGRAM_laa_primary.tsv"
FAM_PRIM = ROOT / "results" / "fig4_v3_FAMILY_laa_primary.tsv"
FAM_GT1M = ROOT / "results" / "fig4_v3_FAMILY_laa_sens_gt1m_eq1.tsv"
FAM_ORD = ROOT / "results" / "fig4_v3_FAMILY_laa_sens_ordinal.tsv"
FAM_LAPV = ROOT / "results" / "fig4_v3_FAMILY_lapv_sensitivity.tsv"
ATLAS_TSV = ROOT / "results" / "fig_integrated_af_rewiring_atlas.tsv"
OUT_PDF = ROOT / "figures" / "figS4_v3_supp_fullmatrix.pdf"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"

CODINGS = ["primary", "gt1m_eq1", "ordinal", "lapv"]
CODING_LABEL = {
    "primary": "Primary\n(months)",
    "gt1m_eq1": "gt1m_eq1",
    "ordinal": "Ordinal",
    "lapv": "LAPV\n(region)",
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
        "rho_div_supp", [(0.0, RA_DEEP), (0.5, "#FFFFFF"), (1.0, LA_DEEP)]
    )
    cmap.set_bad("#E0E0E0")
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    masked = np.ma.masked_invalid(mat)
    n_rows, n_cols = mat.shape
    im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto",
                    extent=[-0.5, n_cols - 0.5, n_rows - 0.5, -0.5])
    for i in range(n_rows):
        for j in range(n_cols):
            v = mat[i, j]; p = pmat[i, j]
            if not np.isfinite(v):
                continue
            txt_col = "#1A1A1A" if abs(v) < vmax * 0.7 else "white"
            ax.text(j, i + 0.05, f"{v:+.2f}", ha="center", va="center",
                     fontsize=5.4, color=txt_col, fontweight="bold")
            if np.isfinite(p):
                if p < 0.05:
                    ax.add_patch(mpatches.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1,
                        facecolor="none", edgecolor="#FFE066", linewidth=1.4,
                    ))
                    ax.text(j, i - 0.30, "*", ha="center", va="center",
                             fontsize=9, color="#FFE066", fontweight="bold",
                             path_effects=[mpe.Stroke(linewidth=1.2, foreground="#1A1A1A"),
                                              mpe.Normal()])
                elif p < 0.10:
                    ax.text(j, i - 0.30, ".", ha="center", va="center",
                             fontsize=9, color="#1A1A1A", fontweight="bold")
    ax.set_xlim(-0.5, n_cols - 0.5); ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(col_labels, fontsize=6.0, fontweight="bold")
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels, fontsize=5.4)
    ax.set_ylabel(ylabel, fontsize=7.0)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def main():
    set_style()
    fam_dfs = {
        "primary": pd.read_csv(FAM_PRIM, sep="\t").set_index("family_id"),
        "gt1m_eq1": pd.read_csv(FAM_GT1M, sep="\t").set_index("family_id"),
        "ordinal": pd.read_csv(FAM_ORD, sep="\t").set_index("family_id"),
        "lapv": pd.read_csv(FAM_LAPV, sep="\t").set_index("family_id"),
    }
    prog = pd.read_csv(PROG_PRIM, sep="\t")
    atlas = pd.read_csv(ATLAS_TSV, sep="\t").sort_values("PriorityScore", ascending=False)
    fam_order = atlas["family_id"].astype(int).tolist()
    short_lookup = dict(zip(atlas["family_id"].astype(int), atlas["short"].astype(str)))

    # FAMILY matrix
    n_fams = len(fam_order)
    fam_rho = np.full((n_fams, len(CODINGS)), np.nan)
    fam_p = np.full((n_fams, len(CODINGS)), np.nan)
    for i, fid in enumerate(fam_order):
        for j, c in enumerate(CODINGS):
            df = fam_dfs[c]
            if fid in df.index:
                fam_rho[i, j] = float(df.loc[fid, "rho_spearman"])
                fam_p[i, j] = float(df.loc[fid, "p_perm"])
    fam_labels = [f"F{fid}  {short_lookup.get(fid, '')}" for fid in fam_order]

    # PROGRAM matrix - sort by family order, only primary coding (others not computed for programs)
    prog["family_id"] = prog["family_id"].astype(int)
    prog["program_order"] = prog["program_order"].astype(int)
    fam_rank_map = {fid: i for i, fid in enumerate(fam_order)}
    prog["fam_rank"] = prog["family_id"].map(fam_rank_map).fillna(99).astype(int)
    prog = prog.sort_values(["fam_rank", "program_order"]).reset_index(drop=True)
    prog["program_short"] = prog["program_id"].map(lambda s: short_pathway_label(str(s), 36))
    prog_labels = [f"F{int(prog.iloc[i]['family_id'])}  {prog.iloc[i]['program_short']}" for i in range(len(prog))]
    prog_rho_arr = prog[["rho_spearman"]].to_numpy(dtype=float)
    prog_p_arr = prog[["p_perm"]].to_numpy(dtype=float)
    # Pad to 4 columns to match family panel? No, keep 1 column primary only.
    prog_col_labels = ["Primary\n(LAA, months)"]

    vmax_fam = max(0.6, float(np.nanmax(np.abs(fam_rho))) * 1.05)
    vmax_prog = max(0.6, float(np.nanmax(np.abs(prog_rho_arr))) * 1.05)

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 280 * MM))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[10, 0.6], height_ratios=[2.4, 5.5],
        left=0.30, right=0.97, top=0.94, bottom=0.06, wspace=0.04, hspace=0.18,
    )
    ax_fam = fig.add_subplot(gs[0, 0])
    ax_fam_cb = fig.add_subplot(gs[0, 1])
    ax_prog = fig.add_subplot(gs[1, 0])
    ax_prog_cb = fig.add_subplot(gs[1, 1])

    im_fam = draw_heatmap(ax_fam, fam_rho, fam_p, vmax_fam,
                            fam_labels, [CODING_LABEL[c] for c in CODINGS],
                            "16 NCP_v3 families  (atlas-priority order)")
    cbar_fam = fig.colorbar(im_fam, cax=ax_fam_cb, orientation="vertical")
    cbar_fam.set_label("Spearman \u03c1", fontsize=6.4)
    cbar_fam.ax.tick_params(labelsize=5.6)

    im_prog = draw_heatmap(ax_prog, prog_rho_arr, prog_p_arr, vmax_prog,
                             prog_labels, prog_col_labels,
                             "30 NCP_v3 programs  (within-family order; LAA primary only)")
    cbar_prog = fig.colorbar(im_prog, cax=ax_prog_cb, orientation="vertical")
    cbar_prog.set_label("Spearman \u03c1", fontsize=6.4)
    cbar_prog.ax.tick_params(labelsize=5.6)

    fig.suptitle(
        "Fig S4  Supplementary  -  full family- and program-level Spearman \u03c1 vs AF duration "
        "(GSE41177 LAA + LAPV)",
        x=0.30, y=0.97, ha="left", fontsize=8.6, fontweight="bold",
    )
    fig.text(
        0.30, 0.025,
        "Adaptive color scale per panel.  Yellow-bordered cells with * mark p < 0.05 (donor permutation, n_perm = 10,000); dots mark p < 0.10.  "
        "Family panel reports 4 sensitivity codings (Primary numeric months / gt1m_eq1 / Ordinal 3-tier / LAPV region).  "
        "Program panel reports LAA primary only (sensitivity codings not recomputed at program level).  "
        "Use results/fig4_v3_*.tsv tables for downstream filtering and bootstrap CI.",
        fontsize=5.7, color="#444444",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
