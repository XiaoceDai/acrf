#!/usr/bin/env python3
"""
Supplementary heatmap: Variance-Heterogeneity Score (VHS) per family x compartment.

Visualises log2 VHS = log2 [SD_AF(family singscore) / SD_Ctrl(family singscore)]
on GSE255612 LA snRNA pseudobulk.  Positive log2_VHS = AF is more dispersed
between donors (candidate sub-population emergence); negative = AF is more
uniform (cell-intrinsic homogenisation or composition lock-in).

Same layout / sort order as Fig 3b for direct visual comparison.
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

VHS_TSV = ROOT / "results" / "fig3_v3_VHS_per_family_ct.tsv"
ATLAS_TSV = ROOT / "results" / "fig_integrated_af_rewiring_atlas.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
DSCORE_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
OUT_PDF = ROOT / "figures" / "figS3_v3_VHS_overlay.pdf"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
SIG_COLOR = "#FFE066"

COMPARTMENTS = [
    "atrialCM", "fibroblast", "endothelial", "macrophage", "lymphocyte",
    "pericyte", "vascular_smooth_muscle", "adipocyte", "neuronal",
]
COMP_SHORT = {
    "atrialCM": "aCM", "fibroblast": "Fibro", "endothelial": "Endo",
    "macrophage": "Macro", "lymphocyte": "Lymph", "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC", "adipocyte": "Adipo", "neuronal": "Neuro",
}
COMPARTMENT_PALETTE = {
    "atrialCM": "#4E79A7", "fibroblast": "#A0CBE8", "endothelial": "#F28E2B",
    "macrophage": "#59A14F", "lymphocyte": "#8CD17D", "pericyte": "#499894",
    "vascular_smooth_muscle": "#86BCB6", "adipocyte": "#E15759", "neuronal": "#79706E",
}
MIXED_COLOR = "#BDBDBD"
MIXED_LABEL = "Mixed/Unresolved"

# Color scale (purple = AF more uniform / homogenised, orange = AF more dispersed)
VMIN, VMAX = -1.5, 1.5


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def main():
    set_style()
    vhs = pd.read_csv(VHS_TSV, sep="\t")
    atlas = pd.read_csv(ATLAS_TSV, sep="\t").sort_values("PriorityScore", ascending=False)
    hca = pd.read_csv(HCA_TSV, sep="\t")
    dscore = pd.read_csv(DSCORE_TSV, sep="\t")

    fam_order = atlas["family_id"].astype(int).tolist()
    short_lookup = dict(zip(atlas["family_id"].astype(int), atlas["short"].astype(str)))
    dir_lookup = dict(zip(atlas["family_id"].astype(int), atlas["HealthyBias"].astype(str)))
    carriers = dict(zip(hca["family_id"].astype(int), hca["LA_carrier"].astype(str)))

    vhs["family_id"] = vhs["family_id"].astype(int)
    vhs_pivot = vhs.pivot(index="family_id", columns="compartment", values="log2_VHS")
    vhs_pivot = vhs_pivot.reindex(index=fam_order, columns=COMPARTMENTS)
    p_pivot = vhs.pivot(index="family_id", columns="compartment", values="p_perm")
    p_pivot = p_pivot.reindex(index=fam_order, columns=COMPARTMENTS)
    dscore["family_id"] = dscore["family_id"].astype(int)
    dscore_p_pivot = dscore.pivot(index="family_id", columns="compartment", values="p_perm")
    dscore_p_pivot = dscore_p_pivot.reindex(index=fam_order, columns=COMPARTMENTS)

    n_rows = len(fam_order)
    n_cols = len(COMPARTMENTS)

    # Diverging cmap distinct from Fig 3b's blue-red: purple <-> orange
    cmap = LinearSegmentedColormap.from_list(
        "VHS_diverge",
        [(0.0, "#4B2E83"), (0.45, "#FFFFFF"), (0.55, "#FFFFFF"), (1.0, "#D2691E")],
    )
    cmap.set_under("#1A1545")
    cmap.set_over("#5C2D04")
    cmap.set_bad("#E0E0E0")
    norm = TwoSlopeNorm(vmin=VMIN, vcenter=0.0, vmax=VMAX)

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 150 * MM))
    gs = fig.add_gridspec(
        1, 6,
        width_ratios=[0.30, 0.42, 3.6, 6.4, 0.50, 1.30],
        left=0.04, right=0.985, top=0.91, bottom=0.18,
        wspace=0.04,
    )
    ax_carrier = fig.add_subplot(gs[0, 0])
    ax_dir = fig.add_subplot(gs[0, 1])
    ax_lbl = fig.add_subplot(gs[0, 2])
    ax_main = fig.add_subplot(gs[0, 3])
    ax_sigtally = fig.add_subplot(gs[0, 4])
    ax_leg = fig.add_subplot(gs[0, 5])

    # ---- carrier strip
    for i, fid in enumerate(fam_order):
        c = carriers.get(fid, MIXED_LABEL)
        if c == MIXED_LABEL:
            ax_carrier.add_patch(mpatches.Rectangle(
                (0, i), 1, 1, facecolor=MIXED_COLOR, edgecolor="white", linewidth=0.4, hatch="///",
            ))
        else:
            ax_carrier.add_patch(mpatches.Rectangle(
                (0, i), 1, 1, facecolor=COMPARTMENT_PALETTE.get(c, MIXED_COLOR), edgecolor="white", linewidth=0.4,
            ))
    ax_carrier.set_xlim(0, 1); ax_carrier.set_ylim(n_rows, 0)
    ax_carrier.set_xticks([0.5]); ax_carrier.set_xticklabels(["LA\nbaseline\ncarrier"], fontsize=5.4, color="#444444")
    ax_carrier.set_yticks([])
    for s in ax_carrier.spines.values():
        s.set_visible(False)

    # ---- direction strip
    for i, fid in enumerate(fam_order):
        d = dir_lookup.get(fid, "?")
        ax_dir.add_patch(mpatches.Rectangle((0, i), 1, 1,
                                              facecolor=LA_DEEP if d == "LA-biased" else RA_DEEP,
                                              edgecolor="white", linewidth=0.4))
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(n_rows, 0)
    ax_dir.set_xticks([]); ax_dir.set_yticks([])
    for s in ax_dir.spines.values():
        s.set_visible(False)

    # ---- family labels
    ax_lbl.set_xlim(0, 1); ax_lbl.set_ylim(n_rows, 0); ax_lbl.axis("off")
    for i, fid in enumerate(fam_order):
        d = dir_lookup.get(fid, "?")
        col = LA_DEEP if d == "LA-biased" else RA_DEEP
        ax_lbl.text(0.99, i + 0.5,
                    f"F{int(fid)}  {short_lookup.get(fid, '')}",
                    fontsize=6.5, color=col, fontweight="bold",
                    ha="right", va="center")

    # ---- main heatmap
    mat = vhs_pivot.to_numpy(dtype=float)
    pmat = p_pivot.to_numpy(dtype=float)
    dscore_p = dscore_p_pivot.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(mat)
    im = ax_main.imshow(masked, cmap=cmap, norm=norm, aspect="auto",
                          extent=[-0.5, n_cols - 0.5, n_rows - 0.5, -0.5],
                          interpolation="nearest")
    for i in range(n_rows):
        for j in range(n_cols):
            p = pmat[i, j]
            v = mat[i, j]
            if np.isfinite(p):
                if p < 0.05:
                    ax_main.add_patch(mpatches.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1,
                        facecolor="none", edgecolor=SIG_COLOR, linewidth=1.6, zorder=4,
                    ))
                    ax_main.text(j, i, "*", ha="center", va="center", fontsize=12,
                                  color=SIG_COLOR, fontweight="bold", zorder=5,
                                  path_effects=[mpe.Stroke(linewidth=1.6, foreground="#1A1A1A"),
                                                  mpe.Normal()])
                elif p < 0.10:
                    ax_main.text(j, i, ".", ha="center", va="center", fontsize=11,
                                  color="#1A1A1A", fontweight="bold", zorder=5)
            # marker for cells where dScore is also significant (cross-info)
            if np.isfinite(dscore_p[i, j]) and dscore_p[i, j] < 0.05:
                ax_main.scatter([j + 0.32], [i - 0.32], s=20,
                                  marker="^", facecolor="#1A1A1A", edgecolor="white",
                                  linewidth=0.4, zorder=6)

    ax_main.set_xlim(-0.5, n_cols - 0.5); ax_main.set_ylim(n_rows - 0.5, -0.5)
    ax_main.set_xticks(np.arange(n_cols))
    ax_main.set_xticklabels([COMP_SHORT[c] for c in COMPARTMENTS], rotation=30, ha="right", fontsize=7.0)
    ax_main.set_yticks([])
    ax_main.set_xlabel("cell compartment (LA snRNA)", fontsize=7.5)
    for s in ax_main.spines.values():
        s.set_visible(False)

    # ---- sig tally per family
    sig_per_row = (pmat < 0.05).sum(axis=1)
    ax_sigtally.set_xlim(0, max(sig_per_row.max() + 0.5, 2.5))
    ax_sigtally.set_ylim(n_rows - 0.5, -0.5)
    for i, n_sig in enumerate(sig_per_row):
        if n_sig > 0:
            ax_sigtally.barh(i, n_sig, height=0.55, color="#1A1A1A", edgecolor="none")
            ax_sigtally.text(n_sig + 0.10, i, str(int(n_sig)),
                              fontsize=6.0, ha="left", va="center", color="#1A1A1A", fontweight="bold")
    ax_sigtally.set_xticks([0, 1, 2])
    ax_sigtally.set_xticklabels(["0", "1", "2"], fontsize=5.5, color="#444444")
    ax_sigtally.set_yticks([])
    ax_sigtally.set_xlabel("# sig\nVHS/family", fontsize=5.5, color="#444444")
    for s in ("top", "right", "left"):
        ax_sigtally.spines[s].set_visible(False)
    ax_sigtally.spines["bottom"].set_linewidth(0.4)

    # ---- legend
    ax_leg.set_xlim(0, 1); ax_leg.set_ylim(1, 0); ax_leg.axis("off")
    cax = fig.add_axes([0.91, 0.74, 0.018, 0.13])
    cbar = fig.colorbar(im, cax=cax, orientation="vertical", extend="both")
    cbar.set_label("log\u2082 VHS\n(SD_AF / SD_Ctrl)", fontsize=5.8)
    cbar.ax.tick_params(labelsize=5.6)
    cbar.set_ticks([VMIN, -0.5, 0, 0.5, VMAX])

    ax_leg.text(0.0, 0.46, "Marker key", fontsize=6.6, fontweight="bold", color="#1A1A1A")
    ax_leg.text(0.0, 0.52, "*  VHS p < 0.05", fontsize=5.8, color="#1A1A1A")
    ax_leg.text(0.0, 0.57, ".  VHS p < 0.10", fontsize=5.8, color="#444444")
    ax_leg.text(0.0, 0.63, "\u25B2  dScore also p<0.05\n     (Fig 3b)", fontsize=5.8, color="#1A1A1A")

    ax_leg.text(0.0, 0.74, "Bias", fontsize=6.6, fontweight="bold")
    ax_leg.add_patch(mpatches.Rectangle((0.0, 0.78), 0.10, 0.026, facecolor=LA_DEEP, edgecolor="none"))
    ax_leg.text(0.14, 0.793, "LA-biased", fontsize=5.8)
    ax_leg.add_patch(mpatches.Rectangle((0.0, 0.82), 0.10, 0.026, facecolor=RA_DEEP, edgecolor="none"))
    ax_leg.text(0.14, 0.833, "RA-biased", fontsize=5.8)

    ax_leg.text(0.0, 0.92, "Carrier (Fig 1F)", fontsize=6.0, fontweight="bold")
    ly = 0.96
    present = list(dict.fromkeys([carriers[fid] for fid in fam_order]))
    for c in present:
        if c == MIXED_LABEL:
            ax_leg.add_patch(mpatches.Rectangle((0.0, ly), 0.08, 0.024,
                                                facecolor=MIXED_COLOR, edgecolor="#888", linewidth=0.4, hatch="///"))
            ax_leg.text(0.12, ly + 0.012, "Mixed", fontsize=5.5, va="center")
        else:
            ax_leg.add_patch(mpatches.Rectangle((0.0, ly), 0.08, 0.024,
                                                facecolor=COMPARTMENT_PALETTE[c], edgecolor="none"))
            ax_leg.text(0.12, ly + 0.012, COMP_SHORT[c], fontsize=5.5, va="center")
        ly += 0.028

    fig.suptitle(
        "Fig S3 (VHS)  Variance-Heterogeneity Score: hidden sub-population signal beyond mean \u0394singscore",
        x=0.04, y=0.965, ha="left", fontsize=8.6, fontweight="bold",
    )
    fig.text(
        0.04, 0.020,
        "VHS = SD_AF(family singscore) / SD_Ctrl(family singscore) per (family x compartment) on GSE255612 LA snRNA pseudobulk.",
        fontsize=5.7, color="#444444",
    )
    fig.text(
        0.04, 0.045,
        "Orange (log\u2082 VHS > 0): AF cohort more dispersed than Ctrl \u2192 candidate sub-population emergence in AF (hidden by mean dScore).",
        fontsize=5.7, color="#444444",
    )
    fig.text(
        0.04, 0.070,
        "Purple (log\u2082 VHS < 0): AF cohort more uniform than Ctrl \u2192 cell-intrinsic homogenisation or composition lock-in.",
        fontsize=5.7, color="#444444",
    )
    fig.text(
        0.04, 0.095,
        "Triangle (\u25B2): same (family, ct) cell also shows significant mean \u0394singscore in Fig 3b.  Independence of VHS and dScore signals supports both as orthogonal evidence axes.",
        fontsize=5.7, color="#444444",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
