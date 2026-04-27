#!/usr/bin/env python3
"""
Figure S5 — Δsingscore (AF − SR) heatmap on GSE238242 LAA snRNA pseudobulks.

6 frozen support programs × 4 main compartments (atrialCM, fibroblast,
endothelial_or_vascular, immune).  Honest negative finding: the
pre-registered support cells (locked program × compartment grid) are
**marked**, and none reaches "moderate / strong" direction support.
This panel is therefore a **Supplementary** sanity table, NOT a main
mechanism figure.

Inputs:
  Project_P5/GSE238242/results/figure5_program_delta_main.tsv
  Project_P5/GSE238242/results/figure5_program_support_check.tsv

Output:
  figures/figS5_multiome_delta_heatmap.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

from plot_style import AFCTRL_CMAP, COL_ONEHALF, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "Project_P5" / "GSE238242" / "results"
OUT = ROOT / "figures" / "figS5_multiome_delta_heatmap.pdf"

# Locked program order (axis groups): STRUCTURAL on top, METABOLIC bottom
PROGRAM_ORDER = [
    ("ECM proteoglycans",       "STRUCTURAL"),
    ("HS-GAG Metabolism",       "STRUCTURAL"),
    ("HS-GAG Biosynthesis",     "STRUCTURAL"),
    ("OXPHOS",                  "METABOLIC"),
    ("Aerobic Resp. & ETC",     "METABOLIC"),
    ("Cytoplasmic Translation", "METABOLIC"),
]
COMPARTMENT_ORDER = ["atrialCM", "fibroblast",
                     "endothelial_or_vascular", "immune"]
COMPARTMENT_SHORT = {
    "atrialCM": "atrialCM",
    "fibroblast": "fibroblast",
    "endothelial_or_vascular": "endo/vasc",
    "immune": "immune",
}
# Locked expected cells (axis -> compartment, sign)
EXPECT = {
    "STRUCTURAL": ("fibroblast", +1),
    "METABOLIC":  ("atrialCM",   -1),
}


def main() -> None:
    set_nature_style()

    delta = pd.read_csv(IN / "figure5_program_delta_main.tsv", sep="\t")
    sup = pd.read_csv(IN / "figure5_program_support_check.tsv", sep="\t")

    # Build (program × compartment) Δ matrix in locked order
    n_prog = len(PROGRAM_ORDER)
    n_cmp = len(COMPARTMENT_ORDER)
    M = np.full((n_prog, n_cmp), np.nan)
    for i, (lab, _ax) in enumerate(PROGRAM_ORDER):
        for j, c in enumerate(COMPARTMENT_ORDER):
            row = delta[(delta["program_label"] == lab) &
                        (delta["compartment"] == c)]
            if len(row) == 1:
                M[i, j] = float(row["delta_score"].iloc[0])

    # Symmetric color scale
    vmax = max(0.05, float(np.nanmax(np.abs(M))) * 1.05)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list("af_ctrl", AFCTRL_CMAP, N=256)

    fig_w = COL_ONEHALF
    fig_h = 0.36 * n_prog + 1.95
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.34, right=0.78, top=0.78, bottom=0.22)

    # Cells
    for i in range(n_prog):
        for j in range(n_cmp):
            d = M[i, j]
            if np.isnan(d):
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                       facecolor="#F2F2F2",
                                       edgecolor="white", linewidth=0.5))
                continue
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                   facecolor=cmap(norm(d)),
                                   edgecolor="white", linewidth=0.5))
            txt_color = "white" if abs(norm(d) - 0.5) > 0.30 else "#222222"
            ax.text(j, i, f"{d:+.3f}",
                    ha="center", va="center", fontsize=5.6,
                    color=txt_color, fontweight="bold")

    # Mark locked-expected cells with thick black border (the support gate cells)
    for ax_name, (cmp_exp, _sgn) in EXPECT.items():
        j = COMPARTMENT_ORDER.index(cmp_exp)
        for i, (lab, ax_lab) in enumerate(PROGRAM_ORDER):
            if ax_lab == ax_name:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                       facecolor="none",
                                       edgecolor="#222222", linewidth=1.2,
                                       zorder=4))

    # Per-cell support_level chip (only for the locked gate cells)
    sup_by_pid = sup.set_index("program_label")
    LEVEL_DOT = {"strong": "●", "moderate": "◐", "weak": "·", "absent": "x"}
    LEVEL_COL = {"strong": "#1B7837", "moderate": "#5AAE61",
                 "weak": "#999999", "absent": "#B2182B"}
    for i, (lab, ax_lab) in enumerate(PROGRAM_ORDER):
        if lab in sup_by_pid.index:
            lvl = sup_by_pid.at[lab, "support_level"]
            cmp_exp = sup_by_pid.at[lab, "compartment_expected"]
            j = COMPARTMENT_ORDER.index(cmp_exp)
            ax.text(j + 0.42, i - 0.42, LEVEL_DOT.get(lvl, "?"),
                    ha="right", va="top", fontsize=7,
                    color=LEVEL_COL.get(lvl, "#444444"), fontweight="bold")

    # Axes
    ax.set_xlim(-0.55, n_cmp - 0.45)
    ax.set_ylim(n_prog - 0.5, -0.5)
    ax.set_xticks(range(n_cmp))
    ax.set_xticklabels([COMPARTMENT_SHORT[c] for c in COMPARTMENT_ORDER],
                       fontsize=6, rotation=30, ha="right")
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_yticks(range(n_prog))
    ylabels = [
        f"{lab}\n{('STRUCT.' if ax_lab=='STRUCTURAL' else 'METAB.')}"
        for lab, ax_lab in PROGRAM_ORDER
    ]
    ax.set_yticklabels(ylabels, fontsize=6, linespacing=1.15)
    ax.tick_params(axis="y", length=0, pad=2)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    # Colorbar
    cbar_ax = fig.add_axes([0.36, 0.085, 0.36, 0.024])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cbar_ax, orientation="horizontal")
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(width=0.4, length=2, labelsize=5.5, pad=1)
    cbar.set_ticks([-vmax, -vmax / 2, 0, vmax / 2, vmax])
    cbar.set_ticklabels([f"{-vmax:+.02f}", f"{-vmax/2:+.02f}", "0",
                         f"{vmax/2:+.02f}", f"{vmax:+.02f}"])
    cbar.ax.set_title(r"$\Delta$ singscore (AF − SR), donor-level",
                      fontsize=6, pad=4)
    fig.text(0.36, 0.142, "AF < SR", fontsize=5.5,
             ha="left", va="bottom", color="#01665E")
    fig.text(0.72, 0.142, "AF > SR", fontsize=5.5,
             ha="right", va="bottom", color="#8C510A")

    # Right-side legend: support gate symbols
    lx, ly0 = 0.795, 0.74
    fig.text(lx, ly0, "support gate",
             fontsize=6, fontweight="bold", color="#222222",
             ha="left", va="bottom")
    fig.text(lx, ly0 - 0.030,
             "thick black box =\nlocked expected cell",
             fontsize=5.2, color="#333333",
             ha="left", va="top", linespacing=1.3)
    fig.text(lx, ly0 - 0.10,
             "support level (in expected cell):",
             fontsize=5.4, color="#333333",
             ha="left", va="top")
    for k, (lvl, mark) in enumerate(LEVEL_DOT.items()):
        y = ly0 - 0.13 - k * 0.028
        fig.text(lx, y, mark, fontsize=7, color=LEVEL_COL[lvl],
                 ha="left", va="center", fontweight="bold")
        fig.text(lx + 0.03, y, lvl, fontsize=5.4, color="#333333",
                 ha="left", va="center")

    # Title + caveat
    fig.suptitle(
        "GSE238242 LAA snRNA pseudobulk — Δsingscore on 6 frozen support "
        "programs × 4 main compartments",
        fontsize=7, y=0.96, fontweight="bold",
    )
    fig.text(0.5, 0.025,
             "AF n=3, SR n=4 (LAA, valvular surgery cohort).  "
             "By the locked support-gate rule (this script), no expected "
             "cell reaches moderate/strong support;\nthe multiome regulatory "
             "layer is therefore not pursued in the main figure.",
             ha="center", va="bottom", fontsize=5.0, color="#444444",
             linespacing=1.4)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
