#!/usr/bin/env python3
"""
Figure S6 — Mural (Pericyte + VSMC) compartment sensitivity for the 6
frozen support programs.  Reported as supplementary; not part of the
4-compartment main analysis (FIGURE5_WORKFLOW.md).

Reads:
  Project_P5/GSE238242/results/figure5_program_delta_mural.tsv
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle, FancyBboxPatch

from plot_style import AFCTRL_CMAP, COL_ONEHALF, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "Project_P5" / "GSE238242" / "results" / "figure5_program_delta_mural.tsv"
OUT = ROOT / "figures" / "figS6_mural_sensitivity.pdf"

PROGRAM_ORDER = [
    ("ECM proteoglycans",       "STRUCTURAL"),
    ("HS-GAG Metabolism",       "STRUCTURAL"),
    ("HS-GAG Biosynthesis",     "STRUCTURAL"),
    ("OXPHOS",                  "METABOLIC"),
    ("Aerobic Resp. & ETC",     "METABOLIC"),
    ("Cytoplasmic Translation", "METABOLIC"),
]


def main() -> None:
    set_nature_style()
    df = pd.read_csv(IN, sep="\t")
    df_idx = df.set_index("program_label")

    n_prog = len(PROGRAM_ORDER)
    n_col = 1
    M = np.full((n_prog, n_col), np.nan)
    for i, (lab, _) in enumerate(PROGRAM_ORDER):
        if lab in df_idx.index:
            M[i, 0] = float(df_idx.at[lab, "delta_score"])

    vmax = max(0.10, float(np.nanmax(np.abs(M))) * 1.10)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list("af_ctrl", AFCTRL_CMAP, N=256)

    fig_w = COL_ONEHALF * 0.85
    fig_h = 0.40 * n_prog + 1.45
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.42, right=0.78, top=0.78, bottom=0.20)

    for i in range(n_prog):
        d = M[i, 0]
        if np.isnan(d):
            ax.add_patch(Rectangle((-0.5, i - 0.5), 1, 1,
                                   facecolor="#F2F2F2",
                                   edgecolor="white", linewidth=0.5))
            continue
        ax.add_patch(Rectangle((-0.5, i - 0.5), 1, 1,
                               facecolor=cmap(norm(d)),
                               edgecolor="white", linewidth=0.5))
        txt_color = "white" if abs(norm(d) - 0.5) > 0.30 else "#222222"
        ax.text(0, i, f"{d:+.3f}",
                ha="center", va="center", fontsize=6.0,
                color=txt_color, fontweight="bold")

    ax.set_xlim(-0.55, 0.55)
    ax.set_ylim(n_prog - 0.5, -0.5)
    ax.set_xticks([0])
    ax.set_xticklabels(["mural\n(PC + SMC)"],
                       fontsize=6, rotation=0, ha="center")
    ax.tick_params(axis="x", length=0, pad=4)
    ax.set_yticks(range(n_prog))
    ax.set_yticklabels([
        f"{lab}\n{('STRUCT.' if ax_lab=='STRUCTURAL' else 'METAB.')}"
        for lab, ax_lab in PROGRAM_ORDER
    ], fontsize=6.0, linespacing=1.15)
    ax.tick_params(axis="y", length=0, pad=2)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    cbar_ax = fig.add_axes([0.44, 0.085, 0.34, 0.024])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cbar_ax, orientation="horizontal")
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(width=0.4, length=2, labelsize=5.5, pad=1)
    cbar.set_ticks([-vmax, 0, vmax])
    cbar.set_ticklabels([f"{-vmax:+.02f}", "0", f"{vmax:+.02f}"])
    cbar.ax.set_title(r"$\Delta$ singscore (AF − SR), donor-level",
                      fontsize=6, pad=4)
    fig.text(0.44, 0.142, "AF < SR", fontsize=5.5,
             ha="left", va="bottom", color="#01665E")
    fig.text(0.78, 0.142, "AF > SR", fontsize=5.5,
             ha="right", va="bottom", color="#8C510A")

    fig.suptitle(
        "Mural sensitivity — Δsingscore on the 6 frozen support programs",
        fontsize=7.0, y=0.96, fontweight="bold",
    )
    fig.text(0.5, 0.025,
             "Supplementary: GSE238242 mural pseudobulk (Pericyte + VSMC, "
             "n_AF = 3, n_SR = 4); not part of the 4-compartment main "
             "analysis.  Directionally consistent with main metabolic "
             "panel; no claim of an independent mural arm.",
             ha="center", va="bottom", fontsize=4.7,
             color="#666666", style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
