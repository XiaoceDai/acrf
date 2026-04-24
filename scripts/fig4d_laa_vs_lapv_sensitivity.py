#!/usr/bin/env python3
"""
Figure 4D - LAA (primary) vs LAPV (sensitivity) ρ side-by-side.

Two-column mini-heatmap, same row order as Fig 4B.
Narrative (user-locked): "appendage-weighted / LAA-dominant trends with
region-specific heterogeneity" - NOT "appendage-specific only" because
several programs carry a same-sign but weaker LAPV trend.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle, Circle

from plot_style import COL_ONEHALF, DURATION_CMAP, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "figures" / "fig4d_laa_vs_lapv_sensitivity.pdf"
IN   = ROOT / "Project_P4" / "GSE41177" / "results"


def main() -> None:
    set_nature_style()
    laa  = pd.read_csv(IN / "duration_trend_laa_primary.tsv", sep="\t")
    lapv = pd.read_csv(IN / "duration_trend_lapv_sensitivity.tsv", sep="\t")
    laa  = laa.sort_values("program_order").reset_index(drop=True)

    lapv_by_id = lapv.set_index("program_id")
    n_row = len(laa)
    n_col = 2

    vmax = max(0.6, float(max(laa["statistic"].abs().max(),
                              lapv["statistic"].abs().max())) * 1.1)
    vmin = -vmax
    norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list("dur_prgn", DURATION_CMAP, N=256)

    fig_w = COL_ONEHALF
    fig_h = 0.30 * n_row + 1.88
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.38, right=0.78, top=0.86, bottom=0.20)

    def _draw_cell(j, i, rho, marker):
        rgba = (1, 1, 1, 0) if pd.isna(rho) else cmap(norm(rho))
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0,
                               facecolor=rgba, edgecolor="white", linewidth=0.5))
        if not pd.isna(rho):
            txt_color = "white" if abs(norm(rho) - 0.5) > 0.30 else "#222222"
            ax.text(j, i - 0.08, f"{rho:+.2f}",
                    ha="center", va="center", fontsize=6.2,
                    color=txt_color, fontweight="bold")
            if isinstance(marker, str) and marker:
                ax.text(j, i + 0.24, marker,
                        ha="center", va="center", fontsize=9,
                        fontweight="bold", color=txt_color)

    # Fill cells + compute sign-concordance per row
    concord_flags = []
    for i, row in laa.iterrows():
        pid = row["program_id"]
        laa_rho = row["statistic"]
        laa_mk  = row.get("sig_marker", "")
        if pid in lapv_by_id.index:
            lapv_rho = lapv_by_id.at[pid, "statistic"]
            lapv_mk  = lapv_by_id.at[pid, "sig_marker"]
        else:
            lapv_rho = float("nan"); lapv_mk = ""

        _draw_cell(0, i, laa_rho,  laa_mk)
        _draw_cell(1, i, lapv_rho, lapv_mk)

        # concordance (|ρ|>0.02 both sides and same sign)
        if any(pd.isna(v) for v in [laa_rho, lapv_rho]):
            concord_flags.append(False)
        elif abs(laa_rho) < 0.02 or abs(lapv_rho) < 0.02:
            concord_flags.append(False)
        else:
            concord_flags.append(np.sign(laa_rho) == np.sign(lapv_rho))

    # sign-concordance chip on the right
    CHIP_X = n_col - 0.5 + 0.60
    for i, ok in enumerate(concord_flags):
        if ok:
            ax.add_patch(Circle((CHIP_X, i), 0.13, facecolor="#444444",
                                edgecolor="#444444", linewidth=0.5))
        else:
            ax.add_patch(Circle((CHIP_X, i), 0.13, facecolor="white",
                                edgecolor="#888888", linewidth=0.5))

    ax.set_xlim(-0.6, CHIP_X + 0.4)
    ax.set_ylim(n_row - 0.5, -0.5)
    ax.set_aspect("auto")

    ax.set_xticks([0, 1, CHIP_X])
    ax.set_xticklabels([f"LAA ρ\n(primary)", f"LAPV ρ\n(sensitivity)",
                        "sign\nconcordant"],
                       fontsize=6, linespacing=1.15)
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_yticks(range(n_row))
    ax.set_yticklabels(laa["program_label"].tolist())
    ax.tick_params(axis="y", length=0, pad=2)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)

    # colorbar
    cbar_ax = fig.add_axes([0.40, 0.095, 0.36, 0.024])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cbar_ax, orientation="horizontal")
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(width=0.4, length=2, labelsize=5.5, pad=1)
    cbar.set_ticks([-0.6, -0.3, 0.0, 0.3, 0.6])
    cbar.ax.set_title(r"Spearman $\rho$  (singscore vs AF duration, months)",
                      fontsize=6, pad=4)
    fig.text(0.40, 0.152, "falls w. duration", fontsize=5.5,
             ha="left", va="bottom", color="#762A83")
    fig.text(0.76, 0.152, "rises w. duration", fontsize=5.5,
             ha="right", va="bottom", color="#1B7837")

    # Right-side legend (aligned with Fig 4B)
    lx, y0 = 0.795, 0.785
    fig.text(lx, y0, "perm p",
             fontsize=6, fontweight="bold", ha="left", va="bottom", color="#333333")
    fig.text(lx, y0 - 0.028, r"$\star$ < 0.05",
             fontsize=5.5, ha="left", va="top", color="#333333")
    fig.text(lx, y0 - 0.052, ".  < 0.10",
             fontsize=5.5, ha="left", va="top", color="#333333")

    y1 = 0.615
    fig.text(lx, y1, "sign",
             fontsize=6, fontweight="bold", ha="left", va="bottom", color="#333333")
    fig.patches.append(Circle((0.803, 0.588), 0.009,
                              transform=fig.transFigure,
                              facecolor="#444444", edgecolor="#444444"))
    fig.text(0.817, 0.588, "same sign",
             fontsize=5.5, ha="left", va="center", color="#333333")
    fig.patches.append(Circle((0.803, 0.562), 0.009,
                              transform=fig.transFigure,
                              facecolor="white", edgecolor="#888888"))
    fig.text(0.817, 0.562, "discordant / ~0",
             fontsize=5.5, ha="left", va="center", color="#333333")

    fig.text(lx, 0.485, "region",
             fontsize=6, fontweight="bold", ha="left", va="bottom", color="#333333")
    fig.text(lx, 0.448,
             "LAA primary\nLAPV sensitivity\n(BH families\nseparate)",
             fontsize=5.3, ha="left", va="top", color="#555555", linespacing=1.35)

    fig.suptitle(
        "LAA-weighted duration trends with region-specific heterogeneity",
        fontsize=7.4, y=0.965, fontweight="bold",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
