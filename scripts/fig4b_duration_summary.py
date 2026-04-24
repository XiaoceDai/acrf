#!/usr/bin/env python3
"""
Figure 4B - Panel program × duration-trend ρ heatmap.

Single main column: LAA primary ρ (n = 12).
Markers: `*` if p_perm < 0.05,  `.` if 0.05 ≤ p_perm < 0.10.
Side indicator: sign-concordance between primary and sens-A (filled dot =
concordant, hollow = discordant).
Row order: frozen `program_order` shared with Fig 1C / Fig 2 / Fig 3.

No LA/RA identity strip (same principle as Fig 3B - Fig 4 axis is duration
trend, not chamber identity).
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
OUT  = ROOT / "figures" / "fig4b_duration_summary.pdf"
IN   = ROOT / "Project_P4" / "GSE41177" / "results"


def main() -> None:
    set_nature_style()
    primary = pd.read_csv(IN / "duration_trend_laa_primary.tsv", sep="\t")
    sens_a  = pd.read_csv(IN / "duration_trend_laa_sens_gt1m_eq1.tsv", sep="\t")

    # align by program_order
    primary = primary.sort_values("program_order").reset_index(drop=True)
    sens_a_by_id = sens_a.set_index("program_id")

    # sign concordance flag
    def _concord(row):
        rho_p = row["statistic"]
        rho_s = sens_a_by_id.loc[row["program_id"], "statistic"]
        if pd.isna(rho_p) or pd.isna(rho_s):
            return False
        if abs(rho_p) < 0.02 or abs(rho_s) < 0.02:
            return False   # both near zero; call it ambiguous rather than concordant
        return np.sign(rho_p) == np.sign(rho_s)
    primary["sign_concordant"] = primary.apply(_concord, axis=1)

    n_row = len(primary)

    vmax = max(0.6, float(primary["statistic"].abs().max()) * 1.1)
    vmin = -vmax
    norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list("dur_prgn", DURATION_CMAP, N=256)

    fig_w = COL_ONEHALF
    fig_h = 0.30 * n_row + 1.78
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.38, right=0.78, top=0.86, bottom=0.20)

    CELL_X = 0          # single main column (LAA ρ)
    CONC_X = 1.00       # sign-concordance chip column

    # Main cells (LAA ρ)
    for i, row in primary.iterrows():
        rho = row["statistic"]
        rgba = (1, 1, 1, 0) if pd.isna(rho) else cmap(norm(rho))
        ax.add_patch(Rectangle((CELL_X - 0.5, i - 0.5), 1.0, 1.0,
                               facecolor=rgba, edgecolor="white", linewidth=0.5))
        if not pd.isna(rho):
            txt_color = "white" if abs(norm(rho) - 0.5) > 0.30 else "#222222"
            ax.text(CELL_X, i - 0.08, f"{rho:+.2f}",
                    ha="center", va="center", fontsize=6.2,
                    color=txt_color, fontweight="bold")
            mk = row["sig_marker"]
            if isinstance(mk, str) and mk:
                ax.text(CELL_X, i + 0.24, mk,
                        ha="center", va="center", fontsize=9,
                        fontweight="bold", color=txt_color)

    # Sign-concordance chips (filled = concordant, hollow circle = discordant)
    for i, row in primary.iterrows():
        y = i
        if row["sign_concordant"]:
            ax.add_patch(Circle((CONC_X, y), 0.13, facecolor="#444444",
                                edgecolor="#444444", linewidth=0.5))
        else:
            ax.add_patch(Circle((CONC_X, y), 0.13, facecolor="white",
                                edgecolor="#888888", linewidth=0.5))

    ax.set_xlim(-0.6, CONC_X + 0.4)
    ax.set_ylim(n_row - 0.5, -0.5)
    ax.set_aspect("auto")

    ax.set_xticks([CELL_X, CONC_X])
    ax.set_xticklabels([f"LAA ρ\n(n=12)", "sign\nconcordance"],
                       fontsize=6, linespacing=1.15)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_yticks(range(n_row))
    ax.set_yticklabels(primary["program_label"].tolist())
    ax.tick_params(axis="y", length=0, pad=2)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)

    # colorbar — slightly higher for clear separation from axis labels
    cbar_ax = fig.add_axes([0.40, 0.075, 0.36, 0.024])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cbar_ax, orientation="horizontal")
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(width=0.4, length=2, labelsize=5.5, pad=1)
    cbar.set_ticks([-0.6, -0.3, 0.0, 0.3, 0.6])
    cbar.ax.set_title(r"Spearman $\rho$  (singscore vs AF duration, months)",
                      fontsize=6, pad=4)
    fig.text(0.40, 0.132, "falls w. duration", fontsize=5.5,
             ha="left", va="bottom", color="#762A83")
    fig.text(0.76, 0.132, "rises w. duration", fontsize=5.5,
             ha="right", va="bottom", color="#1B7837")

    # right-side legend — aligned block with consistent line height
    lx, y0 = 0.795, 0.785
    fig.text(lx, y0, "perm p",
             fontsize=6, fontweight="bold", ha="left", va="bottom",
             color="#333333")
    fig.text(lx, y0 - 0.028, r"$\star$ < 0.05",
             fontsize=5.5, ha="left", va="top", color="#333333")
    fig.text(lx, y0 - 0.052, ".  < 0.10",
             fontsize=5.5, ha="left", va="top", color="#333333")

    y1 = 0.615
    fig.text(lx, y1, "sign",
             fontsize=6, fontweight="bold", ha="left", va="bottom",
             color="#333333")
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

    fig.suptitle(
        "LAA program trends with AF duration show sensitivity-supported directionality",
        fontsize=7.4, y=0.965, fontweight="bold",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
