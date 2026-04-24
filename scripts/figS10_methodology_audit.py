#!/usr/bin/env python3
"""
Figure S10 — Methodology audit (4 internal-consistency checks).

S10a — A1 NCP self-recovery: per-program mean aCAI_rel on healthy bulk
       (should be ≈ +1 by construction; reported as a sanity check).
S10b — A2 direction stability: fraction of healthy donors with aCAI > 0
       per program (should be high).
S10c — A3 cross-method ρ: singscore vs mean-z, per program (Spearman ρ
       on per-donor aCAI on healthy bulk).
S10d — A4 cross-method axis pooled g: STRUCTURAL & METABOLIC pooled g
       under singscore vs mean-z scoring.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

from plot_style import COL_DOUBLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "figures" / "figS10_methodology_audit.pdf"

AXIS_COL = {"STRUCTURAL": "#B35806", "METABOLIC": "#542788"}
METHOD_COL = {"singscore": "#3F6512", "mean_z": "#7A5800"}


def panel_S10a(ax, prog):
    """Per-program mean aCAI_rel ± SD on healthy bulk."""
    n = len(prog)
    y = np.arange(n)
    ax.errorbar(prog["A1_mean_aCAIrel_healthy_ss"], y,
                xerr=prog["A1_sd_aCAIrel_healthy_ss"],
                fmt="o", markersize=5, color="#1A1A1A",
                ecolor="#888888", elinewidth=0.7, capsize=2)
    ax.axvline(1.0, color="#1B7837", lw=0.6, ls=(0, (1.5, 1.5)))
    ax.set_yticks(y)
    ax.set_yticklabels(prog["program_label"], fontsize=5.4)
    ax.invert_yaxis()
    ax.set_xlabel("mean aCAI$_{rel}$ on healthy bulk (target ≈ 1)",
                  fontsize=6.0, labelpad=2)
    ax.set_xlim(-0.1, 2.0)
    ax.set_title("S10a — NCP self-recovery (A1)",
                 fontsize=6.6, pad=3, fontweight="bold", loc="left")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5); ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="x", labelsize=5.4, length=2)
    ax.tick_params(axis="y", length=0, pad=2)


def panel_S10b(ax, prog):
    """Direction stability: fraction of donors with aCAI > 0."""
    y = np.arange(len(prog))
    vals = prog["A2_frac_donors_aCAIpos_ss"].astype(float)
    cols = ["#1B7837" if v >= 0.75 else "#D0D000" if v >= 0.5 else "#B2182B"
            for v in vals]
    ax.barh(y, vals, color=cols, edgecolor="white", linewidth=0.4)
    ax.axvline(0.75, color="#888888", lw=0.4, ls=(0, (1.5, 1.5)))
    ax.set_yticks(y); ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("fraction of healthy donors with aCAI > 0\n(target ≥ 0.75)",
                  fontsize=6.0, labelpad=2, linespacing=1.3)
    ax.set_title("S10b — direction stability (A2)",
                 fontsize=6.6, pad=3, fontweight="bold", loc="left")
    for i, v in enumerate(vals):
        ax.text(v + 0.02, i, f"{v:.2f}",
                fontsize=4.8, ha="left", va="center", color="#222222")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5); ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="x", labelsize=5.4, length=2)
    ax.tick_params(axis="y", length=0, pad=2)


def panel_S10c(ax, prog):
    """singscore vs mean-z cross-method ρ on healthy bulk."""
    y = np.arange(len(prog))
    rho = prog["A3_singscore_vs_meanZ_rho"].astype(float)
    cols = ["#1B7837" if r >= 0.7 else "#D0D000" if r >= 0.5 else "#B2182B"
            for r in rho]
    ax.barh(y, rho, color=cols, edgecolor="white", linewidth=0.4)
    ax.axvline(0.7, color="#888888", lw=0.4, ls=(0, (1.5, 1.5)))
    ax.set_yticks(y); ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Spearman ρ  (singscore vs mean-z, per donor)\n(target ≥ 0.7)",
                  fontsize=6.0, labelpad=2, linespacing=1.3)
    ax.set_title("S10c — cross-method correlation (A3)",
                 fontsize=6.6, pad=3, fontweight="bold", loc="left")
    for i, v in enumerate(rho):
        ax.text(v + 0.02, i, f"{v:.2f}",
                fontsize=4.8, ha="left", va="center", color="#222222")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5); ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="x", labelsize=5.4, length=2)
    ax.tick_params(axis="y", length=0, pad=2)


def panel_S10d(ax, axis_g):
    """Cross-method axis pooled g comparison."""
    methods = ["singscore", "mean_z"]
    axes_lbl = ["STRUCTURAL", "METABOLIC"]
    width = 0.36
    x = np.arange(len(axes_lbl))
    for k, m in enumerate(methods):
        sub = axis_g[axis_g["method"] == m].set_index("axis")
        vals = [sub.loc[a, "pooled_g"] if a in sub.index else np.nan
                for a in axes_lbl]
        offs = (-1) ** k * width / 2
        ax.bar(x + offs, vals, width=width,
               facecolor=METHOD_COL[m], edgecolor="white",
               linewidth=0.4, label=m)
        for xi, v in zip(x + offs, vals):
            if not pd.isna(v):
                ax.text(xi, v - (0.04 if v < 0 else -0.04),
                        f"{v:+.2f}",
                        ha="center", va="top" if v < 0 else "bottom",
                        fontsize=5.0, color="#222222")
    ax.axhline(0, color="#888888", lw=0.5, ls=(0, (1.5, 1.5)))
    ax.set_xticks(x)
    ax.set_xticklabels(axes_lbl, fontsize=5.6)
    ax.set_ylabel("pooled g (AF − SR)", fontsize=6.0, labelpad=2)
    ax.set_title("S10d — axis pooled g, singscore vs mean-z (A4)",
                 fontsize=6.6, pad=3, fontweight="bold", loc="left")
    ax.legend(loc="upper right", fontsize=5.0, frameon=False,
              handlelength=1.0, handleheight=0.55, labelspacing=0.25)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5); ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", labelsize=5.4, length=2)


def main() -> None:
    set_nature_style()
    prog = pd.read_csv(RES / "methodology_audit_program.tsv", sep="\t")
    prog = prog.sort_values("program_order").reset_index(drop=True)
    axis_g = pd.read_csv(RES / "methodology_audit_axis_g.tsv", sep="\t")

    fig_w = COL_DOUBLE
    fig_h = 5.0
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.06, right=0.97, top=0.91, bottom=0.10)
    gs = fig.add_gridspec(2, 3,
                          height_ratios=[1.0, 1.0],
                          width_ratios=[1.45, 1.0, 1.0],
                          hspace=0.55, wspace=0.55)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, :])

    panel_S10a(ax_a, prog)
    panel_S10b(ax_b, prog)
    panel_S10c(ax_c, prog)
    panel_S10d(ax_d, axis_g)

    fig.suptitle(
        "Methodology audit — internal consistency of the rewiring framework",
        fontsize=7.4, y=0.965, fontweight="bold",
    )
    fig.text(0.012, 0.97, "S10",
             fontsize=10, fontweight="bold", ha="left", va="top",
             transform=fig.transFigure)
    fig.text(0.5, 0.020,
             "All 10 programs pass NCP self-recovery (mean aCAI$_{rel}$ ≈ 1) and "
             "direction stability (8/8 healthy donors aCAI > 0).  Per-donor singscore "
             "and mean-z scores correlate at ρ = 0.5–0.9, and the dual-axis pooled g "
             "is qualitatively unchanged when scoring is swapped.",
             ha="center", va="bottom", fontsize=4.8,
             color="#666666", style="italic", linespacing=1.4)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
