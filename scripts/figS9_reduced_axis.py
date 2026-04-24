#!/usr/bin/env python3
"""
Figure S9 — Reduced-axis readout (feasibility check, not a classifier).

3 sub-panels:
  S9a — bar of Spearman ρ (reduced vs full axis composite) per axis × cohort
  S9b — donor-level scatter (full vs reduced) per axis × cohort
  S9c — LAA duration ρ comparison (reduced matrix / metabolic vs full OXPHOS)

Inputs (all from h6_reduced_axis.R):
  results/figure_s9_reduced_panels.tsv
  results/figure_s9_reduced_scores_long.tsv
  results/figure_s9_reduced_correlation.tsv
  results/figure_s9_reduced_centroids.tsv
  results/figure_s9_reduced_duration.tsv
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle

from plot_style import COL_DOUBLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "figures" / "figS9_reduced_axis.pdf"

AXIS_COL = {"matrix": "#B35806", "metabolic": "#542788"}
COHORT_MARKER = {"GSE128188": "o", "GSE79768": "s", "GSE115574": "^"}
COL_AF = "#B2182B"
COL_SR = "#2166AC"


def _round(ax, x, y, w, h, fc, ec="none", lw=0.4, z=1, rounding=0.012):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    ))


def panel_S9a(ax, cor_tbl):
    """Grouped bars: Spearman ρ (reduced vs full) per axis × cohort."""
    cohorts = ["GSE128188", "GSE79768", "GSE115574"]
    axes_lbl = ["matrix", "metabolic"]

    width = 0.36
    x = np.arange(len(cohorts))
    for k, ax_lbl in enumerate(axes_lbl):
        sub = cor_tbl[cor_tbl["axis"] == ax_lbl].set_index("cohort")
        vals = [sub.loc[c, "spearman_rho"] if c in sub.index else np.nan
                for c in cohorts]
        offs = (-1) ** k * width / 2
        ax.bar(x + offs, vals, width=width,
               facecolor=AXIS_COL[ax_lbl], edgecolor="white",
               linewidth=0.4, label=ax_lbl)
        for xi, v in zip(x + offs, vals):
            if not np.isnan(v):
                ax.text(xi, v + 0.02, f"{v:.2f}",
                        ha="center", va="bottom", fontsize=4.6,
                        color="#222222")
    ax.axhline(0, color="#888888", lw=0.5, ls=(0, (1.5, 1.5)), zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(cohorts, fontsize=5.6)
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("Spearman ρ\n(reduced vs full)", fontsize=6.0,
                  labelpad=2)
    ax.set_title("S9a — reduced panel reproduces the full axis composite",
                 fontsize=6.8, pad=3, fontweight="bold", loc="left")
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", labelsize=5.6, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.legend(loc="upper right", fontsize=5.0, frameon=False,
              handlelength=1.0, handleheight=0.55, labelspacing=0.25)


def panel_S9b(ax, scores):
    """Per-axis × cohort scatter: full (x) vs reduced (y) aCAI_rel."""
    cohorts = ["GSE128188", "GSE79768", "GSE115574"]
    for ax_lbl in ("matrix", "metabolic"):
        sub_axis = scores[scores["axis"] == ax_lbl]
        for cohort in cohorts:
            sub = sub_axis[sub_axis["cohort"] == cohort]
            color = AXIS_COL[ax_lbl]
            mark = COHORT_MARKER.get(cohort, "o")
            for _, r in sub.iterrows():
                fc = COL_AF if r["rhythm"] == "AF" else COL_SR
                ax.scatter(r["aCAI_rel_full"], r["aCAI_rel_reduced"],
                           marker=mark, s=20, facecolor=fc,
                           edgecolor=color, linewidth=0.5, alpha=0.85,
                           zorder=4)
    ax.axhline(0, color="#BBBBBB", lw=0.4, zorder=1)
    ax.axvline(0, color="#BBBBBB", lw=0.4, zorder=1)
    # diagonal
    lim = max(abs(scores["aCAI_rel_full"]).max(),
              abs(scores["aCAI_rel_reduced"]).max()) * 1.05
    ax.plot([-lim, lim], [-lim, lim], color="#888888", lw=0.5, ls=":",
            zorder=2)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("full axis composite (aCAI$_{rel}$)",
                  fontsize=6.0, labelpad=2)
    ax.set_ylabel("reduced 5-gene score (aCAI$_{rel}$)",
                  fontsize=6.0, labelpad=2)
    ax.set_title("S9b — full vs reduced, donor-level",
                 fontsize=6.8, pad=3, fontweight="bold", loc="left")
    ax.tick_params(axis="both", labelsize=5.4, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    # legend
    handles = [
        Line2D([], [], marker="o", linestyle="None", markersize=4,
               markerfacecolor=COL_AF, markeredgecolor="#444444",
               markeredgewidth=0.4, label="AF"),
        Line2D([], [], marker="o", linestyle="None", markersize=4,
               markerfacecolor=COL_SR, markeredgecolor="#444444",
               markeredgewidth=0.4, label="SR"),
        Line2D([], [], marker="s", linestyle="None", markersize=3.5,
               markerfacecolor=AXIS_COL["matrix"],
               markeredgecolor="white", markeredgewidth=0.3,
               label="matrix"),
        Line2D([], [], marker="s", linestyle="None", markersize=3.5,
               markerfacecolor=AXIS_COL["metabolic"],
               markeredgecolor="white", markeredgewidth=0.3,
               label="metabolic"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=4.6,
              frameon=False, ncol=2, handlelength=1.0,
              handleheight=0.5, labelspacing=0.2, columnspacing=1.0)


def panel_S9c(ax, dur):
    """LAA duration ρ — reduced matrix / reduced metabolic / full OXPHOS."""
    rows = []
    for _, r in dur.iterrows():
        rows.append({"label": f"reduced {r['axis']}",
                     "rho": r["reduced_rho_LAA"],
                     "color": AXIS_COL[r["axis"]]})
    # full OXPHOS
    full_oxphos = dur[dur["axis"] == "metabolic"]["full_rho_LAA_reference"].dropna()
    if not full_oxphos.empty:
        rows.append({"label": "full OXPHOS\n(Fig 4 reference)",
                     "rho": float(full_oxphos.iloc[0]),
                     "color": "#1A1A1A"})
    df = pd.DataFrame(rows)
    y = np.arange(len(df))
    for i, r in df.iterrows():
        ax.plot([0, r["rho"]], [i, i], color=r["color"], lw=2.5,
                solid_capstyle="round", zorder=3)
        ax.scatter([r["rho"]], [i], s=42, marker="o",
                   facecolor=r["color"], edgecolor="white",
                   linewidth=0.6, zorder=4)
        ax.text(r["rho"] + (0.04 if r["rho"] >= 0 else -0.04), i,
                f"ρ = {r['rho']:+.2f}",
                ha="left" if r["rho"] >= 0 else "right",
                va="center", fontsize=5.0, color="#222222",
                family="monospace")
    ax.axvline(0, color="#888888", lw=0.45, ls=(0, (1.5, 1.5)), zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"].tolist(), fontsize=5.4)
    ax.invert_yaxis()
    ax.set_xlabel("Spearman ρ  (singscore vs AF duration, GSE41177 LAA)",
                  fontsize=6.0, labelpad=2)
    ax.set_xlim(-0.85, 0.55)
    ax.set_title("S9c — reduced metabolic does NOT replace full OXPHOS for chronicity",
                 fontsize=6.4, pad=3, fontweight="bold", loc="left")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="x", labelsize=5.4, length=2)
    ax.tick_params(axis="y", length=0, pad=2)


def main() -> None:
    set_nature_style()
    panels = pd.read_csv(RES / "figure_s9_reduced_panels.tsv", sep="\t")
    cor   = pd.read_csv(RES / "figure_s9_reduced_correlation.tsv", sep="\t")
    sc    = pd.read_csv(RES / "figure_s9_reduced_scores_long.tsv", sep="\t")
    dur   = pd.read_csv(RES / "figure_s9_reduced_duration.tsv", sep="\t")

    fig_w = COL_DOUBLE
    fig_h = 4.0
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.07, right=0.97, top=0.84, bottom=0.16)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0],
                          width_ratios=[1.0, 1.0, 1.0],
                          hspace=0.55, wspace=0.45)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1:])
    ax_c = fig.add_subplot(gs[1, :])

    panel_S9a(ax_a, cor)
    panel_S9b(ax_b, sc)
    panel_S9c(ax_c, dur)

    # Locked panel composition footer
    matrix_genes = panels[panels["axis"] == "matrix"]["gene"].tolist()
    metab_genes  = panels[panels["axis"] == "metabolic"]["gene"].tolist()

    fig.suptitle(
        "Reduced-axis readout — feasibility check, not a classifier",
        fontsize=7.4, y=0.96, fontweight="bold",
    )
    fig.text(0.012, 0.97, "S9",
             fontsize=10, fontweight="bold", ha="left", va="top",
             transform=fig.transFigure)
    fig.text(0.5, 0.045,
             f"matrix reduced panel:  {' · '.join(matrix_genes)}",
             ha="center", va="bottom", fontsize=5.2, color="#444444")
    fig.text(0.5, 0.025,
             f"metabolic reduced panel:  {' · '.join(metab_genes)}",
             ha="center", va="bottom", fontsize=5.2, color="#444444")
    fig.text(0.5, 0.005,
             "Genes drawn only from already-used component sets (Figs 3E + 6C "
             "+ Fig 1 OXPHOS LE).  This is a framework-feasibility check, "
             "not a clinical classifier.",
             ha="center", va="bottom", fontsize=4.6,
             color="#666666", style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
