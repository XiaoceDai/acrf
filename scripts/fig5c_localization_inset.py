#!/usr/bin/env python3
"""
Figure 5C — ROI-level localization of fibroblast-associated ECM signal.

Subordinate spatial inset panel (non-primary):
  - LEFT  : ECM proteoglycans score in AA ROIs, x = area (border / distal),
            colored by rhythm (AF / SR), area effect only.
  - MIDDLE: fibroblast (canonical) vs ECM proteoglycans scatter (24 AA ROIs),
            Spearman ρ printed in panel; rhythm + area encoded.
  - RIGHT : one short callout sentence.

Title constraint: do NOT advertise "spatial" — GeoMx provides ROI-level
segments, not Visium/Xenium-style tissue heatmaps.

Color discipline:
  - AF / SR are encoded as low-saturation tints (not the saturated Fig 3 BrBG)
  - border / distal are encoded as point shapes, not a fifth color axis.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from plot_style import COL_DOUBLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "Project_P5" / "GSE261363" / "roi_scores.tsv"
OUT = ROOT / "figures" / "fig5c_localization_inset.pdf"

# Low-saturation rhythm tints (deliberately muted; do not compete with 5D)
COL_AF = "#9C6B49"   # muted earthy brown
COL_SR = "#65857F"   # muted teal-grey
SHAPE = {"border": "^", "distal": "o"}


def main() -> None:
    set_nature_style()
    df = pd.read_csv(IN, sep="\t")
    elig = df[df["eligible_for_5D"] == True].copy()
    if elig.empty:
        raise SystemExit("No eligible ROIs found")

    elig["color"] = elig["rhythm"].map({"AF": COL_AF, "SR": COL_SR})
    elig["marker_shape"] = elig["area"].map(SHAPE)

    fig_w = COL_DOUBLE * 0.78
    fig_h = 2.35
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.075, right=0.98, top=0.80, bottom=0.22,
                        wspace=0.45)

    # =================== LEFT: ECM by area, color = rhythm ===================
    ax_l = fig.add_subplot(1, 3, 1)
    areas = ["border", "distal"]
    rhythms = ["SR", "AF"]
    jitter_w = 0.085
    rng = np.random.default_rng(7)

    # boxplots per area (area effect is the message)
    box_data = [elig.loc[elig["area"] == a, "ecm_score"].to_numpy()
                for a in areas]
    bp = ax_l.boxplot(
        box_data, positions=[0, 1], widths=0.55, patch_artist=True,
        showcaps=False, whiskerprops=dict(color="#888888", linewidth=0.6),
        medianprops=dict(color="#222222", linewidth=1.0),
        flierprops=dict(marker=""),
        boxprops=dict(facecolor="#F2F2F2", edgecolor="#888888", linewidth=0.6),
    )

    for x_pos, area in enumerate(areas):
        for r in rhythms:
            sub = elig[(elig["area"] == area) & (elig["rhythm"] == r)]
            if sub.empty:
                continue
            xs = x_pos + (rng.random(len(sub)) - 0.5) * 2 * jitter_w
            ax_l.scatter(xs, sub["ecm_score"],
                         marker=SHAPE[area], s=22,
                         facecolor=COL_AF if r == "AF" else COL_SR,
                         edgecolor="white", linewidth=0.4, zorder=4)

    ax_l.set_xticks([0, 1])
    ax_l.set_xticklabels(["border", "distal"], fontsize=6.5)
    ax_l.set_xlabel("AA ROI position", fontsize=6.5, labelpad=2)
    ax_l.set_ylabel("ECM proteoglycans score", fontsize=6.5, labelpad=2)
    ax_l.set_title("ECM activity by AA position",
                   fontsize=7, pad=3, fontweight="bold")
    for s in ("top", "right"):
        ax_l.spines[s].set_visible(False)
    ax_l.tick_params(axis="both", labelsize=6, pad=2)

    # =================== MIDDLE: fibro vs ECM scatter ====================
    ax_m = fig.add_subplot(1, 3, 2)
    for r in rhythms:
        for a in areas:
            sub = elig[(elig["rhythm"] == r) & (elig["area"] == a)]
            if sub.empty:
                continue
            ax_m.scatter(sub["fibro_score"], sub["ecm_score"],
                         marker=SHAPE[a], s=24,
                         facecolor=COL_AF if r == "AF" else COL_SR,
                         edgecolor="white", linewidth=0.4, zorder=3,
                         label=f"{r} · {a}")
    rho, p = stats.spearmanr(elig["fibro_score"], elig["ecm_score"])
    ax_m.text(0.04, 0.96,
              f"Spearman ρ = {rho:+.2f}\n(n = {len(elig)}, descriptive)",
              transform=ax_m.transAxes, ha="left", va="top",
              fontsize=5.4, color="#222222", linespacing=1.25,
              bbox=dict(facecolor="white", alpha=0.9,
                        edgecolor="#E0E0E0", linewidth=0.35, pad=1.2))
    ax_m.set_xlabel("fibroblast (canonical) score", fontsize=6.5, labelpad=2)
    ax_m.set_ylabel("ECM proteoglycans score", fontsize=6.5, labelpad=2)
    ax_m.set_title("fibroblast-rich ROIs co-carry ECM signal",
                   fontsize=7, pad=3, fontweight="bold")
    for s in ("top", "right"):
        ax_m.spines[s].set_visible(False)
    ax_m.tick_params(axis="both", labelsize=6, pad=2)

    # =================== RIGHT: callout ====================
    ax_r = fig.add_subplot(1, 3, 3)
    ax_r.set_axis_off()
    ax_r.text(
        0.0, 0.78,
        "ECM-associated activity is enriched at the appendage "
        "border relative to distal appendage,",
        transform=ax_r.transAxes,
        ha="left", va="top", fontsize=6.0,
        color="#222222", linespacing=1.45,
    )
    ax_r.text(
        0.0, 0.62,
        "in both AF and SR; "
        "fibroblast-rich ROIs also carry higher ECM proteoglycan signal.",
        transform=ax_r.transAxes,
        ha="left", va="top", fontsize=6.0,
        color="#222222", linespacing=1.45,
    )
    ax_r.text(
        0.0, 0.32,
        "Spatial layer provides localization support only;\n"
        "mechanistic inference is not derived from GeoMx and\n"
        "relies primarily on cross-layer synthesis.",
        transform=ax_r.transAxes,
        ha="left", va="top", fontsize=5.5,
        color="#555555", style="italic", linespacing=1.5,
    )

    # =================== Shared legend (bottom) ====================
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", linestyle="None", markersize=4.5,
               markerfacecolor=COL_AF, markeredgecolor="white",
               markeredgewidth=0.4, label="AF"),
        Line2D([], [], marker="o", linestyle="None", markersize=4.5,
               markerfacecolor=COL_SR, markeredgecolor="white",
               markeredgewidth=0.4, label="SR"),
        Line2D([], [], marker="^", linestyle="None", markersize=5.5,
               markerfacecolor="#888888", markeredgecolor="white",
               markeredgewidth=0.4, label="border"),
        Line2D([], [], marker="o", linestyle="None", markersize=4.5,
               markerfacecolor="#888888", markeredgecolor="white",
               markeredgewidth=0.4, label="distal"),
    ]
    fig.legend(
        handles=handles, loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        fontsize=5.8, frameon=False, ncol=4,
        handlelength=1.05, handleheight=0.55, columnspacing=1.6,
    )

    fig.suptitle(
        "ROI-level localization of fibroblast-associated ECM proteoglycan signal "
        "(GSE261363 AA segments, n = 24)",
        fontsize=7.2, y=0.965, fontweight="bold",
    )

    # Panel label
    fig.text(0.012, 0.965, "c",
             fontsize=10, fontweight="bold", ha="left", va="top")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}  (n_ROIs={len(elig)}, ρ={rho:+.3f})")


if __name__ == "__main__":
    main()
