#!/usr/bin/env python3
"""
Figure 3D - donor-level singscore dots for representative programs,
faceted by cell type.

Spec §6.4 selection is relaxed because under Spec A (q_perm < 0.10) no
program has a clear major carrier; we instead pick four programs by
narrative role:

  1. LA-biased metabolic with CM-vs-immune dichotomy
         -> Cytoplasmic Translation
  2. RA-biased ECM, the only * hit (ECM proteoglycans x fibroblast)
         -> ECM proteoglycans
  3. Cell-type-reversed (CM up, immune down)
         -> BMP Signalling
  4. Flat / preserved-ish
         -> OXPHOS
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from plot_style import (
    AFCTRL_CMAP, COL_DOUBLE, COMPARTMENT_COLORS, set_nature_style,
)

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "figures" / "fig3d_representative_dots.pdf"

COMPARTMENTS = ["atrialCM", "fibroblast", "endothelial_or_vascular", "immune"]
COMPARTMENT_LABELS = {
    "atrialCM": "atrialCM",
    "fibroblast": "fibroblast",
    "endothelial_or_vascular": "endo/vasc",
    "immune": "immune",
}

REPRESENTATIVES = [
    ("REACTOME__REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM",  "nominal fibroblast-localized signal"),
    ("REACTOME__REACTOME_SIGNALING_BY_BMP",                                       "compartmentally discordant shift"),
    ("REACTOME__REACTOME_TRANSLATION",                                            "minimal LA-arm shift"),
    ("HALLMARK__HALLMARK_OXIDATIVE_PHOSPHORYLATION",                              "minimal LA-arm shift"),
]

# rhythm cosmetics
RHY_COL = {"Ctrl": "#666666", "AF": AFCTRL_CMAP[-1]}  # teal-ish for Ctrl, brown for AF
RHY_MARK = {"Ctrl": "o", "AF": "s"}


def ci95(values: np.ndarray) -> tuple[float, float]:
    n = len(values)
    if n < 2:
        return (float(values.mean()) if n else float("nan"), float("nan"))
    m = float(values.mean())
    se = float(values.std(ddof=1) / math.sqrt(n))
    half = float(stats.t.ppf(0.975, df=n - 1) * se)
    return m, half


def sig_marker(p: float) -> str:
    if math.isnan(p):
        return ""
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "."
    return ""


def load_data():
    src = pd.read_csv(ROOT / "results" / "figure3_source_data.tsv", sep="\t")
    delta = pd.read_csv(ROOT / "results" / "af_la_program_celltype_delta.tsv",
                        sep="\t")
    return src, delta


def plot_panel(ax_row, pick_id, role, src, delta, y_lim):
    prog = src[src["program_id"] == pick_id].iloc[0]
    prog_label = prog["program_label"]
    base_dir   = int(prog["baseline_direction"])

    sub = src[src["program_id"] == pick_id]
    d_stats = delta[delta["program_id"] == pick_id].set_index("cell_type")

    for idx, ct in enumerate(COMPARTMENTS):
        ax = ax_row[idx]
        ct_data = sub[sub["cell_type"] == ct]
        for rhy_idx, rhy in enumerate(["Ctrl", "AF"]):
            vals = ct_data.loc[ct_data["rhythm"] == rhy, "singscore"].to_numpy()
            if len(vals) == 0:
                continue
            rng = np.random.default_rng(42 + idx * 10 + rhy_idx)
            jitter = rng.uniform(-0.14, 0.14, size=len(vals))
            ax.scatter(rhy_idx + jitter, vals,
                       marker=RHY_MARK[rhy], s=11,
                       facecolor=RHY_COL[rhy], edgecolor="white",
                       linewidth=0.3, alpha=0.85, zorder=3)
            m, half = ci95(vals)
            ax.plot([rhy_idx - 0.28, rhy_idx + 0.28], [m, m],
                    color=RHY_COL[rhy], lw=1.0, zorder=4)
            if not math.isnan(half):
                ax.plot([rhy_idx, rhy_idx], [m - half, m + half],
                        color=RHY_COL[rhy], lw=0.7, zorder=4)

        # x axis
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Ctrl", "AF"], fontsize=6)
        ax.tick_params(axis="x", length=0, pad=2)
        ax.set_xlim(-0.55, 1.55)
        ax.set_ylim(*y_lim)
        ax.tick_params(axis="y", labelsize=6)

        # cell-type label strip above each small axis
        from matplotlib.transforms import blended_transform_factory
        trans = blended_transform_factory(ax.transAxes, ax.transAxes)
        ax.add_patch(plt.Rectangle((0.0, 1.02), 1.0, 0.065,
                                   transform=trans,
                                   facecolor=COMPARTMENT_COLORS[ct],
                                   edgecolor="none", clip_on=False,
                                   zorder=5))
        ax.text(0.5, 1.055, COMPARTMENT_LABELS[ct],
                transform=trans, ha="center", va="center",
                fontsize=6, color="white", fontweight="bold",
                clip_on=False, zorder=6)

        # bracket + sig_marker if any.  Draw marker at the very top of the
        # panel's plot area (in axes fraction coords so it does not fight
        # the shared y-limits, which are dominated by outlier points).
        if ct in d_stats.index:
            p = d_stats.at[ct, "p_perm"]
            marker = sig_marker(p)
            if marker:
                from matplotlib.transforms import blended_transform_factory
                tf = blended_transform_factory(ax.transData, ax.transAxes)
                ax.plot([0, 0, 1, 1], [0.78, 0.82, 0.82, 0.78],
                        color="#222222", lw=0.8, transform=tf,
                        clip_on=False)
                ax.text(0.5, 0.84, marker,
                        ha="center", va="bottom", fontsize=12,
                        fontweight="bold", color="#222222",
                        transform=tf, clip_on=False)

        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.spines["left"].set_linewidth(0.5)
        ax.spines["bottom"].set_linewidth(0.5)

    # leftmost axis gets a y-label showing the program + role
    y_text = f"{prog_label}\n{role}"
    ax_row[0].set_ylabel(y_text, fontsize=6.5, labelpad=8)


def main() -> None:
    set_nature_style()
    src, delta = load_data()

    # validate representatives exist
    available = set(src["program_id"].unique())
    picks = [(pid, role) for pid, role in REPRESENTATIVES if pid in available]
    if len(picks) == 0:
        print("no representatives available; skipping fig3d")
        return

    # shared y range
    rep_src = src[src["program_id"].isin([p for p, _ in picks])]
    y_lo = float(rep_src["singscore"].min()) - 0.01
    y_hi = float(rep_src["singscore"].max()) + 0.02

    n_rows = len(picks)
    n_cols = len(COMPARTMENTS)
    fig_w = COL_DOUBLE
    fig_h = 1.4 * n_rows + 0.7
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h),
                             sharey=True)
    if n_rows == 1:
        axes = np.array([axes])
    fig.subplots_adjust(left=0.13, right=0.96, top=0.92, bottom=0.10,
                        wspace=0.14, hspace=0.55)

    for ax_row, (pid, role) in zip(axes, picks):
        plot_panel(ax_row, pid, role, src, delta, (y_lo, y_hi))

    # legend at the top-right as figure text
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", color=RHY_COL["Ctrl"],
               linestyle="None", markersize=5,
               markerfacecolor=RHY_COL["Ctrl"],
               markeredgecolor="white", markeredgewidth=0.4,
               label="Ctrl (n=16)"),
        Line2D([], [], marker="s", color=RHY_COL["AF"],
               linestyle="None", markersize=5,
               markerfacecolor=RHY_COL["AF"],
               markeredgecolor="white", markeredgewidth=0.4,
               label="AF   (n=18)"),
    ]
    fig.legend(handles=handles, loc="upper right",
               bbox_to_anchor=(0.98, 0.99),
               fontsize=6, frameon=False,
               handlelength=1.0, handleheight=0.8)

    fig.suptitle("Representative programs: donor-level singscore by compartment",
                 fontsize=8, y=0.985)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}  "
          f"({len(picks)} representative programs)")


if __name__ == "__main__":
    main()
