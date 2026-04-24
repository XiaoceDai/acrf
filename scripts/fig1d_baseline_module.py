#!/usr/bin/env python3
"""
Figure 1d (merged baseline module) - HORIZONTAL layout for composite assembly.

Layout (top -> bottom, all sharing the program x-axis):
  row 0: pathway short labels (rotated 90 deg)
  row 1: LA / RA direction strip (with text)
  row 2: donor x program aCAI heatmap (donors on y, programs on x)
  row 3: family color band with F-id
  row 4: mean aCAI bars + n/8 stability annotations

Right side: a small horizontal colorbar at the bottom.

Inputs:
  results/fig1e_baseline_axis_families.tsv

Output:
  figures/fig1d_baseline_module.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import matplotlib as mpl

from plot_style import COL_DOUBLE, set_nature_style
from _pathway_labels import short_pathway_label

set_nature_style()

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

IN_TSV = ROOT / "results" / "fig1e_baseline_axis_families.tsv"
OUT_PDF = ROOT / "figures" / "fig1d_baseline_module.pdf"

LA_COLOR = "#D6604D"
RA_COLOR = "#4393C3"

FAMILY_PALETTE = [
    "#E6AB02", "#66A61E", "#7570B3", "#D95F02",
    "#1B9E77", "#E7298A", "#A6761D", "#666666",
    "#1F78B4", "#FB9A99", "#33A02C", "#B15928",
]


def short_label(s: str) -> str:
    return short_pathway_label(s, 26)


def main() -> None:
    df = pd.read_csv(IN_TSV, sep="\t")
    donor_cols = [c for c in df.columns if c.startswith("aCAI_")]
    donors = [c.replace("aCAI_", "") for c in donor_cols]

    # Canonical row order (LA top -> RA bottom): we want RA on the LEFT in the
    # horizontal layout per user's reference. Reverse the order: RA first, LA next.
    fam_la = df[df["direction"] == "LA-biased"].copy()
    fam_ra = df[df["direction"] == "RA-biased"].copy()
    sorted_df = pd.concat([fam_ra, fam_la], ignore_index=True)
    K = len(sorted_df)
    n_ra = len(fam_ra)

    # rows = donors, cols = programs
    mat = sorted_df[donor_cols].to_numpy(dtype=float).T
    n_donors = len(donors)
    vmax = float(np.nanmax(np.abs(mat))) or 1e-3
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list(
        "LA_RA",
        [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)],
    )

    fam_ids = sorted(sorted_df["family_id"].unique())
    fam_cmap = {fid: FAMILY_PALETTE[i % len(FAMILY_PALETTE)]
                for i, fid in enumerate(fam_ids)}

    # ---- figure layout -------------------------------------------------------
    fig = plt.figure(figsize=(COL_DOUBLE * 1.55, COL_DOUBLE * 0.58))
    gs = fig.add_gridspec(
        5, 1,
        height_ratios=[1.6, 0.30, 1.7, 0.32, 1.4],
        hspace=0.10,
        top=0.97, bottom=0.10, left=0.07, right=0.985,
    )

    # row 0 : angled bold pathway labels (top)
    ax_lbl = fig.add_subplot(gs[0, 0])
    ax_lbl.set_xlim(-0.5, K - 0.5)
    ax_lbl.set_ylim(0, 1.12)
    ax_lbl.axis("off")
    for i, r in sorted_df.iterrows():
        dcol = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        ax_lbl.text(
            i, 0.02, short_label(r["pathway"]),
            rotation=52, ha="left", va="bottom",
            rotation_mode="anchor",
            fontsize=7.3, color=dcol, fontweight="bold",
        )

    # row 1 : LA / RA direction strip (with text)
    ax_dir = fig.add_subplot(gs[1, 0], sharex=ax_lbl)
    for i, d in enumerate(sorted_df["direction"]):
        ax_dir.add_patch(mpatches.Rectangle(
            (i - 0.5, 0), 1, 1,
            facecolor=LA_COLOR if d == "LA-biased" else RA_COLOR,
            edgecolor="white", linewidth=0.2,
        ))
    ax_dir.text(n_ra / 2 - 0.5, 0.5, "RA",
                ha="center", va="center", color="white",
                fontweight="bold", fontsize=8.6)
    ax_dir.text(n_ra + (K - n_ra) / 2 - 0.5, 0.5, "LA",
                ha="center", va="center", color="white",
                fontweight="bold", fontsize=8.6)
    ax_dir.set_xlim(-0.5, K - 0.5)
    ax_dir.set_ylim(0, 1)
    ax_dir.axis("off")

    # row 2 : main donor x program heatmap (each cell as independent vector)
    ax_hm = fig.add_subplot(gs[2, 0], sharex=ax_lbl)
    for i in range(n_donors):
        for j in range(K):
            v = mat[i, j]
            if np.isnan(v):
                fc = "#E5E5E5"
            else:
                fc = cmap(norm(v))
            ax_hm.add_patch(mpatches.Rectangle(
                (j - 0.5, i - 0.5), 1.0, 1.0,
                facecolor=fc, edgecolor="white", linewidth=0.18,
            ))
    ax_hm.set_xlim(-0.5, K - 0.5)
    ax_hm.set_ylim(n_donors - 0.5, -0.5)
    ax_hm.set_yticks(np.arange(n_donors))
    ax_hm.set_yticklabels(donors, fontsize=6.8)
    ax_hm.set_ylabel("donor", fontsize=8.0)
    ax_hm.tick_params(axis="y", length=0, pad=2)
    ax_hm.set_xticks([])
    for s in ("top", "right", "left", "bottom"):
        ax_hm.spines[s].set_visible(False)
    # divider between RA and LA program blocks
    ax_hm.axvline(n_ra - 0.5, color="black", lw=0.8)

    # row 3 : family color band with F-id (along program axis)
    ax_fam = fig.add_subplot(gs[3, 0], sharex=ax_lbl)
    for i, fid in enumerate(sorted_df["family_id"]):
        ax_fam.add_patch(mpatches.Rectangle(
            (i - 0.5, 0), 1, 1,
            facecolor=fam_cmap[fid], edgecolor="white", linewidth=0.2,
        ))
    fam_groups = sorted_df.groupby("family_id", sort=False)
    for fid, gp in fam_groups:
        idxs = sorted_df.index[sorted_df["family_id"] == fid].tolist()
        xc = float(np.mean(idxs))
        ax_fam.text(xc, 0.5, f"F{int(fid)}",
                    ha="center", va="center",
                    fontsize=7.2, color="white", fontweight="bold")
    ax_fam.set_xlim(-0.5, K - 0.5)
    ax_fam.set_ylim(0, 1)
    ax_fam.axis("off")

    # row 4 : mean aCAI bars + n/8 stability
    ax_bar = fig.add_subplot(gs[4, 0], sharex=ax_lbl)
    bar_colors = [LA_COLOR if d == "LA-biased" else RA_COLOR
                  for d in sorted_df["direction"]]
    x = np.arange(K)
    ax_bar.bar(x, sorted_df["mean_aCAI"], color=bar_colors,
               edgecolor="white", linewidth=0.3, width=0.85)
    ax_bar.axhline(0, color="#444", lw=0.4)
    ax_bar.set_xlim(-0.5, K - 0.5)
    ax_bar.set_xticks([])
    ax_bar.tick_params(axis="y", labelsize=6.8)
    ax_bar.set_ylabel("mean aCAI", fontsize=8.0)
    for s in ("top", "right"):
        ax_bar.spines[s].set_visible(False)
    ax_bar.axvline(n_ra - 0.5, color="black", lw=0.8)

    ymax_abs = float(np.nanmax(np.abs(sorted_df["mean_aCAI"]))) or 0.1
    ypad = ymax_abs * 0.10
    for xi, v, sc in zip(x, sorted_df["mean_aCAI"], sorted_df["sign_concordance"]):
        ya = v + (ypad if v >= 0 else -ypad)
        va = "bottom" if v >= 0 else "top"
        ax_bar.text(
            xi, ya, f"{int(sc)}/8",
            fontsize=6.0, ha="center", va=va,
            color="#1B7837" if sc == 8 else "#666",
            fontweight="bold" if sc == 8 else "normal",
            rotation=90,
        )
    # leave a little headroom for rotated n/8 text
    ax_bar.set_ylim(-ymax_abs * 1.45, ymax_abs * 1.45)

    # ---- colorbar at the very bottom right ----
    cax = fig.add_axes([0.78, 0.045, 0.18, 0.018])
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label("aCAI", fontsize=7.4, labelpad=2)
    cbar.ax.tick_params(labelsize=6.6)
    cbar.ax.text(-0.04, 0.5, "RA", transform=cbar.ax.transAxes,
                 ha="right", va="center", fontsize=6.8,
                 color=RA_COLOR, fontweight="bold")
    cbar.ax.text(1.04, 0.5, "LA", transform=cbar.ax.transAxes,
                 ha="left", va="center", fontsize=6.8,
                 color=LA_COLOR, fontweight="bold")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.03)
    plt.close(fig)
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
