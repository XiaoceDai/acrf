#!/usr/bin/env python3
"""
Figure 1d (SUMMARY version, for the main composite Figure 1).

Drops the 30 x 8 per-donor heatmap to reduce visual fatigue across the
Figure 1 panels (Fig 1d / 1e / 1f all used 30-row heatmaps).

What this panel shows, for each NCP_v2 program:
    - mean_aCAI averaged over 8 healthy donors  (bar length / direction)
    - sign-concordance n/8                       (text annotation)
    - chamber direction tag (LA = red, RA = blue)
    - pathway short name on the right

The full per-donor heatmap is kept as a supplementary figure
(fig1b_baseline_donor_aCAI.pdf).

Inputs:
    results/fig1b_baseline_donor_aCAI.tsv   (cached by fig1b_baseline_donor_aCAI.py)

Output:
    figures/fig1d_summary.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import COL_DOUBLE, set_nature_style
from _pathway_labels import short_pathway_label

set_nature_style()

mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
})

LA_COLOR = "#D6604D"   # LA = red (new convention)
RA_COLOR = "#4393C3"   # RA = blue

IN  = ROOT / "results" / "fig1b_baseline_donor_aCAI.tsv"
OUT = ROOT / "figures" / "fig1d_summary.pdf"


def main() -> None:
    df = pd.read_csv(IN, sep="\t")

    # Keep canonical sort: LA first (most LA on top), then RA (most RA on bottom)
    sorted_df = pd.concat([
        df[df["direction"] == "LA-biased"]
          .sort_values("mean_aCAI", ascending=False),
        df[df["direction"] == "RA-biased"]
          .sort_values("mean_aCAI", ascending=True),
    ]).reset_index(drop=True)

    n_la  = (sorted_df["direction"] == "LA-biased").sum()
    n_ra  = (sorted_df["direction"] == "RA-biased").sum()
    n_8of8 = int((sorted_df["sign_concordance"] == 8).sum())
    K = len(sorted_df)

    # Compact landscape figure: direction band | bar chart | pathway name
    fig = plt.figure(figsize=(COL_DOUBLE * 0.62, COL_DOUBLE * 0.85))
    gs = fig.add_gridspec(
        1, 3, width_ratios=[0.20, 1.10, 1.55], wspace=0.06,
        top=0.88, bottom=0.12, left=0.05, right=0.985,
    )

    # ---- col 0 : LA / RA direction band -----------------------------------
    ax_dir = fig.add_subplot(gs[0, 0])
    for i, d in enumerate(sorted_df["direction"]):
        ax_dir.add_patch(mpatches.Rectangle(
            (0.0, i - 0.5), 1.0, 1.0,
            facecolor=LA_COLOR if d == "LA-biased" else RA_COLOR,
            edgecolor="white", linewidth=0.2,
        ))
    ax_dir.text(0.5, n_la / 2 - 0.5, "LA-biased",
                rotation=90, ha="center", va="center",
                color="white", fontweight="bold", fontsize=8.5)
    ax_dir.text(0.5, n_la + n_ra / 2 - 0.5, "RA-biased",
                rotation=90, ha="center", va="center",
                color="white", fontweight="bold", fontsize=8.5)
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(K - 0.5, -0.5)
    ax_dir.set_xticks([]); ax_dir.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        ax_dir.spines[s].set_visible(False)

    # ---- col 1 : mean aCAI bars + n/8 ------------------------------------
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax_dir)
    colors = [LA_COLOR if d == "LA-biased" else RA_COLOR
              for d in sorted_df["direction"]]
    y = np.arange(K)
    ax2.barh(y, sorted_df["mean_aCAI"], color=colors,
             edgecolor="white", linewidth=0.3)
    ax2.axvline(0, color="#444", lw=0.4)
    ax2.set_xlabel("mean aCAI\nover 8 healthy donors", fontsize=8.5)
    ax2.tick_params(axis="x", labelsize=7.5)
    ax2.set_yticks([])
    ax2.axhline(n_la - 0.5, color="black", lw=0.7)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    # NOTE: sharey already inverts; no invert_yaxis() needed.

    # n/8 sign-concordance label outside each bar end
    xpad = 0.02
    for yi, v, sc in zip(y, sorted_df["mean_aCAI"],
                          sorted_df["sign_concordance"]):
        sc_color = "#1B7837" if sc == 8 else "#666"
        sc_weight = "bold" if sc == 8 else "normal"
        ax2.text(v + (xpad if v > 0 else -xpad), yi,
                  f"{int(sc)}/8",
                  fontsize=7, ha="left" if v > 0 else "right",
                  va="center", color=sc_color, fontweight=sc_weight)

    # ---- col 2 : pathway names (coloured by chamber direction) -----------
    ax_lab = fig.add_subplot(gs[0, 2], sharey=ax_dir)
    for i, (_, r) in enumerate(sorted_df.iterrows()):
        col = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        ax_lab.text(0.0, i, short_pathway_label(r["pathway"], 38),
                    ha="left", va="center", fontsize=7.8, color=col)
    ax_lab.axhline(n_la - 0.5, color="black", lw=0.7)
    ax_lab.set_xlim(0, 1)
    ax_lab.set_xticks([])
    ax_lab.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        ax_lab.spines[s].set_visible(False)

    # ---- title + summary stat ---------------------------------------------
    fig.suptitle(
        f"d  Donor-level aCAI primitive — {n_8of8}/{K} programs sign-concordant in 8/8 donors",
        x=0.05, y=0.96, ha="left", va="top",
        fontsize=9.5, fontweight="bold",
    )
    fig.text(
        0.05, 0.92,
        "aCAI = log\u2082(LA) \u2212 log\u2082(RA), averaged across 8 healthy donors  \u00b7  "
        "n/8 = donors sharing the row mean's sign  \u00b7  8/8 (green, bold) = stable",
        fontsize=7.2, color="#444", ha="left", va="top",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
