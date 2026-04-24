#!/usr/bin/env python3
"""
Figure 1A - GSE112339 bulk PCA, LA vs RA.

Rebuilt in Python so typography / colour convention / panel sizing match
the other Figure 1 panels (C, D).

Data:
    Project_P01/qc/qc_pca.tsv   columns: sample_id, PC1, PC2, donor, chamber

Visual rules:
    RA = blue  (#4393C3)   LA = red  (#D6604D)
    LA filled, RA filled; donor pair joined by a thin grey line.
    Small donor id label near each LA point.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_style import COL_SINGLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
IN   = ROOT / "bulk/Project_P01/qc/qc_pca.tsv"
OUT  = ROOT / "figures" / "fig1a_pca_la_vs_ra.pdf"

# Match Figure 6B red/blue tones while swapping chamber mapping for this panel:
# RA -> blue, LA -> red.
RA_COLOR = "#4393C3"
LA_COLOR = "#D6604D"


def load_pca() -> pd.DataFrame:
    df = pd.read_csv(IN, sep="\t")
    # strip chamber prefix from donor-ish sample_id to get a short tag
    df["donor_tag"] = df["sample_id"].str.replace(r"^[LR]A_", "", regex=True)
    return df


def main() -> None:
    set_nature_style()
    df = load_pca()

    # Variance explained isn't in qc_pca.tsv; fall back to showing raw PC scores.
    # We can at least report axis ranges to reviewers via PC1/PC2 labels.

    fig_w = COL_SINGLE
    fig_h = COL_SINGLE * 0.95
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.15)

    # Draw donor-pair connector lines first (under points).
    for _, g in df.groupby("donor_tag"):
        if len(g) == 2:
            ax.plot(g["PC1"], g["PC2"],
                    color="#BEBEBE", linewidth=0.5, zorder=1)

    la = df[df["chamber"] == "LA"]
    ra = df[df["chamber"] == "RA"]

    ax.scatter(ra["PC1"], ra["PC2"],
               s=34, facecolor=RA_COLOR, edgecolor="white",
               linewidth=0.6, label="RA", zorder=3)
    ax.scatter(la["PC1"], la["PC2"],
               s=34, facecolor=LA_COLOR, edgecolor="white",
               linewidth=0.6, label="LA", zorder=3)

    # Donor labels -- one per donor, anchored to the LA point.
    for _, row in la.iterrows():
        ax.annotate(row["donor_tag"],
                    xy=(row["PC1"], row["PC2"]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=5, color="#444444")

    # Axes / spines.
    ax.axhline(0, color="#D0D0D0", lw=0.4, zorder=0)
    ax.axvline(0, color="#D0D0D0", lw=0.4, zorder=0)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # Symmetric legend in upper right, no frame.
    leg = ax.legend(loc="upper right", frameon=False,
                    handletextpad=0.3, borderpad=0.2)
    for handle in leg.legend_handles:
        handle.set_sizes([18])

    ax.set_title("a  PCA (GSE112339)", fontsize=7.6, pad=4)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
