#!/usr/bin/env python3
"""
Figure 6D  Top reversal-gene anchor profiles.

For 4 representative reversal genes, show the per-cohort aCAI profile
across the three paired cohorts (healthy reference + SR + AF) so the
chamber-flip pattern is visible at single-cohort resolution.  Each
panel additionally shows the gene's dominant cellular carrier (from
Fig 6C, n2b_primary refresh).

Picks (locked):
    DCN     RA→LA   fibroblast carrier   strongest fibroblast Δ
    C1S     RA→LA   fibroblast carrier   complement / matrix gene
    LAPTM5  RA→LA   macrophage carrier    immune compartment exemplar
    ACTA1   LA→RA   mast_cell carrier     LA→RA flip exemplar

Inputs:
    results/n2_chamber_reversal_long.tsv
    results/n2_chamber_reversal_genes.tsv
    results/n2b_primary_reversal_per_celltype_long.tsv

Output:
    figures/fig6d_top_anchor_genes.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

LONG_TSV   = ROOT / "results" / "n2_chamber_reversal_long.tsv"
GENES_TSV  = ROOT / "results" / "n2_chamber_reversal_genes.tsv"
N2B_TSV    = ROOT / "results" / "n2b_primary_reversal_per_celltype_long.tsv"
FIG_OUT    = ROOT / "figures" / "fig6d_top_anchor_genes.pdf"

LA_COLOR  = "#C97A74"
RA_COLOR  = "#7D99B8"
GOLD      = "#D4A35E"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"
SR_FILL = "#9DBAA9"
AF_FILL = "#D67373"

COHORT_ORDER = ["GSE128188", "GSE79768", "GSE115574"]

PANELS = [
    {"gene": "DCN",    "carrier_ct": "fibroblast",
     "tag": "RA→LA  ·  fibroblast carrier"},
    {"gene": "C1S",    "carrier_ct": "fibroblast",
     "tag": "RA→LA  ·  fibroblast carrier"},
    {"gene": "LAPTM5", "carrier_ct": "macrophage",
     "tag": "RA→LA  ·  macrophage carrier"},
    {"gene": "ACTA1",  "carrier_ct": "mast_cell",
     "tag": "LA→RA  ·  mast cell carrier"},
]


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "axes.titlesize":    9,
        "axes.labelsize":    8,
        "xtick.labelsize":   7,
        "ytick.labelsize":   7,
        "legend.fontsize":   6.5,
        "axes.linewidth":    0.6,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def main() -> None:
    set_publication_style()
    long  = pd.read_csv(LONG_TSV, sep="\t")
    genes = pd.read_csv(GENES_TSV, sep="\t").set_index("gene")
    n2b   = pd.read_csv(N2B_TSV, sep="\t")

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(200 * MM, 110 * MM))
    grid = fig.add_gridspec(
        nrows=1, ncols=len(PANELS), wspace=0.30,
        left=0.06, right=0.98, top=0.84, bottom=0.16)

    fig.suptitle(
        "Fig 6D  Top reversal-gene anchor profiles  "
        "(per-cohort aCAI, healthy / SR / AF)",
        fontsize=8.8, fontweight="bold", y=0.965, x=0.02, ha="left")

    for idx, p in enumerate(PANELS):
        ax = fig.add_subplot(grid[idx])
        gene = p["gene"]
        sub = long[long["gene"] == gene].copy()
        sub["__co"] = sub["cohort"].map(
            {c: i for i, c in enumerate(COHORT_ORDER)})
        sub = sub.sort_values("__co")

        if gene not in genes.index:
            ax.text(0.5, 0.5, f"NO RECORD: {gene}",
                     transform=ax.transAxes, ha="center", va="center")
            continue
        g_meta = genes.loc[gene]
        healthy = float(g_meta["aCAI_healthy"])
        h_dir   = g_meta["healthy_direction"]
        h_color = LA_COLOR if h_dir == "LA-biased" else RA_COLOR

        # x layout: 0 healthy, 1-3 cohort SR, 5-7 cohort AF
        x_h = 0.0
        x_sr = np.array([1.5, 2.5, 3.5])
        x_af = np.array([5.0, 6.0, 7.0])

        # background blocks
        ax.axvspan(-0.6, 0.6, color=GOLD, alpha=0.10, zorder=0)
        ax.axvspan(1.0, 4.0, color=SR_FILL, alpha=0.13, zorder=0)
        ax.axvspan(4.5, 7.6, color=AF_FILL, alpha=0.13, zorder=0)

        # zero line
        ax.axhline(0, color=GREY_DARK, lw=0.5, ls="--", zorder=1)

        # healthy point
        ax.scatter([x_h], [healthy],
                      marker="D", s=80,
                      facecolor=GOLD, edgecolor="black", linewidth=0.6,
                      zorder=4, label="bulk healthy")
        ax.text(x_h, healthy + 0.05 * np.sign(healthy + 1e-9),
                  f"{healthy:+.2f}",
                  ha="center", va="bottom" if healthy > 0 else "top",
                  fontsize=6.0, color=GREY_DARK)

        # SR + AF dots per cohort + connecting line
        sr_vals = []; af_vals = []; cohort_labels = []
        for ci, cohort in enumerate(COHORT_ORDER):
            row = sub[sub["cohort"] == cohort]
            if row.empty:
                sr_vals.append(np.nan); af_vals.append(np.nan)
                cohort_labels.append(cohort); continue
            r = row.iloc[0]
            sr_vals.append(float(r["aCAI_SR"]) if pd.notna(r["aCAI_SR"]) else np.nan)
            af_vals.append(float(r["aCAI_AF"]) if pd.notna(r["aCAI_AF"]) else np.nan)
            cohort_labels.append(cohort)

        for k in range(3):
            if not pd.isna(sr_vals[k]):
                ax.scatter([x_sr[k]], [sr_vals[k]],
                              marker="o", s=55,
                              facecolor=SR_FILL, edgecolor="black",
                              linewidth=0.5, zorder=4,
                              label="SR" if k == 0 else None)
                ax.text(x_sr[k], sr_vals[k] + 0.05 * np.sign(sr_vals[k] + 1e-9),
                          f"{sr_vals[k]:+.2f}",
                          ha="center", va="bottom" if sr_vals[k] > 0 else "top",
                          fontsize=5.5, color=GREY_DARK)
            if not pd.isna(af_vals[k]):
                ax.scatter([x_af[k]], [af_vals[k]],
                              marker="s", s=55,
                              facecolor=AF_FILL, edgecolor="black",
                              linewidth=0.5, zorder=4,
                              label="AF" if k == 0 else None)
                ax.text(x_af[k], af_vals[k] + 0.05 * np.sign(af_vals[k] + 1e-9),
                          f"{af_vals[k]:+.2f}",
                          ha="center", va="bottom" if af_vals[k] > 0 else "top",
                          fontsize=5.5, color=GREY_DARK)
            # connecting line SR → AF (per cohort)
            if not (pd.isna(sr_vals[k]) or pd.isna(af_vals[k])):
                ax.plot([x_sr[k], x_af[k]],
                          [sr_vals[k], af_vals[k]],
                          color=GREY_MID, lw=0.6, ls=":", zorder=2.5)

        # X tick labels: H | cohorts SR | cohorts AF
        all_x = [x_h] + list(x_sr) + list(x_af)
        all_labels = ["bulk", "1", "2", "3", "1", "2", "3"]
        ax.set_xticks(all_x)
        ax.set_xticklabels(all_labels, fontsize=6.0)
        # group labels under axis
        ax.text(x_h, ax.get_ylim()[0], "",
                  ha="center", va="top", fontsize=6.0)
        # bottom block titles
        ax.annotate("Healthy",
                       (x_h, 0), xytext=(0, -22), textcoords="offset points",
                       ha="center", fontsize=6.5, fontweight="bold",
                       color=GOLD)
        ax.annotate("SR  (cohort 1-3)",
                       ((x_sr[0] + x_sr[-1]) / 2, 0),
                       xytext=(0, -22), textcoords="offset points",
                       ha="center", fontsize=6.5, fontweight="bold",
                       color="#3F6B59")
        ax.annotate("AF  (cohort 1-3)",
                       ((x_af[0] + x_af[-1]) / 2, 0),
                       xytext=(0, -22), textcoords="offset points",
                       ha="center", fontsize=6.5, fontweight="bold",
                       color="#9E3636")

        ax.set_xlim(-0.6, 7.6)
        # Y limits with padding
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin - 0.20, ymax + 0.20)

        # Carrier badge top-right
        cell = p["carrier_ct"]
        sub_ct = n2b[(n2b["gene"] == gene)
                       & (n2b["compartment"] == cell)]
        if not sub_ct.empty:
            d_align = float(sub_ct.iloc[0]["delta_aligned"])
            badge = (f"carrier: {cell.replace('_', ' ')}\n"
                       f"Δ_aligned = {d_align:+.3f}")
            ax.text(0.98, 0.98, badge,
                      transform=ax.transAxes,
                      ha="right", va="top",
                      fontsize=6.0, fontweight="bold",
                      color="#1F5C39",
                      bbox=dict(boxstyle="round,pad=0.30",
                                  facecolor="#F2F8F4",
                                  edgecolor="#3C8C57", linewidth=0.5))

        # Title
        ax.set_title(
            f"{gene}    ({h_dir})",
            loc="left", fontsize=8.0, fontweight="bold",
            color=h_color, pad=4)
        ax.text(0.02, 0.95, p["tag"],
                  transform=ax.transAxes, ha="left", va="top",
                  fontsize=6.5, color=GREY_DARK, style="italic")
        ax.set_ylabel("aCAI  (LA − RA)", fontsize=7)

    fig.text(0.02, 0.025,
              "Diamond = bulk healthy aCAI (8-donor reference, GSE112339).   "
              "Circle = SR per cohort.   Square = AF per cohort.   "
              "Dotted line connects SR → AF within each cohort.   "
              "Reversal pattern:  AF colour-flips against bulk healthy "
              "while SR preserves it.",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
