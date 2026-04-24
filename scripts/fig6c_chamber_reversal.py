#!/usr/bin/env python3
"""
Figure 6c (REPLACEMENT) — Framework-discovered chamber-identity REVERSAL
genes as gene-level positive evidence for the distributed-substrate model.

Background:
    The original Figure 6c reported pre-registered cis-eQTL × AF-GWAS
    colocalisation on four anchor genes (PRELP, OMD, IDH2, NDUFA8) plus
    two canonical positive-control AF GWAS loci (PITX2, KCNN3).  At
    current GTEx Atrial Appendage panel power (n = 372), the conservative
    cis-eQTL gate (lead p < 5e-8) was not passed by either the framework
    anchors or the canonical positive controls; that classical
    single-gene mRNA-mediated regime check is now reported in
    Supplementary Fig. S12.

    Figure 6c instead reports a gene-level positive object class that
    is uniquely enabled by CRF: chamber-identity REVERSAL genes,
    defined per gene as genes whose chamber direction in healthy
    donors is flipped in AF donors of paired-chamber cohorts while SR
    donors of the same cohorts preserve the healthy direction.

Inputs:
    results/n2_chamber_reversal_genes.tsv      (filtered final list)
    results/n2_chamber_reversal_long.tsv       (per-gene per-cohort)

Output:
    figures/fig6c_chamber_reversal.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from plot_style import (
    COL_DOUBLE, LA_COLOR, RA_COLOR, panel_label,
    set_nature_style, strip_spines,
)

ROOT = Path(__file__).resolve().parent.parent
IN_GENES = ROOT / "results" / "n2_chamber_reversal_genes.tsv"
IN_LONG  = ROOT / "results" / "n2_chamber_reversal_long.tsv"
OUT      = ROOT / "figures" / "fig6c_chamber_reversal.pdf"

ANCHOR_HIGHLIGHT = {"PRELP", "OMD", "IDH2", "NDUFA8"}


def main() -> None:
    set_nature_style()
    final = pd.read_csv(IN_GENES, sep="\t")
    long  = pd.read_csv(IN_LONG, sep="\t")

    final["is_anchor"] = final["gene"].isin(ANCHOR_HIGHLIGHT)
    final = final.sort_values(
        ["n_cohorts_call", "mean_abs_aCAI_AF"],
        ascending=[False, False],
    ).reset_index(drop=True)

    cohorts = ["GSE128188", "GSE79768", "GSE115574"]
    n_la = (final["healthy_direction"] == "LA-biased").sum()
    n_ra = (final["healthy_direction"] == "RA-biased").sum()

    fig = plt.figure(figsize=(COL_DOUBLE, 0.55 * COL_DOUBLE))
    gs = fig.add_gridspec(
        2, 2, hspace=0.5, wspace=0.40,
        height_ratios=[1.0, 1.05], width_ratios=[1.5, 1.0],
    )

    # ---- panel left-top : ranked horizontal lollipop of all reversal genes
    ax = fig.add_subplot(gs[0, 0])
    n_show = len(final)
    y = np.arange(n_show)[::-1]
    colors = [LA_COLOR if d == "LA-biased" else RA_COLOR
              for d in final["healthy_direction"]]
    # lollipop stems
    for yi, mag, col in zip(y, final["mean_abs_aCAI_AF"], colors):
        ax.plot([0, mag], [yi, yi], color=col, lw=0.8, alpha=0.65, zorder=2)
    ax.scatter(final["mean_abs_aCAI_AF"], y, s=22, c=colors,
               edgecolor="white", linewidth=0.4, zorder=3)
    # bold-mark anchor-name overlap
    for yi, gene, is_anc in zip(y, final["gene"], final["is_anchor"]):
        weight = "bold" if is_anc else "normal"
        size = 5.6 if is_anc else 5.0
        ax.text(-0.005, yi, gene, fontsize=size, ha="right", va="center",
                color="#1A1A1A", fontweight=weight)
    ax.set_xlabel("mean |aCAI AF| across paired cohorts (log2 LA − RA)")
    ax.set_yticks([]); ax.set_xlim(-0.18, max(final["mean_abs_aCAI_AF"]) * 1.10)
    ax.set_title(
        "Framework-discovered chamber-identity reversal genes "
        f"(n = {len(final)}, ≥ 2/3 paired cohorts; bold = also Fig 6 anchor)",
        loc="left", pad=4, fontsize=6.6,
    )
    strip_spines(ax)
    panel_label(ax, "c", x=-0.15, y=1.04)

    # ---- panel right-top : direction asymmetry mini-bar
    ax = fig.add_subplot(gs[0, 1])
    bars = ax.bar(["LA → RA", "RA → LA"], [n_la, n_ra],
                  color=[LA_COLOR, RA_COLOR],
                  edgecolor="white", linewidth=0.4)
    for b, n in zip(bars, [n_la, n_ra]):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.5,
                f"n = {int(n)}", ha="center", va="bottom",
                fontsize=7, fontweight="bold")
    ax.set_ylabel("reversal genes")
    ax.set_ylim(0, max(n_la, n_ra) * 1.25)
    ax.set_title(
        "Reversal direction is RA → LA dominated\n"
        f"(asymmetry ratio {n_ra}:{n_la})",
        loc="left", pad=4, fontsize=6.6,
    )
    strip_spines(ax)

    # ---- panel bottom : per-cohort scatter of healthy aCAI vs AF aCAI
    inner = gs[1, :].subgridspec(1, 3, wspace=0.30)
    keep_set = set(final["gene"])
    for k, c in enumerate(cohorts):
        ax = fig.add_subplot(inner[0, k])
        sub = long.query("cohort == @c").set_index("gene")
        x = sub["aCAI_healthy"].values
        y = sub["aCAI_AF"].values
        ax.scatter(x, y, s=1.0, color="#CFCFCF", alpha=0.4,
                   edgecolor="none", zorder=1)
        is_kept = sub.index.isin(keep_set)
        is_la = is_kept & (sub["healthy_dir"] > 0)
        is_ra = is_kept & (sub["healthy_dir"] < 0)
        ax.scatter(x[is_la], y[is_la], s=14, color=LA_COLOR, alpha=0.95,
                   edgecolor="white", linewidth=0.3, zorder=3,
                   label=f"LA → RA")
        ax.scatter(x[is_ra], y[is_ra], s=14, color=RA_COLOR, alpha=0.95,
                   edgecolor="white", linewidth=0.3, zorder=3,
                   label=f"RA → LA")
        # highlight PRELP if in panel
        for ag in ANCHOR_HIGHLIGHT:
            if ag in sub.index:
                xa, ya = sub.loc[ag, ["aCAI_healthy", "aCAI_AF"]]
                if ag in keep_set:
                    ax.scatter([xa], [ya], s=42, facecolor="none",
                               edgecolor="black", linewidth=0.9, zorder=4)
                    ax.annotate(
                        ag, xy=(xa, ya),
                        xytext=(8, -8), textcoords="offset points",
                        fontsize=5.2, fontweight="bold", color="black",
                        arrowprops=dict(arrowstyle="-", color="black",
                                        lw=0.4, shrinkA=2),
                    )
        lim = max(np.percentile(np.abs(x), 99.5),
                  np.percentile(np.abs(y), 99.5))
        ax.axhline(0, color="#888", lw=0.4, zorder=0)
        ax.axvline(0, color="#888", lw=0.4, zorder=0)
        ax.plot([-lim, lim], [-lim, lim], color="#888",
                lw=0.4, ls="--", zorder=0)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xlabel("aCAI healthy")
        if k == 0:
            ax.set_ylabel("aCAI AF")
        n_af = int(sub["n_AF_donors"].iloc[0])
        n_sr = int(sub["n_SR_donors"].iloc[0])
        ax.set_title(f"{c}  (AF n={n_af}, SR n={n_sr})", pad=2, fontsize=6.6)
        if k == 2:
            ax.legend(loc="upper left", frameon=False, fontsize=5)
        strip_spines(ax)

    fig.suptitle(
        "Figure 6c — Framework-discovered chamber-identity reversal genes "
        "(positive gene-level evidence)",
        y=1.00, fontsize=7.4, fontweight="bold",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
