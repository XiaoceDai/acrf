#!/usr/bin/env python3
"""
Fig S4b v3 - LAA primary rho vs LAPV sensitivity rho concordance scatter.

Top 6 family by |LAA primary rho| are highlighted; identity line (y=x) and
flipped-sign line (y=-x) drawn for visual reference.  Full-16-family points
shown as background.

Used to address the LAPV sensitivity arm: at n=12 LAPV samples, individual
rho estimates are wide; the panel reports concordance count (number of
families whose LAA and LAPV rho share sign) as the only robust statistic.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

CLS_TSV = ROOT / "results" / "fig4_v3_chronicity_class.tsv"
ATLAS_TSV = ROOT / "results" / "fig_integrated_af_rewiring_atlas.tsv"
OUT_PDF = ROOT / "figures" / "figS4b_LAA_LAPV_concordance.pdf"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def main():
    set_style()
    cls = pd.read_csv(CLS_TSV, sep="\t")
    atlas = pd.read_csv(ATLAS_TSV, sep="\t").set_index("family_id")

    cls["abs_primary"] = cls["rho_primary"].abs()
    top6 = cls.nlargest(6, "abs_primary")["family_id"].tolist()

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(150 * MM, 130 * MM))
    gs = fig.add_gridspec(1, 1, left=0.13, right=0.96, top=0.86, bottom=0.30)
    ax = fig.add_subplot(gs[0, 0])

    # All families as light-grey background
    ax.scatter(cls["rho_primary"], cls["rho_lapv"],
               s=22, c="#CCCCCC", edgecolor="none", alpha=0.7, zorder=1)
    # Top 6 highlighted
    for _, r in cls.iterrows():
        fid = int(r["family_id"])
        if fid not in top6:
            continue
        x = float(r["rho_primary"]); y = float(r["rho_lapv"])
        same = (np.sign(x) == np.sign(y)) and x != 0
        col = "#1B5E20" if same else "#9C2A2A"
        ax.scatter([x], [y], s=80, c=col, edgecolor="white", linewidth=0.6, zorder=4)
        short = atlas.loc[fid, "short"] if fid in atlas.index else f"F{fid}"
        ax.text(x + 0.04, y + 0.04, f"F{fid}  {short[:14]}",
                fontsize=6.2, color=col, fontweight="bold", zorder=5)

    # Reference lines
    lim = 0.7
    ax.plot([-lim, lim], [-lim, lim], color="#1A1A1A", lw=0.7, ls="--", alpha=0.6, zorder=2)
    ax.plot([-lim, lim], [lim, -lim], color="#888", lw=0.6, ls=":", alpha=0.5, zorder=2)
    ax.axhline(0, color="#888", lw=0.5, alpha=0.5)
    ax.axvline(0, color="#888", lw=0.5, alpha=0.5)

    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("LAA primary  Spearman \u03c1", fontsize=7.4)
    ax.set_ylabel("LAPV sensitivity  Spearman \u03c1", fontsize=7.4)
    ax.tick_params(labelsize=6.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # Concordance count
    same_n = int(((np.sign(cls["rho_primary"]) == np.sign(cls["rho_lapv"])) & (cls["rho_primary"] != 0)).sum())
    total_n = int((cls["rho_primary"] != 0).sum())
    ax.text(0.04, 0.96,
            f"LAA-LAPV same-sign concordance:\n{same_n} / {total_n} families\n"
            f"Top 6 by |LAA \u03c1|: green = same sign,  red = opposite sign",
            transform=ax.transAxes,
            fontsize=6.0, color="#1A1A1A", ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#F5F5F5", edgecolor="#888", linewidth=0.4))

    fig.suptitle(
        "Fig S4b  LAA vs LAPV duration-trend concordance (GSE41177)",
        x=0.05, y=0.94, ha="left", fontsize=8.6, fontweight="bold",
    )
    footer_lines = [
        "Each dot: one NCP axis family; x = LAA primary Spearman \u03c1, y = LAPV sensitivity Spearman \u03c1 (both n = 12 AF donors).",
        "Solid diagonal y = x marks perfect LAA-LAPV agreement;  dotted antidiagonal y = -x marks perfect inversion.",
        "",
        "At n = 12 the LAPV per-family \u03c1 is highly uncertain (bootstrap CI is wide, often spanning zero); the only",
        "robust statistic at this resolution is the LAA-LAPV same-sign concordance count.  Individual LAPV \u03c1 values",
        "should not be interpreted as independent confirmation of LAA \u03c1, but as a region sensitivity arm only.",
        "",
        "Top 6 by |LAA \u03c1| are highlighted: green dots = LAA and LAPV same direction;  red dots = opposite direction.",
    ]
    y_top = 0.26
    line_h = 0.022
    for k, line in enumerate(footer_lines):
        fig.text(0.05, y_top - line_h * k, line, fontsize=5.7, color="#444444")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
