#!/usr/bin/env python3
"""
Figure 1e — Data-driven module landscape: healthy chamber bias vs AF
rewiring magnitude for the 12 modules identified in Fig. 1d.

Each module is a bubble located by:
    x = mean aCAI in healthy donors  (signed: + LA, - RA)
    y = mean ΔaCAI(AF − SR)          (+ LA-side expansion in AF,
                                      - LA-side collapse in AF)
    bubble area = number of genes in the module
    bubble colour = annotation class (canonical / novel; here all
                    are novel given Jaccard < 0.30 against any NCP
                    panel program)

Three modules are explicitly labelled:
    M11 (PITX2-anchored LA-identity collapse), bottom-right quadrant
    M7  (adipose program, LA expansion),        top-right quadrant
    M8  (BMP10/HAMP atrial CM, partial RA→LA),  top-left quadrant

Inputs:
    results/n1_modules_summary.tsv
Output:
    figures/fig1e_data_driven_modules.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import COL_SINGLE, panel_label, set_nature_style, strip_spines

set_nature_style()

IN  = ROOT / "results" / "n1_modules_summary.tsv"
OUT = ROOT / "figures" / "fig1e_data_driven_modules.pdf"


def main() -> None:
    df = pd.read_csv(IN, sep="\t")
    fig = plt.figure(figsize=(COL_SINGLE * 1.10, COL_SINGLE * 0.95))
    ax = fig.add_axes([0.18, 0.18, 0.76, 0.74])

    # bubbles
    sizes = (df["n_genes"] / 5.0).clip(8, 220)  # area scaling
    colors = ["#B2182B" if r["module_id"] == 11
              else "#1B7837" if r["module_id"] == 7
              else "#762A83" if r["module_id"] == 8
              else "#999999" for _, r in df.iterrows()]
    edge   = ["black" if r["module_id"] in (11, 7, 8) else "white"
              for _, r in df.iterrows()]
    ax.scatter(df["healthy_aCAI_mean"], df["delta_AF_minus_SR"],
               s=sizes, c=colors, edgecolor=edge, linewidth=0.5,
               alpha=0.85, zorder=3)

    # quadrant guides
    ax.axhline(0, color="#888", lw=0.5, zorder=1)
    ax.axvline(0, color="#888", lw=0.5, zorder=1)

    # quadrant labels
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    pad = 0.02
    ax.text(xmax * 0.96, ymax * 0.96, "LA-biased  AF expand",
            ha="right", va="top", fontsize=4.8, color="#666",
            style="italic")
    ax.text(xmax * 0.96, ymin * 0.96, "LA-biased  AF collapse",
            ha="right", va="bottom", fontsize=4.8, color="#666",
            style="italic")
    ax.text(xmin * 0.96, ymax * 0.96, "RA-biased  AF expand",
            ha="left", va="top", fontsize=4.8, color="#666",
            style="italic")
    ax.text(xmin * 0.96, ymin * 0.96, "RA-biased  AF collapse",
            ha="left", va="bottom", fontsize=4.8, color="#666",
            style="italic")

    # explicit labels for highlighted modules
    label_offsets = {11: (-2.0, -0.06), 7: (1.4, 0.06), 8: (1.5, 0.04)}
    for mod_id, (xoff, yoff) in label_offsets.items():
        row = df[df["module_id"] == mod_id]
        if row.empty:
            continue
        x = row["healthy_aCAI_mean"].iloc[0]
        y = row["delta_AF_minus_SR"].iloc[0]
        n = int(row["n_genes"].iloc[0])
        if mod_id == 11:
            text = f"M11  PITX2 + ITLN1, MSLN\n(LA-identity collapse; n={n})"
        elif mod_id == 7:
            text = f"M7  ADIPOQ + LEP + PLIN1\n(LA adipose expansion; n={n})"
        else:
            text = f"M8  HAMP + BMP10\n(atrial CM partial; n={n})"
        ax.annotate(
            text, xy=(x, y), xytext=(x + xoff, y + yoff),
            fontsize=5.0, fontweight="bold", color="black",
            ha="left" if xoff > 0 else "right",
            va="center",
            arrowprops=dict(arrowstyle="-", color="black", lw=0.4, shrinkA=4),
        )

    ax.set_xlabel("mean aCAI in healthy  (LA − RA, log2)")
    ax.set_ylabel("Δ aCAI  (AF − SR, log2)")
    strip_spines(ax)

    # bubble-size legend (top-right inset)
    handles = []
    for n in [50, 300, 1000]:
        handles.append(plt.scatter([], [], s=n / 5.0, c="#999999",
                                    edgecolor="white", linewidth=0.3,
                                    label=f"{n} genes"))
    ax.legend(handles=handles, loc="upper left", frameon=False,
              fontsize=5, scatterpoints=1, labelspacing=1.0,
              borderpad=0.4, handletextpad=0.6)

    fig.text(0.04, 0.96, "e", fontsize=10, fontweight="bold",
             ha="left", va="top")
    fig.text(0.18, 0.96,
             "Data-driven module landscape: healthy bias vs AF rewiring "
             "(12 modules; bubble = n genes)",
             fontsize=6.6, fontweight="bold", ha="left", va="top")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
