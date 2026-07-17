#!/usr/bin/env python3
"""Figure 1d — Data-driven discovery primitive of the framework."""
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

LOAD_PATH = ROOT / "results" / "n1_pc_loadings.tsv"
MOD_PATH  = ROOT / "results" / "n1_modules_summary.tsv"
LONG_PATH = ROOT / "results" / "n1_modules_long.tsv"
OUT       = ROOT / "figures" / "fig1d_data_driven_pca.pdf"

# Genes to label by name (one anchor per highlighted module)
GENE_LABELS = {
    "PITX2":  "M11 / PITX2",
    "ITLN1":  "M11",
    "MSLN":   "M11",
    "ADIPOQ": "M7 / adipose",
    "LEP":    "M7",
    "PLIN1":  "M7",
    "BMP10":  "M8",
    "HAMP":   "M8",
}
HIGHLIGHT_MODULES = {11, 7, 8}


def main() -> None:
    loadings = pd.read_csv(LOAD_PATH, sep="\t", index_col="gene")
    modlong  = pd.read_csv(LONG_PATH, sep="\t").set_index("gene")
    summary  = pd.read_csv(MOD_PATH,  sep="\t").set_index("module_id")

    # attach module id from long table
    loadings["module"] = modlong["module_id"].reindex(loadings.index)

    # rather than reread, we hard-encode from the reported run:
    pc_var = [36.5, 6.6, 5.5, 5.0, 3.7, 3.1, 2.6, 2.3, 2.0, 1.8]
    pc_var_cum = np.cumsum(pc_var)

    fig = plt.figure(figsize=(COL_SINGLE * 1.10, COL_SINGLE * 0.95))
    ax = fig.add_axes([0.16, 0.18, 0.78, 0.74])

    cmap = plt.colormaps.get_cmap("tab20").resampled(12)

    # Background : non-highlighted modules in muted greys
    bg_mask = ~loadings["module"].isin(HIGHLIGHT_MODULES)
    ax.scatter(loadings.loc[bg_mask, "PC1"], loadings.loc[bg_mask, "PC2"],
               s=1.6, c="#CCCCCC", alpha=0.55, edgecolor="none", zorder=1)

    # Foreground : highlighted modules
    mod_color = {11: "#B2182B", 7: "#1B7837", 8: "#762A83"}
    mod_label = {11: f"M11 (PITX2 collapse, n={int(summary.loc[11, 'n_genes'])})",
                 7:  f"M7 (adipose expansion, n={int(summary.loc[7, 'n_genes'])})",
                 8:  f"M8 (atrial CM, n={int(summary.loc[8, 'n_genes'])})"}
    for m in HIGHLIGHT_MODULES:
        idx = loadings["module"] == m
        ax.scatter(loadings.loc[idx, "PC1"], loadings.loc[idx, "PC2"],
                   s=10, color=mod_color[m], alpha=0.9,
                   edgecolor="white", linewidth=0.25, zorder=3,
                   label=mod_label[m])

    # Anchor gene labels
    for g, _ in GENE_LABELS.items():
        if g not in loadings.index:
            continue
        x, y = loadings.loc[g, "PC1"], loadings.loc[g, "PC2"]
        ax.scatter([x], [y], s=24, facecolor="none",
                   edgecolor="black", linewidth=0.6, zorder=4)
        ax.annotate(
            g, xy=(x, y), xytext=(7, 5), textcoords="offset points",
            fontsize=5.4, fontweight="bold", color="black",
            arrowprops=dict(arrowstyle="-", color="black", lw=0.4, shrinkA=2),
        )

    ax.axhline(0, color="#888", lw=0.4, zorder=0)
    ax.axvline(0, color="#888", lw=0.4, zorder=0)
    ax.set_xlabel(f"PC1 ({pc_var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({pc_var[1]:.1f}%)")
    strip_spines(ax)
    ax.legend(loc="lower left", frameon=False, fontsize=5)

    ax_in = fig.add_axes([0.65, 0.66, 0.27, 0.20])
    x = np.arange(1, len(pc_var) + 1)
    ax_in.bar(x, pc_var, color="#444", edgecolor="white", linewidth=0.3)
    ax_in.plot(x, pc_var_cum, color="#B2182B", lw=0.8, marker="o",
               markersize=2.5, markerfacecolor="#B2182B",
               markeredgecolor="white", markeredgewidth=0.3)
    ax_in.set_xlim(0.4, 10.6); ax_in.set_ylim(0, 100)
    ax_in.set_xticks([1, 5, 10])
    ax_in.set_yticks([0, 50, 100])
    ax_in.tick_params(labelsize=4.5, length=1.5, pad=1)
    ax_in.set_xlabel("PC", fontsize=5, labelpad=1)
    ax_in.set_ylabel("var % (cum.)", fontsize=5, labelpad=1)
    ax_in.text(10.4, pc_var_cum[-1] + 2, f"{pc_var_cum[-1]:.0f}%",
               fontsize=4.5, color="#B2182B", ha="right", va="bottom")
    for s in ("top", "right"):
        ax_in.spines[s].set_visible(False)

    fig.text(0.04, 0.96, "d", fontsize=10, fontweight="bold",
             ha="left", va="top")
    fig.text(0.16, 0.96,
             "Data-driven primitive: PCA on the chamber asymmetry "
             "matrix M_aCAI (3,781 genes × 59 donors)",
             fontsize=6.6, fontweight="bold", ha="left", va="top")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
