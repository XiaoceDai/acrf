#!/usr/bin/env python3
"""
Figure 1f — Program × cell-type attribution at single-nucleus
resolution, across the FULL 13-compartment cell-type taxonomy of
GSE255612 (atrialCM, fibroblast, endothelial, lymphatic_endothelial,
macrophage, lymphocyte, mast_cell, B_cell, pericyte,
vascular_smooth_muscle, adipocyte, epicardial, neuronal).

Two complementary attribution heatmaps:

  Left panel: 10 NCP programs (Fig. 1c, curated primitive) ×
              13 compartments
              Per program, per compartment: mean expression of the
              program's leading-edge genes in that compartment,
              z-scored across compartments.

  Right panel: 12 data-driven modules (Fig. 1d–e, data-driven
               primitive; M0–M11) × 13 compartments
               Same per-compartment z-score, evaluated on the gene
               members of each module.

Validation predictions (text-locked before plotting):
  - Adipose module M7 (ADIPOQ, LEP, PLIN1, PCK1, CIDEC) should peak
    in the adipocyte compartment.
  - Mesothelial / PITX2 module M11 (PITX2, ITLN1, MSLN, KRT19, PRG4,
    SLPI) should peak in the epicardial compartment.
  - Atrial-CM module M8 (HAMP, BMP10) should peak in the atrialCM
    compartment.

Inputs:
  Project_P3/GSE255612/pseudobulk_full/<compartment>/expr.tsv
  bulk/Project_P01/normal_chamber_programs.tsv  (10 NCP programs)
  results/n1_modules_long.tsv                   (gene -> module id)
  results/n1_modules_summary.tsv                (module summary)

Output:
  figures/fig1f_attribution_full.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import COL_DOUBLE, panel_label, set_nature_style, strip_spines

set_nature_style()

PB_DIR  = ROOT / "Project_P3" / "GSE255612" / "pseudobulk_full"
NCP_TSV = ROOT / "bulk" / "Project_P01" / "normal_chamber_programs.tsv"
MOD_LONG = ROOT / "results" / "n1_modules_long.tsv"
MOD_SUMM = ROOT / "results" / "n1_modules_summary.tsv"
OUT      = ROOT / "figures" / "fig1f_attribution_full.pdf"

# Ordering / colour conventions for the 13 compartments.
COMPARTMENTS = [
    "atrialCM", "fibroblast", "endothelial", "lymphatic_endothelial",
    "macrophage", "lymphocyte", "mast_cell", "B_cell",
    "pericyte", "vascular_smooth_muscle",
    "adipocyte", "epicardial", "neuronal",
]

# Highlight modules to call out
HIGHLIGHT_MODULES = {7: "M7 adipose", 11: "M11 PITX2/mesothelial", 8: "M8 BMP10/CM"}


def load_compartment_means() -> pd.DataFrame:
    """Return DataFrame (gene × compartment) of mean expression across all
    pseudobulk samples in each compartment (donors regardless of rhythm)."""
    out = {}
    for c in COMPARTMENTS:
        f = PB_DIR / c / "expr.tsv"
        if not f.exists():
            print(f"  WARN: missing {f}")
            continue
        m = pd.read_csv(f, sep="\t", index_col=0)
        out[c] = m.mean(axis=1)
    df = pd.concat(out, axis=1)
    return df


def parse_ncp_le() -> dict[str, list[str]]:
    progs = pd.read_csv(NCP_TSV, sep="\t")
    progs.columns = [c.strip().strip('"') for c in progs.columns]
    progs["program"] = progs["program"].astype(str).str.strip().str.strip('"')
    progs["leading_edge_genes"] = progs["leading_edge_genes"] \
        .astype(str).str.strip().str.strip('"')
    return {
        row["program"]: [g.strip() for g in row["leading_edge_genes"].split(",")
                         if g.strip()]
        for _, row in progs.iterrows()
    }


def parse_module_genes() -> tuple[dict[int, list[str]], pd.DataFrame]:
    long = pd.read_csv(MOD_LONG, sep="\t")
    summ = pd.read_csv(MOD_SUMM, sep="\t").set_index("module_id")
    out = {}
    for mid, sub in long.groupby("module_id"):
        out[int(mid)] = sub["gene"].tolist()
    return out, summ


def attribution_zmatrix(gene_sets: dict, expr_means: pd.DataFrame,
                        max_genes: int | None = None) -> pd.DataFrame:
    """Per gene-set per compartment: mean expression of gene-set members
    in that compartment.  Then z-score per row across compartments so a
    cell's value indicates relative attribution to that compartment."""
    rows = []
    keys = []
    for k, genes in gene_sets.items():
        if max_genes is not None:
            genes = genes[:max_genes]
        present = [g for g in genes if g in expr_means.index]
        if len(present) < 3:
            continue
        sub = expr_means.loc[present]
        means = sub.mean(axis=0)  # mean across genes per compartment
        rows.append(means.values)
        keys.append(k)
    arr = np.array(rows)
    # row-wise z-score
    z = (arr - arr.mean(axis=1, keepdims=True)) \
        / (arr.std(axis=1, keepdims=True) + 1e-9)
    return pd.DataFrame(z, index=keys, columns=expr_means.columns)


def main() -> None:
    print("loading per-compartment expression means ...")
    expr_means = load_compartment_means()
    print(f"  expr_means shape: {expr_means.shape}")
    expr_means = expr_means[[c for c in COMPARTMENTS if c in expr_means.columns]]

    print("parsing NCP leading-edge gene sets ...")
    ncp_sets = parse_ncp_le()
    # cap LE size for noise control
    ncp_z = attribution_zmatrix(ncp_sets, expr_means, max_genes=50)
    print(f"  NCP attribution: {ncp_z.shape}")

    print("parsing N1 module gene sets ...")
    mod_sets, mod_summary = parse_module_genes()
    mod_z = attribution_zmatrix(mod_sets, expr_means)
    print(f"  Module attribution: {mod_z.shape}")
    # rename module index for display, sorted by module_id
    mod_z = mod_z.loc[sorted(mod_z.index)]
    mod_labels = [
        f"M{int(m)} ({HIGHLIGHT_MODULES[int(m)]})"
        if int(m) in HIGHLIGHT_MODULES
        else f"M{int(m)} (n={int(mod_summary.loc[int(m), 'n_genes'])})"
        for m in mod_z.index
    ]
    mod_z.index = mod_labels

    # ---- figure -----------------------------------------------------------
    fig = plt.figure(figsize=(COL_DOUBLE, 0.75 * COL_DOUBLE))
    gs = fig.add_gridspec(1, 2, wspace=0.10, width_ratios=[1.0, 1.0])

    cmap = "RdBu_r"
    vmax = max(np.nanmax(np.abs(ncp_z.values)),
               np.nanmax(np.abs(mod_z.values)))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # ---- left: NCP attribution
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(ncp_z.values, cmap=cmap, norm=norm, aspect="auto")
    ax.set_yticks(np.arange(len(ncp_z)))
    ax.set_yticklabels(ncp_z.index, fontsize=5.5)
    ax.set_xticks(np.arange(len(ncp_z.columns)))
    ax.set_xticklabels(ncp_z.columns, rotation=45, ha="right", fontsize=5)
    ax.set_title("Curated primitive — NCP programs × compartment",
                 loc="left", fontsize=6.6, pad=4)
    panel_label(ax, "f", x=-0.32, y=1.04)
    # mark dominant compartment per row
    for i in range(len(ncp_z)):
        j = int(np.nanargmax(ncp_z.values[i]))
        ax.scatter([j], [i], s=22, facecolor="none", edgecolor="black",
                   linewidth=0.8, zorder=4)

    # ---- right: N1 module attribution
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(mod_z.values, cmap=cmap, norm=norm, aspect="auto")
    ax2.set_yticks(np.arange(len(mod_z)))
    ax2.set_yticklabels(mod_z.index, fontsize=5.5)
    ax2.set_xticks(np.arange(len(mod_z.columns)))
    ax2.set_xticklabels(mod_z.columns, rotation=45, ha="right", fontsize=5)
    ax2.set_title("Data-driven primitive — N1 modules × compartment",
                  loc="left", fontsize=6.6, pad=4)
    # mark dominant compartment + special star for highlighted predictions
    predictions = {
        "M7 (M7 adipose)":               ("adipocyte", "ok"),
        "M11 (M11 PITX2/mesothelial)":   ("epicardial", "ok"),
        "M8 (M8 BMP10/CM)":              ("atrialCM", "ok"),
    }
    for i, label in enumerate(mod_z.index):
        j = int(np.nanargmax(mod_z.values[i]))
        ax2.scatter([j], [i], s=22, facecolor="none", edgecolor="black",
                    linewidth=0.8, zorder=4)
        if label in predictions:
            pred_comp, marker = predictions[label]
            if pred_comp in mod_z.columns:
                pj = list(mod_z.columns).index(pred_comp)
                if pj == j:
                    # prediction matches data
                    ax2.add_patch(mpatches.Rectangle(
                        (j - 0.42, i - 0.42), 0.84, 0.84,
                        facecolor="none", edgecolor="#1B7837",
                        linewidth=1.5, zorder=5,
                    ))
                else:
                    ax2.text(pj, i, "?", fontsize=7, color="#444",
                             ha="center", va="center", fontweight="bold")

    # shared colorbar
    cax = fig.add_axes([0.93, 0.20, 0.012, 0.55])
    cbar = fig.colorbar(im2, cax=cax)
    cbar.set_label("row z-score (relative attribution)", fontsize=5)
    cbar.ax.tick_params(labelsize=5)

    fig.suptitle(
        "Figure 1f — Cell-compartment attribution of both discovery "
        "primitives across the full 13-compartment GSE255612 taxonomy "
        "(green box = data-driven prediction confirmed by single-nucleus carrier)",
        y=0.99, fontsize=7.0, fontweight="bold",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")

    # ---- console summary --------------------------------------------------
    print("\n=== NCP programs: dominant carrier per program (z-score) ===")
    for prog in ncp_z.index:
        j = int(np.nanargmax(ncp_z.loc[prog].values))
        z = ncp_z.loc[prog].values[j]
        print(f"  {prog:50s} -> {ncp_z.columns[j]:25s}  z = {z:+.2f}")

    print("\n=== N1 modules: dominant carrier per module (z-score) ===")
    for mod in mod_z.index:
        j = int(np.nanargmax(mod_z.loc[mod].values))
        z = mod_z.loc[mod].values[j]
        print(f"  {mod:50s} -> {mod_z.columns[j]:25s}  z = {z:+.2f}")


if __name__ == "__main__":
    main()
