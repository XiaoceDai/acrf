#!/usr/bin/env python3
"""Figure 1f (REDESIGN, scheme 2): family-level cell-compartment attribution."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

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

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
STRONG_Z = 2.0

PROG_TSV = ROOT / "results" / "fig1f_baseline_attribution.tsv"
FAM_TSV  = ROOT / "results" / "fig1e_baseline_axis_families.tsv"
OUT_TSV = ROOT / "results" / "fig1f_family_attribution.tsv"
OUT_PDF = ROOT / "figures" / "fig1f_family_attribution.pdf"

COMPARTMENTS = [
    "atrialCM", "fibroblast", "endothelial", "lymphatic_endothelial",
    "macrophage", "lymphocyte", "mast_cell", "B_cell",
    "pericyte", "vascular_smooth_muscle",
    "adipocyte", "epicardial", "neuronal",
]
COMP_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "lymphatic_endothelial": "Lymph-endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "mast_cell": "Mast",
    "B_cell": "B",
    "pericyte": "Pericyte",
    "vascular_smooth_muscle": "vSMC",
    "adipocyte": "Adipo",
    "epicardial": "Epi",
    "neuronal": "Neuro",
}

FAMILY_PALETTE = [
    "#E6AB02", "#66A61E", "#7570B3", "#D95F02",
    "#1B9E77", "#E7298A", "#A6761D", "#666666",
    "#1F78B4", "#FB9A99", "#33A02C", "#B15928",
]


def short_label(s: str) -> str:
    return short_pathway_label(s, 36)


def build_figure(mat, fam_summary, fam_color_map,
                 mode: str, out_path: Path) -> None:
    K = len(fam_summary)
    n_la = int((fam_summary["direction"] == "LA-biased").sum())
    ncc = len(COMPARTMENTS)
    vmax = float(np.nanmax(np.abs(mat))) or 1e-3
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list(
        "LA_RA", [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)])
    cmap.set_bad("#E5E5E5")

    fig = plt.figure(figsize=(COL_DOUBLE * 1.30, COL_DOUBLE * 0.55))
    gs = fig.add_gridspec(
        1, 5, width_ratios=[0.18, 0.22, 2.5, 4.6, 1.7], wspace=0.05,
        top=0.86, bottom=0.22, left=0.04, right=0.985,
    )

    # col 0: LA / RA strip
    ax_dir = fig.add_subplot(gs[0, 0])
    for i, d in enumerate(fam_summary["direction"]):
        ax_dir.add_patch(mpatches.Rectangle(
            (0, i - 0.5), 1, 1,
            facecolor=LA_COLOR if d == "LA-biased" else RA_COLOR,
            edgecolor="white", linewidth=0.2,
        ))
    ax_dir.text(0.5, n_la / 2 - 0.5, "LA",
                rotation=90, ha="center", va="center",
                color="white", fontweight="bold", fontsize=9.0)
    ax_dir.text(0.5, n_la + (K - n_la) / 2 - 0.5, "RA",
                rotation=90, ha="center", va="center",
                color="white", fontweight="bold", fontsize=9.0)
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(K - 0.5, -0.5)
    ax_dir.axis("off")

    # col 1: family color band with F-id
    ax_fam = fig.add_subplot(gs[0, 1], sharey=ax_dir)
    for i, fid in enumerate(fam_summary["family_id"]):
        ax_fam.add_patch(mpatches.Rectangle(
            (0, i - 0.5), 1, 1,
            facecolor=fam_color_map[fid], edgecolor="white", linewidth=0.2,
        ))
        ax_fam.text(0.5, i, f"F{int(fid)}", ha="center", va="center",
                    fontsize=8.6, color="white", fontweight="bold")
    ax_fam.set_xlim(0, 1); ax_fam.set_ylim(K - 0.5, -0.5)
    ax_fam.axis("off")

    # col 2: representative pathway + n
    ax_lbl = fig.add_subplot(gs[0, 2], sharey=ax_dir)
    ax_lbl.set_xlim(0, 1); ax_lbl.set_ylim(K - 0.5, -0.5)
    ax_lbl.axis("off")
    for i, r in fam_summary.iterrows():
        dcol = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        ax_lbl.text(0.0, i, short_label(r["representative"]),
                    fontsize=8.4, ha="left", va="center",
                    color=dcol, fontweight="bold")
        ax_lbl.text(0.0, i + 0.32,
                    f"n = {int(r['n_programs'])} programs",
                    fontsize=7.0, ha="left", va="center",
                    color="#555", style="italic")

    # col 3: heatmap (mode-dependent)
    ax_hm = fig.add_subplot(gs[0, 3], sharey=ax_dir)
    ax_hm.set_xlim(-0.5, ncc - 0.5)
    ax_hm.set_ylim(K - 0.5, -0.5)
    ax_hm.set_aspect("auto")

    for i in range(K):
        for j in range(ncc):
            v = mat[i, j]
            if np.isnan(v):
                ax_hm.add_patch(mpatches.Rectangle(
                    (j - 0.5, i - 0.5), 1.0, 1.0,
                    facecolor="#E5E5E5",
                    edgecolor="white", linewidth=0.2,
                ))
                continue

            if mode == "size":
                # base white background tile
                ax_hm.add_patch(mpatches.Rectangle(
                    (j - 0.5, i - 0.5), 1.0, 1.0,
                    facecolor="#FAFAFA",
                    edgecolor="white", linewidth=0.2,
                ))
                # size ~ |z|/vmax, capped to 1
                s = max(0.18, min(0.92, abs(v) / vmax))
                ax_hm.add_patch(mpatches.Rectangle(
                    (j - s / 2, i - s / 2), s, s,
                    facecolor=cmap(norm(v)),
                    edgecolor="white", linewidth=0.2,
                ))
            else:
                # plain or dot: full color tile
                ax_hm.add_patch(mpatches.Rectangle(
                    (j - 0.5, i - 0.5), 1.0, 1.0,
                    facecolor=cmap(norm(v)),
                    edgecolor="white", linewidth=0.2,
                ))
                if mode == "dot" and abs(v) >= STRONG_Z:
                    inner = "#1A1A1A" if v > 0 else "#1A1A1A"
                    ax_hm.add_patch(mpatches.Circle(
                        (j, i), 0.18,
                        facecolor=inner, edgecolor="white", linewidth=0.3,
                        zorder=4,
                    ))

    ax_hm.set_yticks([])
    ax_hm.set_xticks(np.arange(ncc))
    ax_hm.set_xticklabels([COMP_SHORT[c] for c in COMPARTMENTS],
                          rotation=30, ha="right", fontsize=7.8)
    ax_hm.set_xlabel("compartment", fontsize=9.0)
    ax_hm.axhline(n_la - 0.5, color="black", lw=0.8)

    # col 4: dominant carrier text
    ax_dom = fig.add_subplot(gs[0, 4], sharey=ax_dir)
    ax_dom.set_xlim(0, 1); ax_dom.set_ylim(K - 0.5, -0.5)
    ax_dom.axis("off")
    ax_dom.set_title("dominant carrier", loc="left", fontsize=8.6,
                     fontweight="bold", pad=2)
    for i in range(K):
        row = mat[i]
        if np.all(np.isnan(row)):
            ax_dom.text(0.0, i, "\u2014", fontsize=8.0, ha="left", va="center",
                        color="#888")
            continue
        j = int(np.nanargmax(row))
        z_val = row[j]
        comp = COMP_SHORT[COMPARTMENTS[j]]
        col = "black" if z_val >= 1.5 else "#666"
        wt = "bold" if z_val >= STRONG_Z else "normal"
        ax_dom.text(0.0, i, comp, fontsize=8.2, ha="left", va="center",
                    color=col, fontweight=wt)
        ax_dom.text(1.0, i, f"z={z_val:+.1f}", fontsize=7.6,
                    ha="right", va="center", color="#666")

    # title and shared colorbar
    title_suffix = {
        "plain": "",
        "dot":   " (strong cells: dot)",
        "size":  " (cell size encodes |z|)",
    }[mode]
    ax_hm.set_title("f  Family-level cell compartment attribution (baseline)"
                    + title_suffix,
                    loc="left", fontsize=10.0, fontweight="bold", pad=6)

    cax = fig.add_axes([0.34, 0.06, 0.22, 0.020])
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label("row z-score (mean across family programs)", fontsize=7.8)
    cbar.ax.tick_params(labelsize=7.2)
    cbar.ax.text(0.02, 0.5, "depleted \u2190", transform=cbar.ax.transAxes,
                 ha="left", va="center", fontsize=7.4,
                 color=RA_COLOR, fontweight="bold")
    cbar.ax.text(0.98, 0.5, "\u2192 enriched", transform=cbar.ax.transAxes,
                 ha="right", va="center", fontsize=7.4,
                 color=LA_COLOR, fontweight="bold")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {out_path.relative_to(ROOT)}")


def main() -> None:
    print("Loading per-program attribution + family table ...")
    prog = pd.read_csv(PROG_TSV, sep="\t")
    fam = pd.read_csv(FAM_TSV, sep="\t")

    fam_summary = (
        fam.assign(abs_mean=fam["mean_aCAI"].abs())
           .sort_values("abs_mean", ascending=False)
           .groupby("family_id", sort=False)
           .agg(direction=("direction", "first"),
                n_programs=("pathway", "size"),
                representative=("pathway", "first"))
           .reset_index()
    )
    fam_summary["dir_rank"] = fam_summary["direction"].map(
        {"LA-biased": 0, "RA-biased": 1})
    fam_summary = fam_summary.sort_values(["dir_rank", "family_id"]).reset_index(drop=True)

    fam_z = (prog.groupby("family_id")[COMPARTMENTS]
                 .mean(numeric_only=True))
    fam_z = fam_z.reindex(fam_summary["family_id"]).reset_index()
    fam_z.merge(fam_summary, on="family_id").to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")

    mat = fam_z[COMPARTMENTS].to_numpy(dtype=float)
    fam_color_map = {fid: FAMILY_PALETTE[i % len(FAMILY_PALETTE)]
                     for i, fid in enumerate(fam_summary["family_id"].tolist())}

    build_figure(mat, fam_summary, fam_color_map, mode="size", out_path=OUT_PDF)


if __name__ == "__main__":
    main()
