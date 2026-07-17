#!/usr/bin/env python3
"""
Figure 1F - 30-program healthy LA/RA cell-type carrier atlas.

Design choices:
  - show all 30 selected programs (no extra filtering)
  - family is used as grouping annotation only
  - heatmap tile color encodes direction-aware z-score sign
  - heatmap tile size encodes |z-score| magnitude
"""
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

from _pathway_labels import short_pathway_label
from plot_style import COL_DOUBLE, FG_COLOR, LA_COLOR, LINE_COLOR, MUTED_COLOR, RA_COLOR, set_nature_style

set_nature_style()
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

META_IN = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
ATTR_IN = ROOT / "results" / "fig1f_v3_baseline_attribution.tsv"
OUT_PDF = ROOT / "figures" / "fig1f_deployed_panel_map.pdf"
OUT_ROWMAP = ROOT / "results" / "fig1f_v3_program_celltype_rowmap.tsv"

COMPARTMENTS = [
    "atrialCM",
    "fibroblast",
    "endothelial",
    "lymphatic_endothelial",
    "macrophage",
    "lymphocyte",
    "mast_cell",
    "B_cell",
    "pericyte",
    "vascular_smooth_muscle",
    "adipocyte",
    "epicardial",
    "neuronal",
]
COMP_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "lymphatic_endothelial": "L-Endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "mast_cell": "Mast",
    "B_cell": "B",
    "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC",
    "adipocyte": "Adipo",
    "epicardial": "Epi",
    "neuronal": "Neuro",
}
FAMILY_PALETTE = [
    "#E6AB02",
    "#66A61E",
    "#7570B3",
    "#D95F02",
    "#1B9E77",
    "#E7298A",
    "#A6761D",
    "#666666",
    "#1F78B4",
    "#FB9A99",
    "#33A02C",
    "#B15928",
]


def build_family_spans(df: pd.DataFrame) -> list[tuple[int, int, int]]:
    spans: list[tuple[int, int, int]] = []
    if df.empty:
        return spans
    start = 0
    current = int(df.loc[0, "family_id"])
    for i in range(1, len(df)):
        fid = int(df.loc[i, "family_id"])
        if fid != current:
            spans.append((start, i - 1, current))
            start = i
            current = fid
    spans.append((start, len(df) - 1, current))
    return spans


def main() -> None:
    meta = pd.read_csv(META_IN, sep="\t")
    attr = pd.read_csv(ATTR_IN, sep="\t")

    comp_cols = [c for c in COMPARTMENTS if c in attr.columns]
    if not comp_cols:
        raise ValueError("No compartment z-score columns found in fig1f_v3_baseline_attribution.tsv.")

    keep_meta = ["program_order", "pathway", "direction", "family_id"]
    keep_meta = [c for c in keep_meta if c in meta.columns]
    df = meta[keep_meta].merge(
        attr[["pathway", "n_genes_used"] + comp_cols],
        on="pathway",
        how="left",
    )

    df["dir_rank"] = df["direction"].map({"LA-biased": 0, "RA-biased": 1}).fillna(2)
    df = df.sort_values(["dir_rank", "family_id", "program_order"]).reset_index(drop=True)
    df["row_id"] = np.arange(1, len(df) + 1)
    df.to_csv(OUT_ROWMAP, sep="\t", index=False)

    n = len(df)
    n_la = int((df["direction"] == "LA-biased").sum())
    family_ids = df["family_id"].tolist()
    fam_color_map = {}
    for i, fid in enumerate(sorted(set(int(x) for x in family_ids))):
        fam_color_map[fid] = FAMILY_PALETTE[i % len(FAMILY_PALETTE)]

    mat = df[comp_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    vmax = float(np.nanmax(np.abs(mat))) if np.isfinite(mat).any() else 1.0
    vmax = max(vmax, 1e-6)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list("LA_RA", [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)])

    fig = plt.figure(figsize=(COL_DOUBLE * 1.34, COL_DOUBLE * 0.86))
    gs = fig.add_gridspec(
        1,
        5,
        width_ratios=[0.16, 0.36, 3.10, 5.95, 1.75],
        wspace=0.03,
        top=0.90,
        bottom=0.16,
        left=0.03,
        right=0.995,
    )

    ax_dir = fig.add_subplot(gs[0, 0])
    for i, d in enumerate(df["direction"]):
        ax_dir.add_patch(
            mpatches.Rectangle(
                (0, i - 0.5),
                1,
                1,
                facecolor=LA_COLOR if d == "LA-biased" else RA_COLOR,
                edgecolor="white",
                linewidth=0.25,
            )
        )
    ax_dir.text(0.5, (n_la - 1) / 2 if n_la > 0 else 0.0, "LA", color="white", fontsize=8.6, fontweight="bold", ha="center", va="center", rotation=90)
    if n - n_la > 0:
        ax_dir.text(0.5, n_la + (n - n_la - 1) / 2, "RA", color="white", fontsize=8.6, fontweight="bold", ha="center", va="center", rotation=90)
    ax_dir.set_xlim(0, 1)
    ax_dir.set_ylim(n - 0.5, -0.5)
    ax_dir.axis("off")

    ax_family = fig.add_subplot(gs[0, 1], sharey=ax_dir)
    ax_family.set_xlim(0, 1)
    ax_family.set_ylim(n - 0.5, -0.5)
    ax_family.axis("off")
    for y0, y1, fid in build_family_spans(df):
        ax_family.add_patch(
            mpatches.Rectangle(
                (0, y0 - 0.5),
                1,
                (y1 - y0 + 1),
                facecolor=fam_color_map[int(fid)],
                edgecolor="white",
                linewidth=0.25,
            )
        )
        ax_family.text(0.5, (y0 + y1) / 2, f"F{int(fid)}", color="white", fontsize=8.0, fontweight="bold", ha="center", va="center")

    ax_lbl = fig.add_subplot(gs[0, 2], sharey=ax_dir)
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(n - 0.5, -0.5)
    ax_lbl.axis("off")
    for i, row in df.iterrows():
        col = LA_COLOR if row["direction"] == "LA-biased" else RA_COLOR
        ax_lbl.text(
            0.0,
            i,
            short_pathway_label(str(row["pathway"]), 42),
            fontsize=7.7,
            color=col,
            fontweight="bold",
            ha="left",
            va="center",
        )

    ax_hm = fig.add_subplot(gs[0, 3], sharey=ax_dir)
    ax_hm.set_xlim(-0.5, len(comp_cols) - 0.5)
    ax_hm.set_ylim(n - 0.5, -0.5)
    ax_hm.set_xticks(np.arange(len(comp_cols)))
    ax_hm.set_xticklabels([COMP_SHORT.get(c, c) for c in comp_cols], rotation=30, ha="right", fontsize=7.2)
    ax_hm.set_yticks([])
    ax_hm.set_xlabel("cell compartment (tile color = signed z, tile size = |z|)", fontsize=8.0)

    # Draw size-encoded heatmap tiles.
    for i in range(n):
        for j in range(len(comp_cols)):
            ax_hm.add_patch(
                mpatches.Rectangle(
                    (j - 0.5, i - 0.5),
                    1.0,
                    1.0,
                    facecolor="#F9F9F9",
                    edgecolor="white",
                    linewidth=0.2,
                )
            )
            v = mat[i, j]
            if not np.isfinite(v):
                continue
            scale = min(abs(v) / vmax, 1.0)
            side = 0.16 + 0.78 * scale
            ax_hm.add_patch(
                mpatches.Rectangle(
                    (j - side / 2, i - side / 2),
                    side,
                    side,
                    facecolor=cmap(norm(v)),
                    edgecolor="white",
                    linewidth=0.25,
                )
            )

    # Draw group separators.
    if n_la > 0 and n_la < n:
        for ax in (ax_lbl, ax_hm):
            ax.axhline(n_la - 0.5, color=LINE_COLOR, lw=0.85)
    for y0, y1, _ in build_family_spans(df):
        if y1 < n - 1:
            ysep = y1 + 0.5
            for ax in (ax_lbl, ax_hm):
                ax.axhline(ysep, color="#E5E5E5", lw=0.6)

    ax_dom = fig.add_subplot(gs[0, 4], sharey=ax_dir)
    ax_dom.set_xlim(0, 1)
    ax_dom.set_ylim(n - 0.5, -0.5)
    ax_dom.axis("off")
    ax_dom.set_title("dominant carrier", fontsize=8.2, fontweight="bold", loc="left", pad=3)
    for i in range(n):
        row = mat[i, :]
        if not np.isfinite(row).any():
            ax_dom.text(0.0, i, "-", fontsize=7.4, color=MUTED_COLOR, ha="left", va="center")
            continue
        j = int(np.nanargmax(row))
        comp = COMP_SHORT.get(comp_cols[j], comp_cols[j])
        val = float(row[j])
        ax_dom.text(
            0.0,
            i,
            comp,
            fontsize=7.9,
            fontweight="bold" if abs(val) >= 1.5 else "normal",
            color=FG_COLOR,
            ha="left",
            va="center",
        )

    ax_hm.set_title(
        "Fig 1e  NCP deployment map (30 programs, 16 families; HCA LA/RA carriers, family-grouped)",
        fontsize=10.0,
        fontweight="bold",
        loc="left",
        pad=6,
    )

    cax = fig.add_axes([0.405, 0.07, 0.22, 0.018])
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("family/program cell-type z-score", fontsize=7.1)
    cb.ax.tick_params(labelsize=6.8)
    cb.ax.text(0.02, 0.5, "RA-like", transform=cb.ax.transAxes, color=RA_COLOR, fontsize=6.8, fontweight="bold", ha="left", va="center")
    cb.ax.text(0.98, 0.5, "LA-like", transform=cb.ax.transAxes, color=LA_COLOR, fontsize=6.8, fontweight="bold", ha="right", va="center")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)

    print(f"wrote {OUT_ROWMAP.relative_to(ROOT)}")
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
