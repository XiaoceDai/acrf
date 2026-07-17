#!/usr/bin/env python3
"""
Fig 1F (v3, HCA) — chamber-resolved baseline cell-compartment attribution.

Visual style mirrors fig1f_deployed_panel_map.py (size-encoded tile heatmap,
family color band, dominant-carrier text column), with the key extension
that LA and RA arms each get their own heatmap and their own dominant-carrier
column.  The consensus class strip on the right tells the reader at a glance
whether the LA and RA arms agree on the dominant carrier.
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label  # noqa: E402
from plot_style import (  # noqa: E402
    COL_DOUBLE,
    FG_COLOR,
    LA_COLOR,
    LINE_COLOR,
    MUTED_COLOR,
    RA_COLOR,
    set_nature_style,
    enforce_font_floor, FS_PANEL_TITLE, FS_TICK, FS_AXIS_LABEL, FS_LEGEND,
)

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
ATTR_IN = ROOT / "results" / "fig1f_v3_HCA_baseline_attribution.tsv"
OUT_PDF = ROOT / "figures" / "fig1f_HCA_baseline_attribution.pdf"
OUT_ROWMAP = ROOT / "results" / "fig1f_v3_HCA_program_celltype_rowmap.tsv"

COMPARTMENTS = [
    "atrialCM",
    "fibroblast",
    "endothelial",
    "macrophage",
    "lymphocyte",
    "pericyte",
    "vascular_smooth_muscle",
    "adipocyte",
    "neuronal",
]
COMP_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC",
    "adipocyte": "Adipo",
    "neuronal": "Neuro",
}
FAMILY_PALETTE = [
    "#E6AB02", "#66A61E", "#7570B3", "#D95F02",
    "#1B9E77", "#E7298A", "#A6761D", "#666666",
    "#1F78B4", "#FB9A99", "#33A02C", "#B15928",
]

MIXED_LABEL = "Mixed/Unresolved"
MIXED_COLOR = "#BDBDBD"
CONSENSUS_COLORS = {
    "SHARED": "#1B5E20",
    "SPLIT": "#B07AA1",
    "LA_ONLY": LA_COLOR,
    "RA_ONLY": RA_COLOR,
    "MIXED": "#9C9C9C",
}
CONSENSUS_SHORT = {"SHARED": "S", "SPLIT": "X", "LA_ONLY": "L", "RA_ONLY": "R", "MIXED": "M"}


def build_family_spans(df: pd.DataFrame) -> list[tuple[int, int, int]]:
    spans: list[tuple[int, int, int]] = []
    if df.empty:
        return spans
    start = 0
    current = int(df.iloc[0]["family_id"])
    for i in range(1, len(df)):
        fid = int(df.iloc[i]["family_id"])
        if fid != current:
            spans.append((start, i - 1, current))
            start = i
            current = fid
    spans.append((start, len(df) - 1, current))
    return spans


def draw_size_tile_heatmap(
    ax: plt.Axes,
    mat: np.ndarray,
    n: int,
    n_cols: int,
    norm: TwoSlopeNorm,
    cmap,
    vmax: float,
) -> None:
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    for i in range(n):
        for j in range(n_cols):
            ax.add_patch(
                mpatches.Rectangle(
                    (j - 0.5, i - 0.5), 1.0, 1.0,
                    facecolor="#F9F9F9", edgecolor="white", linewidth=0.2,
                )
            )
            v = mat[i, j]
            if not np.isfinite(v):
                continue
            scale = min(abs(v) / vmax, 1.0)
            side = 0.16 + 0.78 * scale
            ax.add_patch(
                mpatches.Rectangle(
                    (j - side / 2, i - side / 2), side, side,
                    facecolor=cmap(norm(v)), edgecolor="white", linewidth=0.25,
                )
            )


def main() -> None:
    if not ATTR_IN.exists():
        sys.exit(f"missing {ATTR_IN}")
    if not META_IN.exists():
        sys.exit(f"missing {META_IN}")

    meta = pd.read_csv(META_IN, sep="\t")
    attr = pd.read_csv(ATTR_IN, sep="\t")

    la_cols = [f"LA_{c}" for c in COMPARTMENTS if f"LA_{c}" in attr.columns]
    ra_cols = [f"RA_{c}" for c in COMPARTMENTS if f"RA_{c}" in attr.columns]
    if not la_cols or not ra_cols:
        sys.exit("attribution TSV missing LA_/RA_ columns")
    keep_meta = [c for c in ("program_order", "pathway", "direction", "family_id") if c in meta.columns]

    keep_attr = [
        "pathway",
        "n_genes_used_LA",
        "n_genes_used_RA",
        "LA_carrier",
        "RA_carrier",
        "consensus_class",
    ] + la_cols + ra_cols
    keep_attr = [c for c in keep_attr if c in attr.columns]
    df = meta[keep_meta].merge(attr[keep_attr], on="pathway", how="left")

    df["dir_rank"] = df["direction"].map({"LA-biased": 0, "RA-biased": 1}).fillna(2)
    df = df.sort_values(["dir_rank", "family_id", "program_order"]).reset_index(drop=True)
    df["row_id"] = np.arange(1, len(df) + 1)
    OUT_ROWMAP.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_ROWMAP, sep="\t", index=False)

    n = len(df)
    n_la_dir = int((df["direction"] == "LA-biased").sum())

    family_ids = df["family_id"].tolist()
    fam_color_map: dict[int, str] = {}
    for i, fid in enumerate(sorted({int(x) for x in family_ids})):
        fam_color_map[fid] = FAMILY_PALETTE[i % len(FAMILY_PALETTE)]

    la_mat = df[la_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    ra_mat = df[ra_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    finite = np.concatenate([la_mat[np.isfinite(la_mat)], ra_mat[np.isfinite(ra_mat)]])
    vmax = float(np.nanmax(np.abs(finite))) if finite.size else 1.0
    vmax = max(vmax, 1e-6)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list("LA_RA", [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)])

    n_cols_comp = len(COMPARTMENTS)
    fig = plt.figure(figsize=(COL_DOUBLE * 1.66, COL_DOUBLE * 1.02))  # taller for larger fonts
    gs = fig.add_gridspec(
        1,
        8,
        width_ratios=[0.18, 0.38, 3.10, 4.40, 4.40, 1.20, 1.20, 0.40],
        wspace=0.04,
        top=0.88,
        bottom=0.14,
        left=0.025,
        right=0.997,
    )

    ax_dir = fig.add_subplot(gs[0, 0])
    ax_family = fig.add_subplot(gs[0, 1], sharey=ax_dir)
    ax_lbl = fig.add_subplot(gs[0, 2], sharey=ax_dir)
    ax_la = fig.add_subplot(gs[0, 3], sharey=ax_dir)
    ax_ra = fig.add_subplot(gs[0, 4], sharey=ax_dir)
    ax_dom_la = fig.add_subplot(gs[0, 5], sharey=ax_dir)
    ax_dom_ra = fig.add_subplot(gs[0, 6], sharey=ax_dir)
    ax_cls = fig.add_subplot(gs[0, 7], sharey=ax_dir)

    # Direction strip
    for i, d in enumerate(df["direction"]):
        ax_dir.add_patch(
            mpatches.Rectangle(
                (0, i - 0.5), 1, 1,
                facecolor=LA_COLOR if d == "LA-biased" else RA_COLOR,
                edgecolor="white", linewidth=0.25,
            )
        )
    ax_dir.text(0.5, (n_la_dir - 1) / 2 if n_la_dir > 0 else 0.0, "LA",
                color="white", fontsize=enforce_font_floor(FS_PANEL_TITLE + 0.5), fontweight="bold", ha="center", va="center", rotation=90)
    if n - n_la_dir > 0:
        ax_dir.text(0.5, n_la_dir + (n - n_la_dir - 1) / 2, "RA",
                    color="white", fontsize=enforce_font_floor(FS_PANEL_TITLE + 0.5), fontweight="bold", ha="center", va="center", rotation=90)
    ax_dir.set_xlim(0, 1)
    ax_dir.set_ylim(n - 0.5, -0.5)
    ax_dir.axis("off")

    # Family band
    ax_family.set_xlim(0, 1)
    ax_family.set_ylim(n - 0.5, -0.5)
    ax_family.axis("off")
    for y0, y1, fid in build_family_spans(df):
        ax_family.add_patch(
            mpatches.Rectangle(
                (0, y0 - 0.5), 1, (y1 - y0 + 1),
                facecolor=fam_color_map[int(fid)], edgecolor="white", linewidth=0.25,
            )
        )
        ax_family.text(0.5, (y0 + y1) / 2, f"F{int(fid)}",
                       color="white", fontsize=enforce_font_floor(FS_PANEL_TITLE), fontweight="bold", ha="center", va="center")

    # Pathway labels
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(n - 0.5, -0.5)
    ax_lbl.axis("off")
    for i, row in df.iterrows():
        col = LA_COLOR if row["direction"] == "LA-biased" else RA_COLOR
        ax_lbl.text(
            0.0, i,
            short_pathway_label(str(row["pathway"]), 36),
            fontsize=enforce_font_floor(FS_TICK + 0.5), color=col, fontweight="bold",
            ha="left", va="center",
        )

    # LA + RA size-encoded heatmaps
    draw_size_tile_heatmap(ax_la, la_mat, n, n_cols_comp, norm, cmap, vmax)
    draw_size_tile_heatmap(ax_ra, ra_mat, n, n_cols_comp, norm, cmap, vmax)

    for ax, ttl, color in [(ax_la, "LA  baseline z-score", LA_COLOR),
                           (ax_ra, "RA  baseline z-score", RA_COLOR)]:
        ax.set_xticks(np.arange(n_cols_comp))
        ax.set_xticklabels([COMP_SHORT[c] for c in COMPARTMENTS], rotation=30, ha="right", fontsize=enforce_font_floor(FS_TICK - 0.5))
        ax.set_yticks([])
        ax.set_title(ttl, fontsize=enforce_font_floor(FS_PANEL_TITLE + 0.5), fontweight="bold", color=color, pad=4, loc="left")

    ax_la.set_xlabel("cell compartment  (tile size = |z|, color = sign)", fontsize=enforce_font_floor(FS_AXIS_LABEL))
    ax_ra.set_xlabel("cell compartment  (tile size = |z|, color = sign)", fontsize=enforce_font_floor(FS_AXIS_LABEL))

    # Group separators (LA/RA direction band + family bands)
    if 0 < n_la_dir < n:
        for ax in (ax_lbl, ax_la, ax_ra, ax_dom_la, ax_dom_ra, ax_cls):
            ax.axhline(n_la_dir - 0.5, color=LINE_COLOR, lw=0.85)
    for y0, y1, _ in build_family_spans(df):
        if y1 < n - 1:
            ysep = y1 + 0.5
            for ax in (ax_lbl, ax_la, ax_ra, ax_dom_la, ax_dom_ra, ax_cls):
                ax.axhline(ysep, color="#E5E5E5", lw=0.6)

    # Dominant carrier text columns (LA / RA)
    for ax, col_key, color, title in [
        (ax_dom_la, "LA_carrier", LA_COLOR, "LA carrier"),
        (ax_dom_ra, "RA_carrier", RA_COLOR, "RA carrier"),
    ]:
        ax.set_xlim(0, 1)
        ax.set_ylim(n - 0.5, -0.5)
        ax.axis("off")
        ax.set_title(title, fontsize=enforce_font_floor(FS_PANEL_TITLE), fontweight="bold", color=color, loc="left", pad=3)
        for i, row in df.iterrows():
            carrier = str(row[col_key])
            if carrier == MIXED_LABEL:
                ax.text(0.0, i, "Mix", fontsize=enforce_font_floor(FS_TICK + 0.5), color=MUTED_COLOR, style="italic", ha="left", va="center")
            else:
                ax.text(
                    0.0, i,
                    COMP_SHORT.get(carrier, carrier),
                    fontsize=enforce_font_floor(FS_TICK + 0.5), color=FG_COLOR, fontweight="bold",
                    ha="left", va="center",
                )

    # Consensus class strip
    ax_cls.set_xlim(0, 1)
    ax_cls.set_ylim(n - 0.5, -0.5)
    ax_cls.axis("off")
    ax_cls.set_title("class", fontsize=enforce_font_floor(FS_PANEL_TITLE), fontweight="bold", color=FG_COLOR, loc="left", pad=3)
    for i, row in df.iterrows():
        cls = str(row["consensus_class"])
        ax_cls.add_patch(
            mpatches.Rectangle(
                (0, i - 0.40), 1.0, 0.80,
                facecolor=CONSENSUS_COLORS.get(cls, "#9C9C9C"), edgecolor="white", linewidth=0.25,
            )
        )
        ax_cls.text(0.5, i, CONSENSUS_SHORT.get(cls, "?"), color="white",
                    fontsize=enforce_font_floor(FS_TICK + 0.5), fontweight="bold", ha="center", va="center")

    # Header / colorbar
    ax_la.text(
        -0.5, -1.6,
        "Fig 1e  NCP deployment map (30 programs, 16 families; HCA LA/RA carriers, family-grouped)",
        fontsize=enforce_font_floor(FS_PANEL_TITLE + 2), fontweight="bold", ha="left",
    )

    cax = fig.add_axes([0.395, 0.045, 0.22, 0.018])
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("baseline expression z-score (across compartments)", fontsize=enforce_font_floor(FS_AXIS_LABEL))
    cb.ax.tick_params(labelsize=enforce_font_floor(FS_TICK - 0.5))
    cb.ax.text(0.02, 0.5, "RA-like", transform=cb.ax.transAxes, color=RA_COLOR, fontsize=enforce_font_floor(FS_LEGEND), fontweight="bold", ha="left", va="center")
    cb.ax.text(0.98, 0.5, "LA-like", transform=cb.ax.transAxes, color=LA_COLOR, fontsize=enforce_font_floor(FS_LEGEND), fontweight="bold", ha="right", va="center")

    # Class legend (right of colorbar) - only show classes that actually
    # occur in this panel; SPLIT/MIXED tend to be empty in many runs.
    present_classes = list(dict.fromkeys(df["consensus_class"].astype(str).tolist()))
    full_items = [
        ("SHARED", "S = LA & RA same"),
        ("LA_ONLY", "L = LA only"),
        ("RA_ONLY", "R = RA only"),
        ("SPLIT", "X = LA \u2260 RA"),
        ("MIXED", "M = neither"),
    ]
    items = [(lbl, CONSENSUS_COLORS[k]) for k, lbl in full_items if k in present_classes]

    leg_x = 0.66
    fig.text(leg_x, 0.075, "carrier class:", fontsize=enforce_font_floor(FS_AXIS_LABEL), fontweight="bold", color=FG_COLOR)
    x0 = leg_x + 0.06
    for label, color in items:
        fig.patches.append(mpatches.Rectangle((x0, 0.058), 0.011, 0.018, transform=fig.transFigure, facecolor=color, edgecolor="none"))
        fig.text(x0 + 0.015, 0.067, label, fontsize=enforce_font_floor(FS_TICK - 0.5), color=FG_COLOR, va="center")
        x0 += 0.072

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)

    print(f"wrote {OUT_ROWMAP.relative_to(ROOT)}")
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
