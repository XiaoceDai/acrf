#!/usr/bin/env python3
"""
Fig 3b — family × compartment Δ-singscore matrix with fidelity class column.

- Fill: AF − Ctrl family Δ singscore
- Black filled dot: nominal family-level p_perm < 0.05 (no p < 0.10 trends on main)
- Black outline: healthy HCA LA carrier cell (Fig 1F prior; interpretation only)
- Hatch: ct_included == False
- Right class column + footer: 1 Fidelity | 2 Relay | 1 Extension | 12 Below-detection
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label  # noqa: E402
from plot_style import (  # noqa: E402
    FG_COLOR,
    FAINT_COLOR,
    LA_DEEP,
    RA_DEEP,
    set_nature_style,
    enforce_font_floor,
    FS_PANEL_TITLE,
    FS_TICK,
    FS_AXIS_LABEL,
    FS_FOOTNOTE,
)

DELTA_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
CLS_TSV = ROOT / "results" / "fig3_v3_carrier_fidelity_class.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
CANONICAL_ORDER_TSV = ROOT / "results" / "ncp_v3_canonical_display_order.tsv"
OUT_PDF = ROOT / "figures" / "fig3_b_family_compartment_heatmap.pdf"

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
COMPARTMENT_PALETTE = {
    "atrialCM": "#4E79A7",
    "fibroblast": "#A0CBE8",
    "endothelial": "#F28E2B",
    "macrophage": "#59A14F",
    "lymphocyte": "#8CD17D",
    "pericyte": "#499894",
    "vascular_smooth_muscle": "#86BCB6",
    "adipocyte": "#E15759",
    "neuronal": "#79706E",
}
COMPARTMENT_SHORT = {
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
MIXED_COLOR = "#BDBDBD"
MIXED_LABEL = "Mixed/Unresolved"

CLASS_ORDER = ["Fidelity", "Relay", "Extension", "Below_detection"]
CLASS_COLORS = {
    "Fidelity": "#7FB069",
    "Relay": "#8E6BAA",
    "Extension": "#5F9EA0",
    "Below_detection": "#A9A9A9",
}
CLASS_SHORT = {
    "Fidelity": "Fid",
    "Relay": "Rel",
    "Extension": "Ext",
    "Below_detection": "BD",
}

VMIN, VMAX = -0.04, 0.04


def main():
    set_nature_style()
    delta = pd.read_csv(DELTA_TSV, sep="\t")
    cls = pd.read_csv(CLS_TSV, sep="\t")
    hca = pd.read_csv(HCA_TSV, sep="\t")
    canonical = pd.read_csv(CANONICAL_ORDER_TSV, sep="\t")
    labels = pd.read_csv(ROOT / "results" / "ncp_v3_family_labels.tsv", sep="\t")

    fam_order = (
        canonical[["family_display_order", "family_id", "canonical_direction"]]
        .drop_duplicates("family_id")
        .sort_values("family_display_order")["family_id"]
        .astype(int)
        .tolist()
    )

    carriers = dict(zip(hca["family_id"].astype(int), hca["LA_carrier"].astype(str)))
    class_lookup = dict(
        zip(cls["family_id"].astype(int), cls["fidelity_class"].astype(str))
    )
    short_lookup = {
        int(r.family_id): short_pathway_label(str(r.family_label))
        for r in labels.itertuples(index=False)
    }
    dir_lookup = dict(
        zip(
            canonical.drop_duplicates("family_id")["family_id"].astype(int),
            canonical.drop_duplicates("family_id")["canonical_direction"].astype(str),
        )
    )

    delta["family_id"] = delta["family_id"].astype(int)
    delta_pivot = delta.pivot(index="family_id", columns="compartment", values="dScore")
    delta_pivot = delta_pivot.reindex(index=fam_order, columns=COMPARTMENTS)
    p_pivot = delta.pivot(index="family_id", columns="compartment", values="p_perm")
    p_pivot = p_pivot.reindex(index=fam_order, columns=COMPARTMENTS)
    incl_pivot = delta.pivot(index="family_id", columns="compartment", values="ct_included")
    incl_pivot = incl_pivot.reindex(index=fam_order, columns=COMPARTMENTS)

    n_rows = len(fam_order)
    n_cols = len(COMPARTMENTS)

    cmap = LinearSegmentedColormap.from_list(
        "AF_d", [(0.0, RA_DEEP), (0.45, "#FFFFFF"), (0.55, "#FFFFFF"), (1.0, LA_DEEP)]
    )
    cmap.set_under("#1F3650")
    cmap.set_over("#5C2D2A")
    cmap.set_bad("#E0E0E0")
    norm = TwoSlopeNorm(vmin=VMIN, vcenter=0.0, vmax=VMAX)

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 210 * MM))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(
        1, 6,
        width_ratios=[0.30, 0.42, 3.4, 6.6, 0.55, 1.35],
        left=0.04, right=0.985, top=0.90, bottom=0.10,
        wspace=0.04,
    )
    ax_carrier = fig.add_subplot(gs[0, 0])
    ax_dir = fig.add_subplot(gs[0, 1])
    ax_lbl = fig.add_subplot(gs[0, 2])
    ax_main = fig.add_subplot(gs[0, 3])
    ax_class = fig.add_subplot(gs[0, 4])
    ax_leg = fig.add_subplot(gs[0, 5])

    for i, fid in enumerate(fam_order):
        carrier = carriers.get(fid, MIXED_LABEL)
        if carrier == MIXED_LABEL:
            ax_carrier.add_patch(mpatches.Rectangle(
                (0, i), 1, 1,
                facecolor=MIXED_COLOR, edgecolor="white", linewidth=0.4, hatch="///",
            ))
        else:
            ax_carrier.add_patch(mpatches.Rectangle(
                (0, i), 1, 1,
                facecolor=COMPARTMENT_PALETTE.get(carrier, MIXED_COLOR),
                edgecolor="white", linewidth=0.4,
            ))
    ax_carrier.set_xlim(0, 1)
    ax_carrier.set_ylim(n_rows, 0)
    ax_carrier.set_xticks([0.5])
    ax_carrier.set_xticklabels(
        ["LA\nbaseline\ncarrier"], fontsize=enforce_font_floor(FS_TICK), color=FAINT_COLOR
    )
    ax_carrier.set_yticks([])
    ax_carrier.tick_params(axis="x", length=0, pad=2)
    for s in ax_carrier.spines.values():
        s.set_visible(False)

    for i, fid in enumerate(fam_order):
        d = dir_lookup.get(fid, "?")
        col = LA_DEEP if d == "LA-biased" else RA_DEEP
        ax_dir.add_patch(mpatches.Rectangle(
            (0, i), 1, 1, facecolor=col, edgecolor="white", linewidth=0.4
        ))
    ax_dir.set_xlim(0, 1)
    ax_dir.set_ylim(n_rows, 0)
    ax_dir.set_xticks([])
    ax_dir.set_yticks([])
    for s in ax_dir.spines.values():
        s.set_visible(False)

    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(n_rows, 0)
    ax_lbl.axis("off")
    for i, fid in enumerate(fam_order):
        d = dir_lookup.get(fid, "?")
        col = LA_DEEP if d == "LA-biased" else RA_DEEP
        ax_lbl.text(
            0.99, i + 0.5,
            f"F{int(fid)}  {short_lookup.get(fid, '')}",
            fontsize=enforce_font_floor(FS_TICK), color=col, fontweight="bold",
            ha="right", va="center",
        )

    mat = delta_pivot.to_numpy(dtype=float)
    pmat = p_pivot.to_numpy(dtype=float)
    incl = incl_pivot.to_numpy(dtype=object)
    masked = np.ma.masked_invalid(mat)
    im = ax_main.imshow(
        masked, cmap=cmap, norm=norm, aspect="auto",
        extent=[-0.5, n_cols - 0.5, n_rows - 0.5, -0.5],
        interpolation="nearest",
    )

    for i, fid in enumerate(fam_order):
        carrier = carriers.get(fid, MIXED_LABEL)
        for j, ct in enumerate(COMPARTMENTS):
            included = incl[i, j]
            if included is False or (isinstance(included, (int, float, np.bool_)) and not bool(included)):
                ax_main.add_patch(mpatches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    facecolor="#E8E8E8", edgecolor="#888888", linewidth=0.3,
                    hatch="xxx", zorder=3,
                ))
                continue
            p = pmat[i, j]
            if np.isfinite(p) and p < 0.05:
                ax_main.plot(
                    j, i, "o",
                    markersize=4.2, color="#1A1A1A",
                    markeredgecolor="white", markeredgewidth=0.35, zorder=5,
                )
            if carrier == ct:
                ax_main.add_patch(mpatches.Rectangle(
                    (j - 0.48, i - 0.48), 0.96, 0.96,
                    facecolor="none", edgecolor="#1A1A1A", linewidth=1.15, zorder=4,
                ))

    ax_main.set_xlim(-0.5, n_cols - 0.5)
    ax_main.set_ylim(n_rows - 0.5, -0.5)
    ax_main.set_xticks(np.arange(n_cols))
    ax_main.set_xticklabels(
        [COMPARTMENT_SHORT[c] for c in COMPARTMENTS],
        rotation=30, ha="right", fontsize=enforce_font_floor(FS_TICK),
    )
    ax_main.set_yticks([])
    ax_main.set_xlabel("cell compartment (LA snRNA)", fontsize=enforce_font_floor(FS_AXIS_LABEL))
    for s in ax_main.spines.values():
        s.set_visible(False)

    # class column
    ax_class.set_xlim(0, 1)
    ax_class.set_ylim(n_rows, 0)
    for i, fid in enumerate(fam_order):
        cname = class_lookup.get(fid, "Below_detection")
        ax_class.add_patch(mpatches.Rectangle(
            (0, i), 1, 1,
            facecolor=CLASS_COLORS.get(cname, MIXED_COLOR),
            edgecolor="white", linewidth=0.4,
        ))
        ax_class.text(
            0.5, i + 0.5, CLASS_SHORT.get(cname, "?"),
            ha="center", va="center",
            fontsize=enforce_font_floor(FS_TICK - 0.5),
            color="white" if cname != "Below_detection" else "#333333",
            fontweight="bold",
        )
    ax_class.set_xticks([0.5])
    ax_class.set_xticklabels(["class"], fontsize=enforce_font_floor(FS_TICK), color=FAINT_COLOR)
    ax_class.set_yticks([])
    ax_class.tick_params(axis="x", length=0, pad=2)
    for s in ax_class.spines.values():
        s.set_visible(False)

    counts = cls["fidelity_class"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int)

    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(1, 0)
    ax_leg.axis("off")
    ax_leg.text(0.0, 0.04, "Color scale", fontsize=enforce_font_floor(FS_TICK),
                fontweight="bold", color=FG_COLOR)
    cax = fig.add_axes([0.91, 0.72, 0.018, 0.12])
    cbar = fig.colorbar(im, cax=cax, orientation="vertical")
    cbar.set_label("Δsingscore  (AF − Ctrl)", fontsize=enforce_font_floor(FS_TICK))
    cbar.ax.tick_params(labelsize=enforce_font_floor(FS_TICK - 0.5))
    cbar.set_ticks([VMIN, -0.02, 0, 0.02, VMAX])
    cbar.set_ticklabels([f"{VMIN}", "-0.02", "0", "+0.02", f"+{VMAX:g}"])

    ax_leg.text(0.0, 0.28, "Cell marks", fontsize=enforce_font_floor(FS_TICK),
                fontweight="bold", color=FG_COLOR)
    ax_leg.plot(0.05, 0.34, "o", markersize=5, color="#1A1A1A",
                transform=ax_leg.transData, clip_on=False)
    ax_leg.text(0.12, 0.34, "p < 0.05 (primary 10k)",
                fontsize=enforce_font_floor(FS_TICK - 0.5), va="center", color=FG_COLOR)
    ax_leg.add_patch(mpatches.Rectangle(
        (0.02, 0.38), 0.08, 0.028,
        facecolor="none", edgecolor="#1A1A1A", linewidth=1.1,
    ))
    ax_leg.text(0.12, 0.394, "HCA LA carrier",
                fontsize=enforce_font_floor(FS_TICK - 0.5), va="center")
    ax_leg.add_patch(mpatches.Rectangle(
        (0.02, 0.43), 0.08, 0.028,
        facecolor="#E8E8E8", edgecolor="#888", linewidth=0.4, hatch="xxx",
    ))
    ax_leg.text(0.12, 0.444, "Not included (NA)",
                fontsize=enforce_font_floor(FS_TICK - 0.5), va="center")

    ax_leg.text(0.0, 0.52, "Class", fontsize=enforce_font_floor(FS_TICK),
                fontweight="bold", color=FG_COLOR)
    ly = 0.57
    for cname in CLASS_ORDER:
        ax_leg.add_patch(mpatches.Rectangle(
            (0.0, ly), 0.10, 0.026, facecolor=CLASS_COLORS[cname], edgecolor="none"
        ))
        ax_leg.text(
            0.14, ly + 0.013,
            f"{cname.replace('_', '-')}  (n={int(counts.get(cname, 0))})",
            fontsize=enforce_font_floor(FS_TICK - 0.5), va="center",
        )
        ly += 0.035

    ax_leg.text(0.0, 0.78, "Bias", fontsize=enforce_font_floor(FS_TICK),
                fontweight="bold", color=FG_COLOR)
    ax_leg.add_patch(mpatches.Rectangle((0.0, 0.82), 0.10, 0.026, facecolor=LA_DEEP, edgecolor="none"))
    ax_leg.text(0.14, 0.833, "LA-biased", fontsize=enforce_font_floor(FS_TICK - 0.5))
    ax_leg.add_patch(mpatches.Rectangle((0.0, 0.86), 0.10, 0.026, facecolor=RA_DEEP, edgecolor="none"))
    ax_leg.text(0.14, 0.873, "RA-biased", fontsize=enforce_font_floor(FS_TICK - 0.5))

    fig.suptitle(
        "Fig 3b  AF Δsingscore at LA single-cell resolution "
        "(16 family × 9 compartment, GSE255612)",
        x=0.04, y=0.965, ha="left",
        fontsize=enforce_font_floor(FS_PANEL_TITLE + 1), fontweight="bold", color=FG_COLOR,
    )
    fig.text(
        0.04, 0.055,
        "Summary:  "
        f"{int(counts.get('Fidelity', 0))} Fidelity  |  "
        f"{int(counts.get('Relay', 0))} Relay  |  "
        f"{int(counts.get('Extension', 0))} Extension  |  "
        f"{int(counts.get('Below_detection', 0))} Below-detection",
        fontsize=enforce_font_floor(FS_TICK), fontweight="bold", color=FG_COLOR,
    )
    fig.text(
        0.04, 0.012,
        "Rows: frozen Figure 1D canonical family order.  "
        "Fill: AF−Ctrl family-mean exact singscore on selection-time frozen LE.  "
        "Black dots: nominal p_perm < 0.05 (primary 10,000-perm / seed 20260418).  "
        "Outlined cells: HCA LA carrier prior (interpretation only).  "
        "Hatched: ct_included = False.  "
        "Class column uses family-level nominal p < 0.05 only.  "
        "p < 0.10 trends and gene-level decompositions are Supporting Figures.",
        fontsize=enforce_font_floor(FS_FOOTNOTE), color=FAINT_COLOR,
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
