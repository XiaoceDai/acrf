#!/usr/bin/env python3
"""
Fig 5e v3 - Asymmetric mechanism cartoon (NCP_v3 Atlas-priority families).

v3 honest layout:  ALL atlas-priority *** families are LA-driven, so the
v2-style symmetric "RA-axis | center | LA-axis" cartoon is structurally
broken in v3.  Instead we show:

    Top banner   |  All 4 *** Priority families are LA-driven (4 / 4)
    Left tier    |  RA-axis empty + note (no high-priority families)
                    plus the 3 RA-driven * Silent families listed for
                    transparency
    Right tier   |  6 chips for LA-driven + Bilateral priority families
                    (F15, F11, F9, F10 = ***;  F2, F4 = * Bilateral)

Each chip carries 5 independent visual encodings:
    background  =  PriorityClass  (***=deep green, **=mid, *=light)
    top strip   =  HCA LA carrier cell type
    middle      =  Fig 4 rho_dur  +  Fig 3 fidelity_class
                                  +  Fig 2E replication_class
    bottom      =  PriorityStars + ChamberDriver
    border      =  cross_axis_class colour from 5C

Inputs:
    results/fig_integrated_af_rewiring_atlas.tsv
    results/fig5b_v3_triangulation.tsv
    results/fig5c_v3_cross_axis.tsv

Output:
    figures/fig5e_v3_mechanism_cartoon.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

ATLAS_TSV = ROOT / "results" / "fig_integrated_af_rewiring_atlas.tsv"
TRI_TSV = ROOT / "results" / "fig5b_v3_triangulation.tsv"
CROSS_TSV = ROOT / "results" / "fig5c_v3_cross_axis.tsv"
OUT_PDF = ROOT / "figures" / "fig5e_v3_mechanism_cartoon.pdf"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
ACCENT = "#C8A076"
GREY = "#666666"
LIGHT = "#F2F2F2"

# Cell-type carrier colours (match Fig 1F + 5b)
CARRIER_COLORS = {
    "atrialCM": "#4E79A7",
    "fibroblast": "#A0CBE8",
    "endothelial": "#F28E2B",
    "macrophage": "#59A14F",
    "lymphocyte": "#8CD17D",
    "pericyte": "#499894",
    "vascular_smooth_muscle": "#86BCB6",
    "adipocyte": "#E15759",
    "neuronal": "#79706E",
    "Mixed/Unresolved": "#BDBDBD",
}
CARRIER_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC",
    "adipocyte": "Adipo",
    "neuronal": "Neuro",
    "Mixed/Unresolved": "Mix",
}

# PriorityClass background palette
PRIORITY_BG = {
    "***": "#D9EAD3",  # deep light green
    "**":  "#EAF1E0",  # medium
    "*":   "#F5F5F0",  # very light
}
PRIORITY_BORDER = {
    "***": "#2E7B5E",
    "**":  "#7BB36C",
    "*":   "#BBBBBB",
}

# Cross-axis class border (echoes 5C colours)
CROSS_AXIS_BORDER = {
    "Triple-positive": "#2E7B5E",
    "Double-positive": "#7BB36C",
    "Single-axis":     "#C8A076",
    "Silent":          "#BBBBBB",
}

# Chamber-driver text/icon
DRIVER_COLOR = {
    "LA-driven": LA_DEEP,
    "RA-driven": RA_DEEP,
    "Bilateral": "#7E5BBF",
}


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def draw_chip(ax, x, y, w, h, info):
    """Render a single family chip at axes-coord (x, y) with size (w, h)."""
    pc = info["PriorityStars"]
    cax_class = info["cross_axis_class"]
    border_color = CROSS_AXIS_BORDER.get(cax_class, "#888")

    # outer chip background (PriorityClass tinted, cross-axis border)
    ax.add_patch(mpatches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.005,rounding_size=0.012",
        facecolor=PRIORITY_BG.get(pc, LIGHT),
        edgecolor=border_color,
        linewidth=1.3, zorder=2,
    ))

    # top strip: HCA LA carrier
    carrier = info["LA_carrier"]
    strip_h = h * 0.18
    strip_color = CARRIER_COLORS.get(carrier, "#BDBDBD")
    ax.add_patch(mpatches.FancyBboxPatch(
        (x + 0.003, y + h - strip_h - 0.002),
        w - 0.006, strip_h,
        boxstyle="round,pad=0.0,rounding_size=0.008",
        facecolor=strip_color,
        edgecolor="#888", linewidth=0.3, zorder=3,
    ))
    ax.text(x + w / 2.0, y + h - strip_h / 2.0 - 0.003,
              f"carrier: {CARRIER_SHORT.get(carrier, '?')}",
              ha="center", va="center",
              fontsize=6.0, fontweight="bold", color="#1A1A1A", zorder=4)

    # title row: F-id + short name (separated)
    title_y = y + h - strip_h - 0.020
    ax.text(x + 0.012, title_y, f"F{int(info['family_id'])}",
              ha="left", va="top", fontsize=8.4, fontweight="bold",
              color="#1A1A1A", zorder=4)
    ax.text(x + 0.045, title_y, str(info["short"])[:22],
              ha="left", va="top", fontsize=6.4, color="#1A1A1A", zorder=4)

    # middle stat lines
    mid_y = title_y - 0.024
    line_h = 0.022
    rho = float(info["rho_primary"])
    chro = str(info["chronicity_class"]).replace("_", "-")
    fid_cls = str(info["fidelity_class"]).replace("_", "-")
    repl = str(info["replication_class"])
    ax.text(x + 0.012, mid_y,
              f"L3 dur \u03c1 = {rho:+.2f}  ({chro})",
              fontsize=5.8, color="#1A1A1A", zorder=4)
    ax.text(x + 0.012, mid_y - line_h,
              f"L2 fidelity:  {fid_cls}",
              fontsize=5.8, color="#1A1A1A", zorder=4)
    ax.text(x + 0.012, mid_y - 2 * line_h,
              f"L1 replication:  {repl}",
              fontsize=5.8, color="#1A1A1A", zorder=4)
    ax.text(x + 0.012, mid_y - 3 * line_h,
              f"cross-axis:  {cax_class}",
              fontsize=5.8, fontweight="bold",
              color=border_color, zorder=4)

    # bottom row: PriorityStars + ChamberDriver
    bot_y = y + 0.012
    ax.text(x + 0.012, bot_y, info["PriorityStars"],
              ha="left", va="bottom", fontsize=10.0,
              fontweight="bold", color=PRIORITY_BORDER.get(pc, "#888"),
              zorder=4)
    drv = str(info["ChamberDriver"])
    ax.text(x + w - 0.012, bot_y, drv,
              ha="right", va="bottom", fontsize=6.4, fontweight="bold",
              color=DRIVER_COLOR.get(drv, "#444"),
              zorder=4)


def main():
    set_style()
    atlas = pd.read_csv(ATLAS_TSV, sep="\t")
    tri = pd.read_csv(TRI_TSV, sep="\t")
    cross = pd.read_csv(CROSS_TSV, sep="\t")

    # Merge into a single dict per family.  Triangulation TSV duplicates
    # 'LA_carrier' / 'replication_class'; pull only the unique columns.
    tri_cols = [c for c in ["fidelity_class", "chronicity_class", "rho_primary"]
                  if c in tri.columns]
    info_all = (
        atlas.set_index("family_id")
             .join(tri.set_index("family_id")[tri_cols], how="left")
             .join(cross.set_index("family_id")[
                 ["cross_axis_class", "cross_axis_subtype", "n_supportive"]],
                   how="left")
    )
    # Pull replication_class from atlas/triangulation (Atlas already has it
    # under "replication_class")
    if "replication_class" not in info_all.columns:
        info_all = info_all.join(tri.set_index("family_id")[["replication_class"]],
                                   how="left")
    info_all["family_id"] = info_all.index.astype(int)

    # 6 chip families (atlas top + biology anchors)
    chip_fids = [15, 11, 9, 10, 2, 4]
    chips = [info_all.loc[fid].to_dict() | {"family_id": fid}
              for fid in chip_fids]

    # ---- figure ----
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(220 * MM, 175 * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # ---- title ----
    fig.text(0.025, 0.965,
              "Fig 5e v3   AF rewiring mechanism cartoon  (asymmetric, NCP_v3)",
              fontsize=10.0, fontweight="bold", color="#1A1A1A")

    # ---- BANNER (full width) ----
    banner_y = 0.875
    banner_h = 0.060
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.025, banner_y), 0.95, banner_h,
        boxstyle="round,pad=0.005,rounding_size=0.012",
        facecolor=LA_DEEP, edgecolor="white", linewidth=1.0,
    ))
    ax.text(0.5, banner_y + banner_h / 2.0,
              "All 4  *** Priority families are LA-driven  (4 / 4)",
              ha="center", va="center",
              fontsize=11.2, fontweight="bold", color="white")
    ax.text(0.97, banner_y - 0.012,
              "Source: atlas ChamberDriver column  (Fig 2D pooled g + per-cohort delta_LA / delta_RA)",
              ha="right", va="top", fontsize=5.6, color="#666666", style="italic")

    # ---- LEFT TIER: RA-axis empty + note ----
    left_x, left_w = 0.025, 0.30
    left_y, left_h = 0.10, 0.74
    ax.add_patch(mpatches.FancyBboxPatch(
        (left_x, left_y), left_w, left_h,
        boxstyle="round,pad=0.005,rounding_size=0.012",
        facecolor="#FAFAFA", edgecolor=RA_DEEP, linewidth=0.8,
    ))
    ax.text(left_x + left_w / 2.0, left_y + left_h - 0.030,
              "RA-axis  (left tier)",
              ha="center", va="top",
              fontsize=8.6, fontweight="bold", color=RA_DEEP)
    # Big "empty" callout
    ax.text(left_x + left_w / 2.0, left_y + left_h - 0.110,
              "no high-priority families",
              ha="center", va="center",
              fontsize=8.6, color=RA_DEEP, fontweight="bold", style="italic")
    # The 3 RA-driven Silent families listed for transparency
    ra_silent_fids = [7, 8, 6]
    ra_y = left_y + left_h - 0.18
    ax.text(left_x + 0.018, ra_y,
              "RA-driven families (all  *  PriorityClass):",
              fontsize=6.4, fontweight="bold", color="#1A1A1A")
    ra_y -= 0.025
    for fid in ra_silent_fids:
        info = info_all.loc[fid]
        cax = str(info["cross_axis_class"])
        carrier = str(info["LA_carrier"])
        ax.text(left_x + 0.030, ra_y,
                  f"\u25CF  F{int(fid)}  {str(info['short'])[:18]}",
                  fontsize=6.3, fontweight="bold", color=RA_DEEP)
        ax.text(left_x + 0.030, ra_y - 0.018,
                  f"      carrier: {CARRIER_SHORT.get(carrier, '?')}   "
                  f"\u03c1 = {float(info['rho_primary']):+.2f}   "
                  f"({cax})",
                  fontsize=5.6, color="#444444")
        ra_y -= 0.052
    # interpretation footer
    foot_y = left_y + 0.065
    ax.text(left_x + 0.018, foot_y,
              "Interpretation:",
              fontsize=6.4, fontweight="bold", color="#1A1A1A")
    ax.text(left_x + 0.018, foot_y - 0.020,
              "AF chamber rewiring at the highest-priority\n"
              "tier is concentrated on the LA arm; the RA\n"
              "arm shows only low-priority single-axis\n"
              "trends (mostly metabolic OXPHOS-related).",
              fontsize=5.7, color="#333333", va="top")

    # ---- RIGHT TIER: 6 chips ----
    right_x = 0.355
    right_w = 0.62
    right_y, right_h = 0.10, 0.74
    ax.add_patch(mpatches.FancyBboxPatch(
        (right_x, right_y), right_w, right_h,
        boxstyle="round,pad=0.005,rounding_size=0.012",
        facecolor="#FFFCF6", edgecolor=LA_DEEP, linewidth=0.8,
    ))
    ax.text(right_x + right_w / 2.0, right_y + right_h - 0.030,
              "LA-axis  (right tier)   \u2022   6 priority families  (4 LA-driven  +  2 Bilateral)",
              ha="center", va="top",
              fontsize=8.6, fontweight="bold", color=LA_DEEP)

    # 6 chips arranged in 3 rows x 2 cols
    chip_w = 0.285
    chip_h = 0.198
    chip_xs = [right_x + 0.012, right_x + 0.012 + chip_w + 0.012]
    chip_ys = [right_y + right_h - 0.080 - chip_h,
                 right_y + right_h - 0.080 - 2 * chip_h - 0.012,
                 right_y + right_h - 0.080 - 3 * chip_h - 2 * 0.012]
    for k, info in enumerate(chips):
        col = k % 2
        row = k // 2
        draw_chip(ax, chip_xs[col], chip_ys[row], chip_w, chip_h, info)

    # ---- legend (bottom strip) ----
    leg_y = 0.060
    ax.add_patch(mpatches.Rectangle((0.025, 0.020), 0.95, leg_y,
                                      facecolor="#FFF8F0", edgecolor=ACCENT, linewidth=0.5))
    ax.text(0.040, 0.063,
              "Chip encoding:",
              fontsize=6.5, fontweight="bold", color="#1A1A1A", va="top")
    ax.text(0.135, 0.063,
              "background = PriorityClass (***/**/*);  border colour = cross-axis class (5C);  top strip = HCA LA carrier (cell-type, Fig 1F);  "
              "middle = L3 dur \u03c1 + chronicity / L2 fidelity / L1 replication;  bottom = PriorityStars + ChamberDriver.",
              fontsize=5.7, color="#333333", va="top")
    ax.text(0.040, 0.040,
              "Honest framing:",
              fontsize=6.5, fontweight="bold", color="#1A1A1A", va="top")
    ax.text(0.135, 0.040,
              "Asymmetry is the message:  no atlas Top-priority family is RA-driven.  "
              "Bilateral chips (F2, F4) are placed in the right tier because Bilateral includes the LA arm.",
              fontsize=5.7, color="#333333", va="top")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
