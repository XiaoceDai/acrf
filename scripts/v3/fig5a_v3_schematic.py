#!/usr/bin/env python3
"""
Fig 5a v3 - Cross-figure synthesis schematic.

A pure schematic explaining how Fig 5 integrates Fig 1F / Fig 2 / Fig 3 /
Fig 4 evidence into a cross-axis convergence classification, then anchors
the result with GeoMx ROI spatial support and a mechanism cartoon.

Flow:
    16 NCP axis family anchor
        |
        v
    5 evidence streams in parallel
        L0  HCA LA carrier  (Fig 1F)
        L1  Bulk pooled g + state + replication  (Fig 2D / 2C / 2E)
        L2  snRNA dScore at carrier + fidelity   (Fig 3)
        L3  LAA duration rho + chronicity        (Fig 4)
        Atlas  PriorityScore + Class             (atlas oncoprint)
        |
        v
    Path A strict supportive criteria
    (L1 = repl=high; L2 = fidelity != Below; L3 = chron != Stable)
        |
        v
    Cross-axis convergence class  (Triple / Double / Single / Silent)
        |     |     |
        v     v     v
    5b heatmap  5c bar+matrix  5d GeoMx (6-fam v3 LE)  5e cartoon
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT_PDF = ROOT / "figures" / "fig5a_schematic.pdf"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
ACCENT = "#C8A076"
GREY = "#666666"
LIGHT = "#F2F2F2"
HL_GREEN = "#3C8C57"
BG_GOLD = "#FFFCF4"


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def stage_box(ax, x, y, w, h, label_top, label_body, *,
                color=LIGHT, edgecolor="#888", body_align="center",
                title_fs=7.2, body_fs=5.9, body_color="#333333"):
    ax.add_patch(mpatches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.005,rounding_size=0.012",
        facecolor=color, edgecolor=edgecolor, linewidth=0.6,
    ))
    ax.text(x + w / 2.0, y + h - 0.020, label_top,
              ha="center", va="top",
              fontsize=title_fs, fontweight="bold", color="#1A1A1A")
    if body_align == "left":
        ax.text(x + 0.010, y + h - 0.052, label_body,
                  ha="left", va="top", fontsize=body_fs, color=body_color)
    else:
        ax.text(x + w / 2.0, y + h / 2.0 - 0.005, label_body,
                  ha="center", va="center",
                  fontsize=body_fs, color=body_color)


def arrow(ax, x1, y1, x2, y2, *, label=None, lw=0.7, color="#1A1A1A",
            label_color="#444444", label_fs=5.6):
    ax.annotate(
        "", xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle="->", color=color, lw=lw),
    )
    if label:
        ax.text((x1 + x2) / 2.0, (y1 + y2) / 2.0 + 0.010, label,
                  fontsize=label_fs, ha="center", va="bottom",
                  color=label_color, style="italic")


def main():
    set_style()
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(220 * MM, 175 * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # ---- title ----
    ax.text(0.025, 0.965,
              "Fig 5a  Cross-figure synthesis pipeline (NCP, 16 families)",
              fontsize=9.6, fontweight="bold", color="#1A1A1A")
    ax.text(0.025, 0.940,
              "Integrates Fig 1F (carrier prior), Fig 2 (bulk rewiring), Fig 3 (snRNA carrier exec) "
              "and Fig 4 (duration trend) under a strict per-axis supportive criterion (Path A).",
              fontsize=6.5, color="#444444", style="italic")

    # ---- Stage 0: anchor ----
    stage_box(ax, 0.025, 0.78, 0.18, 0.10,
              "0.  NCP axis family anchor",
              "16 axis families\nfrozen panel\n(Fig 1c-e)",
              color=BG_GOLD, edgecolor=ACCENT, body_align="left",
              title_fs=7.0, body_fs=5.7)

    # ---- Stage 1: 5 evidence streams (horizontal row) ----
    stream_y = 0.55
    stream_h = 0.18
    stream_w = 0.17
    stream_xs = [0.025, 0.220, 0.415, 0.610, 0.805]
    stream_specs = [
        ("L0  HCA LA carrier",
         "Fig 1e\n"
         "HCA Heart Atlas\n"
         "(LA + RA paired)\n"
         "cell-type per family\n"
         "+ consensus class",
         "#FFF4E6", ACCENT),
        ("L1  Bulk rewiring",
         "Fig 2C / 2D / 2E\n"
         "3 cohorts pooled g\n"
         "+ rewiring state\n"
         "+ replication class\n"
         "(low / mod / high)",
         "#EAF1F7", RA_DEEP),
        ("L2  snRNA carrier exec",
         "Fig 3\n"
         "GSE255612 LA snRNA\n"
         "rank singscore (LE)\n"
         "dScore at carrier\n"
         "+ fidelity class",
         "#F2EAEA", LA_DEEP),
        ("L3  LAA duration trend",
         "Fig 4\n"
         "GSE41177 LAA n=12\n"
         "Spearman rho on\n"
         "duration months\n"
         "+ chronicity class",
         "#FFF8E0", "#C8A076"),
        ("Atlas  Integrated",
         "Atlas oncoprint\n"
         "PriorityScore =\n"
         "E + R + C + H\n"
         "PriorityClass\n"
         "(*** / ** / *)",
         "#E8F0E8", HL_GREEN),
    ]
    for x, (top, body, fc, ec) in zip(stream_xs, stream_specs):
        stage_box(ax, x, stream_y, stream_w, stream_h,
                    top, body,
                    color=fc, edgecolor=ec, body_align="left",
                    title_fs=6.8, body_fs=5.5)

    # arrows from anchor box down to each stream box
    for x in stream_xs:
        arrow(ax, 0.115, 0.78, x + stream_w / 2.0, stream_y + stream_h)

    # ---- Stage 2: Path A supportive criteria (centred wide box) ----
    stage_box(ax, 0.150, 0.355, 0.700, 0.135,
              "2.  Path A strict supportive criteria  (per-axis equal threshold)",
              "L1 supportive  =  replication_class = high  (CI excludes 0;  3 cohorts pooled)\n"
              "L2 supportive  =  fidelity_class  in  {Fidelity, Relay, Extension}  "
              "(i.e. \u2260 Below_detection;  significant cell-type-specific signal in \u2265 1 compartment)\n"
              "L3 supportive  =  chronicity_class \u2260 Stable  (|\u03c1| \u2265 0.2 directional trend)\n"
              "                  HCA carrier (L0) provides healthy prior; not used in supportive count.",
              color="#F8F8F8", edgecolor="#888", body_align="left",
              title_fs=7.0, body_fs=5.8)

    # arrows from streams to Path A box
    for x in stream_xs:
        arrow(ax, x + stream_w / 2.0, stream_y, 0.50, 0.49,
                lw=0.5, color="#888", label=None)

    # ---- Stage 3: Cross-axis convergence (4 outcome chips) ----
    chip_y = 0.305
    chip_h = 0.080
    chip_specs = [
        ("Triple-positive",  1, "#2E7B5E", "white"),
        ("Double-positive",  2, "#7BB36C", "#1A1A1A"),
        ("Single-axis",      7, "#C8C39A", "#1A1A1A"),
        ("Silent",           6, "#BBBBBB", "white"),
    ]
    chip_w = 0.16
    chip_xs = [0.085, 0.265, 0.445, 0.625]
    for x, (cls_name, c, fc, txt_color) in zip(chip_xs, chip_specs):
        ax.add_patch(mpatches.FancyBboxPatch(
            (x, chip_y), chip_w, chip_h,
            boxstyle="round,pad=0.006,rounding_size=0.012",
            facecolor=fc, edgecolor="white", linewidth=1.2,
        ))
        ax.text(x + chip_w / 2.0, chip_y + chip_h - 0.022, cls_name,
                  ha="center", va="top", fontsize=7.0,
                  fontweight="bold", color=txt_color)
        ax.text(x + chip_w / 2.0, chip_y + chip_h / 2.0 - 0.008,
                  f"n = {c}",
                  ha="center", va="center", fontsize=11.0,
                  fontweight="bold", color=txt_color)
    # arrow from Path A to chips row
    arrow(ax, 0.50, 0.355, 0.50, chip_y + chip_h, lw=0.7,
            label="cross-axis count")

    # rightmost: output panels (sits beside the chip row at the same y)
    stage_box(ax, 0.810, chip_y, 0.165, chip_h,
                "Output panels",
                "5b  triangulation heatmap\n"
                "5c  cross-axis class\n"
                "5d  GeoMx (5d below)\n"
                "5e  mechanism cartoon",
                color="#FFF4E6", edgecolor=ACCENT, body_align="left",
                title_fs=6.4, body_fs=5.2)
    arrow(ax, 0.785, chip_y + chip_h / 2.0, 0.810, chip_y + chip_h / 2.0,
            lw=0.7)

    # ---- Stage 4: GeoMx spatial corroboration row ----
    stage_box(ax, 0.025, 0.175, 0.45, 0.105,
                "5d  Spatial corroboration  (support-only, NOT discovery)",
                "GSE261363 GeoMx WTA\n"
                "n = 24 atrial-appendage ROIs (border / distal x AF / SR)\n"
                "6 Atlas-priority families x NCP LE rank singscore\n"
                "  -> family x cell-type Spearman matrix\n"
                "  -> AF - SR \u0394, Border - Distal \u0394 (Welch CI)\n"
                "Selection: F15, F11, F9, F10, F2, F4 (covers structure + metabolism axes)",
                color="#EAF1F7", edgecolor=RA_DEEP, body_align="left",
                title_fs=6.6, body_fs=5.3)

    # ---- Stage 5: cartoon / asymmetry banner ----
    stage_box(ax, 0.495, 0.175, 0.480, 0.105,
                "5e  Mechanism cartoon  (asymmetric, integrative layout)",
                "Banner:  All 4 *** Priority families are LA-driven (4 / 4)\n"
                "Right tier (LA-axis):  6 chips for top-priority families\n"
                "  F15 (Triple), F11, F9, F10 (***);  F16, F3 (**);  F2, F4 (Bilateral *)\n"
                "Left tier (RA-axis):  empty + note 'no high-priority families'\n"
                "Each chip: bg = PriorityClass;  top strip = LA carrier (cell-type);\n"
                "  middle = \u03c1_dur + fidelity + repl;  bottom = ***  ChamberDriver",
                color="#F2EAEA", edgecolor=LA_DEEP, body_align="left",
                title_fs=6.6, body_fs=5.3)

    # ---- footer banner (3 framing lines) ----
    ax.add_patch(mpatches.Rectangle((0.025, 0.020), 0.95, 0.135,
                                     facecolor="#FFF8F0", edgecolor=ACCENT, linewidth=0.6))
    ax.text(0.040, 0.140, "Integrative framing:",
              fontsize=6.5, fontweight="bold", color="#1A1A1A", va="top")
    ax.text(0.155, 0.140,
              "no single integrated metric replaces per-axis classes; Fig 5 chains the three already-published "
              "class systems (Fig 2 state, Fig 3 fidelity, Fig 4 chronicity) under a strict equal-threshold supportive criterion.",
              fontsize=5.7, color="#333333", va="top")
    ax.text(0.040, 0.097, "Power note:",
              fontsize=6.5, fontweight="bold", color="#1A1A1A", va="top")
    ax.text(0.155, 0.097,
              "L1 (n = 49 donors pooled across 3 cohorts) > L2 (n = 34 LA snRNA) > L3 (n = 12 LAA bulk).  "
              "GeoMx layer (n = 24 ROIs) provides spatial localization corroboration only; mechanistic causality is not inferred from 5d.",
              fontsize=5.7, color="#333333", va="top")
    ax.text(0.040, 0.054, "Atlas integration:",
              fontsize=6.5, fontweight="bold", color="#1A1A1A", va="top")
    ax.text(0.155, 0.054,
              "PriorityScore  =  E (effect strength)  +  R (replication)  +  C (state concordance)  +  H (heterogeneity penalty); "
              "PriorityScore drives the row order in 5b / 5c / 5e and the family selection in 5d.",
              fontsize=5.7, color="#333333", va="top")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
