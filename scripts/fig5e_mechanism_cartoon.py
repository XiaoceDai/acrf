#!/usr/bin/env python3
"""
Figure 5E  Dual-axis mechanism map  (post-NCP_v2 narrative).

Hybrid 3-column layout:
    LEFT      RA-axis evidence chips (SM-1, SM-2, SM-3)
    CENTER    Atrial-wall tissue cartoon with cell-type icons + directional
              arrows;  serves as the visual hub for the seven sub-mechanisms.
    RIGHT     LA-axis evidence chips (SM-4, SM-5, SM-6/7)

Top banner:    cross-axis directional-asymmetry headline
                (21 RA→LA  vs  3 LA→RA reversal genes,  binomtest p = 2.8e-4)

Output:  figures/fig5e_mechanism_cartoon.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Rectangle

ROOT = Path(__file__).resolve().parent.parent
FIG_OUT = ROOT / "figures" / "fig5e_mechanism_cartoon.pdf"

LA_COLOR    = "#D6604D"
RA_COLOR    = "#4393C3"
BG_LA       = "#E6F0F7"
BG_RA       = "#FBE7E2"
BG_LA_DEEP  = "#C8DCEB"
BG_RA_DEEP  = "#F4D2C8"
GREY_LIGHT  = "#E6E6E6"
GREY_MID    = "#9C9C9C"
GREY_DARK   = "#444444"
BG_GREY     = "#FAFAFA"
BG_GOLD     = "#FFFCF4"
BG_GREEN    = "#F2F8F4"
HL_GREEN    = "#3C8C57"
BRBG_BROWN  = "#8C510A"
BRBG_TEAL   = "#01665E"
PURPLE      = "#7E5BBF"

# Cell-icon colours
COL_FIB   = "#66A61E"
COL_CM    = "#E6AB02"
COL_ADIPO = "#B79268"
COL_LYMPH = "#7DA3B5"
COL_MAC   = "#D95F02"
COL_VENDO = "#7570B3"


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def add_box(ax, x, y, w, h, text="", *,
              fc=BG_GREY, ec=GREY_DARK, text_color="black",
              fontsize=7.0, fontweight="normal",
              radius=0.012, lw=0.7):
    box = mpatches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.005,rounding_size={radius}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2)
    ax.add_patch(box)
    if text:
        ax.text(x + w / 2, y + h / 2, text,
                 ha="center", va="center",
                 fontsize=fontsize, fontweight=fontweight,
                 color=text_color, zorder=3)


def draw_evidence_chip(ax, x, y, w, h, sm_id, title, color,
                          evidence_lines):
    """Narrow evidence chip with SM tag, colored title, evidence lines.
    evidence_lines: list of (label, value, is_star)."""
    add_box(ax, x, y, w, h, "",
              fc="white", ec=color, lw=1.0, radius=0.008)
    # SM badge
    ax.text(x + 0.018, y + h - 0.013, sm_id,
             ha="center", va="center",
             fontsize=5.8, fontweight="bold",
             color="white",
             bbox=dict(boxstyle="round,pad=0.18", facecolor=color,
                         edgecolor="none"),
             zorder=4)
    # title
    ax.text(x + 0.040, y + h - 0.013, title,
             ha="left", va="center",
             fontsize=6.6, fontweight="bold",
             color=color, zorder=4)
    # evidence
    for k, line in enumerate(evidence_lines):
        ax.text(x + 0.010, y + h - 0.030 - k * 0.014,
                 line,
                 ha="left", va="center",
                 fontsize=5.3, color=GREY_DARK, zorder=4)


# ---- cell-type icons -------------------------------------------------------

def draw_fibroblast(ax, x, y, scale=0.022, color=COL_FIB, label=True):
    """Fibroblast: spindle (rotated thin ellipse)."""
    e = mpatches.Ellipse((x, y), scale * 1.7, scale * 0.55,
                            angle=20, facecolor=color, edgecolor="black",
                            linewidth=0.4, zorder=5)
    ax.add_patch(e)
    if label:
        ax.text(x, y + scale * 0.55, "fibroblast",
                  ha="center", va="bottom", fontsize=5.0,
                  color=color, fontweight="bold", zorder=6)


def draw_cm(ax, x, y, scale=0.024, color=COL_CM, label=True):
    """Cardiomyocyte: rectangle with horizontal stripes (sarcomeric)."""
    rect = Rectangle((x - scale * 0.9, y - scale * 0.4),
                       scale * 1.8, scale * 0.8,
                       facecolor=color, edgecolor="black",
                       linewidth=0.5, zorder=5)
    ax.add_patch(rect)
    # internal stripes
    for k in range(1, 4):
        ax.plot([x - scale * 0.9, x + scale * 0.9],
                  [y - scale * 0.4 + k * scale * 0.20,
                   y - scale * 0.4 + k * scale * 0.20],
                  color="white", lw=0.4, zorder=6)
    if label:
        ax.text(x, y + scale * 0.55, "atrialCM",
                  ha="center", va="bottom", fontsize=5.0,
                  color=color, fontweight="bold", zorder=6)


def draw_adipocyte(ax, x, y, scale=0.020, color=COL_ADIPO, label=True):
    """Adipocyte: large circle with inner lipid droplet."""
    c = Circle((x, y), scale, facecolor=color, alpha=0.6,
                  edgecolor="black", linewidth=0.5, zorder=5)
    ax.add_patch(c)
    inner = Circle((x - scale * 0.15, y + scale * 0.15),
                      scale * 0.55, facecolor="white",
                      edgecolor=color, linewidth=0.4, zorder=6)
    ax.add_patch(inner)
    if label:
        ax.text(x, y + scale * 1.25, "adipocyte",
                  ha="center", va="bottom", fontsize=5.0,
                  color=color, fontweight="bold", zorder=6)


def draw_lymph_endo(ax, x, y, scale=0.022, color=COL_LYMPH, label=True):
    """Lymphatic endothelial: thin elongated tube."""
    e = mpatches.Ellipse((x, y), scale * 1.6, scale * 0.45,
                            angle=0, facecolor=color, alpha=0.6,
                            edgecolor="black", linewidth=0.4, zorder=5)
    ax.add_patch(e)
    # inner lumen line
    ax.plot([x - scale * 0.7, x + scale * 0.7], [y, y],
              color="white", lw=0.5, zorder=6)
    if label:
        ax.text(x, y + scale * 0.45, "lymph endo",
                  ha="center", va="bottom", fontsize=5.0,
                  color=color, fontweight="bold", zorder=6)


def draw_macrophage(ax, x, y, scale=0.018, color=COL_MAC, label=True):
    """Macrophage: circle with small dots around."""
    c = Circle((x, y), scale * 0.7, facecolor=color, alpha=0.7,
                  edgecolor="black", linewidth=0.4, zorder=5)
    ax.add_patch(c)
    # dots = macrophage protrusions / bumps
    import numpy as np
    for ang in np.linspace(0, 360, 5, endpoint=False):
        rad = (ang) * 3.14159 / 180
        bx = x + scale * 1.0 * np.cos(rad)
        by = y + scale * 1.0 * np.sin(rad)
        ax.add_patch(Circle((bx, by), scale * 0.18,
                                facecolor=color, edgecolor="none",
                                zorder=5))
    if label:
        ax.text(x, y + scale * 1.4, "macrophage",
                  ha="center", va="bottom", fontsize=5.0,
                  color=color, fontweight="bold", zorder=6)


def draw_arrow(ax, x0, y0, x1, y1, *, color=GREY_DARK, lw=1.0,
                  style="-|>", curve=0.0, label=None,
                  label_color=None, label_offset=(0, 0)):
    """Draw a curved arrow between (x0,y0) and (x1,y1)."""
    arr = FancyArrowPatch((x0, y0), (x1, y1),
                              arrowstyle=style,
                              connectionstyle=f"arc3,rad={curve}",
                              color=color, linewidth=lw,
                              mutation_scale=12, zorder=4)
    ax.add_patch(arr)
    if label:
        ax.text((x0 + x1) / 2 + label_offset[0],
                 (y0 + y1) / 2 + label_offset[1],
                 label, fontsize=5.5, fontweight="bold",
                 color=(label_color or color),
                 ha="center", va="center",
                 bbox=dict(boxstyle="round,pad=0.15",
                              facecolor="white", edgecolor="none"),
                 zorder=5)


def draw_tissue_cartoon(ax, x0, y0, w, h):
    """Center cartoon: stylized atrial cross-section with LA + RA blocks
    + cell icons + sub-mechanism flow arrows."""

    # background panel
    add_box(ax, x0, y0, w, h, "",
              fc="#FFFFFF", ec=GREY_MID, lw=0.6, radius=0.012)

    # title
    ax.text(x0 + w / 2, y0 + h - 0.020, "Atrial-wall mechanism schematic",
             ha="center", va="center", fontsize=7.5, fontweight="bold",
             color=GREY_DARK)

    # ===== LA / RA tissue blocks =====
    # vertical layout: LA block on TOP, arrows going down to LA cells
    # RA block on the LEFT, arrows curving to LA tissue (showing invasion)
    pad = 0.012
    band_h = (h - 0.080) / 2 - 0.005
    band_w = w - pad * 2

    # LA tissue block (lower) - "target" of invasion
    la_x, la_y = x0 + pad, y0 + 0.020
    la_w, la_h = band_w, band_h
    add_box(ax, la_x, la_y, la_w, la_h, "",
              fc=BG_LA, ec=LA_COLOR, lw=0.8, radius=0.008)
    ax.text(la_x + 0.010, la_y + la_h - 0.012,
             "LA  (left atrium)  —  TARGET",
             ha="left", va="center", fontsize=6.5,
             fontweight="bold", color=LA_COLOR)

    # RA tissue block (upper)
    ra_x, ra_y = x0 + pad, y0 + band_h + 0.045
    ra_w, ra_h = band_w, band_h
    add_box(ax, ra_x, ra_y, ra_w, ra_h, "",
              fc=BG_RA, ec=RA_COLOR, lw=0.8, radius=0.008)
    ax.text(ra_x + 0.010, ra_y + ra_h - 0.012,
             "RA  (right atrium)  —  SOURCE",
             ha="left", va="center", fontsize=6.5,
             fontweight="bold", color=RA_COLOR)

    # ===== Cell icons distributed within LA tissue =====
    cy_la = la_y + la_h * 0.45
    # Position cells within LA from left to right
    fib_x  = la_x + la_w * 0.15
    cm_x   = la_x + la_w * 0.40
    adip_x = la_x + la_w * 0.65
    lymph_x = la_x + la_w * 0.85

    draw_fibroblast(ax, fib_x, cy_la, scale=0.015)
    draw_cm(ax, cm_x, cy_la, scale=0.018)
    draw_adipocyte(ax, adip_x, cy_la, scale=0.014)
    draw_lymph_endo(ax, lymph_x, cy_la, scale=0.015)

    # ===== Sub-mechanism arrows in LA block =====

    # SM-1 TGF-β / SLRP fibrosis: RA → LA fibroblast (matrix invasion)
    # curved brown arrow from upper RA tissue middle down to LA fibroblast
    draw_arrow(ax,
                  ra_x + ra_w * 0.20, ra_y + 0.005,
                  fib_x, cy_la + 0.020,
                  color=BRBG_BROWN, lw=1.8, curve=-0.30,
                  label="SM-1\nTGF-β / SLRP", label_color=BRBG_BROWN,
                  label_offset=(-0.04, 0.005))

    # SM-2a BMP gain in adipocyte: small UP arrow on the adipocyte
    ax.annotate("", xy=(adip_x, cy_la - 0.030),
                  xytext=(adip_x, cy_la - 0.045),
                  arrowprops=dict(arrowstyle="-|>", color=PURPLE, lw=1.4,
                                    mutation_scale=10), zorder=5)
    ax.text(adip_x + 0.013, cy_la - 0.038, "SM-2 BMP ↑",
              fontsize=5.5, fontweight="bold", color=PURPLE,
              ha="left", va="center", zorder=5)

    # SM-2b BMP loss in lymph endo: small DOWN arrow
    ax.annotate("", xy=(lymph_x, cy_la - 0.045),
                  xytext=(lymph_x, cy_la - 0.030),
                  arrowprops=dict(arrowstyle="-|>", color=PURPLE, lw=1.4,
                                    mutation_scale=10), zorder=5)
    ax.text(lymph_x + 0.013, cy_la - 0.038, "SM-2 BMP ↓",
              fontsize=5.5, fontweight="bold", color=PURPLE,
              ha="left", va="center", zorder=5)

    # SM-4 + SM-5 LA chronic CM decay: time-curve arrow next to CM
    # downward curving arrow showing chronic decline
    draw_arrow(ax,
                  cm_x - 0.005, cy_la - 0.018,
                  cm_x - 0.025, cy_la - 0.045,
                  color=LA_COLOR, lw=1.4, curve=-0.30,
                  label="SM-4/5\nMito Ca²⁺ + OXPHOS\nchronic decline",
                  label_color=LA_COLOR,
                  label_offset=(-0.040, -0.005))
    # add small "time" arrow icon
    ax.annotate("", xy=(cm_x + 0.030, cy_la - 0.045),
                  xytext=(cm_x + 0.005, cy_la - 0.045),
                  arrowprops=dict(arrowstyle="->", color=GREY_MID,
                                    lw=0.8, mutation_scale=8), zorder=5)
    ax.text(cm_x + 0.040, cy_la - 0.045, "AF months",
              fontsize=5.0, color=GREY_MID, style="italic",
              ha="left", va="center", zorder=5)

    # ===== RA tissue — annotations (just labels, no detailed cells) =====
    # show 21:3 directional asymmetry as the dominant arrow
    ax.text(ra_x + ra_w * 0.55, ra_y + ra_h * 0.45,
              "21 RA-biased reversal genes\n→ invade LA tissue",
              ha="center", va="center",
              fontsize=6.0, color=BRBG_BROWN, fontweight="bold",
              bbox=dict(boxstyle="round,pad=0.30", facecolor=BG_RA,
                          edgecolor=RA_COLOR, linewidth=0.5))

    # SM-3 HS-GAG diffuse - small note in LA tissue corner
    ax.text(la_x + la_w * 0.55, la_y + 0.012,
              "SM-3  HS-GAG metabolism diffuse  (no single carrier)",
              ha="center", va="bottom", fontsize=5.4,
              color=BRBG_BROWN, style="italic")

    # SM-6/7 dediff small note
    ax.text(la_x + la_w * 0.10, la_y + 0.012,
              "SM-6/7  ACTA1 + ZFP57  loss",
              ha="left", va="bottom", fontsize=5.4,
              color=LA_COLOR, style="italic")


def main() -> None:
    set_publication_style()

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(230 * MM, 175 * MM))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # ==== TITLE ====
    fig.text(0.022, 0.965,
              "Fig 5E  Dual-axis mechanism map of AF atrial remodeling",
              fontsize=9.5, fontweight="bold", ha="left", va="center")
    fig.text(0.022, 0.937,
              "Seven sub-mechanisms (SM-1 to SM-7) + cross-axis directional "
              "asymmetry, derived from Fig 2 / 3 / 4 / 5 / 6 + N2d statistics.",
              fontsize=7.5, color=GREY_DARK, ha="left", va="center",
              style="italic")

    # ==== TOP BANNER ====
    bx, by, bw, bh = 0.020, 0.860, 0.960, 0.045
    add_box(ax, bx, by, bw, bh, "",
              fc=BG_GOLD, ec="black", lw=1.2, radius=0.012)
    ax.text(bx + 0.015, by + bh / 2,
             "DIRECTIONAL  ASYMMETRY",
             ha="left", va="center",
             fontsize=7.5, fontweight="bold", color="white",
             bbox=dict(boxstyle="round,pad=0.20", facecolor="black",
                         edgecolor="none"))
    ax.text(bx + 0.150, by + bh * 0.65,
             "AF rewiring is biased toward RA → LA matrix invasion, "
             "not LA → RA homeostatic loss",
             ha="left", va="center",
             fontsize=7.0, fontweight="bold", color=GREY_DARK)
    ax.text(bx + 0.150, by + bh * 0.30,
             "21 RA→LA  vs  3 LA→RA  reversal genes (24 total)        "
             "binomtest p = 2.8 × 10⁻⁴   ★",
             ha="left", va="center",
             fontsize=6.5, color=BRBG_BROWN, fontweight="bold")

    # ==== 3-column layout ====
    # Left: RA chips (3); Center: tissue cartoon; Right: LA chips (3)
    left_x, left_w = 0.020, 0.215
    right_x, right_w = 0.765, 0.215
    cen_x, cen_w = 0.250, 0.500
    body_y = 0.165
    body_h = 0.680

    # RA-axis container header
    ax.text(left_x + left_w / 2, body_y + body_h - 0.008,
             "RA-AXIS  ·  ACUTE  ·  STRUCTURAL",
             ha="center", va="center",
             fontsize=7.2, fontweight="bold", color=RA_COLOR,
             bbox=dict(boxstyle="round,pad=0.22", facecolor=BG_RA,
                         edgecolor=RA_COLOR, linewidth=0.6))

    # LA-axis container header
    ax.text(right_x + right_w / 2, body_y + body_h - 0.008,
             "LA-AXIS  ·  CHRONIC  ·  METABOLIC",
             ha="center", va="center",
             fontsize=7.2, fontweight="bold", color=LA_COLOR,
             bbox=dict(boxstyle="round,pad=0.22", facecolor=BG_LA,
                         edgecolor=LA_COLOR, linewidth=0.6))

    # ===== LEFT chips (RA) =====
    chip_h = 0.180
    chip_pad = 0.014
    chip_top = body_y + body_h - 0.045
    # SM-1
    draw_evidence_chip(ax, left_x, chip_top - chip_h, left_w, chip_h,
                          "SM-1", "TGF-β / SLRP fibrosis",
                          color=RA_COLOR,
                          evidence_lines=[
                              "F10 bulk g = +0.68 ★ (Fig 2D)",
                              "F10 × fibro Δ = +0.010 p=0.026 ★",
                              "GAG × fibro ROI ρ = +0.43 ★ (Fig 5D)",
                              "SLRP enrich  OR=111.6  p=2.4e-4 ★",
                              "(N2d Test 2)  → DCN, PRELP",
                          ])
    # SM-2
    sm2_top = chip_top - chip_h - chip_pad
    draw_evidence_chip(ax, left_x, sm2_top - chip_h, left_w, chip_h,
                          "SM-2", "BMP  ·  bidirectional",
                          color=RA_COLOR,
                          evidence_lines=[
                              "F11 bulk g = +0.85 ★ (Fig 2D)",
                              "F11 × adipo  Δ = +0.032 ★",
                              "F11 × lymph  Δ = -0.068 ★",
                              "Δ_adipo > Δ_lymph 8/9 donors",
                              "Wilcoxon p = 0.008 ★ (N2d T3)",
                          ])
    # SM-3
    sm3_top = sm2_top - chip_h - chip_pad
    draw_evidence_chip(ax, left_x, sm3_top - chip_h, left_w, chip_h,
                          "SM-3", "HS-GAG  ·  diffuse",
                          color=RA_COLOR,
                          evidence_lines=[
                              "F9 bulk g = +0.79 ★ (Fig 2D)",
                              "Fig 5C class:  BULK-LED",
                              "no single LA cell type p<0.05",
                              "→ pan-tissue PG turnover",
                              "(no resolvable carrier)",
                          ])

    # ===== RIGHT chips (LA) =====
    # SM-4
    draw_evidence_chip(ax, right_x, chip_top - chip_h, right_w, chip_h,
                          "SM-4", "Mito Ca²⁺ buffering",
                          color=LA_COLOR,
                          evidence_lines=[
                              "F2 LAA ρ = -0.56 p=0.06 (Fig 4B)",
                              "Fig 5C class:  DURATION-LED",
                              "(only LA family classified)",
                              "→ chronic Ca²⁺ buffer loss",
                              "→ cytosolic Ca²⁺ overload",
                          ])
    # SM-5
    sm5_top = chip_top - chip_h - chip_pad
    draw_evidence_chip(ax, right_x, sm5_top - chip_h, right_w, chip_h,
                          "SM-5", "OXPHOS / TCA decline",
                          color=LA_COLOR,
                          evidence_lines=[
                              "F7 LAA  ρ = -0.58 p=0.06 (Fig 4)",
                              "F8 Complex I  trend (sub-thresh)",
                              "Fig 5D ROI:",
                              "  OXPHOS × atrialCM",
                              "  ρ = +0.45  p = 0.028 ★",
                          ])
    # SM-6/7
    sm67_top = sm5_top - chip_h - chip_pad
    draw_evidence_chip(ax, right_x, sm67_top - chip_h, right_w, chip_h,
                          "SM-6/7", "LA dediff (small)",
                          color=LA_COLOR,
                          evidence_lines=[
                              "3 LA→RA reversal genes (Fig 6B)",
                              "  ACTA1  sarcomeric actin loss",
                              "  ZFP57  imprinting loss",
                              "  PCDH20 adhesion loss",
                              "→ contractile + epigenetic dediff",
                          ])

    # ===== CENTER tissue cartoon =====
    draw_tissue_cartoon(ax, cen_x, body_y, cen_w, body_h)

    # ==== bottom synthesis ====
    fig.text(0.022, 0.130,
              "Synthesis: AF atrial remodeling decomposes into a structural "
              "RA-axis (3 RA-biased families:  TGF-β / SLRP fibrotic invasion, "
              "BMP context-dependent bidirectional signalling, HS-GAG diffuse "
              "turnover) acutely reversing chamber identity in LA tissue, "
              "and a chronic LA-axis (mitochondrial Ca²⁺ buffering decay + "
              "bioenergetic decline + small-but-consistent sarcomeric "
              "dedifferentiation) in atrialCM that intensifies with AF duration.",
              fontsize=6.7, color=GREY_DARK, ha="left", va="top",
              fontstyle="italic", wrap=True)
    fig.text(0.022, 0.038,
              "Cell-icon legend (centre cartoon): ",
              fontsize=6.0, color=GREY_DARK, fontweight="bold",
              ha="left", va="bottom")
    fig.text(0.150, 0.038,
              "fibroblast (green spindle)  ·  atrialCM (yellow striated)  ·  "
              "adipocyte (brown circle + droplet)  ·  lymph endo (cyan tube)",
              fontsize=5.7, color=GREY_DARK, ha="left", va="bottom")
    fig.text(0.022, 0.018,
              "All 7 sub-mechanisms carry explicit per-test p-values "
              "(★ marks p < 0.05).   Source: Fig 2D + 3B/3C + 4B + 5C/5D + "
              "6B/6C + N2d (results/n2d_mechanism_statistics.tsv).",
              fontsize=5.7, color="#444", ha="left", va="bottom",
              fontstyle="italic")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
