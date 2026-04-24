#!/usr/bin/env python3
"""
Figure 6A  Chamber-identity reversal gene discovery primitive (schematic).

Five-stage flow (mirrors Fig 3A / 4A / 5A layout):
    1. Per-gene aCAI primitive   (gene-level analogue of the Fig 2A
                                    program-level aCAI)
    2. Three paired AF cohorts (GSE128188 / GSE79768 / GSE115574) +
       8-donor healthy reference (GSE112339)
    3. Per-gene per-cohort calls:
            af_reversed       =  Healthy_dir × AF_dir < 0
            sr_preserved      =  Healthy_dir × SR_dir > 0
            af_specific_flip  =  af_reversed AND sr_preserved
            above_magnitude   =  |aCAI_healthy| ≥ 0.3 AND |aCAI_AF| ≥ 0.15
    4. Final reversal-gene call:
            af_specific_flip  AND  above_magnitude  in ≥ 2 of 3 cohorts
       → 24 reversal genes  (21 RA→LA  +  3 LA→RA)
    5. Cellular carrier attribution on GSE255612 LA snRNA → feeds 6B–6E.

Output:  figures/fig6a_schematic.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FIG_OUT = ROOT / "figures" / "fig6a_schematic.pdf"

LA_COLOR  = "#D6604D"
RA_COLOR  = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"
BG_GREY    = "#FAFAFA"
BG_GOLD    = "#FFFCF4"
BG_GREEN   = "#F2F8F4"
HL_GREEN   = "#3C8C57"
BRBG_BROWN = "#8C510A"
BRBG_TEAL  = "#01665E"


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def add_box(ax, x, y, w, h, text="", *, fc=BG_GREY, ec=GREY_DARK,
              text_color="black", fontsize=7.0, fontweight="normal",
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


def add_stage_header(ax, x, y, idx, title):
    ax.text(x, y, f"STAGE {idx}",
             fontsize=6.5, fontweight="bold",
             color="white", ha="left", va="center",
             bbox=dict(boxstyle="round,pad=0.20", facecolor=GREY_DARK,
                         edgecolor="none"),
             zorder=4)
    ax.text(x + 0.052, y, title,
             fontsize=7.5, fontweight="bold",
             color=GREY_DARK, ha="left", va="center")


def add_arrow(ax, x0, y0, x1, y1, *, color=GREY_DARK, lw=1.0):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                  arrowprops=dict(arrowstyle="-|>", color=color,
                                    lw=lw, mutation_scale=10), zorder=4)


def main() -> None:
    set_publication_style()
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(200 * MM, 130 * MM))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # title
    fig.text(0.025, 0.965,
              "Fig 6A  Chamber-identity reversal gene discovery primitive",
              fontsize=9.0, fontweight="bold", ha="left", va="center")
    fig.text(0.025, 0.938,
              "Per-gene aCAI applied to 3 paired AF cohorts + 8-donor "
              "healthy reference detects genes whose chamber identity "
              "flips in AF while SR donors preserve it.",
              fontsize=7.5, color=GREY_DARK, ha="left", va="center",
              style="italic")

    # column anchors (top row stages 1-4)
    s1_x, s1_w = 0.030, 0.190
    s2_x, s2_w = 0.245, 0.215
    s3_x, s3_w = 0.485, 0.215
    s4_x, s4_w = 0.725, 0.250
    row1_y_top = 0.880

    add_stage_header(ax, s1_x, row1_y_top, 1, "Per-gene aCAI primitive")
    add_stage_header(ax, s2_x, row1_y_top, 2, "Cohorts + healthy reference")
    add_stage_header(ax, s3_x, row1_y_top, 3, "Per-cohort gene calls")
    add_stage_header(ax, s4_x, row1_y_top, 4, "Reversal gene shortlist")

    # cards
    card_y, card_h = 0.510, 0.340
    for cx, cw, fc, ec, lw in [
        (s1_x, s1_w, "#FCFCFC", GREY_MID, 0.6),
        (s2_x, s2_w, "#FCFCFC", GREY_MID, 0.6),
        (s3_x, s3_w, "#FCFCFC", GREY_MID, 0.6),
        (s4_x, s4_w, "#FFFCF4", "black", 1.1),
    ]:
        add_box(ax, cx - 0.005, card_y, cw + 0.010, card_h, "",
                  fc=fc, ec=ec, lw=lw, radius=0.012)

    # ---- STAGE 1: per-gene aCAI ----
    ax.text(s1_x + s1_w / 2, 0.820,
             "Per-gene aCAI",
             fontsize=8.0, fontweight="bold", ha="center")
    ax.text(s1_x + s1_w / 2, 0.793,
             "(gene-level chamber-asymmetry score)",
             fontsize=6.5, color=GREY_DARK, ha="center", style="italic")
    ax.text(s1_x + s1_w / 2, 0.748,
             r"$\mathrm{aCAI}_g(d)\ =\ "
             r"\log_2 \mathrm{score\_LA}_{g,d}\ -\ "
             r"\log_2 \mathrm{score\_RA}_{g,d}$",
             fontsize=8.0, ha="center", va="center")
    ax.text(s1_x + s1_w / 2, 0.700,
             "donor d, gene g",
             fontsize=6.0, ha="center", color=GREY_DARK, style="italic")
    # quick legend bar
    ax.text(s1_x + s1_w / 2, 0.665,
             "+ → LA-biased    0 → symmetric    − → RA-biased",
             fontsize=6.0, ha="center", color=GREY_DARK)
    ax.text(s1_x + s1_w / 2, 0.610,
             "Reduced from Fig 2A primitive to single-gene granularity",
             fontsize=6.0, ha="center", color=GREY_DARK, style="italic")
    ax.text(s1_x + s1_w / 2, 0.560,
             "limma cannot ask the chamber-direction question;",
             fontsize=5.8, ha="center", color=GREY_MID, style="italic")
    ax.text(s1_x + s1_w / 2, 0.535,
             "GSEA exposes program-level NES, not per-gene sign-flip.",
             fontsize=5.8, ha="center", color=GREY_MID, style="italic")

    # ---- STAGE 2: cohorts ----
    ax.text(s2_x + s2_w / 2, 0.820,
             "3 paired AF cohorts",
             fontsize=8.0, fontweight="bold", ha="center")
    ax.text(s2_x + s2_w / 2, 0.795,
             "+ 8-donor bulk healthy reference",
             fontsize=6.4, ha="center", color=GREY_DARK, style="italic")
    cohort_labels = [
        ("GSE128188",  "5 AF + 5 SR"),
        ("GSE79768",   "7 AF + 6 SR"),
        ("GSE115574", "13 AF + 15 SR"),
    ]
    base_y = 0.745
    for i, (lbl, det) in enumerate(cohort_labels):
        y = base_y - i * 0.038
        ax.text(s2_x + 0.020, y, lbl,
                 fontsize=6.5, fontweight="bold", ha="left", va="center",
                 color=GREY_DARK)
        ax.text(s2_x + 0.118, y, det,
                 fontsize=6.2, ha="left", va="center", color=GREY_DARK)
    # healthy ref box
    add_box(ax, s2_x + 0.020, 0.598, s2_w - 0.040, 0.040,
              "GSE112339 healthy  ·  8 paired LA / RA donors",
              fc="#FFE9B0", ec="#D4A35E",
              fontsize=6.5, fontweight="bold", radius=0.008)
    ax.text(s2_x + s2_w / 2, 0.555,
             "all 4 datasets in gene-symbol space",
             fontsize=6.0, ha="center", color=GREY_MID, style="italic")
    ax.text(s2_x + s2_w / 2, 0.530,
             "(score_LA, score_RA per donor per gene)",
             fontsize=5.8, ha="center", color=GREY_MID, style="italic")

    add_arrow(ax, s1_x + s1_w + 0.005, 0.700,
                s2_x - 0.005, 0.700,
                color=GREY_DARK, lw=0.9)

    # ---- STAGE 3: per-cohort calls ----
    ax.text(s3_x + s3_w / 2, 0.820,
             "Per-cohort flip + magnitude",
             fontsize=7.5, fontweight="bold", ha="center")
    rules = [
        "af_reversed",         "Healthy_dir × AF_dir  <  0",
        "sr_preserved",        "Healthy_dir × SR_dir  >  0",
        "af_specific_flip",    "af_reversed  AND  sr_preserved",
        "above_magnitude",     "|aCAI_healthy| ≥ 0.30  AND  |aCAI_AF| ≥ 0.15",
    ]
    base_y = 0.770
    for k in range(0, len(rules), 2):
        lbl = rules[k]; cond = rules[k + 1]
        y = base_y - (k // 2) * 0.040
        ax.text(s3_x + 0.015, y, lbl,
                 fontsize=6.5, fontweight="bold", ha="left", va="center",
                 color=GREY_DARK)
        ax.text(s3_x + 0.105, y, cond,
                 fontsize=5.8, ha="left", va="center", color=GREY_DARK,
                 style="italic", family="monospace")
    ax.text(s3_x + s3_w / 2, 0.578,
             "applied independently per cohort",
             fontsize=6.0, ha="center", color=GREY_MID, style="italic")
    ax.text(s3_x + s3_w / 2, 0.548,
             "yields per-cohort flag set per gene",
             fontsize=6.0, ha="center", color=GREY_MID, style="italic")

    add_arrow(ax, s2_x + s2_w + 0.005, 0.700,
                s3_x - 0.005, 0.700,
                color=GREY_DARK, lw=0.9)

    # ---- STAGE 4: shortlist (golden) ----
    ax.text(s4_x + s4_w / 2, 0.820,
             "Reversal gene shortlist",
             fontsize=8.0, fontweight="bold", ha="center")
    ax.text(s4_x + s4_w / 2, 0.793,
             "af_specific_flip  AND  above_magnitude",
             fontsize=6.5, ha="center", color=GREY_DARK, style="italic")
    ax.text(s4_x + s4_w / 2, 0.770,
             "in  ≥  2  of  3  cohorts",
             fontsize=6.5, ha="center", color=GREY_DARK, style="italic")
    # The shortlist count — large, bold
    ax.text(s4_x + s4_w / 2, 0.715,
             "n  =  24  genes",
             fontsize=10.0, fontweight="bold", ha="center", va="center",
             color="black")
    # split
    ax.text(s4_x + s4_w / 2, 0.665,
             "21  RA → LA",
             fontsize=7.5, ha="center", va="center",
             color=RA_COLOR, fontweight="bold")
    ax.text(s4_x + s4_w / 2, 0.640,
             "3  LA → RA",
             fontsize=7.5, ha="center", va="center",
             color=LA_COLOR, fontweight="bold")
    ax.text(s4_x + s4_w / 2, 0.598,
             "matrix-invasion bias  (RA → LA dominant)",
             fontsize=6.0, ha="center", color=GREY_DARK, style="italic")
    ax.text(s4_x + s4_w / 2, 0.555,
             "framework-discovered;  unique to CRF",
             fontsize=6.0, ha="center", color=HL_GREEN, fontweight="bold")
    ax.text(s4_x + s4_w / 2, 0.530,
             "(not visible to limma / GSEA)",
             fontsize=5.8, ha="center", color=GREY_MID, style="italic")

    add_arrow(ax, s3_x + s3_w + 0.005, 0.700,
                s4_x - 0.005, 0.700,
                color=GREY_DARK, lw=0.9)

    # ===================================================================
    # ROW 2  STAGE 5  cellular carrier + downstream panels
    # ===================================================================
    add_stage_header(ax, 0.030, 0.470, 5,
                       "Cellular carrier attribution  ·  feeds Fig 6B–6E")

    # Mini reversal-gene heatmap (program × cohort thumb)
    mat_x, mat_y = 0.030, 0.130
    mat_w, mat_h = 0.220, 0.330
    add_box(ax, mat_x, mat_y, mat_w, mat_h, "",
              fc="#FFFFFF", ec=GREY_DARK, lw=0.8, radius=0.012)
    ax.text(mat_x + mat_w / 2, mat_y + mat_h - 0.024,
             "Reversal gene table",
             fontsize=7.5, fontweight="bold", ha="center")
    ax.text(mat_x + mat_w / 2, mat_y + mat_h - 0.048,
             "(24 genes × 7 columns:  healthy + 3 SR + 3 AF)",
             fontsize=5.8, color=GREY_DARK, ha="center", style="italic")
    # mini fake heatmap
    rng = np.random.default_rng(13)
    nrows, ncols = 14, 7
    fake = rng.normal(scale=0.4, size=(nrows, ncols))
    # healthy + SR cols : RA-biased rows negative
    fake[:11, 0:4] -= 0.4
    # AF cols : flipped to positive (chamber-reversal)
    fake[:11, 4:] += 0.6
    # last 3 rows: opposite (LA→RA)
    fake[11:, 0:4] += 0.5
    fake[11:, 4:]  -= 0.5
    cmap_fake = mpl.colors.LinearSegmentedColormap.from_list(
        "h", [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)])
    norm_fake = mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    mini_x0 = mat_x + 0.020; mini_y0 = mat_y + 0.025
    mini_w  = mat_w - 0.040; mini_h  = mat_h - 0.105
    cw = mini_w / ncols; ch = mini_h / nrows
    for r in range(nrows):
        for c in range(ncols):
            ax.add_patch(mpatches.Rectangle(
                (mini_x0 + c * cw, mini_y0 + r * ch),
                cw, ch,
                facecolor=cmap_fake(norm_fake(fake[r, c])),
                edgecolor="white", linewidth=0.2, zorder=3))
    # vertical block separators in mini
    for vx in (mini_x0 + 1 * cw, mini_x0 + 4 * cw):
        ax.plot([vx, vx], [mini_y0, mini_y0 + mini_h],
                  color=GREY_DARK, lw=0.5, zorder=4)

    add_arrow(ax, s4_x + s4_w / 2, 0.508,
                mat_x + mat_w / 2, mat_y + mat_h + 0.005,
                color=GREY_DARK, lw=1.0)

    # Carrier attribution box
    car_x, car_y, car_w, car_h = 0.290, 0.130, 0.290, 0.330
    add_box(ax, car_x, car_y, car_w, car_h, "",
              fc=BG_GREY, ec=GREY_DARK, lw=0.8, radius=0.012)
    ax.text(car_x + car_w / 2, car_y + car_h - 0.024,
             "Cellular carrier (GSE255612)",
             fontsize=7.5, fontweight="bold", ha="center")
    ax.text(car_x + car_w / 2, car_y + car_h - 0.054,
             r"$\Delta_{g,c}\ =\ "
             r"\overline{\mathrm{expr}}^{\mathrm{AF}}_{g,c}\ -\ "
             r"\overline{\mathrm{expr}}^{\mathrm{Ctrl}}_{g,c}$",
             fontsize=8.0, ha="center", va="center")
    ax.text(car_x + car_w / 2, car_y + car_h - 0.105,
             "→ sign-aligned by gene's expected reversal direction",
             fontsize=6.2, ha="center", va="center", color=GREY_DARK)
    ax.text(car_x + car_w / 2, car_y + car_h - 0.135,
             "→ argmax over 11 cell types  =  dominant carrier",
             fontsize=6.2, ha="center", va="center", color=GREY_DARK,
             style="italic")
    ax.text(car_x + car_w / 2, car_y + 0.078,
             "fibroblast  Δ_aligned  =  +0.042  (top carrier)",
             fontsize=6.5, ha="center", va="center",
             color="#1F5C39", fontweight="bold")
    ax.text(car_x + car_w / 2, car_y + 0.050,
             "macrophage  +0.024     mast cell  +0.039",
             fontsize=6.0, ha="center", va="center",
             color=GREY_DARK)
    ax.text(car_x + car_w / 2, car_y + 0.025,
             "(across 24 reversal genes, mean Δ_aligned per cell type)",
             fontsize=5.8, ha="center", va="center",
             color=GREY_MID, style="italic")

    add_arrow(ax, mat_x + mat_w + 0.005, car_y + car_h / 2,
                car_x - 0.005, car_y + car_h / 2,
                color=GREY_DARK, lw=1.0)

    # Downstream panels (right green box)
    pan_x, pan_y, pan_w, pan_h = 0.620, 0.130, 0.360, 0.330
    add_box(ax, pan_x, pan_y, pan_w, pan_h, "",
              fc=BG_GREEN, ec=HL_GREEN, lw=1.0, radius=0.012)
    ax.text(pan_x + pan_w / 2, pan_y + pan_h - 0.024,
             "Downstream Fig 6 panels",
             fontsize=7.5, fontweight="bold", ha="center", color="#1F5C39")
    bullets = [
        ("Fig 6B",
         "24 reversal genes × 3 cohorts heatmap (Healthy / SR / AF)"),
        ("Fig 6C",
         "Cellular carrier matrix (24 genes × 11 cell types)"),
        ("Fig 6D",
         "Top anchor genes' per-cohort aCAI profiles"),
        ("Fig 6E",
         "Classical mRNA cis-eQTL × GWAS regime check"),
    ]
    for k, (tag, desc) in enumerate(bullets):
        y = pan_y + pan_h - 0.075 - k * 0.058
        ax.text(pan_x + 0.020, y, tag,
                 fontsize=7.0, fontweight="bold",
                 color="#1F5C39", ha="left", va="center")
        ax.text(pan_x + 0.075, y, desc,
                 fontsize=6.4, color=GREY_DARK, ha="left", va="center")

    add_arrow(ax, car_x + car_w + 0.005, car_y + car_h / 2,
                pan_x - 0.005, car_y + car_h / 2,
                color=GREY_DARK, lw=1.0)

    # ---- bottom notation ----
    fig.text(0.025, 0.030,
              "Notation:  d = donor,  g = gene,  c = cell type;  "
              "Healthy_dir from GSE112339 (8-donor bulk reference);  "
              "AF_dir / SR_dir from each paired AF cohort.",
              fontsize=6.2, color="#444", ha="left", va="bottom")
    fig.text(0.025, 0.010,
              "The reversal-gene primitive is uniquely enabled by per-gene "
              "donor-level chamber-asymmetry — neither limma DE (no "
              "chamber direction) nor pathway GSEA (no per-gene flip) "
              "can produce this object class.",
              fontsize=6.2, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
