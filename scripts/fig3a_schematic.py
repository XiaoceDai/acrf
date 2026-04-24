#!/usr/bin/env python3
"""
Figure 3A  Single-cell carrier-attribution analytical schematic.

Five-stage flow:
    1. NCP_v2 panel        (30 programs in 11 axis families, from Fig 1c-e)
    2. GSE255612 LA snRNA  (34 donors: 18 AF, 16 Ctrl;  13 cell types
       passed the n ≥ 30 cells gate; 11 used here, B_cell + epicardial
       excluded for n_Ctrl < 4)
    3. Donor × cell-type pseudobulk  (mean log-norm expression per gene)
    4. Per (donor, cell type, program) singscore on the program's
       leading-edge genes  →  Δ singscore = mean_AF − mean_Ctrl
    5. Family-level  &  program-level  Δ + permutation p_perm
       feeding Figs 3B / 3C / 3D / 3E

Output:  figures/fig3a_schematic.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FIG_OUT = ROOT / "figures" / "fig3a_schematic.pdf"

LA_COLOR  = "#D6604D"
RA_COLOR  = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"
BG_GREY    = "#FAFAFA"
BG_AF      = "#444444"
BG_CTRL    = "#FFFFFF"
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


def add_stage_header(ax, x, y, idx, title):
    ax.text(x, y,
             f"STAGE {idx}",
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
                                    lw=lw, mutation_scale=10),
                  zorder=4)


def main() -> None:
    set_publication_style()
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(200 * MM, 130 * MM))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # ---- title ----
    fig.text(0.025, 0.965,
              "Fig 3A  Cell-type carrier attribution of NCP_v2 rewiring "
              "in the LA snRNA arm",
              fontsize=9.0, fontweight="bold", ha="left", va="center")
    fig.text(0.025, 0.938,
              "30 NCP_v2 programs (Fig 1c–e) → per-donor per-cell-type "
              "singscore → Δ singscore (AF − Ctrl) per program × cell type "
              "→ family-level aggregation.",
              fontsize=7.5, color=GREY_DARK, ha="left", va="center",
              style="italic")

    # ===================================================================
    # ROW 1   STAGE 1 → STAGE 4   (per-donor pipeline)
    # ===================================================================
    s1_x, s1_w = 0.030, 0.180
    s2_x, s2_w = 0.235, 0.215
    s3_x, s3_w = 0.475, 0.215
    s4_x, s4_w = 0.715, 0.260
    row1_y_top = 0.880

    add_stage_header(ax, s1_x, row1_y_top, 1, "NCP_v2 panel")
    add_stage_header(ax, s2_x, row1_y_top, 2, "GSE255612 LA snRNA")
    add_stage_header(ax, s3_x, row1_y_top, 3, "Pseudobulk + singscore")
    add_stage_header(ax, s4_x, row1_y_top, 4, "Δ singscore (AF − Ctrl)")

    # ---- Stage cards (taller, with explicit padding so all text sits
    # comfortably inside) ----
    card_y, card_h = 0.510, 0.340
    for cx, cw, fc, ec, ec_lw in [
        (s1_x, s1_w, "#FCFCFC", GREY_MID, 0.6),
        (s2_x, s2_w, "#FCFCFC", GREY_MID, 0.6),
        (s3_x, s3_w, "#FCFCFC", GREY_MID, 0.6),
        (s4_x, s4_w, "#FFFCF4", "black", 1.1),
    ]:
        add_box(ax, cx - 0.005, card_y, cw + 0.010, card_h,
                  "", fc=fc, ec=ec, lw=ec_lw, radius=0.012)

    # Inner content y-bands (top → bottom).
    #   row A : title          y = 0.820
    #   row B : sub-text       y = 0.795
    #   row C : main artwork   y = 0.620–0.770
    #   row D : caption        y = 0.555  (well inside card; card bottom = 0.510)

    # ---- STAGE 1  panel ----
    ax.text(s1_x + s1_w / 2, 0.820,
             "30 NCP_v2 programs",
             fontsize=8.0, fontweight="bold", ha="center", va="center")
    ax.text(s1_x + s1_w / 2, 0.795,
             "in 11 axis families",
             fontsize=7.0, color=GREY_DARK, ha="center", va="center")
    swatch_pal = [
        "#4F86C6", "#6FA8DC", "#9EC4E0",
        "#82B5A1", "#A2C8B0", "#B7D2BC",
        "#D4B58E", "#C8C26C",
        "#E89C8B", "#D86A6A", "#B94B4B",
    ]
    sw, sh = 0.020, 0.024
    rows, cols = 3, 4
    base_x = s1_x + (s1_w - cols * sw - (cols - 1) * 0.005) / 2
    base_y = 0.660
    for k in range(11):
        rr, cc = divmod(k, cols)
        x = base_x + cc * (sw + 0.005)
        y = base_y + (rows - 1 - rr) * (sh + 0.005)
        ax.add_patch(mpatches.Rectangle(
            (x, y), sw, sh, facecolor=swatch_pal[k],
            edgecolor="white", linewidth=0.4, zorder=3))
        ax.text(x + sw / 2, y + sh / 2, f"F{k+1}",
                 ha="center", va="center", fontsize=5.5,
                 fontweight="bold",
                 color="white" if k in (0, 8, 9, 10) else "black",
                 zorder=4)
    ax.text(s1_x + s1_w / 2, 0.555,
             "leading-edge genes from",
             fontsize=6.0, ha="center", va="center", color=GREY_DARK,
             style="italic")
    ax.text(s1_x + s1_w / 2, 0.535,
             "fgsea on bulk LA-vs-RA rank",
             fontsize=6.0, ha="center", va="center", color=GREY_DARK,
             style="italic")

    # ---- STAGE 2  donors / cell types ----
    ax.text(s2_x + s2_w / 2, 0.820,
             "34 donors  ·  13 cell types",
             fontsize=8.0, fontweight="bold", ha="center")
    af_x = s2_x + 0.013; ctl_x = s2_x + 0.115
    ax.add_patch(mpatches.Rectangle(
        (af_x, 0.760), 0.085, 0.030,
        facecolor=BG_AF, edgecolor=GREY_MID, linewidth=0.4, zorder=3))
    ax.text(af_x + 0.0425, 0.775, "AF  n = 18",
             ha="center", va="center", color="white",
             fontsize=6.5, fontweight="bold", zorder=4)
    ax.add_patch(mpatches.Rectangle(
        (ctl_x, 0.760), 0.085, 0.030,
        facecolor=BG_CTRL, edgecolor=GREY_MID, linewidth=0.4, zorder=3))
    ax.text(ctl_x + 0.0425, 0.775, "Ctrl  n = 16",
             ha="center", va="center", color="black",
             fontsize=6.5, fontweight="bold", zorder=4)

    ct_labels = ["atrialCM", "fibro", "endo", "lymph endo",
                  "pericyte", "VSMC", "adipo",
                  "macrophage", "lymphocyte", "mast", "neuronal"]
    base_x2 = s2_x + 0.014
    base_y2 = 0.715
    chip_w, chip_h = 0.044, 0.026
    for k, lbl in enumerate(ct_labels):
        rr, cc = divmod(k, 4)
        x = base_x2 + cc * (chip_w + 0.006)
        y = base_y2 - rr * (chip_h + 0.008)
        ax.add_patch(mpatches.FancyBboxPatch(
            (x, y), chip_w, chip_h,
            boxstyle="round,pad=0.002,rounding_size=0.005",
            facecolor=GREY_LIGHT, edgecolor=GREY_MID, lw=0.4,
            zorder=3))
        ax.text(x + chip_w / 2, y + chip_h / 2, lbl,
                 ha="center", va="center", fontsize=5.4,
                 color=GREY_DARK, zorder=4)
    ax.text(s2_x + s2_w / 2, 0.545,
             "B_cell + epicardial excluded   (n_Ctrl < 4)",
             fontsize=6.0, ha="center", va="center",
             color=GREY_MID, style="italic")

    add_arrow(ax, s1_x + s1_w + 0.005, 0.690, s2_x - 0.005, 0.690,
                color=GREY_DARK, lw=0.9)

    # ---- STAGE 3  pseudobulk + singscore ----
    ax.text(s3_x + s3_w / 2, 0.820,
             "Donor × cell-type pseudobulk",
             fontsize=7.5, fontweight="bold", ha="center")
    ax.text(s3_x + s3_w / 2, 0.793,
             "mean log-norm expression  per gene",
             fontsize=6.6, color=GREY_DARK, ha="center")
    add_box(ax, s3_x + 0.020, 0.720, s3_w - 0.040, 0.045,
              "rankGenes  →  simpleScore (LE genes)",
              fc=BG_GREY, ec=GREY_DARK, fontsize=7.0, fontweight="bold",
              radius=0.008)
    ax.text(s3_x + s3_w / 2, 0.670,
             r"singscore$_{i,\,k,\,c}$",
             fontsize=9.0, ha="center", va="center")
    ax.text(s3_x + s3_w / 2, 0.625,
             "donor i,  program k,  cell type c",
             fontsize=6.2, ha="center", va="center", color=GREY_DARK,
             style="italic")
    ax.text(s3_x + s3_w / 2, 0.583,
             "9,090 (donor × ct × prog) scores",
             fontsize=6.0, ha="center", va="center", color=GREY_DARK)
    ax.text(s3_x + s3_w / 2, 0.545,
             "across all 11 included cell types",
             fontsize=5.8, ha="center", va="center", color=GREY_MID,
             style="italic")

    add_arrow(ax, s2_x + s2_w + 0.005, 0.690, s3_x - 0.005, 0.690,
                color=GREY_DARK, lw=0.9)

    # ---- STAGE 4  Δ formula card (golden) ----
    ax.text(s4_x + s4_w / 2, 0.820,
             "Per-program × cell-type effect",
             fontsize=8.0, fontweight="bold", ha="center")
    ax.text(s4_x + s4_w / 2, 0.760,
             r"$\Delta_{k,c}\ =\ "
             r"\overline{\mathrm{singscore}}^{\mathrm{AF}}_{k,c}\ -\ "
             r"\overline{\mathrm{singscore}}^{\mathrm{Ctrl}}_{k,c}$",
             fontsize=9.0, ha="center", va="center")
    ax.text(s4_x + s4_w / 2, 0.700,
             "+ Δ  →  AF higher",
             fontsize=7.0, color=BRBG_BROWN, ha="center", va="center",
             fontweight="bold")
    ax.text(s4_x + s4_w / 2, 0.674,
             "0  →  no change",
             fontsize=7.0, color=GREY_DARK, ha="center", va="center")
    ax.text(s4_x + s4_w / 2, 0.648,
             "− Δ  →  AF lower",
             fontsize=7.0, color=BRBG_TEAL, ha="center", va="center",
             fontweight="bold")
    ax.text(s4_x + s4_w / 2, 0.598,
             "10 000 MC perm  →  p_perm",
             fontsize=6.5, ha="center", va="center", color=GREY_DARK,
             style="italic")
    ax.text(s4_x + s4_w / 2, 0.560,
             "BH q across programs / families",
             fontsize=6.2, ha="center", va="center", color=GREY_DARK,
             style="italic")
    ax.text(s4_x + s4_w / 2, 0.535,
             "(separate test per cell type)",
             fontsize=5.8, ha="center", va="center", color=GREY_MID,
             style="italic")

    add_arrow(ax, s3_x + s3_w + 0.005, 0.690, s4_x - 0.005, 0.690,
                color=GREY_DARK, lw=0.9)

    # ===================================================================
    # ROW 2  STAGE 5  family aggregation + downstream panels
    # ===================================================================
    add_stage_header(ax, 0.030, 0.470, 5,
                       "Family aggregation  ·  feeds Fig 3B / 3C / 3D / 3E")

    # Mini Δ heatmap thumbnail (representing M_dScore)
    mat_x, mat_y = 0.030, 0.130
    mat_w, mat_h = 0.220, 0.330
    add_box(ax, mat_x, mat_y, mat_w, mat_h, "",
              fc="#FFFFFF", ec=GREY_DARK, lw=0.8, radius=0.012)
    ax.text(mat_x + mat_w / 2, mat_y + mat_h - 0.024,
             "M_Δ singscore", fontsize=7.5, fontweight="bold",
             ha="center")
    ax.text(mat_x + mat_w / 2, mat_y + mat_h - 0.048,
             "(30 programs × 11 cell types)",
             fontsize=6.0, color=GREY_DARK, ha="center", style="italic")
    # tiny BrBG fake heatmap
    rng = np.random.default_rng(11)
    fake = rng.normal(scale=0.5, size=(15, 11))
    fake[12:, :] += 0.5  # RA-biased rows positive shift
    fake[:3, :] -= 0.3   # LA-biased rows attenuated
    cmap_fake = mpl.colors.LinearSegmentedColormap.from_list(
        "h", [(0.0, BRBG_TEAL), (0.5, "#FFFFFF"), (1.0, BRBG_BROWN)])
    norm_fake = mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    mini_x0 = mat_x + 0.015; mini_y0 = mat_y + 0.020
    mini_w  = mat_w - 0.030; mini_h  = mat_h - 0.090
    nrows, ncols = fake.shape
    cw = mini_w / ncols; ch = mini_h / nrows
    for r in range(nrows):
        for c in range(ncols):
            ax.add_patch(mpatches.Rectangle(
                (mini_x0 + c * cw, mini_y0 + r * ch),
                cw, ch,
                facecolor=cmap_fake(norm_fake(fake[r, c])),
                edgecolor="none", zorder=3))

    add_arrow(ax, s4_x + s4_w / 2, 0.508,
                mat_x + mat_w / 2, mat_y + mat_h + 0.005,
                color=GREY_DARK, lw=1.0)

    # Family aggregation card
    fam_x, fam_y, fam_w, fam_h = 0.290, 0.130, 0.290, 0.330
    add_box(ax, fam_x, fam_y, fam_w, fam_h, "",
              fc=BG_GREY, ec=GREY_DARK, lw=0.8, radius=0.012)
    ax.text(fam_x + fam_w / 2, fam_y + fam_h - 0.024,
             "Family aggregation",
             fontsize=7.5, fontweight="bold", ha="center")
    ax.text(fam_x + fam_w / 2, fam_y + fam_h - 0.054,
             r"$\Delta^{\mathrm{family}}_{F,c}\ =\ "
             r"\frac{1}{|\mathrm{F}|} \sum_{k\in F} \Delta_{k,c}$",
             fontsize=8.5, ha="center", va="center")
    ax.text(fam_x + fam_w / 2, fam_y + fam_h - 0.105,
             "→ 11 axis families × 11 cell types  =  121 cells",
             fontsize=6.5, ha="center", va="center", color=GREY_DARK)
    ax.text(fam_x + fam_w / 2, fam_y + fam_h - 0.140,
             "Independent perm test on family-mean singscore",
             fontsize=6.2, ha="center", va="center", color=GREY_DARK,
             style="italic")
    ax.text(fam_x + fam_w / 2, fam_y + 0.030,
             "bridges Fig 2D forest  →  Fig 3B carrier heatmap",
             fontsize=6.5, ha="center", va="center",
             color="#1F5C39", fontweight="bold")

    add_arrow(ax, mat_x + mat_w + 0.005, fam_y + fam_h / 2,
                fam_x - 0.005, fam_y + fam_h / 2,
                color=GREY_DARK, lw=1.0)

    # Output panels card
    pan_x, pan_y, pan_w, pan_h = 0.620, 0.130, 0.360, 0.330
    add_box(ax, pan_x, pan_y, pan_w, pan_h, "",
              fc=BG_GREEN, ec=HL_GREEN, lw=1.0, radius=0.012)
    ax.text(pan_x + pan_w / 2, pan_y + pan_h - 0.024,
             "Downstream Fig 3 panels",
             fontsize=7.5, fontweight="bold", ha="center", color="#1F5C39")
    bullets = [
        ("Fig 3B",
         "family × cell-type Δ heatmap (11 × 11)"),
        ("Fig 3C",
         "RA-biased family decomposition (F9 / F10 / F11)"),
        ("Fig 3D",
         "full program × cell-type heatmap (30 × 11)"),
        ("Fig 3E",
         "anchor genes  +  Fig 6c reversal-gene cross-mark"),
    ]
    for k, (tag, desc) in enumerate(bullets):
        y = pan_y + pan_h - 0.075 - k * 0.058
        ax.text(pan_x + 0.020, y, tag,
                 fontsize=7.0, fontweight="bold",
                 color="#1F5C39", ha="left", va="center")
        ax.text(pan_x + 0.075, y, desc,
                 fontsize=6.5, color=GREY_DARK, ha="left", va="center")

    add_arrow(ax, fam_x + fam_w + 0.005, fam_y + fam_h / 2,
                pan_x - 0.005, fam_y + fam_h / 2,
                color=GREY_DARK, lw=1.0)

    # ---- bottom notation ----
    fig.text(0.025, 0.030,
              "Notation:  i = donor;  k = NCP_v2 program;  "
              "c = cell type;  F = axis family;  LE = leading-edge gene set.",
              fontsize=6.2, color="#444", ha="left", va="bottom")
    fig.text(0.025, 0.010,
              "All inference is on Δ singscore between AF and Ctrl donor "
              "groups within each cell type, never on raw expression — "
              "this controls for cell-type composition and donor batch.",
              fontsize=6.2, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
