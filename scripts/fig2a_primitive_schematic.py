#!/usr/bin/env python3
"""Figure 2A — Bridge from Fig1 locked v3 panel to AF rewiring quantification."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FIG_OUT = ROOT / "figures" / "fig2a_primitive_schematic.pdf"

LA_COLOR  = "#C97A74"
RA_COLOR  = "#7D99B8"
GREY_LIGHT = "#E6E6E6"
GREY_MID  = "#9C9C9C"
GREY_DARK = "#444444"
BG_LA     = "#FBE7E2"
BG_RA     = "#E6F0F7"
BOX_BG    = "#FAFAFA"
HL_GREEN  = "#3C8C57"
BG_GREEN  = "#F2F8F4"


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "axes.titlesize":    9,
        "axes.labelsize":    8,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def add_box(ax, x, y, w, h, text="", *,
              fc=BOX_BG, ec=GREY_DARK, text_color="black",
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


def add_stage_header(ax, x_left, x_right, y, idx, title):
    ax.text(x_left, y,
             f"STAGE {idx}",
             fontsize=6.5, fontweight="bold",
             color="white", ha="left", va="center",
             bbox=dict(boxstyle="round,pad=0.20", facecolor=GREY_DARK,
                         edgecolor="none"),
             zorder=4)
    ax.text(x_left + 0.052, y, title,
             fontsize=7.5, fontweight="bold",
             color=GREY_DARK, ha="left", va="center")


def add_arrow(ax, x0, y0, x1, y1, *, color=GREY_DARK, lw=1.0,
                label=None, label_color=None,
                label_offset=(0, 0.012), label_fontsize=6.5):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                  arrowprops=dict(arrowstyle="-|>", color=color,
                                    lw=lw, mutation_scale=10),
                  zorder=4)
    if label:
        ax.text((x0 + x1) / 2 + label_offset[0],
                 (y0 + y1) / 2 + label_offset[1],
                 label,
                 fontsize=label_fontsize, fontweight="bold",
                 color=(label_color or color),
                 ha="center", va="center", zorder=5,
                 bbox=dict(boxstyle="round,pad=0.10",
                              facecolor="white", edgecolor="none"))


def main() -> None:
    set_publication_style()

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(200 * MM, 130 * MM))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    fig.text(0.025, 0.965,
              "Fig 2A  Fig1->Fig2 bridge: frozen NCP panel to AF rewiring unit",
              fontsize=9.0, fontweight="bold", ha="left", va="center")
    fig.text(0.025, 0.935,
              "No panel re-selection in Fig2: NCP programs are fixed from Fig1, "
              "then projected to donor-paired AF cohorts for rewiring estimation.",
              fontsize=7.5, color=GREY_DARK, ha="left", va="center",
              style="italic")

    # ROW 1  STAGE 1 → STAGE 4   (per-donor primitive)
    row1_top = 0.86; row1_bot = 0.51

    # Stage column anchors
    s1_x = 0.030;  s1_w = 0.180
    s2_x = 0.235;  s2_w = 0.180
    s3_x = 0.440;  s3_w = 0.220
    s4_x = 0.685;  s4_w = 0.290

    # Stage headers
    add_stage_header(ax, s1_x, s1_x + s1_w, row1_top,
                       1, "Per-donor paired tissue")
    add_stage_header(ax, s2_x, s2_x + s2_w, row1_top,
                       2, "Per-chamber expression")
    add_stage_header(ax, s3_x, s3_x + s3_w, row1_top,
                       3, "Frozen program object from Fig1")
    add_stage_header(ax, s4_x, s4_x + s4_w, row1_top,
                       4, "Donor-level rewiring primitive")

    card_y = 0.555; card_h = 0.260
    for cx, cw in [(s1_x, s1_w), (s2_x, s2_w),
                     (s3_x, s3_w), (s4_x, s4_w)]:
        add_box(ax, cx - 0.005, card_y, cw + 0.010, card_h,
                  "", fc="#FCFCFC", ec=GREY_LIGHT, lw=0.6,
                  radius=0.012)

    ax.text(s1_x + s1_w / 2, 0.785,
             "donor i",
             fontsize=7.5, fontweight="bold", ha="center", va="center")
    add_box(ax, s1_x + 0.012, 0.660, 0.073, 0.105,
              "LA", fc=LA_COLOR, ec=LA_COLOR,
              text_color="white", fontsize=8.5, fontweight="bold",
              radius=0.008)
    add_box(ax, s1_x + 0.095, 0.660, 0.073, 0.105,
              "RA", fc=RA_COLOR, ec=RA_COLOR,
              text_color="white", fontsize=8.5, fontweight="bold",
              radius=0.008)
    ax.text(s1_x + s1_w / 2, 0.605,
             "n = 51 paired donors",
             fontsize=6.8, color=GREY_DARK, ha="center", va="center",
             fontweight="bold")
    ax.text(s1_x + s1_w / 2, 0.580,
             "3 paired AF cohorts",
             fontsize=6.3, color=GREY_DARK, ha="center", va="center",
             style="italic")

    add_box(ax, s2_x, 0.700, s2_w, 0.080,
              "expr_LA  (gene-level)",
              fc=BG_LA, ec=LA_COLOR, fontsize=7.0, fontweight="bold")
    add_box(ax, s2_x, 0.605, s2_w, 0.080,
              "expr_RA  (gene-level)",
              fc=BG_RA, ec=RA_COLOR, fontsize=7.0, fontweight="bold")
    ax.text(s2_x + s2_w / 2, 0.575,
             "bulk RNA-seq  /  pseudobulk per chamber",
             fontsize=6.0, color=GREY_DARK, ha="center", va="center",
             style="italic")

    add_arrow(ax, s1_x + s1_w + 0.005, 0.740,
                s2_x - 0.002, 0.740,
                color=LA_COLOR, lw=0.9)
    add_arrow(ax, s1_x + s1_w + 0.005, 0.645,
                s2_x - 0.002, 0.645,
                color=RA_COLOR, lw=0.9)

    add_box(ax, s3_x, 0.700, s3_w, 0.080,
              r"$S_{ik}^{\mathrm{LA}}=\mathrm{NES}(P_k\mid\mathrm{LA}_i)$",
              fc=BG_LA, ec=LA_COLOR, fontsize=7.2)
    add_box(ax, s3_x, 0.605, s3_w, 0.080,
              r"$S_{ik}^{\mathrm{RA}}=\mathrm{NES}(P_k\mid\mathrm{RA}_i)$",
              fc=BG_RA, ec=RA_COLOR, fontsize=7.2)
    ax.text(s3_x + s3_w / 2, 0.575,
             "P_k from locked NCP panel (K = 30; defined in Fig1)",
             fontsize=6.0, color=GREY_DARK, ha="center", va="center",
             style="italic")

    add_arrow(ax, s2_x + s2_w + 0.005, 0.740,
                s3_x - 0.002, 0.740,
                color=LA_COLOR, lw=0.9)
    add_arrow(ax, s2_x + s2_w + 0.005, 0.645,
                s3_x - 0.002, 0.645,
                color=RA_COLOR, lw=0.9)

    add_box(ax, s4_x - 0.005, card_y, s4_w + 0.010, card_h,
              "", fc="#FFFCF4", ec="black", lw=1.2, radius=0.014)
    ax.text(s4_x + s4_w / 2, 0.785,
             "aCAI primitive",
             fontsize=8.0, fontweight="bold", ha="center", va="center")
    ax.text(s4_x + s4_w / 2, 0.734,
             r"$\mathrm{aCAI}_{ik}=S_{ik}^{\mathrm{LA}}-S_{ik}^{\mathrm{RA}}$",
             fontsize=9.2, ha="center", va="center")
    ax.text(s4_x + s4_w / 2, 0.704,
             r"$\mathrm{aCAI}_{ik}^{\mathrm{rel}}="
             r"\frac{\mathrm{aCAI}_{ik}}{\left|\mathrm{aCAI}_{k}^{\mathrm{healthy}}\right|}$",
             fontsize=8.2, ha="center", va="center")
    ax.text(s4_x + s4_w / 2, 0.668,
             "+ value  →  LA-biased",
             fontsize=7.0, color=LA_COLOR, ha="center", va="center",
             fontweight="bold")
    ax.text(s4_x + s4_w / 2, 0.640,
             "0  →  symmetric",
             fontsize=7.0, color=GREY_DARK, ha="center", va="center")
    ax.text(s4_x + s4_w / 2, 0.612,
             "− value  →  RA-biased",
             fontsize=7.0, color=RA_COLOR, ha="center", va="center",
             fontweight="bold")
    ax.text(s4_x + s4_w / 2, 0.578,
             "healthy-baseline normalized, comparable across cohorts",
             fontsize=6.0, color=GREY_DARK, ha="center", va="center",
             style="italic")

    # arrows from S3 boxes into S4 box (no +/− labels — redundant with
    # the equation, and they were colliding with the aCAI formula)
    add_arrow(ax, s3_x + s3_w + 0.005, 0.740,
                s4_x - 0.007, 0.722,
                color=LA_COLOR, lw=0.9)
    add_arrow(ax, s3_x + s3_w + 0.005, 0.645,
                s4_x - 0.007, 0.700,
                color=RA_COLOR, lw=0.9)

    # ROW 2  STAGE 5  (cohort-stratified inference)
    row2_top = 0.46
    add_stage_header(ax, 0.030, 0.985, row2_top,
                       5, "AF-stage inference on fixed NCP object")

    # M_aCAI matrix
    mat_x = 0.030; mat_y = 0.13; mat_w = 0.225; mat_h = 0.30
    add_box(ax, mat_x, mat_y, mat_w, mat_h, "",
              fc="#FFFFFF", ec=GREY_DARK, lw=0.8, radius=0.012)
    ax.text(mat_x + mat_w / 2, mat_y + mat_h - 0.024,
             "M_aCAI^rel  matrix",
             fontsize=7.5, fontweight="bold", ha="center")
    ax.text(mat_x + mat_w / 2, mat_y + mat_h - 0.048,
             "(30 programs  ×  51 cohort-scoped donors)",
             fontsize=6.0, color=GREY_DARK, ha="center", style="italic")
    # mini fake heatmap inside
    mini_x0 = mat_x + 0.015; mini_y0 = mat_y + 0.02
    mini_w  = mat_w - 0.030; mini_h  = mat_h - 0.085
    rng = np.random.default_rng(7)
    nrows, ncols = 14, 16
    cell_w = mini_w / ncols; cell_h = mini_h / nrows
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "h", [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)])
    norm = mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    fake = rng.normal(scale=0.45, size=(nrows, ncols))
    fake[0:5, :]  += 0.45
    fake[-3:, :]  -= 0.55
    for r in range(nrows):
        for c in range(ncols):
            ax.add_patch(mpatches.Rectangle(
                (mini_x0 + c * cell_w, mini_y0 + r * cell_h),
                cell_w, cell_h,
                facecolor=cmap(norm(fake[r, c])),
                edgecolor="none", zorder=3))

    # Cohort × rhythm pool box
    pool_x = 0.300; pool_y = 0.13; pool_w = 0.305; pool_h = 0.30
    add_box(ax, pool_x, pool_y, pool_w, pool_h, "",
              fc=BOX_BG, ec=GREY_DARK, lw=0.8, radius=0.012)
    ax.text(pool_x + pool_w / 2, pool_y + pool_h - 0.024,
             "Block by rhythm,   pool across cohorts",
             fontsize=7.5, fontweight="bold", ha="center")

    cohort_labels = ["GSE128188",  "GSE79768",  "GSE115574"]
    n_sr = [5, 6, 15]; n_af = [5, 7, 13]
    base_y = pool_y + 0.040
    row_h = 0.058
    col_lab = pool_x + 0.012
    col_sr  = pool_x + 0.108
    col_af  = pool_x + 0.208
    block_w = 0.085; block_h = row_h - 0.018
    for j, c in enumerate(cohort_labels):
        y = base_y + (2 - j) * row_h
        ax.text(col_lab, y + row_h / 2 - 0.009, c,
                 fontsize=6.8, fontweight="bold",
                 color=GREY_DARK, ha="left", va="center")
        ax.add_patch(mpatches.Rectangle(
            (col_sr, y + 0.005), block_w, block_h,
            facecolor="#FFFFFF", edgecolor=GREY_MID, linewidth=0.5,
            zorder=3))
        ax.text(col_sr + block_w / 2,
                 y + 0.005 + block_h / 2,
                 f"SR  n = {n_sr[j]}",
                 fontsize=6.5, ha="center", va="center", color="black")
        ax.add_patch(mpatches.Rectangle(
            (col_af, y + 0.005), block_w, block_h,
            facecolor="#444444", edgecolor=GREY_MID, linewidth=0.5,
            zorder=3))
        ax.text(col_af + block_w / 2,
                 y + 0.005 + block_h / 2,
                 f"AF  n = {n_af[j]}",
                 fontsize=6.5, ha="center", va="center", color="white")

    # Hedges' g / pooling box
    hk_x = 0.650; hk_y = 0.13; hk_w = 0.330; hk_h = 0.30
    add_box(ax, hk_x, hk_y, hk_w, hk_h, "",
              fc=BG_GREEN, ec=HL_GREEN, lw=1.0, radius=0.012)
    ax.text(hk_x + hk_w / 2, hk_y + hk_h - 0.024,
             "Per-program AF rewiring effect",
             fontsize=7.5, fontweight="bold", ha="center", color="#1F5C39")
    ax.text(hk_x + hk_w / 2, hk_y + hk_h - 0.090,
             r"$g_k^{(c)}=\mathrm{Hedges'}\,g\left("
             r"\mathrm{aCAI}_{k,c}^{\mathrm{rel,AF}},"
             r"\mathrm{aCAI}_{k,c}^{\mathrm{rel,SR}}\right)$",
             fontsize=8.5, ha="center", va="center")
    ax.text(hk_x + hk_w / 2, hk_y + hk_h - 0.150,
             "pool 3 cohorts  →  fixed-effect IVW 95% CI (family)",
             fontsize=6.8, ha="center", va="center", color="#1F5C39")
    ax.text(hk_x + hk_w / 2, hk_y + hk_h - 0.180,
             "+ leave-one-cohort-out sensitivity",
             fontsize=6.2, ha="center", va="center", color=GREY_DARK,
             style="italic")
    ax.text(hk_x + hk_w / 2, hk_y + hk_h - 0.205,
             "+ paired-only re-fit",
             fontsize=6.2, ha="center", va="center", color=GREY_DARK,
             style="italic")
    ax.text(hk_x + hk_w / 2, hk_y + 0.030,
             "feeds  Fig 2C (states)   &   Fig 2D (effects/consistency)",
             fontsize=6.8, ha="center", va="center",
             color="#1F5C39", fontweight="bold")

    # arrows: STAGE 4 box → M_aCAI ; M_aCAI → pool ; pool → HK
    add_arrow(ax, s4_x + s4_w / 2, 0.580,
                mat_x + mat_w / 2, mat_y + mat_h + 0.005,
                color=GREY_DARK, lw=1.0,
                label="stack across donors × programs",
                label_offset=(-0.18, 0.01))
    add_arrow(ax, mat_x + mat_w + 0.005,
                mat_y + mat_h / 2,
                pool_x - 0.002, pool_y + pool_h / 2,
                color=GREY_DARK, lw=1.0)
    add_arrow(ax, pool_x + pool_w + 0.005,
                pool_y + pool_h / 2,
                hk_x - 0.002, hk_y + hk_h / 2,
                color=GREY_DARK, lw=1.0)

    fig.text(0.025, 0.030,
              "Notation:   i = donor;   k = NCP pathway program;   "
              "c = cohort;   P_k = program gene set;   "
              "expr = per-gene expression vector.",
              fontsize=6.2, color="#444", ha="left", va="bottom")
    fig.text(0.025, 0.010,
              "All inference is on aCAI, never on raw chamber expression "
              "— this controls for donor batch and bulk composition while "
              "preserving sign and magnitude of chamber asymmetry.",
              fontsize=6.2, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
