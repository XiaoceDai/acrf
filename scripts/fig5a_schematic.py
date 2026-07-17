#!/usr/bin/env python3
"""Figure 5A  Triangulation analytical schematic."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FIG_OUT = ROOT / "figures" / "fig5a_schematic.pdf"

LA_COLOR  = "#D6604D"
RA_COLOR  = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"
BG_GREY    = "#FAFAFA"
BG_GREEN   = "#F2F8F4"
HL_GREEN   = "#3C8C57"
BG_GOLD    = "#FFFCF4"
BRBG_BROWN = "#8C510A"
BRBG_TEAL  = "#01665E"
PRGN_GREEN = "#1B7837"
PRGN_PURPLE = "#762A83"


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
    fig = plt.figure(figsize=(200 * MM, 140 * MM))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    fig.text(0.025, 0.965,
              "Fig 5A  Per-family triangulation design",
              fontsize=9.0, fontweight="bold", ha="left", va="center")
    fig.text(0.025, 0.940,
              "Three independent evidence streams — bulk cohort rewiring, "
              "LA snRNA carrier, and LAA duration trend — are gathered "
              "per axis family and assigned to an evidence class.",
              fontsize=7.5, color=GREY_DARK, ha="left", va="center",
              style="italic")

    # ROW 1   STAGE 1 (families)  →  STAGE 2 (three streams in parallel)
    s1_x, s1_w = 0.025, 0.185
    s2_x, s2_w = 0.235, 0.425   # wide three-panel strip
    s3_x, s3_w = 0.680, 0.300
    row_top = 0.895

    add_stage_header(ax, s1_x, row_top, 1, "11 axis families")
    add_stage_header(ax, s2_x, row_top, 2, "Three evidence streams")
    add_stage_header(ax, s3_x, row_top, 3, "Per-family support")

    card_y, card_h = 0.500, 0.360
    # families card
    add_box(ax, s1_x - 0.005, card_y, s1_w + 0.010, card_h, "",
              fc="#FCFCFC", ec=GREY_MID, lw=0.6, radius=0.012)
    # streams card (gold)
    add_box(ax, s2_x - 0.005, card_y, s2_w + 0.010, card_h, "",
              fc=BG_GOLD, ec="black", lw=0.9, radius=0.012)
    # support card
    add_box(ax, s3_x - 0.005, card_y, s3_w + 0.010, card_h, "",
              fc=BG_GREEN, ec=HL_GREEN, lw=0.9, radius=0.012)

    ax.text(s1_x + s1_w / 2, 0.830,
             "11 axis families",
             fontsize=8.0, fontweight="bold", ha="center")
    ax.text(s1_x + s1_w / 2, 0.805,
             "shared across Fig 2 / 3 / 4",
             fontsize=6.5, color=GREY_DARK, ha="center", style="italic")
    swatch_pal = [
        "#4F86C6", "#6FA8DC", "#9EC4E0",
        "#82B5A1", "#A2C8B0", "#B7D2BC",
        "#D4B58E", "#C8C26C",
        "#E89C8B", "#D86A6A", "#B94B4B",
    ]
    sw, sh = 0.022, 0.028
    rows, cols = 3, 4
    base_x = s1_x + (s1_w - cols * sw - (cols - 1) * 0.005) / 2
    base_y = 0.670
    for k in range(11):
        rr, cc = divmod(k, cols)
        x = base_x + cc * (sw + 0.005)
        y = base_y + (rows - 1 - rr) * (sh + 0.005)
        ax.add_patch(mpatches.Rectangle(
            (x, y), sw, sh, facecolor=swatch_pal[k],
            edgecolor="white", linewidth=0.4, zorder=3))
        ax.text(x + sw / 2, y + sh / 2, f"F{k+1}",
                 ha="center", va="center", fontsize=6.0,
                 fontweight="bold",
                 color="white" if k in (0, 8, 9, 10) else "black",
                 zorder=4)
    # LA/RA direction badges below swatches
    ax.text(s1_x + s1_w / 2, 0.580,
             "LA-biased:  F1 – F8",
             fontsize=6.3, ha="center", va="center",
             color=LA_COLOR, fontweight="bold")
    ax.text(s1_x + s1_w / 2, 0.555,
             "RA-biased:  F9 – F11",
             fontsize=6.3, ha="center", va="center",
             color=RA_COLOR, fontweight="bold")
    ax.text(s1_x + s1_w / 2, 0.525,
             "30 NCP_v2 programs",
             fontsize=6.0, ha="center", va="center", color=GREY_MID,
             style="italic")

    # Three vertical sub-columns with their own icon
    sub_w = (s2_w - 0.020) / 3
    sub_xs = [s2_x + 0.010 + k * (sub_w + 0.0) for k in range(3)]

    # column titles
    stream_titles = [
        ("L1  bulk aCAI", "Fig 2D pooled g", BRBG_BROWN),
        ("L2  LA snRNA",  "Fig 3B family Δ",  BRBG_BROWN),
        ("L3  LAA duration", "Fig 4B ρ",     PRGN_GREEN),
    ]
    for k, ((title, sub, col), x0) in enumerate(zip(stream_titles, sub_xs)):
        ax.text(x0 + sub_w / 2, 0.830, title,
                 fontsize=7.4, fontweight="bold", ha="center", va="center",
                 color=col)
        ax.text(x0 + sub_w / 2, 0.808, sub,
                 fontsize=6.3, ha="center", va="center",
                 color=GREY_DARK, style="italic")

    # iconographic data descriptor per stream
    stream_desc = [
        ("3 paired AF cohorts",
         "GSE128188  GSE79768  GSE115574",
         "REML + Hartung-Knapp pooling",
         "pooled g  per family"),
        ("GSE255612 LA snRNA",
         "34 donors × 11 cell types",
         "max-sig cell-type carrier",
         "Δ singscore  per family"),
        ("GSE41177  LAA primary",
         "12 AF donors  ·  10–150 mo",
         "Spearman ρ  vs duration",
         "signed ρ  per family"),
    ]
    for k, (a, b, c, d), x0 in zip(range(3), stream_desc, sub_xs):
        y_start = 0.730
        for idx, txt in enumerate([a, b, c, d]):
            y = y_start - idx * 0.042
            fw = "bold" if idx == 0 else "normal"
            size = 6.5 if idx == 0 else 6.0
            color = GREY_DARK if idx < 3 else ("#1F5C39")
            ax.text(x0 + sub_w / 2, y, txt,
                     fontsize=size, fontweight=fw,
                     color=color, ha="center", va="center")
        # tiny family × cohort / family × ct mini heatmap icon
        mini_x = x0 + 0.01
        mini_y = 0.545
        mini_w = sub_w - 0.02; mini_h = 0.040
        # just a gradient bar
        import numpy as np
        nrows, ncols = 4, 6
        cw = mini_w / ncols; ch = mini_h / nrows
        rng = np.random.default_rng(11 + k)
        fake = rng.normal(scale=0.5, size=(nrows, ncols))
        if k == 0:
            cmap = mpl.colors.LinearSegmentedColormap.from_list(
                "h", [(0.0, BRBG_TEAL), (0.5, "#FFFFFF"),
                         (1.0, BRBG_BROWN)])
            fake[-1, :] += 0.6
        elif k == 1:
            cmap = mpl.colors.LinearSegmentedColormap.from_list(
                "h", [(0.0, BRBG_TEAL), (0.5, "#FFFFFF"),
                         (1.0, BRBG_BROWN)])
        else:
            cmap = mpl.colors.LinearSegmentedColormap.from_list(
                "h", [(0.0, PRGN_PURPLE), (0.5, "#FFFFFF"),
                         (1.0, PRGN_GREEN)])
            fake[0:2, :] -= 0.4
        norm = mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
        for r in range(nrows):
            for c in range(ncols):
                ax.add_patch(mpatches.Rectangle(
                    (mini_x + c * cw, mini_y + r * ch),
                    cw, ch,
                    facecolor=cmap(norm(fake[r, c])),
                    edgecolor="none", zorder=3))

    # arrow from families card into streams card
    add_arrow(ax, s1_x + s1_w + 0.005, 0.700,
                s2_x - 0.005, 0.700,
                color=GREY_DARK, lw=1.0)

    # arrow from streams card into support card (per-family collapse)
    add_arrow(ax, s2_x + s2_w + 0.005, 0.700,
                s3_x - 0.005, 0.700,
                color=GREY_DARK, lw=1.0)

    ax.text(s3_x + s3_w / 2, 0.830,
             "Per-family support rules",
             fontsize=8.0, fontweight="bold", ha="center",
             color="#1F5C39")
    ax.text(s3_x + s3_w / 2, 0.805,
             "applied independently per layer",
             fontsize=6.2, ha="center", va="center",
             color=GREY_DARK, style="italic")
    rules = [
        ("L1 bulk",  "95 % CI excludes 0  (Hartung-Knapp pooled)"),
        ("L2 snRNA", "p_perm < 0.05   at the max-sig cell type"),
        ("L3 duration", "p_perm < 0.10   AND  |ρ| ≥ 0.50  (n = 12)"),
    ]
    for k, (lbl, body) in enumerate(rules):
        y = 0.742 - k * 0.050
        ax.text(s3_x + 0.012, y, lbl,
                 fontsize=6.8, fontweight="bold",
                 color="#1F5C39", ha="left", va="center")
        ax.text(s3_x + 0.085, y, body,
                 fontsize=6.2, color=GREY_DARK, ha="left", va="center",
                 style="italic")

    cls_x = s3_x + 0.010
    cls_y = 0.550
    cls_w = s3_w - 0.020
    cls_h = 0.050
    add_box(ax, cls_x, cls_y, cls_w, cls_h, "",
              fc=BG_GREY, ec="#1F5C39", lw=0.6, radius=0.010)
    ax.text(cls_x + cls_w / 2, cls_y + cls_h - 0.014,
             "Evidence class  =  " +
             " / ".join(["CONVERGENT", "BULK-LED", "SNRNA-LED",
                          "DURATION-LED", "SILENT"]),
             fontsize=5.8, ha="center", va="center",
             color="#1F5C39", fontweight="bold")
    ax.text(cls_x + cls_w / 2, cls_y + 0.014,
             "CONVERGENT = ≥ 2 supported layers with consistent sign",
             fontsize=5.7, ha="center", va="center",
             color=GREY_DARK, style="italic")

    # ROW 2   STAGE 4 + STAGE 5
    add_stage_header(ax, 0.025, 0.455, 4,
                       "Family × class matrix  ·  STAGE 5  downstream panels")

    # Left mini 11 × 3 heatmap (family × layer)
    mat_x = 0.025; mat_y = 0.120; mat_w = 0.220; mat_h = 0.310
    add_box(ax, mat_x, mat_y, mat_w, mat_h, "",
              fc="#FFFFFF", ec=GREY_DARK, lw=0.8, radius=0.012)
    ax.text(mat_x + mat_w / 2, mat_y + mat_h - 0.022,
             "Fig 5B  triangulation",
             fontsize=7.3, fontweight="bold", ha="center")
    ax.text(mat_x + mat_w / 2, mat_y + mat_h - 0.046,
             "11 families × 3 layers",
             fontsize=6.0, color=GREY_DARK, ha="center", style="italic")
    # mini 3-col heatmap
    import numpy as np
    mini_x0 = mat_x + 0.030; mini_y0 = mat_y + 0.020
    mini_w  = mat_w - 0.060; mini_h  = mat_h - 0.090
    nrows, ncols = 11, 3
    cw = mini_w / ncols; ch = mini_h / nrows
    rng = np.random.default_rng(3)
    fake = rng.normal(scale=0.3, size=(nrows, ncols))
    fake[8:, 0] = np.array([0.75, 0.65, 0.85])  # RA bulk_g large
    fake[9, 1]  = 0.20                             # F10 snRNA fib pos
    fake[10, 1] = -0.65                            # F11 snRNA lymph endo neg
    fake[1, 2] = -0.55                              # F2 duration neg
    fake[6, 2] = -0.30                              # F7 duration neg
    cmap_generic = mpl.colors.LinearSegmentedColormap.from_list(
        "h", [(0.0, BRBG_TEAL), (0.5, "#FFFFFF"), (1.0, BRBG_BROWN)])
    norm_generic = mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    for r in range(nrows):
        for c in range(ncols):
            ax.add_patch(mpatches.Rectangle(
                (mini_x0 + c * cw, mini_y0 + (nrows - 1 - r) * ch),
                cw, ch,
                facecolor=cmap_generic(norm_generic(fake[r, c])),
                edgecolor="white", linewidth=0.3, zorder=3))

    add_arrow(ax, s3_x + s3_w / 2, 0.505,
                mat_x + mat_w / 2, mat_y + mat_h + 0.005,
                color=GREY_DARK, lw=1.0)

    # Middle: 5C class table mini
    tab_x = 0.285; tab_y = 0.120; tab_w = 0.290; tab_h = 0.310
    add_box(ax, tab_x, tab_y, tab_w, tab_h, "",
              fc=BG_GREY, ec=GREY_DARK, lw=0.8, radius=0.012)
    ax.text(tab_x + tab_w / 2, tab_y + tab_h - 0.022,
             "Fig 5C  evidence class per family",
             fontsize=7.3, fontweight="bold", ha="center")
    # Class pills
    pill_classes = [
        ("CONVERGENT",   "#2E7B5E", "2 families"),
        ("BULK-LED",     "#8C510A", "1 family"),
        ("SNRNA-LED",    "#5E3C99", "0 families"),
        ("DURATION-LED", "#1F78B4", "1 family"),
        ("SILENT",       "#BBBBBB", "7 families"),
    ]
    py = tab_y + tab_h - 0.055
    for name, col, cnt in pill_classes:
        pill = mpatches.FancyBboxPatch(
            (tab_x + 0.015, py - 0.015),
            0.110, 0.030,
            boxstyle="round,pad=0.003,rounding_size=0.010",
            facecolor=col, edgecolor=col,
            linewidth=0.6, zorder=3)
        ax.add_patch(pill)
        ax.text(tab_x + 0.070, py, name,
                 ha="center", va="center",
                 fontsize=6.3, fontweight="bold", color="white",
                 zorder=4)
        ax.text(tab_x + 0.140, py, cnt,
                 ha="left", va="center",
                 fontsize=6.5, color=GREY_DARK, fontweight="bold")
        py -= 0.040

    add_arrow(ax, mat_x + mat_w + 0.005, tab_y + tab_h / 2,
                tab_x - 0.005, tab_y + tab_h / 2,
                color=GREY_DARK, lw=1.0)

    # Right: Fig 5D + 5E
    pan_x = 0.615; pan_y = 0.120; pan_w = 0.360; pan_h = 0.310
    add_box(ax, pan_x, pan_y, pan_w, pan_h, "",
              fc=BG_GREEN, ec=HL_GREEN, lw=1.0, radius=0.012)
    ax.text(pan_x + pan_w / 2, pan_y + pan_h - 0.022,
             "Downstream panels",
             fontsize=7.3, fontweight="bold", ha="center", color="#1F5C39")
    bullets = [
        ("Fig 5D",
         "GeoMx AA cardiac ROIs (n = 24):  multi-program × cell-type ρ matrix"),
        ("Fig 5E",
         "dual-axis mechanism cartoon  (structural  /  metabolic)"),
    ]
    for k, (tag, desc) in enumerate(bullets):
        y = pan_y + pan_h - 0.080 - k * 0.060
        ax.text(pan_x + 0.020, y, tag,
                 fontsize=7.0, fontweight="bold",
                 color="#1F5C39", ha="left", va="center")
        ax.text(pan_x + 0.075, y, desc,
                 fontsize=6.3, color=GREY_DARK, ha="left", va="center")

    ax.text(pan_x + pan_w / 2, pan_y + 0.060,
             "Fig 5B / 5C classify families  →  Fig 5D / 5E explain the two axes",
             fontsize=6.2, ha="center", va="center",
             color="#1F5C39", fontweight="bold", style="italic")

    add_arrow(ax, tab_x + tab_w + 0.005, tab_y + tab_h / 2,
                pan_x - 0.005, tab_y + tab_h / 2,
                color=GREY_DARK, lw=1.0)

    fig.text(0.025, 0.035,
              "Notation:  L1 / L2 / L3 = evidence layers;  "
              "g = pooled Hedges' g;  Δ = AF − Ctrl singscore change;  "
              "ρ = Spearman with AF duration in months.",
              fontsize=6.2, color="#444", ha="left", va="bottom")
    fig.text(0.025, 0.015,
              "CONVERGENT families carry ≥ 2 supported layers with "
              "consistent sign;  single-layer support gives a *-LED "
              "class;  SILENT families pass no layer and are kept for "
              "transparency.",
              fontsize=6.2, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
