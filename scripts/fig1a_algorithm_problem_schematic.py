#!/usr/bin/env python3
"""Figure 1A — donor-paired chamber-aware prior schematic (methods opening panel).

Layout matches the manuscript reference: Step 1 → Step 2 (core contrast) → Step 3
(prior aggregation), with Step 2b naive baseline off the main path, and a bottom
summary ribbon + LA/RA legend.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import COL_DOUBLE, FG_COLOR, MUTED_COLOR, set_nature_style

set_nature_style()
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.default": "regular",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

# Reference palette (mockup)
LA_HEX = "#D65A47"
RA_HEX = "#4E79A7"
CARD_BG = "#F5F5F5"
CARD_EDGE = "#C8C8C8"
LINE_MUTED = "#888888"

OUT = ROOT / "figures" / "fig1a_algorithm_problem_schematic.pdf"


def step_badge(ax: plt.Axes, cx: float, cy: float, n: str, *, accent: bool = False) -> None:
    fc = "#E8E8E8" if not accent else "#FAD4CE"
    ec = "#B0B0B0" if not accent else LA_HEX
    lw = 0.8 if not accent else 1.2
    circ = plt.Circle((cx, cy), 0.022, transform=ax.transAxes, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=5)
    ax.add_patch(circ)
    ax.text(cx, cy, n, transform=ax.transAxes, fontsize=8.5, fontweight="bold",
            ha="center", va="center", color=FG_COLOR if not accent else LA_HEX, zorder=6)


def rounded_card(ax: plt.Axes, x: float, y: float, w: float, h: float, *, ec: str | None = None, lw: float = 0.8) -> None:
    ax.add_patch(
        mpatches.FancyBboxPatch(
            (x, y), w, h,
            transform=ax.transAxes,
            boxstyle="round,pad=0.008,rounding_size=0.014",
            facecolor=CARD_BG,
            edgecolor=ec or CARD_EDGE,
            linewidth=lw,
            clip_on=False,
        )
    )


def draw_step1(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    rounded_card(ax, x, y, w, h)
    step_badge(ax, x + 0.028, y + h - 0.038, "1")
    ax.text(
        x + 0.055, y + h - 0.028,
        "Paired atrial samples (Healthy donors)",
        transform=ax.transAxes, fontsize=8.2, fontweight="bold", color=FG_COLOR, va="top", ha="left",
    )

    yh = y + h - 0.11
    ax.text(x + 0.04, yh, "donor", fontsize=6.8, color=MUTED_COLOR, fontweight="bold", transform=ax.transAxes)
    ax.text(x + 0.22, yh, "LA", fontsize=6.8, color=LA_HEX, fontweight="bold", transform=ax.transAxes)
    ax.text(x + 0.34, yh, "RA", fontsize=6.8, color=RA_HEX, fontweight="bold", transform=ax.transAxes)
    ax.plot([x + 0.03, x + w - 0.04], [yh - 0.012, yh - 0.012], color=LINE_MUTED, lw=0.55, transform=ax.transAxes)

    bw, bh = 0.048, 0.028
    y0 = yh - 0.028
    step_y = 0.048
    for i in range(8):
        yy = y0 - i * step_y
        ax.text(x + 0.035, yy, f"D{i + 1}", fontsize=7.0, color=FG_COLOR, ha="left", va="center", transform=ax.transAxes)
        ax.add_patch(
            mpatches.Rectangle(
                (x + 0.18, yy - bh / 2), bw, bh, transform=ax.transAxes,
                facecolor=LA_HEX, edgecolor="white", linewidth=0.35, zorder=3,
            )
        )
        ax.add_patch(
            mpatches.Rectangle(
                (x + 0.30, yy - bh / 2), bw, bh, transform=ax.transAxes,
                facecolor=RA_HEX, edgecolor="white", linewidth=0.35, zorder=3,
            )
        )
        ax.text(x + 0.18 + bw / 2, yy, "LA", fontsize=6.0, color="white", ha="center", va="center", transform=ax.transAxes, fontweight="bold")
        ax.text(x + 0.30 + bw / 2, yy, "RA", fontsize=6.0, color="white", ha="center", va="center", transform=ax.transAxes, fontweight="bold")
        # connector
        ax.plot([x + 0.18 + bw, x + 0.30], [yy, yy], color=FG_COLOR, lw=0.65, transform=ax.transAxes, zorder=2)

    yfoot = y + 0.10
    ax.text(x + 0.04, yfoot, "×8 donor pairs total", fontsize=7.0, color=FG_COLOR, fontweight="bold", transform=ax.transAxes)
    # simple "two people" icon
    icx, icy = x + 0.04, yfoot - 0.055
    for dx in (0, 0.022):
        ax.add_patch(plt.Circle((icx + dx, icy), 0.011, transform=ax.transAxes, facecolor=MUTED_COLOR, zorder=4))
    ax.plot([icx - 0.006, icx + 0.028], [icy - 0.018, icy - 0.018], color=MUTED_COLOR, lw=1.0, transform=ax.transAxes)
    ax.text(
        x + 0.095, yfoot - 0.035,
        "Each donor contributes one LA and one RA sample.",
        fontsize=6.5, color=MUTED_COLOR, transform=ax.transAxes, va="center", ha="left",
    )


def draw_step2(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    rounded_card(ax, x, y, w, h, ec=LA_HEX, lw=1.35)
    step_badge(ax, x + 0.032, y + h - 0.038, "2", accent=True)
    ax.text(
        x + 0.062, y + h - 0.028,
        "Within-donor contrast",
        transform=ax.transAxes, fontsize=8.2, fontweight="bold", color=FG_COLOR, va="top", ha="left",
    )
    cy = y + h * 0.48
    cx = x + w / 2
    ax.text(cx - 0.065, cy, r"$LA_d$", fontsize=15, color=LA_HEX, ha="right", va="center", transform=ax.transAxes, fontweight="bold")
    ax.text(cx - 0.038, cy, "−", fontsize=15, color=FG_COLOR, ha="center", va="center", transform=ax.transAxes)
    ax.text(cx + 0.038, cy, r"$RA_d$", fontsize=15, color=RA_HEX, ha="left", va="center", transform=ax.transAxes, fontweight="bold")
    ax.text(
        x + w / 2, y + h * 0.22,
        "Removes donor background and preserves\nchamber-specific signal.",
        fontsize=6.8, color=MUTED_COLOR, ha="center", va="center", transform=ax.transAxes,
    )


def draw_step2b(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    ax.add_patch(
        mpatches.FancyBboxPatch(
            (x, y), w, h,
            transform=ax.transAxes,
            boxstyle="round,pad=0.006,rounding_size=0.012",
            facecolor="#FAFAFA",
            edgecolor="#B8B8B8",
            linewidth=0.9,
            linestyle=(0, (3, 2)),
        )
    )
    ax.text(x + w / 2, y + h - 0.06, "Naive baseline", fontsize=7.0, fontweight="bold",
            color=MUTED_COLOR, ha="center", va="top", transform=ax.transAxes)
    ax.text(x + w / 2, y + h * 0.52, r"$y \sim \mathrm{chamber}$", fontsize=8.5, color=FG_COLOR,
            ha="center", va="center", transform=ax.transAxes)
    ax.text(x + 0.06, y + h * 0.28, "✗  ignores donor pairing", fontsize=6.5, color=MUTED_COLOR,
            ha="left", va="center", transform=ax.transAxes)
    ax.text(x + 0.06, y + h * 0.12, "✗  more sensitive to inter-donor variability", fontsize=6.5, color=MUTED_COLOR,
            ha="left", va="center", transform=ax.transAxes)


def draw_step3(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    rounded_card(ax, x, y, w, h)
    step_badge(ax, x + 0.022, y + h - 0.038, "3")
    ax.text(
        x + 0.050, y + h - 0.028,
        "Healthy chamber prior",
        transform=ax.transAxes, fontsize=8.2, fontweight="bold", color=FG_COLOR, va="top", ha="left",
    )
    ax.text(
        x + 0.022, y + h - 0.072,
        "Aggregate donor-paired contrasts across programs",
        fontsize=6.6, color=MUTED_COLOR, ha="left", va="top", transform=ax.transAxes,
    )

    y0 = y + h - 0.14
    ax.text(
        x + 0.022, y0,
        r"$aCAI_{d,p} = \frac{1}{|G_p|}\sum_{g \in G_p} E^{LA}_{g,d}"
        r" - \frac{1}{|G_p|}\sum_{g \in G_p} E^{RA}_{g,d}$",
        fontsize=8.0, color=FG_COLOR, ha="left", va="top", transform=ax.transAxes,
    )
    ax.text(
        x + 0.022, y0 - 0.10,
        r"$\overline{aCAI}_p = \frac{1}{D}\sum_{d=1}^{D} aCAI_{d,p}$",
        fontsize=8.5, color=FG_COLOR, ha="left", va="top", transform=ax.transAxes,
    )

    yleg = y + h - 0.34
    ax.annotate("", xy=(x + 0.045, yleg), xytext=(x + 0.015, yleg), arrowprops=dict(arrowstyle="->", color=LA_HEX, lw=1.2),
                transform=ax.transAxes)
    ax.text(x + 0.052, yleg, r"$\overline{aCAI}_p > 0$  →  LA-biased", fontsize=6.5, color=FG_COLOR, va="center", transform=ax.transAxes)
    ax.annotate("", xy=(x + 0.015, yleg - 0.045), xytext=(x + 0.045, yleg - 0.045), arrowprops=dict(arrowstyle="->", color=RA_HEX, lw=1.2),
                transform=ax.transAxes)
    ax.text(x + 0.052, yleg - 0.045, r"$\overline{aCAI}_p < 0$  →  RA-biased", fontsize=6.5, color=FG_COLOR, va="center", transform=ax.transAxes)
    ax.text(x + 0.015, yleg - 0.095, r"$\Rightarrow\!\Rightarrow$", fontsize=9, color=MUTED_COLOR, va="center", transform=ax.transAxes)
    ax.text(x + 0.052, yleg - 0.095, r"higher $c_p$  $\rightarrow$  stronger prior (donor sign concordance)", fontsize=6.3, color=MUTED_COLOR,
            va="center", transform=ax.transAxes)

    ax.text(
        x + 0.022, y + 0.06,
        r"$D=8$ donors;  $G_p$ = genes in program $p$;  $E$ = normalized expression.",
        fontsize=6.4, color=MUTED_COLOR, ha="left", va="bottom", transform=ax.transAxes,
    )


def draw_bottom_banner(ax: plt.Axes) -> None:
    bx, by, bw, bh = 0.02, 0.02, 0.96, 0.11
    rounded_card(ax, bx, by, bw, bh, ec="#D0D0D0", lw=0.7)
    # target / bullseye
    tcx, tcy = bx + 0.055, by + bh / 2
    rings = [(0.028, 1.0, LA_HEX), (0.017, 0.75, RA_HEX), (0.007, 0.5, LA_HEX)]
    for r, lw, ec in rings:
        ax.add_patch(plt.Circle((tcx, tcy), r, transform=ax.transAxes, fill=False, edgecolor=ec, linewidth=lw, zorder=4))
    ax.add_patch(plt.Circle((tcx, tcy), 0.0035, transform=ax.transAxes, facecolor=RA_HEX, zorder=5))

    ax.text(
        bx + 0.11, by + bh * 0.58,
        "Pairing controls inter-donor variation and yields a more reliable chamber-aware prior.",
        fontsize=7.2, color=FG_COLOR, ha="left", va="center", transform=ax.transAxes, fontweight="bold",
    )

    lx = bx + 0.11
    ly = by + bh * 0.22
    ax.add_patch(mpatches.Rectangle((lx, ly - 0.008), 0.018, 0.016, transform=ax.transAxes, facecolor=LA_HEX, edgecolor="none"))
    ax.text(lx + 0.026, ly, "LA", fontsize=6.8, color=FG_COLOR, va="center", transform=ax.transAxes)
    ax.add_patch(mpatches.Rectangle((lx + 0.055, ly - 0.008), 0.018, 0.016, transform=ax.transAxes, facecolor=RA_HEX, edgecolor="none"))
    ax.text(lx + 0.083, ly, "RA", fontsize=6.8, color=FG_COLOR, va="center", transform=ax.transAxes)


def main() -> None:
    fig_h = COL_DOUBLE * 0.56
    fig, ax = plt.subplots(figsize=(COL_DOUBLE * 1.05, fig_h))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    ax.text(
        0.0, 0.995,
        "a  Donor-paired contrasts define a chamber-aware healthy asymmetry prior.",
        fontsize=9.4, fontweight="bold", ha="left", va="top", color=FG_COLOR, transform=ax.transAxes,
    )

    # Step 1 (left)
    draw_step1(ax, 0.02, 0.16, 0.24, 0.78)

    # Step 2 + 2b (center column)
    draw_step2(ax, 0.29, 0.38, 0.22, 0.56)
    draw_step2b(ax, 0.29, 0.16, 0.22, 0.19)

    # Step 3 (right)
    draw_step3(ax, 0.545, 0.16, 0.435, 0.78)

    # Arrows: 1 → 2 (solid)
    ax.annotate(
        "",
        xy=(0.29, 0.66), xycoords=ax.transAxes,
        xytext=(0.26, 0.66), textcoords=ax.transAxes,
        arrowprops=dict(arrowstyle="->", color=FG_COLOR, lw=1.35, shrinkA=2, shrinkB=2),
    )
    # 2 → 3 (solid)
    ax.annotate(
        "",
        xy=(0.545, 0.66), xycoords=ax.transAxes,
        xytext=(0.51, 0.66), textcoords=ax.transAxes,
        arrowprops=dict(arrowstyle="->", color=FG_COLOR, lw=1.35, shrinkA=2, shrinkB=2),
    )
    # 1 → 2b (dashed, lower)
    ax.annotate(
        "",
        xy=(0.29, 0.28), xycoords=ax.transAxes,
        xytext=(0.24, 0.32), textcoords=ax.transAxes,
        arrowprops=dict(arrowstyle="->", color=MUTED_COLOR, lw=1.0, linestyle=(0, (4, 3)), shrinkA=2, shrinkB=2),
    )

    draw_bottom_banner(ax)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
