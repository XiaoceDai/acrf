#!/usr/bin/env python3
"""Figure 1A - minimal Nature Methods-style formula panel."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import COL_DOUBLE, FG_COLOR, LA_COLOR, LINE_COLOR, MUTED_COLOR, RA_COLOR, set_nature_style

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

OUT = ROOT / "figures" / "fig1a_algorithm_problem_schematic.pdf"


def card(ax: plt.Axes, x: float, y: float, w: float, h: float, title: str, subtitle: str) -> None:
    ax.add_patch(
        mpatches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.006,rounding_size=0.012",
            facecolor="#FDFDFD",
            edgecolor="#D8D8D8",
            linewidth=0.70,
        )
    )
    ax.text(x + 0.018, y + h - 0.048, title, fontsize=8.9, fontweight="bold", ha="left", va="top", color=FG_COLOR)
    ax.text(x + 0.018, y + h - 0.092, subtitle, fontsize=6.8, ha="left", va="top", color=MUTED_COLOR)


def draw_design(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    card(
        ax,
        x,
        y,
        w,
        h,
        "Input: healthy paired atria",
        "8 donor pairs; each donor contributes one LA and one RA sample",
    )
    yh = y + h - 0.145
    ax.text(x + 0.02, yh, "donor", fontsize=6.9, color=MUTED_COLOR, fontweight="bold")
    ax.text(x + 0.18, yh, "LA", fontsize=6.9, color=MUTED_COLOR, fontweight="bold")
    ax.text(x + 0.26, yh, "RA", fontsize=6.9, color=MUTED_COLOR, fontweight="bold")
    ax.plot([x + 0.02, x + w - 0.02], [yh - 0.010, yh - 0.010], color=LINE_COLOR, lw=0.55)

    bw = 0.060
    bh = 0.030
    y0 = yh - 0.040
    step = 0.050
    for i in range(8):
        yy = y0 - i * step
        ax.text(x + 0.02, yy, f"D{i + 1}", fontsize=7.1, color=FG_COLOR, ha="left", va="center")
        ax.add_patch(mpatches.Rectangle((x + 0.16, yy - bh / 2), bw, bh, facecolor=LA_COLOR, edgecolor="#FFFFFF", linewidth=0.30))
        ax.add_patch(mpatches.Rectangle((x + 0.24, yy - bh / 2), bw, bh, facecolor=RA_COLOR, edgecolor="#FFFFFF", linewidth=0.30))
        ax.text(x + 0.16 + bw / 2, yy, "LA", fontsize=6.0, color="white", ha="center", va="center")
        ax.text(x + 0.24 + bw / 2, yy, "RA", fontsize=6.0, color="white", ha="center", va="center")

    ax.text(
        x + 0.02,
        y + 0.050,
        r"Primitive unit: paired contrast $(LA_d - RA_d)$ within each donor.",
        fontsize=6.9,
        color=MUTED_COLOR,
        ha="left",
    )


def draw_naive_ribbon(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    """Slim middle ribbon for baseline comparator."""
    ax.add_patch(
        mpatches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.005,rounding_size=0.010",
            facecolor="#F8F8F8",
            edgecolor="#DCDCDC",
            linewidth=0.65,
        )
    )
    ax.text(
        x + w / 2,
        y + h * 0.84,
        "Naive baseline",
        fontsize=7.1,
        fontweight="bold",
        color=MUTED_COLOR,
        ha="center",
        va="center",
    )
    ax.text(
        x + w / 2,
        y + h * 0.58,
        r"$y \sim chamber$",
        fontsize=7.8,
        color=FG_COLOR,
        ha="center",
        va="center",
    )
    ax.text(
        x + w / 2,
        y + h * 0.34,
        "drops donor pairing",
        fontsize=6.3,
        color=MUTED_COLOR,
        ha="center",
        va="center",
    )
    ax.text(
        x + w / 2,
        y + h * 0.20,
        "weaker chamber prior",
        fontsize=6.3,
        color=MUTED_COLOR,
        ha="center",
        va="center",
    )


def draw_formula(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    card(
        ax,
        x,
        y,
        w,
        h,
        "Chamber-aware primitive definition",
        "Donor-level equations used to build the healthy chamber prior",
    )
    # Primary equation 1: donor-level paired contrast
    ax.text(
        x + 0.018,
        y + h - 0.18,
        r"$(1)\ \ aCAI_{d,p}=\frac{1}{|G_p|}\sum_{g\in G_p}E^{LA}_{g,d}"
        r"-\frac{1}{|G_p|}\sum_{g\in G_p}E^{RA}_{g,d}$",
        fontsize=9.8,
        color=FG_COLOR,
        ha="left",
        va="top",
    )
    # Primary equation 2: cohort-level program score
    ax.text(
        x + 0.018,
        y + h - 0.33,
        r"$(2)\ \ \bar{aCAI}_p=\frac{1}{D}\sum_{d=1}^{D}aCAI_{d,p}$",
        fontsize=10.2,
        color=FG_COLOR,
        ha="left",
        va="top",
    )
    # Secondary metric: sign concordance across donors (de-emphasized)
    ax.text(
        x + 0.018,
        y + h - 0.47,
        r"$(3)\ \ c_p=\frac{1}{D}\sum_{d=1}^{D}\mathbf{1}\!\left[\mathrm{sign}(aCAI_{d,p})=\mathrm{sign}(\bar{aCAI}_p)\right]$",
        fontsize=7.4,
        color=MUTED_COLOR,
        ha="left",
        va="top",
    )
    ax.text(
        x + 0.018,
        y + 0.12,
        r"$D=8$ donors;  $G_p$ genes in program $p$;  $E$ normalized expression.",
        fontsize=6.7,
        color=MUTED_COLOR,
        ha="left",
    )
    ax.text(
        x + 0.018,
        y + 0.075,
        r"Sign rule: $\bar{aCAI}_p>0$ (LA-biased), $\bar{aCAI}_p<0$ (RA-biased).",
        fontsize=6.7,
        color=MUTED_COLOR,
        ha="left",
    )


def main() -> None:
    fig, ax = plt.subplots(figsize=(COL_DOUBLE * 1.08, COL_DOUBLE * 0.48))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.0,
        1.018,
        "a  Healthy left-right chamber asymmetry defined on donor-paired samples",
        fontsize=9.8,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=FG_COLOR,
    )

    draw_design(ax, 0.02, 0.13, 0.39, 0.76)
    draw_naive_ribbon(ax, 0.435, 0.385, 0.12, 0.25)
    draw_formula(ax, 0.58, 0.13, 0.40, 0.76)

    ax.annotate("", xy=(0.44, 0.51), xytext=(0.41, 0.51), arrowprops=dict(arrowstyle="->", lw=0.9, color=LINE_COLOR))
    ax.annotate("", xy=(0.58, 0.51), xytext=(0.55, 0.51), arrowprops=dict(arrowstyle="->", lw=0.9, color=LINE_COLOR))

    ax.text(
        0.02,
        0.048,
        "LA = red, RA = blue; Figure 1A uses healthy paired donors only (no rhythm split).",
        fontsize=6.5,
        color=MUTED_COLOR,
        ha="left",
        va="bottom",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
