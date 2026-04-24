#!/usr/bin/env python3
"""
Figure 6C (v3) — Convergent-balance synthesis schematic (three-layer vertical).

Layer 1: three equal-width evidence boxes (no decorative icons):
    cross-cohort axis · family-level programs · anchor evidence (10 + PITX2).
Layer 2: two parallel inferential legs (charcoal emphasis), each ending with
    a negation marker (✗) for the tested mediation mode.
Layer 3: unified verdict (thick SUPPORT_COLOR border, light fill).

Pure cartoon: no data coordinates, no causal arrowheads.  Centre line is
chosen from §10.11 via results/figure6_mr_anchors_v3_summary.tsv.

Outputs:
    figures/fig6c_v3_double_triangle.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch

from plot_style import (
    COL_DOUBLE,
    FG_COLOR,
    FAINT_COLOR,
    FS_AXIS_LABEL,
    FS_FOOTNOTE,
    FS_PANEL_LETTER,
    FS_PANEL_TITLE,
    MUTED_COLOR,
    SUPPORT_COLOR,
    set_nature_style,
)

ROOT = Path(__file__).resolve().parent.parent
SUM = ROOT / "results" / "figure6_mr_anchors_v3_summary.tsv"
OUT = ROOT / "figures" / "fig6c_v3_double_triangle.pdf"

CHARCOAL = "#2F2F2F"
LEG_FRAME_BG = "#FAFAFA"
LEG_FRAME_EDGE = CHARCOAL
MERGE_PT = (0.50, 0.725)

CENTRE_TEXTS = {
    "A": "Distributed program-level remodeling",
    "B": "Distributed program-level remodeling",
    "C": "Distributed program-level remodeling",
    "D": "Distributed (with weak hybrid signal at {gene})",
    "E": "Distributed + isolated single-locus driver: {gene}",
    "F": "Distributed + family-level driver axis: {family}",
    "G": "Mixed model: revise framing toward joint mechanism",
}


def _evidence_box(
    ax,
    cx: float,
    cy: float,
    *,
    title: str,
    line1: str,
    line2: str,
    width: float = 0.28,
    height: float = 0.115,
) -> None:
    x0 = cx - width / 2
    y0 = cy - height / 2
    ax.add_patch(
        FancyBboxPatch(
            (x0, y0),
            width,
            height,
            boxstyle="round,pad=0.0,rounding_size=0.01",
            facecolor=LEG_FRAME_BG,
            edgecolor=LEG_FRAME_EDGE,
            linewidth=0.75,
            zorder=2,
        )
    )
    ax.text(
        cx,
        cy + height * 0.22,
        title,
        ha="center",
        va="center",
        fontsize=FS_AXIS_LABEL,
        fontweight="bold",
        color=FG_COLOR,
        zorder=3,
    )
    ax.text(
        cx,
        cy - height * 0.05,
        line1,
        ha="center",
        va="center",
        fontsize=FS_FOOTNOTE,
        color=MUTED_COLOR,
        zorder=3,
    )
    ax.text(
        cx,
        cy - height * 0.32,
        line2,
        ha="center",
        va="center",
        fontsize=FS_FOOTNOTE,
        color=FAINT_COLOR,
        zorder=3,
        style="italic",
    )


def _negation_cross(ax, x: float, y: float, span: float = 0.022, lw: float = 1.15) -> None:
    """Draw a compact × (no special font glyphs)."""
    ax.plot(
        [x - span, x + span],
        [y - span, y + span],
        color=CHARCOAL,
        lw=lw,
        solid_capstyle="round",
        zorder=4,
    )
    ax.plot(
        [x - span, x + span],
        [y + span, y - span],
        color=CHARCOAL,
        lw=lw,
        solid_capstyle="round",
        zorder=4,
    )


def _inferential_leg(
    ax,
    cx: float,
    cy: float,
    *,
    title: str,
    body: str,
    fig_ref: str,
    width: float = 0.34,
    height: float = 0.168,
) -> tuple[float, float, float, float]:
    """Draw leg frame; return (x0, y0, width, height) in axes coords."""
    x0 = cx - width / 2
    y0 = cy - height / 2
    ax.add_patch(
        FancyBboxPatch(
            (x0, y0),
            width,
            height,
            boxstyle="round,pad=0.0,rounding_size=0.012",
            facecolor=LEG_FRAME_BG,
            edgecolor=LEG_FRAME_EDGE,
            linewidth=0.9,
            zorder=2,
        )
    )
    ax.text(
        cx,
        cy + height * 0.30,
        title,
        ha="center",
        va="center",
        fontsize=FS_AXIS_LABEL,
        fontweight="bold",
        color=CHARCOAL,
        zorder=3,
    )
    ax.text(
        cx,
        cy - height * 0.02,
        body,
        ha="center",
        va="center",
        fontsize=FS_FOOTNOTE,
        color=MUTED_COLOR,
        zorder=3,
        linespacing=1.35,
    )
    ax.text(
        cx,
        cy - height * 0.38,
        fig_ref,
        ha="center",
        va="center",
        fontsize=FS_FOOTNOTE,
        color=FAINT_COLOR,
        zorder=3,
        family="monospace",
    )
    # Negation mark (drawn ×; avoids missing-glyph fonts for ✗).
    _negation_cross(ax, cx + width * 0.40, cy + height * 0.34, span=0.019, lw=1.2)
    return x0, y0, width, height


def main() -> None:
    set_nature_style()

    summ = pd.read_csv(SUM, sep="\t").iloc[0]
    scenario_full = str(summ["scenario"])
    scenario_letter = scenario_full.strip()[0]
    centre_text_template = CENTRE_TEXTS.get(scenario_letter, CENTRE_TEXTS["A"])
    centre_text = centre_text_template

    fig_w = COL_DOUBLE
    fig_h = 3.75
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # ---- Layer 1: three evidence boxes (equal width, aligned) ----------
    y_e = 0.865
    gap = 0.024
    total_w = 0.90
    w_box = (total_w - 2 * gap) / 3
    x0_row = 0.5 - total_w / 2
    centers = [x0_row + w_box / 2 + i * (w_box + gap) for i in range(3)]

    _evidence_box(
        ax,
        centers[0],
        y_e,
        title="Cross-cohort axis",
        line1="Random-effects meta · locus architecture",
        line2="Fig 2E · Fig 6A",
        width=w_box,
        height=0.118,
    )
    _evidence_box(
        ax,
        centers[1],
        y_e,
        title="Family-level programs",
        line1="Atlas context · RA-biased families",
        line2="Fig 5",
        width=w_box,
        height=0.118,
    )
    _evidence_box(
        ax,
        centers[2],
        y_e,
        title="Anchor evidence",
        line1="10 pre-registered anchors + PITX2 calibrator",
        line2="Fig 6A · Fig 6B",
        width=w_box,
        height=0.118,
    )

    y_bottom_evidence = y_e - 0.118 / 2

    # ---- Connectors: evidence → merge → legs (solid charcoal) ----------
    leg_cx = (0.29, 0.71)
    leg_cy = 0.548
    y_top_leg = leg_cy + 0.168 / 2

    for xc in centers:
        ax.plot(
            [xc, MERGE_PT[0]],
            [y_bottom_evidence, MERGE_PT[1]],
            color=CHARCOAL,
            lw=0.85,
            solid_capstyle="round",
            zorder=1,
        )
    ax.plot(
        [MERGE_PT[0], leg_cx[0]],
        [MERGE_PT[1], y_top_leg],
        color=CHARCOAL,
        lw=0.85,
        solid_capstyle="round",
        zorder=1,
    )
    ax.plot(
        [MERGE_PT[0], leg_cx[1]],
        [MERGE_PT[1], y_top_leg],
        color=CHARCOAL,
        lw=0.85,
        solid_capstyle="round",
        zorder=1,
    )

    # ---- Layer 2: inferential legs --------------------------------------
    _inferential_leg(
        ax,
        leg_cx[0],
        leg_cy,
        title="Leg I · Single-gene mediation",
        body="Bulk atrial cis-eQTL × AF GWAS\ncis–trans colocalisation",
        fig_ref="Fig 6B (+ sQTL layer)",
    )
    _inferential_leg(
        ax,
        leg_cx[1],
        leg_cy,
        title="Leg II · Single-cell attribution",
        body="Carrier-cell framing for\nprogram drivers",
        fig_ref="Fig 3D–E",
    )

    # ---- Layer 3: verdict -----------------------------------------------
    centre_x, centre_y = 0.50, 0.168
    box_w, box_h = 0.64, 0.148
    box_x = centre_x - box_w / 2
    box_y = centre_y - box_h / 2
    ax.add_patch(
        FancyBboxPatch(
            (box_x, box_y),
            box_w,
            box_h,
            boxstyle="round,pad=0.0,rounding_size=0.018",
            facecolor="#F2F6F2",
            edgecolor=SUPPORT_COLOR,
            linewidth=2.1,
            zorder=3,
        )
    )
    y_leg_bottom = leg_cy - 0.168 / 2
    y_verdict_top = box_y + box_h
    for x_leg in leg_cx:
        ax.plot(
            [x_leg, centre_x],
            [y_leg_bottom, y_verdict_top],
            color=CHARCOAL,
            lw=0.95,
            solid_capstyle="round",
            zorder=1,
        )

    ax.text(
        centre_x,
        box_y + box_h - 0.028,
        f"Unified verdict (scenario {scenario_full.strip()})",
        ha="center",
        va="center",
        fontsize=FS_FOOTNOTE,
        color=MUTED_COLOR,
        fontweight="bold",
        style="italic",
        zorder=4,
    )
    ax.text(
        centre_x,
        box_y + box_h - 0.072,
        centre_text,
        ha="center",
        va="center",
        fontsize=FS_PANEL_TITLE,
        fontweight="bold",
        color=FG_COLOR,
        zorder=4,
    )
    ax.text(
        centre_x,
        box_y + 0.034,
        "Parallel legs do not support single-gene or single-carrier mediation;\n"
        "evidence converges on a distributed program-level substrate.",
        ha="center",
        va="center",
        fontsize=FS_FOOTNOTE,
        color=MUTED_COLOR,
        style="italic",
        zorder=4,
        linespacing=1.45,
    )

    fig.suptitle(
        "Convergent balance: independent evidence → parallel inferential legs → "
        "one verdict",
        fontsize=FS_PANEL_TITLE,
        y=0.972,
        fontweight="bold",
        color=FG_COLOR,
    )
    fig.text(
        0.012,
        0.972,
        "c",
        fontsize=FS_PANEL_LETTER,
        fontweight="bold",
        ha="left",
        va="top",
        transform=fig.transFigure,
        color=FG_COLOR,
    )

    fig.text(
        0.50,
        0.055,
        "Evidence: (i) Fig 2E + Fig 6A   ·   (ii) Fig 5   ·   (iii) Fig 6A–B   ·   "
        "Legs: coloc / carriers",
        ha="center",
        va="center",
        fontsize=FS_FOOTNOTE,
        color=MUTED_COLOR,
    )
    fig.text(
        0.50,
        0.018,
        "Schematic only.  No effect sizes; connectors are non-causal (no arrowheads).  "
        "Verdict template from spec §10.11; scenario rule §10.8.",
        ha="center",
        va="bottom",
        fontsize=FS_FOOTNOTE,
        color=FAINT_COLOR,
        style="italic",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
