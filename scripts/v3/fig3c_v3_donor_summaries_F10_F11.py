#!/usr/bin/env python3
"""Fig 3c — family-level donor summaries for F10 (Fidelity) and F11 (Relay)."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "v3"))

from fig3_donor_summary_common import (  # noqa: E402
    CLASS_COLORS,
    COMPARTMENT_SHORT,
    draw_donor_pair,
    family_donor_frame,
    load_tables,
    lookup_delta,
    lookup_sensitivity,
)
from plot_style import AF_RHYTHM_COLOR, SR_RHYTHM_COLOR, set_nature_style  # noqa: E402

OUT_PDF = ROOT / "figures" / "fig3_c_donor_summaries_F10_F11.pdf"


def main() -> None:
    set_nature_style()
    scores, delta = load_tables()

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 115 * MM))
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[0.22, 1.0],
        left=0.07, right=0.98, top=0.88, bottom=0.12,
        hspace=0.55, wspace=0.28,
    )

    # header badges
    ax_h = fig.add_subplot(gs[0, :])
    ax_h.axis("off")
    ax_h.set_xlim(0, 1)
    ax_h.set_ylim(0, 1)
    ax_h.add_patch(mpatches.FancyBboxPatch(
        (0.02, 0.15), 0.44, 0.70,
        boxstyle="round,pad=0.01,rounding_size=0.02",
        facecolor=CLASS_COLORS["Fidelity"], edgecolor="none",
    ))
    ax_h.text(0.24, 0.55, "F10  Fidelity", ha="center", va="center",
              color="white", fontsize=9, fontweight="bold")
    ax_h.text(0.24, 0.28, "fibroblast (HCA carrier, sig)", ha="center", va="center",
              color="white", fontsize=6.2)
    ax_h.add_patch(mpatches.FancyBboxPatch(
        (0.54, 0.15), 0.44, 0.70,
        boxstyle="round,pad=0.01,rounding_size=0.02",
        facecolor=CLASS_COLORS["Relay"], edgecolor="none",
    ))
    ax_h.text(0.76, 0.55, "F11  Relay", ha="center", va="center",
              color="white", fontsize=9, fontweight="bold")
    ax_h.text(0.76, 0.28, "adipocyte sig; endothelial HCA carrier (flat)",
              ha="center", va="center", color="white", fontsize=6.2)

    ax_f10 = fig.add_subplot(gs[1, 0])
    ax_f11a = fig.add_subplot(gs[1, 1])

    # F10 fibroblast
    fam = family_donor_frame(scores, 10, "fibroblast")
    d, p = lookup_delta(delta, 10, "fibroblast")
    draw_donor_pair(
        ax_f10, fam,
        ctrl_color=SR_RHYTHM_COLOR, af_color=AF_RHYTHM_COLOR,
        title=f"F10 × {COMPARTMENT_SHORT['fibroblast']}",
        subtitle="Fidelity support (HCA carrier)",
        dscore=d, p_perm=p, highlight=True,
        carrier_note="HCA LA carrier = fibroblast",
    )
    ax_f10.set_ylabel("family mean singscore", fontsize=7.5)

    # F11: show adipocyte (sig) as primary; add small twin for endothelial reference via second axes in same cell? 
    # Plan: F11 Relay: adipocyte (sig) + HCA endothelial carrier reference (show carrier flat/non-sig)
    # Use nested gridspec within right column.
    # Rebuild right column as 1x2 nested.
    fig.delaxes(ax_f11a)
    inner = gs[1, 1].subgridspec(1, 2, wspace=0.35)
    ax_ad = fig.add_subplot(inner[0, 0])
    ax_en = fig.add_subplot(inner[0, 1])

    fam_ad = family_donor_frame(scores, 11, "adipocyte")
    d_ad, p_ad = lookup_delta(delta, 11, "adipocyte")
    sens_ad = lookup_sensitivity(11, "adipocyte")
    draw_donor_pair(
        ax_ad, fam_ad,
        ctrl_color=SR_RHYTHM_COLOR, af_color=AF_RHYTHM_COLOR,
        title=f"F11 × {COMPARTMENT_SHORT['adipocyte']}",
        subtitle="Relay support (non-carrier)",
        dscore=d_ad, p_perm=p_ad, highlight=True,
        sensitivity_p=sens_ad,
    )

    fam_en = family_donor_frame(scores, 11, "endothelial")
    d_en, p_en = lookup_delta(delta, 11, "endothelial")
    draw_donor_pair(
        ax_en, fam_en,
        ctrl_color=SR_RHYTHM_COLOR, af_color=AF_RHYTHM_COLOR,
        title=f"F11 × {COMPARTMENT_SHORT['endothelial']}",
        subtitle="HCA carrier reference (non-sig)",
        dscore=d_en, p_perm=p_en, highlight=False,
        carrier_note="HCA LA carrier = endothelial",
    )

    fig.suptitle(
        "Fig 3c  Family-level donor scores — F10 Fidelity and F11 Relay",
        x=0.07, y=0.965, ha="left", fontsize=10, fontweight="bold",
    )
    fig.text(
        0.07, 0.02,
        "Points: donor family-mean exact singscore (program-first mean of frozen-LE program scores).  "
        "Horizontal ticks: group means.  Annotated p_perm are primary 10,000-MC values; "
        "1e6 sensitivity shown for the F11×adipocyte boundary cell when available.  "
        "No gene-level bars.",
        fontsize=5.8, color="#666666",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
