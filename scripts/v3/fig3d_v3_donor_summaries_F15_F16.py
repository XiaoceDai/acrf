#!/usr/bin/env python3
"""Fig 3d — family-level donor summaries for F15 (Relay) and F16 (Extension)."""
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

OUT_PDF = ROOT / "figures" / "fig3_d_donor_summaries_F15_F16.pdf"


def main() -> None:
    set_nature_style()
    scores, delta = load_tables()

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 130 * MM))
    gs = fig.add_gridspec(
        3, 2,
        height_ratios=[0.18, 1.0, 1.0],
        left=0.07, right=0.98, top=0.90, bottom=0.10,
        hspace=0.45, wspace=0.28,
    )

    ax_h = fig.add_subplot(gs[0, :])
    ax_h.axis("off")
    ax_h.set_xlim(0, 1)
    ax_h.set_ylim(0, 1)
    ax_h.add_patch(mpatches.FancyBboxPatch(
        (0.02, 0.15), 0.44, 0.70,
        boxstyle="round,pad=0.01,rounding_size=0.02",
        facecolor=CLASS_COLORS["Relay"], edgecolor="none",
    ))
    ax_h.text(0.24, 0.55, "F15  Relay", ha="center", va="center",
              color="white", fontsize=9, fontweight="bold")
    ax_h.text(0.24, 0.28, "endothelial only (fibroblast = HCA carrier)",
              ha="center", va="center", color="white", fontsize=6.0)
    ax_h.add_patch(mpatches.FancyBboxPatch(
        (0.54, 0.15), 0.44, 0.70,
        boxstyle="round,pad=0.01,rounding_size=0.02",
        facecolor=CLASS_COLORS["Extension"], edgecolor="none",
    ))
    ax_h.text(0.76, 0.55, "F16  Extension", ha="center", va="center",
              color="white", fontsize=9, fontweight="bold")
    ax_h.text(0.76, 0.28, "fibroblast + atrialCM (both sig)",
              ha="center", va="center", color="white", fontsize=6.0)

    # F15 endothelial (sig) + fibroblast carrier reference
    ax_f15e = fig.add_subplot(gs[1, 0])
    ax_f15f = fig.add_subplot(gs[1, 1])
    fam_e = family_donor_frame(scores, 15, "endothelial")
    d_e, p_e = lookup_delta(delta, 15, "endothelial")
    sens_e = lookup_sensitivity(15, "endothelial")
    draw_donor_pair(
        ax_f15e, fam_e,
        ctrl_color=SR_RHYTHM_COLOR, af_color=AF_RHYTHM_COLOR,
        title=f"F15 × {COMPARTMENT_SHORT['endothelial']}",
        subtitle="Relay support (non-carrier)",
        dscore=d_e, p_perm=p_e, highlight=True,
        sensitivity_p=sens_e,
    )
    ax_f15e.set_ylabel("family mean singscore", fontsize=7.5)

    fam_f = family_donor_frame(scores, 15, "fibroblast")
    d_f, p_f = lookup_delta(delta, 15, "fibroblast")
    draw_donor_pair(
        ax_f15f, fam_f,
        ctrl_color=SR_RHYTHM_COLOR, af_color=AF_RHYTHM_COLOR,
        title=f"F15 × {COMPARTMENT_SHORT['fibroblast']}",
        subtitle="HCA carrier reference (non-sig)",
        dscore=d_f, p_perm=p_f, highlight=False,
        carrier_note="HCA LA carrier = fibroblast",
    )

    # F16 fibroblast + atrialCM
    ax_f16f = fig.add_subplot(gs[2, 0])
    ax_f16a = fig.add_subplot(gs[2, 1])
    fam_16f = family_donor_frame(scores, 16, "fibroblast")
    d16f, p16f = lookup_delta(delta, 16, "fibroblast")
    draw_donor_pair(
        ax_f16f, fam_16f,
        ctrl_color=SR_RHYTHM_COLOR, af_color=AF_RHYTHM_COLOR,
        title=f"F16 × {COMPARTMENT_SHORT['fibroblast']}",
        subtitle="Extension support (HCA carrier)",
        dscore=d16f, p_perm=p16f, highlight=True,
        carrier_note="HCA LA carrier = fibroblast",
    )
    ax_f16f.set_ylabel("family mean singscore", fontsize=7.5)

    fam_16a = family_donor_frame(scores, 16, "atrialCM")
    d16a, p16a = lookup_delta(delta, 16, "atrialCM")
    draw_donor_pair(
        ax_f16a, fam_16a,
        ctrl_color=SR_RHYTHM_COLOR, af_color=AF_RHYTHM_COLOR,
        title=f"F16 × {COMPARTMENT_SHORT['atrialCM']}",
        subtitle="Extension support (neighbour)",
        dscore=d16a, p_perm=p16a, highlight=True,
    )

    fig.suptitle(
        "Fig 3d  Family-level donor scores — F15 Relay and F16 Extension",
        x=0.07, y=0.97, ha="left", fontsize=10, fontweight="bold",
    )
    fig.text(
        0.07, 0.015,
        "F15 support is endothelial-only at family nominal p < 0.05; fibroblast is shown as the HCA carrier reference.  "
        "F16 extension = fibroblast + atrialCM.  Primary p_perm annotated; 1e6 sensitivity for F15×endothelial when available.  "
        "F16×adipocyte trend (p = 0.054) is Supporting Figure only.  No gene-level bars.",
        fontsize=5.8, color="#666666",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
