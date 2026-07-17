#!/usr/bin/env python3
"""
Figure 6D — Integrated G-tier traceability / interpretation matrix.

Reads:
  results/figure6_g_tier_traceability.tsv
  results/fig6a_anchor_coloc_calibration.tsv (optional display fields)
  results/fig6a_v3_locus_zoom_summary.tsv (locus callout distances)

Outputs:
  figures/fig6d_traceability_matrix.pdf
  results/fig6d_traceability_matrix.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import COL_DOUBLE, FG_COLOR, FS_LEGEND, FS_TICK, MUTED_COLOR, panel_label, set_nature_style, strip_spines
from _figure_export import pdf_safe_text, save_figure_bundle

G_TIER = ROOT / "results" / "figure6_g_tier_traceability.tsv"
LOCUS = ROOT / "results" / "fig6a_v3_locus_zoom_summary.tsv"
OUT_TSV = ROOT / "results" / "fig6d_traceability_matrix.tsv"
OUT_PDF = ROOT / "figures/fig6d_traceability_matrix.pdf"

TIER_COLORS = {
    "G0": "#E0E0E0",
    "G1": "#B8CCE4",
    "G2": "#F4D4A8",
    "G3": "#C5E0B4",
    "G4": "#6E9B72",
}
CALLOUT_GENES = ["SMAD9", "IDH2", "PITX2"]


def main() -> None:
    set_nature_style()
    gt = pd.read_csv(G_TIER, sep="\t")
    anchors = gt[gt["entity_type"] == "anchor"].copy()
    families = gt[gt["entity_type"] == "family"].copy()

    locus = pd.read_csv(LOCUS, sep="\t") if LOCUS.is_file() else pd.DataFrame()
    locus_map = locus.set_index("anchor") if "anchor" in locus.columns else {}

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    gt.to_csv(OUT_TSV, sep="\t", index=False)

    fig = plt.figure(figsize=(COL_DOUBLE, 0.95 * COL_DOUBLE))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.1, 0.75, 0.35], hspace=0.45,
                          left=0.10, right=0.97, top=0.90, bottom=0.08)

    # --- panel a: anchor G-tier strip ---
    ax_a = fig.add_subplot(gs[0])
    order = anchors["entity_id"].tolist()
    tiers = anchors.set_index("entity_id").loc[order, "G_tier"].tolist()
    x = np.arange(len(order))
    colors = [TIER_COLORS.get(t, "#CCC") for t in tiers]
    ax_a.bar(x, [1] * len(order), color=colors, edgecolor="white", linewidth=0.5)
    for i, g in enumerate(order):
        ax_a.text(i, 0.5, g, ha="center", va="center", fontsize=FS_LEGEND - 1,
                  rotation=90, fontweight="bold")
        ax_a.text(i, 1.08, tiers[i], ha="center", va="bottom", fontsize=FS_LEGEND - 1.5,
                  color=MUTED_COLOR)
    ax_a.set_xlim(-0.5, len(order) - 0.5)
    ax_a.set_ylim(0, 1.25)
    ax_a.set_xticks([])
    ax_a.set_yticks([])
    ax_a.set_title("d  Anchor-level G-tier traceability (12 pre-specified genes)",
                   loc="left", fontsize=FS_TICK, fontweight="bold")
    panel_label(ax_a, "d", x=-0.06, y=1.02)
    strip_spines(ax_a)

    # --- panel b: family G-tier strip ---
    ax_b = fig.add_subplot(gs[1])
    ford = families.sort_values("family_id")
    fx = np.arange(len(ford))
    ft = ford["G_tier"].tolist()
    fcolors = [TIER_COLORS.get(t, "#CCC") for t in ft]
    ax_b.barh(fx, [1] * len(ford), color=fcolors, edgecolor="white", linewidth=0.4)
    for i, (_, row) in enumerate(ford.iterrows()):
        ax_b.text(0.02, i, pdf_safe_text(row["entity_id"]), va="center",
                  fontsize=FS_LEGEND - 1)
        ax_b.text(0.98, i, row["G_tier"], va="center", ha="right",
                  fontsize=FS_LEGEND - 1.5, color=MUTED_COLOR)
    ax_b.set_yticks([])
    ax_b.set_xlim(0, 1)
    ax_b.set_title("Family-level G-tier (NCP leading-edge x external AF gene sets)",
                   loc="left", fontsize=FS_TICK, fontweight="bold")
    strip_spines(ax_b)

    # --- panel c: locus callouts + sQTL null footnote ---
    ax_c = fig.add_subplot(gs[2])
    ax_c.axis("off")
    lines = [
        "sQTL extension: all tested framework anchors - verdict_sqtl = no_eligible_sQTL "
        "(null at current junction-level coverage; not omitted).",
        "Locus callouts (frozen lead-variant offsets):",
    ]
    for g in CALLOUT_GENES:
        if g in locus_map.index:
            dist_kb = round(float(locus_map.loc[g]["delta_bp"]) / 1000, 1)
            lines.append(f"  - {g}: GWAS-eQTL lead distance ~ {dist_kb} kb")
        else:
            lines.append(f"  - {g}: calibration / H3-dominant locus (see Fig. 6a, 6b)")
    y = 0.95
    for ln in lines:
        ax_c.text(0.0, y, ln, transform=ax_c.transAxes, fontsize=FS_LEGEND - 0.5,
                  va="top", color=FG_COLOR if not ln.startswith("sQTL") else MUTED_COLOR)
        y -= 0.22

    # legend
    handles = [mpatches.Patch(facecolor=TIER_COLORS[t], edgecolor="none", label=t)
               for t in ["G0", "G1", "G2", "G3", "G4"]]
    ax_a.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 1.35),
                ncol=5, frameon=False, fontsize=FS_LEGEND - 1)

    fig.suptitle("Genetic traceability tiers - calibration-limited, not genetics-excluding",
                 y=0.98, fontsize=8, fontweight="bold")
    pdf_path, png_path = save_figure_bundle(fig, OUT_PDF)
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")
    print(f"wrote {pdf_path.relative_to(ROOT)}")
    if png_path:
        print(f"wrote {png_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
