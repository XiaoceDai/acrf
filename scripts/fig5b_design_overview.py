#!/usr/bin/env python3
"""
Figure 5B — Slim design / gate overview (was the previous "5A overview").

Compressed to a half-height strip so the main bridge panel (5A
triangulation) carries the visual weight.  Two columns:
  LEFT   : pre-registered support gate
           (6 frozen support programs × 2 expected cells)
  RIGHT  : outcome on GSE238242 LAA snRNA pseudobulks; gate NOT met
           → multiome regulator inference deferred (see Fig S5).

This panel is not a result; it documents the locked decision rule.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from plot_style import COL_DOUBLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "fig5b_design_overview.pdf"


def _box(fig, x, y, w, h, face, edge, lw=0.6):
    fig.patches.append(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.005,rounding_size=0.012",
        transform=fig.transFigure,
        facecolor=face, edgecolor=edge, linewidth=lw, zorder=1,
    ))


def _arrow(fig, a, b):
    fig.patches.append(FancyArrowPatch(
        a, b, transform=fig.transFigure, arrowstyle="-|>", mutation_scale=10,
        linewidth=0.85, color="#666666", shrinkA=2, shrinkB=2, zorder=3,
    ))


def main() -> None:
    set_nature_style()
    fig_w = COL_DOUBLE
    fig_h = 0.95
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.78, bottom=0.10)

    _box(fig, 0.025, 0.10, 0.43, 0.62, "#FFF4E0", "#C58A2D")
    _box(fig, 0.515, 0.10, 0.46, 0.62, "#FCEAEA", "#9C2E2E")

    fig.text(0.045, 0.62, "Pre-registered support gate",
             fontsize=6.0, fontweight="bold", color="#222222",
             ha="left", va="top")
    fig.text(0.045, 0.45,
             "STRUCTURAL → fibroblast (AF↑)   ·   METABOLIC → atrialCM (AF↓)",
             fontsize=5.0, color="#222222",
             ha="left", va="top",
             family="monospace")
    fig.text(0.045, 0.27,
             "Else: do not proceed to ATAC / motif / regulator inference.",
             fontsize=4.9, color="#555555", style="italic",
             ha="left", va="top")

    fig.text(0.530, 0.62, "GSE238242 LAA snRNA pseudobulk (AF n=3, SR n=4)",
             fontsize=6.0, fontweight="bold", color="#222222",
             ha="left", va="top")
    fig.text(0.530, 0.45,
             "Outcome: gate NOT met "
             "(no expected cell reaches moderate / strong support).",
             fontsize=5.1, color="#9C2E2E", fontweight="bold",
             ha="left", va="top")
    fig.text(0.530, 0.27,
             "Multiome regulator-tier inference deferred — see Fig S5.",
             fontsize=4.9, color="#555555", style="italic",
             ha="left", va="top")

    _arrow(fig, (0.460, 0.42), (0.512, 0.42))

    fig.text(0.012, 0.93, "b", fontsize=9, fontweight="bold",
             ha="left", va="top")
    fig.suptitle("Pre-registered support gate (deferred multiome track)",
                 fontsize=6.6, y=0.92, fontweight="bold")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
