#!/usr/bin/env python3
"""
Fig 3c v3 - Cellular fidelity classification (16 family).

Single horizontal stacked bar with 4 categories:
    Cellular fidelity (carrier sig only)
    Cellular relay    (non-carrier sig)
    Cellular extension (carrier + neighbour, same direction)
    Below detection   (no compartment p < 0.05)

Plus methodological footnote and per-class family list.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

CLS_TSV = ROOT / "results" / "fig3_v3_carrier_fidelity_class.tsv"
OUT_PDF = ROOT / "figures" / "fig3_c_fidelity_class.pdf"

CLASS_ORDER = ["Fidelity", "Relay", "Extension", "Below_detection"]
CLASS_LABELS = {
    "Fidelity": "Cellular fidelity\n(baseline carrier sig)",
    "Relay": "Cellular relay\n(non-carrier sig)",
    "Extension": "Cellular extension\n(carrier + neighbour, same direction)",
    "Below_detection": "Below detection\n(no compartment p < 0.05)",
}
CLASS_COLORS = {
    "Fidelity": "#7FB069",          # green (no carrier-color conflict)
    "Relay": "#8E6BAA",             # purple (was orange; avoid Endo orange in 3B)
    "Extension": "#5F9EA0",         # teal/cadet (was blue; avoid Fibro blue in 3B)
    "Below_detection": "#A9A9A9",   # neutral grey (distinct from Mixed grey hatch)
}


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def main():
    set_style()
    cls = pd.read_csv(CLS_TSV, sep="\t")
    counts = cls["fidelity_class"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int)
    total = int(counts.sum())

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 95 * MM))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[3.6, 1.4], height_ratios=[1.5, 1.0],
        left=0.05, right=0.985, top=0.86, bottom=0.26, wspace=0.06, hspace=0.30,
    )
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])
    ax_fams = fig.add_subplot(gs[1, :])

    # ---- stacked horizontal bar ----
    x_left = 0.0
    bar_y = 0.5
    bar_h = 0.55
    for cls_name in CLASS_ORDER:
        n = int(counts.get(cls_name, 0))
        if n == 0:
            continue
        ax_bar.barh(
            bar_y, n, left=x_left, height=bar_h,
            color=CLASS_COLORS[cls_name], edgecolor="white", linewidth=1.5,
        )
        # count text inside bar
        ax_bar.text(
            x_left + n / 2.0, bar_y, str(n),
            ha="center", va="center", fontsize=10, color="white", fontweight="bold",
        )
        x_left += n
    ax_bar.set_xlim(0, total)
    ax_bar.set_ylim(0, 1)
    ax_bar.set_yticks([])
    ax_bar.set_xticks(range(0, total + 1, 2))
    ax_bar.set_xlabel(f"families  ({total} total)", fontsize=7.2)
    ax_bar.tick_params(axis="x", labelsize=6.5)
    for s in ("top", "right", "left"):
        ax_bar.spines[s].set_visible(False)
    ax_bar.spines["bottom"].set_linewidth(0.4)

    # ---- legend ----
    ax_leg.axis("off")
    ax_leg.set_xlim(0, 1); ax_leg.set_ylim(1, 0)
    y = 0.05
    for cls_name in CLASS_ORDER:
        ax_leg.add_patch(mpatches.Rectangle(
            (0.0, y), 0.12, 0.12, facecolor=CLASS_COLORS[cls_name], edgecolor="none"
        ))
        ax_leg.text(
            0.16, y + 0.06, CLASS_LABELS[cls_name],
            fontsize=5.8, ha="left", va="center", color="#1A1A1A",
        )
        y += 0.24

    # ---- per-class family list ----
    ax_fams.axis("off")
    ax_fams.set_xlim(0, 1); ax_fams.set_ylim(1, 0)
    ax_fams.text(0.0, 0.05, "Members per class", fontsize=7.0, fontweight="bold", color="#1A1A1A")
    y = 0.20
    for cls_name in CLASS_ORDER:
        sub = cls[cls["fidelity_class"] == cls_name].sort_values("family_id")
        if sub.empty:
            txt = "(none)"
        else:
            members = []
            for _, r in sub.iterrows():
                fid = int(r["family_id"])
                # Only annotate cell-types when they are real (non-empty,
                # non-NaN); Below_detection rows have no sig_compartments.
                raw = r.get("sig_compartments", "")
                if pd.notna(raw) and str(raw).strip():
                    members.append(f"F{fid} ({raw})")
                else:
                    members.append(f"F{fid}")
            txt = ",  ".join(members)
        ax_fams.text(
            0.0, y, f"\u25A0  {cls_name.replace('_', ' ')}:",
            fontsize=6.2, fontweight="bold", color=CLASS_COLORS[cls_name],
        )
        ax_fams.text(0.18, y, txt, fontsize=6.0, color="#1A1A1A")
        y += 0.13

    fig.suptitle(
        "Fig 3c  Cellular fidelity classification (LA snRNA, 16 families)",
        x=0.05, y=0.96, ha="left", fontsize=8.8, fontweight="bold",
    )
    n_resolved = counts.get('Fidelity', 0) + counts.get('Relay', 0) + counts.get('Extension', 0)
    footer_lines = [
        f"GSE255612 LA snRNA resolves cellular attribution for only {n_resolved}/{total} families at nominal p < 0.05.",
        "Bulk AF rewiring signals (Fig 2) for the remaining 14 families are either",
        "(i) driven by RA-side changes not observable in this LA-only dataset,",
        "(ii) diluted below per-compartment detection limit at compartment-level pseudobulk resolution, or",
        "(iii) confounded by sub-cell-type composition shifts not captured by per-cell averaging.",
    ]
    for k, line in enumerate(footer_lines):
        fig.text(0.05, 0.165 - 0.025 * k, line, fontsize=5.8, color="#444444")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
