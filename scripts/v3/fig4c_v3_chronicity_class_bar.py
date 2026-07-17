#!/usr/bin/env python3
"""Fig 4c v3 - 3-class chronicity classification stacked bar (effect-size driven)."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

CLS_TSV = ROOT / "results" / "fig4_v3_chronicity_class.tsv"
ATLAS_TSV = ROOT / "results" / "fig_integrated_af_rewiring_atlas.tsv"
OUT_PDF = ROOT / "figures" / "fig4_c_chronicity_class.pdf"

CLASS_ORDER = ["Strong_trend", "Weak_trend", "Stable", "Inconclusive"]
CLASS_LABELS = {
    "Strong_trend": "Strong-trend\n(|\u03c1| \u2265 0.4)",
    "Weak_trend": "Weak-trend\n(0.2 \u2264 |\u03c1| < 0.4)",
    "Stable": "Stable\n(|\u03c1| < 0.2)",
    "Inconclusive": "Inconclusive\n(NA)",
}
CLASS_COLORS = {
    "Strong_trend": "#C44569",   # deep rose - distinct from Fig 3c's Relay purple
    "Weak_trend": "#F8B500",     # amber
    "Stable": "#A9A9A9",         # neutral grey
    "Inconclusive": "#E0E0E0",   # light grey
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
    atlas = pd.read_csv(ATLAS_TSV, sep="\t").set_index("family_id")

    counts = cls["chronicity_class"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int)
    total = int(counts.sum())

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(200 * MM, 130 * MM))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[3.6, 1.4], height_ratios=[1.0, 1.4],
        left=0.05, right=0.985, top=0.88, bottom=0.18, wspace=0.06, hspace=0.30,
    )
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])
    ax_fams = fig.add_subplot(gs[1, :])

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

    ax_leg.axis("off")
    ax_leg.set_xlim(0, 1); ax_leg.set_ylim(1, 0)
    y = 0.05
    for cls_name in CLASS_ORDER:
        if counts.get(cls_name, 0) == 0:
            continue
        ax_leg.add_patch(mpatches.Rectangle(
            (0.0, y), 0.12, 0.10, facecolor=CLASS_COLORS[cls_name], edgecolor="none",
        ))
        ax_leg.text(
            0.16, y + 0.05, CLASS_LABELS[cls_name],
            fontsize=5.8, ha="left", va="center", color="#1A1A1A",
        )
        y += 0.20

    ax_fams.axis("off")
    ax_fams.set_xlim(0, 1); ax_fams.set_ylim(1, 0)
    ax_fams.text(0.0, 0.02, "Members per class  (with primary \u03c1 and sensitivity-concordance n / 4)",
                  fontsize=7.0, fontweight="bold", color="#1A1A1A")

    short_overrides = {
        2: "Mito Ca2+ trans", 9: "HS-GAG metab", 15: "NTRK2/TrkB",
        3: "TG metab", 7: "TCA cycle", 8: "Complex I",
        13: "Amino-acid metab", 10: "TGF-beta", 16: "Keratan",
        12: "rRNA proc", 6: "Mito prot import", 4: "BCAA cat",
        14: "EuTrans init", 1: "Ribosome QC", 5: "tRNA aa",
        11: "BMP signalling",
    }

    y = 0.10
    line_h = 0.045
    items_per_line = 3
    for cls_name in CLASS_ORDER:
        sub = cls[cls["chronicity_class"] == cls_name].sort_values(
            "rho_primary", key=lambda s: s.abs(), ascending=False,
        )
        if sub.empty:
            continue
        n_total = len(sub)
        # Header line: class label + count + threshold
        threshold_label = {
            "Strong_trend": "(|\u03c1| \u2265 0.4)",
            "Weak_trend": "(0.2 \u2264 |\u03c1| < 0.4)",
            "Stable": "(|\u03c1| < 0.2)",
            "Inconclusive": "(NA)",
        }[cls_name]
        ax_fams.text(0.0, y, f"\u25A0  {cls_name.replace('_', '-')}  {threshold_label}  n = {n_total}",
                      fontsize=6.6, fontweight="bold", color=CLASS_COLORS[cls_name])
        y += line_h * 0.85
        # Member rows
        members = []
        for _, r in sub.iterrows():
            fid = int(r["family_id"])
            short = short_overrides.get(fid, "")
            rho = float(r["rho_primary"])
            ssn = int(r["sens_same_sign_n"])
            members.append(f"F{fid:2d} {short:16s} (\u03c1={rho:+.2f}, {ssn}/4)")
        # Wrap into chunks of items_per_line
        for k in range(0, len(members), items_per_line):
            chunk = members[k:k + items_per_line]
            txt = "   ".join(chunk)
            ax_fams.text(0.04, y, txt, fontsize=6.0, color="#1A1A1A", family="monospace")
            y += line_h * 0.75
        y += line_h * 0.4

    fig.suptitle(
        "Fig 4c  Chronicity classification (GSE41177 LAA, 16 families)",
        x=0.05, y=0.96, ha="left", fontsize=8.8, fontweight="bold",
    )
    # Compact footnote (full version belongs to the manuscript Figure Legend).
    fig.text(
        0.05, 0.10,
        "GSE41177 LAA n = 12.  Effect-size partition under power limitation.  "
        "Strong / Weak / Stable defined by |\u03c1| threshold (see headers above).  "
        "No family reaches |\u03c1| \u2265 0.4 + p < 0.05.",
        fontsize=5.9, color="#444444",
    )
    fig.text(
        0.05, 0.075,
        "Strong-trend candidates F2 + F9 (and Weak-trend F15 for cross-figure link with Fig 3D) are detailed in Fig 4d / 4e / 4f.  "
        "See Figure Legend for full statistics.",
        fontsize=5.9, color="#444444",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
