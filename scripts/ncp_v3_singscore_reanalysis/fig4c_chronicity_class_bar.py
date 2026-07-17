#!/usr/bin/env python3
"""Corrected Fig 4C — Strong / Weak / Stable class summary (neutral greys)."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
OUT_BASE = FIGURES / "fig4c_chronicity_class_bar"
WIDTH_MM, HEIGHT_MM = 178.0, 35.0

CHRONICITY_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"

CLASS_ORDER = ["Strong_trend", "Weak_trend", "Stable"]
CLASS_LABELS = {
    "Strong_trend": "Strong-trend\n(|ρ| ≥ 0.4)",
    "Weak_trend": "Weak-trend\n(0.2 ≤ |ρ| < 0.4)",
    "Stable": "Stable\n(|ρ| < 0.2)",
}
CLASS_COLORS = {
    "Strong_trend": "#4A4A4A",
    "Weak_trend": "#8A8A8A",
    "Stable": "#C8C8C8",
}


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 6.5,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def save_outputs(fig: plt.Figure) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BASE.with_suffix(".pdf"), bbox_inches=None, pad_inches=0.02)
    fig.savefig(OUT_BASE.with_suffix(".svg"), bbox_inches=None, pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    set_style()
    cls = pd.read_csv(CHRONICITY_TSV, sep="\t")
    labels = pd.read_csv(LABELS_TSV, sep="\t").set_index("family_id")
    counts = cls["chronicity_class"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int)
    assert int(counts["Strong_trend"]) == 2
    assert int(counts["Weak_trend"]) == 6
    assert int(counts["Stable"]) == 8
    total = int(counts.sum())

    mm = 1.0 / 25.4
    fig = plt.figure(figsize=(WIDTH_MM * mm, HEIGHT_MM * mm))
    gs = fig.add_gridspec(1, 3, width_ratios=[2.5, 1.0, 3.5],
                          left=0.05, right=0.98, top=0.82, bottom=0.22, wspace=0.06)
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])
    ax_fams = fig.add_subplot(gs[0, 2])

    x_left = 0.0
    for cls_name in CLASS_ORDER:
        n = int(counts[cls_name])
        ax_bar.barh(0.5, n, left=x_left, height=0.6, color=CLASS_COLORS[cls_name],
                    edgecolor="white", linewidth=1.2)
        ax_bar.text(x_left + n / 2.0, 0.5, str(n), ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold")
        x_left += n
    ax_bar.set_xlim(0, total)
    ax_bar.set_ylim(0, 1)
    ax_bar.set_yticks([])
    ax_bar.set_xticks(range(0, total + 1, 2))
    ax_bar.set_xlabel("families", fontsize=6.0)
    for spine in ("top", "right", "left"):
        ax_bar.spines[spine].set_visible(False)

    ax_leg.axis("off")
    y = 0.1
    for cls_name in CLASS_ORDER:
        ax_leg.add_patch(mpatches.Rectangle((0.0, y), 0.15, 0.12,
                                            facecolor=CLASS_COLORS[cls_name], edgecolor="none"))
        ax_leg.text(0.2, y + 0.06, CLASS_LABELS[cls_name], fontsize=5.5, va="center")
        y += 0.28

    ax_fams.axis("off")
    ax_fams.set_xlim(0, 1)
    ax_fams.set_ylim(1, 0)
    y = 0.05
    line_h = 0.14
    for cls_name in CLASS_ORDER:
        sub = cls[cls["chronicity_class"] == cls_name].sort_values(
            "rho_primary", key=lambda s: s.abs(), ascending=False,
        )
        members = []
        for _, row in sub.iterrows():
            fid = int(row["family_id"])
            short = str(labels.loc[fid, "short_label"])[:18]
            members.append(f"F{fid} {short} (ρ={float(row['rho_primary']):+.2f})")
        ax_fams.text(0.0, y, f"{cls_name.replace('_', '-')} (n={len(sub)}):",
                     fontsize=5.8, fontweight="bold", color=CLASS_COLORS[cls_name])
        ax_fams.text(0.22, y, "  ".join(members), fontsize=5.4, color="#333333", family="monospace")
        y += line_h

    fig.text(0.05, 0.92, "Fig 4c  Chronicity class summary (GSE41177 LAA primary, 16 families)",
             fontsize=7.2, fontweight="bold", color="#1A1A1A")
    fig.text(0.05, 0.08,
             "Classes use |ρ_LAA| only.  L3 layer for Fig 5: chronicity_class ≠ Stable.  "
             "Exemplar donors in Fig 4d (F9, F2, F15).",
             fontsize=5.5, color="#444444")

    save_outputs(fig)
    print(f"wrote {OUT_BASE.with_suffix('.pdf').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
