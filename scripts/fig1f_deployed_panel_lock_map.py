#!/usr/bin/env python3
"""
Figure 1F - baseline v3 panel lock map.

Purpose:
    Freeze the baseline v3 object by showing program-level topology,
    score confidence floor, and perturbation robustness in one compact panel.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label
from plot_style import COL_DOUBLE, FG_COLOR, LA_COLOR, LINE_COLOR, MUTED_COLOR, RA_COLOR, SUPPORT_COLOR, set_nature_style

set_nature_style()
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

PANEL_TSV = ROOT / "results" / "innovation_ncp_v3_panel.tsv"
FREQ_TSV = ROOT / "results" / "innovation_ncp_bootstrap_frequency.tsv"
OUT_TSV = ROOT / "results" / "fig1f_deployed_panel_lock.tsv"
OUT_PDF = ROOT / "figures" / "fig1f_deployed_panel_map.pdf"

TIER_COLORS = {
    "core": "#E9F2EC",
    "supporting": "#F3EFE7",
    "exploratory": "#EFECEF",
}


def main() -> None:
    panel = pd.read_csv(PANEL_TSV, sep="\t").sort_values("program_order").reset_index(drop=True)
    freq = pd.read_csv(FREQ_TSV, sep="\t")[["pathway", "inclusion_freq", "robust_tier"]]

    df = panel.merge(freq, on="pathway", how="left")
    df["inclusion_freq"] = df["inclusion_freq"].fillna(0.0)
    df["robust_tier"] = df["robust_tier"].fillna("sensitive")
    df.to_csv(OUT_TSV, sep="\t", index=False)

    n = len(df)
    y = np.arange(n)
    n_la = int((df["direction"] == "LA-biased").sum())
    n_ra = int((df["direction"] == "RA-biased").sum())

    fig = plt.figure(figsize=(COL_DOUBLE * 1.10, COL_DOUBLE * 0.95))
    gs = fig.add_gridspec(
        1,
        4,
        width_ratios=[0.16, 3.2, 2.45, 1.45],
        wspace=0.03,
        top=0.92,
        bottom=0.10,
        left=0.03,
        right=0.99,
    )

    # Direction strip
    ax_dir = fig.add_subplot(gs[0, 0])
    ax_dir.set_xlim(0, 1)
    ax_dir.set_ylim(n - 0.5, -0.5)
    ax_dir.axis("off")
    for i, r in df.iterrows():
        col = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        ax_dir.add_patch(mpatches.Rectangle((0, i - 0.5), 1, 1, facecolor=col, edgecolor="white", linewidth=0.2))
    if 0 < n_la < n:
        ax_dir.axhline(n_la - 0.5, color=LINE_COLOR, lw=0.8)

    # Labels and topology
    ax_lbl = fig.add_subplot(gs[0, 1], sharey=ax_dir)
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(n - 0.5, -0.5)
    ax_lbl.axis("off")
    for i, r in df.iterrows():
        bg = TIER_COLORS.get(str(r["innovation_tier"]), "#F6F6F6")
        ax_lbl.add_patch(
            mpatches.Rectangle(
                (0.0, i - 0.5),
                1.0,
                1.0,
                facecolor=bg,
                edgecolor="white",
                linewidth=0.2,
                zorder=0,
            )
        )
        col = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        ax_lbl.text(
            0.01,
            i,
            f"{r['cluster_id']}",
            fontsize=7.1,
            color=col,
            ha="left",
            va="center",
            fontweight="bold",
        )
        ax_lbl.text(
            0.20,
            i,
            short_pathway_label(str(r["pathway"]), 42),
            fontsize=7.1,
            color=FG_COLOR,
            ha="left",
            va="center",
        )
    if 0 < n_la < n:
        ax_lbl.axhline(n_la - 0.5, color=LINE_COLOR, lw=0.8)

    # Score/robustness tracks
    ax_trk = fig.add_subplot(gs[0, 2], sharey=ax_dir)
    ax_trk.set_xlim(0, 1.02)
    ax_trk.set_ylim(n - 0.5, -0.5)
    ax_trk.axvline(0.5, color="#9A9A9A", lw=0.7, ls="--")
    ax_trk.axvline(0.8, color=SUPPORT_COLOR, lw=0.8, ls="--")
    for i, r in df.iterrows():
        col = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        floor = float(r["confidence_floor"])
        score = float(r["consensus_score"])
        freqv = float(r["inclusion_freq"])
        # confidence segment
        ax_trk.plot([floor, score], [i, i], color="#7E7E7E", lw=1.4, solid_capstyle="round", zorder=2)
        # score point
        ax_trk.scatter(score, i, s=24, color=col, edgecolor="white", linewidth=0.3, zorder=3)
        # robustness marker (slight vertical offset)
        ax_trk.scatter(freqv, i + 0.18, s=20, marker="s", color="#4D4D4D", edgecolor="white", linewidth=0.25, zorder=3)

    ax_trk.set_yticks([])
    ax_trk.set_xlabel("score / robustness track", fontsize=7.2)
    ax_trk.tick_params(axis="x", labelsize=7.0)
    for s in ("top", "right"):
        ax_trk.spines[s].set_visible(False)
    if 0 < n_la < n:
        ax_trk.axhline(n_la - 0.5, color=LINE_COLOR, lw=0.8)
    ax_trk.text(0.5, -1.1, "stable", fontsize=6.3, color="#7A7A7A", ha="center", va="bottom")
    ax_trk.text(0.8, -1.1, "very stable", fontsize=6.3, color=SUPPORT_COLOR, ha="center", va="bottom")
    ax_trk.set_title("confidence floor -> consensus score; square = inclusion freq", loc="left", fontsize=7.4, pad=2)

    # Right summary block
    ax_sum = fig.add_subplot(gs[0, 3], sharey=ax_dir)
    ax_sum.axis("off")
    ax_sum.set_xlim(0, 1)
    ax_sum.set_ylim(1, 0)
    core_n = int((df["innovation_tier"] == "core").sum())
    sup_n = int((df["innovation_tier"] == "supporting").sum())
    exp_n = int((df["innovation_tier"] == "exploratory").sum())
    vstable = int((df["inclusion_freq"] >= 0.8).sum())
    stable = int((df["inclusion_freq"] >= 0.5).sum())
    med_score = float(df["consensus_score"].median())
    med_floor = float(df["confidence_floor"].median())

    ax_sum.text(0.00, 0.05, "baseline object lock", fontsize=8.4, fontweight="bold", color=FG_COLOR, ha="left")
    ax_sum.text(0.00, 0.13, f"K=30\nLA={n_la}, RA={n_ra}", fontsize=7.3, color=MUTED_COLOR, ha="left", va="top")

    ax_sum.text(0.00, 0.30, "innovation tier", fontsize=7.2, fontweight="bold", color=FG_COLOR, ha="left")
    ax_sum.text(0.02, 0.36, f"core: {core_n}", fontsize=7.0, color=MUTED_COLOR, ha="left")
    ax_sum.text(0.02, 0.41, f"supporting: {sup_n}", fontsize=7.0, color=MUTED_COLOR, ha="left")
    ax_sum.text(0.02, 0.46, f"exploratory: {exp_n}", fontsize=7.0, color=MUTED_COLOR, ha="left")

    ax_sum.text(0.00, 0.58, "perturbation stability", fontsize=7.2, fontweight="bold", color=FG_COLOR, ha="left")
    ax_sum.text(0.02, 0.64, f"freq >= 0.8: {vstable}", fontsize=7.0, color=MUTED_COLOR, ha="left")
    ax_sum.text(0.02, 0.69, f"freq >= 0.5: {stable}", fontsize=7.0, color=MUTED_COLOR, ha="left")

    ax_sum.text(0.00, 0.80, "score summary", fontsize=7.2, fontweight="bold", color=FG_COLOR, ha="left")
    ax_sum.text(0.02, 0.86, f"median consensus: {med_score:.3f}", fontsize=7.0, color=MUTED_COLOR, ha="left")
    ax_sum.text(0.02, 0.91, f"median floor: {med_floor:.3f}", fontsize=7.0, color=MUTED_COLOR, ha="left")

    fig.text(
        0.03,
        0.955,
        "f  Baseline panel lock map and metadata",
        fontsize=9.8,
        fontweight="bold",
        ha="left",
        va="center",
        color=FG_COLOR,
    )
    fig.text(
        0.03,
        0.05,
        "Rows are frozen by program_order in innovation_ncp_v3_panel.tsv; no downstream re-ranking is allowed.",
        fontsize=6.4,
        color=MUTED_COLOR,
        ha="left",
        va="center",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)

    print(f"wrote {OUT_TSV.relative_to(ROOT)}")
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
