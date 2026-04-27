#!/usr/bin/env python3
"""
Figure 1E - baseline panel stability under healthy resampling perturbation.

Uses bootstrap-like perturbation outputs to summarize how stable the
v3 baseline panel remains across score-threshold and weight perturbations.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
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

IN_FREQ = ROOT / "results" / "innovation_ncp_bootstrap_frequency.tsv"
IN_SUM = ROOT / "results" / "innovation_ncp_bootstrap_summary.tsv"
IN_PANEL = ROOT / "results" / "innovation_ncp_v3_panel.tsv"

OUT_FIG1E = ROOT / "figures" / "fig1e_v3_robustness_stability.pdf"
OUT_FIGS14 = ROOT / "figures" / "figS14_v3_bootstrap_stability.pdf"


def main() -> None:
    freq = pd.read_csv(IN_FREQ, sep="\t")
    summ = pd.read_csv(IN_SUM, sep="\t").iloc[0]
    panel = pd.read_csv(IN_PANEL, sep="\t")

    panel_paths = set(panel["pathway"].astype(str))
    freq = freq.copy()
    freq["is_panel"] = freq["pathway"].astype(str).isin(panel_paths)

    panel_df = (
        panel[["program_order", "pathway", "direction"]]
        .merge(freq[["pathway", "inclusion_freq", "robust_tier"]], on="pathway", how="left")
        .sort_values("program_order")
        .reset_index(drop=True)
    )
    panel_df["inclusion_freq"] = panel_df["inclusion_freq"].fillna(0.0)
    panel_df["robust_tier"] = panel_df["robust_tier"].fillna("sensitive")

    nonpanel_df = freq.loc[~freq["is_panel"]].copy()

    fig = plt.figure(figsize=(COL_DOUBLE * 1.08, COL_DOUBLE * 0.92))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[2.8, 1.4],
        height_ratios=[1.0, 1.0],
        wspace=0.22,
        hspace=0.28,
        left=0.05,
        right=0.99,
        top=0.90,
        bottom=0.10,
    )

    # Left: panel-member frequencies
    ax = fig.add_subplot(gs[:, 0])
    d = panel_df.sort_values("inclusion_freq", ascending=True).reset_index(drop=True)
    y = np.arange(len(d))
    colors = [LA_COLOR if x == "LA-biased" else RA_COLOR for x in d["direction"]]
    ax.barh(y, d["inclusion_freq"], color=colors, edgecolor="white", linewidth=0.25)
    ax.axvline(0.5, color="#8B8B8B", lw=0.7, ls="--")
    ax.axvline(0.8, color=SUPPORT_COLOR, lw=0.8, ls="--")
    ax.text(0.5, -1.2, "stable", ha="center", va="bottom", fontsize=6.4, color="#666")
    ax.text(0.8, -1.2, "very stable", ha="center", va="bottom", fontsize=6.4, color=SUPPORT_COLOR)
    ax.set_xlim(0, 1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([short_pathway_label(p, 38) for p in d["pathway"]], fontsize=6.5)
    ax.set_xlabel("Inclusion frequency across perturbation runs", fontsize=7.2)
    ax.set_title("panel-member baseline stability profile (K=30)", loc="left", fontsize=8.6, fontweight="bold", pad=4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", labelsize=7.0)

    # Right top: tier composition among panel members
    ax_tier = fig.add_subplot(gs[0, 1])
    tier_order = ["very-stable", "stable", "sensitive"]
    tier_colors = [SUPPORT_COLOR, "#9C9C9C", "#C8A076"]
    counts = [int((panel_df["robust_tier"] == t).sum()) for t in tier_order]
    ax_tier.bar(np.arange(3), counts, color=tier_colors, width=0.62)
    ax_tier.set_xticks(np.arange(3), ["very\nstable", "stable", "sensitive"], fontsize=6.7)
    ax_tier.set_ylabel("n pathways", fontsize=7.0)
    ax_tier.set_ylim(0, max(1, len(panel_df)))
    ax_tier.set_title("tier composition", loc="left", fontsize=8.4, fontweight="bold", pad=3)
    for i, c in enumerate(counts):
        ax_tier.text(i, c + 0.6, str(c), ha="center", va="bottom", fontsize=6.8, color=FG_COLOR)
    for s in ("top", "right"):
        ax_tier.spines[s].set_visible(False)
    ax_tier.tick_params(axis="y", labelsize=6.6)

    # Right bottom: panel vs non-panel frequency distribution
    ax_dist = fig.add_subplot(gs[1, 1])
    panel_vals = panel_df["inclusion_freq"].to_numpy(float)
    nonpanel_vals = nonpanel_df["inclusion_freq"].to_numpy(float)
    bins = np.linspace(0, 1, 21)
    ax_dist.hist(nonpanel_vals, bins=bins, color="#D8D8D8", alpha=1.0, label="non-panel", edgecolor="white")
    ax_dist.hist(panel_vals, bins=bins, color="#8F9A7A", alpha=0.9, label="panel", edgecolor="white")
    ax_dist.set_xlabel("inclusion frequency", fontsize=7.0)
    ax_dist.set_ylabel("count", fontsize=7.0)
    ax_dist.set_title("panel vs background", loc="left", fontsize=8.4, fontweight="bold", pad=3)
    ax_dist.legend(frameon=False, fontsize=6.5, loc="upper right")
    for s in ("top", "right"):
        ax_dist.spines[s].set_visible(False)
    ax_dist.tick_params(labelsize=6.6)

    fig.text(
        0.05,
        0.95,
        "e  baseline panel stability under healthy resampling perturbation",
        fontsize=9.8,
        fontweight="bold",
        ha="left",
        va="center",
        color=FG_COLOR,
    )
    fig.text(
        0.73,
        0.88,
        f"runs={int(summ['n_runs'])} | target K={int(summ['target_k'])}\n"
        f"mean RA/run={summ['mean_ra_per_run']:.2f} | very-stable={int(summ['n_path_very_stable_freq_ge_0_8'])}",
        fontsize=6.6,
        color=MUTED_COLOR,
        ha="left",
        va="top",
    )

    OUT_FIG1E.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG1E, bbox_inches=None, pad_inches=0.04)
    fig.savefig(OUT_FIGS14, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)

    print(f"wrote {OUT_FIG1E.relative_to(ROOT)}")
    print(f"wrote {OUT_FIGS14.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
