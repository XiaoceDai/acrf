#!/usr/bin/env python3
"""
Figure 1B (algorithm storyline) - baseline aCAI effect + concordance view.

Replaces heatmap-heavy rendering with a method-centric panel:
  1) donor-level sign profile (dots)
  2) program-level effect size (mean with min-max and IQR)
  3) sign concordance summary (k/8)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
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

IN_TSV = ROOT / "results" / "fig1b_baseline_donor_aCAI.tsv"
OUT_TSV = ROOT / "results" / "fig1b_baseline_effect_concordance.tsv"
OUT_PDF = ROOT / "figures" / "fig1b_baseline_effect_concordance.pdf"


def sort_programs(df: pd.DataFrame) -> pd.DataFrame:
    la = df[df["direction"] == "LA-biased"].sort_values("mean_aCAI", ascending=False)
    ra = df[df["direction"] == "RA-biased"].sort_values("mean_aCAI", ascending=True)
    other = df[~df["direction"].isin(["LA-biased", "RA-biased"])].sort_values("mean_aCAI", ascending=False)
    return pd.concat([la, ra, other], ignore_index=True)


def main() -> None:
    df = pd.read_csv(IN_TSV, sep="\t")
    donor_cols = [c for c in df.columns if c.startswith("aCAI_")]
    if not donor_cols:
        raise ValueError("No donor columns found in fig1b_baseline_donor_aCAI.tsv")

    if "mean_aCAI" not in df.columns:
        df["mean_aCAI"] = df[donor_cols].mean(axis=1)
    if "sign_concordance" not in df.columns:
        df["sign_concordance"] = np.sign(df[donor_cols]).eq(np.sign(df["mean_aCAI"]), axis=0).sum(axis=1).astype(int)

    df["q25"] = df[donor_cols].quantile(0.25, axis=1)
    df["q75"] = df[donor_cols].quantile(0.75, axis=1)
    df["vmin"] = df[donor_cols].min(axis=1)
    df["vmax"] = df[donor_cols].max(axis=1)
    df["abs_mean"] = df["mean_aCAI"].abs()
    # Naive y~chamber estimate on balanced paired healthy data:
    # mean(LA) - mean(RA) is algebraically identical to mean paired difference.
    # We keep it as an explicit column for visual comparison with paired estimate.
    df["naive_chamber_est"] = df[donor_cols].mean(axis=1)

    plot_df = sort_programs(df)
    top_idx = set(plot_df["abs_mean"].nlargest(6).index.tolist())
    plot_df.to_csv(
        OUT_TSV,
        sep="\t",
        index=False,
        columns=[
            "program_order",
            "pathway",
            "direction",
            "mean_aCAI",
            "q25",
            "q75",
            "vmin",
            "vmax",
            "naive_chamber_est",
            "sign_concordance",
        ]
        + donor_cols,
    )

    n = len(plot_df)
    n_la = int((plot_df["direction"] == "LA-biased").sum())
    y = np.arange(n)
    donor_names = [c.replace("aCAI_", "") for c in donor_cols]

    fig = plt.figure(figsize=(COL_DOUBLE * 1.04, COL_DOUBLE * 0.96))
    gs = fig.add_gridspec(
        1,
        5,
        width_ratios=[0.18, 3.0, 1.7, 2.55, 0.9],
        wspace=0.03,
        top=0.92,
        bottom=0.12,
        left=0.03,
        right=0.99,
    )

    # Column 0: direction strip
    ax_dir = fig.add_subplot(gs[0, 0])
    ax_dir.set_xlim(0, 1)
    ax_dir.set_ylim(n - 0.5, -0.5)
    ax_dir.axis("off")
    for i, (_, r) in enumerate(plot_df.iterrows()):
        face = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        ax_dir.add_patch(mpatches.Rectangle((0, i - 0.5), 1, 1, facecolor=face, edgecolor="white", linewidth=0.2))
    if 0 < n_la < n:
        ax_dir.axhline(n_la - 0.5, color=LINE_COLOR, lw=0.8)

    # Column 1: program labels
    ax_lbl = fig.add_subplot(gs[0, 1], sharey=ax_dir)
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(n - 0.5, -0.5)
    ax_lbl.axis("off")
    for i, (_, r) in enumerate(plot_df.iterrows()):
        color = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        text = short_pathway_label(r["pathway"], 38)
        if i in top_idx:
            text = f"* {text}"
        ax_lbl.text(
            0.0,
            i,
            text,
            fontsize=7.8,
            color=color,
            ha="left",
            va="center",
            fontweight="bold" if i in top_idx else "normal",
        )
    if 0 < n_la < n:
        ax_lbl.axhline(n_la - 0.5, color=LINE_COLOR, lw=0.8)

    # Column 2: donor-level sign profile dots
    ax_dot = fig.add_subplot(gs[0, 2], sharey=ax_lbl)
    ax_dot.set_xlim(-0.5, len(donor_cols) - 0.5)
    ax_dot.set_ylim(n - 0.5, -0.5)
    for i, row in plot_df.iterrows():
        msign = np.sign(row["mean_aCAI"])
        for j, c in enumerate(donor_cols):
            v = float(row[c])
            dsign = np.sign(v)
            dot_color = LA_COLOR if v >= 0 else RA_COLOR
            ax_dot.scatter(j, i, s=22, color=dot_color, edgecolor="white", linewidth=0.25, zorder=3)
            if msign != 0 and dsign != 0 and dsign != msign:
                ax_dot.scatter(j, i, s=15, marker="x", color=FG_COLOR, linewidth=0.6, zorder=4)
    ax_dot.set_xticks(np.arange(len(donor_cols)))
    ax_dot.set_xticklabels(donor_names, rotation=45, ha="right", fontsize=6.8)
    ax_dot.set_yticks([])
    ax_dot.set_xlabel("donor sign profile", fontsize=7.2)
    for spine in ("top", "right"):
        ax_dot.spines[spine].set_visible(False)
    if 0 < n_la < n:
        ax_dot.axhline(n_la - 0.5, color=LINE_COLOR, lw=0.8)

    # Column 3: effect size distribution summary
    ax_eff = fig.add_subplot(gs[0, 3], sharey=ax_lbl)
    ax_eff.set_ylim(n - 0.5, -0.5)
    ax_eff.axvline(0, color=LINE_COLOR, lw=0.8, zorder=0)
    for i, row in plot_df.iterrows():
        col = LA_COLOR if row["direction"] == "LA-biased" else RA_COLOR
        ax_eff.plot([row["vmin"], row["vmax"]], [i, i], color="#C7C7C7", lw=0.9, zorder=1)
        ax_eff.plot([row["q25"], row["q75"]], [i, i], color="#7F7F7F", lw=2.0, zorder=2)
        # Naive estimate: hollow circle (larger ring so overlap is still visible).
        ax_eff.scatter(
            row["naive_chamber_est"],
            i,
            s=44,
            facecolor="none",
            edgecolor="#4A4A4A",
            linewidth=0.75,
            zorder=3,
        )
        # Paired estimate: filled circle.
        ax_eff.scatter(row["mean_aCAI"], i, s=20, color=col, edgecolor="white", linewidth=0.3, zorder=4)
        if i in top_idx:
            ax_eff.scatter(row["mean_aCAI"], i, s=48, facecolor="none", edgecolor="#6A5A2B", linewidth=0.8, zorder=5)
    xlim = max(0.15, float(np.nanmax(np.abs(plot_df[donor_cols].values))) * 1.08)
    ax_eff.set_xlim(-xlim, xlim)
    ax_eff.set_yticks([])
    ax_eff.set_xlabel("aCAI effect size (mean, IQR, range)", fontsize=7.2)
    ax_eff.tick_params(axis="x", labelsize=7.0)
    for spine in ("top", "right"):
        ax_eff.spines[spine].set_visible(False)
    if 0 < n_la < n:
        ax_eff.axhline(n_la - 0.5, color=LINE_COLOR, lw=0.8)
    ax_eff.set_title(
        "b  Baseline chamber prior: effect size and donor concordance",
        loc="left",
        fontsize=9.8,
        fontweight="bold",
        pad=6,
    )

    # Column 4: sign concordance
    ax_sc = fig.add_subplot(gs[0, 4], sharey=ax_lbl)
    ax_sc.set_xlim(0, 1)
    ax_sc.set_ylim(n - 0.5, -0.5)
    ax_sc.axis("off")
    ax_sc.set_title("k/8", loc="left", fontsize=7.6, fontweight="bold", pad=3)
    for i, row in plot_df.iterrows():
        sc = int(row["sign_concordance"])
        color = SUPPORT_COLOR if sc == len(donor_cols) else MUTED_COLOR
        fw = "bold" if sc == len(donor_cols) else "normal"
        ax_sc.text(0.5, i, f"{sc}/{len(donor_cols)}", ha="center", va="center", fontsize=7.0, color=color, fontweight=fw)
    if 0 < n_la < n:
        ax_sc.axhline(n_la - 0.5, color=LINE_COLOR, lw=0.8)

    # Legend for effect-size summary primitives
    legend_handles = [
        mlines.Line2D([], [], color="#C7C7C7", lw=1.0, label="range (min-max)"),
        mlines.Line2D([], [], color="#7F7F7F", lw=2.2, label="IQR (Q1-Q3)"),
        mlines.Line2D([], [], marker="o", markersize=5.0, linestyle="None", markerfacecolor="none", markeredgecolor="#4A4A4A", label="naive y~chamber"),
        mlines.Line2D([], [], marker="o", markersize=4.4, linestyle="None", markerfacecolor=FG_COLOR, markeredgecolor="white", label="paired aCAI"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.57, 0.045),
        ncol=4,
        frameon=False,
        fontsize=6.4,
        handlelength=1.7,
        columnspacing=0.9,
    )

    fig.text(
        0.03,
        0.025,
        "In this balanced healthy paired cohort, naive and paired estimates can overlap; hollow ring is kept visible around filled point.",
        fontsize=6.5,
        color=MUTED_COLOR,
        ha="left",
        va="bottom",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
