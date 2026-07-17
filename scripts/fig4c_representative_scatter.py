#!/usr/bin/env python3
"""Figure 4C  Representative duration scatter plots."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

SCORES_TSV = ROOT / "Project_P4" / "GSE41177" / "program_scores_primary.tsv"
META_TSV   = ROOT / "Project_P4" / "GSE41177" / "sample_metadata.tsv"
LAA_TSV    = ROOT / "Project_P4" / "GSE41177" / "results_primary" \
                / "duration_trend_laa_primary.tsv"
FIG_OUT    = ROOT / "figures" / "fig4c_representative_scatter.pdf"

LA_COLOR = "#D6604D"
RA_COLOR = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"
PRGN_PURPLE = "#762A83"
PRGN_GREEN  = "#1B7837"

# Top-3 program ids (locked to current LAA primary p < 0.06 hits)
TARGETS = [
    {"program_id": "HALLMARK__HALLMARK_OXIDATIVE_PHOSPHORYLATION",
     "label": "OXPHOS",
     "family": "F7",
     "direction": "LA-biased"},
    {"program_id": "REACTOME__REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT",
     "label": "Mito Ca²⁺ transport",
     "family": "F2",
     "direction": "LA-biased"},
    {"program_id": "REACTOME__REACTOME_CLASS_B_2_SECRETIN_FAMILY_RECEPTORS",
     "label": "Class B/2 secretin receptors",
     "family": "F9",
     "direction": "RA-biased"},
]


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "axes.titlesize":    9,
        "axes.labelsize":    8,
        "xtick.labelsize":   7,
        "ytick.labelsize":   7,
        "legend.fontsize":   6.5,
        "axes.linewidth":    0.6,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def main() -> None:
    set_publication_style()

    scores = pd.read_csv(SCORES_TSV, sep="\t")
    laa_stats = pd.read_csv(LAA_TSV, sep="\t").set_index("program_id")

    # restrict to LAA scores
    scores = scores[scores["tissue"] == "LAA"].copy()

    n = len(TARGETS)
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(180 * MM, 65 * MM))
    grid = fig.add_gridspec(
        nrows=1, ncols=n, wspace=0.30,
        left=0.07, right=0.98, top=0.85, bottom=0.18)

    fig.suptitle(
        "Fig 4C  Representative AF-duration trends in LAA "
        "(GSE41177, n = 12 AF + 7 SR; p < 0.06 program hits)",
        fontsize=8.8, fontweight="bold", y=0.97, x=0.02, ha="left")

    for idx, tgt in enumerate(TARGETS):
        ax = fig.add_subplot(grid[idx])
        ps = scores[scores["program_id"] == tgt["program_id"]].copy()
        if ps.empty:
            ax.text(0.5, 0.5, f"NO DATA: {tgt['label']}",
                     ha="center", va="center",
                     transform=ax.transAxes)
            continue

        af   = ps[(ps["rhythm"] == "AF") & ps["duration_months"].notna()] \
                  .sort_values("duration_months")
        sr   = ps[ps["rhythm"] == "SR"]

        x_af = af["duration_months"].to_numpy(dtype=float)
        y_af = af["score"].to_numpy(dtype=float)

        col = LA_COLOR if tgt["direction"] == "LA-biased" else RA_COLOR

        # SR reference markers at x = 0
        if len(sr) > 0:
            jitter = (np.random.default_rng(7).normal(0, 1.2, len(sr)))
            ax.scatter(jitter, sr["score"],
                          marker="o", s=22, facecolor="white",
                          edgecolor=GREY_MID, linewidth=0.6,
                          alpha=0.7, zorder=2,
                          label=f"SR  n = {len(sr)}")

        # AF dots
        ax.scatter(x_af, y_af,
                      marker="o", s=30, facecolor=col,
                      edgecolor="black", linewidth=0.4,
                      alpha=0.9, zorder=3,
                      label=f"AF  n = {len(x_af)}")

        # Linear fit (least squares on AF only)
        if len(x_af) >= 3 and np.std(x_af) > 0:
            slope, intercept = np.polyfit(x_af, y_af, 1)
            xs = np.linspace(min(x_af), max(x_af), 80)
            ys = slope * xs + intercept
            # CI band via residual std
            resid = y_af - (slope * x_af + intercept)
            sse = np.sum(resid ** 2)
            mse = sse / max(1, len(x_af) - 2)
            x_mean = x_af.mean()
            sx = np.sqrt(np.sum((x_af - x_mean) ** 2))
            se_y = np.sqrt(mse * (1 / len(x_af) + (xs - x_mean) ** 2 / sx ** 2))
            ax.fill_between(xs, ys - 1.96 * se_y, ys + 1.96 * se_y,
                              color=col, alpha=0.15, zorder=1.5,
                              edgecolor="none")
            ax.plot(xs, ys, color=col, lw=1.2, zorder=2.5)

        # Stats annotation
        if tgt["program_id"] in laa_stats.index:
            srow = laa_stats.loc[tgt["program_id"]]
            rho = float(srow["statistic"])
            p   = float(srow["p_perm"])
            sig = ""
            if p < 0.05: sig = " *"
            elif p < 0.10: sig = " ."
            ann = (f"ρ = {rho:+.3f}{sig}\n"
                     f"p_perm = {p:.3f}\n"
                     f"n_AF = {int(srow['n_samples'])}")
            ax.text(0.05, 0.95, ann,
                     transform=ax.transAxes,
                     ha="left", va="top",
                     fontsize=6.5, color=GREY_DARK,
                     bbox=dict(boxstyle="round,pad=0.25",
                                  facecolor="white",
                                  edgecolor=GREY_LIGHT, linewidth=0.5))

        ax.axvline(0, color=GREY_LIGHT, lw=0.5, zorder=1)
        ax.set_xlim(-5, max(x_af.max() if len(x_af) else 100, 100) * 1.05)
        ax.set_xlabel("AF duration  (months)", fontsize=7.5)
        ax.set_ylabel("singscore", fontsize=7.5)

        # Title : F# + label, coloured by direction
        ax.set_title(
            f"{tgt['family']}  {tgt['label']}",
            loc="left", fontsize=8.0, fontweight="bold",
            color=col, pad=4)
        ax.legend(loc="lower left", fontsize=6.0, frameon=False,
                     handletextpad=0.4, borderaxespad=0.3)

    fig.text(0.02, 0.020,
              "AF dots colour = baseline chamber direction (blue LA-biased, "
              "red RA-biased).  Line = least-squares fit on AF only;  "
              "shaded band = 95 % CI of fitted line.   "
              "SR samples plotted as grey hollow markers near x = 0  "
              "(visual reference; not fitted).",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
