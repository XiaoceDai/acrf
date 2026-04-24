#!/usr/bin/env python3
"""
Figure S11 — Method-validity P5 (simulation + bootstrap + classical concordance).

S11a — simulation power: mean pooled g vs true effect size + detection rate
       (|g| > 0.3) on cohort-matched (5+5, 7+6, 13+15) synthetic data.
S11b — bootstrap CI for STRUCTURAL / METABOLIC pooled g vs the locked
       Hartung–Knapp CI from Fig 6A.
S11c — classical limma per-chamber AF vs SR direction concordance with
       the framework's per-program direction (per cohort).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle

from plot_style import COL_DOUBLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "figures" / "figS11_method_p5.pdf"

AXIS_COL = {"STRUCTURAL": "#B35806", "METABOLIC": "#542788"}


def panel_S11a(ax, sim):
    """Simulation: mean pooled g & detection rate vs true effect size."""
    eff = sim["effect_size"].astype(float)
    mean_g = sim["mean_pooled_g"].astype(float)
    sd_g = sim["sd_pooled_g"].astype(float)
    rate = sim["detect_rate"].astype(float)

    ax.errorbar(eff, mean_g, yerr=sd_g,
                fmt="o-", color="#1A1A1A", markersize=4,
                ecolor="#888888", elinewidth=0.6, capsize=2)
    ax.axhline(0, color="#888888", lw=0.4, ls=(0, (1.5, 1.5)))
    ax.set_xlabel("true rewiring effect (aCAI$_{rel}$ units)",
                  fontsize=6.0, labelpad=2)
    ax.set_ylabel("recovered pooled g  ± SD\n(REML / Hartung–Knapp)",
                  fontsize=6.0, labelpad=2, linespacing=1.3)
    ax.set_title(
        "S11a — recovery of axis pooled g (cohort-matched synthetic data)",
        fontsize=6.4, pad=3, fontweight="bold", loc="left")
    ax.tick_params(axis="both", labelsize=5.4, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5); ax.spines["bottom"].set_linewidth(0.5)

    # Right axis: detection rate
    ax2 = ax.twinx()
    ax2.plot(eff, rate, color="#1B7837", lw=1.0, marker="s", markersize=3.5,
             markerfacecolor="#1B7837", markeredgecolor="white",
             markeredgewidth=0.5)
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("detection rate  (|g| > 0.3, 200 reps)",
                   color="#1B7837", fontsize=6.0, labelpad=2)
    ax2.tick_params(axis="y", labelsize=5.4, length=2, colors="#1B7837")
    ax2.spines["top"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.spines["right"].set_linewidth(0.5)
    # legend
    handles = [
        Line2D([], [], marker="o", linestyle="-", color="#1A1A1A",
               markersize=4, label="recovered g"),
        Line2D([], [], marker="s", linestyle="-", color="#1B7837",
               markersize=3.5, label="detection rate"),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=5.0,
              frameon=False, handlelength=1.0, handleheight=0.55,
              labelspacing=0.25)


def panel_S11b(ax, boot, hk_struct, hk_metab):
    """Bootstrap CI vs Hartung–Knapp CI for axis pooled g."""
    rows = [
        {"axis": "STRUCTURAL", "method": "HK (Fig 6A)",
         "g": -0.568, "lo": hk_struct[0], "hi": hk_struct[1]},
        {"axis": "STRUCTURAL", "method": "Bootstrap",
         "g": float(boot.loc[boot["axis"] == "STRUCTURAL", "real_g"].iloc[0]),
         "lo": float(boot.loc[boot["axis"] == "STRUCTURAL", "boot_q025"].iloc[0]),
         "hi": float(boot.loc[boot["axis"] == "STRUCTURAL", "boot_q975"].iloc[0])},
        {"axis": "METABOLIC",  "method": "HK (Fig 6A)",
         "g": -0.268, "lo": hk_metab[0], "hi": hk_metab[1]},
        {"axis": "METABOLIC",  "method": "Bootstrap",
         "g": float(boot.loc[boot["axis"] == "METABOLIC", "real_g"].iloc[0]),
         "lo": float(boot.loc[boot["axis"] == "METABOLIC", "boot_q025"].iloc[0]),
         "hi": float(boot.loc[boot["axis"] == "METABOLIC", "boot_q975"].iloc[0])},
    ]
    df = pd.DataFrame(rows)
    y = np.arange(len(df))
    for i, r in df.iterrows():
        col = AXIS_COL[r["axis"]]
        ax.plot([r["lo"], r["hi"]], [i, i],
                color=col, lw=1.2, solid_capstyle="round", zorder=3)
        ax.scatter([r["g"]], [i], s=42, color=col,
                   edgecolor="white", linewidth=0.6, zorder=4)
        ax.text(r["hi"] + 0.03, i,
                f"{r['g']:+.2f}  [{r['lo']:+.2f}, {r['hi']:+.2f}]",
                fontsize=4.8, ha="left", va="center", color="#222222",
                family="monospace")
    ax.axvline(0, color="#888888", lw=0.5, ls=(0, (1.5, 1.5)), zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['axis']}\n{r['method']}" for _, r in df.iterrows()],
                       fontsize=5.2, linespacing=1.15)
    ax.invert_yaxis()
    ax.set_xlim(-2.0, 1.4)
    ax.set_xlabel("pooled g (AF − SR)", fontsize=6.0, labelpad=2)
    ax.set_title("S11b — bootstrap CI vs Hartung–Knapp CI",
                 fontsize=6.4, pad=3, fontweight="bold", loc="left")
    ax.tick_params(axis="x", labelsize=5.4, length=2)
    ax.tick_params(axis="y", length=0, pad=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5); ax.spines["bottom"].set_linewidth(0.5)


def panel_S11c(ax, classic):
    """Classical limma per-program direction concordance with framework."""
    cohorts = ["GSE128188", "GSE79768", "GSE115574"]

    summary = []
    for cohort in cohorts:
        sub = classic[classic["cohort"] == cohort]
        n_total = len(sub)
        n_concord = int(sub["concordant"].sum())
        pct = n_concord / max(n_total, 1) * 100
        summary.append({"cohort": cohort, "n_concord": n_concord,
                        "n_total": n_total, "pct": pct})
    sumdf = pd.DataFrame(summary)

    y = np.arange(len(sumdf))
    cols = ["#1B7837" if v >= 80 else "#D0D000" if v >= 60 else "#B2182B"
            for v in sumdf["pct"]]
    ax.barh(y, sumdf["pct"], color=cols, edgecolor="white", linewidth=0.4)
    ax.axvline(80, color="#888888", lw=0.4, ls=(0, (1.5, 1.5)))
    for i, r in sumdf.iterrows():
        ax.text(r["pct"] + 1, i,
                f"{int(r['n_concord'])} / {int(r['n_total'])}  ({r['pct']:.0f}%)",
                fontsize=5.0, ha="left", va="center", color="#222222")
    ax.set_yticks(y)
    ax.set_yticklabels(sumdf["cohort"], fontsize=5.6)
    ax.invert_yaxis()
    ax.set_xlim(0, 110)
    ax.set_xlabel("framework-vs-classical direction concordance  (%)",
                  fontsize=6.0, labelpad=2)
    ax.set_title("S11c — classical limma per-chamber AF vs SR concordance",
                 fontsize=6.4, pad=3, fontweight="bold", loc="left")
    ax.tick_params(axis="x", labelsize=5.4, length=2)
    ax.tick_params(axis="y", length=0, pad=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5); ax.spines["bottom"].set_linewidth(0.5)


def main() -> None:
    set_nature_style()
    sim = pd.read_csv(RES / "method_p5_simulation.tsv", sep="\t")
    boot = pd.read_csv(RES / "method_p5_bootstrap_axis.tsv", sep="\t")
    classic = pd.read_csv(RES / "method_p5_classical_vs_framework.tsv",
                          sep="\t")

    # Hartung–Knapp values pulled from Fig 6A meta tsv
    meta = pd.read_csv(RES / "figure6_meta_axis.tsv", sep="\t")
    pool = meta[meta["type"] == "pooled"].set_index("axis")
    hk_struct = (pool.loc["STRUCTURAL", "ci_lo"], pool.loc["STRUCTURAL", "ci_hi"])
    hk_metab  = (pool.loc["METABOLIC",  "ci_lo"], pool.loc["METABOLIC",  "ci_hi"])

    fig_w = COL_DOUBLE
    fig_h = 5.6
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.07, right=0.97, top=0.92, bottom=0.10)
    gs = fig.add_gridspec(2, 2,
                          height_ratios=[1.0, 1.0],
                          width_ratios=[1.4, 1.0],
                          hspace=0.55, wspace=0.45)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    panel_S11a(ax_a, sim)
    panel_S11b(ax_b, boot, hk_struct, hk_metab)
    panel_S11c(ax_c, classic)

    fig.suptitle(
        "Method-validity P5 — simulation, bootstrap, classical-comparison",
        fontsize=7.4, y=0.965, fontweight="bold",
    )
    fig.text(0.012, 0.97, "S11",
             fontsize=10, fontweight="bold", ha="left", va="top",
             transform=fig.transFigure)
    fig.text(0.5, 0.018,
             "Synthetic data: cohort-matched donor counts (5+5, 7+6, 13+15); "
             "200 reps per effect size.  Bootstrap: 2000 cohort-donor "
             "resamples.  Classical: per-cohort, per-chamber limma AF vs SR; "
             "concordance is sign-equality with framework Welch t.",
             ha="center", va="bottom", fontsize=4.6,
             color="#666666", style="italic", linespacing=1.4)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
