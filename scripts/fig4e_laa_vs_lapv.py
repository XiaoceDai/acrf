#!/usr/bin/env python3
"""
Figure 4E  LAA primary vs LAPV sensitivity ρ comparison.

Per program (n = 30) and per family (n = 11) we plot
    x = ρ in LAA primary
    y = ρ in LAPV sensitivity
on a 1:1-anchored scatter, so concordant trends fall on the y = x line
and LAA-specific trends fall on the x-axis.  Diagonal shaded band marks
the |ρ_LAA − ρ_LAPV| ≤ 0.2 "concordant" region.  Markers shaped /
coloured by baseline LA-vs-RA direction; size encodes |ρ_LAA|.

Output:  figures/fig4e_laa_vs_lapv.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from _pathway_labels import short_pathway_label

D4_DIR  = ROOT / "Project_P4" / "GSE41177" / "results_primary"
FIG_OUT = ROOT / "figures" / "fig4e_laa_vs_lapv.pdf"

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"


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
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
        "axes.spines.top":   False,
        "axes.spines.right": False,
    })


def make_axes(ax, x, y, *, marker, color, size, edgecolor="black"):
    ax.scatter(x, y, marker=marker, s=size, color=color,
                  edgecolor=edgecolor, linewidth=0.5,
                  alpha=0.85, zorder=4)


def main() -> None:
    set_publication_style()

    laa  = pd.read_csv(D4_DIR / "duration_trend_laa_primary.tsv", sep="\t")
    lapv = pd.read_csv(D4_DIR / "duration_trend_lapv_sensitivity.tsv",
                          sep="\t")
    fam_laa  = pd.read_csv(D4_DIR / "duration_trend_FAMILY_laa_primary.tsv",
                              sep="\t")
    fam_lapv = pd.read_csv(D4_DIR / "duration_trend_FAMILY_lapv_sensitivity.tsv",
                              sep="\t")

    # program-level join
    prog = laa[["program_id", "program_label", "program_order",
                  "family_id", "baseline_direction", "statistic",
                  "p_perm"]].rename(
            columns={"statistic": "rho_laa", "p_perm": "p_laa"})
    prog = prog.merge(
        lapv[["program_id", "statistic", "p_perm"]].rename(
            columns={"statistic": "rho_lapv", "p_perm": "p_lapv"}),
        on="program_id", how="inner")
    prog["short"] = prog["program_id"].apply(
        lambda s: short_pathway_label(s, 18))

    fam = fam_laa[["family_id", "direction_label", "n_programs",
                       "statistic", "p_perm"]].rename(
            columns={"statistic": "rho_laa", "p_perm": "p_laa"})
    fam = fam.merge(
        fam_lapv[["family_id", "statistic", "p_perm"]].rename(
            columns={"statistic": "rho_lapv", "p_perm": "p_lapv"}),
        on="family_id", how="inner")

    # ---- figure : 2 panels (program + family) ------------------------
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(180 * MM, 100 * MM))
    grid = fig.add_gridspec(
        nrows=1, ncols=2, wspace=0.35,
        left=0.07, right=0.97, top=0.86, bottom=0.16)

    fig.suptitle(
        "Fig 4E  LAA primary  vs  LAPV sensitivity  AF-duration ρ "
        "(n = 12 each)",
        fontsize=8.8, fontweight="bold", y=0.97, x=0.02, ha="left")

    panels = [
        ("program-level   (n = 30 programs)", prog,
         "rho_laa", "rho_lapv", "p_laa", "p_lapv", "baseline_direction",
         True),
        ("family-level   (n = 11 families)", fam,
         "rho_laa", "rho_lapv", "p_laa", "p_lapv", "direction_label",
         False),
    ]

    for idx, (title, df, xc, yc, pxc, pyc, dir_col, is_prog) in enumerate(panels):
        ax = fig.add_subplot(grid[idx])
        # axis range
        all_vals = np.concatenate([df[xc].dropna().to_numpy(),
                                       df[yc].dropna().to_numpy()])
        lim = max(0.6, float(np.max(np.abs(all_vals))) + 0.05)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

        # Concordant band: |x − y| ≤ 0.2
        xs = np.linspace(-lim, lim, 50)
        ax.fill_between(xs, xs - 0.2, xs + 0.2,
                          color=GREY_LIGHT, alpha=0.5, zorder=1,
                          label="concordant band ±0.2")
        # 1:1 line
        ax.plot([-lim, lim], [-lim, lim], color=GREY_MID, lw=0.6,
                  ls="--", zorder=2)
        # axes through 0
        ax.axvline(0, color=GREY_DARK, lw=0.5, zorder=1.5)
        ax.axhline(0, color=GREY_DARK, lw=0.5, zorder=1.5)

        # Quadrant labels
        ax.text(lim * 0.96, lim * 0.96, "both ↑",
                  ha="right", va="top", fontsize=6.0, color=GREY_DARK,
                  alpha=0.5)
        ax.text(-lim * 0.96, -lim * 0.96, "both ↓",
                  ha="left", va="bottom", fontsize=6.0, color=GREY_DARK,
                  alpha=0.5)
        ax.text(lim * 0.96, -lim * 0.96, "LAA ↑  LAPV ↓",
                  ha="right", va="bottom", fontsize=6.0,
                  color=GREY_DARK, alpha=0.5)
        ax.text(-lim * 0.96, lim * 0.96, "LAA ↓  LAPV ↑",
                  ha="left", va="top", fontsize=6.0, color=GREY_DARK,
                  alpha=0.5)

        # Plot points
        for _, r in df.iterrows():
            x, y = r[xc], r[yc]
            if pd.isna(x) or pd.isna(y):
                continue
            if is_prog:
                d = int(r[dir_col])
                col = LA_COLOR if d > 0 else RA_COLOR
                marker = "o" if d > 0 else "s"
            else:
                d = r[dir_col]
                col = LA_COLOR if d == "LA-biased" else RA_COLOR
                marker = "o" if d == "LA-biased" else "s"
            size = 30 + 80 * abs(x)
            edge = "black" if r[pxc] < 0.05 else "#777"
            lw   = 1.0 if r[pxc] < 0.05 else 0.4
            ax.scatter(x, y, marker=marker, s=size,
                          color=col, edgecolor=edge, linewidth=lw,
                          alpha=0.85, zorder=4)
            # Label significant LAA hits at family level always; at
            # program level only the strongest few
            if not is_prog and r[pxc] < 0.10:
                lbl = f"F{int(r['family_id'])}"
                ax.annotate(lbl, (x, y),
                              xytext=(6, 4), textcoords="offset points",
                              fontsize=6.5, fontweight="bold",
                              color=col,
                              zorder=5)
            if is_prog and r[pxc] < 0.10:
                ax.annotate(r["short"], (x, y),
                              xytext=(6, 4), textcoords="offset points",
                              fontsize=6.0, fontweight="bold",
                              color=col,
                              zorder=5)

        ax.set_xlabel("ρ  LAA primary", fontsize=7.5)
        ax.set_ylabel("ρ  LAPV sensitivity", fontsize=7.5)
        ax.set_title(title, loc="left", fontsize=8.0, fontweight="bold",
                       color=GREY_DARK, pad=4)
        ax.set_aspect("equal")

    # legend (right of last panel, in figure coords)
    handles = [
        mlines.Line2D([], [], marker="o", color="white",
                          markerfacecolor=LA_COLOR, markeredgecolor="black",
                          markersize=6, linestyle="",
                          label="LA-biased"),
        mlines.Line2D([], [], marker="s", color="white",
                          markerfacecolor=RA_COLOR, markeredgecolor="black",
                          markersize=6, linestyle="",
                          label="RA-biased"),
        mlines.Line2D([], [], marker="o", color="white",
                          markerfacecolor="#888", markeredgecolor="black",
                          markersize=8, linestyle="",
                          label="bigger = larger |ρ_LAA|"),
        mpatches.Patch(facecolor=GREY_LIGHT, alpha=0.5,
                          label="concordant band ±0.2"),
        mlines.Line2D([], [], color="black", lw=1.0,
                          label="thick edge: p_LAA < 0.05"),
    ]
    fig.legend(handles=handles, loc="lower right",
                  bbox_to_anchor=(0.97, 0.02), ncols=5, frameon=False,
                  fontsize=6.2)

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
