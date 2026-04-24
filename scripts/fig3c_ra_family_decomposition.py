#!/usr/bin/env python3
"""
Figure 3C  RA-biased family AF rewiring decomposition by cell type.

Three rows, one per RA-biased axis family (F9 HS-GAG, F10 TGF-β,
F11 BMP).  In each row, an 11-bar plot of Δ singscore (AF − Ctrl) per
LA-arm cell type, with 95 % CI bars derived from the per-donor SE.

This is the "where in LA tissue does the RA-biased family signal land?"
panel that links Fig 2D (the family-level forest of bulk aCAI rewiring)
to single-cell carriers in GSE255612.

Inputs:
    results/fig3_primary_family_celltype_delta.tsv

Output:
    figures/fig3c_ra_family_decomposition.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

DELTA_TSV = ROOT / "results" / "fig3_primary_family_celltype_delta.tsv"
FIG_OUT   = ROOT / "figures" / "fig3c_ra_family_decomposition.pdf"

LA_COLOR  = "#D6604D"
RA_COLOR  = "#4393C3"
BRBG_BROWN = "#8C510A"
BRBG_TEAL  = "#01665E"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"

CT_LABELS = {
    "atrialCM":                "atrialCM",
    "fibroblast":              "fibroblast",
    "endothelial":             "endothelial",
    "lymphatic_endothelial":   "lymph endo",
    "pericyte":                "pericyte",
    "vascular_smooth_muscle":  "VSMC",
    "adipocyte":               "adipocyte",
    "macrophage":              "macrophage",
    "lymphocyte":              "lymphocyte",
    "mast_cell":               "mast cell",
    "neuronal":                "neuronal",
}
CT_ORDER = list(CT_LABELS.keys())

FAMILY_LABELS = {
    9:  "F9   HS-GAG metabolism family",
    10: "F10  TGF-β signalling family",
    11: "F11  BMP signalling family",
}


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "axes.titlesize":    9,
        "axes.labelsize":    8,
        "xtick.labelsize":   6.5,
        "ytick.labelsize":   7,
        "axes.linewidth":    0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def main() -> None:
    set_publication_style()
    df = pd.read_csv(DELTA_TSV, sep="\t")

    # SE of dScore = sqrt(sd_af^2 / n_af + sd_ctrl^2 / n_ctrl)
    df["se_dscore"] = np.sqrt(
        (df["sd_score_af"] ** 2) / df["n_af_donors"]
        + (df["sd_score_ctrl"] ** 2) / df["n_ctrl_donors"])
    df["ci_lo"] = df["dScore"] - 1.96 * df["se_dscore"]
    df["ci_hi"] = df["dScore"] + 1.96 * df["se_dscore"]

    ra = df[df["family_id"].isin([9, 10, 11])].copy()
    ra = ra[ra["cell_type_included"]]

    # ---- common y-range across the 3 family panels (so bars are comparable) -
    y_min = float(np.nanmin(ra["ci_lo"]))
    y_max = float(np.nanmax(ra["ci_hi"]))
    pad = max(0.005, (y_max - y_min) * 0.10)
    y_lo, y_hi = y_min - pad, y_max + pad

    # ---- figure --------------------------------------------------------
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(180 * MM, 145 * MM))
    grid = fig.add_gridspec(
        nrows=3, ncols=1, hspace=0.35,
        left=0.10, right=0.97, top=0.91, bottom=0.10)

    fig.suptitle(
        "Fig 3C  RA-biased family AF rewiring decomposition across "
        "11 LA snRNA cell types  (chamber-reversal carriers)",
        fontsize=8.8, fontweight="bold", y=0.985, x=0.02, ha="left")

    fam_ids = [9, 10, 11]
    for idx, fid in enumerate(fam_ids):
        ax = fig.add_subplot(grid[idx])
        sub = ra[ra["family_id"] == fid].set_index("cell_type")
        sub = sub.reindex(CT_ORDER)

        x = np.arange(len(CT_ORDER))
        d = sub["dScore"].to_numpy()
        ci_lo = sub["ci_lo"].to_numpy()
        ci_hi = sub["ci_hi"].to_numpy()
        p = sub["p_perm"].to_numpy()

        # Bar colours: brown for + (AF > Ctrl), teal for −, grey for NaN
        colors = []
        for v in d:
            if pd.isna(v):
                colors.append("#DDDDDD")
            elif v >= 0:
                colors.append(BRBG_BROWN)
            else:
                colors.append(BRBG_TEAL)

        # zero reference line
        ax.axhline(0, color=GREY_DARK, lw=0.6, zorder=1)

        # Bars
        ax.bar(x, d, width=0.62, color=colors, edgecolor="black",
                 linewidth=0.4, zorder=3)
        # Error bars
        for xi, lo, hi, val in zip(x, ci_lo, ci_hi, d):
            if pd.isna(val):
                continue
            ax.plot([xi, xi], [lo, hi],
                     color="black", lw=0.7, zorder=4)
            ax.plot([xi - 0.10, xi + 0.10], [lo, lo],
                     color="black", lw=0.7, zorder=4)
            ax.plot([xi - 0.10, xi + 0.10], [hi, hi],
                     color="black", lw=0.7, zorder=4)

        # Significance markers above each bar
        for xi, val, hi, pv in zip(x, d, ci_hi, p):
            if pd.isna(val) or pd.isna(pv) or pv >= 0.10:
                continue
            mark = "*" if pv < 0.05 else "."
            y_anchor = hi + (y_hi - y_lo) * 0.025
            if val < 0:
                # for negative bars, place marker above zero / at top of CI
                y_anchor = max(hi, 0) + (y_hi - y_lo) * 0.025
            ax.text(xi, y_anchor, mark,
                     ha="center", va="bottom",
                     fontsize=11 if mark == "*" else 13,
                     fontweight="bold", color="black",
                     zorder=5)

        # Number of programs in family
        n_progs = int(sub["n_programs"].dropna().iloc[0]) if not sub["n_programs"].dropna().empty else 0

        ax.set_xticks(x)
        ax.set_xticklabels([CT_LABELS[c] for c in CT_ORDER],
                              rotation=30, ha="right",
                              rotation_mode="anchor")
        ax.set_ylim(y_lo, y_hi)
        ax.set_ylabel("Δ singscore  (AF − Ctrl)", fontsize=7)
        ax.tick_params(axis="x", pad=2)
        ax.set_title(
            f"{FAMILY_LABELS[fid]}   ·   "
            f"{n_progs} program{'s' if n_progs != 1 else ''} in family   "
            f"·   pooled bulk aCAI g (Fig 2D) = "
            f"{ {9:'+0.79', 10:'+0.68', 11:'+0.85'}[fid] }",
            loc="left", fontsize=7.5, fontweight="bold",
            color=RA_COLOR, pad=4)

        # Box around significant cells for emphasis
        for xi, val, pv in zip(x, d, p):
            if pd.isna(val) or pd.isna(pv) or pv >= 0.05:
                continue
            ax.add_patch(mpatches.Rectangle(
                (xi - 0.34, min(val, 0) - 0.001),
                0.68, abs(val) + 0.002,
                facecolor="none", edgecolor=RA_COLOR,
                linewidth=1.2, linestyle="--", zorder=6))

    # ---- caption ----
    fig.text(0.02, 0.020,
              "Bars: family-mean Δ singscore (AF − Ctrl).   "
              "Whiskers: 95 % CI from per-donor SE.   "
              "Markers: * p_perm < 0.05,  . p_perm < 0.10  (10,000 MC perms).   "
              "Dashed red box: cells with p_perm < 0.05.   "
              "Brown = AF higher  /  teal = AF lower  (BrBG palette).",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")

    # console summary of the key findings for the manuscript text
    print("\n=== RA-family × cell-type p<0.10 hits ===")
    for fid in fam_ids:
        sub = ra[(ra["family_id"] == fid) & (ra["p_perm"] < 0.10)] \
                  .sort_values("p_perm")
        print(f"\n  F{fid} ({FAMILY_LABELS[fid]}):")
        if sub.empty:
            print("     (none with p<0.10)")
        for _, r in sub.iterrows():
            print(f"     {r['cell_type']:25s}  Δ = {r['dScore']:+.4f}  "
                    f"95% CI [{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}]  "
                    f"p = {r['p_perm']:.3f}")


if __name__ == "__main__":
    main()
