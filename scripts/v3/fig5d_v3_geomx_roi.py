#!/usr/bin/env python3
"""
Fig 5d v3 - GeoMx ROI spatial corroboration  (6 v3-priority families).

A 4-panel composite (mirrors v2 layout but rebuilt with v3 LE signatures):

    a  Heatmap of all signature scores across 24 AA cardiac DNA+ ROIs
        (rows = 6 family LE signatures + 6 cell-type marker panels;
         cols = ROIs ordered by rhythm then area then patient),
        z-scored per row.

    b  Family x cell-type Spearman rho matrix (6 x 6).  Which v3 family
        signature co-varies with which cell type across the 24 ROIs?
        This is the spatial-localisation panel.

    c  AF - SR Delta per signature with 95 % CI (Welch).

    d  Border - Distal Delta per signature with 95 % CI (Welch).

Family signatures = UNION of NCP_v3 LE genes per family (g3b_v3 output).
Cell-type panels = canonical hand-curated markers (NOT data-derived).

Inputs:
    Project_P5/GSE261363/roi_scores_v3.tsv  (g3b_v3 output)

Output:
    figures/fig5d_v3_geomx_roi.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
ROI_TSV = ROOT / "Project_P5" / "GSE261363" / "roi_scores_v3.tsv"
FIG_OUT = ROOT / "figures" / "fig5d_v3_geomx_roi.pdf"

LA_COLOR = "#D6604D"
RA_COLOR = "#4393C3"
GREY_DARK = "#444444"
BRBG_BROWN = "#8C510A"
BRBG_TEAL  = "#01665E"
PRGN_PURPLE = "#762A83"
PRGN_GREEN  = "#1B7837"
AF_COLOR = "#B56F6F"
SR_COLOR = "#6B9AAD"

# Family signatures, in atlas-priority descending order.
# All 6 are LA-driven or Bilateral by ChamberDriver (no RA-driven priority families).
FAMILY_ORDER = [
    "F15  NTRK2 / TrkB",
    "F11  BMP signalling",
    "F9   HS-GAG metabolism",
    "F10  TGF-beta signalling",
    "F2   Mito Ca2+ transport",
    "F4   BCAA catabolism",
]
FAMILY_DRIVER = {
    "F15  NTRK2 / TrkB":           "LA-driven",
    "F11  BMP signalling":         "LA-driven",
    "F9   HS-GAG metabolism":      "LA-driven",
    "F10  TGF-beta signalling":    "LA-driven",
    "F2   Mito Ca2+ transport":    "Bilateral",
    "F4   BCAA catabolism":        "Bilateral",
}
FAMILY_PRIORITY = {
    "F15  NTRK2 / TrkB":           "***",
    "F11  BMP signalling":         "***",
    "F9   HS-GAG metabolism":      "***",
    "F10  TGF-beta signalling":    "***",
    "F2   Mito Ca2+ transport":    "*",
    "F4   BCAA catabolism":        "*",
}
FAMILY_SHORT = {
    "F15  NTRK2 / TrkB":           "F15  NTRK2",
    "F11  BMP signalling":         "F11  BMP",
    "F9   HS-GAG metabolism":      "F9   HS-GAG",
    "F10  TGF-beta signalling":    "F10  TGF-\u03b2",
    "F2   Mito Ca2+ transport":    "F2   Mito Ca\u00b2\u207a",
    "F4   BCAA catabolism":        "F4   BCAA",
}

CELL_ORDER = [
    "fibroblast",
    "lymphatic_endothelial",
    "vascular_endothelial",
    "atrialCM",
    "macrophage",
    "adipocyte",
]
CELL_SHORT = {
    "fibroblast":            "fibroblast",
    "lymphatic_endothelial": "lymph endo",
    "vascular_endothelial":  "vasc endo",
    "atrialCM":              "atrialCM",
    "macrophage":            "macrophage",
    "adipocyte":             "adipocyte",
}

ALL_SIG_ORDER = FAMILY_ORDER + CELL_ORDER


def set_publication_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size":         8,
        "axes.titlesize":    9,
        "axes.labelsize":    8,
        "xtick.labelsize":   6.5,
        "ytick.labelsize":   6.5,
        "legend.fontsize":   6.5,
        "axes.linewidth":    0.6,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def welch_ci(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan"), float("nan")
    diff = a.mean() - b.mean()
    se = float(np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)))
    if se == 0:
        return diff, 0.0, 1.0
    df = (a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) ** 2 / \
         ((a.var(ddof=1) ** 2) / (len(a) ** 2 * (len(a) - 1))
          + (b.var(ddof=1) ** 2) / (len(b) ** 2 * (len(b) - 1)))
    p = 2 * (1 - stats.t.cdf(abs(diff) / se, df))
    tcrit = stats.t.ppf(0.975, df)
    return diff, tcrit * se, p


def main():
    set_publication_style()
    long = pd.read_csv(ROI_TSV, sep="\t")
    long = long[long["eligible_for_5d"] == True].copy()  # noqa: E712
    long = long[long["signature_lab"].isin(ALL_SIG_ORDER)]

    W = (long.pivot_table(index="signature_lab", columns="gsm",
                              values="score").reindex(ALL_SIG_ORDER))
    meta = (long[["gsm", "patient", "rhythm", "area"]]
                .drop_duplicates("gsm")
                .set_index("gsm")
                .loc[W.columns])
    order = meta.sort_values(["rhythm", "area", "patient"]).index.tolist()
    W = W[order]
    meta = meta.loc[order]
    n_roi = len(order)

    # Per-row z-score (heatmap a)
    Wz = W.sub(W.mean(axis=1), axis=0).div(W.std(axis=1, ddof=1), axis=0)

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(210 * MM, 180 * MM))
    grid = fig.add_gridspec(
        nrows=2, ncols=2,
        height_ratios=[1.05, 1.0],
        width_ratios=[1.55, 1.0],
        hspace=0.45, wspace=0.30,
        left=0.16, right=0.96, top=0.91, bottom=0.08)

    fig.suptitle(
        "Fig 5d v3   GeoMx AA ROI spatial corroboration   "
        f"(n = {n_roi} cardiac DNA+ ROIs)",
        fontsize=8.8, fontweight="bold", y=0.97, x=0.02, ha="left")
    fig.text(0.02, 0.945,
              "6 NCP_v3 priority families (rank-singscore on LE union)  +  6 canonical cell-type marker panels.  "
              "Spatial layer is corroboration only; mechanistic causality is NOT inferred from 5d.",
              fontsize=6.0, color="#444444", style="italic")

    # ===================================================================
    # Panel a: signature x ROI heatmap
    # ===================================================================
    ax_a = fig.add_subplot(grid[0, 0])
    K, _ = Wz.shape
    vmax = float(np.nanmax(np.abs(Wz.values)))
    vmax = min(vmax, 2.5)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "BrBG", [(0.0, BRBG_TEAL), (0.5, "#FFFFFF"), (1.0, BRBG_BROWN)])
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    ax_a.imshow(Wz.values, aspect="auto", cmap=cmap, norm=norm,
                  interpolation="nearest")
    ax_a.set_yticks(range(K))
    yticklabels = [FAMILY_SHORT.get(s, CELL_SHORT.get(s, s)) for s in ALL_SIG_ORDER]
    ax_a.set_yticklabels(yticklabels, fontsize=6.5)
    for i, s in enumerate(ALL_SIG_ORDER):
        if s in FAMILY_ORDER:
            c = LA_COLOR if FAMILY_DRIVER[s] == "LA-driven" else "#7E5BBF"
        else:
            c = GREY_DARK
        ax_a.get_yticklabels()[i].set_color(c)

    # gap between families and cells
    n_fam = len(FAMILY_ORDER)
    ax_a.axhline(n_fam - 0.5, color="white", lw=2.5, zorder=4)
    ax_a.axhline(n_fam - 0.5, color=GREY_DARK, lw=0.5, zorder=4.5, ls=":")

    ax_a.set_xticks([])
    band_h = 0.15
    band_y0 = -0.6
    cur_x = 0
    for (rhy, area), grp in meta.groupby(["rhythm", "area"], sort=False):
        n = len(grp)
        rcol = AF_COLOR if rhy == "AF" else SR_COLOR
        ax_a.add_patch(mpatches.Rectangle(
            (cur_x - 0.5, band_y0 - band_h * 1.2),
            n, band_h, facecolor=rcol, alpha=0.6,
            edgecolor="white", linewidth=0.5,
            clip_on=False, zorder=4))
        ax_a.text(cur_x - 0.5 + n / 2, band_y0 - band_h * 0.6,
                    f"{rhy}  \u00b7  {area}  (n={n})",
                    ha="center", va="bottom",
                    fontsize=6.0, fontweight="bold",
                    color="black", clip_on=False, zorder=5)
        cur_x += n

    cax = fig.add_axes([0.555, 0.625, 0.010, 0.180])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
                                       orientation="vertical")
    cb.set_label("z-score (per row)", fontsize=6.5)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)
    ax_a.set_title("a    Signature \u00d7 ROI heatmap   (z-scored per row)",
                     loc="left", fontsize=8.0, fontweight="bold", pad=14)

    # ===================================================================
    # Panel b: family x cell-type Spearman rho
    # ===================================================================
    ax_b = fig.add_subplot(grid[0, 1])
    fams = FAMILY_ORDER
    cells = CELL_ORDER
    R = np.zeros((len(fams), len(cells)))
    P = np.zeros((len(fams), len(cells)))
    for i, p in enumerate(fams):
        for j, c in enumerate(cells):
            r, pv = stats.spearmanr(W.loc[p].values, W.loc[c].values)
            R[i, j] = r; P[i, j] = pv

    cmap2 = mcolors.LinearSegmentedColormap.from_list(
        "PRGn", [(0.0, PRGN_PURPLE), (0.5, "#FFFFFF"), (1.0, PRGN_GREEN)])
    rmax = max(0.5, float(np.max(np.abs(R))) + 0.05)
    norm2 = mcolors.TwoSlopeNorm(vmin=-rmax, vcenter=0, vmax=rmax)
    for i, p in enumerate(fams):
        for j, c in enumerate(cells):
            r = R[i, j]; pv = P[i, j]
            ax_b.add_patch(mpatches.Rectangle(
                (j, i), 1.0, 1.0,
                facecolor=cmap2(norm2(r)),
                edgecolor="white", linewidth=0.6))
            shade = abs(r) / rmax
            txt_color = "white" if shade > 0.65 else "black"
            ax_b.text(j + 0.5, i + 0.55, f"{r:+.2f}",
                         ha="center", va="center",
                         fontsize=6.5, color=txt_color, fontweight="bold")
            if pv < 0.10:
                mark = "*" if pv < 0.05 else "."
                badge_face = "#000000" if mark == "*" else "#FFFFFF"
                badge_text = "white" if mark == "*" else "black"
                ax_b.add_patch(mpatches.Circle(
                    (j + 0.85, i + 0.18), 0.10,
                    facecolor=badge_face, edgecolor="black",
                    linewidth=0.4, zorder=4))
                ax_b.text(j + 0.85, i + 0.20, mark,
                             ha="center", va="center",
                             fontsize=7.5 if mark == "*" else 8.5,
                             fontweight="bold", color=badge_text, zorder=5)
    ax_b.set_xlim(0, len(cells)); ax_b.set_ylim(len(fams), 0)
    ax_b.set_xticks(np.arange(len(cells)) + 0.5)
    ax_b.set_xticklabels([CELL_SHORT[c] for c in cells],
                            rotation=35, ha="right", rotation_mode="anchor",
                            fontsize=6.5)
    ax_b.set_yticks(np.arange(len(fams)) + 0.5)
    ax_b.set_yticklabels([FAMILY_SHORT[p] for p in fams], fontsize=6.5)
    for i, p in enumerate(fams):
        c = LA_COLOR if FAMILY_DRIVER[p] == "LA-driven" else "#7E5BBF"
        ax_b.get_yticklabels()[i].set_color(c)
    for s in ax_b.spines.values():
        s.set_linewidth(0.4); s.set_color("#999")
    cax2 = fig.add_axes([0.945, 0.625, 0.010, 0.180])
    cb2 = mpl.colorbar.ColorbarBase(cax2, cmap=cmap2, norm=norm2,
                                        orientation="vertical")
    cb2.set_label("Spearman \u03c1", fontsize=6.5)
    cb2.outline.set_linewidth(0.4)
    cb2.ax.tick_params(labelsize=6, length=2, width=0.4)
    ax_b.set_title("b    Family \u00d7 cell-type \u03c1 across 24 AA cardiac ROIs",
                     loc="left", fontsize=8.0, fontweight="bold", pad=4)

    # ===================================================================
    # Panel c: AF - SR Delta per signature
    # ===================================================================
    ax_c = fig.add_subplot(grid[1, 0])
    rows_c = []
    for s in ALL_SIG_ORDER:
        af = W.loc[s, meta["rhythm"] == "AF"].dropna().values
        sr = W.loc[s, meta["rhythm"] == "SR"].dropna().values
        d, ci, p = welch_ci(af, sr)
        rows_c.append({"sig": s, "delta": d, "ci": ci, "p": p})
    cdf = pd.DataFrame(rows_c).iloc[::-1].reset_index(drop=True)
    y = np.arange(len(cdf))
    cols = []
    for s in cdf["sig"]:
        if s in FAMILY_ORDER:
            cols.append("#E89C8B" if FAMILY_DRIVER[s] == "LA-driven" else "#C8B8E0")
        else:
            cols.append("#CCCCCC")
    ax_c.barh(y, cdf["delta"], xerr=cdf["ci"],
                color=cols, edgecolor="black", linewidth=0.4,
                error_kw=dict(lw=0.7), height=0.7)
    ax_c.axvline(0, color=GREY_DARK, lw=0.6)
    for yi, p in zip(y, cdf["p"]):
        if p < 0.10:
            mark = "*" if p < 0.05 else "."
            ax_c.text(ax_c.get_xlim()[1] * 0.98, yi, mark,
                         ha="right", va="center",
                         fontsize=10, fontweight="bold")
    ax_c.set_yticks(y)
    ax_c.set_yticklabels([FAMILY_SHORT.get(s, CELL_SHORT.get(s, s))
                              for s in cdf["sig"]], fontsize=6.5)
    for i, s in enumerate(cdf["sig"]):
        if s in FAMILY_ORDER:
            c = LA_COLOR if FAMILY_DRIVER[s] == "LA-driven" else "#7E5BBF"
        else:
            c = GREY_DARK
        ax_c.get_yticklabels()[i].set_color(c)
    ax_c.set_xlabel("\u0394 singscore  (AF \u2212 SR)   \u00b7   Welch  \u00b1 95 % CI",
                       fontsize=7)
    ax_c.set_title("c    Rhythm effect (AF \u2212 SR)",
                     loc="left", fontsize=8.0, fontweight="bold", pad=4)

    # ===================================================================
    # Panel d: Border - Distal Delta
    # ===================================================================
    ax_d = fig.add_subplot(grid[1, 1])
    rows_d = []
    for s in ALL_SIG_ORDER:
        bo = W.loc[s, meta["area"] == "border"].dropna().values
        di = W.loc[s, meta["area"] == "distal"].dropna().values
        d, ci, p = welch_ci(bo, di)
        rows_d.append({"sig": s, "delta": d, "ci": ci, "p": p})
    ddf = pd.DataFrame(rows_d).iloc[::-1].reset_index(drop=True)
    y = np.arange(len(ddf))
    cols = []
    for s in ddf["sig"]:
        if s in FAMILY_ORDER:
            cols.append("#E89C8B" if FAMILY_DRIVER[s] == "LA-driven" else "#C8B8E0")
        else:
            cols.append("#CCCCCC")
    ax_d.barh(y, ddf["delta"], xerr=ddf["ci"],
                color=cols, edgecolor="black", linewidth=0.4,
                error_kw=dict(lw=0.7), height=0.7)
    ax_d.axvline(0, color=GREY_DARK, lw=0.6)
    for yi, p in zip(y, ddf["p"]):
        if p < 0.10:
            mark = "*" if p < 0.05 else "."
            ax_d.text(ax_d.get_xlim()[1] * 0.98, yi, mark,
                         ha="right", va="center",
                         fontsize=10, fontweight="bold")
    ax_d.set_yticks(y)
    ax_d.set_yticklabels([FAMILY_SHORT.get(s, CELL_SHORT.get(s, s))
                              for s in ddf["sig"]], fontsize=6.5)
    for i, s in enumerate(ddf["sig"]):
        if s in FAMILY_ORDER:
            c = LA_COLOR if FAMILY_DRIVER[s] == "LA-driven" else "#7E5BBF"
        else:
            c = GREY_DARK
        ax_d.get_yticklabels()[i].set_color(c)
    ax_d.set_xlabel("\u0394 singscore  (border \u2212 distal)   \u00b7   Welch  \u00b1 95 % CI",
                       fontsize=7)
    ax_d.set_title("d    Area effect (border \u2212 distal)",
                     loc="left", fontsize=8.0, fontweight="bold", pad=4)

    fig.text(0.02, 0.025,
              "Family signatures: union of NCP_v3 leading-edge genes per family (rank-singscore).  "
              "Cell-type panels: hand-curated markers (NOT data-derived).  "
              "Markers: * p_t < 0.05;  . p_t < 0.10  (Welch).  "
              "Tick-label colour: red = LA-driven family, purple = Bilateral family, grey = cell-type.",
              fontsize=5.6, color="#444", ha="left", va="bottom")
    fig.text(0.02, 0.005,
              "Family selection rationale:  F15/F11/F9/F10 = atlas Triple-/Double-/Single-positive (***);  "
              "F2 = Strong-trend duration (Bilateral *);  F4 = Adipo-carrier metabolic (Bilateral *) for Fig 3F sub-family link.",
              fontsize=5.6, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")

    print("\n=== Top |rho| family x cell-type cross-correlations ===")
    pairs = []
    for i, p in enumerate(fams):
        for j, c in enumerate(cells):
            pairs.append((p, c, R[i, j], P[i, j]))
    for p, c, r, pv in sorted(pairs, key=lambda x: -abs(x[2]))[:10]:
        sig = "*" if pv < 0.05 else ("." if pv < 0.10 else "")
        print(f"  {p:30s}  x  {c:24s}  rho={r:+.2f}  p={pv:.3f} {sig}")


if __name__ == "__main__":
    main()
