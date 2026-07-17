#!/usr/bin/env python3
"""Figure 5D  GeoMx AA ROI corroboration of multi-family carrier patterns."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
ROI_TSV = ROOT / "Project_P5" / "GSE261363" / "roi_scores.tsv"
FIG_OUT = ROOT / "figures" / "fig5d_geomx_roi.pdf"

LA_COLOR  = "#D6604D"
RA_COLOR  = "#4393C3"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"
BRBG_BROWN = "#8C510A"
BRBG_TEAL  = "#01665E"
PRGN_PURPLE = "#762A83"
PRGN_GREEN  = "#1B7837"

AF_COLOR = "#B56F6F"
SR_COLOR = "#6B9AAD"
BORDER_COLOR = "#DBA76A"
DISTAL_COLOR = "#7A8FA6"

# Display order for signatures (programs first by family, then cell types)
SIG_ORDER = [
    # programs
    "F10  GAG metabolism",
    "F11  BMP signalling",
    "F9   HS-GAG metabolism",
    "F2   Mito Ca²⁺ transport",
    "F7   OXPHOS (LA control)",
    # cell types
    "fibroblast",
    "lymphatic_endothelial",
    "vascular_endothelial",
    "atrialCM",
    "macrophage",
]
# Program rows are RA-biased (top 3) + LA-biased (next 2)
SIG_TYPES = {
    "F10  GAG metabolism":         "RA-program",
    "F11  BMP signalling":         "RA-program",
    "F9   HS-GAG metabolism":      "RA-program",
    "F2   Mito Ca²⁺ transport":    "LA-program",
    "F7   OXPHOS (LA control)":    "LA-program",
    "fibroblast":                  "cell",
    "lymphatic_endothelial":       "cell",
    "vascular_endothelial":        "cell",
    "atrialCM":                    "cell",
    "macrophage":                  "cell",
}
# Display labels (shorter for plots)
SIG_DISPLAY = {
    "F10  GAG metabolism":         "F10  GAG metab",
    "F11  BMP signalling":         "F11  BMP",
    "F9   HS-GAG metabolism":      "F9   HS-GAG",
    "F2   Mito Ca²⁺ transport":    "F2   Mito Ca²⁺",
    "F7   OXPHOS (LA control)":    "F7   OXPHOS",
    "fibroblast":                  "fibroblast",
    "lymphatic_endothelial":       "lymph endo",
    "vascular_endothelial":        "vasc endo",
    "atrialCM":                    "atrialCM",
    "macrophage":                  "macrophage",
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
        "ytick.labelsize":   6.5,
        "legend.fontsize":   6.5,
        "axes.linewidth":    0.6,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def welch_ci(a, b):
    """Return (mean(a) − mean(b), 95% CI half-width).  Welch."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan"), float("nan")
    diff = a.mean() - b.mean()
    se = float(np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)))
    if se == 0:
        return diff, 0.0, 1.0
    # df via Welch-Satterthwaite
    df = (a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) ** 2 / \
         ((a.var(ddof=1) ** 2) / (len(a) ** 2 * (len(a) - 1))
          + (b.var(ddof=1) ** 2) / (len(b) ** 2 * (len(b) - 1)))
    tcrit = stats.t.ppf(0.975, df)
    p = 2 * (1 - stats.t.cdf(abs(diff) / se, df))
    return diff, tcrit * se, p


def main() -> None:
    set_publication_style()
    long = pd.read_csv(ROI_TSV, sep="\t")
    long = long[long["eligible_for_5D"] == True].copy()
    long = long[long["signature_lab"].isin(SIG_ORDER)]

    # wide: rows = signature, cols = gsm
    W = long.pivot_table(index="signature_lab", columns="gsm",
                            values="score").reindex(SIG_ORDER)
    meta = (long[["gsm", "patient", "rhythm", "area"]]
                .drop_duplicates("gsm")
                .set_index("gsm")
                .loc[W.columns])

    # column order: rhythm (SR before AF), then area (border before distal),
    # then patient ID; this groups paired blocks visually.
    order = meta.sort_values(["rhythm", "area", "patient"]).index.tolist()
    W = W[order]
    meta = meta.loc[order]
    n_roi = len(order)

    # Z-score per row for panel a heatmap (so signatures with different
    # native scales are comparable).
    Wz = W.sub(W.mean(axis=1), axis=0).div(W.std(axis=1, ddof=1), axis=0)

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(200 * MM, 180 * MM))
    grid = fig.add_gridspec(
        nrows=2, ncols=2,
        height_ratios=[1.05, 1.0],
        width_ratios=[1.55, 1.0],
        hspace=0.45, wspace=0.30,
        left=0.16, right=0.96, top=0.92, bottom=0.07)

    fig.suptitle(
        "Fig 5D  GeoMx AA ROI corroboration of multi-family carrier patterns "
        "(n = 24 cardiac DNA+ ROIs)",
        fontsize=8.8, fontweight="bold", y=0.97, x=0.02, ha="left")

    # Panel a: signature × ROI heatmap (z-scored per row)
    ax_a = fig.add_subplot(grid[0, 0])
    K, C = Wz.shape
    vmax = float(np.nanmax(np.abs(Wz.values)))
    vmax = min(vmax, 2.5)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "BrBG", [(0.0, BRBG_TEAL), (0.5, "#FFFFFF"), (1.0, BRBG_BROWN)])
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    ax_a.imshow(Wz.values, aspect="auto",
                  cmap=cmap, norm=norm,
                  interpolation="nearest")
    ax_a.set_yticks(range(K))
    ax_a.set_yticklabels([SIG_DISPLAY[s] for s in SIG_ORDER],
                            fontsize=6.5)
    # tick colour by type
    for i, s in enumerate(SIG_ORDER):
        t = SIG_TYPES[s]
        if t == "RA-program":   c = RA_COLOR
        elif t == "LA-program": c = LA_COLOR
        else:                   c = GREY_DARK
        ax_a.get_yticklabels()[i].set_color(c)

    # gap between programs and cells (row 4 / 5 boundary)
    ax_a.axhline(4.5, color="white", lw=2.5, zorder=4)
    ax_a.axhline(4.5, color=GREY_DARK, lw=0.5, zorder=4.5, ls=":")

    # x-tick: hide individual ROIs but draw rhythm × area band markers
    ax_a.set_xticks([])
    # top band: rhythm + area
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
                    f"{rhy}  ·  {area}  (n={n})",
                    ha="center", va="bottom",
                    fontsize=6.0, fontweight="bold",
                    color="black", clip_on=False, zorder=5)
        cur_x += n

    # colorbar inset (right of heatmap)
    cax = fig.add_axes([0.555, 0.665, 0.010, 0.180])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
                                       orientation="vertical")
    cb.set_label("z-score (per row)", fontsize=6.5)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)

    ax_a.set_title(
        "a    Signature × ROI heatmap   (z-scored per row)",
        loc="left", fontsize=8.0, fontweight="bold", pad=14)

    # Panel b: program × cell-type Spearman ρ matrix  (the core finding)
    ax_b = fig.add_subplot(grid[0, 1])
    progs = [s for s in SIG_ORDER if SIG_TYPES[s].endswith("program")]
    cells = [s for s in SIG_ORDER if SIG_TYPES[s] == "cell"]

    R = np.zeros((len(progs), len(cells)))
    P = np.zeros((len(progs), len(cells)))
    for i, p in enumerate(progs):
        for j, c in enumerate(cells):
            r, pv = stats.spearmanr(W.loc[p].values, W.loc[c].values)
            R[i, j] = r; P[i, j] = pv

    cmap2 = mcolors.LinearSegmentedColormap.from_list(
        "PRGn", [(0.0, PRGN_PURPLE), (0.5, "#FFFFFF"), (1.0, PRGN_GREEN)])
    rmax = max(0.5, float(np.max(np.abs(R))) + 0.05)
    norm2 = mcolors.TwoSlopeNorm(vmin=-rmax, vcenter=0, vmax=rmax)

    for i, p in enumerate(progs):
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
                         fontsize=6.5, color=txt_color,
                         fontweight="bold")
            if pv < 0.10:
                mark = "*" if pv < 0.05 else "."
                badge_face = "#000000" if mark == "*" else "#FFFFFF"
                badge_text = "white"  if mark == "*" else "black"
                ax_b.add_patch(mpatches.Circle(
                    (j + 0.85, i + 0.18), 0.10,
                    facecolor=badge_face, edgecolor="black",
                    linewidth=0.4, zorder=4))
                ax_b.text(j + 0.85, i + 0.20, mark,
                             ha="center", va="center",
                             fontsize=7.5 if mark == "*" else 8.5,
                             fontweight="bold", color=badge_text,
                             zorder=5)
    ax_b.set_xlim(0, len(cells)); ax_b.set_ylim(len(progs), 0)
    ax_b.set_xticks(np.arange(len(cells)) + 0.5)
    ax_b.set_xticklabels([SIG_DISPLAY[c] for c in cells],
                            rotation=35, ha="right",
                            rotation_mode="anchor",
                            fontsize=6.5)
    ax_b.set_yticks(np.arange(len(progs)) + 0.5)
    ax_b.set_yticklabels([SIG_DISPLAY[p] for p in progs],
                            fontsize=6.5)
    for i, p in enumerate(progs):
        c = RA_COLOR if SIG_TYPES[p] == "RA-program" else LA_COLOR
        ax_b.get_yticklabels()[i].set_color(c)
    for s in ax_b.spines.values():
        s.set_linewidth(0.4); s.set_color("#999")

    cax2 = fig.add_axes([0.945, 0.665, 0.010, 0.180])
    cb2 = mpl.colorbar.ColorbarBase(cax2, cmap=cmap2, norm=norm2,
                                        orientation="vertical")
    cb2.set_label("Spearman ρ", fontsize=6.5)
    cb2.outline.set_linewidth(0.4)
    cb2.ax.tick_params(labelsize=6, length=2, width=0.4)

    ax_b.set_title(
        "b    Program × cell-type ρ across 24 AA cardiac ROIs",
        loc="left", fontsize=8.0, fontweight="bold", pad=4)

    # Panel c: AF − SR Δ per signature (Welch t)
    ax_c = fig.add_subplot(grid[1, 0])
    rows_c = []
    for s in SIG_ORDER:
        af = W.loc[s, meta["rhythm"] == "AF"].dropna().values
        sr = W.loc[s, meta["rhythm"] == "SR"].dropna().values
        d, ci, p = welch_ci(af, sr)
        rows_c.append({"sig": s, "delta": d, "ci": ci, "p": p})
    cdf = pd.DataFrame(rows_c)
    y = np.arange(len(cdf))
    cdf = cdf.iloc[::-1].reset_index(drop=True); y = np.arange(len(cdf))
    cols = []
    for s in cdf["sig"]:
        if SIG_TYPES[s] == "RA-program":   cols.append("#E89C8B")
        elif SIG_TYPES[s] == "LA-program": cols.append("#9EC4E0")
        else:                              cols.append("#CCCCCC")
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
    ax_c.set_yticklabels([SIG_DISPLAY[s] for s in cdf["sig"]],
                            fontsize=6.5)
    for i, s in enumerate(cdf["sig"]):
        t = SIG_TYPES[s]
        c = RA_COLOR if t == "RA-program" else LA_COLOR if t == "LA-program" else GREY_DARK
        ax_c.get_yticklabels()[i].set_color(c)
    ax_c.set_xlabel("Δ singscore  (AF − SR)   ·   Welch  ±95 % CI", fontsize=7)
    ax_c.set_title("c    Rhythm effect (AF − SR)",
                     loc="left", fontsize=8.0, fontweight="bold", pad=4)

    # Panel d: Border − Distal Δ per signature
    ax_d = fig.add_subplot(grid[1, 1])
    rows_d = []
    for s in SIG_ORDER:
        bo = W.loc[s, meta["area"] == "border"].dropna().values
        di = W.loc[s, meta["area"] == "distal"].dropna().values
        d, ci, p = welch_ci(bo, di)
        rows_d.append({"sig": s, "delta": d, "ci": ci, "p": p})
    ddf = pd.DataFrame(rows_d)
    ddf = ddf.iloc[::-1].reset_index(drop=True); y = np.arange(len(ddf))
    cols = []
    for s in ddf["sig"]:
        if SIG_TYPES[s] == "RA-program":   cols.append("#E89C8B")
        elif SIG_TYPES[s] == "LA-program": cols.append("#9EC4E0")
        else:                              cols.append("#CCCCCC")
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
    ax_d.set_yticklabels([SIG_DISPLAY[s] for s in ddf["sig"]],
                            fontsize=6.5)
    for i, s in enumerate(ddf["sig"]):
        t = SIG_TYPES[s]
        c = RA_COLOR if t == "RA-program" else LA_COLOR if t == "LA-program" else GREY_DARK
        ax_d.get_yticklabels()[i].set_color(c)
    ax_d.set_xlabel("Δ singscore  (border − distal)   ·   Welch  ±95 % CI",
                       fontsize=7)
    ax_d.set_title("d    Area effect (border − distal)",
                     loc="left", fontsize=8.0, fontweight="bold", pad=4)

    fig.text(0.02, 0.020,
              "Signatures: 5 NCP_v2 programs (LE genes) + 5 canonical "
              "cell-type marker panels.   "
              "Markers: * p_t < 0.05;  . p_t < 0.10  (Welch).   "
              "Colour key (left tick labels): RA-biased program (red), "
              "LA-biased program (blue), cell-type panel (grey).",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")

    # console summary of strongest cross-correlations
    print("\n=== Top |ρ| program × cell-type cross-correlations ===")
    pairs = []
    for i, p in enumerate(progs):
        for j, c in enumerate(cells):
            pairs.append((p, c, R[i, j], P[i, j]))
    for p, c, r, pv in sorted(pairs, key=lambda x: -abs(x[2]))[:10]:
        sig = "★" if pv < 0.05 else ("." if pv < 0.10 else "")
        print(f"  {p:32s}  ×  {c:24s}  ρ={r:+.2f}  p={pv:.3f} {sig}")


if __name__ == "__main__":
    main()
