#!/usr/bin/env python3
"""
Fig 3b v3 - family x compartment Delta-singscore heatmap (LA snRNA).

Layout (left to right):
    LEFT  carrier strip      : HCA baseline LA carrier (Fig 1F)
    LEFT  family / direction : F-id + bias color
    LEFT  family label       : short pathway name
    CENTER heatmap           : 16 family x 9 compartment dScore
    BOTTOM histogram         : distribution of all 144 dScore values

Color scale locked:  vmin = -0.04, vmax = +0.04   (under/over -> deep)
Significance markers: *  p_perm < 0.05 ;  .  p_perm < 0.10
Row order: atlas PriorityScore descending (cross-figure consistency).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label  # noqa: E402

DELTA_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
ATLAS_TSV = ROOT / "results" / "fig_integrated_af_rewiring_atlas.tsv"
FAM_MAP_TSV = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
OUT_PDF = ROOT / "figures" / "fig3_v3_b_family_compartment_heatmap.pdf"

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"

COMPARTMENTS = [
    "atrialCM",
    "fibroblast",
    "endothelial",
    "macrophage",
    "lymphocyte",
    "pericyte",
    "vascular_smooth_muscle",
    "adipocyte",
    "neuronal",
]
COMPARTMENT_PALETTE = {
    "atrialCM": "#4E79A7",
    "fibroblast": "#A0CBE8",
    "endothelial": "#F28E2B",
    "macrophage": "#59A14F",
    "lymphocyte": "#8CD17D",
    "pericyte": "#499894",
    "vascular_smooth_muscle": "#86BCB6",
    "adipocyte": "#E15759",
    "neuronal": "#79706E",
}
COMPARTMENT_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC",
    "adipocyte": "Adipo",
    "neuronal": "Neuro",
}
MIXED_COLOR = "#BDBDBD"
MIXED_LABEL = "Mixed/Unresolved"

VMIN, VMAX = -0.04, 0.04


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 8.5,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def main():
    set_style()
    delta = pd.read_csv(DELTA_TSV, sep="\t")
    hca = pd.read_csv(HCA_TSV, sep="\t")
    atlas = pd.read_csv(ATLAS_TSV, sep="\t")

    # Use atlas-derived family short labels (the same labels shown in
    # Fig 2C / atlas), so cross-figure naming is identical.
    atlas = atlas.sort_values("PriorityScore", ascending=False).reset_index(drop=True)
    fam_order = atlas["family_id"].astype(int).tolist()

    carriers = dict(zip(hca["family_id"].astype(int), hca["LA_carrier"].astype(str)))
    short_lookup = dict(zip(atlas["family_id"].astype(int), atlas["short"].astype(str)))
    dir_lookup = dict(zip(atlas["family_id"].astype(int), atlas["HealthyBias"].astype(str)))

    delta["family_id"] = delta["family_id"].astype(int)
    delta_pivot = delta.pivot(index="family_id", columns="compartment", values="dScore")
    delta_pivot = delta_pivot.reindex(index=fam_order, columns=COMPARTMENTS)
    p_pivot = delta.pivot(index="family_id", columns="compartment", values="p_perm")
    p_pivot = p_pivot.reindex(index=fam_order, columns=COMPARTMENTS)

    n_rows = len(fam_order)
    n_cols = len(COMPARTMENTS)

    cmap = LinearSegmentedColormap.from_list(
        "AF_d", [(0.0, RA_DEEP), (0.45, "#FFFFFF"), (0.55, "#FFFFFF"), (1.0, LA_DEEP)]
    )
    cmap.set_under("#1F3650")  # very deep blue for underflow
    cmap.set_over("#5C2D2A")   # very deep red for overflow
    cmap.set_bad("#E0E0E0")
    norm = TwoSlopeNorm(vmin=VMIN, vcenter=0.0, vmax=VMAX)

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 200 * MM))
    gs = fig.add_gridspec(
        2, 6,
        width_ratios=[0.30, 0.42, 3.6, 6.4, 0.50, 1.30],
        height_ratios=[6.5, 1.5],
        left=0.04, right=0.985, top=0.93, bottom=0.06,
        wspace=0.04, hspace=0.10,
    )
    ax_carrier = fig.add_subplot(gs[0, 0])
    ax_dir = fig.add_subplot(gs[0, 1])
    ax_lbl = fig.add_subplot(gs[0, 2])
    ax_main = fig.add_subplot(gs[0, 3])
    ax_sigtally = fig.add_subplot(gs[0, 4])
    ax_leg = fig.add_subplot(gs[0, 5])
    ax_hist = fig.add_subplot(gs[1, 3])
    ax_hist_leg = fig.add_subplot(gs[1, 5])

    # ---- carrier strip (left) ----
    for i, fid in enumerate(fam_order):
        carrier = carriers.get(fid, MIXED_LABEL)
        if carrier == MIXED_LABEL:
            ax_carrier.add_patch(mpatches.Rectangle(
                (0, i), 1, 1,
                facecolor=MIXED_COLOR, edgecolor="white", linewidth=0.4, hatch="///"
            ))
        else:
            ax_carrier.add_patch(mpatches.Rectangle(
                (0, i), 1, 1,
                facecolor=COMPARTMENT_PALETTE.get(carrier, MIXED_COLOR),
                edgecolor="white", linewidth=0.4
            ))
    ax_carrier.set_xlim(0, 1)
    ax_carrier.set_ylim(n_rows, 0)
    ax_carrier.set_xticks([0.5])
    ax_carrier.set_xticklabels(["LA\nbaseline\ncarrier"], fontsize=5.4, color="#444444")
    ax_carrier.set_yticks([])
    ax_carrier.tick_params(axis="x", length=0, pad=2)
    for s in ax_carrier.spines.values():
        s.set_visible(False)

    # ---- direction strip ----
    for i, fid in enumerate(fam_order):
        d = dir_lookup.get(fid, "?")
        col = LA_DEEP if d == "LA-biased" else RA_DEEP
        ax_dir.add_patch(mpatches.Rectangle((0, i), 1, 1, facecolor=col, edgecolor="white", linewidth=0.4))
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(n_rows, 0); ax_dir.set_xticks([]); ax_dir.set_yticks([])
    for s in ax_dir.spines.values():
        s.set_visible(False)

    # ---- family labels ----
    ax_lbl.set_xlim(0, 1); ax_lbl.set_ylim(n_rows, 0); ax_lbl.axis("off")
    for i, fid in enumerate(fam_order):
        d = dir_lookup.get(fid, "?")
        col = LA_DEEP if d == "LA-biased" else RA_DEEP
        ax_lbl.text(
            0.99, i + 0.5,
            f"F{int(fid)}  {short_lookup.get(fid, '')}",
            fontsize=6.5, color=col, fontweight="bold",
            ha="right", va="center"
        )

    # ---- main heatmap ----
    mat = delta_pivot.to_numpy(dtype=float)
    pmat = p_pivot.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(mat)
    im = ax_main.imshow(
        masked, cmap=cmap, norm=norm, aspect="auto",
        extent=[-0.5, n_cols - 0.5, n_rows - 0.5, -0.5],
        interpolation="nearest",
    )
    # significance markers - use a yellow-on-dark high-contrast asterisk that
    # remains legible against deep blue / red cells; box outline in same hue.
    SIG_COLOR = "#FFE066"   # warm yellow, high contrast on both red and blue
    for i in range(n_rows):
        for j in range(n_cols):
            v = mat[i, j]
            p = pmat[i, j]
            if not np.isfinite(p):
                continue
            if p < 0.05:
                ax_main.add_patch(mpatches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    facecolor="none", edgecolor=SIG_COLOR, linewidth=1.6, zorder=4,
                ))
                ax_main.text(j, i, "*", ha="center", va="center", fontsize=12,
                             color=SIG_COLOR, fontweight="bold", zorder=5,
                             path_effects=[
                                 mpe.Stroke(linewidth=1.6, foreground="#1A1A1A"),
                                 mpe.Normal(),
                             ])
            elif p < 0.10:
                ax_main.text(j, i, ".", ha="center", va="center", fontsize=11,
                             color="#1A1A1A", fontweight="bold", zorder=5)
    ax_main.set_xlim(-0.5, n_cols - 0.5)
    ax_main.set_ylim(n_rows - 0.5, -0.5)
    ax_main.set_xticks(np.arange(n_cols))
    ax_main.set_xticklabels([COMPARTMENT_SHORT[c] for c in COMPARTMENTS], rotation=30, ha="right", fontsize=7.0)
    ax_main.set_yticks([])
    ax_main.set_xlabel("cell compartment (LA snRNA)", fontsize=7.5)
    for s in ax_main.spines.values():
        s.set_visible(False)

    # ---- right: sig tally per family (count of sig cells in this row) ----
    sig_per_row = (pmat < 0.05).sum(axis=1)
    ax_sigtally.set_xlim(0, max(sig_per_row.max() + 0.5, 2.5))
    ax_sigtally.set_ylim(n_rows - 0.5, -0.5)
    for i, n_sig in enumerate(sig_per_row):
        if n_sig > 0:
            ax_sigtally.barh(i, n_sig, height=0.55, color="#1A1A1A", edgecolor="none")
            ax_sigtally.text(n_sig + 0.10, i, str(int(n_sig)),
                             fontsize=6.0, ha="left", va="center", color="#1A1A1A", fontweight="bold")
    ax_sigtally.set_xticks([0, 1, 2])
    ax_sigtally.set_xticklabels(["0", "1", "2"], fontsize=5.5, color="#444444")
    ax_sigtally.set_yticks([])
    ax_sigtally.set_xlabel("# sig\nct/family", fontsize=5.5, color="#444444")
    for s in ("top", "right", "left"):
        ax_sigtally.spines[s].set_visible(False)
    ax_sigtally.spines["bottom"].set_linewidth(0.4)

    # ---- legend (right) ----
    ax_leg.set_xlim(0, 1); ax_leg.set_ylim(1, 0); ax_leg.axis("off")
    ax_leg.text(0.0, 0.04, "Color scale", fontsize=6.6, fontweight="bold", color="#1A1A1A")
    # No values exceed +/- 0.04 in this panel, so colorbar uses no extends.
    cax = fig.add_axes([0.91, 0.78, 0.018, 0.10])
    cbar = fig.colorbar(im, cax=cax, orientation="vertical")
    cbar.set_label("\u0394singscore  (AF \u2212 Ctrl)", fontsize=5.8)
    cbar.ax.tick_params(labelsize=5.6)
    cbar.set_ticks([VMIN, -0.02, 0, 0.02, VMAX])
    cbar.set_ticklabels([f"{VMIN}", "-0.02", "0", "+0.02", f"+{VMAX:g}"])

    ax_leg.text(0.0, 0.46, "Significance", fontsize=6.6, fontweight="bold", color="#1A1A1A")
    ax_leg.text(0.0, 0.52, "*   p < 0.05  (boxed)", fontsize=5.8, color="#1A1A1A")
    ax_leg.text(0.0, 0.58, ".   p < 0.10", fontsize=5.8, color="#444444")

    ax_leg.text(0.0, 0.68, "Bias", fontsize=6.6, fontweight="bold", color="#1A1A1A")
    ax_leg.add_patch(mpatches.Rectangle((0.0, 0.72), 0.10, 0.026, facecolor=LA_DEEP, edgecolor="none"))
    ax_leg.text(0.14, 0.733, "LA-biased", fontsize=5.8)
    ax_leg.add_patch(mpatches.Rectangle((0.0, 0.76), 0.10, 0.026, facecolor=RA_DEEP, edgecolor="none"))
    ax_leg.text(0.14, 0.773, "RA-biased", fontsize=5.8)

    ax_leg.text(0.0, 0.86, "LA baseline carrier\n(Fig 1F HCA)", fontsize=6.0, fontweight="bold", color="#1A1A1A")
    present_carriers = list(dict.fromkeys([carriers[fid] for fid in fam_order]))
    ly = 0.91
    for c in present_carriers:
        if c == MIXED_LABEL:
            ax_leg.add_patch(mpatches.Rectangle((0.0, ly), 0.08, 0.024,
                                                facecolor=MIXED_COLOR, edgecolor="#888", linewidth=0.4, hatch="///"))
            ax_leg.text(0.12, ly + 0.012, "Mixed", fontsize=5.5, va="center")
        else:
            ax_leg.add_patch(mpatches.Rectangle((0.0, ly), 0.08, 0.024,
                                                facecolor=COMPARTMENT_PALETTE[c], edgecolor="none"))
            ax_leg.text(0.12, ly + 0.012, COMPARTMENT_SHORT[c], fontsize=5.5, va="center")
        ly += 0.028

    # ---- bottom histogram of all 144 dScore values ----
    flat = mat[np.isfinite(mat)]
    bins = np.arange(-0.05, 0.0501, 0.005)
    counts, edges = np.histogram(flat, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bar_colors = [cmap(norm(c)) for c in centers]
    ax_hist.bar(centers, counts, width=0.0045, color=bar_colors, edgecolor="#444444", linewidth=0.3)

    sig_mask = (pmat < 0.05) & np.isfinite(mat)
    sig_idx = np.argwhere(sig_mask)
    # Compose human-readable tags like "F15-Lymph", staggered vertically to
    # avoid overlap when multiple hits share a similar x.
    ct_short_for_tag = {c: COMPARTMENT_SHORT[c][:4] for c in COMPARTMENTS}
    sig_records = []
    for ii, jj in sig_idx:
        fid = fam_order[ii]
        ct = COMPARTMENTS[jj]
        v = float(mat[ii, jj])
        sig_records.append((v, f"F{int(fid)}-{ct_short_for_tag[ct]}"))
    sig_records.sort(key=lambda r: r[0])
    y_top = max(counts) * 1.05
    for k, (v, tag) in enumerate(sig_records):
        # Stagger label heights so adjacent values don't overprint.
        y_label = y_top * (0.55 + 0.18 * (k % 2))
        ax_hist.annotate(
            tag,
            xy=(v, 0), xytext=(v, y_label),
            ha="center", va="bottom",
            fontsize=5.6, color="#1A1A1A", fontweight="bold",
            arrowprops=dict(arrowstyle="-", color="#1A1A1A", lw=0.5, shrinkA=0, shrinkB=2),
        )
    ax_hist.set_xlim(-0.05, 0.05)
    ax_hist.set_ylim(0, max(counts) * 1.20)
    ax_hist.axvline(0, color="#888", lw=0.4, ls="--")
    ax_hist.set_xlabel("\u0394singscore  (AF \u2212 Ctrl)", fontsize=7.2)
    ax_hist.set_ylabel("# (family \u00d7 ct) cells", fontsize=6.6)
    ax_hist.tick_params(labelsize=6.0)
    for s in ("top", "right"):
        ax_hist.spines[s].set_visible(False)
    ax_hist.spines["left"].set_linewidth(0.4)
    ax_hist.spines["bottom"].set_linewidth(0.4)

    ax_hist_leg.axis("off")
    ax_hist_leg.set_xlim(0, 1); ax_hist_leg.set_ylim(0, 1)
    n_total = int(np.isfinite(mat).sum())
    n_sig = int(sig_mask.sum())
    pct_zero = 100.0 * (np.abs(flat) < 0.005).sum() / max(n_total, 1)
    ax_hist_leg.text(
        0.0, 0.95,
        "Distribution of all\n"
        f"{n_total} (family \u00d7 ct)\n"
        f"\u0394 values.\n\n"
        f"{pct_zero:.0f}% within \u00b10.005;\n"
        f"only {n_sig} cells reach\n"
        f"p < 0.05 (boxed *).",
        fontsize=5.6, color="#444444", va="top",
    )

    fig.suptitle(
        "Fig 3b  AF \u0394singscore at LA single-cell resolution (16 family \u00d7 9 compartment, GSE255612)",
        x=0.04, y=0.972, ha="left", fontsize=8.8, fontweight="bold",
    )
    fig.text(
        0.04, 0.012,
        "Rows: 16 NCP_v3 families ordered by Atlas PriorityScore (descending).  "
        "Cells: AF\u2212Ctrl difference of family-mean rank-singscore on Fig 1F-derived "
        "leading-edge genes (per-program LE then averaged across programs in family).  "
        "Statistics: per-cell donor-permutation, n = 18 AF + 16 Ctrl LA donors.  "
        "Carrier strip on the left shows the HCA paired-LA baseline executor (Fig 1F).  "
        "Yellow-bordered cells with * mark cell-types where AF significantly modifies "
        "the family score (nominal p < 0.05; \".\" = p < 0.10).  "
        "No |\u0394| value exceeds 0.04 in this dataset (peak |\u0394| = 0.031, F15 \u00d7 lymphocyte); "
        "color scale chosen to maximise visibility of the four sig hits.",
        fontsize=5.6, color="#444444",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
