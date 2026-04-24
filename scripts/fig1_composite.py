#!/usr/bin/env python3
"""
Fig 1 publication composite — Cell / Nature / Science style.
3 × 2 grid; panels D / F aggregated to FAMILY level so heatmap cells
are roughly square and the figure is compact.

Layout:
    a  PCA donors                  b  Volcano + chamber-biased gene cohort
    c  Curated NCP selection       (continued)
    d  Family-level aCAI heatmap   e  Axis family overview (dendrogram + list)
    f  Family-level cell attribution

All panels share fonts, line weights, colour palette.  Heat-maps follow
Cell / Nature / Science conventions: continuous colour, no per-cell
borders, no in-cell text, square-ish aspect, thin colour bars.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster import hierarchy as sch
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label

OUT = ROOT / "figures" / "fig1_composite.pdf"

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
MM = 1.0 / 25.4
COL_DOUBLE = 183 * MM
FIG_HEIGHT = 195 * MM      # 6 panels in 2 rows; row 1 compact, row 2 spacious

LA_COLOR = "#B2182B"
RA_COLOR = "#2166AC"
HL_GREEN = "#1B7837"
GREY_LIGHT = "#D8D8D8"
GREY_MID   = "#888888"
DIVERGING_CMAP = "RdBu_r"

FAMILY_PALETTE = [
    "#E6AB02", "#66A61E", "#7570B3", "#D95F02",
    "#1B9E77", "#E7298A", "#A6761D", "#666666",
    "#1F78B4", "#FB9A99", "#33A02C", "#B15928",
]

PANEL_LETTER_KW = dict(fontsize=11, fontweight="bold",
                        ha="left", va="top")

short = lambda s, n=28: short_pathway_label(s, n)


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family":     "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         7,
        "axes.titlesize":    7.5,
        "axes.labelsize":    7,
        "xtick.labelsize":   6.5,
        "ytick.labelsize":   6.5,
        "legend.fontsize":   6,
        "axes.linewidth":    0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size":  2.5,
        "ytick.major.size":  2.5,
        "lines.linewidth":   0.8,
        "patch.linewidth":   0.4,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "figure.dpi":         150,
        "savefig.dpi":        300,
        "pdf.fonttype":       42,
        "ps.fonttype":        42,
        "svg.fonttype":       "none",
        "legend.frameon":     False,
    })


def strip_axes(ax) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_letter(fig, x, y, label):
    fig.text(x, y, label, transform=fig.transFigure, **PANEL_LETTER_KW)


# ---------------------------------------------------------------------------
# Panel A — PCA donors
# ---------------------------------------------------------------------------
def draw_panel_a(ax) -> None:
    qc = pd.read_csv(ROOT / "bulk" / "Project_P01" / "qc" / "qc_pca.tsv",
                      sep="\t")
    for d, sub in qc.groupby("donor"):
        if len(sub) == 2:
            ax.plot(sub["PC1"], sub["PC2"], color=GREY_MID,
                    lw=0.4, alpha=0.6, zorder=1)
    for _, r in qc.iterrows():
        col = LA_COLOR if r["chamber"] == "LA" else RA_COLOR
        ax.scatter(r["PC1"], r["PC2"], color=col, edgecolor="white",
                    linewidth=0.5, s=24, zorder=3)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    ax.set_title("Healthy LA / RA RNA-seq\n(GSE112339, n = 8 paired)",
                  loc="left", pad=8, fontsize=7, fontweight="bold")
    handles = [
        plt.Line2D([], [], marker="o", color="white",
                    markerfacecolor=LA_COLOR, markeredgecolor="white",
                    markersize=5, label="LA"),
        plt.Line2D([], [], marker="o", color="white",
                    markerfacecolor=RA_COLOR, markeredgecolor="white",
                    markersize=5, label="RA"),
        plt.Line2D([], [], color=GREY_MID, lw=0.6, alpha=0.6,
                    label="paired donor"),
    ]
    ax.legend(handles=handles, loc="best", fontsize=5.5)
    strip_axes(ax)


# ---------------------------------------------------------------------------
# Panel B — Volcano
# ---------------------------------------------------------------------------
def draw_panel_b(ax) -> None:
    de = pd.read_csv(ROOT / "bulk" / "Project_P01" / "bulk_chamber_DE_full.tsv",
                      sep="\t")
    de = de.dropna(subset=["logFC", "FDR"])
    de["nlog10"] = -np.log10(de["FDR"].clip(lower=1e-300))
    fdr_thr, fc_thr = 0.05, np.log2(1.5)
    bg = (de["FDR"] >= fdr_thr) | (de["logFC"].abs() <= fc_thr)
    la = (~bg) & (de["logFC"] > 0)
    ra = (~bg) & (de["logFC"] < 0)

    ax.scatter(de.loc[bg, "logFC"], de.loc[bg, "nlog10"],
                s=0.8, c=GREY_LIGHT, alpha=0.5, edgecolor="none", zorder=1)
    ax.scatter(de.loc[la, "logFC"], de.loc[la, "nlog10"],
                s=1.4, c=LA_COLOR, alpha=0.9, edgecolor="none", zorder=2,
                label=f"LA-biased  n = {int(la.sum())}")
    ax.scatter(de.loc[ra, "logFC"], de.loc[ra, "nlog10"],
                s=1.4, c=RA_COLOR, alpha=0.9, edgecolor="none", zorder=2,
                label=f"RA-biased  n = {int(ra.sum())}")
    ax.axhline(-np.log10(fdr_thr), color=GREY_MID, lw=0.4, ls="--")
    ax.axvline( fc_thr, color=GREY_MID, lw=0.4, ls="--")
    ax.axvline(-fc_thr, color=GREY_MID, lw=0.4, ls="--")

    canonical = ["PITX2", "BMP10", "NPPA", "NPPB", "MYL7"]
    for gene in canonical:
        m = de[de["gene_symbol"].str.strip('"') == gene]
        if len(m) and m["FDR"].iloc[0] < fdr_thr:
            x = float(m["logFC"].iloc[0]); y = float(m["nlog10"].iloc[0])
            ax.scatter([x], [y], s=18, facecolor="none",
                        edgecolor="black", linewidth=0.5, zorder=4)
            ax.annotate(gene, xy=(x, y), xytext=(5, 4),
                         textcoords="offset points",
                         fontsize=6, fontweight="bold", color="black",
                         arrowprops=dict(arrowstyle="-", color="black",
                                          lw=0.4, shrinkA=2))

    ax.set_xlabel("log₂ FC (LA − RA)")
    ax.set_ylabel("−log₁₀ FDR")
    ax.set_title("Paired LA-vs-RA limma\n→ chamber-biased gene cohort",
                  loc="left", pad=8, fontsize=7, fontweight="bold")
    ax.legend(loc="upper right", fontsize=5.5)
    strip_axes(ax)


# ---------------------------------------------------------------------------
# Panel C — Curated NCP_v2 selection (no numbers; colour-only)
# ---------------------------------------------------------------------------
def draw_panel_c(ax_sc, ax_lst) -> None:
    long  = pd.read_csv(ROOT / "results" / "fig1c_pathway_selection_long.tsv",
                         sep="\t")
    panel = pd.read_csv(ROOT / "results" / "fig1c_ncp_panel_v2.tsv",
                         sep="\t")
    panel_set = set(panel["pathway"])
    long["is_panel"] = long["pathway"].isin(panel_set)

    # Background candidates: stability colour ramp (no edges, no numbers)
    stab_colors = ["#E5E5E5", "#C9C9C9", "#AEAEAE", "#959595",
                   "#7E7E7E", "#666666", "#525252", "#3D3D3D", "#222222"]
    bg = long[~long["is_panel"]]
    ax_sc.scatter(bg["NES"], bg["gsva_logFC"], s=10,
                   c=[stab_colors[min(int(s), 8)]
                      for s in bg["stability_consistent_runs"]],
                   edgecolor="none", alpha=0.85, zorder=2)

    # NCP panel: clean colour-only markers (no numbers, no in-plot names)
    pnl = long[long["is_panel"]]
    la_pts = pnl[pnl["NES"] > 0]; ra_pts = pnl[pnl["NES"] < 0]
    ax_sc.scatter(la_pts["NES"], la_pts["gsva_logFC"],
                   s=22, c=LA_COLOR, edgecolor="black", linewidth=0.4,
                   zorder=4, label=f"NCP LA  n = {len(la_pts)}")
    ax_sc.scatter(ra_pts["NES"], ra_pts["gsva_logFC"],
                   s=22, c=RA_COLOR, edgecolor="black", linewidth=0.4,
                   zorder=4, label=f"NCP RA  n = {len(ra_pts)}")

    ax_sc.axhline(0, color=GREY_MID, lw=0.4)
    ax_sc.axvline(0, color=GREY_MID, lw=0.4)
    ax_sc.set_xlabel("signed GSEA NES (LA-biased > 0)")
    ax_sc.set_ylabel("GSVA paired log₂FC (LA − RA)")
    ax_sc.set_title(
        f"Dual-method GSEA + GSVA + stability + Jaccard\n"
        f"{len(long)} candidates → {(long['stability_consistent_runs']>=7).sum()} stable "
        f"→ {len(panel)} NCP panel",
        loc="left", pad=8, fontsize=7, fontweight="bold")

    # Stability legend
    handles = [
        plt.Line2D([], [], marker="o", color="white",
                    markerfacecolor=stab_colors[0], markeredgecolor="none",
                    markersize=4, label="0 / 8 stable"),
        plt.Line2D([], [], marker="o", color="white",
                    markerfacecolor=stab_colors[4], markeredgecolor="none",
                    markersize=4, label="4 / 8"),
        plt.Line2D([], [], marker="o", color="white",
                    markerfacecolor=stab_colors[8], markeredgecolor="none",
                    markersize=4, label="8 / 8"),
    ]
    leg1 = ax_sc.legend(handles=handles, loc="lower left",
                          fontsize=5.5, title="background",
                          title_fontsize=6, labelspacing=0.4)
    ax_sc.add_artist(leg1)
    handles2 = [
        plt.Line2D([], [], marker="o", color="white",
                    markerfacecolor=LA_COLOR, markeredgecolor="black",
                    markersize=6, label=f"NCP LA  ({len(la_pts)})"),
        plt.Line2D([], [], marker="o", color="white",
                    markerfacecolor=RA_COLOR, markeredgecolor="black",
                    markersize=6, label=f"NCP RA  ({len(ra_pts)})"),
    ]
    ax_sc.legend(handles=handles2, loc="upper left",
                  fontsize=5.5, labelspacing=0.4)
    strip_axes(ax_sc)

    # ---- Right panel: clean colour-coded list, LA in two columns on top,
    # RA in single column at bottom.  Better balance for K=24+6.
    la_panel = panel[panel["NES"] > 0].sort_values("NES", ascending=False).reset_index(drop=True)
    ra_panel = panel[panel["NES"] < 0].sort_values("NES", ascending=True).reset_index(drop=True)
    n_la = len(la_panel); n_ra = len(ra_panel)
    half = (n_la + 1) // 2          # split LA into 2 columns
    K_rows = half + 2 + n_ra        # +2 = headers
    ax_lst.set_xlim(0, 1); ax_lst.set_ylim(K_rows + 0.5, -0.5)
    ax_lst.axis("off")
    ax_lst.text(0.0, -0.2,
                  f"NCP_v2 panel  (K = {len(panel)})",
                  fontsize=7, fontweight="bold",
                  ha="left", va="bottom")

    # ---- LA section header ----
    ax_lst.text(0.0,  0.45, f"LA-biased  ({n_la})",
                  fontsize=6.0, fontweight="bold",
                  color=LA_COLOR, ha="left", va="bottom")
    ax_lst.text(0.55, 0.45, "(continued)",
                  fontsize=5.0, fontstyle="italic",
                  color=LA_COLOR, ha="left", va="bottom")

    # LA in 2 columns; label cap at 17 chars to avoid bleeding into right col
    for i in range(n_la):
        in_left = i < half
        x_col = 0.0 if in_left else 0.55
        y_row = (i if in_left else (i - half)) + 1.0
        ax_lst.add_patch(mpatches.Circle(
            (x_col + 0.020, y_row), radius=0.18,
            facecolor=LA_COLOR, edgecolor="none"))
        ax_lst.text(x_col + 0.060, y_row,
                      short(la_panel.iloc[i]["pathway"], 17),
                      ha="left", va="center", fontsize=5.5,
                      color="black")

    # ---- RA section header (below LA block) ----
    ra_y0 = half + 1.5
    ax_lst.text(0.0, ra_y0 - 0.15, f"RA-biased  ({n_ra})",
                  fontsize=6.5, fontweight="bold",
                  color=RA_COLOR, ha="left", va="bottom")
    for i in range(n_ra):
        y_row = ra_y0 + 0.5 + i
        ax_lst.add_patch(mpatches.Circle(
            (0.018, y_row), radius=0.16,
            facecolor=RA_COLOR, edgecolor="none"))
        ax_lst.text(0.05, y_row,
                      short(ra_panel.iloc[i]["pathway"], 26),
                      ha="left", va="center", fontsize=5.8,
                      color="black")


# ---------------------------------------------------------------------------
# Helpers — family-level aggregation
# ---------------------------------------------------------------------------
def load_family_membership():
    """Return df with one row per program: pathway, direction, family_id."""
    fam = pd.read_csv(ROOT / "results" / "fig1e_baseline_axis_families.tsv",
                       sep="\t")
    return fam[["pathway", "direction", "family_id", "program_order",
                  "n_genes_used"]
                + [c for c in fam.columns if c.startswith("aCAI_")]]


def family_summary():
    """Per family: representative pathway, direction, n_programs."""
    fam = load_family_membership()
    donor_cols = [c for c in fam.columns if c.startswith("aCAI_")]
    fam["mean_aCAI"] = fam[donor_cols].mean(axis=1)
    fam["abs_mean_aCAI"] = fam["mean_aCAI"].abs()
    summary = (fam.sort_values("abs_mean_aCAI", ascending=False)
                  .groupby("family_id")
                  .agg(direction=("direction", "first"),
                        n_programs=("pathway", "size"),
                        representative=("pathway", "first"))
                  .reset_index())
    return summary, fam, donor_cols


# ---------------------------------------------------------------------------
# Panel D — Family-level aCAI heatmap (11 families × 8 donors, square cells)
# ---------------------------------------------------------------------------
def draw_panel_d(ax_main, ax_cbar) -> None:
    """Family-level baseline aCAI heatmap.

    Clean Cell/Nature/Science style:
      - no separate side bands (F-id encoded as colored bullet in
        y-tick label text rendered with a small marker prefix)
      - inline vertical colorbar on the right
      - generous font sizes
      - LA/RA boundary as single bold horizontal line
    """
    summary, fam, donor_cols = family_summary()
    donors = [c.replace("aCAI_", "") for c in donor_cols]

    fam_mat = (fam.groupby("family_id")[donor_cols].mean()
                  .reindex(summary["family_id"]))
    summary["row_mean"] = fam_mat.mean(axis=1).values
    la_sorted = summary[summary["direction"] == "LA-biased"] \
                 .sort_values("row_mean", ascending=False)
    ra_sorted = summary[summary["direction"] == "RA-biased"] \
                 .sort_values("row_mean", ascending=True)
    sorted_summary = pd.concat([la_sorted, ra_sorted]).reset_index(drop=True)
    fam_mat = fam_mat.reindex(sorted_summary["family_id"])
    n_la = (sorted_summary["direction"] == "LA-biased").sum()

    # Family colour map (consistent across panels)
    fam_ids = sorted(summary["family_id"])
    fam_color_map = {fid: FAMILY_PALETTE[i % len(FAMILY_PALETTE)]
                      for i, fid in enumerate(fam_ids)}

    mat = fam_mat.values
    vmax = float(np.nanmax(np.abs(mat))) or 0.001
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax_main.imshow(mat, cmap=DIVERGING_CMAP, norm=norm,
                          aspect="auto", interpolation="nearest",
                          extent=[-0.5, len(donors) - 0.5,
                                   len(sorted_summary) - 0.5, -0.5])

    # Y ticks: pathway name only (no F-id - shown as coloured marker)
    ax_main.set_yticks(np.arange(len(sorted_summary)))
    ax_main.set_yticklabels(
        [short(r["representative"], 24)
         for _, r in sorted_summary.iterrows()],
        fontsize=7)

    # Coloured family-id bullets immediately to the LEFT of each y-tick label
    # (drawn in axes coordinates, just outside the heatmap area)
    for i, (_, r) in enumerate(sorted_summary.iterrows()):
        fid = int(r["family_id"])
        ax_main.add_patch(mpatches.Circle(
            (-1.4, i), radius=0.35,
            facecolor=fam_color_map[fid], edgecolor="white",
            linewidth=0.5, clip_on=False, zorder=10))
        ax_main.text(-1.4, i, f"F{fid}",
                       fontsize=5.0, fontweight="bold",
                       color="white", ha="center", va="center",
                       clip_on=False, zorder=11)

    ax_main.set_xticks(np.arange(len(donors)))
    ax_main.set_xticklabels(donors, rotation=45, ha="right", fontsize=6.5)
    ax_main.set_xlabel("healthy donor (GSE112339)", fontsize=7)
    ax_main.tick_params(length=2, pad=2)
    ax_main.axhline(n_la - 0.5, color="black", lw=1.0, zorder=4)
    for s in ("top", "right", "left", "bottom"):
        ax_main.spines[s].set_visible(False)

    ax_main.set_title(
        f"Family-level baseline aCAI · {len(sorted_summary)} families × 8 donors",
        loc="left", pad=8, fontsize=7.5, fontweight="bold")

    # Inline vertical colorbar
    cbar = mpl.colorbar.ColorbarBase(ax_cbar, cmap=plt.colormaps[DIVERGING_CMAP],
                                       norm=norm, orientation="vertical")
    cbar.set_label("aCAI (LA − RA, log₂)", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6, length=2, pad=1)
    cbar.outline.set_linewidth(0.5)


# ---------------------------------------------------------------------------
# Panel E — Axis families overview (dendrogram + family list)
# ---------------------------------------------------------------------------
def draw_panel_e(ax) -> None:
    """Compact axis-family overview: dendrogram with bullet+name leaves.

    Single panel: dendrogram on left, leaf labels on right with
    coloured F# bullets and family representative names + n_programs.
    """
    summary, fam, donor_cols = family_summary()

    fam_mat = (fam.groupby("family_id")[donor_cols].mean()
                  .reindex(summary["family_id"]))
    dir_map = dict(zip(summary["family_id"], summary["direction"]))

    def fam_dendro(direction):
        ids = [fid for fid in fam_mat.index if dir_map[fid] == direction]
        if len(ids) < 2:
            return None, ids
        M = fam_mat.loc[ids].values
        n = M.shape[0]
        D = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                r = np.corrcoef(M[i], M[j])[0, 1]
                D[i, j] = D[j, i] = 1 - r
        Z = sch.linkage(squareform(D, checks=False), method="average")
        return Z, ids

    Z_la, ids_la = fam_dendro("LA-biased")
    Z_ra, ids_ra = fam_dendro("RA-biased")
    n_la = len(ids_la); n_ra = len(ids_ra)

    fam_color_map = {fid: FAMILY_PALETTE[i % len(FAMILY_PALETTE)]
                      for i, fid in enumerate(sorted(summary["family_id"]))}

    # Get leaf order from dendrograms
    leaf_order_la = list(sch.leaves_list(Z_la)) if Z_la is not None else [0]
    leaf_order_ra = list(sch.leaves_list(Z_ra)) if Z_ra is not None else [0]
    ordered_ids = ([ids_la[i] for i in leaf_order_la]
                    + [ids_ra[i] for i in leaf_order_ra])

    # Draw dendrogram on left third of axis; labels on right two-thirds
    DENDRO_X_FRAC = 0.30   # fraction of x-range used by dendrogram
    K = len(ordered_ids)

    # Find max distance to scale dendrogram x
    max_d = 0.0
    if Z_la is not None: max_d = max(max_d, Z_la[:, 2].max())
    if Z_ra is not None: max_d = max(max_d, Z_ra[:, 2].max())
    max_d = max(max_d, 0.5)

    def draw_one(Z, ids, color, y_offset):
        if Z is None or len(ids) < 2:
            return
        sch.set_link_color_palette([color])
        d = sch.dendrogram(Z, no_plot=True, orientation="left",
                            labels=ids)
        for icoord, dcoord in zip(d["icoord"], d["dcoord"]):
            # rescale dcoord (height) to [0, DENDRO_X_FRAC]
            x = [DENDRO_X_FRAC * (1.0 - v / max_d) for v in dcoord]
            y = [(v / 10.0 - 0.5) + y_offset for v in icoord]
            ax.plot(x, y, color=color, lw=0.9)
        sch.set_link_color_palette(None)

    draw_one(Z_la, ids_la, LA_COLOR, 0)
    draw_one(Z_ra, ids_ra, RA_COLOR, n_la)

    # Cut-line at h = 0.35 mapped into the rescaled coordinate
    CUT = 0.35
    cut_x = DENDRO_X_FRAC * (1.0 - CUT / max_d)
    ax.axvline(cut_x, color=HL_GREEN, lw=1.0, ls="--",
                ymin=0.05, ymax=0.95)
    ax.text(cut_x, -1.0, f"cut\nh = {CUT:.2f}",
             fontsize=5.5, color=HL_GREEN,
             ha="center", va="bottom", fontweight="bold")

    # LA / RA boundary horizontal divider
    ax.axhline(n_la - 0.5, color="black", lw=0.6, alpha=0.5,
                xmin=0.0, xmax=1.0)

    # Leaf labels (coloured bullet + family text)
    summary_idx = summary.set_index("family_id")
    for slot, fid in enumerate(ordered_ids):
        col = fam_color_map[fid]
        rep = summary_idx.loc[fid, "representative"]
        n_progs = int(summary_idx.loc[fid, "n_programs"])
        # Bullet
        ax.add_patch(mpatches.Circle(
            (DENDRO_X_FRAC + 0.04, slot), radius=0.32,
            facecolor=col, edgecolor="white", linewidth=0.5))
        ax.text(DENDRO_X_FRAC + 0.04, slot, f"F{int(fid)}",
                  fontsize=5.0, fontweight="bold",
                  color="white", ha="center", va="center")
        # Pathway name + n_programs on a single line
        ax.text(DENDRO_X_FRAC + 0.10, slot,
                  f"{short(rep, 22)}",
                  fontsize=6.5, ha="left", va="center", color="black")
        ax.text(0.99, slot, f"n = {n_progs}",
                  fontsize=5.5, ha="right", va="center",
                  color="#666", style="italic")

    ax.set_xlim(0, 1)
    ax.set_ylim(K + 0.5, -1.6)
    ax.axis("off")
    ax.set_title(
        f"Axis families ({K})  ·  hierarchical clustering of family aCAI profiles",
        fontsize=7.5, pad=4, loc="left", fontweight="bold")


# ---------------------------------------------------------------------------
# Panel F — Family-level cell-compartment attribution
# ---------------------------------------------------------------------------
def draw_panel_f(ax_main, ax_cbar) -> None:
    """Family-level cell-compartment attribution heatmap.

    Same clean style as Panel D: F-id encoded as coloured bullet on
    each y-tick label, inline vertical colorbar on right.
    """
    z = pd.read_csv(ROOT / "results" / "fig1f_baseline_attribution.tsv",
                     sep="\t")
    summary, fam, donor_cols = family_summary()
    z = z.merge(fam[["pathway"]], on="pathway", how="inner")
    COMPARTMENTS = [
        "atrialCM", "fibroblast", "endothelial", "lymphatic_endothelial",
        "macrophage", "lymphocyte", "mast_cell", "B_cell",
        "pericyte", "vascular_smooth_muscle",
        "adipocyte", "epicardial", "neuronal",
    ]
    fam_z = (z.groupby("family_id")[COMPARTMENTS].mean()
                .reindex(summary["family_id"]))

    summary["row_mean"] = fam_z.replace([np.inf, -np.inf], np.nan).mean(axis=1).values
    la_sorted = summary[summary["direction"] == "LA-biased"] \
                 .sort_values("family_id")
    ra_sorted = summary[summary["direction"] == "RA-biased"] \
                 .sort_values("family_id")
    sorted_summary = pd.concat([la_sorted, ra_sorted]).reset_index(drop=True)
    fam_z = fam_z.reindex(sorted_summary["family_id"])
    n_la = (sorted_summary["direction"] == "LA-biased").sum()

    fam_ids = sorted(summary["family_id"])
    fam_color_map = {fid: FAMILY_PALETTE[i % len(FAMILY_PALETTE)]
                      for i, fid in enumerate(fam_ids)}

    mat = fam_z.values.astype(float)
    vmax = float(np.nanmax(np.abs(mat))) or 0.001
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap_div = plt.colormaps[DIVERGING_CMAP].copy()
    cmap_div.set_bad("#D8D8D8")
    mat_ma = np.ma.masked_invalid(mat)

    im = ax_main.imshow(mat_ma, cmap=cmap_div, norm=norm,
                          aspect="auto", interpolation="nearest",
                          extent=[-0.5, len(COMPARTMENTS) - 0.5,
                                   len(sorted_summary) - 0.5, -0.5])

    ax_main.set_yticks(np.arange(len(sorted_summary)))
    ax_main.set_yticklabels(
        [short(r["representative"], 24)
         for _, r in sorted_summary.iterrows()],
        fontsize=7)

    # Coloured family-id bullets to the LEFT of each y-tick label
    for i, (_, r) in enumerate(sorted_summary.iterrows()):
        fid = int(r["family_id"])
        ax_main.add_patch(mpatches.Circle(
            (-2.0, i), radius=0.35,
            facecolor=fam_color_map[fid], edgecolor="white",
            linewidth=0.5, clip_on=False, zorder=10))
        ax_main.text(-2.0, i, f"F{fid}",
                       fontsize=5.0, fontweight="bold",
                       color="white", ha="center", va="center",
                       clip_on=False, zorder=11)

    ax_main.set_xticks(np.arange(len(COMPARTMENTS)))
    short_comp = {
        "atrialCM": "atrialCM",
        "fibroblast": "fibroblast",
        "endothelial": "endothelial",
        "lymphatic_endothelial": "lymphatic-endo",
        "macrophage": "macrophage",
        "lymphocyte": "lymphocyte",
        "mast_cell": "mast cell",
        "B_cell": "B cell",
        "pericyte": "pericyte",
        "vascular_smooth_muscle": "VSMC",
        "adipocyte": "adipocyte",
        "epicardial": "epicardial",
        "neuronal": "neuronal",
    }
    ax_main.set_xticklabels(
        [short_comp.get(c, c) for c in COMPARTMENTS],
        rotation=45, ha="right", fontsize=6.5)
    ax_main.set_xlabel("compartment (GSE255612 Ctrl donors)", fontsize=7)
    ax_main.tick_params(length=2, pad=2)
    ax_main.axhline(n_la - 0.5, color="black", lw=1.0, zorder=4)

    # Dominant carrier black box
    for i in range(len(sorted_summary)):
        row_vals = mat[i]
        if np.all(np.isnan(row_vals)):
            continue
        j = int(np.nanargmax(row_vals))
        ax_main.add_patch(mpatches.Rectangle(
            (j - 0.45, i - 0.45), 0.90, 0.90,
            facecolor="none", edgecolor="black", linewidth=0.9, zorder=5))
    for s in ("top", "right", "left", "bottom"):
        ax_main.spines[s].set_visible(False)
    ax_main.set_title(
        f"Family attribution · {len(sorted_summary)} families × 13 compartments "
        "(black box = dominant carrier)",
        loc="left", pad=8, fontsize=7.5, fontweight="bold")

    cbar = mpl.colorbar.ColorbarBase(ax_cbar, cmap=cmap_div, norm=norm,
                                       orientation="vertical")
    cbar.set_label("row z-score", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6, length=2, pad=1)
    cbar.outline.set_linewidth(0.5)


# ---------------------------------------------------------------------------
# Master composite
# ---------------------------------------------------------------------------
def main() -> None:
    set_publication_style()
    fig = plt.figure(figsize=(COL_DOUBLE, FIG_HEIGHT))
    # Manual margins so leftmost y-tick labels do not get clipped
    # Generous left + bottom margins so y-tick and x-tick labels do not clip
    fig.subplots_adjust(left=0.10, right=0.985, top=0.965, bottom=0.075)

    outer = fig.add_gridspec(
        2, 1,
        height_ratios=[0.55, 1.20],   # row 1 is compact summary, row 2 dominant
        hspace=0.55,
    )

    # ---- ROW 1 : a + b + c ---------------------------------------------
    row1 = outer[0].subgridspec(
        1, 4, width_ratios=[0.95, 1.00, 1.10, 1.30], wspace=0.55,
    )
    ax_a = fig.add_subplot(row1[0, 0])
    ax_b = fig.add_subplot(row1[0, 1])
    ax_c1 = fig.add_subplot(row1[0, 2])
    ax_c2 = fig.add_subplot(row1[0, 3])
    draw_panel_a(ax_a)
    draw_panel_b(ax_b)
    draw_panel_c(ax_c1, ax_c2)
    panel_letter(fig, 0.005, 0.985, "a")
    panel_letter(fig, 0.230, 0.985, "b")
    panel_letter(fig, 0.460, 0.985, "c")

    # ---- ROW 2 : d + e + f, each with inline vertical colour bar -------
    row2 = outer[1].subgridspec(
        1, 3,
        width_ratios=[0.95, 1.00, 1.45],
        wspace=0.55,
    )

    # Panel d : main heatmap + inline vertical colour bar
    d_inner = row2[0, 0].subgridspec(
        1, 2, width_ratios=[1.0, 0.045], wspace=0.06)
    ax_d_main = fig.add_subplot(d_inner[0, 0])
    ax_d_cbar = fig.add_subplot(d_inner[0, 1])
    draw_panel_d(ax_d_main, ax_d_cbar)
    panel_letter(fig, 0.005, 0.555, "d")

    # Panel e : single axis with dendrogram + leaf labels
    ax_e = fig.add_subplot(row2[0, 1])
    draw_panel_e(ax_e)
    panel_letter(fig, 0.345, 0.555, "e")

    # Panel f : main heatmap + inline vertical colour bar
    f_inner = row2[0, 2].subgridspec(
        1, 2, width_ratios=[1.0, 0.030], wspace=0.04)
    ax_f_main = fig.add_subplot(f_inner[0, 0])
    ax_f_cbar = fig.add_subplot(f_inner[0, 1])
    draw_panel_f(ax_f_main, ax_f_cbar)
    panel_letter(fig, 0.620, 0.555, "f")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
