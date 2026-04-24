#!/usr/bin/env python3
"""
Figure 1e (REBUILD, baseline-only) — Hierarchical clustering of the 30
NCP_v2 programs by their cross-donor aCAI patterns identifies the
baseline AXIS FAMILIES of healthy chamber asymmetry.

Input  : results/fig1d_baseline_donor_aCAI.tsv  (30 programs x 8 healthy
         donors aCAI matrix produced by fig1d_baseline_donor_aCAI.py)
Method : pairwise distance = 1 - Pearson correlation of each program's
         cross-donor aCAI vector; agglomerative hierarchical clustering
         with average linkage; cut the tree at a distance threshold to
         define K data-driven axis families.  LA-biased and RA-biased
         programs are clustered separately so chamber direction is
         preserved.
Outputs:
    results/fig1e_baseline_axis_families.tsv   (program -> family id)
    figures/fig1e_baseline_axis_families.pdf   (dendrogram + family heatmap)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster import hierarchy as sch
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import matplotlib as mpl

from plot_style import COL_DOUBLE, panel_label, set_nature_style, strip_spines
from _pathway_labels import short_pathway_label, family_short_label

set_nature_style()

# Lock Arial first so any unicode (β, ⁺) is taken from Arial directly.
mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
})

# Updated chamber convention: RA = blue, LA = red.
LA_COLOR = "#D6604D"
RA_COLOR = "#4393C3"

IN  = ROOT / "results" / "fig1d_baseline_donor_aCAI.tsv"
OUT_TSV = ROOT / "results" / "fig1e_baseline_axis_families.tsv"
FIG_OUT = ROOT / "figures" / "fig1e_baseline_axis_families.pdf"

CUT_HEIGHT = 0.35   # 1 - Pearson correlation threshold


# Family colour palette (12 distinct, colourblind-safe)
FAMILY_PALETTE = [
    "#E6AB02", "#66A61E", "#7570B3", "#D95F02",
    "#1B9E77", "#E7298A", "#A6761D", "#666666",
    "#1F78B4", "#FB9A99", "#33A02C", "#B15928",
]


short_label = lambda s: short_pathway_label(s, 30)


def cluster_one_direction(df_dir: pd.DataFrame, donor_cols: list[str],
                          fam_id_offset: int = 0) -> pd.DataFrame:
    if len(df_dir) < 2:
        df_dir = df_dir.copy()
        df_dir["family_id"] = fam_id_offset + 1
        df_dir["leaf_order"] = 0
        return df_dir

    M = df_dir[donor_cols].values
    # 1 - Pearson correlation
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            r = np.corrcoef(M[i], M[j])[0, 1]
            d = 1 - r
            D[i, j] = d
            D[j, i] = d
    Z = sch.linkage(squareform(D, checks=False), method="average")
    fam = sch.fcluster(Z, t=CUT_HEIGHT, criterion="distance")
    leaf_order = sch.leaves_list(Z)

    df = df_dir.copy().reset_index(drop=True)
    df["family_id"] = fam + fam_id_offset
    df["leaf_order_within_dir"] = np.argsort(np.argsort(leaf_order))
    df["_link_z"] = [Z] * len(df)  # carry linkage for plotting
    return df


def main() -> None:
    print("Loading baseline aCAI matrix ...")
    df = pd.read_csv(IN, sep="\t")
    donor_cols = [c for c in df.columns if c.startswith("aCAI_")]
    donors = [c.replace("aCAI_", "") for c in donor_cols]
    print(f"  programs: {len(df)}  donors: {len(donors)}")

    # Cluster LA and RA separately so chamber direction is preserved
    la = df[df["direction"] == "LA-biased"].copy()
    ra = df[df["direction"] == "RA-biased"].copy()
    la_clu = cluster_one_direction(la, donor_cols, fam_id_offset=0)
    n_la_fam = la_clu["family_id"].nunique()
    ra_clu = cluster_one_direction(ra, donor_cols, fam_id_offset=n_la_fam)
    print(f"  LA families: {n_la_fam}  "
          f"RA families: {ra_clu['family_id'].nunique()}")

    # Order programs : LA first (by family then by leaf order),
    # then RA (by family then by leaf order)
    la_clu = la_clu.sort_values(["family_id", "leaf_order_within_dir"])
    ra_clu = ra_clu.sort_values(["family_id", "leaf_order_within_dir"])

    sorted_df = pd.concat([la_clu, ra_clu], ignore_index=True)

    # Keep linkage objects for dendrograms before they get dropped
    Z_la = la_clu["_link_z"].iloc[0] if len(la_clu) else None
    Z_ra = ra_clu["_link_z"].iloc[0] if len(ra_clu) else None
    sorted_df = sorted_df.drop(columns=["_link_z"])

    # Per-family summary : pick representative = highest |mean_aCAI|
    sorted_df["abs_mean_aCAI"] = sorted_df["mean_aCAI"].abs()
    fam_summary = (sorted_df.sort_values("abs_mean_aCAI", ascending=False)
                            .groupby("family_id")
                            .agg(direction=("direction", "first"),
                                 n_programs=("pathway", "size"),
                                 representative=("pathway", "first"),
                                 mean_aCAI_rep=("mean_aCAI", "first"))
                            .reset_index())
    print("\n=== axis families ===")
    print(fam_summary.to_string(index=False))

    # Save program-to-family table
    out_cols = (["program_order", "pathway", "direction",
                 "family_id", "n_genes_used",
                 "mean_aCAI", "sign_concordance"]
                + donor_cols)
    sorted_df[out_cols].to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\nwrote {OUT_TSV.relative_to(ROOT)}")

    # ---- Figure -----------------------------------------------------------
    # Build a per-row family label (representative pathway short name)
    fam_rep_map = dict(zip(fam_summary["family_id"],
                            fam_summary["representative"]))
    fam_n_map   = dict(zip(fam_summary["family_id"],
                            fam_summary["n_programs"]))
    family_lookup = (sorted_df.groupby("family_id")
                              .agg(rep_pathway=("pathway", "first"),
                                   start_idx=("pathway", lambda s: s.index[0])
                                   if False else ("pathway", "first"))
                              .reset_index())

    fig = plt.figure(figsize=(COL_DOUBLE * 1.35, COL_DOUBLE * 0.95))
    # Drop the redundant direction strip (already encoded by sort + dendrogram colour)
    # 4 columns instead of 5: dendrogram, heatmap, family band, legend
    gs = fig.add_gridspec(1, 4,
                          width_ratios=[1.0, 4.0, 0.45, 2.4],
                          wspace=0.05,
                          top=0.90, bottom=0.13, left=0.03, right=0.985)

    # ---- left : combined dendrogram (LA on top, RA on bottom) -----------
    ax_dendro = fig.add_subplot(gs[0, 0])
    n_la = len(la_clu); n_ra = len(ra_clu)

    def draw_dendro(Z, n, y_offset, color):
        if Z is None:
            return
        sch.set_link_color_palette([color])
        d = sch.dendrogram(Z, no_plot=True, orientation="left",
                            labels=[None] * n)
        for icoord, dcoord in zip(d["icoord"], d["dcoord"]):
            x = [-v for v in dcoord]
            y = [(v / 10.0 - 0.5) + y_offset for v in icoord]
            ax_dendro.plot(x, y, color=color, lw=0.5)
        sch.set_link_color_palette(None)

    draw_dendro(Z_la, n_la, 0, LA_COLOR)
    draw_dendro(Z_ra, n_ra, n_la, RA_COLOR)
    # cut-height line (vertical, since orientation is left)
    ax_dendro.axvline(-CUT_HEIGHT, color="#1B7837", lw=0.9, ls="--")
    ax_dendro.text(-CUT_HEIGHT, -1.4,
                    f"cut h = {CUT_HEIGHT:.2f}",
                    fontsize=6, color="#1B7837",
                    ha="center", va="bottom", fontweight="bold")
    ax_dendro.axhline(n_la - 0.5, color="#888", lw=0.5, ls=":")
    ax_dendro.set_ylim(-2.0, n_la + n_ra - 0.5)
    ax_dendro.invert_yaxis()
    ax_dendro.axis("off")
    ax_dendro.set_title("dendrogram", fontsize=8.8, pad=2)
    # Tiny direction key inside the dendrogram column
    ax_dendro.text(-CUT_HEIGHT, n_la / 2, "LA",
                   fontsize=7.8, color=LA_COLOR, fontweight="bold",
                   ha="left", va="center", rotation=90)
    ax_dendro.text(-CUT_HEIGHT, n_la + n_ra / 2, "RA",
                   fontsize=7.8, color=RA_COLOR, fontweight="bold",
                   ha="left", va="center", rotation=90)

    # ---- centre : aCAI heatmap ----
    # Diverging colormap aligned with new chamber convention
    # (positive aCAI -> LA -> red,  negative -> RA -> blue).
    from matplotlib.colors import LinearSegmentedColormap
    acai_cmap = LinearSegmentedColormap.from_list(
        "LA_RA",
        [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)],
    )

    ax = fig.add_subplot(gs[0, 1], sharey=ax_dendro)
    mat = sorted_df[donor_cols].values
    vmax = float(np.nanmax(np.abs(mat))) or 0.001
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax.imshow(mat, cmap=acai_cmap, norm=norm, aspect="auto",
                    extent=[-0.5, len(donors) - 0.5,
                             len(sorted_df) - 0.5, -0.5])
    ax.set_yticks(np.arange(len(sorted_df)))
    ax.set_yticklabels(
        [short_label(r["pathway"]) for _, r in sorted_df.iterrows()],
        fontsize=8.1,
    )
    # Color y-tick labels by chamber direction for quick orientation
    for tl, direction in zip(ax.get_yticklabels(), sorted_df["direction"]):
        tl.set_color(LA_COLOR if direction == "LA-biased" else RA_COLOR)
    ax.set_xticks(np.arange(len(donors)))
    ax.set_xticklabels(donors, rotation=45, ha="right", fontsize=8.2)
    ax.set_xlabel("donor", fontsize=9.2)
    ax.set_title("e  Axis families", loc="left", pad=6,
                 fontsize=10.0, fontweight="bold")

    # ---- col 2 : family colour with F-id label ----
    ax_fam = fig.add_subplot(gs[0, 2], sharey=ax)
    fam_ids_seen = sorted(sorted_df["family_id"].unique())
    fam_color_map = {fid: FAMILY_PALETTE[i % len(FAMILY_PALETTE)]
                     for i, fid in enumerate(fam_ids_seen)}
    for i in range(len(sorted_df)):
        fid = sorted_df.iloc[i]["family_id"]
        ax_fam.add_patch(mpatches.Rectangle(
            (0.05, i - 0.5), 0.9, 1.0,
            facecolor=fam_color_map[fid], edgecolor="white", linewidth=0.2,
        ))
    fam_groups = sorted_df.groupby("family_id", sort=False)
    for fid, gp in fam_groups:
        idxs = sorted_df.index[sorted_df["family_id"] == fid].tolist()
        y_centre = np.mean(idxs)
        ax_fam.text(0.5, y_centre, f"F{int(fid)}",
                     fontsize=8.1, fontweight="bold",
                     ha="center", va="center", color="white")
    ax_fam.set_xlim(0, 1); ax_fam.set_ylim(len(sorted_df) - 0.5, -0.5)
    ax_fam.set_xticks([])
    ax_fam.set_yticks([])
    ax_fam.set_title("family", fontsize=8.8, pad=2)
    for s in ("top", "right", "left", "bottom"):
        ax_fam.spines[s].set_visible(False)

    # ---- col 3 : family legend (representative names) ----
    ax_leg = fig.add_subplot(gs[0, 3])
    ax_leg.set_xlim(0, 1)
    K_total = len(fam_ids_seen)
    ax_leg.set_ylim(K_total + 0.5, -0.5)
    ax_leg.axis("off")
    ax_leg.set_title("families",
                     fontsize=10.0, pad=2, loc="left", fontweight="bold")
    fam_dir_map = dict(zip(fam_summary["family_id"], fam_summary["direction"]))
    for slot, fid in enumerate(fam_ids_seen):
        rep = fam_rep_map[fid]
        n   = fam_n_map[fid]
        col = fam_color_map[fid]
        chamber_col = LA_COLOR if fam_dir_map[fid] == "LA-biased" else RA_COLOR
        # Large coloured square with F-id
        ax_leg.add_patch(mpatches.Rectangle(
            (0.0, slot - 0.35), 0.10, 0.7,
            facecolor=col, edgecolor="white", linewidth=0.4,
        ))
        ax_leg.text(0.05, slot, f"F{int(fid)}",
                     fontsize=8.1, fontweight="bold",
                     ha="center", va="center", color="white")
        ax_leg.text(0.13, slot, f"{short_label(rep)}",
                     fontsize=8.1, ha="left", va="center",
                     color=chamber_col, fontweight="bold")
        ax_leg.text(0.13, slot + 0.36, f"n = {int(n)} programs",
                     fontsize=6.8, ha="left", va="center",
                     color="#666", style="italic")

    # Family-boundary horizontal lines on the heatmap
    fam_changes = sorted_df["family_id"].ne(sorted_df["family_id"].shift())
    for i, change in enumerate(fam_changes):
        if i > 0 and change:
            ax.axhline(i - 0.5, color="black", lw=0.5, ls="-")
            ax_fam.axhline(i - 0.5, color="black", lw=0.5, ls="-",
                            clip_on=False)
    # LA / RA boundary line (heavier)
    ax.axhline(n_la - 0.5, color="black", lw=1.0, zorder=5)

    # Colorbar for heatmap (use Figure space below)
    cax = fig.add_axes([0.18, 0.05, 0.22, 0.018])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label("aCAI", fontsize=8.4)
    cbar.ax.tick_params(labelsize=7.8)
    cbar.ax.text(
        0.02, 0.5, "RA-biased ←", transform=cbar.ax.transAxes,
        ha="left", va="center", fontsize=7.8,
        color=RA_COLOR, fontweight="bold",
    )
    cbar.ax.text(
        0.98, 0.5, "→ LA-biased", transform=cbar.ax.transAxes,
        ha="right", va="center", fontsize=7.8,
        color=LA_COLOR, fontweight="bold",
    )

    # Detailed method notes are moved to figures/fig1_legend.txt.

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
