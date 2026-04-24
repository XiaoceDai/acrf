#!/usr/bin/env python3
"""
Figure 1e (SUMMARY version, for the main composite Figure 1).

Drops the per-donor aCAI heatmap to reduce visual fatigue across the
Figure 1 panels.  Keeps the structural elements that carry the message:

    - dendrogram (LA branch red, RA branch blue)
    - thin family colour band with F1..FK labels
    - family list with representative pathway and n_programs

The full per-donor heatmap version is kept as supplementary
(fig1e_baseline_axis_families.pdf).

Inputs:
    results/fig1e_baseline_axis_families.tsv   (cached by fig1e_baseline_axis_families.py)
    results/fig1d_baseline_donor_aCAI.tsv      (per-donor aCAI matrix; needed
                                                 to recompute Pearson distance
                                                 -> dendrogram for the panel)

Output:
    figures/fig1e_summary.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy as sch
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import COL_DOUBLE, set_nature_style
from _pathway_labels import short_pathway_label

set_nature_style()

mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
})

LA_COLOR = "#D6604D"
RA_COLOR = "#4393C3"
CUT_HEIGHT = 0.35

FAMILY_PALETTE = [
    "#E6AB02", "#66A61E", "#7570B3", "#D95F02",
    "#1B9E77", "#E7298A", "#A6761D", "#666666",
    "#1F78B4", "#FB9A99", "#33A02C", "#B15928",
]

FAM_TSV  = ROOT / "results" / "fig1e_baseline_axis_families.tsv"
ACAI_TSV = ROOT / "results" / "fig1d_baseline_donor_aCAI.tsv"
OUT      = ROOT / "figures" / "fig1e_summary.pdf"


def short_label(s: str) -> str:
    return short_pathway_label(s, 36)


def cluster_one_direction(df_dir: pd.DataFrame, donor_cols: list[str]):
    """Return (linkage Z, leaf order indices)."""
    if len(df_dir) < 2:
        return None, np.arange(len(df_dir))
    M = df_dir[donor_cols].values
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            r = np.corrcoef(M[i], M[j])[0, 1]
            d = 1 - r
            D[i, j] = d
            D[j, i] = d
    Z = sch.linkage(squareform(D, checks=False), method="average")
    return Z, sch.leaves_list(Z)


def main() -> None:
    fam = pd.read_csv(FAM_TSV, sep="\t")
    acai = pd.read_csv(ACAI_TSV, sep="\t")
    donor_cols = [c for c in acai.columns if c.startswith("aCAI_")]

    # Merge family info onto aCAI table so we have everything in one frame
    df = acai.merge(fam[["pathway", "family_id"]], on="pathway", how="inner")

    la = df[df["direction"] == "LA-biased"].copy().reset_index(drop=True)
    ra = df[df["direction"] == "RA-biased"].copy().reset_index(drop=True)

    Z_la, _ = cluster_one_direction(la, donor_cols)
    Z_ra, _ = cluster_one_direction(ra, donor_cols)

    # Reuse Fig 1e ordering: by family_id then by program_order so families stay grouped
    sorted_df = pd.concat([
        la.sort_values(["family_id", "program_order"]),
        ra.sort_values(["family_id", "program_order"]),
    ], ignore_index=True)

    n_la = len(la); n_ra = len(ra)
    K = len(sorted_df)

    fam_summary = (
        sorted_df.assign(abs_mean=sorted_df["mean_aCAI"].abs())
                 .sort_values("abs_mean", ascending=False)
                 .groupby("family_id", sort=False)
                 .agg(direction=("direction", "first"),
                      n_programs=("pathway", "size"),
                      representative=("pathway", "first"))
                 .reset_index()
    )
    fam_summary = fam_summary.sort_values("family_id").reset_index(drop=True)

    fam_ids_seen = sorted(sorted_df["family_id"].unique())
    fam_color_map = {fid: FAMILY_PALETTE[i % len(FAMILY_PALETTE)]
                     for i, fid in enumerate(fam_ids_seen)}
    fam_dir_map = dict(zip(fam_summary["family_id"], fam_summary["direction"]))

    # ---- Figure layout : dendrogram | family band | family list ----------
    fig = plt.figure(figsize=(COL_DOUBLE * 0.78, COL_DOUBLE * 0.85))
    gs = fig.add_gridspec(
        1, 3, width_ratios=[1.1, 0.32, 2.55], wspace=0.06,
        top=0.88, bottom=0.10, left=0.04, right=0.985,
    )

    # ---- col 0 : combined dendrogram (LA top, RA bottom) -----------------
    ax_dendro = fig.add_subplot(gs[0, 0])

    def draw_dendro(Z, n, y_offset, color):
        if Z is None:
            return
        sch.set_link_color_palette([color])
        d = sch.dendrogram(Z, no_plot=True, orientation="left",
                            labels=[None] * n)
        for icoord, dcoord in zip(d["icoord"], d["dcoord"]):
            x = [-v for v in dcoord]
            y = [(v / 10.0 - 0.5) + y_offset for v in icoord]
            ax_dendro.plot(x, y, color=color, lw=0.6)
        sch.set_link_color_palette(None)

    draw_dendro(Z_la, n_la, 0, LA_COLOR)
    draw_dendro(Z_ra, n_ra, n_la, RA_COLOR)

    ax_dendro.axvline(-CUT_HEIGHT, color="#1B7837", lw=0.9, ls="--")
    ax_dendro.text(-CUT_HEIGHT, -1.2,
                    f"cut h = {CUT_HEIGHT:.2f}",
                    fontsize=7, color="#1B7837",
                    ha="center", va="bottom", fontweight="bold")
    ax_dendro.axhline(n_la - 0.5, color="#888", lw=0.5, ls=":")

    # Direction tags inside the dendrogram column
    ax_dendro.text(-CUT_HEIGHT, n_la / 2, "LA branch",
                   fontsize=7.5, color=LA_COLOR, fontweight="bold",
                   ha="left", va="center", rotation=90)
    ax_dendro.text(-CUT_HEIGHT, n_la + n_ra / 2, "RA branch",
                   fontsize=7.5, color=RA_COLOR, fontweight="bold",
                   ha="left", va="center", rotation=90)

    ax_dendro.set_ylim(-2.0, K - 0.5)
    ax_dendro.invert_yaxis()
    ax_dendro.axis("off")
    ax_dendro.set_title("dendrogram\n(1 − Pearson r,\naverage linkage)",
                        fontsize=8, pad=2)

    # ---- col 1 : family colour band with F-id ----------------------------
    ax_fam = fig.add_subplot(gs[0, 1], sharey=ax_dendro)
    for i in range(K):
        fid = sorted_df.iloc[i]["family_id"]
        ax_fam.add_patch(mpatches.Rectangle(
            (0.05, i - 0.5), 0.9, 1.0,
            facecolor=fam_color_map[fid], edgecolor="white", linewidth=0.2,
        ))
    fam_groups = sorted_df.groupby("family_id", sort=False)
    for fid, gp in fam_groups:
        idxs = sorted_df.index[sorted_df["family_id"] == fid].tolist()
        y_centre = float(np.mean(idxs))
        ax_fam.text(0.5, y_centre, f"F{int(fid)}",
                     fontsize=7.5, fontweight="bold",
                     ha="center", va="center", color="white")
    ax_fam.set_xlim(0, 1); ax_fam.set_ylim(K - 0.5, -0.5)
    ax_fam.set_xticks([]); ax_fam.set_yticks([])
    ax_fam.set_title("family", fontsize=8, pad=2)
    for s in ("top", "right", "left", "bottom"):
        ax_fam.spines[s].set_visible(False)

    # Family-boundary horizontal lines on the family band (not dendrogram)
    fam_changes = sorted_df["family_id"].ne(sorted_df["family_id"].shift())
    for i, change in enumerate(fam_changes):
        if i > 0 and change:
            ax_fam.axhline(i - 0.5, color="black", lw=0.4, ls="-",
                            alpha=0.6, clip_on=False)
    ax_fam.axhline(n_la - 0.5, color="black", lw=1.0, zorder=5)

    # ---- col 2 : family list (representative pathway, n_programs) -------
    ax_leg = fig.add_subplot(gs[0, 2])
    ax_leg.set_xlim(0, 1)
    K_total = len(fam_ids_seen)
    ax_leg.set_ylim(K_total + 0.5, -0.5)
    ax_leg.axis("off")
    ax_leg.set_title("axis families (representative pathway)",
                     fontsize=9, pad=2, loc="left", fontweight="bold")
    for slot, fid in enumerate(fam_ids_seen):
        rep = fam_summary.loc[fam_summary["family_id"] == fid,
                              "representative"].iloc[0]
        n = fam_summary.loc[fam_summary["family_id"] == fid,
                            "n_programs"].iloc[0]
        col = fam_color_map[fid]
        chamber_col = LA_COLOR if fam_dir_map[fid] == "LA-biased" else RA_COLOR
        ax_leg.add_patch(mpatches.Rectangle(
            (0.0, slot - 0.35), 0.10, 0.7,
            facecolor=col, edgecolor="white", linewidth=0.4,
        ))
        ax_leg.text(0.05, slot, f"F{int(fid)}",
                     fontsize=7.5, fontweight="bold",
                     ha="center", va="center", color="white")
        ax_leg.text(0.13, slot, f"{short_label(rep)}",
                     fontsize=7.8, ha="left", va="center",
                     color=chamber_col, fontweight="bold")
        ax_leg.text(0.13, slot + 0.36, f"n = {int(n)} programs",
                     fontsize=6.2, ha="left", va="center",
                     color="#666", style="italic")

    # ---- title + footnote ------------------------------------------------
    n_la_fam = (fam_summary["direction"] == "LA-biased").sum()
    n_ra_fam = (fam_summary["direction"] == "RA-biased").sum()
    fig.suptitle(
        f"e  Hierarchical clustering identifies {n_la_fam + n_ra_fam} baseline axis families "
        f"({n_la_fam} LA + {n_ra_fam} RA)",
        x=0.04, y=0.96, ha="left", va="top",
        fontsize=9.5, fontweight="bold",
    )
    fig.text(
        0.04, 0.92,
        "Programs are clustered by their cross-donor aCAI vector "
        "(1 \u2212 Pearson, average linkage); LA / RA branches clustered separately.",
        fontsize=7.2, color="#444", ha="left", va="top",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
