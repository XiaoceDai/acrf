#!/usr/bin/env python3
"""
Figure 1f (REBUILD, baseline-only) — Cell-compartment attribution of the
30 NCP_v2 programs on the HEALTHY (Ctrl) donor subset of GSE255612, at
the full 13-compartment cell-type taxonomy.

Forward analytical step from Fig. 1e (baseline axis families):
    Question: at single-nucleus resolution, which compartment carries
    each NCP_v2 program at the healthy baseline?
    Method  : per program, per compartment, mean expression of program
              gene members in healthy (Ctrl) donor pseudobulks; row
              z-score across compartments; the compartment with the
              largest row z-score is the candidate carrier.

This panel is strictly healthy-baseline — AF biosamples are dropped.

Inputs:
    results/fig1c_ncp_panel_v2.tsv             (30 NCP_v2 programs)
    results/fig1e_baseline_axis_families.tsv   (program -> family id)
    bulk/Project_P01/gene_sets/h.all.v2026.1.Hs.symbols.gmt
    bulk/Project_P01/gene_sets/c2.cp.reactome.v2026.1.Hs.symbols.gmt
    Project_P3/GSE255612/pseudobulk_full/<compartment>/expr.tsv,samples.tsv
Output:
    figures/fig1f_baseline_attribution.pdf
    results/fig1f_baseline_attribution.tsv  (program x compartment z-score)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import matplotlib as mpl

from plot_style import COL_DOUBLE, panel_label, set_nature_style, strip_spines
from _pathway_labels import short_pathway_label

set_nature_style()

# Lock Arial (for unicode glyphs like β / ⁺) and updated chamber convention.
mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
})

LA_COLOR = "#C97A74"   # LA = red (new convention, low saturation)
RA_COLOR = "#7D99B8"   # RA = blue (low saturation)

NCP_TSV   = ROOT / "results" / "fig1c_ncp_panel_v2.tsv"
FAM_TSV   = ROOT / "results" / "fig1e_baseline_axis_families.tsv"
GMT_HALL  = ROOT / "bulk" / "Project_P01" / "gene_sets" \
                / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC  = ROOT / "bulk" / "Project_P01" / "gene_sets" \
                / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
PB_DIR    = ROOT / "Project_P3" / "GSE255612" / "pseudobulk_full"

OUT_TSV   = ROOT / "results" / "fig1f_baseline_attribution.tsv"
FIG_OUT   = ROOT / "figures" / "fig1f_baseline_attribution.pdf"

COMPARTMENTS = [
    "atrialCM", "fibroblast", "endothelial", "lymphatic_endothelial",
    "macrophage", "lymphocyte", "mast_cell", "B_cell",
    "pericyte", "vascular_smooth_muscle",
    "adipocyte", "epicardial", "neuronal",
]

# Muted per-compartment colors for the dominant-carrier side strip.
COMPARTMENT_PALETTE = [
    "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
    "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
    "#E15759", "#FF9D9A", "#79706E",
]
COMPARTMENT_COLOR_MAP = dict(zip(COMPARTMENTS, COMPARTMENT_PALETTE))

# Per-compartment minimum Ctrl-donor count: compartments with fewer than
# this many Ctrl donors are shown in the heatmap as N/A and excluded from
# the row-z computation and from "dominant carrier" calls.  This protects
# against single-donor compartments (e.g., B_cell n=1, epicardial n=1)
# producing spurious extreme z-scores.
MIN_CTRL_DONORS = 5

# Family palette (matches Fig 1e palette ordering)
FAMILY_PALETTE = [
    "#E6AB02", "#66A61E", "#7570B3", "#D95F02",
    "#1B9E77", "#E7298A", "#A6761D", "#666666",
    "#1F78B4", "#FB9A99", "#33A02C", "#B15928",
]


def load_gmt(path: Path, prefix: str) -> dict[str, list[str]]:
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            out[f"{prefix}__{parts[0]}"] = [g for g in parts[2:] if g]
    return out


short_label = lambda s: short_pathway_label(s, 34)


def main() -> None:
    print("Loading NCP_v2 panel + family assignments ...")
    ncp = pd.read_csv(NCP_TSV, sep="\t")
    fam = pd.read_csv(FAM_TSV, sep="\t")
    ncp = ncp.merge(fam[["pathway", "family_id"]], on="pathway", how="left")
    print(f"  programs: {len(ncp)} (with family assignments: "
          f"{ncp['family_id'].notna().sum()})")

    pathways = {**load_gmt(GMT_HALL, "HALLMARK"),
                **load_gmt(GMT_REAC, "REACTOME")}

    # ---- Per-compartment per-gene MEAN over healthy (Ctrl) donors only ----
    print("Computing per-compartment MEAN expression on Ctrl donors only ...")
    expr_means = {}
    n_donors_per_comp = {}
    for c in COMPARTMENTS:
        f = PB_DIR / c / "expr.tsv"
        m = PB_DIR / c / "samples.tsv"
        if not f.exists():
            print(f"  WARN: missing {f}")
            continue
        expr = pd.read_csv(f, sep="\t", index_col=0)
        meta = pd.read_csv(m, sep="\t")
        ctrl_donors = meta.loc[meta["rhythm"] == "Ctrl", "donor_id"].tolist()
        ctrl_present = [d for d in ctrl_donors if d in expr.columns]
        n_donors_per_comp[c] = len(ctrl_present)
        if len(ctrl_present) == 0:
            expr_means[c] = pd.Series(np.nan, index=expr.index)
        else:
            expr_means[c] = expr[ctrl_present].mean(axis=1)
    means_df = pd.concat(expr_means, axis=1)[COMPARTMENTS]
    print(f"  Ctrl-donor counts per compartment: "
          f"{[(c, n_donors_per_comp[c]) for c in COMPARTMENTS]}")

    # Compartments to include in the z-score computation (n_donors >= MIN_CTRL_DONORS)
    reliable_comps = [c for c in COMPARTMENTS
                      if n_donors_per_comp.get(c, 0) >= MIN_CTRL_DONORS]
    excluded_comps = [c for c in COMPARTMENTS if c not in reliable_comps]
    print(f"  RELIABLE  (n >= {MIN_CTRL_DONORS} Ctrl donors): "
          f"{len(reliable_comps)} compartments")
    print(f"  EXCLUDED  (n < {MIN_CTRL_DONORS}; shown as N/A in heatmap): "
          f"{excluded_comps}")

    # ---- Per program: row-z mean expression of LE genes across reliable compartments ----
    rows = []
    for _, r in ncp.iterrows():
        gset = pathways.get(r["pathway"], [])
        gset = [g for g in gset if g in means_df.index]
        if len(gset) < 5:
            continue
        sub = means_df.loc[gset]
        means_per_comp = sub.mean(axis=0)
        # z-score using only RELIABLE compartments
        rel_means = means_per_comp[reliable_comps]
        mu = rel_means.mean()
        sd = rel_means.std() + 1e-9
        rec = {"program_order": int(r["program_order"]),
                "pathway": r["pathway"],
                "direction": r["direction"],
                "family_id": int(r["family_id"]) if pd.notna(r["family_id"]) else -1,
                "n_genes_used": len(gset)}
        for c in COMPARTMENTS:
            if c in reliable_comps:
                rec[c] = float((means_per_comp[c] - mu) / sd)
            else:
                rec[c] = np.nan      # excluded for low coverage
        rows.append(rec)
    z_df = pd.DataFrame(rows)
    z_df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")

    # ---- Sort programs: by family_id (ordering matches Fig 1e ordering) ----
    fam_order = (fam.sort_values(["direction", "family_id",
                                    "leaf_order_within_dir"]
                                   if "leaf_order_within_dir" in fam.columns
                                   else ["direction", "family_id"])
                       ["pathway"].tolist())
    z_df["sort_idx"] = z_df["pathway"].map({p: i for i, p in enumerate(fam_order)})
    z_df = z_df.sort_values("sort_idx").reset_index(drop=True)

    n_la = (z_df["direction"] == "LA-biased").sum()
    n_ra = (z_df["direction"] == "RA-biased").sum()
    K = len(z_df)

    # Family colour map (consistent with Fig 1e)
    fam_ids_seen = sorted(z_df["family_id"].unique())
    fam_color_map = {fid: FAMILY_PALETTE[i % len(FAMILY_PALETTE)]
                      for i, fid in enumerate(fam_ids_seen)}

    # ---- Figure -----------------------------------------------------------
    fig = plt.figure(figsize=(COL_DOUBLE * 1.40, COL_DOUBLE * 1.18))
    gs = fig.add_gridspec(
        1, 5, width_ratios=[0.28, 0.30, 2.00, 4.90, 1.25], wspace=0.05,
        top=0.92, bottom=0.16, left=0.04, right=0.985,
    )

    mat = z_df[COMPARTMENTS].values.astype(float)
    vmax = float(np.nanmax(np.abs(mat))) or 0.001
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    # Mask N/A cells (low-donor compartments) so they render as grey
    mat_ma = np.ma.masked_invalid(mat)
    # Unified heatmap palette with Fig 6b (RA blue -> white -> LA red).
    from matplotlib.colors import LinearSegmentedColormap
    cmap_div = LinearSegmentedColormap.from_list(
        "LA_RA",
        [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)],
    )
    cmap_div.set_bad("#D8D8D8")

    # ---- col 0 : direction band with text (LA = red, RA = blue) ----
    ax_dir = fig.add_subplot(gs[0, 0])
    for i, d in enumerate(z_df["direction"]):
        ax_dir.add_patch(mpatches.Rectangle(
            (0.0, i - 0.5), 1.0, 1.0,
            facecolor=LA_COLOR if d == "LA-biased" else RA_COLOR,
            edgecolor="white", linewidth=0.2,
        ))
    ax_dir.text(0.5, n_la / 2 - 0.5, "LA-biased",
                rotation=90, ha="center", va="center",
                color="white", fontweight="bold", fontsize=9.2)
    ax_dir.text(0.5, n_la + n_ra / 2 - 0.5, "RA-biased",
                rotation=90, ha="center", va="center",
                color="white", fontweight="bold", fontsize=9.2)
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(K - 0.5, -0.5)
    ax_dir.set_xticks([]); ax_dir.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        ax_dir.spines[s].set_visible(False)

    # ---- col 1 : family colour band with F-id labels ----
    ax_fam = fig.add_subplot(gs[0, 1], sharey=ax_dir)
    for i in range(K):
        fid = z_df.iloc[i]["family_id"]
        ax_fam.add_patch(mpatches.Rectangle(
            (0.0, i - 0.5), 1.0, 1.0,
            facecolor=fam_color_map.get(fid, "#CCCCCC"),
            edgecolor="white", linewidth=0.2,
        ))
    fam_groups = z_df.groupby("family_id", sort=False)
    for fid, gp in fam_groups:
        idxs = z_df.index[z_df["family_id"] == fid].tolist()
        y_centre = np.mean(idxs)
        ax_fam.text(0.5, y_centre, f"F{int(fid)}",
                     fontsize=8.1, fontweight="bold",
                     ha="center", va="center", color="white")
    ax_fam.set_xlim(0, 1); ax_fam.set_ylim(K - 0.5, -0.5)
    ax_fam.set_xticks([]); ax_fam.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        ax_fam.spines[s].set_visible(False)

    # Pre-compute dominant carrier index per row for linked displays.
    dominant_idx = []
    for i in range(K):
        row_vals = mat[i]
        if np.all(np.isnan(row_vals)):
            dominant_idx.append(None)
        else:
            dominant_idx.append(int(np.nanargmax(row_vals)))

    # ---- col 2 : pathway labels + dominant strip (separate from heatmap) ---
    ax_lbl = fig.add_subplot(gs[0, 2], sharey=ax_dir)
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(K - 0.5, -0.5)
    ax_lbl.axis("off")
    for i, (_, r) in enumerate(z_df.iterrows()):
        dcol = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        didx = dominant_idx[i]
        if didx is not None:
            dcomp = COMPARTMENTS[didx]
            ax_lbl.add_patch(mpatches.Rectangle(
                (0.00, i - 0.32), 0.045, 0.64,
                facecolor=COMPARTMENT_COLOR_MAP.get(dcomp, "#999999"),
                edgecolor="white", linewidth=0.2,
            ))
        else:
            ax_lbl.add_patch(mpatches.Rectangle(
                (0.00, i - 0.32), 0.045, 0.64,
                facecolor="#D8D8D8", edgecolor="white", linewidth=0.2,
            ))
        ax_lbl.text(
            0.07, i, short_label(r["pathway"]),
            fontsize=8.3, ha="left", va="center",
            color=dcol, fontweight="normal",
        )

    # ---- col 3 : main attribution heatmap (occupies most of the panel) ---
    ax = fig.add_subplot(gs[0, 3], sharey=ax_dir)
    im = ax.imshow(mat_ma, cmap=cmap_div, norm=norm, aspect="auto",
                    extent=[-0.5, len(COMPARTMENTS) - 0.5, K - 0.5, -0.5])
    ax.set_yticks([])

    ax.set_xticks(np.arange(len(COMPARTMENTS)))
    comp_short = {
        "atrialCM": "aCM",
        "fibroblast": "Fibro",
        "endothelial": "Endo",
        "lymphatic_endothelial": "Lymph-endo",
        "macrophage": "Macro",
        "lymphocyte": "Lymph",
        "mast_cell": "Mast",
        "B_cell": "B",
        "pericyte": "Pericyte",
        "vascular_smooth_muscle": "vSMC",
        "adipocyte": "Adipo",
        "epicardial": "Epi",
        "neuronal": "Neuro",
    }
    xticklabels = []
    for c in COMPARTMENTS:
        n = n_donors_per_comp[c]
        base = comp_short.get(c, c)
        if c in reliable_comps:
            xticklabels.append(f"{base}\n(n={n})")
        else:
            xticklabels.append(f"{base}\n(n={n}, N/A)")
    ax.set_xticklabels(xticklabels, rotation=30, ha="right", fontsize=7.8)
    ax.set_xlabel("compartment", fontsize=9.2)

    # Dominant carrier marker is shown in the left strip (ax_lbl), not on heatmap.
    # LA / RA boundary
    ax.axhline(n_la - 0.5, color="black", lw=0.7, zorder=3)
    # Family boundary lines
    fam_changes = z_df["family_id"].ne(z_df["family_id"].shift())
    for i, change in enumerate(fam_changes):
        if i > 0 and change:
            ax.axhline(i - 0.5, color="black", lw=0.4, ls="-", alpha=0.6)
            ax_fam.axhline(i - 0.5, color="black", lw=0.4, ls="-",
                            alpha=0.6, clip_on=False)

    # ---- col 4 : dominant carrier short-name strip ----
    ax_dom = fig.add_subplot(gs[0, 4], sharey=ax_dir)
    ax_dom.set_xlim(0, 1); ax_dom.set_ylim(K - 0.5, -0.5)
    ax_dom.axis("off")
    ax_dom.set_title("carrier", fontsize=8.8, pad=2,
                      fontweight="bold")
    for i in range(K):
        row_vals = mat[i]
        if np.all(np.isnan(row_vals)):
            ax_dom.text(0.5, i, "\u2014", fontsize=8.0,
                         ha="center", va="center", color="#888")
            continue
        j = dominant_idx[i]
        if j is None:
            ax_dom.text(0.5, i, "\u2014", fontsize=8.0,
                         ha="center", va="center", color="#888")
            continue
        z_val = row_vals[j]
        comp_name = comp_short.get(COMPARTMENTS[j], COMPARTMENTS[j])
        # Color the carrier name dimly if z < 1.5 to flag weak attributions
        z_color = "black" if z_val >= 1.5 else "#888"
        z_weight = "bold" if z_val >= 2.0 else "normal"
        ax_dom.text(0.02, i, f"{comp_name}",
                     fontsize=8.0, ha="left", va="center", color=z_color,
                     fontweight=z_weight)
        ax_dom.text(0.98, i, f"z={z_val:+.1f}",
                     fontsize=7.5, ha="right", va="center", color="#666")

    # ---- shared colour bar at bottom ----
    cax = fig.add_axes([0.33, 0.05, 0.25, 0.017])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label("row z-score", fontsize=8.4)
    cbar.ax.tick_params(labelsize=7.8)
    cbar.ax.text(
        0.02, 0.5, "depleted ←", transform=cbar.ax.transAxes,
        ha="left", va="center", fontsize=7.8,
        color=RA_COLOR, fontweight="bold",
    )
    cbar.ax.text(
        0.98, 0.5, "→ enriched", transform=cbar.ax.transAxes,
        ha="right", va="center", fontsize=7.8,
        color=LA_COLOR, fontweight="bold",
    )
    cbar.ax.text(
        0.50, 1.75, "left strip = dominant carrier",
        transform=cbar.ax.transAxes,
        ha="center", va="bottom", fontsize=7.2, color="#333",
    )

    ax.set_title("f  Compartment attribution", loc="left",
                 fontsize=10.0, fontweight="bold", pad=4)
    # Detailed method notes are moved to figures/fig1_legend.txt.

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")

    # ---- console summary --------------------------------------------------
    print("\n=== Dominant carrier per program (z-score; reliable compartments only) ===")
    for _, r in z_df.iterrows():
        comp_vals = r[COMPARTMENTS].astype(float)
        if comp_vals.isna().all():
            continue
        j = int(np.nanargmax(comp_vals.values))
        print(f"  F{int(r['family_id'])}  {short_label(r['pathway']):42s} "
              f"-> {COMPARTMENTS[j]:25s} z = {comp_vals.iloc[j]:+.2f}")


if __name__ == "__main__":
    main()
