#!/usr/bin/env python3
"""
N1 — Data-driven discovery of REWIRING-DEFINED gene modules via PCA on
the donor-level chamber asymmetry (aCAI) matrix.

Design (CRF-native discovery primitive):
    aCAI_g(d) = log2(score_LA(d, g)) - log2(score_RA(d, g))
    M_aCAI    = matrix indexed by (gene, donor)
                rows  = ~12,000 genes (intersection across 4 cohorts)
                cols  = ~50 donors (8 healthy + 17 SR + 25 AF, all 4 cohorts)

    PCA on the gene-side of M_aCAI returns data-driven AXES OF
    CHAMBER-REWIRING that are independent of MSigDB pathway taxonomy.
    Genes that load strongly together on these axes are coordinately
    rewired, regardless of whether MSigDB groups them into the same
    pathway.  Clustering on the top-K PC loadings therefore yields
    REWIRING-BEHAVIOUR-DEFINED gene modules, the discovery primitive
    that conventional GSEA-on-DE cannot produce.

Pipeline:
    1. Build (gene x donor) aCAI matrix from healthy (GSE112339) +
       3 paired chamber cohorts (GSE128188, GSE79768, GSE115574).
    2. Restrict to genes with non-trivial chamber bias in at least one
       cohort and to donors with both LA + RA samples.
    3. Centre and scale features (gene-wise unit variance), PCA, keep
       top-K (default 10) PCs.
    4. KMeans clustering (default K_MOD = 12) on the gene loadings.
    5. Per module annotate:
         - mean aCAI across healthy donors (signed direction)
         - mean ΔaCAI(AF − SR) across AF + SR paired donors (rewiring)
         - Jaccard overlap with each program in `normal_chamber_programs`
         - 'canonical' (Jaccard >= 0.30 with any panel program) or
           'novel'    (max Jaccard < 0.30)
         - top representative genes (largest absolute PC1+PC2 loading)

Outputs:
    results/n1_modules_long.tsv               (gene -> module id)
    results/n1_modules_summary.tsv            (per-module characterisation)
    results/n1_pc_loadings.tsv                (gene x top-K PC loadings)
    figures/figN1_rewiring_modules.pdf        (PCA scatter + module heatmap)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (
    COL_DOUBLE, LA_COLOR, RA_COLOR, panel_label,
    set_nature_style, strip_spines,
)

# import data loaders from N2
from n2_chamber_reversal_genes import (
    load_healthy_GSE112339,
    load_paired_cohort,
    per_donor_aCAI,
)

set_nature_style()

NCP_PROGRAM_PATH = ROOT / "bulk" / "Project_P01" / "normal_chamber_programs.tsv"

OUT_LONG    = ROOT / "results" / "n1_modules_long.tsv"
OUT_SUMMARY = ROOT / "results" / "n1_modules_summary.tsv"
OUT_LOAD    = ROOT / "results" / "n1_pc_loadings.tsv"
FIG_OUT     = ROOT / "figures" / "figN1_rewiring_modules.pdf"

K_PCS    = 10
K_MOD    = 12
RNG_SEED = 0
MIN_ABS_HEALTHY_aCAI = 0.20   # gene must have |aCAI_healthy| >= this in healthy
MIN_DONORS_PER_GENE  = 30      # gene must be measurable in at least this many donors


def jaccard(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def main() -> None:
    print("Loading healthy GSE112339 ...")
    h_expr, h_pairs = load_healthy_GSE112339()
    h_aCAI = per_donor_aCAI(h_expr, h_pairs)
    h_aCAI.columns = [f"healthy__{c}" for c in h_aCAI.columns]
    print(f"  healthy aCAI: {h_aCAI.shape[0]} genes x {h_aCAI.shape[1]} donors")

    cohort_donor_meta: list[dict] = []
    for d in h_aCAI.columns:
        cohort_donor_meta.append({
            "donor_col": d, "cohort": "GSE112339",
            "rhythm": "healthy",
        })

    paired_aCAI = []
    for cohort in ["GSE128188", "GSE79768", "GSE115574"]:
        print(f"Loading paired cohort {cohort} ...")
        expr, pairs = load_paired_cohort(cohort)
        aCAI = per_donor_aCAI(expr, pairs)
        # rename donor columns with cohort+rhythm prefix
        col_map = {}
        for p in pairs:
            la, ra, rhythm = p
            old_col = f"{la}__{ra}"
            if old_col not in aCAI.columns:
                continue
            new_col = f"{cohort}__{rhythm}__{la}_{ra}"
            col_map[old_col] = new_col
            cohort_donor_meta.append({
                "donor_col": new_col, "cohort": cohort, "rhythm": rhythm,
            })
        aCAI = aCAI.rename(columns=col_map)
        aCAI = aCAI[[c for c in aCAI.columns if c in col_map.values()]]
        print(f"  {cohort}: {aCAI.shape[0]} genes x {aCAI.shape[1]} donors")
        paired_aCAI.append(aCAI)

    # Common gene set across all 4 cohorts
    common = set(h_aCAI.index)
    for a in paired_aCAI:
        common &= set(a.index)
    common = sorted(common)
    print(f"\nCommon gene symbols across 4 cohorts: {len(common)}")

    # Concatenate column-wise into one (gene x donor) matrix
    M = pd.concat(
        [h_aCAI.loc[common]] + [a.loc[common] for a in paired_aCAI],
        axis=1,
    )
    cohort_donor = pd.DataFrame(cohort_donor_meta).set_index("donor_col")
    cohort_donor = cohort_donor.loc[M.columns]
    print(f"M_aCAI matrix: {M.shape}  "
          f"(healthy={(cohort_donor['rhythm']=='healthy').sum()}, "
          f"SR={(cohort_donor['rhythm']=='SR').sum()}, "
          f"AF={(cohort_donor['rhythm']=='AF').sum()})")

    healthy_cols = cohort_donor[cohort_donor["rhythm"] == "healthy"].index.tolist()
    h_mean = M[healthy_cols].mean(axis=1)
    keep_mask = (h_mean.abs() >= MIN_ABS_HEALTHY_aCAI) \
                & (M.notna().sum(axis=1) >= MIN_DONORS_PER_GENE)
    M_keep = M.loc[keep_mask].fillna(0.0)
    print(f"After gene filter (|aCAI_healthy| >= {MIN_ABS_HEALTHY_aCAI} "
          f"AND measured in >= {MIN_DONORS_PER_GENE} donors): "
          f"{M_keep.shape[0]} genes")

    print("Running PCA (gene-side, top-{0} PCs) ...".format(K_PCS))
    scaler = StandardScaler(with_mean=True, with_std=True)
    X = scaler.fit_transform(M_keep.values)        # gene x donor (centred & scaled)
    pca = PCA(n_components=K_PCS, random_state=RNG_SEED)
    pca.fit(X)
    # Per-gene loading on each PC (project gene back into PC space)
    loadings = pca.transform(X)                    # gene x K
    var_ratio = pca.explained_variance_ratio_
    print(f"  variance explained (top {K_PCS} PCs): "
          f"{[f'{v*100:.1f}%' for v in var_ratio]}; "
          f"total = {var_ratio.sum()*100:.1f}%")

    loadings_df = pd.DataFrame(
        loadings,
        index=M_keep.index,
        columns=[f"PC{i+1}" for i in range(K_PCS)],
    )
    loadings_df.to_csv(OUT_LOAD, sep="\t", index_label="gene")
    print(f"wrote {OUT_LOAD.relative_to(ROOT)}")

    print(f"KMeans clustering (K = {K_MOD}) on top-{K_PCS} PC loadings ...")
    km = KMeans(n_clusters=K_MOD, n_init=20, random_state=RNG_SEED)
    cluster = km.fit_predict(loadings)
    loadings_df["module"] = cluster

    progs = pd.read_csv(NCP_PROGRAM_PATH, sep="\t")
    progs.columns = [c.strip().strip('"') for c in progs.columns]
    progs["program_id"] = progs["program_id"].astype(str).str.strip().str.strip('"')
    progs["program"] = progs["program"].astype(str).str.strip().str.strip('"')
    progs["leading_edge_genes"] = progs["leading_edge_genes"] \
        .astype(str).str.strip().str.strip('"')
    program_sets = {
        row["program"]: set(g.strip() for g in row["leading_edge_genes"].split(","))
        for _, row in progs.iterrows()
    }

    af_donor_cols = cohort_donor[cohort_donor["rhythm"] == "AF"].index.tolist()
    sr_donor_cols = cohort_donor[cohort_donor["rhythm"] == "SR"].index.tolist()

    summaries = []
    long_rows = []
    for mod_id in sorted(loadings_df["module"].unique()):
        genes_in_mod = loadings_df.index[loadings_df["module"] == mod_id].tolist()
        sub = M_keep.loc[genes_in_mod]
        # signed module direction in healthy
        h_module = sub[healthy_cols].mean(axis=1).mean()
        af_module = sub[af_donor_cols].mean(axis=1).mean()
        sr_module = sub[sr_donor_cols].mean(axis=1).mean()
        # MSigDB best overlap
        best_program, best_jacc, best_n_overlap = None, 0.0, 0
        for prog, gset in program_sets.items():
            j = jaccard(set(genes_in_mod), gset)
            if j > best_jacc:
                best_jacc = j
                best_program = prog
                best_n_overlap = len(set(genes_in_mod) & gset)
        is_canonical = best_jacc >= 0.30
        # top representative genes (by absolute PC1+PC2 loading)
        top_loading = (loadings_df.loc[genes_in_mod, ["PC1", "PC2", "PC3"]]
                       .abs().sum(axis=1)
                       .sort_values(ascending=False))
        top_genes = top_loading.head(8).index.tolist()
        summaries.append({
            "module_id":            int(mod_id),
            "n_genes":              len(genes_in_mod),
            "healthy_aCAI_mean":    round(h_module, 4),
            "healthy_direction":    "LA-biased" if h_module > 0 else "RA-biased",
            "AF_aCAI_mean":         round(af_module, 4),
            "SR_aCAI_mean":         round(sr_module, 4),
            "delta_AF_minus_SR":    round(af_module - sr_module, 4),
            "best_NCP_program":     best_program if best_jacc > 0 else "",
            "best_NCP_jaccard":     round(best_jacc, 3),
            "best_NCP_n_overlap":   int(best_n_overlap),
            "annotation_class":     "canonical" if is_canonical else "novel",
            "top_genes":            ",".join(top_genes),
        })
        for g in genes_in_mod:
            long_rows.append({
                "gene":      g,
                "module_id": int(mod_id),
                "PC1":       round(float(loadings_df.loc[g, "PC1"]), 4),
                "PC2":       round(float(loadings_df.loc[g, "PC2"]), 4),
                "PC3":       round(float(loadings_df.loc[g, "PC3"]), 4),
            })

    summary = pd.DataFrame(summaries).sort_values(
        ["annotation_class", "best_NCP_jaccard", "n_genes"],
        ascending=[True, False, False],
    )
    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)}")

    long = pd.DataFrame(long_rows)
    long.to_csv(OUT_LONG, sep="\t", index=False)
    print(f"wrote {OUT_LONG.relative_to(ROOT)}")

    print("\n=== module summary ===")
    print(summary.to_string(index=False))
    n_canonical = (summary["annotation_class"] == "canonical").sum()
    n_novel     = (summary["annotation_class"] == "novel").sum()
    print(f"\ncanonical modules: {n_canonical}; novel modules: {n_novel}")

    fig = plt.figure(figsize=(COL_DOUBLE, 0.55 * COL_DOUBLE))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.50,
                          width_ratios=[1.2, 1.2, 1.0],
                          height_ratios=[1.0, 1.05])

    # panel a: gene-loading scatter PC1 vs PC2 colored by module
    ax = fig.add_subplot(gs[0, 0])
    cmap = plt.colormaps.get_cmap("tab20").resampled(K_MOD)
    for k in range(K_MOD):
        idx = loadings_df["module"] == k
        ax.scatter(loadings_df.loc[idx, "PC1"], loadings_df.loc[idx, "PC2"],
                   s=2.5, color=cmap(k), alpha=0.7,
                   edgecolor="none", label=f"M{k}")
    ax.axhline(0, color="#888", lw=0.4)
    ax.axvline(0, color="#888", lw=0.4)
    ax.set_xlabel(f"PC1 ({var_ratio[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_ratio[1]*100:.1f}%)")
    ax.set_title("a  Gene loadings on top-2 rewiring PCs",
                 loc="left", pad=4, fontsize=6.6)
    strip_spines(ax)
    panel_label(ax, "a", x=-0.18, y=1.04)

    # panel b: scree
    ax = fig.add_subplot(gs[0, 1])
    x = np.arange(1, K_PCS + 1)
    ax.bar(x, var_ratio * 100, color="#666666",
           edgecolor="white", linewidth=0.4)
    ax.set_xticks(x); ax.set_xticklabels([f"PC{i}" for i in x],
                                          rotation=0, fontsize=5)
    ax.set_ylabel("variance explained (%)")
    ax.set_title(f"b  Scree (top {K_PCS} PCs; "
                 f"cumulative {var_ratio.sum()*100:.1f}%)",
                 loc="left", pad=4, fontsize=6.6)
    strip_spines(ax)

    # panel c: module direction × novelty mini-bar
    ax = fig.add_subplot(gs[0, 2])
    cls_counts = summary.groupby(["annotation_class", "healthy_direction"]) \
                         .size().unstack(fill_value=0)
    cats = ["canonical", "novel"]
    la_vals = [cls_counts.loc[c, "LA-biased"] if c in cls_counts.index
               and "LA-biased" in cls_counts.columns else 0 for c in cats]
    ra_vals = [cls_counts.loc[c, "RA-biased"] if c in cls_counts.index
               and "RA-biased" in cls_counts.columns else 0 for c in cats]
    width = 0.35
    xpos = np.arange(len(cats))
    ax.bar(xpos - width/2, la_vals, width=width, color=LA_COLOR, label="LA-biased")
    ax.bar(xpos + width/2, ra_vals, width=width, color=RA_COLOR, label="RA-biased")
    for xi, la, ra in zip(xpos, la_vals, ra_vals):
        ax.text(xi - width/2, la + 0.05, str(la), ha="center",
                va="bottom", fontsize=5, fontweight="bold")
        ax.text(xi + width/2, ra + 0.05, str(ra), ha="center",
                va="bottom", fontsize=5, fontweight="bold")
    ax.set_xticks(xpos); ax.set_xticklabels(cats)
    ax.set_ylabel("modules (n)")
    ax.set_title(f"c  {n_canonical} canonical + {n_novel} novel modules",
                 loc="left", pad=4, fontsize=6.6)
    ax.legend(loc="upper right", frameon=False, fontsize=5)
    strip_spines(ax)

    # panel d (bottom span): module x donor heatmap of mean module aCAI
    ax = fig.add_subplot(gs[1, :])
    mod_donor = []
    for mod_id in sorted(loadings_df["module"].unique()):
        genes_in_mod = loadings_df.index[loadings_df["module"] == mod_id]
        mod_donor.append(M_keep.loc[genes_in_mod].mean(axis=0).rename(f"M{mod_id}"))
    mod_donor = pd.concat(mod_donor, axis=1).T  # module x donor
    # order donor columns by rhythm (healthy, SR, AF) then cohort
    rhy_pri = cohort_donor["rhythm"].map({"healthy":0, "SR":1, "AF":2})
    coh_pri = cohort_donor["cohort"].map(
        {"GSE112339":0, "GSE128188":1, "GSE79768":2, "GSE115574":3})
    order = (cohort_donor.assign(_r=rhy_pri, _c=coh_pri)
                          .sort_values(["_r", "_c"]).index.tolist())
    mod_donor = mod_donor[order]
    # order modules: by direction (LA first), then by ΔAF-SR magnitude
    summary_for_order = summary.assign(
        _dir=summary["healthy_direction"].map({"LA-biased":0, "RA-biased":1}),
        _absd=summary["delta_AF_minus_SR"].abs(),
    ).sort_values(["_dir", "_absd"], ascending=[True, False])
    mod_order_ids = [f"M{i}" for i in summary_for_order["module_id"].tolist()]
    mod_donor = mod_donor.loc[mod_order_ids]
    vmax = float(np.nanmax(np.abs(mod_donor.values))) or 0.001
    im = ax.imshow(mod_donor.values, cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(np.arange(len(mod_donor)))
    summary_idx = summary.set_index("module_id")
    ylabels = []
    for m in mod_donor.index:
        mid = int(m[1:])
        ann = summary_idx.loc[mid, "annotation_class"][:4]
        dir_ = summary_idx.loc[mid, "healthy_direction"]
        n   = summary_idx.loc[mid, "n_genes"]
        ylabels.append(f"{m}  ({ann}; {dir_}; n={n})")
    ax.set_yticklabels(ylabels, fontsize=5)
    ax.set_xticks([])
    # add rhythm colour stripe at top of heatmap
    rhy_colors = {"healthy":"#1B7837", "SR":"#2166AC", "AF":"#B2182B"}
    bar_h = 0.55
    y_bar = -0.85
    for xi, donor in enumerate(mod_donor.columns):
        rhy = cohort_donor.loc[donor, "rhythm"]
        ax.add_patch(mpatches.Rectangle(
            (xi - 0.5, y_bar), 1, bar_h,
            facecolor=rhy_colors[rhy], edgecolor="none",
            clip_on=False, zorder=3,
        ))
    ax.set_xlim(-0.5, len(mod_donor.columns) - 0.5)
    ax.set_ylim(len(mod_donor) - 0.5, y_bar - 0.2)
    # legend for rhythm bar
    handles = [
        mpatches.Patch(facecolor=rhy_colors["healthy"], label="healthy"),
        mpatches.Patch(facecolor=rhy_colors["SR"],      label="SR"),
        mpatches.Patch(facecolor=rhy_colors["AF"],      label="AF"),
    ]
    ax.legend(handles=handles, loc="upper right",
              bbox_to_anchor=(1.0, -0.04), ncols=3,
              frameon=False, fontsize=5)
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("module mean aCAI (LA − RA)", fontsize=5)
    cbar.ax.tick_params(labelsize=5)
    ax.set_title(
        f"d  Module × donor heatmap (donors ordered: healthy → SR → AF; "
        f"green/blue/red bar at top of heatmap)",
        loc="left", pad=4, fontsize=6.4,
    )
    panel_label(ax, "d", x=-0.06, y=1.06)

    fig.suptitle(
        "Data-driven discovery of REWIRING-DEFINED gene modules "
        "(PCA on the donor-level chamber asymmetry matrix)",
        y=1.00, fontsize=7.2, fontweight="bold",
    )

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
