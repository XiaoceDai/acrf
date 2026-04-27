"""
S5 — Cell/Nature-style LA vs RA atrial composition figures (PDF).

Purely descriptive (no differential-abundance testing).
Slice: source=="Nuclei" & Used=="Yes" & region in {LA, RA}   (matches DE).
Excludes: doublets, NotAssigned, Ventricular_Cardiomyocyte.

Outputs (figures/, all PDF):
    fig_composition_overall.pdf
    fig_composition_per_donor.pdf
    fig_composition_top4_paired.pdf
    composition_LA_RA_table.tsv
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plot_style import (  # noqa: E402
    COL_DOUBLE,
    MM,
    panel_label,
    set_nature_style,
    strip_spines,
)

H5AD = Path("/home/bio/Project_P02/global_raw.h5ad")
FIG_DIR = Path("/home/bio/Project_P02/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

TOP4 = {
    "Atrial_Cardiomyocyte":    "atrial_cardiomyocyte",
    "Fibroblast":              "fibroblast",
    "Endothelial":             "endothelial_or_vascular",
    "Myeloid":                 "immune",
    "Lymphoid":                "immune",
}
EXCLUDE_CT = {"doublets", "NotAssigned", "Ventricular_Cardiomyocyte"}

CT_PALETTE = [
    "#1f6fb4", "#9ac6e5", "#ef9020", "#f7c286", "#2ca02c", "#98df8a",
    "#d62728", "#ff9896", "#8e5ea2", "#c5b0d5", "#7f7f7f", "#bcbd22", "#17becf",
]
TOP4_PALETTE = {
    "atrial_cardiomyocyte":    "#1f6fb4",
    "fibroblast":              "#9ac6e5",
    "endothelial_or_vascular": "#ef9020",
    "immune":                  "#2ca02c",
}


def _stacked_hbar(ax, order, values_df, palette, label_map=None,
                  edgecolor="white", edgewidth=0.3):
    left = np.zeros(len(values_df))
    y = np.arange(len(values_df))
    for ct in order:
        vals = values_df[ct].values
        ax.barh(y, vals, left=left, color=palette[ct],
                edgecolor=edgecolor, linewidth=edgewidth,
                label=(label_map[ct] if label_map else ct))
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(values_df.index)
    ax.invert_yaxis()


def main() -> None:
    set_nature_style()

    a = ad.read_h5ad(H5AD, backed="r")
    obs = a.obs[["donor", "region", "cell_type", "source", "Used"]].copy()
    obs = obs[
        (obs["source"].astype(str) == "Nuclei")
        & obs["region"].isin(["LA", "RA"])
        & (obs["Used"].astype(str) == "Yes")
        & (~obs["cell_type"].astype(str).isin(EXCLUDE_CT))
    ]
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["region"] = obs["region"].astype(str)
    print(f"n cells after filter: {len(obs):,}")

    tab = (obs.groupby(["donor", "region", "cell_type"], observed=True)
              .size().rename("n").reset_index())
    tab["pct"] = tab["n"] / tab.groupby(["donor", "region"], observed=True)["n"].transform("sum") * 100
    tab.to_csv(FIG_DIR / "composition_LA_RA_table.tsv", sep="\t", index=False)

    ct_order = (obs.groupby("cell_type", observed=True).size()
                   .sort_values(ascending=False).index.tolist())
    palette = {ct: CT_PALETTE[i % len(CT_PALETTE)] for i, ct in enumerate(ct_order)}
    pretty = {ct: ct.replace("_", " ") for ct in ct_order}

    # -------- Fig 1: overall LA vs RA --------
    overall = (obs.groupby(["region", "cell_type"], observed=True).size()
                  .unstack("cell_type", fill_value=0).reindex(columns=ct_order)).loc[["LA", "RA"]]
    overall_pct = overall.div(overall.sum(axis=1), axis=0) * 100

    fig = plt.figure(figsize=(COL_DOUBLE, 55 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.35,
                          left=0.06, right=0.80, top=0.80, bottom=0.22)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    _stacked_hbar(ax_a, ct_order, overall_pct, palette, pretty)
    ax_a.set_xlim(0, 100)
    ax_a.set_xlabel("cells (%)")
    ax_a.set_title("composition")
    strip_spines(ax_a); panel_label(ax_a, "a", x=-0.10)

    _stacked_hbar(ax_b, ct_order, overall, palette, pretty)
    ax_b.set_xlabel("cells (n)")
    ax_b.set_title("absolute counts")
    strip_spines(ax_b); ax_b.ticklabel_format(style="plain", axis="x")
    panel_label(ax_b, "b", x=-0.10)

    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[ct]) for ct in ct_order]
    fig.legend(handles, [pretty[ct] for ct in ct_order],
               loc="center left", bbox_to_anchor=(0.81, 0.55),
               frameon=False, title="cell type", title_fontsize=7,
               fontsize=6, labelspacing=0.35, handlelength=1.0)
    fig.text(0.06, 0.94, "LA vs RA cell-type composition", fontsize=8, fontweight="bold")
    fig.text(0.06, 0.04,
             "Heart Cell Atlas, Nuclei only (Harvard-Nuclei + Sanger-Nuclei), Used = Yes. "
             f"n = {len(obs):,} cells from 13 donors.", fontsize=5.5, color="#444")
    fig.savefig(FIG_DIR / "fig_composition_overall.pdf")
    plt.close(fig)

    # -------- Fig 2: per-donor --------
    donors_all = sorted(obs["donor"].unique().tolist())
    has_la = obs.groupby("donor", observed=True)["region"].apply(lambda r: "LA" in set(r))
    has_ra = obs.groupby("donor", observed=True)["region"].apply(lambda r: "RA" in set(r))
    paired = has_la & has_ra
    d_d = [d for d in donors_all if d.startswith("D") and paired[d]]
    h_d = [d for d in donors_all if d.startswith("H") and paired[d]]
    unpaired = [d for d in donors_all if not paired[d]]
    donors = d_d + h_d + unpaired

    wide = (obs.groupby(["donor", "region", "cell_type"], observed=True).size()
               .unstack("cell_type", fill_value=0).reindex(columns=ct_order))
    wide_pct = wide.div(wide.sum(axis=1), axis=0) * 100

    fig = plt.figure(figsize=(COL_DOUBLE, 72 * MM))
    gs = fig.add_gridspec(1, 1, left=0.05, right=0.82, top=0.82, bottom=0.18)
    ax = fig.add_subplot(gs[0, 0])
    group_width = 0.85; bar_w = group_width / 2.0
    x_centers = np.arange(len(donors), dtype=float)
    for i, region in enumerate(["LA", "RA"]):
        xs = x_centers - group_width / 2 + bar_w * (i + 0.5)
        bottoms = np.zeros(len(donors))
        for ct in ct_order:
            vals = np.array([
                wide_pct.loc[(d, region), ct] if (d, region) in wide_pct.index else 0.0
                for d in donors
            ])
            ax.bar(xs, vals, bottom=bottoms, width=bar_w * 0.95,
                   color=palette[ct], edgecolor="white", linewidth=0.25)
            bottoms += vals
    ax.set_xticks(x_centers); ax.set_xticklabels(donors)
    ax.set_ylabel("cells (%)"); ax.set_ylim(0, 100)
    strip_spines(ax)
    for x in x_centers:
        ax.text(x - bar_w / 2, -3.5, "LA", ha="center", va="top", fontsize=5, color="#555")
        ax.text(x + bar_w / 2, -3.5, "RA", ha="center", va="top", fontsize=5, color="#555")
    if unpaired:
        ax.text(0.99, 1.02, f"unpaired (LA-only): {', '.join(unpaired)}",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=5.5, color="#a00", style="italic")
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[ct]) for ct in ct_order]
    fig.legend(handles, [pretty[ct] for ct in ct_order],
               loc="center left", bbox_to_anchor=(0.83, 0.5),
               frameon=False, title="cell type", title_fontsize=7,
               fontsize=6, labelspacing=0.35, handlelength=1.0)
    fig.text(0.05, 0.94, "Per-donor LA / RA composition", fontsize=8, fontweight="bold")
    fig.text(0.05, 0.03,
             "Left bar per donor = LA, right = RA. Heart Cell Atlas, Nuclei only, Used = Yes.",
             fontsize=5.5, color="#444")
    fig.savefig(FIG_DIR / "fig_composition_per_donor.pdf")
    plt.close(fig)

    # -------- Fig 3: top-4 paired --------
    obs4 = obs.copy()
    obs4["top_ct"] = obs4["cell_type"].map(TOP4)
    obs4 = obs4[obs4["top_ct"].notna()]
    g4 = obs4.groupby(["donor", "region", "top_ct"], observed=True).size().rename("n").reset_index()
    g4_pct = g4.copy()
    g4_pct["pct"] = g4_pct["n"] / g4_pct.groupby(["donor", "region"], observed=True)["n"].transform("sum") * 100

    top4_order = ["atrial_cardiomyocyte", "fibroblast", "endothelial_or_vascular", "immune"]
    top4_pretty = {
        "atrial_cardiomyocyte":    "atrial CM",
        "fibroblast":              "fibroblast",
        "endothelial_or_vascular": "endo / vascular",
        "immune":                  "immune",
    }
    paired_donors = [d for d in donors if d in paired.index and paired[d]]

    fig = plt.figure(figsize=(COL_DOUBLE, 95 * MM))
    gs = fig.add_gridspec(2, 1, hspace=0.35, left=0.07, right=0.82, top=0.88, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0], sharex=ax_a)
    bar_w = 0.38
    x_centers = np.arange(len(paired_donors), dtype=float)

    def _draw(ax, source_df, value_col, ylabel, put_legend=False):
        for i, region in enumerate(["LA", "RA"]):
            xs = x_centers + (i - 0.5) * bar_w
            bottoms = np.zeros(len(paired_donors))
            for ct in top4_order:
                vals = []
                for d in paired_donors:
                    v = source_df[(source_df["donor"] == d) & (source_df["region"] == region)
                                  & (source_df["top_ct"] == ct)][value_col]
                    vals.append(float(v.iloc[0]) if len(v) else 0.0)
                vals = np.array(vals)
                ax.bar(xs, vals, bottom=bottoms, width=bar_w * 0.92,
                       color=TOP4_PALETTE[ct], edgecolor="white", linewidth=0.3,
                       label=top4_pretty[ct] if (region == "LA" and put_legend) else None)
                bottoms += vals
        ax.set_ylabel(ylabel); strip_spines(ax)

    _draw(ax_a, g4, "n", "cells (n)", put_legend=True)
    panel_label(ax_a, "a", x=-0.05)
    ax_a.legend(loc="upper right", frameon=False, ncol=4, fontsize=6,
                bbox_to_anchor=(1.0, 1.08), handlelength=1.0, columnspacing=1.0)

    _draw(ax_b, g4_pct, "pct", "cells (%)")
    ax_b.set_ylim(0, 100); ax_b.set_xticks(x_centers); ax_b.set_xticklabels(paired_donors)
    panel_label(ax_b, "b", x=-0.05)
    for x in x_centers:
        ax_b.text(x - bar_w / 2, -5.5, "LA", ha="center", va="top", fontsize=5, color="#555")
        ax_b.text(x + bar_w / 2, -5.5, "RA", ha="center", va="top", fontsize=5, color="#555")

    fig.text(0.07, 0.955, "Top-4 cell classes: LA (left) vs RA (right) per donor", fontsize=8, fontweight="bold")
    fig.text(0.07, 0.03,
             "Mapping: atrial CM ← Atrial_Cardiomyocyte;  fibroblast ← Fibroblast;  "
             "endo/vascular ← Endothelial;  immune ← Myeloid + Lymphoid. "
             "Nuclei only, Used = Yes, n_cells ≥ 30 per (donor × chamber × class).",
             fontsize=5.5, color="#444")
    fig.savefig(FIG_DIR / "fig_composition_top4_paired.pdf")
    plt.close(fig)

    print("done.")


if __name__ == "__main__":
    main()
