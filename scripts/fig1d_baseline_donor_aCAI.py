#!/usr/bin/env python3
"""
Figure 1d (REBUILD, baseline-only) — Donor-level aCAI primitive on the
healthy reference cohort (GSE112339, 8 paired LA / RA donors).

For each NCP_v2 program (Fig. 1c, n = 30) and each healthy donor d:
    aCAI_d(p) = mean(log2-CPM of program p member genes in LA donor d)
              - mean(log2-CPM of program p member genes in RA donor d)
A positive aCAI_d(p) means program p is LA-biased in donor d at the
healthy baseline; a negative value means RA-biased.

The figure shows a heatmap (30 programs x 8 donors) of aCAI_d(p), with
programs grouped by their NCP_v2 chamber direction.  The right strip
reports each program's mean aCAI across donors and its sign-concordance
across the 8 donors.

This is the framework's atomic-unit panel: it shows what the donor-level
aCAI metric looks like at the healthy baseline before any AF / SR
information enters the analysis.

Inputs:
    results/fig1c_ncp_panel_v2.tsv
    bulk/Project_P01/GSE112339/GSE112339_HTSeq_hg19_Ensembl.txt.gz
    bulk/Project_P01/gene_sets/h.all.v2026.1.Hs.symbols.gmt
    bulk/Project_P01/gene_sets/c2.cp.reactome.v2026.1.Hs.symbols.gmt
Output:
    figures/fig1d_baseline_donor_aCAI.pdf
    results/fig1d_baseline_donor_aCAI.tsv  (program x donor aCAI matrix)
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

# Lock Arial first so glyphs like ⁺ / β are taken from Arial (which has them)
# rather than falling back to Liberation Sans (which doesn't).
mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
})

NCP_TSV   = ROOT / "results" / "fig1c_ncp_panel_v2.tsv"
FAM_TSV   = ROOT / "results" / "fig1e_baseline_axis_families.tsv"
COUNTS    = ROOT / "bulk" / "Project_P01" / "GSE112339" \
                / "GSE112339_HTSeq_hg19_Ensembl.txt.gz"
GMT_HALL  = ROOT / "bulk" / "Project_P01" / "gene_sets" \
                / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC  = ROOT / "bulk" / "Project_P01" / "gene_sets" \
                / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"

OUT_TSV   = ROOT / "results" / "fig1d_baseline_donor_aCAI.tsv"
FIG_OUT   = ROOT / "figures" / "fig1d_baseline_donor_aCAI.pdf"


def load_gmt(path: Path, prefix: str) -> dict[str, list[str]]:
    """Read a MSigDB-style GMT and return {prefix__name: [genes]}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = f"{prefix}__{parts[0]}"
            genes = [g for g in parts[2:] if g]
            out[name] = genes
    return out


def load_log_cpm() -> pd.DataFrame:
    raw = pd.read_csv(COUNTS, sep="\t", index_col=0)
    raw = raw.groupby(raw.index).max()
    keep = [c for c in raw.columns if c.startswith(("LA_", "RA_"))]
    raw = raw[keep]
    cpm = raw.div(raw.sum(axis=0), axis=1) * 1e6
    return np.log2(cpm + 1.0)


def align_to_canonical_order(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep d/e/f row order linked:
    if fig1e family table exists, follow its pathway order exactly.
    Otherwise fall back to stable direction-aware ordering.
    """
    if FAM_TSV.exists():
        fam = pd.read_csv(FAM_TSV, sep="\t")
        order = fam["pathway"].tolist()
        order_map = {p: i for i, p in enumerate(order)}
        tmp = df.copy()
        tmp["sort_idx"] = tmp["pathway"].map(order_map)
        tmp = tmp.sort_values(["sort_idx", "program_order"], na_position="last")
        tmp = tmp.drop(columns=["sort_idx"]).reset_index(drop=True)
        return tmp

    # Fallback (first build): LA on top, RA below.
    return pd.concat([
        df[df["direction"] == "LA-biased"].sort_values("mean_aCAI", ascending=False),
        df[df["direction"] == "RA-biased"].sort_values("mean_aCAI", ascending=True),
    ]).reset_index(drop=True)


def main() -> None:
    print("Loading NCP_v2 panel and reference data ...")
    ncp = pd.read_csv(NCP_TSV, sep="\t")
    print(f"  NCP_v2 panel: {len(ncp)} programs "
          f"({(ncp['direction']=='LA-biased').sum()} LA, "
          f"{(ncp['direction']=='RA-biased').sum()} RA)")

    pathways = {**load_gmt(GMT_HALL, "HALLMARK"),
                **load_gmt(GMT_REAC, "REACTOME")}
    log_cpm = load_log_cpm()
    print(f"  log-CPM: {log_cpm.shape[0]} genes x {log_cpm.shape[1]} samples")

    la_cols = sorted([c for c in log_cpm.columns if c.startswith("LA_")])
    ra_cols = sorted([c for c in log_cpm.columns if c.startswith("RA_")])
    donors = [c[3:] for c in la_cols]
    print(f"  paired donors: {len(donors)}  ({', '.join(donors)})")

    # Compute per-donor per-program aCAI = mean LA expr - mean RA expr
    rows = []
    for _, r in ncp.iterrows():
        gset = pathways.get(r["pathway"], [])
        gset = [g for g in gset if g in log_cpm.index]
        if len(gset) < 5:
            continue
        sub = log_cpm.loc[gset]
        per_donor = []
        for d in donors:
            la_mean = sub[f"LA_{d}"].mean()
            ra_mean = sub[f"RA_{d}"].mean()
            per_donor.append(la_mean - ra_mean)
        rows.append({
            "program_order":  int(r["program_order"]),
            "pathway":        r["pathway"],
            "direction":      r["direction"],
            "n_genes_used":   len(gset),
            **{f"aCAI_{d}": v for d, v in zip(donors, per_donor)},
        })
    aCAI_df = pd.DataFrame(rows)

    donor_cols = [f"aCAI_{d}" for d in donors]
    aCAI_df["mean_aCAI"] = aCAI_df[donor_cols].mean(axis=1)
    aCAI_df["sign_concordance"] = (
        np.sign(aCAI_df[donor_cols]).eq(np.sign(aCAI_df["mean_aCAI"]),
                                          axis=0).sum(axis=1).astype(int)
    )
    aCAI_df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")

    # ---- Figure -----------------------------------------------------------
    sorted_df = align_to_canonical_order(aCAI_df)

    n_la = (sorted_df["direction"] == "LA-biased").sum()
    n_ra = (sorted_df["direction"] == "RA-biased").sum()
    n_8of8 = int((sorted_df["sign_concordance"] == 8).sum())

    mat = sorted_df[donor_cols].values
    vmax = float(np.nanmax(np.abs(mat))) or 0.001
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Wider figure: extra column on the far right for inline pathway labels
    # so the bar-chart strip is self-contained.
    fig = plt.figure(figsize=(COL_DOUBLE * 1.30, COL_DOUBLE * 0.92))
    gs = fig.add_gridspec(
        1, 4, width_ratios=[0.32, 6.2, 1.05, 1.85], wspace=0.04,
        top=0.90, bottom=0.13, left=0.04, right=0.985,
    )

    # ---- col 0 : LA / RA direction band with rotated text label ----
    # Updated convention: RA = blue (#4393C3), LA = red (#D6604D)
    LA_COLOR = "#D6604D"
    RA_COLOR = "#4393C3"
    ax_dir = fig.add_subplot(gs[0, 0])
    for i, d in enumerate(sorted_df["direction"]):
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
    ax_dir.set_xlim(0, 1); ax_dir.set_ylim(len(sorted_df) - 0.5, -0.5)
    ax_dir.set_xticks([]); ax_dir.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        ax_dir.spines[s].set_visible(False)

    # ---- col 1 : main heatmap (no in-cell numbers; clean) ----
    # Diverging colormap aligned with the new convention:
    #   positive aCAI -> LA-biased -> red
    #   negative aCAI -> RA-biased -> blue
    from matplotlib.colors import LinearSegmentedColormap
    acai_cmap = LinearSegmentedColormap.from_list(
        "LA_RA",
        [(0.0, RA_COLOR), (0.5, "#FFFFFF"), (1.0, LA_COLOR)],
    )

    ax = fig.add_subplot(gs[0, 1], sharey=ax_dir)
    im = ax.imshow(mat, cmap=acai_cmap, norm=norm, aspect="auto",
                    extent=[-0.5, len(donors) - 0.5,
                             len(sorted_df) - 0.5, -0.5])

    ax.set_yticks(np.arange(len(sorted_df)))
    ax.set_yticklabels(
        [short_pathway_label(r["pathway"], 36)
         for _, r in sorted_df.iterrows()],
        fontsize=8.1,
    )
    ax.set_xticks(np.arange(len(donors)))
    ax.set_xticklabels(donors, rotation=45, ha="right", fontsize=8.2)
    ax.set_xlabel("donor", fontsize=9.2)
    # Horizontal separator between LA and RA blocks
    ax.axhline(n_la - 0.5, color="black", lw=0.9, zorder=4)

    ax.set_title("d  Baseline donor aCAI (shared row order)", loc="left", pad=8,
                 fontsize=10.0, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, fraction=0.020, pad=0.012)
    cbar.set_label("aCAI = log₂(LA) − log₂(RA)", fontsize=8.8)
    cbar.ax.tick_params(labelsize=7.8)
    # Direction tags on the colorbar to make sign meaning explicit
    cbar.ax.text(
        0.5, 1.04, "LA-biased ↑",
        transform=cbar.ax.transAxes, ha="center", va="bottom",
        fontsize=7.8, color=LA_COLOR, fontweight="bold",
    )
    cbar.ax.text(
        0.5, -0.04, "RA-biased ↓",
        transform=cbar.ax.transAxes, ha="center", va="top",
        fontsize=7.8, color=RA_COLOR, fontweight="bold",
    )

    # ---- col 2 : mean aCAI bar (one summary column instead of two) ----
    ax2 = fig.add_subplot(gs[0, 2], sharey=ax_dir)
    colors = [LA_COLOR if d == "LA-biased" else RA_COLOR
              for d in sorted_df["direction"]]
    y = np.arange(len(sorted_df))
    ax2.barh(y, sorted_df["mean_aCAI"], color=colors,
             edgecolor="white", linewidth=0.3)
    ax2.axvline(0, color="#444", lw=0.4)
    ax2.set_xlabel("mean aCAI", fontsize=8.4)
    ax2.tick_params(labelsize=7.8)
    ax2.set_yticks([])
    ax2.axhline(n_la - 0.5, color="black", lw=0.7)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    # NOTE: ax2 shares y with ax_dir, which is already inverted (row 0 on top).
    # Do NOT call invert_yaxis() here — it would flip the bars relative to the
    # heatmap and the direction band.
    # Annotate sign concordance on each bar (8/8 = green tick)
    xpad = 0.02
    for yi, v, sc in zip(y, sorted_df["mean_aCAI"],
                          sorted_df["sign_concordance"]):
        sc_color = "#1B7837" if sc == 8 else "#666"
        sc_weight = "bold" if sc == 8 else "normal"
        ax2.text(v + (xpad if v > 0 else -xpad), yi,
                  f"{int(sc)}/8",
                  fontsize=7.2, ha="left" if v > 0 else "right",
                  va="center", color=sc_color, fontweight=sc_weight)
    # ---- col 3 : inline pathway labels next to each bar ------------------
    ax3 = fig.add_subplot(gs[0, 3], sharey=ax_dir)
    for yi, (_, r) in enumerate(sorted_df.iterrows()):
        col = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        ax3.text(
            0.0, yi, short_pathway_label(r["pathway"], 36),
            ha="left", va="center", fontsize=8.1, color=col,
        )
    ax3.axhline(n_la - 0.5, color="black", lw=0.7)
    ax3.set_xlim(0, 1)
    ax3.set_xticks([])
    ax3.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        ax3.spines[s].set_visible(False)

    # Detailed method notes are moved to figures/fig1_legend.txt.

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")

    # ---- console summary --------------------------------------------------
    print("\n=== NCP_v2 baseline aCAI summary ===")
    print(f"sign concordance distribution (out of 8 donors):")
    print(sorted_df["sign_concordance"].value_counts().sort_index().to_string())
    print(f"\nprograms with 8/8 donor sign concordance: "
          f"{(sorted_df['sign_concordance']==8).sum()} / {len(sorted_df)}")


if __name__ == "__main__":
    main()
