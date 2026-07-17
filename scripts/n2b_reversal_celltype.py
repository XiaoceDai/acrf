#!/usr/bin/env python3
"""
N2b — Map the 24 chamber-identity reversal genes (Fig. 7, results/
n2_chamber_reversal_genes.tsv) onto the four single-nucleus
compartments of GSE255612 to identify the cellular carrier(s) of the
matrix-invasion signature.

Per gene, per compartment:
    delta_expr = mean(expression in AF biosamples) - mean(expression in
                  Ctrl biosamples)
For RA->LA reversal genes (n = 21; matrix-invasion candidates), a
positive delta_expr in a compartment indicates that the gene is
up-regulated in AF in that compartment (i.e., that compartment is a
candidate carrier of the reversal).  For LA->RA reversal genes (n = 3),
a negative delta_expr indicates the same.  The sign-aligned delta is
reported as `delta_aligned` so that positive = consistent with the
reversal direction in either case.

Outputs:
    results/n2b_reversal_per_celltype_long.tsv   (gene x compartment rows)
    results/n2b_reversal_per_celltype_summary.tsv
    figures/figN2b_reversal_celltype.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (
    COL_DOUBLE, LA_COLOR, RA_COLOR, panel_label,
    set_nature_style, strip_spines,
)

set_nature_style()

REVERSAL_TSV = ROOT / "results" / "n2_chamber_reversal_genes.tsv"

COMPARTMENTS = [
    "atrialCM",
    "fibroblast",
    "endothelial_or_vascular",
    "immune",
]

COMPARTMENT_COLORS = {
    "atrialCM":                "#E6AB02",
    "fibroblast":              "#66A61E",
    "endothelial_or_vascular": "#7570B3",
    "immune":                  "#D95F02",
}

OUT_LONG    = ROOT / "results" / "n2b_reversal_per_celltype_long.tsv"
OUT_SUMMARY = ROOT / "results" / "n2b_reversal_per_celltype_summary.tsv"
FIG_OUT     = ROOT / "figures" / "figN2b_reversal_celltype.pdf"


def main() -> None:
    rev = pd.read_csv(REVERSAL_TSV, sep="\t")
    # sign convention for "aligned" delta:
    #   RA-biased reversal (RA->LA): expected delta_expr POSITIVE in carrier
    #   LA-biased reversal (LA->RA): expected delta_expr NEGATIVE in carrier
    rev["expected_dir_in_AF"] = np.where(
        rev["healthy_direction"] == "RA-biased", +1, -1
    )
    print(f"Reversal genes loaded: {len(rev)} "
          f"(RA->LA = {(rev['healthy_direction']=='RA-biased').sum()}, "
          f"LA->RA = {(rev['healthy_direction']=='LA-biased').sum()})")

    long_rows = []
    summary_rows = []
    per_compartment_aligned: dict[str, list[float]] = {c: [] for c in COMPARTMENTS}

    for compartment in COMPARTMENTS:
        expr_path = ROOT / "Project_P3" / "GSE255612" / "pseudobulk" \
                    / compartment / "expr.tsv"
        meta_path = ROOT / "Project_P3" / "GSE255612" / "pseudobulk" \
                    / compartment / "samples.tsv"
        if not expr_path.exists():
            print(f"  WARN: missing {expr_path}")
            continue
        expr = pd.read_csv(expr_path, sep="\t", index_col=0)
        meta = pd.read_csv(meta_path, sep="\t")
        # expr columns use bare donor_id (e.g. "P1279_3n") while meta sample_id
        # carries the compartment prefix (e.g. "fibroblast__P1279_3n").
        af_samples = meta.loc[meta["rhythm"] == "AF",   "donor_id"].tolist()
        ct_samples = meta.loc[meta["rhythm"] == "Ctrl", "donor_id"].tolist()
        af_present = [s for s in af_samples if s in expr.columns]
        ct_present = [s for s in ct_samples if s in expr.columns]
        print(f"  {compartment}: AF n = {len(af_present)}, Ctrl n = {len(ct_present)}, "
              f"genes = {expr.shape[0]}")

        for _, r in rev.iterrows():
            g = r["gene"]
            if g not in expr.index:
                # gene not measured in this compartment
                long_rows.append({
                    "gene": g, "compartment": compartment,
                    "healthy_direction": r["healthy_direction"],
                    "expected_dir_in_AF": int(r["expected_dir_in_AF"]),
                    "n_AF": len(af_present), "n_Ctrl": len(ct_present),
                    "mean_expr_AF": np.nan, "mean_expr_Ctrl": np.nan,
                    "delta_expr": np.nan, "delta_aligned": np.nan,
                    "fold_AF_vs_Ctrl": np.nan,
                })
                continue
            af_v = expr.loc[g, af_present].dropna()
            ct_v = expr.loc[g, ct_present].dropna()
            if len(af_v) < 3 or len(ct_v) < 3:
                continue
            af_mean = float(af_v.mean()); ct_mean = float(ct_v.mean())
            d = af_mean - ct_mean
            d_aligned = d * r["expected_dir_in_AF"]
            per_compartment_aligned[compartment].append(d_aligned)
            long_rows.append({
                "gene": g, "compartment": compartment,
                "healthy_direction": r["healthy_direction"],
                "expected_dir_in_AF": int(r["expected_dir_in_AF"]),
                "n_AF": len(af_present), "n_Ctrl": len(ct_present),
                "mean_expr_AF": round(af_mean, 5),
                "mean_expr_Ctrl": round(ct_mean, 5),
                "delta_expr": round(d, 5),
                "delta_aligned": round(d_aligned, 5),
                "fold_AF_vs_Ctrl": round((af_mean / ct_mean) if ct_mean > 0 else np.nan, 4),
            })

    long = pd.DataFrame(long_rows)
    OUT_LONG.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(OUT_LONG, sep="\t", index=False)
    print(f"wrote {OUT_LONG.relative_to(ROOT)}")

    for compartment in COMPARTMENTS:
        sub = long[long["compartment"] == compartment].dropna(subset=["delta_aligned"])
        if len(sub) == 0:
            continue
        # split by direction
        ra2la = sub[sub["healthy_direction"] == "RA-biased"]
        la2ra = sub[sub["healthy_direction"] == "LA-biased"]
        # sign-aligned mean: positive = expression goes in the expected
        # reversal direction in this compartment
        summary_rows.append({
            "compartment":             compartment,
            "n_genes_evaluated":       len(sub),
            "mean_delta_aligned":      round(sub["delta_aligned"].mean(), 5),
            "median_delta_aligned":    round(sub["delta_aligned"].median(), 5),
            "frac_genes_aligned_pos":  round((sub["delta_aligned"] > 0).mean(), 3),
            "n_genes_ra2la":           len(ra2la),
            "mean_delta_aligned_ra2la": round(ra2la["delta_aligned"].mean(), 5)
                if len(ra2la) else np.nan,
            "frac_aligned_pos_ra2la":  round((ra2la["delta_aligned"] > 0).mean(), 3)
                if len(ra2la) else np.nan,
            "n_genes_la2ra":           len(la2ra),
            "mean_delta_aligned_la2ra": round(la2ra["delta_aligned"].mean(), 5)
                if len(la2ra) else np.nan,
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)}")

    print("\n=== compartment-level summary ===")
    print(summary.to_string(index=False))

    fig = plt.figure(figsize=(COL_DOUBLE, 0.55 * COL_DOUBLE))
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.40,
                          height_ratios=[1.0, 1.1])

    # panel a: per-compartment mean delta_aligned (bar)
    ax = fig.add_subplot(gs[0, 0])
    cs = list(summary["compartment"])
    means = list(summary["mean_delta_aligned"])
    fracs = list(summary["frac_genes_aligned_pos"])
    colors = [COMPARTMENT_COLORS.get(c, "#888888") for c in cs]
    x = np.arange(len(cs))
    bars = ax.bar(x, means, color=colors, edgecolor="white", linewidth=0.4)
    for xi, m, f in zip(x, means, fracs):
        ax.text(xi, m + (0.0005 if m >= 0 else -0.0005),
                f"{m:+.4f}\n({int(f*100)}% +)",
                ha="center", va="bottom" if m >= 0 else "top",
                fontsize=5, fontweight="bold")
    ax.axhline(0, color="#555", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(cs, rotation=20, ha="right")
    ax.set_ylabel("mean Δ aligned (AF − Ctrl)\n(positive = consistent with reversal)")
    ax.set_title("a  Carrier compartment of the chamber-reversal signature",
                 loc="left", pad=4, fontsize=6.6)
    strip_spines(ax)
    panel_label(ax, "a", x=-0.18, y=1.04)

    # panel b: split by healthy_direction
    ax = fig.add_subplot(gs[0, 1])
    width = 0.35
    ra_means = [s["mean_delta_aligned_ra2la"] for s in summary_rows]
    la_means = [s["mean_delta_aligned_la2ra"] for s in summary_rows]
    ax.bar(x - width/2, ra_means, width=width,
           color=RA_COLOR, edgecolor="white", linewidth=0.4,
           label=f"RA→LA (n=21 matrix)")
    ax.bar(x + width/2, la_means, width=width,
           color=LA_COLOR, edgecolor="white", linewidth=0.4,
           label=f"LA→RA (n=3)")
    ax.axhline(0, color="#555", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(cs, rotation=20, ha="right")
    ax.set_ylabel("mean Δ aligned (AF − Ctrl)")
    ax.set_title("b  Direction-specific carrier signal",
                 loc="left", pad=4, fontsize=6.6)
    ax.legend(loc="upper left", frameon=False, fontsize=5)
    strip_spines(ax)
    panel_label(ax, "b", x=-0.18, y=1.04)

    # panel c-bottom: per-gene heatmap (genes x compartments)
    ax = fig.add_subplot(gs[1, :])
    pivot = (long.dropna(subset=["delta_aligned"])
                  .pivot_table(index="gene", columns="compartment",
                               values="delta_aligned"))
    # order columns
    pivot = pivot[[c for c in COMPARTMENTS if c in pivot.columns]]
    # sort genes: RA-biased first (largest n=21), then LA-biased; within each, by overall fibro carrier signal
    gene_dir = (rev.set_index("gene")["healthy_direction"]
                  .reindex(pivot.index).fillna("RA-biased"))
    # sort by direction then by max across compartments
    pivot["_dir"] = gene_dir
    pivot["_rank"] = pivot[[c for c in COMPARTMENTS if c in pivot.columns]].max(axis=1)
    pivot = pivot.sort_values(["_dir", "_rank"], ascending=[True, False])
    pivot = pivot.drop(columns=["_dir", "_rank"])

    vmax = float(np.nanmax(np.abs(pivot.values))) or 0.001
    im = ax.imshow(pivot.values, cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(np.arange(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=5)
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=20, ha="right")
    # mark direction bar on left
    for i, g in enumerate(pivot.index):
        d = gene_dir.get(g, "RA-biased")
        bar_color = RA_COLOR if d == "RA-biased" else LA_COLOR
        ax.add_patch(mpatches.Rectangle(
            (-0.65, i - 0.5), 0.25, 1,
            facecolor=bar_color, edgecolor="none", clip_on=False,
        ))
    ax.set_xlim(-0.7, pivot.shape[1] - 0.5)
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Δ aligned (AF − Ctrl)", fontsize=5)
    cbar.ax.tick_params(labelsize=5)
    ax.set_title(
        "c  Per-gene reversal signal across compartments "
        "(left bar: red = RA→LA matrix-invasion gene; blue = LA→RA)",
        loc="left", pad=4, fontsize=6.4,
    )
    panel_label(ax, "c", x=-0.06, y=1.04)

    fig.suptitle(
        "Chamber-reversal genes localise to the fibroblast compartment "
        "in single-nucleus AF atrial tissue (GSE255612)",
        y=1.00, fontsize=7.2, fontweight="bold",
    )

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
