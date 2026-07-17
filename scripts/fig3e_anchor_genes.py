#!/usr/bin/env python3
"""Figure 3E  Anchor genes leading the strongest Fig 3D shifts."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

PB_BASE   = ROOT / "Project_P3" / "GSE255612" / "pseudobulk_full"
LE_TSV    = ROOT / "results" / "fig3_primary_program_le.tsv"
DELTA_TSV = ROOT / "results" / "fig3_primary_program_celltype_delta.tsv"
FIG_OUT   = ROOT / "figures" / "fig3e_anchor_genes.pdf"

LA_COLOR  = "#D6604D"
RA_COLOR  = "#4393C3"
BRBG_BROWN = "#8C510A"
BRBG_TEAL  = "#01665E"
GREY_LIGHT = "#E6E6E6"
GREY_MID  = "#9C9C9C"
GREY_DARK = "#444444"

TARGETS = [
    {"program_id": "REACTOME__REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM",
     "title": "F10 GAG metabolism × LA fibroblast",
     "compartment": "fibroblast",
     "subtitle": "RA-biased program (LA gain)  ·  family TGF-β"},
    {"program_id": "REACTOME__REACTOME_SIGNALING_BY_BMP",
     "title": "F11 BMP signalling × LA adipocyte",
     "compartment": "adipocyte",
     "subtitle": "RA-biased program (LA gain)  ·  family BMP"},
    {"program_id": "REACTOME__REACTOME_SIGNALING_BY_BMP",
     "title": "F11 BMP signalling × LA lymphatic endothelial",
     "compartment": "lymphatic_endothelial",
     "subtitle": "RA-biased program (LA loss)  ·  family BMP"},
]

N_TOP = 10


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "axes.titlesize":    9,
        "axes.labelsize":    8,
        "xtick.labelsize":   7,
        "ytick.labelsize":   7,
        "axes.linewidth":    0.6,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
        "axes.spines.top":   False,
        "axes.spines.right": False,
    })


def load_pseudobulk(compartment: str):
    """Returns (expr_df, samples_df).  expr indexed by gene symbol when
    possible, otherwise the original index (Ensembl ID)."""
    gdir = PB_BASE / compartment
    expr = pd.read_csv(gdir / "expr.tsv", sep="\t", index_col=0)
    samps = pd.read_csv(gdir / "samples.tsv", sep="\t")
    return expr, samps


def main() -> None:
    set_publication_style()

    # Load LE genes per program (NCP_v2)
    le_df = pd.read_csv(LE_TSV, sep="\t")
    le_lookup = {row["program_id"]:
                    [g for g in str(row["leading_edge_genes"]).split(",")
                       if g and g != "nan"]
                  for _, row in le_df.iterrows()}
    delta = pd.read_csv(DELTA_TSV, sep="\t")
    delta_idx = delta.set_index(["program_id", "cell_type"])

    n_rows = len(TARGETS)
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(180 * MM, 60 * MM * n_rows))
    grid = fig.add_gridspec(
        nrows=n_rows, ncols=1, hspace=0.55,
        left=0.20, right=0.96, top=0.93, bottom=0.10)

    fig.suptitle(
        "Fig 3E  Anchor genes leading the strongest single-cell rewiring "
        "shifts in LA snRNA  (per program × cell type)",
        fontsize=8.8, fontweight="bold", y=0.985, x=0.02, ha="left")

    for idx, tgt in enumerate(TARGETS):
        ax = fig.add_subplot(grid[idx])
        pid = tgt["program_id"]
        ct  = tgt["compartment"]
        expr, samps = load_pseudobulk(ct)

        # align expr columns to samps rows
        expr = expr[samps["donor_id"]]
        af_mask  = (samps["rhythm"] == "AF").to_numpy()
        ctl_mask = (samps["rhythm"] == "Ctrl").to_numpy()

        le = le_lookup.get(pid, [])
        if not le:
            ax.text(0.5, 0.5, f"NO LE genes for {pid}",
                     ha="center", va="center",
                     transform=ax.transAxes)
            continue
        # restrict to LE genes present in the pseudobulk
        present = [g for g in le if g in expr.index]
        if not present:
            ax.text(0.5, 0.5,
                     f"No LE genes present in {ct} pseudobulk",
                     ha="center", va="center",
                     transform=ax.transAxes)
            continue
        sub_expr = expr.loc[present]
        af_mean   = sub_expr.iloc[:, af_mask].mean(axis=1)
        ctl_mean  = sub_expr.iloc[:, ctl_mask].mean(axis=1)
        delta_g   = af_mean - ctl_mean

        # SE per gene — use across-donor SD / sqrt(n)
        af_sd  = sub_expr.iloc[:, af_mask].std(axis=1, ddof=1)
        ctl_sd = sub_expr.iloc[:, ctl_mask].std(axis=1, ddof=1)
        n_af  = int(af_mask.sum()); n_ctl = int(ctl_mask.sum())
        se_g  = np.sqrt(af_sd**2 / n_af + ctl_sd**2 / n_ctl)

        # rank by direction-aligned shift  (sign matches program Δ)
        prog_dscore = float(delta_idx.loc[(pid, ct), "dScore"])
        sign = np.sign(prog_dscore)
        score_g = delta_g * sign
        ranked = (pd.DataFrame({"gene": present, "delta": delta_g.values,
                                  "se": se_g.values,
                                  "score": score_g.values})
                      .sort_values("score", ascending=False)
                      .head(N_TOP)
                      .reset_index(drop=True))

        # Bar plot horizontal
        ranked = ranked.iloc[::-1].reset_index(drop=True)  # so largest at top
        y = np.arange(len(ranked))
        colors = [BRBG_BROWN if v >= 0 else BRBG_TEAL for v in ranked["delta"]]
        ax.barh(y, ranked["delta"], color=colors,
                  edgecolor="black", linewidth=0.4, height=0.65, zorder=3)
        for yi, val, se in zip(y, ranked["delta"], ranked["se"]):
            ax.plot([val - 1.96 * se, val + 1.96 * se], [yi, yi],
                     color="black", lw=0.6, zorder=4)

        ax.axvline(0, color=GREY_DARK, lw=0.6, zorder=1)
        ax.set_yticks(y)
        ax.set_yticklabels(ranked["gene"].tolist(), fontsize=7)

        ax.set_xlabel("Δ mean expression  (AF − Ctrl)   ·   "
                        "log-norm pseudobulk", fontsize=7)
        ax.set_title(
            f"{tgt['title']}   "
            f"·   program-level Δ = {prog_dscore:+.4f}   "
            f"·   {tgt['subtitle']}",
            loc="left", fontsize=7.5, fontweight="bold",
            color=RA_COLOR, pad=4)
        ax.tick_params(axis="x", pad=2)

    fig.text(0.02, 0.020,
              "Bars: per-gene AF − Ctrl mean of log-normalized snRNA "
              "pseudobulk in the named LA cell type, restricted to that "
              "program's leading-edge genes;  whiskers = 95 % CI from per-donor SE.\n"
              "Brown = AF higher,  teal = AF lower  (BrBG palette).  "
              "Top 10 genes ranked by direction-aligned shift (signed by program-level Δ).",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")

    # Console summary
    print("\n=== Anchor genes per target ===")
    for tgt in TARGETS:
        pid = tgt["program_id"]; ct = tgt["compartment"]
        le = le_lookup.get(pid, [])
        expr, samps = load_pseudobulk(ct)
        expr = expr[samps["donor_id"]]
        af_mask  = (samps["rhythm"] == "AF").to_numpy()
        ctl_mask = (samps["rhythm"] == "Ctrl").to_numpy()
        present = [g for g in le if g in expr.index]
        if not present: continue
        sub_expr = expr.loc[present]
        delta_g  = sub_expr.iloc[:, af_mask].mean(axis=1) - sub_expr.iloc[:, ctl_mask].mean(axis=1)
        prog_dscore = float(delta_idx.loc[(pid, ct), "dScore"])
        sign = np.sign(prog_dscore)
        ranked = pd.DataFrame({"gene": present, "delta": delta_g.values,
                                "score": delta_g.values * sign}) \
                      .sort_values("score", ascending=False).head(N_TOP)
        print(f"\n  {tgt['title']}  (n_LE={len(le)}, present={len(present)})")
        for _, r in ranked.iterrows():
            print(f"     {r['gene']:12s}  Δ = {r['delta']:+.4f}")


if __name__ == "__main__":
    main()
