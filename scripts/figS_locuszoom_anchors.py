#!/usr/bin/env python3
"""
Figure S  LocusZoom-style mini-Manhattan for the 6 cis-eQTL × AF-GWAS
            anchors (4 framework + 2 positive controls).

For each anchor we render a 2-track stacked plot over the ±1 Mb cis-window:
    upper track  : −log10(eQTL p-value) for the anchor gene
    lower track  : −log10(AF GWAS p-value)
    x-axis       : genomic position (Mb)
The anchor gene's TSS / body is annotated with a vertical band; the lead
eQTL and lead GWAS SNP positions are highlighted with vertical lines.
PP3 / PP4 from the original coloc.abf are reported in panel.

Read this figure as visual evidence for the high-PP3 / low-PP4 pattern
in Fig 6E:  if the eQTL and GWAS Manhattan peaks fall at *different*
genomic positions, the two traits have *different* causal variants —
which is exactly what coloc.abf reports as PP3 ≫ PP4.

Inputs:
    Project_P6/cache/eqtl_<anchor>_<window>.tsv.gz
    Project_P6/cache/gwas_<anchor>_<window>.tsv.gz
    results/figure6_mr_anchors.tsv             (PP-values)
    results/figure6_mr_positive_controls.tsv

Output:
    figures/figS_locuszoom_anchors.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "Project_P6" / "cache"
FIG_OUT = ROOT / "figures" / "figS_locuszoom_anchors.pdf"

LA_COLOR  = "#C97A74"
RA_COLOR  = "#7D99B8"
HL_GREEN  = "#3C8C57"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"

EQTL_COLOR = "#7E5BBF"
GWAS_COLOR = "#C25450"

EQTL_NAMES = ["molecular_trait_id", "chromosome", "position", "ref", "alt",
               "variant", "ma_samples", "maf", "pvalue", "beta", "se",
               "type", "ac", "an", "r2", "molecular_trait_object_id",
               "gene_id", "median_tpm", "rsid", "nlog10p"]
GWAS_NAMES = ["chromosome", "base_pair_location", "effect_allele",
               "other_allele", "beta", "standard_error",
               "effect_allele_frequency", "p_value", "rsid", "rs_id",
               "n", "hm_coordinate_conversion", "hm_code", "variant_id"]

# Anchor metadata (must match h4 / h4c definitions exactly).
ANCHORS = [
    {"symbol": "PRELP",  "axis": "matrix",            "group": "anchor",
     "gene_id": "ENSG00000188783", "chrom": 1,
     "tss": 203475828, "tes": 203491352,
     "win": (202475828, 204491352)},
    {"symbol": "OMD",    "axis": "matrix",            "group": "anchor",
     "gene_id": "ENSG00000127083", "chrom": 9,
     "tss": 92414245, "tes": 92424461,
     "win": (91414245, 93424461)},
    {"symbol": "IDH2",   "axis": "metabolic",         "group": "anchor",
     "gene_id": "ENSG00000182054", "chrom": 15,
     "tss": 90083045, "tes": 90102504,
     "win": (89083045, 91102504)},
    {"symbol": "NDUFA8", "axis": "metabolic",         "group": "anchor",
     "gene_id": "ENSG00000119421", "chrom": 9,
     "tss": 122144058, "tes": 122159819,
     "win": (121144058, 123159819)},
    {"symbol": "PITX2",  "axis": "positive_control",  "group": "pos_ctrl",
     "gene_id": "ENSG00000164093", "chrom": 4,
     "tss": 110610001, "tes": 110633149,
     "win": (109610001, 111633149)},
    {"symbol": "KCNN3",  "axis": "positive_control",  "group": "pos_ctrl",
     "gene_id": "ENSG00000143603", "chrom": 1,
     "tss": 154837802, "tes": 154961963,
     "win": (153837802, 155961963)},
]


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "axes.titlesize":    9,
        "axes.labelsize":    8,
        "xtick.labelsize":   6.5,
        "ytick.labelsize":   6.5,
        "legend.fontsize":   6.5,
        "axes.linewidth":    0.6,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def load_eqtl_for_anchor(a):
    win_lo, win_hi = a["win"]
    fn = CACHE / f"eqtl_{a['symbol']}_chr{a['chrom']}_{win_lo}_{win_hi}.tsv.gz"
    if not fn.exists():
        print(f"WARN: missing {fn}")
        return pd.DataFrame()
    df = pd.read_csv(fn, sep="\t", header=None,
                       names=EQTL_NAMES, low_memory=False)
    # filter to anchor gene
    df = df[df["gene_id"].str.startswith(a["gene_id"]) |
              df["molecular_trait_id"].str.startswith(a["gene_id"])]
    df = df.assign(position=df["position"].astype(int),
                     pvalue=df["pvalue"].astype(float))
    df = df[df["pvalue"] > 0]
    return df


def load_gwas_for_anchor(a):
    win_lo, win_hi = a["win"]
    fn = CACHE / f"gwas_{a['symbol']}_chr{a['chrom']}_{win_lo}_{win_hi}.tsv.gz"
    if not fn.exists():
        print(f"WARN: missing {fn}")
        return pd.DataFrame()
    df = pd.read_csv(fn, sep="\t", header=None,
                       names=GWAS_NAMES, low_memory=False)
    df = df.assign(base_pair_location=df["base_pair_location"].astype(int),
                     p_value=df["p_value"].astype(float))
    df = df[df["p_value"] > 0]
    return df


def main() -> None:
    set_publication_style()

    # PP values for annotation
    pp_anchors = pd.read_csv(ROOT / "results" / "figure6_mr_anchors.tsv",
                                  sep="\t")
    pp_pos     = pd.read_csv(ROOT / "results" / "figure6_mr_positive_controls.tsv",
                                  sep="\t")
    pp_pos = pp_pos[pp_pos["tissue"] == "atrial_appendage"]
    pp = pd.concat([pp_anchors[["anchor", "coloc_PP3", "coloc_PP4"]],
                       pp_pos[["anchor", "coloc_PP3", "coloc_PP4"]]],
                      ignore_index=True).set_index("anchor")

    # ---- figure layout : 6 anchor blocks (3 rows × 2 cols), each 2 sub-tracks ----
    n_anchors = len(ANCHORS)
    rows, cols = 3, 2
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(200 * MM, 220 * MM))
    outer = fig.add_gridspec(
        nrows=rows, ncols=cols, hspace=0.55, wspace=0.30,
        left=0.07, right=0.97, top=0.93, bottom=0.05)

    fig.suptitle(
        "Fig S  LocusZoom-style cis-eQTL × AF GWAS mini-Manhattan  "
        "(6 anchors, ±1 Mb cis-window)",
        fontsize=8.8, fontweight="bold", y=0.97, x=0.02, ha="left")

    fig.text(0.02, 0.965,
              "Read each block as: do the eQTL and GWAS −log10(p) peaks "
              "land at the same genomic position?   "
              "If yes → shared causal variant (high PP4).   "
              "If no  → different causal variants (high PP3, low PP4) — "
              "the dominant pattern observed for all 6 anchors.",
              fontsize=6.0, color=GREY_DARK, ha="left", va="top",
              style="italic")

    for k, a in enumerate(ANCHORS):
        rr, cc = divmod(k, cols)
        sub = outer[rr, cc].subgridspec(
            nrows=3, ncols=1, hspace=0.10,
            height_ratios=[3.0, 3.0, 0.55])
        ax_e = fig.add_subplot(sub[0])
        ax_g = fig.add_subplot(sub[1], sharex=ax_e)
        ax_b = fig.add_subplot(sub[2], sharex=ax_e)

        e = load_eqtl_for_anchor(a)
        g = load_gwas_for_anchor(a)

        if e.empty or g.empty:
            for ax in (ax_e, ax_g, ax_b):
                ax.text(0.5, 0.5, f"missing data: {a['symbol']}",
                          transform=ax.transAxes, ha="center", va="center")
                ax.axis("off")
            continue

        # nlog10
        e["nlp"] = -np.log10(e["pvalue"])
        g["nlp"] = -np.log10(g["p_value"])

        # x in Mb
        e["x_mb"] = e["position"] / 1e6
        g["x_mb"] = g["base_pair_location"] / 1e6

        # subsample for plot speed (huge n_GWAS)
        rng = np.random.default_rng(42 + k)
        if len(g) > 6000:
            g_plot = g.sample(6000, random_state=42 + k)
        else:
            g_plot = g
        if len(e) > 6000:
            e_plot = e.sample(6000, random_state=42 + k)
        else:
            e_plot = e

        # eQTL track
        ax_e.scatter(e_plot["x_mb"], e_plot["nlp"],
                       s=2.0, color=EQTL_COLOR, alpha=0.45,
                       edgecolor="none", zorder=2)
        # lead eQTL
        lead_e = e.loc[e["pvalue"].idxmin()]
        ax_e.scatter([lead_e["position"] / 1e6], [-np.log10(lead_e["pvalue"])],
                       marker="D", s=45,
                       facecolor=EQTL_COLOR, edgecolor="black",
                       linewidth=0.6, zorder=4,
                       label=f"lead eQTL  p = {lead_e['pvalue']:.1e}")
        ax_e.axvline(lead_e["position"] / 1e6, color=EQTL_COLOR,
                       lw=0.5, ls=":", alpha=0.6, zorder=1)

        # GWAS track
        ax_g.scatter(g_plot["x_mb"], g_plot["nlp"],
                       s=2.0, color=GWAS_COLOR, alpha=0.45,
                       edgecolor="none", zorder=2)
        lead_g = g.loc[g["p_value"].idxmin()]
        ax_g.scatter([lead_g["base_pair_location"] / 1e6],
                       [-np.log10(lead_g["p_value"])],
                       marker="D", s=45,
                       facecolor=GWAS_COLOR, edgecolor="black",
                       linewidth=0.6, zorder=4,
                       label=f"lead GWAS  p = {lead_g['p_value']:.1e}")
        ax_g.axvline(lead_g["base_pair_location"] / 1e6, color=GWAS_COLOR,
                       lw=0.5, ls=":", alpha=0.6, zorder=1)

        # gene body band
        for ax in (ax_e, ax_g, ax_b):
            ax.axvspan(a["tss"] / 1e6, a["tes"] / 1e6,
                         color=HL_GREEN, alpha=0.20, zorder=0)

        # gene-body track (bottom)
        ax_b.add_patch(mpatches.Rectangle(
            (a["tss"] / 1e6, 0.4), (a["tes"] - a["tss"]) / 1e6, 0.4,
            facecolor=HL_GREEN, edgecolor="none", alpha=0.7,
            zorder=3))
        ax_b.text((a["tss"] + a["tes"]) / 2 / 1e6, 0.6, a["symbol"],
                    ha="center", va="center",
                    fontsize=7.0, fontweight="bold", color="white",
                    zorder=4)
        ax_b.set_ylim(0, 1); ax_b.set_yticks([])
        ax_b.set_xlabel(f"chromosome {a['chrom']}  position (Mb)",
                          fontsize=7)
        ax_b.spines["left"].set_visible(False)

        # title with PP3 / PP4
        if a["symbol"] in pp.index:
            pp3 = float(pp.loc[a["symbol"], "coloc_PP3"])
            pp4 = float(pp.loc[a["symbol"], "coloc_PP4"])
            pp_text = f"PP3 = {pp3:.2f}   PP4 = {pp4:.2f}"
        else:
            pp_text = ""
        axis_col = (RA_COLOR if a["axis"] == "matrix"
                       else (LA_COLOR if a["axis"] == "metabolic" else HL_GREEN))
        title_lbl = ("▶ " if a["group"] == "anchor" else "★ ") + a["symbol"]
        ax_e.set_title(
            f"{title_lbl}    [{a['axis'].replace('_', ' ')}]    "
            f"·    {pp_text}",
            loc="left", fontsize=7.8, fontweight="bold",
            color=axis_col, pad=4)

        ax_e.set_ylabel("−log10 p\neQTL", fontsize=6.5)
        ax_g.set_ylabel("−log10 p\nGWAS", fontsize=6.5)

        # genomic significance reference
        ax_e.axhline(-np.log10(5e-8), color="grey", lw=0.4, ls="--",
                       alpha=0.5)
        ax_g.axhline(-np.log10(5e-8), color="grey", lw=0.4, ls="--",
                       alpha=0.5)

        # tighten the x-range to the actual cache window
        ax_e.set_xlim(a["win"][0] / 1e6 - 0.02, a["win"][1] / 1e6 + 0.02)

        # legends
        ax_e.legend(loc="upper right", frameon=False, fontsize=5.8)
        ax_g.legend(loc="upper right", frameon=False, fontsize=5.8)

        # hide x labels on top tracks
        plt.setp(ax_e.get_xticklabels(), visible=False)
        plt.setp(ax_g.get_xticklabels(), visible=False)

    fig.text(0.02, 0.020,
              "Diamonds: lead variant in each track.   "
              "Green band: anchor gene body  (TSS - TES).   "
              "Dashed grey horizontal: −log10(5e-8) (GWAS-grade).   "
              "▶ = framework anchor;  ★ = positive control.   "
              "Note that for several anchors the lead eQTL and lead GWAS "
              "diamonds fall at clearly distinct genomic positions, "
              "consistent with the high-PP3 / low-PP4 pattern.",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
