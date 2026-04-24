#!/usr/bin/env python3
"""
Figure 6E  Classical single-gene mRNA cis-eQTL × AF GWAS regime check
            with tiered, eQTL-appropriate thresholds + LocusZoom visual.

Three panels:
    a  Per-anchor tiered evidence table (n = 4 framework anchors + 2
        positive controls).  Now includes PP3, PP4, and the physical
        distance (kb) between lead eQTL and lead GWAS variant — the
        single quantity that converts the high-PP3 / low-PP4 pattern
        into a physical, interpretable claim.

    b  PP4 × −log10(eQTL p) tiered scatter (legacy panel).  Background
        shading marks the 4 PP4 zones (FAIL / EXPLORATORY / STANDARD /
        CONFIRMATORY) and 2 eQTL p reference lines.

    c  PITX2 mini-LocusZoom (canonical 4q25 AF locus).  eQTL track
        (atrial appendage GTEx) + GWAS track (AF GWAS) over the ±1 Mb
        cis-window.  Lead variants of the two tracks fall ~ 580 kb
        apart — the visual proof that the high-PP3 / low-PP4 verdict
        for PITX2 is not a power problem but a physical one.

Threshold philosophy (mirrors GTEx / eQTLGen practice):
    eQTL evidence
        STRINGENT  : raw p < 5e-8
        STANDARD   : raw p < 1e-5  (≈ Bonferroni for ~200 effective
                                     tests in cis-window)
        EXPLORATORY: BH-FDR < 0.05 across the 6 anchors
    Coloc evidence
        CONFIRMATORY: PP4 ≥ 0.90
        STANDARD    : PP4 ≥ 0.70
        EXPLORATORY : PP4 ≥ 0.50

Inputs:
    results/figure6_mr_anchors.tsv             (atrial appendage)
    results/figure6_mr_positive_controls.tsv   (PITX2, KCNN3 sensitivity)
    Project_P6/cache/{eqtl,gwas}_<anchor>_<window>.tsv.gz
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "Project_P6" / "cache"

ANCHORS_TSV = ROOT / "results" / "figure6_mr_anchors.tsv"
POS_TSV     = ROOT / "results" / "figure6_mr_positive_controls.tsv"
FIG_OUT     = ROOT / "figures" / "fig6e_classical_regime_check.pdf"

GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"
LA_COLOR   = "#C97A74"
RA_COLOR   = "#7D99B8"
HL_GREEN   = "#3C8C57"
COL_PASS   = "#1B7837"
COL_PARTIAL = "#D4A35E"
COL_FAIL   = "#B2182B"
EQTL_COLOR = "#7E5BBF"
GWAS_COLOR = "#C25450"

EQTL_STRINGENT  = 5e-8
EQTL_STANDARD   = 1e-5
PP4_CONFIRMATORY = 0.90
PP4_STANDARD     = 0.70
PP4_EXPLORATORY  = 0.50

EQTL_NAMES = ["molecular_trait_id", "chromosome", "position", "ref", "alt",
               "variant", "ma_samples", "maf", "pvalue", "beta", "se",
               "type", "ac", "an", "r2", "molecular_trait_object_id",
               "gene_id", "median_tpm", "rsid", "nlog10p"]
GWAS_NAMES = ["chromosome", "base_pair_location", "effect_allele",
               "other_allele", "beta", "standard_error",
               "effect_allele_frequency", "p_value", "rsid", "rs_id",
               "n", "hm_coordinate_conversion", "hm_code", "variant_id"]

ANCHOR_META = {
    "PRELP":  {"gene_id": "ENSG00000188783", "chrom": 1,
                 "tss": 203475828, "tes": 203491352,
                 "win": (202475828, 204491352)},
    "OMD":    {"gene_id": "ENSG00000127083", "chrom": 9,
                 "tss": 92414245, "tes": 92424461,
                 "win": (91414245, 93424461)},
    "IDH2":   {"gene_id": "ENSG00000182054", "chrom": 15,
                 "tss": 90083045, "tes": 90102504,
                 "win": (89083045, 91102504)},
    "NDUFA8": {"gene_id": "ENSG00000119421", "chrom": 9,
                 "tss": 122144058, "tes": 122159819,
                 "win": (121144058, 123159819)},
    "PITX2":  {"gene_id": "ENSG00000164093", "chrom": 4,
                 "tss": 110610001, "tes": 110633149,
                 "win": (109610001, 111633149)},
    "KCNN3":  {"gene_id": "ENSG00000143603", "chrom": 1,
                 "tss": 154837802, "tes": 154961963,
                 "win": (153837802, 155961963)},
}


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
        "legend.fontsize":   6.5,
        "axes.linewidth":    0.6,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def gate_status(value: float, gates: list, ascending: bool):
    if pd.isna(value):
        return -1
    tier = 0
    for i, g in enumerate(gates, start=1):
        if (ascending and value < g) or (not ascending and value >= g):
            tier = i
        else:
            break
    return tier


def load_eqtl(symbol: str):
    a = ANCHOR_META[symbol]
    win = a["win"]
    fn = CACHE / f"eqtl_{symbol}_chr{a['chrom']}_{win[0]}_{win[1]}.tsv.gz"
    if not fn.exists():
        return pd.DataFrame()
    df = pd.read_csv(fn, sep="\t", header=None,
                       names=EQTL_NAMES, low_memory=False)
    df = df[df["gene_id"].astype(str).str.startswith(a["gene_id"])
              | df["molecular_trait_id"].astype(str).str.startswith(a["gene_id"])]
    df = df.assign(position=df["position"].astype(int),
                     pvalue=df["pvalue"].astype(float))
    return df[df["pvalue"] > 0]


def load_gwas(symbol: str):
    a = ANCHOR_META[symbol]
    win = a["win"]
    fn = CACHE / f"gwas_{symbol}_chr{a['chrom']}_{win[0]}_{win[1]}.tsv.gz"
    if not fn.exists():
        return pd.DataFrame()
    df = pd.read_csv(fn, sep="\t", header=None,
                       names=GWAS_NAMES, low_memory=False)
    df = df.assign(base_pair_location=df["base_pair_location"].astype(int),
                     p_value=df["p_value"].astype(float))
    return df[df["p_value"] > 0]


def lead_offset_kb(symbol: str) -> float:
    e = load_eqtl(symbol); g = load_gwas(symbol)
    if e.empty or g.empty:
        return float("nan")
    e_lead = int(e.loc[e["pvalue"].idxmin(), "position"])
    g_lead = int(g.loc[g["p_value"].idxmin(), "base_pair_location"])
    return abs(e_lead - g_lead) / 1e3


def main() -> None:
    set_publication_style()

    anchors = pd.read_csv(ANCHORS_TSV, sep="\t").assign(group="anchor")
    pos     = pd.read_csv(POS_TSV, sep="\t")
    pos     = pos[pos["tissue"] == "atrial_appendage"].assign(group="pos_ctrl")
    df = pd.concat([
        anchors[["anchor", "axis", "n_eqtl", "lead_eqtl_p",
                  "coloc_PP3", "coloc_PP4", "verdict", "group"]],
        pos[["anchor", "axis", "n_eqtl", "lead_eqtl_p",
              "coloc_PP3", "coloc_PP4", "verdict", "group"]],
    ], ignore_index=True)
    df["fdr_bh"] = false_discovery_control(
        df["lead_eqtl_p"].to_numpy(), method="bh")
    df["offset_kb"] = df["anchor"].apply(lead_offset_kb)

    # row order
    axis_order = {"matrix": 0, "metabolic": 1, "positive_control": 2}
    df["__o"] = df["axis"].map(axis_order)
    df = df.sort_values(["__o", "anchor"]).reset_index(drop=True)

    eqtl_tiers_p = [1.0, EQTL_STANDARD, EQTL_STRINGENT]
    pp4_tiers   = [PP4_EXPLORATORY, PP4_STANDARD, PP4_CONFIRMATORY]

    rows = []
    for _, r in df.iterrows():
        t_p   = gate_status(r["lead_eqtl_p"], eqtl_tiers_p, ascending=True)
        t_fdr = 1 if r["fdr_bh"] < 0.05 else 0
        eqtl_tier  = max(t_p, t_fdr)
        c_tier = gate_status(r["coloc_PP4"], pp4_tiers, ascending=False)
        rows.append({"anchor": r["anchor"], "axis": r["axis"],
                       "group": r["group"],
                       "lead_p": r["lead_eqtl_p"], "fdr": r["fdr_bh"],
                       "PP3": r["coloc_PP3"], "PP4": r["coloc_PP4"],
                       "offset_kb": r["offset_kb"],
                       "eqtl_tier": eqtl_tier,
                       "coloc_tier": c_tier})
    tier = pd.DataFrame(rows)

    print("=== anchor tier table (with PP3, offset) ===")
    print(tier.to_string(index=False))

    # ---- figure layout: top (a + b), bottom (c PITX2 + d KCNN3) ----
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(230 * MM, 175 * MM))
    outer = fig.add_gridspec(
        nrows=2, ncols=1,
        height_ratios=[1.05, 1.0],
        hspace=0.34,
        left=0.05, right=0.97, top=0.91, bottom=0.06)
    top = outer[0].subgridspec(nrows=1, ncols=2, wspace=0.18,
                                  width_ratios=[1.55, 1.0])
    # bottom : 2-col grid; each col itself a 3-row LocusZoom block
    bot_outer = outer[1].subgridspec(nrows=1, ncols=2, wspace=0.20)
    bot_pitx2 = bot_outer[0].subgridspec(
        nrows=3, ncols=1, hspace=0.10,
        height_ratios=[3.5, 3.5, 0.6])
    bot_kcnn3 = bot_outer[1].subgridspec(
        nrows=3, ncols=1, hspace=0.10,
        height_ratios=[3.5, 3.5, 0.6])

    fig.suptitle(
        "Fig 6E  Classical single-gene mRNA cis-eQTL × AF GWAS coloc "
        "regime check  (GTEx atrial appendage, n = 372)",
        fontsize=8.7, fontweight="bold", y=0.965, x=0.02, ha="left")

    # ===================================================================
    # PANEL a: anchor table
    # ===================================================================
    ax = fig.add_subplot(top[0])
    ax.axis("off")
    ax.set_title("a    Per-anchor tiered evidence",
                   loc="left", fontsize=8.4, fontweight="bold", pad=4)

    col_xs = {
        "name":   0.025,
        "axis":   0.135,
        "fdr":    0.300,
        "pp3":    0.380,
        "pp4":    0.450,
        "off":    0.530,
        "etier":  0.640,
        "ctier":  0.815,
    }
    headers = [
        ("name",  "Anchor"),
        ("axis",  "Axis"),
        ("fdr",   "BH-FDR"),
        ("pp3",   "PP3"),
        ("pp4",   "PP4"),
        ("off",   "lead-offset (kb)"),
        ("etier", "eQTL tier"),
        ("ctier", "coloc tier"),
    ]
    header_y = 0.94
    for k, label in headers:
        ax.text(col_xs[k], header_y, label,
                  ha="left", va="center", fontsize=7.0, fontweight="bold",
                  color=GREY_DARK, transform=ax.transAxes)
    ax.plot([0.005, 0.985], [header_y - 0.022, header_y - 0.022],
              color=GREY_MID, lw=0.6, transform=ax.transAxes)

    EQTL_PILL_LABELS = {
        0: ("FAIL",         COL_FAIL),
        1: ("EXPLORATORY",  COL_PARTIAL),
        2: ("STANDARD",     COL_PASS),
        3: ("STRINGENT",    "#0F4D2A"),
    }
    COLOC_PILL_LABELS = {
        -1: ("n/a",          GREY_MID),
        0:  ("FAIL",         COL_FAIL),
        1:  ("EXPLORATORY",  COL_PARTIAL),
        2:  ("STANDARD",     COL_PASS),
        3:  ("CONFIRMATORY", "#0F4D2A"),
    }

    row_h = 0.085
    for i, r in tier.iterrows():
        y = header_y - 0.04 - (i + 0.5) * row_h
        if r["axis"] == "matrix":
            axis_col = RA_COLOR; axis_lbl = "matrix (RA)"
        elif r["axis"] == "metabolic":
            axis_col = LA_COLOR; axis_lbl = "metabolic (LA)"
        else:
            axis_col = HL_GREEN; axis_lbl = "positive ctrl"

        # icon
        ax.text(0.005, y,
                  "\u25b6" if r["group"] == "anchor" else "\u2605",
                  ha="center", va="center",
                  fontsize=9 if r["group"] == "anchor" else 11,
                  color=axis_col, transform=ax.transAxes)
        ax.text(col_xs["name"], y, r["anchor"],
                  ha="left", va="center", fontsize=7.6, fontweight="bold",
                  color=axis_col, transform=ax.transAxes)
        ax.text(col_xs["axis"], y, axis_lbl,
                  ha="left", va="center", fontsize=6.4,
                  color=axis_col, transform=ax.transAxes)
        ax.text(col_xs["fdr"], y, f"{r['fdr']:.2e}",
                  ha="left", va="center", fontsize=6.6,
                  color=GREY_DARK, transform=ax.transAxes)

        # PP3 (high)  - bold red if > 0.30
        pp3_color = COL_FAIL if r["PP3"] > 0.30 else GREY_DARK
        ax.text(col_xs["pp3"], y, f"{r['PP3']:.2f}",
                  ha="left", va="center", fontsize=6.6,
                  color=pp3_color, fontweight="bold",
                  transform=ax.transAxes)
        # PP4 (low)
        ax.text(col_xs["pp4"], y, f"{r['PP4']:.2f}",
                  ha="left", va="center", fontsize=6.6,
                  color=GREY_DARK,
                  transform=ax.transAxes)
        # offset (kb) - bold dark if > 200 kb (clear separation)
        off = r["offset_kb"]
        off_color = "black" if off > 200 else GREY_DARK
        off_str = f"{off:,.0f}" if pd.notna(off) else "n/a"
        ax.text(col_xs["off"], y, off_str,
                  ha="left", va="center", fontsize=6.8,
                  color=off_color, fontweight="bold",
                  transform=ax.transAxes)

        # eQTL pill
        e_lbl, e_col = EQTL_PILL_LABELS[r["eqtl_tier"]]
        pill = mpatches.FancyBboxPatch(
            (col_xs["etier"] - 0.005, y - 0.026), 0.155, 0.052,
            boxstyle="round,pad=0.003,rounding_size=0.010",
            facecolor=e_col, edgecolor=e_col,
            linewidth=0.6, transform=ax.transAxes, zorder=3)
        ax.add_patch(pill)
        ax.text(col_xs["etier"] + 0.072, y, e_lbl,
                  ha="center", va="center", fontsize=6.3, fontweight="bold",
                  color="white", transform=ax.transAxes, zorder=4)

        # coloc pill
        c_lbl, c_col = COLOC_PILL_LABELS[r["coloc_tier"]]
        pill = mpatches.FancyBboxPatch(
            (col_xs["ctier"] - 0.005, y - 0.026), 0.165, 0.052,
            boxstyle="round,pad=0.003,rounding_size=0.010",
            facecolor=c_col, edgecolor=c_col,
            linewidth=0.6, transform=ax.transAxes, zorder=3)
        ax.add_patch(pill)
        ax.text(col_xs["ctier"] + 0.077, y, c_lbl,
                  ha="center", va="center", fontsize=6.3, fontweight="bold",
                  color="white", transform=ax.transAxes, zorder=4)

        ax.plot([0.005, 0.985], [y - row_h / 2, y - row_h / 2],
                  color="#EEEEEE", lw=0.3, transform=ax.transAxes)

    # legend
    ax.text(0.005, 0.075,
              "Tier definitions:",
              fontsize=6.5, fontweight="bold", color=GREY_DARK,
              transform=ax.transAxes)
    ax.text(0.005, 0.046,
              f"eQTL  ::  EXPLORATORY (BH-FDR < 0.05)  ·  "
              f"STANDARD (p < {EQTL_STANDARD:.0e})  ·  "
              f"STRINGENT (p < {EQTL_STRINGENT:.0e})",
              fontsize=5.7, color=GREY_DARK,
              transform=ax.transAxes)
    ax.text(0.005, 0.018,
              f"coloc ::  EXPLORATORY (PP4 ≥ {PP4_EXPLORATORY:.1f})  ·  "
              f"STANDARD (PP4 ≥ {PP4_STANDARD:.1f})  ·  "
              f"CONFIRMATORY (PP4 ≥ {PP4_CONFIRMATORY:.1f})",
              fontsize=5.7, color=GREY_DARK,
              transform=ax.transAxes)
    ax.text(0.62, 0.046,
              "PP3 = different causal variants per coloc.abf;  "
              "PP3 > 0.30 marked bold red.",
              fontsize=5.7, color=GREY_DARK,
              transform=ax.transAxes)
    ax.text(0.62, 0.018,
              "lead-offset (kb) = physical distance between the two "
              "tracks' lead variants;  > 200 kb marked bold black.",
              fontsize=5.7, color=GREY_DARK,
              transform=ax.transAxes)

    # ===================================================================
    # PANEL b: PP4 × −log10 p scatter
    # ===================================================================
    axb = fig.add_subplot(top[1])
    axb.set_title("b    Distance to PP4 / eQTL gates",
                    loc="left", fontsize=8.4, fontweight="bold", pad=4)

    NLP_MAX = max(8, -np.log10(tier["lead_p"]).max() + 0.5)
    axb.axvspan(0,                PP4_EXPLORATORY,    alpha=0.06,
                  color=COL_FAIL, zorder=0)
    axb.axvspan(PP4_EXPLORATORY,  PP4_STANDARD,      alpha=0.07,
                  color=COL_PARTIAL, zorder=0)
    axb.axvspan(PP4_STANDARD,     PP4_CONFIRMATORY,  alpha=0.10,
                  color=COL_PASS, zorder=0)
    axb.axvspan(PP4_CONFIRMATORY, 1.0,                alpha=0.18,
                  color=COL_PASS, zorder=0)

    axb.axvline(PP4_EXPLORATORY,  color=COL_PARTIAL, lw=0.6, ls="--", zorder=1)
    axb.axvline(PP4_STANDARD,     color=COL_PASS,    lw=0.7, ls="--", zorder=1)
    axb.axvline(PP4_CONFIRMATORY, color="#0F4D2A",   lw=0.7, ls="--", zorder=1)
    axb.axhline(-np.log10(EQTL_STANDARD),  color=COL_PASS,
                  lw=0.6, ls=":", zorder=1)
    axb.axhline(-np.log10(EQTL_STRINGENT), color="#0F4D2A",
                  lw=0.6, ls=":", zorder=1)

    for _, r in tier.iterrows():
        col = (RA_COLOR if r["axis"] == "matrix"
                else (LA_COLOR if r["axis"] == "metabolic" else HL_GREEN))
        marker = "o" if r["group"] == "anchor" else "*"
        size = 110 if marker == "*" else 70
        axb.scatter(r["PP4"], -np.log10(r["lead_p"]),
                      marker=marker, s=size,
                      facecolor=col, edgecolor="black", linewidth=0.6,
                      zorder=4)
        axb.annotate(r["anchor"],
                       (r["PP4"], -np.log10(r["lead_p"])),
                       xytext=(7, 4), textcoords="offset points",
                       fontsize=7.0, fontweight="bold",
                       color=col, zorder=5)

    axb.set_xlim(-0.02, 1.0)
    axb.set_ylim(0, NLP_MAX)
    axb.set_xlabel("coloc PP4")
    axb.set_ylabel("−log10 (lead cis-eQTL p)")

    axb.text(PP4_EXPLORATORY / 2,           NLP_MAX * 0.97,
                "FAIL", ha="center", va="top",
                fontsize=6.0, color=COL_FAIL, alpha=0.7, fontweight="bold")
    axb.text((PP4_EXPLORATORY + PP4_STANDARD) / 2, NLP_MAX * 0.97,
                "EXPL", ha="center", va="top",
                fontsize=6.0, color=COL_PARTIAL, alpha=0.85, fontweight="bold")
    axb.text((PP4_STANDARD + PP4_CONFIRMATORY) / 2, NLP_MAX * 0.97,
                "STD", ha="center", va="top",
                fontsize=6.0, color=COL_PASS, alpha=0.85, fontweight="bold")
    axb.text((PP4_CONFIRMATORY + 1.0) / 2, NLP_MAX * 0.97,
                "CONF", ha="center", va="top",
                fontsize=6.0, color="#0F4D2A", alpha=0.85, fontweight="bold")

    handles = [
        plt.Line2D([], [], marker="o", linestyle="",
                       markerfacecolor="#888", markeredgecolor="black",
                       markersize=6, markeredgewidth=0.4,
                       label="framework anchor (\u25b6)"),
        plt.Line2D([], [], marker="*", linestyle="",
                       markerfacecolor="#888", markeredgecolor="black",
                       markersize=10, markeredgewidth=0.4,
                       label="positive control (\u2605)"),
    ]
    axb.legend(handles=handles, loc="lower right",
                  frameon=False, fontsize=6.2)

    # ===================================================================
    # PANEL c / d: PITX2 + KCNN3 mini-LocusZooms (positive controls)
    # ===================================================================
    def render_locuszoom(symbol: str, gridspec, panel_letter: str,
                            cytoband: str, pp3: float, pp4: float):
        ax_e = fig.add_subplot(gridspec[0])
        ax_g = fig.add_subplot(gridspec[1], sharex=ax_e)
        ax_b = fig.add_subplot(gridspec[2], sharex=ax_e)

        a = ANCHOR_META[symbol]
        e = load_eqtl(symbol); g = load_gwas(symbol)
        e["nlp"] = -np.log10(e["pvalue"])
        g["nlp"] = -np.log10(g["p_value"])
        e["x_mb"] = e["position"] / 1e6
        g["x_mb"] = g["base_pair_location"] / 1e6

        e_plot = e.sample(min(len(e), 6000), random_state=7)
        g_plot = g.sample(min(len(g), 6000), random_state=7)

        ax_e.scatter(e_plot["x_mb"], e_plot["nlp"],
                       s=2.5, color=EQTL_COLOR, alpha=0.45,
                       edgecolor="none", zorder=2)
        lead_e = e.loc[e["pvalue"].idxmin()]
        ax_e.scatter([lead_e["position"] / 1e6],
                       [-np.log10(lead_e["pvalue"])],
                       marker="D", s=55,
                       facecolor=EQTL_COLOR, edgecolor="black",
                       linewidth=0.7, zorder=4,
                       label=(f"lead eQTL  p = {lead_e['pvalue']:.1e}\n"
                                 f"@ {lead_e['position'] / 1e6:.3f} Mb"))
        ax_e.axvline(lead_e["position"] / 1e6, color=EQTL_COLOR,
                       lw=0.6, ls=":", alpha=0.7, zorder=1)
        ax_e.axhline(-np.log10(5e-8), color="grey", lw=0.4, ls="--",
                       alpha=0.5)

        ax_g.scatter(g_plot["x_mb"], g_plot["nlp"],
                       s=2.5, color=GWAS_COLOR, alpha=0.45,
                       edgecolor="none", zorder=2)
        lead_g = g.loc[g["p_value"].idxmin()]
        ax_g.scatter([lead_g["base_pair_location"] / 1e6],
                       [-np.log10(lead_g["p_value"])],
                       marker="D", s=55,
                       facecolor=GWAS_COLOR, edgecolor="black",
                       linewidth=0.7, zorder=4,
                       label=(f"lead GWAS  p = {lead_g['p_value']:.1e}\n"
                                 f"@ {lead_g['base_pair_location'] / 1e6:.3f} Mb"))
        ax_g.axvline(lead_g["base_pair_location"] / 1e6, color=GWAS_COLOR,
                       lw=0.6, ls=":", alpha=0.7, zorder=1)
        ax_g.axhline(-np.log10(5e-8), color="grey", lw=0.4, ls="--",
                       alpha=0.5)

        # gene body band on all 3 axes
        for ax in (ax_e, ax_g, ax_b):
            ax.axvspan(a["tss"] / 1e6, a["tes"] / 1e6,
                         color=HL_GREEN, alpha=0.20, zorder=0)

        # gene body track
        ax_b.add_patch(mpatches.Rectangle(
            (a["tss"] / 1e6, 0.4), (a["tes"] - a["tss"]) / 1e6, 0.4,
            facecolor=HL_GREEN, edgecolor="none", alpha=0.8, zorder=3))
        ax_b.text((a["tss"] + a["tes"]) / 2 / 1e6, 0.6, symbol,
                    ha="center", va="center", fontsize=7.0,
                    fontweight="bold", color="white", zorder=4)
        ax_b.set_ylim(0, 1); ax_b.set_yticks([])
        ax_b.set_xlabel(f"chromosome {a['chrom']}  position (Mb)",
                          fontsize=6.8)
        ax_b.spines["left"].set_visible(False)

        ax_e.set_ylabel("−log10 p\neQTL", fontsize=6.5)
        ax_g.set_ylabel("−log10 p\nGWAS", fontsize=6.5)

        ax_e.legend(loc="upper right", frameon=False, fontsize=5.8)
        ax_g.legend(loc="upper right", frameon=False, fontsize=5.8)

        plt.setp(ax_e.get_xticklabels(), visible=False)
        plt.setp(ax_g.get_xticklabels(), visible=False)

        # distance between leads (kb)
        e_x = lead_e["position"] / 1e6
        g_x = lead_g["base_pair_location"] / 1e6
        offset_kb = abs(e_x - g_x) * 1000

        ax_e.set_title(
            f"{panel_letter}    {symbol} ({cytoband})  cis-eQTL × AF GWAS    "
            f"·    {offset_kb:.0f} kb apart    ·    "
            f"PP3 = {pp3:.2f}   PP4 = {pp4:.2f}",
            loc="left", fontsize=7.8, fontweight="bold",
            color=HL_GREEN, pad=4)

        # double-headed arrow on gene-body track
        yarr = 0.05
        ax_b.annotate(
            "", xy=(g_x, yarr), xytext=(e_x, yarr),
            arrowprops=dict(arrowstyle="<->", color=GREY_DARK, lw=1.0,
                              shrinkA=0, shrinkB=0),
            zorder=5)
        ax_b.text((e_x + g_x) / 2, yarr - 0.15,
                     f"{offset_kb:.0f} kb",
                     ha="center", va="top", fontsize=6.3,
                     fontweight="bold", color=GREY_DARK)

    pitx2_pp3 = float(tier.loc[tier["anchor"] == "PITX2", "PP3"].iloc[0])
    pitx2_pp4 = float(tier.loc[tier["anchor"] == "PITX2", "PP4"].iloc[0])
    kcnn3_pp3 = float(tier.loc[tier["anchor"] == "KCNN3", "PP3"].iloc[0])
    kcnn3_pp4 = float(tier.loc[tier["anchor"] == "KCNN3", "PP4"].iloc[0])

    render_locuszoom("PITX2", bot_pitx2, "c", "chr 4q25",
                        pitx2_pp3, pitx2_pp4)
    render_locuszoom("KCNN3", bot_kcnn3, "d", "chr 1q21.3",
                        kcnn3_pp3, kcnn3_pp4)

    fig.text(0.02, 0.020,
              "Narrative:  6/6 anchors PASS the eQTL EXPLORATORY tier "
              "(BH-FDR < 0.05);  IDH2 also passes STANDARD.   "
              "0/6 anchors reach the coloc EXPLORATORY tier (PP4 ≥ 0.50);  "
              "all 6 carry high PP3 (different causal variants) and a "
              "large physical lead-variant offset (459-1,214 kb).   "
              "Both canonical AF positive controls — PITX2 (panel c) and "
              "KCNN3 (panel d, smallest offset 459 kb) — fail equally.",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
