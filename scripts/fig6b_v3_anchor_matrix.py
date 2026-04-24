#!/usr/bin/env python3
"""
Figure 6B (v3) — 10-anchor + PITX2-calibrator evidence matrix.

Reads results/figure6_mr_anchors_v3_consolidated.tsv (the per-anchor table
built by build_v3_consolidated.py from the v2 + v3 + PITX2 raw outputs;
verdict / shape_tag / family_verdict already applied per spec §10.5–10.9).

Renders one row per anchor (10 anchors, 4 v2 + 6 v3) plus a divider and
one row for the PITX2 calibrator at the bottom.  Columns:

    Family · Anchor · Source · cis-eQTL primary · cis-eQTL blood (sens) ·
    AF GWAS min p · PP3 · PP4 · Verdict · Shape · Family verdict

Verdict colour (locked, spec §10.10):
    supported          #1B7837 (green)
    weak_hybrid        #E08214 (deep orange)
    not_supported      #B2182B (red)
    no_eligible_eQTL   #9E9E9E (grey)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

from plot_style import (COL_DOUBLE, FG_COLOR, MUTED_COLOR,
                        set_nature_style)

ROOT = Path(__file__).resolve().parent.parent
IN  = ROOT / "results" / "figure6_mr_anchors_v3_consolidated.tsv"
SUM = ROOT / "results" / "figure6_mr_anchors_v3_summary.tsv"
OUT = ROOT / "figures" / "fig6b_v3_anchor_matrix.pdf"

VERDICT_COL = {
    "supported":         "#1B7837",
    "weak_hybrid":       "#E08214",
    "not_supported":     "#B2182B",
    "no_eligible_eQTL":  "#9E9E9E",
}
VERDICT_LBL = {
    "supported":         "supported",
    "weak_hybrid":       "weak hybrid",
    "not_supported":     "not supported",
    "no_eligible_eQTL":  "no eligible eQTL",
}

SOURCE_LBL = {
    "v2_pre_registered":   "v2 pre-reg.",
    "v3_LE_extension":     "v3 LE ext.",
    "v2_positive_control": "v2 calibrator",
}

FAMILY_LBL = {
    "F8":         "F8 RA-driven Metabolic",
    "F10":        "F10 Matrix",
    "F11":        "F11 BMP (RA-biased)",
    "F15":        "F15 NTRK2 (RA-biased)",
    "Metabolic":  "Metabolic",
    "calibrator": "Calibrator (PITX2)",
}

SHAPE_COL = {
    "B":                       "#7B3294",   # PITX2-mode (purple)
    "C":                       "#A6611A",   # eQTL absence (brown)
    "fail-other":              "#BDBDBD",
    "fail-other_B-demoted":    "#BDBDBD",
    "weak":                    "#FDB863",
    "pass":                    "#1B7837",
    "no_data":                 "#FFFFFF",
}
SHAPE_LBL = {
    "B": "B",
    "C": "C",
    "fail-other": "—",
    "fail-other_B-demoted": "—†",
    "weak": "w",
    "pass": "✓",
    "no_data": "n/a",
}


def fmt_p(p) -> str:
    if p is None or pd.isna(p):
        return "n/a"
    p = float(p)
    if p == 0:
        return "<1e-300"
    if p < 1e-3:
        return f"{p:.1e}".replace("e-0", "e-")
    return f"{p:.3f}"


def fmt_pp(pp) -> str:
    if pp is None or pd.isna(pp):
        return "n/a"
    return f"{float(pp):.2f}"


def coloc_color(pp4):
    """For coloring PP4 cell to give a quick visual sense; not driving verdict."""
    if pd.isna(pp4):
        return "#F2F2F2"
    pp4 = float(pp4)
    if pp4 >= 0.7:
        return "#A6DBA0"
    if pp4 >= 0.5:
        return "#F1B377"   # weak hybrid orange
    if pp4 >= 0.3:
        return "#E0E0E0"
    return "#F4ACA0"


def gate_glyph(p_eqtl, threshold=5e-8, marginal=1e-4):
    if pd.isna(p_eqtl):
        return ("—", "#AAAAAA")
    p = float(p_eqtl)
    if p < threshold:
        return ("●", "#1B7837")
    if p < marginal:
        return ("◐", "#888888")
    return ("○", "#B2182B")


def _round_box(ax, x, y, w, h, fc, ec="white", lw=0.4, z=2, rounding=0.012):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    ))


def main() -> None:
    set_nature_style()
    # keep_default_na=False prevents pandas from turning the literal string
    # "NA" (used in family_verdict for the calibrator row) into NaN.
    df = pd.read_csv(IN, sep="\t", keep_default_na=False, na_values=[""])
    # Re-cast numeric columns now that na_values is restricted to "".
    for c in ("lead_eqtl_p", "n_eqtl", "n_gwas", "n_overlap", "gwas_min_p",
              "coloc_PP0", "coloc_PP1", "coloc_PP2", "coloc_PP3", "coloc_PP4",
              "lead_eqtl_p_blood", "coloc_PP4_blood", "coloc_PP3_blood"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    summ = pd.read_csv(SUM, sep="\t").iloc[0]

    # Split anchors vs calibrator; keep canonical order
    anchors = df.query("kind == 'anchor'").copy()
    calib   = df.query("kind == 'calibrator'").copy()
    df = pd.concat([anchors, calib], ignore_index=True)
    n_anchor = len(anchors)
    n_row    = len(df)

    # ---- column geometry (axes 0..1) --------------------------------
    # PP3+PP4 merged into a single wider "coloc" cell with two horizontal
    # bars sharing a 0–1 axis (PP3 in PITX2-mode purple, PP4 in verdict
    # colour) so that B-mode (high PP3, low PP4) vs C-mode (eQTL absence)
    # failure modes are visually separable at a glance.
    # The legacy 5e-8 GWAS column is replaced by an explicit eGene/Bonf
    # PRIMARY/SENSITIVITY chip pair (Plan C, spec §10.3).
    columns = [
        ("family",   "Family",                                       0.140),
        ("anchor",   "Anchor",                                       0.063),
        ("source",   "Source",                                       0.075),
        ("eligibility", "Plan C eligibility\n(eGene  /  Bonf)",      0.115),
        ("eqtl",     "cis-eQTL primary\nlead p (atrial)",            0.105),
        ("eqtl_b",   "cis-eQTL sens.\nlead p (blood)",               0.105),
        ("gwas",     "AF GWAS\nmin p (window)",                      0.080),
        ("coloc",    "Coloc bars: PP3 / PP4ᵉ(eQTL) / PP4ˢ(sQTL)",    0.140),
        ("verdict",  "Verdict",                                      0.090),
        ("shape",    "Shape",                                        0.040),
        ("famv",     "Family\nverdict",                              0.057),
    ]
    LEFT_PAD = 0.0
    xs = np.cumsum([LEFT_PAD] + [c[2] for c in columns])
    centers = (xs[:-1] + xs[1:]) / 2
    col_x  = dict(zip([c[0] for c in columns], centers))
    col_w  = dict(zip([c[0] for c in columns], [c[2] for c in columns]))
    col_lo = dict(zip([c[0] for c in columns], xs[:-1]))
    col_hi = dict(zip([c[0] for c in columns], xs[1:]))

    fig_w = COL_DOUBLE
    row_h = 0.62          # taller rows to fit PP3 + PP4ᵉ + PP4ˢ 3-bar mini chart
    fig_h = row_h * n_row + 3.40   # extra room for 4 legend rows + 3 captions
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.012, right=0.988, top=0.87, bottom=0.25)
    ax.set_axis_off()

    # ---- header -----------------------------------------------------
    HEAD_TOP = n_row + 0.20
    HEAD_BOT = n_row - 0.55
    HEAD_Y   = (HEAD_TOP + HEAD_BOT) / 2
    _round_box(ax, xs[0], HEAD_BOT, xs[-1] - xs[0],
               HEAD_TOP - HEAD_BOT,
               fc="#F0EEEA", ec="none", z=1, rounding=0.014)
    for (key, label, _w), x in zip(columns, centers):
        ax.text(x, HEAD_Y, label, ha="center", va="center",
                fontsize=5.8, fontweight="bold", color=FG_COLOR,
                linespacing=1.2)

    # ---- rows -------------------------------------------------------
    ROW_HALF = 0.43
    for i, r in df.iterrows():
        y = n_row - 1 - i
        # alt zebra
        if i % 2 == 1:
            ax.add_patch(Rectangle((xs[0], y - ROW_HALF),
                                   xs[-1] - xs[0], 2 * ROW_HALF,
                                   facecolor="#F8F8F8", edgecolor="none",
                                   zorder=0))
        # divider before calibrator
        if i == n_anchor:
            ax.plot([xs[0] + 0.005, xs[-1] - 0.005],
                    [y + ROW_HALF + 0.06, y + ROW_HALF + 0.06],
                    color="#666666", lw=1.0, zorder=0)
            ax.text(xs[0] + 0.012, y + ROW_HALF + 0.02,
                    "POSITIVE-CONTROL CALIBRATOR  —  pre-declared expected shape_tag = B (PITX2-mode)",
                    fontsize=4.8, fontweight="bold", color="#444444",
                    style="italic", ha="left", va="bottom", zorder=3)

        # Family
        fam_lbl = FAMILY_LBL.get(r["family"], r["family"])
        ax.text(col_lo["family"] + 0.005, y, fam_lbl,
                ha="left", va="center", fontsize=5.6,
                color="#333333")

        # Anchor (bold)
        ax.text(col_x["anchor"], y, str(r["anchor"]),
                ha="center", va="center",
                fontsize=7.2, fontweight="bold", color=FG_COLOR)

        # Source chip
        src = r.get("source", "")
        src_lbl = SOURCE_LBL.get(src, src)
        ax.text(col_x["source"], y, src_lbl,
                ha="center", va="center", fontsize=5.0,
                color="#444444", style="italic")

        # ---- Plan C eligibility (eGene + Bonf chip pair) ----
        is_eg = bool(r.get("is_eGene_atrial", False))
        is_bf = bool(r.get("bonf_passed", False))
        # Two side-by-side mini chips
        chip_w, chip_h = 0.040, 0.30
        chip_y = y - chip_h / 2
        chip_x_eg = col_x["eligibility"] - 0.025 - chip_w
        chip_x_bf = col_x["eligibility"] + 0.025
        for cx, ok, lbl in [(chip_x_eg, is_eg, "eGene"),
                             (chip_x_bf, is_bf, "Bonf")]:
            cfc = "#1B7837" if ok else "#E0E0E0"
            ctxt = "white" if ok else "#666666"
            _round_box(ax, cx, chip_y, chip_w, chip_h,
                       fc=cfc, ec="white", lw=0.4, z=2, rounding=0.010)
            ax.text(cx + chip_w / 2, y, lbl,
                    ha="center", va="center", fontsize=4.4,
                    color=ctxt, fontweight="bold", zorder=3)
        # Plan-C agreement annotation (✓/✗)
        agree = (is_eg == is_bf)
        ax.text(col_x["eligibility"], y, "=" if agree else "≠",
                ha="center", va="center", fontsize=6.0,
                color=("#1B7837" if agree else "#B2182B"),
                fontweight="bold", zorder=3)

        # ---- cis-eQTL primary lead p ----
        ax.text(col_x["eqtl"], y,
                f"p = {fmt_p(r['lead_eqtl_p'])}",
                ha="center", va="center", fontsize=5.4,
                color="#444444", family="monospace")

        # ---- cis-eQTL blood (sensitivity) ----
        if pd.isna(r.get("lead_eqtl_p_blood")):
            ax.text(col_x["eqtl_b"], y, "n/a (gene absent)",
                    ha="center", va="center", fontsize=5.0,
                    color="#888888", style="italic")
        else:
            ax.text(col_x["eqtl_b"], y,
                    f"p = {fmt_p(r['lead_eqtl_p_blood'])}",
                    ha="center", va="center", fontsize=5.4,
                    color="#444444", family="monospace")

        # ---- GWAS min p ----
        ax.text(col_x["gwas"], y, fmt_p(r["gwas_min_p"]),
                ha="center", va="center", fontsize=5.6,
                color="#222222", family="monospace")

        # ---- Coloc 3-bar mini-chart (PP3, PP4_eqtl, PP4_sqtl) ----
        # Three stacked bars sharing 0..1 axis:
        #   top    : PP3 (different causal variants, purple)
        #   middle : PP4 from cis-eQTL (verdict colour)
        #   bottom : PP4 from cis-sQTL (green family) OR n/a band if no sGene
        # Spec §10.14 RNA-layer extension.
        track_lo = col_lo["coloc"] + 0.005
        track_hi = col_hi["coloc"] - 0.040     # room for right-edge labels
        track_w  = track_hi - track_lo
        bar_h    = 0.13
        gap      = 0.025
        pp3      = float(r["coloc_PP3"]) if not pd.isna(r["coloc_PP3"]) else 0.0
        pp4      = float(r["coloc_PP4"]) if not pd.isna(r["coloc_PP4"]) else 0.0
        sqtl_pp4 = (float(r["sqtl_PP4_max"])
                    if not pd.isna(r.get("sqtl_PP4_max", float("nan")))
                    else float("nan"))
        is_sgene = bool(r.get("is_sGene_atrial", False))
        v        = r["verdict"]
        pp4_col  = VERDICT_COL.get(v, "#BBBBBB")
        SQTL_COL = "#1B7837"   # green: sQTL layer (distinct from purple PP3)

        # Vertical positions (top-down: PP3, PP4_eqtl, PP4_sqtl)
        pp3_y      = y + gap + bar_h
        pp4_y      = y - bar_h / 2
        sqtl_y     = y - gap - 2 * bar_h - bar_h / 2

        # background troughs + 0.5 reference ticks
        for bar_y in (pp3_y, pp4_y, sqtl_y):
            ax.add_patch(Rectangle((track_lo, bar_y), track_w, bar_h,
                                    facecolor="#F0F0F0", edgecolor="#CCCCCC",
                                    linewidth=0.3, zorder=2))
            ax.plot([track_lo + 0.5 * track_w] * 2,
                    [bar_y + 0.015, bar_y + bar_h - 0.015],
                    color="#999999", lw=0.3, zorder=3)

        # PP3 bar (purple)
        ax.add_patch(Rectangle((track_lo, pp3_y), track_w * pp3, bar_h,
                                facecolor="#7B3294", edgecolor="white",
                                linewidth=0.3, zorder=4))
        ax.text(track_hi + 0.003, pp3_y + bar_h / 2,
                f"{pp3:.2f}", ha="left", va="center", fontsize=4.4,
                color="#7B3294", fontweight="bold", family="monospace")
        ax.text(track_lo - 0.005, pp3_y + bar_h / 2,
                "PP3", ha="right", va="center", fontsize=3.8,
                color="#7B3294", fontweight="bold")

        # PP4 eQTL bar (verdict colour)
        ax.add_patch(Rectangle((track_lo, pp4_y), track_w * pp4, bar_h,
                                facecolor=pp4_col, edgecolor="white",
                                linewidth=0.3, zorder=4))
        ax.text(track_hi + 0.003, pp4_y + bar_h / 2,
                f"{pp4:.2f}", ha="left", va="center", fontsize=4.4,
                color=pp4_col, fontweight="bold", family="monospace")
        ax.text(track_lo - 0.005, pp4_y + bar_h / 2,
                "PP4ᵉ", ha="right", va="center", fontsize=3.8,
                color=pp4_col, fontweight="bold")

        # PP4 sQTL bar (green, OR cross-hatched n/a if no sGene)
        if (not is_sgene) or pd.isna(sqtl_pp4):
            # Mark the bar as "no sGene" with a neutral hatched fill
            ax.add_patch(Rectangle((track_lo, sqtl_y), track_w, bar_h,
                                    facecolor="#FAFAFA",
                                    edgecolor="#CCCCCC",
                                    linewidth=0.3, hatch="//////",
                                    zorder=4))
            ax.text(track_hi + 0.003, sqtl_y + bar_h / 2,
                    "n/a", ha="left", va="center", fontsize=4.0,
                    color="#888888", fontweight="bold",
                    style="italic", family="monospace")
        else:
            ax.add_patch(Rectangle((track_lo, sqtl_y),
                                    track_w * sqtl_pp4, bar_h,
                                    facecolor=SQTL_COL, edgecolor="white",
                                    linewidth=0.3, zorder=4))
            ax.text(track_hi + 0.003, sqtl_y + bar_h / 2,
                    f"{sqtl_pp4:.2f}", ha="left", va="center", fontsize=4.4,
                    color=SQTL_COL, fontweight="bold", family="monospace")
        ax.text(track_lo - 0.005, sqtl_y + bar_h / 2,
                "PP4ˢ", ha="right", va="center", fontsize=3.8,
                color=SQTL_COL, fontweight="bold")

        # Verdict pill
        v = r["verdict"]
        v_w = col_w["verdict"] * 0.92
        v_h = 0.46
        v_x = col_x["verdict"] - v_w / 2
        _round_box(ax, v_x, y - v_h / 2, v_w, v_h,
                   fc=VERDICT_COL.get(v, "#BBBBBB"),
                   ec="white", lw=0.5, z=2, rounding=0.020)
        ax.text(col_x["verdict"], y, VERDICT_LBL.get(v, v),
                ha="center", va="center", fontsize=5.4,
                color="white", fontweight="bold", zorder=3)

        # Shape tag chip (small)
        s = r.get("shape_tag", "fail-other")
        s_col = SHAPE_COL.get(s, "#BDBDBD")
        s_lbl = SHAPE_LBL.get(s, s)
        chip_w, chip_h = col_w["shape"] * 0.78, 0.36
        chip_x = col_x["shape"] - chip_w / 2
        _round_box(ax, chip_x, y - chip_h / 2, chip_w, chip_h,
                   fc=s_col, ec="white", lw=0.4, z=2, rounding=0.012)
        s_txt_col = "white" if s in ("B", "C", "pass") else "#1A1A1A"
        ax.text(col_x["shape"], y, s_lbl,
                ha="center", va="center", fontsize=6.0,
                color=s_txt_col, fontweight="bold", zorder=3)

        # Family verdict
        fv = r.get("family_verdict", "NA")
        fv_lbl = {"pass": "pass", "fail": "fail",
                  "NA_underpowered": "n=1\nunder-pwr",
                  "NA": "—"}.get(fv, fv)
        fv_col_map = {"pass": "#1B7837", "fail": "#B2182B",
                      "NA_underpowered": "#9E9E9E", "NA": "#FFFFFF"}
        ax.text(col_x["famv"], y, fv_lbl,
                ha="center", va="center", fontsize=5.0,
                color=fv_col_map.get(fv, "#666666"),
                fontweight="bold", linespacing=1.1)

    # ---- bottom legend strip ----------------------------------------
    LEG_Y = 0.180
    # Plan C eligibility legend (replaces the old gate-glyph row)
    fig.text(0.020, LEG_Y, "Eligibility (Plan C)",
             fontsize=5.4, fontweight="bold", color=FG_COLOR,
             ha="left", va="center")
    items_elig = [
        ("eGene", "#1B7837", "PRIMARY: GTEx atrial eGene q<0.05"),
        ("Bonf",  "#1B7837", "SENSITIVITY: per-gene Bonf p < 0.05/n"),
        ("=",     "#1B7837", "primary == sensitivity (agreement)"),
        ("≠",     "#B2182B", "disagreement (see caption)"),
    ]
    x = 0.155
    for sym, col, lab in items_elig:
        if sym in ("eGene", "Bonf"):
            fig.patches.append(FancyBboxPatch(
                (x, LEG_Y - 0.011), 0.018, 0.020,
                boxstyle="round,pad=0.0,rounding_size=0.004",
                transform=fig.transFigure,
                facecolor=col, edgecolor="white", linewidth=0.4, zorder=10))
            fig.text(x + 0.009, LEG_Y, sym, fontsize=4.0,
                     color="white", fontweight="bold",
                     ha="center", va="center", zorder=11)
            fig.text(x + 0.022, LEG_Y, lab, fontsize=4.4,
                     color="#333333", ha="left", va="center")
        else:
            fig.text(x + 0.005, LEG_Y, sym, fontsize=7,
                     color=col, fontweight="bold",
                     ha="center", va="center")
            fig.text(x + 0.018, LEG_Y, lab, fontsize=4.4,
                     color="#333333", ha="left", va="center")
        x += 0.205

    # Verdict legend
    VERDICT_Y = LEG_Y - 0.030
    fig.text(0.020, VERDICT_Y, "Verdict",
             fontsize=5.4, fontweight="bold", color=FG_COLOR,
             ha="left", va="center")
    x = 0.085
    for v, col in VERDICT_COL.items():
        fig.patches.append(FancyBboxPatch(
            (x, VERDICT_Y - 0.011), 0.016, 0.020,
            boxstyle="round,pad=0.0,rounding_size=0.005",
            transform=fig.transFigure,
            facecolor=col, edgecolor="white", linewidth=0.4, zorder=10))
        fig.text(x + 0.020, VERDICT_Y, VERDICT_LBL[v], fontsize=4.6,
                 color="#333333", ha="left", va="center")
        x += 0.110

    # Coloc 3-bar legend (PP3 / PP4_eQTL / PP4_sQTL)
    COLOC_Y = LEG_Y - 0.060
    fig.text(0.020, COLOC_Y, "Coloc bars",
             fontsize=5.4, fontweight="bold", color=FG_COLOR,
             ha="left", va="center")
    items_coloc = [
        ("PP3",   "#7B3294",
         "PP3 = different causal variants (PITX2-mode)"),
        ("PP4ᵉ", "#1B7837",
         "PP4 from cis-eQTL (eQTL layer)"),
        ("PP4ˢ", "#1B7837",
         "PP4 from cis-sQTL leafcutter (sQTL layer)"),
    ]
    x = 0.085
    for sym, col, lab in items_coloc:
        fig.patches.append(FancyBboxPatch(
            (x, COLOC_Y - 0.011), 0.022, 0.020,
            boxstyle="round,pad=0.0,rounding_size=0.005",
            transform=fig.transFigure,
            facecolor=col, edgecolor="white", linewidth=0.4, zorder=10))
        fig.text(x + 0.011, COLOC_Y, sym, fontsize=4.4, color="white",
                 fontweight="bold", ha="center", va="center", zorder=11)
        fig.text(x + 0.025, COLOC_Y, lab, fontsize=4.0,
                 color="#333333", ha="left", va="center")
        x += 0.245

    # Shape tag legend
    SHAPE_Y = LEG_Y - 0.090
    fig.text(0.020, SHAPE_Y, "Shape tag",
             fontsize=5.4, fontweight="bold", color=FG_COLOR,
             ha="left", va="center")
    items_shape = [("B", "#7B3294", "PITX2-mode (PP3≥0.5, PP4<0.3)"),
                   ("C", "#A6611A", "eQTL absence (PP1+PP0≥0.7)"),
                   ("—†", "#BDBDBD", "B demoted (calibrator failed)")]
    x = 0.085
    for sym, col, lab in items_shape:
        fig.patches.append(FancyBboxPatch(
            (x, SHAPE_Y - 0.011), 0.020, 0.020,
            boxstyle="round,pad=0.0,rounding_size=0.005",
            transform=fig.transFigure,
            facecolor=col, edgecolor="white", linewidth=0.4, zorder=10))
        fig.text(x + 0.010, SHAPE_Y, sym, fontsize=4.6,
                 color=("white" if sym in ("B", "C") else "#1A1A1A"),
                 fontweight="bold",
                 ha="center", va="center", zorder=11)
        fig.text(x + 0.024, SHAPE_Y, lab, fontsize=4.4,
                 color="#333333", ha="left", va="center")
        x += 0.225

    # ---- title ------------------------------------------------------
    fig.suptitle(
        "10 anchors × cis-eQTL × AF GWAS coloc — none reach the "
        "single-gene mediation threshold",
        fontsize=7.4, y=0.965, fontweight="bold",
    )
    fig.text(0.012, 0.965, "b",
             fontsize=10, fontweight="bold", ha="left", va="top",
             transform=fig.transFigure)

    # ---- captions (locked, spec §10.12) -----------------------------
    n_B = int(summ.get("n_shape_B_atrial", 0))
    n_C = int(summ.get("n_shape_C_atrial", 0))
    n_other = (summ["n_anchors_total"] - summ["n_supported"] - n_B - n_C)
    n_primary = int(summ.get("planC_n_primary_eligible", 0))
    n_sens    = int(summ.get("planC_n_sensitivity_eligible", 0))
    n_agree   = int(summ.get("planC_n_agreement", 0))
    disagree  = str(summ.get("planC_disagreement_anchors", "(none)"))
    cap_main = (
        f"AF GWAS = GCST90624411 (Miyazawa 2023, ~1.6M); eQTL primary = "
        f"GTEx v8 atrial appendage (n=372); blood sensitivity = whole blood "
        f"(n=670). \u00b11 Mb cis window, coloc.abf default priors. "
        f"Figure-level scenario = "
        f"\u201c{summ['scenario']}\u201d "
        f"({summ['n_supported']}/{summ['n_anchors_total']} anchors supported, "
        f"{summ['n_fam_pass']}/3 multi-anchor families pass). "
        f"PITX2 calibrator: shape_tag = {summ['PITX2_shape_tag']} "
        f"(expected B; calibrator "
        f"{'PASS' if summ['PITX2_calibrator_validated'] else 'FAIL'})."
    )
    fig.text(0.500, 0.115, cap_main, ha="center", va="bottom",
             fontsize=4.6, color="#666666", style="italic")

    # Plan C eligibility (spec \u00a710.3 / \u00a710.12)
    cap_planC = (
        f"Eligibility (Plan C hybrid): PRIMARY = GTEx v8 atrial appendage "
        f"eGene status (q < 0.05, BH-FDR on permutation-adjusted p_beta); "
        f"SENSITIVITY = per-gene Bonferroni (lead p < 0.05 / n_variants in "
        f"cis window). The two gates AGREE on {n_agree}/{summ['n_anchors_total']} "
        f"anchors ({n_primary}/{summ['n_anchors_total']} eGene-eligible vs "
        f"{n_sens}/{summ['n_anchors_total']} Bonferroni-eligible; "
        f"disagreements: {disagree}). Verdict driven by primary; the legacy "
        f"5e-8 GWAS-significance threshold is methodologically misapplied to "
        f"cis-eQTL discovery and is reported only as `gate_5e8_legacy` for "
        f"transparency \u2014 never as a gate in v3."
    )
    fig.text(0.500, 0.092, cap_planC, ha="center", va="bottom",
             fontsize=4.6, color="#666666", style="italic")

    # Failure-mode description (spec §10.12 / §10.14)
    n_sGene = int(summ.get("n_sGene_atrial", 0))
    n_total = int(summ.get("n_anchors_total", 10))
    cap_failure = (
        f"3-bar coloc: PP3 (different causal variants), PP4ᵉ from cis-eQTL, "
        f"PP4ˢ from cis-sQTL (leafcutter, atrial appendage). "
        f"sQTL layer is structurally inapplicable: only {n_sGene}/{n_total} "
        f"anchors are atrial sGenes (q<0.05; spec §10.14), so the PP4ˢ bar "
        f"is hatched n/a for the rest. Combined RNA-layer verdict: "
        f"{int(summ.get('n_supp_either_RNA_layer', 0))}/{n_total} anchors "
        f"supported across either layer. "
        f"Among failed anchors, n_B = {n_B} reproduce PITX2-mode "
        f"(B; demoted post-calibrator-fail), n_C = {n_C} show eQTL absence; "
        f"the remaining {int(n_other)} are intermediate fail-other."
    )
    fig.text(0.500, 0.045, cap_failure, ha="center", va="bottom",
             fontsize=4.6, color="#666666", style="italic")

    # RA-biased clause (spec §10.12, paragraph 1)
    cap_RA = (
        "Both v3 extension families (F11 BMP, F15 NTRK2/TrkB) are "
        "RA-biased in the healthy chamber prior (Fig 1f); positive coloc "
        "here would mark RA-side individual drivers, complementing the "
        "LA-side ECM (F10/PRELP/OMD) v2 set. This LA-RA symmetry across "
        "anchor cohorts strengthens the distributed-substrate interpretation."
    )
    fig.text(0.500, 0.015, cap_RA, ha="center", va="bottom",
             fontsize=4.6, color="#666666", style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
