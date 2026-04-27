#!/usr/bin/env python3
"""
Figure 6C — Pre-registered anchor evidence matrix.

Reads `results/figure6_mr_anchors.tsv` (primary, GTEx atrial appendage)
and `results/figure6_mr_anchors_blood.tsv` (sensitivity, GTEx whole blood)
and renders a 7-column evidence matrix per anchor:

  Anchor | Axis | cis-eQTL primary | cis-eQTL blood | AF GWAS | coloc PP4 | Verdict

Verdict bins (locked colocalisation thresholds):
  supported     #1B7837 (green)
  inconclusive  #9E9E9E (grey)
  not supported #B2182B (red)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

from plot_style import COL_DOUBLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "results" / "figure6_mr_anchors.tsv"
IN_BLOOD = ROOT / "results" / "figure6_mr_anchors_blood.tsv"
IN_POSITIVE = ROOT / "results" / "figure6_mr_positive_controls.tsv"
OUT = ROOT / "figures" / "fig6c_mr_plausibility.pdf"

VERDICT_COL = {
    "supported":     "#1B7837",
    "inconclusive":  "#9E9E9E",
    "not supported": "#B2182B",
}
AXIS_COL = {
    "matrix":           "#B35806",
    "metabolic":        "#542788",
    "positive_control": "#444444",
}

ANCHOR_ORDER = ["PRELP", "OMD", "IDH2", "NDUFA8"]
POSITIVE_CONTROL_ORDER = ["PITX2", "KCNN3"]


def fmt_p(p):
    if p is None or pd.isna(p):
        return "n/a"
    if p == 0:
        return "<1e-300"
    if p < 1e-3:
        s = f"{p:.1e}".replace("e-0", "e-").replace("e+0", "e+")
        return s
    return f"{p:.3f}"


def gate_glyph(p_eqtl, threshold=5e-8, marginal=1e-5):
    if pd.isna(p_eqtl):
        return ("", "#888888")
    if p_eqtl < threshold:
        return ("●", "#1B7837")
    if p_eqtl < marginal:
        return ("◐", "#888888")
    return ("○", "#B2182B")


def coloc_color(pp4):
    if pd.isna(pp4):
        return "#F2F2F2"
    if pp4 >= 0.7:
        return "#A6DBA0"
    if pp4 >= 0.3:
        return "#E0E0E0"
    return "#F4ACA0"


def _round_box(ax, x, y, w, h, fc, ec="white", lw=0.4, z=2, rounding=0.012):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    ))


def _load_positive_controls() -> pd.DataFrame:
    """Load PITX2 / KCNN3 positive-control results and reshape to one row
    per anchor with both atrial-appendage and blood columns."""
    if not IN_POSITIVE.exists():
        return pd.DataFrame()
    pc = pd.read_csv(IN_POSITIVE, sep="\t")
    primary = (pc.query("tissue == 'atrial_appendage'")
                 .set_index("anchor")[["axis", "lead_eqtl_p",
                                        "gwas_min_p", "coloc_PP4",
                                        "verdict"]])
    blood = (pc.query("tissue == 'whole_blood'")
               .set_index("anchor")[["lead_eqtl_p", "coloc_PP4"]]
               .rename(columns={"lead_eqtl_p": "lead_eqtl_p_blood",
                                "coloc_PP4": "coloc_PP4_blood"}))
    out = primary.join(blood, how="left").reset_index()
    out["is_positive_control"] = True
    out["axis"] = "positive_control"
    out["blood_verdict"] = pd.NA
    out["gate_passed_blood"] = pd.NA
    return out


def main() -> None:
    set_nature_style()
    df = pd.read_csv(IN, sep="\t").set_index("anchor").loc[ANCHOR_ORDER].reset_index()
    if IN_BLOOD.exists():
        bl = pd.read_csv(IN_BLOOD, sep="\t").rename(columns={
            "lead_eqtl_p": "lead_eqtl_p_blood",
            "coloc_PP4":  "coloc_PP4_blood",
            "gate_passed": "gate_passed_blood",
            "blood_verdict": "blood_verdict",
        })[["anchor", "lead_eqtl_p_blood", "coloc_PP4_blood",
            "gate_passed_blood", "blood_verdict"]]
        df = df.merge(bl, on="anchor", how="left")
    else:
        df["lead_eqtl_p_blood"] = pd.NA
        df["coloc_PP4_blood"] = pd.NA
        df["blood_verdict"] = pd.NA
    df["is_positive_control"] = False

    # ---- positive controls (PITX2, KCNN3) ----------------------------------
    pc = _load_positive_controls()
    if not pc.empty:
        pc = pc.set_index("anchor").loc[
            [a for a in POSITIVE_CONTROL_ORDER if a in pc.index]
        ].reset_index()
        # align columns to the main df
        for col in df.columns:
            if col not in pc.columns:
                pc[col] = pd.NA
        df = pd.concat([df, pc[df.columns]], ignore_index=True)
    n_anchor = len(ANCHOR_ORDER)

    # ---- column geometry --------------------------------------------------
    columns = [
        ("anchor",     "Anchor",                 0.11),
        ("axis",       "Axis",                   0.11),
        ("eqtl",       "cis-eQTL primary\nGTEx atrial appendage",   0.20),
        ("eqtl_blood", "cis-eQTL sensitivity\nGTEx whole blood",    0.20),
        ("gwas",       "AF GWAS\nmin p (window)",                   0.13),
        ("coloc",      "coloc PP4\n(primary)",                      0.11),
        ("verdict",    "Verdict",                                   0.14),
    ]
    n_row = len(df)
    fig_w = COL_DOUBLE
    row_h = 0.50
    fig_h = row_h * n_row + 1.85
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.025, right=0.975, top=0.78, bottom=0.20)
    ax.set_axis_off()

    xs = np.cumsum([0.0] + [c[2] for c in columns])
    x_centers = (xs[:-1] + xs[1:]) / 2
    col_x = dict(zip([c[0] for c in columns], x_centers))
    col_w = dict(zip([c[0] for c in columns], [c[2] for c in columns]))

    # ---- header band ------------------------------------------------------
    HEAD_TOP = n_row + 0.20
    HEAD_BOT = n_row - 0.40
    HEAD_Y = (HEAD_TOP + HEAD_BOT) / 2
    _round_box(ax, xs[0], HEAD_BOT, xs[-1] - xs[0], HEAD_TOP - HEAD_BOT,
               fc="#F0EEEA", ec="none", z=1, rounding=0.014)
    for (key, label, w), x in zip(columns, x_centers):
        ax.text(x, HEAD_Y, label, ha="center", va="center",
                fontsize=6.0, fontweight="bold", color="#222222",
                linespacing=1.25)

    # ---- rows -------------------------------------------------------------
    ROW_HALF = 0.42
    for i, r in df.iterrows():
        y = n_row - 1 - i
        # alternating subtle stripe
        if i % 2 == 1:
            ax.add_patch(Rectangle((xs[0], y - ROW_HALF),
                                   xs[-1] - xs[0], 2 * ROW_HALF,
                                   facecolor="#F8F8F8", edgecolor="none",
                                   zorder=0))
        # bold divider between anchors and positive controls
        if i == n_anchor:
            ax.plot([xs[0] + 0.005, xs[-1] - 0.005],
                    [y + ROW_HALF + 0.02, y + ROW_HALF + 0.02],
                    color="#666666", lw=1.0, zorder=0)
            ax.text(xs[0] + 0.012, y + ROW_HALF - 0.02,
                    "Pre-registered POSITIVE CONTROLS — canonical AF GWAS loci",
                    fontsize=4.6, fontweight="bold", color="#444444",
                    style="italic", ha="left", va="top", zorder=3)
        else:
            # row separator (very faint)
            ax.plot([xs[0] + 0.005, xs[-1] - 0.005],
                    [y - ROW_HALF, y - ROW_HALF],
                    color="#E2E2E2", lw=0.4, zorder=0)

        # Anchor name
        ax.text(col_x["anchor"], y, r["anchor"],
                ha="center", va="center", fontsize=7.0, fontweight="bold",
                color="#1A1A1A")

        # Axis chip (rounded)
        chip_w, chip_h = 0.075, 0.34
        _round_box(ax, col_x["axis"] - chip_w / 2, y - chip_h / 2,
                   chip_w, chip_h, fc=AXIS_COL[r["axis"]],
                   ec="white", lw=0.4, z=2, rounding=0.012)
        ax.text(col_x["axis"], y, r["axis"],
                ha="center", va="center", fontsize=5.6, color="white",
                fontweight="bold", zorder=3)

        # cis-eQTL primary: glyph above, monospace p-value below
        glyph, glyph_col = gate_glyph(r["lead_eqtl_p"])
        ax.text(col_x["eqtl"], y + 0.13, glyph,
                ha="center", va="center", fontsize=10.5,
                color=glyph_col, fontweight="bold")
        ax.text(col_x["eqtl"], y - 0.18,
                f"lead p = {fmt_p(r['lead_eqtl_p'])}",
                ha="center", va="center", fontsize=5.4, color="#444444",
                family="monospace")

        # cis-eQTL blood (sensitivity)
        if pd.isna(r.get("lead_eqtl_p_blood")):
            ax.text(col_x["eqtl_blood"], y + 0.13, "—",
                    ha="center", va="center", fontsize=10,
                    color="#AAAAAA", fontweight="bold")
            ax.text(col_x["eqtl_blood"], y - 0.18,
                    "n/a (gene absent in tissue)",
                    ha="center", va="center", fontsize=5.0,
                    color="#888888", style="italic")
        else:
            glyph_b, glyph_col_b = gate_glyph(r["lead_eqtl_p_blood"])
            ax.text(col_x["eqtl_blood"], y + 0.13, glyph_b,
                    ha="center", va="center", fontsize=10.5,
                    color=glyph_col_b, fontweight="bold")
            ax.text(col_x["eqtl_blood"], y - 0.18,
                    f"lead p = {fmt_p(r['lead_eqtl_p_blood'])}",
                    ha="center", va="center", fontsize=5.4, color="#444444",
                    family="monospace")

        # GWAS min p
        ax.text(col_x["gwas"], y, fmt_p(r["gwas_min_p"]),
                ha="center", va="center", fontsize=6.0,
                color="#222222", family="monospace", fontweight="bold")

        # Coloc PP4 (rounded colored cell)
        cell_w, cell_h = 0.078, 0.36
        _round_box(ax, col_x["coloc"] - cell_w / 2, y - cell_h / 2,
                   cell_w, cell_h, fc=coloc_color(r["coloc_PP4"]),
                   ec="#CCCCCC", lw=0.4, z=2, rounding=0.010)
        pp_label = (f"{r['coloc_PP4']:.2f}" if not pd.isna(r["coloc_PP4"])
                    else "n/a")
        ax.text(col_x["coloc"], y, pp_label,
                ha="center", va="center", fontsize=6.4,
                color="#1A1A1A", fontweight="bold", zorder=3)

        # Verdict pill (rounded)
        v = r["verdict"]
        v_w, v_h = 0.115, 0.42
        _round_box(ax, col_x["verdict"] - v_w / 2, y - v_h / 2,
                   v_w, v_h, fc=VERDICT_COL[v],
                   ec="white", lw=0.5, z=2, rounding=0.020)
        ax.text(col_x["verdict"], y, v,
                ha="center", va="center", fontsize=5.6, color="white",
                fontweight="bold", zorder=3)

    # outer rounded border
    _round_box(ax, xs[0], -ROW_HALF - 0.05, xs[-1] - xs[0],
               n_row + 0.65, fc="none", ec="#BBBBBB", lw=0.45, z=4,
               rounding=0.014)

    # ---- bottom legend strip ---------------------------------------------
    LEG_Y = 0.105
    fig.text(0.045, LEG_Y, "cis-eQTL gate",
             fontsize=5.6, fontweight="bold", color="#222222",
             ha="left", va="center")
    items = [("●", "#1B7837", "p < 5e-8"),
             ("◐", "#888888", "5e-8 ≤ p < 1e-5  marginal"),
             ("○", "#B2182B", "p ≥ 1e-5  gate fail")]
    x = 0.16
    for sym, col, lab in items:
        fig.text(x, LEG_Y, sym, fontsize=9, color=col,
                 ha="left", va="center", fontweight="bold")
        fig.text(x + 0.018, LEG_Y, lab, fontsize=5.0, color="#333333",
                 ha="left", va="center")
        x += 0.13

    fig.text(0.575, LEG_Y, "Verdict",
             fontsize=5.6, fontweight="bold", color="#222222",
             ha="left", va="center")
    x = 0.625
    for v, col in VERDICT_COL.items():
        fig.patches.append(FancyBboxPatch(
            (x, LEG_Y - 0.012), 0.018, 0.024,
            boxstyle="round,pad=0.0,rounding_size=0.005",
            transform=fig.transFigure,
            facecolor=col, edgecolor="white", linewidth=0.4,
        ))
        fig.text(x + 0.022, LEG_Y, v, fontsize=5.0, color="#333333",
                 ha="left", va="center")
        x += 0.10

    # ---- caption ----------------------------------------------------------
    fig.text(0.5, 0.025,
             "Plausibility — not causal.  cis-eQTL × AF-GWAS coloc on a "
             "pre-registered anchor set; primary tissue = GTEx atrial "
             "appendage (n = 372); blood sensitivity = GTEx whole blood "
             "(n = 670, eQTL Catalogue substitute for eQTLGen).  Two "
             "canonical AF GWAS loci (PITX2, KCNN3) are shown below the "
             "divider as pre-registered positive controls; both fail the "
             "same conservative gate at this panel power, indicating "
             "that the framework anchor null reflects the calibration of "
             "the cis-eQTL panel rather than absence of any single-gene "
             "mediation in principle.  Verdict bins fixed in spec §6.7; "
             "gate = lead cis-eQTL p < 5e-8 within ±1 Mb of the gene.",
             ha="center", va="bottom", fontsize=4.6,
             color="#666666", style="italic")

    fig.suptitle(
        "Pre-registered anchors do not pass the conservative cis-eQTL "
        "gate at current atrial-appendage panel power",
        fontsize=7.4, y=0.965, fontweight="bold",
    )
    fig.text(0.012, 0.965, "c",
             fontsize=10, fontweight="bold", ha="left", va="top",
             transform=fig.transFigure)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
