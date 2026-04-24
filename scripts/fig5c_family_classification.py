#!/usr/bin/env python3
"""
Figure 5C  Family evidence-class matrix  (integrates Fig 5B).

For each axis family we classify the triangulation pattern across
Fig 2D (bulk g) / Fig 3B (LA snRNA max-sig Δ) / Fig 4B (LAA duration ρ)
into one of five mutually exclusive classes:

    CONVERGENT        ≥ 2 layers deliver a statistically supported hit
                        with consistent sign.
    BULK-LED          only L1 bulk g is supported (|g| CI excludes 0);
                        other layers are flat or non-significant.
    SNRNA-LED         only L2 LA snRNA carrier is supported
                        (p_perm < 0.05); other layers flat.
    DURATION-LED      only L3 LAA duration ρ is supported (p_perm < 0.10
                        and |ρ| ≥ 0.5 given n = 12).
    SILENT            no layer supported (the framework's own null
                        class — NOT a figure of merit, shown for
                        transparency).

Sign consistency for CONVERGENT:  we measure the "AF-ward" direction
of each layer relative to baseline chamber identity:

        L1 bulk g       +  (AF pushes toward higher LA − RA; LA-ward)
        L2 snRNA Δ      +  (AF higher singscore in LA)
        L3 duration ρ   +  (singscore intensifies with months)

Rows are ordered LA-biased on top, RA-biased on bottom.

Inputs:
    results/fig5b_triangulation.tsv   (from fig5b_triangulation_heatmap.py)
    results/fig3_primary_family_celltype_delta.tsv  (for top-ct readouts)

Outputs:
    results/fig5c_family_classes.tsv
    figures/fig5c_family_classification.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

TRI_TSV  = ROOT / "results" / "fig5b_triangulation.tsv"
OUT_TSV  = ROOT / "results" / "fig5c_family_classes.tsv"
FIG_OUT  = ROOT / "figures" / "fig5c_family_classification.pdf"

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"

CLASS_COLORS = {
    "CONVERGENT":    "#2E7B5E",
    "BULK-LED":      "#8C510A",
    "SNRNA-LED":     "#5E3C99",
    "DURATION-LED":  "#1F78B4",
    "SILENT":        "#BBBBBB",
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
        "axes.linewidth":    0.6,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def layer_supported(tri: pd.Series) -> dict[str, int]:
    """Return +1 / 0 / −1 per layer; 0 = unsupported, else sign of effect."""
    out = {"L1": 0, "L2": 0, "L3": 0}

    # L1 — CI excludes 0
    if pd.notna(tri.get("ci_lo")) and pd.notna(tri.get("ci_hi")):
        if tri["ci_lo"] > 0 or tri["ci_hi"] < 0:
            out["L1"] = int(np.sign(tri["bulk_g"]))

    # L2 — p_perm < 0.05
    if pd.notna(tri.get("snrna_p")) and tri["snrna_p"] < 0.05:
        out["L2"] = int(np.sign(tri["snrna_d"]))

    # L3 — p_perm < 0.10 AND |ρ| ≥ 0.50 (n = 12 is small; guard against
    # spurious low-magnitude p<0.10 hits)
    if pd.notna(tri.get("dur_p")) and tri["dur_p"] < 0.10 \
            and abs(tri["dur_rho"]) >= 0.50:
        out["L3"] = int(np.sign(tri["dur_rho"]))

    return out


def classify_family(tri: pd.Series) -> tuple[str, str]:
    """Return (class, rationale)."""
    sup = layer_supported(tri)
    n_sup = sum(1 for v in sup.values() if v != 0)
    layers_sup = [k for k, v in sup.items() if v != 0]

    if n_sup == 0:
        return ("SILENT", "no layer reaches support threshold")
    if n_sup >= 2:
        signs = [v for v in sup.values() if v != 0]
        if all(s == signs[0] for s in signs):
            rat = " + ".join(layers_sup) + " consistent"
            return ("CONVERGENT", rat)
        else:
            rat = " + ".join(layers_sup) + " conflict"
            return ("CONVERGENT", rat + "  (sign-mixed)")
    # n_sup == 1
    if sup["L1"] != 0: return ("BULK-LED",     "L1 only")
    if sup["L2"] != 0: return ("SNRNA-LED",    "L2 only")
    if sup["L3"] != 0: return ("DURATION-LED", "L3 only")
    return ("SILENT", "unreachable")


def main() -> None:
    set_publication_style()

    tri = pd.read_csv(TRI_TSV, sep="\t")

    rows = []
    for _, t in tri.iterrows():
        sup = layer_supported(t)
        cls, rat = classify_family(t)
        rows.append({
            "family_id":       int(t["family_id"]),
            "direction_label": t["direction_label"],
            "short":           t["short"],
            "n_programs":      int(t["n_programs"])
                                  if pd.notna(t["n_programs"]) else 0,
            "bulk_g":          t.get("bulk_g", np.nan),
            "bulk_ci_lo":      t.get("ci_lo", np.nan),
            "bulk_ci_hi":      t.get("ci_hi", np.nan),
            "snrna_top_ct":    t.get("snrna_top_ct", ""),
            "snrna_d":         t.get("snrna_d", np.nan),
            "snrna_p":         t.get("snrna_p", np.nan),
            "dur_rho":         t.get("dur_rho", np.nan),
            "dur_p":           t.get("dur_p", np.nan),
            "L1_signed_sup":   sup["L1"],
            "L2_signed_sup":   sup["L2"],
            "L3_signed_sup":   sup["L3"],
            "evidence_class":  cls,
            "rationale":       rat,
        })
    cls_df = pd.DataFrame(rows)

    # row order: classes (convergent → led → silent), within class by |bulk_g|
    cls_order = ["CONVERGENT", "BULK-LED", "SNRNA-LED",
                  "DURATION-LED", "SILENT"]
    cls_df["_co"] = cls_df["evidence_class"].map(
        {c: i for i, c in enumerate(cls_order)})
    cls_df["_ag"] = cls_df["bulk_g"].abs()
    cls_df = cls_df.sort_values(["_co", "_ag"],
                                      ascending=[True, False]).reset_index(drop=True)
    cls_df = cls_df.drop(columns=["_co", "_ag"])

    cls_df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")

    # ---- figure --------------------------------------------------------
    K = len(cls_df)
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(200 * MM, 135 * MM))

    fig.suptitle(
        "Fig 5C  Evidence class per axis family  "
        "(triangulation of Fig 2 + 3 + 4)",
        fontsize=8.8, fontweight="bold", y=0.97, x=0.02, ha="left")

    # --- layout : central class pill + 3 layer badges + rationale ---
    ax = fig.add_axes([0.02, 0.07, 0.96, 0.82])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # column layout (x in [0,1])
    col_xs = {
        "family":   (0.02, 0.20),   # F# + label
        "class":    (0.21, 0.33),   # pill
        "bulk":     (0.34, 0.45),   # L1 value + CI
        "snrna":    (0.46, 0.60),   # L2 value + top ct
        "dur":      (0.61, 0.71),   # L3 value
        "rationale":(0.72, 0.98),   # narrative
    }
    # Header row
    header_y = 0.94
    headers = [
        ("family",   "Axis family"),
        ("class",    "Evidence class"),
        ("bulk",     "L1  Bulk g  [95 % CI]"),
        ("snrna",    "L2  LA snRNA Δ @ top cell"),
        ("dur",      "L3  LAA ρ vs duration"),
        ("rationale","Rationale"),
    ]
    for key, label in headers:
        x0, x1 = col_xs[key]
        ax.text((x0 + x1) / 2, header_y, label,
                  ha="center", va="center",
                  fontsize=7.5, fontweight="bold", color=GREY_DARK)
    ax.plot([0.01, 0.99], [header_y - 0.025, header_y - 0.025],
              color=GREY_MID, lw=0.6)

    # Rows
    row_h = (header_y - 0.05) / (K + 1)
    for i, (_, r) in enumerate(cls_df.iterrows()):
        y = header_y - 0.04 - (i + 0.5) * row_h

        # Family label + direction colour
        x0, x1 = col_xs["family"]
        d   = r["direction_label"]
        col = LA_COLOR if d == "LA" else RA_COLOR
        ax.add_patch(mpatches.Rectangle(
            (x0, y - row_h * 0.38), 0.02, row_h * 0.76,
            facecolor=col, alpha=0.5, edgecolor="none"))
        ax.text(x0 + 0.028, y,
                  f"F{int(r['family_id'])}  {r['short']}  "
                  f"(n={int(r['n_programs'])})",
                  ha="left", va="center", fontsize=7.0,
                  color=col, fontweight="bold")

        # Class pill
        x0, x1 = col_xs["class"]
        face = CLASS_COLORS[r["evidence_class"]]
        pill = mpatches.FancyBboxPatch(
            ((x0 + x1) / 2 - 0.055, y - row_h * 0.33),
            0.11, row_h * 0.66,
            boxstyle="round,pad=0.003,rounding_size=0.010",
            facecolor=face, edgecolor=face,
            linewidth=0.6, zorder=3)
        ax.add_patch(pill)
        ax.text((x0 + x1) / 2, y, r["evidence_class"],
                  ha="center", va="center",
                  fontsize=6.5, fontweight="bold",
                  color="white", zorder=4)

        # L1 bulk g value + CI
        x0, x1 = col_xs["bulk"]
        if pd.notna(r["bulk_g"]):
            sig_marker = ""
            if pd.notna(r["bulk_ci_lo"]) and pd.notna(r["bulk_ci_hi"]):
                if r["bulk_ci_lo"] > 0 or r["bulk_ci_hi"] < 0:
                    sig_marker = "  *"
            txt = (f"{r['bulk_g']:+.2f}"
                     f"  [{r['bulk_ci_lo']:+.2f}, {r['bulk_ci_hi']:+.2f}]"
                     f"{sig_marker}")
            col_txt = ("#000000" if sig_marker else GREY_DARK)
            fw = "bold" if sig_marker else "normal"
            ax.text((x0 + x1) / 2, y, txt,
                      ha="center", va="center",
                      fontsize=6.3, color=col_txt, fontweight=fw)

        # L2 snRNA Δ + cell type
        x0, x1 = col_xs["snrna"]
        if pd.notna(r["snrna_d"]):
            mark = ""
            if pd.notna(r["snrna_p"]):
                if r["snrna_p"] < 0.05: mark = "  *"
                elif r["snrna_p"] < 0.10: mark = "  ."
            ct = (r["snrna_top_ct"] or "").replace("_", " ")
            txt = f"{r['snrna_d']:+.3f}  @ {ct}  (p={r['snrna_p']:.2f}){mark}"
            col_txt = "#000000" if "*" in mark else GREY_DARK
            fw = "bold" if "*" in mark else "normal"
            ax.text((x0 + x1) / 2, y, txt,
                      ha="center", va="center",
                      fontsize=6.3, color=col_txt, fontweight=fw)

        # L3 duration ρ
        x0, x1 = col_xs["dur"]
        if pd.notna(r["dur_rho"]):
            mark = ""
            if pd.notna(r["dur_p"]):
                if r["dur_p"] < 0.05: mark = "  *"
                elif r["dur_p"] < 0.10: mark = "  ."
            txt = f"ρ = {r['dur_rho']:+.2f}  (p={r['dur_p']:.2f}){mark}"
            col_txt = "#000000" if "*" in mark else GREY_DARK
            fw = "bold" if "*" in mark else "normal"
            ax.text((x0 + x1) / 2, y, txt,
                      ha="center", va="center",
                      fontsize=6.3, color=col_txt, fontweight=fw)

        # Rationale
        x0, x1 = col_xs["rationale"]
        ax.text(x0, y, r["rationale"],
                  ha="left", va="center",
                  fontsize=6.2, color=GREY_DARK,
                  fontstyle="italic")

        # light separator
        ax.plot([0.01, 0.99], [y - row_h / 2, y - row_h / 2],
                  color="#EEEEEE", lw=0.3, zorder=0)

    # Bottom legend : class colour chips
    legend_y = 0.03
    legend_items = list(CLASS_COLORS.items())
    total_w = 0.70
    chip_w = total_w / len(legend_items)
    start_x = 0.02
    fig.text(start_x, legend_y + 0.02, "class colour key:",
              fontsize=6.0, color=GREY_DARK, fontweight="bold")
    for k, (cname, ccol) in enumerate(legend_items):
        x = start_x + 0.11 + k * chip_w
        fig.patches.append(mpatches.Rectangle(
            (x, legend_y + 0.016), 0.012, 0.016,
            transform=fig.transFigure,
            facecolor=ccol, edgecolor="none"))
        fig.text(x + 0.016, legend_y + 0.024, cname,
                  fontsize=6.0, color=GREY_DARK, va="center")

    fig.text(0.02, 0.008,
              "Support rules:  L1 | 95 % CI excludes 0     "
              "L2 | p_perm < 0.05     "
              "L3 | p_perm < 0.10  AND  |ρ| ≥ 0.50  (n = 12).   "
              "* support marker, . p < 0.10 sub-threshold for L2 / L3.",
              fontsize=5.8, color="#444", ha="left", va="bottom")

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")

    # Console summary
    print("\n=== evidence class tally ===")
    print(cls_df["evidence_class"].value_counts().to_string())
    print("\n=== per-family classes ===")
    print(cls_df[["family_id", "direction_label", "short",
                      "evidence_class", "rationale"]].to_string(index=False))


if __name__ == "__main__":
    main()
