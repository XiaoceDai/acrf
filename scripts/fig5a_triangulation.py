#!/usr/bin/env python3
"""Figure 5A — Program triangulation heatmap."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

from plot_style import (
    COL_DOUBLE,
    COMPARTMENT_COLORS,
    DURATION_CMAP,
    LA_COLOR,
    RA_COLOR,
    REWIRING_CLASS_COLORS,
    set_nature_style,
)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "fig5a_triangulation.pdf"

PANEL = ROOT / "results" / "figure3_program_panel.tsv"
CARRIER = ROOT / "results" / "figure3_carrier_shift.tsv"
DURATION = ROOT / "Project_P4" / "GSE41177" / "results" / "duration_trend_laa_primary.tsv"


CLASS_COLORS = {
    "Convergent":      "#1B7837",
    "Structure-led":   "#66A61E",
    "Metabolism-led":  "#E6AB02",
    "Fragmented":      "#9E9E9E",
}


def bulk_score(consensus, trend128):
    """Map class string to a signed strength in {-2,-1,0,+1,+2}."""
    label = consensus if isinstance(consensus, str) and consensus and consensus != '""' else trend128
    label = (label or "").strip().strip('"')
    table = {
        "Reversed":         -2, "Reversed-weak":   -1,
        "Attenuated":       -1, "Attenuated-weak": -1,
        "Preserved":         0, "":                 0, "Low-support": 0,
        "Amplified-weak":   +1, "Amplified":       +2,
    }
    return table.get(label, 0), label or "n/a"


def dur_dir(rho, thr=0.30):
    if pd.isna(rho):
        return 0
    if rho <= -thr:
        return -1
    if rho >= +thr:
        return +1
    return 0


def triangulate(row):
    bulk = row["bulk_strength"]
    dd = row["dur_dir"]
    is_metab = row["baseline_direction"] == 1
    is_struct = row["baseline_direction"] == -1
    loc = (row["AF_major_carrier_CT"] or "").lower()

    # Metabolism-led
    if is_metab and dd == -1:
        return "Metabolism-led"

    # Structure-led
    if is_struct and bulk <= -1:
        if "fibroblast" in loc:
            return "Structure-led"
        # Allow fallback: healthy_top_expressor was fibroblast (Fig 1 expressor)
        if "fibroblast" in (row.get("healthy_top_expressor_CT") or "").lower():
            return "Structure-led"

    # Convergent: bulk and duration agree on weakening/decline
    if bulk != 0 and dd != 0:
        # For LA-high programs, AF-side should decline (bulk -1/-2 and dur -1)
        # For RA-high programs (in LAA arm), AF score may rise as identity blurs (bulk -1/-2 and dur +1).
        if (is_metab and bulk <= -1 and dd == -1) or (is_struct and bulk <= -1 and dd == +1):
            return "Convergent"

    return "Fragmented"


def main() -> None:
    set_nature_style()

    panel = pd.read_csv(PANEL, sep="\t")
    carrier = pd.read_csv(CARRIER, sep="\t")[
        ["program_id", "AF_major_carrier_CT", "healthy_top_expressor_CT"]
    ]
    dur = pd.read_csv(DURATION, sep="\t")[["program_id", "statistic", "p_perm"]]
    dur = dur.rename(columns={"statistic": "rho"})

    df = panel.merge(carrier, on="program_id", how="left").merge(dur, on="program_id", how="left")
    df = df.sort_values("program_order").reset_index(drop=True)

    parsed = df.apply(lambda r: bulk_score(r.get("consensus_class"),
                                           r.get("trend_class_GSE128188")), axis=1)
    df["bulk_strength"] = [p[0] for p in parsed]
    df["bulk_label"] = [p[1] for p in parsed]
    df["dur_dir"] = df["rho"].apply(dur_dir)
    df["AF_major_carrier_CT"] = df["AF_major_carrier_CT"].fillna("Unresolved")
    df["triangulation"] = df.apply(triangulate, axis=1)

    n = len(df)
    fig_w = COL_DOUBLE * 0.95
    fig_h = 0.40 * n + 1.85
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.27, right=0.97, top=0.80, bottom=0.16)

    cmap_dur = LinearSegmentedColormap.from_list("dur", DURATION_CMAP, N=256)
    norm_dur = TwoSlopeNorm(vcenter=0.0, vmin=-0.6, vmax=+0.6)

    col_x = {
        "Healthy bias":       0,
        "Bulk rewiring":      1,
        "LA AF carrier":      2,
        "Duration ρ (LAA)":   3,
        "Triangulation":      4,
    }

    def cell(j, i, face, edge="white", lw=0.6):
        ax.add_patch(Rectangle((j - 0.46, i - 0.46), 0.92, 0.92,
                               facecolor=face, edgecolor=edge, linewidth=lw,
                               zorder=2))

    def text(j, i, s, color="#222222", weight="normal", size=5.7):
        ax.text(j, i, s, ha="center", va="center", fontsize=size,
                color=color, fontweight=weight, zorder=3)

    # ROW LOOP
    for i, row in df.iterrows():
        # 1. Healthy bias
        face = LA_COLOR if row["baseline_direction"] == 1 else RA_COLOR
        cell(0, i, face)
        text(0, i, "LA" if row["baseline_direction"] == 1 else "RA",
             color="white", weight="bold")

        # 2. Bulk rewiring
        lab = row["bulk_label"] if row["bulk_label"] in REWIRING_CLASS_COLORS \
            else "Preserved"
        face = REWIRING_CLASS_COLORS.get(lab, "#FFFFFF")
        cell(1, i, face)
        # short label keeps "Attenuated-weak" readable at small size
        short_bulk = {
            "Amplified": "Amplified",
            "Amplified-weak": "Ampl-w",
            "Preserved": "Preserved",
            "Attenuated-weak": "Atten-w",
            "Attenuated": "Attenuated",
            "Reversed-weak": "Rev-w",
            "Reversed": "Reversed",
            "Low-support": "Low-supp",
            "NA": "n/a",
        }.get(lab, lab)
        text(1, i, short_bulk, color="#222222")

        # 3. LA AF carrier — when all rows are Unresolved we keep a uniform
        # grey tile (no per-row text) and explain in the caption.  Otherwise
        # we render the carrier name on a coloured tile.
        loc = row["AF_major_carrier_CT"]
        if not isinstance(loc, str) or loc == "" or loc == "Unresolved":
            cell(2, i, "#E5E5E5", edge="white")
            # no per-cell text; caption note explains
        else:
            face = COMPARTMENT_COLORS.get(loc, COMPARTMENT_COLORS["unresolved"])
            short = {"atrialCM": "atrialCM", "fibroblast": "fibroblast",
                     "endothelial_or_vascular": "endo/vasc",
                     "immune": "immune"}.get(loc, loc)
            cell(2, i, face, edge="white")
            text(2, i, short, color="white")

        # 4. Duration ρ
        rho = row["rho"]
        if pd.isna(rho):
            cell(3, i, "#F2F2F2"); text(3, i, "n/a", color="#777777")
        else:
            face = cmap_dur(norm_dur(rho))
            cell(3, i, face)
            txtc = "white" if abs(norm_dur(rho) - 0.5) > 0.35 else "#222222"
            text(3, i, f"{rho:+.2f}", color=txtc, weight="bold")

        # 5. Triangulation
        tcl = row["triangulation"]
        cell(4, i, CLASS_COLORS[tcl])
        text(4, i, tcl, color="white", weight="bold", size=5.4)

    # Axes labeling
    ax.set_xlim(-0.55, len(col_x) - 0.45)
    ax.set_ylim(n - 0.4, -0.6)
    ax.set_xticks(list(col_x.values()))
    ax.set_xticklabels(list(col_x.keys()), fontsize=6.4, fontweight="bold")
    ax.tick_params(axis="x", length=0, pad=5)
    ax.xaxis.set_ticks_position("top")
    ax.set_yticks(range(n))
    ax.set_yticklabels(df["program_label"].tolist(),
                       fontsize=6.2)
    ax.tick_params(axis="y", length=0, pad=3)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    # subtle row-group separators between metabolic (LA-high) and structural (RA-high)
    flips = [k for k in range(1, n)
             if df["baseline_direction"].iloc[k] != df["baseline_direction"].iloc[k - 1]]
    for k in flips:
        ax.axhline(k - 0.5, color="#BDBDBD", linewidth=0.5, zorder=1)

    fig.suptitle(
        "Cross-figure triangulation of the 10 frozen chamber programs",
        fontsize=7.6, y=0.96, fontweight="bold",
    )

    # Bottom legend strip — triangulation class colors
    leg_y = 0.060
    items = list(CLASS_COLORS.items())
    x0 = 0.27
    spacing = 0.175
    for k, (lab, col) in enumerate(items):
        x = x0 + k * spacing
        fig.patches.append(Rectangle((x, leg_y), 0.016, 0.022,
                                     transform=fig.transFigure,
                                     facecolor=col, edgecolor="white",
                                     linewidth=0.4))
        fig.text(x + 0.020, leg_y + 0.011, lab,
                 fontsize=5.5, ha="left", va="center", color="#222222")
    fig.text(0.27, leg_y - 0.024,
             "LA AF carrier (Fig 3): all rows Unresolved → uniform grey, no per-row label.",
             fontsize=4.8, ha="left", va="top", color="#777777", style="italic")
    fig.text(0.27, leg_y - 0.044,
             "Rules locked in script header (auditable); duration ρ uses |ρ| ≥ 0.30 cutoff.",
             fontsize=4.8, ha="left", va="top", color="#777777", style="italic")

    fig.text(0.012, 0.965, "a",
             fontsize=10, fontweight="bold", ha="left", va="top")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)

    # Also write the assembled triangulation table for the record
    out_tsv = ROOT / "results" / "figure5_triangulation_table.tsv"
    df[["program_id", "program_label", "program_order", "baseline_direction",
        "bulk_label", "bulk_strength", "AF_major_carrier_CT", "rho", "p_perm",
        "dur_dir", "triangulation"]].to_csv(out_tsv, sep="\t", index=False)
    print(f"wrote {OUT.relative_to(ROOT)}  +  {out_tsv.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
