#!/usr/bin/env python3
"""
Fig 3f v3 mini panel - F4 BCAA catabolism (program-level relay).

Single horizontal bar plot of BCAA program dScore across all 9 LA snRNA
compartments.  Highlights fibroblast as the AF-active site even though
the F4 family carrier (Fig 1F HCA) is adipocyte.  Sub-family signal that
gets diluted at family-mean level.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PROG_TSV = ROOT / "results" / "fig3_v3_program_celltype_delta.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
OUT_PDF = ROOT / "figures" / "fig3_v3_f_F4_BCAA_mini.pdf"

PROGRAM_ID = "REACTOME__REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM"
PROGRAM_LABEL = "F4 BCAA catabolism (program-level)"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
COMPARTMENT_PALETTE = {
    "atrialCM": "#4E79A7",
    "fibroblast": "#A0CBE8",
    "endothelial": "#F28E2B",
    "macrophage": "#59A14F",
    "lymphocyte": "#8CD17D",
    "pericyte": "#499894",
    "vascular_smooth_muscle": "#86BCB6",
    "adipocyte": "#E15759",
    "neuronal": "#79706E",
}
COMPARTMENT_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC",
    "adipocyte": "Adipo",
    "neuronal": "Neuro",
}


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def main():
    set_style()
    df = pd.read_csv(PROG_TSV, sep="\t")
    df = df[df["program_id"] == PROGRAM_ID].copy()
    df = df.sort_values("dScore", ascending=False).reset_index(drop=True)
    cts = df["compartment"].astype(str).tolist()

    hca = pd.read_csv(HCA_TSV, sep="\t")
    family_carrier = hca.loc[hca["family_id"] == 4, "LA_carrier"]
    family_carrier = str(family_carrier.iloc[0]) if not family_carrier.empty else "unknown"

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(170 * MM, 110 * MM))
    gs = fig.add_gridspec(
        1, 1, left=0.18, right=0.96, top=0.83, bottom=0.36,
    )
    ax = fig.add_subplot(gs[0, 0])

    n = len(df)
    y = range(n)
    for i, r in df.iterrows():
        v = float(r["dScore"])
        p = float(r["p_perm"])
        ct = str(r["compartment"])
        color = LA_DEEP if v >= 0 else RA_DEEP
        ax.barh(i, v, height=0.65, color=color, edgecolor="#1A1A1A", linewidth=0.4, alpha=0.92)
        # carrier marker (left margin)
        if ct == family_carrier:
            ax.text(-0.024, i, "\u25C0  family carrier", fontsize=5.6, color="#888888",
                    style="italic", ha="left", va="center")
        # significance text
        if p < 0.05:
            ax.text(v + (0.001 if v >= 0 else -0.001), i, f"  *  p={p:.3f}",
                    fontsize=6.4, fontweight="bold", color="#1A1A1A",
                    ha="left" if v >= 0 else "right", va="center")
        elif p < 0.10:
            ax.text(v + (0.001 if v >= 0 else -0.001), i, f"  .  p={p:.3f}",
                    fontsize=5.8, color="#444444",
                    ha="left" if v >= 0 else "right", va="center")

    ax.axvline(0, color="#1A1A1A", lw=0.5)
    ax.set_yticks(list(y))
    ax.set_yticklabels([COMPARTMENT_SHORT[c] for c in cts], fontsize=7.0, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("\u0394singscore  (AF \u2212 Ctrl)", fontsize=7.4)
    ax.tick_params(axis="x", labelsize=6.4)
    ax.set_xlim(-0.025, 0.025)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(0.4)
    ax.spines["bottom"].set_linewidth(0.4)

    fig.suptitle(
        f"Fig 3f  {PROGRAM_LABEL}",
        x=0.05, y=0.94, ha="left", fontsize=8.4, fontweight="bold",
    )
    fig.text(
        0.18, 0.86,
        "Sub-family relay  -  family carrier (Fig 1F HCA) = "
        f"{COMPARTMENT_SHORT.get(family_carrier, family_carrier)}",
        fontsize=6.4, color="#444444", style="italic",
    )

    footer_lines = [
        "Although F4 (4 metabolic programs) shows no family-level cellular signal when averaged,",
        "its BCAA catabolism component reveals fibroblast as the AF-active site (p = 0.011) \u2014",
        "a sub-family relay from the family's adipocyte baseline carrier.",
        "Consistent with the documented metabolic reprogramming of cardiac fibroblasts during AF-associated remodelling.",
        "",
        "Note on Lymph: lymphocyte shows the largest \u0394Score (+0.023) but does not reach nominal p < 0.05",
        "(p = 0.14) because of higher within-group variance in the smaller-N lymphocyte sample (n = 10 Ctrl).",
        "Fibroblast is therefore the only compartment passing the gene-level significance threshold here.",
        "",
        "Color note: bars are colored by direction (RED = AF > Ctrl upshift, BLUE = AF < Ctrl downshift),",
        "the same convention as Fig 3b, 3d and 3e.  Color encodes direction (up/down), NOT cellular-fidelity",
        "class (relay/extension/fidelity); the latter is shown by the badge color in each case-study panel.",
    ]
    y_top = 0.20
    line_h = 0.018
    for k, line in enumerate(footer_lines):
        fig.text(
            0.05, y_top - line_h * k,
            line, fontsize=5.6, color="#444444",
        )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
