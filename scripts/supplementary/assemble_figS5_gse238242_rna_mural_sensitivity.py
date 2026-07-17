#!/usr/bin/env python3
"""
Assemble final Supplementary Figure S5 — GSE238242 RNA-side and mural-cell sensitivity.

A = RNA-side four-compartment sensitivity (frozen Project_P5 main delta)
B = combined mural-cell sensitivity (frozen Project_P5 mural delta)

No score recomputation. Publication-clean titles; no regulator-tier language.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "supplementary"))

from plot_style import AFCTRL_CMAP  # noqa: E402
from _export import save_fig_bundle, set_pub_style  # noqa: E402

IN = ROOT / "Project_P5" / "GSE238242" / "results"
DELTA_MAIN = IN / "figure5_program_delta_main.tsv"
SUPPORT = IN / "figure5_program_support_check.tsv"
DELTA_MURAL = IN / "figure5_program_delta_mural.tsv"

PROGRAM_ORDER = [
    ("ECM proteoglycans", "STRUCTURAL"),
    ("HS-GAG Metabolism", "STRUCTURAL"),
    ("HS-GAG Biosynthesis", "STRUCTURAL"),
    ("OXPHOS", "METABOLIC"),
    ("Aerobic Resp. & ETC", "METABOLIC"),
    ("Cytoplasmic Translation", "METABOLIC"),
]
COMPARTMENT_ORDER = ["atrialCM", "fibroblast", "endothelial_or_vascular", "immune"]
COMPARTMENT_SHORT = {
    "atrialCM": "atrialCM",
    "fibroblast": "fibroblast",
    "endothelial_or_vascular": "endo/vasc",
    "immune": "immune",
}
EXPECT = {
    "STRUCTURAL": ("fibroblast", +1),
    "METABOLIC": ("atrialCM", -1),
}
LEVEL_COL = {
    "none": "#888888",
    "weak": "#B8860B",
    "moderate": "#D2691E",
    "strong": "#8B0000",
}
LEVEL_DOT = {"none": "○", "weak": "◐", "moderate": "●", "strong": "★"}


def main() -> dict:
    set_pub_style()
    delta = pd.read_csv(DELTA_MAIN, sep="\t")
    sup = pd.read_csv(SUPPORT, sep="\t")
    mural = pd.read_csv(DELTA_MURAL, sep="\t")

    n_prog = len(PROGRAM_ORDER)
    n_cmp = len(COMPARTMENT_ORDER)
    M = np.full((n_prog, n_cmp), np.nan)
    for i, (lab, _ax) in enumerate(PROGRAM_ORDER):
        for j, c in enumerate(COMPARTMENT_ORDER):
            row = delta[(delta["program_label"] == lab) & (delta["compartment"] == c)]
            if len(row) == 1:
                M[i, j] = float(row["delta_score"].iloc[0])

    # mural values aligned to program order
    mural_vals = []
    for lab, _ax in PROGRAM_ORDER:
        row = mural[mural["program_label"] == lab]
        if len(row) == 0 and "program" in mural.columns:
            row = mural[mural["program"] == lab]
        # try flexible
        if len(row) == 0:
            for col in mural.columns:
                if "label" in col.lower() or col.lower() == "program":
                    row = mural[mural[col] == lab]
                    if len(row):
                        break
        val_col = "delta_score" if "delta_score" in mural.columns else (
            "delta" if "delta" in mural.columns else mural.columns[-1]
        )
        mural_vals.append(float(row[val_col].iloc[0]) if len(row) else np.nan)
    mural_vals = np.asarray(mural_vals, dtype=float)

    MM = 1 / 25.4
    fig = plt.figure(figsize=(180 * MM, 120 * MM))
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.35, 0.85],
        left=0.14, right=0.96, top=0.82, bottom=0.18, wspace=0.35,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    vmax = max(0.05, float(np.nanmax(np.abs(M))) * 1.05)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list("af_ctrl", AFCTRL_CMAP, N=256)
    im = ax_a.imshow(M, cmap=cmap, norm=norm, aspect="auto")
    ax_a.set_xticks(range(n_cmp))
    ax_a.set_xticklabels([COMPARTMENT_SHORT[c] for c in COMPARTMENT_ORDER], rotation=30, ha="right", fontsize=6.2)
    ax_a.set_yticks(range(n_prog))
    ax_a.set_yticklabels(
        [f"{lab}\n({axg})" for lab, axg in PROGRAM_ORDER],
        fontsize=5.8,
    )
    # mark expected cells
    for i, (lab, axg) in enumerate(PROGRAM_ORDER):
        exp_ct, _sign = EXPECT[axg]
        if exp_ct in COMPARTMENT_ORDER:
            j = COMPARTMENT_ORDER.index(exp_ct)
            ax_a.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#111111", lw=1.2))
    ax_a.set_title("A  RNA-side four-compartment sensitivity", loc="left", fontsize=7.4, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax_a, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=5.5)
    cbar.set_label("Δsingscore (AF − SR)", fontsize=6.0)

    # Panel B mural bars
    vmax_b = max(0.05, float(np.nanmax(np.abs(mural_vals))) * 1.1)
    colors = ["#8C510A" if v >= 0 else "#01665E" for v in mural_vals]
    y = np.arange(n_prog)
    ax_b.barh(y, mural_vals, color=colors, edgecolor="#222222", linewidth=0.4, height=0.7)
    ax_b.axvline(0, color="#222222", lw=0.5)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([lab for lab, _ in PROGRAM_ORDER], fontsize=5.8)
    ax_b.invert_yaxis()
    ax_b.set_xlim(-vmax_b, vmax_b)
    ax_b.set_xlabel("Δsingscore (AF − SR)", fontsize=6.2)
    ax_b.set_title("B  Combined mural-cell sensitivity", loc="left", fontsize=7.4, fontweight="bold")
    ax_b.text(
        0.98, -0.12, "mural = pericyte + VSMC",
        transform=ax_b.transAxes, ha="right", va="top", fontsize=5.4, color="#555555",
    )
    for s in ("top", "right"):
        ax_b.spines[s].set_visible(False)

    fig.suptitle(
        "GSE238242 RNA-side and mural-cell sensitivity",
        x=0.14, y=0.96, ha="left", fontsize=9.0, fontweight="bold",
    )
    fig.text(
        0.14, 0.90,
        "Supplementary sensitivity only · not a Path A gate · not an independent main-text evidence layer · "
        "does not alter Figure 3 localization classes or Figure 5 support counts · "
        "thick boxes mark locked expected RNA-side cells",
        fontsize=5.5, color="#333333", ha="left",
    )
    fig.text(
        0.50, 0.04,
        "Frozen GSE238242 LAA snRNA pseudobulk Δsingscore (AF − SR). "
        "No expected RNA-side cell reaches moderate/strong support under the locked support-gate rule.",
        ha="center", va="bottom", fontsize=5.3, color="#444444",
    )

    paths = save_fig_bundle(fig, "figS5_gse238242_rna_mural_sensitivity", pad_inches=0.03)
    # verify no regulator-tier in support table usage — text scan later
    _ = sup  # lineage retained
    return {
        "paths": paths,
        "panel_structure": "A=RNA-side 4-compartment; B=mural combined",
        "notes": "Merged old S5+S6 from frozen Project_P5 sources; regulator-tier removed",
    }


if __name__ == "__main__":
    print(main()["paths"]["pdf"])
