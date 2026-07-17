#!/usr/bin/env python3
"""Corrected Fig 4B — 16-family LAA/LAPV Spearman trend heatmap (canonical order)."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
OUT_BASE = FIGURES / "fig4b_family_duration_heatmap"
WIDTH_MM, HEIGHT_MM = 178.0, 90.0

CHRONICITY_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"
CANONICAL_TSV = ROOT / "results" / "ncp_v3_canonical_display_order.tsv"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"

VIOLET = "#6B4C9A"
OCHRE = "#C8A076"
VMIN, VMAX = -1.0, 1.0

COMPARTMENT_SHORT = {
    "atrialCM": "aCM", "fibroblast": "Fibro", "endothelial": "Endo",
    "macrophage": "Macro", "lymphocyte": "Lymph", "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC", "adipocyte": "Adipo", "neuronal": "Neuro",
    "Mixed/Unresolved": "Mix",
}


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 6.5,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def canonical_family_order() -> list[int]:
    order = (
        pd.read_csv(CANONICAL_TSV, sep="\t")
        .sort_values("family_display_order")
        .drop_duplicates("family_id")
    )
    return order["family_id"].astype(int).tolist()


def regional_direction(rho_laa: float, rho_lapv: float) -> str:
    if not np.isfinite(rho_lapv):
        return "NA"
    if rho_laa == 0 or rho_lapv == 0:
        return "same" if np.sign(rho_laa) == np.sign(rho_lapv) else "opposite"
    return "same" if np.sign(rho_laa) == np.sign(rho_lapv) else "opposite"


def save_outputs(fig: plt.Figure) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BASE.with_suffix(".pdf"), bbox_inches=None, pad_inches=0.02)
    fig.savefig(OUT_BASE.with_suffix(".svg"), bbox_inches=None, pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    set_style()
    cls = pd.read_csv(CHRONICITY_TSV, sep="\t").set_index("family_id")
    labels = pd.read_csv(LABELS_TSV, sep="\t").set_index("family_id")
    hca = pd.read_csv(HCA_TSV, sep="\t").set_index("family_id")
    fam_order = canonical_family_order()
    assert len(fam_order) == 16

    rho_laa = np.array([float(cls.loc[fid, "rho_primary"]) for fid in fam_order])
    rho_lapv = np.array([float(cls.loc[fid, "rho_lapv"]) for fid in fam_order])
    directions = [regional_direction(rho_laa[i], rho_lapv[i]) for i in range(16)]

    heat = np.column_stack([rho_laa, rho_lapv])
    cmap = LinearSegmentedColormap.from_list("rho_vwo", [(0.0, VIOLET), (0.5, "#FFFFFF"), (1.0, OCHRE)])
    cmap.set_bad("#E0E0E0")
    norm = TwoSlopeNorm(vmin=VMIN, vcenter=0.0, vmax=VMAX)

    mm = 1.0 / 25.4
    fig = plt.figure(figsize=(WIDTH_MM * mm, HEIGHT_MM * mm))
    gs = fig.add_gridspec(1, 4, width_ratios=[0.8, 2.8, 2.2, 0.7],
                          left=0.06, right=0.94, top=0.90, bottom=0.12, wspace=0.08)
    ax_carrier = fig.add_subplot(gs[0, 0])
    ax_labels = fig.add_subplot(gs[0, 1])
    ax_heat = fig.add_subplot(gs[0, 2])
    ax_dir = fig.add_subplot(gs[0, 3])

    for ax in (ax_carrier, ax_labels, ax_heat, ax_dir):
        ax.set_xlim(0, 1)
        ax.set_ylim(16, 0)
        ax.axis("off")

    for i, fid in enumerate(fam_order):
        carrier = str(hca.loc[fid, "LA_carrier"]) if fid in hca.index else "Mix"
        short = COMPARTMENT_SHORT.get(carrier, carrier[:4])
        ax_carrier.text(0.5, i + 0.5, short, ha="center", va="center", fontsize=5.8, color="#444444")
        direction = str(labels.loc[fid, "axis_group"])
        short_label = str(labels.loc[fid, "short_label"])
        ax_labels.text(0.02, i + 0.5, f"F{fid:2d}", ha="left", va="center", fontsize=6.0, fontweight="bold")
        ax_labels.text(0.12, i + 0.5, f"{direction[0]}  {short_label[:28]}",
                       ha="left", va="center", fontsize=5.6, color="#333333")

    im = ax_heat.imshow(heat, aspect="auto", cmap=cmap, norm=norm)
    ax_heat.set_xticks([0, 1])
    ax_heat.set_xticklabels(["LAA ρ", "LAPV ρ"], fontsize=6.5)
    ax_heat.set_yticks(range(16))
    ax_heat.set_yticklabels([])
    for i in range(16):
        for j in range(2):
            val = heat[i, j]
            if not np.isfinite(val):
                continue
            txt_color = "white" if abs(val) >= 0.35 else "#1A1A1A"
            ax_heat.text(j, i, f"{val:+.2f}", ha="center", va="center",
                         fontsize=5.8, color=txt_color, fontweight="bold")

    for i, direction in enumerate(directions):
        color = "#2E7B5E" if direction == "same" else ("#9C5550" if direction == "opposite" else "#888888")
        ax_dir.text(0.5, i + 0.5, direction, ha="center", va="center", fontsize=5.8, color=color, fontweight="bold")

    cbar_ax = fig.add_axes([0.55, 0.04, 0.30, 0.025])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cbar.ax.tick_params(labelsize=5.5)
    cbar.set_label("Spearman ρ (fixed scale −1 to +1)", fontsize=6.0)

    fig.text(0.06, 0.95, "Fig 4b  Family-level duration trends (canonical Fig 1D order, n = 12 AF-LAA)",
             fontsize=7.5, fontweight="bold", color="#1A1A1A", ha="left")
    fig.text(0.06, 0.91,
             "Program-first family activity scores (exact-reference).  ρ text only; no p-value stars.  "
             "Regional direction: LAA vs LAPV sign concordance (sensitivity only).",
             fontsize=5.5, color="#444444", ha="left")

    save_outputs(fig)
    print(f"wrote {OUT_BASE.with_suffix('.pdf').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
