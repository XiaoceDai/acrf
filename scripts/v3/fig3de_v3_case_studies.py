#!/usr/bin/env python3
"""
Fig 3d / 3e v3 case-study panels.

Fig 3d: F15 NTRK2 (cellular relay)
    Two side-by-side per-gene Delta-rank bar plots:
        left  = lymphocyte (sig)
        right = endothelial (sig)
    Each gene bar is overlaid with a grey vertical reference line marking
    the gene's healthy-LA log1p-CPM in the BASELINE CARRIER (fibroblast,
    HCA paired-LA donors).  This visualises the relay signature: black bars
    drift one direction (AF-driven cellular shift in non-carrier ct) while
    grey lines sit at a different position (the resting executor remains
    the carrier).

Fig 3e: F16 Keratan sulfate (cellular extension)
    Same dual-bar layout, but the sig cell types are fibroblast (the carrier
    itself) and atrialCM (a neighbour) - both upregulated.

Inputs
------
results/fig3_v3_per_gene_delta.tsv      (already computed by step 1)
results/fig1f_v3_HCA_family_carriers.tsv

Outputs
-------
figures/fig3_v3_d_F15_NTRK2_relay.pdf
figures/fig3_v3_e_F16_Keratan_extension.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

DELTA_TSV = ROOT / "results" / "fig3_v3_per_gene_delta.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"

OUT_DIR = ROOT / "figures"

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
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

CASES = {
    "fig3d": {
        "program": "REACTOME__REACTOME_SIGNALING_BY_NTRK2_TRKB",
        "label": "F15 NTRK2 / TrkB",
        "class_name": "Cellular relay",
        "class_color": "#8E6BAA",
        "left_ct": "lymphocyte",
        "right_ct": "endothelial",
        "carrier_ct": "fibroblast",
        "panel_letter": "3d",
        "out_pdf": OUT_DIR / "fig3_v3_d_F15_NTRK2_relay.pdf",
        "interpretation_lines": [
            "Bars in non-carrier compartments (lymphocyte, endothelial) reveal a predominantly AF-time downshift.",
            "The bulk-detected attenuation (Fig 2C, F15 = Att 2/3) therefore originates from cellular relay rather than carrier-intrinsic change.",
            "Note that the relay is NOT a uniform downshift of every LE gene: SRC (in lymphocyte) and PLCG1 (in endothelial) show modest upward",
            "drift, indicating compartment-specific divergence within the same pathway.",
        ],
    },
    "fig3e": {
        "program": "REACTOME__REACTOME_KERATAN_SULFATE_BIOSYNTHESIS",
        "label": "F16 Keratan sulfate biosynthesis",
        "class_name": "Cellular extension",
        "class_color": "#5F9EA0",
        "left_ct": "fibroblast",
        "right_ct": "atrialCM",
        "carrier_ct": "fibroblast",
        "panel_letter": "3e",
        "out_pdf": OUT_DIR / "fig3_v3_e_F16_Keratan_extension.pdf",
        "interpretation_lines": [
            "Both compartments show concordant AF-time upshift (positive RED bars), in contrast to the predominantly",
            "downshift (blue) seen in the relay case (Fig 3d).  The baseline carrier (fibroblast) is co-active with the",
            "neighbouring cardiomyocyte compartment; the bulk-detected change therefore reflects cellular extension",
            "\u2014 the carrier's activity change is amplified by adjacent cell types rather than re-routed.",
            "Note that OGN, OMD and CHST3 show small divergent (negative) drift in atrialCM despite the overall upshift,",
            "indicating that extension is not a perfectly uniform pathway-wide event; the dominant extension signal is",
            "carried by ACAN, FMOD and PRELP (gene-level p < 0.05 in both compartments).",
        ],
    },
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


def render_case(delta, key, cfg):
    sub = delta[delta["program_id"] == cfg["program"]]
    if sub.empty:
        print(f"[SKIP] no data for {cfg['program']}")
        return
    left_ct = cfg["left_ct"]; right_ct = cfg["right_ct"]; carrier = cfg["carrier_ct"]

    left_df = sub[sub["compartment"] == left_ct].copy()
    right_df = sub[sub["compartment"] == right_ct].copy()
    carrier_df = sub[sub["compartment"] == carrier].copy()

    # Gene order: by absolute LA-baseline (carrier) expression descending,
    # so genes with the strongest baseline carrier signal appear at top.
    gene_order = (
        carrier_df.assign(abs_ref=lambda d: d["baseline_LA_logcpm"].abs())
        .sort_values("abs_ref", ascending=False)["gene"].tolist()
    )
    if not gene_order:
        gene_order = sorted(sub["gene"].unique())

    # Filter to genes present in both panels
    left_df = left_df.set_index("gene").reindex(gene_order)
    right_df = right_df.set_index("gene").reindex(gene_order)
    carrier_df = carrier_df.set_index("gene").reindex(gene_order)

    # Reference line position: scale carrier baseline to share x-range
    # with dRank.  We map [0, max_baseline] -> [0, abs(max(dRank))] visually.
    all_dr = pd.concat([left_df["dRank"], right_df["dRank"]]).abs()
    max_dr = max(0.05, float(all_dr.max() * 1.15))

    n_genes = len(gene_order)
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(170 * MM, max(85, 38 + 8.0 * n_genes) * MM))
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1, 1],
        left=0.10, right=0.96, top=0.84, bottom=0.30, wspace=0.18,
    )
    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1])

    def _draw_panel(ax, df, ct, title):
        y = np.arange(n_genes)
        d_vals = df["dRank"].to_numpy(dtype=float)
        p_vals = df["p_perm"].to_numpy(dtype=float)
        baseline = df["baseline_LA_logcpm"].to_numpy(dtype=float)
        # Map baseline log-CPM (positive scalar) onto a horizontal grey marker
        # rescaled into +x_axis units.  We use baseline / global_max_baseline * max_dr
        # so the line position is comparable across genes.
        if np.isfinite(baseline).any():
            ref_max = float(np.nanmax(baseline)) or 1.0
            ref_scaled = baseline / ref_max * max_dr * 0.85  # leave headroom
        else:
            ref_scaled = np.zeros_like(d_vals)

        for i in range(n_genes):
            v = d_vals[i] if np.isfinite(d_vals[i]) else 0.0
            color = LA_DEEP if v >= 0 else RA_DEEP
            ax.barh(
                i, v, height=0.65,
                color=color, edgecolor="#1A1A1A", linewidth=0.4, alpha=0.9, zorder=3,
            )
            # significance asterisk
            p = p_vals[i] if np.isfinite(p_vals[i]) else 1.0
            if p < 0.05:
                ax.text(
                    v + (max_dr * 0.04 if v >= 0 else -max_dr * 0.04),
                    i, "*", fontsize=10, fontweight="bold", color="#1A1A1A",
                    ha="left" if v >= 0 else "right", va="center", zorder=5,
                )
            elif p < 0.10:
                ax.text(
                    v + (max_dr * 0.04 if v >= 0 else -max_dr * 0.04),
                    i, ".", fontsize=10, fontweight="bold", color="#444444",
                    ha="left" if v >= 0 else "right", va="center", zorder=5,
                )
            # Carrier reference line (vertical short tick at scaled position)
            if np.isfinite(ref_scaled[i]):
                ax.plot(
                    [ref_scaled[i], ref_scaled[i]], [i - 0.32, i + 0.32],
                    color="#888888", linewidth=1.2, zorder=4,
                )
        ax.axvline(0, color="#1A1A1A", lw=0.5, zorder=2)
        ax.set_xlim(-max_dr, max_dr)
        ax.set_ylim(-0.5, n_genes - 0.5)
        ax.invert_yaxis()
        ax.set_yticks(np.arange(n_genes))
        ax.set_yticklabels(gene_order, fontsize=6.6, fontweight="bold")
        ax.set_xlabel("\u0394 norm-rank  (AF \u2212 Ctrl)", fontsize=7.0)
        ax.tick_params(axis="x", labelsize=6.2)
        ax.set_title(
            title,
            fontsize=7.6, fontweight="bold",
            color=COMPARTMENT_PALETTE.get(ct, "#1A1A1A"), loc="left", pad=4,
        )
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.spines["left"].set_linewidth(0.4)
        ax.spines["bottom"].set_linewidth(0.4)

    _draw_panel(ax_left, left_df, left_ct, COMPARTMENT_SHORT[left_ct] + "  AF\u2212Ctrl")
    _draw_panel(ax_right, right_df, right_ct, COMPARTMENT_SHORT[right_ct] + "  AF\u2212Ctrl")

    # ---- header strip with class badge ----
    fig.text(
        0.10, 0.92,
        f"Fig {cfg['panel_letter']}  {cfg['label']}  \u2014  {cfg['class_name']}",
        fontsize=9.2, fontweight="bold", color="#1A1A1A",
    )
    fig.add_artist(
        mpatches.FancyBboxPatch(
            (0.78, 0.905), 0.16, 0.040,
            boxstyle="round,pad=0.01,rounding_size=0.012",
            facecolor=cfg["class_color"], edgecolor="none",
            transform=fig.transFigure,
        )
    )
    fig.text(0.86, 0.925, cfg["class_name"], fontsize=6.6, fontweight="bold",
             color="white", ha="center", va="center")

    # ---- footer with full caption (each line manually wrapped to fit 170mm width) ----
    footer_lines = [
        "Bars: per-gene \u0394 norm-rank (AF \u2212 Ctrl) on GSE255612 LA snRNA pseudobulk (n = 18 AF + 16 Ctrl donors).",
        "Red = AF > Ctrl (upshift), Blue = AF < Ctrl (downshift) \u2014 same convention as Fig 3b cells.",
        f"Grey vertical ticks: per-gene healthy log1p-CPM in the baseline carrier ({COMPARTMENT_SHORT[carrier]}, HCA paired-LA donors),",
        "rescaled to panel x-range.",
        "",
    ]
    footer_lines.extend(cfg["interpretation_lines"])
    footer_lines.extend([
        "",
        "Significance markers: * = gene-level p < 0.05;  . = gene-level p < 0.10  (donor permutation, n_perm = 5000).",
        "Color note: cellular-fidelity badge colors (purple = relay, teal = extension) are decoupled",
        "from the carrier palette in Fig 3b to avoid cross-figure confusion.",
    ])
    y_top = 0.27
    line_h = 0.018
    for k, line in enumerate(footer_lines):
        fig.text(
            0.05, y_top - line_h * k,
            line, fontsize=5.6, color="#444444",
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(cfg["out_pdf"], bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(cfg["out_pdf"].relative_to(ROOT)))


def main():
    set_style()
    delta = pd.read_csv(DELTA_TSV, sep="\t")
    for key, cfg in CASES.items():
        render_case(delta, key, cfg)


if __name__ == "__main__":
    main()
