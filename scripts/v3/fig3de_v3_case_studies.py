#!/usr/bin/env python3
"""
Supporting gene-level decompositions for F15 (Relay) and F16 (Extension).

Formerly main-text Fig 3d/3e; main Figure 3 now uses family-level donor
summaries (3C/3D). These panels remain as Supporting Figures.

F15 NTRK2 (cellular relay; endothelial-only family support):
    left  = endothelial (sig)
    right = fibroblast (HCA carrier reference; non-sig at family level)

F16 Keratan sulfate (cellular extension):
    fibroblast (carrier) + atrialCM (neighbour) — both upregulated.

Outputs
-------
figures/figS_gene_decomp_F15_NTRK2_relay.pdf
figures/figS_gene_decomp_F16_Keratan_extension.pdf
(also mirrors legacy fig3_d_* / fig3_e_* paths for backwards compatibility)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_style import enforce_font_floor, FS_PANEL_TITLE, FS_TICK, FS_AXIS_LABEL, FS_FOOTNOTE

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
        "left_ct": "endothelial",
        "right_ct": "fibroblast",
        "carrier_ct": "fibroblast",
        "panel_letter": "S (F15 gene)",
        "out_pdf": OUT_DIR / "figS_gene_decomp_F15_NTRK2_relay.pdf",
        "legacy_pdf": OUT_DIR / "fig3_d_F15_NTRK2_relay.pdf",
        "interpretation_lines": [
            "Family-level support after exact-singscore correction is endothelial-only (Relay); fibroblast is the HCA carrier reference.",
            "Gene bars in endothelial show a predominantly AF-time downshift; fibroblast carrier bars are shown for context.",
            "Lymphocyte is no longer a family-level support compartment and is omitted from this supporting panel.",
            "Gene-level marks are descriptive and do not redefine family classification.",
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
        "panel_letter": "S (F16 gene)",
        "out_pdf": OUT_DIR / "figS_gene_decomp_F16_Keratan_extension.pdf",
        "legacy_pdf": OUT_DIR / "fig3_e_F16_Keratan_extension.pdf",
        "interpretation_lines": [
            "Both compartments show concordant AF-time upshift (positive RED bars), matching family-level Extension.",
            "The baseline carrier (fibroblast) is co-active with neighbouring atrialCM; this is cellular extension,",
            "not relay. Gene-level marks are Supporting only; main Fig 3D uses family-level donor summaries.",
        ],
    },
}


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": enforce_font_floor(FS_TICK),
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
    fig = plt.figure(figsize=(170 * MM, max(95, 44 + 8.5 * n_genes) * MM))  # taller for larger fonts
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1, 1],
        left=0.10, right=0.96, top=0.82, bottom=0.28, wspace=0.18,
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
                    i, "*", fontsize=enforce_font_floor(FS_PANEL_TITLE + 2), fontweight="bold", color="#1A1A1A",
                    ha="left" if v >= 0 else "right", va="center", zorder=5,
                )
            elif p < 0.10:
                ax.text(
                    v + (max_dr * 0.04 if v >= 0 else -max_dr * 0.04),
                    i, ".", fontsize=enforce_font_floor(FS_PANEL_TITLE + 2), fontweight="bold", color="#444444",
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
        ax.set_yticklabels(gene_order, fontsize=enforce_font_floor(FS_TICK - 0.2), fontweight="bold")
        ax.set_xlabel("\u0394 norm-rank  (AF \u2212 Ctrl)", fontsize=enforce_font_floor(FS_AXIS_LABEL))
        ax.tick_params(axis="x", labelsize=enforce_font_floor(FS_TICK - 0.5))
        ax.set_title(
            title,
            fontsize=enforce_font_floor(FS_PANEL_TITLE), fontweight="bold",
            color=COMPARTMENT_PALETTE.get(ct, "#1A1A1A"), loc="left", pad=4,
        )
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.spines["left"].set_linewidth(0.4)
        ax.spines["bottom"].set_linewidth(0.4)

    _draw_panel(ax_left, left_df, left_ct, COMPARTMENT_SHORT[left_ct] + "  AF\u2212Ctrl")
    _draw_panel(ax_right, right_df, right_ct, COMPARTMENT_SHORT[right_ct] + "  AF\u2212Ctrl")

    fig.text(
        0.10, 0.92,
        f"Fig {cfg['panel_letter']}  {cfg['label']}  \u2014  {cfg['class_name']}",
        fontsize=enforce_font_floor(FS_PANEL_TITLE + 2), fontweight="bold", color="#1A1A1A",
    )
    fig.add_artist(
        mpatches.FancyBboxPatch(
            (0.78, 0.905), 0.16, 0.040,
            boxstyle="round,pad=0.01,rounding_size=0.012",
            facecolor=cfg["class_color"], edgecolor="none",
            transform=fig.transFigure,
        )
    )
    fig.text(0.86, 0.925, cfg["class_name"], fontsize=enforce_font_floor(FS_TICK), fontweight="bold",
             color="white", ha="center", va="center")

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
    line_h = 0.020
    for k, line in enumerate(footer_lines):
        fig.text(
            0.05, y_top - line_h * k,
            line, fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5), color="#444444",
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(cfg["out_pdf"], bbox_inches=None, pad_inches=0.04)
    print("wrote " + str(cfg["out_pdf"].relative_to(ROOT)))
    legacy = cfg.get("legacy_pdf")
    if legacy is not None:
        fig.savefig(legacy, bbox_inches=None, pad_inches=0.04)
        print("wrote " + str(legacy.relative_to(ROOT)) + " (legacy mirror)")
    plt.close(fig)


def main():
    set_style()
    delta = pd.read_csv(DELTA_TSV, sep="\t")
    for key, cfg in CASES.items():
        render_case(delta, key, cfg)


if __name__ == "__main__":
    main()
