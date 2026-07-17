#!/usr/bin/env python3
"""
Fig 4b v3 - 16 family x 4 sensitivity coding Spearman rho heatmap (LA chronicity).

Layout:
    [carrier] [state] [family label] [4-coding rho heatmap with rho text]
    [sensitivity tally column] [legend column]

4 sensitivity codings (left -> right):
    Primary    : numeric duration_months (>1M dropped)
    gt1m_eq1   : ">1M" coded as 1 month
    Ordinal    : 3-tier (SR=0, short_AF=1, long_AF=2)
    LAPV       : LAPV region sensitivity (numeric duration)

Row order: atlas PriorityScore descending (cross-figure consistency with Fig 3b).
Sig markers: yellow * for p<0.05, "." for p<0.10  (donor-permutation, n_perm = 10,000).
Each cell shows Spearman rho text (2 decimals) for direct readability.
Right tally: number of codings with same sign as primary.
"""
from __future__ import annotations

import sys
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_style import enforce_font_floor, FS_PANEL_TITLE, FS_TICK, FS_AXIS_LABEL, FS_LEGEND, FS_FOOTNOTE

CLS_TSV = ROOT / "results" / "fig4_v3_chronicity_class.tsv"
ATLAS_TSV = ROOT / "results" / "fig_integrated_af_rewiring_atlas.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
PRIMARY_TSV = ROOT / "results" / "fig4_v3_FAMILY_laa_primary.tsv"
GT1M_TSV = ROOT / "results" / "fig4_v3_FAMILY_laa_sens_gt1m_eq1.tsv"
ORD_TSV = ROOT / "results" / "fig4_v3_FAMILY_laa_sens_ordinal.tsv"
LAPV_TSV = ROOT / "results" / "fig4_v3_FAMILY_lapv_sensitivity.tsv"
OUT_PDF = ROOT / "figures" / "fig4_b_family_duration_heatmap.pdf"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
SIG_COLOR = "#FFE066"

CODINGS = ["primary", "gt1m_eq1", "ordinal", "lapv"]
CODING_LABEL = {
    "primary": "Primary\n(numeric mo)",
    "gt1m_eq1": "gt1m_eq1\n(>1M = 1)",
    "ordinal": "Ordinal\n(3-tier)",
    "lapv": "LAPV\n(region sens)",
}
COMPARTMENT_PALETTE = {
    "atrialCM": "#4E79A7", "fibroblast": "#A0CBE8", "endothelial": "#F28E2B",
    "macrophage": "#59A14F", "lymphocyte": "#8CD17D", "pericyte": "#499894",
    "vascular_smooth_muscle": "#86BCB6", "adipocyte": "#E15759", "neuronal": "#79706E",
}
COMPARTMENT_SHORT = {
    "atrialCM": "aCM", "fibroblast": "Fibro", "endothelial": "Endo",
    "macrophage": "Macro", "lymphocyte": "Lymph", "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC", "adipocyte": "Adipo", "neuronal": "Neuro",
}
MIXED_COLOR = "#BDBDBD"
MIXED_LABEL = "Mixed/Unresolved"
STATE_COLORS = {
    "Preserved": "#A9C9B5", "Attenuated": "#E9C27B",
    "Reversed": "#C6A2D7", "Amplified": "#E59595",
}
STATE_SHORT = {"Preserved": "Pre", "Attenuated": "Att", "Reversed": "Rev", "Amplified": "Amp"}

VMIN, VMAX = -0.6, 0.6


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": enforce_font_floor(FS_TICK),
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def main():
    set_style()
    cls = pd.read_csv(CLS_TSV, sep="\t")
    atlas = pd.read_csv(ATLAS_TSV, sep="\t").sort_values("PriorityScore", ascending=False)
    hca = pd.read_csv(HCA_TSV, sep="\t")

    fam_order = atlas["family_id"].astype(int).tolist()
    short_lookup = dict(zip(atlas["family_id"].astype(int), atlas["short"].astype(str)))
    dir_lookup = dict(zip(atlas["family_id"].astype(int), atlas["HealthyBias"].astype(str)))
    state_lookup = dict(zip(atlas["family_id"].astype(int), atlas["state_consensus"].astype(str)))
    carriers = dict(zip(hca["family_id"].astype(int), hca["LA_carrier"].astype(str)))

    # Build rho + p matrices: 4 codings
    coding_dfs = {
        "primary": pd.read_csv(PRIMARY_TSV, sep="\t").set_index("family_id"),
        "gt1m_eq1": pd.read_csv(GT1M_TSV, sep="\t").set_index("family_id"),
        "ordinal": pd.read_csv(ORD_TSV, sep="\t").set_index("family_id"),
        "lapv": pd.read_csv(LAPV_TSV, sep="\t").set_index("family_id"),
    }
    n_rows = len(fam_order)
    n_cols = len(CODINGS)
    rho_mat = np.full((n_rows, n_cols), np.nan)
    p_mat = np.full((n_rows, n_cols), np.nan)
    for j, c in enumerate(CODINGS):
        df = coding_dfs[c]
        for i, fid in enumerate(fam_order):
            if fid in df.index:
                rho_mat[i, j] = float(df.loc[fid, "rho_spearman"])
                p_mat[i, j] = float(df.loc[fid, "p_perm"])

    # Sensitivity tally per row
    cls_idx = cls.set_index("family_id")
    sens_tally = np.array([int(cls_idx.loc[fid, "sens_same_sign_n"]) if fid in cls_idx.index else 0 for fid in fam_order])

    # Color map blue (negative rho) to red (positive rho)
    cmap = LinearSegmentedColormap.from_list(
        "rho_div", [(0.0, RA_DEEP), (0.5, "#FFFFFF"), (1.0, LA_DEEP)]
    )
    cmap.set_bad("#E0E0E0")
    norm = TwoSlopeNorm(vmin=VMIN, vcenter=0.0, vmax=VMAX)
    tally_cmap = LinearSegmentedColormap.from_list(
        "tally", [(0.0, "#FFFFFF"), (1.0, "#9C2A2A")]
    )

    # Layout simplified to 2 left strips (Carrier 5mm + State 3mm) + family label.
    # LA/RA bias is encoded as L/R prefix in the family label (no separate strip).
    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 220 * MM))  # taller for larger fonts
    gs = fig.add_gridspec(
        1, 6,
        width_ratios=[0.55, 0.30, 3.20, 4.40, 0.55, 1.30],
        left=0.04, right=0.985, top=0.91, bottom=0.16,
        wspace=0.035,
    )
    ax_carrier = fig.add_subplot(gs[0, 0])
    ax_state = fig.add_subplot(gs[0, 1])
    ax_lbl = fig.add_subplot(gs[0, 2])
    ax_main = fig.add_subplot(gs[0, 3])
    ax_tally = fig.add_subplot(gs[0, 4])
    ax_leg = fig.add_subplot(gs[0, 5])
    # ax_dir kept as None to avoid breaking later code paths
    ax_dir = None

    # ---- carrier strip
    for i, fid in enumerate(fam_order):
        c = carriers.get(fid, MIXED_LABEL)
        if c == MIXED_LABEL:
            ax_carrier.add_patch(mpatches.Rectangle(
                (0, i), 1, 1, facecolor=MIXED_COLOR, edgecolor="white", linewidth=0.4, hatch="///",
            ))
        else:
            ax_carrier.add_patch(mpatches.Rectangle(
                (0, i), 1, 1, facecolor=COMPARTMENT_PALETTE.get(c, MIXED_COLOR),
                edgecolor="white", linewidth=0.4,
            ))
    ax_carrier.set_xlim(0, 1); ax_carrier.set_ylim(n_rows, 0)
    ax_carrier.set_xticks([0.5]); ax_carrier.set_xticklabels(["LA\nbaseline\ncarrier"], fontsize=enforce_font_floor(FS_TICK - 0.5), color="#444444")
    ax_carrier.set_yticks([]); ax_carrier.tick_params(axis="x", length=0, pad=2)
    for s in ax_carrier.spines.values():
        s.set_visible(False)

    # ---- state strip (Atlas Fig 2C state)
    for i, fid in enumerate(fam_order):
        st = state_lookup.get(fid, "")
        c = STATE_COLORS.get(st, "#E8E8E8")
        ax_state.add_patch(mpatches.Rectangle((0, i), 1, 1, facecolor=c, edgecolor="white", linewidth=0.4))
        ax_state.text(0.5, i + 0.5, STATE_SHORT.get(st, "?"), ha="center", va="center",
                      fontsize=enforce_font_floor(FS_TICK - 0.5), color="#1A1A1A", fontweight="bold")
    ax_state.set_xlim(0, 1); ax_state.set_ylim(n_rows, 0)
    ax_state.set_xticks([0.5]); ax_state.set_xticklabels(["AF\nstate\n(Fig 2C)"], fontsize=enforce_font_floor(FS_TICK - 0.5), color="#444444")
    ax_state.set_yticks([]); ax_state.tick_params(axis="x", length=0, pad=2)
    for s in ax_state.spines.values():
        s.set_visible(False)

    # ---- family labels (LA/RA bias encoded as colored "L" / "R" prefix)
    ax_lbl.set_xlim(0, 1); ax_lbl.set_ylim(n_rows, 0); ax_lbl.axis("off")
    for i, fid in enumerate(fam_order):
        d = dir_lookup.get(fid, "?")
        col = LA_DEEP if d == "LA-biased" else RA_DEEP
        prefix = "[L]" if d == "LA-biased" else "[R]"
        ax_lbl.text(0.99, i + 0.5,
                    f"{prefix}  F{int(fid)}  {short_lookup.get(fid, '')}",
                    fontsize=enforce_font_floor(FS_TICK), color=col, fontweight="bold",
                    ha="right", va="center")

    # ---- main heatmap
    masked = np.ma.masked_invalid(rho_mat)
    im = ax_main.imshow(masked, cmap=cmap, norm=norm, aspect="auto",
                          extent=[-0.5, n_cols - 0.5, n_rows - 0.5, -0.5],
                          interpolation="nearest")
    primary_idx = CODINGS.index("primary")
    for i in range(n_rows):
        v_primary = rho_mat[i, primary_idx]
        for j in range(n_cols):
            v = rho_mat[i, j]; p = p_mat[i, j]
            if not np.isfinite(v):
                continue
            txt_col = "#1A1A1A" if abs(v) < 0.45 else "white"
            ax_main.text(j, i + 0.05, f"{v:+.2f}", ha="center", va="center",
                          fontsize=enforce_font_floor(FS_TICK - 0.5), color=txt_col, fontweight="bold")
            # opposite-sign marker (only for non-primary columns)
            if j != primary_idx and np.isfinite(v_primary) and v != 0 and v_primary != 0:
                if np.sign(v) != np.sign(v_primary):
                    ax_main.add_patch(mpatches.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1,
                        facecolor="none", edgecolor="#888888", linewidth=0.8,
                        linestyle=(0, (3, 2)), zorder=3,
                    ))
            if np.isfinite(p):
                if p < 0.05:
                    ax_main.add_patch(mpatches.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1,
                        facecolor="none", edgecolor=SIG_COLOR, linewidth=1.6, zorder=4,
                    ))
                    ax_main.text(j, i - 0.30, "*", ha="center", va="center",
                                  fontsize=enforce_font_floor(FS_PANEL_TITLE + 2), color=SIG_COLOR, fontweight="bold", zorder=5,
                                  path_effects=[mpe.Stroke(linewidth=1.4, foreground="#1A1A1A"),
                                                  mpe.Normal()])
                elif p < 0.10:
                    ax_main.text(j, i - 0.30, ".", ha="center", va="center",
                                  fontsize=enforce_font_floor(FS_PANEL_TITLE + 2), color="#1A1A1A", fontweight="bold", zorder=5)
    ax_main.set_xlim(-0.5, n_cols - 0.5); ax_main.set_ylim(n_rows - 0.5, -0.5)
    ax_main.set_xticks(np.arange(n_cols))
    ax_main.set_xticklabels([CODING_LABEL[c] for c in CODINGS], fontsize=enforce_font_floor(FS_TICK), fontweight="bold")
    ax_main.set_yticks([])
    ax_main.set_xlabel("Spearman \u03c1 vs AF duration  (LAA primary, sensitivity codings, LAPV region)", fontsize=enforce_font_floor(FS_AXIS_LABEL))
    for s in ax_main.spines.values():
        s.set_visible(False)

    # ---- sensitivity tally column with opposite-sign indicator
    tally_norm = mpl.colors.Normalize(vmin=0, vmax=4)
    for i, n_sens in enumerate(sens_tally):
        v_primary = rho_mat[i, primary_idx]
        # count number of non-primary codings flipped sign vs primary
        n_flip = 0
        for j in range(n_cols):
            if j == primary_idx:
                continue
            v = rho_mat[i, j]
            if np.isfinite(v) and np.isfinite(v_primary) and v != 0 and v_primary != 0:
                if np.sign(v) != np.sign(v_primary):
                    n_flip += 1
        ax_tally.add_patch(mpatches.Rectangle(
            (0, i), 1, 1, facecolor=tally_cmap(tally_norm(n_sens)),
            edgecolor="white", linewidth=0.4,
        ))
        txt_col = "#1A1A1A" if n_sens < 3 else "white"
        if n_flip > 0:
            ax_tally.text(0.5, i + 0.42, f"{int(n_sens)}/4", ha="center", va="center",
                          fontsize=enforce_font_floor(FS_TICK - 0.5), color=txt_col, fontweight="bold")
            ax_tally.text(0.5, i + 0.74, f"\u2193{n_flip}", ha="center", va="center",
                          fontsize=enforce_font_floor(FS_TICK - 1.5), color=txt_col)
        else:
            ax_tally.text(0.5, i + 0.5, f"{int(n_sens)}/4", ha="center", va="center",
                          fontsize=enforce_font_floor(FS_TICK - 0.5), color=txt_col, fontweight="bold")
    ax_tally.set_xlim(0, 1); ax_tally.set_ylim(n_rows, 0)
    ax_tally.set_xticks([0.5])
    ax_tally.set_xticklabels(["sens\n(\u2193 = flips)"], fontsize=enforce_font_floor(FS_TICK - 0.5), color="#444444")
    ax_tally.set_yticks([])
    ax_tally.tick_params(axis="x", length=0, pad=2)
    for s in ax_tally.spines.values():
        s.set_visible(False)

    # ---- legend column
    ax_leg.set_xlim(0, 1); ax_leg.set_ylim(1, 0); ax_leg.axis("off")
    cax = fig.add_axes([0.92, 0.74, 0.018, 0.14])
    cbar = fig.colorbar(im, cax=cax, orientation="vertical")
    cbar.set_label("Spearman \u03c1", fontsize=enforce_font_floor(FS_TICK))
    cbar.ax.tick_params(labelsize=enforce_font_floor(FS_TICK - 0.5))
    cbar.set_ticks([VMIN, -0.3, 0, 0.3, VMAX])
    cbar.ax.axhline(0, color="#1A1A1A", lw=0.8, zorder=10)
    # Add a small triangle marker on the colorbar at rho=0 for unambiguous zero
    cbar.ax.annotate("", xy=(1.6, 0), xycoords=("axes fraction", "data"),
                     xytext=(2.3, 0), textcoords=("axes fraction", "data"),
                     arrowprops=dict(arrowstyle="-|>", color="#1A1A1A", lw=0.7),
                     annotation_clip=False)

    ax_leg.text(0.0, 0.42, "Markers", fontsize=enforce_font_floor(FS_TICK + 0.5), fontweight="bold", color="#1A1A1A")
    ax_leg.text(0.0, 0.475, "*  p < 0.05  (yellow box)", fontsize=enforce_font_floor(FS_TICK - 0.5), color="#1A1A1A")
    ax_leg.text(0.0, 0.515, ".  p < 0.10", fontsize=enforce_font_floor(FS_TICK - 0.5), color="#444444")
    ax_leg.text(0.0, 0.555, "dashed border:\nsign-flip vs primary", fontsize=enforce_font_floor(FS_TICK - 0.5), color="#666666")

    ax_leg.text(0.0, 0.63, "Bias prefix", fontsize=enforce_font_floor(FS_TICK + 0.5), fontweight="bold")
    ax_leg.text(0.0, 0.665, "[L] = LA-biased family", fontsize=enforce_font_floor(FS_TICK - 0.5), color=LA_DEEP, fontweight="bold")
    ax_leg.text(0.0, 0.695, "[R] = RA-biased family", fontsize=enforce_font_floor(FS_TICK - 0.5), color=RA_DEEP, fontweight="bold")

    ax_leg.text(0.0, 0.76, "AF state (Fig 2C)", fontsize=enforce_font_floor(FS_TICK), fontweight="bold")
    ly = 0.79
    for st in ["Preserved", "Attenuated", "Reversed", "Amplified"]:
        ax_leg.add_patch(mpatches.Rectangle((0.0, ly), 0.08, 0.024, facecolor=STATE_COLORS[st], edgecolor="none"))
        ax_leg.text(0.12, ly + 0.012, f"{STATE_SHORT[st]} = {st}", fontsize=enforce_font_floor(FS_TICK - 0.5), va="center")
        ly += 0.028

    ax_leg.text(0.0, 0.92, "Carrier (Fig 1F)", fontsize=enforce_font_floor(FS_TICK), fontweight="bold")
    ly = 0.95
    present = list(dict.fromkeys([carriers[fid] for fid in fam_order]))
    for c in present:
        if c == MIXED_LABEL:
            ax_leg.add_patch(mpatches.Rectangle((0.0, ly), 0.08, 0.024,
                                                facecolor=MIXED_COLOR, edgecolor="#888", linewidth=0.4, hatch="///"))
            ax_leg.text(0.12, ly + 0.012, "Mixed", fontsize=enforce_font_floor(FS_TICK - 0.5), va="center")
        else:
            ax_leg.add_patch(mpatches.Rectangle((0.0, ly), 0.08, 0.024,
                                                facecolor=COMPARTMENT_PALETTE[c], edgecolor="none"))
            ax_leg.text(0.12, ly + 0.012, COMPARTMENT_SHORT[c], fontsize=enforce_font_floor(FS_TICK - 0.5), va="center")
        ly += 0.025

    fig.suptitle(
        "Fig 4b  Family-level Spearman \u03c1 vs AF duration across 4 sensitivity codings (GSE41177 LAA + LAPV)",
        x=0.04, y=0.965, ha="left", fontsize=enforce_font_floor(FS_PANEL_TITLE + 1), fontweight="bold",
    )
    fig.text(
        0.04, 0.110,
        "Spearman \u03c1 vs AF duration (months), n = 12.  * p < 0.05, . p < 0.10  (10,000 donor permutations).  "
        "Right: same-sign concordance / 4 sensitivity codings (\u2193 = number of codings flipping sign vs primary).",
        fontsize=enforce_font_floor(FS_FOOTNOTE), color="#444444",
    )
    fig.text(
        0.04, 0.075,
        "Row order: atlas PriorityScore descending (cross-figure consistency with Fig 3b).  Family bias is encoded as [L] / [R] "
        "prefix in the family label; carrier strip = HCA baseline (Fig 1F);  state strip = Fig 2C consensus.",
        fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5), color="#444444",
    )
    fig.text(
        0.04, 0.040,
        "Caveats:  (i) F15 \u03c1 = -0.385 in primary (rounded -0.39 < 0.4 \u2192 Weak-trend, NOT Strong-trend).  "
        "(ii) F11 BMP gt1m_eq1 \u03c1 = -0.43 (p<0.05) but primary \u03c1 = -0.20 (NS) -- single-coding sig is not robust here.",
        fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5), color="#5A5A5A", style="italic",
    )
    fig.text(
        0.04, 0.020,
        "(iii) F3 TG metabolism shows 4/4 concordant positive \u03c1 (lipid metabolic adaptation strengthening with chronicity, distinct from F9/F15 attenuation).",
        fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5), color="#5A5A5A", style="italic",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(OUT_PDF.relative_to(ROOT)))


if __name__ == "__main__":
    main()
