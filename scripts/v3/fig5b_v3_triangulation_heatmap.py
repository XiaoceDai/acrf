#!/usr/bin/env python3
"""
Fig 5b v3 - Cross-figure triangulation heatmap (16 NCP_v3 families).

Stage-grouped layout (left -> right = narrative time arrow):

    [F-id + short]
      L0 stage:  HCA LA carrier  (cell-type colour band)
      L1 stage:  bulk pooled g (continuous, BrBG)  +  state band (Pre/Att/Rev/Amp)
      L2 stage:  snRNA dScore at carrier compartment (continuous, PRGn)
                 +  fidelity band (Fid/Relay/Ext/Below)
      L3 stage:  duration rho_primary (continuous, BrBG)
                 +  chronicity band (Strong/Weak/Stable)
      Atlas:     PriorityScore (horizontal bar, with star marker)

Continuous cells are labelled with the numeric value; discrete bands are
solid colour with no in-cell text (legend on the right).  Cells that are
"supportive" by the strict v3 criteria carry a yellow border:
    L1 supportive  =  replication_class == high       (CI excludes 0)
    L2 supportive  =  fidelity_class != Below_detection
    L3 supportive  =  chronicity_class != Stable

Row order:  Atlas PriorityScore descending  (cross-figure consistency).

Inputs:
    results/fig_integrated_af_rewiring_atlas.tsv     (row order, anchor)
    results/fig1f_v3_HCA_family_carriers.tsv         (LA carrier)
    results/fig2e_replication_lollipop_v3.tsv        (pooled g, repl class, state)
    results/fig3_v3_family_celltype_delta.tsv        (snRNA dScore per ct)
    results/fig3_v3_carrier_fidelity_class.tsv       (fidelity class)
    results/fig4_v3_chronicity_class.tsv             (rho_primary, chronicity)

Outputs:
    results/fig5b_v3_triangulation.tsv               (per-family rolled-up table)
    figures/fig5b_v3_triangulation_heatmap.pdf       (the heatmap)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

ATLAS_TSV = ROOT / "results" / "fig_integrated_af_rewiring_atlas.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
REPL_TSV = ROOT / "results" / "fig2e_replication_lollipop_v3.tsv"
DELTA_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
FID_TSV = ROOT / "results" / "fig3_v3_carrier_fidelity_class.tsv"
CHRO_TSV = ROOT / "results" / "fig4_v3_chronicity_class.tsv"

OUT_TSV = ROOT / "results" / "fig5b_v3_triangulation.tsv"
OUT_PDF = ROOT / "figures" / "fig5b_v3_triangulation_heatmap.pdf"

# --- shared cell-type colour palette (matches Fig 1F-v3 / Fig 3B) ---
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
    "Mixed/Unresolved": "#BDBDBD",
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
    "Mixed/Unresolved": "Mix",
}

# --- discrete band palettes ---
STATE_COLORS = {
    "Preserved": "#9C9C9C",
    "Attenuated": "#A6CEE3",
    "Reversed": "#762A83",
    "Amplified": "#1B7837",
}
FIDELITY_COLORS = {
    "Fidelity": "#7FB069",
    "Relay": "#8E6BAA",
    "Extension": "#5F9EA0",
    "Below_detection": "#F0F0F0",
}
CHRONICITY_COLORS = {
    "Strong_trend": "#C44569",
    "Weak_trend": "#F8B500",
    "Stable": "#E8E8E8",
}

# --- continuous palettes (commensurate within stage, not across) ---
G_VMAX = 1.5      # |pooled g|
DSCORE_VMAX = 0.04
RHO_VMAX = 0.6

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
SUPPORT_BORDER = "#F5C518"  # yellow border for supportive cells


def diverging_cmap():
    return mcolors.LinearSegmentedColormap.from_list(
        "diverging_v3",
        [(0.0, RA_DEEP), (0.45, "#FFFFFF"), (0.55, "#FFFFFF"), (1.0, LA_DEEP)],
    )


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def value_to_color(v: float, vmax: float, cmap) -> tuple:
    if not np.isfinite(v):
        return (0.92, 0.92, 0.92, 1.0)
    norm = mcolors.Normalize(vmin=-vmax, vmax=vmax, clip=True)
    return cmap(norm(v))


def build_table() -> pd.DataFrame:
    atlas = pd.read_csv(ATLAS_TSV, sep="\t")
    hca = pd.read_csv(HCA_TSV, sep="\t")
    repl = pd.read_csv(REPL_TSV, sep="\t")
    delta = pd.read_csv(DELTA_TSV, sep="\t")
    fid = pd.read_csv(FID_TSV, sep="\t")
    chro = pd.read_csv(CHRO_TSV, sep="\t")

    atlas = atlas.sort_values("PriorityScore", ascending=False).reset_index(drop=True)
    rows = []
    for _, ar in atlas.iterrows():
        fid_id = int(ar["family_id"])
        # L0 baseline carrier (cell-type)
        hca_row = hca[hca["family_id"] == fid_id].iloc[0]
        la_carrier = str(hca_row["LA_carrier"])
        # L1 bulk
        repl_row = repl[repl["family_id"] == fid_id]
        if len(repl_row):
            repl_row = repl_row.iloc[0]
            pooled_g = float(repl_row["pooled_g"])
            repl_class = str(repl_row["replication_class"])
            state_consensus = str(repl_row["dominant_state"])
        else:
            pooled_g = float("nan"); repl_class = "low"; state_consensus = "Preserved"
        # L2 snRNA dScore at LA carrier compartment
        if la_carrier == "Mixed/Unresolved":
            d_at_carrier = float("nan")
            p_at_carrier = float("nan")
        else:
            sub = delta[(delta["family_id"] == fid_id) & (delta["compartment"] == la_carrier)]
            if len(sub):
                d_at_carrier = float(sub["dScore"].iloc[0])
                p_at_carrier = float(sub["p_perm"].iloc[0])
            else:
                d_at_carrier = float("nan"); p_at_carrier = float("nan")
        # L2' fidelity
        fid_row = fid[fid["family_id"] == fid_id]
        fidelity_class = str(fid_row.iloc[0]["fidelity_class"]) if len(fid_row) else "Below_detection"
        # L3 chronicity
        chro_row = chro[chro["family_id"] == fid_id]
        if len(chro_row):
            chro_row = chro_row.iloc[0]
            rho_primary = float(chro_row["rho_primary"])
            p_chro = float(chro_row["p_perm_primary"])
            chronicity_class = str(chro_row["chronicity_class"])
        else:
            rho_primary = float("nan"); p_chro = float("nan"); chronicity_class = "Stable"
        rows.append({
            "family_id": fid_id,
            "short": str(ar["short"]),
            "HealthyBias": str(ar["HealthyBias"]),
            "PriorityScore": float(ar["PriorityScore"]),
            "PriorityClass": str(ar["PriorityClass"]),
            "PriorityStars": str(ar["PriorityStars"]),
            "ChamberDriver": str(ar["ChamberDriver"]),
            # L0
            "LA_carrier": la_carrier,
            # L1
            "pooled_g": pooled_g,
            "replication_class": repl_class,
            "L1_supportive": int(repl_class == "high"),
            "state_consensus": state_consensus,
            # L2
            "dScore_at_carrier": d_at_carrier,
            "p_perm_carrier": p_at_carrier,
            "fidelity_class": fidelity_class,
            "L2_supportive": int(fidelity_class != "Below_detection"),
            # L3
            "rho_primary": rho_primary,
            "p_perm_primary": p_chro,
            "chronicity_class": chronicity_class,
            "L3_supportive": int(chronicity_class != "Stable"),
        })
    df = pd.DataFrame(rows)
    df["n_supportive"] = df[["L1_supportive", "L2_supportive", "L3_supportive"]].sum(axis=1)
    return df


def render_pdf(df: pd.DataFrame) -> None:
    set_style()
    n = len(df)
    MM = 1.0 / 25.4
    fig_w_mm, fig_h_mm = 230, 155
    fig = plt.figure(figsize=(fig_w_mm * MM, fig_h_mm * MM))
    ax = fig.add_axes((0.005, 0.04, 0.99, 0.84))
    ax.set_xlim(0, 100)
    # Reserve 2 axes-units above row 0 for two header rows (stage band + col names)
    ax.set_ylim(-2.2, n + 1)
    ax.invert_yaxis()
    ax.axis("off")

    # Column geometry (in axes-coord units, 0..100)
    X = {
        "fid":      (1.0,  4.0),     # F15 etc
        "short":    (4.5,  20.0),    # short label
        "L0":       (20.5, 23.5),    # carrier strip
        "L1_g":     (24.0, 31.5),    # bulk g cell
        "L1_state": (31.7, 34.7),    # state strip
        "L2_d":     (35.2, 42.7),    # snRNA dScore cell
        "L2_fid":   (42.9, 45.9),    # fidelity strip
        "L3_rho":   (46.4, 53.9),    # duration rho cell
        "L3_chro":  (54.1, 57.1),    # chronicity strip
        "atlas":    (57.6, 80.0),    # priority score bar
        "stars":    (80.4, 84.0),    # priority stars
        "n_sup":    (84.5, 87.5),    # number of supportive layers
    }

    cmap = diverging_cmap()
    row_h = 0.92  # height of each row inside the unit-spaced grid

    # ---- header rows ----
    # Row A (stage banner) at y = -1.55, Row B (column names) at y = -0.55
    stage_y = -1.55
    col_y = -0.30

    def hdr_col(x_centre, txt, *, fontsize=5.7, weight="bold"):
        ax.text(x_centre, col_y, txt, ha="center", va="bottom",
                fontsize=fontsize, fontweight=weight, color="#1A1A1A")

    # Stage banner backgrounds (light pastel rectangles to make stages visually distinct)
    stage_bg = "#F5F5F5"

    def stage_band(x0, x1, label):
        ax.add_patch(mpatches.FancyBboxPatch(
            (x0, stage_y - 0.40), x1 - x0, 0.85,
            boxstyle="round,pad=0.02,rounding_size=0.10",
            facecolor=stage_bg, edgecolor="#999", linewidth=0.4, zorder=0.5,
        ))
        ax.text((x0 + x1) / 2.0, stage_y, label,
                ha="center", va="center", fontsize=6.6,
                fontweight="bold", color="#1A1A1A")

    # Family / L0 stays as a "anchor" header (no band)
    ax.text((X["fid"][0] + X["short"][1]) / 2.0, stage_y,
              "Family / Atlas anchor", ha="center", va="center",
              fontsize=6.6, fontweight="bold", color="#1A1A1A")

    stage_band(X["L0"][0] - 0.2, X["L0"][1] + 0.2, "L0  Carrier")
    stage_band(X["L1_g"][0] - 0.2, X["L1_state"][1] + 0.2, "L1  Bulk (Fig 2)")
    stage_band(X["L2_d"][0] - 0.2, X["L2_fid"][1] + 0.2, "L2  snRNA (Fig 3)")
    stage_band(X["L3_rho"][0] - 0.2, X["L3_chro"][1] + 0.2, "L3  Duration (Fig 4)")
    stage_band(X["atlas"][0] - 0.2, X["n_sup"][1] + 0.2, "Atlas integration")

    # Column names (Row B)
    hdr_col((X["fid"][0] + X["short"][1]) / 2.0, "Family")
    hdr_col((X["L0"][0] + X["L0"][1]) / 2.0, "HCA carrier")
    hdr_col((X["L1_g"][0] + X["L1_g"][1]) / 2.0, "pooled g")
    hdr_col((X["L1_state"][0] + X["L1_state"][1]) / 2.0, "state")
    hdr_col((X["L2_d"][0] + X["L2_d"][1]) / 2.0, "\u0394 at carrier")
    hdr_col((X["L2_fid"][0] + X["L2_fid"][1]) / 2.0, "fidelity")
    hdr_col((X["L3_rho"][0] + X["L3_rho"][1]) / 2.0, "\u03c1 (LAA)")
    hdr_col((X["L3_chro"][0] + X["L3_chro"][1]) / 2.0, "chron")
    hdr_col((X["atlas"][0] + X["atlas"][1]) / 2.0, "PriorityScore")
    hdr_col((X["stars"][0] + X["stars"][1]) / 2.0, "Class")
    hdr_col((X["n_sup"][0] + X["n_sup"][1]) / 2.0, "n.sup")

    # ---- vertical light separators between stages ----
    for x in [X["L0"][1] + 0.25, X["L1_state"][1] + 0.25,
              X["L2_fid"][1] + 0.25, X["L3_chro"][1] + 0.25,
              X["atlas"][1] + 0.20]:
        ax.plot([x, x], [-0.1, n + 0.05], color="#CCCCCC", lw=0.5, zorder=1)

    # ---- per-row rendering ----
    for i, row in df.iterrows():
        y = i + 0.5  # cell centre
        y0 = i + 0.5 - row_h / 2.0
        # --- F-id + short label
        bias_color = LA_DEEP if row["HealthyBias"] == "LA-biased" else RA_DEEP
        ax.text(X["fid"][0], y, f"F{int(row['family_id'])}", ha="left", va="center",
                fontsize=6.8, fontweight="bold", color=bias_color)
        ax.text(X["short"][0] + 0.5, y, str(row["short"])[:24],
                ha="left", va="center", fontsize=6.4, color="#1A1A1A")

        # --- L0 carrier strip (cell-type colour) ---
        carrier = row["LA_carrier"]
        c_color = COMPARTMENT_PALETTE.get(carrier, "#BDBDBD")
        rect = mpatches.Rectangle((X["L0"][0], y0), X["L0"][1] - X["L0"][0], row_h,
                                    facecolor=c_color, edgecolor="#888", linewidth=0.3)
        ax.add_patch(rect)
        if carrier == "Mixed/Unresolved":
            # diagonal hatch overlay
            rect.set_hatch("///")
            rect.set_edgecolor("#666")
        # carrier short label
        ax.text((X["L0"][0] + X["L0"][1]) / 2.0, y,
                COMPARTMENT_SHORT.get(carrier, "?"),
                ha="center", va="center", fontsize=5.4, color="#1A1A1A")

        # --- L1 g value cell ---
        g = row["pooled_g"]
        gc = value_to_color(g, G_VMAX, cmap)
        edge = SUPPORT_BORDER if row["L1_supportive"] else "#888"
        lw_edge = 1.2 if row["L1_supportive"] else 0.3
        rect = mpatches.Rectangle((X["L1_g"][0], y0), X["L1_g"][1] - X["L1_g"][0], row_h,
                                    facecolor=gc, edgecolor=edge, linewidth=lw_edge)
        ax.add_patch(rect)
        if np.isfinite(g):
            txt_color = "white" if abs(g) > G_VMAX * 0.6 else "#1A1A1A"
            ax.text((X["L1_g"][0] + X["L1_g"][1]) / 2.0, y, f"{g:+.2f}",
                    ha="center", va="center", fontsize=6.0, color=txt_color)

        # --- L1' state strip ---
        s_color = STATE_COLORS.get(row["state_consensus"], "#FFFFFF")
        rect = mpatches.Rectangle((X["L1_state"][0], y0),
                                    X["L1_state"][1] - X["L1_state"][0], row_h,
                                    facecolor=s_color, edgecolor="#888", linewidth=0.3)
        ax.add_patch(rect)

        # --- L2 dScore cell ---
        d = row["dScore_at_carrier"]
        if not np.isfinite(d):
            rect = mpatches.Rectangle((X["L2_d"][0], y0),
                                        X["L2_d"][1] - X["L2_d"][0], row_h,
                                        facecolor="#F0F0F0",
                                        edgecolor="#AAA", linewidth=0.3, hatch="///")
            ax.add_patch(rect)
            ax.text((X["L2_d"][0] + X["L2_d"][1]) / 2.0, y, "—",
                    ha="center", va="center", fontsize=6.0, color="#888")
        else:
            dc = value_to_color(d, DSCORE_VMAX, cmap)
            edge = SUPPORT_BORDER if row["L2_supportive"] else "#888"
            lw_edge = 1.2 if row["L2_supportive"] else 0.3
            rect = mpatches.Rectangle((X["L2_d"][0], y0),
                                        X["L2_d"][1] - X["L2_d"][0], row_h,
                                        facecolor=dc, edgecolor=edge, linewidth=lw_edge)
            ax.add_patch(rect)
            txt_color = "white" if abs(d) > DSCORE_VMAX * 0.6 else "#1A1A1A"
            ax.text((X["L2_d"][0] + X["L2_d"][1]) / 2.0, y, f"{d:+.3f}",
                    ha="center", va="center", fontsize=6.0, color=txt_color)
            # significance dot/star inside the cell, top-right
            if np.isfinite(row["p_perm_carrier"]):
                if row["p_perm_carrier"] < 0.05:
                    ax.text(X["L2_d"][1] - 0.4, y0 + 0.15, "*",
                            ha="right", va="top", fontsize=7.0, color="#1A1A1A")
                elif row["p_perm_carrier"] < 0.10:
                    ax.text(X["L2_d"][1] - 0.4, y0 + 0.15, ".",
                            ha="right", va="top", fontsize=7.0, color="#1A1A1A")

        # --- L2' fidelity strip ---
        f_color = FIDELITY_COLORS.get(row["fidelity_class"], "#F0F0F0")
        edge = SUPPORT_BORDER if row["L2_supportive"] else "#888"
        lw_edge = 1.2 if row["L2_supportive"] else 0.3
        rect = mpatches.Rectangle((X["L2_fid"][0], y0),
                                    X["L2_fid"][1] - X["L2_fid"][0], row_h,
                                    facecolor=f_color, edgecolor=edge, linewidth=lw_edge)
        ax.add_patch(rect)
        # If fidelity is Relay or Extension AND carrier-Δ is small, show
        # a tiny annotation in the L2 dScore cell pointing to where the
        # actual signal lives, so readers don't misread Δ ≈ 0.
        fid_cls = str(row["fidelity_class"])
        if fid_cls in ("Relay", "Extension") and np.isfinite(row["dScore_at_carrier"]):
            if abs(row["dScore_at_carrier"]) < DSCORE_VMAX * 0.30:
                tag = "Relay\u2192non-carrier" if fid_cls == "Relay" else "Extension"
                ax.text(X["L2_d"][1] - 0.15, y0 + row_h - 0.05, tag,
                          ha="right", va="top",
                          fontsize=4.6, color="#5E3C99", style="italic",
                          fontweight="bold")

        # --- L3 rho cell ---
        rho = row["rho_primary"]
        rc = value_to_color(rho, RHO_VMAX, cmap)
        edge = SUPPORT_BORDER if row["L3_supportive"] else "#888"
        lw_edge = 1.2 if row["L3_supportive"] else 0.3
        rect = mpatches.Rectangle((X["L3_rho"][0], y0),
                                    X["L3_rho"][1] - X["L3_rho"][0], row_h,
                                    facecolor=rc, edgecolor=edge, linewidth=lw_edge)
        ax.add_patch(rect)
        if np.isfinite(rho):
            txt_color = "white" if abs(rho) > RHO_VMAX * 0.6 else "#1A1A1A"
            ax.text((X["L3_rho"][0] + X["L3_rho"][1]) / 2.0, y, f"{rho:+.2f}",
                    ha="center", va="center", fontsize=6.0, color=txt_color)
            if np.isfinite(row["p_perm_primary"]):
                if row["p_perm_primary"] < 0.05:
                    ax.text(X["L3_rho"][1] - 0.4, y0 + 0.15, "*",
                            ha="right", va="top", fontsize=7.0, color="#1A1A1A")
                elif row["p_perm_primary"] < 0.10:
                    ax.text(X["L3_rho"][1] - 0.4, y0 + 0.15, ".",
                            ha="right", va="top", fontsize=7.0, color="#1A1A1A")

        # --- L3' chronicity strip ---
        chro_color = CHRONICITY_COLORS.get(row["chronicity_class"], "#F0F0F0")
        rect = mpatches.Rectangle((X["L3_chro"][0], y0),
                                    X["L3_chro"][1] - X["L3_chro"][0], row_h,
                                    facecolor=chro_color, edgecolor="#888", linewidth=0.3)
        ax.add_patch(rect)

        # --- Atlas PriorityScore bar ---
        ps = float(row["PriorityScore"])
        bar_w_full = X["atlas"][1] - X["atlas"][0]
        bar_w = bar_w_full * (ps / 100.0)
        # background track
        ax.add_patch(mpatches.Rectangle((X["atlas"][0], y0 + row_h * 0.25),
                                          bar_w_full, row_h * 0.5,
                                          facecolor="#EEEEEE", edgecolor="none"))
        # filled bar coloured by priority stars
        stars = str(row["PriorityStars"])
        bar_color = {"***": "#2E7B5E", "**": "#7BB36C", "*": "#BBBBBB"}.get(stars, "#BBBBBB")
        ax.add_patch(mpatches.Rectangle((X["atlas"][0], y0 + row_h * 0.25),
                                          bar_w, row_h * 0.5,
                                          facecolor=bar_color, edgecolor="none"))
        ax.text(X["atlas"][0] + bar_w + 0.4, y, f"{ps:.0f}",
                ha="left", va="center", fontsize=5.8, color="#1A1A1A")

        # --- stars (★ glyphs) ---
        ax.text((X["stars"][0] + X["stars"][1]) / 2.0, y,
                stars,
                ha="center", va="center",
                fontsize=8.4, fontweight="bold",
                color={"***": "#2E7B5E", "**": "#7BB36C", "*": "#888888"}.get(stars, "#888"))
        # --- n.sup ---
        ax.text((X["n_sup"][0] + X["n_sup"][1]) / 2.0, y,
                f"{int(row['n_supportive'])}/3",
                ha="center", va="center",
                fontsize=6.2, fontweight="bold",
                color="#1A1A1A")

    # ---- Title block ----
    fig.text(0.005, 0.965,
             "Fig 5b  Cross-figure triangulation  (16 NCP_v3 families)",
             fontsize=10.0, fontweight="bold", color="#1A1A1A")
    fig.text(0.005, 0.940,
             "Stage-grouped layout:  L0 HCA carrier  \u2192  L1 bulk  \u2192  L2 snRNA  \u2192  L3 duration  \u2192  Atlas",
             fontsize=6.8, color="#444444", style="italic")

    # ---- Compact bottom caption (full legend in figure caption / local manuscript) ----
    fig.text(0.005, 0.022,
             "Stage-grouped layout: L0 HCA carrier  \u2192  L1 bulk  \u2192  L2 snRNA  \u2192  L3 duration  \u2192  Atlas.  "
             "Continuous cells: BrBG palette (white = 0; deep red = positive / LA-ward; deep blue = negative).  "
             "Yellow border = layer-supportive.   *  p_perm < 0.05;   .  p_perm < 0.10.",
             fontsize=5.6, color="#444444")
    fig.text(0.005, 0.005,
             "L2 \u0394 = signed \u0394Score at the healthy LA carrier cell type;  "
             "L2 supportive (yellow border on fidelity strip) = fidelity_class \u2208 {Fidelity, Relay, Extension} "
             "(signal may reside in a non-carrier compartment, e.g. F15 Relay \u2192 lymph + endothelial).  "
             "Hatched L2 cells = LA carrier is Mixed/Unresolved (no single carrier-specific \u0394).  "
             "See Figure Legend for full colour key.",
             fontsize=5.6, color="#444444")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


def main():
    df = build_table()
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.4f")
    print(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(df)} rows)")
    render_pdf(df)


if __name__ == "__main__":
    main()
