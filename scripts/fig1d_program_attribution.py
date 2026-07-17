#!/usr/bin/env python3
"""
Figure 1D - Program x cell-type attribution heatmap (signed + significance).

Design (agreed final rules)
---------------------------
Color (signed attribution score per program x cell type):
    signed_score = sign(sc_direction vs bulk_direction)
                 * |sc_NES|
                 * LE_Jaccard
  * positive  -> this cell type supports the bulk program's direction
  * negative  -> this cell type reproduces the program in the OPPOSITE direction
  * zero      -> missing value (not a significance gate; significance is
                 overlaid separately)

Significance marker (overlaid text):
    padj_ct < 0.05  -> '*'
    padj_ct < 0.10  -> '.'
    else            -> ''

Rows match 1C exactly (signed bulk NES descending).
Columns: Atrial CM, Fibroblast, Endothelial / Vascular, Immune.

Diverging colormap (#3B73B9 - white - #C74A44) with a TwoSlopeNorm centered
at 0.  Cells are drawn as vector Rectangle patches so Adobe Illustrator
renders RGB fills faithfully.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

from plot_style import (
    COL_ONEHALF, LA_COLOR, RA_COLOR, set_nature_style,
)

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "figures" / "fig1d_program_attribution.pdf"

CT_KEYS   = ["atrialCM", "fibro", "endo", "immune"]
CT_LABELS = ["Atrial CM", "Fibroblast", "Endothelial /\nVascular", "Immune"]

DIV_CMAP = LinearSegmentedColormap.from_list(
    "signed_attr", ["#3B73B9", "#FFFFFF", "#C74A44"]
)


def load_row_order() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "normal_chamber_programs.tsv", sep="\t")
    for c in ("program", "program_id", "direction"):
        df[c] = df[c].astype(str).str.strip('"')
    df["signed_NES"] = df["NES"]
    df = df.sort_values("signed_NES", ascending=False).reset_index(drop=True)
    return df[["program", "program_id", "direction", "signed_NES", "NES"]]


def load_attribution() -> pd.DataFrame:
    return pd.read_csv(ROOT / "results/normal_chamber_program_attribution.tsv",
                       sep="\t")


def compute_signed_score_matrix(
    row_order: pd.DataFrame, attr: pd.DataFrame, ct_keys: list[str],
) -> np.ndarray:
    """signed_score = sign(sc_dir == bulk_dir) * |NES| * jac."""
    attr_by_id = attr.set_index("program_id")
    mat = np.zeros((len(row_order), len(ct_keys)), dtype=float)

    for i, prow in row_order.iterrows():
        arow = attr_by_id.loc[prow["program_id"]]
        bulk_dir = "LA_up" if str(prow["direction"]).startswith("LA") else "RA_up"

        for j, ct in enumerate(ct_keys):
            sc_dir = arow.get(f"sc_dir_{ct}", np.nan)
            nes    = arow.get(f"NES_{ct}", np.nan)
            jac    = arow.get(f"LE_overlap_jac_{ct}", np.nan)

            if pd.isna(sc_dir) or pd.isna(nes) or pd.isna(jac):
                continue

            sign = 1.0 if str(sc_dir) == bulk_dir else -1.0
            mat[i, j] = sign * abs(float(nes)) * float(jac)

    return mat


def compute_sig_label_matrix(
    row_order: pd.DataFrame, attr: pd.DataFrame, ct_keys: list[str],
) -> np.ndarray:
    """'*' padj<0.05 ; '.' padj<0.10 ; '' otherwise."""
    attr_by_id = attr.set_index("program_id")
    lab = np.empty((len(row_order), len(ct_keys)), dtype=object)
    lab[:] = ""

    for i, prow in row_order.iterrows():
        arow = attr_by_id.loc[prow["program_id"]]
        for j, ct in enumerate(ct_keys):
            padj = arow.get(f"padj_{ct}", np.nan)
            if pd.isna(padj):
                continue
            p = float(padj)
            if p < 0.05:
                lab[i, j] = "*"
            elif p < 0.10:
                lab[i, j] = "."
    return lab


def main() -> None:
    set_nature_style()
    row_order = load_row_order()
    attr      = load_attribution()
    M    = compute_signed_score_matrix(row_order, attr, CT_KEYS)
    SIG  = compute_sig_label_matrix(row_order, attr, CT_KEYS)
    n_rows, n_cols = M.shape

    fig_w = COL_ONEHALF
    fig_h = 0.28 * n_rows + 0.9
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.44, right=0.94, top=0.92, bottom=0.24)

    # Symmetric two-slope norm centered at 0 so white = 0.
    amax = float(np.max(np.abs(M)))
    amax = max(amax, 1e-6)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-amax, vmax=+amax)

    # Vector cell rectangles.
    for i in range(n_rows):
        for j in range(n_cols):
            rgba = DIV_CMAP(norm(M[i, j]))
            ax.add_patch(Rectangle(
                (j - 0.5, i - 0.5), 1.0, 1.0,
                facecolor=rgba, edgecolor="white", linewidth=0.5,
            ))

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_aspect("auto")

    # Rows / columns.
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(CT_LABELS)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_order["program"].tolist())
    ax.tick_params(axis="both", length=0)
    ax.tick_params(axis="x", pad=2)
    ax.tick_params(axis="y", pad=3)

    # Color row labels by bulk direction.
    for tick_label, direction in zip(
        ax.get_yticklabels(), row_order["direction"]
    ):
        tick_label.set_color(LA_COLOR if direction.startswith("LA") else RA_COLOR)

    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    # Divider between LA and RA blocks.
    n_LA = int(row_order["direction"].str.startswith("LA").sum())
    if 0 < n_LA < n_rows:
        ax.axhline(n_LA - 0.5, color="black", lw=0.8)

    # Significance markers (no numeric scores).
    for i in range(n_rows):
        for j in range(n_cols):
            s = SIG[i, j]
            if not s:
                continue
            # '.' rides a bit higher than '*' so it looks centered.
            va_offset = -0.08 if s == "*" else -0.22
            ax.text(j, i + va_offset, s,
                    ha="center", va="center",
                    fontsize=9, fontweight="bold", color="black")

    ax.set_title("Cell-compartment attribution of healthy chamber programs",
                 fontsize=7.6, pad=4)

    # Colorbar.
    cax = fig.add_axes([0.48, 0.07, 0.42, 0.025])
    sm = ScalarMappable(norm=norm, cmap=DIV_CMAP)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label("signed attribution score", labelpad=2)
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(width=0.4, length=2, pad=1, labelsize=6)

    # Flank labels for the colorbar.
    cb_y_mid = 0.07 + 0.025 / 2
    fig.text(0.48 - 0.013, cb_y_mid, "opposite",
             color="#3B73B9", fontsize=6, fontweight="bold",
             ha="right", va="center")
    fig.text(0.48 + 0.42 + 0.013, cb_y_mid, "concordant",
             color="#C74A44", fontsize=6, fontweight="bold",
             ha="left", va="center")

    # Significance legend + direction legend (row-colour key).
    fig.text(0.06, 0.10, "* padj < 0.05",
             fontsize=6, ha="left", va="center")
    fig.text(0.06, 0.075, ". padj < 0.10",
             fontsize=6, ha="left", va="center")

    fig.text(0.06, 0.035, "LA-high",
             color=LA_COLOR, fontsize=6, fontweight="bold", ha="left")
    fig.text(0.17, 0.035, "RA-high",
             color=RA_COLOR, fontsize=6, fontweight="bold", ha="left")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")

    # Diagnostics: print the signed-score matrix + sig labels
    disp_score = pd.DataFrame(M,   index=row_order["program"], columns=CT_KEYS)
    disp_sig   = pd.DataFrame(SIG, index=row_order["program"], columns=CT_KEYS)
    print("\nSigned attribution score matrix:")
    print(disp_score.round(3).to_string())
    print("\nSignificance marker matrix:")
    print(disp_sig.to_string())


if __name__ == "__main__":
    main()
