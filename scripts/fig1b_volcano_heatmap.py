#!/usr/bin/env python3
"""
Figure 1B - GSE112339 bulk volcano + top-gene heatmap.

Rebuilt in Python to match the Cell/Nature-style format used across
Figure 1 panels (A, C, D).

Inputs:
    Project_P01/normal_chamber_genes_full.tsv   DE result (one row per gene)
    bulk/logcpm_voom.tsv                        voom-normalized logCPM matrix
    bulk/sample_meta.tsv                        curated sample metadata

Visual rules:
    LA = blue (#2166AC)   RA = red (#B2182B)
    Volcano: x = logFC (LA - RA),  y = -log10(FDR)
             LA-high dots = blue, RA-high dots = red, NS = grey
             dashed guides at |logFC|=0.5 and FDR=0.05
    Heatmap: top 20 LA + top 20 RA stable, paired-donor-ordered, row-z
             standard blue-white-red z-score palette
             column chamber strip colored LA=blue / RA=red
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
    COL_SINGLE, COL_ONEHALF, LA_COLOR, LA_FILL, RA_COLOR, RA_FILL, set_nature_style,
)

ROOT = Path(__file__).resolve().parent.parent
DE_TSV     = ROOT / "bulk/Project_P01/normal_chamber_genes_full.tsv"
EXPR_TSV   = ROOT / "bulk/logcpm_voom.tsv"
META_TSV   = ROOT / "bulk/sample_meta.tsv"
OUT_VOLC   = ROOT / "figures" / "fig1b_volcano.pdf"
OUT_HMAP   = ROOT / "figures" / "fig1b_heatmap.pdf"

FDR_CUT     = 0.05
LOGFC_CUT   = 0.5
TOP_EACH    = 15
LABEL_EACH  = 4

# z-score colormap for the heatmap: unambiguous blue->white->red (low->high).
ZCMAP = LinearSegmentedColormap.from_list(
    "zscore_bwr", ["#2166AC", "#4393C3", "#D1E5F0",
                   "#F7F7F7",
                   "#FDDBC7", "#D6604D", "#B2182B"]
)

NS_COLOR = "#BEBEBE"


def _strip_quotes(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip('"')


def load_de() -> pd.DataFrame:
    df = pd.read_csv(DE_TSV, sep="\t")
    for c in ("ensembl_id", "gene_symbol", "direction"):
        if c in df.columns:
            df[c] = _strip_quotes(df[c])
    df["label"] = df["gene_symbol"].where(
        df["gene_symbol"].notna() & (df["gene_symbol"] != ""),
        df["ensembl_id"],
    )
    df["neglog10fdr"] = -np.log10(np.clip(df["FDR"], 1e-300, None))
    df["sig_class"] = np.select(
        [
            (df["FDR"] < FDR_CUT) & (df["logFC"] > 0),
            (df["FDR"] < FDR_CUT) & (df["logFC"] < 0),
        ],
        ["LA", "RA"],
        default="NS",
    )
    return df


def load_expr_meta():
    expr = pd.read_csv(EXPR_TSV, sep="\t")
    for c in ("gene_key", "gene_symbol"):
        expr[c] = _strip_quotes(expr[c])
    expr = expr.set_index("gene_key")

    meta = pd.read_csv(META_TSV, sep="\t")
    for c in meta.columns:
        meta[c] = _strip_quotes(meta[c])
    return expr, meta


def draw_volcano(ax: plt.Axes, de: pd.DataFrame) -> None:
    color_map = {"LA": LA_FILL, "RA": RA_FILL, "NS": NS_COLOR}
    for cls in ("NS", "LA", "RA"):  # draw NS first so sig points sit on top
        sub = de[de["sig_class"] == cls]
        ax.scatter(
            sub["logFC"], sub["neglog10fdr"],
            s=16, facecolor=color_map[cls], edgecolor="none",
            alpha=0.85 if cls != "NS" else 0.55,
            linewidth=0,
        )

    ax.axvline(-LOGFC_CUT, color="#808080", lw=0.4, ls=(0, (3, 2)))
    ax.axvline(+LOGFC_CUT, color="#808080", lw=0.4, ls=(0, (3, 2)))
    ax.axhline(-np.log10(FDR_CUT), color="#808080", lw=0.4, ls=(0, (3, 2)))

    # Top genes (by moderated_t) on each side get text labels, with a
    # tiny vertical stagger so they don't print on top of each other.
    la_lab = (
        de[(de["FDR"] < FDR_CUT) & (de["logFC"] > LOGFC_CUT)]
        .nlargest(LABEL_EACH, "moderated_t")
        .reset_index(drop=True)
    )
    ra_lab = (
        de[(de["FDR"] < FDR_CUT) & (de["logFC"] < -LOGFC_CUT)]
        .nsmallest(LABEL_EACH, "moderated_t")
        .reset_index(drop=True)
    )
    for i, row in la_lab.iterrows():
        dx = 6
        dy = 4 + (i % 2) * 5   # zig-zag
        ax.annotate(row["label"],
                    xy=(row["logFC"], row["neglog10fdr"]),
                    xytext=(dx, dy), textcoords="offset points",
                    fontsize=5, ha="left", va="bottom",
                    color="#222222",
                    arrowprops=dict(arrowstyle="-", color="#808080",
                                    lw=0.3, shrinkA=0, shrinkB=1))
    for i, row in ra_lab.iterrows():
        dx = -6
        dy = 4 + (i % 2) * 5
        ax.annotate(row["label"],
                    xy=(row["logFC"], row["neglog10fdr"]),
                    xytext=(dx, dy), textcoords="offset points",
                    fontsize=5, ha="right", va="bottom",
                    color="#222222",
                    arrowprops=dict(arrowstyle="-", color="#808080",
                                    lw=0.3, shrinkA=0, shrinkB=1))

    ax.set_xlabel("logFC (LA - RA)")
    ax.set_ylabel(r"$-\log_{10}(\mathrm{FDR})$")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # Legend swatches.
    handles = [
        plt.Line2D([0], [0], marker='o', linestyle='', color=LA_FILL,
                   markersize=5.5, label="LA-high"),
        plt.Line2D([0], [0], marker='o', linestyle='', color=RA_FILL,
                   markersize=5.5, label="RA-high"),
        plt.Line2D([0], [0], marker='o', linestyle='', color=NS_COLOR,
                   markersize=5.5, label="NS"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False,
              handletextpad=0.3, borderpad=0.2, fontsize=5.5)

    ax.set_title("b1  Volcano", fontsize=7.6, pad=4)


def pick_top_keys(de: pd.DataFrame, expr_index: pd.Index) -> pd.DataFrame:
    """Pick top LA + RA genes (by |moderated_t|, FDR-filtered) that are
    present in the expression matrix."""
    la = de[(de["FDR"] < FDR_CUT) & (de["logFC"] > LOGFC_CUT)]
    ra = de[(de["FDR"] < FDR_CUT) & (de["logFC"] < -LOGFC_CUT)]
    la = la.sort_values("moderated_t", ascending=False)
    ra = ra.sort_values("moderated_t", ascending=True)
    la = la[la["ensembl_id"].isin(expr_index)]
    ra = ra[ra["ensembl_id"].isin(expr_index)]
    la = la.head(TOP_EACH)
    ra = ra.head(TOP_EACH)
    top = pd.concat([la, ra], ignore_index=True)
    top = top.sort_values("logFC", ascending=False).reset_index(drop=True)
    return top


def row_zscore(mat: np.ndarray) -> np.ndarray:
    m = mat.mean(axis=1, keepdims=True)
    s = mat.std(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return (mat - m) / s


def draw_heatmap(ax: plt.Axes, strip_ax: plt.Axes,
                 de: pd.DataFrame, expr: pd.DataFrame, meta: pd.DataFrame,
                 cbar_cax: plt.Axes) -> None:
    top = pick_top_keys(de, expr.index)
    if top.empty:
        ax.text(0.5, 0.5, "No genes passed selection",
                ha="center", va="center", transform=ax.transAxes)
        return

    meta = meta.copy()
    meta["chamber"] = pd.Categorical(meta["chamber"],
                                     categories=["RA", "LA"], ordered=True)
    meta = meta.sort_values(["donor", "chamber", "sample_id"]).reset_index(drop=True)

    sample_cols = meta["sample_id"].tolist()
    row_keys = top["ensembl_id"].tolist()
    hm = expr.loc[row_keys, sample_cols].to_numpy(dtype=float)
    hm = row_zscore(hm)

    lim = float(np.max(np.abs(hm)))
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=+lim)

    n_rows, n_cols = hm.shape
    for i in range(n_rows):
        for j in range(n_cols):
            rgba = ZCMAP(norm(hm[i, j]))
            ax.add_patch(Rectangle(
                (j - 0.5, i - 0.5), 1.0, 1.0,
                facecolor=rgba, edgecolor="none",
            ))

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_aspect("auto")

    # Row labels
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(top["label"].tolist(), fontsize=6)
    ax.tick_params(axis="y", length=0, pad=2)

    # Color row labels by direction.
    for tl, direction in zip(ax.get_yticklabels(), top["direction"]):
        tl.set_color(LA_COLOR if direction.startswith("LA") else RA_COLOR)

    # Column ticks: donor tags (once per pair).
    donor_tags = meta["donor"].tolist()
    chamber_tags = meta["chamber"].astype(str).tolist()
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([f"{c}\n{d}" for c, d in zip(chamber_tags, donor_tags)],
                       fontsize=5, rotation=0)
    ax.tick_params(axis="x", length=0, pad=2)

    # Donor separators between pairs.
    donor_breaks = [j + 0.5 for j in range(n_cols - 1)
                    if donor_tags[j] != donor_tags[j + 1]]
    for xb in donor_breaks:
        ax.axvline(xb, color="#404040", lw=0.5)

    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_title("b2  Top DEG heatmap", fontsize=8, pad=16)

    for j, ch in enumerate(chamber_tags):
        rgba = LA_COLOR if ch == "LA" else RA_COLOR
        strip_ax.add_patch(Rectangle(
            (j - 0.5, 0), 1.0, 1.0,
            facecolor=rgba, edgecolor="white", linewidth=0.5,
        ))
    strip_ax.set_xlim(-0.5, n_cols - 0.5)
    strip_ax.set_ylim(0, 1)
    strip_ax.set_aspect("auto")
    strip_ax.set_xticks([])
    strip_ax.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        strip_ax.spines[s].set_visible(False)
    # Keep strip unlabeled to reduce panel text load.

    sm = ScalarMappable(norm=norm, cmap=ZCMAP)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=cbar_cax, orientation="horizontal")
    cbar.set_label("row z-score", fontsize=6, labelpad=2)
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(width=0.4, length=2, pad=1, labelsize=6)


def export_volcano(de: pd.DataFrame) -> None:
    """Export volcano as a standalone panel for downstream composite assembly."""
    fig_w = COL_SINGLE
    fig_h = COL_SINGLE * 0.95
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0.16, 0.15, 0.80, 0.76])
    draw_volcano(ax, de)
    OUT_VOLC.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_VOLC, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT_VOLC.relative_to(ROOT)}  (size {fig_w:.2f} x {fig_h:.2f} in)")


def export_heatmap(de: pd.DataFrame, expr: pd.DataFrame, meta: pd.DataFrame) -> None:
    """Export heatmap as a standalone panel for downstream composite assembly."""
    n_rows_hm = 2 * TOP_EACH
    # Keep row labels readable for publication assembly.
    row_in = 0.14
    hm_h_in = row_in * n_rows_hm
    top_pad = 0.55
    bot_pad = 1.05
    fig_h = hm_h_in + top_pad + bot_pad
    fig_w = COL_ONEHALF
    fig = plt.figure(figsize=(fig_w, fig_h))

    def fr_y(inches: float) -> float:
        return inches / fig_h

    # Allocate extra left margin for gene labels; this panel is often shrunk in montage.
    hm_left, hm_right = 0.28, 0.98
    hm_bot = fr_y(bot_pad)
    hm_top = hm_bot + fr_y(hm_h_in)
    hm_ax = fig.add_axes([hm_left, hm_bot, hm_right - hm_left, hm_top - hm_bot])

    strip_h = fr_y(0.14)
    strip_ax = fig.add_axes([hm_left, hm_top + fr_y(0.05), hm_right - hm_left, strip_h])

    cbar_w = 0.22
    cbar_cax = fig.add_axes([hm_left + (hm_right - hm_left - cbar_w) / 2, fr_y(0.30), cbar_w, fr_y(0.15)])

    draw_heatmap(hm_ax, strip_ax, de, expr, meta, cbar_cax)

    OUT_HMAP.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_HMAP, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT_HMAP.relative_to(ROOT)}  (size {fig_w:.2f} x {fig_h:.2f} in)")


def main() -> None:
    set_nature_style()
    de = load_de()
    expr, meta = load_expr_meta()
    export_volcano(de)
    export_heatmap(de, expr, meta)


if __name__ == "__main__":
    main()
