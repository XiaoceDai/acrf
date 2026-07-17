#!/usr/bin/env python3
"""Supplementary Figure S4 — Per-cohort PCA mini-grid."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (
    COL_DOUBLE,
    LA_COLOR,
    RA_COLOR,
    panel_label,
    set_nature_style,
    strip_spines,
)

set_nature_style()

FIG_OUT  = ROOT / "figures" / "figS4_per_cohort_pca.pdf"
TBL_OUT  = ROOT / "results" / "figS4_per_cohort_pca.tsv"

AF_COLOR = "#B2182B"
SR_COLOR = "#2166AC"
HEALTHY_LA = LA_COLOR
HEALTHY_RA = RA_COLOR


# PCA helper.
def top_variance_pca(expr: pd.DataFrame, n_features: int = 2_000,
                     log_transform: bool = False) -> tuple[pd.DataFrame, np.ndarray]:
    """Return DataFrame with PC1, PC2 per sample (columns of expr)."""
    from sklearn.decomposition import PCA

    x = expr.copy()
    x = x.replace([np.inf, -np.inf], np.nan).dropna(how="any", axis=0)
    if log_transform:
        x = np.log2(x.clip(lower=0) + 1.0)
    # drop zero-variance features then take top variance
    rv = x.var(axis=1)
    rv = rv[rv > 1e-8].sort_values(ascending=False)
    keep = rv.head(min(n_features, len(rv))).index
    x = x.loc[keep]

    M = x.values.T  # samples × features
    M = M - M.mean(axis=0, keepdims=True)
    pca = PCA(n_components=2, svd_solver="randomized", random_state=0)
    pcs = pca.fit_transform(M)
    var_explained = pca.explained_variance_ratio_
    return (pd.DataFrame(pcs, index=x.columns, columns=["PC1", "PC2"]),
            np.asarray(var_explained))


# Cohort plotters.
def plot_gse112339(ax: plt.Axes) -> tuple[str, list[dict]]:
    """Use precomputed PCA in qc_pca.tsv."""
    qc = pd.read_csv(ROOT / "bulk" / "Project_P01" / "qc" / "qc_pca.tsv",
                     sep="\t")
    rows = []
    for _, r in qc.iterrows():
        col = HEALTHY_LA if r["chamber"] == "LA" else HEALTHY_RA
        ax.scatter(r["PC1"], r["PC2"], color=col, edgecolor="white",
                   linewidth=0.3, s=12, zorder=3)
        rows.append({
            "cohort": "GSE112339", "sample_id": r["sample_id"],
            "PC1": r["PC1"], "PC2": r["PC2"],
            "label": r["chamber"], "category": "healthy_chamber",
        })
    # connect paired LA<->RA samples
    for d, sub in qc.groupby("donor"):
        if len(sub) == 2:
            ax.plot(sub["PC1"], sub["PC2"], color="#888888",
                    linewidth=0.3, alpha=0.6, zorder=1)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("GSE112339 (n = 16; healthy LA / RA)", pad=2)
    return "PC1 / PC2 from qc_pca.tsv", rows


def plot_paired_chamber_cohort(ax: plt.Axes, cohort: str, expr_path: Path,
                               md_path: Path, log_transform: bool
                               ) -> tuple[str, list[dict]]:
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    md   = pd.read_csv(md_path, sep="\t")
    common = [c for c in expr.columns if c in set(md["gsm"])]
    expr = expr[common]
    md   = md.set_index("gsm").loc[common]

    pcs, var = top_variance_pca(expr, n_features=2_000,
                                log_transform=log_transform)
    rows = []
    for sid, row in pcs.iterrows():
        chamber = md.loc[sid, "chamber"]
        rhythm  = md.loc[sid, "rhythm"]
        col = AF_COLOR if rhythm == "AF" else SR_COLOR
        marker = "o" if chamber == "LA" else "^"
        ax.scatter(row["PC1"], row["PC2"], color=col, marker=marker,
                   edgecolor="white", linewidth=0.3, s=14, zorder=3)
        rows.append({
            "cohort": cohort, "sample_id": sid,
            "PC1": float(row["PC1"]), "PC2": float(row["PC2"]),
            "label": f"{rhythm}-{chamber}", "category": "AF_vs_SR",
        })
    ax.set_xlabel(f"PC1 ({var[0]*100:.0f} %)")
    ax.set_ylabel(f"PC2 ({var[1]*100:.0f} %)")
    n_la = int((md["chamber"] == "LA").sum())
    n_ra = int((md["chamber"] == "RA").sum())
    ax.set_title(f"{cohort} (LA = {n_la}, RA = {n_ra})", pad=2)
    return f"top-2000 var PCA, {'log2(x+1)' if log_transform else 'as-is'}", rows


def plot_gse41177(ax: plt.Axes) -> tuple[str, list[dict]]:
    expr_path = ROOT / "Project_P4" / "GSE41177" / "expr_symbol_LAA.tsv"
    md_path   = ROOT / "Project_P4" / "GSE41177" / "sample_metadata.tsv"
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    md   = pd.read_csv(md_path, sep="\t")
    common = [c for c in expr.columns if c in set(md["gsm"])]
    expr = expr[common]
    md   = md.set_index("gsm").loc[common]

    pcs, var = top_variance_pca(expr, n_features=2_000, log_transform=False)
    rows = []
    for sid, row in pcs.iterrows():
        rhythm = md.loc[sid, "rhythm"]
        col = AF_COLOR if rhythm == "AF" else SR_COLOR
        ax.scatter(row["PC1"], row["PC2"], color=col, marker="o",
                   edgecolor="white", linewidth=0.3, s=14, zorder=3)
        rows.append({
            "cohort": "GSE41177", "sample_id": sid,
            "PC1": float(row["PC1"]), "PC2": float(row["PC2"]),
            "label": rhythm, "category": "AF_vs_SR_LAA",
        })
    ax.set_xlabel(f"PC1 ({var[0]*100:.0f} %)")
    ax.set_ylabel(f"PC2 ({var[1]*100:.0f} %)")
    n_af = int((md["rhythm"] == "AF").sum())
    n_sr = int((md["rhythm"] == "SR").sum())
    ax.set_title(f"GSE41177 LAA (AF = {n_af}, SR = {n_sr})", pad=2)
    return "top-2000 var PCA on LAA samples only", rows


def main() -> None:
    fig = plt.figure(figsize=(COL_DOUBLE, 0.5 * COL_DOUBLE))
    gs = fig.add_gridspec(2, 3, hspace=0.65, wspace=0.55)

    all_rows: list[dict] = []
    panels = [
        ("a", "GSE112339", plot_gse112339, None),
        ("b", "GSE128188", "paired", True),
        ("c", "GSE79768",  "paired", False),
        ("d", "GSE115574", "paired", False),
        ("e", "GSE41177",  plot_gse41177, None),
    ]

    coords = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
    for (lab, name, fn_or_kind, log_flag), (r, c) in zip(panels, coords):
        ax = fig.add_subplot(gs[r, c])
        if callable(fn_or_kind):
            note, rows = fn_or_kind(ax)
        else:
            cohort_dir = "Project_P1"
            note, rows = plot_paired_chamber_cohort(
                ax, name,
                ROOT / cohort_dir / name / "expr_symbol.tsv",
                ROOT / cohort_dir / name / "sample_metadata.tsv",
                log_transform=log_flag,
            )
        strip_spines(ax)
        panel_label(ax, lab, x=-0.18, y=1.06)
        all_rows.extend(rows)

    legend_ax = fig.add_subplot(gs[1, 2])
    legend_ax.axis("off")
    handles = [
        plt.Line2D([], [], marker="o", color="white",
                   markerfacecolor=HEALTHY_LA, markeredgecolor="white",
                   markersize=6, label="healthy LA (Fig 1)"),
        plt.Line2D([], [], marker="o", color="white",
                   markerfacecolor=HEALTHY_RA, markeredgecolor="white",
                   markersize=6, label="healthy RA (Fig 1)"),
        plt.Line2D([], [], marker="o", color="white",
                   markerfacecolor=SR_COLOR, markeredgecolor="white",
                   markersize=6, label="SR — LA"),
        plt.Line2D([], [], marker="^", color="white",
                   markerfacecolor=SR_COLOR, markeredgecolor="white",
                   markersize=6, label="SR — RA"),
        plt.Line2D([], [], marker="o", color="white",
                   markerfacecolor=AF_COLOR, markeredgecolor="white",
                   markersize=6, label="AF — LA"),
        plt.Line2D([], [], marker="^", color="white",
                   markerfacecolor=AF_COLOR, markeredgecolor="white",
                   markersize=6, label="AF — RA"),
    ]
    legend_ax.legend(handles=handles, loc="center", frameon=False,
                     ncols=1, fontsize=6, handlelength=1.5)
    legend_ax.text(0.5, -0.1,
                   "Each panel : top-variance 2,000 features, sample-centred PCA",
                   transform=legend_ax.transAxes,
                   ha="center", va="top", fontsize=5, color="#444444")

    fig.suptitle("Per-cohort principal-component overview",
                 y=1.00, fontsize=8)

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")

    pd.DataFrame(all_rows).to_csv(TBL_OUT, sep="\t", index=False)
    print(f"wrote {TBL_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
