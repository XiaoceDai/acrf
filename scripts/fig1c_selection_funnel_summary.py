#!/usr/bin/env python3
"""Final Figure 1C: healthy-only selection provenance for the frozen panel."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import innovation_ncp_consensus as inc  # noqa: E402
from plot_style import (  # noqa: E402
    COL_DOUBLE,
    FG_COLOR,
    LA_COLOR,
    LINE_COLOR,
    MUTED_COLOR,
    RA_COLOR,
    set_nature_style,
)

LONG = ROOT / "results" / "fig1c_pathway_selection_long.tsv"
SCORED = ROOT / "results" / "innovation_ncp_consensus_long.tsv"
PANEL = ROOT / "results" / "innovation_ncp_v3_panel.tsv"
FAMILY_MAP = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
OUT_BASE = ROOT / "figures" / "fig1c_selection_provenance"
OUT_SOURCE = ROOT / "results" / "fig1c_selection_provenance_source.tsv"

set_nature_style()
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def build_source() -> tuple[pd.DataFrame, dict[str, int]]:
    long_df = pd.read_csv(LONG, sep="\t")
    scored = pd.read_csv(SCORED, sep="\t")
    panel = pd.read_csv(PANEL, sep="\t")
    fam = pd.read_csv(FAMILY_MAP, sep="\t")

    if len(long_df) != 80 or len(scored) != 80:
        raise AssertionError("Expected 80 candidate pathways")
    if len(panel) != 30:
        raise AssertionError("Expected 30 deployed programs")

    strict = (
        scored["direction_agree"].astype(bool)
        & (scored["stability_consistent_runs"].astype(float) >= inc.MIN_STAB)
        & (scored["consensus_score"].astype(float) >= inc.MIN_SCORE)
        & (scored["padj"].astype(float) < 0.05)
        & (scored["gsva_padj"].astype(float) < 0.05)
    )
    scored = scored.copy()
    scored["strict_pool"] = strict
    deployed = set(panel["pathway"].astype(str))
    scored["deployed"] = scored["pathway"].astype(str).isin(deployed)
    scored["deployed_tier"] = np.select(
        [scored["deployed"] & scored["strict_pool"], scored["deployed"] & ~scored["strict_pool"]],
        ["strict", "fallback"],
        default="not_deployed",
    )

    counts = {
        "n_candidates": int(len(scored)),
        "n_strict": int(scored["strict_pool"].sum()),
        "n_deployed": int(scored["deployed"].sum()),
        "n_deployed_strict": int((scored["deployed_tier"] == "strict").sum()),
        "n_deployed_fallback": int((scored["deployed_tier"] == "fallback").sum()),
        "n_la": int((panel["direction"] == "LA-biased").sum()),
        "n_ra": int((panel["direction"] == "RA-biased").sum()),
        "n_families": int(fam["family_id"].nunique()),
    }
    expected = {
        "n_candidates": 80,
        "n_strict": 38,
        "n_deployed": 30,
        "n_deployed_strict": 25,
        "n_deployed_fallback": 5,
        "n_la": 22,
        "n_ra": 8,
        "n_families": 16,
    }
    if counts != expected:
        raise AssertionError(f"Selection provenance changed: {counts}")

    source_cols = [
        "pathway",
        "direction",
        "NES",
        "padj",
        "gsva_logFC",
        "gsva_padj",
        "stability_consistent_runs",
        "consensus_score",
        "direction_agree",
        "strict_pool",
        "deployed",
        "deployed_tier",
    ]
    OUT_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    scored[source_cols].to_csv(OUT_SOURCE, sep="\t", index=False)
    return scored, counts


def _node(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    subtitle: str,
    *,
    face: str = "#F7F7F7",
    edge: str = LINE_COLOR,
) -> None:
    ax.add_patch(
        patches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.008,rounding_size=0.018",
            facecolor=face,
            edgecolor=edge,
            linewidth=0.9,
            transform=ax.transAxes,
        )
    )
    ax.text(x + w / 2, y + h * 0.64, title, fontsize=9.0, fontweight="bold", ha="center", va="center", transform=ax.transAxes)
    ax.text(x + w / 2, y + h * 0.25, subtitle, fontsize=6.7, color=MUTED_COLOR, ha="center", va="center", linespacing=1.15, transform=ax.transAxes)


def draw_provenance(ax: plt.Axes, counts: dict[str, int]) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Healthy-only selection provenance", loc="left", fontsize=9.2, fontweight="bold", pad=7)

    _node(ax, 0.05, 0.75, 0.90, 0.14, "80 candidates", "healthy paired GSEA shortlist")
    _node(ax, 0.05, 0.48, 0.57, 0.14, "38 strict pool", "dual evidence\n+ stability + score")
    _node(ax, 0.67, 0.48, 0.28, 0.14, "5 fallback", "selected RA\narm", face="#FAFAFA")
    _node(
        ax,
        0.05,
        0.12,
        0.90,
        0.19,
        "30 frozen programs",
        "25 strict + 5 fallback\n diversity and direction lock",
        face="#F3EFE9",
        edge="#C8A076",
    )
    ax.annotate("", xy=(0.335, 0.63), xytext=(0.50, 0.75), xycoords=ax.transAxes, arrowprops={"arrowstyle": "-|>", "color": MUTED_COLOR, "lw": 1.0})
    ax.annotate("", xy=(0.81, 0.63), xytext=(0.50, 0.75), xycoords=ax.transAxes, arrowprops={"arrowstyle": "-|>", "color": MUTED_COLOR, "lw": 1.0})
    ax.annotate("", xy=(0.42, 0.32), xytext=(0.335, 0.48), xycoords=ax.transAxes, arrowprops={"arrowstyle": "-|>", "color": MUTED_COLOR, "lw": 1.0})
    ax.annotate("", xy=(0.70, 0.32), xytext=(0.81, 0.48), xycoords=ax.transAxes, arrowprops={"arrowstyle": "-|>", "color": MUTED_COLOR, "lw": 0.9})
    ax.text(
        0.05,
        0.025,
        "AF cohorts were not used for candidate scoring,\npruning, balancing, or reselection.",
        fontsize=6.8,
        color=FG_COLOR,
        fontweight="bold",
        ha="left",
        va="bottom",
        transform=ax.transAxes,
    )


def draw_scatter(ax: plt.Axes, df: pd.DataFrame) -> None:
    x = df["NES"].to_numpy(float)
    y = df["gsva_logFC"].to_numpy(float)
    strict = df["strict_pool"].to_numpy(bool)
    deployed = df["deployed"].to_numpy(bool)
    dirs = df["direction"].astype(str).to_numpy()

    ax.set_facecolor("#FAFAFA")
    ax.axhline(0, color=LINE_COLOR, lw=0.8, ls="--", zorder=0)
    ax.axvline(0, color=LINE_COLOR, lw=0.8, ls="--", zorder=0)
    ax.scatter(x, y, s=18, color="#C9C9C9", edgecolor="white", linewidth=0.25, zorder=1)
    for direction, color in (("LA-biased", LA_COLOR), ("RA-biased", RA_COLOR)):
        mask = strict & (dirs == direction)
        ax.scatter(x[mask], y[mask], s=25, color=color, edgecolor="white", linewidth=0.3, zorder=2)
    ax.scatter(x[deployed], y[deployed], s=48, facecolor="none", edgecolor=FG_COLOR, linewidth=0.75, zorder=3)

    lim_x = max(0.8, float(np.nanmax(np.abs(x))) * 1.08)
    lim_y = max(0.4, float(np.nanmax(np.abs(y))) * 1.08)
    ax.set_xlim(-lim_x, lim_x)
    ax.set_ylim(-lim_y, lim_y)
    ax.set_xlabel("signed GSEA NES", fontsize=7.6)
    ax.set_ylabel("paired GSVA logFC (LA − RA)", fontsize=7.6)
    ax.set_title("Healthy dual-evidence landscape", loc="left", fontsize=9.2, fontweight="bold", pad=7)
    ax.tick_params(labelsize=6.8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    handles = [
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="#C9C9C9", markeredgecolor="white", markersize=4.5, label="candidate"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor=LA_COLOR, markeredgecolor="white", markersize=4.5, label="strict LA"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor=RA_COLOR, markeredgecolor="white", markersize=4.5, label="strict RA"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="none", markeredgecolor=FG_COLOR, markersize=5.0, label="deployed"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=6.7, ncol=2, handletextpad=0.3, columnspacing=0.7)


def draw_lock(ax: plt.Axes, counts: dict[str, int]) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Locked deployed object", loc="left", fontsize=9.2, fontweight="bold", pad=7)
    items = [
        ("Programs", f"{counts['n_deployed']}"),
        ("Direction", f"{counts['n_la']} LA  /  {counts['n_ra']} RA"),
        ("Families", f"{counts['n_families']}"),
    ]
    y0 = 0.78
    for i, (label, value) in enumerate(items):
        yy = y0 - i * 0.22
        ax.text(0.06, yy + 0.035, label, fontsize=7.3, color=MUTED_COLOR, fontweight="bold", ha="left", va="center")
        ax.text(0.94, yy - 0.030, value, fontsize=11.0 if i != 1 else 8.5, color=FG_COLOR, fontweight="bold", ha="right", va="center")
        ax.plot([0.06, 0.94], [yy - 0.105, yy - 0.105], color=LINE_COLOR, lw=0.7)
    ax.text(
        0.06,
        0.12,
        "Membership and direction\nwere frozen before AF analyses.",
        fontsize=7.6,
        color=FG_COLOR,
        fontweight="bold",
        ha="left",
        va="bottom",
        linespacing=1.3,
    )


def render(df: pd.DataFrame, counts: dict[str, int]) -> None:
    fig = plt.figure(figsize=(COL_DOUBLE * 1.25, COL_DOUBLE * 0.67))
    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[1.12, 1.85, 0.78],
        left=0.045,
        right=0.985,
        bottom=0.17,
        top=0.82,
        wspace=0.20,
    )
    ax_flow = fig.add_subplot(gs[0, 0])
    ax_scatter = fig.add_subplot(gs[0, 1])
    ax_lock = fig.add_subplot(gs[0, 2])

    fig.text(0.045, 0.965, "c", fontsize=13, fontweight="bold", ha="left", va="top")
    fig.text(
        0.070,
        0.965,
        "Selection provenance: healthy candidates to the frozen deployed panel",
        fontsize=10.8,
        fontweight="bold",
        ha="left",
        va="top",
    )
    fig.text(
        0.070,
        0.910,
        "Selection used paired healthy atrial evidence only; downstream AF cohorts were evaluation data.",
        fontsize=7.8,
        color=MUTED_COLOR,
        ha="left",
        va="top",
    )

    draw_provenance(ax_flow, counts)
    draw_scatter(ax_scatter, df)
    draw_lock(ax_lock, counts)

    OUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BASE.with_suffix(".pdf"), bbox_inches=None, facecolor="white")
    fig.savefig(OUT_BASE.with_suffix(".svg"), bbox_inches=None, facecolor="white")
    fig.savefig(OUT_BASE.with_suffix(".png"), dpi=600, bbox_inches=None, facecolor="white")
    plt.close(fig)


def main() -> None:
    df, counts = build_source()
    render(df, counts)
    print(f"wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"summary: {counts}")


if __name__ == "__main__":
    main()
