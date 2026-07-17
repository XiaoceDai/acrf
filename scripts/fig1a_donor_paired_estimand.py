#!/usr/bin/env python3
"""Final Figure 1A: donor-paired chamber estimand."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import FG_COLOR, LA_COLOR, LINE_COLOR, MUTED_COLOR, RA_COLOR, set_nature_style  # noqa: E402

META_IN = ROOT / "bulk" / "sample_meta.tsv"
SOURCE_OUT = ROOT / "results" / "fig1a_donor_paired_estimand_source.tsv"
SUMMARY_OUT = ROOT / "results" / "fig1a_donor_paired_estimand_summary.tsv"
FIG_BASE = ROOT / "figures" / "fig1a_donor_paired_estimand"

MM = 1 / 25.4
PANEL_WIDTH_MM = 68.2
PANEL_HEIGHT_MM = 105

set_nature_style()
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "mathtext.default": "regular",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.bbox": "standard",
    }
)


def build_provenance() -> tuple[pd.DataFrame, int]:
    meta = pd.read_csv(META_IN, sep="\t", dtype=str)
    required = {"sample_id", "chamber", "donor"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Missing metadata columns: {sorted(missing)}")

    meta = meta[list(required)].copy()
    meta["chamber"] = meta["chamber"].str.upper()
    unexpected = sorted(set(meta["chamber"].dropna()) - {"LA", "RA"})
    if unexpected:
        raise ValueError(f"Unexpected chamber labels: {unexpected}")

    counts = (
        meta.groupby(["donor", "chamber"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=["LA", "RA"], fill_value=0)
        .sort_index()
    )
    source = pd.DataFrame(
        {
            "donor_id": counts.index.astype(str),
            "has_LA": counts["LA"].to_numpy() > 0,
            "has_RA": counts["RA"].to_numpy() > 0,
            "complete_pair": (counts["LA"].to_numpy() > 0) & (counts["RA"].to_numpy() > 0),
        }
    )
    if source.empty:
        raise AssertionError("No healthy donors found in donor/chamber metadata")
    if not source["complete_pair"].all():
        bad = source.loc[~source["complete_pair"], "donor_id"].tolist()
        raise AssertionError(f"Analysis-level healthy donors missing a paired chamber profile: {bad}")

    n_pairs = int(source["complete_pair"].sum())
    SOURCE_OUT.parent.mkdir(parents=True, exist_ok=True)
    source.to_csv(SOURCE_OUT, sep="\t", index=False)
    pd.DataFrame(
        [
            {"metric": "analysis_level_healthy_donors", "value": len(source)},
            {"metric": "complete_LA_RA_donor_pairs", "value": n_pairs},
            {"metric": "all_analysis_donors_complete", "value": bool(source["complete_pair"].all())},
        ]
    ).to_csv(SUMMARY_OUT, sep="\t", index=False)
    return source, n_pairs


def draw_pair(ax: plt.Axes, y: float, donor_label: str) -> None:
    ax.text(0.12, y, donor_label, fontsize=7.1, color=FG_COLOR, ha="left", va="center")
    ax.scatter(0.42, y, s=80, marker="s", color=LA_COLOR, edgecolor="white", linewidth=0.45, zorder=3)
    ax.scatter(0.74, y, s=80, marker="s", color=RA_COLOR, edgecolor="white", linewidth=0.45, zorder=3)
    ax.plot([0.45, 0.71], [y, y], color=MUTED_COLOR, lw=0.75, zorder=2)
    ax.text(0.42, y, "LA", fontsize=6.5, color="white", fontweight="bold", ha="center", va="center", zorder=4)
    ax.text(0.74, y, "RA", fontsize=6.5, color="white", fontweight="bold", ha="center", va="center", zorder=4)


def render(source: pd.DataFrame, n_pairs: int) -> None:
    fig, ax = plt.subplots(figsize=(PANEL_WIDTH_MM * MM, PANEL_HEIGHT_MM * MM))
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.03, top=0.97)
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.04, 0.975, "A", fontsize=11.0, fontweight="bold", color=FG_COLOR, ha="left", va="top")
    ax.text(
        0.125,
        0.975,
        "Donor-paired chamber estimand",
        fontsize=8.5,
        fontweight="bold",
        color=FG_COLOR,
        ha="left",
        va="top",
    )
    ax.plot([0.04, 0.96], [0.925, 0.925], color=LINE_COLOR, lw=0.65)

    ax.text(0.58, 0.885, "Healthy reference donors", fontsize=7.3, fontweight="bold", color=FG_COLOR, ha="center")
    ax.text(0.58, 0.847, "sample-level program activity", fontsize=6.8, color=MUTED_COLOR, ha="center")
    ax.text(0.42, 0.812, "LA score", fontsize=6.7, color=LA_COLOR, fontweight="bold", ha="center")
    ax.text(0.74, 0.812, "RA score", fontsize=6.7, color=RA_COLOR, fontweight="bold", ha="center")

    representative = source["donor_id"].head(3).tolist()
    for i, _ in enumerate(representative, start=1):
        draw_pair(ax, 0.765 - (i - 1) * 0.072, f"Donor {i}")
    ax.text(0.58, 0.555, "⋮", fontsize=11, color=MUTED_COLOR, ha="center", va="center")
    ax.text(
        0.58,
        0.510,
        f"n = {n_pairs} complete LA–RA donor pairs",
        fontsize=7.2,
        fontweight="bold",
        color=FG_COLOR,
        ha="center",
        va="center",
    )

    ax.text(0.50, 0.468, "↓", fontsize=10.0, color=MUTED_COLOR, ha="center", va="center")
    ax.text(
        0.50,
        0.422,
        r"$\Delta_{d,p}=S^{LA}_{d,p}-S^{RA}_{d,p}$",
        fontsize=9.2,
        color=FG_COLOR,
        ha="center",
        va="center",
    )
    ax.text(
        0.50,
        0.380,
        r"$S$: sample-level program activity",
        fontsize=6.8,
        color=MUTED_COLOR,
        ha="center",
        va="center",
    )
    ax.text(
        0.50,
        0.323,
        "one donor × one program\n× one paired contrast",
        fontsize=7.1,
        fontweight="bold",
        color=FG_COLOR,
        ha="center",
        va="center",
        linespacing=1.18,
    )

    ax.text(
        0.08,
        0.245,
        "Within-donor subtraction reduces\n"
        "donor-specific baseline variation.",
        fontsize=7.1,
        color=FG_COLOR,
        ha="left",
        va="center",
        linespacing=1.25,
    )
    ax.plot([0.08, 0.92], [0.200, 0.200], color=LINE_COLOR, lw=0.6)
    ax.text(
        0.08,
        0.165,
        "Independent-sample pooling can mix\n"
        "donor and chamber effects.",
        fontsize=6.5,
        color=MUTED_COLOR,
        ha="left",
        va="center",
        linespacing=1.2,
    )
    ax.text(
        0.50,
        0.078,
        "Direction alignment and reference scaling → panel B",
        fontsize=6.5,
        color=MUTED_COLOR,
        ha="center",
        va="center",
    )

    FIG_BASE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_BASE.with_suffix(".pdf"), bbox_inches=None, facecolor="white")
    fig.savefig(FIG_BASE.with_suffix(".svg"), bbox_inches=None, facecolor="white")
    fig.savefig(FIG_BASE.with_suffix(".png"), dpi=600, bbox_inches=None, facecolor="white")
    plt.close(fig)


def main() -> None:
    source, n_pairs = build_provenance()
    render(source, n_pairs)
    print(f"donor/chamber source: {META_IN.relative_to(ROOT)}")
    print(f"complete paired donors: {n_pairs}")
    print(f"wrote {SOURCE_OUT.relative_to(ROOT)}")
    print(f"wrote {SUMMARY_OUT.relative_to(ROOT)}")
    print(f"wrote {FIG_BASE.relative_to(ROOT)}.[pdf|svg|png]")


if __name__ == "__main__":
    main()
