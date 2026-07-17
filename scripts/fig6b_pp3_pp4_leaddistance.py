#!/usr/bin/env python3
"""
Figure 6B — Aligned colocalisation posterior and lead-distance diagnostics.

Left: paired PP3 / PP4 points for each anchor on a shared 0-1 scale.
Right: lead GWAS-lead eQTL distance derived from absolute delta_bp / 1000.
KCNN3 is displayed as N/A by panel specification.

FROZEN INPUTS (see build_figure6_mr_v3_manuscript_tables.py):
  results/figure6_mr_anchors_v3_pipeline_wide.tsv
  results/fig6a_v3_locus_zoom_summary.tsv

Filenames containing "mr" denote coloc-based diagnostics; no Mendelian
randomization is performed.

Outputs:
  figures/fig6b_pp3_pp4_leaddistance.pdf
  figures/fig6b_pp3_pp4_leaddistance.svg
  figures/fig6b_pp3_pp4_leaddistance.png  (600 dpi)
  results/fig6b_pp3_pp4_leaddistance.tsv
  results/fig6b_pp3_pp4_leaddistance_qc.txt
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (  # noqa: E402
    COL_DOUBLE,
    FG_COLOR,
    HIGHLIGHT_COLOR,
    LINE_COLOR,
    MUTED_COLOR,
    set_nature_style,
    strip_spines,
)
from _figure_export import save_figure_bundle  # noqa: E402

PIPELINE_WIDE = ROOT / "results" / "figure6_mr_anchors_v3_pipeline_wide.tsv"
LOCUS_SUMMARY = ROOT / "results" / "fig6a_v3_locus_zoom_summary.tsv"

OUT_PDF = ROOT / "figures" / "fig6b_pp3_pp4_leaddistance.pdf"
OUT_SVG = ROOT / "figures" / "fig6b_pp3_pp4_leaddistance.svg"
OUT_TSV = ROOT / "results" / "fig6b_pp3_pp4_leaddistance.tsv"
OUT_QC = ROOT / "results" / "fig6b_pp3_pp4_leaddistance_qc.txt"

ANCHOR_ORDER = [
    "PRELP", "OMD", "IDH2", "NDUFA8", "BMP10", "SMAD9", "SMAD6",
    "NTRK2", "NTF3", "FRS2", "PITX2", "KCNN3",
]
CALIBRATORS = {"PITX2", "KCNN3"}
PP4_SUPPORT = 0.70
PP3_COLOR = "#526D82"
PP4_COLOR = HIGHLIGHT_COLOR
KEY_ANCHORS = {"SMAD9", "PITX2", "KCNN3"}
KCNN3_DISTANCE_NOTE = "N/A by panel specification"


def build_source_table() -> pd.DataFrame:
    wide = pd.read_csv(PIPELINE_WIDE, sep="\t").set_index("anchor")
    loc = pd.read_csv(LOCUS_SUMMARY, sep="\t").set_index("anchor")
    missing_wide = [g for g in ANCHOR_ORDER if g not in wide.index]
    if missing_wide:
        raise ValueError(f"Anchors missing from frozen coloc table: {missing_wide}")
    rows = []
    for g in ANCHOR_ORDER:
        r = wide.loc[g]
        raw_delta = (
            float(loc.loc[g, "delta_bp"])
            if g in loc.index and pd.notna(loc.loc[g, "delta_bp"])
            else np.nan
        )
        computed_kb = abs(raw_delta) / 1000.0 if pd.notna(raw_delta) else np.nan
        display_kb = np.nan if g == "KCNN3" else computed_kb
        rows.append({
            "gene": g,
            "anchor_type": "calibration" if g in CALIBRATORS else "framework",
            "coloc_PP3": float(r["coloc_PP3"]),
            "coloc_PP4": float(r["coloc_PP4"]),
            "PP4_supported": bool(float(r["coloc_PP4"]) >= PP4_SUPPORT),
            "PP3_gt_PP4": bool(float(r["coloc_PP3"]) > float(r["coloc_PP4"])),
            "source_delta_bp": raw_delta,
            "computed_distance_kb": (
                round(computed_kb, 1) if pd.notna(computed_kb) else np.nan
            ),
            "display_distance_kb": (
                round(display_kb, 1) if pd.notna(display_kb) else np.nan
            ),
            "distance_display_status": (
                KCNN3_DISTANCE_NOTE
                if g == "KCNN3"
                else ("observed" if pd.notna(display_kb) else "N/A in source")
            ),
        })
    return pd.DataFrame(rows)


def _computed_callout(df: pd.DataFrame) -> str:
    n_supported = int(df["PP4_supported"].sum())
    n_total = len(df)
    n_mismatch = int(df["PP3_gt_PP4"].sum())
    pattern = (
        "dominant pattern consistent with H3 signal mismatch"
        if n_mismatch > n_total / 2
        else "no dominant H3 signal-mismatch pattern"
    )
    return f"{n_supported}/{n_total} PP4-supported; {pattern}"


def _assert_source(df: pd.DataFrame) -> None:
    assert list(df["gene"]) == ANCHOR_ORDER
    assert len(df) == 12
    kcnn3 = df.loc[df["gene"] == "KCNN3"].iloc[0]
    assert pd.isna(kcnn3["display_distance_kb"])
    assert kcnn3["distance_display_status"] == KCNN3_DISTANCE_NOTE
    observed = df["display_distance_kb"].dropna()
    assert len(observed) == 11
    assert np.isfinite(observed).all()


def render(df: pd.DataFrame, callout: str):
    set_nature_style()
    fig = plt.figure(figsize=(COL_DOUBLE, COL_DOUBLE * 0.75), facecolor="white")
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.65, 1.0],
        left=0.13, right=0.98, top=0.78, bottom=0.14, wspace=0.16,
    )
    ax_pp = fig.add_subplot(gs[0, 0])
    ax_dist = fig.add_subplot(gs[0, 1], sharey=ax_pp)

    y = np.arange(len(df))
    key_mask = df["gene"].isin(KEY_ANCHORS)

    # Left: paired posterior points and connectors.
    for i, r in df.iterrows():
        ax_pp.plot(
            [r["coloc_PP4"], r["coloc_PP3"]], [i, i],
            color=LINE_COLOR, lw=1.0, zorder=1,
        )
    ax_pp.scatter(
        df["coloc_PP3"], y, marker="o", s=42,
        facecolor=PP3_COLOR, edgecolor="#334A5C", linewidth=0.6,
        label="PP3", zorder=3,
    )
    ax_pp.scatter(
        df["coloc_PP4"], y, marker="D", s=34,
        facecolor=PP4_COLOR, edgecolor="#785B2C", linewidth=0.6,
        label="PP4", zorder=4,
    )
    ax_pp.axvline(
        PP4_SUPPORT, color="#666666", lw=0.9, ls="--", zorder=0,
    )
    ax_pp.text(
        PP4_SUPPORT + 0.012, -0.72, "PP4 support threshold",
        fontsize=6.5, color=MUTED_COLOR, ha="left", va="center",
    )
    ax_pp.set_xlim(0, 1)
    ax_pp.set_xticks(np.linspace(0, 1, 6))
    ax_pp.set_xlabel("Colocalisation posterior probability", fontsize=8)
    ax_pp.set_yticks(y)
    labels = []
    for gene, is_key in zip(df["gene"], key_mask):
        labels.append(gene)
    ax_pp.set_yticklabels(labels, fontsize=7)
    for tick, is_key in zip(ax_pp.get_yticklabels(), key_mask):
        if is_key:
            tick.set_fontweight("bold")
            tick.set_color("#6B3A1F")
    ax_pp.invert_yaxis()
    ax_pp.legend(
        loc="lower right", frameon=False, ncol=2, fontsize=7,
        handletextpad=0.4, columnspacing=1.0,
    )
    ax_pp.set_title(
        "Paired PP3 and PP4 by anchor",
        loc="left", fontsize=8, fontweight="bold", pad=8,
    )
    strip_spines(ax_pp)

    # Right: aligned lead-distance lollipops.
    observed = df["display_distance_kb"].dropna()
    dmin = float(observed.min())
    dmax = float(observed.max())
    x_pad = max(50.0, 0.08 * (dmax - dmin))
    for i, r in df.iterrows():
        d = r["display_distance_kb"]
        if pd.notna(d):
            ax_dist.hlines(i, 0, d, color=LINE_COLOR, lw=1.3, zorder=1)
            ax_dist.scatter(
                d, i, s=38, marker="o",
                facecolor="#777777", edgecolor="#444444", linewidth=0.5,
                zorder=3,
            )
        else:
            ax_dist.text(
                0, i, "N/A", ha="left", va="center",
                fontsize=7, fontweight="bold", color=MUTED_COLOR,
            )
    ax_dist.set_xlim(0, dmax + x_pad)
    ax_dist.set_xlabel("Lead GWAS-lead eQTL distance (kb)", fontsize=8)
    ax_dist.set_title(
        "Lead-variant separation",
        loc="left", fontsize=8, fontweight="bold", pad=8,
    )
    ax_dist.tick_params(axis="y", left=False, labelleft=False)
    strip_spines(ax_dist)

    # Computed annotation and range note (title remains range-free).
    fig.text(
        0.02, 0.965,
        "Aligned colocalisation geometry and lead-variant separation",
        fontsize=9, fontweight="bold", ha="left", va="top", color=FG_COLOR,
    )
    fig.text(
        0.02, 0.925, callout,
        fontsize=7.5, fontweight="bold", ha="left", va="top", color=FG_COLOR,
    )
    fig.text(
        0.02, 0.892,
        f"Observed lead-distance range: {dmin:.1f}-{dmax:.1f} kb "
        "(KCNN3 displayed as N/A).",
        fontsize=7, ha="left", va="top", color=MUTED_COLOR,
    )
    fig.text(
        0.02, 0.045,
        "Rows follow Figure 6A order. Circle = PP3; diamond = PP4. "
        "Distances are abs(delta_bp)/1000 from the frozen locus table.",
        fontsize=6.5, ha="left", va="bottom", color=MUTED_COLOR,
        style="italic",
    )
    return fig


def main() -> None:
    df = build_source_table()
    _assert_source(df)
    callout = _computed_callout(df)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    observed = df["display_distance_kb"].dropna()
    dmin, dmax = float(observed.min()), float(observed.max())
    qc_lines = [
        "Figure 6B QC",
        f"anchors: {len(df)}",
        f"PP4_supported: {int(df['PP4_supported'].sum())}/{len(df)}",
        f"observed_distance_min_kb: {dmin:.1f}",
        f"observed_distance_max_kb: {dmax:.1f}",
        "KCNN3_display_distance_kb: NA",
        f"KCNN3_status: {KCNN3_DISTANCE_NOTE}",
        "KCNN3_NA_check: PASS",
    ]
    OUT_QC.write_text("\n".join(qc_lines) + "\n", encoding="utf-8")

    fig = render(df, callout)
    pdf_path, png_path = save_figure_bundle(
        fig, OUT_PDF, dpi=600, bbox_inches=None,
    )
    fig = render(df, callout)
    fig.savefig(
        OUT_SVG, format="svg", bbox_inches=None, pad_inches=0.04,
        facecolor="white", edgecolor="none",
    )
    plt.close(fig)

    print(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(df)} anchors)")
    print(f"wrote {OUT_QC.relative_to(ROOT)}")
    print(f"callout: {callout}")
    print(f"lead-distance range: {dmin:.1f}-{dmax:.1f} kb")
    print(f"wrote {pdf_path.relative_to(ROOT)}")
    if png_path:
        print(f"wrote {png_path.relative_to(ROOT)}")
    print(f"wrote {OUT_SVG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
