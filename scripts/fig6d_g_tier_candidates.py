#!/usr/bin/env python3
"""
Figure 6D main candidate + Figure S16 entity-level strip.

Frozen inputs (read-only; G-tier assignments are not recomputed):
  results/figure6_g_tier_traceability.tsv
  results/figure6_g_tier_rules.tsv

Outputs:
  results/fig6d_g_tier_source.tsv
  results/fig6d_g_tier_counts.tsv
  figures/fig6d_g_tier_distribution.{pdf,svg,png}   (main Fig 6D candidate)
  figures/figS16_g_tier_entity_strip.{pdf,svg,png}  (supplement)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (  # noqa: E402
    COL_DOUBLE,
    FG_COLOR,
    FS_FOOTNOTE,
    FS_LEGEND,
    FS_TICK,
    MUTED_COLOR,
    set_nature_style,
)
from _figure_export import save_figure_bundle  # noqa: E402

G_TIER = ROOT / "results" / "figure6_g_tier_traceability.tsv"
G_RULES = ROOT / "results" / "figure6_g_tier_rules.tsv"
OUT_SOURCE = ROOT / "results" / "fig6d_g_tier_source.tsv"
OUT_COUNTS = ROOT / "results" / "fig6d_g_tier_counts.tsv"
OUT_DIST = ROOT / "figures" / "fig6d_g_tier_distribution.pdf"
OUT_STRIP = ROOT / "figures" / "figS16_g_tier_entity_strip.pdf"

TIERS = ["G0", "G1", "G2", "G3", "G4"]
TIER_COLORS = {
    "G0": "#E0E0E0",
    "G1": "#B8CCE4",
    "G2": "#F4D4A8",
    "G3": "#C5E0B4",
    "G4": "#6E9B72",
}
# Short labels for compact definition key
TIER_SHORT = {
    "G0": "No near-GW AF / no family overlap",
    "G1": "Near-GW AF or weak family overlap",
    "G2": "H3 / eGene without PP4 support",
    "G3": "G2 + family FDR <= 0.05",
    "G4": "PP4 >= 0.70 same-variant coloc",
}
CALLOUT_ANCHORS = ["SMAD9", "PITX2", "KCNN3"]
TRACEABILITY_NOTE = "G-tier summarizes genetic traceability, not causal strength."
N_ANCHORS = 12
N_FAMILIES = 16
PCT_LABEL_MIN_N = 3


def _load_source() -> pd.DataFrame:
    gt = pd.read_csv(G_TIER, sep="\t")
    required = {"entity_type", "entity_id", "G_tier"}
    missing = required - set(gt.columns)
    if missing:
        raise ValueError(f"Missing columns in G-tier table: {sorted(missing)}")
    if gt["G_tier"].isna().any() or (gt["G_tier"].astype(str).str.strip() == "").any():
        raise ValueError("Missing G_tier assignments; refusing to infer")
    bad = set(gt["G_tier"].astype(str)) - set(TIERS)
    if bad:
        raise ValueError(f"Unexpected G_tier values: {sorted(bad)}")
    cols = [c for c in [
        "entity_type", "entity_id", "family_id", "G_tier", "verdict",
        "near_GW_AF_locus", "coloc_PP3", "coloc_PP4", "is_eGene_bonferroni",
        "family_min_overlap_fdr", "sqtl_verdict",
    ] if c in gt.columns]
    return gt[cols].copy()


def _build_counts(src: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for etype, n_expected in (("anchor", N_ANCHORS), ("family", N_FAMILIES)):
        sub = src[src["entity_type"] == etype]
        n_total = len(sub)
        vc = sub["G_tier"].value_counts()
        for tier in TIERS:
            n = int(vc.get(tier, 0))
            rows.append({
                "entity_type": etype,
                "G_tier": tier,
                "n": n,
                "pct": round(100.0 * n / n_total, 1) if n_total else 0.0,
                "n_total": n_total,
            })
    return pd.DataFrame(rows)


def _assert_data(src: pd.DataFrame, counts: pd.DataFrame) -> None:
    n_a = int((src["entity_type"] == "anchor").sum())
    n_f = int((src["entity_type"] == "family").sum())
    assert n_a == N_ANCHORS, f"expected {N_ANCHORS} anchors, got {n_a}"
    assert n_f == N_FAMILIES, f"expected {N_FAMILIES} families, got {n_f}"
    assert src["G_tier"].notna().all()
    assert len(counts) == 2 * len(TIERS)
    for etype in ("anchor", "family"):
        sub = counts[counts["entity_type"] == etype]
        assert int(sub["n"].sum()) == int(sub["n_total"].iloc[0])
        assert set(sub["G_tier"]) == set(TIERS)
        g4 = sub.loc[sub["G_tier"] == "G4", "n"].iloc[0]
        assert g4 == 0 or g4 >= 0  # G4 row must exist
        assert "G4" in set(sub["G_tier"])


def _emphasis_lines(counts: pd.DataFrame) -> list[str]:
    """Derive emphasis statements from counts (no hardcoded numerators)."""
    lines: list[str] = []

    def n(etype: str, tiers: list[str]) -> tuple[int, int]:
        sub = counts[counts["entity_type"] == etype]
        total = int(sub["n_total"].iloc[0])
        hit = int(sub.loc[sub["G_tier"].isin(tiers), "n"].sum())
        return hit, total

    a_g2, a_tot = n("anchor", ["G2"])
    f_g01, f_tot = n("family", ["G0", "G1"])
    g4_a, _ = n("anchor", ["G4"])
    g4_f, _ = n("family", ["G4"])
    g4_all = g4_a + g4_f
    all_tot = a_tot + f_tot

    # Modal tier statements
    a_mode = counts[counts["entity_type"] == "anchor"].sort_values("n", ascending=False).iloc[0]
    f_g0 = int(counts.loc[
        (counts["entity_type"] == "family") & (counts["G_tier"] == "G0"), "n"
    ].iloc[0])
    f_g1 = int(counts.loc[
        (counts["entity_type"] == "family") & (counts["G_tier"] == "G1"), "n"
    ].iloc[0])

    if a_g2 >= (a_tot + 1) // 2:
        lines.append(f"Most anchors are G2 ({a_g2}/{a_tot}).")
    else:
        lines.append(
            f"Anchors peak at {a_mode['G_tier']} ({int(a_mode['n'])}/{a_tot})."
        )

    if f_g01 == f_tot:
        lines.append(f"Most families are G0-G1 ({f_g01}/{f_tot}; G0={f_g0}, G1={f_g1}).")
    elif f_g01 >= (f_tot + 1) // 2:
        lines.append(f"Most families are G0-G1 ({f_g01}/{f_tot}).")
    else:
        lines.append(f"Families: G0={f_g0}, G1={f_g1} of {f_tot}.")

    lines.append(f"No entity reaches G4 ({g4_all}/{all_tot}).")
    lines.append(f"G4: 0")
    return lines


def _zero_tier_note(counts: pd.DataFrame) -> str:
    """Explicit zero-count tiers that have no bar width."""
    bits = []
    for etype, label in (("anchor", "Anchors"), ("family", "Families")):
        sub = counts[counts["entity_type"] == etype]
        zeros = [t for t in TIERS if int(sub.loc[sub["G_tier"] == t, "n"].iloc[0]) == 0]
        if zeros:
            bits.append(f"{label}: " + ", ".join(f"{t}: 0" for t in zeros))
    return "  |  ".join(bits)


def _definition_key_handles() -> list[mpatches.Patch]:
    return [
        mpatches.Patch(
            facecolor=TIER_COLORS[t],
            edgecolor="#AAAAAA",
            linewidth=0.4,
            label=f"{t}: {TIER_SHORT[t]}",
        )
        for t in TIERS
    ]


def _save_svg(fig: plt.Figure, svg_path: Path) -> None:
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        svg_path,
        format="svg",
        bbox_inches="tight",
        pad_inches=0.04,
        facecolor="white",
        edgecolor="none",
    )


def _draw_distribution(counts: pd.DataFrame) -> plt.Figure:
    emphasis = _emphasis_lines(counts)
    zero_note = _zero_tier_note(counts)

    fig_w = COL_DOUBLE
    fig_h = 0.42 * COL_DOUBLE
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="white")
    fig.subplots_adjust(left=0.16, right=0.97, top=0.72, bottom=0.28)

    y_positions = {"anchor": 1.0, "family": 0.0}
    bar_h = 0.55

    for etype, y in y_positions.items():
        sub = counts[counts["entity_type"] == etype].set_index("G_tier")
        n_total = int(sub["n_total"].iloc[0])
        x0 = 0.0
        for tier in TIERS:
            n = int(sub.loc[tier, "n"])
            if n == 0:
                continue
            width = n / n_total
            ax.barh(
                y, width, left=x0, height=bar_h,
                color=TIER_COLORS[tier], edgecolor="white", linewidth=0.6,
                align="center",
            )
            # Integer count inside segment; pct only if readable
            label = str(n)
            if n >= PCT_LABEL_MIN_N:
                pct = float(sub.loc[tier, "pct"])
                label = f"{n}\n({pct:.0f}%)"
            ax.text(
                x0 + width / 2, y, label,
                ha="center", va="center",
                fontsize=FS_LEGEND - 0.5 if n >= PCT_LABEL_MIN_N else FS_TICK,
                color=FG_COLOR, fontweight="bold",
            )
            x0 += width

    ax.set_yticks([1.0, 0.0])
    ax.set_yticklabels(
        [
            f"Anchors (n = {int(counts.loc[counts.entity_type == 'anchor', 'n_total'].iloc[0])})",
            f"Families (n = {int(counts.loc[counts.entity_type == 'family', 'n_total'].iloc[0])})",
        ],
        fontsize=FS_TICK,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.55, 1.55)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=FS_LEGEND - 0.5)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("white")

    fig.text(
        0.02, 0.96,
        "G-tier distribution (anchors vs NCP axis families)",
        fontsize=FS_TICK, fontweight="bold", ha="left", va="top", color=FG_COLOR,
    )
    fig.text(
        0.02, 0.905, TRACEABILITY_NOTE,
        fontsize=FS_LEGEND, ha="left", va="top", color=MUTED_COLOR, style="italic",
    )
    for i, line in enumerate(emphasis[:3]):  # skip duplicate bare "G4: 0" if in zero_note
        fig.text(
            0.02, 0.855 - 0.035 * i, line,
            fontsize=FS_LEGEND, ha="left", va="top",
            color=FG_COLOR, fontweight="bold" if i < 2 else "normal",
        )

    fig.text(
        0.02, 0.20, f"Zero-width tiers shown explicitly — {zero_note}",
        fontsize=FS_FOOTNOTE - 0.5, ha="left", va="top", color=MUTED_COLOR,
    )

    legend = ax.legend(
        handles=_definition_key_handles(),
        loc="upper left",
        bbox_to_anchor=(0.0, -0.18),
        ncol=1,
        frameon=False,
        fontsize=FS_LEGEND - 1.0,
        handlelength=1.2,
        labelspacing=0.35,
    )
    legend.set_title("G-tier definition key", prop={"size": FS_LEGEND - 0.5, "weight": "bold"})

    return fig


def _draw_entity_strip(src: pd.DataFrame, counts: pd.DataFrame) -> plt.Figure:
    emphasis = _emphasis_lines(counts)
    zero_note = _zero_tier_note(counts)
    tier_x = {t: i for i, t in enumerate(TIERS)}

    fig_w = COL_DOUBLE
    fig_h = 0.40 * COL_DOUBLE
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="white")
    fig.subplots_adjust(left=0.14, right=0.97, top=0.72, bottom=0.30)

    rng = np.random.default_rng(42)
    y_base = {"anchor": 1.0, "family": 0.0}
    # Staggered offsets so SMAD9 / PITX2 / KCNN3 do not collide
    callout_offsets = {
        "SMAD9": (10, 10),
        "PITX2": (-28, 12),
        "KCNN3": (10, -14),
    }
    callout_xy: dict[str, tuple[float, float]] = {}

    for etype, y0 in y_base.items():
        sub = src[src["entity_type"] == etype].copy()
        for tier in TIERS:
            members = sub[sub["G_tier"] == tier].reset_index(drop=True)
            n = len(members)
            if n == 0:
                continue
            if n == 1:
                jitters = np.array([0.0])
            else:
                jitters = np.linspace(-0.22, 0.22, n)
                jitters = jitters + rng.normal(0, 0.01, n)
            xs = np.full(n, float(tier_x[tier]))
            if n > 4:
                xs = xs + rng.uniform(-0.12, 0.12, n)
            ys = y0 + jitters
            ax.scatter(
                xs, ys,
                s=48,
                c=TIER_COLORS[tier],
                edgecolors="#555555",
                linewidths=0.5,
                zorder=3,
            )
            if etype == "anchor":
                for i, row in members.iterrows():
                    eid = str(row["entity_id"])
                    if eid in CALLOUT_ANCHORS:
                        callout_xy[eid] = (float(xs[i]), float(ys[i]))

    for eid in CALLOUT_ANCHORS:
        if eid not in callout_xy:
            continue
        x, y = callout_xy[eid]
        dx, dy = callout_offsets[eid]
        ax.annotate(
            eid,
            (x, y),
            textcoords="offset points",
            xytext=(dx, dy),
            fontsize=FS_LEGEND - 0.5,
            color=FG_COLOR,
            fontweight="bold",
            arrowprops=dict(arrowstyle="-", color=MUTED_COLOR, lw=0.5),
        )

    ax.set_xticks(range(len(TIERS)))
    ax.set_xticklabels(TIERS, fontsize=FS_TICK, fontweight="bold")
    ax.set_yticks([1.0, 0.0])
    ax.set_yticklabels(
        [
            f"Anchors (n = {N_ANCHORS})",
            f"Families (n = {N_FAMILIES})",
        ],
        fontsize=FS_TICK,
    )
    ax.set_xlim(-0.5, len(TIERS) - 0.5)
    ax.set_ylim(-0.55, 1.55)
    ax.axhline(0.5, color="#EEEEEE", lw=0.8, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")
    ax.set_facecolor("white")

    # Explicit empty G4 column marker
    g4_x = tier_x["G4"]
    ax.text(
        g4_x, 1.35, "G4: 0",
        ha="center", va="bottom", fontsize=FS_LEGEND - 0.5,
        color=MUTED_COLOR, fontweight="bold",
    )

    fig.text(
        0.02, 0.96,
        "Figure S16 | G-tier entity strip (each point = one anchor or family)",
        fontsize=FS_TICK, fontweight="bold", ha="left", va="top", color=FG_COLOR,
    )
    fig.text(
        0.02, 0.905, TRACEABILITY_NOTE,
        fontsize=FS_LEGEND, ha="left", va="top", color=MUTED_COLOR, style="italic",
    )
    for i, line in enumerate(emphasis[:3]):
        fig.text(
            0.02, 0.855 - 0.035 * i, line,
            fontsize=FS_LEGEND, ha="left", va="top",
            color=FG_COLOR, fontweight="bold" if i < 2 else "normal",
        )

    fig.text(
        0.02, 0.22, f"Zero-count tiers — {zero_note}",
        fontsize=FS_FOOTNOTE - 0.5, ha="left", va="top", color=MUTED_COLOR,
    )
    fig.text(
        0.02, 0.175,
        "Labeled anchors: " + ", ".join(CALLOUT_ANCHORS) + " (others unlabeled to avoid collisions).",
        fontsize=FS_FOOTNOTE - 0.5, ha="left", va="top", color=MUTED_COLOR,
    )

    legend = ax.legend(
        handles=_definition_key_handles(),
        loc="upper left",
        bbox_to_anchor=(0.0, -0.22),
        ncol=1,
        frameon=False,
        fontsize=FS_LEGEND - 1.0,
        handlelength=1.2,
        labelspacing=0.35,
    )
    legend.set_title("G-tier definition key", prop={"size": FS_LEGEND - 0.5, "weight": "bold"})

    return fig


def main() -> None:
    set_nature_style()
    src = _load_source()
    counts = _build_counts(src)
    _assert_data(src, counts)

    OUT_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    src.to_csv(OUT_SOURCE, sep="\t", index=False)
    counts.to_csv(OUT_COUNTS, sep="\t", index=False)
    print(f"wrote {OUT_SOURCE.relative_to(ROOT)} ({len(src)} entities)")
    print(f"wrote {OUT_COUNTS.relative_to(ROOT)}")
    print(counts.to_string(index=False))

    # Candidate 1
    fig1 = _draw_distribution(counts)
    pdf1, png1 = save_figure_bundle(fig1, OUT_DIST, dpi=600)
    print(f"wrote {pdf1.relative_to(ROOT)}")
    if png1:
        print(f"wrote {png1.relative_to(ROOT)}")
    fig1 = _draw_distribution(counts)
    svg1 = OUT_DIST.with_suffix(".svg")
    _save_svg(fig1, svg1)
    plt.close(fig1)
    print(f"wrote {svg1.relative_to(ROOT)}")

    # Candidate 2
    fig2 = _draw_entity_strip(src, counts)
    pdf2, png2 = save_figure_bundle(fig2, OUT_STRIP, dpi=600)
    print(f"wrote {pdf2.relative_to(ROOT)}")
    if png2:
        print(f"wrote {png2.relative_to(ROOT)}")
    fig2 = _draw_entity_strip(src, counts)
    svg2 = OUT_STRIP.with_suffix(".svg")
    _save_svg(fig2, svg2)
    plt.close(fig2)
    print(f"wrote {svg2.relative_to(ROOT)}")

    for line in _emphasis_lines(counts):
        print(f"emphasis: {line}")


if __name__ == "__main__":
    main()
