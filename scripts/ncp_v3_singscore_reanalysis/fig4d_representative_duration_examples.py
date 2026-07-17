#!/usr/bin/env python3
"""Corrected Figure 4D — representative LAA duration trends and regional sensitivity.

Full-width three-facet panel (F9, F2, F15). Scatter only; no locked trend-guide
exists in the corrected pipeline, so no fit line is drawn.
"""
from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
OUT_NAME = "fig4d_representative_duration_examples"

SCORES_TSV = RESULTS / "gse41177_family_sample_scores.tsv"
LAA_TSV = RESULTS / "gse41177_family_duration_laa_primary.tsv"
LAPV_TSV = RESULTS / "gse41177_family_duration_lapv_sensitivity.tsv"
CHRONICITY_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
# Corrected eligibility freeze (duration_months + primary_duration_eligible).
# Required to select the 12 AF-LAA primary samples; not a legacy duration table.
META_TSV = RESULTS / "gse41177_duration_metadata.tsv"

WIDTH_MM = 178.0
HEIGHT_MM = 72.0
MARGIN_L_MM = 16.0
MARGIN_R_MM = 5.0
MARGIN_TOP_MM = 1.5
MARGIN_BOT_MM = 1.5
GUTTER_MM = 4.0
N_FACETS = 3

TITLE_H_MM = 10.0
HEADER_H_MM = 9.5
STATS_H_MM = 7.5
FOOTER_H_MM = 6.5
XLABEL_H_MM = 4.0
TICK_GAP_MM = 3.5
PLOT_BOTTOM_MM = FOOTER_H_MM + XLABEL_H_MM + STATS_H_MM + TICK_GAP_MM
PLOT_TOP_MM = HEIGHT_MM - TITLE_H_MM - HEADER_H_MM - 1.0
PLOT_H_MM = PLOT_TOP_MM - PLOT_BOTTOM_MM

# Fig 4C greys for Strong / Weak badges
STRONG_COLOR = "#2E2E2E"
WEAK_COLOR = "#7A7A7A"
STABLE_COLOR = "#C5C5C5"
TAG_OUTLINE = "#8A8A8A"
TAG_TEXT = "#444444"
# Fig 4B LA/RA strip colors
LA_COLOR = "#9C5550"
RA_COLOR = "#4E6B89"
POINT_COLOR = "#2E2E2E"
POINT_EDGE = "#D8D8D8"
POINT_SIZE = 28.0 * 0.95  # ~5% smaller
ZERO_LINE_COLOR = "#E8E8E8"
FG_COLOR = "#1A1A1A"
MUTED_COLOR = "#555555"
FAINT_COLOR = "#777777"
SPINE_COLOR = "#666666"
WHITE = "#FFFFFF"

FACET_ORDER = [9, 2, 15]
EXPECTED = {
    9: {
        "class": "Strong_trend",
        "rho_laa": -0.420,
        "rho_lapv": -0.161,
        "regional": "Same sign",
        "label": "HS-GAG / bone / GPCR",
        "axis": "RA",
    },
    2: {
        "class": "Strong_trend",
        "rho_laa": -0.564,
        "rho_lapv": 0.196,
        "regional": "Opposite sign",
        "label": "Mito Ca²⁺ transport",
        "axis": "LA",
    },
    15: {
        "class": "Weak_trend",
        "rho_laa": -0.385,
        "rho_lapv": -0.186,
        "regional": "Same sign",
        "label": "NTRK2 / TrkB",
        "axis": "RA",
    },
}

MIN_FONT_PT = 5.0
BBOX_TOL_MM = 0.35
OVERLAP_TOL_MM2 = 0.25
Y_PAD_FRAC = 0.08
JITTER_HALFWIDTH = 1.8  # months; only for duplicate durations


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "font.size": 6.5,
        "axes.linewidth": 0.45,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })


def minus_fmt(value: float, decimals: int) -> str:
    """Format with a true minus sign for negative values."""
    mag = f"{abs(value):.{decimals}f}"
    return f"−{mag}" if value < 0 else mag


def regional_label(agreement: bool) -> str:
    return "Same sign" if bool(agreement) else "Opposite sign"


def class_from_rho(rho: float) -> str:
    abs_rho = abs(float(rho))
    if abs_rho >= 0.4:
        return "Strong_trend"
    if abs_rho >= 0.2:
        return "Weak_trend"
    return "Stable"


def class_badge_color(chronicity_class: str) -> str:
    if chronicity_class == "Strong_trend":
        return STRONG_COLOR
    if chronicity_class == "Weak_trend":
        return WEAK_COLOR
    return STABLE_COLOR


def class_display(chronicity_class: str) -> str:
    return chronicity_class.replace("_", "-")


def deterministic_x_jitter(sample_ids: list[str], durations: np.ndarray) -> np.ndarray:
    """Small horizontal jitter only where duration values are duplicated."""
    jitter = np.zeros(len(durations), dtype=float)
    for duration in np.unique(durations):
        idx = np.where(durations == duration)[0]
        if len(idx) < 2:
            continue
        for rank, i in enumerate(sorted(idx, key=lambda j: sample_ids[j])):
            # Symmetric offsets for 2 points; spread for more.
            if len(idx) == 2:
                offset = -JITTER_HALFWIDTH if rank == 0 else JITTER_HALFWIDTH
            else:
                offset = ((rank / (len(idx) - 1)) * 2.0 - 1.0) * JITTER_HALFWIDTH
            jitter[i] = offset
    return jitter


def load_inputs() -> dict[str, pd.DataFrame]:
    for path in (SCORES_TSV, LAA_TSV, LAPV_TSV, CHRONICITY_TSV, LABELS_TSV, META_TSV):
        require(path.is_file(), f"missing required input: {path}")

    scores = pd.read_csv(SCORES_TSV, sep="\t")
    laa = pd.read_csv(LAA_TSV, sep="\t").set_index("family_id")
    lapv = pd.read_csv(LAPV_TSV, sep="\t").set_index("family_id")
    chronicity = pd.read_csv(CHRONICITY_TSV, sep="\t").set_index("family_id")
    labels = pd.read_csv(LABELS_TSV, sep="\t").set_index("family_id")
    meta = pd.read_csv(META_TSV, sep="\t")

    require("python_exact_reference" in set(scores["scorer"]), "unexpected scorer lineage")
    require(int(meta["primary_duration_eligible"].sum()) == 12, "expected 12 primary-eligible samples")
    return {
        "scores": scores,
        "laa": laa,
        "lapv": lapv,
        "chronicity": chronicity,
        "labels": labels,
        "meta": meta,
    }


def build_source(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    primary = data["meta"].loc[data["meta"]["primary_duration_eligible"]].copy()
    require(len(primary) == 12, "primary subset must have 12 rows")
    require((primary["rhythm"] == "AF").all(), "primary samples must be AF")
    require((primary["region"] == "LAA").all(), "primary samples must be LAA")
    require(primary["duration_months"].notna().all(), "primary durations must be numeric")
    require((primary["duration_months"] > 0).all(), "SR/zero-duration samples must not enter primary set")

    # Assert no SR samples leak in
    sr_ids = set(data["meta"].loc[data["meta"]["rhythm"] == "SR", "sample_id"])
    require(sr_ids.isdisjoint(set(primary["sample_id"])), "SR sample in primary-eligible set")

    rows = []
    for facet_rank, fid in enumerate(FACET_ORDER, start=1):
        sub = data["scores"][
            (data["scores"]["family_id"] == fid)
            & (data["scores"]["region"] == "LAA")
            & data["scores"]["sample_id"].isin(primary["sample_id"])
        ].merge(
            primary[["sample_id", "donor_id", "duration_months", "rhythm", "exclusion_reason"]],
            on="sample_id",
            how="inner",
        )
        require(len(sub) == 12, f"F{fid}: expected 12 AF-LAA points, got {len(sub)}")
        require(not sub["sample_id"].isin(sr_ids).any(), f"F{fid}: SR sample rendered")
        require((sub["rhythm"] == "AF").all(), f"F{fid}: non-AF sample rendered")
        require(sub["duration_months"].notna().all(), f"F{fid}: missing duration")

        laa_row = data["laa"].loc[fid]
        lapv_row = data["lapv"].loc[fid]
        cls_row = data["chronicity"].loc[fid]
        label_row = data["labels"].loc[fid]

        for _, r in sub.iterrows():
            rows.append({
                "facet_rank": facet_rank,
                "family_id": fid,
                "short_label": str(label_row["short_label"]),
                "axis_group": str(label_row["axis_group"]),
                "sample_id": str(r["sample_id"]),
                "donor_id": str(r["donor_id"]),
                "region": "LAA",
                "rhythm": "AF",
                "duration_months": float(r["duration_months"]),
                "score": float(r["score"]),
                "chronicity_class": str(cls_row["chronicity_class"]),
                "rho_laa": float(laa_row["rho_spearman"]),
                "ci_lo_laa": float(laa_row["bootstrap_ci_lo"]),
                "ci_hi_laa": float(laa_row["bootstrap_ci_hi"]),
                "p_perm_laa": float(laa_row["p_perm"]),
                "rho_lapv": float(lapv_row["rho_spearman"]),
                "regional_label": regional_label(cls_row["regional_sign_agreement"]),
                "n_samples": 12,
            })

    source = pd.DataFrame(rows)
    require(len(source) == 36, "source must have 36 rows (3 × 12)")
    return source


def assert_scientific_contract(source: pd.DataFrame, data: dict[str, pd.DataFrame]) -> None:
    for fid, exp in EXPECTED.items():
        cls_row = data["chronicity"].loc[fid]
        laa_row = data["laa"].loc[fid]
        lapv_row = data["lapv"].loc[fid]
        rho_laa = float(laa_row["rho_spearman"])
        recomputed = class_from_rho(rho_laa)
        require(recomputed == exp["class"], f"F{fid} class from unrounded rho mismatch")
        require(str(cls_row["chronicity_class"]) == recomputed,
                f"F{fid} chronicity_class does not match unrounded |rho_LAA| bins")
        require(abs(rho_laa - exp["rho_laa"]) < 0.002, f"F{fid} rho_LAA mismatch")
        require(abs(float(lapv_row["rho_spearman"]) - exp["rho_lapv"]) < 0.002, f"F{fid} rho_LAPV mismatch")
        require(regional_label(cls_row["regional_sign_agreement"]) == exp["regional"], f"F{fid} regional mismatch")
        require(int(laa_row["n_samples"]) == 12, f"F{fid} LAA n != 12")
        require(int(lapv_row["n_samples"]) == 12, f"F{fid} LAPV n != 12")
        facet = source[source["family_id"] == fid]
        require(len(facet) == 12, f"F{fid} facet n != 12")
        require(facet["sample_id"].nunique() == 12, f"F{fid} duplicate sample_ids")


def build_summary(source: pd.DataFrame, data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for facet_rank, fid in enumerate(FACET_ORDER, start=1):
        facet = source[source["family_id"] == fid]
        cls_row = data["chronicity"].loc[fid]
        laa_row = data["laa"].loc[fid]
        lapv_row = data["lapv"].loc[fid]
        label_row = data["labels"].loc[fid]
        rows.append({
            "facet_rank": facet_rank,
            "family_id": fid,
            "short_label": str(label_row["short_label"]),
            "axis_group": str(label_row["axis_group"]),
            "chronicity_class": str(cls_row["chronicity_class"]),
            "regional_label": regional_label(cls_row["regional_sign_agreement"]),
            "n_samples": 12,
            "rho_laa": float(laa_row["rho_spearman"]),
            "ci_lo_laa": float(laa_row["bootstrap_ci_lo"]),
            "ci_hi_laa": float(laa_row["bootstrap_ci_hi"]),
            "p_perm_laa": float(laa_row["p_perm"]),
            "rho_lapv": float(lapv_row["rho_spearman"]),
            "score_min": float(facet["score"].min()),
            "score_max": float(facet["score"].max()),
            "duration_min": float(facet["duration_months"].min()),
            "duration_max": float(facet["duration_months"].max()),
        })
    return pd.DataFrame(rows)


def build_axis_contract(source: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    y_min_raw = float(source["score"].min())
    y_max_raw = float(source["score"].max())
    pad = max((y_max_raw - y_min_raw) * Y_PAD_FRAC, 0.02)
    y_min = y_min_raw - pad
    y_max = y_max_raw + pad
    # Nice 0.1-spaced ticks strictly inside the locked limits
    y_lo_tick = np.ceil(y_min * 10.0 - 1e-9) / 10.0
    y_hi_tick = np.floor(y_max * 10.0 + 1e-9) / 10.0
    if y_hi_tick < y_lo_tick:
        y_ticks = np.array([round((y_min + y_max) / 2.0, 2)])
    else:
        y_ticks = np.round(np.arange(y_lo_tick, y_hi_tick + 0.05, 0.1), 2)

    x_max = float(source["duration_months"].max())
    x_lo = 0.0
    x_hi = x_max + max(x_max * 0.03, 3.0)
    x_ticks = np.array([0.0, 50.0, 100.0, 150.0])
    x_ticks = x_ticks[(x_ticks >= x_lo - 1e-9) & (x_ticks <= x_hi + 1e-9)]

    return pd.DataFrame([{
        "contract_id": "figure4d_v1",
        "panel_width_mm": WIDTH_MM,
        "panel_height_mm": HEIGHT_MM,
        "n_facets": N_FACETS,
        "facet_order": ",".join(f"F{f}" for f in FACET_ORDER),
        "gutter_mm": GUTTER_MM,
        "x_label": "AF duration (months)",
        "y_label": "Family activity score (centered rank units)",
        "x_min": x_lo,
        "x_max": x_hi,
        "x_ticks": ",".join(str(int(t)) for t in x_ticks),
        "y_min": y_min,
        "y_max": y_max,
        "y_ticks": ",".join(f"{t:.2f}" for t in y_ticks),
        "point_color": POINT_COLOR,
        "jitter_halfwidth_months": JITTER_HALFWIDTH,
        "trend_guide": "omitted_no_locked_implementation",
        "n_points_per_facet": 12,
        "score_union_min": y_min_raw,
        "score_union_max": y_max_raw,
    }])


def write_tables(source: pd.DataFrame, summary: pd.DataFrame, axis: pd.DataFrame) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    source_path = RESULTS / f"{OUT_NAME}_source.tsv"
    summary_path = RESULTS / f"{OUT_NAME}_summary.tsv"
    axis_path = RESULTS / "figure4d_axis_contract.tsv"
    source.to_csv(source_path, sep="\t", index=False, float_format="%.6f")
    summary.to_csv(summary_path, sep="\t", index=False, float_format="%.6f")
    axis.to_csv(axis_path, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {source_path.relative_to(ROOT)}")
    print(f"wrote {summary_path.relative_to(ROOT)}")
    print(f"wrote {axis_path.relative_to(ROOT)}")


def mm_axes(fig: plt.Figure, x_mm: float, y_mm: float, w_mm: float, h_mm: float) -> plt.Axes:
    return fig.add_axes([x_mm / WIDTH_MM, y_mm / HEIGHT_MM, w_mm / WIDTH_MM, h_mm / HEIGHT_MM])


def facet_geometry() -> list[tuple[float, float]]:
    usable = WIDTH_MM - MARGIN_L_MM - MARGIN_R_MM
    facet_w = (usable - (N_FACETS - 1) * GUTTER_MM) / N_FACETS
    return [(MARGIN_L_MM + i * (facet_w + GUTTER_MM), facet_w) for i in range(N_FACETS)]


def stats_text(row: pd.Series) -> str:
    return (
        f"LAA ρ = {minus_fmt(float(row['rho_laa']), 3)} "
        f"[95% CI {minus_fmt(float(row['ci_lo_laa']), 3)}, {minus_fmt(float(row['ci_hi_laa']), 3)}]\n"
        f"permutation p = {float(row['p_perm_laa']):.4f} · "
        f"LAPV ρ = {minus_fmt(float(row['rho_lapv']), 3)}"
    )


def draw_tag(
    ax: plt.Axes,
    x: float,
    y: float,
    text: str,
    *,
    facecolor: str,
    edgecolor: str,
    textcolor: str,
    fontsize: float = 5.2,
) -> float:
    """Draw a compact rounded tag; return right edge x in axis mm coords."""
    # Approximate width from character count (mm at ~5.2 pt)
    width = max(14.0, 1.55 * len(text) + 2.8)
    height = 2.6
    ax.add_patch(mpatches.FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.12,rounding_size=0.45",
        facecolor=facecolor, edgecolor=edgecolor, linewidth=0.55,
        mutation_aspect=1.0, clip_on=False,
    ))
    ax.text(x + width / 2, y + height / 2, text, fontsize=fontsize,
            ha="center", va="center", color=textcolor, fontweight="bold",
            clip_on=False)
    return x + width


def render(source: pd.DataFrame, summary: pd.DataFrame, axis: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    fig.patch.set_facecolor(WHITE)
    fig.set_dpi(300)

    ax_row = axis.iloc[0]
    x_lo, x_hi = float(ax_row["x_min"]), float(ax_row["x_max"])
    y_lo, y_hi = float(ax_row["y_min"]), float(ax_row["y_max"])
    x_ticks = [float(t) for t in str(ax_row["x_ticks"]).split(",") if t]
    y_ticks = [float(t) for t in str(ax_row["y_ticks"]).split(",") if t]
    geoms = facet_geometry()

    # Title band
    title_ax = mm_axes(fig, MARGIN_L_MM, HEIGHT_MM - TITLE_H_MM - MARGIN_TOP_MM,
                       WIDTH_MM - MARGIN_L_MM - MARGIN_R_MM, TITLE_H_MM)
    title_ax.set_xlim(0, 1)
    title_ax.set_ylim(0, 1)
    title_ax.axis("off")
    title_ax.text(0, 0.92, "D", fontsize=10.5, fontweight="bold", va="top", color=FG_COLOR)
    title_ax.text(0.035, 0.92, "Representative LAA duration trends and regional sensitivity",
                  fontsize=8.8, fontweight="bold", va="top", color=FG_COLOR)
    title_ax.text(0.035, 0.35,
                  "Strong and Weak LAA trends with concordant or discordant LAPV signs",
                  fontsize=6.5, va="top", color=MUTED_COLOR)

    plot_axes: list[plt.Axes] = []
    n_points_drawn: dict[int, int] = {}
    for i, fid in enumerate(FACET_ORDER):
        x0, fw = geoms[i]
        sum_row = summary.loc[summary["family_id"] == fid].iloc[0]
        facet = source[source["family_id"] == fid].copy().reset_index(drop=True)

        # Header: direction strip + title + two separate tags
        head = mm_axes(fig, x0, PLOT_TOP_MM, fw, HEADER_H_MM)
        head.set_xlim(0, fw)
        head.set_ylim(0, HEADER_H_MM)
        head.axis("off")
        strip_color = LA_COLOR if sum_row["axis_group"] == "LA" else RA_COLOR
        head.add_patch(mpatches.Rectangle((0, 1.8), 2.2, 5.8, facecolor=strip_color, edgecolor="none"))
        head.text(1.1, 4.7, sum_row["axis_group"], fontsize=5.5, color=WHITE,
                  ha="center", va="center", fontweight="bold", rotation=90)
        head.text(3.2, 7.6, f"F{fid}  {display_label(fid, str(sum_row['short_label']))}",
                  fontsize=7.0, fontweight="bold", va="top", color=FG_COLOR)

        class_text = class_display(str(sum_row["chronicity_class"]))
        regional_text = str(sum_row["regional_label"])
        tag_y = 1.4
        tag_x = 3.2
        tag_x = draw_tag(
            head, tag_x, tag_y, class_text,
            facecolor=class_badge_color(str(sum_row["chronicity_class"])),
            edgecolor="none",
            textcolor=WHITE,
        )
        tag_x += 1.2
        draw_tag(
            head, tag_x, tag_y, regional_text,
            facecolor=WHITE,
            edgecolor=TAG_OUTLINE,
            textcolor=TAG_TEXT,
        )

        # Scatter
        ax = mm_axes(fig, x0, PLOT_BOTTOM_MM, fw, PLOT_H_MM)
        durations = facet["duration_months"].to_numpy(dtype=float)
        scores = facet["score"].to_numpy(dtype=float)
        sample_ids = facet["sample_id"].astype(str).tolist()
        require(len(sample_ids) == 12, f"F{fid}: expected 12 observations, got {len(sample_ids)}")
        x_plot = durations + deterministic_x_jitter(sample_ids, durations)
        ax.axhline(0.0, color=ZERO_LINE_COLOR, linewidth=0.45, zorder=1)
        ax.scatter(
            x_plot, scores,
            s=POINT_SIZE, c=POINT_COLOR, edgecolors=POINT_EDGE, linewidths=0.28,
            zorder=3, alpha=0.95,
        )
        n_points_drawn[fid] = len(scores)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.tick_params(labelsize=5.5, length=2.0, width=0.4, colors=SPINE_COLOR,
                       pad=1.5)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(SPINE_COLOR)
            ax.spines[spine].set_linewidth(0.45)
        if i == 0:
            ax.set_ylabel("Family activity score\n(centered rank units)",
                          fontsize=6.5, color=FG_COLOR, labelpad=4)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("")
        plot_axes.append(ax)

        # Stats under facet (exactly two lines)
        stats_ax = mm_axes(fig, x0, FOOTER_H_MM + XLABEL_H_MM, fw, STATS_H_MM)
        stats_ax.set_xlim(0, 1)
        stats_ax.set_ylim(0, 1)
        stats_ax.axis("off")
        stats_ax.text(0.0, 0.88, stats_text(sum_row), fontsize=5.2, va="top", ha="left",
                      color=MUTED_COLOR, linespacing=1.45, family="sans-serif")

    require(all(n == 12 for n in n_points_drawn.values()),
            f"point counts drifted: {n_points_drawn}")

    # Shared x-axis label centered under all facets
    xlab = mm_axes(fig, MARGIN_L_MM, FOOTER_H_MM, WIDTH_MM - MARGIN_L_MM - MARGIN_R_MM, XLABEL_H_MM)
    xlab.set_xlim(0, 1)
    xlab.set_ylim(0, 1)
    xlab.axis("off")
    xlab.text(0.5, 0.55, "AF duration (months)", fontsize=6.5, ha="center", va="center",
              color=FG_COLOR)

    # Shared-axis contract assertion at render time
    for ax in plot_axes:
        xl, xh = ax.get_xlim()
        yl, yh = ax.get_ylim()
        require(abs(xl - x_lo) < 1e-9 and abs(xh - x_hi) < 1e-9, "facet x limits drifted")
        require(abs(yl - y_lo) < 1e-6 and abs(yh - y_hi) < 1e-6, "facet y limits drifted")
        require(list(ax.get_xticks()) == list(x_ticks) or np.allclose(ax.get_xticks(), x_ticks),
                "facet x ticks drifted")
        # y ticks may include extras from MaxNLocator; force exact shared ticks
        ax.set_yticks(y_ticks)
        require(np.allclose(ax.get_yticks(), y_ticks), "facet y ticks drifted")

    # Footer
    foot = mm_axes(fig, MARGIN_L_MM, 0.6, WIDTH_MM - MARGIN_L_MM - MARGIN_R_MM, FOOTER_H_MM - 0.8)
    foot.set_xlim(0, 1)
    foot.set_ylim(0, 1)
    foot.axis("off")
    foot.text(
        0, 0.70,
        "12 AF-LAA samples per facet; classes use full-precision |ρ_LAA|;\n"
        "LAPV is regional sensitivity only.",
        fontsize=5.0, va="top", color=FAINT_COLOR, linespacing=1.35,
    )
    return fig


def text_bbox_mm(artist: mpl.text.Text, fig: plt.Figure) -> tuple[float, float, float, float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    return (
        float(bbox.x0) * WIDTH_MM,
        float(bbox.y0) * HEIGHT_MM,
        float(bbox.x1) * WIDTH_MM,
        float(bbox.y1) * HEIGHT_MM,
    )


def display_label(fid: int, short_label: str) -> str:
    """Publication label with mathtext for Ca2+."""
    if fid == 2:
        return r"Mito Ca$^{2+}$ transport"
    return short_label


def validate_layout(fig: plt.Figure) -> dict[str, float | int]:
    entries: list[tuple[str, mpl.text.Text, tuple[float, float, float, float]]] = []
    for ax in fig.axes:
        for child in ax.get_children():
            if isinstance(child, mpl.text.Text) and child.get_text().strip():
                if not child.get_visible():
                    continue
                entries.append((child.get_text().strip(), child, text_bbox_mm(child, fig)))
        if not ax.axison:
            continue
        for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            if lab.get_text().strip() and lab.get_visible():
                entries.append((lab.get_text().strip(), lab, text_bbox_mm(lab, fig)))

    require(entries, "no visible text for layout QC")
    min_font = min(child.get_fontsize() for _, child, _ in entries)
    require(min_font >= MIN_FONT_PT - 0.05, f"min font {min_font:.2f} < {MIN_FONT_PT}")

    violations = 0
    for label, _, bbox in entries:
        if bbox[0] < -BBOX_TOL_MM or bbox[2] > WIDTH_MM + BBOX_TOL_MM:
            violations += 1
            require(False, f"text clipped horizontally: {label[:60]!r} bbox={bbox}")
        if bbox[1] < -BBOX_TOL_MM or bbox[3] > HEIGHT_MM + BBOX_TOL_MM:
            violations += 1
            require(False, f"text clipped vertically: {label[:60]!r} bbox={bbox}")

    return {
        "n_text_objects": len(entries),
        "min_font_pt": float(min_font),
        "n_bbox_violations": violations,
    }


def save_outputs(fig: plt.Figure) -> dict[str, Path]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    paths = {
        "pdf": FIGURES / f"{OUT_NAME}.pdf",
        "svg": FIGURES / f"{OUT_NAME}.svg",
        "png": FIGURES / f"{OUT_NAME}.png",
        "preview": FIGURES / f"{OUT_NAME}_preview.png",
    }
    fig.savefig(paths["pdf"], bbox_inches=None, pad_inches=0.02)
    fig.savefig(paths["svg"], bbox_inches=None, pad_inches=0.02)
    fig.savefig(paths["png"], dpi=600, bbox_inches=None, pad_inches=0.02)
    fig.savefig(paths["preview"], dpi=150, bbox_inches=None, pad_inches=0.02)
    plt.close(fig)
    for path in paths.values():
        print(f"wrote {path.relative_to(ROOT)}")
    return paths


def write_qc(
    summary: pd.DataFrame,
    axis: pd.DataFrame,
    layout_qc: dict[str, float | int],
    output_paths: dict[str, Path],
) -> None:
    qc_path = RESULTS / f"{OUT_NAME}_qc.md"
    inputs = {
        str(SCORES_TSV.relative_to(ROOT)): sha256(SCORES_TSV),
        str(LAA_TSV.relative_to(ROOT)): sha256(LAA_TSV),
        str(LAPV_TSV.relative_to(ROOT)): sha256(LAPV_TSV),
        str(CHRONICITY_TSV.relative_to(ROOT)): sha256(CHRONICITY_TSV),
        str(LABELS_TSV.relative_to(ROOT)): sha256(LABELS_TSV),
        str(META_TSV.relative_to(ROOT)): sha256(META_TSV),
    }
    outputs = {str(p.relative_to(ROOT)): sha256(p) for p in output_paths.values()}
    for name in (f"{OUT_NAME}_source.tsv", f"{OUT_NAME}_summary.tsv", "figure4d_axis_contract.tsv"):
        path = RESULTS / name
        outputs[str(path.relative_to(ROOT))] = sha256(path)

    ax = axis.iloc[0]
    lines = [
        "# Figure 4D representative duration examples QC",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "Status: **PASS**",
        "",
        "## Scientific contract",
        "- Facet order: F9, F2, F15",
        "- exactly 12 AF-LAA primary-eligible samples in every facet",
        "- No SR samples rendered; no censored/nonnumeric durations; no LAPV donor points",
        "- No trend-guide / regression line / confidence band",
        "- Chronicity class and regional sign shown as separate tags",
        "- Classes reproduce full-precision |rho_LAA| bins; p/CI/LAPV do not enter class assignment",
        "- Footer: 12 AF-LAA samples per facet; classes use full-precision |ρ_LAA|; LAPV is regional sensitivity only",
        "- Statistics read from corrected LAA primary / LAPV sensitivity / chronicity tables",
        "",
        "## Spot checks",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- F{int(row['family_id'])}: class={row['chronicity_class']}; "
            f"ρ_LAA={float(row['rho_laa']):+.4f}; ρ_LAPV={float(row['rho_lapv']):+.4f}; "
            f"regional={row['regional_label']}; n={int(row['n_samples'])}"
        )
    lines += [
        "",
        "## Axis contract",
        f"- Panel: {float(ax['panel_width_mm']):.1f} × {float(ax['panel_height_mm']):.1f} mm",
        f"- Shared x: [{float(ax['x_min']):.2f}, {float(ax['x_max']):.2f}]; ticks={ax['x_ticks']}",
        f"- Shared y: [{float(ax['y_min']):.4f}, {float(ax['y_max']):.4f}]; ticks={ax['y_ticks']}",
        f"- Identical x/y limits and ticks across all three facets",
        f"- Trend guide: {ax['trend_guide']}",
        "",
        "## Layout QC",
        f"- visible text objects: {layout_qc['n_text_objects']}",
        f"- minimum font size: {layout_qc['min_font_pt']:.2f} pt",
        f"- bbox violations: {layout_qc['n_bbox_violations']}",
        "",
        "## Allowed inputs",
        *[f"- `{p}` — `{d}`" for p, d in inputs.items()],
        "",
        "## Outputs",
        *[f"- `{p}` — `{d}`" for p, d in outputs.items()],
        "",
    ]
    qc_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {qc_path.relative_to(ROOT)}")


def main() -> None:
    set_style()
    data = load_inputs()
    source = build_source(data)
    assert_scientific_contract(source, data)
    summary = build_summary(source, data)
    axis = build_axis_contract(source, summary)
    write_tables(source, summary, axis)
    fig = render(source, summary, axis)
    layout_qc = validate_layout(fig)
    paths = save_outputs(fig)
    write_qc(summary, axis, layout_qc, paths)


if __name__ == "__main__":
    main()
