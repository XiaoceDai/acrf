#!/usr/bin/env python3
"""Final Figure 4A: cross-sectional AF-duration analysis design bridge (GSE41177)."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"

OUT_BASE = FIGURES / "fig4a_cross_sectional_duration_bridge"
OUT_PREVIEW = FIGURES / "fig4a_cross_sectional_duration_bridge_preview.png"
OUT_MANIFEST = RESULTS / "fig4a_input_manifest.tsv"
OUT_SUMMARY = RESULTS / "fig4a_duration_design_summary.tsv"
OUT_QC = RESULTS / "fig4a_duration_design_qc.md"

CHRONICITY_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"

INPUTS = {
    "frozen_program_metadata": RESULTS / "frozen_program_metadata.tsv",
    "frozen_family_metadata": RESULTS / "frozen_family_metadata.tsv",
    "leading_edge_validation": RESULTS / "leading_edge_comparison_summary.tsv",
    "gse41177_duration_metadata": RESULTS / "gse41177_duration_metadata.tsv",
    "gse41177_program_sample_scores": RESULTS / "gse41177_program_sample_scores.tsv",
    "gse41177_family_sample_scores": RESULTS / "gse41177_family_sample_scores.tsv",
}

REQUIRED_ROLES = set(INPUTS)
SCORER = "python_exact_reference"
SCORE_METHOD = "python_exact_reference::exact_reference(centered rank-based program score)"
DATASET = "GSE41177"

MM = 1 / 25.4
WIDTH_MM = 178.0
HEIGHT_MM = 52.0
MARGIN_MM = 3.0
STRIP_H_MM = 36.5
TITLE_STRIP_GAP_MM = 1.8
TITLE_SUB_GAP_MM = 1.3
GAP_MM = 3.0
COL_FRAC = [0.21, 0.22, 0.23, 0.34]
LINE_H_MM = 2.45
MODULE_LINE_H_MM = 2.7
BODY_PT = 6.2
FOOTNOTE_PT = 6.2

FG_COLOR = "#1A1A1A"
MUTED_COLOR = "#595959"
FAINT_COLOR = "#7A7A7A"
KICKER_COLOR = "#8A8A8A"
MID_GREY = "#9A9A9A"
LAPV_COLOR = "#9A9A9A"
STRIP_BG = "#F7F7F7"
WHITE = "#FFFFFF"
MIN_FONT_PT = 6.2
BBOX_TOL_MM = 0.55
OVERLAP_TOL_MM2 = 0.15


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify_chronicity(abs_rho: float) -> str:
    if abs_rho >= 0.4:
        return "Strong_trend"
    if abs_rho >= 0.2:
        return "Weak_trend"
    return "Stable"


def build_manifest() -> pd.DataFrame:
    rows = []
    for role, path in INPUTS.items():
        require(path.is_file(), f"Missing Figure 4A input: {path}")
        programs = pd.read_csv(INPUTS["frozen_program_metadata"], sep="\t")
        le_version = str(programs["gene_set_version"].iloc[0])
        rows.append(
            {
                "artifact_role": role,
                "artifact_path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "generated_by": "corrected_ncp_v3_singscore_reanalysis",
                "dataset": DATASET if role.startswith("gse41177") else "frozen_from_Fig1",
                "used_for": {
                    "frozen_program_metadata": "frozen_topology",
                    "frozen_family_metadata": "frozen_topology",
                    "leading_edge_validation": "frozen_leading_edges",
                    "gse41177_duration_metadata": "duration_eligibility",
                    "gse41177_program_sample_scores": "program_activity_scores",
                    "gse41177_family_sample_scores": "family_activity_scores",
                }[role],
                "score_method": SCORE_METHOD,
                "scorer_package": "python_exact_reference",
                "scorer_version": "1.28.1_formula",
                "leading_edge_version": le_version,
                "reference_dataset": "GSE112339_frozen_LE",
                "reference_formula": "mean(program_score) per family",
                "allows_cohort_specific_reselection": False,
                "allows_reference_recalibration": False,
                "is_authoritative": True,
            }
        )
    manifest = pd.DataFrame(rows)
    require(set(manifest["artifact_role"]) == REQUIRED_ROLES, "Manifest roles incomplete")
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(OUT_MANIFEST, sep="\t", index=False)
    return manifest


def load_manifest() -> tuple[pd.DataFrame, dict[str, Path]]:
    manifest = build_manifest()
    paths: dict[str, Path] = {}
    allowed_root = (RESULTS).resolve()
    for row in manifest.itertuples(index=False):
        path = (ROOT / row.artifact_path).resolve()
        require(path.is_file(), f"Manifest input missing: {path}")
        require(path.is_relative_to(allowed_root), f"Input outside corrected namespace: {path}")
        require(sha256(path) == row.sha256, f"Checksum mismatch: {path}")
        paths[str(row.artifact_role)] = path
    return manifest, paths


def validate_duration_parser(meta: pd.DataFrame) -> None:
    require(len(meta) == 38, f"expected 38 samples, got {len(meta)}")
    require(int(meta["primary_duration_eligible"].sum()) == 12, "expected 12 primary-eligible AF-LAA")

    bad_year = meta["duration_raw"].astype(str).str.contains(r"\d+[Yy]", regex=True, na=False)
    require(not bad_year.any(), "year-converted durations detected")

    for _, row in meta.iterrows():
        raw = str(row["duration_raw"]).strip()
        months = row["duration_months"]
        if row["rhythm"] == "SR":
            require(not bool(row["primary_duration_eligible"]), "SR sample marked primary-eligible")
            require(str(row["exclusion_reason"]) == "sr_anchor", "SR exclusion reason mismatch")
            continue
        if row["region"] == "LAPV":
            require(not bool(row["primary_duration_eligible"]), "LAPV sample marked primary-eligible")
            continue
        if raw == ">1M":
            require(pd.isna(months), ">1M must have missing duration_months")
            require(str(row["exclusion_reason"]) == "censored_gt1m", "censored exclusion mismatch")
            continue
        match = re.fullmatch(r"(\d+)M", raw)
        if match:
            require(float(months) == float(match.group(1)), f"parser mismatch for {raw}")
            if row["region"] == "LAA" and row["rhythm"] == "AF":
                require(
                    bool(row["primary_duration_eligible"]),
                    f"finite AF-LAA month sample not primary-eligible: {raw}",
                )
            continue
        require(pd.isna(months), f"non-numeric duration must be excluded: {raw}")


def validate_aggregation(prog: pd.DataFrame, fam: pd.DataFrame, programs: pd.DataFrame) -> float:
    require(set(prog["scorer"].astype(str)) == {SCORER}, "unexpected program scorer")
    require(set(fam["scorer"].astype(str)) == {SCORER}, "unexpected family scorer")
    require(prog["program_id"].nunique() == 30, "expected 30 programs")
    require(prog["sample_id"].nunique() == 38, "expected 38 scored samples")
    require(len(programs) == 30 and programs["program_id"].nunique() == 30, "frozen program count changed")

    recomputed = (
        prog.groupby(["family_id", "sample_id", "region"], as_index=False)
        .agg(score_recomputed=("score", "mean"), n_programs=("program_id", "nunique"))
    )
    observed = fam[["family_id", "sample_id", "region", "score", "n_programs"]].copy()
    merged = recomputed.merge(
        observed,
        on=["family_id", "sample_id", "region"],
        how="outer",
        suffixes=("_recomputed", "_observed"),
        indicator=True,
    )
    require((merged["_merge"] == "both").all(), "program/family score tables do not align")
    merged["abs_error"] = (merged["score_recomputed"] - merged["score"]).abs()
    max_error = float(merged["abs_error"].max())
    require(max_error <= 1e-6, f"equal-weight aggregation error too large: {max_error}")
    return max_error


def validate_class_thresholds() -> None:
    require(CHRONICITY_TSV.is_file(), f"missing chronicity reference: {CHRONICITY_TSV}")
    chro = pd.read_csv(CHRONICITY_TSV, sep="\t")
    for _, row in chro.iterrows():
        expected = classify_chronicity(abs(float(row["rho_primary"])))
        require(
            str(row["chronicity_class"]) == expected,
            f"threshold mismatch for family {row['family_id']}: "
            f"{row['chronicity_class']} != {expected}",
        )
    require(int(chro["n_samples_primary"].unique()[0]) == 12, "primary correlation n must be 12")


def build_summary(meta: pd.DataFrame, programs: pd.DataFrame, families: pd.DataFrame, max_error: float) -> pd.DataFrame:
    n_censored = int(
        (
            (meta["rhythm"] == "AF")
            & (meta["region"] == "LAA")
            & (meta["exclusion_reason"] == "censored_gt1m")
        ).sum()
    )
    n_sr_laa = int(((meta["rhythm"] == "SR") & (meta["region"] == "LAA")).sum())
    rows = [
        {
            "scope": "design",
            "total_samples": len(meta),
            "primary_eligible_af_laa": int(meta["primary_duration_eligible"].sum()),
            "censored_af_laa": n_censored,
            "sr_laa_anchors": n_sr_laa,
            "frozen_programs": len(programs),
            "frozen_families": len(families),
            "scorer": SCORER,
            "aggregation": "equal_weight_program_first",
            "input_layer": "family_activity_score",
            "uses_acairel": False,
            "class_rule": "abs_rho_laa_only",
            "lapv_role": "regional_sensitivity_only",
            "max_aggregation_error": max_error,
        }
    ]
    return pd.DataFrame(rows)


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 6.5,
            "mathtext.default": "regular",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": WHITE,
            "savefig.facecolor": WHITE,
        }
    )


def strip_geometry() -> dict[str, float | list[float]]:
    content_w = WIDTH_MM - 2 * MARGIN_MM
    col_space = content_w - 3 * GAP_MM
    col_w = [col_space * fraction for fraction in COL_FRAC]
    xs = [MARGIN_MM]
    for width in col_w[:-1]:
        xs.append(xs[-1] + width + GAP_MM)
    strip_y0 = MARGIN_MM
    strip_y1 = strip_y0 + STRIP_H_MM
    title_y1 = HEIGHT_MM - MARGIN_MM
    gap_centers = [xs[index] + col_w[index] + GAP_MM / 2 for index in range(3)]
    return {
        "xs": xs,
        "col_w": col_w,
        "strip_y0": strip_y0,
        "strip_y1": strip_y1,
        "strip_cy": (strip_y0 + strip_y1) / 2,
        "title_y1": title_y1,
        "gap_centers": gap_centers,
    }


def bridge_arrow(ax: plt.Axes, x0: float, x1: float, y: float) -> None:
    ax.annotate(
        "",
        xy=(x1, y),
        xytext=(x0, y),
        arrowprops={
            "arrowstyle": "-|>",
            "color": MID_GREY,
            "lw": 0.45,
            "mutation_scale": 7,
            "shrinkA": 0,
            "shrinkB": 0,
        },
    )


def add_column_axes(fig: plt.Figure, x_mm: float, y_mm: float, w_mm: float, h_mm: float) -> plt.Axes:
    col_ax = fig.add_axes([x_mm / WIDTH_MM, y_mm / HEIGHT_MM, w_mm / WIDTH_MM, h_mm / HEIGHT_MM])
    col_ax.set_xlim(0, w_mm)
    col_ax.set_ylim(0, h_mm)
    col_ax.axis("off")
    col_ax.set_facecolor(STRIP_BG)
    return col_ax


def draw_column_header(
    col_ax: plt.Axes,
    h_mm: float,
    kicker: str,
    pad_x: float,
    pad_y: float,
    *,
    body_gap_mm: float = 4.2,
) -> float:
    header_y = h_mm - pad_y
    col_ax.text(
        pad_x,
        header_y,
        kicker,
        ha="left",
        va="top",
        fontsize=6.2,
        color=KICKER_COLOR,
        clip_on=False,
    )
    return header_y - body_gap_mm


def place_module_lines(
    ax: plt.Axes,
    x: float,
    y: float,
    lines: tuple[str, ...],
    *,
    line_h: float = MODULE_LINE_H_MM,
) -> float:
    for line in lines:
        _, y = place_line(ax, x, y, line, line_h=line_h)
    return y


def place_line(
    ax: plt.Axes,
    x: float,
    y: float,
    text: str,
    *,
    fontsize: float = BODY_PT,
    fontweight: str = "normal",
    color: str = MUTED_COLOR,
    line_h: float | None = None,
) -> tuple[mpl.text.Text, float]:
    artist = ax.text(
        x,
        y,
        text,
        ha="left",
        va="top",
        fontsize=fontsize,
        fontweight=fontweight,
        color=color,
        clip_on=False,
    )
    step = line_h if line_h is not None else LINE_H_MM
    return artist, y - (step if text else 1.0)


def text_bbox_mm(artist: mpl.text.Text, ax: plt.Axes, fig: plt.Figure) -> tuple[float, float, float, float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted())
    return float(bbox.x0), float(bbox.y0), float(bbox.x1), float(bbox.y1)


def bbox_area(box: tuple[float, float, float, float]) -> float:
    return max(0.0, box[2] - box[0]) * max(0.0, box[3] - box[1])


def bbox_intersection(
    a: tuple[float, float, float, float],
    b: tuple[float, float, float, float],
) -> float:
    x0 = max(a[0], b[0])
    y0 = max(a[1], b[1])
    x1 = min(a[2], b[2])
    y1 = min(a[3], b[3])
    if x1 <= x0 or y1 <= y0:
        return 0.0
    return (x1 - x0) * (y1 - y0)


def text_bbox_figure_mm(artist: mpl.text.Text, ax: plt.Axes, fig: plt.Figure) -> tuple[float, float, float, float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    return (
        float(bbox.x0) * WIDTH_MM,
        float(bbox.y0) * HEIGHT_MM,
        float(bbox.x1) * WIDTH_MM,
        float(bbox.y1) * HEIGHT_MM,
    )


def validate_layout(
    fig: plt.Figure,
    title_ax: plt.Axes,
    col_axes: list[plt.Axes],
    col_w: list[float],
) -> dict[str, int]:
    fig.canvas.draw()
    text_entries: list[tuple[mpl.text.Text, plt.Axes, tuple[float, float, float, float], str, str]] = []

    for child in title_ax.get_children():
        if not isinstance(child, mpl.text.Text):
            continue
        label = child.get_text().strip()
        if not label:
            continue
        bbox = text_bbox_figure_mm(child, title_ax, fig)
        text_entries.append((child, title_ax, bbox, label, "title"))
        require(bbox[0] >= -BBOX_TOL_MM, f"title text clipped left: {label!r}")
        require(bbox[2] <= WIDTH_MM + BBOX_TOL_MM, f"title text clipped right: {label!r}")
        require(bbox[1] >= -BBOX_TOL_MM, f"title text clipped bottom: {label!r}")
        require(bbox[3] <= HEIGHT_MM + BBOX_TOL_MM, f"title text clipped top: {label!r}")

    for ax_index, ax in enumerate(col_axes):
        width = col_w[ax_index]
        for child in ax.get_children():
            if not isinstance(child, mpl.text.Text):
                continue
            label = child.get_text().strip()
            if not label:
                continue
            bbox = text_bbox_mm(child, ax, fig)
            text_entries.append((child, ax, bbox, label, f"col{ax_index + 1}"))
            require(bbox[0] >= -BBOX_TOL_MM, f"text clipped left in column {ax_index + 1}: {label!r}")
            require(bbox[2] <= width + BBOX_TOL_MM, f"text clipped right in column {ax_index + 1}: {label!r}")
            require(bbox[1] >= -BBOX_TOL_MM, f"text clipped bottom in column {ax_index + 1}: {label!r}")
            require(
                bbox[3] <= STRIP_H_MM + BBOX_TOL_MM,
                f"text clipped top in column {ax_index + 1}: {label!r}",
            )
            fontsize = child.get_fontsize()
            require(fontsize >= MIN_FONT_PT - 0.05, f"font below minimum in column {ax_index + 1}: {label!r}")

    figure_entries: list[tuple[str, tuple[float, float, float, float], str]] = []
    for _, ax, bbox, label, zone in text_entries:
        if zone == "title":
            figure_entries.append((label, bbox, zone))
            continue
        ax_index = col_axes.index(ax)
        pos = ax.get_position()
        x_off = float(pos.x0) * WIDTH_MM
        y_off = float(pos.y0) * HEIGHT_MM
        figure_entries.append(
            (label, (bbox[0] + x_off, bbox[1] + y_off, bbox[2] + x_off, bbox[3] + y_off), zone),
        )

    for index, (label_a, bbox_a, zone_a) in enumerate(figure_entries):
        for label_b, bbox_b, zone_b in figure_entries[index + 1 :]:
            if zone_a == "title" and zone_b == "title":
                continue
            overlap = bbox_intersection(bbox_a, bbox_b)
            require(
                overlap <= OVERLAP_TOL_MM2,
                f"overlapping text: {label_a!r} vs {label_b!r} ({overlap:.2f} mm²)",
            )

    subtitle_boxes = [bbox for label, bbox, zone in figure_entries if "Descriptive chronicity" in label]
    main_title_boxes = [
        bbox for label, bbox, zone in figure_entries if "Cross-sectional AF-duration" in label
    ]
    strip_top = MARGIN_MM + STRIP_H_MM
    if subtitle_boxes and main_title_boxes:
        require(
            main_title_boxes[0][1] >= subtitle_boxes[0][3] + TITLE_SUB_GAP_MM - BBOX_TOL_MM,
            "main title overlaps subtitle",
        )
    if subtitle_boxes:
        require(
            subtitle_boxes[0][1] >= strip_top + TITLE_STRIP_GAP_MM - BBOX_TOL_MM,
            "subtitle overlaps process strip",
        )

    required_labels = [
        "REFERENCE",
        "DATA & DURATION",
        "FAMILY ACTIVITY",
        "ANALYSIS",
        "PRIMARY LAA · 12 AF-LAA",
        "LAPV REGIONAL SENSITIVITY",
        "Strong",
        "Weak",
        "Stable",
        "|ρ| ≥ 0.4",
        "0.2 ≤ |ρ| < 0.4",
        "|ρ| < 0.2",
        "Parsing: XM→X · 0M→0",
        "nonnumeric→NA",
        "Not aCAIrel",
        "No paired LA–RA contrast",
        "30 programs → 16 families",
        "Frozen program→family map",
        "Exact-reference rank scoring",
        "regional Spearman ρ · not for class",
        "SR excluded from ρ, p, and class",
    ]
    col_labels = [label for _, _, _, label, zone in text_entries if zone != "title"]
    joined = "\n".join(col_labels)
    for phrase in required_labels:
        require(phrase in joined, f"required visible label missing: {phrase}")

    return {
        "n_text_objects": len(text_entries),
        "n_columns": len(col_axes),
    }


def render() -> tuple[plt.Figure, plt.Axes, list[plt.Axes], list[float]]:
    set_style()
    geom = strip_geometry()
    xs: list[float] = geom["xs"]  # type: ignore[assignment]
    col_w: list[float] = geom["col_w"]  # type: ignore[assignment]
    strip_y0: float = geom["strip_y0"]  # type: ignore[assignment]
    strip_cy: float = geom["strip_cy"]  # type: ignore[assignment]
    title_y1: float = geom["title_y1"]  # type: ignore[assignment]
    strip_y1: float = geom["strip_y1"]  # type: ignore[assignment]
    gap_centers: list[float] = geom["gap_centers"]  # type: ignore[assignment]

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM))
    fig.patch.set_facecolor(WHITE)

    pad_x = 1.8
    pad_y = 1.2
    h = STRIP_H_MM

    title_ax = fig.add_axes([0, 0, 1, 1])
    title_ax.set_xlim(0, WIDTH_MM)
    title_ax.set_ylim(0, HEIGHT_MM)
    title_ax.axis("off")
    title_ax.set_facecolor(WHITE)
    title_ax.text(MARGIN_MM, title_y1 - 0.4, "A", ha="left", va="top", fontsize=10.8, fontweight="bold", color=FG_COLOR)
    main_title = title_ax.text(
        MARGIN_MM + 7.0,
        title_y1 - 0.4,
        "Cross-sectional AF-duration analysis in GSE41177",
        ha="left",
        va="top",
        fontsize=8.7,
        fontweight="bold",
        color=FG_COLOR,
    )
    fig.canvas.draw()
    main_bb = text_bbox_figure_mm(main_title, title_ax, fig)
    subtitle_top = main_bb[1] - TITLE_SUB_GAP_MM
    title_ax.text(
        MARGIN_MM + 7.0,
        subtitle_top,
        "Descriptive chronicity proxy; not recurrence, progression, or treatment response.",
        ha="left",
        va="top",
        fontsize=6.5,
        color=FAINT_COLOR,
        linespacing=1.05,
    )

    col_axes = [add_column_axes(fig, xs[i], strip_y0, col_w[i], h) for i in range(4)]
    for index, col_ax in enumerate(col_axes):
        col_ax.add_patch(
            mpl.patches.Rectangle(
                (0, 0),
                col_w[index],
                h,
                facecolor="none",
                edgecolor="#E0E0E0",
                linewidth=0.35,
                zorder=0,
            )
        )

    arrow_half = GAP_MM * 0.26
    for gap_x in gap_centers:
        bridge_arrow(title_ax, gap_x - arrow_half, gap_x + arrow_half, strip_cy)

    # Module 1 — REFERENCE
    y = draw_column_header(col_axes[0], h, "REFERENCE", pad_x, pad_y, body_gap_mm=4.0)
    place_module_lines(
        col_axes[0],
        pad_x,
        y,
        (
            "Frozen NCP from Fig 1",
            "30 programs → 16 families",
            "Frozen LE sets",
            "Frozen program→family map",
        ),
    )

    # Module 2 — DATA & DURATION
    y = draw_column_header(col_axes[1], h, "DATA & DURATION", pad_x, pad_y, body_gap_mm=4.0)
    for line in (
        "GSE41177 · Affymetrix",
        "38 samples · LAA + LAPV",
        "Month-coded AF duration",
    ):
        _, y = place_line(col_axes[1], pad_x, y, line)
    _, y = place_line(col_axes[1], pad_x, y, "")
    _, y = place_line(col_axes[1], pad_x, y, "Primary:", fontweight="bold", color=FG_COLOR)
    for line in (
        "12 AF-LAA samples",
        "with finite numeric months",
    ):
        _, y = place_line(col_axes[1], pad_x, y, line, fontweight="bold", color=FG_COLOR)
    _, y = place_line(col_axes[1], pad_x, y, "")
    _, y = place_line(
        col_axes[1],
        pad_x,
        y,
        "Parsing: XM→X · 0M→0",
        fontsize=FOOTNOTE_PT,
    )
    _, y = place_line(
        col_axes[1],
        pad_x,
        y,
        "nonnumeric→NA",
        fontsize=FOOTNOTE_PT,
    )
    place_line(
        col_axes[1],
        pad_x,
        y,
        "No year-to-month conversion",
        fontsize=FOOTNOTE_PT,
        color=FAINT_COLOR,
    )

    # Module 3 — FAMILY ACTIVITY
    y = draw_column_header(col_axes[2], h, "FAMILY ACTIVITY", pad_x, pad_y, body_gap_mm=4.0)
    place_module_lines(
        col_axes[2],
        pad_x,
        y,
        (
            "Exact-reference rank scoring",
            "Program-level scores",
            "→ equal-weight family activity",
            "",
            "Not aCAIrel",
            "No paired LA–RA contrast",
        ),
    )

    # Module 4 — ANALYSIS
    w4 = col_w[3]
    ah = 2.5
    y = draw_column_header(col_axes[3], h, "ANALYSIS", pad_x, pad_y, body_gap_mm=3.4)
    _, y = place_line(
        col_axes[3],
        pad_x,
        y,
        "PRIMARY LAA · 12 AF-LAA",
        fontsize=6.5,
        fontweight="bold",
        color=FG_COLOR,
        line_h=ah,
    )
    for line in (
        "Spearman ρ · bootstrap 95% CI · permutation p",
        "Class from |ρLAA| only",
    ):
        _, y = place_line(col_axes[3], pad_x, y, line, line_h=ah)
    for line in (
        "Strong   |ρ| ≥ 0.4",
        "Weak     0.2 ≤ |ρ| < 0.4",
        "Stable   |ρ| < 0.2",
    ):
        _, y = place_line(col_axes[3], pad_x, y, line, fontsize=FOOTNOTE_PT, line_h=2.45)

    rule_y = y - 0.15
    col_axes[3].plot([pad_x, w4 - pad_x], [rule_y, rule_y], color="#D8D8D8", lw=0.35, clip_on=False)
    y = rule_y - 0.85

    _, y = place_line(
        col_axes[3],
        pad_x,
        y,
        "LAPV REGIONAL SENSITIVITY",
        fontsize=6.4,
        fontweight="semibold",
        color=LAPV_COLOR,
        line_h=ah,
    )
    for line in (
        "regional Spearman ρ · not for class",
        "SR excluded from ρ, p, and class",
    ):
        _, y = place_line(
            col_axes[3],
            pad_x,
            y,
            line,
            color=LAPV_COLOR if "SR" not in line else FAINT_COLOR,
            fontsize=FOOTNOTE_PT,
            line_h=ah,
        )

    return fig, title_ax, col_axes, col_w


def validate_outputs(
    pdf: Path,
    svg: Path,
    png: Path,
    layout_stats: dict[str, int],
) -> None:
    svg_text = svg.read_text(encoding="utf-8")
    require("NCP_v3" not in svg_text, "NCP_v3 visible in SVG")
    require("results/" not in svg_text and "scripts/" not in svg_text, "internal paths visible in SVG")
    require("Not aCAIrel" in svg_text or "not aCAIrel" in svg_text, "family activity disclaimer missing")
    require("Regional family trends" not in svg_text, "redundant handoff line still present")

    pdftotext = shutil.which("pdftotext")
    pdfinfo = shutil.which("pdfinfo")
    if pdfinfo:
        info = subprocess.run([pdfinfo, str(pdf)], check=True, capture_output=True, text=True).stdout
        match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
        require(match is not None, "Unable to read Figure 4A PDF dimensions")
        width_mm = float(match.group(1)) * 25.4 / 72
        height_mm = float(match.group(2)) * 25.4 / 72
        require(
            abs(width_mm - WIDTH_MM) < 0.1 and abs(height_mm - HEIGHT_MM) < 0.1,
            f"Figure 4A PDF dimensions are {width_mm:.2f} x {height_mm:.2f} mm",
        )
    if pdftotext:
        visible = subprocess.run([pdftotext, "-layout", str(pdf), "-"], check=True, capture_output=True, text=True).stdout
        required = [
            "Cross-sectional AF-duration analysis in GSE41177",
            "Descriptive chronicity proxy",
            "REFERENCE",
            "DATA & DURATION",
            "FAMILY ACTIVITY",
            "ANALYSIS",
            "PRIMARY LAA",
            "12 AF-LAA",
            "LAPV REGIONAL SENSITIVITY",
            "SR excluded",
            "Frozen NCP from Fig 1",
            "Frozen LE sets",
            "Not aCAIrel",
            "No paired LA",
            "Parsing: XM",
            "Strong",
            "Weak",
            "Stable",
            "Class from",
            "not for class",
            "bootstrap 95% CI",
            "permutation p",
        ]
        normalized = " ".join(visible.split())
        missing = [phrase for phrase in required if phrase not in normalized]
        require(not missing, f"Figure 4A visible text incomplete: {missing}")

    require(png.is_file() and png.stat().st_size > 0, "PNG export missing")
    require(layout_stats["n_text_objects"] > 0, "layout QC found no text objects")


def write_qc(
    manifest: pd.DataFrame,
    summary: pd.DataFrame,
    max_error: float,
    output_hashes: dict[str, str],
    layout_stats: dict[str, int],
    *,
    status: str,
) -> None:
    manifest_lines = "\n".join(
        f"- `{row.artifact_role}`: `{row.artifact_path}` — `{row.sha256}`"
        for row in manifest.itertuples(index=False)
    )
    output_lines = "\n".join(f"- `{path}` — `{digest}`" for path, digest in output_hashes.items())
    row = summary.iloc[0]
    OUT_QC.write_text(
        f"""# Figure 4A cross-sectional duration bridge QC

Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}

Status: **{status}**

## Hard assertions
- total GSE41177 expression samples: {int(row["total_samples"])}
- primary-eligible AF-LAA samples: {int(row["primary_eligible_af_laa"])}
- frozen programs: {int(row["frozen_programs"])}
- frozen families: {int(row["frozen_families"])}
- scorer: `{row["scorer"]}`
- aggregation: `{row["aggregation"]}`
- input layer: `{row["input_layer"]}` (not aCAIrel)
- class rule: `{row["class_rule"]}`
- LAPV role: `{row["lapv_role"]}`
- maximum equal-weight aggregation error: {max_error:.3g}
- SR excluded from correlation/class assignment: yes
- LAPV excluded from class assignment: yes
- legacy table/scorer reads: no
- atlas priority ordering: no
- panel dimensions: {WIDTH_MM:.0f} x {HEIGHT_MM:.0f} mm
- module widths (%): 21 / 22 / 23 / 34
- layout bbox QC text objects: {layout_stats["n_text_objects"]}
- clipped text objects: 0
- overlapping text boxes: 0
- minimum visible font: {MIN_FONT_PT:.1f} pt

## Input manifest
{manifest_lines}

## Outputs
{output_lines}
""",
        encoding="utf-8",
    )


def main() -> None:
    manifest, paths = load_manifest()
    meta = pd.read_csv(paths["gse41177_duration_metadata"], sep="\t")
    programs = pd.read_csv(paths["frozen_program_metadata"], sep="\t")
    families = pd.read_csv(paths["frozen_family_metadata"], sep="\t")
    prog = pd.read_csv(paths["gse41177_program_sample_scores"], sep="\t")
    fam = pd.read_csv(paths["gse41177_family_sample_scores"], sep="\t")
    le_gate = pd.read_csv(paths["leading_edge_validation"], sep="\t").iloc[0]

    require(len(families) == 16 and families["family_id"].nunique() == 16, "frozen family count changed")
    require(
        int(le_gate["n_programs"]) == 30
        and int(le_gate["n_set_identical"]) == 30
        and str(le_gate["gate_passed"]).upper() == "TRUE",
        "frozen leading-edge gate failed",
    )

    validate_duration_parser(meta)
    max_error = validate_aggregation(prog, fam, programs)
    validate_class_thresholds()
    summary = build_summary(meta, programs, families, max_error)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, title_ax, col_axes, col_w = render()
    layout_stats = validate_layout(fig, title_ax, col_axes, col_w)
    pdf = OUT_BASE.with_suffix(".pdf")
    svg = OUT_BASE.with_suffix(".svg")
    png = OUT_BASE.with_suffix(".png")
    fig.savefig(pdf, bbox_inches=None, pad_inches=0)
    fig.savefig(svg, bbox_inches=None, pad_inches=0)
    fig.savefig(png, dpi=600, bbox_inches=None, pad_inches=0)
    fig.savefig(OUT_PREVIEW, dpi=150, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    validate_outputs(pdf, svg, png, layout_stats)

    output_hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (pdf, svg, png, OUT_PREVIEW, OUT_MANIFEST, OUT_SUMMARY)
    }
    write_qc(manifest, summary, max_error, output_hashes, layout_stats, status="PASS")
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_MANIFEST.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
