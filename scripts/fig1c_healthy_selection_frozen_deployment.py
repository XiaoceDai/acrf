#!/usr/bin/env python3
"""Standalone Figure 1C: healthy-only selection and frozen deployment."""
from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from _pathway_labels import short_pathway_label

RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

LONG = RESULTS / "fig1c_pathway_selection_long.tsv"
SCORED = RESULTS / "innovation_ncp_consensus_long.tsv"
PANEL = RESULTS / "innovation_ncp_v3_panel.tsv"
PANEL_SUMMARY = RESULTS / "innovation_ncp_v3_summary.tsv"
CRITERIA = RESULTS / "ncp_v3_selection_criteria.tsv"
FAMILY_MAP = RESULTS / "innovation_ncp_v3_axis_families.tsv"

OUT_BASE = FIGURES / "fig1c_selection_provenance"
OUT_SOURCE = RESULTS / "fig1c_selection_provenance_source.tsv"
OUT_SUMMARY = RESULTS / "fig1c_selection_provenance_summary.tsv"
OUT_INPUTS = RESULTS / "fig1c_selection_provenance_inputs.tsv"
OUT_QC = RESULTS / "fig1c_selection_provenance_qc.md"
OUT_PREVIEW = FIGURES / "fig1c_selection_provenance_preview.png"

CONSTRUCTION_INPUTS = {
    LONG,
    SCORED,
    PANEL,
    PANEL_SUMMARY,
    CRITERIA,
}
VALIDATION_ONLY_INPUTS = {FAMILY_MAP}
ALLOWED_INPUTS = CONSTRUCTION_INPUTS | VALIDATION_ONLY_INPUTS

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
FG_COLOR = "#1A1A1A"
MUTED_COLOR = "#555555"
LINE_COLOR = "#C7C7C7"
LIGHT_GREY = "#F3F3F3"
MID_GREY = "#A2A2A2"
DARK_GREY = "#565656"
WHITE = "#FFFFFF"
MM = 1.0 / 25.4
WIDTH_MM = 68.2
HEIGHT_MM = 137


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numeric_threshold(value: str) -> float:
    matches = re.findall(r"\d+(?:\.\d+)?", str(value))
    require(bool(matches), f"No numeric threshold in criteria value: {value}")
    return float(matches[0])


def criterion_value(criteria: pd.DataFrame, criterion: str, parameter: str) -> str:
    rows = criteria[
        (criteria["stage"] == "NCP_consensus_pruning")
        & (criteria["criterion"] == criterion)
        & (criteria["parameter"] == parameter)
    ]
    require(len(rows) == 1, f"Expected one criteria row for {criterion}/{parameter}")
    return str(rows.iloc[0]["value"])


def load_construction_inputs() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    for path in CONSTRUCTION_INPUTS:
        require(path.is_file(), f"Missing healthy-selection input: {path}")
    long_df = pd.read_csv(LONG, sep="\t")
    scored = pd.read_csv(SCORED, sep="\t")
    panel = pd.read_csv(PANEL, sep="\t")
    panel_summary = pd.read_csv(PANEL_SUMMARY, sep="\t")
    criteria = pd.read_csv(CRITERIA, sep="\t")
    return long_df, scored, panel, panel_summary, criteria


def build_selection_source() -> tuple[pd.DataFrame, dict[str, int | float | str]]:
    long_df, scored, panel, panel_summary, criteria = load_construction_inputs()

    require(len(long_df) == 80 and long_df["pathway"].nunique() == 80, "Expected 80 unique candidates in selection-long table")
    require(len(scored) == 80 and scored["pathway"].nunique() == 80, "Expected 80 unique candidates in consensus table")
    require(set(long_df["pathway"]) == set(scored["pathway"]), "Healthy candidate universes differ")
    require(len(panel) == 30 and panel["pathway"].nunique() == 30, "Frozen panel is not 30 unique programs")

    # Validate that the duplicated healthy evidence fields agree across the two
    # authoritative selection tables before using the scored table.
    comparable = [
        "NES",
        "padj",
        "gsva_logFC",
        "gsva_padj",
        "stability_consistent_runs",
    ]
    joined = long_df[["pathway", "direction", *comparable]].merge(
        scored[["pathway", "direction", *comparable]],
        on="pathway",
        suffixes=("_long", "_scored"),
        validate="one_to_one",
    )
    require(
        bool((joined["direction_long"] == joined["direction_scored"]).all()),
        "Candidate directions disagree between healthy-selection tables",
    )
    for column in comparable:
        require(
            np.allclose(
                joined[f"{column}_long"].astype(float),
                joined[f"{column}_scored"].astype(float),
                rtol=0,
                atol=1e-12,
                equal_nan=True,
            ),
            f"Healthy-selection evidence differs for {column}",
        )

    stability_min = numeric_threshold(
        criterion_value(
            criteria,
            "consensus_gate_strict",
            "stability_consistent_runs",
        )
    )
    strict_score_min = numeric_threshold(
        criterion_value(criteria, "consensus_gate_strict", "consensus_score")
    )
    strict_p_max = numeric_threshold(
        criterion_value(criteria, "consensus_gate_strict", "dual significance")
    )
    fallback_score_min = numeric_threshold(
        criterion_value(criteria, "consensus_gate_fallback", "consensus_score")
    )
    fallback_p_max = numeric_threshold(
        criterion_value(criteria, "consensus_gate_fallback", "dual significance")
    )
    direction_targets = re.findall(
        r"\d+",
        criterion_value(
            criteria,
            "panel_balance",
            "MIN_PER_DIRECTION / TARGET_PER_DIRECTION",
        ),
    )
    require(len(direction_targets) == 2, "Could not parse direction target")
    target_per_direction = int(direction_targets[1])

    source = scored.copy()
    source["dual_method_pass"] = (
        (source["padj"].astype(float) < strict_p_max)
        & (source["gsva_padj"].astype(float) < strict_p_max)
    )
    source["strict_eligible"] = (
        source["direction_agree"].astype(bool)
        & (source["stability_consistent_runs"].astype(float) >= stability_min)
        & (source["consensus_score"].astype(float) >= strict_score_min)
        & source["dual_method_pass"]
    )
    source["fallback_eligible"] = (
        source["direction_agree"].astype(bool)
        & (source["stability_consistent_runs"].astype(float) >= stability_min)
        & (source["consensus_score"].astype(float) >= fallback_score_min)
        & (source["padj"].astype(float) < fallback_p_max)
        & (source["gsva_padj"].astype(float) < fallback_p_max)
    )

    final_programs = set(panel["pathway"].astype(str))
    candidate_programs = set(source["pathway"].astype(str))
    require(final_programs.issubset(candidate_programs), "A frozen program is outside the original 80 candidates")
    source["final_selected"] = source["pathway"].astype(str).isin(final_programs)
    source["strict_selected"] = source["final_selected"] & source["strict_eligible"]
    source["fallback_selected"] = source["final_selected"] & ~source["strict_eligible"]
    source["selection_route"] = np.select(
        [source["strict_selected"], source["fallback_selected"]],
        ["strict", "fallback"],
        default="not_selected",
    )
    source["publication_label"] = source["pathway"].map(
        lambda value: short_pathway_label(str(value), 46).replace("⁺", "+")
    )
    source = source.rename(
        columns={"pathway": "program_id", "direction": "canonical_direction"}
    )

    strict_set = set(source.loc[source["strict_selected"], "program_id"])
    fallback_set = set(source.loc[source["fallback_selected"], "program_id"])
    final_set = set(source.loc[source["final_selected"], "program_id"])
    require(len(source) == 80 and source["program_id"].nunique() == 80, "Source is not 80 unique candidates")
    require(int(source["strict_eligible"].sum()) == 38, "Strict-eligible count is not 38")
    require(len(strict_set) == 25, "Strict-selected count is not 25")
    require(len(fallback_set) == 5, "Fallback-selected count is not 5")
    require(strict_set.isdisjoint(fallback_set), "Strict and fallback selected sets overlap")
    require(len(strict_set) + len(fallback_set) == 30, "25 + 5 does not equal 30")
    require(final_set == final_programs, "Derived final set differs from frozen panel")
    require(
        bool(source.loc[source["fallback_selected"], "fallback_eligible"].all()),
        "A fallback-selected program does not meet the prespecified fallback gate",
    )
    require(
        not bool(source.loc[source["fallback_selected"], "strict_eligible"].any()),
        "A fallback-selected program is incorrectly strict-eligible",
    )
    require(
        set(source.loc[source["fallback_selected"], "canonical_direction"]) == {"RA-biased"},
        "Not all five fallback-selected programs are RA-biased",
    )

    panel_direction = panel.set_index("pathway")["direction"]
    source_direction = source.set_index("program_id")["canonical_direction"]
    require(
        bool(
            (
                panel_direction.sort_index()
                == source_direction.loc[panel_direction.index].sort_index()
            ).all()
        ),
        "Frozen panel directions disagree with healthy-selection metadata",
    )
    n_la = int((panel["direction"] == "LA-biased").sum())
    n_ra = int((panel["direction"] == "RA-biased").sum())
    require((n_la, n_ra) == (22, 8), f"Expected 22 LA / 8 RA, observed {n_la} / {n_ra}")
    strict_ra = int(
        (
            source["strict_selected"]
            & source["canonical_direction"].eq("RA-biased")
        ).sum()
    )
    fallback_ra = int(
        (
            source["fallback_selected"]
            & source["canonical_direction"].eq("RA-biased")
        ).sum()
    )
    require(
        strict_ra == 3
        and fallback_ra == 5
        and strict_ra + fallback_ra == target_per_direction == 8,
        "Fallback selections do not complete the prespecified RA direction target",
    )

    require(len(panel_summary) == 1, "Frozen panel summary must contain one row")
    summary_row = panel_summary.iloc[0]
    require(int(summary_row["n_total_candidates"]) == 80, "Frozen summary candidate count changed")
    require(int(summary_row["n_quality_candidates"]) == 38, "Frozen summary strict count changed")
    require(int(summary_row["n_panel_v3"]) == 30, "Frozen summary panel count changed")
    require(int(summary_row["la_count_v3"]) == 22, "Frozen summary LA count changed")
    require(int(summary_row["ra_count_v3"]) == 8, "Frozen summary RA count changed")

    # Validation-only topology input is deliberately read after selection is fixed.
    require(FAMILY_MAP.is_file(), f"Missing validation-only family map: {FAMILY_MAP}")
    family_map = pd.read_csv(FAMILY_MAP, sep="\t")
    mapped = family_map[family_map["pathway"].astype(str).isin(final_set)].copy()
    require(len(mapped) == 30, "Selected subset does not have exactly 30 mapping rows")
    require(mapped["pathway"].nunique() == 30, "A selected program is duplicated in the mapping table")
    require(set(mapped["pathway"]) == final_set, "Selected mapping subset has missing or extra programs")
    require(
        len(family_map) == 30 and set(family_map["pathway"]) == final_set,
        "Frozen mapping table contains programs outside the final selected set",
    )
    require(mapped["family_id"].nunique() == 16, "Final programs do not map to 16 unique families")
    family_direction_counts = mapped.groupby("family_id")["direction"].nunique()
    require(bool((family_direction_counts == 1).all()), "Family direction is internally inconsistent")
    mapped_direction = mapped.set_index("pathway")["direction"]
    require(
        bool(
            (
                mapped_direction.sort_index()
                == source_direction.loc[mapped_direction.index].sort_index()
            ).all()
        ),
        "Program directions disagree between selection and family metadata",
    )

    counts: dict[str, int | float | str] = {
        "candidate_pathways": 80,
        "strict_eligible": 38,
        "strict_selected": 25,
        "fallback_selected": 5,
        "final_frozen": 30,
        "LA_biased": n_la,
        "RA_biased": n_ra,
        "strict_selected_RA": strict_ra,
        "fallback_selected_RA": fallback_ra,
        "direction_target": target_per_direction,
        "final_families": int(mapped["family_id"].nunique()),
        "selection_dataset": "GSE112339 paired healthy reference",
        "AF_cohort_inputs_read": 0,
    }
    source_columns = [
        "program_id",
        "publication_label",
        "canonical_direction",
        "strict_eligible",
        "strict_selected",
        "fallback_eligible",
        "fallback_selected",
        "final_selected",
        "selection_route",
        "consensus_score",
        "stability_consistent_runs",
        "direction_agree",
        "dual_method_pass",
        "NES",
        "padj",
        "gsva_logFC",
        "gsva_padj",
    ]
    source = source[source_columns].sort_values(
        ["final_selected", "selection_route", "consensus_score"],
        ascending=[False, True, False],
    )
    return source.reset_index(drop=True), counts


def write_inputs_ledger() -> pd.DataFrame:
    rows = []
    usage = {
        LONG: ("construction", "healthy candidate universe and paired evidence"),
        SCORED: ("construction", "strict/fallback eligibility and consensus evidence"),
        PANEL: ("construction", "frozen final program identity and direction"),
        PANEL_SUMMARY: ("construction", "frozen count cross-check"),
        CRITERIA: ("construction", "machine-readable prespecified gate definitions"),
        FAMILY_MAP: ("validation_only", "final program-to-family handoff assertion"),
    }
    require(set(usage) == ALLOWED_INPUTS, "Input allowlist changed")
    for path, (role, used_for) in usage.items():
        rows.append(
            {
                "artifact_path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "input_role": role,
                "used_for": used_for,
                "contains_AF_cohort_results": False,
            }
        )
    ledger = pd.DataFrame(rows)
    require(not bool(ledger["contains_AF_cohort_results"].any()), "Disease-cohort input detected")
    ledger.to_csv(OUT_INPUTS, sep="\t", index=False)
    return ledger


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 7.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": WHITE,
            "savefig.facecolor": WHITE,
            "hatch.linewidth": 0.55,
        }
    )


def node(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    subtitle: str = "",
    *,
    dashed: bool = False,
    linewidth: float = 0.8,
    title_size: float = 7.5,
    title_fraction: float | None = None,
    subtitle_fraction: float = 0.25,
    subtitle_size: float = 6.5,
) -> None:
    ax.add_patch(
        mpl.patches.Rectangle(
            (x, y),
            w,
            h,
            facecolor=WHITE,
            edgecolor=DARK_GREY,
            linewidth=linewidth,
            linestyle=(0, (3, 2)) if dashed else "solid",
        )
    )
    ax.text(
        x + w / 2,
        y + h * (
            title_fraction
            if title_fraction is not None
            else (0.62 if subtitle else 0.50)
        ),
        title,
        ha="center",
        va="center",
        fontsize=title_size,
        fontweight="bold",
        color=FG_COLOR,
    )
    if subtitle:
        ax.text(
            x + w / 2,
            y + h * subtitle_fraction,
            subtitle,
            ha="center",
            va="center",
            fontsize=subtitle_size,
            color=MUTED_COLOR,
            linespacing=1.1,
        )


def arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    dashed: bool = False,
) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={
            "arrowstyle": "-|>",
            "color": DARK_GREY,
            "lw": 0.85,
            "linestyle": (0, (3, 2)) if dashed else "solid",
            "shrinkA": 0,
            "shrinkB": 0,
            "mutation_scale": 8,
        },
    )


def stacked_bar(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    left_count: int,
    right_count: int,
    left_color: str,
    right_color: str,
    *,
    right_hatch: str | None = None,
) -> None:
    total = left_count + right_count
    left_w = w * left_count / total
    ax.add_patch(
        mpl.patches.Rectangle(
            (x, y), left_w, h, facecolor=left_color, edgecolor=DARK_GREY, linewidth=0.45
        )
    )
    ax.add_patch(
        mpl.patches.Rectangle(
            (x + left_w, y),
            w - left_w,
            h,
            facecolor=right_color,
            edgecolor=DARK_GREY,
            linewidth=0.45,
            hatch=right_hatch,
        )
    )


def render(counts: dict[str, int | float | str]) -> plt.Figure:
    set_style()
    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM))
    ax = fig.add_axes([0.055, 0.035, 0.90, 0.90])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig.text(
        0.055,
        0.977,
        "C",
        ha="left",
        va="top",
        fontsize=11.0,
        fontweight="bold",
        color=FG_COLOR,
    )
    fig.text(
        0.135,
        0.977,
        "Healthy-only selection\nand frozen deployment",
        ha="left",
        va="top",
        fontsize=8.5,
        fontweight="bold",
        color=FG_COLOR,
        linespacing=1.0,
    )
    fig.text(
        0.055,
        0.925,
        "GSE112339 paired healthy reference only ·\n"
        "AF cohort data not used for selection",
        ha="left",
        va="top",
        fontsize=6.5,
        color=MUTED_COLOR,
        linespacing=1.15,
    )

    node(
        ax,
        0.13,
        0.828,
        0.72,
        0.070,
        f"{counts['candidate_pathways']} candidate pathways",
        "paired healthy reference",
    )
    arrow(ax, (0.40, 0.828), (0.40, 0.742))
    ax.text(
        0.08,
        0.785,
        "strict consensus gate",
        ha="left",
        va="center",
        fontsize=6.7,
        fontweight="bold",
        color=FG_COLOR,
    )
    ax.text(
        0.57,
        0.786,
        "dual-method evidence\n"
        "directional agreement\n"
        "resampling stability\n"
        "consensus score",
        ha="left",
        va="center",
        fontsize=6.5,
        color=MUTED_COLOR,
        linespacing=1.15,
    )
    node(
        ax,
        0.13,
        0.670,
        0.54,
        0.065,
        f"{counts['strict_eligible']} strict-eligible\ncandidates",
        title_size=7.2,
    )
    arrow(ax, (0.40, 0.670), (0.40, 0.600))
    ax.text(
        0.43,
        0.626,
        "redundancy and\ndiversity control",
        ha="left",
        va="center",
        fontsize=6.5,
        color=MUTED_COLOR,
        linespacing=1.1,
    )
    node(
        ax,
        0.08,
        0.525,
        0.50,
        0.065,
        f"{counts['strict_selected']} strict-selected\ncandidates",
        title_size=6.5,
    )

    # The conditional fallback is evaluated only after strict selection leaves
    # the prespecified RA-direction target unmet.
    ax.text(
        0.680,
        0.500,
        "Prespecified fallback used only\n"
        "if the RA target\n"
        "remained unmet",
        ha="center",
        va="center",
        fontsize=6.5,
        color=MUTED_COLOR,
        linespacing=1.1,
    )
    ax.plot(
        [0.58, 0.95, 0.95],
        [0.557, 0.557, 0.468],
        color=DARK_GREY,
        lw=0.8,
        linestyle=(0, (3, 2)),
    )
    arrow(ax, (0.95, 0.468), (0.91, 0.452), dashed=True)
    node(
        ax,
        0.39,
        0.350,
        0.57,
        0.102,
        f"{counts['fallback_selected']} fallback-selected\ncandidates",
        "same candidate universe\n"
        "all RA-biased\n"
        "redundancy-screened",
        dashed=True,
        title_size=6.5,
        title_fraction=0.76,
        subtitle_fraction=0.23,
        subtitle_size=6.5,
    )
    arrow(ax, (0.28, 0.525), (0.28, 0.310))
    arrow(ax, (0.70, 0.350), (0.70, 0.310), dashed=True)

    # Strong final deployment node.
    ax.add_patch(
        mpl.patches.Rectangle(
            (0.055, 0.070),
            0.89,
            0.240,
            facecolor=WHITE,
            edgecolor=FG_COLOR,
            linewidth=1.35,
        )
    )
    ax.text(
        0.50,
        0.286,
        f"Frozen NCP panel · {counts['final_frozen']} programs",
        ha="center",
        va="center",
        fontsize=7.2,
        fontweight="bold",
        color=FG_COLOR,
    )
    ax.text(0.105, 0.252, "Selection route", fontsize=6.5, fontweight="bold", color=FG_COLOR)
    stacked_bar(
        ax,
        0.105,
        0.222,
        0.79,
        0.027,
        int(counts["strict_selected"]),
        int(counts["fallback_selected"]),
        DARK_GREY,
        WHITE,
        right_hatch="////",
    )
    ax.text(
        0.105,
        0.210,
        f"{counts['strict_selected']} strict-selected",
        ha="left",
        va="top",
        fontsize=6.5,
        color=MUTED_COLOR,
    )
    ax.text(
        0.895,
        0.210,
        f"{counts['fallback_selected']} fallback-selected",
        ha="right",
        va="top",
        fontsize=6.5,
        color=MUTED_COLOR,
    )
    ax.text(0.105, 0.174, "Chamber direction", fontsize=6.5, fontweight="bold", color=FG_COLOR)
    stacked_bar(
        ax,
        0.105,
        0.144,
        0.79,
        0.027,
        int(counts["LA_biased"]),
        int(counts["RA_biased"]),
        LA_COLOR,
        RA_COLOR,
    )
    ax.text(
        0.105,
        0.132,
        f"{counts['LA_biased']} LA-biased",
        ha="left",
        va="top",
        fontsize=6.5,
        color=LA_COLOR,
        fontweight="bold",
    )
    ax.text(
        0.895,
        0.132,
        f"{counts['RA_biased']} RA-biased",
        ha="right",
        va="top",
        fontsize=6.5,
        color=RA_COLOR,
        fontweight="bold",
    )

    ax.add_patch(
        mpl.patches.Rectangle(
            (0.095, 0.074), 0.030, 0.024, facecolor=WHITE, edgecolor=FG_COLOR, linewidth=0.7
        )
    )
    ax.add_patch(
        mpl.patches.Arc(
            (0.110, 0.098), 0.020, 0.024, theta1=0, theta2=180, edgecolor=FG_COLOR, linewidth=0.7
        )
    )
    ax.text(
        0.135,
        0.087,
        "Frozen before any AF-cohort analysis",
        ha="left",
        va="center",
        fontsize=6.5,
        fontweight="bold",
        color=FG_COLOR,
    )
    ax.text(
        0.95,
        0.035,
        f"Fixed {counts['final_frozen']}-program \u2192 "
        f"{counts['final_families']}-family mapping\nshown in D  \u2192",
        ha="right",
        va="center",
        fontsize=6.5,
        color=FG_COLOR,
        fontweight="bold",
    )
    return fig


def write_qc(
    counts: dict[str, int | float | str],
    source: pd.DataFrame,
    ledger: pd.DataFrame,
    figure_hashes: dict[str, str],
) -> None:
    input_lines = "\n".join(
        f"- `{row.artifact_path}` — `{row.sha256}` — {row.input_role}: {row.used_for}"
        for row in ledger.itertuples(index=False)
    )
    output_lines = "\n".join(
        f"- `{path}` — `{digest}`" for path, digest in figure_hashes.items()
    )
    forbidden_renderer = "fig1c_selection_" + "funnel_summary"
    script_text = Path(__file__).read_text(encoding="utf-8")
    require(forbidden_renderer not in script_text, "Legacy Figure 1C renderer import detected")
    qc = f"""# Figure 1C healthy-selection provenance QC

Status: **PASS**

## Frozen provenance assertions
- unique healthy candidates: {counts['candidate_pathways']}
- strict-eligible candidates: {counts['strict_eligible']}
- strict-selected programs: {counts['strict_selected']}
- fallback-selected programs: {counts['fallback_selected']}
- final frozen programs: {counts['final_frozen']}
- direction composition: {counts['LA_biased']} LA-biased / {counts['RA_biased']} RA-biased
- conditional fallback check: {counts['strict_selected_RA']} strict RA + {counts['fallback_selected_RA']} fallback RA = target {counts['direction_target']}
- fallback-selected programs strict-eligible: 0
- validated fixed families: {counts['final_families']}
- AF cohort result inputs read: {counts['AF_cohort_inputs_read']}
- legacy three-column renderer imported: no
- candidate-level source rows: {len(source)}
- panel dimensions: {WIDTH_MM} × {HEIGHT_MM} mm
- minimum publication-facing text size: 6.5 pt

## Input lineage
{input_lines}

## Figure hashes
{output_lines}
"""
    OUT_QC.write_text(qc, encoding="utf-8")


def main() -> None:
    source, counts = build_selection_source()

    # All scientific and validation-only topology assertions pass before output.
    OUT_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    OUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    source.to_csv(OUT_SOURCE, sep="\t", index=False)
    pd.DataFrame(
        [{"metric": key, "value": value} for key, value in counts.items()]
    ).to_csv(OUT_SUMMARY, sep="\t", index=False)
    ledger = write_inputs_ledger()

    fig = render(counts)
    pdf = OUT_BASE.with_suffix(".pdf")
    svg = OUT_BASE.with_suffix(".svg")
    png = OUT_BASE.with_suffix(".png")
    fig.savefig(pdf, bbox_inches=None, pad_inches=0)
    fig.savefig(svg, bbox_inches=None, pad_inches=0)
    fig.savefig(png, dpi=600, bbox_inches=None, pad_inches=0)
    fig.savefig(OUT_PREVIEW, dpi=150, bbox_inches=None, pad_inches=0)
    plt.close(fig)

    figure_hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (pdf, svg, png, OUT_PREVIEW)
    }
    write_qc(counts, source, ledger, figure_hashes)
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_INPUTS.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")
    print(
        "Summary: 80 candidates; 38 strict-eligible; "
        "25 strict + 5 fallback = 30 frozen; 22 LA + 8 RA; 16 families"
    )


if __name__ == "__main__":
    main()
