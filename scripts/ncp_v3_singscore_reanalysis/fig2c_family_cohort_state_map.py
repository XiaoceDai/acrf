#!/usr/bin/env python3
"""Final Figure 2C: descriptive family-by-cohort rewiring states."""
from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fig2b_family_donor_heatmap as fig2b


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"

OUT_SOURCE = RESULTS / "fig2c_family_cohort_states_source.tsv"
OUT_SUMMARY = RESULTS / "fig2c_family_cohort_states_summary.tsv"
OUT_QC = RESULTS / "fig2c_family_cohort_states_qc.md"
OUT_GEOMETRY = RESULTS / "figure2_family_row_geometry.tsv"
OUT_GEOMETRY_COMPOSITE = RESULTS / "figure2_family_row_geometry_composite.tsv"
OUT_BASE = FIGURES / "fig2c_family_cohort_state_map"
OUT_PREVIEW = FIGURES / "fig2c_family_cohort_state_map_preview.png"
OUT_BASE_COMPOSITE = FIGURES / "fig2c_family_cohort_state_map_composite"
OUT_PREVIEW_COMPOSITE = FIGURES / "fig2c_family_cohort_state_map_preview_composite.png"
FIG2B_SOURCE = RESULTS / "fig2b_family_donor_source.tsv"

COHORT_ORDER = ["GSE128188", "GSE79768", "GSE115574"]
EXPECTED_COUNTS = {
    "GSE128188": {"SR": 5, "AF": 5},
    "GSE79768": {"SR": 6, "AF": 7},
    "GSE115574": {"SR": 15, "AF": 13},
}
EXPECTED_FAMILY_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]
STATE_ORDER = ["Amplified", "Preserved", "Attenuated", "Reversed"]
STATE_ABBREVIATIONS = {
    "Amplified": "Amp",
    "Preserved": "Pres",
    "Attenuated": "Att",
    "Reversed": "Rev",
}
STATE_COLORS = {
    "Amplified": "#D58A45",
    "Preserved": "#659478",
    "Attenuated": "#D8C27E",
    "Reversed": "#B5A3CC",
}
LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
FG = "#202020"
MUTED = "#606060"
LIGHT_LINE = "#EEEEEE"
MARKER_COLOR = "#6F6F6F"
WEDGE_FRACTION = 0.11

SOURCE_COLUMNS = [
    "family_id",
    "family_label",
    "family_direction",
    "canonical_row_order",
    "cohort",
    "n_SR",
    "n_AF",
    "mean_SR_aCAIrel",
    "mean_AF_aCAIrel",
    "sr_reference",
    "floor_applied",
    "af_sr_ratio",
    "state",
    "state_order",
]
GEOMETRY_COLUMNS = [
    "family_id",
    "canonical_order",
    "row_center",
    "row_height",
    "direction_block",
    "LA_RA_separator_y",
    "matrix_body_top",
    "matrix_body_bottom",
]

MM = 1 / 25.4
PANEL_WIDTH_MM = 88.0

# Shared composite row contract with Figure 2D (side-by-side collage only).
COMPOSITE_PANEL_HEIGHT_MM = 134.0
COMPOSITE_TABLE_TOP_MM = 27.0
COMPOSITE_ROW_HEIGHT_MM = 5.60


@dataclass(frozen=True)
class PanelLayout:
    """Physical millimetre contract for a Figure 2C render variant."""

    panel_height_mm: float
    table_top_mm: float
    row_height_mm: float
    header_top_mm: float
    label_header_y_mm: float
    legend_top_mm: float
    floor_top_mm: float

    @property
    def table_bottom_mm(self) -> float:
        return self.table_top_mm + 16 * self.row_height_mm

    @property
    def separator_mm(self) -> float:
        return self.table_top_mm + 11 * self.row_height_mm

    @property
    def header_height_mm(self) -> float:
        return self.table_top_mm - self.header_top_mm


STANDALONE_LAYOUT = PanelLayout(
    panel_height_mm=102.0,
    table_top_mm=29.0,
    row_height_mm=3.35,
    header_top_mm=21.0,
    label_header_y_mm=25.8,
    legend_top_mm=86.5,
    floor_top_mm=94.0,
)
COMPOSITE_LAYOUT = PanelLayout(
    panel_height_mm=COMPOSITE_PANEL_HEIGHT_MM,
    table_top_mm=COMPOSITE_TABLE_TOP_MM,
    row_height_mm=COMPOSITE_ROW_HEIGHT_MM,
    header_top_mm=19.0,
    label_header_y_mm=23.8,
    legend_top_mm=119.0,
    floor_top_mm=126.5,
)

# Standalone aliases retained for importers / QC text.
PANEL_HEIGHT_MM = STANDALONE_LAYOUT.panel_height_mm
TABLE_TOP_MM = STANDALONE_LAYOUT.table_top_mm
ROW_HEIGHT_MM = STANDALONE_LAYOUT.row_height_mm
TABLE_BOTTOM_MM = STANDALONE_LAYOUT.table_bottom_mm
SEPARATOR_MM = STANDALONE_LAYOUT.separator_mm
HEADER_TOP_MM = STANDALONE_LAYOUT.header_top_mm
HEADER_HEIGHT_MM = STANDALONE_LAYOUT.header_height_mm

DIR_X_MM = 2.5
DIR_W_MM = 3.2
ID_X_MM = 6.2
ID_W_MM = 5.5
LABEL_X_MM = 12.2
LABEL_W_MM = 36.3
CELL_X_MM = 49.0
CELL_W_MM = 38.0
COHORT_CELL_W_MM = CELL_W_MM / 3


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def direction_short(value: object) -> str:
    text = str(value)
    if text.startswith("LA"):
        return "LA"
    if text.startswith("RA"):
        return "RA"
    raise AssertionError(f"Unexpected frozen direction: {value}")


def classify_state(mean_af: float, sr_reference: float) -> str:
    if mean_af <= -0.02:
        return "Reversed"
    ratio = mean_af / sr_reference
    if ratio < 2 / 3:
        return "Attenuated"
    if ratio <= 4 / 3:
        return "Preserved"
    return "Amplified"


def validate_classifier_boundaries() -> None:
    reference = 1.0
    require(classify_state(-0.02, reference) == "Reversed", "Reversal boundary changed")
    require(
        classify_state(np.nextafter(-0.02, np.inf), reference) == "Attenuated",
        "Reversal classification order changed",
    )
    require(classify_state(2 / 3, reference) == "Preserved", "Lower preserved boundary changed")
    require(
        classify_state(np.nextafter(2 / 3, -np.inf), reference) == "Attenuated",
        "Attenuated boundary changed",
    )
    require(classify_state(4 / 3, reference) == "Preserved", "Upper preserved boundary changed")
    require(
        classify_state(np.nextafter(4 / 3, np.inf), reference) == "Amplified",
        "Amplified boundary changed",
    )


def derive_states(
    donor_source: pd.DataFrame,
    families: pd.DataFrame,
    paths: dict[str, Path],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float | int]]:
    validate_classifier_boundaries()
    require(donor_source.shape == (816, 10), "Corrected donor-family source is not 16 × 51")
    require(
        not any("fig2c" in str(path).lower() for path in paths.values()),
        "A Figure 2C state output was introduced as an input",
    )

    count_table = (
        donor_source[["dataset", "donor_id", "rhythm"]]
        .drop_duplicates()
        .groupby(["dataset", "rhythm"])["donor_id"]
        .nunique()
    )
    for cohort, counts in EXPECTED_COUNTS.items():
        for rhythm, expected in counts.items():
            require(
                int(count_table.loc[(cohort, rhythm)]) == expected,
                f"{cohort} {rhythm} donor count changed",
            )

    means = (
        donor_source.groupby(["family_id", "dataset", "rhythm"], as_index=False)
        .agg(
            mean_aCAIrel=("family_aCAIrel", "mean"),
            n_donors=("donor_id", "nunique"),
        )
        .pivot(index=["family_id", "dataset"], columns="rhythm")
    )
    means.columns = [f"{metric}_{rhythm}" for metric, rhythm in means.columns]
    means = means.reset_index().rename(columns={"dataset": "cohort"})

    metadata = families[
        ["family_id", "short_label", "canonical_direction", "canonical_family_display_order"]
    ].rename(
        columns={
            "short_label": "family_label",
            "canonical_family_display_order": "canonical_row_order",
        }
    )
    states = means.merge(metadata, on="family_id", how="left", validate="many_to_one")
    states["family_direction"] = states["canonical_direction"].map(direction_short)
    states["n_SR"] = states["n_donors_SR"].astype(int)
    states["n_AF"] = states["n_donors_AF"].astype(int)
    states["mean_SR_aCAIrel"] = states["mean_aCAIrel_SR"].astype(float)
    states["mean_AF_aCAIrel"] = states["mean_aCAIrel_AF"].astype(float)
    states["sr_reference"] = states["mean_SR_aCAIrel"].clip(lower=0.02)
    states["floor_applied"] = states["mean_SR_aCAIrel"].lt(0.02)
    states["af_sr_ratio"] = states["mean_AF_aCAIrel"] / states["sr_reference"]
    states["state"] = [
        classify_state(float(mean_af), float(reference))
        for mean_af, reference in zip(states["mean_AF_aCAIrel"], states["sr_reference"])
    ]
    state_rank = {state: index + 1 for index, state in enumerate(STATE_ORDER)}
    states["state_order"] = states["state"].map(state_rank).astype(int)
    cohort_rank = {cohort: index for index, cohort in enumerate(COHORT_ORDER)}
    states["_cohort_order"] = states["cohort"].map(cohort_rank)
    states = states.sort_values(["canonical_row_order", "_cohort_order"]).reset_index(drop=True)
    states = states[SOURCE_COLUMNS]

    require(states.shape == (48, len(SOURCE_COLUMNS)), "Expected exactly 48 family-cohort calls")
    require(states["family_id"].nunique() == 16, "Expected exactly 16 families")
    require(states["cohort"].nunique() == 3, "Expected exactly 3 cohorts")
    require(states["state"].notna().all(), "A state call is missing")
    require(
        states.groupby("family_id")["cohort"].nunique().eq(3).all(),
        "A family is missing a cohort call",
    )
    require(
        states.drop_duplicates("canonical_row_order")["family_id"].tolist()
        == EXPECTED_FAMILY_ORDER,
        "Family rows differ from the locked order",
    )
    for cohort, counts in EXPECTED_COUNTS.items():
        cohort_rows = states.loc[states["cohort"].eq(cohort)]
        require(cohort_rows["n_SR"].eq(counts["SR"]).all(), f"{cohort} SR count changed")
        require(cohort_rows["n_AF"].eq(counts["AF"]).all(), f"{cohort} AF count changed")

    authoritative = pd.read_csv(paths["family_acairel"], sep="\t")
    direct = (
        authoritative.loc[authoritative["dataset"].isin(COHORT_ORDER)]
        .groupby(["family_id", "dataset", "rhythm"])["family_aCAIrel"]
        .mean()
        .unstack("rhythm")
        .reset_index()
        .rename(columns={"dataset": "cohort", "SR": "direct_SR", "AF": "direct_AF"})
    )
    checked = states.merge(direct, on=["family_id", "cohort"], how="outer", validate="one_to_one")
    mean_error = float(
        np.max(
            np.abs(
                np.concatenate(
                    [
                        checked["mean_SR_aCAIrel"].to_numpy() - checked["direct_SR"].to_numpy(),
                        checked["mean_AF_aCAIrel"].to_numpy() - checked["direct_AF"].to_numpy(),
                    ]
                )
            )
        )
    )
    require(mean_error <= 1e-12, f"Cohort means do not reproduce corrected donor values: {mean_error}")

    runtime_counts = states["state"].value_counts()
    summary = pd.DataFrame(
        {
            "state": STATE_ORDER,
            "state_abbreviation": [STATE_ABBREVIATIONS[state] for state in STATE_ORDER],
            "state_order": np.arange(1, len(STATE_ORDER) + 1),
            "n_calls": [int(runtime_counts.get(state, 0)) for state in STATE_ORDER],
            "n_floor_applied": [
                int(states.loc[states["state"].eq(state), "floor_applied"].sum())
                for state in STATE_ORDER
            ],
        }
    )
    require(int(summary["n_calls"].sum()) == 48, "Runtime state counts do not sum to 48")
    require(
        int(summary["n_floor_applied"].sum()) == int(states["floor_applied"].sum()),
        "Runtime floor counts disagree",
    )
    stats: dict[str, float | int] = {
        "n_families": int(states["family_id"].nunique()),
        "n_cohorts": int(states["cohort"].nunique()),
        "n_calls": int(len(states)),
        "n_floor_applied": int(states["floor_applied"].sum()),
        "max_mean_reproduction_error": mean_error,
    }
    return states, summary, stats


def build_geometry(
    families: pd.DataFrame,
    donor_source: pd.DataFrame,
    layout: PanelLayout = STANDALONE_LAYOUT,
) -> pd.DataFrame:
    geometry = (
        families[
            ["family_id", "canonical_direction", "canonical_family_display_order"]
        ]
        .rename(
            columns={
                "canonical_family_display_order": "canonical_order",
            }
        )
        .sort_values("canonical_order")
        .reset_index(drop=True)
    )
    geometry["direction_block"] = geometry["canonical_direction"].map(direction_short)
    geometry["row_center"] = (
        layout.table_top_mm + (geometry["canonical_order"] - 0.5) * layout.row_height_mm
    )
    geometry["row_height"] = layout.row_height_mm
    geometry["LA_RA_separator_y"] = layout.separator_mm
    geometry["matrix_body_top"] = layout.table_top_mm
    geometry["matrix_body_bottom"] = layout.table_bottom_mm
    geometry = geometry[GEOMETRY_COLUMNS]

    require(len(geometry) == 16, "Shared row geometry is not 16 rows")
    require(geometry["family_id"].tolist() == EXPECTED_FAMILY_ORDER, "Geometry order changed")
    require(geometry["canonical_order"].tolist() == list(range(1, 17)), "Row ranks changed")
    require(geometry["direction_block"].tolist() == ["LA"] * 11 + ["RA"] * 5, "Direction blocks changed")
    require(geometry["row_height"].eq(layout.row_height_mm).all(), "Row heights changed")
    require(geometry["LA_RA_separator_y"].eq(layout.separator_mm).all(), "LA/RA separator changed")
    require(geometry["matrix_body_top"].eq(layout.table_top_mm).all(), "Matrix top changed")
    require(geometry["matrix_body_bottom"].eq(layout.table_bottom_mm).all(), "Matrix bottom changed")
    require(FIG2B_SOURCE.is_file(), "Corrected Figure 2B source is required for row-order validation")
    fig2b_rows = (
        pd.read_csv(FIG2B_SOURCE, sep="\t")
        [["family_id", "canonical_row_order"]]
        .drop_duplicates()
        .sort_values("canonical_row_order")
    )
    require(
        fig2b_rows["family_id"].astype(int).tolist() == geometry["family_id"].tolist(),
        "Shared geometry does not match Figure 2B row order",
    )
    return geometry


def load_row_geometry(path: Path = OUT_GEOMETRY) -> pd.DataFrame:
    """Load the mandatory physical row contract for Figure 2C/2D renderers."""
    require(path.is_file(), f"Missing shared Figure 2 row geometry: {path}")
    geometry = pd.read_csv(path, sep="\t")
    require(geometry.columns.tolist() == GEOMETRY_COLUMNS, "Shared geometry schema changed")
    require(len(geometry) == 16, "Shared geometry row count changed")
    require(geometry["family_id"].astype(int).tolist() == EXPECTED_FAMILY_ORDER, "Shared geometry order changed")
    require(geometry["canonical_order"].astype(int).tolist() == list(range(1, 17)), "Shared canonical order changed")
    require(geometry["direction_block"].tolist() == ["LA"] * 11 + ["RA"] * 5, "Shared direction blocks changed")
    return geometry


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 6.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def add_axes_mm(
    fig: plt.Figure,
    left: float,
    top: float,
    width: float,
    height: float,
    *,
    panel_height_mm: float,
) -> plt.Axes:
    return fig.add_axes(
        [
            left / PANEL_WIDTH_MM,
            (panel_height_mm - top - height) / panel_height_mm,
            width / PANEL_WIDTH_MM,
            height / panel_height_mm,
        ]
    )


def render(
    states: pd.DataFrame,
    summary: pd.DataFrame,
    geometry: pd.DataFrame,
    layout: PanelLayout = STANDALONE_LAYOUT,
) -> plt.Figure:
    set_style()
    panel_h = layout.panel_height_mm
    table_top = layout.table_top_mm
    row_h = layout.row_height_mm
    fig = plt.figure(figsize=(PANEL_WIDTH_MM * MM, panel_h * MM))
    fig.text(2.5 / PANEL_WIDTH_MM, 1 - 2.5 / panel_h, "C", fontsize=11, fontweight="bold", ha="left", va="top", color=FG)
    fig.text(
        8.0 / PANEL_WIDTH_MM,
        1 - 2.5 / panel_h,
        "Cohort-specific family rewiring states",
        fontsize=9.0,
        fontweight="bold",
        ha="left",
        va="top",
        color=FG,
    )
    fig.text(
        8.0 / PANEL_WIDTH_MM,
        1 - 7.3 / panel_h,
        "16 frozen NCP families × 3 cohorts = 48 descriptive state calls",
        fontsize=6.5,
        ha="left",
        va="top",
        color=MUTED,
    )
    fig.text(
        2.5 / PANEL_WIDTH_MM,
        1 - 12.5 / panel_h,
        "Descriptive AF-relative-to-SR states; not significance tests.",
        fontsize=6.5,
        ha="left",
        va="top",
        color=MUTED,
    )

    fig.text(
        ID_X_MM / PANEL_WIDTH_MM,
        1 - layout.label_header_y_mm / panel_h,
        "ID",
        fontsize=6.5,
        ha="left",
        va="center",
        color=MUTED,
    )
    fig.text(
        LABEL_X_MM / PANEL_WIDTH_MM,
        1 - layout.label_header_y_mm / panel_h,
        "Frozen NCP family",
        fontsize=6.5,
        ha="left",
        va="center",
        color=MUTED,
    )

    ax_header = add_axes_mm(
        fig,
        CELL_X_MM,
        layout.header_top_mm,
        CELL_W_MM,
        layout.header_height_mm,
        panel_height_mm=panel_h,
    )
    ax_header.set_xlim(0, 3)
    ax_header.set_ylim(0, 1)
    ax_header.axis("off")
    for index, cohort in enumerate(COHORT_ORDER):
        ax_header.add_patch(
            mpatches.Rectangle(
                (index, 0),
                1,
                1,
                facecolor="#C8C8C8",
                edgecolor="white",
                linewidth=0.6,
            )
        )
        counts = EXPECTED_COUNTS[cohort]
        ax_header.text(
            index + 0.5,
            0.65,
            cohort,
            ha="center",
            va="center",
            fontsize=6.5,
            fontfamily="Roboto Condensed",
            color=FG,
        )
        ax_header.text(
            index + 0.5,
            0.27,
            f'{counts["SR"]} SR / {counts["AF"]} AF',
            ha="center",
            va="center",
            fontsize=6.5,
            fontfamily="Roboto Condensed",
            color=FG,
        )

    axes = {
        "direction": add_axes_mm(fig, DIR_X_MM, table_top, DIR_W_MM, 16 * row_h, panel_height_mm=panel_h),
        "ids": add_axes_mm(fig, ID_X_MM, table_top, ID_W_MM, 16 * row_h, panel_height_mm=panel_h),
        "labels": add_axes_mm(fig, LABEL_X_MM, table_top, LABEL_W_MM, 16 * row_h, panel_height_mm=panel_h),
        "cells": add_axes_mm(fig, CELL_X_MM, table_top, CELL_W_MM, 16 * row_h, panel_height_mm=panel_h),
    }
    matrix_top = float(geometry["matrix_body_top"].iloc[0])
    matrix_bottom = float(geometry["matrix_body_bottom"].iloc[0])
    separator_y = float(geometry["LA_RA_separator_y"].iloc[0])
    for ax in axes.values():
        ax.set_ylim(matrix_bottom, matrix_top)
        ax.axis("off")

    ax_dir = axes["direction"]
    ax_dir.set_xlim(0, 1)
    for label, color in [("LA", LA_COLOR), ("RA", RA_COLOR)]:
        block = geometry.loc[geometry["direction_block"].eq(label)]
        start = float((block["row_center"] - block["row_height"] / 2).min())
        end = float((block["row_center"] + block["row_height"] / 2).max())
        ax_dir.add_patch(
            mpatches.Rectangle(
                (0, start),
                1,
                end - start,
                facecolor=mcolors.to_rgba(color, 0.22),
                edgecolor="none",
            )
        )
        ax_dir.text(
            0.5,
            (start + end) / 2,
            label,
            ha="center",
            va="center",
            rotation=90,
            fontsize=6.5,
            fontweight="bold",
            color=color,
        )

    ax_ids = axes["ids"]
    ax_ids.set_xlim(0, 1)
    ax_labels = axes["labels"]
    ax_labels.set_xlim(0, 1)
    row_meta = (
        states[
            [
                "family_id",
                "family_label",
                "family_direction",
                "canonical_row_order",
            ]
        ]
        .drop_duplicates()
        .merge(
            geometry,
            left_on=["family_id", "canonical_row_order", "family_direction"],
            right_on=["family_id", "canonical_order", "direction_block"],
            how="left",
            validate="one_to_one",
        )
        .sort_values("canonical_order")
    )
    require(row_meta["row_center"].notna().all(), "State rows do not resolve through shared geometry")
    for row in row_meta.itertuples(index=False):
        y = float(row.row_center)
        ax_ids.text(0.96, y, f"F{int(row.family_id)}", ha="right", va="center", fontsize=6.5, color=FG)
        ax_labels.text(0.01, y, str(row.family_label), ha="left", va="center", fontsize=6.5, color=FG)

    state_lookup = states.set_index(["family_id", "cohort"])
    ax_cells = axes["cells"]
    ax_cells.set_xlim(0, 3)
    for row in row_meta.itertuples(index=False):
        row_top = float(row.row_center - row.row_height / 2)
        row_height = float(row.row_height)
        for cohort_index, cohort in enumerate(COHORT_ORDER):
            call = state_lookup.loc[(int(row.family_id), cohort)]
            state = str(call["state"])
            ax_cells.add_patch(
                mpatches.Rectangle(
                    (cohort_index, row_top),
                    1,
                    row_height,
                    facecolor=STATE_COLORS[state],
                    edgecolor="white",
                    linewidth=0.6,
                )
            )
            ax_cells.text(
                cohort_index + 0.5,
                float(row.row_center),
                STATE_ABBREVIATIONS[state],
                ha="center",
                va="center",
                fontsize=6.5,
                fontweight="semibold",
                color=FG,
            )
            if bool(call["floor_applied"]):
                ax_cells.add_patch(
                    mpatches.Polygon(
                        [
                            (cohort_index + 1, row_top),
                            (cohort_index + 1 - WEDGE_FRACTION, row_top),
                            (cohort_index + 1, row_top + WEDGE_FRACTION * row_height),
                        ],
                        closed=True,
                        facecolor=MARKER_COLOR,
                        edgecolor="none",
                    )
                )

    ordinary_boundaries = (
        geometry["row_center"] + geometry["row_height"] / 2
    ).iloc[:-1]
    for ax in axes.values():
        for boundary in ordinary_boundaries:
            if not np.isclose(float(boundary), separator_y):
                ax.axhline(boundary, color=LIGHT_LINE, lw=0.25, zorder=5)
        ax.axhline(separator_y, color="#454545", lw=1.0, zorder=6)
    ax_cells.add_patch(
        mpatches.Rectangle(
            (0, matrix_top),
            3,
            matrix_bottom - matrix_top,
            fill=False,
            edgecolor="#A8A8A8",
            linewidth=0.45,
        )
    )

    ax_legend = add_axes_mm(
        fig, 2.5, layout.legend_top_mm, 83.0, 5.5, panel_height_mm=panel_h
    )
    ax_legend.set_xlim(0, 83)
    ax_legend.set_ylim(0, 1)
    ax_legend.axis("off")
    summary_lookup = summary.set_index("state")
    for index, state in enumerate(STATE_ORDER):
        x_start = index * 20.75
        ax_legend.add_patch(
            mpatches.Rectangle(
                (x_start + 0.1, 0.27),
                2.2,
                0.46,
                facecolor=STATE_COLORS[state],
                edgecolor="#777777",
                linewidth=0.25,
            )
        )
        ax_legend.text(
            x_start + 3.0,
            0.5,
            f"{state} (n={int(summary_lookup.loc[state, 'n_calls'])})",
            ha="left",
            va="center",
            fontsize=6.5,
            fontfamily="Roboto Condensed",
            color=FG,
        )

    floor_count = int(states["floor_applied"].sum())
    ax_floor = add_axes_mm(
        fig, 2.5, layout.floor_top_mm, 83.0, 4.0, panel_height_mm=panel_h
    )
    ax_floor.set_xlim(0, 1)
    ax_floor.set_ylim(0, 1)
    ax_floor.axis("off")
    sample_x, sample_y = 0.0, 0.21
    sample_width, sample_height = 0.10, 0.58
    ax_floor.add_patch(
        mpatches.Rectangle(
            (sample_x, sample_y),
            sample_width,
            sample_height,
            facecolor="white",
            edgecolor="#A8A8A8",
            linewidth=0.35,
        )
    )
    ax_floor.add_patch(
        mpatches.Polygon(
            [
                (sample_x + sample_width, sample_y + sample_height),
                (
                    sample_x + sample_width * (1 - WEDGE_FRACTION),
                    sample_y + sample_height,
                ),
                (
                    sample_x + sample_width,
                    sample_y + sample_height * (1 - WEDGE_FRACTION),
                ),
            ],
            closed=True,
            facecolor=MARKER_COLOR,
            edgecolor="none",
        )
    )
    ax_floor.text(
        0.115,
        0.53,
        f"corner wedge: SR-reference floor applied ({floor_count}/{len(states)} calls)",
        ha="left",
        va="center",
        fontsize=6.5,
        color=MUTED,
    )
    return fig


def save_outputs(
    fig: plt.Figure,
    *,
    out_base: Path = OUT_BASE,
    out_preview: Path = OUT_PREVIEW,
) -> dict[str, Path]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    outputs = {
        "pdf": out_base.with_suffix(".pdf"),
        "svg": out_base.with_suffix(".svg"),
        "png_600dpi": out_base.with_suffix(".png"),
        "preview": out_preview,
    }
    fig.savefig(outputs["pdf"])
    fig.savefig(outputs["svg"])
    fig.savefig(outputs["png_600dpi"], dpi=600)
    fig.savefig(outputs["preview"], dpi=200)
    plt.close(fig)
    return outputs


def validate_palette() -> dict[str, float]:
    direction_rgb = [np.asarray(mcolors.to_rgb(color)) for color in (LA_COLOR, RA_COLOR)]
    state_rgb = [np.asarray(mcolors.to_rgb(color)) for color in STATE_COLORS.values()]
    min_direction_distance = min(
        float(np.linalg.norm(state_color - direction_color))
        for state_color in state_rgb
        for direction_color in direction_rgb
    )
    luminance = [
        float(np.dot(color, np.asarray([0.2126, 0.7152, 0.0722])))
        for color in state_rgb
    ]
    min_grayscale_distance = min(
        abs(luminance[left] - luminance[right])
        for left in range(len(luminance))
        for right in range(left)
    )
    simulations = {
        "protanopia": np.asarray(
            [[0.56667, 0.43333, 0.0], [0.55833, 0.44167, 0.0], [0.0, 0.24167, 0.75833]]
        ),
        "deuteranopia": np.asarray(
            [[0.625, 0.375, 0.0], [0.700, 0.300, 0.0], [0.0, 0.300, 0.700]]
        ),
        "tritanopia": np.asarray(
            [[0.950, 0.050, 0.0], [0.0, 0.43333, 0.56667], [0.0, 0.475, 0.525]]
        ),
    }
    min_cvd_state_distance = np.inf
    min_cvd_direction_distance = np.inf
    for matrix in simulations.values():
        simulated_states = [matrix @ color for color in state_rgb]
        simulated_directions = [matrix @ color for color in direction_rgb]
        min_cvd_state_distance = min(
            min_cvd_state_distance,
            min(
                float(np.linalg.norm(simulated_states[left] - simulated_states[right]))
                for left in range(len(simulated_states))
                for right in range(left)
            ),
        )
        min_cvd_direction_distance = min(
            min_cvd_direction_distance,
            min(
                float(np.linalg.norm(state_color - direction_color))
                for state_color in simulated_states
                for direction_color in simulated_directions
            ),
        )
    require(min_direction_distance >= 0.15, "A state fill is too similar to a frozen-direction color")
    require(min_grayscale_distance >= 0.045, "State fills are insufficiently separated in grayscale")
    require(min_cvd_state_distance >= 0.18, "State fills collapse under a common color-vision simulation")
    require(
        min_cvd_direction_distance >= 0.10,
        "State and frozen-direction colors collapse under color-vision simulation",
    )
    require(set(STATE_COLORS.values()).isdisjoint({LA_COLOR, RA_COLOR}), "Direction color reused as state fill")
    return {
        "normal_state_direction": min_direction_distance,
        "grayscale_state": min_grayscale_distance,
        "cvd_state": float(min_cvd_state_distance),
        "cvd_state_direction": float(min_cvd_direction_distance),
    }


def validate_outputs(
    outputs: dict[str, Path],
    layout: PanelLayout = STANDALONE_LAYOUT,
) -> None:
    pdfinfo = shutil.which("pdfinfo")
    pdftotext = shutil.which("pdftotext")
    require(pdfinfo is not None and pdftotext is not None, "PDF QC tools are unavailable")
    info = subprocess.run(
        [pdfinfo, str(outputs["pdf"])], check=True, capture_output=True, text=True
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read Figure 2C PDF dimensions")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    require(
        abs(width_mm - PANEL_WIDTH_MM) < 0.05
        and abs(height_mm - layout.panel_height_mm) < 0.05,
        f"Figure 2C PDF dimensions are {width_mm:.2f} × {height_mm:.2f} mm",
    )
    visible = subprocess.run(
        [pdftotext, "-layout", str(outputs["pdf"]), "-"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    normalized = " ".join(visible.split())
    required = [
        "Cohort-specific family rewiring states",
        "16 frozen NCP families × 3 cohorts = 48 descriptive state calls",
        "Descriptive AF-relative-to-SR states; not significance tests.",
        "GSE128188",
        "GSE79768",
        "GSE115574",
        "Amplified",
        "Preserved",
        "Attenuated",
        "Reversed",
        "corner wedge: SR-reference floor applied (5/48 calls)",
    ]
    require(all(phrase in normalized for phrase in required), "Figure 2C visible text is incomplete")
    forbidden = [
        "*",
        "p value",
        "p-value",
        "Hedges",
        "AF/SR",
        "AF mean",
        "SR mean",
        "significance star",
    ]
    require(not any(token in visible for token in forbidden), "Prohibited inferential or numeric cell text is visible")
    require(outputs["svg"].is_file() and outputs["svg"].stat().st_size > 0, "SVG export is missing")
    expected_width = round(PANEL_WIDTH_MM / 25.4 * 600)
    expected_height = round(layout.panel_height_mm / 25.4 * 600)
    image = plt.imread(outputs["png_600dpi"])
    require(
        abs(image.shape[1] - expected_width) <= 1
        and abs(image.shape[0] - expected_height) <= 1,
        f"600-dpi PNG dimensions changed: {image.shape[1]} × {image.shape[0]}",
    )


def write_qc(
    manifest: pd.DataFrame,
    paths: dict[str, Path],
    states: pd.DataFrame,
    summary: pd.DataFrame,
    stats: dict[str, float | int],
    outputs: dict[str, Path],
    palette_stats: dict[str, float],
    composite_outputs: dict[str, Path],
) -> None:
    input_lines = [
        f"- `{path.relative_to(ROOT)}`: `{fig2b.sha256(path)}`"
        for path in sorted(paths.values(), key=lambda item: str(item))
    ]
    output_paths = {
        "source": OUT_SOURCE,
        "summary": OUT_SUMMARY,
        "geometry": OUT_GEOMETRY,
        "geometry_composite": OUT_GEOMETRY_COMPOSITE,
        **outputs,
        **{f"composite_{k}": v for k, v in composite_outputs.items()},
    }
    output_lines = [
        f"- `{path.relative_to(ROOT)}`: `{fig2b.sha256(path)}`"
        for path in output_paths.values()
    ]
    state_lines = [
        f"- {row.state}: {int(row.n_calls)} calls; {int(row.n_floor_applied)} floor-applied"
        for row in summary.itertuples(index=False)
    ]
    cohort_lines = [
        f"- {cohort}: {counts['SR']} SR / {counts['AF']} AF"
        for cohort, counts in EXPECTED_COUNTS.items()
    ]
    text = "\n".join(
        [
            "# Figure 2C corrected family-cohort state QC",
            "",
            "Status: **PASS**",
            "",
            "## Matrix and derivation",
            f"- Families: {stats['n_families']}",
            f"- Cohorts: {stats['n_cohorts']}",
            f"- Descriptive calls: {stats['n_calls']}",
            f"- SR-reference floor applied: {stats['n_floor_applied']} calls",
            f"- Maximum cohort-mean reproduction error: {float(stats['max_mean_reproduction_error']):.3e}",
            "- Donor-level family aCAIrel independently verified against equal-weight corrected constituent-program aCAIrel.",
            "- State rule order: AF mean <= -0.02 first; otherwise cohort SR reference = max(SR mean, 0.02), then exact 2/3 and 4/3 boundaries.",
            "- No pooled SR reference, sample-level scoring, healthy-reference reconstruction, family reordering, or significance testing was used.",
            "",
            "## Donor counts",
            *cohort_lines,
            "",
            "## Runtime state counts",
            *state_lines,
            "",
            "## Geometry and rendering",
            "- Shared panel width: 88 mm",
            f"- Standalone row height: {STANDALONE_LAYOUT.row_height_mm:.2f} mm; panel {STANDALONE_LAYOUT.panel_height_mm:.0f} mm.",
            f"- Composite row height: {COMPOSITE_LAYOUT.row_height_mm:.2f} mm; panel {COMPOSITE_LAYOUT.panel_height_mm:.0f} mm (shared with Figure 2D collage).",
            f"- Standalone table top/bottom: {STANDALONE_LAYOUT.table_top_mm:.2f} / {STANDALONE_LAYOUT.table_bottom_mm:.2f} mm",
            f"- Composite table top/bottom: {COMPOSITE_LAYOUT.table_top_mm:.2f} / {COMPOSITE_LAYOUT.table_bottom_mm:.2f} mm",
            "- Locked rows: 11 LA followed by 5 RA; identical to corrected Figure 2B.",
            "- Standalone geometry units: millimetres from the 102 mm panel top.",
            "- Composite geometry written to `figure2_family_row_geometry_composite.tsv` for side-by-side alignment with Figure 2D.",
            f"- Minimum state-to-direction RGB distance: {palette_stats['normal_state_direction']:.3f}",
            f"- Minimum grayscale state luminance distance: {palette_stats['grayscale_state']:.3f}",
            f"- Minimum simulated-CVD state distance: {palette_stats['cvd_state']:.3f}",
            f"- Minimum simulated-CVD state-to-direction distance: {palette_stats['cvd_state_direction']:.3f}",
            "- All cells retain state abbreviations for grayscale interpretation.",
            "- No significance star, p value, Hedges' g, AF/SR ratio, AF mean, or SR mean is rendered.",
            "",
            "## Dependency audit",
            f"- Corrected manifest rows used: {len(manifest)}",
            "- No Figure 2C state table appears in the runtime input set.",
            "- Legacy Figure 2C renderers and state tables were not read.",
            "",
            "## Input SHA256",
            *input_lines,
            "",
            "## Output SHA256",
            *output_lines,
            "",
        ]
    )
    OUT_QC.write_text(text, encoding="utf-8")


def _assert_render_text(fig: plt.Figure, states: pd.DataFrame, summary: pd.DataFrame) -> None:
    fig.canvas.draw()
    visible_text = [
        text
        for text in fig.texts
        + [item for ax in fig.axes for item in ax.texts]
        if text.get_text().strip()
    ]
    require(min(text.get_fontsize() for text in visible_text) >= 6.5, "Visible text below 6.5 pt")
    require(not any("*" in text.get_text() for text in visible_text), "A significance star is rendered")
    abbreviation_counts = {
        abbreviation: sum(text.get_text() == abbreviation for text in visible_text)
        for abbreviation in STATE_ABBREVIATIONS.values()
    }
    expected_abbreviation_counts = {
        STATE_ABBREVIATIONS[row.state]: int(row.n_calls)
        for row in summary.itertuples(index=False)
    }
    require(
        abbreviation_counts == expected_abbreviation_counts,
        "Rendered cell abbreviations do not match runtime state counts",
    )
    require(not any(text.get_text() == "Pre" for text in visible_text), "Legacy Pre abbreviation remains")
    preserved_labels = [text for text in visible_text if text.get_text() == "Pres"]
    require(
        all(
            text.get_window_extent().width
            < abs(
                text.axes.transData.transform((1, 0))[0]
                - text.axes.transData.transform((0, 0))[0]
            )
            * 0.8
            for text in preserved_labels
        ),
        "Pres does not fit comfortably within a state cell",
    )
    wedge_count = sum(
        isinstance(patch, mpatches.Polygon)
        for ax in fig.axes
        for patch in ax.patches
    )
    require(
        wedge_count == int(states["floor_applied"].sum()) + 1,
        "Rendered floor wedges do not match floor flags plus legend sample",
    )


def main() -> None:
    manifest, paths = fig2b.load_manifest()
    donor_source, families, donor_stats = fig2b.build_source(paths)
    states, summary, stats = derive_states(donor_source, families, paths)
    require(
        float(donor_stats["max_verification_error"]) <= 1e-12,
        "Donor-family validation failed upstream",
    )
    geometry = build_geometry(families, donor_source, STANDALONE_LAYOUT)
    geometry_composite = build_geometry(families, donor_source, COMPOSITE_LAYOUT)

    RESULTS.mkdir(parents=True, exist_ok=True)
    states.to_csv(OUT_SOURCE, sep="\t", index=False, float_format="%.15g")
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    geometry.to_csv(OUT_GEOMETRY, sep="\t", index=False, float_format="%.6f")
    geometry_composite.to_csv(OUT_GEOMETRY_COMPOSITE, sep="\t", index=False, float_format="%.6f")
    geometry = load_row_geometry()
    geometry_composite = load_row_geometry(OUT_GEOMETRY_COMPOSITE)

    fig = render(states, summary, geometry, STANDALONE_LAYOUT)
    _assert_render_text(fig, states, summary)
    outputs = save_outputs(fig)
    validate_outputs(outputs, STANDALONE_LAYOUT)

    fig_composite = render(states, summary, geometry_composite, COMPOSITE_LAYOUT)
    _assert_render_text(fig_composite, states, summary)
    composite_outputs = save_outputs(
        fig_composite,
        out_base=OUT_BASE_COMPOSITE,
        out_preview=OUT_PREVIEW_COMPOSITE,
    )
    validate_outputs(composite_outputs, COMPOSITE_LAYOUT)

    palette_stats = validate_palette()
    write_qc(
        manifest,
        paths,
        states,
        summary,
        stats,
        outputs,
        palette_stats,
        composite_outputs,
    )

    counts = ", ".join(
        f"{row.state}={int(row.n_calls)}" for row in summary.itertuples(index=False)
    )
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)} ({len(states)} calls)")
    print(f"Wrote {OUT_GEOMETRY.relative_to(ROOT)} ({len(geometry)} locked rows)")
    print(f"Wrote {OUT_GEOMETRY_COMPOSITE.relative_to(ROOT)} ({len(geometry_composite)} composite rows)")
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_BASE_COMPOSITE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW_COMPOSITE.relative_to(ROOT)}")
    print(f"Runtime counts: {counts}; floor applied={int(stats['n_floor_applied'])}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
