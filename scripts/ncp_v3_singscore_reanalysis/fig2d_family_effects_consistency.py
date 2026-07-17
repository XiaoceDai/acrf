#!/usr/bin/env python3
"""Final Figure 2D: cross-cohort AF-SR effects and consistency.

This panel replaces the legacy main-text chamber-decomposition Figure 2D and
the legacy separate Figure 2E consistency panel. It shows, for each of the 16
frozen NCP families, cohort-specific and pooled Hedges' g (AF - SR) together
with the corrected High/Moderate/Low consistency class.

Inputs are read through the dedicated corrected manifest used by Figures 2A-C.
No legacy effect or consistency table is read as a plotting input. The family-
level pooling model is audited at runtime: only fixed-effect inverse-variance
weighting (no tau2) is accepted for the pooled estimates.
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fig2b_family_donor_heatmap as fig2b
import fig2c_family_cohort_state_map as fig2c


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"

OUT_SOURCE = RESULTS / "fig2d_family_effects_consistency_source.tsv"
OUT_SUMMARY = RESULTS / "fig2d_family_effects_consistency_summary.tsv"
OUT_QC = RESULTS / "fig2d_family_effects_consistency_qc.md"
OUT_AUDIT = RESULTS / "legacy_vs_corrected_consistency_changes.tsv"
OUT_BASE = FIGURES / "fig2d_family_effects_consistency"
OUT_PREVIEW = FIGURES / "fig2d_family_effects_consistency_preview.png"

COHORT_ORDER = ["GSE128188", "GSE79768", "GSE115574"]
EXPECTED_COUNTS = {
    "GSE128188": {"SR": 5, "AF": 5},
    "GSE79768": {"SR": 6, "AF": 7},
    "GSE115574": {"SR": 15, "AF": 13},
}
EXPECTED_FAMILY_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]

MM = 1 / 25.4
PANEL_WIDTH_MM = 88.0
# Taller panel: row height unlocked from standalone Fig 2C so cohort CIs can be shown.
# Matches fig2c COMPOSITE_LAYOUT (shared C|D collage contract).
PANEL_HEIGHT_MM = fig2c.COMPOSITE_PANEL_HEIGHT_MM
TABLE_TOP_MM = fig2c.COMPOSITE_TABLE_TOP_MM
ROW_HEIGHT_MM = fig2c.COMPOSITE_ROW_HEIGHT_MM
TABLE_BOTTOM_MM = TABLE_TOP_MM + 16 * ROW_HEIGHT_MM
SEPARATOR_MM = TABLE_TOP_MM + 11 * ROW_HEIGHT_MM
# Shared Fig 2C standalone geometry used only for family order / LA-RA block membership.
FIG2C_ROW_HEIGHT_MM = fig2c.STANDALONE_LAYOUT.row_height_mm
OUT_GEOMETRY_COMPOSITE = fig2c.OUT_GEOMETRY_COMPOSITE

DIR_X_MM = 2.5
DIR_W_MM = 3.2
ID_X_MM = 6.2
ID_W_MM = 5.5
FOREST_X_MM = 14.0
FOREST_W_MM = 52.0
COL_X_MM = 66.0
COL_W_MM = 22.0
HEADER_TOP_MM = 19.0
HEADER_HEIGHT_MM = TABLE_TOP_MM - HEADER_TOP_MM

FG = "#202020"
MUTED = "#606060"
LIGHT_LINE = "#EEEEEE"
LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
DIRECTION_ALPHA = 0.22

MARKER_COLOR = "#333333"
COHORT_MARKERS = {
    "GSE128188": "o",
    "GSE79768": "^",
    "GSE115574": "s",
}
COHORT_COLORS = {
    "GSE128188": "#A0A0A0",
    "GSE79768": "#808080",
    "GSE115574": "#606060",
}
CONSISTENCY_COLORS = {
    "High": "#3C6E47",
    "Moderate": "#C69C3E",
    "Low": "#BFBFBF",
}
CONSISTENCY_SYMBOLS = {
    "High": "●",
    "Moderate": "▲",
    "Low": "■",
}

SOURCE_COLUMNS = [
    "family_id",
    "family_label",
    "family_direction",
    "canonical_order",
    "estimate_type",
    "cohort",
    "n_AF",
    "n_SR",
    "hedges_g",
    "sampling_variance",
    "standard_error",
    "ci_lower",
    "ci_upper",
    "pooling_model",
    "tau2",
    "I2",
    "pooled_supported",
    "pooled_sign",
    "cohort_sign",
    "sign_matches_pooled",
    "n_sign_match",
    "consistency_class",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def direction_short(value: object) -> str:
    text = str(value)
    if text.startswith("LA"):
        return "LA"
    if text.startswith("RA"):
        return "RA"
    raise AssertionError(f"Unexpected frozen direction: {value}")


def hedges_g(af: np.ndarray, sr: np.ndarray) -> tuple[float, float]:
    """Small-sample-corrected Hedges' g and its sampling variance."""
    nx, ny = len(af), len(sr)
    sx, sy = np.std(af, ddof=1), np.std(sr, ddof=1)
    pooled_sd = np.sqrt(((nx - 1) * sx**2 + (ny - 1) * sy**2) / (nx + ny - 2))
    correction = 1 - 3 / (4 * (nx + ny - 2) - 1)
    g = correction * (np.mean(af) - np.mean(sr)) / pooled_sd
    variance = (nx + ny) / (nx * ny) + g**2 / (2 * (nx + ny))
    return float(g), float(variance)


def round_axis_bound(value: float) -> float:
    """Round a positive bound outward to a readable symmetric axis limit."""
    if value <= 0:
        return 1.0
    # Choose a human-friendly step based on magnitude.
    if value <= 0.5:
        step = 0.1
    elif value <= 1.0:
        step = 0.25
    elif value <= 2.0:
        step = 0.5
    elif value <= 3.0:
        step = 0.5
    elif value <= 4.0:
        step = 1.0
    elif value <= 8.0:
        step = 2.0
    else:
        step = 5.0
    return float(np.ceil(value / step) * step)


def plotted_axis_bound(effects: pd.DataFrame) -> float:
    """Symmetric limit from all drawn CI bounds (cohort + pooled; no clipping)."""
    max_abs = float(effects[["ci_lower", "ci_upper"]].abs().max().max())
    return round_axis_bound(max_abs)


def expand_row_geometry(geometry_2c: pd.DataFrame) -> pd.DataFrame:
    """Rebuild row centers with taller rows while preserving Fig 2C order/blocks.

    Family order and LA/RA membership stay locked to Figure 2C. Physical row
    height is increased so cohort CIs can be shown without stacking. The
    resulting centres must match `figure2_family_row_geometry_composite.tsv`.
    """
    geo = geometry_2c.copy().sort_values("canonical_order").reset_index(drop=True)
    require(
        geo["family_id"].astype(int).tolist() == EXPECTED_FAMILY_ORDER,
        "Shared geometry order differs from expected",
    )
    require(len(geo) == 16, "Expected 16 families in shared geometry")

    geo["row_height"] = ROW_HEIGHT_MM
    geo["row_center"] = TABLE_TOP_MM + (np.arange(16) + 0.5) * ROW_HEIGHT_MM
    geo["matrix_body_top"] = TABLE_TOP_MM
    geo["matrix_body_bottom"] = TABLE_BOTTOM_MM
    geo["LA_RA_separator_y"] = SEPARATOR_MM
    # Keep direction_block from Fig 2C; assert it still splits after family 11 in order.
    la_count = int((geo["direction_block"] == "LA").sum())
    require(la_count == 11, f"Expected 11 LA families, got {la_count}")
    require(
        float(geo.loc[geo["direction_block"].eq("LA"), "row_center"].max())
        < SEPARATOR_MM
        < float(geo.loc[geo["direction_block"].eq("RA"), "row_center"].min()),
        "LA/RA separator does not fall between direction blocks",
    )
    return geo


def assert_matches_composite_geometry(geometry: pd.DataFrame) -> None:
    """Require 2D row centres equal the shared composite C|D geometry TSV."""
    composite = fig2c.load_row_geometry(OUT_GEOMETRY_COMPOSITE)
    require(
        geometry["family_id"].astype(int).tolist()
        == composite["family_id"].astype(int).tolist(),
        "Figure 2D family order diverges from composite geometry",
    )
    require(
        np.allclose(
            geometry["row_center"].to_numpy(dtype=float),
            composite["row_center"].to_numpy(dtype=float),
            atol=1e-6,
        ),
        "Figure 2D row centers diverge from composite geometry contract",
    )
    require(
        np.allclose(
            geometry["row_height"].to_numpy(dtype=float),
            composite["row_height"].to_numpy(dtype=float),
            atol=1e-6,
        ),
        "Figure 2D row heights diverge from composite geometry contract",
    )
    require(
        abs(float(geometry["matrix_body_top"].iloc[0]) - float(composite["matrix_body_top"].iloc[0]))
        < 1e-6
        and abs(
            float(geometry["matrix_body_bottom"].iloc[0])
            - float(composite["matrix_body_bottom"].iloc[0])
        )
        < 1e-6,
        "Figure 2D matrix body bounds diverge from composite geometry",
    )


def compute_effects(
    donor_source: pd.DataFrame,
    families: pd.DataFrame,
) -> pd.DataFrame:
    """Compute cohort-specific and pooled Hedges' g per family."""
    rows: list[dict[str, Any]] = []
    for family_id in EXPECTED_FAMILY_ORDER:
        family_group = donor_source.loc[donor_source["family_id"] == family_id]
        family_meta = families.loc[families["family_id"] == family_id].iloc[0]
        family_label = str(family_meta["short_label"])
        family_direction = direction_short(family_meta["canonical_direction"])
        canonical_order = int(family_meta["canonical_family_display_order"])

        cohort_records: list[dict[str, Any]] = []
        for cohort in COHORT_ORDER:
            cohort_group = family_group.loc[family_group["dataset"] == cohort]
            af = cohort_group.loc[cohort_group["rhythm"] == "AF", "family_aCAIrel"].to_numpy(dtype=float)
            sr = cohort_group.loc[cohort_group["rhythm"] == "SR", "family_aCAIrel"].to_numpy(dtype=float)
            g, variance = hedges_g(af, sr)
            record = {
                "family_id": int(family_id),
                "family_label": family_label,
                "family_direction": family_direction,
                "canonical_order": canonical_order,
                "estimate_type": "cohort-specific",
                "cohort": cohort,
                "n_AF": len(af),
                "n_SR": len(sr),
                "hedges_g": g,
                "sampling_variance": variance,
                "standard_error": float(np.sqrt(variance)),
                "ci_lower": g - 1.96 * np.sqrt(variance),
                "ci_upper": g + 1.96 * np.sqrt(variance),
                "pooling_model": "NA",
                "tau2": "NA",
                "I2": "NA",
            }
            rows.append(record)
            cohort_records.append(record)

        valid = [r for r in cohort_records if np.isfinite(r["sampling_variance"]) and r["sampling_variance"] > 0]
        weights = np.array([1 / r["sampling_variance"] for r in valid])
        pooled_g = float(np.sum(weights * [r["hedges_g"] for r in valid]) / weights.sum())
        pooled_var = float(1 / weights.sum())
        pooled_se = float(np.sqrt(pooled_var))
        ci_lower = pooled_g - 1.96 * pooled_se
        ci_upper = pooled_g + 1.96 * pooled_se

        # I2 based on the fixed-effect Q statistic (Higgins-Thompson).
        subset = pd.DataFrame(valid)
        q = float(np.sum(weights * (subset["hedges_g"].to_numpy() - pooled_g) ** 2))
        df = len(valid) - 1
        i2 = float(max(0.0, (q - df) / q) if q > 0 else 0.0)

        pooled_record = {
            "family_id": int(family_id),
            "family_label": family_label,
            "family_direction": family_direction,
            "canonical_order": canonical_order,
            "estimate_type": "pooled",
            "cohort": "POOLED",
            "n_AF": sum(r["n_AF"] for r in valid),
            "n_SR": sum(r["n_SR"] for r in valid),
            "hedges_g": pooled_g,
            "sampling_variance": pooled_var,
            "standard_error": pooled_se,
            "ci_lower": ci_lower,
            "ci_upper": ci_upper,
            "pooling_model": "fixed-effect_inverse_variance",
            "tau2": "NA",
            "I2": i2,
        }
        rows.append(pooled_record)

    effects = pd.DataFrame(rows)
    require(effects.shape[0] == 64, f"Expected 64 effect rows, got {effects.shape[0]}")
    require(
        effects.groupby("family_id").size().eq(4).all(),
        "Each family must have 3 cohort + 1 pooled estimate",
    )
    require(
        effects.groupby("family_id")["cohort"].apply(lambda x: set(x) == set(COHORT_ORDER + ["POOLED"])).all(),
        "Cohort set per family is incomplete",
    )
    return effects


def compute_consistency(effects: pd.DataFrame) -> pd.DataFrame:
    """Assign corrected High/Moderate/Low consistency classes."""
    summary_rows: list[dict[str, Any]] = []
    for family_id in EXPECTED_FAMILY_ORDER:
        subset = effects.loc[
            (effects["family_id"] == family_id) & (effects["cohort"].isin(COHORT_ORDER))
        ].copy()
        pooled = effects.loc[
            (effects["family_id"] == family_id) & (effects["cohort"] == "POOLED")
        ].iloc[0]

        pooled_sign = int(np.sign(pooled["hedges_g"]))
        pooled_supported = bool((pooled["ci_lower"] > 0) or (pooled["ci_upper"] < 0))
        sign_matches = (np.sign(subset["hedges_g"]) == pooled_sign).astype(int)
        n_sign_match = int(sign_matches.sum())

        if pooled_supported and n_sign_match >= 2:
            consistency_class = "High"
        elif n_sign_match >= 2:
            consistency_class = "Moderate"
        else:
            consistency_class = "Low"

        summary_rows.append(
            {
                "family_id": int(family_id),
                "pooled_g": pooled["hedges_g"],
                "pooled_ci_lower": pooled["ci_lower"],
                "pooled_ci_upper": pooled["ci_upper"],
                "pooled_supported": pooled_supported,
                "pooled_sign": pooled_sign,
                "n_sign_match": n_sign_match,
                "consistency_class": consistency_class,
                "I2": pooled["I2"],
            }
        )

        # Write per-row classification metadata back into the effects frame.
        effects.loc[pooled.name, "pooled_supported"] = pooled_supported
        effects.loc[pooled.name, "pooled_sign"] = pooled_sign
        effects.loc[pooled.name, "n_sign_match"] = n_sign_match
        effects.loc[pooled.name, "consistency_class"] = consistency_class

        for _, row in subset.iterrows():
            effects.loc[row.name, "pooled_supported"] = pooled_supported
            effects.loc[row.name, "pooled_sign"] = pooled_sign
            effects.loc[row.name, "cohort_sign"] = int(np.sign(row["hedges_g"]))
            effects.loc[row.name, "sign_matches_pooled"] = bool(int(np.sign(row["hedges_g"])) == pooled_sign)
            effects.loc[row.name, "n_sign_match"] = n_sign_match
            effects.loc[row.name, "consistency_class"] = consistency_class

    summary = pd.DataFrame(summary_rows)
    return summary


def build_legacy_audit(
    effects: pd.DataFrame,
    families: pd.DataFrame,
) -> pd.DataFrame:
    """Build the 7-row legacy-to-corrected consistency migration audit."""
    legacy_path = ROOT / "results" / "fig2e_replication_lollipop_v3.tsv"
    require(legacy_path.is_file(), f"Legacy consistency table missing: {legacy_path}")
    legacy = pd.read_csv(legacy_path, sep="\t")
    required_legacy_cols = {
        "family_id",
        "replication_class",
        "same_sign_n",
        "pooled_sig",
        "pooled_g",
        "pooled_ci_lo",
        "pooled_ci_hi",
    }
    require(
        required_legacy_cols.issubset(legacy.columns),
        f"Legacy table missing columns: {required_legacy_cols - set(legacy.columns)}",
    )

    corrected = effects.loc[effects["cohort"] == "POOLED"].copy()
    corrected["consistency_class"] = corrected["consistency_class"].str.capitalize()
    legacy["replication_class"] = legacy["replication_class"].str.capitalize()

    # Rename legacy columns to avoid collisions and produce readable output.
    legacy_renamed = legacy[
        [
            "family_id",
            "replication_class",
            "same_sign_n",
            "pooled_sig",
            "pooled_g",
            "pooled_ci_lo",
            "pooled_ci_hi",
        ]
    ].rename(
        columns={
            "replication_class": "legacy_consistency_class",
            "same_sign_n": "legacy_n_sign_match",
            "pooled_sig": "legacy_pooled_supported",
            "pooled_g": "legacy_pooled_g",
            "pooled_ci_lo": "legacy_pooled_ci_lower",
            "pooled_ci_hi": "legacy_pooled_ci_upper",
        }
    )
    corrected_renamed = corrected[
        [
            "family_id",
            "consistency_class",
            "n_sign_match",
            "pooled_supported",
            "hedges_g",
            "ci_lower",
            "ci_upper",
        ]
    ].rename(
        columns={
            "consistency_class": "corrected_consistency_class",
            "n_sign_match": "corrected_n_sign_match",
            "pooled_supported": "corrected_pooled_supported",
            "hedges_g": "corrected_pooled_g",
            "ci_lower": "corrected_pooled_ci_lower",
            "ci_upper": "corrected_pooled_ci_upper",
        }
    )

    merged = corrected_renamed.merge(legacy_renamed, on="family_id", how="inner")
    changed = merged.loc[merged["corrected_consistency_class"].ne(merged["legacy_consistency_class"])].copy()

    label_map = families.set_index("family_id")["short_label"].to_dict()
    changed["family_label"] = changed["family_id"].map(label_map)
    audit = changed[
        [
            "family_id",
            "family_label",
            "legacy_consistency_class",
            "corrected_consistency_class",
            "legacy_n_sign_match",
            "corrected_n_sign_match",
            "legacy_pooled_supported",
            "corrected_pooled_supported",
            "legacy_pooled_g",
            "corrected_pooled_g",
            "legacy_pooled_ci_lower",
            "legacy_pooled_ci_upper",
            "corrected_pooled_ci_lower",
            "corrected_pooled_ci_upper",
        ]
    ]
    require(len(audit) == 7, f"Expected exactly 7 changed families, got {len(audit)}")
    return audit


def validate_model_lock(effects: pd.DataFrame) -> None:
    """Assert the pooled estimates come from the locked fixed-effect IVW model."""
    pooled = effects.loc[effects["cohort"] == "POOLED"]
    require(
        (pooled["pooling_model"] == "fixed-effect_inverse_variance").all(),
        "Pooled estimates do not use the locked fixed-effect IVW model",
    )
    require(
        (pooled["tau2"] == "NA").all(),
        "Family-level tau2 is not NA; a random-effects model may have been used",
    )
    require(
        pooled["hedges_g"].notna().all() and np.isfinite(pooled["hedges_g"]).all(),
        "Pooled Hedges' g contains non-finite values",
    )
    require(
        pooled["ci_lower"].notna().all() and pooled["ci_upper"].notna().all(),
        "Pooled CIs contain missing values",
    )
    require(
        pooled["I2"].notna().all() and (pooled["I2"] >= 0).all() and (pooled["I2"] <= 1).all(),
        "I2 is outside [0, 1]",
    )


def add_axes_mm(
    fig: plt.Figure,
    left: float,
    top: float,
    width: float,
    height: float,
) -> plt.Axes:
    return fig.add_axes(
        [
            left / PANEL_WIDTH_MM,
            (PANEL_HEIGHT_MM - top - height) / PANEL_HEIGHT_MM,
            width / PANEL_WIDTH_MM,
            height / PANEL_HEIGHT_MM,
        ]
    )


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


def render(
    effects: pd.DataFrame,
    summary: pd.DataFrame,
    geometry: pd.DataFrame,
    *,
    composite_mode: bool = False,
) -> plt.Figure:
    set_style()
    fig = plt.figure(figsize=(PANEL_WIDTH_MM * MM, PANEL_HEIGHT_MM * MM))

    # Panel label and title (compressed, no explanatory note in figure body)
    fig.text(
        2.5 / PANEL_WIDTH_MM,
        1 - 2.5 / PANEL_HEIGHT_MM,
        "D",
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="top",
        color=FG,
    )
    fig.text(
        8.0 / PANEL_WIDTH_MM,
        1 - 2.5 / PANEL_HEIGHT_MM,
        "Cross-cohort AF-SR effects and consistency",
        fontsize=9.0,
        fontweight="bold",
        ha="left",
        va="top",
        color=FG,
    )
    fig.text(
        8.0 / PANEL_WIDTH_MM,
        1 - 7.3 / PANEL_HEIGHT_MM,
        "Cohort-specific and pooled Hedges' g across 16 frozen NCP families",
        fontsize=6.5,
        ha="left",
        va="top",
        color=MUTED,
    )

    matrix_top = float(geometry["matrix_body_top"].iloc[0])
    matrix_bottom = float(geometry["matrix_body_bottom"].iloc[0])
    separator_y = float(geometry["LA_RA_separator_y"].iloc[0])
    row_lookup = geometry.set_index("family_id")
    summary_lookup = summary.set_index("family_id")

    # Determine symmetric x-axis from all CI bounds (no clipping).
    x_limit = plotted_axis_bound(effects)

    # Layout differs slightly between standalone and composite versions.
    if composite_mode:
        dir_x_mm = None
        dir_w_mm = None
        id_x_mm = 2.5
        id_w_mm = 5.5
        forest_x_mm = 9.0
        forest_w_mm = 57.0
    else:
        dir_x_mm = DIR_X_MM
        dir_w_mm = DIR_W_MM
        id_x_mm = ID_X_MM
        id_w_mm = ID_W_MM
        forest_x_mm = FOREST_X_MM
        forest_w_mm = FOREST_W_MM

    col_x_mm = COL_X_MM
    col_w_mm = COL_W_MM

    # Direction strip (standalone only)
    if not composite_mode and dir_x_mm is not None and dir_w_mm is not None:
        ax_dir = add_axes_mm(fig, dir_x_mm, TABLE_TOP_MM, dir_w_mm, 16 * ROW_HEIGHT_MM)
        ax_dir.set_ylim(matrix_bottom, matrix_top)
        ax_dir.axis("off")
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
                    facecolor=mcolors.to_rgba(color, DIRECTION_ALPHA),
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
    else:
        ax_dir = None

    # Family IDs
    ax_ids = add_axes_mm(fig, id_x_mm, TABLE_TOP_MM, id_w_mm, 16 * ROW_HEIGHT_MM)
    ax_ids.set_ylim(matrix_bottom, matrix_top)
    ax_ids.axis("off")
    ax_ids.set_xlim(0, 1)
    for row in geometry.itertuples(index=False):
        ax_ids.text(
            0.96,
            float(row.row_center),
            f"F{int(row.family_id)}",
            ha="right",
            va="center",
            fontsize=6.5,
            color=FG,
        )

    # Forest plot axes (single axes with visible x-axis)
    ax_forest = add_axes_mm(fig, forest_x_mm, TABLE_TOP_MM, forest_w_mm, 16 * ROW_HEIGHT_MM)
    ax_forest.set_xlim(-x_limit, x_limit)
    ax_forest.set_ylim(matrix_bottom, matrix_top)
    ax_forest.spines["top"].set_visible(False)
    ax_forest.spines["right"].set_visible(False)
    ax_forest.spines["left"].set_visible(False)
    ax_forest.spines["bottom"].set_visible(True)
    ax_forest.spines["bottom"].set_linewidth(0.5)
    ax_forest.spines["bottom"].set_color(FG)
    ax_forest.tick_params(axis="y", left=False, labelleft=False)
    ax_forest.tick_params(axis="x", direction="out", length=2, width=0.5, colors=FG, labelsize=6.5)
    if abs(x_limit - 4.0) < 1e-9:
        xticks = np.array([-4.0, -2.0, 0.0, 2.0, 4.0])
    elif abs(x_limit - 2.5) < 1e-9:
        xticks = np.array([-2.5, -1.0, 0.0, 1.0, 2.5])
    else:
        xticks = np.linspace(-x_limit, x_limit, 5)
    ax_forest.set_xticks(xticks)
    ax_forest.set_xticklabels([f"{x:g}" for x in xticks])
    ax_forest.axvline(0, color="#888888", linewidth=0.8, zorder=1)

    # Cohort-marker legend in the upper forest region (one clean line)
    ax_legend_top = add_axes_mm(fig, 2.5, 14.5, 83.0, 5.0)
    ax_legend_top.set_xlim(0, 1)
    ax_legend_top.set_ylim(0, 1)
    ax_legend_top.axis("off")
    legend_items = [
        ("GSE128188", "o"),
        ("GSE79768", "^"),
        ("GSE115574", "s"),
        ("pooled", "D"),
    ]
    x_positions = [0.06, 0.30, 0.55, 0.79]
    for x, (label, marker) in zip(x_positions, legend_items):
        if marker == "D":
            ax_legend_top.plot(
                x,
                0.45,
                marker=marker,
                markersize=3.4,
                markerfacecolor=MARKER_COLOR,
                markeredgecolor=MARKER_COLOR,
                markeredgewidth=0.35,
                zorder=4,
            )
        else:
            cohort_color = COHORT_COLORS[label]
            ax_legend_top.plot(
                x,
                0.45,
                marker=marker,
                markersize=2.4,
                markerfacecolor="white",
                markeredgecolor=cohort_color,
                markeredgewidth=0.8,
                zorder=4,
            )
        ax_legend_top.text(x + 0.03, 0.45, label, ha="left", va="center", fontsize=6.5, color=FG)

    # Four vertically dodged bands within each taller row; cohort CIs restored.
    usable = ROW_HEIGHT_MM * 0.82
    step = usable / 3.0
    offsets = {
        "GSE128188": -1.5 * step,
        "GSE79768": -0.5 * step,
        "GSE115574": 0.5 * step,
        "POOLED": 1.5 * step,
    }
    ci_lw_cohort = 0.55
    ci_lw_pooled = 1.55
    cohort_marker_size = 2.2
    pooled_marker_size = 3.6

    for family_id in EXPECTED_FAMILY_ORDER:
        y_center = float(row_lookup.loc[family_id, "row_center"])
        family_effects = effects.loc[effects["family_id"] == family_id]

        for cohort in COHORT_ORDER:
            est = family_effects.loc[family_effects["cohort"] == cohort].iloc[0]
            y = y_center + offsets[cohort]
            cohort_color = COHORT_COLORS[cohort]
            ax_forest.plot(
                [est["ci_lower"], est["ci_upper"]],
                [y, y],
                color=cohort_color,
                linewidth=ci_lw_cohort,
                solid_capstyle="butt",
                alpha=0.70,
                zorder=2,
            )
            ax_forest.plot(
                est["hedges_g"],
                y,
                marker=COHORT_MARKERS[cohort],
                markersize=cohort_marker_size,
                markerfacecolor="white",
                markeredgecolor=cohort_color,
                markeredgewidth=0.85,
                zorder=4,
            )

        pooled = family_effects.loc[family_effects["cohort"] == "POOLED"].iloc[0]
        y_pooled = y_center + offsets["POOLED"]
        ax_forest.plot(
            [pooled["ci_lower"], pooled["ci_upper"]],
            [y_pooled, y_pooled],
            color=MARKER_COLOR,
            linewidth=ci_lw_pooled,
            solid_capstyle="butt",
            zorder=5,
        )
        ax_forest.plot(
            pooled["hedges_g"],
            y_pooled,
            marker="D",
            markersize=pooled_marker_size,
            markerfacecolor=MARKER_COLOR,
            markeredgecolor=MARKER_COLOR,
            markeredgewidth=0.35,
            zorder=6,
        )

    # Right-side column headers (edge-aligned to avoid overlap in narrow column)
    fig.text(
        (col_x_mm + 0.5) / PANEL_WIDTH_MM,
        1 - (HEADER_TOP_MM + HEADER_HEIGHT_MM / 2) / PANEL_HEIGHT_MM,
        "same sign",
        ha="left",
        va="center",
        fontsize=5.0,
        color=MUTED,
    )
    fig.text(
        (col_x_mm + col_w_mm - 0.5) / PANEL_WIDTH_MM,
        1 - (HEADER_TOP_MM + HEADER_HEIGHT_MM / 2) / PANEL_HEIGHT_MM,
        "consistency",
        ha="right",
        va="center",
        fontsize=5.0,
        color=MUTED,
    )

    # Right-side columns: same sign and consistency
    ax_cols = add_axes_mm(fig, col_x_mm, TABLE_TOP_MM, col_w_mm, 16 * ROW_HEIGHT_MM)
    ax_cols.set_xlim(0, 1)
    ax_cols.set_ylim(matrix_bottom, matrix_top)
    ax_cols.axis("off")

    for family_id in EXPECTED_FAMILY_ORDER:
        y = float(row_lookup.loc[family_id, "row_center"])
        s = summary_lookup.loc[family_id]
        same_sign_text = f"{int(s['n_sign_match'])}/3"
        consistency_class = str(s["consistency_class"])
        color = CONSISTENCY_COLORS[consistency_class]
        symbol = CONSISTENCY_SYMBOLS[consistency_class]
        ax_cols.text(0.25, y, same_sign_text, ha="center", va="center", fontsize=6.5, color=FG)
        ax_cols.text(0.58, y, symbol, ha="center", va="center", fontsize=7.0, color=color, fontweight="bold")
        ax_cols.text(0.84, y, consistency_class[0], ha="center", va="center", fontsize=6.5, color=color, fontweight="bold")

    # Row dividers and separator (light ordinary lines, strong separator)
    ordinary_boundaries = (geometry["row_center"] + geometry["row_height"] / 2).iloc[:-1]
    axes_to_line = [ax_ids, ax_forest, ax_cols]
    if ax_dir is not None:
        axes_to_line.append(ax_dir)
    for ax in axes_to_line:
        for boundary in ordinary_boundaries:
            if not np.isclose(float(boundary), separator_y):
                ax.axhline(boundary, color="#F5F5F5", lw=0.2, zorder=5)
        ax.axhline(separator_y, color="#454545", lw=1.0, zorder=6)

    # X-axis label moved below ticks; SMD (Hedges' g).
    fig.text(
        (forest_x_mm + forest_w_mm / 2) / PANEL_WIDTH_MM,
        1 - (TABLE_BOTTOM_MM + 5.5) / PANEL_HEIGHT_MM,
        "Hedges' g (SMD)",
        ha="center",
        va="top",
        fontsize=6.5,
        color=FG,
    )
    fig.text(
        (forest_x_mm + forest_w_mm / 2) / PANEL_WIDTH_MM,
        1 - (TABLE_BOTTOM_MM + 9.0) / PANEL_HEIGHT_MM,
        "aligned family aCAIrel (AF − SR)",
        ha="center",
        va="top",
        fontsize=6.0,
        color=MUTED,
    )

    # Consistency legend below the axis labels
    class_counts = summary["consistency_class"].value_counts().to_dict()
    ax_cons_legend = add_axes_mm(fig, 5.0, PANEL_HEIGHT_MM - 5.2, 78.0, 4.0)
    ax_cons_legend.set_xlim(0, 1)
    ax_cons_legend.set_ylim(0, 1)
    ax_cons_legend.axis("off")
    x_positions = [0.10, 0.46, 0.82]
    for x, cls in zip(x_positions, ["High", "Moderate", "Low"]):
        count = int(class_counts.get(cls, 0))
        ax_cons_legend.text(
            x + 0.02,
            0.55,
            CONSISTENCY_SYMBOLS[cls],
            ha="center",
            va="center",
            fontsize=6.5,
            color=CONSISTENCY_COLORS[cls],
            fontweight="bold",
        )
        ax_cons_legend.text(
            x + 0.05,
            0.55,
            f"{cls[0]} {cls} (n = {count})",
            ha="left",
            va="center",
            fontsize=6.0,
            color=FG,
        )

    return fig


def save_figure(fig: plt.Figure, base: Path, preview: Path | None = None) -> dict[str, Path]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    outputs = {
        "pdf": base.with_suffix(".pdf"),
        "svg": base.with_suffix(".svg"),
        "png_600dpi": base.with_suffix(".png"),
    }
    if preview is not None:
        outputs["preview"] = preview
    fig.savefig(outputs["pdf"])
    fig.savefig(outputs["svg"])
    fig.savefig(outputs["png_600dpi"], dpi=600)
    if preview is not None:
        fig.savefig(preview, dpi=200)
    plt.close(fig)
    return outputs


def validate_outputs(outputs: dict[str, Path]) -> None:
    pdfinfo = shutil.which("pdfinfo")
    pdftotext = shutil.which("pdftotext")
    require(pdfinfo is not None and pdftotext is not None, "PDF QC tools are unavailable")
    info = subprocess.run(
        [pdfinfo, str(outputs["pdf"])], check=True, capture_output=True, text=True
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read Figure 2D PDF dimensions")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    require(
        abs(width_mm - PANEL_WIDTH_MM) < 0.05 and abs(height_mm - PANEL_HEIGHT_MM) < 0.05,
        f"Figure 2D PDF dimensions are {width_mm:.2f} × {height_mm:.2f} mm",
    )
    visible = subprocess.run(
        [pdftotext, "-layout", str(outputs["pdf"]), "-"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    normalized = " ".join(visible.split())
    required = [
        "Cross-cohort AF-SR effects and consistency",
        "Cohort-specific and pooled Hedges' g across 16 frozen NCP families",
        "Hedges' g (SMD)",
        "aligned family aCAIrel (AF",
        "GSE128188",
        "GSE79768",
        "GSE115574",
        "consistency",
        "pooled",
    ]
    require(all(phrase in normalized for phrase in required), "Figure 2D visible text is incomplete")
    # "same sign" may be split across lines by pdftotext; check both words are present.
    require("same" in normalized and "sign" in normalized, "Right-side 'same sign' header is missing")
    # Ensure GSE115574 is not truncated in the visible text stream.
    require("GSE115574" in visible, "GSE115574 accession is missing or truncated in extracted text")
    forbidden = [
        "REML",
        "Hartung",
        "Knapp",
        "fixed-effect_inverse_variance",
    ]
    require(not any(token in visible for token in forbidden), "Prohibited model terminology visible in panel")

    require(outputs["svg"].is_file() and outputs["svg"].stat().st_size > 0, "SVG export is missing")
    expected_width = round(PANEL_WIDTH_MM / 25.4 * 600)
    expected_height = round(PANEL_HEIGHT_MM / 25.4 * 600)
    image = plt.imread(outputs["png_600dpi"])
    require(
        abs(image.shape[1] - expected_width) <= 1 and abs(image.shape[0] - expected_height) <= 1,
        f"600-dpi PNG dimensions changed: {image.shape[1]} × {image.shape[0]}",
    )


def write_qc(
    manifest: pd.DataFrame,
    paths: dict[str, Path],
    effects: pd.DataFrame,
    summary: pd.DataFrame,
    audit: pd.DataFrame,
    outputs: dict[str, Path],
    composite_outputs: dict[str, Path],
) -> None:
    input_lines = [
        f"- `{path.relative_to(ROOT)}`: `{sha256(path)}`"
        for path in sorted(paths.values(), key=lambda item: str(item))
    ]
    output_paths = {
        "source": OUT_SOURCE,
        "summary": OUT_SUMMARY,
        "audit": OUT_AUDIT,
        **outputs,
        **{f"composite_{k}": v for k, v in composite_outputs.items()},
    }
    output_lines = [
        f"- `{path.relative_to(ROOT)}`: `{sha256(path)}`"
        for path in output_paths.values()
    ]

    max_all_ci = float(effects[["ci_lower", "ci_upper"]].abs().max().max())
    x_limit = plotted_axis_bound(effects)
    max_drawn = max_all_ci

    class_counts = summary["consistency_class"].value_counts().to_dict()
    class_lines = [
        f"- {cls}: {int(class_counts.get(cls, 0))} families"
        for cls in ["High", "Moderate", "Low"]
    ]

    text = "\n".join(
        [
            "# Figure 2D corrected family effects and consistency QC",
            "",
            "Status: **PASS**",
            "",
            "## Visual design",
            "- Title area compressed to title + subtitle only.",
            "- Explanatory note ('Consistency combines pooled-CI support and cohort sign agreement; descriptive framework summary.') is reserved for the figure caption, not shown in the panel body.",
            "- Cohort-marker legend is placed in the upper forest region; consistency class legend is at the bottom.",
            "- Pooled estimates are emphasised with a larger filled diamond and thicker CI line; cohort-specific estimates remain open and lighter.",
            "- Ordinary row separators are lightened; LA/RA separator remains strong.",
            "- Two panel variants are exported: standalone (with LA/RA direction strip) and composite (direction strip removed for use beside Figure 2C).",
            "- Row height is unlocked from Fig 2C mm lock (taller rows) so cohort CIs can be drawn without stacking; family order and LA/RA block membership remain locked to Fig 2C.",
            "- Within each family row, three cohort markers + thin CIs and the pooled diamond + thick CI occupy four vertically dodged bands.",
            "- X-axis is scaled to the full cohort+pooled CI envelope; effect measure is labelled as Hedges' g (SMD).",
            "",
            "## Pooling model audit",
            "- Family-level pooling model: **fixed-effect inverse-variance weighting**.",
            "- Weight formula: w_i = 1 / Var(g_i).",
            "- Pooled variance: 1 / Σ w_i.",
            "- 95% CI: Hedges' g ± 1.96 × SE (normal approximation).",
            "- tau² estimator: **not applicable** (no between-study variance estimated).",
            "- I² implementation: Higgins–Thompson on fixed-effect Q statistic, max(0, (Q - df) / Q).",
            "- Pooled estimates verified to come from the locked model (`pooling_model == 'fixed-effect_inverse_variance'`; `tau2 == 'NA'`).",
            "- No random-effects / REML / Hartung–Knapp model was used for family-level pooling.",
            "- No legacy effect or consistency table was read as a plotting input.",
            "",
            "## Classification",
            "- High: pooled 95% CI excludes 0 AND n_sign_match >= 2.",
            "- Moderate: not High AND n_sign_match >= 2.",
            "- Low: all remaining cases.",
            "",
            "## Runtime consistency classes",
            *class_lines,
            "",
            f"## Geometry",
            f"- Panel size: {PANEL_WIDTH_MM:.1f} × {PANEL_HEIGHT_MM:.1f} mm.",
            f"- Row height: {ROW_HEIGHT_MM:.2f} mm.",
            f"- Table top/bottom: {TABLE_TOP_MM:.2f} / {TABLE_BOTTOM_MM:.2f} mm.",
            f"- LA/RA separator: {SEPARATOR_MM:.2f} mm from panel top.",
            f"- X-axis symmetric bound: ±{x_limit:.2f} Hedges' g (from all drawn CI bounds; max |CI| = {max_drawn:.3f}).",
            "- Family order / LA–RA membership from shared `figure2_family_row_geometry.tsv`; physical row centers rebuilt via `expand_row_geometry()` and asserted against `figure2_family_row_geometry_composite.tsv`.",
            "- Forest encoding: cohort open markers + thin CI lines and pooled filled diamond + thick CI on four vertically dodged bands.",
            "- X-axis label: \"Hedges' g (SMD)\" with secondary line for aligned family aCAIrel (AF − SR), placed below the forest ticks.",
            f"- Composite panel contract: {PANEL_WIDTH_MM:.0f} × {PANEL_HEIGHT_MM:.0f} mm; TABLE_TOP={TABLE_TOP_MM:.1f}; ROW_HEIGHT={ROW_HEIGHT_MM:.2f} (shared with Figure 2C composite).",
            "",
            "## Migration audit",
            f"- Families with changed consistency class: {len(audit)}.",
            f"- Family IDs: {sorted(audit['family_id'].tolist())}.",
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


def main() -> None:
    manifest, paths = fig2b.load_manifest()
    donor_source, families, donor_stats = fig2b.build_source(paths)
    require(
        float(donor_stats["max_verification_error"]) <= 1e-12,
        "Donor-family validation failed upstream",
    )

    geometry = expand_row_geometry(fig2c.load_row_geometry())
    assert_matches_composite_geometry(geometry)

    effects = compute_effects(donor_source, families)
    validate_model_lock(effects)
    summary = compute_consistency(effects)
    audit = build_legacy_audit(effects, families)

    # Final source schema ordering and assertion
    source = effects[SOURCE_COLUMNS].copy()
    require(
        source.groupby("family_id").size().eq(4).all(),
        "Source table does not have 4 rows per family",
    )
    require(
        source.loc[source["cohort"] == "POOLED", "pooling_model"].eq("fixed-effect_inverse_variance").all(),
        "Pooled source rows have wrong pooling_model",
    )
    require(
        source.loc[source["cohort"] != "POOLED", "pooling_model"].eq("NA").all(),
        "Cohort-specific source rows have wrong pooling_model",
    )
    require(
        source["tau2"].eq("NA").all(),
        "tau2 column contains non-NA values",
    )
    require(
        source[["ci_lower", "ci_upper"]].notna().all().all(),
        "Source CIs contain missing values",
    )
    require(
        np.isfinite(source[["ci_lower", "ci_upper"]]).all().all(),
        "Source CIs contain non-finite values",
    )

    # Donor count assertions against Figs. 2A and 2B
    for cohort in COHORT_ORDER:
        cohort_rows = source.loc[source["cohort"] == cohort]
        require(cohort_rows["n_SR"].eq(EXPECTED_COUNTS[cohort]["SR"]).all(), f"{cohort} SR count changed")
        require(cohort_rows["n_AF"].eq(EXPECTED_COUNTS[cohort]["AF"]).all(), f"{cohort} AF count changed")

    # n_sign_match bounds
    require(
        source["n_sign_match"].dropna().astype(int).between(0, 3).all(),
        "n_sign_match outside 0-3",
    )

    # No state/chamber-decomposition field in classifier
    for col in ["state", "dominant_state", "quadrant", "driver", "delta_LA", "delta_RA"]:
        require(col not in source.columns, f"State/chamber-decomposition field '{col}' entered the classifier")

    # Consistency rule verification
    for _, row in summary.iterrows():
        if row["pooled_supported"] and row["n_sign_match"] >= 2:
            expected = "High"
        elif row["n_sign_match"] >= 2:
            expected = "Moderate"
        else:
            expected = "Low"
        require(row["consistency_class"] == expected, f"Consistency rule violated for family {row['family_id']}")

    RESULTS.mkdir(parents=True, exist_ok=True)
    source.to_csv(OUT_SOURCE, sep="\t", index=False, float_format="%.15g")
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False, float_format="%.15g")
    audit.to_csv(OUT_AUDIT, sep="\t", index=False, float_format="%.15g")

    # Standalone Figure 2D (with LA/RA direction strip)
    fig_standalone = render(source, summary, geometry, composite_mode=False)
    outputs = save_figure(fig_standalone, OUT_BASE, OUT_PREVIEW)
    validate_outputs(outputs)

    # Composite-panel version (direction strip removed, forest widened)
    composite_base = OUT_BASE.with_name(OUT_BASE.name + "_composite")
    composite_preview = OUT_PREVIEW.with_name(OUT_PREVIEW.stem + "_composite.png")
    fig_composite = render(source, summary, geometry, composite_mode=True)
    composite_outputs = save_figure(fig_composite, composite_base, composite_preview)
    validate_outputs(composite_outputs)

    write_qc(manifest, paths, source, summary, audit, outputs, composite_outputs)

    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)} ({len(source)} rows)")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)} ({len(summary)} rows)")
    print(f"Wrote {OUT_AUDIT.relative_to(ROOT)} ({len(audit)} rows)")
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {composite_base.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {composite_preview.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")
    print(
        "Runtime counts: "
        + ", ".join(f"{cls}={int(summary['consistency_class'].eq(cls).sum())}" for cls in ["High", "Moderate", "Low"])
    )


if __name__ == "__main__":
    main()
