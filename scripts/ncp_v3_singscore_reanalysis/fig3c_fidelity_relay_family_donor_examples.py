#!/usr/bin/env python3
"""Final Figure 3C: family-level donor examples of Fidelity and Relay."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"

BRIDGE = RESULTS / "fig3c_fidelity_relay_donor_source.tsv"
FIG3B_SOURCE = RESULTS / "fig3b_family_compartment_localization_source.tsv"
CLASS_TABLE = ROOT / "results" / "fig3_v3_carrier_fidelity_class.tsv"
HCA_TABLE = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
AUTHORITY_LOCK = ROOT / "results" / "fig3_publication_authority_lock.md"
SENSITIVITY_TABLE = ROOT / "results" / "fig3_boundary_permutation_sensitivity.tsv"

OUT_SUMMARY = RESULTS / "fig3c_fidelity_relay_summary.tsv"
OUT_QC = RESULTS / "fig3c_fidelity_relay_qc.md"
OUT_GEOMETRY = RESULTS / "figure3cd_case_study_geometry.tsv"
AXIS_CONTRACT = RESULTS / "figure3cd_case_study_axis_contract.tsv"
OUT_BASE = FIGURES / "fig3c_fidelity_relay_family_donor_examples"
OUT_PREVIEW = FIGURES / "fig3c_fidelity_relay_family_donor_examples_preview.png"

APPROVED_BRIDGE_SHA256 = "b2f4cadbc19e0501879f85dbaf5caf19d5c6f879038551ad67cc5d4cdbcebb7e"
APPROVED_GEOMETRY_SHA256 = "4197ec0db7e4b0e6296b0b3912ffa441ee8b3ee14d96dd1d6f9d2b9456ec3a40"
APPROVED_AXIS_CONTRACT_SHA256 = "061566fa14fd27513c6731162300a26db32ddba5dbda35a7d1a5a89eddc83cc6"
TARGETS = [(10, "fibroblast"), (11, "endothelial"), (11, "adipocyte")]
FIDELITY_COLOR = "#8FAF7E"
RELAY_COLOR = "#A898B8"
FG = "#1A1A1A"
MUTED = "#5D5D5D"
LIGHT = "#777777"
AF_COLOR = "#303030"
ZERO_COLOR = "#E3E3E3"
BADGE_WIDTH_MM = 12.2
BADGE_HEIGHT_MM = 4.2
POINT_SIZE = 8.55
MEAN_LINE_WIDTH = 0.90
ZERO_LINE_WIDTH = 0.50
ANNOTATION_FONT_PT = 5.3

MM = 1 / 25.4
WIDTH_MM = 88.0
HEIGHT_MM = 82.0
MIN_FONT_PT = 5.3

GEOMETRY_COLUMNS = [
    "record_type",
    "target_panel",
    "block_id",
    "family_id",
    "compartment",
    "left_mm",
    "top_mm",
    "width_mm",
    "height_mm",
    "separator_y_mm",
    "panel_width_mm",
    "panel_height_mm",
    "shared_y_min",
    "shared_y_max",
    "status",
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


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": MIN_FONT_PT,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def load_family_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    runtime_inputs = [BRIDGE, FIG3B_SOURCE, CLASS_TABLE, HCA_TABLE, AUTHORITY_LOCK]
    require(all(path.is_file() for path in runtime_inputs), "A family-level Figure 3C input is missing")
    require(sha256(BRIDGE) == APPROVED_BRIDGE_SHA256, "Frozen family-level donor bridge changed")
    require("LOCKED" in AUTHORITY_LOCK.read_text(encoding="utf-8"), "Figure 3 authority is not locked")
    require(
        not any(
            token in str(path).lower()
            for path in runtime_inputs
            for token in ("program", "gene", "pre_correction", "legacy")
        ),
        "Renderer runtime includes a forbidden program-, gene-, pre-correction, or legacy input",
    )
    bridge = pd.read_csv(BRIDGE, sep="\t")
    fig3b = pd.read_csv(FIG3B_SOURCE, sep="\t")
    classes = pd.read_csv(CLASS_TABLE, sep="\t")
    hca = pd.read_csv(HCA_TABLE, sep="\t")
    return bridge, fig3b, classes, hca


def load_sensitivity_qc(summary: pd.DataFrame) -> dict[str, float | int | str]:
    """Load family-level precision sensitivity for QC only, never for rendering."""
    require(SENSITIVITY_TABLE.is_file(), "Family-level permutation sensitivity table is missing")
    sensitivity = pd.read_csv(SENSITIVITY_TABLE, sep="\t")
    row = sensitivity.loc[
        sensitivity["family_id"].astype(int).eq(11)
        & sensitivity["compartment"].eq("adipocyte")
    ]
    require(len(row) == 1, "Expected one F11 Adipocyte sensitivity record")
    row = row.iloc[0]
    primary = summary.loc[
        summary["family_id"].eq(11) & summary["compartment"].eq("adipocyte")
    ].iloc[0]
    require(int(row["n_perm_primary"]) == 10000, "Sensitivity record primary count changed")
    require(int(row["n_perm_sensitivity"]) == 1_000_000, "Targeted sensitivity count changed")
    require(
        abs(float(row["p_perm_primary_10k"]) - float(primary["p_perm"])) <= 1e-12,
        "Sensitivity record no longer reproduces the locked primary p value",
    )
    require(
        str(row["publication_primary"]) == "p_perm_primary_10k",
        "Sensitivity result incorrectly replaced the publication-primary result",
    )
    return {
        "p_primary": float(row["p_perm_primary_10k"]),
        "n_primary": int(row["n_perm_primary"]),
        "p_sensitivity": float(row["p_perm_sensitivity_1e6"]),
        "n_sensitivity": int(row["n_perm_sensitivity"]),
        "note": str(row["note"]),
    }


def build_summary(
    bridge: pd.DataFrame,
    fig3b: pd.DataFrame,
    classes: pd.DataFrame,
    hca: pd.DataFrame,
) -> pd.DataFrame:
    require(len(bridge) == 99, "Family-level donor bridge is not 99 rows")
    provenance_columns = [
        "source_donor_scores_sha256",
        "source_family_delta_sha256",
        "source_class_table_sha256",
        "source_hca_table_sha256",
        "source_fig3b_sha256",
        "source_authority_lock_sha256",
    ]
    require(all(column in bridge for column in provenance_columns), "Bridge provenance hashes are missing")
    require(
        all(bridge[column].nunique() == 1 for column in provenance_columns),
        "Bridge provenance hashes are not frozen",
    )
    require(
        set(map(tuple, bridge[["family_id", "compartment"]].drop_duplicates().to_numpy()))
        == set(TARGETS),
        "Figure 3C facets changed",
    )
    class_lookup = dict(zip(classes["family_id"].astype(int), classes["fidelity_class"].astype(str)))
    carrier_lookup = dict(zip(hca["family_id"].astype(int), hca["LA_carrier"].astype(str)))
    support_lookup = {
        int(row.family_id): set(str(row.sig_compartments).split(","))
        if pd.notna(row.sig_compartments) and str(row.sig_compartments)
        else set()
        for row in classes.itertuples(index=False)
    }
    require(class_lookup[10] == "Fidelity" and class_lookup[11] == "Relay", "F10/F11 classes changed")
    require(carrier_lookup[10] == "fibroblast", "F10 carrier is not Fibroblast")
    require(carrier_lookup[11] == "endothelial", "F11 carrier is not Endothelial")
    require(support_lookup[10] == {"fibroblast"}, "F10 support is not Fibroblast only")
    require(support_lookup[11] == {"adipocyte"}, "F11 support is not Adipocyte only")

    rows = []
    for order, (family_id, compartment) in enumerate(TARGETS, start=1):
        donor = bridge.loc[
            bridge["family_id"].eq(family_id) & bridge["compartment"].eq(compartment)
        ].copy()
        b3 = fig3b.loc[
            fig3b["family_id"].eq(family_id) & fig3b["compartment"].eq(compartment)
        ]
        require(len(b3) == 1, f"Missing Figure 3B comparator for F{family_id} {compartment}")
        b3 = b3.iloc[0]
        counts = donor["group"].value_counts().to_dict()
        means = donor.groupby("group")["family_score"].mean().to_dict()
        delta = float(means["AF"] - means["Ctrl"])
        p_perm = float(donor["p_perm"].iloc[0])
        require(counts == {"AF": int(b3["n_af"]), "Ctrl": int(b3["n_ctrl"])}, "Donor counts do not reproduce Figure 3B")
        require(abs(delta - float(b3["dScore"])) <= 1e-9, "Donor delta does not reproduce Figure 3B")
        require(abs(p_perm - float(b3["p_perm"])) <= 1e-9, "p_perm does not reproduce Figure 3B")
        require(donor["n_perm"].eq(10000).all(), "Primary permutation count changed")
        require(donor["perm_strategy"].eq("mc").all(), "Primary permutation strategy changed")
        rows.append(
            {
                "facet_order": order,
                "family_id": family_id,
                "family_label": str(donor["family_label"].iloc[0]),
                "fidelity_class": class_lookup[family_id],
                "compartment": compartment,
                "compartment_label": str(donor["compartment_label"].iloc[0]),
                "healthy_carrier": carrier_lookup[family_id],
                "is_carrier": compartment == carrier_lookup[family_id],
                "is_supported": compartment in support_lookup[family_id],
                "n_ctrl": int(counts["Ctrl"]),
                "n_af": int(counts["AF"]),
                "mean_ctrl": float(means["Ctrl"]),
                "mean_af": float(means["AF"]),
                "delta": delta,
                "p_perm": p_perm,
                "n_perm": 10000,
                "perm_strategy": "mc",
            }
        )
    summary = pd.DataFrame(rows)
    f11_endo = summary.loc[
        summary["family_id"].eq(11) & summary["compartment"].eq("endothelial")
    ].iloc[0]
    f11_adipo = summary.loc[
        summary["family_id"].eq(11) & summary["compartment"].eq("adipocyte")
    ].iloc[0]
    require(not bool(f11_endo["is_supported"]), "F11 Endothelial incorrectly supported")
    require(bool(f11_adipo["is_supported"]) and not bool(f11_adipo["is_carrier"]), "F11 Adipocyte logic changed")
    require(abs(float(f11_adipo["p_perm"]) - 0.0431) < 5e-5, "F11 Adipocyte p_perm does not reproduce 0.0431")
    return summary


def shared_y_range(bridge: pd.DataFrame) -> tuple[float, float]:
    values = bridge["family_score"].to_numpy(dtype=float)
    span = float(values.max() - values.min())
    padding = max(0.02, span * 0.08)
    y_min = min(0.0, float(values.min()) - padding)
    y_max = float(values.max()) + padding
    require(y_min < y_max and y_min <= 0 <= y_max, "Shared y-axis does not include zero")
    return y_min, y_max


def build_geometry(y_min: float, y_max: float) -> pd.DataFrame:
    rows = [
        {
            "record_type": "panel",
            "target_panel": "Fig3C",
            "block_id": "",
            "family_id": np.nan,
            "compartment": "",
            "left_mm": 0.0,
            "top_mm": 0.0,
            "width_mm": WIDTH_MM,
            "height_mm": HEIGHT_MM,
            "separator_y_mm": 40.8,
            "panel_width_mm": WIDTH_MM,
            "panel_height_mm": HEIGHT_MM,
            "shared_y_min": y_min,
            "shared_y_max": y_max,
            "status": "active",
        },
        {
            "record_type": "block",
            "target_panel": "Fig3C",
            "block_id": "block1",
            "family_id": 10,
            "compartment": "",
            "left_mm": 2.5,
            "top_mm": 10.5,
            "width_mm": 83.0,
            "height_mm": 29.8,
            "separator_y_mm": 40.8,
            "panel_width_mm": WIDTH_MM,
            "panel_height_mm": HEIGHT_MM,
            "shared_y_min": y_min,
            "shared_y_max": y_max,
            "status": "active",
        },
        {
            "record_type": "facet",
            "target_panel": "Fig3C",
            "block_id": "block1",
            "family_id": 10,
            "compartment": "fibroblast",
            "left_mm": 8.0,
            "top_mm": 17.0,
            "width_mm": 35.0,
            "height_mm": 13.5,
            "separator_y_mm": 40.8,
            "panel_width_mm": WIDTH_MM,
            "panel_height_mm": HEIGHT_MM,
            "shared_y_min": y_min,
            "shared_y_max": y_max,
            "status": "active",
        },
        {
            "record_type": "logic",
            "target_panel": "Fig3C",
            "block_id": "block1",
            "family_id": 10,
            "compartment": "",
            "left_mm": 47.0,
            "top_mm": 17.0,
            "width_mm": 38.0,
            "height_mm": 20.5,
            "separator_y_mm": 40.8,
            "panel_width_mm": WIDTH_MM,
            "panel_height_mm": HEIGHT_MM,
            "shared_y_min": y_min,
            "shared_y_max": y_max,
            "status": "active",
        },
        {
            "record_type": "block",
            "target_panel": "Fig3C",
            "block_id": "block2",
            "family_id": 11,
            "compartment": "",
            "left_mm": 2.5,
            "top_mm": 41.3,
            "width_mm": 83.0,
            "height_mm": 31.7,
            "separator_y_mm": 40.8,
            "panel_width_mm": WIDTH_MM,
            "panel_height_mm": HEIGHT_MM,
            "shared_y_min": y_min,
            "shared_y_max": y_max,
            "status": "active",
        },
        {
            "record_type": "facet",
            "target_panel": "Fig3C",
            "block_id": "block2",
            "family_id": 11,
            "compartment": "endothelial",
            "left_mm": 8.0,
            "top_mm": 48.0,
            "width_mm": 21.5,
            "height_mm": 12.0,
            "separator_y_mm": 40.8,
            "panel_width_mm": WIDTH_MM,
            "panel_height_mm": HEIGHT_MM,
            "shared_y_min": y_min,
            "shared_y_max": y_max,
            "status": "active",
        },
        {
            "record_type": "facet",
            "target_panel": "Fig3C",
            "block_id": "block2",
            "family_id": 11,
            "compartment": "adipocyte",
            "left_mm": 31.5,
            "top_mm": 48.0,
            "width_mm": 22.5,
            "height_mm": 12.0,
            "separator_y_mm": 40.8,
            "panel_width_mm": WIDTH_MM,
            "panel_height_mm": HEIGHT_MM,
            "shared_y_min": y_min,
            "shared_y_max": y_max,
            "status": "active",
        },
        {
            "record_type": "logic",
            "target_panel": "Fig3C",
            "block_id": "block2",
            "family_id": 11,
            "compartment": "",
            "left_mm": 57.0,
            "top_mm": 48.0,
            "width_mm": 28.0,
            "height_mm": 22.0,
            "separator_y_mm": 40.8,
            "panel_width_mm": WIDTH_MM,
            "panel_height_mm": HEIGHT_MM,
            "shared_y_min": y_min,
            "shared_y_max": y_max,
            "status": "active",
        },
        {
            "record_type": "footer",
            "target_panel": "Fig3C",
            "block_id": "",
            "family_id": np.nan,
            "compartment": "",
            "left_mm": 2.5,
            "top_mm": 74.0,
            "width_mm": 83.0,
            "height_mm": 6.0,
            "separator_y_mm": 40.8,
            "panel_width_mm": WIDTH_MM,
            "panel_height_mm": HEIGHT_MM,
            "shared_y_min": y_min,
            "shared_y_max": y_max,
            "status": "active",
        },
    ]
    for block_id, top_mm, height_mm in [("block1", 10.5, 29.8), ("block2", 41.3, 31.7)]:
        rows.append(
            {
                "record_type": "reserved_block",
                "target_panel": "Fig3D",
                "block_id": block_id,
                "family_id": np.nan,
                "compartment": "",
                "left_mm": 2.5,
                "top_mm": top_mm,
                "width_mm": 83.0,
                "height_mm": height_mm,
                "separator_y_mm": 40.8,
                "panel_width_mm": WIDTH_MM,
                "panel_height_mm": HEIGHT_MM,
                "shared_y_min": y_min,
                "shared_y_max": y_max,
                "status": "reserved",
            }
        )
    geometry = pd.DataFrame(rows)[GEOMETRY_COLUMNS]
    require(len(geometry) == 11, "Unexpected case-study geometry row count")
    return geometry


def load_geometry() -> pd.DataFrame:
    require(OUT_GEOMETRY.is_file(), "Shared Figure 3C/3D geometry is missing")
    require(
        sha256(OUT_GEOMETRY) == APPROVED_GEOMETRY_SHA256,
        "Shared Figure 3C/3D geometry contract changed",
    )
    geometry = pd.read_csv(OUT_GEOMETRY, sep="\t")
    require(geometry.columns.tolist() == GEOMETRY_COLUMNS, "Shared case-study geometry schema changed")
    require(len(geometry) == 11, "Shared case-study geometry row count changed")
    require(
        geometry.loc[geometry["target_panel"].eq("Fig3D"), "status"].eq("reserved").all(),
        "Figure 3D reserved geometry changed",
    )
    return geometry


def load_axis_contract() -> pd.Series:
    require(AXIS_CONTRACT.is_file(), "Shared Figure 3C/3D axis contract is missing")
    require(
        sha256(AXIS_CONTRACT) == APPROVED_AXIS_CONTRACT_SHA256,
        "Shared Figure 3C/3D axis contract changed",
    )
    table = pd.read_csv(AXIS_CONTRACT, sep="\t")
    rows = table.loc[table["contract_id"].eq("figure3cd_v1")]
    require(len(rows) == 1, "Expected one standalone Figure 3C/3D axis contract")
    axis = rows.iloc[0]
    require(float(axis["point_size_pt2"]) == POINT_SIZE, "Point-size contract mismatch")
    require(float(axis["mean_line_width_pt"]) == MEAN_LINE_WIDTH, "Mean-line contract mismatch")
    require(str(axis["control_edge"]) == LIGHT, "Control-edge contract mismatch")
    require(str(axis["zero_line_color"]) == ZERO_COLOR, "Zero-line color contract mismatch")
    require(float(axis["zero_line_width_pt"]) == ZERO_LINE_WIDTH, "Zero-line contract mismatch")
    require(float(axis["badge_width_mm"]) == BADGE_WIDTH_MM, "Badge-width contract mismatch")
    require(float(axis["badge_height_mm"]) == BADGE_HEIGHT_MM, "Badge-height contract mismatch")
    require(int(axis["annotation_lines"]) == 2, "Annotation-line contract mismatch")
    require(float(axis["annotation_font_pt"]) == ANNOTATION_FONT_PT, "Annotation-font contract mismatch")
    return axis


def mm_axes(
    fig: plt.Figure,
    left: float,
    top: float,
    width: float,
    height: float,
) -> plt.Axes:
    return fig.add_axes(
        [
            left / WIDTH_MM,
            (HEIGHT_MM - top - height) / HEIGHT_MM,
            width / WIDTH_MM,
            height / HEIGHT_MM,
        ]
    )


def stable_jitter(family_id: int, compartment: str, group: str, n: int) -> np.ndarray:
    token = f"{family_id}|{compartment}|{group}".encode("utf-8")
    seed = int.from_bytes(hashlib.sha256(token).digest()[:4], "little")
    return np.random.default_rng(seed).uniform(-0.10, 0.10, size=n)


def format_delta(value: float) -> str:
    magnitude = f"{abs(value):.3f}"
    return f"−{magnitude}" if value < 0 else magnitude


def draw_facet(
    ax: plt.Axes,
    donor: pd.DataFrame,
    summary_row: pd.Series,
    y_min: float,
    y_max: float,
    show_y_labels: bool,
    role: str,
) -> None:
    ax.axhline(0, color=ZERO_COLOR, lw=ZERO_LINE_WIDTH, zorder=0)
    for group, x in [("Ctrl", 0.0), ("AF", 1.0)]:
        values = donor.loc[donor["group"].eq(group), "family_score"].to_numpy(dtype=float)
        jitter = stable_jitter(
            int(summary_row["family_id"]),
            str(summary_row["compartment"]),
            group,
            len(values),
        )
        if group == "Ctrl":
            ax.scatter(
                x + jitter,
                values,
                s=POINT_SIZE,
                facecolors="white",
                edgecolors=LIGHT,
                linewidths=0.6,
                zorder=3,
            )
            mean_color = "#777777"
        else:
            ax.scatter(
                x + jitter,
                values,
                s=POINT_SIZE,
                facecolors=AF_COLOR,
                edgecolors=AF_COLOR,
                linewidths=0.3,
                zorder=3,
            )
            mean_color = AF_COLOR
        ax.hlines(
            float(values.mean()),
            x - 0.20,
            x + 0.20,
            color=mean_color,
            linewidth=MEAN_LINE_WIDTH,
            zorder=4,
        )
    ax.set_xlim(-0.42, 1.42)
    ax.set_ylim(y_min, y_max)
    ax.set_yticks([tick for tick in (0.0, 0.2, 0.4) if y_min <= tick <= y_max])
    ax.set_xticks([0, 1], ["Control", "AF"])
    ax.tick_params(axis="both", labelsize=6.5, length=2, width=0.4, pad=1.5)
    if not show_y_labels:
        ax.tick_params(axis="y", labelleft=False)
    ax.text(
        0.5,
        1.10,
        str(summary_row["compartment_label"]),
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.5,
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.91,
        role,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.0,
        fontfamily="Roboto Condensed",
        color=MUTED,
    )
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.45)
        ax.spines[side].set_color("#777777")


def draw_logic(
    ax: plt.Axes,
    fields: list[tuple[str, str]],
    mode: str,
) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    label_y = [0.96, 0.69, 0.42]
    value_y = [0.84, 0.57, 0.30]
    for (label, value), y_label, y_value in zip(fields, label_y, value_y):
        ax.text(
            0.0,
            y_label,
            label,
            ha="left",
            va="top",
            fontsize=6.0,
            color=MUTED,
            fontweight="semibold",
        )
        ax.text(
            0.09,
            y_value,
            value,
            ha="left",
            va="top",
            fontsize=6.3,
            color=FG,
        )
    glyph_y = 0.055
    if mode == "Fidelity":
        ax.add_patch(
            mpatches.Rectangle(
                (0.02, glyph_y),
                0.065,
                0.085,
                facecolor="white",
                edgecolor="#555555",
                linewidth=0.45,
            )
        )
        ax.plot(0.0525, glyph_y + 0.0425, "o", ms=2.0, color=FG)
        ax.text(
            0.135,
            glyph_y + 0.0425,
            "→ Fidelity",
            ha="left",
            va="center",
            fontsize=6.2,
            fontweight="semibold",
        )
    else:
        ax.add_patch(
            mpatches.Rectangle(
                (0.02, glyph_y),
                0.065,
                0.085,
                facecolor="white",
                edgecolor="#555555",
                linewidth=0.45,
            )
        )
        ax.text(
            0.125,
            glyph_y + 0.0425,
            "+",
            ha="center",
            va="center",
            fontsize=6.0,
        )
        ax.plot(0.205, glyph_y + 0.0425, "o", ms=2.0, color=FG)
        ax.text(
            0.285,
            glyph_y + 0.0425,
            "→ Relay",
            ha="left",
            va="center",
            fontsize=6.2,
            fontweight="semibold",
        )


def render(
    bridge: pd.DataFrame,
    summary: pd.DataFrame,
    geometry: pd.DataFrame,
    *,
    composite: bool = False,
    panel_width_mm: float = WIDTH_MM,
    panel_height_mm: float = HEIGHT_MM,
    annotation_font_pt: float = ANNOTATION_FONT_PT,
) -> plt.Figure:
    set_style()
    y_min = float(geometry["shared_y_min"].iloc[0])
    y_max = float(geometry["shared_y_max"].iloc[0])
    width_mm = float(panel_width_mm)
    height_mm = float(panel_height_mm)
    fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
    def panel_axes(left: float, top: float, width: float, height: float) -> plt.Axes:
        return fig.add_axes(
            [
                left / width_mm,
                (height_mm - top - height) / height_mm,
                width / width_mm,
                height / height_mm,
            ]
        )

    fig.text(2.5 / width_mm, 1 - 2.3 / height_mm, "C", fontsize=11, fontweight="bold", ha="left", va="top", color=FG)
    fig.text(
        8.0 / width_mm,
        1 - 2.3 / height_mm,
        "Fidelity and Relay examples",
        fontsize=8.6,
        fontweight="bold",
        ha="left",
        va="top",
        color=FG,
    )
    fig.text(
        8.0 / width_mm,
        1 - 7.0 / height_mm,
        "F10 GAG disease / TGF-β and F11 BMP signalling",
        fontsize=6.5,
        ha="left",
        va="top",
        color=MUTED,
    )
    fig.text(
        0.012,
        0.47,
        "Family score (centered rank units)",
        fontsize=6.5,
        rotation=90,
        ha="center",
        va="center",
        color=FG,
    )

    fig.text(2.5 / width_mm, 1 - 10.4 / height_mm, "F10  GAG disease / TGF-β", fontsize=7.0, fontweight="bold", ha="left", va="top")
    badge1 = panel_axes(85.0 - BADGE_WIDTH_MM, 10.7, BADGE_WIDTH_MM, BADGE_HEIGHT_MM)
    badge1.set_gid("fidelity-badge-axis")
    badge1.axis("off")
    badge1_patch = mpatches.FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0.02", facecolor=FIDELITY_COLOR, edgecolor="none")
    badge1_patch.set_gid("fidelity-badge")
    badge1.add_patch(badge1_patch)
    badge1_text = badge1.text(0.5, 0.5, "Fidelity", ha="center", va="center", fontsize=6.5, fontweight="bold", color="white")
    badge1_text.set_gid("fidelity-badge-text")

    f10 = summary.loc[(summary["family_id"].eq(10)) & summary["compartment"].eq("fibroblast")].iloc[0]
    f10_donor = bridge.loc[(bridge["family_id"].eq(10)) & bridge["compartment"].eq("fibroblast")]
    f10_geo = geometry.loc[(geometry["record_type"].eq("facet")) & geometry["family_id"].eq(10)].iloc[0]
    ax_f10 = panel_axes(f10_geo.left_mm, f10_geo.top_mm, f10_geo.width_mm, f10_geo.height_mm)
    draw_facet(
        ax_f10,
        f10_donor,
        f10,
        y_min,
        y_max,
        True,
        "carrier · supported",
    )
    fig.text(
        8.0 / width_mm,
        1 - 34.8 / height_mm,
        (
            f"Ctrl/AF {int(f10.n_ctrl)}/{int(f10.n_af)} · Δ {format_delta(float(f10.delta))}\n"
            f"permutation p {float(f10.p_perm):.4f}"
            if composite
            else f"Ctrl/AF n = {int(f10.n_ctrl)}/{int(f10.n_af)} · "
            f"Δ = {format_delta(float(f10.delta))}\n"
            f"permutation p = {float(f10.p_perm):.4f}"
        ),
        fontsize=annotation_font_pt,
        fontfamily="Roboto Condensed",
        linespacing=1.05,
        ha="left",
        va="top",
        color=MUTED,
    )
    f10_logic_geo = geometry.loc[(geometry["record_type"].eq("logic")) & geometry["family_id"].eq(10)].iloc[0]
    ax_f10_logic = panel_axes(f10_logic_geo.left_mm, f10_logic_geo.top_mm, f10_logic_geo.width_mm, f10_logic_geo.height_mm)
    draw_logic(
        ax_f10_logic,
        [
            ("Healthy LA carrier", "Fibroblast"),
            ("Nominal support", "Fibroblast only"),
            ("Additional support", "None"),
        ],
        "Fidelity",
    )

    separator_y = float(geometry["separator_y_mm"].iloc[0])
    ax_rule = panel_axes(2.5, separator_y, width_mm - 5.0, 0.3)
    ax_rule.axhline(0.5, color="#C8C8C8", lw=0.55)
    ax_rule.axis("off")

    fig.text(2.5 / width_mm, 1 - 41.2 / height_mm, "F11  BMP signalling", fontsize=7.0, fontweight="bold", ha="left", va="top")
    badge2 = panel_axes(85.0 - BADGE_WIDTH_MM, 41.5, BADGE_WIDTH_MM, BADGE_HEIGHT_MM)
    badge2.set_gid("relay-badge-axis")
    badge2.axis("off")
    badge2_patch = mpatches.FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0.02", facecolor=RELAY_COLOR, edgecolor="none")
    badge2_patch.set_gid("relay-badge")
    badge2.add_patch(badge2_patch)
    badge2_text = badge2.text(0.5, 0.5, "Relay", ha="center", va="center", fontsize=6.5, fontweight="bold", color="white")
    badge2_text.set_gid("relay-badge-text")

    for compartment, show_y in [("endothelial", True), ("adipocyte", False)]:
        row = summary.loc[
            summary["family_id"].eq(11) & summary["compartment"].eq(compartment)
        ].iloc[0]
        donor = bridge.loc[
            bridge["family_id"].eq(11) & bridge["compartment"].eq(compartment)
        ]
        geo = geometry.loc[
            geometry["record_type"].eq("facet")
            & geometry["family_id"].eq(11)
            & geometry["compartment"].eq(compartment)
        ].iloc[0]
        ax = panel_axes(geo.left_mm, geo.top_mm, geo.width_mm, geo.height_mm)
        role = (
            "carrier · not supported"
            if compartment == "endothelial"
            else "non-carrier · supported"
        )
        draw_facet(ax, donor, row, y_min, y_max, show_y, role)
        fig.text(
            float(geo.left_mm) / width_mm,
            1 - 64.0 / height_mm,
            (
                f"Ctrl/AF {int(row.n_ctrl)}/{int(row.n_af)} · Δ {format_delta(float(row.delta))}\n"
                f"permutation p {float(row.p_perm):.4f}"
                if composite
                else f"Ctrl/AF n = {int(row.n_ctrl)}/{int(row.n_af)} · "
                f"Δ = {format_delta(float(row.delta))}\n"
                f"permutation p = {float(row.p_perm):.4f}"
            ),
            fontsize=annotation_font_pt,
            fontfamily="Roboto Condensed",
            linespacing=1.05,
            ha="left",
            va="top",
            color=MUTED,
        )

    f11_logic_geo = geometry.loc[(geometry["record_type"].eq("logic")) & geometry["family_id"].eq(11)].iloc[0]
    ax_f11_logic = panel_axes(f11_logic_geo.left_mm, f11_logic_geo.top_mm, f11_logic_geo.width_mm, f11_logic_geo.height_mm)
    draw_logic(
        ax_f11_logic,
        [
            ("Healthy LA carrier", "Endothelial"),
            ("Nominal support", "Adipocyte only"),
            ("Carrier compartment", "No nominal support"),
        ],
        "Relay",
    )

    if not composite:
        fig.text(
            2.5 / width_mm,
            1 - 74.2 / height_mm,
            "Donor×compartment family scores; nominal permutation p values are\n"
            "localization cues, not FDR-controlled discoveries.",
            fontsize=6.0,
            linespacing=1.05,
            ha="left",
            va="top",
            color=MUTED,
        )
    return fig


def validate_figure(fig: plt.Figure, summary: pd.DataFrame) -> dict[str, float | int]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [
        text
        for text in fig.findobj(mpl.text.Text)
        if text.get_visible() and text.get_text().strip()
    ]
    min_font = min(float(text.get_fontsize()) for text in texts)
    require(min_font >= MIN_FONT_PT, f"Visible text below {MIN_FONT_PT} pt")
    outside = []
    for text in texts:
        bbox = text.get_window_extent(renderer)
        if (
            bbox.x0 < fig.bbox.x0 - 1
            or bbox.y0 < fig.bbox.y0 - 1
            or bbox.x1 > fig.bbox.x1 + 1
            or bbox.y1 > fig.bbox.y1 + 1
        ):
            outside.append(text.get_text())
    require(not outside, f"Visible text outside canvas: {outside}")
    overlaps: list[tuple[str, str]] = []
    for index, first in enumerate(texts):
        first_box = first.get_window_extent(renderer)
        for second in texts[index + 1 :]:
            second_box = second.get_window_extent(renderer)
            overlap_width = max(
                0.0,
                min(first_box.x1, second_box.x1) - max(first_box.x0, second_box.x0),
            )
            overlap_height = max(
                0.0,
                min(first_box.y1, second_box.y1) - max(first_box.y0, second_box.y0),
            )
            if overlap_width * overlap_height > 1.0:
                overlaps.append((first.get_text(), second.get_text()))
    require(not overlaps, f"Visible text overlaps detected: {overlaps}")
    visible = "\n".join(text.get_text() for text in texts)
    forbidden = [
        "*",
        "q=",
        "q_perm",
        "p<0.10",
        "Path A",
        "bulk consistency",
        "duration class",
        "gene-level",
        "program-level",
        "p_perm",
        "carrier unsupported",
    ]
    require(not any(token in visible for token in forbidden), "Forbidden Figure 3C text is visible")
    require(
        "permutation p = 0.0431" in visible,
        "F11 Adipocyte primary permutation p is not visible",
    )
    require("carrier · supported" in visible, "F10 facet role is missing")
    require("carrier · not supported" in visible, "F11 Endothelial facet role is missing")
    require("non-carrier · supported" in visible, "F11 Adipocyte facet role is missing")
    require(len(summary) == 3, "Expected three facet annotations")
    return {
        "min_font_pt": min_font,
        "n_visible_text": len(texts),
        "n_bbox_violations": len(outside),
        "n_overlap_violations": len(overlaps),
    }


def save_outputs(fig: plt.Figure) -> dict[str, Path]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    outputs = {
        "pdf": OUT_BASE.with_suffix(".pdf"),
        "svg": OUT_BASE.with_suffix(".svg"),
        "png_600dpi": OUT_BASE.with_suffix(".png"),
        "preview": OUT_PREVIEW,
    }
    fig.savefig(outputs["pdf"], bbox_inches=None, pad_inches=0)
    fig.savefig(outputs["svg"], bbox_inches=None, pad_inches=0)
    fig.savefig(outputs["png_600dpi"], dpi=600, bbox_inches=None, pad_inches=0)
    fig.savefig(outputs["preview"], dpi=200, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    return outputs


def validate_exports(outputs: dict[str, Path]) -> None:
    pdftotext = shutil.which("pdftotext")
    pdfinfo = shutil.which("pdfinfo")
    require(pdftotext is not None and pdfinfo is not None, "PDF QC tools unavailable")
    info = subprocess.run(
        [pdfinfo, str(outputs["pdf"])], check=True, capture_output=True, text=True
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read Figure 3C PDF dimensions")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    require(
        abs(width_mm - WIDTH_MM) < 0.05 and abs(height_mm - HEIGHT_MM) < 0.05,
        f"Figure 3C dimensions changed: {width_mm:.2f} × {height_mm:.2f} mm",
    )
    visible = subprocess.run(
        [pdftotext, "-layout", str(outputs["pdf"]), "-"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    normalized = " ".join(visible.split())
    required = [
        "Fidelity and Relay examples",
        "F10 GAG disease / TGF-β and F11 BMP signalling",
        "Healthy LA carrier",
        "Fibroblast only",
        "Endothelial",
        "Adipocyte only",
        "Carrier compartment",
        "No nominal support",
        "carrier · supported",
        "carrier · not supported",
        "non-carrier · supported",
        "permutation p = 0.0431",
        "not FDR-controlled discoveries",
    ]
    require(all(phrase in normalized for phrase in required), "Figure 3C PDF text is incomplete")
    image = plt.imread(outputs["png_600dpi"])
    expected_width = round(WIDTH_MM / 25.4 * 600)
    expected_height = round(HEIGHT_MM / 25.4 * 600)
    require(
        abs(image.shape[1] - expected_width) <= 1
        and abs(image.shape[0] - expected_height) <= 1,
        "Figure 3C 600-dpi PNG dimensions changed",
    )
    require(outputs["svg"].is_file() and outputs["svg"].stat().st_size > 0, "SVG export missing")


def write_qc(
    bridge: pd.DataFrame,
    summary: pd.DataFrame,
    geometry: pd.DataFrame,
    sensitivity_qc: dict[str, float | int | str],
    layout_qc: dict[str, float | int],
    outputs: dict[str, Path],
) -> None:
    authoritative = bridge[
        [
            "family_id",
            "compartment",
            "mean_ctrl",
            "mean_af",
            "delta",
            "p_perm",
        ]
    ].drop_duplicates()
    comparison = summary.merge(
        authoritative,
        on=["family_id", "compartment"],
        suffixes=("_derived", "_authoritative"),
        validate="one_to_one",
    )
    reproduction_errors = {
        metric: float(
            np.max(
                np.abs(
                    comparison[f"{metric}_derived"].to_numpy(dtype=float)
                    - comparison[f"{metric}_authoritative"].to_numpy(dtype=float)
                )
            )
        )
        for metric in ("mean_ctrl", "mean_af", "delta", "p_perm")
    }
    provenance_columns = [
        "source_donor_scores_sha256",
        "source_family_delta_sha256",
        "source_class_table_sha256",
        "source_hca_table_sha256",
        "source_fig3b_sha256",
        "source_authority_lock_sha256",
    ]
    bridge_provenance = {
        column: str(bridge[column].iloc[0]) for column in provenance_columns
    }
    hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in [
            BRIDGE,
            FIG3B_SOURCE,
            CLASS_TABLE,
            HCA_TABLE,
            SENSITIVITY_TABLE,
            AXIS_CONTRACT,
            OUT_SUMMARY,
            OUT_GEOMETRY,
            *outputs.values(),
        ]
    }
    y_min = float(geometry["shared_y_min"].iloc[0])
    y_max = float(geometry["shared_y_max"].iloc[0])
    lines = [
        "# Figure 3C Fidelity/Relay family-donor QC",
        "",
        "Status: **PASS**",
        "",
        "## Scientific assertions",
        "- All plotted values originate from the SHA256-locked family-level donor bridge.",
        "- Renderer opens no gene-level or program-level input; the sensitivity table is family-level and QC-only.",
        "- F10 class = Fidelity; healthy LA carrier = Fibroblast; nominal support = Fibroblast only; additional support = none.",
        "- F11 class = Relay; healthy LA carrier = Endothelial; nominal support = Adipocyte only.",
        "- F11 Endothelial is not supported; F11 Adipocyte is supported.",
        "- All displayed tests use the primary 10,000-permutation family-level result.",
        f"- Family-level donor rows: {len(bridge)}",
        f"- Shared y-axis: [{y_min:.6g}, {y_max:.6g}]",
        f"- Maximum Ctrl-mean reproduction error: {reproduction_errors['mean_ctrl']:.3e}",
        f"- Maximum AF-mean reproduction error: {reproduction_errors['mean_af']:.3e}",
        f"- Maximum delta reproduction error: {reproduction_errors['delta']:.3e}",
        f"- Maximum p_perm reproduction error: {reproduction_errors['p_perm']:.3e}",
        "",
        "## Frozen bridge provenance",
        "- Builder: `scripts/ncp_v3_singscore_reanalysis/build_fig3c_family_donor_bridge.py`",
        *[
            f"- `{column}`: `{digest}`"
            for column, digest in bridge_provenance.items()
        ],
        "",
        "## Targeted precision sensitivity",
        (
            "- F11 Adipocyte primary: "
            f"p={float(sensitivity_qc['p_primary']):.10f}, "
            f"n_perm={int(sensitivity_qc['n_primary']):,}."
        ),
        (
            "- F11 Adipocyte targeted sensitivity: "
            f"p={float(sensitivity_qc['p_sensitivity']):.10f}, "
            f"n_perm={int(sensitivity_qc['n_sensitivity']):,}."
        ),
        f"- Interpretation: {sensitivity_qc['note']}.",
        "",
        "## Facet reproduction",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"- F{row.family_id} {row.compartment}: Ctrl n={row.n_ctrl}, AF n={row.n_af}, "
            f"mean Ctrl={row.mean_ctrl:.8f}, mean AF={row.mean_af:.8f}, "
            f"delta={row.delta:+.8f}, p_perm={row.p_perm:.8f}, n_perm={row.n_perm}"
        )
    lines.extend(
        [
            "",
            "## Layout/export QC",
            f"- Minimum visible font: {float(layout_qc['min_font_pt']):.2f} pt",
            f"- Visible-text bbox violations: {int(layout_qc['n_bbox_violations'])}",
            f"- Visible-text overlap violations: {int(layout_qc['n_overlap_violations'])}",
            "- Panel dimensions: 88 × 82 mm",
            "- PDF, SVG, 600-dpi PNG, and 200-dpi reduced preview validated.",
            "",
            "## SHA256",
            *[f"- `{path}`: `{digest}`" for path, digest in hashes.items()],
            "",
        ]
    )
    OUT_QC.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    bridge, fig3b, classes, hca = load_family_inputs()
    summary = build_summary(bridge, fig3b, classes, hca)
    sensitivity_qc = load_sensitivity_qc(summary)
    axis_contract = load_axis_contract()
    y_min, y_max = shared_y_range(bridge)
    RESULTS.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False, float_format="%.10g")
    geometry = load_geometry()
    require(
        abs(float(geometry["shared_y_min"].iloc[0]) - y_min) <= 5e-7
        and abs(float(geometry["shared_y_max"].iloc[0]) - y_max) <= 5e-7,
        "Locked shared y-axis range no longer reproduces the bridge",
    )
    require(
        abs(float(axis_contract["y_min"]) - float(geometry["shared_y_min"].iloc[0])) <= 5e-7
        and abs(float(axis_contract["y_max"]) - float(geometry["shared_y_max"].iloc[0])) <= 5e-7,
        "Shared Figure 3C/3D y-axis contract changed",
    )
    fig = render(bridge, summary, geometry)
    layout_qc = validate_figure(fig, summary)
    outputs = save_outputs(fig)
    validate_exports(outputs)
    write_qc(bridge, summary, geometry, sensitivity_qc, layout_qc, outputs)
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Validated {OUT_GEOMETRY.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
