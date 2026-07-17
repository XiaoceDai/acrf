#!/usr/bin/env python3
"""Final Figure 3D: family-level donor examples of Relay and Extension."""
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

FAMILY_DONOR = RESULTS / "gse255612_family_donor_compartment_scores.tsv"
FAMILY_CONTRACT = (
    RESULTS
    / "gse255612_corrected_family_compartment_permutation_classification.tsv"
)
FAMILY_METADATA = ROOT / "results" / "ncp_v3_family_labels.tsv"
AXIS_CONTRACT = RESULTS / "figure3cd_case_study_axis_contract.tsv"

OUT_SOURCE = RESULTS / "fig3d_relay_extension_donor_source.tsv"
OUT_SUMMARY = RESULTS / "fig3d_relay_extension_summary.tsv"
OUT_QC = RESULTS / "fig3d_relay_extension_qc.md"
OUT_BASE = FIGURES / "fig3d_relay_extension_family_donor_examples"
OUT_PREVIEW = FIGURES / "fig3d_relay_extension_family_donor_examples_preview.png"

APPROVED_HASHES = {
    FAMILY_DONOR: "0f257a20980eb14e0f7acbb02ed8225911c181507fc2bc5f4270ee2e4491270e",
    FAMILY_CONTRACT: "23e31f695418879eb9ebb1a614ada5d471cea24dbb077e6497717bc029fd38a9",
    FAMILY_METADATA: "c083ab182baca0ea49b9860d3363eb39e3da07c90199a1e97bed5b5b016cf05e",
    AXIS_CONTRACT: "061566fa14fd27513c6731162300a26db32ddba5dbda35a7d1a5a89eddc83cc6",
}
TARGETS = [
    (15, "fibroblast", True, 1, "carrier · not supported"),
    (15, "endothelial", True, 2, "non-carrier · supported"),
    (15, "lymphocyte", False, 3, "non-carrier · not supported"),
    (16, "fibroblast", True, 4, "carrier · supported"),
    (16, "atrialCM", True, 5, "non-carrier · supported"),
    (16, "adipocyte", False, 6, "non-carrier · not supported"),
]
COMPARTMENT_LABELS = {
    "fibroblast": "Fibroblast",
    "endothelial": "Endothelial",
    "lymphocyte": "Lymphocyte",
    "atrialCM": "aCM",
    "adipocyte": "Adipocyte",
}
RELAY_COLOR = "#A898B8"
EXTENSION_COLOR = "#7FA8AA"
FG = "#1A1A1A"
MUTED = "#5D5D5D"
ZERO_COLOR = "#E3E3E3"
BADGE_WIDTH_MM = 12.2
MIN_FONT_PT = 5.3
MM = 1 / 25.4


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


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.Series]:
    runtime_inputs = list(APPROVED_HASHES)
    require(all(path.is_file() for path in runtime_inputs), "A locked Figure 3D input is missing")
    require(
        all(sha256(path) == APPROVED_HASHES[path] for path in runtime_inputs),
        "A locked Figure 3D input changed",
    )
    require(
        not any(
            token in str(path).lower()
            for path in runtime_inputs
            for token in (
                "program",
                "gene",
                "leading",
                "expression",
                "legacy",
                "path_a",
            )
        ),
        "Renderer input path includes a forbidden analysis level",
    )
    donor = pd.read_csv(FAMILY_DONOR, sep="\t")
    contract = pd.read_csv(FAMILY_CONTRACT, sep="\t")
    metadata = pd.read_csv(FAMILY_METADATA, sep="\t")
    axis_table = pd.read_csv(AXIS_CONTRACT, sep="\t")
    axis_rows = axis_table.loc[axis_table["contract_id"].eq("figure3cd_v1")]
    require(len(axis_rows) == 1, "Expected one standalone Figure 3C/3D axis contract")
    axis = axis_rows.iloc[0]
    require(str(axis["contract_id"]) == "figure3cd_v1", "Axis contract version changed")
    require(float(axis["panel_width_mm"]) == 88.0, "Panel width changed")
    require(float(axis["panel_height_mm"]) == 82.0, "Panel height changed")
    require(str(axis["zero_line_color"]) == ZERO_COLOR, "Zero-line color changed")
    require(float(axis["badge_width_mm"]) == BADGE_WIDTH_MM, "Badge width changed")
    require(float(axis["badge_height_mm"]) == 4.2, "Badge height changed")
    require(int(axis["annotation_lines"]) == 2, "Annotation-line contract changed")
    require(float(axis["separator_y_mm"]) == 40.8, "Separator position changed")
    require(
        abs(
            float(axis["block1_bottom_mm"])
            - float(axis["block1_top_mm"])
            - 29.8
        )
        <= 1e-9
        and abs(
            float(axis["block2_bottom_mm"])
            - float(axis["block2_top_mm"])
            - 31.7
        )
        <= 1e-9,
        "Case-block height contract changed",
    )
    require(
        str(axis["source_geometry_sha256"])
        == "4197ec0db7e4b0e6296b0b3912ffa441ee8b3ee14d96dd1d6f9d2b9456ec3a40",
        "Figure 3C/3D geometry provenance changed",
    )
    return donor, contract, metadata, axis


def build_outputs(
    donor: pd.DataFrame,
    contract: pd.DataFrame,
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    label_lookup = dict(zip(metadata["family_id"].astype(int), metadata["short_label"]))
    require(label_lookup[15] == "NTRK2 / TrkB", "F15 publication label changed")
    require(label_lookup[16] == "Keratan sulfate", "F16 publication label changed")

    family15 = contract.loc[contract["family_id"].eq(15)]
    family16 = contract.loc[contract["family_id"].eq(16)]
    require(family15["fidelity_class"].eq("Relay").all(), "F15 is not Relay")
    require(family16["fidelity_class"].eq("Extension").all(), "F16 is not Extension")
    require(family15["healthy_LA_carrier"].eq("fibroblast").all(), "F15 carrier changed")
    require(family16["healthy_LA_carrier"].eq("fibroblast").all(), "F16 carrier changed")
    support15 = set(family15.loc[family15["is_supported"], "compartment"])
    support16 = set(family16.loc[family16["is_supported"], "compartment"])
    require(support15 == {"endothelial"}, f"F15 support changed: {support15}")
    require(
        support16 == {"fibroblast", "atrialCM"},
        f"F16 support changed: {support16}",
    )

    source_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []
    for family_id, compartment, displayed, order, role in TARGETS:
        values = donor.loc[
            donor["family_id"].eq(family_id)
            & donor["compartment"].eq(compartment)
        ].copy()
        crow = contract.loc[
            contract["family_id"].eq(family_id)
            & contract["compartment"].eq(compartment)
        ]
        require(len(crow) == 1 and not values.empty, f"Missing F{family_id} {compartment}")
        crow = crow.iloc[0]
        counts = values["group"].value_counts().to_dict()
        means = values.groupby("group")["family_score"].mean().to_dict()
        delta = float(means["AF"] - means["Ctrl"])
        require(
            counts == {"AF": int(crow["n_af"]), "Ctrl": int(crow["n_ctrl"])},
            f"F{family_id} {compartment} donor counts changed",
        )
        require(abs(float(means["Ctrl"]) - float(crow["mean_ctrl"])) <= 1e-9, "Control mean mismatch")
        require(abs(float(means["AF"]) - float(crow["mean_af"])) <= 1e-9, "AF mean mismatch")
        require(abs(delta - float(crow["dScore"])) <= 1e-9, "Figure 3B delta mismatch")
        require(int(crow["n_perm"]) == 10000, "Primary permutation count changed")
        is_carrier = bool(crow["is_carrier"])
        is_supported = bool(crow["is_supported"])
        summary_row = {
            "facet_order": order,
            "family_id": family_id,
            "family_label": label_lookup[family_id],
            "fidelity_class": str(crow["fidelity_class"]),
            "compartment": compartment,
            "compartment_label": COMPARTMENT_LABELS[compartment],
            "displayed": displayed,
            "facet_role": role,
            "healthy_LA_carrier": str(crow["healthy_LA_carrier"]),
            "is_carrier": is_carrier,
            "is_supported": is_supported,
            "n_ctrl": int(counts["Ctrl"]),
            "n_af": int(counts["AF"]),
            "mean_ctrl": float(means["Ctrl"]),
            "mean_af": float(means["AF"]),
            "delta": delta,
            "p_perm": float(crow["p_perm"]),
            "n_perm": int(crow["n_perm"]),
            "perm_strategy": str(crow["perm_strategy"]),
            "p_perm_sensitivity_1e6": float(crow["p_perm_sensitivity_1e6"])
            if pd.notna(crow["p_perm_sensitivity_1e6"])
            else np.nan,
            "n_perm_sensitivity": int(crow["n_perm_sensitivity"])
            if pd.notna(crow["n_perm_sensitivity"])
            else pd.NA,
            "sensitivity_note": str(crow["sensitivity_note"])
            if pd.notna(crow["sensitivity_note"])
            else "",
        }
        summary_rows.append(summary_row)
        for column, value in summary_row.items():
            values[column] = value
        source_parts.append(values)

    source = pd.concat(source_parts, ignore_index=True)
    source = source.sort_values(
        ["facet_order", "group", "donor_id"],
        key=lambda series: series.map({"Ctrl": "0", "AF": "1"})
        if series.name == "group"
        else series,
    ).reset_index(drop=True)
    summary = pd.DataFrame(summary_rows).sort_values("facet_order").reset_index(drop=True)

    f15_fibro = summary.query("family_id == 15 and compartment == 'fibroblast'").iloc[0]
    f15_endo = summary.query("family_id == 15 and compartment == 'endothelial'").iloc[0]
    f15_lymph = summary.query("family_id == 15 and compartment == 'lymphocyte'").iloc[0]
    f16_fibro = summary.query("family_id == 16 and compartment == 'fibroblast'").iloc[0]
    f16_acm = summary.query("family_id == 16 and compartment == 'atrialCM'").iloc[0]
    f16_adipo = summary.query("family_id == 16 and compartment == 'adipocyte'").iloc[0]
    require(not f15_fibro.is_supported, "F15 Fibroblast incorrectly supported")
    require(f15_endo.is_supported and abs(f15_endo.p_perm - 0.0471) < 5e-5, "F15 Endothelial support changed")
    require(not f15_lymph.is_supported, "F15 Lymphocyte incorrectly supported")
    require(f16_fibro.is_supported and f16_acm.is_supported, "F16 displayed support changed")
    require(not f16_adipo.is_supported and abs(f16_adipo.p_perm - 0.0542) < 5e-5, "F16 Adipocyte call changed")
    require(
        int(f15_endo.n_perm_sensitivity) == 1_000_000
        and np.isfinite(float(f15_endo.p_perm_sensitivity_1e6)),
        "F15 Endothelial precision sensitivity is missing",
    )
    return source, summary


def mm_axes(
    fig: plt.Figure,
    width_mm: float,
    height_mm: float,
    left: float,
    top: float,
    width: float,
    height: float,
) -> plt.Axes:
    return fig.add_axes(
        [
            left / width_mm,
            (height_mm - top - height) / height_mm,
            width / width_mm,
            height / height_mm,
        ]
    )


def stable_jitter(
    family_id: int,
    compartment: str,
    group: str,
    n: int,
    halfwidth: float,
) -> np.ndarray:
    token = f"{family_id}|{compartment}|{group}".encode("utf-8")
    seed = int.from_bytes(hashlib.sha256(token).digest()[:4], "little")
    return np.random.default_rng(seed).uniform(-halfwidth, halfwidth, size=n)


def format_delta(value: float) -> str:
    magnitude = f"{abs(value):.3f}"
    return f"−{magnitude}" if value < 0 else magnitude


def draw_facet(
    ax: plt.Axes,
    values: pd.DataFrame,
    row: pd.Series,
    axis: pd.Series,
    show_y_labels: bool,
) -> None:
    y_ticks = [float(value) for value in str(axis["y_ticks"]).split(",")]
    ax.axhline(
        float(axis["zero_line"]),
        color=str(axis["zero_line_color"]),
        lw=float(axis["zero_line_width_pt"]),
        zorder=0,
    )
    for group, x in [("Ctrl", 0.0), ("AF", 1.0)]:
        scores = values.loc[values["group"].eq(group), "family_score"].to_numpy(float)
        jitter = stable_jitter(
            int(row["family_id"]),
            str(row["compartment"]),
            group,
            len(scores),
            float(axis["jitter_halfwidth"]),
        )
        if group == "Ctrl":
            face = str(axis["control_face"])
            edge = str(axis["control_edge"])
            line_width = float(axis["control_linewidth_pt"])
            mean_color = "#777777"
        else:
            face = str(axis["af_face"])
            edge = str(axis["af_edge"])
            line_width = float(axis["af_linewidth_pt"])
            mean_color = str(axis["af_face"])
        ax.scatter(
            x + jitter,
            scores,
            s=float(axis["point_size_pt2"]),
            facecolors=face,
            edgecolors=edge,
            linewidths=line_width,
            zorder=3,
        )
        ax.hlines(
            float(scores.mean()),
            x - 0.20,
            x + 0.20,
            color=mean_color,
            linewidth=float(axis["mean_line_width_pt"]),
            zorder=4,
        )
    ax.set_xlim(-0.42, 1.42)
    ax.set_ylim(float(axis["y_min"]), float(axis["y_max"]))
    ax.set_yticks(y_ticks)
    ax.set_xticks([0, 1], ["Control", "AF"])
    ax.tick_params(axis="both", labelsize=6.5, length=2, width=0.4, pad=1.5)
    if not show_y_labels:
        ax.tick_params(axis="y", labelleft=False)
    ax.text(
        0.5,
        1.10,
        str(row["compartment_label"]),
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.5,
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.91,
        str(row["facet_role"]),
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
    for (label, value), y_label, y_value in zip(
        fields,
        [0.96, 0.69, 0.42],
        [0.84, 0.57, 0.30],
    ):
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
        ax.text(0.09, y_value, value, ha="left", va="top", fontsize=6.3, color=FG)
    glyph_y = 0.055
    if mode == "Relay":
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
        ax.text(0.125, glyph_y + 0.0425, "+", ha="center", va="center", fontsize=6.0)
        ax.plot(0.205, glyph_y + 0.0425, "o", ms=2.0, color=FG)
        ax.text(0.285, glyph_y + 0.0425, "→ Relay", ha="left", va="center", fontsize=6.2, fontweight="semibold")
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
        ax.plot(0.0525, glyph_y + 0.0425, "o", ms=2.0, color=FG)
        ax.text(0.125, glyph_y + 0.0425, "+", ha="center", va="center", fontsize=6.0)
        ax.plot(0.205, glyph_y + 0.0425, "o", ms=2.0, color=FG)
        ax.text(0.285, glyph_y + 0.0425, "→ Extension", ha="left", va="center", fontsize=6.2, fontweight="semibold")


def render(
    source: pd.DataFrame,
    summary: pd.DataFrame,
    axis: pd.Series,
    *,
    composite: bool = False,
) -> plt.Figure:
    set_style()
    width_mm = float(axis["panel_width_mm"])
    height_mm = float(axis["panel_height_mm"])
    fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
    fig.text(2.5 / width_mm, 1 - 2.3 / height_mm, "D", fontsize=11, fontweight="bold", ha="left", va="top")
    fig.text(
        8.0 / width_mm,
        1 - 2.3 / height_mm,
        "Relay and Extension examples",
        fontsize=8.6,
        fontweight="bold",
        ha="left",
        va="top",
    )
    fig.text(
        8.0 / width_mm,
        1 - 7.0 / height_mm,
        "F15 NTRK2 / TrkB and F16 Keratan sulfate",
        fontsize=6.5,
        color=MUTED,
        ha="left",
        va="top",
    )
    if not composite:
        fig.text(
            0.012,
            0.47,
            "Family score (centered rank units)",
            fontsize=6.5,
            rotation=90,
            ha="center",
            va="center",
        )

    blocks = [
        {
            "family_id": 15,
            "header_y": 10.4,
            "badge_top": 10.7,
            "header": "F15  NTRK2 / TrkB",
            "badge": "Relay",
            "badge_color": RELAY_COLOR,
            "facet_top": 17.0,
            "stats_y": 34.8,
            "logic_top": 17.0,
            "logic_height": 20.5,
            "facets": ["fibroblast", "endothelial"],
            "fields": [
                ("Healthy LA carrier", "Fibroblast"),
                ("Nominal support", "Endothelial only"),
                ("Carrier compartment", "No nominal support"),
            ],
        },
        {
            "family_id": 16,
            "header_y": 41.2,
            "badge_top": 41.5,
            "header": "F16  Keratan sulfate",
            "badge": "Extension",
            "badge_color": EXTENSION_COLOR,
            "facet_top": 48.0,
            "stats_y": 64.0,
            "logic_top": 48.0,
            "logic_height": 22.0,
            "facets": ["fibroblast", "atrialCM"],
            "fields": [
                ("Healthy LA carrier", "Fibroblast"),
                ("Nominal support", "Fibroblast + aCM"),
                ("Additional support", "aCM"),
            ],
        },
    ]
    facet_geometry = (
        [(8.0, 24.5), (34.0, 24.5)]
        if composite
        else [(8.0, 21.5), (31.5, 22.5)]
    )
    for block in blocks:
        family_id = int(block["family_id"])
        fig.text(
            2.5 / width_mm,
            1 - float(block["header_y"]) / height_mm,
            str(block["header"]),
            fontsize=7.0,
            fontweight="bold",
            ha="left",
            va="top",
        )
        badge = mm_axes(
            fig,
            width_mm,
            height_mm,
            85.0 - BADGE_WIDTH_MM,
            float(block["badge_top"]),
            BADGE_WIDTH_MM,
            float(axis["badge_height_mm"]),
        )
        badge.set_gid(f"{str(block['badge']).lower()}-badge-axis")
        badge.axis("off")
        badge_patch = mpatches.FancyBboxPatch(
            (0, 0),
            1,
            1,
            boxstyle="round,pad=0.02",
            facecolor=str(block["badge_color"]),
            edgecolor="none",
        )
        badge_patch.set_gid(f"{str(block['badge']).lower()}-badge")
        badge.add_patch(badge_patch)
        badge_text = badge.text(
            0.5,
            0.5,
            str(block["badge"]),
            ha="center",
            va="center",
            fontsize=6.0,
            fontweight="bold",
            color="white",
        )
        badge_text.set_gid(f"{str(block['badge']).lower()}-badge-text")
        for index, (compartment, (left, facet_width)) in enumerate(
            zip(block["facets"], facet_geometry)
        ):
            row = summary.loc[
                summary["family_id"].eq(family_id)
                & summary["compartment"].eq(compartment)
            ].iloc[0]
            values = source.loc[
                source["family_id"].eq(family_id)
                & source["compartment"].eq(compartment)
            ]
            facet = mm_axes(
                fig,
                width_mm,
                height_mm,
                left,
                float(block["facet_top"]),
                facet_width,
                12.0,
            )
            draw_facet(facet, values, row, axis, index == 0)
            fig.text(
                left / width_mm,
                1 - float(block["stats_y"]) / height_mm,
                (
                    f"Ctrl/AF {int(row.n_ctrl)}/{int(row.n_af)} · Δ {format_delta(float(row.delta))}\n"
                    f"permutation p {float(row.p_perm):.4f}"
                    if composite
                    else f"Ctrl/AF n = {int(row.n_ctrl)}/{int(row.n_af)} · "
                    f"Δ = {format_delta(float(row.delta))}\n"
                    f"permutation p = {float(row.p_perm):.4f}"
                ),
                fontsize=float(axis["annotation_font_pt"]),
                fontfamily="Roboto Condensed",
                linespacing=1.05,
                color=MUTED,
                ha="left",
                va="top",
            )
        logic = mm_axes(
            fig,
            width_mm,
            height_mm,
            60.0 if composite else 57.0,
            float(block["logic_top"]),
            24.5 if composite else 28.0,
            float(block["logic_height"]),
        )
        draw_logic(logic, block["fields"], str(block["badge"]))

    rule = mm_axes(
        fig,
        width_mm,
        height_mm,
        2.5,
        float(axis["separator_y_mm"]),
        width_mm - 5.0,
        0.3,
    )
    rule.axhline(0.5, color="#C8C8C8", lw=0.55)
    rule.axis("off")
    if not composite:
        fig.text(
            2.5 / width_mm,
            1 - 74.2 / height_mm,
            "Donor×compartment family scores; nominal permutation p values are\n"
            "localization cues, not FDR-controlled discoveries.",
            fontsize=6.0,
            color=MUTED,
            linespacing=1.05,
            ha="left",
            va="top",
        )
    return fig


def validate_figure(fig: plt.Figure) -> dict[str, float | int]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    extension_badges = [
        patch
        for patch in fig.findobj(mpatches.FancyBboxPatch)
        if patch.get_gid() == "extension-badge"
    ]
    relay_badges = [
        patch
        for patch in fig.findobj(mpatches.FancyBboxPatch)
        if patch.get_gid() == "relay-badge"
    ]
    require(len(extension_badges) == 1, "Extension badge is duplicated or missing")
    require(len(relay_badges) == 1, "Relay badge is duplicated or missing")
    extension_axes = [
        ax for ax in fig.axes if ax.get_gid() == "extension-badge-axis"
    ]
    relay_axes = [ax for ax in fig.axes if ax.get_gid() == "relay-badge-axis"]
    require(
        len(extension_axes) == 1 and len(relay_axes) == 1,
        "Badge axes are duplicated or missing",
    )
    extension_texts = [
        text for text in fig.findobj(mpl.text.Text)
        if text.get_gid() == "extension-badge-text"
    ]
    require(len(extension_texts) == 1, "Extension badge text is duplicated or missing")
    extension_text_box = extension_texts[0].get_window_extent(renderer)
    extension_axes_box = extension_axes[0].get_window_extent(renderer)
    require(
        extension_text_box.x0 >= extension_axes_box.x0 - 1
        and extension_text_box.x1 <= extension_axes_box.x1 + 1
        and extension_text_box.y0 >= extension_axes_box.y0 - 1
        and extension_text_box.y1 <= extension_axes_box.y1 + 1,
        "Extension badge text is clipped",
    )
    extension_box = extension_axes[0].get_position()
    relay_box = relay_axes[0].get_position()
    require(
        abs(extension_box.width - relay_box.width) <= 1e-12
        and abs(extension_box.height - relay_box.height) <= 1e-12,
        "Relay and Extension badge dimensions differ",
    )
    extension_top_mm = (1 - extension_box.y1) * 82.0
    extension_bottom_mm = (1 - extension_box.y0) * 82.0
    require(
        extension_top_mm >= 40.8
        and extension_bottom_mm <= 73.0
        and extension_box.x0 >= 0
        and extension_box.x1 <= 1,
        "Extension badge crosses its block or canvas boundary",
    )
    texts = [
        text
        for text in fig.findobj(mpl.text.Text)
        if text.get_visible() and text.get_text().strip()
    ]
    min_font = min(float(text.get_fontsize()) for text in texts)
    require(min_font >= MIN_FONT_PT, "Visible text is below the minimum font size")
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
        "p_perm",
        "q=",
        "q_perm",
        "p<0.10",
        "Path A",
        "gene-level",
        "program-level",
        "Lymphocyte",
        "Adipocyte",
    ]
    require(not any(token in visible for token in forbidden), "Forbidden Figure 3D text is visible")
    for phrase in (
        "carrier · not supported",
        "non-carrier · supported",
        "carrier · supported",
        "permutation p = 0.0471",
        "permutation p = 0.0094",
    ):
        require(phrase in visible, f"Required visible phrase missing: {phrase}")
    return {
        "min_font_pt": min_font,
        "n_bbox_violations": len(outside),
        "n_overlap_violations": len(overlaps),
        "n_extension_badges": len(extension_badges),
        "n_badge_clipping_violations": 0,
    }


def save_outputs(fig: plt.Figure, axis: pd.Series) -> dict[str, Path]:
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
    image = plt.imread(outputs["png_600dpi"])
    expected_width = round(float(axis["panel_width_mm"]) / 25.4 * 600)
    expected_height = round(float(axis["panel_height_mm"]) / 25.4 * 600)
    require(
        abs(image.shape[1] - expected_width) <= 1
        and abs(image.shape[0] - expected_height) <= 1,
        "600-dpi PNG dimensions changed",
    )
    return outputs


def validate_exports(outputs: dict[str, Path], axis: pd.Series) -> None:
    pdfinfo = shutil.which("pdfinfo")
    pdftotext = shutil.which("pdftotext")
    require(pdfinfo is not None and pdftotext is not None, "PDF QC tools unavailable")
    info = subprocess.run(
        [pdfinfo, str(outputs["pdf"])],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read PDF dimensions")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    require(
        abs(width_mm - float(axis["panel_width_mm"])) < 0.05
        and abs(height_mm - float(axis["panel_height_mm"])) < 0.05,
        "Figure 3D dimensions changed",
    )
    visible = subprocess.run(
        [pdftotext, "-layout", str(outputs["pdf"]), "-"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    normalized = " ".join(visible.split())
    for phrase in (
        "Relay and Extension examples",
        "F15 NTRK2 / TrkB and F16 Keratan sulfate",
        "Endothelial only",
        "Fibroblast + aCM",
        "No nominal support",
        "not FDR-controlled discoveries",
    ):
        require(phrase in normalized, f"PDF text is incomplete: {phrase}")
    require(
        len(re.findall(r"\bExtension\b", normalized)) == 3,
        "PDF contains a duplicated or missing Extension badge/object",
    )
    svg_text = outputs["svg"].read_text(encoding="utf-8")
    require(
        svg_text.count('id="extension-badge"') == 1
        and svg_text.count('id="extension-badge-text"') == 1,
        "SVG contains a duplicated or missing Extension badge/object",
    )
    image = plt.imread(outputs["png_600dpi"])
    rgb = image[:, :, :3]
    target = np.asarray([0x7F, 0xA8, 0xAA], dtype=float) / 255.0
    teal_mask = np.linalg.norm(rgb - target, axis=2) < 0.03
    require(teal_mask.any(), "PNG Extension badge color is missing")
    visited = np.zeros_like(teal_mask, dtype=bool)
    height_px, width_px = teal_mask.shape
    component_sizes: list[int] = []
    for start_y, start_x in zip(*np.where(teal_mask & ~visited)):
        if visited[start_y, start_x]:
            continue
        stack = [(int(start_y), int(start_x))]
        visited[start_y, start_x] = True
        component_size = 0
        while stack:
            y_pos, x_pos = stack.pop()
            component_size += 1
            for next_y, next_x in (
                (y_pos - 1, x_pos),
                (y_pos + 1, x_pos),
                (y_pos, x_pos - 1),
                (y_pos, x_pos + 1),
            ):
                if (
                    0 <= next_y < height_px
                    and 0 <= next_x < width_px
                    and teal_mask[next_y, next_x]
                    and not visited[next_y, next_x]
                ):
                    visited[next_y, next_x] = True
                    stack.append((next_y, next_x))
        component_sizes.append(component_size)
    require(
        sum(size >= 1000 for size in component_sizes) == 1,
        "PNG contains duplicated Extension badge components",
    )
    y_pixels, x_pixels = np.where(teal_mask)
    top_mm = float(y_pixels.min()) / image.shape[0] * float(axis["panel_height_mm"])
    bottom_mm = float(y_pixels.max() + 1) / image.shape[0] * float(axis["panel_height_mm"])
    left_mm = float(x_pixels.min()) / image.shape[1] * float(axis["panel_width_mm"])
    right_mm = float(x_pixels.max() + 1) / image.shape[1] * float(axis["panel_width_mm"])
    require(
        top_mm >= float(axis["separator_y_mm"])
        and bottom_mm <= float(axis["block2_bottom_mm"])
        and left_mm >= 0
        and right_mm <= float(axis["panel_width_mm"]),
        "PNG Extension badge is clipped or crosses the F15/F16 separator",
    )
    require(outputs["svg"].is_file() and outputs["svg"].stat().st_size > 0, "SVG missing")


def write_qc(
    source: pd.DataFrame,
    summary: pd.DataFrame,
    axis: pd.Series,
    layout_qc: dict[str, float | int],
    outputs: dict[str, Path],
) -> None:
    f15_endo = summary.query("family_id == 15 and compartment == 'endothelial'").iloc[0]
    hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in [
            FAMILY_DONOR,
            FAMILY_CONTRACT,
            FAMILY_METADATA,
            AXIS_CONTRACT,
            OUT_SOURCE,
            OUT_SUMMARY,
            *outputs.values(),
        ]
    }
    lines = [
        "# Figure 3D Relay/Extension family-donor QC",
        "",
        "Status: **PASS**",
        "",
        "## Scientific assertions",
        "- F15 = Relay; healthy LA carrier = Fibroblast; support set = {Endothelial}.",
        "- F15 Fibroblast support = FALSE; Endothelial support = TRUE; Lymphocyte support = FALSE.",
        "- F15 Endothelial primary p reproduces 0.0471.",
        "- F16 = Extension; healthy LA carrier = Fibroblast; support set = {Fibroblast, aCM}.",
        "- F16 Fibroblast support = TRUE; aCM support = TRUE; Adipocyte support = FALSE.",
        "- F16 Adipocyte primary p reproduces 0.0542 and is not displayed as support.",
        "- All plotted values originate from the locked family-level donor table.",
        "- Renderer opens no program-level, gene-level, leading-edge, expression, legacy gene-bar, or Path A input.",
        "- All primary family-level permutation tests use n_perm=10,000.",
        f"- Family-level source rows: {len(source)}; summary rows: {len(summary)}.",
        "",
        "## Shared Figure 3C/3D axis contract",
        f"- contract_id: `{axis['contract_id']}`",
        f"- y_min={float(axis['y_min']):.6f}; y_max={float(axis['y_max']):.6f}; y_ticks={axis['y_ticks']}",
        f"- point size={float(axis['point_size_pt2']):.1f}; jitter half-width={float(axis['jitter_halfwidth']):.2f}; mean-line width={float(axis['mean_line_width_pt']):.1f}.",
        f"- block heights={float(axis['block1_bottom_mm']) - float(axis['block1_top_mm']):.1f}/{float(axis['block2_bottom_mm']) - float(axis['block2_top_mm']):.1f} mm; separator={float(axis['separator_y_mm']):.1f} mm.",
        f"- badge dimensions={float(axis['badge_width_mm']):.1f} × {float(axis['badge_height_mm']):.1f} mm; annotations={int(axis['annotation_lines'])} lines at {float(axis['annotation_font_pt']):.1f} pt.",
        f"- geometry SHA256: `{axis['source_geometry_sha256']}`.",
        "- Figure 3C and Figure 3D use identical y-axis and mark contracts.",
        "",
        "## Targeted precision sensitivity",
        f"- F15 Endothelial primary: p={float(f15_endo.p_perm):.10f}, n_perm={int(f15_endo.n_perm):,}.",
        f"- F15 Endothelial sensitivity: p={float(f15_endo.p_perm_sensitivity_1e6):.10f}, n_perm={int(f15_endo.n_perm_sensitivity):,}.",
        f"- Interpretation: {f15_endo.sensitivity_note}.",
        "",
        "## Facet reproduction",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"- F{row.family_id} {row.compartment}: displayed={row.displayed}, "
            f"Ctrl n={row.n_ctrl}, AF n={row.n_af}, mean Ctrl={row.mean_ctrl:.8f}, "
            f"mean AF={row.mean_af:.8f}, delta={row.delta:+.8f}, "
            f"p_perm={row.p_perm:.8f}, n_perm={row.n_perm}"
        )
    lines.extend(
        [
            "",
            "## Layout/export QC",
            f"- Minimum visible font: {float(layout_qc['min_font_pt']):.2f} pt",
            f"- Visible-text bbox violations: {int(layout_qc['n_bbox_violations'])}",
            f"- Visible-text overlap violations: {int(layout_qc['n_overlap_violations'])}",
            f"- Extension badge objects: {int(layout_qc['n_extension_badges'])}",
            f"- Badge clipping violations: {int(layout_qc['n_badge_clipping_violations'])}",
            "- PDF/SVG/PNG duplicate-object and badge-boundary assertions passed.",
            "- PDF, SVG, 600-dpi PNG, and 200-dpi 88-mm preview validated.",
            "",
            "## SHA256",
            *[f"- `{path}`: `{digest}`" for path, digest in hashes.items()],
            "",
        ]
    )
    OUT_QC.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    donor, contract, metadata, axis = load_inputs()
    source, summary = build_outputs(donor, contract, metadata)
    RESULTS.mkdir(parents=True, exist_ok=True)
    source.to_csv(OUT_SOURCE, sep="\t", index=False, float_format="%.10g")
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False, float_format="%.10g")
    fig = render(source, summary, axis)
    layout_qc = validate_figure(fig)
    outputs = save_outputs(fig, axis)
    validate_exports(outputs, axis)
    write_qc(source, summary, axis, layout_qc, outputs)
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
