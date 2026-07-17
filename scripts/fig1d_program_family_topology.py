#!/usr/bin/env python3
"""Standalone Figure 1D: frozen 30-program to 16-family topology."""
from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

FAMILY_MAP = RESULTS / "innovation_ncp_v3_axis_families.tsv"
FAMILY_LABELS = RESULTS / "ncp_v3_family_labels.tsv"
CANONICAL_ORDER = RESULTS / "ncp_v3_canonical_display_order.tsv"
PROGRAM_METADATA = (
    RESULTS
    / "ncp_v3_singscore_reanalysis"
    / "fig1b_canonical_display_order.tsv"
)
ALLOWED_INPUTS = {
    FAMILY_MAP,
    FAMILY_LABELS,
    CANONICAL_ORDER,
    PROGRAM_METADATA,
}

OUT_BASE = FIGURES / "fig1d_program_family_topology"
OUT_SOURCE = RESULTS / "fig1d_program_family_topology_source.tsv"
OUT_SUMMARY = RESULTS / "fig1d_program_family_topology_summary.tsv"
OUT_QC = RESULTS / "fig1d_program_family_topology_qc.md"
OUT_PREVIEW = FIGURES / "fig1d_program_family_topology_preview.png"

FAMILY_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]
LA_FAMILIES = FAMILY_ORDER[:11]
RA_FAMILIES = FAMILY_ORDER[11:]

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
FG_COLOR = "#1A1A1A"
MUTED_COLOR = "#555555"
LINE_COLOR = "#D1D1D1"
WHITE = "#FFFFFF"
MM = 1.0 / 25.4
WIDTH_MM = 102.3
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


def wrap_programs(labels: list[str], max_chars: int = 42) -> str:
    joined = " · ".join(labels)
    if len(joined) <= max_chars:
        return joined
    if len(labels) >= 5:
        lines = [
            " · ".join(labels[:2]),
            " · ".join(labels[2:4]),
            " · ".join(labels[4:]),
        ]
        require(
            all(len(line) <= max_chars for line in lines),
            f"Constituent program labels cannot be compactly wrapped: {joined}",
        )
        return "\n".join(lines)
    candidates: list[tuple[int, str, str]] = []
    for split_at in range(1, len(labels)):
        left = " · ".join(labels[:split_at])
        right = " · ".join(labels[split_at:])
        penalty = abs(len(left) - len(right))
        if len(left) > max_chars or len(right) > max_chars:
            penalty += 1000
        candidates.append((penalty, left, right))
    require(bool(candidates), "Unable to wrap constituent program labels")
    _, left, right = min(candidates, key=lambda item: item[0])
    require(
        len(left) <= max_chars and len(right) <= max_chars,
        f"Constituent program labels require more than two lines: {joined}",
    )
    return f"{left}\n{right}"


def wrap_family_label(label: str, max_chars: int = 20) -> str:
    if len(label) <= max_chars:
        return label
    words = label.split()
    candidates = []
    for split_at in range(1, len(words)):
        left = " ".join(words[:split_at])
        right = " ".join(words[split_at:])
        penalty = abs(len(left) - len(right))
        if len(left) > max_chars or len(right) > max_chars:
            penalty += 1000
        candidates.append((penalty, left, right))
    require(bool(candidates), f"Unable to wrap family label: {label}")
    _, left, right = min(candidates, key=lambda item: item[0])
    require(
        len(left) <= max_chars and len(right) <= max_chars,
        f"Family label requires more than two lines: {label}",
    )
    return f"{left}\n{right}"


def display_family_label(value: str) -> str:
    return str(value)


def build_source() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int | str]]:
    for path in ALLOWED_INPUTS:
        require(path.is_file(), f"Missing frozen topology input: {path}")

    mapping = pd.read_csv(FAMILY_MAP, sep="\t").rename(
        columns={"pathway": "program_id", "direction": "mapping_direction"}
    )
    labels = pd.read_csv(FAMILY_LABELS, sep="\t")
    canonical = pd.read_csv(CANONICAL_ORDER, sep="\t")
    program_meta = pd.read_csv(PROGRAM_METADATA, sep="\t")

    require(len(mapping) == 30, "Family mapping does not contain exactly 30 rows")
    require(mapping["program_id"].nunique() == 30, "Every program must map exactly once")
    require(len(canonical) == 30 and canonical["program_id"].nunique() == 30, "Canonical order is not 30 unique programs")
    require(len(program_meta) == 30 and program_meta["program_id"].nunique() == 30, "Figure 1B metadata is not 30 unique programs")
    require(labels["family_id"].nunique() == 16 and len(labels) == 16, "Family labels are not 16 unique families")

    program_sets = [
        set(mapping["program_id"]),
        set(canonical["program_id"]),
        set(program_meta["program_id"]),
    ]
    require(program_sets[0] == program_sets[1] == program_sets[2], "Frozen program sets differ across topology inputs")

    source = canonical.merge(
        mapping[
            [
                "program_id",
                "program_order",
                "family_id",
                "mapping_direction",
                "family_source",
            ]
        ],
        on="program_id",
        suffixes=("_canonical", "_mapping"),
        validate="one_to_one",
    )
    source = source.merge(
        program_meta[
            [
                "program_id",
                "publication_label",
                "canonical_direction",
                "program_order",
                "family_id",
                "canonical_display_order",
            ]
        ],
        on="program_id",
        suffixes=("", "_fig1b"),
        validate="one_to_one",
    )
    source = source.merge(
        labels[["family_id", "family_label", "short_label", "axis_group"]],
        on="family_id",
        validate="many_to_one",
    )

    require(
        bool(
            (
                source["program_order_canonical"]
                == source["program_order_mapping"]
            ).all()
        )
        and bool((source["program_order_canonical"] == source["program_order"]).all()),
        "Frozen program_order differs across topology inputs",
    )
    require(
        bool((source["family_id_canonical"] == source["family_id_mapping"]).all())
        and bool((source["family_id_canonical"] == source["family_id"]).all()),
        "Frozen family assignment differs across topology inputs",
    )
    require(
        bool(
            (
                source["canonical_display_order"]
                == source["canonical_display_order_fig1b"]
            ).all()
        ),
        "Canonical display order differs from Figure 1B metadata",
    )
    require(
        bool((source["canonical_direction"] == source["mapping_direction"]).all())
        and bool(
            (
                source["canonical_direction"]
                == source["canonical_direction_fig1b"]
            ).all()
        ),
        "Frozen chamber direction differs across topology inputs",
    )

    source = source.sort_values("canonical_display_order").reset_index(drop=True)
    require(
        source["canonical_display_order"].tolist() == list(range(1, 31)),
        "Canonical program order is not contiguous 1..30",
    )
    observed_family_order = source["family_id"].drop_duplicates().astype(int).tolist()
    require(
        observed_family_order == FAMILY_ORDER,
        f"Family order is not the locked direction-blocked order: {observed_family_order}",
    )
    require(
        source["family_display_order"].drop_duplicates().astype(int).tolist()
        == list(range(1, 17)),
        "Family display order is not contiguous 1..16",
    )
    for _, group in source.groupby("family_id", sort=False):
        require(
            group["program_order"].tolist()
            == sorted(group["program_order"].tolist()),
            f"Within-family order changed for F{int(group.iloc[0]['family_id'])}",
        )

    require(source["family_id"].nunique() == 16, "Expected exactly 16 families")
    family_counts = source.groupby("family_id").size()
    require(int(family_counts.sum()) == 30, "Family program counts do not sum to 30")
    require(bool((family_counts >= 1).all()), "A family contains no programs")
    family_direction_count = source.groupby("family_id")["canonical_direction"].nunique()
    require(bool((family_direction_count == 1).all()), "A family mixes LA and RA directions")
    n_la_programs = int((source["canonical_direction"] == "LA-biased").sum())
    n_ra_programs = int((source["canonical_direction"] == "RA-biased").sum())
    require((n_la_programs, n_ra_programs) == (22, 8), "Expected 22 LA and 8 RA programs")
    family_directions = source.groupby("family_id", sort=False)["canonical_direction"].first()
    n_la_families = int((family_directions == "LA-biased").sum())
    n_ra_families = int((family_directions == "RA-biased").sum())
    require((n_la_families, n_ra_families) == (11, 5), "Expected 11 LA and 5 RA families")
    require(
        source[
            ["publication_label", "family_label", "short_label"]
        ].notna().all().all(),
        "A publication-facing label is missing",
    )
    require(
        bool(
            source[
                ["publication_label", "family_label", "short_label"]
            ].apply(lambda column: column.astype(str).str.strip().ne("")).all().all()
        ),
        "A publication-facing label is blank",
    )
    expected_axis = source["canonical_direction"].map(
        {"LA-biased": "LA", "RA-biased": "RA"}
    )
    require(bool((expected_axis == source["axis_group"]).all()), "Family axis_group disagrees with frozen direction")

    family_rows = []
    for family_index, (family_id, group) in enumerate(
        source.groupby("family_id", sort=False), start=1
    ):
        group = group.sort_values("program_order")
        family_rows.append(
            {
                "family_display_order": family_index,
                "family_id": int(family_id),
                "canonical_direction": str(group.iloc[0]["canonical_direction"]),
                "family_label": str(group.iloc[0]["short_label"]),
                "family_label_wrapped": wrap_family_label(
                    str(group.iloc[0]["short_label"])
                ),
                "constituent_program_labels": " · ".join(
                    group["publication_label"].astype(str)
                ),
                "constituent_program_labels_wrapped": wrap_programs(
                    group["publication_label"].astype(str).tolist()
                ),
                "n_programs": int(len(group)),
                "program_ids": ";".join(group["program_id"].astype(str)),
                "program_orders": ";".join(
                    group["program_order"].astype(int).astype(str)
                ),
            }
        )
    family_table = pd.DataFrame(family_rows)
    require(len(family_table) == 16, "Family-level display table is not 16 rows")

    summary: dict[str, int | str] = {
        "frozen_programs": 30,
        "fixed_families": 16,
        "LA_biased_families": n_la_families,
        "RA_biased_families": n_ra_families,
        "LA_biased_programs": n_la_programs,
        "RA_biased_programs": n_ra_programs,
        "aggregation_rule": "program-first equal-weight mean",
        "disallowed_HCA_or_disease_inputs_read": 0,
    }
    output_columns = [
        "canonical_display_order",
        "family_display_order",
        "family_id",
        "canonical_direction",
        "program_order",
        "program_id",
        "publication_label",
        "short_label",
        "family_label",
        "axis_group",
        "family_source",
    ]
    return source[output_columns], family_table, summary


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
        }
    )


def render(
    family_table: pd.DataFrame, summary: dict[str, int | str]
) -> plt.Figure:
    set_style()
    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM))
    fig.text(
        0.040,
        0.978,
        "D",
        ha="left",
        va="top",
        fontsize=11.0,
        fontweight="bold",
        color=FG_COLOR,
    )
    fig.text(
        0.090,
        0.978,
        "Frozen NCP program-to-family topology",
        ha="left",
        va="top",
        fontsize=8.5,
        fontweight="bold",
        color=FG_COLOR,
    )
    fig.text(
        0.040,
        0.937,
        f"{summary['frozen_programs']} frozen programs \u2192 "
        f"{summary['fixed_families']} families",
        ha="left",
        va="top",
        fontsize=7.7,
        fontweight="bold",
        color=FG_COLOR,
    )
    fig.text(
        0.040,
        0.908,
        f"{summary['LA_biased_families']} LA-biased families | "
        f"{summary['RA_biased_families']} RA-biased families",
        ha="left",
        va="top",
        fontsize=7.3,
        fontweight="bold",
        color=FG_COLOR,
    )
    fig.text(
        0.040,
        0.879,
        "Fixed membership, direction and display order used in Figs 2\u20136",
        ha="left",
        va="top",
        fontsize=6.7,
        color=MUTED_COLOR,
    )

    ax = fig.add_axes([0.040, 0.095, 0.92, 0.745])
    ax.set_xlim(0, 1)
    ax.set_ylim(16.85, -0.80)
    ax.axis("off")

    x_dir0, x_dir1 = 0.005, 0.050
    x_family_id = 0.085
    x_family_label = 0.155
    x_programs = 0.470
    x_count = 0.995
    header_y = -0.34
    ax.text(x_family_id, header_y, "ID", ha="center", va="center", fontsize=6.6, fontweight="bold")
    ax.text(x_family_label, header_y, "family", ha="left", va="center", fontsize=6.6, fontweight="bold")
    ax.text(x_programs, header_y, "constituent programs", ha="left", va="center", fontsize=6.6, fontweight="bold")
    ax.text(x_count, header_y, "n", ha="right", va="center", fontsize=6.6, fontweight="bold")
    ax.plot([0.0, 1.0], [0.08, 0.08], color=FG_COLOR, lw=0.75)

    la_top, la_bottom = 0.19, 11.05
    ra_top, ra_bottom = 11.19, 16.05
    ax.add_patch(
        mpl.patches.Rectangle(
            (x_dir0, la_top),
            x_dir1 - x_dir0,
            la_bottom - la_top,
            facecolor=LA_COLOR,
            edgecolor="none",
        )
    )
    ax.add_patch(
        mpl.patches.Rectangle(
            (x_dir0, ra_top),
            x_dir1 - x_dir0,
            ra_bottom - ra_top,
            facecolor=RA_COLOR,
            edgecolor="none",
        )
    )

    for row_index, row in family_table.iterrows():
        y = row_index + 0.62
        ax.text(
            x_family_id,
            y,
            f"F{int(row['family_id'])}",
            ha="center",
            va="center",
            fontsize=6.7,
            color=MUTED_COLOR,
        )
        ax.text(
            x_family_label,
            y,
            display_family_label(str(row["family_label_wrapped"])),
            ha="left",
            va="center",
            fontsize=6.5,
            fontweight="bold",
            color=FG_COLOR,
            fontfamily=(
                "DejaVu Sans"
                if "⁺" in str(row["family_label_wrapped"])
                else "sans-serif"
            ),
        )
        ax.text(
            x_programs,
            y,
            display_family_label(
                str(row["constituent_program_labels_wrapped"])
            ),
            ha="left",
            va="center",
            fontsize=6.5,
            color=MUTED_COLOR,
            linespacing=0.90,
            fontfamily=(
                "DejaVu Sans"
                if "⁺" in str(row["constituent_program_labels_wrapped"])
                else "sans-serif"
            ),
        )
        ax.text(
            x_count,
            y,
            str(int(row["n_programs"])),
            ha="right",
            va="center",
            fontsize=6.7,
            color=FG_COLOR,
            fontweight="bold",
        )
        separator_lw = 1.15 if row_index == 10 else 0.42
        separator_color = FG_COLOR if row_index == 10 else LINE_COLOR
        ax.plot(
            [0.060 if row_index != 10 else 0.0, 1.0],
            [y + 0.49, y + 0.49],
            color=separator_color,
            lw=separator_lw,
        )

    ax.text(
        (x_dir0 + x_dir1) / 2,
        (0.62 + 10.62) / 2,
        "LA-biased",
        ha="center",
        va="center",
        rotation=90,
        fontsize=6.5,
        fontweight="bold",
        color=WHITE,
    )
    ax.text(
        (x_dir0 + x_dir1) / 2,
        (11.62 + 15.62) / 2,
        "RA-biased",
        ha="center",
        va="center",
        rotation=90,
        fontsize=6.5,
        fontweight="bold",
        color=WHITE,
    )

    fig.text(
        0.040,
        0.055,
        "Primary transcriptomic family summaries use equal-weight program-first\n"
        "aggregation.",
        ha="left",
        va="top",
        fontsize=6.7,
        color=FG_COLOR,
        linespacing=1.15,
    )
    return fig


def write_qc(
    source: pd.DataFrame,
    family_table: pd.DataFrame,
    summary: dict[str, int | str],
    figure_hashes: dict[str, str],
) -> None:
    input_lines = "\n".join(
        f"- `{path.relative_to(ROOT)}` — `{sha256(path)}`"
        for path in sorted(ALLOWED_INPUTS)
    )
    output_lines = "\n".join(
        f"- `{path}` — `{digest}`" for path, digest in figure_hashes.items()
    )
    qc = f"""# Figure 1D frozen topology QC

Status: **PASS**

## Topology assertions
- unique frozen programs: {summary['frozen_programs']}
- unique fixed families: {summary['fixed_families']}
- mapping rows: {len(source)}
- family display rows: {len(family_table)}
- family program-count sum: {int(family_table['n_programs'].sum())}
- direction composition: {summary['LA_biased_programs']} LA / {summary['RA_biased_programs']} RA programs
- family direction composition: {summary['LA_biased_families']} LA / {summary['RA_biased_families']} RA families
- family order: {','.join('F' + str(value) for value in family_table['family_id'])}
- all family and program publication labels present: yes
- within-family frozen program_order preserved: yes
- program-first equal-weight aggregation stated: yes
- disallowed HCA or disease inputs read: {summary['disallowed_HCA_or_disease_inputs_read']}
- panel dimensions: {WIDTH_MM} × {HEIGHT_MM} mm
- minimum publication-facing text size: 6.5 pt

## Authoritative input hashes
{input_lines}

## Figure hashes
{output_lines}
"""
    OUT_QC.write_text(qc, encoding="utf-8")


def main() -> None:
    source, family_table, summary = build_source()

    # All frozen-topology assertions pass before any output is created.
    source.to_csv(OUT_SOURCE, sep="\t", index=False)
    pd.DataFrame(
        [{"metric": key, "value": value} for key, value in summary.items()]
    ).to_csv(OUT_SUMMARY, sep="\t", index=False)

    fig = render(family_table, summary)
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
    write_qc(source, family_table, summary, figure_hashes)
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")
    print(
        "Summary: 30 frozen programs -> 16 families; "
        "11 LA / 5 RA families; 22 LA / 8 RA programs"
    )


if __name__ == "__main__":
    main()
