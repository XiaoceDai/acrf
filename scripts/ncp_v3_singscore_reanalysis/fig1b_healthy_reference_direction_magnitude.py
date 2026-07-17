#!/usr/bin/env python3
"""Standalone corrected Figure 1B renderer.

All inputs are resolved and checksum-verified through fig1b_input_manifest.tsv.
No scoring, reselection, or legacy Figure 1B code is invoked.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
MANIFEST = RESULTS / "fig1b_input_manifest.tsv"

OUT_BASE = FIGURES / "fig1b_healthy_reference_direction_magnitude"
OUT_SOURCE = RESULTS / "fig1b_healthy_reference_direction_magnitude_source.tsv"
OUT_SUMMARY = RESULTS / "fig1b_healthy_reference_direction_magnitude_summary.tsv"
OUT_QC = RESULTS / "fig1b_healthy_reference_direction_magnitude_qc.md"
OUT_PREVIEW = FIGURES / "fig1b_healthy_reference_direction_magnitude_preview.png"

SCORER_PACKAGE = "singscore"
SCORER_VERSION = "1.28.1"
REFERENCE_FORMULA = "mean(abs(direction_sign * delta))"
EXPECTED_MEDIAN_M = 0.0283547
NUMERICAL_TOLERANCE = 1e-12
MEDIAN_TOLERANCE = 1e-7
FAMILY_BLOCKS = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
DONOR_COLOR = "#858585"
RANGE_COLOR = "#A8A8A8"
TEXT_COLOR = "#222222"
MUTED_COLOR = "#5F5F5F"
GRID_COLOR = "#E4E4E4"
WHITE = "#FFFFFF"
MM = 1.0 / 25.4
FIGURE_WIDTH_MM = 102.3
FIGURE_HEIGHT_MM = 105


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def load_manifest() -> tuple[pd.DataFrame, dict[str, Path]]:
    require(MANIFEST.is_file(), f"Missing Figure 1B input manifest: {MANIFEST}")
    manifest = pd.read_csv(MANIFEST, sep="\t")
    required_columns = {
        "artifact_role",
        "artifact_path",
        "sha256",
        "generated_by",
        "dataset",
        "score_method",
        "scorer_package",
        "scorer_version",
        "leading_edge_version",
        "reference_formula",
        "is_authoritative",
    }
    require(
        required_columns.issubset(manifest.columns),
        f"Manifest missing columns: {sorted(required_columns - set(manifest.columns))}",
    )
    require(manifest["artifact_role"].is_unique, "Manifest artifact roles are not unique")
    require(manifest["is_authoritative"].astype(bool).all(), "Manifest has non-authoritative rows")
    require(
        set(manifest["scorer_package"]) == {SCORER_PACKAGE},
        "Manifest scorer package is not singscore",
    )
    require(
        set(manifest["scorer_version"].astype(str)) == {SCORER_VERSION},
        "Manifest scorer version is not singscore 1.28.1",
    )
    require(
        set(manifest["reference_formula"]) == {REFERENCE_FORMULA},
        "Manifest reference formula is not the approved formula",
    )
    require(
        manifest["leading_edge_version"]
        .str.contains("exact_match_selection_time_fgsea_30_of_30", regex=False)
        .all(),
        "Manifest does not certify a 30/30 selection-time fgsea match",
    )

    resolved: dict[str, Path] = {}
    namespace = RESULTS.resolve()
    for row in manifest.itertuples(index=False):
        path = (ROOT / row.artifact_path).resolve()
        require(path.is_relative_to(namespace), f"Input outside corrected namespace: {path}")
        require(path.is_file(), f"Manifest input does not exist: {path}")
        require(sha256(path) == row.sha256, f"Checksum mismatch: {path}")
        resolved[str(row.artifact_role)] = path
    expected_roles = {
        "healthy_donor_program_contrasts",
        "healthy_program_reference_summary",
        "frozen_program_metadata",
        "frozen_family_metadata",
        "canonical_display_order",
        "leading_edge_validation",
        "scorer_version_record",
    }
    require(set(resolved) == expected_roles, "Manifest roles differ from the approved input set")
    return manifest, resolved


def validate_sources(
    manifest: pd.DataFrame, paths: dict[str, Path]
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float | int | str], list[str]]:
    contrasts = pd.read_csv(paths["healthy_donor_program_contrasts"], sep="\t")
    reference = pd.read_csv(paths["healthy_program_reference_summary"], sep="\t")
    order = pd.read_csv(paths["canonical_display_order"], sep="\t")
    leading_edge = pd.read_csv(paths["leading_edge_validation"], sep="\t")

    require(
        int(leading_edge.loc[0, "n_programs"]) == 30
        and int(leading_edge.loc[0, "n_set_identical"]) == 30
        and bool(leading_edge.loc[0, "gate_passed"]),
        "Leading-edge sets are not an exact 30/30 selection-time fgsea match",
    )
    require(len(order) == 30 and order["program_id"].nunique() == 30, "Order table is not 30 programs")
    require(
        order["figure_display_order"].tolist() == list(range(1, 31)),
        "Figure display order is not contiguous 1..30",
    )
    observed_family_blocks = order["family_id"].drop_duplicates().astype(int).tolist()
    require(
        observed_family_blocks == FAMILY_BLOCKS,
        f"Canonical family blocks differ: {observed_family_blocks}",
    )
    require(
        order["canonical_display_order"].tolist() == order["figure_display_order"].tolist(),
        "Figure row order differs from frozen canonical display order",
    )

    require(contrasts["program_id"].nunique() == 30, "Contrast table is not 30 programs")
    donor_order = sorted(contrasts["donor_id"].astype(str).unique())
    require(len(donor_order) == 8, f"Expected 8 healthy donors, found {len(donor_order)}")
    require(len(contrasts) == 240, f"Expected 240 donor-program rows, found {len(contrasts)}")
    require(
        not contrasts.isna().any().any(),
        "Missing values detected; no missing value may be converted to zero",
    )
    require(
        int((contrasts["aCAI"] == 0).sum()) == 0,
        "Zero aligned aCAI detected; missing-to-zero conversion cannot be excluded",
    )
    require(
        bool((contrasts["aCAI"] > 0).all()),
        "At least one donor-program aligned aCAI is not positive",
    )
    require(
        np.allclose(
            contrasts["aCAI"],
            contrasts["direction_sign"] * contrasts["delta"],
            rtol=0,
            atol=NUMERICAL_TOLERANCE,
        ),
        "aCAI does not equal direction_sign * delta",
    )
    require(
        np.allclose(
            contrasts["aCAIrel"],
            contrasts["aCAI"] / contrasts["M"],
            rtol=0,
            atol=NUMERICAL_TOLERANCE,
        ),
        "aCAIrel does not equal aCAI / M",
    )
    per_program_n = contrasts.groupby("program_id")["donor_id"].nunique()
    require(bool((per_program_n == 8).all()), "A program does not have exactly 8 donors")

    recomputed = (
        contrasts.groupby("program_id")["aCAI"]
        .apply(lambda values: float(np.mean(np.abs(values))))
        .rename("M_recomputed")
    )
    reference_indexed = reference.set_index("program_id")
    require(len(reference_indexed) == 30, "Reference summary is not 30 programs")
    require(
        set(recomputed.index) == set(reference_indexed.index),
        "Contrast and reference program sets differ",
    )
    m_error = float(
        np.max(
            np.abs(
                recomputed.loc[reference_indexed.index].to_numpy()
                - reference_indexed["M"].to_numpy()
            )
        )
    )
    require(
        m_error <= NUMERICAL_TOLERANCE,
        f"M[p] differs from mean(abs(aCAI)); max error={m_error}",
    )
    concordance = contrasts.groupby("program_id")["aCAI"].apply(
        lambda values: int((values > 0).sum())
    )
    require(bool((concordance == 8).all()), "Not every program is 8/8 concordant")
    require(
        bool((reference_indexed["directional_concordance_k"] == 8).all()),
        "Reference summary does not report 8/8 for all programs",
    )
    direction_counts = order["canonical_direction"].value_counts().to_dict()
    require(
        direction_counts == {"LA-biased": 22, "RA-biased": 8},
        f"Expected 22 LA-biased and 8 RA-biased programs, found {direction_counts}",
    )
    median_m = float(reference["M"].median())
    require(
        np.isclose(median_m, EXPECTED_MEDIAN_M, rtol=0, atol=MEDIAN_TOLERANCE),
        f"Median M[p] {median_m} differs from approximately {EXPECTED_MEDIAN_M}",
    )

    source = contrasts.merge(
        order[
            [
                "program_id",
                "figure_display_order",
                "family_block_order",
                "family_id",
                "program_order",
                "publication_label",
                "canonical_direction",
            ]
        ],
        on=[
            "program_id",
            "family_id",
            "program_order",
            "canonical_direction",
        ],
        validate="many_to_one",
    )
    source = source.merge(
        recomputed.reset_index(), on="program_id", validate="many_to_one"
    )
    source["donor_min_aCAI"] = source.groupby("program_id")["aCAI"].transform("min")
    source["donor_max_aCAI"] = source.groupby("program_id")["aCAI"].transform("max")
    source = source.sort_values(["figure_display_order", "donor_id"]).reset_index(drop=True)
    require(
        source["program_id"].drop_duplicates().tolist() == order["program_id"].tolist(),
        "Figure source rows do not follow authoritative order",
    )

    summary: dict[str, float | int | str] = {
        "n_programs": 30,
        "n_healthy_donors": 8,
        "n_donor_program_values": 240,
        "n_la_biased": 22,
        "n_ra_biased": 8,
        "n_programs_8_of_8_concordant": 30,
        "median_M": median_m,
        "minimum_M": float(reference["M"].min()),
        "maximum_M": float(reference["M"].max()),
        "maximum_M_recomputation_error": m_error,
        "score_method": str(manifest["score_method"].iloc[0]),
        "scorer_package": SCORER_PACKAGE,
        "scorer_version": SCORER_VERSION,
        "reference_formula": REFERENCE_FORMULA,
        "leading_edge_match": "30/30 exact selection-time fgsea match",
    }
    return source, order, summary, donor_order


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 7.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.6,
            "savefig.facecolor": WHITE,
            "figure.facecolor": WHITE,
        }
    )


def publication_label(value: str) -> str:
    """Return the exact frozen publication-facing label."""
    return str(value)


def render(
    source: pd.DataFrame,
    order: pd.DataFrame,
    summary: dict[str, float | int | str],
    donor_order: list[str],
) -> plt.Figure:
    set_style()
    fig = plt.figure(
        figsize=(FIGURE_WIDTH_MM * MM, FIGURE_HEIGHT_MM * MM),
        constrained_layout=False,
    )
    grid = fig.add_gridspec(
        1,
        4,
        width_ratios=[0.70, 0.50, 2.60, 3.20],
        left=0.035,
        right=0.965,
        bottom=0.110,
        top=0.695,
        wspace=0.045,
    )
    ax_direction = fig.add_subplot(grid[0, 0])
    ax_family = fig.add_subplot(grid[0, 1], sharey=ax_direction)
    ax_label = fig.add_subplot(grid[0, 2], sharey=ax_direction)
    ax_value = fig.add_subplot(grid[0, 3], sharey=ax_direction)

    fig.text(
        0.035,
        0.975,
        "B",
        ha="left",
        va="top",
        fontsize=11.0,
        fontweight="bold",
        color=TEXT_COLOR,
    )
    fig.text(
        0.085,
        0.975,
        "Healthy-reference direction and magnitude\nacross 30 programs",
        ha="left",
        va="top",
        fontsize=8.5,
        fontweight="bold",
        color=TEXT_COLOR,
    )
    fig.text(
        0.035,
        0.905,
        r"$\mathrm{aCAI}_{d,p} = s_p\,\Delta_{d,p}$"
        r"     $M_p = \mathrm{mean}_{d\in H}\,|\mathrm{aCAI}_{d,p}|$",
        ha="left",
        va="top",
        fontsize=6.5,
        color=TEXT_COLOR,
    )
    fig.text(
        0.035,
        0.870,
        r"$\mathrm{aCAIrel}_{d,p} = \mathrm{aCAI}_{d,p}\,/\,M_p$",
        ha="left",
        va="top",
        fontsize=6.5,
        color=TEXT_COLOR,
    )
    fig.text(
        0.035,
        0.835,
        "All values are direction-aligned; color denotes the frozen LA/RA bias.",
        ha="left",
        va="top",
        fontsize=6.5,
        color=MUTED_COLOR,
    )
    ribbon = mpl.patches.FancyBboxPatch(
        (0.030, 0.725),
        0.935,
        0.082,
        boxstyle="round,pad=0.004,rounding_size=0.004",
        transform=fig.transFigure,
        facecolor="#F4F4F4",
        edgecolor="#D7D7D7",
        linewidth=0.55,
        zorder=0,
    )
    fig.add_artist(ribbon)
    fig.text(
        0.045,
        0.788,
        f"{summary['n_la_biased']} LA-biased | {summary['n_ra_biased']} RA-biased",
        ha="left",
        va="top",
        fontsize=6.7,
        fontweight="bold",
        color=TEXT_COLOR,
    )
    fig.text(
        0.045,
        0.753,
        f"{summary['n_programs_8_of_8_concordant']}/"
        f"{summary['n_programs']} programs concordant in all 8 reference donors",
        ha="left",
        va="top",
        fontsize=6.5,
        fontweight="bold",
        color=TEXT_COLOR,
    )
    fig.text(
        0.950,
        0.788,
        rf"median $M_p$ = {float(summary['median_M']):.4f}",
        ha="right",
        va="top",
        fontsize=6.7,
        fontweight="bold",
        color=TEXT_COLOR,
    )

    n_rows = len(order)
    for axis in (ax_direction, ax_family, ax_label, ax_value):
        axis.set_ylim(n_rows - 0.5, -0.5)
        axis.set_yticks([])
        axis.set_facecolor(WHITE)

    ax_direction.set_xlim(0, 1)
    ax_family.set_xlim(0, 1)
    ax_label.set_xlim(0, 1)
    for axis in (ax_direction, ax_family, ax_label):
        axis.axis("off")

    ax_direction.set_title("")
    ax_family.set_title("family", fontsize=6.5, fontweight="bold", pad=3)
    ax_label.set_title("program", fontsize=6.5, fontweight="bold", pad=3)
    ax_value.set_title(r"Aligned donor aCAI and $M_p$", fontsize=6.5, fontweight="bold", pad=3)

    donor_offsets = np.linspace(-0.13, 0.13, len(donor_order))
    ax_value.set_xlim(0, 0.40)
    ax_value.set_xticks(np.arange(0, 0.401, 0.10))
    ax_value.set_xticks(np.arange(0.05, 0.401, 0.10), minor=True)
    ax_value.axvline(0, color="#777777", lw=0.75, zorder=0)
    ax_value.grid(axis="x", color=GRID_COLOR, lw=0.45, zorder=0)
    ax_value.spines[["top", "right", "left"]].set_visible(False)
    ax_value.tick_params(axis="x", labelsize=6.6, length=2.5, width=0.5)
    ax_value.set_xlabel(
        "Aligned healthy chamber contrast\n(centered singscore units)",
        fontsize=6.5,
        labelpad=4,
    )

    grouped = source.groupby("program_id", sort=False)
    display_order = order.reset_index(drop=True)
    previous_family: int | None = None
    for row_index, row in display_order.iterrows():
        program_id = str(row["program_id"])
        family_id = int(row["family_id"])
        direction = str(row["canonical_direction"])
        color = LA_COLOR if direction == "LA-biased" else RA_COLOR
        values_frame = grouped.get_group(program_id).set_index("donor_id").loc[donor_order]
        values = values_frame["aCAI"].to_numpy(dtype=float)
        m_value = float(values_frame["M_recomputed"].iloc[0])

        ax_direction.add_patch(
            mpl.patches.Rectangle(
                (0, row_index - 0.48),
                1,
                0.96,
                facecolor=color,
                edgecolor=WHITE,
                linewidth=0.25,
            )
        )
        ax_label.text(
            0.01,
            row_index,
            publication_label(str(row["publication_label"])),
            ha="left",
            va="center",
            fontsize=6.6,
            color=TEXT_COLOR,
            fontfamily=(
                "DejaVu Sans"
                if "⁺" in str(row["publication_label"])
                else "sans-serif"
            ),
        )
        ax_value.plot(
            [float(values.min()), float(values.max())],
            [row_index, row_index],
            color=RANGE_COLOR,
            lw=1.0,
            solid_capstyle="round",
            zorder=1,
        )
        ax_value.scatter(
            values,
            row_index + donor_offsets,
            s=11.0,
            marker="o",
            facecolor=DONOR_COLOR,
            edgecolor=WHITE,
            linewidth=0.25,
            alpha=0.95,
            zorder=2,
        )
        ax_value.scatter(
            [m_value],
            [row_index],
            s=27,
            marker="D",
            facecolor=color,
            edgecolor=WHITE,
            linewidth=0.45,
            zorder=3,
        )

        if previous_family is not None and family_id != previous_family:
            y = row_index - 0.5
            for axis in (ax_direction, ax_family, ax_label, ax_value):
                axis.axhline(y, color="#D2D2D2", lw=0.45, zorder=0)
        previous_family = family_id

    for family_id, family_rows in display_order.groupby("family_id", sort=False):
        family_center = float(family_rows.index.to_numpy(dtype=float).mean())
        ax_family.text(
            0.5,
            family_center,
            f"F{int(family_id)}",
            ha="center",
            va="center",
            fontsize=6.7,
            color=TEXT_COLOR,
        )

    la_rows = int((order["canonical_direction"] == "LA-biased").sum())
    ra_rows = n_rows - la_rows
    ax_direction.text(
        0.5,
        (la_rows - 1) / 2,
        "LA",
        rotation=90,
        ha="center",
        va="center",
        color=WHITE,
        fontweight="bold",
        fontsize=6.5,
    )
    ax_direction.text(
        0.5,
        la_rows + (ra_rows - 1) / 2,
        "RA",
        rotation=90,
        ha="center",
        va="center",
        color=WHITE,
        fontweight="bold",
        fontsize=6.5,
    )
    donor_handle = mpl.lines.Line2D(
        [],
        [],
        linestyle="none",
        marker="o",
        markersize=3.4,
        markerfacecolor=DONOR_COLOR,
        markeredgecolor=WHITE,
    )
    magnitude_handle = mpl.lines.Line2D(
        [],
        [],
        linestyle="none",
        marker="D",
        markersize=4.6,
        markerfacecolor=LA_COLOR,
        markeredgecolor=WHITE,
    )
    direction_handle = (
        mpl.lines.Line2D(
            [], [], linestyle="none", marker="s", markersize=4.5,
            markerfacecolor=LA_COLOR, markeredgecolor="none"
        ),
        mpl.lines.Line2D(
            [], [], linestyle="none", marker="s", markersize=4.5,
            markerfacecolor=RA_COLOR, markeredgecolor="none"
        ),
    )
    ax_value.legend(
        handles=[donor_handle, magnitude_handle, direction_handle],
        labels=[
            "donor-level\naligned aCAI",
            "program-specific\n" + r"$M_p$",
            "color = frozen\nchamber direction",
        ],
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.25)},
        loc="center right",
        bbox_to_anchor=(1.0, 0.60),
        frameon=False,
        fontsize=6.5,
        ncol=1,
        borderaxespad=0.3,
        handletextpad=0.5,
        labelspacing=0.35,
    )
    return fig


def write_qc(
    manifest: pd.DataFrame,
    source: pd.DataFrame,
    summary: dict[str, float | int | str],
    figure_hashes: dict[str, str],
) -> None:
    manifest_rows = "\n".join(
        f"- `{row.artifact_role}`: `{row.artifact_path}` — `{row.sha256}`"
        for row in manifest.itertuples(index=False)
    )
    figure_rows = "\n".join(
        f"- `{path}` — `{digest}`" for path, digest in figure_hashes.items()
    )
    forbidden_import = "fig1b_baseline_" + "effect_concordance"
    source_code = Path(__file__).read_text(encoding="utf-8")
    require(forbidden_import not in source_code, "Legacy renderer import detected")
    qc = f"""# Corrected Figure 1B QC

Status: **PASS**

## Scientific assertions
- scorer: `{SCORER_PACKAGE} {SCORER_VERSION}`
- leading edges: 30/30 exact match to selection-time fgsea
- programs: {summary['n_programs']}
- healthy donors: {summary['n_healthy_donors']}
- aligned donor-program values: {summary['n_donor_program_values']}, all positive
- direction split: {summary['n_la_biased']} LA-biased / {summary['n_ra_biased']} RA-biased
- concordance: {summary['n_programs_8_of_8_concordant']}/30 programs at 8/8
- median M[p]: {float(summary['median_M']):.10f}
- M[p] range: {float(summary['minimum_M']):.10f}–{float(summary['maximum_M']):.10f}
- maximum M[p] recomputation error: {float(summary['maximum_M_recomputation_error']):.3g}
- reference formula: `{summary['reference_formula']}`
- missing values converted to zero: no
- legacy renderer imported: no
- legacy numerical table read: no
- figure width: {FIGURE_WIDTH_MM} mm
- figure height: {FIGURE_HEIGHT_MM} mm
- source rows: {len(source)}

## Authoritative input hashes
{manifest_rows}

## Figure hashes
{figure_rows}

## Preview
`{OUT_PREVIEW.relative_to(ROOT)}` is a reduced-resolution preview rendered at the
same {FIGURE_WIDTH_MM}-mm final panel width.
"""
    OUT_QC.write_text(qc, encoding="utf-8")


def main() -> None:
    manifest, paths = load_manifest()
    source, order, summary, donor_order = validate_sources(manifest, paths)

    # All hard assertions have passed before any figure or source output is created.
    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    source_columns = [
        "figure_display_order",
        "family_block_order",
        "family_id",
        "program_order",
        "program_id",
        "publication_label",
        "canonical_direction",
        "direction_sign",
        "dataset",
        "donor_id",
        "score_LA",
        "score_RA",
        "delta",
        "aCAI",
        "M",
        "M_recomputed",
        "aCAIrel",
        "donor_min_aCAI",
        "donor_max_aCAI",
        "score_method",
        "gene_set_version",
    ]
    source[source_columns].to_csv(OUT_SOURCE, sep="\t", index=False)
    pd.DataFrame([summary]).to_csv(OUT_SUMMARY, sep="\t", index=False)

    fig = render(source, order, summary, donor_order)
    pdf = OUT_BASE.with_suffix(".pdf")
    svg = OUT_BASE.with_suffix(".svg")
    png = OUT_BASE.with_suffix(".png")
    fig.savefig(pdf, bbox_inches=None, pad_inches=0)
    fig.savefig(svg, bbox_inches=None, pad_inches=0)
    fig.savefig(png, dpi=600, bbox_inches=None, pad_inches=0)
    fig.savefig(OUT_PREVIEW, dpi=120, bbox_inches=None, pad_inches=0)
    plt.close(fig)

    figure_hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (pdf, svg, png, OUT_PREVIEW)
    }
    write_qc(manifest, source, summary, figure_hashes)
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")
    print(
        f"Summary: 22 LA-biased | 8 RA-biased; 30/30 at 8/8; "
        f"median M[p]={float(summary['median_M']):.10f}; "
        f"range={float(summary['minimum_M']):.10f}–"
        f"{float(summary['maximum_M']):.10f}"
    )


if __name__ == "__main__":
    main()
