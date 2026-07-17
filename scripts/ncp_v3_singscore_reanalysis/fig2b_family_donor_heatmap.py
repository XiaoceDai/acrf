#!/usr/bin/env python3
"""Corrected Figure 2B: frozen-family aCAIrel across 51 paired AF donors."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures"
MANIFEST = RESULTS / "fig2a_input_manifest.tsv"
OUT_SOURCE = RESULTS / "fig2b_family_donor_source.tsv"
OUT_QC = RESULTS / "fig2b_family_donor_qc.md"
OUT_BASE = FIGURES / "fig2b_donor_acai_heatmap"
OUT_PREVIEW = FIGURES / "fig2b_donor_acai_heatmap_preview.png"

COHORT_ORDER = ["GSE128188", "GSE79768", "GSE115574"]
RHYTHM_ORDER = ["SR", "AF"]
EXPECTED_COUNTS = {
    "GSE128188": {"SR": 5, "AF": 5},
    "GSE79768": {"SR": 6, "AF": 7},
    "GSE115574": {"SR": 15, "AF": 13},
}
REQUIRED_ROLES = {
    "frozen_program_metadata",
    "frozen_family_metadata",
    "disease_program_contrasts",
    "family_acairel",
    "GSE128188_sample_metadata",
    "GSE79768_sample_metadata",
    "GSE115574_sample_metadata",
}
SOURCE_COLUMNS = [
    "family_id",
    "family_label",
    "family_direction",
    "canonical_row_order",
    "dataset",
    "donor_id",
    "rhythm",
    "canonical_column_order",
    "n_constituent_programs",
    "family_aCAIrel",
]

MM = 1 / 25.4
WIDTH_MM = 178
HEIGHT_MM = 105
LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
FG = "#1A1A1A"
MUTED = "#595959"
LIGHT_LINE = "#D8D8D8"
COHORT_COLORS = {
    "GSE128188": "#D8D8D8",
    "GSE79768": "#BEBEBE",
    "GSE115574": "#A4A4A4",
}
RHYTHM_COLORS = {"SR": "#F5F5F5", "AF": "#555555"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def natural_key(value: object) -> tuple[object, ...]:
    return tuple(
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", str(value))
    )


def load_manifest() -> tuple[pd.DataFrame, dict[str, Path]]:
    require(MANIFEST.is_file(), f"Missing corrected manifest: {MANIFEST}")
    manifest = pd.read_csv(MANIFEST, sep="\t")
    required_columns = {
        "artifact_role",
        "artifact_path",
        "sha256",
        "is_authoritative",
        "allows_cohort_specific_reselection",
        "allows_reference_recalibration",
    }
    require(required_columns.issubset(manifest.columns), "Corrected manifest schema is incomplete")
    selected = manifest.loc[manifest["artifact_role"].isin(REQUIRED_ROLES)].copy()
    require(set(selected["artifact_role"]) == REQUIRED_ROLES, "Required Figure 2B manifest roles are missing")
    require(selected["artifact_role"].is_unique, "Figure 2B manifest roles are not unique")
    require(
        selected["is_authoritative"].astype(str).str.lower().eq("true").all(),
        "A Figure 2B input is not authoritative",
    )
    require(
        selected["allows_cohort_specific_reselection"].astype(str).str.lower().eq("false").all(),
        "Manifest permits cohort-specific re-selection",
    )
    require(
        selected["allows_reference_recalibration"].astype(str).str.lower().eq("false").all(),
        "Manifest permits reference recalibration",
    )

    paths: dict[str, Path] = {}
    for row in selected.itertuples(index=False):
        path = (ROOT / row.artifact_path).resolve()
        require(path.is_file(), f"Missing manifest input: {path}")
        require(sha256(path) == row.sha256, f"Checksum mismatch: {path}")
        paths[str(row.artifact_role)] = path
    return selected, paths


def validate_paired_metadata(paths: dict[str, Path]) -> pd.DataFrame:
    donor_rows: list[pd.DataFrame] = []
    for cohort in COHORT_ORDER:
        metadata = pd.read_csv(paths[f"{cohort}_sample_metadata"], sep="\t")
        included = metadata.loc[
            metadata["analysis_inclusion"].astype(str).str.lower().eq("true")
        ].copy()
        included["donor_id"] = included["donor_id"].astype(str)
        chambers = included.groupby("donor_id")["chamber"].agg(lambda values: set(values))
        require(chambers.map(lambda value: value == {"LA", "RA"}).all(), f"Incomplete pair in {cohort}")
        rhythms = included.groupby("donor_id")["rhythm"].agg(lambda values: set(values))
        require(rhythms.map(len).eq(1).all(), f"Mixed rhythm within donor in {cohort}")
        donors = included.drop_duplicates("donor_id")[["dataset", "donor_id", "rhythm"]]
        counts = donors["rhythm"].value_counts().to_dict()
        require(counts == EXPECTED_COUNTS[cohort], f"{cohort} counts changed: {counts}")
        donor_rows.append(donors)
    donors = pd.concat(donor_rows, ignore_index=True)
    require(len(donors) == 51, "Expected 51 paired disease donors")
    require((donors["rhythm"] == "SR").sum() == 26, "Expected 26 SR donors")
    require((donors["rhythm"] == "AF").sum() == 25, "Expected 25 AF donors")
    return donors


def build_source(
    paths: dict[str, Path],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float | int | str]]:
    programs = pd.read_csv(paths["frozen_program_metadata"], sep="\t")
    families = pd.read_csv(paths["frozen_family_metadata"], sep="\t")
    program_values = pd.read_csv(paths["disease_program_contrasts"], sep="\t")
    family_values = pd.read_csv(paths["family_acairel"], sep="\t")
    donors = validate_paired_metadata(paths)

    for frame in (program_values, family_values, donors):
        frame["donor_id"] = frame["donor_id"].astype(str)

    require(len(programs) == 30 and programs["program_id"].nunique() == 30, "Frozen program count changed")
    require(len(families) == 16 and families["family_id"].nunique() == 16, "Frozen family count changed")
    require(programs["family_id"].nunique() == 16, "Frozen program-family topology changed")
    expected_sizes = programs.groupby("family_id")["program_id"].nunique().rename("expected_programs")

    program_values = program_values.loc[program_values["dataset"].isin(COHORT_ORDER)].copy()
    family_values = family_values.loc[family_values["dataset"].isin(COHORT_ORDER)].copy()
    require(len(program_values) == 51 * 30, "Corrected program table is not 51 × 30")
    require(len(family_values) == 51 * 16, "Authoritative family slice is not 51 × 16")
    require(program_values["aCAIrel"].notna().all(), "Corrected program aCAIrel contains missing values")
    require(family_values["family_aCAIrel"].notna().all(), "Family aCAIrel contains missing values")

    recomputed = (
        program_values.groupby(
            ["dataset", "donor_id", "rhythm", "family_id"], as_index=False
        )
        .agg(
            family_aCAIrel_recomputed=("aCAIrel", "mean"),
            n_constituent_programs=("program_id", "nunique"),
        )
        .merge(expected_sizes, on="family_id", how="left", validate="many_to_one")
    )
    require(
        (recomputed["n_constituent_programs"] == recomputed["expected_programs"]).all(),
        "A donor-family cell omits a frozen constituent program",
    )

    authoritative = family_values[
        [
            "dataset",
            "donor_id",
            "rhythm",
            "family_id",
            "family_aCAIrel",
            "n_programs_observed",
            "family_label",
            "canonical_direction",
            "canonical_family_display_order",
        ]
    ].copy()
    checked = authoritative.merge(
        recomputed,
        on=["dataset", "donor_id", "rhythm", "family_id"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    require((checked["_merge"] == "both").all(), "Authoritative and recomputed family cells do not align")
    checked["absolute_verification_error"] = (
        checked["family_aCAIrel"] - checked["family_aCAIrel_recomputed"]
    ).abs()
    max_error = float(checked["absolute_verification_error"].max())
    require(max_error <= 1e-12, f"Family values fail equal-weight verification: {max_error}")
    require(
        (checked["n_programs_observed"] == checked["n_constituent_programs"]).all(),
        "Authoritative constituent-program counts disagree",
    )

    donor_keys = donors[["dataset", "donor_id", "rhythm"]].drop_duplicates()
    checked_keys = checked[["dataset", "donor_id", "rhythm"]].drop_duplicates()
    require(
        set(map(tuple, donor_keys.to_numpy())) == set(map(tuple, checked_keys.to_numpy())),
        "Family table donors differ from corrected paired metadata",
    )

    ordered_donors: list[tuple[str, str, str]] = []
    for cohort in COHORT_ORDER:
        for rhythm in RHYTHM_ORDER:
            subset = donor_keys.loc[
                donor_keys["dataset"].eq(cohort) & donor_keys["rhythm"].eq(rhythm)
            ]
            for donor_id in sorted(subset["donor_id"].tolist(), key=natural_key):
                ordered_donors.append((cohort, str(donor_id), rhythm))
    require(len(ordered_donors) == 51, "Canonical donor order is not length 51")
    column_order = {
        key: index + 1 for index, key in enumerate(ordered_donors)
    }

    source = checked.copy()
    source["family_direction"] = source["canonical_direction"]
    source["canonical_row_order"] = source["canonical_family_display_order"].astype(int)
    source["canonical_column_order"] = [
        column_order[(row.dataset, str(row.donor_id), row.rhythm)]
        for row in source.itertuples(index=False)
    ]
    source = source.rename(columns={"n_programs_observed": "authoritative_n_programs"})
    source = source[SOURCE_COLUMNS].sort_values(
        ["canonical_row_order", "canonical_column_order"]
    )

    require(source.shape == (816, 10), "Figure 2B source table has unexpected dimensions")
    require(source["family_id"].nunique() == 16, "Figure 2B source does not have 16 families")
    require(source["canonical_column_order"].nunique() == 51, "Figure 2B source does not have 51 donors")
    require(
        source.groupby("canonical_column_order")["family_id"].nunique().eq(16).all(),
        "A donor is missing a frozen family",
    )
    require(
        source.groupby("canonical_row_order")["canonical_column_order"].nunique().eq(51).all(),
        "A family is missing a donor",
    )

    values = source["family_aCAIrel"].to_numpy(dtype=float)
    quantiles = np.quantile(values, [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99])
    stats: dict[str, float | int | str] = {
        "n_families": 16,
        "n_donors": 51,
        "n_SR": 26,
        "n_AF": 25,
        "min": float(values.min()),
        "max": float(values.max()),
        "q01": float(quantiles[0]),
        "q05": float(quantiles[1]),
        "q25": float(quantiles[2]),
        "median": float(quantiles[3]),
        "q75": float(quantiles[4]),
        "q95": float(quantiles[5]),
        "q99": float(quantiles[6]),
        "max_verification_error": max_error,
        "clipping_policy": "none",
    }
    return source, families, stats


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
    source: pd.DataFrame,
    families: pd.DataFrame,
    stats: dict[str, float | int | str],
) -> plt.Figure:
    set_style()
    row_meta = (
        source[
            [
                "family_id",
                "family_label",
                "family_direction",
                "canonical_row_order",
            ]
        ]
        .drop_duplicates()
        .sort_values("canonical_row_order")
    )
    short_labels = families[
        ["family_id", "short_label", "canonical_family_display_order"]
    ].rename(columns={"canonical_family_display_order": "canonical_row_order"})
    row_meta = row_meta.merge(
        short_labels,
        on=["family_id", "canonical_row_order"],
        how="left",
        validate="one_to_one",
    )
    donors = (
        source[
            ["dataset", "donor_id", "rhythm", "canonical_column_order"]
        ]
        .drop_duplicates()
        .sort_values("canonical_column_order")
    )
    source = source.assign(
        sample_key=source["dataset"] + "::" + source["donor_id"].astype(str)
    )
    donors = donors.assign(
        sample_key=donors["dataset"] + "::" + donors["donor_id"].astype(str)
    )
    matrix = source.pivot(
        index="canonical_row_order",
        columns="sample_key",
        values="family_aCAIrel",
    ).reindex(
        index=row_meta["canonical_row_order"],
        columns=donors["sample_key"],
    )
    require(matrix.shape == (16, 51) and matrix.notna().all().all(), "Rendered matrix is not complete 16 × 51")

    vmin = float(stats["min"])
    vmax = float(stats["max"])
    require(vmin < 0 < 1 < vmax, "Corrected full range must contain negative, zero, one, and >1 values")
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "aligned_asymmetry",
        ["#756F91", "#D9D5E2", "#FFFFFF", "#F3E4BD", "#A98245"],
        N=256,
    )
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM))
    ax_labels = fig.add_axes([0.030, 0.165, 0.245, 0.620])
    ax_dir = fig.add_axes([0.286, 0.165, 0.026, 0.620])
    ax_heat = fig.add_axes([0.318, 0.165, 0.600, 0.620])
    ax_cohort = fig.add_axes([0.318, 0.805, 0.600, 0.055])
    ax_rhythm = fig.add_axes([0.318, 0.785, 0.600, 0.020])
    ax_cbar = fig.add_axes([0.318, 0.080, 0.600, 0.022])

    fig.text(0.030, 0.965, "B", fontsize=11.0, fontweight="bold", ha="left", va="top", color=FG)
    fig.text(
        0.066,
        0.965,
        "Donor-level family aCAIrel across three paired-chamber bulk cohorts",
        fontsize=9.0,
        fontweight="bold",
        ha="left",
        va="top",
        color=FG,
    )
    fig.text(
        0.066,
        0.925,
        "16 frozen NCP families · 51 donors with complete LA–RA pairs · 26 SR · 25 AF",
        fontsize=6.5,
        ha="left",
        va="top",
        color=MUTED,
    )
    fig.text(
        0.910,
        0.925,
        "Columns: cohort → rhythm → donor",
        fontsize=6.5,
        ha="right",
        va="top",
        color=MUTED,
    )

    image = ax_heat.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
        origin="upper",
    )
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    for spine in ax_heat.spines.values():
        spine.set_linewidth(0.4)
        spine.set_color("#A0A0A0")

    cohort_boundaries: list[int] = []
    rhythm_boundaries: list[int] = []
    for index in range(1, len(donors)):
        prev = donors.iloc[index - 1]
        curr = donors.iloc[index]
        if curr["dataset"] != prev["dataset"]:
            cohort_boundaries.append(index)
        elif curr["rhythm"] != prev["rhythm"]:
            rhythm_boundaries.append(index)
    for boundary in rhythm_boundaries:
        ax_heat.axvline(boundary - 0.5, color="#777777", lw=0.35)
    for boundary in cohort_boundaries:
        ax_heat.axvline(boundary - 0.5, color="white", lw=3.2)
        ax_heat.axvline(boundary - 0.5, color="#555555", lw=0.50)

    ax_cohort.set_xlim(-0.5, 50.5)
    ax_cohort.set_ylim(0, 1)
    ax_cohort.axis("off")
    for cohort in COHORT_ORDER:
        positions = np.flatnonzero(donors["dataset"].eq(cohort).to_numpy())
        x0, x1 = positions[0] - 0.5, positions[-1] + 0.5
        ax_cohort.add_patch(
            mpatches.Rectangle(
                (x0, 0),
                x1 - x0,
                1,
                facecolor="#B8B8B8",
                edgecolor="white",
                linewidth=0.5,
            )
        )
        ax_cohort.text(
            (x0 + x1) / 2,
            0.5,
            f"{cohort} (n={len(positions)})",
            ha="center",
            va="center",
            fontsize=6.5,
            fontweight="normal",
            fontstretch="condensed",
            color=FG,
        )

    ax_rhythm.set_xlim(-0.5, 50.5)
    ax_rhythm.set_ylim(0, 1)
    ax_rhythm.axis("off")
    for cohort in COHORT_ORDER:
        for rhythm in RHYTHM_ORDER:
            positions = np.flatnonzero(
                (donors["dataset"].eq(cohort) & donors["rhythm"].eq(rhythm)).to_numpy()
            )
            x0, x1 = positions[0] - 0.5, positions[-1] + 0.5
            ax_rhythm.add_patch(
                mpatches.Rectangle(
                    (x0, 0),
                    x1 - x0,
                    1,
                    facecolor=RHYTHM_COLORS[rhythm],
                    edgecolor="#8A8A8A",
                    linewidth=0.25,
                )
            )
            ax_rhythm.text(
                (x0 + x1) / 2,
                0.5,
                rhythm,
                ha="center",
                va="center",
                fontsize=6.5,
                fontweight="bold",
                color="white" if rhythm == "AF" else FG,
            )

    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(16, 0)
    ax_labels.axis("off")
    ax_labels.text(
        0.98,
        -0.45,
        "Frozen NCP family",
        ha="right",
        va="bottom",
        fontsize=6.5,
        color=MUTED,
    )
    for i, row in enumerate(row_meta.itertuples(index=False)):
        ax_labels.text(
            0.98,
            i + 0.5,
            f"F{int(row.family_id)}  {row.short_label}",
            ha="right",
            va="center",
            fontsize=6.5,
            color=FG,
        )
        if i:
            ax_labels.axhline(i, color="#F1F1F1", lw=0.25)

    ax_dir.set_xlim(0, 1)
    ax_dir.set_ylim(16, 0)
    ax_dir.axis("off")
    ax_dir.text(
        0.5,
        -0.45,
        "",
        ha="center",
        va="bottom",
        fontsize=6.5,
        color=MUTED,
    )
    direction_blocks = [
        (0, 11, "LA", LA_COLOR),
        (11, 16, "RA", RA_COLOR),
    ]
    for start, stop, label, color in direction_blocks:
        ax_dir.add_patch(
            mpatches.Rectangle(
                (0, start),
                1,
                stop - start,
                facecolor=mpl.colors.to_rgba(color, 0.20),
                edgecolor="none",
            )
        )
        ax_dir.text(
            0.5,
            (start + stop) / 2,
            label,
            ha="center",
            va="center",
            fontsize=6.5,
            fontweight="bold",
            color=color,
        )
    for axis in (ax_labels, ax_dir):
        axis.axhline(11, color="#4F4F4F", lw=1.1)
    ax_heat.axhline(10.5, color="#4F4F4F", lw=1.1)

    colorbar = fig.colorbar(image, cax=ax_cbar, orientation="horizontal")
    ticks = [vmin, 0.0, 1.0, vmax]
    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels(
        [
            f"{vmin:.2f}",
            "",
            "",
            f"{vmax:.2f}",
        ]
    )
    colorbar.ax.tick_params(labelsize=6.5, length=2, width=0.4, pad=2)
    colorbar.outline.set_linewidth(0.4)
    colorbar.ax.axvline(norm(0.0), color=FG, lw=0.7, ymin=0, ymax=1)
    colorbar.ax.axvline(norm(1.0), color=FG, lw=0.7, ymin=0, ymax=1)
    fig.text(
        0.318,
        0.118,
        "Aligned family aCAIrel",
        ha="left",
        va="center",
        fontsize=6.5,
        color=FG,
    )

    fig.text(
        0.918,
        0.118,
        "Global scale; no row or cohort standardization",
        ha="right",
        va="center",
        fontsize=6.5,
        color=MUTED,
    )
    zero_x = 0.318 + 0.600 * norm(0.0)
    one_x = 0.318 + 0.600 * norm(1.0)
    fig.text(
        zero_x,
        0.050,
        "0 = chamber symmetry",
        ha="right",
        va="center",
        fontsize=6.5,
        fontweight="bold",
        color=FG,
    )
    fig.text(
        one_x,
        0.050,
        "1 = healthy reference",
        ha="left",
        va="center",
        fontsize=6.5,
        fontweight="bold",
        color=FG,
    )
    fig.text(
        0.318,
        0.018,
        "<0 reversal",
        ha="left",
        va="center",
        fontsize=6.5,
        color=MUTED,
    )
    fig.text(
        0.918,
        0.018,
        ">1 greater-than-reference asymmetry",
        ha="right",
        va="center",
        fontsize=6.5,
        color=MUTED,
    )
    return fig


def validate_outputs(pdf: Path, svg: Path) -> None:
    pdfinfo = shutil.which("pdfinfo")
    pdftotext = shutil.which("pdftotext")
    require(pdfinfo is not None and pdftotext is not None, "PDF QC tools are unavailable")
    info = subprocess.run([pdfinfo, str(pdf)], check=True, capture_output=True, text=True).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read Figure 2B PDF dimensions")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    require(
        abs(width_mm - WIDTH_MM) < 0.05 and abs(height_mm - HEIGHT_MM) < 0.05,
        f"Figure 2B PDF dimensions are {width_mm:.2f} × {height_mm:.2f} mm",
    )
    visible = subprocess.run(
        [pdftotext, "-layout", str(pdf), "-"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    normalized = " ".join(visible.split())
    required = [
        "Donor-level family aCAIrel across three paired-chamber bulk cohorts",
        "16 frozen NCP families",
        "51 donors with complete LA–RA pairs",
        "Columns: cohort → rhythm → donor",
        "GSE128188 (n=10)",
        "GSE79768 (n=13)",
        "GSE115574 (n=28)",
        "Aligned family aCAIrel",
        "Global scale; no row or cohort standardization",
        "chamber symmetry",
        "healthy reference",
        "greater-than-reference asymmetry",
    ]
    require(all(phrase in normalized for phrase in required), "Figure 2B visible text is incomplete")
    forbidden = ["raw aCAI", "30 programs", "row z-score", "cohort centered"]
    require(not any(token in visible for token in forbidden), "Legacy/prohibited Figure 2B text is visible")
    svg_text = svg.read_text(encoding="utf-8")
    require("raw aCAI" not in svg_text and "30 programs" not in svg_text, "Legacy text is present in SVG")


def write_qc(
    manifest: pd.DataFrame,
    source: pd.DataFrame,
    stats: dict[str, float | int | str],
    output_hashes: dict[str, str],
) -> None:
    manifest_lines = "\n".join(
        f"- `{row.artifact_role}`: `{row.artifact_path}` — `{row.sha256}`"
        for row in manifest.itertuples(index=False)
    )
    output_lines = "\n".join(f"- `{path}` — `{digest}`" for path, digest in output_hashes.items())
    OUT_QC.write_text(
        f"""# Figure 2B corrected family-donor heatmap QC

Status: **PASS**

## Matrix and provenance
- matrix: 16 frozen NCP families × 51 paired disease donors ({len(source)} cells)
- rhythm composition: {stats['n_SR']} SR / {stats['n_AF']} AF
- GSE128188: 5 SR / 5 AF
- GSE79768: 6 SR / 7 AF
- GSE115574: 15 SR / 13 AF
- authoritative family values verified against equal-weight corrected constituent-program means
- maximum verification error: {float(stats['max_verification_error']):.3g}
- column order: cohort → rhythm (SR, AF) → natural donor ID
- row order: frozen canonical family display order

## Scale statistics
- global minimum: {float(stats['min']):.8g}
- q01: {float(stats['q01']):.8g}
- q05: {float(stats['q05']):.8g}
- q25: {float(stats['q25']):.8g}
- median: {float(stats['median']):.8g}
- q75: {float(stats['q75']):.8g}
- q95: {float(stats['q95']):.8g}
- q99: {float(stats['q99']):.8g}
- global maximum: {float(stats['max']):.8g}
- clipping policy: none; complete corrected range shown
- semantic anchors: <0 reversal; 0 chamber symmetry; 1 healthy-reference magnitude; >1 greater-than-reference asymmetry

## Prohibited transformations
- upstream sample-level scoring rerun: no
- row z-scoring: no
- cohort centering: no
- donor clustering: no
- imputation: no
- family gene-union scoring: no

## Manifest-resolved authoritative inputs
{manifest_lines}

## Outputs
{output_lines}
""",
        encoding="utf-8",
    )


def main() -> None:
    manifest, paths = load_manifest()
    source, families, stats = build_source(paths)
    OUT_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    OUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    source.to_csv(OUT_SOURCE, sep="\t", index=False)

    fig = render(source, families, stats)
    pdf = OUT_BASE.with_suffix(".pdf")
    svg = OUT_BASE.with_suffix(".svg")
    png = OUT_BASE.with_suffix(".png")
    fig.savefig(pdf, bbox_inches=None, pad_inches=0)
    fig.savefig(svg, bbox_inches=None, pad_inches=0)
    fig.savefig(png, dpi=600, bbox_inches=None, pad_inches=0)
    fig.savefig(OUT_PREVIEW, dpi=150, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    validate_outputs(pdf, svg)

    output_hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (pdf, svg, png, OUT_PREVIEW, OUT_SOURCE)
    }
    write_qc(manifest, source, stats, output_hashes)
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)} ({len(source)} family-donor rows)")
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(
        "Range: "
        f"{float(stats['min']):.4f} to {float(stats['max']):.4f}; "
        f"verification error={float(stats['max_verification_error']):.3g}"
    )
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
