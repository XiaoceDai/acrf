#!/usr/bin/env python3
"""Final Figure 2A: frozen NCP deployment into three paired AF cohorts."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
MANIFEST = RESULTS / "fig2a_input_manifest.tsv"

OUT_BASE = FIGURES / "fig2a_frozen_ncp_bulk_deployment"
OUT_PREVIEW = FIGURES / "fig2a_frozen_ncp_bulk_deployment_preview.png"
OUT_SOURCE = RESULTS / "fig2a_bulk_deployment_source.tsv"
OUT_SUMMARY = RESULTS / "fig2a_bulk_deployment_summary.tsv"
OUT_QC = RESULTS / "fig2a_bulk_deployment_qc.md"

SCORER_PACKAGE = "singscore"
SCORER_VERSION = "1.28.1"
SCORE_METHOD = (
    "singscore::rankGenes(tiesMethod=min)+"
    "simpleScore(upSet=measured_frozen_LE,centerScore=TRUE)"
)
COHORT_ORDER = ["GSE128188", "GSE79768", "GSE115574"]
EXPECTED_COUNTS = {
    "GSE128188": {"SR": 5, "AF": 5},
    "GSE79768": {"SR": 6, "AF": 7},
    "GSE115574": {"SR": 15, "AF": 13},
}
REQUIRED_ROLES = {
    "frozen_program_metadata",
    "frozen_family_metadata",
    "leading_edge_validation",
    "scorer_version_record",
    "healthy_program_reference_summary",
    "disease_program_contrasts",
    "family_acairel",
    "cohort_pairing_summary",
    "GSE128188_sample_metadata",
    "GSE79768_sample_metadata",
    "GSE115574_sample_metadata",
}

MM = 1 / 25.4
WIDTH_MM = 178
HEIGHT_MM = 55
MARGIN_MM = 4
TITLE_H_MM = 11
STRIP_H_MM = 36
GAP_MM = 3.0
COL_W_MM = [32, 37, 58, 34]
LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
FG_COLOR = "#1A1A1A"
MUTED_COLOR = "#595959"
LINE_COLOR = "#E2E2E2"
STRIP_BG = "#F7F7F7"
QUANT_BG = "#F2F2F2"
MID_GREY = "#9A9A9A"
FAINT_COLOR = "#7A7A7A"
KICKER_COLOR = "#8A8A8A"
WHITE = "#FFFFFF"
MIN_FONT_PT = 6.0


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest() -> tuple[pd.DataFrame, dict[str, Path]]:
    require(MANIFEST.is_file(), f"Missing Figure 2A manifest: {MANIFEST}")
    manifest = pd.read_csv(MANIFEST, sep="\t")
    required_columns = {
        "artifact_role",
        "artifact_path",
        "sha256",
        "scorer_package",
        "scorer_version",
        "reference_dataset",
        "allows_cohort_specific_reselection",
        "allows_reference_recalibration",
        "is_authoritative",
    }
    require(required_columns.issubset(manifest.columns), "Figure 2A manifest schema is incomplete")
    require(set(manifest["artifact_role"]) == REQUIRED_ROLES, "Figure 2A manifest roles differ from the approved set")
    require(manifest["artifact_role"].is_unique, "Figure 2A manifest roles are not unique")
    require(set(manifest["scorer_package"].astype(str)) == {SCORER_PACKAGE}, "Unexpected scorer package")
    require(set(manifest["scorer_version"].astype(str)) == {SCORER_VERSION}, "Unexpected scorer version")
    require(set(manifest["reference_dataset"].astype(str)) == {"GSE112339"}, "Reference dataset is not exclusively GSE112339")
    for column in (
        "allows_cohort_specific_reselection",
        "allows_reference_recalibration",
    ):
        require(
            manifest[column].astype(str).str.lower().eq("false").all(),
            f"Manifest permits forbidden operation: {column}",
        )
    require(
        manifest["is_authoritative"].astype(str).str.lower().eq("true").all(),
        "Non-authoritative Figure 2A input detected",
    )

    paths: dict[str, Path] = {}
    allowed_roots = (
        (ROOT / "results" / "ncp_v3_singscore_reanalysis").resolve(),
        (ROOT / "data" / "ncp_v3_singscore_reanalysis").resolve(),
    )
    for row in manifest.itertuples(index=False):
        path = (ROOT / row.artifact_path).resolve()
        require(path.is_file(), f"Manifest input is missing: {path}")
        require(any(path.is_relative_to(base) for base in allowed_roots), f"Input outside corrected namespaces: {path}")
        require(sha256(path) == row.sha256, f"Checksum mismatch: {path}")
        paths[str(row.artifact_role)] = path
    return manifest, paths


def included_metadata(path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(path, sep="\t")
    included = metadata.loc[
        metadata["analysis_inclusion"].astype(str).str.lower().eq("true")
    ].copy()
    included["donor_id"] = included["donor_id"].astype(str)
    return included


def build_source(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float | int | str]]:
    programs = pd.read_csv(paths["frozen_program_metadata"], sep="\t")
    families = pd.read_csv(paths["frozen_family_metadata"], sep="\t")
    contrasts = pd.read_csv(paths["disease_program_contrasts"], sep="\t")
    contrasts["donor_id"] = contrasts["donor_id"].astype(str)
    observed_family = pd.read_csv(paths["family_acairel"], sep="\t")
    observed_family["donor_id"] = observed_family["donor_id"].astype(str)
    reference = pd.read_csv(paths["healthy_program_reference_summary"], sep="\t")
    versions = pd.read_csv(paths["scorer_version_record"], sep="\t")
    leading_edge = pd.read_csv(paths["leading_edge_validation"], sep="\t").iloc[0]

    require(len(programs) == 30 and programs["program_id"].nunique() == 30, "Frozen program count changed")
    require(len(families) == 16 and families["family_id"].nunique() == 16, "Frozen family count changed")
    require(programs["family_id"].nunique() == 16, "Program topology does not contain 16 families")
    require(
        versions.loc[versions["component"] == SCORER_PACKAGE, "version"].astype(str).tolist()
        == [SCORER_VERSION],
        "singscore version is not 1.28.1",
    )
    require(
        int(leading_edge["n_programs"]) == 30
        and int(leading_edge["n_set_identical"]) == 30
        and str(leading_edge["gate_passed"]).lower() == "true",
        "Frozen leading-edge gate failed",
    )
    require(
        len(reference) == 30
        and reference["program_id"].nunique() == 30
        and set(reference["reference_dataset"].astype(str)) == {"GSE112339"},
        "Healthy M_p provenance changed",
    )
    require(len(contrasts) == 1530, "Disease program contrasts are not 51 x 30")
    require(set(contrasts["dataset"]) == set(COHORT_ORDER), "Unexpected cohort in disease contrasts")
    require(set(contrasts["reference_dataset"].astype(str)) == {"GSE112339"}, "Disease denominator was recalibrated")
    require(set(contrasts["score_method"].astype(str)) == {SCORE_METHOD}, "Disease contrasts use an unexpected scorer")
    require(contrasts["program_id"].nunique() == 30, "Disease contrasts do not contain 30 programs")

    donor_rows = []
    for cohort in COHORT_ORDER:
        metadata = included_metadata(paths[f"{cohort}_sample_metadata"])
        donor_chambers = metadata.groupby("donor_id")["chamber"].agg(lambda values: set(values))
        require(donor_chambers.map(lambda value: value == {"LA", "RA"}).all(), f"Incomplete pair in {cohort}")
        donor_rhythm = metadata.groupby("donor_id")["rhythm"].agg(lambda values: set(values))
        require(donor_rhythm.map(len).eq(1).all(), f"Mixed rhythm within donor in {cohort}")
        donors = metadata.drop_duplicates("donor_id")[["dataset", "donor_id", "rhythm"]]
        counts = donors["rhythm"].value_counts().to_dict()
        require(counts == EXPECTED_COUNTS[cohort], f"{cohort} counts changed: {counts}")
        donor_rows.append(donors)
    donor_table = pd.concat(donor_rows, ignore_index=True)
    require(len(donor_table) == 51, "Expected 51 included paired donors")
    require(
        donor_table[["dataset", "donor_id"]].drop_duplicates().shape[0] == 51,
        "Not 51 unique cohort-scoped donors",
    )
    require((donor_table["rhythm"] == "SR").sum() == 26, "SR donor count is not 26")
    require((donor_table["rhythm"] == "AF").sum() == 25, "AF donor count is not 25")

    contrast_keys = contrasts[["dataset", "donor_id", "rhythm"]].drop_duplicates()
    require(
        set(map(tuple, contrast_keys.to_numpy())) == set(map(tuple, donor_table.to_numpy())),
        "Contrast donors do not exactly match included paired metadata",
    )
    per_donor_programs = contrasts.groupby(["dataset", "donor_id"])["program_id"].nunique()
    require((per_donor_programs == 30).all(), "A disease donor is missing a frozen program")

    expected_family_sizes = programs.groupby("family_id")["program_id"].nunique().rename("expected_n_programs")
    recomputed = (
        contrasts.groupby(["dataset", "donor_id", "rhythm", "family_id"], as_index=False)
        .agg(
            family_aCAIrel_recomputed=("aCAIrel", "mean"),
            n_programs_recomputed=("program_id", "nunique"),
        )
        .merge(expected_family_sizes, on="family_id", how="left", validate="many_to_one")
    )
    require(
        (recomputed["n_programs_recomputed"] == recomputed["expected_n_programs"]).all(),
        "Family aggregation omitted a constituent program",
    )
    disease_observed = observed_family.loc[
        observed_family["dataset"].isin(COHORT_ORDER),
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
        ],
    ].copy()
    source = recomputed.merge(
        disease_observed,
        on=["dataset", "donor_id", "rhythm", "family_id"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    require((source["_merge"] == "both").all(), "Observed and recomputed family tables do not align")
    source["absolute_aggregation_error"] = (
        source["family_aCAIrel"] - source["family_aCAIrel_recomputed"]
    ).abs()
    max_error = float(source["absolute_aggregation_error"].max())
    require(max_error <= 1e-12, f"Family aCAIrel is not equal-weight program mean; max error={max_error}")
    require(len(source) == 51 * 16, "Figure 2A source is not 51 donors x 16 families")
    source = source.drop(columns="_merge").sort_values(
        ["dataset", "rhythm", "donor_id", "canonical_family_display_order"]
    )

    summary_rows = []
    for cohort in COHORT_ORDER:
        cohort_donors = donor_table.loc[donor_table["dataset"] == cohort]
        summary_rows.append(
            {
                "scope": cohort,
                "paired_donors": len(cohort_donors),
                "SR_donors": int((cohort_donors["rhythm"] == "SR").sum()),
                "AF_donors": int((cohort_donors["rhythm"] == "AF").sum()),
                "frozen_programs": 30,
                "frozen_families": 16,
                "scorer_package": SCORER_PACKAGE,
                "scorer_version": SCORER_VERSION,
                "reference_dataset": "GSE112339",
                "cohort_specific_reselection": False,
                "reference_recalibration": False,
            }
        )
    summary_rows.append(
        {
            "scope": "total",
            "paired_donors": 51,
            "SR_donors": 26,
            "AF_donors": 25,
            "frozen_programs": 30,
            "frozen_families": 16,
            "scorer_package": SCORER_PACKAGE,
            "scorer_version": SCORER_VERSION,
            "reference_dataset": "GSE112339",
            "cohort_specific_reselection": False,
            "reference_recalibration": False,
        }
    )
    summary_table = pd.DataFrame(summary_rows)
    facts: dict[str, float | int | str] = {
        "n_cohorts": 3,
        "n_paired_donors": 51,
        "n_SR": 26,
        "n_AF": 25,
        "n_programs": 30,
        "n_families": 16,
        "max_aggregation_error": max_error,
        "leading_edge_version": programs["gene_set_version"].iloc[0],
    }
    return source, summary_table, facts


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
    """Return mm geometry for the 178 x 55 mm publication layout."""
    content_w = WIDTH_MM - 2 * MARGIN_MM
    require(
        sum(COL_W_MM) + 3 * GAP_MM == content_w,
        "Column widths and gaps do not fill the content area",
    )
    xs = [MARGIN_MM]
    for width in COL_W_MM[:-1]:
        xs.append(xs[-1] + width + GAP_MM)
    strip_y0 = MARGIN_MM
    strip_y1 = MARGIN_MM + STRIP_H_MM
    title_y1 = HEIGHT_MM - MARGIN_MM
    title_y0 = title_y1 - TITLE_H_MM
    gap_centers = [
        xs[index] + COL_W_MM[index] + GAP_MM / 2
        for index in range(3)
    ]
    return {
        "xs": xs,
        "strip_y0": strip_y0,
        "strip_y1": strip_y1,
        "strip_cy": (strip_y0 + strip_y1) / 2,
        "title_y0": title_y0,
        "title_y1": title_y1,
        "gap_centers": gap_centers,
        "content_right": MARGIN_MM + content_w,
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


def draw_lock_icon(ax: plt.Axes, x: float, y: float, size: float = 2.2) -> None:
    """Low-prominence line-art lock in the module corner."""
    body_w = size * 0.72
    body_h = size * 0.55
    shackle_r = size * 0.18
    body_x = x - body_w
    body_y = y - body_h
    ax.add_patch(
        mpl.patches.Rectangle(
            (body_x, body_y),
            body_w,
            body_h,
            fill=False,
            edgecolor="#B8B8B8",
            linewidth=0.35,
            zorder=3,
        )
    )
    ax.add_patch(
        mpl.patches.Arc(
            (body_x + body_w / 2, body_y + body_h),
            width=body_w * 0.72,
            height=body_h * 1.05,
            angle=0,
            theta1=0,
            theta2=180,
            fill=False,
            edgecolor="#B8B8B8",
            linewidth=0.35,
            zorder=3,
        )
    )
    ax.plot(
        [body_x + body_w / 2 - shackle_r, body_x + body_w / 2 + shackle_r],
        [body_y + body_h * 0.42, body_y + body_h * 0.42],
        color="#B8B8B8",
        lw=0.35,
        zorder=3,
    )


def draw_output_badge(
    ax: plt.Axes,
    x: float,
    y: float,
    letter: str,
    badge_w: float = 3.2,
    badge_h: float = 3.2,
) -> float:
    ax.add_patch(
        mpl.patches.Rectangle(
            (x, y - badge_h / 2),
            badge_w,
            badge_h,
            facecolor=WHITE,
            edgecolor="#C8C8C8",
            linewidth=0.35,
            zorder=2,
        )
    )
    ax.text(
        x + badge_w / 2,
        y,
        letter,
        ha="center",
        va="center",
        fontsize=6.8,
        fontweight="bold",
        color=FG_COLOR,
        zorder=3,
    )
    return badge_w


def add_column_axes(
    fig: plt.Figure,
    x_mm: float,
    y_mm: float,
    w_mm: float,
    h_mm: float,
) -> plt.Axes:
    left = x_mm / WIDTH_MM
    bottom = y_mm / HEIGHT_MM
    width = w_mm / WIDTH_MM
    height = h_mm / HEIGHT_MM
    col_ax = fig.add_axes([left, bottom, width, height])
    col_ax.set_xlim(0, w_mm)
    col_ax.set_ylim(0, h_mm)
    col_ax.axis("off")
    col_ax.set_facecolor(WHITE)
    return col_ax


def draw_column_header(
    col_ax: plt.Axes,
    w_mm: float,
    h_mm: float,
    kicker: str,
    title: str | None,
    pad_x: float,
    pad_y: float,
) -> float:
    header_y = h_mm - pad_y
    col_ax.text(
        pad_x,
        header_y,
        kicker,
        ha="left",
        va="top",
        fontsize=6.0,
        color=KICKER_COLOR,
        clip_on=True,
    )
    if title:
        col_ax.text(
            pad_x,
            header_y - 2.8,
            title,
            ha="left",
            va="top",
            fontsize=7.0,
            fontweight="semibold",
            color=FG_COLOR,
            linespacing=1.05,
            clip_on=True,
        )
        return header_y - 9.0
    return header_y - 4.5


def formula_bbox_bottom(col_ax: plt.Axes, text_artist: mpl.text.Text) -> float:
    """Return the bottom edge (mm) of a rendered text/math artist."""
    fig = col_ax.figure
    fig.canvas.draw()
    bbox = text_artist.get_window_extent().transformed(col_ax.transData.inverted())
    return float(bbox.y0)


def formula_bottom_mm(col_ax: plt.Axes, text_artist: mpl.text.Text, pad_mm: float = 1.6) -> float:
    """Return y (mm) for a divider placed below a rendered formula block."""
    return formula_bbox_bottom(col_ax, text_artist) - pad_mm


def render(summary: pd.DataFrame) -> plt.Figure:
    set_style()
    geom = strip_geometry()
    xs: list[float] = geom["xs"]  # type: ignore[assignment]
    strip_y0: float = geom["strip_y0"]  # type: ignore[assignment]
    strip_y1: float = geom["strip_y1"]  # type: ignore[assignment]
    strip_cy: float = geom["strip_cy"]  # type: ignore[assignment]
    title_y1: float = geom["title_y1"]  # type: ignore[assignment]
    gap_centers: list[float] = geom["gap_centers"]  # type: ignore[assignment]

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM))
    fig.patch.set_facecolor(WHITE)

    pad_x = 2.0
    pad_y = 1.8
    body_size = 6.2

    # Title band on full-figure axes.
    title_ax = fig.add_axes([0, 0, 1, 1])
    title_ax.set_xlim(0, WIDTH_MM)
    title_ax.set_ylim(0, HEIGHT_MM)
    title_ax.axis("off")
    title_ax.set_facecolor(WHITE)
    title_ax.text(
        MARGIN_MM,
        title_y1 - 0.6,
        "A",
        ha="left",
        va="top",
        fontsize=11.0,
        fontweight="bold",
        color=FG_COLOR,
    )
    title_ax.text(
        MARGIN_MM + 7.5,
        title_y1 - 0.6,
        "Frozen NCP deployment in three paired-chamber AF cohorts",
        ha="left",
        va="top",
        fontsize=8.8,
        fontweight="bold",
        color=FG_COLOR,
    )
    title_ax.text(
        MARGIN_MM + 7.5,
        title_y1 - 5.0,
        "Frozen Fig 1 definitions were applied unchanged across cohorts; "
        "no cohort-specific re-selection or recalibration.",
        ha="left",
        va="top",
        fontsize=6.5,
        color=FAINT_COLOR,
        linespacing=1.12,
    )

    panel_faces = [WHITE, WHITE, QUANT_BG, WHITE]
    col_axes = [
        add_column_axes(fig, xs[index], strip_y0, COL_W_MM[index], STRIP_H_MM)
        for index in range(4)
    ]
    for index, col_ax in enumerate(col_axes):
        col_ax.set_facecolor(panel_faces[index])
        col_ax.add_patch(
            mpl.patches.Rectangle(
                (0, 0),
                COL_W_MM[index],
                STRIP_H_MM,
                facecolor="none",
                edgecolor="#E0E0E0",
                linewidth=0.35,
                zorder=0,
            )
        )

    arrow_half = GAP_MM * 0.26
    for gap_x in gap_centers:
        bridge_arrow(title_ax, gap_x - arrow_half, gap_x + arrow_half, strip_cy)

    # Module 1 — INPUT
    w1 = COL_W_MM[0]
    h = STRIP_H_MM
    body_top = draw_column_header(
        col_axes[0], w1, h, "INPUT", "Frozen NCP\nfrom Fig 1", pad_x, pad_y
    )
    col_axes[0].text(
        pad_x,
        body_top,
        "30 programs  →\n16 families",
        fontsize=6.6,
        fontweight="bold",
        color=FG_COLOR,
        ha="left",
        va="top",
        linespacing=1.05,
        clip_on=True,
    )
    module1_lines = [
        "leading-edge gene sets",
        r"canonical direction  $s_p$",
        r"reference magnitude  $M_p$",
        "family topology and\ndisplay order",
    ]
    line_y = body_top - 7.0
    for text in module1_lines:
        col_axes[0].text(
            pad_x,
            line_y,
            text,
            fontsize=body_size,
            color=MUTED_COLOR,
            va="top",
            linespacing=1.04,
            clip_on=True,
        )
        line_y -= 3.0 if "\n" in text else 2.8
    draw_lock_icon(col_axes[0], w1 - 1.0, h - 1.0)

    # Module 2 — COHORTS
    w2 = COL_W_MM[1]
    body_top = draw_column_header(
        col_axes[1], w2, h, "COHORTS", "Three paired\nLA–RA cohorts", pad_x, pad_y
    )
    rows = summary.set_index("scope")
    cohort_x = pad_x
    sr_x = 24.5
    af_x = w2 - pad_x
    table_top = body_top - 0.2
    row_step = 3.5
    for index, cohort in enumerate(COHORT_ORDER):
        row = rows.loc[cohort]
        row_y = table_top - index * row_step
        col_axes[1].text(cohort_x, row_y, cohort, fontsize=body_size, color=FG_COLOR, va="top", clip_on=True)
        col_axes[1].text(
            sr_x,
            row_y,
            f"{int(row['SR_donors'])} SR",
            fontsize=body_size,
            color=MUTED_COLOR,
            ha="right",
            va="top",
            clip_on=True,
        )
        col_axes[1].text(
            af_x,
            row_y,
            f"{int(row['AF_donors'])} AF",
            fontsize=body_size,
            color=MUTED_COLOR,
            ha="right",
            va="top",
            clip_on=True,
        )

    total = rows.loc["total"]
    sep_y = table_top - len(COHORT_ORDER) * row_step - 1.0
    col_axes[1].plot(
        [cohort_x, w2 - pad_x],
        [sep_y, sep_y],
        color="#D8D8D8",
        lw=0.35,
        clip_on=True,
    )
    col_axes[1].text(
        cohort_x,
        sep_y - 2.8,
        f"{int(total['paired_donors'])} paired donors",
        fontsize=6.6,
        fontweight="bold",
        color=FG_COLOR,
        va="top",
        clip_on=True,
    )
    col_axes[1].text(
        cohort_x,
        sep_y - 5.4,
        f"{int(total['SR_donors'])} SR  |  {int(total['AF_donors'])} AF",
        fontsize=body_size,
        color=MUTED_COLOR,
        ha="left",
        va="top",
        clip_on=True,
    )
    col_axes[1].text(
        cohort_x,
        1.0,
        "One complete LA–RA pair\nper donor",
        fontsize=body_size,
        color=MUTED_COLOR,
        ha="left",
        va="bottom",
        linespacing=1.04,
        clip_on=True,
    )

    # Module 3 — QUANTIFICATION
    w3 = COL_W_MM[2]
    body_top = draw_column_header(
        col_axes[2],
        w3,
        h,
        "QUANTIFICATION",
        "Donor-level chamber quantification",
        pad_x,
        pad_y,
    )
    # The title fits on one line in the 58 mm column, recovering vertical
    # space for the formula and note blocks.
    body_top += 2.6
    formula_top = body_top - 0.4
    formula_gap_mm = 1.2
    formula1_text = col_axes[2].text(
        w3 / 2,
        formula_top,
        r"$\mathrm{aCAIrel}_{d,p}="
        r"\dfrac{s_p\left(S^{LA}_{d,p}-S^{RA}_{d,p}\right)}{M_p}$",
        fontsize=6.6,
        color=FG_COLOR,
        ha="center",
        va="top",
        clip_on=True,
    )
    formula2_top = formula_bbox_bottom(col_axes[2], formula1_text) - formula_gap_mm
    formula2_text = col_axes[2].text(
        w3 / 2,
        formula2_top,
        r"$\mathrm{aCAIrel}_{d,f}="
        r"\dfrac{1}{|P_f|}\sum_{p\in P_f}\mathrm{aCAIrel}_{d,p}$",
        fontsize=6.6,
        color=FG_COLOR,
        ha="center",
        va="top",
        clip_on=True,
    )
    divider_y = formula_bottom_mm(col_axes[2], formula2_text)
    col_axes[2].plot(
        [pad_x, w3 - pad_x],
        [divider_y, divider_y],
        color="#D8D8D8",
        lw=0.35,
        clip_on=True,
    )
    notes = [
        "Same singscore implementation",
        r"$M_p$ fixed from GSE112339",
        "Equal-weight program-first aggregation",
        "No disease-cohort re-selection or recalibration",
    ]
    note_size = 6.0
    note_y = divider_y - 1.4
    for text in notes:
        col_axes[2].text(pad_x, note_y, "•", fontsize=note_size, color=MID_GREY, va="top", clip_on=True)
        col_axes[2].text(
            pad_x + 1.8,
            note_y,
            text,
            fontsize=note_size,
            color=MUTED_COLOR,
            va="top",
            linespacing=1.04,
            clip_on=True,
        )
        note_y -= 2.3

    # Module 4 — OUTPUTS
    w4 = COL_W_MM[3]
    body_top = draw_column_header(col_axes[3], w4, h, "OUTPUTS", None, pad_x, pad_y)
    outputs = [
        ("B", "Donor-level family\naCAIrel"),
        ("C", "Cohort-specific\nrewiring states"),
        ("D", "AF–SR effects and\nconsistency"),
    ]
    row_centers = [body_top - 1.5, body_top - 9.5, body_top - 17.5]
    for index, (letter, label) in enumerate(outputs):
        row_y = row_centers[index]
        if index:
            col_axes[3].plot(
                [pad_x, w4 - pad_x],
                [row_y + 3.6, row_y + 3.6],
                color="#D8D8D8",
                lw=0.35,
                clip_on=True,
            )
        badge_end = draw_output_badge(col_axes[3], pad_x, row_y, letter)
        col_axes[3].text(
            pad_x + badge_end + 1.0,
            row_y + 1.4,
            label,
            fontsize=body_size,
            color=MUTED_COLOR,
            va="top",
            linespacing=1.04,
            clip_on=True,
        )

    return fig


def write_qc(
    manifest: pd.DataFrame,
    source: pd.DataFrame,
    summary: pd.DataFrame,
    facts: dict[str, float | int | str],
    output_hashes: dict[str, str],
) -> None:
    manifest_lines = "\n".join(
        f"- `{row.artifact_role}`: `{row.artifact_path}` — `{row.sha256}`"
        for row in manifest.itertuples(index=False)
    )
    output_lines = "\n".join(f"- `{path}` — `{digest}`" for path, digest in output_hashes.items())
    script_text = Path(__file__).read_text(encoding="utf-8")
    forbidden_fragments = [
        "fig2_primary_" + "acai",
        "fig1b_baseline_" + "donor_aCAI",
        "innovation_ncp_v3_" + "acai_long",
        "fig2a_" + "schematic",
    ]
    require(not any(fragment in script_text for fragment in forbidden_fragments), "Legacy Figure 2A dependency detected")
    font_sizes = [
        float(value)
        for value in re.findall(r"fontsize\s*=\s*([0-9]+(?:\.[0-9]+)?)", script_text)
    ]
    require(font_sizes and min(font_sizes) >= MIN_FONT_PT, "Publication-facing text below 6 pt detected")
    require(len(source) == 816, "Unexpected Figure 2A source row count")
    require(len(summary) == 4, "Unexpected Figure 2A summary row count")
    OUT_QC.write_text(
        f"""# Figure 2A frozen bulk deployment QC

Status: **PASS**

## Hard assertions
- paired disease cohorts: {facts['n_cohorts']}
- complete paired donors: {facts['n_paired_donors']}
- rhythm composition: {facts['n_SR']} SR / {facts['n_AF']} AF
- frozen programs: {facts['n_programs']}
- frozen families: {facts['n_families']}
- scorer: {SCORER_PACKAGE} {SCORER_VERSION}
- score method: `{SCORE_METHOD}`
- healthy magnitude source: GSE112339 only
- maximum equal-weight family aggregation error: {float(facts['max_aggregation_error']):.3g}
- frozen leading-edge version: `{facts['leading_edge_version']}`
- cohort-specific pathway re-selection: no
- disease-cohort reference recalibration: no
- legacy numerical table read: no
- panel dimensions: {WIDTH_MM} x {HEIGHT_MM} mm
- minimum publication-facing text size: {MIN_FONT_PT} pt
- visible title, subtitle and approved module content: checked
- visible internal paths or `NCP_v3`: none

## Corrected authoritative inputs
{manifest_lines}

## Outputs
{output_lines}
""",
        encoding="utf-8",
    )


def validate_outputs(pdf: Path, svg: Path) -> None:
    pdfinfo = shutil.which("pdfinfo")
    pdftotext = shutil.which("pdftotext")
    require(pdfinfo is not None and pdftotext is not None, "PDF QC tools are unavailable")
    info = subprocess.run(
        [pdfinfo, str(pdf)], check=True, capture_output=True, text=True
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read Figure 2A PDF dimensions")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    require(
        abs(width_mm - WIDTH_MM) < 0.05 and abs(height_mm - HEIGHT_MM) < 0.05,
        f"Figure 2A PDF dimensions are {width_mm:.2f} x {height_mm:.2f} mm",
    )
    visible = subprocess.run(
        [pdftotext, "-layout", str(pdf), "-"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    normalized = " ".join(visible.split())
    required_phrases = [
        "Frozen NCP deployment in three paired-chamber AF cohorts",
        "Frozen Fig 1 definitions were applied unchanged across cohorts",
        "no cohort-specific re-selection or recalibration.",
        "leading-edge gene sets",
        "reference magnitude",
        "family topology and",
        "display order",
        "Three paired",
        "LA–RA cohorts",
        "One complete LA–RA pair",
        "per donor",
        "51 paired donors",
        "Donor-level chamber",
        "quantification",
        "Same singscore implementation",
        "fixed from GSE112339",
        "Equal-weight program-first aggregation",
        "No disease-cohort re-selection",
        "or recalibration",
        "Donor-level family",
        "aCAIrel",
        "Cohort-specific",
        "rewiring states",
        "AF–SR effects and",
        "consistency",
        "INPUT",
        "COHORTS",
        "QUANTIFICATION",
        "OUTPUTS",
    ]
    require(all(phrase in normalized for phrase in required_phrases), "Figure 2A visible text is incomplete")
    forbidden = [
        "NCP_v3",
        "results/",
        "scripts/",
        "STAGE 1",
        "NES",
        "REML",
        "Hartung",
        "leave-one-cohort-out",
    ]
    require(not any(token in visible for token in forbidden), "Forbidden legacy/internal text is visible")
    svg_text = svg.read_text(encoding="utf-8")
    require(
        "NCP_v3" not in svg_text and "results/" not in svg_text and "scripts/" not in svg_text,
        "Figure 2A SVG exposes internal nomenclature or paths",
    )


def main() -> None:
    manifest, paths = load_manifest()
    source, summary, facts = build_source(paths)

    OUT_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    OUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    source.to_csv(OUT_SOURCE, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    fig = render(summary)
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
        for path in (pdf, svg, png, OUT_PREVIEW, OUT_SOURCE, OUT_SUMMARY)
    }
    write_qc(manifest, source, summary, facts, output_hashes)
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)} ({len(source)} donor-family rows)")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
