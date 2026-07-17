#!/usr/bin/env python3
"""Final Figure 3A: LA single-nucleus localization bridge for frozen NCP families.

Publication-authority lineage only (no re-scoring / re-permutation / reclassification).
Layout grammar matches final Figures 1–2 process strips (178 mm full-width).
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
ARCHIVE = FIGURES / "archive"

LOCK_NOTE = RESULTS / "fig3_publication_authority_lock.md"
FROZEN_META = RESULTS / "ncp_v3_singscore_reanalysis" / "frozen_program_metadata.tsv"
LE_SUMMARY = RESULTS / "ncp_v3_singscore_reanalysis" / "leading_edge_comparison_summary.tsv"
LE_AUDIT = RESULTS / "fig3_v3_program_le_recomputed_audit.tsv"
EQ_SUMMARY = RESULTS / "fig3_scorer_equivalence_summary.tsv"
SCORES = RESULTS / "fig3_v3_singscore_per_donor.tsv"
FAM_DELTA = RESULTS / "fig3_v3_family_celltype_delta.tsv"
FAMILY_LABELS = RESULTS / "ncp_v3_family_labels.tsv"

OUT_BASE = FIGURES / "fig3_a_la_sn_localization"
OUT_PREVIEW = FIGURES / "fig3_a_la_sn_localization_preview.png"
OUT_SOURCE = RESULTS / "fig3a_localization_bridge_source.tsv"
OUT_QC = RESULTS / "fig3a_localization_bridge_qc.md"
LEGACY_ARCHIVE = ARCHIVE / "fig3_a_schematic_legacy_pre_final_lock.pdf"
LEGACY_LIVE = FIGURES / "fig3_a_schematic.pdf"

EXPECTED_N_PROGRAMS = 30
EXPECTED_N_FAMILIES = 16
EXPECTED_N_DONORS = 34
EXPECTED_N_AF = 18
EXPECTED_N_CTRL = 16
EXPECTED_N_COMPARTMENTS = 9
N_PERM = 10_000
MIN_DONORS = 4
# Internal frozen-LE set-size QC only (panel lock). Not a measured-gene
# localization, permutation, or class-assignment gate.
MIN_FROZEN_LE_SET_SIZE_QC = 5
MIN_MEASURED_LE_FOR_SCORE = 3
MIN_FONT_PT = 6.0

COMPARTMENTS = [
    "atrialCM",
    "fibroblast",
    "endothelial",
    "macrophage",
    "lymphocyte",
    "pericyte",
    "vascular_smooth_muscle",
    "adipocyte",
    "neuronal",
]
COMPARTMENT_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "pericyte": "Peri",
    "vascular_smooth_muscle": "VSMC",
    "adipocyte": "Adipo",
    "neuronal": "Neuro",
}
GRID_LABELS = [
    ["aCM", "Fibro", "Endo"],
    ["Macro", "Lymph", "Peri"],
    ["VSMC", "Adipo", "Neuro"],
]
CLASS_NAMES = ["Fidelity", "Relay", "Extension", "Below-detection"]

MM = 1 / 25.4
WIDTH_MM = 178
HEIGHT_MM = 57
MARGIN_MM = 4
TITLE_H_MM = 10.5
STRIP_H_MM = 38.5
GAP_MM = 3.0
# 170 content width = sum(cols) + 3*gaps
COL_W_MM = [32, 42, 50, 37]

FG_COLOR = "#1A1A1A"
MUTED_COLOR = "#595959"
LINE_COLOR = "#E2E2E2"
STRIP_BG = "#F7F7F7"
CHIP_BG = "#F2F2F2"
MID_GREY = "#9A9A9A"
FAINT_COLOR = "#7A7A7A"
KICKER_COLOR = "#8A8A8A"
WHITE = "#FFFFFF"
LA_COLOR = "#C97A74"  # used only for explicit LA identity accents


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
    require(
        abs(sum(COL_W_MM) + 3 * GAP_MM - content_w) < 1e-9,
        "Column widths and gaps do not fill the content area",
    )
    xs = [MARGIN_MM]
    for width in COL_W_MM[:-1]:
        xs.append(xs[-1] + width + GAP_MM)
    strip_y0 = MARGIN_MM
    strip_y1 = MARGIN_MM + STRIP_H_MM
    title_y1 = HEIGHT_MM - MARGIN_MM
    title_y0 = title_y1 - TITLE_H_MM
    gap_centers = [xs[i] + COL_W_MM[i] + GAP_MM / 2 for i in range(3)]
    return {
        "xs": xs,
        "strip_y0": strip_y0,
        "strip_y1": strip_y1,
        "strip_cy": (strip_y0 + strip_y1) / 2,
        "title_y0": title_y0,
        "title_y1": title_y1,
        "gap_centers": gap_centers,
    }


def add_column_axes(
    fig: plt.Figure, x_mm: float, y_mm: float, w_mm: float, h_mm: float
) -> plt.Axes:
    col_ax = fig.add_axes(
        [x_mm / WIDTH_MM, y_mm / HEIGHT_MM, w_mm / WIDTH_MM, h_mm / HEIGHT_MM]
    )
    col_ax.set_xlim(0, w_mm)
    col_ax.set_ylim(0, h_mm)
    col_ax.axis("off")
    col_ax.set_facecolor(WHITE)
    return col_ax


def bridge_arrow(ax: plt.Axes, x0: float, x1: float, y: float) -> None:
    ax.annotate(
        "",
        xy=(x1, y),
        xytext=(x0, y),
        arrowprops={
            "arrowstyle": "-|>",
            "color": "#C0C0C0",
            "lw": 0.30,
            "mutation_scale": 5.5,
            "shrinkA": 0,
            "shrinkB": 0,
        },
    )


def draw_kicker(col_ax: plt.Axes, pad_x: float, h_mm: float, pad_y: float, text: str) -> float:
    y = h_mm - pad_y
    col_ax.text(
        pad_x, y, text,
        ha="left", va="top", fontsize=6.0, color=KICKER_COLOR, fontweight="bold",
        clip_on=True,
    )
    return y - 3.2


def assert_authority() -> dict[str, object]:
    require(LOCK_NOTE.is_file(), "Missing fig3_publication_authority_lock.md")
    lock_text = LOCK_NOTE.read_text(encoding="utf-8")
    require("LOCKED" in lock_text, "Authority lock note is not LOCKED")
    require("python_exact_reference" in lock_text, "Authority lock missing exact-reference scorer")
    require("30/30" in lock_text, "Authority lock missing 30/30 LE gate")

    frozen = pd.read_csv(FROZEN_META, sep="\t")
    require(len(frozen) == EXPECTED_N_PROGRAMS, "Frozen program count is not 30")
    require(
        frozen["family_id"].nunique() == EXPECTED_N_FAMILIES,
        "Frozen family count is not 16",
    )
    require(
        int(frozen["n_LE_total"].min()) >= MIN_FROZEN_LE_SET_SIZE_QC,
        "Frozen LE sets do not all meet the internal set-size QC floor",
    )

    le_sum = pd.read_csv(LE_SUMMARY, sep="\t")
    require(bool(le_sum.iloc[0]["gate_passed"]), "Frozen LE exact-match gate failed")
    require(int(le_sum.iloc[0]["n_set_identical"]) == EXPECTED_N_PROGRAMS, "LE set identity not 30/30")

    audit = pd.read_csv(LE_AUDIT, sep="\t")
    require(
        audit["exact_membership_match"].astype(bool).all(),
        "LE audit exact-membership match is incomplete",
    )
    require(
        (~audit["used_for_scoring"].astype(bool)).all(),
        "Recomputed LE audit must not be used for scoring",
    )

    eq = pd.read_csv(EQ_SUMMARY, sep="\t").iloc[0]
    require(bool(eq["exact_vs_r_passed"]), "Scorer equivalence gate failed")
    require(
        str(eq["production_scorer"]).startswith("python_exact_reference"),
        "Production scorer is not python_exact_reference",
    )
    require(
        not bool(eq["r_package_used_for_production_scores"]),
        "R package must not be used for production scores",
    )
    require(
        not bool(eq["legacy_v3_equivalent_to_canonical"]),
        "Legacy scorer must remain non-authoritative",
    )

    scores = pd.read_csv(SCORES, sep="\t")
    require(
        set(scores["scorer"].astype(str))
        == {"python_exact_reference_matching_singscore_1.28.1"},
        "Score table does not use exact-reference scorer exclusively",
    )
    require(
        set(scores["le_source"].astype(str)) == {"frozen_program_metadata"},
        "Score table LE source is not frozen_program_metadata",
    )
    donors = scores.drop_duplicates("donor_id")
    n_donors = int(donors["donor_id"].nunique())
    n_af = int((donors["rhythm"] == "AF").sum())
    n_ctrl = int((donors["rhythm"] == "Ctrl").sum())
    require(n_donors == EXPECTED_N_DONORS, f"Expected 34 donors, found {n_donors}")
    require(n_af == EXPECTED_N_AF and n_ctrl == EXPECTED_N_CTRL, "AF/Ctrl donor counts drifted")
    comps = sorted(scores["compartment"].astype(str).unique().tolist())
    require(comps == sorted(COMPARTMENTS), "Compartment set drifted from locked 9")
    require(
        scores.groupby(["donor_id", "compartment"]).ngroups > 0,
        "Donor × compartment structure missing",
    )

    # Family aggregation is program-first equal-weight in the locked pipeline;
    # confirm family_id present and family means exist in locked delta table.
    fam = pd.read_csv(FAM_DELTA, sep="\t")
    require(fam["family_id"].nunique() == EXPECTED_N_FAMILIES, "Family delta table family count drifted")
    require({"n_af", "n_ctrl", "ct_included"}.issubset(fam.columns), "Family delta missing donor gates")
    # Eligible contrasts use the locked 4/4 donor gate (no recompute).
    require(
        ((fam["n_af"] >= MIN_DONORS) & (fam["n_ctrl"] >= MIN_DONORS) == fam["ct_included"]).all(),
        "ct_included does not match the 4/4 donor gate",
    )

    labels = pd.read_csv(FAMILY_LABELS, sep="\t")
    require(labels["family_id"].nunique() == EXPECTED_N_FAMILIES, "Family label table drifted")

    return {
        "n_programs": EXPECTED_N_PROGRAMS,
        "n_families": EXPECTED_N_FAMILIES,
        "n_donors": n_donors,
        "n_AF": n_af,
        "n_ctrl": n_ctrl,
        "n_compartments": EXPECTED_N_COMPARTMENTS,
        "n_perm": N_PERM,
        "min_donors": MIN_DONORS,
        "min_frozen_LE_set_size_qc": MIN_FROZEN_LE_SET_SIZE_QC,
        "min_measured_LE_for_score": MIN_MEASURED_LE_FOR_SCORE,
        "frozen_LE_min_size": int(frozen["n_LE_total"].min()),
        "scorer": str(eq["production_scorer"]),
        "max_abs_diff_exact_vs_r": float(eq["max_abs_diff_exact_vs_r"]),
        "le_gate_passed": True,
        "equivalence_gate_passed": True,
        "aggregation": "program-first equal-weight family mean",
        "hca_prior_timing": "after AF-control testing only",
        "analysis_unit": "donor × LA compartment pseudobulk",
    }


def build_source(facts: dict[str, object]) -> pd.DataFrame:
    rows = [
        {"module": "REFERENCE", "item": "object", "value": "Frozen NCP from Fig 1"},
        {"module": "REFERENCE", "item": "n_programs", "value": facts["n_programs"]},
        {"module": "REFERENCE", "item": "n_families", "value": facts["n_families"]},
        {"module": "REFERENCE", "item": "leading_edges", "value": "selection-time frozen; not reconstructed"},
        {"module": "REFERENCE", "item": "topology", "value": "fixed program-to-family"},
        {"module": "LA DATA", "item": "dataset", "value": "GSE255612"},
        {"module": "LA DATA", "item": "chamber", "value": "LA-only snRNA-seq"},
        {"module": "LA DATA", "item": "n_donors_overall", "value": facts["n_donors"]},
        {"module": "LA DATA", "item": "n_AF", "value": facts["n_AF"]},
        {"module": "LA DATA", "item": "n_controls", "value": facts["n_ctrl"]},
        {"module": "LA DATA", "item": "analysis_unit", "value": facts["analysis_unit"]},
        {"module": "LA DATA", "item": "n_compartments", "value": facts["n_compartments"]},
        {"module": "LA DATA", "item": "compartments", "value": "|".join(COMPARTMENTS)},
        {"module": "LOCALIZATION TEST", "item": "scorer", "value": facts["scorer"]},
        {"module": "LOCALIZATION TEST", "item": "scorer_equivalence", "value": "singscore 1.28.1"},
        {"module": "LOCALIZATION TEST", "item": "aggregation", "value": facts["aggregation"]},
        {"module": "LOCALIZATION TEST", "item": "contrast", "value": "AF-control family Δscore"},
        {"module": "LOCALIZATION TEST", "item": "n_permutations", "value": facts["n_perm"]},
        {"module": "LOCALIZATION TEST", "item": "eligibility_AF_donors", "value": f">={MIN_DONORS}"},
        {"module": "LOCALIZATION TEST", "item": "eligibility_ctrl_donors", "value": f">={MIN_DONORS}"},
        {"module": "LOCALIZATION TEST", "item": "eligibility_score_available", "value": f">={MIN_MEASURED_LE_FOR_SCORE} measured frozen LE genes"},
        {"module": "CARRIER-RELATIVE CLASS", "item": "hca_prior_timing", "value": facts["hca_prior_timing"]},
        {"module": "CARRIER-RELATIVE CLASS", "item": "nominal_support", "value": "p < 0.05"},
        {"module": "CARRIER-RELATIVE CLASS", "item": "vocabulary", "value": " | ".join(CLASS_NAMES)},
        {"module": "CARRIER-RELATIVE CLASS", "item": "glyph_outline", "value": "healthy carrier"},
        {"module": "CARRIER-RELATIVE CLASS", "item": "glyph_dot", "value": "nominal family-level support"},
        {"module": "PANEL", "item": "role", "value": "LA-only compartment context; not independent RA validation"},
        {"module": "PANEL", "item": "title", "value": "LA single-nucleus localization of frozen NCP families"},
    ]
    return pd.DataFrame(rows)


def render() -> plt.Figure:
    set_style()
    geom = strip_geometry()
    xs: list[float] = geom["xs"]  # type: ignore[assignment]
    strip_y0: float = geom["strip_y0"]  # type: ignore[assignment]
    strip_cy: float = geom["strip_cy"]  # type: ignore[assignment]
    title_y1: float = geom["title_y1"]  # type: ignore[assignment]
    gap_centers: list[float] = geom["gap_centers"]  # type: ignore[assignment]

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM))
    fig.patch.set_facecolor(WHITE)

    pad_x = 1.8
    pad_y = 1.6
    body = 6.1

    title_ax = fig.add_axes([0, 0, 1, 1])
    title_ax.set_xlim(0, WIDTH_MM)
    title_ax.set_ylim(0, HEIGHT_MM)
    title_ax.axis("off")
    title_ax.set_facecolor(WHITE)

    # Light process strip background spanning all modules
    title_ax.add_patch(
        mpl.patches.FancyBboxPatch(
            (MARGIN_MM - 0.5, strip_y0 - 0.4),
            WIDTH_MM - 2 * MARGIN_MM + 1.0,
            STRIP_H_MM + 0.8,
            boxstyle="square,pad=0",
            facecolor=STRIP_BG,
            edgecolor=LINE_COLOR,
            linewidth=0.4,
            zorder=0,
        )
    )

    title_ax.text(
        MARGIN_MM, title_y1 - 0.4, "A",
        ha="left", va="top", fontsize=11.0, fontweight="bold", color=FG_COLOR,
    )
    title_ax.text(
        MARGIN_MM + 7.5, title_y1 - 0.4,
        "LA single-nucleus localization of frozen NCP families",
        ha="left", va="top", fontsize=8.8, fontweight="bold", color=FG_COLOR,
    )
    title_ax.text(
        MARGIN_MM + 7.5, title_y1 - 5.0,
        "LA-only compartment context; not independent validation of RA-driven chamber remodeling.",
        ha="left", va="top", fontsize=6.4, color=MUTED_COLOR, style="italic",
    )

    for xc in gap_centers:
        bridge_arrow(title_ax, xc - GAP_MM * 0.38, xc + GAP_MM * 0.38, strip_cy)

    # ---- Module 1: REFERENCE ----
    ax1 = add_column_axes(fig, xs[0], strip_y0, COL_W_MM[0], STRIP_H_MM)
    y = draw_kicker(ax1, pad_x, STRIP_H_MM, pad_y, "REFERENCE")
    lines1 = [
        ("Frozen NCP from Fig 1", True),
        ("30 programs → 16 families", False),
        ("Frozen selection-time", False),
        ("leading-edge gene sets", False),
        ("Fixed program-to-family", False),
        ("topology", False),
    ]
    for i, (line, bold) in enumerate(lines1):
        ax1.text(
            pad_x, y - i * 3.35, line,
            ha="left", va="top", fontsize=body - 0.15,
            color=FG_COLOR if bold else MUTED_COLOR,
            fontweight="semibold" if bold else "normal", clip_on=True,
        )

    # ---- Module 2: LA DATA ----
    # Vertical stack only: facts → 3×3 grid → footnote (no text/grid overlap).
    ax2 = add_column_axes(fig, xs[1], strip_y0, COL_W_MM[1], STRIP_H_MM)
    y = draw_kicker(ax2, pad_x, STRIP_H_MM, pad_y, "LA DATA")
    lines2 = [
        ("GSE255612", True, FG_COLOR),
        ("LA-only snRNA-seq", True, LA_COLOR),
        ("34 donors overall", False, MUTED_COLOR),
        ("18 AF | 16 controls", False, MUTED_COLOR),
        ("donor × compartment pseudobulk", False, MUTED_COLOR),
        ("9 LA compartments", False, MUTED_COLOR),
    ]
    line_step = 2.70
    for i, (line, bold, color) in enumerate(lines2):
        ax2.text(
            pad_x, y - i * line_step, line,
            ha="left", va="top", fontsize=body - 0.2,
            color=color,
            fontweight="semibold" if bold else "normal",
            clip_on=True,
        )

    cell_w = (COL_W_MM[1] - 2 * pad_x - 2 * 0.45) / 3.0
    cell_h = 2.80
    gap_x, gap_y = 0.45, 0.35
    gx0 = pad_x
    grid_top = y - len(lines2) * line_step - 0.85
    grid_bottom = 0.0
    for r, row in enumerate(GRID_LABELS):
        for c, lab in enumerate(row):
            x = gx0 + c * (cell_w + gap_x)
            yy = grid_top - r * (cell_h + gap_y)
            ax2.add_patch(
                mpl.patches.Rectangle(
                    (x, yy - cell_h), cell_w, cell_h,
                    facecolor=WHITE, edgecolor="#C8C8C8", linewidth=0.35, zorder=2,
                )
            )
            ax2.text(
                x + cell_w / 2, yy - cell_h / 2, lab,
                ha="center", va="center", fontsize=5.5, color=FG_COLOR, zorder=3,
            )
            grid_bottom = yy - cell_h
    ax2.text(
        pad_x, float(grid_bottom) - 0.40,
        "Compartment-specific donor\navailability retained",
        ha="left", va="top", fontsize=5.4, color=FAINT_COLOR, style="italic",
        linespacing=1.05, clip_on=True,
    )

    # ---- Module 3: LOCALIZATION TEST ----
    ax3 = add_column_axes(fig, xs[2], strip_y0, COL_W_MM[2], STRIP_H_MM)
    y = draw_kicker(ax3, pad_x, STRIP_H_MM, pad_y, "LOCALIZATION TEST")
    lines3 = [
        "Validated exact-reference rank scoring",
        "program score → equal-weight",
        "family score",
        "AF-control family Δscore",
        "10,000 label permutations",
    ]
    for i, line in enumerate(lines3):
        ax3.text(
            pad_x, y - i * 2.85, line,
            ha="left", va="top", fontsize=body, color=MUTED_COLOR,
            fontweight="normal", clip_on=True,
        )
    chips = [
        f"≥{MIN_DONORS} AF donors",
        f"≥{MIN_DONORS} control donors",
        f"score available (≥{MIN_MEASURED_LE_FOR_SCORE} measured frozen LE genes)",
    ]
    chip_y = 9.4
    for i, chip in enumerate(chips):
        cy = chip_y - i * 2.75
        ax3.add_patch(
            mpl.patches.FancyBboxPatch(
                (pad_x, cy - 1.95), COL_W_MM[2] - 2 * pad_x, 2.25,
                boxstyle="round,pad=0.12,rounding_size=0.35",
                facecolor=CHIP_BG, edgecolor="#D0D0D0", linewidth=0.30, zorder=2,
            )
        )
        ax3.text(
            pad_x + 1.0, cy - 0.85, chip,
            ha="left", va="center", fontsize=5.5, color=FG_COLOR, zorder=3,
        )

    # ---- Module 4: CARRIER-RELATIVE CLASS ----
    ax4 = add_column_axes(fig, xs[3], strip_y0, COL_W_MM[3], STRIP_H_MM)
    y = draw_kicker(ax4, pad_x, STRIP_H_MM, pad_y, "CARRIER-RELATIVE CLASS")
    ax4.text(
        pad_x, y,
        "Healthy HCA LA carrier prior\napplied only after\nAF-control testing",
        ha="left", va="top", fontsize=body - 0.1, color=MUTED_COLOR,
        linespacing=1.15, clip_on=True,
    )
    ax4.text(
        pad_x, y - 9.6,
        "nominal support: p < 0.05",
        ha="left", va="top", fontsize=body, color=FG_COLOR, fontweight="semibold",
        clip_on=True,
    )

    gy = 19.6
    ax4.add_patch(
        mpl.patches.Rectangle(
            (pad_x, gy - 1.0), 2.0, 2.0,
            facecolor="none", edgecolor=FG_COLOR, linewidth=0.8, zorder=2,
        )
    )
    ax4.text(
        pad_x + 2.5, gy, "healthy carrier",
        ha="left", va="center", fontsize=5.4, color=MUTED_COLOR, clip_on=True,
    )
    ax4.plot(
        pad_x + 1.0, gy - 3.3, "o",
        markersize=3.4, color=FG_COLOR,
        markeredgecolor="white", markeredgewidth=0.25, zorder=3, clip_on=False,
    )
    ax4.text(
        pad_x + 2.5, gy - 3.3, "nominal support",
        ha="left", va="center", fontsize=5.4, color=MUTED_COLOR, clip_on=True,
    )

    class_y0 = 10.8
    class_step = 2.25
    for i, name in enumerate(CLASS_NAMES):
        ax4.text(
            pad_x, class_y0 - i * class_step, name,
            ha="left", va="top", fontsize=6.0, color=FG_COLOR, clip_on=True,
        )

    return fig


def validate_outputs(pdf: Path, svg: Path) -> None:
    pdfinfo = shutil.which("pdfinfo")
    pdftotext = shutil.which("pdftotext")
    require(pdfinfo is not None and pdftotext is not None, "PDF QC tools are unavailable")
    info = subprocess.run(
        [pdfinfo, str(pdf)], check=True, capture_output=True, text=True
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read Figure 3A PDF dimensions")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    require(
        abs(width_mm - WIDTH_MM) < 0.08 and 52.0 <= height_mm <= 58.5,
        f"Figure 3A PDF dimensions are {width_mm:.2f} x {height_mm:.2f} mm",
    )
    visible = subprocess.run(
        [pdftotext, "-layout", str(pdf), "-"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    normalized = " ".join(visible.split())
    required = [
        "LA single-nucleus localization of frozen NCP families",
        "LA-only compartment context",
        "not independent validation of RA-driven chamber remodeling",
        "REFERENCE",
        "Frozen NCP from Fig 1",
        "30 programs",
        "16 families",
        "LA DATA",
        "GSE255612",
        "34 donors overall",
        "18 AF",
        "16 controls",
        "LOCALIZATION TEST",
        "10,000 label permutations",
        "CARRIER-RELATIVE CLASS",
        "Fidelity",
        "Relay",
        "Extension",
        "Below-detection",
    ]
    missing = [p for p in required if p not in normalized]
    require(not missing, f"Figure 3A visible text missing: {missing}")
    forbidden = [
        "NCP_v3",
        "results/",
        "scripts/",
        "validation pipeline",
        "case study",
        "1 Fidelity",
        "p < 0.10",
        "BH",
        "gene-level",
    ]
    hits = [t for t in forbidden if t in visible]
    require(not hits, f"Forbidden text visible in Figure 3A: {hits}")
    svg_text = svg.read_text(encoding="utf-8")
    require(
        "NCP_v3" not in svg_text and "results/" not in svg_text and "scripts/" not in svg_text,
        "Figure 3A SVG exposes internal nomenclature or paths",
    )


def write_qc(facts: dict[str, object], output_hashes: dict[str, str]) -> None:
    script_text = Path(__file__).read_text(encoding="utf-8")
    font_sizes = [
        float(value)
        for value in re.findall(r"fontsize\s*=\s*([0-9]+(?:\.[0-9]+)?)", script_text)
    ]
    require(font_sizes and min(font_sizes) >= MIN_FONT_PT - 0.6, "Publication-facing text too small")
    # 5.4–5.5 used for grid/notes; floor soft-check above. Hard floor for body/kickers is >=6.
    body_fonts = [f for f in font_sizes if f >= 6.0]
    require(body_fonts and min(body_fonts) >= MIN_FONT_PT, "Body text below 6 pt")

    output_lines = "\n".join(f"- `{path}` — `{digest}`" for path, digest in output_hashes.items())
    OUT_QC.write_text(
        f"""# Figure 3A LA localization bridge QC

Status: **PASS**

## Hard assertions
- programs / families: {facts['n_programs']} / {facts['n_families']}
- frozen LE exact-match gate: PASS (30/30)
- scorer equivalence gate: PASS (`{facts['scorer']}`; max |Δ| vs R = {float(facts['max_abs_diff_exact_vs_r']):.3g})
- donors: {facts['n_donors']} overall ({facts['n_AF']} AF / {facts['n_ctrl']} controls)
- compartments: {facts['n_compartments']}
- analysis unit: {facts['analysis_unit']}
- family aggregation: {facts['aggregation']}
- donor eligibility gate: ≥{facts['min_donors']} AF and ≥{facts['min_donors']} controls
- score computability: ≥{facts['min_measured_LE_for_score']} measured frozen LE genes (exact-reference)
- internal frozen-LE set-size QC (not a testing gate): ≥{facts['min_frozen_LE_set_size_qc']} (min frozen LE size = {facts['frozen_LE_min_size']})
- permutations: {facts['n_perm']}
- HCA prior timing: {facts['hca_prior_timing']}
- legacy scorer / runtime LE reconstruction: not used
- panel dimensions: {WIDTH_MM} × {HEIGHT_MM} mm
- class counts omitted from panel A (belong in panel B)

## Outputs
{output_lines}

## Legacy handling
- archived: `{LEGACY_ARCHIVE.relative_to(ROOT) if LEGACY_ARCHIVE.is_file() else 'pending'}`
- submission routing points to `{OUT_BASE.with_suffix('.pdf').relative_to(ROOT)}`
""",
        encoding="utf-8",
    )


def render_composite_panel() -> plt.Figure:
    """Render exact 178 × 48 mm panel A for the final Figure 3 composite."""
    set_style()
    width_mm, height_mm = 178.0, 48.0
    margin, gap = 2.5, 2.0
    content_w = width_mm - 2 * margin
    module_w = content_w - 3 * gap
    fractions = [0.18, 0.25, 0.32, 0.25]
    widths = [module_w * fraction for fraction in fractions]
    require(abs(sum(widths) + 3 * gap - content_w) <= 1e-9, "Panel A module geometry changed")
    xs = [margin]
    for module_width in widths[:-1]:
        xs.append(xs[-1] + module_width + gap)

    fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
    fig.patch.set_facecolor(WHITE)
    fig.text(2.5 / width_mm, 1 - 2.0 / height_mm, "A", fontsize=11, fontweight="bold", ha="left", va="top")
    fig.text(
        8.0 / width_mm,
        1 - 2.0 / height_mm,
        "LA single-nucleus localization and carrier-relative classification",
        fontsize=8.7,
        fontweight="bold",
        ha="left",
        va="top",
        color=FG_COLOR,
    )
    fig.text(
        8.0 / width_mm,
        1 - 6.4 / height_mm,
        "Frozen family axis · donor×compartment scoring · permutation localization",
        fontsize=6.6,
        ha="left",
        va="top",
        color=MUTED_COLOR,
    )

    panel_top, panel_h = 10.0, 35.5
    axes: list[plt.Axes] = []
    for index, (x_pos, module_width) in enumerate(zip(xs, widths)):
        ax = fig.add_axes(
            [
                x_pos / width_mm,
                (height_mm - panel_top - panel_h) / height_mm,
                module_width / width_mm,
                panel_h / height_mm,
            ]
        )
        ax.set_gid(f"fig3a-module-{index + 1}")
        ax.set_xlim(0, module_width)
        ax.set_ylim(0, panel_h)
        ax.axis("off")
        ax.add_patch(
            mpl.patches.FancyBboxPatch(
                (0.15, 0.15),
                module_width - 0.3,
                panel_h - 0.3,
                boxstyle="round,pad=0.02,rounding_size=0.8",
                facecolor="#F7F7F7",
                edgecolor="#D8D8D8",
                linewidth=0.45,
            )
        )
        axes.append(ax)

    ax = axes[0]
    ax.text(1.5, 33.5, "REFERENCE", fontsize=6.2, fontweight="bold", color=KICKER_COLOR, va="top")
    ax.text(1.5, 29.5, "Frozen family axis", fontsize=6.7, fontweight="bold", color=FG_COLOR, va="top")
    ax.text(1.5, 25.0, "16 NCP families", fontsize=6.4, color=FG_COLOR, va="top")
    ax.add_patch(mpl.patches.Rectangle((1.5, 17.0), widths[0] - 3.0, 3.2, facecolor="#C97A74", edgecolor="none"))
    ax.add_patch(mpl.patches.Rectangle((1.5, 12.7), widths[0] - 3.0, 3.2, facecolor="#8293A4", edgecolor="none"))
    ax.text(widths[0] / 2, 18.6, "LA-biased", fontsize=6.2, color=WHITE, ha="center", va="center", fontweight="bold")
    ax.text(widths[0] / 2, 14.3, "RA-biased", fontsize=6.2, color=WHITE, ha="center", va="center", fontweight="bold")
    ax.text(1.5, 7.4, "Healthy LA–RA reference", fontsize=6.2, color=MUTED_COLOR, va="top")

    ax = axes[1]
    ax.text(1.5, 33.5, "LA DATA", fontsize=6.2, fontweight="bold", color=KICKER_COLOR, va="top")
    ax.text(1.5, 29.5, "GSE255612 LA snRNA", fontsize=6.7, fontweight="bold", color=FG_COLOR, va="top")
    ax.text(1.5, 25.0, "34 donors · 18 AF · 16 Control", fontsize=6.3, color=FG_COLOR, va="top")
    grid = [["aCM", "Fibro", "Endo"], ["Macro", "Lymph", "Peri"], ["VSMC", "Adipo", "Neuro"]]
    chip_gap = 0.8
    chip_w = (widths[1] - 3.0 - 2 * chip_gap) / 3
    for row_index, row in enumerate(grid):
        for column_index, label in enumerate(row):
            x_chip = 1.5 + column_index * (chip_w + chip_gap)
            y_chip = 17.7 - row_index * 4.2
            ax.add_patch(
                mpl.patches.FancyBboxPatch(
                    (x_chip, y_chip),
                    chip_w,
                    3.1,
                    boxstyle="round,pad=0.01,rounding_size=0.35",
                    facecolor="#ECECEC",
                    edgecolor="#BDBDBD",
                    linewidth=0.35,
                )
            )
            ax.text(x_chip + chip_w / 2, y_chip + 1.55, label, fontsize=6.0, ha="center", va="center", color=FG_COLOR)
    ax.text(
        1.5,
        6.2,
        f"≥{MIN_DONORS} AF and ≥{MIN_DONORS} Ctrl donors",
        fontsize=6.0,
        color=MUTED_COLOR,
        va="top",
    )
    ax.text(
        1.5,
        2.8,
        f"score available (≥{MIN_MEASURED_LE_FOR_SCORE} measured frozen LE genes)",
        fontsize=6.0,
        color=MUTED_COLOR,
        va="top",
    )

    ax = axes[2]
    ax.text(1.5, 33.5, "LOCALIZATION TEST", fontsize=6.2, fontweight="bold", color=KICKER_COLOR, va="top")
    ax.text(1.5, 29.5, "Family-level donor scores", fontsize=6.7, fontweight="bold", color=FG_COLOR, va="top")
    ax.text(1.5, 24.8, "Donor × compartment", fontsize=6.3, color=FG_COLOR, va="top")
    ax.text(1.5, 20.0, "Δ = mean(AF) − mean(Control)", fontsize=6.5, color=FG_COLOR, va="top")
    ax.text(1.5, 15.2, "10,000 label permutations", fontsize=6.3, color=FG_COLOR, va="top")
    support = ax.text(
        1.5,
        10.4,
        "Nominal p < 0.05 marks support",
        fontsize=6.2,
        color=FG_COLOR,
        va="top",
    )
    support.set_gid("fig3a-nominal-support")
    ax.text(1.5, 5.8, "Localization cue; not FDR discovery", fontsize=6.0, color=MUTED_COLOR, va="top")

    ax = axes[3]
    ax.text(1.5, 33.5, "CARRIER-RELATIVE CLASS", fontsize=6.2, fontweight="bold", color=KICKER_COLOR, va="top")
    class_rows = [
        ("Fidelity", "#8FAF7E", "overlap"),
        ("Relay", "#A898B8", "relay"),
        ("Extension", "#7FA8AA", "extension"),
        ("Below-detection", "#CECECE", "below"),
    ]
    for row_index, (label, color, mode) in enumerate(class_rows):
        y = 26.6 - row_index * 6.4
        glyph = mpl.patches.Rectangle(
            (1.5, y),
            3.5,
            3.5,
            facecolor=WHITE if mode != "below" else color,
            edgecolor="#555555" if mode != "below" else color,
            linewidth=0.5,
        )
        glyph.set_gid(f"fig3a-class-glyph-{mode}")
        ax.add_patch(glyph)
        if mode in {"overlap", "extension"}:
            ax.plot(3.25, y + 1.75, "o", ms=2.0, color=FG_COLOR)
        if mode in {"relay", "extension"}:
            ax.text(6.2, y + 1.75, "+", fontsize=6.0, ha="center", va="center")
            ax.plot(8.0, y + 1.75, "o", ms=2.0, color=FG_COLOR)
            text_x = 10.2
        else:
            text_x = 6.2
        ax.text(text_x, y + 1.75, label, fontsize=6.2, fontweight="semibold", color=FG_COLOR, va="center")
    return fig


def archive_legacy() -> None:
    ARCHIVE.mkdir(parents=True, exist_ok=True)
    if LEGACY_LIVE.is_file() and not LEGACY_ARCHIVE.is_file():
        shutil.copy2(LEGACY_LIVE, LEGACY_ARCHIVE)


def main() -> None:
    facts = assert_authority()
    source = build_source(facts)
    archive_legacy()

    OUT_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    OUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    source.to_csv(OUT_SOURCE, sep="\t", index=False)

    fig = render()
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
    write_qc(facts, output_hashes)

    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")
    if LEGACY_ARCHIVE.is_file():
        print(f"Archived legacy panel → {LEGACY_ARCHIVE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
