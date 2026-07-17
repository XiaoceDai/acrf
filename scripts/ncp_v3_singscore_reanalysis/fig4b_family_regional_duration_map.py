#!/usr/bin/env python3
"""Figure 4B — family-level regional duration-trend map (corrected namespace)."""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"

OUT_BASE = FIGURES / "fig4b_family_regional_duration_map"
OUT_PREVIEW = FIGURES / "fig4b_family_regional_duration_map_preview.png"
OUT_SOURCE = RESULTS / "fig4b_family_regional_duration_source.tsv"
OUT_SUMMARY = RESULTS / "fig4b_family_regional_duration_summary.tsv"
OUT_SCALE_QC = RESULTS / "fig4b_family_regional_duration_scale_qc.tsv"
OUT_QC = RESULTS / "fig4b_family_regional_duration_qc.md"

LAA_TSV = RESULTS / "gse41177_family_duration_laa_primary.tsv"
LAPV_TSV = RESULTS / "gse41177_family_duration_lapv_sensitivity.tsv"
CHRONICITY_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"
CANONICAL_TSV = ROOT / "results" / "ncp_v3_canonical_display_order.tsv"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"

ALLOWED_INPUTS = {
    "laa_primary": LAA_TSV,
    "lapv_sensitivity": LAPV_TSV,
    "chronicity_class": CHRONICITY_TSV,
    "canonical_order": CANONICAL_TSV,
    "family_labels": LABELS_TSV,
}

FORBIDDEN_INPUT_PATTERNS = (
    re.compile(r"results/fig4_v3_"),
    re.compile(r"results/fig4b_family_duration_correlations"),
    re.compile(r"fig_integrated_af_rewiring_atlas"),
    re.compile(r"fig1f_v3_HCA_family_carriers"),
    re.compile(r"acairel", re.I),
    re.compile(r"gse41177_duration_metadata"),
)

PUBLICATION_LABELS: dict[int, str] = {
    1: "Translation / RQC",
    2: "Mito Ca²⁺ transport",
    3: "Adipo–lipid axis",
    4: "Mito proteostasis / FA–BCAA",
    5: "tRNA aminoacylation",
    6: "Mito translation / import",
    7: "OXPHOS / ETC",
    8: "Complex I biogenesis",
    12: "rRNA processing",
    13: "Amino-acid metabolism",
    14: "Translation initiation",
    9: "HS-GAG / bone / GPCR",
    10: "GAG disease / TGF-β",
    11: "BMP signalling",
    15: "NTRK2 / TrkB",
    16: "Keratan sulfate",
}

EXPECTED_LA_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14]
EXPECTED_RA_ORDER = [9, 10, 11, 15, 16]
EXPECTED_FAM_ORDER = EXPECTED_LA_ORDER + EXPECTED_RA_ORDER

SCORER = "python_exact_reference"
# Minimum uniform display precision that eliminates class-boundary ambiguity
# (F11 rho=-0.19965 is Stable; 3dp shows -0.200 → reads Weak).
RHO_DISPLAY_DECIMALS = 4
VMIN, VMAX = -1.0, 1.0

MM = 1 / 25.4
# Content ends ~94.5 mm; keep a tight right margin so C can sit close in the composite.
WIDTH_MM = 98.0
HEIGHT_MM = 82.0
PREVIEW_WIDTH_MM = 178.0
MARGIN_MM = 2.5
TITLE_H_MM = 9.5
TITLE_HEADER_GAP_MM = 1.0
HEADER_H_MM = 8.5
ROW_H_MM = 3.0
FOOTER_H_MM = 15.0
N_ROWS = 16

W_DIR_MM = 3.0
W_FID_MM = 6.5
W_LBL_MM = 33.0
W_RHO_COL_MM = 13.0
W_RHO_GAP_MM = 2.0
W_REG_COL_MM = 13.5
GAP_MM = 1.0
GAP_FID_LBL_MM = 2.0
GAP_LBL_HEAT_MM = 3.0
GAP_REG_MM = 2.0

FG_COLOR = "#1A1A1A"
MUTED_COLOR = "#595959"
FAINT_COLOR = "#7A7A7A"
WHITE = "#FFFFFF"
ROW_EDGE_COLOR = "#F5F5F5"
LA_COLOR = "#9C5550"
RA_COLOR = "#4E6B89"
VIOLET = "#6B4C9A"
OCHRE = "#C8A076"

MIN_FONT_PT = 6.5
BBOX_TOL_MM = 0.55
OVERLAP_TOL_MM2 = 0.15

FORBIDDEN_RENDER_TOKENS = (
    "p <",
    "p=",
    "q <",
    "PriorityScore",
    "Path A",
    "aCAIrel",
    "✓",
    "✗",
    "★",
    "☆",
)


@dataclass(frozen=True)
class PanelGeometry:
    x_dir: float
    x_fid: float
    x_lbl: float
    x_laa: float
    x_lapv: float
    x_reg: float
    matrix_bottom: float
    matrix_top: float
    title_bottom: float


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
            "axes.linewidth": 0.6,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": WHITE,
            "savefig.facecolor": WHITE,
        }
    )


def assert_no_forbidden_reads() -> None:
    allowed = {str(path) for path in ALLOWED_INPUTS.values()}
    for line in Path(__file__).read_text(encoding="utf-8").splitlines():
        if "read_csv(" not in line:
            continue
        for pattern in FORBIDDEN_INPUT_PATTERNS:
            match = pattern.search(line)
            if match is None:
                continue
            require(
                any(path in line for path in allowed),
                f"Forbidden read_csv path reference: {match.group()}",
            )


def svg_dimensions_mm(svg_text: str) -> tuple[float, float]:
    width = re.search(r'\bwidth="([0-9.]+)(mm|pt)"', svg_text)
    height = re.search(r'\bheight="([0-9.]+)(mm|pt)"', svg_text)
    require(width is not None and height is not None, "SVG dimensions missing")
    scale = 1.0 if width.group(2) == "mm" else 25.4 / 72.0
    return float(width.group(1)) * scale, float(height.group(1)) * scale


def regional_direction(rho_laa: float, rho_lapv: float) -> str:
    if not np.isfinite(rho_lapv):
        return "NA"
    if np.sign(rho_laa) == np.sign(rho_lapv):
        return "Same"
    return "Opp."


def classify_from_abs_rho(abs_rho: float) -> str:
    if abs_rho >= 0.4:
        return "Strong_trend"
    if abs_rho >= 0.2:
        return "Weak_trend"
    return "Stable"


def displayed_rho_string(rho: float) -> str:
    return f"{rho:.{RHO_DISPLAY_DECIMALS}f}"


def assert_display_class_consistency(source: pd.DataFrame) -> None:
    """Displayed rounded |ρ| must not imply a different chronicity class than full precision."""
    boundary_strings = {
        f"{0.2:.{RHO_DISPLAY_DECIMALS}f}",
        f"{0.4:.{RHO_DISPLAY_DECIMALS}f}",
    }
    for row in source.itertuples(index=False):
        rho = float(row.rho_laa)
        abs_full = abs(rho)
        full_class = classify_from_abs_rho(abs_full)
        require(
            full_class == row.chronicity_class,
            f"F{row.family_id} chronicity_class != class from full-precision |rho_LAA|",
        )
        rounded = round(rho, RHO_DISPLAY_DECIMALS)
        display_class = classify_from_abs_rho(abs(rounded))
        disp = displayed_rho_string(rho)
        abs_disp = f"{abs_full:.{RHO_DISPLAY_DECIMALS}f}"
        require(
            display_class == full_class,
            f"F{row.family_id}: displayed {disp} implies {display_class} but full-precision class is {full_class}",
        )
        if abs_disp in boundary_strings and abs_full not in (0.2, 0.4):
            require(
                False,
                f"F{row.family_id}: displayed |rho|={abs_disp} sits on a class boundary "
                f"while full-precision |rho|={abs_full:.10g} does not",
            )


def canonical_family_order() -> list[int]:
    order = (
        pd.read_csv(CANONICAL_TSV, sep="\t")
        .sort_values("family_display_order")
        .drop_duplicates("family_id")
    )
    return order["family_id"].astype(int).tolist()


def panel_geometry() -> PanelGeometry:
    matrix_h = N_ROWS * ROW_H_MM
    require(
        abs(TITLE_H_MM + TITLE_HEADER_GAP_MM + HEADER_H_MM + matrix_h + FOOTER_H_MM - HEIGHT_MM) < 0.01,
        "Vertical geometry does not sum to panel height",
    )
    matrix_bottom = FOOTER_H_MM
    matrix_top = matrix_bottom + matrix_h
    title_bottom = matrix_top + HEADER_H_MM + TITLE_HEADER_GAP_MM

    x_dir = MARGIN_MM
    x_fid = x_dir + W_DIR_MM + GAP_MM
    x_lbl = x_fid + W_FID_MM + GAP_FID_LBL_MM
    x_laa = x_lbl + W_LBL_MM + GAP_LBL_HEAT_MM
    x_lapv = x_laa + W_RHO_COL_MM + W_RHO_GAP_MM
    x_reg = x_lapv + W_RHO_COL_MM + GAP_REG_MM
    require(x_reg + W_REG_COL_MM <= WIDTH_MM - MARGIN_MM + 0.01, "Right edge would clip")

    return PanelGeometry(
        x_dir=x_dir,
        x_fid=x_fid,
        x_lbl=x_lbl,
        x_laa=x_laa,
        x_lapv=x_lapv,
        x_reg=x_reg,
        matrix_bottom=matrix_bottom,
        matrix_top=matrix_top,
        title_bottom=title_bottom,
    )


def load_and_build_source() -> pd.DataFrame:
    for role, path in ALLOWED_INPUTS.items():
        require(path.is_file(), f"Missing allowed input ({role}): {path}")

    laa = pd.read_csv(LAA_TSV, sep="\t")
    lapv = pd.read_csv(LAPV_TSV, sep="\t")
    chronicity = pd.read_csv(CHRONICITY_TSV, sep="\t")
    canonical = pd.read_csv(CANONICAL_TSV, sep="\t")
    labels = pd.read_csv(LABELS_TSV, sep="\t")

    fam_order = canonical_family_order()
    require(len(fam_order) == 16, f"Expected 16 families, got {len(fam_order)}")
    require(fam_order == EXPECTED_FAM_ORDER, "Canonical order does not match Figure 1D")
    missing_labels = [fid for fid in fam_order if fid not in PUBLICATION_LABELS]
    require(not missing_labels, f"Missing publication labels for families: {missing_labels}")

    laa = laa.set_index("family_id")
    lapv = lapv.set_index("family_id")
    chronicity = chronicity.set_index("family_id")
    labels = labels.set_index("family_id")
    dir_lookup = dict(
        zip(
            canonical.drop_duplicates("family_id")["family_id"].astype(int),
            canonical.drop_duplicates("family_id")["canonical_direction"].astype(str),
        )
    )

    rows = []
    for rank, fid in enumerate(fam_order, start=1):
        require(fid in laa.index, f"Missing LAA row for family {fid}")
        require(fid in lapv.index, f"Missing LAPV row for family {fid}")
        require(fid in chronicity.index, f"Missing chronicity row for family {fid}")
        require(fid in labels.index, f"Missing label for family {fid}")

        rho_laa = float(laa.loc[fid, "rho_spearman"])
        rho_lapv = float(lapv.loc[fid, "rho_spearman"])
        n_laa = int(laa.loc[fid, "n_samples"])
        n_lapv = int(lapv.loc[fid, "n_samples"])

        rows.append(
            {
                "canonical_order": rank,
                "family_id": fid,
                "canonical_direction": dir_lookup[fid],
                "family_label": str(labels.loc[fid, "short_label"]),
                "publication_label": PUBLICATION_LABELS[fid],
                "rho_laa": rho_laa,
                "n_laa": n_laa,
                "rho_lapv": rho_lapv,
                "n_lapv": n_lapv,
                "regional_direction": regional_direction(rho_laa, rho_lapv),
                "chronicity_class": str(chronicity.loc[fid, "chronicity_class"]),
                "laa_analysis": str(laa.loc[fid, "analysis"]),
                "lapv_analysis": str(lapv.loc[fid, "analysis"]),
            }
        )

    source = pd.DataFrame(rows)
    assert_scientific_contract(source, chronicity)
    return source


def assert_scientific_contract(source: pd.DataFrame, chronicity: pd.DataFrame) -> None:
    require(len(source) == 16, "Source table must have 16 rows")
    require(source["family_id"].tolist() == EXPECTED_FAM_ORDER, "Source row order drifted")
    require((source["n_laa"] == 12).all(), "LAA primary n must be 12 for all families")
    require((source["n_lapv"] == 12).all(), "LAPV sensitivity n must be 12 for all families")
    require((source["laa_analysis"] == "laa_primary").all(), "LAA analysis tag must be laa_primary")
    require((source["lapv_analysis"] == "lapv_sensitivity").all(), "LAPV analysis tag must be lapv_sensitivity")

    rhos = pd.concat([source["rho_laa"], source["rho_lapv"]])
    require(rhos.notna().all(), "Missing rho values")
    require((rhos >= VMIN).all() and (rhos <= VMAX).all(), "rho values outside [-1, 1]")

    counts = chronicity["chronicity_class"].value_counts()
    require(int(counts.get("Strong_trend", 0)) == 2, "Strong_trend count != 2")
    require(int(counts.get("Weak_trend", 0)) == 6, "Weak_trend count != 6")
    require(int(counts.get("Stable", 0)) == 8, "Stable count != 8")

    f2 = source.loc[source["family_id"] == 2].iloc[0]
    require(f2["chronicity_class"] == "Strong_trend", "F2 must be Strong_trend")
    require(f2["rho_laa"] < 0, "F2 LAA rho must be negative")
    require(f2["rho_lapv"] > 0, "F2 LAPV rho must be positive")
    require(f2["regional_direction"] == "Opp.", "F2 regional direction must be Opp.")

    f9 = source.loc[source["family_id"] == 9].iloc[0]
    require(f9["chronicity_class"] == "Strong_trend", "F9 must be Strong_trend")
    require(f9["rho_laa"] < 0, "F9 LAA rho must be negative")
    require(f9["rho_lapv"] < 0, "F9 LAPV rho must be negative")
    require(f9["regional_direction"] == "Same", "F9 regional direction must be Same")

    f10 = source.loc[source["family_id"] == 10].iloc[0]
    require(f10["chronicity_class"] == "Weak_trend", "F10 must be Weak_trend")

    f15 = source.loc[source["family_id"] == 15].iloc[0]
    require(f15["chronicity_class"] == "Weak_trend", "F15 must be Weak_trend")

    assert_display_class_consistency(source)


def text_bbox_figure_mm(artist: mpl.text.Text, fig: plt.Figure) -> tuple[float, float, float, float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    return (
        float(bbox.x0) * WIDTH_MM,
        float(bbox.y0) * HEIGHT_MM,
        float(bbox.x1) * WIDTH_MM,
        float(bbox.y1) * HEIGHT_MM,
    )


def bbox_intersection(
    a: tuple[float, float, float, float],
    b: tuple[float, float, float, float],
) -> float:
    x0 = max(a[0], b[0])
    y0 = max(a[1], b[1])
    x1 = min(a[2], b[2])
    y1 = min(a[3], b[3])
    if x1 <= x0 or y1 <= y0:
        return 0.0
    return (x1 - x0) * (y1 - y0)


def collect_visible_text(fig: plt.Figure) -> list[tuple[str, mpl.text.Text, tuple[float, float, float, float], str]]:
    fig.canvas.draw()
    geom = panel_geometry()
    entries: list[tuple[str, mpl.text.Text, tuple[float, float, float, float], str]] = []
    for ax in fig.axes:
        for child in ax.get_children():
            if not isinstance(child, mpl.text.Text):
                continue
            label = child.get_text().strip()
            if not label:
                continue
            bbox = text_bbox_figure_mm(child, fig)
            y_mid = 0.5 * (bbox[1] + bbox[3])
            if y_mid >= geom.title_bottom:
                zone = "title"
            elif y_mid >= geom.matrix_top:
                zone = "header"
            elif y_mid >= FOOTER_H_MM:
                zone = "matrix"
            else:
                zone = "footer"
            entries.append((label, child, bbox, zone))
    return entries


def validate_layout(fig: plt.Figure) -> dict[str, float | int]:
    entries = collect_visible_text(fig)
    require(entries, "layout QC found no visible text")

    min_font = min(child.get_fontsize() for _, child, _, _ in entries)
    violations = 0
    overlaps = 0
    skip_zone_pairs = {("title", "header"), ("title", "title")}

    for label, _, bbox, _ in entries:
        if bbox[0] < -BBOX_TOL_MM:
            violations += 1
            require(False, f"text clipped left: {label!r}")
        if bbox[2] > WIDTH_MM + BBOX_TOL_MM:
            violations += 1
            require(False, f"text clipped right: {label!r}")
        if bbox[1] < -BBOX_TOL_MM:
            violations += 1
            require(False, f"text clipped bottom: {label!r}")
        if bbox[3] > HEIGHT_MM + BBOX_TOL_MM:
            violations += 1
            require(False, f"text clipped top: {label!r}")

    for index, (label_a, child_a, bbox_a, zone_a) in enumerate(entries):
        fontsize = child_a.get_fontsize()
        require(fontsize >= MIN_FONT_PT - 0.05, f"font below minimum: {label_a!r} ({fontsize:.2f} pt)")
        for label_b, child_b, bbox_b, zone_b in entries[index + 1 :]:
            if (zone_a, zone_b) in skip_zone_pairs or (zone_b, zone_a) in skip_zone_pairs:
                continue
            if zone_a == "header" and zone_b == "header" and child_a.axes is child_b.axes:
                continue
            if zone_a == "matrix" and zone_b == "matrix" and child_a.axes is not child_b.axes:
                overlap = bbox_intersection(bbox_a, bbox_b)
                if overlap > OVERLAP_TOL_MM2:
                    overlaps += 1
                    require(False, f"overlapping matrix text: {label_a!r} vs {label_b!r} ({overlap:.2f} mm²)")
                continue
            overlap = bbox_intersection(bbox_a, bbox_b)
            if overlap > OVERLAP_TOL_MM2:
                overlaps += 1
                require(False, f"overlapping text: {label_a!r} vs {label_b!r} ({overlap:.2f} mm²)")

    primary_boxes = [bbox for label, _, bbox, zone in entries if label == "Primary" and zone == "header"]
    sensitivity_boxes = [
        bbox
        for label, _, bbox, zone in entries
        if label in {"Regional", "sensitivity", "Regional\nsensitivity", "Regional\nsign"} and zone == "header"
    ]
    require(primary_boxes and sensitivity_boxes, "header labels missing for overlap QC")
    for sens_bbox in sensitivity_boxes:
        overlap = bbox_intersection(primary_boxes[0], sens_bbox)
        require(
            overlap <= OVERLAP_TOL_MM2,
            f"Primary and Regional sensitivity headers overlap ({overlap:.2f} mm²)",
        )

    cbar_title_left = [bbox for label, _, bbox, zone in entries if "Spearman" in label and "duration" in label]
    cbar_title_right = [bbox for label, _, bbox, zone in entries if "Global scale" in label]
    if cbar_title_left and cbar_title_right:
        require(
            bbox_intersection(cbar_title_left[0], cbar_title_right[0]) <= OVERLAP_TOL_MM2,
            "colorbar title lines overlap",
        )

    return {
        "n_text_objects": len(entries),
        "min_font_pt": float(min_font),
        "n_bbox_violations": violations,
        "n_overlaps": overlaps,
    }


def render(source: pd.DataFrame) -> tuple[plt.Figure, LinearSegmentedColormap, TwoSlopeNorm]:
    geom = panel_geometry()
    matrix_h = N_ROWS * ROW_H_MM

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM))
    fig.patch.set_facecolor(WHITE)

    def mm_axes(x_mm: float, y_mm: float, w_mm: float, h_mm: float) -> plt.Axes:
        return fig.add_axes([x_mm / WIDTH_MM, y_mm / HEIGHT_MM, w_mm / WIDTH_MM, h_mm / HEIGHT_MM])

    title_ax = mm_axes(MARGIN_MM, geom.title_bottom, WIDTH_MM - 2 * MARGIN_MM, TITLE_H_MM)
    title_ax.axis("off")
    title_ax.set_xlim(0, WIDTH_MM - 2 * MARGIN_MM)
    title_ax.set_ylim(0, TITLE_H_MM)
    title_ax.text(0, 8.8, "B", fontsize=10.5, fontweight="bold", va="top", color=FG_COLOR)
    title_ax.text(
        6.0,
        8.8,
        "Family-level duration trends across atrial regions",
        fontsize=8.5,
        fontweight="bold",
        va="top",
        color=FG_COLOR,
    )
    title_ax.text(
        6.0,
        5.0,
        "Primary LAA associations and LAPV regional sensitivity",
        fontsize=MIN_FONT_PT,
        va="top",
        color=MUTED_COLOR,
    )
    title_ax.text(
        6.0,
        1.8,
        "16 frozen NCP families · LAA primary: 12 AF samples · Spearman ρ",
        fontsize=MIN_FONT_PT,
        va="top",
        color=FAINT_COLOR,
    )

    n_la = int((source["canonical_direction"] == "LA-biased").sum())
    n_laa = int(source["n_laa"].iloc[0])
    n_lapv = int(source["n_lapv"].iloc[0])

    ax_dir = mm_axes(geom.x_dir, geom.matrix_bottom, W_DIR_MM, matrix_h)
    ax_dir.set_xlim(0, 1)
    ax_dir.set_ylim(N_ROWS, 0)
    ax_dir.axis("off")
    for row_index, row in enumerate(source.itertuples(index=False)):
        color = LA_COLOR if row.canonical_direction == "LA-biased" else RA_COLOR
        ax_dir.add_patch(
            mpatches.Rectangle((0, row_index), 1, 1, facecolor=color, edgecolor=ROW_EDGE_COLOR, linewidth=0.12)
        )
    if n_la > 0:
        ax_dir.text(
            0.5,
            n_la / 2,
            "LA",
            fontsize=MIN_FONT_PT,
            color=WHITE,
            fontweight="bold",
            ha="center",
            va="center",
            rotation=90,
            clip_on=True,
        )
    ax_dir.text(
        0.5,
        n_la + (N_ROWS - n_la) / 2,
        "RA",
        fontsize=MIN_FONT_PT,
        color=WHITE,
        fontweight="bold",
        ha="center",
        va="center",
        rotation=90,
        clip_on=True,
    )
    ax_dir.plot([0, 1], [n_la, n_la], color=WHITE, linewidth=0.6, clip_on=False, zorder=5)

    ax_fid = mm_axes(geom.x_fid, geom.matrix_bottom, W_FID_MM, matrix_h)
    ax_fid.set_xlim(0, 1)
    ax_fid.set_ylim(N_ROWS, 0)
    ax_fid.axis("off")

    ax_lbl = mm_axes(geom.x_lbl, geom.matrix_bottom, W_LBL_MM, matrix_h)
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(N_ROWS, 0)
    ax_lbl.axis("off")

    heat = source[["rho_laa", "rho_lapv"]].to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list("rho_vwo", [(0.0, VIOLET), (0.5, WHITE), (1.0, OCHRE)])
    norm = TwoSlopeNorm(vmin=VMIN, vcenter=0.0, vmax=VMAX)

    for row_index, row in enumerate(source.itertuples(index=False)):
        ax_fid.text(
            0.92,
            row_index + 0.5,
            f"F{row.family_id}",
            fontsize=MIN_FONT_PT,
            fontweight="bold",
            ha="right",
            va="center",
            color=FG_COLOR,
            clip_on=True,
        )
        ax_lbl.text(
            0.04,
            row_index + 0.5,
            row.publication_label,
            fontsize=MIN_FONT_PT,
            ha="left",
            va="center",
            color=FG_COLOR,
            clip_on=True,
        )

    def draw_rho_column(ax: plt.Axes, column_index: int, x_mm: float) -> None:
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(N_ROWS - 0.5, -0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
        for row_index in range(N_ROWS):
            val = heat[row_index, column_index]
            if not np.isfinite(val):
                continue
            face = cmap(norm(val))
            ax.add_patch(
                mpatches.Rectangle(
                    (-0.5, row_index - 0.5),
                    1,
                    1,
                    facecolor=face,
                    edgecolor=ROW_EDGE_COLOR,
                    linewidth=0.12,
                    zorder=1,
                )
            )
            txt_color = WHITE if abs(val) >= 0.35 else FG_COLOR
            ax.text(
                0,
                row_index,
                f"{val:.{RHO_DISPLAY_DECIMALS}f}",
                ha="center",
                va="center",
                fontsize=MIN_FONT_PT,
                color=txt_color,
                fontweight="normal",
                zorder=2,
            )

    ax_laa = mm_axes(geom.x_laa, geom.matrix_bottom, W_RHO_COL_MM, matrix_h)
    ax_lapv = mm_axes(geom.x_lapv, geom.matrix_bottom, W_RHO_COL_MM, matrix_h)
    draw_rho_column(ax_laa, 0, geom.x_laa)
    draw_rho_column(ax_lapv, 1, geom.x_lapv)

    ax_reg = mm_axes(geom.x_reg, geom.matrix_bottom, W_REG_COL_MM, matrix_h)
    ax_reg.set_xlim(0, 1)
    ax_reg.set_ylim(N_ROWS, 0)
    ax_reg.axis("off")
    for row_index, row in enumerate(source.itertuples(index=False)):
        ax_reg.text(
            0.5,
            row_index + 0.5,
            row.regional_direction,
            ha="center",
            va="center",
            fontsize=MIN_FONT_PT,
            color=FG_COLOR,
        )

    def draw_header(
        ax: plt.Axes,
        lines: tuple[tuple[str, dict], ...],
        y_positions: tuple[float, ...] | None = None,
    ) -> None:
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, HEADER_H_MM)
        default_y = (7.0, 4.0, 1.2)
        positions = y_positions if y_positions is not None else default_y[: len(lines)]
        for (text, style), y in zip(lines, positions):
            style_copy = dict(style)
            fontsize = style_copy.pop("fontsize", MIN_FONT_PT)
            linespacing = style_copy.pop("linespacing", None)
            va = "center" if "\n" in text else "baseline"
            ax.text(
                0.5,
                y,
                text,
                ha="center",
                va=va,
                fontsize=fontsize,
                linespacing=linespacing if linespacing is not None else 1.0,
                **style_copy,
            )

    hdr_laa = mm_axes(geom.x_laa, geom.matrix_top, W_RHO_COL_MM, HEADER_H_MM)
    draw_header(
        hdr_laa,
        (
            ("LAA", {"fontweight": "bold", "fontsize": 7.0, "color": FG_COLOR}),
            ("Primary", {"fontweight": "semibold", "color": MUTED_COLOR}),
            (f"n = {n_laa}", {"color": FAINT_COLOR}),
        ),
    )

    hdr_lapv = mm_axes(geom.x_lapv, geom.matrix_top, W_RHO_COL_MM, HEADER_H_MM)
    draw_header(
        hdr_lapv,
        (
            ("LAPV", {"fontweight": "semibold", "color": MUTED_COLOR}),
            ("Regional\nsensitivity", {"fontweight": "normal", "color": MUTED_COLOR, "linespacing": 1.1}),
            (f"n = {n_lapv}", {"color": FAINT_COLOR}),
        ),
        y_positions=(7.0, 4.0, 1.2),
    )

    hdr_reg = mm_axes(geom.x_reg, geom.matrix_top, W_REG_COL_MM, HEADER_H_MM)
    hdr_reg.axis("off")
    hdr_reg.set_xlim(0, 1)
    hdr_reg.set_ylim(0, HEADER_H_MM)
    hdr_reg.text(
        0.5,
        4.5,
        "Regional\nsign",
        ha="center",
        va="center",
        fontsize=MIN_FONT_PT,
        color=FG_COLOR,
        linespacing=1.1,
    )

    footer_ax = mm_axes(MARGIN_MM, 0.0, WIDTH_MM - 2 * MARGIN_MM, FOOTER_H_MM)
    footer_ax.axis("off")
    footer_ax.set_xlim(0, WIDTH_MM - 2 * MARGIN_MM)
    footer_ax.set_ylim(0, FOOTER_H_MM)
    footer_w = WIDTH_MM - 2 * MARGIN_MM
    cbar_width = 2 * W_RHO_COL_MM + W_RHO_GAP_MM
    cbar_left = geom.x_laa - MARGIN_MM
    cbar_right = cbar_left + cbar_width

    # Compact footer: notes under labels (left); colorbar under LAA+LAPV (right)
    note_lines = (
        (13.4, "No significance encoding."),
        (10.8, "Regional sign: Same = shared"),
        (8.2, "LAA/LAPV sign; Opp. = opposite."),
        (5.6, "Chronicity classes use"),
        (3.0, "unrounded LAA Spearman ρ."),
    )
    for y_mm, text in note_lines:
        footer_ax.text(
            0.0,
            y_mm,
            text,
            ha="left",
            va="top",
            fontsize=MIN_FONT_PT,
            color=FAINT_COLOR,
            clip_on=False,
        )

    cbar_y = 11.8
    cbar_h = 1.15
    cbar_mm = mm_axes(geom.x_laa, cbar_y, cbar_width, cbar_h)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_mm, orientation="horizontal")
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cbar.ax.tick_params(labelsize=MIN_FONT_PT, length=1.8, width=0.35, pad=0.8)
    cbar.outline.set_linewidth(0.35)

    cbar_center = cbar_left + cbar_width / 2
    assoc_y = cbar_y - 3.4
    footer_ax.text(
        cbar_left,
        assoc_y,
        "negative",
        ha="left",
        va="top",
        fontsize=MIN_FONT_PT,
        color=FAINT_COLOR,
    )
    footer_ax.text(
        cbar_right,
        assoc_y,
        "positive",
        ha="right",
        va="top",
        fontsize=MIN_FONT_PT,
        color=FAINT_COLOR,
    )
    footer_ax.text(
        cbar_center,
        assoc_y - 2.6,
        "Spearman ρ with AF duration (−1 to +1)",
        ha="center",
        va="top",
        fontsize=MIN_FONT_PT,
        color=FG_COLOR,
    )

    return fig, cmap, norm


def save_preview(fig: plt.Figure, path: Path) -> None:
    scale = WIDTH_MM / PREVIEW_WIDTH_MM
    orig_size = fig.get_size_inches().copy()
    fig.set_size_inches(PREVIEW_WIDTH_MM * MM, HEIGHT_MM * MM)
    for ax in fig.axes:
        pos = ax.get_position()
        ax.set_position([pos.x0 * scale, pos.y0, pos.width * scale, pos.height])
    fig.savefig(path, dpi=150, bbox_inches=None, pad_inches=0)
    fig.set_size_inches(orig_size[0], orig_size[1])
    for ax in fig.axes:
        pos = ax.get_position()
        ax.set_position([pos.x0 / scale, pos.y0, pos.width / scale, pos.height])


def validate_outputs(pdf: Path, svg: Path, png: Path) -> dict[str, object]:
    require(pdf.is_file() and pdf.stat().st_size > 0, "PDF export missing")
    require(svg.is_file() and svg.stat().st_size > 0, "SVG export missing")
    require(png.is_file() and png.stat().st_size > 0, "PNG export missing")

    svg_text = svg.read_text(encoding="utf-8")
    svg_width_mm, svg_height_mm = svg_dimensions_mm(svg_text)
    require(abs(svg_width_mm - WIDTH_MM) < 0.5, f"SVG width {svg_width_mm:.2f} mm not ~{WIDTH_MM:.0f} mm")
    require(abs(svg_height_mm - HEIGHT_MM) < 0.5, f"SVG height {svg_height_mm:.2f} mm not ~{HEIGHT_MM:.0f} mm")

    text_chunks = re.findall(r"<text[^>]*>([^<]*)</text>", svg_text)
    combined = "\n".join(text_chunks)
    for token in FORBIDDEN_RENDER_TOKENS:
        require(token not in combined, f"Forbidden token rendered: {token}")
    star_texts = [text for text in text_chunks if "*" in text]
    require(not star_texts, f"Star marker found in SVG text: {star_texts}")

    for phrase in ("Ca²⁺", "FA–BCAA", "TGF-β"):
        require(phrase in combined, f"Publication label glyph missing in SVG: {phrase}")
    # F11 is the known near-boundary Stable case (rho≈−0.19965 → 4dp −0.1996)
    require("-0.1996" in combined, "F11 4-decimal LAA display (-0.1996) missing from SVG")
    require(not re.search(r"(?<![\d.])-0\.200(?![\d])", combined), "Ambiguous −0.200 boundary display still present")

    return {"svg_width_mm": svg_width_mm, "svg_height_mm": svg_height_mm}


def write_artifacts(
    source: pd.DataFrame,
    output_hashes: dict[str, str],
    validation: dict[str, object],
    layout_stats: dict[str, float | int],
) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    export_cols = [
        "canonical_order",
        "family_id",
        "canonical_direction",
        "family_label",
        "publication_label",
        "rho_laa",
        "n_laa",
        "rho_lapv",
        "n_lapv",
        "regional_direction",
        "chronicity_class",
    ]
    source[export_cols].to_csv(OUT_SOURCE, sep="\t", index=False, float_format="%.10g")

    summary = pd.DataFrame(
        [
            {"metric": "n_families", "value": 16},
            {"metric": "canonical_order", "value": ",".join(str(fid) for fid in EXPECTED_FAM_ORDER)},
            {"metric": "order_provenance", "value": "ncp_v3_canonical_display_order family_display_order"},
            {"metric": "laa_primary_n", "value": int(source["n_laa"].iloc[0])},
            {"metric": "lapv_sensitivity_n", "value": int(source["n_lapv"].iloc[0])},
            {"metric": "scorer", "value": SCORER},
            {"metric": "rho_display_decimals", "value": RHO_DISPLAY_DECIMALS},
            {"metric": "rho_source_laa", "value": str(LAA_TSV.relative_to(ROOT))},
            {"metric": "rho_source_lapv", "value": str(LAPV_TSV.relative_to(ROOT))},
            {"metric": "n_strong_trend", "value": int((source["chronicity_class"] == "Strong_trend").sum())},
            {"metric": "n_weak_trend", "value": int((source["chronicity_class"] == "Weak_trend").sum())},
            {"metric": "n_stable", "value": int((source["chronicity_class"] == "Stable").sum())},
            {"metric": "n_same_direction", "value": int((source["regional_direction"] == "Same").sum())},
            {"metric": "n_opposite_direction", "value": int((source["regional_direction"] == "Opp.").sum())},
            {"metric": "legacy_tables_read", "value": False},
            {"metric": "acairel_tables_read", "value": False},
        ]
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    scale_qc = pd.DataFrame(
        [
            {"metric": "vmin", "value": VMIN},
            {"metric": "vcenter", "value": 0.0},
            {"metric": "vmax", "value": VMAX},
            {"metric": "row_standardized", "value": False},
            {"metric": "clipping", "value": False},
            {"metric": "significance_encoding", "value": False},
            {"metric": "rho_display_decimals", "value": RHO_DISPLAY_DECIMALS},
            {"metric": "rho_min", "value": float(pd.concat([source["rho_laa"], source["rho_lapv"]]).min())},
            {"metric": "rho_max", "value": float(pd.concat([source["rho_laa"], source["rho_lapv"]]).max())},
            {"metric": "panel_width_mm", "value": WIDTH_MM},
            {"metric": "panel_height_mm", "value": HEIGHT_MM},
            {"metric": "preview_width_mm", "value": PREVIEW_WIDTH_MM},
        ]
    )
    scale_qc.to_csv(OUT_SCALE_QC, sep="\t", index=False)

    input_lines = "\n".join(
        f"- `{path.relative_to(ROOT)}` — `{sha256(path)}`"
        for path in ALLOWED_INPUTS.values()
    )
    output_lines = "\n".join(f"- `{path}` — `{digest}`" for path, digest in sorted(output_hashes.items()))

    OUT_QC.write_text(
        f"""# Figure 4B family regional duration map QC

Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}

## Gate: **PASS**

## Scientific contract
- 16 frozen NCP families in Figure 1D canonical order
- LAA primary n = 12 AF samples for all families
- LAPV regional sensitivity n = 12 for all families
- Spearman ρ displayed from corrected duration tables only (not chronicity ρ columns)
- Display precision: {RHO_DISPLAY_DECIMALS} decimals (render-only; source rho unrounded)
- Display-class consistency: rounded display |ρ_LAA| never implies a different class than full precision
- Scorer lineage: `{SCORER}` (upstream corrected pipeline; not recomputed here)
- SR samples excluded via upstream `laa_primary` / `lapv_sensitivity` cohort tags
- No legacy Figure 4 tables read
- No aCAIrel tables read
- No PriorityScore ordering
- No p/q/stars/significance colors rendered

## Class counts (assertion-only from chronicity table)
- Strong_trend: 2
- Weak_trend: 6
- Stable: 8

## Spot checks
- F2: LAA negative Strong_trend; LAPV positive; regional Opp.
- F9: LAA negative Strong_trend; LAPV negative; regional Same
- F10: Weak_trend
- F15: Weak_trend

## Layout QC
- visible text objects: {layout_stats["n_text_objects"]}
- minimum font size: {layout_stats["min_font_pt"]:.1f} pt
- bbox violations: {layout_stats["n_bbox_violations"]}
- text overlaps: {layout_stats["n_overlaps"]}
- Primary vs Regional sensitivity headers: separated
- colorbar aligned under LAA+LAPV columns only

## Panel geometry
- Target: {WIDTH_MM:.1f} mm × {HEIGHT_MM:.1f} mm
- SVG export: {validation["svg_width_mm"]:.1f} mm × {validation["svg_height_mm"]:.1f} mm
- Preview width: {PREVIEW_WIDTH_MM:.1f} mm (B|C composite slot)

## Allowed inputs
{input_lines}

## Outputs
{output_lines}
""",
        encoding="utf-8",
    )


def main() -> None:
    set_style()
    assert_no_forbidden_reads()

    source = load_and_build_source()

    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, _, _ = render(source)
    layout_stats = validate_layout(fig)

    pdf = OUT_BASE.with_suffix(".pdf")
    svg = OUT_BASE.with_suffix(".svg")
    png = OUT_BASE.with_suffix(".png")

    fig.savefig(pdf, bbox_inches=None, pad_inches=0)
    fig.savefig(svg, bbox_inches=None, pad_inches=0)
    fig.savefig(png, dpi=600, bbox_inches=None, pad_inches=0)
    save_preview(fig, OUT_PREVIEW)
    plt.close(fig)

    validation = validate_outputs(pdf, svg, png)
    output_hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (pdf, svg, png, OUT_PREVIEW)
    }
    write_artifacts(source, output_hashes, validation, layout_stats)

    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_SCALE_QC.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
