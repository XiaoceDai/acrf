#!/usr/bin/env python3
"""Figure 5A — cross-layer triangulation contract on the frozen NCP topology.

Writes authoritative manifest + canonical triangulation table, then renders
panel A only. Does not modify Figure 5B–D or formal submission routing.
"""
from __future__ import annotations

import hashlib
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
RESULTS_ROOT = ROOT / "results"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"

OUT_BASE = FIGURES / "fig5a_cross_layer_triangulation_contract"
OUT_PREVIEW = FIGURES / "fig5a_cross_layer_triangulation_contract_preview.png"
OUT_MANIFEST = RESULTS / "figure5_authoritative_input_manifest.tsv"
OUT_CANONICAL = RESULTS / "figure5_triangulation_canonical.tsv"
OUT_SOURCE = RESULTS / "fig5a_triangulation_contract_source.tsv"
OUT_QC = RESULTS / "fig5a_triangulation_contract_qc.md"

CANONICAL_TSV = RESULTS_ROOT / "ncp_v3_canonical_display_order.tsv"
LABELS_TSV = RESULTS_ROOT / "ncp_v3_family_labels.tsv"
FROZEN_META_TSV = RESULTS / "frozen_family_metadata.tsv"
HCA_TSV = RESULTS_ROOT / "fig1f_v3_HCA_family_carriers.tsv"
BULK_FE_TSV = RESULTS / "fig2e_cross_cohort_consistency.tsv"
BULK_SUMMARY_TSV = RESULTS / "fig2d_family_effects_consistency_summary.tsv"
SN_TSV = RESULTS_ROOT / "fig3_v3_carrier_fidelity_class.tsv"
DURATION_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"
GEOMX_ROI_TSV = ROOT / "Project_P5" / "GSE261363" / "roi_scores_v3.tsv"
GEOMX_CORR_TSV = RESULTS_ROOT / "fig5d_geomx_roi_correlations.tsv"

FORBIDDEN_READ_PATTERNS = (
    re.compile(r"results/fig5b_v3_triangulation\.tsv"),
    re.compile(r"fig2e_replication_lollipop_v3"),
)

EXPECTED_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]
EXPECTED_PATH_A = {10, 15}
EXPECTED_COUNTS = {"triple": 2, "double": 1, "single": 10, "silent": 3}

MM = 1 / 25.4
WIDTH_MM = 178.0
HEIGHT_MM = 56.0
MARGIN_MM = 3.0
GAP_MM = 2.0
COL_FRAC = [0.19, 0.17, 0.19, 0.17, 0.28]
TITLE_H_MM = 11.0
GEOMX_H_MM = 6.5
STRIP_GAP_MM = 1.5
STRIP_H_MM = HEIGHT_MM - MARGIN_MM - TITLE_H_MM - STRIP_GAP_MM - GEOMX_H_MM - MARGIN_MM

BODY_PT = 6.2
FOOTNOTE_PT = 5.8
MIN_FONT_PT = 5.8
PANEL_LETTER_PT = 10.8
TITLE_TEXT_PT = 8.6
SUBTITLE_PT = 6.4
KICKER_PT = 6.0
GEOMX_PT = 5.8

_STANDALONE = {
    "WIDTH_MM": 178.0,
    "HEIGHT_MM": 56.0,
    "MARGIN_MM": 3.0,
    "TITLE_H_MM": 11.0,
    "GEOMX_H_MM": 6.5,
    "STRIP_GAP_MM": 1.5,
    "PANEL_LETTER_PT": 10.8,
    "TITLE_TEXT_PT": 8.6,
    "SUBTITLE_PT": 6.4,
    "BODY_PT": 6.2,
    "FOOTNOTE_PT": 5.8,
    "MIN_FONT_PT": 5.8,
    "KICKER_PT": 6.0,
    "GEOMX_PT": 5.8,
    "LINE_H_MM": 2.35,
}
_COMPOSITE = {
    "WIDTH_MM": 178.0,
    "HEIGHT_MM": 48.0,
    "MARGIN_MM": 2.5,
    "TITLE_H_MM": 8.5,  # free strip height for one-step type bump
    "GEOMX_H_MM": 5.2,
    "STRIP_GAP_MM": 1.0,
    # One-step type bump for composite readability
    "PANEL_LETTER_PT": 11.0,
    "TITLE_TEXT_PT": 9.0,
    "SUBTITLE_PT": 6.8,
    "BODY_PT": 6.7,
    "FOOTNOTE_PT": 6.2,
    "MIN_FONT_PT": 5.8,
    "KICKER_PT": 6.4,
    "GEOMX_PT": 6.2,
    "LINE_H_MM": 2.35,
}


def apply_layout(mode: str = "standalone") -> None:
    """Switch between standalone and Figure 5 composite canvas geometry."""
    global WIDTH_MM, HEIGHT_MM, MARGIN_MM, TITLE_H_MM, GEOMX_H_MM, STRIP_GAP_MM, STRIP_H_MM
    global PANEL_LETTER_PT, TITLE_TEXT_PT, SUBTITLE_PT, BODY_PT, FOOTNOTE_PT, MIN_FONT_PT
    global KICKER_PT, GEOMX_PT, LINE_H_MM
    cfg = _COMPOSITE if mode == "composite" else _STANDALONE
    WIDTH_MM = cfg["WIDTH_MM"]
    HEIGHT_MM = cfg["HEIGHT_MM"]
    MARGIN_MM = cfg["MARGIN_MM"]
    TITLE_H_MM = cfg["TITLE_H_MM"]
    GEOMX_H_MM = cfg["GEOMX_H_MM"]
    STRIP_GAP_MM = cfg["STRIP_GAP_MM"]
    PANEL_LETTER_PT = cfg["PANEL_LETTER_PT"]
    TITLE_TEXT_PT = cfg["TITLE_TEXT_PT"]
    SUBTITLE_PT = cfg["SUBTITLE_PT"]
    BODY_PT = cfg["BODY_PT"]
    FOOTNOTE_PT = cfg["FOOTNOTE_PT"]
    MIN_FONT_PT = cfg["MIN_FONT_PT"]
    KICKER_PT = cfg["KICKER_PT"]
    GEOMX_PT = cfg["GEOMX_PT"]
    LINE_H_MM = cfg["LINE_H_MM"]
    STRIP_H_MM = HEIGHT_MM - MARGIN_MM - TITLE_H_MM - STRIP_GAP_MM - GEOMX_H_MM - MARGIN_MM


LINE_H_MM = 2.35

FG_COLOR = "#1A1A1A"
MUTED_COLOR = "#595959"
FAINT_COLOR = "#7A7A7A"
KICKER_COLOR = "#8A8A8A"
MID_GREY = "#9A9A9A"
STRIP_BG = "#F7F7F7"
WHITE = "#FFFFFF"
GATE_EDGE = "#8A8A8A"
MODULE_EDGE = "#EBEBEB"
TRIAGE_EDGE = "#C8C8C8"
GEOMX_EDGE = "#A0A0A0"
BBOX_TOL_MM = 0.55
OVERLAP_TOL_MM2 = 0.15
PATH_A_GATE_PHRASES = (
    "High consistency",
    "non-Below-detection",
    "non-Stable",
)
MINI_GATE_LABELS = (
    "Bulk\nHigh",
    "snRNA\nnon-Below",
    "Duration\nnon-Stable",
)
RHO_LAA_MATH = r"$|\rho_{\mathrm{LAA}}|$"

INPUT_ROLES = {
    "canonical_order": CANONICAL_TSV,
    "family_labels": LABELS_TSV,
    "frozen_family_metadata": FROZEN_META_TSV,
    "hca_carrier": HCA_TSV,
    "bulk_consistency_ivw": BULK_FE_TSV,
    "bulk_consistency_summary": BULK_SUMMARY_TSV,
    "snrna_fidelity_class": SN_TSV,
    "duration_chronicity_class": DURATION_TSV,
    "geomx_roi_scores_context": GEOMX_ROI_TSV,
    "geomx_correlations_context": GEOMX_CORR_TSV,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def title_consistency(value: str) -> str:
    text = str(value).strip()
    low = text.lower()
    if low == "high":
        return "High"
    if low == "moderate":
        return "Moderate"
    if low == "low":
        return "Low"
    return text


def classify_duration(abs_rho: float) -> str:
    if abs_rho >= 0.4:
        return "Strong_trend"
    if abs_rho >= 0.2:
        return "Weak_trend"
    return "Stable"


def support_tier(count: int) -> str:
    return {3: "triple", 2: "double", 1: "single", 0: "silent"}[count]


def assert_no_forbidden_reads() -> None:
    for line in Path(__file__).read_text(encoding="utf-8").splitlines():
        if "read_csv(" not in line:
            continue
        if line.lstrip().startswith("#"):
            continue
        for pattern in FORBIDDEN_READ_PATTERNS:
            require(pattern.search(line) is None, f"Forbidden read: {line.strip()}")


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": BODY_PT,
            "mathtext.default": "regular",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": WHITE,
            "savefig.facecolor": WHITE,
        }
    )


def write_manifest() -> pd.DataFrame:
    notes = {
        "canonical_order": "Figure_1D_frozen_order",
        "family_labels": "publication_short_labels",
        "frozen_family_metadata": "direction_and_labels",
        "hca_carrier": "not_a_gate",
        "bulk_consistency_ivw": "fixed-effect_IVW; dominant_state context only",
        "bulk_consistency_summary": "High/Moderate/Low for bulk_gate",
        "snrna_fidelity_class": "non-Below gate",
        "duration_chronicity_class": "non-Stable gate; unrounded |rho_LAA|",
        "geomx_roi_scores_context": "24_eligible_ROIs; not_in_support_count",
        "geomx_correlations_context": "shown_in_panel_D; not_in_Path_A",
    }
    used_for = {
        "canonical_order": "frozen_topology",
        "family_labels": "frozen_topology",
        "frozen_family_metadata": "frozen_topology",
        "hca_carrier": "context_only_carrier_prior",
        "bulk_consistency_ivw": "bulk_consistency_and_dominant_state_context",
        "bulk_consistency_summary": "bulk_consistency_gate",
        "snrna_fidelity_class": "snRNA_class_gate",
        "duration_chronicity_class": "duration_class_gate",
        "geomx_roi_scores_context": "context_only_spatial_corroboration",
        "geomx_correlations_context": "context_only_spatial_corroboration",
    }
    gate_roles = {
        "bulk_consistency_ivw",
        "bulk_consistency_summary",
        "snrna_fidelity_class",
        "duration_chronicity_class",
    }
    rows = []
    for role, path in INPUT_ROLES.items():
        require(path.is_file(), f"Missing authoritative input ({role}): {path}")
        rows.append(
            {
                "artifact_role": role,
                "artifact_path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "used_for": used_for[role],
                "enters_Path_A_gate": role in gate_roles,
                "notes": notes[role],
            }
        )
    manifest = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(OUT_MANIFEST, sep="\t", index=False)
    return manifest


def build_canonical_triangulation() -> pd.DataFrame:
    fams = (
        pd.read_csv(CANONICAL_TSV, sep="\t")
        .sort_values("family_display_order")
        .drop_duplicates("family_id")["family_id"]
        .astype(int)
        .tolist()
    )
    require(fams == EXPECTED_ORDER, f"Canonical order drifted: {fams}")

    labels = pd.read_csv(LABELS_TSV, sep="\t").set_index("family_id")
    frozen = pd.read_csv(FROZEN_META_TSV, sep="\t").set_index("family_id")
    hca = pd.read_csv(HCA_TSV, sep="\t").set_index("family_id")
    bulk_fe = pd.read_csv(BULK_FE_TSV, sep="\t").set_index("family_id")
    bulk_sum = pd.read_csv(BULK_SUMMARY_TSV, sep="\t").set_index("family_id")
    sn = pd.read_csv(SN_TSV, sep="\t").set_index("family_id")
    duration = pd.read_csv(DURATION_TSV, sep="\t").set_index("family_id")

    require(
        set(bulk_fe["pooling_model"].astype(str)) == {"fixed-effect_inverse_variance"},
        "Bulk pooling model must be fixed-effect_inverse_variance",
    )

    rows = []
    for rank, fid in enumerate(fams, start=1):
        consistency = title_consistency(bulk_sum.loc[fid, "consistency_class"])
        consistency_fe = title_consistency(bulk_fe.loc[fid, "replication_class"])
        require(consistency == consistency_fe, f"F{fid} bulk consistency mismatch")
        rho = float(duration.loc[fid, "rho_primary"])
        dur_class = str(duration.loc[fid, "chronicity_class"])
        require(classify_duration(abs(rho)) == dur_class, f"F{fid} duration class != unrounded |rho|")
        sn_class = str(sn.loc[fid, "fidelity_class"])
        bulk_gate = consistency == "High"
        sn_gate = sn_class != "Below_detection"
        duration_gate = dur_class != "Stable"
        support = int(bulk_gate) + int(sn_gate) + int(duration_gate)
        rows.append(
            {
                "family_id": fid,
                "family_label": str(labels.loc[fid, "short_label"]),
                "family_direction": str(frozen.loc[fid, "canonical_direction"]),
                "canonical_order": rank,
                "hca_carrier": str(hca.loc[fid, "LA_carrier"]),
                "bulk_consistency": consistency,
                "bulk_dominant_state": str(bulk_fe.loc[fid, "dominant_state"]),
                "bulk_gate": bulk_gate,
                "snRNA_class": sn_class,
                "sn_gate": sn_gate,
                "duration_class": dur_class,
                "duration_gate": duration_gate,
                "support_count": support,
                "Path_A": support == 3,
            }
        )

    out = pd.DataFrame(rows)
    path_a = set(out.loc[out["Path_A"], "family_id"].astype(int))
    require(path_a == EXPECTED_PATH_A, f"Path A {path_a} != {EXPECTED_PATH_A}")
    tiers = Counter(support_tier(int(c)) for c in out["support_count"])
    for tier, n in EXPECTED_COUNTS.items():
        require(tiers[tier] == n, f"{tier} count {tiers[tier]} != {n}")
    out.to_csv(OUT_CANONICAL, sep="\t", index=False)
    return out


def strip_geometry() -> dict[str, float | list[float]]:
    content_w = WIDTH_MM - 2 * MARGIN_MM
    col_space = content_w - 4 * GAP_MM
    col_w = [col_space * fraction for fraction in COL_FRAC]
    require(abs(sum(COL_FRAC) - 1.0) < 1e-9, "COL_FRAC must sum to 1")
    xs = [MARGIN_MM]
    for width in col_w[:-1]:
        xs.append(xs[-1] + width + GAP_MM)
    geomx_y0 = MARGIN_MM
    strip_y0 = geomx_y0 + GEOMX_H_MM + STRIP_GAP_MM
    title_y1 = HEIGHT_MM - MARGIN_MM
    gap_centers = [xs[i] + col_w[i] + GAP_MM / 2 for i in range(4)]
    return {
        "xs": xs,
        "col_w": col_w,
        "strip_y0": strip_y0,
        "strip_cy": strip_y0 + STRIP_H_MM / 2,
        "geomx_y0": geomx_y0,
        "title_y1": title_y1,
        "gap_centers": gap_centers,
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


def add_axes(fig: plt.Figure, x_mm: float, y_mm: float, w_mm: float, h_mm: float) -> plt.Axes:
    ax = fig.add_axes([x_mm / WIDTH_MM, y_mm / HEIGHT_MM, w_mm / WIDTH_MM, h_mm / HEIGHT_MM])
    ax.set_xlim(0, w_mm)
    ax.set_ylim(0, h_mm)
    ax.axis("off")
    return ax


def place_line(
    ax: plt.Axes,
    x: float,
    y: float,
    text: str,
    *,
    fontsize: float = BODY_PT,
    fontweight: str = "normal",
    color: str = MUTED_COLOR,
    line_h: float | None = None,
) -> tuple[mpl.text.Text, float]:
    artist = ax.text(
        x,
        y,
        text,
        ha="left",
        va="top",
        fontsize=fontsize,
        fontweight=fontweight,
        color=color,
        clip_on=False,
    )
    step = line_h if line_h is not None else LINE_H_MM
    blank = 0.75 if HEIGHT_MM <= 48.5 else 0.9
    return artist, y - (step if text else blank)


def draw_header(ax: plt.Axes, h: float, title: str, pad_x: float, pad_y: float) -> float:
    y = h - pad_y
    header_step = 2.45 if HEIGHT_MM <= 48.5 else 2.55
    header_gap = 0.7 if HEIGHT_MM <= 48.5 else 0.85
    for line in title.split("\n"):
        ax.text(pad_x, y, line, ha="left", va="top", fontsize=KICKER_PT, color=KICKER_COLOR, fontweight="bold")
        y -= header_step
    return y - header_gap


def draw_mini_gate(ax: plt.Axes, x: float, y: float, w: float, h: float, label: str) -> None:
    ax.add_patch(
        mpl.patches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.08,rounding_size=0.35",
            facecolor=WHITE,
            edgecolor=GATE_EDGE,
            linewidth=0.55,
            clip_on=False,
        )
    )
    ax.text(
        x + w / 2,
        y + h / 2,
        label,
        ha="center",
        va="center",
        fontsize=max(5.2, FOOTNOTE_PT - 0.4),
        color=MUTED_COLOR,
        linespacing=1.15,
        clip_on=False,
    )


def text_bbox_figure_mm(artist: mpl.text.Text, ax: plt.Axes, fig: plt.Figure) -> tuple[float, float, float, float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    return (
        float(bbox.x0) * WIDTH_MM,
        float(bbox.y0) * HEIGHT_MM,
        float(bbox.x1) * WIDTH_MM,
        float(bbox.y1) * HEIGHT_MM,
    )


def text_bbox_mm(artist: mpl.text.Text, ax: plt.Axes, fig: plt.Figure) -> tuple[float, float, float, float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted())
    return float(bbox.x0), float(bbox.y0), float(bbox.x1), float(bbox.y1)


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


def render() -> tuple[plt.Figure, plt.Axes, list[plt.Axes], list[float], plt.Axes]:
    set_style()
    geom = strip_geometry()
    xs: list[float] = geom["xs"]  # type: ignore[assignment]
    col_w: list[float] = geom["col_w"]  # type: ignore[assignment]
    strip_y0: float = geom["strip_y0"]  # type: ignore[assignment]
    strip_cy: float = geom["strip_cy"]  # type: ignore[assignment]
    title_y1: float = geom["title_y1"]  # type: ignore[assignment]
    geomx_y0: float = geom["geomx_y0"]  # type: ignore[assignment]
    gap_centers: list[float] = geom["gap_centers"]  # type: ignore[assignment]

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM))
    fig.patch.set_facecolor(WHITE)
    pad_x = 1.4
    pad_y = 1.1
    h = STRIP_H_MM

    title_ax = fig.add_axes([0, 0, 1, 1])
    title_ax.set_xlim(0, WIDTH_MM)
    title_ax.set_ylim(0, HEIGHT_MM)
    title_ax.axis("off")
    title_ax.set_facecolor(WHITE)

    title_ax.text(MARGIN_MM, title_y1 - 0.3, "A", ha="left", va="top", fontsize=PANEL_LETTER_PT, fontweight="bold", color=FG_COLOR)
    main = title_ax.text(
        MARGIN_MM + 6.8,
        title_y1 - 0.3,
        "Cross-layer triangulation on the frozen NCP topology",
        ha="left",
        va="top",
        fontsize=TITLE_TEXT_PT,
        fontweight="bold",
        color=FG_COLOR,
    )
    fig.canvas.draw()
    main_bb = text_bbox_figure_mm(main, title_ax, fig)
    title_ax.text(
        MARGIN_MM + 6.8,
        main_bb[1] - 1.2,
        "Layer-specific summaries are juxtaposed without fitting a merged latent model.",
        ha="left",
        va="top",
        fontsize=SUBTITLE_PT,
        color=FAINT_COLOR,
    )

    col_axes = [add_axes(fig, xs[i], strip_y0, col_w[i], h) for i in range(5)]
    for index, ax in enumerate(col_axes):
        ax.set_facecolor(STRIP_BG)
        edge = TRIAGE_EDGE if index == 4 else MODULE_EDGE
        lw = 0.55 if index == 4 else 0.28
        ax.add_patch(
            mpl.patches.Rectangle(
                (0, 0),
                col_w[index],
                h,
                facecolor="none",
                edgecolor=edge,
                linewidth=lw,
                zorder=0,
            )
        )

    arrow_half = GAP_MM * 0.28
    for gap_x in gap_centers:
        bridge_arrow(title_ax, gap_x - arrow_half, gap_x + arrow_half, strip_cy)

    y = draw_header(col_axes[0], h, "FROZEN BACKBONE", pad_x, pad_y)
    for line in (
        "16 frozen NCP families",
        "Canonical membership,",
        "direction & display order",
        "from Fig 1D",
        "",
        "HCA LA carrier prior:",
        "context only",
    ):
        color = FAINT_COLOR if "context" in line else MUTED_COLOR
        _, y = place_line(col_axes[0], pad_x, y, line, color=color)

    y = draw_header(col_axes[1], h, "PAIRED-CHAMBER\nBULK", pad_x, pad_y)
    for line in (
        "Cross-cohort",
        "consistency",
        "Dominant rewiring",
        "state",
        "",
        "Path A gate:",
        "High consistency",
        "",
        "Dominant state:",
        "context only",
    ):
        weight = "bold" if line in {"Path A gate:", "High consistency"} else "normal"
        if weight == "bold":
            color = FG_COLOR
        elif "context" in line or line == "Dominant state:":
            color = FAINT_COLOR
        else:
            color = MUTED_COLOR
        _, y = place_line(col_axes[1], pad_x, y, line, fontweight=weight, color=color)

    y = draw_header(col_axes[2], h, "LA SINGLE-NUCLEUS", pad_x, pad_y)
    for line in (
        "Fidelity",
        "Relay",
        "Extension",
        "Below-detection",
        "",
        "Path A gate:",
        "non-Below-detection",
    ):
        weight = "bold" if line in {"Path A gate:", "non-Below-detection"} else "normal"
        color = FG_COLOR if weight == "bold" else MUTED_COLOR
        _, y = place_line(col_axes[2], pad_x, y, line, fontweight=weight, color=color)

    y = draw_header(col_axes[3], h, "LAA DURATION", pad_x, pad_y)
    for line in (
        "Strong-trend",
        "Weak-trend",
        "Stable",
        "",
        "Path A gate:",
        "non-Stable",
        "",
        "Class uses unrounded",
        RHO_LAA_MATH,
    ):
        weight = "bold" if line in {"Path A gate:", "non-Stable"} else "normal"
        if weight == "bold":
            color = FG_COLOR
        elif "unrounded" in line or line == RHO_LAA_MATH:
            color = FAINT_COLOR
        else:
            color = MUTED_COLOR
        fs = FOOTNOTE_PT if color == FAINT_COLOR else BODY_PT
        line_h = 2.7 if line == RHO_LAA_MATH else None
        _, y = place_line(col_axes[3], pad_x, y, line, fontweight=weight, color=color, fontsize=fs, line_h=line_h)

    y = draw_header(col_axes[4], h, "RULE-BASED TRIAGE", pad_x, pad_y)
    _, y = place_line(col_axes[4], pad_x, y, "Three Path A gates", fontweight="bold", color=FG_COLOR)
    gate_y = y - 5.0
    gate_w = (col_w[4] - 2 * pad_x - 2 * 1.0) / 3
    gate_h = 4.6
    for i, lab in enumerate(MINI_GATE_LABELS):
        draw_mini_gate(col_axes[4], pad_x + i * (gate_w + 1.0), gate_y, gate_w, gate_h, lab)
    y = gate_y - 1.2
    for line in (
        "support count = 0–3",
        "Path A = 3/3 gates",
        "",
        "Rule-based prioritization",
        "Not a unified significance",
        "test or causal validation",
    ):
        weight = "bold" if line.startswith("Path A") else "normal"
        if weight == "bold":
            color = FG_COLOR
        elif line.startswith("Not") or "causal" in line or line.startswith("test"):
            color = FAINT_COLOR
        else:
            color = MUTED_COLOR
        fs = FOOTNOTE_PT if color == FAINT_COLOR else BODY_PT
        _, y = place_line(col_axes[4], pad_x, y, line, fontweight=weight, color=color, fontsize=fs)

    geomx_ax = add_axes(fig, MARGIN_MM, geomx_y0, WIDTH_MM - 2 * MARGIN_MM, GEOMX_H_MM)
    geomx_ax.add_patch(
        mpl.patches.FancyBboxPatch(
            (0.3, 0.35),
            WIDTH_MM - 2 * MARGIN_MM - 0.6,
            GEOMX_H_MM - 0.7,
            boxstyle="round,pad=0.15,rounding_size=0.4",
            facecolor=WHITE,
            edgecolor=GEOMX_EDGE,
            linewidth=0.7,
            linestyle=(0, (2.2, 1.6)),
            clip_on=False,
        )
    )
    geomx_ax.text(
        1.6,
        GEOMX_H_MM / 2,
        "GeoMx spatial context · 24 eligible ROIs · 4 patients · ROI-level only · "
        "not included in support count or Path A · shown in D",
        ha="left",
        va="center",
        fontsize=GEOMX_PT,
        color=FAINT_COLOR,
        clip_on=False,
    )

    return fig, title_ax, col_axes, col_w, geomx_ax


def validate_layout(
    fig: plt.Figure,
    title_ax: plt.Axes,
    col_axes: list[plt.Axes],
    col_w: list[float],
    geomx_ax: plt.Axes,
) -> dict[str, float | int]:
    fig.canvas.draw()
    entries: list[tuple[str, tuple[float, float, float, float], str, float]] = []

    def collect(ax: plt.Axes, zone: str, width_limit: float | None = None, height_limit: float | None = None) -> None:
        for child in ax.get_children():
            if not isinstance(child, mpl.text.Text):
                continue
            label = child.get_text().strip()
            if not label:
                continue
            if zone in {"title", "geomx"}:
                bbox = text_bbox_figure_mm(child, ax, fig)
            else:
                local = text_bbox_mm(child, ax, fig)
                pos = ax.get_position()
                bbox = (
                    local[0] + float(pos.x0) * WIDTH_MM,
                    local[1] + float(pos.y0) * HEIGHT_MM,
                    local[2] + float(pos.x0) * WIDTH_MM,
                    local[3] + float(pos.y0) * HEIGHT_MM,
                )
                require(local[0] >= -BBOX_TOL_MM, f"clipped left {zone}: {label!r}")
                if width_limit is not None:
                    require(local[2] <= width_limit + BBOX_TOL_MM, f"clipped right {zone}: {label!r}")
                require(local[1] >= -BBOX_TOL_MM, f"clipped bottom {zone}: {label!r}")
                if height_limit is not None:
                    require(local[3] <= height_limit + BBOX_TOL_MM, f"clipped top {zone}: {label!r}")
            require(bbox[0] >= -BBOX_TOL_MM, f"figure clipped left: {label!r}")
            require(bbox[2] <= WIDTH_MM + BBOX_TOL_MM, f"figure clipped right: {label!r}")
            require(bbox[1] >= -BBOX_TOL_MM, f"figure clipped bottom: {label!r}")
            require(bbox[3] <= HEIGHT_MM + BBOX_TOL_MM, f"figure clipped top: {label!r}")
            fontsize = float(child.get_fontsize())
            is_mini_gate = label in MINI_GATE_LABELS
            min_pt = 5.15 if is_mini_gate else MIN_FONT_PT
            require(fontsize >= min_pt - 0.05, f"font below minimum: {label!r}")
            entries.append((label, bbox, zone, fontsize))

    collect(title_ax, "title")
    for i, ax in enumerate(col_axes):
        collect(ax, f"col{i + 1}", width_limit=col_w[i], height_limit=STRIP_H_MM)
    collect(geomx_ax, "geomx")

    for i, (label_a, bbox_a, zone_a, _) in enumerate(entries):
        for label_b, bbox_b, zone_b, _ in entries[i + 1 :]:
            if zone_a == "title" and zone_b == "title":
                continue
            if zone_a.startswith("col") and zone_b.startswith("col") and zone_a != zone_b:
                continue
            overlap = bbox_intersection(bbox_a, bbox_b)
            require(overlap <= OVERLAP_TOL_MM2, f"overlapping text: {label_a!r} vs {label_b!r} ({overlap:.2f})")

    joined = "\n".join(label for label, *_ in entries)
    joined_plain = joined.replace(RHO_LAA_MATH, "|rho_LAA|")
    for phrase in (
        "Cross-layer triangulation on the frozen NCP topology",
        "FROZEN BACKBONE",
        "PAIRED-CHAMBER",
        "LA SINGLE-NUCLEUS",
        "LAA DURATION",
        "RULE-BASED TRIAGE",
        "High consistency",
        "non-Below-detection",
        "non-Stable",
        "Path A = 3/3 gates",
        "GeoMx spatial context",
        "not included in support count or Path A",
        "context only",
        "Cross-cohort",
        "Dominant rewiring",
        "Class uses unrounded",
        "Not a unified significance",
        "test or causal validation",
        "Three Path A gates",
        "Bulk\nHigh",
        "snRNA\nnon-Below",
        "Duration\nnon-Stable",
    ):
        require(phrase in joined, f"required phrase missing: {phrase}")
    require(RHO_LAA_MATH in joined, "required formal |rho_LAA| mathtext missing")
    require("BULK" in joined, "required phrase missing: BULK")

    path_a_gate_hits = sum(1 for label, *_ in entries if label == "Path A gate:")
    require(path_a_gate_hits == 3, f"expected exactly 3 Path A gate labels, got {path_a_gate_hits}")
    for phrase in PATH_A_GATE_PHRASES:
        require(phrase in joined, f"Path A gate condition missing: {phrase}")
    for lab in MINI_GATE_LABELS:
        require(lab in joined, f"mini-gate missing: {lab!r}")

    require("NCP_v3" not in joined_plain, "NCP_v3 visible")
    require("PriorityScore" not in joined_plain, "PriorityScore visible")
    require("PASS" not in joined_plain and "FAIL" not in joined_plain, "PASS/FAIL checks visible")
    require(re.search(r"\bF(?:1[0-6]|[1-9])\b", joined_plain) is None, "family IDs visible in panel A")
    require(re.search(r"support\s*=\s*\d", joined_plain, flags=re.IGNORECASE) is None, "result support counts visible")
    require(
        all(label != "non-Below" for label, *_ in entries),
        "stale Path A gate wording non-Below (use non-Below-detection)",
    )
    require("Not unified significance" not in joined_plain, "stale significance disclaimer")
    require("LAA |ρ|" not in joined_plain and "LAA |rho|" not in joined_plain.lower(), "stale LAA |rho| wording")

    # HCA / dominant-state / GeoMx must remain context-only, never Path A conditions.
    gate_block = "\n".join(
        label
        for label, *_ in entries
        if label.startswith("Path A")
        or label in PATH_A_GATE_PHRASES
        or label in MINI_GATE_LABELS
        or label == "Three Path A gates"
    )
    for banned in ("HCA", "Dominant", "GeoMx", "carrier", "rewiring"):
        require(banned not in gate_block, f"{banned} appears inside Path A gate block")

    return {
        "n_text_objects": len(entries),
        "min_font_pt": float(min(fs for *_, fs in entries)),
        "n_path_a_gate_labels": path_a_gate_hits,
    }


def assert_export_geometry(pdf: Path, svg: Path, png: Path) -> dict[str, str]:
    """Confirm PDF/SVG/PNG share the same physical canvas geometry."""
    from PIL import Image

    svg_text = svg.read_text(encoding="utf-8")
    require("NCP_v3" not in svg_text, "NCP_v3 in SVG")
    require("PriorityScore" not in svg_text, "PriorityScore in SVG")
    require("GeoMx spatial context" in svg_text, "GeoMx strip missing in SVG")
    require("non-Below-detection" in svg_text, "polished snRNA gate missing in SVG")
    require("Not a unified significance" in svg_text, "polished disclaimer missing in SVG")
    require("Bulk" in svg_text and "High" in svg_text, "Bulk/High mini-gate missing in SVG")
    require("snRNA" in svg_text and "non-Below" in svg_text, "snRNA/non-Below mini-gate missing in SVG")
    require("Duration" in svg_text and "non-Stable" in svg_text, "Duration/non-Stable mini-gate missing in SVG")
    require(re.search(r"\bF(?:1[0-6]|[1-9])\b", svg_text) is None, "family IDs in SVG")

    # SVG viewBox / width-height in mm or pt
    vb = re.search(r'viewBox="([0-9.eE+-]+) ([0-9.eE+-]+) ([0-9.eE+-]+) ([0-9.eE+-]+)"', svg_text)
    require(vb is not None, "SVG viewBox missing")
    svg_w, svg_h = float(vb.group(3)), float(vb.group(4))
    # Matplotlib SVG often uses points (1 pt = 1/72 in). Convert to mm.
    svg_w_mm = svg_w * 25.4 / 72.0
    svg_h_mm = svg_h * 25.4 / 72.0
    require(abs(svg_w_mm - WIDTH_MM) < 0.8, f"SVG width {svg_w_mm:.2f} mm != {WIDTH_MM}")
    require(abs(svg_h_mm - HEIGHT_MM) < 0.8, f"SVG height {svg_h_mm:.2f} mm != {HEIGHT_MM}")

    img = Image.open(png)
    png_w_px, png_h_px = img.size
    png_w_mm = png_w_px * 25.4 / 600.0
    png_h_mm = png_h_px * 25.4 / 600.0
    require(abs(png_w_mm - WIDTH_MM) < 0.8, f"PNG width {png_w_mm:.2f} mm != {WIDTH_MM}")
    require(abs(png_h_mm - HEIGHT_MM) < 0.8, f"PNG height {png_h_mm:.2f} mm != {HEIGHT_MM}")

    # PDF page size from MediaBox (points).
    pdf_bytes = pdf.read_bytes()
    media = re.search(
        rb"/MediaBox\s*\[\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*\]",
        pdf_bytes,
    )
    require(media is not None, "PDF MediaBox missing")
    pdf_w_mm = (float(media.group(3)) - float(media.group(1))) * 25.4 / 72.0
    pdf_h_mm = (float(media.group(4)) - float(media.group(2))) * 25.4 / 72.0
    require(abs(pdf_w_mm - WIDTH_MM) < 0.8, f"PDF width {pdf_w_mm:.2f} mm != {WIDTH_MM}")
    require(abs(pdf_h_mm - HEIGHT_MM) < 0.8, f"PDF height {pdf_h_mm:.2f} mm != {HEIGHT_MM}")
    require(abs(pdf_w_mm - svg_w_mm) < 0.8, "PDF/SVG width mismatch")
    require(abs(pdf_h_mm - svg_h_mm) < 0.8, "PDF/SVG height mismatch")
    require(abs(pdf_w_mm - png_w_mm) < 0.8, "PDF/PNG width mismatch")
    require(abs(pdf_h_mm - png_h_mm) < 0.8, "PDF/PNG height mismatch")

    return {
        "pdf_mm": f"{pdf_w_mm:.2f}x{pdf_h_mm:.2f}",
        "svg_mm": f"{svg_w_mm:.2f}x{svg_h_mm:.2f}",
        "png_mm": f"{png_w_mm:.2f}x{png_h_mm:.2f}",
        "png_px": f"{png_w_px}x{png_h_px}",
    }


def write_source_and_qc(
    manifest: pd.DataFrame,
    canonical: pd.DataFrame,
    layout_stats: dict[str, float | int],
    output_hashes: dict[str, str],
    geometry: dict[str, str],
) -> None:
    source = pd.DataFrame(
        [
            {"metric": "panel", "value": "Figure_5A"},
            {"metric": "canvas_mm", "value": f"{WIDTH_MM:.0f}x{HEIGHT_MM:.0f}"},
            {"metric": "module_widths_frac", "value": "/".join(str(int(f * 100)) for f in COL_FRAC)},
            {"metric": "canonical_order_source", "value": "Figure_1D_ncp_v3_canonical_display_order"},
            {"metric": "n_Path_A_gates", "value": 3},
            {"metric": "Path_A_gate_conditions", "value": "High consistency; non-Below-detection; non-Stable"},
            {"metric": "mini_gate_labels", "value": "Bulk/High; snRNA/non-Below; Duration/non-Stable"},
            {"metric": "Path_A_families", "value": "F10,F15"},
            {"metric": "triple_double_single_silent", "value": "2/1/10/3"},
            {"metric": "hca_is_gate", "value": False},
            {"metric": "dominant_state_is_gate", "value": False},
            {"metric": "geomx_is_gate", "value": False},
            {"metric": "genetics_is_gate", "value": False},
            {"metric": "priorityscore_opened", "value": False},
            {"metric": "legacy_fig5_table_read", "value": False},
            {"metric": "family_ids_visible", "value": False},
            {"metric": "n_text_objects", "value": layout_stats["n_text_objects"]},
            {"metric": "min_font_pt", "value": layout_stats["min_font_pt"]},
            {"metric": "pdf_geometry_mm", "value": geometry["pdf_mm"]},
            {"metric": "svg_geometry_mm", "value": geometry["svg_mm"]},
            {"metric": "png_geometry_mm", "value": geometry["png_mm"]},
            {"metric": "canonical_table", "value": str(OUT_CANONICAL.relative_to(ROOT))},
            {"metric": "manifest_table", "value": str(OUT_MANIFEST.relative_to(ROOT))},
        ]
    )
    source.to_csv(OUT_SOURCE, sep="\t", index=False)

    path_a = sorted(canonical.loc[canonical["Path_A"], "family_id"].astype(int))
    OUT_QC.write_text(
        f"""# Figure 5A triangulation contract QC

Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}

Status: **PASS**

## Scope
- Panel A only. Figure 5B–D and formal submission routing were not modified.
- No legacy `results/fig5b_v3_triangulation.tsv` read.
- No Atlas PriorityScore field opened.

## Final terminology / publication polish
- snRNA Path A gate wording: `non-Below-detection`
- Duration class note: formal mathtext `|ρ_LAA|`
- Disclaimer: `Not a unified significance test or causal validation`
- Mini-gates: Bulk/High; snRNA/non-Below; Duration/non-Stable (neutral outlines)
- Ordinary module outlines lightened; RULE-BASED TRIAGE modestly stronger
- GeoMx dashed context strip retained

## Final QC checks
- Exactly three Path A gates displayed: **PASS** (`{layout_stats.get("n_path_a_gate_labels", 3)}` labels)
- HCA / dominant-state / GeoMx never presented as gates: **PASS**
- No PriorityScore text: **PASS**
- No result family IDs or per-family support counts in panel A: **PASS**
- All visible text inside canvas: **PASS**
- No visible `NCP_v3`: **PASS**
- PDF / SVG / PNG geometries match: **PASS** (PDF {geometry["pdf_mm"]} mm; SVG {geometry["svg_mm"]} mm; PNG {geometry["png_mm"]} mm / {geometry["png_px"]} px @ 600 dpi)

## Canonical triangulation
- Order source: Figure 1D (`ncp_v3_canonical_display_order.tsv`)
- Path A gates: exactly 3 (bulk High; snRNA non-Below-detection; duration non-Stable)
- Non-gates: HCA carrier, dominant rewiring state, GeoMx, genetics
- Path A families: {path_a}
- Counts triple/double/single/silent: 2/1/10/3

## Layout
- Canvas: {WIDTH_MM:.0f} × {HEIGHT_MM:.0f} mm
- Module fractions: {COL_FRAC}
- Visible text objects: {layout_stats["n_text_objects"]}
- Minimum font: {float(layout_stats["min_font_pt"]):.1f} pt
- GeoMx strip: dashed context only

## Inputs
{chr(10).join(f"- `{r.artifact_role}`: `{r.artifact_path}` — `{r.sha256}`" for r in manifest.itertuples(index=False))}

## Outputs
{chr(10).join(f"- `{path}` — `{digest}`" for path, digest in sorted(output_hashes.items()))}
""",
        encoding="utf-8",
    )


def main() -> None:
    assert_no_forbidden_reads()
    require(STRIP_H_MM > 25.0, f"Strip height too small: {STRIP_H_MM:.2f}")

    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    manifest = write_manifest()
    canonical = build_canonical_triangulation()

    fig, title_ax, col_axes, col_w, geomx_ax = render()
    layout_stats = validate_layout(fig, title_ax, col_axes, col_w, geomx_ax)

    pdf = OUT_BASE.with_suffix(".pdf")
    svg = OUT_BASE.with_suffix(".svg")
    png = OUT_BASE.with_suffix(".png")
    fig.savefig(pdf, bbox_inches=None, pad_inches=0)
    fig.savefig(svg, bbox_inches=None, pad_inches=0)
    fig.savefig(png, dpi=600, bbox_inches=None, pad_inches=0)
    fig.savefig(OUT_PREVIEW, dpi=150, bbox_inches=None, pad_inches=0)
    plt.close(fig)

    geometry = assert_export_geometry(pdf, svg, png)

    output_hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (pdf, svg, png, OUT_PREVIEW, OUT_MANIFEST, OUT_CANONICAL)
    }
    write_source_and_qc(manifest, canonical, layout_stats, output_hashes, geometry)
    output_hashes[str(OUT_SOURCE.relative_to(ROOT))] = sha256(OUT_SOURCE)
    write_source_and_qc(manifest, canonical, layout_stats, output_hashes, geometry)

    print(f"Wrote {OUT_MANIFEST.relative_to(ROOT)}")
    print(f"Wrote {OUT_CANONICAL.relative_to(ROOT)}")
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")
    print(f"Geometry match: PDF={geometry['pdf_mm']} SVG={geometry['svg_mm']} PNG={geometry['png_mm']} mm")


if __name__ == "__main__":
    main()