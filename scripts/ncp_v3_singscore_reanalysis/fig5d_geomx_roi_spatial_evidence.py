#!/usr/bin/env python3
"""Figure 5D: GeoMx ROI-level spatial evidence table (three-row aligned layout).

Three example rows link NCP family signatures to GeoMx spatial context across the
24 eligible cardiac AA DNA+ ROIs. GeoMx is descriptive corroboration only and is
not used as a Path A gate.

Input:
    Project_P5/GSE261363/roi_scores_v3.tsv

Outputs:
    figures/ncp_v3_singscore_reanalysis/fig5d_geomx_roi_spatial_evidence.{pdf,svg,png}
    results/ncp_v3_singscore_reanalysis/fig5d_geomx_roi_spatial_evidence_source.tsv
    results/ncp_v3_singscore_reanalysis/fig5d_geomx_roi_spatial_evidence_summary.tsv
    results/ncp_v3_singscore_reanalysis/fig5d_geomx_roi_spatial_evidence_qc.md
"""
from __future__ import annotations

import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
IN_GEO = ROOT / "Project_P5" / "GSE261363" / "roi_scores_v3.tsv"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
OUT_BASE = FIGURES / "fig5d_geomx_roi_spatial_evidence"
OUT_SOURCE = RESULTS / "fig5d_geomx_roi_spatial_evidence_source.tsv"
OUT_SUMMARY = RESULTS / "fig5d_geomx_roi_spatial_evidence_summary.tsv"
OUT_QC = RESULTS / "fig5d_geomx_roi_spatial_evidence_qc.md"
OUT_PREVIEW = OUT_BASE.with_name(OUT_BASE.name + "_preview.png")

WIDTH_MM = 98.0
HEIGHT_MM = 88.0
MARGIN_MM = 2.0
MM = 1.0 / 25.4

# Absolute type system shared with Figure 5C in the composite.
PANEL_LETTER_PT = 10.8
TITLE_PT = 8.6
SUBTITLE_PT = 6.5
BODY_PT = 6.5
NOTE_PT = 6.0
TICK_PT = 6.0
FOOTER_PT = 5.5

_STANDALONE = {
    "WIDTH_MM": 98.0,
    "HEIGHT_MM": 88.0,
    "MARGIN_MM": 2.0,
    "PANEL_LETTER_PT": 10.8,
    "TITLE_PT": 8.6,
    "SUBTITLE_PT": 6.5,
    "BODY_PT": 6.5,
    "NOTE_PT": 6.0,
    "TICK_PT": 6.0,
    "FOOTER_PT": 5.5,
}
_COMPOSITE = {
    "WIDTH_MM": 98.0,
    "HEIGHT_MM": 88.0,
    "MARGIN_MM": 2.5,  # match C/D letter baseline with composite C margin_top
    # One-step type bump for composite readability
    "PANEL_LETTER_PT": 11.0,
    "TITLE_PT": 9.0,
    "SUBTITLE_PT": 7.0,
    "BODY_PT": 6.7,
    "NOTE_PT": 6.5,
    "TICK_PT": 6.5,
    "FOOTER_PT": 5.8,
}


def apply_layout(mode: str = "standalone") -> None:
    """Switch between standalone and Figure 5 composite canvas geometry."""
    global WIDTH_MM, HEIGHT_MM, MARGIN_MM
    global PANEL_LETTER_PT, TITLE_PT, SUBTITLE_PT, BODY_PT, NOTE_PT, TICK_PT, FOOTER_PT
    cfg = _COMPOSITE if mode == "composite" else _STANDALONE
    WIDTH_MM = cfg["WIDTH_MM"]
    HEIGHT_MM = cfg["HEIGHT_MM"]
    MARGIN_MM = cfg["MARGIN_MM"]
    PANEL_LETTER_PT = cfg["PANEL_LETTER_PT"]
    TITLE_PT = cfg["TITLE_PT"]
    SUBTITLE_PT = cfg["SUBTITLE_PT"]
    BODY_PT = cfg["BODY_PT"]
    NOTE_PT = cfg["NOTE_PT"]
    TICK_PT = cfg["TICK_PT"]
    FOOTER_PT = cfg["FOOTER_PT"]

FG = "#222222"
MUTED = "#666666"
FAINT = "#888888"
WHITE = "#FFFFFF"
BG_ALT = "#F7F7F7"
ROW_RULE = "#E8E8E8"
MEAN_LINE = "#444444"
AXIS_COLOR = "#B0B0B0"
AXIS_LW = 0.45
ZERO_LINE = "#C8C8C8"
ZERO_LW = 0.4
GRID_COLOR = "#EEEEEE"
SR_EDGE_LW = 0.80  # ~14% above prior 0.70 for grayscale readability
AF_EDGE_LW = 0.70

# Local Figure 5D semantics only — do not reuse for evidence classes elsewhere.
# Softened Okabe–Ito-style patient palette, distinct from Fig 5B LA/RA,
# bulk-state, snRNA-class and Path A greys.
PATIENT_ORDER = ["P161", "P167", "P175", "P187"]
PATIENT_COLORS = {
    "P161": "#4A8BB8",  # soft sky blue
    "P167": "#C9A03A",  # soft gold
    "P175": "#3A8A72",  # soft teal
    "P187": "#B06A90",  # soft rose
}
PATIENT_SHAPES = {
    "P161": "o",
    "P167": "^",
    "P175": "s",
    "P187": "D",
}
PATIENT_SHORT = {
    "P161": "P1",
    "P167": "P2",
    "P175": "P3",
    "P187": "P4",
}

FORBIDDEN_READ_PATTERNS = [
    re.compile(r"PriorityScore"),
    re.compile(r"figure5_triangulation_(canonical|legacy)"),
    re.compile(r"g3b|g3a|fig2|fig3|fig4|fig6", re.I),
    re.compile(r"genetics|genotype|variant", re.I),
    re.compile(r"hca.*carrier|carrier.*hca", re.I),
    re.compile(r"bulk.*state|dominant.*state", re.I),
]


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_no_forbidden_reads() -> None:
    script_text = Path(__file__).read_text(encoding="utf-8")
    for i, line in enumerate(script_text.splitlines(), 1):
        if "pd.read_csv(" not in line:
            continue
        for pat in FORBIDDEN_READ_PATTERNS:
            if pat.search(line):
                raise AssertionError(f"Forbidden read on line {i}: {line.strip()}")


def darken(hex_color: str, factor: float = 0.72) -> str:
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    r = max(0, min(255, int(r * factor)))
    g = max(0, min(255, int(g * factor)))
    b = max(0, min(255, int(b * factor)))
    return f"#{r:02X}{g:02X}{b:02X}"


def patient_style(patient: str, rhythm: str) -> tuple[str, str, str, float]:
    color = PATIENT_COLORS[patient]
    marker = PATIENT_SHAPES[patient]
    if rhythm == "SR":
        return marker, WHITE, color, SR_EDGE_LW
    return marker, color, darken(color), AF_EDGE_LW


def load_and_assert() -> pd.DataFrame:
    assert_no_forbidden_reads()
    require(IN_GEO.is_file(), f"GeoMx source not found: {IN_GEO}")
    df = pd.read_csv(IN_GEO, sep="\t")
    require("eligible_for_5d" in df.columns, "missing eligible_for_5d column")
    df = df[df["eligible_for_5d"] == True].copy()  # noqa: E712

    n_roi = df["gsm"].nunique()
    require(n_roi == 24, f"expected 24 eligible ROIs, got {n_roi}")

    pats = sorted(df["patient"].unique())
    require(len(pats) == 4, f"expected 4 patients, got {len(pats)}")
    require(
        set(pats) == set(PATIENT_ORDER),
        f"unexpected patients: {pats}",
    )

    markers = set(df["marker"].unique())
    require(
        markers <= {"DNA+", "DNA"},
        f"all eligible ROIs must be DNA-positive, got markers: {markers}",
    )
    require(
        all(df["tissue"] == "Atrial appendage"),
        "all eligible ROIs must be atrial appendage",
    )
    return df


def welch_ci(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan"), float("nan")
    diff = a.mean() - b.mean()
    se = float(np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)))
    if se == 0:
        return diff, 0.0, 1.0
    df = (a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) ** 2 / (
        (a.var(ddof=1) ** 2) / (len(a) ** 2 * (len(a) - 1))
        + (b.var(ddof=1) ** 2) / (len(b) ** 2 * (len(b) - 1))
    )
    p = 2 * (1 - stats.t.cdf(abs(diff) / se, df))
    tcrit = stats.t.ppf(0.975, df)
    return diff, tcrit * se, p


def build_sources(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, Any]], dict[str, str]]:
    W = df.pivot_table(index="signature_lab", columns="gsm", values="score")
    meta = (
        df[["gsm", "patient", "rhythm", "area"]]
        .drop_duplicates("gsm")
        .set_index("gsm")
        .loc[W.columns]
    )

    markers = {p: PATIENT_SHAPES[p] for p in PATIENT_ORDER}

    r1_rho, r1_p = stats.spearmanr(
        W.loc["F10  TGF-beta signalling"], W.loc["fibroblast"]
    )
    r2_rho, r2_p = stats.spearmanr(
        W.loc["F15  NTRK2 / TrkB"], W.loc["macrophage"]
    )

    af = W.loc["F9   HS-GAG metabolism", meta["rhythm"] == "AF"].dropna().values
    sr = W.loc["F9   HS-GAG metabolism", meta["rhythm"] == "SR"].dropna().values
    r3_delta, r3_ci, r3_p = welch_ci(af, sr)

    require(abs(r1_rho - 0.439130) < 0.001, f"F10 rho mismatch: {r1_rho}")
    require(abs(r1_p - 0.031804) < 0.001, f"F10 p mismatch: {r1_p}")
    require(abs(r2_rho - 0.431304) < 0.001, f"F15 rho mismatch: {r2_rho}")
    require(abs(r2_p - 0.035350) < 0.001, f"F15 p mismatch: {r2_p}")
    require(abs(r3_delta - (-0.033333)) < 0.001, f"F9 delta mismatch: {r3_delta}")
    require(abs(r3_p - 0.048130) < 0.001, f"F9 p mismatch: {r3_p}")

    rows: list[dict[str, Any]] = [
        {
            "name": "F10",
            "family_id": "F10",
            "family_name": "GAG disease /\nTGF-β",
            "context": "Fibroblast-panel\nassociation",
            "plot_type": "scatter",
            "x": W.loc["fibroblast"],
            "y": W.loc["F10  TGF-beta signalling"],
            "x_context": "x = fibroblast-panel score",
            "stat_type": "Spearman",
            "rho": r1_rho,
            "p": r1_p,
            "n": 24,
        },
        {
            "name": "F15",
            "family_id": "F15",
            "family_name": "NTRK2 / TrkB",
            "context": "Macrophage-panel\nassociation",
            "plot_type": "scatter",
            "x": W.loc["macrophage"],
            "y": W.loc["F15  NTRK2 / TrkB"],
            "x_context": "x = macrophage-panel score",
            "stat_type": "Spearman",
            "rho": r2_rho,
            "p": r2_p,
            "n": 24,
        },
        {
            "name": "F9",
            "family_id": "F9",
            "family_name": "HS-GAG / bone /\nGPCR",
            "context": "AF–SR ROI\ncontrast",
            "plot_type": "strip",
            "y": W.loc["F9   HS-GAG metabolism"],
            "x_context": "",
            "stat_type": "Welch",
            "delta": r3_delta,
            "ci": r3_ci,
            "p": r3_p,
            "n": 24,
            "sr_n": len(sr),
            "af_n": len(af),
        },
    ]
    return W, meta, rows, markers


def write_source(W: pd.DataFrame, meta: pd.DataFrame) -> None:
    source = W.T.reset_index().rename(columns={"gsm": "roi_id"})
    source["patient"] = source["roi_id"].map(meta["patient"])
    source["rhythm"] = source["roi_id"].map(meta["rhythm"])
    source["area"] = source["roi_id"].map(meta["area"])
    source.to_csv(OUT_SOURCE, sep="\t", index=False)


def write_summary(
    W: pd.DataFrame, meta: pd.DataFrame, rows: list[dict[str, Any]]
) -> pd.DataFrame:
    out_rows: list[dict[str, Any]] = []
    for row in rows:
        d: dict[str, Any] = {
            "name": row["name"],
            "family_id": row["family_id"],
            "family_name": row["family_name"],
            "context": row["context"],
            "stat_type": row["stat_type"],
            "n": row["n"],
        }
        if "rho" in row:
            d["rho"] = row["rho"]
            d["p"] = row["p"]
        else:
            d["delta"] = row["delta"]
            d["ci"] = row["ci"]
            d["p"] = row["p"]
            d["sr_n"] = row.get("sr_n", 0)
            d["af_n"] = row.get("af_n", 0)
        out_rows.append(d)
    summary = pd.DataFrame(out_rows)
    summary["n_rois"] = len(W.columns)
    summary["n_patients"] = meta["patient"].nunique()
    summary["assertion_status"] = "PASS"
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    return summary


def _draw_summary_table(
    ax: plt.Axes,
    x_label: float,
    x_value: float,
    y_top: float,
    header: str,
    pairs: list[tuple[str, str]],
) -> None:
    line_h = 3.15
    ax.text(
        x_label,
        y_top,
        header,
        ha="left",
        va="center",
        fontsize=NOTE_PT,
        color=FG,
        fontweight="semibold",
    )
    for i, (lab, val) in enumerate(pairs):
        y = y_top - (i + 1) * line_h
        ax.text(
            x_label,
            y,
            lab,
            ha="left",
            va="center",
            fontsize=NOTE_PT,
            color=MUTED,
        )
        ax.text(
            x_value,
            y,
            val,
            ha="right",
            va="center",
            fontsize=BODY_PT,
            color=FG,
            fontfamily="monospace",
        )


def render(
    W: pd.DataFrame, meta: pd.DataFrame, rows: list[dict[str, Any]], markers: dict[str, str]
) -> tuple[plt.Figure, plt.Axes]:
    del markers  # canonical mapping is PATIENT_SHAPES
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Liberation Sans"]
    mpl.rcParams["axes.unicode_minus"] = True
    mpl.rcParams["pdf.fonttype"] = 42

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM), facecolor=WHITE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, WIDTH_MM)
    ax.set_ylim(0, HEIGHT_MM)
    ax.axis("off")
    ax.set_facecolor(WHITE)

    # Column geometry (approx 25 / 48 / 21 mm with fixed gutters)
    gutter = 0.75
    content_w = WIDTH_MM - 2 * MARGIN_MM
    label_w = 23.0
    summary_w = 25.5
    plot_w = content_w - label_w - summary_w - 2 * gutter
    require(plot_w > 40.0, f"plot column too narrow: {plot_w:.2f}")
    plot_left = MARGIN_MM + label_w + gutter
    summary_left = plot_left + plot_w + gutter
    label_x = MARGIN_MM + 0.8
    plot_center_x = plot_left + plot_w / 2
    summary_label_x = summary_left + 0.3
    summary_value_x = summary_left + summary_w - 0.4

    # Title block — panel letter / title sizes match Figure 5C
    y_title = HEIGHT_MM - MARGIN_MM - 0.2
    # Composite type bump: keep subtitle below title with va=top (not center).
    subtitle_dy = 3.8 if MARGIN_MM >= 2.5 and TITLE_PT >= 8.9 else 5.6
    subtitle_va = "top" if subtitle_dy < 5.0 else "center"
    header_dy = 12.2 if subtitle_va == "top" else 11.0
    ax.text(
        MARGIN_MM,
        y_title,
        "D",
        ha="left",
        va="top",
        fontsize=PANEL_LETTER_PT,
        fontweight="bold",
        color=FG,
    )
    ax.text(
        MARGIN_MM + 6.5,
        y_title,
        "GeoMx ROI spatial context",
        ha="left",
        va="top",
        fontsize=TITLE_PT,
        fontweight="bold",
        color=FG,
    )
    ax.text(
        MARGIN_MM + 6.5,
        y_title - subtitle_dy,
        "24 eligible cardiac AA DNA-positive ROIs\n· 4 patients · descriptive context only",
        ha="left",
        va=subtitle_va,
        fontsize=SUBTITLE_PT,
        color=FAINT,
        linespacing=0.9,
    )
    # Path A exclusion is stated in the Figure 5 legend, not as an in-panel badge.

    # Column headers
    header_y = y_title - header_dy
    ax.text(
        label_x,
        header_y,
        "Family /\nspatial context",
        ha="left",
        va="center",
        fontsize=NOTE_PT,
        color=MUTED,
        fontweight="semibold",
        linespacing=0.9,
    )
    ax.text(
        plot_center_x,
        header_y,
        "ROI-level data",
        ha="center",
        va="center",
        fontsize=NOTE_PT,
        color=MUTED,
        fontweight="semibold",
    )
    ax.text(
        summary_label_x,
        header_y,
        "ROI-level\nsummary",
        ha="left",
        va="center",
        fontsize=NOTE_PT,
        color=MUTED,
        fontweight="semibold",
        linespacing=0.9,
    )
    ax.text(
        plot_center_x,
        header_y - 3.0,
        "y = GeoMx family-signature score",
        ha="center",
        va="center",
        fontsize=NOTE_PT,
        color=FAINT,
    )

    # Row band
    legend_top = 13.0
    row_top = header_y - 4.2
    n_rows = 3
    row_gap = 0.8
    available = row_top - legend_top
    row_h = (available - (n_rows - 1) * row_gap) / n_rows
    patients = PATIENT_ORDER

    plot_axes: list[plt.Axes] = []
    for idx, row in enumerate(rows):
        row_y = row_top - idx * (row_h + row_gap)
        row_bottom = row_y - row_h
        row_center = (row_y + row_bottom) / 2

        if idx % 2 == 0:
            ax.add_patch(
                Rectangle(
                    (MARGIN_MM, row_bottom),
                    WIDTH_MM - 2 * MARGIN_MM,
                    row_h,
                    facecolor=BG_ALT,
                    edgecolor="none",
                    zorder=0,
                )
            )
        if idx < n_rows - 1:
            sep_y = row_bottom - row_gap / 2
            ax.plot(
                [MARGIN_MM, WIDTH_MM - MARGIN_MM],
                [sep_y, sep_y],
                color=ROW_RULE,
                lw=0.4,
                zorder=0,
            )

        # Family column: ID / name / context (kept inside the 25 mm column)
        ax.text(
            label_x,
            row_center + 5.0,
            row["family_id"],
            ha="left",
            va="center",
            fontsize=BODY_PT,
            color=FG,
            fontweight="semibold",
        )
        ax.text(
            label_x,
            row_center + 1.6,
            row["family_name"],
            ha="left",
            va="center",
            fontsize=NOTE_PT,
            color=FG,
            fontweight="semibold",
            linespacing=0.9,
        )
        ax.text(
            label_x,
            row_center - 3.2,
            row["context"],
            ha="left",
            va="center",
            fontsize=NOTE_PT,
            color=MUTED,
            linespacing=0.9,
        )

        # Small x-context above the plot (not an axis title)
        if row["x_context"]:
            ax.text(
                plot_center_x,
                row_y - 1.6,
                row["x_context"],
                ha="center",
                va="center",
                fontsize=NOTE_PT,
                color=FAINT,
            )

        # Plot axes — leave headroom for x-context; reserve inset for tick labels
        # so they stay inside the ROI-level plot column / row band.
        plot_pad_top = 3.4 if row["x_context"] else 1.8
        plot_pad_bottom = 3.6 if row["plot_type"] == "strip" else 3.4
        ytick_inset = 4.0
        ax_left = (plot_left + ytick_inset) / WIDTH_MM
        ax_bottom = (row_bottom + plot_pad_bottom) / HEIGHT_MM
        ax_width = (plot_w - ytick_inset - 0.5) / WIDTH_MM
        ax_height = (row_h - plot_pad_top - plot_pad_bottom) / HEIGHT_MM
        pax = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])
        plot_axes.append(pax)

        n_points = 0
        if row["plot_type"] == "scatter":
            x = row["x"]
            y = row["y"]
            for patient in patients:
                pmask = meta["patient"] == patient
                for rhythm in ("SR", "AF"):
                    mask = pmask & (meta["rhythm"] == rhythm)
                    if not mask.any():
                        continue
                    marker, face, edge, lw = patient_style(patient, rhythm)
                    pax.scatter(
                        x[mask],
                        y[mask],
                        marker=marker,
                        s=14,
                        facecolors=face,
                        edgecolors=edge,
                        linewidths=lw,
                        zorder=2,
                    )
                    n_points += int(mask.sum())
            pax.set_xlabel("")
            pax.set_ylabel("")
        else:
            y_all = row["y"]
            group_ys: dict[str, list[float]] = {"SR": [], "AF": []}
            for g_idx, rhythm in enumerate(("SR", "AF")):
                rmask = meta["rhythm"] == rhythm
                for p_idx, patient in enumerate(patients):
                    mask = rmask & (meta["patient"] == patient)
                    if not mask.any():
                        continue
                    jitter = (p_idx - 1.5) * 0.11
                    x_vals = np.full(mask.sum(), float(g_idx) + jitter)
                    y_vals = y_all[mask].values
                    marker, face, edge, lw = patient_style(patient, rhythm)
                    pax.scatter(
                        x_vals,
                        y_vals,
                        marker=marker,
                        s=14,
                        facecolors=face,
                        edgecolors=edge,
                        linewidths=lw,
                        zorder=2,
                    )
                    group_ys[rhythm].extend(y_vals.tolist())
                    n_points += int(mask.sum())
            for g_idx, rhythm in enumerate(("SR", "AF")):
                ys = group_ys[rhythm]
                if ys:
                    mean_y = float(np.mean(ys))
                    pax.plot(
                        [g_idx - 0.22, g_idx + 0.22],
                        [mean_y, mean_y],
                        color=MEAN_LINE,
                        lw=1.0,
                        zorder=1,
                    )
            pax.set_xticks([0, 1])
            pax.set_xticklabels(["SR", "AF"], fontsize=TICK_PT, color=FAINT)
            pax.set_xlim(-0.55, 1.55)
            pax.set_xlabel("")
            pax.set_ylabel("")

        require(n_points == 24, f"{row['name']}: expected 24 points, got {n_points}")

        # Standardized axis typography / weights across all rows
        pax.spines["top"].set_visible(False)
        pax.spines["right"].set_visible(False)
        pax.spines["left"].set_color(AXIS_COLOR)
        pax.spines["bottom"].set_color(AXIS_COLOR)
        pax.spines["left"].set_linewidth(AXIS_LW)
        pax.spines["bottom"].set_linewidth(AXIS_LW)
        pax.tick_params(
            axis="both",
            labelsize=TICK_PT,
            length=1.8,
            width=0.35,
            colors=FAINT,
            pad=1.2,
        )
        y0, y1 = pax.get_ylim()
        if y0 < 0.0 < y1:
            pax.axhline(0.0, color=ZERO_LINE, lw=ZERO_LW, zorder=0)
        pax.yaxis.grid(True, color=GRID_COLOR, lw=0.35, zorder=0)
        pax.set_axisbelow(True)
        pax.set_facecolor("none")

        # Simplified aligned summary mini-table
        if row["plot_type"] == "scatter":
            header = "Spearman"
            pairs = [
                ("ρ", f"{row['rho']:.3f}"),
                ("nominal p", f"{row['p']:.4f}"),
                ("ROIs", f"{row['n']}"),
            ]
        else:
            delta_txt = f"−{abs(row['delta']):.3f}" if row["delta"] < 0 else f"{row['delta']:.3f}"
            header = "Welch contrast"
            pairs = [
                ("ΔAF−SR", delta_txt),
                ("nominal p", f"{row['p']:.4f}"),
                ("SR / AF n", f"{row['sr_n']} / {row['af_n']}"),
            ]
        _draw_summary_table(
            ax,
            summary_label_x,
            summary_value_x,
            row_center + 5.0,
            header,
            pairs,
        )

    # Bottom legend: Rhythm and Patient as separate dimensions
    legend_y = 10.2
    ax.text(
        MARGIN_MM,
        legend_y + 1.6,
        "Rhythm",
        ha="left",
        va="center",
        fontsize=NOTE_PT,
        color=MUTED,
        fontweight="semibold",
    )
    # Neutral open / filled circles for rhythm only
    ax.scatter(
        [MARGIN_MM + 12.0],
        [legend_y + 1.6],
        marker="o",
        s=18,
        facecolors=WHITE,
        edgecolors=FG,
        linewidths=SR_EDGE_LW,
        transform=ax.transData,
        zorder=3,
        clip_on=False,
    )
    ax.text(
        MARGIN_MM + 13.8,
        legend_y + 1.6,
        "SR ROI",
        ha="left",
        va="center",
        fontsize=NOTE_PT,
        color=MUTED,
    )
    ax.scatter(
        [MARGIN_MM + 28.0],
        [legend_y + 1.6],
        marker="o",
        s=18,
        facecolors=FG,
        edgecolors=FG,
        linewidths=AF_EDGE_LW,
        transform=ax.transData,
        zorder=3,
        clip_on=False,
    )
    ax.text(
        MARGIN_MM + 29.8,
        legend_y + 1.6,
        "AF ROI",
        ha="left",
        va="center",
        fontsize=NOTE_PT,
        color=MUTED,
    )

    ax.text(
        MARGIN_MM + 48.0,
        legend_y + 1.6,
        "Patient",
        ha="left",
        va="center",
        fontsize=NOTE_PT,
        color=MUTED,
        fontweight="semibold",
    )
    px0 = MARGIN_MM + 60.0
    for i, patient in enumerate(patients):
        px = px0 + i * 9.0
        color = PATIENT_COLORS[patient]
        marker = PATIENT_SHAPES[patient]
        ax.scatter(
            [px],
            [legend_y + 1.6],
            marker=marker,
            s=18,
            facecolors=color,
            edgecolors=darken(color),
            linewidths=0.7,
            transform=ax.transData,
            zorder=3,
            clip_on=False,
        )
        ax.text(
            px + 1.8,
            legend_y + 1.6,
            PATIENT_SHORT[patient],
            ha="left",
            va="center",
            fontsize=NOTE_PT,
            color=MUTED,
        )

    # Footer
    ax.text(
        MARGIN_MM,
        5.6,
        "Patient = color/shape; rhythm = open/filled; panel-specific axes.",
        ha="left",
        va="bottom",
        fontsize=FOOTER_PT,
        color=MUTED,
    )
    ax.text(
        MARGIN_MM,
        3.0,
        "ROI-level summaries; within-patient ROI clustering was not modeled.",
        ha="left",
        va="bottom",
        fontsize=FOOTER_PT,
        color=MUTED,
    )

    return fig, ax


def validate_layout(fig: plt.Figure, ax: plt.Axes) -> dict[str, float | int]:
    fig.canvas.draw()
    texts = []
    renderer = fig.canvas.get_renderer()
    for candidate_ax in [ax] + list(fig.axes):
        for artist in candidate_ax.get_children():
            if not isinstance(artist, mpl.text.Text):
                continue
            label = artist.get_text().strip()
            if not label:
                continue
            bbox = artist.get_window_extent(renderer).transformed(ax.transData.inverted())
            x0, y0, x1, y1 = float(bbox.x0), float(bbox.y0), float(bbox.x1), float(bbox.y1)
            require(x0 >= -0.5, f"text clipped left: {label!r}")
            require(x1 <= WIDTH_MM + 0.5, f"text clipped right: {label!r}")
            require(y0 >= -0.5, f"text clipped bottom: {label!r}")
            require(y1 <= HEIGHT_MM + 0.5, f"text clipped top: {label!r}")
            texts.append((label, (x0, y0, x1, y1), float(artist.get_fontsize()), candidate_ax is ax))

    min_font = min(fs for _, _, fs, _ in texts)
    require(min_font >= 5.5 - 0.05, f"text below 5.5 pt: {min_font:.2f}")

    # Overlap checks among main-canvas labels, and against plot-axis labels
    # (including ticks) so x-context labels cannot collide with neighboring axes.
    main_texts = [(lab, bb) for lab, bb, _, is_main in texts if is_main]
    axis_texts = [(lab, bb) for lab, bb, _, is_main in texts if not is_main]
    for i, (label_a, bbox_a) in enumerate(main_texts):
        for label_b, bbox_b in main_texts[i + 1 :]:
            if label_a == label_b:
                continue
            x0 = max(bbox_a[0], bbox_b[0])
            y0 = max(bbox_a[1], bbox_b[1])
            x1 = min(bbox_a[2], bbox_b[2])
            y1 = min(bbox_a[3], bbox_b[3])
            if x1 > x0 and y1 > y0:
                overlap = (x1 - x0) * (y1 - y0)
                require(
                    overlap < 0.15,
                    f"text overlap: {label_a!r} vs {label_b!r} ({overlap:.2f} mm²)",
                )
        for label_b, bbox_b in axis_texts:
            x0 = max(bbox_a[0], bbox_b[0])
            y0 = max(bbox_a[1], bbox_b[1])
            x1 = min(bbox_a[2], bbox_b[2])
            y1 = min(bbox_a[3], bbox_b[3])
            if x1 > x0 and y1 > y0:
                overlap = (x1 - x0) * (y1 - y0)
                require(
                    overlap < 0.15,
                    f"axis/label overlap: {label_a!r} vs {label_b!r} ({overlap:.2f} mm²)",
                )

    # Forbidden axis labels
    all_labels = " ".join(lab for lab, *_ in texts)
    for forbidden in ("F10 score", "F15 score", "F9 score", "Family score"):
        require(forbidden not in all_labels, f"forbidden axis label present: {forbidden}")

    return {"n_text_objects": len(texts), "min_font_pt": min_font}


def write_qc(
    hashes: dict[str, str], layout_stats: dict[str, float | int], input_hash: str
) -> None:
    OUT_QC.write_text(
        f"""# Figure 5D GeoMx ROI spatial context QC

Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}

Status: **PASS**

## Layout checks
- Canvas: {WIDTH_MM:.0f} × {HEIGHT_MM:.0f} mm
- Text objects: {layout_stats["n_text_objects"]}
- Minimum font: {layout_stats["min_font_pt"]:.1f} pt
- No text overlap, no canvas clipping (validated)
- No repeated vertical family-score axis labels
- Patient encoding: color + shape; rhythm encoding: open / filled
- Patient palette is local to Figure 5D only

## Inputs
- `{IN_GEO.relative_to(ROOT)}` — `{input_hash}`

## Outputs
{chr(10).join(f"- `{path}` — `{digest}`" for path, digest in sorted(hashes.items()))}
""",
        encoding="utf-8",
    )


def main() -> None:
    assert_no_forbidden_reads()
    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    df = load_and_assert()
    W, meta, rows, markers = build_sources(df)
    write_source(W, meta)
    write_summary(W, meta, rows)

    fig, ax = render(W, meta, rows, markers)
    layout_stats = validate_layout(fig, ax)

    pdf = OUT_BASE.with_suffix(".pdf")
    svg = OUT_BASE.with_suffix(".svg")
    png = OUT_BASE.with_suffix(".png")
    fig.savefig(pdf, bbox_inches=None, pad_inches=0)
    fig.savefig(svg, bbox_inches=None, pad_inches=0)
    fig.savefig(png, dpi=600, bbox_inches=None, pad_inches=0)
    fig.savefig(OUT_PREVIEW, dpi=150, bbox_inches=None, pad_inches=0)
    plt.close(fig)

    svg_text = svg.read_text(encoding="utf-8")
    require("GeoMx ROI spatial context" in svg_text, "title missing")
    require("Not a Path A gate" not in svg_text, "badge must be removed")
    require("descriptive context only" in svg_text, "subtitle boundary missing")
    require("y = GeoMx family-signature score" in svg_text, "y-context missing")
    require("x = fibroblast-panel score" in svg_text, "F10 x-context missing")
    require("x = macrophage-panel score" in svg_text, "F15 x-context missing")
    require("F10" in svg_text and "GAG disease" in svg_text and "TGF" in svg_text, "F10 label missing")
    require("F15" in svg_text and "NTRK2" in svg_text and "TrkB" in svg_text, "F15 label missing")
    require("F9" in svg_text and "HS-GAG" in svg_text and "GPCR" in svg_text, "F9 label missing")
    require("Spearman" in svg_text and "0.439" in svg_text, "F10 rho missing")
    require("Spearman" in svg_text and "0.431" in svg_text, "F15 rho missing")
    require("Welch contrast" in svg_text and "0.033" in svg_text, "F9 delta missing")
    require("SR ROI" in svg_text and "AF ROI" in svg_text, "rhythm legend missing")
    require("P1" in svg_text and "P4" in svg_text, "patient legend missing")
    require(
        "Patient = color/shape; rhythm = open/filled; panel-specific axes." in svg_text,
        "footer missing",
    )
    require("Family score" not in svg_text, "Family score label must be removed")
    require(
        "F10 score" not in svg_text
        and "F15 score" not in svg_text
        and "F9 score" not in svg_text,
        "family score axis labels must be removed",
    )

    from PIL import Image

    img = Image.open(png)
    w_mm = img.size[0] * 25.4 / 600.0
    h_mm = img.size[1] * 25.4 / 600.0
    require(abs(w_mm - WIDTH_MM) < 1.0, f"PNG width {w_mm:.2f} != {WIDTH_MM}")
    require(abs(h_mm - HEIGHT_MM) < 1.0, f"PNG height {h_mm:.2f} != {HEIGHT_MM}")

    hashes = {
        str(p.relative_to(ROOT)): sha256(p)
        for p in (pdf, svg, png, OUT_PREVIEW, OUT_SOURCE, OUT_SUMMARY)
    }
    write_qc(hashes, layout_stats, sha256(IN_GEO))

    print(f"Wrote {OUT_SOURCE}")
    print(f"Wrote {OUT_SUMMARY}")
    print(f"Wrote {OUT_QC}")
    print(f"Wrote {OUT_BASE}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW}")


if __name__ == "__main__":
    main()
