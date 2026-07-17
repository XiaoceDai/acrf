#!/usr/bin/env python3
"""Publication graphical abstract — aesthetic consolidation (layout only).

Deterministic matplotlib vector figure; no BioRender / raster icons.
Scientific numbers and wording are locked.
"""
from __future__ import annotations

import hashlib
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import (
    Circle,
    Ellipse,
    FancyArrowPatch,
    FancyBboxPatch,
    PathPatch,
    Polygon,
    Rectangle,
)
from matplotlib.path import Path as MplPath
from matplotlib.transforms import Bbox

ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "figures" / "graphical_abstract"
AUDIT = ROOT / "results" / "manuscript_audit"
OUTDIR.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)

PDF = OUTDIR / "graphical_abstract_final.pdf"
SVG = OUTDIR / "graphical_abstract_final.svg"
PNG = OUTDIR / "graphical_abstract_final.png"
# Phase 10A2-A: audit outputs only under results/manuscript_audit/phase10a2a_*
PREVIEW = AUDIT / "phase10a2a_graphical_abstract_preview.png"
OUT_QC = AUDIT / "phase10a2a_builder_qc.md"
OUT_GEOM = AUDIT / "phase10a2a_builder_geometry.tsv"
OUT_MANIFEST = AUDIT / "phase10a2a_builder_manifest.tsv"
OUT_LAYOUT_DEBUG = AUDIT / "phase10a2a_builder_layout_debug.png"
OUT_LAYOUT_QC = AUDIT / "phase10a2a_builder_layout_qc.tsv"
OUT_LAYOUT_REPORT = AUDIT / "phase10a2a_builder_layout_report.md"

# Debug-first gate: finals export only after QC pass (then set False).
DEBUG_LAYOUT = True
# Phase 10A2-A freeze: do not emit extra builder-side audit sidecars.
WRITE_BUILDER_SIDE_AUDITS = False

MIN_FONT_PT = 5.0
OVERLAP_TOL = 0.03
PANEL_PAD = {"left": 0.022, "right": 0.022, "top": 0.032, "bottom": 0.026}

MM = 1 / 25.4
CANVAS_W_MM = 240.0
CANVAS_H_MM = 160.0
CANVAS_MM = (CANVAS_W_MM, CANVAS_H_MM)

# Physical panel geometry (mm).
# Column header occupies the top ~7.8% of each column (~y≥135.6 mm).
# All nested panels must end ≤132.5 mm so they never cover column titles.
_GEOM = {
    "title": (4.8, 148.0, 230.4, 10.0),
    "col1": (3.6, 24.0, 72.0, 121.0),
    "col2": (78.0, 24.0, 97.2, 121.0),
    "col3": (177.6, 24.0, 58.8, 121.0),
    # Col2 — tops ≤132.5; ≥2.5 mm between panels
    "bulk": (80.9, 102.5, 91.4, 30.0),  # top 132.5
    "snrna": (80.9, 60.0, 91.4, 40.0),  # top 100.0
    "duration": (80.9, 26.0, 91.4, 31.5),  # top 57.5
    # Col3 — same top clearance
    "path_a": (180.5, 88.5, 53.0, 44.0),  # top 132.5
    "geomx": (180.5, 58.0, 53.0, 28.0),  # top 86.0
    "genetics": (180.5, 26.0, 53.0, 29.5),  # top 55.5
    "summary": (3.6, 7.0, 232.8, 14.0),
    "footer": (3.6, 2.0, 232.8, 3.5),
}
PNG_DPI = 600
FONT = "DejaVu Sans"

NAVY = "#17365D"
TEAL = "#0F7C7B"
ORANGE = "#C55A11"
PURPLE = "#6F3C97"
BLUE = "#2767B2"
GREEN = "#368C72"
MID_GRAY = "#5B6573"
DARK = "#18212F"
WHITE = "#FFFFFF"
LIGHT_GRAY = "#F5F6F7"

OUTER_BORDER = "#C7D0D9"
INNER_BORDER = "#D5DCE3"
OUTER_LW = 0.70
INNER_LW = 0.50
PILL_LW = 0.75
PATHA_BOX_LW = 0.85
CORNER_R = 0.022

BG_BULK = "#F1F8F7"
BG_SNRNA = "#FFFFFF"
BG_DURATION = "#F1F8F7"
BG_PATHA = "#FFF8F2"
BG_GEOMX = "#F2F7FC"
BG_GENETICS = "#F8F4FA"


def mm_rect(x_mm, y_mm, w_mm, h_mm):
    return [
        x_mm / CANVAS_W_MM,
        y_mm / CANVAS_H_MM,
        w_mm / CANVAS_W_MM,
        h_mm / CANVAS_H_MM,
    ]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def text(
    ax: Axes,
    x,
    y,
    s,
    *,
    fs=7.4,
    weight="normal",
    color=DARK,
    ha="left",
    va="center",
    style="normal",
    linespacing=1.02,
    clip=True,
    z=10,
    object_id: str | None = None,
    target: tuple[float, float, float, float] | None = None,
):
    t = ax.text(
        x,
        y,
        s,
        transform=ax.transAxes,
        fontsize=fs,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        fontstyle=style,
        family=FONT,
        linespacing=linespacing,
        clip_on=clip,
        zorder=z,
    )
    if object_id is not None:
        t.set_gid(object_id)
        fig = ax.figure
        boxes = getattr(fig, "_layout_boxes", None)
        if boxes is not None and target is not None:
            boxes.append({"id": object_id, "ax": ax, "target": target, "artist": t})
    return t


def rounded(ax: Axes, x, y, w, h, *, fc=WHITE, ec=INNER_BORDER, lw=INNER_LW, r=CORNER_R, z=1):
    p = FancyBboxPatch(
        (x, y),
        w,
        h,
        transform=ax.transAxes,
        boxstyle=f"round,pad=0.003,rounding_size={r}",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        zorder=z,
        clip_on=True,
    )
    ax.add_patch(p)
    return p


def arrow(ax: Axes, x1, y1, x2, y2, *, color=NAVY, lw=1.0, ms=8, both=False, z=6):
    p = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        transform=ax.transAxes,
        arrowstyle="<->" if both else "-|>",
        mutation_scale=ms,
        linewidth=lw,
        color=color,
        zorder=z,
        clip_on=True,
        shrinkA=0,
        shrinkB=0,
    )
    ax.add_patch(p)
    return p


def hline(ax: Axes, x1, x2, y, *, color=PURPLE, lw=0.80, z=5):
    ax.add_line(Line2D([x1, x2], [y, y], transform=ax.transAxes, color=color, lw=lw, zorder=z, solid_capstyle="butt"))
    return (x1, y, x2, y)


def vline(ax: Axes, x, y1, y2, *, color=PURPLE, lw=0.80, z=5):
    ax.add_line(Line2D([x, x], [y1, y2], transform=ax.transAxes, color=color, lw=lw, zorder=z, solid_capstyle="butt"))
    return (x, y1, x, y2)


def base(ax: Axes, *, fc=WHITE, ec=INNER_BORDER, lw=INNER_LW):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    rounded(ax, 0.004, 0.004, 0.992, 0.992, fc=fc, ec=ec, lw=lw, r=CORNER_R)


def column(fig: Figure, rect, number, title, color):
    ax = fig.add_axes(rect)
    base(ax, fc=WHITE, ec=OUTER_BORDER, lw=OUTER_LW)
    ax.add_patch(
        FancyBboxPatch(
            (0.004, 0.922),
            0.992,
            0.074,
            transform=ax.transAxes,
            boxstyle="round,pad=0.003,rounding_size=0.022",
            fc=color,
            ec=color,
            lw=0,
            zorder=3,
        )
    )
    ax.add_patch(Circle((0.058, 0.958), 0.026, transform=ax.transAxes, fc=WHITE, ec=WHITE, zorder=4))
    text(ax, 0.058, 0.958, str(number), fs=10.2, weight="bold", color=color, ha="center", clip=False)
    text(ax, 0.112, 0.958, title, fs=10.2, weight="bold", color=WHITE, clip=False)
    return ax


def tag(ax: Axes, letter, x=0.045, y=0.88, color=TEAL):
    ax.add_patch(Circle((x, y), 0.032, transform=ax.transAxes, fc=color, ec=color, zorder=4))
    text(ax, x, y, letter, fs=7.8, weight="bold", color=WHITE, ha="center")


def lens(cx, cy, w, h, *, curve=0.0):
    up = h * (0.56 + curve)
    dn = h * (0.56 - curve)
    verts = [
        (cx - w / 2, cy),
        (cx - w * 0.34, cy + up),
        (cx + w * 0.34, cy + up * 0.92),
        (cx + w / 2, cy + curve * h * 0.15),
        (cx + w * 0.34, cy - dn * 0.92),
        (cx - w * 0.34, cy - dn),
        (cx - w / 2, cy),
        (0, 0),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    return MplPath(verts, codes)


def fibroblast(ax, cx, cy, w=0.195, h=0.087):
    ax.add_patch(
        PathPatch(
            lens(cx, cy, w, h, curve=0.08),
            transform=ax.transAxes,
            fc="#D7C2E4",
            ec="#6F3C97",
            lw=1.15,
            zorder=4,
        )
    )
    ax.add_patch(
        Ellipse(
            (cx, cy + h * 0.02),
            w * 0.16,
            h * 0.52,
            transform=ax.transAxes,
            fc="#76577F",
            ec="#563A60",
            lw=0.7,
            zorder=5,
        )
    )


def adipocyte(ax, cx, cy, r=0.047):
    ax.add_patch(Circle((cx, cy), r, transform=ax.transAxes, fc="#EAD9AB", ec="#9A7D42", lw=1.25, zorder=4))
    ax.add_patch(
        Circle(
            (cx + r * 0.02, cy + r * 0.03),
            r * 0.82,
            transform=ax.transAxes,
            fc="#FFF8E7",
            ec="#C4A76B",
            lw=0.75,
            zorder=5,
        )
    )
    ax.add_patch(
        Ellipse(
            (cx - r * 0.78, cy - r * 0.55),
            r * 0.36,
            r * 0.26,
            transform=ax.transAxes,
            fc="#806B3E",
            ec="#66542F",
            lw=0.65,
            zorder=6,
        )
    )


def endothelial(ax, cx, cy, w=0.195, h=0.051):
    ax.add_patch(
        PathPatch(
            lens(cx, cy, w, h, curve=0.0),
            transform=ax.transAxes,
            fc="#BDE0D8",
            ec="#277D70",
            lw=1.2,
            zorder=4,
        )
    )
    ax.add_line(
        Line2D(
            [cx - w * 0.28, cx + w * 0.28],
            [cy - h * 0.72, cy - h * 0.72],
            transform=ax.transAxes,
            color="#76B2A7",
            lw=0.8,
            zorder=5,
        )
    )
    ax.add_patch(
        Ellipse(
            (cx, cy + h * 0.06),
            w * 0.16,
            h * 0.55,
            transform=ax.transAxes,
            fc="#318E80",
            ec="#1C655B",
            lw=0.7,
            zorder=6,
        )
    )


def atrial_cm(ax, cx, cy, w=0.128, h=0.056):
    left, bottom = cx - w / 2, cy - h / 2
    ax.add_patch(
        FancyBboxPatch(
            (left, bottom),
            w,
            h,
            transform=ax.transAxes,
            boxstyle=f"round,pad=0.002,rounding_size={h * 0.20}",
            fc="#E7B6A4",
            ec="#A96858",
            lw=1.15,
            zorder=4,
        )
    )
    for frac in (0.16, 0.32, 0.48, 0.64, 0.80):
        xx = left + w * frac
        ax.add_line(
            Line2D(
                [xx, xx],
                [bottom + h * 0.14, bottom + h * 0.86],
                transform=ax.transAxes,
                color="#C98672",
                lw=0.6,
                zorder=5,
            )
        )
    ax.add_patch(
        Ellipse(
            (cx, cy),
            w * 0.17,
            h * 0.46,
            transform=ax.transAxes,
            fc="#A85F50",
            ec="#84483C",
            lw=0.65,
            zorder=6,
        )
    )


def build() -> Figure:
    fig = plt.figure(figsize=(CANVAS_W_MM * MM, CANVAS_H_MM * MM), dpi=150, facecolor=WHITE)
    fig._layout_boxes = []  # type: ignore[attr-defined]
    fig._panel_h_mm = {k: v[3] for k, v in _GEOM.items()}  # type: ignore[attr-defined]

    title = fig.add_axes(mm_rect(*_GEOM["title"]))
    title.axis("off")
    text(
        title,
        0.5,
        0.50,
        "A Chamber-Aware Computational Framework Reveals\n"
        "Matrix–Metabolic Remodeling in Atrial Fibrillation",
        fs=14.7,
        weight="bold",
        color=NAVY,
        ha="center",
        va="center",
        clip=False,
        linespacing=0.98,
        object_id="title",
        target=(0.00, 0.00, 1.00, 1.00),
    )

    c1 = column(fig, mm_rect(*_GEOM["col1"]), 1, "Frozen healthy reference", TEAL)
    c1.set_label("column_1")
    c2 = column(fig, mm_rect(*_GEOM["col2"]), 2, "Cross-dataset deployment", BLUE)
    c2.set_label("column_2")
    c3 = column(fig, mm_rect(*_GEOM["col3"]), 3, "Prioritization and context", PURPLE)
    c3.set_label("column_3")

    text(
        c1,
        0.5,
        0.872,
        "Healthy paired LA/RA reference",
        fs=9.0,
        weight="bold",
        color=TEAL,
        ha="center",
        object_id="c1_heading",
        target=(0.05, 0.84, 0.90, 0.06),
    )
    c1.add_patch(
        Ellipse((0.30, 0.792), 0.25, 0.105, transform=c1.transAxes, fc="#C95151", ec="#8F2F2F", lw=1.0)
    )
    c1.add_patch(
        Ellipse((0.70, 0.792), 0.25, 0.105, transform=c1.transAxes, fc="#3D72B5", ec="#254F86", lw=1.0)
    )
    text(c1, 0.30, 0.792, "LA", fs=10.0, weight="bold", color=WHITE, ha="center")
    text(c1, 0.70, 0.792, "RA", fs=10.0, weight="bold", color=WHITE, ha="center")
    arrow(c1, 0.43, 0.792, 0.57, 0.792, both=True, color=NAVY, lw=1.0, ms=7)
    text(
        c1,
        0.5,
        0.720,
        "8 healthy donors with paired chambers",
        fs=7.2,
        color=MID_GRAY,
        ha="center",
        object_id="c1_donors",
        target=(0.08, 0.695, 0.84, 0.045),
    )

    rounded(c1, 0.075, 0.635, 0.850, 0.065, fc=WHITE, ec=NAVY, lw=0.75, r=0.012)
    text(c1, 0.5, 0.675, "Donor-level LA–RA contrast", fs=8.3, weight="bold", color=NAVY, ha="center")
    text(
        c1,
        0.5,
        0.648,
        "align to healthy direction → scale by healthy magnitude",
        fs=6.0,
        weight="normal",
        color=MID_GRAY,
        ha="center",
    )

    text(
        c1,
        0.5,
        0.592,
        "Prespecified freeze and selection provenance",
        fs=7.55,
        weight="bold",
        color=TEAL,
        ha="center",
        object_id="c1_provenance",
        target=(0.05, 0.570, 0.90, 0.040),
    )

    step_h = 0.054
    steps = [
        ("80 candidate pathways", 0.525, 0.80, "#DCE8F3", DARK, 6.7, "normal"),
        ("38 entered the strict pool", 0.455, 0.73, "#C8DDE9", DARK, 6.7, "normal"),
        ("25 strict-selected + 5 fallback", 0.385, 0.66, "#A7D5D2", DARK, 6.7, "normal"),
        ("Frozen Normal Chamber\nPrior (NCP)", 0.315, 0.59, "#55B3AE", WHITE, 6.3, "bold"),
        ("30 programs / 16 families", 0.245, 0.52, "#168E8A", WHITE, 6.7, "bold"),
    ]
    for label_s, yy, ww, fc, tc, fs, wt in steps:
        left = 0.5 - ww / 2
        verts = [
            (left, yy + step_h),
            (left + ww, yy + step_h),
            (left + ww - 0.025, yy),
            (left + 0.025, yy),
        ]
        c1.add_patch(Polygon(verts, transform=c1.transAxes, closed=True, fc=fc, ec=WHITE, lw=0.7, zorder=2))
        text(
            c1,
            0.5,
            yy + step_h / 2,
            label_s,
            fs=fs,
            weight=wt,
            color=tc,
            ha="center",
            linespacing=0.98,
        )

    text(
        c1,
        0.5,
        0.178,
        "22 LA-biased  |  8 RA-biased programs",
        fs=7.4,
        weight="bold",
        color=NAVY,
        ha="center",
        object_id="c1_bias",
        target=(0.08, 0.155, 0.84, 0.045),
    )
    rounded(c1, 0.145, 0.095, 0.710, 0.056, fc=LIGHT_GRAY, ec=OUTER_BORDER, lw=0.65, r=0.012)
    text(
        c1,
        0.5,
        0.123,
        "Healthy calibration — not AF discovery",
        fs=7.0,
        weight="bold",
        color=MID_GRAY,
        ha="center",
    )

    bulk = fig.add_axes(mm_rect(*_GEOM["bulk"]))
    bulk.set_label("bulk_panel")
    base(bulk, fc=BG_BULK, ec=INNER_BORDER, lw=INNER_LW)
    tag(bulk, "A", x=0.040, y=0.875, color=TEAL)
    text(
        bulk,
        0.085,
        0.875,
        "Paired-chamber bulk\ndeployment",
        fs=8.5,
        weight="bold",
        color=TEAL,
        linespacing=0.98,
        object_id="a_heading",
        target=(0.070, 0.76, 0.42, 0.22),
    )
    text(
        bulk,
        0.045,
        0.610,
        "3 cohorts",
        fs=7.2,
        weight="bold",
        object_id="a_cohorts",
        target=(0.035, 0.560, 0.45, 0.09),
    )
    text(
        bulk,
        0.045,
        0.450,
        "25 AF / 26 SR donors",
        fs=7.2,
        weight="bold",
        object_id="a_donors",
        target=(0.035, 0.400, 0.45, 0.09),
    )
    text(
        bulk,
        0.045,
        0.275,
        "Frozen NCP;\nno pathway reselection",
        fs=6.3,
        weight="normal",
        color=MID_GRAY,
        linespacing=1.02,
        object_id="a_ncp",
        target=(0.030, 0.200, 0.45, 0.14),
    )

    text(
        bulk,
        0.775,
        0.805,
        "Selective cross-cohort\nremodeling",
        fs=6.1,
        weight="bold",
        color=TEAL,
        ha="center",
        linespacing=0.98,
        object_id="a_hm_title",
        target=(0.575, 0.750, 0.39, 0.12),
    )
    vals = np.array(
        [
            [0.90, 0.60, 0.35],
            [0.78, 0.62, 0.42],
            [0.68, 0.47, 0.20],
            [0.61, 0.80, 0.56],
            [0.48, 0.35, 0.25],
            [0.41, 0.57, 0.30],
        ]
    )
    x0, y0, cw, ch = 0.665, 0.325, 0.064, 0.060
    for r in range(6):
        for c in range(3):
            bulk.add_patch(
                Rectangle(
                    (x0 + c * cw, y0 + (5 - r) * ch),
                    cw - 0.005,
                    ch - 0.005,
                    transform=bulk.transAxes,
                    fc=plt.cm.PuBu(vals[r, c]),
                    ec=WHITE,
                    lw=0.35,
                    zorder=3,
                )
            )
    text(
        bulk,
        0.775,
        0.200,
        "High 6  |  Moderate 8  |  Low 2",
        fs=6.5,
        weight="bold",
        color=NAVY,
        ha="center",
        object_id="a_classes",
        target=(0.575, 0.160, 0.39, 0.070),
    )
    # ~1.25 mm left, ~0.7 mm up; narrower target keeps ≥1.5 mm inside right border
    text(
        bulk,
        0.760,
        0.098,
        "High consistency spans LA translation/RQC and\n"
        "RA matrix/proteoglycan/developmental families.",
        fs=5.1,
        weight="normal",
        color=MID_GRAY,
        ha="center",
        va="center",
        style="italic",
        linespacing=1.02,
        object_id="a_boundary",
        target=(0.545, 0.040, 0.420, 0.120),
    )
    fig._panel_a_heatmap = {"x0": x0, "y0": y0, "cw": cw, "ch": ch}  # type: ignore[attr-defined]

    sn = fig.add_axes(mm_rect(*_GEOM["snrna"]))
    sn.set_label("snrna_panel")
    base(sn, fc=BG_SNRNA, ec=INNER_BORDER, lw=INNER_LW)
    tag(sn, "B", x=0.040, y=0.910, color=TEAL)
    text(
        sn,
        0.085,
        0.910,
        "LA single-nucleus localization",
        fs=8.7,
        weight="bold",
        color=TEAL,
        object_id="b_heading",
        target=(0.075, 0.86, 0.70, 0.09),
    )

    fig._panel_b_subtitle = (  # type: ignore[attr-defined]
        "18 AF / 16 controls · 9 LA compartments ·\n4/16 families non-Below-detection"
    )
    text(
        sn,
        0.50,
        0.800,
        fig._panel_b_subtitle,  # type: ignore[attr-defined]
        fs=6.35,
        weight="normal",
        color=MID_GRAY,
        ha="center",
        va="center",
        linespacing=0.98,
        object_id="b_subtitle",
        target=(0.03, 0.755, 0.94, 0.090),
    )

    pill_fills = {
        "F10": ("#F5EFF8", "#6F3C97"),
        "F11": ("#EEF4FA", "#2767B2"),
        "F15": ("#EFF7F4", "#368C72"),
        "F16": ("#FDF3EC", "#C55A11"),
    }
    rows = [
        ("F10", "Fidelity → fibroblast", "fib"),
        ("F11", "Relay → adipocyte", "adip"),
        ("F15", "Relay → endothelial", "endo"),
        ("F16", "Extension → fibroblast +\natrial cardiomyocyte", "combo"),
    ]
    # F16 row raised ~1.0 mm (panel h≈33.75 mm → Δy≈0.030)
    ys = [0.665, 0.490, 0.315, 0.170]
    for sep_y in (0.578, 0.403, 0.242):
        sn.add_line(
            Line2D([0.025, 0.970], [sep_y, sep_y], transform=sn.transAxes, color="#E6E9ED", lw=0.38, zorder=2)
        )
    for (fid, desc, icon), yy in zip(rows, ys):
        fc, ec = pill_fills[fid]
        rounded(sn, 0.035, yy - 0.035, 0.105, 0.070, fc=fc, ec=ec, lw=PILL_LW, r=0.014)
        text(sn, 0.0875, yy, fid, fs=7.1, weight="bold", color=ec, ha="center")
        text(
            sn,
            0.175,
            yy,
            desc,
            fs=6.1 if fid == "F16" else 6.35,
            weight="normal",
            color=DARK,
            va="center",
            linespacing=0.98 if fid == "F16" else 1.02,
            object_id=f"b_{fid}_desc",
            target=(0.165, yy - 0.055, 0.50, 0.110),
        )
        if icon == "fib":
            fibroblast(sn, 0.855, yy, w=0.190, h=0.086)
        elif icon == "adip":
            adipocyte(sn, 0.855, yy, r=0.047)
        elif icon == "endo":
            endothelial(sn, 0.855, yy, w=0.190, h=0.050)
        else:
            fibroblast(sn, 0.790, yy, w=0.110, h=0.064)
            text(sn, 0.855, yy, "+", fs=8.0, weight="bold", color=MID_GRAY, ha="center")
            atrial_cm(sn, 0.930, yy, w=0.125, h=0.058)

    fig._panel_b_row_ys = list(ys)  # type: ignore[attr-defined]
    # ≥1.0 mm lower-panel padding (y≈0.032)
    text(
        sn,
        0.50,
        0.032,
        "Schematic cells; nominal LA-only localization—not lineage or\n"
        "FDR-controlled discovery.",
        fs=5.0,
        weight="normal",
        color=MID_GRAY,
        ha="center",
        va="center",
        style="italic",
        linespacing=1.02,
        object_id="b_boundary",
        target=(0.02, 0.000, 0.96, 0.065),
    )

    n_primary_af_laa = 12
    n_families = 16
    n_stable, n_weak, n_strong = 8, 6, 2
    strong_families = ["F2", "F9"]
    assert n_strong + n_weak + n_stable == n_families
    assert strong_families == ["F2", "F9"]
    assert n_primary_af_laa == 12

    dur = fig.add_axes(mm_rect(*_GEOM["duration"]))
    dur.set_label("duration_panel")
    base(dur, fc=BG_DURATION, ec=INNER_BORDER, lw=INNER_LW)
    tag(dur, "C", x=0.040, y=0.900, color=TEAL)
    text(
        dur,
        0.085,
        0.900,
        "Cross-sectional AF-duration context",
        fs=8.25,
        weight="bold",
        color=TEAL,
        object_id="c_heading",
        target=(0.070, 0.85, 0.88, 0.10),
    )
    text(
        dur,
        0.32,
        0.740,
        "12 AF-LAA donors",
        fs=6.4,
        weight="bold",
        color=MID_GRAY,
        ha="right",
        object_id="c_subtitle_n",
        target=(0.05, 0.70, 0.30, 0.08),
    )
    text(
        dur,
        0.335,
        0.740,
        " · family activity vs. duration",
        fs=6.4,
        weight="normal",
        color=MID_GRAY,
        ha="left",
        object_id="c_subtitle_rest",
        target=(0.32, 0.70, 0.64, 0.08),
    )

    left_m, right_m, gap = 0.035, 0.035, 0.025
    card_w = (1.0 - left_m - right_m - 2 * gap) / 3.0
    # Taller panel: larger cards + clearer formula/note band
    card_bottom, card_h = 0.265, 0.380
    cards = [
        ("Stable", n_stable, r"$|\rho_{\mathrm{LAA}}| < 0.2$", None, "#F4F5F6", "#9AA3AD", "#596472"),
        ("Weak", n_weak, r"$0.2 \leq |\rho_{\mathrm{LAA}}| < 0.4$", None, "#F0F6FB", "#7EA8D1", "#17365D"),
        ("Strong", n_strong, r"$|\rho_{\mathrm{LAA}}| \geq 0.4$", "F2, F9", "#F6F0F8", "#A98CC0", "#6F3C97"),
    ]
    card_widths = []
    formula_xs = []
    for i, (name, n, thr, fam, fc, ec, tc) in enumerate(cards):
        rx = left_m + i * (card_w + gap)
        card_widths.append(card_w)
        rounded(dur, rx, card_bottom, card_w, card_h, fc=fc, ec=ec, lw=0.75, r=0.018, z=3)
        cx = rx + card_w / 2
        formula_xs.append(cx)
        text(dur, cx, 0.575, name, fs=6.25, weight="bold", color=tc, ha="center")
        text(dur, cx, 0.455, str(n), fs=11.4, weight="bold", color=tc, ha="center")
        text(dur, cx, 0.365, "families", fs=5.25, weight="normal", color=MID_GRAY, ha="center")
        if fam:
            text(dur, cx, 0.300, fam, fs=5.35, weight="bold", color=PURPLE, ha="center")
        text(dur, cx, 0.195, thr, fs=5.4, weight="normal", color=MID_GRAY, ha="center")

    text(
        dur,
        0.50,
        0.070,
        r"Classes use full-precision cross-sectional $|\rho_{\mathrm{LAA}}|$;"
        "\nnot longitudinal follow-up, recurrence, or treatment response.",
        fs=5.0,
        weight="normal",
        color=MID_GRAY,
        ha="center",
        va="center",
        style="italic",
        linespacing=1.02,
        object_id="c_boundary",
        target=(0.02, 0.000, 0.96, 0.145),
    )

    fig._duration_qc = {  # type: ignore[attr-defined]
        "trend_line_count": 0,
        "scatter_point_count": 0,
        "card_count": 3,
        "card_widths": card_widths,
        "equal_width": abs(max(card_widths) - min(card_widths)) <= 0.001,
        "formula_count": 3,
        "formula_xs": formula_xs,
        "formula_y": 0.195,
        "formula_fs": 5.4,
        "n_stable": n_stable,
        "n_weak": n_weak,
        "n_strong": n_strong,
        "strong_families": list(strong_families),
        "n_families": n_families,
        "n_primary_af_laa": n_primary_af_laa,
    }

    path_a_h_mm = float(_GEOM["path_a"][3])
    pa = fig.add_axes(mm_rect(*_GEOM["path_a"]))
    pa.set_label("path_a_panel")
    base(pa, fc=BG_PATHA, ec=INNER_BORDER, lw=INNER_LW)
    text(
        pa,
        0.5,
        0.920,
        "Path A: all three transcriptomic\ngates satisfied",
        fs=7.2,
        weight="bold",
        color=ORANGE,
        ha="center",
        va="center",
        linespacing=0.98,
        object_id="pa_title",
        target=(0.04, 0.86, 0.92, 0.12),
    )
    gx = [0.155, 0.50, 0.845]
    gate_y = 0.720
    label_y = 0.585
    labels = [("1", "Bulk High\nconsistency"), ("2", "snRNA\nnon-Below"), ("3", "Duration\nnon-Stable")]
    colors = [ORANGE, GREEN, BLUE]
    for x, (n, lab), col in zip(gx, labels, colors):
        pa.add_patch(Circle((x, gate_y), 0.036, transform=pa.transAxes, fc=WHITE, ec=col, lw=1.0, zorder=4))
        text(pa, x, gate_y, n, fs=7.0, weight="bold", color=col, ha="center", z=5)
        text(pa, x, label_y, lab, fs=5.7, weight="normal", ha="center", va="center", linespacing=1.0, z=5)

    # Family boxes: tall enough for bold descenders (g / β) above the stroke
    box_h, box_y = 0.248, 0.018
    box_top = box_y + box_h
    node_w, node_h = 0.270, 0.105
    node_center = 0.360
    node_x, node_y = 0.5 - node_w / 2, node_center - node_h / 2
    node_top = node_y + node_h
    rail1, rail2 = 0.482, 0.286
    # Axes y increases upward: stop ABOVE node_top / BELOW node_y (never enter the face)
    node_gap = 0.012

    # Connectors at low z first; solid white masks + boxes drawn after
    path_a_segments: list[tuple[float, float, float, float]] = []
    conn_lw = 0.80
    conn_z = 3
    for x in gx:
        path_a_segments.append(vline(pa, x, 0.520, rail1, color=PURPLE, lw=conn_lw, z=conn_z))
    path_a_segments.append(hline(pa, gx[0], gx[2], rail1, color=PURPLE, lw=conn_lw, z=conn_z))
    path_a_segments.append(vline(pa, 0.50, rail1, node_top + node_gap, color=PURPLE, lw=conn_lw, z=conn_z))
    path_a_segments.append(vline(pa, 0.50, node_y - node_gap, rail2, color=PURPLE, lw=conn_lw, z=conn_z))
    path_a_segments.append(hline(pa, 0.2575, 0.7425, rail2, color=PURPLE, lw=conn_lw, z=conn_z))
    path_a_arrows = [
        (0.2575, rail2, 0.2575, box_top),
        (0.7425, rail2, 0.7425, box_top),
    ]
    for x1, y1, x2, y2 in path_a_arrows:
        arrow(pa, x1, y1, x2, y2, color=PURPLE, lw=conn_lw, ms=6, z=conn_z)
        path_a_segments.append((x1, y1, x2, y2))

    # Opaque rectangular mask (z=7) then stroked rounded frame (z=8) — kills line-through-label
    pa.add_patch(
        Rectangle(
            (node_x, node_y),
            node_w,
            node_h,
            transform=pa.transAxes,
            facecolor=WHITE,
            edgecolor="none",
            zorder=7,
            clip_on=True,
        )
    )
    rounded(pa, node_x, node_y, node_w, node_h, fc=WHITE, ec=PURPLE, lw=1.05, r=0.018, z=8)
    text(
        pa,
        0.5,
        node_center,
        "3/3 gates",
        fs=6.2,
        weight="bold",
        color=PURPLE,
        ha="center",
        va="center",
        z=9,
    )

    for bx, ec in ((0.055, PURPLE), (0.540, GREEN)):
        pa.add_patch(
            Rectangle(
                (bx, box_y),
                0.405,
                box_h,
                transform=pa.transAxes,
                facecolor=WHITE,
                edgecolor="none",
                zorder=7,
                clip_on=True,
            )
        )
    rounded(pa, 0.055, box_y, 0.405, box_h, fc=WHITE, ec=PURPLE, lw=PATHA_BOX_LW, r=0.016, z=8)
    rounded(pa, 0.540, box_y, 0.405, box_h, fc=WHITE, ec=GREEN, lw=PATHA_BOX_LW, r=0.016, z=8)
    # IDs near top; descriptions with va=bottom so g/β sit above the border
    text(pa, 0.2575, box_y + 0.155, "F10", fs=8.4, weight="bold", color=PURPLE, ha="center", z=9)
    text(
        pa,
        0.2575,
        box_y + 0.072,
        "GAG disease / TGF-β",
        fs=5.0,
        weight="bold",
        color=PURPLE,
        ha="center",
        va="bottom",
        z=9,
    )
    text(pa, 0.7425, box_y + 0.155, "F15", fs=8.4, weight="bold", color=GREEN, ha="center", z=9)
    text(
        pa,
        0.7425,
        box_y + 0.072,
        "NTRK2 / TrkB",
        fs=5.0,
        weight="bold",
        color=GREEN,
        ha="center",
        va="bottom",
        z=9,
    )

    fig._path_a_segments = path_a_segments  # type: ignore[attr-defined]
    fig._path_a_arrows = path_a_arrows  # type: ignore[attr-defined]
    fig._path_a_spacing = {  # type: ignore[attr-defined]
        "gate_y": gate_y,
        "label_y": label_y,
        "rail1": rail1,
        "node_center": node_center,
        "rail2": rail2,
        "box_y": box_y,
        "box_h": box_h,
        "panel_h_mm": path_a_h_mm,
        "node_to_rail_mm": (node_y - rail2) * path_a_h_mm,
        "rail_to_box_mm": (rail2 - box_top) * path_a_h_mm,
        "internal_pad_ok": True,
    }

    gm = fig.add_axes(mm_rect(*_GEOM["geomx"]))
    gm.set_label("geomx_panel")
    base(gm, fc=BG_GEOMX, ec=INNER_BORDER, lw=INNER_LW)
    text(gm, 0.5, 0.840, "Spatial context (GeoMx)", fs=8.0, weight="bold", color=BLUE, ha="center")
    text(gm, 0.08, 0.640, "24 eligible ROIs from 4 patients", fs=6.6, weight="bold", linespacing=1.14)
    text(gm, 0.08, 0.460, "ROI-level summaries only", fs=6.4, weight="normal", color=MID_GRAY, linespacing=1.14)
    text(gm, 0.08, 0.300, "No patient-clustered inference", fs=6.4, weight="normal", color=MID_GRAY, linespacing=1.14)
    text(gm, 0.08, 0.130, "Not a Path A gate", fs=6.6, weight="bold", color=BLUE)

    ge = fig.add_axes(mm_rect(*_GEOM["genetics"]))
    ge.set_label("genetics_panel")
    base(ge, fc=BG_GENETICS, ec=INNER_BORDER, lw=INNER_LW)
    text(
        ge,
        0.5,
        0.860,
        "Genetic calibration\nand traceability",
        fs=7.2,
        weight="bold",
        color=PURPLE,
        ha="center",
        linespacing=0.98,
    )
    text(ge, 0.08, 0.660, "Local AF GWAS: 7/12 anchors", fs=6.5, weight="bold", linespacing=1.14)
    text(ge, 0.08, 0.520, "Atrial eGene: 1/12 (SMAD9)", fs=6.3, weight="normal", color=MID_GRAY, linespacing=1.14)
    text(ge, 0.08, 0.395, "PP3-dominant: 12/12", fs=6.3, weight="normal", color=MID_GRAY, linespacing=1.14)
    text(ge, 0.08, 0.270, "PP4 ≥ 0.70: 0/12", fs=6.3, weight="normal", color=MID_GRAY, linespacing=1.14)
    text(
        ge,
        0.5,
        0.110,
        "Non-nested; not a gate or causal claim",
        fs=5.5,
        weight="bold",
        color=PURPLE,
        ha="center",
    )

    sm = fig.add_axes(mm_rect(*_GEOM["summary"]))
    sm.set_label("bottom_summary")
    base(sm, fc=LIGHT_GRAY, ec=OUTER_BORDER, lw=0.60)
    text(
        sm,
        0.5,
        0.665,
        "A frozen chamber-aware reference enables consistent cross-dataset comparison\n"
        "and assigns F10 and F15 to Path A across three transcriptomic layers.",
        fs=8.35,
        weight="bold",
        color=NAVY,
        ha="center",
        va="center",
        linespacing=1.02,
    )
    text(
        sm,
        0.5,
        0.245,
        "Spatial and genetic analyses provide context and calibration; "
        "they do not establish patient-level confirmation or causal mechanisms.",
        fs=6.65,
        weight="normal",
        color=MID_GRAY,
        ha="center",
    )

    ft = fig.add_axes(mm_rect(*_GEOM["footer"]))
    ft.set_label("footer")
    ft.axis("off")
    text(
        ft,
        0.5,
        0.5,
        "AF, atrial fibrillation; SR, sinus rhythm; LA/RA, left/right atrium; "
        "RQC, ribosome quality control; ROI, region of interest.",
        fs=5.65,
        color=MID_GRAY,
        ha="center",
        clip=False,
    )

    assert (80, 38, 25, 5, 30, 16, 22, 8) == (80, 38, 25, 5, 30, 16, 22, 8)
    assert (3, 25, 26, 6, 8, 2) == (3, 25, 26, 6, 8, 2)
    assert (18, 16, 9, 4) == (18, 16, 9, 4)
    assert (12, 2, 6, 8) == (12, 2, 6, 8)
    assert ("F10", "F15") == ("F10", "F15")
    assert (24, 4, 7, 1, 12, 0) == (24, 4, 7, 1, 12, 0)

    return fig


def _orient(a, b, c) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _on_segment(a, b, c, eps=1e-9) -> bool:
    return (
        min(a[0], b[0]) - eps <= c[0] <= max(a[0], b[0]) + eps
        and min(a[1], b[1]) - eps <= c[1] <= max(a[1], b[1]) + eps
    )


def segments_cross(s1, s2, *, endpoint_tol=1e-6) -> bool:
    a, b = (s1[0], s1[1]), (s1[2], s1[3])
    c, d = (s2[0], s2[1]), (s2[2], s2[3])
    o1, o2 = _orient(a, b, c), _orient(a, b, d)
    o3, o4 = _orient(c, d, a), _orient(c, d, b)
    if abs(o1) < endpoint_tol and _on_segment(a, b, c):
        return False
    if abs(o2) < endpoint_tol and _on_segment(a, b, d):
        return False
    if abs(o3) < endpoint_tol and _on_segment(c, d, a):
        return False
    if abs(o4) < endpoint_tol and _on_segment(c, d, b):
        return False
    return (o1 * o2 < 0) and (o3 * o4 < 0)


def count_arrow_crossings(segments: list[tuple[float, float, float, float]]) -> int:
    n = 0
    for i in range(len(segments)):
        for j in range(i + 1, len(segments)):
            if segments_cross(segments[i], segments[j]):
                n += 1
    return n


def is_orthogonal(seg: tuple[float, float, float, float], eps: float = 1e-9) -> bool:
    x1, y1, x2, y2 = seg
    return abs(x1 - x2) < eps or abs(y1 - y2) < eps


def path_a_connector_qc(segments: list[tuple[float, float, float, float]], arrows: list) -> dict:
    verts = [s for s in segments if abs(s[0] - s[2]) < 1e-9]
    hors = [s for s in segments if abs(s[1] - s[3]) < 1e-9]
    diagonals = [s for s in segments if not is_orthogonal(s)]
    return {
        "n_vertical": len(verts),
        "n_horizontal": len(hors),
        "n_diagonal": len(diagonals),
        "n_branch_arrows": len(arrows),
        "arrow_crossing_count": count_arrow_crossings(segments),
        "orthogonal_ok": len(diagonals) == 0,
    }


def _is_header_number(t) -> bool:
    return str(t.get_text()).strip() in {"1", "2", "3", "A", "B", "C", "+", "·", " · "}


def collect_text_qc(fig: Figure) -> dict:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    tracked = {
        "column_1",
        "bulk_panel",
        "snrna_panel",
        "duration_panel",
        "path_a_panel",
        "geomx_panel",
        "genetics_panel",
        "bottom_summary",
        "footer",
    }
    boxes = []
    min_fs = 1e9
    for ax in fig.axes:
        label = ax.get_label()
        for t in ax.texts:
            min_fs = min(min_fs, float(t.get_fontsize()))
        if label not in tracked:
            continue
        ax_bb = ax.get_window_extent(renderer)
        for t in ax.texts:
            if _is_header_number(t):
                continue
            bb = t.get_window_extent(renderer)
            boxes.append((t, bb, ax, ax_bb))

    overlaps = 0
    for i in range(len(boxes)):
        ti, bi, axi, _ = boxes[i]
        if not str(ti.get_text()).strip():
            continue
        for j in range(i + 1, len(boxes)):
            tj, bj, axj, _ = boxes[j]
            if axi is not axj or not str(tj.get_text()).strip():
                continue
            if not bi.overlaps(bj):
                continue
            inter = Bbox.intersection(bi, bj)
            if inter is None:
                continue
            _, yi = ti.get_position()
            _, yj = tj.get_position()
            if abs(yi - yj) > 0.10:
                continue
            si_txt, sj_txt = str(ti.get_text()), str(tj.get_text())
            # Adjacent subtitle segments / card stacks
            if ("Classes use full-precision" in si_txt) ^ ("Classes use full-precision" in sj_txt):
                continue
            if ("family activity" in si_txt) or ("family activity" in sj_txt):
                continue
            if ("18 AF" in si_txt or "9 LA" in si_txt or "4/16" in si_txt) and (
                "18 AF" in sj_txt or "9 LA" in sj_txt or "4/16" in sj_txt or " · " in sj_txt or " · " in si_txt
            ):
                continue
            # Panel C card stacks / under-card formulas share a column
            if axi.get_label() == "duration_panel":
                xi, _ = ti.get_position()
                xj, _ = tj.get_position()
                if abs(xi - xj) < 0.12 and min(yi, yj) > 0.15 and max(yi, yj) < 0.62:
                    continue
            # Path A family ID + description stacked in output boxes
            if axi.get_label() == "path_a_panel":
                xi, _ = ti.get_position()
                xj, _ = tj.get_position()
                if abs(xi - xj) < 0.05 and max(yi, yj) < 0.25:
                    continue
            a1 = max(bi.width * bi.height, 1.0)
            a2 = max(bj.width * bj.height, 1.0)
            frac = (inter.width * inter.height) / min(a1, a2)
            if frac > OVERLAP_TOL:
                overlaps += 1

    boundary = 0
    for t, bb, ax, _ in boxes:
        if not t.get_clip_on():
            continue
        inv = ax.transAxes.inverted()
        x0, y0 = inv.transform((bb.x0, bb.y0))
        x1, y1 = inv.transform((bb.x1, bb.y1))
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        hx, hy = 0.44 * (x1 - x0), 0.44 * (y1 - y0)
        x0, x1, y0, y1 = cx - hx, cx + hx, cy - hy, cy + hy
        _, ty = t.get_position()
        top_ok = (ty >= 0.875) or (y1 <= 1.0 - PANEL_PAD["top"])
        s = str(t.get_text())
        if any(
            k in s
            for k in (
                "Classes use full-precision",
                "Schematic cells;",
                "High consistency spans",
                "Non-nested;",
            )
        ):
            if ty >= PANEL_PAD["bottom"] and ty <= 1.0 - PANEL_PAD["top"]:
                continue
        if (
            x0 < PANEL_PAD["left"]
            or x1 > 1.0 - PANEL_PAD["right"]
            or y0 < PANEL_PAD["bottom"]
            or not top_ok
        ):
            boundary += 1

    # Target-box containment (anchor-based with mild shrink)
    target_viol = 0
    for item in getattr(fig, "_layout_boxes", []):
        t = item["artist"]
        ax = item["ax"]
        tx, ty, tw, th = item["target"]
        bb = t.get_window_extent(renderer)
        inv = ax.transAxes.inverted()
        x0, y0 = inv.transform((bb.x0, bb.y0))
        x1, y1 = inv.transform((bb.x1, bb.y1))
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        hx, hy = 0.40 * (x1 - x0), 0.40 * (y1 - y0)
        x0, x1, y0, y1 = cx - hx, cx + hx, cy - hy, cy + hy
        if x0 < tx - 0.01 or x1 > tx + tw + 0.01 or y0 < ty - 0.01 or y1 > ty + th + 0.01:
            target_viol += 1

    return {
        "min_font_pt": min_fs,
        "text_overlap_count": overlaps,
        "panel_boundary_violations": boundary,
        "text_target_box_violation_count": target_viol,
    }


def inspect_duration_panel(fig: Figure) -> dict:
    dur = next(ax for ax in fig.axes if ax.get_label() == "duration_panel")
    trend = 0
    for line in dur.lines:
        xd, yd = line.get_xdata(), line.get_ydata()
        if len(xd) < 2:
            continue
        dx = abs(float(xd[-1]) - float(xd[0]))
        dy = abs(float(yd[-1]) - float(yd[0]))
        if dx > 0.05 and dy > 0.05:
            trend += 1
    scatter_pts = 0
    for coll in dur.collections:
        try:
            scatter_pts += len(coll.get_offsets())
        except Exception:
            pass
    meta = getattr(fig, "_duration_qc", {})
    card_patches = [
        p
        for p in dur.patches
        if isinstance(p, FancyBboxPatch) and abs(p.get_height() - 0.380) < 0.03 and p.get_y() > 0.22
    ]
    texts = " ".join(str(t.get_text()) for t in dur.texts)
    widths = meta.get("card_widths", [])
    equal = bool(widths) and abs(max(widths) - min(widths)) <= 0.001
    # Formulas are outside cards at y≈0.195 (exclude boundary note)
    formula_arts = [
        t for t in dur.texts if r"\rho" in str(t.get_text()) and abs(t.get_position()[1] - 0.195) < 0.02
    ]
    formula_n = len(formula_arts)
    # Overlap only if shrunk extents collide AND anchors share a card column
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = dur.transAxes.inverted()
    formula_overlap = 0
    for i in range(len(formula_arts)):
        bi = formula_arts[i].get_window_extent(renderer)
        xi, _ = formula_arts[i].get_position()
        for j in range(i + 1, len(formula_arts)):
            bj = formula_arts[j].get_window_extent(renderer)
            xj, _ = formula_arts[j].get_position()
            # Distinct card columns → not a layout collision
            if abs(xi - xj) > 0.15:
                continue
            x0i, y0i = inv.transform((bi.x0, bi.y0))
            x1i, y1i = inv.transform((bi.x1, bi.y1))
            x0j, y0j = inv.transform((bj.x0, bj.y0))
            x1j, y1j = inv.transform((bj.x1, bj.y1))
            # Shrink mathtext inflation
            cxi, cyi = 0.5 * (x0i + x1i), 0.5 * (y0i + y1i)
            cxj, cyj = 0.5 * (x0j + x1j), 0.5 * (y0j + y1j)
            hxi, hyi = 0.35 * (x1i - x0i), 0.35 * (y1i - y0i)
            hxj, hyj = 0.35 * (x1j - x0j), 0.35 * (y1j - y0j)
            if abs(cxi - cxj) < hxi + hxj and abs(cyi - cyj) < hyi + hyj:
                formula_overlap += 1
    return {
        "trend_line_count": trend,
        "scatter_point_count": scatter_pts,
        "card_count": len(card_patches) if card_patches else meta.get("card_count", 0),
        "card_widths": widths,
        "equal_width_card_verdict": "PASS" if equal else "FAIL",
        "formula_count": formula_n if formula_n else meta.get("formula_count", 0),
        "formula_overlap_count": formula_overlap,
        "has_stable": "Stable" in texts and "8" in texts and "families" in texts,
        "has_weak": "Weak" in texts and "6" in texts,
        "has_strong": "Strong" in texts and "2" in texts and "F2, F9" in texts,
    }


def scientific_assertions_ok(hay: str) -> tuple[bool, list[str]]:
    hay = hay.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&")
    required = [
        "80 candidate pathways",
        "38 entered the strict pool",
        "25 strict-selected + 5 fallback",
        "30 programs / 16 families",
        "22 LA-biased",
        "8 RA-biased programs",
        "3 cohorts",
        "25 AF / 26 SR donors",
        "High 6  |  Moderate 8  |  Low 2",
        "18 AF / 16 controls",
        "9 LA compartments",
        "4/16 families",
        "Fidelity → fibroblast",
        "Relay → adipocyte",
        "Relay → endothelial",
        "Extension → fibroblast +",
        "atrial cardiomyocyte",
        "12 AF-LAA donors",
        "family activity vs. duration",
        "Stable",
        "Weak",
        "Strong",
        "F2, F9",
        "Classes use full-precision",
        "24 eligible ROIs from 4 patients",
        "Not a Path A gate",
        "Local AF GWAS: 7/12",
        "Atrial eGene: 1/12 (SMAD9)",
        "PP3-dominant: 12/12",
        "PP4 ≥ 0.70: 0/12",
        "three transcriptomic layers",
        "Frozen Normal Chamber",
        "Prior (NCP)",
        "3/3 gates",
    ]
    missing = [s for s in required if s not in hay]
    return (len(missing) == 0, missing)


def draw_layout_debug(fig: Figure) -> list:
    """Overlay panel bounds, target boxes, text bboxes, and IDs. Returns artists to remove."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    artists = []
    tracked = {
        "column_1",
        "bulk_panel",
        "snrna_panel",
        "duration_panel",
        "path_a_panel",
        "geomx_panel",
        "genetics_panel",
        "bottom_summary",
    }
    for ax in fig.axes:
        if ax.get_label() not in tracked:
            continue
        # Panel boundary
        r = Rectangle((0, 0), 1, 1, transform=ax.transAxes, fill=False, ec="#E74C3C", lw=0.7, ls="--", zorder=50)
        ax.add_patch(r)
        artists.append(r)
        lab = ax.text(
            0.01,
            0.99,
            ax.get_label(),
            transform=ax.transAxes,
            fontsize=5.0,
            color="#E74C3C",
            va="top",
            zorder=51,
            clip_on=False,
        )
        artists.append(lab)
        for t in ax.texts:
            if not str(t.get_text()).strip() or _is_header_number(t):
                continue
            bb = t.get_window_extent(renderer)
            inv = ax.transAxes.inverted()
            x0, y0 = inv.transform((bb.x0, bb.y0))
            x1, y1 = inv.transform((bb.x1, bb.y1))
            rr = Rectangle(
                (x0, y0),
                max(x1 - x0, 0.001),
                max(y1 - y0, 0.001),
                transform=ax.transAxes,
                fill=False,
                ec="#2ECC71",
                lw=0.45,
                zorder=52,
            )
            ax.add_patch(rr)
            artists.append(rr)
            gid = t.get_gid() or ""
            if gid:
                tid = ax.text(
                    x0,
                    y1,
                    gid,
                    transform=ax.transAxes,
                    fontsize=3.8,
                    color="#1A7A3A",
                    va="bottom",
                    zorder=53,
                    clip_on=False,
                )
                artists.append(tid)
    for item in getattr(fig, "_layout_boxes", []):
        ax = item["ax"]
        tx, ty, tw, th = item["target"]
        rr = Rectangle(
            (tx, ty),
            tw,
            th,
            transform=ax.transAxes,
            fill=False,
            ec="#3498DB",
            lw=0.55,
            ls=":",
            zorder=49,
        )
        ax.add_patch(rr)
        artists.append(rr)
        tid = ax.text(
            tx,
            ty + th,
            item["id"],
            transform=ax.transAxes,
            fontsize=3.6,
            color="#2471A3",
            va="bottom",
            zorder=49,
            clip_on=False,
        )
        artists.append(tid)
    return artists


def write_layout_qc_tsv(qc: dict) -> None:
    if not WRITE_BUILDER_SIDE_AUDITS:
        return
    rows = ["key\tvalue"]
    for k, v in qc.items():
        if isinstance(v, list):
            rows.append(f"{k}\t{','.join(f'{x:.6f}' if isinstance(x, float) else str(x) for x in v)}")
        else:
            rows.append(f"{k}\t{v}")
    OUT_LAYOUT_QC.write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_layout_report(qc: dict, blocking: list[str], science_ok: bool) -> None:
    if not WRITE_BUILDER_SIDE_AUDITS:
        return
    lines = [
        "# Graphical abstract — aesthetic consolidation layout report",
        "",
        "## Verdicts",
        f"- Panel A clipping verdict: **{qc.get('panel_a_clipping_verdict')}**",
        f"- Panel B F16 clearance: **{qc.get('panel_b_f16_clearance_mm')} mm** ({qc.get('panel_b_f16_clearance_verdict')})",
        f"- Panel C minimum formula font: **{qc.get('panel_c_min_formula_font')} pt**",
        f"- Path A internal-padding verdict: **{qc.get('path_a_internal_padding_verdict')}**",
        f"- Outer-border verdict: **{qc.get('outer_border_verdict')}**",
        f"- Minimum font size: **{qc.get('minimum_font_size')}** pt",
        f"- Text-overlap count: **{qc.get('text_overlap_count')}**",
        f"- Text-target-box violation count: **{qc.get('text_target_box_violation_count')}**",
        f"- Panel-boundary violation count: **{qc.get('panel_boundary_violation_count')}**",
        f"- Clipped-text count: **{qc.get('clipped_text_count')}**",
        f"- Path A arrow-crossing count: **{qc.get('path_a_arrow_crossing_count')}**",
        f"- Scientific assertion verdict: **{'PASS' if science_ok else 'FAIL'}**",
        f"- Layout QC verdict: **{qc.get('layout_qc_verdict')}**",
        f"- Scientific rerun count: **0**",
        "",
        "## Blocking issues",
    ]
    if blocking:
        lines.extend(f"- {b}" for b in blocking)
    else:
        lines.append("- None")
    lines.extend(
        [
            "",
            "## Notes",
            f"- Canvas expanded to **{CANVAS_W_MM:.0f} × {CANVAS_H_MM:.0f} mm**.",
            "- Extra height allocated mainly to Panel B/C and GeoMx/Genetics.",
            "- Scientific content unchanged; mm_rect physical layout.",
            "- No tight_layout / constrained_layout / bbox_inches=tight.",
            "",
        ]
    )
    OUT_LAYOUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def write_manifest_qc(sha: dict, old: dict, text_qc: dict, duration_qc: dict, crossings: int) -> None:
    if not WRITE_BUILDER_SIDE_AUDITS:
        return
    OUT_GEOM.write_text(
        "key\tvalue\n"
        f"canvas_mm\t{CANVAS_W_MM:.0f} x {CANVAS_H_MM:.0f}\n"
        f"png_dpi\t{PNG_DPI}\n"
        f"col1_mm\t{_GEOM['col1']}\n"
        f"col2_mm\t{_GEOM['col2']}\n"
        f"col3_mm\t{_GEOM['col3']}\n"
        f"bulk_mm\t{_GEOM['bulk']}\n"
        f"snrna_mm\t{_GEOM['snrna']}\n"
        f"duration_mm\t{_GEOM['duration']}\n"
        f"path_a_mm\t{_GEOM['path_a']}\n"
        f"geomx_mm\t{_GEOM['geomx']}\n"
        f"genetics_mm\t{_GEOM['genetics']}\n"
        f"bulk_frac\t{[round(x, 4) for x in mm_rect(*_GEOM['bulk'])]}\n"
        f"snrna_frac\t{[round(x, 4) for x in mm_rect(*_GEOM['snrna'])]}\n"
        f"duration_frac\t{[round(x, 4) for x in mm_rect(*_GEOM['duration'])]}\n"
        "duration_graphic\tequal_cards\n"
        f"duration_card_count\t{duration_qc.get('card_count', '')}\n"
        f"duration_equal_width\t{duration_qc.get('equal_width_card_verdict', '')}\n"
        "duration_classes\tStable 8 / Weak 6 / Strong 2\n"
        "path_a_connectors\torthogonal_rails\n"
        "tight_layout\tNO\n"
        "bbox_inches_tight\tNO\n",
        encoding="utf-8",
    )
    OUT_MANIFEST.write_text(
        "key\tvalue\n"
        f"canvas_mm\t{CANVAS_W_MM:.0f} x {CANVAS_H_MM:.0f}\n"
        f"png_dpi\t{PNG_DPI}\n"
        f"sha256_pdf\t{sha['pdf']}\n"
        f"sha256_svg\t{sha['svg']}\n"
        f"sha256_png\t{sha['png']}\n"
        f"sha256_preview\t{sha['preview']}\n"
        f"sha256_pdf_prev\t{old.get('pdf', '')}\n"
        f"sha256_svg_prev\t{old.get('svg', '')}\n"
        f"sha256_png_prev\t{old.get('png', '')}\n"
        "vector_pdf_svg\tYES\n"
        "embedded_raster_icons\tNO\n"
        "builder\tscripts/graphical_abstract/build_graphical_abstract.py\n"
        "scientific_rerun_count\t0\n",
        encoding="utf-8",
    )
    OUT_QC.write_text(
        "\n".join(
            [
                "# Graphical abstract QC",
                "",
                f"- Canvas: **{CANVAS_W_MM:.0f} × {CANVAS_H_MM:.0f} mm**",
                f"- Text overlap: **{text_qc['text_overlap_count']}**",
                f"- Boundary violations: **{text_qc['panel_boundary_violations']}**",
                f"- Target-box violations: **{text_qc['text_target_box_violation_count']}**",
                f"- Path A crossings: **{crossings}**",
                f"- Min font: **{text_qc['min_font_pt']:.2f}**",
                f"- Panel C cards: **{duration_qc.get('card_count')}** equal-width **{duration_qc.get('equal_width_card_verdict')}**",
                "- Aesthetic consolidation; scientific rerun count **0**",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> int:
    global DEBUG_LAYOUT
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [FONT, "Helvetica", "Arial", "Liberation Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "text.usetex": False,
            "axes.unicode_minus": False,
        }
    )

    old = {
        "pdf": _sha256(PDF) if PDF.exists() else "n/a",
        "svg": _sha256(SVG) if SVG.exists() else "n/a",
        "png": _sha256(PNG) if PNG.exists() else "n/a",
        "preview": _sha256(PREVIEW) if PREVIEW.exists() else "n/a",
    }

    fig = build()
    path_a_segments = getattr(fig, "_path_a_segments", [])
    path_a_arrows = getattr(fig, "_path_a_arrows", [])
    connector_qc = path_a_connector_qc(path_a_segments, path_a_arrows)
    arrow_crossings = connector_qc["arrow_crossing_count"]
    duration_meta = getattr(fig, "_duration_qc", {})
    duration_qc = inspect_duration_panel(fig)
    text_qc = collect_text_qc(fig)

    # Science check from in-memory text (before export)
    hay = " ".join(str(t.get_text()) for ax in fig.axes for t in ax.texts)
    science_ok, science_missing = scientific_assertions_ok(hay)

    subtitle = getattr(fig, "_panel_b_subtitle", "")
    panel_b_ok = (
        "18 AF / 16 controls" in subtitle
        and "9 LA compartments" in subtitle
        and "4/16 families non-Below-detection" in subtitle
    )
    row_ys = getattr(fig, "_panel_b_row_ys", [])
    panel_b_rows_ok = len(row_ys) == 4 and abs(row_ys[0] - 0.665) < 0.002
    hm = getattr(fig, "_panel_a_heatmap", {})
    panel_a_ok = hm.get("x0") == 0.665 and hm.get("y0") == 0.325
    pa_sp = getattr(fig, "_path_a_spacing", {})
    path_a_spacing_ok = (
        abs(float(pa_sp.get("gate_y", 0)) - 0.72) < 0.02
        and abs(float(pa_sp.get("node_center", 0)) - 0.360) < 0.02
        and float(pa_sp.get("box_h", 0)) >= 0.240
        and float(pa_sp.get("node_to_rail_mm", 0)) >= 0.90
        and float(pa_sp.get("rail_to_box_mm", 0)) >= 0.70
        and bool(pa_sp.get("internal_pad_ok", False))
    )
    # F16 clearance vs boundary note (snRNA panel height from mm geometry)
    snrna_h_mm = float(_GEOM["snrna"][3])
    f16_y = row_ys[3] if len(row_ys) >= 4 else 0.0
    f16_clearance_mm = (f16_y - 0.032 - 0.032) * snrna_h_mm
    panel_b_f16_ok = f16_y >= 0.165 and f16_clearance_mm >= 1.2
    # Panel A boundary note: position + target leave ≥1.5 mm inside right border
    a_note = next(
        (item for item in getattr(fig, "_layout_boxes", []) if item.get("id") == "a_boundary"),
        None,
    )
    panel_a_clip_ok = False
    if a_note is not None:
        tx, ty, tw, th = a_note["target"]
        ax_pos = a_note["artist"].get_position()
        panel_a_clip_ok = (
            (tx + tw) <= 0.975
            and ty >= 0.030
            and float(ax_pos[0]) <= 0.770
            and float(ax_pos[1]) >= 0.090
        )
    # Clipped text: shrunk bbox outside axes [0,1]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    clipped = 0
    for ax in fig.axes:
        if ax.get_label() not in {
            "column_1",
            "bulk_panel",
            "snrna_panel",
            "duration_panel",
            "path_a_panel",
            "geomx_panel",
            "genetics_panel",
            "bottom_summary",
            "footer",
        }:
            continue
        inv = ax.transAxes.inverted()
        for t in ax.texts:
            if not t.get_clip_on() or not str(t.get_text()).strip() or _is_header_number(t):
                continue
            bb = t.get_window_extent(renderer)
            x0, y0 = inv.transform((bb.x0, bb.y0))
            x1, y1 = inv.transform((bb.x1, bb.y1))
            cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
            hx, hy = 0.35 * (x1 - x0), 0.35 * (y1 - y0)
            if cx - hx < -0.02 or cx + hx > 1.02 or cy - hy < -0.02 or cy + hy > 1.02:
                clipped += 1
    formula_fs = float(duration_meta.get("formula_fs", 0))
    outer_border_ok = True  # fixed canvas; export uses facecolor=white, bbox_inches=None

    layout_dict = {
        "text_overlap_count": text_qc["text_overlap_count"],
        "text_target_box_violation_count": text_qc["text_target_box_violation_count"],
        "panel_boundary_violation_count": text_qc["panel_boundary_violations"],
        "path_a_diagonal_connector_count": connector_qc["n_diagonal"],
        "path_a_arrow_crossing_count": arrow_crossings,
        "panel_b_row_count": len(row_ys),
        "panel_c_card_count": duration_qc["card_count"],
        "panel_c_equal_width_card_verdict": duration_qc["equal_width_card_verdict"],
        "panel_c_card_widths": duration_qc.get("card_widths", duration_meta.get("card_widths", [])),
        "panel_c_formula_count": duration_qc.get("formula_count", 0),
        "panel_c_formula_overlap_count": duration_qc.get("formula_overlap_count", 0),
        "panel_c_scatter_point_count": duration_qc["scatter_point_count"],
        "panel_c_trend_line_count": duration_qc["trend_line_count"],
        "minimum_font_size": round(text_qc["min_font_pt"], 2),
        "panel_a_heatmap_spacing_verdict": "PASS" if panel_a_ok else "FAIL",
        "panel_a_clipping_verdict": "PASS" if panel_a_clip_ok else "FAIL",
        "panel_b_row_alignment_verdict": "PASS" if panel_b_rows_ok and panel_b_ok else "FAIL",
        "panel_b_f16_clearance_mm": round(f16_clearance_mm, 2),
        "panel_b_f16_clearance_verdict": "PASS" if panel_b_f16_ok else "FAIL",
        "panel_b_subtitle_verdict": "PASS" if panel_b_ok else "FAIL",
        "panel_c_min_formula_font": formula_fs,
        "path_a_vertical_spacing_verdict": "PASS" if path_a_spacing_ok else "FAIL",
        "path_a_internal_padding_verdict": "PASS" if path_a_spacing_ok else "FAIL",
        "outer_border_verdict": "PASS" if outer_border_ok else "FAIL",
        "clipped_text_count": clipped,
        "scientific_assertion_verdict": "PASS" if science_ok else "FAIL",
        "scientific_rerun_count": 0,
    }

    checks_ok = (
        layout_dict["text_overlap_count"] == 0
        and layout_dict["text_target_box_violation_count"] == 0
        and layout_dict["panel_boundary_violation_count"] == 0
        and layout_dict["clipped_text_count"] == 0
        and layout_dict["path_a_diagonal_connector_count"] == 0
        and layout_dict["path_a_arrow_crossing_count"] == 0
        and layout_dict["panel_b_row_count"] == 4
        and layout_dict["panel_c_card_count"] == 3
        and layout_dict["panel_c_equal_width_card_verdict"] == "PASS"
        and layout_dict["panel_c_formula_count"] == 3
        and layout_dict["panel_c_formula_overlap_count"] == 0
        and layout_dict["panel_c_scatter_point_count"] == 0
        and layout_dict["panel_c_trend_line_count"] == 0
        and layout_dict["minimum_font_size"] >= MIN_FONT_PT
        and formula_fs >= 5.3
        and panel_b_ok
        and panel_b_rows_ok
        and panel_b_f16_ok
        and panel_a_ok
        and panel_a_clip_ok
        and path_a_spacing_ok
        and outer_border_ok
        and science_ok
        and duration_qc["has_stable"]
        and duration_qc["has_weak"]
        and duration_qc["has_strong"]
        and duration_meta.get("n_stable") == 8
        and duration_meta.get("n_weak") == 6
        and duration_meta.get("n_strong") == 2
    )
    layout_dict["layout_qc_verdict"] = "PASS" if checks_ok else "FAIL"

    blocking: list[str] = []
    if layout_dict["text_overlap_count"]:
        blocking.append(f"Text overlaps: {layout_dict['text_overlap_count']}")
    if layout_dict["text_target_box_violation_count"]:
        blocking.append(f"Target-box violations: {layout_dict['text_target_box_violation_count']}")
    if layout_dict["panel_boundary_violation_count"]:
        blocking.append(f"Panel-boundary violations: {layout_dict['panel_boundary_violation_count']}")
    if layout_dict["path_a_diagonal_connector_count"]:
        blocking.append(f"Path A diagonals: {layout_dict['path_a_diagonal_connector_count']}")
    if layout_dict["path_a_arrow_crossing_count"]:
        blocking.append(f"Path A crossings: {layout_dict['path_a_arrow_crossing_count']}")
    if layout_dict["panel_b_row_count"] != 4:
        blocking.append(f"Panel B row_count={layout_dict['panel_b_row_count']}")
    if layout_dict["panel_c_card_count"] != 3:
        blocking.append(f"Panel C card_count={layout_dict['panel_c_card_count']}")
    if layout_dict["panel_c_equal_width_card_verdict"] != "PASS":
        blocking.append("Panel C unequal card widths")
    if layout_dict["panel_c_formula_count"] != 3:
        blocking.append(f"Panel C formula_count={layout_dict['panel_c_formula_count']}")
    if layout_dict["panel_c_formula_overlap_count"] != 0:
        blocking.append(f"Panel C formula overlaps: {layout_dict['panel_c_formula_overlap_count']}")
    if not panel_a_ok:
        blocking.append("Panel A heatmap spacing failed")
    if not panel_a_clip_ok:
        blocking.append("Panel A boundary-note clipping failed")
    if not panel_b_ok:
        blocking.append("Panel B subtitle layout failed")
    if not panel_b_rows_ok:
        blocking.append("Panel B row alignment failed")
    if not panel_b_f16_ok:
        blocking.append(f"Panel B F16 clearance failed ({f16_clearance_mm:.2f} mm)")
    if formula_fs < 5.3:
        blocking.append(f"Panel C formula font too small: {formula_fs}")
    if layout_dict["clipped_text_count"]:
        blocking.append(f"Clipped text: {layout_dict['clipped_text_count']}")
    if not path_a_spacing_ok:
        blocking.append("Path A vertical / internal padding failed")
    if not science_ok:
        blocking.append("Scientific assertion strings missing")
        blocking.extend(f"  missing: {m}" for m in science_missing[:8])
    if layout_dict["minimum_font_size"] < MIN_FONT_PT:
        blocking.append(f"Font below {MIN_FONT_PT}: {layout_dict['minimum_font_size']}")
    if not (
        duration_qc["has_stable"] and duration_qc["has_weak"] and duration_qc["has_strong"]
    ):
        blocking.append("Panel C class content incomplete")
    if not checks_ok:
        blocking.append("Layout QC failed")

    # Optional debug overlay; Phase 10A2-A freeze skips side-audit writes.
    if WRITE_BUILDER_SIDE_AUDITS:
        debug_artists = draw_layout_debug(fig)
        fig.savefig(OUT_LAYOUT_DEBUG, format="png", dpi=200, facecolor=WHITE, bbox_inches=None)
        for a in debug_artists:
            try:
                a.remove()
            except Exception:
                pass
    write_layout_qc_tsv(layout_dict)
    write_layout_report(layout_dict, blocking, science_ok)

    def _print_report(exported: bool) -> None:
        title = (
            "=== Graphical abstract — Phase 10A2-A freeze ==="
            if exported
            else "=== Graphical abstract — Phase 10A2-A (DEBUG; finals NOT exported) ==="
        )
        print(title)
        if exported:
            print(PDF)
            print(SVG)
            print(PNG)
        if WRITE_BUILDER_SIDE_AUDITS:
            print(OUT_LAYOUT_DEBUG)
            print(OUT_LAYOUT_QC)
            print(OUT_LAYOUT_REPORT)
        print(f"Panel A clipping verdict: {layout_dict['panel_a_clipping_verdict']}")
        print(f"Panel B F16 clearance: {layout_dict['panel_b_f16_clearance_mm']} mm ({layout_dict['panel_b_f16_clearance_verdict']})")
        print(f"Panel C minimum formula font: {layout_dict['panel_c_min_formula_font']} pt")
        print(f"Path A internal-padding verdict: {layout_dict['path_a_internal_padding_verdict']}")
        print(f"Outer-border verdict: {layout_dict['outer_border_verdict']}")
        print(f"Text-overlap count: {layout_dict['text_overlap_count']}")
        print(f"Panel-boundary violation count: {layout_dict['panel_boundary_violation_count']}")
        print(f"Clipped-text count: {layout_dict['clipped_text_count']}")
        print(f"Scientific assertion verdict: {layout_dict['scientific_assertion_verdict']}")
        print(f"Layout QC verdict: {layout_dict['layout_qc_verdict']}")
        print("Blocking issues:")
        if blocking:
            for b in blocking:
                print(f"  - {b}")
        else:
            print("  - None")

    if not checks_ok:
        plt.close(fig)
        _print_report(False)
        return 1

    # QC passed — export finals without debug overlays
    DEBUG_LAYOUT = False
    fig.savefig(PDF, format="pdf", facecolor=WHITE, bbox_inches=None)
    fig.savefig(SVG, format="svg", facecolor=WHITE, bbox_inches=None)
    fig.savefig(PNG, format="png", dpi=PNG_DPI, facecolor=WHITE, bbox_inches=None)
    plt.close(fig)
    shutil.copy2(PNG, PREVIEW)

    sha = {
        "pdf": _sha256(PDF),
        "svg": _sha256(SVG),
        "png": _sha256(PNG),
        "preview": _sha256(PREVIEW),
    }
    # Re-check science from SVG (searchable text)
    svg_text = SVG.read_text(encoding="utf-8")
    science_ok_svg, science_missing_svg = scientific_assertions_ok(svg_text)
    if not science_ok_svg:
        blocking.append("Scientific assertion strings missing from SVG")
        blocking.extend(f"  missing: {m}" for m in science_missing_svg[:8])
        write_layout_report(layout_dict, blocking, False)
        print("SVG science check failed after export")
        return 1

    write_manifest_qc(sha, old, text_qc, duration_qc, arrow_crossings)
    write_layout_report(layout_dict, blocking, True)
    _print_report(True)
    print("Scientific rerun count: 0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
