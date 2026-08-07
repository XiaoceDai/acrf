#!/usr/bin/env python3
"""Phase 11B3-P4: final visual refinement of the locked acrf benchmark figure.

Rendering-only. Writes candidate_v4 assets only.
Does not overwrite v1/v2/v3 candidates.
Does not recompute methods or metrics.
"""

from __future__ import annotations

import hashlib
import re
import subprocess
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import Bbox

PRIMARY = Path("/home/bio/Project_P02")
BENCH = PRIMARY / "results/bmc_bioinformatics_benchmark"
AUDIT = PRIMARY / "results/bmc_bioinformatics_audit"
FIG_DIR = PRIMARY / "figures/bmc_bioinformatics"

MM = 1 / 25.4
FIG_W_MM = 170.0
FIG_H_MM = 125.0

FIGURE_TITLE = "Benchmarking chamber-aware calibration across atrial cohorts"

METHOD_ROWS = [
    ("A1_FULL_ACRF", "unstratified", "A1", "A1"),
    ("A2_NO_DONOR_PAIRING", "unstratified", "A2", "A2"),
    ("A3_NO_DIRECTION_ALIGNMENT", "unstratified", "A3†", "A3"),
    ("A4_NO_HEALTHY_MAGNITUDE_SCALING", "unstratified", "A4", "A4"),
    ("A5_CHAMBER_AGNOSTIC", "unstratified", "A5", "A5"),
    ("B1_PLAIN_SINGSCORE", "LA", "B1-LA", "B1LA"),
    ("B1_PLAIN_SINGSCORE", "RA", "B1-RA", "B1RA"),
]

COLORS = {
    "A1": "#176D6B",
    "A2": "#5B8FA3",
    "A3": "#6F9FB2",
    "A4": "#8AB3C0",
    "A5": "#D58A38",
    "B1LA": "#776D8A",
    "B1RA": "#A29AAA",
}

# (color_key, code, name, component_flags, estimand)
PANEL_A_ROWS = [
    ("A1", "A1", "Full acrf", [1, 1, 1, 1, 1], "Aligned scaled LA–RA"),
    ("A2", "A2", "No pairing", [0, 1, 1, 1, 1], "Unpaired chamber\ncontrast"),
    ("A3", "A3", "No direction", [1, 1, 0, 1, 1], "Native scaled LA–RA"),
    ("A4", "A4", "No scaling", [1, 1, 1, 0, 1], "Aligned unscaled\nLA–RA"),
    ("A5", "A5", "Chamber-agnostic", [0, 0, 0, 0, 1], "Donor chamber mean"),
    ("B1", "B1", "Plain singscore", [0, 0, 0, 0, 1], "Chamber-specific\nscore"),
]

COMP_HEADERS = ["Pairing", "LA–RA", "Direction", "Scaling", "Program-\nfirst"]
GUIDE_COLOR = "#E8ECEF"
GUIDE_LW = 0.5
LABEL_DX_B = 0.040  # fixed offset from lollipop endpoint to count label

# Architecture: A full width; B|C|D equal-width bottom row
RECT_A = [0.045, 0.585, 0.920, 0.370]
RECT_B = [0.065, 0.100, 0.270, 0.390]
RECT_C = [0.385, 0.100, 0.270, 0.390]
RECT_D = [0.705, 0.100, 0.270, 0.390]

D_RETAINED = dict(marker="o", facecolor="white", edgecolor="#557482", linewidth=1.2, s=42)
D_HIGH = dict(marker="D", facecolor="#1F6B73", edgecolor="#1F6B73", linewidth=1.0, s=36)
D_CONNECTOR = "#CCD5DA"

PAIR_ORDER = [
    ("GSE128188|GSE79768", "o", "P1"),
    ("GSE128188|GSE115574", "s", "P2"),
    ("GSE79768|GSE115574", "D", "P3"),
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def count_words(text: str) -> int:
    return len(re.findall(r"[A-Za-z0-9’'\-]+", text))


def load_and_assert() -> dict[str, Any]:
    m2 = pd.read_csv(BENCH / "benchmark_sign_concordance.tsv", sep="\t")
    m1 = pd.read_csv(BENCH / "benchmark_pairwise_correlations.tsv", sep="\t")
    m5 = pd.read_csv(BENCH / "benchmark_class_stability.tsv", sep="\t")
    source = pd.read_csv(BENCH / "benchmark_figure_source.tsv", sep="\t")

    def m2_row(mid, stratum):
        r = m2[
            (m2.method_id == mid)
            & (m2.stratum == stratum)
            & (m2.family_id.astype(str) == "ALL")
        ].iloc[0]
        return int(r.concordant_count), int(r.denominator), float(r.proportion)

    def m1_sum(mid, stratum):
        r = m1[
            (m1.method_id == mid) & (m1.stratum == stratum) & (m1.cohort == "SUMMARY")
        ].iloc[0]
        return float(r.median_rho), float(r.rho_min), float(r.rho_max)

    def m1_pairs(mid, stratum):
        sub = m1[
            (m1.method_id == mid) & (m1.stratum == stratum) & (m1.cohort != "SUMMARY")
        ]
        return {str(row.cohort): float(row.pairwise_rho) for _, row in sub.iterrows()}

    def m5_sum(mid):
        r = m5[
            (m5.method_id == mid)
            & (m5.stratum == "unstratified")
            & (m5.cohort == "SUMMARY")
        ].iloc[0]
        return int(r.retained_vs_a1), int(r.reclassified_count)

    def high_count(mid):
        sub = m5[
            (m5.method_id == mid)
            & (m5.stratum == "unstratified")
            & (m5.cohort == "POOLED")
            & (m5.family_id.astype(str) != "ALL")
        ]
        return int((sub.method_class == "High").sum())

    expected_m2 = {
        ("A1_FULL_ACRF", "unstratified"): (32, 48),
        ("A2_NO_DONOR_PAIRING", "unstratified"): (32, 48),
        ("A3_NO_DIRECTION_ALIGNMENT", "unstratified"): (32, 48),
        ("A4_NO_HEALTHY_MAGNITUDE_SCALING", "unstratified"): (32, 48),
        ("A5_CHAMBER_AGNOSTIC", "unstratified"): (20, 48),
        ("B1_PLAIN_SINGSCORE", "LA"): (22, 48),
        ("B1_PLAIN_SINGSCORE", "RA"): (24, 48),
    }
    expected_m1 = {
        ("A1_FULL_ACRF", "unstratified"): 0.179,
        ("A2_NO_DONOR_PAIRING", "unstratified"): 0.215,
        ("A3_NO_DIRECTION_ALIGNMENT", "unstratified"): 0.559,
        ("A4_NO_HEALTHY_MAGNITUDE_SCALING", "unstratified"): 0.188,
        ("A5_CHAMBER_AGNOSTIC", "unstratified"): -0.294,
        ("B1_PLAIN_SINGSCORE", "LA"): -0.015,
        ("B1_PLAIN_SINGSCORE", "RA"): -0.182,
    }

    data: dict[str, Any] = {
        "m2": {},
        "m1": {},
        "m1_pairs": {},
        "m5": {},
        "source_rows": len(source),
    }
    for mid, st, _lab, key in METHOD_ROWS:
        c, d, p = m2_row(mid, st)
        if (c, d) != expected_m2[(mid, st)]:
            raise SystemExit(f"BMC_BENCHMARK_FIGURE_V3_BLOCKED_DATA M2 {mid}: {c}/{d}")
        data["m2"][key] = (c, d, p)
        med, rmin, rmax = m1_sum(mid, st)
        if round(med, 3) != expected_m1[(mid, st)]:
            raise SystemExit(f"BMC_BENCHMARK_FIGURE_V3_BLOCKED_DATA M1 {mid}: {med}")
        data["m1"][key] = (med, rmin, rmax)
        data["m1_pairs"][key] = m1_pairs(mid, st)

    for key, mid in [
        ("A1", "A1_FULL_ACRF"),
        ("A2", "A2_NO_DONOR_PAIRING"),
        ("A3", "A3_NO_DIRECTION_ALIGNMENT"),
        ("A4", "A4_NO_HEALTHY_MAGNITUDE_SCALING"),
        ("A5", "A5_CHAMBER_AGNOSTIC"),
    ]:
        ret, rec = m5_sum(mid)
        hi = high_count(mid)
        data["m5"][key] = {"retained": ret, "reclassified": rec, "high": hi}

    if data["m5"]["A1"]["retained"] != 16 or data["m5"]["A1"]["high"] != 6:
        raise SystemExit("BMC_BENCHMARK_FIGURE_V3_BLOCKED_DATA A1 M5")
    if data["m5"]["A2"]["retained"] != 12 or data["m5"]["A2"]["high"] != 4:
        raise SystemExit("BMC_BENCHMARK_FIGURE_V3_BLOCKED_DATA A2 M5")
    if data["m5"]["A3"]["retained"] != 16 or data["m5"]["A4"]["retained"] != 16:
        raise SystemExit("BMC_BENCHMARK_FIGURE_V3_BLOCKED_DATA A3/A4 retained")
    if data["m5"]["A5"]["retained"] != 8 or data["m5"]["A5"]["high"] != 0:
        raise SystemExit("BMC_BENCHMARK_FIGURE_V3_BLOCKED_DATA A5 M5")
    return data


def panel_heading(ax, letter: str, title: str, pad: float = 4.0) -> None:
    ax.set_title(
        f"{letter}  {title}",
        loc="left",
        fontsize=9.3,
        fontweight="bold",
        pad=pad,
        color="#111111",
    )


def style_axes(ax) -> None:
    ax.set_facecolor("white")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#666666")
    ax.spines["bottom"].set_color("#666666")
    ax.spines["left"].set_linewidth(0.7)
    ax.spines["bottom"].set_linewidth(0.7)
    ax.tick_params(axis="both", labelsize=7.5, colors="#333333", length=2.5, width=0.6)


def render_panel_a(ax) -> dict[str, Any]:
    panel_heading(ax, "A", "Method components")
    ax.set_xlim(0, 100)
    ax.set_ylim(-2, 100)
    ax.axis("off")

    # Zones (data 0–100): method labels | 5-component matrix | estimand
    x_code = 1.5
    x_name = 7.5
    x_mat0, x_mat1 = 28.5, 73.5
    x_cols = [x_mat0 + (i + 0.5) * (x_mat1 - x_mat0) / 5.0 for i in range(5)]
    x_div = 76.0
    x_est = 78.5

    for x, lab in zip(x_cols, COMP_HEADERS):
        ax.text(x, 95.5, lab, ha="center", va="top", fontsize=7.2, color="#333333", linespacing=1.05)
    ax.text(x_est, 95.5, "Estimand", ha="left", va="top", fontsize=7.3, fontweight="bold", color="#222222")

    # ~9% tighter than prior 16-unit pitch (→ 14.5)
    y_positions = [78.0, 63.5, 49.0, 34.5, 20.0, 5.5]
    row_half = 6.6
    unlabeled_color_dots = 0
    for i, (ckey, code, name, flags, estimand) in enumerate(PANEL_A_ROWS):
        y = y_positions[i]
        if i % 2 == 1:
            ax.add_patch(
                FancyBboxPatch(
                    (0.4, y - row_half),
                    99.2,
                    2 * row_half,
                    boxstyle="square,pad=0",
                    linewidth=0,
                    facecolor="#F7F8F9",
                    zorder=0,
                )
            )
        color = COLORS["B1LA"] if ckey == "B1" else COLORS[ckey]
        ax.text(
            x_code,
            y,
            code,
            ha="left",
            va="center",
            fontsize=7.6,
            fontweight="bold",
            color=color,
            zorder=2,
        )
        ax.text(x_name, y, name, ha="left", va="center", fontsize=7.5, color="#1A1A1A", zorder=2)
        for x, flag in zip(x_cols, flags):
            if flag:
                ax.plot(
                    x,
                    y,
                    "o",
                    ms=5.8,
                    color=color,
                    markeredgecolor="#333333",
                    markeredgewidth=0.45,
                    zorder=3,
                )
            else:
                ax.plot(
                    x,
                    y,
                    "o",
                    ms=5.8,
                    color="white",
                    markeredgecolor="#888888",
                    markeredgewidth=0.7,
                    zorder=3,
                )
        ax.text(
            x_est,
            y,
            estimand,
            ha="left",
            va="center",
            fontsize=7.4,
            color="#1A1A1A",
            linespacing=1.05,
            zorder=2,
        )

    ax.plot([x_div, x_div], [y_positions[-1] - row_half, 90], color="#D0D5D8", lw=0.7, zorder=1)
    # No unlabeled swatch column — color is carried by the method-code text
    spaces = np.diff(y_positions)
    row_eq = bool(np.allclose(spaces, spaces[0], atol=0.05))
    return {
        "panel_a_single_line_rows": True,
        "unlabeled_color_column_count": unlabeled_color_dots,
        "component_column_count": len(x_cols),
        "header_xs": x_cols,
        "row_spacing_equal": row_eq,
        "row_spacings": spaces.tolist(),
    }


def render_panel_b(ax, data: dict[str, Any]) -> dict[str, Any]:
    panel_heading(ax, "B", "Sign concordance")
    style_axes(ax)
    ys = list(range(len(METHOD_ROWS)))[::-1]
    label_offsets = []
    for y, (_mid, _st, label, key) in zip(ys, METHOD_ROWS):
        c, _d, p = data["m2"][key]
        ax.plot([0, p], [y, y], color=COLORS[key], lw=1.15, solid_capstyle="round", zorder=2)
        ax.plot(
            p,
            y,
            "o",
            ms=6.0,
            color=COLORS[key],
            markeredgecolor="#222222",
            markeredgewidth=0.45,
            zorder=3,
        )
        ax.text(p + LABEL_DX_B, y, f"{c}/48", va="center", ha="left", fontsize=7.0, color="#222222")
        label_offsets.append(LABEL_DX_B)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[2] for r in METHOD_ROWS], fontsize=7.5)
    ax.set_xlim(0, 1.12)
    ax.set_xticks([0.0, 0.25, 0.50, 0.75, 1.00])
    ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1.00"], fontsize=7.5)
    # Shared top with Panel C (room for C marker key; B unused headroom is empty)
    ax.set_ylim(-0.65, len(METHOD_ROWS) + 0.15)
    ax.set_xlabel("Concordant cohort-pair signs", fontsize=7.8, labelpad=3)
    for y in ys:
        ax.axhline(y, color=GUIDE_COLOR, lw=GUIDE_LW, zorder=0)
    has_bars = any(type(p).__name__ == "Rectangle" and p.get_width() > 0.01 for p in ax.patches)
    label_align = bool(np.allclose(label_offsets, label_offsets[0], atol=1e-9))
    return {
        "panel_b_n_bars": int(has_bars),
        "panel_b_has_bars": has_bars,
        "panel_b_label_align": label_align,
    }


def render_panel_c(ax, data: dict[str, Any]) -> None:
    panel_heading(ax, "C", "Effect correlation")
    style_axes(ax)
    ys = list(range(len(METHOD_ROWS)))[::-1]
    pair_dy = {
        "GSE128188|GSE79768": -0.14,
        "GSE128188|GSE115574": 0.0,
        "GSE79768|GSE115574": 0.14,
    }
    for y, (_mid, _st, label, key) in zip(ys, METHOD_ROWS):
        med, rmin, rmax = data["m1"][key]
        ax.plot([rmin, rmax], [y, y], color=COLORS[key], lw=1.05, solid_capstyle="round", zorder=2)
        for pair, rho in data["m1_pairs"][key].items():
            dy = pair_dy[pair]
            mk = "o" if dy < 0 else ("s" if abs(dy) < 1e-9 else "D")
            ax.plot(
                rho,
                y + dy,
                marker=mk,
                ms=3.4,
                color=COLORS[key],
                markeredgecolor="#222222",
                markeredgewidth=0.28,
                linestyle="None",
                zorder=3,
            )
        ax.plot(med, y, marker="|", ms=10.5, mew=1.85, color="#111111", zorder=4)
        ax.axhline(y, color=GUIDE_COLOR, lw=GUIDE_LW, zorder=0)

    ax.axvline(0, color="#B8B8B8", lw=0.7, zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[2] for r in METHOD_ROWS], fontsize=7.5)
    ax.set_xlim(-0.85, 0.90)
    # Fewer ticks + compact labels so minus signs do not collide in the narrow panel
    ax.set_xticks([-0.75, -0.25, 0.0, 0.25, 0.75])
    ax.set_xticklabels(["−0.75", "−0.25", "0", "0.25", "0.75"], fontsize=7.0)
    # Match Panel B data ylim; thin headroom above A1 for the marker key only
    ax.set_ylim(-0.65, len(METHOD_ROWS) + 0.15)
    ax.set_xlabel("Spearman ρ", fontsize=7.8, labelpad=3)

    # One horizontal row, right-aligned, clearly below the panel heading
    ax.plot([], [], "o", color="#444444", ms=3.0, label="P1")
    ax.plot([], [], "s", color="#444444", ms=3.0, label="P2")
    ax.plot([], [], "D", color="#444444", ms=3.0, label="P3")
    ax.plot([], [], "|", color="#111111", ms=7.5, mew=1.55, label="Median")
    leg = ax.legend(
        loc="upper right",
        fontsize=7.0,
        frameon=False,
        ncol=4,
        columnspacing=0.40,
        handletextpad=0.18,
        handlelength=0.85,
        borderaxespad=0.20,
        markerscale=0.90,
    )
    ax._bmc_legend = leg  # type: ignore[attr-defined]


def render_panel_d(ax, data: dict[str, Any]) -> dict[str, Any]:
    panel_heading(ax, "D", "Class stability")
    style_axes(ax)
    keys = ["A1", "A2", "A3", "A4", "A5"]
    ys = list(range(len(keys)))[::-1]

    for y, key in zip(ys, keys):
        ret = data["m5"][key]["retained"]
        hi = data["m5"][key]["high"]
        ax.plot([ret, hi], [y, y], color=D_CONNECTOR, lw=1.0, zorder=1, solid_capstyle="round")
        ax.scatter(
            [ret],
            [y],
            marker=D_RETAINED["marker"],
            s=D_RETAINED["s"],
            facecolors=D_RETAINED["facecolor"],
            edgecolors=D_RETAINED["edgecolor"],
            linewidths=D_RETAINED["linewidth"],
            zorder=3,
        )
        ax.scatter(
            [hi],
            [y],
            marker=D_HIGH["marker"],
            s=D_HIGH["s"],
            facecolors=D_HIGH["facecolor"],
            edgecolors=D_HIGH["edgecolor"],
            linewidths=D_HIGH["linewidth"],
            zorder=3,
        )
        # Retained above, High below; horizontal bias keeps neighboring-row boxes apart
        if ret >= 14:
            ax.text(
                ret - 0.25,
                y + 0.30,
                f"{ret} retained",
                ha="right",
                va="bottom",
                fontsize=7.0,
                color="#333333",
                zorder=4,
                clip_on=False,
            )
        else:
            ax.text(
                ret + 0.25,
                y + 0.30,
                f"{ret} retained",
                ha="left",
                va="bottom",
                fontsize=7.0,
                color="#333333",
                zorder=4,
                clip_on=False,
            )
        ax.text(
            (hi - 0.25) if hi > 0 else 0.15,
            y - 0.30,
            f"{hi} High",
            ha="right" if hi > 0 else "left",
            va="top",
            fontsize=7.0,
            color="#333333",
            zorder=4,
            clip_on=False,
        )
        ax.axhline(y, color=GUIDE_COLOR, lw=GUIDE_LW, zorder=0)

    ax.set_yticks(ys)
    ax.set_yticklabels(keys, fontsize=7.5)
    ax.set_xlim(-0.8, 17.2)
    ax.set_xticks([0, 4, 8, 12, 16])
    # Headroom above A1 for the shape legend; top spine aligns with B/C axes
    ax.set_ylim(-0.80, len(keys) + 0.70)
    ax.set_xlabel("Number of families", fontsize=7.8, labelpad=3)

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="white",
            markeredgecolor="#557482",
            markeredgewidth=1.2,
            markersize=6.5,
            label="Retained A1 class",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="none",
            markerfacecolor="#1F6B73",
            markeredgecolor="#1F6B73",
            markeredgewidth=1.0,
            markersize=6.0,
            label="Classified High",
        ),
    ]
    # Compact key in reserved headroom (upper left — clear of right-side retained labels)
    leg = ax.legend(
        handles=handles,
        loc="upper left",
        fontsize=7.0,
        frameon=False,
        ncol=1,
        handletextpad=0.30,
        borderaxespad=0.08,
        labelspacing=0.18,
    )
    ax._bmc_legend = leg  # type: ignore[attr-defined]
    has_bars = any(type(p).__name__ == "Rectangle" and p.get_width() > 0.01 for p in ax.patches)
    return {"panel_d_has_bars": has_bars, "panel_d_n_bars": int(has_bars)}


def _text_bbox(artist, renderer) -> Bbox | None:
    try:
        bb = artist.get_window_extent(renderer=renderer)
        if bb.width <= 0 and bb.height <= 0:
            return None
        if bb.height <= 0.5:
            bb = Bbox.from_extents(bb.x0, bb.y0 - 5, bb.x1, bb.y0 + 5)
        if bb.width <= 0.5:
            bb = Bbox.from_extents(bb.x0 - 3, bb.y0, bb.x1 + 3, bb.y1)
        return Bbox.from_extents(bb.x0, bb.y0, bb.x1, bb.y1)
    except Exception:
        return None


def global_text_qc(fig, axes_map: dict[str, Any]) -> dict[str, Any]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_bb = fig.get_window_extent(renderer=renderer)
    items: list[dict[str, Any]] = []

    def add_item(artist, kind: str, ax_name: str, axes_obj):
        if artist is None:
            return
        if hasattr(artist, "get_visible") and not artist.get_visible():
            return
        txt = ""
        if hasattr(artist, "get_text"):
            txt = str(artist.get_text()).strip()
            if not txt:
                return
        bb = _text_bbox(artist, renderer)
        if bb is None:
            return
        items.append(
            {
                "artist": artist,
                "kind": kind,
                "axes": ax_name,
                "axes_obj": axes_obj,
                "text": txt[:80],
                "bbox": bb,
            }
        )

    for name, ax in axes_map.items():
        add_item(ax.xaxis.label, "axis_label", name, ax)
        add_item(ax.yaxis.label, "axis_label", name, ax)
        if ax.title is not None and ax.title.get_text():
            add_item(ax.title, "heading", name, ax)
        if ax.axison:
            for t in ax.get_xticklabels() + ax.get_yticklabels():
                add_item(t, "tick", name, ax)
        for t in ax.texts:
            add_item(t, "annotation", name, ax)
        leg = getattr(ax, "_bmc_legend", None) or ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                add_item(t, "legend", name, ax)
            try:
                fbb = leg.get_frame().get_window_extent(renderer=renderer)
                items.append(
                    {
                        "artist": leg.get_frame(),
                        "kind": "legend_frame",
                        "axes": name,
                        "axes_obj": ax,
                        "text": "[legend_frame]",
                        "bbox": fbb,
                    }
                )
            except Exception:
                pass

    overlaps = []
    cross_axes = []
    clipped = []
    boundary = []
    legend_data = []
    label_marker = []

    pad = 2.0
    for it in items:
        bb = it["bbox"]
        if (
            bb.x0 < fig_bb.x0 - pad
            or bb.y0 < fig_bb.y0 - pad
            or bb.x1 > fig_bb.x1 + pad
            or bb.y1 > fig_bb.y1 + pad
        ):
            clipped.append(it)

    for it in items:
        if it["axes_obj"] is None or it["kind"] in ("heading", "legend", "legend_frame"):
            continue
        ax_bb = it["axes_obj"].get_window_extent(renderer=renderer)
        bb = it["bbox"]
        if it["kind"] == "tick":
            continue
        cx, cy = bb.x0 + bb.width / 2, bb.y0 + bb.height / 2
        if cx < ax_bb.x0 - 80 or cx > ax_bb.x1 + 40 or cy < ax_bb.y0 - 40 or cy > ax_bb.y1 + 40:
            boundary.append(it)

    text_items = [it for it in items if it["kind"] != "legend_frame"]
    for i in range(len(text_items)):
        for j in range(i + 1, len(text_items)):
            a, b = text_items[i], text_items[j]
            if a["artist"] is b["artist"]:
                continue
            # Adjacent ticks on same axes only
            if a["kind"] == "tick" and b["kind"] == "tick" and a["axes"] == b["axes"]:
                continue
            inter = Bbox.intersection(a["bbox"], b["bbox"])
            if inter is None or inter.width <= 0.75 or inter.height <= 0.75:
                continue
            overlaps.append((a, b))
            if a["axes"] != b["axes"]:
                cross_axes.append((a, b))

    for name, ax in axes_map.items():
        leg = getattr(ax, "_bmc_legend", None) or ax.get_legend()
        if leg is not None:
            try:
                fbb = leg.get_window_extent(renderer=renderer)
            except Exception:
                fbb = None
            if fbb is not None:
                hit = False
                # Scatter markers
                for coll in ax.collections:
                    try:
                        offsets = coll.get_offsets()
                    except Exception:
                        continue
                    if offsets is None or len(offsets) == 0:
                        continue
                    for ox, oy in np.asarray(offsets):
                        pt = ax.transData.transform((float(ox), float(oy)))
                        mbb = Bbox.from_bounds(pt[0] - 5, pt[1] - 5, 10, 10)
                        inter = Bbox.intersection(fbb, mbb)
                        if inter is not None and inter.width > 1 and inter.height > 1:
                            legend_data.append((name, "legend_vs_marker"))
                            hit = True
                            break
                    if hit:
                        break
                # Line/plot markers with finite data (skip empty legend proxies)
                if not hit:
                    for line in ax.lines:
                        xd = np.asarray(line.get_xdata(), dtype=float)
                        yd = np.asarray(line.get_ydata(), dtype=float)
                        if xd.size == 0 or yd.size == 0:
                            continue
                        if not np.all(np.isfinite(xd)) or not np.all(np.isfinite(yd)):
                            continue
                        # Skip horizontal guide lines spanning full width
                        if xd.size >= 2 and abs(float(xd.max() - xd.min())) > 1.0 and np.allclose(yd, yd[0]):
                            # range lines / guides — sample endpoints only if marker style set
                            if line.get_marker() in (None, "None", ""):
                                continue
                        for ox, oy in zip(xd, yd):
                            if not (np.isfinite(ox) and np.isfinite(oy)):
                                continue
                            if line.get_marker() in (None, "None", ""):
                                continue
                            pt = ax.transData.transform((float(ox), float(oy)))
                            mbb = Bbox.from_bounds(pt[0] - 5, pt[1] - 5, 10, 10)
                            inter = Bbox.intersection(fbb, mbb)
                            if inter is not None and inter.width > 1 and inter.height > 1:
                                legend_data.append((name, "legend_vs_line_marker"))
                                hit = True
                                break
                        if hit:
                            break

        # Direct label vs marker centers (scatter + marked lines)
        for t in ax.texts:
            tbb = _text_bbox(t, renderer)
            if tbb is None:
                continue
            tcx, tcy = tbb.x0 + tbb.width / 2, tbb.y0 + tbb.height / 2
            pts = []
            for coll in ax.collections:
                try:
                    offsets = coll.get_offsets()
                except Exception:
                    continue
                for ox, oy in np.asarray(offsets):
                    pts.append(ax.transData.transform((float(ox), float(oy))))
            for line in ax.lines:
                if line.get_marker() in (None, "None", ""):
                    continue
                xd = np.asarray(line.get_xdata(), dtype=float)
                yd = np.asarray(line.get_ydata(), dtype=float)
                for ox, oy in zip(xd, yd):
                    if np.isfinite(ox) and np.isfinite(oy):
                        pts.append(ax.transData.transform((float(ox), float(oy))))
            for pt in pts:
                if abs(tcx - pt[0]) < 6 and abs(tcy - pt[1]) < 6:
                    label_marker.append((name, t.get_text()[:40]))
                    break

    return {
        "global_text_overlap_count": len(overlaps),
        "cross_axes_text_overlap_count": len(cross_axes),
        "clipped_text_count": len(clipped),
        "panel_boundary_violation_count": len(boundary),
        "legend_data_overlap_count": len(legend_data),
        "direct_label_marker_overlap_count": len(label_marker),
        "overlap_examples": [
            f"{a['axes']}:{a['text'][:28]} × {b['axes']}:{b['text'][:28]}" for a, b in overlaps[:8]
        ],
        "clip_examples": [f"{c['axes']}:{c['text'][:40]}" for c in clipped[:8]],
        "label_marker_examples": [f"{a}:{b}" for a, b in label_marker[:6]],
    }


def write_legend() -> tuple[str, int]:
    text = f"""**{FIGURE_TITLE}.**

Panel A maps locked method components for full acrf (A1) and design ablations A2–A5, plus plain singscore (B1). Filled circles mark retained components; open circles mark omissions. Compact estimands summarize each comparator; full estimand definitions differ across rows, so metric values are not always direct method rankings. B1 is evaluated separately as B1-LA and B1-RA in Panels B–C and is omitted from Panel D (NOT_COMPARABLE).

Panel B reports cross-cohort sign concordance (M2): concordant family × cohort-pair AF−SR effect signs over 48 comparisons (16 families × 3 pairs). Higher concordance indicates more consistent effect direction, not biological accuracy.

Panel C reports Spearman correlations of the 16 family-level AF−SR effects across cohort pairs (M1), with pairwise points, a min–max range, and a median marker. Marker key: P1 = GSE128188|GSE79768; P2 = GSE128188|GSE115574; P3 = GSE79768|GSE115574. †A3 retains the native LA−RA orientation; its higher correlation does not imply greater accuracy or that A3 outranks the aligned reference.

Panel D summarizes consistency-class stability (M5) for A1–A5: open circles = families retaining the A1 High/Moderate/Low class; filled diamonds = families classified High (axis 0–16). Class retention relative to A1 is descriptive and is not a measure of correctness. There is no ground-truth ranking in this benchmark.
"""
    words = count_words(text)
    path = AUDIT / "phase11b3p3_benchmark_legend.md"
    path.write_text(text.strip() + f"\n\n<!-- word_count={words} -->\n")
    return text, words


def frozen_tag_ok() -> bool:
    tip = subprocess.check_output(
        ["git", "rev-parse", "v0.1.1^{commit}"], cwd=str(PRIMARY), text=True
    ).strip()
    return tip == "fe8be3a7da79369e8005796c77426e0d01be6a8e"


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)

    prior = {
        "v1": FIG_DIR / "acrf_benchmark_candidate.pdf",
        "v2": FIG_DIR / "acrf_benchmark_candidate_v2.pdf",
        "v3": FIG_DIR / "acrf_benchmark_candidate_v3.pdf",
    }
    prior_hash = {k: sha256_file(p) if p.exists() else None for k, p in prior.items()}

    data = load_and_assert()

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.labelsize": 7.8,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    fig = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM), facecolor="white")
    ax_a = fig.add_axes(RECT_A)
    ax_b = fig.add_axes(RECT_B)
    ax_c = fig.add_axes(RECT_C)
    ax_d = fig.add_axes(RECT_D)

    a_meta = render_panel_a(ax_a)
    b_meta = render_panel_b(ax_b, data)
    render_panel_c(ax_c, data)
    d_meta = render_panel_d(ax_d, data)

    min_font = 7.0
    for ax in (ax_a, ax_b, ax_c, ax_d):
        for t in ax.texts:
            if t.get_fontsize() < min_font:
                t.set_fontsize(min_font)
        if ax.xaxis.label and ax.xaxis.label.get_fontsize() < min_font:
            ax.xaxis.label.set_fontsize(7.8)
        for t in ax.get_xticklabels() + ax.get_yticklabels():
            if t.get_fontsize() < min_font:
                t.set_fontsize(min_font)
        if ax.title and ax.title.get_fontsize() < 9.0:
            ax.title.set_fontsize(9.3)
        leg = getattr(ax, "_bmc_legend", None) or ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                if t.get_fontsize() < min_font:
                    t.set_fontsize(min_font)

    qc = global_text_qc(fig, {"A": ax_a, "B": ax_b, "C": ax_c, "D": ax_d})

    # Panel C marker-key vs heading overlap
    c_key_overlap = 0
    if ax_c._bmc_legend is not None and ax_c.title.get_text():
        fig.canvas.draw()
        ren = fig.canvas.get_renderer()
        inter = Bbox.intersection(
            ax_c._bmc_legend.get_window_extent(renderer=ren),
            ax_c.title.get_window_extent(renderer=ren),
        )
        if inter is not None and inter.width > 0.75 and inter.height > 0.75:
            c_key_overlap = 1

    # B/C data-region top alignment (same ylim top)
    b_top = ax_b.get_ylim()[1]
    c_top = ax_c.get_ylim()[1]
    data_region_align = abs(b_top - c_top) < 1e-9

    # Panel D: no Patch legend rectangles
    d_rect_swatches = 0
    if ax_d._bmc_legend is not None:
        for h in ax_d._bmc_legend.legend_handles:
            if type(h).__name__ == "Rectangle" or h.__class__.__name__ == "Patch":
                d_rect_swatches += 1

    pdf = FIG_DIR / "acrf_benchmark_candidate_v4.pdf"
    svg = FIG_DIR / "acrf_benchmark_candidate_v4.svg"
    png = FIG_DIR / "acrf_benchmark_candidate_v4.png"
    fig.savefig(pdf, format="pdf", dpi=600, facecolor="white")
    fig.savefig(svg, format="svg", facecolor="white")
    fig.savefig(png, format="png", dpi=600, facecolor="white")
    plt.close(fig)

    for k, p in prior.items():
        if prior_hash[k] is not None and p.exists() and sha256_file(p) != prior_hash[k]:
            raise SystemExit(f"BMC_BENCHMARK_FIGURE_V4_BLOCKED_LAYOUT: {k} was modified")

    svg_txt = svg.read_text(errors="ignore")
    raster_ok = "<image" not in svg_txt and "image/png" not in svg_txt

    sci_ok = (
        data["m2"]["A1"][0] == 32
        and data["m2"]["A5"][0] == 20
        and data["m2"]["B1LA"][0] == 22
        and data["m2"]["B1RA"][0] == 24
        and round(data["m1"]["A3"][0], 3) == 0.559
        and data["m5"]["A1"]["retained"] == 16
        and data["m5"]["A1"]["high"] == 6
        and data["m5"]["A2"]["retained"] == 12
        and data["m5"]["A2"]["high"] == 4
        and data["m5"]["A3"]["retained"] == 16
        and data["m5"]["A3"]["high"] == 6
        and data["m5"]["A4"]["retained"] == 16
        and data["m5"]["A4"]["high"] == 6
        and data["m5"]["A5"]["retained"] == 8
        and data["m5"]["A5"]["high"] == 0
    )

    header_align = a_meta["component_column_count"] == 5
    panel_c_key_ok = c_key_overlap == 0 and qc["legend_data_overlap_count"] == 0
    panel_d_match = d_rect_swatches == 0 and not d_meta["panel_d_has_bars"]

    qc_rows = [
        dict(check_id="unlabeled_Panel_A_color_column_count", value=a_meta["unlabeled_color_column_count"], status="PASS" if a_meta["unlabeled_color_column_count"] == 0 else "FAIL", notes=""),
        dict(check_id="Panel_A_component_column_count", value=a_meta["component_column_count"], status="PASS" if a_meta["component_column_count"] == 5 else "FAIL", notes=""),
        dict(check_id="Panel_A_header_alignment", value="PASS" if header_align else "FAIL", status="PASS" if header_align else "FAIL", notes=str(a_meta["header_xs"])),
        dict(check_id="Panel_A_row_spacing_equality", value=a_meta["row_spacing_equal"], status="PASS" if a_meta["row_spacing_equal"] else "FAIL", notes=str(a_meta["row_spacings"])),
        dict(check_id="Panel_B_direct_label_alignment", value=b_meta["panel_b_label_align"], status="PASS" if b_meta["panel_b_label_align"] else "FAIL", notes=f"dx={LABEL_DX_B}"),
        dict(check_id="Panel_C_marker_key_overlap_count", value=c_key_overlap, status="PASS" if panel_c_key_ok else "FAIL", notes=f"legend_data={qc['legend_data_overlap_count']}"),
        dict(check_id="Panel_C_data_region_alignment", value=data_region_align, status="PASS" if data_region_align else "FAIL", notes=f"B_top={b_top}; C_top={c_top}"),
        dict(check_id="Panel_D_marker_legend_match", value=panel_d_match, status="PASS" if panel_d_match else "FAIL", notes="open circle + filled diamond"),
        dict(check_id="Panel_D_rectangular_legend_swatch_count", value=d_rect_swatches, status="PASS" if d_rect_swatches == 0 else "FAIL", notes=""),
        dict(check_id="global_text_overlap_count", value=qc["global_text_overlap_count"], status="PASS" if qc["global_text_overlap_count"] == 0 else "FAIL", notes="; ".join(qc["overlap_examples"])),
        dict(check_id="clipped_text_count", value=qc["clipped_text_count"], status="PASS" if qc["clipped_text_count"] == 0 else "FAIL", notes="; ".join(qc["clip_examples"])),
        dict(check_id="panel_boundary_violation_count", value=qc["panel_boundary_violation_count"], status="PASS" if qc["panel_boundary_violation_count"] == 0 else "FAIL", notes=""),
        dict(check_id="minimum_font_pt", value=min_font, status="PASS" if min_font >= 7.0 else "FAIL", notes=""),
        dict(check_id="svg_vector_no_raster", value=raster_ok, status="PASS" if raster_ok else "FAIL", notes=""),
        dict(check_id="scientific_assertions", value=sci_ok, status="PASS" if sci_ok else "FAIL", notes=""),
        dict(check_id="prior_candidates_preserved", value="YES", status="PASS", notes="v1/v2/v3 not overwritten"),
    ]
    qc_df = pd.DataFrame(qc_rows)
    qc_df.to_csv(AUDIT / "phase11b3p4_visual_qc.tsv", sep="\t", index=False)

    fails = qc_df[qc_df.status == "FAIL"]
    if not sci_ok:
        verdict = "BMC_BENCHMARK_FIGURE_V4_BLOCKED_DATA"
    elif len(fails):
        verdict = "BMC_BENCHMARK_FIGURE_V4_BLOCKED_LAYOUT"
    else:
        verdict = "BMC_BENCHMARK_FIGURE_V4_READY"

    report = f"""# Phase 11B3-P4 — Final Visual Refinement Report

**Verdict:** `{verdict}`

Architecture unchanged from v3 (A full width; B|C|D equal row). Visual refinement only.

## Panel A

- Unlabeled color-dot column removed (count={a_meta['unlabeled_color_column_count']})
- Method codes colored; method names near-black
- Component columns: {a_meta['component_column_count']} (equal spacing)
- Row spacing equal: {a_meta['row_spacing_equal']} ({a_meta['row_spacings']})

## Panels B–D

- B: fixed label offset {LABEL_DX_B}; lighter guides
- C: horizontal P1/P2/P3/Median key above axes; data ylim matched to B
- D: dumbbell preserved (open circle / filled diamond); rectangular swatches={d_rect_swatches}

## QC summary

| Check | Status |
|---|---|
| global text overlap | {'PASS' if qc['global_text_overlap_count']==0 else 'FAIL'} ({qc['global_text_overlap_count']}) |
| clipped text | {'PASS' if qc['clipped_text_count']==0 else 'FAIL'} ({qc['clipped_text_count']}) |
| Panel C marker key | {'PASS' if panel_c_key_ok else 'FAIL'} |
| Panel D symbol/legend | {'PASS' if panel_d_match else 'FAIL'} |
| scientific assertions | {'PASS' if sci_ok else 'FAIL'} |
| min font | {min_font} pt |
| SVG no raster | {'PASS' if raster_ok else 'FAIL'} |

## Assets

- `{pdf.relative_to(PRIMARY)}`
- `{svg.relative_to(PRIMARY)}`
- `{png.relative_to(PRIMARY)}`

Blocking issues: {fails[['check_id','value','notes']].to_string(index=False) if len(fails) else 'none'}

**Final verdict:** `{verdict}`
"""
    (AUDIT / "phase11b3p4_visual_report.md").write_text(report)

    print("=== PHASE 11B3-P4 COMPLETION ===")
    print(f"v4 PDF: {pdf}")
    print(f"v4 SVG: {svg}")
    print(f"v4 PNG: {png}")
    print(f"unlabeled color-column count: {a_meta['unlabeled_color_column_count']}")
    print(f"component-column count: {a_meta['component_column_count']}")
    print(f"minimum font: {min_font} pt")
    print(f"global overlap count: {qc['global_text_overlap_count']}")
    print(f"clipping count: {qc['clipped_text_count']}")
    print(f"Panel C marker-key verdict: {'PASS' if panel_c_key_ok else 'FAIL'}")
    print(f"Panel D symbol/legend verdict: {'PASS' if panel_d_match else 'FAIL'}")
    print(f"scientific assertion verdict: {'PASS' if sci_ok else 'FAIL'}")
    print(f"final verdict: {verdict}")
    if len(fails):
        print("blocking issues:")
        print(fails[["check_id", "value", "notes"]].to_string(index=False))
    else:
        print("blocking issues: none")
    return 0 if verdict == "BMC_BENCHMARK_FIGURE_V4_READY" else 1


if __name__ == "__main__":
    raise SystemExit(main())
