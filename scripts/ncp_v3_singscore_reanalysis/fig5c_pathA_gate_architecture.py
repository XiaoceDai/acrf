#!/usr/bin/env python3
"""Figure 5C — rule-based Path A prioritization.

Reads only the Figure 5 canonical table. Does not modify 5A/B/D or formal
routing. Does not open PriorityScore, GeoMx, genetics, or Figures 2-4 sources.
"""
from __future__ import annotations

import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"

IN_CANONICAL = RESULTS / "figure5_triangulation_canonical.tsv"

OUT_BASE = FIGURES / "fig5c_pathA_gate_architecture"
OUT_PREVIEW = FIGURES / "fig5c_pathA_gate_architecture_preview.png"
OUT_SOURCE = RESULTS / "fig5c_pathA_gate_architecture_source.tsv"
OUT_SUMMARY = RESULTS / "fig5c_pathA_gate_architecture_summary.tsv"
OUT_FAMILY = RESULTS / "fig5c_pathA_family_detail.tsv"
OUT_QC = RESULTS / "fig5c_pathA_gate_architecture_qc.md"

FORBIDDEN_READ_PATTERNS = (
    re.compile(r"PriorityScore"),
    re.compile(r"results/fig5b_v3_triangulation\.tsv"),
    re.compile(r"fig2e_cross_cohort_consistency"),
    re.compile(r"fig2d_family_effects"),
    re.compile(r"fig3_v3_carrier_fidelity"),
    re.compile(r"fig4_v3_chronicity"),
    re.compile(r"fig1f_v3_HCA"),
    re.compile(r"geomx", re.IGNORECASE),
    re.compile(r"genetics"),
    re.compile(r"roi_scores_v3"),
)

MM = 1 / 25.4
WIDTH_MM = 84.0
HEIGHT_MM = 88.0
MARGIN_MM = 4.0
COL_WIDTHS = [12.3, 8.0, 23.0, 15.5, 6.0, 11.2]

FG = "#1A1A1A"
MUTED = "#595959"
FAINT = "#7A7A7A"
WHITE = "#FFFFFF"
BG = "#FAFAFA"
PATH_A_EDGE = "#8A8A8A"

GREY_SEQ = {
    3: "#4A4A4A",
    2: "#6E6E6E",
    1: "#9A9A9A",
    0: "#C6C6C6",
}

# Absolute type system (shared with Figure 5D in composite)
PANEL_LETTER_PT = 10.8
TITLE_PT = 10.8
TITLE_TEXT_PT = 8.6
SUBTITLE_PT = 6.5
BODY_PT = 6.2
NOTE_PT = 6.0
FOOTER_PT = 5.5

_STANDALONE = {
    "WIDTH_MM": 84.0,
    "HEIGHT_MM": 88.0,
    "MARGIN_MM": 4.0,
    "COL_WIDTHS": [12.3, 8.0, 23.0, 15.5, 6.0, 11.2],
    "PANEL_LETTER_PT": 10.8,
    "TITLE_TEXT_PT": 8.6,
    "SUBTITLE_PT": 6.3,
    "BODY_PT": 6.2,
    "NOTE_PT": 6.0,
    "FOOTER_PT": 6.0,
}
_COMPOSITE = {
    "WIDTH_MM": 76.0,
    "HEIGHT_MM": 88.0,
    "MARGIN_MM": 2.5,  # shared top margin with composite D letter baseline
    "COL_WIDTHS": [11.5, 7.5, 21.5, 14.5, 5.6, 10.4],  # sum = 71 (= 76 - 2*2.5)
    # One-step type bump for composite readability
    "PANEL_LETTER_PT": 11.0,
    "TITLE_TEXT_PT": 9.0,
    "SUBTITLE_PT": 7.0,
    "BODY_PT": 6.7,
    "NOTE_PT": 6.5,
    "FOOTER_PT": 5.8,
}


def apply_layout(mode: str = "standalone") -> None:
    """Switch between standalone and Figure 5 composite canvas geometry."""
    global WIDTH_MM, HEIGHT_MM, MARGIN_MM, COL_WIDTHS
    global PANEL_LETTER_PT, TITLE_TEXT_PT, SUBTITLE_PT, BODY_PT, NOTE_PT, FOOTER_PT
    cfg = _COMPOSITE if mode == "composite" else _STANDALONE
    WIDTH_MM = cfg["WIDTH_MM"]
    HEIGHT_MM = cfg["HEIGHT_MM"]
    MARGIN_MM = cfg["MARGIN_MM"]
    COL_WIDTHS = list(cfg["COL_WIDTHS"])
    PANEL_LETTER_PT = cfg["PANEL_LETTER_PT"]
    TITLE_TEXT_PT = cfg["TITLE_TEXT_PT"]
    SUBTITLE_PT = cfg["SUBTITLE_PT"]
    BODY_PT = cfg["BODY_PT"]
    NOTE_PT = cfg["NOTE_PT"]
    FOOTER_PT = cfg["FOOTER_PT"]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_no_forbidden_reads() -> None:
    src = Path(__file__).read_text(encoding="utf-8")
    for line in src.splitlines():
        if not line.strip().startswith("#") and "pd.read_csv(" in line:
            for pattern in FORBIDDEN_READ_PATTERNS:
                require(pattern.search(line) is None, f"Forbidden read pattern in: {line.strip()}")


def load_and_assert() -> pd.DataFrame:
    df = pd.read_csv(IN_CANONICAL, sep="\t")
    require(set(df.columns) >= {"family_id", "family_label", "bulk_gate", "sn_gate", "duration_gate", "support_count", "Path_A"}, "missing columns")

    # Canonical assertions
    require(len(df) == 16, f"expected 16 families, got {len(df)}")
    require(df["family_id"].nunique() == 16, "duplicate family_id")

    path_a = set(df.loc[df["Path_A"], "family_id"])
    require(path_a == {10, 15}, f"Path A families {path_a} != {{10,15}}")

    # Compute combinations
    combo = df.copy()
    combo["gate_sum"] = combo["bulk_gate"].astype(int) + combo["sn_gate"].astype(int) + combo["duration_gate"].astype(int)
    require(set(combo["gate_sum"]) == {0, 1, 2, 3}, f"unexpected gate sums {set(combo['gate_sum'])}")
    require(combo["gate_sum"].value_counts().to_dict() == {3: 2, 2: 1, 1: 10, 0: 3}, "support-count distribution mismatch")

    # Observed combination buckets
    bulk = combo["bulk_gate"].astype(bool)
    sn = combo["sn_gate"].astype(bool)
    dur = combo["duration_gate"].astype(bool)

    def families(where: pd.Series) -> set[int]:
        return set(combo.loc[where, "family_id"])

    all_three = families(bulk & sn & dur)
    bulk_dur = families(bulk & ~sn & dur)
    bulk_only = families(bulk & ~sn & ~dur)
    sn_only = families(~bulk & sn & ~dur)
    dur_only = families(~bulk & ~sn & dur)
    no_gates = families(~bulk & ~sn & ~dur)
    bulk_sn = families(bulk & sn & ~dur)
    sn_dur = families(~bulk & sn & dur)

    require(all_three == {10, 15}, f"all-three families {all_three}")
    require(bulk_dur == {9}, f"bulk+duration families {bulk_dur}")
    require(bulk_only == {1, 12, 14}, f"bulk-only families {bulk_only}")
    require(sn_only == {11, 16}, f"snRNA-only families {sn_only}")
    require(dur_only == {2, 3, 7, 8, 13}, f"duration-only families {dur_only}")
    require(no_gates == {4, 5, 6}, f"no-gate families {no_gates}")
    require(bulk_sn == set(), "unexpected bulk+snRNA-only families")
    require(sn_dur == set(), "unexpected snRNA+duration-only families")

    all_observed = all_three | bulk_dur | bulk_only | sn_only | dur_only | no_gates
    require(all_observed == set(df["family_id"]), "observed combinations do not cover all families")
    require(len(all_observed) == 16, "observed combination count mismatch")

    # No forbidden fields used
    for col in ("hca_carrier", "bulk_dominant_state", "bulk_consistency", "snRNA_class", "duration_class"):
        require(col in combo.columns, f"context column {col} missing")
    for forbidden in ("PriorityScore", "GeoMx", "genetics"):
        require(forbidden not in " ".join(combo.columns), f"forbidden column {forbidden}")

    return combo


def build_sources(combo: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Support-count composition
    counts = combo["gate_sum"].value_counts().sort_index(ascending=False).rename_axis("support_count").reset_index(name="n_families")
    source = counts.copy()
    source["fraction"] = (source["n_families"] / 16.0).round(4)

    # Observed combinations
    rows = []
    for bulk, sn, dur, name in [
        (True, True, True, "all_three"),
        (True, False, True, "bulk_duration"),
        (True, False, False, "bulk_only"),
        (False, True, False, "snRNA_only"),
        (False, False, True, "duration_only"),
        (False, False, False, "no_gates"),
    ]:
        mask = (combo["bulk_gate"] == bulk) & (combo["sn_gate"] == sn) & (combo["duration_gate"] == dur)
        fams = sorted(combo.loc[mask, "family_id"])
        rows.append({
            "combination_name": name,
            "bulk_gate": bulk,
            "sn_gate": sn,
            "duration_gate": dur,
            "n_families": len(fams),
            "family_ids": ",".join(f"F{i}" for i in fams),
        })
    combo_source = pd.DataFrame(rows)

    # Path A family detail
    detail = combo.loc[combo["Path_A"], ["family_id", "family_label"]].copy()
    detail["family_id"] = detail["family_id"].apply(lambda x: f"F{x}")
    detail = detail.rename(columns={"family_id": "family", "family_label": "label"})

    return source, combo_source, detail



def format_family_ids(family_ids: str) -> str:
    """Return middle-dot separated family IDs, wrapped to fit narrow column."""
    if not family_ids:
        return "—"
    ids = [x.strip() for x in family_ids.split(",")]
    if len(ids) == 1:
        return ids[0]
    if len(ids) == 2:
        return " · ".join(ids)
    if len(ids) == 3:
        return f"{ids[0]} · {ids[1]}\n{ids[2]}"
    if len(ids) == 5:
        return f"{ids[0]} · {ids[1]}\n{ids[2]} · {ids[3]}\n{ids[4]}"
    return " · ".join(ids)


def render(source: pd.DataFrame, combo_source: pd.DataFrame, detail: pd.DataFrame) -> tuple[plt.Figure, plt.Axes]:
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Liberation Sans"]
    mpl.rcParams["axes.unicode_minus"] = True

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM), facecolor=WHITE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, WIDTH_MM)
    ax.set_ylim(0, HEIGHT_MM)
    ax.axis("off")
    ax.set_facecolor(WHITE)

    # Title / subtitle (compact 5.5 mm title block)
    y_title = HEIGHT_MM - MARGIN_MM - 0.2
    ax.text(MARGIN_MM, y_title, "C", ha="left", va="top", fontsize=PANEL_LETTER_PT, fontweight="bold", color=FG)
    ax.text(
        MARGIN_MM + 6.5,
        y_title,
        "Rule-based Path A prioritization",
        ha="left",
        va="top",
        fontsize=TITLE_TEXT_PT,
        fontweight="bold",
        color=FG,
    )
    subtitle = (
        "Observed gate combinations\nacross 16 frozen NCP families"
        if WIDTH_MM <= 78.0
        else "Observed gate combinations across 16 frozen NCP families"
    )
    ax.text(
        MARGIN_MM + 6.5,
        y_title - 3.5,
        subtitle,
        ha="left",
        va="top",
        fontsize=SUBTITLE_PT,
        color=FAINT,
        linespacing=0.9,
    )

    # Top section: composition bar for half-panel width
    top_y = y_title - (8.2 if WIDTH_MM <= 78.0 else 5.5)
    bar_h = 6.0
    bar_w = WIDTH_MM - 2 * MARGIN_MM
    counts = {int(row.support_count): int(row.n_families) for row in source.itertuples(index=False)}
    total = sum(counts.values())

    seg_x: dict[int, tuple[float, float]] = {}
    x = MARGIN_MM
    for n in (3, 2, 1, 0):
        w = bar_w * counts[n] / total
        seg_x[n] = (x, x + w)
        ax.add_patch(Rectangle((x, top_y - bar_h), w, bar_h, facecolor=GREY_SEQ[n], edgecolor="none"))
        ax.text(x + w / 2, top_y - bar_h / 2, str(counts[n]), ha="center", va="center", fontsize=BODY_PT, color=WHITE, fontweight="bold")
        x += w

    # Fraction labels outside the bar; leader for the narrow 2/3 segment
    label_y = top_y - bar_h - 1.4
    for n, frac in [(3, "3/3"), (2, "2/3"), (1, "1/3"), (0, "0/3")]:
        x0, x1 = seg_x[n]
        cx = (x0 + x1) / 2
        if n == 2:
            lx = x1 + 1.0
            ax.text(lx, label_y, frac, ha="left", va="center", fontsize=NOTE_PT, color=MUTED, fontweight="bold")
            ax.plot([x1, lx - 0.5], [top_y - bar_h / 2, label_y], color=MUTED, lw=0.5, zorder=1)
        else:
            ax.text(cx, label_y, frac, ha="center", va="center", fontsize=NOTE_PT, color=MUTED, fontweight="bold")

    # Small Path A label beneath the 3/3 fraction label
    ax.text(
        (seg_x[3][0] + seg_x[3][1]) / 2,
        label_y - 2.8,
        "Path A",
        ha="center",
        va="center",
        fontsize=NOTE_PT,
        color=MUTED,
        fontweight="bold",
    )

    # Bottom section: gate combination matrix
    matrix_top = top_y - bar_h - 6.0
    matrix_left = MARGIN_MM
    matrix_right = WIDTH_MM - MARGIN_MM
    matrix_w = matrix_right - matrix_left

    col_labels = ["Bulk\nHigh", "LA snRNA\ngate:\nnon-\nBelow-detection", "LAA\nduration\nnon-\nStable", "n", "Families"]
    col_widths = list(COL_WIDTHS)
    require(abs(sum(col_widths) - matrix_w) < 0.1, f"column widths sum {sum(col_widths)} != {matrix_w}")
    col_x = [matrix_left]
    for w in col_widths[:-1]:
        col_x.append(col_x[-1] + w)
    col_cx = [x + w / 2 for x, w in zip(col_x, col_widths)]

    matrix_rows = [
        ("all_three", "All three\n/ Path A", True, True, True, True),
        ("bulk_duration", "Bulk +\nduration", True, False, True, False),
        ("bulk_only", "Bulk only", True, False, False, False),
        ("snRNA_only", "snRNA\nonly", False, True, False, False),
        ("duration_only", "Duration\nonly", False, False, True, False),
        ("no_gates", "No gates", False, False, False, False),
    ]

    compact = WIDTH_MM <= 78.0
    row_h = 5.7 if compact else 6.5
    row_gap = 0.35 if compact else 0.5
    header_h = 9.5 if compact else 10.5
    n_matrix_rows = len(matrix_rows)
    matrix_h = header_h + 0.4 + row_h / 2 + (n_matrix_rows - 1) * (row_h + row_gap) + row_h / 2

    # Draw matrix background
    ax.add_patch(
        Rectangle(
            (matrix_left, matrix_top - matrix_h),
            matrix_w,
            matrix_h,
            facecolor=BG,
            edgecolor="#E0D0D0",
            lw=0.5,
            zorder=0,
        )
    )

    # Header row (combination column has no header; labels cover the remaining columns)
    header_y = matrix_top - header_h / 2
    for cx, label in zip(col_cx[1:], col_labels):
        ax.text(cx, header_y, label, ha="center", va="center", fontsize=BODY_PT, color=MUTED, fontweight="bold", linespacing=0.85)

    # Header rule
    header_rule_y = matrix_top - header_h
    ax.plot([matrix_left, matrix_right], [header_rule_y, header_rule_y], color="#D0D0D0", lw=0.55, zorder=1)

    # Row lines
    y_cursor = header_rule_y - row_h / 2 - 0.4
    for idx, (name, combo_label, bulk, sn, dur, is_path_a) in enumerate(matrix_rows):
        row = combo_source.loc[combo_source["combination_name"] == name].iloc[0]
        fams = row["family_ids"]
        n = int(row["n_families"])

        # Alternate row background
        if idx % 2 == 1:
            ax.add_patch(
                Rectangle(
                    (matrix_left, y_cursor - row_h / 2),
                    matrix_w,
                    row_h,
                    facecolor="#F3F3F3",
                    edgecolor="none",
                    zorder=0,
                )
            )

        # Observed combination label
        ax.text(col_x[0] + 0.5, y_cursor, combo_label, ha="left", va="center", fontsize=BODY_PT, color=FG, linespacing=0.85)

        # Gate symbols in the three gate columns
        for col_idx, satisfied in enumerate([bulk, sn, dur]):
            cx = col_cx[col_idx + 1]
            if satisfied:
                ax.add_patch(plt.Circle((cx, y_cursor), 1.15, facecolor="#4A4A4A", edgecolor="none", zorder=2))
            else:
                ax.add_patch(plt.Circle((cx, y_cursor), 1.15, facecolor=WHITE, edgecolor="#4A4A4A", linewidth=0.8, zorder=2))

        # n
        ax.text(col_cx[4], y_cursor, str(n), ha="center", va="center", fontsize=BODY_PT, color=FG, fontweight="bold")

        # Families
        if is_path_a:
            fam_text = "F10 · F15"
        else:
            fam_text = format_family_ids(fams)
        ax.text(col_x[5] + 0.3, y_cursor, fam_text, ha="left", va="center", fontsize=BODY_PT, color=FG, linespacing=0.85)

        y_cursor -= row_h + row_gap

    # Path A callout below matrix (full family names with badge)
    matrix_bottom = matrix_top - matrix_h
    callout_top = matrix_bottom - (0.8 if compact else 0.2)
    callout_h = 2.6 if compact else 2.8
    badge_w, badge_h = 13.0, 2.6
    badge_x = matrix_left
    badge_y = callout_top - callout_h / 2
    ax.add_patch(
        FancyBboxPatch(
            (badge_x, badge_y - badge_h / 2),
            badge_w,
            badge_h,
            boxstyle="round,pad=0.06,rounding_size=0.3",
            facecolor=WHITE,
            edgecolor=PATH_A_EDGE,
            linewidth=0.7,
            zorder=1,
        )
    )
    ax.text(badge_x + badge_w / 2, badge_y, "Path A", ha="center", va="center", fontsize=BODY_PT, color=FG, fontweight="bold", zorder=2)

    f10_label = detail.loc[detail["family"] == "F10", "label"].iloc[0]
    f15_label = detail.loc[detail["family"] == "F15", "label"].iloc[0]
    ax.text(
        badge_x + badge_w + 1.0,
        badge_y,
        f"F10 {f10_label}    F15 {f15_label}",
        ha="left",
        va="center",
        fontsize=BODY_PT,
        color=FG,
    )

    # Footer: two required lines; optional third only when vertical space allows
    footer_y = 2.8 if compact else 4.0
    footer_gap = 2.8 if compact else 2.3
    ax.text(
        MARGIN_MM,
        footer_y,
        "● gate met; ○ gate unmet. Path A requires 3/3 gates.",
        ha="left",
        va="bottom",
        fontsize=FOOTER_PT,
        color=FAINT,
    )
    ax.text(
        MARGIN_MM,
        footer_y + footer_gap,
        "Context annotations, GeoMx, and genetics are not counted.",
        ha="left",
        va="bottom",
        fontsize=FOOTER_PT,
        color=FAINT,
    )
    if not compact:
        ax.text(
            MARGIN_MM,
            footer_y + 2 * footer_gap,
            "No bulk+snRNA-only or snRNA+duration-only families.",
            ha="left",
            va="bottom",
            fontsize=FOOTER_PT,
            color=FAINT,
        )

    return fig, ax


def validate_layout(fig: plt.Figure, ax: plt.Axes) -> dict[str, float | int]:
    fig.canvas.draw()
    texts = []
    renderer = fig.canvas.get_renderer()
    for artist in ax.get_children():
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
        texts.append((label, (x0, y0, x1, y1), float(artist.get_fontsize())))

    min_font = min(fs for _, _, fs in texts)
    require(min_font >= 5.5 - 0.05, f"text below 5.5 pt: {min_font:.2f}")

    for i, (label_a, bbox_a, _) in enumerate(texts):
        for label_b, bbox_b, _ in texts[i + 1 :]:
            if label_a == label_b:
                continue
            x0 = max(bbox_a[0], bbox_b[0])
            y0 = max(bbox_a[1], bbox_b[1])
            x1 = min(bbox_a[2], bbox_b[2])
            y1 = min(bbox_a[3], bbox_b[3])
            if x1 > x0 and y1 > y0:
                overlap = (x1 - x0) * (y1 - y0)
                require(overlap < 0.15, f"text overlap: {label_a!r} vs {label_b!r} ({overlap:.2f} mm²)")

    return {"n_text_objects": len(texts), "min_font_pt": min_font}


def write_summary(source: pd.DataFrame, combo_source: pd.DataFrame, detail: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {"metric": "n_families", "value": 16},
        {"metric": "Path_A_families", "value": ",".join(detail["family"])},
        {"metric": "support_3", "value": int(source.loc[source["support_count"] == 3, "n_families"].iloc[0])},
        {"metric": "support_2", "value": int(source.loc[source["support_count"] == 2, "n_families"].iloc[0])},
        {"metric": "support_1", "value": int(source.loc[source["support_count"] == 1, "n_families"].iloc[0])},
        {"metric": "support_0", "value": int(source.loc[source["support_count"] == 0, "n_families"].iloc[0])},
        {"metric": "combination_all_three", "value": combo_source.loc[combo_source["combination_name"] == "all_three", "family_ids"].iloc[0]},
        {"metric": "combination_bulk_duration", "value": combo_source.loc[combo_source["combination_name"] == "bulk_duration", "family_ids"].iloc[0]},
        {"metric": "combination_bulk_only", "value": combo_source.loc[combo_source["combination_name"] == "bulk_only", "family_ids"].iloc[0]},
        {"metric": "combination_snRNA_only", "value": combo_source.loc[combo_source["combination_name"] == "snRNA_only", "family_ids"].iloc[0]},
        {"metric": "combination_duration_only", "value": combo_source.loc[combo_source["combination_name"] == "duration_only", "family_ids"].iloc[0]},
        {"metric": "combination_no_gates", "value": combo_source.loc[combo_source["combination_name"] == "no_gates", "family_ids"].iloc[0]},
        {"metric": "forbidden_fields_used", "value": False},
        {"metric": "canvas_mm", "value": f"{WIDTH_MM:.0f}x{HEIGHT_MM:.0f}"},
        {"metric": "assertion_status", "value": "PASS"},
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    return summary


def write_qc(hashes: dict[str, str], layout_stats: dict[str, float | int]) -> None:
    OUT_QC.write_text(
        f"""# Figure 5C Path A gate architecture QC

Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}

Status: **PASS**

## Layout checks
- Text objects: {layout_stats["n_text_objects"]}
- Minimum font: {layout_stats["min_font_pt"]:.1f} pt
- No text overlap, no canvas clipping (validated)

## Inputs
- `{IN_CANONICAL.relative_to(ROOT)}` — `{sha256(IN_CANONICAL)}`

## Outputs
{chr(10).join(f"- `{path}` — `{digest}`" for path, digest in sorted(hashes.items()))}
""",
        encoding="utf-8",
    )


def main() -> None:
    assert_no_forbidden_reads()
    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    combo = load_and_assert()
    source, combo_source, detail = build_sources(combo)
    source.to_csv(OUT_SOURCE, sep="\t", index=False)
    detail.to_csv(OUT_FAMILY, sep="\t", index=False)
    write_summary(source, combo_source, detail)

    fig, ax = render(source, combo_source, detail)
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
    require("Rule-based Path A prioritization" in svg_text, "title missing")
    require("Observed gate combinations across 16 frozen NCP families" in svg_text, "subtitle missing")
    require("Path A" in svg_text, "Path A tag missing")
    require(
        "F10" in svg_text
        and "GAG disease" in svg_text
        and "TGF-β" in svg_text
        and "F15" in svg_text
        and "NTRK2" in svg_text
        and "TrkB" in svg_text,
        "Path A family labels missing",
    )
    require("● gate met" in svg_text, "footer missing")
    require("Context annotations, GeoMx" in svg_text, "footer missing")
    require("No bulk+snRNA-only" in svg_text, "note missing")

    from PIL import Image

    img = Image.open(png)
    w_mm = img.size[0] * 25.4 / 600.0
    h_mm = img.size[1] * 25.4 / 600.0
    require(abs(w_mm - WIDTH_MM) < 1.0, f"PNG width {w_mm:.2f} != {WIDTH_MM}")
    require(abs(h_mm - HEIGHT_MM) < 1.0, f"PNG height {h_mm:.2f} != {HEIGHT_MM}")

    hashes = {
        str(p.relative_to(ROOT)): sha256(p)
        for p in (pdf, svg, png, OUT_PREVIEW, OUT_SOURCE, OUT_SUMMARY, OUT_FAMILY, IN_CANONICAL)
    }
    write_qc(hashes, layout_stats)

    print(f"Wrote {OUT_SOURCE}")
    print(f"Wrote {OUT_SUMMARY}")
    print(f"Wrote {OUT_FAMILY}")
    print(f"Wrote {OUT_QC}")
    print(f"Wrote {OUT_BASE}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW}")


if __name__ == "__main__":
    main()
