#!/usr/bin/env python3
"""Corrected Figure 4C — chronicity class composition summary (16 frozen NCP families).

Layout-polish edition: consistent class blocks, separator clearance, Fig 4B typography,
mathematical |ρ_LAA| footer, 43 × 82 mm canvas.
"""
from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
OUT_SCRIPT_NAME = "fig4c_chronicity_class_summary"

CHRONICITY_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
CANONICAL_TSV = ROOT / "results" / "ncp_v3_canonical_display_order.tsv"

WIDTH_MM = 43.0
HEIGHT_MM = 82.0
MARGIN_MM = 1.5

STRONG_COLOR = "#2E2E2E"
WEAK_COLOR = "#7A7A7A"
STABLE_COLOR = "#C5C5C5"
FG_COLOR = "#1A1A1A"
MUTED_COLOR = "#555555"
FAINT_COLOR = "#777777"
RULE_COLOR = "#D0D0D0"
WHITE = "#FFFFFF"

CLASS_ORDER = ["Strong_trend", "Weak_trend", "Stable"]
CLASS_LABELS = {
    "Strong_trend": "Strong-trend",
    "Weak_trend": "Weak-trend",
    "Stable": "Stable",
}
CLASS_THRESHOLDS = {
    "Strong_trend": r"$|\rho_{\mathrm{LAA}}| \geq 0.4$",
    "Weak_trend": r"$0.2 \leq |\rho_{\mathrm{LAA}}| < 0.4$",
    "Stable": r"$|\rho_{\mathrm{LAA}}| < 0.2$",
}
CLASS_COLORS = {
    "Strong_trend": STRONG_COLOR,
    "Weak_trend": WEAK_COLOR,
    "Stable": STABLE_COLOR,
}
HEADER_COLORS = {
    "Strong_trend": STRONG_COLOR,
    "Weak_trend": WEAK_COLOR,
    "Stable": "#888888",
}
EXPECTED_SETS = {
    "Strong_trend": {2, 9},
    "Weak_trend": {3, 7, 8, 10, 13, 15},
    "Stable": {1, 4, 5, 6, 11, 12, 14, 16},
}

# Typography matched to Figure 4B at final physical size
PANEL_LETTER_PT = 10.5
PANEL_TITLE_PT = 8.8
SUBTITLE_PT = 6.5
CLASS_HEADING_PT = 6.5
THRESHOLD_PT = 5.5
FAMILY_PT = 6.0
BAR_COUNT_PT = 7.0
FOOTER_PT = 5.0
MIN_FONT_PT = 5.0

BBOX_TOL_MM = 0.2
OVERLAP_TOL_MM2 = 0.2
RULE_CLEARANCE_MM = 2.0

# Layout (mm from figure bottom; y increases upward)
# Title band near top; composition bar; three equal body blocks; footer band.
TITLE_TOP_MM = 80.5
TITLE_LETTER_Y = 80.0
TITLE_TEXT_Y = 80.0
SUBTITLE_Y = 76.2
SUBTITLE_Y2 = 73.8

BAR_H_MM = 7.0  # ~12.5% shorter than prior 8 mm
BAR_TOP_MM = 71.5
BAR_CENTER_MM = BAR_TOP_MM - BAR_H_MM / 2.0

# Three equal body blocks between bar bottom and footer top
BODY_TOP_MM = 63.0
FOOTER_TOP_MM = 14.0
N_BLOCKS = 3
BODY_SPAN_MM = BODY_TOP_MM - FOOTER_TOP_MM
BLOCK_PITCH_MM = BODY_SPAN_MM / N_BLOCKS  # ~16.33 mm per slot
BLOCK_CONTENT_H_MM = 12.0  # title + threshold + up to 2 family lines
BLOCK_TOPS = [BODY_TOP_MM - i * BLOCK_PITCH_MM for i in range(N_BLOCKS)]
# Separators sit in the residual gap below each block's content, with ≥2 mm
# clearance above the next block top (higher y = toward top of figure).
SEPARATOR_YS = [
    BLOCK_TOPS[i] - BLOCK_CONTENT_H_MM - (BLOCK_PITCH_MM - BLOCK_CONTENT_H_MM) / 2.0
    for i in range(N_BLOCKS - 1)
]

# Within-block vertical offsets from block top (positive = downward)
PAD_TOP_MM = 0.6
TITLE_OFFSET_MM = PAD_TOP_MM
THRESH_OFFSET_MM = 3.8
FAM1_OFFSET_MM = 7.4
FAM2_OFFSET_MM = 10.8

FOOTER_LINES = [
    (r"Classes use unrounded $|\rho_{\mathrm{LAA}}|$.", 11.5),
    ("Direction shown in B.", 8.0),
    ("Descriptive bins, not significance calls.", 4.5),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "font.size": 6.5,
        "axes.linewidth": 0.5,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def canonical_family_order() -> list[int]:
    order = (
        pd.read_csv(CANONICAL_TSV, sep="\t")
        .sort_values("family_display_order")
        .drop_duplicates("family_id")
    )
    return order["family_id"].astype(int).tolist()


def recompute_class(rho: float) -> str:
    abs_rho = abs(float(rho))
    if abs_rho >= 0.4:
        return "Strong_trend"
    if abs_rho >= 0.2:
        return "Weak_trend"
    return "Stable"


def load_and_assert() -> tuple[pd.DataFrame, pd.DataFrame]:
    cls = pd.read_csv(CHRONICITY_TSV, sep="\t")
    labels = pd.read_csv(LABELS_TSV, sep="\t")
    _ = pd.read_csv(CANONICAL_TSV, sep="\t")  # assert readable; used via canonical_family_order

    require(set(cls.columns) >= {"family_id", "rho_primary", "chronicity_class"},
            "chronicity table missing required columns")
    require(len(cls) == 16, f"expected 16 families, got {len(cls)}")
    require(cls["family_id"].nunique() == 16, "family_id values are not unique")

    cls = cls.copy()
    cls["recomputed_class"] = cls["rho_primary"].apply(recompute_class)
    cls["class_matches"] = cls["recomputed_class"] == cls["chronicity_class"]
    require(cls["class_matches"].all(),
            "recomputed classes do not match input chronicity_class")

    counts = cls["chronicity_class"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int)
    require(int(counts["Strong_trend"]) == 2, "Strong_trend count != 2")
    require(int(counts["Weak_trend"]) == 6, "Weak_trend count != 6")
    require(int(counts["Stable"]) == 8, "Stable count != 8")
    require(int(counts.sum()) == 16, "class totals must sum to 16")

    for class_name, expected in EXPECTED_SETS.items():
        actual = set(cls.loc[cls["chronicity_class"] == class_name, "family_id"])
        require(actual == expected, f"{class_name} set {actual} != expected {expected}")

    all_sets = [set(cls.loc[cls["chronicity_class"] == c, "family_id"]) for c in CLASS_ORDER]
    require(
        all_sets[0].isdisjoint(all_sets[1])
        and all_sets[0].isdisjoint(all_sets[2])
        and all_sets[1].isdisjoint(all_sets[2]),
        "chronicity classes are not disjoint",
    )
    require(all_sets[0] | all_sets[1] | all_sets[2] == set(cls["family_id"]),
            "classes do not partition all 16 families")
    return cls, labels


def build_source(cls: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    fam_order = canonical_family_order()
    labels = labels.set_index("family_id")
    order_lookup = {fid: rank for rank, fid in enumerate(fam_order, start=1)}

    rows = []
    for _, row in cls.iterrows():
        fid = int(row["family_id"])
        rho = float(row["rho_primary"])
        rows.append({
            "canonical_order": order_lookup[fid],
            "family_id": fid,
            "family_label": str(labels.loc[fid, "family_label"]),
            "short_label": str(labels.loc[fid, "short_label"]),
            "axis_group": str(labels.loc[fid, "axis_group"]),
            "rho_laa": rho,
            "abs_rho_laa": abs(rho),
            "chronicity_class": str(row["chronicity_class"]),
            "recomputed_class": str(row["recomputed_class"]),
            "class_matches_input": bool(row["class_matches"]),
        })
    source = pd.DataFrame(rows).sort_values("canonical_order")
    require(len(source) == 16, "source table must have 16 rows")
    require(source["class_matches_input"].all(), "source table contains mismatched recomputed class")
    return source


def build_summary(source: pd.DataFrame) -> pd.DataFrame:
    fam_order = canonical_family_order()
    order_lookup = {fid: rank for rank, fid in enumerate(fam_order, start=1)}

    rows = []
    for class_name in CLASS_ORDER:
        sub = source[source["chronicity_class"] == class_name].copy()
        sub["sort_key"] = sub["family_id"].map(order_lookup)
        sub = sub.sort_values("sort_key")
        members = [f"F{fid}" for fid in sub["family_id"].tolist()]
        rows.append({
            "chronicity_class": CLASS_LABELS[class_name],
            "count": len(sub),
            "fraction": f"{len(sub)}/16",
            "threshold": CLASS_THRESHOLDS[class_name],
            "members": " · ".join(members),
            "member_family_ids": ",".join(str(fid) for fid in sub["family_id"].tolist()),
        })
    return pd.DataFrame(rows)


def write_tables(source: pd.DataFrame, summary: pd.DataFrame) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    source_path = RESULTS / f"{OUT_SCRIPT_NAME}_source.tsv"
    summary_path = RESULTS / f"{OUT_SCRIPT_NAME}_summary.tsv"
    source.to_csv(source_path, sep="\t", index=False, float_format="%.6f")
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"wrote {source_path.relative_to(ROOT)}")
    print(f"wrote {summary_path.relative_to(ROOT)}")


def mm_axes(fig: plt.Figure, x_mm: float, y_mm: float, w_mm: float, h_mm: float) -> plt.Axes:
    return fig.add_axes([x_mm / WIDTH_MM, y_mm / HEIGHT_MM, w_mm / WIDTH_MM, h_mm / HEIGHT_MM])


def member_lines(class_name: str, members: str) -> tuple[str, str]:
    tokens = members.split(" · ")
    if class_name == "Weak_trend":
        return " · ".join(tokens[:3]), " · ".join(tokens[3:])
    if class_name == "Stable":
        return " · ".join(tokens[:4]), " · ".join(tokens[4:])
    return members, ""


def render(summary: pd.DataFrame) -> tuple[plt.Figure, list[Line2D]]:
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    fig.patch.set_facecolor("#FFFFFF")
    fig.set_dpi(300)
    usable_w = WIDTH_MM - 2 * MARGIN_MM

    ax = mm_axes(fig, MARGIN_MM, 0.0, usable_w, HEIGHT_MM)
    ax.set_xlim(0, usable_w)
    ax.set_ylim(0, HEIGHT_MM)
    ax.axis("off")

    # Title
    ax.text(0, TITLE_LETTER_Y, "C", fontsize=PANEL_LETTER_PT, fontweight="bold",
            va="top", color=FG_COLOR)
    ax.text(4.2, TITLE_TEXT_Y, "Chronicity classes", fontsize=PANEL_TITLE_PT,
            fontweight="bold", va="top", color=FG_COLOR)
    ax.text(0, SUBTITLE_Y, "LAA |ρ| bins across 16 frozen",
            fontsize=SUBTITLE_PT, va="top", color=MUTED_COLOR)
    ax.text(0, SUBTITLE_Y2, "NCP families",
            fontsize=SUBTITLE_PT, va="top", color=MUTED_COLOR)

    # Composition bar (data coords: x in [0, usable_w] mapped via secondary transform)
    ax_bar = mm_axes(fig, MARGIN_MM, 0.0, usable_w, HEIGHT_MM)
    ax_bar.set_xlim(0, 1)
    ax_bar.set_ylim(0, HEIGHT_MM)
    ax_bar.axis("off")

    x_left = 0.0
    for class_name in CLASS_ORDER:
        count = int(summary.loc[summary["chronicity_class"] == CLASS_LABELS[class_name], "count"].iloc[0])
        width = count / 16.0
        color = CLASS_COLORS[class_name]
        ax_bar.barh(
            BAR_CENTER_MM, width, left=x_left, height=BAR_H_MM,
            color=color, edgecolor=WHITE, linewidth=0.8,
        )
        text_color = WHITE if class_name in ("Strong_trend", "Weak_trend") else FG_COLOR
        ax_bar.text(
            x_left + width / 2, BAR_CENTER_MM, str(count),
            ha="center", va="center", fontsize=BAR_COUNT_PT,
            color=text_color, fontweight="bold",
        )
        x_left += width

    # Separator rules (stored for QC intersection checks)
    separator_lines: list[Line2D] = []
    for y in SEPARATOR_YS:
        (line,) = ax.plot(
            [0, usable_w], [y, y],
            color=RULE_COLOR, linewidth=0.6, clip_on=False, solid_capstyle="butt",
        )
        separator_lines.append(line)

    # Class blocks with consistent internal spacing
    for i, class_name in enumerate(CLASS_ORDER):
        top = BLOCK_TOPS[i]
        count = int(summary.loc[summary["chronicity_class"] == CLASS_LABELS[class_name], "count"].iloc[0])
        members = summary.loc[summary["chronicity_class"] == CLASS_LABELS[class_name], "members"].iloc[0]
        line1, line2 = member_lines(class_name, members)

        ax.text(
            0, top - TITLE_OFFSET_MM,
            f"{CLASS_LABELS[class_name]}  {count}/16",
            fontsize=CLASS_HEADING_PT, fontweight="bold", va="top",
            color=HEADER_COLORS[class_name],
        )
        ax.text(
            0, top - THRESH_OFFSET_MM,
            CLASS_THRESHOLDS[class_name],
            fontsize=THRESHOLD_PT, va="top", color=MUTED_COLOR,
        )
        ax.text(
            0, top - FAM1_OFFSET_MM,
            line1,
            fontsize=FAMILY_PT, va="top", color=FG_COLOR,
        )
        if line2:
            ax.text(
                0, top - FAM2_OFFSET_MM,
                line2,
                fontsize=FAMILY_PT, va="top", color=FG_COLOR,
            )

    # Footer
    for text, y in FOOTER_LINES:
        ax.text(0, y, text, fontsize=FOOTER_PT, va="top", color=FAINT_COLOR)

    return fig, separator_lines


def text_bbox_mm(artist: mpl.text.Text, fig: plt.Figure) -> tuple[float, float, float, float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    return (
        float(bbox.x0) * WIDTH_MM,
        float(bbox.y0) * HEIGHT_MM,
        float(bbox.x1) * WIDTH_MM,
        float(bbox.y1) * HEIGHT_MM,
    )


def collect_visible_text(fig: plt.Figure) -> list[tuple[str, mpl.text.Text, tuple[float, float, float, float]]]:
    entries: list[tuple[str, mpl.text.Text, tuple[float, float, float, float]]] = []
    for ax in fig.axes:
        for child in ax.get_children():
            if not isinstance(child, mpl.text.Text):
                continue
            label = child.get_text().strip()
            if not label:
                continue
            entries.append((label, child, text_bbox_mm(child, fig)))
    return entries


def bbox_overlap(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> float:
    x0 = max(a[0], b[0])
    y0 = max(a[1], b[1])
    x1 = min(a[2], b[2])
    y1 = min(a[3], b[3])
    if x1 <= x0 or y1 <= y0:
        return 0.0
    return (x1 - x0) * (y1 - y0)


def line_intersects_bbox(y_mm: float, bbox: tuple[float, float, float, float], tol_mm: float = 0.15) -> bool:
    """True if a horizontal separator at y_mm intersects the text bbox (with tolerance)."""
    return (bbox[1] - tol_mm) <= y_mm <= (bbox[3] + tol_mm)


def validate_layout(fig: plt.Figure, separator_lines: list[Line2D]) -> dict[str, float | int]:
    entries = collect_visible_text(fig)
    require(entries, "layout QC found no visible text")

    min_font = min(child.get_fontsize() for _, child, _ in entries)
    require(min_font >= MIN_FONT_PT - 0.05, f"minimum font size {min_font:.2f} pt < {MIN_FONT_PT} pt")

    violations = 0
    for label, _, bbox in entries:
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

    overlaps = 0
    for idx, (label_a, _, bbox_a) in enumerate(entries):
        for label_b, _, bbox_b in entries[idx + 1:]:
            overlap = bbox_overlap(bbox_a, bbox_b)
            if overlap > OVERLAP_TOL_MM2:
                overlaps += 1
                require(False, f"overlapping text: {label_a!r} vs {label_b!r} ({overlap:.2f} mm²)")

    # Hard assertion: no separator intersects any text bbox
    rule_hits = 0
    for line in separator_lines:
        ydata = line.get_ydata()
        y_mm = float(ydata[0])
        for label, _, bbox in entries:
            if line_intersects_bbox(y_mm, bbox):
                rule_hits += 1
                require(False, f"separator at y={y_mm:.2f} mm intersects text {label!r}")

    # Clearance: separator must sit above the next block top by ≥ RULE_CLEARANCE_MM
    for sep_y, next_top in zip(SEPARATOR_YS, BLOCK_TOPS[1:]):
        clearance_to_next = sep_y - next_top
        require(
            clearance_to_next >= RULE_CLEARANCE_MM - 0.05,
            f"separator clearance {clearance_to_next:.2f} mm < {RULE_CLEARANCE_MM} mm to next block top",
        )

    return {
        "n_text_objects": len(entries),
        "min_font_pt": float(min_font),
        "n_bbox_violations": violations,
        "n_overlaps": overlaps,
        "n_rule_text_hits": rule_hits,
        "separator_ys": SEPARATOR_YS,
        "block_tops": BLOCK_TOPS,
    }


def save_outputs(fig: plt.Figure) -> dict[str, Path]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    paths = {
        "pdf": FIGURES / f"{OUT_SCRIPT_NAME}.pdf",
        "svg": FIGURES / f"{OUT_SCRIPT_NAME}.svg",
        "png": FIGURES / f"{OUT_SCRIPT_NAME}.png",
        "preview": FIGURES / f"{OUT_SCRIPT_NAME}_preview.png",
    }
    for key, path in paths.items():
        dpi = 600 if key == "png" else (150 if key == "preview" else None)
        kwargs = {"bbox_inches": None, "pad_inches": 0.02}
        if dpi is not None:
            kwargs["dpi"] = dpi
        fig.savefig(path, **kwargs)
        print(f"wrote {path.relative_to(ROOT)}")
    plt.close(fig)
    return paths


def write_qc(
    summary: pd.DataFrame,
    layout_qc: dict[str, float | int],
    output_paths: dict[str, Path],
) -> None:
    qc_path = RESULTS / f"{OUT_SCRIPT_NAME}_qc.md"
    input_hashes = {
        str(CHRONICITY_TSV.relative_to(ROOT)): sha256(CHRONICITY_TSV),
        str(LABELS_TSV.relative_to(ROOT)): sha256(LABELS_TSV),
        str(CANONICAL_TSV.relative_to(ROOT)): sha256(CANONICAL_TSV),
    }
    output_hashes = {str(p.relative_to(ROOT)): sha256(p) for p in output_paths.values()}
    source_path = RESULTS / f"{OUT_SCRIPT_NAME}_source.tsv"
    summary_path = RESULTS / f"{OUT_SCRIPT_NAME}_summary.tsv"
    output_hashes[str(source_path.relative_to(ROOT))] = sha256(source_path)
    output_hashes[str(summary_path.relative_to(ROOT))] = sha256(summary_path)

    lines = [
        "# Figure 4C chronicity class summary QC",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "Status: **PASS**",
        "",
        "## Scientific contract",
        "- 16 unique frozen NCP families",
        f"- Strong-trend: {int(summary.iloc[0]['count'])} ({summary.iloc[0]['fraction']})",
        f"- Weak-trend: {int(summary.iloc[1]['count'])} ({summary.iloc[1]['fraction']})",
        f"- Stable: {int(summary.iloc[2]['count'])} ({summary.iloc[2]['fraction']})",
        f"- Strong set: {{{summary.iloc[0]['member_family_ids']}}}",
        f"- Weak set: {{{summary.iloc[1]['member_family_ids']}}}",
        f"- Stable set: {{{summary.iloc[2]['member_family_ids']}}}",
        "- Sets are disjoint and complete; totals sum to 16",
        "- Classes recomputed from full-precision |rho_laa|; all match input chronicity_class",
        "- No p, CI, q, LAPV or regional-sign field used for class assignment or display",
        "- No legacy Figure 4 tables read",
        "",
        "## Typography (Fig 4B matched)",
        f"- panel letter: {PANEL_LETTER_PT} pt",
        f"- panel title: {PANEL_TITLE_PT} pt",
        f"- subtitle: {SUBTITLE_PT} pt",
        f"- class headings: {CLASS_HEADING_PT} pt",
        f"- footer: {FOOTER_PT} pt",
        "",
        "## Layout QC",
        f"- visible text objects: {layout_qc['n_text_objects']}",
        f"- minimum font size: {layout_qc['min_font_pt']:.2f} pt",
        f"- bbox violations: {layout_qc['n_bbox_violations']}",
        f"- text overlaps: {layout_qc['n_overlaps']}",
        f"- separator–text intersections: {layout_qc['n_rule_text_hits']}",
        f"- separator y (mm): {layout_qc['separator_ys']}",
        f"- block tops (mm): {layout_qc['block_tops']}",
        f"- panel geometry: {WIDTH_MM:.1f} mm × {HEIGHT_MM:.1f} mm",
        f"- composition bar height: {BAR_H_MM:.1f} mm",
        "",
        "## Allowed inputs",
        *[f"- `{path}` — `{digest}`" for path, digest in input_hashes.items()],
        "",
        "## Outputs",
        *[f"- `{path}` — `{digest}`" for path, digest in output_hashes.items()],
        "",
    ]
    qc_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {qc_path.relative_to(ROOT)}")


def main() -> None:
    set_style()
    cls, labels = load_and_assert()
    source = build_source(cls, labels)
    summary = build_summary(source)
    write_tables(source, summary)
    fig, separator_lines = render(summary)
    layout_qc = validate_layout(fig, separator_lines)
    output_paths = save_outputs(fig)
    write_qc(summary, layout_qc, output_paths)


if __name__ == "__main__":
    main()
