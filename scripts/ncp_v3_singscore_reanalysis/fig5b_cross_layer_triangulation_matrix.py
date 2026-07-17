#!/usr/bin/env python3
"""Figure 5B — family-level cross-layer triangulation matrix.

Reads only the two Figure 5 authoritative tables. Does not modify 5A/C/D
or formal submission routing. Does not open PriorityScore or legacy Fig 5.
"""
from __future__ import annotations

import hashlib
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"

IN_MANIFEST = RESULTS / "figure5_authoritative_input_manifest.tsv"
IN_CANONICAL = RESULTS / "figure5_triangulation_canonical.tsv"

OUT_BASE = FIGURES / "fig5b_cross_layer_triangulation_matrix"
OUT_PREVIEW = FIGURES / "fig5b_cross_layer_triangulation_matrix_preview.png"
OUT_SOURCE = RESULTS / "fig5b_cross_layer_triangulation_source.tsv"
OUT_SUMMARY = RESULTS / "fig5b_cross_layer_triangulation_summary.tsv"
OUT_QC = RESULTS / "fig5b_cross_layer_triangulation_qc.md"

FORBIDDEN_READ_PATTERNS = (
    re.compile(r"PriorityScore"),
    re.compile(r"results/fig5b_v3_triangulation\.tsv"),
    re.compile(r"fig2e_cross_cohort_consistency"),
    re.compile(r"fig2d_family_effects"),
    re.compile(r"fig3_v3_carrier_fidelity"),
    re.compile(r"fig4_v3_chronicity"),
    re.compile(r"fig1f_v3_HCA"),
    re.compile(r"geomx", re.IGNORECASE),
    re.compile(r"roi_scores_v3"),
)

EXPECTED_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]
EXPECTED_PATH_A = {10, 15}
EXPECTED_DOUBLE = {9}
EXPECTED_SILENT = {4, 5, 6}
EXPECTED_COUNTS = {"triple": 2, "double": 1, "single": 10, "silent": 3}

MM = 1 / 25.4
WIDTH_MM = 178.0
HEIGHT_MM = 96.0
MARGIN_MM = 3.0
TITLE_H_MM = 12.5
LEGEND_H_MM = 20.0
TABLE_TOP = HEIGHT_MM - MARGIN_MM - TITLE_H_MM
TABLE_BOTTOM = MARGIN_MM + LEGEND_H_MM
TABLE_H = TABLE_TOP - TABLE_BOTTOM

BODY_PT = 6.2
HEADER_PT = 6.0
GATE_CRIT_PT = 5.2
GROUP_PT = 6.0
NOTE_PT = 6.0
MIN_FONT_PT = 5.2
PANEL_LETTER_PT = 10.8
TITLE_TEXT_PT = 8.6
SUBTITLE_PT = 6.3

_STANDALONE = {
    "WIDTH_MM": 178.0,
    "HEIGHT_MM": 96.0,
    "MARGIN_MM": 3.0,
    "TITLE_H_MM": 12.5,
    "LEGEND_H_MM": 20.0,
    "PANEL_LETTER_PT": 10.8,
    "TITLE_TEXT_PT": 8.6,
    "SUBTITLE_PT": 6.3,
    "BODY_PT": 6.2,
    "HEADER_PT": 6.0,
    "GATE_CRIT_PT": 5.2,
    "GROUP_PT": 6.0,
    "NOTE_PT": 6.0,
    "MIN_FONT_PT": 5.2,
}
_COMPOSITE = {
    "WIDTH_MM": 178.0,
    "HEIGHT_MM": 96.0,  # +4 mm vs prior composite so body 6.7 pt clears 16 rows
    "MARGIN_MM": 3.0,
    "TITLE_H_MM": 10.5,  # free body row height for type bump
    "LEGEND_H_MM": 20.0,  # full two-line legend + footer note
    # One-step type bump for composite readability
    "PANEL_LETTER_PT": 11.0,
    "TITLE_TEXT_PT": 9.0,
    "SUBTITLE_PT": 6.8,
    "BODY_PT": 6.7,
    "HEADER_PT": 6.5,
    "GATE_CRIT_PT": 5.6,
    "GROUP_PT": 6.5,
    "NOTE_PT": 6.5,
    "MIN_FONT_PT": 5.6,
}

LAYOUT_MODE = "standalone"


def apply_layout(mode: str = "standalone") -> None:
    """Switch between standalone and Figure 5 composite canvas geometry."""
    global WIDTH_MM, HEIGHT_MM, MARGIN_MM, TITLE_H_MM, LEGEND_H_MM
    global TABLE_TOP, TABLE_BOTTOM, TABLE_H, LAYOUT_MODE
    global PANEL_LETTER_PT, TITLE_TEXT_PT, SUBTITLE_PT
    global BODY_PT, HEADER_PT, GATE_CRIT_PT, GROUP_PT, NOTE_PT, MIN_FONT_PT
    LAYOUT_MODE = mode
    cfg = _COMPOSITE if mode == "composite" else _STANDALONE
    WIDTH_MM = cfg["WIDTH_MM"]
    HEIGHT_MM = cfg["HEIGHT_MM"]
    MARGIN_MM = cfg["MARGIN_MM"]
    TITLE_H_MM = cfg["TITLE_H_MM"]
    LEGEND_H_MM = cfg["LEGEND_H_MM"]
    PANEL_LETTER_PT = cfg["PANEL_LETTER_PT"]
    TITLE_TEXT_PT = cfg["TITLE_TEXT_PT"]
    SUBTITLE_PT = cfg["SUBTITLE_PT"]
    BODY_PT = cfg["BODY_PT"]
    HEADER_PT = cfg["HEADER_PT"]
    GATE_CRIT_PT = cfg["GATE_CRIT_PT"]
    GROUP_PT = cfg["GROUP_PT"]
    NOTE_PT = cfg["NOTE_PT"]
    MIN_FONT_PT = cfg["MIN_FONT_PT"]
    TABLE_TOP = HEIGHT_MM - MARGIN_MM - TITLE_H_MM
    TABLE_BOTTOM = MARGIN_MM + LEGEND_H_MM
    TABLE_H = TABLE_TOP - TABLE_BOTTOM

FG = "#1A1A1A"
MUTED = "#595959"
FAINT = "#7A7A7A"
GATE_CRIT_COLOR = "#8A8A8A"
WHITE = "#FFFFFF"
ROW_LINE = "#F3F3F3"
LA_RA_LINE = "#8A8A8A"
GROUP_LINE = "#BDBDBD"
GROUP_RULE = "#999999"
PATH_A_EDGE = "#8A8A8A"
CONTEXT_BG = "#EFEFEF"
GATES_BG = "#FFFFFF"
TRIAGE_BG = "#EFEFEF"
STATE_ALPHA = 0.88  # context-only desaturation

LA_COLOR = "#9C5550"
RA_COLOR = "#4E6B89"

STATE_COLORS = {
    "Amplified": "#D58A45",
    "Preserved": "#659478",
    "Attenuated": "#D8C27E",
    "Reversed": "#B5A3CC",
}
STATE_ABBREV = {
    "Amplified": "Amp",
    "Preserved": "Pres",
    "Attenuated": "Att",
    "Reversed": "Rev",
}
STATE_LEGEND = {
    "Amplified": "Amplified",
    "Preserved": "Preserved",
    "Attenuated": "Attenuated",
    "Reversed": "Reversed",
}

CONS_SHAPE = {
    "High": ("circle", "#2E7B5E", "H"),
    "Moderate": ("triangle", "#D4A017", "M"),
    "Low": ("square", "#9A9A9A", "L"),
}

SN_COLORS = {
    "Fidelity": "#8FAF7E",
    "Relay": "#A898B8",
    "Extension": "#7FA8AA",
    "Below_detection": "#DBDBDB",  # ~12% lighter than #CECECE
}
SN_LABEL = {
    "Fidelity": "Fidelity",
    "Relay": "Relay",
    "Extension": "Extension",
    "Below_detection": "Below",
}
SN_LEGEND = {
    "Fidelity": "Fidelity",
    "Relay": "Relay",
    "Extension": "Extension",
    "Below_detection": "Below-detection",
}

DUR_COLORS = {
    "Strong_trend": "#2E2E2E",
    "Weak_trend": "#7A7A7A",
    "Stable": "#CECECE",  # ~7% lighter than #C5C5C5
}
DUR_LABEL = {
    "Strong_trend": "Strong",
    "Weak_trend": "Weak",
    "Stable": "Stable",
}

HCA_ABBREV = {
    "endothelial": "Endo",
    "atrialCM": "aCM",
    "adipocyte": "Adipo",
    "fibroblast": "Fibro",
    "lymphocyte": "Lymph",
    "Mixed/Unresolved": "Mixed",
}

# Publication-facing label corrections (canonical values may already match).
LABEL_FIXES = {
    2: "Mito Ca²⁺ transport",
    4: "Mito proteostasis / FA–BCAA",
    10: "GAG disease / TGF-β",
    15: "NTRK2 / TrkB",
}

GATE_CRITERIA = {
    "cons": "gate: High",
    "sn": "gate: non-Below-detection",
    "dur": "gate: non-Stable",
}

COL_SPEC = [
    ("strip", 3.5),
    ("fid", 7.0),
    ("label", 34.0),  # room for longest family label at composite body 6.7 pt
    ("hca", 17.0),
    ("state", 22.0),
    ("cons", 18.0),
    ("sn", 21.0),
    ("dur", 20.0),
    ("supp", 15.0),
    ("patha", 14.5),
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def support_tier(count: int) -> str:
    return {3: "triple", 2: "double", 1: "single", 0: "silent"}[count]


def as_bool(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value) == "True"


def assert_no_forbidden_reads() -> None:
    text = Path(__file__).read_text(encoding="utf-8")
    for line in text.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("#"):
            continue
        # Match actual pandas reads only (ignore this guard itself).
        if not re.search(r"\bpd\.read_csv\s*\(", line):
            continue
        require(
            "IN_MANIFEST" in line or "IN_CANONICAL" in line,
            f"Unexpected read_csv: {line.strip()}",
        )
        for pattern in FORBIDDEN_READ_PATTERNS:
            require(pattern.search(line) is None, f"Forbidden read pattern in: {line.strip()}")


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": BODY_PT,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": WHITE,
            "savefig.facecolor": WHITE,
        }
    )


def load_and_assert() -> tuple[pd.DataFrame, pd.DataFrame]:
    require(IN_MANIFEST.is_file(), f"Missing manifest: {IN_MANIFEST}")
    require(IN_CANONICAL.is_file(), f"Missing canonical: {IN_CANONICAL}")
    manifest = pd.read_csv(IN_MANIFEST, sep="\t")
    df = pd.read_csv(IN_CANONICAL, sep="\t")
    df = df.sort_values("canonical_order").reset_index(drop=True)

    require(list(df["family_id"].astype(int)) == EXPECTED_ORDER, "Canonical order drifted")
    require("PriorityScore" not in df.columns, "PriorityScore column present")
    joined_cols = ",".join(df.columns)
    for banned in ("GeoMx", "geomx", "genetics", "PriorityScore"):
        require(banned not in joined_cols, f"Banned column marker: {banned}")

    recomputed_support = []
    recomputed_path_a = []
    for row in df.itertuples(index=False):
        bulk_gate2 = str(row.bulk_consistency) == "High"
        sn_gate2 = str(row.snRNA_class) != "Below_detection"
        duration_gate2 = str(row.duration_class) != "Stable"
        require(as_bool(row.bulk_gate) == bulk_gate2, f"F{row.family_id} bulk_gate mismatch")
        require(as_bool(row.sn_gate) == sn_gate2, f"F{row.family_id} sn_gate mismatch")
        require(as_bool(row.duration_gate) == duration_gate2, f"F{row.family_id} duration_gate mismatch")
        support = int(bulk_gate2) + int(sn_gate2) + int(duration_gate2)
        require(int(row.support_count) == support, f"F{row.family_id} support_count mismatch")
        path_a = support == 3
        require(as_bool(row.Path_A) == path_a, f"F{row.family_id} Path_A mismatch")
        recomputed_support.append(support)
        recomputed_path_a.append(path_a)

    df = df.copy()
    df["support_recomputed"] = recomputed_support
    df["Path_A_recomputed"] = recomputed_path_a

    path_a_ids = set(df.loc[df["Path_A_recomputed"], "family_id"].astype(int))
    require(path_a_ids == EXPECTED_PATH_A, f"Path A {path_a_ids} != {EXPECTED_PATH_A}")
    for fid in EXPECTED_PATH_A:
        row = df.loc[df["family_id"] == fid].iloc[0]
        require(int(row["support_recomputed"]) == 3, f"F{fid} not 3/3")
    for fid in EXPECTED_DOUBLE:
        row = df.loc[df["family_id"] == fid].iloc[0]
        require(int(row["support_recomputed"]) == 2, f"F{fid} not 2/3")
        require(not bool(row["Path_A_recomputed"]), f"F{fid} wrongly Path A")
    for fid in EXPECTED_SILENT:
        row = df.loc[df["family_id"] == fid].iloc[0]
        require(int(row["support_recomputed"]) == 0, f"F{fid} not 0/3")

    tiers = Counter(support_tier(int(c)) for c in df["support_recomputed"])
    for tier, n in EXPECTED_COUNTS.items():
        require(tiers[tier] == n, f"{tier} count {tiers[tier]} != {n}")

    return manifest, df


def build_source(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in df.itertuples(index=False):
        hca = str(row.hca_carrier)
        state = str(row.bulk_dominant_state)
        sn = str(row.snRNA_class)
        dur = str(row.duration_class)
        cons = str(row.bulk_consistency)
        require(hca in HCA_ABBREV, f"Unknown HCA carrier: {hca}")
        require(state in STATE_ABBREV, f"Unknown dominant state: {state}")
        require(cons in CONS_SHAPE, f"Unknown consistency: {cons}")
        require(sn in SN_LABEL, f"Unknown snRNA class: {sn}")
        require(dur in DUR_LABEL, f"Unknown duration class: {dur}")
        support = int(row.support_recomputed)
        fid = int(row.family_id)
        label = LABEL_FIXES.get(fid, str(row.family_label))
        rows.append(
            {
                "family_id": fid,
                "canonical_order": int(row.canonical_order),
                "family_label": label,
                "family_direction": str(row.family_direction),
                "hca_carrier": hca,
                "hca_abbrev": HCA_ABBREV[hca],
                "bulk_dominant_state": state,
                "state_abbrev": STATE_ABBREV[state],
                "bulk_consistency": cons,
                "consistency_glyph": CONS_SHAPE[cons][2],
                "snRNA_class": sn,
                "snRNA_display": SN_LABEL[sn],
                "duration_class": dur,
                "duration_display": DUR_LABEL[dur],
                "support_count": support,
                "support_display": f"{support}/3",
                "Path_A": bool(row.Path_A_recomputed),
                "path_a_display": "Path A" if row.Path_A_recomputed else "–",
            }
        )
    return pd.DataFrame(rows)


def column_geometry() -> dict:
    content_w = WIDTH_MM - 2 * MARGIN_MM
    widths = [w for _, w in COL_SPEC]
    require(abs(sum(widths) - content_w) < 0.05, f"Column widths sum {sum(widths)} != {content_w}")
    xs = []
    x = MARGIN_MM
    for name, w in COL_SPEC:
        xs.append((name, x, w, x + w / 2))
        x += w
    return {"cols": xs, "content_w": content_w}


def draw_consistency_glyph(
    ax: plt.Axes,
    cx: float,
    cy: float,
    consistency: str,
    *,
    with_letter: bool = True,
    glyph_x: float | None = None,
) -> None:
    shape, color, letter = CONS_SHAPE[consistency]
    gx = glyph_x if glyph_x is not None else (cx - 2.0 if with_letter else cx)
    size = 1.15
    if shape == "circle":
        ax.add_patch(Circle((gx, cy), size * 0.85, facecolor=color, edgecolor="none", zorder=3))
    elif shape == "triangle":
        h = size * 1.55
        ax.add_patch(
            Polygon(
                [
                    (gx, cy + h * 0.55),
                    (gx - h * 0.55, cy - h * 0.45),
                    (gx + h * 0.55, cy - h * 0.45),
                ],
                closed=True,
                facecolor=color,
                edgecolor="none",
                zorder=3,
            )
        )
    else:
        s = size * 1.35
        ax.add_patch(Rectangle((gx - s / 2, cy - s / 2), s, s, facecolor=color, edgecolor="none", zorder=3))
    if with_letter:
        ax.text(cx + 1.3, cy, letter, ha="center", va="center", fontsize=BODY_PT, color=FG, fontweight="bold", zorder=4)


def support_fontweight(count: int) -> str:
    if count >= 3:
        return "bold"
    if count == 2:
        return "demibold"
    return "normal"


def render(source: pd.DataFrame) -> tuple[plt.Figure, plt.Axes]:
    set_style()
    geom = column_geometry()
    cols = {name: (x0, w, xc) for name, x0, w, xc in geom["cols"]}

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM), facecolor=WHITE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, WIDTH_MM)
    ax.set_ylim(0, HEIGHT_MM)
    ax.axis("off")
    ax.set_facecolor(WHITE)

    y_title = HEIGHT_MM - MARGIN_MM - 0.2
    title_step = 3.3 if LAYOUT_MODE == "composite" else 3.6
    note_step = 6.1 if LAYOUT_MODE == "composite" else 6.6
    ax.text(MARGIN_MM, y_title, "B", ha="left", va="top", fontsize=PANEL_LETTER_PT, fontweight="bold", color=FG)
    ax.text(
        MARGIN_MM + 6.5,
        y_title,
        "Family-level cross-layer triangulation matrix",
        ha="left",
        va="top",
        fontsize=TITLE_TEXT_PT,
        fontweight="bold",
        color=FG,
    )
    ax.text(
        MARGIN_MM + 6.5,
        y_title - title_step,
        "Context annotations and three Path A gates in the frozen Fig 1D order",
        ha="left",
        va="top",
        fontsize=SUBTITLE_PT,
        color=FAINT,
    )
    ax.text(
        MARGIN_MM + 6.5,
        y_title - note_step,
        "HCA carrier and bulk dominant state are context only; support count uses bulk, snRNA and duration gates.",
        ha="left",
        va="top",
        fontsize=NOTE_PT,
        color=FAINT,
    )

    n_rows = len(source)
    # Three-level header: group text -> group rule -> column headers -> gate defs -> column rule.
    # Composite type bump: taller header than standalone, but leave >= ~2.4 mm/row for body 6.7 pt.
    header_h = 19.5 if LAYOUT_MODE == "composite" else 18.5
    body_h = TABLE_H - header_h
    row_h = body_h / n_rows
    header_y0 = TABLE_TOP - header_h
    body_y0 = TABLE_BOTTOM

    hca_x0, _, _ = cols["hca"]
    state_x0, state_w, _ = cols["state"]
    cons_x0, _, _ = cols["cons"]
    dur_x0, dur_w, _ = cols["dur"]
    supp_x0, _, _ = cols["supp"]
    patha_x0, patha_w, _ = cols["patha"]
    label_x0, _, _ = cols["label"]
    strip_x0, strip_w, strip_xc = cols["strip"]

    # Group backgrounds spanning header + body
    context_w = (state_x0 + state_w) - hca_x0
    gates_w = (dur_x0 + dur_w) - cons_x0
    triage_w = (patha_x0 + patha_w) - supp_x0
    band_top = TABLE_TOP
    band_h = TABLE_TOP - body_y0
    ax.add_patch(Rectangle((hca_x0, body_y0), context_w, band_h, facecolor=CONTEXT_BG, edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((cons_x0, body_y0), gates_w, band_h, facecolor=GATES_BG, edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((supp_x0, body_y0), triage_w, band_h, facecolor=TRIAGE_BG, edgecolor="none", zorder=0))

    # Stronger vertical rules between groups (full height of table band)
    ax.plot([hca_x0, hca_x0], [body_y0, band_top], color=GROUP_RULE, lw=1.05, zorder=1)
    ax.plot([cons_x0, cons_x0], [body_y0, band_top], color=GROUP_RULE, lw=1.05, zorder=1)
    ax.plot([supp_x0, supp_x0], [body_y0, band_top], color=GROUP_RULE, lw=1.05, zorder=1)
    ax.plot([patha_x0 + patha_w, patha_x0 + patha_w], [body_y0, band_top], color=GROUP_RULE, lw=0.55, zorder=1)

    # LEVEL 1 — group headers (fixed baseline, independent from column headers)
    group_y = TABLE_TOP - 2.5
    ax.text((hca_x0 + state_x0 + state_w) / 2, group_y, "CONTEXT", ha="center", va="center", fontsize=GROUP_PT, color=MUTED, fontweight="bold", zorder=2)
    ax.text((cons_x0 + dur_x0 + dur_w) / 2, group_y, "PATH A GATES", ha="center", va="center", fontsize=GROUP_PT, color=MUTED, fontweight="bold", zorder=2)
    ax.text((supp_x0 + patha_x0 + patha_w) / 2, group_y, "TRIAGE", ha="center", va="center", fontsize=GROUP_PT, color=MUTED, fontweight="bold", zorder=2)

    group_rule_y = TABLE_TOP - 5.5
    ax.plot([hca_x0 + 0.4, state_x0 + state_w - 0.4], [group_rule_y, group_rule_y], color=GROUP_LINE, lw=0.7, zorder=2)
    ax.plot([cons_x0 + 0.4, dur_x0 + dur_w - 0.4], [group_rule_y, group_rule_y], color=GROUP_LINE, lw=0.7, zorder=2)
    ax.plot([supp_x0 + 0.4, patha_x0 + patha_w - 0.4], [group_rule_y, group_rule_y], color=GROUP_LINE, lw=0.7, zorder=2)

    # LEVEL 2 — column headers (fixed baseline, two-line where needed)
    col_headers = {
        "hca": "HCA LA\ncarrier",
        "state": "Bulk dominant\nstate",
        "cons": "Bulk\nconsistency",
        "sn": "LA snRNA\nclass",
        "dur": "LAA duration\nclass",
        "supp": "Support\ncount",
        "patha": "Path A",
    }
    col_header_y = TABLE_TOP - (9.6 if LAYOUT_MODE == "composite" else 9.5)
    for key, text in col_headers.items():
        _, _, xc = cols[key]
        ax.text(
            xc,
            col_header_y,
            text,
            ha="center",
            va="center",
            fontsize=HEADER_PT,
            color=MUTED,
            linespacing=0.95 if LAYOUT_MODE == "composite" else 1.0,
            zorder=2,
        )

    # LEVEL 3 — gate definitions (fixed baseline, only under gate columns)
    gate_def_y = TABLE_TOP - (15.6 if LAYOUT_MODE == "composite" else 14.5)
    gate_defs = {
        "cons": "gate: High",
        "sn": "gate: non-\nBelow-detection",
        "dur": "gate: non-\nStable",
    }
    for key, crit in gate_defs.items():
        _, _, xc = cols[key]
        ax.text(
            xc,
            gate_def_y,
            crit,
            ha="center",
            va="center",
            fontsize=GATE_CRIT_PT,
            color=GATE_CRIT_COLOR,
            linespacing=0.95 if LAYOUT_MODE == "composite" else 1.0,
            zorder=2,
        )

    # Column rule at bottom of header region (raised slightly so first-row text stays below it).
    header_rule_y = header_y0 + 0.3
    ax.plot([MARGIN_MM, WIDTH_MM - MARGIN_MM], [header_rule_y, header_rule_y], color="#D0D0D0", lw=0.55, zorder=2)

    # Record rule segments for layout validation.
    rule_segments = [
        (hca_x0, body_y0, hca_x0, band_top, "vertical"),
        (cons_x0, body_y0, cons_x0, band_top, "vertical"),
        (supp_x0, body_y0, supp_x0, band_top, "vertical"),
        (patha_x0 + patha_w, body_y0, patha_x0 + patha_w, band_top, "vertical"),
        (hca_x0 + 0.4, group_rule_y, state_x0 + state_w - 0.4, group_rule_y, "horizontal"),
        (cons_x0 + 0.4, group_rule_y, dur_x0 + dur_w - 0.4, group_rule_y, "horizontal"),
        (supp_x0 + 0.4, group_rule_y, patha_x0 + patha_w - 0.4, group_rule_y, "horizontal"),
        (MARGIN_MM, header_rule_y, WIDTH_MM - MARGIN_MM, header_rule_y, "horizontal"),
    ]

    # Continuous LA / RA direction strips
    n_la = int((source["family_direction"] == "LA-biased").sum())
    n_ra = n_rows - n_la
    require(n_la == 11 and n_ra == 5, f"Expected 11 LA + 5 RA rows, got {n_la}/{n_ra}")
    la_top = header_y0
    la_bot = header_y0 - n_la * row_h
    ra_top = la_bot
    ra_bot = body_y0
    strip_draw_w = strip_w * 0.72
    ax.add_patch(
        Rectangle((strip_x0, la_bot), strip_draw_w, la_top - la_bot, facecolor=LA_COLOR, edgecolor="none", zorder=2)
    )
    ax.add_patch(
        Rectangle((strip_x0, ra_bot), strip_draw_w, ra_top - ra_bot, facecolor=RA_COLOR, edgecolor="none", zorder=2)
    )
    ax.text(
        strip_x0 + strip_draw_w / 2,
        (la_top + la_bot) / 2,
        "LA",
        ha="center",
        va="center",
        fontsize=6.0,
        color=WHITE,
        fontweight="bold",
        rotation=90,
        zorder=3,
    )
    ax.text(
        strip_x0 + strip_draw_w / 2,
        (ra_top + ra_bot) / 2,
        "RA",
        ha="center",
        va="center",
        fontsize=6.0,
        color=WHITE,
        fontweight="bold",
        rotation=90,
        zorder=3,
    )

    for i, row in enumerate(source.itertuples(index=False)):
        y_top = header_y0 - i * row_h
        y_bot = y_top - row_h
        cy = (y_top + y_bot) / 2

        if i > 0:
            prev = source.iloc[i - 1]
            stronger = str(prev["family_direction"]) == "LA-biased" and str(row.family_direction) == "RA-biased"
            if stronger:
                ax.plot([label_x0, WIDTH_MM - MARGIN_MM], [y_top, y_top], color=LA_RA_LINE, lw=1.15, zorder=2)
            else:
                ax.plot([label_x0, WIDTH_MM - MARGIN_MM], [y_top, y_top], color=ROW_LINE, lw=0.28, zorder=2)

        _, _, fid_xc = cols["fid"]
        ax.text(fid_xc, cy, f"F{row.family_id}", ha="center", va="center", fontsize=BODY_PT, color=FG, fontweight="bold", zorder=3)
        ax.text(label_x0 + 0.4, cy, row.family_label, ha="left", va="center", fontsize=BODY_PT, color=FG, zorder=3)

        _, _, hca_xc = cols["hca"]
        ax.text(hca_xc, cy, row.hca_abbrev, ha="center", va="center", fontsize=BODY_PT, color=MUTED, zorder=3)

        _, _, state_xc = cols["state"]
        state_color = STATE_COLORS[row.bulk_dominant_state]
        chip_w, chip_h = 9.5, min(3.2, row_h * 0.72)
        ax.add_patch(
            FancyBboxPatch(
                (state_xc - chip_w / 2, cy - chip_h / 2),
                chip_w,
                chip_h,
                boxstyle="round,pad=0.08,rounding_size=0.35",
                facecolor=state_color,
                edgecolor="none",
                alpha=STATE_ALPHA,
                zorder=2,
            )
        )
        text_color = FG if row.bulk_dominant_state == "Attenuated" else WHITE
        ax.text(
            state_xc,
            cy,
            row.state_abbrev,
            ha="center",
            va="center",
            fontsize=BODY_PT,
            color=text_color,
            fontweight="bold",
            zorder=3,
        )

        _, _, cons_xc = cols["cons"]
        draw_consistency_glyph(ax, cons_xc, cy, row.bulk_consistency, with_letter=True)

        _, _, sn_xc = cols["sn"]
        sn_color = SN_COLORS[row.snRNA_class]
        sn_w, sn_h = 14.5, min(3.2, row_h * 0.72)
        ax.add_patch(
            FancyBboxPatch(
                (sn_xc - sn_w / 2, cy - sn_h / 2),
                sn_w,
                sn_h,
                boxstyle="round,pad=0.08,rounding_size=0.35",
                facecolor=sn_color,
                edgecolor="none",
                zorder=2,
            )
        )
        sn_alpha_note = 0.92 if row.snRNA_class == "Below_detection" else 1.0
        ax.text(
            sn_xc,
            cy,
            row.snRNA_display,
            ha="center",
            va="center",
            fontsize=BODY_PT,
            color=FG,
            alpha=sn_alpha_note,
            zorder=3,
        )

        _, _, dur_xc = cols["dur"]
        dur_color = DUR_COLORS[row.duration_class]
        dur_w, dur_h = 11.5, min(3.2, row_h * 0.72)
        ax.add_patch(
            FancyBboxPatch(
                (dur_xc - dur_w / 2, cy - dur_h / 2),
                dur_w,
                dur_h,
                boxstyle="round,pad=0.08,rounding_size=0.35",
                facecolor=dur_color,
                edgecolor="none",
                zorder=2,
            )
        )
        dur_text = WHITE if row.duration_class in {"Strong_trend", "Weak_trend"} else FG
        ax.text(dur_xc, cy, row.duration_display, ha="center", va="center", fontsize=BODY_PT, color=dur_text, zorder=3)

        _, _, supp_xc = cols["supp"]
        ax.text(
            supp_xc,
            cy,
            row.support_display,
            ha="center",
            va="center",
            fontsize=BODY_PT,
            color=FG,
            fontweight=support_fontweight(int(row.support_count)),
            zorder=3,
        )

        _, _, path_xc = cols["patha"]
        if row.Path_A:
            tw, th = 14.5, min(3.0, row_h * 0.68)
            ax.add_patch(
                FancyBboxPatch(
                    (path_xc - tw / 2, cy - th / 2),
                    tw,
                    th,
                    boxstyle="round,pad=0.06,rounding_size=0.3",
                    facecolor=WHITE,
                    edgecolor=PATH_A_EDGE,
                    linewidth=0.7,
                    zorder=2,
                )
            )
            ax.text(path_xc, cy, "Path A", ha="center", va="center", fontsize=BODY_PT, color=FG, fontweight="bold", zorder=3)
        else:
            ax.text(path_xc, cy, "–", ha="center", va="center", fontsize=BODY_PT, color="#B0B0B0", zorder=3)

    ax.plot([MARGIN_MM, WIDTH_MM - MARGIN_MM], [body_y0, body_y0], color="#D0D0D0", lw=0.55, zorder=2)

    # Legend — cleaned consistency; full formal names for state / snRNA
    ly = MARGIN_MM + LEGEND_H_MM - 1.5
    ax.text(MARGIN_MM, ly, "Legend", ha="left", va="top", fontsize=NOTE_PT, color=MUTED, fontweight="bold")

    item_y = ly - 5.5  # swatches/items sit well below group headings
    h_sw = 2.6
    h_gap = 3.0

    # --- Row 1: Bulk consistency + Bulk state
    gx = MARGIN_MM + 12.0
    ax.text(gx, ly, "Bulk consistency:", ha="left", va="top", fontsize=NOTE_PT, color=FAINT)
    ix = gx + 1.0
    for cons, full in (("High", "High"), ("Moderate", "Moderate"), ("Low", "Low")):
        draw_consistency_glyph(ax, ix + h_sw / 2, item_y + h_sw / 2, cons, with_letter=False, glyph_x=ix + h_sw / 2)
        ax.text(ix + h_sw + 0.8, item_y + h_sw / 2, full, ha="left", va="center", fontsize=NOTE_PT, color=FAINT)
        ix += h_sw + 0.8 + len(full) * 1.35 + h_gap

    gx2 = ix + 2.0
    ax.text(gx2, ly, "Bulk state:", ha="left", va="top", fontsize=NOTE_PT, color=FAINT)
    ix2 = gx2 + 1.0
    for state, full in STATE_LEGEND.items():
        ax.add_patch(Rectangle((ix2, item_y), h_sw, h_sw, facecolor=STATE_COLORS[state], edgecolor="none"))
        ax.text(ix2 + h_sw + 0.8, item_y + h_sw / 2, full, ha="left", va="center", fontsize=NOTE_PT, color=FAINT)
        ix2 += h_sw + 0.8 + len(full) * 1.22 + h_gap

    # --- Row 2: LA snRNA + LAA duration
    ly2 = ly - 8.5
    item_y2 = ly2 - 6.0
    gx = MARGIN_MM + 12.0
    ax.text(gx, ly2, "LA snRNA:", ha="left", va="top", fontsize=NOTE_PT, color=FAINT)
    ix = gx + 1.0
    for sn, lab in SN_LEGEND.items():
        ax.add_patch(Rectangle((ix, item_y2), h_sw, h_sw, facecolor=SN_COLORS[sn], edgecolor="none"))
        ax.text(ix + h_sw + 0.8, item_y2 + h_sw / 2, lab, ha="left", va="center", fontsize=NOTE_PT, color=FAINT)
        ix += h_sw + 0.8 + len(lab) * 1.22 + h_gap

    gx2 = ix + 2.0
    ax.text(gx2, ly2, "LAA duration:", ha="left", va="top", fontsize=NOTE_PT, color=FAINT)
    ix2 = gx2 + 1.0
    for dur, lab in DUR_LABEL.items():
        ax.add_patch(Rectangle((ix2, item_y2), h_sw, h_sw, facecolor=DUR_COLORS[dur], edgecolor="none"))
        ax.text(ix2 + h_sw + 0.8, item_y2 + h_sw / 2, lab, ha="left", va="center", fontsize=NOTE_PT, color=FAINT)
        ix2 += h_sw + 0.8 + len(lab) * 1.22 + h_gap

    ax.text(
        MARGIN_MM,
        MARGIN_MM + 0.4,
        "Context columns do not enter support count. Path A requires all three gate columns.",
        ha="left",
        va="bottom",
        fontsize=NOTE_PT,
        color=FAINT,
    )

    return fig, ax, rule_segments


def validate_layout(fig: plt.Figure, ax: plt.Axes, rule_segments: list[tuple[float, float, float, float, str]]) -> dict[str, float | int]:
    """Collect all text bboxes and assert no overlap, no clipping, font floor, and rule-text separation."""
    fig.canvas.draw()
    texts: list[tuple[str, tuple[float, float, float, float], float]] = []
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
        fontsize = float(artist.get_fontsize())
        texts.append((label, (x0, y0, x1, y1), fontsize))

    # Allow gate-criteria text to be 5.2 pt; everything else must be >= 6.0 pt.
    min_other = min((fs for label, _, fs in texts if "gate:" not in label), default=99.0)
    require(min_other >= 6.0 - 0.05, f"non-gate text below 6.0 pt: {min_other:.2f}")
    min_gate = min((fs for label, _, fs in texts if "gate:" in label), default=99.0)
    require(min_gate >= GATE_CRIT_PT - 0.05, f"gate text below {GATE_CRIT_PT} pt: {min_gate:.2f}")

    # Header baselines (same as render) for position-based classification.
    table_top = HEIGHT_MM - MARGIN_MM - TITLE_H_MM
    header_h = 19.5 if LAYOUT_MODE == "composite" else 18.5
    group_y = table_top - 2.5
    col_header_y = table_top - (9.6 if LAYOUT_MODE == "composite" else 9.5)
    gate_def_y = table_top - (15.6 if LAYOUT_MODE == "composite" else 14.5)
    header_y0 = table_top - header_h
    title_labels = {
        "B",
        "Family-level cross-layer triangulation matrix",
        "Context annotations and three Path A gates in the frozen Fig 1D order",
        "HCA carrier and bulk dominant state are context only; support count uses bulk, snRNA and duration gates.",
    }

    def at_baseline(y: float, target: float, tol: float = 0.8) -> bool:
        return abs(y - target) <= tol

    def y_range_for(target_y: float, tol: float = 0.8) -> tuple[float, float]:
        bboxes = [bbox for label, bbox, fs in texts if at_baseline(bbox[1] + (bbox[3] - bbox[1]) / 2, target_y, tol) and label not in title_labels]
        require(bboxes, f"no header text found near y={target_y}")
        return min(b[1] for b in bboxes), max(b[3] for b in bboxes)

    group_y_min, group_y_max = y_range_for(group_y)
    col_y_min, col_y_max = y_range_for(col_header_y)
    gate_y_min, gate_y_max = y_range_for(gate_def_y)

    # Strict three-level header assertions (larger y == higher on page).
    require(group_y_min > col_y_max + 0.5, f"group header intersects column header: group bottom {group_y_min:.2f} <= col top {col_y_max:.2f}")
    require(col_y_min > gate_y_max + 0.5, f"column header intersects gate definition: col bottom {col_y_min:.2f} <= gate top {gate_y_max:.2f}")

    # No horizontal or vertical rule crosses visible text.
    for rx0, ry0, rx1, ry1, kind in rule_segments:
        for label, bbox, _ in texts:
            if label in title_labels:
                continue
            # Header rules are not allowed to cross any text in the header region.
            # Vertical group rules run the full table height and must stay in gutters (no text within tolerance).
            if kind == "horizontal":
                # Rule is a horizontal line segment at y = ry0 = ry1.
                if bbox[1] - 0.2 <= ry0 <= bbox[3] + 0.2 and not (bbox[2] < rx0 - 0.2 or bbox[0] > rx1 + 0.2):
                    require(False, f"horizontal rule at y={ry0:.2f} crosses text {label!r}")
            else:
                # Rule is a vertical line segment at x = rx0 = rx1.
                # Only flag if the text is horizontally centered on the rule (i.e. the rule cuts through the text box).
                if bbox[0] - 0.2 <= rx0 <= bbox[2] + 0.2 and not (bbox[3] < ry0 - 0.2 or bbox[1] > ry1 + 0.2):
                    require(False, f"vertical rule at x={rx0:.2f} crosses text {label!r}")

    # Overlap check between any two text boxes (excluding the title block itself, which is intentionally stacked).
    for i, (label_a, bbox_a, _) in enumerate(texts):
        for label_b, bbox_b, _ in texts[i + 1 :]:
            if label_a in title_labels and label_b in title_labels:
                continue
            x0 = max(bbox_a[0], bbox_b[0])
            y0 = max(bbox_a[1], bbox_b[1])
            x1 = min(bbox_a[2], bbox_b[2])
            y1 = min(bbox_a[3], bbox_b[3])
            if x1 > x0 and y1 > y0:
                overlap = (x1 - x0) * (y1 - y0)
                require(overlap < 0.15, f"text overlap: {label_a!r} vs {label_b!r} ({overlap:.2f} mm²)")

    return {
        "n_text_objects": len(texts),
        "min_font_pt": min(fs for _, _, fs in texts),
        "min_non_gate_font_pt": min_other,
        "min_gate_font_pt": min_gate,
    }


def write_summary(source: pd.DataFrame) -> pd.DataFrame:
    tiers = Counter(support_tier(int(c)) for c in source["support_count"])
    path_a = sorted(source.loc[source["Path_A"], "family_id"].astype(int))
    rows = [
        {"metric": "n_families", "value": len(source)},
        {"metric": "order", "value": ",".join(f"F{i}" for i in source["family_id"])},
        {"metric": "Path_A_families", "value": ",".join(f"F{i}" for i in path_a)},
        {"metric": "triple", "value": tiers["triple"]},
        {"metric": "double", "value": tiers["double"]},
        {"metric": "single", "value": tiers["single"]},
        {"metric": "silent", "value": tiers["silent"]},
        {"metric": "silent_families", "value": ",".join(f"F{i}" for i in sorted(EXPECTED_SILENT))},
        {"metric": "double_families", "value": "F9"},
        {"metric": "hca_enters_support", "value": False},
        {"metric": "dominant_state_enters_support", "value": False},
        {"metric": "geomx_present", "value": False},
        {"metric": "genetics_present", "value": False},
        {"metric": "priorityscore_read", "value": False},
        {"metric": "legacy_fig5_read", "value": False},
        {"metric": "dominant_state_counts_shown", "value": False},
        {"metric": "canvas_mm", "value": f"{WIDTH_MM:.0f}x{HEIGHT_MM:.0f}"},
        {"metric": "assertion_status", "value": "PASS"},
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    return summary


def write_qc(manifest: pd.DataFrame, source: pd.DataFrame, hashes: dict[str, str], layout_stats: dict[str, float | int]) -> None:
    path_a = sorted(source.loc[source["Path_A"], "family_id"].astype(int))
    OUT_QC.write_text(
        f"""# Figure 5B cross-layer triangulation matrix QC

Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}

Status: **PASS**

## Scope
- Panel B only. Figure 5A / 5C / 5D and formal routing were not modified.
- Renderer read only `figure5_authoritative_input_manifest.tsv` and `figure5_triangulation_canonical.tsv`.
- No Atlas PriorityScore field opened.
- No root/legacy Figure 5 table read.
- No Fig 2–4 source tables reopened in the renderer.
- Dominant-state concordance counts omitted (not present in allowed inputs).

## Assertions
- Order: F1–F8, F12–F14, F9–F11, F15–F16 — **PASS**
- Path A = F10, F15 at 3/3 — **PASS** ({path_a})
- F9 = 2/3, not Path A — **PASS**
- Silent F4,F5,F6 = 0/3 — **PASS**
- Counts triple/double/single/silent = 2/1/10/3 — **PASS**
- HCA and dominant state do not enter support_count — **PASS**
- GeoMx and genetics absent from panel — **PASS**

## Layout
- Canvas: {WIDTH_MM:.0f} × {HEIGHT_MM:.0f} mm
- Header: strict three-level (group / column / gate definitions), no overlap
- Group headers: CONTEXT / PATH A GATES / TRIAGE with light grey / white / light grey bands
- Gate criteria under headers: High / non-Below-detection / non-Stable
- Continuous LA (11) and RA (5) direction strips with single LA/RA labels
- Path A: outlined "Path A" tags only; no PASS/FAIL or full-row highlight
- Legend groups: Bulk consistency / Bulk state / LA snRNA / LAA duration
- Legend: formal Amplified / Preserved / Attenuated / Reversed and Below-detection
- Visible text objects: {layout_stats["n_text_objects"]}
- Minimum non-gate font: {float(layout_stats["min_non_gate_font_pt"]):.1f} pt
- Minimum gate-criteria font: {float(layout_stats["min_gate_font_pt"]):.1f} pt
- No text overlap, no canvas clipping (validated)

## Inputs
- `{IN_MANIFEST.relative_to(ROOT)}` — `{sha256(IN_MANIFEST)}`
- `{IN_CANONICAL.relative_to(ROOT)}` — `{sha256(IN_CANONICAL)}`

Manifest roles recorded (context only; not reopened here):
{chr(10).join(f"- `{r.artifact_role}`: `{r.artifact_path}`" for r in manifest.itertuples(index=False))}

## Outputs
{chr(10).join(f"- `{path}` — `{digest}`" for path, digest in sorted(hashes.items()))}
""",
        encoding="utf-8",
    )


def main() -> None:
    assert_no_forbidden_reads()
    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    manifest, canonical = load_and_assert()
    source = build_source(canonical)
    source.to_csv(OUT_SOURCE, sep="\t", index=False)
    write_summary(source)

    fig, ax, rule_segments = render(source)
    layout_stats = validate_layout(fig, ax, rule_segments)
    pdf = OUT_BASE.with_suffix(".pdf")
    svg = OUT_BASE.with_suffix(".svg")
    png = OUT_BASE.with_suffix(".png")
    fig.savefig(pdf, bbox_inches=None, pad_inches=0)
    fig.savefig(svg, bbox_inches=None, pad_inches=0)
    fig.savefig(png, dpi=600, bbox_inches=None, pad_inches=0)
    fig.savefig(OUT_PREVIEW, dpi=150, bbox_inches=None, pad_inches=0)
    plt.close(fig)

    svg_text = svg.read_text(encoding="utf-8")
    require("PriorityScore" not in svg_text, "PriorityScore in SVG")
    require("NCP_v3" not in svg_text, "NCP_v3 in SVG")
    require("GeoMx" not in svg_text, "GeoMx in SVG")
    require("PASS" not in svg_text and "FAIL" not in svg_text, "PASS/FAIL in SVG")
    require("Family-level cross-layer triangulation matrix" in svg_text, "title missing")
    require("CONTEXT" in svg_text and "PATH A GATES" in svg_text and "TRIAGE" in svg_text, "group headers missing")
    require("Context columns do not enter support count" in svg_text, "legend note missing")
    require(
        "gate: High" in svg_text
        and "gate: non-" in svg_text
        and "Below-detection" in svg_text
        and "Stable" in svg_text,
        "gate criteria missing",
    )
    require("Below-detection" in svg_text, "formal Below-detection legend missing")
    require("Amplified" in svg_text and "Preserved" in svg_text and "Attenuated" in svg_text and "Reversed" in svg_text, "formal state legend missing")
    require("Path A" in svg_text, "Path A tag text missing")
    for label in LABEL_FIXES.values():
        require(label in svg_text, f"publication label missing: {label}")
    require("PriorityScore" not in "".join(source.columns.astype(str)), "PriorityScore in source")
    require(layout_stats["min_non_gate_font_pt"] >= 6.0, "non-gate font below 6.0 pt")
    from PIL import Image

    img = Image.open(png)
    w_mm = img.size[0] * 25.4 / 600.0
    h_mm = img.size[1] * 25.4 / 600.0
    require(abs(w_mm - WIDTH_MM) < 1.0, f"PNG width {w_mm:.2f} != {WIDTH_MM}")
    require(abs(h_mm - HEIGHT_MM) < 1.0, f"PNG height {h_mm:.2f} != {HEIGHT_MM}")

    hashes = {
        str(p.relative_to(ROOT)): sha256(p)
        for p in (pdf, svg, png, OUT_PREVIEW, OUT_SOURCE, OUT_SUMMARY, IN_MANIFEST, IN_CANONICAL)
    }
    write_qc(manifest, source, hashes, layout_stats)

    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
