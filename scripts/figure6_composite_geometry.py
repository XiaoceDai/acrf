#!/usr/bin/env python3
"""
Shared Figure 6 A–D composite geometry — coordinated 2 × 2 layout.

Does not recompute genetic statistics.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

# Composite canvas (178 × 174 mm)
GUTTER_H_MM = 4.0
GUTTER_V_MM = 4.0

PANEL_A = {"x": 0.0, "y": 0.0, "w": 72.0, "h": 92.0}
PANEL_B = {
    "x": PANEL_A["w"] + GUTTER_H_MM,
    "y": 0.0,
    "w": 102.0,
    "h": 92.0,
}
PANEL_C = {
    "x": 0.0,
    "y": PANEL_A["h"] + GUTTER_V_MM,
    "w": 70.0,
    "h": 78.0,
}
PANEL_D = {
    "x": PANEL_C["w"] + GUTTER_H_MM,
    "y": PANEL_C["y"],
    "w": 104.0,
    "h": 78.0,
}

COMPOSITE_W_MM = 178.0
COMPOSITE_H_MM = PANEL_C["y"] + PANEL_C["h"]  # 174.0

# Typography (absolute pt; one system across A–D)
FS_LETTER = 10.5
FS_TITLE = 8.5
FS_SUB = 6.3
FS_SUMMARY = 5.9
FS_HEADER = 5.9
FS_BODY = 5.9
FS_FOOT = 5.15

FONT_FAMILY = "DejaVu Sans"

# Evidence colors (type only — not validation)
ACCENT_GWAS = "#5C6B7A"
ACCENT_EGENE = "#7A6554"
ACCENT_PP3 = "#8A7348"
ACCENT_PP4 = "#6B5A6E"
ACCENT_DIST = "#5C6570"
ACCENT_RISK = "#5C6B7A"
ACCENT_REG = "#5E8A84"
ACCENT_RARE = "#9A7A8A"
OPEN_EDGE = "#6E6E6E"
OPEN_FACE = "#FFFFFF"
NA_COLOR = "#A8A8A8"
BAND_CAL = "#F0F0EE"
SMAD9_HL = "#F7F0E8"
ROW_LINE = "#EAEAEA"
RULE = "#D0CDC8"
MUTED = "#555555"
FG = "#222222"

ANCHOR_ORDER = [
    "PRELP", "OMD", "IDH2", "NDUFA8", "BMP10", "SMAD9", "SMAD6",
    "NTRK2", "NTF3", "FRS2", "PITX2", "KCNN3",
]

# A/B shared horizontal row geometry (mm from panel bottom-left)
# Shared height 92 mm; A and B row centers match exactly.
AB_HEIGHT = 92.0
AB_A_WIDTH = 72.0
AB_B_WIDTH = 102.0

AB_TITLE_Y = AB_HEIGHT - 2.4
AB_SUB_Y = AB_TITLE_Y - 3.15
AB_HEAD_H = 8.2
AB_BODY_GAP = 1.1
AB_FOOTER_BASELINE = 5.6
AB_NOTES_H = 6.8  # space for A compact notes above footer
AB_BODY_BOT = AB_FOOTER_BASELINE + AB_NOTES_H  # A uses notes; B uses same body_bot
AB_BODY_TOP = AB_SUB_Y - AB_BODY_GAP - AB_HEAD_H
AB_N_ROWS = 12
AB_ROW_H = (AB_BODY_TOP - AB_BODY_BOT) / AB_N_ROWS
AB_ROW_CENTERS = [AB_BODY_TOP - (i + 0.5) * AB_ROW_H for i in range(AB_N_ROWS)]
AB_CAL_SEP_Y = AB_BODY_TOP - 10 * AB_ROW_H  # before PITX2 (index 10)

OUT_AB_GEOM = ROOT / "results" / "figure6ab_horizontal_geometry.tsv"
OUT_CD_GEOM = ROOT / "results" / "figure6cd_horizontal_geometry.tsv"
OUT_COMPOSITE_GEOM = ROOT / "results" / "figure6_composite_geometry.tsv"

# Legacy paths kept for older readers (written as aliases of horizontal)
OUT_AB_GEOM_LEGACY = ROOT / "results" / "figure6ab_anchor_row_geometry.tsv"
OUT_CD_GEOM_LEGACY = ROOT / "results" / "figure6cd_composite_geometry.tsv"

# C/D shared outer contract (height 78 mm)
CD_HEIGHT = 78.0
CD_C_WIDTH = 70.0
CD_D_WIDTH = 104.0
CD_TOP = CD_HEIGHT
CD_LETTER_Y = CD_HEIGHT - 2.2
CD_SUB_Y = CD_LETTER_Y - 3.0
CD_BODY_GAP = 1.0
CD_FOOTER_BASELINE = 5.4
CD_BODY_BOT = CD_FOOTER_BASELINE + 3.8
CD_BODY_TOP = CD_SUB_Y - CD_BODY_GAP - 7.5  # header band under subtitle


def write_ab_geometry() -> Path:
    rows = [
        {"metric": "panel_a_width", "value": AB_A_WIDTH, "unit": "mm",
         "notes": "labels + core calibration"},
        {"metric": "panel_b_width", "value": AB_B_WIDTH, "unit": "mm",
         "notes": "posterior + distance; no anchor labels"},
        {"metric": "height", "value": AB_HEIGHT, "unit": "mm", "notes": "shared"},
        {"metric": "title_baseline", "value": AB_TITLE_Y, "unit": "mm", "notes": ""},
        {"metric": "subtitle_baseline", "value": AB_SUB_Y, "unit": "mm", "notes": ""},
        {"metric": "body_top", "value": round(AB_BODY_TOP, 4), "unit": "mm",
         "notes": "top of row-0 cell"},
        {"metric": "body_bottom", "value": round(AB_BODY_BOT, 4), "unit": "mm",
         "notes": ""},
        {"metric": "footer_baseline", "value": AB_FOOTER_BASELINE, "unit": "mm",
         "notes": ""},
        {"metric": "head_h", "value": AB_HEAD_H, "unit": "mm", "notes": ""},
        {"metric": "row_h", "value": round(AB_ROW_H, 6), "unit": "mm", "notes": ""},
        {"metric": "cal_sep_y", "value": round(AB_CAL_SEP_Y, 4), "unit": "mm",
         "notes": "separator before PITX2"},
        {"metric": "fs_letter", "value": FS_LETTER, "unit": "pt", "notes": ""},
        {"metric": "fs_title", "value": FS_TITLE, "unit": "pt", "notes": ""},
        {"metric": "fs_sub", "value": FS_SUB, "unit": "pt", "notes": ""},
        {"metric": "fs_body", "value": FS_BODY, "unit": "pt", "notes": ""},
        {"metric": "fs_foot", "value": FS_FOOT, "unit": "pt", "notes": ""},
    ]
    for i, (gene, y) in enumerate(zip(ANCHOR_ORDER, AB_ROW_CENTERS), start=1):
        rows.append({
            "metric": f"row_center_{i:02d}_{gene}",
            "value": round(y, 6),
            "unit": "mm",
            "notes": f"anchor_order={i}",
        })
    OUT_AB_GEOM.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_AB_GEOM, sep="\t", index=False)
    pd.DataFrame(rows).to_csv(OUT_AB_GEOM_LEGACY, sep="\t", index=False)
    return OUT_AB_GEOM


def write_cd_geometry() -> Path:
    rows = [
        {"metric": "panel_top", "panel_C": CD_TOP, "panel_D": CD_TOP,
         "unit": "mm", "notes": ""},
        {"metric": "panel_bottom", "panel_C": 0.0, "panel_D": 0.0,
         "unit": "mm", "notes": ""},
        {"metric": "height", "panel_C": CD_HEIGHT, "panel_D": CD_HEIGHT,
         "unit": "mm", "notes": "shared; do not stretch"},
        {"metric": "width", "panel_C": CD_C_WIDTH, "panel_D": CD_D_WIDTH,
         "unit": "mm", "notes": "C+gutter4+D = 178"},
        {"metric": "panel_letter_baseline", "panel_C": CD_LETTER_Y,
         "panel_D": CD_LETTER_Y, "unit": "mm", "notes": ""},
        {"metric": "title_baseline", "panel_C": CD_LETTER_Y,
         "panel_D": CD_LETTER_Y, "unit": "mm", "notes": "aligned with letter"},
        {"metric": "summary_subtitle_baseline", "panel_C": CD_SUB_Y,
         "panel_D": CD_SUB_Y, "unit": "mm", "notes": ""},
        {"metric": "body_top", "panel_C": round(CD_BODY_TOP, 4),
         "panel_D": round(CD_BODY_TOP, 4), "unit": "mm", "notes": ""},
        {"metric": "body_bottom", "panel_C": round(CD_BODY_BOT, 4),
         "panel_D": round(CD_BODY_BOT, 4), "unit": "mm", "notes": ""},
        {"metric": "footer_baseline", "panel_C": CD_FOOTER_BASELINE,
         "panel_D": CD_FOOTER_BASELINE, "unit": "mm", "notes": ""},
        {"metric": "fs_letter", "panel_C": FS_LETTER, "panel_D": FS_LETTER,
         "unit": "pt", "notes": ""},
        {"metric": "fs_title", "panel_C": FS_TITLE, "panel_D": FS_TITLE,
         "unit": "pt", "notes": ""},
        {"metric": "fs_sub", "panel_C": FS_SUB, "panel_D": FS_SUB,
         "unit": "pt", "notes": ""},
        {"metric": "fs_foot", "panel_C": FS_FOOT, "panel_D": FS_FOOT,
         "unit": "pt", "notes": ""},
    ]
    OUT_CD_GEOM.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_CD_GEOM, sep="\t", index=False)
    pd.DataFrame(rows).to_csv(OUT_CD_GEOM_LEGACY, sep="\t", index=False)
    return OUT_CD_GEOM


def write_composite_geometry() -> Path:
    rows = []
    for name, p in (("A", PANEL_A), ("B", PANEL_B), ("C", PANEL_C), ("D", PANEL_D)):
        rows.append({
            "panel": name,
            "x_mm": p["x"],
            "y_mm": p["y"],
            "width_mm": p["w"],
            "height_mm": p["h"],
            "scale_percent": 100,
        })
    rows.append({
        "panel": "canvas",
        "x_mm": 0.0,
        "y_mm": 0.0,
        "width_mm": COMPOSITE_W_MM,
        "height_mm": COMPOSITE_H_MM,
        "scale_percent": 100,
    })
    OUT_COMPOSITE_GEOM.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_COMPOSITE_GEOM, sep="\t", index=False)
    return OUT_COMPOSITE_GEOM


# Back-compat aliases used by older imports
def write_ab_geometry_alias() -> Path:
    return write_ab_geometry()
