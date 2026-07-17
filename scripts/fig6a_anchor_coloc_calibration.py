#!/usr/bin/env python3
"""Figure 6A — Two-column anchor-level AF GWAS / atrial cis-eQTL calibration."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (  # noqa: E402
    FG_COLOR,
    FS_FOOTNOTE,
    FS_LEGEND,
    FS_PANEL_TITLE,
    FS_TICK,
    MUTED_COLOR,
    set_nature_style,
)
from _figure_export import pdf_safe_text, save_figure_bundle  # noqa: E402

PIPELINE_WIDE = ROOT / "results" / "figure6_mr_anchors_v3_pipeline_wide.tsv"
FAMILY_LABELS = ROOT / "results" / "ncp_v3_family_labels.tsv"
PROVENANCE_TSV = (
    ROOT / "results" / "provenance" / "gwas_gcst90624411"
    / "gwas_gcst90624411_provenance.tsv"
)

OUT_PDF = ROOT / "figures" / "fig6a_anchor_coloc_calibration.pdf"
OUT_SVG = ROOT / "figures" / "fig6a_anchor_coloc_calibration.svg"
OUT_PNG = ROOT / "figures" / "fig6a_anchor_coloc_calibration.png"
OUT_SOURCE = ROOT / "results" / "fig6a_anchor_coloc_calibration_source.tsv"
OUT_MANIFEST = ROOT / "results" / "fig6a_anchor_coloc_calibration_manifest.tsv"
OUT_QC = ROOT / "results" / "fig6a_anchor_coloc_calibration_qc.md"
OUT_TSV_LEGACY = ROOT / "results" / "fig6a_anchor_coloc_calibration.tsv"

TITLE = "A  Anchor-level calibration against AF GWAS and atrial cis-eQTL"
SOURCE_LINE = (
    "Yuan et al. 2025 cross-population AF GWAS · GCST90624411 | "
    "GTEx LAA cis-eQTL"
)

# Fixed publication geometry (mm) — composite contract with Figure 6B
from figure6_composite_geometry import (  # noqa: E402
    AB_ANCHOR_LEFT,
    AB_BODY_BOT,
    AB_BODY_TOP,
    AB_CAL_SEP_Y,
    AB_FAMILY_LEFT,
    AB_FAMILY_W,
    AB_HEAD_H,
    AB_HEIGHT,
    AB_ROW_CENTERS,
    AB_ROW_H,
    AB_WIDTH,
    FS_BODY,
    FS_FOOT,
    FS_LETTER,
    FS_SOURCE,
    FS_TITLE,
    write_ab_geometry,
)

WIDTH_MM = AB_WIDTH
HEIGHT_MM = AB_HEIGHT
MARGIN_L = AB_FAMILY_LEFT
LEFT_W = 118.0          # compact matrix (116–120)
GUTTER = 4.5
RIGHT_W = 54.0          # calibration readout (53–57)
# trailing right pad = 178 - 1.5 - 118 - 4.5 - 54 = 0.0

EXPECTED_MD5 = "04da4c8c82f9d29da802924a02a47441"
EXPECTED_SHA256 = (
    "d883b43c524ec179232d055435071d4125f2bb1abbeffa19d3bb3d924d2ddeeb"
)
EXPECTED_FILE_USED = "harmonised/GCST90624411.h.tsv.gz"

ANCHOR_ORDER = [
    "PRELP", "OMD",
    "IDH2", "NDUFA8",
    "BMP10", "SMAD9", "SMAD6",
    "NTRK2", "NTF3", "FRS2",
    "PITX2", "KCNN3",
]
CALIBRATORS = {"PITX2", "KCNN3"}
PP4_SUPPORT = 0.70
GATE_P = 5e-8
N_ANCHORS = 12

# Compact group keys → full publication names (legend / QC only)
GROUP_COMPACT = {
    "F10": "F10",
    "Metabolic": "Metabolic",
    "F8": "F8",
    "F11": "F11",
    "F15": "F15",
    "AF calibration controls": "AF calibration controls",
}
GROUP_FULL = {
    "F10": "F10 / GAG disease / TGF-beta",
    "Metabolic": "Metabolic",
    "F8": "F8 / Complex I biogenesis",
    "F11": "F11 / BMP signalling",
    "F15": "F15 / NTRK2 / TrkB",
    "AF calibration controls": "AF calibration controls",
}

# Evidence accents (headers + filled symbols + right-column headings only)
ACCENT_GWAS = "#5C6B7A"      # muted slate blue
ACCENT_EGENE = "#7A6554"     # muted warm brown
ACCENT_COLOC = "#6B5A6E"     # muted plum
OPEN_EDGE = "#B5B5B5"
NA_COLOR = "#A8A8A8"
BAND_CAL = "#F2F0EC"
SMAD9_HL = "#F7F0E8"
ROW_LINE = "#E8E8E8"
RULE = "#D0CDC8"

SYM_YES = "●"
SYM_NO = "○"
SYM_NA = "—"

FORBIDDEN_VISIBLE = [
    "Miyazawa",
    "miyazawa",
    "mr_",
    "figure6_mr_anchors_v3_pipeline_wide",
    "NCP_v3",
    "ncp_v3",
]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _bool(v) -> bool:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return False
    if isinstance(v, str):
        s = v.strip().lower()
        if s in ("", "nan", "na", "n/a"):
            return False
        return s in ("true", "1", "yes")
    if pd.isna(v):
        return False
    return bool(v)


def _family_short() -> dict[int, str]:
    if not FAMILY_LABELS.is_file():
        return {}
    fb = pd.read_csv(FAMILY_LABELS, sep="\t")
    return dict(zip(fb["family_id"].astype(int), fb["short_label"].astype(str)))


def _family_display(raw: str, fam_short: dict[int, str]) -> str:
    raw = str(raw) if pd.notna(raw) else ""
    if raw.startswith("F") and raw[1:].isdigit():
        fid = int(raw[1:])
        lab = fam_short.get(fid, "")
        return f"{raw} / {pdf_safe_text(lab)}" if lab else raw
    return pdf_safe_text(raw)


def _compact_group(anchor_family: str, anchor_type: str) -> str:
    if anchor_type == "calibration":
        return "AF calibration controls"
    s = str(anchor_family)
    if s.startswith("F") and s[1:2].isdigit():
        # F10 / ..., F8 / ..., F11 / ..., F15 / ...
        return s.split()[0].split("/")[0]
    if s.lower().startswith("metabolic"):
        return "Metabolic"
    return s.split("/")[0].strip()


def _verdict(pp4_ok: bool, egene: bool) -> str:
    if pp4_ok:
        return "supported"
    if egene:
        return "eGene, PP4-low"
    return "not supported"


def assert_gwas_provenance() -> dict[str, str]:
    if not PROVENANCE_TSV.is_file():
        raise FileNotFoundError(f"Missing provenance freeze: {PROVENANCE_TSV}")
    prov = dict(
        zip(
            *pd.read_csv(PROVENANCE_TSV, sep="\t").iloc[:, :2].astype(str).values.T
        )
    )
    assert prov.get("accession") == "GCST90624411", prov.get("accession")
    assert prov.get("file_used") == EXPECTED_FILE_USED, prov.get("file_used")
    assert prov.get("local_md5") == EXPECTED_MD5, prov.get("local_md5")
    assert prov.get("local_sha256") == EXPECTED_SHA256, prov.get("local_sha256")
    assert str(prov.get("coloc_rerun_required", "")).upper() == "FALSE"
    return prov


def build_source_table() -> pd.DataFrame:
    wide = pd.read_csv(PIPELINE_WIDE, sep="\t").set_index("anchor")
    fam_short = _family_short()
    missing = [g for g in ANCHOR_ORDER if g not in wide.index]
    if missing:
        raise ValueError(f"Anchors missing from frozen table: {missing}")

    rows = []
    for g in ANCHOR_ORDER:
        r = wide.loc[g]
        gwas_p = float(r["gwas_min_p"]) if pd.notna(r["gwas_min_p"]) else np.nan
        near_gw = bool(gwas_p < GATE_P) if pd.notna(gwas_p) else False
        egene = _bool(r.get("bonf_passed")) or _bool(r.get("is_eGene_atrial"))
        legacy = _bool(r.get("gate_5e8_legacy"))
        pp3 = float(r["coloc_PP3"]) if pd.notna(r["coloc_PP3"]) else np.nan
        pp4 = float(r["coloc_PP4"]) if pd.notna(r["coloc_PP4"]) else np.nan
        pp4_ok = bool(pp4 >= PP4_SUPPORT) if pd.notna(pp4) else False
        pp3_dom = (
            bool(pp3 > pp4 and pp4 < PP4_SUPPORT)
            if pd.notna(pp3) and pd.notna(pp4)
            else False
        )
        is_cal = g in CALIBRATORS
        fam_full = _family_display(r.get("family", ""), fam_short)
        rows.append({
            "gene": g,
            "anchor_family": fam_full,
            "group_compact": _compact_group(fam_full, "calibration" if is_cal else "framework"),
            "anchor_type": "calibration" if is_cal else "framework",
            "local_AF_GWAS_support": near_gw,
            "atrial_eGene_bonferroni": egene,
            "legacy_eQTL_p_lt_5e8": legacy,
            "PP3_dominant": pp3_dom,
            "PP4_ge_0.70": pp4_ok,
            "calibration_verdict": _verdict(pp4_ok, egene),
            "gwas_min_p": gwas_p,
            "coloc_PP3": pp3,
            "coloc_PP4": pp4,
        })
    return pd.DataFrame(rows)


def _summary_counts(df: pd.DataFrame) -> dict[str, int]:
    return {
        "n_gwas": int(df["local_AF_GWAS_support"].sum()),
        "n_egene": int(df["atrial_eGene_bonferroni"].sum()),
        "n_pp3": int(df["PP3_dominant"].sum()),
        "n_pp4": int(df["PP4_ge_0.70"].sum()),
        "n": len(df),
    }


def _assert_table(df: pd.DataFrame, counts: dict[str, int]) -> None:
    assert len(df) == N_ANCHORS
    assert list(df["gene"]) == ANCHOR_ORDER
    assert set(df.loc[df["anchor_type"] == "calibration", "gene"]) == CALIBRATORS
    for col in (
        "local_AF_GWAS_support", "atrial_eGene_bonferroni",
        "legacy_eQTL_p_lt_5e8", "PP3_dominant", "PP4_ge_0.70",
    ):
        assert df[col].notna().all()
    assert counts["n_gwas"] == 7
    assert counts["n_egene"] == 1
    assert counts["n_pp3"] == 12
    assert counts["n_pp4"] == 0
    smad = df.loc[df["gene"] == "SMAD9"].iloc[0]
    assert bool(smad["atrial_eGene_bonferroni"]) is True
    assert bool(smad["PP4_ge_0.70"]) is False
    assert smad["calibration_verdict"] == "eGene, PP4-low"


def _family_groups(df: pd.DataFrame) -> list[tuple[str, list[int]]]:
    groups: list[tuple[str, list[int]]] = []
    cur_lab: str | None = None
    cur_idx: list[int] = []
    for i, r in df.iterrows():
        lab = str(r["group_compact"])
        if lab != cur_lab:
            if cur_lab is not None:
                groups.append((cur_lab, cur_idx))
            cur_lab = lab
            cur_idx = [int(i)]
        else:
            cur_idx.append(int(i))
    if cur_lab is not None:
        groups.append((cur_lab, cur_idx))
    return groups


def _draw_symbol(ax, x: float, y: float, state: str, accent: str) -> None:
    if state == "yes":
        ax.text(
            x, y, SYM_YES, ha="center", va="center",
            fontsize=8.2, color=accent, zorder=4, clip_on=False,
        )
    elif state == "no":
        ax.text(
            x, y, SYM_NO, ha="center", va="center",
            fontsize=8.2, color=OPEN_EDGE, zorder=4, clip_on=False,
        )
    else:
        ax.text(
            x, y, SYM_NA, ha="center", va="center",
            fontsize=7.5, color=NA_COLOR, zorder=4, clip_on=False,
        )


def _wrap_right(ax, x: float, y: float, text: str, width_mm: float,
                fontsize: float, color: str, *, weight: str = "normal",
                style: str = "normal", line_mm: float = 3.0) -> float:
    """Simple word-wrap; returns y after last line."""
    words = text.split()
    lines: list[str] = []
    cur = ""
    # ~0.48 mm per char at 6.5 pt (approx)
    max_chars = max(18, int(width_mm / (fontsize * 0.28)))
    for w in words:
        trial = (cur + " " + w).strip()
        if len(trial) > max_chars and cur:
            lines.append(cur)
            cur = w
        else:
            cur = trial
    if cur:
        lines.append(cur)
    for line in lines:
        ax.text(
            x, y, line, ha="left", va="top",
            fontsize=fontsize, color=color, fontweight=weight,
            style=style, clip_on=False, zorder=4,
        )
        y -= line_mm
    return y


def render(df: pd.DataFrame, counts: dict[str, int]) -> plt.Figure:
    set_nature_style()
    fig = plt.figure(
        figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4),
        facecolor="white", dpi=300,
    )
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, WIDTH_MM)
    ax.set_ylim(0, HEIGHT_MM)
    ax.set_axis_off()

    x_left = MARGIN_L
    x_mat_r = x_left + LEFT_W
    x_right = x_mat_r + GUTTER
    x_right_r = x_right + RIGHT_W

    from figure6_composite_geometry import AB_LINE2_Y, AB_TITLE_Y
    ax.text(
        x_left, AB_TITLE_Y, TITLE, fontsize=FS_TITLE, fontweight="bold",
        ha="left", va="top", color=FG_COLOR, clip_on=False,
    )
    ax.text(
        x_left, AB_LINE2_Y, SOURCE_LINE, fontsize=FS_SOURCE,
        ha="left", va="top", color=MUTED_COLOR, clip_on=False,
    )

    n_row = len(df)
    n_framework = int((df["anchor_type"] == "framework").sum())
    head_h = AB_HEAD_H
    head_bot = AB_BODY_TOP  # top of first data row / bottom of header
    matrix_top = head_bot + head_h
    matrix_bot = AB_BODY_BOT
    row_h = AB_ROW_H
    row_centers = list(AB_ROW_CENTERS)

    # Column layout: locked family/anchor left edges; evidence cols fill LEFT_W
    fam_w = AB_FAMILY_W
    gene_w = 22.1
    rest = LEFT_W - fam_w - gene_w
    ev = rest / 3.0
    xs = [
        AB_FAMILY_LEFT,
        AB_FAMILY_LEFT + fam_w,
        AB_ANCHOR_LEFT + gene_w,
        AB_ANCHOR_LEFT + gene_w + ev,
        AB_ANCHOR_LEFT + gene_w + 2 * ev,
        AB_FAMILY_LEFT + LEFT_W,
    ]
    col = {
        "family": (xs[0], xs[1]),
        "gene": (xs[1], xs[2]),
        "gwas": (xs[2], xs[3]),
        "egene": (xs[3], xs[4]),
        "pp4": (xs[4], xs[5]),
    }

    def _cx(key: str) -> float:
        a, b = col[key]
        return 0.5 * (a + b)

    # Light header wash (not a card)
    ax.add_patch(Rectangle(
        (x_left, head_bot), LEFT_W, head_h,
        facecolor="#F6F5F3", edgecolor="none", zorder=1,
    ))
    ax.plot(
        [x_left, x_mat_r], [head_bot, head_bot],
        color=RULE, lw=0.6, zorder=2,
    )

    # Static labels
    ax.text(
        _cx("family"), head_bot + head_h * 0.55, "Family / role",
        ha="center", va="center", fontsize=FS_TICK - 0.4,
        fontweight="bold", color=FG_COLOR, zorder=3,
    )
    ax.text(
        _cx("gene"), head_bot + head_h * 0.55, "Anchor",
        ha="center", va="center", fontsize=FS_TICK - 0.4,
        fontweight="bold", color=FG_COLOR, zorder=3,
    )
    # Count-integrated evidence headers
    header_ev = [
        ("gwas", "Local AF GWAS", f"{counts['n_gwas']}/{counts['n']}", ACCENT_GWAS),
        ("egene", "Atrial eGene", f"{counts['n_egene']}/{counts['n']}", ACCENT_EGENE),
        ("pp4", "Colocalization\nPP4 ≥ 0.70", f"{counts['n_pp4']}/{counts['n']}", ACCENT_COLOC),
    ]
    for key, lab, cnt, accent in header_ev:
        ax.text(
            _cx(key), head_bot + head_h * 0.72, lab,
            ha="center", va="center", fontsize=FS_TICK - 0.6,
            fontweight="bold", color=accent, linespacing=1.05, zorder=3,
        )
        ax.text(
            _cx(key), head_bot + head_h * 0.28, cnt,
            ha="center", va="center", fontsize=FS_TICK,
            fontweight="bold", color=FG_COLOR, zorder=3,
        )

    # Rows — shared A/B centers
    row_y: dict[int, float] = {}
    for i in range(n_row):
        row_y[i] = row_centers[i]
        y_cursor = row_y[i] - row_h / 2
        r = df.iloc[i]

        # Calibration-control block shading
        if r["anchor_type"] == "calibration":
            ax.add_patch(Rectangle(
                (x_left, y_cursor), LEFT_W, row_h,
                facecolor=BAND_CAL, edgecolor="none", zorder=0,
            ))
        elif r["gene"] == "SMAD9":
            ax.add_patch(Rectangle(
                (x_left, y_cursor), LEFT_W, row_h,
                facecolor=SMAD9_HL, edgecolor="none", zorder=0,
            ))

        # Light row separator
        ax.plot(
            [x_left + 0.4, x_mat_r - 0.4], [y_cursor, y_cursor],
            color=ROW_LINE, lw=0.35, zorder=1,
        )

        ax.text(
            _cx("gene"), row_y[i], r["gene"],
            ha="center", va="center",
            fontsize=FS_TICK - 0.1, fontweight="bold",
            color="#5C4030" if r["gene"] == "SMAD9" else FG_COLOR,
            zorder=4,
        )
        _draw_symbol(
            ax, _cx("gwas"), row_y[i],
            "yes" if r["local_AF_GWAS_support"] else "no", ACCENT_GWAS,
        )
        _draw_symbol(
            ax, _cx("egene"), row_y[i],
            "yes" if r["atrial_eGene_bonferroni"] else "no", ACCENT_EGENE,
        )
        _draw_symbol(
            ax, _cx("pp4"), row_y[i],
            "yes" if r["PP4_ge_0.70"] else "no", ACCENT_COLOC,
        )

    # Compact group labels + brackets (once per group)
    for lab, idxs in _family_groups(df):
        y0 = row_y[idxs[-1]] - row_h / 2 + 0.2
        y1 = row_y[idxs[0]] + row_h / 2 - 0.2
        bx = col["family"][0] + 1.0
        ax.plot([bx, bx], [y0, y1], color="#C9C4BC", lw=1.1,
                solid_capstyle="round", zorder=2)
        # end ticks
        ax.plot([bx, bx + 1.0], [y0, y0], color="#C9C4BC", lw=1.0, zorder=2)
        ax.plot([bx, bx + 1.0], [y1, y1], color="#C9C4BC", lw=1.0, zorder=2)
        # label once at group center; calibrators get two-line compact label
        y_mid = 0.5 * (row_y[idxs[0]] + row_y[idxs[-1]])
        if lab == "AF calibration controls":
            ax.text(
                _cx("family") + 0.8, y_mid, "AF\ncalibration\ncontrols",
                ha="center", va="center", fontsize=FS_FOOTNOTE - 0.6,
                color=MUTED_COLOR, linespacing=1.05, zorder=4,
            )
        else:
            ax.text(
                _cx("family") + 0.8, y_mid, lab,
                ha="center", va="center", fontsize=FS_TICK - 0.3,
                fontweight="bold", color=MUTED_COLOR, zorder=4,
            )

    # Thin rule above calibration block (shared A/B y)
    if n_framework < n_row:
        y_div = AB_CAL_SEP_Y
        ax.plot(
            [x_left + 0.3, x_mat_r - 0.3], [y_div, y_div],
            color="#8A8580", lw=0.7, zorder=3,
        )

    ry = matrix_top
    ax.text(
        x_right, ry, "CALIBRATION READOUT",
        ha="left", va="top", fontsize=FS_TICK - 0.2,
        fontweight="bold", color=FG_COLOR, clip_on=False,
    )
    ry -= 3.6
    ax.plot([x_right, x_right_r], [ry, ry], color=RULE, lw=0.7, zorder=2)
    ry -= 1.8

    sections = [
        (
            ACCENT_GWAS,
            "Locus relevance ≠ gene mediation",
            f"{counts['n_gwas']}/{counts['n']} local GWAS windows; "
            f"{counts['n_pp4']}/{counts['n']} coloc support.",
        ),
        (
            ACCENT_EGENE,
            "eGene ≠ colocalization",
            "SMAD9 is atrial-eGene positive but PP4-low.",
        ),
        (
            ACCENT_COLOC,
            "Calibration-limited non-support",
            "PITX2 and KCNN3 also did not meet the same conservative gate.",
        ),
        (
            MUTED_COLOR,
            f"PP3-dominant: {counts['n_pp3']}/{counts['n']}",
            "Continuous posteriors and lead-variant separation are shown in B.",
        ),
    ]
    body_fs = FS_FOOTNOTE - 0.4
    head_fs = FS_FOOTNOTE - 0.1
    for i, (accent, head, body) in enumerate(sections):
        ax.text(
            x_right, ry, head,
            ha="left", va="top", fontsize=head_fs, fontweight="bold",
            color=accent, clip_on=False, zorder=4,
        )
        ry -= 3.1
        ry = _wrap_right(
            ax, x_right, ry, body, RIGHT_W - 1.0,
            fontsize=body_fs, color=MUTED_COLOR, line_mm=2.7,
        )
        ry -= 1.2
        if i < len(sections) - 1:
            ax.plot(
                [x_right, x_right_r], [ry + 0.4, ry + 0.4],
                color=RULE, lw=0.45, zorder=2,
            )
            ry -= 1.4

    ax.text(
        x_left, 4.4,
        f"{SYM_YES} criterion met     {SYM_NO} evaluated, gate not met",
        ha="left", va="center", fontsize=FS_FOOTNOTE - 0.5,
        color=MUTED_COLOR, clip_on=False,
    )
    ax.text(
        x_left, 2.2,
        "Open symbols do not exclude biological involvement.  "
        "Evidence dimensions are non-nested.",
        ha="left", va="center", fontsize=FS_FOOTNOTE - 0.7,
        color=MUTED_COLOR, style="italic", clip_on=False,
    )
    fig._fig6a_layout = {  # type: ignore[attr-defined]
        "row_y": [row_y[i] for i in range(n_row)],
        "row_h": row_h,
        "family_left": AB_FAMILY_LEFT,
        "anchor_left": AB_ANCHOR_LEFT,
        "body_top": AB_BODY_TOP,
        "body_bottom": AB_BODY_BOT,
        "cal_sep_y": AB_CAL_SEP_Y,
    }
    return fig


def _pdf_page_mm(pdf: Path) -> tuple[float, float]:
    info = subprocess.run(
        ["pdfinfo", str(pdf)], check=True, capture_output=True, text=True,
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    if not match:
        raise RuntimeError("Unable to read PDF page size")
    return float(match.group(1)) * 25.4 / 72, float(match.group(2)) * 25.4 / 72


def _svg_mm(svg: Path) -> tuple[float, float]:
    text = svg.read_text(encoding="utf-8")
    m = re.search(
        r'width="([0-9.]+)pt".*?height="([0-9.]+)pt"', text, flags=re.DOTALL,
    )
    if m:
        return float(m.group(1)) * 25.4 / 72, float(m.group(2)) * 25.4 / 72
    m = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', text)
    if not m:
        raise RuntimeError("Unable to read SVG geometry")
    return float(m.group(1)) * 25.4 / 72, float(m.group(2)) * 25.4 / 72


def _png_mm(png: Path, dpi: int = 600) -> tuple[float, float, int, int]:
    with Image.open(png) as im:
        w_px, h_px = im.size
    return w_px * 25.4 / dpi, h_px * 25.4 / dpi, w_px, h_px


def _extract_visible_text(pdf: Path, svg: Path) -> str:
    chunks = [svg.read_text(encoding="utf-8")]
    if shutil.which("pdftotext"):
        out = subprocess.run(
            ["pdftotext", "-layout", str(pdf), "-"],
            check=True, capture_output=True, text=True,
        ).stdout
        chunks.append(out)
    return "\n".join(chunks)


def run_qc(
    df: pd.DataFrame,
    counts: dict[str, int],
    prov: dict[str, str],
) -> dict[str, object]:
    for path in (OUT_PDF, OUT_SVG, OUT_PNG, OUT_SOURCE):
        assert path.is_file(), path

    src = pd.read_csv(OUT_SOURCE, sep="\t")
    assert len(src) == N_ANCHORS
    assert list(src["gene"]) == ANCHOR_ORDER
    assert int(src["local_AF_GWAS_support"].astype(str).str.lower().isin(
        ["true", "1"]).sum()) == 7
    assert int(src["atrial_eGene_bonferroni"].astype(str).str.lower().isin(
        ["true", "1"]).sum()) == 1
    assert int(src["PP3_dominant"].astype(str).str.lower().isin(
        ["true", "1"]).sum()) == 12
    assert int(src["PP4_ge_0.70"].astype(str).str.lower().isin(
        ["true", "1"]).sum()) == 0
    assert "legacy_eQTL_p_lt_5e8" in src.columns
    for col in ("coloc_PP3", "coloc_PP4", "gwas_min_p"):
        assert np.allclose(
            src[col].astype(float).values,
            df[col].astype(float).values,
            equal_nan=True,
        )

    visible = _extract_visible_text(OUT_PDF, OUT_SVG)
    for bad in FORBIDDEN_VISIBLE:
        assert bad not in visible, f"forbidden visible text: {bad}"
    assert not re.search(r"(?<![A-Za-z])MR(?![A-Za-z])", visible)
    assert not re.search(r"\bmr_", visible)
    assert "Yuan et al. 2025" in visible
    assert "GCST90624411" in visible
    assert "GTEx LAA" in visible
    assert "CALIBRATION READOUT" in visible
    assert "7/12" in visible and "1/12" in visible
    assert "0/12" in visible
    assert "PP3-dominant: 12/12" in visible or "PP3-dominant: 12/12" in OUT_SVG.read_text()
    assert "SMAD9" in visible and "PP4-low" in visible
    assert "PITX2" in visible and "KCNN3" in visible
    assert "AF calibration controls" in visible or "calibration" in visible.lower()
    assert visible.count("not supported") == 0
    assert "Legacy eQTL" not in visible
    assert "Calibration verdict" not in visible
    # Long family names should not spam the panel (full names are legend-only)
    assert "GAG disease" not in visible
    assert "Complex I biogenesis" not in visible
    assert "BMP signalling" not in visible
    assert ".tsv" not in visible
    for g in ANCHOR_ORDER:
        assert g in visible, f"missing gene {g}"
    for lab in ("F10", "F8", "F11", "F15", "Metabolic"):
        assert lab in visible, f"missing compact group {lab}"

    pdf_w, pdf_h = _pdf_page_mm(OUT_PDF)
    svg_w, svg_h = _svg_mm(OUT_SVG)
    png_w, png_h, png_px_w, png_px_h = _png_mm(OUT_PNG, dpi=600)
    assert abs(svg_w - png_w) < 0.50, (svg_w, png_w)
    assert abs(svg_h - png_h) < 0.50, (svg_h, png_h)
    assert abs(pdf_w - svg_w) < 2.0, (pdf_w, svg_w)
    assert abs(pdf_h - svg_h) < 2.0, (pdf_h, svg_h)
    assert abs(pdf_w - WIDTH_MM) < 2.0, pdf_w
    assert 68.0 <= pdf_h <= 78.0, pdf_h

    svg_text = OUT_SVG.read_text(encoding="utf-8")
    sizes = [float(x) for x in re.findall(r"font-size:\s*([0-9.]+)px", svg_text)]
    sizes += [float(x) for x in re.findall(r"font-size:\s*([0-9.]+)", svg_text)]
    if sizes:
        assert min(sizes) >= 5.5, min(sizes)

    return {
        "counts": counts,
        "pdf_mm": (round(pdf_w, 2), round(pdf_h, 2)),
        "svg_mm": (round(svg_w, 2), round(svg_h, 2)),
        "png_mm": (round(png_w, 2), round(png_h, 2)),
        "png_px": (png_px_w, png_px_h),
        "provenance_accession": prov["accession"],
        "provenance_md5": prov["local_md5"],
        "coloc_rerun_required": prov["coloc_rerun_required"],
        "min_font_px": min(sizes) if sizes else None,
        "group_full_legend": GROUP_FULL,
    }


def write_manifest(qc: dict[str, object]) -> None:
    rows = [
        {
            "artifact": "figure_pdf",
            "path": str(OUT_PDF.relative_to(ROOT)),
            "sha256": _sha256(OUT_PDF),
            "geometry_mm": f"{qc['pdf_mm'][0]}x{qc['pdf_mm'][1]}",
            "notes": "two-column publication panel",
        },
        {
            "artifact": "figure_svg",
            "path": str(OUT_SVG.relative_to(ROOT)),
            "sha256": _sha256(OUT_SVG),
            "geometry_mm": f"{qc['svg_mm'][0]}x{qc['svg_mm'][1]}",
            "notes": "two-column publication panel",
        },
        {
            "artifact": "figure_png_600dpi",
            "path": str(OUT_PNG.relative_to(ROOT)),
            "sha256": _sha256(OUT_PNG),
            "geometry_mm": f"{qc['png_mm'][0]}x{qc['png_mm'][1]}",
            "notes": f"600 dpi; {qc['png_px'][0]}x{qc['png_px'][1]} px",
        },
        {
            "artifact": "source_table",
            "path": str(OUT_SOURCE.relative_to(ROOT)),
            "sha256": _sha256(OUT_SOURCE),
            "geometry_mm": "",
            "notes": "full frozen fields retained; matrix shows compact subset",
        },
        {
            "artifact": "frozen_anchor_wide_input",
            "path": str(PIPELINE_WIDE.relative_to(ROOT)),
            "sha256": _sha256(PIPELINE_WIDE),
            "geometry_mm": "",
            "notes": "read-only; not shown on panel",
        },
        {
            "artifact": "gwas_provenance",
            "path": str(PROVENANCE_TSV.relative_to(ROOT)),
            "sha256": _sha256(PROVENANCE_TSV),
            "geometry_mm": "",
            "notes": "asserted off-panel",
        },
    ]
    pd.DataFrame(rows).to_csv(OUT_MANIFEST, sep="\t", index=False)


def write_qc(qc: dict[str, object]) -> None:
    c = qc["counts"]
    legend = qc["group_full_legend"]
    hashes = {
        str(p.relative_to(ROOT)): _sha256(p)
        for p in (OUT_PDF, OUT_SVG, OUT_PNG, OUT_SOURCE, OUT_MANIFEST)
    }
    OUT_QC.write_text(
        "\n".join(
            [
                "# Figure 6A anchor coloc calibration QC",
                "",
                f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
                "",
                "Status: **PASS**",
                "",
                "## Scope",
                "- Two-column publication redesign; no numerical changes.",
                "- No coloc recomputation; no frozen-table rebuild.",
                "- No Figure 6 A–D composite routing change.",
                "",
                "## Layout",
                f"- Canvas: {WIDTH_MM:.0f} × {HEIGHT_MM:.0f} mm.",
                f"- Left matrix ≈ {LEFT_W:.1f} mm; gutter {GUTTER:.1f} mm; "
                f"right readout ≈ {RIGHT_W:.1f} mm.",
                "- No KPI cards; counts integrated into matrix column headers.",
                "",
                "## Counts",
                f"- Local AF GWAS = {c['n_gwas']}/{c['n']}",
                f"- Atrial eGene = {c['n_egene']}/{c['n']}",
                f"- PP3-dominant = {c['n_pp3']}/{c['n']}",
                f"- PP4 ≥ 0.70 = {c['n_pp4']}/{c['n']}",
                "",
                "## Compact family / role labels (panel)",
                "- F10, Metabolic, F8, F11, F15, AF calibration controls",
                "",
                "## Full publication-facing family names (legend)",
                *[f"- `{k}`: {v}" for k, v in legend.items()],
                "",
                "## Scientific assertions",
                "- 12 anchors in locked order.",
                "- SMAD9 eGene-positive / PP4-low (subtle row tint; no matrix prose).",
                "- PITX2 / KCNN3 in shaded calibration-control block.",
                "- No repeated per-row verdict text.",
                "",
                "## Forbidden visible text",
                "- zero Miyazawa / MR / mr_ / internal filenames / NCP_v3",
                "- zero long repeated family phrases on panel",
                "",
                "## GWAS provenance (off-panel)",
                f"- accession = `{qc['provenance_accession']}`",
                f"- file used = `{EXPECTED_FILE_USED}`",
                f"- local MD5 = `{qc['provenance_md5']}`",
                f"- coloc_rerun_required = `{qc['coloc_rerun_required']}`",
                "",
                "## Geometry / typography",
                f"- PDF: {qc['pdf_mm'][0]} × {qc['pdf_mm'][1]} mm",
                f"- SVG: {qc['svg_mm'][0]} × {qc['svg_mm'][1]} mm",
                f"- PNG: {qc['png_mm'][0]} × {qc['png_mm'][1]} mm "
                f"({qc['png_px'][0]} × {qc['png_px'][1]} px @ 600 dpi)",
                f"- Minimum detected SVG font-size: {qc['min_font_px']} px",
                "",
                "## SHA256",
                *[f"- `{path}`: `{digest}`" for path, digest in sorted(hashes.items())],
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    prov = assert_gwas_provenance()
    df = build_source_table()
    counts = _summary_counts(df)
    _assert_table(df, counts)

    OUT_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_SOURCE, sep="\t", index=False)
    df.to_csv(OUT_TSV_LEGACY, sep="\t", index=False)
    print(f"wrote {OUT_SOURCE.relative_to(ROOT)} ({len(df)} anchors)")
    print(
        "counts: "
        f"GWAS {counts['n_gwas']}/12 · eGene {counts['n_egene']}/12 · "
        f"PP3 {counts['n_pp3']}/12 · PP4 {counts['n_pp4']}/12"
    )

    write_ab_geometry()
    print(f"wrote results/figure6ab_anchor_row_geometry.tsv")

    from matplotlib.transforms import Bbox
    fixed = Bbox.from_bounds(0.0, 0.0, WIDTH_MM / 25.4, HEIGHT_MM / 25.4)
    fig = render(df, counts)
    pdf_path, png_path = save_figure_bundle(
        fig, OUT_PDF, dpi=600, bbox_inches=fixed, pad_inches=0,
    )
    print(f"wrote {pdf_path.relative_to(ROOT)}")
    if png_path:
        print(f"wrote {png_path.relative_to(ROOT)}")

    fig = render(df, counts)
    OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        OUT_SVG, format="svg",
        bbox_inches=fixed, pad_inches=0,
        facecolor="white", edgecolor="none",
    )
    plt.close(fig)
    print(f"wrote {OUT_SVG.relative_to(ROOT)}")

    qc = run_qc(df, counts, prov)
    write_manifest(qc)
    write_qc(qc)
    print(f"wrote {OUT_MANIFEST.relative_to(ROOT)}")
    print(f"wrote {OUT_QC.relative_to(ROOT)}")
    print("QC: PASS")


if __name__ == "__main__":
    main()
