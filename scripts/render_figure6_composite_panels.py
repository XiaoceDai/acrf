#!/usr/bin/env python3
"""
Composite-specific Figure 6 A–D vector renders (2 × 2 layout).

Reads frozen source tables only — no coloc or external-set rematch.
Writes exact-size panel SVG/PDF/PNG for nesting at 100% scale.
"""
from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle
from matplotlib.transforms import Bbox

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from figure6_composite_geometry import (  # noqa: E402
    AB_A_WIDTH,
    AB_B_WIDTH,
    AB_BODY_BOT,
    AB_BODY_TOP,
    AB_CAL_SEP_Y,
    AB_FOOTER_BASELINE,
    AB_HEAD_H,
    AB_HEIGHT,
    AB_ROW_CENTERS,
    AB_ROW_H,
    AB_SUB_Y,
    AB_TITLE_Y,
    ACCENT_DIST,
    ACCENT_EGENE,
    ACCENT_GWAS,
    ACCENT_PP3,
    ACCENT_PP4,
    ACCENT_RARE,
    ACCENT_REG,
    ACCENT_RISK,
    ANCHOR_ORDER,
    BAND_CAL,
    CD_BODY_BOT,
    CD_BODY_TOP,
    CD_C_WIDTH,
    CD_D_WIDTH,
    CD_FOOTER_BASELINE,
    CD_HEIGHT,
    CD_LETTER_Y,
    CD_SUB_Y,
    FG,
    FONT_FAMILY,
    FS_BODY,
    FS_FOOT,
    FS_HEADER,
    FS_LETTER,
    FS_SUMMARY,
    FS_TITLE,
    MUTED,
    NA_COLOR,
    OPEN_EDGE,
    OPEN_FACE,
    ROW_LINE,
    SMAD9_HL,
    write_ab_geometry,
    write_cd_geometry,
)

FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
PANEL_DIR = FIGURES / "_composite_fig6"

SRC_A = RESULTS / "fig6a_anchor_coloc_calibration_source.tsv"
SRC_B = RESULTS / "fig6b_posterior_lead_distance_source.tsv"
SRC_C = RESULTS / "fig6c_family_af_overlap_source.tsv"
SRC_D = RESULTS / "fig6d_non_nested_traceability_source.tsv"

SHA_A = "6a688697ead4b86c65551a4c6eb709e89b259f3ad90ca36555df4ac87329261d"
SHA_B = "be0751c6c969335c66db56db52dfa854bb71e44835be37a79a637ad82023dd8f"
SHA_C = "aa52c0d960f284d4c52f77b56a8c69f744dd32d57493f35a50ea7d51684593e7"
SHA_D = "9a8c186d927e2f66996172e80ffb3220d89c05b20766aac3eb9106f5397f2bce"

SET_ORDER = [
    ("yuan2025_cross_population", "Yuan 2025", 29),
    ("roselli2018_loci", "Roselli 2018", 21),
    ("selewa2023_regulatory", "Selewa 2023", 20),
    ("choi_mendelian_rare", "Choi 2025", 22),
]
FAMILY_ORDER = [3, 13, 9]
FAMILY_LABELS = {
    3: "Adipo-lipid axis",
    13: "Amino-acid metabolism",
    9: "HS-GAG / bone / GPCR",
}
EXPECTED_CELLS = {
    (3, "yuan2025_cross_population"): "CAV1",
    (3, "roselli2018_loci"): "CAV1",
    (13, "yuan2025_cross_population"): "CAV1",
    (13, "roselli2018_loci"): "CAV1",
    (9, "selewa2023_regulatory"): "HEY2",
}

STATUS_PRESENT = "present"
STATUS_NO = "evaluated_no_support"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _bool_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def _fixed_bbox(w_mm: float, h_mm: float) -> Bbox:
    return Bbox.from_bounds(0.0, 0.0, w_mm / 25.4, h_mm / 25.4)


def _setup(w_mm: float, h_mm: float):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT_FAMILY, "Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4), facecolor="white", dpi=300)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, w_mm)
    ax.set_ylim(0, h_mm)
    ax.set_axis_off()
    return fig, ax


def _save(fig, stem: str, w_mm: float, h_mm: float) -> None:
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    targets = {
        "A": FIGURES / "fig6a_anchor_coloc_calibration",
        "B": FIGURES / "fig6b_posterior_lead_distance",
        "C": FIGURES / "fig6c_family_af_overlap",
        "D": FIGURES / "fig6d_traceability_matrix",
    }
    base = targets[stem]
    fixed = _fixed_bbox(w_mm, h_mm)
    for ext, dpi in (("pdf", None), ("svg", None), ("png", 600)):
        path = base.with_suffix(f".{ext}")
        kw = dict(facecolor="white", edgecolor="none", bbox_inches=fixed, pad_inches=0)
        if ext == "png":
            fig.savefig(path, dpi=dpi, **kw)
        else:
            fig.savefig(path, **kw)
        mirror = PANEL_DIR / f"panel_{stem.lower()}.{ext}"
        mirror.write_bytes(path.read_bytes())
    plt.close(fig)


def _draw_symbol(ax, x: float, y: float, state: str, accent: str, r: float = 0.95) -> None:
    if state in ("yes", STATUS_PRESENT, True):
        ax.add_patch(Circle(
            (x, y), r, facecolor=accent, edgecolor=accent,
            linewidth=0.55, zorder=4, clip_on=False,
        ))
    elif state in ("no", STATUS_NO, False):
        ax.add_patch(Circle(
            (x, y), r, facecolor=OPEN_FACE, edgecolor=OPEN_EDGE,
            linewidth=0.85, zorder=4, clip_on=False,
        ))
    else:
        ax.text(
            x, y, "—", ha="center", va="center", fontsize=FS_BODY,
            color=NA_COLOR, zorder=4, clip_on=False,
        )


def _title_block(ax, letter_title: str, line2: str, *, x0: float = 1.2) -> None:
    letter = letter_title.split("  ", 1)[0]
    rest = letter_title.split("  ", 1)[1] if "  " in letter_title else ""
    ax.text(
        x0, AB_TITLE_Y, letter,
        fontsize=FS_LETTER, fontweight="bold", color=FG,
        ha="left", va="center", clip_on=False,
    )
    ax.text(
        x0 + 4.0, AB_TITLE_Y, rest,
        fontsize=FS_TITLE, fontweight="bold", color=FG,
        ha="left", va="center", clip_on=False,
    )
    ax.text(
        x0, AB_SUB_Y, line2,
        fontsize=FS_SUMMARY, color=MUTED,
        ha="left", va="center", clip_on=False,
    )


def _title_block_cd(ax, letter_title: str, line2: str, *, x0: float = 1.2) -> None:
    letter = letter_title.split("  ", 1)[0]
    rest = letter_title.split("  ", 1)[1] if "  " in letter_title else letter_title
    ax.text(
        x0, CD_LETTER_Y, letter,
        fontsize=FS_LETTER, fontweight="bold", color=FG,
        ha="left", va="center", clip_on=False,
    )
    ax.text(
        x0 + 4.2, CD_LETTER_Y, rest,
        fontsize=FS_TITLE, fontweight="bold", color=FG,
        ha="left", va="center", clip_on=False,
    )
    ax.text(
        x0, CD_SUB_Y, line2,
        fontsize=FS_SUMMARY, color=MUTED,
        ha="left", va="center", clip_on=False,
    )


def _family_groups(labels: list[str]) -> list[tuple[str, list[int]]]:
    groups: list[tuple[str, list[int]]] = []
    cur: str | None = None
    idxs: list[int] = []
    for i, lab in enumerate(labels):
        if lab != cur:
            if cur is not None:
                groups.append((cur, idxs))
            cur = lab
            idxs = [i]
        else:
            idxs.append(i)
    if cur is not None:
        groups.append((cur, idxs))
    return groups


def render_a(df: pd.DataFrame) -> None:
    w, h = AB_A_WIDTH, AB_HEIGHT
    fig, ax = _setup(w, h)
    _title_block(ax, "A  Anchor calibration", "Yuan 2025 GCST90624411 · GTEx LAA")

    x0 = 1.0
    fam_w = 14.5
    anc_w = 12.5
    ev_w = (w - x0 - fam_w - anc_w - 1.0) / 3.0
    xs = [
        x0,
        x0 + fam_w,
        x0 + fam_w + anc_w,
        x0 + fam_w + anc_w + ev_w,
        x0 + fam_w + anc_w + 2 * ev_w,
        w - 0.6,
    ]
    head_bot = AB_BODY_TOP
    ax.add_patch(Rectangle(
        (x0, head_bot), xs[-1] - x0, AB_HEAD_H,
        facecolor="#F7F7F6", edgecolor="none", zorder=0,
    ))
    headers = [
        (0, "Family / role", FG),
        (1, "Anchor", FG),
        (2, "Local GWAS\n7/12", ACCENT_GWAS),
        (3, "Atrial eGene\n1/12", ACCENT_EGENE),
        (4, "PP4 ≥ 0.70\n0/12", ACCENT_PP4),
    ]
    for i, text, color in headers:
        cx = 0.5 * (xs[i] + xs[i + 1])
        ax.text(
            cx, head_bot + AB_HEAD_H * 0.52, text,
            ha="center", va="center", fontsize=FS_HEADER_SAFE(),
            color=color, fontweight="semibold", linespacing=1.15, clip_on=False,
        )

    df = df.set_index("gene").loc[ANCHOR_ORDER].reset_index()
    groups = _family_groups(df["group_compact"].astype(str).tolist())

    for i, r in df.iterrows():
        y = AB_ROW_CENTERS[i]
        y0 = y - 0.5 * AB_ROW_H
        if r["gene"] == "SMAD9":
            ax.add_patch(Rectangle(
                (x0, y0), xs[-1] - x0, AB_ROW_H,
                facecolor=SMAD9_HL, edgecolor="none", zorder=0,
            ))
        if i >= 10:
            ax.add_patch(Rectangle(
                (x0, y0), xs[-1] - x0, AB_ROW_H,
                facecolor=BAND_CAL, edgecolor="none", zorder=0, alpha=0.85,
            ))
        ax.plot([x0, xs[-1]], [y0, y0], color=ROW_LINE, lw=0.4, zorder=1)
        ax.text(
            xs[1] + 0.4, y, r["gene"],
            ha="left", va="center", fontsize=FS_BODY, color=FG,
            fontstyle="italic", clip_on=False,
        )
        gwas = "yes" if bool(r["local_AF_GWAS_support"]) else "no"
        egene = "yes" if bool(r["atrial_eGene_bonferroni"]) else "no"
        pp4 = "yes" if bool(r["PP4_ge_0.70"]) else "no"
        _draw_symbol(ax, 0.5 * (xs[2] + xs[3]), y, gwas, ACCENT_GWAS, r=0.85)
        _draw_symbol(ax, 0.5 * (xs[3] + xs[4]), y, egene, ACCENT_EGENE, r=0.85)
        _draw_symbol(ax, 0.5 * (xs[4] + xs[5]), y, pp4, ACCENT_PP4, r=0.85)

    for lab, idxs in groups:
        y_top = AB_ROW_CENTERS[idxs[0]] + 0.42 * AB_ROW_H
        y_bot = AB_ROW_CENTERS[idxs[-1]] - 0.42 * AB_ROW_H
        bx = xs[0] + 0.7
        ax.plot([bx, bx], [y_bot, y_top], color="#9A9A9A", lw=0.7, zorder=2)
        ax.plot([bx, bx + 0.6], [y_top, y_top], color="#9A9A9A", lw=0.7, zorder=2)
        ax.plot([bx, bx + 0.6], [y_bot, y_bot], color="#9A9A9A", lw=0.7, zorder=2)
        y_mid = 0.5 * (y_top + y_bot)
        text = "AF\ncontrols" if lab == "AF calibration controls" else lab
        ax.text(
            bx + 1.0, y_mid, text,
            ha="left", va="center", fontsize=FS_BODY - 0.3, color=MUTED,
            linespacing=1.05, clip_on=False,
        )

    ax.plot(
        [x0, xs[-1]], [AB_CAL_SEP_Y, AB_CAL_SEP_Y],
        color="#B8B5B0", lw=0.85, zorder=3,
    )
    note_y = AB_FOOTER_BASELINE + 4.2
    ax.text(
        x0, note_y, "SMAD9: atrial eGene+ · PP4-low",
        ha="left", va="center", fontsize=FS_FOOT, color=MUTED, clip_on=False,
    )
    ax.text(
        x0, note_y - 2.4,
        "PITX2 / KCNN3: calibration controls also unsupported",
        ha="left", va="center", fontsize=FS_FOOT, color=MUTED, clip_on=False,
    )
    _save(fig, "A", w, h)


def FS_HEADER_SAFE() -> float:
    return FS_HEADER - 0.2


def render_b(df: pd.DataFrame) -> None:
    w, h = AB_B_WIDTH, AB_HEIGHT
    fig, ax = _setup(w, h)
    _title_block(
        ax,
        "B  Coloc posterior and lead separation",
        "PP3-dominant 12/12 · PP4 support 0/12 · 0.46–1.92 Mb",
    )

    df = df.set_index("gene").loc[ANCHOR_ORDER].reset_index()
    x0 = 1.2
    mb_w = 8.5
    x_mb = w - 1.0 - mb_w
    dist_w = 28.0
    x_dist0 = x_mb - 1.2 - dist_w
    post_w = x_dist0 - x0 - 2.0
    track_gap = 1.4
    half = (post_w - track_gap) / 2.0
    x_pp3_0, x_pp3_1 = x0, x0 + half
    x_pp4_0, x_pp4_1 = x_pp3_1 + track_gap, x_pp3_1 + track_gap + half

    head_bot = AB_BODY_TOP
    ax.add_patch(Rectangle(
        (x0, head_bot), w - x0 - 0.6, AB_HEAD_H,
        facecolor="#F7F7F6", edgecolor="none", zorder=0,
    ))
    ax.text(
        0.5 * (x_pp3_0 + x_pp3_1), head_bot + AB_HEAD_H * 0.55, "PP3",
        ha="center", va="center", fontsize=FS_HEADER, color=ACCENT_PP3,
        fontweight="semibold", clip_on=False,
    )
    ax.text(
        0.5 * (x_pp4_0 + x_pp4_1), head_bot + AB_HEAD_H * 0.55, "PP4",
        ha="center", va="center", fontsize=FS_HEADER, color=ACCENT_PP4,
        fontweight="semibold", clip_on=False,
    )
    ax.text(
        0.5 * (x_dist0 + x_mb - 1.2), head_bot + AB_HEAD_H * 0.55,
        "Lead separation",
        ha="center", va="center", fontsize=FS_HEADER, color=ACCENT_DIST,
        fontweight="semibold", clip_on=False,
    )
    ax.text(
        x_mb + mb_w * 0.55, head_bot + AB_HEAD_H * 0.55, "Mb",
        ha="center", va="center", fontsize=FS_HEADER, color=ACCENT_DIST,
        fontweight="semibold", clip_on=False,
    )

    band_x0 = x_pp4_0 + 0.70 * (x_pp4_1 - x_pp4_0)
    ax.add_patch(Rectangle(
        (band_x0, AB_BODY_BOT), x_pp4_1 - band_x0, AB_BODY_TOP - AB_BODY_BOT,
        facecolor="#E8DFE9", edgecolor="none", alpha=0.85, zorder=0,
    ))
    ax.text(
        0.5 * (band_x0 + x_pp4_1), AB_BODY_TOP - 0.9, "≥0.70",
        ha="center", va="top", fontsize=FS_FOOT, color=ACCENT_PP4,
        clip_on=False,
    )

    dmax = float(df["delta_mb"].max())
    dmax_plot = max(2.0, dmax * 1.02)

    for i, r in df.iterrows():
        y = AB_ROW_CENTERS[i]
        y0 = y - 0.5 * AB_ROW_H
        if r["gene"] == "SMAD9":
            ax.add_patch(Rectangle(
                (x0, y0), w - x0 - 0.6, AB_ROW_H,
                facecolor=SMAD9_HL, edgecolor="none", zorder=0,
            ))
        if i >= 10:
            ax.add_patch(Rectangle(
                (x0, y0), w - x0 - 0.6, AB_ROW_H,
                facecolor=BAND_CAL, edgecolor="none", zorder=0, alpha=0.85,
            ))
        ax.plot([x0, w - 0.6], [y0, y0], color=ROW_LINE, lw=0.4, zorder=1)

        pp3 = float(r["coloc_PP3"])
        pp4 = float(r["coloc_PP4"])
        ax.plot([x_pp3_0, x_pp3_1], [y, y], color="#DDD8D0", lw=0.7, zorder=1)
        ax.plot([x_pp4_0, x_pp4_1], [y, y], color="#DDD8D0", lw=0.7, zorder=1)
        x3 = x_pp3_0 + np.clip(pp3, 0, 1) * (x_pp3_1 - x_pp3_0)
        x4 = x_pp4_0 + np.clip(pp4, 0, 1) * (x_pp4_1 - x_pp4_0)
        ax.plot([x_pp3_0, x3], [y, y], color=ACCENT_PP3, lw=1.35, zorder=2)
        ax.plot([x_pp4_0, x4], [y, y], color=ACCENT_PP4, lw=1.35, zorder=2)
        ax.plot(x3, y, "o", ms=3.2, color=ACCENT_PP3, zorder=3)
        ax.plot(x4, y, "o", ms=3.2, color=ACCENT_PP4, zorder=3)

        dmb = float(r["delta_mb"])
        xd = x_dist0 + (dmb / dmax_plot) * dist_w
        ax.plot([x_dist0, xd], [y, y], color=ACCENT_DIST, lw=1.1, zorder=2)
        ax.plot(xd, y, "o", ms=3.4, color=ACCENT_DIST, zorder=3)
        ax.text(
            x_mb + mb_w - 0.4, y, f"{dmb:.2f}",
            ha="right", va="center", fontsize=FS_BODY, color=FG,
            clip_on=False, family="monospace",
        )

    ax.plot(
        [x0, w - 0.6], [AB_CAL_SEP_Y, AB_CAL_SEP_Y],
        color="#B8B5B0", lw=0.85, zorder=3,
    )
    ax.text(
        x0, AB_FOOTER_BASELINE,
        "Lead separation is descriptive and is not used in the coloc posterior.",
        ha="left", va="center", fontsize=FS_FOOT, color=MUTED, clip_on=False,
    )
    _save(fig, "B", w, h)


def render_c(df: pd.DataFrame) -> None:
    w, h = CD_C_WIDTH, CD_HEIGHT
    fig, ax = _setup(w, h)
    _title_block_cd(
        ax,
        "C  External AF-gene overlap",
        "5 cells · 2 genes · 3/16 families · 0/64 FDR",
    )

    x0 = 1.2
    fam_w = 22.0
    x_mat = x0 + fam_w
    n_sets = 4
    col_w = (w - x_mat - 1.0) / n_sets
    body_top = CD_BODY_TOP
    body_bot = CD_BODY_BOT + 6.5
    row_h = (body_top - body_bot) / len(FAMILY_ORDER)

    header_y0 = CD_SUB_Y - 1.0
    for j, (sid, lab, n) in enumerate(SET_ORDER):
        cx = x_mat + (j + 0.5) * col_w
        ax.text(
            cx, header_y0 - 1.4, lab.replace(" ", "\n"),
            ha="center", va="top", fontsize=FS_HEADER - 0.35, color=FG,
            fontweight="semibold", linespacing=1.1, clip_on=False,
        )
        ax.text(
            cx, header_y0 - 5.0, f"n = {n}",
            ha="center", va="top", fontsize=FS_FOOT, color=MUTED, clip_on=False,
        )

    lookup = {
        (int(r.family_id), r.external_set_id): r
        for r in df.itertuples(index=False)
    }
    role_edge = {
        "yuan2025_cross_population": ACCENT_RISK,
        "roselli2018_loci": ACCENT_RISK,
        "selewa2023_regulatory": ACCENT_REG,
        "choi_mendelian_rare": ACCENT_RARE,
    }
    role_face = {
        "yuan2025_cross_population": "#EDF1F6",
        "roselli2018_loci": "#EDF1F6",
        "selewa2023_regulatory": "#E8F3F1",
        "choi_mendelian_rare": "#F4EBEF",
    }

    for i, fid in enumerate(FAMILY_ORDER):
        y = body_top - (i + 0.5) * row_h
        ax.text(
            x0, y + 0.9, f"F{fid}",
            ha="left", va="center", fontsize=FS_BODY, fontweight="bold",
            color=FG, clip_on=False,
        )
        ax.text(
            x0, y - 1.1, FAMILY_LABELS[fid],
            ha="left", va="center", fontsize=FS_FOOT, color=MUTED, clip_on=False,
        )
        for j, (sid, lab, n) in enumerate(SET_ORDER):
            cx = x_mat + (j + 0.5) * col_w
            row = lookup.get((fid, sid))
            genes = ""
            if row is not None and int(row.n_overlap) > 0:
                genes = str(row.overlap_genes).split(";")[0]
            if genes:
                pw, ph = 7.6, 3.0
                ax.add_patch(FancyBboxPatch(
                    (cx - pw / 2, y - ph / 2), pw, ph,
                    boxstyle="round,pad=0.15,rounding_size=0.8",
                    facecolor=role_face[sid], edgecolor=role_edge[sid],
                    linewidth=0.7, zorder=2, clip_on=False,
                ))
                ax.text(
                    cx, y, genes,
                    ha="center", va="center", fontsize=FS_BODY - 0.2,
                    color=FG, fontstyle="italic", zorder=3, clip_on=False,
                )

    sum_y = body_bot - 1.6
    ax.text(
        x0, sum_y, "Zero overlap: 13/16 families",
        ha="left", va="center", fontsize=FS_FOOT, color=MUTED, clip_on=False,
    )
    ax.text(
        x0, sum_y - 2.2, "Unique genes: CAV1 · HEY2",
        ha="left", va="center", fontsize=FS_FOOT, color=MUTED, clip_on=False,
    )
    ax.text(
        x0, CD_FOOTER_BASELINE,
        "Labels = member overlap; no enrichment passed BH-FDR.",
        ha="left", va="center", fontsize=FS_FOOT, color=MUTED, clip_on=False,
    )
    _save(fig, "C", w, h)


def render_d(df: pd.DataFrame) -> None:
    w, h = CD_D_WIDTH, CD_HEIGHT
    fig, ax = _setup(w, h)
    _title_block_cd(
        ax,
        "D  Non-nested genetic traceability",
        "Independent flags across 12 anchors",
    )

    tag = "No score · not Path A"
    tag_x = w - 1.4
    ax.add_patch(FancyBboxPatch(
        (tag_x - 28.5, CD_LETTER_Y - 1.55), 28.0, 3.1,
        boxstyle="round,pad=0.12,rounding_size=0.6",
        facecolor="#FAFAFA", edgecolor="#B0B0B0", linewidth=0.7,
        zorder=2, clip_on=False,
    ))
    ax.text(
        tag_x - 14.5, CD_LETTER_Y, tag,
        ha="center", va="center", fontsize=FS_FOOT, color=MUTED,
        zorder=3, clip_on=False,
    )

    df = df.set_index("gene").loc[ANCHOR_ORDER].reset_index()
    x0 = 1.2
    fam_w = 16.0
    anc_w = 13.0
    rem = w - x0 - fam_w - anc_w - 1.0
    core_w = rem * 0.36 / 3
    ext_w = rem * 0.64 / 3
    widths = [fam_w, anc_w, core_w, core_w, core_w, ext_w, ext_w, ext_w]
    xs = [x0]
    for ww in widths:
        xs.append(xs[-1] + ww)

    body_top = CD_BODY_TOP
    body_bot = CD_BODY_BOT
    n = 12
    row_h = (body_top - body_bot) / n
    row_centers = [body_top - (i + 0.5) * row_h for i in range(n)]
    cal_sep = body_top - 10 * row_h

    head_h = 7.2
    ax.add_patch(Rectangle(
        (x0, body_top), xs[-1] - x0, head_h,
        facecolor="#F7F7F6", edgecolor="none", zorder=0,
    ))
    ax.text(
        0.5 * (xs[2] + xs[5]), body_top + head_h * 0.78, "CORE",
        ha="center", va="center", fontsize=FS_FOOT - 0.2, color=MUTED,
        fontweight="semibold", clip_on=False,
    )
    ax.text(
        0.5 * (xs[5] + xs[8]), body_top + head_h * 0.78, "EXTERNAL",
        ha="center", va="center", fontsize=FS_FOOT - 0.2, color=MUTED,
        fontweight="semibold", clip_on=False,
    )
    headers = [
        (0, "Family / role", FG),
        (1, "Anchor", FG),
        (2, "G\n7/12", ACCENT_GWAS),
        (3, "E\n1/12", ACCENT_EGENE),
        (4, "C\n0/12", ACCENT_PP4),
        (5, "Risk /\nprioritized\n2/12", ACCENT_RISK),
        (6, "Regulatory\n2/12", ACCENT_REG),
        (7, "Rare /\nfunctional\n0/12", ACCENT_RARE),
    ]
    for i, text, color in headers:
        cx = 0.5 * (xs[i] + xs[i + 1])
        ax.text(
            cx, body_top + head_h * 0.32, text,
            ha="center", va="center", fontsize=FS_HEADER - 0.35, color=color,
            fontweight="semibold", linespacing=1.05, clip_on=False,
        )

    if "group_compact" in df.columns:
        groups = _family_groups(df["group_compact"].astype(str).tolist())
    else:
        groups = _family_groups([
            "F10", "F10", "Metabolic", "F8", "F11", "F11", "F11",
            "F15", "F15", "F15", "AF calibration controls", "AF calibration controls",
        ])

    status_cols = [
        ("local_gwas_status", ACCENT_GWAS),
        ("atrial_egene_status", ACCENT_EGENE),
        ("coloc_status", ACCENT_PP4),
        ("common_risk_status", ACCENT_RISK),
        ("regulatory_status", ACCENT_REG),
        ("rare_functional_status", ACCENT_RARE),
    ]

    for i, r in df.iterrows():
        y = row_centers[i]
        y0 = y - 0.5 * row_h
        if r["gene"] == "SMAD9":
            ax.add_patch(Rectangle(
                (x0, y0), xs[-1] - x0, row_h,
                facecolor=SMAD9_HL, edgecolor="none", zorder=0,
            ))
        if i >= 10:
            ax.add_patch(Rectangle(
                (x0, y0), xs[-1] - x0, row_h,
                facecolor=BAND_CAL, edgecolor="none", zorder=0, alpha=0.85,
            ))
        ax.plot([x0, xs[-1]], [y0, y0], color=ROW_LINE, lw=0.35, zorder=1)
        ax.text(
            xs[1] + 0.35, y, r["gene"],
            ha="left", va="center", fontsize=FS_BODY - 0.1, color=FG,
            fontstyle="italic", clip_on=False,
        )
        for j, (col, accent) in enumerate(status_cols):
            cx = 0.5 * (xs[2 + j] + xs[3 + j])
            _draw_symbol(ax, cx, y, str(r[col]), accent, r=0.78)

    for lab, idxs in groups:
        y_top = row_centers[idxs[0]] + 0.42 * row_h
        y_bot = row_centers[idxs[-1]] - 0.42 * row_h
        bx = xs[0] + 0.55
        ax.plot([bx, bx], [y_bot, y_top], color="#9A9A9A", lw=0.65, zorder=2)
        ax.plot([bx, bx + 0.5], [y_top, y_top], color="#9A9A9A", lw=0.65, zorder=2)
        ax.plot([bx, bx + 0.5], [y_bot, y_bot], color="#9A9A9A", lw=0.65, zorder=2)
        y_mid = 0.5 * (y_top + y_bot)
        text = "AF\ncontrols" if lab == "AF calibration controls" else lab
        ax.text(
            bx + 0.85, y_mid, text,
            ha="left", va="center", fontsize=FS_BODY - 0.4, color=MUTED,
            linespacing=1.05, clip_on=False,
        )

    ax.plot([x0, xs[-1]], [cal_sep, cal_sep], color="#B8B5B0", lw=0.8, zorder=3)
    ax.text(
        x0, CD_FOOTER_BASELINE + 2.0,
        "G = local GWAS · E = atrial eGene · C = PP4 ≥ 0.70",
        ha="left", va="center", fontsize=FS_FOOT, color=MUTED, clip_on=False,
    )
    ax.text(
        x0, CD_FOOTER_BASELINE,
        "● present · ○ evaluated, not supported · — not evaluated",
        ha="left", va="center", fontsize=FS_FOOT, color=MUTED, clip_on=False,
    )
    _save(fig, "D", w, h)


def load_and_assert() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path, expected in (
        (SRC_A, SHA_A), (SRC_B, SHA_B), (SRC_C, SHA_C), (SRC_D, SHA_D),
    ):
        got = sha256(path)
        require(got == expected, f"frozen SHA mismatch {path.name}: {got}")

    a = pd.read_csv(SRC_A, sep="\t")
    b = pd.read_csv(SRC_B, sep="\t")
    c = pd.read_csv(SRC_C, sep="\t")
    d = pd.read_csv(SRC_D, sep="\t")

    for col in ("local_AF_GWAS_support", "atrial_eGene_bonferroni", "PP4_ge_0.70"):
        a[col] = _bool_series(a[col])
    for col in ("pp3_dominant", "pp4_support"):
        b[col] = _bool_series(b[col])

    require(list(a["gene"]) == ANCHOR_ORDER, "A order")
    require(list(b["gene"]) == ANCHOR_ORDER, "B order")
    require(list(d["gene"]) == ANCHOR_ORDER, "D order")
    require(int(a["local_AF_GWAS_support"].sum()) == 7, "GWAS 7/12")
    require(int(a["atrial_eGene_bonferroni"].sum()) == 1, "eGene 1/12")
    require(int(a["PP4_ge_0.70"].sum()) == 0, "PP4 0/12")
    require(int(b["pp3_dominant"].sum()) == 12, "PP3 12/12")
    require(int(b["pp4_support"].sum()) == 0, "B PP4 0/12")
    require(len(b) == 12 and b["delta_mb"].notna().all(), "12 distances")
    dmin, dmax = float(b["delta_mb"].min()), float(b["delta_mb"].max())
    require(abs(dmin - 0.46) < 0.02, f"dmin {dmin}")
    require(abs(dmax - 1.92) < 0.02, f"dmax {dmax}")
    k = b.loc[b["gene"] == "KCNN3", "delta_mb"].iloc[0]
    require(abs(float(k) - 0.46) < 0.02, "KCNN3")
    bmp = b.loc[b["gene"] == "BMP10", "delta_mb"].iloc[0]
    require(abs(float(bmp) - 1.92) < 0.02, "BMP10")

    require(len(c) == 64, "64 tests")
    require(int((c["n_overlap"] > 0).sum()) == 5, "5 nonzero")
    require(int((c["fdr_bh"] < 0.05).sum()) == 0, "FDR 0")
    fams = sorted(c.loc[c["n_overlap"] > 0, "family_id"].astype(int).unique())
    require(fams == [3, 9, 13], f"families {fams}")
    genes: set[str] = set()
    for r in c.loc[c["n_overlap"] > 0].itertuples(index=False):
        genes.update(str(r.overlap_genes).split(";"))
    require(genes == {"CAV1", "HEY2"}, str(genes))
    for (fid, sid), gene in EXPECTED_CELLS.items():
        row = c[(c["family_id"] == fid) & (c["external_set_id"] == sid)].iloc[0]
        require(int(row["n_overlap"]) == 1, f"{fid}x{sid}")
        require(str(row["overlap_genes"]).startswith(gene), f"{fid}x{sid} gene")

    require(int((d["local_gwas_status"] == STATUS_PRESENT).sum()) == 7, "D G")
    require(int((d["atrial_egene_status"] == STATUS_PRESENT).sum()) == 1, "D E")
    require(int((d["coloc_status"] == STATUS_PRESENT).sum()) == 0, "D C")
    require(int((d["common_risk_status"] == STATUS_PRESENT).sum()) == 2, "D risk")
    require(int((d["regulatory_status"] == STATUS_PRESENT).sum()) == 2, "D reg")
    require(int((d["rare_functional_status"] == STATUS_PRESENT).sum()) == 0, "D rare")
    smad = d.loc[d["gene"] == "SMAD9"].iloc[0]
    require(smad["local_gwas_status"] == STATUS_NO, "SMAD9 G")
    require(smad["atrial_egene_status"] == STATUS_PRESENT, "SMAD9 E")
    require(smad["coloc_status"] == STATUS_NO, "SMAD9 C")
    return a, b, c, d


def panel_text_qc() -> None:
    for name, path, ban_labels in (
        ("A", FIGURES / "fig6a_anchor_coloc_calibration.svg", False),
        ("B", FIGURES / "fig6b_posterior_lead_distance.svg", True),
        ("C", FIGURES / "fig6c_family_af_overlap.svg", False),
        ("D", FIGURES / "fig6d_traceability_matrix.svg", False),
    ):
        text = path.read_text(encoding="utf-8")
        vis = re.sub(r"<[^>]+>", " ", text)
        for bad in (
            "NCP_v3", "Miyazawa", "mr_", "G0", "G1", "G2", "G3", "G4",
            "evidence score", "highest tier", ".tsv", "pipeline_wide",
        ):
            require(bad not in vis, f"{name} forbidden {bad}")
        if ban_labels:
            for gene in ANCHOR_ORDER:
                require(
                    not re.search(rf"\b{gene}\b", vis),
                    f"B must not repeat anchor label {gene}",
                )
        sizes = [float(x) for x in re.findall(r"font-size:\s*([0-9.]+)", text)]
        if sizes:
            require(min(sizes) >= FS_FOOT - 0.25, f"{name} min font {min(sizes)}")


def main() -> None:
    write_ab_geometry()
    write_cd_geometry()
    a, b, c, d = load_and_assert()
    # D needs group labels — merge from A if absent
    if "group_compact" not in d.columns:
        d = d.merge(a[["gene", "group_compact"]], on="gene", how="left")
    render_a(a)
    render_b(b)
    render_c(c)
    render_d(d)
    panel_text_qc()
    print("composite panels: PASS")
    print(f"wrote {FIGURES / 'fig6a_anchor_coloc_calibration.svg'}")
    print(f"wrote {FIGURES / 'fig6b_posterior_lead_distance.svg'}")
    print(f"wrote {FIGURES / 'fig6c_family_af_overlap.svg'}")
    print(f"wrote {FIGURES / 'fig6d_traceability_matrix.svg'}")


if __name__ == "__main__":
    main()
