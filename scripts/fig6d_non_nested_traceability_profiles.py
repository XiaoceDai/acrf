#!/usr/bin/env python3
"""
Figure 6D — Non-nested genetic traceability profiles (publication panel).

Synthesis panel of independent evidence flags for 10 framework anchors and
2 AF calibration controls. Not a new statistical analysis; not a Path A gate.

Authoritative inputs (frozen; read-only):
  results/fig6a_anchor_coloc_calibration_source.tsv
  results/fig6b_posterior_lead_distance_source.tsv
  results/provenance/figure6_external_gene_sets_manifest.tsv
  results/external_af_gene_sets.tsv
  results/ncp_v3_gene_universe.tsv

Optional p-value enrichment (source table only; never display logic):
  results/figure6_mr_anchors_v3_pipeline_wide.tsv  (lead_eqtl_p column)

Outputs:
  results/fig6d_non_nested_traceability_source.tsv
  results/fig6d_non_nested_traceability_summary.tsv
  results/fig6d_non_nested_traceability_manifest.tsv
  results/fig6d_non_nested_traceability_qc.md
  figures/fig6d_traceability_matrix.{pdf,svg,png}
  figures/_pdf_preview/fig6d_traceability_matrix.png
"""
from __future__ import annotations

import hashlib
import re
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Rectangle
from matplotlib.transforms import Bbox

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (  # noqa: E402
    FG_COLOR,
    FONT_FLOOR_HARD,
    MUTED_COLOR,
    set_nature_style,
)
from _figure_export import save_figure_bundle  # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
FIG6A = ROOT / "results" / "fig6a_anchor_coloc_calibration_source.tsv"
FIG6B = ROOT / "results" / "fig6b_posterior_lead_distance_source.tsv"
MANIFEST_EXT = ROOT / "results" / "provenance" / "figure6_external_gene_sets_manifest.tsv"
EXT = ROOT / "results" / "external_af_gene_sets.tsv"
UNI = ROOT / "results" / "ncp_v3_gene_universe.tsv"
# Optional: lead_eqtl_p only (never used for status / display)
PIPELINE_WIDE = ROOT / "results" / "figure6_mr_anchors_v3_pipeline_wide.tsv"

OUT_SOURCE = ROOT / "results" / "fig6d_non_nested_traceability_source.tsv"
OUT_SUMMARY = ROOT / "results" / "fig6d_non_nested_traceability_summary.tsv"
OUT_MANIFEST = ROOT / "results" / "fig6d_non_nested_traceability_manifest.tsv"
OUT_QC = ROOT / "results" / "fig6d_non_nested_traceability_qc.md"
OUT_PDF = ROOT / "figures" / "fig6d_traceability_matrix.pdf"
OUT_SVG = ROOT / "figures" / "fig6d_traceability_matrix.svg"
OUT_PNG = ROOT / "figures" / "fig6d_traceability_matrix.png"
OUT_PREVIEW = ROOT / "figures" / "_pdf_preview" / "fig6d_traceability_matrix.png"
OUT_GRAY = ROOT / "figures" / "_pdf_preview" / "fig6d_traceability_matrix_grayscale.png"

# ---------------------------------------------------------------------------
# Geometry / style
# ---------------------------------------------------------------------------
WIDTH_MM = 90.0
HEIGHT_MM = 74.0

ANCHOR_ORDER = [
    "PRELP", "OMD",
    "IDH2", "NDUFA8",
    "BMP10", "SMAD9", "SMAD6",
    "NTRK2", "NTF3", "FRS2",
    "PITX2", "KCNN3",
]
CALIBRATORS = {"PITX2", "KCNN3"}
N_ANCHORS = 12
PP4_GATE = 0.70

STATUS_PRESENT = "present"
STATUS_NO = "evaluated_no_support"
STATUS_NA = "not_evaluated"
ALLOWED_STATUS = {STATUS_PRESENT, STATUS_NO, STATUS_NA}

# Charcoal core calibration; muted category colors for external only
CORE_FILL = "#3A3A3A"
OPEN_EDGE = "#8E8E8E"       # ~10–15% stronger open-circle contrast
NA_COLOR = "#A8A8A8"
ACCENT_COMMON = "#7A8FA3"   # muted slate blue
ACCENT_REG = "#6A9A94"      # muted teal
ACCENT_RARE = "#A08090"     # muted mauve
BAND_CAL = "#F7F6F3"        # lighter calibration-control wash
SMAD9_HL = "#FBF7F1"        # ~10–15% lighter than prior warm tint
ROW_LINE = "#EDEDED"
SEP_CAL = "#9A9590"
RULE = "#D8D5D0"
BRACKET = "#D8D4CC"
GROUP_RULE = "#B8B3AC"
CORE_EXT_SEP = "#E2DFD9"

SYM_YES = "●"
SYM_NO = "○"
SYM_NA = "—"
SYM_R_MM = 0.92             # identical diameter across six evidence columns

# Shared Figure 6 A–D typography
FS_LETTER = 10.8
FS_TITLE = 8.7
FS_SUB = 6.6
FS_BOUND = 6.2
FS_GROUP = 5.35
FS_COL = 6.0
FS_COUNT = 5.9
FS_FAM = 5.9
FS_GENE = 6.55
FS_SYM = 7.2
FS_FOOT = 5.5
FS_LEG = 5.5
FS_CAL_LAB = 5.45
from figure6_composite_geometry import (  # noqa: E402
    CD_BOUND_DY,
    CD_FOOTER_BASELINE,
    CD_LETTER_Y,
    CD_SUB_DY,
)
EXPECTED_SOURCE_SHA256 = (
    "9a8c186d927e2f66996172e80ffb3220d89c05b20766aac3eb9106f5397f2bce"
)

FORBIDDEN_READ_COLS = {
    "Path_A", "path_a", "PriorityScore", "priority_score",
    "support_count", "G_tier", "g_tier", "AtlasPriority",
    "highest_tier", "evidence_score", "evidence_count",
}
FORBIDDEN_VISIBLE = [
    "G0", "G1", "G2", "G3", "G4",
    "highest tier", "evidence score", "evidence count",
    "support count", "Miyazawa", "miyazawa", "mr_",
    "NCP_v3", "ncp_v3",
    "fig6a_anchor", "fig6b_posterior", "external_af_gene_sets",
    "figure6_mr", "pipeline_wide",
]

COMMON_ROLE_TOKENS = (
    "common-risk", "common risk", "prioritized", "risk-gene", "risk gene",
    "common-risk / prioritized",
)


def _stop(msg: str) -> None:
    raise SystemExit(f"STOP: {msg}")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _fixed_bbox(w: float = WIDTH_MM, h: float = HEIGHT_MM) -> Bbox:
    return Bbox.from_bounds(0.0, 0.0, w / 25.4, h / 25.4)


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


def _status_from_bool(present: bool) -> str:
    return STATUS_PRESENT if present else STATUS_NO


# ---------------------------------------------------------------------------
# Load / build source
# ---------------------------------------------------------------------------
def load_manifest() -> pd.DataFrame:
    if not MANIFEST_EXT.is_file():
        _stop(f"missing external-set manifest: {MANIFEST_EXT}")
    m = pd.read_csv(MANIFEST_EXT, sep="\t")
    m = m[m["included_in_6C"] == True].copy()  # noqa: E712
    m = m.sort_values("display_order")
    if len(m) != 4:
        _stop(f"expected 4 included external sets, got {len(m)}")
    for col in (
        "source_id", "evidence_role", "display_author_year",
        "gene_list_sha256", "n_genes_normalized", "citation",
    ):
        if col not in m.columns:
            _stop(f"manifest missing column: {col}")
    if m["gene_list_sha256"].isna().any() or (
        m["gene_list_sha256"].astype(str).str.len() < 32
    ).any():
        _stop("any external set lacks provenance SHA256")
    yuan = m[m["source_id"].str.contains("yuan", case=False)]
    roselli = m[m["source_id"].str.contains("roselli", case=False)]
    if len(yuan) != 1 or len(roselli) != 1:
        _stop("Yuan and Roselli cannot be separated")
    if yuan.iloc[0]["source_id"] == roselli.iloc[0]["source_id"]:
        _stop("Yuan and Roselli conflated to one source_id")
    return m.reset_index(drop=True)


def _role_category(evidence_role: str) -> str:
    role = str(evidence_role).strip().lower()
    if "rare" in role or "functional" in role or "mendelian" in role:
        return "rare_functional"
    if "regulatory" in role:
        return "regulatory"
    if any(tok in role for tok in COMMON_ROLE_TOKENS) or "priorit" in role:
        return "common_risk"
    _stop(f"unclassified evidence_role: {evidence_role}")


def load_external_sets(manifest: pd.DataFrame, universe: set[str]) -> dict[str, set[str]]:
    ext = pd.read_csv(EXT, sep="\t")
    for col in ext.columns:
        if col in FORBIDDEN_READ_COLS:
            _stop(f"forbidden field in external gene sets: {col}")
    ext_u = ext[ext["in_universe"] == True]  # noqa: E712
    out: dict[str, set[str]] = {}
    for row in manifest.itertuples(index=False):
        sid = str(row.source_id)
        genes_all = sorted(
            ext.loc[ext["source_id"] == sid, "gene_symbol"].astype(str).unique()
        )
        if not genes_all:
            _stop(f"external set missing: {sid}")
        gene_sha = _sha256_text("\n".join(genes_all))
        if gene_sha != str(row.gene_list_sha256):
            _stop(
                f"gene-list SHA256 mismatch for {sid}: "
                f"manifest={row.gene_list_sha256} live={gene_sha}"
            )
        genes_u = set(
            ext_u.loc[ext_u["source_id"] == sid, "gene_symbol"].astype(str)
        ) & universe
        if len(genes_u) != int(row.n_genes_normalized):
            _stop(
                f"normalized gene count mismatch for {sid}: "
                f"manifest={row.n_genes_normalized} live={len(genes_u)}"
            )
        out[sid] = genes_u
    return out


def _optional_egene_p() -> dict[str, float]:
    """Source-table enrichment only; never used for status."""
    if not PIPELINE_WIDE.is_file():
        return {}
    wide = pd.read_csv(PIPELINE_WIDE, sep="\t", usecols=lambda c: c in {"anchor", "lead_eqtl_p"})
    if "anchor" not in wide.columns or "lead_eqtl_p" not in wide.columns:
        return {}
    out: dict[str, float] = {}
    for r in wide.itertuples(index=False):
        g = str(r.anchor)
        if g in ANCHOR_ORDER and pd.notna(r.lead_eqtl_p):
            out[g] = float(r.lead_eqtl_p)
    return out


def build_source(
    manifest: pd.DataFrame,
    ext_sets: dict[str, set[str]],
) -> tuple[pd.DataFrame, dict[str, int]]:
    if not FIG6A.is_file():
        _stop(f"missing Figure 6A source: {FIG6A}")
    if not FIG6B.is_file():
        _stop(f"missing Figure 6B source: {FIG6B}")

    a = pd.read_csv(FIG6A, sep="\t")
    b = pd.read_csv(FIG6B, sep="\t")
    for col in a.columns.tolist() + b.columns.tolist():
        if col in FORBIDDEN_READ_COLS:
            _stop(f"forbidden field in frozen source: {col}")

    if list(a["gene"]) != ANCHOR_ORDER:
        _stop(f"Figure 6A row order differs from locked order: {list(a['gene'])}")
    if list(b["gene"]) != ANCHOR_ORDER:
        _stop(f"Figure 6B row order differs from locked order: {list(b['gene'])}")
    if len(a) != N_ANCHORS or a["gene"].nunique() != N_ANCHORS:
        _stop("anchor duplication or missing rows in Figure 6A")
    if len(b) != N_ANCHORS or b["gene"].nunique() != N_ANCHORS:
        _stop("anchor duplication or missing rows in Figure 6B")

    # Cross-table PP conflict check
    for g in ANCHOR_ORDER:
        pa = float(a.loc[a["gene"] == g, "coloc_PP4"].iloc[0])
        pb = float(b.loc[b["gene"] == g, "coloc_PP4"].iloc[0])
        if abs(pa - pb) > 1e-9:
            _stop(f"PP4 conflict for {g}: 6A={pa} 6B={pb}")
        p3a = float(a.loc[a["gene"] == g, "coloc_PP3"].iloc[0])
        p3b = float(b.loc[b["gene"] == g, "coloc_PP3"].iloc[0])
        if abs(p3a - p3b) > 1e-9:
            _stop(f"PP3 conflict for {g}: 6A={p3a} 6B={p3b}")

    # Category → source_ids
    cat_ids: dict[str, list[str]] = {
        "common_risk": [], "regulatory": [], "rare_functional": [],
    }
    id_meta = {}
    for row in manifest.itertuples(index=False):
        cat = _role_category(str(row.evidence_role))
        cat_ids[cat].append(str(row.source_id))
        id_meta[str(row.source_id)] = row

    if len(cat_ids["common_risk"]) < 2:
        _stop("common-risk category must retain separate Yuan and Roselli sets")
    egene_p = _optional_egene_p()

    sha_a = _sha256(FIG6A)
    sha_b = _sha256(FIG6B)
    sha_m = _sha256(MANIFEST_EXT)
    sha_e = _sha256(EXT)
    source_paths = ";".join([
        str(FIG6A.relative_to(ROOT)),
        str(FIG6B.relative_to(ROOT)),
        str(MANIFEST_EXT.relative_to(ROOT)),
        str(EXT.relative_to(ROOT)),
    ])
    source_shas = ";".join([sha_a, sha_b, sha_m, sha_e])

    rows = []
    for order, g in enumerate(ANCHOR_ORDER, start=1):
        ra = a.loc[a["gene"] == g].iloc[0]
        rb = b.loc[b["gene"] == g].iloc[0]
        is_cal = g in CALIBRATORS
        fam = str(ra["group_compact"])

        gwas_present = _bool(ra["local_AF_GWAS_support"])
        egene_present = _bool(ra["atrial_eGene_bonferroni"])
        pp4 = float(ra["coloc_PP4"])
        pp3 = float(ra["coloc_PP3"])
        # Prefer 6A boolean; assert against PP4 gate
        coloc_present = _bool(ra["PP4_ge_0.70"])
        if coloc_present != (pp4 >= PP4_GATE):
            _stop(f"coloc status conflict for {g}")
        if _bool(rb["pp4_support"]) != coloc_present:
            _stop(f"Figure 6B pp4_support conflicts for {g}")

        def _ext_status(cat: str) -> tuple[str, str, str]:
            matched = []
            for sid in cat_ids[cat]:
                if g in ext_sets[sid]:
                    meta = id_meta[sid]
                    matched.append(
                        f"{sid}|{meta.display_author_year}|{g}|{meta.gene_list_sha256}"
                    )
            if matched:
                return STATUS_PRESENT, ";".join(sid for sid in cat_ids[cat] if g in ext_sets[sid]), ";".join(matched)
            # Successfully compared against all category sets
            return STATUS_NO, ";".join(cat_ids[cat]), ""

        cr_st, cr_src, cr_match = _ext_status("common_risk")
        rg_st, rg_src, rg_match = _ext_status("regulatory")
        rf_st, rf_src, rf_match = _ext_status("rare_functional")

        note_parts = []
        if g == "SMAD9":
            note_parts.append(
                "atrial eGene present without same-variant colocalization"
            )
        if is_cal:
            note_parts.append("AF calibration control")
        if coloc_present is False and pd.notna(pp4):
            note_parts.append(f"PP4={pp4:.6f} < {PP4_GATE:g}")

        rows.append({
            "anchor_order": order,
            "family_or_role": fam,
            "gene": g,
            "is_calibration_control": is_cal,
            "local_gwas_status": _status_from_bool(gwas_present),
            "local_gwas_definition": (
                "frozen ±1 Mb window contains genome-wide significant AF "
                "association (Figure 6A local_AF_GWAS_support)"
            ),
            "local_gwas_sources": str(FIG6A.relative_to(ROOT)),
            "atrial_egene_status": _status_from_bool(egene_present),
            "atrial_egene_definition": (
                "Bonferroni atrial eGene (Figure 6A atrial_eGene_bonferroni); "
                "not the legacy lead-p < 5e-8 rule"
            ),
            "atrial_egene_p": egene_p.get(g, np.nan),
            "coloc_status": _status_from_bool(coloc_present),
            "coloc_PP3": pp3,
            "coloc_PP4": pp4,
            "coloc_definition": f"PP4 >= {PP4_GATE:g} (frozen; Figure 6A/6B)",
            "common_risk_status": cr_st,
            "common_risk_sources": cr_src,
            "common_risk_matched_sets": cr_match,
            "regulatory_status": rg_st,
            "regulatory_sources": rg_src,
            "regulatory_matched_sets": rg_match,
            "rare_functional_status": rf_st,
            "rare_functional_sources": rf_src,
            "rare_functional_matched_sets": rf_match,
            "evaluation_note": "; ".join(note_parts),
            "source_paths": source_paths,
            "source_sha256s": source_shas,
        })

    df = pd.DataFrame(rows)

    # Hard scientific assertions
    if len(df) != N_ANCHORS or df["gene"].nunique() != N_ANCHORS:
        _stop("exactly 12 unique anchors required")
    if list(df["gene"]) != ANCHOR_ORDER:
        _stop("row order must match Figures 6A/6B")
    n_fw = int((~df["is_calibration_control"]).sum())
    n_cal = int(df["is_calibration_control"].sum())
    if n_fw != 10 or n_cal != 2:
        _stop(f"framework={n_fw} calibrators={n_cal}")

    counts = {
        "local_gwas": int((df["local_gwas_status"] == STATUS_PRESENT).sum()),
        "atrial_egene": int((df["atrial_egene_status"] == STATUS_PRESENT).sum()),
        "coloc": int((df["coloc_status"] == STATUS_PRESENT).sum()),
        "common_risk": int((df["common_risk_status"] == STATUS_PRESENT).sum()),
        "regulatory": int((df["regulatory_status"] == STATUS_PRESENT).sum()),
        "rare_functional": int((df["rare_functional_status"] == STATUS_PRESENT).sum()),
    }
    if counts["local_gwas"] != 7:
        _stop(f"local AF GWAS present != 7/12 (got {counts['local_gwas']})")
    if counts["atrial_egene"] != 1:
        _stop(f"atrial eGene present != 1/12 (got {counts['atrial_egene']})")
    if counts["coloc"] != 0:
        _stop(f"coloc present != 0/12 (got {counts['coloc']})")

    smad = df.loc[df["gene"] == "SMAD9"].iloc[0]
    if smad["atrial_egene_status"] != STATUS_PRESENT:
        _stop("SMAD9 atrial eGene must be present")
    if smad["coloc_status"] != STATUS_NO:
        _stop("SMAD9 coloc must be evaluated_no_support")
    for g in ("PITX2", "KCNN3"):
        r = df.loc[df["gene"] == g].iloc[0]
        if not r["is_calibration_control"]:
            _stop(f"{g} must be calibration control")
        if r["coloc_status"] != STATUS_NO:
            _stop(f"{g} coloc must be evaluated_no_support")
        if r["family_or_role"] != "AF calibration controls":
            _stop(f"{g} family_or_role must be AF calibration controls")

    for r in df.itertuples(index=False):
        for st_col, match_col in (
            ("common_risk_status", "common_risk_matched_sets"),
            ("regulatory_status", "regulatory_matched_sets"),
            ("rare_functional_status", "rare_functional_matched_sets"),
        ):
            st = getattr(r, st_col)
            if st not in ALLOWED_STATUS:
                _stop(f"invalid status {st} for {r.gene}.{st_col}")
            if st == STATUS_PRESENT and not str(getattr(r, match_col)).strip():
                _stop(f"present flag without matched set: {r.gene}.{st_col}")

    return df, counts


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------
def _family_groups(df: pd.DataFrame) -> list[tuple[str, list[int]]]:
    groups: list[tuple[str, list[int]]] = []
    cur_lab: str | None = None
    cur_idx: list[int] = []
    for i, r in df.iterrows():
        lab = str(r["family_or_role"])
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


def _draw_status(ax, x: float, y: float, status: str, fill: str) -> None:
    """Fixed-diameter circle symbols for grayscale-stable present/open states."""
    if status == STATUS_PRESENT:
        ax.add_patch(Circle(
            (x, y), SYM_R_MM,
            facecolor=fill, edgecolor=fill, linewidth=0.0,
            zorder=5, clip_on=False,
        ))
    elif status == STATUS_NO:
        ax.add_patch(Circle(
            (x, y), SYM_R_MM,
            facecolor="white", edgecolor=OPEN_EDGE, linewidth=0.95,
            zorder=5, clip_on=False,
        ))
    elif status == STATUS_NA:
        ax.text(
            x, y, SYM_NA, ha="center", va="center",
            fontsize=FS_SYM, color=NA_COLOR, zorder=5, clip_on=False,
        )
    else:
        _stop(f"unrenderable status: {status}")


def _artist_bbox_mm(ax, fig, artist) -> tuple[float, float, float, float]:
    fig.canvas.draw()
    bb = artist.get_window_extent(renderer=fig.canvas.get_renderer())
    inv = ax.transData.inverted()
    x0, y0 = inv.transform((bb.x0, bb.y0))
    x1, y1 = inv.transform((bb.x1, bb.y1))
    return min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)


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

    x0 = 0.85
    # Fixed column allocation for 90 mm composite width
    fam_w = 12.8
    gene_w = 13.5
    core_w = 23.0
    gap_ce = 3.0
    cr_w = 12.2
    reg_w = 11.3
    rare_w = 11.8
    gap_ext = 0.4

    x_fam = x0
    x_gene = x_fam + fam_w
    x_core = x_gene + gene_w
    x_ext0 = x_core + core_w + gap_ce
    x_cr = x_ext0
    x_reg = x_cr + cr_w + gap_ext
    x_rare = x_reg + reg_w + gap_ext
    x_right = x_rare + rare_w
    assert x_right <= WIDTH_MM - 0.05, x_right
    ext_w = x_right - x_ext0
    x_core_ext_sep = x_core + core_w + gap_ce / 2

    cx_fam = x_fam + fam_w / 2
    cx_gene = x_gene + gene_w / 2
    core_pad = 1.8
    g_span = core_w - 2 * core_pad
    cx_g, cx_e, cx_c = [
        x_core + core_pad + (i + 0.5) * (g_span / 3) for i in range(3)
    ]
    cx_cr = x_cr + cr_w / 2
    cx_reg = x_reg + reg_w / 2
    cx_rare = x_rare + rare_w / 2

    # Compact title block — shared C/D baselines
    y = CD_LETTER_Y
    ax.text(
        x0, y, "D", fontsize=FS_LETTER, fontweight="bold",
        ha="left", va="top", color=FG_COLOR, clip_on=True,
    )
    ax.text(
        x0 + 3.4, y, "Non-nested genetic traceability",
        fontsize=FS_TITLE, fontweight="bold",
        ha="left", va="top", color=FG_COLOR, clip_on=True,
    )
    y -= CD_SUB_DY
    ax.text(
        x0, y, "Independent flags across 12 prespecified anchors",
        fontsize=FS_SUB, ha="left", va="top", color="#5A5A5A", clip_on=True,
    )
    y -= CD_BOUND_DY
    ax.text(
        x0, y, "No composite score or ordinal evidence tier.",
        fontsize=FS_BOUND, ha="left", va="top", color=MUTED_COLOR, clip_on=True,
    )

    footer_h = 5.9
    head_h = 12.4   # three distinct baselines + underlines
    cal_label_h = 1.95
    body_top = y - 1.0 - head_h
    body_bot = footer_h
    usable = body_top - body_bot - cal_label_h
    if usable < 36:
        _stop(f"row region too short after compact title/header: {usable:.1f} mm")
    row_h = usable / N_ANCHORS
    n_framework = 10
    row_y: list[float] = []
    for i in range(n_framework):
        row_y.append(body_top - (i + 0.5) * row_h)
    y_fw_bot = body_top - n_framework * row_h
    y_cal_label = y_fw_bot - cal_label_h / 2
    y_cal_top = y_fw_bot - cal_label_h
    for i in range(2):
        row_y.append(y_cal_top - (i + 0.5) * row_h)
    matrix_bot = row_y[-1] - row_h / 2

    # Header wash
    ax.add_patch(Rectangle(
        (x0, body_top), x_right - x0, head_h,
        facecolor="#F8F7F5", edgecolor="none", zorder=1,
    ))
    ax.plot(
        [x0, x_right], [body_top, body_top],
        color=RULE, lw=0.45, zorder=2,
    )

    # Level 1 — lighter group headings
    y_l1 = body_top + head_h - 2.0
    t_core = ax.text(
        x_core + core_w / 2, y_l1, "CORE CALIBRATION",
        ha="center", va="center", fontsize=FS_GROUP, fontweight="semibold",
        color="#4A4A4A", zorder=4, clip_on=True,
    )
    t_ext = ax.text(
        x_ext0 + ext_w / 2, y_l1, "EXTERNAL TRACEABILITY",
        ha="center", va="center", fontsize=FS_GROUP, fontweight="semibold",
        color="#4A4A4A", zorder=4, clip_on=True,
    )
    # Thinner / lighter underline per group
    y_rule = y_l1 - 1.4
    ax.plot(
        [x_core + 0.8, x_core + core_w - 0.8], [y_rule, y_rule],
        color=GROUP_RULE, lw=0.55, solid_capstyle="round", zorder=3,
    )
    ax.plot(
        [x_ext0 + 0.6, x_right - 0.6], [y_rule, y_rule],
        color=GROUP_RULE, lw=0.55, solid_capstyle="round", zorder=3,
    )

    # Level 2 — evidence labels (own baseline)
    y_l2 = body_top + 5.65
    t_fam = ax.text(
        cx_fam, y_l2, "Family /\nrole",
        ha="center", va="center", fontsize=FS_COL - 0.2, fontweight="bold",
        color=MUTED_COLOR, linespacing=0.92, zorder=4,
    )
    ax.text(
        cx_gene, y_l2, "Anchor",
        ha="center", va="center", fontsize=FS_COL - 0.1, fontweight="bold",
        color=MUTED_COLOR, zorder=4,
    )
    l2_artists = []
    for cx, lab in ((cx_g, "G"), (cx_e, "E"), (cx_c, "C")):
        l2_artists.append(ax.text(
            cx, y_l2, lab, ha="center", va="center",
            fontsize=FS_COL + 0.2, fontweight="bold", color=CORE_FILL, zorder=4,
        ))
    for cx, lab, col in (
        (cx_cr, "Risk /\nprioritized", ACCENT_COMMON),
        (cx_reg, "Regulatory", ACCENT_REG),
        (cx_rare, "Rare /\nfunctional", ACCENT_RARE),
    ):
        l2_artists.append(ax.text(
            cx, y_l2, lab, ha="center", va="center",
            fontsize=FS_COL - 0.3, fontweight="bold", color=col,
            linespacing=0.92, zorder=4,
        ))

    # Level 3 — counts (own baseline)
    y_l3 = body_top + 2.05
    for cx, n in (
        (cx_g, counts["local_gwas"]),
        (cx_e, counts["atrial_egene"]),
        (cx_c, counts["coloc"]),
        (cx_cr, counts["common_risk"]),
        (cx_reg, counts["regulatory"]),
        (cx_rare, counts["rare_functional"]),
    ):
        ax.text(
            cx, y_l3, f"{n}/{N_ANCHORS}",
            ha="center", va="center", fontsize=FS_COUNT,
            color=FG_COLOR, zorder=4,
        )

    # Very light vertical separator between C and Risk / prioritized
    ax.plot(
        [x_core_ext_sep, x_core_ext_sep], [matrix_bot, body_top + head_h - 0.3],
        color=CORE_EXT_SEP, lw=0.45, zorder=2, solid_capstyle="butt",
    )

    # Measure group-heading bboxes for hard layout assertions
    core_bb = _artist_bbox_mm(ax, fig, t_core)
    ext_bb = _artist_bbox_mm(ax, fig, t_ext)
    fam_bb = _artist_bbox_mm(ax, fig, t_fam)
    header_spans = {
        "core_group": core_bb,
        "ext_group": ext_bb,
        "family_role": fam_bb,
        "l2": [_artist_bbox_mm(ax, fig, t) for t in l2_artists],
        "core_span": (x_core, x_core + core_w),
        "ext_span": (x_ext0, x_right),
        "fam_span": (x_fam, x_fam + fam_w),
        "col_bounds": [
            (x_core + core_pad, x_core + core_pad + g_span / 3),
            (x_core + core_pad + g_span / 3, x_core + core_pad + 2 * g_span / 3),
            (x_core + core_pad + 2 * g_span / 3, x_core + core_pad + g_span),
            (x_cr, x_cr + cr_w),
            (x_reg, x_reg + reg_w),
            (x_rare, x_rare + rare_w),
        ],
        "core_ext_sep_x": x_core_ext_sep,
    }

    # Data rows
    for i, r in df.iterrows():
        idx = int(i)
        yi = row_y[idx]
        y_bot = yi - row_h / 2
        if r["is_calibration_control"]:
            ax.add_patch(Rectangle(
                (x0, y_bot), x_right - x0, row_h,
                facecolor=BAND_CAL, edgecolor="none", zorder=0,
            ))
        elif r["gene"] == "SMAD9":
            ax.add_patch(Rectangle(
                (x0, y_bot), x_right - x0, row_h,
                facecolor=SMAD9_HL, edgecolor="none", zorder=0,
            ))
        ax.plot(
            [x0 + 0.25, x_right - 0.25], [y_bot, y_bot],
            color=ROW_LINE, lw=0.28, zorder=1,
        )
        ax.text(
            cx_gene, yi, r["gene"],
            ha="center", va="center",
            fontsize=FS_GENE, fontweight="semibold",
            color="#4A4038" if r["gene"] == "SMAD9" else FG_COLOR,
            zorder=4,
        )
        _draw_status(ax, cx_g, yi, r["local_gwas_status"], CORE_FILL)
        _draw_status(ax, cx_e, yi, r["atrial_egene_status"], CORE_FILL)
        _draw_status(ax, cx_c, yi, r["coloc_status"], CORE_FILL)
        _draw_status(ax, cx_cr, yi, r["common_risk_status"], ACCENT_COMMON)
        _draw_status(ax, cx_reg, yi, r["regulatory_status"], ACCENT_REG)
        _draw_status(ax, cx_rare, yi, r["rare_functional_status"], ACCENT_RARE)

    # Stronger separator + calibration-control label (not a data row)
    y_div = y_fw_bot
    ax.plot(
        [x0 + 0.2, x_right - 0.2], [y_div, y_div],
        color=SEP_CAL, lw=0.75, zorder=3,
    )
    ax.add_patch(Rectangle(
        (x0, y_cal_top), x_gene + gene_w - x0, cal_label_h,
        facecolor=BAND_CAL, edgecolor="none", zorder=0,
    ))
    ax.text(
        (x_fam + x_gene + gene_w) / 2, y_cal_label, "AF calibration controls",
        ha="center", va="center", fontsize=FS_CAL_LAB, fontweight="semibold",
        color=MUTED_COLOR, zorder=4, clip_on=True,
    )

    # Family brackets for framework groups only (light grouping guides)
    for lab, idxs in _family_groups(df):
        if lab == "AF calibration controls":
            continue
        y0 = row_y[idxs[-1]] - row_h / 2 + 0.2
        y1 = row_y[idxs[0]] + row_h / 2 - 0.2
        bx = x_fam + 0.95
        ax.plot([bx, bx], [y0, y1], color=BRACKET, lw=0.64,
                solid_capstyle="round", zorder=2)
        ax.plot([bx, bx + 0.55], [y0, y0], color=BRACKET, lw=0.6, zorder=2)
        ax.plot([bx, bx + 0.55], [y1, y1], color=BRACKET, lw=0.6, zorder=2)
        y_mid = 0.5 * (row_y[idxs[0]] + row_y[idxs[-1]])
        ax.text(
            cx_fam + 0.55, y_mid, lab,
            ha="center", va="center", fontsize=FS_FAM,
            fontweight="semibold", color=MUTED_COLOR, zorder=4,
        )

    # Footer — two panel-internal lines only (shared C/D baseline)
    foot_y = CD_FOOTER_BASELINE
    ax.text(
        x0, foot_y,
        "G = local AF GWAS · E = atrial eGene · C = PP4 ≥ 0.70",
        ha="left", va="top", fontsize=FS_FOOT, color=MUTED_COLOR, clip_on=True,
    )
    ax.text(
        x0, foot_y - 2.2,
        f"{SYM_YES} present · {SYM_NO} evaluated, criterion not met/"
        f"no membership · {SYM_NA} not evaluated",
        ha="left", va="top", fontsize=FS_FOOT, color=MUTED_COLOR, clip_on=True,
    )

    fig._fig6d_layout = {  # type: ignore[attr-defined]
        "row_y": row_y,
        "row_h": row_h,
        "centers": {
            "G": cx_g, "E": cx_e, "C": cx_c,
            "common": cx_cr, "reg": cx_reg, "rare": cx_rare,
            "gene": cx_gene, "fam": cx_fam,
        },
        "y_div": y_div,
        "y_cal_label": y_cal_label,
        "body_top": body_top,
        "body_bot": body_bot,
        "matrix_bot": matrix_bot,
        "x_right": x_right,
        "header_spans": header_spans,
        "baselines": {"l1": y_l1, "l2": y_l2, "l3": y_l3, "rule": y_rule},
        "sym_r_mm": SYM_R_MM,
    }
    return fig


# ---------------------------------------------------------------------------
# QC / outputs
# ---------------------------------------------------------------------------
def _pdf_page_mm(pdf: Path) -> tuple[float, float]:
    import subprocess
    info = subprocess.run(
        ["pdfinfo", str(pdf)], check=True, capture_output=True, text=True,
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    if not match:
        raise RuntimeError("Unable to read PDF page size")
    return float(match.group(1)) * 25.4 / 72, float(match.group(2)) * 25.4 / 72


def _svg_mm(svg: Path) -> tuple[float, float]:
    text = svg.read_text(encoding="utf-8")
    m = re.search(r'width="([0-9.]+)(?:mm|pt)?"', text)
    h = re.search(r'height="([0-9.]+)(?:mm|pt)?"', text)
    if not m or not h:
        # viewBox fallback
        vb = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', text)
        if not vb:
            raise RuntimeError("Unable to read SVG size")
        return float(vb.group(1)) * 25.4 / 72, float(vb.group(2)) * 25.4 / 72
    w, hh = float(m.group(1)), float(h.group(1))
    # matplotlib often writes inches as unitless pt-like numbers in svg
    if w < 20:  # inches
        return w * 25.4, hh * 25.4
    if "mm" in (m.group(0) + h.group(0)):
        return w, hh
    return w * 25.4 / 72, hh * 25.4 / 72


def _png_mm(png: Path, dpi: int = 600) -> tuple[float, float, int, int]:
    with png.open("rb") as fh:
        sig = fh.read(8)
        if sig[:4] != b"\x89PNG":
            raise RuntimeError(f"not a PNG: {png}")
        fh.read(8)  # length + IHDR
        wh = fh.read(8)
    ww, hh = struct.unpack(">II", wh)
    return ww * 25.4 / dpi, hh * 25.4 / dpi, ww, hh


def _min_svg_font(svg: Path) -> float:
    text = svg.read_text(encoding="utf-8")
    sizes = [float(x) for x in re.findall(r"font-size:\s*([0-9.]+)px", text)]
    sizes += [float(x) for x in re.findall(r"font-size:\s*([0-9.]+)pt", text)]
    return min(sizes) if sizes else 99.0


def _assert_header_layout(layout: dict) -> None:
    hs = layout["header_spans"]
    core = hs["core_group"]
    ext = hs["ext_group"]
    fam = hs["family_role"]
    c0, c1 = hs["core_span"]
    e0, e1 = hs["ext_span"]
    f0, f1 = hs["fam_span"]
    # Fully visible within group spans; no clipping at panel edge
    assert core[0] >= c0 - 0.2 and core[1] <= c1 + 0.2, core
    assert ext[0] >= e0 - 0.2 and ext[1] <= e1 + 0.2, ext
    assert ext[1] <= WIDTH_MM - 0.2, ext
    assert fam[0] >= f0 - 0.25 and fam[1] <= f1 + 0.25, fam
    # Group headings do not overlap; extra core↔external gutter
    assert core[1] + 0.35 <= ext[0], (core, ext)
    assert (e0 - c1) >= 1.4, (c1, e0)
    # Distinct baselines
    bl = layout["baselines"]
    assert bl["l1"] > bl["rule"] > bl["l2"] > bl["l3"]
    assert abs(bl["l1"] - bl["l2"]) >= 2.0
    assert abs(bl["l2"] - bl["l3"]) >= 2.0
    # Level-2 labels stay inside their columns and do not overlap neighbors
    spans = hs["l2"]
    bounds = hs["col_bounds"]
    assert len(spans) == 6 and len(bounds) == 6
    for i, ((x0, x1, _y0, _y1), (b0, b1)) in enumerate(zip(spans, bounds)):
        assert x0 >= b0 - 0.4, (i, x0, b0)
        assert x1 <= b1 + 0.4, (i, x1, b1)
        if i > 0:
            assert x0 >= spans[i - 1][1] + 0.12, (i, spans[i - 1], (x0, x1))
    # Core/external separator sits in the gutter, not on symbols
    sep_x = hs["core_ext_sep_x"]
    assert c1 < sep_x < e0, (c1, sep_x, e0)
    assert abs(sep_x - spans[2][1]) > 0.8
    assert abs(spans[3][0] - sep_x) > 0.8


def run_qc(
    df: pd.DataFrame,
    counts: dict[str, int],
    layout: dict,
    manifest: pd.DataFrame,
) -> dict:
    assert list(df["gene"]) == ANCHOR_ORDER
    assert len(df) == 12
    assert int((~df["is_calibration_control"]).sum()) == 10
    assert int(df["is_calibration_control"].sum()) == 2
    assert counts["local_gwas"] == 7
    assert counts["atrial_egene"] == 1
    assert counts["coloc"] == 0
    assert counts["common_risk"] == 2
    assert counts["regulatory"] == 2
    assert counts["rare_functional"] == 0

    smad = df.loc[df["gene"] == "SMAD9"].iloc[0]
    assert smad["local_gwas_status"] == STATUS_NO
    assert smad["atrial_egene_status"] == STATUS_PRESENT
    assert smad["coloc_status"] == STATUS_NO
    for g in ("PITX2", "KCNN3"):
        r = df.loc[df["gene"] == g].iloc[0]
        assert bool(r["is_calibration_control"])
        assert r["local_gwas_status"] == STATUS_PRESENT
        assert r["atrial_egene_status"] == STATUS_NO
        assert r["coloc_status"] == STATUS_NO
        assert r["common_risk_status"] == STATUS_PRESENT
        assert r["regulatory_status"] == STATUS_PRESENT
        assert r["rare_functional_status"] == STATUS_NO

    svg_text = OUT_SVG.read_text(encoding="utf-8")
    visible = re.sub(r"<[^>]+>", " ", svg_text)
    for bad in FORBIDDEN_VISIBLE:
        assert bad not in visible, bad
    for bad in (
        "highest tier", "evidence score", "Path A family", "support count",
        "PUBLISHED", "Common-risk", "G-tier", "G0", "G1", "G2", "G3", "G4",
        "prior.", "Open symbols do not exclude",
    ):
        assert bad not in visible, bad
    assert "Evidence categories are non-nested and are not combined" not in visible

    for r in df.itertuples(index=False):
        for col in (
            "local_gwas_status", "atrial_egene_status", "coloc_status",
            "common_risk_status", "regulatory_status", "rare_functional_status",
        ):
            assert getattr(r, col) in ALLOWED_STATUS

    for g in ANCHOR_ORDER:
        assert g in visible
    assert "CORE CALIBRATION" in visible
    assert "EXTERNAL TRACEABILITY" in visible
    assert "Family" in visible and "role" in visible
    assert "prioritized" in visible
    assert "Regulatory" in visible
    assert "AF calibration controls" in visible
    assert "Non-nested genetic traceability" in visible
    assert "No composite score or ordinal evidence tier." in visible
    assert "Path A" not in visible
    assert "G = local AF GWAS" in visible
    assert f"{SYM_YES} present" in visible or "present ·" in visible
    for name in ("Yuan", "Roselli", "Selewa", "Choi"):
        assert name not in visible, name
    # Matrix headers use G/E/C, not full-word GWAS/eGene/Coloc triad
    assert "eGene" not in visible.replace("atrial eGene", "")
    assert "Coloc" not in visible
    assert abs(layout["sym_r_mm"] - SYM_R_MM) < 1e-9

    pdf_w, pdf_h = _pdf_page_mm(OUT_PDF)
    svg_w, svg_h = _svg_mm(OUT_SVG)
    png_w, png_h, px_w, px_h = _png_mm(OUT_PNG, dpi=600)
    assert abs(pdf_w - WIDTH_MM) < 1.5, pdf_w
    assert 72.0 <= pdf_h <= 76.5, pdf_h
    assert abs(pdf_h - HEIGHT_MM) < 1.5, pdf_h
    assert abs(svg_w - pdf_w) < 1.8 and abs(svg_h - pdf_h) < 1.8
    assert abs(png_w - pdf_w) < 1.8 and abs(png_h - pdf_h) < 1.8

    min_font = _min_svg_font(OUT_SVG)
    assert min_font >= FS_FOOT - 0.15, min_font

    centers = layout["centers"]
    assert abs(centers["G"] - centers["E"]) > 3.0
    assert abs(centers["E"] - centers["C"]) > 3.0
    assert len(layout["row_y"]) == 12
    # Framework rows equal; calibrator rows equal (gap for label between)
    fw = np.array(layout["row_y"][:10])
    cal = np.array(layout["row_y"][10:])
    assert np.allclose(np.diff(fw), np.diff(fw)[0], atol=0.05)
    assert np.allclose(np.diff(cal), np.diff(cal)[0], atol=0.05)
    assert abs(layout["row_h"] - abs(np.diff(fw)[0])) < 0.05
    # Calibration label is not a data-row center
    assert layout["y_cal_label"] > layout["row_y"][10]
    assert layout["y_cal_label"] < layout["row_y"][9]

    _assert_header_layout(layout)

    source_sha = _sha256(OUT_SOURCE)
    assert source_sha == EXPECTED_SOURCE_SHA256, source_sha

    na_cells = []
    for r in df.itertuples(index=False):
        for col in (
            "local_gwas_status", "atrial_egene_status", "coloc_status",
            "common_risk_status", "regulatory_status", "rare_functional_status",
        ):
            if getattr(r, col) == STATUS_NA:
                na_cells.append(f"{r.gene}:{col}")

    return {
        "counts": counts,
        "na_cells": na_cells,
        "pdf_mm": (round(pdf_w, 2), round(pdf_h, 2)),
        "svg_mm": (round(svg_w, 2), round(svg_h, 2)),
        "png_mm": (round(png_w, 2), round(png_h, 2)),
        "png_px": (px_w, px_h),
        "min_font": min_font,
        "source_sha256": source_sha,
        "pdf_sha256": _sha256(OUT_PDF),
        "svg_sha256": _sha256(OUT_SVG),
        "png_sha256": _sha256(OUT_PNG),
        "manifest_sha256": _sha256(MANIFEST_EXT),
        "layout_qc": "PASS",
        "scientific_qc": "PASS",
        "provenance_qc": "PASS",
    }


def write_outputs(
    df: pd.DataFrame,
    counts: dict[str, int],
    manifest: pd.DataFrame,
    qc: dict,
) -> None:
    df.to_csv(OUT_SOURCE, sep="\t", index=False)

    summary = pd.DataFrame([
        {"dimension": "local_AF_GWAS", "present": counts["local_gwas"], "n": N_ANCHORS,
         "definition": "Figure 6A local_AF_GWAS_support"},
        {"dimension": "atrial_eGene", "present": counts["atrial_egene"], "n": N_ANCHORS,
         "definition": "Figure 6A atrial_eGene_bonferroni"},
        {"dimension": "same_variant_coloc", "present": counts["coloc"], "n": N_ANCHORS,
         "definition": f"PP4 >= {PP4_GATE:g}"},
        {"dimension": "common_risk", "present": counts["common_risk"], "n": N_ANCHORS,
         "definition": "membership in COMMON-RISK / PRIORITIZED sets"},
        {"dimension": "regulatory", "present": counts["regulatory"], "n": N_ANCHORS,
         "definition": "membership in REGULATORY sets"},
        {"dimension": "rare_functional", "present": counts["rare_functional"], "n": N_ANCHORS,
         "definition": "membership in RARE / FUNCTIONAL sets"},
    ])
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    man_rows = [
        {"artifact": "publication_pdf", "path": str(OUT_PDF.relative_to(ROOT)),
         "sha256": qc["pdf_sha256"], "geometry_mm": f"{qc['pdf_mm'][0]}x{qc['pdf_mm'][1]}"},
        {"artifact": "publication_svg", "path": str(OUT_SVG.relative_to(ROOT)),
         "sha256": qc["svg_sha256"], "geometry_mm": f"{qc['svg_mm'][0]}x{qc['svg_mm'][1]}"},
        {"artifact": "publication_png_600dpi", "path": str(OUT_PNG.relative_to(ROOT)),
         "sha256": qc["png_sha256"], "geometry_mm": f"{qc['png_mm'][0]}x{qc['png_mm'][1]}"},
        {"artifact": "print_preview_92mm", "path": str(OUT_PREVIEW.relative_to(ROOT)),
         "sha256": _sha256(OUT_PREVIEW), "geometry_mm": f"{WIDTH_MM:.0f}x{HEIGHT_MM:.0f}"},
        {"artifact": "grayscale_preview", "path": str(OUT_GRAY.relative_to(ROOT)),
         "sha256": _sha256(OUT_GRAY), "geometry_mm": ""},
        {"artifact": "source_table", "path": str(OUT_SOURCE.relative_to(ROOT)),
         "sha256": qc["source_sha256"], "geometry_mm": ""},
        {"artifact": "summary_table", "path": str(OUT_SUMMARY.relative_to(ROOT)),
         "sha256": _sha256(OUT_SUMMARY), "geometry_mm": ""},
        {"artifact": "external_set_manifest", "path": str(MANIFEST_EXT.relative_to(ROOT)),
         "sha256": qc["manifest_sha256"], "geometry_mm": ""},
    ]
    pd.DataFrame(man_rows).to_csv(OUT_MANIFEST, sep="\t", index=False)

    matrix_lines = []
    for r in df.itertuples(index=False):
        matrix_lines.append(
            f"| {r.anchor_order} | {r.gene} | {r.family_or_role} | "
            f"{r.local_gwas_status} | {r.atrial_egene_status} | {r.coloc_status} | "
            f"{r.common_risk_status} | {r.regulatory_status} | {r.rare_functional_status} |"
        )

    set_lines = [
        f"- `{r.source_id}` / {r.display_author_year} / {r.evidence_role} / "
        f"n={int(r.n_genes_normalized)} / SHA256=`{r.gene_list_sha256}` / "
        f"citation: {r.citation}"
        for r in manifest.itertuples(index=False)
    ]
    na_lines = qc["na_cells"] or ["(none)"]

    OUT_QC.write_text(
        "\n".join([
            "# Figure 6D non-nested genetic traceability QC",
            "",
            f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
            "",
            "Status: **PASS**",
            "",
            "## Canvas",
            f"- {WIDTH_MM:.0f} × {HEIGHT_MM:.0f} mm (target 90–92 × 72–76)",
            "- layout polish only; frozen source SHA256 unchanged",
            f"- PDF: {qc['pdf_mm'][0]} × {qc['pdf_mm'][1]} mm",
            f"- SVG: {qc['svg_mm'][0]} × {qc['svg_mm'][1]} mm",
            f"- PNG: {qc['png_mm'][0]} × {qc['png_mm'][1]} mm "
            f"({qc['png_px'][0]} × {qc['png_px'][1]} px @ 600 dpi)",
            f"- minimum SVG font-size: {qc['min_font']:.2f} px",
            "",
            "## Runtime present counts",
            f"- Local AF GWAS (G): {counts['local_gwas']}/12",
            f"- Atrial eGene (E): {counts['atrial_egene']}/12",
            f"- Same-variant coloc (C): {counts['coloc']}/12",
            f"- Common-risk: {counts['common_risk']}/12",
            f"- Regulatory: {counts['regulatory']}/12",
            f"- Rare / functional: {counts['rare_functional']}/12",
            "",
            "## 12-anchor status matrix",
            "| order | gene | role | G | E | C | common-risk | regulatory | rare/functional |",
            "|---|---|---|---|---|---|---|---|---|",
            *matrix_lines,
            "",
            "## not_evaluated cells",
            *[f"- {x}" for x in na_lines],
            "",
            "## External source IDs and hashes",
            *set_lines,
            "",
            "## SHA256",
            f"- source: `{qc['source_sha256']}`",
            f"- pdf: `{qc['pdf_sha256']}`",
            f"- svg: `{qc['svg_sha256']}`",
            f"- png: `{qc['png_sha256']}`",
            "",
            "## QC verdicts",
            f"- scientific_qc = {qc['scientific_qc']}",
            f"- provenance_qc = {qc['provenance_qc']}",
            f"- layout_qc = {qc['layout_qc']}",
            "",
            "## Scope notes",
            "- No G-tier / highest-tier / evidence-score fields.",
            "- No Path A, PriorityScore, or support-count fields read.",
            "- No coloc rerun; PP4 from frozen Figure 6A/6B sources.",
            "- Yuan and Roselli retained as separate common-risk source records.",
            "- Formal Figure 6 A–D composite routing not updated.",
            "",
        ]),
        encoding="utf-8",
    )


def main() -> None:
    universe = set(pd.read_csv(UNI, sep="\t")["gene"].astype(str))
    manifest = load_manifest()
    ext_sets = load_external_sets(manifest, universe)
    df, counts = build_source(manifest, ext_sets)

    fig = render(df, counts)
    layout = dict(fig._fig6d_layout)  # type: ignore[attr-defined]
    save_figure_bundle(
        fig, OUT_PDF, dpi=600,
        bbox_inches=_fixed_bbox(), pad_inches=0,
    )
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")
    print(f"wrote {OUT_PNG.relative_to(ROOT)}")

    fig = render(df, counts)
    layout = dict(fig._fig6d_layout)  # type: ignore[attr-defined]
    OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
    OUT_PREVIEW.parent.mkdir(parents=True, exist_ok=True)
    bb = _fixed_bbox()
    fig.savefig(
        OUT_SVG, format="svg", bbox_inches=bb, pad_inches=0,
        facecolor="white", edgecolor="none",
    )
    fig.savefig(
        OUT_PREVIEW, format="png", dpi=150, bbox_inches=bb, pad_inches=0,
        facecolor="white", edgecolor="none",
    )
    fig.savefig(
        OUT_GRAY, format="png", dpi=150, bbox_inches=bb, pad_inches=0,
        facecolor="white", edgecolor="none",
    )
    plt.close(fig)
    try:
        from PIL import Image
        Image.open(OUT_GRAY).convert("L").save(OUT_GRAY)
    except Exception:
        pass

    # Write source before QC (QC hashes source)
    df.to_csv(OUT_SOURCE, sep="\t", index=False)
    qc = run_qc(df, counts, layout, manifest)
    write_outputs(df, counts, manifest, qc)

    print(f"wrote {OUT_SVG.relative_to(ROOT)}")
    print(f"wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"wrote {OUT_MANIFEST.relative_to(ROOT)}")
    print(f"wrote {OUT_QC.relative_to(ROOT)}")
    print(f"canvas_mm: {WIDTH_MM:.0f} x {HEIGHT_MM:.0f}")
    print(
        "counts: "
        f"G={counts['local_gwas']}/12 "
        f"E={counts['atrial_egene']}/12 "
        f"C={counts['coloc']}/12 "
        f"common={counts['common_risk']}/12 "
        f"reg={counts['regulatory']}/12 "
        f"rare={counts['rare_functional']}/12"
    )
    print("QC: PASS")
    print(f"SOURCE_SHA256={qc['source_sha256']}")
    print(f"PDF_SHA256={qc['pdf_sha256']}")
    print(f"SVG_SHA256={qc['svg_sha256']}")
    print(f"PNG_SHA256={qc['png_sha256']}")


if __name__ == "__main__":
    main()
