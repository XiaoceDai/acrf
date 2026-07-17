#!/usr/bin/env python3
"""
Figure 6C — Family overlap with external AF gene sets (publication panel).

Sparse member-level overlap between the 16 frozen NCP families (Figure 1D
order) and prespecified external AF gene sets. Genetic-traceability context
layer only — not a Path A gate.

Inputs (frozen; read-only):
  results/provenance/figure6_external_gene_sets_manifest.tsv
  results/ncp_v3_family_leading_edge_genes.tsv
  results/ncp_v3_gene_universe.tsv
  results/ncp_v3_family_labels.tsv
  results/ncp_v3_canonical_display_order.tsv
  results/external_af_gene_sets.tsv
  results/fig6c_family_af_overlap_source.tsv  (locked stats; remapped if needed)

Outputs:
  results/fig6c_family_af_overlap_manifest.tsv
  results/fig6c_family_af_overlap_qc.md
  figures/fig6c_family_af_overlap.{pdf,svg,png}
  figures/figS_family_af_overlap_full_matrix.{pdf,svg,png}  (optional full 16×4)
  figures/_pdf_preview/fig6c_family_af_overlap.png

Frozen source (not rewritten):
  results/fig6c_family_af_overlap_source.tsv
"""
from __future__ import annotations

import hashlib
import re
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.transforms import Bbox

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (  # noqa: E402
    FG_COLOR,
    FONT_FLOOR_HARD,
    FS_TICK,
    LA_DEEP,
    MUTED_COLOR,
    RA_DEEP,
    set_nature_style,
)
from _figure_export import save_figure_bundle  # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
MANIFEST_EXT = ROOT / "results" / "provenance" / "figure6_external_gene_sets_manifest.tsv"
FLE = ROOT / "results" / "ncp_v3_family_leading_edge_genes.tsv"
UNI = ROOT / "results" / "ncp_v3_gene_universe.tsv"
LABELS = ROOT / "results" / "ncp_v3_family_labels.tsv"
ORDER = ROOT / "results" / "ncp_v3_canonical_display_order.tsv"
EXT = ROOT / "results" / "external_af_gene_sets.tsv"
LOCKED_SOURCE = ROOT / "results" / "fig6c_family_af_overlap_source.tsv"

OUT_SOURCE = ROOT / "results" / "fig6c_family_af_overlap_source.tsv"
OUT_MANIFEST = ROOT / "results" / "fig6c_family_af_overlap_manifest.tsv"
OUT_QC = ROOT / "results" / "fig6c_family_af_overlap_qc.md"
OUT_PDF = ROOT / "figures" / "fig6c_family_af_overlap.pdf"
OUT_SVG = ROOT / "figures" / "fig6c_family_af_overlap.svg"
OUT_PNG = ROOT / "figures" / "fig6c_family_af_overlap.png"
OUT_PREVIEW = ROOT / "figures" / "_pdf_preview" / "fig6c_family_af_overlap.png"
OUT_GRAY = ROOT / "figures" / "_pdf_preview" / "fig6c_family_af_overlap_grayscale.png"
OUT_SUPP_PDF = ROOT / "figures" / "figS_family_af_overlap_full_matrix.pdf"
OUT_SUPP_SVG = ROOT / "figures" / "figS_family_af_overlap_full_matrix.svg"
OUT_SUPP_PNG = ROOT / "figures" / "figS_family_af_overlap_full_matrix.png"

# ---------------------------------------------------------------------------
# Geometry — sparse main-text panel (3 nonzero families only)
# ---------------------------------------------------------------------------
WIDTH_MM = 84.0
HEIGHT_MM = 74.0  # composite C/D shared outer height
MARGIN_L = 1.3
MARGIN_R = 1.0
DIR_W = 3.6          # small LA/RA tag column
ID_W = 5.2
LABEL_W = 24.0
GUTTER_W = 2.0
COL_GAP = 1.4
ROW_H_MM = 7.5       # compact fixed row height (~20–25% shorter)
HEADER_FONT = "Liberation Sans Narrow"
RULE = "#D0CDC8"
ROW_RULE = "#EDEBE7"
COL_RULE = "#E8E6E2"
AXIS_SEP = "#B8B4AC"   # modest separator between F13 (LA) and F9 (RA)
ID_COLOR = "#4A4A4A"
EXPECTED_ROSELLI_YEAR = 2018
EXPECTED_UNIQUE_GENES = {"CAV1", "HEY2"}

# Typography (pt) — shared Figure 6 A–D system
FS_LETTER = 10.8
FS_TITLE = 8.7
FS_SUB = 6.6
FS_RUNTIME = 6.2
FS_HEADER = 6.2
FS_LABEL = 6.0
FS_FOOT = 5.5
FS_SUMMARY = 6.2
HEADER_PAD_MM = 0.35
TITLE_TO_HEADER_GAP_MM = 3.4  # keep role rules below runtime summary
# Shared C/D title baselines (mm from panel bottom)
from figure6_composite_geometry import (  # noqa: E402
    CD_BOUND_DY,
    CD_FOOTER_BASELINE,
    CD_LETTER_Y,
    CD_SUB_DY,
)


def _fixed_bbox(width_mm: float = WIDTH_MM, height_mm: float = HEIGHT_MM) -> Bbox:
    return Bbox.from_bounds(0.0, 0.0, width_mm / 25.4, height_mm / 25.4)

N_FAMILIES = 16
UNIVERSE_N = 14635
FDR_ALPHA = 0.05

# Main-panel rows: families with ≥1 nonzero overlap, Figure 1D order
DISPLAY_FAMILY_IDS = [3, 13, 9]

# Legacy → current source_id (frozen tables may still use the old Yuan id)
ID_REMAP = {"cross_population_gwas2025": "yuan2025_cross_population"}

PUB_FAMILY_LABEL = {
    1: "Translation / RQC",
    2: "Mito Ca²⁺ transport",
    3: "Adipo-lipid axis",
    4: "Mito proteostasis / FA–BCAA",
    5: "tRNA aminoacylation",
    6: "Mito translation / import",
    7: "OXPHOS / ETC",
    8: "Complex I biogenesis",
    12: "rRNA processing",
    13: "Amino-acid metabolism",
    14: "Translation initiation",
    9: "HS-GAG / bone / GPCR",
    10: "GAG disease / TGF-β",
    11: "BMP signalling",
    15: "NTRK2 / TrkB",
    16: "Keratan sulfate",
}

EXPECTED_NONZERO = {
    (3, "yuan2025_cross_population"): ("CAV1",),
    (3, "roselli2018_loci"): ("CAV1",),
    (13, "yuan2025_cross_population"): ("CAV1",),
    (13, "roselli2018_loci"): ("CAV1",),
    (9, "selewa2023_regulatory"): ("HEY2",),
}

# Evidence-role fills / outlines (identity only — not significance)
ROLE_STYLE = {
    "COMMON-RISK / PRIORITIZED": {
        "fill": "#EDF1F6",
        "edge": "#7A8FA3",
        "rule": "#7A8FA3",
    },
    "REGULATORY": {
        "fill": "#E8F3F1",
        "edge": "#6A9A94",
        "rule": "#6A9A94",
    },
    "RARE / FUNCTIONAL": {
        "fill": "#F4EBEF",
        "edge": "#A08090",
        "rule": "#A08090",
    },
}

FORBIDDEN_READ_COLS = {
    "Path_A", "path_a", "PriorityScore", "priority_score",
    "support_count", "G_tier", "g_tier", "AtlasPriority",
}


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


# ---------------------------------------------------------------------------
# Load / preserve locked statistics
# ---------------------------------------------------------------------------
def load_manifest() -> pd.DataFrame:
    if not MANIFEST_EXT.is_file():
        _stop(f"missing external-set manifest: {MANIFEST_EXT}")
    m = pd.read_csv(MANIFEST_EXT, sep="\t")
    m = m[m["included_in_6C"].astype(str).str.upper().isin(["TRUE", "1", "YES"])].copy()
    m = m.sort_values("display_order").reset_index(drop=True)
    if m.empty:
        _stop("no external sets with included_in_6C = TRUE")
    for col in (
        "source_id", "display_author_year", "display_evidence_label",
        "display_panel_label", "evidence_role", "gene_list_sha256",
        "n_genes_normalized", "source_publication", "publication_year",
    ):
        if col not in m.columns:
            _stop(f"manifest missing column: {col}")
    # Panel must not use long evidence phrases
    forbidden_panel = {
        "Cross-population GWAS", "Gene prioritization",
        "Regulatory mapping", "Rare coding and SV",
    }
    if set(m["display_panel_label"].astype(str)) & forbidden_panel:
        _stop("display_panel_label still contains long legend phrases")
    if m["gene_list_sha256"].isna().any() or (m["gene_list_sha256"].astype(str).str.len() < 32).any():
        _stop("any external set lacks a source hash")
    if m["source_id"].duplicated().any():
        _stop("duplicate source_id in manifest")
    # Yuan and Roselli must remain separable
    yuan = m[m["source_id"].str.contains("yuan", case=False)]
    roselli = m[m["source_id"].str.contains("roselli", case=False)]
    if len(yuan) != 1 or len(roselli) != 1:
        _stop("Yuan and Roselli cannot be separated as independent included columns")
    if yuan.iloc[0]["display_author_year"] == roselli.iloc[0]["display_author_year"]:
        _stop("Yuan and Roselli share the same display label")
    for role in m["evidence_role"]:
        if role not in ROLE_STYLE:
            _stop(f"unknown evidence_role: {role}")
    return m


def family_order_1d() -> list[int]:
    ord_df = pd.read_csv(ORDER, sep="\t")
    fam = (
        ord_df[["family_display_order", "family_id"]]
        .drop_duplicates()
        .sort_values("family_display_order")
    )
    order = fam["family_id"].astype(int).tolist()
    expected = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]
    if order != expected:
        _stop(f"family membership / order differs from Figure 1D: {order}")
    if len(order) != N_FAMILIES:
        _stop(f"expected {N_FAMILIES} families, got {len(order)}")
    return order


def load_family_sets(universe: set[str]) -> dict[int, set[str]]:
    fle = pd.read_csv(FLE, sep="\t")
    for col in fle.columns:
        if col in FORBIDDEN_READ_COLS:
            _stop(f"forbidden field present in leading-edge table: {col}")
    family_sets = {
        int(fid): {g for g in grp["gene"].astype(str) if g in universe}
        for fid, grp in fle.groupby("family_id")
    }
    if set(family_sets) != set(range(1, N_FAMILIES + 1)):
        _stop("family leading-edge membership does not cover families 1–16")
    return family_sets


def load_external_sets(manifest: pd.DataFrame, universe: set[str]) -> dict[str, set[str]]:
    ext = pd.read_csv(EXT, sep="\t")
    for col in ext.columns:
        if col in FORBIDDEN_READ_COLS:
            _stop(f"forbidden field present in external gene sets: {col}")
    ext_u = ext[ext["in_universe"] == True]  # noqa: E712
    out: dict[str, set[str]] = {}
    for row in manifest.itertuples(index=False):
        sid = str(row.source_id)
        genes_all = sorted(ext.loc[ext["source_id"] == sid, "gene_symbol"].astype(str).unique())
        if not genes_all:
            _stop(f"external set missing from gene-set table: {sid}")
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


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        val = pvals[idx] * n / (i + 1)
        prev = min(prev, val)
        ranked[idx] = min(prev, 1.0)
    return ranked


def load_locked_stats(
    manifest: pd.DataFrame,
    family_sets: dict[int, set[str]],
    ext_sets: dict[str, set[str]],
    labels: pd.DataFrame,
    fam_order: list[int],
) -> pd.DataFrame:
    """Preserve locked hypergeometric OR / p / BH-q; assert against live membership."""
    if not LOCKED_SOURCE.is_file():
        _stop(f"missing locked source table: {LOCKED_SOURCE}")

    locked = pd.read_csv(LOCKED_SOURCE, sep="\t")
    locked["external_set_id"] = locked["external_set_id"].replace(ID_REMAP)

    required = {
        "family_id", "external_set_id", "n_overlap", "overlap_genes",
        "odds_ratio", "hypergeom_p", "fdr_bh", "n_family_le",
        "n_external_in_universe",
    }
    missing = required - set(locked.columns)
    if missing:
        _stop(f"locked source missing columns: {sorted(missing)}")

    set_ids = manifest["source_id"].tolist()
    label_map = labels.set_index("family_id")
    pub = {
        r.source_id: (r.display_author_year, r.display_evidence_label, r.evidence_role)
        for r in manifest.itertuples(index=False)
    }

    rows = []
    for fid in fam_order:
        fgenes = family_sets[fid]
        short = str(label_map.loc[fid, "short_label"])
        axis = str(label_map.loc[fid, "axis_group"])
        for sid in set_ids:
            live_overlap = sorted(fgenes & ext_sets[sid])
            k_live = len(live_overlap)
            hit = locked[
                (locked["family_id"].astype(int) == fid)
                & (locked["external_set_id"] == sid)
            ]
            if hit.empty:
                _stop(f"family–set pair not evaluated in locked source: F{fid} × {sid}")
            if len(hit) != 1:
                _stop(f"duplicate locked rows for F{fid} × {sid}")
            r = hit.iloc[0]
            genes_locked = [
                g.strip() for g in str(r["overlap_genes"]).split(";") if g.strip()
            ] if pd.notna(r["overlap_genes"]) and str(r["overlap_genes"]).strip() else []
            if k_live != int(r["n_overlap"]):
                _stop(
                    f"overlap count disagrees with audited source: "
                    f"F{fid}×{sid} live={k_live} locked={int(r['n_overlap'])}"
                )
            if genes_locked != live_overlap:
                _stop(
                    f"overlap genes disagree with audited source: "
                    f"F{fid}×{sid} live={live_overlap} locked={genes_locked}"
                )
            if int(r["n_external_in_universe"]) != len(ext_sets[sid]):
                _stop(f"external set size drift for {sid}")
            if int(r["n_family_le"]) != len(fgenes):
                _stop(f"family leading-edge size drift for F{fid}")

            line1, line2, role = pub[sid]
            rows.append({
                "family_id": fid,
                "family_display_order": fam_order.index(fid) + 1,
                "family_short_label": short,
                "publication_family_label": f"F{fid} {short}",
                "axis_group": axis,
                "external_set_id": sid,
                "evidence_role": role,
                "publication_label_line1": line1,
                "publication_label_line2": line2,
                "n_family_le": int(r["n_family_le"]),
                "n_external_in_universe": int(r["n_external_in_universe"]),
                "n_overlap": int(r["n_overlap"]),
                "overlap_genes": ";".join(live_overlap),
                "expected_overlap": float(r["expected_overlap"]),
                "odds_ratio": float(r["odds_ratio"]),
                "hypergeom_p": float(r["hypergeom_p"]),
                "fdr_bh": float(r["fdr_bh"]),
                "fdr_significant": bool(float(r["fdr_bh"]) < FDR_ALPHA),
                "evaluated": True,
                "universe_n": UNIVERSE_N,
            })

    df = pd.DataFrame(rows)
    n_sets = len(set_ids)
    if len(df) != N_FAMILIES * n_sets:
        _stop(f"test count != 16 × n_sets ({N_FAMILIES * n_sets}); got {len(df)}")

    # Recompute BH only as a consistency check — do not replace locked q values
    recomputed_q = _bh_fdr(df["hypergeom_p"].to_numpy())
    if not np.allclose(recomputed_q, df["fdr_bh"].to_numpy(), atol=1e-12, rtol=0):
        # Locked tables store BH after full-matrix correction; allow tiny float noise
        if (recomputed_q - df["fdr_bh"].to_numpy()).max() > 1e-10:
            _stop("locked BH-FDR values disagree with recomputed BH from locked p values")

    if int((df["fdr_bh"] < FDR_ALPHA).sum()) != 0:
        _stop("any q < 0.05 appears unexpectedly")
    if df["universe_n"].nunique() != 1 or int(df["universe_n"].iloc[0]) != UNIVERSE_N:
        _stop("universe_n differs across source tables")
    return df


# ---------------------------------------------------------------------------
# Cell text
# ---------------------------------------------------------------------------
def cell_label(n_overlap: int, genes: list[str], evaluated: bool) -> str:
    if not evaluated:
        return "—"
    if n_overlap == 0:
        return ""
    if n_overlap == 1:
        return genes[0]
    if n_overlap == 2:
        return f"{genes[0]}\n{genes[1]}"
    return f"k = {n_overlap}"


def unique_overlap_genes(df: pd.DataFrame) -> list[str]:
    genes: set[str] = set()
    for raw in df.loc[df["n_overlap"] > 0, "overlap_genes"]:
        genes.update(g.strip() for g in str(raw).split(";") if g.strip())
    return sorted(genes)


def panel_strings(
    df: pd.DataFrame, n_sets: int, labels: pd.DataFrame,
) -> dict[str, str]:
    n_tests = N_FAMILIES * n_sets
    n_nonzero = int((df["n_overlap"] > 0).sum())
    n_fam_overlap = int(df.loc[df["n_overlap"] > 0, "family_id"].nunique())
    n_zero_fam = N_FAMILIES - n_fam_overlap
    n_sig = int((df["fdr_bh"] < FDR_ALPHA).sum())
    uniq = unique_overlap_genes(df)
    la_ra = _zero_overlap_la_ra(df, labels)
    if la_ra is None:
        _stop("could not verify LA/RA breakdown of zero-overlap families")
    n_la_zero, n_ra_zero = la_ra

    if n_nonzero != 5:
        _stop(f"n_nonzero_cells != 5 (got {n_nonzero})")
    if len(uniq) != 2 or set(uniq) != EXPECTED_UNIQUE_GENES:
        _stop(f"unique overlap genes mismatch: {uniq}")
    if n_fam_overlap != 3:
        _stop(f"n_families_with_overlap != 3 (got {n_fam_overlap})")
    if n_zero_fam != 13:
        _stop(f"n_zero_overlap_families != 13 (got {n_zero_fam})")
    if n_sig != 0:
        _stop(f"n_FDR_significant != 0 (got {n_sig})")

    gene_clause = " · ".join(uniq)
    return {
        "title_body": "Sparse overlap with external AF gene sets",
        "subtitle": (
            f"{n_nonzero} non-zero cells involving {len(uniq)} unique genes "
            f"across {n_fam_overlap}/16 families"
        ),
        "runtime": (
            f"{n_sig}/{n_tests} enrichments at BH-FDR < {FDR_ALPHA:g} · "
            f"universe N = {UNIVERSE_N:,}"
        ),
        "zero_summary": (
            f"Zero overlap: {n_zero_fam}/16 families—"
            f"{n_la_zero} LA-biased and {n_ra_zero} RA-biased."
        ),
        "gene_summary": (
            f"Unique overlap genes: {gene_clause} · "
            f"FDR-significant enrichments: {n_sig}/{n_tests}."
        ),
        "unique_genes": uniq,
        "n_nonzero": n_nonzero,
        "n_fam_overlap": n_fam_overlap,
        "n_zero_fam": n_zero_fam,
        "n_la_zero": n_la_zero,
        "n_ra_zero": n_ra_zero,
        "n_tests": n_tests,
        "n_sig": n_sig,
    }


def assert_manifest_display(manifest: pd.DataFrame) -> None:
    """Verify displayed author/year/label/size/hash against frozen manifest."""
    for row in manifest.itertuples(index=False):
        ay = str(row.display_author_year).strip()
        if " " not in ay:
            _stop(f"manifest display_author_year lacks year: {ay}")
        author, year_s = ay.rsplit(" ", 1)
        if not year_s.isdigit():
            _stop(f"manifest year not numeric: {ay}")
        year = int(year_s)
        if int(row.publication_year) != year:
            _stop(
                f"manifest publication_year ({row.publication_year}) "
                f"disagrees with display_author_year ({ay})"
            )
        if not str(row.display_panel_label).strip():
            _stop(f"missing panel label for {row.source_id}")
        if int(row.n_genes_normalized) <= 0:
            _stop(f"invalid set size for {row.source_id}")
        if len(str(row.gene_list_sha256)) < 32:
            _stop(f"missing SHA256 for {row.source_id}")

    roselli = manifest[manifest["source_id"].str.contains("roselli", case=False)]
    if len(roselli) != 1:
        _stop("Roselli must appear as exactly one manifest row")
    r = roselli.iloc[0]
    visible = f"{str(r['display_author_year']).rsplit(' ', 1)[0]} {int(r['publication_year'])}"
    expected = f"Roselli {EXPECTED_ROSELLI_YEAR}"
    if visible != expected:
        _stop(
            f"manifest Roselli year/label is '{visible}', "
            f"expected '{expected}' — STOP before render"
        )
    if int(r["publication_year"]) != EXPECTED_ROSELLI_YEAR:
        _stop(
            f"manifest Roselli publication_year={r['publication_year']} "
            f"differs from expected {EXPECTED_ROSELLI_YEAR}"
        )


def _zero_overlap_la_ra(df: pd.DataFrame, labels: pd.DataFrame) -> tuple[int, int] | None:
    """Return (n_LA_zero, n_RA_zero) among families with zero member overlap."""
    axis = dict(zip(labels["family_id"].astype(int), labels["axis_group"].astype(str)))
    overlap_fams = set(df.loc[df["n_overlap"] > 0, "family_id"].astype(int))
    zero_fams = [fid for fid in family_order_1d() if fid not in overlap_fams]
    if len(zero_fams) != 13:
        return None
    n_la = sum(1 for fid in zero_fams if axis.get(fid) == "LA")
    n_ra = sum(1 for fid in zero_fams if axis.get(fid) == "RA")
    if n_la + n_ra != 13:
        return None
    return n_la, n_ra


def _column_geometry(n_cols: int) -> tuple[float, list[float], float, float]:
    x_dir = MARGIN_L
    x_id = x_dir + DIR_W + 0.35
    x_lab = x_id + ID_W
    x_mat = x_lab + LABEL_W + GUTTER_W
    mat_w = WIDTH_MM - MARGIN_R - x_mat
    col_w = (mat_w - COL_GAP * (n_cols - 1)) / n_cols
    col_x = [x_mat + j * (col_w + COL_GAP) for j in range(n_cols)]
    return mat_w, col_x, col_w, x_mat


# ---------------------------------------------------------------------------
# Render — sparse main-text panel
# ---------------------------------------------------------------------------
def _header_font(fontsize: float, weight: str = "normal"):
    from matplotlib.font_manager import FontProperties
    return FontProperties(family=HEADER_FONT, size=fontsize, weight=weight)


def _text_width_mm(ax, fig, text: str, fontproperties) -> float:
    """Measure rendered text width in data/mm units (1 data unit = 1 mm)."""
    t = ax.text(0, 0, text, fontproperties=fontproperties, alpha=0.0)
    fig.canvas.draw()
    bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
    inv = ax.transData.inverted()
    x0, _ = inv.transform((bb.x0, bb.y0))
    x1, _ = inv.transform((bb.x1, bb.y1))
    t.remove()
    return abs(x1 - x0)


def _wrap_family_label(label: str, max_w: float, ax, fig, fp) -> str:
    """Wrap family label to ≤2 lines if wider than max_w."""
    if _text_width_mm(ax, fig, label, fp) <= max_w:
        return label
    # Prefer break at spaces / slashes
    tokens = re.split(r"(\s+|/)", label)
    parts: list[str] = []
    cur = ""
    for tok in tokens:
        trial = f"{cur}{tok}"
        if cur and _text_width_mm(ax, fig, trial, fp) > max_w:
            parts.append(cur.rstrip())
            cur = tok.lstrip() if tok.startswith(" ") else tok
        else:
            cur = trial
    if cur:
        parts.append(cur.rstrip())
    if len(parts) <= 1:
        return label
    # Keep at most two lines
    if len(parts) == 2:
        return f"{parts[0]}\n{parts[1]}"
    return f"{parts[0]}\n{' '.join(parts[1:]).replace(' /', '/').strip()}"


def _header_plan(
    ax, fig, manifest: pd.DataFrame, col_w: float,
) -> tuple[bool, list[bool], list[float]]:
    """Decide wrap per column; return (any_wrap, wraps, half_widths)."""
    hdr = _header_font(FS_HEADER, "normal")
    hdr_bold = _header_font(FS_HEADER, "semibold")
    hdr_n = _header_font(FS_HEADER - 0.15, "normal")
    max_w = col_w - 2 * HEADER_PAD_MM
    wraps: list[bool] = []
    halves: list[float] = []
    for row in manifest.itertuples(index=False):
        author = str(row.display_author_year).rsplit(" ", 1)[0]
        year = int(row.publication_year)
        line1 = f"{author} {year}"
        panel_label = str(row.display_panel_label)
        n_genes = int(row.n_genes_normalized)
        line2 = f"{panel_label} · n = {n_genes}"
        wrap = _text_width_mm(ax, fig, line2, hdr) > max_w
        wraps.append(wrap)
        w1 = _text_width_mm(ax, fig, line1, hdr_bold)
        if wrap:
            w2 = _text_width_mm(ax, fig, panel_label, hdr)
            w3 = _text_width_mm(ax, fig, f"· n = {n_genes}", hdr_n)
            halves.append(max(w1, w2, w3) / 2)
        else:
            halves.append(max(w1, _text_width_mm(ax, fig, line2, hdr)) / 2)
    return any(wraps), wraps, halves


def _draw_set_headers(
    ax, fig, manifest: pd.DataFrame, col_x: list[float], col_w: float,
    body_top: float, head_h: float, wraps: list[bool],
) -> list[tuple[float, float]]:
    """
    Draw set headers from the frozen manifest.

    Specified content (two lines):
      Author Year
      PanelLabel · n = N

    The second content line wraps visually when wider than the column.
    """
    hdr = _header_font(FS_HEADER, "normal")
    hdr_bold = _header_font(FS_HEADER, "semibold")
    hdr_n = _header_font(FS_HEADER - 0.15, "normal")
    spans: list[tuple[float, float]] = []

    for j, row in enumerate(manifest.itertuples(index=False)):
        role = str(row.evidence_role)
        style = ROLE_STYLE[role]
        cx = col_x[j]
        mid = cx + col_w / 2
        ax.plot(
            [cx + 0.4, cx + col_w - 0.4],
            [body_top + head_h - 0.35, body_top + head_h - 0.35],
            color=style["rule"], lw=1.25, solid_capstyle="round", zorder=3,
        )
        author = str(row.display_author_year).rsplit(" ", 1)[0]
        year = int(row.publication_year)
        line1 = f"{author} {year}"
        panel_label = str(row.display_panel_label)
        n_genes = int(row.n_genes_normalized)
        line2 = f"{panel_label} · n = {n_genes}"

        if wraps[j]:
            y1 = body_top + head_h - 2.15
            y2 = body_top + head_h - 4.35
            y3 = body_top + head_h - 6.45
            t1 = ax.text(
                mid, y1, line1, ha="center", va="center",
                fontproperties=hdr_bold, color=FG_COLOR, zorder=3, clip_on=True,
            )
            t2 = ax.text(
                mid, y2, panel_label, ha="center", va="center",
                fontproperties=hdr, color=MUTED_COLOR, zorder=3, clip_on=True,
            )
            t3 = ax.text(
                mid, y3, f"· n = {n_genes}", ha="center", va="center",
                fontproperties=hdr_n, color=MUTED_COLOR, zorder=3, clip_on=True,
            )
            fig.canvas.draw()
            xs: list[float] = []
            for t in (t1, t2, t3):
                bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
                inv = ax.transData.inverted()
                x0, _ = inv.transform((bb.x0, bb.y0))
                x1, _ = inv.transform((bb.x1, bb.y1))
                xs.extend([x0, x1])
            spans.append((min(xs), max(xs)))
        else:
            y1 = body_top + head_h - 2.35
            y2 = body_top + head_h - 4.85
            t1 = ax.text(
                mid, y1, line1, ha="center", va="center",
                fontproperties=hdr_bold, color=FG_COLOR, zorder=3, clip_on=True,
            )
            t2 = ax.text(
                mid, y2, line2, ha="center", va="center",
                fontproperties=hdr, color=MUTED_COLOR, zorder=3, clip_on=True,
            )
            fig.canvas.draw()
            xs = []
            for t in (t1, t2):
                bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
                inv = ax.transData.inverted()
                x0, _ = inv.transform((bb.x0, bb.y0))
                x1, _ = inv.transform((bb.x1, bb.y1))
                xs.extend([x0, x1])
            spans.append((min(xs), max(xs)))

    return spans


def render(df: pd.DataFrame, manifest: pd.DataFrame, strings: dict[str, str]) -> plt.Figure:
    set_nature_style()
    fig = plt.figure(
        figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4),
        facecolor="white", dpi=300,
    )
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, WIDTH_MM)
    ax.set_ylim(0, HEIGHT_MM)
    ax.set_axis_off()

    labels = pd.read_csv(LABELS, sep="\t")
    label_map = {
        int(r.family_id): r for r in labels.itertuples(index=False)
    }
    set_ids = manifest["source_id"].tolist()
    n_cols = len(set_ids)
    display_ids = list(DISPLAY_FAMILY_IDS)
    n_rows = len(display_ids)
    mat_w, col_x, col_w, x_mat = _column_geometry(n_cols)

    x0 = MARGIN_L
    x_dir = MARGIN_L
    x_id = x_dir + DIR_W + 0.35
    x_lab = x_id + ID_W

    # Compact title block — shared C/D baselines
    y = CD_LETTER_Y
    ax.text(
        x0, y, "C",
        fontsize=FS_LETTER, fontweight="bold",
        ha="left", va="top", color=FG_COLOR, clip_on=True,
    )
    ax.text(
        x0 + 3.4, y, strings["title_body"],
        fontsize=FS_TITLE, fontweight="bold",
        ha="left", va="top", color=FG_COLOR, clip_on=True,
    )
    y -= CD_SUB_DY
    ax.text(
        x0, y, strings["subtitle"],
        fontsize=FS_SUB, ha="left", va="top", color="#3A3A3A", clip_on=True,
    )
    y -= CD_BOUND_DY
    ax.text(
        x0, y, strings["runtime"],
        fontsize=FS_RUNTIME, ha="left", va="top", color=MUTED_COLOR, clip_on=True,
    )

    any_wrap, wraps, _halves = _header_plan(ax, fig, manifest, col_w)
    head_h = 8.6 if any_wrap else 6.8
    # Role-rule y = body_top + head_h − 0.35 must sit below runtime text
    body_top = y - TITLE_TO_HEADER_GAP_MM - head_h
    row_h = ROW_H_MM
    matrix_bot = body_top - n_rows * row_h
    # Packed summary + footer height under the matrix
    summary_block_h = 5.0
    footer_block_h = 5.2
    bottom_margin = 1.0
    content_bot = bottom_margin + footer_block_h + summary_block_h
    if matrix_bot < content_bot + 0.6:
        _stop(
            f"matrix collides with summary/footer "
            f"(matrix_bot={matrix_bot:.1f}, content_bot={content_bot:.1f}); "
            f"increase HEIGHT_MM or reduce ROW_H_MM/head_h"
        )

    lookup = {
        (int(r.family_id), r.external_set_id): r
        for r in df.itertuples(index=False)
    }

    header_spans = _draw_set_headers(
        ax, fig, manifest, col_x, col_w, body_top, head_h, wraps,
    )

    # Light guides only — no heavy outer grid; do not stretch empty cells
    for i in range(n_rows + 1):
        yy = body_top - i * row_h
        lw = 0.55 if i == 2 else 0.28
        color = AXIS_SEP if i == 2 else ROW_RULE
        ax.plot(
            [x_mat, x_mat + mat_w], [yy, yy],
            color=color, lw=lw, zorder=1,
        )
    for j in range(n_cols):
        ax.plot(
            [col_x[j], col_x[j]], [matrix_bot, body_top],
            color=COL_RULE, lw=0.28, zorder=1,
        )
    ax.plot(
        [col_x[-1] + col_w, col_x[-1] + col_w], [matrix_bot, body_top],
        color=COL_RULE, lw=0.28, zorder=1,
    )

    lab_fp = _header_font(FS_LABEL, "normal")
    pill_boxes: list[tuple[float, float, float, float, float, float]] = []

    # Rows: only families with ≥1 overlap (reduced row / pill height)
    for i, fid in enumerate(display_ids):
        y_top = body_top - i * row_h
        y_bot = y_top - row_h
        y_mid = (y_top + y_bot) / 2
        axis = str(label_map[fid].axis_group)
        dir_col = LA_DEEP if axis == "LA" else RA_DEEP

        tag_h = min(2.15, row_h * 0.52)
        ax.add_patch(FancyBboxPatch(
            (x_dir, y_mid - tag_h / 2), DIR_W - 0.25, tag_h,
            boxstyle="round,pad=0.015,rounding_size=0.12",
            facecolor=dir_col, edgecolor="none", zorder=2,
        ))
        ax.text(
            x_dir + (DIR_W - 0.25) / 2, y_mid, axis,
            ha="center", va="center",
            fontproperties=_header_font(FS_HEADER - 0.3, "bold"),
            color="white", zorder=3, clip_on=True,
        )

        ax.text(
            x_id, y_mid, f"F{fid}",
            ha="left", va="center",
            fontproperties=_header_font(FS_LABEL, "semibold"),
            color=ID_COLOR, zorder=3, clip_on=True,
        )
        fam_txt = _wrap_family_label(
            PUB_FAMILY_LABEL[fid], LABEL_W - 0.4, ax, fig, lab_fp,
        )
        ax.text(
            x_lab, y_mid, fam_txt,
            ha="left", va="center",
            fontproperties=lab_fp, color=FG_COLOR, zorder=3, clip_on=True,
            linespacing=0.95,
        )

        for j, sid in enumerate(set_ids):
            cell = lookup[(fid, sid)]
            evaluated = bool(cell.evaluated)
            k = int(cell.n_overlap)
            genes = [
                g.strip() for g in str(cell.overlap_genes).split(";") if g.strip()
            ] if k else []
            txt = cell_label(k, genes, evaluated)
            cx = col_x[j]
            if not evaluated:
                ax.text(
                    cx + col_w / 2, y_mid, "—",
                    ha="center", va="center",
                    fontsize=FS_LABEL, color=MUTED_COLOR, zorder=3,
                )
                continue
            if k == 0:
                continue
            role = str(cell.evidence_role)
            style = ROLE_STYLE[role]
            pill_w = min(col_w - 2.4, 7.2)
            pill_h = min(row_h * 0.48, 3.15)  # ~10–15% shorter gene pills
            pill_x = cx + (col_w - pill_w) / 2
            pill_y = y_mid - pill_h / 2
            ax.add_patch(FancyBboxPatch(
                (pill_x, pill_y), pill_w, pill_h,
                boxstyle="round,pad=0.015,rounding_size=0.14",
                facecolor=style["fill"], edgecolor=style["edge"],
                linewidth=0.65, zorder=2,
            ))
            ax.text(
                cx + col_w / 2, y_mid, txt,
                ha="center", va="center",
                fontproperties=_header_font(FS_LABEL, "semibold"),
                color=FG_COLOR, zorder=3, clip_on=True,
            )
            pill_boxes.append((pill_x, pill_y, pill_w, pill_h, y_bot, y_top))

    # Compact bottom summary + footer packed under the matrix (no large box)
    sum_y = matrix_bot - 1.25
    ax.text(
        x0, sum_y, strings["zero_summary"],
        fontsize=FS_SUMMARY, ha="left", va="top", color=FG_COLOR, clip_on=True,
    )
    ax.text(
        x0, sum_y - 2.3, strings["gene_summary"],
        fontsize=FS_SUMMARY, ha="left", va="top", color="#3A3A3A", clip_on=True,
    )

    foot_items = [
        "Labels denote overlapping frozen family genes; blanks denote tested zero overlap.",
        "Member overlap is a traceability cue, not family-level enrichment or validation.",
    ]
    # Shared C/D footer baseline (second line below)
    foot_y = CD_FOOTER_BASELINE
    for line in foot_items:
        ax.text(
            x0, foot_y, line,
            fontsize=FS_FOOT, ha="left", va="top",
            color=MUTED_COLOR, clip_on=True,
        )
        foot_y -= 2.0

    # Layout geometry for QC (mm)
    rule_y = body_top + head_h - 0.35
    fig._fig6c_layout = {  # type: ignore[attr-defined]
        "col_x": col_x,
        "col_w": col_w,
        "header_spans": header_spans,
        "pill_boxes": pill_boxes,
        "label_w": LABEL_W,
        "row_h": row_h,
        "head_h": head_h,
        "body_top": body_top,
        "matrix_bot": matrix_bot,
        "runtime_y": y,
        "rule_y": rule_y,
        "sum_y": sum_y,
    }
    return fig


def render_supplement_full(
    df: pd.DataFrame, manifest: pd.DataFrame,
) -> plt.Figure:
    """Optional 16-row full matrix for the supplement (all tested pairs)."""
    from matplotlib.font_manager import FontProperties

    w_mm, h_mm = 120.0, 110.0
    set_nature_style()
    fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4), facecolor="white", dpi=300)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, w_mm)
    ax.set_ylim(0, h_mm)
    ax.set_axis_off()

    fam_order = family_order_1d()
    set_ids = manifest["source_id"].tolist()
    labels = pd.read_csv(LABELS, sep="\t")
    label_map = {int(r.family_id): r for r in labels.itertuples(index=False)}
    lookup = {
        (int(r.family_id), r.external_set_id): r
        for r in df.itertuples(index=False)
    }

    x0, x_id, x_lab, x_mat = 2.0, 4.0, 10.0, 42.0
    mat_w = w_mm - x_mat - 2.0
    col_w = mat_w / len(set_ids)
    col_x = [x_mat + j * col_w for j in range(len(set_ids))]

    ax.text(
        x0, h_mm - 3.0,
        "Supplementary — full 16 × 4 family × external AF gene-set overlap",
        fontsize=8.0, fontweight="bold", ha="left", va="top", color=FG_COLOR,
    )
    ax.text(
        x0, h_mm - 6.0,
        "All 64 tested pairs; main-text Figure 6C shows only the three families with nonzero overlap.",
        fontsize=6.5, ha="left", va="top", color=MUTED_COLOR,
    )

    body_top = h_mm - 16.0
    tracks_bot = 6.0
    row_h = (body_top - tracks_bot) / N_FAMILIES
    uni = FontProperties(family="DejaVu Sans", size=6.0)

    for j, row in enumerate(manifest.itertuples(index=False)):
        mid = col_x[j] + col_w / 2
        ax.plot(
            [col_x[j] + 0.5, col_x[j] + col_w - 0.5],
            [body_top + 8.5, body_top + 8.5],
            color=ROLE_STYLE[str(row.evidence_role)]["rule"], lw=1.2, zorder=3,
        )
        ax.text(mid, body_top + 6.2, str(row.display_author_year),
                ha="center", va="center", fontsize=6.0, fontweight="bold")
        ax.text(mid, body_top + 3.8, str(row.display_panel_label),
                ha="center", va="center", fontsize=5.8, color=MUTED_COLOR)
        ax.text(mid, body_top + 1.6, f"n = {int(row.n_genes_normalized)}",
                ha="center", va="center", fontsize=5.8, color=MUTED_COLOR)

    for i, fid in enumerate(fam_order):
        y_mid = body_top - (i + 0.5) * row_h
        ax.text(x_id, y_mid, f"F{fid}", ha="left", va="center",
                fontsize=6.0, fontweight="bold", color=ID_COLOR)
        ax.text(x_lab, y_mid, PUB_FAMILY_LABEL[fid], ha="left", va="center",
                fontproperties=uni, color=FG_COLOR)
        ax.plot([x_mat, x_mat + mat_w], [body_top - i * row_h, body_top - i * row_h],
                color=ROW_RULE, lw=0.3, zorder=1)
        for j, sid in enumerate(set_ids):
            cell = lookup[(fid, sid)]
            k = int(cell.n_overlap)
            if k == 0:
                continue
            genes = [g for g in str(cell.overlap_genes).split(";") if g]
            style = ROLE_STYLE[str(cell.evidence_role)]
            cx = col_x[j] + col_w / 2
            ax.add_patch(FancyBboxPatch(
                (cx - 3.2, y_mid - 1.3), 6.4, 2.6,
                boxstyle="round,pad=0.02,rounding_size=0.12",
                facecolor=style["fill"], edgecolor=style["edge"], linewidth=0.6,
            ))
            ax.text(cx, y_mid, genes[0] if k == 1 else f"k={k}",
                    ha="center", va="center", fontsize=6.0, fontweight="bold")

    ax.text(
        x0, 3.0,
        "Blank cells = tested zero member overlap. No family–set enrichment survived BH-FDR.",
        fontsize=5.5, ha="left", va="center", color=MUTED_COLOR,
    )
    fig._supp_size_mm = (w_mm, h_mm)  # type: ignore[attr-defined]
    return fig


# ---------------------------------------------------------------------------
# QC / export helpers
# ---------------------------------------------------------------------------
def _pdf_page_mm(pdf: Path) -> tuple[float, float]:
    raw = pdf.read_bytes()
    m = re.search(rb"/MediaBox\s*\[\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*\]", raw)
    if not m:
        _stop("could not parse PDF MediaBox")
    x0, y0, x1, y1 = map(float, m.groups())
    return (x1 - x0) * 25.4 / 72.0, (y1 - y0) * 25.4 / 72.0


def _svg_mm(svg: Path) -> tuple[float, float]:
    root = ET.parse(svg).getroot()
    w = root.attrib.get("width", "")
    h = root.attrib.get("height", "")

    def _to_mm(val: str) -> float:
        if val.endswith("pt"):
            return float(val[:-2]) * 25.4 / 72.0
        if val.endswith("in"):
            return float(val[:-2]) * 25.4
        if val.endswith("mm"):
            return float(val[:-2])
        if val.endswith("px"):
            return float(val[:-2]) * 25.4 / 72.0
        return float(val) * 25.4 / 72.0

    return _to_mm(w), _to_mm(h)


def _png_mm(png: Path, dpi: int = 600) -> tuple[float, float, int, int]:
    data = png.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    w, h = struct.unpack(">II", data[16:24])
    return w * 25.4 / dpi, h * 25.4 / dpi, w, h


def _min_svg_font(svg: Path) -> float:
    """Return minimum font-size in pt (SVG may use px or bare numbers)."""
    text = svg.read_text(encoding="utf-8")
    sizes = [
        float(x) for x in re.findall(r"font-size:\s*([0-9.]+)px", text)
    ]
    sizes += [
        float(x) for x in re.findall(r"font-size:\s*([0-9.]+)pt", text)
    ]
    sizes += [float(x) for x in re.findall(r'font-size="([0-9.]+)"', text)]
    # Matplotlib SVG often emits style font-size in px ≈ pt at 72 dpi export
    return min(sizes) if sizes else 99.0


def _assert_layout_qc(layout: dict) -> None:
    """No header collisions; pills stay inside their row cells."""
    col_x = layout["col_x"]
    col_w = layout["col_w"]
    spans = layout["header_spans"]
    assert len(spans) == len(col_x) == 4
    for j, (x0, x1) in enumerate(spans):
        assert x0 >= col_x[j] - 0.15, (j, x0, col_x[j])
        assert x1 <= col_x[j] + col_w + 0.15, (j, x1, col_x[j] + col_w)
        if j > 0:
            prev = spans[j - 1]
            assert x0 >= prev[1] + 0.2, (
                f"header overlap cols {j-1}/{j}: {prev} vs {(x0, x1)}"
            )
    for pill_x, pill_y, pill_w, pill_h, y_bot, y_top in layout["pill_boxes"]:
        assert pill_y >= y_bot - 0.05
        assert pill_y + pill_h <= y_top + 0.05
        assert pill_w > 0 and pill_h > 0
    # Role rules must sit below the runtime summary line
    assert layout["rule_y"] <= layout["runtime_y"] - 2.4, (
        layout["rule_y"], layout["runtime_y"],
    )
    assert layout["sum_y"] <= layout["matrix_bot"] - 0.8
    assert 72.0 <= HEIGHT_MM <= 76.5
    assert abs(layout["row_h"] - ROW_H_MM) < 0.05


def run_qc(
    df: pd.DataFrame,
    manifest: pd.DataFrame,
    strings: dict[str, str],
    layout: dict | None = None,
) -> dict[str, object]:
    n_sets = len(manifest)
    n_tests = N_FAMILIES * n_sets
    n_sig = int((df["fdr_bh"] < FDR_ALPHA).sum())
    min_q = float(df["fdr_bh"].min())

    assert len(df["family_id"].unique()) == N_FAMILIES
    assert len(df) == n_tests == 64
    assert n_sets == 4
    assert n_sig == 0
    assert min_q >= FDR_ALPHA - 1e-15
    assert int(df["universe_n"].iloc[0]) == UNIVERSE_N
    assert set(manifest["source_id"]) == set(df["external_set_id"])
    assert (df["evaluated"] == True).all()  # noqa: E712 — no not-evaluated as zero

    for r in df.itertuples(index=False):
        if not r.evaluated:
            continue
        if int(r.n_overlap) == 0:
            assert (not r.overlap_genes) or str(r.overlap_genes).strip() == ""
        else:
            genes = [g for g in str(r.overlap_genes).split(";") if g]
            assert len(genes) == int(r.n_overlap)

    svg_text = OUT_SVG.read_text(encoding="utf-8")
    visible = re.sub(r"<[^>]+>", " ", svg_text)
    for bad in (
        "NCP_v3", "choi_mendelian_rare", "yuan2025_cross_population",
        "roselli2018_loci", "selewa2023_regulatory", "PriorityScore",
        "Path_A", "G-tier", "−log10 FDR", "-log10 FDR", "colorbar",
        "COMMON-RISK / PRIORITIZED", "REGULATORY", "RARE / FUNCTIONAL",
    ):
        assert bad not in visible, bad
    assert "★" not in visible

    pdf_w, pdf_h = _pdf_page_mm(OUT_PDF)
    svg_w, svg_h = _svg_mm(OUT_SVG)
    png_w, png_h, png_px_w, png_px_h = _png_mm(OUT_PNG, dpi=600)
    assert abs(pdf_w - WIDTH_MM) < 1.5, pdf_w
    assert 72.0 <= pdf_h <= 76.5, pdf_h
    assert abs(pdf_h - HEIGHT_MM) < 1.5, pdf_h
    assert abs(svg_w - pdf_w) < 1.5 and abs(svg_h - pdf_h) < 1.5
    assert abs(png_w - pdf_w) < 1.5 and abs(png_h - pdf_h) < 1.5

    min_font = _min_svg_font(OUT_SVG)
    assert min_font >= FS_FOOT - 0.15, min_font

    n_nonzero = int((df["n_overlap"] > 0).sum())
    overlap_fams = sorted(df.loc[df["n_overlap"] > 0, "family_id"].astype(int).unique())
    n_fam_overlap = len(overlap_fams)
    n_zero_fam = N_FAMILIES - n_fam_overlap
    uniq = unique_overlap_genes(df)
    assert n_nonzero == 5
    assert len(uniq) == 2 and set(uniq) == EXPECTED_UNIQUE_GENES
    assert set(overlap_fams) == set(DISPLAY_FAMILY_IDS) == {3, 13, 9}
    assert n_zero_fam == 13
    assert DISPLAY_FAMILY_IDS == [3, 13, 9]

    got = {}
    for r in df.itertuples(index=False):
        if int(r.n_overlap) > 0:
            genes = tuple(
                g.strip() for g in str(r.overlap_genes).split(";") if g.strip()
            )
            got[(int(r.family_id), r.external_set_id)] = genes
    assert got == EXPECTED_NONZERO, got

    choi = df[df["external_set_id"] == "choi_mendelian_rare"]
    assert (choi["n_overlap"] == 0).all()
    assert len(choi) == N_FAMILIES

    # All 13 omitted families have overlap_n = 0 across all four sets
    omitted = [fid for fid in family_order_1d() if fid not in set(DISPLAY_FAMILY_IDS)]
    assert len(omitted) == 13
    for fid in omitted:
        sub = df[df["family_id"].astype(int) == fid]
        assert len(sub) == 4
        assert (sub["n_overlap"] == 0).all(), fid

    for fid in DISPLAY_FAMILY_IDS:
        # Wrapped labels may insert newlines; check token presence
        for tok in re.split(r"[\s/]+", PUB_FAMILY_LABEL[fid]):
            if tok:
                assert tok in visible, tok
    for fid in (1, 2, 4, 16):
        assert PUB_FAMILY_LABEL[fid] not in visible, PUB_FAMILY_LABEL[fid]

    assert "Sparse overlap with external AF gene sets" in visible
    assert "2 unique genes" in visible
    assert "0/64 enrichments" in visible
    assert "universe N = 14,635" in visible
    assert "Zero overlap: 13/16 families" in visible
    assert "9 LA-biased and 4 RA-biased" in visible
    assert "Unique overlap genes: CAV1 · HEY2" in visible
    assert "CAV1" in visible and "HEY2" in visible
    assert "Roselli 2018" in visible
    assert "All 64 family–set pairs tested" not in visible
    assert "source table" not in visible

    # Manifest-derived headers visible; SHA256 recorded in QC/manifest
    for row in manifest.itertuples(index=False):
        author = str(row.display_author_year).rsplit(" ", 1)[0]
        year = int(row.publication_year)
        assert f"{author} {year}" in visible
        assert str(row.display_panel_label) in visible
        assert f"n = {int(row.n_genes_normalized)}" in visible
        assert len(str(row.gene_list_sha256)) == 64

    for phrase in (
        "Cross-population GWAS", "Gene prioritization",
        "Regulatory mapping", "Rare coding and SV",
    ):
        assert phrase not in visible, phrase

    if layout is not None:
        _assert_layout_qc(layout)

    assert OUT_SOURCE.is_file()
    src = pd.read_csv(OUT_SOURCE, sep="\t")
    assert len(src) == 64

    nonzero = df[df["n_overlap"] > 0].copy()
    return {
        "n_sets": n_sets,
        "n_tests": n_tests,
        "n_sig": n_sig,
        "n_nonzero": n_nonzero,
        "n_unique_genes": len(uniq),
        "unique_genes": uniq,
        "n_fam_overlap": n_fam_overlap,
        "n_zero_fam": n_zero_fam,
        "min_q": min_q,
        "strings": strings,
        "nonzero": nonzero,
        "pdf_mm": (round(pdf_w, 2), round(pdf_h, 2)),
        "svg_mm": (round(svg_w, 2), round(svg_h, 2)),
        "png_mm": (round(png_w, 2), round(png_h, 2)),
        "png_px": (png_px_w, png_px_h),
        "min_font_px": min_font,
        "source_sha256": _sha256(OUT_SOURCE),
        "pdf_sha256": _sha256(OUT_PDF),
        "svg_sha256": _sha256(OUT_SVG),
        "png_sha256": _sha256(OUT_PNG),
        "manifest_sha256": _sha256(MANIFEST_EXT),
        "layout_qc": "PASS",
        "scientific_qc": "PASS",
        "col_w_mm": round(float(layout["col_w"]), 2) if layout else None,
        "row_h_mm": round(float(layout["row_h"]), 2) if layout else None,
    }


def write_outputs(df: pd.DataFrame, manifest: pd.DataFrame, strings: dict[str, str], qc: dict) -> None:
    """Update figure manifest + QC only; leave frozen source tables untouched."""
    man_rows = [
        {
            "artifact": "publication_pdf",
            "path": str(OUT_PDF.relative_to(ROOT)),
            "sha256": qc["pdf_sha256"],
            "geometry_mm": f"{qc['pdf_mm'][0]}x{qc['pdf_mm'][1]}",
            "notes": "sparse main-text Figure 6C (3 overlap families)",
        },
        {
            "artifact": "publication_svg",
            "path": str(OUT_SVG.relative_to(ROOT)),
            "sha256": qc["svg_sha256"],
            "geometry_mm": f"{qc['svg_mm'][0]}x{qc['svg_mm'][1]}",
            "notes": "sparse main-text Figure 6C",
        },
        {
            "artifact": "publication_png_600dpi",
            "path": str(OUT_PNG.relative_to(ROOT)),
            "sha256": qc["png_sha256"],
            "geometry_mm": f"{qc['png_mm'][0]}x{qc['png_mm'][1]}",
            "notes": "600 dpi",
        },
        {
            "artifact": "final_width_preview",
            "path": str(OUT_PREVIEW.relative_to(ROOT)),
            "sha256": _sha256(OUT_PREVIEW),
            "geometry_mm": "",
            "notes": "150 dpi print preview",
        },
        {
            "artifact": "grayscale_preview",
            "path": str(OUT_GRAY.relative_to(ROOT)),
            "sha256": _sha256(OUT_GRAY),
            "geometry_mm": "",
            "notes": "grayscale interpretability check",
        },
        {
            "artifact": "supplement_full_matrix_pdf",
            "path": str(OUT_SUPP_PDF.relative_to(ROOT)),
            "sha256": _sha256(OUT_SUPP_PDF) if OUT_SUPP_PDF.is_file() else "",
            "geometry_mm": "",
            "notes": "optional full 16×4 matrix",
        },
        {
            "artifact": "source_table",
            "path": str(OUT_SOURCE.relative_to(ROOT)),
            "sha256": qc["source_sha256"],
            "geometry_mm": "",
            "notes": "frozen 64-test hypergeometric statistics",
        },
        {
            "artifact": "external_set_manifest",
            "path": str(MANIFEST_EXT.relative_to(ROOT)),
            "sha256": qc["manifest_sha256"],
            "geometry_mm": "",
            "notes": "frozen included_in_6C sets",
        },
    ]
    pd.DataFrame(man_rows).to_csv(OUT_MANIFEST, sep="\t", index=False)

    # Order nonzero lines in display order
    order_key = {fid: i for i, fid in enumerate(DISPLAY_FAMILY_IDS)}
    nz = qc["nonzero"].copy()
    nz["_ord"] = nz["family_id"].map(order_key)
    nz = nz.sort_values(["_ord", "external_set_id"])
    nonzero_lines = [
        f"- F{int(r.family_id)} × {r.publication_label_line1}: {r.overlap_genes}"
        for r in nz.itertuples(index=False)
    ]

    set_lines = [
        f"- {r.display_author_year} / {r.display_panel_label} / "
        f"n={int(r.n_genes_normalized)} [{r.evidence_role}]; "
        f"SHA256=`{r.gene_list_sha256}`"
        for r in manifest.itertuples(index=False)
    ]

    OUT_QC.write_text(
        "\n".join(
            [
                "# Figure 6C sparse external AF gene-set overlap QC",
                "",
                f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
                "",
                "Status: **PASS**",
                "",
                "## Scope",
                "- Main-text panel shows only families with ≥1 nonzero member overlap.",
                "- Full 16 × 4 results retained in frozen source table.",
                "- No change to OR / p / BH-q / universe / gene sets.",
                "- No formal Figure 6 A–D composite routing update.",
                "",
                "## Panel copy",
                f"- title: `C  {strings['title_body']}`",
                f"- subtitle: `{strings['subtitle']}`",
                f"- runtime: `{strings['runtime']}`",
                f"- zero summary: `{strings['zero_summary']}`",
                f"- gene summary: `{strings['gene_summary']}`",
                "",
                "## Displayed rows (Figure 1D order among nonzero families)",
                "- F3 Adipo-lipid axis (LA)",
                "- F13 Amino-acid metabolism (LA)",
                "- F9 HS-GAG / bone / GPCR (RA)",
                "",
                "## Displayed external sets (manifest-derived)",
                *set_lines,
                "",
                "## Exact visible nonzero cells",
                *nonzero_lines,
                "",
                "## Counts",
                f"- families tested = {N_FAMILIES}",
                f"- external sets = {qc['n_sets']}",
                f"- tests = {qc['n_tests']}",
                f"- non-zero cells = {qc['n_nonzero']}",
                f"- unique overlap genes = {qc['n_unique_genes']} "
                f"({', '.join(qc['unique_genes'])})",
                f"- families with overlap = {qc['n_fam_overlap']}",
                f"- families with zero overlap = {qc['n_zero_fam']}",
                f"- FDR-significant (BH < 0.05) = {qc['n_sig']}",
                f"- minimum q = {qc['min_q']}",
                "",
                "## Layout QC",
                f"- canvas {WIDTH_MM:.0f} × {HEIGHT_MM:.0f} mm (target 84 × 64–70)",
                f"- PDF: {qc['pdf_mm'][0]} × {qc['pdf_mm'][1]} mm",
                f"- SVG: {qc['svg_mm'][0]} × {qc['svg_mm'][1]} mm",
                f"- PNG: {qc['png_mm'][0]} × {qc['png_mm'][1]} mm "
                f"({qc['png_px'][0]} × {qc['png_px'][1]} px @ 600 dpi)",
                f"- column width: {qc.get('col_w_mm')} mm; row height: {qc.get('row_h_mm')} mm",
                f"- minimum SVG font-size: {qc['min_font_px']:.2f} px",
                "- no set-header overlap; gene pills inside cell bounds",
                f"- layout_qc = {qc['layout_qc']}; scientific_qc = {qc['scientific_qc']}",
                "",
                "## Provenance (manifest-matched)",
                *set_lines,
                "",
                "## SHA256",
                f"- source (frozen): `{qc['source_sha256']}`",
                f"- external-set manifest: `{qc['manifest_sha256']}`",
                f"- pdf: `{qc['pdf_sha256']}`",
                f"- svg: `{qc['svg_sha256']}`",
                f"- png: `{qc['png_sha256']}`",
                "",
                "## Legend note (not on panel)",
                "- Repeated CAV1 labels reflect shared membership across frozen "
                "family sets and external gene sets; matrix cells are not "
                "independent genetic observations.",
                "- Full 16 × 4 results are available in the figure source table.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    manifest = load_manifest()
    assert_manifest_display(manifest)
    fam_order = family_order_1d()
    labels = pd.read_csv(LABELS, sep="\t")
    universe_df = pd.read_csv(UNI, sep="\t")
    universe = set(universe_df["gene"].astype(str))
    if len(universe) != UNIVERSE_N:
        _stop(f"universe_n differs: live={len(universe)} expected={UNIVERSE_N}")

    family_sets = load_family_sets(universe)
    ext_sets = load_external_sets(manifest, universe)
    df = load_locked_stats(manifest, family_sets, ext_sets, labels, fam_order)

    strings = panel_strings(df, len(manifest), labels)

    fig = render(df, manifest, strings)
    layout = dict(fig._fig6c_layout)  # type: ignore[attr-defined]
    save_figure_bundle(
        fig, OUT_PDF, dpi=600,
        bbox_inches=_fixed_bbox(), pad_inches=0,
    )
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")
    print(f"wrote {OUT_PNG.relative_to(ROOT)}")

    fig = render(df, manifest, strings)
    layout = dict(fig._fig6c_layout)  # type: ignore[attr-defined]
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

    # Optional supplementary full 16×4 matrix
    fig_s = render_supplement_full(df, manifest)
    w_mm, h_mm = fig_s._supp_size_mm  # type: ignore[attr-defined]
    bb_s = _fixed_bbox(w_mm, h_mm)
    save_figure_bundle(fig_s, OUT_SUPP_PDF, dpi=600, bbox_inches=bb_s, pad_inches=0)
    fig_s = render_supplement_full(df, manifest)
    fig_s.savefig(
        OUT_SUPP_SVG, format="svg", bbox_inches=bb_s, pad_inches=0,
        facecolor="white", edgecolor="none",
    )
    fig_s.savefig(
        OUT_SUPP_PNG, format="png", dpi=300, bbox_inches=bb_s, pad_inches=0,
        facecolor="white", edgecolor="none",
    )
    plt.close(fig_s)
    print(f"wrote {OUT_SUPP_PDF.relative_to(ROOT)}")

    print(f"wrote {OUT_SVG.relative_to(ROOT)}")
    print(f"wrote {OUT_PREVIEW.relative_to(ROOT)}")

    qc = run_qc(df, manifest, strings, layout=layout)
    write_outputs(df, manifest, strings, qc)

    print(f"wrote {OUT_MANIFEST.relative_to(ROOT)}")
    print(f"wrote {OUT_QC.relative_to(ROOT)}")
    print(f"subtitle: {strings['subtitle']}")
    print(f"canvas_mm: {WIDTH_MM:.0f} x {HEIGHT_MM:.0f}")
    print("QC: PASS")
    print(f"SOURCE_SHA256={qc['source_sha256']}")
    print(f"PDF_SHA256={qc['pdf_sha256']}")
    print(f"SVG_SHA256={qc['svg_sha256']}")
    print(f"PNG_SHA256={qc['png_sha256']}")
    print(f"MIN_Q={qc['min_q']}")
    print(f"N_NONZERO={qc['n_nonzero']}")
    print(f"N_ZERO_FAM={qc['n_zero_fam']}")
    print(f"N_TESTS={qc['n_tests']}")


if __name__ == "__main__":
    main()
