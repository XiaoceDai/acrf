#!/usr/bin/env python3
"""
Figure 6B — Posterior configuration and lead-variant separation (publication).

Continuous diagnostic companion to Figure 6A. Frozen inputs only; no coloc
recomputation and no alteration of PP3 / PP4 values.

Authoritative inputs:
  results/figure6_mr_anchors_v3_pipeline_wide.tsv
  results/fig6a_v3_locus_zoom_summary.tsv
  results/fig6a_anchor_coloc_calibration_source.tsv
  results/fig6b_kcnn3_distance_recovery_audit.tsv

Outputs:
  figures/fig6b_posterior_lead_distance.{pdf,svg,png}
  figures/_pdf_preview/fig6b_posterior_lead_distance.png
  results/fig6b_posterior_lead_distance_{source,manifest,qc}.*
"""
from __future__ import annotations

import hashlib
import math
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
    FONT_FLOOR_HARD,
    FS_AXIS_LABEL,
    FS_FOOTNOTE,
    FS_LEGEND,
    FS_PANEL_TITLE,
    FS_TICK,
    MUTED_COLOR,
    set_nature_style,
)
from _figure_export import save_figure_bundle  # noqa: E402
from compute_anchor_lead_distance import (  # noqa: E402
    ANCHORS,
    EQTL_NAMES,
    GWAS_NAMES,
    _cache_paths,
    _read_slice,
)

PIPELINE_WIDE = ROOT / "results" / "figure6_mr_anchors_v3_pipeline_wide.tsv"
LOCUS_SUMMARY = ROOT / "results" / "fig6a_v3_locus_zoom_summary.tsv"
FIG6A_SOURCE = ROOT / "results" / "fig6a_anchor_coloc_calibration_source.tsv"
KCNN3_AUDIT = ROOT / "results" / "fig6b_kcnn3_distance_recovery_audit.tsv"

OUT_PDF = ROOT / "figures" / "fig6b_posterior_lead_distance.pdf"
OUT_SVG = ROOT / "figures" / "fig6b_posterior_lead_distance.svg"
OUT_PNG = ROOT / "figures" / "fig6b_posterior_lead_distance.png"
OUT_PREVIEW = ROOT / "figures" / "_pdf_preview" / "fig6b_posterior_lead_distance.png"
OUT_SOURCE = ROOT / "results" / "fig6b_posterior_lead_distance_source.tsv"
OUT_MANIFEST = ROOT / "results" / "fig6b_posterior_lead_distance_manifest.tsv"
OUT_QC = ROOT / "results" / "fig6b_posterior_lead_distance_qc.md"

TITLE = "B  Posterior configuration and lead-variant separation"
SUBTITLE = (
    "High PP3, low PP4, and widely separated lead variants across "
    "prespecified anchors"
)
FOOTER_1 = (
    "Lead separation is descriptive and is not used in the coloc posterior."
)
FOOTER_2 = (
    "High PP3 with low PP4 favors distinct local signals over supported "
    "same-variant mediation."
)

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
    FS_FOOT,
    FS_SOURCE,
    FS_TITLE,
    write_ab_geometry,
)

WIDTH_MM = AB_WIDTH
HEIGHT_MM = AB_HEIGHT
MARGIN_L = AB_FAMILY_LEFT
FAMILY_W = AB_FAMILY_W
ANCHOR_W = 22.1
# Remaining width after locked family/anchor columns + right margin
_RIGHT_MARGIN = 1.5
_REMAIN = WIDTH_MM - AB_ANCHOR_LEFT - ANCHOR_W - _RIGHT_MARGIN
# Three tracks + gutters between gene|pp3|pp4|dist
N_GUTTERS = 3
GUTTER = 1.2
PP3_W = (_REMAIN - N_GUTTERS * GUTTER) * 0.28
PP4_W = (_REMAIN - N_GUTTERS * GUTTER) * 0.28
DIST_W = (_REMAIN - N_GUTTERS * GUTTER) * 0.44

PP4_SUPPORT = 0.70
POINT_SIZE = 16.0
STEM_LW = 0.50

ACCENT_PP3 = "#8A7348"
ACCENT_PP4 = "#6B5A6E"
ACCENT_DIST = "#5C6B7A"
BAND_CAL = "#F3F2EF"
SMAD9_HL = "#FAF6F0"   # ~10% lighter than prior warm tint
ROW_LINE = "#F0F0F0"   # ~15–20% lighter separators
RULE = "#D0CDC8"
SUPPORT_BAND = "#EFEAF1"
OPEN_NA = "#8A8A8A"
FAMILY_LABEL = "#6E6E6E"
BRACKET = "#D4CFC7"

ANCHOR_ORDER = [
    "PRELP", "OMD", "IDH2", "NDUFA8", "BMP10", "SMAD9", "SMAD6",
    "NTRK2", "NTF3", "FRS2", "PITX2", "KCNN3",
]
CALIBRATORS = {"PITX2", "KCNN3"}

FORBIDDEN_VISIBLE = [
    "Miyazawa", "miyazawa", "mr_", "figure6_mr", "NCP_v3", "ncp_v3",
    "459–1214", "459-1214", "point size",
]

GROUP_FULL = {
    "F10": "F10 / GAG disease / TGF-beta",
    "Metabolic": "Metabolic",
    "F8": "F8 / Complex I biogenesis",
    "F11": "F11 / BMP signalling",
    "F15": "F15 / NTRK2 / TrkB",
    "AF calibration controls": "AF calibration controls",
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _stop(msg: str) -> None:
    raise SystemExit(f"STOP: {msg}")


def _lead_from_cache(symbol: str) -> dict[str, object]:
    """Documented min-p lead rule on local cache only (no remote)."""
    meta = ANCHORS[symbol]
    eq_path, gw_path = _cache_paths(symbol, meta)
    if not eq_path.is_file() or not gw_path.is_file():
        return {}
    eq = _read_slice(eq_path, EQTL_NAMES)
    gw = _read_slice(gw_path, GWAS_NAMES)
    eq = eq[eq["gene_id"].astype(str).str.startswith(meta["gene_id"])].copy()
    eq["position"] = pd.to_numeric(eq["position"], errors="coerce")
    eq["pvalue"] = pd.to_numeric(eq["pvalue"], errors="coerce")
    eq = eq.dropna(subset=["position", "pvalue"])
    if eq.empty:
        return {}
    eq_lead = eq.loc[eq["pvalue"].idxmin()]
    gw["base_pair_location"] = pd.to_numeric(gw["base_pair_location"], errors="coerce")
    gw["p_value"] = pd.to_numeric(gw["p_value"], errors="coerce")
    gw = gw.dropna(subset=["base_pair_location", "p_value"])
    if gw.empty:
        return {}
    gw_lead = gw.loc[gw["p_value"].idxmin()]
    return {
        "lead_gwas_variant": str(
            gw_lead.get("rsid") or gw_lead.get("rs_id") or gw_lead.get("variant_id") or ""
        ),
        "lead_gwas_bp": int(gw_lead["base_pair_location"]),
        "lead_eqtl_variant": str(eq_lead.get("rsid") or eq_lead.get("variant") or ""),
        "lead_eqtl_bp": int(eq_lead["position"]),
    }


def _read_kcnn3_audit() -> pd.Series:
    if not KCNN3_AUDIT.is_file():
        _stop(f"missing KCNN3 recovery audit: {KCNN3_AUDIT}")
    audit = pd.read_csv(KCNN3_AUDIT, sep="\t").iloc[0]
    if str(audit["gene"]) != "KCNN3":
        _stop("KCNN3 audit gene mismatch")
    status = str(audit["distance_status"])
    if status not in ("frozen_existing", "frozen_reconstructed", "unavailable_with_reason"):
        _stop(f"invalid KCNN3 distance_status: {status}")
    return audit


def build_source_table() -> tuple[pd.DataFrame, dict[str, object]]:
    for path in (PIPELINE_WIDE, LOCUS_SUMMARY, FIG6A_SOURCE):
        if not path.is_file():
            _stop(f"missing required source: {path}")

    wide = pd.read_csv(PIPELINE_WIDE, sep="\t").set_index("anchor")
    loc = pd.read_csv(LOCUS_SUMMARY, sep="\t").set_index("anchor")
    fig6a = pd.read_csv(FIG6A_SOURCE, sep="\t")
    kcnn3_audit = _read_kcnn3_audit()

    a6 = list(fig6a["gene"].astype(str))
    if a6 != ANCHOR_ORDER:
        _stop(f"Figure 6A order mismatch: {a6}")
    if len(set(a6)) != 12:
        _stop("Figure 6A source does not contain 12 unique anchors")
    missing_pp = [g for g in ANCHOR_ORDER if g not in wide.index]
    if missing_pp:
        _stop(f"PP3/PP4 anchors missing: {missing_pp}")
    missing_loc = [g for g in ANCHOR_ORDER if g not in loc.index]
    if missing_loc:
        _stop(f"lead-distance anchors missing from locus summary: {missing_loc}")

    source_sha = _sha256(PIPELINE_WIDE)
    locus_sha = _sha256(LOCUS_SUMMARY)

    rows = []
    for order_i, g in enumerate(ANCHOR_ORDER, start=1):
        r = wide.loc[g]
        arow = fig6a.loc[fig6a["gene"] == g].iloc[0]
        if pd.isna(r["coloc_PP3"]) or pd.isna(r["coloc_PP4"]):
            _stop(f"non-finite PP3/PP4 for {g}")
        pp3 = float(r["coloc_PP3"])
        pp4 = float(r["coloc_PP4"])
        if not (0.0 <= pp3 <= 1.0 and 0.0 <= pp4 <= 1.0):
            _stop(f"PP3/PP4 out of [0,1] for {g}")
        if not np.isclose(pp3, float(arow["coloc_PP3"]), atol=1e-12) or not np.isclose(
            pp4, float(arow["coloc_PP4"]), atol=1e-12
        ):
            _stop(f"PP3/PP4 disagree across sources for {g}")

        lr = loc.loc[g]
        if pd.isna(lr["delta_bp"]):
            _stop(f"frozen delta_bp missing for {g}")
        delta_bp = float(lr["delta_bp"])
        meta = ANCHORS[g]
        ref = (meta["start"] + meta["end"]) / 2.0
        lead_gwas_kb = float(lr["lead_gwas_kb"])
        lead_eqtl_kb = float(lr["lead_eqtl_kb"])
        gwas_bp = int(round(ref + lead_gwas_kb * 1000.0))
        eqtl_bp = int(round(ref + lead_eqtl_kb * 1000.0))
        if abs(abs(gwas_bp - eqtl_bp) - delta_bp) > 1.0:
            _stop(
                f"delta_bp inconsistent with lead_*_kb for {g}: "
                f"|{gwas_bp}-{eqtl_bp}| vs {delta_bp}"
            )

        # Prefer absolute coordinates from local cache when the documented
        # min-p rule reproduces the frozen delta_bp (kb fields are rounded).
        cache_lead = _lead_from_cache(g)
        gwas_var = ""
        eqtl_var = ""
        deriv = (
            "results/fig6a_v3_locus_zoom_summary.tsv; absolute bp from "
            "frozen lead_*_kb + locked GRCh38 gene midpoint"
        )
        if cache_lead:
            cache_delta = abs(
                int(cache_lead["lead_gwas_bp"]) - int(cache_lead["lead_eqtl_bp"])
            )
            if abs(cache_delta - delta_bp) <= 1:
                gwas_bp = int(cache_lead["lead_gwas_bp"])
                eqtl_bp = int(cache_lead["lead_eqtl_bp"])
                gwas_var = str(cache_lead["lead_gwas_variant"])
                eqtl_var = str(cache_lead["lead_eqtl_variant"])
                deriv = (
                    "results/fig6a_v3_locus_zoom_summary.tsv; absolute bp / "
                    "variants from local Project_P6/cache (documented min-p "
                    "rule; delta_bp matches freeze)"
                )
            else:
                _stop(
                    f"cache-derived delta conflicts with frozen delta_bp for {g}: "
                    f"{cache_delta} vs {delta_bp}"
                )

        if g == "KCNN3":
            status = str(kcnn3_audit["distance_status"])
            if status == "unavailable_with_reason":
                available = False
                delta_mb = np.nan
            else:
                available = True
                delta_mb = float(kcnn3_audit["delta_mb"])
                # Prefer audit absolute bp / variants when present
                if pd.notna(kcnn3_audit.get("lead_gwas_bp")):
                    gwas_bp = int(kcnn3_audit["lead_gwas_bp"])
                if pd.notna(kcnn3_audit.get("lead_eqtl_bp")):
                    eqtl_bp = int(kcnn3_audit["lead_eqtl_bp"])
                if str(kcnn3_audit.get("lead_gwas_variant", "")):
                    gwas_var = str(kcnn3_audit["lead_gwas_variant"])
                if str(kcnn3_audit.get("lead_eqtl_variant", "")):
                    eqtl_var = str(kcnn3_audit["lead_eqtl_variant"])
                if abs(float(kcnn3_audit["delta_bp"]) - delta_bp) > 1e-6:
                    _stop("KCNN3 audit delta_bp disagrees with locus summary")
                deriv = str(kcnn3_audit["derivation_source"])
            if str(kcnn3_audit.get("genome_build", "")) not in ("", "GRCh38", "nan"):
                if str(kcnn3_audit["genome_build"]) != "GRCh38":
                    _stop(f"KCNN3 genome_build not GRCh38: {kcnn3_audit['genome_build']}")
        else:
            status = "frozen_existing"
            available = True
            delta_mb = abs(delta_bp) / 1e6

        rows.append({
            "anchor_order": order_i,
            "family_or_role": str(arow["group_compact"]),
            "gene": g,
            "coloc_PP3": pp3,
            "coloc_PP4": pp4,
            "pp3_dominant": bool(pp3 > pp4 and pp4 < PP4_SUPPORT),
            "pp4_support": bool(pp4 >= PP4_SUPPORT),
            "lead_gwas_variant": gwas_var,
            "lead_gwas_bp": gwas_bp,
            "lead_eqtl_variant": eqtl_var,
            "lead_eqtl_bp": eqtl_bp,
            "delta_bp": delta_bp,
            "delta_mb": delta_mb if available else np.nan,
            "distance_status": status if g == "KCNN3" else "frozen_existing",
            "distance_derivation_source": deriv,
            "distance_available": available,
            "source_table": str(PIPELINE_WIDE.relative_to(ROOT)),
            "source_sha256": source_sha,
            "locus_source_sha256": locus_sha,
            "anchor_type": "calibration" if g in CALIBRATORS else "framework",
        })

    df = pd.DataFrame(rows)
    stats = _assert_source(df, kcnn3_audit)
    return df, stats


def _assert_source(df: pd.DataFrame, kcnn3_audit: pd.Series) -> dict[str, object]:
    if len(df) != 12 or df["gene"].nunique() != 12:
        _stop("expected exactly 12 unique anchors")
    if list(df["gene"]) != ANCHOR_ORDER:
        _stop("anchor order does not match Figure 6A")

    pp3 = df["coloc_PP3"].astype(float)
    pp4 = df["coloc_PP4"].astype(float)
    if not np.isfinite(pp3).all() or not np.isfinite(pp4).all():
        _stop("non-finite PP3/PP4")
    if ((pp3 < 0) | (pp3 > 1) | (pp4 < 0) | (pp4 > 1)).any():
        _stop("PP3/PP4 outside [0,1]")

    n_pp3_dom = int(df["pp3_dominant"].sum())
    n_pp4 = int(df["pp4_support"].sum())
    if n_pp3_dom != 12:
        _stop(f"PP3-dominant count {n_pp3_dom}/12 != 12/12")
    if n_pp4 != 0:
        _stop(f"PP4 support count {n_pp4}/12 != 0/12")

    avail = df.loc[df["distance_available"]]
    missing = df.loc[~df["distance_available"], "gene"].tolist()
    k_status = str(kcnn3_audit["distance_status"])
    if k_status in ("frozen_existing", "frozen_reconstructed"):
        if "KCNN3" in missing or not bool(df.loc[df["gene"] == "KCNN3", "distance_available"].iloc[0]):
            _stop("KCNN3 recovered but still marked unavailable")
        if len(avail) != 12:
            _stop(f"expected 12 available distances after recovery, got {len(avail)}")
    else:
        if missing != ["KCNN3"]:
            _stop(f"unexpected missing distances: {missing}")

    dmin = float(avail["delta_mb"].min())
    dmax = float(avail["delta_mb"].max())
    raw_min = float(avail["delta_bp"].abs().min()) / 1e6
    raw_max = float(avail["delta_bp"].abs().max()) / 1e6
    if not np.isclose(dmin, raw_min) or not np.isclose(dmax, raw_max):
        _stop("distance min/max disagree with frozen delta_bp")

    upper = max(2.0, math.ceil(dmax / 0.5) * 0.5)
    return {
        "n_pp3_dom": n_pp3_dom,
        "n_pp4": n_pp4,
        "n_dist_available": int(len(avail)),
        "dist_min_mb": dmin,
        "dist_max_mb": dmax,
        "dist_upper_mb": upper,
        "unavailable": missing,
        "kcnn3_status": k_status,
        "kcnn3_delta_mb": (
            float(df.loc[df["gene"] == "KCNN3", "delta_mb"].iloc[0])
            if bool(df.loc[df["gene"] == "KCNN3", "distance_available"].iloc[0])
            else None
        ),
        "kcnn3_gwas": str(df.loc[df["gene"] == "KCNN3", "lead_gwas_variant"].iloc[0]),
        "kcnn3_eqtl": str(df.loc[df["gene"] == "KCNN3", "lead_eqtl_variant"].iloc[0]),
    }


def runtime_summary(stats: dict[str, object]) -> str:
    dmin = stats["dist_min_mb"]
    dmax = stats["dist_max_mb"]
    base = (
        f"PP3-dominant {stats['n_pp3_dom']}/12 · "
        f"PP4 ≥ 0.70 support {stats['n_pp4']}/12 · "
        f"lead separation {dmin:.2f}–{dmax:.2f} Mb"
    )
    if stats["unavailable"]:
        base += " · KCNN3 unavailable"
    return base


def _family_groups(df: pd.DataFrame) -> list[tuple[str, list[int]]]:
    groups: list[tuple[str, list[int]]] = []
    cur: str | None = None
    idxs: list[int] = []
    for i, r in df.iterrows():
        lab = str(r["family_or_role"])
        if lab != cur:
            if cur is not None:
                groups.append((cur, idxs))
            cur = lab
            idxs = [int(i)]
        else:
            idxs.append(int(i))
    if cur is not None:
        groups.append((cur, idxs))
    return groups


def render(df: pd.DataFrame, summary: str, stats: dict[str, object]) -> plt.Figure:
    set_nature_style()
    fig = plt.figure(
        figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4),
        facecolor="white", dpi=300,
    )
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, WIDTH_MM)
    ax.set_ylim(0, HEIGHT_MM)
    ax.set_axis_off()

    x0 = MARGIN_L
    x_fam = AB_FAMILY_LEFT
    x_gene = AB_ANCHOR_LEFT
    x_pp3 = x_gene + ANCHOR_W + GUTTER
    x_pp4 = x_pp3 + PP3_W + GUTTER
    x_dist = x_pp4 + PP4_W + GUTTER
    span_w = (x_dist + DIST_W) - x_fam

    # Shared A/B two-line title contract: title + runtime summary
    from figure6_composite_geometry import AB_LINE2_Y, AB_TITLE_Y
    ax.text(
        x0, AB_TITLE_Y, TITLE, fontsize=FS_TITLE, fontweight="bold",
        ha="left", va="top", color=FG_COLOR, clip_on=False,
    )
    ax.text(
        x0, AB_LINE2_Y, summary, fontsize=FS_SOURCE,
        ha="left", va="top", color="#3A3A3A", clip_on=False,
    )

    # Shared A/B body geometry (do not recompute from free title pack)
    head_h = AB_HEAD_H
    body_top = AB_BODY_TOP
    tracks_bot = AB_BODY_BOT
    tracks_top = body_top + head_h
    n_row = len(df)
    row_h = AB_ROW_H
    n_framework = int((df["anchor_type"] == "framework").sum())
    dist_upper = float(stats["dist_upper_mb"])
    row_centers = list(AB_ROW_CENTERS)

    # Strict two-level header
    ax.add_patch(Rectangle(
        (x_fam, body_top), span_w, head_h,
        facecolor="#F7F6F4", edgecolor="none", zorder=0,
    ))
    ax.plot(
        [x_fam, x_dist + DIST_W], [body_top, body_top],
        color=RULE, lw=0.45, zorder=2,
    )

    # Level 1 — primary column titles (upper band)
    y_l1 = body_top + 7.15
    ax.text(
        x_fam + FAMILY_W / 2, y_l1, "Family / role",
        ha="center", va="center", fontsize=FS_TICK - 0.3,
        fontweight="bold", color=FG_COLOR, zorder=3,
    )
    ax.text(
        x_gene + ANCHOR_W / 2, y_l1, "Anchor",
        ha="center", va="center", fontsize=FS_TICK - 0.3,
        fontweight="bold", color=FG_COLOR, zorder=3,
    )
    ax.text(
        x_pp3 + PP3_W / 2, y_l1, "PP3",
        ha="center", va="center", fontsize=FS_TICK,
        fontweight="bold", color=ACCENT_PP3, zorder=3,
    )
    ax.text(
        x_pp4 + PP4_W / 2, y_l1, "PP4",
        ha="center", va="center", fontsize=FS_TICK,
        fontweight="bold", color=ACCENT_PP4, zorder=3,
    )
    ax.text(
        x_dist + DIST_W / 2, y_l1, "Lead GWAS–eQTL separation",
        ha="center", va="center", fontsize=FS_TICK - 0.2,
        fontweight="bold", color=ACCENT_DIST, zorder=3,
    )

    # Level 2 — subtitles under numerical tracks (clear of threshold label)
    y_l2 = body_top + 4.55
    ax.text(
        x_pp3 + PP3_W / 2, y_l2, "distinct local signals",
        ha="center", va="center", fontsize=FONT_FLOOR_HARD,
        color=MUTED_COLOR, zorder=3,
    )
    # Keep PP4 subtitle left of the support band
    ax.text(
        x_pp4 + PP4_W * 0.38, y_l2, "shared causal signal",
        ha="center", va="center", fontsize=FONT_FLOOR_HARD,
        color=MUTED_COLOR, zorder=3,
    )
    ax.text(
        x_dist + DIST_W / 2, y_l2, "Mb",
        ha="center", va="center", fontsize=FONT_FLOOR_HARD,
        color=MUTED_COLOR, zorder=3,
    )

    def x_pp(track_x: float, track_w: float, v: float) -> float:
        return track_x + 3.2 + (track_w - 5.8) * float(v)

    def x_mb(v: float) -> float:
        # Leave right pad for left-placed high-distance labels
        return x_dist + 5.5 + (DIST_W - 12.0) * (float(v) / dist_upper)

    # PP4 support band in body only
    band_x0 = x_pp(x_pp4, PP4_W, PP4_SUPPORT)
    band_x1 = x_pp(x_pp4, PP4_W, 1.0)
    ax.add_patch(Rectangle(
        (band_x0, tracks_bot), band_x1 - band_x0, body_top - tracks_bot,
        facecolor=SUPPORT_BAND, edgecolor="none", zorder=0, alpha=0.85,
    ))
    ax.plot(
        [band_x0, band_x0], [tracks_bot, body_top],
        color=ACCENT_PP4, lw=0.5, linestyle=(0, (2, 1.4)), zorder=2,
    )
    # Threshold annotation: bottom of header, centered on the support band only.
    ax.text(
        (band_x0 + band_x1) / 2, body_top + 2.35,
        "PP4 ≥ 0.70",
        ha="center", va="center", fontsize=FONT_FLOOR_HARD,
        color=ACCENT_PP4, zorder=5, clip_on=False,
    )
    ax.text(
        (band_x0 + band_x1) / 2, body_top + 0.85,
        "support region",
        ha="center", va="center", fontsize=FONT_FLOOR_HARD,
        color="#8A7A8C", zorder=5, clip_on=False,
    )

    row_y: dict[int, float] = {}
    dist_label_artists: list[tuple[str, float, float, str]] = []
    for i in range(n_row):
        row_y[i] = row_centers[i]
        y_cursor = row_y[i] - row_h / 2
        r = df.iloc[i]

        if r["anchor_type"] == "calibration":
            ax.add_patch(Rectangle(
                (x_fam, y_cursor), span_w, row_h,
                facecolor=BAND_CAL, edgecolor="none", zorder=0,
            ))
        elif r["gene"] == "SMAD9":
            ax.add_patch(Rectangle(
                (x_fam, y_cursor), span_w, row_h,
                facecolor=SMAD9_HL, edgecolor="none", zorder=0,
            ))

        ax.plot(
            [x_fam + 0.3, x_dist + DIST_W - 0.3], [y_cursor, y_cursor],
            color=ROW_LINE, lw=0.24, zorder=1,
        )

        ax.text(
            x_gene + 1.2, row_y[i], r["gene"],
            ha="left", va="center",
            fontsize=FS_TICK - 0.05, fontweight="bold",
            color="#5C4030" if r["gene"] == "SMAD9" else FG_COLOR,
            zorder=4,
        )

        # PP3
        x3 = x_pp(x_pp3, PP3_W, r["coloc_PP3"])
        x3_0 = x_pp(x_pp3, PP3_W, 0.0)
        ax.plot([x3_0, x3], [row_y[i], row_y[i]], color=ACCENT_PP3, lw=STEM_LW, zorder=3)
        ax.scatter(
            [x3], [row_y[i]], s=POINT_SIZE, c=ACCENT_PP3,
            linewidths=0, zorder=4, clip_on=False,
        )

        # PP4
        x4 = x_pp(x_pp4, PP4_W, r["coloc_PP4"])
        x4_0 = x_pp(x_pp4, PP4_W, 0.0)
        ax.plot([x4_0, x4], [row_y[i], row_y[i]], color=ACCENT_PP4, lw=STEM_LW, zorder=3)
        ax.scatter(
            [x4], [row_y[i]], s=POINT_SIZE, c=ACCENT_PP4,
            linewidths=0, zorder=4, clip_on=False,
        )

        # Distance — every available value labeled; dynamic side placement
        if bool(r["distance_available"]):
            dmb = float(r["delta_mb"])
            xd = x_mb(dmb)
            xd0 = x_mb(0.0)
            ax.plot([xd0, xd], [row_y[i], row_y[i]], color=ACCENT_DIST, lw=STEM_LW, zorder=3)
            ax.scatter(
                [xd], [row_y[i]], s=POINT_SIZE, c=ACCENT_DIST,
                linewidths=0, zorder=4, clip_on=False,
            )
            lab = f"{dmb:.2f}"
            if dmb < 1.70:
                lx, ha = xd + 1.35, "left"
            else:
                lx, ha = xd - 1.35, "right"
            ax.text(
                lx, row_y[i], lab,
                ha=ha, va="center", fontsize=FS_AXIS_LABEL,
                color="#7A7A7A", zorder=5, clip_on=False,
            )
            dist_label_artists.append((str(r["gene"]), lx, row_y[i], lab))
        else:
            ax.text(
                x_dist + DIST_W / 2, row_y[i], "N/A",
                ha="center", va="center", fontsize=FS_TICK - 0.2,
                color=OPEN_NA, zorder=4,
            )

    # Store for QC via figure canvas metadata
    fig._fig6b_dist_labels = dist_label_artists  # type: ignore[attr-defined]
    fig._fig6b_row_y = row_y  # type: ignore[attr-defined]

    # Compact family brackets (lighter weight)
    for lab, idxs in _family_groups(df):
        y0 = row_y[idxs[-1]] - row_h / 2 + 0.12
        y1 = row_y[idxs[0]] + row_h / 2 - 0.12
        bx = x_fam + 1.0
        ax.plot([bx, bx], [y0, y1], color=BRACKET, lw=0.75, solid_capstyle="round", zorder=2)
        ax.plot([bx, bx + 0.7], [y0, y0], color=BRACKET, lw=0.65, zorder=2)
        ax.plot([bx, bx + 0.7], [y1, y1], color=BRACKET, lw=0.65, zorder=2)
        y_mid = 0.5 * (row_y[idxs[0]] + row_y[idxs[-1]])
        if lab == "AF calibration controls":
            ax.text(
                x_fam + FAMILY_W / 2 + 0.5, y_mid + 1.05, "AF calibration",
                ha="center", va="center", fontsize=FONT_FLOOR_HARD,
                color=FAMILY_LABEL, zorder=4,
            )
            ax.text(
                x_fam + FAMILY_W / 2 + 0.5, y_mid - 1.05, "controls",
                ha="center", va="center", fontsize=FONT_FLOOR_HARD,
                color=FAMILY_LABEL, zorder=4,
            )
        else:
            ax.text(
                x_fam + FAMILY_W / 2 + 0.5, y_mid, lab,
                ha="center", va="center", fontsize=FS_TICK - 0.2,
                fontweight="bold", color=FAMILY_LABEL, zorder=4,
            )

    if n_framework < n_row:
        y_div = AB_CAL_SEP_Y
        ax.plot(
            [x_fam + 0.2, x_dist + DIST_W - 0.2], [y_div, y_div],
            color="#8A8580", lw=0.8, zorder=3,
        )

    # Axes — shared baseline; ≥4.5 mm to first footer line
    axis_y = tracks_bot + 0.35
    footer1_y = 4.7
    footer2_y = 2.4

    def _pp_axis(track_x: float, track_w: float) -> None:
        x_lo = x_pp(track_x, track_w, 0.0)
        x_hi = x_pp(track_x, track_w, 1.0)
        ax.plot([x_lo, x_hi], [axis_y, axis_y], color=RULE, lw=0.4, zorder=2)
        for t, lab in [(0.0, "0"), (0.5, "0.5"), (1.0, "1.0")]:
            xt = x_pp(track_x, track_w, t)
            ax.plot([xt, xt], [axis_y, axis_y + 0.45], color=RULE, lw=0.35, zorder=2)
            ax.text(
                xt, axis_y - 0.15, lab, ha="center", va="top",
                fontsize=FS_AXIS_LABEL, color=MUTED_COLOR, zorder=3,
            )

    _pp_axis(x_pp3, PP3_W)
    _pp_axis(x_pp4, PP4_W)

    ticks = [0.0, 0.5, 1.0, 1.5, 2.0]
    x_lo = x_mb(0.0)
    x_hi = x_mb(dist_upper)
    ax.plot([x_lo, x_hi], [axis_y, axis_y], color=RULE, lw=0.4, zorder=2)
    for t in ticks:
        if t > dist_upper + 1e-9:
            continue
        xt = x_mb(t)
        ax.plot([xt, xt], [axis_y, axis_y + 0.45], color=RULE, lw=0.35, zorder=2)
        ax.text(
            xt, axis_y - 0.15, f"{t:g}", ha="center", va="top",
            fontsize=FS_AXIS_LABEL, color=MUTED_COLOR, zorder=3,
        )

    ax.text(
        x0, footer1_y, FOOTER_1,
        ha="left", va="center", fontsize=FONT_FLOOR_HARD,
        color=MUTED_COLOR, clip_on=False,
    )
    ax.text(
        x0, footer2_y, FOOTER_2,
        ha="left", va="center", fontsize=FONT_FLOOR_HARD,
        color=MUTED_COLOR, clip_on=False,
    )
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


def run_qc(df: pd.DataFrame, stats: dict[str, object], summary: str) -> dict[str, object]:
    for path in (OUT_PDF, OUT_SVG, OUT_PNG, OUT_SOURCE, OUT_PREVIEW):
        assert path.is_file(), path

    src = pd.read_csv(OUT_SOURCE, sep="\t")
    assert list(src["gene"]) == ANCHOR_ORDER
    assert int(src["pp3_dominant"].astype(str).str.lower().isin(["true", "1"]).sum()) == 12
    assert int(src["pp4_support"].astype(str).str.lower().isin(["true", "1"]).sum()) == 0
    assert np.allclose(src["coloc_PP3"].astype(float), df["coloc_PP3"].astype(float), rtol=0, atol=0)
    assert np.allclose(src["coloc_PP4"].astype(float), df["coloc_PP4"].astype(float), rtol=0, atol=0)

    avail_mask = src["distance_status"].astype(str).isin(
        ["frozen_existing", "frozen_reconstructed"]
    )
    # Also accept distance_available column if present
    if "distance_available" in src.columns:
        avail_mask = src["distance_available"].astype(str).str.lower().isin(["true", "1"])
    avail = src.loc[avail_mask]
    assert len(avail) == stats["n_dist_available"]
    dmin = float(avail["delta_mb"].min())
    dmax = float(avail["delta_mb"].max())
    assert abs(dmin - stats["dist_min_mb"]) < 1e-12
    assert abs(dmax - stats["dist_max_mb"]) < 1e-12
    assert f"{dmin:.2f}" in summary and f"{dmax:.2f}" in summary
    assert "PP4 ≥ 0.70" in summary
    if stats["kcnn3_status"] in ("frozen_existing", "frozen_reconstructed"):
        assert "KCNN3 unavailable" not in summary
        pdf_txt = subprocess.run(
            ["pdftotext", "-layout", str(OUT_PDF), "-"],
            check=True, capture_output=True, text=True,
        ).stdout
        for line in pdf_txt.splitlines():
            if "KCNN3" in line:
                assert "N/A" not in line, line
        krow = src.loc[src["gene"] == "KCNN3"].iloc[0]
        assert float(krow["delta_mb"]) > 0
    # Publication typography: SVG must contain ≥, not ASCII >=
    svg_text = OUT_SVG.read_text(encoding="utf-8")
    assert "≥" in svg_text
    assert "PP4 &gt;= 0.70" not in svg_text
    assert "PP4 >= 0.70" not in svg_text

    visible = _extract_visible_text(OUT_PDF, OUT_SVG)
    for bad in FORBIDDEN_VISIBLE:
        assert bad not in visible, f"forbidden visible text: {bad}"
    assert not re.search(r"(?<![A-Za-z])MR(?![A-Za-z])", visible)
    assert "PP4 ≥ 0.70" in visible
    assert "support region" in visible
    assert "AF calibration" in visible
    assert "same-variant mediation" in visible
    assert "GWAS–eQTL" in visible or "GWAS" in visible
    assert "459–1214" not in visible
    assert "point size" not in visible.lower()
    for g in ANCHOR_ORDER:
        assert g in visible, f"missing gene {g}"
    assert "GAG disease" not in visible
    assert ".tsv" not in visible

    # 12 visible numeric distance labels (two-decimal Mb)
    expected_labs = [f"{float(v):.2f}" for v in avail["delta_mb"].tolist()]
    assert len(expected_labs) == 12
    for lab in expected_labs:
        assert lab in visible, f"missing distance label {lab}"
    assert "1.92" in visible  # BMP10
    assert "0.46" in visible  # KCNN3
    assert "N/A" not in visible

    # Label positions from last render metadata (regenerated below in main)
    # Approximate in-canvas check via SVG text content positions if present
    svg_text = OUT_SVG.read_text(encoding="utf-8")
    # Count distance-label occurrences (unique values may collide; require ≥12
    # matches of xx.xx pattern near distance track is soft — hard count via
    # expected_labs uniqueness-aware check above).
    assert expected_labs.count("1.92") >= 1
    assert expected_labs.count("0.46") >= 1

    pdf_w, pdf_h = _pdf_page_mm(OUT_PDF)
    svg_w, svg_h = _svg_mm(OUT_SVG)
    png_w, png_h, png_px_w, png_px_h = _png_mm(OUT_PNG, dpi=600)
    assert abs(svg_w - png_w) < 0.50
    assert abs(svg_h - png_h) < 0.50
    assert abs(pdf_w - svg_w) < 2.0
    assert abs(pdf_h - svg_h) < 2.0
    assert abs(pdf_w - WIDTH_MM) < 2.0
    assert abs(pdf_h - HEIGHT_MM) < 1.5, pdf_h

    # Composite typography floor is FS_FOOT (5.2–5.7 pt contract), not the
    # legacy FONT_FLOOR_HARD = 6.5 used for standalone denser panels.
    sizes = [float(x) for x in re.findall(r"font-size:\s*([0-9.]+)px", svg_text)]
    sizes += [float(x) for x in re.findall(r"font-size:\s*([0-9.]+)", svg_text)]
    if sizes:
        assert min(sizes) >= FS_FOOT - 0.05, min(sizes)

    return {
        "summary": summary,
        "stats": stats,
        "pdf_mm": (round(pdf_w, 2), round(pdf_h, 2)),
        "svg_mm": (round(svg_w, 2), round(svg_h, 2)),
        "png_mm": (round(png_w, 2), round(png_h, 2)),
        "png_px": (png_px_w, png_px_h),
        "min_font_px": min(sizes) if sizes else None,
        "n_distance_labels": len(expected_labs),
        "distance_labels": expected_labs,
        "source_sha256": _sha256(OUT_SOURCE),
        "pdf_sha256": _sha256(OUT_PDF),
        "svg_sha256": _sha256(OUT_SVG),
        "png_sha256": _sha256(OUT_PNG),
    }


def write_manifest(qc: dict[str, object]) -> None:
    rows = [
        {"artifact": "figure_pdf", "path": str(OUT_PDF.relative_to(ROOT)),
         "sha256": qc["pdf_sha256"],
         "geometry_mm": f"{qc['pdf_mm'][0]}x{qc['pdf_mm'][1]}", "notes": "publication PDF"},
        {"artifact": "figure_svg", "path": str(OUT_SVG.relative_to(ROOT)),
         "sha256": qc["svg_sha256"],
         "geometry_mm": f"{qc['svg_mm'][0]}x{qc['svg_mm'][1]}", "notes": "publication SVG"},
        {"artifact": "figure_png_600dpi", "path": str(OUT_PNG.relative_to(ROOT)),
         "sha256": qc["png_sha256"],
         "geometry_mm": f"{qc['png_mm'][0]}x{qc['png_mm'][1]}",
         "notes": f"600 dpi; {qc['png_px'][0]}x{qc['png_px'][1]} px"},
        {"artifact": "manuscript_preview_png", "path": str(OUT_PREVIEW.relative_to(ROOT)),
         "sha256": _sha256(OUT_PREVIEW), "geometry_mm": "", "notes": "150 dpi preview"},
        {"artifact": "source_table", "path": str(OUT_SOURCE.relative_to(ROOT)),
         "sha256": qc["source_sha256"], "geometry_mm": "",
         "notes": "includes KCNN3 frozen distance after recovery audit"},
        {"artifact": "kcnn3_recovery_audit", "path": str(KCNN3_AUDIT.relative_to(ROOT)),
         "sha256": _sha256(KCNN3_AUDIT), "geometry_mm": "",
         "notes": str(qc["stats"]["kcnn3_status"])},
        {"artifact": "pp_source", "path": str(PIPELINE_WIDE.relative_to(ROOT)),
         "sha256": _sha256(PIPELINE_WIDE), "geometry_mm": "",
         "notes": "read-only; filename not shown on panel"},
        {"artifact": "locus_source", "path": str(LOCUS_SUMMARY.relative_to(ROOT)),
         "sha256": _sha256(LOCUS_SUMMARY), "geometry_mm": "",
         "notes": "frozen lead-distance including KCNN3"},
    ]
    pd.DataFrame(rows).to_csv(OUT_MANIFEST, sep="\t", index=False)


def write_qc(qc: dict[str, object]) -> None:
    s = qc["stats"]
    OUT_QC.write_text(
        "\n".join(
            [
                "# Figure 6B posterior / lead-distance QC",
                "",
                f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
                "",
                "Status: **PASS**",
                "",
                "## Scope",
                "- Publication Figure 6B polish only; Figure 6A/C/D untouched.",
                "- No coloc recomputation; no PP3/PP4 alteration.",
                "- No Figure 6 A–D composite routing change.",
                "",
                "## KCNN3 recovery",
                f"- distance_status = `{s['kcnn3_status']}`",
                f"- lead GWAS variant = `{s['kcnn3_gwas']}`",
                f"- lead eQTL variant = `{s['kcnn3_eqtl']}`",
                f"- delta_mb = `{s['kcnn3_delta_mb']}`",
                "- Prior N/A was an obsolete panel-policy omission; frozen locus "
                "summary already contained KCNN3.",
                "",
                "## Runtime summary",
                f"- `{qc['summary']}`",
                "",
                "## Counts",
                f"- PP3-dominant = {s['n_pp3_dom']}/12",
                f"- PP4 ≥ 0.70 = {s['n_pp4']}/12",
                f"- Available lead distances = {s['n_dist_available']}/12",
                f"- Visible numeric distance labels = {qc['n_distance_labels']}/12",
                f"- Distance range = {s['dist_min_mb']:.2f}–{s['dist_max_mb']:.2f} Mb",
                f"- Distance axis upper = {s['dist_upper_mb']:.1f} Mb",
                f"- Missing-distance rows = {s['unavailable'] or 'none'}",
                f"- Labels: {', '.join(qc['distance_labels'])}",
                "",
                "## Hard assertions",
                "- 12 unique anchors; order matches Figure 6A.",
                "- All PP3/PP4 finite in [0,1]; values unchanged.",
                "- 12/12 distance labels visible (BMP10 1.92; KCNN3 0.46).",
                "- Fixed point size; no area encoding; no 459–1214 kb text.",
                "- No visible MR / mr_ / Miyazawa / internal filename / NCP_v3 / N/A.",
                "",
                "## Layout polish",
                "- Two-level track header; PP4 threshold label = `PP4 ≥ 0.70`.",
                "- Dynamic distance-label side placement (left if ≥1.70 Mb).",
                "- Reduced separator / bracket / stem weight; lighter SMAD9 tint.",
                "",
                "## Full family names (legend)",
                *[f"- `{k}`: {v}" for k, v in GROUP_FULL.items()],
                "",
                "## Geometry / typography",
                f"- PDF: {qc['pdf_mm'][0]} × {qc['pdf_mm'][1]} mm",
                f"- SVG: {qc['svg_mm'][0]} × {qc['svg_mm'][1]} mm",
                f"- PNG: {qc['png_mm'][0]} × {qc['png_mm'][1]} mm "
                f"({qc['png_px'][0]} × {qc['png_px'][1]} px @ 600 dpi)",
                f"- Minimum detected SVG font-size: {qc['min_font_px']} px",
                "",
                "## SHA256",
                f"- source: `{qc['source_sha256']}`",
                f"- pdf: `{qc['pdf_sha256']}`",
                f"- svg: `{qc['svg_sha256']}`",
                f"- png: `{qc['png_sha256']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    df, stats = build_source_table()
    summary = runtime_summary(stats)

    OUT_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    out_cols = [
        "anchor_order", "family_or_role", "gene",
        "coloc_PP3", "coloc_PP4", "pp3_dominant", "pp4_support",
        "lead_gwas_variant", "lead_gwas_bp", "lead_eqtl_variant", "lead_eqtl_bp",
        "delta_bp", "delta_mb", "distance_status", "distance_derivation_source",
        "source_table", "source_sha256",
    ]
    df[out_cols].to_csv(OUT_SOURCE, sep="\t", index=False)
    print(f"wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"summary: {summary}")
    print(
        f"KCNN3 recovery: {stats['kcnn3_status']} "
        f"delta_mb={stats['kcnn3_delta_mb']} "
        f"gwas={stats['kcnn3_gwas']} eqtl={stats['kcnn3_eqtl']}"
    )

    # fontweight medium may not exist — map via matplotlib
    write_ab_geometry()
    from matplotlib.transforms import Bbox
    fixed = Bbox.from_bounds(0.0, 0.0, WIDTH_MM / 25.4, HEIGHT_MM / 25.4)
    fig = render(df, summary, stats)
    save_figure_bundle(fig, OUT_PDF, dpi=600, bbox_inches=fixed, pad_inches=0)
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")
    print(f"wrote {OUT_PNG.relative_to(ROOT)}")

    fig = render(df, summary, stats)
    OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        OUT_SVG, format="svg", bbox_inches=fixed, pad_inches=0,
        facecolor="white", edgecolor="none",
    )
    OUT_PREVIEW.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        OUT_PREVIEW, format="png", dpi=150, bbox_inches=fixed, pad_inches=0,
        facecolor="white", edgecolor="none",
    )
    plt.close(fig)
    print(f"wrote {OUT_SVG.relative_to(ROOT)}")
    print(f"wrote {OUT_PREVIEW.relative_to(ROOT)}")

    qc = run_qc(df, stats, summary)
    write_manifest(qc)
    write_qc(qc)
    print(f"wrote {OUT_MANIFEST.relative_to(ROOT)}")
    print(f"wrote {OUT_QC.relative_to(ROOT)}")
    print("QC: PASS")
    print(f"SOURCE_SHA256={qc['source_sha256']}")
    print(f"PDF_SHA256={qc['pdf_sha256']}")
    print(f"SVG_SHA256={qc['svg_sha256']}")
    print(f"PNG_SHA256={qc['png_sha256']}")


if __name__ == "__main__":
    main()
