#!/usr/bin/env python3
"""
Figure 6A (v3 follow-up) — Locus architecture: AF GWAS vs atrial cis-eQTL
                            for each of 10 anchors + 1 PITX2 calibrator.

Replaces the deprecated Fig 6A heatmap, which was structurally redundant
with Fig 2E (same families × cohorts × replication class).

Layout: 3 rows × 4 cols = 12 panels
        - 11 panels: one per (anchor + PITX2 calibrator) — locus-zoom mini
        - 12th panel: legend + global locus-architecture summary

Per-panel design (greyscale; no LA/RA red/blue reuse):
  * Top half:  AF GWAS  -log10(p) scatter, lead = filled diamond
  * Bottom half: atrial cis-eQTL -log10(p) scatter (mirror, downward),
                 lead = open diamond
  * Vertical dashed grey line at gene TSS
  * Horizontal Δ_lead arrow between GWAS-lead and eQTL-lead positions:
      Δ < 50 kb  → green  #1B7837
      50–250 kb  → yellow #FDB863
      ≥ 250 kb   → red    #B2182B   (PITX2-mode physical-distance failure)
  * Verdict pill (top-right) + eGene badge (bottom-right)

Inputs (read-only):
  results/figure6_mr_anchors_v3_consolidated.tsv
      anchor metadata + verdict + eligibility
  Project_P6/cache/gwas_<ANCHOR>_chr<N>_<lo>_<hi>.tsv.gz   × 11
  Project_P6/cache/eqtl_<ANCHOR>_chr<N>_<lo>_<hi>.tsv.gz   × 11
      previously cached by h4_mr_coloc_anchors.R / h4d_*.R / h4c_*.R

Output:
  figures/fig6a_v3_locus_zoom.pdf
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import (Arrow, FancyArrowPatch, FancyBboxPatch,
                                 Rectangle)

from plot_style import COL_DOUBLE, FG_COLOR, MUTED_COLOR, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
RES  = ROOT / "results"
CACHE = ROOT / "Project_P6" / "cache"

IN_CONSOL = RES / "figure6_mr_anchors_v3_consolidated.tsv"
OUT       = ROOT / "figures" / "fig6a_v3_locus_zoom.pdf"

# Cache slice column names (mirror h4d header definition)
EQTL_NAMES = ["molecular_trait_id", "chromosome", "position", "ref", "alt",
              "variant", "ma_samples", "maf", "pvalue", "beta", "se",
              "type", "ac", "an", "r2", "molecular_trait_object_id",
              "gene_id", "median_tpm", "rsid", "nlog10p"]
GWAS_NAMES = ["chromosome", "base_pair_location", "effect_allele",
              "other_allele", "beta", "standard_error",
              "effect_allele_frequency", "p_value", "rsid", "rs_id",
              "n", "hm_coordinate_conversion", "hm_code", "variant_id"]

# --- v3 follow-up #3 (visual refresh, 2026-04-23) -------------------------
# Style policy for this panel only (deliberate exemption from plot_style.py
# 6-pt floor): with 11 mini-locus-zooms in a 3x4 grid on a double-column
# page, no element exceeds ~45x35 mm so axis ticks must drop to 4-5 pt.
# All other figures retain the project floor.
#
# Colour palette switched to grayscale + a single muted brick-red accent
# for the Delta-lead arrow.  Verdict pill colours (4-class) are no longer
# used inside the panel data area: with 10/10 anchors landing in
# {no_eligible_eQTL, not_supported}, multi-colour pills carried no
# information and competed with the scatter for attention.  Single dark-
# grey pill is now drawn OUTSIDE the panel data area (figure coords) just
# above the title bar.
#
GRAY_GWAS_PT      = "#1A1A1A"   # GWAS scatter
GRAY_GWAS_LEAD    = "#000000"   # GWAS lead diamond (slightly darker)
GRAY_EQTL_PT      = "#B0B0B0"   # eQTL scatter (mid-grey)
GRAY_EQTL_EDGE    = "#666666"   # eQTL lead diamond edge
DELTA_BRICK       = "#C44E52"   # Delta arrow + label, single accent
PILL_BG           = "#3A3A3A"   # verdict pill background (single colour)

VERDICT_SHORT = {
    "supported":         "supp.",
    "weak_hybrid":       "weak",
    "not_supported":     "not supp.",
    "no_eligible_eQTL":  "no elig.",
}

# Δ_lead colour BANDS (kept in spec; current data lands in 250K+ band)
DELTA_50K  = 50_000
DELTA_250K = 250_000


# ---------------------------------------------------------------------------
# Cache slice loaders (find by anchor name prefix)
# ---------------------------------------------------------------------------
def find_cache_file(prefix: str) -> Path | None:
    candidates = sorted(CACHE.glob(f"{prefix}_chr*.tsv.gz"))
    return candidates[0] if candidates else None


def load_eqtl_slice(anchor: str, gene_id: str) -> pd.DataFrame:
    """Read cached eQTL slice for anchor; filter to gene_id rows."""
    f = find_cache_file(f"eqtl_{anchor}")
    if f is None:
        return pd.DataFrame()
    n_cols = pd.read_csv(f, sep="\t", header=None, nrows=1).shape[1]
    df = pd.read_csv(f, sep="\t", header=None,
                     names=EQTL_NAMES[:n_cols],
                     low_memory=False)
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["pvalue"]   = pd.to_numeric(df["pvalue"], errors="coerce")
    if "gene_id" in df.columns:
        gid_re = rf"^{re.escape(gene_id)}(?:\.|$)"
        mask = (df["gene_id"].astype(str).str.contains(gid_re, regex=True)
                | df["molecular_trait_id"].astype(str)
                                          .str.contains(gid_re, regex=True))
        df = df[mask].copy()
    df = df.dropna(subset=["position", "pvalue"])
    return df


def load_gwas_slice(anchor: str) -> pd.DataFrame:
    f = find_cache_file(f"gwas_{anchor}")
    if f is None:
        return pd.DataFrame()
    n_cols = pd.read_csv(f, sep="\t", header=None, nrows=1).shape[1]
    df = pd.read_csv(f, sep="\t", header=None,
                     names=GWAS_NAMES[:n_cols],
                     low_memory=False)
    df["base_pair_location"] = pd.to_numeric(df["base_pair_location"],
                                              errors="coerce")
    df["p_value"]            = pd.to_numeric(df["p_value"], errors="coerce")
    df = df.dropna(subset=["base_pair_location", "p_value"])
    return df


# ---------------------------------------------------------------------------
# Anchor coordinates (from h4 / h4d ANCHORS list, locked)
# ---------------------------------------------------------------------------
ANCHOR_COORD = {
    # v2 pre-registered
    "PRELP":  (1,  203475828, 203491352),
    "OMD":    (9,  92414245,  92424461),
    "IDH2":   (15, 90083045,  90102504),
    "NDUFA8": (9,  122144058, 122159819),
    # v3 LE extension
    "NTRK2":  (9,  84668522,  85029699),
    "NTF3":   (12, 5404888,   5495026),
    "FRS2":   (12, 69464533,  69513650),
    "BMP10":  (2,  68869054,  68886842),
    "SMAD9":  (13, 36879326,  36925087),
    "SMAD6":  (15, 66703041,  66782082),
    # PITX2 calibrator
    "PITX2":  (4,  110610001, 110633149),
}

GENE_ID_MAP = {
    "PRELP":  "ENSG00000188783", "OMD":    "ENSG00000127083",
    "IDH2":   "ENSG00000182054", "NDUFA8": "ENSG00000119421",
    "NTRK2":  "ENSG00000148053", "NTF3":   "ENSG00000185652",
    "FRS2":   "ENSG00000166225", "BMP10":  "ENSG00000163217",
    "SMAD9":  "ENSG00000120693", "SMAD6":  "ENSG00000137834",
    "PITX2":  "ENSG00000164093",
}

# 11 panels in this canonical order (4 v2 → 6 v3 → PITX2)
ANCHOR_ORDER = ["PRELP", "OMD", "IDH2", "NDUFA8",
                "NTRK2", "NTF3", "FRS2",
                "BMP10", "SMAD9", "SMAD6",
                "PITX2"]


# ---------------------------------------------------------------------------
# Per-panel rendering
# ---------------------------------------------------------------------------
def delta_color(delta_bp: float) -> str:
    """Δ-arrow colour band — kept for spec compliance; in current data
    100% of anchors land in the >=250 kb (brick-red) band.  All three
    bands resolve to the same single accent colour DELTA_BRICK so that
    the visual is monochrome; the threshold logic is preserved (and
    surfaced in the figure legend) so future runs that produce green/
    yellow Δs render with distinct colours without code change.
    """
    if pd.isna(delta_bp):
        return "#BDBDBD"
    return DELTA_BRICK   # single muted brick-red for the v3 dataset


def render_panel(ax, anchor: str, meta: dict) -> dict:
    """Render one mini locus-zoom; return summary stats for legend panel."""
    chrom, gene_start, gene_end = ANCHOR_COORD[anchor]
    gene_id = GENE_ID_MAP[anchor]
    tss = gene_start

    eqtl = load_eqtl_slice(anchor, gene_id)
    gwas = load_gwas_slice(anchor)

    # --- coordinates: convert positions to kb relative to TSS ---
    win_lo_bp = tss - 1_000_000
    win_hi_bp = tss + 1_000_000

    summary = {
        "anchor":       anchor,
        "delta_bp":     np.nan,
        "lead_gwas_kb": np.nan,
        "lead_eqtl_kb": np.nan,
        "lead_gwas_p":  np.nan,
        "lead_eqtl_p":  np.nan,
    }

    # GWAS
    g_kb = (gwas["base_pair_location"].to_numpy(float) - tss) / 1000.0 \
           if not gwas.empty else np.array([])
    g_log = -np.log10(np.clip(gwas["p_value"].to_numpy(float),
                               1e-300, 1.0)) if not gwas.empty else np.array([])
    if not gwas.empty:
        ig = int(np.argmin(gwas["p_value"].to_numpy(float)))
        lead_g_pos = float(gwas["base_pair_location"].iloc[ig])
        lead_g_p   = float(gwas["p_value"].iloc[ig])
        lead_g_kb  = (lead_g_pos - tss) / 1000.0
    else:
        lead_g_pos = np.nan; lead_g_p = np.nan; lead_g_kb = np.nan

    # eQTL
    e_kb = (eqtl["position"].to_numpy(float) - tss) / 1000.0 \
           if not eqtl.empty else np.array([])
    e_log = -np.log10(np.clip(eqtl["pvalue"].to_numpy(float),
                               1e-300, 1.0)) if not eqtl.empty else np.array([])
    if not eqtl.empty:
        ie = int(np.argmin(eqtl["pvalue"].to_numpy(float)))
        lead_e_pos = float(eqtl["position"].iloc[ie])
        lead_e_p   = float(eqtl["pvalue"].iloc[ie])
        lead_e_kb  = (lead_e_pos - tss) / 1000.0
    else:
        lead_e_pos = np.nan; lead_e_p = np.nan; lead_e_kb = np.nan

    # ---- shared y range (mirror) ----
    ymax = max(2.0,
               (np.nanmax(g_log) if g_log.size else 0),
               (np.nanmax(e_log) if e_log.size else 0))
    # Round up to nicer grid value
    ymax = float(np.ceil(ymax / 5.0) * 5.0)

    # ---- background card ----
    ax.set_xlim(-1000, 1000)               # ±1 Mb in kb
    ax.set_ylim(-ymax * 1.10, ymax * 1.10)
    ax.axhline(0, color="#888888", lw=0.45, zorder=2)
    ax.axvline(0, color="#888888", lw=0.4,  ls=(0, (1.5, 1.5)), zorder=2)

    # gene body strip on the centre line
    gene_lo_kb = (gene_start - tss) / 1000.0
    gene_hi_kb = (gene_end   - tss) / 1000.0
    ax.add_patch(Rectangle((gene_lo_kb, -ymax * 0.025),
                           max(gene_hi_kb - gene_lo_kb, 4.0),
                           ymax * 0.05,
                           facecolor="#1A1A1A", edgecolor="none", zorder=4))

    # ---- GWAS scatter (top, filled) ----
    if g_kb.size:
        ax.scatter(g_kb, g_log, s=3.0, c=GRAY_GWAS_PT,
                    alpha=0.45, edgecolors="none", zorder=3,
                    rasterized=True)
        # lead diamond — slightly larger, true black, anchors the eye
        ax.scatter([lead_g_kb], [-np.log10(max(lead_g_p, 1e-300))],
                    s=50, marker="D", facecolor=GRAY_GWAS_LEAD,
                    edgecolor="white", linewidth=0.6, zorder=6)
        ax.axvline(lead_g_kb, color=GRAY_GWAS_LEAD, lw=0.45,
                    ls=(0, (3, 2)), alpha=0.45, zorder=2)

    # ---- eQTL scatter (bottom, mid-grey filled, mirrored about y=0) ----
    if e_kb.size:
        ax.scatter(e_kb, -e_log, s=3.0, c=GRAY_EQTL_PT,
                    alpha=0.45, edgecolors="none", zorder=3,
                    rasterized=True)
        # lead: hollow white diamond with mid-grey edge
        ax.scatter([lead_e_kb], [np.log10(max(lead_e_p, 1e-300))],
                    s=50, marker="D", facecolor="white",
                    edgecolor=GRAY_EQTL_EDGE, linewidth=0.7, zorder=6)
        ax.axvline(lead_e_kb, color=GRAY_EQTL_EDGE, lw=0.45,
                    ls=(0, (1, 2)), alpha=0.55, zorder=2)

    # ---- Δ_lead arrow (single muted brick, thin) ----
    if not (np.isnan(lead_g_pos) or np.isnan(lead_e_pos)):
        delta_bp = abs(lead_g_pos - lead_e_pos)
        col = delta_color(delta_bp)
        y_arrow = ymax * 0.12
        ax.annotate("",
                    xy=(lead_g_kb, y_arrow),
                    xytext=(lead_e_kb, y_arrow),
                    arrowprops=dict(arrowstyle="<|-|>",
                                     mutation_scale=4.5,
                                     lw=0.6, color=col,
                                     shrinkA=0, shrinkB=0))
        mid_kb = (lead_g_kb + lead_e_kb) / 2
        delta_label = (f"\u0394 = {int(round(delta_bp/1000))} kb"
                       if delta_bp >= 1000
                       else f"\u0394 = {int(round(delta_bp))} bp")
        ax.text(mid_kb, y_arrow + ymax * 0.05, delta_label,
                ha="center", va="bottom", fontsize=4.0,
                color=col, style="italic", fontweight="normal")
        summary["delta_bp"]     = delta_bp
        summary["lead_gwas_kb"] = lead_g_kb
        summary["lead_eqtl_kb"] = lead_e_kb
        summary["lead_gwas_p"]  = lead_g_p
        summary["lead_eqtl_p"]  = lead_e_p

    # ---- title bar (anchor + family + role) ----
    # Two-tone: anchor name in 7pt bold black; family/role in 5pt grey.
    fam_lbl  = meta.get("family", "")
    role_lbl = meta.get("role", "")
    sub_lbl  = fam_lbl
    if role_lbl and role_lbl not in ("anchor", "calibrator"):
        sub_lbl = f"{fam_lbl} \u00b7 {role_lbl}"
    ax.set_title("", pad=2.0)  # clear default title slot
    # Custom two-tone title via figure-coords text just above axes.
    bbox = ax.get_position()
    fig_obj = ax.figure
    title_y = bbox.y1 + 0.012
    # Anchor name (bold black, left-aligned).  Sub-label (family/role)
    # offset by enough to clear up to 7-char anchor names (NDUFA8 = 6).
    fig_obj.text(bbox.x0 + 0.002, title_y, anchor,
                  fontsize=7.0, fontweight="bold", color=FG_COLOR,
                  ha="left", va="bottom")
    fig_obj.text(bbox.x0 + 0.002 + 0.055, title_y, sub_lbl,
                  fontsize=5.0, color="#666666",
                  ha="left", va="bottom")

    # ---- verdict pill (figure coords, OUTSIDE axes data area) ----
    # Single dark-grey pill at the right end of the title bar so it
    # never competes with scatter or Δ arrow inside the panel.
    verdict = str(meta.get("verdict", "no_eligible_eQTL"))
    v_short = VERDICT_SHORT.get(verdict, verdict)
    pill_w  = 0.060   # figure-coord width
    pill_h  = 0.018
    pill_x  = bbox.x1 - pill_w - 0.001
    pill_y  = title_y - 0.002
    fig_obj.patches.append(FancyBboxPatch(
        (pill_x, pill_y), pill_w, pill_h,
        boxstyle="round,pad=0.0,rounding_size=0.005",
        transform=fig_obj.transFigure,
        facecolor=PILL_BG, edgecolor="white", linewidth=0.4, zorder=10,
    ))
    fig_obj.text(pill_x + pill_w / 2, pill_y + pill_h / 2, v_short,
                  fontsize=4.6, color="white", fontweight="bold",
                  ha="center", va="center", zorder=11)

    # ---- eGene / Bonf inline next to TSS (replaces bottom-right badge) ----
    is_eg = bool(meta.get("is_eGene_atrial", False))
    bonf  = bool(meta.get("bonf_passed", False))
    if is_eg and bonf:
        eg_label = "eGene\u00b7Bonf"
    elif is_eg:
        eg_label = "eGene"
    elif bonf:
        eg_label = "Bonf"
    else:
        eg_label = "no eGene"
    eg_col = FG_COLOR if (is_eg or bonf) else "#999999"
    ax.text(0, ymax * 0.05, eg_label,
             ha="center", va="bottom", fontsize=4.0,
             color=eg_col, style="italic")

    # ---- y-axis MAX corner annotation (right corners) ----
    ax.text(990, ymax * 0.93, f"max={int(ymax)}",
             ha="right", va="top", fontsize=3.8,
             color="#888888", style="italic")
    ax.text(990, -ymax * 0.93, f"max={int(ymax)}",
             ha="right", va="bottom", fontsize=3.8,
             color="#888888", style="italic")

    # ---- y-axis ticks ----
    ax.set_yticks([-ymax, -ymax/2, 0, ymax/2, ymax])
    yticklabels = [f"{int(ymax)}", f"{int(ymax/2)}", "0",
                   f"{int(ymax/2)}", f"{int(ymax)}"]
    ax.set_yticklabels(yticklabels, fontsize=4.4, color="#888888")
    ax.tick_params(axis="y", length=1.5, pad=1.0,
                    colors="#888888")

    # ---- x-axis ticks (kb relative to TSS) ----
    ax.set_xticks([-1000, -500, 0, 500, 1000])
    ax.set_xticklabels(["-1 Mb", "-500", "0\nTSS", "500", "+1 Mb"],
                        fontsize=4.4, color="#888888")
    ax.tick_params(axis="x", length=1.5, pad=1.0,
                    colors="#888888")

    # ---- spines ----
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.4)
        ax.spines[s].set_color("#888888")

    # ---- track sub-labels (small italic, left edge) ----
    ax.text(-980, ymax * 0.92, "GWAS", ha="left", va="top",
             fontsize=4.0, color="#444444",
             style="italic")
    ax.text(-980, -ymax * 0.92, "atrial eQTL", ha="left", va="bottom",
             fontsize=4.0, color="#444444",
             style="italic")

    return summary


# ---------------------------------------------------------------------------
# Summary / legend panel (12th cell)
# ---------------------------------------------------------------------------
def render_summary_panel(ax, summaries: list[dict],
                          consol: pd.DataFrame) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # Title
    ax.text(0.0, 0.96, "ANCHOR LOCUS ARCHITECTURE",
            fontsize=6.4, fontweight="bold", color=FG_COLOR,
            ha="left", va="top")

    # Distance distribution
    deltas = pd.Series([s["delta_bp"] for s in summaries
                         if s["anchor"] != "PITX2"])
    n_anchor = len(deltas)
    n_lt50   = int((deltas < DELTA_50K).sum())
    n_lt250  = int(((deltas >= DELTA_50K) & (deltas < DELTA_250K)).sum())
    n_ge250  = int((deltas >= DELTA_250K).sum())

    rows_dist = [
        ("Δ < 50 kb (green)",
         f"{n_lt50:>2d} / {n_anchor}", "#1B7837"),
        ("50 ≤ Δ < 250 kb (yellow)",
         f"{n_lt250:>2d} / {n_anchor}", "#E08214"),
        ("Δ ≥ 250 kb  (red, PITX2-mode)",
         f"{n_ge250:>2d} / {n_anchor}", "#B2182B"),
    ]
    y0 = 0.86
    for lab, val, col in rows_dist:
        ax.text(0.00, y0, lab, ha="left", va="center",
                fontsize=4.8, color=col, fontweight="bold")
        ax.text(0.78, y0, val, ha="right", va="center",
                fontsize=4.8, color=col, fontweight="bold",
                family="monospace")
        y0 -= 0.060

    # Eligibility distribution (anchors only, excludes PITX2)
    anchors = consol.query("kind == 'anchor'")
    n_eg   = int(anchors["is_eGene_atrial"].sum())
    n_bonf = int(anchors["bonf_passed"].sum())
    n_5e8  = int(anchors["gate_5e8_legacy"].sum())
    n_agree = int(anchors["eligibility_agreement"].sum())

    y0 -= 0.02
    ax.text(0.0, y0, "Plan C eligibility (anchors only)",
            fontsize=5.2, fontweight="bold", color=FG_COLOR,
            ha="left", va="center")
    y0 -= 0.060
    egene_anchors = anchors.query("is_eGene_atrial")["anchor"].tolist()
    bonf_anchors  = anchors.query("bonf_passed")["anchor"].tolist()
    rows_elig = [
        ("PRIMARY (atrial eGene q<0.05)",
         f"{n_eg:>2d} / {len(anchors)}",
         f"{', '.join(egene_anchors) if egene_anchors else '—'}"),
        ("SENSITIVITY (per-gene Bonf)",
         f"{n_bonf:>2d} / {len(anchors)}",
         f"{', '.join(bonf_anchors) if bonf_anchors else '—'}"),
        ("Plan-C agreement",
         f"{n_agree:>2d} / {len(anchors)}",
         "all agree" if n_agree == len(anchors) else "(see caption)"),
        ("Legacy 5e-8 (NOT used)",
         f"{n_5e8:>2d} / {len(anchors)}",
         "transparency only"),
    ]
    for lab, val, note in rows_elig:
        ax.text(0.00, y0, lab, ha="left", va="center",
                fontsize=4.6, color=FG_COLOR)
        ax.text(0.55, y0, val, ha="right", va="center",
                fontsize=4.6, color=FG_COLOR,
                family="monospace", fontweight="bold")
        ax.text(0.58, y0, note, ha="left", va="center",
                fontsize=4.2, color=MUTED_COLOR, style="italic")
        y0 -= 0.045

    # Verdict legend (now monochrome — see header comment)
    y0 -= 0.020
    ax.text(0.0, y0, "Verdict pill (monochrome)",
            fontsize=5.0, fontweight="bold", color=FG_COLOR,
            ha="left", va="center")
    y0 -= 0.055
    n_supp = int((anchors["verdict"] == "supported").sum())
    n_weak = int((anchors["verdict"] == "weak_hybrid").sum())
    n_not  = int((anchors["verdict"] == "not_supported").sum())
    n_no   = int((anchors["verdict"] == "no_eligible_eQTL").sum())
    legend_rows = [
        ("supported",         n_supp),
        ("weak_hybrid",       n_weak),
        ("not_supported",     n_not),
        ("no_eligible_eQTL",  n_no),
    ]
    for lab, n in legend_rows:
        # All pills use the single PILL_BG colour (verdict variation
        # adds no information for this dataset; spec lock 2026-04-23).
        ax.add_patch(FancyBboxPatch(
            (0.00, y0 - 0.012), 0.045, 0.024,
            boxstyle="round,pad=0.0,rounding_size=0.005",
            facecolor=PILL_BG, edgecolor="white", linewidth=0.4))
        ax.text(0.058, y0, lab, fontsize=4.6, color=FG_COLOR,
                ha="left", va="center")
        ax.text(0.78, y0, f"{n:>2d} / {len(anchors)}",
                fontsize=4.6, color=FG_COLOR,
                ha="right", va="center", family="monospace")
        y0 -= 0.045

    # Footer note (updated for v3 visual refresh)
    ax.text(0.0, 0.04,
            "Top: AF GWAS \u2212log10(p)   Bottom: atrial cis-eQTL \u2212log10(p)\n"
            "Filled diamond = GWAS lead; hollow = eQTL lead.\n"
            "Brick \u0394 arrow = physical distance between leads.\n"
            "Verdict pill above each panel; eGene/Bonf badge at TSS.",
            fontsize=3.8, color="#666666", style="italic",
            ha="left", va="bottom", linespacing=1.4)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    set_nature_style()
    consol = pd.read_csv(IN_CONSOL, sep="\t",
                          keep_default_na=False, na_values=[""])
    for c in ("lead_eqtl_p", "n_eqtl", "coloc_PP4", "coloc_PP3",
              "bonf_threshold"):
        if c in consol.columns:
            consol[c] = pd.to_numeric(consol[c], errors="coerce")
    consol_idx = consol.set_index("anchor")

    # Build figure
    n_cols = 4
    n_rows = 3
    fig_w = COL_DOUBLE
    fig_h = 6.6
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(fig_w, fig_h))
    # Extra wspace/hspace + a bit more top so the figure-coord
    # title bar + verdict pill above each axes have room to render.
    fig.subplots_adjust(left=0.045, right=0.985, top=0.91, bottom=0.05,
                         wspace=0.25, hspace=0.42)

    summaries = []
    panel_idx = 0
    for r in range(n_rows):
        for c in range(n_cols):
            ax = axes[r, c]
            if panel_idx < len(ANCHOR_ORDER):
                anchor = ANCHOR_ORDER[panel_idx]
                meta = {
                    "family":           consol_idx.loc[anchor, "family"],
                    "role":             consol_idx.loc[anchor, "role"],
                    "verdict":          consol_idx.loc[anchor, "verdict"],
                    "is_eGene_atrial":  consol_idx.loc[anchor, "is_eGene_atrial"],
                    "bonf_passed":      consol_idx.loc[anchor, "bonf_passed"],
                }
                # bool conversion (TSV may have stringified)
                for k in ("is_eGene_atrial", "bonf_passed"):
                    v = meta[k]
                    meta[k] = (v is True) or (str(v).lower() == "true")
                summary = render_panel(ax, anchor, meta)
                summaries.append(summary)
                panel_idx += 1
            else:
                render_summary_panel(ax, summaries, consol)

    fig.suptitle(
        "Locus architecture: AF GWAS lead and atrial cis-eQTL lead are "
        "physically distinct across all anchors",
        fontsize=8.0, y=0.978, fontweight="bold",
    )
    fig.text(0.012, 0.97, "a",
             fontsize=10, fontweight="bold", ha="left", va="top",
             transform=fig.transFigure)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")

    # Also write the per-anchor summary so 6A figure note can quote it
    summary_df = pd.DataFrame(summaries)
    summary_tsv = RES / "fig6a_v3_locus_zoom_summary.tsv"
    summary_df.to_csv(summary_tsv, sep="\t", index=False)
    print(f"wrote {summary_tsv.relative_to(ROOT)}")
    # Print the distance distribution for the caption-injection step
    da = summary_df.query("anchor != 'PITX2'")["delta_bp"]
    print(f"Anchor Δ_lead distribution (n=10): "
          f"<50kb = {int((da<50_000).sum())}, "
          f"50-250kb = {int(((da>=50_000)&(da<250_000)).sum())}, "
          f"≥250kb = {int((da>=250_000).sum())}")


if __name__ == "__main__":
    main()
