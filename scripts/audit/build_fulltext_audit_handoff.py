#!/usr/bin/env python3
"""Build read-only fulltext audit handoff for whole-manuscript reconciliation.

Default output namespace (v1_1):
  results/handoff/fulltext_audit_v1_1/
  results/handoff/fulltext_audit_v1_1.zip(+.sha256)

Does not overwrite results/handoff/fulltext_audit_v1/.
Does not modify manuscript, figures, result tables, other scripts, routing,
or submission-pack assets. Rerunnable and deterministic for a fixed snapshot.
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from PIL import Image


def detect_root() -> Path:
    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        if (parent / ".git").exists():
            return parent
    raise SystemExit("FATAL: project root with .git not found")


ROOT = detect_root()
BUNDLE_NAME = "fulltext_audit_v1_1"
OUT = ROOT / "results" / "handoff" / BUNDLE_NAME
PREV = OUT / "figure_previews"
MSRC = OUT / "manuscript_sources"
AUTH = OUT / "07_authoritative_source_tables"
ZIP_PATH = ROOT / "results" / "handoff" / f"{BUNDLE_NAME}.zip"
ZIP_SHA = ROOT / "results" / "handoff" / f"{BUNDLE_NAME}.zip.sha256"
MAX_COPY_BYTES = 20 * 1024 * 1024
UTC = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
# Protect prior frozen handoff namespace
LEGACY_V1 = ROOT / "results" / "handoff" / "fulltext_audit_v1"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args: str) -> str:
    r = subprocess.run(["git", *args], cwd=str(ROOT), capture_output=True, text=True, check=False)
    return (r.stdout or "").strip()


def run_version(cmd: list[str]) -> str:
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, check=False)
        out = (r.stdout or r.stderr or "").strip().splitlines()
        return out[0] if out else "unavailable"
    except Exception:
        return "unavailable"


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def mtime_utc(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def kv_tsv(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t")
    if set(df.columns) >= {"metric", "value"}:
        return dict(zip(df["metric"].astype(str), df["value"].astype(str)))
    if set(df.columns) >= {"key", "value"}:
        return dict(zip(df["key"].astype(str), df["value"].astype(str)))
    raise ValueError(f"not a key/value TSV: {path}")


def tracked_set() -> set[str]:
    return set(git("ls-files").splitlines())


def status_map() -> dict[str, str]:
    m: dict[str, str] = {}
    for line in git("status", "--porcelain").splitlines():
        if len(line) < 4:
            continue
        st, rest = line[:2].strip(), line[3:]
        if " -> " in rest:
            rest = rest.split(" -> ", 1)[1]
        m[rest] = st
    return m


class ClaimRegistry:
    def __init__(self) -> None:
        self.rows: list[dict] = []
        self._by_id: dict[str, dict] = {}

    def add(self, **kw) -> None:
        cid = kw["claim_id"]
        if cid in self._by_id:
            prev = self._by_id[cid]
            if str(prev["value_text"]) != str(kw["value_text"]):
                raise SystemExit(
                    f"FATAL conflicting claim {cid}: "
                    f"{prev['value_text']!r} ({prev['authoritative_source_path']}) vs "
                    f"{kw['value_text']!r} ({kw['authoritative_source_path']})"
                )
            return
        row = {
            "claim_id": cid,
            "figure_id": kw["figure_id"],
            "panel_id": kw.get("panel_id", ""),
            "analysis_layer": kw.get("analysis_layer", ""),
            "claim_name": kw.get("claim_name", ""),
            "value_text": kw["value_text"],
            "value_numeric": kw.get("value_numeric", ""),
            "denominator": kw.get("denominator", ""),
            "unit": kw.get("unit", ""),
            "entity_ids": kw.get("entity_ids", ""),
            "definition": kw.get("definition", ""),
            "interpretation_boundary": kw.get("interpretation_boundary", ""),
            "authoritative_source_path": kw["authoritative_source_path"],
            "authoritative_source_sha256": kw["authoritative_source_sha256"],
            "qc_status": kw.get("qc_status", "PASS"),
        }
        self._by_id[cid] = row
        self.rows.append(row)

    def to_df(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)


class BundleFiles:
    def __init__(self, tracked: set[str], status: dict[str, str]) -> None:
        self.tracked = tracked
        self.status = status
        self.rows: list[dict] = []

    def record(
        self,
        bundle_path: Path,
        original: Path | None,
        role: str,
        figure_id: str = "",
        authoritative: bool = False,
        notes: str = "",
    ) -> None:
        orig_rel = rel(original) if original and original.is_file() else ""
        git_tracked = orig_rel in self.tracked if orig_rel else False
        git_st = self.status.get(
            orig_rel,
            "clean" if git_tracked else ("generated" if not orig_rel else "untracked"),
        )
        self.rows.append({
            "bundle_path": str(bundle_path.relative_to(OUT)),
            "original_project_path": orig_rel,
            "artifact_role": role,
            "figure_id": figure_id,
            "bytes": bundle_path.stat().st_size if bundle_path.is_file() else 0,
            "mtime_utc": mtime_utc(bundle_path) if bundle_path.is_file() else "",
            "sha256": sha256(bundle_path) if bundle_path.is_file() else "",
            "git_tracked": git_tracked,
            "git_status": git_st,
            "authoritative": authoritative,
            "notes": notes,
        })

    def record_skipped_large(self, src: Path, role: str, figure_id: str = "") -> None:
        orig_rel = rel(src)
        self.rows.append({
            "bundle_path": "",
            "original_project_path": orig_rel,
            "artifact_role": role + "_too_large_not_copied",
            "figure_id": figure_id,
            "bytes": src.stat().st_size,
            "mtime_utc": mtime_utc(src),
            "sha256": sha256(src),
            "git_tracked": orig_rel in self.tracked,
            "git_status": self.status.get(orig_rel, ""),
            "authoritative": True,
            "notes": f"exceeds {MAX_COPY_BYTES} bytes; path+sha only",
        })

    def to_df(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)


def copy_small(
    src: Path,
    dest: Path,
    bundle: BundleFiles,
    role: str,
    figure_id: str = "",
    authoritative: bool = True,
    notes: str = "",
) -> bool:
    if not src.is_file():
        return False
    if src.stat().st_size > MAX_COPY_BYTES:
        bundle.record_skipped_large(src, role, figure_id)
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    bundle.record(dest, src, role, figure_id, authoritative, notes)
    return True


def preview_png(src_png: Path, dest: Path, bundle: BundleFiles, figure_id: str, max_w: int = 1100) -> None:
    if not src_png.is_file():
        return
    im = Image.open(src_png).convert("RGB")
    w, h = im.size
    if w > max_w:
        im = im.resize((max_w, int(h * max_w / w)), Image.Resampling.LANCZOS)
    dest.parent.mkdir(parents=True, exist_ok=True)
    im.save(dest, format="PNG", optimize=True)
    bundle.record(dest, src_png, "figure_preview", figure_id, False, "reduced preview")


FIGURES = {
    "1": {
        "title": "Figure 1 | Healthy atrial chamber-identity NCP definition",
        "panels": "A–D",
        "pdf": ROOT / "figures/figure1_composite.pdf",
        "svg": ROOT / "figures/figure1_composite.svg",
        "png": ROOT / "figures/figure1_composite.png",
        "qc": ROOT / "results/figure1_composite_qc.md",
        "renderer": "scripts/assemble_figure1.py",
        "sources": [
            ROOT / "results/ncp_v3_singscore_reanalysis/fig1b_healthy_reference_direction_magnitude_summary.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig1b_healthy_reference_concordance_summary.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/frozen_program_metadata.tsv",
            ROOT / "results/ncp_v3_family_labels.tsv",
            ROOT / "results/ncp_v3_canonical_display_order.tsv",
            ROOT / "results/innovation_ncp_v3_summary.tsv",
            ROOT / "results/fig1c_selection_provenance_qc.md",
        ],
    },
    "2": {
        "title": "Figure 2 | Bulk disease deployment of the frozen NCP",
        "panels": "A–D",
        "pdf": ROOT / "figures/figure2_composite.pdf",
        "svg": ROOT / "figures/figure2_composite.svg",
        "png": ROOT / "figures/figure2_composite.png",
        "qc": ROOT / "results/figure2_composite_qc.md",
        "renderer": "scripts/ncp_v3_singscore_reanalysis/assemble_figure2_composite.py",
        "sources": [
            ROOT / "results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_summary.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_summary.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_summary.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv",
        ],
    },
    "3": {
        "title": "Figure 3 | LA single-nucleus localization / carrier fidelity",
        "panels": "A–D",
        "pdf": ROOT / "figures/ncp_v3_singscore_reanalysis/figure3_la_single_nucleus_localization.pdf",
        "svg": ROOT / "figures/ncp_v3_singscore_reanalysis/figure3_la_single_nucleus_localization.svg",
        "png": ROOT / "figures/ncp_v3_singscore_reanalysis/figure3_la_single_nucleus_localization.png",
        "qc": ROOT / "results/ncp_v3_singscore_reanalysis/figure3_composite_qc.md",
        "renderer": "scripts/ncp_v3_singscore_reanalysis/assemble_figure3_la_single_nucleus_localization.py",
        "sources": [
            ROOT / "results/fig3c_family_fidelity_classes.tsv",
            ROOT / "results/fig3_v3_carrier_fidelity_class.tsv",
            ROOT / "results/fig3_v3_family_celltype_delta.tsv",
            ROOT / "results/fig3a_localization_bridge_source.tsv",
        ],
    },
    "4": {
        "title": "Figure 4 | Cross-sectional duration / chronicity (GSE41177)",
        "panels": "A–D",
        "pdf": ROOT / "figures/ncp_v3_singscore_reanalysis/figure4_cross_sectional_duration.pdf",
        "svg": ROOT / "figures/ncp_v3_singscore_reanalysis/figure4_cross_sectional_duration.svg",
        "png": ROOT / "figures/ncp_v3_singscore_reanalysis/figure4_cross_sectional_duration.png",
        "qc": ROOT / "results/ncp_v3_singscore_reanalysis/figure4_composite_qc.md",
        "renderer": "scripts/ncp_v3_singscore_reanalysis/assemble_figure4_chronicity_composite.py",
        "sources": [
            ROOT / "results/ncp_v3_singscore_reanalysis/fig4a_duration_design_summary.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig4c_chronicity_class_summary_summary.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig4_lapv_sensitivity.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig4c_chronicity_class_summary_qc.md",
        ],
    },
    "5": {
        "title": "Figure 5 | Cross-layer triangulation (Path A)",
        "panels": "A–D",
        "pdf": ROOT / "figures/ncp_v3_singscore_reanalysis/figure5_cross_layer_triangulation.pdf",
        "svg": ROOT / "figures/ncp_v3_singscore_reanalysis/figure5_cross_layer_triangulation.svg",
        "png": ROOT / "figures/ncp_v3_singscore_reanalysis/figure5_cross_layer_triangulation.png",
        "qc": ROOT / "results/ncp_v3_singscore_reanalysis/figure5_composite_qc.md",
        "renderer": "scripts/ncp_v3_singscore_reanalysis/assemble_figure5_cross_layer_triangulation.py",
        "composite_manifest": ROOT / "results/ncp_v3_singscore_reanalysis/figure5_composite_manifest.tsv",
        "sources": [
            ROOT / "results/ncp_v3_singscore_reanalysis/figure5_triangulation_canonical.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig5b_cross_layer_triangulation_summary.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig5b_cross_layer_triangulation_source.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig5c_pathA_gate_architecture_summary.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/fig5d_geomx_roi_spatial_evidence_summary.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/figure5_geomx_lineage_audit.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/figure5_authoritative_input_manifest.tsv",
            ROOT / "results/ncp_v3_singscore_reanalysis/figure5_composite_manifest.tsv",
        ],
    },
    "6": {
        "title": "Figure 6 | Genetic traceability and colocalisation calibration",
        "panels": "A–D",
        "pdf": ROOT / "figures/figure6_genetic_traceability.pdf",
        "svg": ROOT / "figures/figure6_genetic_traceability.svg",
        "png": ROOT / "figures/figure6_genetic_traceability.png",
        "qc": ROOT / "results/figure6_composite_qc.md",
        "renderer": "scripts/assemble_figure6_genetic_traceability.py",
        "sources": [
            ROOT / "results/fig6a_anchor_coloc_calibration_source.tsv",
            ROOT / "results/fig6b_posterior_lead_distance_source.tsv",
            ROOT / "results/fig6c_family_af_overlap_source.tsv",
            ROOT / "results/fig6c_family_af_overlap_summary.tsv",
            ROOT / "results/fig6d_non_nested_traceability_source.tsv",
            ROOT / "results/fig6d_non_nested_traceability_summary.tsv",
            ROOT / "results/provenance/figure6_external_gene_sets_manifest.tsv",
            ROOT / "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv",
            ROOT / "results/figure6_composite_geometry.tsv",
        ],
    },
}


def qc_pass(path: Path) -> str:
    if not path.is_file():
        return "MISSING_QC"
    text = path.read_text(encoding="utf-8", errors="ignore")
    if re.search(r"Status:\s*\*\*PASS\*\*|assertion_status\s+PASS|QC:\s*PASS|\*\*PASS\*\*", text):
        return "PASS"
    if "PASS" in text and "FAIL" not in text:
        return "PASS"
    return "UNKNOWN"


def canvas_mm(svg: Path) -> tuple[str, str]:
    if not svg.is_file():
        return "", ""
    t = svg.read_text(encoding="utf-8", errors="ignore")
    w = re.search(r'width="([0-9.]+)mm"', t)
    h = re.search(r'height="([0-9.]+)mm"', t)
    return (w.group(1) if w else ""), (h.group(1) if h else "")


def build_claims(claims: ClaimRegistry) -> None:
    # ----- FIGURE 1 -----
    s1 = ROOT / "results/ncp_v3_singscore_reanalysis/fig1b_healthy_reference_direction_magnitude_summary.tsv"
    c1 = ROOT / "results/ncp_v3_singscore_reanalysis/fig1b_healthy_reference_concordance_summary.tsv"
    r = pd.read_csv(s1, sep="\t").iloc[0]
    rc = pd.read_csv(c1, sep="\t").iloc[0]
    h1, hc1 = sha256(s1), sha256(c1)
    if int(r["n_programs"]) != 30 or int(rc["n_programs"]) != 30:
        raise SystemExit("FATAL: Fig1 n_programs conflict")
    if int(r["n_la_biased"]) != int(rc["n_la"]) or int(r["n_ra_biased"]) != int(rc["n_ra"]):
        raise SystemExit("FATAL: Fig1 LA/RA counts disagree")
    if int(r["n_programs_8_of_8_concordant"]) != int(rc["n_full_concordance"]):
        raise SystemExit("FATAL: Fig1 concordance disagree")
    claims.add(claim_id="F1_N_PROGRAMS", figure_id="1", panel_id="B", analysis_layer="healthy_reference",
               claim_name="n_frozen_programs", value_text="30", value_numeric=30, unit="programs",
               entity_ids="frozen_panel", definition="Number of frozen programs",
               interpretation_boundary="Panel size lock",
               authoritative_source_path=rel(s1), authoritative_source_sha256=h1)
    claims.add(claim_id="F1_LA_BIASED", figure_id="1", panel_id="B", analysis_layer="healthy_reference",
               claim_name="n_la_biased", value_text="22", value_numeric=22, denominator=30, unit="programs",
               entity_ids="LA", definition="LA-biased programs", interpretation_boundary="Healthy composition",
               authoritative_source_path=rel(s1), authoritative_source_sha256=h1)
    claims.add(claim_id="F1_RA_BIASED", figure_id="1", panel_id="B", analysis_layer="healthy_reference",
               claim_name="n_ra_biased", value_text="8", value_numeric=8, denominator=30, unit="programs",
               entity_ids="RA", definition="RA-biased programs", interpretation_boundary="Healthy composition",
               authoritative_source_path=rel(s1), authoritative_source_sha256=h1)
    claims.add(claim_id="F1_CONCORDANT_8OF8", figure_id="1", panel_id="B", analysis_layer="healthy_reference",
               claim_name="donor_concordance_8of8",
               value_text="30/30 programs concordant in all 8 healthy donor pairs",
               value_numeric=30, denominator=30, unit="programs", entity_ids="8_healthy_donors",
               definition="Full 8/8 healthy-donor concordance", interpretation_boundary="Healthy concordance gate",
               authoritative_source_path=rel(s1), authoritative_source_sha256=h1)
    claims.add(claim_id="F1_MEDIAN_M", figure_id="1", panel_id="B", analysis_layer="healthy_reference",
               claim_name="median_M", value_text=str(r["median_M"]), value_numeric=float(r["median_M"]),
               unit="M", entity_ids="30_programs", definition="Median healthy program magnitude M",
               interpretation_boundary="Healthy-reference magnitude",
               authoritative_source_path=rel(s1), authoritative_source_sha256=h1)
    claims.add(claim_id="F1_PROGRAMS_TO_FAMILIES", figure_id="1", panel_id="D", analysis_layer="topology",
               claim_name="programs_to_families", value_text="30 programs mapped to 16 frozen families",
               value_numeric=16, denominator=30, unit="families", entity_ids="NCP_topology",
               definition="Fixed program→family topology", interpretation_boundary="Structural mapping only",
               authoritative_source_path=rel(c1), authoritative_source_sha256=hc1)
    sel = ROOT / "results/fig1c_selection_provenance_qc.md"
    txt = sel.read_text(encoding="utf-8")
    m_strict = re.search(r"strict-selected programs:\s*(\d+)", txt)
    m_fallback = re.search(r"fallback-selected programs:\s*(\d+)", txt)
    if not (m_strict and m_fallback):
        raise SystemExit("FATAL: cannot parse strict/fallback from fig1c_selection_provenance_qc.md")
    if (int(m_strict.group(1)), int(m_fallback.group(1))) != (25, 5):
        raise SystemExit(f"FATAL: expected 25/5 strict/fallback, got {m_strict.group(1)}/{m_fallback.group(1)}")
    claims.add(claim_id="F1_STRICT", figure_id="1", panel_id="C", analysis_layer="selection",
               claim_name="strict_selected", value_text="25", value_numeric=25, denominator=30, unit="programs",
               entity_ids="strict", definition="Strict-selected programs", interpretation_boundary="Selection provenance",
               authoritative_source_path=rel(sel), authoritative_source_sha256=sha256(sel))
    claims.add(claim_id="F1_FALLBACK", figure_id="1", panel_id="C", analysis_layer="selection",
               claim_name="fallback_selected", value_text="5", value_numeric=5, denominator=30, unit="programs",
               entity_ids="fallback", definition="Fallback-selected programs", interpretation_boundary="Selection provenance",
               authoritative_source_path=rel(sel), authoritative_source_sha256=sha256(sel))

    # ----- FIGURE 2 -----
    f2a = ROOT / "results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_summary.tsv"
    df2a = pd.read_csv(f2a, sep="\t")
    tot = df2a.loc[df2a["scope"].astype(str) == "total"].iloc[0]
    h2a = sha256(f2a)
    if (int(tot["paired_donors"]), int(tot["SR_donors"]), int(tot["AF_donors"])) != (51, 26, 25):
        raise SystemExit(f"FATAL: Fig2 donor totals unexpected: {dict(tot)}")
    claims.add(claim_id="F2_PAIRED_DONORS", figure_id="2", panel_id="A", analysis_layer="bulk_deployment",
               claim_name="paired_donors_total", value_text="51", value_numeric=51, unit="donors",
               entity_ids="paired_cohort", definition="Paired donors total (scope=total)",
               interpretation_boundary="Design N",
               authoritative_source_path=rel(f2a), authoritative_source_sha256=h2a)
    claims.add(claim_id="F2_SR_DONORS", figure_id="2", panel_id="A", analysis_layer="bulk_deployment",
               claim_name="SR_donors", value_text="26", value_numeric=26, unit="donors", entity_ids="SR",
               definition="SR paired donors", interpretation_boundary="Design N",
               authoritative_source_path=rel(f2a), authoritative_source_sha256=h2a)
    claims.add(claim_id="F2_AF_DONORS", figure_id="2", panel_id="A", analysis_layer="bulk_deployment",
               claim_name="AF_donors", value_text="25", value_numeric=25, unit="donors", entity_ids="AF",
               definition="AF paired donors", interpretation_boundary="Design N",
               authoritative_source_path=rel(f2a), authoritative_source_sha256=h2a)
    f2c = ROOT / "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_summary.tsv"
    f2cs = ROOT / "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv"
    s = pd.read_csv(f2c, sep="\t")
    src = pd.read_csv(f2cs, sep="\t")
    total_calls = int(s["n_calls"].sum())
    if total_calls != 48 or len(src) != 48:
        raise SystemExit(f"FATAL: Fig2C calls={total_calls} rows={len(src)}; expected 48")
    parts = [f"{row.state} {int(row.n_calls)}" for row in s.itertuples(index=False)]
    claims.add(claim_id="F2_STATE_COUNTS", figure_id="2", panel_id="C", analysis_layer="family_state",
               claim_name="corrected_state_counts", value_text=" / ".join(parts),
               value_numeric=48, denominator=48, unit="family_x_cohort_calls",
               entity_ids="Amplified,Preserved,Attenuated,Reversed",
               definition="Corrected family×cohort state call counts (sum to 48)",
               interpretation_boundary="Do not use legacy Figure 2 state counts",
               authoritative_source_path=rel(f2c), authoritative_source_sha256=sha256(f2c))
    f2d = ROOT / "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_summary.tsv"
    d = pd.read_csv(f2d, sep="\t")
    vc = d["consistency_class"].value_counts()
    high, mod, low = int(vc.get("High", 0)), int(vc.get("Moderate", 0)), int(vc.get("Low", 0))
    if (high, mod, low) != (6, 8, 2):
        raise SystemExit(f"FATAL: Fig2D classes {high}/{mod}/{low} != 6/8/2")
    claims.add(claim_id="F2_CONSISTENCY", figure_id="2", panel_id="D", analysis_layer="consistency",
               claim_name="High_Moderate_Low", value_text="High 6, Moderate 8, Low 2",
               denominator=16, unit="families", entity_ids="16_families",
               definition="Bulk cross-cohort consistency classes", interpretation_boundary="Descriptive tiers",
               authoritative_source_path=rel(f2d), authoritative_source_sha256=sha256(f2d))
    f2ds = ROOT / "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv"
    ds = pd.read_csv(f2ds, sep="\t")
    models = ds.loc[ds["estimate_type"].astype(str) == "pooled", "pooling_model"].dropna().astype(str).unique().tolist()
    if models != ["fixed-effect_inverse_variance"]:
        raise SystemExit(f"FATAL: unexpected pooling models {models}")
    claims.add(claim_id="F2_POOLING", figure_id="2", panel_id="D", analysis_layer="consistency",
               claim_name="pooling_model",
               value_text="fixed-effect inverse-variance weighting",
               entity_ids="family_effects", definition="Primary family pooling model",
               interpretation_boundary="Fixed-effect IVW; tau2 not estimated",
               authoritative_source_path=rel(f2ds), authoritative_source_sha256=sha256(f2ds))

    # ----- FIGURE 3 -----
    f3 = ROOT / "results/fig3c_family_fidelity_classes.tsv"
    f3b = ROOT / "results/fig3_v3_carrier_fidelity_class.tsv"
    df3 = pd.read_csv(f3, sep="\t")
    df3b = pd.read_csv(f3b, sep="\t")
    h3 = sha256(f3)

    def fams(frame, label):
        return sorted(
            frame.loc[frame["fidelity_class"].astype(str).str.contains(label, case=False, na=False), "family_id"]
            .astype(int).unique().tolist()
        )

    for label in ("Fidelity", "Relay", "Extension"):
        if fams(df3, label) != fams(df3b, label):
            raise SystemExit(f"FATAL: Fig3 {label} class disagreement")
    if len(fams(df3, "Below")) != len(fams(df3b, "Below")):
        raise SystemExit("FATAL: Fig3 Below-detection count disagreement")
    fid, relay, ext, bel = fams(df3, "Fidelity"), fams(df3, "Relay"), fams(df3, "Extension"), fams(df3, "Below")
    if fid != [10] or relay != [11, 15] or ext != [16] or len(bel) != 12:
        raise SystemExit(f"FATAL: Fig3 unexpected classes {fid}/{relay}/{ext}/{len(bel)}")
    claims.add(claim_id="F3_FIDELITY", figure_id="3", panel_id="C", analysis_layer="carrier_fidelity",
               claim_name="Fidelity", value_text="F10", value_numeric=1, denominator=16, unit="families",
               entity_ids="F10", definition="Fidelity-class families", interpretation_boundary="LA snRNA localization",
               authoritative_source_path=rel(f3), authoritative_source_sha256=h3)
    claims.add(claim_id="F3_RELAY", figure_id="3", panel_id="C", analysis_layer="carrier_fidelity",
               claim_name="Relay", value_text="F11, F15", value_numeric=2, denominator=16, unit="families",
               entity_ids="F11,F15", definition="Relay-class families", interpretation_boundary="LA snRNA localization",
               authoritative_source_path=rel(f3), authoritative_source_sha256=h3)
    claims.add(claim_id="F3_EXTENSION", figure_id="3", panel_id="C", analysis_layer="carrier_fidelity",
               claim_name="Extension", value_text="F16", value_numeric=1, denominator=16, unit="families",
               entity_ids="F16", definition="Extension-class families", interpretation_boundary="LA snRNA localization",
               authoritative_source_path=rel(f3), authoritative_source_sha256=h3)
    claims.add(claim_id="F3_BELOW", figure_id="3", panel_id="C", analysis_layer="carrier_fidelity",
               claim_name="Below_detection", value_text="12", value_numeric=12, denominator=16, unit="families",
               entity_ids=",".join(f"F{i}" for i in bel), definition="Below-detection families",
               interpretation_boundary="Not biological exclusion",
               authoritative_source_path=rel(f3), authoritative_source_sha256=h3)
    row11 = df3.loc[df3["family_id"] == 11].iloc[0]
    row15 = df3.loc[df3["family_id"] == 15].iloc[0]
    row16 = df3.loc[df3["family_id"] == 16].iloc[0]
    if str(row11["snRNA_supported_compartment"]).lower() != "adipocyte":
        raise SystemExit("FATAL: F11 compartment")
    if str(row15["snRNA_supported_compartment"]).lower() != "endothelial":
        raise SystemExit("FATAL: F15 compartment")
    if str(row16["snRNA_supported_compartment"]) != "atrialCM,fibroblast":
        raise SystemExit("FATAL: F16 compartments")
    claims.add(claim_id="F3_F11_COMPARTMENT", figure_id="3", panel_id="C", analysis_layer="carrier_fidelity",
               claim_name="F11_supported_compartment", value_text="Adipocyte", entity_ids="F11",
               definition="F11 supported snRNA compartment", interpretation_boundary="LA-only localization",
               authoritative_source_path=rel(f3), authoritative_source_sha256=h3)
    claims.add(claim_id="F3_F15_COMPARTMENT", figure_id="3", panel_id="C", analysis_layer="carrier_fidelity",
               claim_name="F15_supported_compartment", value_text="Endothelial only", entity_ids="F15",
               definition="F15 supported snRNA compartment", interpretation_boundary="LA-only localization",
               authoritative_source_path=rel(f3), authoritative_source_sha256=h3)
    claims.add(claim_id="F3_F16_COMPARTMENTS", figure_id="3", panel_id="C", analysis_layer="carrier_fidelity",
               claim_name="F16_supported_compartments", value_text="Fibroblast and aCM",
               entity_ids="F16", definition="F16 supported snRNA compartments (atrialCM,fibroblast)",
               interpretation_boundary="LA-only localization",
               authoritative_source_path=rel(f3), authoritative_source_sha256=h3)
    lock = ROOT / "results/fig3_publication_authority_lock.md"
    claims.add(claim_id="F3_LA_ONLY_BOUNDARY", figure_id="3", panel_id="A/C", analysis_layer="carrier_fidelity",
               claim_name="LA_only_localization_boundary", value_text="LA-only localization boundary",
               entity_ids="LA", definition="Localization restricted to LA single-nucleus context",
               interpretation_boundary="Not a bi-chamber discovery screen",
               authoritative_source_path=rel(lock), authoritative_source_sha256=sha256(lock))
    claims.add(claim_id="F3_PERM_NOT_FDR", figure_id="3", panel_id="C", analysis_layer="carrier_fidelity",
               claim_name="nominal_perm_not_FDR_discovery",
               value_text="nominal permutation localization is not FDR-controlled discovery",
               entity_ids="p_perm", definition="Permutation p for class evidence; not BH-FDR discovery",
               interpretation_boundary="Do not interpret as genome-wide FDR hits",
               authoritative_source_path=rel(f3), authoritative_source_sha256=h3)

    # ----- FIGURE 4 -----
    f4a = ROOT / "results/ncp_v3_singscore_reanalysis/fig4a_duration_design_summary.tsv"
    r4 = pd.read_csv(f4a, sep="\t").iloc[0]
    h4 = sha256(f4a)
    if int(r4["primary_eligible_af_laa"]) != 12:
        raise SystemExit("FATAL: Fig4 primary n != 12")
    claims.add(claim_id="F4_PRIMARY_N", figure_id="4", panel_id="A", analysis_layer="duration_design",
               claim_name="primary_AF_LAA_n", value_text="12", value_numeric=12, unit="donors",
               entity_ids="AF-LAA", definition="Primary-eligible AF-LAA donors",
               interpretation_boundary="SR not used for ρ/class",
               authoritative_source_path=rel(f4a), authoritative_source_sha256=h4)
    claims.add(claim_id="F4_LAPV_SENSITIVITY", figure_id="4", panel_id="A/B", analysis_layer="duration_design",
               claim_name="LAPV_regional_sensitivity_only", value_text=str(r4["lapv_role"]),
               entity_ids="LAPV", definition="LAPV role in duration design",
               interpretation_boundary="Sensitivity only; does not redefine primary classes",
               authoritative_source_path=rel(f4a), authoritative_source_sha256=h4)
    claims.add(claim_id="F4_NO_SR_ANCHORS", figure_id="4", panel_id="A", analysis_layer="duration_design",
               claim_name="SR_not_duration_anchors",
               value_text="SR samples are not used as duration anchors in the final main figure",
               entity_ids="SR", definition="SR excluded from correlation/class assignment",
               interpretation_boundary=f"sr_laa_anchors design field={int(r4['sr_laa_anchors'])} (non-analytic)",
               authoritative_source_path=rel(f4a), authoritative_source_sha256=h4)
    f4c = ROOT / "results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv"
    d4 = pd.read_csv(f4c, sep="\t")
    h4c = sha256(f4c)

    def cmem(label):
        return sorted(d4.loc[d4["chronicity_class"].astype(str).str.contains(label, case=False), "family_id"].astype(int).tolist())

    strong, weak, stable = cmem("Strong"), cmem("Weak"), cmem("Stable")
    if (strong, weak, stable) != ([2, 9], [3, 7, 8, 10, 13, 15], [1, 4, 5, 6, 11, 12, 14, 16]):
        raise SystemExit(f"FATAL: Fig4 class membership unexpected {strong}/{weak}/{stable}")
    claims.add(claim_id="F4_STRONG", figure_id="4", panel_id="C", analysis_layer="chronicity",
               claim_name="Strong_trend", value_text="F2, F9", value_numeric=2, denominator=16, unit="families",
               entity_ids="F2,F9", definition="Strong-trend families",
               interpretation_boundary="Cross-sectional AF-LAA duration",
               authoritative_source_path=rel(f4c), authoritative_source_sha256=h4c)
    claims.add(claim_id="F4_WEAK", figure_id="4", panel_id="C", analysis_layer="chronicity",
               claim_name="Weak_trend", value_text="F3, F7, F8, F10, F13, F15", value_numeric=6, denominator=16,
               unit="families", entity_ids="F3,F7,F8,F10,F13,F15", definition="Weak-trend families",
               interpretation_boundary="Cross-sectional AF-LAA duration",
               authoritative_source_path=rel(f4c), authoritative_source_sha256=h4c)
    claims.add(claim_id="F4_STABLE", figure_id="4", panel_id="C", analysis_layer="chronicity",
               claim_name="Stable", value_text="F1, F4, F5, F6, F11, F12, F14, F16", value_numeric=8,
               denominator=16, unit="families", entity_ids="F1,F4,F5,F6,F11,F12,F14,F16",
               definition="Stable families", interpretation_boundary="Cross-sectional AF-LAA duration",
               authoritative_source_path=rel(f4c), authoritative_source_sha256=h4c)
    claims.add(claim_id="F4_CLASS_COUNTS", figure_id="4", panel_id="C", analysis_layer="chronicity",
               claim_name="Strong_Weak_Stable_counts", value_text="2 / 6 / 8", denominator=16, unit="families",
               entity_ids="16_families", definition="Chronicity class counts",
               interpretation_boundary="Uses full-precision abs(rho_LAA)",
               authoritative_source_path=rel(f4c), authoritative_source_sha256=h4c)
    f4qc = ROOT / "results/ncp_v3_singscore_reanalysis/fig4c_chronicity_class_summary_qc.md"
    claims.add(claim_id="F4_FULL_PRECISION_RHO", figure_id="4", panel_id="C", analysis_layer="chronicity",
               claim_name="class_uses_full_precision_abs_rho_LAA",
               value_text="class uses full-precision abs(rho_LAA)", entity_ids="rho_LAA",
               definition="Chronicity classes from full-precision |rho_LAA|",
               interpretation_boundary="Display rounding does not redefine classes",
               authoritative_source_path=rel(f4qc), authoritative_source_sha256=sha256(f4qc))

    # ----- FIGURE 5 -----
    f5canon = ROOT / "results/ncp_v3_singscore_reanalysis/figure5_triangulation_canonical.tsv"
    f5s = ROOT / "results/ncp_v3_singscore_reanalysis/fig5b_cross_layer_triangulation_summary.tsv"
    f5a = ROOT / "results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv"
    f5c = ROOT / "results/ncp_v3_singscore_reanalysis/fig5c_pathA_gate_architecture_summary.tsv"
    r5, ra, rc5 = kv_tsv(f5s), kv_tsv(f5a), kv_tsv(f5c)
    if r5["Path_A_families"] != ra["Path_A_families"] or r5["Path_A_families"] != rc5["Path_A_families"]:
        raise SystemExit("FATAL: Path A families disagree across Fig5 sources")
    if r5["Path_A_families"] != "F10,F15":
        raise SystemExit(f"FATAL: Path A != F10,F15: {r5['Path_A_families']}")
    h5 = sha256(f5s)
    canon = pd.read_csv(f5canon, sep="\t")
    h5c = sha256(f5canon)
    if "support_count" not in canon.columns:
        raise SystemExit("FATAL: figure5_triangulation_canonical.tsv missing support_count")

    def fams_at(sc: int) -> list[int]:
        return sorted(canon.loc[canon["support_count"].astype(int) == sc, "family_id"].astype(int).tolist())

    fams_3, fams_2, fams_1, fams_0 = fams_at(3), fams_at(2), fams_at(1), fams_at(0)
    expected_1 = [1, 2, 3, 7, 8, 11, 12, 13, 14, 16]
    if fams_1 != expected_1:
        raise SystemExit(f"FATAL: support_count==1 set {fams_1} != expected {expected_1}")
    if (len(fams_3), len(fams_2), len(fams_1), len(fams_0)) != (2, 1, 10, 3):
        raise SystemExit(
            f"FATAL: support counts n={len(fams_3)}/{len(fams_2)}/{len(fams_1)}/{len(fams_0)} != 2/1/10/3"
        )
    # Cross-check summary metric names (internal 'silent' may remain in frozen summary)
    if int(float(r5["triple"])) != 2 or int(float(r5["double"])) != 1:
        raise SystemExit("FATAL: summary triple/double disagree with canonical support_count")
    if int(float(r5["single"])) != 10 or int(float(r5["silent"])) != 3:
        raise SystemExit("FATAL: summary single/silent disagree with canonical support_count")

    claims.add(claim_id="F5_PATH_A", figure_id="5", panel_id="B", analysis_layer="Path_A",
               claim_name="Path_A_families", value_text="F10, F15", value_numeric=2, denominator=16,
               unit="families", entity_ids="F10,F15", definition="Path A triple-positive families",
               interpretation_boundary="Triangulation support; genetics/GeoMx excluded",
               authoritative_source_path=rel(f5s), authoritative_source_sha256=h5)
    claims.add(claim_id="F5_3OF3", figure_id="5", panel_id="B", analysis_layer="triangulation",
               claim_name="support_3of3", value_text="2", value_numeric=2, denominator=16,
               unit="families", entity_ids=",".join(f"F{i}" for i in fams_3),
               definition="Families with 3/3 gates", interpretation_boundary="Support-count distribution",
               authoritative_source_path=rel(f5canon), authoritative_source_sha256=h5c)
    claims.add(claim_id="F5_2OF3_N", figure_id="5", panel_id="B", analysis_layer="triangulation",
               claim_name="support_2of3_n", value_text="1", value_numeric=1, denominator=16,
               unit="families", entity_ids=",".join(f"F{i}" for i in fams_2),
               definition="Count of families with 2/3 gates", interpretation_boundary="Support-count distribution",
               authoritative_source_path=rel(f5canon), authoritative_source_sha256=h5c)
    # F5_1OF3: membership from canonical support_count == 1 (never fall back to double_families)
    entity_1 = ",".join(f"F{i}" for i in fams_1)
    claims.add(claim_id="F5_1OF3", figure_id="5", panel_id="B", analysis_layer="triangulation",
               claim_name="support_1of3", value_text="10", value_numeric=10, denominator=16,
               unit="families", entity_ids=entity_1,
               definition="Families with exactly 1/3 gates (support_count == 1)",
               interpretation_boundary="Support-count distribution",
               authoritative_source_path=rel(f5canon), authoritative_source_sha256=h5c)
    zero_def = "Families satisfying zero of the three Path A gates"
    zero_bound = (
        "No predefined Path A gate was satisfied; this is not biological exclusion."
    )
    claims.add(claim_id="F5_0OF3_N", figure_id="5", panel_id="B", analysis_layer="triangulation",
               claim_name="support_0of3_n", value_text="3", value_numeric=3, denominator=16,
               unit="families", entity_ids=",".join(f"F{i}" for i in fams_0),
               definition=zero_def,
               interpretation_boundary=zero_bound,
               authoritative_source_path=rel(f5canon), authoritative_source_sha256=h5c)
    claims.add(claim_id="F5_2OF3_FAM", figure_id="5", panel_id="B", analysis_layer="triangulation",
               claim_name="support_2of3_families", value_text="F9", value_numeric=1, denominator=16,
               unit="families", entity_ids="F9", definition="Families with 2/3 gates",
               interpretation_boundary="Not Path A",
               authoritative_source_path=rel(f5canon), authoritative_source_sha256=h5c)
    claims.add(claim_id="F5_0OF3_FAM", figure_id="5", panel_id="B", analysis_layer="triangulation",
               claim_name="support_0of3_families", value_text="F4,F5,F6", value_numeric=3, denominator=16,
               unit="families", entity_ids="F4,F5,F6", definition=zero_def,
               interpretation_boundary=zero_bound,
               authoritative_source_path=rel(f5canon), authoritative_source_sha256=h5c)
    if int(float(ra["n_Path_A_gates"])) != 3:
        raise SystemExit("FATAL: n_Path_A_gates != 3")
    claims.add(claim_id="F5_N_GATES", figure_id="5", panel_id="A", analysis_layer="Path_A",
               claim_name="n_Path_A_gates",
               value_text="Path A has exactly three gates: bulk High, snRNA non-Below-detection, duration non-Stable",
               value_numeric=3, unit="gates", entity_ids=ra["Path_A_gate_conditions"],
               definition="Exactly three Path A gates", interpretation_boundary=ra["Path_A_gate_conditions"],
               authoritative_source_path=rel(f5a), authoritative_source_sha256=sha256(f5a))
    for cid, key, label in [
        ("F5_HCA_NOT_GATE", "hca_is_gate", "HCA carrier does not enter support count"),
        ("F5_STATE_NOT_GATE", "dominant_state_is_gate", "bulk dominant state does not enter support count"),
        ("F5_GEOMX_NOT_GATE", "geomx_is_gate", "GeoMx does not enter support count or Path A"),
        ("F5_GENETICS_NOT_GATE", "genetics_is_gate", "genetics does not enter support count or Path A"),
    ]:
        if str(ra[key]).lower() not in {"false", "0"}:
            raise SystemExit(f"FATAL: expected {key}=False, got {ra[key]}")
        claims.add(claim_id=cid, figure_id="5", panel_id="A", analysis_layer="Path_A_contract",
                   claim_name=key, value_text=label, value_numeric=0, entity_ids=key,
                   definition=label, interpretation_boundary="Non-gate evidence layers",
                   authoritative_source_path=rel(f5a), authoritative_source_sha256=sha256(f5a))
    gsum = ROOT / "results/ncp_v3_singscore_reanalysis/fig5d_geomx_roi_spatial_evidence_summary.tsv"
    gs = pd.read_csv(gsum, sep="\t")
    hg = sha256(gsum)
    n_rois = int(gs["n_rois"].iloc[0])
    n_pat = int(gs["n_patients"].iloc[0])
    if (n_rois, n_pat) != (24, 4):
        raise SystemExit(f"FATAL: GeoMx eligible {n_rois}/{n_pat} != 24/4")
    claims.add(claim_id="F5_GEOMX_ELIGIBLE", figure_id="5", panel_id="D", analysis_layer="GeoMx",
               claim_name="eligible_ROIs", value_text="24 ROIs from 4 patients", value_numeric=24,
               unit="ROIs", entity_ids="GeoMx", definition="Eligible cardiac AA DNA-positive ROI subset",
               interpretation_boundary="Outside Path A",
               authoritative_source_path=rel(gsum), authoritative_source_sha256=hg)
    for fam, cid, use_delta in [("F10", "F5_GEOMX_F10", False), ("F15", "F5_GEOMX_F15", False), ("F9", "F5_GEOMX_F9", True)]:
        row = gs.loc[gs["name"].astype(str) == fam].iloc[0]
        if use_delta:
            val = f"delta={float(row['delta']):.6f}; p={float(row['p']):.6f}"
        else:
            val = f"rho={float(row['rho']):.6f}; p={float(row['p']):.6f}"
        claims.add(claim_id=cid, figure_id="5", panel_id="D", analysis_layer="GeoMx",
                   claim_name=f"{fam}_GeoMx_stats", value_text=val, entity_ids=fam,
                   definition=f"GeoMx summary for {fam}", interpretation_boundary="Spatial evidence; not a Path A gate",
                   authoritative_source_path=rel(gsum), authoritative_source_sha256=hg)

    # ----- FIGURE 6 -----
    gwas = ROOT / "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv"
    gv = kv_tsv(gwas)
    hgwas = sha256(gwas)
    if gv["accession"] != "GCST90624411" or gv["locked_af_cases"] != "252438" or gv["locked_controls"] != "1959739":
        raise SystemExit(f"FATAL: GWAS provenance unexpected {gv}")
    claims.add(claim_id="F6_GWAS_PUB", figure_id="6", panel_id="A", analysis_layer="GWAS_provenance",
               claim_name="Yuan_2025_cross_population", value_text=gv["publication"], entity_ids="Yuan2025",
               definition="Publication for cross-population AF GWAS",
               interpretation_boundary="Catalog accession identifies deposited sumstats",
               authoritative_source_path=rel(gwas), authoritative_source_sha256=hgwas)
    claims.add(claim_id="F6_GCST", figure_id="6", panel_id="A", analysis_layer="GWAS_provenance",
               claim_name="accession", value_text="GCST90624411", entity_ids="GCST90624411",
               definition="GWAS Catalog accession", interpretation_boundary="File identity",
               authoritative_source_path=rel(gwas), authoritative_source_sha256=hgwas)
    claims.add(claim_id="F6_N_CASES", figure_id="6", panel_id="A", analysis_layer="GWAS_provenance",
               claim_name="locked_af_cases", value_text="252438", value_numeric=252438, unit="cases",
               entity_ids="AF_cases", definition="Corrected AF case N",
               interpretation_boundary="Replaces erroneous 168007",
               authoritative_source_path=rel(gwas), authoritative_source_sha256=hgwas)
    claims.add(claim_id="F6_N_CONTROLS", figure_id="6", panel_id="A", analysis_layer="GWAS_provenance",
               claim_name="locked_controls", value_text="1959739", value_numeric=1959739, unit="controls",
               entity_ids="controls", definition="Corrected control N",
               interpretation_boundary="Author-correction locked N",
               authoritative_source_path=rel(gwas), authoritative_source_sha256=hgwas)
    a_path = ROOT / "results/fig6a_anchor_coloc_calibration_source.tsv"
    b_path = ROOT / "results/fig6b_posterior_lead_distance_source.tsv"
    a = pd.read_csv(a_path, sep="\t")
    b = pd.read_csv(b_path, sep="\t")
    ha, hb = sha256(a_path), sha256(b_path)
    for col in ("local_AF_GWAS_support", "atrial_eGene_bonferroni", "PP4_ge_0.70"):
        a[col] = a[col].astype(str).str.lower().isin(["true", "1"])
    for col in ("pp3_dominant", "pp4_support"):
        b[col] = b[col].astype(str).str.lower().isin(["true", "1"])
    gwas_n = int(a["local_AF_GWAS_support"].sum())
    egene_n = int(a["atrial_eGene_bonferroni"].sum())
    pp4_a = int(a["PP4_ge_0.70"].sum())
    pp3_n = int(b["pp3_dominant"].sum())
    pp4_b = int(b["pp4_support"].sum())
    if pp4_a != pp4_b:
        raise SystemExit(f"FATAL: Fig6 PP4 A={pp4_a} vs B={pp4_b}")
    if (gwas_n, egene_n, pp3_n, pp4_a) != (7, 1, 12, 0):
        raise SystemExit(f"FATAL: Fig6 A/B counts {(gwas_n, egene_n, pp3_n, pp4_a)}")
    claims.add(claim_id="F6_GWAS_7OF12", figure_id="6", panel_id="A", analysis_layer="genetic_traceability",
               claim_name="local_AF_GWAS", value_text="7/12", value_numeric=7, denominator=12, unit="anchors",
               entity_ids="12_anchors", definition="Local AF GWAS support", interpretation_boundary="Traceability flag; not coloc",
               authoritative_source_path=rel(a_path), authoritative_source_sha256=ha)
    claims.add(claim_id="F6_EGENE_1OF12", figure_id="6", panel_id="A", analysis_layer="genetic_traceability",
               claim_name="atrial_eGene", value_text="1/12", value_numeric=1, denominator=12, unit="anchors",
               entity_ids="SMAD9", definition="Bonferroni atrial eGene", interpretation_boundary="eGene ≠ coloc",
               authoritative_source_path=rel(a_path), authoritative_source_sha256=ha)
    claims.add(claim_id="F6_PP3_12OF12", figure_id="6", panel_id="B", analysis_layer="coloc",
               claim_name="PP3_dominant", value_text="12/12", value_numeric=12, denominator=12, unit="anchors",
               entity_ids="12_anchors", definition="PP3-dominant configuration",
               interpretation_boundary="Distinct-local-signals posterior shape",
               authoritative_source_path=rel(b_path), authoritative_source_sha256=hb)
    claims.add(claim_id="F6_PP4_0OF12", figure_id="6", panel_id="A/B", analysis_layer="coloc",
               claim_name="PP4_ge_0.70", value_text="0/12", value_numeric=0, denominator=12, unit="anchors",
               entity_ids="12_anchors", definition="PP4 ≥ 0.70 support",
               interpretation_boundary="Calibration-limited non-support",
               authoritative_source_path=rel(a_path), authoritative_source_sha256=ha)
    dmin, dmax = float(b["delta_mb"].min()), float(b["delta_mb"].max())
    k = float(b.loc[b["gene"] == "KCNN3", "delta_mb"].iloc[0])
    claims.add(claim_id="F6_LEAD_RANGE", figure_id="6", panel_id="B", analysis_layer="lead_distance",
               claim_name="lead_separation_range_Mb", value_text=f"{dmin:.2f}–{dmax:.2f} Mb", unit="Mb",
               entity_ids="12_anchors", definition="Lead GWAS–eQTL separation range",
               interpretation_boundary="Descriptive; not in coloc posterior",
               authoritative_source_path=rel(b_path), authoritative_source_sha256=hb)
    claims.add(claim_id="F6_KCNN3_DIST", figure_id="6", panel_id="B", analysis_layer="lead_distance",
               claim_name="KCNN3_lead_distance_Mb", value_text=str(k), value_numeric=k, unit="Mb",
               entity_ids="KCNN3", definition="KCNN3 absolute lead separation", interpretation_boundary="Descriptive",
               authoritative_source_path=rel(b_path), authoritative_source_sha256=hb)
    c_path = ROOT / "results/fig6c_family_af_overlap_source.tsv"
    c = pd.read_csv(c_path, sep="\t")
    hc = sha256(c_path)
    n_tests, n_nz = len(c), int((c["n_overlap"] > 0).sum())
    n_fdr = int((c["fdr_bh"] < 0.05).sum())
    fams_c = sorted(c.loc[c["n_overlap"] > 0, "family_id"].astype(int).unique().tolist())
    genes = sorted({g for x in c.loc[c["n_overlap"] > 0, "overlap_genes"] for g in str(x).split(";") if g})
    univ = int(c["universe_n"].iloc[0])
    if (n_tests, n_nz, n_fdr, fams_c, set(genes), univ) != (64, 5, 0, [3, 9, 13], {"CAV1", "HEY2"}, 14635):
        raise SystemExit(f"FATAL: Fig6C unexpected {(n_tests, n_nz, n_fdr, fams_c, genes, univ)}")
    for cid, vt, vn, den, unit, ent, name in [
        ("F6C_N_TESTS", "64", 64, "", "tests", "family_x_set", "n_tested_pairs"),
        ("F6C_N_NONZERO", "5", 5, 64, "cells", "F3,F9,F13", "nonzero_cells"),
        ("F6C_UNIQUE_GENES", "CAV1, HEY2", 2, "", "genes", "CAV1,HEY2", "unique_genes"),
        ("F6C_OVERLAP_FAMS", "F3, F13, F9", 3, 16, "families", "F3,F13,F9", "overlap_families"),
        ("F6C_FDR", "0/64", 0, 64, "tests", "64_pairs", "BH_FDR_significant"),
        ("F6C_UNIVERSE", "14635", 14635, "", "genes", "universe", "universe_n"),
    ]:
        claims.add(claim_id=cid, figure_id="6", panel_id="C", analysis_layer="external_overlap",
                   claim_name=name, value_text=vt, value_numeric=vn, denominator=den, unit=unit, entity_ids=ent,
                   definition=name, interpretation_boundary="Member-overlap cue; not enrichment validation",
                   authoritative_source_path=rel(c_path), authoritative_source_sha256=hc)
    dsum_path = ROOT / "results/fig6d_non_nested_traceability_summary.tsv"
    dsum = pd.read_csv(dsum_path, sep="\t")
    hd = sha256(dsum_path)
    dmap = {str(r.dimension): (int(r.present), int(r.n)) for r in dsum.itertuples(index=False)}
    expected = {
        "local_AF_GWAS": (7, 12), "atrial_eGene": (1, 12), "same_variant_coloc": (0, 12),
        "common_risk": (2, 12), "regulatory": (2, 12), "rare_functional": (0, 12),
    }
    for k, exp in expected.items():
        if dmap.get(k) != exp:
            raise SystemExit(f"FATAL: Fig6D {k}={dmap.get(k)} != {exp}")
    claims.add(claim_id="F6D_CORE", figure_id="6", panel_id="D", analysis_layer="non_nested_traceability",
               claim_name="G_E_C_counts", value_text="7/12, 1/12, 0/12", denominator=12, unit="anchors",
               entity_ids="G,E,C", definition="Core calibration present counts",
               interpretation_boundary="Non-nested; no G-tier",
               authoritative_source_path=rel(dsum_path), authoritative_source_sha256=hd)
    claims.add(claim_id="F6D_RISK", figure_id="6", panel_id="D", analysis_layer="non_nested_traceability",
               claim_name="risk_prioritized", value_text="2/12", value_numeric=2, denominator=12, unit="anchors",
               entity_ids="risk", definition="Common-risk/prioritized membership", interpretation_boundary="External flag",
               authoritative_source_path=rel(dsum_path), authoritative_source_sha256=hd)
    claims.add(claim_id="F6D_REG", figure_id="6", panel_id="D", analysis_layer="non_nested_traceability",
               claim_name="regulatory", value_text="2/12", value_numeric=2, denominator=12, unit="anchors",
               entity_ids="regulatory", definition="Regulatory-set membership", interpretation_boundary="External flag",
               authoritative_source_path=rel(dsum_path), authoritative_source_sha256=hd)
    claims.add(claim_id="F6D_RARE", figure_id="6", panel_id="D", analysis_layer="non_nested_traceability",
               claim_name="rare_functional", value_text="0/12", value_numeric=0, denominator=12, unit="anchors",
               entity_ids="rare", definition="Rare/functional-set membership", interpretation_boundary="External flag",
               authoritative_source_path=rel(dsum_path), authoritative_source_sha256=hd)
    for cid, vt in [
        ("F6_NON_NESTED", "genetic evidence categories are non-nested"),
        ("F6_NO_GTIER", "no G0–G4 tier"),
        ("F6_NO_SCORE", "no composite genetic score"),
        ("F6_NOT_PATHA", "genetics is not a Path A gate"),
    ]:
        claims.add(claim_id=cid, figure_id="6", panel_id="D", analysis_layer="interpretation_boundary",
                   claim_name=cid, value_text=vt, entity_ids="NA", definition=vt,
                   interpretation_boundary="Figure 6 boundary language",
                   authoritative_source_path=rel(dsum_path), authoritative_source_sha256=hd)


def build_conflicts() -> tuple[list[dict], dict[str, int]]:
    rows: list[dict] = []

    def add(severity, file, location, observed, expected, auth, action):
        rows.append({
            "severity": severity, "file": file, "location": location,
            "observed": observed, "expected": expected,
            "authoritative_source": auth, "proposed_migration_action": action,
        })

    sub = ROOT / "submission/PLOS_CB_submission_pack/02_Main_figures"
    # Figure 5: distinguish legacy submission panel filenames (A–E) from final scientific A–D structure.
    if not list(sub.glob("Figure_5_figure5*")) and list(sub.glob("Figure_5A*")):
        add(
            "MAJOR",
            "submission/PLOS_CB_submission_pack/02_Main_figures/",
            "Figure 5 routing (legacy panel filenames vs final A–D structure)",
            "Submission pack contains legacy standalone Figure 5 panel assets A–E; "
            "the final A–D composite is not routed.",
            "The authoritative final Figure 5 structure is A–D. At routing freeze, ship the "
            "canonical A–D composite or explicit A–D panel assets only; Figure 5E must not "
            "remain in formal routing.",
            "figures/ncp_v3_singscore_reanalysis/figure5_cross_layer_triangulation.pdf",
            "Defer asset changes to the final routing phase.",
        )

    c6 = ROOT / "figures/figure6_genetic_traceability.pdf"
    s6 = sub / "Figure_6_figure6_genetic_traceability.pdf"
    if c6.is_file() and s6.is_file() and sha256(c6) != sha256(s6):
        add("BLOCKER", rel(s6), "Figure 6 composite hash",
            f"submission sha {sha256(s6)[:16]}…", f"canonical sha {sha256(c6)[:16]}…",
            rel(c6), "Rebuild submission pack after Figure 6 2×2 QC sign-off (do not edit bundle)")

    if (ROOT / "figures/fig6b_pp3_pp4_leaddistance.pdf").is_file():
        add("MAJOR", "figures/fig6b_pp3_pp4_leaddistance.pdf", "legacy Fig6B",
            "legacy file still on disk", "canonical fig6b_posterior_lead_distance.pdf",
            "figures/fig6b_posterior_lead_distance.pdf",
            "Keep archived; ensure no routing points to legacy")

    layout = ROOT / "results/figure_layout_audit.tsv"
    if layout.is_file():
        t = layout.read_text(encoding="utf-8", errors="ignore")
        if "pp3_pp4_leaddistance" in t or "g_tier" in t:
            add("MAJOR", rel(layout), "Figure 6 rows",
                "stale legacy 6B/G-tier paths", "posterior_lead_distance + non-nested D",
                "results/figure6_composite_manifest.tsv",
                "Regenerate layout audit after routing freeze")

    repro = ROOT / "results/reproducibility_table.tsv"
    if repro.is_file():
        t = repro.read_text(encoding="utf-8", errors="ignore")
        if "assemble_figure5_cross_layer" not in t and "fig5a_triangulation.py" in t:
            add("MAJOR", rel(repro), "Fig 5 scripts",
                "legacy panel drivers listed", "v3 assemble_figure5 composite path",
                "results/ncp_v3_singscore_reanalysis/figure5_composite_manifest.tsv",
                "Update reproducibility table rows after manuscript freeze")

    m4 = ROOT / "scripts/m4_results_docx.py"
    if m4.is_file():
        t = m4.read_text(encoding="utf-8", errors="ignore")
        if "G-tier" in t or re.search(r"\bG[0-4]\b", t):
            add("MAJOR", rel(m4), "Results generator",
                "G-tier/G0–G4 language present", "non-nested Fig6D wording only",
                "results/fig6d_non_nested_traceability_summary.tsv",
                "Regenerate Results.docx after clearing legacy G-tier text")

    m2 = ROOT / "scripts/m2_figure_captions_docx.py"
    if m2.is_file():
        t = m2.read_text(encoding="utf-8", errors="ignore")
        # crude panel-d check for Figure 6 legend builder
        if "Figure 6" in t and re.search(r"figure.?6", t, re.I):
            chunk = t.split("Figure 6", 1)[-1][:4000]
            if '"(d)"' not in chunk and "(\"d\")" not in chunk and "Panel D" not in chunk and "panel D" not in chunk:
                add("MAJOR", rel(m2), "Figure 6 legend",
                    "possible A–C-only structure", "final A–D legend",
                    "manuscript_drafts/Figure_Captions.docx",
                    "Verify captions generator encodes panels A–D")

    if (ROOT / "results/figure6_g_tier_traceability.tsv").is_file():
        add("MINOR", "results/figure6_g_tier_traceability.tsv", "legacy G-tier table",
            "file still present", "non-nested source only for publication",
            "results/fig6d_non_nested_traceability_source.tsv",
            "Do not cite as final; archive already exists")

    for fid, meta in FIGURES.items():
        svg = meta["svg"]
        if svg.is_file():
            vis = re.sub(r"<[^>]+>", " ", svg.read_text(encoding="utf-8", errors="ignore"))
            if "NCP_v3" in vis or "ncp_v3" in vis:
                add("BLOCKER", rel(svg), "publication SVG text",
                    "NCP_v3 visible", "no publication-facing NCP_v3",
                    rel(meta["qc"]), "Remove/replace visible string in panel renderer")
            if "PriorityScore" in vis or "AtlasPriority" in vis:
                add("BLOCKER", rel(svg), "publication SVG text",
                    "Atlas PriorityScore visible", "not used for ordering",
                    rel(meta["qc"]), "Remove PriorityScore from publication assets")
            if re.search(r"\bmr_|\bMiyazawa\b", vis):
                add("BLOCKER", rel(svg), "publication SVG text",
                    "MR/Miyazawa visible", "absent from publication assets",
                    rel(meta["qc"]), "Strip internal genetics filenames from figure text")

    f5s = ROOT / "results/ncp_v3_singscore_reanalysis/fig5b_cross_layer_triangulation_summary.tsv"
    if f5s.is_file():
        pa = kv_tsv(f5s).get("Path_A_families", "")
        if pa.strip() in {"F15", "15"}:
            add("BLOCKER", rel(f5s), "Path_A_families", pa, "F10,F15", rel(f5s), "Restore F10,F15 Path A")

    for gen in [ROOT / "scripts/m4_results_docx.py", ROOT / "scripts/m1_methods_docx.py",
                ROOT / "scripts/m2_figure_captions_docx.py"]:
        if gen.is_file() and "168007" in gen.read_text(encoding="utf-8", errors="ignore"):
            add("BLOCKER", rel(gen), "GWAS N", "168007 present", "252438 cases / 1959739 controls",
                "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv",
                "Replace legacy case N before regenerating manuscript")

    # Yuan/Roselli conflation scan in generators
    for gen in [ROOT / "scripts/m4_results_docx.py", ROOT / "scripts/m1_methods_docx.py",
                ROOT / "scripts/m2_figure_captions_docx.py"]:
        if not gen.is_file():
            continue
        t = gen.read_text(encoding="utf-8", errors="ignore")
        if "Roselli" in t and ("GCST90624411" in t or "Yuan" in t):
            # only flag if they appear near each other
            if re.search(r"Roselli.{0,200}Yuan|Yuan.{0,200}Roselli", t, re.S):
                add("MAJOR", rel(gen), "GWAS citation",
                    "Yuan/Roselli appear near each other", "Yuan 2025 for GCST90624411; Roselli separate if cited",
                    "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv",
                    "Disambiguate citations in manuscript generators")

    counts = {"BLOCKER": 0, "MAJOR": 0, "MINOR": 0, "EDITORIAL": 0}
    for r in rows:
        counts[r["severity"]] = counts.get(r["severity"], 0) + 1
    return rows, counts



def write_conflicts_md(conflicts: list[dict], counts: dict[str, int]) -> str:
    lines = [
        f"# Unresolved conflicts — {BUNDLE_NAME}",
        "",
        "Read-only audit. No manuscript/routing/figure edits were performed by the builder.",
        "",
        f"| Severity | Count |",
        f"|----------|------:|",
        f"| BLOCKER | {counts.get('BLOCKER', 0)} |",
        f"| MAJOR | {counts.get('MAJOR', 0)} |",
        f"| MINOR | {counts.get('MINOR', 0)} |",
        f"| EDITORIAL | {counts.get('EDITORIAL', 0)} |",
        "",
    ]
    order = {"BLOCKER": 0, "MAJOR": 1, "MINOR": 2, "EDITORIAL": 3}
    for i, c in enumerate(sorted(conflicts, key=lambda x: (order.get(x["severity"], 9), x["file"])), 1):
        lines.extend([
            f"## {i}. [{c['severity']}] {c['location']}",
            "",
            f"- **file**: `{c['file']}`",
            f"- **observed**: {c['observed']}",
            f"- **expected**: {c['expected']}",
            f"- **authoritative_source**: `{c['authoritative_source']}`",
            f"- **proposed_migration_action**: {c['proposed_migration_action']}",
            "",
        ])
    if not conflicts:
        lines.append("No unresolved conflicts detected by automated checks.")
        lines.append("")
    return "\n".join(lines)


def assert_figure5_panels_ad_only() -> None:
    """Assert final Figure 5 composite is A–D only (no panel E)."""
    meta = FIGURES["5"]
    if meta["panels"] != "A–D":
        raise SystemExit(f"FATAL: Figure 5 panel_structure must be A–D, got {meta['panels']!r}")
    man = meta.get("composite_manifest")
    if man is None or not Path(man).is_file():
        raise SystemExit("FATAL: Figure 5 composite_manifest missing")
    mdf = pd.read_csv(man, sep="\t")
    panels = sorted(mdf["panel"].astype(str).unique().tolist())
    if panels != ["A", "B", "C", "D"]:
        raise SystemExit(f"FATAL: Figure 5 composite_manifest panels={panels}, expected A–D only")
    if "E" in panels:
        raise SystemExit("FATAL: Figure 5 composite_manifest contains panel E")
    svg = meta["svg"]
    if svg.is_file():
        text = svg.read_text(encoding="utf-8", errors="ignore")
        if re.search(r'id="panel-E"|Panel E|panel E', text):
            raise SystemExit("FATAL: Figure 5 SVG contains panel E markers")
        ids = sorted(set(re.findall(r'id="panel-([A-E])"', text)))
        if ids and ids != ["A", "B", "C", "D"]:
            raise SystemExit(f"FATAL: Figure 5 SVG panel ids={ids}, expected A–D only")


def main() -> None:
    if OUT.resolve() == LEGACY_V1.resolve():
        raise SystemExit("FATAL: refusing to overwrite fulltext_audit_v1; use v1_1 namespace")
    if BUNDLE_NAME == "fulltext_audit_v1":
        raise SystemExit("FATAL: BUNDLE_NAME must not be fulltext_audit_v1")
    assert_figure5_panels_ad_only()
    if OUT.exists():
        shutil.rmtree(OUT)
    global UTC
    UTC = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    OUT.mkdir(parents=True)
    PREV.mkdir(parents=True)
    MSRC.mkdir(parents=True)
    AUTH.mkdir(parents=True)

    branch = git("rev-parse", "--abbrev-ref", "HEAD") or "NA"
    commit = git("rev-parse", "HEAD") or "NA"
    porcelain = git("status", "--porcelain")
    tracked = tracked_set()
    status = status_map()
    py_ver = sys.version.split()[0]
    r_ver = run_version(["R", "--version"])
    if r_ver == "unavailable":
        r_ver = run_version(["Rscript", "--version"])

    bundle = BundleFiles(tracked, status)
    claims = ClaimRegistry()
    build_claims(claims)

    fig_rows = []
    routing_rows = []
    sub_dir = ROOT / "submission/PLOS_CB_submission_pack/02_Main_figures"
    builder = ROOT / "scripts/build_plos_submission_folder.py"
    builder_txt = builder.read_text(encoding="utf-8") if builder.is_file() else ""

    submission_map = {
        "1": "Figure_1_figure1_composite.pdf",
        "2": "Figure_2_figure2_composite.pdf",
        "3": "Figure_3_figure3_la_single_nucleus_localization.pdf",
        "4": "Figure_4_figure4_cross_sectional_duration.pdf",
        "5": None,
        "6": "Figure_6_figure6_genetic_traceability.pdf",
    }

    for fid, meta in FIGURES.items():
        pdf, svg, png, qc = meta["pdf"], meta["svg"], meta["png"], meta["qc"]
        if not pdf.is_file():
            raise SystemExit(f"FATAL: missing final composite PDF for Figure {fid}: {pdf}")
        # uniqueness: FIGURES defines exactly one composite per figure_id
        qstat = qc_pass(qc)
        srcs = [p for p in meta["sources"] if p.is_file()]
        for p in srcs:
            dest = AUTH / Path(*p.relative_to(ROOT).parts)
            copy_small(p, dest, bundle, "authoritative_source_table", fid, True)
        if png.is_file():
            preview_png(png, PREV / f"figure{fid}_preview.png", bundle, fid)
        elif pdf.is_file() and shutil.which("mutool"):
            tmp = PREV / f"_tmp{fid}.png"
            subprocess.run(
                ["mutool", "draw", "-q", "-r", "100", "-o", str(tmp), str(pdf)],
                check=False, capture_output=True,
            )
            if tmp.is_file():
                preview_png(tmp, PREV / f"figure{fid}_preview.png", bundle, fid)
                tmp.unlink(missing_ok=True)
        w_mm, h_mm = canvas_mm(svg)

        # routing
        sub_name = submission_map[fid]
        if fid == "5":
            panels = sorted(sub_dir.glob("Figure_5*.pdf")) if sub_dir.is_dir() else []
            sub_path = None
            for p in panels:
                if "figure5_cross_layer" in p.name:
                    sub_path = p
                    break
            if sub_path is None:
                route_status = "partially_routed"
                hash_match = False
                sub_rel = (
                    "legacy standalone Figure_5A–Figure_5E panel assets "
                    "(final A–D composite not routed)"
                )
                sub_sha = ""
                legacy = "yes"
            else:
                hash_match = sha256(sub_path) == sha256(pdf)
                route_status = "fully_routed" if hash_match else "hash_mismatched"
                sub_rel = rel(sub_path)
                sub_sha = sha256(sub_path)
                legacy = "no"
        else:
            sub_path = sub_dir / sub_name if sub_name and sub_dir.is_dir() else None
            if sub_path and sub_path.is_file():
                hash_match = sha256(sub_path) == sha256(pdf)
                route_status = "fully_routed" if hash_match else "hash_mismatched"
                sub_rel = rel(sub_path)
                sub_sha = sha256(sub_path)
            else:
                hash_match = False
                route_status = "not_routed"
                sub_rel = sub_name or "MISSING"
                sub_sha = ""
            legacy = "yes" if fid == "6" and (ROOT / "figures/fig6b_pp3_pp4_leaddistance.pdf").is_file() else "no"

        fig_rows.append({
            "figure_id": fid,
            "panel_structure": meta["panels"],
            "final_title": meta["title"],
            "pdf_path": rel(pdf),
            "svg_path": rel(svg) if svg.is_file() else "MISSING",
            "png_path": rel(png) if png.is_file() else "MISSING",
            "pdf_sha256": sha256(pdf),
            "svg_sha256": sha256(svg) if svg.is_file() else "",
            "png_sha256": sha256(png) if png.is_file() else "",
            "canvas_width_mm": w_mm,
            "canvas_height_mm": h_mm,
            "png_dpi": 600 if png.is_file() else "",
            "renderer_script": meta["renderer"],
            "source_table_paths": ";".join(rel(p) for p in srcs),
            "source_table_sha256s": ";".join(sha256(p) for p in srcs),
            "qc_report_path": rel(qc) if qc.is_file() else "MISSING",
            "freeze_status": "scientifically_frozen" if qstat == "PASS" else qstat,
            "routing_status": route_status,
        })
        routing_rows.append({
            "figure_id": fid,
            "canonical_asset_path": rel(pdf),
            "canonical_sha256": sha256(pdf),
            "submission_asset_path": sub_rel,
            "submission_sha256": sub_sha,
            "builder_source": "scripts/build_plos_submission_folder.py",
            "routing_status": route_status,
            "legacy_asset_still_active": legacy,
            "hash_match": hash_match,
            "notes": (
                "builder includes figure6_genetic_traceability"
                if fid == "6" and "figure6_genetic_traceability" in builder_txt
                else (
                    "Final scientific structure is A–D; submission currently retains "
                    "legacy standalone panel filenames A–E without the routed A–D composite."
                    if fid == "5"
                    else ""
                )
            ),
        })

    # Provenance (copy small + manifest)
    prov_items = [
        ("frozen_program_metadata", ROOT / "results/ncp_v3_singscore_reanalysis/frozen_program_metadata.tsv", True, False, "Frozen 30-program metadata"),
        ("frozen_family_metadata", ROOT / "results/ncp_v3_family_labels.tsv", True, False, "Frozen 16-family labels"),
        ("figure1d_canonical_family_order", ROOT / "results/ncp_v3_canonical_display_order.tsv", True, False, "Figure 1D canonical family order"),
        ("selection_time_leading_edges", ROOT / "results/ncp_v3_singscore_reanalysis/selection_time_leading_edges.tsv", True, False, "Selection-time leading-edge provenance"),
        ("family_leading_edge", ROOT / "results/ncp_v3_family_leading_edge_genes.tsv", True, False, "Program-first family aggregation LE genes"),
        ("scorer_equivalence_summary", ROOT / "results/fig3_scorer_equivalence_summary.tsv", True, False, "R singscore 1.28.1 equivalence evidence"),
        ("innovation_ncp_v3_summary", ROOT / "results/innovation_ncp_v3_summary.tsv", True, False, "Program-first family aggregation contract summary"),
        ("geomx_gene_union_exception", ROOT / "results/ncp_v3_singscore_reanalysis/figure5_geomx_lineage_audit.tsv", True, False, "GeoMx gene-union scoring exception / lineage"),
        ("gse255612_lineage", ROOT / "results/fig3_publication_authority_lock.md", True, False, "GSE255612 corrected lineage / Fig3 authority lock"),
        ("gse41177_lineage", ROOT / "results/ncp_v3_singscore_reanalysis/fig4a_duration_design_summary.tsv", True, False, "GSE41177 corrected duration design lineage"),
        ("gwas_gcst90624411_provenance", ROOT / "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv", True, False, "GCST90624411 provenance package"),
        ("figure6_external_af_gene_sets_manifest", ROOT / "results/provenance/figure6_external_gene_sets_manifest.tsv", True, False, "Figure 6 external AF gene-set manifest"),
        ("geomx_eligibility_source", ROOT / "results/ncp_v3_singscore_reanalysis/fig5d_geomx_roi_spatial_evidence_summary.tsv", True, False, "GeoMx eligibility/source summary"),
        ("path_a_gate_definition", ROOT / "results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv", True, False, "Path A gate definition"),
        ("final_family_label_mapping", ROOT / "results/ncp_v3_family_labels.tsv", True, False, "Final family-label mapping"),
        ("reproducibility_table", ROOT / "results/reproducibility_table.tsv", True, False, "Reproducibility table (routing audit input)"),
        ("figure_layout_audit", ROOT / "results/figure_layout_audit.tsv", True, False, "Figure layout audit (routing audit input)"),
        ("submission_readme_index", ROOT / "submission/PLOS_CB_submission_pack/00_README_SUBMISSION_INDEX.txt", True, True, "Submission README/index"),
    ]
    prov_rows = []
    for role, p, auth, vis, notes in prov_items:
        if p.is_file():
            dest = AUTH / Path(*p.relative_to(ROOT).parts)
            copy_small(p, dest, bundle, "provenance_source", "", auth, notes)
            # scorer version from F1 summary if available
            ver = "current_snapshot"
            if "scorer" in role or role == "innovation_ncp_v3_summary":
                ver = "singscore_1.28.1"
            prov_rows.append({
                "artifact_role": role,
                "artifact_path": rel(p),
                "sha256": sha256(p),
                "version": ver,
                "authoritative": auth,
                "publication_visible": vis,
                "notes": notes,
            })
        else:
            prov_rows.append({
                "artifact_role": role,
                "artifact_path": rel(p) if p else "",
                "sha256": "MISSING",
                "version": "current_snapshot",
                "authoritative": auth,
                "publication_visible": vis,
                "notes": notes + "; missing_at_snapshot",
            })

    # Manuscript sources
    ms_spec = [
        ("Abstract", ROOT / "manuscript_drafts/Abstract_AuthorSummary_plos_cb.md", "", "scripts/", True, "current", "editable Abstract+Author Summary MD"),
        ("Author Summary", ROOT / "manuscript_drafts/Abstract_AuthorSummary_plos_cb.md", "", "scripts/", True, "current", "same file as Abstract"),
        ("Introduction", ROOT / "manuscript_drafts/Introduction.docx", "scripts/m6_introduction_docx.py", "", True, "current", "editable generated DOCX"),
        ("Results", ROOT / "manuscript_drafts/Results.docx", "scripts/m4_results_docx.py", "", True, "current", "editable generated DOCX"),
        ("Methods", ROOT / "manuscript_drafts/Methods.docx", "scripts/m1_methods_docx.py", "", True, "current", "editable generated DOCX"),
        ("Discussion", ROOT / "manuscript_drafts/Discussion.docx", "scripts/m5_discussion_docx.py", "", True, "current", "editable generated DOCX"),
        ("Figure 1–6 legends", ROOT / "manuscript_drafts/Figure_Captions.docx", "scripts/m2_figure_captions_docx.py", "", True, "current", "figure legends"),
        ("Supporting Information legends", ROOT / "manuscript_drafts/Figure_Captions.docx", "scripts/m2_figure_captions_docx.py", "", True, "current", "SI legends share captions doc if present"),
        ("Data availability", ROOT / "manuscript_drafts/Data_Code_Availability.docx", "scripts/m3_data_availability_docx.py", "", True, "current", "data+code availability"),
        ("Code availability", ROOT / "manuscript_drafts/Data_Code_Availability.docx", "scripts/m3_data_availability_docx.py", "", True, "current", "data+code availability"),
        ("latest manuscript DOCX", ROOT / "submission/PLOS_CB_submission_pack/01_Manuscript_Word/Manuscript_main.docx", "scripts/build_plos_submission_folder.py", "", True, "current", "submission-pack main manuscript DOCX"),
        ("latest manuscript Markdown", ROOT / "manuscript_drafts/RESULTS_draft_v3.md", "", "", False, "legacy", "draft MD; not necessarily publication-synced"),
        ("Results draft MD", ROOT / "manuscript_drafts/RESULTS_draft_v3.md", "", "", False, "legacy", "legacy draft"),
        ("Discussion draft MD", ROOT / "manuscript_drafts/Discussion_draft_v3.md", "", "", False, "legacy", "legacy draft"),
    ]
    # PDF if any
    pdf_candidates = list((ROOT / "submission/PLOS_CB_submission_pack").rglob("*.pdf"))
    ms_pdf = None
    for p in pdf_candidates:
        if "Manuscript" in p.name or "manuscript" in p.name.lower():
            ms_pdf = p
            break
    if ms_pdf is None:
        # no full manuscript PDF in pack — record MISSING
        ms_spec.append(("latest manuscript PDF", ROOT / "MISSING_manuscript.pdf", "", "", True, "current", "no full manuscript PDF found in submission pack"))

    ms_rows = []
    for section, src, gen, outp, pub, cur, notes in ms_spec:
        if not src.is_file():
            ms_rows.append({
                "section": section,
                "source_path": rel(src) if src.name != "MISSING_manuscript.pdf" else "MISSING",
                "source_sha256": "",
                "generator_script": gen,
                "generated_output_path": outp,
                "generated_output_sha256": "",
                "publication_facing": pub,
                "current_or_legacy": cur,
                "notes": notes + "; missing",
            })
            continue
        # copy small editable sources (md/docx/py) under manuscript_sources/
        if src.suffix.lower() in {".md", ".docx", ".txt"} and src.stat().st_size <= MAX_COPY_BYTES:
            dest = MSRC / Path(*src.relative_to(ROOT).parts)
            copy_small(src, dest, bundle, "manuscript_source", "", False, notes)
        gen_sha = ""
        gen_path = ROOT / gen if gen and (ROOT / gen).is_file() else None
        if gen_path and gen_path.is_file():
            gdest = MSRC / Path(*gen_path.relative_to(ROOT).parts)
            copy_small(gen_path, gdest, bundle, "manuscript_generator", "", False, "generator copy")
            gen_sha = sha256(gen_path)
        ms_rows.append({
            "section": section,
            "source_path": rel(src),
            "source_sha256": sha256(src),
            "generator_script": gen,
            "generated_output_path": outp or rel(src),
            "generated_output_sha256": sha256(src),
            "publication_facing": pub,
            "current_or_legacy": cur,
            "notes": notes,
        })

    # Also copy abstract md generator if exists (none) — done

    conflicts, counts = build_conflicts()

    # Determine routing switched?
    fig6_routed = any(r["figure_id"] == "6" and r["routing_status"] in {"fully_routed", "hash_mismatched"} for r in routing_rows)
    fig6_hash_ok = any(r["figure_id"] == "6" and r["hash_match"] for r in routing_rows)
    submission_switched = fig6_routed and "figure6_genetic_traceability" in builder_txt
    formal_routing_switched = submission_switched and fig6_hash_ok

    # Write manifests
    claims_df = claims.to_df()
    claims_df.to_csv(OUT / "02_final_claims_registry.tsv", sep="\t", index=False)
    bundle.record(OUT / "02_final_claims_registry.tsv", None, "claims_registry", notes="generated")

    fig_df = pd.DataFrame(fig_rows)
    fig_df.to_csv(OUT / "03_final_main_figure_manifest.tsv", sep="\t", index=False)
    bundle.record(OUT / "03_final_main_figure_manifest.tsv", None, "figure_manifest", notes="generated")

    ms_df = pd.DataFrame(ms_rows)
    ms_df.to_csv(OUT / "04_manuscript_source_manifest.tsv", sep="\t", index=False)
    bundle.record(OUT / "04_manuscript_source_manifest.tsv", None, "manuscript_manifest", notes="generated")

    route_df = pd.DataFrame(routing_rows)
    route_df.to_csv(OUT / "05_routing_manifest.tsv", sep="\t", index=False)
    bundle.record(OUT / "05_routing_manifest.tsv", None, "routing_manifest", notes="generated")

    prov_df = pd.DataFrame(prov_rows)
    prov_df.to_csv(OUT / "06_provenance_manifest.tsv", sep="\t", index=False)
    bundle.record(OUT / "06_provenance_manifest.tsv", None, "provenance_manifest", notes="generated")

    conflicts_md = write_conflicts_md(conflicts, counts)
    (OUT / "08_unresolved_conflicts.md").write_text(conflicts_md, encoding="utf-8")
    bundle.record(OUT / "08_unresolved_conflicts.md", None, "conflicts_report", notes="generated")

    fig_ns = "\n".join(
        f"| {fid} | `{meta['pdf'].relative_to(ROOT)}` | {next(r['freeze_status'] for r in fig_rows if r['figure_id']==fid)} | {next(r['routing_status'] for r in fig_rows if r['figure_id']==fid)} |"
        for fid, meta in FIGURES.items()
    )
    latest_ms = "submission/PLOS_CB_submission_pack/01_Manuscript_Word/Manuscript_main.docx"
    if (ROOT / latest_ms).is_file():
        latest_ms_line = f"`{latest_ms}` (sha256 `{sha256(ROOT / latest_ms)[:16]}…`)"
    else:
        latest_ms_line = "MISSING"

    readme = f"""# Full-text audit handoff bundle v1

Read-only snapshot for whole-manuscript reconciliation.

## Project snapshot

| Field | Value |
|-------|-------|
| project_root | `{ROOT}` |
| snapshot_time_utc | `{UTC}` |
| git_branch | `{branch}` |
| git_commit | `{commit}` |
| Python version | `{py_ver}` |
| R version | `{r_ver}` |
| git_status_porcelain lines | `{len(porcelain.splitlines()) if porcelain else 0}` |

Porcelain (first 40 lines):

```
{chr(10).join(porcelain.splitlines()[:40])}
```

## Final Figures 1–6 namespace

Scientifically frozen final composites live under the paths below (not standalone panels).

| Figure | Canonical composite | Freeze | Routing |
|--------|---------------------|--------|---------|
{fig_ns}

- **Figures 1–6 scientifically frozen?** Yes for source-locked composites with QC PASS (see `03_final_main_figure_manifest.tsv` `freeze_status`).
- **Formal submission routing switched?** {"Yes (builder points to final Fig6 composite and hashes match)" if formal_routing_switched else "No / incomplete — builder may reference finals but hash or packaging mismatches remain (see `05_routing_manifest.tsv`)"}.

## Latest manuscript output

- Main DOCX: {latest_ms_line}
- Section DOCX/MD: see `04_manuscript_source_manifest.tsv` and copies under `manuscript_sources/`
- Full manuscript PDF: not found in submission pack at snapshot (if absent, regenerate outside this bundle)

## Known unresolved manuscript/routing conflicts

See `08_unresolved_conflicts.md` (severity-tagged).

Counts at snapshot: BLOCKER={counts.get('BLOCKER',0)}, MAJOR={counts.get('MAJOR',0)}, MINOR={counts.get('MINOR',0)}, EDITORIAL={counts.get('EDITORIAL',0)}.

## Bundle guarantees

- **Read-only**: this builder writes only under `results/handoff/{BUNDLE_NAME}/` and the zip sidecars.
- **Does not overwrite** `results/handoff/fulltext_audit_v1/`.
- **No raw human-level data**: no expression matrices, CEL/FASTQ, Seurat/AnnData, full GWAS sumstats, or large intermediates.
- Claims in `02_final_claims_registry.tsv` are derived only from final authoritative TSVs/manifests/QC — not from manuscript prose, PDF text, or PNG/SVG appearance.
- Full-resolution figure PDFs are **not** copied into the bundle; compact previews are under `figure_previews/`.
- Figure 5 final composite panel structure is **A–D** (not A–E).

## Contents

1. `00_README.md` — this file
2. `01_bundle_file_manifest.tsv` — every file copied into the bundle
3. `02_final_claims_registry.tsv`
4. `03_final_main_figure_manifest.tsv`
5. `04_manuscript_source_manifest.tsv` + `manuscript_sources/`
6. `05_routing_manifest.tsv`
7. `06_provenance_manifest.tsv`
8. `07_authoritative_source_tables/`
9. `08_unresolved_conflicts.md`
10. `figure_previews/figure{{1..6}}_preview.png`

Builder: `scripts/audit/build_fulltext_audit_handoff.py`
"""
    (OUT / "00_README.md").write_text(readme, encoding="utf-8")
    bundle.record(OUT / "00_README.md", None, "readme", notes="generated")

    # Post-write Figure 5 / registry / routing assertions
    fig_man_path = OUT / "03_final_main_figure_manifest.tsv"
    f5_row = pd.read_csv(fig_man_path, sep="\t")
    f5_ps = str(f5_row.loc[f5_row["figure_id"].astype(str) == "5", "panel_structure"].iloc[0])
    if f5_ps != "A–D":
        raise SystemExit(f"FATAL: bundle figure manifest Fig5 panel_structure={f5_ps!r}, expected A–D")
    prev5 = PREV / "figure5_preview.png"
    if not prev5.is_file():
        raise SystemExit("FATAL: missing figure5_preview.png")
    # Scientific structure field must remain A–D (not described as A–E)
    if f5_ps != "A–D" or "A–E" in f5_ps:
        raise SystemExit("FATAL: Figure 5 panel_structure must be A–D only")

    claims_path = OUT / "02_final_claims_registry.tsv"
    cdf = pd.read_csv(claims_path, sep="\t")
    pub_cols = [
        "claim_id", "figure_id", "panel_id", "analysis_layer", "claim_name",
        "value_text", "value_numeric", "denominator", "unit", "entity_ids",
        "definition", "interpretation_boundary", "qc_status",
    ]
    for col in pub_cols:
        if col not in cdf.columns:
            continue
        bad = cdf[cdf[col].astype(str).str.contains(r"\bsilent\b", case=False, na=False, regex=True)]
        if len(bad):
            raise SystemExit(
                f"FATAL: publication-facing registry field {col} contains \\bsilent\\b: "
                f"{bad['claim_id'].tolist()}"
            )
    # Internal frozen summary metric 'silent' must remain readable and equal to 3
    f5s_chk = ROOT / "results/ncp_v3_singscore_reanalysis/fig5b_cross_layer_triangulation_summary.tsv"
    silent_metric = kv_tsv(f5s_chk).get("silent")
    if silent_metric is None or int(float(silent_metric)) != 3:
        raise SystemExit(f"FATAL: internal source metric silent changed/missing: {silent_metric!r}")

    r0n = cdf.loc[cdf["claim_id"] == "F5_0OF3_N"].iloc[0]
    r0f = cdf.loc[cdf["claim_id"] == "F5_0OF3_FAM"].iloc[0]
    r1 = cdf.loc[cdf["claim_id"] == "F5_1OF3"].iloc[0]
    if int(float(r0n["value_numeric"])) != 3:
        raise SystemExit("FATAL: F5_0OF3_N != 3")
    if str(r0f["entity_ids"]) != "F4,F5,F6" and str(r0f["value_text"]).replace(" ", "") != "F4,F5,F6":
        raise SystemExit(f"FATAL: F5_0OF3_FAM unexpected {r0f['entity_ids']!r} / {r0f['value_text']!r}")
    if str(r1["entity_ids"]) != "F1,F2,F3,F7,F8,F11,F12,F13,F14,F16":
        raise SystemExit(f"FATAL: F5_1OF3 entity_ids={r1['entity_ids']!r}")

    route_path = OUT / "05_routing_manifest.tsv"
    rdf = pd.read_csv(route_path, sep="\t")
    r5 = rdf.loc[rdf["figure_id"].astype(str) == "5"].iloc[0]
    if str(r5["routing_status"]) != "partially_routed":
        raise SystemExit(f"FATAL: Fig5 routing_status={r5['routing_status']!r}")
    if str(r5["legacy_asset_still_active"]).lower() not in {"yes", "true", "1"}:
        raise SystemExit(f"FATAL: Fig5 legacy_asset_still_active={r5['legacy_asset_still_active']!r}")
    conf_path = OUT / "08_unresolved_conflicts.md"
    conf = conf_path.read_text(encoding="utf-8")
    if "legacy standalone Figure 5 panel assets A–E" not in conf:
        raise SystemExit("FATAL: routing conflict missing legacy A–E panel-filename wording")
    if "authoritative final Figure 5 structure is A–D" not in conf:
        raise SystemExit("FATAL: routing conflict missing final A–D scientific-structure wording")
    if "Defer asset changes to the final routing phase" not in conf:
        raise SystemExit("FATAL: routing conflict missing deferred migration action")
    # Must not claim the scientific figure itself "is A–E"
    if re.search(r"Figure 5 is A–E|final Figure 5.*A–E structure", conf, re.I):
        raise SystemExit("FATAL: conflict text incorrectly describes final Figure 5 as A–E")

    # Finalize bundle file manifest (all files under OUT)
    # Re-scan to ensure every file appears
    existing = {r["bundle_path"] for r in bundle.rows if r["bundle_path"]}
    for p in sorted(OUT.rglob("*")):
        if not p.is_file():
            continue
        bp = str(p.relative_to(OUT))
        if bp in existing or bp == "01_bundle_file_manifest.tsv":
            continue
        bundle.record(p, None, "bundle_other", notes="present_at_finalize")

    man_df = bundle.to_df()
    # sort deterministically
    man_df = man_df.sort_values(["bundle_path", "original_project_path"], kind="mergesort")
    man_path = OUT / "01_bundle_file_manifest.tsv"
    man_df.to_csv(man_path, sep="\t", index=False)
    # append self-row
    self_row = {
        "bundle_path": "01_bundle_file_manifest.tsv",
        "original_project_path": "",
        "artifact_role": "bundle_file_manifest",
        "figure_id": "",
        "bytes": man_path.stat().st_size,
        "mtime_utc": mtime_utc(man_path),
        "sha256": sha256(man_path),
        "git_tracked": False,
        "git_status": "generated",
        "authoritative": False,
        "notes": "manifest self-entry",
    }
    man_df = pd.concat([man_df, pd.DataFrame([self_row])], ignore_index=True)
    man_df = man_df.sort_values(["bundle_path", "original_project_path"], kind="mergesort")
    man_df.to_csv(man_path, sep="\t", index=False)

    # ZIP with top-level BUNDLE_NAME/
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(OUT.rglob("*")):
            if p.is_file():
                arc = Path(BUNDLE_NAME) / p.relative_to(OUT)
                zf.write(p, arcname=str(arc))

    zip_digest = sha256(ZIP_PATH)
    ZIP_SHA.write_text(f"{zip_digest}  {ZIP_PATH.name}\n", encoding="utf-8")

    file_count = sum(1 for p in OUT.rglob("*") if p.is_file())
    uncompressed = sum(p.stat().st_size for p in OUT.rglob("*") if p.is_file())
    compressed = ZIP_PATH.stat().st_size

    print(f"=== {BUNDLE_NAME} handoff complete ===")
    print(f"bundle_absolute_path\t{OUT}")
    print(f"bundle_relative_path\tresults/handoff/{BUNDLE_NAME}/")
    print(f"zip_absolute_path\t{ZIP_PATH}")
    print(f"bundle_SHA256\t{zip_digest}")
    print(f"file_count\t{file_count}")
    print(f"uncompressed_size\t{uncompressed}")
    print(f"compressed_size\t{compressed}")
    print(f"git_commit\t{commit}")
    print(f"BLOCKER\t{counts.get('BLOCKER', 0)}")
    print(f"MAJOR\t{counts.get('MAJOR', 0)}")
    print(f"MINOR\t{counts.get('MINOR', 0)}")
    print(f"EDITORIAL\t{counts.get('EDITORIAL', 0)}")
    print("STOP — read-only bundle generated; no manuscript/routing/figure edits.")


if __name__ == "__main__":
    main()
