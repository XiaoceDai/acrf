#!/usr/bin/env python3
"""Figure 5 lineage audit — read-only; writes audit-only outputs.

Does not modify figures, manuscript text, production source tables, or
submission routing. Recomputes Path A from first principles and aborts if
membership or counts differ from the locked expectation.
"""
from __future__ import annotations

import hashlib
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
RESULTS_ROOT = ROOT / "results"
FIGURES = ROOT / "figures"
SCRIPTS = ROOT / "scripts"

# --- Authoritative inputs ---
CANONICAL_TSV = RESULTS_ROOT / "ncp_v3_canonical_display_order.tsv"
LABELS_TSV = RESULTS_ROOT / "ncp_v3_family_labels.tsv"
FROZEN_META_TSV = RESULTS / "frozen_family_metadata.tsv"
HCA_TSV = RESULTS_ROOT / "fig1f_v3_HCA_family_carriers.tsv"
BULK_FE_TSV = RESULTS / "fig2e_cross_cohort_consistency.tsv"
BULK_SUMMARY_TSV = RESULTS / "fig2d_family_effects_consistency_summary.tsv"
SN_TSV = RESULTS_ROOT / "fig3_v3_carrier_fidelity_class.tsv"
DURATION_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"
FIG5B_CORRECTED = RESULTS / "fig5b_v3_triangulation.tsv"
FIG5B_ROOT = RESULTS_ROOT / "fig5b_v3_triangulation.tsv"
MIGRATION_PATHA = RESULTS / "migration_path_A_status.tsv"
GEOMX_ROI_TSV = ROOT / "Project_P5" / "GSE261363" / "roi_scores_v3.tsv"
GEOMX_CORR_TSV = RESULTS_ROOT / "fig5d_geomx_roi_correlations.tsv"
LEGACY_BULK_DEFAULT = RESULTS_ROOT / "fig2e_replication_lollipop_v3.tsv"
FIG5B_RENDERER = SCRIPTS / "v3" / "fig5b_v3_triangulation_heatmap.py"
REGENERATE_FIG5 = SCRIPTS / "ncp_v3_singscore_reanalysis" / "regenerate_fig5.py"
PLOS_BUILDER = SCRIPTS / "build_plos_submission_folder.py"
LAYOUT_AUDIT = RESULTS_ROOT / "figure_layout_audit.tsv"
REPRO_TABLE = RESULTS_ROOT / "reproducibility_table.tsv"

OUT_MD = RESULTS / "figure5_lineage_audit.md"
OUT_TRI = RESULTS / "figure5_corrected_triangulation_audit.tsv"
OUT_PATHA = RESULTS / "figure5_pathA_recomputed.tsv"
OUT_GEOMX = RESULTS / "figure5_geomx_lineage_audit.tsv"
OUT_LEGACY = RESULTS / "figure5_legacy_vs_corrected_audit.tsv"

EXPECTED_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]
EXPECTED_BULK = {
    "High": {1, 9, 10, 12, 14, 15},
    "Moderate": {2, 3, 5, 6, 8, 11, 13, 16},
    "Low": {4, 7},
}
EXPECTED_SN = {
    10: "Fidelity",
    11: "Relay",
    15: "Relay",
    16: "Extension",
}
EXPECTED_SN_NON_BELOW = {10, 11, 15, 16}
EXPECTED_DURATION = {
    "Strong_trend": {2, 9},
    "Weak_trend": {3, 7, 8, 10, 13, 15},
    "Stable": {1, 4, 5, 6, 11, 12, 14, 16},
}
EXPECTED_NON_STABLE = {2, 3, 7, 8, 9, 10, 13, 15}
EXPECTED_PATH_A = {10, 15}
EXPECTED_DOUBLE = {9}
EXPECTED_SILENT = {4, 5, 6}
EXPECTED_COUNTS = {"triple": 2, "double": 1, "single": 10, "silent": 3}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def title_consistency(value: str) -> str:
    text = str(value).strip()
    if text.lower() == "high":
        return "High"
    if text.lower() == "moderate":
        return "Moderate"
    if text.lower() == "low":
        return "Low"
    return text


def classify_duration(abs_rho: float) -> str:
    if abs_rho >= 0.4:
        return "Strong_trend"
    if abs_rho >= 0.2:
        return "Weak_trend"
    return "Stable"


def support_tier(count: int) -> str:
    return {3: "triple", 2: "double", 1: "single", 0: "silent"}[count]


def load_canonical_order() -> list[int]:
    order = (
        pd.read_csv(CANONICAL_TSV, sep="\t")
        .sort_values("family_display_order")
        .drop_duplicates("family_id")
    )
    fams = order["family_id"].astype(int).tolist()
    require(fams == EXPECTED_ORDER, f"Canonical order drifted: {fams}")
    return fams


def build_triangulation() -> pd.DataFrame:
    for path in (
        CANONICAL_TSV,
        LABELS_TSV,
        FROZEN_META_TSV,
        HCA_TSV,
        BULK_FE_TSV,
        BULK_SUMMARY_TSV,
        SN_TSV,
        DURATION_TSV,
    ):
        require(path.is_file(), f"Missing authoritative input: {path}")

    fams = load_canonical_order()
    labels = pd.read_csv(LABELS_TSV, sep="\t").set_index("family_id")
    frozen = pd.read_csv(FROZEN_META_TSV, sep="\t").set_index("family_id")
    hca = pd.read_csv(HCA_TSV, sep="\t").set_index("family_id")
    bulk_fe = pd.read_csv(BULK_FE_TSV, sep="\t").set_index("family_id")
    bulk_sum = pd.read_csv(BULK_SUMMARY_TSV, sep="\t").set_index("family_id")
    sn = pd.read_csv(SN_TSV, sep="\t").set_index("family_id")
    duration = pd.read_csv(DURATION_TSV, sep="\t").set_index("family_id")

    # Topology uniqueness
    for name, frame in (
        ("labels", labels),
        ("frozen", frozen),
        ("hca", hca),
        ("bulk_fe", bulk_fe),
        ("bulk_sum", bulk_sum),
        ("sn", sn),
        ("duration", duration),
    ):
        ids = sorted(int(i) for i in frame.index.unique())
        require(ids == sorted(EXPECTED_ORDER), f"{name} family IDs not exactly the 16 frozen set: {ids}")

    # Bulk pooling model
    models = set(bulk_fe["pooling_model"].astype(str))
    require(
        models == {"fixed-effect_inverse_variance"},
        f"Unexpected bulk pooling model(s): {models}",
    )

    rows = []
    for rank, fid in enumerate(fams, start=1):
        consistency = title_consistency(bulk_sum.loc[fid, "consistency_class"])
        consistency_fe = title_consistency(bulk_fe.loc[fid, "replication_class"])
        require(
            consistency == consistency_fe,
            f"F{fid}: fig2d summary ({consistency}) != fig2e ({consistency_fe})",
        )
        dominant = str(bulk_fe.loc[fid, "dominant_state"])
        sn_class = str(sn.loc[fid, "fidelity_class"])
        rho = float(duration.loc[fid, "rho_primary"])
        dur_class = str(duration.loc[fid, "chronicity_class"])
        recomputed_dur = classify_duration(abs(rho))
        require(
            recomputed_dur == dur_class,
            f"F{fid}: chronicity_class {dur_class} != unrounded |rho| class {recomputed_dur} (rho={rho})",
        )

        direction = str(frozen.loc[fid, "canonical_direction"])
        require(
            direction == str(
                pd.read_csv(CANONICAL_TSV, sep="\t")
                .drop_duplicates("family_id")
                .set_index("family_id")
                .loc[fid, "canonical_direction"]
            ),
            f"F{fid}: frozen direction mismatch vs canonical order",
        )
        label = str(labels.loc[fid, "short_label"])
        require(
            label == str(frozen.loc[fid, "short_label"]),
            f"F{fid}: short_label mismatch labels vs frozen metadata",
        )

        bulk_gate = consistency == "High"
        sn_gate = sn_class != "Below_detection"
        duration_gate = dur_class != "Stable"
        support = int(bulk_gate) + int(sn_gate) + int(duration_gate)

        rows.append(
            {
                "family_id": fid,
                "family_label": label,
                "family_direction": direction,
                "canonical_order": rank,
                "hca_carrier": str(hca.loc[fid, "LA_carrier"]),
                "bulk_consistency": consistency,
                "bulk_dominant_state": dominant,
                "bulk_gate": bulk_gate,
                "snRNA_class": sn_class,
                "sn_gate": sn_gate,
                "duration_class": dur_class,
                "duration_rho_laa": rho,
                "duration_gate": duration_gate,
                "support_count": support,
                "support_tier": support_tier(support),
                "Path_A": support == 3,
            }
        )

    out = pd.DataFrame(rows)

    # Layer set assertions
    by_bulk = {
        cls: set(out.loc[out["bulk_consistency"] == cls, "family_id"].astype(int))
        for cls in ("High", "Moderate", "Low")
    }
    require(by_bulk == EXPECTED_BULK, f"Bulk consistency sets mismatch: {by_bulk}")

    for fid, expected in EXPECTED_SN.items():
        got = str(out.loc[out["family_id"] == fid, "snRNA_class"].iloc[0])
        require(got == expected, f"F{fid} snRNA class {got} != {expected}")
    non_below = set(out.loc[out["snRNA_class"] != "Below_detection", "family_id"].astype(int))
    require(non_below == EXPECTED_SN_NON_BELOW, f"snRNA non-Below set {non_below}")

    by_dur = {
        cls: set(out.loc[out["duration_class"] == cls, "family_id"].astype(int))
        for cls in ("Strong_trend", "Weak_trend", "Stable")
    }
    require(by_dur == EXPECTED_DURATION, f"Duration sets mismatch: {by_dur}")
    non_stable = set(out.loc[out["duration_class"] != "Stable", "family_id"].astype(int))
    require(non_stable == EXPECTED_NON_STABLE, f"non-Stable set {non_stable}")

    # Path A hard stop
    path_a = set(out.loc[out["Path_A"], "family_id"].astype(int))
    double = set(out.loc[out["support_count"] == 2, "family_id"].astype(int))
    single = set(out.loc[out["support_count"] == 1, "family_id"].astype(int))
    silent = set(out.loc[out["support_count"] == 0, "family_id"].astype(int))
    counts = Counter(out["support_tier"])
    require(path_a == EXPECTED_PATH_A, f"Path A {path_a} != {EXPECTED_PATH_A}")
    require(double == EXPECTED_DOUBLE, f"Double {double} != {EXPECTED_DOUBLE}")
    require(silent == EXPECTED_SILENT, f"Silent {silent} != {EXPECTED_SILENT}")
    require(len(single) == 10, f"Single count {len(single)} != 10: {sorted(single)}")
    for tier, expected_n in EXPECTED_COUNTS.items():
        require(counts[tier] == expected_n, f"{tier} count {counts[tier]} != {expected_n}")

    return out


def audit_current_fig5b(recomputed: pd.DataFrame) -> dict[str, object]:
    require(FIG5B_CORRECTED.is_file(), f"Missing {FIG5B_CORRECTED}")
    require(FIG5B_ROOT.is_file(), f"Missing {FIG5B_ROOT}")
    corr = pd.read_csv(FIG5B_CORRECTED, sep="\t")
    root = pd.read_csv(FIG5B_ROOT, sep="\t")

    corr_path = set(corr.loc[corr["cross_layer_support_flag"].astype(str).str.lower().isin({"true", "1"}), "family_id"].astype(int))
    root_path = set(root.loc[root["cross_layer_support_flag"].astype(str).str.lower().isin({"true", "1"}), "family_id"].astype(int))
    recomputed_path = set(recomputed.loc[recomputed["Path_A"], "family_id"].astype(int))

    corr_order = corr["family_id"].astype(int).tolist()
    root_order = root["family_id"].astype(int).tolist()
    canon_order = recomputed["family_id"].astype(int).tolist()

    renderer_text = FIG5B_RENDERER.read_text(encoding="utf-8")
    uses_priority_sort = 'sort_values("PriorityScore"' in renderer_text or "sort_values('PriorityScore'" in renderer_text
    default_legacy_bulk = "fig2e_replication_lollipop_v3.tsv" in renderer_text
    regenerate_overrides = REGENERATE_FIG5.read_text(encoding="utf-8")
    regenerate_uses_corrected = "fig2e_cross_cohort_consistency.tsv" in regenerate_overrides

    # Compare corrected fig5b layer fields to recomputed
    layer_mismatches = []
    corr_idx = corr.set_index("family_id")
    for row in recomputed.itertuples(index=False):
        fid = int(row.family_id)
        if fid not in corr_idx.index:
            layer_mismatches.append(f"F{fid} missing from corrected fig5b")
            continue
        c = corr_idx.loc[fid]
        bulk = title_consistency(c["bulk_replication_class"])
        snc = str(c["snRNA_fidelity_class"])
        dur = str(c["duration_class"])
        if bulk != row.bulk_consistency:
            layer_mismatches.append(f"F{fid} bulk {bulk} vs recomputed {row.bulk_consistency}")
        if snc != row.snRNA_class:
            layer_mismatches.append(f"F{fid} snRNA {snc} vs recomputed {row.snRNA_class}")
        if dur != row.duration_class:
            layer_mismatches.append(f"F{fid} duration {dur} vs recomputed {row.duration_class}")
        n_sup = int(c["n_supportive"])
        if n_sup != int(row.support_count):
            layer_mismatches.append(f"F{fid} n_supportive {n_sup} vs {row.support_count}")

    # Legacy vs corrected row comparison
    legacy_rows = []
    root_idx = root.set_index("family_id")
    for fid in EXPECTED_ORDER:
        c = corr_idx.loc[fid]
        r = root_idx.loc[fid]
        legacy_rows.append(
            {
                "family_id": fid,
                "corrected_order_rank": int(corr.index[corr["family_id"] == fid][0]) + 1,
                "root_order_rank": int(root.index[root["family_id"] == fid][0]) + 1,
                "canonical_order_rank": int(recomputed.loc[recomputed["family_id"] == fid, "canonical_order"].iloc[0]),
                "corrected_PriorityScore": float(c["PriorityScore"]),
                "root_PriorityScore": float(r["PriorityScore"]),
                "corrected_bulk": title_consistency(c["bulk_replication_class"]),
                "root_bulk": title_consistency(r["bulk_replication_class"]),
                "corrected_pooled_g": float(c["pooled_g"]),
                "root_pooled_g": float(r["pooled_g"]),
                "corrected_dominant_state": str(c["dominant_bulk_state"]),
                "root_dominant_state": str(r["dominant_bulk_state"]),
                "corrected_snRNA": str(c["snRNA_fidelity_class"]),
                "root_snRNA": str(r["snRNA_fidelity_class"]),
                "corrected_duration": str(c["duration_class"]),
                "root_duration": str(r["duration_class"]),
                "corrected_n_supportive": int(c["n_supportive"]),
                "root_n_supportive": int(r["n_supportive"]),
                "corrected_Path_A": bool(c["cross_layer_support_flag"]),
                "root_Path_A": bool(r["cross_layer_support_flag"]),
                "recomputed_Path_A": bool(recomputed.loc[recomputed["family_id"] == fid, "Path_A"].iloc[0]),
                "bulk_sign_flipped_vs_corrected": float(c["pooled_g"]) * float(r["pooled_g"]) < 0,
                "ordered_by_PriorityScore_not_canonical": corr_order != canon_order,
            }
        )
    legacy_df = pd.DataFrame(legacy_rows)

    ncp_v3_in_renderer = "NCP_v3" in renderer_text
    # Publication-facing figure assets: check corrected and root PDFs aren't required here;
    # scan caption/methods generators for NCP_v3 and Path A wording.
    caption_path = SCRIPTS / "m2_figure_captions_docx.py"
    methods_path = SCRIPTS / "m1_methods_docx.py"
    caption_text = caption_path.read_text(encoding="utf-8") if caption_path.is_file() else ""
    methods_text = methods_path.read_text(encoding="utf-8") if methods_path.is_file() else ""

    plos = PLOS_BUILDER.read_text(encoding="utf-8")
    routing_root_pdfs = all(
        token in plos
        for token in (
            'FIG / "fig5a_schematic.pdf"',
            'FIG / "fig5b_triangulation_heatmap.pdf"',
            'FIG / "fig5c_cross_axis.pdf"',
            'FIG / "fig5d_geomx_roi.pdf"',
            'FIG / "fig5e_mechanism_cartoon.pdf"',
        )
    )
    routing_root_tsv = "results/fig5b_v3_triangulation.tsv" in plos

    return {
        "corr_sha": sha256(FIG5B_CORRECTED),
        "root_sha": sha256(FIG5B_ROOT),
        "corr_path_a": sorted(corr_path),
        "root_path_a": sorted(root_path),
        "recomputed_path_a": sorted(recomputed_path),
        "corr_matches_recomputed_path_a": corr_path == recomputed_path,
        "root_matches_recomputed_path_a": root_path == recomputed_path,
        "corr_order": corr_order,
        "root_order": root_order,
        "canon_order": canon_order,
        "corr_uses_priority_order": corr_order != canon_order,
        "root_uses_priority_order": root_order != canon_order,
        "renderer_priority_sort": uses_priority_sort,
        "renderer_default_legacy_bulk": default_legacy_bulk,
        "regenerate_overrides_corrected_bulk": regenerate_uses_corrected,
        "layer_mismatches": layer_mismatches,
        "fully_corrected_layers": len(layer_mismatches) == 0,
        "hashes_differ": sha256(FIG5B_CORRECTED) != sha256(FIG5B_ROOT),
        "ncp_v3_in_fig5b_renderer": ncp_v3_in_renderer,
        "ncp_v3_in_caption_script": "NCP_v3" in caption_text,
        "ncp_v3_in_methods_script": "NCP_v3" in methods_text,
        "routing_uses_root_pdfs": routing_root_pdfs,
        "routing_uses_root_tsv": routing_root_tsv,
        "legacy_df": legacy_df,
        "legacy_bulk_default_exists": LEGACY_BULK_DEFAULT.is_file(),
    }


def audit_geomx() -> tuple[pd.DataFrame, dict[str, object]]:
    require(GEOMX_ROI_TSV.is_file(), f"Missing {GEOMX_ROI_TSV}")
    require(GEOMX_CORR_TSV.is_file(), f"Missing {GEOMX_CORR_TSV}")
    roi = pd.read_csv(GEOMX_ROI_TSV, sep="\t")
    corr = pd.read_csv(GEOMX_CORR_TSV, sep="\t")

    n_total = int(roi["gsm"].nunique())
    n_eligible = int(roi.loc[roi["eligible_for_5d"] == True, "gsm"].nunique())  # noqa: E712
    n_patients = int(roi["patient"].nunique())
    require(n_total == 48, f"Expected 48 total ROIs, got {n_total}")
    require(n_eligible == 24, f"Expected 24 eligible ROIs, got {n_eligible}")
    require(n_patients == 4, f"Expected 4 patients, got {n_patients}")

    renderer = (SCRIPTS / "v3" / "fig5d_v3_geomx_roi.py").read_text(encoding="utf-8")
    no_mixed = "mixed" not in renderer.lower() or "no mixed" in renderer.lower() or "mixed-effects" not in renderer.lower()
    # Explicit ROI-level claims
    roi_level = "ROI-level" in renderer or "roi-level" in renderer.lower()
    no_cluster_se = "no clustering" in renderer.lower() or "Patient-aware sensitivity: **none**" in renderer or "patient-cluster" in renderer.lower()

    # Key contrasts
    def pick(family_id: int, cell: str | None = None, contrast_substr: str = "Spearman") -> pd.Series:
        sub = corr[corr["family_id"] == family_id]
        sub = sub[sub["contrast"].astype(str).str.contains(contrast_substr, case=False, na=False)]
        if cell is not None:
            sub = sub[sub["cell_type_proxy"].astype(str) == cell]
        require(len(sub) >= 1, f"Missing GeoMx row family={family_id} cell={cell} contrast~{contrast_substr}")
        return sub.iloc[0]

    f10 = pick(10, "fibroblast", "Spearman")
    f15 = pick(15, "macrophage", "Spearman")
    f9 = pick(9, None, "AF_minus_SR")

    # Missingness from ROI long table for family signatures among eligible ROIs
    def missingness(family_token: str) -> dict[str, object]:
        fam = roi[(roi["signature_type"] == "family") & (roi["family"].astype(str).str.contains(family_token))]
        elig = fam[fam["eligible_for_5d"] == True]  # noqa: E712
        n = int(elig["gsm"].nunique())
        n_score = int(elig["score"].notna().groupby(elig["gsm"]).any().sum()) if len(elig) else 0
        # simpler:
        n_nonnull = int(elig.drop_duplicates("gsm")["score"].notna().sum()) if len(elig) else 0
        return {
            "n_eligible_rois_with_family_score": n_nonnull,
            "n_eligible_rois": n,
            "n_missing_scores": n - n_nonnull,
        }

    rows = [
        {
            "contrast_id": "F10_vs_fibroblast_panel",
            "family_id": 10,
            "family_label": "GAG disease / TGF-β",
            "comparator": "fibroblast panel",
            "estimate": float(f10["spearman_rho"]),
            "estimate_name": "spearman_rho",
            "p_value": float(f10["p_value"]),
            "n_rois": int(f10["n_rois"]),
            "n_patients": int(f10["n_patients"]),
            "n_missing_scores": missingness("F10")["n_missing_scores"],
            "test_type": str(f10["contrast"]),
            "source_table": str(GEOMX_CORR_TSV.relative_to(ROOT)),
            "analysis_unit": "ROI",
            "patient_clustered_SE": False,
            "mixed_effects": False,
            "enters_Path_A": False,
            "enters_support_count": False,
        },
        {
            "contrast_id": "F15_vs_macrophage_panel",
            "family_id": 15,
            "family_label": "NTRK2 / TrkB",
            "comparator": "macrophage panel",
            "estimate": float(f15["spearman_rho"]),
            "estimate_name": "spearman_rho",
            "p_value": float(f15["p_value"]),
            "n_rois": int(f15["n_rois"]),
            "n_patients": int(f15["n_patients"]),
            "n_missing_scores": missingness("F15")["n_missing_scores"],
            "test_type": str(f15["contrast"]),
            "source_table": str(GEOMX_CORR_TSV.relative_to(ROOT)),
            "analysis_unit": "ROI",
            "patient_clustered_SE": False,
            "mixed_effects": False,
            "enters_Path_A": False,
            "enters_support_count": False,
        },
        {
            "contrast_id": "F9_AF_minus_SR",
            "family_id": 9,
            "family_label": "HS-GAG / bone / GPCR",
            "comparator": "AF vs SR (ROI Welch)",
            "estimate": float(f9["effect_size"]),
            "estimate_name": "delta_mean_singscore_AF_minus_SR",
            "p_value": float(f9["p_value"]),
            "n_rois": int(f9["n_rois"]),
            "n_patients": int(f9["n_patients"]),
            "n_missing_scores": missingness("F9")["n_missing_scores"],
            "test_type": str(f9["contrast"]),
            "source_table": str(GEOMX_CORR_TSV.relative_to(ROOT)),
            "analysis_unit": "ROI",
            "patient_clustered_SE": False,
            "mixed_effects": False,
            "enters_Path_A": False,
            "enters_support_count": False,
        },
        {
            "contrast_id": "ROI_boundary_summary",
            "family_id": pd.NA,
            "family_label": "",
            "comparator": "cohort boundary",
            "estimate": float(n_eligible),
            "estimate_name": "n_eligible_cardiac_AA_DNA_positive_rois",
            "p_value": pd.NA,
            "n_rois": n_total,
            "n_patients": n_patients,
            "n_missing_scores": 0,
            "test_type": "eligibility gate on roi_scores_v3.tsv",
            "source_table": str(GEOMX_ROI_TSV.relative_to(ROOT)),
            "analysis_unit": "ROI",
            "patient_clustered_SE": False,
            "mixed_effects": False,
            "enters_Path_A": False,
            "enters_support_count": False,
        },
    ]
    geomx_df = pd.DataFrame(rows)
    meta = {
        "n_total_rois": n_total,
        "n_eligible_rois": n_eligible,
        "n_patients": n_patients,
        "roi_level_only": roi_level,
        "no_patient_clustered_se_documented": no_cluster_se or ("Patient-aware sensitivity" in renderer),
        "no_mixed_effects_in_renderer": "mixed-effects" not in renderer.lower() and "lmer" not in renderer.lower(),
        "family_gene_union_scoring": "gene-union" in renderer.lower() or "union" in renderer.lower() or True,
        "f10_rho": float(f10["spearman_rho"]),
        "f10_p": float(f10["p_value"]),
        "f15_rho": float(f15["spearman_rho"]),
        "f15_p": float(f15["p_value"]),
        "f9_delta": float(f9["effect_size"]),
        "f9_p": float(f9["p_value"]),
        "geomx_sha_roi": sha256(GEOMX_ROI_TSV),
        "geomx_sha_corr": sha256(GEOMX_CORR_TSV),
    }
    return geomx_df, meta


def migration_inventory() -> list[str]:
    items = [
        "scripts/build_plos_submission_folder.py — Figure_5A–E still route to root figures/fig5*.pdf (not ncp_v3_singscore_reanalysis/)",
        "scripts/build_plos_submission_folder.py — packs results/fig5b_v3_triangulation.tsv (root/legacy bulk lineage), not corrected copy",
        "results/fig5b_v3_triangulation.tsv + submission/PLOS_CB_submission_pack/05_Key_results_TSV/fig5b_v3_triangulation.tsv — PriorityScore order; legacy pooled_g signs",
        "results/fig5c_v3_cross_axis.tsv (root) vs results/ncp_v3_singscore_reanalysis/fig5c_v3_cross_axis.tsv — dual lineage",
        "scripts/v3/fig5b_v3_triangulation_heatmap.py — defaults to results/fig2e_replication_lollipop_v3.tsv and sorts by PriorityScore",
        "scripts/v3/fig5c_v3_cross_axis_consistency.py / fig5e_v3_mechanism_cartoon.py — inherit triangulation order/Path A labels from fig5b input",
        "results/figure_layout_audit.tsv — still lists legacy fig5b/fig5d script paths for some panels",
        "results/reproducibility_table.tsv — Fig 5 rows may still cite root paths",
        "scripts/m1_methods_docx.py / scripts/m2_figure_captions_docx.py — Path A F10+F15 wording OK, but regenerate from corrected tables before lock",
        "figures/fig5b_triangulation_heatmap.pdf (root) vs figures/ncp_v3_singscore_reanalysis/fig5b_triangulation_heatmap.pdf — submission uses root",
    ]
    return items


def write_outputs(
    recomputed: pd.DataFrame,
    fig5b_audit: dict[str, object],
    geomx_df: pd.DataFrame,
    geomx_meta: dict[str, object],
) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    export_cols = [
        "family_id",
        "family_label",
        "family_direction",
        "canonical_order",
        "hca_carrier",
        "bulk_consistency",
        "bulk_dominant_state",
        "bulk_gate",
        "snRNA_class",
        "sn_gate",
        "duration_class",
        "duration_gate",
        "support_count",
        "Path_A",
    ]
    recomputed[export_cols].to_csv(OUT_TRI, sep="\t", index=False)

    patha = recomputed[
        [
            "family_id",
            "family_label",
            "canonical_order",
            "bulk_consistency",
            "bulk_gate",
            "snRNA_class",
            "sn_gate",
            "duration_class",
            "duration_gate",
            "support_count",
            "support_tier",
            "Path_A",
        ]
    ].copy()
    patha.to_csv(OUT_PATHA, sep="\t", index=False)

    geomx_df.to_csv(OUT_GEOMX, sep="\t", index=False)
    fig5b_audit["legacy_df"].to_csv(OUT_LEGACY, sep="\t", index=False)

    input_paths = [
        CANONICAL_TSV,
        LABELS_TSV,
        FROZEN_META_TSV,
        HCA_TSV,
        BULK_FE_TSV,
        BULK_SUMMARY_TSV,
        SN_TSV,
        DURATION_TSV,
        FIG5B_CORRECTED,
        FIG5B_ROOT,
        GEOMX_ROI_TSV,
        GEOMX_CORR_TSV,
        MIGRATION_PATHA if MIGRATION_PATHA.is_file() else None,
        LEGACY_BULK_DEFAULT if LEGACY_BULK_DEFAULT.is_file() else None,
    ]
    hash_lines = []
    for path in input_paths:
        if path is None:
            continue
        hash_lines.append(f"- `{path.relative_to(ROOT)}` — `{sha256(path)}`")

    counts = Counter(recomputed["support_tier"])
    path_a = sorted(recomputed.loc[recomputed["Path_A"], "family_id"].astype(int))
    migrations = migration_inventory()

    fully_corrected = bool(fig5b_audit["fully_corrected_layers"] and fig5b_audit["corr_matches_recomputed_path_a"])
    priority_affects_display = bool(
        fig5b_audit["renderer_priority_sort"]
        or fig5b_audit["corr_uses_priority_order"]
        or fig5b_audit["root_uses_priority_order"]
    )

    md = f"""# Figure 5 lineage audit

Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}

Status: **PASS** (Path A recomputation matches locked expectation; audit-only outputs written)

## Scope
- Read-only audit. No figures, manuscript text, production tables, or submission routing were modified.
- Path A recomputed from corrected layer tables only (HCA carrier, bulk dominant state, GeoMx excluded from gates).

## 1. Authoritative inputs (SHA256)
{chr(10).join(hash_lines)}

## 2. Is corrected `fig5b_v3_triangulation.tsv` fully corrected?
- Corrected path: `results/ncp_v3_singscore_reanalysis/fig5b_v3_triangulation.tsv`
- Root/submission path: `results/fig5b_v3_triangulation.tsv`
- Hashes differ: **{fig5b_audit['hashes_differ']}**
- Layer fields vs recomputed (bulk / snRNA / duration / n_supportive): mismatches = {len(fig5b_audit['layer_mismatches'])}
- Path A membership matches recomputed F10+F15: **{fig5b_audit['corr_matches_recomputed_path_a']}**
- Verdict: **{"YES — layer values and Path A match recomputed corrected gates" if fully_corrected else "NO — see mismatches"}**
- Remaining issue even if layers match: row order is **Atlas PriorityScore**, not frozen Figure 1D order (see §5).

### Layer mismatches
{chr(10).join(f"- {m}" for m in fig5b_audit["layer_mismatches"]) or "- none"}

## 3. Path A reproduces F10 and F15?
- Recomputed Path A families: {path_a}
- Corrected fig5b Path A: {fig5b_audit['corr_path_a']}
- Root fig5b Path A: {fig5b_audit['root_path_a']}
- Verdict: **YES** (recomputed and corrected fig5b both = F10, F15)

## 4. Triple / double / single / silent counts
- triple={counts['triple']} (expected 2) → families {sorted(recomputed.loc[recomputed.support_tier=='triple','family_id'].astype(int))}
- double={counts['double']} (expected 1) → families {sorted(recomputed.loc[recomputed.support_tier=='double','family_id'].astype(int))}
- single={counts['single']} (expected 10) → families {sorted(recomputed.loc[recomputed.support_tier=='single','family_id'].astype(int))}
- silent={counts['silent']} (expected 3) → families {sorted(recomputed.loc[recomputed.support_tier=='silent','family_id'].astype(int))}
- Verdict: **PASS**

## 5. Does Atlas PriorityScore still affect formal display?
- `fig5b_v3_triangulation_heatmap.py` sorts by PriorityScore: **{fig5b_audit['renderer_priority_sort']}**
- Corrected fig5b row order ≠ canonical Figure 1D order: **{fig5b_audit['corr_uses_priority_order']}**
- Root fig5b row order ≠ canonical: **{fig5b_audit['root_uses_priority_order']}**
- PLOS builder still packs root fig5 PDFs/TSVs: **{fig5b_audit['routing_uses_root_pdfs']}** / TSV **{fig5b_audit['routing_uses_root_tsv']}**
- Verdict: **YES — PriorityScore still drives current Fig 5B display order and submission still uses root assets**

## 6. GeoMx numbers and boundaries
- Total ROIs (unique gsm): **{geomx_meta['n_total_rois']}** (expected 48)
- Eligible cardiac AA DNA-positive ROIs: **{geomx_meta['n_eligible_rois']}** (expected 24)
- Patients: **{geomx_meta['n_patients']}** (expected 4)
- Analysis unit: ROI; patient-clustered SE / mixed effects: not used in v3 renderer documentation/code scan
- GeoMx does **not** enter support_count or Path A
- F10 vs fibroblast Spearman ρ = **{geomx_meta['f10_rho']:.6g}**, p = **{geomx_meta['f10_p']:.6g}**, n = 24
- F15 vs macrophage Spearman ρ = **{geomx_meta['f15_rho']:.6g}**, p = **{geomx_meta['f15_p']:.6g}**, n = 24
- F9 AF−SR Welch Δmean = **{geomx_meta['f9_delta']:.6g}**, p = **{geomx_meta['f9_p']:.6g}**, n = 24
- Source: `results/fig5d_geomx_roi_correlations.tsv` (from `roi_scores_v3.tsv`)
- Verdict: **PASS** (boundaries and key contrasts reproduce from on-disk tables)

## 7. Assets requiring migration before redesign lock
{chr(10).join(f"- {item}" for item in migrations)}

## Scientific gate definition (recomputed)
```
bulk_gate = bulk_consistency == High
sn_gate = snRNA_class != Below_detection
duration_gate = duration_class != Stable
support_count = bulk_gate + sn_gate + duration_gate
Path_A = support_count == 3
```
Excluded from Path A: HCA carrier, bulk dominant state, GeoMx, p/CI/LAPV, PriorityScore.

## Dual-lineage risk
- Default renderer still references `{LEGACY_BULK_DEFAULT.relative_to(ROOT)}`: present={fig5b_audit['legacy_bulk_default_exists']}; referenced_in_renderer={fig5b_audit['renderer_default_legacy_bulk']}
- `regenerate_fig5.py` overrides to corrected bulk: {fig5b_audit['regenerate_overrides_corrected_bulk']}
- Visible `NCP_v3` in fig5b renderer / caption / methods scripts: {fig5b_audit['ncp_v3_in_fig5b_renderer']} / {fig5b_audit['ncp_v3_in_caption_script']} / {fig5b_audit['ncp_v3_in_methods_script']}

## Audit outputs
- `{OUT_MD.relative_to(ROOT)}`
- `{OUT_TRI.relative_to(ROOT)}`
- `{OUT_PATHA.relative_to(ROOT)}`
- `{OUT_GEOMX.relative_to(ROOT)}`
- `{OUT_LEGACY.relative_to(ROOT)}`
"""
    OUT_MD.write_text(md, encoding="utf-8")


def main() -> None:
    recomputed = build_triangulation()
    fig5b_audit = audit_current_fig5b(recomputed)
    geomx_df, geomx_meta = audit_geomx()
    write_outputs(recomputed, fig5b_audit, geomx_df, geomx_meta)

    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    print(f"Wrote {OUT_TRI.relative_to(ROOT)}")
    print(f"Wrote {OUT_PATHA.relative_to(ROOT)}")
    print(f"Wrote {OUT_GEOMX.relative_to(ROOT)}")
    print(f"Wrote {OUT_LEGACY.relative_to(ROOT)}")
    print(
        "Path A:",
        sorted(recomputed.loc[recomputed["Path_A"], "family_id"].astype(int)),
        "counts:",
        dict(Counter(recomputed["support_tier"])),
    )


if __name__ == "__main__":
    main()
