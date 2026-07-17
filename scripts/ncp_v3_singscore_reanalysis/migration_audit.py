#!/usr/bin/env python3
"""Generate old-versus-corrected tables and audit legacy consumers."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
SCRIPT_ROOT = ROOT / "scripts"


def write_m_comparison() -> None:
    old = pd.read_csv(ROOT / "results" / "fig1b_baseline_donor_aCAI.tsv", sep="\t")
    new = pd.read_csv(RESULTS / "healthy_reference_program_summary.tsv", sep="\t")
    old = old[["pathway", "mean_aCAI", "sign_concordance"]].rename(
        columns={
            "pathway": "program_id",
            "mean_aCAI": "legacy_mean_delta",
            "sign_concordance": "legacy_sign_concordance",
        }
    )
    old["legacy_M_abs_mean_delta"] = old["legacy_mean_delta"].abs()
    merged = new.merge(old, on="program_id", how="outer", indicator=True)
    merged["delta_corrected_minus_legacy_M"] = (
        merged["M"] - merged["legacy_M_abs_mean_delta"]
    )
    merged["concordance_changed"] = (
        merged["directional_concordance_k"] != merged["legacy_sign_concordance"]
    )
    merged.to_csv(RESULTS / "migration_M_program.tsv", sep="\t", index=False)


def write_donor_comparisons() -> None:
    old_program = pd.read_csv(
        ROOT / "results" / "fig2b_program_donor_acairel.tsv", sep="\t"
    ).rename(columns={"cohort": "dataset", "aCAIrel": "legacy_aCAIrel"})
    new_program = pd.read_csv(
        RESULTS / "disease_program_contrasts_long.tsv", sep="\t"
    ).rename(columns={"aCAIrel": "corrected_aCAIrel"})
    key = ["dataset", "donor_id", "rhythm", "program_id"]
    program = new_program[key + ["family_id", "corrected_aCAIrel"]].merge(
        old_program[key + ["legacy_aCAIrel"]],
        on=key,
        how="outer",
        indicator=True,
    )
    program["delta_corrected_minus_legacy"] = (
        program["corrected_aCAIrel"] - program["legacy_aCAIrel"]
    )
    program.to_csv(
        RESULTS / "migration_program_donor_acairel.tsv", sep="\t", index=False
    )

    old_family = pd.read_csv(
        ROOT / "results" / "fig2b_family_donor_acairel.tsv", sep="\t"
    ).rename(columns={"cohort": "dataset", "aCAIrel": "legacy_family_aCAIrel"})
    new_family = pd.read_csv(RESULTS / "family_acairel_long.tsv", sep="\t")
    new_family = new_family[new_family["dataset"] != "GSE112339"].rename(
        columns={"family_aCAIrel": "corrected_family_aCAIrel"}
    )
    family_key = ["dataset", "donor_id", "rhythm", "family_id"]
    family = new_family[family_key + ["corrected_family_aCAIrel"]].merge(
        old_family[family_key + ["legacy_family_aCAIrel"]],
        on=family_key,
        how="outer",
        indicator=True,
    )
    family["delta_corrected_minus_legacy"] = (
        family["corrected_family_aCAIrel"] - family["legacy_family_aCAIrel"]
    )
    family.to_csv(
        RESULTS / "migration_family_donor_acairel.tsv", sep="\t", index=False
    )


def write_classification_comparisons() -> pd.DataFrame:
    old_states = pd.read_csv(
        ROOT / "results" / "fig2c_chamber_decomposition_v3.tsv", sep="\t"
    )[["family_id", "cohort", "state"]].rename(columns={"state": "legacy_state"})
    new_states = pd.read_csv(
        RESULTS / "fig2c_family_cohort_states.tsv", sep="\t"
    )[["family_id", "cohort", "state"]].rename(columns={"state": "corrected_state"})
    states = new_states.merge(
        old_states, on=["family_id", "cohort"], how="outer", indicator=True
    )
    states["classification_changed"] = (
        states["corrected_state"] != states["legacy_state"]
    )
    states.to_csv(
        RESULTS / "migration_family_cohort_states.tsv", sep="\t", index=False
    )

    old_consistency = pd.read_csv(
        ROOT / "results" / "fig2e_replication_lollipop_v3.tsv", sep="\t"
    )
    new_consistency = pd.read_csv(
        RESULTS / "fig2e_cross_cohort_consistency.tsv", sep="\t"
    )
    keep = [
        "family_id", "pooled_g", "pooled_ci_lo", "pooled_ci_hi",
        "replication_class", "dominant_state", "state_concordance_n",
    ]
    pooled = new_consistency[keep].merge(
        old_consistency[keep],
        on="family_id",
        how="outer",
        suffixes=("_corrected", "_legacy"),
        indicator=True,
    )
    pooled["replication_class_changed"] = (
        pooled["replication_class_corrected"] != pooled["replication_class_legacy"]
    )
    pooled["dominant_state_changed"] = (
        pooled["dominant_state_corrected"] != pooled["dominant_state_legacy"]
    )
    pooled.to_csv(
        RESULTS / "migration_pooled_effects_and_consistency.tsv",
        sep="\t",
        index=False,
    )

    old_path = pd.read_csv(ROOT / "results" / "fig5b_v3_triangulation.tsv", sep="\t")
    new_path = pd.read_csv(RESULTS / "fig5b_v3_triangulation.tsv", sep="\t")
    path_keep = [
        "family_id", "L1_supportive", "L2_supportive", "L3_supportive",
        "n_supportive", "cross_layer_support_flag",
    ]
    path_a = new_path[path_keep].merge(
        old_path[path_keep],
        on="family_id",
        how="outer",
        suffixes=("_corrected", "_legacy"),
        indicator=True,
    )
    old_cross = pd.read_csv(ROOT / "results" / "fig5c_v3_cross_axis.tsv", sep="\t")
    new_cross = pd.read_csv(RESULTS / "fig5c_v3_cross_axis.tsv", sep="\t")
    path_a = path_a.merge(
        new_cross[["family_id", "support_pattern", "cross_axis_class"]].rename(
            columns={
                "support_pattern": "support_pattern_corrected",
                "cross_axis_class": "cross_axis_class_corrected",
            }
        ),
        on="family_id",
        how="left",
    ).merge(
        old_cross[["family_id", "support_pattern", "cross_axis_class"]].rename(
            columns={
                "support_pattern": "support_pattern_legacy",
                "cross_axis_class": "cross_axis_class_legacy",
            }
        ),
        on="family_id",
        how="left",
    )
    path_a["path_A_changed"] = (
        path_a["support_pattern_corrected"] != path_a["support_pattern_legacy"]
    ) | (
        path_a["cross_axis_class_corrected"] != path_a["cross_axis_class_legacy"]
    )
    path_a.to_csv(RESULTS / "migration_path_A_status.tsv", sep="\t", index=False)

    changes = []
    for row in states[states["classification_changed"]].itertuples():
        changes.append(
            {
                "classification_type": "family_cohort_state",
                "family_id": row.family_id,
                "cohort": row.cohort,
                "legacy": row.legacy_state,
                "corrected": row.corrected_state,
            }
        )
    for row in pooled[pooled["replication_class_changed"]].itertuples():
        changes.append(
            {
                "classification_type": "replication_class",
                "family_id": row.family_id,
                "cohort": "POOLED",
                "legacy": row.replication_class_legacy,
                "corrected": row.replication_class_corrected,
            }
        )
    for row in pooled[pooled["dominant_state_changed"]].itertuples():
        changes.append(
            {
                "classification_type": "dominant_state",
                "family_id": row.family_id,
                "cohort": "cross_cohort",
                "legacy": row.dominant_state_legacy,
                "corrected": row.dominant_state_corrected,
            }
        )
    for row in path_a[path_a["path_A_changed"]].itertuples():
        changes.append(
            {
                "classification_type": "Path_A",
                "family_id": row.family_id,
                "cohort": "cross_layer",
                "legacy": f"{row.support_pattern_legacy}:{row.cross_axis_class_legacy}",
                "corrected": (
                    f"{row.support_pattern_corrected}:"
                    f"{row.cross_axis_class_corrected}"
                ),
            }
        )
    changed = pd.DataFrame(
        changes,
        columns=["classification_type", "family_id", "cohort", "legacy", "corrected"],
    )
    changed.to_csv(
        RESULTS / "migration_changed_classifications.tsv", sep="\t", index=False
    )
    return changed


def dependency_audit() -> None:
    patterns = {
        "legacy_program_scores": re.compile(r"Project_P1.*program_scores_long|program_scores_long\\.tsv"),
        "legacy_primary_alias": re.compile(r"fig2_primary_(?:acai|axis)"),
        "legacy_pooled_effect": re.compile(r"fig2d_family_forest(?:_v3)?\\.tsv"),
        "legacy_consistency": re.compile(r"fig2e_replication_lollipop(?:_v3)?\\.tsv"),
        "legacy_acairel_export": re.compile(r"fig2b_(?:program|family)_donor_acairel\\.tsv"),
    }
    rows = []
    for path in sorted(SCRIPT_ROOT.rglob("*")):
        if path.suffix not in {".py", ".R"}:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        hits = [name for name, pattern in patterns.items() if pattern.search(text)]
        if not hits:
            continue
        relative = str(path.relative_to(ROOT))
        is_reanalysis = "ncp_v3_singscore_reanalysis" in relative
        rows.append(
            {
                "script": relative,
                "legacy_dependency_types": ";".join(hits),
                "classification": (
                    "reanalysis_audit_or_migration_only"
                    if is_reanalysis
                    else "legacy_or_primary_consumer_quarantined_from_corrected_run"
                ),
                "included_in_corrected_execution": False if not is_reanalysis else True,
            }
        )
    pd.DataFrame(rows).to_csv(
        RESULTS / "supplementary_dependency_audit.tsv", sep="\t", index=False
    )

    corrected_runtime = [
        "score_programs.R",
        "compute_acairel.py",
        "regenerate_fig2.py",
        "regenerate_fig5.py",
    ]
    forbidden = [
        "Project_P1/program_scores_long.tsv",
        "fig2_primary_acai",
        "fig2d_family_forest_v3.tsv",
        "fig2e_replication_lollipop_v3.tsv",
    ]
    lineage = []
    for filename in corrected_runtime:
        path = SCRIPT_ROOT / "ncp_v3_singscore_reanalysis" / filename
        text = path.read_text(encoding="utf-8")
        hits = [token for token in forbidden if token in text]
        lineage.append(
            {
                "script": str(path.relative_to(ROOT)),
                "forbidden_legacy_bulk_hits": ";".join(hits),
                "passed": len(hits) == 0,
            }
        )
    lineage_df = pd.DataFrame(lineage)
    lineage_df.to_csv(RESULTS / "corrected_lineage_audit.tsv", sep="\t", index=False)
    if not lineage_df["passed"].all():
        raise AssertionError("Corrected runtime mixes legacy bulk inputs")


def main() -> None:
    write_m_comparison()
    write_donor_comparisons()
    changed = write_classification_comparisons()
    dependency_audit()
    summary = (
        changed.groupby("classification_type").size().rename("n_changed").reset_index()
        if not changed.empty
        else pd.DataFrame(columns=["classification_type", "n_changed"])
    )
    summary.to_csv(
        RESULTS / "migration_changed_classification_summary.tsv", sep="\t", index=False
    )
    print(f"Migration audit complete; changed classification rows={len(changed)}")


if __name__ == "__main__":
    main()
