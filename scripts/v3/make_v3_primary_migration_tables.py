#!/usr/bin/env python3
"""Build migration checklists for the v3-primary switch.

Outputs
-------
results/v3_primary_replacement_checklist.tsv
results/v3_affected_results_comparison.tsv
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"


def build_replacement_checklist() -> pd.DataFrame:
    rows = [
        {
            "layer": "fig2_input_alias",
            "old_or_virtual_target": "results/fig2_primary_acai_long.tsv",
            "new_source": "results/innovation_ncp_v3_acai_long.tsv",
            "action": "copied_from_v3_source",
        },
        {
            "layer": "fig2_input_alias",
            "old_or_virtual_target": "results/fig2_primary_acai_summary.tsv",
            "new_source": "results/innovation_ncp_v3_acai_summary.tsv",
            "action": "copied_from_v3_source",
        },
        {
            "layer": "fig2_input_alias",
            "old_or_virtual_target": "results/fig2_primary_axis_families.tsv",
            "new_source": "results/innovation_ncp_v3_axis_families.tsv",
            "action": "copied_from_v3_source",
        },
        {
            "layer": "fig3_input_alias",
            "old_or_virtual_target": "results/fig3_primary_program_celltype_delta.tsv",
            "new_source": "results/fig3_v2_program_celltype_delta.tsv",
            "action": "copied_before_v2_cleanup",
        },
        {
            "layer": "fig3_input_alias",
            "old_or_virtual_target": "results/fig3_primary_family_celltype_delta.tsv",
            "new_source": "results/fig3_v2_family_celltype_delta.tsv",
            "action": "copied_before_v2_cleanup",
        },
        {
            "layer": "fig3_input_alias",
            "old_or_virtual_target": "results/fig3_primary_program_carrier.tsv",
            "new_source": "results/fig3_v2_program_carrier.tsv",
            "action": "copied_before_v2_cleanup",
        },
        {
            "layer": "fig3_input_alias",
            "old_or_virtual_target": "results/fig3_primary_source_data.tsv",
            "new_source": "results/fig3_v2_source_data.tsv",
            "action": "copied_before_v2_cleanup",
        },
        {
            "layer": "fig3_input_alias",
            "old_or_virtual_target": "results/fig3_primary_program_le.tsv",
            "new_source": "results/fig3_ncpv2_program_le.tsv",
            "action": "copied_before_v2_cleanup",
        },
        {
            "layer": "fig4_input_alias",
            "old_or_virtual_target": "Project_P4/GSE41177/program_scores_primary.tsv",
            "new_source": "Project_P4/GSE41177/program_scores_v2.tsv",
            "action": "copied_before_v2_cleanup",
        },
        {
            "layer": "fig4_input_alias",
            "old_or_virtual_target": "Project_P4/GSE41177/results_primary",
            "new_source": "Project_P4/GSE41177/results_v2",
            "action": "copied_before_v2_cleanup",
        },
        {
            "layer": "fig6_input_alias",
            "old_or_virtual_target": "results/n2b_primary_reversal_per_celltype_long.tsv",
            "new_source": "results/n2b_v2_reversal_per_celltype_long.tsv",
            "action": "copied_before_v2_cleanup",
        },
        {
            "layer": "fig6_input_alias",
            "old_or_virtual_target": "results/n2b_primary_reversal_per_celltype_summary.tsv",
            "new_source": "results/n2b_v2_reversal_per_celltype_summary.tsv",
            "action": "copied_before_v2_cleanup",
        },
        {
            "layer": "fig3_input_alias",
            "old_or_virtual_target": "Project_P3/GSE255612/program_scores_primary.tsv",
            "new_source": "Project_P3/GSE255612/program_scores_v2.tsv",
            "action": "copied_before_v2_cleanup",
        },
        {
            "layer": "fig2_main_results",
            "old_or_virtual_target": "results/fig2c_chamber_decomposition.tsv",
            "new_source": "results/fig2c_chamber_decomposition_v3.tsv",
            "action": "promoted_to_primary_name",
        },
        {
            "layer": "fig2_main_results",
            "old_or_virtual_target": "results/fig2d_family_forest.tsv",
            "new_source": "results/fig2d_family_forest_v3.tsv",
            "action": "promoted_to_primary_name",
        },
        {
            "layer": "fig2_main_results",
            "old_or_virtual_target": "results/fig2e_replication_lollipop.tsv",
            "new_source": "results/fig2e_replication_lollipop_v3.tsv",
            "action": "promoted_to_primary_name",
        },
        {
            "layer": "fig2_main_figure",
            "old_or_virtual_target": "figures/fig2b_donor_acai_heatmap.pdf",
            "new_source": "figures/fig2b_donor_acai_heatmap_v3.pdf",
            "action": "promoted_to_primary_name",
        },
        {
            "layer": "fig2_main_figure",
            "old_or_virtual_target": "figures/fig2c_chamber_decomposition.pdf",
            "new_source": "figures/fig2c_chamber_decomposition_v3.pdf",
            "action": "promoted_to_primary_name",
        },
        {
            "layer": "fig2_main_figure",
            "old_or_virtual_target": "figures/fig2d_family_forest.pdf",
            "new_source": "figures/fig2d_family_forest_v3.pdf",
            "action": "promoted_to_primary_name",
        },
        {
            "layer": "fig2_main_figure",
            "old_or_virtual_target": "figures/fig2e_replication_lollipop.pdf",
            "new_source": "figures/fig2e_replication_lollipop_v3.pdf",
            "action": "promoted_to_primary_name",
        },
        {
            "layer": "downstream_script_default",
            "old_or_virtual_target": "scripts/fig3b_family_celltype_heatmap.py",
            "new_source": "results/fig2_primary_acai_long.tsv + results/fig2_primary_axis_families.tsv",
            "action": "default_input_switched_to_primary",
        },
        {
            "layer": "downstream_script_default",
            "old_or_virtual_target": "scripts/fig3d_program_celltype_heatmap.py",
            "new_source": "results/fig2_primary_acai_long.tsv + results/fig2_primary_axis_families.tsv",
            "action": "default_input_switched_to_primary",
        },
        {
            "layer": "downstream_script_default",
            "old_or_virtual_target": "scripts/fig4b_family_duration_heatmap.py",
            "new_source": "results/fig2_primary_acai_long.tsv + results/fig2_primary_axis_families.tsv",
            "action": "default_input_switched_to_primary",
        },
        {
            "layer": "downstream_script_default",
            "old_or_virtual_target": "scripts/fig5b_triangulation_heatmap.py",
            "new_source": "results/fig2_primary_acai_long.tsv + results/fig2_primary_axis_families.tsv",
            "action": "default_input_switched_to_primary",
        },
    ]
    out = pd.DataFrame(rows)
    out["exists_target"] = out["old_or_virtual_target"].apply(
        lambda p: (ROOT / p).exists()
    )
    out["exists_source"] = out["new_source"].apply(
        lambda p: all((ROOT / piece.strip()).exists() for piece in p.split("+"))
    )
    out["status"] = out.apply(
        lambda r: "done" if bool(r["exists_target"]) and bool(r["exists_source"]) else "missing",
        axis=1,
    )
    return out


def build_affected_results_comparison() -> pd.DataFrame:
    metrics = pd.read_csv(RESULTS / "innovation_v2_vs_v3_metrics.tsv", sep="\t").iloc[0]
    decision = pd.read_csv(RESULTS / "innovation_v2_vs_v3_decision_table.tsv", sep="\t")

    rows = [
        {
            "domain": "panel",
            "metric": "panel_size",
            "v2_value": float(metrics["panel_size_v2"]),
            "v3_value": float(metrics["panel_size_v3"]),
            "delta_v3_minus_v2": float(metrics["panel_size_v3"] - metrics["panel_size_v2"]),
            "rule_or_note": "matched panel size",
        },
        {
            "domain": "panel",
            "metric": "overlap_fraction_vs_v2",
            "v2_value": 1.0,
            "v3_value": float(metrics["overlap_fraction_vs_v2"]),
            "delta_v3_minus_v2": float(metrics["overlap_fraction_vs_v2"] - 1.0),
            "rule_or_note": "continuity >= 0.80",
        },
        {
            "domain": "direction_balance",
            "metric": "ra_count",
            "v2_value": float(metrics["ra_count_v2"]),
            "v3_value": float(metrics["ra_count_v3"]),
            "delta_v3_minus_v2": float(metrics["ra_count_v3"] - metrics["ra_count_v2"]),
            "rule_or_note": "RA coverage >= 8",
        },
        {
            "domain": "redundancy",
            "metric": "mean_pairwise_jaccard",
            "v2_value": float(metrics["mean_pairwise_jaccard_v2"]),
            "v3_value": float(metrics["mean_pairwise_jaccard_v3"]),
            "delta_v3_minus_v2": float(
                metrics["mean_pairwise_jaccard_v3"] - metrics["mean_pairwise_jaccard_v2"]
            ),
            "rule_or_note": "increase <= 0.005",
        },
    ]

    out = pd.DataFrame(rows)
    rule_map = {
        "Panel continuity": "overlap_fraction_vs_v2",
        "Direction balance (RA coverage)": "ra_count",
        "Redundancy control": "mean_pairwise_jaccard",
        "No cohort sign flip (all-program delta)": "n_sign_flip",
        "Transferability (delta stability)": "mean_abs_delta_shift",
        "Transferability (effect-size stability)": "mean_abs_g_shift",
    }

    decisions = decision.copy()
    decisions["metric"] = decisions["criterion"].map(rule_map)
    decisions = decisions[["criterion", "metric", "rule", "observed", "pass", "overall_decision"]]
    decisions = decisions.rename(
        columns={
            "criterion": "domain",
            "rule": "rule_or_note",
            "observed": "v3_value",
            "pass": "rule_pass",
        }
    )
    decisions["v2_value"] = pd.NA
    decisions["delta_v3_minus_v2"] = pd.NA
    decisions["metric"] = decisions["metric"].fillna("decision_gate")
    decisions = decisions[
        ["domain", "metric", "v2_value", "v3_value", "delta_v3_minus_v2", "rule_or_note", "rule_pass", "overall_decision"]
    ]

    out["rule_pass"] = pd.NA
    out["overall_decision"] = "upgrade-ready"
    out = out[
        ["domain", "metric", "v2_value", "v3_value", "delta_v3_minus_v2", "rule_or_note", "rule_pass", "overall_decision"]
    ]
    decisions = decisions.assign(
        v2_value="",
        delta_v3_minus_v2="",
        rule_pass=decisions["rule_pass"].astype("string"),
    )
    out = out.assign(
        rule_pass=out["rule_pass"].astype("string"),
    )
    return pd.concat([out, decisions], ignore_index=True)


def build_v2_cleanup_plan() -> pd.DataFrame:
    rows = [
        {
            "item_path": "results/fig2_ncp_v2_acai_long.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by results/fig2_primary_acai_long.tsv",
        },
        {
            "item_path": "results/fig2_ncp_v2_acai_summary.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by results/fig2_primary_acai_summary.tsv",
        },
        {
            "item_path": "results/fig3_v2_program_celltype_delta.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by results/fig3_primary_program_celltype_delta.tsv",
        },
        {
            "item_path": "results/fig3_v2_family_celltype_delta.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by results/fig3_primary_family_celltype_delta.tsv",
        },
        {
            "item_path": "results/fig3_v2_program_carrier.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by results/fig3_primary_program_carrier.tsv",
        },
        {
            "item_path": "results/fig3_v2_source_data.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by results/fig3_primary_source_data.tsv",
        },
        {
            "item_path": "results/fig3_ncpv2_program_le.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by results/fig3_primary_program_le.tsv",
        },
        {
            "item_path": "results/n2b_v2_reversal_per_celltype_long.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by results/n2b_primary_reversal_per_celltype_long.tsv",
        },
        {
            "item_path": "results/n2b_v2_reversal_per_celltype_summary.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by results/n2b_primary_reversal_per_celltype_summary.tsv",
        },
        {
            "item_path": "Project_P3/GSE255612/program_scores_v2.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by Project_P3/GSE255612/program_scores_primary.tsv",
        },
        {
            "item_path": "Project_P4/GSE41177/program_scores_v2.tsv",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by Project_P4/GSE41177/program_scores_primary.tsv",
        },
        {
            "item_path": "Project_P4/GSE41177/results_v2",
            "kind": "legacy_v2_runtime",
            "action": "delete",
            "target_section": "none",
            "replacement_or_note": "replaced by Project_P4/GSE41177/results_primary",
        },
        {
            "item_path": "figures/innovation_v2_vs_v3_benchmark.pdf",
            "kind": "method_comparison",
            "action": "keep_and_relabel",
            "target_section": "Figure S13",
            "replacement_or_note": "copied to figures/figS13_v2_vs_v3_benchmark.pdf",
        },
        {
            "item_path": "figures/innovation_ncp_bootstrap_frequency.pdf",
            "kind": "method_validation",
            "action": "keep_and_relabel",
            "target_section": "Figure S14",
            "replacement_or_note": "copied to figures/figS14_v3_bootstrap_stability.pdf",
        },
        {
            "item_path": "results/innovation_v2_vs_v3_metrics.tsv",
            "kind": "method_comparison_source",
            "action": "keep",
            "target_section": "Supplementary methods benchmark",
            "replacement_or_note": "source for Fig S13",
        },
        {
            "item_path": "results/innovation_v2_vs_v3_axis_effects.tsv",
            "kind": "method_comparison_source",
            "action": "keep",
            "target_section": "Supplementary methods benchmark",
            "replacement_or_note": "source for Fig S13",
        },
        {
            "item_path": "results/innovation_v2_vs_v3_decision_table.tsv",
            "kind": "method_comparison_source",
            "action": "keep",
            "target_section": "Supplementary methods benchmark",
            "replacement_or_note": "source for Fig S13",
        },
        {
            "item_path": "results/innovation_ncp_v2_acai_long.tsv",
            "kind": "method_comparison_source",
            "action": "keep",
            "target_section": "Supplementary methods benchmark",
            "replacement_or_note": "needed for reproducible v2-v3 comparison",
        },
        {
            "item_path": "results/innovation_ncp_v2_acai_summary.tsv",
            "kind": "method_comparison_source",
            "action": "keep",
            "target_section": "Supplementary methods benchmark",
            "replacement_or_note": "needed for reproducible v2-v3 comparison",
        },
        {
            "item_path": "results/fig1c_ncp_panel_v2.tsv",
            "kind": "baseline_panel_definition",
            "action": "keep",
            "target_section": "Supplementary methods benchmark",
            "replacement_or_note": "frozen baseline panel used for all v2 comparisons",
        },
    ]
    out = pd.DataFrame(rows)
    out["exists_now"] = out["item_path"].apply(lambda p: (ROOT / p).exists())
    return out


def main() -> None:
    checklist = build_replacement_checklist()
    affected = build_affected_results_comparison()
    cleanup = build_v2_cleanup_plan()

    out_a = RESULTS / "v3_primary_replacement_checklist.tsv"
    out_b = RESULTS / "v3_affected_results_comparison.tsv"
    out_c = RESULTS / "v2_cleanup_and_figureS_plan.tsv"
    checklist.to_csv(out_a, sep="\t", index=False)
    affected.to_csv(out_b, sep="\t", index=False)
    cleanup.to_csv(out_c, sep="\t", index=False)
    print(f"wrote {out_a.relative_to(ROOT)} ({len(checklist)} rows)")
    print(f"wrote {out_b.relative_to(ROOT)} ({len(affected)} rows)")
    print(f"wrote {out_c.relative_to(ROOT)} ({len(cleanup)} rows)")


if __name__ == "__main__":
    main()
