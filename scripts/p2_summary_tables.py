#!/usr/bin/env python3
"""
Build three summary tables from existing `results/` and `scripts/` artifacts:

  P-2 Power table              results/power_table.tsv
  P-3 Multi-test correction    results/mtc_summary_pipeline.tsv
     (Supplementary-style 8-column table: run scripts/p3_mtc_summary.py → results/mtc_summary.tsv)
  P-4 Reproducibility table    results/reproducibility_table.tsv

No re-analysis; tabling only. `spec_source` points at implementing scripts.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"


# ---------------------------------------------------------------------------
# B2 Power table  (collapses h8 simulation + cohort sample sizes)
# ---------------------------------------------------------------------------
def build_power_table() -> None:
    sim = pd.read_csv(RES / "method_p5_simulation.tsv", sep="\t")
    rows = []
    for _, r in sim.iterrows():
        rows.append({
            "true_effect": r["effect_size"],
            "recovered_pooled_g_mean": r["mean_pooled_g"],
            "recovered_pooled_g_sd":   r["sd_pooled_g"],
            "detection_rate_g_gt_0p3": r["detect_rate"],
            "n_simulations":           r["reps"],
            "cohort_donor_structure":  "5+5, 7+6, 13+15  (matched to "
                                        "GSE128188, GSE79768, GSE115574)",
        })
    df = pd.DataFrame(rows)
    df.to_csv(RES / "power_table.tsv", sep="\t", index=False)
    print(f"wrote results/power_table.tsv  ({len(df)} rows)")
    print(df.to_string(index=False))


# ---------------------------------------------------------------------------
# B3 Multi-test correction summary  (one table; per figure where it applies)
# ---------------------------------------------------------------------------
def build_mtc_table() -> None:
    rows = [
        # (figure, statistic, count, family, correction, threshold, source spec)
        ("Fig 2 rewiring class",
         "permutation gate on |ΔaCAI_rel|",
         "16 axis families × 3 cohorts (gate per family × cohort)",
         "per-cohort, per-family",
         "permutation; q_perm BH within cohort",
         "q_perm < 0.10 for class assignment",
         "scripts/p4_acai_classify.py"),
        ("Fig 3 LA-arm Δsingscore",
         "permutation Δsingscore (AF − Ctrl)",
         "16 families × 4 main compartments = 64 cells",
         "per family (across the 4 compartments)",
         "BH-FDR within family",
         "q_perm < 0.10 to assign carrier-fidelity class (Below-detection otherwise)",
         "scripts/r4_delta_classify.py"),
        ("Fig 4 LAA duration ρ",
         "Spearman ρ + Monte-Carlo permutation",
         "16 axis families × 4 statistical conditions "
         "(LAA-primary, LAA-sens-A, LAA-ordinal, LAPV-sens)",
         "region-specific (LAA family separate from LAPV family); "
         "16 families per region",
         "BH-FDR per region per condition",
         "q_perm < 0.10 for `sig_marker = .`; q_perm < 0.05 for `*`",
         "scripts/d4_duration_correlate.py"),
        ("Fig 5A triangulation class",
         "rule-based composite (no p-value)",
         "16 axis families",
         "n/a (deterministic rule)",
         "n/a",
         "auditable in script header",
         "scripts/fig5a_triangulation.py"),
        ("Fig 5B multiome support gate",
         "Δsingscore (AF − SR) on snRNA pseudobulk",
         "6 programs × 4 main compartments = 24 cells",
         "n/a (gate is a pre-registered direction call)",
         "no FDR (effect-size + direction)",
         "moderate / strong support tier (locked support-gate)",
         "scripts/figS5_multiome_delta_heatmap.py"),
        ("Fig 6A axis pooled meta",
         "Hedges' g (AF − SR) per cohort + REML / Hartung–Knapp pool",
         "2 axes (STRUCTURAL, METABOLIC); k = 3 cohorts",
         "no multiple-testing correction "
         "(framework reports g, CI, I², τ², p_HK individually)",
         "Hartung–Knapp variance adjustment",
         "consistency-check framing; not a high-powered pooled discovery",
         "scripts/h2_meta_axis.R"),
        ("Fig 6B per-program meta",
         "Hedges' g per program + REML / HK pool",
         "30 NCP pathway programs",
         "no MTC (descriptive consistency summary)",
         "Hartung–Knapp",
         "I² + p_HK reported per program",
         "scripts/h3_meta_program.R"),
        ("Fig 6A/6B cis-eQTL × AF GWAS coloc",
         "coloc.abf default priors (p1=p2=1e-4, p12=1e-5)",
         "12 pre-specified genes (GTEx Atrial Appendage primary; "
         "GTEx Whole Blood sensitivity)",
         "n/a (small pre-specified anchor set; support flags, not BH)",
         "local AF GWAS / atrial eGene / PP4 ≥ 0.70 support flags; "
         "PP3-dominant configuration readout",
         "PP4 ≥ 0.70 support; PP3-dominant when PP3 > PP4; "
         "lead separation descriptive (0.46–1.92 Mb)",
         "scripts/fig6a_anchor_coloc_calibration.py + "
         "scripts/fig6b_posterior_lead_distance.py"),
        ("Fig S8a size-matched perms",
         "size-matched random gene-set pooled g",
         "6 frozen support programs × 1000 perms each",
         "per program",
         "two-sided empirical p_emp",
         "p_emp reported; no multi-program adjustment "
         "(specificity check, not a discovery screen)",
         "scripts/figS8_falsification.py (size-matched)"),
        ("Fig S8b axis-membership perm",
         "axis-membership permutation pooled |g|",
         "6 frozen IDs randomly partitioned 3+3 × 1000 perms",
         "per axis",
         "two-sided empirical p_emp",
         "p_emp reported per axis",
         "scripts/figS8_falsification.py (axis-membership)"),
        ("Fig S8c negative controls",
         "pooled g + I² + LAA ρ for 6 pre-registered control pathways",
         "6 controls × 3 statistics",
         "n/a (descriptive comparison)",
         "no MTC; comparison is qualitative",
         "values reported alongside real frozen programs",
         "scripts/figS8_falsification.py (negative controls)"),
        ("Fig S8d LOPO axis robustness",
         "pooled g for axis with each component program removed",
         "2 axes × 3 LOPO removals",
         "per axis",
         "no MTC (direction-preservation tag)",
         "label: robust / partially-sensitive / single-program-dependent",
         "scripts/h9_lopo_axis.R + scripts/figS8_falsification.py"),
        ("Fig S10 methodology audit",
         "NCP self-recovery; donor-level direction; cross-method ρ; "
         "cross-method axis pooled g",
         "30 NCP pathway programs × 4 tests",
         "n/a (internal consistency)",
         "no MTC",
         "self-recovery + 8/8 direction stability + "
         "cross-method ρ ≥ 0.5",
         "scripts/h7_methodology_audit.R"),
        ("Fig S11 P5 method validity",
         "simulation; bootstrap; classical limma comparison",
         "various",
         "n/a",
         "no MTC; cross-validation framing",
         "see scripts/h8_methodology_p5.R",
         "scripts/h8_methodology_p5.R"),
    ]
    cols = ["figure", "statistic", "count", "family", "correction",
            "threshold_or_call", "spec_source"]
    df = pd.DataFrame(rows, columns=cols)
    out = RES / "mtc_summary_pipeline.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out.relative_to(ROOT)}  ({len(df)} rows)")


# ---------------------------------------------------------------------------
# B4 Reproducibility table  (figure → script → input → output)
# ---------------------------------------------------------------------------
def build_reproducibility_table() -> None:
    rows = [
        ("Fig 1", "scripts/assemble_figure1.py",
         "figures/fig1a_donor_paired_estimand.pdf + "
         "figures/ncp_v3_singscore_reanalysis/"
         "fig1b_healthy_reference_direction_magnitude.pdf + "
         "figures/fig1c_selection_provenance.pdf + "
         "figures/fig1d_program_family_topology.pdf",
         "figures/figure1_composite.pdf"),
        ("Fig 1A", "scripts/fig1a_donor_paired_estimand.py",
         "results/fig1a_donor_paired_estimand_source.tsv",
         "figures/fig1a_donor_paired_estimand.pdf"),
        ("Fig 1B",
         "scripts/ncp_v3_singscore_reanalysis/"
         "fig1b_healthy_reference_direction_magnitude.py",
         "results/ncp_v3_singscore_reanalysis/fig1b_input_manifest.tsv",
         "figures/ncp_v3_singscore_reanalysis/"
         "fig1b_healthy_reference_direction_magnitude.pdf"),
        ("Fig 1C", "scripts/fig1c_healthy_selection_frozen_deployment.py",
         "results/fig1c_pathway_selection_long.tsv + "
         "results/innovation_ncp_consensus_long.tsv + "
         "results/innovation_ncp_v3_panel.tsv + "
         "results/innovation_ncp_v3_summary.tsv + "
         "results/ncp_v3_selection_criteria.tsv + "
         "results/innovation_ncp_v3_axis_families.tsv (validation-only)",
         "figures/fig1c_selection_provenance.pdf"),
        ("Fig 1D", "scripts/fig1d_program_family_topology.py",
         "results/innovation_ncp_v3_axis_families.tsv + "
         "results/ncp_v3_family_labels.tsv + "
         "results/ncp_v3_canonical_display_order.tsv + "
         "results/ncp_v3_singscore_reanalysis/fig1b_canonical_display_order.tsv",
         "figures/fig1d_program_family_topology.pdf"),
        ("Fig 2A",
         "scripts/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_bridge.py",
         "results/ncp_v3_singscore_reanalysis/fig2a_input_manifest.tsv",
         "figures/ncp_v3_singscore_reanalysis/"
         "fig2a_frozen_ncp_bulk_deployment.pdf"),
        ("Fig 2B",
         "scripts/ncp_v3_singscore_reanalysis/fig2b_family_donor_heatmap.py",
         "results/ncp_v3_singscore_reanalysis/fig2a_input_manifest.tsv",
         "figures/fig2b_donor_acai_heatmap.pdf + "
         "figures/fig2b_donor_acai_heatmap.svg + "
         "results/ncp_v3_singscore_reanalysis/fig2b_family_donor_source.tsv"),
        ("Fig 2C",
         "scripts/ncp_v3_singscore_reanalysis/fig2c_family_cohort_state_map.py",
         "results/ncp_v3_singscore_reanalysis/fig2a_input_manifest.tsv",
         "figures/ncp_v3_singscore_reanalysis/fig2c_family_cohort_state_map.pdf + "
         "figures/ncp_v3_singscore_reanalysis/fig2c_family_cohort_state_map.svg + "
         "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv + "
         "results/ncp_v3_singscore_reanalysis/figure2_family_row_geometry.tsv"),
        ("Fig 2D", "scripts/fig2d_representative_dots.py",
         "Project_P1/program_scores_long.tsv",
         "figures/fig2d_representative_dots.pdf"),
        ("Fig 2E", "scripts/fig2e_replication_lollipop.py",
         "results/af_chamber_rewiring_consensus.tsv",
         "figures/fig2e_replication_lollipop.pdf"),
        ("Fig 3",
         "scripts/ncp_v3_singscore_reanalysis/"
         "assemble_figure3_la_single_nucleus_localization.py",
         "results/ncp_v3_singscore_reanalysis/"
         "fig3b_family_compartment_localization_source.tsv + "
         "results/ncp_v3_singscore_reanalysis/fig3c_fidelity_relay_donor_source.tsv + "
         "results/ncp_v3_singscore_reanalysis/gse255612_family_donor_compartment_scores.tsv + "
         "results/ncp_v3_singscore_reanalysis/"
         "gse255612_corrected_family_compartment_permutation_classification.tsv + "
         "results/ncp_v3_singscore_reanalysis/figure3cd_case_study_axis_contract.tsv + "
         "results/ncp_v3_singscore_reanalysis/figure3cd_case_study_composite_geometry.tsv",
         "figures/ncp_v3_singscore_reanalysis/"
         "figure3_la_single_nucleus_localization.pdf + "
         "figures/ncp_v3_singscore_reanalysis/"
         "figure3_la_single_nucleus_localization.svg + "
         "results/ncp_v3_singscore_reanalysis/figure3_composite_manifest.tsv + "
         "results/ncp_v3_singscore_reanalysis/figure3_composite_qc.md"),
        ("Fig 4", "scripts/ncp_v3_singscore_reanalysis/regenerate_figure4.py",
         "results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv + "
         "gse41177_family_sample_scores.tsv + gse41177_duration_metadata.tsv",
         "figures/ncp_v3_singscore_reanalysis/figure4_cross_sectional_duration.pdf + "
         "figures/ncp_v3_singscore_reanalysis/figure4_cross_sectional_duration.svg + "
         "results/ncp_v3_singscore_reanalysis/figure4_composite_manifest.tsv + "
         "results/ncp_v3_singscore_reanalysis/figure4_composite_qc.md"),
        ("Fig 5A", "scripts/fig5a_triangulation.py",
         "results/figure3_program_panel.tsv + figure3_carrier_shift.tsv + "
         "Project_P4/GSE41177/results/duration_trend_laa_primary.tsv",
         "figures/fig5a_triangulation.pdf + "
         "results/figure5_triangulation_table.tsv"),
        ("Fig 5B", "scripts/fig5b_design_overview.py",
         "(schematic, no data)",
         "figures/fig5b_design_overview.pdf"),
        ("Fig 5C", "scripts/fig5c_localization_inset.py",
         "Project_P5/GSE261363/roi_scores.tsv",
         "figures/fig5c_localization_inset.pdf"),
        ("Fig 5D", "scripts/fig5d_dual_axis_mechanism.py",
         "(schematic, no data)",
         "figures/fig5d_dual_axis_mechanism.pdf"),
        ("Fig 6A", "scripts/fig6a_anchor_coloc_calibration.py",
         "results/figure6_mr_anchors_v3_pipeline_wide.tsv + "
         "results/fig6a_v3_locus_zoom_summary.tsv (frozen coloc artifacts; "
         "KCNN3 lead dist via scripts/compute_anchor_lead_distance.py)",
         "figures/fig6a_anchor_coloc_calibration.pdf + "
         "results/fig6a_anchor_coloc_calibration_source.tsv"),
        ("Fig 6B", "scripts/fig6b_posterior_lead_distance.py",
         "results/figure6_mr_anchors_v3_pipeline_wide.tsv + "
         "results/fig6a_v3_locus_zoom_summary.tsv (frozen coloc artifacts)",
         "figures/fig6b_posterior_lead_distance.pdf + "
         "results/fig6b_posterior_lead_distance_source.tsv"),
        ("Fig 6C", "scripts/fetch_external_af_gene_sets.py + "
         "scripts/fig6c_family_af_overlap.py",
         "results/ncp_v3_family_leading_edge_genes.tsv + "
         "results/ncp_v3_gene_universe.tsv + "
         "results/external_af_gene_sets.tsv",
         "figures/fig6c_family_af_overlap.pdf + "
         "results/fig6c_family_af_overlap_source.tsv"),
        ("Fig 6D", "scripts/fig6d_non_nested_traceability_profiles.py",
         "results/fig6a_anchor_coloc_calibration_source.tsv + "
         "results/fig6b_posterior_lead_distance_source.tsv + "
         "results/external_af_gene_sets.tsv (frozen; non-nested flags)",
         "figures/fig6d_traceability_matrix.pdf + "
         "results/fig6d_non_nested_traceability_source.tsv"),
        ("Fig 6 composite", "scripts/assemble_figure6_genetic_traceability.py",
         "figures/fig6a_anchor_coloc_calibration.svg + "
         "figures/fig6b_posterior_lead_distance.svg + "
         "figures/fig6c_family_af_overlap.svg + "
         "figures/fig6d_traceability_matrix.svg (exact-size nest)",
         "figures/figure6_genetic_traceability.pdf + "
         "results/figure6_composite_manifest.tsv"),
        ("Fig 7", "scripts/n2_chamber_reversal_genes.py + "
         "scripts/n2b_reversal_celltype.py + "
         "scripts/fig7_reversal_carrier_composite.py",
         "results/n2_chamber_reversal_long.tsv + n2_chamber_reversal_genes.tsv + "
         "n2b_reversal_per_celltype_long.tsv + n2b_reversal_per_celltype_summary.tsv "
         "(frozen reversal / carrier artifacts)",
         "figures/fig7_reversal_carrier_composite.pdf"),
        ("Graphical abstract", "scripts/graphical_abstract.py",
         "(schematic)",
         "figures/graphical_abstract.pdf"),
        ("Fig S5",  "scripts/figS5_multiome_delta_heatmap.py",
         "Project_P5/GSE238242/results/figure5_program_delta_main.tsv",
         "figures/figS5_multiome_delta_heatmap.pdf"),
        ("Fig S6",  "scripts/f5b_program_singscore_mural.R + "
         "scripts/figS6_mural_sensitivity.py",
         "Project_P5/GSE238242/pseudobulk/gse238242_pseudobulk_counts_supp_mural.rds",
         "figures/figS6_mural_sensitivity.pdf"),
        ("Fig S7",  "scripts/figS7_loco_sensitivity.py",
         "results/figure6_meta_axis_loco.tsv",
         "figures/figS7_loco_sensitivity.pdf"),
        ("Fig S8",  "scripts/h5_falsification.R + h9_lopo_axis.R + "
         "scripts/figS8_falsification.py",
         "Project_P1/program_scores_long.tsv + cohort expression matrices + "
         "MSigDB GMT files",
         "figures/figS8_falsification.pdf + "
         "results/figure_s8a_size_matched.tsv + figure_s8b_axis_perm.tsv + "
         "figure_s8c_negative_controls.tsv + figure_s8d_lopo_axis.tsv"),
        ("Fig S9",  "scripts/h6_reduced_axis.R + scripts/figS9_reduced_axis.py",
         "Project_P1/<cohort>/expr_symbol.tsv + sample_metadata.tsv + "
         "Project_P4/GSE41177/expr_symbol_LAA.tsv + figure5_program_sets.tsv",
         "figures/figS9_reduced_axis.pdf + results/figure_s9_*.tsv"),
        ("Fig S10", "scripts/h7_methodology_audit.R + "
         "scripts/figS10_methodology_audit.py",
         "Project_P1/<cohort>/expr_symbol.tsv + bulk/logcpm_voom.tsv",
         "figures/figS10_methodology_audit.pdf + "
         "results/methodology_audit_program.tsv + methodology_audit_axis_g.tsv"),
        ("Fig S11", "scripts/h8_methodology_p5.R + "
         "scripts/figS11_method_p5.py",
         "Project_P1/program_scores_long.tsv + cohort expression matrices",
         "figures/figS11_method_p5.pdf + "
         "results/method_p5_simulation.tsv + method_p5_bootstrap_axis.tsv + "
         "method_p5_classical_vs_framework.tsv"),
    ]
    cols = ["figure", "script", "input", "output"]
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(RES / "reproducibility_table.tsv", sep="\t", index=False)
    print(f"\nwrote results/reproducibility_table.tsv  ({len(df)} rows)")


def main() -> None:
    build_power_table()
    build_mtc_table()
    build_reproducibility_table()


if __name__ == "__main__":
    main()
