#!/usr/bin/env python3
"""Assemble a PLOS CB–style submission folder (Word + PDF + key TSV)."""
from __future__ import annotations

import shutil
from pathlib import Path

from audit_publication_ncp_nomenclature import run_qc as run_ncp_nomenclature_qc
from migrate_publication_ncp_nomenclature import migrate_docx

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
DST = ROOT / "submission" / "PLOS_CB_submission_pack"

MAIN_FIGURES: list[tuple[str, Path]] = [
    ("Figure_1", FIG / "figure1_composite.pdf"),
    ("Figure_2", FIG / "figure2_composite.pdf"),
    (
        "Figure_2A",
        FIG
        / "ncp_v3_singscore_reanalysis"
        / "fig2a_frozen_ncp_bulk_deployment.pdf",
    ),
    ("Figure_2B", FIG / "fig2b_donor_acai_heatmap.pdf"),
    (
        "Figure_2C",
        FIG
        / "ncp_v3_singscore_reanalysis"
        / "fig2c_family_cohort_state_map.pdf",
    ),
    (
        "Figure_2D",
        FIG
        / "ncp_v3_singscore_reanalysis"
        / "fig2d_family_effects_consistency.pdf",
    ),
    (
        "Figure_3",
        FIG
        / "ncp_v3_singscore_reanalysis"
        / "figure3_la_single_nucleus_localization.pdf",
    ),
    (
        "Figure_4",
        FIG
        / "ncp_v3_singscore_reanalysis"
        / "figure4_cross_sectional_duration.pdf",
    ),
    ("Figure_5A", FIG / "fig5a_schematic.pdf"),
    ("Figure_5B", FIG / "fig5b_triangulation_heatmap.pdf"),
    ("Figure_5C", FIG / "fig5c_cross_axis.pdf"),
    ("Figure_5D", FIG / "fig5d_geomx_roi.pdf"),
    ("Figure_5E", FIG / "fig5e_mechanism_cartoon.pdf"),
    ("Figure_6", FIG / "figure6_genetic_traceability.pdf"),
    ("Figure_6A", FIG / "fig6a_anchor_coloc_calibration.pdf"),
    ("Figure_6B", FIG / "fig6b_posterior_lead_distance.pdf"),
    ("Figure_6C", FIG / "fig6c_family_af_overlap.pdf"),
    ("Figure_6D", FIG / "fig6d_traceability_matrix.pdf"),
    ("Figure_7", FIG / "fig7_reversal_carrier_composite.pdf"),
]

SUPPORTING_FIGURES: list[tuple[str, Path]] = [
    ("Supporting_Figure_S1", FIG / "figS1_cohort_overview.pdf"),
    ("Supporting_Figure_S2", FIG / "figS2_pairing_attribute_heatmap.pdf"),
    ("Supporting_Figure_S3a_legacy_panel_stability", FIG / "figS3_ncp_panel_stability.pdf"),
    ("Supporting_Figure_S3b_fullmatrix", FIG / "figS3_supp_fullmatrix.pdf"),
    ("Supporting_Figure_S3c_VHS_overlay", FIG / "figS3_VHS_overlay.pdf"),
    ("Supporting_Figure_S4a_PCA", FIG / "figS4_per_cohort_pca.pdf"),
    ("Supporting_Figure_S4b_fullmatrix", FIG / "figS4_supp_fullmatrix.pdf"),
    ("Supporting_Figure_S4c_LAA_LAPV_concordance", FIG / "figS4b_LAA_LAPV_concordance.pdf"),
    ("Supporting_Figure_S5_multiome_gate", FIG / "figS5_multiome_delta_heatmap.pdf"),
    ("Supporting_Figure_S6_mural", FIG / "figS6_mural_sensitivity.pdf"),
    ("Supporting_Figure_S7_LOCO", FIG / "figS7_loco_sensitivity.pdf"),
    ("Supporting_Figure_S8_falsification", FIG / "figS8_falsification.pdf"),
    ("Supporting_Figure_S9_reduced_axis", FIG / "figS9_reduced_axis.pdf"),
    ("Supporting_Figure_S10_audit", FIG / "figS10_methodology_audit.pdf"),
    ("Supporting_Figure_S11_P5", FIG / "figS11_method_p5.pdf"),
    ("Supporting_Figure_S13_prototype_vs_ncp", FIG / "figS13_prototype_vs_ncp_benchmark.pdf"),
    ("Supporting_Figure_S14_bootstrap", FIG / "figS14_ncp_bootstrap_stability.pdf"),
    ("Supporting_Figure_S15_pairing_stress", FIG / "figS15_pairing_stress_test.pdf"),
    (
        "Supporting_Figure_S16_program_donor_acairel",
        FIG / "figS16_program_donor_acairel_heatmap.pdf",
    ),
    (
        "Supporting_Figure_chamber_contribution_to_af_ncp_family_rewiring",
        FIG
        / "ncp_v3_singscore_reanalysis"
        / "supp_chamber_contribution_to_af_ncp_family_rewiring.pdf",
    ),
    ("Supporting_Figure_extra_n2c_sensitivity", FIG / "figS_n2c_panel_sensitivity.pdf"),
    ("Supporting_Figure_extra_bootstrap_frequency", FIG / "innovation_ncp_bootstrap_frequency.pdf"),
    ("Supporting_Figure_fig3_gene_decomp_F15", FIG / "figS_gene_decomp_F15_NTRK2_relay.pdf"),
    ("Supporting_Figure_fig3_gene_decomp_F16", FIG / "figS_gene_decomp_F16_Keratan_extension.pdf"),
    ("Supporting_Figure_fig3_F4_BCAA", FIG / "figS_F4_BCAA_program_exception.pdf"),
    ("Supporting_Figure_fig3_p010_trends", FIG / "figS_fig3_p010_trends.pdf"),
    ("Supporting_Figure_fig3_boundary_perm_hist", FIG / "figS_fig3_boundary_perm_histograms.pdf"),
]

# TSV files explicitly named in manuscript (dxc_benchmark) — reproducibility bundle
RESULTS_TSV = [
    "results/cohort_demographics_long.tsv",
    "results/cohort_demographics_table.tsv",
    "results/fig1c_ncp_panel_v2.tsv",
    "results/fig1c_pathway_selection_long.tsv",
    "results/fig1f_v3_HCA_family_carriers.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_summary.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_qc.md",
    "results/ncp_v3_singscore_reanalysis/figure2_family_row_geometry.tsv",
    "results/ncp_v3_singscore_reanalysis/figure2_family_row_geometry_composite.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_summary.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_qc.md",
    "results/ncp_v3_singscore_reanalysis/legacy_vs_corrected_consistency_changes.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2d_family_effects.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2e_cross_cohort_consistency.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2_replication_classes.tsv",
    "results/ncp_v3_singscore_reanalysis/supp_chamber_contribution_to_af_ncp_family_rewiring_source.tsv",
    "results/ncp_v3_singscore_reanalysis/supp_chamber_contribution_to_af_ncp_family_rewiring_qc.md",
    "results/fig3_v3_carrier_fidelity_class.tsv",
    "results/fig3_v3_singscore_per_donor.tsv",
    "results/fig3_v3_family_celltype_delta.tsv",
    "results/fig3_boundary_permutation_sensitivity.tsv",
    "results/fig3_publication_authority_lock.md",
    "results/ncp_v3_singscore_reanalysis/fig3b_family_compartment_localization_source.tsv",
    "results/ncp_v3_singscore_reanalysis/fig3b_family_compartment_localization_summary.tsv",
    "results/ncp_v3_singscore_reanalysis/fig3b_family_compartment_localization_scale_qc.tsv",
    "results/ncp_v3_singscore_reanalysis/fig3b_family_compartment_localization_qc.md",
    "results/fig3b_family_compartment_dscore.tsv",
    "results/ncp_v3_singscore_reanalysis/fig3c_fidelity_relay_donor_source.tsv",
    "results/ncp_v3_singscore_reanalysis/fig3c_fidelity_relay_summary.tsv",
    "results/ncp_v3_singscore_reanalysis/fig3c_fidelity_relay_qc.md",
    "results/ncp_v3_singscore_reanalysis/figure3cd_case_study_geometry.tsv",
    "results/ncp_v3_singscore_reanalysis/figure3cd_case_study_axis_contract.tsv",
    "results/ncp_v3_singscore_reanalysis/gse255612_family_donor_compartment_scores.tsv",
    "results/ncp_v3_singscore_reanalysis/gse255612_corrected_family_compartment_permutation_classification.tsv",
    "results/ncp_v3_singscore_reanalysis/fig3d_relay_extension_donor_source.tsv",
    "results/ncp_v3_singscore_reanalysis/fig3d_relay_extension_summary.tsv",
    "results/ncp_v3_singscore_reanalysis/fig3d_relay_extension_qc.md",
    "results/ncp_v3_singscore_reanalysis/fig3cd_family_donor_examples_composite_qc.md",
    "results/ncp_v3_singscore_reanalysis/figure3cd_case_study_composite_geometry.tsv",
    "results/ncp_v3_singscore_reanalysis/figure3_composite_manifest.tsv",
    "results/ncp_v3_singscore_reanalysis/figure3_composite_qc.md",
    "results/fig3_v3_p010_trends.tsv",
    "results/ncp_v3_singscore_reanalysis/gse41177_family_duration_laa_primary.tsv",
    "results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv",
    "results/ncp_v3_singscore_reanalysis/gse41177_family_duration_lapv_sensitivity.tsv",
    "results/ncp_v3_singscore_reanalysis/fig4b_family_duration_correlations.tsv",
    "results/ncp_v3_singscore_reanalysis/fig4c_chronicity_classes.tsv",
    "results/ncp_v3_singscore_reanalysis/fig4c_chronicity_aggregate.tsv",
    "results/ncp_v3_singscore_reanalysis/fig4def_exemplar_scatter.tsv",
    "results/ncp_v3_singscore_reanalysis/fig4def_exemplar_summary.tsv",
    "results/ncp_v3_singscore_reanalysis/fig4_lapv_sensitivity.tsv",
    "results/ncp_v3_singscore_reanalysis/fig4_duration_scores_long.tsv",
    "results/ncp_v3_singscore_reanalysis/figure4_composite_manifest.tsv",
    "results/ncp_v3_singscore_reanalysis/figure4_composite_qc.md",
    "results/fig5b_v3_triangulation.tsv",
    "results/fig5c_v3_cross_axis.tsv",
    "results/fig5d_geomx_roi_correlations.tsv",
    "results/figS1_cohort_overview_summary.tsv",
    "results/figS2_pairing_summary.tsv",
    "results/figS3_ncp_stability_summary.tsv",
    "results/figS4_per_cohort_pca.tsv",
    "results/fig_integrated_af_rewiring_atlas.tsv",
    "results/figure6_meta_axis_loco.tsv",
    "results/figure6_mr_anchors_v3_pipeline_wide.tsv",
    "results/figure6_mr_anchors_v3_consolidated.tsv",
    "results/figure6_mr_coloc_gate_rules.tsv",
    "results/figure6_mr_calibration_genetics.tsv",
    "results/fig6a_v3_locus_plot_windows.tsv",
    "results/external_af_gene_sets.tsv",
    "results/external_af_gene_set_summary.tsv",
    "results/fig6a_anchor_coloc_calibration_source.tsv",
    "results/fig6b_posterior_lead_distance_source.tsv",
    "results/fig6c_family_af_overlap_source.tsv",
    "results/fig6c_family_af_overlap.tsv",
    "results/fig6d_non_nested_traceability_source.tsv",
    "results/fig6d_non_nested_traceability_summary.tsv",
    "results/figure6_composite_manifest.tsv",
    "results/figure6_composite_qc.md",
    "results/figure6_composite_geometry.tsv",
    "results/figure6ab_anchor_row_geometry.tsv",
    "results/figure6cd_composite_geometry.tsv",
    # Legacy G-tier tables retained for provenance only (not main-text routing)
    "results/figure6_g_tier_traceability.tsv",
    "results/figure6_g_tier_rules.tsv",
    "results/ncp_v3_gene_universe.tsv",
    "results/ncp_v3_family_leading_edge_genes.tsv",
    "results/coloc_rerun_diff_report.tsv",
    "results/innovation_ncp_bootstrap_frequency.tsv",
    "results/innovation_ncp_bootstrap_summary.tsv",
    "results/innovation_ncp_v3_axis_families.tsv",
    "results/innovation_ncp_v3_panel.tsv",
    "results/innovation_ncp_v3_summary.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2b_family_donor_source.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2b_family_donor_qc.md",
    "results/ncp_v3_singscore_reanalysis/figS16_program_donor_source.tsv",
    "results/ncp_v3_singscore_reanalysis/figS16_program_donor_qc.md",
    "results/fig1c_ncp_v3_program_plotting_table.tsv",
    "results/figure1_publication_asset_map.tsv",
    "results/figure1_composite_manifest.tsv",
    "results/figure2_composite_manifest.tsv",
    "results/figure2_composite_qc.md",
    "results/reproducibility_table.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2a_input_manifest.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_source.tsv",
    "results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_summary.tsv",
    "results/ncp_v3_family_labels.tsv",
    "results/ncp_v3_selection_criteria.tsv",
    "results/innovation_v2_vs_v3_decision_table.tsv",
    "results/innovation_v2_vs_v3_metrics.tsv",
    "results/methodology_audit_program.tsv",
    "results/v3_primary_replacement_checklist.tsv",
]


def copy_if_exists(src: Path, dst: Path) -> str:
    if not src.is_file():
        return f"MISSING: {src.relative_to(ROOT)}"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return f"OK: {dst.relative_to(DST)}"


def main() -> None:
    for source_docx in (
        ROOT / "manuscript_drafts" / "dxc_benchmark.docx",
        ROOT / "manuscript_drafts" / "Figure_Captions.docx",
        ROOT / "manuscript_drafts" / "Data_Code_Availability.docx",
    ):
        if source_docx.is_file():
            migrate_docx(source_docx)
    if DST.exists():
        shutil.rmtree(DST)
    (DST / "01_Manuscript_Word").mkdir(parents=True)
    (DST / "02_Main_figures").mkdir(parents=True)
    (DST / "03_Supporting_figures").mkdir(parents=True)
    (DST / "04_Supporting_tables_and_drafts").mkdir(parents=True)
    (DST / "05_Key_results_TSV").mkdir(parents=True)
    (DST / "00_SCAN").mkdir(parents=True)

    log_lines: list[str] = []

    m1 = ROOT / "manuscript_drafts" / "dxc_benchmark.docx"
    m2 = ROOT / "manuscript_drafts" / "Figure_Captions.docx"
    m3 = ROOT / "manuscript_drafts" / "Data_Code_Availability.docx"
    log_lines.append(copy_if_exists(m1, DST / "01_Manuscript_Word" / "Manuscript_main.docx"))
    log_lines.append(copy_if_exists(m2, DST / "01_Manuscript_Word" / "Figure_legends.docx"))
    log_lines.append(copy_if_exists(m3, DST / "01_Manuscript_Word" / "Data_Code_Availability.docx"))
    leg = FIG / "fig1_legend.txt"
    log_lines.append(copy_if_exists(leg, DST / "01_Manuscript_Word" / "fig1_legend.txt"))

    scan = ROOT / "submission" / "SCAN_v2_panel_refs_in_docx.txt"
    if scan.is_file():
        shutil.copy2(scan, DST / "00_SCAN" / scan.name)

    for label, path in MAIN_FIGURES:
        out = DST / "02_Main_figures" / f"{label}_{path.name}"
        log_lines.append(copy_if_exists(path, out))

    for label, path in SUPPORTING_FIGURES:
        out = DST / "03_Supporting_figures" / f"{label}_{path.name}"
        log_lines.append(copy_if_exists(path, out))

    drafts = ROOT / "docs" / "manuscript" / "drafts"
    for name in [
        "Supplementary_Table_S1_cohort_demographics.docx",
        "Supplementary_Table_S2_MTC.docx",
        "Supplementary_Table_S3_chamber_reversal_genes.docx",
        "Supplementary_Table_S4_reversal_carrier.docx",
        "Supplementary_Table_S5_data_driven_modules.docx",
        "Supplementary_Table_S6_anchor_matrix.docx",
    ]:
        p = drafts / name
        log_lines.append(copy_if_exists(p, DST / "04_Supporting_tables_and_drafts" / name))

    sub_mech = ROOT / "Manuscript_outputs" / "Supplementary_Table_Sub_Mechanisms.docx"
    log_lines.append(
        copy_if_exists(sub_mech, DST / "04_Supporting_tables_and_drafts" / sub_mech.name)
    )

    for rel in RESULTS_TSV:
        src = ROOT / rel
        dst = DST / "05_Key_results_TSV" / Path(rel).name
        log_lines.append(copy_if_exists(src, dst))

    index = DST / "00_README_SUBMISSION_INDEX.txt"
    lines = [
        "PLOS Computational Biology — submission pack (assembled from repository)",
        "Root: submission/PLOS_CB_submission_pack/",
        "",
        "01_Manuscript_Word/",
        "  Manuscript_main.docx       — full article (source: manuscript_drafts/dxc_benchmark.docx)",
        "  Figure_legends.docx        — optional separate upload if journal splits legends",
        "  fig1_legend.txt            — text legend for Fig. 1 if needed",
        "",
        "02_Main_figures/",
        "  Figure_1 … Figure_7     — one PDF per assembled main figure",
        "  Figure 1 title: Definition and freezing of the healthy paired-chamber NCP coordinate system.",
        "  Figure 1 contains four panels (A-D); prototype and robustness extensions are S13/S14.",
        "  HCA carrier heatmaps are not part of Figure 1.",
        "",
        "03_Supporting_figures/",
        "  Supporting_Figure_S*       — S3/S4 have multiple PDFs (suffix a/b/c)",
        "",
        "04_Supporting_tables_and_drafts/",
        "  Legacy Word tables — cross-check table numbers vs Manuscript_main.docx (S1–S6 v3 text may differ).",
        "",
        "05_Key_results_TSV/",
        "  Tab-separated outputs referenced in the manuscript Methods / legends.",
        "",
        "00_SCAN/",
        "  SCAN_v2_panel_refs_in_docx.txt — automated scan for v2 / panel / results paths in Word",
        "",
        "Journal portal: usually also need Cover letter, CRediT, conflict, funding — not copied here.",
        "See docs/submission/PLOS_CB_SUBMISSION_PACK.md for the authoritative checklist.",
        "",
        "--- copy log ---",
        *log_lines,
    ]
    index.write_text("\n".join(lines) + "\n", encoding="utf-8")
    run_ncp_nomenclature_qc(ROOT)
    print(f"Wrote {index.relative_to(ROOT)} ({len(log_lines)} copy operations)")


if __name__ == "__main__":
    main()
