#!/usr/bin/env python3
"""
Build three end-of-study summary tables for the supplement / Methods:

  P-2 Power table              results/power_table.tsv
  P-3 Multi-test correction    results/mtc_summary.tsv
  P-4 Reproducibility table    results/reproducibility_table.tsv

All sources are existing files in `results/` and `scripts/`.  No
re-analysis; this is documentation tabling only.
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
         "10 programs × 3 cohorts (gate per program × cohort)",
         "per-cohort, per-program",
         "permutation; q_perm BH within cohort",
         "q_perm < 0.10 for class assignment",
         "FIGURE2_WORKFLOW.md"),
        ("Fig 3 LA-arm Δsingscore",
         "permutation Δsingscore (AF − Ctrl)",
         "10 programs × 4 main compartments = 40 cells",
         "per program (across the 4 compartments)",
         "BH-FDR within program",
         "q_perm < 0.10 to assign AF major carrier (Unresolved otherwise)",
         "FIGURE3_WORKFLOW.md"),
        ("Fig 4 LAA duration ρ",
         "Spearman ρ + Monte-Carlo permutation",
         "10 programs × 4 statistical conditions "
         "(LAA-primary, LAA-sens-A, LAA-ordinal, LAPV-sens)",
         "region-specific (LAA family separate from LAPV family); "
         "10 programs per family",
         "BH-FDR per region per condition",
         "q_perm < 0.10 for `sig_marker = .`; q_perm < 0.05 for `*`",
         "FIGURE4_WORKFLOW.md"),
        ("Fig 5A triangulation class",
         "rule-based composite (no p-value)",
         "10 programs",
         "n/a (deterministic rule)",
         "n/a",
         "auditable in script header",
         "scripts/fig5a_triangulation.py"),
        ("Fig 5B multiome support gate",
         "Δsingscore (AF − SR) on snRNA pseudobulk",
         "6 programs × 4 main compartments = 24 cells",
         "n/a (gate is a pre-registered direction call)",
         "no FDR (effect-size + direction)",
         "moderate / strong support tier (FIGURE5_WORKFLOW.md §5A)",
         "FIGURE5_WORKFLOW.md"),
        ("Fig 6A axis pooled meta",
         "Hedges' g (AF − SR) per cohort + REML / Hartung–Knapp pool",
         "2 axes (STRUCTURAL, METABOLIC); k = 3 cohorts",
         "no multiple-testing correction "
         "(framework reports g, CI, I², τ², p_HK individually)",
         "Hartung–Knapp variance adjustment",
         "consistency-check framing; not a high-powered pooled discovery",
         "FIGURE6_WORKFLOW.md §3"),
        ("Fig 6B per-program meta",
         "Hedges' g per program + REML / HK pool",
         "10 programs",
         "no MTC (descriptive consistency summary)",
         "Hartung–Knapp",
         "I² + p_HK reported per program",
         "FIGURE6_WORKFLOW.md §5"),
        ("Fig 6C cis-eQTL × AF GWAS coloc",
         "coloc.abf default priors (p1=p2=1e-4, p12=1e-5)",
         "4 anchors × 2 tissues (GTEx Atrial Appendage primary; "
         "GTEx Whole Blood sensitivity)",
         "n/a (small pre-registered anchor set; verdict bins, not BH)",
         "verdict bins: supported / inconclusive / not supported",
         "PP4 ≥ 0.7 supported; PP4 ∈ [0.3, 0.7) inconclusive; "
         "PP4 < 0.3 or cis-eQTL gate p ≥ 5e-8 → not supported",
         "FIGURE6_WORKFLOW.md §6.7"),
        ("Fig S8a size-matched perms",
         "size-matched random gene-set pooled g",
         "6 frozen support programs × 1000 perms each",
         "per program",
         "two-sided empirical p_emp",
         "p_emp reported; no multi-program adjustment "
         "(specificity check, not a discovery screen)",
         "FIGURE_S8_FALSIFICATION_SPEC.md §1"),
        ("Fig S8b axis-membership perm",
         "axis-membership permutation pooled |g|",
         "6 frozen IDs randomly partitioned 3+3 × 1000 perms",
         "per axis",
         "two-sided empirical p_emp",
         "p_emp reported per axis",
         "FIGURE_S8_FALSIFICATION_SPEC.md §2"),
        ("Fig S8c negative controls",
         "pooled g + I² + LAA ρ for 6 pre-registered control pathways",
         "6 controls × 3 statistics",
         "n/a (descriptive comparison)",
         "no MTC; comparison is qualitative",
         "values reported alongside real frozen programs",
         "FIGURE_S8_FALSIFICATION_SPEC.md §3"),
        ("Fig S8d LOPO axis robustness",
         "pooled g for axis with each component program removed",
         "2 axes × 3 LOPO removals",
         "per axis",
         "no MTC (direction-preservation tag)",
         "label: robust / partially-sensitive / single-program-dependent",
         "FIGURE_S8_FALSIFICATION_SPEC.md §4"),
        ("Fig S10 methodology audit",
         "NCP self-recovery; donor-level direction; cross-method ρ; "
         "cross-method axis pooled g",
         "10 programs × 4 tests",
         "n/a (internal consistency)",
         "no MTC",
         "self-recovery + 8/8 direction stability + "
         "cross-method ρ ≥ 0.5",
         "METHODOLOGY_AUDIT_REPORT.md"),
        ("Fig S11 P5 method validity",
         "simulation; bootstrap; classical limma comparison",
         "various",
         "n/a",
         "no MTC; cross-validation framing",
         "see METHOD_VALIDITY_P5_REPORT.md",
         "METHOD_VALIDITY_P5_REPORT.md"),
    ]
    cols = ["figure", "statistic", "count", "family", "correction",
            "threshold_or_call", "spec_source"]
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(RES / "mtc_summary.tsv", sep="\t", index=False)
    print(f"\nwrote results/mtc_summary.tsv  ({len(df)} rows)")


# ---------------------------------------------------------------------------
# B4 Reproducibility table  (figure → script → input → output)
# ---------------------------------------------------------------------------
def build_reproducibility_table() -> None:
    rows = [
        ("Fig 1A", "scripts/fig1a_pca.py",
         "bulk/Project_P01/qc/qc_pca.tsv",
         "figures/fig1a_pca_la_vs_ra.pdf"),
        ("Fig 1B", "scripts/fig1b_volcano_heatmap.py",
         "bulk/Project_P01/normal_chamber_genes_full.tsv + "
         "bulk/logcpm_voom.tsv + bulk/sample_meta.tsv",
         "figures/fig1b_volcano_heatmap.pdf"),
        ("Fig 1C", "scripts/fig1c_bulk_programs.py",
         "normal_chamber_programs.tsv",
         "figures/fig1c_bulk_programs.pdf"),
        ("Fig 1D", "scripts/fig1d_program_attribution.py",
         "results/normal_chamber_program_attribution.tsv",
         "figures/fig1d_program_attribution.pdf"),
        ("Fig 2A", "scripts/fig2a_schematic.py",
         "(schematic, no data input)",
         "figures/fig2a_schematic.pdf"),
        ("Fig 2B", "scripts/fig2b_rewiring_heatmap.py",
         "Project_P1/program_scores_long.tsv",
         "figures/fig2b_rewiring_heatmap.pdf"),
        ("Fig 2C", "scripts/fig2c_rewiring_classes.py",
         "results/af_chamber_rewiring_programs.tsv",
         "figures/fig2c_rewiring_classes.pdf"),
        ("Fig 2D", "scripts/fig2d_representative_dots.py",
         "Project_P1/program_scores_long.tsv",
         "figures/fig2d_representative_dots.pdf"),
        ("Fig 2E", "scripts/fig2e_replication_lollipop.py",
         "results/af_chamber_rewiring_consensus.tsv",
         "figures/fig2e_replication_lollipop.pdf"),
        ("Fig 3A", "scripts/fig3a_schematic.py",
         "(schematic)",
         "figures/fig3a_schematic.pdf"),
        ("Fig 3B", "scripts/fig3b_delta_heatmap.py",
         "results/af_la_program_celltype_delta.tsv",
         "figures/fig3b_delta_heatmap.pdf"),
        ("Fig 3C", "scripts/fig3c_carrier_shift.py",
         "results/figure3_carrier_shift.tsv",
         "figures/fig3c_carrier_shift.pdf"),
        ("Fig 3D", "scripts/fig3d_representative_dots.py",
         "results/af_la_program_celltype_delta.tsv",
         "figures/fig3d_representative_dots.pdf"),
        ("Fig 3E", "scripts/fig3e_anchor_genes.py",
         "Project_P3/GSE255612/pseudobulk/<compartment>/expr.tsv",
         "figures/fig3e_anchor_genes.pdf"),
        ("Fig 4A", "scripts/fig4a_schematic.py",
         "(schematic)",
         "figures/fig4a_schematic.pdf"),
        ("Fig 4B", "scripts/fig4b_duration_summary.py",
         "Project_P4/GSE41177/results/duration_trend_laa_primary.tsv",
         "figures/fig4b_duration_summary.pdf"),
        ("Fig 4C", "scripts/fig4c_representative_scatter.py",
         "Project_P4/GSE41177/program_scores.tsv + sample_metadata.tsv",
         "figures/fig4c_representative_scatter.pdf"),
        ("Fig 4D", "scripts/fig4d_laa_vs_lapv_sensitivity.py",
         "Project_P4/GSE41177/results/duration_trend_*.tsv",
         "figures/fig4d_laa_vs_lapv_sensitivity.pdf"),
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
        ("Fig 6A", "scripts/h2_meta_axis.R + scripts/fig6a_axis_forest.py",
         "results/figure6_axis_composite_long.tsv (from h1_axis_composite.R)",
         "figures/fig6a_axis_forest.pdf + results/figure6_meta_axis.tsv"),
        ("Fig 6B", "scripts/h3_meta_program.R + scripts/fig6b_program_pooled.py",
         "Project_P1/program_scores_long.tsv + figure3_program_panel.tsv",
         "figures/fig6b_program_pooled.pdf + results/figure6_meta_program.tsv"),
        ("Fig 6C", "scripts/h4_mr_coloc_anchors.R + h4b_mr_coloc_blood.R + "
         "scripts/fig6c_mr_plausibility.py",
         "AF GWAS GCST90624411 (tabix HTTP) + "
         "eQTL Catalogue QTD000251/QTD000356 (tabix HTTP)",
         "figures/fig6c_mr_plausibility.pdf + "
         "results/figure6_mr_anchors.tsv + figure6_mr_anchors_blood.tsv"),
        ("Fig 6D", "scripts/fig6d_synthesis.py",
         "(schematic, references Fig 6A/B/C narratively)",
         "figures/fig6d_synthesis.pdf"),
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
