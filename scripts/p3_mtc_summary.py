#!/usr/bin/env python3
"""
B3 — Multi-test correction summary (Supplementary Table S2).

Outputs:
  results/mtc_summary.tsv               (raw 8-column table)
  manuscript_drafts/Supplementary_Table_S2_MTC.docx
                                       (publication-style version)
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from docx import Document

from m0_doc_style import apply_base_style, add_body, add_table


# ---------------------------------------------------------------------------
# Data — single source of truth for MTC strategy across the manuscript.
# Each row = one statistical test family, NOT one p-value.
# ---------------------------------------------------------------------------
COLS = [
    "figure", "panel", "test_family", "test_method",
    "n_tests", "mtc_method", "threshold", "outcome_summary",
]

ROWS: list[list[str]] = [
    # ---- Figure 1 ----
    ["Fig. 1", "b",
     "Paired LA-vs-RA differential expression (GSE112339)",
     "donor-adjusted limma, moderated t",
     "≈ 22,000 genes",
     "Benjamini–Hochberg (per-test, FDR)",
     "FDR < 0.05; |log2FC| ≥ log2(1.5)",
     "Defines the chamber-biased gene pool used to build NCP programs."],

    ["Fig. 1", "c",
     "NCP pathway-program selection (consensus pruning)",
     "dual GSEA + GSVA evidence fusion, LOO stability, Jaccard pruning",
     "80 candidate pathways → 30 deployed programs",
     "Consensus gates + diversity pruning (no FDR)",
     "consensus_score ≥ 0.62; stability ≥ 7/8; dual padj < 0.05",
     "Thirty pathway programs frozen as the NCP panel "
     "(results/innovation_ncp_v3_panel.tsv)."],

    ["Fig. 1", "d",
     "Dual-method prototype panel versus consensus-pruned NCP",
     "side-by-side benchmark on healthy reference and cross-cohort transfer",
     "30 prototype programs vs 30 NCP programs",
     "None (methodological benchmark; no discovery claim)",
     "descriptive continuity / balance / redundancy metrics",
     "Prototype retained only for benchmarking (Supplementary Fig. S13)."],

    ["Fig. 1", "e",
     "NCP selection robustness under evidence-weight perturbation",
     "bootstrap-style resampling of consensus weights / thresholds",
     "30 pathway programs",
     "None (stability panel)",
     "structural preservation of panel size and LA/RA balance",
     "Panel composition remains stable under perturbation "
     "(Supplementary Fig. S14)."],

    ["Fig. 1", "f",
     "Frozen NCP deployment map and HCA LA-carrier priors",
     "family-level carrier assignment from healthy cardiac atlas",
     "16 axis families",
     "None (descriptive deployment panel)",
     "n / a",
     "Sixteen families and LA baseline carriers locked for Figs. 2–5."],

    # ---- Figure 2 ----
    ["Fig. 2", "b–c",
     "Per family × cohort AF rewiring class",
     "donor-level aCAIrel, label-permutation gate",
     "16 axis families × 3 cohorts = 48 family × cohort cells",
     "Empirical permutation p, BH within cohort",
     "qperm < 0.10 to call a strict class",
     "Trend-class overlay shown for context only and not used "
     "for the main-text claim."],

    ["Fig. 2", "e",
     "Replication of bulk rewiring (GSE115574)",
     "donor-level aCAIrel, lollipop direction",
     "16 axis families",
     "None (descriptive replication panel)",
     "Direction concordance with discovery cohorts",
     "Direction-concordance check; not a discovery sample."],

    # ---- Figure 3 ----
    ["Fig. 3", "b–c",
     "AF carrier-fidelity call per NCP axis family (LA-arm, GSE255612)",
     "per family × compartment Δsingscore vs perm null",
     "16 families × 4 compartments = 64 tests",
     "Empirical permutation qperm (BH within family)",
     "qperm < 0.10 to assign Fidelity / Relay / Extension",
     "Most families do not pass and are tagged Below-detection."],

    ["Fig. 3", "e",
     "Anchor-gene support",
     "leading-edge component genes, descriptive",
     "≈ 30 anchor genes",
     "None (descriptive support, not an inference test)",
     "n / a",
     "Anchor genes are descriptive component support; no causal "
     "claim."],

    # ---- Figure 4 ----
    ["Fig. 4", "a–c",
     "LAA AF-only Spearman ρ vs numeric AF duration (12 donors)",
     "Spearman correlation, label permutation (10,000 perms)",
     "16 NCP families (LAA primary)",
     "BH within region (LAA separate from LAPV)",
     "|ρ_LAA| effect-size class (Strong / Weak / Stable)",
     "Descriptive chronicity proxy; SR excluded; not recurrence analysis."],

    ["Fig. 4", "d",
     "Donor-level exemplar scatter (F9, F2, F15)",
     "Descriptive display only",
     "3 families × 12 AF-LAA donors",
     "None",
     "No inferential fit; LAPV ρ shown in facet annotation only",
     "Exemplar panel for narrative linkage to Fig. 2–3."],

    # ---- Figure 5 ----
    ["Fig. 5", "a",
     "Cross-figure triangulation classification",
     "deterministic rule on Fig. 1–4 outputs",
     "16 axis families",
     "None (composite of upstream MTC-controlled tests)",
     "Rule-based class assignment",
     "Each upstream input retains its own MTC; this panel "
     "composites them."],

    ["Fig. 5", "b",
     "Pre-registered multiome support gate (GSE238242)",
     "Δsingscore on 6 support programs × 2 expected compartments",
     "12 pre-registered (program × compartment) cells",
     "Pre-registered gate (no FDR); decision is binary",
     "≥ 4 / 12 cells reaching moderate / strong support",
     "Gate not met; multiome regulator-tier inference deferred "
     "(Supplementary Fig. S5)."],

    ["Fig. 5", "c",
     "ROI-level ECM × fibroblast localisation (GSE261363)",
     "ROI-level Spearman ρ on Q3-normalised WTA matrix",
     "≈ 24 ROIs",
     "None (descriptive localisation panel; no AF / SR primary "
     "test)",
     "n / a",
     "Localisation support only; mechanistic inference relies on "
     "cross-layer synthesis."],

    # ---- Figure 6 ----
    ["Fig. 6", "a",
     "Cross-cohort dual-axis pooled effect",
     "Hedges' g (AF − SR) per cohort, REML + Hartung–Knapp pool",
     "2 axes (k = 3 cohorts each)",
     "None (consistency check across pre-specified axes)",
     "Hartung–Knapp 95 % CI excluding 0",
     "k = 3 → pooled estimate is interpreted as a consistency "
     "check, not a high-powered pooled discovery."],

    ["Fig. 6", "b",
     "Per-program pooled effect (NCP pathway programs)",
     "Hedges' g per cohort, REML + Hartung–Knapp pool",
     "30 NCP pathway programs (k = 3 cohorts each)",
     "BH across the 30 pathway programs (within Fig. 6b family)",
     "Hartung–Knapp 95 % CI; q < 0.10 highlighted",
     "Structural programs the most consistently attenuated; "
     "metabolic programs directionally consistent but more "
     "heterogeneous."],

    ["Fig. 6", "c",
     "Framework-discovered chamber-identity REVERSAL gene class",
     "per-gene aCAI sign-flip with SR-preservation control + "
     "magnitude gate, cross-cohort consistency requirement",
     "12,309 common genes × 3 paired-chamber cohorts",
     "Pre-registered call: AF-specific reversal in ≥ 2 of 3 cohorts "
     "(no FDR; structural call)",
     "|aCAI_healthy| ≥ 0.3 log2 AND |aCAI_AF| ≥ 0.15 log2",
     "24 reversal genes identified (21 RA → LA + 3 LA → RA); "
     "PRELP independently recovered. Positive gene-level evidence "
     "for the distributed substrate model."],

    ["Complementary", "M_aCAI modules",
     "Data-driven rewiring-defined gene modules (PCA + k-means on M_aCAI)",
     "PCA on (gene × donor) aCAI matrix + k-means clustering on "
     "top-10 PC loadings; modules characterised by direction, "
     "ΔaCAI(AF − SR), Jaccard against NCP pathway programs, and "
     "top-loading genes",
     "3,781 genes (filtered) × 59 donors → 12 modules",
     "None (unsupervised discovery; canonical / novel call by "
     "Jaccard ≥ 0.30 against any of the 30 NCP pathway programs)",
     "PC1 = 36.5 %; top-10 PCs = 69.1 % cumulative",
     "12 / 12 modules are novel (no canonical NCP match exceeds "
     "Jaccard 0.30); PITX2 module (M11) collapses in AF, adipose "
     "module (M7) expands toward LA in AF, BMP10 module (M8) "
     "partially equalised — both primitives independently recover "
     "PITX2."],

    ["Fig. 6", "d",
     "Cellular carrier of the chamber-reversal signature "
     "(single-nucleus GSE255612 projection)",
     "per gene per compartment Δ_expr = mean(AF) − mean(Ctrl); "
     "sign-aligned by expected AF direction; reported as "
     "compartment mean Δ_aligned across the 24 reversal genes",
     "24 genes × 4 main compartments = 96 cells",
     "None (descriptive cell-type localisation; no inferential "
     "FDR claim)",
     "Largest mean Δ_aligned identifies the candidate carrier",
     "Fibroblast carries the matrix-invasion signature: mean "
     "Δ_aligned = +0.042 vs +0.009–0.023 in other compartments; "
     "71 % of the 21 RA → LA matrix-invasion genes up-regulated in "
     "AF in fibroblasts."],

    # ---- Supplementary Figures S1–S4 (cohort QC, no inferential MTC) ----
    ["Supp. Fig. S1", "a–d",
     "Cohort overview (donor counts, age, sex, disease background)",
     "descriptive cross-tabulation",
     "8 cohorts × 4 attributes",
     "None (descriptive)",
     "n / a",
     "Reports the study population summarised in Table 1."],

    ["Supp. Fig. S2", "a–d",
     "Donor × chamber / region pairing structure and covariates",
     "descriptive presence / absence + covariate display",
     "4 cohorts × 2 axis levels",
     "None (descriptive)",
     "n / a",
     "Documents the paired structure that motivates the donor-paired "
     "aCAIrel analysis."],

    ["Supp. Fig. S3", "a–c",
     "NCP panel and gene-level stability across donor resamplings",
     "stability_consistent_runs (out of 8) and per-gene "
     "stability_score (out of 8) on GSE112339",
     "30 NCP pathway programs (panel a–b); ≈ 14,600 genes (panel c)",
     "None (consistency check)",
     "Programs retained when stable in 8 / 8 resamplings",
     "All thirty deployed NCP pathway programs are direction-stable "
     "in 8 / 8 resamplings; LE-gene stability is consistently high."],

    ["Supp. Fig. S4", "a–e",
     "Per-cohort top-2,000-feature PCA",
     "PCA on log2-CPM (RNA-seq) or RMA-normalised intensity (Affy)",
     "5 bulk cohorts",
     "None (descriptive QC)",
     "n / a",
     "Donor / chamber dominates the healthy reference; rhythm is a "
     "smaller axis on top of donor variation in the AF cohorts."],

    # ---- Supplementary Figure S5 ----
    ["Supp. Fig. S5", "all",
     "Multiome RNA-side support (full 4-compartment heatmap)",
     "Δsingscore on 6 support programs × 4 compartments",
     "24 cells",
     "None (companion display to the locked gate in Fig. 5b)",
     "Same gate as Fig. 5b",
     "Gate not met; HS-GAG programs trend in the opposite "
     "direction."],

    # ---- Supplementary Figure S6 ----
    ["Supp. Fig. S6", "all",
     "Mural (pericyte + VSMC) compartment sensitivity (GSE238242)",
     "Δsingscore on 6 support programs × mural pseudobulk",
     "6 cells",
     "None (sensitivity panel; not part of the 4-compartment main "
     "analysis)",
     "Direction concordance with main metabolic panel",
     "Directionally consistent with the 4-compartment main panel; "
     "no claim of an independent mural arm."],

    # ---- Supplementary Figure S7 ----
    ["Supp. Fig. S7", "all",
     "Leave-one-cohort-out sensitivity for the dual-axis pooled g",
     "REML + Hartung–Knapp on 2 cohorts at a time",
     "2 axes × 3 LOCO removals = 6 sub-analyses",
     "None (sensitivity check)",
     "Direction stability of structural axis",
     "Structural axis remains direction-preserved across all 3 "
     "LOCO removals."],

    # ---- Supplementary Figure S9 ----
    ["Supp. Fig. S9", "all",
     "Reduced 5-gene panel readout per axis",
     "Spearman ρ between reduced-axis composite and full-axis "
     "composite per donor",
     "six axis-support programs (3 structural + 3 metabolic)",
     "None (framework-feasibility check; not a clinical classifier)",
     "ρ = 0.55–0.90 cohort-level reproducibility",
     "Reduced panel reproduces full-axis composite and AF − SR "
     "centroid direction; reduced metabolic panel does not replace "
     "OXPHOS for chronicity."],

    # ---- Supplementary Figure S8 ----
    ["Supp. Fig. S8", "a",
     "Size-matched random gene-set permutation",
     "permutation null per program (1,000 random sets)",
     "6 support programs",
     "Empirical permutation p (no BH; specificity panel)",
     "Empirical p < 0.10",
     "ECM proteoglycans empirical p = 0.10 (specificity claim "
     "supported jointly by g + I² + LAA-ρ, not by single test)."],

    ["Supp. Fig. S8", "b",
     "Axis-membership permutation",
     "1,000 random 3:3 partitions of the 6 frozen support programs",
     "1 partition test",
     "Empirical permutation p",
     "Empirical p < 0.10",
     "|g_STRUCT| empirical p = 0.33; both axes pool negatively "
     "under any 3:3 partition."],

    ["Supp. Fig. S8", "c",
     "Pre-registered negative-control pathways",
     "axis-style pooled g on 6 control pathways",
     "6 pre-registered controls (3 axis-specificity, 3 unrelated)",
     "None (specificity panel; pre-registered pathway list)",
     "Direction and heterogeneity vs real support programs",
     "Specificity arises from the joint combination of pooled g, "
     "low I², and a coherent LAA-duration ρ — not from any "
     "single statistic."],

    ["Supp. Fig. S8", "d",
     "Leave-one-program-out (LOPO) axis robustness",
     "axis-pooled g after removing each support program in turn",
     "2 axes × 3 LOPO removals = 6 tests",
     "None (binary direction-preservation check)",
     "Direction preserved across all LOPO removals",
     "6 / 6 LOPO removals direction-preserved for both axes."],

    # ---- Supplementary Figure S10 / S11 ----
    ["Supp. Fig. S10", "a–d",
     "Methodology audit (NCP self-recovery, donor stability, "
     "cross-method correlation, cross-method pooled g)",
     "internal-consistency tests (no DE / association call)",
     "4 tests",
     "None (consistency check)",
     "Pre-stated thresholds in scripts/h7_methodology_audit.R",
     "All 4 tests pass; framework internally consistent."],

    ["Supp. Fig. S11", "a–c",
     "Method-validity supplement (simulation power, cohort-donor "
     "bootstrap, classical limma comparison)",
     "external-style validations",
     "3 sub-panels",
     "None (validation, not discovery)",
     "Pre-stated thresholds in scripts/h8_methodology_p5.R",
     "Power ≥ 92 % at true effect = 0.3; bootstrap CIs ⊂ HK CIs; "
     "limma direction concordance 25 / 30 cells."],

    ["Fig. 6c", "all",
     "NCP family × external AF risk-gene overlap (genetic traceability)",
     "One-sided hypergeometric enrichment per family × reference set",
     "16 NCP families × 4 external AF gene sets; background universe = "
     "bulk chamber DE genes (n = 14,635)",
     "BH-FDR across all family×set tests",
     "FDR ≤ 0.05 for enrichment claim",
     "Family-level overlap supports traceability at program level; "
     "not interpreted as single-gene causal validation."],
    ["Fig. 6a / 6b", "all",
     "Single-gene cis-eQTL × AF-GWAS coloc traceability and calibration "
     "check (genetic layer)",
     "coloc.abf (default priors p1 = p2 = 1e-4, p12 = 1e-5)",
     "12 genes × primary tissue (GTEx atrial appendage): 10 pre-registered "
     "framework anchors + PITX2 + KCNN3 calibration genes; optional whole-blood "
     "sensitivity columns where tabulated",
     "No FDR — pre-registered gate-and-coloc",
     "lead cis-eQTL p < 5e-8 AND PP4 ≥ 0.7 to support",
     "No tested anchor reaches the conservative supporting gate (including "
     "PITX2 / KCNN3); absence of PP4 ≥ 0.7 support is calibration-limited at "
     "current GTEx atrial-appendage panel resolution — not informative for "
     "excluding AF-relevant genetics. Anchors not replaced post hoc."],
]


def main() -> None:
    df = pd.DataFrame(ROWS, columns=COLS)

    out_tsv = ROOT / "results" / "mtc_summary.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"wrote {out_tsv.relative_to(ROOT)}  ({len(df)} rows × {len(df.columns)} cols)")

    # ---- Word version ----
    doc = Document()
    apply_base_style(doc)

    doc.add_heading(
        "Supplementary Table S2 | Multiple-testing correction "
        "strategy across figures",
        level=1,
    )

    add_body(
        doc,
        "Each row corresponds to one statistical test family rather than "
        "one p-value. Where an analysis is descriptive or pre-registered "
        "as a binary gate, no formal multiple-testing correction is "
        "applied and the rationale is stated in the outcome column. The "
        "Hartung–Knapp variance adjustment is the default for the "
        "random-effects meta-analyses in Fig. 6 (k = 3 cohorts).",
    )

    header = [c.replace("_", " ") for c in COLS]
    rows_for_doc = [header] + [list(r) for r in ROWS]
    add_table(
        doc, rows_for_doc,
        col_widths_cm=[1.6, 0.9, 3.2, 2.6, 1.8, 2.0, 2.4, 4.0],
    )

    out_doc = ROOT / "manuscript_drafts" / "Supplementary_Table_S2_MTC.docx"
    out_doc.parent.mkdir(parents=True, exist_ok=True)
    doc.save(out_doc)
    print(f"wrote {out_doc.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
