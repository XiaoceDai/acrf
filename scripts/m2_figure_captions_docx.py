#!/usr/bin/env python3
"""Build manuscript_drafts/Figure_Captions.docx (matches Methods.docx style)."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt, RGBColor

from m0_doc_style import (
    apply_base_style,
    add_body,
)

OUT = ROOT / "manuscript_drafts" / "Figure_Captions.docx"
OUT.parent.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Caption content. Each entry: (figure_title, [(panel_letter, body), ...],
# trailing_italic_note_or_None)
# ---------------------------------------------------------------------------
MAIN_FIGURES: list[tuple[str, list[str]]] = [
    (
        "Fig 1. A frozen healthy chamber reference defines the Normal Chamber Prior.",
        [
            (
                "(A) Donor-paired chamber estimand. In paired atrial cohorts, rhythm is assigned "
                "at donor level and chamber is measured within donor. Donor-level LA–RA program "
                "contrasts are aligned to the frozen healthy chamber direction and scaled by the "
                "corresponding healthy reference magnitude to generate the relative aligned "
                "chamber asymmetry index used downstream. This panel defines the measurement "
                "framework and is not an AF association result."
            ),
            (
                "(B) Healthy-reference calibration in eight healthy donors with paired LA and RA "
                "samples from GSE112339. The frozen panel contains 30 programs, of which 22 are "
                "LA-biased and eight RA-biased. All 30 programs matched the canonical chamber "
                "direction in every donor. Healthy reference magnitude varied across programs "
                "(median M = 0.0284; range approximately 0.00940–0.311), motivating "
                "program-specific scaling."
            ),
            (
                "(C) Selection and freeze provenance. Eighty candidate pathways were evaluated; "
                "38 entered the strict pool. Diversity-aware redundancy control retained 25 "
                "strict-selected programs, and five additional programs were admitted under the "
                "prespecified fallback rule to complete the frozen 30-program panel. Disease "
                "cohorts did not contribute to candidate scoring, reselection, program membership, "
                "or chamber direction."
            ),
            (
                "(D) Frozen program-to-family topology. The 30 programs map to 16 NCP families in "
                "the canonical direction-blocked order used throughout subsequent analyses. The "
                "topology retains 22 LA-biased and eight RA-biased programs and spans translation "
                "and ribosome quality control, extracellular-matrix and glycosaminoglycan biology, "
                "growth-factor and receptor signaling, metabolism, and related atrial-remodeling "
                "themes."
            ),
            (
                "LA, left atrium; RA, right atrium; AF, atrial fibrillation; NCP, Normal Chamber "
                "Prior; M, healthy program-specific reference magnitude."
            ),
        ],
    ),
    (
        "Fig 2. Paired-chamber bulk cohorts reveal selective cross-cohort NCP remodeling "
        "in atrial fibrillation.",
        [
            (
                "(A) Bulk deployment of the frozen NCP. Fixed selection-time leading-edge gene "
                "sets were scored at program level in three paired-chamber bulk cohorts. "
                "Sample-level scores were converted to donor-level LA–RA contrasts, aligned to the "
                "healthy chamber direction, scaled by the program-specific healthy magnitude, and "
                "summarized at family level without pathway reselection."
            ),
            (
                "(B) Donor-level family heatmap across 51 paired donors: 25 with AF and 26 in "
                "sinus rhythm. Rows show the 16 frozen families and columns show donors, with "
                "cohort and rhythm annotations. The heatmap displays AF-associated "
                "shifts together with donor, cohort, and family heterogeneity."
            ),
            (
                "(C) Descriptive family-by-cohort remodeling states. Across 16 families and three "
                "cohorts, the 48 family–cohort calls comprised five Amplified, 16 Preserved, 15 "
                "Attenuated, and 12 Reversed states. These labels summarize cohort-level remodeling "
                "patterns and are not independent significance tests or unique causal trajectories."
            ),
            (
                "(D) Cohort-specific and pooled family effects with cross-cohort consistency. "
                "Points and intervals show cohort-specific AF-minus-SR Hedges’ g estimates and "
                "their uncertainty; pooled effects use fixed-effect inverse-variance weighting. "
                "Six families were classified as High consistency, eight as Moderate, and two as "
                "Low. The High set comprised three LA-biased translation/RQC families "
                "(F1, F12, and F14) and three RA-biased matrix/developmental families "
                "(F9, F10, and F15). All six High-consistency pooled effects were negative on the "
                "aligned family scale. Consistency classes summarize pooled-effect uncertainty and "
                "cross-cohort direction agreement and are not causal claims."
            ),
            (
                "AF, atrial fibrillation; SR, sinus rhythm; LA, left atrium; RA, right atrium; "
                "NCP, Normal Chamber Prior; RQC, ribosome quality control; Hedges’ g, standardized "
                "mean difference corrected for small-sample bias."
            ),
        ],
    ),
    (
        "Fig 3. LA single-nucleus data resolve selective NCP family localization.",
        [
            (
                "(A) LA single-nucleus localization design. GSE255612 contributed 18 AF and "
                "16 non-fibrillating control donors. Frozen NCP programs were scored in "
                "donor-aware compartment pseudobulks across nine LA compartments: atrial "
                "cardiomyocytes, fibroblasts, endothelial cells, macrophages, lymphocytes, "
                "pericytes, vascular smooth muscle cells, adipocytes, and neuronal cells. Because "
                "the resource contains LA samples only, this layer provides LA compartment context "
                "rather than an independent test of RA-driven remodeling."
            ),
            (
                "(B) Family-by-compartment AF-minus-control localization map. Cell values show "
                "family-score differences across included compartments. Nominal permutation "
                "markers identify compartments meeting the preset donor-count, detectable-gene, "
                "and permutation-support criteria. These markers are localization cues and are not "
                "FDR-controlled discovery."
            ),
            (
                "(C) Carrier-relative class architecture. Final classes comprised one Fidelity "
                "family (F10 GAG disease / TGF-β), two Relay families "
                "(F11 BMP signaling and F15 NTRK2 / TrkB), one Extension family "
                "(F16 Keratan sulfate), and 12 Below-detection families. F10 showed nominal "
                "Fidelity localization in fibroblasts, F11 Relay localization in adipocytes, F15 "
                "Relay localization in endothelial cells only, and F16 Extension involving "
                "fibroblasts and atrial cardiomyocytes. Below-detection indicates that no included "
                "LA compartment met the preset criteria and does not imply absence of biological "
                "remodeling."
            ),
            (
                "(D) Donor-level examples of Relay and Extension. For F15, the healthy carrier "
                "prior was fibroblast, whereas nominal AF-associated localization was confined to "
                "endothelial cells, illustrating Relay. F16 illustrates Extension, with nominal "
                "support in the fibroblast carrier compartment and atrial cardiomyocytes. The "
                "displayed examples do not expand the four-family non-Below set."
            ),
            (
                "AF, atrial fibrillation; LA, left atrium; RA, right atrium; NCP, Normal Chamber "
                "Prior; FDR, false discovery rate."
            ),
        ],
    ),
    (
        "Fig 4. Cross-sectional AF-duration analysis identifies selective chronicity "
        "trends.",
        [
            (
                "(A) Duration-analysis design in GSE41177. The primary analysis used LAA samples "
                "from 12 AF donors with finite numeric month-coded duration. Because this was not "
                "a paired LA–RA duration design, frozen family activity scores rather than "
                "relative chamber-asymmetry scores were analyzed. LAPV samples provided regional "
                "sensitivity only. Sinus-rhythm samples were not included as duration anchors in "
                "the final main-figure panels."
            ),
            (
                "(B) Family-level LAA duration associations with LAPV sensitivity. Spearman "
                "correlations summarize the relationship between family activity and AF duration "
                "in the 12 primary AF-LAA samples; corresponding LAPV estimates provide a "
                "regional-direction sensitivity view."
            ),
            (
                "(C) Chronicity classes across the 16 families. Classes were assigned from the "
                "full-precision absolute LAA Spearman correlation: Strong-trend for |ρ| ≥ 0.4, "
                "Weak-trend for 0.2 ≤ |ρ| < 0.4, and Stable for |ρ| < 0.2. Two families were "
                "Strong-trend, six Weak-trend, and eight Stable. Display rounding, permutation "
                "p values, bootstrap intervals, LAPV estimates, and regional sign did not determine "
                "class. The classes are descriptive trend-magnitude summaries rather than "
                "multiple-testing-controlled significance calls."
            ),
            (
                "(D) Representative duration patterns. F9 HS-GAG / bone / GPCR showed a Strong "
                "negative LAA trend (ρ = −0.420) and a same-sign negative LAPV estimate "
                "(ρ = −0.161). F2 Mito Ca²⁺ transport showed a Strong negative LAA trend "
                "(ρ = −0.564) but an opposite-sign positive LAPV estimate (ρ = 0.196). "
                "F15 NTRK2 / TrkB showed a Weak negative LAA trend (ρ = −0.385) and a "
                "same-sign negative LAPV estimate (ρ = −0.186). Regional sign is a sensitivity "
                "descriptor and is not interpreted as replication or failure."
            ),
            (
                "Duration provides cross-sectional chronicity context and does not establish "
                "temporal progression, recurrence risk, or treatment response. AF, atrial "
                "fibrillation; LAA, left-atrial appendage; LAPV, left atrial–pulmonary vein "
                "junction; ρ, Spearman correlation coefficient."
            ),
        ],
    ),
    (
        "Fig 5. Rule-based cross-layer triangulation prioritizes two NCP families and "
        "provides ROI-level spatial context.",
        [
            (
                "(A) Cross-layer triangulation contract. Corrected paired-bulk, LA "
                "single-nucleus, and LAA duration summaries were placed on the frozen 16-family "
                "NCP topology without fitting a merged latent model. High bulk consistency, "
                "non-Below-detection LA single-nucleus classification, and non-Stable duration "
                "classification defined the three Path A gates. The healthy HCA LA carrier prior "
                "and dominant bulk state were contextual annotations; GeoMx and genetics were "
                "also excluded from the support count."
            ),
            (
                "(B) Family-level triangulation matrix in the canonical Figure 1D order. Context "
                "columns show the healthy HCA LA carrier prior and corrected dominant bulk state. "
                "Gate columns show bulk consistency, LA single-nucleus class, and LAA duration "
                "class. Support count is the number of satisfied gates, and Path A requires all "
                "three. F10 GAG disease / TGF-β and F15 NTRK2 / TrkB satisfied 3/3 gates and "
                "were assigned Path A. F9 HS-GAG / bone / GPCR satisfied 2/3 gates."
            ),
            (
                "(C) Observed gate architecture. Across 16 families, two satisfied 3/3 gates, one "
                "satisfied 2/3, ten satisfied 1/3, and three satisfied 0/3. F1, F12, and F14 were "
                "bulk-only; F11 and F16 were single-nucleus-only; F2, F3, F7, F8, and F13 were "
                "duration-only; and F4, F5, and F6 met none of the three gates. Path A is a "
                "rule-based triage label across non-exchangeable evidence layers, not a joint "
                "statistical test or causal-validation claim."
            ),
            (
                "(D) GeoMx ROI-level spatial context. The full GSE261363 resource contained "
                "48 ROIs from four patients; the panel displays 24 eligible cardiac "
                "atrial-appendage DNA-positive ROIs. F10 family-signature scores showed a nominal "
                "positive correlation with fibroblast-panel scores "
                "(Spearman ρ = 0.439, nominal p = 0.0318), and F15 scores showed a nominal "
                "positive correlation with macrophage-panel scores "
                "(ρ = 0.431, nominal p = 0.0354). F9 showed a nominal AF-minus-SR ROI-level "
                "difference (Δ = −0.033, Welch p = 0.0481). Marker color and shape identify "
                "patients, and open versus filled markers identify rhythm; axes are panel-specific. "
                "Multiple ROIs arose from each patient, and no patient-clustered or mixed-effects "
                "inference was applied. GeoMx family signatures used the approved frozen "
                "gene-union scoring exception and were not numerically merged with "
                "the program-first family scores used in Figures 2–4. These summaries provide "
                "spatial context only and did not modify Path A."
            ),
            (
                "HCA, Human Cell Atlas; LA, left atrium; LAA, left-atrial appendage; ROI, region "
                "of interest; AF, atrial fibrillation; SR, sinus rhythm; NCP, Normal Chamber Prior."
            ),
        ],
    ),
    (
        "Fig 6. Genetic calibration reveals limited, non-nested traceability of "
        "prespecified NCP anchors and families.",
        [
            (
                "(A) Anchor-level calibration against AF GWAS and atrial cis-eQTL evidence. "
                "Twelve prespecified anchors—10 framework anchors and the PITX2 and KCNN3 AF "
                "calibration controls—were evaluated using the Yuan et al. 2025 cross-population "
                "AF GWAS deposited as GWAS Catalog accession GCST90624411 and GTEx "
                "left-atrial-appendage cis-eQTL data. The corrected GWAS publication reports "
                "252,438 AF cases and 1,959,739 controls. Seven anchor windows contained local "
                "genome-wide AF associations, one anchor, SMAD9, met the corrected atrial eGene "
                "criterion, and no anchor reached the prespecified same-variant colocalization "
                "threshold of PP4 ≥ 0.70. SMAD9 illustrates that atrial eGene evidence does not "
                "imply GWAS–eQTL colocalization. The lack of colocalization support for PITX2 and "
                "KCNN3 under the same regime helps calibrate the overall non-support rather than "
                "exclude gene-level mechanisms."
            ),
            (
                "(B) Continuous colocalization posterior configuration and lead-variant "
                "separation. All "
                "12 anchors were PP3-dominant, and none reached PP4 ≥ 0.70. The adjacent "
                "lollipop track shows the absolute genomic separation between the frozen lead AF "
                "GWAS and atrial eQTL variants, ranging from 0.46 to 1.92 Mb. Lead separation is "
                "a descriptive diagnostic and was not used to calculate the colocalization "
                "posterior. The "
                "combined high-PP3, low-PP4 pattern is consistent with distinct local association "
                "signals rather than supported same-variant bulk-LAA total-expression mediation."
            ),
            (
                "(C) Sparse overlap between frozen NCP family leading-edge genes and four "
                "prespecified external AF gene sets. Sixty-four family–set pairs were tested for "
                "enrichment by a hypergeometric test against a background of 14,635 NCP-eligible "
                "genes. "
                "Five cells contained member overlap, involving two unique genes across three "
                "families. Four cells reflected CAV1 membership in F3 Adipo-lipid axis and "
                "F13 Amino-acid metabolism across two common-risk or prioritized sets; the fifth "
                "reflected HEY2 membership in F9 HS-GAG / bone / GPCR and a regulatory set. "
                "No family–set enrichment survived BH-FDR correction. Repeated CAV1 cells reflect "
                "shared set membership rather than independent genetic observations; member "
                "overlap is therefore shown as a traceability cue rather than evidence of "
                "family-level enrichment."
            ),
            (
                "(D) Non-nested genetic traceability profiles across the 12 anchors. The compact "
                "core profile shows local AF GWAS, atrial eGene, and same-variant colocalization "
                "flags, present for 7/12, 1/12, and 0/12 anchors. External columns show risk or "
                "prioritized gene-set membership, regulatory annotation, and rare or functional "
                "membership, present for 2/12, 2/12, and 0/12 anchors. Filled symbols indicate "
                "evidence present, open symbols indicate that an evaluated criterion was not met "
                "or no membership was found, and dashes indicate evidence not evaluated. The "
                "categories are non-nested and were not combined into an ordinal tier, composite "
                "score, or causal verdict, and they were not used to exclude genetic "
                "mechanisms or as a Figure 5 Path A gate."
            ),
            (
                "AF, atrial fibrillation; GWAS, genome-wide association study; eQTL, expression "
                "quantitative trait locus; PP3, posterior probability favoring distinct local "
                "signals; PP4, posterior probability favoring a shared causal signal; BH-FDR, "
                "Benjamini–Hochberg false discovery rate; LAA, left-atrial appendage; NCP, Normal "
                "Chamber Prior."
            ),
        ],
    ),
]


# Phase 5B4 — final formal SI figure legends (S1–S8 only; no S9+ / no S12)
SUPP_FIGURES: list[tuple[str, list[str]]] = [
    (
        "S1 Fig. Cohort overview and analysis roles of public datasets.",
        [
            (
                "(A) Dataset-specific donor, sample, patient, region-of-interest, or summary-"
                "statistics counts, stratified by rhythm or control group where applicable. "
                "(B) Age distributions by cohort and rhythm or control status where available. "
                "(C) Sex composition by cohort. (D) Reported disease-background annotations. "
                "Counts use the analysis unit appropriate to each resource and should not be "
                "interpreted as a single pooled donor total. The analysis role, primary unit, and "
                "source provenance for each dataset are provided in S1 Table."
            ),
            (
                "AF, atrial fibrillation; SR, sinus rhythm; ROI, region of interest."
            ),
        ],
    ),
    (
        "S2 Fig. Donor pairing and sample-attribute structure.",
        [
            (
                "Pairing and sample-availability map for the three paired-chamber bulk cohorts and "
                "the GSE41177 regional-duration cohort. For GSE128188, GSE79768, and GSE115574, "
                "rows represent donors and columns represent atrial chambers. For GSE41177, rows "
                "represent donors and columns represent left-atrial-appendage and left "
                "atrial–pulmonary vein junction samples. Filled cells indicate available samples; "
                "rhythm, age, and sex annotations are shown where available. GSE41177 "
                "sinus-rhythm samples document cohort structure only and were not used as duration "
                "anchors in the final correlation, classification, or Figure 4D analyses. This "
                "figure documents sample structure and is not an independent AF-remodeling "
                "evidence layer."
            ),
            (
                "AF, atrial fibrillation; SR, sinus rhythm."
            ),
        ],
    ),
    (
        "S3 Fig. Full LA single-nucleus family and program localization.",
        [
            (
                "Donor-aware compartment pseudobulk profiles from GSE255612 were analyzed using "
                "the frozen exact-reference NCP scoring lineage."
            ),
            (
                "(A) Complete 16-family × 9-compartment AF-minus-control localization map. "
                "(B) Complete 30-program × 9-compartment localization map. "
                "(C) Corrected family class and supported-compartment summary."
            ),
            (
                "Final family classes comprised one Fidelity family "
                "(F10 GAG disease / TGF-β), two Relay families "
                "(F11 BMP signaling and F15 NTRK2 / TrkB), one Extension family "
                "(F16 Keratan sulfate), and 12 Below-detection families. Supported compartments "
                "were fibroblast for F10, adipocyte for F11, endothelial only for F15, and "
                "fibroblast plus atrial cardiomyocyte for F16. Markers denote nominal permutation "
                "localization under the preset donor-count and detectable-gene criteria and are not "
                "FDR-controlled discovery. Because the resource contains LA samples only, this "
                "figure supplies LA compartment context rather than an independent test of "
                "RA-driven remodeling. Program-level values do not redefine the family classes; "
                "in particular, any F4 BCAA program signal does not make F4 a fifth non-Below "
                "family."
            ),
            (
                "AF, atrial fibrillation; LA, left atrium; RA, right atrium; FDR, false discovery "
                "rate."
            ),
        ],
    ),
    (
        "S4 Fig. Per-cohort principal-component overview of expression datasets.",
        [
            (
                "Principal-component analysis was performed separately within each expression "
                "cohort using the fixed preprocessing recorded in the figure source lineage. "
                "Panels show the GSE112339 healthy paired LA/RA reference, the three paired-"
                "chamber AF/SR bulk cohorts GSE128188, GSE79768, and GSE115574, and the GSE41177 "
                "duration cohort. Points are annotated by rhythm or control status and by chamber "
                "or region where applicable; paired healthy donors are connected in the reference "
                "panel. This figure provides cohort-level quality control and sample context only "
                "and is not an independent statistical test of AF remodeling."
            ),
            (
                "PCA, principal-component analysis; AF, atrial fibrillation; SR, sinus rhythm; "
                "LA, left atrium; RA, right atrium."
            ),
        ],
    ),
    (
        "S5 Fig. GSE238242 RNA-side and mural-cell sensitivity.",
        [
            (
                "(A) RNA-side AF-minus-SR score differences for six prespecified support programs "
                "across the four primary compartments in the GSE238242 left-atrial-appendage "
                "multiome cohort. (B) Sensitivity analysis after combining pericyte and vascular "
                "smooth muscle cell profiles into a mural-cell compartment. The combined mural-"
                "cell result was directionally consistent with the four-compartment RNA-side "
                "analysis. No expected program–compartment pair met the prespecified moderate or "
                "strong support criterion. These panels are supplementary sensitivity analyses "
                "only: they are not among the three Path A gates, do not alter the Figure 3 family "
                "classes, and do not constitute an independent main-text evidence layer."
            ),
            (
                "AF, atrial fibrillation; SR, sinus rhythm."
            ),
        ],
    ),
    (
        "S6 Fig. Chamber contribution to AF-associated NCP family rewiring.",
        [
            (
                "Chamber-contribution view of AF-associated family shifts across the three "
                "paired-chamber bulk cohorts GSE128188, GSE79768, and GSE115574. For each cohort, "
                "families are plotted by the AF-minus-SR change in relative left-atrial "
                "(ΔaCAI_LA) and right-atrial (ΔaCAI_RA) scores; bubble size encodes "
                "|ΔaCAIrel|, the absolute AF-minus-SR change in the donor-paired relative aligned "
                "chamber asymmetry index. Points are colored by chamber source "
                "(LA-dominant, RA-dominant, or Bilateral), indicating whether the paired-chamber "
                "family shift is driven primarily by the LA component, the RA component, or both. "
                "Family rows follow the final canonical NCP topology and display order. The figure "
                "decomposes chamber contributions to the observed AF-associated family rewiring "
                "and does not constitute mechanism-direction validation or a separate Path A gate."
            ),
            (
                "AF, atrial fibrillation; LA, left atrium; RA, right atrium; NCP, Normal Chamber "
                "Prior; aCAIrel, relative aligned chamber asymmetry index."
            ),
        ],
    ),
    (
        "S7 Fig. F4 BCAA constituent-program fibroblast support.",
        [
            (
                "The parent F4 Mito proteostasis / FA–BCAA family remained Below-detection in the "
                "final family-level LA single-nucleus classification. Within that family, the "
                "branched-chain amino-acid catabolism constituent program showed nominal "
                "fibroblast support under the corrected exact-reference scorer "
                "(Δsingscore = 0.0125; p_perm = 0.0115; n_control = 16; n_AF = 18). The figure "
                "displays the program-level donor and compartment evidence underlying this nominal "
                "result. This is a constituent-program example only: it is not a fifth non-Below "
                "family, does not assign F4 to Fidelity, Relay, or Extension, is not FDR-controlled "
                "discovery, and does not change the Figure 5 support count or Path A classification."
            ),
            (
                "AF, atrial fibrillation; BCAA, branched-chain amino-acid."
            ),
        ],
    ),
    (
        "S8 Fig. Full NCP family overlap with external AF gene sets.",
        [
            (
                "Complete 16-family × 4-set matrix of member overlap between frozen NCP family "
                "leading-edge genes and the prespecified external AF common-risk or prioritized, "
                "regulatory, and rare or functional gene sets. All 64 family–set pairs were tested "
                "against the frozen background universe of 14,635 NCP-eligible genes. Five cells "
                "contained member overlap, involving two unique genes across three families: "
                "CAV1 occurred in F3 Adipo-lipid axis and F13 Amino-acid metabolism across the "
                "Yuan and Roselli sets, and HEY2 occurred in F9 HS-GAG / bone / GPCR and the "
                "Selewa regulatory set. The Choi rare/structural-variant set contained no overlap "
                "with the frozen families. No family–set enrichment survived BH-FDR correction. "
                "Repeated CAV1 cells reflect shared set membership rather than independent genetic "
                "observations. Cell labels indicate member overlap; blank cells indicate tested zero "
                "overlap. Member overlap is not interpreted as family-level enrichment or a Path A "
                "gate."
            ),
            (
                "AF, atrial fibrillation; BH-FDR, Benjamini–Hochberg false discovery rate."
            ),
        ],
    ),
]


# Phase 5B4 — final formal SI table legends (S1–S6 only)
SUPP_TABLES: list[tuple[str, list[str]]] = [
    (
        "S1 Table. Public datasets, cohort metadata, and analysis roles.",
        [
            (
                "The workbook lists the public molecular, spatial, eQTL, GWAS, and external "
                "gene-set resources used in the study. For each resource, it reports accession or "
                "resource identifier, platform, tissue or region, chamber information, rhythm or "
                "control groups, dataset-specific counts, primary analysis unit, analysis role, "
                "main-figure or supplementary link, source publication, source path, and "
                "provenance hash. Counts use the unit appropriate to each dataset—donor, sample, "
                "patient, ROI, or summary-statistics resource—and should not be interpreted as a "
                "single pooled donor total. The workbook contains README, datasets, and "
                "Provenance sheets."
            ),
        ],
    ),
    (
        "S2 Table. Multiple-testing, uncertainty, and evidence-gating strategy.",
        [
            (
                "The workbook summarizes the tested or summarized unit, estimand, statistic, "
                "uncertainty measure, resampling or multiplicity treatment, classification rule, "
                "Path A contribution, and interpretation boundary for each analysis layer. "
                "Exactly three rows contribute to Path A: High bulk cross-cohort consistency, "
                "non-Below-detection LA single-nucleus classification, and non-Stable LAA duration "
                "classification. The healthy carrier prior, dominant bulk state, GeoMx, genetics, "
                "GSE238242 sensitivity, and the F4 BCAA constituent-program example are context or "
                "sensitivity layers and do not enter the Path A support count. The workbook "
                "distinguishes formal effect estimates, nominal localization cues, descriptive "
                "effect-size bins, ROI-level summaries, and non-nested genetic diagnostics. It "
                "contains README, evidence_gating, and Provenance sheets."
            ),
        ],
    ),
    (
        "S3 Table. Corrected bulk paired-chamber effects and cross-cohort consistency.",
        [
            (
                "The family_summary sheet reports the 16 frozen NCP families in canonical order, "
                "including healthy direction, pooled AF-minus-SR Hedges’ g, standard error, "
                "95% confidence interval, fixed-effect inverse-variance pooling metadata, cohort-"
                "sign agreement, and final High, Moderate, or Low consistency class. The final "
                "classes comprise six High, eight Moderate, and two Low families. The "
                "cohort_effects sheet reports the family-by-cohort estimates and the 48 descriptive "
                "state calls: five Amplified, 16 Preserved, 15 Attenuated, and 12 Reversed. State "
                "calls summarize remodeling patterns and are not independent significance tests. "
                "The workbook contains README, family_summary, cohort_effects, and Provenance "
                "sheets."
            ),
        ],
    ),
    (
        "S4 Table. Corrected LA single-nucleus localization and carrier-relative classes.",
        [
            (
                "The family_summary sheet reports the 16 frozen families, healthy HCA LA carrier "
                "prior, final Fidelity, Relay, Extension, or Below-detection class, supported "
                "compartments, and evidence summary. Final classes comprise F10 Fidelity, "
                "F11 and F15 Relay, F16 Extension, and 12 Below-detection families. Supported "
                "compartments are fibroblast for F10, adipocyte for F11, endothelial only for F15, "
                "and fibroblast plus atrial cardiomyocyte for F16. The family_compartment_detail "
                "sheet reports donor counts, detectable leading-edge genes, AF-minus-control family "
                "score differences, permutation method, nominal permutation p values, testing "
                "eligibility, carrier relationship, and localization status for all 16 families "
                "across nine LA compartments. Localization is LA-only and nominal; it is not "
                "FDR-controlled discovery, and Below-detection does not imply biological absence. "
                "The workbook contains README, family_summary, family_compartment_detail, and "
                "Provenance sheets."
            ),
            (
                "HCA, Human Cell Atlas; LA, left atrium; FDR, false discovery rate."
            ),
        ],
    ),
    (
        "S5 Table. Corrected cross-sectional AF-duration trends and chronicity classes.",
        [
            (
                "The family_summary sheet reports the 16 frozen families, full-precision and "
                "displayed LAA Spearman correlations, absolute LAA correlation, final Strong-, "
                "Weak-, or Stable-trend class, permutation p value, bootstrap interval, LAPV "
                "sensitivity estimate, and regional-sign descriptor. Final classes comprise two "
                "Strong-trend, six Weak-trend, and eight Stable families. Classification uses only "
                "the full-precision absolute LAA Spearman correlation; displayed rounding, "
                "p values, confidence intervals, LAPV estimates, and regional sign do not determine "
                "class. The design_summary sheet records the primary 12 AF-LAA samples, duration "
                "parsing, scoring design, LAPV sensitivity role, and the exclusion of SR samples "
                "from duration correlation and class assignment. The analysis provides "
                "cross-sectional chronicity context and is not a progression, recurrence, or "
                "treatment-response model. The workbook contains README, family_summary, "
                "design_summary, and Provenance sheets."
            ),
            (
                "AF, atrial fibrillation; LAA, left-atrial appendage; LAPV, left atrial–pulmonary "
                "vein junction; SR, sinus rhythm."
            ),
        ],
    ),
    (
        "S6 Table. Anchor-level genetic calibration and non-nested traceability.",
        [
            (
                "The anchor_summary sheet reports the 12 prespecified anchors in the final Figure "
                "6 order, including calibration-control role, local AF GWAS status, lead GWAS and "
                "atrial eQTL variants and positions, atrial eGene status, PP3, PP4, same-variant "
                "colocalization status, lead-variant separation, and membership in common-risk or "
                "prioritized, regulatory, and rare or functional external evidence sets. The final "
                "summary contains local AF GWAS support for 7/12 anchors, atrial eGene support for "
                "1/12 (SMAD9), PP3-dominant configuration for 12/12, and PP4 ≥ 0.70 support for "
                "0/12. Risk or prioritized and regulatory membership are each present for 2/12 "
                "anchors, PITX2 and KCNN3, whereas rare or functional membership is present for "
                "0/12. The dimension_summary sheet provides definitions, present counts, thresholds, "
                "source lineage, and interpretation boundaries for the seven evidence dimensions. "
                "The categories are non-nested, are not combined into a score or ordinal tier, and "
                "do not enter Path A. GWAS provenance is locked to Yuan et al. 2025, "
                "GCST90624411, with 252,438 AF cases and 1,959,739 controls. The workbook "
                "contains README, anchor_summary, dimension_summary, and Provenance sheets."
            ),
            (
                "AF, atrial fibrillation; GWAS, genome-wide association study; eQTL, expression "
                "quantitative trait locus; PP3, posterior probability favoring distinct local "
                "signals; PP4, posterior probability favoring a shared causal signal."
            ),
        ],
    ),
]


GRAPHICAL_ABSTRACT = (
    "Three-band schematic. Top: environmental modifiers (air pollution, "
    "noise, heat, built and social environment) as dashed upstream "
    "connectors. Centre: healthy chamber identity → AF chamber rewiring "
    "→ distributed matrix–metabolic substrate. Bottom: clinical "
    "translation tagline — substrate phenotyping beyond a binary AF/SR "
    "state. Dashed arrows denote hypothesised relationships and convey "
    "no causal claim; solid arrows denote relationships established "
    "within this study."
)


def main() -> None:
    doc = Document()
    apply_base_style(doc)

    doc.add_heading("Figure legends", level=1)

    intro = doc.add_paragraph()
    intro.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    intro_run = intro.add_run(
        "All numerical claims in the figure legends below trace to the "
        "tab-separated deliverable files released alongside the "
        "manuscript (see Code availability). Each legend ends with the "
        "pre-registered boundary statement when applicable."
    )
    intro_run.italic = True
    intro_run.font.size = Pt(10)
    intro_run.font.color.rgb = RGBColor(0x66, 0x66, 0x66)

    # Main figures (Phase 5A: six A–D legends; no Figure 7)
    doc.add_heading("Main figures", level=2)
    for title, paragraphs in MAIN_FIGURES:
        doc.add_heading(title, level=3)
        for para in paragraphs:
            add_body(doc, para)

    # Supporting Information figures (Phase 5B4: S1–S8 only)
    doc.add_heading("Supporting Information figures", level=2)
    for title, paragraphs in SUPP_FIGURES:
        doc.add_heading(title, level=3)
        for para in paragraphs:
            add_body(doc, para)

    # Supporting Information tables (Phase 5B4: S1–S6 only)
    doc.add_heading("Supporting Information tables", level=2)
    for title, paragraphs in SUPP_TABLES:
        doc.add_heading(title, level=3)
        for para in paragraphs:
            add_body(doc, para)

    # Graphical abstract
    doc.add_heading("Graphical abstract", level=2)
    add_body(doc, GRAPHICAL_ABSTRACT)

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
