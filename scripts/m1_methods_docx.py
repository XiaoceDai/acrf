#!/usr/bin/env python3
"""
Build the Nature Communications-style Methods Word document.

Output: manuscript_drafts/Methods.docx

Style:
  - Times New Roman 11pt body
  - Arial 11.5pt bold subheadings
  - Double line-spacing
  - Wide margins (2.5 cm)
  - Hanging-indent for the data-availability list

Source: METHODS_TEXT (single Python multiline string composed below).
"""
from __future__ import annotations

from pathlib import Path

from docx import Document
from docx.enum.text import WD_LINE_SPACING, WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "manuscript_drafts"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "Methods.docx"


# ---------------------------------------------------------------------------
# Methods text (publication-ready narrative).
# ---------------------------------------------------------------------------
SECTIONS: list[tuple[str, list[str]]] = [
    (
        "Datasets and cohorts",
        [
            "Eight publicly available human atrial transcriptomic and "
            "multi-omic resources were analysed in this study; no new "
            "sequencing data were generated. Cohort identifiers, "
            "platforms, sample numbers, and analytical roles are "
            "summarised in Supplementary Table 1, and per-cohort "
            "donor counts, age, sex, and disease background are shown "
            "in Supplementary Fig. S1. The donor × chamber (or LAA / "
            "LAPV) pairing structure of the four cohorts that drive "
            "Figs 2 and 4 is shown in Supplementary Fig. S2.",

            "The Normal Chamber Prior was defined in GSE112339 (paired "
            "left and right atrial bulk RNA-seq from 8 healthy donors). "
            "AF chamber rewiring was discovered in two paired-chamber "
            "cohorts—GSE128188 (RNA-seq; 5 sinus rhythm and 5 atrial "
            "fibrillation donors) and GSE79768 (Affymetrix HG-U133 Plus "
            "2.0; 6 sinus rhythm and 7 atrial fibrillation donors)—and "
            "replicated in GSE115574 (Affymetrix HG-U133 Plus 2.0; 16 "
            "sinus rhythm and 15 atrial fibrillation donors with severe "
            "mitral regurgitation as the disease background). The "
            "left-atrial single-nucleus arm was supplied by GSE255612 "
            "(10x snRNA-seq from 34 donors with paired LA + RA "
            "sampling). The chronicity analysis used GSE41177 "
            "(Affymetrix HG-U133 Plus 2.0; 16 valvular AF or sinus "
            "rhythm donors with paired LAA and LAPV sampling and 12 AF "
            "donors with numeric AF duration in months). The "
            "pre-registered multiome support gate used GSE238242 (10x "
            "Multiome of LAA from 7 donors; 3 AF and 4 sinus rhythm). "
            "ROI-level spatial localisation used GSE261363 (NanoString "
            "GeoMx Digital Spatial Profiling of 4 patients across 48 "
            "ROIs in atrial appendage and epicardial adipose tissue).",

            "Genetic plausibility analyses (Fig. 6c) were carried out "
            "against the multi-ancestry atrial fibrillation GWAS "
            "GCST90624411 retrieved from the EMBL-EBI GWAS Catalog, and "
            "against cis-eQTL summary statistics from GTEx v8 Heart "
            "Atrial Appendage (eQTL Catalogue dataset QTD000251; "
            "n = 372). A blood-tissue sensitivity analysis used GTEx v8 "
            "Whole Blood (eQTL Catalogue QTD000356; n = 670) as a "
            "cross-tissue substitute for eQTLGen, which is not "
            "redistributed through the eQTL Catalogue tabix "
            "infrastructure used here.",
        ],
    ),
    (
        "Pre-processing and normalisation",
        [
            "Bulk RNA-seq counts were filtered with edgeR::filterByExpr "
            "(minimum count 5 per group), TMM-normalised, and "
            "transformed to log2-CPM with prior count 1. Affymetrix "
            "microarray CEL files were RMA-normalised with the affy "
            "package and collapsed from probe to gene symbol by "
            "selecting the probe with the largest interquartile range "
            "per symbol. Single-nucleus RNA-seq counts (GSE255612 and "
            "GSE238242) were aggregated to donor × compartment "
            "pseudobulks by summing raw counts; pseudobulks were then "
            "TMM-normalised and transformed to log2-CPM exactly as for "
            "bulk RNA-seq. Cell-to-compartment mapping for GSE238242 "
            "followed RULES.md (atrialCM, fibroblast, "
            "endothelial_or_vascular, immune); pericyte and smooth "
            "muscle cells were retained as a mural-sensitivity "
            "compartment (Supplementary Fig. S6) and not entered in the "
            "primary 4-compartment analysis.",
        ],
    ),
    (
        "Chamber-aware Rewiring Framework",
        [
            "Rhythm is a donor-level attribute in paired chamber "
            "cohorts, so a sample-level donor × chamber × rhythm "
            "interaction term is rank-deficient. We therefore quantified "
            "AF chamber rewiring at the donor level using the Aligned "
            "Chamber Asymmetry Index, aCAIrel.",

            "For each donor d and program p, we computed CAId,p = "
            "score(LA)d,p − score(RA)d,p, where score(·) is the "
            "singscore of the program's leading-edge gene set "
            "(singscore::simpleScore with knownDirection = TRUE and "
            "centerScore = TRUE). The aligned chamber asymmetry was "
            "aCAId,p = sp · CAId,p with sp = +1 for LA-biased and "
            "sp = −1 for RA-biased programs in the Normal Chamber Prior, "
            "and the relative aligned index was "
            "aCAIrel,d,p = aCAId,p / |aCAIbulk_healthy(p)|, where "
            "|aCAIbulk_healthy(p)| was estimated from the GSE112339 "
            "healthy bulk reference. Values of aCAIrel near +1 indicate "
            "preservation of the healthy chamber direction; values "
            "between 0 and 1 indicate attenuation; negative values "
            "indicate reversal. Normalising by the healthy bulk reference "
            "magnitude removes the dependence of singscore magnitudes on "
            "leading-edge set size and therefore allows programs of "
            "different sizes to share a common scale.",

            "The 10 frozen panel programs (program_order locked at "
            "1–10), and the assignment of each program to the structural "
            "or metabolic axis used in Fig. 6, are recorded in the "
            "single-source-of-truth file results/figure5_program_sets.tsv. "
            "Two anchor families were locked before all downstream "
            "analyses: STRUCTURAL_SUPPORT_IDS = {ECM proteoglycans, "
            "HS-GAG Metabolism, HS-GAG Biosynthesis} and "
            "METABOLIC_SUPPORT_IDS = {OXPHOS, Aerobic Resp & ETC, "
            "Cytoplasmic Translation}.",
        ],
    ),
    (
        "Per-figure analyses",
        [
            "Figure 1 establishes the framework's two complementary "
            "discovery primitives on the donor-level chamber asymmetry "
            "matrix M_aCAI (gene × donor; cell value = log2(LA) − "
            "log2(RA) per donor per gene, computed across all 59 "
            "paired-chamber donors of the four cohorts: 8 healthy + 26 "
            "SR + 25 AF; common-gene intersection 12,309). The "
            "curated discovery primitive (Fig. 1a–c) is the Normal "
            "Chamber Prior: paired LA-vs-RA limma analysis of "
            "GSE112339 followed by signed gene-set enrichment over "
            "MSigDB Hallmark and Reactome, retaining the 10 stable "
            "programs that satisfied a cross-resampling stability "
            "threshold (8 / 8 donor resamplings; per-program direction "
            "and leading-edge gene stability are shown in "
            "Supplementary Fig. S3, and the per-cohort PCA structure "
            "of GSE112339 and the four other bulk cohorts is shown in "
            "Supplementary Fig. S4). The data-driven discovery "
            "primitive (Fig. 1d–e) operates on the same matrix without "
            "an external pathway database: M_aCAI was filtered to the "
            "3,781 genes with |aCAI_healthy| ≥ 0.20 log2 (top ~30 % of "
            "chamber-biased genes) and measurable aCAI in at least 30 "
            "donors, gene-wise centred and scaled to unit variance, "
            "and decomposed by principal-component analysis "
            "(scikit-learn 1.5; PC1 = 36.5 %, top 10 PCs = 69.1 % of "
            "variance). Genes were clustered into K = 12 modules by "
            "k-means (k-means++ initialisation, 20 restarts, "
            "random_state = 0) on the 10-PC loading space. Each module "
            "was characterised by its mean aCAI in healthy donors "
            "(signed direction), its mean ΔaCAI between AF and SR "
            "donors (rewiring magnitude), the maximum Jaccard overlap "
            "with the leading-edge set of any of the 10 NCP panel "
            "programs (the canonical threshold was Jaccard ≥ 0.30), "
            "and the genes that contribute most strongly to its "
            "top-3 PC loading. Modules with maximum Jaccard < 0.30 "
            "against any NCP program were called novel "
            "framework-discovered modules. The two primitives are "
            "validated against each other in Fig. 1e: the canonical "
            "GWAS gene PITX2 is independently recovered inside an "
            "unsupervised data-driven module (M11), confirming that "
            "the data-driven primitive does not contradict the "
            "curated primitive on canonical AF biology. "
            "Cell-compartment attribution at single-nucleus resolution "
            "(Fig. 1f) was performed on the GSE255612 snRNA-seq "
            "pseudobulks, rebuilt across the full 13-cell-type "
            "taxonomy of the dataset (atrialCM, fibroblast, vascular "
            "endothelial, lymphatic endothelial, macrophage, "
            "lymphocyte, mast cell, B cell, pericyte, vascular "
            "smooth muscle, adipocyte, epicardial, neuronal; per-donor "
            "pseudobulk samples retained when n_cells ≥ 30). For each "
            "NCP program and each data-driven module, we computed the "
            "mean expression of the gene members in each compartment "
            "and z-scored the resulting vector across compartments; "
            "the compartment with the largest row z-score was reported "
            "as the candidate single-cell carrier. The original "
            "4-compartment lock from RULES.md is preserved as the set "
            "used by Fig. 3 (which was pre-registered on that lock); "
            "the extended 13-compartment scheme is used in Fig. 1f "
            "specifically because the data-driven primitive of "
            "Fig. 1d–e predicted carrier compartments (adipocyte and "
            "epicardial cells) that the original lock had excluded.",

            "Figure 2. AF chamber rewiring was quantified at the donor "
            "level by aCAIrel and classified per program × cohort into "
            "Preserved, Amplified, Attenuated, Reversed, or Low-support "
            "by a permutation-tested |ΔaCAIrel| gate (FIGURE2_WORKFLOW.md). "
            "Discovery used GSE128188 and GSE79768; replication used "
            "GSE115574 (Fig. 2e). Permutation gates and trend-class "
            "overlays are pre-specified.",

            "Figure 3. GSE255612 cells were aggregated by donor × "
            "compartment with raw-count summation. AF major carrier per "
            "program was assigned only when one of the four main "
            "compartments passed qperm < 0.10; otherwise the program was "
            "tagged Unresolved. Anchor genes shown in Fig. 3e are "
            "component genes drawn from each program's leading-edge set "
            "and are reported as descriptive support (no causal-driver "
            "claim).",

            "Figure 4. GSE41177 LAA samples were singscored for each of "
            "the 10 frozen panel programs and Spearman ρ with numeric AF "
            "duration in months was computed on AF samples only "
            "(n = 12). Sinus rhythm samples were displayed as visual "
            "anchors at 0 months and excluded from the primary "
            "statistical model. Sensitivity analyses recoded the \">1 "
            "month\" entries as 1 month or used a three-tier ordinal "
            "duration coding; LAPV samples were carried as a region "
            "sensitivity panel only (Fig. 4d). Multiple-testing "
            "correction was per region (LAA family separate from LAPV "
            "family), with BH-adjusted permutation q values.",

            "Figure 5. Cross-figure triangulation (Fig. 5a) integrated "
            "the per-program healthy bias, bulk rewiring class, AF major "
            "carrier (or Unresolved), and LAA duration ρ into a "
            "rule-based Triangulation class (Convergent, Structure-led, "
            "Metabolism-led, Fragmented), encoded in the script header "
            "of scripts/fig5a_triangulation.py for auditability. The "
            "pre-registered multiome support gate (Fig. 5b) tested the "
            "six axis-support programs in their expected compartments "
            "on the GSE238242 LAA snRNA pseudobulk; the gate was not "
            "met (Supplementary Fig. S5) and the regulator-tier inference "
            "was therefore deferred. ROI-level localisation (Fig. 5c) "
            "computed singscore on the GSE261363 GeoMx Q3-normalised "
            "WTA matrix for an ECM proteoglycans module and a canonical "
            "fibroblast marker panel; ROI-level Spearman ρ is reported "
            "as descriptive (no AF / SR primary statistic). Figure 5d "
            "is a schematic.",

            "Figure 6. For Fig. 6a, donor-level axis composites were "
            "computed as the equal-weight mean of aCAIrel over the three "
            "support programs in each axis; Hedges' g (AF − SR) was "
            "computed per cohort with the standard small-sample "
            "correction; effect sizes were combined across cohorts by "
            "REML random-effects meta-analysis with the Hartung–Knapp "
            "variance adjustment (metafor 4.8.0). With k = 3 cohorts the "
            "pooled estimate is interpreted as a consistency check "
            "across discovery cohorts, not a high-powered pooled "
            "discovery. Figure 6b applied the same recipe at the "
            "individual-program level for the 10 frozen panel programs.",

            "Figures 6c–d report the framework-discovered chamber-"
            "identity REVERSAL gene class together with its single-cell "
            "compartment carrier signature. Figure 6c reports the "
            "reversal gene class itself. A reversal gene is defined per gene as "
            "a gene whose chamber direction in healthy donors (LA-biased "
            "or RA-biased) is flipped in AF donors of paired chamber "
            "cohorts, while sinus-rhythm donors of the same cohorts "
            "preserve the healthy direction. Per cohort, per donor, we "
            "computed aCAI per gene as log2(LA) − log2(RA); per cohort, "
            "we computed AF_dir(g) and SR_dir(g) as the sign of the "
            "donor-mean of aCAI in AF and SR donors respectively, and "
            "Healthy_dir(g) as the sign of the donor-mean across the "
            "eight GSE112339 healthy paired donors. A gene was called "
            "AF-specific reversed in a cohort when Healthy_dir × AF_dir "
            "< 0, Healthy_dir × SR_dir > 0, |aCAI_healthy| ≥ 0.3 log2 "
            "(top ~15 % of chamber-biased genes in the healthy "
            "reference), and |aCAI_AF| ≥ 0.15 log2 (~80th percentile of "
            "|aCAI_AF|). The final reversal gene list required "
            "AF-specific reversal in at least 2 of the 3 paired chamber "
            "cohorts (GSE128188, GSE79768, GSE115574). For Figure 6d "
            "(carrier compartment of the reversal signature), each of "
            "the 24 reversal genes was projected onto the four main "
            "GSE255612 single-nucleus compartments (atrialCM, "
            "fibroblast, endothelial / vascular, immune; RULES.md). "
            "Per gene per compartment we computed Δexpr = "
            "mean(expression in AF biosamples) − mean(expression in "
            "Ctrl biosamples), then sign-aligned by the gene's expected "
            "AF direction (positive Δ for RA → LA reversal genes, "
            "negative Δ for LA → RA reversal genes) so that a positive "
            "Δ_aligned indicates that the gene's expression shifts in "
            "the direction consistent with its chamber-reversal call. "
            "The compartment with the largest mean Δ_aligned across the "
            "24 reversal genes was reported as the candidate cellular "
            "carrier of the reversal signature.",

            "A complementary, classical-regime check by single-gene "
            "cis-eQTL × AF-GWAS colocalisation is reported in "
            "Supplementary Fig. S12. The same four pre-registered "
            "anchor genes (PRELP, OMD, IDH2, NDUFA8), together with two "
            "pre-registered canonical AF GWAS positive controls "
            "(PITX2 at the 4q25 locus and KCNN3 at 1q21), were tested "
            "by ±1 Mb tabix slicing of the harmonised AF GWAS "
            "(GCST90624411) and the GTEx Atrial Appendage cis-eQTL "
            "dataset (QTD000251), allele-alignment, and coloc.abf with "
            "default priors (p1 = p2 = 1e-4, p12 = 1e-5). Verdict bins "
            "were locked before analysis: an anchor was supported when "
            "its lead cis-eQTL p in the primary tissue was below 5e-8 "
            "and PP4 was at least 0.7; inconclusive when the gate was "
            "met and PP4 was between 0.3 and 0.7; and not supported "
            "otherwise. Anchors that failed the gate were not replaced "
            "after the fact. A blood-tissue sensitivity analysis (GTEx "
            "Whole Blood, QTD000356) was carried out with the same "
            "pipeline.",
        ],
    ),
    (
        "Statistical methods",
        [
            "Effect sizes followed the Hedges' g small-sample-corrected "
            "standardised mean difference, with the AF − SR sign "
            "convention so that negative g indicates AF lower than SR "
            "(attenuation of the healthy axis direction). Random-effects "
            "meta-analysis used the REML estimator with the Hartung–"
            "Knapp variance adjustment as the default, implemented in "
            "metafor 4.8.0. Multiple-testing correction families and "
            "thresholds are specified per figure in Supplementary "
            "Table S2 (mtc_summary.tsv) and in the per-figure workflow "
            "documents. Colocalisation analyses used coloc 5.2.3 with "
            "default priors (Fig. 6c). All software versions are listed "
            "in the Code Availability section.",
        ],
    ),
    (
        "Methodology validation",
        [
            "We verified the framework's pipeline by four "
            "internal-consistency tests (Supplementary Fig. S10): NCP "
            "self-recovery (mean aCAIrel ≈ 1 for all 10 panel programs "
            "on the healthy bulk reference); donor-level direction "
            "stability (8/8 healthy donors with aCAI > 0 for every one "
            "of the 10 programs); cross-method per-donor correlation "
            "between singscore and a mean-z scoring (Spearman ρ between "
            "0.52 and 0.93 across programs); and cross-method robustness "
            "of the dual-axis pooled g (the structural axis remained the "
            "more stable arm under either scoring method, with pooled "
            "g = −0.57 / −0.79 and I² = 0%).",

            "We further verified analytical validity (Supplementary "
            "Fig. S11) by a cohort-matched simulation in which the "
            "pipeline recovered known axis effects with at least 92% "
            "detection at a true effect size of 0.3 in aCAIrel units; a "
            "cohort-donor bootstrap (B = 2000) of the dual-axis pooled "
            "g whose 95% intervals were consistent with and tighter than "
            "the Hartung–Knapp intervals reported in Fig. 6a; and a "
            "classical limma per-chamber AF-vs-SR analysis whose "
            "per-program direction agreed with the framework in 25 of "
            "30 program × cohort cells (83%, 100% in GSE115574). "
            "Specificity of the dual-axis story was further checked by "
            "a four-component falsification supplement (Supplementary "
            "Fig. S8) comprising 1,000 size-matched random gene-set "
            "permutations per support program, a 1,000-iteration "
            "axis-membership permutation, six pre-registered "
            "negative-control pathways (axis-specificity controls "
            "INFLAMMATORY_RESPONSE, KRAS_SIGNALING_UP, E2F_TARGETS, and "
            "biologically unrelated controls OLFACTORY_SIGNALING, "
            "KERATINIZATION, SENSORY_PERCEPTION_OF_TASTE), and a "
            "leave-one-program-out (LOPO) axis-attenuation test in "
            "which both axes remained direction-preserved across all "
            "six removals.",

            "Where individual control pathways achieved comparable "
            "single-statistic magnitudes to the real frozen support "
            "programs (for example, REACTOME_KERATINIZATION pooled g = "
            "−1.05), they were accompanied by high cross-cohort "
            "heterogeneity (I² = 78%) and a much weaker LAA-duration ρ "
            "(|ρ| ≤ 0.29), so the framework's specificity is in the "
            "joint combination of pooled effect direction, low "
            "cross-cohort heterogeneity, and a coherent LAA-duration "
            "association rather than in any single statistic.",
        ],
    ),
    (
        "Sex-stratified sensitivity",
        [
            "Per-donor sex was unavailable in the GEO record of "
            "GSE115574; GSE128188 contained only male donors; only "
            "GSE79768 contained both male (n = 5) and female (n = 8) "
            "donors with sufficient AF and SR representation for "
            "sex-stratified analysis. Within GSE79768 the structural "
            "axis was directionally consistent in both sexes (Hedges' g "
            "of −0.76 in male donors and −0.27 in female donors), "
            "whereas the metabolic axis was too small to call (g = "
            "+1.29 in 3 + 2 male donors and −0.40 in 4 + 4 female "
            "donors). Cross-cohort sex-stratified meta-analysis was "
            "therefore underpowered (k ≤ 2 per stratum) and is not "
            "presented in the main text; per-cohort, per-stratum values "
            "are reported in Supplementary Table "
            "(sex_sensitivity_axis_long.tsv).",
        ],
    ),
    (
        "Software, data, and code availability",
        [
            "All computational analyses were carried out in R 4.4 "
            "(packages: singscore, metafor 4.8.0, limma, edgeR, coloc "
            "5.2.3, Rsamtools, data.table) and Python 3.12 (pandas, "
            "numpy, scipy, matplotlib). Tabix HTTP slicing of remote "
            "GWAS and eQTL summary statistics was performed via "
            "Bioconductor Rsamtools::TabixFile. All analysed datasets "
            "are public; accession identifiers are listed in "
            "Supplementary Table 1 and reproduced in the Data "
            "Availability statement. All analysis code, intermediate "
            "tables, and figure-generation scripts are available at the "
            "project repository (URL to be assigned) and will be "
            "permanently archived on Zenodo upon publication "
            "(DOI to be assigned).",
        ],
    ),
]


# ---------------------------------------------------------------------------
# Build the document.
# ---------------------------------------------------------------------------
def main() -> None:
    doc = Document()

    # ---- Page setup ----
    for section in doc.sections:
        section.top_margin    = Cm(2.5)
        section.bottom_margin = Cm(2.5)
        section.left_margin   = Cm(2.5)
        section.right_margin  = Cm(2.5)

    # ---- Default body style ----
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(11)
    style.font.color.rgb = RGBColor(0x10, 0x10, 0x10)
    pf = style.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    pf.space_before = Pt(0)
    pf.space_after = Pt(6)

    # Heading 1 — section title only
    h1 = doc.styles["Heading 1"]
    h1.font.name = "Arial"
    h1.font.size = Pt(13)
    h1.font.bold = True
    h1.font.color.rgb = RGBColor(0x1A, 0x1A, 0x1A)

    # Heading 2 — subsection headers
    h2 = doc.styles["Heading 2"]
    h2.font.name = "Arial"
    h2.font.size = Pt(11.5)
    h2.font.bold = True
    h2.font.color.rgb = RGBColor(0x1A, 0x1A, 0x1A)
    pf2 = h2.paragraph_format
    pf2.space_before = Pt(14)
    pf2.space_after = Pt(4)

    # ---- Title ----
    title = doc.add_heading("Methods", level=1)
    title.alignment = WD_ALIGN_PARAGRAPH.LEFT

    # ---- Sections ----
    for header, paragraphs in SECTIONS:
        doc.add_heading(header, level=2)
        for p in paragraphs:
            para = doc.add_paragraph(p)
            # justify body text for journal-like look
            para.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    # ---- Trailing note (kept off the publication-ready text) ----
    note = doc.add_paragraph()
    note.paragraph_format.space_before = Pt(18)
    run = note.add_run(
        "End-of-Methods note (not for typesetter): the full "
        "single-source-of-truth specifications, frozen anchor lists, "
        "permutation details, and per-figure deliverable paths are "
        "released alongside the manuscript at the project repository "
        "(see Code availability)."
    )
    run.italic = True
    run.font.size = Pt(9)
    run.font.color.rgb = RGBColor(0x66, 0x66, 0x66)

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
