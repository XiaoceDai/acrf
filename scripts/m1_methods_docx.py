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
        'Ethical approval and public human data',
        [
            'All analyses used publicly available, de-identified human molecular profiles, spatial ROI summaries, and genetic summary statistics. No new participants or biospecimens were enrolled, and no participant-level identifiers were accessed. Each resource was analyzed with the dataset-specific unit appropriate to its design (donor, sample, patient, ROI, or summary-statistic resource); resources were not treated as a single pooled donor denominator.',
        ],
    ),
    (
        'Public datasets and analysis roles',
        [
            'Public atrial transcriptomic, spatial, eQTL, GWAS, and external gene-set resources used in this study are inventoried in S1 Table, which reports accession or resource identifier, platform, tissue or region, chamber information, rhythm or control groups, dataset-specific counts, primary analysis unit, analysis role, and source provenance. Cohort structure for the paired bulk and duration resources is further illustrated in S1 Fig and S2 Fig.',
            'GSE112339 provided the healthy paired left-atrium (LA) and right-atrium (RA) bulk RNA-seq reference: eight healthy donors, each with paired LA and RA samples. This resource defines the frozen healthy chamber direction and program-specific reference magnitude for NCP calibration and is not used to select AF disease signals.',
            'Paired-chamber bulk deployment used three cohorts analyzed in parallel, with no cohort designated as a discovery or replication arm: GSE128188 (10 paired donors; 5 AF and 5 SR), GSE79768 (13 paired donors; 7 AF and 6 SR), and GSE115574 (28 paired donors; 13 AF and 15 SR). Across the three cohorts, the analysis included 51 paired donors, comprising 25 AF and 26 SR donors.',
            'The Human Cell Atlas heart resource supplied healthy LA compartment carrier priors for localization context only. It was not used to select AF disease signals and is not a Path A gate; implementation details are given with the single-nucleus and Path A Methods.',
            'GSE255612 provided LA single-nucleus RNA-seq from 18 AF and 16 non-fibrillating control donors for donor×compartment pseudobulk localization context. The primary localization analysis is LA only. Detailed testing procedures are described with the single-nucleus Methods.',
            'GSE41177 provided the cross-sectional AF-duration resource. The primary duration analysis uses LAA samples from 12 AF donors with numeric duration; LAPV samples are regional sensitivity only. SR samples document cohort structure but are not duration anchors in final inference or main panels. Detailed duration Methods are described later.',
            'GSE238242 was used only for supplementary RNA-side and mural-cell sensitivity (S5 Fig). It is not a Path A gate, is not an independent main-text evidence layer, and does not alter Figure 3 family classes.',
            'GSE261363 provided GeoMx spatial context: 48 ROIs from four patients in the full resource, of which 24 eligible ROIs from four patients enter the final ROI-level descriptive analysis. GeoMx is not a Path A gate; detailed spatial Methods are described later.',
            'Genetic calibration used the Yuan et al. 2025 cross-population AF GWAS (GWAS Catalog accession GCST90624411; corrected publication size 252,438 AF cases and 1,959,739 controls), GTEx left-atrial-appendage cis-eQTL summary statistics through the eQTL Catalogue, and the external AF gene sets used in Figure 6C and S8 Fig. The Yuan GWAS and the Roselli external gene set were retained as distinct prespecified resources. Detailed genetic Methods are described later.',
        ],
    ),
    (
        'Preprocessing and normalization',
        [
            'For the paired bulk RNA-seq resources GSE112339 and GSE128188, gene-level count matrices were filtered with edgeR::filterByExpr, TMM-normalized, and transformed with limma::voom. Duplicate gene symbols were collapsed after voom by retaining the row with maximum variance per symbol.',
            'For the paired Affymetrix resources GSE79768 and GSE115574, deposited GEO Series Matrix probe intensities were mapped to gene symbols and collapsed by selecting the probe with the largest interquartile range per symbol.',
            'For the duration microarray resource GSE41177, CEL files were RMA-normalized with the affy package and collapsed from probe to gene symbol by the maximum-IQR rule.',
            'For GSE255612 single-nucleus localization inputs, deposited expression values were aggregated to donor × compartment pseudobulks by the mean across cells, converted from Ensembl identifiers to gene symbols by mean collapse of duplicates, and transformed to log1p counts per million before program scoring. Detailed localization testing is described with the single-nucleus Methods.',
            'For the supplementary GSE238242 resource, raw counts were summed to donor × compartment pseudobulks using locked compartment labels (atrialCM, fibroblast, endothelial_or_vascular, immune), with mural cells retained as a sensitivity compartment. Within each compartment, pseudobulks were filtered with edgeR::filterByExpr (minimum count 5), TMM-normalized, and transformed to log2-CPM with prior count 1.',
            'For GSE261363 GeoMx spatial context, analyses used the author-provided Q3-normalized whole-transcriptome atlas matrix. Detailed GeoMx scoring, including the family gene-union exception, is described with the GeoMx Methods.',
        ],
    ),
    (
        'Frozen NCP reference construction',
        [
            'The frozen NCP is a single healthy-reference panel comprising 30 pathway programs mapped to 16 families (22 LA-biased and eight RA-biased programs). The panel was frozen before AF-cohort analysis. Disease cohorts did not contribute to candidate scoring, pathway reselection, leading-edge membership, chamber direction, or family assignment.',
            'Healthy-reference selection evaluated 80 candidate pathways. Of these, 38 entered the strict eligibility pool. Diversity-aware pruning retained 25 programs from the strict pool, and five additional programs entered through the prespecified fallback rule to complete the frozen 30-program panel. Fallback admission used the locked fallback consensus and significance thresholds and did not imply that all 30 final programs were drawn from the strict pool.',
            'Prespecified gates included dual-method directional coherence (sign(NES) matching sign(GSVA logFC)), leave-one-donor-out stability (at least seven of eight consistent runs for the strict pool), a consensus-score threshold of at least 0.62 for the strict pool (at least 0.50 for fallback), and dual significance thresholds (strict: adjusted P value < 0.05 and GSVA adjusted P value < 0.05; fallback: adjusted P value < 0.10 and GSVA adjusted P value < 0.10).',
            'Greedy Jaccard redundancy control used a maximum similarity of 0.5 for LA-biased programs and 0.75 for RA-biased programs, with a locked relaxation ladder if needed to complete the panel. Direction-balance targets required a minimum of four and a target of eight programs per direction, and fallback eligibility followed the locked fallback consensus and significance thresholds. Selection-time leading-edge gene sets were frozen at panel lock and reused unchanged across disease datasets; downstream cohorts did not recalculate or replace those leading edges.',
        ],
    ),
    (
        'Exact-reference scoring and family aggregation',
        [
            'Production program scoring used a shared Python exact-reference implementation of the one-direction scoring formula in Bioconductor singscore 1.28.1. Within each sample, measured genes were ranked using the minimum-ties method. For a frozen program with k measured leading-edge genes among N ranked genes, the centered score was',
            '(mean rank − (k + 1)/2)/(N − k) − 0.5,',
            'matching singscore::simpleScore with centerScore = TRUE for a one-direction up-set. Scoring required at least three measured leading-edge genes. An equivalence audit against Bioconductor singscore 1.28.1 compared 8,310 scores and found a maximum absolute difference of 5.55 × 10⁻16. Production scores were generated with the Python implementation rather than by invoking the R package for each sample.',
            'For every sample and program, scoring used the intersection of the frozen selection-time leading-edge set and measured genes; program membership was not re-estimated per cohort. Family aggregation was program-first and equal-weight: each constituent program contributed one score, and the family score was the mean of those program scores. Family scores were not formed by pooling all family genes into one union before scoring for Figures 2–4. The explicit family-gene-union exception used for GeoMx is described with the GeoMx Methods.',
        ],
    ),
    (
        'Donor-level chamber asymmetry and healthy-reference scaling',
        [
            'Because rhythm was constant within donor, disease comparisons were performed on donor-level LA–RA contrasts rather than by treating chamber samples as independent observations. AF-associated chamber remodeling was therefore quantified at the donor level using the relative aligned chamber asymmetry index (aCAIrel). Healthy-reference calibration defines the measurement scale and is not AF evidence.',
            'For donor d and frozen program p, the raw paired contrast is',
            'Δ_{d,p} = score_{d,LA,p} − score_{d,RA,p}.',
            'The frozen healthy-direction sign is s_p = +1 for LA-biased programs and s_p = −1 for RA-biased programs. The aligned chamber asymmetry is',
            'aCAI_{d,p} = s_p × Δ_{d,p}.',
            'The healthy reference magnitude is',
            'M_p = mean over the eight healthy donors of |aCAI_{d,p}|,',
            'that is, the mean of absolute aligned contrasts rather than the absolute value of the mean contrast. The relative index is',
            'aCAIrel_{d,p} = aCAI_{d,p} / M_p.',
            'An aCAIrel value of 1 corresponds to the mean healthy aligned magnitude for that program. Values greater than 1 indicate amplification beyond the healthy magnitude, values between 0 and 1 indicate attenuation while preserving the healthy direction, a value of 0 indicates no aligned chamber asymmetry, and negative values indicate reversal. Program-specific M_p enables comparison across programs with different healthy magnitudes. For paired-bulk family summaries, the production order was: program scores; paired donor contrasts; direction alignment; program-specific scaling; then equal-weight program-first family aggregation of aCAIrel.',
        ],
    ),
    (
        'Paired-chamber bulk remodeling and cross-cohort consistency',
        [
            'Paired-chamber bulk remodeling was analyzed in GSE128188, GSE79768, and GSE115574 (51 paired donors; 25 AF and 26 SR) using frozen NCP definitions and selection-time leading edges without reselection.',
            'Figure 2C reports descriptive family-by-cohort remodeling states. For each family and cohort, mean AF and SR family aCAIrel values were computed. Reversed was assigned first when mean AF aCAIrel was ≤ −0.02. Remaining calls used the AF-to-SR aligned ratio with the SR mean clipped by a denominator floor of 0.02: Attenuated if the ratio was < 2/3, Preserved if 2/3 ≤ ratio ≤ 4/3, and Amplified if the ratio was > 4/3. The classifier was applied to all 16 families in each of the three cohorts, yielding 48 family–cohort calls; runtime assertions verified complete, non-duplicated classification. These labels summarize remodeling patterns and are not independent significance tests or unique causal trajectories.',
            'For each family and cohort, the AF-minus-SR contrast in donor-level family aCAIrel was summarized as Hedges’ g with the standard small-sample correction J = 1 − 3/(4(n_AF + n_SR − 2) − 1), using the pooled within-group standard deviation and sampling variance (n_AF + n_SR)/(n_AF·n_SR) + g²/(2(n_AF + n_SR)). Cohort-specific estimates and normal-approximation 95% confidence intervals were retained for pooling.',
            'Primary family pooling used fixed-effect inverse-variance weighting. For cohort i, the weight was w_i = 1/v_i, where v_i was the sampling variance of Hedges’ g. The pooled estimate was Σw_i g_i/Σw_i, with standard error sqrt(1/Σw_i) and a normal-approximation 95% confidence interval. Between-study variance was not estimated and tau² was recorded as NA. REML with Hartung–Knapp adjustment was restricted to designated archived axis/program-level sensitivity analyses and was not used for any primary family-level estimate.',
            'Consistency classes combined pooled uncertainty and direction agreement. High required that the pooled 95% CI exclude zero and that at least two of three cohort-effect signs match the pooled sign. Moderate required the same sign rule without High status. All remaining families were Low. Consistency classes are evidence summaries and are not causal claims or multiplicity corrections.',
        ],
    ),
    (
        'LA single-nucleus localization and carrier-relative classification',
        [
            'LA single-nucleus localization used GSE255612 atrial snRNA-seq from 18 AF and 16 non-fibrillating control donors, restricted to left-atrial (LA) samples. The analysis unit was donor × compartment pseudobulk across nine compartments: atrial cardiomyocytes, fibroblasts, endothelial cells, macrophages, lymphocytes, pericytes, vascular smooth muscle cells, adipocytes, and neuronal cells. Pseudobulks used donor-aware mean aggregation with log1p counts-per-million representation. Program scores used frozen selection-time leading-edge gene sets and the shared Python exact-reference scorer described above, followed by equal-weight program-first family aggregation.',
            'Program scoring required at least three measured frozen leading-edge genes. Compartments entered localization testing when they contained at least four AF and four control donors and an available exact-reference score. No additional measured-gene threshold was imposed at the localization-testing or class-assignment stage.',
            'For each eligible family × compartment, and for the parallel program-level series, the test statistic was the AF-minus-control difference in mean donor scores. Donor AF/control labels were permuted within compartment, and support was evaluated with a two-sided absolute-value statistic. When the number of distinct label assignments was at most 20,000, all assignments were enumerated and the exact P value was the proportion whose absolute statistic was at least as large as the observed value. Otherwise, 10,000 Monte Carlo label permutations were performed using base seed 20260418 and a deterministic per-test seed. The Monte Carlo exceedance count was initialized to one to account for the observed configuration; for B permutations, P = (1 + number of permuted statistics at least as extreme as observed)/(B + 1). Donors without an available score were omitted from the corresponding compartment contrast. Benjamini–Hochberg q values were calculated within compartment among included tests for reporting, but final class-support calls used nominal permutation support (P < 0.05) and were not FDR-controlled discovery.',
            'Healthy LA carrier compartments were assigned from the Human Cell Atlas heart resource using donor/region/compartment pseudobulk. Full-pathway gene-mean expression was converted to a CPM-log1p scale, row-standardized across reliable compartments, and summarized at family level by equal-weight program aggregation before assigning the dominant LA carrier. A single dominant compartment was required (top standardized value at least 1.0 with gap to the second compartment at least 0.25); otherwise the prior was recorded as unresolved. The carrier prior was used to interpret the relationship between healthy compartment identity and disease localization. It did not enter the AF-minus-control test statistic.',
            'Carrier-relative classes used the production order after testing. Fidelity: the healthy carrier compartment is nominally supported and the Extension rule is not met. Relay: nominal support occurs outside the healthy carrier while the carrier is not supported. Extension: the healthy carrier is supported together with at least one additional compartment that shows the same AF-minus-control sign. Below-detection: no eligible LA compartment meets the preset nominal-support rule. Below-detection means no included LA compartment met the preset support criteria; it does not imply biological absence.',
            'Program-level localization used the same exact-reference and permutation framework. A constituent program may meet nominal localization support even when its parent family remains Below-detection; this does not alter the family class or Path A gate.',
        ],
    ),
    (
        'Cross-sectional AF-duration and chronicity classification',
        [
            'Cross-sectional AF-duration analysis used GSE41177 with frozen exact-reference program scores and equal-weight program-first family activity scores (family activity, not aCAIrel). The primary analysis was restricted to left atrial appendage (LAA) samples from 12 AF donors with finite numeric month-coded duration. Left atrial pulmonary vein (LAPV) samples were analyzed only as regional sensitivity. Sinus-rhythm samples were not used in correlations, fitted trends, class assignment, or final representative panels.',
            'Duration labels of the form XM were parsed to integer months, with 0M retained as zero. Censored or non-numeric values such as >1M were treated as missing for the primary analysis; no unsupported year-to-month conversion was applied. Primary correlations used complete cases with finite duration and finite family activity.',
            'For each family, the primary association was the Spearman correlation between family activity and duration in months, using average ranks for ties and complete observations with finite duration and family activity. At least three complete observations were required. Two-sided permutation P values were obtained by shuffling duration ranks for B = 10,000 iterations using seed 20260418. The exceedance count was initialized to one, giving P = (1 + number of permuted absolute correlations at least as large as the observed absolute correlation)/(B + 1). Donor-level bootstrap resampling with 5,000 replicates from the same seed stream provided percentile 95% confidence intervals for Spearman rho. Asymptotic Spearman P values and Benjamini–Hochberg q values were retained for reporting but did not determine chronicity class.',
            'Chronicity classes used full-precision unrounded |ρ_LAA|: Strong-trend if |ρ_LAA| ≥ 0.4; Weak-trend if 0.2 ≤ |ρ_LAA| < 0.4; and Stable if |ρ_LAA| < 0.2. Permutation P values, bootstrap intervals, LAPV estimates, and regional sign did not determine class.',
            'LAPV sensitivity computed the equivalent family activity–duration Spearman correlation in eligible AF-LAPV samples. Regional sign was classified as Same when both LAA and LAPV correlations were finite and nonzero with matching signs, Opposite when both were finite and nonzero with different signs, and NA when either estimate was zero or missing. Regional sign provided sensitivity context only. Cross-sectional duration supplies chronicity context and does not establish progression, recurrence, treatment response, or temporal causal order.',
        ],
    ),
    (
        'Rule-based cross-layer triangulation and Path A',
        [
            'Cross-layer triangulation evaluated the 16 frozen NCP families in canonical Figure 1D order. No merged latent model was fit; corrected layer summaries retained their own units and uncertainty structures. Healthy HCA LA carrier labels and dominant bulk rewiring states were retained as context-only fields.',
            'Exactly three binary gates defined support. The bulk gate was met when the cross-cohort consistency class was High; the single-nucleus gate was met when the carrier-relative class was anything other than Below-detection; and the duration gate was met when the chronicity class was anything other than Stable. The support count was the number of these three gates that were met, and Path A was assigned when all three gates were satisfied.',
            'HCA carrier status, dominant bulk rewiring state, GeoMx ROI summaries, genetic calibration outputs, GSE238242 sensitivity contrasts, and F4 constituent-program evidence did not enter the support count. Path A is a prespecified rule-based triage label across non-exchangeable evidence layers. It is not a joint statistical test, multiplicity-controlled discovery criterion, or causal-validation claim.',
        ],
    ),
    (
        'GeoMx ROI-level spatial context',
        [
            'GeoMx spatial context used GSE261363. The full resource comprises 48 ROIs from four patients; the final eligible subset comprised 24 cardiac atrial-appendage DNA-positive ROIs from four patients, with ROI as the descriptive analysis unit. Eligibility required atrial-appendage tissue, DNA-positive marking (DNA or DNA+), an available rhythm label, and area restricted to distal or border regions. ROIs lacking required scores were excluded from the corresponding contrast.',
            'GeoMx family signatures used the frozen union of leading-edge genes assigned to each family. This is an explicit exception to program-first aggregation and reflects the reduced spatial expression panel and data structure. GeoMx scores were not numerically merged with transcriptomic family scores. Analyses used the author-provided Q3-normalized whole-transcriptome matrix. Within each ROI, genes were ranked and scored with singscore::simpleScore (knownDirection = TRUE; centerScore = TRUE). Scoring required at least three measured genes from the signature or cell-panel set. Prespecified cell-marker panels were scored with the same rank-based procedure; family–panel gene overlap was handled by independent scoring of each signature rather than by subtracting shared genes.',
            'Prespecified ROI-level summaries comprised two-sided Spearman correlation tests for the F10 family score versus the fibroblast-panel score and the F15 family score versus the macrophage-panel score, together with a two-sided Welch comparison of the F9 family score between AF and SR ROIs. Missing scores were omitted pairwise. Multiple ROIs originated from each patient; inference remained ROI-level without patient-clustered standard errors or mixed-effects models. GeoMx did not enter Path A, and the outputs provide ROI-level spatial context only.',
        ],
    ),
    (
        'Genetic calibration, colocalization, external overlap, and non-nested traceability',
        [
            'Genetic calibration used the Yuan et al. 2025 cross-population AF GWAS (GWAS Catalog accession GCST90624411; 252,438 AF cases and 1,959,739 controls) under frozen harmonized GRCh38 summary-statistic provenance. The Yuan GWAS resource was retained as distinct from the Roselli external gene set.',
            'Ten framework anchors were analyzed: PRELP, OMD, IDH2, NDUFA8, BMP10, SMAD9, SMAD6, NTRK2, NTF3, and FRS2. PITX2 and KCNN3 were included as AF calibration controls. Anchors were prespecified; controls were not selected from the colocalization results; genetics did not enter Path A.',
            'For each anchor, local AF GWAS evidence used a ±1 Mb window about the gene coordinates on GRCh38. Variants were harmonized to the frozen GWAS build, and the lead AF GWAS variant was the first row attaining the minimum P value in the frozen locus slice; exact P-value ties therefore followed the frozen source-row order. Local genome-wide significance used the conventional 5 × 10⁻⁸ threshold for descriptive local-GWAS support flags. Anchors without usable local variants were recorded as not evaluated for that dimension.',
            'Atrial eGene status used the corrected Bonferroni atrial eGene definition. GTEx left-atrial-appendage cis-eQTL summary statistics were obtained through the eQTL Catalogue. For each anchor, the gene-level testing unit was that gene’s cis-eQTL variants in the analysis window; the primary threshold was α = 0.05 / n, where n is the number of cis variants tested for the gene, applied to the lead cis-eQTL P value. The primary atrial eGene call was whether the lead cis-eQTL P value passed this gene-specific Bonferroni threshold. The lead atrial eQTL variant was selected by the same first-minimum-P rule within the frozen gene-specific slice. Raw lead cis-eQTL P values are descriptive only and do not define the primary eGene gate.',
            'Colocalization used coloc.abf on the shared-variant intersection of the local AF GWAS and atrial cis-eQTL summary statistics within the same ±1 Mb window, after allele alignment, with quantitative eQTL and case–control GWAS trait assumptions, eQTL/GWAS sample-size inputs from the summary-statistic fields, and MAF taken from the eQTL table. Analyses required at least 50 shared variants. Priors were p1 = 1 × 10⁻⁴, p2 = 1 × 10⁻⁴, and p12 = 1 × 10⁻⁵. PP4 is the posterior probability of a shared causal variant (hypothesis H4), and PP3 is the posterior probability of distinct local signals (hypothesis H3). Same-variant colocalization support required PP4 ≥ 0.70. PP3-dominant status was assigned when PP3 > PP4 and PP4 < 0.70. Calibration controls were processed with the same pipeline and were not used to redefine thresholds. PP4 is a same-variant colocalization criterion. PP3 dominance is consistent with distinct local signals but does not prove distinct causal variants.',
            'Lead separation was defined as the absolute GRCh38 distance between the frozen lead AF GWAS variant and the frozen lead atrial eQTL variant. Exact P-value ties retained the first matching row in the frozen source order, and missing leads were recorded as not evaluated. Distances were also expressed in megabases. Lead separation was descriptive and did not enter the colocalization posterior.',
            'External AF gene-set overlap used four separate frozen sources: Yuan, Roselli, Selewa, and Choi. Gene symbols were normalized to the frozen NCP-eligible universe of 14,635 genes. For each of the 16 families and four external sets, overlap with frozen family leading-edge genes was tested using a one-sided hypergeometric test. Benjamini–Hochberg correction was applied across all 64 family–set tests, and overlapping genes were retained for display. Member overlap is distinct from family-level enrichment.',
            'Non-nested genetic traceability treated local AF GWAS, atrial eGene, same-variant colocalization, risk/prioritized membership, regulatory annotation, and rare/functional evidence as independent dimensions. Each dimension was recorded as evidence present, evaluated but unsupported, or not evaluated; the corresponding machine-readable values were present, evaluated_no_support, and not_evaluated. The evidence categories are non-nested: no ordinal tier, evidence count, composite score, or causal verdict was formed. Genetics is not a Path A gate, and evaluated non-support does not exclude genetic involvement.',
        ],
    ),
    (
        'Formal supplementary analyses and reporting assets',
        [
            'Formal Supporting Information comprises exactly eight figures (S1–S8 Fig) and exactly six tables (S1–S6 Table). All formal SI items were generated from the frozen final Figures 1–6 source lineage. No formal SI item introduced an additional Path A gate. Archived development assets were retained only in the repository archive and were not part of the formal SI scope.',
            'S1 Fig summarizes dataset-specific cohort structure and analysis units from frozen metadata. It is descriptive only; reported counts are not a pooled donor total across resources. S2 Fig documents donor, chamber, and region availability together with pairing and sample-attribute structure. In GSE41177, sinus-rhythm samples illustrate cohort structure only and do not enter duration correlations or chronicity-class assignment. S4 Fig shows per-cohort principal-component overviews for sample context. For GSE112339, PCA reused the frozen healthy-reference QC coordinates computed with R prcomp on all genes after log counts-per-million transformation with gene centering and scaling. For GSE128188, GSE79768, GSE115574, and GSE41177 (LAA expression), PCA used the top min(2,000, number of nonzero-variance features) features by variance after dropping non-finite values, with feature-mean centering and without unit-variance scaling, via scikit-learn PCA (two components; randomized SVD; random_state = 0); GSE128188 intensities were log2-transformed before feature selection, whereas the Affymetrix cohorts and GSE41177 LAA matrix were analyzed on the stored intensity scale. S4 provides QC and sample context only and does not constitute an independent AF-remodeling test.',
            'S3 Fig expands LA single-nucleus localization without scientific recomputation. S3A displays the complete 16-family × 9-compartment matrix; S3B displays the complete 30-program × 9-compartment matrix; and S3C summarizes corrected family classes with supported compartments. All panels use the final exact-reference scorer lineage and the same donor eligibility, permutation framework, nominal-support rule, and LA-only boundary as Figure 3. No new test or class was introduced by the S3 assembly, and constituent-program values do not redefine family classes. The final S3 composite was assembled from frozen source panels.',
            'S5 Fig reports GSE238242 supplementary sensitivity. Donor × compartment pseudobulks were formed by summing raw counts across cells within each donor and compartment, requiring at least 30 cells; duplicate gene symbols were collapsed by summation. S5A scored six prespecified programs in four primary compartments (atrial cardiomyocytes, fibroblasts, endothelial or vascular cells, and immune cells). S5B combined pericytes and vascular smooth muscle cells into a single mural compartment for the same six programs. Within each analysis matrix, genes were filtered with edgeR::filterByExpr, TMM-normalized, converted to log-CPM, ranked, and scored with singscore::simpleScore (knownDirection = TRUE; centerScore = TRUE); scores with fewer than three matched leading-edge genes were set to missing. The displayed contrast was the AF-minus-SR difference in mean donor scores. Support labels used absolute-delta thresholds of 0.010 and 0.030 with expected-compartment sign rules for display only. S5 is supplementary sensitivity only: it is not one of the three Path A gates, does not change Figure 3 family classes, and is not an independent main-text evidence layer.',
            'S6 Fig decomposes AF-associated family rewiring into chamber-specific contributions. For each family and paired-bulk cohort, the AF-minus-SR difference in the aligned relative family score was calculated separately in LA and RA:',
            'Δ_LA = mean_AF(A_LA) − mean_SR(A_LA),',
            'Δ_RA = mean_AF(A_RA) − mean_SR(A_RA),',
            'where A_LA and A_RA denote the donor-level chamber-specific aligned relative family scores. Families used the corrected equal-weight program-first aggregation and canonical display order. Chamber-source labels were assigned by the frozen multiplicative 1.25 rule on absolute magnitudes (|Δ_LA| ≥ 1.25 × |Δ_RA| for LA-dominant; |Δ_RA| ≥ 1.25 × |Δ_LA| for RA-dominant); otherwise the pattern was classified as bilateral. Because the rule compares absolute magnitudes rather than a reciprocal ratio, a zero opposite-chamber difference does not create a division-by-zero case. No confidence intervals or hypothesis tests were calculated. The decomposition did not validate mechanism direction, define a separate Path A gate, or alter the pooled family estimates in Figure 2.',
            'S7 Fig shows the F4 BCAA constituent-program example under the same corrected exact-reference scorer and program-level Figure 3 permutation framework. Program-level nominal fibroblast support does not redefine the parent family class, does not create a fifth non-Below family, does not enter Path A, and is not FDR-controlled discovery.',
            'S8 Fig displays the complete 16 × 4 family–set overlap matrix corresponding to the same 64 locked comparisons used in Figure 6C. Analyses used the frozen NCP-eligible universe of 14,635 genes, one-sided hypergeometric tests, and Benjamini–Hochberg correction across all 64 tests. No additional tests and no continuous significance scale were introduced. Overlap labels indicate member identity only; member overlap is not family-level enrichment.',
            'Supplementary Tables S1–S6 were built deterministically from frozen final source tables. Each workbook contains a README sheet, visible data sheets, and a Provenance sheet, with machine-readable TSV exports accompanying each visible data sheet. Workbooks contain no hidden analysis formulas; workbook and TSV values were checked for equality; rebuilding tables did not rerun the scientific analyses. S1 Table inventories datasets, metadata, and analysis roles. S2 Table records testing, uncertainty, and evidence-gating strategy. S3 Table reports bulk family and cohort effects with consistency classes. S4 Table reports LA single-nucleus localization and carrier-relative classes. S5 Table reports duration and chronicity. S6 Table reports anchor-level genetic calibration and non-nested traceability.',
        ],
    ),
    (
        'Multiple-testing, uncertainty, and evidence-gating strategy',
        [
            'The study combined heterogeneous analysis layers with different analysis units, estimands, and uncertainty measures. No single global P value or universal multiple-testing correction was applied across all layers. Layer-specific testing, uncertainty, and gating rules are summarized in S2 Table.',
            'Healthy-reference analyses provided calibration and selection provenance and were not an AF hypothesis-testing layer. Figure 2C family-by-cohort states are descriptive and are not independent significance tests. Figure 2D used cohort-specific Hedges’ g with fixed-effect inverse-variance pooling; consistency classes used the pooled confidence interval and cross-cohort direction agreement and are not FDR categories. The High consistency class is one of the three Path A gates.',
            'Figure 3 localization used two-sided nominal permutation support, with exact enumeration where applicable and 10,000 Monte Carlo permutations otherwise. Benjamini–Hochberg q values were reported within compartment, but final carrier-relative classes used nominal P < 0.05 and were not FDR-controlled discovery. Figure 4 chronicity classes used Spearman correlation with permutation and bootstrap uncertainty reporting; Strong-, Weak-, and Stable-trend classes were based only on full-precision absolute ρ_LAA. P values, q values, confidence intervals, LAPV estimates, and regional sign did not define class.',
            'Path A used exactly three rule-based gates: bulk High consistency, single-nucleus non-Below-detection, and duration non-Stable. The support count is not a joint statistical test, and no multiplicity correction is applied to the three-gate rule. Healthy HCA carrier labels, dominant bulk rewiring state, GeoMx ROI summaries, genetic calibration outputs, GSE238242 sensitivity, and F4 constituent-program evidence are not Path A gates.',
            'GeoMx summaries were two-sided nominal ROI-level Spearman and Welch contrasts without patient-clustered or mixed-effects inference and without a multiplicity-controlled discovery claim; GeoMx is not a Path A gate. Genetic calibration separated local AF GWAS genome-wide thresholding, gene-specific Bonferroni atrial eGene calls, colocalization posterior support at PP4 ≥ 0.70 (not a conventional P-value correction), Benjamini–Hochberg correction across the 64 external family–set hypergeometric tests, and non-nested descriptive traceability statuses. Genetics is not a Path A gate, and evaluated non-support does not exclude genetic involvement.',
            'Among formal SI assets, S1–S4 Fig, S6 Fig, and S1–S6 Table primarily expose design, expanded matrices, decomposition, and provenance rather than additional independent discoveries. S5 Fig and S7 Fig retain their explicitly bounded sensitivity and program-level nominal frameworks. S8 Fig displays the same 64 tests as Figure 6C and does not create a second family of tests.',
        ],
    ),
    (
        'Software versions and reproducibility',
        [
            'Primary analyses used R and Python. Frozen production records identified R 4.5.1 with Bioconductor singscore 1.28.1 (matrixStats 1.5.0; GSEABase 1.70.1) for reference-scorer comparisons and R 4.5.1 with edgeR 4.6.3 and limma 3.64.3 for RNA-seq filtering, TMM normalization, and related linear-model workflows. Colocalization used coloc.abf from the R coloc package; the exact production coloc version could not be recovered from the frozen artifacts and is therefore not reported as a production version. Python-based workflows required Python 3.10 or later and used pandas, NumPy, SciPy, matplotlib, and scikit-learn; openpyxl and python-docx supported workbook and manuscript generation. An exact production Python environment lock was not available.',
            'Program and family scores in the corrected production lineage used the shared Python exact-reference implementation. Equivalence to Bioconductor singscore 1.28.1 was verified across 8,310 score comparisons with a maximum absolute difference of 5.55 × 10⁻¹⁶.',
            'Resampling used base seed 20260418 for LA single-nucleus and AF-duration permutations, with deterministic per-test seed derivation where implemented and exact enumeration where the distinct label space was at most 20,000. No broader claim of bit-identical determinism is made for library routines whose implementations are not guaranteed deterministic.',
            'Analyses used frozen selection-time leading-edge gene sets, frozen program and family metadata, the canonical family display order, corrected figure-specific source tables, GWAS provenance for GCST90624411, and the external AF gene-set manifest.',
            'The analysis implementation is distributed as the acrf package (import name crf). Repository, archival, and final release identifiers are reported in Code Availability. Final figures, formal Supporting Information assets, and supplementary workbooks were generated from manifest-listed frozen source tables with accompanying source hashes and QC reports; visible workbook sheets have exact machine-readable TSV exports, and manuscript sections are generated from versioned source scripts.',
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

    def _is_equation_paragraph(text: str) -> bool:
        t = text.strip()
        if not t:
            return False
        if t.startswith('Δ_') or t.startswith('aCAI') or t.startswith('M_p'):
            return True
        if t.startswith('(mean rank') and t.endswith(','):
            return True
        return False

    # ---- Sections ----
    for header, paragraphs in SECTIONS:
        doc.add_heading(header, level=2)
        for p in paragraphs:
            para = doc.add_paragraph(p)
            if _is_equation_paragraph(p):
                para.alignment = WD_ALIGN_PARAGRAPH.CENTER
            else:
                # justify body text for journal-like look
                para.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
