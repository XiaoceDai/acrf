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
    add_caption_paragraph,
)

OUT = ROOT / "manuscript_drafts" / "Figure_Captions.docx"
OUT.parent.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Caption content. Each entry: (figure_title, [(panel_letter, body), ...],
# trailing_italic_note_or_None)
# ---------------------------------------------------------------------------
MAIN_FIGURES: list[tuple[str, list[tuple[str, str]], str | None]] = [
    (
        "Figure 1 | Framework construction: dual discovery primitives on the "
        "donor-level chamber asymmetry matrix",
        [
            ("a", "Principal-component analysis of GSE112339 healthy bulk "
             "RNA-seq (n = 16 samples; 8 paired donors). Donor-level "
             "variation dominates the global transcriptome, motivating the "
             "paired chamber-aware analysis used downstream and the "
             "framework's atomic unit M_aCAI(d, g) = log2(LA donor d, "
             "gene g) − log2(RA donor d, gene g)."),
            ("b", "Volcano plot and heatmap of paired LA-vs-RA differential "
             "expression (donor-adjusted limma). Structured chamber-biased "
             "genes constitute the substrate for both discovery primitives "
             "below."),
            ("c", "Curated discovery primitive — Normal Chamber Prior (NCP): "
             "10 stable LA- or RA-biased gene programs derived by signed "
             "GSEA over MSigDB Hallmark and Reactome. LA-biased programs "
             "are metabolic and translational; RA-biased programs are "
             "matrix and glycan. Program order is locked as program_order "
             "1–10 and used unchanged downstream as the calibration set."),
            ("d", "Data-driven discovery primitive — gene loadings on the "
             "top-2 principal components of M_aCAI (3,781 chamber-biased "
             "genes × 59 paired chamber donors across 4 cohorts; gene-wise "
             "centred and scaled), coloured by k-means module assignment "
             "(k = 12). PC1 explains 36.5 % of the chamber-asymmetry "
             "variance; the top-10 PCs cumulatively explain 69.1 %. All "
             "12 modules are framework-discovered novel modules (Jaccard "
             "≥ 0.30 against any NCP program is not reached by any of "
             "the 12 modules)."),
            ("e", "Module × donor heatmap of mean module aCAI (donors "
             "ordered healthy → SR → AF; modules ordered LA-biased → "
             "RA-biased and within direction by |ΔaCAI(AF − SR)|). The "
             "PITX2-anchored mesothelial module (M11; ITLN1, MSLN, KRT19, "
             "PRG4, SLPI, HP) collapses toward RA in AF "
             "(ΔaCAI(AF − SR) = −0.33); the adipose module (M7; ADIPOQ, "
             "LEP, PLIN1, PCK1, CIDEC, RBP4, TF) expands toward LA in AF "
             "(ΔaCAI(AF − SR) = +0.53); the small atrial-cardiomyocyte "
             "module (M8; HAMP, BMP10) is partially equalised toward LA "
             "(ΔaCAI(AF − SR) = +0.50). The independent recovery of "
             "PITX2 inside the unsupervised M11 cluster validates the "
             "data-driven primitive against canonical AF biology."),
            ("f", "Cell-compartment attribution at single-nucleus "
             "resolution across the full 13-cell-type taxonomy of "
             "GSE255612 (atrialCM, fibroblast, vascular endothelial, "
             "lymphatic endothelial, macrophage, lymphocyte, mast "
             "cell, B cell, pericyte, vascular smooth muscle, "
             "adipocyte, epicardial, neuronal). Left, 10 NCP programs "
             "(curated primitive). Right, 12 data-driven modules "
             "(data-driven primitive). Each cell shows the row "
             "z-score of mean program / module gene expression across "
             "compartments; the compartment with the largest row "
             "z-score is circled. Green box = predicted single-cell "
             "carrier confirmed: the M7 adipose module peaks in "
             "adipocytes (z = +3.45), and the M11 PITX2 / mesothelial "
             "module peaks in epicardial cells (z = +3.46). Both "
             "carrier compartments had been excluded from the "
             "original 4-compartment lock used by Fig. 3 and are "
             "added here precisely because the data-driven primitive "
             "predicted them."),
        ],
        "Cell-compartment attribution in panel f is descriptive and is not "
        "interpreted as a causal-driver claim. The original "
        "4-compartment lock (atrialCM, fibroblast, "
        "endothelial / vascular, immune) is preserved for the "
        "pre-registered Fig. 3 analysis; the extended 13-compartment "
        "scheme is used in panel f only.",
    ),
    (
        "Figure 2 | Atrial-fibrillation chamber rewiring",
        [
            ("a", "Analytical schematic of the Chamber-aware Rewiring "
             "Framework. Because rhythm is a donor-level attribute in "
             "paired chamber cohorts, a sample-level donor × chamber × "
             "rhythm interaction is rank-deficient; rewiring is therefore "
             "quantified at the donor level using "
             "aCAIrel = aCAI / |aCAI healthy_bulk(p)|, which aligns each "
             "program to its NCP direction and normalises by the healthy "
             "bulk magnitude."),
            ("b", "AF chamber-rewiring landscape (mean aCAIrel; programs × "
             "cohort × rhythm) across paired atrial cohorts. A diverging "
             "PuOr palette is used: positive values indicate preserved or "
             "amplified healthy direction; negative values indicate "
             "reversal."),
            ("c", "Per-program rewiring class (Preserved, Attenuated, "
             "Reversed, Amplified, or Low-support) by cohort, called by a "
             "strict permutation gate. The trend-class overlay is shown "
             "for context only and is not used for the main-text claim."),
            ("d", "Representative donor-level aCAIrel distributions for "
             "example programs."),
            ("e", "Replication in the GSE115574 mitral-regurgitation "
             "cohort. The replication cohort provides directional support "
             "and is not the discovery sample."),
        ],
        None,
    ),
    (
        "Figure 3 | LA-arm cell-compartment localisation",
        [
            ("a", "Pseudobulk schematic. GSE255612 cells were mapped to "
             "four main compartments (atrialCM, fibroblast, "
             "endothelial_or_vascular, immune; locked cell-type buckets) and aggregated by "
             "donor × compartment."),
            ("b", "AF-associated LA-arm Δsingscore is modest in absolute "
             "magnitude and distributed across compartments (per program × "
             "compartment Δsingscore (AF − Ctrl))."),
            ("c", "Carrier-shift table: dominant carrier in Fig. 1, top "
             "expressor in Fig. 1, and AF major carrier in Fig. 3. Most "
             "programs do not pass the qperm < 0.10 gate and are reported "
             "as Unresolved (uniform grey tile)."),
            ("d", "Representative programs at the donor × compartment "
             "level."),
            ("e", "Anchor genes underlying compartment-localised program "
             "shifts. Anchor genes are component genes drawn from each "
             "program's leading-edge set and are reported as descriptive "
             "support."),
        ],
        "Figure 3 quantifies the LA arm only; RA-arm-driven rewiring "
        "cannot be resolved in GSE255612. The framework claim is "
        "distributed but structured rewiring rather than a single "
        "compartment-of-origin claim.",
    ),
    (
        "Figure 4 | AF-duration chronicity in the LA appendage",
        [
            ("a", "Schematic. The frozen Fig. 1–3 program panel is "
             "projected to GSE41177 LAA samples and Spearman ρ is computed "
             "on AF-only LAA samples with numeric AF duration in months "
             "(n = 12). Sinus-rhythm samples are visual anchors at month = "
             "0 and are not in the primary statistical model. The duration "
             "covariate represents disease duration in valvular persistent "
             "AF and is not a recurrence covariate."),
            ("b", "LAA program trends with AF duration: heatmap of ρ, "
             "permutation p, and sign concordance under sensitivity "
             "coding A. Sensitivity-supported directionality is "
             "emphasised over a single-cell FDR call."),
            ("c", "Representative programs (linear fit on AF-only samples; "
             "SR anchors are displayed but not fitted): OXPHOS shows the "
             "strongest LAA duration-associated decline."),
            ("d", "LAA-weighted duration trends with region-specific "
             "heterogeneity (LAA primary vs LAPV sensitivity)."),
        ],
        "Power statement: the primary numeric-duration LAA sample is "
        "n = 12. Figure 4 is interpreted as a direction- and "
        "sensitivity-based progression analysis rather than a "
        "high-powered association screen.",
    ),
    (
        "Figure 5 | Cross-figure synthesis and dual-axis cartoon",
        [
            ("a", "Cross-figure triangulation of the 10 frozen chamber "
             "programs. Columns: NCP healthy bias, Fig. 2 bulk rewiring "
             "class, Fig. 3 AF major carrier (uniform grey when "
             "Unresolved), Fig. 4 LAA duration ρ, and triangulation "
             "class. Anchor calls: OXPHOS = Metabolism-led; ECM "
             "proteoglycans = Structure-led."),
            ("b", "Pre-registered multiome support gate. Six frozen "
             "support programs × two expected compartments are tested on "
             "the GSE238242 LAA snRNA pseudobulk (AF n = 3, SR n = 4). "
             "Outcome: the gate is not met and multiome regulator-tier "
             "inference is therefore deferred (see Supplementary "
             "Fig. S5)."),
            ("c", "ROI-level localisation of fibroblast-associated ECM "
             "proteoglycan signal (GSE261363 GeoMx atrial-appendage "
             "segments, n = 24). The spatial layer provides localisation "
             "support only; mechanistic inference is not derived from "
             "GeoMx and relies on cross-layer synthesis."),
            ("d", "Distributed dual-axis mechanism cartoon. Healthy "
             "chamber identity → AF rewiring (matrix arm and metabolic "
             "arm) → distributed substrate. There is no single FDR-robust "
             "dominant carrier; remodelling is distributed but structured "
             "(Fig. 3 Unresolved carriers). Schematic."),
        ],
        None,
    ),
    (
        "Figure 6 | Robustness, plausibility, and synthesis",
        [
            ("a", "Cross-cohort consistency of dual-axis chamber "
             "rewiring. Random-effects meta-analysis on the donor-level "
             "axis composite (equal-weight mean of aCAIrel over the six "
             "frozen support programs per axis). Per-cohort Hedges' g "
             "(AF − SR) and the Hartung–Knapp pooled estimate (REML, "
             "k = 3 cohorts) are reported. With k = 3, the pooled "
             "estimate is interpreted as a cross-cohort consistency "
             "check, not a high-powered pooled discovery."),
            ("b", "Per-program pooled effects for the 10 frozen panel "
             "programs. Structural programs show the most consistent "
             "pooled attenuation; metabolic programs are directionally "
             "consistent but more heterogeneous."),
            ("c", "Framework-discovered chamber-identity REVERSAL "
             "genes — gene-level positive evidence for the distributed "
             "substrate. A reversal gene is defined per gene as a gene "
             "whose chamber direction in healthy donors is flipped in "
             "AF donors of paired-chamber cohorts while sinus-rhythm "
             "donors of the same cohorts preserve the healthy direction "
             "(Methods). Twenty-four reversal genes were identified at "
             "AF-specific reversal in ≥ 2 of 3 paired chamber cohorts; "
             "21 of 24 are RA → LA reversals (RA-biased genes pulled "
             "into the LA in AF), and only 3 are LA → RA reversals. "
             "PRELP (the pre-registered structural-axis anchor of "
             "Fig. 6) is independently recovered by the chamber-"
             "reversal screen. A complementary single-gene cis-eQTL "
             "regime check, including two canonical AF GWAS positive "
             "controls (PITX2, KCNN3), is reported in Supplementary "
             "Fig. S12."),
            ("d", "Cellular carrier of the chamber-identity reversal "
             "signature in single-nucleus AF atrial tissue (GSE255612, "
             "AF n = 18, Ctrl n = 16). Per gene, Δ_expr = "
             "mean(expression in AF) − mean(expression in Ctrl) is "
             "computed in each of the four main compartments and "
             "sign-aligned by the gene's expected AF direction; "
             "positive Δ_aligned indicates expression shift consistent "
             "with the chamber-reversal call. The 24 reversal genes "
             "localise quantitatively to fibroblasts (mean Δ_aligned = "
             "+0.042, ~5-fold higher than any other compartment; 71 % "
             "of the 21 RA → LA matrix-invasion genes up-regulated in "
             "AF in fibroblasts), with smaller contributions from "
             "atrialCM and immune cells. Per-gene heatmap is given in "
             "panel c. The synthesis cartoon previously shown here is "
             "moved to the Graphical Abstract."),
        ],
        None,
    ),
]


SUPP_FIGURES: list[tuple[str, str]] = [
    (
        "Supplementary Figure S1 | Cohort overview",
        "Per-cohort summary of the 71 donors that compose the eight "
        "publicly available datasets analysed in this study. a Donor "
        "counts per cohort, stacked by rhythm (atrial fibrillation in "
        "red, sinus rhythm in blue). b Age distribution per cohort, "
        "split by rhythm (boxes: median and IQR; points: individual "
        "donors). c Sex composition per cohort (male, female, or not "
        "reported). d Disease background composition per cohort. The "
        "underlying donor-level table is released as "
        "results/cohort_demographics_long.tsv and the figure-side "
        "summary as results/figS1_cohort_overview_summary.tsv.",
    ),
    (
        "Supplementary Figure S2 | Donor pairing and sample-attribute "
        "structure",
        "Donor × chamber (GSE128188, GSE79768, GSE115574) or donor × "
        "region (GSE41177 LAA / LAPV) pairing matrix for the four "
        "cohorts that drive Figs 2 and 4. Filled cells indicate that "
        "the corresponding sample is present in the dataset and are "
        "coloured by donor rhythm (AF = red, SR = blue); empty cells "
        "indicate a missing sample. Right-hand side panels report "
        "per-donor sex (green = male, purple = female, grey = not "
        "reported) and age. Companion table: "
        "results/figS2_pairing_summary.tsv.",
    ),
    (
        "Supplementary Figure S3 | Normal Chamber Prior — panel and "
        "gene-level stability",
        "a Program-level direction stability across the 8 donor "
        "resamplings of GSE112339 for each of the 10 panel programs "
        "(stability_consistent_runs, maximum 8). b Mean stability of "
        "the leading-edge genes of each program across the same 8 "
        "resamplings; the right-hand annotation reports the fraction "
        "of leading-edge genes that are significant in all 8 "
        "resamplings. c Gene-level stability distribution split by "
        "direction (LA-high in blue, RA-high in red), supporting the "
        "per-program summaries in panels a–b. Companion deliverable: "
        "results/figS3_ncp_stability_summary.tsv.",
    ),
    (
        "Supplementary Figure S4 | Per-cohort principal-component "
        "overview",
        "Top-variance 2,000-feature, sample-centred PCA on each of the "
        "five bulk cohorts. a GSE112339 (healthy LA / RA reference; "
        "PCA reused from bulk/Project_P01/qc/qc_pca.tsv); paired "
        "donors are connected. b–d Paired-chamber AF-vs-SR cohorts "
        "(GSE128188, GSE79768, GSE115574). e GSE41177 LAA samples "
        "only. In b–e, points are coloured by rhythm (AF = red, SR = "
        "blue) and shaped by chamber (circle = LA, triangle = RA where "
        "applicable). The healthy donor / chamber axis dominates the "
        "global structure in panel a, whereas the rhythm signal is a "
        "smaller-magnitude axis on top of donor variation in the AF "
        "cohorts — motivating the donor-paired aCAIrel framework. "
        "Companion deliverable: results/figS4_per_cohort_pca.tsv.",
    ),
    (
        "Supplementary Figure S5 | Multiome RNA-side support gate (negative)",
        "GSE238242 LAA snRNA pseudobulk Δsingscore (AF − SR) on the six "
        "frozen support programs × four main compartments. The "
        "locked multiome support-gate rule is not met: "
        "no expected (program × compartment) cell reaches moderate or "
        "strong support, and the HS-GAG programs even trend in the "
        "opposite direction (Δscore = −0.03 to −0.04). Multiome "
        "regulator-tier inference is therefore deferred to an independent "
        "cohort.",
    ),
    (
        "Supplementary Figure S6 | Mural sensitivity (pericyte + VSMC)",
        "Δsingscore on the six frozen support programs in the mural "
        "(pericyte + vascular smooth-muscle cell) pseudobulk of GSE238242. "
        "Direction is consistent with the four-compartment main panel and "
        "is not a separate finding.",
    ),
    (
        "Supplementary Figure S7 | Leave-one-cohort-out sensitivity (Fig. 6a)",
        "Pooled Hedges' g for the structural and metabolic axes when each "
        "of the three discovery cohorts is excluded in turn. Direction "
        "stability of the structural axis is the intended check.",
    ),
    (
        "Supplementary Figure S8 | Specificity and falsification supplement",
        "Four-component specificity supplement "
        "(locked falsification protocol; scripts/figS8_falsification.py). "
        "a Size-matched random "
        "gene-set permutations (1,000 per program; ECM proteoglycans "
        "empirical p = 0.10). b Axis-membership permutation of the six "
        "frozen support programs (1,000 partitions; |g_STRUCT| empirical "
        "p = 0.33, consistent with both axes pooling negatively under any "
        "3:3 partition). c Six pre-registered negative-control pathways "
        "(axis-specificity controls: HALLMARK_INFLAMMATORY_RESPONSE, "
        "KRAS_SIGNALING_UP, E2F_TARGETS; biologically unrelated controls: "
        "REACTOME_OLFACTORY_SIGNALING_PATHWAY, KERATINIZATION, "
        "SENSORY_PERCEPTION_OF_TASTE). d Leave-one-program-out (LOPO) "
        "axis robustness — both axes remain direction-preserved (6/6 LOPO "
        "removals), confirming that the dual-axis signal is not driven by "
        "any single component program.",
    ),
    (
        "Supplementary Figure S9 | Reduced-axis readout",
        "Five-gene reduced panel per axis (matrix: PRELP, OMD, FMOD, "
        "CHST3, B4GALT1; metabolic: IDH2, NDUFA8, NDUFS2, COX5A, UQCRC1). "
        "The reduced panel reproduces the full-axis donor-level composite "
        "(Spearman ρ = 0.55–0.90) and the cohort-level SR→AF centroid "
        "direction. The reduced metabolic panel does not replace OXPHOS "
        "for chronicity (ρ_LAA = −0.06 vs full OXPHOS ρ_LAA = −0.58). "
        "This panel is a framework-feasibility check and not a clinical "
        "classifier.",
    ),
    (
        "Supplementary Figure S10 | Methodology audit",
        "Internal-consistency tests (scripts/h7_methodology_audit.R). a NCP "
        "self-recovery (mean aCAIrel ≈ 1 for 10/10 programs in the "
        "healthy bulk reference). b Donor-level direction stability (8/8 "
        "healthy donors with aCAI > 0 for every program). c singscore vs "
        "mean-z cross-method correlation (Spearman ρ = 0.52–0.93 across "
        "programs). d Axis-pooled g unchanged when singscore is replaced "
        "with mean-z (structural g = −0.57 / −0.79; both with I² = 0%).",
    ),
    (
        "Supplementary Figure S11 | Method-validity supplement (P5)",
        "External-style validations (scripts/h8_methodology_p5.R). a "
        "Simulation power (≥ 92% detection at a true effect of 0.3 in "
        "aCAIrel units, cohort-matched donor counts). b Cohort-donor "
        "bootstrap (B = 2,000) of the dual-axis pooled g; intervals are "
        "consistent with and tighter than the Hartung–Knapp pooled "
        "estimates of Fig. 6a. c Classical limma per-chamber AF-vs-SR "
        "comparison; framework-direction concordance = 25/30 program × "
        "cohort cells (83%; 100% in GSE115574); discordant cells were "
        "near-zero classical effects.",
    ),
    (
        "Supplementary Figure S12 | Single-gene cis-eQTL × AF-GWAS "
        "regime check (calibration)",
        "Pre-registered cis-eQTL × AF-GWAS colocalisation on the four "
        "framework anchor genes (PRELP, OMD, IDH2, NDUFA8) and the two "
        "pre-registered canonical AF GWAS positive controls (PITX2 at "
        "the 4q25 locus and KCNN3 at the 1q21 locus). At current GTEx "
        "Atrial Appendage panel power (n = 372), neither the framework "
        "anchors nor the canonical positive controls pass the "
        "conservative cis-eQTL gate (lead p < 5e-8 within ± 1 Mb), and "
        "the colocalisation posterior PP4 is below the supporting "
        "threshold for all six anchors. The failure of the canonical "
        "positive controls indicates that the single-gene mRNA-mediated "
        "regime is calibration-limited at this panel power and is "
        "therefore not informative against a single-gene mediation "
        "argument; the framework-level positive gene-level evidence is "
        "supplied by the chamber-identity reversal genes shown in "
        "Fig. 6c. Anchors are not replaced post hoc. A blood-tissue "
        "sensitivity column (GTEx Whole Blood, eQTL Catalogue "
        "QTD000356) is shown for the framework anchors; the two "
        "positive-control anchors are heart-tissue-specific and have "
        "no whole-blood expression sufficient for cis-eQTL testing.",
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

    # Main figures
    doc.add_heading("Main figures", level=2)
    for title, panels, trailing in MAIN_FIGURES:
        doc.add_heading(title, level=3)
        add_caption_paragraph(doc, panels, trailing_note=trailing)

    # Supplementary figures
    doc.add_heading("Supplementary figures", level=2)
    for title, body in SUPP_FIGURES:
        doc.add_heading(title, level=3)
        add_body(doc, body)

    # Graphical abstract
    doc.add_heading("Graphical abstract", level=2)
    add_body(doc, GRAPHICAL_ABSTRACT)

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
