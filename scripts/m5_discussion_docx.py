#!/usr/bin/env python3
"""
Build manuscript_drafts/Discussion.docx (NC style; ~1,200 words; 6 paragraphs).

Tone : framework-first / conservative / no clinical-target claim.
References are placeholders [N] to be back-filled at submission.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt, RGBColor

from m0_doc_style import apply_base_style

OUT = ROOT / "manuscript_drafts" / "Discussion.docx"
OUT.parent.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Discussion text. NC style: no internal sub-headings; paragraph breaks only.
# ---------------------------------------------------------------------------
PARAGRAPHS: list[str] = [
    # ---- 1. Opening / restated main finding ---------------------------------
    "Across eight publicly available human atrial transcriptomic and "
    "multi-omic resources, a single chamber-aware framework anchored "
    "on a healthy LA-vs-RA program panel resolved atrial fibrillation "
    "(AF) into a distributed, dual-axis substrate rather than a "
    "master-gene or master-compartment phenotype. The matrix arm "
    "(ECM proteoglycans and heparan-sulfate / glycosaminoglycan "
    "metabolism) and the metabolic arm (oxidative phosphorylation, "
    "aerobic respiration and electron transport, and cytoplasmic "
    "translation) were each identifiable as cross-cohort consistent "
    "attenuations of the healthy chamber direction; neither arm was "
    "carried by a single FDR-robust dominant cell carrier in the "
    "left-atrial single-nucleus arm; and neither arm was mediated "
    "through a pre-registered single-gene cis-eQTL anchor in the "
    "tissue-matched GTEx Atrial Appendage panel. Read together, these "
    "three observations argue that AF substrate phenotyping is best "
    "framed at the level of cross-compartment transcriptional "
    "programs rather than at the level of any single causal gene or "
    "any single carrier cell type.",

    # ---- 2. Methodological contribution -------------------------------------
    "A central methodological observation is that paired chamber × "
    "rhythm cohorts are intrinsically rank-deficient at the "
    "sample-level interaction term: rhythm is a donor-level "
    "attribute, so a sample-level donor × chamber × rhythm model "
    "cannot identify how AF rewires the chamber axis. The "
    "donor-level aligned chamber asymmetry index, aCAIrel, was "
    "designed to recover that interaction by aligning each program "
    "to its healthy chamber direction and then normalising by the "
    "healthy bulk magnitude. The normalisation step removes the "
    "dependence of singscore magnitudes on leading-edge set size, "
    "which is what allows programs of very different sizes to share "
    "a single, interpretable scale on which preservation, "
    "attenuation, and reversal of the healthy chamber direction can "
    "be read off directly. Internal-consistency checks "
    "(Supplementary Fig. S10) and external-style validations "
    "(Supplementary Fig. S11) confirmed that this metric recovers "
    "known direction in healthy bulk, agrees with a classical "
    "limma chamber-stratified analysis in 25 / 30 program × cohort "
    "cells, and matches its REML / Hartung–Knapp meta-analytic "
    "intervals against an independent cohort-donor bootstrap. "
    "We anticipate that the same chamber-aware donor-level recipe "
    "will be useful in any setting where a tissue is structurally "
    "asymmetric and the perturbation is donor-level rather than "
    "sample-level. A second design feature of the framework is "
    "that, on the same donor-level chamber asymmetry matrix M_aCAI "
    "established at framework construction (Fig. 1), it provides "
    "two complementary discovery primitives: a curated primitive — "
    "signed GSEA over MSigDB / Reactome on the chamber-biased gene "
    "cohort, yielding the 10 canonical NCP programs that constitute "
    "the locked calibration set used downstream — and a data-driven "
    "primitive — PCA of M_aCAI followed by k-means on the top-PC "
    "loadings (Fig. 1d–e), yielding twelve rewiring-defined modules "
    "none of which exceeded a 30 % Jaccard overlap with any "
    "canonical NCP program. Among the twelve data-driven modules a "
    "tight PITX2-anchored mesothelial module (M11) collapses in AF "
    "and a 27-gene adipose-program module (M7; ADIPOQ, LEP, PLIN1, "
    "PCK1) expands toward the LA in AF, both observations that the "
    "GSEA-based primitive cannot reach. The independent recovery of "
    "PITX2 inside an unsupervised cluster validates the data-driven "
    "primitive against canonical AF biology and demonstrates that "
    "the two primitives are mutually consistent on the locus that "
    "matters most; the LA expansion of the adipose module is a "
    "framework-specific finding that motivates further investigation "
    "of epicardial-adipose-tissue involvement in AF substrate. Both "
    "data-driven module predictions are independently validated at "
    "single-nucleus resolution in Fig. 1f: the M7 adipose module "
    "peaks in adipocytes (row z = +3.45) and the M11 PITX2 / "
    "mesothelial module peaks in epicardial cells (row z = +3.46), "
    "in both cases on cell types that had been excluded from the "
    "original 4-compartment analysis lock and that were added back "
    "specifically because the data-driven primitive identified them "
    "as biologically central.",

    # ---- 3. Relation to the prior literature --------------------------------
    "The framework's findings sit naturally alongside three strands "
    "of the existing AF literature. First, multi-ancestry AF GWAS "
    "have identified hundreds of common-variant signals dominated by "
    "the 4q25 / PITX2 locus and a small number of "
    "transcription-factor and ion-channel loci [1, 2]; targeted "
    "single-gene mediation analyses have nevertheless struggled to "
    "produce a generalisable, tissue-matched cis-eQTL × AF-GWAS "
    "colocalisation outside that locus [3]. In this study, even the "
    "two pre-registered canonical AF GWAS positive controls (PITX2 "
    "and KCNN3) failed the conservative single-gene cis-eQTL gate at "
    "current GTEx Atrial Appendage panel power (Supplementary "
    "Fig. S12), confirming that the classical single-gene mRNA-"
    "mediated regime is presently calibration-limited rather than "
    "diagnostically informative for distributed substrate models. "
    "Against that classical-regime null, our framework provides a "
    "gene-level positive signature uniquely visible to CRF: a class "
    "of 24 chamber-identity REVERSAL genes (Fig. 6c) that are "
    "reproducibly RA → LA reversed across at least two of the three "
    "paired-chamber AF cohorts. The 21:3 RA → LA asymmetry of the "
    "reversal class, dominated by extracellular-matrix proteoglycans "
    "and fibroblast genes (PRELP, OMD, DCN, FGFBP2, C1S, CTSK), "
    "describes a previously unreported pattern of chamber-identity "
    "invasion in AF — RA-biased matrix programs are pulled into the "
    "LA — and the independent recovery of PRELP through both the "
    "framework's pre-registered structural anchor list and the "
    "unbiased chamber-reversal screen provides cross-method "
    "convergence on the same gene. The cellular carrier of the "
    "reversal signature, identified by projecting the 24 reversal "
    "genes onto the GSE255612 single-nucleus compartments (Fig. 6d), "
    "is the fibroblast (mean Δ_aligned ≈ 5-fold higher than any other "
    "compartment), with smaller contributions from atrialCM and "
    "immune cells; this gene-level carrier resolution is itself a "
    "framework-specific finding that the prior literature has not "
    "reported. Second, single-nucleus surveys of human atrial tissue "
    "have reported compartment-level transcriptional differences in "
    "AF without resolving a single dominant cellular carrier at the "
    "program level [4–7]; the program-level Fig. 3 finding "
    "(unresolved at the strict permutation gate) recapitulates that "
    "lack of program-level carrier resolution, while the gene-level "
    "Fig. 6d analysis identifies the chamber-reversal carrier "
    "precisely. Third, atrial fibrosis and mitochondrial "
    "dysfunction have been described separately in AF [8–11]; the "
    "dual-axis synthesis in Figs 5–6 places them within the same "
    "chamber-aware framework and shows that they pool negatively as "
    "cross-cohort consistent axes — most consistently for the "
    "structural arm, more heterogeneously for the metabolic arm.",

    # ---- 4. Translational implications (conservative) -----------------------
    "The translational implication that we feel can be defended at "
    "this stage is therefore framework-level rather than "
    "target-level. A donor-level aCAIrel readout, or the reduced "
    "five-gene-per-axis panel that reproduces the full-axis "
    "donor composite at Spearman ρ = 0.55–0.90 (Supplementary "
    "Fig. S9), provides a candidate molecular phenotyping tool that "
    "operates beyond a binary AF / SR state and could in principle "
    "be applied to bulk transcriptomic data from clinical atrial "
    "biopsy material. The OXPHOS attenuation that tracks AF duration "
    "in the LA appendage (Fig. 4) is consistent with a metabolic "
    "progression component that warrants prospective testing in "
    "cohorts with paired clinical follow-up, and the distributed "
    "nature of the substrate suggests that combination strategies "
    "that engage both the matrix and metabolic arms are conceptually "
    "more defensible than single-target approaches. We deliberately "
    "stop short of nominating actionable molecular targets, because "
    "the evidence assembled here does not include independent "
    "regulator-tier or causal-perturbation data; the framework is "
    "intended as a substrate-phenotyping scaffold, not as a target "
    "discovery study.",

    # ---- 5. Limitations / honest negatives ----------------------------------
    "Several limitations frame the strength of the conclusions. The "
    "primary chronicity panel (Fig. 4) used n = 12 AF donors with "
    "numeric AF duration in months and is interpreted as a direction- "
    "and sensitivity-supported progression analysis rather than a "
    "high-powered association screen; the duration covariate is a "
    "disease-duration proxy in valvular persistent AF, not a "
    "recurrence covariate. The pre-registered multiome regulator-tier "
    "support gate (Fig. 5b) was not met on the GSE238242 LAA snRNA "
    "pseudobulk (AF n = 3, SR n = 4), so all chromatin- and "
    "regulator-tier inference was deferred to an independent cohort; "
    "the field will need a multiome cohort that satisfies the "
    "support gate before the dual-axis substrate can be extended to "
    "transcription-factor-level claims. The cis-eQTL × AF-GWAS coloc "
    "(Fig. 6c) used four pre-registered anchor genes and was not "
    "replaced post hoc; broader anchor sets, larger or "
    "atrial-tissue-specific cis-eQTL resources [12], and rare-variant "
    "association would together provide the more sensitive single-gene "
    "test that this study deliberately did not perform. Sex-stratified "
    "meta-analysis was constrained by data availability — only "
    "GSE79768 contained both sexes with sufficient AF and SR "
    "representation — and is therefore not presented as a main-text "
    "result. Spatial localisation (Fig. 5c) is at the GeoMx ROI level "
    "and would benefit from single-cell-resolution spatial "
    "confirmation in independent tissue. Finally, all cohorts "
    "analysed here were drawn from public archives; recruitment and "
    "selection biases of those original studies propagate into the "
    "framework's effect sizes.",

    # ---- 6. Conclusion / outlook --------------------------------------------
    "In sum, a chamber-aware framework based on a locked panel of ten "
    "stable LA-vs-RA programs and a donor-level aCAIrel metric "
    "yields a reproducible, distributed, dual-axis substrate model "
    "of AF that is supported by cross-cohort consistency and is not "
    "supported by a single-gene mediation argument at the resolution "
    "tested. The framework, the locked specifications, and the "
    "intermediate deliverables that ground every numerical claim in "
    "this paper are released openly to permit independent replication "
    "and extension. The most informative follow-up experiments would "
    "be (i) an independent multiome cohort sized to satisfy the "
    "pre-registered regulator-tier support gate, (ii) "
    "single-cell-resolution spatial profiling of AF and SR atrial "
    "tissue to test whether the LA-arm Unresolved carrier signature "
    "becomes resolved at sufficient compartment depth, and (iii) "
    "prospective phenotyping cohorts to test whether the dual-axis "
    "donor-level readout adds clinically actionable information "
    "beyond the binary AF / SR distinction.",
]


def main() -> None:
    doc = Document()
    apply_base_style(doc)

    doc.add_heading("Discussion", level=1)

    for p in PARAGRAPHS:
        para = doc.add_paragraph(p)
        para.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    note = doc.add_paragraph()
    note.paragraph_format.space_before = Pt(18)
    run = note.add_run(
        "End-of-Discussion note (not for typesetter): citation "
        "placeholders [1]–[12] in the Discussion text are to be "
        "back-filled at submission against the bibliography file. "
        "All numerical claims trace to the deliverables released "
        "alongside the manuscript (see Code availability)."
    )
    run.italic = True
    run.font.size = Pt(9)
    run.font.color.rgb = RGBColor(0x66, 0x66, 0x66)

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
