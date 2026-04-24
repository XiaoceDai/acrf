#!/usr/bin/env python3
"""
Build manuscript_drafts/Introduction.docx

Framing: methodology-paper route (Genome Biology / Bioinformatics class).
   - The Chamber-aware Rewiring Framework (CRF) and its donor-level aCAIrel
     metric are the primary deliverables.
   - Atrial fibrillation is the application case.
   - Pre-registration discipline is presented as a distinguishing feature.
   - Honest negatives are framed as falsifiability evidence.

Length target ~650 words / 4 paragraphs (NC Introduction-style; suitable
for Genome Biology / Bioinformatics introductions).
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

OUT = ROOT / "manuscript_drafts" / "Introduction.docx"
OUT.parent.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Introduction paragraphs.
#   P1 : motivating biology — chamber asymmetry of the atria
#   P2 : analytical gap — paired chamber × rhythm rank-deficiency
#   P3 : the methodological deliverable + its falsifiability discipline
#   P4 : "here we show" — concrete application + open release
# ---------------------------------------------------------------------------
PARAGRAPHS: list[str] = [
    # ---- P1 ---------------------------------------------------------------
    "The two atria of the human heart are not symmetric organs. The "
    "left and right atria arise from distinct embryological "
    "progenitors, retain distinct chamber-specifying transcription "
    "factor programs in the adult, and contribute differently to the "
    "substrate of atrial fibrillation (AF), the most common "
    "sustained cardiac arrhythmia [1–3]. Common-variant association "
    "studies have repeatedly highlighted chamber-axis biology as the "
    "dominant heritable contribution to AF risk — most prominently "
    "the 4q25 locus around PITX2, which is itself a transcription "
    "factor of the left-atrial gene program [2, 4]. Despite this "
    "biological backdrop, transcriptomic analyses of AF have "
    "typically pooled left and right atrial samples, or analysed "
    "each chamber separately as if the chamber axis were a nuisance "
    "covariate, rather than treating the chamber axis itself as the "
    "analytical scaffold along which AF rewiring should be "
    "quantified.",

    # ---- P2 ---------------------------------------------------------------
    "There is a concrete statistical reason why the chamber axis is "
    "rarely treated as the scaffold. In paired-chamber cohorts where "
    "every donor contributes both an LA and an RA sample, rhythm "
    "(AF or sinus rhythm) is a donor-level attribute, not a "
    "sample-level attribute: a donor × chamber × rhythm interaction "
    "term in a sample-level linear model is therefore rank-deficient "
    "and cannot identify how AF rewires the chamber axis. As a "
    "consequence, prior chamber-aware AF transcriptomic studies "
    "have most often pooled chambers (losing the chamber axis), "
    "stratified by chamber and reported per-chamber AF-vs-SR "
    "differential expression (losing the cross-chamber comparison), "
    "or relied on aggregate sample-level interactions whose "
    "identifiability conditions are not made explicit. The absence "
    "of a donor-level metric that quantifies how AF preserves, "
    "attenuates, or reverses the *healthy* chamber direction has "
    "left a gap between the chamber-axis biology that motivates AF "
    "research and the analytical pipelines that summarise AF "
    "transcriptomic data.",

    # ---- P3 ---------------------------------------------------------------
    "Here we introduce the Chamber-aware Rewiring Framework (CRF), a "
    "general donor-level pipeline for quantifying transcriptional "
    "rewiring of a chamber-asymmetric tissue under a donor-level "
    "perturbation. CRF first defines a Normal Chamber Prior (NCP) "
    "from a healthy paired-chamber cohort and freezes a panel of "
    "stable LA-vs-RA gene programs. It then quantifies rewiring at "
    "the donor level using an aligned chamber asymmetry index, "
    "aCAIrel, which aligns each program to its healthy chamber "
    "direction and normalises by the healthy bulk magnitude so that "
    "programs of arbitrary leading-edge size share a single, "
    "interpretable scale: aCAIrel ≈ +1 indicates donor-level "
    "preservation, intermediate values quantify attenuation, and "
    "negative values quantify reversal of the healthy direction. "
    "Cohort-level effect sizes (Hedges' g) are pooled across "
    "datasets by REML random-effects meta-analysis with the "
    "Hartung–Knapp variance adjustment [5], and program panels, "
    "axis-membership assignments, single-cell compartment rules, "
    "multiome support gates, and the small set of anchor genes "
    "tested by cis-eQTL × GWAS colocalisation are pre-registered as "
    "text specifications before downstream analyses. Any "
    "pre-registered gate that fails is reported as a finding rather "
    "than replaced post hoc, and the framework's specificity is "
    "checked through a four-component falsification supplement "
    "(size-matched random gene-set permutation, axis-membership "
    "permutation, six pre-registered negative-control pathways, and "
    "a leave-one-program-out test).",

    # ---- P4 ---------------------------------------------------------------
    "We apply CRF to atrial fibrillation across eight publicly "
    "available human atrial transcriptomic and multi-omic datasets "
    "(71 donors; bulk RNA-seq, microarray, single-nucleus RNA-seq, "
    "single-nucleus multiome, and GeoMx spatial; Supplementary "
    "Fig. S1, Supplementary Table 1) [6–9]. We show that AF "
    "consistently attenuates the healthy chamber direction across "
    "three independent paired-chamber cohorts and that the "
    "attenuation is itself distributed across cell compartments "
    "(≥ 7 / 10 panel programs Unresolved at the strict permutation "
    "gate), that mitochondrial OXPHOS attenuation tracks AF "
    "duration in the left atrial appendage, and that a "
    "pre-registered cis-eQTL × AF-GWAS colocalisation on four "
    "canonical anchor genes is uniformly negative in the "
    "tissue-matched GTEx Atrial Appendage panel [10] — a result "
    "consistent with, and not contradictory to, a distributed "
    "dual-axis substrate. All locked specifications, intermediate "
    "deliverables, and figure-generation scripts are released "
    "openly so that independent investigators can reanalyse the "
    "data, change individual gates, or apply CRF to other "
    "chamber-asymmetric tissues without re-implementing the "
    "framework.",
]


def main() -> None:
    doc = Document()
    apply_base_style(doc)

    doc.add_heading("Introduction", level=1)

    for p in PARAGRAPHS:
        para = doc.add_paragraph(p)
        para.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    note = doc.add_paragraph()
    note.paragraph_format.space_before = Pt(18)
    run = note.add_run(
        "End-of-Introduction note (not for typesetter): citation "
        "placeholders [1]–[10] are to be back-filled at submission "
        "against the bibliography file. The framework name 'CRF', the "
        "donor-level metric 'aCAIrel', and the four-component "
        "falsification supplement are introduced here and used "
        "throughout Methods and Results without further definition."
    )
    run.italic = True
    run.font.size = Pt(9)
    run.font.color.rgb = RGBColor(0x66, 0x66, 0x66)

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
