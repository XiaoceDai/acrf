#!/usr/bin/env python3
"""Build manuscript_drafts/Data_Code_Availability.docx (matches Methods.docx style)."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH

from m0_doc_style import apply_base_style, add_body, add_table

OUT = ROOT / "manuscript_drafts" / "Data_Code_Availability.docx"
OUT.parent.mkdir(parents=True, exist_ok=True)


DATA_TABLE: list[list[str]] = [
    ["Resource", "Accession", "Source", "Role"],
    ["Healthy bulk RNA-seq, paired LA/RA, 8 donors",
     "GSE112339", "NCBI GEO",
     "Fig. 1 — Normal Chamber Prior reference"],
    ["Bulk RNA-seq, paired LA/RA, 5 SR + 5 AF",
     "GSE128188", "NCBI GEO",
     "Fig. 2 — paired chamber-rewiring discovery"],
    ["Affymetrix HG-U133 Plus 2.0, paired LA/RA, 6 SR + 7 AF",
     "GSE79768", "NCBI GEO",
     "Fig. 2 — paired chamber-rewiring discovery"],
    ["Affymetrix HG-U133 Plus 2.0, paired LA/RA, 16 SR + 15 AF "
     "(severe MR background)",
     "GSE115574", "NCBI GEO",
     "Fig. 2 — paired-cohort replication"],
    ["10x snRNA-seq, paired LA/RA, 34 donors",
     "GSE255612", "NCBI GEO",
     "Fig. 3 — LA-arm cell-compartment localisation"],
    ["Affymetrix HG-U133 Plus 2.0, LAA + LAPV, valvular AF + SR donors "
     "with numeric AF duration",
     "GSE41177", "NCBI GEO",
     "Fig. 4 — chronicity / duration analyses"],
    ["10x Multiome (snRNA + snATAC), LAA, 7 donors",
     "GSE238242", "NCBI GEO",
     "Supplementary Fig. S5 — pre-registered multiome support gate"],
    ["NanoString GeoMx DSP, 4 patients, 48 ROIs",
     "GSE261363", "NCBI GEO",
     "Fig. 5c — ROI-level ECM × fibroblast localisation"],
    ["Atrial-fibrillation GWAS (multi-ancestry)",
     "GCST90624411", "EMBL-EBI GWAS Catalog",
     "Fig. 6c — cis-eQTL × AF GWAS coloc"],
    ["cis-eQTL, GTEx v8 Heart Atrial Appendage (n = 372)",
     "QTD000251", "EMBL-EBI eQTL Catalogue",
     "Fig. 6c — primary tissue cis-eQTL"],
    ["cis-eQTL, GTEx v8 Whole Blood (n = 670)",
     "QTD000356", "EMBL-EBI eQTL Catalogue",
     "Fig. 6c — blood-tissue sensitivity"],
]

SOFTWARE_TABLE: list[list[str]] = [
    ["Software", "Version", "Use"],
    ["R", "4.4", "analysis pipeline, meta-analysis, coloc"],
    ["Python", "3.12",
     "figure generation, demographics tabling, schematics"],
    ["singscore (Bioconductor)", "latest",
     "per-sample gene-set scoring"],
    ["metafor", "4.8.0",
     "random-effects meta-analysis (REML / Hartung–Knapp)"],
    ["limma", "latest",
     "classical chamber-stratified differential expression "
     "(Supplementary Fig. S11c)"],
    ["edgeR", "latest", "TMM normalisation, log-CPM"],
    ["coloc", "5.2.3", "coloc.abf colocalisation analysis (Fig. 6c)"],
    ["Rsamtools", "latest",
     "tabix HTTP-range slicing of remote GWAS / eQTL files"],
    ["data.table", "latest", "tabular operations"],
    ["pandas, numpy, scipy, matplotlib", "latest",
     "tabling and figure generation"],
]


def main() -> None:
    doc = Document()
    apply_base_style(doc)

    doc.add_heading("Data and Code availability", level=1)

    # ---- Data availability ----
    doc.add_heading("Data availability", level=2)
    add_body(
        doc,
        "All datasets analysed in this study are publicly available; no "
        "new sequencing or spatial data were generated. The accessions "
        "and analytical roles are listed in the table below. NCBI GEO "
        "accessions resolve to "
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<ID>; GWAS "
        "Catalog accessions to https://www.ebi.ac.uk/gwas/studies/<ID>; "
        "and eQTL Catalogue dataset identifiers to "
        "https://www.ebi.ac.uk/eqtl/api/v2/datasets/<ID>.",
    )
    add_table(doc, DATA_TABLE, col_widths_cm=[6.0, 3.0, 2.8, 4.5])

    # ---- Code availability ----
    doc.add_heading("Code availability", level=2)
    add_body(
        doc,
        "All analysis code, intermediate tables, figure-generation "
        "scripts, and the specification documents that lock the "
        "analytical decisions are available at the project repository "
        "(URL: to be assigned) and will be archived on Zenodo upon "
        "publication (DOI: to be assigned).",
    )
    add_body(
        doc,
        "The repository organises analyses under the locked specification "
        "documents. MANUSCRIPT_FRAMEWORK.md is the single source of truth "
        "for the Chamber-aware Rewiring Framework, the aCAIrel formula, "
        "the rewiring and triangulation classes, and the two "
        "program-defined phenotypes (PRFS, MA-aCM). FIGURE1_WORKFLOW.md "
        "through FIGURE6_WORKFLOW.md are the per-figure pre-registered "
        "specifications (panel maps, colour conventions, forbidden items, "
        "and gates). FIGURE_S8_FALSIFICATION_SPEC.md is the pre-registered "
        "specificity / falsification supplement (negative-control "
        "pathway list and leave-one-program-out rules). "
        "SUPPLEMENTARY_INDEX.md links each main and supplementary figure "
        "to its source script and to its results/*.tsv deliverables. "
        "METHODOLOGY_AUDIT_REPORT.md and METHOD_VALIDITY_P5_REPORT.md "
        "are the internal-consistency and external-style validation "
        "reports that support Supplementary Fig. S10 and Supplementary "
        "Fig. S11, respectively. ILLUSTRATOR_MUSTFIX.md lists the "
        "post-matplotlib editorial polish to be applied in vector-editing "
        "software; that step has no effect on numerical results.",
    )
    add_body(
        doc,
        "Each main or supplementary figure can be reproduced by "
        "executing the corresponding panel script (python3 "
        "scripts/<fig_panel_script>.py) after running its upstream "
        "analysis script (scripts/h*.R, scripts/p*.py, scripts/r*.R, "
        "scripts/d*.R, scripts/f*.R, or scripts/g*.R). The required "
        "input files for each script are listed in its docstring header. "
        "A consolidated figure → script → deliverable mapping is provided "
        "in Supplementary Table (reproducibility_table.tsv).",
    )

    # ---- Software versions ----
    doc.add_heading("Software versions", level=2)
    add_table(doc, SOFTWARE_TABLE, col_widths_cm=[6.0, 2.5, 7.8])
    add_body(
        doc,
        "External tools used: tabix HTTP-range slicing via Rsamtools "
        "(no local index download required); the GTEx Portal API for "
        "gene coordinates; and the GWAS Catalog FTP for harmonised "
        "summary statistics.",
    )

    # ---- Reporting summary ----
    doc.add_heading("Reporting summary", level=2)
    add_body(
        doc,
        "All journal-required reporting (sample sizes, replicates, "
        "statistical tests, software, randomisation, and blinding) is "
        "documented in the Methods section, in the per-figure "
        "specification documents released with the code repository, and "
        "in the consolidated multiple-testing-correction summary "
        "(Supplementary Table mtc_summary.tsv). The Nature Research "
        "Reporting Summary will be completed and attached at submission.",
    )

    # ---- Author contact ----
    doc.add_heading("Author contact", level=2)
    add_body(
        doc,
        "Correspondence and requests for materials should be addressed "
        "to the corresponding author (to be assigned).",
    )

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
