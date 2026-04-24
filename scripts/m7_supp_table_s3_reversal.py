#!/usr/bin/env python3
"""Build manuscript_drafts/Supplementary_Table_S3_chamber_reversal_genes.docx
listing the 24 framework-discovered chamber-identity reversal genes."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import pandas as pd
from docx import Document

from m0_doc_style import apply_base_style, add_body, add_table

OUT = ROOT / "manuscript_drafts" / "Supplementary_Table_S3_chamber_reversal_genes.docx"
OUT.parent.mkdir(parents=True, exist_ok=True)

IN = ROOT / "results" / "n2_chamber_reversal_genes.tsv"


def main() -> None:
    df = pd.read_csv(IN, sep="\t")
    df = df.sort_values(
        ["healthy_direction", "n_cohorts_call", "mean_abs_aCAI_AF"],
        ascending=[True, False, False],
    ).reset_index(drop=True)

    doc = Document()
    apply_base_style(doc)

    doc.add_heading(
        "Supplementary Table S3 | Framework-discovered "
        "chamber-identity reversal genes",
        level=1,
    )

    add_body(
        doc,
        "Twenty-four genes that are AF-specifically reversed in their "
        "chamber identity in at least 2 of the 3 paired chamber AF "
        "cohorts (GSE128188, GSE79768, GSE115574). Per gene, per "
        "cohort, AF-specific reversal requires "
        "Healthy_direction × AF_direction < 0 (sign flip), "
        "Healthy_direction × SR_direction > 0 (SR donors preserve the "
        "healthy direction), |aCAI_healthy| ≥ 0.3 log2 (top ~15 % of "
        "chamber-biased genes in the healthy reference), and "
        "|aCAI_AF| ≥ 0.15 log2 (~80th percentile of |aCAI_AF|). All "
        "aCAI values are reported in log2 (LA − RA) units. Healthy "
        "direction is computed from the 8 paired GSE112339 healthy "
        "donors. PRELP is the framework's pre-registered structural-"
        "axis anchor (Fig. 6) and is independently recovered by this "
        "screen.",
    )

    rows = [[
        "gene", "healthy direction", "n cohorts",
        "aCAI healthy", "mean |aCAI AF|",
        "aCAI AF GSE128188", "aCAI AF GSE79768", "aCAI AF GSE115574",
    ]]
    for _, r in df.iterrows():
        rows.append([
            r["gene"],
            r["healthy_direction"],
            str(int(r["n_cohorts_call"])),
            f"{r['aCAI_healthy']:+.3f}",
            f"{r['mean_abs_aCAI_AF']:.3f}",
            f"{r['aCAI_AF_GSE128188']:+.3f}",
            f"{r['aCAI_AF_GSE79768']:+.3f}",
            f"{r['aCAI_AF_GSE115574']:+.3f}",
        ])

    add_table(doc, rows, col_widths_cm=[2.0, 2.4, 1.2, 1.8, 1.8, 2.2, 2.2, 2.2])

    add_body(
        doc,
        "Direction asymmetry: 21 of 24 reversal genes are RA → LA "
        "(RA-biased in healthy → LA-biased in AF), and only 3 are "
        "LA → RA. The 21 RA → LA genes are dominated by canonical "
        "extracellular-matrix and fibroblast genes (PRELP, OMD, DCN, "
        "FGFBP2, C1S, CTSK, SUSD5, SCARA5) together with stromal-"
        "signalling genes (DLL1, F2R, DNER) — a coherent "
        "matrix-program identity-invasion signature in the LA "
        "chamber under AF.",
    )

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}  ({len(df)} reversal genes)")


if __name__ == "__main__":
    main()
