#!/usr/bin/env python3
"""Build Supplementary Table S5 — full characterisation of the 12
data-driven rewiring-defined gene modules (Fig. 7 source data)."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import pandas as pd
from docx import Document

from m0_doc_style import apply_base_style, add_body, add_table

OUT = ROOT / "manuscript_drafts" / "Supplementary_Table_S5_data_driven_modules.docx"
OUT.parent.mkdir(parents=True, exist_ok=True)
IN = ROOT / "results" / "n1_modules_summary.tsv"


def main() -> None:
    df = pd.read_csv(IN, sep="\t")
    df = df.sort_values(
        ["healthy_direction", "delta_AF_minus_SR"],
        ascending=[True, False],
    ).reset_index(drop=True)

    doc = Document()
    apply_base_style(doc)
    doc.add_heading(
        "Supplementary Table S5 | Data-driven rewiring-defined gene "
        "modules (Fig. 7 source data)",
        level=1,
    )
    add_body(
        doc,
        "Twelve gene modules identified by k-means clustering "
        "(k = 12) on the top-10 PC loadings of the donor-level "
        "chamber asymmetry matrix M_aCAI (3,781 chamber-biased genes "
        "× 59 paired chamber donors across 4 cohorts). Per module "
        "the table reports: number of genes; mean aCAI in healthy "
        "donors and the implied chamber direction; mean aCAI in AF "
        "and SR donors and the rewiring magnitude ΔaCAI(AF − SR); "
        "the maximum Jaccard overlap against any of the 10 NCP "
        "panel programs and the corresponding NCP program identity; "
        "an annotation class (canonical when Jaccard ≥ 0.30, novel "
        "otherwise); and the top eight genes that contribute most to "
        "the module's PC1 + PC2 + PC3 loading.",
    )

    rows = [["module", "n genes", "healthy direction", "mean aCAI healthy",
             "ΔaCAI (AF − SR)", "best NCP program",
             "best NCP Jaccard", "annotation", "top genes"]]
    for _, r in df.iterrows():
        rows.append([
            f"M{int(r['module_id'])}",
            str(int(r["n_genes"])),
            str(r["healthy_direction"]),
            f"{r['healthy_aCAI_mean']:+.3f}",
            f"{r['delta_AF_minus_SR']:+.3f}",
            str(r["best_NCP_program"]) or "—",
            f"{r['best_NCP_jaccard']:.3f}",
            str(r["annotation_class"]),
            str(r["top_genes"]),
        ])
    add_table(doc, rows,
              col_widths_cm=[1.4, 1.2, 2.0, 2.2, 2.4, 4.4,
                             1.6, 1.4, 5.0])

    add_body(
        doc,
        "All 12 modules are framework-discovered novel modules "
        "(maximum Jaccard against any canonical NCP program is below "
        "0.30 in every case). The two main-text-highlighted modules "
        "are M11 (PITX2-anchored mesothelial / LA-identity module, "
        "13 genes; collapses in AF) and M7 (adipose-program module, "
        "27 genes; expands toward LA in AF). The independent recovery "
        "of PITX2 inside the unsupervised M11 cluster validates the "
        "data-driven discovery primitive against the canonical 4q25 "
        "AF GWAS biology.",
    )

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}  ({len(df)} modules)")


if __name__ == "__main__":
    main()
