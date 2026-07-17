#!/usr/bin/env python3
"""Build G-tier genetic traceability table (anchors + families)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PIPELINE = ROOT / "results" / "figure6_mr_anchors_v3_pipeline_wide.tsv"
OVERLAP = ROOT / "results" / "fig6c_family_af_overlap.tsv"
LABELS = ROOT / "results" / "ncp_v3_family_labels.tsv"
FAMILY_MAP = {
    "PRELP": 10, "OMD": 10, "IDH2": 7, "NDUFA8": 8,
    "BMP10": 11, "SMAD9": 11, "SMAD6": 11,
    "NTRK2": 15, "NTF3": 15, "FRS2": 15,
}

OUT_TABLE = ROOT / "results" / "figure6_g_tier_traceability.tsv"
OUT_RULES = ROOT / "results" / "figure6_g_tier_rules.tsv"

G_TIER_RULES = [
    ("G0", "No genome-wide AF signal in ±1 Mb (gwas_min_p ≥ 5e-8); family overlap FDR > 0.10 for all external sets"),
    ("G1", "Near-GW AF locus (gwas_min_p < 5e-8) but no strong family overlap (min FDR > 0.10)"),
    ("G2", "Shared-locus / H3-dominant (PP3 high, PP4 < 0.70) OR Bonferroni eGene without PP4 support; OR family overlap FDR ≤ 0.10"),
    ("G3", "G2 anchor criteria + family overlap FDR ≤ 0.05 for anchor's family; OR family-only strong enrichment FDR ≤ 0.05"),
    ("G4", "Same-variant colocalisation PP4 ≥ 0.70 (expected 0/12 at current atrial eQTL power)"),
]


def _anchor_tier(row: pd.Series, fam_min_fdr: float) -> str:
    gwas_p = float(row["gwas_min_p"]) if pd.notna(row["gwas_min_p"]) else 1.0
    pp3 = float(row["coloc_PP3"]) if pd.notna(row["coloc_PP3"]) else 0.0
    pp4 = float(row["coloc_PP4"]) if pd.notna(row["coloc_PP4"]) else 0.0
    is_egene = bool(row.get("is_eGene_atrial_bonferroni", row.get("bonf_passed", False)))
    if pp4 >= 0.70:
        return "G4"
    near_gw = gwas_p < 5e-8
    h3 = pp3 >= 0.30 and pp4 < 0.70
    if near_gw and (h3 or is_egene or fam_min_fdr <= 0.05):
        return "G3" if fam_min_fdr <= 0.05 else "G2"
    if near_gw:
        return "G1"
    if h3 or is_egene or fam_min_fdr <= 0.10:
        return "G2"
    return "G0"


def _family_tier(min_fdr: float, best_overlap: int) -> str:
    if min_fdr <= 0.05 and best_overlap >= 1:
        return "G3"
    if min_fdr <= 0.10 and best_overlap >= 1:
        return "G2"
    if best_overlap >= 1:
        return "G1"
    return "G0"


def main() -> None:
    pw = pd.read_csv(PIPELINE, sep="\t")
    ov = pd.read_csv(OVERLAP, sep="\t")
    labels = pd.read_csv(LABELS, sep="\t")

    fam_fdr = ov.groupby("family_id")["fdr_bh"].min()
    fam_best_n = ov.groupby("family_id")["n_overlap"].max()

    anchor_rows = []
    for _, r in pw.iterrows():
        g = r["anchor"]
        fid = FAMILY_MAP.get(g)
        fdr = float(fam_fdr.get(fid, 1.0)) if fid else 1.0
        tier = _anchor_tier(r, fdr)
        gwas_p = r.get("gwas_min_p", np.nan)
        near_gw = pd.notna(gwas_p) and float(gwas_p) < 5e-8
        sqtl_null = str(r.get("verdict_sqtl", "no_eligible_sQTL"))
        anchor_rows.append({
            "entity_type": "anchor",
            "entity_id": g,
            "family_id": fid if fid else "",
            "G_tier": tier,
            "near_GW_AF_locus": near_gw,
            "coloc_PP3": r.get("coloc_PP3"),
            "coloc_PP4": r.get("coloc_PP4"),
            "is_eGene_bonferroni": r.get("bonf_passed", False),
            "family_min_overlap_fdr": fdr,
            "sqtl_verdict": sqtl_null,
            "verdict": (
                f"{tier}: calibration-limited coloc" if tier in ("G1", "G2", "G3") and float(r.get("coloc_PP4", 0) or 0) < 0.70
                else f"{tier}: traceability tier"
            ),
        })

    family_rows = []
    for _, lab in labels.sort_values("family_id").iterrows():
        fid = int(lab["family_id"])
        min_fdr = float(fam_fdr.get(fid, 1.0))
        best_n = int(fam_best_n.get(fid, 0))
        tier = _family_tier(min_fdr, best_n)
        family_rows.append({
            "entity_type": "family",
            "entity_id": lab["short_label"],
            "family_id": fid,
            "G_tier": tier,
            "near_GW_AF_locus": "",
            "coloc_PP3": "",
            "coloc_PP4": "",
            "is_eGene_bonferroni": "",
            "family_min_overlap_fdr": min_fdr,
            "sqtl_verdict": "n/a (family-level)",
            "verdict": f"{tier}: family×AF-set overlap (best n={best_n}, min FDR={min_fdr:.3g})",
        })

    out = pd.DataFrame(anchor_rows + family_rows)
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TABLE, sep="\t", index=False)

    rules = pd.DataFrame(G_TIER_RULES, columns=["G_tier", "rule"])
    rules.to_csv(OUT_RULES, sep="\t", index=False)
    print(f"wrote {OUT_TABLE.relative_to(ROOT)} ({len(out)} rows)")
    print(f"wrote {OUT_RULES.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
