#!/usr/bin/env python3
"""Build manuscript-facing tables for NCP:"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label

import innovation_ncp_consensus as inc

AXIS_FAM = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
FIG1B = ROOT / "results" / "fig1b_baseline_donor_aCAI.tsv"
FAM_LAB = ROOT / "results" / "ncp_v3_family_labels.tsv"

OUT_PROGRAM = ROOT / "results" / "fig1c_ncp_v3_program_plotting_table.tsv"
OUT_CRIT = ROOT / "results" / "ncp_v3_selection_criteria.tsv"

DISPLAY_MAX_CHARS = 36


def main() -> None:
    axis = pd.read_csv(AXIS_FAM, sep="\t")
    fb = pd.read_csv(FIG1B, sep="\t")[["pathway", "n_genes_used"]]
    labs = pd.read_csv(FAM_LAB, sep="\t")

    merged = axis.merge(fb, on="pathway", how="left")
    merged = merged.merge(labs, on="family_id", how="left")
    miss_lab = merged["family_label"].isna().sum()
    if miss_lab:
        raise ValueError(f"{miss_lab} rows missing family_label — update ncp_v3_family_labels.tsv")

    merged["display_label"] = merged["pathway"].map(
        lambda p: short_pathway_label(p, DISPLAY_MAX_CHARS)
    )

    out_cols = [
        "program_order",
        "pathway",
        "display_label",
        "direction",
        "family_id",
        "family_label",
        "n_genes_used",
        "family_source",
    ]
    merged[out_cols].to_csv(OUT_PROGRAM, sep="\t", index=False)
    print(f"wrote {OUT_PROGRAM.relative_to(ROOT)} ({len(merged)} rows)")

    w = inc.ScoreWeights()
    rows = [
        # --- R discovery (fig1c_dual_method_selection.R) ---
        ("R_discovery", "candidate_gene_sets", "Hallmark + Reactome MSigDB",
         "gmtPathways on h.all.v2026.1.Hs + c2.cp.reactome.v2026.1.Hs", ""),
        ("R_discovery", "fgsea_set_size", "minSize / maxSize",
         "15 / 500", "Gene sets outside range excluded from fgsea / GSVA."),
        ("R_discovery", "candidate_pool_gate", "fgsea padj",
         "< 0.25", "Among significant pathways, ranked by |NES|."),
        ("R_discovery", "candidate_pool_size", "TOP_N_TO_SHOW",
         "80", "Cap on pathways carried into GSVA + LOO."),
        ("R_discovery", "loo_stability_hit", "fgsea padj (per LOO fold)",
         "< 0.05", "Same-direction NES as full-sample run."),
        ("R_discovery", "expression_filter", "edgeR filterByExpr",
         "min.count = 5", "Applied before log-CPM for GSVA."),
        ("prototype_panel_lock", "dual_method_gate", "GSEA padj / GSVA padj / LOO count",
         "< 0.01 / < 0.01 / 8 of 8",
         "Dual-method prototype panel (fig1c_ncp_panel_v2.tsv); benchmark comparator only."),
        ("prototype_panel_lock", "redundancy", "per-direction clustering",
         "Jaccard distance ≥ 0.5",
         "Representative = max |NES| per cluster (see fig1c_dual_method_selection.R)."),
        # --- NCP consensus pruning (innovation_ncp_consensus.py) ---
        ("NCP_consensus_pruning", "input_table", "rows",
         "results/fig1c_pathway_selection_long.tsv",
         "Same 80-row candidate × donor-stability table from R."),
        ("NCP_consensus_pruning", "consensus_gate_strict", "direction_agree",
         "sign(NES) == sign(gsva_logFC)", ""),
        ("NCP_consensus_pruning", "consensus_gate_strict", "stability_consistent_runs",
         f">= {inc.MIN_STAB}", "Minimum donors reproducing GSEA hit under LOO."),
        ("NCP_consensus_pruning", "consensus_gate_strict", "consensus_score",
         f">= {inc.MIN_SCORE}", "Weighted fusion of fgsea/GSVA significance, effect, stability."),
        ("NCP_consensus_pruning", "consensus_gate_strict", "dual significance",
         "padj < 0.05 AND gsva_padj < 0.05", "Strict pool for diversity pruning."),
        ("NCP_consensus_pruning", "consensus_gate_fallback", "consensus_score",
         ">= 0.50", "Used only if strict pool under-fills."),
        ("NCP_consensus_pruning", "consensus_gate_fallback", "dual significance",
         "padj < 0.10 AND gsva_padj < 0.10", ""),
        ("NCP_consensus_pruning", "innovation_tier_core", "consensus / stability / agree",
         ">= 0.78 / 8 / yes", "Tagged core in innovation_ncp_consensus_long.tsv."),
        ("NCP_consensus_pruning", "jaccard_prune_LA", "max Jaccard similarity",
         str(inc.MAX_JACCARD_SIM), "Greedy keep by consensus_score."),
        ("NCP_consensus_pruning", "jaccard_prune_RA", "max Jaccard similarity",
         str(inc.MAX_JACCARD_SIM_RA), "RA arm allowance."),
        ("NCP_consensus_pruning", "jaccard_relaxation", "similarity ladder",
         "0.60, 0.70, 0.80, 0.95",
         "Only if panel still short after greedy selection."),
        ("NCP_consensus_pruning", "panel_balance", "MIN_PER_DIRECTION / TARGET_PER_DIRECTION",
         f"{inc.MIN_PER_DIRECTION} / {inc.TARGET_PER_DIRECTION}", ""),
        ("NCP_consensus_pruning", "panel_size", "target_k",
         "30 (matches dual-method prototype panel when present)",
         "Locks deployed NCP breadth to prototype panel cardinality for benchmarking parity."),
        ("NCP_consensus_pruning", "scoring_norm", "FGSEA_Q_REF / GSVA_Q_REF / NES_REF / LOGFC_REF / STAB_MAX",
         f"{inc.FGSEA_Q_REF} / {inc.GSVA_Q_REF} / {inc.NES_REF} / {inc.LOGFC_REF} / {inc.STAB_MAX}",
         "Clip-to-[0,1] pieces before weighted sum."),
        ("NCP_consensus_pruning", "score_weights", "fgsea_sig / gsva_sig / fgsea_eff / gsva_eff / stability",
         f"{w.fgsea_sig} / {w.gsva_sig} / {w.fgsea_eff} / {w.gsva_eff} / {w.stability}",
         "Sum = 1.0"),
        ("NCP_consensus_pruning", "coherence_penalty", "NES vs GSVA sign mismatch",
         "× 0.35", "Otherwise × 1.0."),
        # --- Axis families ---
        ("Axis_families", "mapping", "script",
         "scripts/v3/innovation_build_v3_axis_families.py",
         "NCP pathway programs inherit fig1e baseline family_id where possible; max(family_id)+1 for new."),
        ("Fig1B_primitive", "baseline_aCAI_gene_overlap",
         "genes intersecting GSE112339 counts",
         ">= 5",
         "Per-program mean log2-CPM(LA)−log2-CPM(RA); see scripts/fig1b_baseline_donor_aCAI.py."),
    ]
    crit = pd.DataFrame(
        rows,
        columns=["stage", "criterion", "parameter", "value", "notes"],
    )
    crit.to_csv(OUT_CRIT, sep="\t", index=False)
    print(f"wrote {OUT_CRIT.relative_to(ROOT)} ({len(crit)} rows)")


if __name__ == "__main__":
    main()
