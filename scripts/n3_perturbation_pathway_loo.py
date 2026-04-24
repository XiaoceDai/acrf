#!/usr/bin/env python3
"""
N3  Layer 1 perturbation analysis  ·  Pathway leave-one-out.

For each chamber-reversal gene g  (n = 24)  and each NCP_v2 program p
(n = 30):

    score_with    = singscore(LE_p)   on each (donor, cell type)
    score_without = singscore(LE_p − {g})  on the same samples
    Δsing_g,p     = score_with − score_without

    keystone score per (gene, program) =
        |mean over (AF donors × LA fibroblasts)| of Δsing
        (relative to the program's overall  Δsing variance)

This quantifies how much each individual reversal gene contributes to the
program's signal in LA fibroblasts (the dominant Fig 6C carrier);  large
|keystone| genes are candidate hubs whose silencing collapses program-
level scoring.

Inputs:
    results/n2_chamber_reversal_genes.tsv          (24 reversal genes)
    results/fig3_primary_program_le.tsv            (NCP primary LE genes)
    Project_P3/GSE255612/pseudobulk_full/<ct>/{expr.tsv, samples.tsv}

Outputs:
    results/n3_perturbation_loo_per_celltype.tsv   long: gene × program × ct
    results/n3_perturbation_loo_summary.tsv         per (gene, program) summary
    results/n3_keystone_genes.tsv                   per gene keystone ranking
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PB_BASE     = ROOT / "Project_P3" / "GSE255612" / "pseudobulk_full"
REV_TSV     = ROOT / "results" / "n2_chamber_reversal_genes.tsv"
LE_TSV      = ROOT / "results" / "fig3_primary_program_le.tsv"

OUT_LONG    = ROOT / "results" / "n3_perturbation_loo_per_celltype.tsv"
OUT_SUMMARY = ROOT / "results" / "n3_perturbation_loo_summary.tsv"
OUT_KEY     = ROOT / "results" / "n3_keystone_genes.tsv"

# 11 cell types matching Fig 3 (B_cell + epicardial excluded)
COMPARTMENTS = [
    "atrialCM", "fibroblast", "endothelial",
    "lymphatic_endothelial",
    "pericyte", "vascular_smooth_muscle",
    "adipocyte",
    "macrophage", "lymphocyte", "mast_cell",
    "neuronal",
]


def singscore_centered(expr_sub: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Centred-singscore equivalent on a (genes × samples) sub-matrix.
    For each sample: mean of within-sample percentile-rank for the LE
    genes minus 0.5 (so 0 = symmetric, +0.5 = top of every gene).
    """
    if expr_sub.empty or not genes:
        return pd.Series(dtype=float)
    present = [g for g in genes if g in expr_sub.index]
    if len(present) < 3:
        return pd.Series(np.nan,
                            index=expr_sub.columns, dtype=float)
    # within-sample rank (descending → high expr = high rank)
    ranks = expr_sub.rank(axis=0, pct=True, ascending=True, method="average")
    score = ranks.loc[present].mean(axis=0) - 0.5
    return score


def main() -> None:
    rev = pd.read_csv(REV_TSV, sep="\t")
    rev_genes = rev["gene"].astype(str).tolist()
    print(f"loaded reversal panel: {len(rev_genes)} genes")

    le = pd.read_csv(LE_TSV, sep="\t")
    print(f"loaded NCP_v2 LE table: {len(le)} programs")

    long_rows = []

    for compartment in COMPARTMENTS:
        gdir = PB_BASE / compartment
        if not (gdir / "expr.tsv").exists():
            print(f"  SKIP {compartment}: no pseudobulk")
            continue
        expr = pd.read_csv(gdir / "expr.tsv", sep="\t", index_col=0)
        meta = pd.read_csv(gdir / "samples.tsv", sep="\t")
        meta_idx = meta.set_index("donor_id")
        # restrict expr to donors with metadata
        common = [d for d in expr.columns if d in meta_idx.index]
        expr = expr[common]
        meta_idx = meta_idx.loc[common]
        print(f"  {compartment}: {expr.shape[0]} genes × {len(common)} samples")

        # per-program loop
        for _, p in le.iterrows():
            program_id    = p["program_id"]
            program_order = int(p["program_order"])
            family_id     = int(p["family_id"]) if not pd.isna(p["family_id"]) else -1
            le_genes      = [g.strip() for g in str(p["leading_edge_genes"]).split(",")]
            le_genes      = [g for g in le_genes if g]

            # baseline score (with all LE genes)
            base_score = singscore_centered(expr, le_genes)
            if base_score.dropna().empty:
                continue

            # split AF / Ctrl per donor
            rhythm = meta_idx.loc[base_score.index, "rhythm"].values
            af_mask = (rhythm == "AF")
            ct_mask = (rhythm == "Ctrl")

            base_af   = base_score[af_mask].mean()
            base_ctrl = base_score[ct_mask].mean()
            base_delta = base_af - base_ctrl   # baseline Δ singscore

            # leave-one-out for each reversal gene that's in this program
            for g in rev_genes:
                if g not in le_genes:
                    # gene not in this program's LE → no leave-one-out signal
                    continue
                le_minus = [x for x in le_genes if x != g]
                pert_score = singscore_centered(expr, le_minus)
                if pert_score.dropna().empty:
                    continue
                pert_af   = pert_score[af_mask].mean()
                pert_ctrl = pert_score[ct_mask].mean()
                pert_delta = pert_af - pert_ctrl

                # impact of removing gene g: how much does program-level
                # Δ change?  positive = gene g contributes to AF rewiring
                impact_delta = base_delta - pert_delta
                # also AF-only baseline shift (drop in AF score when g
                # silenced) - more interpretable for "if we silence g in AF":
                impact_af_score = base_af - pert_af

                long_rows.append({
                    "gene":           g,
                    "program_id":     program_id,
                    "program_order":  program_order,
                    "family_id":      family_id,
                    "compartment":    compartment,
                    "n_AF":           int(af_mask.sum()),
                    "n_Ctrl":         int(ct_mask.sum()),
                    "n_LE_total":     len(le_genes),
                    "base_score_AF":  round(float(base_af), 5),
                    "base_score_Ctrl": round(float(base_ctrl), 5),
                    "base_delta":     round(float(base_delta), 5),
                    "pert_score_AF":  round(float(pert_af), 5),
                    "pert_score_Ctrl": round(float(pert_ctrl), 5),
                    "pert_delta":     round(float(pert_delta), 5),
                    "impact_delta":   round(float(impact_delta), 5),
                    "impact_af_score": round(float(impact_af_score), 5),
                })

    long = pd.DataFrame(long_rows)
    OUT_LONG.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(OUT_LONG, sep="\t", index=False)
    print(f"\nwrote {OUT_LONG.relative_to(ROOT)} ({len(long)} rows)")

    # per (gene, program) summary across the 11 cell types
    summary = (long.groupby(["gene", "program_id", "program_order",
                                  "family_id"])
                    .agg(
                        n_celltypes_present=("compartment", "nunique"),
                        mean_impact_delta=("impact_delta", "mean"),
                        max_abs_impact_delta=("impact_delta",
                                                  lambda s: s.abs().max()),
                        argmax_celltype=("compartment",
                                            lambda s: s.iloc[s.reset_index(drop=True)
                                                                .index[long.loc[s.index, "impact_delta"].abs().argmax()]]),
                    )
                    .reset_index())
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)} ({len(summary)} rows)")

    # per gene keystone ranking - LA fibroblast is the dominant carrier so
    # keystone score = mean |impact_delta| over LA fibroblast across all the
    # programs the gene appears in
    fib_long = long[long["compartment"] == "fibroblast"]
    keystone = (fib_long.groupby("gene")
                          .agg(n_programs_member=("program_id", "nunique"),
                                 max_abs_impact_in_fib=("impact_delta",
                                                          lambda s: s.abs().max()),
                                 mean_abs_impact_in_fib=("impact_delta",
                                                           lambda s: s.abs().mean()),
                                 sum_impact_in_fib=("impact_delta", "sum"))
                          .reset_index()
                          .sort_values("max_abs_impact_in_fib",
                                          ascending=False))

    # also overall (any cell type) max |impact_delta|
    overall = (long.groupby("gene")
                    .agg(any_max_abs_impact=("impact_delta",
                                                  lambda s: s.abs().max()),
                          any_argmax_program=("program_id",
                                                  lambda s: long.loc[s.index]
                                                               .iloc[long.loc[s.index, "impact_delta"]
                                                                          .abs().argmax()]
                                                               ["program_id"]),
                          any_argmax_celltype=("compartment",
                                                  lambda s: long.loc[s.index]
                                                               .iloc[long.loc[s.index, "impact_delta"]
                                                                          .abs().argmax()]
                                                               ["compartment"]))
                    .reset_index())
    keystone = keystone.merge(overall, on="gene", how="left")

    keystone.to_csv(OUT_KEY, sep="\t", index=False)
    print(f"wrote {OUT_KEY.relative_to(ROOT)} ({len(keystone)} rows)")

    # console preview
    print("\n=== TOP 10 keystone reversal genes (by |impact| in LA fibroblast) ===")
    print(keystone.head(10)[["gene", "n_programs_member",
                                  "max_abs_impact_in_fib",
                                  "any_argmax_program",
                                  "any_argmax_celltype"]].to_string(index=False))

    print("\n=== TOP 10 (gene, program) pairs by |impact_delta| in fibroblast ===")
    fib_pairs = fib_long.assign(abs_impact=fib_long["impact_delta"].abs()) \
                          .sort_values("abs_impact", ascending=False)
    print(fib_pairs.head(10)[["gene", "program_id", "n_LE_total",
                                   "base_delta", "pert_delta",
                                   "impact_delta"]].to_string(index=False))


if __name__ == "__main__":
    main()
