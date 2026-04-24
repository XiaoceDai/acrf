#!/usr/bin/env python3
"""
N2b v2 — Map the 24 chamber-identity reversal genes onto the
13-compartment GSE255612 LA snRNA pseudobulk (Project_P3/GSE255612/
pseudobulk_full).  Mirrors Fig 3's 11 informative cell types.

Per gene, per compartment:
    delta_expr     = mean(expr in AF biosamples) - mean(expr in Ctrl)
    delta_aligned  = delta_expr * sign(expected_AF_dir)

For RA-biased reversal genes (RA→LA, n = 21), positive delta_aligned
in a compartment means the gene is up-regulated in AF in that
compartment, consistent with matrix-invasion direction.  For LA-biased
reversal genes (LA→RA, n = 3), positive delta_aligned means the gene
is down-regulated in AF.  delta_aligned > 0  ⇒  carrier candidate.

Outputs:
    results/n2b_primary_reversal_per_celltype_long.tsv     (gene × ct rows)
    results/n2b_primary_reversal_per_celltype_summary.tsv  (ct rows)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

REVERSAL_TSV = ROOT / "results" / "n2_chamber_reversal_genes.tsv"
PB_BASE      = ROOT / "Project_P3" / "GSE255612" / "pseudobulk_full"

# 11 cell types passing donor gates in Fig 3 (B_cell + epicardial excluded
# for n_Ctrl < 4); we still score them all here and mark which pass.
COMPARTMENTS_ALL = [
    "atrialCM", "fibroblast", "endothelial",
    "lymphatic_endothelial",
    "pericyte", "vascular_smooth_muscle",
    "adipocyte",
    "macrophage", "lymphocyte", "mast_cell",
    "neuronal",
    # excluded from Fig 3 main analysis but keep here for transparency
    "B_cell", "epicardial",
]
MIN_DONORS = 4   # match R4b convention

OUT_LONG    = ROOT / "results" / "n2b_primary_reversal_per_celltype_long.tsv"
OUT_SUMMARY = ROOT / "results" / "n2b_primary_reversal_per_celltype_summary.tsv"


def main() -> None:
    rev = pd.read_csv(REVERSAL_TSV, sep="\t")
    rev["expected_dir_in_AF"] = np.where(
        rev["healthy_direction"] == "RA-biased", +1, -1)
    print(f"Reversal genes: {len(rev)}  "
          f"(RA→LA n={(rev['healthy_direction']=='RA-biased').sum()}, "
          f"LA→RA n={(rev['healthy_direction']=='LA-biased').sum()})")

    long_rows = []
    summary_rows = []
    for compartment in COMPARTMENTS_ALL:
        gdir = PB_BASE / compartment
        if not (gdir / "expr.tsv").exists():
            print(f"  skip {compartment}: no pseudobulk")
            continue
        expr = pd.read_csv(gdir / "expr.tsv", sep="\t", index_col=0)
        meta = pd.read_csv(gdir / "samples.tsv", sep="\t")

        af_donors = meta.loc[meta["rhythm"] == "AF",   "donor_id"].tolist()
        ct_donors = meta.loc[meta["rhythm"] == "Ctrl", "donor_id"].tolist()
        af_present = [s for s in af_donors if s in expr.columns]
        ct_present = [s for s in ct_donors if s in expr.columns]
        included = (len(af_present) >= MIN_DONORS) and \
                     (len(ct_present) >= MIN_DONORS)
        print(f"  {compartment:24s}  AF={len(af_present)}  Ctrl={len(ct_present)}  "
                f"{'INCLUDED' if included else 'excluded (<4)'}")

        for _, r in rev.iterrows():
            g = str(r["gene"])
            row = {
                "gene": g, "compartment": compartment,
                "compartment_included": included,
                "healthy_direction": r["healthy_direction"],
                "expected_dir_in_AF": int(r["expected_dir_in_AF"]),
                "n_AF": len(af_present), "n_Ctrl": len(ct_present),
            }
            if g not in expr.index or not included:
                row.update(mean_expr_AF=np.nan, mean_expr_Ctrl=np.nan,
                              delta_expr=np.nan, delta_aligned=np.nan,
                              fold_AF_vs_Ctrl=np.nan)
                long_rows.append(row); continue
            af_v = expr.loc[g, af_present].dropna()
            ct_v = expr.loc[g, ct_present].dropna()
            if len(af_v) < 3 or len(ct_v) < 3:
                continue
            af_mean = float(af_v.mean()); ct_mean = float(ct_v.mean())
            d = af_mean - ct_mean
            d_aligned = d * r["expected_dir_in_AF"]
            row.update(
                mean_expr_AF=round(af_mean, 5),
                mean_expr_Ctrl=round(ct_mean, 5),
                delta_expr=round(d, 5),
                delta_aligned=round(d_aligned, 5),
                fold_AF_vs_Ctrl=round((af_mean / ct_mean)
                                          if ct_mean > 0 else np.nan, 4))
            long_rows.append(row)

    long = pd.DataFrame(long_rows)
    long.to_csv(OUT_LONG, sep="\t", index=False)
    print(f"wrote {OUT_LONG.relative_to(ROOT)} ({len(long)} rows)")

    # per-compartment summary
    for compartment in COMPARTMENTS_ALL:
        sub = long[(long["compartment"] == compartment)
                       & long["delta_aligned"].notna()]
        if len(sub) == 0:
            continue
        ra2la = sub[sub["healthy_direction"] == "RA-biased"]
        la2ra = sub[sub["healthy_direction"] == "LA-biased"]
        included = bool(sub["compartment_included"].iloc[0])
        summary_rows.append({
            "compartment":              compartment,
            "compartment_included":     included,
            "n_genes_evaluated":        len(sub),
            "mean_delta_aligned":       round(sub["delta_aligned"].mean(), 5),
            "median_delta_aligned":     round(sub["delta_aligned"].median(), 5),
            "frac_genes_aligned_pos":   round((sub["delta_aligned"] > 0).mean(), 3),
            "n_genes_ra2la":            len(ra2la),
            "mean_delta_aligned_ra2la": (round(ra2la["delta_aligned"].mean(), 5)
                                          if len(ra2la) else np.nan),
            "frac_aligned_pos_ra2la":   (round((ra2la["delta_aligned"] > 0).mean(), 3)
                                          if len(ra2la) else np.nan),
            "n_genes_la2ra":            len(la2ra),
            "mean_delta_aligned_la2ra": (round(la2ra["delta_aligned"].mean(), 5)
                                          if len(la2ra) else np.nan),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)} ({len(summary)} rows)")

    print("\n=== per-compartment carrier summary (sorted by mean_delta_aligned) ===")
    print(summary.sort_values("mean_delta_aligned", ascending=False)
                  [["compartment", "compartment_included",
                    "n_genes_evaluated", "mean_delta_aligned",
                    "frac_genes_aligned_pos",
                    "mean_delta_aligned_ra2la"]].to_string(index=False))


if __name__ == "__main__":
    main()
