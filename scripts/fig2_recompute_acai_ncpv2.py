#!/usr/bin/env python3
"""
Fig 2 prep — recompute donor-level aCAI per NCP_v2 program × 3 paired AF
cohorts, using the new 30-program NCP_v2 panel (results/fig1c_ncp_panel_v2.tsv).

aCAI(d, p) = mean log expression of program p member genes in LA donor d
           - mean log expression in RA donor d

Cohorts:
    GSE128188  (5 AF + 5 SR donors,  paired LA/RA, RNA-seq → log2-CPM)
    GSE79768   (7 AF + 6 SR donors,  paired LA/RA, RMA-normalised microarray)
    GSE115574  (15 AF + 16 SR donors, paired LA/RA, RMA-normalised microarray)

Outputs:
    results/fig2_primary_acai_long.tsv     (long format: donor x program)
    results/fig2_primary_acai_summary.tsv  (per program × cohort: mean / sd
                                            in AF and SR + Δ aCAI(AF − SR))
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

NCP_TSV = ROOT / "results" / "fig1c_ncp_panel_v2.tsv"
GMT_HALL = ROOT / "bulk" / "Project_P01" / "gene_sets" \
                 / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC = ROOT / "bulk" / "Project_P01" / "gene_sets" \
                 / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"

OUT_LONG    = ROOT / "results" / "fig2_primary_acai_long.tsv"
OUT_SUMMARY = ROOT / "results" / "fig2_primary_acai_summary.tsv"


def load_gmt(path: Path, prefix: str) -> dict[str, list[str]]:
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            out[f"{prefix}__{parts[0]}"] = [g for g in parts[2:] if g]
    return out


def load_paired_cohort(cohort: str
                       ) -> tuple[pd.DataFrame,
                                    list[tuple[str, str, str, str]]]:
    """Return log-scale expression matrix and paired-sample list:
    [(LA_gsm, RA_gsm, donor_id, rhythm), ...]."""
    expr = pd.read_csv(ROOT / "Project_P1" / cohort / "expr_symbol.tsv",
                        sep="\t", index_col=0)
    expr = expr.groupby(expr.index).max()
    md = pd.read_csv(ROOT / "Project_P1" / cohort / "sample_metadata.tsv",
                     sep="\t")
    md = md[md["gsm"].isin(expr.columns)]

    if expr.values.max() > 100:
        expr = np.log2(expr.clip(lower=0) + 1.0)

    pairs = []
    for donor, sub in md.groupby("donor_id"):
        la_rows = sub[sub["chamber"] == "LA"]
        ra_rows = sub[sub["chamber"] == "RA"]
        if len(la_rows) >= 1 and len(ra_rows) >= 1:
            rhythm = la_rows["rhythm"].iloc[0]
            pairs.append((la_rows["gsm"].iloc[0],
                          ra_rows["gsm"].iloc[0],
                          str(donor), rhythm))
    return expr, pairs


def main() -> None:
    print("Loading NCP_v2 panel + pathway gene sets ...")
    ncp = pd.read_csv(NCP_TSV, sep="\t")
    pathways = {**load_gmt(GMT_HALL, "HALLMARK"),
                **load_gmt(GMT_REAC, "REACTOME")}
    print(f"  NCP_v2 programs: {len(ncp)}")

    rows = []
    for cohort in ["GSE128188", "GSE79768", "GSE115574"]:
        print(f"Loading paired cohort {cohort} ...")
        expr, pairs = load_paired_cohort(cohort)
        n_af = sum(1 for p in pairs if p[3] == "AF")
        n_sr = sum(1 for p in pairs if p[3] == "SR")
        print(f"  {cohort}: {len(pairs)} paired donors "
              f"({n_af} AF + {n_sr} SR), {expr.shape[0]} genes")

        for _, p in ncp.iterrows():
            gset = pathways.get(p["pathway"], [])
            present = [g for g in gset if g in expr.index]
            if len(present) < 5:
                continue
            sub = expr.loc[present]
            for la_gsm, ra_gsm, donor, rhythm in pairs:
                la_mean = sub[la_gsm].mean()
                ra_mean = sub[ra_gsm].mean()
                aCAI = la_mean - ra_mean
                rows.append({
                    "cohort":         cohort,
                    "donor_id":       donor,
                    "rhythm":         rhythm,
                    "program_id":     p["pathway"],
                    "program_order":  int(p["program_order"]),
                    "direction":      p["direction"],
                    "score_LA":       round(float(la_mean), 4),
                    "score_RA":       round(float(ra_mean), 4),
                    "aCAI":           round(float(aCAI), 4),
                    "n_genes_used":   len(present),
                })
    long = pd.DataFrame(rows)
    OUT_LONG.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(OUT_LONG, sep="\t", index=False)
    print(f"\nwrote {OUT_LONG.relative_to(ROOT)}  ({len(long)} rows)")

    # Per program × cohort × rhythm aggregate
    summ = (long.groupby(["program_id", "program_order", "direction",
                            "cohort", "rhythm"])
                 .agg(n_donors=("donor_id", "nunique"),
                       aCAI_mean=("aCAI", "mean"),
                       aCAI_sd=("aCAI", "std"))
                 .reset_index())
    # Pivot to AF / SR side by side, compute Δ
    pv = (summ.pivot_table(index=["program_id", "program_order",
                                    "direction", "cohort"],
                            columns="rhythm",
                            values=["n_donors", "aCAI_mean",
                                     "aCAI_sd"]))
    pv.columns = [f"{a}_{b}" for a, b in pv.columns]
    pv = pv.reset_index()
    pv["delta_aCAI_AF_minus_SR"] = (pv["aCAI_mean_AF"]
                                       - pv["aCAI_mean_SR"]).round(4)
    pv["aCAI_mean_AF"] = pv["aCAI_mean_AF"].round(4)
    pv["aCAI_mean_SR"] = pv["aCAI_mean_SR"].round(4)
    pv["aCAI_sd_AF"]   = pv["aCAI_sd_AF"].round(4)
    pv["aCAI_sd_SR"]   = pv["aCAI_sd_SR"].round(4)
    pv = pv.sort_values(["cohort", "program_order"]).reset_index(drop=True)
    pv.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)}  ({len(pv)} rows)")

    print("\n=== summary preview ===")
    print(pv[["cohort", "program_order", "direction",
              "n_donors_AF", "n_donors_SR",
              "aCAI_mean_AF", "aCAI_mean_SR",
              "delta_aCAI_AF_minus_SR"]].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
