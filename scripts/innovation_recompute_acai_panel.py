#!/usr/bin/env python3
"""
Recompute donor-level aCAI for an arbitrary pathway panel.

Usage
-----
python3 scripts/innovation_recompute_acai_panel.py \
  --panel results/fig1c_ncp_panel_v2.tsv \
  --tag innovation_ncp_v2

python3 scripts/innovation_recompute_acai_panel.py \
  --panel results/innovation_ncp_v3_panel.tsv \
  --tag innovation_ncp_v3

Outputs
-------
results/{tag}_acai_long.tsv
results/{tag}_acai_summary.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

GMT_HALL = ROOT / "bulk" / "Project_P01" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC = ROOT / "bulk" / "Project_P01" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
COHORTS = ["GSE128188", "GSE79768", "GSE115574"]


def load_gmt(path: Path, prefix: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            out[f"{prefix}__{parts[0]}"] = [g for g in parts[2:] if g]
    return out


def load_paired_cohort(cohort: str) -> tuple[pd.DataFrame, list[tuple[str, str, str, str]]]:
    expr = pd.read_csv(
        ROOT / "Project_P1" / cohort / "expr_symbol.tsv",
        sep="\t",
        index_col=0,
    )
    expr = expr.groupby(expr.index).max()
    md = pd.read_csv(ROOT / "Project_P1" / cohort / "sample_metadata.tsv", sep="\t")
    md = md[md["gsm"].isin(expr.columns)]

    if expr.values.max() > 100:
        expr = np.log2(expr.clip(lower=0) + 1.0)

    pairs: list[tuple[str, str, str, str]] = []
    for donor, sub in md.groupby("donor_id"):
        la_rows = sub[sub["chamber"] == "LA"]
        ra_rows = sub[sub["chamber"] == "RA"]
        if len(la_rows) >= 1 and len(ra_rows) >= 1:
            rhythm = str(la_rows["rhythm"].iloc[0])
            pairs.append(
                (
                    str(la_rows["gsm"].iloc[0]),
                    str(ra_rows["gsm"].iloc[0]),
                    str(donor),
                    rhythm,
                )
            )
    return expr, pairs


def main() -> None:
    parser = argparse.ArgumentParser(description="Recompute donor-level aCAI for a panel.")
    parser.add_argument("--panel", required=True, help="Path to panel TSV.")
    parser.add_argument("--tag", required=True, help="Output tag, e.g. innovation_ncp_v3.")
    args = parser.parse_args()

    panel_path = Path(args.panel)
    if not panel_path.is_absolute():
        panel_path = ROOT / panel_path
    if not panel_path.exists():
        raise FileNotFoundError(f"Panel file not found: {panel_path}")

    ncp = pd.read_csv(panel_path, sep="\t")
    if "pathway" not in ncp.columns:
        raise ValueError(f"Panel file missing 'pathway' column: {panel_path}")
    if "program_order" not in ncp.columns:
        ncp = ncp.copy()
        ncp["program_order"] = np.arange(1, len(ncp) + 1)
    if "direction" not in ncp.columns:
        raise ValueError(f"Panel file missing 'direction' column: {panel_path}")

    out_long = ROOT / "results" / f"{args.tag}_acai_long.tsv"
    out_summary = ROOT / "results" / f"{args.tag}_acai_summary.tsv"

    print(f"Loading panel: {panel_path.relative_to(ROOT)} (K={len(ncp)})")
    pathways = {
        **load_gmt(GMT_HALL, "HALLMARK"),
        **load_gmt(GMT_REAC, "REACTOME"),
    }

    rows: list[dict[str, object]] = []
    for cohort in COHORTS:
        print(f"Loading paired cohort {cohort} ...")
        expr, pairs = load_paired_cohort(cohort)
        n_af = sum(1 for p in pairs if p[3] == "AF")
        n_sr = sum(1 for p in pairs if p[3] == "SR")
        print(f"  {cohort}: {len(pairs)} paired donors ({n_af} AF + {n_sr} SR), {expr.shape[0]} genes")

        for _, p in ncp.iterrows():
            pathway = str(p["pathway"])
            gset = pathways.get(pathway, [])
            present = [g for g in gset if g in expr.index]
            if len(present) < 5:
                continue
            sub = expr.loc[present]
            for la_gsm, ra_gsm, donor, rhythm in pairs:
                la_mean = float(sub[la_gsm].mean())
                ra_mean = float(sub[ra_gsm].mean())
                acai = la_mean - ra_mean
                rows.append(
                    {
                        "cohort": cohort,
                        "donor_id": donor,
                        "rhythm": rhythm,
                        "program_id": pathway,
                        "program_order": int(p["program_order"]),
                        "direction": str(p["direction"]),
                        "score_LA": round(la_mean, 4),
                        "score_RA": round(ra_mean, 4),
                        "aCAI": round(acai, 4),
                        "n_genes_used": len(present),
                    }
                )

    long = pd.DataFrame(rows)
    out_long.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(out_long, sep="\t", index=False)
    print(f"wrote {out_long.relative_to(ROOT)} ({len(long)} rows)")

    summ = (
        long.groupby(["program_id", "program_order", "direction", "cohort", "rhythm"])
        .agg(n_donors=("donor_id", "nunique"), aCAI_mean=("aCAI", "mean"), aCAI_sd=("aCAI", "std"))
        .reset_index()
    )
    pv = (
        summ.pivot_table(
            index=["program_id", "program_order", "direction", "cohort"],
            columns="rhythm",
            values=["n_donors", "aCAI_mean", "aCAI_sd"],
        )
    )
    pv.columns = [f"{a}_{b}" for a, b in pv.columns]
    pv = pv.reset_index()
    pv["delta_aCAI_AF_minus_SR"] = pv["aCAI_mean_AF"] - pv["aCAI_mean_SR"]
    for c in ["aCAI_mean_AF", "aCAI_mean_SR", "aCAI_sd_AF", "aCAI_sd_SR", "delta_aCAI_AF_minus_SR"]:
        pv[c] = pv[c].round(4)
    pv = pv.sort_values(["cohort", "program_order"]).reset_index(drop=True)
    pv.to_csv(out_summary, sep="\t", index=False)
    print(f"wrote {out_summary.relative_to(ROOT)} ({len(pv)} rows)")


if __name__ == "__main__":
    main()
