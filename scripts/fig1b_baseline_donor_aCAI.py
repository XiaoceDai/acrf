#!/usr/bin/env python3
"""Healthy-baseline donor-level aCAI table (GSE112339; not a main-text figure)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

AXIS_FAM_TSV = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
COUNTS = ROOT / "bulk" / "Project_P01" / "GSE112339" / "GSE112339_HTSeq_hg19_Ensembl.txt.gz"
GMT_HALL = ROOT / "bulk" / "Project_P01" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC = ROOT / "bulk" / "Project_P01" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
OUT_TSV = ROOT / "results" / "fig1b_baseline_donor_aCAI.tsv"


def load_gmt(path: Path, prefix: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = f"{prefix}__{parts[0]}"
            genes = [g for g in parts[2:] if g]
            out[name] = genes
    return out


def load_log_cpm() -> pd.DataFrame:
    raw = pd.read_csv(COUNTS, sep="\t", index_col=0)
    raw = raw.groupby(raw.index).max()
    keep = [c for c in raw.columns if c.startswith(("LA_", "RA_"))]
    raw = raw[keep]
    cpm = raw.div(raw.sum(axis=0), axis=1) * 1e6
    return np.log2(cpm + 1.0)


def main() -> None:
    print("Loading NCP deployed panel and reference data ...")
    if not AXIS_FAM_TSV.exists():
        raise FileNotFoundError(
            f"Missing {AXIS_FAM_TSV.relative_to(ROOT)} — build it first "
            "(e.g. scripts/v3/innovation_build_v3_axis_families.py)."
        )
    ncp = pd.read_csv(AXIS_FAM_TSV, sep="\t")[["program_order", "pathway", "direction"]].copy()
    print(
        f"  NCP panel: {len(ncp)} programs "
        f"({(ncp['direction'] == 'LA-biased').sum()} LA, "
        f"{(ncp['direction'] == 'RA-biased').sum()} RA)"
    )

    pathways = {**load_gmt(GMT_HALL, "HALLMARK"), **load_gmt(GMT_REAC, "REACTOME")}
    log_cpm = load_log_cpm()
    print(f"  log-CPM: {log_cpm.shape[0]} genes x {log_cpm.shape[1]} samples")

    la_cols = sorted(c for c in log_cpm.columns if c.startswith("LA_"))
    donors = [c[3:] for c in la_cols]
    print(f"  paired donors: {len(donors)}  ({', '.join(donors)})")

    rows = []
    for _, r in ncp.iterrows():
        gset = [g for g in pathways.get(r["pathway"], []) if g in log_cpm.index]
        if len(gset) < 5:
            raise ValueError(
                f"{r['pathway']}: only {len(gset)} genes overlap counts (require ≥5)."
            )
        sub = log_cpm.loc[gset]
        per_donor = [
            sub[f"LA_{d}"].mean() - sub[f"RA_{d}"].mean()
            for d in donors
        ]
        rows.append({
            "program_order": int(r["program_order"]),
            "pathway": r["pathway"],
            "direction": r["direction"],
            "n_genes_used": len(gset),
            **{f"aCAI_{d}": v for d, v in zip(donors, per_donor)},
        })

    out = pd.DataFrame(rows)
    donor_cols = [f"aCAI_{d}" for d in donors]
    out["mean_aCAI"] = out[donor_cols].mean(axis=1)
    out["sign_concordance"] = (
        np.sign(out[donor_cols]).eq(np.sign(out["mean_aCAI"]), axis=0).sum(axis=1).astype(int)
    )

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")

    print("\n=== NCP baseline aCAI summary ===")
    print("sign concordance distribution (out of 8 donors):")
    print(out["sign_concordance"].value_counts().sort_index().to_string())
    print(
        f"\nprograms with 8/8 donor sign concordance: "
        f"{(out['sign_concordance'] == 8).sum()} / {len(out)}"
    )


if __name__ == "__main__":
    main()
