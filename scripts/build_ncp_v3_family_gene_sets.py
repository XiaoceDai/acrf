#!/usr/bin/env python3
"""Prep for Fig 6C (NCP_v3 family x published AF risk-gene overlap)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROGRAM_LE = ROOT / "results" / "fig3_v3_program_le.tsv"
FAMILY_LABELS = ROOT / "results" / "ncp_v3_family_labels.tsv"
BULK_DE_FULL = ROOT / "bulk" / "Project_P01" / "bulk_chamber_DE_full.tsv"

OUT_LONG = ROOT / "results" / "ncp_v3_family_leading_edge_genes.tsv"
OUT_SUMMARY = ROOT / "results" / "ncp_v3_family_gene_set_summary.tsv"
OUT_UNIVERSE = ROOT / "results" / "ncp_v3_gene_universe.tsv"


def load_universe() -> set[str]:
    de = pd.read_csv(BULK_DE_FULL, sep="\t")
    genes = de["gene_symbol"].astype(str).str.strip().str.strip('"')
    return set(g for g in genes if g and g.lower() != "nan")


def load_family_labels() -> dict[int, dict[str, str]]:
    fb = pd.read_csv(FAMILY_LABELS, sep="\t")
    out: dict[int, dict[str, str]] = {}
    for _, r in fb.iterrows():
        out[int(r["family_id"])] = {
            "family_label": str(r.get("family_label", "")),
            "short_label": str(r.get("short_label", "")),
            "axis_group": str(r.get("axis_group", "")),
        }
    return out


def main() -> None:
    universe = load_universe()
    labels = load_family_labels()

    le = pd.read_csv(PROGRAM_LE, sep="\t")

    long_rows: list[dict] = []
    for _, r in le.iterrows():
        fid = int(r["family_id"])
        genes = [g.strip() for g in str(r["leading_edge_genes"]).split(",") if g.strip()]
        lab = labels.get(fid, {})
        for g in genes:
            long_rows.append({
                "family_id": fid,
                "family_label": lab.get("family_label", ""),
                "short_label": lab.get("short_label", ""),
                "axis_group": lab.get("axis_group", ""),
                "program_id": r["program_id"],
                "direction": r["direction"],
                "gene": g,
                "in_universe": g in universe,
            })

    long_df = pd.DataFrame(long_rows)
    # De-duplicate to unique family x gene (a gene may recur across programs).
    fam_gene = (
        long_df.sort_values(["family_id", "gene"])
        .drop_duplicates(subset=["family_id", "gene"])
        .reset_index(drop=True)
    )
    fam_gene.to_csv(OUT_LONG, sep="\t", index=False)

    # Per-family summary.
    summary_rows = []
    for fid, grp in fam_gene.groupby("family_id"):
        lab = labels.get(int(fid), {})
        n_genes = grp["gene"].nunique()
        n_in_univ = int(grp["in_universe"].sum())
        n_programs = long_df.loc[long_df["family_id"] == fid, "program_id"].nunique()
        summary_rows.append({
            "family_id": int(fid),
            "family_label": lab.get("family_label", ""),
            "short_label": lab.get("short_label", ""),
            "axis_group": lab.get("axis_group", ""),
            "n_programs": n_programs,
            "n_leading_edge_genes": n_genes,
            "n_in_universe": n_in_univ,
        })
    summary_df = pd.DataFrame(summary_rows).sort_values("family_id")
    summary_df.to_csv(OUT_SUMMARY, sep="\t", index=False)

    # Background universe file.
    uni_df = pd.DataFrame({"gene": sorted(universe)})
    uni_df["source"] = "bulk_chamber_DE_full (expressed / FGSEA-ranked universe)"
    uni_df.to_csv(OUT_UNIVERSE, sep="\t", index=False)

    total_le = fam_gene["gene"].nunique()
    total_in_univ = fam_gene.loc[fam_gene["in_universe"], "gene"].nunique()
    print(f"wrote {OUT_LONG.relative_to(ROOT)} "
          f"({len(fam_gene)} family x gene rows)")
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)} "
          f"({len(summary_df)} families)")
    print(f"wrote {OUT_UNIVERSE.relative_to(ROOT)} "
          f"({len(universe)} background genes)")
    print(f"leading-edge gene coverage in universe: "
          f"{total_in_univ}/{total_le} unique genes")


if __name__ == "__main__":
    main()
