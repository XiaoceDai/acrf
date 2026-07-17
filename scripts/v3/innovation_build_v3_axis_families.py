#!/usr/bin/env python3
"""
Build a family-mapping table for the innovation NCP panel.

Policy
------
- If a v3 pathway exists in the baseline family table, keep its family_id.
- If missing from baseline, assign a new family_id starting from max+1.

Output
------
results/innovation_ncp_v3_axis_families.tsv
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PANEL_V3 = ROOT / "results" / "innovation_ncp_v3_panel.tsv"
FAM_BASE = ROOT / "results" / "fig1e_baseline_axis_families.tsv"
OUT = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"


def main() -> None:
    v3 = pd.read_csv(PANEL_V3, sep="\t")
    fam = pd.read_csv(FAM_BASE, sep="\t")

    fam_sub = fam[["pathway", "family_id"]].drop_duplicates()
    merged = v3[["program_order", "pathway", "direction"]].merge(
        fam_sub, on="pathway", how="left"
    )

    next_id = int(fam_sub["family_id"].max()) + 1
    assigned = []
    for i, row in merged.iterrows():
        if pd.isna(row["family_id"]):
            merged.loc[i, "family_id"] = next_id
            assigned.append((row["pathway"], next_id))
            next_id += 1

    merged["family_id"] = merged["family_id"].astype(int)
    merged["family_source"] = merged["pathway"].isin(set(fam_sub["pathway"])).map(
        {True: "baseline", False: "new_v3"}
    )
    merged = merged.sort_values(["family_id", "program_order"]).reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT, sep="\t", index=False)

    print(f"wrote {OUT.relative_to(ROOT)} ({len(merged)} rows)")
    print(f"new families assigned: {len(assigned)}")
    for p, fid in assigned:
        print(f"  F{fid}: {p}")


if __name__ == "__main__":
    main()
