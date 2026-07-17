#!/usr/bin/env python3
"""Compare frozen coloc table vs latest h4d rerun (diff-only; no unfreeze)."""
from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
FROZEN = ROOT / "results" / "figure6_mr_anchors_v3_pipeline_wide.tsv"
OUT = ROOT / "results" / "coloc_rerun_diff_report.tsv"

COLS = ["lead_eqtl_p", "gwas_min_p", "coloc_PP3", "coloc_PP4", "gate_passed"]


def main() -> None:
    reruns = sorted(glob.glob(str(ROOT / "results" / "figure6_mr_anchors_v3_rerun_*.tsv")))
    if not reruns:
        print("No rerun TSV found — run scripts/h4d_mr_coloc_v3_all_anchors.R first.")
        return
    rerun_path = reruns[-1]
    frozen = pd.read_csv(FROZEN, sep="\t")
    rerun = pd.read_csv(rerun_path, sep="\t")
    merged = frozen.merge(rerun, on="anchor", suffixes=("_frozen", "_rerun"), how="outer")

    rows = []
    for _, r in merged.iterrows():
        for c in COLS:
            a = r.get(f"{c}_frozen")
            b = r.get(f"{c}_rerun")
            if pd.isna(a) and pd.isna(b):
                continue
            try:
                diff = float(b) - float(a) if pd.notna(a) and pd.notna(b) else np.nan
            except (TypeError, ValueError):
                diff = np.nan
            rows.append({
                "anchor": r["anchor"],
                "field": c,
                "frozen": a,
                "rerun": b,
                "abs_diff": abs(diff) if pd.notna(diff) else np.nan,
                "rerun_file": Path(rerun_path).name,
            })
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"wrote {OUT.relative_to(ROOT)} from {Path(rerun_path).name}")
    if len(df):
        print(f"  max |diff| PP4: {df.loc[df['field']=='coloc_PP4','abs_diff'].max():.4g}")


if __name__ == "__main__":
    main()
