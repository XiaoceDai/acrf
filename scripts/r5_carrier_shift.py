#!/usr/bin/env python3
"""
R5 - 3C carrier shift table (three columns).

Spec: carrier-shift summary (paired with r4_delta_classify outputs).

Column 1  healthy_dominant_CT         (Fig 1D signed-attribution argmax)
Column 2  healthy_top_expressor_CT    (argmax over cell types of
                                        mean singscore on Fig 1 pseudobulk,
                                        AF-naive i.e. all 16-24 healthy LA+RA
                                        pseudobulks)
Column 3  AF_major_carrier_CT         (from R4, §3.7 with q_perm < 0.10)

These three measure different quantities and must NOT be collapsed into
one "dominant cell type" summary; R5 lays them side-by-side so the Fig
3C panel can render them as three columns of a single row per program.

Inputs
    results/normal_chamber_program_attribution.tsv   (Fig 1D)
    Project_P3/fig1_healthy_program_scores.tsv       (R3)
    results/af_la_major_carrier_af.tsv               (R4)
    results/figure3_program_panel.tsv                (R3)

Output
    results/figure3_carrier_shift.tsv                (per-program, §6.3 schema)
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    panel = pd.read_csv(ROOT / "results" / "figure3_program_panel.tsv", sep="\t")
    panel = panel.sort_values("program_order").reset_index(drop=True)

    # ---- col 1: healthy_dominant_CT from Fig 1D attribution ----
    attr = pd.read_csv(ROOT / "results" / "normal_chamber_program_attribution.tsv",
                       sep="\t")
    attr_col1 = attr[["program_id", "dominant_cell_type", "attribution_strength"]].copy()
    attr_col1 = attr_col1.rename(columns={
        "dominant_cell_type":   "healthy_dominant_CT",
        "attribution_strength": "healthy_dominant_attribution_strength",
    })

    # ---- col 2: healthy_top_expressor_CT from Fig 1 pseudobulk singscore ----
    fig1 = pd.read_csv(ROOT / "Project_P3" / "fig1_healthy_program_scores.tsv", sep="\t")
    # mean singscore per (program_id, compartment4) across all healthy samples
    mean_by_ct = (fig1.groupby(["program_id", "compartment4"], as_index=False)
                      ["score"].mean())
    # for each program_id, pick ct with max mean score
    def _argmax_ct(grp):
        i = grp["score"].idxmax()
        return pd.Series({
            "healthy_top_expressor_CT":    grp.loc[i, "compartment4"],
            "healthy_top_expressor_score": float(grp.loc[i, "score"]),
        })
    col2 = mean_by_ct.groupby("program_id", as_index=False).apply(_argmax_ct)
    # reshape: apply may return MultiIndex in recent pandas, normalise
    col2 = col2.reset_index(drop=True)

    # ---- col 3: AF_major_carrier_CT from R4 ----
    carrier = pd.read_csv(ROOT / "results" / "af_la_major_carrier_af.tsv", sep="\t")
    col3 = carrier[["program_id",
                    "AF_major_carrier_CT",
                    "AF_major_carrier_dScore",
                    "AF_major_carrier_q_perm"]].copy()

    # ---- assemble ----
    out = panel[["program_id", "program_label", "program_order",
                 "baseline_direction"]].copy()
    out = (out.merge(attr_col1, on="program_id", how="left")
              .merge(col2,      on="program_id", how="left")
              .merge(col3,      on="program_id", how="left"))

    # carrier_shift_flag: top expressor -> AF major carrier
    def _flag(row):
        af = row["AF_major_carrier_CT"]
        hop = row["healthy_top_expressor_CT"]
        if pd.isna(af) or af == "Unresolved" or pd.isna(hop):
            return pd.NA
        return af != hop
    out["carrier_shift_flag_vs_top_expressor"] = out.apply(_flag, axis=1)

    def _flag_dominant(row):
        af  = row["AF_major_carrier_CT"]
        dom = row["healthy_dominant_CT"]
        if pd.isna(af) or af == "Unresolved" or pd.isna(dom) or dom == "unresolved":
            return pd.NA
        return af != dom
    out["carrier_shift_flag_vs_dominant"] = out.apply(_flag_dominant, axis=1)

    out = out.sort_values("program_order").reset_index(drop=True)

    out_path = ROOT / "results" / "figure3_carrier_shift.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"wrote {out_path.relative_to(ROOT)}  ({len(out)} rows)")

    # stdout preview
    print("\n=== carrier shift table (3C source) ===")
    print(f"{'#':>2}  {'program':28s}  {'dir':3s}  "
          f"{'healthy_dominant (Fig1D)':25s}  {'top_expressor (Fig1 ss)':25s}  {'AF_major (Fig3 Δ)':25s}")
    for _, r in out.iterrows():
        print(f"{r['program_order']:>2}  {r['program_label']:28s}  "
              f"{r['baseline_direction']:+d}   "
              f"{str(r['healthy_dominant_CT']):25s}  "
              f"{str(r['healthy_top_expressor_CT']):25s}  "
              f"{str(r['AF_major_carrier_CT']):25s}")


if __name__ == "__main__":
    main()
