#!/usr/bin/env python3
"""Freeze GSE41177 duration metadata for corrected Figure 4 analysis."""
from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
META_IN = ROOT / "Project_P4" / "GSE41177" / "sample_metadata.tsv"
OUT = ROOT / "results" / "ncp_v3_singscore_reanalysis" / "gse41177_duration_metadata.tsv"


def exclusion_reason(row: pd.Series) -> str:
    if row["rhythm"] == "SR":
        return "sr_anchor"
    if row["tissue"] == "LAPV":
        return "lapv_only"
    if row["rhythm"] == "AF" and row["tissue"] == "LAA":
        if pd.isna(row["duration_months"]) or str(row["duration_raw"]).strip() == ">1M":
            return "censored_gt1m"
        return "eligible_primary"
    return "other"


def main() -> None:
    meta = pd.read_csv(META_IN, sep="\t")
    meta["duration_months"] = pd.to_numeric(meta["duration_months"], errors="coerce")
    out = pd.DataFrame(
        {
            "sample_id": meta["gsm"],
            "donor_id": meta["donor_id"],
            "region": meta["tissue"],
            "rhythm": meta["rhythm"],
            "duration_raw": meta["duration_raw"],
            "duration_months": meta["duration_months"],
            "primary_duration_eligible": (
                (meta["rhythm"] == "AF")
                & (meta["tissue"] == "LAA")
                & meta["duration_months"].notna()
            ),
            "exclusion_reason": meta.apply(exclusion_reason, axis=1),
        }
    )

    assert len(out) == 38, f"expected 38 samples, got {len(out)}"
    n_primary = int(out["primary_duration_eligible"].sum())
    assert n_primary == 12, f"expected 12 primary-eligible AF-LAA, got {n_primary}"

    n_censored = int(
        (
            (out["rhythm"] == "AF")
            & (out["region"] == "LAA")
            & (out["exclusion_reason"] == "censored_gt1m")
        ).sum()
    )
    assert n_censored == 4, f"expected 4 censored AF-LAA, got {n_censored}"

    n_sr_laa = int(((out["rhythm"] == "SR") & (out["region"] == "LAA")).sum())
    assert n_sr_laa == 3, f"expected 3 SR-LAA, got {n_sr_laa}"

    # Locked D1 contract: no year-converted durations (all end with M or are 0M/>1M/SR)
    bad = out[
        out["duration_raw"].astype(str).str.contains(r"\d+[Yy]", regex=True, na=False)
    ]
    assert bad.empty, "year-converted durations detected"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(f"wrote {OUT.relative_to(ROOT)} ({len(out)} rows, {n_primary} primary-eligible)")


if __name__ == "__main__":
    main()
