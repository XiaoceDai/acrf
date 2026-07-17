#!/usr/bin/env python3
"""Freeze corrected GSE255612 family-level donor×compartment scores."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
SCORES = RESULTS / "fig3_v3_singscore_per_donor.tsv"
FAMILY_METADATA = RESULTS / "ncp_v3_family_labels.tsv"
AUTHORITY_LOCK = RESULTS / "fig3_publication_authority_lock.md"
OUT = (
    RESULTS
    / "ncp_v3_singscore_reanalysis"
    / "gse255612_family_donor_compartment_scores.tsv"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def natural_key(value: object) -> tuple[object, ...]:
    return tuple(
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", str(value))
    )


def main() -> None:
    require(
        SCORES.is_file() and FAMILY_METADATA.is_file() and AUTHORITY_LOCK.is_file(),
        "A corrected family-score bridge input is missing",
    )
    require(
        "LOCKED" in AUTHORITY_LOCK.read_text(encoding="utf-8"),
        "Figure 3 publication authority is not locked",
    )
    scores = pd.read_csv(SCORES, sep="\t")
    metadata = pd.read_csv(FAMILY_METADATA, sep="\t")
    require(
        scores["scorer"].eq("python_exact_reference_matching_singscore_1.28.1").all(),
        "Corrected scorer changed",
    )
    require(
        scores["le_source"].eq("frozen_program_metadata").all(),
        "Frozen family-program source changed",
    )
    require(set(scores["rhythm"]) == {"Ctrl", "AF"}, "Unexpected rhythm groups")
    expected_programs = (
        scores.groupby("family_id")["program_id"].nunique().astype(int).to_dict()
    )
    family = (
        scores.loc[scores["score"].notna()]
        .groupby(
            ["family_id", "direction", "compartment", "donor_id", "rhythm"],
            as_index=False,
        )
        .agg(
            family_score=("score", "mean"),
            n_constituent_programs=("program_id", "nunique"),
        )
    )
    require(
        family.apply(
            lambda row: int(row["n_constituent_programs"])
            == expected_programs[int(row["family_id"])],
            axis=1,
        ).all(),
        "A donor×family×compartment row omits a constituent family program",
    )
    labels = metadata[["family_id", "family_label", "short_label"]].copy()
    family = family.merge(labels, on="family_id", how="left", validate="many_to_one")
    require(
        family[["family_label", "short_label"]].notna().all().all(),
        "Publication-facing family metadata is incomplete",
    )
    family = family.rename(columns={"rhythm": "group"})
    family["source_donor_scores_sha256"] = sha256(SCORES)
    family["source_family_metadata_sha256"] = sha256(FAMILY_METADATA)
    family["source_authority_lock_sha256"] = sha256(AUTHORITY_LOCK)
    family["_group_order"] = family["group"].map({"Ctrl": 0, "AF": 1})
    family["_donor_key"] = family["donor_id"].map(natural_key)
    family = family.sort_values(
        ["family_id", "compartment", "_group_order", "_donor_key"]
    ).reset_index(drop=True)
    family = family[
        [
            "family_id",
            "family_label",
            "short_label",
            "direction",
            "compartment",
            "donor_id",
            "group",
            "family_score",
            "n_constituent_programs",
            "source_donor_scores_sha256",
            "source_family_metadata_sha256",
            "source_authority_lock_sha256",
        ]
    ]
    require(family["family_id"].nunique() == 16, "Expected 16 frozen families")
    require(family["compartment"].nunique() == 9, "Expected nine compartments")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    family.to_csv(OUT, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {OUT.relative_to(ROOT)} ({len(family)} family-level rows)")
    print(f"SHA256 {sha256(OUT)}")


if __name__ == "__main__":
    main()
