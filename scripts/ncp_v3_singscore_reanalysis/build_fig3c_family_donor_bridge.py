#!/usr/bin/env python3
"""Freeze corrected family-level donor scores for final Figure 3C."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
OUT = RESULTS / "ncp_v3_singscore_reanalysis" / "fig3c_fidelity_relay_donor_source.tsv"

SCORES = RESULTS / "fig3_v3_singscore_per_donor.tsv"
FAMILY_DELTA = RESULTS / "fig3_v3_family_celltype_delta.tsv"
CLASS_TABLE = RESULTS / "fig3_v3_carrier_fidelity_class.tsv"
HCA_TABLE = RESULTS / "fig1f_v3_HCA_family_carriers.tsv"
FIG3B_SOURCE = (
    RESULTS
    / "ncp_v3_singscore_reanalysis"
    / "fig3b_family_compartment_localization_source.tsv"
)
AUTHORITY_LOCK = RESULTS / "fig3_publication_authority_lock.md"

TARGETS = [
    (10, "fibroblast", "GAG disease / TGF-β"),
    (11, "endothelial", "BMP signalling"),
    (11, "adipocyte", "BMP signalling"),
]
EXPECTED_CLASS = {10: "Fidelity", 11: "Relay"}
EXPECTED_CARRIER = {10: "fibroblast", 11: "endothelial"}
EXPECTED_SUPPORT = {10: {"fibroblast"}, 11: {"adipocyte"}}

OUTPUT_COLUMNS = [
    "family_id",
    "family_label",
    "fidelity_class",
    "compartment",
    "compartment_label",
    "donor_id",
    "group",
    "family_score",
    "n_family_programs",
    "is_carrier",
    "is_supported",
    "n_ctrl",
    "n_af",
    "mean_ctrl",
    "mean_af",
    "delta",
    "p_perm",
    "n_perm",
    "perm_strategy",
    "source_donor_scores_sha256",
    "source_family_delta_sha256",
    "source_class_table_sha256",
    "source_hca_table_sha256",
    "source_fig3b_sha256",
    "source_authority_lock_sha256",
]


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


def build_bridge() -> pd.DataFrame:
    allowed_inputs = [
        SCORES,
        FAMILY_DELTA,
        CLASS_TABLE,
        HCA_TABLE,
        FIG3B_SOURCE,
        AUTHORITY_LOCK,
    ]
    require(all(path.is_file() for path in allowed_inputs), "A corrected bridge input is missing")
    require("LOCKED" in AUTHORITY_LOCK.read_text(encoding="utf-8"), "Figure 3 authority is not locked")
    require(
        not any(
            token in path.name.lower()
            for path in allowed_inputs
            for token in ("gene", "pre_correction", "legacy")
        ),
        "Forbidden gene-level, pre-correction, or legacy bridge input",
    )

    scores = pd.read_csv(SCORES, sep="\t")
    delta = pd.read_csv(FAMILY_DELTA, sep="\t")
    classes = pd.read_csv(CLASS_TABLE, sep="\t")
    hca = pd.read_csv(HCA_TABLE, sep="\t")
    fig3b = pd.read_csv(FIG3B_SOURCE, sep="\t")
    source_hashes = {
        "source_donor_scores_sha256": sha256(SCORES),
        "source_family_delta_sha256": sha256(FAMILY_DELTA),
        "source_class_table_sha256": sha256(CLASS_TABLE),
        "source_hca_table_sha256": sha256(HCA_TABLE),
        "source_fig3b_sha256": sha256(FIG3B_SOURCE),
        "source_authority_lock_sha256": sha256(AUTHORITY_LOCK),
    }

    require(scores["scorer"].eq("python_exact_reference_matching_singscore_1.28.1").all(), "Corrected scorer changed")
    require(scores["le_source"].eq("frozen_program_metadata").all(), "Frozen leading-edge source changed")
    class_lookup = dict(zip(classes["family_id"].astype(int), classes["fidelity_class"]))
    carrier_lookup = dict(zip(hca["family_id"].astype(int), hca["LA_carrier"].astype(str)))
    support_lookup = {
        int(row.family_id): set(str(row.sig_compartments).split(","))
        if pd.notna(row.sig_compartments) and str(row.sig_compartments)
        else set()
        for row in classes.itertuples(index=False)
    }
    for family_id in (10, 11):
        require(class_lookup[family_id] == EXPECTED_CLASS[family_id], f"F{family_id} class changed")
        require(carrier_lookup[family_id] == EXPECTED_CARRIER[family_id], f"F{family_id} carrier changed")
        require(support_lookup[family_id] == EXPECTED_SUPPORT[family_id], f"F{family_id} support changed")

    rows: list[pd.DataFrame] = []
    for target_order, (family_id, compartment, family_label) in enumerate(TARGETS, start=1):
        subset = scores.loc[
            scores["family_id"].astype(int).eq(family_id)
            & scores["compartment"].eq(compartment)
            & scores["score"].notna()
        ].copy()
        require(not subset.empty, f"Missing corrected donor scores for F{family_id} {compartment}")
        expected_programs = int(subset["program_id"].nunique())
        donor = (
            subset.groupby(["donor_id", "rhythm"], as_index=False)
            .agg(
                family_score=("score", "mean"),
                n_family_programs=("program_id", "nunique"),
            )
        )
        require(
            donor["n_family_programs"].eq(expected_programs).all(),
            f"F{family_id} {compartment} omits a constituent family program",
        )
        require(set(donor["rhythm"]) == {"Ctrl", "AF"}, "Unexpected donor group")

        delta_row = delta.loc[
            delta["family_id"].astype(int).eq(family_id)
            & delta["compartment"].eq(compartment)
        ]
        require(len(delta_row) == 1, f"Expected one family permutation row for F{family_id} {compartment}")
        drow = delta_row.iloc[0]
        fig3b_row = fig3b.loc[
            fig3b["family_id"].astype(int).eq(family_id)
            & fig3b["compartment"].eq(compartment)
        ]
        require(len(fig3b_row) == 1, f"Expected one Figure 3B row for F{family_id} {compartment}")
        brow = fig3b_row.iloc[0]

        counts = donor["rhythm"].value_counts().to_dict()
        means = donor.groupby("rhythm")["family_score"].mean().to_dict()
        observed_delta = float(means["AF"] - means["Ctrl"])
        require(counts == {"AF": int(drow["n_af"]), "Ctrl": int(drow["n_ctrl"])}, "Donor counts do not reproduce family summary")
        require(abs(float(means["Ctrl"]) - float(drow["mean_ctrl"])) <= 1e-10, "Control mean mismatch")
        require(abs(float(means["AF"]) - float(drow["mean_af"])) <= 1e-10, "AF mean mismatch")
        require(abs(observed_delta - float(drow["dScore"])) <= 1e-10, "Family delta mismatch")
        require(abs(float(drow["dScore"]) - float(brow["dScore"])) <= 1e-10, "Figure 3B delta mismatch")
        require(abs(float(drow["p_perm"]) - float(brow["p_perm"])) <= 1e-10, "Figure 3B p mismatch")
        require(int(drow["n_perm"]) == 10000 and str(drow["perm_strategy"]) == "mc", "Primary permutation definition changed")

        donor["family_id"] = family_id
        donor["family_label"] = family_label
        donor["fidelity_class"] = class_lookup[family_id]
        donor["compartment"] = compartment
        donor["compartment_label"] = {
            "fibroblast": "Fibroblast",
            "endothelial": "Endothelial",
            "adipocyte": "Adipocyte",
        }[compartment]
        donor["group"] = donor["rhythm"]
        donor["is_carrier"] = compartment == carrier_lookup[family_id]
        donor["is_supported"] = compartment in support_lookup[family_id]
        donor["n_ctrl"] = int(drow["n_ctrl"])
        donor["n_af"] = int(drow["n_af"])
        donor["mean_ctrl"] = float(drow["mean_ctrl"])
        donor["mean_af"] = float(drow["mean_af"])
        donor["delta"] = float(drow["dScore"])
        donor["p_perm"] = float(drow["p_perm"])
        donor["n_perm"] = int(drow["n_perm"])
        donor["perm_strategy"] = str(drow["perm_strategy"])
        for column, digest in source_hashes.items():
            donor[column] = digest
        donor["_target_order"] = target_order
        rows.append(donor)

    bridge = pd.concat(rows, ignore_index=True)
    bridge["_group_order"] = bridge["group"].map({"Ctrl": 0, "AF": 1})
    bridge["_donor_key"] = bridge["donor_id"].map(natural_key)
    bridge = bridge.sort_values(
        ["_target_order", "_group_order", "_donor_key"]
    ).reset_index(drop=True)
    bridge = bridge[OUTPUT_COLUMNS]

    require(len(bridge) == 99, f"Expected 99 family-level donor rows, got {len(bridge)}")
    require(
        bridge.groupby(["family_id", "compartment"])["donor_id"].nunique().to_dict()
        == {(10, "fibroblast"): 34, (11, "adipocyte"): 31, (11, "endothelial"): 34},
        "Facet donor totals changed",
    )
    f11_ad = bridge.loc[
        bridge["family_id"].eq(11) & bridge["compartment"].eq("adipocyte"),
        "p_perm",
    ].iloc[0]
    require(abs(float(f11_ad) - 0.0431) < 5e-5, "F11 adipocyte p_perm no longer reproduces 0.0431")
    return bridge


def main() -> None:
    bridge = build_bridge()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    bridge.to_csv(OUT, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {OUT.relative_to(ROOT)} ({len(bridge)} family-level donor rows)")
    print(f"SHA256 {sha256(OUT)}")


if __name__ == "__main__":
    main()
