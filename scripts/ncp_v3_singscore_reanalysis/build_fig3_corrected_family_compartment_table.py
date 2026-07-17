#!/usr/bin/env python3
"""Freeze the corrected family-compartment permutation/classification contract."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
FIG3B = (
    RESULTS
    / "ncp_v3_singscore_reanalysis"
    / "fig3b_family_compartment_localization_source.tsv"
)
DELTA = RESULTS / "fig3_v3_family_celltype_delta.tsv"
SENSITIVITY = RESULTS / "fig3_boundary_permutation_sensitivity.tsv"
OUT = (
    RESULTS
    / "ncp_v3_singscore_reanalysis"
    / "gse255612_corrected_family_compartment_permutation_classification.tsv"
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


def main() -> None:
    require(
        FIG3B.is_file() and DELTA.is_file() and SENSITIVITY.is_file(),
        "A corrected family-level contract input is missing",
    )
    fig3b = pd.read_csv(FIG3B, sep="\t")
    delta = pd.read_csv(DELTA, sep="\t")
    sensitivity = pd.read_csv(SENSITIVITY, sep="\t")
    delta_contract = delta[
        [
            "family_id",
            "compartment",
            "mean_ctrl",
            "mean_af",
            "n_perm",
            "perm_strategy",
        ]
    ]
    table = fig3b.merge(
        delta_contract,
        on=["family_id", "compartment"],
        how="left",
        validate="one_to_one",
    )
    require(table["n_perm"].eq(10000).all(), "Primary permutation count changed")
    require(table["perm_strategy"].eq("mc").all(), "Permutation strategy changed")
    table = table.merge(
        sensitivity[
            [
                "family_id",
                "compartment",
                "p_perm_sensitivity_1e6",
                "n_perm_sensitivity",
                "note",
            ]
        ],
        on=["family_id", "compartment"],
        how="left",
        validate="one_to_one",
    )
    table["p_perm_sensitivity_1e6"] = table[
        "p_perm_sensitivity_1e6"
    ].astype(float)
    table["n_perm_sensitivity"] = table["n_perm_sensitivity"].astype("Int64")
    table["source_fig3b_sha256"] = sha256(FIG3B)
    table["source_family_delta_sha256"] = sha256(DELTA)
    table["source_sensitivity_sha256"] = sha256(SENSITIVITY)
    table = table.rename(
        columns={
            "canonical_direction": "direction",
            "is_carrier_cell": "is_carrier",
            "is_sig_p_lt_0_05": "is_supported",
            "hca_LA_carrier": "healthy_LA_carrier",
            "note": "sensitivity_note",
        }
    )
    f15_endo = table.loc[
        table["family_id"].eq(15) & table["compartment"].eq("endothelial")
    ].iloc[0]
    f16_adipo = table.loc[
        table["family_id"].eq(16) & table["compartment"].eq("adipocyte")
    ].iloc[0]
    require(abs(float(f15_endo["p_perm"]) - 0.0471) < 5e-5, "F15 Endothelial p changed")
    require(
        int(f15_endo["n_perm_sensitivity"]) == 1_000_000
        and np.isfinite(float(f15_endo["p_perm_sensitivity_1e6"])),
        "F15 Endothelial sensitivity is missing",
    )
    require(abs(float(f16_adipo["p_perm"]) - 0.0542) < 5e-5, "F16 Adipocyte p changed")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {OUT.relative_to(ROOT)} ({len(table)} family-compartment rows)")
    print(f"SHA256 {sha256(OUT)}")


if __name__ == "__main__":
    main()
