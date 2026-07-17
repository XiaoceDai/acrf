#!/usr/bin/env python3
"""Build the dedicated, checksummed input resolver for corrected Figure 1B."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from _pathway_labels import short_pathway_label

RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
MANIFEST = RESULTS / "fig1b_input_manifest.tsv"
ORDER = RESULTS / "fig1b_canonical_display_order.tsv"

CONTRASTS = RESULTS / "healthy_program_contrasts_long.tsv"
REFERENCE = RESULTS / "healthy_reference_program_summary.tsv"
PROGRAMS = RESULTS / "frozen_program_metadata.tsv"
FAMILIES = RESULTS / "frozen_family_metadata.tsv"
LE_VALIDATION = RESULTS / "leading_edge_comparison_summary.tsv"
SCORER_VERSIONS = RESULTS / "singscore_software_versions.tsv"

SCORER_PACKAGE = "singscore"
SCORER_VERSION = "1.28.1"
SCORE_METHOD = (
    "singscore::rankGenes(tiesMethod=min)+"
    "simpleScore(upSet=measured_frozen_LE,centerScore=TRUE)"
)
REFERENCE_FORMULA = "mean(abs(direction_sign * delta))"
FAMILY_BLOCKS = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_order_table() -> pd.DataFrame:
    programs = pd.read_csv(PROGRAMS, sep="\t")
    families = pd.read_csv(FAMILIES, sep="\t")
    if len(programs) != 30 or programs["program_id"].nunique() != 30:
        raise AssertionError("Frozen program metadata is not exactly 30 programs")
    if families["family_id"].nunique() != 16:
        raise AssertionError("Frozen family metadata is not exactly 16 families")
    if set(programs["family_id"]) != set(FAMILY_BLOCKS):
        raise AssertionError("Frozen family IDs do not match the approved blocks")

    block_rank = {family_id: rank for rank, family_id in enumerate(FAMILY_BLOCKS)}
    order = programs.copy()
    order["family_block_order"] = order["family_id"].map(block_rank)
    order["publication_label"] = order["program_id"].map(
        lambda value: short_pathway_label(str(value), 42)
    )
    order = order.sort_values(
        ["family_block_order", "program_order", "program_id"]
    ).reset_index(drop=True)
    order["figure_display_order"] = range(1, len(order) + 1)
    order["approved_direction_block"] = order["family_id"].map(
        lambda value: "LA" if value in FAMILY_BLOCKS[:11] else "RA"
    )
    expected_direction = order["approved_direction_block"].map(
        {"LA": "LA-biased", "RA": "RA-biased"}
    )
    if not (order["canonical_direction"] == expected_direction).all():
        raise AssertionError("Canonical directions disagree with approved family blocks")
    if order["canonical_display_order"].tolist() != list(range(1, 31)):
        raise AssertionError("Frozen canonical display order is not contiguous 1..30")
    if order["program_id"].tolist() != programs.sort_values(
        "canonical_display_order"
    )["program_id"].tolist():
        raise AssertionError(
            "Approved family-block order differs from frozen canonical program order"
        )
    columns = [
        "figure_display_order",
        "canonical_display_order",
        "family_block_order",
        "approved_direction_block",
        "family_id",
        "program_order",
        "program_id",
        "publication_label",
        "canonical_direction",
        "direction_sign",
        "gene_set_version",
    ]
    order[columns].to_csv(ORDER, sep="\t", index=False)
    return order[columns]


def main() -> None:
    required = [
        CONTRASTS,
        REFERENCE,
        PROGRAMS,
        FAMILIES,
        LE_VALIDATION,
        SCORER_VERSIONS,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing authoritative Figure 1B inputs: {missing}")

    leading_edge = pd.read_csv(LE_VALIDATION, sep="\t")
    if not (
        int(leading_edge.loc[0, "n_programs"]) == 30
        and int(leading_edge.loc[0, "n_set_identical"]) == 30
        and bool(leading_edge.loc[0, "gate_passed"])
    ):
        raise AssertionError("Leading-edge validation is not an exact 30/30 match")
    versions = pd.read_csv(SCORER_VERSIONS, sep="\t")
    observed_version = versions.loc[
        versions["component"] == SCORER_PACKAGE, "version"
    ]
    if observed_version.tolist() != [SCORER_VERSION]:
        raise AssertionError(
            f"Expected {SCORER_PACKAGE} {SCORER_VERSION}, "
            f"observed {observed_version.tolist()}"
        )

    order = build_order_table()
    leading_edge_sha = order["gene_set_version"].drop_duplicates().tolist()
    if len(leading_edge_sha) != 1:
        raise AssertionError("Programs do not share one frozen leading-edge version")
    leading_edge_version = (
        f"{leading_edge_sha[0]}; exact_match_selection_time_fgsea_30_of_30"
    )
    artifacts = [
        (
            "healthy_donor_program_contrasts",
            CONTRASTS,
            "scripts/ncp_v3_singscore_reanalysis/compute_acairel.py",
            "GSE112339",
        ),
        (
            "healthy_program_reference_summary",
            REFERENCE,
            "scripts/ncp_v3_singscore_reanalysis/compute_acairel.py",
            "GSE112339",
        ),
        (
            "frozen_program_metadata",
            PROGRAMS,
            "scripts/ncp_v3_singscore_reanalysis/freeze_inputs.py",
            "frozen_NCP_v3",
        ),
        (
            "frozen_family_metadata",
            FAMILIES,
            "scripts/ncp_v3_singscore_reanalysis/freeze_inputs.py",
            "frozen_NCP_v3",
        ),
        (
            "canonical_display_order",
            ORDER,
            "scripts/ncp_v3_singscore_reanalysis/build_fig1b_input_manifest.py",
            "frozen_NCP_v3",
        ),
        (
            "leading_edge_validation",
            LE_VALIDATION,
            "scripts/ncp_v3_singscore_reanalysis/compare_leading_edges.R",
            "GSE112339_selection_rank",
        ),
        (
            "scorer_version_record",
            SCORER_VERSIONS,
            "scripts/ncp_v3_singscore_reanalysis/validate_scorer_equivalence.R",
            "software_environment",
        ),
    ]
    rows = []
    namespace = (RESULTS.resolve(),)
    for role, path, generated_by, dataset in artifacts:
        resolved = path.resolve()
        if not any(resolved.is_relative_to(base) for base in namespace):
            raise AssertionError(f"Input outside corrected namespace: {path}")
        rows.append(
            {
                "artifact_role": role,
                "artifact_path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "generated_by": generated_by,
                "dataset": dataset,
                "score_method": SCORE_METHOD,
                "scorer_package": SCORER_PACKAGE,
                "scorer_version": SCORER_VERSION,
                "leading_edge_version": leading_edge_version,
                "reference_formula": REFERENCE_FORMULA,
                "is_authoritative": True,
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(MANIFEST, sep="\t", index=False)
    print(f"Wrote {MANIFEST.relative_to(ROOT)} ({len(manifest)} artifacts)")
    print(f"Wrote {ORDER.relative_to(ROOT)} ({len(order)} programs)")


if __name__ == "__main__":
    main()
