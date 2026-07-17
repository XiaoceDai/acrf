#!/usr/bin/env python3
"""Build the checksummed corrected-input resolver for final Figure 2A."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
DATA = ROOT / "data" / "ncp_v3_singscore_reanalysis"
MANIFEST = RESULTS / "fig2a_input_manifest.tsv"

PROGRAMS = RESULTS / "frozen_program_metadata.tsv"
FAMILIES = RESULTS / "frozen_family_metadata.tsv"
LE_VALIDATION = RESULTS / "leading_edge_comparison_summary.tsv"
SCORER_VERSIONS = RESULTS / "singscore_software_versions.tsv"
REFERENCE = RESULTS / "healthy_reference_program_summary.tsv"
DISEASE_CONTRASTS = RESULTS / "disease_program_contrasts_long.tsv"
FAMILY_VALUES = RESULTS / "family_acairel_long.tsv"
PAIRING_SUMMARY = RESULTS / "cohort_pairing_summary.tsv"
METADATA = {
    cohort: DATA / "metadata" / f"{cohort}_sample_metadata.tsv"
    for cohort in ("GSE128188", "GSE79768", "GSE115574")
}

SCORER_PACKAGE = "singscore"
SCORER_VERSION = "1.28.1"
SCORE_METHOD = (
    "singscore::rankGenes(tiesMethod=min)+"
    "simpleScore(upSet=measured_frozen_LE,centerScore=TRUE)"
)
REFERENCE_FORMULA = "mean(abs(direction_sign * delta))"
EXPECTED_COUNTS = {
    "GSE128188": {"SR": 5, "AF": 5},
    "GSE79768": {"SR": 6, "AF": 7},
    "GSE115574": {"SR": 15, "AF": 13},
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def included_metadata(path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(path, sep="\t")
    required = {
        "dataset",
        "sample_id",
        "donor_id",
        "chamber",
        "rhythm",
        "analysis_inclusion",
    }
    require(required.issubset(metadata.columns), f"Metadata columns missing from {path}")
    included = metadata.loc[
        metadata["analysis_inclusion"].astype(str).str.lower().eq("true")
    ].copy()
    return included


def validate_sources() -> str:
    required_paths = [
        PROGRAMS,
        FAMILIES,
        LE_VALIDATION,
        SCORER_VERSIONS,
        REFERENCE,
        DISEASE_CONTRASTS,
        FAMILY_VALUES,
        PAIRING_SUMMARY,
        *METADATA.values(),
    ]
    missing = [str(path) for path in required_paths if not path.is_file()]
    require(not missing, f"Missing corrected Figure 2A inputs: {missing}")

    programs = pd.read_csv(PROGRAMS, sep="\t")
    families = pd.read_csv(FAMILIES, sep="\t")
    require(len(programs) == 30 and programs["program_id"].nunique() == 30, "Not 30 frozen programs")
    require(families["family_id"].nunique() == 16 and len(families) == 16, "Not 16 frozen families")
    require(programs["family_id"].nunique() == 16, "Frozen programs do not span 16 families")
    versions = programs["gene_set_version"].drop_duplicates().tolist()
    require(len(versions) == 1, "Frozen programs do not share one leading-edge version")

    leading_edge = pd.read_csv(LE_VALIDATION, sep="\t").iloc[0]
    require(
        int(leading_edge["n_programs"]) == 30
        and int(leading_edge["n_set_identical"]) == 30
        and str(leading_edge["gate_passed"]).lower() == "true",
        "Leading-edge validation is not an exact 30/30 match",
    )
    software = pd.read_csv(SCORER_VERSIONS, sep="\t")
    observed = software.loc[software["component"] == SCORER_PACKAGE, "version"].astype(str).tolist()
    require(observed == [SCORER_VERSION], f"Expected singscore {SCORER_VERSION}, observed {observed}")

    reference = pd.read_csv(REFERENCE, sep="\t")
    require(len(reference) == 30 and reference["program_id"].nunique() == 30, "Reference M_p table is not 30 programs")
    require(set(reference["reference_dataset"].astype(str)) == {"GSE112339"}, "M_p is not exclusively GSE112339")

    contrasts = pd.read_csv(DISEASE_CONTRASTS, sep="\t")
    require(len(contrasts) == 51 * 30, "Disease contrast table is not 51 donors x 30 programs")
    require(set(contrasts["dataset"]) == set(EXPECTED_COUNTS), "Unexpected disease cohort in contrasts")
    require(set(contrasts["reference_dataset"].astype(str)) == {"GSE112339"}, "Disease aCAIrel uses a non-healthy denominator")
    require(set(contrasts["score_method"].astype(str)) == {SCORE_METHOD}, "Disease contrasts use an unexpected scorer")
    require(contrasts["program_id"].nunique() == 30, "Disease contrasts do not contain 30 programs")

    for cohort, expected in EXPECTED_COUNTS.items():
        metadata = included_metadata(METADATA[cohort])
        require(set(metadata["dataset"].astype(str)) == {cohort}, f"Dataset mismatch in {cohort} metadata")
        donor_chambers = metadata.groupby("donor_id")["chamber"].agg(lambda values: set(values))
        require(donor_chambers.map(lambda value: value == {"LA", "RA"}).all(), f"Incomplete included pair in {cohort}")
        donor_rhythm = metadata.groupby("donor_id")["rhythm"].nunique()
        require((donor_rhythm == 1).all(), f"Donor rhythm mismatch in {cohort}")
        counts = metadata.drop_duplicates("donor_id")["rhythm"].value_counts().to_dict()
        require(counts == expected, f"{cohort} donor counts differ: {counts} != {expected}")

    return versions[0]


def main() -> None:
    leading_edge_version = validate_sources()
    artifacts = [
        ("frozen_program_metadata", PROGRAMS, "frozen_topology", "all"),
        ("frozen_family_metadata", FAMILIES, "frozen_topology", "all"),
        ("leading_edge_validation", LE_VALIDATION, "frozen_leading_edges", "GSE112339_selection"),
        ("scorer_version_record", SCORER_VERSIONS, "software_version", "all"),
        ("healthy_program_reference_summary", REFERENCE, "healthy_M_p_only", "GSE112339"),
        ("disease_program_contrasts", DISEASE_CONTRASTS, "program_aCAIrel", "three_AF_cohorts"),
        ("family_acairel", FAMILY_VALUES, "family_aCAIrel_validation", "healthy_and_three_AF_cohorts"),
        ("cohort_pairing_summary", PAIRING_SUMMARY, "pairing_audit", "all"),
    ]
    artifacts.extend(
        (f"{cohort}_sample_metadata", path, "paired_donor_counts", cohort)
        for cohort, path in METADATA.items()
    )

    allowed_roots = (RESULTS.resolve(), DATA.resolve())
    rows = []
    for role, path, used_for, dataset in artifacts:
        resolved = path.resolve()
        require(
            any(resolved.is_relative_to(base) for base in allowed_roots),
            f"Input outside corrected namespaces: {path}",
        )
        rows.append(
            {
                "artifact_role": role,
                "artifact_path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "generated_by": {
                    "frozen_program_metadata": "scripts/ncp_v3_singscore_reanalysis/freeze_inputs.py",
                    "frozen_family_metadata": "scripts/ncp_v3_singscore_reanalysis/freeze_inputs.py",
                    "disease_program_contrasts": "scripts/ncp_v3_singscore_reanalysis/compute_acairel.py",
                    "family_acairel": "scripts/ncp_v3_singscore_reanalysis/compute_acairel.py",
                }.get(role, "corrected_ncp_v3_singscore_reanalysis"),
                "dataset": dataset,
                "used_for": used_for,
                "score_method": SCORE_METHOD,
                "scorer_package": SCORER_PACKAGE,
                "scorer_version": SCORER_VERSION,
                "leading_edge_version": leading_edge_version,
                "reference_dataset": "GSE112339",
                "reference_formula": REFERENCE_FORMULA,
                "allows_cohort_specific_reselection": False,
                "allows_reference_recalibration": False,
                "is_authoritative": True,
            }
        )

    manifest = pd.DataFrame(rows)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(MANIFEST, sep="\t", index=False)
    print(f"Wrote {MANIFEST.relative_to(ROOT)} ({len(manifest)} corrected artifacts)")


if __name__ == "__main__":
    main()
