#!/usr/bin/env python3
"""Validate final gates and checksum all namespaced reanalysis outputs."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
MANIFEST = RESULTS / "analysis_output_manifest.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    assertions = pd.read_csv(RESULTS / "assertion_report.tsv", sep="\t")
    leading_edge = pd.read_csv(
        RESULTS / "leading_edge_comparison_summary.tsv", sep="\t"
    )
    scorer = pd.read_csv(RESULTS / "scorer_equivalence_summary.tsv", sep="\t")
    lineage = pd.read_csv(RESULTS / "corrected_lineage_audit.tsv", sep="\t")
    if not assertions["passed"].all():
        raise AssertionError("One or more scientific assertions failed")
    if not bool(leading_edge.loc[0, "gate_passed"]):
        raise AssertionError("Leading-edge gate failed")
    if not bool(scorer.loc[0, "exact_reference_gate_passed"]):
        raise AssertionError("Scorer equivalence gate failed")
    if not lineage["passed"].all():
        raise AssertionError("Corrected/legacy lineage gate failed")

    figure_summary = pd.read_csv(
        RESULTS / "fig1b_healthy_reference_concordance_summary.tsv", sep="\t"
    ).iloc[0]
    program_scores = pd.read_csv(RESULTS / "program_scores_long.tsv", sep="\t")
    changed = pd.read_csv(
        RESULTS / "migration_changed_classifications.tsv", sep="\t"
    )
    path_a = pd.read_csv(RESULTS / "fig5c_v3_cross_axis.tsv", sep="\t")
    summary = pd.DataFrame(
        [
            {
                "n_programs": int(figure_summary["n_programs"]),
                "n_families": int(figure_summary["n_families"]),
                "n_sample_program_scores": len(program_scores),
                "n_missing_program_scores": int(
                    program_scores["program_score"].isna().sum()
                ),
                "n_full_8_of_8_concordance": int(
                    figure_summary["n_full_concordance"]
                ),
                "n_ge_6_of_8_concordance": int(
                    figure_summary["n_ge6_concordance"]
                ),
                "median_M": float(figure_summary["median_M"]),
                "n_changed_classification_rows": len(changed),
                "n_Path_A_triple_positive": int(
                    (path_a["cross_axis_class"] == "Triple-positive").sum()
                ),
                "all_gates_passed": True,
            }
        ]
    )
    summary.to_csv(
        RESULTS / "reanalysis_completion_summary.tsv", sep="\t", index=False
    )

    paths = [
        path
        for directory in (RESULTS, FIGURES)
        for path in directory.iterdir()
        if path.is_file() and path != MANIFEST
    ]
    rows = [
        {
            "namespace": "results" if path.parent == RESULTS else "figures",
            "path": str(path.relative_to(ROOT)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(paths)
    ]
    pd.DataFrame(rows).to_csv(MANIFEST, sep="\t", index=False)
    print(f"All gates passed; checksummed {len(rows)} namespaced outputs")


if __name__ == "__main__":
    main()
