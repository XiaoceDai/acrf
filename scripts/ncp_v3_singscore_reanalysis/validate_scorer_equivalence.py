#!/usr/bin/env python3
"""Validate an exact Python reference against the authoritative R scorer."""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
DATA = ROOT / "data" / "ncp_v3_singscore_reanalysis"
R_SCRIPT = (
    ROOT
    / "scripts"
    / "ncp_v3_singscore_reanalysis"
    / "validate_scorer_equivalence.R"
)
TOLERANCE = 1e-12


def exact_reference(matrix: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Match the installed singscore rankGenes/simpleScore for one up-set.

    rankGenes uses minimum ranks for ties. For k measured genes among N ranked
    genes, simpleScore(centerScore=TRUE) computes:
      (mean_rank - (k + 1) / 2) / (N - k) - 0.5
    """
    measured = [gene for gene in genes if gene in matrix.index]
    if len(measured) < 3:
        return pd.Series(np.nan, index=matrix.columns, dtype=float)
    ranks = matrix.rank(axis=0, method="min")
    k = len(measured)
    n = len(matrix)
    return (ranks.loc[measured].mean(axis=0) - (k + 1) / 2) / (n - k) - 0.5


def legacy_python_approximation(matrix: pd.DataFrame, genes: list[str]) -> pd.Series:
    measured = [gene for gene in genes if gene in matrix.index]
    ranks = matrix.rank(axis=0, method="average")
    return ranks.loc[measured].mean(axis=0) / len(matrix) - 0.5


def write_fixture() -> None:
    fixture = pd.DataFrame(
        {
            "sample_ties": [1, 1, 2, 4, 4, 5, 7, 8, 9, 9, 10, 12],
            "sample_reverse": [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
            "sample_mixed": [0, 3, 3, 2, 8, 5, 5, 1, 9, 7, 11, 10],
        },
        index=[f"G{i}" for i in range(1, 13)],
    )
    fixture.index.name = "gene_symbol"
    fixture.to_csv(RESULTS / "scorer_validation_fixture.tsv", sep="\t")


def load_matrix(dataset: str) -> pd.DataFrame:
    return pd.read_csv(
        DATA / "processed" / dataset / "expr_symbol.tsv", sep="\t", index_col=0
    )


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    write_fixture()
    subprocess.run(["Rscript", str(R_SCRIPT)], cwd=ROOT, check=True)

    programs = pd.read_csv(RESULTS / "frozen_program_metadata.tsv", sep="\t")
    matrices = {
        "fixture": pd.read_csv(
            RESULTS / "scorer_validation_fixture.tsv", sep="\t", index_col=0
        ),
        "GSE112339": load_matrix("GSE112339"),
    }
    gene_sets = {
        "fixture_ties_and_missing": ["G1", "G2", "G3", "G_missing"],
        **{
            row.program_id: [
                gene.strip()
                for gene in str(row.leading_edge_genes).split(",")
                if gene.strip()
            ]
            for row in programs.itertuples()
        },
    }

    r_scores = pd.read_csv(RESULTS / "r_scorer_validation_scores.tsv", sep="\t")
    rows = []
    for r in r_scores.itertuples():
        matrix = matrices[r.dataset]
        genes = gene_sets[r.program_id]
        exact = float(exact_reference(matrix, genes).loc[r.sample_id])
        approximate = float(
            legacy_python_approximation(matrix, genes).loc[r.sample_id]
        )
        rows.append(
            {
                "dataset": r.dataset,
                "sample_id": r.sample_id,
                "program_id": r.program_id,
                "n_measured_genes": r.n_measured_genes,
                "n_ranked_genes": r.n_ranked_genes,
                "r_singscore": r.r_singscore,
                "python_exact_reference": exact,
                "absolute_difference_exact": abs(exact - r.r_singscore),
                "legacy_python_approximation": approximate,
                "absolute_difference_legacy_approximation": abs(
                    approximate - r.r_singscore
                ),
            }
        )
    comparison = pd.DataFrame(rows)
    comparison.to_csv(
        RESULTS / "scorer_equivalence_comparison.tsv", sep="\t", index=False
    )
    summary = pd.DataFrame(
        [
            {
                "n_scores_compared": len(comparison),
                "tolerance": TOLERANCE,
                "max_abs_difference_exact_reference": comparison[
                    "absolute_difference_exact"
                ].max(),
                "exact_reference_gate_passed": bool(
                    (comparison["absolute_difference_exact"] <= TOLERANCE).all()
                ),
                "max_abs_difference_legacy_python_approximation": comparison[
                    "absolute_difference_legacy_approximation"
                ].max(),
                "legacy_python_is_equivalent": bool(
                    (
                        comparison["absolute_difference_legacy_approximation"]
                        <= TOLERANCE
                    ).all()
                ),
                "authoritative_scorer": (
                    "R singscore::rankGenes(tiesMethod=min) + "
                    "simpleScore(centerScore=TRUE)"
                ),
            }
        ]
    )
    summary.to_csv(
        RESULTS / "scorer_equivalence_summary.tsv", sep="\t", index=False
    )

    audit = """# Rank-based scorer audit

Authoritative implementation: `singscore::rankGenes` followed by
`singscore::simpleScore(upSet=measured_LE, centerScore=TRUE)`.

- Rank direction: ascending expression within each sample.
- Tie handling: minimum rank (`rankGenes` default `tiesMethod=\"min\"`).
- Normalization: for `k` measured set genes among `N` ranked genes,
  `(mean_rank - (k+1)/2)/(N-k) - 0.5`.
- Gene identifiers: exact HGNC-symbol intersection.
- Frozen set duplicates: removed before scoring.
- Expression duplicates: resolved upstream; maximum-variance count row for
  RNA-seq and maximum-IQR probe per symbol for GPL570 arrays.
- Measured intersection: only frozen LE symbols present in a cohort matrix.
- Minimum coverage: at least 3 measured LE genes, matching `p3_program_singscore.R`.
- Missing behavior: fewer than 3 genes or an unavailable chamber score yields
  NA; missing values are never converted to zero.

The prior v3 Python approximation uses average tie ranks and `mean(rank/N)-0.5`;
it is retained only in the comparison report and is not authoritative.
"""
    (RESULTS / "rank_based_scorer_audit.md").write_text(audit, encoding="utf-8")

    passed = bool(summary.loc[0, "exact_reference_gate_passed"])
    print(
        f"Exact R/Python scorer equivalence: {passed}; "
        f"max difference={summary.loc[0, 'max_abs_difference_exact_reference']:.3g}"
    )
    if not passed:
        raise SystemExit(42)


if __name__ == "__main__":
    main()
