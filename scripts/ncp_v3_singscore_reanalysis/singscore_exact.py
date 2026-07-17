#!/usr/bin/env python3
"""Shared exact singscore reference matching Bioconductor singscore 1.28.1.

Authoritative formula (validated against R):
  rankGenes(tiesMethod=\"min\") + simpleScore(upSet=measured, centerScore=TRUE)
  => (mean_rank - (k + 1) / 2) / (N - k) - 0.5

Production bulk scoring remains R (`score_programs.R`). This module is the
validated Python equivalent for pipelines that must score in Python (e.g. Fig 3).
"""
from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd

MIN_MEASURED_GENES = 3
TOLERANCE = 1e-12


def normalize_gene_list(genes: Iterable[str]) -> list[str]:
    """Deduplicate and trim LE symbols (order-preserving)."""
    seen: set[str] = set()
    out: list[str] = []
    for gene in genes:
        symbol = str(gene).strip()
        if not symbol or symbol in seen:
            continue
        seen.add(symbol)
        out.append(symbol)
    return out


def parse_leading_edge_field(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    return normalize_gene_list(str(value).split(","))


def exact_reference(matrix: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Match singscore::rankGenes(min) + simpleScore(centerScore=TRUE).

    matrix : genes × samples (index = gene symbols, columns = sample IDs)
    genes  : frozen leading-edge symbols (already normalized)
    """
    measured = [gene for gene in genes if gene in matrix.index]
    if len(measured) < MIN_MEASURED_GENES:
        return pd.Series(np.nan, index=matrix.columns, dtype=float)
    ranks = matrix.rank(axis=0, method="min")
    k = len(measured)
    n = len(matrix)
    return (ranks.loc[measured].mean(axis=0) - (k + 1) / 2) / (n - k) - 0.5


def legacy_v3_rank_singscore(matrix: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Prior Fig 3 Python approximation (average ties; mean(rank/N)-0.5).

    Retained only for equivalence audits. Not authoritative.
    """
    measured = [gene for gene in genes if gene in matrix.index]
    if len(measured) < 5:
        return pd.Series(np.nan, index=matrix.columns, dtype=float)
    ranks = matrix.rank(axis=0, method="average")
    return ranks.loc[measured].mean(axis=0) / len(matrix) - 0.5


def scores_equivalent(
    a: pd.Series,
    b: pd.Series,
    *,
    tolerance: float = TOLERANCE,
) -> tuple[bool, dict[str, float | int | bool]]:
    """Strict equivalence: missingness, eligibility, signs, max |Δ| <= tolerance."""
    a = a.reindex(b.index)
    a_na = a.isna()
    b_na = b.isna()
    missingness_identical = bool((a_na == b_na).all())
    both = ~(a_na | b_na)
    n_both = int(both.sum())
    if n_both == 0:
        signs_identical = True
        max_abs = 0.0
    else:
        aa = a.loc[both].to_numpy(dtype=float)
        bb = b.loc[both].to_numpy(dtype=float)
        signs_identical = bool(np.all(np.sign(aa) == np.sign(bb)))
        max_abs = float(np.max(np.abs(aa - bb)))
    passed = (
        missingness_identical
        and signs_identical
        and max_abs <= tolerance
    )
    return passed, {
        "n_compared": int(len(a)),
        "n_finite_both": n_both,
        "missingness_identical": missingness_identical,
        "signs_identical": signs_identical,
        "max_abs_difference": max_abs,
        "tolerance": tolerance,
        "passed": passed,
    }
