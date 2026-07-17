#!/usr/bin/env python3
"""Freeze internal inputs and compare reconstructed matrices with legacy copies."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
DATA = ROOT / "data" / "ncp_v3_singscore_reanalysis"
PROCESSED = DATA / "processed"
METADATA = DATA / "metadata"
LEGACY = ROOT / "Project_P2" / "Project_P1"

PANEL = ROOT / "results" / "innovation_ncp_v3_panel.tsv"
LEADING_EDGE = ROOT / "results" / "fig3_v3_program_le.tsv"
FAMILY_MAP = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
FAMILY_LABELS = ROOT / "results" / "ncp_v3_family_labels.tsv"
GMT_HALLMARK = ROOT / "bulk" / "Project_P01" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REACTOME = ROOT / "bulk" / "Project_P01" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def freeze_program_tables() -> None:
    panel = pd.read_csv(PANEL, sep="\t")
    le = pd.read_csv(LEADING_EDGE, sep="\t")
    families = pd.read_csv(FAMILY_MAP, sep="\t")
    labels = pd.read_csv(FAMILY_LABELS, sep="\t")
    panel = panel.rename(columns={"pathway": "program_id"})
    families = families.rename(columns={"pathway": "program_id"})

    if len(panel) != 30 or panel["program_id"].nunique() != 30:
        raise AssertionError("Deployed panel is not exactly 30 unique programs")
    merged = panel[["program_id", "program_order", "direction"]].merge(
        families[["program_id", "family_id", "family_source"]],
        on="program_id",
        validate="one_to_one",
    )
    merged = merged.merge(
        le[["program_id", "n_LE_total", "leading_edge_genes"]],
        on="program_id",
        validate="one_to_one",
    )
    if merged["family_id"].nunique() != 16:
        raise AssertionError("Frozen mapping is not exactly 16 families")
    counts = merged["direction"].value_counts().to_dict()
    if counts != {"LA-biased": 22, "RA-biased": 8}:
        raise AssertionError(f"Unexpected direction counts: {counts}")
    merged["direction_sign"] = merged["direction"].map(
        {"LA-biased": 1, "RA-biased": -1}
    )
    merged["canonical_direction"] = merged["direction"]
    merged["direction_rank"] = merged["direction"].map(
        {"LA-biased": 0, "RA-biased": 1}
    )
    merged = merged.sort_values(
        ["direction_rank", "family_id", "program_order"]
    ).reset_index(drop=True)
    merged["canonical_display_order"] = np.arange(1, len(merged) + 1)
    merged["gene_set_version"] = f"sha256:{sha256(LEADING_EDGE)}"
    merged.to_csv(RESULTS / "frozen_program_metadata.tsv", sep="\t", index=False)

    family_direction = (
        merged.groupby("family_id")["direction"].agg(lambda x: sorted(set(x)))
    )
    mixed = family_direction[family_direction.map(len) != 1]
    if not mixed.empty:
        raise AssertionError(f"Families with mixed canonical directions: {mixed.index.tolist()}")
    fam = labels.drop_duplicates("family_id").copy()
    fam["canonical_direction"] = fam["family_id"].map(
        family_direction.map(lambda x: x[0])
    )
    fam["direction_rank"] = fam["canonical_direction"].map(
        {"LA-biased": 0, "RA-biased": 1}
    )
    fam = fam.sort_values(["direction_rank", "family_id"]).reset_index(drop=True)
    fam["canonical_family_display_order"] = np.arange(1, len(fam) + 1)
    if len(fam) != 16 or fam["family_id"].nunique() != 16:
        raise AssertionError("Family label table does not resolve to 16 families")
    fam.to_csv(RESULTS / "frozen_family_metadata.tsv", sep="\t", index=False)


def load_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    return matrix


def legacy_matrix(dataset: str) -> pd.DataFrame | None:
    if dataset == "GSE112339":
        path = ROOT / "bulk" / "logcpm_voom.tsv"
        frame = pd.read_csv(path, sep="\t")
        if not {"gene_key", "gene_symbol"}.issubset(frame.columns):
            raise ValueError(f"Unexpected healthy legacy matrix schema: {path}")
        sample_columns = [c for c in frame.columns if c not in {"gene_key", "gene_symbol"}]
        metadata = pd.read_csv(
            METADATA / "GSE112339_sample_metadata.tsv", sep="\t"
        )
        rename = dict(
            zip(metadata["source_expression_column"], metadata["sample_id"])
        )
        frame = frame.set_index("gene_symbol")[sample_columns].rename(columns=rename)
        frame = frame[~frame.index.duplicated(keep="first")]
        return frame
    path = LEGACY / dataset / "expr_symbol.tsv"
    return load_matrix(path) if path.is_file() else None


def compare_matrices() -> None:
    rows = []
    for dataset in ("GSE112339", "GSE128188", "GSE79768", "GSE115574"):
        corrected = load_matrix(PROCESSED / dataset / "expr_symbol.tsv")
        old = legacy_matrix(dataset)
        if old is None:
            rows.append(
                {
                    "dataset": dataset,
                    "legacy_matrix_found": False,
                    "n_corrected_genes": len(corrected),
                    "n_legacy_genes": np.nan,
                    "n_shared_genes": 0,
                    "n_shared_samples": 0,
                    "median_sample_spearman": np.nan,
                    "minimum_sample_spearman": np.nan,
                    "maximum_absolute_difference": np.nan,
                    "matrix_identical_within_1e-10": False,
                }
            )
            continue
        shared_genes = corrected.index.intersection(old.index)
        shared_samples = corrected.columns.intersection(old.columns)
        correlations = []
        max_difference = 0.0
        for sample in shared_samples:
            x = corrected.loc[shared_genes, sample].to_numpy(dtype=float)
            y = old.loc[shared_genes, sample].to_numpy(dtype=float)
            finite = np.isfinite(x) & np.isfinite(y)
            correlations.append(float(spearmanr(x[finite], y[finite]).statistic))
            max_difference = max(
                max_difference, float(np.max(np.abs(x[finite] - y[finite])))
            )
        rows.append(
            {
                "dataset": dataset,
                "legacy_matrix_found": True,
                "n_corrected_genes": len(corrected),
                "n_legacy_genes": len(old),
                "n_shared_genes": len(shared_genes),
                "n_shared_samples": len(shared_samples),
                "median_sample_spearman": float(np.median(correlations)),
                "minimum_sample_spearman": float(np.min(correlations)),
                "maximum_absolute_difference": max_difference,
                "matrix_identical_within_1e-10": (
                    len(corrected) == len(old)
                    and len(corrected.columns) == len(old.columns)
                    and max_difference <= 1e-10
                ),
            }
        )
    pd.DataFrame(rows).to_csv(
        RESULTS / "public_vs_legacy_matrix_comparison.tsv", sep="\t", index=False
    )


def write_manifest() -> None:
    paths = [
        PANEL,
        LEADING_EDGE,
        FAMILY_MAP,
        FAMILY_LABELS,
        GMT_HALLMARK,
        GMT_REACTOME,
        RESULTS / "public_input_manifest.tsv",
        METADATA / "unified_sample_metadata.tsv",
        METADATA / "gene_identifier_mapping.tsv",
    ]
    paths.extend(
        METADATA / f"{dataset}_sample_metadata.tsv"
        for dataset in ("GSE112339", "GSE128188", "GSE79768", "GSE115574")
    )
    paths.extend(
        PROCESSED / dataset / "expr_symbol.tsv"
        for dataset in ("GSE112339", "GSE128188", "GSE79768", "GSE115574")
    )
    rows = []
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        rows.append(
            {
                "kind": (
                    "frozen_internal_input"
                    if str(path).startswith(str(ROOT / "results"))
                    or str(path).startswith(str(ROOT / "bulk"))
                    else "reconstructed_input"
                ),
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    pd.DataFrame(rows).to_csv(
        RESULTS / "analysis_input_manifest.tsv", sep="\t", index=False
    )


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    freeze_program_tables()
    compare_matrices()
    write_manifest()
    print("Wrote frozen program/family tables, matrix comparison, and input manifest")


if __name__ == "__main__":
    main()
