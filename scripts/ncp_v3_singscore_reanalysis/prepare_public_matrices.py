#!/usr/bin/env python3
"""Reconstruct frozen sample metadata and processed array expression matrices."""
from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "ncp_v3_singscore_reanalysis"
META_DIR = DATA / "metadata"
PROCESSED = DATA / "processed"
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"

SERIES = {
    "GSE112339": META_DIR / "GSE112339_series_matrix.txt.gz",
    "GSE128188": META_DIR / "GSE128188_series_matrix.txt.gz",
    "GSE79768": META_DIR / "GSE79768_series_matrix.txt.gz",
    "GSE115574": META_DIR / "GSE115574_series_matrix.txt.gz",
}
EXPECTED_SAMPLES = {
    "GSE112339": 40,
    "GSE128188": 20,
    "GSE79768": 26,
    "GSE115574": 59,
}
EXPECTED_PAIRS = {
    "GSE112339": 8,
    "GSE128188": 10,
    "GSE79768": 13,
    "GSE115574": 28,
}


def clean_field(value: str | None) -> str:
    return "" if value is None else value.strip().strip('"')


def read_series(path: Path) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    metadata_rows: dict[str, list[list[str]]] = {}
    table_header: list[str] | None = None
    table_rows: list[list[str]] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                fields = [clean_field(v) for v in line.rstrip("\n").split("\t")]
                if table_header is None:
                    table_header = fields
                else:
                    table_rows.append(fields)
                continue
            if line.startswith("!Sample_"):
                fields = [clean_field(v) for v in line.rstrip("\n").split("\t")]
                key = fields[0][len("!Sample_") :]
                metadata_rows.setdefault(key, []).append(fields[1:])

    if "geo_accession" not in metadata_rows:
        raise ValueError(f"No GEO accessions found in {path}")
    gsms = metadata_rows["geo_accession"][0]
    n = len(gsms)
    titles = metadata_rows.get("title", [[""] * n])[0]
    sources = metadata_rows.get("source_name_ch1", [[""] * n])[0]
    characteristics = metadata_rows.get("characteristics_ch1", [])
    meta = pd.DataFrame(
        {
            "sample_id": gsms,
            "title": titles,
            "source_name": sources,
            "characteristics": [
                " | ".join(row[i] for row in characteristics) for i in range(n)
            ],
        }
    )
    matrix = None
    if table_header is not None and table_rows:
        matrix = pd.DataFrame(table_rows, columns=table_header)
        id_column = matrix.columns[0]
        matrix = matrix.set_index(id_column)
        matrix = matrix.apply(pd.to_numeric, errors="coerce")
        matrix.index = matrix.index.astype(str)
    return meta, matrix


def classify_sample(accession: str, title: str) -> dict[str, str]:
    title = title.strip()
    if accession == "GSE112339":
        base = re.sub(r"\s*\[.*$", "", title)
        donor_raw, tissue = base.rsplit("_", 1)
        chamber = {"SEPTUM": "IVS"}.get(tissue.upper(), tissue.upper())
        return {
            "donor_id": donor_raw.replace("-", "_"),
            "chamber": chamber,
            "rhythm": "Healthy",
            "source_expression_column": (
                f"{chamber}_{donor_raw.replace('-', '_')}"
            ),
        }
    if accession == "GSE128188":
        match = re.fullmatch(r"(SR|AF)_(LA|RA)_(\S+)", title)
        if not match:
            raise ValueError(f"Unparsed GSE128188 title: {title}")
        rhythm, chamber, donor = match.groups()
        return {
            "donor_id": donor,
            "chamber": chamber,
            "rhythm": rhythm,
            "source_expression_column": f"{donor}-{chamber}",
        }
    if accession == "GSE79768":
        match = re.fullmatch(
            r"(Atrial Fibrillation|Sinus Rhythm)_(left atrial|right atrial)_(\S+)",
            title,
        )
        if not match:
            raise ValueError(f"Unparsed GSE79768 title: {title}")
        rhythm_raw, chamber_raw, donor = match.groups()
        return {
            "donor_id": donor,
            "chamber": {"left atrial": "LA", "right atrial": "RA"}[chamber_raw],
            "rhythm": {
                "Atrial Fibrillation": "AF",
                "Sinus Rhythm": "SR",
            }[rhythm_raw],
            "source_expression_column": "",
        }
    if accession == "GSE115574":
        match = re.fullmatch(r"(AF|SR)(\d{2})([LR])", title)
        if not match:
            raise ValueError(f"Unparsed GSE115574 title: {title}")
        rhythm, donor_number, side = match.groups()
        return {
            "donor_id": f"{rhythm}{donor_number}",
            "chamber": {"L": "LA", "R": "RA"}[side],
            "rhythm": rhythm,
            "source_expression_column": "",
        }
    raise ValueError(accession)


def build_metadata(accession: str, raw: pd.DataFrame) -> pd.DataFrame:
    classified = pd.DataFrame(
        [classify_sample(accession, title) for title in raw["title"]]
    )
    meta = pd.concat([raw.reset_index(drop=True), classified], axis=1)
    meta.insert(0, "dataset", accession)
    meta["platform"] = (
        "RNA-seq" if accession in {"GSE112339", "GSE128188"} else "GPL570 microarray"
    )
    meta["analysis_inclusion"] = False
    meta["exclusion_reason"] = ""

    if accession == "GSE112339":
        eligible = meta["chamber"].isin(["LA", "RA"])
        meta.loc[~eligible, "exclusion_reason"] = "non_atrial_chamber"
    else:
        eligible = pd.Series(True, index=meta.index)

    for donor_id, group in meta[eligible].groupby("donor_id"):
        complete = (
            len(group) == 2
            and set(group["chamber"]) == {"LA", "RA"}
            and group["sample_id"].nunique() == 2
        )
        if complete:
            meta.loc[group.index, "analysis_inclusion"] = True
        else:
            meta.loc[group.index, "exclusion_reason"] = "incomplete_LA_RA_pair"

    included = meta[meta["analysis_inclusion"]]
    pairs = included["donor_id"].nunique()
    if pairs != EXPECTED_PAIRS[accession]:
        raise AssertionError(
            f"{accession}: expected {EXPECTED_PAIRS[accession]} pairs, found {pairs}"
        )
    for _, group in included.groupby("donor_id"):
        if len(group) != 2 or set(group["chamber"]) != {"LA", "RA"}:
            raise AssertionError(f"{accession}: invalid included pair")
    return meta


def read_gpl570_mapping(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = []
        table_started = False
        for line in handle:
            if line.startswith("ID\t"):
                table_started = True
            if table_started:
                lines.append(line.rstrip("\n"))
    reader = csv.DictReader(lines, delimiter="\t")
    rows = []
    for row in reader:
        probe = clean_field(row.get("ID", ""))
        raw_symbol = clean_field(row.get("Gene symbol", ""))
        symbols = [
            symbol.strip()
            for symbol in re.split(r"\s*///\s*", raw_symbol)
            if symbol.strip() and symbol.strip() != "---"
        ]
        if probe and symbols:
            rows.append(
                {
                    "source_id": probe,
                    "gene_symbol": symbols[0],
                    "mapping_source": "NCBI GEO GPL570.annot.gz; first listed symbol",
                }
            )
    mapping = pd.DataFrame(rows).drop_duplicates("source_id", keep="first")
    if len(mapping) < 20_000:
        raise AssertionError(f"Unexpectedly small GPL570 mapping: {len(mapping)}")
    return mapping


def collapse_array(
    accession: str,
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    mapping: pd.DataFrame,
) -> pd.DataFrame:
    included = metadata.loc[metadata["analysis_inclusion"], "sample_id"].tolist()
    missing_samples = sorted(set(included) - set(matrix.columns))
    if missing_samples:
        raise AssertionError(f"{accession}: matrix missing samples {missing_samples}")
    values = matrix[included].replace([np.inf, -np.inf], np.nan)
    if values.isna().any().any():
        raise AssertionError(f"{accession}: non-finite processed matrix values")

    common = values.index.intersection(mapping["source_id"])
    values = values.loc[common]
    map_index = mapping.set_index("source_id")
    symbols = map_index.loc[values.index, "gene_symbol"]
    iqr = values.quantile(0.75, axis=1) - values.quantile(0.25, axis=1)
    selector = pd.DataFrame(
        {"probe_id": values.index, "gene_symbol": symbols.to_numpy(), "iqr": iqr.to_numpy()}
    )
    selector = selector.sort_values(
        ["gene_symbol", "iqr", "probe_id"], ascending=[True, False, True]
    ).drop_duplicates("gene_symbol", keep="first")
    collapsed = values.loc[selector["probe_id"]].copy()
    collapsed.index = selector["gene_symbol"].to_numpy()
    collapsed.index.name = "gene_symbol"
    collapsed = collapsed.sort_index()
    if collapsed.index.duplicated().any():
        raise AssertionError(f"{accession}: duplicate symbols after collapse")
    return collapsed


def main() -> None:
    metadata_parts = []
    series_matrices: dict[str, pd.DataFrame | None] = {}
    for accession, path in SERIES.items():
        raw_meta, matrix = read_series(path)
        if len(raw_meta) != EXPECTED_SAMPLES[accession]:
            raise AssertionError(
                f"{accession}: expected {EXPECTED_SAMPLES[accession]} samples, "
                f"found {len(raw_meta)}"
            )
        metadata = build_metadata(accession, raw_meta)
        metadata.to_csv(
            META_DIR / f"{accession}_sample_metadata.tsv", sep="\t", index=False
        )
        metadata_parts.append(metadata)
        series_matrices[accession] = matrix

    unified = pd.concat(metadata_parts, ignore_index=True)
    required = [
        "dataset",
        "sample_id",
        "donor_id",
        "chamber",
        "rhythm",
        "platform",
        "analysis_inclusion",
        "exclusion_reason",
    ]
    unified[required + ["title", "source_expression_column"]].to_csv(
        META_DIR / "unified_sample_metadata.tsv", sep="\t", index=False
    )

    mapping = read_gpl570_mapping(META_DIR / "GPL570.annot.gz")
    mapping.to_csv(META_DIR / "gene_identifier_mapping.tsv", sep="\t", index=False)
    for accession in ("GSE79768", "GSE115574"):
        matrix = series_matrices[accession]
        if matrix is None:
            raise AssertionError(f"{accession}: Series Matrix has no expression table")
        metadata = unified[unified["dataset"] == accession]
        collapsed = collapse_array(accession, matrix, metadata, mapping)
        out_dir = PROCESSED / accession
        out_dir.mkdir(parents=True, exist_ok=True)
        collapsed.to_csv(out_dir / "expr_symbol.tsv", sep="\t")

    summary = (
        unified.groupby("dataset", as_index=False)
        .agg(
            n_registry_samples=("sample_id", "nunique"),
            n_analysis_samples=("analysis_inclusion", "sum"),
            n_registry_donors=("donor_id", "nunique"),
        )
    )
    included = unified[unified["analysis_inclusion"]]
    pair_counts = included.groupby("dataset")["donor_id"].nunique()
    summary["n_complete_pairs"] = summary["dataset"].map(pair_counts).astype(int)
    RESULTS.mkdir(parents=True, exist_ok=True)
    summary.to_csv(RESULTS / "cohort_pairing_summary.tsv", sep="\t", index=False)
    print(f"Wrote frozen metadata for {len(unified)} registry samples")
    print("Wrote GPL570 symbol matrices for GSE79768 and GSE115574")


if __name__ == "__main__":
    main()
