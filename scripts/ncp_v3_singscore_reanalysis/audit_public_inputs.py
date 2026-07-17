#!/usr/bin/env python3
"""Audit minimum public inputs for the NCP_v3 singscore reanalysis.

This stage performs metadata-only HTTP HEAD requests. It intentionally does not
download expression matrices, CEL archives, or FASTQ files. The generated table
is the approval gate for the later immutable-input ingestion stage.
"""
from __future__ import annotations

import csv
import datetime as dt
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "ncp_v3_singscore_reanalysis" / "cohort_input_audit.tsv"
USER_AGENT = "NCP-v3-singscore-input-audit/1.0"

GPL570_ANNOT = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/"
    "GPL570.annot.gz"
)


@dataclass(frozen=True)
class RemoteFile:
    name: str
    url: str
    role: str
    recommended: bool
    large_download: bool = False


@dataclass(frozen=True)
class Cohort:
    accession: str
    platform: str
    data_level: str
    gene_identifier_type: str
    sample_count: int
    analysis_sample_count: int
    expected_paired_donor_count: int
    files: tuple[RemoteFile, ...]
    raw_reprocessing_required: str
    recommendation_reason: str


COHORTS = (
    Cohort(
        accession="GSE112339",
        platform="GPL16791; Illumina HiSeq 2500; RNA-seq",
        data_level="deposited gene-level HTSeq counts",
        gene_identifier_type=(
            "HGNC symbols in the deposited hg19 matrix despite the Ensembl filename; "
            "must be re-verified after download"
        ),
        sample_count=40,
        analysis_sample_count=16,
        expected_paired_donor_count=8,
        files=(
            RemoteFile(
                "GSE112339_HTSeq_hg19_Ensembl.txt.gz",
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/"
                "GSE112339/suppl/GSE112339_HTSeq_hg19_Ensembl.txt.gz",
                "recommended expression",
                True,
            ),
            RemoteFile(
                "GSE112339_series_matrix.txt.gz",
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/"
                "GSE112339/matrix/GSE112339_series_matrix.txt.gz",
                "recommended metadata",
                True,
            ),
            RemoteFile(
                "GSE112339_HTSeq_hg38_NONCODE2016.txt.gz",
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/"
                "GSE112339/suppl/GSE112339_HTSeq_hg38_NONCODE2016.txt.gz",
                "alternate processed matrix; not recommended",
                False,
            ),
        ),
        raw_reprocessing_required="no",
        recommendation_reason=(
            "The deposited hg19 gene-level count matrix is sufficient to reconstruct "
            "within-sample ranks; the small Series Matrix supplies authoritative GSM, "
            "donor, and tissue metadata. FASTQ reprocessing is unnecessary."
        ),
    ),
    Cohort(
        accession="GSE128188",
        platform="GPL18573; Illumina NextSeq 500; RNA-seq",
        data_level="deposited gene-level HTSeq raw counts",
        gene_identifier_type="HGNC gene symbols",
        sample_count=20,
        analysis_sample_count=20,
        expected_paired_donor_count=10,
        files=(
            RemoteFile(
                "GSE128188_CountsV2_EdgeR.txt.gz",
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE128nnn/"
                "GSE128188/suppl/GSE128188_CountsV2_EdgeR.txt.gz",
                "recommended expression",
                True,
            ),
            RemoteFile(
                "GSE128188_series_matrix.txt.gz",
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE128nnn/"
                "GSE128188/matrix/GSE128188_series_matrix.txt.gz",
                "recommended metadata",
                True,
            ),
        ),
        raw_reprocessing_required="no",
        recommendation_reason=(
            "The author-deposited complete gene-count matrix covers all 20 paired "
            "samples and preserves sample-level ranks; SRA/FASTQ reprocessing adds no "
            "required information for this estimand."
        ),
    ),
    Cohort(
        accession="GSE79768",
        platform="GPL570; Affymetrix Human Genome U133 Plus 2.0 Array",
        data_level="deposited normalized probe-level Series Matrix",
        gene_identifier_type="Affymetrix GPL570 probe set ID; map to HGNC symbol",
        sample_count=26,
        analysis_sample_count=26,
        expected_paired_donor_count=13,
        files=(
            RemoteFile(
                "GSE79768_series_matrix.txt.gz",
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/"
                "GSE79768/matrix/GSE79768_series_matrix.txt.gz",
                "recommended expression and metadata",
                True,
            ),
            RemoteFile(
                "GPL570.annot.gz",
                GPL570_ANNOT,
                "recommended probe annotation",
                True,
            ),
            RemoteFile(
                "GSE79768_RAW.tar",
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/"
                "GSE79768/suppl/GSE79768_RAW.tar",
                "CEL fallback only",
                False,
                True,
            ),
        ),
        raw_reprocessing_required="no; conditional CEL fallback only",
        recommendation_reason=(
            "GEO reports processed values in the complete Series Matrix for all 26 "
            "samples. It is adequate for within-sample ranking if post-download checks "
            "confirm finite values, complete sample metadata, and resolvable GPL570 "
            "probe IDs; otherwise the CEL archive becomes required."
        ),
    ),
    Cohort(
        accession="GSE115574",
        platform="GPL570; Affymetrix Human Genome U133 Plus 2.0 Array",
        data_level="deposited normalized probe-level Series Matrix",
        gene_identifier_type="Affymetrix GPL570 probe set ID; map to HGNC symbol",
        sample_count=59,
        analysis_sample_count=56,
        expected_paired_donor_count=28,
        files=(
            RemoteFile(
                "GSE115574_series_matrix.txt.gz",
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/"
                "GSE115574/matrix/GSE115574_series_matrix.txt.gz",
                "recommended expression and metadata",
                True,
            ),
            RemoteFile(
                "GPL570.annot.gz",
                GPL570_ANNOT,
                "recommended probe annotation",
                True,
            ),
            RemoteFile(
                "GSE115574_RAW.tar",
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/"
                "GSE115574/suppl/GSE115574_RAW.tar",
                "CEL fallback only",
                False,
                True,
            ),
        ),
        raw_reprocessing_required="no; conditional CEL fallback only",
        recommendation_reason=(
            "GEO reports processed values in the 59-sample Series Matrix. The 28 "
            "complete LA/RA donor pairs can be selected from metadata; AF14, AF15, and "
            "SR16 remain excluded as unpaired. CEL reprocessing is only needed if the "
            "processed matrix fails completeness, reproducibility, or probe-map checks."
        ),
    ),
)


def head_size(url: str) -> tuple[int | None, str]:
    request = urllib.request.Request(
        url,
        method="HEAD",
        headers={"User-Agent": USER_AGENT, "Accept-Encoding": "identity"},
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            raw_size = response.headers.get("Content-Length")
            size = int(raw_size) if raw_size is not None else None
            return size, f"HTTP {response.status}"
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        return None, f"HEAD failed: {type(exc).__name__}: {exc}"


def format_size(size: int | None) -> str:
    if size is None:
        return "unavailable"
    units = ("B", "KiB", "MiB", "GiB")
    value = float(size)
    unit = units[0]
    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            break
        value /= 1024.0
    return f"{size} B ({value:.2f} {unit})"


def main() -> None:
    checked: dict[str, tuple[int | None, str]] = {}
    for cohort in COHORTS:
        for remote in cohort.files:
            if remote.url not in checked:
                checked[remote.url] = head_size(remote.url)

    audit_date = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    rows: list[dict[str, object]] = []
    for cohort in COHORTS:
        recommended = [item for item in cohort.files if item.recommended]
        available_parts = []
        all_recommended_reachable = True
        for item in cohort.files:
            size, status = checked[item.url]
            available_parts.append(
                f"{item.name} [{item.role}; {format_size(size)}; {status}]"
            )
            if item.recommended and size is None:
                all_recommended_reachable = False

        rows.append(
            {
                "accession": cohort.accession,
                "platform": cohort.platform,
                "available_files": "; ".join(available_parts),
                "recommended_authoritative_file": "; ".join(
                    item.name for item in recommended
                ),
                "data_level": cohort.data_level,
                "gene_identifier_type": cohort.gene_identifier_type,
                "sample_count": cohort.sample_count,
                "analysis_sample_count": cohort.analysis_sample_count,
                "expected_paired_donor_count": cohort.expected_paired_donor_count,
                "url": "; ".join(item.url for item in recommended),
                "file_size": "; ".join(
                    f"{item.name}={format_size(checked[item.url][0])}"
                    for item in recommended
                ),
                "raw_reprocessing_required": cohort.raw_reprocessing_required,
                "reason_for_recommendation": cohort.recommendation_reason,
                "remote_audit_status": (
                    "recommended files reachable by HEAD"
                    if all_recommended_reachable
                    else "BLOCKED: one or more recommended file sizes unavailable"
                ),
                "large_download_deferred": True,
                "audit_date_utc": audit_date,
            }
        )

    required_columns = [
        "accession",
        "platform",
        "available_files",
        "recommended_authoritative_file",
        "data_level",
        "gene_identifier_type",
        "sample_count",
        "analysis_sample_count",
        "expected_paired_donor_count",
        "url",
        "file_size",
        "raw_reprocessing_required",
        "reason_for_recommendation",
        "remote_audit_status",
        "large_download_deferred",
        "audit_date_utc",
    ]
    assert len(rows) == 4
    assert {row["accession"] for row in rows} == {
        "GSE112339",
        "GSE128188",
        "GSE79768",
        "GSE115574",
    }
    assert {row["expected_paired_donor_count"] for row in rows} == {8, 10, 13, 28}
    assert all(row["large_download_deferred"] is True for row in rows)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=required_columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    blocked = [row["accession"] for row in rows if str(row["remote_audit_status"]).startswith("BLOCKED")]
    print(f"Wrote {OUT.relative_to(ROOT)} ({len(rows)} cohorts)")
    print("No expression matrix, CEL archive, or FASTQ file was downloaded.")
    if blocked:
        raise SystemExit(f"Remote audit incomplete for: {', '.join(blocked)}")
    print("APPROVAL REQUIRED: stop before downloading the recommended files.")


if __name__ == "__main__":
    main()
