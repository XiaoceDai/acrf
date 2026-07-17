#!/usr/bin/env python3
"""Download the approved minimum GEO inputs into an immutable namespace."""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import shutil
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "results" / "ncp_v3_singscore_reanalysis" / "cohort_input_audit.tsv"
DATA = ROOT / "data" / "ncp_v3_singscore_reanalysis"
RAW = DATA / "raw"
META = DATA / "metadata"
OUT = ROOT / "results" / "ncp_v3_singscore_reanalysis" / "public_input_manifest.tsv"
LEGACY_ROOT = ROOT / "Project_P2" / "Project_P1"
USER_AGENT = "NCP-v3-singscore-input-ingestion/1.0"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        with temporary.open("wb") as handle:
            shutil.copyfileobj(response, handle, length=1024 * 1024)
    temporary.replace(destination)


def destination_for(filename: str) -> Path:
    if "series_matrix" in filename or filename.startswith("GPL"):
        return META / filename
    return RAW / filename


def legacy_candidate(accession: str, filename: str) -> Path | None:
    candidates = [
        LEGACY_ROOT / "_meta" / filename,
        LEGACY_ROOT / accession / filename,
        ROOT / "bulk" / "Project_P01" / accession / filename,
    ]
    return next((path for path in candidates if path.is_file()), None)


def main() -> None:
    if not AUDIT.is_file():
        raise FileNotFoundError(f"Missing approved audit table: {AUDIT}")
    audit_rows = list(csv.DictReader(AUDIT.open(encoding="utf-8"), delimiter="\t"))
    if len(audit_rows) != 4:
        raise ValueError(f"Expected four audit rows, found {len(audit_rows)}")
    if not all(row["remote_audit_status"].startswith("recommended files reachable") for row in audit_rows):
        raise RuntimeError("Remote audit did not pass for every cohort")

    requests: dict[str, dict[str, str]] = {}
    for row in audit_rows:
        names = [item.strip() for item in row["recommended_authoritative_file"].split(";")]
        urls = [item.strip() for item in row["url"].split(";")]
        if len(names) != len(urls):
            raise ValueError(f"Filename/URL count mismatch for {row['accession']}")
        for filename, url in zip(names, urls):
            existing = requests.get(url)
            record = {
                "accession": row["accession"],
                "filename": filename,
                "url": url,
                "platform": row["platform"],
                "data_level": row["data_level"],
            }
            if existing is None:
                requests[url] = record
            elif filename != existing["filename"]:
                raise ValueError(f"One URL maps to multiple filenames: {url}")
            else:
                existing["accession"] = f"{existing['accession']};{row['accession']}"

    download_date = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    manifest: list[dict[str, object]] = []
    for record in requests.values():
        destination = destination_for(record["filename"])
        download(record["url"], destination)
        digest = sha256(destination)
        accessions = record["accession"].split(";")
        comparisons = []
        for accession in accessions:
            legacy = legacy_candidate(accession, record["filename"])
            if legacy is None:
                comparisons.append(f"{accession}:not_found")
            else:
                legacy_digest = sha256(legacy)
                comparisons.append(
                    f"{accession}:{'identical' if legacy_digest == digest else 'different'}:"
                    f"{legacy.relative_to(ROOT)}:{legacy_digest}"
                )
        manifest.append(
            {
                "accession": record["accession"],
                "source_url": record["url"],
                "download_date_utc": download_date,
                "original_filename": record["filename"],
                "local_path": str(destination.relative_to(ROOT)),
                "file_size_bytes": destination.stat().st_size,
                "sha256": digest,
                "platform": record["platform"],
                "annotation_build": (
                    "GPL570 annotation snapshot"
                    if record["filename"].startswith("GPL570")
                    else "GEO record / filename"
                ),
                "data_level": record["data_level"],
                "preprocessing_status": "downloaded_unmodified",
                "legacy_comparison": ";".join(comparisons),
            }
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(manifest[0])
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Wrote {OUT.relative_to(ROOT)} ({len(manifest)} unique files)")
    print(f"Downloaded bytes: {sum(int(row['file_size_bytes']) for row in manifest)}")


if __name__ == "__main__":
    main()
