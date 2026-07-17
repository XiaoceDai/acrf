#!/usr/bin/env python3
"""P1 — Sample metadata 'identity cards' for the three AF chamber-rewiring cohorts."""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path
from typing import Dict, List

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "Project_P1" / "_meta"
OUTDIR = ROOT / "Project_P1"

COLUMNS = [
    "gsm", "cohort", "cohort_role", "platform",
    "donor_id", "chamber", "rhythm",
    "sex", "age", "disease_background",
    "cel_file", "count_column",
    "is_paired", "donor_pair_mate_gsm",
]


# series_matrix parser (minimal, only what we need)

def parse_series_matrix(path: Path) -> pd.DataFrame:
    """Return a DataFrame with one row per sample, indexed by GSM.

    Keeps the headers we care about: title, geo_accession, source_name,
    and every characteristics_ch1 row (joined into `characteristics`).
    """
    rows: Dict[str, List[str]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0][len("!Sample_"):]
            values = [v.strip().strip('"') for v in parts[1:]]
            rows.setdefault(key, []).append(values)

    if "geo_accession" not in rows:
        raise RuntimeError(f"no geo_accession row in {path}")
    gsms = rows["geo_accession"][0]
    n = len(gsms)

    df = pd.DataFrame({"gsm": gsms})
    df["title"] = rows["title"][0] if "title" in rows else ["" for _ in range(n)]
    df["source_name"] = (
        rows["source_name_ch1"][0] if "source_name_ch1" in rows else ["" for _ in range(n)]
    )
    # Concatenate all characteristics_ch1 rows with " | " so we can regex into them.
    char_rows = rows.get("characteristics_ch1", [])
    joined = []
    for i in range(n):
        joined.append(" | ".join(cr[i] for cr in char_rows))
    df["characteristics"] = joined
    return df


def _extract(pattern: str, text: str, default: str | None = None) -> str | None:
    m = re.search(pattern, text, flags=re.IGNORECASE)
    return m.group(1).strip() if m else default


# cohort-specific parsers
# Each returns a DataFrame with the full locked schema.

def build_gse128188(series_matrix: Path) -> pd.DataFrame:
    """GSE128188 — titles are 'SR_LA_7', 'AF_RA_5' etc. Paired LAA/RAA bulk RNA-seq."""
    df = parse_series_matrix(series_matrix)
    if len(df) != 20:
        raise RuntimeError(f"GSE128188 expected 20 samples, got {len(df)}")

    out = pd.DataFrame()
    out["gsm"] = df["gsm"]
    out["cohort"] = "GSE128188"
    out["cohort_role"] = "discovery"
    out["platform"] = "bulk_rnaseq"

    pat = re.compile(r"^(SR|AF)_(LA|RA)_(\S+)$")
    parsed = df["title"].map(lambda t: pat.match(t))
    if parsed.isna().any() or (parsed == None).any():  # noqa: E711
        bad = df.loc[parsed.isna() | (parsed == None), "title"].tolist()  # noqa: E711
        raise RuntimeError(f"GSE128188 title regex failed on: {bad}")

    out["rhythm"] = [m.group(1) for m in parsed]
    out["chamber"] = [m.group(2) for m in parsed]
    out["donor_id"] = [m.group(3) for m in parsed]

    # Counts column header in GSE128188_CountsV2_EdgeR.txt.gz is "<donor>-<chamber>".
    out["count_column"] = out["donor_id"] + "-" + out["chamber"]
    out["cel_file"] = pd.NA

    # Optional demographics (128188 is all-male; age not in this series_matrix).
    out["sex"] = df["characteristics"].map(
        lambda s: _extract(r"Sex:\s*(\w+)", s)
    ).str.lower()
    out["age"] = pd.to_numeric(
        df["characteristics"].map(lambda s: _extract(r"age:\s*(\d+)", s)),
        errors="coerce",
    )
    out["disease_background"] = pd.NA
    return out


def build_gse79768(series_matrix: Path) -> pd.DataFrame:
    """GSE79768 — titles are 'Atrial Fibrillation_left atrial_14' etc. 13 paired donors."""
    df = parse_series_matrix(series_matrix)
    if len(df) != 26:
        raise RuntimeError(f"GSE79768 expected 26 samples, got {len(df)}")

    out = pd.DataFrame()
    out["gsm"] = df["gsm"]
    out["cohort"] = "GSE79768"
    out["cohort_role"] = "discovery"
    out["platform"] = "affy_hgu133plus2"

    # Title format: "<Atrial Fibrillation|Sinus Rhythm>_<left atrial|right atrial>_<donor>"
    pat = re.compile(
        r"^(Atrial Fibrillation|Sinus Rhythm)_"
        r"(left atrial|right atrial)_"
        r"(\S+)$"
    )
    parsed = df["title"].map(lambda t: pat.match(t))
    if parsed.isna().any() or (parsed == None).any():  # noqa: E711
        bad = df.loc[parsed.isna() | (parsed == None), "title"].tolist()  # noqa: E711
        raise RuntimeError(f"GSE79768 title regex failed on: {bad}")

    rhythm_map = {"Atrial Fibrillation": "AF", "Sinus Rhythm": "SR"}
    chamber_map = {"left atrial": "LA", "right atrial": "RA"}
    out["rhythm"] = [rhythm_map[m.group(1)] for m in parsed]
    out["chamber"] = [chamber_map[m.group(2)] for m in parsed]
    out["donor_id"] = [m.group(3) for m in parsed]

    # CEL filename inside the tar is "<GSM>_<donor><L|R>.CEL.gz".
    cel_suffix = {"LA": "L", "RA": "R"}
    out["cel_file"] = [
        f"{g}_{d}{cel_suffix[c]}.CEL.gz"
        for g, d, c in zip(out["gsm"], out["donor_id"], out["chamber"])
    ]
    out["count_column"] = pd.NA

    out["sex"] = df["characteristics"].map(
        lambda s: _extract(r"gender:\s*(\w+)", s)
    ).str.lower()
    out["age"] = pd.to_numeric(
        df["characteristics"].map(lambda s: _extract(r"age:\s*(\d+)", s)),
        errors="coerce",
    )
    out["disease_background"] = pd.NA
    return out


def build_gse115574(series_matrix: Path) -> pd.DataFrame:
    """GSE115574 — titles are 'AF01L', 'SR16R' etc.; GEO registry 15 AF + 16 SR donors.

    Analytical Fig 2 donor totals (13 AF + 15 SR) are resolved downstream against
    innovation_ncp_v3_acai_long.tsv / p1_demographics_table.py, not in this parser.
    """
    df = parse_series_matrix(series_matrix)
    if len(df) != 59:
        raise RuntimeError(f"GSE115574 expected 59 samples, got {len(df)}")

    out = pd.DataFrame()
    out["gsm"] = df["gsm"]
    out["cohort"] = "GSE115574"
    out["cohort_role"] = "replication"
    out["platform"] = "affy_hgu133plus2"

    # Title format: "<AF|SR><nn><L|R>"
    pat = re.compile(r"^(AF|SR)(\d{2})([LR])$")
    parsed = df["title"].map(lambda t: pat.match(t))
    if parsed.isna().any() or (parsed == None).any():  # noqa: E711
        bad = df.loc[parsed.isna() | (parsed == None), "title"].tolist()  # noqa: E711
        raise RuntimeError(f"GSE115574 title regex failed on: {bad}")

    chamber_map = {"L": "LA", "R": "RA"}
    out["rhythm"] = [m.group(1) for m in parsed]
    # Donor id keeps the rhythm prefix so "AF01" and "SR01" don't collide.
    out["donor_id"] = [f"{m.group(1)}{m.group(2)}" for m in parsed]
    out["chamber"] = [chamber_map[m.group(3)] for m in parsed]

    out["cel_file"] = [
        f"{g}_{t}.CEL.gz" for g, t in zip(out["gsm"], df["title"])
    ]
    out["count_column"] = pd.NA

    # 115574 series_matrix does not systematically expose sex/age; leave NA
    # unless we later decide we need it for Supp.
    out["sex"] = pd.NA
    out["age"] = pd.NA
    # Every sample in this series is either AF+MR or SR+MR; the disease background
    # is the MR, which is the whole reason this cohort is separated from discovery.
    out["disease_background"] = df["source_name"].map(
        lambda s: "severe_mitral_regurgitation" if "patient" in s.lower() else pd.NA
    )
    return out


# pairing annotation

def annotate_pairing(df: pd.DataFrame) -> pd.DataFrame:
    """Fill is_paired and donor_pair_mate_gsm by (cohort, donor_id) aggregation."""
    df = df.copy()
    df["is_paired"] = False
    df["donor_pair_mate_gsm"] = pd.NA

    for (cohort, donor), grp in df.groupby(["cohort", "donor_id"]):
        chambers = set(grp["chamber"])
        if chambers == {"LA", "RA"} and len(grp) == 2:
            # complete pair — annotate each row with the other's GSM
            la_gsm = grp.loc[grp["chamber"] == "LA", "gsm"].iloc[0]
            ra_gsm = grp.loc[grp["chamber"] == "RA", "gsm"].iloc[0]
            df.loc[(df["cohort"] == cohort) & (df["donor_id"] == donor) & (df["chamber"] == "LA"),
                   ["is_paired", "donor_pair_mate_gsm"]] = [True, ra_gsm]
            df.loc[(df["cohort"] == cohort) & (df["donor_id"] == donor) & (df["chamber"] == "RA"),
                   ["is_paired", "donor_pair_mate_gsm"]] = [True, la_gsm]
    return df


# reporting

def summarise(df: pd.DataFrame) -> str:
    lines = [f"--- {df['cohort'].iloc[0]} ({df['cohort_role'].iloc[0]}, {df['platform'].iloc[0]}) ---"]
    lines.append(f"total samples: {len(df)}")

    donors = df.groupby("donor_id")
    paired = sum(1 for _, g in donors if set(g["chamber"]) == {"LA", "RA"})
    lines.append(f"unique donors: {donors.ngroups}  (paired LA+RA: {paired})")

    by_rhythm = df.groupby("rhythm")["donor_id"].nunique()
    lines.append(f"donors by rhythm: {by_rhythm.to_dict()}")

    paired_by_rhythm: Dict[str, int] = {"SR": 0, "AF": 0}
    for _, g in donors:
        if set(g["chamber"]) == {"LA", "RA"}:
            paired_by_rhythm[g["rhythm"].iloc[0]] += 1
    lines.append(f"paired donors by rhythm: {paired_by_rhythm}")

    unpaired = [
        (r["donor_id"], r["chamber"], r["rhythm"], r["gsm"])
        for _, r in df[~df["is_paired"]].iterrows()
    ]
    if unpaired:
        lines.append(f"unpaired samples ({len(unpaired)}): " +
                     ", ".join(f"{d}/{c}/{rh}/{g}" for d, c, rh, g in unpaired))
    return "\n".join(lines)


# main

def main() -> None:
    if not META.exists():
        sys.exit(f"metadata dir not found: {META}")

    builders = {
        "GSE128188": (build_gse128188, META / "GSE128188_series_matrix.txt.gz"),
        "GSE79768":  (build_gse79768,  META / "GSE79768_series_matrix.txt.gz"),
        "GSE115574": (build_gse115574, META / "GSE115574_series_matrix.txt.gz"),
    }

    per_cohort: Dict[str, pd.DataFrame] = {}
    for cohort, (fn, path) in builders.items():
        if not path.exists():
            sys.exit(f"missing series_matrix: {path}")
        df = fn(path)
        df = annotate_pairing(df)
        df = df[COLUMNS]
        out = OUTDIR / cohort / "sample_metadata.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, sep="\t", index=False)
        per_cohort[cohort] = df
        print(summarise(df))
        print(f"  -> {out.relative_to(ROOT)}\n")

    union = pd.concat(per_cohort.values(), ignore_index=True)
    union_out = OUTDIR / "samples_union.tsv"
    union.to_csv(union_out, sep="\t", index=False)
    print(f"union: {len(union)} rows -> {union_out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
