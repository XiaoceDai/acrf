#!/usr/bin/env python3
"""
Build the cohort-demographics master table (Supplementary Table S1 / Table 1).

Includes GSE112339 healthy bulk reference plus discovery, replication, scRNA,
multiome, and GeoMx cohorts — eight GEO-linked resources when all rows present.

Output: results/cohort_demographics_table.tsv  (cohort-level; versioned source)
        results/cohort_demographics_long.tsv   (per-donor; Fig S1 input)

For GSE128188 / GSE79768 / GSE115574, donor lists and counts are intersected with
results/innovation_ncp_v3_acai_long.tsv when present so Table S1 / Fig S1 match Fig 2D
sample sizes (n_AF / n_SR). GEO metadata may list additional donors omitted from the
v3 NCP axis panel.

DOCX: run scripts/m_supp_table_s1_cohort_demographics.py after this script.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "results" / "cohort_demographics_table.tsv"
OUT_LONG = ROOT / "results" / "cohort_demographics_long.tsv"

# Same ordering as scripts/figS1_cohort_overview.cohort_order (manuscript flow).
COHORT_SORT_ORDER = [
    "GSE112339", "GSE128188", "GSE79768", "GSE115574",
    "GSE255612", "GSE41177", "GSE238242", "GSE261363",
]

FIG2_COHORT_IDS = frozenset({"GSE128188", "GSE79768", "GSE115574"})
V3_ACAI_LONG = ROOT / "results" / "innovation_ncp_v3_acai_long.tsv"


def ncp_v3_panel_donor_ids(cohort_id: str) -> frozenset[str] | None:
    """Donors present in v3 axis aCAI long table (Fig 2 input); None if path missing."""
    if cohort_id not in FIG2_COHORT_IDS:
        return None
    if not V3_ACAI_LONG.exists():
        return None
    chunk = pd.read_csv(V3_ACAI_LONG, sep="\t", usecols=["cohort", "donor_id"])
    s = chunk.loc[chunk["cohort"] == cohort_id, "donor_id"].dropna().astype(str)
    return frozenset(s.unique())


COHORTS = [
    {"id": "GSE128188", "role": "Fig 2 paired discovery (RNA-seq)",
     "path": "Project_P1/GSE128188/sample_metadata.tsv",
     "platform": "Illumina HiSeq, paired chamber",
     "tissue": "LAA + RAA"},
    {"id": "GSE79768",  "role": "Fig 2 paired discovery (microarray)",
     "path": "Project_P1/GSE79768/sample_metadata.tsv",
     "platform": "Affy HG-U133 Plus 2.0, paired chamber",
     "tissue": "LAA + RAA"},
    {"id": "GSE115574", "role": "Fig 2 replication (MR background)",
     "path": "Project_P1/GSE115574/sample_metadata.tsv",
     "platform": "Affy HG-U133 Plus 2.0, paired chamber",
     "tissue": "LAA + RAA"},
    {"id": "GSE41177",  "role": "Fig 4 LAA + LAPV duration cohort",
     "path": "Project_P4/GSE41177/sample_metadata.tsv",
     "platform": "Affy HG-U133 Plus 2.0, valvular AF",
     "tissue": "LAA + LAPV"},
]


def fmt_age(values):
    arr = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if len(arr) == 0:
        return "n/a"
    return f"{arr.mean():.0f} ± {arr.std(ddof=0):.0f}  (range {int(arr.min())}–{int(arr.max())})"


def fmt_sex(values):
    s = pd.Series(values).astype(str).str.lower().str.strip()
    s = s.replace({"nan": "", "none": ""}).fillna("")
    n_m = s.isin({"male", "m"}).sum()
    n_f = s.isin({"female", "f"}).sum()
    n_u = len(s) - n_m - n_f
    parts = [f"M={n_m}"]
    if n_f: parts.append(f"F={n_f}")
    if n_u: parts.append(f"NA={n_u}")
    return "  ".join(parts)


rows = []
long_rows = []


def append_gse112339_healthy_reference() -> None:
    """Healthy paired-chamber bulk reference (Normal Chamber Prior)."""
    bulk_meta = ROOT / "bulk" / "sample_meta.tsv"
    alt_meta = ROOT / "Project_P1" / "GSE112339" / "sample_metadata.tsv"
    n_donors = 8
    donor_ids: list[str] = []

    if bulk_meta.exists():
        sm = pd.read_csv(bulk_meta, sep="\t")
        if "donor" in sm.columns:
            donor_ids = sm["donor"].drop_duplicates().astype(str).tolist()
            n_donors = len(donor_ids)
    elif alt_meta.exists():
        md = pd.read_csv(alt_meta, sep="\t")
        if "donor_id" in md.columns:
            donor_ids = md["donor_id"].drop_duplicates().astype(str).tolist()
            n_donors = len(donor_ids)

    rows.append({
        "cohort": "GSE112339",
        "role": "Fig 1 Normal Chamber Prior (healthy bulk reference)",
        "platform": "Illumina HiSeq RNA-seq, paired-end (GEO GSE112339)",
        "chamber": "LAA + RAA",
        "n_donors": n_donors,
        "n_AF_donors": 0,
        "n_SR_donors": 0,
        "age": "n/a",
        "sex": "n/a",
        "background": (
            "healthy paired LA/RA donors; chamber reference for aCAI_bulk and "
            "NCP panel construction"
        ),
    })
    for did in donor_ids:
        long_rows.append({
            "cohort": "GSE112339",
            "donor_id": did,
            "rhythm": "healthy_chamber",
            "sex": None,
            "age": np.nan,
            "disease_background": "healthy reference",
        })


append_gse112339_healthy_reference()

for c in COHORTS:
    p = ROOT / c["path"]
    if not p.exists():
        rows.append({"cohort": c["id"], "role": c["role"],
                     "platform": c["platform"], "chamber": c["tissue"],
                     "n_donors": "n/a", "n_AF_donors": "n/a",
                     "n_SR_donors": "n/a", "age": "n/a", "sex": "n/a",
                     "background": "n/a"})
        continue
    md = pd.read_csv(p, sep="\t")
    if "donor_id" not in md.columns or "rhythm" not in md.columns:
        donor_meta = pd.DataFrame()
    else:
        agg_dict = {}
        for col, default in [("rhythm", "n/a"), ("sex", None),
                              ("age", np.nan),
                              ("disease_background", None)]:
            if col in md.columns:
                agg_dict[col] = (lambda x: x.dropna().iloc[0]
                                            if x.dropna().any() else default)
        donor_meta = (md.groupby("donor_id")
                        .agg(agg_dict)
                        .reset_index())

    panel_ids = ncp_v3_panel_donor_ids(c["id"])
    n_geo = len(donor_meta)
    if panel_ids is not None and len(donor_meta) > 0:
        filt = donor_meta[
            donor_meta["donor_id"].astype(str).isin(panel_ids)
        ].copy().reset_index(drop=True)
        if filt.empty:
            raise ValueError(
                f"{c['id']}: no GEO donors overlap "
                f"{V3_ACAI_LONG.relative_to(ROOT)} — check cohort IDs / metadata."
            )
        donor_meta = filt

    n_donors = len(donor_meta)
    n_af = (donor_meta["rhythm"] == "AF").sum() if "rhythm" in donor_meta else 0
    n_sr = (donor_meta["rhythm"] == "SR").sum() if "rhythm" in donor_meta else 0
    age_str = fmt_age(donor_meta["age"]) if "age" in donor_meta else "n/a"
    sex_str = fmt_sex(donor_meta["sex"]) if "sex" in donor_meta else "n/a"
    bg_vals = donor_meta["disease_background"].dropna().unique() \
        if "disease_background" in donor_meta else []
    bg_str = "; ".join(bg_vals) if len(bg_vals) > 0 else "primary AF"

    # add duration info for GSE41177
    duration_str = ""
    if c["id"] == "GSE41177" and "duration_months" in md.columns:
        af_rows = md[(md["rhythm"] == "AF") & (md["tissue"] == "LAA")]
        af_rows = af_rows.dropna(subset=["duration_months"])
        af_rows = af_rows.drop_duplicates("donor_id")
        if len(af_rows) > 0:
            durs = pd.to_numeric(af_rows["duration_months"], errors="coerce").dropna()
            duration_str = f"; AF duration n={len(durs)}, median={durs.median():.0f} mo (range {int(durs.min())}–{int(durs.max())})"

    panel_note = ""
    if panel_ids is not None and n_geo > len(donor_meta):
        panel_note = (
            f" GEO sample registry lists {n_geo} donor(s); Fig 2 / Fig S1 / Table S1 report "
            f"{len(donor_meta)} paired donors ({int(n_af)} AF + {int(n_sr)} SR) retained in "
            f"{V3_ACAI_LONG.name}. Exclusions reflect paired LAA/RAA completeness after "
            "preprocessing, overlap/probe coverage against frozen NCP programmes, and "
            "donor-level QC (Methods)."
        )

    rows.append({
        "cohort": c["id"], "role": c["role"], "platform": c["platform"],
        "chamber": c["tissue"], "n_donors": n_donors,
        "n_AF_donors": n_af, "n_SR_donors": n_sr,
        "age": age_str, "sex": sex_str,
        "background": bg_str + duration_str + panel_note,
    })
    for _, r in donor_meta.iterrows():
        long_rows.append({"cohort": c["id"], "donor_id": r["donor_id"],
                          "rhythm": r.get("rhythm", "n/a"),
                          "sex": r.get("sex"),
                          "age": r.get("age"),
                          "disease_background": r.get("disease_background")})

# Add GSE255612 (single-nucleus) row from snRNA metadata if available
sn_path = ROOT / "Project_P3" / "GSE255612" / "derived" / "gse255612_obs_4compartment.tsv"
if sn_path.exists():
    sn = pd.read_csv(sn_path, sep="\t")
    sn_donors = sn[["donor_id", "rhythm", "sex"]].drop_duplicates() \
        if "donor_id" in sn.columns else pd.DataFrame()
    n_d = len(sn_donors)
    n_af = (sn_donors["rhythm"] == "AF").sum() if "rhythm" in sn_donors else 0
    n_sr = (sn_donors["rhythm"] == "Ctrl").sum() if "rhythm" in sn_donors else 0
    rows.append({
        "cohort": "GSE255612", "role": "Fig 3 LA-arm snRNA pseudobulk",
        "platform": "10x Chromium snRNA-seq",
        "chamber": "LA + RA", "n_donors": n_d,
        "n_AF_donors": n_af, "n_SR_donors": n_sr,
        "age": "n/a (per-cell metadata)", "sex": fmt_sex(sn_donors["sex"]) if "sex" in sn_donors else "n/a",
        "background": "primary AF; donor-level QC retained",
    })

# Add GSE238242 (snRNA + snATAC multiome) from existing pseudobulk meta
sn_meta = ROOT / "Project_P5" / "GSE238242" / "pseudobulk" / "gse238242_pseudobulk_meta_main.tsv"
if sn_meta.exists():
    pm = pd.read_csv(sn_meta, sep="\t")
    donors = pm[["donor", "rhythm", "sex"]].drop_duplicates()
    n_d = len(donors)
    n_af = (donors["rhythm"] == "AF").sum()
    n_sr = (donors["rhythm"] == "SR").sum()
    rows.append({
        "cohort": "GSE238242", "role": "Fig S5 multiome RNA pseudobulk (gate)",
        "platform": "10x Multiome (snRNA + snATAC)",
        "chamber": "LAA",
        "n_donors": n_d, "n_AF_donors": n_af, "n_SR_donors": n_sr,
        "age": "n/a (per-cell metadata)",
        "sex": fmt_sex(donors["sex"]),
        "background": "persistent AF (NCR 2025); AF all male, SR mixed",
    })

# Add GSE261363 (GeoMx spatial) — counts ROIs
roi_meta = ROOT / "Project_P5" / "GSE261363" / "roi_metadata.tsv"
if roi_meta.exists():
    rm = pd.read_csv(roi_meta, sep="\t")
    rm = rm[rm["is_control"] == False] if "is_control" in rm.columns else rm
    n_pat = rm["patient"].nunique() if "patient" in rm.columns else 0
    n_af = (rm["rhythm"] == "AF").sum()
    n_sr = (rm["rhythm"] == "SR").sum()
    n_aa = (rm["tissue"] == "Atrial appendage").sum() if "tissue" in rm.columns else 0
    rows.append({
        "cohort": "GSE261363", "role": "Fig 5C ROI-level localization (GeoMx)",
        "platform": "NanoString GeoMx DSP",
        "chamber": "atrial appendage + EAT",
        "n_donors": n_pat, "n_AF_donors": "—", "n_SR_donors": "—",
        "age": "n/a",
        "sex": "n/a",
        "background": f"valvular AF + SR controls; n_ROIs (AF / SR) = {n_af} / {n_sr}; n_AA segments = {n_aa}",
    })

order_map = {c: i for i, c in enumerate(COHORT_SORT_ORDER)}
table = pd.DataFrame(rows)
table["_sort"] = table["cohort"].map(lambda x: order_map.get(x, 999))
table = table.sort_values("_sort").drop(columns=["_sort"])

OUT.parent.mkdir(parents=True, exist_ok=True)
table.to_csv(OUT, sep="\t", index=False, na_rep="n/a")

long_df = pd.DataFrame(long_rows)
if len(long_df) > 0:
    long_df["_sort"] = long_df["cohort"].map(lambda x: order_map.get(x, 999))
    long_df = long_df.sort_values("_sort").drop(columns=["_sort"])
long_df.to_csv(OUT_LONG, sep="\t", index=False, na_rep="n/a")

print(f"wrote {OUT.relative_to(ROOT)}  ({len(table)} cohort rows)")
print(f"wrote {OUT_LONG.relative_to(ROOT)}  ({len(long_df)} donor rows)")
print()
print(table.to_string(index=False))
