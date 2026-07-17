#!/usr/bin/env python3
"""
Compute GWAS vs cis-eQTL lead-variant distance for Fig 6A/6B locus summaries.

For each anchor, within the locked ±1 Mb cis window:
  - GWAS lead  = variant with minimum AF-GWAS p-value
  - eQTL lead  = variant with minimum cis-eQTL p-value for the anchor gene_id
  - delta_bp   = |GWAS lead pos − eQTL lead pos|
  - lead_*_kb  = offset from gene midpoint (hg38), for locus-zoom display

Primary data source: on-disk tabix cache slices under Project_P6/cache/
(``eqtl_<SYMBOL>_chr*`` / ``gwas_<SYMBOL>_chr*``). If cache is missing, the
script aborts with instructions to run the h4* R drivers first.

Usage
-----
  python3 scripts/compute_anchor_lead_distance.py KCNN3
  python3 scripts/compute_anchor_lead_distance.py --all-missing
  python3 scripts/compute_anchor_lead_distance.py KCNN3 --write

With ``--write``, appends/updates rows in results/fig6a_v3_locus_zoom_summary.tsv.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "Project_P6" / "cache"
OUT_TSV = ROOT / "results" / "fig6a_v3_locus_zoom_summary.tsv"

EQTL_NAMES = [
    "molecular_trait_id", "chromosome", "position", "ref", "alt", "variant",
    "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2",
    "molecular_trait_object_id", "gene_id", "median_tpm", "rsid", "nlog10p",
]
GWAS_NAMES = [
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid",
    "rs_id", "n", "hm_coordinate_conversion", "hm_code", "variant_id",
]

# Locked anchor definitions (gene body + cache window; mirrors h4d / LOCUS_WINDOWS).
ANCHORS: dict[str, dict] = {
    "PRELP":  {"gene_id": "ENSG00000188783", "chr": 1,  "start": 203475828, "end": 203491352,
               "win_start": 202475828, "win_end": 204491352},
    "OMD":    {"gene_id": "ENSG00000127083", "chr": 9,  "start": 92414245,  "end": 92424461,
               "win_start": 91414245,  "win_end": 93424461},
    "IDH2":   {"gene_id": "ENSG00000182054", "chr": 15, "start": 90083045,  "end": 90102504,
               "win_start": 89083045,  "win_end": 91102504},
    "NDUFA8": {"gene_id": "ENSG00000119421", "chr": 9,  "start": 122144058, "end": 122159819,
               "win_start": 121144058, "win_end": 123159819},
    "NTRK2":  {"gene_id": "ENSG00000148053", "chr": 9,  "start": 84668522,  "end": 85029699,
               "win_start": 83668522,  "win_end": 86029699},
    "NTF3":   {"gene_id": "ENSG00000185652", "chr": 12, "start": 5404888,   "end": 5495026,
               "win_start": 4404888,   "win_end": 6495026},
    "FRS2":   {"gene_id": "ENSG00000166225", "chr": 12, "start": 69464533,  "end": 69513650,
               "win_start": 68464533,  "win_end": 70513650},
    "BMP10":  {"gene_id": "ENSG00000163217", "chr": 2,  "start": 68869054,  "end": 68886842,
               "win_start": 67869054,  "win_end": 69886842},
    "SMAD9":  {"gene_id": "ENSG00000120693", "chr": 13, "start": 36879326,  "end": 36925087,
               "win_start": 35879326,  "win_end": 37925087},
    "SMAD6":  {"gene_id": "ENSG00000137834", "chr": 15, "start": 66703041,  "end": 66782082,
               "win_start": 65703041,  "win_end": 67782082},
    "PITX2":  {"gene_id": "ENSG00000164093", "chr": 4,  "start": 110610001, "end": 110633149,
               "win_start": 109610001, "win_end": 111633149},
    "KCNN3":  {"gene_id": "ENSG00000143603", "chr": 1,  "start": 154837802, "end": 154961963,
               "win_start": 153837802, "win_end": 155961963},
}


def _cache_paths(symbol: str, meta: dict) -> tuple[Path, Path]:
    tag = f"chr{meta['chr']}_{meta['win_start']}_{meta['win_end']}"
    eq = CACHE / f"eqtl_{symbol}_{tag}.tsv.gz"
    gw = CACHE / f"gwas_{symbol}_{tag}.tsv.gz"
    return eq, gw


def _read_slice(path: Path, names: list[str]) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty cache slice: {path.relative_to(ROOT)}")
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, sep="\t", header=None, names=names[:20])


def compute_lead_distance(symbol: str) -> dict:
    if symbol not in ANCHORS:
        raise KeyError(f"Unknown anchor {symbol!r}; known: {sorted(ANCHORS)}")
    meta = ANCHORS[symbol]
    eq_path, gw_path = _cache_paths(symbol, meta)

    eq = _read_slice(eq_path, EQTL_NAMES)
    gw = _read_slice(gw_path, GWAS_NAMES)

    gene_id = meta["gene_id"]
    eq = eq[eq["gene_id"].astype(str).str.startswith(gene_id)].copy()
    if eq.empty:
        raise ValueError(f"No eQTL rows for {gene_id} in {eq_path.name}")

    eq["position"] = pd.to_numeric(eq["position"], errors="coerce")
    eq["pvalue"] = pd.to_numeric(eq["pvalue"], errors="coerce")
    eq = eq.dropna(subset=["position", "pvalue"])
    eq_lead = eq.loc[eq["pvalue"].idxmin()]

    gw["base_pair_location"] = pd.to_numeric(gw["base_pair_location"], errors="coerce")
    gw["p_value"] = pd.to_numeric(gw["p_value"], errors="coerce")
    gw = gw.dropna(subset=["base_pair_location", "p_value"])
    if gw.empty:
        raise ValueError(f"No GWAS rows in {gw_path.name}")
    gw_lead = gw.loc[gw["p_value"].idxmin()]

    gwas_pos = int(gw_lead["base_pair_location"])
    eqtl_pos = int(eq_lead["position"])
    ref = (meta["start"] + meta["end"]) / 2.0
    delta_bp = abs(gwas_pos - eqtl_pos)

    return {
        "anchor": symbol,
        "delta_bp": float(delta_bp),
        "lead_gwas_kb": round((gwas_pos - ref) / 1000.0, 3),
        "lead_eqtl_kb": round((eqtl_pos - ref) / 1000.0, 3),
        "lead_gwas_p": float(gw_lead["p_value"]),
        "lead_eqtl_p": float(eq_lead["pvalue"]),
        "gwas_lead_pos": gwas_pos,
        "eqtl_lead_pos": eqtl_pos,
        "eqtl_cache": str(eq_path.relative_to(ROOT)),
        "gwas_cache": str(gw_path.relative_to(ROOT)),
    }


def write_summary(rows: list[dict]) -> None:
    existing = pd.read_csv(OUT_TSV, sep="\t") if OUT_TSV.is_file() else pd.DataFrame()
    new_df = pd.DataFrame([{
        "anchor": r["anchor"],
        "delta_bp": r["delta_bp"],
        "lead_gwas_kb": r["lead_gwas_kb"],
        "lead_eqtl_kb": r["lead_eqtl_kb"],
        "lead_gwas_p": r["lead_gwas_p"],
        "lead_eqtl_p": r["lead_eqtl_p"],
    } for r in rows])

    if not existing.empty:
        existing = existing[~existing["anchor"].isin(new_df["anchor"])]
        out = pd.concat([existing, new_df], ignore_index=True)
    else:
        out = new_df

    # Stable order: framework anchors then calibrators.
    order = list(ANCHORS.keys())
    out["_ord"] = out["anchor"].map({g: i for i, g in enumerate(order)})
    out = out.sort_values("_ord").drop(columns="_ord")
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("anchors", nargs="*", help="Anchor symbol(s), e.g. KCNN3")
    parser.add_argument(
        "--all-missing",
        action="store_true",
        help="Compute for anchors absent from fig6a_v3_locus_zoom_summary.tsv",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Update results/fig6a_v3_locus_zoom_summary.tsv",
    )
    args = parser.parse_args()

    targets: list[str]
    if args.all_missing:
        present = set()
        if OUT_TSV.is_file():
            present = set(pd.read_csv(OUT_TSV, sep="\t")["anchor"].astype(str))
        targets = [g for g in ANCHORS if g not in present]
    elif args.anchors:
        targets = args.anchors
    else:
        targets = ["KCNN3"]

    results: list[dict] = []
    for sym in targets:
        row = compute_lead_distance(sym)
        results.append(row)
        print(
            f"{sym}: delta_bp={row['delta_bp']:.0f} "
            f"({row['delta_bp']/1000:.1f} kb)  "
            f"lead_gwas_kb={row['lead_gwas_kb']}  lead_eqtl_kb={row['lead_eqtl_kb']}  "
            f"lead_eqtl_p={row['lead_eqtl_p']:.3e}"
        )

    if args.write:
        write_summary(results)
        print(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(results)} row(s) updated)")


if __name__ == "__main__":
    main()
