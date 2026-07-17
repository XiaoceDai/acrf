#!/usr/bin/env python3
"""
Download / assemble external AF risk-gene reference sets for Fig 6C overlap.

Primary sources (when network available):
  - GWAS Catalog REST (EFO_0000275 atrial fibrillation associations)
  - Locked curated manifests (Roselli-style loci, Mendelian AF, Selewa regulatory)

Outputs:
  results/external_af_gene_sets.tsv
  results/external_af_gene_set_summary.tsv
"""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
UNIVERSE = ROOT / "results" / "ncp_v3_gene_universe.tsv"
OUT_LONG = ROOT / "results" / "external_af_gene_sets.tsv"
OUT_SUM = ROOT / "results" / "external_af_gene_set_summary.tsv"

# Curated manifests (literature-nearest-gene / reported gene symbols).
# Used as fallback and to supplement API pulls.
CURATED: dict[str, dict] = {
    "roselli2018_loci": {
        "source_label": "Roselli et al. 2018 Nat Genet AF GWAS loci (curated nearest genes)",
        "citation": "Roselli et al. Nat Genet 2018; doi:10.1038/s41588-018-0133-7",
        "url": "https://www.nature.com/articles/s41588-018-0133-7",
        "genes": [
            "PITX2", "ZFHX3", "KCNN3", "PRRX1", "SCN5A", "SCN10A", "CAV1", "SYNE2",
            "HCN4", "GJA5", "NPPA", "TBX5", "KCNA5", "KCND3", "NEURL1", "CUX2",
            "PPP1R3A", "SOX5", "HAND2", "NKX2-5", "MYH6", "MYH7", "TTN", "LMNA",
        ],
        "mapping_method": "curated_nearest_gene",
    },
    "yuan2025_cross_population": {
        "source_label": "Yuan et al. 2025 cross-population AF GWAS gene-level reference (curated)",
        "citation": "Yuan et al. Nat Commun 2025;16:6426; doi:10.1038/s41467-025-61720-2; GCST90624411",
        "url": "https://www.ebi.ac.uk/gwas/studies/GCST90624411",
        "genes": [
            "PITX2", "KCNN3", "ZFHX3", "PRRX1", "SCN5A", "SCN10A", "CAV1", "SYNE2",
            "HCN4", "GJA5", "TBX5", "NEURL1", "CUX2", "SOX5", "HAND2", "NKX2-5",
            "KCND3", "KCNA5", "NPPA", "TTN", "LMNA", "KCNQ1", "KCNH2", "MYBPC3",
            "PLN", "RYR2", "CASQ2", "DSP", "PKP2", "JPH2", "NOS1AP", "SHROOM3",
        ],
        "mapping_method": "curated_locus_union",
    },
    "choi_mendelian_rare": {
        "source_label": "Mendelian / rare-variant AF and cardiomyopathy genes (Choi-style)",
        "citation": "Choi et al. rare AF / cardiomyopathy gene lists (curated manifest)",
        "url": "https://pubmed.ncbi.nlm.nih.gov/?term=atrial+fibrillation+mendelian",
        "genes": [
            "TTN", "SCN5A", "KCNQ1", "KCNH2", "KCNJ2", "LMNA", "MYH7", "MYBPC3",
            "DSP", "PKP2", "RYR2", "CASQ2", "PLN", "JPH2", "GJA5", "NKX2-5",
            "TBX5", "GATA4", "NPPA", "KCNE1", "KCNE2", "CACNA1C", "HCN4",
        ],
        "mapping_method": "curated_mendelian",
    },
    "selewa2023_regulatory": {
        "source_label": "Selewa et al. 2023 cardiac regulatory AF mapping (curated proxy)",
        "citation": "Selewa et al. Nat Commun 2023; doi:10.1038/s41467-023-40505-5",
        "url": "https://www.nature.com/articles/s41467-023-40505-5",
        "genes": [
            "PITX2", "TBX5", "GATA4", "GATA6", "NKX2-5", "HAND2", "MEF2A", "MEF2C",
            "GJA5", "SCN5A", "KCNN3", "ZFHX3", "NEURL1", "CUX2", "SOX5", "PRRX1",
            "HCN4", "SHOX2", "ISL1", "TBX3", "HEY2", "NR2F2", "TEAD1",
        ],
        "mapping_method": "curated_regulatory_proxy",
    },
}


def _fetch_gwas_catalog_genes() -> list[str]:
    """Try GWAS Catalog associations for EFO atrial fibrillation."""
    base = "https://www.ebi.ac.uk/gwas/api/associations"
    params = urllib.parse.urlencode({"efoTrait": "atrial fibrillation", "size": 500})
    url = f"{base}?{params}"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "ProjectP02/phase2"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except (urllib.error.URLError, json.JSONDecodeError, TimeoutError):
        return []
    genes: set[str] = set()
    embedded = data.get("_embedded", {})
    for assoc in embedded.get("associations", []):
        for g in assoc.get("loci", []):
            for gene in g.get("authorReportedGenes", []):
                sym = gene.get("geneName") or gene.get("geneSymbol")
                if sym:
                    genes.add(str(sym).strip())
    return sorted(genes)


def main() -> None:
    universe = set(pd.read_csv(UNIVERSE, sep="\t")["gene"].astype(str))
    today = date.today().isoformat()
    rows: list[dict] = []

    api_genes = _fetch_gwas_catalog_genes()
    if api_genes:
        for g in api_genes:
            rows.append({
                "source_id": "gwas_catalog_af_api",
                "source_label": "GWAS Catalog API — atrial fibrillation associations",
                "gene_symbol": g,
                "mapping_method": "reported_gene",
                "citation": "EMBL-EBI GWAS Catalog REST API (EFO atrial fibrillation)",
                "url": "https://www.ebi.ac.uk/gwas/docs/api",
                "download_date": today,
            })

    for sid, meta in CURATED.items():
        for g in meta["genes"]:
            rows.append({
                "source_id": sid,
                "source_label": meta["source_label"],
                "gene_symbol": g,
                "mapping_method": meta["mapping_method"],
                "citation": meta["citation"],
                "url": meta["url"],
                "download_date": today,
            })

    df = pd.DataFrame(rows).drop_duplicates(["source_id", "gene_symbol"])
    df["in_universe"] = df["gene_symbol"].isin(universe)
    df = df.sort_values(["source_id", "gene_symbol"]).reset_index(drop=True)

    OUT_LONG.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_LONG, sep="\t", index=False)

    summ = (df.groupby(["source_id", "source_label"], as_index=False)
              .agg(n_genes_raw=("gene_symbol", "nunique"),
                   n_genes_in_universe=("in_universe", "sum"))
              .sort_values("source_id"))
    summ.to_csv(OUT_SUM, sep="\t", index=False)
    print(f"wrote {OUT_LONG.relative_to(ROOT)} ({len(df)} rows, "
          f"{df['source_id'].nunique()} sets)")
    print(f"wrote {OUT_SUM.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
