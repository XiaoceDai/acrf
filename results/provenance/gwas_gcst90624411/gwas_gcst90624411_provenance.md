# GWAS provenance — GCST90624411 (Figure 6)

## Authority

| Item | Value |
|------|-------|
| Publication | Yuan et al. *Nat Commun* 2025;16:6426 |
| DOI | 10.1038/s41467-025-61720-2 |
| Author correction | 10.1038/s41467-025-64765-5 |
| PubMed | 40645996 |
| Accession | **GCST90624411** (cross-population AF GWAS) |
| Locked sample size | **252,438** AF cases / **1,959,739** controls |

The author correction replaced the erroneous **168,007**-case figure with **252,438** and corrected a European-discovery case count separately. Manuscript wording uses the corrected N.

## File frozen for Figure 6

| Role | Path |
|------|------|
| Harmonised sumstats | `data/gwas/GCST90624411/harmonised/GCST90624411.h.tsv.gz` |
| Tabix index | `data/gwas/GCST90624411/harmonised/GCST90624411.h.tsv.gz.tbi` |
| FTP meta.yaml (snapshot) | `results/provenance/gwas_gcst90624411/GCST90624411-meta.yaml` |

Remote (EBI FTP harmonised directory):

`https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/GCST90624411/harmonised/`

## Integrity

| File | Bytes | MD5 | SHA256 |
|------|------:|-----|--------|
| `GCST90624411.h.tsv.gz` | 1,720,310,142 | `04da4c8c82f9d29da802924a02a47441` | `d883b43c524ec179232d055435071d4125f2bb1abbeffa19d3bb3d924d2ddeeb` |
| `GCST90624411.h.tsv.gz.tbi` | 1,965,432 | `a77b28bae2c9a307db45cc95c3f7729e` | `6616a7db24b18c28ac3b9904f7244254754ae33e3397052476819679f7b1f97b` |
| `GCST90624411-meta.yaml` | 1,218 | `ab9314bcfa9ffdbd6162823f9b96639f` | `a59ef12e5a53c2aeab8ebdaf6a23741bc620ac3477a643c16aec4afbd39d156f` |

**MD5 assert:** local sumstats MD5 **==** FTP `md5sum.txt` expected (`04da4c8c…`) → **PASS**

- `.tbi` and meta.yaml MD5s also match the FTP `md5sum.txt` snapshot (`ftp_md5sum_snapshot.txt`).
- Download UTC (sumstats complete): **2026-07-14T11:16:22Z**

## Catalog note (provenance only — not main manuscript)

> The GWAS Catalog study page displayed incomplete publication/sample metadata and a "not yet associated" banner on 2026-07-14. Publication attribution, accession role and sample size were therefore taken from the corrected article and its Data availability statement; the Catalog accession and FTP identify the deposited summary-statistics file.

Details: `catalog_check_2026-07-14.tsv`.

## Coloc / Figure 6 rerun gate

**coloc_rerun_required = FALSE**

Do **not** rerun Figure 6 / coloc unless:

1. local harmonised-file MD5 fails,
2. the actual accession differs, or
3. the tabix file does not match the frozen harmonised file.

Analysis scripts continue to use the remote tabix URL for the same accession/file; local freeze is the integrity anchor.

## Manuscript wording (locked)

> AF GWAS summary statistics were obtained from the harmonised NHGRI–EBI GWAS Catalog release GCST90624411, identified by Yuan et al. as the cross-population AF GWAS. The corrected publication reports 252,438 AF cases and 1,959,739 controls.

## Figure 6A visible source label (locked)

> Yuan et al. 2025 cross-population AF GWAS · GCST90624411

Do not display Miyazawa, Roselli-as-GWAS-source, internal `mr_` filenames, or checksums on the panel.
