# GCST90624411 frozen GWAS provenance

Frozen integrity package for the Yuan et al. 2025 cross-population AF GWAS
harmonised summary statistics used by Figure 6.

## Contents

| File | Purpose |
|------|---------|
| `gwas_gcst90624411_provenance.tsv` | Machine-readable key/value freeze |
| `gwas_gcst90624411_provenance.md` | Human-readable provenance + hashes |
| `ftp_md5sum_snapshot.txt` | FTP `md5sum.txt` as retrieved |
| `GCST90624411-meta.yaml` | FTP harmonised meta.yaml snapshot |
| `catalog_check_2026-07-14.tsv` | Catalog UI/REST incompleteness note |
| `download_manifest.tsv` | Local paths, URLs, sizes, hashes |
| `README.md` | This file |

## Local data files (outside this directory)

- `data/gwas/GCST90624411/harmonised/GCST90624411.h.tsv.gz`
- `data/gwas/GCST90624411/harmonised/GCST90624411.h.tsv.gz.tbi`

## Validation (2026-07-14)

- Expected FTP MD5: `04da4c8c82f9d29da802924a02a47441`
- Local MD5: `04da4c8c82f9d29da802924a02a47441` → **PASS**
- Tabix / meta.yaml MD5s match FTP snapshot → **PASS**
- Coloc rerun required → **FALSE**
