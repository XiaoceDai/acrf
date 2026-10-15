# How the 15 October 2026 checks were run

These scripts were executed once and then copied here unchanged. They contain absolute paths under `/home/bio/Project_P02`. They do not download data. Running them again requires the same local files, which are not in this repository.

## eGene row count

`scripts/recount_eqtl_and_hca.py` reads `Project_P6/cache/eqtl_<GENE>_*.tsv.gz` and applies the filter in `scripts/h4_mr_coloc_anchors.R`: the gene id matches, P, beta, and standard error are finite, standard error is above 0, and minor-allele frequency is between 0 and 0.5. `n_eqtl` is that row count. Unique variant ids are fewer because some positions have more than one rsID.

The stored thresholds for 11 genes are in `results/figure6_mr_anchors_v3_consolidated.tsv` at tag v0.2.0 and equal 0.05 divided by `n_eqtl`. KCNN3’s count, 4498, is in `results/figure6_mr_positive_controls.tsv`. That table does not store a Bonferroni threshold or a Benjamini–Hochberg q. Applying 0.05/4498 to the recounted lead P does not pass. The completed comparison is `results/egene_cache_recount.tsv` and `results/egene_definition_reconstruction.tsv`.

No script in tag v0.2.0 assigns `bonf_threshold` or `bonf_passed`. The figure code reads `bonf_passed OR is_eGene_atrial`. This check did not recover the writer of those columns.

## coloc recheck

`scripts/coloc_recheck.R` reads the existing atrial eQTL cache and the matching GWAS slices. Settings, taken from the v0.2.0 driver and not adjusted after the comparison:

- priors p1 = p2 = 1e-4, p12 = 1e-5
- quantitative eQTL with sdY = 1
- case-control GWAS with case fraction 0.05
- allele-matched overlap, at least 50 shared variants
- absolute tolerance for PP0–PP4 fixed at 1e-6 before looking at the deltas

The run used coloc 5.2.3. `results/coloc_sessionInfo.txt` records the package path. `results/coloc_recheck.tsv` compares the new posteriors with the stored table. All 12 genes matched within the tolerance. PP4 remains below 0.70 for every gene. A contemporaneous log of the historical coloc version was not found, so 5.2.3 is this recheck only.

## Heart Cell Atlas file

The local matrix is `/home/bio/Project_P02/global_raw.h5ad`, 5,068,859,764 bytes. SHA256 `80917c4e67e581ee2696691ddd4b031bf45cbc7ba77757c8a203cd88c15f9c3f` was computed for this check. The file modification time is not a download date. Shape 486,134 by 33,538, 14 donors, and six regions match `reports/s1_schema.txt`. No official checksum was found, so identity with the public v1 file remains unverified.

The nuclei, LA/RA, Used = Yes, and mapped-type filter still includes all 14 donors. Paired-donor counts stay in `Project_P3/HCA/pseudobulk_chamber/paired_donors.tsv` and were not recomputed.

## Not found

- the script that first wrote `bonf_passed`
- `results/figure6_mr_anchors_v3_pipeline_wide.tsv`
- a historical coloc sessionInfo tied to the stored posteriors
- an official Heart Cell Atlas checksum or download record
- a Frontiers redraw script for Figure 3 or for the final Supplementary Figures S1–S12 PDFs
