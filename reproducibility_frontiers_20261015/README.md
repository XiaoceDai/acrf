# Frontiers submission materials, 15 October 2026

This directory documents the Frontiers in Bioinformatics Technology and Code submission. It does not change the acrf package. The installable version remains 0.2.0.

The analyses reported in the manuscript were run from Git tag `v0.2.0`, commit `f10135200ca295af2d9e23c5392c3771c532ba87`. The commit that adds this directory is a later documentation commit on the default branch. It is not tag `v0.2.0`.

## Install and run

From the repository root:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e '.[dev]'
pytest -q tests/test_crf_core.py
python -c "import crf; print(crf.__version__)"
```

Python 3.10 or later is required. The import name is `crf`. Program scores used by the package are calculated with singscore in R; `crf` applies the reference direction and the reference magnitude.

## Figures in the submission

| Submission file | Content | Script used for the accepted file |
|---|---|---|
| Figure 1 | Paired contrast and the healthy reference magnitude | `figure_scripts/redraw_frontiers_figure1.py` |
| Figure 2 | Cross-cohort sign agreement, correlation, and class stability | `figure_scripts/redraw_frontiers_figure2.py` |
| Figure 3 | Family effects in the three paired AF cohorts | No Frontiers redraw script was found. The accepted `Figure_3.pdf` is the baseline file. |
| Supplementary Figures S1–S8 | Cohort, selection, localization, and sensitivity displays | The accepted PDF is the baseline file. An earlier assembly script does not regenerate that PDF. |
| Supplementary Figures S9–S12 | Duration, genetics, and related panels | The accepted PDF is the baseline file. Page 4 of S9–S12 was edited in the PDF content stream; it was not redrawn by a plotting script. |

Figure 2 reads the locked benchmark tables already in `results/bmc_bioinformatics_benchmark/`. Those ten `benchmark_*.tsv` files and `phase11b1_software_environment.tsv` match the submission archive byte for byte. Figure 1 also uses `results/ncp_v3_singscore_reanalysis/migration_M_program.tsv` through the program-definition archive.

`figure_scripts/build_frontiers_docx.py` builds the manuscript DOCX from Markdown. It is not a statistical script.

## Source tables

`source_tables/` contains the 13 tables restored for Supplementary Tables S3–S7. They were copied from the accepted submission archive and were not recalculated here.

Localization classes are assigned by `scripts/v3/fig3_v3_carrier_fidelity_class.py`. When the carrier and at least one other compartment are supported, a shared AF-minus-control sign is Extension. Mixed direction is used only when the carrier is supported and every other supported compartment has the opposite sign.

## Provenance checks

`provenance/` holds the scripts and tables from the 15 October 2026 checks. Read `provenance/RUN.md` before treating those outputs as a new analysis. In short:

- The Bonferroni denominator is the count of cis-eQTL rows that pass the historical filter. It is not a unique-variant count. The script that first wrote `bonf_passed` was not found.
- Only SMAD9 of the 12 anchors passes. Local GWAS support remains 7/12 and PP4 ≥ 0.70 remains 0/12.
- coloc 5.2.3 is the version used for this recheck. It is not a recovered historical session.
- The Heart Cell Atlas SHA256 in `provenance/results/global_raw.h5ad.sha256` was computed on the local file during that check. It is not an official checksum. A download date was not recovered, and byte identity with the public v1 file was not verified.

## Not included

Cover letters, suggested reviewers, and internal author-confirmation notes are not in this repository. External expression matrices and the eQTL cache are not redistributed; the checks read those local files and store the resulting tables here.
