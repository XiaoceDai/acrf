# Project P02 — Bi-atrial chamber program rewiring in atrial fibrillation

Computational pipeline linking **healthy LA/RA program asymmetry** (`aCAI`), **AF bulk rewiring**, **single-nucleus pseudobulk carriers**, **duration associations**, **cross-layer synthesis**, and **genetic colocalization** on pre-registered anchors. Figure builds follow strict **forward-only dependencies** (`RULES.md`).

## Biology analysis design (summary)

| Stage | Question | Design highlights |
|--------|-----------|-------------------|
| **Fig 1** | How is a reproducible healthy chamber-asymmetry reference built? | Donor-level `aCAI`, frozen gene-set panel **NCP v3** (vs internal v2 benchmark), robustness checks |
| **Fig 2** | How are those programs rewired in AF in paired bi-atrial bulk? | Same LE genes as Fig 1; multi-cohort; rewiring classes + replication lollipop |
| **Fig 3** | Which LA cell classes carry AF-associated shifts? | Donor-aware **pseudobulk** (no new clustering); scope = LA arm (explicit boundary) |
| **Fig 4** | Do shifts relate to AF duration in persistent AF? | Spearman vs duration; cohort-specific caveats documented |
| **Fig 5** | Do layers triangulate to a narratable dual-axis model? | Triangulation table; multiome gate **honoured when not met** (supplementary negative panel) |
| **Fig 6** | External robustness + genetic plausibility | Meta / coloc; Plan C eGene eligibility; locus architecture; schematic synthesis |

**Strengths:** pre-registered gates, explicit negatives, orthogonal colour semantics per figure, stage lock against circular reasoning.  
**Main caveats:** bulk eQTL vs cell-type-specific regulation; modest n in some layers; genetic layer is associational (coloc PP interpretation as in frozen spec).

Full machine-readable specs: `FIGURE1_WORKFLOW.md` … `FIGURE6_WORKFLOW.md`, `INNOVATION_WORKFLOW.md`, `RULES.md`.

## Install (Python package)

See [`PACKAGING.md`](PACKAGING.md). Quick start:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
chamber-rewiring --help
```

## Data

Raw arrays and large matrices are **not** in this repository. Restore cohorts under `Project_P1/` (and related paths) per each figure workflow (GEO accessions listed there). Optional bulk mirror layout is described in `FIGURE2_WORKFLOW.md` (`bulk/Project_P01/`, …).

## Reproducibility

- Run order and deliverables are frozen in the `FIGURE*_WORKFLOW.md` files.
- Generated tables and figures: `results/`, `figures/` (included in git unless you uncomment those lines in `.gitignore` for a slimmer mirror).

## Third-party tools

- **LDSC**: the `tools/ldsc/` tree is **not** tracked (large upstream test data). Clone [bulik/ldsc](https://github.com/bulik/ldsc) into `tools/ldsc` if you run stratified LD-score steps. Reference panels under `ldsc_data/` may be empty placeholders; download per your `_ldsc_download.sh` / Zenodo mirrors.

## GitHub

Repository was initialised with `git`; large paths are listed in `.gitignore` (see below). To publish:

```bash
# on GitHub: create an empty repo (no README), then:
git remote add origin https://github.com/<you>/<repo>.git
git push -u origin main
```

## License

See [`LICENSE`](LICENSE) (MIT for code in this repo). **GEO / eQTL Catalogue / GWAS** data stay under their original terms; this repository does not redistribute raw cohort matrices.

### What is excluded from git (must be restored locally)

`Project_P1`–`P6`, `bulk/`, `global_raw.h5ad`, root-level `QTD*` / `GCST*` caches, `tools/ldsc/`. Paths match your frozen workflows and `PACKAGING.md`.
