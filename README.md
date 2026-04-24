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

## License

Specify a `LICENSE` file before public release (e.g. MIT for code; data remain under original GEO terms).
