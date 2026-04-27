# acrf

Atrial **Chamber-aware Rewiring Framework (CRF)** analysis repository.

**GitHub:** [github.com/XiaoceDai/acrf](https://github.com/XiaoceDai/acrf)

## Repository layout

| Path | Purpose |
|------|---------|
| **`docs/`** | Specs + audit (**`docs/INDEX.md`**); manuscript and submission materials are **local-only** (not in this public repo) |
| `RESULTS_OUTLINE_V3.md` | Symlink → **`docs/specs/`**; Results ↔ `figures/*_v3*.pdf` |
| `src/crf/` | Installable **acrf** method core (`import crf`) |
| `scripts/` | Pipelines, legacy figures, helpers (`plot_style.py`) |
| **`scripts/v3/`** | **NCP_v3** primary figure drivers |
| `results/` | TSV outputs consumed by figures and manuscript |
| `figures/` | PDFs and TeX compositing |
| `FIGURE*_WORKFLOW.md` | Symlinks → **`docs/specs/`**；**从哪读起：** **`WORKFLOWS_INDEX.md`** |
| `Project_P3/`, `Project_P4/`, `Project_P5/`, `bulk/` | Cohort-specific working trees |

## Python package (method core)

The installable distribution **acrf** (*Atrial Chamber Rewiring Framework*; `pip install acrf` or editable `-e .`; import **`crf`**) exposes portable definitions, scales,
and classification rules used in the paper pipeline:

```bash
python3 -m venv .venv && source .venv/bin/activate   # if pip refuses system-wide (PEP 668)
pip install -e '.[dev]'
pytest -q tests/test_crf_core.py
python -c "import crf; print(crf.__version__)"
```

Import `crf` for `aCAI` / `aCAI_rel`, Fig 2C trend classes, and legacy
ten-program triangulation helpers. **Fig. 5B–5C v3** integration is implemented
in `scripts/v3/fig5b_v3_triangulation_heatmap.py` and
`scripts/v3/fig5c_v3_cross_axis_consistency.py`. Column name contracts live in
`crf.schema`.

Project-specific figure scripts remain under `scripts/`; they can gradually
switch to `import crf` for duplicated logic.
