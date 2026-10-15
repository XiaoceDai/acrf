# acrf

Atrial **Chamber-aware Rewiring Framework (CRF)** — analysis code and installable method core.

**GitHub:** [github.com/XiaoceDai/acrf](https://github.com/XiaoceDai/acrf)

The package version is 0.2.0. The analyses in the Frontiers submission used tag `v0.2.0` (`f10135200ca295af2d9e23c5392c3771c532ba87`). The default branch also contains a later documentation commit. That commit is not tag `v0.2.0` and does not change the scoring code.

Current submission entry point: [`reproducibility_frontiers_20261015/README.md`](reproducibility_frontiers_20261015/README.md). It maps Figure 1–3 and Supplementary Figures S1–S12, lists the 13 supplementary source tables, and records the eGene, coloc, and Heart Cell Atlas checks.

Publication-facing text uses NCP. Internal artifact names retain the historical
ncp_v3 identifier to preserve frozen lineage, checksums and
backward-compatible workflow paths.

## Layout

| Path | Purpose |
|------|---------|
| `src/crf/` | Installable **acrf** package (`import crf`); column contracts in `crf.schema` |
| `scripts/` | Cohort pipelines, figure drivers, shared `plot_style.py` |
| `scripts/v3/` | NCP figure drivers (`*_v3*.py`); thirty pathway programs / sixteen axis families |
| `tests/` | Package tests (`pytest`) |
| `results/` | TSV outputs consumed by figure scripts |
| `figures/` | PDF outputs |
| `Project_P3/` … `bulk/` | Cohort working trees (often gitignored locally) |

## Python package

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e '.[dev]'
pytest -q tests/test_crf_core.py
python -c "import crf; print(crf.__version__)"
```

Use `crf` for `aCAI` / `aCAI_rel`, Fig 2C-style trend classes, and triangulation helpers. v3 integrative panels are implemented in `scripts/v3/fig5b_v3_triangulation_heatmap.py` and `scripts/v3/fig5c_v3_cross_axis_consistency.py`.
