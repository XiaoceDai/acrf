# acrf

Atrial **Chamber-aware Rewiring Framework (CRF)** — analysis code and installable method core.

**GitHub:** [github.com/XiaoceDai/acrf](https://github.com/XiaoceDai/acrf)

## Layout

| Path | Purpose |
|------|---------|
| `src/crf/` | Installable **acrf** package (`import crf`); column contracts in `crf.schema` |
| `scripts/` | Cohort pipelines, figure drivers, shared `plot_style.py` |
| `scripts/v3/` | NCP_v3 figure drivers (`*_v3*.py`) |
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
