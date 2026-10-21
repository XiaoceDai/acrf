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

## Historical documentation

Commit `4b62226dbee44052eb649bf854773ef30725b76c` replaced five root documents with symlinks into `docs/specs/` and `docs/reports/`. Those target files were not committed, so the symlinks did not resolve. The five paths below are ordinary Markdown entry points again. Each one links to its own original text in commit `1f21c505c8448dcf0dccd5f396233b762381aff6`.

These documents are early development notes from before v0.2.0. Restoring the reference does not re-validate the audit conclusions in the original text. Figure numbers and method descriptions in those files may not apply to the current submission. They are not missing files.

| Entry | Original |
|---|---|
| [RULES.md](RULES.md) | https://github.com/XiaoceDai/acrf/blob/1f21c505c8448dcf0dccd5f396233b762381aff6/RULES.md |
| [ILLUSTRATOR_MUSTFIX.md](ILLUSTRATOR_MUSTFIX.md) | https://github.com/XiaoceDai/acrf/blob/1f21c505c8448dcf0dccd5f396233b762381aff6/ILLUSTRATOR_MUSTFIX.md |
| [INNOVATION_WORKFLOW.md](INNOVATION_WORKFLOW.md) | https://github.com/XiaoceDai/acrf/blob/1f21c505c8448dcf0dccd5f396233b762381aff6/INNOVATION_WORKFLOW.md |
| [METHODOLOGY_AUDIT_REPORT.md](METHODOLOGY_AUDIT_REPORT.md) | https://github.com/XiaoceDai/acrf/blob/1f21c505c8448dcf0dccd5f396233b762381aff6/METHODOLOGY_AUDIT_REPORT.md |
| [METHOD_VALIDITY_P5_REPORT.md](METHOD_VALIDITY_P5_REPORT.md) | https://github.com/XiaoceDai/acrf/blob/1f21c505c8448dcf0dccd5f396233b762381aff6/METHOD_VALIDITY_P5_REPORT.md |

Current submission materials are described in [`reproducibility_frontiers_20261015/README.md`](reproducibility_frontiers_20261015/README.md). Items that are still missing are listed in [`reproducibility_frontiers_20261015/MISSING.md`](reproducibility_frontiers_20261015/MISSING.md).
