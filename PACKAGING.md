# Installable workflow (Project P02)

This repository can be installed as a small Python package so that **PLOS Computational Biology** (and other journals) can cite a single install step and a stable CLI.

## Prerequisites

- Python 3.10+ (use a **virtual environment** if your OS uses PEP 668 / “externally managed” system Python).
- Full reproduction also requires **R** and external data under `bulk/`, `Project_P1/`, etc., as documented in each `FIGURE*_WORKFLOW.md`.

## Editable install (recommended)

From the repository root (the directory that contains `FIGURE1_WORKFLOW.md`):

```bash
pip install -e .
```

Optional extras:

```bash
pip install -e ".[full]"    # scikit-learn, statsmodels (some scripts)
pip install -e ".[dev]"    # pytest, ruff
```

## CLI

- **Innovation release bundle** (manifest + core Fig 2 v3 pipeline steps):

  ```bash
  chamber-rewiring release-bundle
  ```

- **Run any script** with `PYTHONPATH` set to the repo root:

  ```bash
  chamber-rewiring run scripts/fig6b_v3_anchor_matrix.py
  ```

If the checkout is not discoverable (e.g. non-standard layout), set:

```bash
export PROJECT_P02_ROOT=/path/to/Project_P02
```

## Scope

The package exposes **entry points and path resolution**; the scientific code remains in `scripts/` to avoid a disruptive refactor. For publication, archive a **tagged commit** (e.g. Zenodo DOI) together with this `pyproject.toml`.
