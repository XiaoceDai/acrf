"""Resolve the project root (checkout) for subprocess runs and PYTHONPATH."""

from __future__ import annotations

import os
from pathlib import Path

_ANCHOR = "FIGURE1_WORKFLOW.md"


def project_root() -> Path:
    """Return absolute path to the Project_P02 repository root.

    Resolution order:
    1. Environment variable ``PROJECT_P02_ROOT`` (recommended for CI).
    2. Walk parents from this file (editable install inside the repo).
    3. Walk parents from the current working directory.
    """
    env = os.environ.get("PROJECT_P02_ROOT")
    if env:
        p = Path(env).expanduser().resolve()
        if not (p / _ANCHOR).is_file():
            raise FileNotFoundError(
                f"PROJECT_P02_ROOT={p} does not contain {_ANCHOR}"
            )
        return p

    here = Path(__file__).resolve()
    for base in (here.parent, *here.parents):
        if (base / _ANCHOR).is_file():
            return base

    cwd = Path.cwd().resolve()
    for base in (cwd, *cwd.parents):
        if (base / _ANCHOR).is_file():
            return base

    raise FileNotFoundError(
        "Cannot locate Project P02 root (missing FIGURE1_WORKFLOW.md). "
        "Clone the repository and run from its root, or set PROJECT_P02_ROOT."
    )
