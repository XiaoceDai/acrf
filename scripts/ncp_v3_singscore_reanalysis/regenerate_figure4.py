#!/usr/bin/env python3
"""Orchestrate corrected Figure 4 panel rendering and final composite assembly."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts" / "ncp_v3_singscore_reanalysis"


def run_module(name: str) -> None:
    path = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    if hasattr(module, "main"):
        module.main()


def main() -> None:
    for script in (
        "fig4a_cross_sectional_duration_bridge",
        "fig4b_family_regional_duration_map",
        "fig4c_chronicity_class_summary",
        "fig4d_representative_duration_examples",
        "assemble_figure4_cross_sectional_duration",
    ):
        print(f"Running {script}...")
        run_module(script)
    print("Figure 4 cross-sectional duration composite complete.")


if __name__ == "__main__":
    main()
