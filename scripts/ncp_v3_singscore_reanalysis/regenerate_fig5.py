#!/usr/bin/env python3
"""Regenerate corrected bulk-dependent Figure 5 assets in isolation."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
V3 = ROOT / "scripts" / "v3"
sys.path.insert(0, str(ROOT / "scripts"))


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    atlas = load_module(
        "corrected_integrated_atlas", V3 / "fig_integrated_af_rewiring_atlas.py"
    )
    atlas.ATTR_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
    atlas.DRIVER_TSV = RESULTS / "fig2d_family_effects.tsv"
    atlas.CONS_TSV = RESULTS / "fig2e_cross_cohort_consistency.tsv"
    atlas.OUT_TSV = RESULTS / "fig_integrated_af_rewiring_atlas.tsv"
    atlas.OUT_PDF = FIGURES / "fig_integrated_af_rewiring_atlas.pdf"
    atlas.LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
    atlas.main()

    fig5b = load_module(
        "corrected_fig5b", V3 / "fig5b_v3_triangulation_heatmap.py"
    )
    fig5b.ATLAS_TSV = RESULTS / "fig_integrated_af_rewiring_atlas.tsv"
    fig5b.HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
    fig5b.REPL_TSV = RESULTS / "fig2e_cross_cohort_consistency.tsv"
    fig5b.DELTA_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
    fig5b.FID_TSV = ROOT / "results" / "fig3_v3_carrier_fidelity_class.tsv"
    fig5b.CHRO_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"
    fig5b.LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
    fig5b.OUT_TSV = RESULTS / "fig5b_v3_triangulation.tsv"
    fig5b.OUT_PDF = FIGURES / "fig5b_triangulation_heatmap.pdf"
    fig5b.main()

    fig5c = load_module(
        "corrected_fig5c", V3 / "fig5c_v3_cross_axis_consistency.py"
    )
    fig5c.IN_TSV = RESULTS / "fig5b_v3_triangulation.tsv"
    fig5c.LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
    fig5c.OUT_TSV = RESULTS / "fig5c_v3_cross_axis.tsv"
    fig5c.OUT_PDF = FIGURES / "fig5c_cross_axis.pdf"
    fig5c.main()

    fig5e = load_module(
        "corrected_fig5e", V3 / "fig5e_v3_mechanism_cartoon.py"
    )
    fig5e.ATLAS_TSV = RESULTS / "fig_integrated_af_rewiring_atlas.tsv"
    fig5e.TRI_TSV = RESULTS / "fig5b_v3_triangulation.tsv"
    fig5e.CROSS_TSV = RESULTS / "fig5c_v3_cross_axis.tsv"
    fig5e.OUT_PDF = FIGURES / "fig5e_mechanism_cartoon.pdf"
    fig5e.main()
    print("Regenerated corrected Figure 5 bulk-dependent fields and Path A labels")


if __name__ == "__main__":
    main()
