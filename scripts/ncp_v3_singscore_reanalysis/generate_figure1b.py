#!/usr/bin/env python3
"""Generate final corrected Figure 1B from the approved reanalysis tables."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"


def load_renderer():
    path = ROOT / "scripts" / "fig1b_baseline_effect_concordance.py"
    spec = importlib.util.spec_from_file_location("corrected_fig1b_renderer", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    contrasts = pd.read_csv(
        RESULTS / "healthy_program_contrasts_long.tsv", sep="\t"
    )
    reference = pd.read_csv(
        RESULTS / "healthy_reference_program_summary.tsv", sep="\t"
    )
    programs = pd.read_csv(RESULTS / "frozen_program_metadata.tsv", sep="\t")
    families = pd.read_csv(RESULTS / "frozen_family_metadata.tsv", sep="\t")
    renderer = load_renderer()

    donor_order = sorted(contrasts["donor_id"].astype(str).unique())
    if len(donor_order) != 8:
        raise AssertionError(f"Figure 1B requires 8 healthy donors, found {len(donor_order)}")
    donor_matrix = contrasts.pivot(
        index="program_id", columns="donor_id", values="aCAI"
    ).reindex(columns=donor_order)
    aligned_cols = [f"aligned_aCAI_{donor}" for donor in donor_order]
    donor_matrix.columns = aligned_cols
    donor_matrix = donor_matrix.reset_index()

    df = (
        reference.merge(
            programs[
                [
                    "program_id", "program_order", "family_id",
                    "canonical_direction", "canonical_display_order",
                ]
            ],
            on=["program_id", "program_order", "family_id", "canonical_direction", "canonical_display_order"],
            validate="one_to_one",
        )
        .merge(
            families[["family_id", "short_label"]],
            on="family_id",
            validate="many_to_one",
        )
        .merge(donor_matrix, on="program_id", validate="one_to_one")
        .sort_values("canonical_display_order")
        .reset_index(drop=True)
    )
    df["pathway"] = df["program_id"]
    df["direction"] = df["canonical_direction"]
    df["display_order"] = df["canonical_display_order"]
    df["display_label"] = df["program_id"].map(
        lambda value: renderer.short_pathway_label(value, 35)
    )

    summary = {
        "n_programs": int(len(df)),
        "n_la": int((df["direction"] == "LA-biased").sum()),
        "n_ra": int((df["direction"] == "RA-biased").sum()),
        "n_full_concordance": int((df["directional_concordance_k"] == 8).sum()),
        "n_ge6_concordance": int((df["directional_concordance_k"] >= 6).sum()),
        "median_M": float(df["M"].median()),
        "n_families": int(df["family_id"].nunique()),
    }
    if summary["n_programs"] != 30 or summary["n_families"] != 16:
        raise AssertionError(summary)

    FIGURES.mkdir(parents=True, exist_ok=True)
    renderer.OUT_BASE = FIGURES / "fig1b_healthy_reference_concordance"
    renderer.OUT_SOURCE = RESULTS / "fig1b_healthy_reference_concordance_source.tsv"
    source_columns = [
        "display_order", "program_order", "program_id", "display_label",
        "direction", "family_id", "short_label", "M",
        "directional_concordance_k", *aligned_cols,
    ]
    df[source_columns].to_csv(renderer.OUT_SOURCE, sep="\t", index=False)
    renderer.render(df, aligned_cols, summary)
    pd.DataFrame([summary]).to_csv(
        RESULTS / "fig1b_healthy_reference_concordance_summary.tsv",
        sep="\t",
        index=False,
    )
    if not np.isclose(
        summary["median_M"], reference["M"].median(), rtol=0, atol=1e-15
    ):
        raise AssertionError("Displayed median M does not match corrected source")
    print(
        "Figure 1B summary: "
        f"{summary['n_la']} LA / {summary['n_ra']} RA; "
        f"8/8={summary['n_full_concordance']}/{summary['n_programs']}; "
        f">=6/8={summary['n_ge6_concordance']}/{summary['n_programs']}; "
        f"median M={summary['median_M']:.10g}"
    )


if __name__ == "__main__":
    main()
