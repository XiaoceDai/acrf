#!/usr/bin/env python3
"""Compute aligned donor contrasts, healthy M[p], and family aCAIrel."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
SCORES = RESULTS / "program_scores_long.tsv"
PROGRAMS = RESULTS / "frozen_program_metadata.tsv"
FAMILIES = RESULTS / "frozen_family_metadata.tsv"
HEALTHY = "GSE112339"
TOLERANCE = 1e-12


def assertion(
    rows: list[dict[str, object]],
    assertion_id: str,
    passed: bool,
    observed: object,
    expected: object,
) -> None:
    rows.append(
        {
            "assertion_id": assertion_id,
            "passed": bool(passed),
            "observed": observed,
            "expected": expected,
        }
    )
    if not passed:
        raise AssertionError(f"{assertion_id}: observed={observed}, expected={expected}")


def build_contrasts(scores: pd.DataFrame) -> pd.DataFrame:
    key = ["dataset", "donor_id", "rhythm", "program_id"]
    pair_counts = (
        scores.groupby(key)["chamber"].agg(lambda x: tuple(sorted(x))).reset_index()
    )
    if not pair_counts["chamber"].map(lambda x: x == ("LA", "RA")).all():
        raise AssertionError("At least one scored donor-program is not exactly LA/RA")

    value_columns = ["program_score"]
    wide = scores.pivot(index=key, columns="chamber", values=value_columns)
    wide.columns = [f"score_{chamber}" for _, chamber in wide.columns]
    wide = wide.reset_index()
    metadata = scores[
        [
            "dataset",
            "program_id",
            "family_id",
            "canonical_direction",
            "direction_sign",
            "n_frozen_leading_edge_genes",
            "n_measured_genes",
            "coverage_fraction",
            "score_method",
            "gene_set_version",
        ]
    ].drop_duplicates()
    contrasts = wide.merge(
        metadata, on=["dataset", "program_id"], validate="many_to_one"
    )
    contrasts["delta"] = contrasts["score_LA"] - contrasts["score_RA"]
    contrasts["aCAI"] = contrasts["direction_sign"] * contrasts["delta"]
    contrasts["concordance"] = np.select(
        [
            contrasts["aCAI"].isna(),
            contrasts["aCAI"] > 0,
            contrasts["aCAI"] < 0,
            contrasts["aCAI"] == 0,
        ],
        ["missing", "concordant", "discordant", "neutral"],
        default="missing",
    )
    return contrasts


def main() -> None:
    scores = pd.read_csv(SCORES, sep="\t")
    programs = pd.read_csv(PROGRAMS, sep="\t")
    families = pd.read_csv(FAMILIES, sep="\t")
    assertion_rows: list[dict[str, object]] = []

    assertion(
        assertion_rows,
        "exactly_30_deployed_programs",
        programs["program_id"].nunique() == 30,
        programs["program_id"].nunique(),
        30,
    )
    assertion(
        assertion_rows,
        "exactly_16_frozen_families",
        families["family_id"].nunique() == 16,
        families["family_id"].nunique(),
        16,
    )
    direction_counts = programs["canonical_direction"].value_counts().to_dict()
    assertion(
        assertion_rows,
        "exactly_22_LA_8_RA_programs",
        direction_counts == {"LA-biased": 22, "RA-biased": 8},
        direction_counts,
        {"LA-biased": 22, "RA-biased": 8},
    )

    contrasts = build_contrasts(scores)
    healthy = contrasts[contrasts["dataset"] == HEALTHY].copy()
    n_healthy_pairs = healthy["donor_id"].nunique()
    assertion(
        assertion_rows,
        "exactly_8_complete_healthy_pairs",
        n_healthy_pairs == 8,
        n_healthy_pairs,
        8,
    )
    reference = (
        healthy.groupby("program_id", as_index=False)
        .agg(
            M=("aCAI", lambda x: float(np.mean(np.abs(x)))),
            directional_concordance_k=(
                "concordance",
                lambda x: int((x == "concordant").sum()),
            ),
            discordant_n=("concordance", lambda x: int((x == "discordant").sum())),
            neutral_n=("concordance", lambda x: int((x == "neutral").sum())),
            missing_n=("concordance", lambda x: int((x == "missing").sum())),
            n_healthy_donors=("donor_id", "nunique"),
        )
    )
    reference["reference_dataset"] = HEALTHY
    all_finite_positive = bool(
        np.isfinite(reference["M"]).all() and (reference["M"] > 0).all()
    )
    assertion(
        assertion_rows,
        "all_M_finite_and_positive",
        all_finite_positive,
        int(all_finite_positive),
        1,
    )
    contrasts = contrasts.merge(
        reference[["program_id", "M", "reference_dataset"]],
        on="program_id",
        validate="many_to_one",
    )
    contrasts["aCAIrel"] = contrasts["aCAI"] / contrasts["M"]
    contrasts["aligned_score_LA_rel"] = (
        contrasts["direction_sign"] * contrasts["score_LA"] / contrasts["M"]
    )
    contrasts["aligned_score_RA_rel"] = (
        contrasts["direction_sign"] * contrasts["score_RA"] / contrasts["M"]
    )

    assertion(
        assertion_rows,
        "no_disease_cohort_contributes_to_M",
        set(reference["reference_dataset"]) == {HEALTHY},
        sorted(reference["reference_dataset"].unique()),
        [HEALTHY],
    )
    healthy_scaled = contrasts[contrasts["dataset"] == HEALTHY]
    scale_check = healthy_scaled.groupby("program_id")["aCAIrel"].apply(
        lambda x: float(np.mean(np.abs(x)))
    )
    max_scale_error = float(np.max(np.abs(scale_check - 1.0)))
    assertion(
        assertion_rows,
        "mean_abs_healthy_aCAIrel_equals_1",
        max_scale_error <= TOLERANCE,
        max_scale_error,
        f"<= {TOLERANCE}",
    )
    component_error = float(
        np.max(
            np.abs(
                contrasts["aligned_score_LA_rel"]
                - contrasts["aligned_score_RA_rel"]
                - contrasts["aCAIrel"]
            )
        )
    )
    assertion(
        assertion_rows,
        "aligned_components_difference_equals_aCAIrel",
        component_error <= TOLERANCE,
        component_error,
        f"<= {TOLERANCE}",
    )
    swap_delta_error = float(
        np.max(np.abs((contrasts["score_RA"] - contrasts["score_LA"]) + contrasts["delta"]))
    )
    swap_acai_error = float(
        np.max(
            np.abs(
                contrasts["direction_sign"]
                * (contrasts["score_RA"] - contrasts["score_LA"])
                + contrasts["aCAI"]
            )
        )
    )
    assertion(
        assertion_rows,
        "swapping_LA_RA_reverses_delta_and_aCAI",
        max(swap_delta_error, swap_acai_error) <= TOLERANCE,
        max(swap_delta_error, swap_acai_error),
        f"<= {TOLERANCE}",
    )
    missing_scores = int(scores["program_score"].isna().sum())
    zero_scores = int((scores["program_score"] == 0).sum())
    assertion(
        assertion_rows,
        "no_missing_score_converted_to_zero",
        missing_scores == 0,
        {"missing": missing_scores, "exact_zero": zero_scores},
        "missing retained as NA; this run has zero missing scores",
    )

    family_keys = ["dataset", "donor_id", "rhythm", "family_id"]
    family = (
        contrasts.groupby(family_keys, as_index=False)
        .agg(
            family_aCAIrel=("aCAIrel", "mean"),
            family_aligned_score_LA_rel=("aligned_score_LA_rel", "mean"),
            family_aligned_score_RA_rel=("aligned_score_RA_rel", "mean"),
            n_programs_observed=("program_id", "nunique"),
        )
    )
    family_expected = contrasts.groupby(family_keys)["aCAIrel"].mean()
    family_observed = family.set_index(family_keys)["family_aCAIrel"]
    family_error = float(np.max(np.abs(family_observed - family_expected)))
    assertion(
        assertion_rows,
        "family_aCAIrel_equal_weight_program_mean",
        family_error <= TOLERANCE,
        family_error,
        f"<= {TOLERANCE}",
    )
    family = family.merge(
        families[
            [
                "family_id",
                "family_label",
                "canonical_direction",
                "canonical_family_display_order",
            ]
        ],
        on="family_id",
        validate="many_to_one",
    )

    canonical_program_order = programs.sort_values("canonical_display_order")[
        "program_id"
    ].tolist()
    displayed_program_order = programs["program_id"].tolist()
    assertion(
        assertion_rows,
        "programs_follow_canonical_display_order",
        displayed_program_order == canonical_program_order,
        len(displayed_program_order),
        "exact canonical order",
    )
    canonical_family_order = families.sort_values("canonical_family_display_order")[
        "family_id"
    ].tolist()
    displayed_family_order = families["family_id"].tolist()
    assertion(
        assertion_rows,
        "families_follow_canonical_display_order",
        displayed_family_order == canonical_family_order,
        len(displayed_family_order),
        "exact canonical order",
    )

    reference = reference.merge(
        programs[
            [
                "program_id",
                "program_order",
                "family_id",
                "canonical_direction",
                "direction_sign",
                "canonical_display_order",
                "gene_set_version",
            ]
        ],
        on="program_id",
        validate="one_to_one",
    ).sort_values("canonical_display_order")
    reference["concordance_fraction"] = (
        reference["directional_concordance_k"] / reference["n_healthy_donors"]
    )
    contrasts = contrasts.merge(
        programs[["program_id", "program_order", "canonical_display_order"]],
        on="program_id",
        validate="many_to_one",
    ).sort_values(["dataset", "donor_id", "canonical_display_order"])
    family = family.sort_values(
        ["dataset", "donor_id", "canonical_family_display_order"]
    )

    healthy_out = contrasts[contrasts["dataset"] == HEALTHY].copy()
    disease_out = contrasts[contrasts["dataset"] != HEALTHY].copy()
    healthy_out.to_csv(
        RESULTS / "healthy_program_contrasts_long.tsv", sep="\t", index=False
    )
    reference.to_csv(
        RESULTS / "healthy_reference_program_summary.tsv", sep="\t", index=False
    )
    disease_out.to_csv(
        RESULTS / "disease_program_contrasts_long.tsv", sep="\t", index=False
    )
    contrasts.to_csv(
        RESULTS / "all_program_contrasts_long.tsv", sep="\t", index=False
    )
    family.to_csv(RESULTS / "family_acairel_long.tsv", sep="\t", index=False)
    pd.DataFrame(assertion_rows).to_csv(
        RESULTS / "assertion_report.tsv", sep="\t", index=False
    )
    print(
        f"Wrote {len(healthy_out)} healthy and {len(disease_out)} disease "
        f"program contrasts; {len(family)} family rows"
    )


if __name__ == "__main__":
    main()
