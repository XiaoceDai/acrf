#!/usr/bin/env python3
"""Regenerate corrected Figure 2 source tables and namespaced panels."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from scipy.stats import chi2, ttest_ind


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
PROGRAMS = pd.read_csv(RESULTS / "frozen_program_metadata.tsv", sep="\t")
FAMILIES = pd.read_csv(RESULTS / "frozen_family_metadata.tsv", sep="\t")
PROGRAM = pd.read_csv(RESULTS / "disease_program_contrasts_long.tsv", sep="\t")
FAMILY = pd.read_csv(RESULTS / "family_acairel_long.tsv", sep="\t")
COHORTS = ["GSE128188", "GSE79768", "GSE115574"]
STATE_ORDER = ["Preserved", "Attenuated", "Reversed", "Amplified"]
STATE_COLORS = ["#A9C9B5", "#E9C27B", "#C6A2D7", "#E59595"]
EPS_REVERSED = 0.02
EPS_FLOOR = 0.02
RATIO_ATTENUATED = 0.67
RATIO_AMPLIFIED = 1.33


def hedges_g(af: np.ndarray, sr: np.ndarray) -> tuple[float, float]:
    nx, ny = len(af), len(sr)
    sx, sy = np.std(af, ddof=1), np.std(sr, ddof=1)
    pooled_sd = np.sqrt(((nx - 1) * sx**2 + (ny - 1) * sy**2) / (nx + ny - 2))
    if not np.isfinite(pooled_sd) or pooled_sd <= 0:
        return np.nan, np.nan
    correction = 1 - 3 / (4 * (nx + ny - 2) - 1)
    g = correction * (np.mean(af) - np.mean(sr)) / pooled_sd
    variance = (nx + ny) / (nx * ny) + g**2 / (2 * (nx + ny))
    return float(g), float(variance)


def classify_state(af: float, reference: float) -> str:
    if af <= -EPS_REVERSED:
        return "Reversed"
    ratio = af / max(reference, EPS_FLOOR)
    if ratio < RATIO_ATTENUATED:
        return "Attenuated"
    if ratio <= RATIO_AMPLIFIED:
        return "Preserved"
    return "Amplified"


def driver(delta_la: float, delta_ra: float) -> str:
    if abs(delta_la) >= 1.25 * abs(delta_ra):
        return "LA-driven"
    if abs(delta_ra) >= 1.25 * abs(delta_la):
        return "RA-driven"
    return "Bilateral"


def quadrant(delta_la: float, delta_ra: float) -> str:
    if delta_la >= 0 and delta_ra >= 0:
        return "Q1 (+,+)"
    if delta_la < 0 <= delta_ra:
        return "Q2 (-,+)"
    if delta_la < 0 and delta_ra < 0:
        return "Q3 (-,-)"
    return "Q4 (+,-)"


def save_all(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg"):
        fig.savefig(FIGURES / f"{stem}.{suffix}", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def donor_heatmap() -> None:
    out = PROGRAM[
        [
            "dataset", "donor_id", "rhythm", "family_id", "program_id",
            "program_order", "canonical_direction", "aCAI", "M", "aCAIrel",
            "n_measured_genes", "coverage_fraction", "canonical_display_order",
        ]
    ].copy()
    out.to_csv(RESULTS / "fig2b_program_donor_acairel.tsv", sep="\t", index=False)
    FAMILY.to_csv(RESULTS / "fig2b_family_donor_acairel.tsv", sep="\t", index=False)
    out["sample_key"] = out["dataset"] + ":" + out["donor_id"].astype(str)
    sample_meta = (
        out[["sample_key", "dataset", "rhythm", "donor_id"]]
        .drop_duplicates()
        .assign(
            cohort_rank=lambda x: x["dataset"].map({c: i for i, c in enumerate(COHORTS)}),
            rhythm_rank=lambda x: x["rhythm"].map({"SR": 0, "AF": 1}),
        )
        .sort_values(["cohort_rank", "rhythm_rank", "donor_id"])
    )
    program_order = PROGRAMS.sort_values("canonical_display_order")["program_id"]
    matrix = out.pivot(index="program_id", columns="sample_key", values="aCAIrel")
    matrix = matrix.reindex(index=program_order, columns=sample_meta["sample_key"])
    vmax = max(1.0, float(np.nanpercentile(np.abs(matrix.to_numpy()), 97.5)))
    fig, ax = plt.subplots(figsize=(12, 8))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(len(matrix)))
    ax.set_yticklabels([f"{i:02d}" for i in range(1, 31)], fontsize=6)
    ax.set_xticks(range(len(sample_meta)))
    ax.set_xticklabels(sample_meta["rhythm"], rotation=90, fontsize=5)
    ax.set_title("Corrected donor-level aligned aCAIrel (30 frozen programs)", loc="left")
    ax.set_xlabel("Donors grouped by cohort and rhythm")
    ax.set_ylabel("Canonical program display order")
    fig.colorbar(image, ax=ax, label="aligned aCAIrel")
    save_all(fig, "fig2b_donor_acairel_heatmap")


def states_table() -> pd.DataFrame:
    pooled_sr = (
        FAMILY[FAMILY["rhythm"] == "SR"].groupby("family_id")["family_aCAIrel"].mean()
    )
    rows = []
    for (family_id, dataset), group in FAMILY.groupby(["family_id", "dataset"]):
        af = group.loc[group["rhythm"] == "AF", "family_aCAIrel"]
        sr = group.loc[group["rhythm"] == "SR", "family_aCAIrel"]
        if af.empty or sr.empty:
            continue
        af_mean, sr_mean = float(af.mean()), float(sr.mean())
        reference = sr_mean
        source = "cohort_SR"
        if not np.isfinite(reference) or reference <= EPS_FLOOR:
            reference = float(pooled_sr.loc[family_id])
            source = "pooled_SR"
        if not np.isfinite(reference) or reference <= EPS_FLOOR:
            reference = EPS_FLOOR
            source = "floor"
        rows.append(
            {
                "family_id": int(family_id),
                "cohort": dataset,
                "direction": group["canonical_direction"].iloc[0],
                "n_SR": int(len(sr)),
                "n_AF": int(len(af)),
                "mean_aCAIrel_SR": sr_mean,
                "mean_aCAIrel_AF": af_mean,
                "reference_aCAIrel_SR": reference,
                "reference_source": source,
                "af_over_ref": af_mean / reference,
                "state": classify_state(af_mean, reference),
            }
        )
    states = pd.DataFrame(rows).merge(
        FAMILIES[["family_id", "family_label", "short_label", "canonical_family_display_order"]],
        on="family_id",
        validate="many_to_one",
    )
    states = states.sort_values(["canonical_family_display_order", "cohort"])
    states.to_csv(RESULTS / "fig2c_family_cohort_states.tsv", sep="\t", index=False)
    return states


def plot_states(states: pd.DataFrame) -> None:
    family_ids = FAMILIES.sort_values("canonical_family_display_order")["family_id"]
    matrix = states.pivot(index="family_id", columns="cohort", values="state").reindex(
        index=family_ids, columns=COHORTS
    )
    codes = matrix.apply(
        lambda column: column.map({state: i for i, state in enumerate(STATE_ORDER)})
    ).to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(6.5, 7))
    ax.imshow(codes, aspect="auto", cmap=ListedColormap(STATE_COLORS), vmin=-0.5, vmax=3.5)
    ax.set_xticks(range(3), [c.replace("GSE", "") for c in COHORTS])
    labels = FAMILIES.set_index("family_id").loc[family_ids, "short_label"]
    ax.set_yticks(range(16), [f"F{fid} {label}" for fid, label in zip(family_ids, labels)], fontsize=7)
    for y in range(16):
        for x in range(3):
            ax.text(x, y, STATE_ORDER[int(codes[y, x])][:3], ha="center", va="center", fontsize=6)
    ax.set_title("Corrected family-by-cohort rewiring states", loc="left")
    save_all(fig, "fig2c_family_cohort_states")


def effects_table() -> pd.DataFrame:
    rows = []
    for family_id in FAMILIES["family_id"]:
        family_group = FAMILY[FAMILY["family_id"] == family_id]
        cohort_rows = []
        for cohort in COHORTS:
            group = family_group[family_group["dataset"] == cohort]
            af = group[group["rhythm"] == "AF"]
            sr = group[group["rhythm"] == "SR"]
            g, variance = hedges_g(
                af["family_aCAIrel"].to_numpy(), sr["family_aCAIrel"].to_numpy()
            )
            delta_la = float(
                af["family_aligned_score_LA_rel"].mean()
                - sr["family_aligned_score_LA_rel"].mean()
            )
            delta_ra = float(
                af["family_aligned_score_RA_rel"].mean()
                - sr["family_aligned_score_RA_rel"].mean()
            )
            p = float(
                ttest_ind(
                    af["family_aCAIrel"], sr["family_aCAIrel"], equal_var=False
                ).pvalue
            )
            record = {
                "family_id": int(family_id),
                "cohort": cohort,
                "hedges_g": g,
                "ci_low": g - 1.96 * np.sqrt(variance),
                "ci_high": g + 1.96 * np.sqrt(variance),
                "p_value": p,
                "n_AF": len(af),
                "n_SR": len(sr),
                "direction": family_group["canonical_direction"].iloc[0],
                "delta_LA": delta_la,
                "delta_RA": delta_ra,
                "delta_aCAI": delta_la - delta_ra,
                "abs_delta_aCAI": abs(delta_la - delta_ra),
                "dist_to_y_eq_x": abs(delta_ra - delta_la) / np.sqrt(2),
                "dist_to_y_eq_negx": abs(delta_ra + delta_la) / np.sqrt(2),
                "quadrant": quadrant(delta_la, delta_ra),
                "driver": driver(delta_la, delta_ra),
                "g": g,
                "var_g": variance,
                "ci_lo": g - 1.96 * np.sqrt(variance),
                "ci_hi": g + 1.96 * np.sqrt(variance),
            }
            rows.append(record)
            cohort_rows.append(record)
        valid = [row for row in cohort_rows if np.isfinite(row["var_g"]) and row["var_g"] > 0]
        weights = np.array([1 / row["var_g"] for row in valid])
        pooled_g = float(np.sum(weights * [row["g"] for row in valid]) / weights.sum())
        pooled_var = float(1 / weights.sum())
        af_all = family_group[family_group["rhythm"] == "AF"]
        sr_all = family_group[family_group["rhythm"] == "SR"]
        delta_la = float(
            af_all["family_aligned_score_LA_rel"].mean()
            - sr_all["family_aligned_score_LA_rel"].mean()
        )
        delta_ra = float(
            af_all["family_aligned_score_RA_rel"].mean()
            - sr_all["family_aligned_score_RA_rel"].mean()
        )
        rows.append(
            {
                **valid[0],
                "cohort": "POOLED",
                "hedges_g": pooled_g,
                "ci_low": pooled_g - 1.96 * np.sqrt(pooled_var),
                "ci_high": pooled_g + 1.96 * np.sqrt(pooled_var),
                "p_value": float(
                    ttest_ind(
                        af_all["family_aCAIrel"],
                        sr_all["family_aCAIrel"],
                        equal_var=False,
                    ).pvalue
                ),
                "n_AF": len(af_all),
                "n_SR": len(sr_all),
                "delta_LA": delta_la,
                "delta_RA": delta_ra,
                "delta_aCAI": delta_la - delta_ra,
                "abs_delta_aCAI": abs(delta_la - delta_ra),
                "dist_to_y_eq_x": abs(delta_ra - delta_la) / np.sqrt(2),
                "dist_to_y_eq_negx": abs(delta_ra + delta_la) / np.sqrt(2),
                "quadrant": quadrant(delta_la, delta_ra),
                "driver": driver(delta_la, delta_ra),
                "g": pooled_g,
                "var_g": pooled_var,
                "ci_lo": pooled_g - 1.96 * np.sqrt(pooled_var),
                "ci_hi": pooled_g + 1.96 * np.sqrt(pooled_var),
            }
        )
    effects = pd.DataFrame(rows).merge(
        FAMILIES[["family_id", "family_label", "short_label", "canonical_family_display_order"]],
        on="family_id",
        validate="many_to_one",
    )
    effects["short"] = effects["short_label"]
    effects["estimate_type"] = effects["cohort"].apply(
        lambda x: "pooled" if x == "POOLED" else "cohort-specific"
    )
    effects["pooling_model"] = effects["cohort"].apply(
        lambda x: "fixed-effect_inverse_variance" if x == "POOLED" else "NA"
    )
    effects["tau2"] = "NA"
    effects["sampling_variance"] = effects["var_g"]
    effects["standard_error"] = np.sqrt(effects["var_g"])
    effects["ci_lower"] = effects["ci_low"]
    effects["ci_upper"] = effects["ci_high"]
    effects = effects.sort_values(["canonical_family_display_order", "cohort"])
    effects.to_csv(RESULTS / "fig2d_family_effects.tsv", sep="\t", index=False)
    return effects


def plot_effects(effects: pd.DataFrame) -> None:
    pooled = effects[effects["cohort"] == "POOLED"].sort_values("canonical_family_display_order")
    y = np.arange(len(pooled))
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.axvline(0, color="#888", lw=0.8)
    ax.errorbar(
        pooled["g"], y,
        xerr=[pooled["g"] - pooled["ci_lo"], pooled["ci_hi"] - pooled["g"]],
        fmt="o", color="#333333", ecolor="#777777", capsize=2,
    )
    ax.set_yticks(y, [f"F{fid} {label}" for fid, label in zip(pooled["family_id"], pooled["short"])], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Pooled Hedges' g (family aCAIrel, AF − SR)")
    ax.set_title("Corrected pooled family effects", loc="left")
    save_all(fig, "fig2d_family_effects")


def consistency_table(effects: pd.DataFrame, states: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for family_id in FAMILIES["family_id"]:
        subset = effects[(effects["family_id"] == family_id) & (effects["cohort"].isin(COHORTS))]
        pooled = effects[(effects["family_id"] == family_id) & (effects["cohort"] == "POOLED")].iloc[0]
        pooled_sign = int(np.sign(pooled["g"]))
        same_sign = int((np.sign(subset["g"]) == pooled_sign).sum())
        cohort_sig = int(((subset["ci_lo"] > 0) | (subset["ci_hi"] < 0)).sum())
        pooled_sig = bool((pooled["ci_lo"] > 0) or (pooled["ci_hi"] < 0))
        replication = "high" if pooled_sig and same_sign >= 2 else ("moderate" if same_sign >= 2 else "low")
        weights = 1 / subset["var_g"].to_numpy()
        q = float(np.sum(weights * (subset["g"].to_numpy() - pooled["g"]) ** 2))
        df = len(subset) - 1
        i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
        p_het = float(1 - chi2.cdf(q, df))
        state_map = (
            states[states["family_id"] == family_id]
            .set_index("cohort")["state"]
            .reindex(COHORTS)
        )
        counts = state_map.value_counts()
        dominant = str(counts.index[0])
        concordance_n = int(counts.iloc[0])
        record = {
            "family_id": int(family_id),
            "short": pooled["short"],
            "direction": pooled["direction"],
            "pooled_g": pooled["g"],
            "pooled_ci_lo": pooled["ci_lo"],
            "pooled_ci_hi": pooled["ci_hi"],
            "pooled_standard_error": (pooled["ci_hi"] - pooled["ci_lo"]) / (2 * 1.96),
            "pooled_sig": pooled_sig,
            "same_sign_n": same_sign,
            "cohort_sig_n": cohort_sig,
            "replication_class": replication,
            "pooling_model": "fixed-effect_inverse_variance",
            "tau2": "NA",
            "Q": q,
            "I2": i2,
            "p_het": p_het,
            "state_128188": state_map["GSE128188"],
            "state_79768": state_map["GSE79768"],
            "state_115574": state_map["GSE115574"],
            "dominant_state": dominant,
            "state_concordance_ratio": concordance_n / 3,
            "state_concordance_n": concordance_n,
            "state_concordance_class": "stable" if concordance_n >= 2 else "heterogeneous",
        }
        for cohort in COHORTS:
            suffix = cohort.replace("GSE", "")
            cohort_row = subset[subset["cohort"] == cohort].iloc[0]
            record[f"g_{suffix}"] = cohort_row["g"]
            record[f"ci_lo_{suffix}"] = cohort_row["ci_lo"]
            record[f"ci_hi_{suffix}"] = cohort_row["ci_hi"]
        rows.append(record)
    consistency = pd.DataFrame(rows).merge(
        FAMILIES[["family_id", "canonical_family_display_order"]],
        on="family_id",
        validate="one_to_one",
    ).sort_values("canonical_family_display_order")
    consistency.to_csv(
        RESULTS / "fig2e_cross_cohort_consistency.tsv", sep="\t", index=False
    )
    consistency[
        [
            "family_id", "direction", "pooled_g", "pooled_ci_lo", "pooled_ci_hi",
            "I2", "replication_class", "dominant_state", "state_concordance_n",
        ]
    ].to_csv(RESULTS / "fig2_replication_classes.tsv", sep="\t", index=False)
    return consistency


def plot_consistency(consistency: pd.DataFrame) -> None:
    codes = consistency[["state_128188", "state_79768", "state_115574"]].apply(
        lambda column: column.map({state: i for i, state in enumerate(STATE_ORDER)})
    ).to_numpy(dtype=float)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9, 7), gridspec_kw={"width_ratios": [2, 1]})
    ax.imshow(codes, aspect="auto", cmap=ListedColormap(STATE_COLORS), vmin=-0.5, vmax=3.5)
    ax.set_xticks(range(3), [c.replace("GSE", "") for c in COHORTS])
    ax.set_yticks(range(16), [f"F{x}" for x in consistency["family_id"]])
    ax.set_title("Cross-cohort states", loc="left")
    colors = consistency["replication_class"].map({"high": "#2E7B5E", "moderate": "#E08214", "low": "#999999"})
    ax2.barh(np.arange(16), consistency["same_sign_n"], color=colors)
    ax2.set_xlim(0, 3)
    ax2.set_yticks([])
    ax2.set_xlabel("Cohorts matching pooled sign")
    ax2.set_title("Replication class", loc="left")
    save_all(fig, "fig2e_cross_cohort_consistency")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    donor_heatmap()
    states = states_table()
    plot_states(states)
    effects = effects_table()
    plot_effects(effects)
    consistency = consistency_table(effects, states)
    plot_consistency(consistency)
    print("Regenerated corrected Figure 2B-E tables and panels")


if __name__ == "__main__":
    main()
