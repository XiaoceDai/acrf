"""Locked estimands for Phase 11B1 BMC acrf benchmark ablations."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

COHORT_ORDER = ["GSE128188", "GSE79768", "GSE115574"]
FAMILY_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]
DISEASE_RHYTHMS = ("AF", "SR")


def hedges_g(af: np.ndarray, sr: np.ndarray) -> tuple[float, float]:
    """Small-sample-corrected Hedges' g (AF−SR) and sampling variance.

    Matches production fig2d_family_effects_consistency.hedges_g.
    """
    nx, ny = len(af), len(sr)
    sx, sy = np.std(af, ddof=1), np.std(sr, ddof=1)
    pooled_sd = np.sqrt(((nx - 1) * sx**2 + (ny - 1) * sy**2) / (nx + ny - 2))
    correction = 1 - 3 / (4 * (nx + ny - 2) - 1)
    g = correction * (np.mean(af) - np.mean(sr)) / pooled_sd
    variance = (nx + ny) / (nx * ny) + g**2 / (2 * (nx + ny))
    return float(g), float(variance)


def welch_mean_diff_var(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Welch unpaired mean difference, variance, and Welch–Satterthwaite df.

    Formula_id: WELCH_UNPAIRED_MEAN_DIFF_V1 (Phase 11B0 freeze).
    """
    nx, ny = len(x), len(y)
    vx = float(np.var(x, ddof=1)) / nx
    vy = float(np.var(y, ddof=1)) / ny
    diff = float(np.mean(x) - np.mean(y))
    var = vx + vy
    if vx == 0.0 and vy == 0.0:
        df = float(nx + ny - 2)
    else:
        df = (vx + vy) ** 2 / (vx**2 / (nx - 1) + vy**2 / (ny - 1))
    return diff, var, float(df)


def fe_ivw_pool(cohort_rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Fixed-effect inverse-variance pool with Higgins–Thompson I²."""
    valid = [
        r
        for r in cohort_rows
        if np.isfinite(r["variance"]) and r["variance"] > 0 and r.get("valid", True)
    ]
    if not valid:
        return {
            "effect": np.nan,
            "variance": np.nan,
            "standard_error": np.nan,
            "lower_CI": np.nan,
            "upper_CI": np.nan,
            "I2": np.nan,
            "n_AF": 0,
            "n_SR": 0,
            "valid": False,
            "invalid_reason": "no_valid_cohort_effects",
        }
    weights = np.array([1.0 / r["variance"] for r in valid], dtype=float)
    effects = np.array([r["effect"] for r in valid], dtype=float)
    pooled = float(np.sum(weights * effects) / weights.sum())
    pooled_var = float(1.0 / weights.sum())
    se = float(np.sqrt(pooled_var))
    q = float(np.sum(weights * (effects - pooled) ** 2))
    df = len(valid) - 1
    i2 = float(max(0.0, (q - df) / q) if q > 0 else 0.0)
    return {
        "effect": pooled,
        "variance": pooled_var,
        "standard_error": se,
        "lower_CI": pooled - 1.96 * se,
        "upper_CI": pooled + 1.96 * se,
        "I2": i2,
        "n_AF": int(sum(r["n_AF"] for r in valid)),
        "n_SR": int(sum(r["n_SR"] for r in valid)),
        "valid": True,
        "invalid_reason": "",
        "pooling_model": "fixed-effect_inverse_variance",
        "Q": q,
    }


def consistency_class_from_effects(
    cohort_effects: list[float],
    pooled_ci_lower: float,
    pooled_ci_upper: float,
    pooled_effect: float,
) -> tuple[str, int, int, bool]:
    """Apply locked High/Moderate/Low rule (identical to fig2d)."""
    pooled_sign = int(np.sign(pooled_effect))
    pooled_supported = bool((pooled_ci_lower > 0) or (pooled_ci_upper < 0))
    n_sign_match = int(sum(int(np.sign(e)) == pooled_sign for e in cohort_effects))
    if pooled_supported and n_sign_match >= 2:
        klass = "High"
    elif n_sign_match >= 2:
        klass = "Moderate"
    else:
        klass = "Low"
    return klass, n_sign_match, pooled_sign, pooled_supported


def effect_sign(value: float) -> int:
    """Sign convention: sign(0)=0."""
    if not np.isfinite(value):
        return 0
    return int(np.sign(value))


def program_first_family_mean(
    program_values: pd.Series, family_id: int, program_to_family: pd.Series
) -> float:
    """Equal-weight mean of constituent program values (program-first aggregation)."""
    members = program_to_family.index[program_to_family == family_id]
    vals = program_values.reindex(members).dropna()
    if vals.empty:
        return float("nan")
    return float(vals.mean())


def load_disease_program_table(path: str | Any) -> pd.DataFrame:
    """Load locked disease program contrast/score cache (paired bulk only)."""
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["rhythm"].isin(DISEASE_RHYTHMS)].copy()
    return df


def build_donor_family_scores(
    disease: pd.DataFrame,
    score_column: str,
) -> pd.DataFrame:
    """Aggregate program-level donor scores to family means (equal weight)."""
    rows: list[dict[str, Any]] = []
    for (dataset, donor_id, rhythm, family_id), grp in disease.groupby(
        ["dataset", "donor_id", "rhythm", "family_id"], sort=False
    ):
        rows.append(
            {
                "dataset": dataset,
                "donor_id": donor_id,
                "rhythm": rhythm,
                "family_id": int(family_id),
                "family_score": float(grp[score_column].mean()),
                "n_programs_observed": int(grp[score_column].notna().sum()),
            }
        )
    return pd.DataFrame(rows)


def hedges_effects_from_donor_family(
    donor_family: pd.DataFrame,
    family_meta: pd.DataFrame,
    *,
    method_id: str,
    stratum: str,
    estimand: str,
) -> pd.DataFrame:
    """Cohort Hedges' g + FE-IVW pool + consistency class from donor-level family scores."""
    label_map = family_meta.set_index("family_id")["short_label"].to_dict()
    rows: list[dict[str, Any]] = []
    for family_id in FAMILY_ORDER:
        cohort_records: list[dict[str, Any]] = []
        for cohort in COHORT_ORDER:
            sub = donor_family.loc[
                (donor_family["dataset"] == cohort)
                & (donor_family["family_id"] == family_id)
            ]
            af = sub.loc[sub["rhythm"] == "AF", "family_score"].to_numpy(dtype=float)
            sr = sub.loc[sub["rhythm"] == "SR", "family_score"].to_numpy(dtype=float)
            valid = len(af) >= 2 and len(sr) >= 2
            if valid:
                g, var = hedges_g(af, sr)
                se = float(np.sqrt(var))
                rec = {
                    "method_id": method_id,
                    "stratum": stratum,
                    "cohort": cohort,
                    "family_id": int(family_id),
                    "family_label": str(label_map.get(family_id, "")),
                    "estimand": estimand,
                    "n_AF": int(len(af)),
                    "n_SR": int(len(sr)),
                    "effect": float(g),
                    "standard_error": se,
                    "variance": float(var),
                    "lower_CI": float(g - 1.96 * se),
                    "upper_CI": float(g + 1.96 * se),
                    "effect_sign": effect_sign(g),
                    "valid": True,
                    "invalid_reason": "",
                    "estimate_type": "cohort-specific",
                    "I2": np.nan,
                    "pooling_model": "NA",
                }
            else:
                rec = {
                    "method_id": method_id,
                    "stratum": stratum,
                    "cohort": cohort,
                    "family_id": int(family_id),
                    "family_label": str(label_map.get(family_id, "")),
                    "estimand": estimand,
                    "n_AF": int(len(af)),
                    "n_SR": int(len(sr)),
                    "effect": np.nan,
                    "standard_error": np.nan,
                    "variance": np.nan,
                    "lower_CI": np.nan,
                    "upper_CI": np.nan,
                    "effect_sign": 0,
                    "valid": False,
                    "invalid_reason": "insufficient_donors",
                    "estimate_type": "cohort-specific",
                    "I2": np.nan,
                    "pooling_model": "NA",
                }
            rows.append(rec)
            cohort_records.append(rec)

        pooled = fe_ivw_pool(cohort_records)
        cohort_effects = [r["effect"] for r in cohort_records if r["valid"]]
        if pooled["valid"] and cohort_effects:
            klass, n_match, psign, psup = consistency_class_from_effects(
                cohort_effects,
                pooled["lower_CI"],
                pooled["upper_CI"],
                pooled["effect"],
            )
        else:
            klass, n_match, psign, psup = "Low", 0, 0, False
        pooled_row = {
            "method_id": method_id,
            "stratum": stratum,
            "cohort": "POOLED",
            "family_id": int(family_id),
            "family_label": str(label_map.get(family_id, "")),
            "estimand": estimand,
            "n_AF": pooled["n_AF"],
            "n_SR": pooled["n_SR"],
            "effect": pooled["effect"],
            "standard_error": pooled["standard_error"],
            "variance": pooled["variance"],
            "lower_CI": pooled["lower_CI"],
            "upper_CI": pooled["upper_CI"],
            "effect_sign": effect_sign(pooled["effect"]) if pooled["valid"] else 0,
            "valid": pooled["valid"],
            "invalid_reason": pooled.get("invalid_reason", ""),
            "estimate_type": "pooled",
            "I2": pooled.get("I2", np.nan),
            "pooling_model": pooled.get("pooling_model", "NA"),
            "consistency_class": klass,
            "n_sign_match": n_match,
            "pooled_sign": psign,
            "pooled_supported": psup,
            "Q": pooled.get("Q", np.nan),
        }
        rows.append(pooled_row)
        for r in rows:
            if r["family_id"] == family_id and r["estimate_type"] == "cohort-specific":
                r["consistency_class"] = klass
                r["n_sign_match"] = n_match
                r["pooled_sign"] = psign
                r["pooled_supported"] = psup
    return pd.DataFrame(rows)


def compute_a3_scores(disease: pd.DataFrame) -> pd.DataFrame:
    """A3: paired CAI_rel without direction alignment (CAI/|M|)."""
    out = disease.copy()
    out["method_score"] = out["delta"] / out["M"].abs()
    return out


def compute_a4_scores(disease: pd.DataFrame) -> pd.DataFrame:
    """A4: paired aCAI without healthy-magnitude scaling."""
    out = disease.copy()
    out["method_score"] = out["direction_sign"] * out["delta"]
    return out


def compute_a5_scores(disease: pd.DataFrame) -> pd.DataFrame:
    """A5: within-donor mean(LA, RA); donor remains biological unit."""
    out = disease.copy()
    out["method_score"] = (out["score_LA"] + out["score_RA"]) / 2.0
    return out


def compute_b1_scores(disease: pd.DataFrame, chamber: str) -> pd.DataFrame:
    """B1: plain singscore within one chamber stratum (LA or RA)."""
    out = disease.copy()
    col = "score_LA" if chamber == "LA" else "score_RA"
    out["method_score"] = out[col]
    return out


def compute_a2_unpaired_effects(
    disease: pd.DataFrame,
    family_meta: pd.DataFrame,
) -> pd.DataFrame:
    """A2 Welch unpaired chamber-mean contrast; AF−SR of aligned/scaled family contrasts.

    Estimand: within rhythm r, C_r = mean(score_LA) − mean(score_RA) treating chambers
    as independent samples of size n_r; aC_rel = s_p * C_r / |M|; program-first family
    mean; effect = aC_rel,AF − aC_rel,SR. Variance: Welch at program level, equal-weight
    propagation across programs, independent AF/SR addition.
    """
    label_map = family_meta.set_index("family_id")["short_label"].to_dict()
    rows: list[dict[str, Any]] = []
    sample_notes: list[dict[str, Any]] = []

    for family_id in FAMILY_ORDER:
        cohort_records: list[dict[str, Any]] = []
        for cohort in COHORT_ORDER:
            sub = disease.loc[disease["dataset"] == cohort]
            prog_effects_af: list[float] = []
            prog_vars_af: list[float] = []
            prog_effects_sr: list[float] = []
            prog_vars_sr: list[float] = []
            n_af = n_sr = 0
            dfs: list[float] = []

            for program_id, pgrp in sub.groupby("program_id", sort=False):
                if int(pgrp["family_id"].iloc[0]) != family_id:
                    continue
                s_p = float(pgrp["direction_sign"].iloc[0])
                m_abs = float(np.abs(pgrp["M"].iloc[0]))
                scale = s_p / m_abs
                af = pgrp.loc[pgrp["rhythm"] == "AF"]
                sr = pgrp.loc[pgrp["rhythm"] == "SR"]
                n_af = len(af)
                n_sr = len(sr)
                # Unpaired chamber means within rhythm
                c_af, v_af, df_af = welch_mean_diff_var(
                    af["score_LA"].to_numpy(dtype=float),
                    af["score_RA"].to_numpy(dtype=float),
                )
                c_sr, v_sr, df_sr = welch_mean_diff_var(
                    sr["score_LA"].to_numpy(dtype=float),
                    sr["score_RA"].to_numpy(dtype=float),
                )
                prog_effects_af.append(scale * c_af)
                prog_vars_af.append((scale**2) * v_af)
                prog_effects_sr.append(scale * c_sr)
                prog_vars_sr.append((scale**2) * v_sr)
                dfs.extend([df_af, df_sr])

            k = len(prog_effects_af)
            if k == 0 or n_af < 2 or n_sr < 2:
                rec = {
                    "method_id": "A2_NO_DONOR_PAIRING",
                    "stratum": "unstratified",
                    "cohort": cohort,
                    "family_id": int(family_id),
                    "family_label": str(label_map.get(family_id, "")),
                    "estimand": "E_A2_unpaired_chamber_mean_contrast_AFSR",
                    "n_AF": int(n_af),
                    "n_SR": int(n_sr),
                    "effect": np.nan,
                    "standard_error": np.nan,
                    "variance": np.nan,
                    "lower_CI": np.nan,
                    "upper_CI": np.nan,
                    "effect_sign": 0,
                    "valid": False,
                    "invalid_reason": "insufficient_programs_or_donors",
                    "estimate_type": "cohort-specific",
                    "I2": np.nan,
                    "pooling_model": "NA",
                    "formula_id": "WELCH_UNPAIRED_MEAN_DIFF_V1",
                    "df_rule": "Welch-Satterthwaite per chamber contrast; family via program mean",
                }
            else:
                fam_af = float(np.mean(prog_effects_af))
                fam_sr = float(np.mean(prog_effects_sr))
                # Equal-weight mean variance under program independence
                var_af = float(np.sum(prog_vars_af) / (k**2))
                var_sr = float(np.sum(prog_vars_sr) / (k**2))
                effect = fam_af - fam_sr
                var = var_af + var_sr
                se = float(np.sqrt(var))
                rec = {
                    "method_id": "A2_NO_DONOR_PAIRING",
                    "stratum": "unstratified",
                    "cohort": cohort,
                    "family_id": int(family_id),
                    "family_label": str(label_map.get(family_id, "")),
                    "estimand": "E_A2_unpaired_chamber_mean_contrast_AFSR",
                    "n_AF": int(n_af),
                    "n_SR": int(n_sr),
                    "effect": effect,
                    "standard_error": se,
                    "variance": var,
                    "lower_CI": effect - 1.96 * se,
                    "upper_CI": effect + 1.96 * se,
                    "effect_sign": effect_sign(effect),
                    "valid": True,
                    "invalid_reason": "",
                    "estimate_type": "cohort-specific",
                    "I2": np.nan,
                    "pooling_model": "NA",
                    "formula_id": "WELCH_UNPAIRED_MEAN_DIFF_V1",
                    "df_rule": "Welch-Satterthwaite per chamber contrast; family via program mean",
                    "variance_rule": "Var(C)=s_LA^2/n+s_RA^2/n; scale by (s_p/|M|)^2; family mean Var=sum(Var_p)/K^2; AF-SR add",
                    "n_programs": k,
                    "min_welch_df": float(np.min(dfs)) if dfs else np.nan,
                }
            rows.append(rec)
            cohort_records.append(rec)
            sample_notes.append(
                {
                    "cohort": cohort,
                    "family_id": family_id,
                    "n_AF": n_af,
                    "n_SR": n_sr,
                    "n_LA_AF": n_af,
                    "n_RA_AF": n_af,
                    "n_LA_SR": n_sr,
                    "n_RA_SR": n_sr,
                }
            )

        pooled = fe_ivw_pool(cohort_records)
        cohort_effects = [r["effect"] for r in cohort_records if r["valid"]]
        if pooled["valid"] and cohort_effects:
            klass, n_match, psign, psup = consistency_class_from_effects(
                cohort_effects,
                pooled["lower_CI"],
                pooled["upper_CI"],
                pooled["effect"],
            )
        else:
            klass, n_match, psign, psup = "Low", 0, 0, False
        rows.append(
            {
                "method_id": "A2_NO_DONOR_PAIRING",
                "stratum": "unstratified",
                "cohort": "POOLED",
                "family_id": int(family_id),
                "family_label": str(label_map.get(family_id, "")),
                "estimand": "E_A2_unpaired_chamber_mean_contrast_AFSR",
                "n_AF": pooled["n_AF"],
                "n_SR": pooled["n_SR"],
                "effect": pooled["effect"],
                "standard_error": pooled["standard_error"],
                "variance": pooled["variance"],
                "lower_CI": pooled["lower_CI"],
                "upper_CI": pooled["upper_CI"],
                "effect_sign": effect_sign(pooled["effect"]) if pooled["valid"] else 0,
                "valid": pooled["valid"],
                "invalid_reason": pooled.get("invalid_reason", ""),
                "estimate_type": "pooled",
                "I2": pooled.get("I2", np.nan),
                "pooling_model": pooled.get("pooling_model", "NA"),
                "consistency_class": klass,
                "n_sign_match": n_match,
                "pooled_sign": psign,
                "pooled_supported": psup,
                "Q": pooled.get("Q", np.nan),
                "formula_id": "WELCH_UNPAIRED_MEAN_DIFF_V1",
                "df_rule": "Welch-Satterthwaite per chamber contrast; family via program mean",
            }
        )
        for r in rows:
            if r["family_id"] == family_id and r["estimate_type"] == "cohort-specific":
                r["consistency_class"] = klass
                r["n_sign_match"] = n_match
                r["pooled_sign"] = psign
                r["pooled_supported"] = psup

    effects = pd.DataFrame(rows)
    effects.attrs["a2_sample_counts"] = pd.DataFrame(sample_notes)
    return effects


def locked_a1_from_fig2d(path: str | Any, family_meta: pd.DataFrame) -> pd.DataFrame:
    """Load locked A1 family×cohort effects; authority_status=LOCKED_REFERENCE."""
    src = pd.read_csv(path, sep="\t")
    label_map = family_meta.set_index("family_id")["short_label"].to_dict()
    rows: list[dict[str, Any]] = []
    for _, r in src.iterrows():
        rows.append(
            {
                "method_id": "A1_FULL_ACRF",
                "stratum": "unstratified",
                "cohort": r["cohort"],
                "family_id": int(r["family_id"]),
                "family_label": str(label_map.get(int(r["family_id"]), r["family_label"])),
                "estimand": "E_A1_paired_aCAI_rel_Hedges_g",
                "n_AF": int(r["n_AF"]),
                "n_SR": int(r["n_SR"]),
                "effect": float(r["hedges_g"]),
                "standard_error": float(r["standard_error"]),
                "variance": float(r["sampling_variance"]),
                "lower_CI": float(r["ci_lower"]),
                "upper_CI": float(r["ci_upper"]),
                "effect_sign": effect_sign(float(r["hedges_g"])),
                "valid": True,
                "invalid_reason": "",
                "estimate_type": r["estimate_type"],
                "I2": r["I2"] if r["I2"] != "NA" and pd.notna(r["I2"]) else np.nan,
                "pooling_model": r["pooling_model"],
                "consistency_class": r["consistency_class"],
                "n_sign_match": r["n_sign_match"] if pd.notna(r["n_sign_match"]) else np.nan,
                "pooled_sign": r["pooled_sign"] if pd.notna(r["pooled_sign"]) else np.nan,
                "pooled_supported": r["pooled_supported"]
                if pd.notna(r["pooled_supported"])
                else np.nan,
                "authority_status": "LOCKED_REFERENCE",
            }
        )
    return pd.DataFrame(rows)


def validate_a1_replay(
    disease_family_acairel: pd.DataFrame,
    locked: pd.DataFrame,
    atol_effect: float = 1e-10,
    atol_var: float = 1e-10,
) -> dict[str, Any]:
    """Replay A1 Hedges' g from family_aCAIrel; require agreement with locked authority."""
    mismatches: list[dict[str, Any]] = []
    cohort_locked = locked.loc[locked["estimate_type"] == "cohort-specific"]
    for _, row in cohort_locked.iterrows():
        sub = disease_family_acairel.loc[
            (disease_family_acairel["dataset"] == row["cohort"])
            & (disease_family_acairel["family_id"] == row["family_id"])
        ]
        af = sub.loc[sub["rhythm"] == "AF", "family_aCAIrel"].to_numpy(dtype=float)
        sr = sub.loc[sub["rhythm"] == "SR", "family_aCAIrel"].to_numpy(dtype=float)
        g, var = hedges_g(af, sr)
        if abs(g - row["effect"]) > atol_effect or abs(var - row["variance"]) > atol_var:
            mismatches.append(
                {
                    "cohort": row["cohort"],
                    "family_id": int(row["family_id"]),
                    "locked_effect": row["effect"],
                    "replay_effect": g,
                    "locked_variance": row["variance"],
                    "replay_variance": var,
                }
            )
    return {
        "status": "PASS" if not mismatches else "FAIL",
        "n_compared": int(len(cohort_locked)),
        "n_mismatch": len(mismatches),
        "atol_effect": atol_effect,
        "atol_var": atol_var,
        "mismatches": mismatches,
    }
