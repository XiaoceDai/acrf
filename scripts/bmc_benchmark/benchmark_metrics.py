"""Primary and secondary metrics for Phase 11B1 BMC acrf benchmark."""

from __future__ import annotations

from itertools import combinations
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from benchmark_estimands import COHORT_ORDER, FAMILY_ORDER, effect_sign

COHORT_PAIRS = list(combinations(COHORT_ORDER, 2))


def _method_arms(effects: pd.DataFrame) -> list[tuple[str, str]]:
    """Unique (method_id, stratum) arms with cohort-specific rows."""
    sub = effects.loc[effects["estimate_type"] == "cohort-specific"]
    return sorted(
        {(str(m), str(s)) for m, s in zip(sub["method_id"], sub["stratum"])}
    )


def _cohort_effect_vector(
    effects: pd.DataFrame, method_id: str, stratum: str, cohort: str
) -> pd.Series:
    """Family-indexed effect vector for one method arm and cohort."""
    sub = effects.loc[
        (effects["method_id"] == method_id)
        & (effects["stratum"] == stratum)
        & (effects["cohort"] == cohort)
        & (effects["estimate_type"] == "cohort-specific")
    ].set_index("family_id")
    return sub.reindex(FAMILY_ORDER)["effect"]


def metric_m1_correlations(
    effects: pd.DataFrame,
    *,
    source_input_hash: str,
    software_version: str,
    seed: str = "NA",
) -> pd.DataFrame:
    """M1: Spearman correlation of 16 family effects across cohort pairs."""
    rows: list[dict[str, Any]] = []
    for method_id, stratum in _method_arms(effects):
        pair_rhos: list[float] = []
        pair_ns: list[int] = []
        for c1, c2 in COHORT_PAIRS:
            v1 = _cohort_effect_vector(effects, method_id, stratum, c1)
            v2 = _cohort_effect_vector(effects, method_id, stratum, c2)
            mask = v1.notna() & v2.notna() & np.isfinite(v1) & np.isfinite(v2)
            n_complete = int(mask.sum())
            if n_complete >= 3:
                rho, _ = spearmanr(v1[mask], v2[mask])
                rho = float(rho)
            else:
                rho = float("nan")
            pair_rhos.append(rho)
            pair_ns.append(n_complete)
            rows.append(
                {
                    "metric_id": "M1_CROSS_COHORT_EFFECT_CORRELATION",
                    "method_id": method_id,
                    "stratum": stratum,
                    "cohort": f"{c1}|{c2}",
                    "family_id": "ALL",
                    "estimand": _estimand(effects, method_id, stratum),
                    "pairwise_rho": rho,
                    "n_complete_families": n_complete,
                    "median_rho": np.nan,
                    "rho_min": np.nan,
                    "rho_max": np.nan,
                    "source_input_hash": source_input_hash,
                    "software_version": software_version,
                    "seed": seed,
                    "valid": bool(np.isfinite(rho)),
                    "notes": "Spearman; no P-value as primary",
                }
            )
        finite = [r for r in pair_rhos if np.isfinite(r)]
        med = float(np.median(finite)) if finite else float("nan")
        rmin = float(np.min(finite)) if finite else float("nan")
        rmax = float(np.max(finite)) if finite else float("nan")
        rows.append(
            {
                "metric_id": "M1_CROSS_COHORT_EFFECT_CORRELATION",
                "method_id": method_id,
                "stratum": stratum,
                "cohort": "SUMMARY",
                "family_id": "ALL",
                "estimand": _estimand(effects, method_id, stratum),
                "pairwise_rho": np.nan,
                "n_complete_families": int(min(pair_ns)) if pair_ns else 0,
                "median_rho": med,
                "rho_min": rmin,
                "rho_max": rmax,
                "source_input_hash": source_input_hash,
                "software_version": software_version,
                "seed": seed,
                "valid": bool(finite),
                "notes": f"pair_rhos={pair_rhos}",
            }
        )
    return pd.DataFrame(rows)


def metric_m2_sign_concordance(
    effects: pd.DataFrame,
    *,
    source_input_hash: str,
    software_version: str,
    seed: str = "NA",
) -> pd.DataFrame:
    """M2: family×cohort-pair sign concordance; denominator 48 when complete."""
    rows: list[dict[str, Any]] = []
    for method_id, stratum in _method_arms(effects):
        denom_target = 16 * 3
        concordant = 0
        missing = 0
        family_hits = {fid: 0 for fid in FAMILY_ORDER}
        family_den = {fid: 0 for fid in FAMILY_ORDER}
        for fid in FAMILY_ORDER:
            signs = {}
            for cohort in COHORT_ORDER:
                v = _cohort_effect_vector(effects, method_id, stratum, cohort)
                val = v.loc[fid]
                signs[cohort] = (
                    None
                    if (pd.isna(val) or not np.isfinite(val))
                    else effect_sign(float(val))
                )
            for c1, c2 in COHORT_PAIRS:
                s1, s2 = signs[c1], signs[c2]
                if s1 is None or s2 is None:
                    missing += 1
                    continue
                family_den[fid] += 1
                # sign(0)=0; zero does not concord with nonzero (equality handles this)
                if s1 == s2:
                    concordant += 1
                    family_hits[fid] += 1
        report_denom = denom_target if missing == 0 else (denom_target - missing)
        prop = concordant / report_denom if report_denom else float("nan")
        rows.append(
            {
                "metric_id": "M2_CROSS_COHORT_SIGN_CONCORDANCE",
                "method_id": method_id,
                "stratum": stratum,
                "cohort": "ALL_PAIRS",
                "family_id": "ALL",
                "estimand": _estimand(effects, method_id, stratum),
                "concordant_count": concordant,
                "denominator": report_denom,
                "proportion": prop,
                "missing_comparisons": missing,
                "source_input_hash": source_input_hash,
                "software_version": software_version,
                "seed": seed,
                "valid": True,
                "notes": "sign(0)=0; zero nonconcordant with nonzero",
            }
        )
        for fid in FAMILY_ORDER:
            rows.append(
                {
                    "metric_id": "M2_CROSS_COHORT_SIGN_CONCORDANCE",
                    "method_id": method_id,
                    "stratum": stratum,
                    "cohort": "ALL_PAIRS",
                    "family_id": fid,
                    "estimand": _estimand(effects, method_id, stratum),
                    "concordant_count": family_hits[fid],
                    "denominator": family_den[fid],
                    "proportion": (family_hits[fid] / family_den[fid])
                    if family_den[fid]
                    else np.nan,
                    "missing_comparisons": 3 - family_den[fid],
                    "source_input_hash": source_input_hash,
                    "software_version": software_version,
                    "seed": seed,
                    "valid": True,
                    "notes": "family_level_breakdown",
                }
            )
    return pd.DataFrame(rows)


def metric_m3_rank_stability(
    effects: pd.DataFrame,
    *,
    source_input_hash: str,
    software_version: str,
    seed: str = "NA",
) -> pd.DataFrame:
    """M3: Spearman of family ranks across cohort pairs; median across pairs."""
    rows: list[dict[str, Any]] = []
    for method_id, stratum in _method_arms(effects):
        pair_rhos: list[float] = []
        for c1, c2 in COHORT_PAIRS:
            v1 = _cohort_effect_vector(effects, method_id, stratum, c1)
            v2 = _cohort_effect_vector(effects, method_id, stratum, c2)
            mask = v1.notna() & v2.notna() & np.isfinite(v1) & np.isfinite(v2)
            if mask.sum() >= 3:
                r1 = v1[mask].rank(method="average")
                r2 = v2[mask].rank(method="average")
                rho, _ = spearmanr(r1, r2)
                rho = float(rho)
            else:
                rho = float("nan")
            pair_rhos.append(rho)
            rows.append(
                {
                    "metric_id": "M3_FAMILY_RANK_STABILITY",
                    "method_id": method_id,
                    "stratum": stratum,
                    "cohort": f"{c1}|{c2}",
                    "family_id": "ALL",
                    "estimand": _estimand(effects, method_id, stratum),
                    "rank_spearman_rho": rho,
                    "median_rho": np.nan,
                    "source_input_hash": source_input_hash,
                    "software_version": software_version,
                    "seed": seed,
                    "valid": bool(np.isfinite(rho)),
                    "notes": "",
                }
            )
        finite = [r for r in pair_rhos if np.isfinite(r)]
        rows.append(
            {
                "metric_id": "M3_FAMILY_RANK_STABILITY",
                "method_id": method_id,
                "stratum": stratum,
                "cohort": "SUMMARY",
                "family_id": "ALL",
                "estimand": _estimand(effects, method_id, stratum),
                "rank_spearman_rho": np.nan,
                "median_rho": float(np.median(finite)) if finite else np.nan,
                "source_input_hash": source_input_hash,
                "software_version": software_version,
                "seed": seed,
                "valid": bool(finite),
                "notes": f"pair_rhos={pair_rhos}",
            }
        )
    return pd.DataFrame(rows)


def metric_m4_heterogeneity(
    effects: pd.DataFrame,
    *,
    source_input_hash: str,
    software_version: str,
    seed: str = "NA",
) -> pd.DataFrame:
    """M4: family-level Q and I² from FE/variance framework (descriptive; n=3)."""
    rows: list[dict[str, Any]] = []
    pooled = effects.loc[effects["estimate_type"] == "pooled"]
    for _, r in pooled.iterrows():
        method_id = str(r["method_id"])
        stratum = str(r["stratum"])
        fid = int(r["family_id"])
        cohort_rows = effects.loc[
            (effects["method_id"] == method_id)
            & (effects["stratum"] == stratum)
            & (effects["family_id"] == fid)
            & (effects["estimate_type"] == "cohort-specific")
            & (effects["valid"])
        ]
        if len(cohort_rows) < 2 or not bool(r["valid"]):
            q = np.nan
            i2 = np.nan
            valid = False
            reason = "insufficient_cohorts_or_invalid_pool"
        else:
            weights = 1.0 / cohort_rows["variance"].to_numpy(dtype=float)
            g = cohort_rows["effect"].to_numpy(dtype=float)
            pooled_g = float(r["effect"])
            q = float(np.sum(weights * (g - pooled_g) ** 2))
            df = len(cohort_rows) - 1
            i2 = float(max(0.0, (q - df) / q) if q > 0 else 0.0)
            valid = True
            reason = "n_cohorts=3_descriptive_only"
        rows.append(
            {
                "metric_id": "M4_EFFECT_HETEROGENEITY",
                "method_id": method_id,
                "stratum": stratum,
                "cohort": "POOLED",
                "family_id": fid,
                "estimand": r["estimand"],
                "Q": q,
                "I2": i2 if pd.isna(r.get("I2")) else r.get("I2"),
                "pooled_effect": r["effect"],
                "source_input_hash": source_input_hash,
                "software_version": software_version,
                "seed": seed,
                "valid": valid,
                "notes": reason,
            }
        )
    return pd.DataFrame(rows)


def metric_m5_class_stability(
    effects: pd.DataFrame,
    *,
    source_input_hash: str,
    software_version: str,
    seed: str = "NA",
) -> pd.DataFrame:
    """M5: High/Moderate/Low retention vs A1; B1 marked NOT_COMPARABLE."""
    a1 = effects.loc[
        (effects["method_id"] == "A1_FULL_ACRF")
        & (effects["estimate_type"] == "pooled")
    ].set_index("family_id")
    rows: list[dict[str, Any]] = []
    classes = ["High", "Moderate", "Low"]

    for method_id, stratum in _method_arms(effects):
        comparable = method_id != "B1_PLAIN_SINGSCORE"
        pooled = effects.loc[
            (effects["method_id"] == method_id)
            & (effects["stratum"] == stratum)
            & (effects["estimate_type"] == "pooled")
        ].set_index("family_id")
        retained = 0
        reclassified = 0
        transitions: dict[tuple[str, str], int] = {
            (a, b): 0 for a in classes for b in classes
        }
        for fid in FAMILY_ORDER:
            a1_class = str(a1.loc[fid, "consistency_class"])
            if not comparable:
                m_class = "NOT_COMPARABLE"
                notes = "B1 chamber-stratified class comparison not locked as comparable"
            else:
                m_class = str(pooled.loc[fid, "consistency_class"])
                notes = "A1 classes are not ground truth"
                if m_class == a1_class:
                    retained += 1
                else:
                    reclassified += 1
                if a1_class in classes and m_class in classes:
                    transitions[(a1_class, m_class)] += 1
            rows.append(
                {
                    "metric_id": "M5_CONSISTENCY_CLASS_STABILITY",
                    "method_id": method_id,
                    "stratum": stratum,
                    "cohort": "POOLED",
                    "family_id": fid,
                    "estimand": _estimand(effects, method_id, stratum),
                    "a1_class": a1_class,
                    "method_class": m_class,
                    "retained_vs_a1": (m_class == a1_class) if comparable else "NOT_COMPARABLE",
                    "source_input_hash": source_input_hash,
                    "software_version": software_version,
                    "seed": seed,
                    "valid": True,
                    "notes": notes,
                }
            )
        rows.append(
            {
                "metric_id": "M5_CONSISTENCY_CLASS_STABILITY",
                "method_id": method_id,
                "stratum": stratum,
                "cohort": "SUMMARY",
                "family_id": "ALL",
                "estimand": _estimand(effects, method_id, stratum),
                "a1_class": "NA",
                "method_class": "NA",
                "retained_vs_a1": retained if comparable else "NOT_COMPARABLE",
                "reclassified_count": reclassified if comparable else "NOT_COMPARABLE",
                "transition_matrix": str(transitions) if comparable else "NOT_COMPARABLE",
                "source_input_hash": source_input_hash,
                "software_version": software_version,
                "seed": seed,
                "valid": True,
                "notes": "retained_count relative to A1; not correctness",
            }
        )
    return pd.DataFrame(rows)


def metric_m6_pooled_sign(
    effects: pd.DataFrame,
    *,
    source_input_hash: str,
    software_version: str,
    seed: str = "NA",
) -> pd.DataFrame:
    """M6: pooled-effect sign agreement with A1 (descriptive)."""
    a1 = effects.loc[
        (effects["method_id"] == "A1_FULL_ACRF")
        & (effects["estimate_type"] == "pooled")
    ].set_index("family_id")
    rows: list[dict[str, Any]] = []
    for method_id, stratum in _method_arms(effects):
        if method_id == "A1_FULL_ACRF":
            continue
        pooled = effects.loc[
            (effects["method_id"] == method_id)
            & (effects["stratum"] == stratum)
            & (effects["estimate_type"] == "pooled")
        ].set_index("family_id")
        agree = 0
        for fid in FAMILY_ORDER:
            s_a1 = effect_sign(float(a1.loc[fid, "effect"]))
            s_m = effect_sign(float(pooled.loc[fid, "effect"]))
            hit = s_a1 == s_m
            agree += int(hit)
            rows.append(
                {
                    "metric_id": "M6_POOLED_EFFECT_DIRECTION_STABILITY",
                    "method_id": method_id,
                    "stratum": stratum,
                    "cohort": "POOLED",
                    "family_id": fid,
                    "estimand": _estimand(effects, method_id, stratum),
                    "a1_pooled_sign": s_a1,
                    "method_pooled_sign": s_m,
                    "sign_agreement": hit,
                    "source_input_hash": source_input_hash,
                    "software_version": software_version,
                    "seed": seed,
                    "valid": True,
                    "notes": "descriptive_only",
                }
            )
        rows.append(
            {
                "metric_id": "M6_POOLED_EFFECT_DIRECTION_STABILITY",
                "method_id": method_id,
                "stratum": stratum,
                "cohort": "SUMMARY",
                "family_id": "ALL",
                "estimand": _estimand(effects, method_id, stratum),
                "a1_pooled_sign": "NA",
                "method_pooled_sign": "NA",
                "sign_agreement": agree,
                "agreement_proportion": agree / 16,
                "source_input_hash": source_input_hash,
                "software_version": software_version,
                "seed": seed,
                "valid": True,
                "notes": "descriptive_only",
            }
        )
    return pd.DataFrame(rows)


def metric_m7_completeness(
    effects: pd.DataFrame,
    *,
    source_input_hash: str,
    software_version: str,
    seed: str = "NA",
) -> pd.DataFrame:
    """M7: estimand completeness / missingness by method arm."""
    rows: list[dict[str, Any]] = []
    cohort_eff = effects.loc[effects["estimate_type"] == "cohort-specific"]
    for method_id, stratum in _method_arms(effects):
        sub = cohort_eff.loc[
            (cohort_eff["method_id"] == method_id) & (cohort_eff["stratum"] == stratum)
        ]
        n_expected = 16 * 3
        n_valid = int(sub["valid"].sum())
        n_missing = n_expected - n_valid
        reasons = (
            sub.loc[~sub["valid"], "invalid_reason"].value_counts().to_dict()
            if n_missing
            else {}
        )
        rows.append(
            {
                "metric_id": "M7_ESTIMAND_COMPLETENESS",
                "method_id": method_id,
                "stratum": stratum,
                "cohort": "ALL",
                "family_id": "ALL",
                "estimand": _estimand(effects, method_id, stratum),
                "n_expected": n_expected,
                "n_valid": n_valid,
                "n_missing": n_missing,
                "missing_reasons": str(reasons),
                "source_input_hash": source_input_hash,
                "software_version": software_version,
                "seed": seed,
                "valid": True,
                "notes": "",
            }
        )
    return pd.DataFrame(rows)


def _estimand(effects: pd.DataFrame, method_id: str, stratum: str) -> str:
    sub = effects.loc[
        (effects["method_id"] == method_id) & (effects["stratum"] == stratum)
    ]
    if sub.empty:
        return "NA"
    return str(sub["estimand"].iloc[0])


def scan_unsupported_claims(texts: list[str]) -> list[str]:
    """Flag unsupported positive claim phrases (ignore explicit negations)."""
    banned = [
        "true biology",
        "ground truth",
        "universal superiority",
        "higher accuracy",
        "clinical prediction",
        "causal mechanism",
        "therapeutic target",
        "biomarker validation",
        "patient-level prediction",
    ]
    negation_prefixes = (
        "not ",
        "no ",
        "without ",
        "aren't ",
        "are not ",
        "is not ",
        "isn't ",
        "never ",
        "nor ",
    )
    hits: list[str] = []
    for text in texts:
        low = text.lower()
        for phrase in banned:
            start = 0
            while True:
                idx = low.find(phrase, start)
                if idx < 0:
                    break
                window = low[max(0, idx - 12) : idx]
                negated = any(window.rstrip().endswith(p.rstrip()) for p in negation_prefixes)
                # also catch "not ground truth" / "are not ground truth"
                if "not " + phrase in low or "no " + phrase in low:
                    negated = True
                if not negated:
                    hits.append(f"{phrase} :: {text[:120]}")
                start = idx + len(phrase)
    return hits
