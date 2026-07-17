#!/usr/bin/env python3
"""Corrected GSE41177 duration trend analysis with exact-reference scores."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "ncp_v3_singscore_reanalysis"
META_TSV = OUT_DIR / "gse41177_duration_metadata.tsv"
FAM_SCORES = OUT_DIR / "gse41177_family_sample_scores.tsv"
OUT_META = META_TSV
OUT_LAA = OUT_DIR / "gse41177_family_duration_laa_primary.tsv"
OUT_LAPV = OUT_DIR / "gse41177_family_duration_lapv_sensitivity.tsv"
OUT_CHRO = OUT_DIR / "fig4_v3_chronicity_class.tsv"

N_PERM = 10_000
N_BOOT = 5_000
SEED = 20_260_418


def _rank_avg(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty_like(values, dtype=float)
    n = len(values)
    i = 0
    while i < n:
        j = i + 1
        while j < n and values[order[j]] == values[order[i]]:
            j += 1
        avg = 0.5 * (i + 1 + j)
        for k in range(i, j):
            ranks[order[k]] = avg
        i = j
    return ranks


def _spearman_fast(x_rank: np.ndarray, y_rank: np.ndarray) -> float:
    n = len(x_rank)
    if n < 3:
        return float("nan")
    xc = x_rank - x_rank.mean()
    yc = y_rank - y_rank.mean()
    denom = np.sqrt((xc * xc).sum() * (yc * yc).sum())
    if denom == 0:
        return float("nan")
    return float((xc * yc).sum() / denom)


def perm_spearman(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    n = len(x)
    if n < 3:
        return float("nan"), float("nan")
    xr = _rank_avg(x)
    yr = _rank_avg(y)
    obs = _spearman_fast(xr, yr)
    if not np.isfinite(obs):
        return float("nan"), float("nan")
    hits = 1
    abs_obs = abs(obs)
    for _ in range(N_PERM):
        perm = rng.permutation(n)
        s = _spearman_fast(xr[perm], yr)
        if abs(s) >= abs_obs - 1e-12:
            hits += 1
    return obs, hits / (N_PERM + 1)


def bootstrap_ci(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    n = len(x)
    if n < 4:
        return float("nan"), float("nan")
    samples = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        xr = _rank_avg(x[idx])
        yr = _rank_avg(y[idx])
        s = _spearman_fast(xr, yr)
        if np.isfinite(s):
            samples.append(s)
    if not samples:
        return float("nan"), float("nan")
    return float(np.percentile(samples, 2.5)), float(np.percentile(samples, 97.5))


def bh_qvalues(p_values: pd.Series) -> pd.Series:
    p = p_values.astype(float)
    n = p.notna().sum()
    if n == 0:
        return pd.Series(np.nan, index=p.index)
    order = p.sort_values()
    ranks = np.arange(1, len(order) + 1, dtype=float)
    q = order * n / ranks
    q = q[::-1].cummin()[::-1]
    out = pd.Series(np.nan, index=p.index, dtype=float)
    out.loc[order.index] = q.values
    return out.clip(upper=1.0)


def classify(rho: float) -> str:
    if not np.isfinite(rho):
        return "Inconclusive"
    a = abs(rho)
    if a >= 0.4:
        return "Strong_trend"
    if a >= 0.2:
        return "Weak_trend"
    return "Stable"


def compute_region_trends(
    scores: pd.DataFrame,
    meta: pd.DataFrame,
    region: str,
    *,
    tag: str,
    rng: np.random.Generator,
) -> pd.DataFrame:
    sub_meta = meta[meta["region"] == region].copy()
    sub_scores = scores[scores["region"] == region].copy()
    merged = sub_scores.merge(
        sub_meta[
            ["sample_id", "donor_id", "rhythm", "duration_months", "primary_duration_eligible"]
        ],
        on="sample_id",
        how="left",
    )
    rows = []
    for (fid, direction, family_source), group in merged.groupby(
        ["family_id", "direction", "family_source"]
    ):
        af = group[
            (group["rhythm"] == "AF") & group["duration_months"].notna()
        ].copy()
        if len(af) < 4:
            continue
        x = af["duration_months"].astype(float).to_numpy()
        y = af["score"].astype(float).to_numpy()
        rho, p_perm = perm_spearman(x, y, rng)
        ci_lo, ci_hi = bootstrap_ci(x, y, rng)
        p_asymp = float(spearmanr(x, y).pvalue) if len(af) >= 3 else float("nan")
        rows.append(
            {
                "analysis": tag,
                "region": region,
                "family_id": int(fid),
                "direction": direction,
                "family_source": family_source,
                "n_samples": int(len(af)),
                "rho_spearman": rho,
                "p_perm": p_perm,
                "p_spearman": p_asymp,
                "bootstrap_ci_lo": ci_lo,
                "bootstrap_ci_hi": ci_hi,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_spearman_bh"] = bh_qvalues(out["p_spearman"])
    return out


def main() -> None:
    rng = np.random.default_rng(SEED)
    meta = pd.read_csv(META_TSV, sep="\t")
    scores = pd.read_csv(FAM_SCORES, sep="\t")

    laa = compute_region_trends(scores, meta, "LAA", tag="laa_primary", rng=rng)
    lapv = compute_region_trends(scores, meta, "LAPV", tag="lapv_sensitivity", rng=rng)
    assert len(laa) == 16, f"expected 16 LAA family trends, got {len(laa)}"
    assert len(lapv) == 16, f"expected 16 LAPV family trends, got {len(lapv)}"
    assert (laa["n_samples"] == 12).all(), "LAA primary must use 12 AF donors"

    laa.to_csv(OUT_LAA, sep="\t", index=False, float_format="%.6f")
    lapv.to_csv(OUT_LAPV, sep="\t", index=False, float_format="%.6f")

    laa_idx = laa.set_index("family_id")
    lapv_idx = lapv.set_index("family_id")
    rows = []
    for fid in sorted(laa_idx.index):
        prow = laa_idx.loc[fid]
        rho = float(prow["rho_spearman"])
        cls = classify(rho)
        rho_lapv = float(lapv_idx.loc[fid, "rho_spearman"]) if fid in lapv_idx.index else np.nan
        rho_sign = np.sign(rho) if rho != 0 else 0
        if np.isfinite(rho_lapv):
            lapv_sign = np.sign(rho_lapv)
            sens_same_sign_n = int(lapv_sign == rho_sign or rho_sign == 0)
            sens_total_n = 2
            regional_sign_agreement = bool(lapv_sign == rho_sign)
        else:
            sens_same_sign_n = 0
            sens_total_n = 1
            regional_sign_agreement = False

        rows.append(
            {
                "family_id": int(fid),
                "direction": str(prow["direction"]),
                "rho_primary": rho,
                "p_perm_primary": float(prow["p_perm"]),
                "ci_lo_primary": float(prow["bootstrap_ci_lo"]),
                "ci_hi_primary": float(prow["bootstrap_ci_hi"]),
                "n_samples_primary": int(prow["n_samples"]),
                "chronicity_class": cls,
                "sens_same_sign_n": sens_same_sign_n,
                "sens_total_n": sens_total_n,
                "rho_gt1m_eq1": np.nan,
                "rho_ordinal": np.nan,
                "rho_lapv": rho_lapv,
                "regional_sign_agreement": regional_sign_agreement,
                "p_spearman_primary": float(prow["p_spearman"]),
                "q_spearman_bh_primary": float(prow["q_spearman_bh"]),
            }
        )

    chro = pd.DataFrame(rows)
    cls_order = {"Strong_trend": 0, "Weak_trend": 1, "Stable": 2, "Inconclusive": 3}
    chro["cls_rank"] = chro["chronicity_class"].map(cls_order)
    chro = chro.sort_values(["cls_rank", "family_id"]).drop(columns=["cls_rank"])
    chro.to_csv(OUT_CHRO, sep="\t", index=False, float_format="%.4f")

    print(f"wrote {OUT_META.relative_to(ROOT)}")
    print(f"wrote {OUT_LAA.relative_to(ROOT)} ({len(laa)} rows)")
    print(f"wrote {OUT_LAPV.relative_to(ROOT)} ({len(lapv)} rows)")
    print(f"wrote {OUT_CHRO.relative_to(ROOT)} ({len(chro)} families)")
    print("\n=== corrected class counts ===")
    print(chro["chronicity_class"].value_counts().to_string())


if __name__ == "__main__":
    main()
