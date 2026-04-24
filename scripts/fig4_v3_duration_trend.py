#!/usr/bin/env python3
"""
Fig 4 v3 - duration trend on GSE41177 (LAA + LAPV) for the NCP_v3 panel.

Pipeline:
    GSE41177 expr  ->  rank singscore on Fig3-derived LE genes per program
                  ->  family-level aggregation (mean of program-level scores)
                  ->  Spearman rho vs numeric AF duration (months)
                       (AF-only, n=12 LAA donors with numeric duration; same for LAPV)
                  ->  permutation p (10,000 sample-label shuffles)
                  ->  bootstrap 95% CI (5,000 resamples)
                  ->  3 sensitivity codings on LAA:
                         primary    : numeric months (drop ">1M" rows)
                         gt1m_eq1   : ">1M" -> 1 month
                         ordinal    : 3-tier (SR=0, short<24, long>=24)

Inputs
------
    Project_P4/GSE41177/expr_symbol_LAA.tsv
    Project_P4/GSE41177/expr_symbol_LAPV.tsv
    Project_P4/GSE41177/sample_metadata.tsv
    results/fig3_v3_program_le.tsv      (LE genes per program; same as Fig 3)
    results/innovation_ncp_v3_axis_families.tsv

Outputs
-------
    results/fig4_v3_singscore_per_sample.tsv
    results/fig4_v3_PROGRAM_laa_primary.tsv
    results/fig4_v3_FAMILY_laa_primary.tsv
    results/fig4_v3_FAMILY_laa_sens_gt1m_eq1.tsv
    results/fig4_v3_FAMILY_laa_sens_ordinal.tsv
    results/fig4_v3_FAMILY_lapv_sensitivity.tsv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent

EXPR_LAA = ROOT / "Project_P4" / "GSE41177" / "expr_symbol_LAA.tsv"
EXPR_LAPV = ROOT / "Project_P4" / "GSE41177" / "expr_symbol_LAPV.tsv"
META_TSV = ROOT / "Project_P4" / "GSE41177" / "sample_metadata.tsv"
LE_TSV = ROOT / "results" / "fig3_v3_program_le.tsv"
FAM_TSV = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"

OUT_LONG = ROOT / "results" / "fig4_v3_singscore_per_sample.tsv"
OUT_PROG_PRIMARY = ROOT / "results" / "fig4_v3_PROGRAM_laa_primary.tsv"
OUT_FAM_PRIMARY = ROOT / "results" / "fig4_v3_FAMILY_laa_primary.tsv"
OUT_FAM_GT1M = ROOT / "results" / "fig4_v3_FAMILY_laa_sens_gt1m_eq1.tsv"
OUT_FAM_ORD = ROOT / "results" / "fig4_v3_FAMILY_laa_sens_ordinal.tsv"
OUT_FAM_LAPV = ROOT / "results" / "fig4_v3_FAMILY_lapv_sensitivity.tsv"

N_PERM = 10_000
N_BOOT = 5_000
SEED = 20_260_418


def _rank_avg(a: np.ndarray) -> np.ndarray:
    """Average-rank within a 1D array.  Faster than scipy.stats.rankdata for
    small n by avoiding kwarg overhead."""
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty_like(a, dtype=float)
    n = len(a)
    i = 0
    while i < n:
        j = i + 1
        while j < n and a[order[j]] == a[order[i]]:
            j += 1
        # average rank for tied block (1-based)
        avg = 0.5 * (i + 1 + j)
        for k in range(i, j):
            ranks[order[k]] = avg
        i = j
    return ranks


def _spearman_fast(x_rank: np.ndarray, y_rank: np.ndarray) -> float:
    """Spearman = Pearson on ranks; no scipy call, ~ 30x faster for n=12."""
    n = len(x_rank)
    if n < 3:
        return float("nan")
    xc = x_rank - x_rank.mean()
    yc = y_rank - y_rank.mean()
    denom = np.sqrt((xc * xc).sum() * (yc * yc).sum())
    if denom == 0:
        return float("nan")
    return float((xc * yc).sum() / denom)


def per_sample_rank(expr):
    """Within each sample (column), rank gene values (avg ties); norm in (0,1]."""
    ranks = expr.rank(axis=0, method="average")
    return ranks / ranks.shape[0]


def singscore(norm_ranks: pd.DataFrame, gene_set: list[str]) -> pd.Series:
    gs = [g for g in gene_set if g in norm_ranks.index]
    if len(gs) < 5:
        return pd.Series(dtype=float)
    return norm_ranks.loc[gs].mean(axis=0) - 0.5


def perm_spearman(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """Permutation Spearman: pre-compute ranks once, then permute one rank
    vector and recompute Pearson on ranks (= Spearman)."""
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
    """Bootstrap 95% CI for Spearman; rank inside each resample (ties may
    differ across resamples so cannot pre-rank)."""
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


def compute_trend(scores_long: pd.DataFrame, level: str, duration_col: str,
                   meta: pd.DataFrame, tag: str, rng: np.random.Generator) -> pd.DataFrame:
    """level in {'program', 'family'}; duration_col is the metadata column to
    use as x (numeric).  tag is a string written to the output table."""
    # Average to family level if needed.
    score_col = "score"
    keys = ["family_id", "direction", "family_source"] if level == "family" else \
            ["program_id", "program_order", "direction", "family_id", "family_source"]
    long = scores_long.copy()
    if level == "family":
        long = (
            long.groupby(["family_id", "direction", "family_source", "gsm"], as_index=False)
            ["score"].mean()
        )
    # Attach metadata duration column
    long = long.merge(meta[["gsm", "donor_id", "rhythm", "tissue", duration_col, "keep_for_primary"]], on="gsm", how="left")

    rows = []
    if level == "family":
        group_iter = long.groupby(["family_id", "direction", "family_source"])
    else:
        group_iter = long.groupby(["program_id", "program_order", "direction", "family_id", "family_source"])

    for keys_vals, sub in group_iter:
        # AF-only with numeric duration
        af = sub[(sub["rhythm"] == "AF") & sub[duration_col].notna()]
        # For sensitivity codings the input column is already numeric.
        af = af[pd.to_numeric(af[duration_col], errors="coerce").notna()].copy()
        if len(af) < 4:
            continue
        x = af[duration_col].astype(float).to_numpy()
        y = af[score_col].astype(float).to_numpy()

        rho, p = perm_spearman(x, y, rng)
        ci_lo, ci_hi = bootstrap_ci(x, y, rng)

        rec = {"analysis": tag, "n_samples": int(len(af)), "rho_spearman": rho,
                "p_perm": p, "bootstrap_ci_lo": ci_lo, "bootstrap_ci_hi": ci_hi}
        if level == "family":
            rec.update({"family_id": int(keys_vals[0]),
                        "direction": keys_vals[1],
                        "family_source": keys_vals[2]})
        else:
            rec.update({"program_id": keys_vals[0],
                        "program_order": int(keys_vals[1]),
                        "direction": keys_vals[2],
                        "family_id": int(keys_vals[3]),
                        "family_source": keys_vals[4]})
        rows.append(rec)
    return pd.DataFrame(rows)


def score_one_tissue(expr_path: Path, meta_tissue: pd.DataFrame, le_lookup: dict[str, list[str]],
                       panel: pd.DataFrame, tissue_label: str) -> pd.DataFrame:
    print(f"  -> reading {expr_path.name}")
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    print(f"     expr shape: {expr.shape}")
    samples_keep = [s for s in meta_tissue["gsm"].tolist() if s in expr.columns]
    expr = expr[samples_keep]
    norm_ranks = per_sample_rank(expr)
    rows = []
    for _, r in panel.iterrows():
        pw = str(r["pathway"])
        le = le_lookup.get(pw, [])
        s = singscore(norm_ranks, le)
        if s.empty:
            continue
        for gsm, score in s.items():
            rows.append({
                "gsm": gsm,
                "tissue": tissue_label,
                "program_id": pw,
                "program_order": int(r["program_order"]),
                "direction": r["direction"],
                "family_id": int(r["family_id"]) if pd.notna(r["family_id"]) else -1,
                "family_source": r["family_source"] if pd.notna(r["family_source"]) else "unknown",
                "score": float(score),
            })
    return pd.DataFrame(rows)


def main():
    rng = np.random.default_rng(SEED)
    panel = pd.read_csv(FAM_TSV, sep="\t")
    le_df = pd.read_csv(LE_TSV, sep="\t")
    le_lookup = {
        str(r["program_id"]): [g for g in str(r["leading_edge_genes"]).split(",") if g]
        for _, r in le_df.iterrows()
    }
    meta = pd.read_csv(META_TSV, sep="\t")
    meta["duration_months"] = pd.to_numeric(meta["duration_months"], errors="coerce")

    # Sensitivity coding: ">1M" rows (duration_months NaN, rhythm AF) -> coded 1
    meta_gt1m = meta.copy()
    meta_gt1m["duration_gt1m_eq1"] = meta_gt1m["duration_months"]
    fillmask = (meta_gt1m["duration_months"].isna()) & (meta_gt1m["rhythm"] == "AF") & (meta_gt1m["duration_raw"].astype(str).str.contains(">1M", regex=False))
    meta_gt1m.loc[fillmask, "duration_gt1m_eq1"] = 1.0

    # Sensitivity coding: ordinal 3-tier  (SR=0, short_AF=1, long_AF=2)
    meta_ord = meta.copy()
    tier_map = {"SR": 0, "short_AF": 1, "long_AF": 2}
    meta_ord["duration_ordinal"] = meta_ord["duration_tier_B"].map(tier_map)

    print("Scoring LAA samples...")
    laa_meta = meta[meta["tissue"] == "LAA"].copy()
    scores_laa = score_one_tissue(EXPR_LAA, laa_meta, le_lookup, panel, "LAA")
    print(f"     LAA scores: {len(scores_laa)} (sample x program)")

    print("Scoring LAPV samples...")
    lapv_meta = meta[meta["tissue"] == "LAPV"].copy()
    scores_lapv = score_one_tissue(EXPR_LAPV, lapv_meta, le_lookup, panel, "LAPV")
    print(f"     LAPV scores: {len(scores_lapv)} (sample x program)")

    long_all = pd.concat([scores_laa, scores_lapv], axis=0, ignore_index=True)
    OUT_LONG.parent.mkdir(parents=True, exist_ok=True)
    long_all.to_csv(OUT_LONG, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_LONG.relative_to(ROOT)} ({len(long_all)} rows)")

    # ---- LAA primary (numeric months) ----
    laa_scores = scores_laa
    df_prog = compute_trend(laa_scores, "program", "duration_months", meta, "laa_primary", rng)
    df_prog.to_csv(OUT_PROG_PRIMARY, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_PROG_PRIMARY.relative_to(ROOT)} ({len(df_prog)} rows)")

    df_fam = compute_trend(laa_scores, "family", "duration_months", meta, "laa_primary", rng)
    df_fam.to_csv(OUT_FAM_PRIMARY, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_FAM_PRIMARY.relative_to(ROOT)} ({len(df_fam)} rows)")

    # ---- LAA sensitivity: gt1m_eq1 ----
    df_fam_gt1m = compute_trend(laa_scores, "family", "duration_gt1m_eq1", meta_gt1m, "laa_gt1m_eq1", rng)
    df_fam_gt1m.to_csv(OUT_FAM_GT1M, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_FAM_GT1M.relative_to(ROOT)} ({len(df_fam_gt1m)} rows)")

    # ---- LAA sensitivity: ordinal ----
    df_fam_ord = compute_trend(laa_scores, "family", "duration_ordinal", meta_ord, "laa_ordinal", rng)
    df_fam_ord.to_csv(OUT_FAM_ORD, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_FAM_ORD.relative_to(ROOT)} ({len(df_fam_ord)} rows)")

    # ---- LAPV sensitivity (region) ----
    df_fam_lapv = compute_trend(scores_lapv, "family", "duration_months", meta, "lapv_sensitivity", rng)
    df_fam_lapv.to_csv(OUT_FAM_LAPV, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_FAM_LAPV.relative_to(ROOT)} ({len(df_fam_lapv)} rows)")


if __name__ == "__main__":
    main()
