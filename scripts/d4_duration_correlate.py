#!/usr/bin/env python3
"""
D4 - AF-duration correlation with sensitivity sidecars.

Writes four TSVs to Project_P4/GSE41177/results/ :

  1. duration_trend_laa_primary.tsv
       LAA, AF samples with numeric duration (n=12 expected)
       statistic: Spearman rho
  2. duration_trend_laa_sens_gt1m_eq1.tsv
       LAA, AF samples with ">1M" coded as 1  (n=16 expected)
       statistic: Spearman rho
  3. duration_trend_laa_sens_ordinal.tsv
       LAA, ALL samples; 3-tier ordinal SR / short_AF / long_AF
       statistic: Kendall tau
  4. duration_trend_lapv_sensitivity.tsv
       LAPV, AF samples with numeric duration (n=12 expected)
       statistic: Spearman rho   (mirrors analysis 1 on LAPV)

Permutation: MC, 10 000 draws, seed 20260418; BH within each file across
the sixteen NCP axis families (spec §3.5: LAA / LAPV families are separated; the
three LAA analyses above are also BH'd independently so their decision
boundaries do not leak into each other).

Spec: duration vs program score trends (LAA-centric primary + sensitivities).
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
IN_DIR  = ROOT / "Project_P4" / "GSE41177"
OUT_DIR = IN_DIR / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_MC_PERM = 10_000
MC_SEED   = 20_260_418


# ----------------------------------------------------------------------
# permutation helpers
# ----------------------------------------------------------------------

def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    r, _ = stats.spearmanr(x, y)
    return float(r)


def _kendall(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    t, _ = stats.kendalltau(x, y, variant="b")
    return float(t)


def permute_corr(scores: np.ndarray, labels: np.ndarray,
                 method: str, rng: np.random.Generator,
                 n_perm: int = N_MC_PERM) -> tuple[float, float, int]:
    """Return (statistic, p_two_sided, n_perm_used) for MC permutation."""
    stat_fn = _spearman if method == "spearman" else _kendall
    observed = stat_fn(scores, labels)
    if math.isnan(observed):
        return (observed, float("nan"), 0)
    hits = 1  # observed counted
    for _ in range(n_perm):
        perm = rng.permutation(len(labels))
        s = stat_fn(scores, labels[perm])
        if not math.isnan(s) and abs(s) >= abs(observed) - 1e-12:
            hits += 1
    p = hits / (n_perm + 1)
    return (observed, p, n_perm)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    mask = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    if mask.sum() == 0:
        return out
    pv = p[mask]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(pv)
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    qv = np.empty_like(pv)
    qv[order] = q
    out[mask] = qv
    return out


def sig_marker(p: float) -> str:
    if math.isnan(p):
        return ""
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "."
    return ""


# ----------------------------------------------------------------------
# per-analysis engine
# ----------------------------------------------------------------------

def run_one(scores_df: pd.DataFrame, sample_df: pd.DataFrame,
            analysis_name: str,
            tissue: str,
            include_mask_col: str,
            label_col: str,
            method: str = "spearman",
            label_transform=None) -> pd.DataFrame:
    """Core per-(program, tissue) loop; returns a per-program DataFrame."""
    rng = np.random.default_rng(MC_SEED)

    # samples entering this analysis
    sub_samples = sample_df[(sample_df["tissue"] == tissue) &
                            sample_df[include_mask_col]].copy()
    if label_transform is not None:
        sub_samples["_label"] = label_transform(sub_samples)
    else:
        sub_samples["_label"] = sub_samples[label_col]

    sub_samples = sub_samples.dropna(subset=["_label"])
    n = len(sub_samples)
    gsm_list = sub_samples["gsm"].tolist()
    labels = np.asarray(sub_samples["_label"], dtype=float)

    programs = (scores_df[["program_id", "program_label", "program_order",
                           "baseline_direction"]]
                .drop_duplicates()
                .sort_values("program_order")
                .reset_index(drop=True))

    rows = []
    for _, prog in programs.iterrows():
        prog_scores = (scores_df[(scores_df["tissue"] == tissue) &
                                 (scores_df["program_id"] == prog["program_id"])]
                       .set_index("gsm")["score"])
        sc = prog_scores.reindex(gsm_list).to_numpy(dtype=float)
        if math.isnan(np.sum(sc)) or n < 3:
            rho, p, _ = (float("nan"), float("nan"), 0)
        else:
            rho, p, _ = permute_corr(sc, labels, method, rng)
        rows.append(dict(
            program_id         = prog["program_id"],
            program_label      = prog["program_label"],
            program_order      = int(prog["program_order"]),
            baseline_direction = int(prog["baseline_direction"]),
            selection_pool     = "LA-biased" if prog["baseline_direction"] == 1 else "RA-biased",
            tissue             = tissue,
            analysis           = analysis_name,
            n_samples          = n,
            statistic_name     = method,
            statistic          = rho,
            p_perm             = p,
            n_perm             = N_MC_PERM,
        ))
    df = pd.DataFrame(rows)
    df["q_perm"]      = bh(df["p_perm"].to_numpy())
    df["sig_marker"]  = df["p_perm"].map(sig_marker)
    df["direction"]   = df["statistic"].map(
        lambda r: "NA" if math.isnan(r) else ("pos" if r > 0 else "neg")
    )
    return df.sort_values("program_order").reset_index(drop=True)


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main() -> None:
    scores = pd.read_csv(IN_DIR / "program_scores.tsv", sep="\t")
    meta   = pd.read_csv(IN_DIR / "sample_metadata.tsv", sep="\t")
    # data.table wrote True/False as TRUE/FALSE strings; coerce
    for c in ("keep_for_primary",):
        if meta[c].dtype == object:
            meta[c] = meta[c].astype(str).str.upper().isin(["TRUE", "T"])

    # --- auxiliary masks ---
    meta["_af_numeric"]  = (meta["rhythm"] == "AF") & meta["duration_months"].notna()
    meta["_af_sensA"]    = (meta["rhythm"] == "AF") & meta["duration_sensitivity_A"].notna()
    meta["_any_tier"]    = meta["duration_tier_B"].notna()

    analyses = [
        dict(
            name         = "laa_primary",
            tissue       = "LAA",
            include_col  = "_af_numeric",
            label_col    = "duration_months",
            method       = "spearman",
            out          = "duration_trend_laa_primary.tsv",
        ),
        dict(
            name         = "laa_sens_gt1m_eq1",
            tissue       = "LAA",
            include_col  = "_af_sensA",
            label_col    = "duration_sensitivity_A",
            method       = "spearman",
            out          = "duration_trend_laa_sens_gt1m_eq1.tsv",
        ),
        dict(
            name         = "laa_sens_ordinal",
            tissue       = "LAA",
            include_col  = "_any_tier",
            label_col    = "duration_tier_B",
            method       = "kendall",
            label_transform = lambda df: df["duration_tier_B"].map(
                {"SR": 0, "short_AF": 1, "long_AF": 2}),
            out          = "duration_trend_laa_sens_ordinal.tsv",
        ),
        dict(
            name         = "lapv_primary",
            tissue       = "LAPV",
            include_col  = "_af_numeric",
            label_col    = "duration_months",
            method       = "spearman",
            out          = "duration_trend_lapv_sensitivity.tsv",
        ),
    ]

    collected = []
    for a in analyses:
        df = run_one(scores, meta, a["name"],
                     tissue=a["tissue"],
                     include_mask_col=a["include_col"],
                     label_col=a["label_col"],
                     method=a["method"],
                     label_transform=a.get("label_transform"))
        out_path = OUT_DIR / a["out"]
        df.to_csv(out_path, sep="\t", index=False)
        print(f"wrote {out_path.relative_to(ROOT)}  (n_programs={len(df)}, "
              f"n_samples={int(df['n_samples'].iloc[0]) if len(df) else 0})")
        collected.append((a["name"], df))

    # stdout preview: rho across analyses, per program
    print("\n=== per-program ρ / τ across analyses (LAA only) ===")
    pivot = (pd.concat([d.assign(analysis=n) for n, d in collected
                        if d["tissue"].iloc[0] == "LAA"])
               .pivot(index=["program_order", "program_label", "baseline_direction"],
                      columns="analysis", values="statistic")
               .round(3))
    print(pivot.to_string())

    print("\n=== per-program LAPV ρ (sensitivity for 4D) ===")
    lapv = [d for n, d in collected if d["tissue"].iloc[0] == "LAPV"][0]
    print(lapv[["program_order", "program_label", "baseline_direction",
                "statistic", "p_perm", "q_perm", "sig_marker"]].round(3).to_string(index=False))

    print("\n=== significance tally per analysis ===")
    for n, d in collected:
        n_star = (d["p_perm"] < 0.05).sum()
        n_dot  = ((d["p_perm"] >= 0.05) & (d["p_perm"] < 0.10)).sum()
        n_q10  = (d["q_perm"] < 0.10).sum()
        print(f"  {n:22s}  n={int(d['n_samples'].iloc[0])}  "
              f"*={n_star}  .={n_dot}  (q<0.10)={n_q10}")


if __name__ == "__main__":
    main()
