#!/usr/bin/env python3
"""D4b - NCP_v2 program / family AF-duration correlation on GSE41177."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT    = Path(__file__).resolve().parent.parent
IN_DIR  = ROOT / "Project_P4" / "GSE41177"
OUT_DIR = IN_DIR / "results_primary"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_MC_PERM = 10_000
MC_SEED   = 20_260_418


def _spearman(x, y):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    r, _ = stats.spearmanr(x, y)
    return float(r)


def _kendall(x, y):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    t, _ = stats.kendalltau(x, y, variant="b")
    return float(t)


def permute_corr(scores, labels, method, rng, n_perm=N_MC_PERM):
    fn = _spearman if method == "spearman" else _kendall
    obs = fn(scores, labels)
    if math.isnan(obs):
        return (obs, float("nan"), 0)
    hits = 1
    for _ in range(n_perm):
        perm = rng.permutation(len(labels))
        s = fn(scores, labels[perm])
        if not math.isnan(s) and abs(s) >= abs(obs) - 1e-12:
            hits += 1
    return (obs, hits / (n_perm + 1), n_perm)


def bh(p):
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
    qv = np.empty_like(pv); qv[order] = q
    out[mask] = qv
    return out


def sig_marker(p):
    if math.isnan(p): return ""
    if p < 0.05:      return "*"
    if p < 0.10:      return "."
    return ""


def run_program_level(scores, samples, *, analysis_name, tissue,
                        include_mask, label_array, method):
    rng = np.random.default_rng(MC_SEED)
    sub_samples = samples[include_mask].copy()
    sub_samples["_label"] = label_array[include_mask.values]
    sub_samples = sub_samples.dropna(subset=["_label"])
    n = len(sub_samples)
    gsm_list = sub_samples["gsm"].tolist()
    labels = np.asarray(sub_samples["_label"], dtype=float)

    programs = (scores[["program_id", "program_label", "program_order",
                            "family_id", "baseline_direction"]]
                  .drop_duplicates("program_id")
                  .sort_values("program_order")
                  .reset_index(drop=True))

    rows = []
    for _, prog in programs.iterrows():
        ps = (scores[(scores["tissue"] == tissue)
                       & (scores["program_id"] == prog["program_id"])]
                  .set_index("gsm")["score"])
        sc = ps.reindex(gsm_list).to_numpy(dtype=float)
        if math.isnan(np.sum(sc)) or n < 3:
            rho, p, _ = (float("nan"), float("nan"), 0)
        else:
            rho, p, _ = permute_corr(sc, labels, method, rng)
        rows.append({
            "program_id":         prog["program_id"],
            "program_label":      prog["program_label"],
            "program_order":      int(prog["program_order"]),
            "family_id":          int(prog["family_id"])
                                       if not pd.isna(prog["family_id"]) else -1,
            "baseline_direction": int(prog["baseline_direction"]),
            "selection_pool":     "LA-biased" if prog["baseline_direction"] == 1
                                                  else "RA-biased",
            "tissue":             tissue,
            "analysis":           analysis_name,
            "n_samples":          n,
            "statistic_name":     method,
            "statistic":          rho,
            "p_perm":             p,
            "n_perm":             N_MC_PERM,
        })
    df = pd.DataFrame(rows)
    df["q_perm"]      = bh(df["p_perm"].to_numpy())
    df["sig_marker"]  = df["p_perm"].map(sig_marker)
    df["direction"]   = df["statistic"].map(
        lambda r: "NA" if math.isnan(r) else ("pos" if r > 0 else "neg"))
    return df.sort_values("program_order").reset_index(drop=True)


def run_family_level(scores, samples, *, analysis_name, tissue,
                       include_mask, label_array, method):
    """Score = mean of programs within each family per sample."""
    rng = np.random.default_rng(MC_SEED + 1)
    sub_samples = samples[include_mask].copy()
    sub_samples["_label"] = label_array[include_mask.values]
    sub_samples = sub_samples.dropna(subset=["_label"])
    n = len(sub_samples)
    gsm_list = sub_samples["gsm"].tolist()
    labels = np.asarray(sub_samples["_label"], dtype=float)

    families = (scores[["family_id"]]
                  .drop_duplicates()
                  .dropna()
                  .assign(family_id=lambda d: d["family_id"].astype(int))
                  .sort_values("family_id")
                  .reset_index(drop=True))

    # Build family meta: representative direction (mean), n_programs
    fam_meta = (scores.dropna(subset=["family_id"])
                       .assign(family_id=lambda d: d["family_id"].astype(int))
                       .groupby("family_id")
                       .agg(direction=("baseline_direction", "mean"),
                              n_programs=("program_id", "nunique"))
                       .reset_index())
    fam_meta["direction_label"] = fam_meta["direction"].apply(
        lambda x: "LA-biased" if x > 0 else "RA-biased")

    rows = []
    for _, fr in families.iterrows():
        fid = int(fr["family_id"])
        # mean score per sample across programs in family
        sub = scores[(scores["tissue"] == tissue)
                       & (scores["family_id"] == fid)]
        agg = (sub.groupby("gsm")["score"].mean())
        sc = agg.reindex(gsm_list).to_numpy(dtype=float)
        if math.isnan(np.sum(sc)) or n < 3:
            rho, p, _ = (float("nan"), float("nan"), 0)
        else:
            rho, p, _ = permute_corr(sc, labels, method, rng)
        meta_row = fam_meta[fam_meta["family_id"] == fid].iloc[0]
        rows.append({
            "family_id":       fid,
            "direction_label": meta_row["direction_label"],
            "n_programs":      int(meta_row["n_programs"]),
            "tissue":          tissue,
            "analysis":        analysis_name,
            "n_samples":       n,
            "statistic_name":  method,
            "statistic":       rho,
            "p_perm":          p,
            "n_perm":          N_MC_PERM,
        })
    df = pd.DataFrame(rows)
    df["q_perm"]     = bh(df["p_perm"].to_numpy())
    df["sig_marker"] = df["p_perm"].map(sig_marker)
    df["direction"]  = df["statistic"].map(
        lambda r: "NA" if math.isnan(r) else ("pos" if r > 0 else "neg"))
    return df.sort_values("family_id").reset_index(drop=True)


def main() -> None:
    scores = pd.read_csv(IN_DIR / "program_scores_primary.tsv", sep="\t")
    meta   = pd.read_csv(IN_DIR / "sample_metadata.tsv", sep="\t")
    if meta["keep_for_primary"].dtype == object:
        meta["keep_for_primary"] = meta["keep_for_primary"].astype(str) \
                                          .str.upper().isin(["TRUE", "T"])

    # auxiliary masks
    af_numeric = (meta["rhythm"] == "AF") & meta["duration_months"].notna()
    af_sensA   = (meta["rhythm"] == "AF") & meta["duration_sensitivity_A"].notna()
    any_tier   = meta["duration_tier_B"].notna()
    tier_label = meta["duration_tier_B"].map(
        {"SR": 0, "short_AF": 1, "long_AF": 2}).to_numpy()

    analyses = [
        dict(name="laa_primary",      tissue="LAA",
              mask=af_numeric & (meta["tissue"] == "LAA"),
              labels=meta["duration_months"].to_numpy(dtype=float),
              method="spearman",
              out="duration_trend_laa_primary.tsv"),
        dict(name="laa_sens_gt1m_eq1", tissue="LAA",
              mask=af_sensA & (meta["tissue"] == "LAA"),
              labels=meta["duration_sensitivity_A"].to_numpy(dtype=float),
              method="spearman",
              out="duration_trend_laa_sens_gt1m_eq1.tsv"),
        dict(name="laa_sens_ordinal", tissue="LAA",
              mask=any_tier & (meta["tissue"] == "LAA"),
              labels=tier_label.astype(float),
              method="kendall",
              out="duration_trend_laa_sens_ordinal.tsv"),
        dict(name="lapv_primary",      tissue="LAPV",
              mask=af_numeric & (meta["tissue"] == "LAPV"),
              labels=meta["duration_months"].to_numpy(dtype=float),
              method="spearman",
              out="duration_trend_lapv_sensitivity.tsv"),
    ]

    # program-level
    collected = []
    for a in analyses:
        df = run_program_level(scores, meta,
                                  analysis_name=a["name"], tissue=a["tissue"],
                                  include_mask=a["mask"],
                                  label_array=a["labels"], method=a["method"])
        out = OUT_DIR / a["out"]
        df.to_csv(out, sep="\t", index=False)
        n_star = int((df["p_perm"] < 0.05).sum())
        n_dot  = int(((df["p_perm"] >= 0.05) & (df["p_perm"] < 0.10)).sum())
        n_q10  = int((df["q_perm"] < 0.10).sum())
        print(f"wrote {out.relative_to(ROOT)}  "
                f"(n_progs={len(df)} n_samples={int(df['n_samples'].iloc[0])} "
                f"*={n_star} .={n_dot} q<0.10={n_q10})")
        collected.append((a["name"], df))

    # family-level (only for primary LAA + LAPV; sensitivities stay program-level)
    fam_analyses = [a for a in analyses if a["name"] in ("laa_primary",
                                                                "lapv_primary")]
    for a in fam_analyses:
        df = run_family_level(scores, meta,
                                 analysis_name=a["name"], tissue=a["tissue"],
                                 include_mask=a["mask"],
                                 label_array=a["labels"], method=a["method"])
        out = OUT_DIR / a["out"].replace("duration_trend_",
                                                "duration_trend_FAMILY_")
        df.to_csv(out, sep="\t", index=False)
        n_star = int((df["p_perm"] < 0.05).sum())
        n_dot  = int(((df["p_perm"] >= 0.05) & (df["p_perm"] < 0.10)).sum())
        print(f"wrote {out.relative_to(ROOT)}  "
                f"(n_families={len(df)} n_samples={int(df['n_samples'].iloc[0])} "
                f"*={n_star} .={n_dot})")

    # console summary
    print("\n=== program-level LAA primary - sig hits ===")
    laa = collected[0][1]
    sig = laa[laa["p_perm"] < 0.10].sort_values("p_perm")
    for _, r in sig.iterrows():
        print(f"  F{int(r['family_id'])} {r['program_label']:30s} "
                f"baseline={'LA' if r['baseline_direction']==1 else 'RA':2s}  "
                f"ρ={r['statistic']:+.3f}  p={r['p_perm']:.3f}  q={r['q_perm']:.3f}")

    print("\n=== family-level LAA primary - sorted by |ρ| ===")
    fam_laa = pd.read_csv(OUT_DIR / "duration_trend_FAMILY_laa_primary.tsv",
                              sep="\t")
    fam_laa["abs"] = fam_laa["statistic"].abs()
    print(fam_laa.sort_values("abs", ascending=False)
                  [["family_id", "direction_label", "n_programs",
                    "statistic", "p_perm", "q_perm", "sig_marker"]]
                  .to_string(index=False))


if __name__ == "__main__":
    main()
