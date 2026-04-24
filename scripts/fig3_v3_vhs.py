#!/usr/bin/env python3
"""
Variance-Heterogeneity Score (VHS) on GSE255612 LA snRNA pseudobulk.

Concept
-------
A pseudobulk score is a single scalar per (donor, compartment, family).
If a compartment in AF contains an emerging sub-population that the bulk-level
score does NOT show as a mean shift, it should still leak into the
between-donor variance: AF donors will spread more widely than Ctrl donors
because the unobserved subpopulation fraction varies between donors.

VHS captures this:

    VHS(family, compartment) = SD_AF(family_singscore) / SD_Ctrl(family_singscore)
    log2_VHS                  = log2(VHS)

Significance: permutation of (Ctrl/AF) labels (5000 shuffles, two-sided
permutation p on |log2_VHS|).

Interpretation
--------------
    log2_VHS > +0.5   AF spread > Ctrl spread by >= ~1.4x
                       --> candidate sub-population emergence in AF
    log2_VHS < -0.5   AF spread < Ctrl spread (homogenisation in AF)
    |log2_VHS| < 0.3  no detectable heterogeneity change

Inputs
------
    Project_P3/GSE255612/pseudobulk_full/<comp>/expr.tsv + samples.tsv
    results/fig3_v3_program_le.tsv
    results/innovation_ncp_v3_axis_families.tsv

Outputs
-------
    results/fig3_v3_VHS_per_family_ct.tsv
        family_id, compartment, n_LE_genes, n_ctrl, n_af,
        sd_af, sd_ctrl, VHS, log2_VHS, p_perm, q_perm, sig_marker
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

PB_DIR = ROOT / "Project_P3" / "GSE255612" / "pseudobulk_full"
LE_TSV = ROOT / "results" / "fig3_v3_program_le.tsv"
FAM_TSV = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
OUT_TSV = ROOT / "results" / "fig3_v3_VHS_per_family_ct.tsv"

COMPARTMENTS = [
    "atrialCM", "fibroblast", "endothelial", "macrophage", "lymphocyte",
    "pericyte", "vascular_smooth_muscle", "adipocyte", "neuronal",
]

N_PERM = 5000
SEED = 20_260_418
MIN_DONORS = 4
MIN_LE = 5


def normalize_expr(expr):
    col_sums = expr.sum(axis=0).replace(0.0, np.nan)
    cpm = expr.div(col_sums, axis=1) * 1e6
    return np.log1p(cpm)


def per_sample_rank(expr_log):
    ranks = expr_log.rank(axis=0, method="average")
    return ranks / ranks.shape[0]


def bh(p):
    p = np.asarray(p, dtype=float)
    finite = np.isfinite(p)
    out = np.full_like(p, np.nan)
    if finite.sum() == 0:
        return out
    pf = p[finite]
    m = len(pf)
    order = np.argsort(pf)
    ranked = pf[order]
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    res = np.empty_like(q)
    res[order] = q
    out[finite] = np.clip(res, 0, 1)
    return out


def sig_marker(p):
    if pd.isna(p):
        return ""
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "."
    return ""


def main():
    rng = np.random.default_rng(SEED)
    le_df = pd.read_csv(LE_TSV, sep="\t")
    le_df["family_id"] = le_df["family_id"].astype(int)

    # Family-level union of LE genes across programs in family
    family_le: dict[int, list[str]] = {}
    for fid, g in le_df.groupby("family_id"):
        all_genes: set[str] = set()
        for _, r in g.iterrows():
            for gn in str(r["leading_edge_genes"]).split(","):
                gn = gn.strip()
                if gn:
                    all_genes.add(gn)
        family_le[int(fid)] = sorted(all_genes)
    print(f"Loaded LE for {len(family_le)} families "
          f"(median LE size = {int(np.median([len(v) for v in family_le.values()]))})")

    rows = []
    for comp in COMPARTMENTS:
        expr_path = PB_DIR / comp / "expr.tsv"
        samp_path = PB_DIR / comp / "samples.tsv"
        if not expr_path.exists():
            print(f"[SKIP] {comp} no expr")
            continue
        expr = pd.read_csv(expr_path, sep="\t", index_col=0)
        meta = pd.read_csv(samp_path, sep="\t")
        expr_log = normalize_expr(expr)
        ranks = per_sample_rank(expr_log)

        ctrl_donors = meta.loc[meta["rhythm"] == "Ctrl", "donor_id"].astype(str).tolist()
        af_donors = meta.loc[meta["rhythm"] == "AF", "donor_id"].astype(str).tolist()
        ctrl_cols = [d for d in ctrl_donors if d in ranks.columns]
        af_cols = [d for d in af_donors if d in ranks.columns]
        if len(ctrl_cols) < MIN_DONORS or len(af_cols) < MIN_DONORS:
            continue

        n_a, n_c = len(af_cols), len(ctrl_cols)
        all_cols = af_cols + ctrl_cols  # fixed order

        for fid, le_genes in family_le.items():
            le_in = [g for g in le_genes if g in ranks.index]
            if len(le_in) < MIN_LE:
                continue
            le_ranks = ranks.loc[le_in]
            # family-level singscore per donor: mean of LE gene norm-ranks - 0.5
            sample_scores = le_ranks.mean(axis=0) - 0.5
            af_vals = sample_scores[af_cols].dropna().values.astype(float)
            ctrl_vals = sample_scores[ctrl_cols].dropna().values.astype(float)
            if len(af_vals) < MIN_DONORS or len(ctrl_vals) < MIN_DONORS:
                continue
            sd_af = float(np.std(af_vals, ddof=1))
            sd_ctrl = float(np.std(ctrl_vals, ddof=1))
            if sd_af <= 0 or sd_ctrl <= 0:
                continue
            vhs = sd_af / sd_ctrl
            log2_vhs = float(np.log2(vhs))
            obs = abs(log2_vhs)

            # Permutation
            scores_arr = sample_scores[all_cols].values.astype(float)
            hits = 1
            for _ in range(N_PERM):
                perm = rng.permutation(len(all_cols))
                af_idx = perm[:n_a]
                ct_idx = perm[n_a:]
                a_v = scores_arr[af_idx]
                c_v = scores_arr[ct_idx]
                sa = float(np.std(a_v, ddof=1))
                sc = float(np.std(c_v, ddof=1))
                if sa <= 0 or sc <= 0:
                    continue
                lv = float(np.log2(sa / sc))
                if abs(lv) >= obs - 1e-12:
                    hits += 1
            p_perm = hits / (N_PERM + 1)

            rows.append({
                "family_id": int(fid),
                "compartment": comp,
                "n_LE_genes": int(len(le_in)),
                "n_ctrl": int(len(ctrl_vals)),
                "n_af": int(len(af_vals)),
                "sd_ctrl": round(sd_ctrl, 6),
                "sd_af": round(sd_af, 6),
                "VHS": round(vhs, 4),
                "log2_VHS": round(log2_vhs, 4),
                "p_perm": p_perm,
            })

    out = pd.DataFrame(rows)
    # BH q within compartment
    out["q_perm"] = np.nan
    for comp in COMPARTMENTS:
        mask = out["compartment"] == comp
        if mask.any():
            out.loc[mask, "q_perm"] = bh(out.loc[mask, "p_perm"].to_numpy())
    out["sig_marker"] = out["p_perm"].map(sig_marker)
    out = out.sort_values(["family_id", "compartment"]).reset_index(drop=True)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(out)} rows)")
    print(f"\n=== Significant VHS (p<0.05) ===")
    sig = out[out["p_perm"] < 0.05]
    if sig.empty:
        print("  none")
    else:
        for _, r in sig.iterrows():
            direction = "AF more variable" if r["log2_VHS"] > 0 else "AF more uniform"
            print(f"  F{r['family_id']:2d} x {r['compartment']:25s}  log2 VHS = {r['log2_VHS']:+.2f}  p={r['p_perm']:.3f}  ({direction})")
    print(f"\n=== Trending VHS (0.05<=p<0.10) ===")
    trend = out[(out["p_perm"] >= 0.05) & (out["p_perm"] < 0.10)]
    if trend.empty:
        print("  none")
    else:
        for _, r in trend.iterrows():
            direction = "AF more variable" if r["log2_VHS"] > 0 else "AF more uniform"
            print(f"  F{r['family_id']:2d} x {r['compartment']:25s}  log2 VHS = {r['log2_VHS']:+.2f}  p={r['p_perm']:.3f}  ({direction})")


if __name__ == "__main__":
    main()
