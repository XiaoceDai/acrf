#!/usr/bin/env python3
"""
Emit manuscript-ready Fig 4 TSV bundles from locked v3 Fig 4 intermediates.

Prerequisite pipeline:
    scripts/v3/fig4_v3_duration_trend.py    → fig4_v3_* TSV + singscore long
    scripts/v3/fig4_v3_chronicity_class.py  → fig4_v3_chronicity_class.tsv

Outputs (results/)
------------------
fig4b_family_duration_correlations.tsv   family × region (LAA + LAPV): ρ, bootstrap CI,
                                         permutation p, scipy Spearman p (nominal), BH q within region,
                                         chronicity class from LAA primary thresholds.

fig4c_chronicity_classes.tsv             per-family summary for Fig 4C prose counts.

fig4c_chronicity_aggregate.tsv           counts + folded member lists per chronicity_class.

fig4def_exemplar_scatter.tsv              Point-level data matching Fig 4d–f panels (PDF): SR anchors +
                                         AF with finite duration_months only (excludes GEO '>1M' censored AF).
                                         Column used_in_duration_correlation marks the n=12 AF primary subset.

fig4def_exemplar_summary.tsv             selection rationale + ρ linkage (Fig 2–5 narrative hooks).

fig4_lapv_sensitivity.tsv                LAA vs LAPV ρ/p comparison per family.

fig4_duration_scores_long.tsv           sample_id × family mean singscore on Fig 3 LE sets (family-mean rank singscore).

Notes
-----
Scores are rank-based singscores on pathway LE genes (Fig 3 pipeline); aggregated as mean across programs → family-mean singscore (not chamber aCAIrel — cohort GSE41177 is single-region bulk samples).

SciPy Spearman p-values are asymptotic two-sided tests (nominal); permutation p_perm matches Fig 4 scripts (donor-label shuffle).

q_value: Benjamini–Hochberg FDR across families within each region (LAA vs LAPV separately), applied to nominal p_spearman.

fig4_lapv_sensitivity.tsv uses p_LAA / p_LAPV as the same permutation p-values as fig4_v3_FAMILY_* tables (not asymptotic Spearman).

chronicity_class on LAPV rows in fig4b repeats the LAA-primary classification (thresholds are defined on LAA).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]

LONG_TSV = ROOT / "results" / "fig4_v3_singscore_per_sample.tsv"
META_TSV = ROOT / "Project_P4" / "GSE41177" / "sample_metadata.tsv"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
LAA_FAM_TSV = ROOT / "results" / "fig4_v3_FAMILY_laa_primary.tsv"
LAPV_FAM_TSV = ROOT / "results" / "fig4_v3_FAMILY_lapv_sensitivity.tsv"
CHRONICITY_TSV = ROOT / "results" / "fig4_v3_chronicity_class.tsv"

OUT_FIG4B = ROOT / "results" / "fig4b_family_duration_correlations.tsv"
OUT_FIG4C_FAM = ROOT / "results" / "fig4c_chronicity_classes.tsv"
OUT_FIG4C_AGG = ROOT / "results" / "fig4c_chronicity_aggregate.tsv"
OUT_DEF_SCATTER = ROOT / "results" / "fig4def_exemplar_scatter.tsv"
OUT_DEF_SUM = ROOT / "results" / "fig4def_exemplar_summary.tsv"
OUT_LAPV = ROOT / "results" / "fig4_lapv_sensitivity.tsv"
OUT_LONG = ROOT / "results" / "fig4_duration_scores_long.tsv"


CASE_SUMMARY = {
    "Fig 4D": {
        "family_id": 9,
        "reason_selected": (
            "Structural-axis HS-GAG / ECM example; strongest robustness (Strong_trend, sens_same_sign 4/4)."
        ),
        "relationship_to_fig2_fig3": (
            "F9 RA-biased ECM programme arm in Fig 2; complements Fig 3 HS-GAG RA programmes where LA snRNA "
            "shows relay/extension — here bulk LAA links programme activity decline with AF duration."
        ),
    },
    "Fig 4E": {
        "family_id": 2,
        "reason_selected": (
            "Metabolic / calcium-handling anchor (Strong_trend); literature-facing progression narrative."
        ),
        "relationship_to_fig2_fig3": (
            "F2 LA-biased Mito Ca²⁺ transport in Fig 2 bulk discovery; duration slope matches mitochondrial "
            "calcium remodelling canon (Voigt / Heijman / Sun cited in figure footer)."
        ),
    },
    "Fig 4F": {
        "family_id": 15,
        "reason_selected": (
            "Weak-trend duration coupling chosen to bridge layers: same F15 NTRK2 family as Fig 3D cellular relay."
        ),
        "relationship_to_fig2_fig3": (
            "Explicit cross-link to Fig 3D relay (lymphocyte + endothelial nominal shifts); Fig 4f shows parallel "
            "bulk LAA duration decline — atlas prioritisation narrative despite marginal single-layer p-values."
        ),
    },
}


def bh_q(p: np.ndarray) -> np.ndarray:
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
    res[order] = np.clip(q, 0, 1)
    out[finite] = res
    return out


def family_xy_duration(long_df: pd.DataFrame, meta: pd.DataFrame, fid: int, tissue: str):
    """AF donors with numeric duration (matches Fig 4 pipeline)."""
    sub = long_df[(long_df["family_id"] == fid) & (long_df["tissue"] == tissue)].copy()
    if sub.empty:
        return None, None
    fam = sub.groupby("gsm", as_index=False)["score"].mean()
    m = fam.merge(meta[["gsm", "rhythm", "duration_months"]], on="gsm", how="left")
    m["duration_months"] = pd.to_numeric(m["duration_months"], errors="coerce")
    af = m[(m["rhythm"] == "AF") & m["duration_months"].notna()].copy()
    af = af[np.isfinite(af["duration_months"]) & np.isfinite(af["score"])]
    if len(af) < 4:
        return None, None
    x = af["duration_months"].to_numpy(dtype=float)
    y = af["score"].to_numpy(dtype=float)
    return x, y


def spearman_nominal(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if x is None or len(x) < 4:
        return float("nan"), float("nan")
    r = spearmanr(x, y)
    return float(r.correlation), float(r.pvalue)


def linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 3:
        return float("nan"), float("nan")
    coef = np.polyfit(x, y, 1)
    return float(coef[0]), float(coef[1])


def main() -> None:
    missing = [p for p in [LONG_TSV, META_TSV, LABELS_TSV, LAA_FAM_TSV, LAPV_FAM_TSV, CHRONICITY_TSV] if not p.exists()]
    if missing:
        sys.exit(f"Missing inputs (run Fig 4 pipeline first): {missing}")

    long_df = pd.read_csv(LONG_TSV, sep="\t")
    meta = pd.read_csv(META_TSV, sep="\t")
    meta["duration_months"] = pd.to_numeric(meta["duration_months"], errors="coerce")
    labs = pd.read_csv(LABELS_TSV, sep="\t")[["family_id", "family_label", "axis_group"]].drop_duplicates(
        subset=["family_id"], keep="first"
    )

    laa = pd.read_csv(LAA_FAM_TSV, sep="\t")
    lapv = pd.read_csv(LAPV_FAM_TSV, sep="\t")
    chronicity = pd.read_csv(CHRONICITY_TSV, sep="\t")

    laa = laa.merge(labs, on="family_id", how="left")
    lapv = lapv.merge(labs, on="family_id", how="left")
    chronicity = chronicity.merge(labs, on="family_id", how="left")

    # --- Spearman nominal p per region for BH ---
    p_nom_laa: dict[int, float] = {}
    p_nom_lapv: dict[int, float] = {}
    for fid in laa["family_id"].unique():
        x, y = family_xy_duration(long_df, meta, int(fid), "LAA")
        _, pv = spearman_nominal(x, y) if x is not None else (np.nan, np.nan)
        p_nom_laa[int(fid)] = pv
    for fid in lapv["family_id"].unique():
        x, y = family_xy_duration(long_df, meta, int(fid), "LAPV")
        _, pv = spearman_nominal(x, y) if x is not None else (np.nan, np.nan)
        p_nom_lapv[int(fid)] = pv

    fids_laa_sorted = sorted({int(x) for x in laa["family_id"].tolist()})
    fids_lapv_sorted = sorted({int(x) for x in lapv["family_id"].tolist()})
    q_vec_laa = bh_q(np.array([p_nom_laa.get(f, np.nan) for f in fids_laa_sorted]))
    q_vec_lapv = bh_q(np.array([p_nom_lapv.get(f, np.nan) for f in fids_lapv_sorted]))
    q_map_laa = dict(zip(fids_laa_sorted, q_vec_laa))
    q_map_lapv = dict(zip(fids_lapv_sorted, q_vec_lapv))

    cls_map = chronicity.set_index("family_id")["chronicity_class"].to_dict()
    label_map = labs.set_index("family_id")["family_label"].to_dict()

    rows_b = []
    for i, r in laa.iterrows():
        fid = int(r["family_id"])
        rows_b.append({
            "family_id": fid,
            "family_label": r["family_label"],
            "axis_group": r["axis_group"],
            "region": "LAA",
            "n_donors": int(r["n_samples"]),
            "rho": float(r["rho_spearman"]),
            "ci_low": float(r["bootstrap_ci_lo"]),
            "ci_high": float(r["bootstrap_ci_hi"]),
            "p_perm": float(r["p_perm"]),
            "p_spearman": p_nom_laa.get(fid, float("nan")),
            "q_value": float(q_map_laa.get(fid, float("nan"))),
            "chronicity_class": cls_map.get(fid, ""),
        })
    lapv_lookup = lapv.set_index("family_id")
    for fid in sorted(lapv_lookup.index.astype(int)):
        r = lapv_lookup.loc[fid]
        rows_b.append({
            "family_id": int(fid),
            "family_label": r["family_label"],
            "axis_group": r["axis_group"],
            "region": "LAPV",
            "n_donors": int(r["n_samples"]),
            "rho": float(r["rho_spearman"]),
            "ci_low": float(r["bootstrap_ci_lo"]),
            "ci_high": float(r["bootstrap_ci_hi"]),
            "p_perm": float(r["p_perm"]),
            "p_spearman": p_nom_lapv.get(int(fid), float("nan")),
            "q_value": float(q_map_lapv.get(int(fid), float("nan"))),
            "chronicity_class": cls_map.get(int(fid), ""),
        })

    fig4b = pd.DataFrame(rows_b).sort_values(["region", "family_id"]).reset_index(drop=True)
    OUT_FIG4B.parent.mkdir(parents=True, exist_ok=True)
    fig4b.to_csv(OUT_FIG4B, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {OUT_FIG4B.relative_to(ROOT)} ({len(fig4b)} rows)")

    # --- Fig 4C per-family ---
    c_rows = []
    for _, r in chronicity.iterrows():
        fid = int(r["family_id"])
        prow = laa[laa["family_id"] == fid].iloc[0]
        c_rows.append({
            "family_id": fid,
            "family_label": r["family_label"],
            "rho": float(prow["rho_spearman"]),
            "p_perm": float(prow["p_perm"]),
            "ci_low": float(prow["bootstrap_ci_lo"]),
            "ci_high": float(prow["bootstrap_ci_hi"]),
            "chronicity_class": r["chronicity_class"],
            "direction": r["direction"],
            "p_spearman": p_nom_laa.get(fid, float("nan")),
            "q_value": float(q_map_laa.get(fid, float("nan"))),
        })
    fig4c_fam = pd.DataFrame(c_rows).sort_values(["chronicity_class", "family_id"]).reset_index(drop=True)
    fig4c_fam.to_csv(OUT_FIG4C_FAM, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {OUT_FIG4C_FAM.relative_to(ROOT)} ({len(fig4c_fam)} rows)")

    # aggregate chronicity summary
    agg_rows = []
    for cls_name, grp in fig4c_fam.groupby("chronicity_class"):
        ids = ",".join(str(int(x)) for x in sorted(grp["family_id"].tolist()))
        flabs = "; ".join(grp.sort_values("family_id")["family_label"].astype(str).tolist())
        dir_vc = grp["direction"].value_counts()
        dir_sum = "; ".join(f"{k}: {int(v)}" for k, v in dir_vc.items())
        agg_rows.append({
            "chronicity_class": cls_name,
            "n_families": int(len(grp)),
            "family_ids": ids,
            "family_labels": flabs,
            "direction_summary": dir_sum,
        })
    pd.DataFrame(agg_rows).to_csv(OUT_FIG4C_AGG, sep="\t", index=False)
    print(f"wrote {OUT_FIG4C_AGG.relative_to(ROOT)} ({len(agg_rows)} rows)")

    # --- LAPV sensitivity wide ---
    lapv_rows = []
    for fid in sorted(laa["family_id"].unique()):
        pr = laa[laa["family_id"] == fid].iloc[0]
        lr = lapv[lapv["family_id"] == fid].iloc[0]
        rho_a, rho_v = float(pr["rho_spearman"]), float(lr["rho_spearman"])
        conc = np.sign(rho_a) == np.sign(rho_v)
        if not np.isfinite(rho_a) or not np.isfinite(rho_v):
            sclass = "incomplete"
        elif conc:
            sclass = "direction_concordant"
        else:
            sclass = "LAPV_opposite_sign"
        lapv_rows.append({
            "family_id": int(fid),
            "family_label": pr["family_label"],
            "rho_LAA": rho_a,
            "p_LAA": float(pr["p_perm"]),
            "rho_LAPV": rho_v,
            "p_LAPV": float(lr["p_perm"]),
            "direction_concordant": conc,
            "sensitivity_class": sclass,
        })
    pd.DataFrame(lapv_rows).to_csv(OUT_LAPV, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {OUT_LAPV.relative_to(ROOT)} ({len(lapv_rows)} rows)")

    # --- duration scores long (family mean per gsm) ---
    # singscore long already has `tissue` (LAA / LAPV); avoid duplicate tissue_* merge columns.
    meta_sub = meta[["gsm", "donor_id", "rhythm", "duration_months"]].copy()
    prog_long = long_df.merge(meta_sub, on="gsm", how="left")
    fam_long = (
        prog_long.groupby(["gsm", "tissue", "donor_id", "rhythm", "duration_months", "family_id"], as_index=False)
        ["score"].mean()
    )
    fam_long = fam_long.merge(labs[["family_id", "family_label"]], on="family_id", how="left")
    fam_long["region"] = fam_long["tissue"]
    fam_long["score_type"] = "family_mean_LE_singscore"
    fam_long = fam_long.rename(columns={"gsm": "sample_id"})
    out_long = fam_long[
        ["donor_id", "sample_id", "region", "duration_months", "family_id", "family_label", "score_type", "score"]
    ].sort_values(["family_id", "region", "sample_id"])
    out_long.to_csv(OUT_LONG, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_LONG.relative_to(ROOT)} ({len(out_long)} rows)")

    # --- Fig 4def scatter + summary ---
    scatter_rows = []
    sum_rows = []
    cls_idx = chronicity.set_index("family_id")

    for panel, info in CASE_SUMMARY.items():
        fid = info["family_id"]

        sub = long_df[(long_df["family_id"] == fid) & (long_df["tissue"] == "LAA")].copy()
        fam_scores = sub.groupby("gsm", as_index=False)["score"].mean()
        merged = fam_scores.merge(meta[["gsm", "donor_id", "rhythm", "duration_months"]], on="gsm", how="left")
        merged["duration_months"] = pd.to_numeric(merged["duration_months"], errors="coerce")
        af = merged[(merged["rhythm"] == "AF") & merged["duration_months"].notna()].copy()
        if len(af) >= 3:
            slope, icept = linear_fit(af["duration_months"].to_numpy(float), af["score"].to_numpy(float))
        else:
            slope, icept = float("nan"), float("nan")

        cr = cls_idx.loc[fid]
        plot_mask = (merged["rhythm"] == "SR") | (
            (merged["rhythm"] == "AF") & merged["duration_months"].notna()
        )
        plot_df = merged.loc[plot_mask].copy()
        for _, row in plot_df.iterrows():
            d = row["duration_months"]
            is_af_numeric = (row["rhythm"] == "AF") and pd.notna(d)
            fitted = (
                slope * float(d) + icept
                if is_af_numeric and np.isfinite(slope)
                else float("nan")
            )
            scatter_rows.append({
                "panel": panel,
                "sample_id": row["gsm"],
                "donor_id": row["donor_id"],
                "region": "LAA",
                "rhythm": row["rhythm"],
                "duration_months": d if pd.notna(d) else "",
                "family_id": fid,
                "family_label": label_map[fid],
                "family_score": float(row["score"]),
                "fitted_value": fitted,
                "used_in_duration_correlation": is_af_numeric,
            })

        sum_rows.append({
            "panel": panel,
            "family_id": fid,
            "family_label": label_map[fid],
            "reason_selected": info["reason_selected"],
            "rho": float(cr["rho_primary"]),
            "ci_low": float(cr["ci_lo_primary"]),
            "ci_high": float(cr["ci_hi_primary"]),
            "p_perm": float(cr["p_perm_primary"]),
            "chronicity_class": str(cr["chronicity_class"]),
            "relationship_to_fig2_fig3": info["relationship_to_fig2_fig3"],
        })

    pd.DataFrame(scatter_rows).to_csv(OUT_DEF_SCATTER, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {OUT_DEF_SCATTER.relative_to(ROOT)} ({len(scatter_rows)} rows)")
    pd.DataFrame(sum_rows).to_csv(OUT_DEF_SUM, sep="\t", index=False)
    print(f"wrote {OUT_DEF_SUM.relative_to(ROOT)} ({len(sum_rows)} rows)")


if __name__ == "__main__":
    main()
