#!/usr/bin/env python3
"""R4b - NCP_v2 Δsinglescore on GSE255612 LA snRNA, with per-program and"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

PERM_EXACT_LIMIT = 20_000
N_MC_PERM        = 10_000
MC_SEED          = 20_260_418
PVAL_THRESHOLD   = 0.10
MIN_DONORS       = 4

# 13-compartment pseudobulk (matches Fig 1f).  Two compartments
# (B_cell, epicardial) only have 1 Ctrl donor and are excluded.
MAIN_CTS = [
    "atrialCM", "fibroblast", "endothelial",
    "lymphatic_endothelial",
    "macrophage", "lymphocyte", "mast_cell",
    "pericyte", "vascular_smooth_muscle",
    "adipocyte", "neuronal",
]


def _stat(values: np.ndarray, is_af: np.ndarray) -> float:
    return float(values[is_af].mean() - values[~is_af].mean())


def permutation_p(values: np.ndarray, is_af: np.ndarray,
                    rng: np.random.Generator) -> tuple[float, str, int]:
    n = len(values)
    n_af = int(is_af.sum())
    if n_af == 0 or n_af == n:
        return (float("nan"), "none", 0)
    observed = _stat(values, is_af)
    n_combinations = math.comb(n, n_af)
    if n_combinations <= PERM_EXACT_LIMIT:
        from itertools import combinations
        hits = 0
        for af_idx in combinations(range(n), n_af):
            mask = np.zeros(n, dtype=bool)
            mask[list(af_idx)] = True
            if abs(_stat(values, mask)) >= abs(observed) - 1e-12:
                hits += 1
        return (hits / n_combinations, "exact", n_combinations)
    hits = 1
    for _ in range(N_MC_PERM):
        perm = rng.permutation(n)
        mask = np.zeros(n, dtype=bool)
        mask[perm[:n_af]] = True
        if abs(_stat(values, mask)) >= abs(observed) - 1e-12:
            hits += 1
    return (hits / (N_MC_PERM + 1), "mc", N_MC_PERM)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = q
    return np.clip(out, 0, 1)


def sig_marker(p: float) -> str:
    if math.isnan(p):
        return ""
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "."
    return ""


def per_unit_table(scores: pd.DataFrame, panel: pd.DataFrame, *,
                     unit_key: str, unit_label_col: str,
                     unit_order_col: str, extra_cols: list[str],
                     rng: np.random.Generator) -> pd.DataFrame:
    """Compute per-(unit × ct) Δ + perm p.  unit_key is either
    'program_id' or 'family_id'.  scores must contain a column matching
    unit_key plus 'cell_type', 'rhythm', 'donor_id', 'score'."""
    rows = []
    cols = [unit_key, unit_label_col]
    if unit_order_col not in cols:
        cols.append(unit_order_col)
    cols.extend(c for c in extra_cols if c not in cols)
    units = panel[cols].drop_duplicates(unit_key)
    units = units.sort_values(unit_order_col).reset_index(drop=True)
    for _, urow in units.iterrows():
        for ct in MAIN_CTS:
            sub = scores[(scores[unit_key] == urow[unit_key])
                          & (scores["cell_type"] == ct)]
            if unit_key == "family_id":
                # average across programs within the family for each donor
                sub = (sub.groupby(["donor_id", "rhythm"])["score"]
                          .mean().reset_index())
            ctrl = sub.loc[sub["rhythm"] == "Ctrl", "score"].to_numpy()
            af   = sub.loc[sub["rhythm"] == "AF",   "score"].to_numpy()
            n_ctrl, n_af = len(ctrl), len(af)
            included = n_ctrl >= MIN_DONORS and n_af >= MIN_DONORS
            row = {unit_key: urow[unit_key],
                    "label":      urow[unit_label_col],
                    "order":      int(urow[unit_order_col]),
                    "cell_type":  ct,
                    "n_ctrl_donors": n_ctrl,
                    "n_af_donors":   n_af,
                    "mean_score_ctrl": float(ctrl.mean()) if n_ctrl else float("nan"),
                    "mean_score_af":   float(af.mean())   if n_af   else float("nan"),
                    "sd_score_ctrl":   float(ctrl.std(ddof=1)) if n_ctrl > 1 else float("nan"),
                    "sd_score_af":     float(af.std(ddof=1))   if n_af   > 1 else float("nan"),
                    "dScore":          float(af.mean() - ctrl.mean()) if (n_ctrl and n_af) else float("nan"),
                    "cell_type_included": included,
                  }
            for c in extra_cols:
                row[c] = urow[c]
            if included:
                values = np.concatenate([ctrl, af])
                is_af = np.concatenate([np.zeros(n_ctrl, dtype=bool),
                                          np.ones(n_af, dtype=bool)])
                p, strat, n_perm = permutation_p(values, is_af, rng)
            else:
                p, strat, n_perm = float("nan"), "none", 0
            row.update(p_perm=p, perm_strategy=strat, n_perm=n_perm)
            rows.append(row)
    tbl = pd.DataFrame(rows)
    tbl["q_perm"] = float("nan")
    for ct in MAIN_CTS:
        mask = (tbl["cell_type"] == ct) & tbl["cell_type_included"]
        if mask.sum() > 0:
            tbl.loc[mask, "q_perm"] = bh(tbl.loc[mask, "p_perm"].to_numpy())
    tbl["sig_marker"] = tbl["p_perm"].map(sig_marker)
    return tbl


def pick_carrier(grp: pd.DataFrame) -> tuple[str, float, float]:
    cand = grp[grp["cell_type_included"] & (grp["q_perm"] < PVAL_THRESHOLD)]
    if cand.empty:
        return ("Unresolved", float("nan"), float("nan"))
    amax = cand["dScore"].abs().idxmax()
    winner = cand.loc[amax]
    wabs = abs(winner["dScore"])
    ties = cand[abs(cand["dScore"]) >= 0.9 * wabs]
    if len(ties) > 1:
        return ("/".join(sorted(ties["cell_type"])),
                  float(winner["dScore"]),
                  float(winner["q_perm"]))
    return (winner["cell_type"],
              float(winner["dScore"]),
              float(winner["q_perm"]))


def main() -> None:
    out_dir = ROOT / "results"
    score_path = ROOT / "Project_P3" / "GSE255612" / "program_scores_primary.tsv"
    panel_path = ROOT / "results" / "fig3_primary_program_le.tsv"
    fam_path   = ROOT / "results" / "fig1e_baseline_axis_families.tsv"

    scores = pd.read_csv(score_path, sep="\t").rename(
        columns={"compartment": "cell_type"})
    avail = set(scores["cell_type"].unique())
    missing = set(MAIN_CTS) - avail
    if missing:
        print(f"WARN: requested cell types not in scores: {sorted(missing)}")
    scores = scores[scores["cell_type"].isin(MAIN_CTS)].copy()
    print(f"scoring across {scores['cell_type'].nunique()} cell types: "
            f"{sorted(scores['cell_type'].unique())}")

    panel = pd.read_csv(panel_path, sep="\t")
    # drop any existing family_id, then re-attach from authoritative source
    panel = panel.drop(columns=[c for c in ("family_id",) if c in panel.columns])
    fam = pd.read_csv(fam_path, sep="\t")[["pathway", "family_id"]] \
              .rename(columns={"pathway": "program_id"})
    panel = panel.merge(fam, on="program_id", how="left")
    # attach program_label from scores (R3b stores it per-row)
    label_map = (scores.drop_duplicates("program_id")
                          [["program_id", "program_label"]])
    panel = panel.merge(label_map, on="program_id", how="left")
    if "family_id" not in scores.columns:
        scores = scores.merge(panel[["program_id", "family_id"]],
                                  on="program_id", how="left")

    # Source data — long per-donor singscores
    src = scores[["cell_type", "donor_id", "rhythm", "n_cells",
                    "program_id", "program_label", "program_order",
                    "family_id", "baseline_direction", "score"]].copy()
    src.columns = ["cell_type", "donor_id", "rhythm", "n_cells",
                     "program_id", "program_label", "program_order",
                     "family_id", "baseline_direction", "singscore"]
    src.to_csv(out_dir / "fig3_primary_source_data.tsv", sep="\t", index=False)
    print(f"wrote results/fig3_primary_source_data.tsv  ({len(src)} rows)")

    # per-program × cell-type
    rng = np.random.default_rng(MC_SEED)
    panel_for_prog = panel[["program_id", "program_label", "program_order",
                                 "family_id", "baseline_direction"]] \
                        .drop_duplicates("program_id")
    prog_tbl = per_unit_table(
        scores.rename(columns={}), panel_for_prog,
        unit_key="program_id",
        unit_label_col="program_label",
        unit_order_col="program_order",
        extra_cols=["family_id", "baseline_direction"],
        rng=rng)
    prog_tbl = prog_tbl.rename(columns={"label": "program_label",
                                              "order": "program_order"})
    prog_tbl.to_csv(out_dir / "fig3_primary_program_celltype_delta.tsv",
                       sep="\t", index=False)
    print(f"wrote results/fig3_primary_program_celltype_delta.tsv "
            f"({len(prog_tbl)} rows)")

    # Per-program major carrier
    car_rows = []
    for pid, grp in prog_tbl.groupby("program_id", sort=False):
        ct, dsc, qp = pick_carrier(grp)
        prow = grp.iloc[0]
        car_rows.append(dict(
            program_id=pid, program_label=prow["program_label"],
            program_order=int(prow["program_order"]),
            family_id=int(prow["family_id"]) if not pd.isna(prow["family_id"]) else -1,
            baseline_direction=int(prow["baseline_direction"]),
            AF_major_carrier_CT=ct,
            AF_major_carrier_dScore=dsc,
            AF_major_carrier_q_perm=qp))
    car_df = pd.DataFrame(car_rows).sort_values("program_order")
    car_df.to_csv(out_dir / "fig3_primary_program_carrier.tsv",
                     sep="\t", index=False)
    print(f"wrote results/fig3_primary_program_carrier.tsv  ({len(car_df)} rows)")

    # per-family × cell-type (NEW — bridge to Fig 2D)
    fam_panel = (panel.dropna(subset=["family_id"])
                       .groupby("family_id", as_index=False)
                       .agg(family_label=("program_label",
                                           lambda s: " / ".join(sorted(set(s))[:2])
                                           + (" ..." if len(set(s)) > 2 else "")),
                              direction=("baseline_direction", "mean"),
                              n_programs=("program_id", "nunique")))
    fam_panel["family_id"] = fam_panel["family_id"].astype(int)
    fam_panel["family_label"] = fam_panel.apply(
        lambda r: f"F{int(r['family_id'])}  {r['family_label']}", axis=1)
    fam_panel["direction_label"] = fam_panel["direction"].apply(
        lambda x: "LA" if x > 0 else "RA")
    fam_panel = fam_panel.sort_values("family_id").reset_index(drop=True)

    rng2 = np.random.default_rng(MC_SEED + 1)
    fam_tbl = per_unit_table(
        scores, fam_panel,
        unit_key="family_id",
        unit_label_col="family_label",
        unit_order_col="family_id",
        extra_cols=["direction_label", "n_programs"],
        rng=rng2)
    fam_tbl = fam_tbl.rename(columns={"label": "family_label",
                                            "order": "family_order"})
    fam_tbl.to_csv(out_dir / "fig3_primary_family_celltype_delta.tsv",
                      sep="\t", index=False)
    print(f"wrote results/fig3_primary_family_celltype_delta.tsv "
            f"({len(fam_tbl)} rows)")

    # stdout preview
    print("\n=== family × cell-type Δsinglescore (q_perm < 0.10) ===")
    for ct in MAIN_CTS:
        sub = fam_tbl[(fam_tbl["cell_type"] == ct)
                       & (fam_tbl["q_perm"] < 0.10)] \
                   .sort_values("dScore", key=lambda s: s.abs(),
                                  ascending=False)
        print(f"  {ct}:")
        if sub.empty:
            print("    (none)")
        else:
            for _, r in sub.iterrows():
                print(f"    F{int(r['family_id'])} {r['direction_label']} "
                        f"({r['n_programs']} progs)  "
                        f"Δ = {r['dScore']:+.4f}  "
                        f"p = {r['p_perm']:.3f}  q = {r['q_perm']:.3f}")

    print("\n=== per-program × ct  significance tally ===")
    for ct in MAIN_CTS:
        sub = prog_tbl[prog_tbl["cell_type"] == ct]
        n_star = int((sub["p_perm"] < 0.05).sum())
        n_dot  = int(((sub["p_perm"] >= 0.05) & (sub["p_perm"] < 0.10)).sum())
        print(f"  {ct:25s}  * = {n_star}, . = {n_dot}, total = {len(sub)}")


if __name__ == "__main__":
    main()
