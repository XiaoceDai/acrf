#!/usr/bin/env python3
"""
R4 - Δsinglescore per (program, cell_type) with donor-level permutation.

Spec: LA pseudobulk Δsingscore carrier classification (this module).

For GSE255612 (18 AF + 16 Ctrl donors), exact permutation requires
C(34, 18) ≈ 2.3e9 enumerations — far above the PERM_EXACT_LIMIT = 20,000,
so we always use Monte Carlo 10,000 draws with a fixed seed (per Fig 2
convention).

Inputs
    Project_P3/GSE255612/program_scores.tsv     (R3)
    results/figure3_program_panel.tsv           (R3)

Outputs
    results/af_la_program_celltype_delta.tsv    main stats table
    results/figure3_source_data.tsv             long per-donor scores
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

# Locked constants (carrier / permutation gates).
PERM_EXACT_LIMIT = 20_000
N_MC_PERM        = 10_000
MC_SEED          = 20_260_418
PVAL_THRESHOLD   = 0.10
MIN_DONORS       = 4        # spec §3.6: ≥ 4 donors on each side per ct


# ----------------------------------------------------------------------
# permutation
# ----------------------------------------------------------------------

def _stat(values: np.ndarray, is_af: np.ndarray) -> float:
    return float(values[is_af].mean() - values[~is_af].mean())


def permutation_p(values: np.ndarray, is_af: np.ndarray,
                  rng: np.random.Generator) -> Tuple[float, str, int]:
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

    # Monte Carlo
    hits = 1  # observed counts itself
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


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main() -> None:
    score_path = ROOT / "Project_P3" / "GSE255612" / "program_scores.tsv"
    panel_path = ROOT / "results" / "figure3_program_panel.tsv"
    out_dir    = ROOT / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    scores = pd.read_csv(score_path, sep="\t")
    panel  = pd.read_csv(panel_path, sep="\t")

    # Main compartments only (spec §3.2).
    main_cts = ["atrialCM", "fibroblast", "endothelial_or_vascular", "immune"]
    scores = scores[~scores["is_supp"]].copy()
    stopifnot_cts = set(scores["group"].unique())
    assert stopifnot_cts == set(main_cts), (
        f"unexpected compartments in program_scores: {stopifnot_cts}")

    # ------------------------------------------------------------------
    # figure3_source_data.tsv: long per-donor scores
    # ------------------------------------------------------------------
    src_cols = [
        "donor_id", "biosample_id", "rhythm", "sex", "n_cells",
        "program_id", "program_label", "program_order", "baseline_direction",
    ]
    src = scores.rename(columns={"group": "cell_type"})
    src = src[["cell_type"] + src_cols + ["score"]]
    src.columns = ["cell_type"] + src_cols + ["singscore"]
    src_out = out_dir / "figure3_source_data.tsv"
    src.to_csv(src_out, sep="\t", index=False)
    print(f"wrote {src_out.relative_to(ROOT)}  ({len(src)} rows)")

    # ------------------------------------------------------------------
    # per (program, cell_type) stats
    # ------------------------------------------------------------------
    panel = panel.sort_values("program_order").reset_index(drop=True)

    # attach Fig 2 provenance so the wide table carries selection_reason
    rewiring = pd.read_csv(ROOT / "results" / "af_chamber_rewiring_programs.tsv",
                           sep="\t")
    panel = panel.merge(
        rewiring[["program_id"]],
        on="program_id", how="left",
    )

    rng = np.random.default_rng(MC_SEED)
    rows = []
    for _, prog in panel.iterrows():
        for ct in main_cts:
            sub = scores[(scores["program_id"] == prog["program_id"]) &
                         (scores["group"] == ct)]
            ctrl = sub.loc[sub["rhythm"] == "Ctrl", "score"].to_numpy()
            af   = sub.loc[sub["rhythm"] == "AF",   "score"].to_numpy()
            n_ctrl, n_af = len(ctrl), len(af)

            row = dict(
                program_id         = prog["program_id"],
                program_label      = prog["program_label"],
                program_order      = int(prog["program_order"]),
                baseline_direction = int(prog["baseline_direction"]),
                selection_reason   = prog["selection_reason"],
                cell_type          = ct,
                n_ctrl_donors      = n_ctrl,
                n_af_donors        = n_af,
                mean_score_ctrl    = float(ctrl.mean()) if n_ctrl else float("nan"),
                mean_score_af      = float(af.mean())   if n_af   else float("nan"),
                sd_score_ctrl      = float(ctrl.std(ddof=1)) if n_ctrl > 1 else float("nan"),
                sd_score_af        = float(af.std(ddof=1))   if n_af   > 1 else float("nan"),
                dScore             = (
                    float(af.mean() - ctrl.mean()) if n_ctrl and n_af else float("nan")
                ),
            )
            # median cell count per donor in this program × ct (for QC)
            if n_ctrl and n_af:
                row["n_ctrl_cells_median"] = int(np.median(
                    sub.loc[sub["rhythm"] == "Ctrl", "n_cells"]))
                row["n_af_cells_median"] = int(np.median(
                    sub.loc[sub["rhythm"] == "AF", "n_cells"]))
            else:
                row["n_ctrl_cells_median"] = 0
                row["n_af_cells_median"]   = 0

            row["cell_type_included"] = (n_ctrl >= MIN_DONORS and n_af >= MIN_DONORS)

            if row["cell_type_included"]:
                values = np.concatenate([ctrl, af])
                is_af  = np.concatenate([np.zeros(n_ctrl, dtype=bool),
                                          np.ones(n_af, dtype=bool)])
                p, strategy, n_perm = permutation_p(values, is_af, rng)
            else:
                p, strategy, n_perm = float("nan"), "none", 0
            row.update(p_perm=p, perm_strategy=strategy, n_perm=n_perm)
            rows.append(row)

    tbl = pd.DataFrame(rows)

    # BH q-value within each cell type, across panel programs
    tbl["q_perm"] = float("nan")
    for ct in main_cts:
        mask = (tbl["cell_type"] == ct) & tbl["cell_type_included"]
        ps = tbl.loc[mask, "p_perm"].to_numpy()
        qs = bh(ps)
        tbl.loc[mask, "q_perm"] = qs

    tbl["sig_marker"] = tbl["p_perm"].map(sig_marker)

    # Add ratio for downstream reading (may be NaN when sign differs or near-zero)
    # ratio_af_over_ctrl = mean_af / mean_ctrl  (singscore can be negative; use when both same sign)
    def _ratio(r):
        c, a = r["mean_score_ctrl"], r["mean_score_af"]
        if pd.isna(c) or pd.isna(a) or abs(c) < 1e-6 or np.sign(c) != np.sign(a):
            return float("nan")
        return a / c
    tbl["ratio_af_over_ctrl"] = tbl.apply(_ratio, axis=1)

    # ------------------------------------------------------------------
    # major_carrier (spec §3.7) -- computed at the program level
    # ------------------------------------------------------------------
    def _pick_carrier(grp: pd.DataFrame) -> Tuple[str, float, float]:
        cand = grp[grp["cell_type_included"] & (grp["q_perm"] < PVAL_THRESHOLD)]
        if cand.empty:
            return ("Unresolved", float("nan"), float("nan"))
        amax = cand["dScore"].abs().idxmax()
        winner = cand.loc[amax]
        # tie check: any other ct within 10 % of winner's |dScore|?
        wabs = abs(winner["dScore"])
        ties = cand[abs(cand["dScore"]) >= 0.9 * wabs]
        if len(ties) > 1:
            return ("/".join(sorted(ties["cell_type"])),
                    float(winner["dScore"]),
                    float(winner["q_perm"]))
        return (winner["cell_type"],
                float(winner["dScore"]),
                float(winner["q_perm"]))

    carrier_rows = []
    for pid, grp in tbl.groupby("program_id", sort=False):
        ct, dsc, qp = _pick_carrier(grp)
        prog = grp.iloc[0]
        carrier_rows.append(dict(
            program_id          = pid,
            program_label       = prog["program_label"],
            program_order       = int(prog["program_order"]),
            baseline_direction  = int(prog["baseline_direction"]),
            selection_reason    = prog["selection_reason"],
            AF_major_carrier_CT = ct,
            AF_major_carrier_dScore = dsc,
            AF_major_carrier_q_perm = qp,
        ))
    carrier_df = pd.DataFrame(carrier_rows).sort_values("program_order")

    # ------------------------------------------------------------------
    # write outputs
    # ------------------------------------------------------------------
    col_order = [
        "program_id", "program_label", "program_order", "baseline_direction",
        "selection_reason", "cell_type",
        "n_ctrl_donors", "n_af_donors",
        "n_ctrl_cells_median", "n_af_cells_median",
        "cell_type_included",
        "mean_score_ctrl", "mean_score_af", "sd_score_ctrl", "sd_score_af",
        "dScore", "ratio_af_over_ctrl",
        "perm_strategy", "n_perm", "p_perm", "q_perm", "sig_marker",
    ]
    tbl = tbl[col_order].sort_values(["program_order", "cell_type"])

    delta_out = out_dir / "af_la_program_celltype_delta.tsv"
    tbl.to_csv(delta_out, sep="\t", index=False)
    print(f"wrote {delta_out.relative_to(ROOT)}  ({len(tbl)} rows, {tbl.shape[1]} cols)")

    carrier_out = out_dir / "af_la_major_carrier_af.tsv"
    carrier_df.to_csv(carrier_out, sep="\t", index=False)
    print(f"wrote {carrier_out.relative_to(ROOT)}  ({len(carrier_df)} rows)")

    # ------------------------------------------------------------------
    # stdout preview
    # ------------------------------------------------------------------
    print("\n=== major carrier per program ===")
    for _, r in carrier_df.iterrows():
        print(f"  {r['program_order']:2d} {r['program_label']:28s} "
              f"{r['baseline_direction']:+d}  -> {r['AF_major_carrier_CT']:25s} "
              f"(dScore={r['AF_major_carrier_dScore']:+.4f}, q={r['AF_major_carrier_q_perm']:.3f})")

    print("\n=== per-(program × ct) significance tally ===")
    for ct in main_cts:
        sub = tbl[tbl["cell_type"] == ct]
        n_star = int((sub["p_perm"] < 0.05).sum())
        n_dot  = int(((sub["p_perm"] >= 0.05) & (sub["p_perm"] < 0.10)).sum())
        print(f"  {ct:25s}  * = {n_star}, . = {n_dot}, total = {len(sub)}")


def stopifnot(condition):
    if not condition:
        raise AssertionError()


if __name__ == "__main__":
    main()
