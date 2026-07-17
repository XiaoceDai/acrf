#!/usr/bin/env python3
"""High-precision permutation sensitivity for three Figure 3 boundary cells.

Prespecified primary results remain the 10,000-MC (or exact) p_perm values in
fig3_v3_family_celltype_delta.tsv. This script does NOT overwrite them.

Targets:
  F11 × adipocyte
  F15 × endothelial
  F16 × adipocyte

Uses the same per-test rng_seed from fig3_v3_permutation_test_order.tsv, then
continues Monte Carlo to N_SENS = 1_000_000 for a sensitivity p and MC SE.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SCORES = ROOT / "results" / "fig3_v3_singscore_per_donor.tsv"
FAM_DELTA = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
TEST_ORDER = ROOT / "results" / "fig3_v3_permutation_test_order.tsv"
OUT_TSV = ROOT / "results" / "fig3_boundary_permutation_sensitivity.tsv"
OUT_QC = ROOT / "results" / "fig3_boundary_permutation_sensitivity_qc.md"

N_PRIMARY = 10_000
N_SENS = 1_000_000
TARGETS = [
    (11, "adipocyte"),
    (15, "endothelial"),
    (16, "adipocyte"),
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def family_donor_scores(scores: pd.DataFrame, family_id: int, compartment: str) -> pd.DataFrame:
    g = scores.loc[
        (scores["family_id"].astype(int) == family_id)
        & (scores["compartment"] == compartment)
        & scores["score"].notna()
    ]
    fam = (
        g.groupby(["donor_id", "rhythm"], as_index=False)["score"]
        .mean()
        .sort_values(["rhythm", "donor_id"])
        .reset_index(drop=True)
    )
    return fam


def mc_p(values: np.ndarray, is_af: np.ndarray, rng: np.random.Generator, n_perm: int) -> tuple[float, float, int]:
    n_af = int(is_af.sum())
    obs = float(values[is_af].mean() - values[~is_af].mean())
    hits = 1
    for _ in range(n_perm):
        perm = rng.permutation(len(values))
        mask = np.zeros(len(values), dtype=bool)
        mask[perm[:n_af]] = True
        stat = float(values[mask].mean() - values[~mask].mean())
        if abs(stat) >= abs(obs) - 1e-12:
            hits += 1
    p = hits / (n_perm + 1)
    se = float(math.sqrt(max(p * (1.0 - p), 0.0) / n_perm))
    return p, se, hits


def main() -> None:
    scores = pd.read_csv(SCORES, sep="\t")
    primary = pd.read_csv(FAM_DELTA, sep="\t")
    order = pd.read_csv(TEST_ORDER, sep="\t")
    order = order.loc[order["test_level"] == "family"].copy()

    rows = []
    for family_id, compartment in TARGETS:
        test_id = f"family|{family_id}|{compartment}"
        tor = order.loc[order["test_id"] == test_id]
        require(len(tor) == 1, f"Missing test order row for {test_id}")
        seed = int(tor.iloc[0]["rng_seed"])

        prim = primary.loc[
            (primary["family_id"].astype(int) == family_id)
            & (primary["compartment"] == compartment)
        ]
        require(len(prim) == 1, f"Missing primary row for {test_id}")
        p_row = prim.iloc[0]
        p_primary = float(p_row["p_perm"])
        se_primary = float(p_row["p_perm_mc_se"]) if pd.notna(p_row["p_perm_mc_se"]) else float("nan")
        n_af_tab = int(p_row["n_af"])
        n_ctrl_tab = int(p_row["n_ctrl"])

        fam = family_donor_scores(scores, family_id, compartment)
        ctrl = fam.loc[fam["rhythm"] == "Ctrl", "score"].to_numpy(dtype=float)
        af = fam.loc[fam["rhythm"] == "AF", "score"].to_numpy(dtype=float)
        require(len(ctrl) == n_ctrl_tab and len(af) == n_af_tab, f"Donor counts drifted for {test_id}")

        values = np.concatenate([ctrl, af])
        is_af = np.concatenate([np.zeros(len(ctrl), dtype=bool), np.ones(len(af), dtype=bool)])
        rng = np.random.default_rng(seed)
        p_sens, se_sens, hits = mc_p(values, is_af, rng, N_SENS)

        rows.append(
            {
                "family_id": family_id,
                "compartment": compartment,
                "test_id": test_id,
                "rng_seed": seed,
                "n_AF": n_af_tab,
                "n_ctrl": n_ctrl_tab,
                "dScore": float(p_row["dScore"]),
                "p_perm_primary_10k": p_primary,
                "p_perm_mc_se_primary": se_primary,
                "n_perm_primary": int(p_row["n_perm"]),
                "perm_strategy_primary": str(p_row["perm_strategy"]),
                "p_perm_sensitivity_1e6": p_sens,
                "p_perm_mc_se_sensitivity": se_sens,
                "n_perm_sensitivity": N_SENS,
                "sensitivity_hits_plus_one": hits,
                "crosses_0_05_primary": bool(p_primary < 0.05),
                "crosses_0_05_sensitivity": bool(p_sens < 0.05),
                "publication_primary": "p_perm_primary_10k",
                "note": "sensitivity does not replace primary p_perm",
            }
        )
        print(
            f"{test_id}: primary={p_primary:.6g} sens_1e6={p_sens:.6g} "
            f"(seed={seed})"
        )

    out = pd.DataFrame(rows)
    out.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.10g")
    OUT_QC.write_text(
        "\n".join(
            [
                "# Figure 3 boundary permutation sensitivity QC",
                "",
                f"- Targets: {TARGETS}",
                f"- Primary N: {N_PRIMARY} (from locked family delta table)",
                f"- Sensitivity N: {N_SENS}",
                "- Same per-test rng_seed as primary schedule; sensitivity continues MC only.",
                "- **Primary publication p-values are unchanged.**",
                "",
                out.to_string(index=False),
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")
    print(f"wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
