#!/usr/bin/env python3
"""
P4 — donor-level aCAI_rel classification into Figure 2 rewiring classes.

Reads the unified long singscore table produced by P3, the bulk reference
from P3, and produces:

  results/af_chamber_rewiring_programs.tsv    per-program × per-cohort stats
  results/af_chamber_rewiring_consensus.tsv   consensus subset + 115574 rep
  results/figure2_source_data.tsv             long, one row per (cohort, donor, program)

All statistics are computed on aCAI_rel (fraction-of-healthy), per
FIGURE2_WORKFLOW.md §3.

Permutation strategy (spec §3.3):
  * ≤ 20 000 combinations -> exact enumeration
  * else                  -> Monte Carlo, 10 000 draws, seed 20260418
"""
from __future__ import annotations

import math
import sys
from itertools import combinations
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

# Locked thresholds (see FIGURE2_WORKFLOW.md §3.4/§3.5).
EPS = 0.30          # baseline support threshold on aCAI_rel
DELTA = 0.30        # shift gate on |dCAI_rel|
AMP_RATIO = 1.25    # aCAI_rel_AF / aCAI_rel_SR lower bound for Amplified
ATT_RATIO = 0.75    # upper bound for Attenuated
PERM_EXACT_LIMIT = 20000
N_MC_PERM = 10000
MC_SEED = 20260418
PVAL_THRESHOLD = 0.10   # perm gate

# Expected dCAI_rel sign per rewiring class (for 115574 replication check).
CLASS_EXPECTED_SIGN = {
    "Amplified":       +1,
    "Amplified-weak":  +1,
    "Preserved":        0,   # around zero
    "Attenuated":      -1,
    "Attenuated-weak": -1,
    "Reversed":        -1,
    "Reversed-weak":   -1,
}

# §3.6 consensus lookup table (discovery = 128188 × 79768).
# Key is (class_128188, class_79768).
CONSENSUS_TABLE = {
    ("Amplified",    "Amplified"):    "Amplified",
    ("Preserved",    "Preserved"):    "Preserved",
    ("Attenuated",   "Attenuated"):   "Attenuated",
    ("Reversed",     "Reversed"):     "Reversed",
    ("Amplified",    "Preserved"):    "Amplified-weak",
    ("Preserved",    "Amplified"):    "Amplified-weak",
    ("Attenuated",   "Preserved"):    "Attenuated-weak",
    ("Preserved",    "Attenuated"):   "Attenuated-weak",
    ("Reversed",     "Preserved"):    "Reversed-weak",
    ("Preserved",    "Reversed"):     "Reversed-weak",
    ("Reversed",     "Attenuated"):   "Reversed-weak",
    ("Attenuated",   "Reversed"):     "Reversed-weak",
    # explicit NA pairs (kept for audit; any missing key also defaults to NA)
    ("Amplified",    "Attenuated"):   pd.NA,
    ("Attenuated",   "Amplified"):    pd.NA,
    ("Amplified",    "Reversed"):     pd.NA,
    ("Reversed",     "Amplified"):    pd.NA,
}


# ----------------------------------------------------------------------
# permutation on rhythm labels (donor-level, balanced by rhythm)
# ----------------------------------------------------------------------

def _stat(values: np.ndarray, is_af: np.ndarray) -> float:
    """mean(AF) - mean(SR)."""
    return float(values[is_af].mean() - values[~is_af].mean())


def permutation_p(values: np.ndarray, is_af: np.ndarray,
                  rng: np.random.Generator) -> Tuple[float, str, int]:
    """Two-sided permutation p-value for mean(AF) - mean(SR).

    Labels are shuffled while keeping n_AF fixed.  Returns (p, strategy, n_perm)
    where strategy ∈ {"exact", "mc"} and n_perm is the number of
    permutations actually evaluated (including the observed statistic's
    canonical count for exact enumeration).
    """
    n = len(values)
    n_af = int(is_af.sum())
    n_sr = n - n_af
    if n_af == 0 or n_sr == 0:
        return (float("nan"), "none", 0)

    observed = _stat(values, is_af)
    n_combinations = math.comb(n, n_af)

    if n_combinations <= PERM_EXACT_LIMIT:
        # exact enumeration over all ways to pick the AF subset
        count_as_extreme = 0
        for af_idx in combinations(range(n), n_af):
            mask = np.zeros(n, dtype=bool)
            mask[list(af_idx)] = True
            if abs(_stat(values, mask)) >= abs(observed) - 1e-12:
                count_as_extreme += 1
        p = count_as_extreme / n_combinations
        return (p, "exact", n_combinations)

    # Monte Carlo
    count_as_extreme = 1  # include observed
    for _ in range(N_MC_PERM):
        perm = rng.permutation(n)
        mask = np.zeros(n, dtype=bool)
        mask[perm[:n_af]] = True
        if abs(_stat(values, mask)) >= abs(observed) - 1e-12:
            count_as_extreme += 1
    p = count_as_extreme / (N_MC_PERM + 1)
    return (p, "mc", N_MC_PERM)


# ----------------------------------------------------------------------
# classification
# ----------------------------------------------------------------------

def _apply_effect_gates(mean_rel_sr: float, mean_rel_af: float,
                        dcai_rel: float) -> dict:
    """Compute the three effect/shift gate booleans used by §3.5 and §3.5b."""
    reversed_effect = (mean_rel_af < -EPS)
    amp_effect = (mean_rel_af > AMP_RATIO * mean_rel_sr)
    amp_shift  = (dcai_rel > DELTA)
    att_effect = (mean_rel_af < ATT_RATIO * mean_rel_sr) and (mean_rel_af >= -EPS)
    att_shift  = (-dcai_rel > DELTA)
    return dict(
        reversed_effect=reversed_effect,
        amp_effect=amp_effect, amp_shift=amp_shift,
        att_effect=att_effect, att_shift=att_shift,
    )


def _label_from_gates(g: dict, perm_pass: bool) -> str:
    """Apply class precedence Reversed > Amplified > Attenuated > Preserved.

    `perm_pass=True` implements spec §3.5 (strict). `perm_pass=False` gives
    Preserved for everything baseline-supported; §3.5b uses a separate call
    with perm_pass implicitly True (the permutation gate is simply ignored).
    """
    if not perm_pass:
        return "Preserved"
    if g["reversed_effect"]:
        return "Reversed"
    if g["amp_effect"] and g["amp_shift"]:
        return "Amplified"
    if g["att_effect"] and g["att_shift"]:
        return "Attenuated"
    return "Preserved"


def classify(mean_rel_sr: float, mean_rel_af: float, dcai_rel: float,
             p_perm: float, baseline_supported: bool) -> dict:
    """Return strict §3.5 classification + §3.5b trend_class + gate flags."""
    if not baseline_supported:
        return dict(
            effect_gate_pass = False,
            shift_gate_pass  = False,
            perm_gate_pass   = False,
            rewiring_class   = "Low-support",
            trend_class      = "Low-support",
        )

    gates = _apply_effect_gates(mean_rel_sr, mean_rel_af, dcai_rel)
    perm_pass = (p_perm is not None) and (not math.isnan(p_perm)) and (p_perm < PVAL_THRESHOLD)

    # Strict class (perm gate enforced)
    strict = _label_from_gates(gates, perm_pass=perm_pass)
    # Trend class (perm gate ignored)
    trend = _label_from_gates(gates, perm_pass=True)

    return dict(
        effect_gate_pass = any([gates["reversed_effect"], gates["amp_effect"], gates["att_effect"]]),
        shift_gate_pass  = abs(dcai_rel) > DELTA,
        perm_gate_pass   = perm_pass,
        rewiring_class   = strict,
        trend_class      = trend,
    )


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main() -> None:
    long_path = ROOT / "Project_P1" / "program_scores_long.tsv"
    bulk_path = ROOT / "Project_P1" / "bulk_reference_aCAI.tsv"
    out_dir = ROOT / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    long = pd.read_csv(long_path, sep="\t")
    bulk = pd.read_csv(bulk_path, sep="\t")[[
        "program_id", "mean_aCAI_bulk", "abs_aCAI_bulk",
        "n_LE_total", "n_LE_in_bulk", "le_coverage_bulk",
        "n_donors_bulk",
    ]]
    bulk = bulk.rename(columns={"mean_aCAI_bulk": "aCAI_bulk"})

    # ------------------------------------------------------------------
    # figure2_source_data.tsv — long format, one row per (cohort, donor, program)
    # ------------------------------------------------------------------
    src = long.merge(bulk[["program_id", "aCAI_bulk"]], on="program_id", how="left")
    # reorder columns to match spec §6.3
    col_order = [
        "cohort", "donor_id", "rhythm",
        "program_id", "program_label", "program_order", "baseline_direction",
        "gsm_LA", "gsm_RA",
        "score_LA", "score_RA",
        "CAI", "aCAI",
        "aCAI_bulk", "abs_aCAI_bulk", "aCAI_rel",
    ]
    src = src[col_order]
    src_out = out_dir / "figure2_source_data.tsv"
    src.to_csv(src_out, sep="\t", index=False)
    print(f"wrote {src_out.relative_to(ROOT)}  ({len(src)} rows)")

    # ------------------------------------------------------------------
    # per-program × per-cohort summary
    # ------------------------------------------------------------------
    rng = np.random.default_rng(MC_SEED)
    cohorts = ["GSE128188", "GSE79768", "GSE115574"]
    programs = long[["program_id", "program_label", "program_order", "baseline_direction"]].drop_duplicates()
    programs = programs.sort_values("program_order").reset_index(drop=True)

    rows = []
    for _, prog in programs.iterrows():
        row = dict(
            program_id         = prog["program_id"],
            program_label      = prog["program_label"],
            program_order      = int(prog["program_order"]),
            baseline_direction = int(prog["baseline_direction"]),
        )
        # join bulk reference
        b = bulk[bulk["program_id"] == prog["program_id"]].iloc[0]
        row["aCAI_bulk"]      = float(b["aCAI_bulk"])
        row["abs_aCAI_bulk"]  = float(b["abs_aCAI_bulk"])
        row["n_LE_total"]     = int(b["n_LE_total"])
        row["n_LE_in_bulk"]   = int(b["n_LE_in_bulk"])
        row["le_coverage_bulk"] = float(b["le_coverage_bulk"])
        row["n_donors_bulk"]  = int(b["n_donors_bulk"])

        for c in cohorts:
            sub = long[(long["cohort"] == c) & (long["program_id"] == prog["program_id"])]
            rel_sr = sub.loc[sub["rhythm"] == "SR", "aCAI_rel"].to_numpy()
            rel_af = sub.loc[sub["rhythm"] == "AF", "aCAI_rel"].to_numpy()
            raw_sr = sub.loc[sub["rhythm"] == "SR", "aCAI"].to_numpy()
            raw_af = sub.loc[sub["rhythm"] == "AF", "aCAI"].to_numpy()
            n_sr, n_af = len(rel_sr), len(rel_af)

            row[f"n_SR_{c}"]            = n_sr
            row[f"n_AF_{c}"]            = n_af
            row[f"mean_aCAI_SR_{c}"]    = float(raw_sr.mean()) if n_sr else float("nan")
            row[f"mean_aCAI_AF_{c}"]    = float(raw_af.mean()) if n_af else float("nan")
            m_rel_sr = float(rel_sr.mean()) if n_sr else float("nan")
            m_rel_af = float(rel_af.mean()) if n_af else float("nan")
            row[f"mean_aCAI_rel_SR_{c}"] = m_rel_sr
            row[f"mean_aCAI_rel_AF_{c}"] = m_rel_af
            # sign check is on raw aCAI sign (equivalent to aCAI_rel sign since abs_aCAI_bulk>0)
            n_pos = int((raw_sr > 0).sum())
            row[f"frac_SR_positive_{c}"] = n_pos / n_sr if n_sr else float("nan")
            row[f"n_SR_positive_{c}"]    = n_pos
            n_req = int(math.ceil(0.8 * n_sr)) if n_sr else 0
            baseline_supported = (m_rel_sr > EPS) and (n_pos >= n_req) if n_sr else False
            row[f"sr_baseline_supported_{c}"] = bool(baseline_supported)

            dcai_rel = m_rel_af - m_rel_sr
            row[f"dCAI_rel_{c}"] = dcai_rel
            row[f"ratio_AF_over_SR_{c}"] = (
                (m_rel_af / m_rel_sr) if baseline_supported and m_rel_sr != 0 else float("nan")
            )

            # permutation: pool SR/AF aCAI_rel values with matched label vector
            values = np.concatenate([rel_sr, rel_af])
            is_af  = np.concatenate([np.zeros(n_sr, dtype=bool), np.ones(n_af, dtype=bool)])
            p_perm, strategy, n_perm = permutation_p(values, is_af, rng)
            row[f"p_perm_{c}"]   = p_perm
            row[f"n_perm_{c}"]   = n_perm
            row[f"perm_strategy_{c}"] = strategy

            cls = classify(m_rel_sr, m_rel_af, dcai_rel, p_perm, baseline_supported)
            for k, v in cls.items():
                row[f"{k}_{c}"] = v

        rows.append(row)

    tbl = pd.DataFrame(rows)

    # BH adjustment per cohort across the 10 programs
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

    for c in cohorts:
        tbl[f"q_perm_{c}"] = bh(tbl[f"p_perm_{c}"].to_numpy())

    # ------------------------------------------------------------------
    # consensus (128188 × 79768) and 115574 replication
    # ------------------------------------------------------------------
    def lookup_consensus(c1: str, c2: str) -> object:
        if c1 == "Low-support" or c2 == "Low-support":
            return pd.NA
        return CONSENSUS_TABLE.get((c1, c2), pd.NA)

    tbl["consensus_class"] = [
        lookup_consensus(a, b)
        for a, b in zip(tbl["rewiring_class_GSE128188"], tbl["rewiring_class_GSE79768"])
    ]

    # 115574 replication: sign of dCAI_rel_115574 against consensus class expectation
    def rep_flag(row) -> dict:
        cc = row["consensus_class"]
        out = dict(
            replicated_in_115574 = False,
            replication_sign_concordance_115574 = False,
            replication_dCAI_rel_115574 = row["dCAI_rel_GSE115574"],
        )
        if pd.isna(cc) or not row["sr_baseline_supported_GSE115574"]:
            return out
        exp_sign = CLASS_EXPECTED_SIGN.get(cc, None)
        if exp_sign is None:
            return out
        d = row["dCAI_rel_GSE115574"]
        if math.isnan(d):
            return out
        if exp_sign == 0:  # Preserved: expect small |d|
            out["replication_sign_concordance_115574"] = abs(d) <= DELTA
        else:
            out["replication_sign_concordance_115574"] = (np.sign(d) == exp_sign)
        out["replicated_in_115574"] = (
            out["replication_sign_concordance_115574"]
            and (row["p_perm_GSE115574"] < PVAL_THRESHOLD)
        )
        return out

    rep_df = pd.DataFrame([rep_flag(r) for _, r in tbl.iterrows()])
    tbl = pd.concat([tbl, rep_df], axis=1)

    # Column ordering: program meta first, then per-cohort blocks, then consensus/replication
    cohort_block = lambda c: [
        f"n_SR_{c}", f"n_AF_{c}",
        f"mean_aCAI_SR_{c}", f"mean_aCAI_AF_{c}",
        f"mean_aCAI_rel_SR_{c}", f"mean_aCAI_rel_AF_{c}",
        f"n_SR_positive_{c}", f"frac_SR_positive_{c}", f"sr_baseline_supported_{c}",
        f"dCAI_rel_{c}", f"ratio_AF_over_SR_{c}",
        f"perm_strategy_{c}", f"n_perm_{c}",
        f"p_perm_{c}", f"q_perm_{c}",
        f"effect_gate_pass_{c}", f"shift_gate_pass_{c}", f"perm_gate_pass_{c}",
        f"rewiring_class_{c}", f"trend_class_{c}",
    ]
    meta_cols = [
        "program_id", "program_label", "program_order", "baseline_direction",
        "aCAI_bulk", "abs_aCAI_bulk", "n_LE_total", "n_LE_in_bulk",
        "le_coverage_bulk", "n_donors_bulk",
    ]
    consensus_cols = [
        "consensus_class",
        "replication_sign_concordance_115574",
        "replication_dCAI_rel_115574",
        "replicated_in_115574",
    ]
    ordered = meta_cols.copy()
    for c in cohorts:
        ordered.extend(cohort_block(c))
    ordered.extend(consensus_cols)
    tbl = tbl[ordered]

    progs_out = out_dir / "af_chamber_rewiring_programs.tsv"
    tbl.to_csv(progs_out, sep="\t", index=False)
    print(f"wrote {progs_out.relative_to(ROOT)}  ({len(tbl)} rows, {tbl.shape[1]} cols)")

    # consensus-only subset
    tbl_consensus = tbl[tbl["consensus_class"].notna()].copy()
    cons_out = out_dir / "af_chamber_rewiring_consensus.tsv"
    tbl_consensus.to_csv(cons_out, sep="\t", index=False)
    print(f"wrote {cons_out.relative_to(ROOT)}  ({len(tbl_consensus)} rows)")

    # ---- stdout summary ----
    print("\n=== per-cohort rewiring class tally (strict §3.5) ===")
    for c in cohorts:
        col = f"rewiring_class_{c}"
        print(f"  {c:10s}  " + ", ".join(
            f"{k}={v}" for k, v in tbl[col].value_counts().sort_index().items()
        ))
    print("\n=== per-cohort trend_class tally (exploratory §3.5b, perm gate ignored) ===")
    for c in cohorts:
        col = f"trend_class_{c}"
        print(f"  {c:10s}  " + ", ".join(
            f"{k}={v}" for k, v in tbl[col].value_counts().sort_index().items()
        ))
    print(f"\n  consensus: " + ", ".join(
        f"{k}={v}" for k, v in tbl["consensus_class"].fillna("NA").value_counts().sort_index().items()
    ))
    rep_pos = int(tbl["replicated_in_115574"].sum())
    print(f"  replicated in 115574: {rep_pos} / {tbl['consensus_class'].notna().sum()} consensus programs")


if __name__ == "__main__":
    main()
