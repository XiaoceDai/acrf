#!/usr/bin/env python3
"""
Build the consolidated v3 anchor evidence table from:

  - results/figure6_mr_anchors.tsv               (4 v2 anchors, atrial)
  - results/figure6_mr_anchors_blood.tsv         (4 v2 anchors, blood)
  - results/figure6_mr_anchors_v3.tsv            (6 v3 anchors, atrial + blood)
  - results/figure6_mr_positive_controls.tsv     (PITX2 calibrator,
                                                  atrial + blood)

Applies the FIGURE6_WORKFLOW.md §10 verdict logic:

  §10.5  verdict (4-class)            — primary tissue
  §10.6  shape_tag (B/C/pass/weak/...) — primary tissue, decoupled from gate
  §10.7  family_verdict                — multi-anchor families only
  §10.8  figure-level scenario         — top-down 5-class

Output:
  results/figure6_mr_anchors_v3_consolidated.tsv
  results/figure6_mr_anchors_v3_summary.tsv     (one-row scenario summary)
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES  = ROOT / "results"

V2_ATRIAL = RES / "figure6_mr_anchors.tsv"
V2_BLOOD  = RES / "figure6_mr_anchors_blood.tsv"
V3        = RES / "figure6_mr_anchors_v3.tsv"
PC        = RES / "figure6_mr_positive_controls.tsv"
EGENE     = RES / "figure6_egene_status.tsv"
SQTL      = RES / "figure6_mr_anchors_v3_sqtl.tsv"      # h4f output
SGENE     = RES / "figure6_sgene_status.tsv"             # h4g output

OUT       = RES / "figure6_mr_anchors_v3_consolidated.tsv"
OUT_SUM   = RES / "figure6_mr_anchors_v3_summary.tsv"

# ---------------------------------------------------------------------------
# Anchor metadata (from FIGURE6_WORKFLOW.md §10.2, locked)
# ---------------------------------------------------------------------------
ANCHOR_META = pd.DataFrame([
    # v2 pre-registered (4)
    ("PRELP",  "F10",       "matrix",     "v2_pre_registered", "anchor"),
    ("OMD",    "F10",       "matrix",     "v2_pre_registered", "anchor"),
    ("IDH2",   "Metabolic", "metabolic",  "v2_pre_registered", "anchor"),
    ("NDUFA8", "F8",        "metabolic",  "v2_pre_registered", "anchor"),
    # v3 LE extension (6)
    ("NTRK2",  "F15",       "receptor",   "v3_LE_extension",   "anchor"),
    ("NTF3",   "F15",       "ligand",     "v3_LE_extension",   "anchor"),
    ("FRS2",   "F15",       "adaptor",    "v3_LE_extension",   "anchor"),
    ("BMP10",  "F11",       "ligand",     "v3_LE_extension",   "anchor"),
    ("SMAD9",  "F11",       "R-SMAD",     "v3_LE_extension",   "anchor"),
    ("SMAD6",  "F11",       "i-SMAD",     "v3_LE_extension",   "anchor"),
    # PITX2 calibrator
    ("PITX2",  "calibrator","calibrator", "v2_positive_control", "calibrator"),
], columns=["anchor", "family", "role", "source", "kind"])

# Anchors-only ordering on the 6B matrix (calibrator on its own row at bottom)
ANCHOR_ORDER = list(ANCHOR_META.query("kind=='anchor'")["anchor"])
CALIBRATOR_ORDER = list(ANCHOR_META.query("kind=='calibrator'")["anchor"])

# Multi-anchor families (only these can earn a family_verdict per §10.7)
MULTI_ANCHOR_FAMS = {
    "F10": 2,
    "F11": 3,
    "F15": 3,
}

# ---------------------------------------------------------------------------
# §10.5 verdict (4-class), §10.6 shape_tag — primary tissue only
# ---------------------------------------------------------------------------
# Plan C eligibility (spec §10.3 v3 follow-up, 2026-04-23):
#   PRIMARY     : GTEx v8 atrial appendage eGene (q < 0.05, BH-FDR on p_beta)
#   SENSITIVITY : per-gene Bonferroni (lead p < 0.05 / n_variants_in_window)
#   LEGACY (off): lead p < 5e-8  — methodologically misapplied, kept only
#                 as a transparency column.
# ---------------------------------------------------------------------------
def compute_verdict(row) -> str:
    """Verdict driven by PRIMARY eligibility (eGene status)."""
    if not bool(row.get("eligibility_primary", False)):
        return "no_eligible_eQTL"
    pp4 = row.get("coloc_PP4")
    if pd.isna(pp4):
        return "no_eligible_eQTL"
    if pp4 >= 0.7:
        return "supported"
    if pp4 >= 0.5:
        return "weak_hybrid"
    return "not_supported"


def compute_shape_tag(row) -> str:
    """Independent of eligibility gate (§10.6). Operates on PP3/PP4 pattern only.
    PITX2 itself is allowed to be no_eligible_eQTL while still earning tag B.
    """
    pp0 = row.get("coloc_PP0"); pp1 = row.get("coloc_PP1")
    pp3 = row.get("coloc_PP3"); pp4 = row.get("coloc_PP4")
    if pd.isna(pp4):
        return "no_data"
    if not pd.isna(pp4) and pp4 >= 0.7:
        return "pass"
    if not pd.isna(pp4) and pp4 >= 0.5:
        return "weak"
    if not pd.isna(pp3) and pp3 >= 0.5 and pp4 < 0.3:
        return "B"  # PITX2-mode (high PP3, low PP4)
    if (not pd.isna(pp0)) and (not pd.isna(pp1)) and (pp0 + pp1) >= 0.7:
        return "C"  # eQTL absence
    return "fail-other"


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def load_primary() -> pd.DataFrame:
    """Build one-row-per-anchor primary-tissue table (atrial appendage)."""
    # v2 atrial
    v2a = pd.read_csv(V2_ATRIAL, sep="\t")
    v2a["tissue"] = "atrial_appendage"

    # v3 atrial
    v3 = pd.read_csv(V3, sep="\t")
    v3a = v3.query("tissue == 'atrial_appendage'").copy()

    # PITX2 atrial
    pc = pd.read_csv(PC, sep="\t")
    pc_a = pc.query("anchor == 'PITX2' and tissue == 'atrial_appendage'").copy()

    # Harmonise columns: pick the union we need
    keep_cols = ["anchor", "tissue", "n_eqtl", "lead_eqtl_p", "gate_passed",
                 "n_gwas", "n_overlap", "gwas_min_p",
                 "coloc_PP0", "coloc_PP1", "coloc_PP2",
                 "coloc_PP3", "coloc_PP4"]
    # v2 atrial has no `prefilter_passed`; impute from lead_eqtl_p
    for df in (v2a, v3a, pc_a):
        if "prefilter_passed" not in df.columns:
            df["prefilter_passed"] = df["lead_eqtl_p"].apply(
                lambda p: (not pd.isna(p)) and (p < 1e-4))
    keep_cols.append("prefilter_passed")
    v2a = v2a[keep_cols].copy()
    v3a = v3a[keep_cols].copy()
    pc_a = pc_a[keep_cols].copy()

    out = pd.concat([v2a, v3a, pc_a], ignore_index=True)
    return out


def load_blood() -> pd.DataFrame:
    """Build one-row-per-anchor blood-sensitivity table (whole blood).

    Returns columns:
      anchor, lead_eqtl_p_blood, gate_passed_blood, prefilter_passed_blood,
      coloc_PP4_blood, coloc_PP3_blood
    """
    rows = []
    # v2 blood
    if V2_BLOOD.exists():
        v2b = pd.read_csv(V2_BLOOD, sep="\t")
        for _, r in v2b.iterrows():
            rows.append({
                "anchor": r["anchor"],
                "lead_eqtl_p_blood": r.get("lead_eqtl_p"),
                "gate_passed_blood": r.get("gate_passed"),
                "prefilter_passed_blood": (
                    (not pd.isna(r.get("lead_eqtl_p"))) and r.get("lead_eqtl_p") < 1e-4
                ),
                "coloc_PP4_blood": r.get("coloc_PP4"),
                "coloc_PP3_blood": pd.NA,  # v2 blood doesn't store PP3
            })

    # v3 blood
    if V3.exists():
        v3 = pd.read_csv(V3, sep="\t")
        v3b = v3.query("tissue == 'whole_blood'")
        for _, r in v3b.iterrows():
            rows.append({
                "anchor": r["anchor"],
                "lead_eqtl_p_blood": r.get("lead_eqtl_p"),
                "gate_passed_blood": r.get("gate_passed"),
                "prefilter_passed_blood": r.get("prefilter_passed"),
                "coloc_PP4_blood": r.get("coloc_PP4"),
                "coloc_PP3_blood": r.get("coloc_PP3"),
            })

    # PITX2 blood
    if PC.exists():
        pc = pd.read_csv(PC, sep="\t")
        pc_b = pc.query("anchor == 'PITX2' and tissue == 'whole_blood'")
        for _, r in pc_b.iterrows():
            rows.append({
                "anchor": r["anchor"],
                "lead_eqtl_p_blood": r.get("lead_eqtl_p"),
                "gate_passed_blood": r.get("gate_passed"),
                "prefilter_passed_blood": (
                    (not pd.isna(r.get("lead_eqtl_p"))) and r.get("lead_eqtl_p") < 1e-4
                ),
                "coloc_PP4_blood": r.get("coloc_PP4"),
                "coloc_PP3_blood": r.get("coloc_PP3"),
            })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Family + figure-level scenario
# ---------------------------------------------------------------------------
def compute_family_verdicts(consolidated: pd.DataFrame) -> dict[str, str]:
    """For each multi-anchor family, apply ⌈2n/3⌉ rule to count of supported anchors.
    Uses verdict_combined (eQTL OR sQTL) so an anchor counts as 'supported' if
    EITHER molecular layer reaches PP4>=0.7 (spec §10.14 RNA-layer aggregate).
    Single-anchor families and the PITX2 calibrator return NA_underpowered / NA.
    """
    out: dict[str, str] = {}
    anchors = consolidated.query("kind == 'anchor'")
    v_col = ("verdict_combined" if "verdict_combined" in anchors.columns
             else "verdict")
    for fam, grp in anchors.groupby("family"):
        n = len(grp)
        n_pass = int((grp[v_col] == "supported").sum())
        if fam not in MULTI_ANCHOR_FAMS:
            out[fam] = "NA_underpowered"
            continue
        threshold = math.ceil(2 * n / 3)
        out[fam] = "pass" if n_pass >= threshold else "fail"
    return out


def compute_scenario(consolidated: pd.DataFrame, fam_verdicts: dict) -> str:
    """§10.8 + §10.14 figure-level scenario, top-down first-match.
    PITX2 (calibrator) does not contribute to anchor counts.
    Uses verdict_combined which aggregates eQTL and sQTL layers.
    """
    anchors = consolidated.query("kind == 'anchor'")
    # Use combined RNA-layer verdict if present, else fall back to eqtl-only
    v_col = ("verdict_combined" if "verdict_combined" in anchors.columns
             else "verdict")
    n_pass = int((anchors[v_col] == "supported").sum())
    n_weak = int((anchors[v_col] == "weak_hybrid").sum())
    n_fam_pass = sum(1 for v in fam_verdicts.values() if v == "pass")
    if n_fam_pass >= 2:
        return "G  multi-family hybrid"
    if n_fam_pass == 1:
        return "F  family-level hybrid"
    if n_pass >= 1:
        return "E  single anchor hybrid"
    if n_weak >= 1:
        return "D  weak hybrid (suggestive)"
    return "A  pure distributed"


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
def load_egene_status() -> dict[str, dict]:
    """Load Plan C eligibility data from h4e_egene_status.R output.
    Returns dict[anchor] -> {atrial_q, atrial_eGene, blood_q, blood_eGene}."""
    if not EGENE.exists():
        raise FileNotFoundError(
            f"Missing {EGENE}; run scripts/h4e_egene_status.R first")
    eg = pd.read_csv(EGENE, sep="\t")
    out: dict[str, dict] = {}
    for anchor, sub in eg.groupby("anchor"):
        atr = sub.query("tissue == 'atrial_appendage'")
        bld = sub.query("tissue == 'whole_blood'")
        out[anchor] = {
            "atrial_q":     float(atr["q_BH"].iloc[0]) if len(atr)
                            and not pd.isna(atr["q_BH"].iloc[0]) else float("nan"),
            "atrial_eGene": bool(atr["is_eGene_q05"].iloc[0]) if len(atr)
                            else False,
            "blood_q":      float(bld["q_BH"].iloc[0]) if len(bld)
                            and not pd.isna(bld["q_BH"].iloc[0]) else float("nan"),
            "blood_eGene":  bool(bld["is_eGene_q05"].iloc[0]) if len(bld)
                            else False,
        }
    return out


def load_sqtl_layer() -> pd.DataFrame:
    """Load atrial sQTL coloc results (h4f) + sGene FDR status (h4g)
    and return one row per anchor with the sQTL-layer fields needed by
    the consolidated table.  If either file is missing, return an empty
    skeleton so the rest of the pipeline still runs (sQTL columns NaN)."""
    cols = ["anchor", "sqtl_n_clusters_tested", "sqtl_n_clusters_eligible",
            "sqtl_lead_junction", "sqtl_lead_eqtl_p",
            "sqtl_PP4_max", "sqtl_PP3_at_max",
            "is_sGene_atrial", "sGene_q_min"]
    if not SQTL.exists():
        return pd.DataFrame(columns=cols)
    sq = pd.read_csv(SQTL, sep="\t")
    sq = sq.rename(columns={
        "n_clusters_tested":   "sqtl_n_clusters_tested",
        "n_clusters_eligible": "sqtl_n_clusters_eligible",
        "lead_junction_id":    "sqtl_lead_junction",
        "lead_eqtl_p":         "sqtl_lead_eqtl_p",
        "PP4_max":             "sqtl_PP4_max",
        "PP3_at_max":          "sqtl_PP3_at_max",
    })[["anchor",
        "sqtl_n_clusters_tested", "sqtl_n_clusters_eligible",
        "sqtl_lead_junction", "sqtl_lead_eqtl_p",
        "sqtl_PP4_max", "sqtl_PP3_at_max"]]
    if SGENE.exists():
        sg = pd.read_csv(SGENE, sep="\t")[["anchor", "is_sGene_q05",
                                            "q_BH_min"]]
        sg = sg.rename(columns={"is_sGene_q05": "is_sGene_atrial",
                                "q_BH_min": "sGene_q_min"})
        sq = sq.merge(sg, on="anchor", how="left")
    else:
        sq["is_sGene_atrial"] = False
        sq["sGene_q_min"]     = float("nan")
    return sq


def main() -> None:
    primary = load_primary()
    blood   = load_blood()
    eg_map  = load_egene_status()
    sqtl    = load_sqtl_layer()

    # ---- Plan C eligibility (spec §10.3 v3 follow-up) ----
    # PRIMARY = GTEx atrial appendage eGene status (q < 0.05).
    # SENSITIVITY = per-gene Bonferroni (lead_eqtl_p < 0.05/n_eqtl).
    # Legacy 5e-8 retained as `gate_5e8_legacy` for transparency only.
    primary["eGene_atrial_q"]     = primary["anchor"].map(
        lambda a: eg_map.get(a, {}).get("atrial_q", float("nan")))
    primary["is_eGene_atrial"]    = primary["anchor"].map(
        lambda a: bool(eg_map.get(a, {}).get("atrial_eGene", False)))
    primary["eGene_blood_q"]      = primary["anchor"].map(
        lambda a: eg_map.get(a, {}).get("blood_q", float("nan")))
    primary["is_eGene_blood"]     = primary["anchor"].map(
        lambda a: bool(eg_map.get(a, {}).get("blood_eGene", False)))

    # Per-gene Bonferroni
    primary["bonf_threshold"] = primary.apply(
        lambda r: (0.05 / float(r["n_eqtl"]))
                  if (not pd.isna(r["n_eqtl"]) and r["n_eqtl"] > 0)
                  else float("nan"),
        axis=1,
    )
    primary["bonf_passed"] = primary.apply(
        lambda r: (not pd.isna(r["lead_eqtl_p"])
                   and not pd.isna(r["bonf_threshold"])
                   and float(r["lead_eqtl_p"]) < float(r["bonf_threshold"])),
        axis=1,
    )

    # Legacy 5e-8 (NEVER used as gate; transparency only)
    primary["gate_5e8_legacy"] = primary["lead_eqtl_p"].apply(
        lambda p: (not pd.isna(p)) and float(p) < 5e-8)

    # Plan C primary / sensitivity / agreement
    primary["eligibility_primary"]     = primary["is_eGene_atrial"]
    primary["eligibility_sensitivity"] = primary["bonf_passed"]
    primary["eligibility_agreement"]   = (
        primary["eligibility_primary"] == primary["eligibility_sensitivity"])

    # Verdict + shape tag are computed on PRIMARY tissue only
    primary["verdict"]   = primary.apply(compute_verdict, axis=1)
    primary["shape_tag"] = primary.apply(compute_shape_tag, axis=1)

    # ---- sQTL layer join (spec §10.14, sQTL RNA-layer extension) ----
    primary = primary.merge(sqtl, on="anchor", how="left")
    primary["is_sGene_atrial"] = primary["is_sGene_atrial"].fillna(False)
    # Per-RNA-layer verdicts.
    # eQTL verdict uses eligibility_primary (eGene); sQTL verdict uses
    # sGene status as the parallel eligibility gate.
    def verdict_eqtl(row) -> str:
        return row["verdict"]      # already computed above

    def verdict_sqtl(row) -> str:
        if not bool(row.get("is_sGene_atrial", False)):
            return "no_eligible_sQTL"
        pp4 = row.get("sqtl_PP4_max")
        if pd.isna(pp4):
            return "no_eligible_sQTL"
        pp4 = float(pp4)
        if pp4 >= 0.7:  return "supported"
        if pp4 >= 0.5:  return "weak_hybrid"
        return "not_supported"

    primary["verdict_eqtl"] = primary.apply(verdict_eqtl, axis=1)
    primary["verdict_sqtl"] = primary.apply(verdict_sqtl, axis=1)

    def verdict_combined(row) -> str:
        """Combined RNA-layer verdict: any layer 'supported' wins; otherwise
        the weakest 'no_eligible' if both layers are no_eligible; else the
        strictly stronger of the two."""
        rank = {"supported": 4, "weak_hybrid": 3, "not_supported": 2,
                "no_eligible_eQTL": 1, "no_eligible_sQTL": 1,
                "no_eligible": 1}
        ve = row["verdict_eqtl"]
        vs = row["verdict_sqtl"]
        # If both no_eligible, return generic no_eligible
        if rank.get(ve, 0) == 1 and rank.get(vs, 0) == 1:
            return "no_eligible"
        # Else pick the stronger
        return ve if rank.get(ve, 0) >= rank.get(vs, 0) else vs

    primary["verdict_combined"] = primary.apply(verdict_combined, axis=1)

    # Merge anchor metadata + blood sensitivity
    df = ANCHOR_META.merge(primary, on="anchor", how="left")
    df = df.merge(blood, on="anchor", how="left")

    # §10.9 PITX2 calibrator contingency:
    # PITX2 was pre-declared to reach shape_tag = B.  If it does not, every
    # anchor-level B tag is demoted to fail-other and the discrepancy is
    # surfaced in `calibrator_status`.
    pitx2_row = df.query("anchor == 'PITX2'")
    pitx2_tag = pitx2_row["shape_tag"].iloc[0] if len(pitx2_row) else "NA"
    calibrator_validated = bool(pitx2_tag == "B")
    if not calibrator_validated:
        n_demoted = int(((df["kind"] == "anchor")
                         & (df["shape_tag"] == "B")).sum())
        df.loc[(df["kind"] == "anchor") & (df["shape_tag"] == "B"),
               "shape_tag"] = "fail-other_B-demoted"
        calibrator_status = (
            f"FAIL: PITX2 reached shape_tag={pitx2_tag} "
            f"(PP3={float(pitx2_row['coloc_PP3'].iloc[0]):.3f}, "
            f"threshold>=0.50); "
            f"{n_demoted} anchor B tag(s) demoted to fail-other "
            f"per spec \u00a710.9"
        )
    else:
        calibrator_status = "PASS: PITX2 reproduces shape_tag=B as expected"
    df["calibrator_status"] = calibrator_status

    # Family verdict + scenario
    fam_verdicts = compute_family_verdicts(df)
    df["family_verdict"] = df["family"].map(fam_verdicts).fillna("NA")
    scenario = compute_scenario(df, fam_verdicts)

    # Order: anchors in canonical order, calibrator last
    order_map = {a: i for i, a in enumerate(ANCHOR_ORDER + CALIBRATOR_ORDER)}
    df["__ord"] = df["anchor"].map(order_map)
    df = df.sort_values("__ord").drop(columns="__ord").reset_index(drop=True)

    # Reorder columns to a stable schema
    col_order = [
        "anchor", "family", "role", "source", "kind",
        "tissue",  # always atrial_appendage in this primary table
        "n_eqtl", "lead_eqtl_p",
        # Plan C eligibility columns (PRIMARY = eGene; SENSITIVITY = Bonferroni)
        "eGene_atrial_q", "is_eGene_atrial",
        "bonf_threshold", "bonf_passed",
        "eligibility_primary", "eligibility_sensitivity",
        "eligibility_agreement",
        "gate_5e8_legacy",        # transparency only — never used as gate
        "gate_passed", "prefilter_passed",  # legacy v3 columns retained
        "n_gwas", "n_overlap", "gwas_min_p",
        "coloc_PP0", "coloc_PP1", "coloc_PP2",
        "coloc_PP3", "coloc_PP4",
        # eQTL-layer verdict (single layer)
        "verdict", "shape_tag",
        # sQTL layer (RNA-layer extension, spec §10.14)
        "is_sGene_atrial", "sGene_q_min",
        "sqtl_n_clusters_tested", "sqtl_n_clusters_eligible",
        "sqtl_lead_junction", "sqtl_lead_eqtl_p",
        "sqtl_PP4_max", "sqtl_PP3_at_max",
        "verdict_eqtl", "verdict_sqtl", "verdict_combined",
        # family + calibrator
        "family_verdict", "calibrator_status",
        # blood sensitivity (right-side columns)
        "eGene_blood_q", "is_eGene_blood",
        "lead_eqtl_p_blood", "gate_passed_blood", "prefilter_passed_blood",
        "coloc_PP4_blood", "coloc_PP3_blood",
    ]
    for c in col_order:
        if c not in df.columns:
            df[c] = pd.NA
    df = df[col_order]

    # Write consolidated long table
    df.to_csv(OUT, sep="\t", index=False)

    # Write one-row scenario summary
    anchors = df.query("kind == 'anchor'")
    # Plan C summary across the 10 anchors only (excluding calibrator)
    plan_c_primary = int(anchors["eligibility_primary"].sum())
    plan_c_sens    = int(anchors["eligibility_sensitivity"].sum())
    plan_c_agree   = int(anchors["eligibility_agreement"].sum())
    plan_c_disagree_anchors = ", ".join(
        sorted(anchors.query("not eligibility_agreement")["anchor"].tolist())) or "(none)"

    # Per-RNA-layer counts (spec §10.14)
    n_supp_eqtl = int((anchors["verdict_eqtl"] == "supported").sum())
    n_supp_sqtl = int((anchors["verdict_sqtl"] == "supported").sum())
    n_supp_either = int((anchors["verdict_combined"] == "supported").sum())
    n_sGene = int(anchors["is_sGene_atrial"].sum())

    summary = pd.DataFrame([{
        "scenario": scenario,
        "n_anchors_total": len(anchors),
        # eQTL-layer (existing fields, unchanged for back-compat)
        "n_supported": int((anchors["verdict"] == "supported").sum()),
        "n_weak_hybrid": int((anchors["verdict"] == "weak_hybrid").sum()),
        "n_not_supported": int((anchors["verdict"] == "not_supported").sum()),
        "n_no_eligible_eQTL": int((anchors["verdict"] == "no_eligible_eQTL").sum()),
        # sQTL layer
        "n_sGene_atrial":           n_sGene,
        "n_supp_sqtl":              n_supp_sqtl,
        "n_supp_eqtl":              n_supp_eqtl,
        "n_supp_either_RNA_layer":  n_supp_either,
        # Shape tags + family
        "n_shape_B_atrial":  int((anchors["shape_tag"] == "B").sum()),
        "n_shape_C_atrial":  int((anchors["shape_tag"] == "C").sum()),
        "n_fam_pass": sum(1 for v in fam_verdicts.values() if v == "pass"),
        "fam_verdicts": "; ".join(f"{k}={v}" for k, v in fam_verdicts.items()),
        # Plan C agreement audit
        "planC_n_primary_eligible":     plan_c_primary,
        "planC_n_sensitivity_eligible": plan_c_sens,
        "planC_n_agreement":            plan_c_agree,
        "planC_disagreement_anchors":   plan_c_disagree_anchors,
        "PITX2_shape_tag": (df.query("anchor == 'PITX2'")["shape_tag"].iloc[0]
                            if "PITX2" in df["anchor"].values else "NA"),
        "PITX2_calibrator_validated": calibrator_validated,
        "calibrator_status": calibrator_status,
    }])
    summary.to_csv(OUT_SUM, sep="\t", index=False)

    # Console echo
    print(f"=== Consolidated v3 anchor verdict — {scenario} ===")
    cols_show = ["anchor", "family", "kind",
                 "is_eGene_atrial", "is_sGene_atrial",
                 "coloc_PP4", "sqtl_PP4_max",
                 "verdict_eqtl", "verdict_sqtl", "verdict_combined"]
    cols_show = [c for c in cols_show if c in df.columns]
    print(df[cols_show].to_string(index=False))
    print()
    print(f"Family verdicts: {fam_verdicts}")
    print(f"Figure-level scenario: {scenario}")
    print(f"Plan C eligibility: PRIMARY (eGene)={plan_c_primary}/{len(anchors)}, "
          f"SENSITIVITY (Bonf)={plan_c_sens}/{len(anchors)}, "
          f"AGREEMENT={plan_c_agree}/{len(anchors)} "
          f"[disagreements: {plan_c_disagree_anchors}]")
    print(f"sQTL RNA-layer (\u00a710.14): n_sGene={n_sGene}/{len(anchors)}; "
          f"supp_eqtl={n_supp_eqtl}; supp_sqtl={n_supp_sqtl}; "
          f"supp_either={n_supp_either}")
    print()
    print(f"Wrote {OUT.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUM.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
