#!/usr/bin/env python3
"""Fig 3 v3 carrier-fidelity classification at family level."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

DELTA_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
OUT_TSV = ROOT / "results" / "fig3_v3_carrier_fidelity_class.tsv"
OUT_MANUSCRIPT_TSV = ROOT / "results" / "fig3c_family_fidelity_classes.tsv"

P_THRESH = 0.05


def evidence_strength_and_below_reason(cls: str, baseline: str, sub: pd.DataFrame) -> tuple[str, str]:
    """Return (evidence_strength, below_detection_reason). Reason empty unless Below_detection."""
    sig = sub[(sub["p_perm"].astype(float) < P_THRESH) & sub["ct_included"]]
    if cls == "Below_detection":
        parts: list[str] = []
        if baseline == "Mixed/Unresolved":
            parts.append("HCA_LA_carrier_prior_unresolved")
        inc = sub[sub["ct_included"]]
        if len(inc) == 0:
            parts.append("no_compartment_met_min_donors_per_rhythm_arm")
        else:
            bp = float(inc["p_perm"].astype(float).min())
            parts.append(f"best_p_perm={bp:.4g}")
            parts.append("no_cell_nominal_p_lt_0.05")
        return "none", "; ".join(parts)
    if sig.empty:
        return "nominal_sig", ""
    pmin = float(sig["p_perm"].astype(float).min())
    return f"p_perm_min={pmin:.4g}", ""


def classify(family_id: int, baseline_carrier: str, sub: pd.DataFrame) -> tuple[str, list[str]]:
    sig = sub[(sub["p_perm"].astype(float) < P_THRESH) & sub["ct_included"]]
    if sig.empty:
        return "Below_detection", []
    sig_cts = sig["compartment"].astype(str).tolist()
    carrier_sig_row = sig[sig["compartment"] == baseline_carrier]
    carrier_sig = not carrier_sig_row.empty
    other_sig = [c for c in sig_cts if c != baseline_carrier]

    if carrier_sig and other_sig:
        carrier_sign = float(carrier_sig_row["dScore"].iloc[0])
        same_sign = [
            c for c in other_sig
            if (sig.loc[sig["compartment"] == c, "dScore"].iloc[0] * carrier_sign) > 0
        ]
        if same_sign:
            return "Extension", sig_cts
        return "Mixed_direction", sig_cts
    if carrier_sig:
        return "Fidelity", sig_cts
    if other_sig:
        return "Relay", other_sig
    return "Below_detection", []


def main():
    if not LABELS_TSV.exists():
        raise FileNotFoundError(f"Missing {LABELS_TSV}")

    delta = pd.read_csv(DELTA_TSV, sep="\t")
    delta["family_id"] = delta["family_id"].astype(int)

    hca = pd.read_csv(HCA_TSV, sep="\t")
    hca["family_id"] = hca["family_id"].astype(int)
    carrier_la = dict(zip(hca["family_id"], hca["LA_carrier"].astype(str)))

    labs = pd.read_csv(LABELS_TSV, sep="\t")[["family_id", "family_label"]].drop_duplicates(
        subset=["family_id"], keep="first"
    )

    rows: list[dict] = []
    ms_rows: list[dict] = []
    for fid, sub in delta.groupby("family_id"):
        baseline = carrier_la.get(fid, "Mixed/Unresolved")
        cls, sig_cts = classify(int(fid), baseline, sub)
        ev, bd_reason = evidence_strength_and_below_reason(cls, baseline, sub)
        sig_joined = ",".join(sig_cts)
        rows.append({
            "family_id": int(fid),
            "direction": str(sub["direction"].iloc[0]),
            "baseline_LA_carrier": baseline,
            "n_sig_cts": int(len(sig_cts)),
            "sig_compartments": sig_joined,
            "fidelity_class": cls,
        })
        ms_rows.append({
            "family_id": int(fid),
            "HCA_carrier_prior": baseline,
            "snRNA_supported_compartment": sig_joined,
            "fidelity_class": cls,
            "evidence_strength": ev,
            "below_detection_reason": bd_reason if cls == "Below_detection" else "",
        })

    out = (
        pd.DataFrame(rows)
        .sort_values(["fidelity_class", "family_id"])
        .reset_index(drop=True)
    )
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print("wrote " + str(OUT_TSV.relative_to(ROOT)) + " (" + str(len(out)) + " rows)")
    print("\n=== class counts ===")
    print(out["fidelity_class"].value_counts().to_string())
    print("\n=== full table ===")
    print(out.to_string(index=False))

    ms = pd.DataFrame(ms_rows).merge(labs, on="family_id", how="left")
    if ms["family_label"].isna().any():
        miss = ms.loc[ms["family_label"].isna(), "family_id"].tolist()
        raise ValueError(f"Missing family_label for family_id (check {LABELS_TSV.name}): {miss}")
    ms = ms[
        [
            "family_id",
            "family_label",
            "HCA_carrier_prior",
            "snRNA_supported_compartment",
            "fidelity_class",
            "evidence_strength",
            "below_detection_reason",
        ]
    ].sort_values("family_id").reset_index(drop=True)
    ms.to_csv(OUT_MANUSCRIPT_TSV, sep="\t", index=False)
    print("\nwrote " + str(OUT_MANUSCRIPT_TSV.relative_to(ROOT)) + " (" + str(len(ms)) + " rows)")


if __name__ == "__main__":
    main()
