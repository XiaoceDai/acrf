#!/usr/bin/env python3
"""
Fig 3 v3 carrier-fidelity classification at family level.

Strict thresholds (locked):
    sig             := family-level p_perm < 0.05  (BH q not enforced;
                       Fig 3 single-cell power is limited and the panel uses
                       nominal p as the visual indicator).
    Fidelity        := the family's HCA baseline LA-carrier compartment is
                       sig; AND no other compartment is sig.
    Relay           := baseline carrier is NOT sig; at least one OTHER
                       compartment is sig.
    Extension       := baseline carrier IS sig AND at least one other
                       compartment is sig with the SAME sign.
    Below detection := no compartment is sig.

Inputs
------
results/fig3_v3_family_celltype_delta.tsv
results/fig1f_v3_HCA_family_carriers.tsv

Outputs
-------
results/fig3_v3_carrier_fidelity_class.tsv
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

DELTA_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
OUT_TSV = ROOT / "results" / "fig3_v3_carrier_fidelity_class.tsv"

P_THRESH = 0.05


def classify(family_id, baseline_carrier, sub):
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
    delta = pd.read_csv(DELTA_TSV, sep="\t")
    delta["family_id"] = delta["family_id"].astype(int)

    hca = pd.read_csv(HCA_TSV, sep="\t")
    hca["family_id"] = hca["family_id"].astype(int)
    carrier_la = dict(zip(hca["family_id"], hca["LA_carrier"].astype(str)))

    rows = []
    for fid, sub in delta.groupby("family_id"):
        baseline = carrier_la.get(fid, "Mixed/Unresolved")
        cls, sig_cts = classify(fid, baseline, sub)
        rows.append({
            "family_id": int(fid),
            "direction": str(sub["direction"].iloc[0]),
            "baseline_LA_carrier": baseline,
            "n_sig_cts": int(len(sig_cts)),
            "sig_compartments": ",".join(sig_cts),
            "fidelity_class": cls,
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


if __name__ == "__main__":
    main()
