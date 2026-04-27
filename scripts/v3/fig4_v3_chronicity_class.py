#!/usr/bin/env python3
"""
Fig 4 v3 - 3-class chronicity classification (effect-size driven).

Strict thresholds (locked, no power-dependent gating):
    Strong-trend := |rho| >= 0.4   (in LAA primary)
    Weak-trend   := 0.2 <= |rho| < 0.4
    Stable       := |rho| < 0.2

Sensitivity concordance := count of {primary, gt1m_eq1, ordinal, LAPV}
                          codings whose rho has the same sign as primary.
                          Range 1-4 (primary itself counts).

Output
------
results/fig4_v3_chronicity_class.tsv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PRIMARY_TSV = ROOT / "results" / "fig4_v3_FAMILY_laa_primary.tsv"
GT1M_TSV = ROOT / "results" / "fig4_v3_FAMILY_laa_sens_gt1m_eq1.tsv"
ORD_TSV = ROOT / "results" / "fig4_v3_FAMILY_laa_sens_ordinal.tsv"
LAPV_TSV = ROOT / "results" / "fig4_v3_FAMILY_lapv_sensitivity.tsv"
OUT_TSV = ROOT / "results" / "fig4_v3_chronicity_class.tsv"


def classify(rho: float) -> str:
    if not np.isfinite(rho):
        return "Inconclusive"
    a = abs(rho)
    if a >= 0.4:
        return "Strong_trend"
    if a >= 0.2:
        return "Weak_trend"
    return "Stable"


def main() -> None:
    primary = pd.read_csv(PRIMARY_TSV, sep="\t").set_index("family_id")
    gt1m = pd.read_csv(GT1M_TSV, sep="\t").set_index("family_id")
    ordn = pd.read_csv(ORD_TSV, sep="\t").set_index("family_id")
    lapv = pd.read_csv(LAPV_TSV, sep="\t").set_index("family_id")

    rows = []
    for fid, prow in primary.iterrows():
        rho = float(prow["rho_spearman"])
        cls = classify(rho)
        rho_sign = np.sign(rho) if rho != 0 else 0

        sens_codings = []
        for label, df in [("primary", primary), ("gt1m_eq1", gt1m),
                          ("ordinal", ordn), ("lapv", lapv)]:
            if fid not in df.index:
                sens_codings.append((label, np.nan))
                continue
            sens_codings.append((label, float(df.loc[fid, "rho_spearman"])))

        same_sign_n = sum(
            1 for _, r in sens_codings
            if np.isfinite(r) and (np.sign(r) == rho_sign or rho_sign == 0)
        )
        sens_dict = {f"rho_{lbl}": r for lbl, r in sens_codings}

        rec = {
            "family_id": int(fid),
            "direction": str(prow["direction"]),
            "rho_primary": rho,
            "p_perm_primary": float(prow["p_perm"]),
            "ci_lo_primary": float(prow["bootstrap_ci_lo"]),
            "ci_hi_primary": float(prow["bootstrap_ci_hi"]),
            "n_samples_primary": int(prow["n_samples"]),
            "chronicity_class": cls,
            "sens_same_sign_n": int(same_sign_n),
            "sens_total_n": int(sum(1 for _, r in sens_codings if np.isfinite(r))),
            **sens_dict,
        }
        rows.append(rec)

    out = pd.DataFrame(rows)
    cls_order = {"Strong_trend": 0, "Weak_trend": 1, "Stable": 2, "Inconclusive": 3}
    out["cls_rank"] = out["chronicity_class"].map(cls_order)
    out = out.sort_values(["cls_rank", "family_id"]).drop(columns=["cls_rank"]).reset_index(drop=True)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.4f")
    print(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(out)} families)")
    print("\n=== class counts ===")
    print(out["chronicity_class"].value_counts().to_string())
    print("\n=== Strong-trend members + sensitivity concordance ===")
    for _, r in out[out["chronicity_class"] == "Strong_trend"].iterrows():
        print(f"  F{int(r['family_id']):2d}  rho={r['rho_primary']:+.2f}  "
              f"sens_same_sign={int(r['sens_same_sign_n'])}/{int(r['sens_total_n'])}  "
              f"(prim={r['rho_primary']:+.2f}, gt1m={r['rho_gt1m_eq1']:+.2f}, "
              f"ord={r['rho_ordinal']:+.2f}, lapv={r['rho_lapv']:+.2f})")


if __name__ == "__main__":
    main()
