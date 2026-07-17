#!/usr/bin/env python3
"""Emit corrected-namespace Fig 4 manuscript TSV bundles."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
SCORES_TSV = RESULTS / "gse41177_family_sample_scores.tsv"
META_TSV = RESULTS / "gse41177_duration_metadata.tsv"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
LAA_TSV = RESULTS / "gse41177_family_duration_laa_primary.tsv"
LAPV_TSV = RESULTS / "gse41177_family_duration_lapv_sensitivity.tsv"
CHRONICITY_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"

OUT_FIG4B = RESULTS / "fig4b_family_duration_correlations.tsv"
OUT_FIG4C_FAM = RESULTS / "fig4c_chronicity_classes.tsv"
OUT_FIG4C_AGG = RESULTS / "fig4c_chronicity_aggregate.tsv"
OUT_DEF_SCATTER = RESULTS / "fig4def_exemplar_scatter.tsv"
OUT_DEF_SUM = RESULTS / "fig4def_exemplar_summary.tsv"
OUT_LAPV = RESULTS / "fig4_lapv_sensitivity.tsv"
OUT_LONG = RESULTS / "fig4_duration_scores_long.tsv"

FACETS = {
    "Fig 4D F9": 9,
    "Fig 4D F2": 2,
    "Fig 4D F15": 15,
}


def main() -> None:
    scores = pd.read_csv(SCORES_TSV, sep="\t")
    meta = pd.read_csv(META_TSV, sep="\t")
    labels = pd.read_csv(LABELS_TSV, sep="\t")
    laa = pd.read_csv(LAA_TSV, sep="\t").merge(labels, on="family_id", how="left")
    lapv = pd.read_csv(LAPV_TSV, sep="\t").merge(labels, on="family_id", how="left")
    chronicity = pd.read_csv(CHRONICITY_TSV, sep="\t").merge(labels, on="family_id", how="left")
    cls_map = chronicity.set_index("family_id")["chronicity_class"].to_dict()

    rows_b = []
    for region, df in [("LAA", laa), ("LAPV", lapv)]:
        for _, row in df.iterrows():
            fid = int(row["family_id"])
            rows_b.append({
                "family_id": fid,
                "family_label": row["family_label"],
                "axis_group": row.get("axis_group", ""),
                "region": region,
                "n_donors": int(row["n_samples"]),
                "rho": float(row["rho_spearman"]),
                "ci_low": float(row["bootstrap_ci_lo"]),
                "ci_high": float(row["bootstrap_ci_hi"]),
                "p_perm": float(row["p_perm"]),
                "p_spearman": float(row.get("p_spearman", np.nan)),
                "q_value": float(row.get("q_spearman_bh", np.nan)),
                "chronicity_class": cls_map.get(fid, ""),
                "regional_direction": (
                    "same" if bool(chronicity.loc[chronicity["family_id"] == fid, "regional_sign_agreement"].iloc[0])
                    else "opposite"
                ) if region == "LAA" else "",
            })
    pd.DataFrame(rows_b).to_csv(OUT_FIG4B, sep="\t", index=False, float_format="%.6g")

    c_rows = []
    for _, row in chronicity.iterrows():
        fid = int(row["family_id"])
        prow = laa.loc[laa["family_id"] == fid].iloc[0]
        c_rows.append({
            "family_id": fid,
            "family_label": row["family_label"],
            "rho": float(prow["rho_spearman"]),
            "p_perm": float(prow["p_perm"]),
            "ci_low": float(prow["bootstrap_ci_lo"]),
            "ci_high": float(prow["bootstrap_ci_hi"]),
            "chronicity_class": row["chronicity_class"],
            "direction": row["direction"],
            "p_spearman": float(prow.get("p_spearman", np.nan)),
            "q_value": float(prow.get("q_spearman_bh", np.nan)),
        })
    fig4c = pd.DataFrame(c_rows)
    fig4c.to_csv(OUT_FIG4C_FAM, sep="\t", index=False, float_format="%.6g")

    agg = (
        fig4c.groupby("chronicity_class", as_index=False)
        .agg(
            n_families=("family_id", "count"),
            family_ids=("family_id", lambda s: ",".join(str(int(x)) for x in sorted(s))),
            family_labels=("family_label", lambda s: "; ".join(s.astype(str))),
        )
    )
    agg.to_csv(OUT_FIG4C_AGG, sep="\t", index=False)

    lapv_rows = []
    for fid in sorted(laa["family_id"].unique()):
        pr = laa.loc[laa["family_id"] == fid].iloc[0]
        lr = lapv.loc[lapv["family_id"] == fid].iloc[0]
        rho_a, rho_v = float(pr["rho_spearman"]), float(lr["rho_spearman"])
        conc = np.sign(rho_a) == np.sign(rho_v) if np.isfinite(rho_v) else False
        lapv_rows.append({
            "family_id": int(fid),
            "family_label": pr["family_label"],
            "rho_LAA": rho_a,
            "p_LAA": float(pr["p_perm"]),
            "rho_LAPV": rho_v,
            "p_LAPV": float(lr["p_perm"]),
            "direction_concordant": conc,
            "sensitivity_class": "direction_concordant" if conc else "LAPV_opposite_sign",
        })
    pd.DataFrame(lapv_rows).to_csv(OUT_LAPV, sep="\t", index=False, float_format="%.6g")

    primary_ids = set(meta.loc[meta["primary_duration_eligible"], "sample_id"])
    long_df = scores.merge(
        meta[["sample_id", "donor_id", "rhythm", "duration_months", "region"]],
        on=["sample_id", "region"],
        how="left",
    ).merge(labels[["family_id", "family_label"]], on="family_id", how="left")
    long_df["score_type"] = "family_exact_reference_score"
    long_df.to_csv(
        OUT_LONG,
        sep="\t",
        index=False,
        columns=["donor_id", "sample_id", "region", "duration_months", "family_id",
                 "family_label", "score_type", "score"],
        float_format="%.6f",
    )

    scatter_rows = []
    summary_rows = []
    cls_idx = chronicity.set_index("family_id")
    for panel, fid in FACETS.items():
        sub = long_df[
            (long_df["family_id"] == fid)
            & (long_df["region"] == "LAA")
            & long_df["sample_id"].isin(primary_ids)
        ]
        for _, row in sub.iterrows():
            scatter_rows.append({
                "panel": panel,
                "sample_id": row["sample_id"],
                "donor_id": row["donor_id"],
                "region": "LAA",
                "rhythm": row["rhythm"],
                "duration_months": row["duration_months"],
                "family_id": fid,
                "family_label": row["family_label"],
                "family_score": float(row["score"]),
                "used_in_duration_correlation": True,
            })
        cr = cls_idx.loc[fid]
        summary_rows.append({
            "panel": panel,
            "family_id": fid,
            "family_label": cr["family_label"],
            "rho": float(cr["rho_primary"]),
            "ci_low": float(cr["ci_lo_primary"]),
            "ci_high": float(cr["ci_hi_primary"]),
            "p_perm": float(cr["p_perm_primary"]),
            "chronicity_class": str(cr["chronicity_class"]),
            "rho_lapv": float(cr["rho_lapv"]),
            "regional_direction": "same" if bool(cr["regional_sign_agreement"]) else "opposite",
        })

    pd.DataFrame(scatter_rows).to_csv(OUT_DEF_SCATTER, sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(summary_rows).to_csv(OUT_DEF_SUM, sep="\t", index=False, float_format="%.6g")

    for path in (OUT_FIG4B, OUT_FIG4C_FAM, OUT_FIG4C_AGG, OUT_DEF_SCATTER, OUT_DEF_SUM, OUT_LAPV, OUT_LONG):
        print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
