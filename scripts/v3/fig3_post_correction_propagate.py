#!/usr/bin/env python3
"""Post-correction audits, Path A impact, and fidelity refresh for Figure 3."""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"

FAM_NEW = RESULTS / "fig3_v3_family_celltype_delta.tsv"
PROG_NEW = RESULTS / "fig3_v3_program_celltype_delta.tsv"
FAM_OLD = RESULTS / "fig3_v3_family_celltype_delta_pre_correction.tsv"
PROG_OLD = RESULTS / "fig3_v3_program_celltype_delta_pre_correction.tsv"
FID_OLD = RESULTS / "fig3_v3_carrier_fidelity_class.tsv"
FID_OLD_SNAP = RESULTS / "fig3_v3_carrier_fidelity_class_pre_correction.tsv"
FID_SCRIPT = ROOT / "scripts" / "v3" / "fig3_v3_carrier_fidelity_class.py"
FID_NEW = RESULTS / "fig3_v3_carrier_fidelity_class.tsv"
FID_MS = RESULTS / "fig3c_family_fidelity_classes.tsv"
FIG5B = RESULTS / "fig5b_v3_triangulation.tsv"
FIG5B_OLD = RESULTS / "fig5b_v3_triangulation_pre_fig3_correction.tsv"

OUT_FAM_AUDIT = RESULTS / "fig3_family_compartment_old_vs_corrected.tsv"
OUT_PROG_AUDIT = RESULTS / "fig3_program_compartment_old_vs_corrected.tsv"
OUT_FID_AUDIT = RESULTS / "fig3_fidelity_old_vs_corrected.tsv"
OUT_PATHA = RESULTS / "fig3_figure5_path_a_impact.tsv"
OUT_BOUNDARY = RESULTS / "fig3_nominal_p_boundary_audit.tsv"
OUT_REPORT = RESULTS / "fig3_equivalence_gated_correction_report.md"

P_LO = 0.03
P_HI = 0.07


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    require(FAM_NEW.is_file() and PROG_NEW.is_file(), "Corrected delta tables missing")
    if FID_OLD.is_file() and not FID_OLD_SNAP.is_file():
        shutil.copy2(FID_OLD, FID_OLD_SNAP)

    # Refresh fidelity from corrected family table
    completed = subprocess.run(
        [sys.executable, str(FID_SCRIPT)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    require(completed.returncode == 0, f"Fidelity classifier failed:\n{completed.stdout}\n{completed.stderr}")
    print(completed.stdout.strip())

    fam_new = pd.read_csv(FAM_NEW, sep="\t")
    prog_new = pd.read_csv(PROG_NEW, sep="\t")
    fid_new = pd.read_csv(FID_NEW, sep="\t")
    fid_ms = pd.read_csv(FID_MS, sep="\t")

    # Family old vs corrected
    if FAM_OLD.is_file():
        fam_old = pd.read_csv(FAM_OLD, sep="\t")
        merged = fam_old.merge(
            fam_new,
            on=["family_id", "compartment"],
            suffixes=("_old", "_new"),
            how="outer",
        )
        merged["dScore_delta"] = merged["dScore_new"] - merged["dScore_old"]
        merged["p_perm_delta"] = merged["p_perm_new"] - merged["p_perm_old"]
        merged["sig_old"] = merged["p_perm_old"] < 0.05
        merged["sig_new"] = merged["p_perm_new"] < 0.05
        merged["sig_changed"] = merged["sig_old"] != merged["sig_new"]
        merged.to_csv(OUT_FAM_AUDIT, sep="\t", index=False, float_format="%.10g")
    else:
        merged = fam_new.copy()
        merged.to_csv(OUT_FAM_AUDIT, sep="\t", index=False, float_format="%.10g")

    if PROG_OLD.is_file():
        prog_old = pd.read_csv(PROG_OLD, sep="\t")
        pm = prog_old.merge(
            prog_new,
            on=["program_id", "compartment"],
            suffixes=("_old", "_new"),
            how="outer",
        )
        pm["dScore_delta"] = pm["dScore_new"] - pm["dScore_old"]
        pm["p_perm_delta"] = pm["p_perm_new"] - pm["p_perm_old"]
        pm["sig_changed"] = (pm["p_perm_old"] < 0.05) != (pm["p_perm_new"] < 0.05)
        pm.to_csv(OUT_PROG_AUDIT, sep="\t", index=False, float_format="%.10g")

    # Fidelity audit
    if FID_OLD_SNAP.is_file():
        fid_old = pd.read_csv(FID_OLD_SNAP, sep="\t")
        fa = fid_old.merge(fid_new, on="family_id", suffixes=("_old", "_new"), how="outer")
        fa["class_changed"] = fa["fidelity_class_old"] != fa["fidelity_class_new"]
        fa.to_csv(OUT_FID_AUDIT, sep="\t", index=False)
        class_changed_ids = fa.loc[fa["class_changed"], "family_id"].astype(int).tolist()
    else:
        fa = fid_new.copy()
        fa["class_changed"] = False
        fa.to_csv(OUT_FID_AUDIT, sep="\t", index=False)
        class_changed_ids = []

    f15_new = fid_ms.loc[fid_ms["family_id"] == 15].iloc[0]
    f16_new = fid_ms.loc[fid_ms["family_id"] == 16].iloc[0]
    f15_old_class = f16_old_class = None
    if FID_OLD_SNAP.is_file():
        fo = pd.read_csv(FID_OLD_SNAP, sep="\t").set_index("family_id")
        f15_old_class = str(fo.loc[15, "fidelity_class"]) if 15 in fo.index else None
        f16_old_class = str(fo.loc[16, "fidelity_class"]) if 16 in fo.index else None

    f15_changed = f15_old_class is not None and f15_old_class != str(f15_new["fidelity_class"])
    f16_changed = f16_old_class is not None and f16_old_class != str(f16_new["fidelity_class"])

    # Path A impact from current Fig 5B vs refreshed L2 definition
    if FIG5B.is_file() and not FIG5B_OLD.is_file():
        shutil.copy2(FIG5B, FIG5B_OLD)

    path_rows = []
    for fid in sorted(fid_new["family_id"].astype(int).unique()):
        cls = str(fid_new.loc[fid_new["family_id"] == fid, "fidelity_class"].iloc[0])
        l2_new = int(cls != "Below_detection")
        l2_old = None
        path_a_old = None
        if FIG5B_OLD.is_file():
            old5 = pd.read_csv(FIG5B_OLD, sep="\t")
            if fid in set(old5["family_id"].astype(int)):
                row = old5.loc[old5["family_id"].astype(int) == fid].iloc[0]
                l2_old = int(row["L2_supportive"])
                path_a_old = bool(row["cross_layer_support_flag"])
        # Recompute Path A eligibility for L2 only; full Path A needs L1+L3 refresh
        path_rows.append(
            {
                "family_id": fid,
                "fidelity_class_new": cls,
                "L2_supportive_new": l2_new,
                "L2_supportive_old": l2_old,
                "L2_changed": (l2_old is not None) and (l2_old != l2_new),
                "cross_layer_support_flag_old": path_a_old,
                "needs_fig5_regen": fid in (15, 16) and (
                    (fid == 15 and f15_changed) or (fid == 16 and f16_changed)
                ),
            }
        )
    path_df = pd.DataFrame(path_rows)
    path_df.to_csv(OUT_PATHA, sep="\t", index=False)

    # Boundary p values
    boundary = fam_new.loc[
        fam_new["ct_included"]
        & fam_new["p_perm"].between(P_LO, P_HI, inclusive="both")
    ][
        [
            "family_id",
            "compartment",
            "dScore",
            "p_perm",
            "p_perm_mc_se",
            "sig_marker",
            "n_af",
            "n_ctrl",
        ]
    ].sort_values("p_perm")
    boundary.to_csv(OUT_BOUNDARY, sep="\t", index=False, float_format="%.10g")

    class_counts = fid_new["fidelity_class"].value_counts().to_dict()

    # Regenerate Fig 5 if F15/F16 class changed
    regenerated_fig5 = False
    if f15_changed or f16_changed:
        regen = ROOT / "scripts" / "ncp_v3_singscore_reanalysis" / "regenerate_fig5.py"
        if regen.is_file():
            completed = subprocess.run(
                [sys.executable, str(regen)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            require(
                completed.returncode == 0,
                f"Fig 5 regeneration failed:\n{completed.stdout}\n{completed.stderr}",
            )
            regenerated_fig5 = True
            print(completed.stdout.strip())
        else:
            # Fall back to fig5b script directly
            fig5b_script = ROOT / "scripts" / "v3" / "fig5b_v3_triangulation_heatmap.py"
            completed = subprocess.run(
                [sys.executable, str(fig5b_script)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            require(completed.returncode == 0, f"fig5b failed:\n{completed.stderr}")
            regenerated_fig5 = True

    # Path A after optional regen
    path_a_families = []
    if FIG5B.is_file():
        fig5 = pd.read_csv(FIG5B, sep="\t")
        path_a_families = (
            fig5.loc[fig5["cross_layer_support_flag"].astype(bool), "family_id"]
            .astype(int)
            .tolist()
        )

    eq = pd.read_csv(RESULTS / "fig3_scorer_equivalence_summary.tsv", sep="\t").iloc[0]

    report = "\n".join(
        [
            "# Figure 3 equivalence-gated correction report",
            "",
            "## Scorer equivalence",
            f"- Exact Python ↔ R singscore 1.28.1: **{'PASS' if bool(eq['exact_vs_r_passed']) else 'FAIL'}**",
            f"- max |Δ| exact vs R: `{float(eq['max_abs_diff_exact_vs_r']):.3g}`",
            f"- Legacy v3 ↔ canonical: **{'PASS' if bool(eq['legacy_v3_equivalent_to_canonical']) else 'FAIL'}**",
            f"- max |Δ| legacy vs exact: `{float(eq['max_abs_diff_legacy_vs_exact']):.3g}`",
            f"- Production scorer: `{eq['production_scorer']}`",
            f"- R used for production scores: `{bool(eq['r_package_used_for_production_scores'])}`",
            "",
            "## Final class counts",
            *[f"- {k}: {int(v)}" for k, v in sorted(class_counts.items())],
            "",
            "## F15 / F16",
            f"- F15 class: `{f15_new['fidelity_class']}` (was `{f15_old_class}`); "
            f"supported: `{f15_new['snRNA_supported_compartment']}`",
            f"- F16 class: `{f16_new['fidelity_class']}` (was `{f16_old_class}`); "
            f"supported: `{f16_new['snRNA_supported_compartment']}`",
            f"- F15 class changed: `{f15_changed}`",
            f"- F16 class changed: `{f16_changed}`",
            "",
            "## Path A impact",
            f"- Families with class change: `{class_changed_ids}`",
            f"- Figure 5 regenerated: `{regenerated_fig5}`",
            f"- Path A families (cross_layer_support_flag): `{path_a_families}`",
            "",
            "## Nominal p near 0.05 boundary "
            f"([{P_LO}, {P_HI}])",
            f"- n boundary family×compartment cells: `{len(boundary)}`",
            "",
        ]
    )
    if len(boundary):
        report += boundary.to_string(index=False) + "\n"
    OUT_REPORT.write_text(report, encoding="utf-8")
    print(report)
    print(f"Wrote {OUT_REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
