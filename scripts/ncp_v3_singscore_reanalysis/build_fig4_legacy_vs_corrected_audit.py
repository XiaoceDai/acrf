#!/usr/bin/env python3
"""Audit legacy Figure 4 duration outputs against corrected exact-reference pipeline."""
from __future__ import annotations

import importlib.util
import sys
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "ncp_v3_singscore_reanalysis"
LEGACY_CHRO = ROOT / "results" / "fig4_v3_chronicity_class.tsv"
CORR_CHRO = OUT_DIR / "fig4_v3_chronicity_class.tsv"
OUT_AUDIT = OUT_DIR / "fig4_legacy_vs_corrected_audit.tsv"
OUT_QC = OUT_DIR / "fig4_corrected_duration_qc.md"
V3 = ROOT / "scripts" / "v3"
KEY_FAMILIES = [2, 9, 10, 15]


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def recompute_triangulation(chro_path: Path) -> pd.DataFrame:
    fig5b = load_module("audit_fig5b", V3 / "fig5b_v3_triangulation_heatmap.py")
    fig5b.ATLAS_TSV = OUT_DIR / "fig_integrated_af_rewiring_atlas.tsv"
    fig5b.HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
    fig5b.REPL_TSV = OUT_DIR / "fig2e_cross_cohort_consistency.tsv"
    fig5b.DELTA_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
    fig5b.FID_TSV = ROOT / "results" / "fig3_v3_carrier_fidelity_class.tsv"
    fig5b.CHRO_TSV = chro_path
    fig5b.LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
    fig5b.OUT_TSV = OUT_DIR / "_tmp_fig5b_audit.tsv"
    fig5b.OUT_PDF = OUT_DIR / "_tmp_fig5b_audit.pdf"
    fig5b.main()
    return pd.read_csv(fig5b.OUT_TSV, sep="\t")


def family_audit_rows(legacy: pd.DataFrame, corrected: pd.DataFrame) -> list[dict]:
    rows = []
    merged = legacy.merge(
        corrected,
        on="family_id",
        how="outer",
        suffixes=("_legacy", "_corrected"),
        indicator=True,
    )
    for row in merged.itertuples():
        rho_laa_l = getattr(row, "rho_primary_legacy", np.nan)
        rho_laa_c = getattr(row, "rho_primary_corrected", np.nan)
        rho_lapv_l = getattr(row, "rho_lapv_legacy", np.nan)
        rho_lapv_c = getattr(row, "rho_lapv_corrected", np.nan)
        cls_l = getattr(row, "chronicity_class_legacy", "")
        cls_c = getattr(row, "chronicity_class_corrected", "")
        rows.append(
            {
                "section": "family",
                "family_id": int(row.family_id),
                "rho_LAA_legacy": rho_laa_l,
                "rho_LAA_corrected": rho_laa_c,
                "delta_rho_LAA": (
                    float(rho_laa_c - rho_laa_l)
                    if np.isfinite(rho_laa_l) and np.isfinite(rho_laa_c)
                    else np.nan
                ),
                "rho_LAPV_legacy": rho_lapv_l,
                "rho_LAPV_corrected": rho_lapv_c,
                "delta_rho_LAPV": (
                    float(rho_lapv_c - rho_lapv_l)
                    if np.isfinite(rho_lapv_l) and np.isfinite(rho_lapv_c)
                    else np.nan
                ),
                "chronicity_class_legacy": cls_l,
                "chronicity_class_corrected": cls_c,
                "class_changed": cls_l != cls_c,
            }
        )
    return rows


def summary_row(label: str, **values: object) -> dict:
    return {"section": "summary", "family_id": label, **values}


def gate_status(legacy: pd.DataFrame, corrected: pd.DataFrame, meta: pd.DataFrame) -> str:
    checks = [
        len(meta) == 38,
        int(meta["primary_duration_eligible"].sum()) == 12,
        len(corrected) == 16,
        corrected["rho_gt1m_eq1"].isna().all(),
        corrected["rho_ordinal"].isna().all(),
        (corrected["chronicity_class"] != "Inconclusive").all(),
    ]
    return "PASS" if all(checks) else "FAIL"


def main() -> None:
    for path in (LEGACY_CHRO, CORR_CHRO, OUT_DIR / "gse41177_duration_metadata.tsv"):
        assert path.is_file(), f"missing: {path}"

    legacy = pd.read_csv(LEGACY_CHRO, sep="\t")
    corrected = pd.read_csv(CORR_CHRO, sep="\t")
    meta = pd.read_csv(OUT_DIR / "gse41177_duration_metadata.tsv", sep="\t")

    rows = family_audit_rows(legacy, corrected)

    for cls in ["Strong_trend", "Weak_trend", "Stable"]:
        rows.append(
            summary_row(
                f"class_count_{cls}",
                chronicity_class_legacy=int((legacy["chronicity_class"] == cls).sum()),
                chronicity_class_corrected=int((corrected["chronicity_class"] == cls).sum()),
            )
        )

    for fid in KEY_FAMILIES:
        leg = legacy.loc[legacy["family_id"] == fid].iloc[0]
        cor = corrected.loc[corrected["family_id"] == fid].iloc[0]
        rows.append(
            summary_row(
                f"F{fid}",
                rho_LAA_legacy=float(leg["rho_primary"]),
                rho_LAA_corrected=float(cor["rho_primary"]),
                rho_LAPV_legacy=float(leg["rho_lapv"]),
                rho_LAPV_corrected=float(cor["rho_lapv"]),
                chronicity_class_legacy=str(leg["chronicity_class"]),
                chronicity_class_corrected=str(cor["chronicity_class"]),
                class_changed=leg["chronicity_class"] != cor["chronicity_class"],
            )
        )

    legacy_tri = pd.read_csv(ROOT / "results" / "fig5b_v3_triangulation.tsv", sep="\t")
    corrected_tri = recompute_triangulation(CORR_CHRO)
    path_a_legacy = legacy_tri.loc[
        legacy_tri["cross_layer_support_flag"], "family_id"
    ].tolist()
    path_a_corrected = corrected_tri.loc[
        corrected_tri["cross_layer_support_flag"], "family_id"
    ].tolist()
    rows.append(
        summary_row(
            "path_A_members",
            chronicity_class_legacy=",".join(f"F{int(x)}" for x in sorted(path_a_legacy)),
            chronicity_class_corrected=",".join(f"F{int(x)}" for x in sorted(path_a_corrected)),
            class_changed=path_a_legacy != path_a_corrected,
        )
    )

    for fid in KEY_FAMILIES:
        leg = legacy_tri.loc[legacy_tri["family_id"] == fid].iloc[0]
        cor = corrected_tri.loc[corrected_tri["family_id"] == fid].iloc[0]
        rows.append(
            summary_row(
                f"pathA_F{fid}",
                rho_LAA_legacy=int(leg["L3_supportive"]),
                rho_LAA_corrected=int(cor["L3_supportive"]),
                rho_LAPV_legacy=int(leg["n_supportive"]),
                rho_LAPV_corrected=int(cor["n_supportive"]),
                chronicity_class_legacy=str(leg["cross_layer_support_flag"]),
                chronicity_class_corrected=str(cor["cross_layer_support_flag"]),
                class_changed=bool(leg["cross_layer_support_flag"] != cor["cross_layer_support_flag"]),
            )
        )

    audit = pd.DataFrame(rows)
    audit.to_csv(OUT_AUDIT, sep="\t", index=False, float_format="%.4f")

    status = gate_status(legacy, corrected, meta)
    class_changes = audit[
        (audit["section"] == "family") & audit["class_changed"].fillna(False)
    ]
    f10 = corrected.loc[corrected["family_id"] == 10].iloc[0]
    f15 = corrected.loc[corrected["family_id"] == 15].iloc[0]

    qc_lines = [
        "# Figure 4 corrected duration QC",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        f"## Gate: **{status}**",
        "",
        "## Scorer proof",
        "- Scorer: `python_exact_reference` (Bioconductor singscore 1.28.1 formula)",
        "- Inputs recovered from GEO: `GSE41177_series_matrix.txt.gz`, `GSE41177_RAW.tar`",
        "- Expression: `expr_symbol_LAA.tsv`, `expr_symbol_LAPV.tsv` (38 GSM columns)",
        "- LE source: `results/fig3_v3_program_le.tsv` (30 programs)",
        "",
        "## Cohort counts",
        f"- Total samples: {len(meta)}",
        f"- Primary duration eligible (AF-LAA, finite months): {int(meta['primary_duration_eligible'].sum())}",
        f"- Censored AF-LAA (>1M): {int((meta['exclusion_reason'] == 'censored_gt1m').sum())}",
        f"- SR-LAA anchors: {int(((meta['rhythm'] == 'SR') & (meta['region'] == 'LAA')).sum())}",
        "",
        "## Class counts",
    ]
    for cls in ["Strong_trend", "Weak_trend", "Stable"]:
        n_leg = int((legacy["chronicity_class"] == cls).sum())
        n_cor = int((corrected["chronicity_class"] == cls).sum())
        qc_lines.append(f"- {cls}: legacy {n_leg} → corrected {n_cor}")

    qc_lines.extend(
        [
            "",
            "## Key families (F2 / F9 / F10 / F15)",
            "",
            "| Family | ρ LAA (leg→cor) | ρ LAPV (leg→cor) | Class (leg→cor) |",
            "|--------|-----------------|------------------|-----------------|",
        ]
    )
    for fid in KEY_FAMILIES:
        leg = legacy.loc[legacy["family_id"] == fid].iloc[0]
        cor = corrected.loc[corrected["family_id"] == fid].iloc[0]
        qc_lines.append(
            f"| F{fid} | {leg['rho_primary']:+.3f}→{cor['rho_primary']:+.3f} | "
            f"{leg['rho_lapv']:+.3f}→{cor['rho_lapv']:+.3f} | "
            f"{leg['chronicity_class']}→{cor['chronicity_class']} |"
        )

    qc_lines.extend(
        [
            "",
            "## Path A delta",
            f"- Legacy triple-positive: {', '.join(f'F{x}' for x in sorted(path_a_legacy)) or 'none'}",
            f"- Corrected triple-positive: {', '.join(f'F{x}' for x in sorted(path_a_corrected)) or 'none'}",
            f"- Path A membership changed: {path_a_legacy != path_a_corrected}",
            "",
            "## F10 / F15 triple-positive status (corrected)",
            f"- F10: L3_supportive={int(corrected_tri.loc[corrected_tri['family_id']==10,'L3_supportive'].iloc[0])}, "
            f"n_supportive={int(corrected_tri.loc[corrected_tri['family_id']==10,'n_supportive'].iloc[0])}, "
            f"class={f10['chronicity_class']}",
            f"- F15: L3_supportive={int(corrected_tri.loc[corrected_tri['family_id']==15,'L3_supportive'].iloc[0])}, "
            f"n_supportive={int(corrected_tri.loc[corrected_tri['family_id']==15,'n_supportive'].iloc[0])}, "
            f"class={f15['chronicity_class']}",
            "",
            "## Legacy column notes",
            "- `rho_gt1m_eq1` and `rho_ordinal` set to NA (sensitivity codings not recomputed).",
            "- `sens_same_sign_n` / `sens_total_n` redefined as primary+LAPV concordance (max 2).",
            "",
            f"## Families with class change: {len(class_changes)}",
        ]
    )
    if not class_changes.empty:
        for row in class_changes.itertuples():
            qc_lines.append(
                f"- F{int(row.family_id)}: {row.chronicity_class_legacy} → {row.chronicity_class_corrected}"
            )

    OUT_QC.write_text("\n".join(qc_lines) + "\n", encoding="utf-8")
    for tmp in (OUT_DIR / "_tmp_fig5b_audit.tsv", OUT_DIR / "_tmp_fig5b_audit.pdf"):
        if tmp.exists():
            tmp.unlink()

    print(f"wrote {OUT_AUDIT.relative_to(ROOT)}")
    print(f"wrote {OUT_QC.relative_to(ROOT)} (gate={status})")
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
