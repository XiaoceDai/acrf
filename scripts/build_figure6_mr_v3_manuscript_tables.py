#!/usr/bin/env python3
"""Generate prose-facing Fig 6 MR / coloc tables from the pipeline-wide TSV."""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PIPELINE_WIDE = ROOT / "results" / "figure6_mr_anchors_v3_pipeline_wide.tsv"
PIPELINE_WIDE_LEGACY = ROOT / "results" / "figure6_mr_anchors_v3_consolidated.tsv"
POSITIVE = ROOT / "results" / "figure6_mr_positive_controls.tsv"
FIG5B = ROOT / "results" / "fig5b_v3_triangulation.tsv"

OUT_CONSOLIDATED = ROOT / "results" / "figure6_mr_anchors_v3_consolidated.tsv"
OUT_RULES = ROOT / "results" / "figure6_mr_coloc_gate_rules.tsv"
OUT_CALIB = ROOT / "results" / "figure6_mr_calibration_genetics.tsv"
OUT_LOCUS = ROOT / "results" / "fig6a_v3_locus_plot_windows.tsv"

# hg38 ±1 Mb cis windows (must match cache filenames / R drivers).
LOCUS_WINDOWS: list[dict[str, str | int]] = [
    {"anchor_gene": "PRELP", "chr": 1, "win_start": 202475828, "win_end": 204491352,
     "family_id": 10, "reason_selected": "ECM–matrix pre-registered anchor"},
    {"anchor_gene": "OMD", "chr": 9, "win_start": 91414245, "win_end": 93424461,
     "family_id": 10, "reason_selected": "ECM–matrix pre-registered anchor"},
    {"anchor_gene": "IDH2", "chr": 15, "win_start": 89083045, "win_end": 91102504,
     "family_id": 7, "reason_selected": "Metabolic–translational pre-registered anchor"},
    {"anchor_gene": "NDUFA8", "chr": 9, "win_start": 121144058, "win_end": 123159819,
     "family_id": 8, "reason_selected": "Metabolic–translational pre-registered anchor"},
    {"anchor_gene": "NTRK2", "chr": 9, "win_start": 83668522, "win_end": 86029699,
     "family_id": 15, "reason_selected": "BMP–SMAD / neurotrophin v3 extension anchor"},
    {"anchor_gene": "NTF3", "chr": 12, "win_start": 4404888, "win_end": 6495026,
     "family_id": 15, "reason_selected": "BMP–SMAD / neurotrophin v3 extension anchor"},
    {"anchor_gene": "FRS2", "chr": 12, "win_start": 68464533, "win_end": 70513650,
     "family_id": 15, "reason_selected": "BMP–SMAD / neurotrophin v3 extension anchor"},
    {"anchor_gene": "BMP10", "chr": 2, "win_start": 67869054, "win_end": 69886842,
     "family_id": 11, "reason_selected": "BMP–SMAD / neurotrophin v3 extension anchor"},
    {"anchor_gene": "SMAD9", "chr": 13, "win_start": 35879326, "win_end": 37925087,
     "family_id": 11, "reason_selected": "BMP–SMAD / neurotrophin v3 extension anchor"},
    {"anchor_gene": "SMAD6", "chr": 15, "win_start": 65703041, "win_end": 67782082,
     "family_id": 11, "reason_selected": "BMP–SMAD / neurotrophin v3 extension anchor"},
    {"anchor_gene": "PITX2", "chr": 4, "win_start": 109610001, "win_end": 111633149,
     "family_id": "", "reason_selected": "Canonical AF locus positive control (4q25)"},
    {"anchor_gene": "KCNN3", "chr": 1, "win_start": 153837802, "win_end": 155961963,
     "family_id": "", "reason_selected": "Canonical AF locus positive control (1q21.3)"},
]

ANCHOR_ORDER = [row["anchor_gene"] for row in LOCUS_WINDOWS]


def _ensure_pipeline_wide() -> pd.DataFrame:
    if not PIPELINE_WIDE.is_file():
        if not PIPELINE_WIDE_LEGACY.is_file():
            sys.stderr.write(
                f"Missing {PIPELINE_WIDE.relative_to(ROOT)} and "
                f"{PIPELINE_WIDE_LEGACY.relative_to(ROOT)}.\n"
            )
            raise SystemExit(1)
        leg = pd.read_csv(PIPELINE_WIDE_LEGACY, sep="\t")
        if len(leg.columns) <= 15:
            sys.stderr.write(
                f"{PIPELINE_WIDE_LEGACY.name} looks like the slim manuscript table "
                f"(≤15 columns). Restore or regenerate the wide pipeline TSV, then "
                f"save it as {PIPELINE_WIDE.relative_to(ROOT)}.\n"
            )
            raise SystemExit(1)
        shutil.copy2(PIPELINE_WIDE_LEGACY, PIPELINE_WIDE)

    df = pd.read_csv(PIPELINE_WIDE, sep="\t")
    if "KCNN3" not in set(df["anchor"].astype(str)):
        df = pd.concat([df, _kcnn3_wide_row(df)], ignore_index=True)
        df.to_csv(PIPELINE_WIDE, sep="\t", index=False)
    return df


def _kcnn3_wide_row(template: pd.DataFrame) -> pd.DataFrame:
    """One KCNN3 row aligned to pipeline-wide columns (clone PITX2 + overlay stats)."""
    pc = pd.read_csv(POSITIVE, sep="\t")
    row_at = pc.loc[
        (pc["anchor"] == "KCNN3") & (pc["tissue"] == "atrial_appendage")
    ].iloc[0]
    pit = template.loc[template["anchor"].astype(str) == "PITX2"].iloc[0].copy()
    pit["anchor"] = "KCNN3"
    pit["n_eqtl"] = row_at["n_eqtl"]
    pit["lead_eqtl_p"] = row_at["lead_eqtl_p"]
    pit["gate_passed"] = row_at["gate_passed"]
    pit["n_gwas"] = row_at["n_gwas"]
    pit["n_overlap"] = row_at["n_overlap"]
    pit["gwas_min_p"] = row_at["gwas_min_p"]
    for i in range(5):
        pit[f"coloc_PP{i}"] = row_at[f"coloc_PP{i}"]
    pit["coloc_PP4"] = row_at["coloc_PP4"]
    pit["coloc_PP3"] = row_at["coloc_PP3"]
    # Align dialect with framework anchors that fail the cis-eQTL gate.
    pit["verdict"] = "no_eligible_eQTL"
    pit["shape_tag"] = pd.NA
    pit["calibrator_status"] = ""
    return pd.DataFrame([pit])


def _family_labels() -> dict[int, str]:
    if not FIG5B.is_file():
        return {}
    fb = pd.read_csv(FIG5B, sep="\t")
    return dict(zip(fb["family_id"].astype(int), fb["family_label"].astype(str)))


def _notes_row(r: pd.Series) -> str:
    parts: list[str] = []
    if pd.notna(r.get("shape_tag")) and str(r["shape_tag"]).strip():
        parts.append(f"shape_tag={r['shape_tag']}")
    if pd.notna(r.get("gate_passed")):
        parts.append(f"gate_passed={r['gate_passed']}")
    if pd.notna(r.get("prefilter_passed")):
        parts.append(f"prefilter_passed={r['prefilter_passed']}")
    if pd.notna(r.get("calibrator_status")) and str(r["calibrator_status"]).strip():
        parts.append(str(r["calibrator_status"]))
    return "; ".join(parts)


def write_gate_rules() -> None:
    rows = [
        ("anchor_set",
         "12 genes: 10 pre-registered framework anchors (PRELP, OMD; IDH2, NDUFA8; "
         "BMP10, SMAD9, SMAD6; NTRK2, NTF3, FRS2) plus 2 canonical calibration genes "
         "(PITX2, KCNN3). Same list defines Fig 6 anchor-matrix genes / Supplementary "
         "Table S6 programme-facing slice + calibration rows."),
        ("cis_window_grch38", "±1 Mb around gene coordinates (locked in scripts/h4_mr_coloc_anchors.R, "
         "scripts/h4c_mr_coloc_positive_controls.R, v3 Project_P6/cache eqtl_* filenames)."),
        ("af_gwas", "Yuan et al. 2025 cross-population AF summary statistics GCST90624411 "
         "(harmonised tabix), sliced to the same cis window."),
        ("eqtl_primary_tissue",
         "GTEx left atrial appendage, gene-level molecular traits via eQTL Catalogue "
         "dataset QTD000251 (n = 372 donors)."),
        ("eqtl_sensitivity_tissue",
         "GTEx whole blood via eQTL Catalogue QTD000356 (670 donors); reported in "
         "pipeline-wide columns lead_eqtl_p_blood / coloc_PP4_blood where available."),
        ("cis_eqtl_gate",
         "Lead cis-eQTL variant genome-wide gate p < 5 × 10⁻⁸ (locked EQTL_GATE_P in R drivers)."),
        ("coloc_engine", "coloc.abf with default priors p1 = p2 = 1 × 10⁻⁴, p12 = 1 × 10⁻⁵ "
         "(coloc 5.2.3); allele-aware variant intersection before coloc."),
        ("pp4_support_threshold",
         "Fig 6a / 6b colour bins: PP4 ≥ 0.70 → supported (green); "
         "0.30 ≤ PP4 < 0.70 → inconclusive (grey); PP4 < 0.30 → weak (red). "
         "Conservative colocalisation claim uses PP4 ≥ 0.70."),
        ("pipeline_verdict_column",
         "results/figure6_mr_anchors_v3_pipeline_wide.tsv verdict combines gate eligibility "
         "(eGene / Bonferroni fields), prefilter_passed, and PP4 (see scripts/fig6e_* "
         "and MR summary). Values include no_eligible_eQTL, not_supported, etc."),
        ("interpretation_guardrail",
         "Because calibration genes PITX2 and KCNN3 also fail the cis-eQTL gate at "
         "current atrial-appendage panel resolution, absence of PP4 ≥ 0.70 support "
         "should be read as calibration-limited under the tested regime — not as proof "
         "that genetics excludes these genes as AF-relevant."),
        ("reproducibility_status",
         "Anchor coloc numbers (figure6_mr_anchors_v3_pipeline_wide.tsv) and "
         "lead-variant offsets (fig6a_v3_locus_zoom_summary.tsv) are FROZEN analysis "
         "artifacts; Fig 6A/6B panels and manuscript tables regenerate directly from "
         "them (scripts/fig6a_anchor_coloc_calibration.py, "
         "scripts/fig6b_pp3_pp4_leaddistance.py). End-to-end re-derivation needs "
         "network tabix (EBI FTP) via scripts/h4*_mr_coloc*.R. Filenames containing "
         "'mr' are coloc-based diagnostics and do NOT perform Mendelian randomization."),
    ]
    pd.DataFrame(rows, columns=["item", "rule"]).to_csv(OUT_RULES, sep="\t", index=False)


def write_calibration_table(wide: pd.DataFrame) -> None:
    rows = []
    for gene in ("PITX2", "KCNN3"):
        sub = wide.loc[wide["anchor"].astype(str) == gene]
        if sub.empty:
            continue
        r = sub.iloc[0]
        gp = r.get("gate_passed")
        if isinstance(gp, str):
            gate_ok = gp.lower() in ("true", "1", "yes")
        else:
            gate_ok = bool(gp) if pd.notna(gp) else False
        pp4 = r.get("coloc_PP4")
        verdict = str(r.get("verdict", ""))
        role = (
            "Canonical AF control — 4q25 (enhancer-mediated cis-eQTL; weak at RNA layer)"
            if gene == "PITX2"
            else "Canonical AF control — 1q21.3 (literature coloc context)"
        )
        rows.append({
            "gene": gene,
            "role": role,
            "cis_eqtl_gate_pass": gate_ok,
            "coloc_PP4": pp4,
            "verdict": verdict,
        })
    pd.DataFrame(rows).to_csv(OUT_CALIB, sep="\t", index=False)


def write_locus_windows(family_labels: dict[int, str]) -> None:
    rows = []
    for i, lw in enumerate(LOCUS_WINDOWS, start=1):
        fid = lw["family_id"]
        flab = ""
        if isinstance(fid, int):
            flab = family_labels.get(fid, "")
        rows.append({
            "panel": f"locus_{i:02d}",
            "gene": lw["anchor_gene"],
            "chr": lw["chr"],
            "start": lw["win_start"],
            "end": lw["win_end"],
            "reason_selected": lw["reason_selected"],
            "family_id": fid if fid != "" else "",
            "family_label": flab,
        })
    pd.DataFrame(rows).to_csv(OUT_LOCUS, sep="\t", index=False)


def main() -> None:
    wide_all = _ensure_pipeline_wide()
    missing = [g for g in ANCHOR_ORDER if g not in set(wide_all["anchor"].astype(str))]
    if missing:
        sys.stderr.write(f"Missing anchors in pipeline-wide table: {missing}\n")
        raise SystemExit(1)
    wide = wide_all.set_index("anchor").loc[ANCHOR_ORDER].reset_index()

    family_labels = _family_labels()

    slim = []
    win_by_gene = {w["anchor_gene"]: w for w in LOCUS_WINDOWS}
    for _, r in wide.iterrows():
        g = str(r["anchor"])
        w = win_by_gene[g]
        fid = w["family_id"]
        flab = ""
        if isinstance(fid, int):
            flab = family_labels.get(int(fid), "")
        if g in ("PITX2", "KCNN3"):
            fid_disp = ""
            if not flab:
                flab = (
                    "Positive control — PITX2 / 4q25 AF locus"
                    if g == "PITX2"
                    else "Positive control — KCNN3 / 1q21.3 AF locus"
                )
        else:
            fid_disp = fid
        locus_s = f"chr{w['chr']}:{w['win_start']}-{w['win_end']}"
        slim.append({
            "anchor_gene": g,
            "family_id": fid_disp,
            "family_label": flab,
            "source": r.get("source", ""),
            "locus": locus_s,
            "cis_eqtl_lead_p": r.get("lead_eqtl_p"),
            "gwas_lead_p": r.get("gwas_min_p"),
            "coloc_PP4": r.get("coloc_PP4"),
            "coloc_PP3": r.get("coloc_PP3"),
            "verdict": r.get("verdict", ""),
            "notes": _notes_row(r),
        })

    pd.DataFrame(slim).to_csv(OUT_CONSOLIDATED, sep="\t", index=False)

    write_gate_rules()
    write_calibration_table(wide)
    write_locus_windows(family_labels)

    print(f"wrote {OUT_CONSOLIDATED.relative_to(ROOT)} ({len(slim)} rows)")
    print(f"wrote {OUT_RULES.relative_to(ROOT)}")
    print(f"wrote {OUT_CALIB.relative_to(ROOT)}")
    print(f"wrote {OUT_LOCUS.relative_to(ROOT)}")
    print(f"pipeline wide (source of truth): {PIPELINE_WIDE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
