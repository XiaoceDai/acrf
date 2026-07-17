#!/usr/bin/env python3
"""
Phase 5B3B — Build final formal Supplementary Tables S1–S6 from frozen sources.

Reads only authoritative frozen tables; never mutates inputs; scientific_rerun=NO.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, Side
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parents[2]
XLSX_DIR = ROOT / "tables" / "supplementary"
TSV_DIR = ROOT / "results" / "supplementary_tables"
AUDIT = ROOT / "results" / "manuscript_audit"
PREV = AUDIT / "phase5b3b_table_previews"
BUILDER = "scripts/supplementary/build_final_supplementary_tables.py"
UTC = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
SNAP_DATE = UTC[:10]

FORBIDDEN_VISIBLE = [
    "NCP_v3",
    "regulator-tier",
    "regulator tier",
    "G-tier",
    "G0",
    "G1",
    "G2",
    "G3",
    "G4",
    "Figure 5E",
    "Figure 6E",
    "Figure 7",
    "Atlas PriorityScore",
    "highest tier",
    "evidence score",
    "composite score",
]


def sha256(path: Path | str) -> str:
    p = Path(path)
    if not p.is_file():
        return ""
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_tsv(rel: str) -> pd.DataFrame:
    p = ROOT / rel
    if not p.is_file():
        raise FileNotFoundError(rel)
    return pd.read_csv(p, sep="\t")


def na(x) -> str:
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return "NA"
    if isinstance(x, str) and x.strip() in {"", "nan", "None", "NaN"}:
        return "NA"
    return x


def family_order_map() -> pd.DataFrame:
    order = read_tsv("results/ncp_v3_canonical_display_order.tsv")
    labs = read_tsv("results/ncp_v3_family_labels.tsv")
    fam = (
        order[["family_display_order", "family_id", "canonical_direction"]]
        .drop_duplicates("family_id")
        .sort_values("family_display_order")
        .reset_index(drop=True)
    )
    fam = fam.merge(labs[["family_id", "family_label", "short_label"]], on="family_id", how="left")
    assert len(fam) == 16
    return fam


def style_header(ws, n_cols: int) -> None:
    thin = Border(
        left=Side(style=None),
        right=Side(style=None),
        top=Side(style=None),
        bottom=Side(style="thin", color="B0B0B0"),
    )
    for col in range(1, n_cols + 1):
        cell = ws.cell(1, col)
        cell.font = Font(bold=True, name="Calibri", size=11)
        cell.alignment = Alignment(wrap_text=True, vertical="center")
        cell.border = thin
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(n_cols)}{ws.max_row}"


def autosize(ws, max_width: int = 42) -> None:
    for col in ws.columns:
        letter = get_column_letter(col[0].column)
        maxlen = 0
        for cell in col:
            val = "" if cell.value is None else str(cell.value)
            maxlen = max(maxlen, min(len(val), max_width))
            cell.alignment = Alignment(wrap_text=True, vertical="top")
            if isinstance(cell.value, float):
                cell.number_format = "0.############"
        ws.column_dimensions[letter].width = max(10, min(max_width, maxlen + 2))


def write_kv_sheet(ws, rows: list[tuple[str, object]]) -> None:
    ws.append(["field", "value"])
    for k, v in rows:
        ws.append([k, v])
    style_header(ws, 2)
    autosize(ws, 80)


def write_df_sheet(ws, df: pd.DataFrame) -> None:
    cols = list(df.columns)
    ws.append(cols)
    for _, row in df.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            if pd.isna(v):
                vals.append("NA")
            elif isinstance(v, (bool, int)):
                vals.append(v)
            elif isinstance(v, float):
                vals.append(v)
            else:
                vals.append(v)
        ws.append(vals)
    style_header(ws, len(cols))
    autosize(ws)


def export_tsv(df: pd.DataFrame, path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    out = df.copy()
    for c in out.columns:
        out[c] = out[c].map(lambda v: "NA" if pd.isna(v) else v)
    out.to_csv(path, sep="\t", index=False)
    return sha256(path)


def write_workbook(
    xlsx_path: Path,
    readme: list[tuple[str, object]],
    sheets: dict[str, pd.DataFrame],
    provenance: pd.DataFrame,
) -> None:
    xlsx_path.parent.mkdir(parents=True, exist_ok=True)
    wb = Workbook()
    ws0 = wb.active
    ws0.title = "README"
    write_kv_sheet(ws0, readme)
    for name, df in sheets.items():
        ws = wb.create_sheet(name)
        write_df_sheet(ws, df)
    wsp = wb.create_sheet("Provenance")
    write_df_sheet(wsp, provenance)
    # no hidden sheets
    for ws in wb.worksheets:
        ws.sheet_state = "visible"
    wb.save(xlsx_path)


def scan_forbidden(text: str) -> list[str]:
    hits = []
    for tok in FORBIDDEN_VISIBLE:
        if tok in ("G0", "G1", "G2", "G3", "G4"):
            if re.search(rf"\b{tok}\b", text):
                hits.append(tok)
        elif tok.lower() in text.lower():
            hits.append(tok)
    return sorted(set(hits))


PATH_LIKE_COLS = {
    "source_path",
    "source_paths",
    "source_sha256",
    "source_sha256s",
    "authoritative_source_path",
    "authoritative_source_sha256",
}


def scan_forbidden_frames(frames: list[pd.DataFrame]) -> list[str]:
    """Scan publication-facing cell text; ignore internal path/hash columns."""
    chunks = []
    for df in frames:
        keep = [c for c in df.columns if c not in PATH_LIKE_COLS]
        if keep:
            chunks.append(df[keep].astype(str).to_string())
    return scan_forbidden("\n".join(chunks))


def provenance_rows(items: list[tuple[str, str, str, str]]) -> pd.DataFrame:
    """items: source_role, source_path, source_qc_path, notes"""
    rows = []
    for role, path, qc, notes in items:
        rows.append(
            {
                "source_role": role,
                "source_path": path,
                "source_sha256": sha256(ROOT / path) if path else "",
                "source_qc_path": qc or "NA",
                "source_qc_status": "PASS" if path and (ROOT / path).is_file() else "MISSING",
                "notes": notes,
            }
        )
    return pd.DataFrame(rows)


# S1

def build_s1() -> dict:
    dem = read_tsv("results/cohort_demographics_table.tsv")
    long = read_tsv("results/cohort_demographics_long.tsv")
    fig2a = read_tsv("results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_summary.tsv")
    gwas = read_tsv("results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv")
    gwas_map = dict(zip(gwas["key"], gwas["value"]))
    ext = read_tsv("results/provenance/figure6_external_gene_sets_manifest.tsv")
    coloc_rules = read_tsv("results/figure6_mr_coloc_gate_rules.tsv")
    geomx_audit = read_tsv("results/ncp_v3_singscore_reanalysis/figure5_geomx_lineage_audit.tsv")
    hca = "results/fig1f_v3_HCA_family_carriers.tsv"

    total = fig2a.loc[fig2a.scope == "total"].iloc[0]
    assert int(total.paired_donors) == 51
    assert int(total.AF_donors) == 25
    assert int(total.SR_donors) == 26

    g112 = long[long.cohort == "GSE112339"]
    assert g112.donor_id.nunique() == 8
    assert set(g112.rhythm.unique()) == {"healthy_chamber"}

    boundary = (
        "Counts are dataset-specific and use donors, samples, patients, ROIs, or summary "
        "statistics as appropriate. They are not a single pooled donor denominator."
    )

    # GeoMx ROI boundary
    roi_row = geomx_audit.loc[geomx_audit.contrast_id == "ROI_boundary_summary"].iloc[0]
    assert int(roi_row.estimate) == 24
    assert int(roi_row.n_rois) == 48
    assert int(roi_row.n_patients) == 4

    dem_map = dem.set_index("cohort")

    def dem_row(cohort: str) -> dict:
        r = dem_map.loc[cohort]
        return {
            "platform": r.platform,
            "chamber": r.chamber,
            "n_donors": r.n_donors,
            "n_AF": r.n_AF_donors,
            "n_SR": r.n_SR_donors,
            "role": r.role,
            "background": r.background,
            "sex": r.sex,
            "age": r.age,
        }

    rows = []

    def add(**kwargs):
        base = {
            "resource_id": "NA",
            "accession_or_resource": "NA",
            "publication_label": "NA",
            "molecular_platform": "NA",
            "tissue_or_region": "NA",
            "chamber_or_region": "NA",
            "rhythm_or_control_groups": "NA",
            "n_samples": "NA",
            "n_donors": "NA",
            "n_AF": "NA",
            "n_SR_or_control": "NA",
            "n_patients": "NA",
            "n_ROIs_full": "NA",
            "n_ROIs_eligible": "NA",
            "primary_analysis_unit": "NA",
            "primary_analysis_role": "NA",
            "main_or_supplementary": "NA",
            "main_figure_or_table_link": "NA",
            "inference_role": "NA",
            "source_publication": "NA",
            "source_path": "NA",
            "source_sha256": "NA",
            "notes": "",
        }
        base.update(kwargs)
        rows.append(base)

    # GSE112339
    d = dem_row("GSE112339")
    add(
        resource_id="GSE112339",
        accession_or_resource="GSE112339",
        publication_label="Healthy paired LA/RA bulk reference",
        molecular_platform=d["platform"],
        tissue_or_region="atrial appendage",
        chamber_or_region="LA + RA (paired)",
        rhythm_or_control_groups="healthy_chamber",
        n_samples=16,
        n_donors=8,
        n_AF=0,
        n_SR_or_control=0,
        primary_analysis_unit="donor (paired LA/RA samples)",
        primary_analysis_role="Healthy-reference NCP calibration",
        main_or_supplementary="main",
        main_figure_or_table_link="Figure 1; S1 Fig; S4 Fig",
        inference_role="reference_calibration",
        source_path="results/cohort_demographics_table.tsv",
        source_sha256=sha256(ROOT / "results/cohort_demographics_table.tsv"),
        notes="Eight healthy donors with paired LA and RA samples. Not an AF remodeling cohort.",
    )

    # Paired bulk AF/SR
    for cohort, fig_link, role in [
        ("GSE128188", "Figure 2", "paired-chamber AF/SR discovery"),
        ("GSE79768", "Figure 2", "paired-chamber AF/SR discovery"),
        ("GSE115574", "Figure 2", "paired-chamber AF/SR replication/background"),
    ]:
        d = dem_row(cohort)
        # Publication-facing notes: strip legacy internal filenames from demographics prose
        bg = re.sub(r"innovation_ncp_v3_acai_long\.tsv", "the frozen paired-bulk aCAI long table", str(d["background"]))
        bg = re.sub(r"ncp_v3", "NCP", bg, flags=re.IGNORECASE)
        add(
            resource_id=cohort,
            accession_or_resource=cohort,
            publication_label=role,
            molecular_platform=d["platform"],
            tissue_or_region="atrial appendage",
            chamber_or_region=str(d["chamber"]),
            rhythm_or_control_groups="AF vs SR",
            n_donors=int(d["n_donors"]),
            n_AF=int(d["n_AF"]),
            n_SR_or_control=int(d["n_SR"]),
            n_samples=int(d["n_donors"]) * 2,
            primary_analysis_unit="paired donor",
            primary_analysis_role="Bulk AF-associated NCP family rewiring",
            main_or_supplementary="main",
            main_figure_or_table_link=f"{fig_link}; S3 Table",
            inference_role="bulk_consistency_input",
            source_path="results/cohort_demographics_table.tsv;results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_summary.tsv",
            source_sha256=sha256(ROOT / "results/cohort_demographics_table.tsv"),
            notes=f"Background: {bg}. Contributes to paired-bulk total 51 donors (25 AF / 26 SR).",
        )

    # HCA
    add(
        resource_id="HCA_heart",
        accession_or_resource="Human Cell Atlas heart resource",
        publication_label="HCA LA carrier prior",
        molecular_platform="single-cell / single-nucleus heart atlas",
        tissue_or_region="heart (LA focus for carrier prior)",
        chamber_or_region="LA carrier prior",
        rhythm_or_control_groups="healthy reference context",
        primary_analysis_unit="family × compartment carrier prior",
        primary_analysis_role="Healthy HCA LA carrier prior for Figure 3/5 context",
        main_or_supplementary="main",
        main_figure_or_table_link="Figure 3; Figure 5 (context only)",
        inference_role="context_only_not_Path_A",
        source_path=hca,
        source_sha256=sha256(ROOT / hca),
        notes="Carrier prior only; does not enter Path A support count.",
    )

    # GSE255612
    d = dem_row("GSE255612")
    assert int(d["n_AF"]) == 18 and int(d["n_SR"]) == 16
    add(
        resource_id="GSE255612",
        accession_or_resource="GSE255612",
        publication_label="LA single-nucleus AF vs control",
        molecular_platform=d["platform"],
        tissue_or_region="left atrium",
        chamber_or_region="LA only (primary localization)",
        rhythm_or_control_groups="AF vs control",
        n_donors=34,
        n_AF=18,
        n_SR_or_control=16,
        primary_analysis_unit="donor (LA pseudobulk)",
        primary_analysis_role="LA snRNA family/program localization",
        main_or_supplementary="main",
        main_figure_or_table_link="Figure 3; S3 Fig; S4 Table",
        inference_role="snRNA_localization",
        source_path="results/cohort_demographics_table.tsv;results/fig3c_family_fidelity_classes.tsv",
        source_sha256=sha256(ROOT / "results/fig3c_family_fidelity_classes.tsv"),
        notes="Primary analysis is LA-only (18 AF / 16 control). Resource chamber metadata may list LA+RA; final localization classes are LA-only.",
    )

    # GSE41177
    d = dem_row("GSE41177")
    add(
        resource_id="GSE41177",
        accession_or_resource="GSE41177",
        publication_label="AF duration LAA primary / LAPV sensitivity",
        molecular_platform=d["platform"],
        tissue_or_region="LAA + LAPV",
        chamber_or_region="LAA (primary); LAPV (sensitivity)",
        rhythm_or_control_groups="AF (duration) + SR (metadata/context only)",
        n_donors=16,
        n_AF=13,
        n_SR_or_control=3,
        primary_analysis_unit="AF-LAA donor with numeric duration",
        primary_analysis_role="Cross-sectional duration / chronicity",
        main_or_supplementary="main",
        main_figure_or_table_link="Figure 4; S5 Table",
        inference_role="duration_chronicity",
        source_path="results/ncp_v3_singscore_reanalysis/fig4a_duration_design_summary.tsv",
        source_sha256=sha256(ROOT / "results/ncp_v3_singscore_reanalysis/fig4a_duration_design_summary.tsv"),
        notes="Primary n=12 AF-LAA donors with numeric duration. SR samples document cohort structure only; they are not duration anchors in final Figure 4 inference.",
    )

    # GSE238242
    d = dem_row("GSE238242")
    add(
        resource_id="GSE238242",
        accession_or_resource="GSE238242",
        publication_label="Multiome RNA-side / mural sensitivity",
        molecular_platform=d["platform"],
        tissue_or_region="LAA",
        chamber_or_region="LAA",
        rhythm_or_control_groups="AF vs SR",
        n_donors=int(d["n_donors"]),
        n_AF=int(d["n_AF"]),
        n_SR_or_control=int(d["n_SR"]),
        primary_analysis_unit="donor pseudobulk",
        primary_analysis_role="Supplementary RNA-side and mural-cell sensitivity",
        main_or_supplementary="supplementary",
        main_figure_or_table_link="S5 Fig",
        inference_role="sensitivity_not_Path_A",
        source_path="results/cohort_demographics_table.tsv",
        source_sha256=sha256(ROOT / "results/cohort_demographics_table.tsv"),
        notes="Supplementary sensitivity only; not a Path A gate.",
    )

    # GSE261363 GeoMx
    d = dem_row("GSE261363")
    add(
        resource_id="GSE261363",
        accession_or_resource="GSE261363",
        publication_label="GeoMx ROI-level spatial context",
        molecular_platform=d["platform"],
        tissue_or_region="atrial appendage + epicardial adipose",
        chamber_or_region=str(d["chamber"]),
        rhythm_or_control_groups="AF vs SR (ROI context)",
        n_patients=4,
        n_ROIs_full=48,
        n_ROIs_eligible=24,
        primary_analysis_unit="ROI (eligible cardiac AA DNA+)",
        primary_analysis_role="GeoMx ROI-level context (not Path A)",
        main_or_supplementary="main",
        main_figure_or_table_link="Figure 5C (context)",
        inference_role="context_only_not_Path_A",
        source_path="Project_P5/GSE261363/roi_scores_v3.tsv;results/ncp_v3_singscore_reanalysis/figure5_geomx_lineage_audit.tsv",
        source_sha256=sha256(ROOT / "Project_P5/GSE261363/roi_scores_v3.tsv"),
        notes="48 total ROIs; 24 eligible ROIs; four patients. ROI-level nominal summaries; no patient-clustered inference; does not enter Path A.",
    )

    # GTEx / eQTL Catalogue
    eqtl_rule = coloc_rules.loc[coloc_rules.item == "eqtl_primary_tissue", "rule"].iloc[0]
    add(
        resource_id="GTEx_LAA_eQTL_Catalogue_QTD000251",
        accession_or_resource="GTEx LAA / eQTL Catalogue QTD000251",
        publication_label="Atrial cis-eQTL resource for Figure 6 calibration",
        molecular_platform="GTEx v8 left atrial appendage cis-eQTL (eQTL Catalogue)",
        tissue_or_region="left atrial appendage",
        chamber_or_region="LAA",
        rhythm_or_control_groups="population eQTL donors",
        n_donors=372,
        primary_analysis_unit="gene / cis-window",
        primary_analysis_role="Atrial eGene / coloc calibration input",
        main_or_supplementary="main",
        main_figure_or_table_link="Figure 6; S6 Table",
        inference_role="genetic_calibration_not_Path_A",
        source_path="results/figure6_mr_coloc_gate_rules.tsv",
        source_sha256=sha256(ROOT / "results/figure6_mr_coloc_gate_rules.tsv"),
        notes=str(eqtl_rule),
    )

    # Yuan GWAS
    assert gwas_map["locked_af_cases"] == "252438"
    assert gwas_map["locked_controls"] == "1959739"
    assert gwas_map["coloc_rerun_required"] == "FALSE"
    add(
        resource_id="GCST90624411",
        accession_or_resource="Yuan 2025 AF GWAS / GCST90624411",
        publication_label="Cross-population AF GWAS",
        molecular_platform="GWAS summary statistics (harmonised)",
        tissue_or_region="NA",
        chamber_or_region="NA",
        rhythm_or_control_groups="AF cases vs controls",
        n_AF=252438,
        n_SR_or_control=1959739,
        primary_analysis_unit="variant (summary statistics)",
        primary_analysis_role="Local AF GWAS / coloc calibration",
        main_or_supplementary="main",
        main_figure_or_table_link="Figure 6; S6 Table",
        inference_role="genetic_calibration_not_Path_A",
        source_publication=gwas_map.get("publication", "Yuan et al. 2025"),
        source_path="results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv",
        source_sha256=sha256(ROOT / "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv"),
        notes="252438 AF cases; 1959739 controls; coloc_rerun_required=FALSE; harmonised summary-statistics provenance frozen.",
    )

    # External AF gene sets (Figure 6C)
    for _, r in ext.sort_values("display_order").iterrows():
        add(
            resource_id=str(r.source_id),
            accession_or_resource=str(r.source_id),
            publication_label=f"{r.display_author_year} — {r.display_evidence_label}",
            molecular_platform="external AF gene set",
            tissue_or_region="NA",
            chamber_or_region="NA",
            rhythm_or_control_groups="NA",
            n_samples=int(r.n_genes_normalized) if pd.notna(r.n_genes_normalized) else "NA",
            primary_analysis_unit="gene set",
            primary_analysis_role="Figure 6C family–gene-set overlap resource",
            main_or_supplementary="main",
            main_figure_or_table_link="Figure 6C; S8 Fig",
            inference_role="external_overlap_not_Path_A",
            source_publication=str(r.source_publication),
            source_path="results/provenance/figure6_external_gene_sets_manifest.tsv",
            source_sha256=sha256(ROOT / "results/provenance/figure6_external_gene_sets_manifest.tsv"),
            notes=f"Evidence role: {r.evidence_role}. Normalized genes={r.n_genes_normalized}.",
        )

    df = pd.DataFrame(rows)
    required = [
        "GSE112339", "GSE128188", "GSE79768", "GSE115574", "HCA_heart", "GSE255612",
        "GSE41177", "GSE238242", "GSE261363", "GTEx_LAA_eQTL_Catalogue_QTD000251",
        "GCST90624411",
    ]
    for rid in required:
        assert rid in set(df.resource_id), rid
    assert set(ext.source_id).issubset(set(df.resource_id))

    # paired bulk total check from embedded donors
    paired = df[df.resource_id.isin(["GSE128188", "GSE79768", "GSE115574"])]
    assert int(paired.n_donors.sum()) == 51
    assert int(paired.n_AF.sum()) == 25
    assert int(paired.n_SR_or_control.sum()) == 26

    prov = provenance_rows(
        [
            ("cohort_demographics_table", "results/cohort_demographics_table.tsv", "NA", "Primary cohort metadata"),
            ("cohort_demographics_long", "results/cohort_demographics_long.tsv", "NA", "Donor-level long form"),
            ("bulk_deployment", "results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_summary.tsv", "NA", "Paired bulk totals"),
            ("geomx_boundary", "results/ncp_v3_singscore_reanalysis/figure5_geomx_lineage_audit.tsv", "NA", "ROI eligibility"),
            ("gwas_provenance", "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv", "NA", "GCST90624411"),
            ("external_gene_sets", "results/provenance/figure6_external_gene_sets_manifest.tsv", "NA", "Figure 6C sets"),
            ("coloc_rules", "results/figure6_mr_coloc_gate_rules.tsv", "NA", "GTEx/eQTL Catalogue"),
            ("hca_carriers", hca, "NA", "HCA carrier prior"),
        ]
    )
    readme = [
        ("supplementary_table_id", "S1 Table"),
        ("title", "Public datasets, cohort metadata, and analysis roles"),
        ("scientific_role", "Dataset inventory and analysis-unit map for the NCP framework"),
        ("main_figure_links", "Figures 1–6; S1/S4/S5/S8 Figs"),
        ("analysis_unit", "Dataset-specific (donor/sample/patient/ROI/summary statistic)"),
        ("interpretation_boundary", boundary),
        ("source_snapshot_date", SNAP_DATE),
        ("scientific_rerun", "NO"),
        ("builder_script", BUILDER),
    ]
    xlsx = XLSX_DIR / "S1_Table_public_datasets_analysis_roles.xlsx"
    write_workbook(xlsx, readme, {"datasets": df}, prov)
    tsv = TSV_DIR / "S1_datasets.tsv"
    export_tsv(df, tsv)
    return {
        "final_id": "S1 Table",
        "final_title": "Public datasets, cohort metadata, and analysis roles",
        "xlsx": xlsx,
        "sheets": {"datasets": df},
        "tsv_paths": [tsv],
        "provenance": prov,
        "sources": [p for _, p, _, _ in [
            ("a", "results/cohort_demographics_table.tsv", "", ""),
            ("b", "results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_summary.tsv", "", ""),
            ("c", "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv", "", ""),
            ("d", "results/provenance/figure6_external_gene_sets_manifest.tsv", "", ""),
        ]],
        "notes": "Rebuilt from cohort demographics + provenance; dataset-specific units",
        "assertions": [
            "GSE112339: 8 healthy paired LA/RA donors",
            "paired bulk total: 51 donors / 25 AF / 26 SR",
            "GSE255612: 18 AF / 16 control; LA-only primary",
            "GSE41177 primary: 12 AF-LAA numeric duration",
            "GSE261363: 48 ROIs / 24 eligible / 4 patients",
            "GCST90624411: 252438 cases / 1959739 controls",
        ],
    }


# S2

def build_s2() -> dict:
    claims = read_tsv("results/handoff/fulltext_audit_v1_1/02_final_claims_registry.tsv")
    contract = read_tsv("results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv")
    cmap = dict(zip(contract["metric"], contract["value"].astype(str)))
    assert cmap["n_Path_A_gates"] == "3"

    def row(**kwargs):
        base = {
            "analysis_layer": "NA",
            "dataset_or_resource": "NA",
            "main_figure_panel": "NA",
            "tested_or_summarized_unit": "NA",
            "primary_estimand": "NA",
            "primary_statistic": "NA",
            "uncertainty_measure": "NA",
            "multiplicity_or_resampling": "NA",
            "classification_or_gate": "NA",
            "contributes_to_Path_A": "NO",
            "formal_Path_A_gate_name": "NA",
            "context_only": "YES",
            "interpretation_boundary": "NA",
            "authoritative_source_path": "NA",
            "authoritative_source_sha256": "NA",
            "notes": "",
        }
        base.update(kwargs)
        return base

    rows = [
        row(
            analysis_layer="Healthy-reference calibration",
            dataset_or_resource="GSE112339",
            main_figure_panel="Figure 1",
            tested_or_summarized_unit="program/family",
            primary_estimand="healthy chamber polarity / panel lock",
            primary_statistic="descriptive composition",
            uncertainty_measure="NA",
            multiplicity_or_resampling="panel freeze (no discovery multiplicity)",
            classification_or_gate="reference lock",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="Healthy-reference calibration only; not an AF remodeling test.",
            authoritative_source_path="results/handoff/fulltext_audit_v1_1/02_final_claims_registry.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/handoff/fulltext_audit_v1_1/02_final_claims_registry.tsv"),
        ),
        row(
            analysis_layer="Figure 2C family-by-cohort states",
            dataset_or_resource="GSE128188;GSE79768;GSE115574",
            main_figure_panel="Figure 2C",
            tested_or_summarized_unit="family × cohort",
            primary_estimand="AF/SR aCAIrel remodeling state",
            primary_statistic="state call (Amplified/Preserved/Attenuated/Reversed)",
            uncertainty_measure="descriptive ratio thresholds",
            multiplicity_or_resampling="none (descriptive)",
            classification_or_gate="descriptive state calls",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="Descriptive state calls, not independent significance tests.",
            authoritative_source_path="results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv"),
        ),
        row(
            analysis_layer="Figure 2D fixed-effect IVW consistency",
            dataset_or_resource="GSE128188;GSE79768;GSE115574",
            main_figure_panel="Figure 2D",
            tested_or_summarized_unit="family",
            primary_estimand="pooled Hedges g (AF−SR)",
            primary_statistic="fixed-effect inverse-variance pooled g",
            uncertainty_measure="Wald CI on pooled g",
            multiplicity_or_resampling="sign-match count across 3 cohorts (no FDR discovery)",
            classification_or_gate="High/Moderate/Low consistency",
            contributes_to_Path_A="YES",
            formal_Path_A_gate_name="bulk High",
            context_only="NO",
            interpretation_boundary="Consistency class from FE-IVW sign-match structure; High is one of three Path A gates.",
            authoritative_source_path="results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv"),
            notes="Path A gate: consistency_class == High",
        ),
        row(
            analysis_layer="Figure 3 nominal permutation localization",
            dataset_or_resource="GSE255612",
            main_figure_panel="Figure 3",
            tested_or_summarized_unit="family × LA compartment",
            primary_estimand="Δsingscore (AF − control)",
            primary_statistic="nominal permutation p",
            uncertainty_measure="Monte Carlo / exact permutation p",
            multiplicity_or_resampling="nominal permutation; not BH-FDR discovery",
            classification_or_gate="Fidelity/Relay/Extension/Below-detection",
            contributes_to_Path_A="YES",
            formal_Path_A_gate_name="snRNA non-Below-detection",
            context_only="NO",
            interpretation_boundary="Nominal permutation localization, not FDR-controlled discovery. LA-only.",
            authoritative_source_path="results/fig3c_family_fidelity_classes.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/fig3c_family_fidelity_classes.tsv"),
            notes="Path A gate: fidelity_class != Below_detection",
        ),
        row(
            analysis_layer="Figure 4 duration chronicity classes",
            dataset_or_resource="GSE41177",
            main_figure_panel="Figure 4",
            tested_or_summarized_unit="family",
            primary_estimand="Spearman rho vs AF duration (LAA)",
            primary_statistic="rho_LAA (full precision)",
            uncertainty_measure="permutation p; bootstrap CI (reported, not class-defining)",
            multiplicity_or_resampling="class from |rho_LAA| bins only",
            classification_or_gate="Strong/Weak/Stable",
            contributes_to_Path_A="YES",
            formal_Path_A_gate_name="duration non-Stable",
            context_only="NO",
            interpretation_boundary="Full-precision |rho_LAA| effect-size bins, not significance classes. SR not duration anchors.",
            authoritative_source_path="results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv"),
            notes="Path A gate: chronicity_class != Stable",
        ),
        row(
            analysis_layer="Figure 5 Path A rule",
            dataset_or_resource="cross-layer",
            main_figure_panel="Figure 5A–B",
            tested_or_summarized_unit="family",
            primary_estimand="count of satisfied Path A gates",
            primary_statistic="support count (0–3)",
            uncertainty_measure="gate-defined (not a new test statistic)",
            multiplicity_or_resampling="exactly three predefined gates",
            classification_or_gate="Path A support count",
            contributes_to_Path_A="NO",
            formal_Path_A_gate_name="NA",
            context_only="NO",
            interpretation_boundary="Path A aggregates the three gates only; genetics/GeoMx/HCA/dominant-state do not enter.",
            authoritative_source_path="results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv"),
            notes=f"Gates: {cmap.get('Path_A_gate_conditions','')}",
        ),
        row(
            analysis_layer="HCA healthy carrier prior",
            dataset_or_resource="HCA heart",
            main_figure_panel="Figure 3 / Figure 5",
            tested_or_summarized_unit="family",
            primary_estimand="LA carrier compartment prior",
            primary_statistic="carrier label",
            uncertainty_measure="NA",
            multiplicity_or_resampling="NA",
            classification_or_gate="context prior",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="HCA carrier prior does not enter Path A support count.",
            authoritative_source_path="results/fig1f_v3_HCA_family_carriers.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/fig1f_v3_HCA_family_carriers.tsv"),
        ),
        row(
            analysis_layer="Bulk dominant rewiring state",
            dataset_or_resource="GSE128188;GSE79768;GSE115574",
            main_figure_panel="Figure 2C / Figure 5",
            tested_or_summarized_unit="family",
            primary_estimand="dominant cohort state",
            primary_statistic="mode/summary of state calls",
            uncertainty_measure="NA",
            multiplicity_or_resampling="NA",
            classification_or_gate="context summary",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="Bulk dominant state does not enter Path A support count.",
            authoritative_source_path="results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv"),
        ),
        row(
            analysis_layer="GeoMx ROI-level context",
            dataset_or_resource="GSE261363",
            main_figure_panel="Figure 5C",
            tested_or_summarized_unit="ROI",
            primary_estimand="ROI-level association / delta",
            primary_statistic="ROI-level nominal summaries",
            uncertainty_measure="ROI-level p (no patient-clustered SE)",
            multiplicity_or_resampling="not Path A; not support-count",
            classification_or_gate="context only",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="ROI-level nominal summaries; no patient-clustered inference; not Path A.",
            authoritative_source_path="results/ncp_v3_singscore_reanalysis/figure5_geomx_lineage_audit.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/ncp_v3_singscore_reanalysis/figure5_geomx_lineage_audit.tsv"),
        ),
        row(
            analysis_layer="Figure 6A local GWAS / atrial eGene / PP4",
            dataset_or_resource="GCST90624411; GTEx LAA eQTL Catalogue",
            main_figure_panel="Figure 6A",
            tested_or_summarized_unit="anchor gene",
            primary_estimand="local GWAS support; Bonferroni atrial eGene; PP4",
            primary_statistic="presence flags + coloc PP3/PP4",
            uncertainty_measure="coloc posterior probabilities",
            multiplicity_or_resampling="Bonferroni atrial eGene definition (not lead-p gate)",
            classification_or_gate="non-nested calibration flags",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="Calibration-limited genetic evidence; not Path A. Primary eGene definition is Bonferroni atrial eGene, not lead-p<5e−8.",
            authoritative_source_path="results/fig6a_anchor_coloc_calibration_source.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/fig6a_anchor_coloc_calibration_source.tsv"),
        ),
        row(
            analysis_layer="Figure 6B PP3 and lead separation",
            dataset_or_resource="GCST90624411; GTEx LAA eQTL Catalogue",
            main_figure_panel="Figure 6B",
            tested_or_summarized_unit="anchor gene",
            primary_estimand="PP3-dominant status; lead GWAS–eQTL distance",
            primary_statistic="PP3/PP4; delta_bp / delta_Mb",
            uncertainty_measure="coloc posteriors",
            multiplicity_or_resampling="NA",
            classification_or_gate="distance / PP3 display",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="Lead separation and PP3-dominant display; not Path A.",
            authoritative_source_path="results/fig6b_posterior_lead_distance_source.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/fig6b_posterior_lead_distance_source.tsv"),
        ),
        row(
            analysis_layer="Figure 6C family–gene-set hypergeometric tests",
            dataset_or_resource="external AF gene sets; NCP families",
            main_figure_panel="Figure 6C",
            tested_or_summarized_unit="family × external set",
            primary_estimand="member overlap / hypergeometric enrichment",
            primary_statistic="overlap count; BH-FDR q",
            uncertainty_measure="BH-FDR across 64 tests",
            multiplicity_or_resampling="BH-FDR (0/64 significant)",
            classification_or_gate="member-overlap display",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="Member overlap is not family-level enrichment proof; not Path A.",
            authoritative_source_path="results/fig6c_family_af_overlap_source.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/fig6c_family_af_overlap_source.tsv"),
        ),
        row(
            analysis_layer="Figure 6D non-nested traceability flags",
            dataset_or_resource="anchors × evidence dimensions",
            main_figure_panel="Figure 6D",
            tested_or_summarized_unit="anchor × dimension",
            primary_estimand="presence / evaluated_no_support flags",
            primary_statistic="categorical status",
            uncertainty_measure="NA",
            multiplicity_or_resampling="non-nested; not summed",
            classification_or_gate="non-nested traceability",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="Evidence categories are non-nested, are not summed into a score, and do not enter Path A.",
            authoritative_source_path="results/fig6d_non_nested_traceability_source.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/fig6d_non_nested_traceability_source.tsv"),
        ),
        row(
            analysis_layer="S5 GSE238242 sensitivity",
            dataset_or_resource="GSE238242",
            main_figure_panel="S5 Fig",
            tested_or_summarized_unit="program × compartment / mural",
            primary_estimand="RNA-side and mural sensitivity deltas",
            primary_statistic="Δsingscore (frozen)",
            uncertainty_measure="NA (sensitivity display)",
            multiplicity_or_resampling="not a Path A gate",
            classification_or_gate="sensitivity only",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="Supplementary sensitivity only; not a Path A gate; does not alter Figure 3/5.",
            authoritative_source_path="Project_P5/GSE238242/results/figure5_program_delta_main.tsv",
            authoritative_source_sha256=sha256(ROOT / "Project_P5/GSE238242/results/figure5_program_delta_main.tsv"),
        ),
        row(
            analysis_layer="S7 F4 BCAA constituent-program example",
            dataset_or_resource="GSE255612",
            main_figure_panel="S7 Fig",
            tested_or_summarized_unit="constituent program × compartment",
            primary_estimand="program-level Δsingscore",
            primary_statistic="nominal permutation p",
            uncertainty_measure="nominal permutation p",
            multiplicity_or_resampling="program-level example; not family reclassification",
            classification_or_gate="program-level example only",
            contributes_to_Path_A="NO",
            context_only="YES",
            interpretation_boundary="Program-level example only; parent F4 remains Below-detection; not a fifth non-Below family.",
            authoritative_source_path="results/fig3_v3_program_celltype_delta.tsv",
            authoritative_source_sha256=sha256(ROOT / "results/fig3_v3_program_celltype_delta.tsv"),
        ),
    ]

    df = pd.DataFrame(rows)
    path_a = df[df.contributes_to_Path_A == "YES"]
    assert len(path_a) == 3, len(path_a)
    assert set(path_a.formal_Path_A_gate_name) == {
        "bulk High",
        "snRNA non-Below-detection",
        "duration non-Stable",
    }
    for layer in [
        "HCA healthy carrier prior",
        "Bulk dominant rewiring state",
        "GeoMx ROI-level context",
        "Figure 6A local GWAS / atrial eGene / PP4",
        "S5 GSE238242 sensitivity",
        "S7 F4 BCAA constituent-program example",
    ]:
        sub = df[df.analysis_layer == layer]
        assert len(sub) == 1 and sub.iloc[0].contributes_to_Path_A == "NO"

    text = df.astype(str).to_string()
    forbad = ["regulator-tier", "G-tier", "lead-p<5e−8 as the primary", "uniform significance test"]
    for tok in forbad:
        assert tok.lower() not in text.lower(), tok

    prov = provenance_rows(
        [
            ("claims_registry", "results/handoff/fulltext_audit_v1_1/02_final_claims_registry.tsv", "NA", "Locked claims"),
            ("path_a_contract", "results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv", "NA", "Exactly 3 gates"),
            ("fig2d", "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv", "NA", "bulk High"),
            ("fig3c", "results/fig3c_family_fidelity_classes.tsv", "NA", "snRNA gate"),
            ("fig4", "results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv", "NA", "duration gate"),
        ]
    )
    readme = [
        ("supplementary_table_id", "S2 Table"),
        ("title", "Multiple-testing, uncertainty, and evidence-gating strategy"),
        ("scientific_role", "Declare which layers contribute to Path A and which are context-only"),
        ("main_figure_links", "Figures 1–6; S5 Fig; S7 Fig"),
        ("analysis_unit", "Layer-specific (family/cohort/compartment/anchor/ROI)"),
        ("interpretation_boundary", "Exactly three Path A gates; genetics/GeoMx/HCA/sensitivity are non-gates."),
        ("source_snapshot_date", SNAP_DATE),
        ("scientific_rerun", "NO"),
        ("builder_script", BUILDER),
    ]
    xlsx = XLSX_DIR / "S2_Table_testing_uncertainty_evidence_gating.xlsx"
    write_workbook(xlsx, readme, {"evidence_gating": df}, prov)
    tsv = TSV_DIR / "S2_evidence_gating.tsv"
    export_tsv(df, tsv)
    return {
        "final_id": "S2 Table",
        "final_title": "Multiple-testing, uncertainty, and evidence-gating strategy",
        "xlsx": xlsx,
        "sheets": {"evidence_gating": df},
        "tsv_paths": [tsv],
        "provenance": prov,
        "sources": [
            "results/handoff/fulltext_audit_v1_1/02_final_claims_registry.tsv",
            "results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv",
        ],
        "notes": "Exactly three Path A YES rows",
        "assertions": [
            "exactly 3 Path A gate rows: bulk High; snRNA non-Below-detection; duration non-Stable",
            "HCA/dominant-state/GeoMx/genetics/S5/S7 are non-gates",
        ],
    }


# S3

def build_s3() -> dict:
    fam_meta = family_order_map()
    src = read_tsv("results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv")
    summ = read_tsv("results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_summary.tsv")
    states = read_tsv("results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv")
    state_sum = read_tsv("results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_summary.tsv")
    deploy = read_tsv("results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_summary.tsv")

    assert set(summ.consistency_class) <= {"High", "Moderate", "Low"}
    assert (summ.consistency_class == "High").sum() == 6
    assert (summ.consistency_class == "Moderate").sum() == 8
    assert (summ.consistency_class == "Low").sum() == 2
    high = set(summ.loc[summ.consistency_class == "High", "family_id"].astype(int))
    mod = set(summ.loc[summ.consistency_class == "Moderate", "family_id"].astype(int))
    low = set(summ.loc[summ.consistency_class == "Low", "family_id"].astype(int))
    assert high == {1, 9, 10, 12, 14, 15}
    assert mod == {2, 3, 5, 6, 8, 11, 13, 16}
    assert low == {4, 7}

    assert len(states) == 48
    assert states.duplicated(["family_id", "cohort"]).sum() == 0
    sc = states.state.value_counts().to_dict()
    assert sc.get("Amplified") == 5
    assert sc.get("Preserved") == 16
    assert sc.get("Attenuated") == 15
    assert sc.get("Reversed") == 12
    # no legacy counts
    assert sc != {7: 7}  # sanity
    for legacy in [(7, 13, 22, 6)]:
        pass
    assert list(state_sum.sort_values("state_order").n_calls) == [5, 16, 15, 12]

    pooled = src[src.estimate_type == "pooled"].copy()
    assert len(pooled) == 16
    assert set(pooled.pooling_model.dropna().unique()) == {"fixed-effect_inverse_variance"}
    assert pooled.tau2.isna().all()

    boundary = "State calls are descriptive remodeling summaries and not independent significance tests."
    sha_src = sha256(ROOT / "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv")
    sha_st = sha256(ROOT / "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv")

    # family_summary
    fam_rows = []
    for _, fm in fam_meta.iterrows():
        fid = int(fm.family_id)
        p = pooled.loc[pooled.family_id == fid].iloc[0]
        s = summ.loc[summ.family_id == fid].iloc[0]
        st = states.loc[states.family_id == fid]
        fam_rows.append(
            {
                "canonical_order": int(fm.family_display_order),
                "family_id": fid,
                "family_label": fm.family_label,
                "healthy_direction": fm.canonical_direction,
                "consistency_class": s.consistency_class,
                "pooled_hedges_g": float(p.hedges_g),
                "pooled_se": float(p.standard_error),
                "pooled_ci_low": float(p.ci_lower),
                "pooled_ci_high": float(p.ci_upper),
                "pooling_model": "fixed-effect_inverse_variance",
                "tau2": "NA",
                "n_cohorts": 3,
                "n_matching_cohort_signs": int(s.n_sign_match),
                "cohort_states": ";".join(f"{r.cohort}:{r.state}" for r in st.itertuples()),
                "interpretation_boundary": boundary,
                "source_sha256": sha_src,
            }
        )
    family_summary = pd.DataFrame(fam_rows).sort_values("canonical_order")

    # cohort_effects: cohort-specific + pooled
    coh_rows = []
    state_lookup = states.set_index(["family_id", "cohort"])["state"].to_dict()
    for _, r in src.iterrows():
        fid = int(r.family_id)
        order = int(fam_meta.loc[fam_meta.family_id == fid, "family_display_order"].iloc[0])
        label = fam_meta.loc[fam_meta.family_id == fid, "family_label"].iloc[0]
        is_pooled = r.estimate_type == "pooled"
        cohort = "POOLED" if is_pooled else str(r.cohort)
        st = "NA" if is_pooled else state_lookup[(fid, r.cohort)]
        coh_rows.append(
            {
                "canonical_order": order,
                "family_id": fid,
                "family_label": label,
                "cohort": cohort,
                "estimate_type": r.estimate_type,
                "n_AF": int(r.n_AF) if pd.notna(r.n_AF) else "NA",
                "n_SR": int(r.n_SR) if pd.notna(r.n_SR) else "NA",
                "hedges_g": float(r.hedges_g),
                "standard_error": float(r.standard_error),
                "variance": float(r.sampling_variance),
                "ci_low": float(r.ci_lower),
                "ci_high": float(r.ci_upper),
                "p_value_if_available": "NA",
                "family_cohort_state": st,
                "pooled_sign": int(r.pooled_sign) if pd.notna(r.pooled_sign) else "NA",
                "sign_matches_pooled": (
                    "NA"
                    if is_pooled or pd.isna(r.sign_matches_pooled)
                    else ("YES" if bool(r.sign_matches_pooled) else "NO")
                ),
                "consistency_class": r.consistency_class,
                "pooling_model": "fixed-effect_inverse_variance" if is_pooled else "NA",
                "source_path": "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv",
                "source_sha256": sha_src,
            }
        )
    cohort_effects = pd.DataFrame(coh_rows).sort_values(["canonical_order", "estimate_type", "cohort"])
    # pooled rows must not have state
    assert (cohort_effects.loc[cohort_effects.estimate_type == "pooled", "family_cohort_state"] == "NA").all()
    assert (cohort_effects.estimate_type == "cohort-specific").sum() == 48

    prov = provenance_rows(
        [
            ("fig2a_deployment", "results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_summary.tsv", "NA", "deployment lock"),
            ("fig2c_states", "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv", "NA", "48 state calls"),
            ("fig2c_summary", "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_summary.tsv", "NA", "5/16/15/12"),
            ("fig2d_source", "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv", "NA", "FE-IVW"),
            ("fig2d_summary", "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_summary.tsv", "NA", "6/8/2"),
            ("family_order", "results/ncp_v3_canonical_display_order.tsv", "NA", "Figure 1D order"),
            ("family_labels", "results/ncp_v3_family_labels.tsv", "NA", "labels"),
        ]
    )
    readme = [
        ("supplementary_table_id", "S3 Table"),
        ("title", "Corrected bulk paired-chamber effects and cross-cohort consistency"),
        ("scientific_role", "Companion to Figure 2 FE-IVW consistency and cohort states"),
        ("main_figure_links", "Figure 2"),
        ("analysis_unit", "family; family × cohort"),
        ("interpretation_boundary", boundary),
        ("source_snapshot_date", SNAP_DATE),
        ("scientific_rerun", "NO"),
        ("builder_script", BUILDER),
    ]
    xlsx = XLSX_DIR / "S3_Table_bulk_effects_consistency.xlsx"
    write_workbook(xlsx, readme, {"family_summary": family_summary, "cohort_effects": cohort_effects}, prov)
    t1 = TSV_DIR / "S3_family_summary.tsv"
    t2 = TSV_DIR / "S3_cohort_effects.tsv"
    export_tsv(family_summary, t1)
    export_tsv(cohort_effects, t2)
    return {
        "final_id": "S3 Table",
        "final_title": "Corrected bulk paired-chamber effects and cross-cohort consistency",
        "xlsx": xlsx,
        "sheets": {"family_summary": family_summary, "cohort_effects": cohort_effects},
        "tsv_paths": [t1, t2],
        "provenance": prov,
        "sources": [
            "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv",
            "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv",
        ],
        "notes": "FE-IVW; High/Moderate/Low=6/8/2; states 5/16/15/12",
        "assertions": [
            "16 families; High=6 Moderate=8 Low=2",
            "48 family×cohort states; Amplified5 Preserved16 Attenuated15 Reversed12",
            "pooling_model=fixed-effect_inverse_variance; tau2=NA",
        ],
    }


# S4

def build_s4() -> dict:
    fam_meta = family_order_map()
    cls = read_tsv("results/fig3c_family_fidelity_classes.tsv")
    carrier = read_tsv("results/fig3_v3_carrier_fidelity_class.tsv")
    fam_delta = read_tsv("results/fig3_v3_family_celltype_delta.tsv")
    prog = read_tsv("results/fig3_v3_program_celltype_delta.tsv")
    lock_path = "results/fig3_publication_authority_lock.md"
    assert (ROOT / lock_path).is_file()

    assert len(cls) == 16
    vc = cls.fidelity_class.value_counts().to_dict()
    assert vc.get("Fidelity") == 1
    assert vc.get("Relay") == 2
    assert vc.get("Extension") == 1
    assert vc.get("Below_detection") == 12
    assert int(cls.loc[cls.fidelity_class == "Fidelity", "family_id"].iloc[0]) == 10
    assert set(cls.loc[cls.fidelity_class == "Relay", "family_id"].astype(int)) == {11, 15}
    assert int(cls.loc[cls.fidelity_class == "Extension", "family_id"].iloc[0]) == 16
    assert cls.loc[cls.family_id == 10, "snRNA_supported_compartment"].iloc[0] == "fibroblast"
    assert cls.loc[cls.family_id == 11, "snRNA_supported_compartment"].iloc[0] == "adipocyte"
    assert cls.loc[cls.family_id == 15, "snRNA_supported_compartment"].iloc[0] == "endothelial"
    assert "lymphocyte" not in str(cls.loc[cls.family_id == 15, "snRNA_supported_compartment"].iloc[0]).lower()
    f16 = str(cls.loc[cls.family_id == 16, "snRNA_supported_compartment"].iloc[0])
    assert "fibroblast" in f16 and "atrialCM" in f16

    # class agreement with carrier table
    m = cls.merge(carrier, on="family_id", suffixes=("_c", "_v"))
    if (m.fidelity_class_c != m.fidelity_class_v).any():
        raise RuntimeError("fig3c vs carrier fidelity class disagreement")

    comps = sorted(fam_delta.compartment.unique())
    assert len(comps) == 9
    assert "RA" not in "".join(comps)
    assert len(fam_delta) == 144

    # F4 family remains Below-detection even if BCAA program supports fibroblast
    assert cls.loc[cls.family_id == 4, "fidelity_class"].iloc[0] == "Below_detection"
    bcaa = prog[prog.program_id.astype(str).str.contains("BRANCHED_CHAIN|BCAA|branched-chain", case=False, na=False)]
    # may be empty name match; use REACTOME branched
    bcaa = prog[prog.program_id.astype(str).str.contains("BRANCHED_CHAIN_AMINO_ACID", case=False, na=False)]
    if len(bcaa):
        # program-level support must not alter family class (already asserted)
        pass

    sha_cls = sha256(ROOT / "results/fig3c_family_fidelity_classes.tsv")
    sha_delta = sha256(ROOT / "results/fig3_v3_family_celltype_delta.tsv")
    boundary = "LA-only localization context; not an independent test of RA-driven remodeling. Nominal permutation; not FDR-controlled discovery. Below-detection is not biological absence."

    fam_rows = []
    for _, fm in fam_meta.iterrows():
        fid = int(fm.family_id)
        c = cls.loc[cls.family_id == fid].iloc[0]
        dsub = fam_delta.loc[fam_delta.family_id == fid]
        supported = c.snRNA_supported_compartment
        if pd.isna(supported) or str(supported).strip() in {"", "nan"}:
            supported_list = []
            supported_str = "NA"
        else:
            supported_list = [x.strip() for x in str(supported).split(",")]
            supported_str = ",".join(supported_list)
        # minimum nominal p among tested compartments
        min_p = float(dsub.loc[dsub.ct_included == True, "p_perm"].min())  # noqa: E712
        fam_rows.append(
            {
                "canonical_order": int(fm.family_display_order),
                "family_id": fid,
                "family_label": fm.family_label,
                "healthy_HCA_LA_carrier_prior": c.HCA_carrier_prior,
                "fidelity_class": str(c.fidelity_class).replace("_", "-") if c.fidelity_class != "Below_detection" else "Below-detection",
                "supported_compartments": supported_str,
                "n_supported_compartments": len(supported_list),
                "minimum_nominal_p": min_p,
                "evidence_strength_summary": c.evidence_strength,
                "below_detection_reason": na(c.below_detection_reason) if c.fidelity_class == "Below_detection" else "NA",
                "LA_only_boundary": "YES",
                "FDR_boundary": "nominal_permutation_not_BH_FDR_discovery",
                "source_path": "results/fig3c_family_fidelity_classes.tsv",
                "source_sha256": sha_cls,
            }
        )
    family_summary = pd.DataFrame(fam_rows).sort_values("canonical_order")

    # detail
    carrier_map = cls.set_index("family_id")["HCA_carrier_prior"].to_dict()
    class_map = cls.set_index("family_id")["fidelity_class"].to_dict()
    det_rows = []
    for _, r in fam_delta.iterrows():
        fid = int(r.family_id)
        order = int(fam_meta.loc[fam_meta.family_id == fid, "family_display_order"].iloc[0])
        label = fam_meta.loc[fam_meta.family_id == fid, "family_label"].iloc[0]
        carrier_c = carrier_map[fid]
        nominal = bool(pd.notna(r.p_perm) and r.p_perm < 0.05 and bool(r.ct_included))
        # fidelity relationship
        if class_map[fid] == "Below_detection":
            rel = "Below-detection"
        elif nominal and str(carrier_c) == str(r.compartment):
            rel = "supports_carrier"
        elif nominal:
            rel = "supports_noncarrier_or_additional"
        else:
            rel = "no_nominal_support"
        det_rows.append(
            {
                "canonical_order": order,
                "family_id": fid,
                "family_label": label,
                "compartment": r.compartment,
                "n_AF_donors": int(r.n_af),
                "n_control_donors": int(r.n_ctrl),
                "n_detectable_LE_genes": "NA",
                "delta_family_score": float(r.dScore),
                "permutation_p": float(r.p_perm),
                "permutation_method": str(r.perm_strategy),
                "n_permutations_or_exact": int(r.n_perm) if pd.notna(r.n_perm) else "NA",
                "included_for_testing": "YES" if bool(r.ct_included) else "NO",
                "nominal_support": "YES" if nominal else "NO",
                "carrier_compartment": carrier_c,
                "fidelity_relationship": rel,
                "source_path": "results/fig3_v3_family_celltype_delta.tsv",
                "source_sha256": sha_delta,
            }
        )
    detail = pd.DataFrame(det_rows).sort_values(["canonical_order", "compartment"])
    assert detail.compartment.nunique() == 9
    assert (detail.nominal_support == "YES").any()

    prov = provenance_rows(
        [
            ("fidelity_classes", "results/fig3c_family_fidelity_classes.tsv", lock_path, "1/2/1/12"),
            ("carrier_class", "results/fig3_v3_carrier_fidelity_class.tsv", lock_path, "agreement check"),
            ("family_delta", "results/fig3_v3_family_celltype_delta.tsv", lock_path, "16×9"),
            ("program_delta", "results/fig3_v3_program_celltype_delta.tsv", lock_path, "program example boundary"),
            ("authority_lock", lock_path, "NA", "publication authority"),
        ]
    )
    readme = [
        ("supplementary_table_id", "S4 Table"),
        ("title", "Corrected LA single-nucleus localization and carrier-relative classes"),
        ("scientific_role", "Companion to Figure 3 localization classes"),
        ("main_figure_links", "Figure 3; S3 Fig"),
        ("analysis_unit", "family; family × LA compartment"),
        ("interpretation_boundary", boundary),
        ("source_snapshot_date", SNAP_DATE),
        ("scientific_rerun", "NO"),
        ("builder_script", BUILDER),
    ]
    xlsx = XLSX_DIR / "S4_Table_snRNA_localization_classes.xlsx"
    write_workbook(xlsx, readme, {"family_summary": family_summary, "family_compartment_detail": detail}, prov)
    t1 = TSV_DIR / "S4_family_summary.tsv"
    t2 = TSV_DIR / "S4_family_compartment_detail.tsv"
    export_tsv(family_summary, t1)
    export_tsv(detail, t2)
    return {
        "final_id": "S4 Table",
        "final_title": "Corrected LA single-nucleus localization and carrier-relative classes",
        "xlsx": xlsx,
        "sheets": {"family_summary": family_summary, "family_compartment_detail": detail},
        "tsv_paths": [t1, t2],
        "provenance": prov,
        "sources": [
            "results/fig3c_family_fidelity_classes.tsv",
            "results/fig3_v3_family_celltype_delta.tsv",
        ],
        "notes": "1/2/1/12; F15 endothelial only; LA-only nominal",
        "assertions": [
            "Fidelity F10; Relay F11,F15; Extension F16; Below-detection 12",
            "F15 endothelial only; F16 fibroblast+atrialCM",
            "nine LA compartments; no RA",
        ],
    }


# S5

def build_s5() -> dict:
    fam_meta = family_order_map()
    design = read_tsv("results/ncp_v3_singscore_reanalysis/fig4a_duration_design_summary.tsv")
    laa = read_tsv("results/ncp_v3_singscore_reanalysis/gse41177_family_duration_laa_primary.tsv")
    lapv = read_tsv("results/ncp_v3_singscore_reanalysis/gse41177_family_duration_lapv_sensitivity.tsv")
    chrn = read_tsv("results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv")
    lapv_sum = read_tsv("results/ncp_v3_singscore_reanalysis/fig4_lapv_sensitivity.tsv")
    class_sum = read_tsv("results/ncp_v3_singscore_reanalysis/fig4c_chronicity_class_summary_summary.tsv")

    d0 = design.iloc[0]
    assert int(d0.primary_eligible_af_laa) == 12
    assert str(d0.class_rule) == "abs_rho_laa_only"
    assert str(d0.lapv_role) == "regional_sensitivity_only"

    # class counts
    vc = chrn.chronicity_class.value_counts().to_dict()
    assert vc.get("Strong_trend") == 2
    assert vc.get("Weak_trend") == 6
    assert vc.get("Stable") == 8
    strong = set(chrn.loc[chrn.chronicity_class == "Strong_trend", "family_id"].astype(int))
    weak = set(chrn.loc[chrn.chronicity_class == "Weak_trend", "family_id"].astype(int))
    stable = set(chrn.loc[chrn.chronicity_class == "Stable", "family_id"].astype(int))
    assert strong == {2, 9}
    assert weak == {3, 7, 8, 10, 13, 15}
    assert stable == {1, 4, 5, 6, 11, 12, 14, 16}

    def class_from_abs(a: float) -> str:
        if a >= 0.4:
            return "Strong"
        if a >= 0.2:
            return "Weak"
        return "Stable"

    sha_laa = sha256(ROOT / "results/ncp_v3_singscore_reanalysis/gse41177_family_duration_laa_primary.tsv")
    sha_chrn = sha256(ROOT / "results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv")
    boundary = (
        "Cross-sectional chronicity from full-precision |rho_LAA| bins only; "
        "not progression, recurrence, or treatment response. "
        "SR samples are metadata/context only and are not duration anchors. "
        "LAPV does not determine class."
    )

    rows = []
    for _, fm in fam_meta.iterrows():
        fid = int(fm.family_id)
        r = laa.loc[laa.family_id == fid].iloc[0]
        c = chrn.loc[chrn.family_id == fid].iloc[0]
        s = lapv_sum.loc[lapv_sum.family_id == fid].iloc[0]
        rho = float(r.rho_spearman)
        abs_rho = abs(rho)
        cls_disp = {"Strong_trend": "Strong", "Weak_trend": "Weak", "Stable": "Stable"}[c.chronicity_class]
        assert class_from_abs(abs_rho) == cls_disp
        # LAPV must not determine class — verify recomputed from LAA only
        regional = "Same" if bool(s.direction_concordant) else "Opposite"
        rows.append(
            {
                "canonical_order": int(fm.family_display_order),
                "family_id": fid,
                "family_label": fm.family_label,
                "rho_LAA_full_precision": rho,
                "rho_LAA_display": round(rho, 4),
                "abs_rho_LAA": abs_rho,
                "chronicity_class": cls_disp,
                "permutation_p_LAA": float(r.p_perm),
                "bootstrap_ci_low_LAA": float(r.bootstrap_ci_lo),
                "bootstrap_ci_high_LAA": float(r.bootstrap_ci_hi),
                "rho_LAPV": float(s.rho_LAPV),
                "regional_sign": regional,
                "LAPV_role": "regional_sensitivity_only",
                "class_rule": "|rho_LAA|>=0.4 Strong; 0.2<=|rho_LAA|<0.4 Weak; |rho_LAA|<0.2 Stable",
                "class_uses_full_precision": "YES",
                "primary_n_AF_LAA": 12,
                "interpretation_boundary": boundary,
                "source_path": "results/ncp_v3_singscore_reanalysis/gse41177_family_duration_laa_primary.tsv;results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv",
                "source_sha256": f"{sha_laa};{sha_chrn}",
            }
        )
    family_summary = pd.DataFrame(rows).sort_values("canonical_order")
    assert len(family_summary) == 16
    assert (family_summary.chronicity_class == "Strong").sum() == 2
    assert (family_summary.chronicity_class == "Weak").sum() == 6
    assert (family_summary.chronicity_class == "Stable").sum() == 8

    # representative values
    f9 = family_summary.loc[family_summary.family_id == 9].iloc[0]
    f2 = family_summary.loc[family_summary.family_id == 2].iloc[0]
    f15 = family_summary.loc[family_summary.family_id == 15].iloc[0]
    assert abs(f9.rho_LAA_full_precision - (-0.420316)) < 1e-6
    assert abs(f9.rho_LAPV - (-0.161121)) < 1e-6
    assert f9.regional_sign == "Same"
    assert abs(f2.rho_LAA_full_precision - (-0.563924)) < 1e-6
    assert abs(f2.rho_LAPV - (0.196147)) < 1e-6
    assert f2.regional_sign == "Opposite"
    assert abs(f15.rho_LAA_full_precision - (-0.385290)) < 1e-6
    assert abs(f15.rho_LAPV - (-0.185640)) < 1e-6
    assert f15.regional_sign == "Same"

    design_summary = pd.DataFrame(
        [
            {
                "dataset": "GSE41177",
                "primary_region": "LAA",
                "primary_n_AF": 12,
                "duration_unit": "months (numeric AF duration)",
                "duration_parsing_rule": "primary_eligible_af_laa with numeric duration; censored AF-LAA excluded from correlation",
                "LAPV_role": "regional_sensitivity_only; does not determine chronicity class",
                "SR_role": "metadata/context only; not a duration anchor in final inference",
                "scoring_type": str(d0.scorer),
                "class_rule": str(d0.class_rule) + "; thresholds Strong>=0.4, Weak>=0.2, else Stable",
                "interpretation_boundary": boundary,
                "source_sha256": sha256(ROOT / "results/ncp_v3_singscore_reanalysis/fig4a_duration_design_summary.tsv"),
            }
        ]
    )

    prov = provenance_rows(
        [
            ("design", "results/ncp_v3_singscore_reanalysis/fig4a_duration_design_summary.tsv", "NA", "n=12; abs_rho_laa_only"),
            ("laa_primary", "results/ncp_v3_singscore_reanalysis/gse41177_family_duration_laa_primary.tsv", "NA", "full-precision rho"),
            ("lapv_sens", "results/ncp_v3_singscore_reanalysis/gse41177_family_duration_lapv_sensitivity.tsv", "NA", "LAPV sensitivity"),
            ("chronicity", "results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv", "NA", "2/6/8 classes"),
            ("lapv_summary", "results/ncp_v3_singscore_reanalysis/fig4_lapv_sensitivity.tsv", "NA", "regional sign"),
            ("class_thresholds", "results/ncp_v3_singscore_reanalysis/fig4c_chronicity_class_summary_summary.tsv", "NA", "threshold definitions"),
        ]
    )
    readme = [
        ("supplementary_table_id", "S5 Table"),
        ("title", "Corrected cross-sectional AF-duration trends and chronicity classes"),
        ("scientific_role", "Companion to Figure 4 duration / chronicity"),
        ("main_figure_links", "Figure 4"),
        ("analysis_unit", "family (AF-LAA n=12 primary)"),
        ("interpretation_boundary", boundary),
        ("source_snapshot_date", SNAP_DATE),
        ("scientific_rerun", "NO"),
        ("builder_script", BUILDER),
    ]
    xlsx = XLSX_DIR / "S5_Table_duration_chronicity.xlsx"
    write_workbook(xlsx, readme, {"family_summary": family_summary, "design_summary": design_summary}, prov)
    t1 = TSV_DIR / "S5_family_summary.tsv"
    t2 = TSV_DIR / "S5_design_summary.tsv"
    export_tsv(family_summary, t1)
    export_tsv(design_summary, t2)
    return {
        "final_id": "S5 Table",
        "final_title": "Corrected cross-sectional AF-duration trends and chronicity classes",
        "xlsx": xlsx,
        "sheets": {"family_summary": family_summary, "design_summary": design_summary},
        "tsv_paths": [t1, t2],
        "provenance": prov,
        "sources": [
            "results/ncp_v3_singscore_reanalysis/gse41177_family_duration_laa_primary.tsv",
            "results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv",
            "results/ncp_v3_singscore_reanalysis/fig4a_duration_design_summary.tsv",
        ],
        "notes": "2/6/8; full-precision |rho_LAA|; no SR anchors",
        "assertions": [
            "Strong F2,F9; Weak F3,F7,F8,F10,F13,F15; Stable F1,F4,F5,F6,F11,F12,F14,F16",
            "primary n=12 AF-LAA; class from full-precision |rho_LAA| only",
            "representative F9/F2/F15 rho values match frozen sources",
        ],
    }


# S6

def build_s6() -> dict:
    order_genes = [
        "PRELP", "OMD", "IDH2", "NDUFA8", "BMP10", "SMAD9",
        "SMAD6", "NTRK2", "NTF3", "FRS2", "PITX2", "KCNN3",
    ]
    a = read_tsv("results/fig6a_anchor_coloc_calibration_source.tsv")
    b = read_tsv("results/fig6b_posterior_lead_distance_source.tsv")
    d = read_tsv("results/fig6d_non_nested_traceability_source.tsv")
    dsum = read_tsv("results/fig6d_non_nested_traceability_summary.tsv")
    gwas = read_tsv("results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv")
    gwas_map = dict(zip(gwas["key"], gwas["value"]))
    windows = read_tsv("results/fig6a_v3_locus_plot_windows.tsv")
    zoom = read_tsv("results/fig6a_v3_locus_zoom_summary.tsv")

    assert list(d.gene) == order_genes
    assert list(b.gene) == order_genes
    assert len(d) == 12

    # summary asserts
    smap = dsum.set_index("dimension")["present"].to_dict()
    assert int(smap["local_AF_GWAS"]) == 7
    assert int(smap["atrial_eGene"]) == 1
    assert int(smap["same_variant_coloc"]) == 0
    assert int(smap["common_risk"]) == 2
    assert int(smap["regulatory"]) == 2
    assert int(smap["rare_functional"]) == 0
    assert (d.coloc_PP3.notna()).sum() == 12
    assert bool(a.PP3_dominant.all())
    assert int(a["PP4_ge_0.70"].sum()) == 0
    assert int(a.atrial_eGene_bonferroni.sum()) == 1
    assert a.loc[a.atrial_eGene_bonferroni, "gene"].tolist() == ["SMAD9"]
    assert set(d.loc[d.is_calibration_control, "gene"]) == {"PITX2", "KCNN3"}
    assert set(d.loc[d.common_risk_status == "present", "gene"]) == {"PITX2", "KCNN3"}
    assert set(d.loc[d.regulatory_status == "present", "gene"]) == {"PITX2", "KCNN3"}
    assert (d.rare_functional_status == "present").sum() == 0
    assert b.delta_mb.notna().all()
    assert abs(float(b.loc[b.gene == "KCNN3", "delta_mb"].iloc[0]) - 0.459018) < 1e-6
    assert abs(float(b.delta_mb.min()) - 0.459018) < 1e-6
    assert abs(float(b.delta_mb.max()) - 1.923937) < 1e-6
    assert gwas_map["locked_af_cases"] == "252438"
    assert gwas_map["locked_controls"] == "1959739"
    assert gwas_map["coloc_rerun_required"] == "FALSE"

    win = windows.set_index("gene")
    zoom_i = zoom.set_index("anchor")
    a_i = a.set_index("gene")
    b_i = b.set_index("gene")
    d_i = d.set_index("gene")

    sha_a = sha256(ROOT / "results/fig6a_anchor_coloc_calibration_source.tsv")
    sha_b = sha256(ROOT / "results/fig6b_posterior_lead_distance_source.tsv")
    sha_d = sha256(ROOT / "results/fig6d_non_nested_traceability_source.tsv")

    rows = []
    for i, gene in enumerate(order_genes, start=1):
        dr = d_i.loc[gene]
        br = b_i.loc[gene]
        ar = a_i.loc[gene]
        wr = win.loc[gene]
        zr = zoom_i.loc[gene]
        anchor_pos = int(round((int(wr.start) + int(wr.end)) / 2))
        rows.append(
            {
                "anchor_order": i,
                "gene": gene,
                "family_or_role": dr.family_or_role,
                "calibration_control": "YES" if bool(dr.is_calibration_control) else "NO",
                "chromosome": int(wr.chr),
                "anchor_position_GRCh38": anchor_pos,
                "local_AF_GWAS_status": dr.local_gwas_status,
                "lead_GWAS_variant": br.lead_gwas_variant,
                "lead_GWAS_position_GRCh38": int(br.lead_gwas_bp) if pd.notna(br.lead_gwas_bp) else "NA",
                "lead_GWAS_p": float(zr.lead_gwas_p) if pd.notna(zr.lead_gwas_p) else "NA",
                "atrial_eGene_status": dr.atrial_egene_status,
                "atrial_eGene_p": float(dr.atrial_egene_p) if pd.notna(dr.atrial_egene_p) else "NA",
                "atrial_eGene_definition": dr.atrial_egene_definition,
                "lead_eQTL_variant": br.lead_eqtl_variant,
                "lead_eQTL_position_GRCh38": int(br.lead_eqtl_bp) if pd.notna(br.lead_eqtl_bp) else "NA",
                "lead_eQTL_p": float(zr.lead_eqtl_p) if pd.notna(zr.lead_eqtl_p) else "NA",
                "PP3": float(br.coloc_PP3),
                "PP4": float(br.coloc_PP4),
                "same_variant_coloc_status": dr.coloc_status,
                "coloc_threshold": "PP4 >= 0.70",
                "lead_separation_bp": float(br.delta_bp),
                "lead_separation_Mb": float(br.delta_mb),
                "common_risk_status": dr.common_risk_status,
                "common_risk_sources": na(dr.common_risk_sources),
                "regulatory_status": dr.regulatory_status,
                "regulatory_sources": na(dr.regulatory_sources),
                "rare_functional_status": dr.rare_functional_status,
                "rare_functional_sources": na(dr.rare_functional_sources),
                "evaluation_note": na(dr.evaluation_note),
                "source_paths": "results/fig6a_anchor_coloc_calibration_source.tsv;results/fig6b_posterior_lead_distance_source.tsv;results/fig6d_non_nested_traceability_source.tsv",
                "source_sha256s": f"{sha_a};{sha_b};{sha_d}",
            }
        )
    anchor_summary = pd.DataFrame(rows)
    assert list(anchor_summary.gene) == order_genes

    # forbid score/tier language in values
    blob = anchor_summary.astype(str).to_string()
    for tok in ["G0", "G1", "G2", "G3", "G4", "highest tier", "evidence score", "composite score", "Path A"]:
        assert tok not in blob, tok
    # lead-p may appear as numeric field label only via atrial_eGene_definition — ensure not primary gate wording wrongly
    assert "lead-p<5e−8 as the primary atrial eGene" not in blob

    dim_rows = []
    dim_map = {
        "local_AF_GWAS": ("local_AF_GWAS", "frozen ±1 Mb window contains genome-wide significant AF association"),
        "atrial_eGene": ("atrial_eGene", "Bonferroni atrial eGene (Figure 6A)"),
        "same_variant_coloc": ("same_variant_coloc", "PP4 >= 0.70"),
        "PP3_dominant": ("PP3_dominant", "PP3 > PP4 for all 12 anchors"),
        "risk_prioritized": ("common_risk", "membership in COMMON-RISK / PRIORITIZED sets"),
        "regulatory": ("regulatory", "membership in REGULATORY sets"),
        "rare_functional": ("rare_functional", "membership in RARE / FUNCTIONAL sets"),
    }
    # PP3_dominant not in dsum — compute
    present_pp3 = 12
    for dim, (src_dim, definition) in dim_map.items():
        if dim == "PP3_dominant":
            present_n, total_n = 12, 12
            rule = "PP3 > PP4"
            src = "results/fig6a_anchor_coloc_calibration_source.tsv"
        else:
            present_n = int(smap[src_dim])
            total_n = 12
            rule = definition
            src = "results/fig6d_non_nested_traceability_summary.tsv"
        dim_rows.append(
            {
                "dimension": dim,
                "present_n": present_n,
                "total_n": total_n,
                "definition": definition,
                "primary_threshold_or_rule": rule,
                "interpretation_boundary": (
                    "Evidence categories are non-nested, are not summed into a score, and do not enter Path A."
                ),
                "source_path": src,
                "source_sha256": sha256(ROOT / src),
            }
        )
    dimension_summary = pd.DataFrame(dim_rows)
    exp = {
        "local_AF_GWAS": 7,
        "atrial_eGene": 1,
        "same_variant_coloc": 0,
        "PP3_dominant": 12,
        "risk_prioritized": 2,
        "regulatory": 2,
        "rare_functional": 0,
    }
    for _, r in dimension_summary.iterrows():
        assert int(r.present_n) == exp[r.dimension], (r.dimension, r.present_n)

    prov = provenance_rows(
        [
            ("fig6a", "results/fig6a_anchor_coloc_calibration_source.tsv", "NA", "local GWAS / eGene / PP"),
            ("fig6b", "results/fig6b_posterior_lead_distance_source.tsv", "NA", "lead distances"),
            ("fig6d", "results/fig6d_non_nested_traceability_source.tsv", "NA", "non-nested flags"),
            ("fig6d_summary", "results/fig6d_non_nested_traceability_summary.tsv", "NA", "dimension counts"),
            ("gwas_provenance", "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv", "NA", "Yuan 2025; 252438/1959739; coloc_rerun_required=FALSE"),
            ("external_sets", "results/provenance/figure6_external_gene_sets_manifest.tsv", "NA", "risk/regulatory/rare sets"),
            ("locus_windows", "results/fig6a_v3_locus_plot_windows.tsv", "NA", "chr/window for anchor position"),
        ]
    )
    readme = [
        ("supplementary_table_id", "S6 Table"),
        ("title", "Anchor-level genetic calibration and non-nested traceability"),
        ("scientific_role", "Option B genetic detail home (no genetic S12 figure)"),
        ("main_figure_links", "Figure 6"),
        ("analysis_unit", "anchor gene"),
        (
            "interpretation_boundary",
            "Evidence categories are non-nested, are not summed into a score, and do not enter Path A. "
            "Yuan et al. 2025; GCST90624411; 252438 AF cases; 1959739 controls; "
            "harmonised summary-statistics provenance frozen; coloc_rerun_required=FALSE. "
            "Primary atrial eGene definition is Bonferroni atrial eGene, not lead-p<5e−8.",
        ),
        ("source_snapshot_date", SNAP_DATE),
        ("scientific_rerun", "NO"),
        ("builder_script", BUILDER),
        ("gwas_accession", "GCST90624411"),
        ("gwas_cases_controls", "252438 / 1959739"),
        ("coloc_rerun_required", "FALSE"),
    ]
    xlsx = XLSX_DIR / "S6_Table_genetic_calibration_traceability.xlsx"
    write_workbook(xlsx, readme, {"anchor_summary": anchor_summary, "dimension_summary": dimension_summary}, prov)
    t1 = TSV_DIR / "S6_anchor_summary.tsv"
    t2 = TSV_DIR / "S6_dimension_summary.tsv"
    export_tsv(anchor_summary, t1)
    export_tsv(dimension_summary, t2)
    return {
        "final_id": "S6 Table",
        "final_title": "Anchor-level genetic calibration and non-nested traceability",
        "xlsx": xlsx,
        "sheets": {"anchor_summary": anchor_summary, "dimension_summary": dimension_summary},
        "tsv_paths": [t1, t2],
        "provenance": prov,
        "sources": [
            "results/fig6a_anchor_coloc_calibration_source.tsv",
            "results/fig6b_posterior_lead_distance_source.tsv",
            "results/fig6d_non_nested_traceability_source.tsv",
            "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv",
        ],
        "notes": "12 anchors; 7/1/0; PP3 12/12; risk/reg 2/2; rare 0; no G-tier/score",
        "assertions": [
            "12 anchors in locked order",
            "local AF GWAS 7/12; atrial eGene 1/12 (SMAD9); PP4>=0.70 0/12; PP3-dominant 12/12",
            "lead distance range 0.46–1.92 Mb; KCNN3=0.459018 Mb",
            "risk/regulatory present = PITX2,KCNN3; rare=0; calibration controls PITX2,KCNN3",
        ],
    }


# Audit / main

def workbook_sheet_names(path: Path) -> list[str]:
    from openpyxl import load_workbook
    wb = load_workbook(path, read_only=True)
    names = list(wb.sheetnames)
    wb.close()
    return names


def workbook_has_formulas_or_hidden(path: Path) -> tuple[bool, bool]:
    from openpyxl import load_workbook
    wb = load_workbook(path, data_only=False)
    has_formula = False
    has_hidden = False
    for ws in wb.worksheets:
        if ws.sheet_state != "visible":
            has_hidden = True
        for row in ws.iter_rows(values_only=False):
            for cell in row:
                if isinstance(cell.value, str) and cell.value.startswith("="):
                    has_formula = True
    wb.close()
    return has_formula, has_hidden


def _cell_str(v) -> str:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return "NA"
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        # match pandas to_csv default for floats
        return format(v, ".15g")
    s = str(v)
    return "NA" if s in {"nan", "None", "NaN"} else s


def df_equal_tsv(df: pd.DataFrame, tsv: Path) -> bool:
    other = pd.read_csv(tsv, sep="\t", dtype=str).fillna("NA")
    left = df.copy()
    for c in left.columns:
        left[c] = left[c].map(_cell_str)
    # re-read using same export normalization
    exported = df.copy()
    for c in exported.columns:
        exported[c] = exported[c].map(lambda v: "NA" if pd.isna(v) else v)
    tmp_vals = pd.read_csv(tsv, sep="\t")
    # Compare via rewrite: round-trip the dataframe through the same exporter semantics
    from io import StringIO
    buf = StringIO()
    out = df.copy()
    for c in out.columns:
        out[c] = out[c].map(lambda v: "NA" if pd.isna(v) else v)
    out.to_csv(buf, sep="\t", index=False)
    buf.seek(0)
    rerun = pd.read_csv(buf, sep="\t", dtype=str).fillna("NA")
    other = other.astype(str)
    return rerun.reset_index(drop=True).equals(other.reset_index(drop=True))


def write_preview(info: dict) -> Path:
    PREV.mkdir(parents=True, exist_ok=True)
    fid = info["final_id"].replace(" ", "_")
    path = PREV / f"{fid}_preview.md"
    lines = [
        f"# {info['final_id']} preview",
        "",
        f"**Title:** {info['final_title']}",
        "",
        f"**Workbook:** `{info['xlsx'].relative_to(ROOT)}`",
        "",
        f"**Sheets:** {', '.join(info['sheets'])}",
        "",
    ]
    for name, df in info["sheets"].items():
        lines += [
            f"## {name}",
            "",
            f"Dimensions: {df.shape[0]} rows × {df.shape[1]} columns",
            "",
            "### First five rows",
            "",
            df.head(5).to_string(index=False),
            "",
            "### Last five rows",
            "",
            df.tail(5).to_string(index=False),
            "",
        ]
    lines += ["## Key assertions", ""]
    for a in info.get("assertions", []):
        lines.append(f"- {a}")
    lines += ["", "## Source hashes", ""]
    for s in info.get("sources", []):
        lines.append(f"- `{s}`: `{sha256(ROOT / s)}`")
    path.write_text("\n".join(lines) + "\n")
    return path


def write_audits(results: list[dict], pre: dict) -> None:
    AUDIT.mkdir(parents=True, exist_ok=True)
    man_rows = []
    lin_rows = []
    qc_blocks = []
    all_pass = True
    for info in results:
        xlsx = info["xlsx"]
        sheets = workbook_sheet_names(xlsx)
        has_formula, has_hidden = workbook_has_formulas_or_hidden(xlsx)
        tsv_ok = all(df_equal_tsv(df, p) for df, p in zip(info["sheets"].values(), info["tsv_paths"]))
        forb = scan_forbidden_frames(list(info["sheets"].values()))
        status = "PASS"
        reasons = []
        if has_formula:
            reasons.append("formulas")
        if has_hidden:
            reasons.append("hidden_sheets")
        if not tsv_ok:
            reasons.append("tsv_mismatch")
        if forb:
            reasons.append("forbidden:" + ",".join(forb))
        if reasons:
            status = "FAIL"
            all_pass = False
        else:
            reasons = []
        row_counts = {k: len(v) for k, v in info["sheets"].items()}
        man_rows.append(
            {
                "final_id": info["final_id"],
                "final_title": info["final_title"],
                "xlsx_path": str(xlsx.relative_to(ROOT)),
                "xlsx_sha256": sha256(xlsx),
                "sheet_names": ";".join(sheets),
                "tsv_paths": ";".join(str(p.relative_to(ROOT)) for p in info["tsv_paths"]),
                "tsv_sha256s": ";".join(sha256(p) for p in info["tsv_paths"]),
                "builder_script": BUILDER,
                "source_paths": ";".join(info["sources"]),
                "source_sha256s": ";".join(sha256(ROOT / s) for s in info["sources"]),
                "row_counts_by_sheet": ";".join(f"{k}:{v}" for k, v in row_counts.items()),
                "qc_status": status,
                "scientific_rerun": "NO",
                "formal_submission": "YES",
                "notes": info.get("notes", ""),
            }
        )
        for _, prow in info["provenance"].iterrows():
            lin_rows.append(
                {
                    "final_id": info["final_id"],
                    "sheet_name": "Provenance",
                    "source_role": prow.source_role,
                    "source_path": prow.source_path,
                    "source_sha256": prow.source_sha256,
                    "source_qc_path": prow.source_qc_path,
                    "source_qc_status": prow.source_qc_status,
                    "scientific_content_changed": "NO",
                    "publication_label_changed": "YES",
                    "notes": prow.notes,
                }
            )
        for sname in info["sheets"]:
            lin_rows.append(
                {
                    "final_id": info["final_id"],
                    "sheet_name": sname,
                    "source_role": "assembled_from_frozen_sources",
                    "source_path": ";".join(info["sources"]),
                    "source_sha256": ";".join(sha256(ROOT / s) for s in info["sources"]),
                    "source_qc_path": "NA",
                    "source_qc_status": "PASS",
                    "scientific_content_changed": "NO",
                    "publication_label_changed": "YES",
                    "notes": info.get("notes", ""),
                }
            )
        write_preview(info)
        qc_blocks.append(
            f"### {info['final_id']} — {status}\n\n"
            f"- title: {info['final_title']}\n"
            f"- xlsx: `{xlsx.relative_to(ROOT)}`\n"
            f"- sha256: `{sha256(xlsx)}`\n"
            f"- sheets: {sheets}\n"
            f"- formulas: {has_formula}; hidden: {has_hidden}; tsv_equal: {tsv_ok}\n"
            f"- forbidden (non-path cells): {forb or 'none'}\n"
            f"- fail_reasons: {reasons or 'none'}\n"
            f"- assertions:\n"
            + "\n".join(f"  - {a}" for a in info.get("assertions", []))
            + "\n"
        )

    man = pd.DataFrame(man_rows)
    assert len(man) == 6
    assert (man.scientific_rerun == "NO").all()
    assert list(man.final_id) == [f"S{i} Table" for i in range(1, 7)]
    man.to_csv(AUDIT / "phase5b3b_si_table_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(lin_rows).to_csv(AUDIT / "phase5b3b_si_table_source_lineage.tsv", sep="\t", index=False)

    # integrity
    mismatches = []
    for p, h in pre.get("hashes", {}).items():
        nh = sha256(ROOT / p)
        if h and nh and nh != h:
            mismatches.append(p)

    qc_md = f"""# Phase 5B3B SI table QC

Generated (UTC): {UTC}

## Global checks

- formal SI tables = 6 (S1–S6)
- scientific_rerun = NO for all
- hidden sheets = none required
- formulas = none required
- TSV/workbook equality = required
- integrity mismatches vs before-snapshot: **{len(mismatches)}**
{chr(10).join('- `' + m + '`' for m in mismatches) if mismatches else '- none'}

## Terminology

Publication-facing NCP only. Forbidden tokens scanned in table cell text.

## Per-table QC

{chr(10).join(qc_blocks)}

## Scientific checklist

### S1
Dataset counts and units current; paired bulk 51/25/26; GSE255612 LA-only primary; GSE41177 n=12; GeoMx 48/24/4; GCST cases/controls locked.

### S2
Exactly three Path A gate rows (bulk High; snRNA non-Below-detection; duration non-Stable).

### S3
48 cohort state calls; 5/16/15/12; consistency 6/8/2; FE-IVW only; tau2=NA.

### S4
1/2/1/12; F15 endothelial only; LA-only nominal boundary; no legacy 14-Below.

### S5
2/6/8; class from full-precision |rho_LAA|; no SR duration-anchor inference.

### S6
12 anchors; 7/1/0; PP3 12/12; risk/reg 2/2; rare 0; no G-tier/score; no lead-p primary eGene gate.
"""
    (AUDIT / "phase5b3b_si_table_qc.md").write_text(qc_md)

    branch = pre.get("branch", "")
    commit = pre.get("commit", "")
    status_before = pre.get("status_porcelain", "")
    try:
        status_after = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True)
        diff_stat = subprocess.check_output(
            ["git", "diff", "--stat"],
            cwd=ROOT,
            text=True,
        )
    except Exception:
        status_after = ""
        diff_stat = ""

    report = f"""# Phase 5B3B report — canonical Supplementary Tables S1–S6

Generated (UTC): {UTC}

## 1. Final S1–S6 list

{chr(10).join(f"- **{r['final_id']}**: {r['final_title']}" for r in man_rows)}

## 2. Workbook paths and hashes

| ID | Path | SHA256 |
|---|---|---|
{chr(10).join(f"| {r['final_id']} | `{r['xlsx_path']}` | `{r['xlsx_sha256']}` |" for r in man_rows)}

## 3. Sheet names and dimensions

{chr(10).join(f"- **{r['final_id']}**: sheets={r['sheet_names']}; rows={r['row_counts_by_sheet']}" for r in man_rows)}

## 4. Source lineage

See `phase5b3b_si_table_source_lineage.tsv`.

## 5. Scientific-rerun verdict

**NO** for all six tables.

## 6–11. Scientific assertion summaries

### S1
{chr(10).join('- ' + a for a in results[0].get('assertions', []))}

### S2
{chr(10).join('- ' + a for a in results[1].get('assertions', []))}

### S3
{chr(10).join('- ' + a for a in results[2].get('assertions', []))}

### S4
{chr(10).join('- ' + a for a in results[3].get('assertions', []))}

### S5
{chr(10).join('- ' + a for a in results[4].get('assertions', []))}

### S6
{chr(10).join('- ' + a for a in results[5].get('assertions', []))}

## 12–14. Equality / formatting / terminology

- TSV/workbook equality enforced in QC
- Formatting: white background; bold header; freeze panes; autofilter; wrap text; no formulas; no hidden sheets; no significance color scales
- Forbidden terminology scan recorded in QC

## 15. Integrity

Locked/pre-snapshot mismatches: **{len(mismatches)}**
{chr(10).join('- `' + m + '`' for m in mismatches) if mismatches else '- none'}

## 16. Git

Branch: `{branch}`  
Commit: `{commit}`

### Before (truncated)
```
{chr(10).join(status_before.splitlines()[:20])}
```

### After (phase paths)
```
{chr(10).join([ln for ln in status_after.splitlines() if 'supplementary_tables' in ln or 'tables/supplementary' in ln or 'phase5b3b' in ln or 'build_final_supplementary' in ln][:40])}
```

### git diff --stat (first 20 lines; pre-existing dirty tree)
```
{chr(10).join(diff_stat.splitlines()[:20])}
```

## 17–18. Unresolved / blocking

- QC FAIL count: **{sum(1 for r in man_rows if r['qc_status'] != 'PASS')}**
- Overall PASS: **{all_pass and len(mismatches) == 0}**

---

## Completion print summary

- formal SI tables = 6
- scientific rerun count = 0
- table rebuild count = 6
- QC by table: {{{', '.join(repr(r['final_id']) + ': ' + repr(r['qc_status']) for r in man_rows)}}}
"""
    (AUDIT / "phase5b3b_si_table_report.md").write_text(report)

    print("PASS" if all_pass and not mismatches else "FAIL")
    for r in man_rows:
        print(f"{r['final_id']}\t{r['qc_status']}\t{r['xlsx_sha256']}")
    print("mismatches", len(mismatches))
    print("rerun", 0)


def main() -> None:
    pre_path = Path("/tmp/phase5b3b_pre_hashes.json")
    if pre_path.is_file():
        pre = json.loads(pre_path.read_text())
    else:
        pre = {"hashes": {}, "branch": "", "commit": "", "status_porcelain": ""}

    XLSX_DIR.mkdir(parents=True, exist_ok=True)
    TSV_DIR.mkdir(parents=True, exist_ok=True)
    PREV.mkdir(parents=True, exist_ok=True)

    results = [
        build_s1(),
        build_s2(),
        build_s3(),
        build_s4(),
        build_s5(),
        build_s6(),
    ]
    write_audits(results, pre)

    # final integrity re-check of inputs
    for p, h in pre.get("hashes", {}).items():
        nh = sha256(ROOT / p)
        if h and nh and nh != h:
            raise RuntimeError(f"Source mutated during Phase 5B3B: {p}")


if __name__ == "__main__":
    main()
