#!/usr/bin/env python3
"""Score GSE41177 samples with python_exact_reference and frozen LE sets."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "ncp_v3_singscore_reanalysis"))
from singscore_exact import exact_reference, parse_leading_edge_field  # noqa: E402

EXPR_LAA = ROOT / "Project_P4" / "GSE41177" / "expr_symbol_LAA.tsv"
EXPR_LAPV = ROOT / "Project_P4" / "GSE41177" / "expr_symbol_LAPV.tsv"
META_TSV = ROOT / "Project_P4" / "GSE41177" / "sample_metadata.tsv"
LE_TSV = ROOT / "results" / "fig3_v3_program_le.tsv"
LE_ALT = ROOT / "results" / "ncp_v3_singscore_reanalysis" / "selection_time_leading_edges.tsv"
FAM_TSV = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
OUT_DIR = ROOT / "results" / "ncp_v3_singscore_reanalysis"
OUT_PROG = OUT_DIR / "gse41177_program_sample_scores.tsv"
OUT_FAM = OUT_DIR / "gse41177_family_sample_scores.tsv"
SCORER = "python_exact_reference"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def score_tissue(
    expr_path: Path,
    tissue: str,
    panel: pd.DataFrame,
    le_lookup: dict[str, list[str]],
    meta_gsms: list[str],
    *,
    le_sha: str,
    expr_sha: str,
) -> pd.DataFrame:
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    samples = [gsm for gsm in meta_gsms if gsm in expr.columns]
    expr = expr[samples]
    rows: list[dict] = []
    for _, row in panel.iterrows():
        program_id = str(row["pathway"])
        genes = le_lookup.get(program_id, [])
        scores = exact_reference(expr, genes)
        for gsm, score in scores.items():
            rows.append(
                {
                    "sample_id": gsm,
                    "region": tissue,
                    "program_id": program_id,
                    "program_order": int(row["program_order"]),
                    "direction": row["direction"],
                    "family_id": int(row["family_id"]),
                    "family_source": row["family_source"],
                    "score": float(score) if pd.notna(score) else np.nan,
                    "scorer": SCORER,
                    "le_manifest_sha256": le_sha,
                    "expr_matrix_sha256": expr_sha,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    for path in (EXPR_LAA, EXPR_LAPV, META_TSV, LE_TSV, FAM_TSV):
        assert path.is_file(), f"missing input: {path}"

    le_df = pd.read_csv(LE_TSV, sep="\t")
    if LE_ALT.is_file():
        alt = pd.read_csv(LE_ALT, sep="\t")
        assert set(le_df["program_id"]) == set(alt["program_id"]), "LE program_id mismatch"
    le_lookup = {
        str(row["program_id"]): parse_leading_edge_field(row["leading_edge_genes"])
        for _, row in le_df.iterrows()
    }
    le_sha = sha256(LE_TSV)

    panel = pd.read_csv(FAM_TSV, sep="\t")
    meta = pd.read_csv(META_TSV, sep="\t")
    laa_gsms = meta.loc[meta["tissue"] == "LAA", "gsm"].tolist()
    lapv_gsms = meta.loc[meta["tissue"] == "LAPV", "gsm"].tolist()

    laa_sha = sha256(EXPR_LAA)
    lapv_sha = sha256(EXPR_LAPV)
    prog_laa = score_tissue(
        EXPR_LAA, "LAA", panel, le_lookup, laa_gsms, le_sha=le_sha, expr_sha=laa_sha
    )
    prog_lapv = score_tissue(
        EXPR_LAPV, "LAPV", panel, le_lookup, lapv_gsms, le_sha=le_sha, expr_sha=lapv_sha
    )
    prog = pd.concat([prog_laa, prog_lapv], ignore_index=True)

    assert prog["program_id"].nunique() == 30, "expected 30 programs"
    assert prog["sample_id"].nunique() == 38, "expected 38 samples"
    assert len(prog) == 30 * 38, f"expected 1140 program rows, got {len(prog)}"
    assert prog["scorer"].eq(SCORER).all()

    eligible = prog.groupby("program_id")["score"].apply(lambda s: s.notna().all())
    assert eligible.all(), "NaN scores for programs with sufficient LE coverage"

    fam = (
        prog.groupby(
            ["family_id", "direction", "family_source", "sample_id", "region"],
            as_index=False,
        )
        .agg(
            score=("score", "mean"),
            n_programs=("program_id", "nunique"),
            scorer=("scorer", "first"),
            le_manifest_sha256=("le_manifest_sha256", "first"),
            expr_matrix_sha256=("expr_matrix_sha256", "first"),
        )
    )
    expected_programs = panel.groupby("family_id")["pathway"].nunique().to_dict()
    for fid, n_prog in expected_programs.items():
        sub = fam[fam["family_id"] == fid]
        assert (sub["n_programs"] == n_prog).all(), f"family {fid} program count mismatch"
    assert len(fam) == 16 * 38, f"expected 608 family rows, got {len(fam)}"

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    prog.to_csv(OUT_PROG, sep="\t", index=False, float_format="%.8f")
    fam.to_csv(OUT_FAM, sep="\t", index=False, float_format="%.8f")
    print(f"wrote {OUT_PROG.relative_to(ROOT)} ({len(prog)} rows)")
    print(f"wrote {OUT_FAM.relative_to(ROOT)} ({len(fam)} rows)")


if __name__ == "__main__":
    main()
