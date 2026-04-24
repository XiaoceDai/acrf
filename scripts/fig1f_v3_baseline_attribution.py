#!/usr/bin/env python3
"""
Build v3 baseline program-by-cell-type attribution table for Fig1F.

This uses:
  - v3 locked panel pathways
  - healthy (Ctrl-only) pseudobulk expression from GSE255612
  - row z-score across reliable compartments (min Ctrl donors threshold)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

PANEL_TSV = ROOT / "results" / "innovation_ncp_v3_panel.tsv"
FAM_TSV = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
GMT_HALL = ROOT / "bulk" / "Project_P01" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC = ROOT / "bulk" / "Project_P01" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
PB_DIR = ROOT / "Project_P3" / "GSE255612" / "pseudobulk_full"

OUT_TSV = ROOT / "results" / "fig1f_v3_baseline_attribution.tsv"

COMPARTMENTS = [
    "atrialCM",
    "fibroblast",
    "endothelial",
    "lymphatic_endothelial",
    "macrophage",
    "lymphocyte",
    "mast_cell",
    "B_cell",
    "pericyte",
    "vascular_smooth_muscle",
    "adipocyte",
    "epicardial",
    "neuronal",
]
MIN_CTRL_DONORS = 5


def load_gmt(path: Path, prefix: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            out[f"{prefix}__{parts[0]}"] = [g for g in parts[2:] if g]
    return out


def main() -> None:
    for p in [PANEL_TSV, FAM_TSV, GMT_HALL, GMT_REAC, PB_DIR]:
        if not p.exists():
            raise FileNotFoundError(f"Missing input: {p}")

    panel = pd.read_csv(PANEL_TSV, sep="\t")
    fam = pd.read_csv(FAM_TSV, sep="\t")
    panel = panel.merge(
        fam[["pathway", "family_id", "family_source"]].drop_duplicates(),
        on="pathway",
        how="left",
    )

    pathways = {}
    pathways.update(load_gmt(GMT_HALL, "HALLMARK"))
    pathways.update(load_gmt(GMT_REAC, "REACTOME"))

    expr_means: dict[str, pd.Series] = {}
    n_ctrl: dict[str, int] = {}
    for comp in COMPARTMENTS:
        expr_path = PB_DIR / comp / "expr.tsv"
        samp_path = PB_DIR / comp / "samples.tsv"
        if not expr_path.exists() or not samp_path.exists():
            continue
        expr = pd.read_csv(expr_path, sep="\t", index_col=0)
        meta = pd.read_csv(samp_path, sep="\t")
        donors = meta.loc[meta["rhythm"] == "Ctrl", "donor_id"].astype(str).tolist()
        keep = [d for d in donors if d in expr.columns]
        n_ctrl[comp] = len(keep)
        if keep:
            expr_means[comp] = expr[keep].mean(axis=1)
        else:
            expr_means[comp] = pd.Series(np.nan, index=expr.index)

    means_df = pd.concat(expr_means, axis=1)
    means_df = means_df[[c for c in COMPARTMENTS if c in means_df.columns]]
    comp_cols = means_df.columns.tolist()
    reliable = [c for c in comp_cols if n_ctrl.get(c, 0) >= MIN_CTRL_DONORS]

    rows: list[dict[str, object]] = []
    for _, r in panel.iterrows():
        pw = str(r["pathway"])
        genes = [g for g in pathways.get(pw, []) if g in means_df.index]
        if len(genes) < 5:
            continue
        sub = means_df.loc[genes]
        mean_vec = sub.mean(axis=0)
        rel = mean_vec[reliable]
        mu = float(rel.mean())
        sd = float(rel.std()) + 1e-9

        rec: dict[str, object] = {
            "program_order": int(r["program_order"]),
            "pathway": pw,
            "direction": r["direction"],
            "family_id": int(r["family_id"]) if pd.notna(r["family_id"]) else -1,
            "family_source": r["family_source"] if pd.notna(r["family_source"]) else "unknown",
            "n_genes_used": int(len(genes)),
        }
        for c in comp_cols:
            rec[c] = float((mean_vec[c] - mu) / sd) if c in reliable else np.nan
        rows.append(rec)

    out = pd.DataFrame(rows).sort_values("program_order").reset_index(drop=True)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(out)} rows)")
    print(f"reliable compartments (n_ctrl >= {MIN_CTRL_DONORS}): {reliable}")


if __name__ == "__main__":
    main()
