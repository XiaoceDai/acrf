#!/usr/bin/env python3
"""
Fig 1F (v3, HCA) — chamber-resolved baseline cell-compartment attribution.

Source data: Project_P3/HCA/pseudobulk_chamber/<comp>/{expr.tsv,samples.tsv}
             (built by scripts/p_hca_pseudobulk_chamber.py from
              global_raw.h5ad / Litvinukova et al. 2020).

For every NCP_v3 program / family we compute, separately for LA and RA:
    (1) per-compartment mean expression averaged across paired donors,
    (2) row z-score across reliable compartments,
    (3) per-side dominant carrier (top z) and Mixed/Unresolved flag.

The two sides are then combined into a per-family carrier consensus:
    SHARED   - LA and RA both call the same dominant carrier
    LA_ONLY  - only LA passes the dominance test
    RA_ONLY  - only RA passes the dominance test
    SPLIT    - LA and RA each call a different dominant carrier
    MIXED    - neither side passes the dominance test

Outputs
-------
results/fig1f_v3_HCA_baseline_attribution.tsv
    per-program (and per-family) z scores for every (compartment, region).
results/fig1f_v3_HCA_family_carriers.tsv
    one row per family: LA_carrier, RA_carrier, consensus class,
    used by Fig 1F figure and Integrated Atlas.
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
PB_DIR = ROOT / "Project_P3" / "HCA" / "pseudobulk_chamber"

OUT_PROGRAM_TSV = ROOT / "results" / "fig1f_v3_HCA_baseline_attribution.tsv"
OUT_FAMILY_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"

COMPARTMENTS = [
    "atrialCM",
    "fibroblast",
    "endothelial",
    "macrophage",
    "lymphocyte",
    "pericyte",
    "vascular_smooth_muscle",
    "adipocyte",
    "neuronal",
]
REGIONS = ("LA", "RA")
MIN_PAIRED_DONORS = 4

CARRIER_MIN_Z = 1.0
CARRIER_MIN_GAP = 0.25
MIXED_LABEL = "Mixed/Unresolved"


def load_gmt(path: Path, prefix: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            out[f"{prefix}__{parts[0]}"] = [g for g in parts[2:] if g]
    return out


def normalize_expr(expr: pd.DataFrame) -> pd.DataFrame:
    """Per-sample CPM-like + log1p normalization to make compartment-level
    expression comparable.  HCA pseudobulk is the per-cell mean of raw counts,
    so adipocyte / atrialCM with high mRNA-per-nucleus would otherwise win
    every program by sheer library size."""
    col_sums = expr.sum(axis=0).replace(0.0, np.nan)
    cpm = expr.div(col_sums, axis=1) * 1e6
    return np.log1p(cpm)


def per_region_mean(
    expr_log: pd.DataFrame, meta: pd.DataFrame, donors: list[str], region: str
) -> pd.Series | None:
    sub = meta[(meta["region"] == region) & (meta["donor"].isin(donors))]
    if sub.empty:
        return None
    cols = [c for c in sub["sample_id"].tolist() if c in expr_log.columns]
    if not cols:
        return None
    return expr_log[cols].mean(axis=1)


def call_carrier(z_vec: pd.Series) -> tuple[str, float, float]:
    finite = z_vec[np.isfinite(z_vec.values)]
    if finite.empty:
        return MIXED_LABEL, np.nan, np.nan
    ranked = finite.sort_values(ascending=False)
    top = str(ranked.index[0])
    top_val = float(ranked.iloc[0])
    second_val = float(ranked.iloc[1]) if len(ranked) > 1 else -np.inf
    if (top_val < CARRIER_MIN_Z) or ((top_val - second_val) < CARRIER_MIN_GAP):
        return MIXED_LABEL, top_val, second_val
    return top, top_val, second_val


def consensus_class(la_carrier: str, ra_carrier: str) -> str:
    la_ok = la_carrier != MIXED_LABEL
    ra_ok = ra_carrier != MIXED_LABEL
    if la_ok and ra_ok and la_carrier == ra_carrier:
        return "SHARED"
    if la_ok and ra_ok and la_carrier != ra_carrier:
        return "SPLIT"
    if la_ok and not ra_ok:
        return "LA_ONLY"
    if ra_ok and not la_ok:
        return "RA_ONLY"
    return "MIXED"


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

    # Load per-compartment expression and samples; pre-compute per-(comp,region)
    # mean expression across PAIRED donors.
    comp_la_means: dict[str, pd.Series] = {}
    comp_ra_means: dict[str, pd.Series] = {}
    n_paired_per_comp: dict[str, int] = {}

    for comp in COMPARTMENTS:
        expr_path = PB_DIR / comp / "expr.tsv"
        samp_path = PB_DIR / comp / "samples.tsv"
        if not expr_path.exists() or not samp_path.exists():
            print(f"[WARN] missing {comp}; skipped")
            continue
        expr_raw = pd.read_csv(expr_path, sep="\t", index_col=0)
        expr = normalize_expr(expr_raw)
        meta = pd.read_csv(samp_path, sep="\t")
        donors_la = set(meta.loc[meta["region"] == "LA", "donor"])
        donors_ra = set(meta.loc[meta["region"] == "RA", "donor"])
        paired = sorted(donors_la & donors_ra)
        n_paired_per_comp[comp] = len(paired)
        if len(paired) < MIN_PAIRED_DONORS:
            print(f"[SKIP] {comp}: paired_donors={len(paired)} < {MIN_PAIRED_DONORS}")
            continue
        la_mean = per_region_mean(expr, meta, paired, "LA")
        ra_mean = per_region_mean(expr, meta, paired, "RA")
        if la_mean is not None:
            comp_la_means[comp] = la_mean
        if ra_mean is not None:
            comp_ra_means[comp] = ra_mean

    if not comp_la_means or not comp_ra_means:
        raise RuntimeError("no compartments passed paired-donor gate")

    la_df = pd.concat(comp_la_means, axis=1)
    ra_df = pd.concat(comp_ra_means, axis=1)
    la_df = la_df[[c for c in COMPARTMENTS if c in la_df.columns]]
    ra_df = ra_df[[c for c in COMPARTMENTS if c in ra_df.columns]]
    print(f"reliable compartments LA: {list(la_df.columns)}")
    print(f"reliable compartments RA: {list(ra_df.columns)}")

    program_rows: list[dict[str, object]] = []
    family_rows: dict[int, dict[str, object]] = {}

    for _, r in panel.iterrows():
        pw = str(r["pathway"])
        family_id = int(r["family_id"]) if pd.notna(r["family_id"]) else -1
        genes_all = pathways.get(pw, [])
        genes_la = [g for g in genes_all if g in la_df.index]
        genes_ra = [g for g in genes_all if g in ra_df.index]
        if len(genes_la) < 5 or len(genes_ra) < 5:
            continue
        la_mean = la_df.loc[genes_la].mean(axis=0)
        ra_mean = ra_df.loc[genes_ra].mean(axis=0)

        la_z = (la_mean - la_mean.mean()) / (la_mean.std() + 1e-9)
        ra_z = (ra_mean - ra_mean.mean()) / (ra_mean.std() + 1e-9)

        la_carrier, la_top, la_second = call_carrier(la_z)
        ra_carrier, ra_top, ra_second = call_carrier(ra_z)
        cls = consensus_class(la_carrier, ra_carrier)

        rec: dict[str, object] = {
            "program_order": int(r["program_order"]),
            "pathway": pw,
            "direction": r["direction"],
            "family_id": family_id,
            "family_source": r["family_source"] if pd.notna(r["family_source"]) else "unknown",
            "n_genes_used_LA": len(genes_la),
            "n_genes_used_RA": len(genes_ra),
            "LA_carrier": la_carrier,
            "RA_carrier": ra_carrier,
            "consensus_class": cls,
            "LA_top_z": round(la_top, 4) if np.isfinite(la_top) else np.nan,
            "RA_top_z": round(ra_top, 4) if np.isfinite(ra_top) else np.nan,
        }
        for comp in COMPARTMENTS:
            rec[f"LA_{comp}"] = float(la_z.get(comp, np.nan)) if comp in la_df.columns else np.nan
            rec[f"RA_{comp}"] = float(ra_z.get(comp, np.nan)) if comp in ra_df.columns else np.nan
        program_rows.append(rec)

        # Aggregate per family by averaging program-level z then re-classifying.
        fr = family_rows.setdefault(
            family_id,
            {
                "family_id": family_id,
                "direction": r["direction"],
                "members": [],
                "la_z_sum": pd.Series(0.0, index=la_df.columns),
                "ra_z_sum": pd.Series(0.0, index=ra_df.columns),
                "n": 0,
            },
        )
        fr["members"].append(pw)
        fr["la_z_sum"] = fr["la_z_sum"].add(la_z.reindex(la_df.columns).fillna(0.0), fill_value=0.0)
        fr["ra_z_sum"] = fr["ra_z_sum"].add(ra_z.reindex(ra_df.columns).fillna(0.0), fill_value=0.0)
        fr["n"] += 1

    program_df = pd.DataFrame(program_rows).sort_values("program_order").reset_index(drop=True)
    OUT_PROGRAM_TSV.parent.mkdir(parents=True, exist_ok=True)
    program_df.to_csv(OUT_PROGRAM_TSV, sep="\t", index=False, float_format="%.4f")
    print(f"wrote {OUT_PROGRAM_TSV.relative_to(ROOT)} ({len(program_df)} rows)")

    fam_records = []
    for fid, fr in family_rows.items():
        if fr["n"] == 0:
            continue
        la_mean_z = fr["la_z_sum"] / fr["n"]
        ra_mean_z = fr["ra_z_sum"] / fr["n"]
        la_carrier, la_top, _ = call_carrier(la_mean_z)
        ra_carrier, ra_top, _ = call_carrier(ra_mean_z)
        cls = consensus_class(la_carrier, ra_carrier)
        rec = {
            "family_id": fid,
            "direction": fr["direction"],
            "n_programs": fr["n"],
            "LA_carrier": la_carrier,
            "RA_carrier": ra_carrier,
            "consensus_class": cls,
            "LA_top_z": round(la_top, 4) if np.isfinite(la_top) else np.nan,
            "RA_top_z": round(ra_top, 4) if np.isfinite(ra_top) else np.nan,
        }
        for comp in COMPARTMENTS:
            rec[f"LA_{comp}"] = float(la_mean_z.get(comp, np.nan))
            rec[f"RA_{comp}"] = float(ra_mean_z.get(comp, np.nan))
        fam_records.append(rec)

    family_df = pd.DataFrame(fam_records).sort_values("family_id").reset_index(drop=True)
    OUT_FAMILY_TSV.parent.mkdir(parents=True, exist_ok=True)
    family_df.to_csv(OUT_FAMILY_TSV, sep="\t", index=False, float_format="%.4f")
    print(f"wrote {OUT_FAMILY_TSV.relative_to(ROOT)} ({len(family_df)} families)")


if __name__ == "__main__":
    main()
