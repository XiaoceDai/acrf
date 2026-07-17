#!/usr/bin/env python3
"""Fig 3 v3 per-gene Delta helper."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PB_DIR = ROOT / "Project_P3" / "GSE255612" / "pseudobulk_full"
LE_TSV = ROOT / "results" / "fig3_v3_program_le.tsv"
HCA_PB_DIR = ROOT / "Project_P3" / "HCA" / "pseudobulk_chamber"

OUT_TSV = ROOT / "results" / "fig3_v3_per_gene_delta.tsv"

CASE_PROGRAMS = [
    "REACTOME__REACTOME_SIGNALING_BY_NTRK2_TRKB",
    "REACTOME__REACTOME_KERATAN_SULFATE_BIOSYNTHESIS",
    "REACTOME__REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM",
]
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
N_MC = 5000
SEED = 20_260_418


def normalize_expr(expr):
    col_sums = expr.sum(axis=0).replace(0.0, np.nan)
    cpm = expr.div(col_sums, axis=1) * 1e6
    return np.log1p(cpm)


def per_sample_rank(expr_log):
    ranks = expr_log.rank(axis=0, method="average")
    return ranks / ranks.shape[0]


def perm_p(values, is_af, rng):
    n = len(values)
    n_af = int(is_af.sum())
    if n_af == 0 or n_af == n:
        return float("nan")
    obs = float(values[is_af].mean() - values[~is_af].mean())
    n_combo = math.comb(n, n_af)
    if n_combo <= 20_000:
        from itertools import combinations
        hits = 0
        for af_idx in combinations(range(n), n_af):
            mask = np.zeros(n, dtype=bool)
            mask[list(af_idx)] = True
            stat = float(values[mask].mean() - values[~mask].mean())
            if abs(stat) >= abs(obs) - 1e-12:
                hits += 1
        return hits / n_combo
    hits = 1
    for _ in range(N_MC):
        perm = rng.permutation(n)
        mask = np.zeros(n, dtype=bool)
        mask[perm[:n_af]] = True
        stat = float(values[mask].mean() - values[~mask].mean())
        if abs(stat) >= abs(obs) - 1e-12:
            hits += 1
    return hits / (N_MC + 1)


def hca_la_per_gene_logcpm(comp, donors_paired):
    expr_path = HCA_PB_DIR / comp / "expr.tsv"
    samp_path = HCA_PB_DIR / comp / "samples.tsv"
    if not expr_path.exists() or not samp_path.exists():
        return pd.Series(dtype=float)
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    samp = pd.read_csv(samp_path, sep="\t")
    keep = samp[(samp["region"] == "LA") & (samp["donor"].isin(donors_paired))]
    if keep.empty:
        return pd.Series(dtype=float)
    cols = [c for c in keep["sample_id"].tolist() if c in expr.columns]
    if not cols:
        return pd.Series(dtype=float)
    expr_log = normalize_expr(expr)
    return expr_log[cols].mean(axis=1)


def main():
    rng = np.random.default_rng(SEED)
    le_df = pd.read_csv(LE_TSV, sep="\t")
    le_lookup = {
        str(r["program_id"]): str(r["leading_edge_genes"]).split(",")
        for _, r in le_df.iterrows()
    }

    paired_path = HCA_PB_DIR / "paired_donors.tsv"
    paired = pd.read_csv(paired_path, sep="\t") if paired_path.exists() else pd.DataFrame()
    paired_donors_per_comp = {
        str(r["compartment"]): str(r["paired_donors"]).split(",")
        for _, r in paired.iterrows()
    }

    rows = []
    for program in CASE_PROGRAMS:
        le_genes = [g for g in le_lookup.get(program, []) if g]
        if not le_genes:
            print(f"[SKIP] no LE for {program}")
            continue
        for comp in COMPARTMENTS:
            expr_path = PB_DIR / comp / "expr.tsv"
            samp_path = PB_DIR / comp / "samples.tsv"
            if not (expr_path.exists() and samp_path.exists()):
                continue
            expr_raw = pd.read_csv(expr_path, sep="\t", index_col=0)
            meta = pd.read_csv(samp_path, sep="\t")
            expr_log = normalize_expr(expr_raw)
            ranks = per_sample_rank(expr_log)

            hca_baseline = (
                hca_la_per_gene_logcpm(comp, paired_donors_per_comp.get(comp, []))
                if comp in paired_donors_per_comp
                else pd.Series(dtype=float)
            )

            ctrl_donors = meta[meta["rhythm"] == "Ctrl"]["donor_id"].tolist()
            af_donors = meta[meta["rhythm"] == "AF"]["donor_id"].tolist()

            for gene in le_genes:
                if gene not in ranks.index:
                    continue
                row = ranks.loc[gene]
                ctrl_vals = row.reindex(ctrl_donors).dropna().to_numpy()
                af_vals = row.reindex(af_donors).dropna().to_numpy()
                if len(ctrl_vals) < 4 or len(af_vals) < 4:
                    continue
                d_rank = float(af_vals.mean() - ctrl_vals.mean())
                values = np.concatenate([ctrl_vals, af_vals])
                is_af = np.concatenate(
                    [np.zeros(len(ctrl_vals), dtype=bool),
                     np.ones(len(af_vals), dtype=bool)]
                )
                p = perm_p(values, is_af, rng)
                rows.append({
                    "program_id": program,
                    "gene": gene,
                    "compartment": comp,
                    "n_ctrl": int(len(ctrl_vals)),
                    "n_af": int(len(af_vals)),
                    "mean_norm_rank_ctrl": float(ctrl_vals.mean()),
                    "mean_norm_rank_af": float(af_vals.mean()),
                    "dRank": d_rank,
                    "p_perm": p,
                    "baseline_LA_logcpm": float(hca_baseline.get(gene, np.nan)),
                })

    out = pd.DataFrame(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.4f")
    print("wrote " + str(OUT_TSV.relative_to(ROOT)) + " (" + str(len(out)) + " rows)")


if __name__ == "__main__":
    main()
