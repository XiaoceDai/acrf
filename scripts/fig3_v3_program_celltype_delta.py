#!/usr/bin/env python3
"""
Fig 3 v3 — recompute program / family x compartment Delta-singscore on the
GSE255612 LA snRNA pseudobulk against the locked NCP_v3 panel
(30 programs / 16 families).

Algorithm (mirrors v2 spirit but updated to v3 panel):

    For each compartment c (9 reliable HCA-style cell types):
        expr (genes x samples, raw per-cell mean of counts)
            -> per-sample CPM normalization
            -> log1p
            -> "score(program, donor)" = mean log1p-CPM across program genes
        per donor scores partitioned into Ctrl (n=16) and AF (n=18)
        dScore = mean(AF) - mean(Ctrl)
        permutation p:  Monte Carlo, N = 10,000, seed = 20260418
        BH q across programs within compartment

    Family-level aggregation = mean of program-level (donor, score)
    then same Ctrl/AF dScore + permutation.

Outputs
-------
    results/fig3_v3_singscore_per_donor.tsv       long (program, donor, ct, score)
    results/fig3_v3_program_celltype_delta.tsv    program x ct stats
    results/fig3_v3_family_celltype_delta.tsv     family x ct stats
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

PANEL_TSV = ROOT / "results" / "innovation_ncp_v3_panel.tsv"
FAM_TSV = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
GMT_HALL = ROOT / "bulk" / "Project_P01" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC = ROOT / "bulk" / "Project_P01" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
DE_TSV = ROOT / "bulk" / "Project_P01" / "bulk_chamber_DE_full.tsv"
PB_DIR = ROOT / "Project_P3" / "GSE255612" / "pseudobulk_full"

OUT_LE_TSV = ROOT / "results" / "fig3_v3_program_le.tsv"
OUT_PROG_TSV = ROOT / "results" / "fig3_v3_program_celltype_delta.tsv"
OUT_FAM_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
OUT_SCORES_TSV = ROOT / "results" / "fig3_v3_singscore_per_donor.tsv"

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
MIN_DONORS = 4
N_MC_PERM = 10_000
PERM_EXACT_LIMIT = 20_000
MC_SEED = 20_260_418


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
    """Per-sample CPM + log1p (used only as input to rank-based singscore)."""
    col_sums = expr.sum(axis=0).replace(0.0, np.nan)
    cpm = expr.div(col_sums, axis=1) * 1e6
    return np.log1p(cpm)


def gsea_leading_edge(
    ranked_genes: list[str],
    ranked_stats: np.ndarray,
    gene_set: set[str],
    *,
    weighted: bool = True,
) -> list[str]:
    """Standard GSEA running-sum leading-edge extraction.

    Implements the same logic as fgsea / Subramanian 2005 with
    weight = |stats| when weighted=True (matches fgsea default `scoreType="std"`).

    Returns the list of gene-set members with rank position on the side of
    the running-sum peak (positive ES -> hits before peak; negative ES ->
    hits after trough).
    """
    n = len(ranked_genes)
    in_set = np.array([g in gene_set for g in ranked_genes], dtype=bool)
    n_h = int(in_set.sum())
    if n_h == 0 or n_h == n:
        return []
    weights = np.abs(ranked_stats) if weighted else np.ones(n)
    sum_hit = float(weights[in_set].sum())
    if sum_hit <= 0:
        return []
    inc_hit = np.where(in_set, weights / sum_hit, 0.0)
    inc_miss = np.where(in_set, 0.0, 1.0 / (n - n_h))
    running = np.cumsum(inc_hit - inc_miss)
    pos_peak = int(np.argmax(running))
    neg_peak = int(np.argmin(running))
    if abs(running[pos_peak]) >= abs(running[neg_peak]):
        peak = pos_peak
        return [ranked_genes[i] for i in range(peak + 1) if in_set[i]]
    peak = neg_peak
    return [ranked_genes[i] for i in range(peak, n) if in_set[i]]


def rank_singscore(expr_sample_genes: pd.DataFrame, gene_set: list[str]) -> pd.Series:
    """Rank-based singscore (Foroutan et al. 2018, simpleScore default).

    For each sample (column), genes are ranked by expression (average ties),
    normalized to (0, 1].  The score is the mean normalized rank of the
    gene-set members minus 0.5 -> range [-0.5, +0.5], centered at 0 = the
    expectation under the null hypothesis "gene set is a random subset of
    the genome".  Positive score = gene set sits in the upper rank half;
    negative = lower half.

    expr_sample_genes : (n_genes x n_samples) DataFrame, rows = gene symbols,
                        columns = sample IDs.
    gene_set          : list of gene symbols (subset of expr index).
    Returns           : pd.Series of scores indexed by sample ID.
    """
    if not gene_set:
        return pd.Series(dtype=float)
    # Rank within each sample (column).
    ranks = expr_sample_genes.rank(axis=0, method="average")
    n = ranks.shape[0]
    norm_ranks = ranks / n
    set_idx = [g for g in gene_set if g in norm_ranks.index]
    if len(set_idx) < 5:
        return pd.Series(dtype=float)
    return norm_ranks.loc[set_idx].mean(axis=0) - 0.5


def permutation_p(values: np.ndarray, is_af: np.ndarray, rng: np.random.Generator) -> tuple[float, str, int]:
    n = len(values)
    n_af = int(is_af.sum())
    if n_af == 0 or n_af == n:
        return float("nan"), "none", 0
    obs = float(values[is_af].mean() - values[~is_af].mean())
    n_combos = math.comb(n, n_af)
    if n_combos <= PERM_EXACT_LIMIT:
        from itertools import combinations
        hits = 0
        idx = np.arange(n)
        for af_idx in combinations(range(n), n_af):
            mask = np.zeros(n, dtype=bool)
            mask[list(af_idx)] = True
            stat = float(values[mask].mean() - values[~mask].mean())
            if abs(stat) >= abs(obs) - 1e-12:
                hits += 1
        return hits / n_combos, "exact", n_combos
    hits = 1
    for _ in range(N_MC_PERM):
        perm = rng.permutation(n)
        mask = np.zeros(n, dtype=bool)
        mask[perm[:n_af]] = True
        stat = float(values[mask].mean() - values[~mask].mean())
        if abs(stat) >= abs(obs) - 1e-12:
            hits += 1
    return hits / (N_MC_PERM + 1), "mc", N_MC_PERM


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    finite = np.isfinite(p)
    out = np.full_like(p, np.nan)
    if finite.sum() == 0:
        return out
    pf = p[finite]
    m = len(pf)
    order = np.argsort(pf)
    ranked = pf[order]
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    res = np.empty_like(q)
    res[order] = q
    out[finite] = np.clip(res, 0, 1)
    return out


def sig_marker(p: float) -> str:
    if pd.isna(p):
        return ""
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "."
    return ""


def main() -> None:
    panel = pd.read_csv(PANEL_TSV, sep="\t")
    fam = pd.read_csv(FAM_TSV, sep="\t")
    panel = panel.merge(
        fam[["pathway", "family_id", "family_source"]].drop_duplicates(),
        on="pathway", how="left",
    )

    pathways: dict[str, list[str]] = {}
    pathways.update(load_gmt(GMT_HALL, "HALLMARK"))
    pathways.update(load_gmt(GMT_REAC, "REACTOME"))

    # ---- Step A: leading-edge extraction on GSE112339 LA/RA bulk DE rank ----
    # This matches the v2 pipeline (r3b_ncpv2_singscore.R) so single-cell
    # singscore is computed on the focal LE subset rather than the full
    # pathway, which dramatically improves cell-type specificity (e.g. GAG
    # metabolism diseases LE = 14 ECM-proteoglycan genes vs 43 in the full
    # pathway).
    de = pd.read_csv(DE_TSV, sep="\t")
    de = de.dropna(subset=["gene_symbol", "t_stat"]).copy()
    de["abs_t"] = de["t_stat"].abs()
    de = de.sort_values("abs_t", ascending=False).drop_duplicates("gene_symbol")
    de = de.sort_values("t_stat", ascending=False)
    ranked_genes = de["gene_symbol"].tolist()
    ranked_stats = de["t_stat"].to_numpy(dtype=float)
    print(f"DE rank: {len(ranked_genes)} unique genes")

    le_per_program: dict[str, list[str]] = {}
    le_rows: list[dict] = []
    for _, r in panel.iterrows():
        pw = str(r["pathway"])
        full_set = set(pathways.get(pw, []))
        le = gsea_leading_edge(ranked_genes, ranked_stats, full_set)
        le_per_program[pw] = le
        le_rows.append({
            "program_id": pw,
            "program_order": int(r["program_order"]),
            "direction": r["direction"],
            "family_id": int(r["family_id"]) if pd.notna(r["family_id"]) else -1,
            "n_full_pathway": len(full_set),
            "n_LE_total": len(le),
            "leading_edge_genes": ",".join(le),
        })
    le_df = pd.DataFrame(le_rows).sort_values("program_order").reset_index(drop=True)
    OUT_LE_TSV.parent.mkdir(parents=True, exist_ok=True)
    le_df.to_csv(OUT_LE_TSV, sep="\t", index=False)
    print(f"wrote {OUT_LE_TSV.relative_to(ROOT)} ({len(le_df)} programs)")

    long_rows: list[dict] = []  # per (program, donor, ct, score)

    for comp in COMPARTMENTS:
        expr_path = PB_DIR / comp / "expr.tsv"
        samp_path = PB_DIR / comp / "samples.tsv"
        if not expr_path.exists() or not samp_path.exists():
            print(f"[SKIP] missing pseudobulk for {comp}")
            continue
        expr_raw = pd.read_csv(expr_path, sep="\t", index_col=0)
        meta = pd.read_csv(samp_path, sep="\t")
        expr_log = normalize_expr(expr_raw)

        # per donor x program: rank-based singscore on the program gene set
        # (Foroutan et al. 2018 simpleScore equivalent).  Operates on the
        # CPM+log1p matrix; rank conversion happens inside rank_singscore.
        meta_by_donor = meta.set_index("donor_id")
        for _, r in panel.iterrows():
            pw = str(r["pathway"])
            le_genes = le_per_program.get(pw, [])
            if len(le_genes) < 5:
                continue
            scores = rank_singscore(expr_log, le_genes)
            if scores.empty:
                continue
            for col_name, score in scores.items():
                donor = str(col_name)
                if donor not in meta_by_donor.index:
                    continue
                rm = meta_by_donor.loc[donor]
                long_rows.append({
                    "compartment": comp,
                    "program_id": pw,
                    "program_order": int(r["program_order"]),
                    "direction": r["direction"],
                    "family_id": int(r["family_id"]) if pd.notna(r["family_id"]) else -1,
                    "family_source": r["family_source"] if pd.notna(r["family_source"]) else "unknown",
                    "donor_id": donor,
                    "rhythm": str(rm["rhythm"]),
                    "n_cells": int(rm["n_cells"]),
                    "score": float(score),
                })

    scores_df = pd.DataFrame(long_rows)
    OUT_SCORES_TSV.parent.mkdir(parents=True, exist_ok=True)
    scores_df.to_csv(OUT_SCORES_TSV, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_SCORES_TSV.relative_to(ROOT)} ({len(scores_df)} rows)")

    rng = np.random.default_rng(MC_SEED)

    # ---- program × compartment ----
    prog_rows = []
    for (pid, comp), g in scores_df.groupby(["program_id", "compartment"]):
        ctrl = g.loc[g["rhythm"] == "Ctrl", "score"].to_numpy()
        af = g.loc[g["rhythm"] == "AF", "score"].to_numpy()
        n_ctrl, n_af = len(ctrl), len(af)
        prog_meta = g.iloc[0]
        rec = {
            "program_id": pid,
            "program_order": int(prog_meta["program_order"]),
            "direction": prog_meta["direction"],
            "family_id": int(prog_meta["family_id"]),
            "family_source": prog_meta["family_source"],
            "compartment": comp,
            "n_ctrl": n_ctrl,
            "n_af": n_af,
            "mean_ctrl": float(ctrl.mean()) if n_ctrl else float("nan"),
            "mean_af": float(af.mean()) if n_af else float("nan"),
            "sd_ctrl": float(ctrl.std(ddof=1)) if n_ctrl > 1 else float("nan"),
            "sd_af": float(af.std(ddof=1)) if n_af > 1 else float("nan"),
            "dScore": float(af.mean() - ctrl.mean()) if (n_ctrl and n_af) else float("nan"),
            "ct_included": (n_ctrl >= MIN_DONORS) and (n_af >= MIN_DONORS),
        }
        if rec["ct_included"]:
            values = np.concatenate([ctrl, af])
            mask_af = np.concatenate([np.zeros(n_ctrl, dtype=bool), np.ones(n_af, dtype=bool)])
            p, strat, n_perm = permutation_p(values, mask_af, rng)
        else:
            p, strat, n_perm = float("nan"), "none", 0
        rec["p_perm"] = p
        rec["perm_strategy"] = strat
        rec["n_perm"] = n_perm
        prog_rows.append(rec)

    prog_df = pd.DataFrame(prog_rows)
    # BH q within compartment
    prog_df["q_perm"] = float("nan")
    for comp in COMPARTMENTS:
        mask = (prog_df["compartment"] == comp) & prog_df["ct_included"]
        prog_df.loc[mask, "q_perm"] = bh(prog_df.loc[mask, "p_perm"].to_numpy())
    prog_df["sig_marker"] = prog_df["p_perm"].map(sig_marker)
    prog_df = prog_df.sort_values(["program_order", "compartment"]).reset_index(drop=True)
    prog_df.to_csv(OUT_PROG_TSV, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_PROG_TSV.relative_to(ROOT)} ({len(prog_df)} rows)")

    # ---- family × compartment ----
    fam_long = (
        scores_df.groupby(["family_id", "direction", "family_source", "compartment", "donor_id", "rhythm"], as_index=False)
        ["score"].mean()
    )
    fam_rows = []
    for (fid, comp), g in fam_long.groupby(["family_id", "compartment"]):
        ctrl = g.loc[g["rhythm"] == "Ctrl", "score"].to_numpy()
        af = g.loc[g["rhythm"] == "AF", "score"].to_numpy()
        n_ctrl, n_af = len(ctrl), len(af)
        fmeta = g.iloc[0]
        rec = {
            "family_id": int(fid),
            "direction": fmeta["direction"],
            "family_source": fmeta["family_source"],
            "compartment": comp,
            "n_ctrl": n_ctrl,
            "n_af": n_af,
            "mean_ctrl": float(ctrl.mean()) if n_ctrl else float("nan"),
            "mean_af": float(af.mean()) if n_af else float("nan"),
            "sd_ctrl": float(ctrl.std(ddof=1)) if n_ctrl > 1 else float("nan"),
            "sd_af": float(af.std(ddof=1)) if n_af > 1 else float("nan"),
            "dScore": float(af.mean() - ctrl.mean()) if (n_ctrl and n_af) else float("nan"),
            "ct_included": (n_ctrl >= MIN_DONORS) and (n_af >= MIN_DONORS),
        }
        if rec["ct_included"]:
            values = np.concatenate([ctrl, af])
            mask_af = np.concatenate([np.zeros(n_ctrl, dtype=bool), np.ones(n_af, dtype=bool)])
            p, strat, n_perm = permutation_p(values, mask_af, rng)
        else:
            p, strat, n_perm = float("nan"), "none", 0
        rec["p_perm"] = p
        rec["perm_strategy"] = strat
        rec["n_perm"] = n_perm
        fam_rows.append(rec)

    fam_df = pd.DataFrame(fam_rows)
    fam_df["q_perm"] = float("nan")
    for comp in COMPARTMENTS:
        mask = (fam_df["compartment"] == comp) & fam_df["ct_included"]
        fam_df.loc[mask, "q_perm"] = bh(fam_df.loc[mask, "p_perm"].to_numpy())
    fam_df["sig_marker"] = fam_df["p_perm"].map(sig_marker)
    fam_df = fam_df.sort_values(["family_id", "compartment"]).reset_index(drop=True)
    fam_df.to_csv(OUT_FAM_TSV, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {OUT_FAM_TSV.relative_to(ROOT)} ({len(fam_df)} rows)")


if __name__ == "__main__":
    main()
