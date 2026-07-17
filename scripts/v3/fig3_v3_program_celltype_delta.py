#!/usr/bin/env python3
"""Fig 3 v3 — program / family × compartment Δ-singscore on GSE255612 LA snRNA."""
from __future__ import annotations

import hashlib
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "ncp_v3_singscore_reanalysis"))

from singscore_exact import (  # noqa: E402
    TOLERANCE,
    exact_reference,
    legacy_v3_rank_singscore,
    normalize_gene_list,
    parse_leading_edge_field,
    scores_equivalent,
)

PANEL_TSV = ROOT / "results" / "innovation_ncp_v3_panel.tsv"
FAM_TSV = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
FROZEN_LE_TSV = (
    ROOT / "results" / "ncp_v3_singscore_reanalysis" / "frozen_program_metadata.tsv"
)
GMT_HALL = ROOT / "bulk" / "Project_P01" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC = ROOT / "bulk" / "Project_P01" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
DE_TSV = ROOT / "bulk" / "Project_P01" / "bulk_chamber_DE_full.tsv"
PB_DIR = ROOT / "Project_P3" / "GSE255612" / "pseudobulk_full"
R_EQ_SCRIPT = ROOT / "scripts" / "v3" / "validate_fig3_scorer_equivalence.R"

OUT_LE_AUDIT = ROOT / "results" / "fig3_v3_program_le_recomputed_audit.tsv"
OUT_LE_TSV = ROOT / "results" / "fig3_v3_program_le.tsv"  # frozen mirror for consumers
OUT_PROG_TSV = ROOT / "results" / "fig3_v3_program_celltype_delta.tsv"
OUT_FAM_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
OUT_FIG3B_TABLE_TSV = ROOT / "results" / "fig3b_family_compartment_dscore.tsv"
OUT_SCORES_TSV = ROOT / "results" / "fig3_v3_singscore_per_donor.tsv"
OUT_TEST_ORDER = ROOT / "results" / "fig3_v3_permutation_test_order.tsv"
OUT_EQ_SUMMARY = ROOT / "results" / "fig3_scorer_equivalence_summary.tsv"
OUT_EQ_DETAIL = ROOT / "results" / "fig3_scorer_equivalence_comparison.tsv"
OUT_EQ_QC = ROOT / "results" / "fig3_scorer_equivalence_qc.md"

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
EXPECTED_N_PROGRAMS = 30


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


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
    """Audit-only recomputation; not used for production scores."""
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


def load_frozen_programs() -> pd.DataFrame:
    require(FROZEN_LE_TSV.is_file(), f"Missing frozen LE manifest: {FROZEN_LE_TSV}")
    frozen = pd.read_csv(FROZEN_LE_TSV, sep="\t")
    require(len(frozen) == EXPECTED_N_PROGRAMS, f"Expected 30 frozen programs, got {len(frozen)}")
    panel = pd.read_csv(PANEL_TSV, sep="\t")
    fam = pd.read_csv(FAM_TSV, sep="\t")
    panel = panel.merge(
        fam[["pathway", "family_id", "family_source"]].drop_duplicates(),
        on="pathway",
        how="left",
    )
    require(len(panel) == EXPECTED_N_PROGRAMS, "Panel does not have 30 programs")
    panel_ids = panel["pathway"].astype(str).tolist()
    frozen_ids = frozen["program_id"].astype(str).tolist()
    require(set(panel_ids) == set(frozen_ids), "Frozen LE program IDs do not match panel")
    require(
        sorted(frozen_ids) == sorted(panel_ids),
        "Frozen LE / panel program-ID sets diverge",
    )

    rows = []
    for _, pr in panel.iterrows():
        pid = str(pr["pathway"])
        fr = frozen.loc[frozen["program_id"] == pid].iloc[0]
        genes = parse_leading_edge_field(fr["leading_edge_genes"])
        require(len(genes) >= 1, f"Empty frozen LE for {pid}")
        rows.append(
            {
                "program_id": pid,
                "program_order": int(pr["program_order"]),
                "direction": pr["direction"],
                "family_id": int(pr["family_id"]) if pd.notna(pr["family_id"]) else -1,
                "family_source": pr["family_source"] if pd.notna(pr["family_source"]) else "unknown",
                "n_LE_total": len(genes),
                "leading_edge_genes": ",".join(genes),
                "le_genes": genes,
            }
        )
    out = pd.DataFrame(rows).sort_values("program_order").reset_index(drop=True)
    require(len(out) == EXPECTED_N_PROGRAMS, "Frozen program table size changed")
    # Exact membership match against frozen file (normalized)
    for _, row in out.iterrows():
        fr = frozen.loc[frozen["program_id"] == row["program_id"]].iloc[0]
        frozen_set = set(parse_leading_edge_field(fr["leading_edge_genes"]))
        require(
            set(row["le_genes"]) == frozen_set,
            f"Normalized LE membership mismatch for {row['program_id']}",
        )
    return out


def write_recomputed_le_audit(programs: pd.DataFrame) -> None:
    """Retain prior bulk-DE/GSEA LE only as an audit mirror (not used for scoring)."""
    pathways: dict[str, list[str]] = {}
    pathways.update(load_gmt(GMT_HALL, "HALLMARK"))
    pathways.update(load_gmt(GMT_REAC, "REACTOME"))
    de = pd.read_csv(DE_TSV, sep="\t")
    de = de.dropna(subset=["gene_symbol", "t_stat"]).copy()
    de["abs_t"] = de["t_stat"].abs()
    de = de.sort_values("abs_t", ascending=False).drop_duplicates("gene_symbol")
    de = de.sort_values("t_stat", ascending=False)
    ranked_genes = de["gene_symbol"].tolist()
    ranked_stats = de["t_stat"].to_numpy(dtype=float)

    rows = []
    for _, r in programs.iterrows():
        pw = str(r["program_id"])
        full_set = set(pathways.get(pw, []))
        le = normalize_gene_list(gsea_leading_edge(ranked_genes, ranked_stats, full_set))
        frozen_set = set(r["le_genes"])
        rows.append(
            {
                "program_id": pw,
                "program_order": int(r["program_order"]),
                "direction": r["direction"],
                "family_id": int(r["family_id"]),
                "n_full_pathway": len(full_set),
                "n_LE_recomputed": len(le),
                "n_LE_frozen": len(frozen_set),
                "jaccard_vs_frozen": (
                    len(frozen_set & set(le)) / len(frozen_set | set(le))
                    if (frozen_set or le)
                    else float("nan")
                ),
                "exact_membership_match": set(le) == frozen_set,
                "leading_edge_genes_recomputed": ",".join(le),
                "leading_edge_genes_frozen": ",".join(r["le_genes"]),
                "used_for_scoring": False,
                "note": "audit_mirror_only; production scores use frozen LE",
            }
        )
    audit = pd.DataFrame(rows).sort_values("program_order").reset_index(drop=True)
    OUT_LE_AUDIT.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(OUT_LE_AUDIT, sep="\t", index=False)
    print(f"wrote {OUT_LE_AUDIT.relative_to(ROOT)} (audit mirror; not used for scoring)")


def deterministic_seed(base_seed: int, test_id: str) -> int:
    digest = hashlib.sha256(f"{base_seed}|{test_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % (2**32 - 1) + 1


def permutation_p(
    values: np.ndarray,
    is_af: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, str, int, float]:
    """Return (p, strategy, n_perm, monte_carlo_se)."""
    n = len(values)
    n_af = int(is_af.sum())
    if n_af == 0 or n_af == n:
        return float("nan"), "none", 0, float("nan")
    obs = float(values[is_af].mean() - values[~is_af].mean())
    n_combos = math.comb(n, n_af)
    if n_combos <= PERM_EXACT_LIMIT:
        from itertools import combinations

        hits = 0
        for af_idx in combinations(range(n), n_af):
            mask = np.zeros(n, dtype=bool)
            mask[list(af_idx)] = True
            stat = float(values[mask].mean() - values[~mask].mean())
            if abs(stat) >= abs(obs) - 1e-12:
                hits += 1
        p = hits / n_combos
        return p, "exact", n_combos, 0.0
    hits = 1
    for _ in range(N_MC_PERM):
        perm = rng.permutation(n)
        mask = np.zeros(n, dtype=bool)
        mask[perm[:n_af]] = True
        stat = float(values[mask].mean() - values[~mask].mean())
        if abs(stat) >= abs(obs) - 1e-12:
            hits += 1
    p = hits / (N_MC_PERM + 1)
    se = float(math.sqrt(max(p * (1.0 - p), 0.0) / N_MC_PERM))
    return p, "mc", N_MC_PERM, se


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
    """Display marker only; classification uses p < 0.05 exclusively."""
    if pd.isna(p):
        return ""
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "."
    return ""


def run_scorer_equivalence(programs: pd.DataFrame) -> dict:
    """Compare legacy v3, exact Python reference, and R singscore 1.28.1."""
    detail_rows: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="fig3-scorer-eq-") as tmp:
        tmp_path = Path(tmp)
        prog_export = programs[["program_id", "leading_edge_genes"]].copy()
        prog_path = tmp_path / "programs.tsv"
        prog_export.to_csv(prog_path, sep="\t", index=False)

        for comp in COMPARTMENTS:
            expr_path = PB_DIR / comp / "expr.tsv"
            samp_path = PB_DIR / comp / "samples.tsv"
            require(expr_path.is_file() and samp_path.is_file(), f"Missing PB for {comp}")
            expr_raw = pd.read_csv(expr_path, sep="\t", index_col=0)
            expr_log = normalize_expr(expr_raw)
            # Drop non-finite rows for ranking stability (same for all scorers).
            expr_log = expr_log.replace([np.inf, -np.inf], np.nan).dropna(how="any")
            require(expr_log.shape[0] > 0, f"Empty expression after cleaning: {comp}")

            matrix_path = tmp_path / f"{comp}_expr_log1p.tsv"
            expr_log.to_csv(matrix_path, sep="\t")
            r_out = tmp_path / f"{comp}_r_scores.tsv"
            completed = subprocess.run(
                ["Rscript", str(R_EQ_SCRIPT), str(matrix_path), str(prog_path), str(r_out)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            require(
                completed.returncode == 0 and r_out.is_file(),
                f"R singscore failed for {comp}:\n{completed.stdout}\n{completed.stderr}",
            )
            r_scores = pd.read_csv(r_out, sep="\t")

            for _, prog in programs.iterrows():
                genes = list(prog["le_genes"])
                legacy = legacy_v3_rank_singscore(expr_log, genes)
                exact = exact_reference(expr_log, genes)
                sub_r = r_scores.loc[r_scores["program_id"] == prog["program_id"]]
                require(len(sub_r) == len(expr_log.columns), f"R row count mismatch {comp}/{prog['program_id']}")
                r_series = pd.Series(
                    sub_r.set_index("sample_id")["r_singscore"],
                    index=expr_log.columns,
                    dtype=float,
                )
                for sample_id in expr_log.columns:
                    detail_rows.append(
                        {
                            "compartment": comp,
                            "sample_id": sample_id,
                            "program_id": prog["program_id"],
                            "legacy_v3_python": float(legacy.loc[sample_id])
                            if sample_id in legacy.index
                            else float("nan"),
                            "python_exact_reference": float(exact.loc[sample_id]),
                            "r_singscore_1_28_1": float(r_series.loc[sample_id]),
                        }
                    )

    detail = pd.DataFrame(detail_rows)
    detail["abs_diff_exact_vs_r"] = (
        detail["python_exact_reference"] - detail["r_singscore_1_28_1"]
    ).abs()
    detail["abs_diff_legacy_vs_r"] = (
        detail["legacy_v3_python"] - detail["r_singscore_1_28_1"]
    ).abs()
    detail["abs_diff_legacy_vs_exact"] = (
        detail["legacy_v3_python"] - detail["python_exact_reference"]
    ).abs()

    # Pairwise gates
    exact_vs_r = scores_equivalent(
        detail.set_index(["compartment", "sample_id", "program_id"])["python_exact_reference"],
        detail.set_index(["compartment", "sample_id", "program_id"])["r_singscore_1_28_1"],
    )[1]
    legacy_vs_exact = scores_equivalent(
        detail.set_index(["compartment", "sample_id", "program_id"])["legacy_v3_python"],
        detail.set_index(["compartment", "sample_id", "program_id"])["python_exact_reference"],
    )[1]
    legacy_vs_r = scores_equivalent(
        detail.set_index(["compartment", "sample_id", "program_id"])["legacy_v3_python"],
        detail.set_index(["compartment", "sample_id", "program_id"])["r_singscore_1_28_1"],
    )[1]

    exact_vs_r_passed = bool(exact_vs_r["passed"])
    legacy_equivalent = bool(legacy_vs_exact["passed"] and legacy_vs_r["passed"])
    # Production uses exact_reference iff it matches R.
    require(exact_vs_r_passed, "Python exact_reference is not equivalent to R singscore 1.28.1")

    summary = pd.DataFrame(
        [
            {
                "n_scores_compared": len(detail),
                "tolerance": TOLERANCE,
                "exact_vs_r_passed": exact_vs_r_passed,
                "max_abs_diff_exact_vs_r": exact_vs_r["max_abs_difference"],
                "legacy_vs_exact_passed": bool(legacy_vs_exact["passed"]),
                "max_abs_diff_legacy_vs_exact": legacy_vs_exact["max_abs_difference"],
                "legacy_vs_r_passed": bool(legacy_vs_r["passed"]),
                "max_abs_diff_legacy_vs_r": legacy_vs_r["max_abs_difference"],
                "legacy_v3_equivalent_to_canonical": legacy_equivalent,
                "production_scorer": "python_exact_reference_matching_singscore_1.28.1",
                "r_package_used_for_production_scores": False,
                "authoritative_formula": (
                    "rankGenes(tiesMethod=min)+simpleScore(centerScore=TRUE) "
                    "=> (mean_rank-(k+1)/2)/(N-k)-0.5; min measured genes = 3"
                ),
            }
        ]
    )
    OUT_EQ_DETAIL.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(OUT_EQ_DETAIL, sep="\t", index=False, float_format="%.15g")
    summary.to_csv(OUT_EQ_SUMMARY, sep="\t", index=False, float_format="%.15g")
    OUT_EQ_QC.write_text(
        "\n".join(
            [
                "# Figure 3 scorer equivalence QC",
                "",
                f"- Scores compared: {len(detail)}",
                f"- Tolerance: {TOLERANCE:g}",
                f"- Exact Python vs R 1.28.1: **{'PASS' if exact_vs_r_passed else 'FAIL'}** "
                f"(max |Δ| = {exact_vs_r['max_abs_difference']:.3g})",
                f"- Legacy v3 Python vs exact: **{'PASS' if legacy_vs_exact['passed'] else 'FAIL'}** "
                f"(max |Δ| = {legacy_vs_exact['max_abs_difference']:.3g})",
                f"- Legacy v3 Python vs R: **{'PASS' if legacy_vs_r['passed'] else 'FAIL'}** "
                f"(max |Δ| = {legacy_vs_r['max_abs_difference']:.3g})",
                f"- Production scorer: `{summary.loc[0, 'production_scorer']}`",
                "- R package is used for the equivalence reference only; production Fig 3 scores are Python exact_reference.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(
        f"Scorer equivalence: exact↔R={'PASS' if exact_vs_r_passed else 'FAIL'} "
        f"(max|{exact_vs_r['max_abs_difference']:.3g}|); "
        f"legacy↔canonical={'PASS' if legacy_equivalent else 'FAIL'}"
    )
    return {
        "exact_vs_r_passed": exact_vs_r_passed,
        "legacy_equivalent": legacy_equivalent,
        "max_abs_exact_vs_r": float(exact_vs_r["max_abs_difference"]),
        "max_abs_legacy_vs_exact": float(legacy_vs_exact["max_abs_difference"]),
        "summary": summary,
    }


def snapshot_prior_tables() -> None:
    for path in (OUT_PROG_TSV, OUT_FAM_TSV, OUT_FIG3B_TABLE_TSV, OUT_SCORES_TSV):
        if path.is_file():
            dest = path.with_name(path.stem + "_pre_correction" + path.suffix)
            if not dest.exists():
                shutil.copy2(path, dest)


def main() -> None:
    snapshot_prior_tables()
    programs = load_frozen_programs()
    # Production LE mirror for downstream consumers expecting fig3_v3_program_le.tsv
    programs.drop(columns=["le_genes"]).assign(
        n_full_pathway=np.nan,
        note="frozen_selection_time_LE_from_Figures_1_2",
    ).to_csv(OUT_LE_TSV, sep="\t", index=False)
    write_recomputed_le_audit(programs)

    eq = run_scorer_equivalence(programs)
    require(eq["exact_vs_r_passed"], "Cannot proceed without exact↔R equivalence")

    long_rows: list[dict] = []
    for comp in COMPARTMENTS:
        expr_path = PB_DIR / comp / "expr.tsv"
        samp_path = PB_DIR / comp / "samples.tsv"
        require(expr_path.is_file() and samp_path.is_file(), f"Missing PB for {comp}")
        expr_raw = pd.read_csv(expr_path, sep="\t", index_col=0)
        meta = pd.read_csv(samp_path, sep="\t")
        expr_log = normalize_expr(expr_raw)
        expr_log = expr_log.replace([np.inf, -np.inf], np.nan).dropna(how="any")
        meta_by_donor = meta.set_index("donor_id")

        for _, r in programs.iterrows():
            genes = list(r["le_genes"])
            scores = exact_reference(expr_log, genes)
            for col_name, score in scores.items():
                donor = str(col_name)
                if donor not in meta_by_donor.index:
                    continue
                rm = meta_by_donor.loc[donor]
                long_rows.append(
                    {
                        "compartment": comp,
                        "program_id": r["program_id"],
                        "program_order": int(r["program_order"]),
                        "direction": r["direction"],
                        "family_id": int(r["family_id"]),
                        "family_source": r["family_source"],
                        "donor_id": donor,
                        "rhythm": str(rm["rhythm"]),
                        "n_cells": int(rm["n_cells"]),
                        "score": float(score) if pd.notna(score) else float("nan"),
                        "scorer": "python_exact_reference_matching_singscore_1.28.1",
                        "le_source": "frozen_program_metadata",
                    }
                )

    scores_df = pd.DataFrame(long_rows)
    # Drop NA scores from aggregation inputs (ineligible measured sets)
    scores_df_valid = scores_df.loc[scores_df["score"].notna()].copy()
    OUT_SCORES_TSV.parent.mkdir(parents=True, exist_ok=True)
    scores_df.to_csv(OUT_SCORES_TSV, sep="\t", index=False, float_format="%.10g")
    print(f"wrote {OUT_SCORES_TSV.relative_to(ROOT)} ({len(scores_df)} rows)")

    test_order_rows: list[dict] = []

    # Program × compartment tests
    prog_keys = (
        scores_df_valid.groupby(["program_id", "compartment"], as_index=False)
        .size()
        .sort_values(["program_id", "compartment"])
        .reset_index(drop=True)
    )
    prog_rows = []
    for idx, key in prog_keys.iterrows():
        pid, comp = str(key["program_id"]), str(key["compartment"])
        g = scores_df_valid.loc[
            (scores_df_valid["program_id"] == pid) & (scores_df_valid["compartment"] == comp)
        ]
        ctrl = g.loc[g["rhythm"] == "Ctrl", "score"].to_numpy()
        af = g.loc[g["rhythm"] == "AF", "score"].to_numpy()
        n_ctrl, n_af = len(ctrl), len(af)
        prog_meta = g.iloc[0]
        test_id = f"program|{pid}|{comp}"
        seed = deterministic_seed(MC_SEED, test_id)
        test_order_rows.append(
            {
                "test_index": len(test_order_rows) + 1,
                "test_level": "program",
                "program_id": pid,
                "family_id": int(prog_meta["family_id"]),
                "compartment": comp,
                "test_id": test_id,
                "rng_seed": seed,
                "n_ctrl": n_ctrl,
                "n_af": n_af,
            }
        )
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
            "test_id": test_id,
            "rng_seed": seed,
        }
        if rec["ct_included"]:
            values = np.concatenate([ctrl, af])
            mask_af = np.concatenate(
                [np.zeros(n_ctrl, dtype=bool), np.ones(n_af, dtype=bool)]
            )
            rng = np.random.default_rng(seed)
            p, strat, n_perm, se = permutation_p(values, mask_af, rng)
        else:
            p, strat, n_perm, se = float("nan"), "none", 0, float("nan")
        rec["p_perm"] = p
        rec["perm_strategy"] = strat
        rec["n_perm"] = n_perm
        rec["p_perm_mc_se"] = se
        prog_rows.append(rec)

    prog_df = pd.DataFrame(prog_rows)
    prog_df["q_perm"] = float("nan")
    for comp in COMPARTMENTS:
        mask = (prog_df["compartment"] == comp) & prog_df["ct_included"]
        prog_df.loc[mask, "q_perm"] = bh(prog_df.loc[mask, "p_perm"].to_numpy())
    prog_df["sig_marker"] = prog_df["p_perm"].map(sig_marker)
    prog_df = prog_df.sort_values(["program_order", "compartment"]).reset_index(drop=True)
    prog_df.to_csv(OUT_PROG_TSV, sep="\t", index=False, float_format="%.10g")
    print(f"wrote {OUT_PROG_TSV.relative_to(ROOT)} ({len(prog_df)} rows)")

    # Family × compartment (program-first mean)
    fam_long = (
        scores_df_valid.groupby(
            ["family_id", "direction", "family_source", "compartment", "donor_id", "rhythm"],
            as_index=False,
        )["score"].mean()
    )
    fam_keys = (
        fam_long.groupby(["family_id", "compartment"], as_index=False)
        .size()
        .sort_values(["family_id", "compartment"])
        .reset_index(drop=True)
    )
    fam_rows = []
    for _, key in fam_keys.iterrows():
        fid, comp = int(key["family_id"]), str(key["compartment"])
        g = fam_long.loc[(fam_long["family_id"] == fid) & (fam_long["compartment"] == comp)]
        ctrl = g.loc[g["rhythm"] == "Ctrl", "score"].to_numpy()
        af = g.loc[g["rhythm"] == "AF", "score"].to_numpy()
        n_ctrl, n_af = len(ctrl), len(af)
        fmeta = g.iloc[0]
        test_id = f"family|{fid}|{comp}"
        seed = deterministic_seed(MC_SEED, test_id)
        test_order_rows.append(
            {
                "test_index": len(test_order_rows) + 1,
                "test_level": "family",
                "program_id": "",
                "family_id": fid,
                "compartment": comp,
                "test_id": test_id,
                "rng_seed": seed,
                "n_ctrl": n_ctrl,
                "n_af": n_af,
            }
        )
        rec = {
            "family_id": fid,
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
            "test_id": test_id,
            "rng_seed": seed,
        }
        if rec["ct_included"]:
            values = np.concatenate([ctrl, af])
            mask_af = np.concatenate(
                [np.zeros(n_ctrl, dtype=bool), np.ones(n_af, dtype=bool)]
            )
            rng = np.random.default_rng(seed)
            p, strat, n_perm, se = permutation_p(values, mask_af, rng)
        else:
            p, strat, n_perm, se = float("nan"), "none", 0, float("nan")
        rec["p_perm"] = p
        rec["perm_strategy"] = strat
        rec["n_perm"] = n_perm
        rec["p_perm_mc_se"] = se
        fam_rows.append(rec)

    fam_df = pd.DataFrame(fam_rows)
    fam_df["q_perm"] = float("nan")
    for comp in COMPARTMENTS:
        mask = (fam_df["compartment"] == comp) & fam_df["ct_included"]
        fam_df.loc[mask, "q_perm"] = bh(fam_df.loc[mask, "p_perm"].to_numpy())
    fam_df["sig_marker"] = fam_df["p_perm"].map(sig_marker)
    fam_df = fam_df.sort_values(["family_id", "compartment"]).reset_index(drop=True)
    fam_df.to_csv(OUT_FAM_TSV, sep="\t", index=False, float_format="%.10g")
    print(f"wrote {OUT_FAM_TSV.relative_to(ROOT)} ({len(fam_df)} rows)")

    pd.DataFrame(test_order_rows).to_csv(OUT_TEST_ORDER, sep="\t", index=False)
    print(f"wrote {OUT_TEST_ORDER.relative_to(ROOT)} ({len(test_order_rows)} tests)")

    # Fig 3B manuscript sidecar
    labels = pd.read_csv(ROOT / "results" / "ncp_v3_family_labels.tsv", sep="\t")[
        ["family_id", "family_label"]
    ].drop_duplicates()
    le_n = programs.groupby("family_id")["n_LE_total"].mean().rename("n_genes_detected")
    fig3b = fam_df.merge(labels, on="family_id", how="left").merge(
        le_n, on="family_id", how="left"
    )
    fig3b_out = pd.DataFrame(
        {
            "family_id": fig3b["family_id"],
            "family_label": fig3b["family_label"],
            "compartment": fig3b["compartment"],
            "dScore_AF_minus_ctrl": fig3b["dScore"],
            "p_perm": fig3b["p_perm"],
            "q_value": fig3b["q_perm"],
            "n_genes_detected": fig3b["n_genes_detected"],
            "n_AF": fig3b["n_af"],
            "n_ctrl": fig3b["n_ctrl"],
            "p_perm_mc_se": fig3b["p_perm_mc_se"],
            "perm_strategy": fig3b["perm_strategy"],
            "n_perm": fig3b["n_perm"],
            "sig_marker": fig3b["sig_marker"],
            "ct_included": fig3b["ct_included"],
        }
    )
    fig3b_out.to_csv(OUT_FIG3B_TABLE_TSV, sep="\t", index=False, float_format="%.10g")
    print(f"wrote {OUT_FIG3B_TABLE_TSV.relative_to(ROOT)}")

    print(
        "Fig 3 corrected scoring complete; "
        f"legacy_equivalent={eq['legacy_equivalent']}; "
        f"max|exact-R|={eq['max_abs_exact_vs_r']:.3g}"
    )


if __name__ == "__main__":
    main()
