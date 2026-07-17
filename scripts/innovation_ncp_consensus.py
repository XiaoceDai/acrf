#!/usr/bin/env python3
"""
Method innovation prototype:
Direction-aware consensus scoring + diversity-aware panel selection.

Motivation
----------
The current NCP selection uses hard gates (dual significance + 8/8 stability)
and then picks cluster representatives by |NES|. This script proposes a more
continuous evidence fusion strategy:

1) Consensus score (0-1)
   - fgsea significance
   - GSVA significance
   - fgsea effect size (|NES|)
   - GSVA effect size (|logFC|)
   - LOO stability fraction
   - directional coherence penalty when NES and GSVA disagree

2) Diversity-aware panel pruning
   - greedy selection by consensus score
   - reject pathways with high Jaccard overlap to already selected pathways
   - run per direction (LA/RA), then merge with minimum per-direction quota

Inputs
------
results/fig1c_pathway_selection_long.tsv
results/fig1c_ncp_panel_v2.tsv                    (optional, for comparison)
bulk/Project_P01/gene_sets/h.all.v2026.1.Hs.symbols.gmt
bulk/Project_P01/gene_sets/c2.cp.reactome.v2026.1.Hs.symbols.gmt

Outputs
-------
results/innovation_ncp_consensus_long.tsv
results/innovation_ncp_v3_panel.tsv
results/innovation_ncp_v3_summary.tsv
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

IN_LONG = ROOT / "results" / "fig1c_pathway_selection_long.tsv"
IN_V2 = ROOT / "results" / "fig1c_ncp_panel_v2.tsv"
GMT_HALL = ROOT / "bulk" / "Project_P01" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC = ROOT / "bulk" / "Project_P01" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"

OUT_LONG = ROOT / "results" / "innovation_ncp_consensus_long.tsv"
OUT_PANEL = ROOT / "results" / "innovation_ncp_v3_panel.tsv"
OUT_SUMMARY = ROOT / "results" / "innovation_ncp_v3_summary.tsv"

# Scoring hyper-parameters (explicitly written for method transparency)
FGSEA_Q_REF = 8.0
GSVA_Q_REF = 6.0
NES_REF = 3.0
LOGFC_REF = 0.60
STAB_MAX = 8.0

# Diversity + panel controls
MIN_STAB = 7
MIN_SCORE = 0.62
MAX_JACCARD_SIM = 0.50  # if similarity > 0.50, considered redundant
MAX_JACCARD_SIM_RA = 0.75  # RA arm is smaller; allow moderately higher overlap
MIN_PER_DIRECTION = 4
TARGET_PER_DIRECTION = 8


@dataclass
class ScoreWeights:
    fgsea_sig: float = 0.30
    gsva_sig: float = 0.20
    fgsea_eff: float = 0.20
    gsva_eff: float = 0.10
    stability: float = 0.20


def load_gmt(path: Path, prefix: str) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            key = f"{prefix}__{parts[0]}"
            out[key] = {g for g in parts[2:] if g}
    return out


def clip01(x: np.ndarray) -> np.ndarray:
    return np.clip(x, 0.0, 1.0)


def score_candidates(df: pd.DataFrame, w: ScoreWeights) -> pd.DataFrame:
    out = df.copy()

    padj = out["padj"].to_numpy(dtype=float)
    gsva_padj = out["gsva_padj"].to_numpy(dtype=float)
    nes = out["NES"].to_numpy(dtype=float)
    lfc = out["gsva_logFC"].to_numpy(dtype=float)
    stab = out["stability_consistent_runs"].to_numpy(dtype=float)

    fgsea_sig = clip01(-np.log10(np.clip(padj, 1e-300, 1.0)) / FGSEA_Q_REF)
    gsva_sig = clip01(-np.log10(np.clip(gsva_padj, 1e-300, 1.0)) / GSVA_Q_REF)
    fgsea_eff = clip01(np.abs(nes) / NES_REF)
    gsva_eff = clip01(np.abs(lfc) / LOGFC_REF)
    stab_frac = clip01(stab / STAB_MAX)

    direction_agree = np.sign(nes) == np.sign(lfc)
    coherence_penalty = np.where(direction_agree, 1.0, 0.35)

    base_score = (
        w.fgsea_sig * fgsea_sig
        + w.gsva_sig * gsva_sig
        + w.fgsea_eff * fgsea_eff
        + w.gsva_eff * gsva_eff
        + w.stability * stab_frac
    )
    consensus = base_score * coherence_penalty
    confidence_floor = np.minimum.reduce([fgsea_sig, gsva_sig, stab_frac]) * coherence_penalty

    out["e_fgsea_sig"] = fgsea_sig
    out["e_gsva_sig"] = gsva_sig
    out["e_fgsea_eff"] = fgsea_eff
    out["e_gsva_eff"] = gsva_eff
    out["e_stability"] = stab_frac
    out["direction_agree"] = direction_agree
    out["coherence_penalty"] = coherence_penalty
    out["consensus_score"] = consensus
    out["confidence_floor"] = confidence_floor

    out["innovation_tier"] = np.where(
        (out["consensus_score"] >= 0.78)
        & (out["stability_consistent_runs"] >= 8)
        & out["direction_agree"],
        "core",
        np.where(
            (out["consensus_score"] >= MIN_SCORE)
            & (out["stability_consistent_runs"] >= MIN_STAB)
            & out["direction_agree"],
            "supporting",
            "exploratory",
        ),
    )
    out = out.sort_values("consensus_score", ascending=False).reset_index(drop=True)
    return out


def jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 0.0
    u = a | b
    if not u:
        return 0.0
    return len(a & b) / len(u)


def select_diverse(
    sub: pd.DataFrame,
    gsets: Dict[str, Set[str]],
    max_sim: float,
) -> pd.DataFrame:
    selected_idx: List[int] = []
    for idx, row in sub.sort_values("consensus_score", ascending=False).iterrows():
        p = row["pathway"]
        gs = gsets.get(p, set())
        redundant = False
        for sidx in selected_idx:
            sp = sub.loc[sidx, "pathway"]
            sgs = gsets.get(sp, set())
            if jaccard(gs, sgs) > max_sim:
                redundant = True
                break
        if not redundant:
            selected_idx.append(idx)
    if not selected_idx:
        return sub.iloc[0:0].copy()
    return sub.loc[selected_idx].sort_values("consensus_score", ascending=False).copy()


def _direction_pool(scored: pd.DataFrame, direction: str) -> pd.DataFrame:
    d = scored[scored["direction"] == direction].copy()
    # Two-stage pool:
    #   strict   : high confidence core/supporting
    #   fallback : keep stability + coherence, relax significance/score slightly
    strict = d[
        d["direction_agree"]
        & (d["stability_consistent_runs"] >= MIN_STAB)
        & (d["consensus_score"] >= MIN_SCORE)
        & (d["padj"] < 0.05)
        & (d["gsva_padj"] < 0.05)
    ].copy()
    fallback = d[
        d["direction_agree"]
        & (d["stability_consistent_runs"] >= MIN_STAB)
        & (d["consensus_score"] >= 0.50)
        & (d["padj"] < 0.10)
        & (d["gsva_padj"] < 0.10)
    ].copy()
    strict["pool_tier"] = "strict"
    fallback["pool_tier"] = "fallback"
    both = pd.concat([strict, fallback], ignore_index=True)
    both = both.sort_values(["consensus_score", "stability_consistent_runs"], ascending=[False, False])
    both = both.drop_duplicates(subset=["pathway"], keep="first")
    return both.reset_index(drop=True)


def _fill_with_relaxed_similarity(
    current: pd.DataFrame,
    candidates: pd.DataFrame,
    gsets: Dict[str, Set[str]],
    target_k: int,
) -> pd.DataFrame:
    if len(current) >= target_k:
        return current.head(target_k).copy()
    # If strict diversity leaves empty slots, progressively relax redundancy cap.
    for max_sim in (0.60, 0.70, 0.80, 0.95):
        if len(current) >= target_k:
            break
        chosen = set(current["pathway"].tolist())
        remain = candidates[~candidates["pathway"].isin(chosen)].copy()
        add = select_diverse(remain, gsets, max_sim=max_sim)
        if add.empty:
            continue
        current = pd.concat([current, add], ignore_index=True)
        current = current.drop_duplicates(subset=["pathway"], keep="first")
        current = current.sort_values("consensus_score", ascending=False)
    return current.head(target_k).copy()


def pick_panel(scored: pd.DataFrame, gsets: Dict[str, Set[str]], target_k: int) -> pd.DataFrame:
    la_pool = _direction_pool(scored, "LA-biased")
    ra_pool = _direction_pool(scored, "RA-biased")

    la_sel = select_diverse(la_pool, gsets, MAX_JACCARD_SIM)
    ra_sel = select_diverse(ra_pool, gsets, MAX_JACCARD_SIM_RA)

    # Lock a minimum (safety) and a target (balanced panel narrative).
    min_la = min(MIN_PER_DIRECTION, len(la_sel))
    min_ra = min(MIN_PER_DIRECTION, len(ra_sel))
    tgt_la = min(TARGET_PER_DIRECTION, len(la_sel))
    tgt_ra = min(TARGET_PER_DIRECTION, len(ra_sel))

    locked_min = pd.concat([la_sel.head(min_la), ra_sel.head(min_ra)], ignore_index=False)
    locked_tgt = pd.concat([la_sel.head(tgt_la), ra_sel.head(tgt_ra)], ignore_index=False)
    locked = pd.concat([locked_min, locked_tgt], ignore_index=False)
    locked = locked.drop_duplicates(subset=["pathway"], keep="first")

    # If one direction is still below target after diversity pruning,
    # top-up from the direction pool by score (allowing higher redundancy).
    need_la = max(0, TARGET_PER_DIRECTION - int((locked["direction"] == "LA-biased").sum()))
    need_ra = max(0, TARGET_PER_DIRECTION - int((locked["direction"] == "RA-biased").sum()))
    if need_la > 0:
        extra_la = la_pool[~la_pool["pathway"].isin(set(locked["pathway"]))].head(need_la)
        locked = pd.concat([locked, extra_la], ignore_index=False)
    if need_ra > 0:
        extra_ra = ra_pool[~ra_pool["pathway"].isin(set(locked["pathway"]))].head(need_ra)
        locked = pd.concat([locked, extra_ra], ignore_index=False)
    locked = locked.drop_duplicates(subset=["pathway"], keep="first")

    locked_paths = set(locked["pathway"].tolist())
    rest = pd.concat([la_sel, ra_sel], ignore_index=False)
    rest = rest[~rest["pathway"].isin(locked_paths)]
    rest = rest.sort_values("consensus_score", ascending=False)

    locked = locked.sort_values("consensus_score", ascending=False)
    if len(locked) >= target_k:
        panel = locked.head(target_k).copy()
    else:
        n_left = target_k - len(locked)
        panel = pd.concat([locked, rest.head(n_left)], ignore_index=False)
        panel = panel.drop_duplicates(subset=["pathway"], keep="first")
        panel = panel.sort_values("consensus_score", ascending=False)

        # Global backfill with relaxed similarity if still short of target K.
        full_pool = pd.concat([la_pool, ra_pool], ignore_index=True)
        full_pool = full_pool.sort_values("consensus_score", ascending=False)
        panel = _fill_with_relaxed_similarity(panel, full_pool, gsets, target_k=target_k)

    panel = panel.sort_values(["direction", "consensus_score"], ascending=[True, False]).copy()
    panel["cluster_id"] = panel.groupby("direction").cumcount() + 1
    panel["cluster_id"] = np.where(
        panel["direction"].eq("LA-biased"),
        "LA_" + panel["cluster_id"].astype(str),
        "RA_" + panel["cluster_id"].astype(str),
    )
    panel = panel.sort_values("consensus_score", ascending=False).reset_index(drop=True)
    panel["program_order"] = np.arange(1, len(panel) + 1)
    return panel


def main() -> None:
    print("Loading pathway table ...")
    df = pd.read_csv(IN_LONG, sep="\t")

    required = {
        "pathway",
        "NES",
        "padj",
        "direction",
        "stability_consistent_runs",
        "gsva_logFC",
        "gsva_padj",
    }
    miss = required - set(df.columns)
    if miss:
        raise ValueError(f"Missing required columns in {IN_LONG.name}: {sorted(miss)}")

    print("Scoring candidates with consensus model ...")
    scored = score_candidates(df, ScoreWeights())
    scored.to_csv(OUT_LONG, sep="\t", index=False)
    print(f"wrote {OUT_LONG.relative_to(ROOT)} ({len(scored)} rows)")

    strict_cand = scored[
        scored["direction_agree"]
        & (scored["stability_consistent_runs"] >= MIN_STAB)
        & (scored["consensus_score"] >= MIN_SCORE)
        & (scored["padj"] < 0.05)
        & (scored["gsva_padj"] < 0.05)
    ].copy()
    print(f"Panel candidates after strict gate: {len(strict_cand)}")

    print("Loading gene sets for redundancy pruning ...")
    gsets: Dict[str, Set[str]] = {}
    gsets.update(load_gmt(GMT_HALL, "HALLMARK"))
    gsets.update(load_gmt(GMT_REAC, "REACTOME"))

    target_k = 30
    if IN_V2.exists():
        v2 = pd.read_csv(IN_V2, sep="\t")
        target_k = int(len(v2))
        v2_paths = set(v2["pathway"].tolist())
    else:
        v2 = None
        v2_paths = set()

    panel = pick_panel(scored, gsets, target_k=target_k)
    panel_cols = [
        "program_order",
        "cluster_id",
        "pathway",
        "direction",
        "consensus_score",
        "confidence_floor",
        "NES",
        "padj",
        "gsva_logFC",
        "gsva_padj",
        "stability_consistent_runs",
        "innovation_tier",
    ]
    panel[panel_cols].to_csv(OUT_PANEL, sep="\t", index=False)
    print(f"wrote {OUT_PANEL.relative_to(ROOT)} ({len(panel)} rows)")

    new_paths = set(panel["pathway"].tolist())
    overlap = len(new_paths & v2_paths) if v2_paths else np.nan
    jaccards = []
    panel_sorted = panel.sort_values("consensus_score", ascending=False)
    p_list = panel_sorted["pathway"].tolist()
    for i in range(len(p_list)):
        for j in range(i + 1, len(p_list)):
            a = gsets.get(p_list[i], set())
            b = gsets.get(p_list[j], set())
            jaccards.append(jaccard(a, b))
    mean_pairwise_jaccard = float(np.mean(jaccards)) if jaccards else np.nan

    summary = pd.DataFrame(
        [
            {
                "n_total_candidates": int(len(scored)),
                "n_quality_candidates": int(len(strict_cand)),
                "n_panel_v3": int(len(panel)),
                "n_panel_v2": int(len(v2_paths)) if v2_paths else np.nan,
                "n_overlap_with_v2": overlap,
                "overlap_fraction_vs_v2": (
                    round(overlap / len(v2_paths), 4) if v2_paths else np.nan
                ),
                "la_count_v3": int((panel["direction"] == "LA-biased").sum()),
                "ra_count_v3": int((panel["direction"] == "RA-biased").sum()),
                "median_consensus_v3": round(float(panel["consensus_score"].median()), 4),
                "mean_pairwise_jaccard_v3": (
                    round(mean_pairwise_jaccard, 4)
                    if not np.isnan(mean_pairwise_jaccard)
                    else np.nan
                ),
            }
        ]
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print("\nTop ten pathway programs by consensus score:")
    print(
        panel.sort_values("consensus_score", ascending=False)
        .head(10)[["program_order", "direction", "pathway", "consensus_score", "stability_consistent_runs"]]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
