#!/usr/bin/env python3
"""
Robustness check for innovation NCP selection.

Goal
----
Quantify how stable each pathway's inclusion is under perturbations of:
  - evidence fusion weights
  - score threshold
  - stability threshold
  - significance thresholds

Inputs
------
results/innovation_ncp_consensus_long.tsv

Outputs
-------
results/innovation_ncp_bootstrap_frequency.tsv
results/innovation_ncp_bootstrap_summary.tsv
figures/innovation_ncp_bootstrap_frequency.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label
from plot_style import COL_DOUBLE, LA_COLOR, RA_COLOR, set_nature_style

IN_TSV = ROOT / "results" / "innovation_ncp_consensus_long.tsv"
OUT_FREQ = ROOT / "results" / "innovation_ncp_bootstrap_frequency.tsv"
OUT_SUM = ROOT / "results" / "innovation_ncp_bootstrap_summary.tsv"
OUT_PDF = ROOT / "figures" / "innovation_ncp_bootstrap_frequency.pdf"

N_RUNS = 400
TARGET_K = 30
TARGET_PER_DIRECTION = 8
RNG_SEED = 20260421

# Base weights follow innovation_ncp_consensus.py
BASE_W = np.array([0.30, 0.20, 0.20, 0.10, 0.20], dtype=float)


def sanitize_label(text: str) -> str:
    return (
        text.replace("²", "2")
        .replace("³", "3")
        .replace("⁺", "+")
        .replace("⁻", "-")
        .replace("β", "beta")
        .replace("α", "alpha")
    )


def direction_pool(
    s: np.ndarray,
    coh: np.ndarray,
    stab: np.ndarray,
    p_adj: np.ndarray,
    g_adj: np.ndarray,
    dir_arr: np.ndarray,
    direction: str,
    score_cut: float,
    padj_cut: float,
    gpadj_cut: float,
) -> np.ndarray:
    strict = (
        (dir_arr == direction)
        & (coh >= 0.99)
        & (stab >= 7)
        & (s >= score_cut)
        & (p_adj < padj_cut)
        & (g_adj < gpadj_cut)
    )
    fallback = (
        (dir_arr == direction)
        & (coh >= 0.99)
        & (stab >= 7)
        & (s >= 0.50)
        & (p_adj < 0.10)
        & (g_adj < 0.10)
    )
    idx = np.where(strict | fallback)[0]
    if len(idx) == 0:
        return idx
    return idx[np.argsort(-s[idx])]


def main() -> None:
    if not IN_TSV.exists():
        raise FileNotFoundError(
            f"Missing {IN_TSV}. Run scripts/innovation_ncp_consensus.py first."
        )
    df = pd.read_csv(IN_TSV, sep="\t")

    req = {
        "pathway",
        "direction",
        "coherence_penalty",
        "e_fgsea_sig",
        "e_gsva_sig",
        "e_fgsea_eff",
        "e_gsva_eff",
        "e_stability",
        "stability_consistent_runs",
        "padj",
        "gsva_padj",
        "consensus_score",
    }
    miss = req - set(df.columns)
    if miss:
        raise ValueError(f"Missing required columns: {sorted(miss)}")

    rng = np.random.default_rng(RNG_SEED)
    comp = df[["e_fgsea_sig", "e_gsva_sig", "e_fgsea_eff", "e_gsva_eff", "e_stability"]].to_numpy(float)
    coh = df["coherence_penalty"].to_numpy(float)

    path = df["pathway"].tolist()
    direction = df["direction"].tolist()
    dir_arr = np.array(direction, dtype=object)
    p_adj = df["padj"].to_numpy(float)
    g_adj = df["gsva_padj"].to_numpy(float)
    stab = df["stability_consistent_runs"].to_numpy(float)
    base = df["consensus_score"].to_numpy(float)

    n = len(df)
    counts = np.zeros(n, dtype=int)
    score_sum = np.zeros(n, dtype=float)
    size_runs = []
    ra_runs = []

    # Dirichlet concentration controls perturbation amplitude.
    # Higher => close to BASE_W; lower => wider perturbation.
    alpha = BASE_W * 60.0

    for _ in range(N_RUNS):
        w = rng.dirichlet(alpha)
        s = (comp @ w) * coh

        score_cut = float(rng.uniform(0.56, 0.68))
        padj_cut = float(rng.uniform(0.03, 0.08))
        gpadj_cut = float(rng.uniform(0.03, 0.08))

        la_pool = direction_pool(
            s, coh, stab, p_adj, g_adj, dir_arr, "LA-biased",
            score_cut=score_cut, padj_cut=padj_cut, gpadj_cut=gpadj_cut
        )
        ra_pool = direction_pool(
            s, coh, stab, p_adj, g_adj, dir_arr, "RA-biased",
            score_cut=score_cut, padj_cut=padj_cut, gpadj_cut=gpadj_cut
        )

        la_take = la_pool[: min(TARGET_PER_DIRECTION, len(la_pool))]
        ra_take = ra_pool[: min(TARGET_PER_DIRECTION, len(ra_pool))]
        chosen = list(la_take) + list(ra_take)

        if len(chosen) < TARGET_K:
            remain_mask = np.ones(n, dtype=bool)
            remain_mask[np.array(chosen, dtype=int)] = False
            remain_mask &= (coh >= 0.99) & (stab >= 7) & (s >= 0.45)
            remain_idx = np.where(remain_mask)[0]
            remain_idx = remain_idx[np.argsort(-s[remain_idx])]
            need = TARGET_K - len(chosen)
            chosen += remain_idx[:need].tolist()

        if not chosen:
            continue
        take = np.array(chosen[:TARGET_K], dtype=int)
        counts[take] += 1
        score_sum[take] += s[take]
        size_runs.append(len(take))
        ra_runs.append(int(np.sum(np.array([direction[i] for i in take]) == "RA-biased")))

    freq = counts / max(N_RUNS, 1)
    mean_sel_score = np.divide(score_sum, np.maximum(counts, 1))

    out = pd.DataFrame(
        {
            "pathway": path,
            "direction": direction,
            "inclusion_freq": np.round(freq, 4),
            "mean_score_when_selected": np.round(mean_sel_score, 4),
            "base_consensus_score": np.round(base, 4),
            "stability_consistent_runs": stab.astype(int),
            "padj": p_adj,
            "gsva_padj": g_adj,
        }
    ).sort_values(["inclusion_freq", "base_consensus_score"], ascending=[False, False])
    out["robust_tier"] = np.where(
        out["inclusion_freq"] >= 0.80,
        "very-stable",
        np.where(out["inclusion_freq"] >= 0.50, "stable", "sensitive"),
    )
    out.to_csv(OUT_FREQ, sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "n_runs": N_RUNS,
                "target_k": TARGET_K,
                "mean_selected_per_run": round(float(np.mean(size_runs)) if size_runs else 0.0, 3),
                "mean_ra_per_run": round(float(np.mean(ra_runs)) if ra_runs else 0.0, 3),
                "n_path_very_stable_freq_ge_0_8": int((out["inclusion_freq"] >= 0.80).sum()),
                "n_path_stable_freq_ge_0_5": int((out["inclusion_freq"] >= 0.50).sum()),
                "top_path_freq": round(float(out["inclusion_freq"].max()), 4),
                "median_path_freq": round(float(out["inclusion_freq"].median()), 4),
            }
        ]
    )
    summary.to_csv(OUT_SUM, sep="\t", index=False)

    # Plot top 36 by inclusion frequency
    set_nature_style()
    top = out.head(36).copy().iloc[::-1]
    y = np.arange(len(top))
    colors = [LA_COLOR if d == "LA-biased" else RA_COLOR for d in top["direction"]]

    fig, ax = plt.subplots(figsize=(COL_DOUBLE * 1.05, COL_DOUBLE * 0.52))
    ax.barh(y, top["inclusion_freq"], color=colors, edgecolor="white", linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [sanitize_label(short_pathway_label(p, 46)) for p in top["pathway"]],
        fontsize=6.6,
    )
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Inclusion frequency across perturbed runs")
    ax.set_title("Innovation robustness: pathway selection frequency", loc="left", fontsize=10, fontweight="bold")
    ax.axvline(0.5, color="#888888", linewidth=0.6, linestyle="--")
    ax.axvline(0.8, color="#666666", linewidth=0.6, linestyle="--")
    ax.text(0.5, -1.2, "stable threshold", ha="center", va="bottom", fontsize=6.5, color="#666")
    ax.text(0.8, -1.2, "very-stable threshold", ha="center", va="bottom", fontsize=6.5, color="#666")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.text(0.85, 0.08, "LA", color=LA_COLOR, fontsize=7, fontweight="bold")
    fig.text(0.89, 0.08, "/", color="#666666", fontsize=7)
    fig.text(0.90, 0.08, "RA", color=RA_COLOR, fontsize=7, fontweight="bold")
    fig.savefig(OUT_PDF, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    print(f"wrote {OUT_FREQ.relative_to(ROOT)} ({len(out)} rows)")
    print(f"wrote {OUT_SUM.relative_to(ROOT)}")
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")
    print("\nTop 12 by inclusion frequency:")
    print(
        out.head(12)[
            ["pathway", "direction", "inclusion_freq", "base_consensus_score", "robust_tier"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
