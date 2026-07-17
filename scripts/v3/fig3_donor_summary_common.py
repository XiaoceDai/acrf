#!/usr/bin/env python3
"""Shared helpers for Fig 3c/3d family-level donor summary panels."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

SCORES = ROOT / "results" / "fig3_v3_singscore_per_donor.tsv"
FAM_DELTA = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
SENS = ROOT / "results" / "fig3_boundary_permutation_sensitivity.tsv"

CLASS_COLORS = {
    "Fidelity": "#7FB069",
    "Relay": "#8E6BAA",
    "Extension": "#5F9EA0",
    "Below_detection": "#A9A9A9",
}
COMPARTMENT_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC",
    "adipocyte": "Adipo",
    "neuronal": "Neuro",
}


def family_donor_frame(scores: pd.DataFrame, family_id: int, compartment: str) -> pd.DataFrame:
    g = scores.loc[
        (scores["family_id"].astype(int) == family_id)
        & (scores["compartment"] == compartment)
        & scores["score"].notna()
    ]
    return (
        g.groupby(["donor_id", "rhythm"], as_index=False)["score"]
        .mean()
        .sort_values(["rhythm", "donor_id"])
        .reset_index(drop=True)
    )


def lookup_delta(delta: pd.DataFrame, family_id: int, compartment: str) -> tuple[float, float]:
    row = delta.loc[
        (delta["family_id"].astype(int) == family_id)
        & (delta["compartment"] == compartment)
    ].iloc[0]
    return float(row["dScore"]), float(row["p_perm"])


def lookup_sensitivity(family_id: int, compartment: str) -> float | None:
    if not SENS.is_file():
        return None
    sens = pd.read_csv(SENS, sep="\t")
    hit = sens.loc[
        (sens["family_id"].astype(int) == family_id)
        & (sens["compartment"] == compartment)
    ]
    if hit.empty:
        return None
    return float(hit.iloc[0]["p_perm_sensitivity_1e6"])


def draw_donor_pair(
    ax,
    fam: pd.DataFrame,
    *,
    ctrl_color: str,
    af_color: str,
    title: str,
    subtitle: str,
    dscore: float,
    p_perm: float,
    carrier_note: str | None = None,
    sensitivity_p: float | None = None,
    highlight: bool = False,
):
    """Jittered Ctrl / AF donor points with group means."""
    rng = np.random.default_rng(abs(hash(title)) % (2**31))
    groups = [("Ctrl", 0.0, ctrl_color), ("AF", 1.0, af_color)]
    for rhythm, x, color in groups:
        vals = fam.loc[fam["rhythm"] == rhythm, "score"].to_numpy(dtype=float)
        if len(vals) == 0:
            continue
        jitter = rng.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(
            x + jitter, vals,
            s=28, color=color, edgecolor="white", linewidth=0.5, zorder=3, alpha=0.9,
        )
        ax.hlines(
            float(vals.mean()), x - 0.22, x + 0.22,
            colors=color, linewidths=1.6, zorder=4,
        )

    ax.set_xlim(-0.55, 1.55)
    ax.set_xticks([0.0, 1.0])
    ax.set_xticklabels(["Ctrl", "AF"])
    ax.set_title(title, fontsize=8.5, fontweight="bold", loc="left", pad=3)
    ax.text(0.0, 1.02, subtitle, transform=ax.transAxes, fontsize=6.2, color="#555555", va="bottom")

    annot = f"Δ = {dscore:+.4f}   p_perm = {p_perm:.4f}"
    if sensitivity_p is not None:
        annot += f"\n(sens 1e6 = {sensitivity_p:.4f})"
    ax.text(
        0.98, 0.98, annot,
        transform=ax.transAxes, ha="right", va="top", fontsize=6.0,
        color="#1A1A1A" if highlight else "#444444",
        fontweight="bold" if highlight else "normal",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#F7F7F7", edgecolor="#CCCCCC", linewidth=0.4),
    )
    if carrier_note:
        ax.text(
            0.02, 0.02, carrier_note,
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=5.8, color="#666666", style="italic",
        )
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(labelsize=7)


def load_tables():
    scores = pd.read_csv(SCORES, sep="\t")
    delta = pd.read_csv(FAM_DELTA, sep="\t")
    return scores, delta
