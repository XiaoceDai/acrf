#!/usr/bin/env python3
"""Supporting Figure — permutation null histograms for three boundary family tests.

Annotates primary 10k p_perm vs 1e6 sensitivity p (does not replace primary).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "v3"))
from fig3_donor_summary_common import family_donor_frame  # noqa: E402
from plot_style import set_nature_style, enforce_font_floor, FS_TICK, FS_PANEL_TITLE  # noqa: E402

SCORES = ROOT / "results" / "fig3_v3_singscore_per_donor.tsv"
DELTA = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
ORDER = ROOT / "results" / "fig3_v3_permutation_test_order.tsv"
SENS = ROOT / "results" / "fig3_boundary_permutation_sensitivity.tsv"
OUT_PDF = ROOT / "figures" / "figS_fig3_boundary_perm_histograms.pdf"

TARGETS = [(11, "adipocyte"), (15, "endothelial"), (16, "adipocyte")]
N_NULL = 10_000


def null_stats(values: np.ndarray, is_af: np.ndarray, seed: int, n_perm: int) -> tuple[float, np.ndarray]:
    n_af = int(is_af.sum())
    obs = float(values[is_af].mean() - values[~is_af].mean())
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        perm = rng.permutation(len(values))
        mask = np.zeros(len(values), dtype=bool)
        mask[perm[:n_af]] = True
        null[i] = float(values[mask].mean() - values[~mask].mean())
    return obs, null


def main() -> None:
    set_nature_style()
    scores = pd.read_csv(SCORES, sep="\t")
    delta = pd.read_csv(DELTA, sep="\t")
    order = pd.read_csv(ORDER, sep="\t")
    order = order.loc[order["test_level"] == "family"]
    sens = pd.read_csv(SENS, sep="\t") if SENS.is_file() else pd.DataFrame()

    MM = 1.0 / 25.4
    fig, axes = plt.subplots(1, 3, figsize=(190 * MM, 70 * MM), sharey=True)
    for ax, (fid, ct) in zip(axes, TARGETS):
        test_id = f"family|{fid}|{ct}"
        seed = int(order.loc[order["test_id"] == test_id, "rng_seed"].iloc[0])
        prim = delta.loc[(delta.family_id == fid) & (delta.compartment == ct)].iloc[0]
        fam = family_donor_frame(scores, fid, ct)
        ctrl = fam.loc[fam.rhythm == "Ctrl", "score"].to_numpy(float)
        af = fam.loc[fam.rhythm == "AF", "score"].to_numpy(float)
        values = np.concatenate([ctrl, af])
        is_af = np.concatenate([np.zeros(len(ctrl), bool), np.ones(len(af), bool)])
        obs, null = null_stats(values, is_af, seed, N_NULL)

        ax.hist(null, bins=40, color="#B0B0B0", edgecolor="white", linewidth=0.3)
        ax.axvline(obs, color="#B2182B", lw=1.4, label="observed Δ")
        ax.axvline(-obs, color="#B2182B", lw=0.8, ls="--", alpha=0.7)
        p_prim = float(prim["p_perm"])
        sens_row = sens.loc[(sens.family_id == fid) & (sens.compartment == ct)]
        p_sens = float(sens_row.iloc[0]["p_perm_sensitivity_1e6"]) if len(sens_row) else float("nan")
        ax.set_title(f"F{fid} × {ct}", fontsize=enforce_font_floor(FS_PANEL_TITLE), fontweight="bold")
        ax.text(
            0.98, 0.96,
            f"primary 10k p = {p_prim:.4f}\nsens 1e6 p = {p_sens:.4f}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=enforce_font_floor(FS_TICK - 0.5),
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#CCC", lw=0.4),
        )
        ax.set_xlabel("null ΔScore", fontsize=enforce_font_floor(FS_TICK))
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    axes[0].set_ylabel("count (of 10k null draws)", fontsize=enforce_font_floor(FS_TICK))
    fig.suptitle(
        "Supporting  Boundary permutation nulls (first 10k draws; primary vs 1e6 sensitivity annotated)",
        x=0.06, y=0.98, ha="left", fontsize=enforce_font_floor(FS_PANEL_TITLE), fontweight="bold",
    )
    fig.text(
        0.06, 0.02,
        "Sensitivity p-values do not replace publication primary p_perm.  "
        "Null histograms use the same per-test rng_seed as the locked primary schedule.",
        fontsize=5.8, color="#666666",
    )
    fig.tight_layout(rect=[0.02, 0.08, 0.98, 0.90])
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
