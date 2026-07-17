#!/usr/bin/env python3
"""Corrected Fig 4D — three-facet donor-level exemplars (F9, F2, F15)."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
OUT_BASE = FIGURES / "fig4d_case_study_exemplars"
WIDTH_MM, HEIGHT_MM = 178.0, 70.0

SCORES_TSV = RESULTS / "gse41177_family_sample_scores.tsv"
META_TSV = RESULTS / "gse41177_duration_metadata.tsv"
CHRONICITY_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"

AF_COLOR = "#9C5550"
FACETS = [
    {"family_id": 9, "title": "F9  HS-GAG / bone / GPCR"},
    {"family_id": 2, "title": "F2  Mito Ca2+ transport"},
    {"family_id": 15, "title": "F15  NTRK2 / TrkB"},
]


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 6.0,
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def regional_direction(rho_laa: float, rho_lapv: float) -> str:
    if not np.isfinite(rho_lapv):
        return "NA"
    return "same" if np.sign(rho_laa) == np.sign(rho_lapv) or rho_laa == 0 else "opposite"


def save_outputs(fig: plt.Figure) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BASE.with_suffix(".pdf"), bbox_inches=None, pad_inches=0.02)
    fig.savefig(OUT_BASE.with_suffix(".svg"), bbox_inches=None, pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    set_style()
    scores = pd.read_csv(SCORES_TSV, sep="\t")
    meta = pd.read_csv(META_TSV, sep="\t")
    cls = pd.read_csv(CHRONICITY_TSV, sep="\t").set_index("family_id")

    primary = meta.loc[meta["primary_duration_eligible"]].copy()
    assert len(primary) == 12

    mm = 1.0 / 25.4
    fig = plt.figure(figsize=(WIDTH_MM * mm, HEIGHT_MM * mm))
    gs = fig.add_gridspec(1, 3, left=0.06, right=0.98, top=0.78, bottom=0.18, wspace=0.28)

    for col, facet in enumerate(FACETS):
        fid = facet["family_id"]
        ax = fig.add_subplot(gs[0, col])
        sub = scores[
            (scores["family_id"] == fid)
            & (scores["region"] == "LAA")
            & scores["sample_id"].isin(primary["sample_id"])
        ].merge(
            primary[["sample_id", "donor_id", "duration_months"]],
            on="sample_id",
            how="inner",
        )
        assert len(sub) == 12, f"F{fid}: expected 12 primary AF-LAA points, got {len(sub)}"

        ax.scatter(sub["duration_months"], sub["score"], s=36, c=AF_COLOR,
                   edgecolor="#1A1A1A", linewidth=0.5, alpha=0.9, zorder=3)
        ax.set_xlabel("AF duration (mo)", fontsize=6.0)
        if col == 0:
            ax.set_ylabel("Family activity score", fontsize=6.0)
        ax.set_title(facet["title"], fontsize=6.5, fontweight="bold", loc="left")
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.tick_params(labelsize=5.5)

        row = cls.loc[fid]
        rho = float(row["rho_primary"])
        rho_lapv = float(row["rho_lapv"])
        direction = regional_direction(rho, rho_lapv)
        stats = (
            f"LAA ρ = {rho:+.2f}\n"
            f"95% CI [{float(row['ci_lo_primary']):+.2f}, {float(row['ci_hi_primary']):+.2f}]\n"
            f"p_perm = {float(row['p_perm_primary']):.3f}\n"
            f"Class: {row['chronicity_class'].replace('_', '-')}\n"
            f"LAPV ρ = {rho_lapv:+.2f}\n"
            f"Regional: {direction}"
        )
        ax.text(0.97, 0.97, stats, transform=ax.transAxes, ha="right", va="top",
                fontsize=5.2, bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                        edgecolor="#888", linewidth=0.4))

    fig.text(0.06, 0.90,
             "Fig 4d  Donor-level exemplars (12 primary AF-LAA only; no SR anchors; descriptive scatter only)",
             fontsize=7.0, fontweight="bold", color="#1A1A1A")
    fig.text(0.06, 0.06,
             "Family activity scores from exact-reference pipeline.  Points are not inferential regression fits.  "
             "F10 and F15 are Path A triple-positive families in Fig 5 (with F10 GAG/TGF-β).",
             fontsize=5.2, color="#444444")

    save_outputs(fig)
    print(f"wrote {OUT_BASE.with_suffix('.pdf').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
