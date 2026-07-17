#!/usr/bin/env python3
"""Supporting Figure — family × compartment cells with primary p_perm < 0.10 trends.

Main Fig 3B shows only p < 0.05 dots; this panel retains '.' trend marks.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_style import set_nature_style, enforce_font_floor, FS_TICK, FS_PANEL_TITLE  # noqa: E402

DELTA = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
OUT_PDF = ROOT / "figures" / "figS_fig3_p010_trends.pdf"
OUT_TSV = ROOT / "results" / "fig3_v3_p010_trends.tsv"


def main() -> None:
    set_nature_style()
    d = pd.read_csv(DELTA, sep="\t")
    trend = d.loc[
        d["p_perm"].notna() & (d["p_perm"] < 0.10) & (d["p_perm"] >= 0.05)
    ].copy()
    trend = trend.sort_values(["family_id", "compartment"])
    cols = [
        "family_id", "compartment", "dScore", "p_perm", "p_perm_mc_se",
        "n_af", "n_ctrl", "sig_marker",
    ]
    trend[cols].to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")

    MM = 1.0 / 25.4
    fig, ax = plt.subplots(figsize=(140 * MM, max(50, 28 + 8 * len(trend)) * MM))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(
        0.0, 0.98,
        "Supporting  Family × compartment nominal 0.05 ≤ p_perm < 0.10 (trends)",
        fontsize=enforce_font_floor(FS_PANEL_TITLE), fontweight="bold", va="top",
        transform=ax.transAxes,
    )
    ax.text(
        0.0, 0.90,
        "Not used for Fidelity / Relay / Extension classification (main Fig 3B).",
        fontsize=enforce_font_floor(FS_TICK), color="#555555", va="top",
        transform=ax.transAxes,
    )
    header = f"{'Family':<8}{'Compartment':<22}{'ΔScore':>10}{'p_perm':>12}"
    ax.text(0.0, 0.78, header, fontsize=7.0, family="monospace", fontweight="bold",
            transform=ax.transAxes, va="top")
    y = 0.70
    if trend.empty:
        ax.text(0.0, y, "(none)", fontsize=7.0, family="monospace", transform=ax.transAxes)
    else:
        for _, r in trend.iterrows():
            line = (
                f"F{int(r.family_id):<7}{str(r.compartment):<22}"
                f"{float(r.dScore):>+10.4f}{float(r.p_perm):>12.4f}"
            )
            ax.text(0.0, y, line, fontsize=7.0, family="monospace", transform=ax.transAxes, va="top")
            y -= 0.08

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
