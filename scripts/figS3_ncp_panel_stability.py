#!/usr/bin/env python3
"""Supplementary Figure S3 — NCP panel stability."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (
    COL_DOUBLE,
    LA_COLOR,
    RA_COLOR,
    panel_label,
    set_nature_style,
    strip_spines,
)

set_nature_style()

PROGS_IN  = ROOT / "bulk" / "Project_P01" / "normal_chamber_programs.tsv"
STAB_IN   = ROOT / "bulk" / "Project_P01" / "bulk_chamber_stability.tsv"
PSET_IN   = ROOT / "results" / "figure5_program_sets.tsv"

FIG_OUT   = ROOT / "figures" / "figS3_ncp_panel_stability.pdf"
SUMM_OUT  = ROOT / "results" / "figS3_ncp_stability_summary.tsv"


def short_label(s: str) -> str:
    s = s.strip().replace(" ", " ")
    return s if len(s) <= 28 else s[:25] + "…"


def main() -> None:
    progs = pd.read_csv(PROGS_IN, sep="\t")
    progs.columns = [c.strip().strip('"') for c in progs.columns]

    progs["program"]    = progs["program"].astype(str).str.strip().str.strip('"')
    progs["direction"]  = progs["direction"].astype(str).str.strip().str.strip('"')
    progs["program_id"] = progs["program_id"].astype(str).str.strip().str.strip('"')
    progs["leading_edge_genes"] = progs["leading_edge_genes"].astype(str).str.strip().str.strip('"')
    progs["stability_consistent_runs"] = pd.to_numeric(
        progs["stability_consistent_runs"], errors="coerce")

    progs = progs.reset_index(drop=True)
    progs["program_order"] = np.arange(1, len(progs) + 1)
    progs["short"] = progs["program"].apply(short_label)

    stab = pd.read_csv(STAB_IN, sep="\t")
    stab.columns = [c.strip().strip('"') for c in stab.columns]
    stab = stab.rename(columns={"gene_symbol": "gene", "stability_score": "stability"})
    stab["gene"]      = stab["gene"].astype(str).str.strip().str.strip('"')
    stab["direction"] = stab["direction"].astype(str).str.strip().str.strip('"')
    stab["stability"] = pd.to_numeric(stab["stability"], errors="coerce")
    stab_lookup = dict(zip(stab["gene"], stab["stability"]))

    # axis assignment (optional)
    axis_lookup: dict[str, str] = {}
    if PSET_IN.exists():
        ps = pd.read_csv(PSET_IN, sep="\t")
        if {"program_id", "axis"}.issubset(ps.columns):
            axis_lookup = dict(zip(
                ps["program_id"].astype(str).str.strip(),
                ps["axis"].astype(str).str.strip(),
            ))

    # per-program LE-gene stability
    summary_rows = []
    for _, row in progs.iterrows():
        le = [g.strip() for g in row["leading_edge_genes"].split(",")
              if g.strip()]
        scores = [stab_lookup.get(g) for g in le]
        scores = [s for s in scores if not pd.isna(s)]
        n_le = len(le)
        n_le_in = len(scores)
        mean_stab = float(np.mean(scores)) if scores else np.nan
        frac_full = float(np.mean(np.array(scores) == 8)) if scores else np.nan
        summary_rows.append({
            "program_order": row["program_order"],
            "program_id":    row["program_id"],
            "program":       row["program"],
            "short":         row["short"],
            "direction":     row["direction"],
            "axis":          axis_lookup.get(row["program_id"], ""),
            "n_LE_total":    n_le,
            "n_LE_in_stab":  n_le_in,
            "LE_mean_stability_runs":   round(mean_stab, 2) if not pd.isna(mean_stab) else None,
            "LE_frac_stable_in_all_runs": round(frac_full, 2) if not pd.isna(frac_full) else None,
            "program_consistent_runs":  int(row["stability_consistent_runs"])
                if not pd.isna(row["stability_consistent_runs"]) else None,
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(SUMM_OUT, sep="\t", index=False)
    print(f"wrote {SUMM_OUT.relative_to(ROOT)}")

    fig = plt.figure(figsize=(COL_DOUBLE, 0.5 * COL_DOUBLE))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.4, 1.2, 1.0], wspace=0.55)

    # Panel a — per-program direction-stability bar (program_consistent_runs / 8)
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(summary))[::-1]
    bar_colors = [LA_COLOR if d == "LA-biased" else RA_COLOR
                  for d in summary["direction"]]
    ax.barh(y, summary["program_consistent_runs"], color=bar_colors,
            edgecolor="black", linewidth=0.3)
    ax.axvline(8, color="#444444", linestyle="--", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(summary["short"])
    ax.set_xlim(0, 8.4)
    ax.set_xlabel("direction-stable resamplings (max = 8)")
    ax.set_title("a  Program-level stability", loc="left", pad=4)
    ax.set_xticks([0, 2, 4, 6, 8])
    strip_spines(ax)
    panel_label(ax, "a", x=-0.55, y=1.02)

    # legend
    handles = [
        mpatches.Patch(facecolor=LA_COLOR, label="LA-biased"),
        mpatches.Patch(facecolor=RA_COLOR, label="RA-biased"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False)

    # Panel b — per-program LE-gene mean stability
    ax = fig.add_subplot(gs[0, 1])
    ax.barh(y, summary["LE_mean_stability_runs"], color=bar_colors,
            edgecolor="black", linewidth=0.3)
    ax.axvline(8, color="#444444", linestyle="--", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlim(0, 8.4)
    ax.set_xlabel("mean LE-gene stability (runs / 8)")
    ax.set_title("b  Leading-edge gene stability", loc="left", pad=4)
    ax.set_xticks([0, 2, 4, 6, 8])
    strip_spines(ax)
    for yi, frac in zip(y, summary["LE_frac_stable_in_all_runs"]):
        if not pd.isna(frac):
            ax.text(8.5, yi, f"{frac:.2f}", fontsize=5,
                    ha="left", va="center", color="#444444")
    panel_label(ax, "b", x=-0.10, y=1.02)
    ax.text(8.5, len(y), "fraction\nLE @ 8/8",
            fontsize=4.5, ha="left", va="bottom", color="#444444")

    # Panel c — overall gene-level stability distribution by direction
    ax = fig.add_subplot(gs[0, 2])
    bins = np.arange(0.5, 9.5, 1)
    for d, col in [("LA-high", LA_COLOR), ("RA-high", RA_COLOR)]:
        sub = stab.query("direction == @d and stability == stability")
        ax.hist(sub["stability"], bins=bins, alpha=0.55, color=col,
                edgecolor="white", linewidth=0.3, label=d)
    ax.set_xlabel("stability score (resamplings significant)")
    ax.set_ylabel("genes (n)")
    ax.set_xticks(range(1, 9))
    ax.set_title("c  Gene-level stability across resamplings", loc="left", pad=4)
    ax.legend(loc="upper left", frameon=False)
    strip_spines(ax)
    panel_label(ax, "c", x=-0.20, y=1.02)

    fig.suptitle(
        "Normal Chamber Prior — panel and gene-level stability "
        "across donor resamplings (GSE112339, 8 paired donors)",
        y=1.02, fontsize=8,
    )

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
