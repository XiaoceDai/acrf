#!/usr/bin/env python3
"""Supplementary Figure S2 — Donor pairing & sample-attribute heatmap."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import COL_DOUBLE, panel_label, set_nature_style, strip_spines

set_nature_style()

OUT  = ROOT / "figures" / "figS2_pairing_attribute_heatmap.pdf"
SUMM = ROOT / "results" / "figS2_pairing_summary.tsv"


AF_COLOR    = "#B2182B"
SR_COLOR    = "#2166AC"
MALE_COLOR  = "#5AAE61"
FEMALE_COLOR = "#9970AB"
NA_COLOR    = "#BDBDBD"
PRESENT_BORDER = "#222222"


COHORTS = [
    {
        "name": "GSE128188",
        "path": ROOT / "Project_P1" / "GSE128188" / "sample_metadata.tsv",
        "axis": "chamber",
        "axis_levels": ["LA", "RA"],
    },
    {
        "name": "GSE79768",
        "path": ROOT / "Project_P1" / "GSE79768" / "sample_metadata.tsv",
        "axis": "chamber",
        "axis_levels": ["LA", "RA"],
    },
    {
        "name": "GSE115574",
        "path": ROOT / "Project_P1" / "GSE115574" / "sample_metadata.tsv",
        "axis": "chamber",
        "axis_levels": ["LA", "RA"],
    },
    {
        "name": "GSE41177",
        "path": ROOT / "Project_P4" / "GSE41177" / "sample_metadata.tsv",
        "axis": "tissue",
        "axis_levels": ["LAA", "LAPV"],
    },
]


def normalise_sex(s: str | float) -> str:
    if pd.isna(s):
        return "NA"
    s = str(s).strip().lower()
    return {"m": "male", "f": "female"}.get(s, s if s in ("male", "female") else "NA")


def draw_one(ax: plt.Axes, ax_cov: plt.Axes, cfg: dict,
             panel_letter: str | None) -> dict:
    md = pd.read_csv(cfg["path"], sep="\t")
    axis = cfg["axis"]
    axis_levels = cfg["axis_levels"]

    md = md.copy()
    md["_rhythm"] = md["rhythm"].fillna("NA")
    md["_sex"]    = md["sex"].apply(normalise_sex) if "sex" in md.columns else "NA"
    md["_age"]    = pd.to_numeric(md.get("age", pd.Series(np.nan, index=md.index)),
                                  errors="coerce")

    rhythm_priority = {"SR": 0, "AF": 1}
    md["_rhythm_pri"] = md["_rhythm"].map(rhythm_priority).fillna(99)
    donors = (md.drop_duplicates("donor_id")
                .sort_values(["_rhythm_pri", "donor_id"])
                ["donor_id"].tolist())
    donor_to_idx = {d: i for i, d in enumerate(donors)}

    n_d = len(donors)
    n_l = len(axis_levels)

    # ---- main matrix : presence (1) / absence (0) per (donor, axis_level) --
    matrix = np.zeros((n_d, n_l))
    rhythm_per_donor = {}
    sex_per_donor = {}
    age_per_donor = {}
    for _, row in md.iterrows():
        d  = row["donor_id"]
        al = row[axis]
        if d not in donor_to_idx or al not in axis_levels:
            continue
        i = donor_to_idx[d]
        j = axis_levels.index(al)
        matrix[i, j] = 1
        rhythm_per_donor[d] = row["_rhythm"]
        sex_per_donor[d]    = row["_sex"]
        if not pd.isna(row["_age"]):
            age_per_donor[d] = row["_age"]

    for i, d in enumerate(donors):
        rh_col = AF_COLOR if rhythm_per_donor.get(d) == "AF" \
                 else SR_COLOR if rhythm_per_donor.get(d) == "SR" else NA_COLOR
        for j in range(n_l):
            facecolor = rh_col if matrix[i, j] == 1 else "#FFFFFF"
            edgecolor = PRESENT_BORDER if matrix[i, j] == 1 else "#DDDDDD"
            ax.add_patch(mpatches.Rectangle(
                (j - 0.5, i - 0.5), 1, 1,
                facecolor=facecolor, edgecolor=edgecolor, linewidth=0.4,
            ))
    ax.set_xlim(-0.6, n_l - 0.4)
    ax.set_ylim(n_d - 0.5, -0.5)
    ax.set_xticks(range(n_l))
    ax.set_xticklabels(axis_levels)
    ax.set_yticks(range(n_d))
    ax.set_yticklabels(donors, fontsize=4)
    ax.tick_params(length=0)
    ax.set_title(cfg["name"], pad=2)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    if panel_letter:
        panel_label(ax, panel_letter, x=-0.45, y=1.04)

    for i, d in enumerate(donors):
        # sex tile
        sx = sex_per_donor.get(d, "NA")
        sx_col = MALE_COLOR if sx == "male" else \
                 FEMALE_COLOR if sx == "female" else NA_COLOR
        ax_cov.add_patch(mpatches.Rectangle(
            (-0.5, i - 0.5), 1, 1,
            facecolor=sx_col, edgecolor="#DDDDDD", linewidth=0.3,
        ))
        # age tile (continuous greyscale)
        if d in age_per_donor:
            a = age_per_donor[d]
            shade = 1.0 - min(max((a - 30) / 50, 0.05), 0.95)
            ax_cov.text(0.5 + 0.5, i, f"{int(a)}", fontsize=4,
                        ha="center", va="center", color="black")
            ax_cov.add_patch(mpatches.Rectangle(
                (0.5, i - 0.5), 1, 1,
                facecolor=(shade, shade, shade), edgecolor="#DDDDDD", linewidth=0.3,
            ))
        else:
            ax_cov.add_patch(mpatches.Rectangle(
                (0.5, i - 0.5), 1, 1,
                facecolor="#FFFFFF", edgecolor="#DDDDDD", linewidth=0.3,
            ))
            ax_cov.text(1.0, i, "·", fontsize=5, ha="center", va="center")
    ax_cov.set_xlim(-0.6, 1.6)
    ax_cov.set_ylim(n_d - 0.5, -0.5)
    ax_cov.set_xticks([0, 1])
    ax_cov.set_xticklabels(["sex", "age"])
    ax_cov.set_yticks([])
    ax_cov.tick_params(length=0)
    for s in ("top", "right", "left", "bottom"):
        ax_cov.spines[s].set_visible(False)

    n_paired = int(np.sum(matrix.sum(axis=1) == n_l))
    n_unpaired = int(np.sum((matrix.sum(axis=1) > 0) & (matrix.sum(axis=1) < n_l)))
    n_af = int(sum(1 for d in donors if rhythm_per_donor.get(d) == "AF"))
    n_sr = int(sum(1 for d in donors if rhythm_per_donor.get(d) == "SR"))
    n_other = n_d - n_af - n_sr

    return {
        "cohort": cfg["name"],
        "axis": axis,
        "axis_levels": ",".join(axis_levels),
        "n_donors": n_d,
        "n_paired_complete": n_paired,
        "n_unpaired": n_unpaired,
        "n_AF_donors": n_af,
        "n_SR_donors": n_sr,
        "n_other_donors": n_other,
    }


def main() -> None:
    fig = plt.figure(figsize=(COL_DOUBLE, 0.55 * COL_DOUBLE))
    panel_letters = ["a", "b", "c", "d"]

    # 4 cohorts laid out in a 1x4 grid; each cohort has main + covariate stripes
    outer = fig.add_gridspec(1, 4, wspace=0.7)
    summary_rows = []
    for k, cfg in enumerate(COHORTS):
        inner = outer[0, k].subgridspec(1, 2, width_ratios=[2.0, 1.4], wspace=0.05)
        ax_main = fig.add_subplot(inner[0, 0])
        ax_cov  = fig.add_subplot(inner[0, 1])
        summary_rows.append(
            draw_one(ax_main, ax_cov, cfg, panel_letters[k])
        )

    legend_handles = [
        mpatches.Patch(facecolor=SR_COLOR, edgecolor=PRESENT_BORDER, label="SR sample present"),
        mpatches.Patch(facecolor=AF_COLOR, edgecolor=PRESENT_BORDER, label="AF sample present"),
        mpatches.Patch(facecolor="#FFFFFF", edgecolor="#999999", label="missing"),
        mpatches.Patch(facecolor=MALE_COLOR,   label="male"),
        mpatches.Patch(facecolor=FEMALE_COLOR, label="female"),
        mpatches.Patch(facecolor=NA_COLOR,     label="not reported"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncols=6,
               frameon=False, bbox_to_anchor=(0.5, -0.04), fontsize=6,
               handlelength=1.2, columnspacing=1.2)

    fig.suptitle(
        "Donor × chamber / region pairing structure and covariates",
        y=1.00, fontsize=8,
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")

    pd.DataFrame(summary_rows).to_csv(SUMM, sep="\t", index=False)
    print(f"wrote {SUMM.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
