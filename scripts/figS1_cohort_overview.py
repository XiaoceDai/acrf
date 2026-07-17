#!/usr/bin/env python3
"""
Supplementary Figure S1 — Cohort overview.

Inputs : results/cohort_demographics_long.tsv   (one row per donor)
Output : figures/figS1_cohort_overview.pdf
         results/figS1_cohort_overview_summary.tsv

Panels:
  a — donor counts per cohort, stacked AF / SR
  b — age distribution per cohort, split AF / SR
  c — sex composition per cohort, stacked male / female / NA
  d — disease background composition per cohort

This figure summarises the study population across eight GEO cohorts
(including GSE112339 healthy reference) and matches
results/cohort_demographics_table.tsv (Supplementary Table S1 source). Chamber
column matches paired LA/RA (or LA+LAPV, etc.) sampling layout; Fig 2 GEO cohorts are
restricted to donors present in innovation_ncp_v3_acai_long.tsv when tables are
regenerated via scripts/p1_demographics_table.py.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (
    AF_RHYTHM_COLOR,
    COL_DOUBLE,
    FS_ANNOT,
    HEALTHY_REF_COLOR,
    QUALITATIVE_CATEGORY_COLORS,
    SR_RHYTHM_COLOR,
    SEX_FEMALE_COLOR,
    SEX_MALE_COLOR,
    SEX_NA_COLOR,
    cohort_axis_labels,
    panel_label,
    set_nature_style,
    strip_spines,
)

set_nature_style()

IN  = ROOT / "results" / "cohort_demographics_long.tsv"
OUT = ROOT / "figures"  / "figS1_cohort_overview.pdf"
OUT_SUMMARY = ROOT / "results" / "figS1_cohort_overview_summary.tsv"


def cohort_order(df: pd.DataFrame) -> list[str]:
    """Sort cohorts by manuscript role: discovery -> replication -> sn -> chronicity."""
    role_order = [
        "GSE112339", "GSE128188", "GSE79768", "GSE115574",
        "GSE255612", "GSE41177", "GSE238242", "GSE261363",
    ]
    seen = list(df["cohort"].drop_duplicates())
    return [c for c in role_order if c in seen] + [c for c in seen if c not in role_order]


def main() -> None:
    df = pd.read_csv(IN, sep="\t")
    cohorts = cohort_order(df)

    fig = plt.figure(figsize=(COL_DOUBLE, 0.55 * COL_DOUBLE))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.42)

    # --- Panel a : donor counts stacked AF/SR --------------------------------
    ax = fig.add_subplot(gs[0, 0])
    counts = (df.groupby(["cohort", "rhythm"]).size().unstack().fillna(0)
                .reindex(cohorts))
    af = counts.get("AF", pd.Series(0, index=cohorts))
    sr = counts.get("SR", pd.Series(0, index=cohorts))
    other = counts.drop(columns=["AF", "SR"], errors="ignore").sum(axis=1) \
            if any(c not in {"AF", "SR"} for c in counts.columns) else pd.Series(0, index=cohorts)
    x = np.arange(len(cohorts))
    ax.bar(x, sr, color=SR_RHYTHM_COLOR, label="SR", width=0.68, linewidth=0, zorder=2)
    ax.bar(x, af, bottom=sr, color=AF_RHYTHM_COLOR, label="AF", width=0.68, linewidth=0, zorder=2)
    if other.sum() > 0:
        ax.bar(x, other, bottom=sr + af, color=HEALTHY_REF_COLOR,
               label="healthy ref. / other", width=0.68, linewidth=0, zorder=2)
    for xi, total in zip(x, sr + af + other):
        ax.text(xi, total + 0.5, f"n = {int(total)}",
                ha="center", va="bottom", fontsize=FS_ANNOT, color="#333333")
    cohort_axis_labels(ax, cohorts)
    ax.set_ylabel("donors (n)")
    ax.set_ylim(0, max((sr + af + other).max() * 1.18, 5))
    ax.legend(loc="upper right", frameon=False)
    strip_spines(ax)
    panel_label(ax, "a")

    # --- Panel b : age distribution per cohort, split AF/SR ------------------
    ax = fig.add_subplot(gs[0, 1])
    plot_data = []
    plot_pos = []
    plot_colors = []
    plot_labels = []
    for i, c in enumerate(cohorts):
        for j, rhythm in enumerate(["SR", "AF"]):
            ages = pd.to_numeric(
                df.query("cohort == @c and rhythm == @rhythm")["age"],
                errors="coerce",
            ).dropna().tolist()
            if not ages:
                continue
            plot_data.append(ages)
            plot_pos.append(i + (j - 0.5) * 0.32)
            plot_colors.append(
                SR_RHYTHM_COLOR if rhythm == "SR" else AF_RHYTHM_COLOR)
            plot_labels.append(rhythm)

    if plot_data:
        bp = ax.boxplot(plot_data, positions=plot_pos, widths=0.28,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", linewidth=0.8),
                        whiskerprops=dict(linewidth=0.6),
                        capprops=dict(linewidth=0.6),
                        boxprops=dict(linewidth=0.5))
        for patch, col in zip(bp["boxes"], plot_colors):
            patch.set_facecolor(col)
            patch.set_alpha(0.55)
        for pos, ages, col in zip(plot_pos, plot_data, plot_colors):
            ax.scatter(np.full(len(ages), pos)
                       + np.random.RandomState(0).uniform(-0.05, 0.05, len(ages)),
                       ages, s=6, color=col, edgecolor="white",
                       linewidth=0.25, zorder=3, alpha=0.85)
    cohort_axis_labels(ax, cohorts)
    ax.set_ylabel("age (years)")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=SR_RHYTHM_COLOR, alpha=0.55, label="SR"),
        plt.Rectangle((0, 0), 1, 1, color=AF_RHYTHM_COLOR, alpha=0.55, label="AF"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False)
    strip_spines(ax)
    panel_label(ax, "b")

    # --- Panel c : sex composition stacked ----------------------------------
    ax = fig.add_subplot(gs[1, 0])
    sex_norm = df["sex"].fillna("NA").astype(str).str.strip().str.lower()
    sex_norm = sex_norm.replace({
        "m": "male", "f": "female", "": "NA", "nan": "NA",
        "unknown": "NA", "na": "NA",
    })
    df = df.assign(_sex=sex_norm)
    sex_counts = (df.groupby(["cohort", "_sex"]).size().unstack(fill_value=0)
                    .reindex(cohorts).fillna(0))
    male   = sex_counts.get("male",   pd.Series(0, index=cohorts))
    female = sex_counts.get("female", pd.Series(0, index=cohorts))
    na     = sex_counts.get("NA",     pd.Series(0, index=cohorts))
    totals = (male + female + na).replace(0, np.nan)
    male_p   = (male   / totals * 100).fillna(0)
    female_p = (female / totals * 100).fillna(0)
    na_p     = (na     / totals * 100).fillna(0)
    x = np.arange(len(cohorts))
    ax.bar(x, male_p, color=SEX_MALE_COLOR, width=0.68, linewidth=0, label="male")
    ax.bar(x, female_p, bottom=male_p, color=SEX_FEMALE_COLOR,
           width=0.68, linewidth=0, label="female")
    ax.bar(x, na_p, bottom=male_p + female_p, color=SEX_NA_COLOR,
           width=0.68, linewidth=0, label="not reported")
    cohort_axis_labels(ax, cohorts)
    ax.set_ylabel("sex composition (%)")
    ax.set_ylim(0, 105)
    ax.legend(loc="lower right", frameon=False, ncols=3, columnspacing=0.6)
    strip_spines(ax)
    panel_label(ax, "c")

    # --- Panel d : disease background composition --------------------------
    ax = fig.add_subplot(gs[1, 1])
    db_col = df.get("disease_background", pd.Series(["unknown"] * len(df))).fillna("not reported")
    db_norm = db_col.astype(str).str.strip().replace({"": "not reported", "nan": "not reported"})
    df = df.assign(_db=db_norm)
    db_counts = (df.groupby(["cohort", "_db"]).size().unstack(fill_value=0)
                   .reindex(cohorts).fillna(0))

    backgrounds = list(db_counts.columns)
    bottoms = pd.Series(0.0, index=cohorts)
    x = np.arange(len(cohorts))
    totals = db_counts.sum(axis=1).replace(0, np.nan)
    for i, b in enumerate(backgrounds):
        share = (db_counts[b] / totals * 100).fillna(0)
        ax.bar(
            x, share, bottom=bottoms,
            color=QUALITATIVE_CATEGORY_COLORS[i % len(QUALITATIVE_CATEGORY_COLORS)],
            width=0.68, linewidth=0, label=str(b),
        )
        bottoms = bottoms + share
    cohort_axis_labels(ax, cohorts)
    ax.set_ylabel("disease background (%)")
    ax.set_ylim(0, 105)
    ax.legend(
        loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False,
        fontsize=FS_ANNOT, handlelength=1.0,
    )
    strip_spines(ax)
    panel_label(ax, "d")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")

    # ---- summary deliverable ------------------------------------------------
    summary_rows = []
    for c in cohorts:
        sub = df.query("cohort == @c")
        ages = pd.to_numeric(sub["age"], errors="coerce").dropna()
        sex_n = sub["_sex"].value_counts().to_dict()
        summary_rows.append({
            "cohort": c,
            "n_donors": len(sub),
            "n_AF": int((sub["rhythm"] == "AF").sum()),
            "n_SR": int((sub["rhythm"] == "SR").sum()),
            "age_mean": round(ages.mean(), 1) if len(ages) else None,
            "age_sd":   round(ages.std(),  1) if len(ages) else None,
            "age_min":  int(ages.min()) if len(ages) else None,
            "age_max":  int(ages.max()) if len(ages) else None,
            "n_male":          sex_n.get("male", 0),
            "n_female":        sex_n.get("female", 0),
            "n_sex_NA":        sex_n.get("NA", 0),
            "disease_backgrounds": " | ".join(
                f"{b}:{int(n)}" for b, n in
                sub["_db"].value_counts().items()
            ),
        })
    pd.DataFrame(summary_rows).to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
