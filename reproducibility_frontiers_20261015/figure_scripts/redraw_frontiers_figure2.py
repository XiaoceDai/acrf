#!/usr/bin/env python3
"""Redraw Frontiers Figure 2 from the locked benchmark tables. Numbers are not recomputed."""

import io
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path("/home/bio/Project_P02")
ZIP = ROOT / "submission/frontiers_bioinformatics/supplementary/Supplementary_benchmark_tables.zip"
OUT = ROOT / "submission/frontiers_bioinformatics/figures/Figure_2.pdf"

METHODS = [
    ("A1_FULL_ACRF", None, "Full acrf", "#0F6E6E"),
    ("A2_NO_DONOR_PAIRING", None, "Unpaired\nchamber-mean contrast", "#3E8F89"),
    ("A3_NO_DIRECTION_ALIGNMENT", None, "No direction alignment", "#7FB3AD"),
    ("A4_NO_HEALTHY_MAGNITUDE_SCALING", None, "No scaling", "#5C7C8A"),
    ("A5_CHAMBER_AGNOSTIC", None, "Average of both chambers", "#D4762C"),
    ("B1_PLAIN_SINGSCORE", "LA", "singscore, left atrium", "#5C4B8A"),
    ("B1_PLAIN_SINGSCORE", "RA", "singscore, right atrium", "#8E82B8"),
]
PAIRS = [
    ("GSE128188|GSE79768", "o", "GSE128188–GSE79768"),
    ("GSE128188|GSE115574", "s", "GSE128188–GSE115574"),
    ("GSE79768|GSE115574", "D", "GSE79768–GSE115574"),
]


def load(name: str) -> pd.DataFrame:
    with zipfile.ZipFile(ZIP) as zf:
        return pd.read_csv(io.BytesIO(zf.read(name)), sep="\t")


def main() -> None:
    summary = load("benchmark_summary.tsv")
    corr = load("benchmark_pairwise_correlations.tsv")
    rows = []
    for method_id, stratum, label, color in METHODS:
        hit = summary[summary.method_id == method_id]
        if stratum is not None:
            hit = hit[hit.stratum == stratum]
        else:
            hit = hit[hit.stratum == "unstratified"]
        if len(hit) != 1:
            raise SystemExit(f"expected one summary row for {method_id} {stratum}, got {len(hit)}")
        rec = hit.iloc[0]
        pair = corr[(corr.method_id == method_id) & (corr.family_id.astype(str) == "ALL") & (corr.cohort != "SUMMARY")]
        if stratum is not None:
            pair = pair[pair.stratum == stratum]
        rhos = {}
        for cohort, _marker, _name in PAIRS:
            val = pair.loc[pair.cohort == cohort, "pairwise_rho"]
            if len(val) != 1:
                raise SystemExit(f"missing rho {method_id} {stratum} {cohort}")
            rhos[cohort] = float(val.iloc[0])
        rows.append(
            dict(
                label=label,
                color=color,
                concordant=int(rec.M2_concordant_count),
                denominator=int(rec.M2_denominator),
                median=float(rec.M1_median_rho),
                rhos=rhos,
                retained=rec.M5_retained_vs_a1,
                comparable=str(rec.M5_retained_vs_a1) != "NOT_COMPARABLE",
            )
        )

    expected = {
        "Full acrf": (32, 0.179412, "16"),
        "Unpaired\nchamber-mean contrast": (32, 0.214706, "12"),
        "No direction alignment": (32, 0.558824, "16"),
        "No scaling": (32, 0.188235, "16"),
        "Average of both chambers": (20, -0.294118, "8"),
        "singscore, left atrium": (22, -0.014706, "NOT_COMPARABLE"),
        "singscore, right atrium": (24, -0.182353, "NOT_COMPARABLE"),
    }
    for row in rows:
        count, median, retained = expected[row["label"]]
        if row["concordant"] != count or abs(row["median"] - median) > 1e-5 or str(row["retained"]) != retained:
            raise SystemExit(f"locked number mismatch: {row}")
        if row["denominator"] != 48:
            raise SystemExit(row)

    high = {
        "Full acrf": 6,
        "Unpaired\nchamber-mean contrast": 4,
        "No direction alignment": 6,
        "No scaling": 6,
        "Average of both chambers": 0,
    }
    cls = load("benchmark_class_stability.tsv")
    pooled = cls[(cls.cohort == "POOLED") & (cls.family_id.astype(str) != "ALL")]
    for method_id, stratum, label, _color in METHODS:
        if label not in high:
            continue
        sub = pooled[pooled.method_id == method_id]
        n_high = int((sub.method_class == "High").sum())
        if n_high != high[label] or len(sub) != 16:
            raise SystemExit(f"high-count mismatch {label}: {n_high}")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.6,
        }
    )
    fig = plt.figure(figsize=(7.05, 4.15))
    ax_b = fig.add_axes([0.30, 0.24, 0.22, 0.66])
    ax_c = fig.add_axes([0.55, 0.24, 0.20, 0.66])
    # Panel C covers only the five acrf rows, aligned with the top of A and B.
    ax_d = fig.add_axes([0.80, 0.429, 0.17, 0.471])

    y = list(range(len(rows) - 1, -1, -1))
    for ax in (ax_b, ax_c, ax_d):
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=7, length=3, width=0.6)

    ax_b.barh(
        y,
        [r["concordant"] / r["denominator"] for r in rows],
        color=[r["color"] for r in rows],
        height=0.55,
        zorder=2,
    )
    for yi, row in zip(y, rows):
        ax_b.text(
            row["concordant"] / row["denominator"] + 0.02,
            yi,
            f"{row['concordant']}/48",
            va="center",
            ha="left",
            fontsize=6.5,
            color="#333333",
        )
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([r["label"] for r in rows], fontsize=7.5)
    ax_b.tick_params(axis="y", pad=9)
    ax_b.set_xlim(0, 1.18)
    ax_b.set_xlabel("Family-effect sign agreement", fontsize=7.5)
    ax_b.set_title("A", loc="left", fontsize=11, fontweight="bold", pad=4)

    for cohort, marker, _name in PAIRS:
        ax_c.scatter(
            [r["rhos"][cohort] for r in rows],
            y,
            marker=marker,
            s=22,
            c=[r["color"] for r in rows],
            zorder=3,
            linewidths=0,
        )
    ax_c.scatter(
        [r["median"] for r in rows],
        y,
        marker="|",
        s=90,
        c="black",
        linewidths=1.2,
        zorder=4,
    )
    ax_c.axvline(0, color="#DDDDDD", lw=0.6, zorder=1)
    ax_c.set_yticks(y)
    ax_c.set_yticklabels([])
    ax_c.set_xlim(-0.9, 0.95)
    ax_c.set_xlabel("Spearman correlation", fontsize=7.5)
    ax_c.set_title("B", loc="left", fontsize=11, fontweight="bold", pad=4)

    d_rows = [(yi, row) for yi, row in zip(y, rows) if row["label"] in high]
    for yi, row in d_rows:
        retained = int(row["retained"])
        n_high = high[row["label"]]
        ax_d.plot([n_high, retained], [yi, yi], color="#D0D0D0", lw=1, zorder=1)
        ax_d.scatter([retained], [yi], s=28, facecolors="white", edgecolors=row["color"], linewidths=1.1, zorder=3)
        ax_d.scatter([n_high], [yi], s=26, marker="D", c=row["color"], zorder=4, linewidths=0)
    ax_d.set_ylim(1.4, 6.4)
    ax_d.set_yticks([yi for yi, row in d_rows])
    ax_d.set_yticklabels([])
    ax_d.set_xlim(-0.8, 17.2)
    ax_d.set_xlabel("Families", fontsize=7.5)
    ax_d.set_title("C", loc="left", fontsize=11, fontweight="bold", pad=4)

    pair_handles = [
        plt.Line2D([0], [0], marker=marker, color="none", markerfacecolor="#555555", markeredgecolor="none", markersize=5, label=name)
        for _cohort, marker, name in PAIRS
    ]
    pair_handles.append(
        plt.Line2D([0], [0], marker="|", color="black", linestyle="None", markersize=8, markeredgewidth=1.2, label="Median")
    )
    class_handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#555555", markersize=5.5, label="Same class as full acrf"),
        plt.Line2D([0], [0], marker="D", color="none", markerfacecolor="#555555", markeredgecolor="none", markersize=4.5, label="High"),
    ]
    fig.legend(
        handles=pair_handles + class_handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=6.4,
        bbox_to_anchor=(0.58, 0.01),
    )
    fig.savefig(OUT, dpi=300)
    print(OUT)


if __name__ == "__main__":
    main()
