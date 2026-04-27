#!/usr/bin/env python3
"""
Benchmark v2 vs v3 NCP panels.

Inputs
------
results/fig1c_ncp_panel_v2.tsv
results/innovation_ncp_v3_panel.tsv
results/innovation_ncp_v2_acai_long.tsv
results/innovation_ncp_v3_acai_long.tsv

Outputs
-------
results/innovation_v2_vs_v3_metrics.tsv
results/innovation_v2_vs_v3_axis_effects.tsv
results/innovation_v2_vs_v3_decision_table.tsv
figures/innovation_v2_vs_v3_benchmark.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plot_style import COL_DOUBLE, LA_COLOR, RA_COLOR, set_nature_style

ROOT = Path(__file__).resolve().parents[2]

PANEL_V2 = ROOT / "results" / "fig1c_ncp_panel_v2.tsv"
PANEL_V3 = ROOT / "results" / "innovation_ncp_v3_panel.tsv"
LONG_V2 = ROOT / "results" / "innovation_ncp_v2_acai_long.tsv"
LONG_V3 = ROOT / "results" / "innovation_ncp_v3_acai_long.tsv"

GMT_HALL = ROOT / "bulk" / "Project_P01" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC = ROOT / "bulk" / "Project_P01" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"

OUT_METRICS = ROOT / "results" / "innovation_v2_vs_v3_metrics.tsv"
OUT_AXIS = ROOT / "results" / "innovation_v2_vs_v3_axis_effects.tsv"
OUT_DECISION = ROOT / "results" / "innovation_v2_vs_v3_decision_table.tsv"
OUT_PDF = ROOT / "figures" / "innovation_v2_vs_v3_benchmark.pdf"
OUT_FIG1D = ROOT / "figures" / "fig1d_v2_to_v3_upgrade_benchmark.pdf"
OUT_FIGS13 = ROOT / "figures" / "figS13_v2_vs_v3_benchmark.pdf"


def load_gmt(path: Path, prefix: str) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            out[f"{prefix}__{parts[0]}"] = {g for g in parts[2:] if g}
    return out


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 0.0
    u = a | b
    if not u:
        return 0.0
    return len(a & b) / len(u)


def mean_pairwise_jaccard(paths: list[str], gsets: dict[str, set[str]]) -> float:
    vals: list[float] = []
    for i in range(len(paths)):
        for j in range(i + 1, len(paths)):
            vals.append(jaccard(gsets.get(paths[i], set()), gsets.get(paths[j], set())))
    if not vals:
        return float("nan")
    return float(np.mean(vals))


def hedges_g(af: np.ndarray, sr: np.ndarray) -> float:
    n1 = len(af)
    n0 = len(sr)
    if n1 < 2 or n0 < 2:
        return float("nan")
    s1 = float(np.std(af, ddof=1))
    s0 = float(np.std(sr, ddof=1))
    pooled_var = ((n1 - 1) * (s1**2) + (n0 - 1) * (s0**2)) / (n1 + n0 - 2)
    if pooled_var <= 0:
        return float("nan")
    d = (float(np.mean(af)) - float(np.mean(sr))) / np.sqrt(pooled_var)
    j = 1 - (3 / (4 * (n1 + n0) - 9))
    return j * d


def build_axis_effects(long_df: pd.DataFrame, version: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort, dfc in long_df.groupby("cohort"):
        for axis_name, mask in [
            ("all-programs", np.ones(len(dfc), dtype=bool)),
            ("LA-biased", dfc["direction"].eq("LA-biased").to_numpy()),
            ("RA-biased", dfc["direction"].eq("RA-biased").to_numpy()),
        ]:
            sub = dfc.loc[mask].copy()
            donor = (
                sub.groupby(["cohort", "donor_id", "rhythm"], as_index=False)["aCAI"]
                .mean()
                .rename(columns={"aCAI": "axis_aCAI"})
            )
            af = donor.loc[donor["rhythm"] == "AF", "axis_aCAI"].to_numpy(float)
            sr = donor.loc[donor["rhythm"] == "SR", "axis_aCAI"].to_numpy(float)
            rows.append(
                {
                    "version": version,
                    "cohort": cohort,
                    "axis": axis_name,
                    "n_af": len(af),
                    "n_sr": len(sr),
                    "mean_af": float(np.mean(af)) if len(af) else np.nan,
                    "mean_sr": float(np.mean(sr)) if len(sr) else np.nan,
                    "delta_af_minus_sr": float(np.mean(af) - np.mean(sr)) if len(af) and len(sr) else np.nan,
                    "hedges_g_af_minus_sr": hedges_g(af, sr),
                }
            )
    out = pd.DataFrame(rows)
    for c in ["mean_af", "mean_sr", "delta_af_minus_sr", "hedges_g_af_minus_sr"]:
        out[c] = out[c].round(4)
    return out


def build_decision_table(metrics: pd.DataFrame, axis: pd.DataFrame) -> pd.DataFrame:
    m = metrics.iloc[0]
    all_axis = axis[axis["axis"] == "all-programs"].copy()
    d2 = all_axis[all_axis["version"] == "v2"].set_index("cohort")
    d3 = all_axis[all_axis["version"] == "v3"].set_index("cohort")
    cohorts = sorted(set(d2.index) & set(d3.index))

    delta_abs = []
    g_abs = []
    sign_flip = 0
    for c in cohorts:
        dv2 = float(d2.loc[c, "delta_af_minus_sr"])
        dv3 = float(d3.loc[c, "delta_af_minus_sr"])
        gv2 = float(d2.loc[c, "hedges_g_af_minus_sr"])
        gv3 = float(d3.loc[c, "hedges_g_af_minus_sr"])
        delta_abs.append(abs(dv3 - dv2))
        g_abs.append(abs(gv3 - gv2))
        if dv2 * dv3 < 0:
            sign_flip += 1

    mean_abs_delta_shift = float(np.mean(delta_abs)) if delta_abs else np.nan
    mean_abs_g_shift = float(np.mean(g_abs)) if g_abs else np.nan
    jaccard_increase = float(m["mean_pairwise_jaccard_v3"] - m["mean_pairwise_jaccard_v2"])

    rows = [
        {
            "criterion": "Panel continuity",
            "rule": "overlap_fraction_vs_v2 >= 0.80",
            "observed": round(float(m["overlap_fraction_vs_v2"]), 4),
            "pass": bool(float(m["overlap_fraction_vs_v2"]) >= 0.80),
        },
        {
            "criterion": "Direction balance (RA coverage)",
            "rule": "ra_count_v3 >= 8",
            "observed": int(m["ra_count_v3"]),
            "pass": bool(int(m["ra_count_v3"]) >= 8),
        },
        {
            "criterion": "Redundancy control",
            "rule": "jaccard_increase <= 0.005",
            "observed": round(jaccard_increase, 4),
            "pass": bool(jaccard_increase <= 0.005),
        },
        {
            "criterion": "No cohort sign flip (all-program delta)",
            "rule": "n_sign_flip == 0",
            "observed": int(sign_flip),
            "pass": bool(sign_flip == 0),
        },
        {
            "criterion": "Transferability (delta stability)",
            "rule": "mean_abs_delta_shift <= 0.02",
            "observed": round(mean_abs_delta_shift, 4),
            "pass": bool(mean_abs_delta_shift <= 0.02),
        },
        {
            "criterion": "Transferability (effect-size stability)",
            "rule": "mean_abs_g_shift <= 0.20",
            "observed": round(mean_abs_g_shift, 4),
            "pass": bool(mean_abs_g_shift <= 0.20),
        },
    ]
    out = pd.DataFrame(rows)
    out["overall_decision"] = np.where(out["pass"].all(), "upgrade-ready", "keep-as-exploratory")
    return out


def main() -> None:
    for p in [PANEL_V2, PANEL_V3, LONG_V2, LONG_V3]:
        if not p.exists():
            raise FileNotFoundError(f"Missing input: {p.relative_to(ROOT)}")

    v2 = pd.read_csv(PANEL_V2, sep="\t")
    v3 = pd.read_csv(PANEL_V3, sep="\t")
    long2 = pd.read_csv(LONG_V2, sep="\t")
    long3 = pd.read_csv(LONG_V3, sep="\t")

    gsets = {}
    gsets.update(load_gmt(GMT_HALL, "HALLMARK"))
    gsets.update(load_gmt(GMT_REAC, "REACTOME"))

    set2 = set(v2["pathway"].tolist())
    set3 = set(v3["pathway"].tolist())
    overlap = len(set2 & set3)

    metrics = pd.DataFrame(
        [
            {
                "panel_size_v2": len(set2),
                "panel_size_v3": len(set3),
                "overlap_n": overlap,
                "overlap_fraction_vs_v2": round(overlap / len(set2), 4),
                "overlap_fraction_vs_v3": round(overlap / len(set3), 4),
                "la_count_v2": int((v2["direction"] == "LA-biased").sum()),
                "ra_count_v2": int((v2["direction"] == "RA-biased").sum()),
                "la_count_v3": int((v3["direction"] == "LA-biased").sum()),
                "ra_count_v3": int((v3["direction"] == "RA-biased").sum()),
                "mean_pairwise_jaccard_v2": round(mean_pairwise_jaccard(v2["pathway"].tolist(), gsets), 4),
                "mean_pairwise_jaccard_v3": round(mean_pairwise_jaccard(v3["pathway"].tolist(), gsets), 4),
            }
        ]
    )
    metrics.to_csv(OUT_METRICS, sep="\t", index=False)

    axis2 = build_axis_effects(long2, "v2")
    axis3 = build_axis_effects(long3, "v3")
    axis = pd.concat([axis2, axis3], ignore_index=True)
    axis.to_csv(OUT_AXIS, sep="\t", index=False)
    decision = build_decision_table(metrics, axis)
    decision.to_csv(OUT_DECISION, sep="\t", index=False)

    # --- plotting ---
    set_nature_style()
    fig, axs = plt.subplots(2, 2, figsize=(COL_DOUBLE * 1.15, COL_DOUBLE * 0.9))
    ax1, ax2, ax3, ax4 = axs.ravel()

    # A. overlap + redundancy
    m = metrics.iloc[0]
    ax1.bar(
        ["overlap/v2", "overlap/v3"],
        [m["overlap_fraction_vs_v2"], m["overlap_fraction_vs_v3"]],
        color=["#4D4D4D", "#7F7F7F"],
    )
    ax1.set_ylim(0, 1)
    ax1.set_ylabel("fraction")
    ax1.set_title("A  Panel continuity")
    ax1.grid(axis="y", alpha=0.2)
    for i, val in enumerate([m["overlap_fraction_vs_v2"], m["overlap_fraction_vs_v3"]]):
        ax1.text(i, val + 0.03, f"{val:.2f}", ha="center", va="bottom", fontsize=9)

    # B. direction composition
    x = np.array([0, 1])
    la = np.array([m["la_count_v2"], m["la_count_v3"]], dtype=float)
    ra = np.array([m["ra_count_v2"], m["ra_count_v3"]], dtype=float)
    ax2.bar(x, la, color=LA_COLOR, label="LA-biased")
    ax2.bar(x, ra, bottom=la, color=RA_COLOR, label="RA-biased")
    ax2.set_xticks(x, ["v2", "v3"])
    ax2.set_ylabel("program count")
    ax2.set_title("B  Direction balance")
    ax2.legend(frameon=False, loc="upper right")

    # C. cohort-level all-program delta
    all_axis = axis[axis["axis"] == "all-programs"].copy()
    cohorts = ["GSE128188", "GSE79768", "GSE115574"]
    xv = np.arange(len(cohorts))
    d2 = all_axis[all_axis["version"] == "v2"].set_index("cohort").reindex(cohorts)["delta_af_minus_sr"].to_numpy(float)
    d3 = all_axis[all_axis["version"] == "v3"].set_index("cohort").reindex(cohorts)["delta_af_minus_sr"].to_numpy(float)
    ax3.plot(xv, d2, marker="o", color="#777777", label="v2")
    ax3.plot(xv, d3, marker="o", color="#111111", label="v3")
    ax3.axhline(0, color="#BBBBBB", linewidth=0.8)
    ax3.set_xticks(xv, cohorts, rotation=20)
    ax3.set_ylabel("AF-SR delta (donor axis aCAI)")
    ax3.set_title("C  Cohort-level overall rewiring shift")
    ax3.legend(frameon=False)

    # D. Hedges g by axis, v2 vs v3
    axes_order = ["LA-biased", "RA-biased"]
    ybase = np.arange(len(axes_order))
    width = 0.24
    for ci, cohort in enumerate(cohorts):
        sub = axis[axis["cohort"] == cohort].set_index(["version", "axis"])
        g2 = np.array([sub.loc[("v2", a), "hedges_g_af_minus_sr"] for a in axes_order], dtype=float)
        g3 = np.array([sub.loc[("v3", a), "hedges_g_af_minus_sr"] for a in axes_order], dtype=float)
        off = (ci - 1) * width
        ax4.scatter(g2, ybase + off, color="#999999", s=24, marker="o")
        ax4.scatter(g3, ybase + off, color="#111111", s=24, marker="x")
    ax4.axvline(0, color="#BBBBBB", linewidth=0.8)
    ax4.set_yticks(ybase, axes_order)
    ax4.set_xlabel("Hedges g (AF-SR)")
    ax4.set_title("D  Effect size by direction arm")
    ax4.text(0.98, 0.03, "o:v2  x:v3", transform=ax4.transAxes, ha="right", va="bottom", fontsize=8)
    ax4.text(
        0.02,
        0.03,
        f"decision: {decision['overall_decision'].iloc[0]}",
        transform=ax4.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        color="#333333",
        fontweight="bold",
    )

    for ax in axs.ravel():
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    fig.suptitle("Innovation benchmark: v2 vs v3 panel transferability", y=0.995, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT_FIG1D, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT_FIGS13, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    print(f"wrote {OUT_METRICS.relative_to(ROOT)}")
    print(f"wrote {OUT_AXIS.relative_to(ROOT)}")
    print(f"wrote {OUT_DECISION.relative_to(ROOT)}")
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")
    print(f"wrote {OUT_FIG1D.relative_to(ROOT)}")
    print(f"wrote {OUT_FIGS13.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
