#!/usr/bin/env python3
"""
Figure 1D (healthy-only): v2 -> v3 upgrade benchmark and gate audit.

This panel is intentionally restricted to baseline construction evidence:
  - panel continuity
  - direction balance
  - redundancy control
  - bootstrap robustness (same resampling universe)

It does not use AF/SR disease-stage comparisons.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_style import (
    COL_DOUBLE,
    FG_COLOR,
    HIGHLIGHT_COLOR,
    LA_COLOR,
    LINE_COLOR,
    MUTED_COLOR,
    RA_COLOR,
    SUPPORT_COLOR,
    set_nature_style,
    enforce_font_floor, FS_PANEL_TITLE, FS_TICK, FS_AXIS_LABEL, FS_FOOTNOTE, FS_LEGEND,
)

ROOT = Path(__file__).resolve().parent.parent

PANEL_V2 = ROOT / "results" / "fig1c_ncp_panel_v2.tsv"
PANEL_V3 = ROOT / "results" / "innovation_ncp_v3_panel.tsv"
BOOT_FREQ = ROOT / "results" / "innovation_ncp_bootstrap_frequency.tsv"

GMT_HALL = ROOT / "bulk" / "Project_P01" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC = ROOT / "bulk" / "Project_P01" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"

OUT_METRICS = ROOT / "results" / "fig1d_healthy_upgrade_metrics.tsv"
OUT_GATES = ROOT / "results" / "fig1d_healthy_upgrade_gates.tsv"
OUT_MEMBERSHIP = ROOT / "results" / "fig1d_healthy_upgrade_membership.tsv"
OUT_PDF = ROOT / "figures" / "fig1d_prototype_vs_ncp_benchmark.pdf"


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
    uni = a | b
    if not uni:
        return 0.0
    return len(a & b) / len(uni)


def mean_pairwise_jaccard(paths: list[str], gsets: dict[str, set[str]]) -> float:
    vals: list[float] = []
    for i in range(len(paths)):
        for j in range(i + 1, len(paths)):
            vals.append(jaccard(gsets.get(paths[i], set()), gsets.get(paths[j], set())))
    if not vals:
        return float("nan")
    return float(np.mean(vals))


def bootstrap_panel_stats(pathways: set[str], boot_df: pd.DataFrame) -> dict[str, float]:
    sub = boot_df[boot_df["pathway"].isin(pathways)].copy()
    freq = sub["inclusion_freq"].to_numpy(float)
    return {
        "n_mapped": int(len(sub)),
        "median_freq": float(np.median(freq)) if len(freq) else np.nan,
        "q1_freq": float(np.quantile(freq, 0.25)) if len(freq) else np.nan,
        "q3_freq": float(np.quantile(freq, 0.75)) if len(freq) else np.nan,
        "share_ge_0_8": float(np.mean(freq >= 0.8)) if len(freq) else np.nan,
        "share_ge_0_5": float(np.mean(freq >= 0.5)) if len(freq) else np.nan,
    }


def main() -> None:
    for p in [PANEL_V2, PANEL_V3, BOOT_FREQ, GMT_HALL, GMT_REAC]:
        if not p.exists():
            raise FileNotFoundError(f"Missing input: {p.relative_to(ROOT)}")

    v2 = pd.read_csv(PANEL_V2, sep="\t")
    v3 = pd.read_csv(PANEL_V3, sep="\t")
    boot = pd.read_csv(BOOT_FREQ, sep="\t")

    set2 = set(v2["pathway"].tolist())
    set3 = set(v3["pathway"].tolist())
    shared = sorted(set2 & set3)
    v2_only = sorted(set2 - set3)
    v3_only = sorted(set3 - set2)

    gsets = {}
    gsets.update(load_gmt(GMT_HALL, "HALLMARK"))
    gsets.update(load_gmt(GMT_REAC, "REACTOME"))

    stats2 = bootstrap_panel_stats(set2, boot)
    stats3 = bootstrap_panel_stats(set3, boot)

    mean_j2 = mean_pairwise_jaccard(v2["pathway"].tolist(), gsets)
    mean_j3 = mean_pairwise_jaccard(v3["pathway"].tolist(), gsets)

    metrics = pd.DataFrame(
        [
            {
                "panel_size_v2": len(set2),
                "panel_size_v3": len(set3),
                "shared_n": len(shared),
                "v2_only_n": len(v2_only),
                "v3_only_n": len(v3_only),
                "overlap_fraction_vs_v2": len(shared) / len(set2),
                "overlap_fraction_vs_v3": len(shared) / len(set3),
                "la_count_v2": int((v2["direction"] == "LA-biased").sum()),
                "ra_count_v2": int((v2["direction"] == "RA-biased").sum()),
                "la_count_v3": int((v3["direction"] == "LA-biased").sum()),
                "ra_count_v3": int((v3["direction"] == "RA-biased").sum()),
                "ra_gain": int((v3["direction"] == "RA-biased").sum() - (v2["direction"] == "RA-biased").sum()),
                "mean_pairwise_jaccard_v2": mean_j2,
                "mean_pairwise_jaccard_v3": mean_j3,
                "jaccard_delta_v3_minus_v2": mean_j3 - mean_j2,
                "median_inclusion_freq_v2set": stats2["median_freq"],
                "median_inclusion_freq_v3set": stats3["median_freq"],
                "share_ge_0_8_v2set": stats2["share_ge_0_8"],
                "share_ge_0_8_v3set": stats3["share_ge_0_8"],
                "share_ge_0_5_v2set": stats2["share_ge_0_5"],
                "share_ge_0_5_v3set": stats3["share_ge_0_5"],
            }
        ]
    )
    metrics.to_csv(OUT_METRICS, sep="\t", index=False, float_format="%.6f")

    m = metrics.iloc[0]
    gates = pd.DataFrame(
        [
            {
                "criterion": "Panel continuity",
                "rule": "overlap_fraction_vs_v2 >= 0.80",
                "observed": float(m["overlap_fraction_vs_v2"]),
                "pass": bool(float(m["overlap_fraction_vs_v2"]) >= 0.80),
            },
            {
                "criterion": "Direction balance",
                "rule": "ra_count_v3 >= 8 and ra_gain >= 2",
                "observed": f"RA {int(m['ra_count_v3'])} (gain {int(m['ra_gain']):+d})",
                "pass": bool(int(m["ra_count_v3"]) >= 8 and int(m["ra_gain"]) >= 2),
            },
            {
                "criterion": "Redundancy control",
                "rule": "mean_pairwise_jaccard_v3 <= 0.025 and delta <= 0.005",
                "observed": f"{float(m['mean_pairwise_jaccard_v3']):.4f} (delta {float(m['jaccard_delta_v3_minus_v2']):+.4f})",
                "pass": bool(float(m["mean_pairwise_jaccard_v3"]) <= 0.025 and float(m["jaccard_delta_v3_minus_v2"]) <= 0.005),
            },
            {
                "criterion": "Bootstrap robustness",
                "rule": "median inclusion freq >= 0.99",
                "observed": float(m["median_inclusion_freq_v3set"]),
                "pass": bool(float(m["median_inclusion_freq_v3set"]) >= 0.99),
            },
            {
                "criterion": "High-stability share",
                "rule": "share(inclusion_freq >= 0.8) >= 0.70",
                "observed": float(m["share_ge_0_8_v3set"]),
                "pass": bool(float(m["share_ge_0_8_v3set"]) >= 0.70),
            },
            {
                "criterion": "Upgrade gain",
                "rule": "v3 median and >=0.8 share >= v2",
                "observed": f"median {float(m['median_inclusion_freq_v2set']):.3f}->{float(m['median_inclusion_freq_v3set']):.3f}; share {float(m['share_ge_0_8_v2set']):.3f}->{float(m['share_ge_0_8_v3set']):.3f}",
                "pass": bool(
                    float(m["median_inclusion_freq_v3set"]) >= float(m["median_inclusion_freq_v2set"])
                    and float(m["share_ge_0_8_v3set"]) >= float(m["share_ge_0_8_v2set"])
                ),
            },
        ]
    )
    gates["overall_decision"] = np.where(gates["pass"].all(), "upgrade-ready", "hold-and-iterate")
    gates.to_csv(OUT_GATES, sep="\t", index=False)

    membership_rows: list[dict[str, str]] = []
    for p in shared:
        membership_rows.append({"pathway": p, "membership": "shared"})
    for p in v2_only:
        membership_rows.append({"pathway": p, "membership": "v2-only"})
    for p in v3_only:
        membership_rows.append({"pathway": p, "membership": "v3-only"})
    pd.DataFrame(membership_rows).sort_values(["membership", "pathway"]).to_csv(
        OUT_MEMBERSHIP, sep="\t", index=False
    )

    # ---------------- plot ----------------
    set_nature_style()
    fig, axs = plt.subplots(2, 2, figsize=(COL_DOUBLE * 1.14, COL_DOUBLE * 1.00))
    ax_a, ax_b, ax_c, ax_d = axs.ravel()

    # a) continuity with controlled replacement
    share_vals = np.array([len(shared), len(v3_only)], dtype=float) / len(set3)
    ax_a.bar(["shared", "NCP-only"], share_vals, color=[MUTED_COLOR, HIGHLIGHT_COLOR], width=0.62)
    ax_a.axhline(0.80, color=LINE_COLOR, linestyle="--", linewidth=0.9)
    ax_a.text(
        0.02,
        0.95,
        "gate: overlap >= 0.80",
        transform=ax_a.transAxes,
        ha="left",
        va="top",
        fontsize=enforce_font_floor(FS_TICK),
        color=MUTED_COLOR,
    )
    ax_a.text(0, share_vals[0] + 0.03, f"{len(shared)}/{len(set3)}", ha="center", va="bottom", fontsize=enforce_font_floor(FS_TICK + 1.2))
    ax_a.text(1, share_vals[1] + 0.03, f"{len(v3_only)}/{len(set3)}", ha="center", va="bottom", fontsize=enforce_font_floor(FS_TICK + 1.2))
    ax_a.set_ylim(0, 1.0)
    ax_a.set_ylabel("fraction of NCP panel", fontsize=enforce_font_floor(FS_AXIS_LABEL))
    ax_a.set_title("a  Panel continuity with controlled replacement", loc="left", fontsize=enforce_font_floor(FS_PANEL_TITLE + 1.5), fontweight="bold")
    ax_a.grid(axis="y", alpha=0.16)

    # b) direction balance (inside chamber-paired family)
    x = np.array([0, 1], dtype=float)
    la = np.array([m["la_count_v2"], m["la_count_v3"]], dtype=float)
    ra = np.array([m["ra_count_v2"], m["ra_count_v3"]], dtype=float)
    ax_b.bar(x, la, width=0.62, color=LA_COLOR, label="LA-biased")
    ax_b.bar(x, ra, width=0.62, bottom=la, color=RA_COLOR, label="RA-biased")
    ax_b.set_xticks(x, ["prototype", "NCP"])
    ax_b.set_ylim(0, max(la + ra) + 2)
    ax_b.set_ylabel("program count", fontsize=enforce_font_floor(FS_AXIS_LABEL))
    ax_b.set_title("b  Direction rebalance in paired framework", loc="left", fontsize=enforce_font_floor(FS_PANEL_TITLE + 1.5), fontweight="bold")
    ax_b.text(
        0.98,
        0.93,
        f"RA gain: {int(m['ra_gain']):+d}",
        transform=ax_b.transAxes,
        ha="right",
        va="top",
        fontsize=enforce_font_floor(FS_TICK + 1),
        color=FG_COLOR,
        fontweight="bold",
    )
    ax_b.legend(frameon=False, loc="upper left", fontsize=enforce_font_floor(FS_LEGEND))

    # c) bootstrap robustness distribution + redundancy note
    freq2 = boot.loc[boot["pathway"].isin(set2), "inclusion_freq"].to_numpy(float)
    freq3 = boot.loc[boot["pathway"].isin(set3), "inclusion_freq"].to_numpy(float)
    bp = ax_c.boxplot(
        [freq2, freq3],
        positions=[0, 1],
        widths=0.42,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": FG_COLOR, "linewidth": 1.2},
        boxprops={"linewidth": 0.8},
        whiskerprops={"linewidth": 0.8},
        capprops={"linewidth": 0.8},
    )
    bp["boxes"][0].set_facecolor("#DADADA")
    bp["boxes"][1].set_facecolor(SUPPORT_COLOR)
    rng = np.random.default_rng(20260418)
    for i, arr in enumerate([freq2, freq3]):
        jitter = rng.uniform(-0.07, 0.07, size=len(arr))
        ax_c.scatter(np.full_like(arr, i, dtype=float) + jitter, arr, s=12, color=FG_COLOR, alpha=0.36, linewidths=0)
    ax_c.axhline(0.8, color=LINE_COLOR, linestyle="--", linewidth=0.9)
    ax_c.set_xticks([0, 1], ["prototype set", "NCP set"])
    ax_c.set_ylim(0, 1.04)
    ax_c.set_ylabel("bootstrap inclusion frequency", fontsize=enforce_font_floor(FS_AXIS_LABEL))
    ax_c.set_title("c  Healthy resampling robustness", loc="left", fontsize=enforce_font_floor(FS_PANEL_TITLE + 1.5), fontweight="bold")
    ax_c.text(
        0.02,
        0.96,
        f"median: {np.median(freq2):.3f} -> {np.median(freq3):.3f}\n"
        f"share >=0.8: {np.mean(freq2 >= 0.8):.3f} -> {np.mean(freq3 >= 0.8):.3f}\n"
        f"mean pairwise Jaccard: {mean_j2:.4f} -> {mean_j3:.4f}",
        transform=ax_c.transAxes,
        ha="left",
        va="top",
        fontsize=enforce_font_floor(FS_TICK),
        color=MUTED_COLOR,
    )

    # d) gate audit summary
    ax_d.set_axis_off()
    ax_d.set_title("d  Healthy-only upgrade gate audit", loc="left", fontsize=enforce_font_floor(FS_PANEL_TITLE + 1.5), fontweight="bold")
    y = 0.95
    row_h = 0.125
    for _, row in gates.iterrows():
        y0 = y - row_h + 0.01
        face = "#EAF4EC" if bool(row["pass"]) else "#F7EAEA"
        rect = mpatches.FancyBboxPatch(
            (0.01, y0),
            0.98,
            row_h - 0.02,
            boxstyle="round,pad=0.006,rounding_size=0.012",
            transform=ax_d.transAxes,
            facecolor=face,
            edgecolor=LINE_COLOR,
            linewidth=0.7,
        )
        ax_d.add_patch(rect)
        ax_d.text(
            0.03,
            y0 + (row_h - 0.02) * 0.63,
            str(row["criterion"]),
            transform=ax_d.transAxes,
            ha="left",
            va="center",
            fontsize=enforce_font_floor(FS_TICK + 0.5),
            color=FG_COLOR,
            fontweight="bold",
        )
        ax_d.text(
            0.03,
            y0 + (row_h - 0.02) * 0.24,
            f"rule: {row['rule']}",
            transform=ax_d.transAxes,
            ha="left",
            va="center",
            fontsize=enforce_font_floor(FS_TICK),
            color=MUTED_COLOR,
        )
        ax_d.text(
            0.74,
            y0 + (row_h - 0.02) * 0.5,
            f"{row['observed']}",
            transform=ax_d.transAxes,
            ha="left",
            va="center",
            fontsize=enforce_font_floor(FS_TICK),
            color=FG_COLOR,
        )
        ax_d.text(
            0.97,
            y0 + (row_h - 0.02) * 0.5,
            "PASS" if bool(row["pass"]) else "FAIL",
            transform=ax_d.transAxes,
            ha="right",
            va="center",
            fontsize=enforce_font_floor(FS_TICK + 0.5),
            color=SUPPORT_COLOR if bool(row["pass"]) else HIGHLIGHT_COLOR,
            fontweight="bold",
        )
        y -= row_h
    ax_d.text(
        0.01,
        0.02,
        f"overall decision: {gates['overall_decision'].iloc[0]}",
        transform=ax_d.transAxes,
        ha="left",
        va="bottom",
        fontsize=enforce_font_floor(FS_TICK + 0.5),
        color=FG_COLOR,
        fontweight="bold",
    )

    fig.suptitle(
        "Fig 1c  Prototype-to-NCP upgrade benchmark (healthy-only, paired framework)",
        y=0.995,
        fontsize=enforce_font_floor(FS_PANEL_TITLE + 2.5),
        fontweight="bold",
    )
    fig.text(
        0.012,
        0.008,
        "Metrics use only baseline panel tables, gene-set overlap, and healthy bootstrap resampling.",
        fontsize=enforce_font_floor(FS_FOOTNOTE),
        color=MUTED_COLOR,
        ha="left",
        va="bottom",
    )

    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    print(f"wrote {OUT_METRICS.relative_to(ROOT)}")
    print(f"wrote {OUT_GATES.relative_to(ROOT)}")
    print(f"wrote {OUT_MEMBERSHIP.relative_to(ROOT)}")
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
