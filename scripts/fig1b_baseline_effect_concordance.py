#!/usr/bin/env python3
"""Final Figure 1B: healthy reference magnitude and donor concordance."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label  # noqa: E402
from plot_style import (  # noqa: E402
    COL_DOUBLE,
    FG_COLOR,
    LA_COLOR,
    LINE_COLOR,
    MUTED_COLOR,
    RA_COLOR,
    SUPPORT_COLOR,
    set_nature_style,
)

SOURCE = ROOT / "results" / "fig1b_baseline_donor_aCAI.tsv"
FAMILY_MAP = ROOT / "results" / "innovation_ncp_v3_axis_families.tsv"
FAMILY_LABELS = ROOT / "results" / "ncp_v3_family_labels.tsv"
OUT_BASE = ROOT / "figures" / "fig1b_healthy_reference_concordance"
OUT_SOURCE = ROOT / "results" / "fig1b_healthy_reference_concordance_source.tsv"

set_nature_style()
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def build_source() -> tuple[pd.DataFrame, list[str], dict[str, float | int]]:
    raw = pd.read_csv(SOURCE, sep="\t")
    fam = pd.read_csv(FAMILY_MAP, sep="\t")[["pathway", "family_id"]]
    labels = pd.read_csv(FAMILY_LABELS, sep="\t")[["family_id", "short_label", "axis_group"]]
    donor_cols = [c for c in raw.columns if c.startswith("aCAI_")]
    if len(donor_cols) != 8:
        raise AssertionError(f"Expected 8 donor columns, found {len(donor_cols)}")

    df = raw.merge(fam, on="pathway", how="left", validate="one_to_one")
    df = df.merge(labels, on="family_id", how="left", validate="many_to_one")
    if df[["family_id", "short_label", "axis_group"]].isna().any().any():
        raise AssertionError("Missing family metadata in Figure 1B join")

    direction_sign = df["direction"].map({"LA-biased": 1.0, "RA-biased": -1.0})
    if direction_sign.isna().any():
        raise AssertionError("Unexpected chamber direction")

    aligned_cols: list[str] = []
    for c in donor_cols:
        new = c.replace("aCAI_", "aligned_aCAI_")
        df[new] = df[c].astype(float) * direction_sign
        aligned_cols.append(new)
    df["M"] = df[aligned_cols].abs().mean(axis=1)
    df["directional_concordance_k"] = (df[aligned_cols] > 0).sum(axis=1).astype(int)
    df["display_label"] = df["pathway"].map(lambda p: short_pathway_label(p, 35))
    df["direction_rank"] = df["direction"].map({"LA-biased": 0, "RA-biased": 1}).astype(int)
    df = df.sort_values(["direction_rank", "family_id", "program_order"]).reset_index(drop=True)
    df["display_order"] = np.arange(1, len(df) + 1)

    summary: dict[str, float | int] = {
        "n_programs": int(len(df)),
        "n_la": int((df["direction"] == "LA-biased").sum()),
        "n_ra": int((df["direction"] == "RA-biased").sum()),
        "n_full_concordance": int((df["directional_concordance_k"] == 8).sum()),
        "n_ge6_concordance": int((df["directional_concordance_k"] >= 6).sum()),
        "median_M": float(df["M"].median()),
        "n_families": int(df["family_id"].nunique()),
    }
    expected = {
        "n_programs": 30,
        "n_la": 22,
        "n_ra": 8,
        "n_full_concordance": 21,
        "n_ge6_concordance": 28,
        "n_families": 16,
    }
    for key, value in expected.items():
        if summary[key] != value:
            raise AssertionError(f"{key}: expected {value}, observed {summary[key]}")

    source_cols = [
        "display_order",
        "program_order",
        "pathway",
        "display_label",
        "direction",
        "family_id",
        "short_label",
        "M",
        "directional_concordance_k",
    ] + donor_cols + aligned_cols
    OUT_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    df[source_cols].to_csv(OUT_SOURCE, sep="\t", index=False, float_format="%.10g")
    return df, aligned_cols, summary


def render(df: pd.DataFrame, aligned_cols: list[str], summary: dict[str, float | int]) -> None:
    n = len(df)
    y = np.arange(n)
    fig = plt.figure(figsize=(COL_DOUBLE * 1.32, COL_DOUBLE * 1.22))
    gs = fig.add_gridspec(
        1,
        6,
        width_ratios=[0.13, 0.72, 2.65, 1.65, 2.05, 0.62],
        left=0.045,
        right=0.985,
        bottom=0.105,
        top=0.835,
        wspace=0.025,
    )
    ax_dir = fig.add_subplot(gs[0, 0])
    ax_fam = fig.add_subplot(gs[0, 1], sharey=ax_dir)
    ax_lbl = fig.add_subplot(gs[0, 2], sharey=ax_dir)
    ax_don = fig.add_subplot(gs[0, 3], sharey=ax_dir)
    ax_mag = fig.add_subplot(gs[0, 4], sharey=ax_dir)
    ax_k = fig.add_subplot(gs[0, 5], sharey=ax_dir)

    fig.text(0.045, 0.965, "b", fontsize=13, fontweight="bold", ha="left", va="top")
    fig.text(
        0.070,
        0.965,
        "Healthy reference magnitude and donor concordance for 30 deployed programs",
        fontsize=10.8,
        fontweight="bold",
        ha="left",
        va="top",
    )
    summary_text = (
        f"{summary['n_la']} LA / {summary['n_ra']} RA  |  "
        f"8/8 concordant: {summary['n_full_concordance']}/{summary['n_programs']}  |  "
        f"≥6/8 concordant: {summary['n_ge6_concordance']}/{summary['n_programs']}  |  "
        f"median M[p] = {summary['median_M']:.3f}"
    )
    fig.text(0.070, 0.915, summary_text, fontsize=8.3, color=FG_COLOR, ha="left", va="top")
    fig.text(
        0.070,
        0.882,
        "Rows are grouped by canonical chamber direction, then family ID and frozen program order.",
        fontsize=7.4,
        color=MUTED_COLOR,
        ha="left",
        va="top",
    )

    family_starts: list[int] = []
    last = None
    for i, r in df.iterrows():
        key = (r["direction"], int(r["family_id"]))
        if key != last:
            family_starts.append(i)
            last = key
    boundaries = [i - 0.5 for i in family_starts[1:]]

    for ax in (ax_dir, ax_fam, ax_lbl, ax_don, ax_mag, ax_k):
        ax.set_ylim(n - 0.5, -0.5)
        for boundary in boundaries:
            ax.axhline(boundary, color=LINE_COLOR, lw=0.6, zorder=0)

    ax_dir.set_xlim(0, 1)
    ax_dir.axis("off")
    for i, direction in enumerate(df["direction"]):
        ax_dir.add_patch(
            patches.Rectangle(
                (0.15, i - 0.47),
                0.70,
                0.94,
                facecolor=LA_COLOR if direction == "LA-biased" else RA_COLOR,
                edgecolor="white",
                linewidth=0.25,
            )
        )
    ax_dir.set_title("dir.", fontsize=6.7, fontweight="bold", pad=7)

    ax_fam.set_xlim(0, 1)
    ax_fam.axis("off")
    for start in family_starts:
        fid = int(df.loc[start, "family_id"])
        members = df.index[df["family_id"].eq(fid)].tolist()
        center = float(np.mean(members))
        ax_fam.text(0.94, center, f"F{fid}", fontsize=7.1, fontweight="bold", ha="right", va="center")
    ax_fam.set_title("family", fontsize=6.7, fontweight="bold", pad=7)

    ax_lbl.set_xlim(0, 1)
    ax_lbl.axis("off")
    for i, r in df.iterrows():
        ax_lbl.text(0.98, i, r["display_label"], fontsize=6.7, color=FG_COLOR, ha="right", va="center")
    ax_lbl.set_title("deployed program", fontsize=6.9, fontweight="bold", pad=7)

    ax_don.set_xlim(-0.55, 7.55)
    ax_don.set_facecolor("#FAFAFA")
    for i, r in df.iterrows():
        vals = r[aligned_cols].to_numpy(dtype=float)
        concordant = vals > 0
        ax_don.scatter(
            np.arange(8)[concordant],
            np.full(int(concordant.sum()), i),
            s=17,
            color=SUPPORT_COLOR,
            edgecolor="white",
            linewidth=0.3,
            zorder=3,
        )
        ax_don.scatter(
            np.arange(8)[~concordant],
            np.full(int((~concordant).sum()), i),
            s=17,
            marker="x",
            color=MUTED_COLOR,
            linewidth=0.8,
            zorder=4,
        )
    ax_don.set_xticks(np.arange(8), [f"d{i}" for i in range(1, 9)], fontsize=6.4)
    ax_don.set_yticks([])
    ax_don.set_title("donor\nconcordance", fontsize=6.9, fontweight="bold", pad=5)
    for side in ("top", "right", "left"):
        ax_don.spines[side].set_visible(False)
    ax_don.spines["bottom"].set_color(LINE_COLOR)

    xmax = float(df["M"].max()) * 1.12
    ax_mag.set_xlim(0, xmax)
    for i, r in df.iterrows():
        color = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        ax_mag.plot([0, r["M"]], [i, i], color=LINE_COLOR, lw=1.0, zorder=1)
        ax_mag.scatter(r["M"], i, s=21, color=color, edgecolor="white", linewidth=0.35, zorder=2)
    ax_mag.set_yticks([])
    ax_mag.set_xlabel(r"$M[p]=mean_d(|aCAI[d,p]|)$", fontsize=7.2, labelpad=5)
    ax_mag.set_title("healthy reference\nmagnitude", fontsize=6.9, fontweight="bold", pad=5)
    ax_mag.tick_params(axis="x", labelsize=6.5)
    for side in ("top", "right", "left"):
        ax_mag.spines[side].set_visible(False)
    ax_mag.spines["bottom"].set_color(LINE_COLOR)

    ax_k.set_xlim(0, 1)
    ax_k.axis("off")
    for i, k in enumerate(df["directional_concordance_k"]):
        ax_k.text(
            0.5,
            i,
            f"{int(k)}/8",
            fontsize=7.0,
            color=SUPPORT_COLOR if int(k) == 8 else MUTED_COLOR,
            fontweight="bold" if int(k) == 8 else "normal",
            ha="center",
            va="center",
        )
    ax_k.set_title("k/8", fontsize=7.2, fontweight="bold", pad=7)

    handles = [
        mlines.Line2D([], [], marker="o", linestyle="none", color=SUPPORT_COLOR, markersize=4.5, label="canonical direction"),
        mlines.Line2D([], [], marker="x", linestyle="none", color=MUTED_COLOR, markersize=4.5, label="opposite direction"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=7.0, bbox_to_anchor=(0.5, 0.025))

    OUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BASE.with_suffix(".pdf"), bbox_inches=None, facecolor="white")
    fig.savefig(OUT_BASE.with_suffix(".svg"), bbox_inches=None, facecolor="white")
    fig.savefig(OUT_BASE.with_suffix(".png"), dpi=600, bbox_inches=None, facecolor="white")
    plt.close(fig)


def main() -> None:
    df, aligned_cols, summary = build_source()
    render(df, aligned_cols, summary)
    print(f"wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(
        f"summary: {summary['n_la']}/{summary['n_ra']} LA/RA; "
        f"{summary['n_full_concordance']}/30 at 8/8; "
        f"{summary['n_ge6_concordance']}/30 at >=6/8; "
        f"median M={summary['median_M']:.10f}"
    )


if __name__ == "__main__":
    main()
