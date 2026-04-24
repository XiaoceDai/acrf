#!/usr/bin/env python3
"""
Figure 1C (alternative A): gate funnel + unlabeled dual-evidence scatter + compact panel summary.

This rendering is designed to reduce visual overload from 30 pathway labels while
keeping the method logic explicit:
  - left: selection funnel (candidate -> dual-significant -> stability gate -> final K)
  - middle: NES vs GSVA logFC scatter without text labels
  - right: final panel composition summary
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label
from plot_style import COL_DOUBLE, FG_COLOR, LA_COLOR, LINE_COLOR, MUTED_COLOR, RA_COLOR, set_nature_style

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

IN_LONG = ROOT / "results" / "fig1c_pathway_selection_long.tsv"
IN_PANEL = ROOT / "results" / "fig1c_ncp_panel_v2.tsv"
OUT_PDF = ROOT / "figures" / "fig1c_dual_method_selection.pdf"
OUT_SUMMARY = ROOT / "results" / "fig1c_selection_summary.tsv"


def draw_funnel(ax: plt.Axes, counts: list[tuple[str, int]], n0: int) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    y_positions = [0.82, 0.62, 0.42, 0.22]
    widths = [0.86, 0.72, 0.58, 0.44]
    h = 0.12

    for i, ((label, n), y, w) in enumerate(zip(counts, y_positions, widths)):
        x = 0.5 - w / 2
        face = "#F8F8F8" if i < 3 else "#EFE9E2"
        edge = "#D9D9D9" if i < 3 else "#CDBEA8"
        ax.add_patch(
            mpatches.FancyBboxPatch(
                (x, y - h / 2),
                w,
                h,
                boxstyle="round,pad=0.01,rounding_size=0.015",
                facecolor=face,
                edgecolor=edge,
                linewidth=0.8,
            )
        )
        pct = 100.0 * n / n0 if n0 > 0 else 0.0
        ax.text(0.5, y + 0.015, label, ha="center", va="center", fontsize=7.1, color=MUTED_COLOR)
        ax.text(0.5, y - 0.018, f"n = {n} ({pct:.1f}%)", ha="center", va="center", fontsize=8.0, color=FG_COLOR, fontweight="bold")
        if i < len(counts) - 1:
            ax.annotate(
                "",
                xy=(0.5, y - h / 2 - 0.02),
                xytext=(0.5, y_positions[i + 1] + h / 2 + 0.01),
                arrowprops=dict(arrowstyle="->", lw=0.8, color=LINE_COLOR),
            )

    ax.set_title("selection gates", fontsize=8.4, fontweight="bold", loc="left", pad=2)


def draw_scatter(ax: plt.Axes, long_df: pd.DataFrame, final_mask: np.ndarray) -> None:
    x = long_df["NES"].to_numpy(float)
    y = long_df["gsva_logFC"].to_numpy(float)
    dirs = long_df["direction"].astype(str).to_numpy()

    ax.axhline(0, color=LINE_COLOR, lw=0.8, zorder=0)
    ax.axvline(0, color=LINE_COLOR, lw=0.8, zorder=0)

    # Background: all candidates
    ax.scatter(x, y, s=20, color="#D0D0D0", edgecolor="white", linewidth=0.25, zorder=1)

    # Gate-passing points by direction color
    dual_mask = (
        (long_df["padj"].to_numpy(float) < 0.01)
        & (long_df["gsva_padj"].to_numpy(float) < 0.01)
        & (np.sign(x) == np.sign(y))
    )
    la_mask = dual_mask & (dirs == "LA-biased")
    ra_mask = dual_mask & (dirs == "RA-biased")
    ax.scatter(x[la_mask], y[la_mask], s=26, color=LA_COLOR, edgecolor="white", linewidth=0.25, zorder=2)
    ax.scatter(x[ra_mask], y[ra_mask], s=26, color=RA_COLOR, edgecolor="white", linewidth=0.25, zorder=2)

    # Final panel ring
    ax.scatter(
        x[final_mask],
        y[final_mask],
        s=52,
        facecolor="none",
        edgecolor="#6B5A44",
        linewidth=0.8,
        zorder=3,
    )

    lim_x = float(max(0.8, np.nanmax(np.abs(x)) * 1.10))
    lim_y = float(max(0.4, np.nanmax(np.abs(y)) * 1.10))
    ax.set_xlim(-lim_x, lim_x)
    ax.set_ylim(-lim_y, lim_y)
    ax.set_xlabel("signed GSEA NES", fontsize=7.4)
    ax.set_ylabel("GSVA paired logFC (LA - RA)", fontsize=7.4)
    ax.tick_params(labelsize=7.0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("dual-evidence consistency (unlabeled)", fontsize=8.4, fontweight="bold", loc="left", pad=2)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#D0D0D0", markeredgecolor="white", markersize=4.6, label="all candidates"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=LA_COLOR, markeredgecolor="white", markersize=4.8, label="dual-significant LA"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=RA_COLOR, markeredgecolor="white", markersize=4.8, label="dual-significant RA"),
        Line2D([0], [0], marker="o", color="#6B5A44", markerfacecolor="none", markersize=5.2, label="final panel (K=30)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=6.2, frameon=False, handletextpad=0.4, borderpad=0.2)


def draw_summary(ax: plt.Axes, panel_df: pd.DataFrame) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("final panel summary", fontsize=8.4, fontweight="bold", loc="left", pad=2)

    k = len(panel_df)
    n_la = int((panel_df["direction"] == "LA-biased").sum())
    n_ra = int((panel_df["direction"] == "RA-biased").sum())
    med_stab = float(panel_df["stability_consistent_runs"].median()) if "stability_consistent_runs" in panel_df.columns else np.nan
    med_abs_nes = float(np.median(np.abs(panel_df["NES"].to_numpy(float))))

    # Simple stacked bar for LA/RA composition
    bar_x, bar_y, bar_w, bar_h = 0.05, 0.81, 0.90, 0.07
    frac_la = n_la / k if k > 0 else 0
    ax.add_patch(mpatches.Rectangle((bar_x, bar_y), bar_w * frac_la, bar_h, facecolor=LA_COLOR, edgecolor="none"))
    ax.add_patch(mpatches.Rectangle((bar_x + bar_w * frac_la, bar_y), bar_w * (1 - frac_la), bar_h, facecolor=RA_COLOR, edgecolor="none"))
    ax.add_patch(mpatches.Rectangle((bar_x, bar_y), bar_w, bar_h, facecolor="none", edgecolor="#BBBBBB", linewidth=0.6))
    ax.text(bar_x, bar_y + 0.10, f"K={k} (LA={n_la}, RA={n_ra})", fontsize=7.1, color=FG_COLOR, ha="left", va="center")

    ax.text(0.05, 0.72, f"median stability: {med_stab:.1f} / 8", fontsize=6.9, color=MUTED_COLOR, ha="left")
    ax.text(0.05, 0.66, f"median |NES|: {med_abs_nes:.2f}", fontsize=6.9, color=MUTED_COLOR, ha="left")

    # Top representatives (no full 30 list)
    top = panel_df.copy()
    top["abs_nes"] = top["NES"].abs()
    top = top.sort_values("abs_nes", ascending=False).head(6)
    ax.text(0.05, 0.57, "representative pathways:", fontsize=6.9, color=MUTED_COLOR, ha="left")
    y0 = 0.52
    for i, (_, r) in enumerate(top.iterrows()):
        yy = y0 - i * 0.075
        col = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        ax.text(0.05, yy, f"- {short_pathway_label(r['pathway'], 30)}", fontsize=6.8, color=col, ha="left")


def main() -> None:
    long_df = pd.read_csv(IN_LONG, sep="\t")
    panel_df = pd.read_csv(IN_PANEL, sep="\t")

    final_pathways = set(panel_df["pathway"].astype(str).tolist())
    final_mask = long_df["pathway"].astype(str).isin(final_pathways).to_numpy()

    gate_dual = (
        (long_df["padj"].to_numpy(float) < 0.01)
        & (long_df["gsva_padj"].to_numpy(float) < 0.01)
        & (np.sign(long_df["NES"].to_numpy(float)) == np.sign(long_df["gsva_logFC"].to_numpy(float)))
    )
    gate_stab = gate_dual & (long_df["stability_consistent_runs"].to_numpy(int) >= 7)

    counts = [
        ("candidate pathways", int(len(long_df))),
        ("dual-significant", int(gate_dual.sum())),
        ("+ stability >= 7/8", int(gate_stab.sum())),
        ("final panel", int(len(panel_df))),
    ]

    fig = plt.figure(figsize=(COL_DOUBLE * 1.05, COL_DOUBLE * 0.62))
    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[1.2, 2.4, 1.25],
        wspace=0.20,
        top=0.90,
        bottom=0.12,
        left=0.04,
        right=0.99,
    )

    ax_f = fig.add_subplot(gs[0, 0])
    ax_s = fig.add_subplot(gs[0, 1])
    ax_r = fig.add_subplot(gs[0, 2])

    draw_funnel(ax_f, counts, n0=counts[0][1])
    draw_scatter(ax_s, long_df, final_mask)
    draw_summary(ax_r, panel_df)

    fig.text(
        0.04,
        0.95,
        "c  Internal baseline implementation (v2): gate logic and final panel definition",
        fontsize=9.8,
        fontweight="bold",
        ha="left",
        va="center",
        color=FG_COLOR,
    )

    fig.text(
        0.04,
        0.05,
        "No pathway labels are overlaid in the scatter to preserve readability; full member list is retained in results tables.",
        fontsize=6.5,
        color=MUTED_COLOR,
        ha="left",
        va="center",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)

    pd.DataFrame(
        {
            "stage": [s for s, _ in counts],
            "n": [n for _, n in counts],
        }
    ).to_csv(OUT_SUMMARY, sep="\t", index=False)

    print(f"wrote {OUT_PDF.relative_to(ROOT)}")
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
