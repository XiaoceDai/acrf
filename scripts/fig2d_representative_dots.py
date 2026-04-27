#!/usr/bin/env python3
"""
Figure 2D - donor-level aCAI_rel dot plots for representative programs.

Selection rule (representative donor dots; script header):
  For each class in {Amplified, Attenuated, Reversed, Preserved}:
    candidates = programs whose consensus_class == class or class-weak
    - class != Preserved:  pick argmin of min(q_perm_128188, q_perm_79768)
    - class == Preserved:  pick argmax of
                           max_c( min(mean_aCAI_rel_SR_c, mean_aCAI_rel_AF_c) )
    - empty -> skip that panel.

Panels are laid out left-to-right; y-axis is shared and each panel is
titled "<short_label> . <consensus_class>".

Within each panel columns are:
    SR_128188   AF_128188  |  SR_79768   AF_79768  |  SR_115574   AF_115574
with a visible gap between cohorts.  Dots are individual donors.
Horizontal marks are mean +- 95 % t-interval.
Reference lines: y = 0 (chamber bias absent) and y = +1 (Fig 1 healthy).
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from plot_style import (
    COL_DOUBLE, REWIRING_CLASS_COLORS,
    set_nature_style,
)

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "figures" / "fig2d_representative_dots.pdf"

COHORTS = ["GSE128188", "GSE79768", "GSE115574"]
COHORT_SHORT = {"GSE128188": "128188", "GSE79768": "79768", "GSE115574": "115574"}
# Cohort colour = muted; strict class colour is kept for panel title only.
COHORT_COLORS = {"GSE128188": "#4C72B0",  # blue-grey
                 "GSE79768":  "#8E6C8A",  # grey-purple
                 "GSE115574": "#666666"}  # neutral dark grey (replication)
COHORT_MARKERS = {"GSE128188": "o", "GSE79768": "^", "GSE115574": "D"}

CLASS_ORDER = ["Amplified", "Attenuated", "Reversed", "Preserved"]
CLASS_TO_WEAK = {
    "Amplified":  ("Amplified",  "Amplified-weak"),
    "Attenuated": ("Attenuated", "Attenuated-weak"),
    "Reversed":   ("Reversed",   "Reversed-weak"),
    "Preserved":  ("Preserved",),
}


def load_data():
    progs = pd.read_csv(ROOT / "results" / "af_chamber_rewiring_programs.tsv",
                        sep="\t")
    src   = pd.read_csv(ROOT / "results" / "figure2_source_data.tsv", sep="\t")
    return progs, src


def select_representatives(progs: pd.DataFrame) -> list[dict]:
    """Return a list of {class, program_id, program_label, consensus_class}."""
    picks: list[dict] = []
    for cls in CLASS_ORDER:
        accepted = CLASS_TO_WEAK[cls]
        cand = progs[progs["consensus_class"].isin(accepted)].copy()
        if cand.empty:
            continue
        if cls == "Preserved":
            cand["tiebreak"] = cand.apply(lambda r: max(
                min(r["mean_aCAI_rel_SR_GSE128188"], r["mean_aCAI_rel_AF_GSE128188"]),
                min(r["mean_aCAI_rel_SR_GSE79768"],  r["mean_aCAI_rel_AF_GSE79768"]),
                min(r["mean_aCAI_rel_SR_GSE115574"], r["mean_aCAI_rel_AF_GSE115574"]),
            ), axis=1)
            pick = cand.sort_values("tiebreak", ascending=False).iloc[0]
        else:
            cand["minq"] = cand[["q_perm_GSE128188", "q_perm_GSE79768"]].min(axis=1)
            pick = cand.sort_values("minq", ascending=True).iloc[0]
        picks.append(dict(
            panel_class=cls,
            program_id=pick["program_id"],
            program_label=pick["program_label"],
            consensus_class=pick["consensus_class"],
        ))
    return picks


def ci95(values: np.ndarray) -> tuple[float, float]:
    n = len(values)
    if n < 2:
        return (float(values.mean()) if n else float("nan"),
                float("nan"))
    m = float(values.mean())
    se = float(values.std(ddof=1) / math.sqrt(n))
    half = float(stats.t.ppf(0.975, df=n - 1) * se)
    return m, half


def perm_marker(p: float) -> str:
    if math.isnan(p):
        return ""
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "."
    return ""


def plot_panel(ax, pick: dict, src: pd.DataFrame, y_lim: tuple[float, float],
               prog_row: pd.Series):
    """One panel for a representative program."""
    df = src[src["program_id"] == pick["program_id"]]

    x_pos = {}   # (cohort, rhythm) -> x
    base = 0
    centres = []
    for c in COHORTS:
        x_pos[(c, "SR")] = base + 0.0
        x_pos[(c, "AF")] = base + 1.0
        centres.append((c, base + 0.5))
        base += 2 + 0.8                 # cohort gap

    # Reference horizontal lines.
    ax.axhline(0.0, color="#888888", lw=0.5, ls="--", zorder=1)
    ax.axhline(1.0, color="#444444", lw=0.6, ls=":", zorder=1)

    # Points + mean / CI per group.
    rng = np.random.default_rng(42)
    for (c, rhy), x in x_pos.items():
        sub = df[(df["cohort"] == c) & (df["rhythm"] == rhy)]
        vals = sub["aCAI_rel"].to_numpy()
        if len(vals) == 0:
            continue
        jitter = rng.uniform(-0.16, 0.16, size=len(vals))
        ax.scatter(x + jitter, vals,
                   marker=COHORT_MARKERS[c],
                   s=12, facecolor=COHORT_COLORS[c],
                   edgecolor="white", linewidth=0.4,
                   alpha=0.85, zorder=3)
        m, half = ci95(vals)
        # mean bar
        ax.plot([x - 0.30, x + 0.30], [m, m],
                color=COHORT_COLORS[c], lw=1.0, zorder=4)
        if not math.isnan(half):
            ax.plot([x, x], [m - half, m + half],
                    color=COHORT_COLORS[c], lw=0.7, zorder=4)

    # x-axis: SR / AF tick labels
    ticks = [x for x in x_pos.values()]
    labels = [rhy for (_, rhy) in x_pos.keys()]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=6)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_xlim(-0.6, base - 0.8)
    ax.set_ylim(*y_lim)
    # Cohort labels placed in axes-fraction coords, below the x-tick labels.
    for c, cx in centres:
        ax.text(cx, -0.14, COHORT_SHORT[c],
                ha="center", va="top",
                fontsize=6.5, color=COHORT_COLORS[c], fontweight="bold",
                transform=ax.get_xaxis_transform(), clip_on=False)

    # -- significance bracket per cohort (no numeric p, same '*'/'.' grammar
    # as Figure 2C).  Silent when neither threshold is crossed, to avoid
    # clutter.
    span = y_lim[1] - y_lim[0]
    y_brk = y_lim[1] - 0.10 * span
    y_txt = y_brk + 0.02 * span
    for c, cx in centres:
        sr_x = x_pos[(c, "SR")]
        af_x = x_pos[(c, "AF")]
        p = prog_row.get(f"p_perm_{c}", float("nan"))
        marker = perm_marker(p)
        if not marker:
            continue
        ax.plot([sr_x, sr_x, af_x, af_x],
                [y_brk - 0.015 * span, y_brk, y_brk,
                 y_brk - 0.015 * span],
                color=COHORT_COLORS[c], lw=0.6, clip_on=False)
        ax.text(cx, y_txt, marker,
                ha="center", va="bottom", fontsize=9, fontweight="bold",
                color=COHORT_COLORS[c], clip_on=False)

    # y-axis ticks: only for leftmost panel (handled by caller)
    ax.tick_params(axis="y", labelsize=6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)

    # Title: program . consensus_class, with a small coloured band above
    # the panel indicating the consensus class (drawn in axes fraction coords
    # so it never spills into the legend).
    title_col = REWIRING_CLASS_COLORS.get(pick["consensus_class"], "#333333")
    from matplotlib.transforms import blended_transform_factory
    trans = blended_transform_factory(ax.transAxes, ax.transAxes)
    ax.add_patch(plt.Rectangle((0.0, 1.01), 1.0, 0.045,
                               transform=trans,
                               facecolor=title_col, edgecolor="none",
                               clip_on=False, zorder=5))
    ax.set_title(f"{pick['program_label']}  ·  {pick['consensus_class']}",
                 fontsize=7, color="#222222", pad=9)


def main() -> None:
    set_nature_style()
    progs, src = load_data()
    picks = select_representatives(progs)
    if not picks:
        print("no representative programs selected; nothing to plot")
        return

    # Shared y range covers all selected programs with a small margin.
    # Modest headroom at the top for optional significance brackets.
    sel_ids = [p["program_id"] for p in picks]
    sel_src = src[src["program_id"].isin(sel_ids)]
    y_lo = min(-0.4, float(sel_src["aCAI_rel"].min()) - 0.1)
    y_hi = max( 1.9, float(sel_src["aCAI_rel"].max()) + 0.25)

    n_panels = len(picks)
    # Reserve ~1.6 in on the right edge for a dedicated legend column.
    per_panel_w = 1.90   # inches each
    legend_w    = 1.40
    panel_gap   = 0.20   # inches between panels
    fig_w = per_panel_w * n_panels + panel_gap * (n_panels - 1) + legend_w + 0.7
    fig_h = 2.8
    fig = plt.figure(figsize=(fig_w, fig_h))
    left_pad = 0.55        # inches for ylabel
    axes = []
    for k in range(n_panels):
        x_in = left_pad + k * (per_panel_w + panel_gap)
        ax = fig.add_axes([x_in / fig_w, 0.22, per_panel_w / fig_w, 0.60])
        if k > 0:
            ax.sharey(axes[0])
        axes.append(ax)
    panel_right_in = left_pad + n_panels * per_panel_w + (n_panels - 1) * panel_gap
    panel_right    = panel_right_in / fig_w

    for ax, pick in zip(axes, picks):
        prog_row = progs[progs["program_id"] == pick["program_id"]].iloc[0]
        plot_panel(ax, pick, src, (y_lo, y_hi), prog_row)
    for ax in axes[1:]:
        ax.tick_params(axis="y", labelleft=False)

    axes[0].set_ylabel(r"$\mathrm{aCAI}_{rel}$", fontsize=7)

    # Dedicated legend axes on the right.
    leg_ax = fig.add_axes([panel_right + 0.02, 0.22, legend_w / fig_w - 0.03, 0.60])
    leg_ax.axis("off")
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", color=COHORT_COLORS["GSE128188"],
               linestyle="None", markersize=5, label="GSE128188 (RNA-seq)"),
        Line2D([], [], marker="^", color=COHORT_COLORS["GSE79768"],
               linestyle="None", markersize=5, label="GSE79768 (array)"),
        Line2D([], [], marker="D", color=COHORT_COLORS["GSE115574"],
               linestyle="None", markersize=5, label="GSE115574 (array, MR)"),
        Line2D([], [], color="#444444", linestyle=":",  lw=0.8,
               label="Fig 1 healthy (+1)"),
        Line2D([], [], color="#888888", linestyle="--", lw=0.6,
               label="chamber-bias absent"),
    ]
    leg_ax.legend(handles=handles, loc="center left",
                  fontsize=6, frameon=False,
                  handlelength=1.1, handleheight=0.7, labelspacing=0.6,
                  borderaxespad=0)
    # Significance marker key under the main legend.
    leg_ax.text(0.0, 0.18, "perm p", fontsize=6, fontweight="bold",
                ha="left", va="bottom", transform=leg_ax.transAxes,
                color="#333333")
    leg_ax.text(0.0, 0.10, r"$\star$  < 0.05",
                fontsize=6, ha="left", va="center",
                transform=leg_ax.transAxes, color="#333333")
    leg_ax.text(0.0, 0.03, r".  < 0.10",
                fontsize=6, ha="left", va="center",
                transform=leg_ax.transAxes, color="#333333")

    # footer if some classes were skipped
    seen = {p["panel_class"] for p in picks}
    skipped = [c for c in CLASS_ORDER if c not in seen]
    if skipped:
        fig.text(0.5, 0.03,
                 "Classes with no consensus candidate (skipped): " + ", ".join(skipped),
                 fontsize=6, ha="center", va="bottom", color="#777777")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}  ({n_panels} panels: "
          + ", ".join(f"{p['panel_class']}={p['program_label']}" for p in picks) + ")")


if __name__ == "__main__":
    main()
