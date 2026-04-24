#!/usr/bin/env python3
"""
Figure 2E (v3) — Cross-cohort consolidation.

Left panel: family-level pooled effect forest (from Fig 2D outputs).
Middle panel: Fig 2C rewiring-state mini-strip across 3 cohorts + state concordance.
Right panel: replication support (lollipop + external support point + heterogeneity):
    - lollipop x: same_sign_n (0..3): number of cohorts whose effect sign
      matches pooled sign.
    - external square: cohort_sig_n (0..3): number of cohorts with CI excluding zero.
    - external text: I2 (%): between-cohort heterogeneity.

Outputs
-------
results/fig2e_replication_lollipop_v3.tsv
figures/fig2e_replication_lollipop_v3.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from scipy.stats import chi2

from plot_style import LA_COLOR, RA_COLOR, set_nature_style

ROOT = Path(__file__).resolve().parent.parent

FOREST_TSV = ROOT / "results" / "fig2d_family_forest_v3.tsv"
STATE_TSV = ROOT / "results" / "fig2c_chamber_decomposition_v3.tsv"
OUT_TSV = ROOT / "results" / "fig2e_replication_lollipop_v3.tsv"
OUT_PDF = ROOT / "figures" / "fig2e_replication_lollipop_v3.pdf"

COHORTS = ["GSE128188", "GSE79768", "GSE115574"]
COHORT_COLORS = {"GSE128188": "#A6A6A6", "GSE79768": "#9DBAA9", "GSE115574": "#A09BC9"}
COHORT_MARKERS = {"GSE128188": "o", "GSE79768": "^", "GSE115574": "s"}

REP_CLASS_COLORS = {
    "high": "#1B5E20",
    "moderate": "#E08214",
    "low": "#8E8E8E",
}
SIG_COUNT_COLORS = {0: "#F0F0F0", 1: "#D9D9D9", 2: "#9ECAE1", 3: "#3182BD"}
STATE_COLORS = {
    "Preserved": "#A9C9B5",
    "Attenuated": "#E9C27B",
    "Reversed": "#C6A2D7",
    "Amplified": "#E59595",
}
CONC_COLORS = {"stable": "#4E79A7", "heterogeneous": "#B07AA1"}


def sign_nonzero(x: float, eps: float = 1e-8) -> int:
    if x > eps:
        return 1
    if x < -eps:
        return -1
    return 0


def classify_replication(pooled_sig: bool, same_sign_n: int) -> str:
    if pooled_sig and same_sign_n >= 2:
        return "high"
    if same_sign_n >= 2:
        return "moderate"
    return "low"


def heterogeneity_metrics(g: np.ndarray, v: np.ndarray) -> tuple[float, float, float]:
    mask = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g = g[mask]
    v = v[mask]
    k = len(g)
    if k < 2:
        return np.nan, np.nan, np.nan
    w = 1.0 / v
    g_bar = float(np.sum(w * g) / np.sum(w))
    q = float(np.sum(w * (g - g_bar) ** 2))
    df = k - 1
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    p_het = float(1.0 - chi2.cdf(q, df))
    return q, i2, p_het


def state_concordance(states: list[str]) -> tuple[str, float, int]:
    s = [x for x in states if isinstance(x, str) and x]
    if not s:
        return "NA", np.nan, 0
    vc = pd.Series(s).value_counts()
    dominant = str(vc.index[0])
    n = int(vc.iloc[0])
    return dominant, n / len(states), n


def main() -> None:
    set_nature_style()

    forest = pd.read_csv(FOREST_TSV, sep="\t")
    state = pd.read_csv(STATE_TSV, sep="\t")

    req = {"family_id", "short", "direction", "cohort", "g", "ci_lo", "ci_hi", "var_g"}
    missing = sorted(req - set(forest.columns))
    if missing:
        raise ValueError(f"Missing required columns in {FOREST_TSV.name}: {missing}")
    state_req = {"family_id", "cohort", "state"}
    state_missing = sorted(state_req - set(state.columns))
    if state_missing:
        raise ValueError(f"Missing required columns in {STATE_TSV.name}: {state_missing}")

    pooled = forest[forest["cohort"] == "POOLED"].copy()
    if pooled.empty:
        raise ValueError("No POOLED rows found in fig2d_family_forest_v3.tsv")

    rows = []
    for _, rp in pooled.iterrows():
        fid = int(rp["family_id"])
        direction = str(rp["direction"])
        sub = forest[(forest["family_id"] == fid) & (forest["cohort"].isin(COHORTS))].copy()
        if sub.empty:
            continue

        pooled_g = float(rp["g"])
        pooled_lo = float(rp["ci_lo"])
        pooled_hi = float(rp["ci_hi"])
        pooled_sig = bool((pooled_lo > 0) or (pooled_hi < 0))
        pooled_sign = sign_nonzero(pooled_g)

        same_sign_n = 0
        cohort_sig_n = 0
        cohort_metrics: dict[str, dict[str, float | bool]] = {}
        for c in COHORTS:
            rc = sub[sub["cohort"] == c]
            if rc.empty:
                cohort_metrics[c] = {"g": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "var_g": np.nan, "same_sign": False, "sig": False}
                continue
            rr = rc.iloc[0]
            g = float(rr["g"])
            lo = float(rr["ci_lo"])
            hi = float(rr["ci_hi"])
            var = float(rr["var_g"]) if np.isfinite(rr["var_g"]) else np.nan
            sig = bool((lo > 0) or (hi < 0))
            same = (pooled_sign != 0) and (sign_nonzero(g) == pooled_sign)
            same_sign_n += int(same)
            cohort_sig_n += int(sig)
            cohort_metrics[c] = {"g": g, "ci_lo": lo, "ci_hi": hi, "var_g": var, "same_sign": same, "sig": sig}

        rep_class = classify_replication(pooled_sig, same_sign_n)
        g_arr = np.array([cohort_metrics[c]["g"] for c in COHORTS], dtype=float)
        v_arr = np.array([cohort_metrics[c]["var_g"] for c in COHORTS], dtype=float)
        q, i2, p_het = heterogeneity_metrics(g_arr, v_arr)

        states = []
        state_by_cohort: dict[str, str] = {}
        for c in COHORTS:
            rs = state[(state["family_id"] == fid) & (state["cohort"] == c)]
            s = str(rs["state"].iloc[0]) if not rs.empty else "NA"
            state_by_cohort[c] = s
            states.append(s)
        dominant_state, st_ratio, st_n = state_concordance(states)
        st_class = "stable" if st_n >= 2 else "heterogeneous"

        rows.append({
            "family_id": fid,
            "short": str(rp["short"]),
            "direction": direction,
            "pooled_g": pooled_g,
            "pooled_ci_lo": pooled_lo,
            "pooled_ci_hi": pooled_hi,
            "pooled_sig": pooled_sig,
            "same_sign_n": int(same_sign_n),
            "cohort_sig_n": int(cohort_sig_n),
            "replication_class": rep_class,
            "Q": q,
            "I2": i2,
            "p_het": p_het,
            "state_128188": state_by_cohort["GSE128188"],
            "state_79768": state_by_cohort["GSE79768"],
            "state_115574": state_by_cohort["GSE115574"],
            "dominant_state": dominant_state,
            "state_concordance_ratio": st_ratio,
            "state_concordance_n": st_n,
            "state_concordance_class": st_class,
            "g_128188": cohort_metrics["GSE128188"]["g"],
            "ci_lo_128188": cohort_metrics["GSE128188"]["ci_lo"],
            "ci_hi_128188": cohort_metrics["GSE128188"]["ci_hi"],
            "g_79768": cohort_metrics["GSE79768"]["g"],
            "ci_lo_79768": cohort_metrics["GSE79768"]["ci_lo"],
            "ci_hi_79768": cohort_metrics["GSE79768"]["ci_hi"],
            "g_115574": cohort_metrics["GSE115574"]["g"],
            "ci_lo_115574": cohort_metrics["GSE115574"]["ci_lo"],
            "ci_hi_115574": cohort_metrics["GSE115574"]["ci_hi"],
        })

    out = pd.DataFrame(rows)
    out["dir_rank"] = out["direction"].map({"LA-biased": 0, "RA-biased": 1}).fillna(2)
    la = out[out["direction"] == "LA-biased"].sort_values("pooled_g", ascending=True)
    ra = out[out["direction"] == "RA-biased"].sort_values("pooled_g", ascending=False)
    out = pd.concat([la, ra], axis=0).reset_index(drop=True)
    out = out.drop(columns=["dir_rank"])
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(out)} rows)")

    # --- Figure ---
    n = len(out)
    fig_h = 0.33 * n + 1.55
    fig = plt.figure(figsize=(8.65, fig_h))
    gs = fig.add_gridspec(
        1, 4, width_ratios=[6.2, 2.35, 3.25, 1.7],
        left=0.35, right=0.985, top=0.94, bottom=0.09, wspace=0.06
    )
    ax_forest = fig.add_subplot(gs[0, 0])
    ax_state = fig.add_subplot(gs[0, 1])
    ax_sup = fig.add_subplot(gs[0, 2])
    ax_leg = fig.add_subplot(gs[0, 3])
    ax_leg.axis("off")

    g_vals = np.r_[
        out["pooled_ci_lo"].to_numpy(float),
        out["pooled_ci_hi"].to_numpy(float),
        out["g_128188"].to_numpy(float),
        out["g_79768"].to_numpy(float),
        out["g_115574"].to_numpy(float),
    ]
    max_abs = float(np.nanmax(np.abs(g_vals))) if np.isfinite(np.nanmax(np.abs(g_vals))) else 0.8
    x_lim = max(0.6, max_abs * 1.18)

    ax_forest.set_xlim(-x_lim, x_lim)
    ax_forest.set_ylim(n - 0.5, -0.5)
    ax_state.set_ylim(n - 0.5, -0.5)
    ax_sup.set_ylim(n - 0.5, -0.5)
    ax_forest.axvline(0, color="#9C9C9C", lw=0.7, ls="--", zorder=1)
    for i in range(n):
        y = i + 0.5
        if i < n - 1:
            ax_forest.axhline(y, color="#ECECEC", lw=0.5, zorder=0)
            ax_state.axhline(y, color="#F0F0F0", lw=0.5, zorder=0)
            ax_sup.axhline(y, color="#F0F0F0", lw=0.5, zorder=0)

    offsets = {"GSE128188": +0.20, "GSE79768": 0.0, "GSE115574": -0.20}
    for i, r in out.iterrows():
        # cohort estimates
        for c in COHORTS:
            g = float(r[f"g_{c.split('GSE')[1]}"]) if c != "GSE115574" else float(r["g_115574"])
            lo = float(r[f"ci_lo_{c.split('GSE')[1]}"]) if c != "GSE115574" else float(r["ci_lo_115574"])
            hi = float(r[f"ci_hi_{c.split('GSE')[1]}"]) if c != "GSE115574" else float(r["ci_hi_115574"])
            if np.isnan(g):
                continue
            y = i + offsets[c]
            ax_forest.plot([lo, hi], [y, y], color=COHORT_COLORS[c], lw=0.9, alpha=0.85, zorder=2)
            ax_forest.scatter(g, y, s=20, marker=COHORT_MARKERS[c], facecolor="white",
                              edgecolor=COHORT_COLORS[c], linewidth=0.9, zorder=3)

        # pooled diamond
        pg = float(r["pooled_g"])
        plo = float(r["pooled_ci_lo"])
        phi = float(r["pooled_ci_hi"])
        col = LA_COLOR if r["direction"] == "LA-biased" else RA_COLOR
        ax_forest.plot([plo, phi], [i, i], color=col, lw=1.5, zorder=4)
        d = Polygon([(pg - 0.055, i), (pg, i + 0.13), (pg + 0.055, i), (pg, i - 0.13)],
                    closed=True, facecolor=col, edgecolor="#333333", linewidth=0.6, zorder=5)
        ax_forest.add_patch(d)
        ax_forest.text(
            x_lim * 0.965,
            i,
            f"{pg:+.2f}",
            ha="right",
            va="center",
            fontsize=5.4,
            color="#505050",
        )

    # y labels
    labels = [f"F{int(r.family_id)}  {r.short}" for r in out.itertuples()]
    ax_forest.set_yticks(np.arange(n))
    ax_forest.set_yticklabels(labels, fontsize=6.4)
    for tick, d in zip(ax_forest.get_yticklabels(), out["direction"]):
        tick.set_color(LA_COLOR if d == "LA-biased" else RA_COLOR)
    ax_forest.tick_params(axis="y", length=0, pad=2)
    ax_forest.tick_params(axis="x", labelsize=6.5)
    ax_forest.set_xlabel("Pooled Hedges' g (AF - SR)", fontsize=7.0)
    ax_forest.set_title("Pooled effect forest", fontsize=7.4, fontweight="bold", pad=4)
    for s in ("top", "right"):
        ax_forest.spines[s].set_visible(False)
    ax_forest.spines["left"].set_visible(False)
    ax_forest.spines["bottom"].set_linewidth(0.6)

    # Fig 2C state mini-strip + concordance
    ax_state.set_xlim(-0.05, 5.2)
    for i, r in out.iterrows():
        svals = [r["state_128188"], r["state_79768"], r["state_115574"]]
        for j, sv in enumerate(svals):
            c = STATE_COLORS.get(str(sv), "#EEEEEE")
            ax_state.add_patch(
                plt.Rectangle((j + 0.05, i - 0.34), 0.86, 0.68,
                              facecolor=c, edgecolor="white", linewidth=0.6, zorder=2)
            )
        conc = int(r["state_concordance_n"])
        conc_col = CONC_COLORS.get(str(r["state_concordance_class"]), "#8E8E8E")
        dom = str(r["dominant_state"])
        dom_short = {"Preserved": "Pre", "Attenuated": "Att", "Reversed": "Rev", "Amplified": "Amp"}.get(dom, dom[:3])
        ax_state.text(3.18, i, f"{conc}/3", ha="center", va="center",
                      fontsize=5.7, color=conc_col, fontweight="bold", zorder=3)
        ax_state.text(4.28, i, f"{dom_short}", ha="center", va="center",
                      fontsize=5.6, color="#333333", fontweight="bold", zorder=3)
    ax_state.set_yticks(np.arange(n))
    ax_state.set_yticklabels([])
    ax_state.tick_params(axis="y", left=False, labelleft=False)
    ax_state.set_xticks([0.48, 1.48, 2.48, 3.18, 4.28])
    ax_state.set_xticklabels(["128", "797", "115", "conc", "state"], fontsize=5.8)
    ax_state.set_title("Fig2C state", fontsize=7.0, fontweight="bold", pad=4)
    for s in ("top", "right", "left"):
        ax_state.spines[s].set_visible(False)
    ax_state.spines["bottom"].set_linewidth(0.6)

    # Replication support lollipop + external support point + I2
    ax_sup.set_xlim(-0.2, 6.45)
    ax_sup.axvline(3.55, color="#DDDDDD", lw=0.8, zorder=1)
    ax_sup.axvline(4.75, color="#E6E6E6", lw=0.8, zorder=1)
    ax_sup.axvline(5.55, color="#E6E6E6", lw=0.8, zorder=1)
    for i, r in out.iterrows():
        same_n = int(r["same_sign_n"])
        sig_n = int(r["cohort_sig_n"])
        cls = str(r["replication_class"])
        col = REP_CLASS_COLORS[cls]

        ax_sup.plot([0, same_n], [i, i], color="#BDBDBD", lw=1.2, zorder=2)
        ax_sup.scatter([same_n], [i], s=30, marker="o", facecolor=col,
                       edgecolor="#333333", linewidth=0.4, zorder=3)

        ax_sup.scatter([4.35], [i], s=48, marker="s",
                       facecolor=SIG_COUNT_COLORS[sig_n], edgecolor="#666666", linewidth=0.4, zorder=3)
        ax_sup.text(4.35, i, f"{sig_n}", ha="center", va="center", fontsize=5.4, color="#222222", zorder=4)

        i2 = float(r["I2"]) if np.isfinite(r["I2"]) else np.nan
        if np.isfinite(i2):
            ax_sup.text(5.10, i, f"{100.0 * i2:.0f}", ha="center", va="center",
                        fontsize=5.2, color="#444444", zorder=4)
        else:
            ax_sup.text(5.10, i, "NA", ha="center", va="center",
                        fontsize=5.0, color="#999999", zorder=4)
        cls_tag = {"high": "H", "moderate": "M", "low": "L"}.get(cls, "L")
        ax_sup.text(5.95, i, cls_tag, ha="center", va="center",
                    fontsize=5.6, color=col, fontweight="bold", zorder=4)

    ax_sup.set_yticks(np.arange(n))
    ax_sup.set_yticklabels([])
    ax_sup.tick_params(axis="y", left=False, labelleft=False)
    ax_sup.set_xticks([0, 1, 2, 3, 4.35, 5.10, 5.95])
    ax_sup.set_xticklabels(["0", "1", "2", "3", "n_sig", "I²", "cls"], fontsize=6.1)
    ax_sup.set_xlabel("same-sign cohort count", fontsize=7.0)
    ax_sup.set_title("Replication support", fontsize=7.4, fontweight="bold", pad=4)
    for s in ("top", "right", "left"):
        ax_sup.spines[s].set_visible(False)
    ax_sup.spines["bottom"].set_linewidth(0.6)

    # legends
    cohort_handles = [
        Line2D([], [], marker="o", linestyle="None", markersize=4.8,
               markerfacecolor="white", markeredgecolor=COHORT_COLORS["GSE128188"], label="GSE128188"),
        Line2D([], [], marker="^", linestyle="None", markersize=4.8,
               markerfacecolor="white", markeredgecolor=COHORT_COLORS["GSE79768"], label="GSE79768"),
        Line2D([], [], marker="s", linestyle="None", markersize=4.8,
               markerfacecolor="white", markeredgecolor=COHORT_COLORS["GSE115574"], label="GSE115574"),
    ]
    fig.legend(
        handles=cohort_handles, loc="upper center", bbox_to_anchor=(0.66, 0.982),
        ncol=3, frameon=False, fontsize=5.8, columnspacing=1.0, handletextpad=0.3
    )

    handles = [
        Line2D([], [], marker="D", linestyle="None", markersize=5,
               markerfacecolor=LA_COLOR, markeredgecolor="#333333", label="pooled LA-biased"),
        Line2D([], [], marker="D", linestyle="None", markersize=5,
               markerfacecolor=RA_COLOR, markeredgecolor="#333333", label="pooled RA-biased"),
    ]
    ax_leg.legend(handles=handles, loc="upper left", frameon=False, fontsize=5.7, handlelength=1.0, labelspacing=0.35)
    ax_leg.text(0.00, 0.73, "Support class", fontsize=6.2, fontweight="bold", color="#333333", transform=ax_leg.transAxes)
    ax_leg.text(0.00, 0.68, "high: pooled sig + >=2 same-sign", fontsize=5.5, color=REP_CLASS_COLORS["high"], transform=ax_leg.transAxes)
    ax_leg.text(0.00, 0.63, "moderate: >=2 same-sign,\n          pooled CI may cross 0", fontsize=5.3, color=REP_CLASS_COLORS["moderate"], transform=ax_leg.transAxes)
    ax_leg.text(0.00, 0.56, "low: otherwise", fontsize=5.5, color=REP_CLASS_COLORS["low"], transform=ax_leg.transAxes)
    ax_leg.text(0.00, 0.48, "n_sig: # cohorts with CI excluding 0", fontsize=5.4, color="#555555", transform=ax_leg.transAxes)
    ax_leg.text(0.00, 0.43, "I²: between-cohort heterogeneity (%)", fontsize=5.4, color="#555555", transform=ax_leg.transAxes)

    ax_leg.text(0.00, 0.34, "State colors (from Fig 2C)", fontsize=6.0, fontweight="bold", color="#333333", transform=ax_leg.transAxes)
    y0 = 0.29
    for st in ["Amplified", "Preserved", "Attenuated", "Reversed"]:
        ax_leg.add_patch(
            plt.Rectangle((0.00, y0 - 0.015), 0.07, 0.026,
                          transform=ax_leg.transAxes, facecolor=STATE_COLORS[st], edgecolor="none")
        )
        ax_leg.text(0.10, y0 - 0.002, st, fontsize=5.4, color="#444444", transform=ax_leg.transAxes, va="center")
        y0 -= 0.048

    fig.suptitle(
        "Fig 2E (v3)  Cross-cohort consolidation: pooled effect + replication support",
        fontsize=8.4, y=0.985
    )
    fig.text(
        0.35, 0.965,
        "State strip shows Fig 2C rewiring states across cohorts (128188/79768/115574); "
        "concordance reports dominant-state frequency (e.g., 2/3). "
        "Forest right text shows pooled g; support panel adds class code (H/M/L).",
        fontsize=5.8, color="#555555",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
