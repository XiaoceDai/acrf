#!/usr/bin/env python3
"""
Figure 2D — Chamber decomposition (mechanism direction).

This panel visualizes family-level AF rewiring as chamber-component shifts:
    x-axis : Δscore_LA = mean(score_LA | AF) - mean(score_LA | SR)
    y-axis : Δscore_RA = mean(score_RA | AF) - mean(score_RA | SR)

Compared with a 1D effect-size forest, this 2D view explains *where* the
rewiring comes from: LA-dominant, RA-dominant, or bilateral changes.

Inputs:
    results/fig2_primary_acai_long.tsv
    results/fig2_primary_axis_families.tsv
Outputs:
    results/fig2d_family_forest.tsv
    figures/fig2d_family_forest.pdf
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from _pathway_labels import short_pathway_label

DEFAULT_ACAI_TSV = ROOT / "results" / "fig2_primary_acai_long.tsv"
DEFAULT_FAM_TSV  = ROOT / "results" / "fig2_primary_axis_families.tsv"
DEFAULT_OUT_TSV  = ROOT / "results" / "fig2d_family_forest.tsv"
DEFAULT_FIG_OUT  = ROOT / "figures" / "fig2d_family_forest.pdf"

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
GREY_LIGHT = "#E6E6E6"
GREY_MID   = "#9C9C9C"
GREY_DARK  = "#444444"
COHORT_COLORS = {
    "GSE128188": "#A6A6A6",
    "GSE79768":  "#9DBAA9",
    "GSE115574": "#A09BC9",
}
DIR_MARKERS = {"LA-biased": "o", "RA-biased": "s"}

DRIVER_COLORS = {
    "LA-driven": LA_COLOR,
    "RA-driven": RA_COLOR,
    "Bilateral": "#A6A6A6",
}

RATIO_THRESH = 1.25
SQRT2 = np.sqrt(2.0)
BUBBLE_SCALE = 0.55

QUADRANT_BIO_LABELS = {
    "Q1": "Bilateral activation\nLA\u2191 RA\u2191",
    "Q2": "RA-compensatory shift\nLA\u2193 RA\u2191",
    "Q3": "Bilateral attenuation\nLA\u2193 RA\u2193",
    "Q4": "LA-dominant divergence\nLA\u2191 RA\u2193",
}


def sanitize_label(text: str) -> str:
    return (
        text.replace("²", "2")
        .replace("³", "3")
        .replace("⁺", "+")
        .replace("⁻", "-")
        .replace("β", "beta")
        .replace("α", "alpha")
    )


def set_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                             "DejaVu Sans"],
        "font.size":         8,
        "axes.titlesize":    9,
        "axes.labelsize":    8,
        "xtick.labelsize":   7,
        "ytick.labelsize":   7,
        "legend.fontsize":   6.5,
        "axes.linewidth":    0.7,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth":   0.8,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "figure.dpi":        150,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,
    })


def hedges_g(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Return (Hedges' g, variance of g).  AF − SR convention with
    standardisation by the pooled SD; small-sample correction J applied.
    """
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan, np.nan
    mx, my = np.mean(x), np.mean(y)
    sx, sy = np.std(x, ddof=1), np.std(y, ddof=1)
    sp = np.sqrt(((nx - 1) * sx**2 + (ny - 1) * sy**2) / (nx + ny - 2))
    if sp <= 0:
        return np.nan, np.nan
    d = (mx - my) / sp
    df = nx + ny - 2
    J = 1 - (3 / (4 * df - 1))
    g = J * d
    var_g = (nx + ny) / (nx * ny) + g**2 / (2 * (nx + ny))
    return float(g), float(var_g)


def classify_quadrant(delta_la: float, delta_ra: float) -> str:
    if delta_la >= 0 and delta_ra >= 0:
        return "Q1 (+,+)"
    if delta_la < 0 <= delta_ra:
        return "Q2 (-,+)"
    if delta_la < 0 and delta_ra < 0:
        return "Q3 (-,-)"
    return "Q4 (+,-)"


def classify_driver(delta_la: float, delta_ra: float, *, ratio_thresh: float = RATIO_THRESH) -> str:
    ala = abs(delta_la)
    ara = abs(delta_ra)
    if ala >= ara * ratio_thresh:
        return "LA-driven"
    if ara >= ala * ratio_thresh:
        return "RA-driven"
    return "Bilateral"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Draw family-level chamber decomposition map.")
    p.add_argument("--acai-tsv", default=str(DEFAULT_ACAI_TSV),
                   help="Input donor-level aCAI table.")
    p.add_argument("--fam-tsv", default=str(DEFAULT_FAM_TSV),
                   help="Input pathway-family mapping table.")
    p.add_argument("--out-tsv", default=str(DEFAULT_OUT_TSV),
                   help="Output decomposition TSV path.")
    p.add_argument("--out-pdf", default=str(DEFAULT_FIG_OUT),
                   help="Output decomposition PDF path.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    acai_tsv = Path(args.acai_tsv)
    fam_tsv = Path(args.fam_tsv)
    out_tsv = Path(args.out_tsv)
    fig_out = Path(args.out_pdf)
    if not acai_tsv.is_absolute():
        acai_tsv = ROOT / acai_tsv
    if not fam_tsv.is_absolute():
        fam_tsv = ROOT / fam_tsv
    if not out_tsv.is_absolute():
        out_tsv = ROOT / out_tsv
    if not fig_out.is_absolute():
        fig_out = ROOT / fig_out

    set_publication_style()

    long_df = pd.read_csv(acai_tsv, sep="\t")
    fam = pd.read_csv(fam_tsv, sep="\t")[["pathway", "family_id"]].rename(
        columns={"pathway": "program_id"})
    long_df = long_df.merge(fam, on="program_id", how="left")
    long_df = long_df.dropna(subset=["family_id"])
    long_df["family_id"] = long_df["family_id"].astype(int)

    # Per donor, per family, per cohort: family-level aCAI + chamber scores
    fam_donor = (long_df.groupby(["cohort", "family_id", "donor_id",
                                     "rhythm", "direction"], as_index=False)
                          [["aCAI", "score_LA", "score_RA"]].mean()
                          .rename(columns={"aCAI": "family_aCAI"}))

    cohorts = ["GSE128188", "GSE79768", "GSE115574"]

    # Representative family label
    fam_rep = (long_df.assign(absACAI=long_df["aCAI"].abs())
                      .sort_values("absACAI", ascending=False)
                      .groupby("family_id")
                      .agg(rep=("program_id", "first"))
                      .reset_index())
    fam_rep["short"] = fam_rep["rep"].apply(
        lambda s: sanitize_label(short_pathway_label(s, 22)))

    # Compute decomposition + g per (family, cohort)
    rows = []
    for fid in sorted(fam_donor["family_id"].unique()):
        fam_sub = fam_donor[fam_donor["family_id"] == fid]
        direction = fam_sub["direction"].iloc[0]
        for cohort in cohorts:
            cs = fam_sub[fam_sub["cohort"] == cohort]
            af = cs[cs["rhythm"] == "AF"]
            sr = cs[cs["rhythm"] == "SR"]
            af_vals = af["family_aCAI"].values
            sr_vals = sr["family_aCAI"].values
            if len(af_vals) == 0 or len(sr_vals) == 0:
                continue

            delta_la = float(af["score_LA"].mean() - sr["score_LA"].mean())
            delta_ra = float(af["score_RA"].mean() - sr["score_RA"].mean())
            delta_acai = float(delta_la - delta_ra)
            quad = classify_quadrant(delta_la, delta_ra)
            driver = classify_driver(delta_la, delta_ra)

            g, var_g = hedges_g(af_vals, sr_vals)
            rows.append({
                "family_id": fid, "direction": direction,
                "cohort": cohort,
                "n_AF": len(af_vals), "n_SR": len(sr_vals),
                "delta_LA": round(delta_la, 4),
                "delta_RA": round(delta_ra, 4),
                "delta_aCAI": round(delta_acai, 4),
                "abs_delta_aCAI": round(abs(delta_acai), 4),
                "dist_to_y_eq_x": round(abs(delta_ra - delta_la) / SQRT2, 4),
                "dist_to_y_eq_negx": round(abs(delta_ra + delta_la) / SQRT2, 4),
                "quadrant": quad,
                "driver": driver,
                "g": round(g, 4) if not np.isnan(g) else np.nan,
                "var_g": round(var_g, 4) if not np.isnan(var_g) else np.nan,
                "ci_lo": round(g - 1.96 * np.sqrt(var_g), 4) if not np.isnan(g) else np.nan,
                "ci_hi": round(g + 1.96 * np.sqrt(var_g), 4) if not np.isnan(g) else np.nan,
            })

        # Pooled row retained for downstream compatibility (e.g. Fig 5B).
        per_coh = [r for r in rows[-3:] if not np.isnan(r["g"])]
        af_all = fam_sub[fam_sub["rhythm"] == "AF"]
        sr_all = fam_sub[fam_sub["rhythm"] == "SR"]
        if len(per_coh) >= 2 and len(af_all) > 0 and len(sr_all) > 0:
            ws = [1 / r["var_g"] for r in per_coh]
            sum_w = sum(ws)
            pooled_g = sum(r["g"] * w for r, w in zip(per_coh, ws)) / sum_w
            pooled_var = 1 / sum_w
            p_delta_la = float(af_all["score_LA"].mean() - sr_all["score_LA"].mean())
            p_delta_ra = float(af_all["score_RA"].mean() - sr_all["score_RA"].mean())
            p_delta_acai = float(p_delta_la - p_delta_ra)
            rows.append({
                "family_id": fid, "direction": direction,
                "cohort": "POOLED",
                "n_AF": int(len(af_all)),
                "n_SR": int(len(sr_all)),
                "delta_LA": round(p_delta_la, 4),
                "delta_RA": round(p_delta_ra, 4),
                "delta_aCAI": round(p_delta_acai, 4),
                "abs_delta_aCAI": round(abs(p_delta_acai), 4),
                "dist_to_y_eq_x": round(abs(p_delta_ra - p_delta_la) / SQRT2, 4),
                "dist_to_y_eq_negx": round(abs(p_delta_ra + p_delta_la) / SQRT2, 4),
                "quadrant": classify_quadrant(p_delta_la, p_delta_ra),
                "driver": classify_driver(p_delta_la, p_delta_ra),
                "g": round(pooled_g, 4),
                "var_g": round(pooled_var, 4),
                "ci_lo": round(pooled_g - 1.96 * np.sqrt(pooled_var), 4),
                "ci_hi": round(pooled_g + 1.96 * np.sqrt(pooled_var), 4),
            })

    forest = pd.DataFrame(rows).merge(
        fam_rep[["family_id", "short"]], on="family_id", how="left"
    )
    cohort_rank = {c: i for i, c in enumerate(cohorts + ["POOLED"])}
    forest["_rank"] = forest["cohort"].map(cohort_rank).fillna(99)
    forest = forest.sort_values(["direction", "family_id", "_rank"]).drop(columns=["_rank"])
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    forest.to_csv(out_tsv, sep="\t", index=False)
    print(f"wrote {out_tsv.relative_to(ROOT)}  ({len(forest)} rows)")

    # ---- Figure: 3 cohort decomposition panels + summary ------------------
    plot_df = forest[forest["cohort"].isin(cohorts)].copy()
    max_abs = float(np.nanmax(np.abs(
        np.r_[plot_df["delta_LA"].values, plot_df["delta_RA"].values]
    )))
    lim = max(0.08, max_abs * 1.22)
    label_dx = lim * 0.014
    label_dy = lim * 0.010

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(188 * MM, 120 * MM))
    gs = fig.add_gridspec(
        1, 4, width_ratios=[1.0, 1.0, 1.0, 0.90],
        left=0.07, right=0.98, top=0.90, bottom=0.13, wspace=0.14
    )
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax_leg = fig.add_subplot(gs[0, 3])

    # Quadrant backgrounds (very light, consistent across cohorts).
    quad_bg = {
        "Q1": "#FAFAFA",   # +,+
        "Q2": "#FFF4F0",   # -,+
        "Q3": "#F7FAFD",   # -,-
        "Q4": "#F2F8F4",   # +,-
    }
    quad_labels = {
        "Q1": QUADRANT_BIO_LABELS["Q1"],
        "Q2": QUADRANT_BIO_LABELS["Q2"],
        "Q3": QUADRANT_BIO_LABELS["Q3"],
        "Q4": QUADRANT_BIO_LABELS["Q4"],
    }
    quad_pos = {"Q1": (0.03, 0.97), "Q2": (0.03, 0.03), "Q3": (0.60, 0.03), "Q4": (0.60, 0.97)}

    for ax, cohort in zip(axes, cohorts):
        sub = plot_df[plot_df["cohort"] == cohort].copy()

        ax.add_patch(mpatches.Rectangle((-lim, 0), lim, lim, facecolor=quad_bg["Q2"], edgecolor="none", zorder=0))
        ax.add_patch(mpatches.Rectangle((0, 0), lim, lim, facecolor=quad_bg["Q1"], edgecolor="none", zorder=0))
        ax.add_patch(mpatches.Rectangle((-lim, -lim), lim, lim, facecolor=quad_bg["Q3"], edgecolor="none", zorder=0))
        ax.add_patch(mpatches.Rectangle((0, -lim), lim, lim, facecolor=quad_bg["Q4"], edgecolor="none", zorder=0))

        ax.axhline(0, color=GREY_MID, lw=0.7, zorder=1)
        ax.axvline(0, color=GREY_MID, lw=0.7, zorder=1)
        ax.plot([-lim, lim], [-lim, lim], color=GREY_LIGHT, lw=0.6, ls="--", zorder=1)
        ax.plot([-lim, lim], [lim, -lim], color=GREY_LIGHT, lw=0.6, ls="--", zorder=1)

        # Quadrant corner labels.
        for q, (xp, yp) in quad_pos.items():
            ax.text(xp, yp, quad_labels[q], transform=ax.transAxes,
                    ha="left", va="top" if yp > 0.5 else "bottom",
                    fontsize=5.0, color=GREY_DARK, linespacing=1.0)

        # Family points.
        for _, r in sub.iterrows():
            fid = int(r["family_id"])
            x = float(r["delta_LA"])
            y = float(r["delta_RA"])
            s = (35 + 260 * (abs(float(r["delta_aCAI"])) / (max_abs + 1e-9))) * BUBBLE_SCALE
            face = DRIVER_COLORS.get(r["driver"], "#A6A6A6")
            m = DIR_MARKERS.get(r["direction"], "o")
            ax.scatter(x, y, s=s, marker=m, facecolor=face, edgecolor="black",
                       linewidth=0.35, alpha=0.92, zorder=3)
            ax.text(x + label_dx, y + label_dy, f"F{fid}", fontsize=5.3, color="#222222",
                    ha="left", va="center", zorder=4)

        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{cohort}  (n families = {len(sub)})", fontsize=7.2, fontweight="bold", pad=4)
        ax.tick_params(axis="both", labelsize=6.0, length=2.5)

    axes[0].set_ylabel(r"$\Delta\mathrm{score}_{RA}$  (AF - SR)", fontsize=7.2, color=RA_COLOR)
    for ax in axes:
        ax.set_xlabel(r"$\Delta\mathrm{score}_{LA}$  (AF - SR)", fontsize=7.2, color=LA_COLOR)

    # Right summary and legends.
    ax_leg.axis("off")
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(1.22, 0)
    ax_leg.text(0.0, 0.02, "Mechanism legend", fontsize=7.2, fontweight="bold", ha="left", va="top")

    y = 0.09
    for d in ["LA-driven", "RA-driven", "Bilateral"]:
        ax_leg.scatter([0.04], [y - 0.005], s=40, marker="o",
                       facecolor=DRIVER_COLORS[d], edgecolor="#555555", linewidth=0.45)
        ax_leg.text(0.11, y - 0.005, d, fontsize=6.3, ha="left", va="center", color=GREY_DARK)
        y += 0.06

    ax_leg.plot([0.0, 0.88], [y - 0.01, y - 0.01], color=GREY_LIGHT, lw=0.6)
    y += 0.01
    ax_leg.text(0.0, y, r"Bubble size = $|\Delta aCAI|$", fontsize=6.8, fontweight="bold", ha="left", va="top")
    y += 0.05
    for frac in [0.25, 0.60, 1.00]:
        size = (35 + 260 * frac) * BUBBLE_SCALE
        val = frac * max_abs
        ax_leg.scatter([0.04], [y], s=size, marker="o", facecolor="#F2F2F2", edgecolor="#7F7F7F", linewidth=0.5)
        ax_leg.text(0.11, y, f"{val:.2f}", fontsize=6.0, ha="left", va="center", color=GREY_DARK)
        y += 0.05

    ax_leg.plot([0.0, 0.88], [y - 0.01, y - 0.01], color=GREY_LIGHT, lw=0.6)
    ax_leg.text(0.0, y, "Baseline direction marker", fontsize=6.8, fontweight="bold", ha="left", va="top")
    y += 0.05
    ax_leg.scatter([0.04], [y], s=32, marker="o", facecolor="white", edgecolor="black", linewidth=0.5)
    ax_leg.text(0.11, y, "LA-biased family", fontsize=6.2, ha="left", va="center")
    y += 0.05
    ax_leg.scatter([0.04], [y], s=32, marker="s", facecolor="white", edgecolor="black", linewidth=0.5)
    ax_leg.text(0.11, y, "RA-biased family", fontsize=6.2, ha="left", va="center")

    y += 0.05
    ax_leg.plot([0.0, 0.88], [y - 0.01, y - 0.01], color=GREY_LIGHT, lw=0.6)
    ax_leg.text(0.0, y, "Diagonal guides", fontsize=6.8, fontweight="bold", ha="left", va="top")
    y += 0.04
    ax_leg.text(0.0, y, "y = x  : symmetric co-shift", fontsize=5.8, ha="left", va="top", color=GREY_DARK)
    y += 0.035
    ax_leg.text(0.0, y, "y = -x : antagonistic shift", fontsize=5.8, ha="left", va="top", color=GREY_DARK)

    y += 0.05
    ax_leg.text(0.0, y, "Driver counts by cohort", fontsize=6.8, fontweight="bold", ha="left", va="top")
    y += 0.04
    ctab = pd.crosstab(plot_df["cohort"], plot_df["driver"]).reindex(cohorts).fillna(0).astype(int)
    for cohort in cohorts:
        row = ctab.loc[cohort]
        txt = (
            f"{cohort}: "
            f"LA {row.get('LA-driven', 0)} | "
            f"RA {row.get('RA-driven', 0)} | "
            f"Bi {row.get('Bilateral', 0)}"
        )
        ax_leg.text(0.0, y, txt, fontsize=5.8, ha="left", va="top", color=GREY_DARK)
        y += 0.055

    fig.suptitle(
        "Fig 2D  Chamber decomposition of AF rewiring (family-level mechanism direction)",
        x=0.07, y=0.965, ha="left", fontsize=8.6, fontweight="bold"
    )
    fig.text(
        0.07, 0.03,
        "Each point is one family in one cohort. "
        "x = \u0394score_LA, y = \u0394score_RA (AF - SR). "
        "Point color denotes dominant chamber source "
        f"(|LA| >= {RATIO_THRESH:.2f}|RA|: LA-driven; |RA| >= {RATIO_THRESH:.2f}|LA|: RA-driven; else bilateral). "
        "Dashed diagonals are y = x (co-shift) and y = -x (antagonistic shift).",
        fontsize=5.8, color=GREY_DARK, ha="left", va="bottom"
    )

    fig_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_out, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {fig_out.relative_to(ROOT)}")

    # ---- console summary ----------------------------------------------
    print("\n=== driver counts by cohort ===")
    print(pd.crosstab(plot_df["cohort"], plot_df["driver"]).to_string())


if __name__ == "__main__":
    main()
