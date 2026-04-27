#!/usr/bin/env python3
"""
Figure S8 — Specificity / Falsification (4-panel supplement).

S8a — Component A: size-matched random gene-set permutations.
S8b — Component B: axis-membership permutation density of |g|.
S8c — Component C: real frozen support programs vs pre-registered
      negative controls on (pooled g, I², LAA duration ρ).
S8d — Leave-one-program-out (LOPO) axis robustness:
      remove each support program in turn and re-pool.
      Inputs from `scripts/h9_lopo_axis.R`.

Components A–C precomputed by `scripts/h5_falsification.R`.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle

from plot_style import COL_DOUBLE, set_nature_style

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "figures" / "figS8_falsification.pdf"

# Frozen support programs in canonical order
PROGRAM_ORDER = [
    ("REACTOME__REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM",
     "ECM proteoglycans", "STRUCTURAL"),
    ("REACTOME__REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM",
     "HS-GAG Metabolism", "STRUCTURAL"),
    ("REACTOME__REACTOME_HS_GAG_BIOSYNTHESIS",
     "HS-GAG Biosynthesis", "STRUCTURAL"),
    ("HALLMARK__HALLMARK_OXIDATIVE_PHOSPHORYLATION",
     "OXPHOS", "METABOLIC"),
    ("REACTOME__REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
     "Aerobic Resp. & ETC", "METABOLIC"),
    ("REACTOME__REACTOME_TRANSLATION",
     "Cytoplasmic Translation", "METABOLIC"),
]

AXIS_COL = {"STRUCTURAL": "#B35806", "METABOLIC": "#542788"}
CTRL_AS_COL = "#5A8DAD"      # axis-specificity controls
CTRL_UN_COL = "#9E9E9E"      # truly unrelated controls
REAL_COL    = "#1A1A1A"      # real-value markers


def _round(ax, x, y, w, h, fc, ec="none", lw=0.4, z=1, rounding=0.012):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    ))


# ---------------------------------------------------------------------
# Panel S8a — size-matched random gene-set permutations (component A)
# ---------------------------------------------------------------------
def panel_S8a(ax):
    summary = pd.read_csv(RES / "figure_s8a_size_matched.tsv", sep="\t")
    summary["program_label"] = summary["program_label"].astype(str)
    # Order programs per PROGRAM_ORDER
    label_to_idx = {lab: i for i, (_, lab, _) in enumerate(PROGRAM_ORDER)}
    summary["sort"] = summary["program_label"].map(label_to_idx)
    summary = summary.sort_values("sort").reset_index(drop=True)

    n = len(summary)
    rng = np.random.default_rng(7)

    # Read each per-program perm distribution
    perm_files = {
        row["program_id"]: RES / f"figure_s8a_perm_dist_{row['program_id']}.tsv.gz"
        for _, row in summary.iterrows()
    }
    perm_data = {pid: pd.read_csv(p, sep="\t")["g"].to_numpy()
                 for pid, p in perm_files.items() if p.exists()}

    for i, row in summary.iterrows():
        axis = next(a for pid, lab, a in PROGRAM_ORDER if lab == row["program_label"])
        col = AXIS_COL[axis]

        # violin/strip: subsample 200 of 1000 perms for clarity
        gs = perm_data.get(row["program_id"], np.array([]))
        if len(gs) > 0:
            idx = rng.choice(len(gs), size=min(200, len(gs)), replace=False)
            jitter = (rng.random(len(idx)) - 0.5) * 0.3
            ax.scatter(i + jitter, gs[idx], s=4,
                       facecolor="#9E9E9E", edgecolor="none",
                       alpha=0.45, zorder=2)
            # 95% CI of the permutation distribution
            ax.plot([i - 0.27, i + 0.27],
                    [row["perm_q025"], row["perm_q025"]],
                    color="#666666", lw=0.7, zorder=3)
            ax.plot([i - 0.27, i + 0.27],
                    [row["perm_q975"], row["perm_q975"]],
                    color="#666666", lw=0.7, zorder=3)
            ax.plot([i, i], [row["perm_q025"], row["perm_q975"]],
                    color="#999999", lw=0.5, zorder=2, ls=":")

        # real value as a coloured diamond
        ax.scatter([i], [row["real_g"]], marker="D", s=44,
                   facecolor=col, edgecolor="white", linewidth=0.8,
                   zorder=5)
        # p_emp text
        ax.text(i, row["real_g"] + 0.10,
                f"p$_\\mathrm{{emp}}$ = {row['p_emp_two_sided']:.2f}",
                fontsize=4.6, ha="center", va="bottom", color="#222222")

    ax.axhline(0, color="#888888", lw=0.5, ls=(0, (1.5, 1.5)), zorder=1)
    ax.set_xticks(range(n))
    ax.set_xticklabels(summary["program_label"].tolist(),
                       fontsize=5.6, rotation=22, ha="right")
    ax.set_ylabel("pooled g  (AF − SR)", fontsize=6.2, labelpad=2)
    ax.set_title(
        "S8a — size-matched random gene-set permutations",
        fontsize=6.8, pad=3, fontweight="bold")
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", labelsize=5.6, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)


# ---------------------------------------------------------------------
# Panel S8b — axis-membership permutation density (component B)
# ---------------------------------------------------------------------
def panel_S8b(ax):
    summary = pd.read_csv(RES / "figure_s8b_axis_perm.tsv", sep="\t")
    perm = pd.read_csv(RES / "figure_s8b_axis_perm_distribution.tsv.gz",
                       sep="\t")
    real_S = float(summary.loc[summary["axis"] == "STRUCTURAL", "real_g"].iloc[0])
    real_M = float(summary.loc[summary["axis"] == "METABOLIC", "real_g"].iloc[0])
    p_S = float(summary.loc[summary["axis"] == "STRUCTURAL",
                            "p_emp_two_sided"].iloc[0])
    p_M = float(summary.loc[summary["axis"] == "METABOLIC",
                            "p_emp_two_sided"].iloc[0])

    # Density on |g|
    abs_S = np.abs(perm["gS"].dropna())
    abs_M = np.abs(perm["gM"].dropna())

    bins = np.linspace(0, max(abs_S.max(), abs_M.max(),
                              abs(real_S), abs(real_M)) * 1.05, 28)
    ax.hist(abs_S, bins=bins, density=True, alpha=0.55,
            facecolor=AXIS_COL["STRUCTURAL"], edgecolor="white",
            linewidth=0.3, label="perm |g| STRUCTURAL")
    ax.hist(abs_M, bins=bins, density=True, alpha=0.55,
            facecolor=AXIS_COL["METABOLIC"], edgecolor="white",
            linewidth=0.3, label="perm |g| METABOLIC")
    ax.axvline(abs(real_S), color=AXIS_COL["STRUCTURAL"], lw=1.0)
    ax.axvline(abs(real_M), color=AXIS_COL["METABOLIC"], lw=1.0)
    ax.text(abs(real_S), ax.get_ylim()[1] * 0.92,
            f" real STRUCT |g| = {abs(real_S):.2f}\n p$_\\mathrm{{emp}}$ = {p_S:.2f}",
            fontsize=4.6, color=AXIS_COL["STRUCTURAL"],
            ha="left", va="top")
    ax.text(abs(real_M), ax.get_ylim()[1] * 0.55,
            f" real METAB |g| = {abs(real_M):.2f}\n p$_\\mathrm{{emp}}$ = {p_M:.2f}",
            fontsize=4.6, color=AXIS_COL["METABOLIC"],
            ha="left", va="top")
    ax.set_xlabel("|pooled g|  (axis composite)",
                  fontsize=6.2, labelpad=2)
    ax.set_ylabel("density", fontsize=6.2, labelpad=2)
    ax.set_title(
        "S8b — axis-membership permutation (B = 1000)",
        fontsize=6.8, pad=3, fontweight="bold")
    ax.tick_params(axis="both", labelsize=5.6, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.legend(loc="upper right", fontsize=4.8, frameon=False,
              handlelength=1.0, handleheight=0.5)


# ---------------------------------------------------------------------
# Panel S8c — real vs negative controls on (g, I², ρ)
# ---------------------------------------------------------------------
def panel_S8c(ax_g, ax_I2, ax_rho):
    real = pd.read_csv(RES / "figure6_meta_program.tsv", sep="\t")
    real = real[(real["type"] == "pooled") &
                (real["axis"].isin(["STRUCTURAL", "METABOLIC"]))].copy()
    real = real[["program_label", "axis", "g", "I2"]].rename(
        columns={"g": "pooled_g"})
    # LAA ρ
    laa = pd.read_csv(ROOT / "Project_P4" / "GSE41177" / "results"
                       / "duration_trend_laa_primary.tsv", sep="\t")
    laa = laa[["program_label", "statistic"]].rename(
        columns={"statistic": "duration_rho_LAA"})
    real = real.merge(laa, on="program_label", how="left")

    ctrl = pd.read_csv(RES / "figure_s8c_negative_controls.tsv", sep="\t")
    # Display the I² in percent and 0–100 scale
    real["I2_pct"] = real["I2"].astype(float)
    ctrl["I2_pct"] = ctrl["I2"].astype(float) * 100  # rma reports proportion

    ctrl["pretty"] = ctrl["pathway"].str.replace("HALLMARK_", "").str.replace(
        "REACTOME_", "").str.replace("_", " ").str.title()
    ctrl_color = ctrl["category"].map(
        {"axis-specificity": CTRL_AS_COL, "unrelated": CTRL_UN_COL})

    # All rows, ordered: real STRUCTURAL, real METABOLIC, ctrl A-spec, ctrl unrelated
    real_S = real[real["axis"] == "STRUCTURAL"].copy()
    real_M = real[real["axis"] == "METABOLIC"].copy()
    ctrl_A = ctrl[ctrl["category"] == "axis-specificity"].copy()
    ctrl_U = ctrl[ctrl["category"] == "unrelated"].copy()

    blocks = [
        ("real · STRUCTURAL", real_S, AXIS_COL["STRUCTURAL"], real_S["program_label"]),
        ("real · METABOLIC",  real_M, AXIS_COL["METABOLIC"],  real_M["program_label"]),
        ("ctrl · axis-spec",  ctrl_A, CTRL_AS_COL,            ctrl_A["pretty"]),
        ("ctrl · unrelated",  ctrl_U, CTRL_UN_COL,            ctrl_U["pretty"]),
    ]

    rows = []
    sep_idx = []
    cur = 0
    for label, df, col, names in blocks:
        for i, (_, r) in enumerate(df.iterrows()):
            rows.append({
                "y": cur,
                "g": r["pooled_g"], "I2": r.get("I2_pct", np.nan),
                "rho": r.get("duration_rho_LAA", np.nan),
                "name": list(names)[i] if i < len(names) else "",
                "color": col,
                "block": label,
            })
            cur += 1
        sep_idx.append(cur - 0.5)
        cur += 0.4

    rows = pd.DataFrame(rows)
    yticks = rows["y"].tolist()
    yticklabels = [f"{r['name']}" for _, r in rows.iterrows()]

    def _draw(ax, key, x_label, x_min=None, x_max=None,
              format_str="{:+.2f}", center0=True):
        if x_min is None: x_min = float(np.nanmin(rows[key]))
        if x_max is None: x_max = float(np.nanmax(rows[key]))
        pad = max(0.10 * (x_max - x_min), 0.05)
        ax.set_xlim(x_min - pad, x_max + pad)
        if center0:
            ax.axvline(0, color="#888888", lw=0.5,
                       ls=(0, (1.5, 1.5)), zorder=1)
        for sep in sep_idx[:-1]:
            ax.axhline(sep, color="#CCCCCC", lw=0.4, zorder=1)
        for _, r in rows.iterrows():
            if pd.isna(r[key]):
                continue
            ax.scatter(r[key], r["y"], marker="o", s=22,
                       facecolor=r["color"], edgecolor="white",
                       linewidth=0.4, zorder=4)
            ax.text(r[key], r["y"] + 0.30,
                    format_str.format(r[key]),
                    fontsize=4.4, ha="center", va="bottom",
                    color="#444444")
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=5.4)
        ax.set_xlabel(x_label, fontsize=6.0, labelpad=2)
        ax.tick_params(axis="x", labelsize=5.4, length=2)
        ax.tick_params(axis="y", length=0, pad=2)
        ax.set_ylim(rows["y"].max() + 1.0, -1.0)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
        ax.spines["left"].set_linewidth(0.5)
        ax.spines["bottom"].set_linewidth(0.5)

    _draw(ax_g,   "g",   "pooled g  (AF − SR)")
    _draw(ax_I2,  "I2",  "I² (%)", x_min=-5, x_max=100, center0=False,
          format_str="{:.0f}%")
    _draw(ax_rho, "rho", "LAA duration ρ")

    # Hide y-tick labels on the I² and ρ panels (shared with first)
    ax_I2.set_yticklabels([])
    ax_rho.set_yticklabels([])

    ax_g.set_title(
        "S8c — real frozen support programs vs pre-registered negative controls",
        fontsize=6.8, pad=3, fontweight="bold", loc="left")


# ---------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------
def panel_S8d(ax):
    """Leave-one-program-out (LOPO) axis robustness."""
    df = pd.read_csv(ROOT / "results" / "figure_s8d_lopo_axis.tsv", sep="\t")

    # Stack rows: STRUCTURAL FULL + 3 LOPO, then METABOLIC FULL + 3 LOPO
    rows = []
    for axis_lab in ["STRUCTURAL", "METABOLIC"]:
        sub = df[df["axis"] == axis_lab].copy()
        # FULL row first, then LOPO rows
        full_row = sub[sub["removed_program"] == "(none)"]
        lopo_rows = sub[sub["removed_program"] != "(none)"]
        for _, r in full_row.iterrows():
            rows.append({"axis": axis_lab, "label": "FULL axis",
                         "g": r["pooled_g"], "lo": r["ci_lo"],
                         "hi": r["ci_hi"], "I2": r["I2"],
                         "is_full": True})
        for _, r in lopo_rows.iterrows():
            rows.append({"axis": axis_lab,
                         "label": f"− {r['removed_label']}",
                         "g": r["pooled_g"], "lo": r["ci_lo"],
                         "hi": r["ci_hi"], "I2": r["I2"],
                         "is_full": False})
    rows = pd.DataFrame(rows)

    n = len(rows)
    y_axis = []
    cur = 0
    sep_at = None
    for i, r in rows.iterrows():
        y_axis.append(cur)
        if not r["is_full"]:
            cur += 1
        else:
            if i > 0:
                sep_at = cur - 0.5
                cur += 0.6
            cur += 1
    rows["y"] = y_axis

    for _, r in rows.iterrows():
        col = AXIS_COL[r["axis"]]
        ax.plot([r["lo"], r["hi"]], [r["y"], r["y"]],
                color=col, lw=1.0, alpha=0.95, zorder=3,
                solid_capstyle="round")
        if r["is_full"]:
            ax.scatter([r["g"]], [r["y"]], s=55,
                       facecolor=col, edgecolor="white", linewidth=0.7,
                       zorder=4)
        else:
            ax.scatter([r["g"]], [r["y"]], s=32,
                       facecolor="white", edgecolor=col, linewidth=1.0,
                       zorder=4)
        ax.text(r["hi"] + 0.04, r["y"],
                f"{r['g']:+.2f}  [{r['lo']:+.2f},{r['hi']:+.2f}]  I²={r['I2']:.0f}%",
                fontsize=4.6, ha="left", va="center",
                color="#222222", family="monospace")
    if sep_at is not None:
        ax.axhline(sep_at, color="#BBBBBB", lw=0.4, zorder=1)
    ax.axvline(0, color="#888888", lw=0.5, ls=(0, (1.5, 1.5)), zorder=1)

    ax.set_yticks(rows["y"].tolist())
    ax.set_yticklabels([
        f"{r['axis']}  ·  {r['label']}" if r["is_full"]
        else f"   {r['label']}"
        for _, r in rows.iterrows()
    ], fontsize=5.2)
    ax.invert_yaxis()
    ax.set_xlim(-2.0, 1.6)
    ax.set_xlabel("pooled g (AF − SR), REML / Hartung–Knapp",
                  fontsize=6.0, labelpad=2)
    ax.set_title(
        "S8d — leave-one-program-out (LOPO) axis robustness",
        fontsize=6.6, pad=3, fontweight="bold", loc="left")
    ax.tick_params(axis="x", labelsize=5.4, length=2)
    ax.tick_params(axis="y", length=0, pad=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["left"].set_linewidth(0.5); ax.spines["bottom"].set_linewidth(0.5)


def main() -> None:
    set_nature_style()
    fig_w = COL_DOUBLE
    fig_h = 9.6
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.07, right=0.97, top=0.955, bottom=0.05)

    gs = fig.add_gridspec(4, 3,
                          height_ratios=[0.95, 0.78, 1.45, 1.20],
                          width_ratios=[1.6, 1.0, 1.0],
                          hspace=0.55, wspace=0.40)

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :])
    ax_c1 = fig.add_subplot(gs[2, 0])
    ax_c2 = fig.add_subplot(gs[2, 1])
    ax_c3 = fig.add_subplot(gs[2, 2])
    ax_d = fig.add_subplot(gs[3, :])

    panel_S8a(ax_a)
    panel_S8b(ax_b)
    panel_S8c(ax_c1, ax_c2, ax_c3)
    panel_S8d(ax_d)

    fig.suptitle(
        "Specificity / falsification of the chamber-aware rewiring framework",
        fontsize=7.6, y=0.97, fontweight="bold",
    )
    fig.text(0.012, 0.97, "S8",
             fontsize=10, fontweight="bold", ha="left", va="top",
             transform=fig.transFigure)
    fig.text(0.5, 0.012,
             "Pre-registered controls (falsification pipeline; this script):  "
             "A) size-matched random gene sets · B) axis-membership permutation · "
             "C) six negative-control pathways (axis-specificity vs unrelated) · "
             "D) leave-one-program-out (LOPO) — both axes remain direction-preserved.",
             ha="center", va="bottom", fontsize=4.6,
             color="#666666", style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
