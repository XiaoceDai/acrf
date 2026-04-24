#!/usr/bin/env python3
"""
Figure 4 (REPLACEMENT) — Matrix-invasion chamber-reversal signature
intensifies with atrial fibrillation duration in the left atrial appendage.

Concept:
    The 24 chamber-identity reversal genes identified in Fig. 6c
    define a gene-level signature whose AF-specific reversal direction
    is known per gene (RA -> LA for 21 genes, LA -> RA for 3 genes).
    If matrix-invasion is a progressive process, the LA expression of
    RA -> LA reversal genes should INCREASE with AF duration in LAA.

Procedure:
    1. Take GSE41177 LAA AF samples with numeric AF duration in months
       (n = 12; duration 10-150 months).
    2. Per gene, compute Spearman rho of LA expression vs AF duration.
    3. Sign-align by expected AF direction:
           RA -> LA (n = 21): positive rho = invasion intensifies
           LA -> RA (n = 3) : negative rho = invasion intensifies
       Aligned rho is positive when the chamber-reversal direction
       intensifies with duration.
    4. Composite signature rho = mean of sign-aligned rhos across the
       24 reversal genes.
    5. Empirical null: 1000 size-matched random gene-set permutations.
    6. Report signature rho, empirical p, and per-gene contributions.

Outputs:
    results/figure4_reversal_chronicity_per_gene.tsv
    results/figure4_reversal_chronicity_summary.tsv
    figures/fig4_reversal_chronicity.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (
    COL_DOUBLE, LA_COLOR, RA_COLOR, panel_label,
    set_nature_style, strip_spines,
)

set_nature_style()

EXPR_PATH = ROOT / "Project_P4" / "GSE41177" / "expr_symbol_LAA.tsv"
META_PATH = ROOT / "Project_P4" / "GSE41177" / "sample_metadata.tsv"
REV_PATH  = ROOT / "results" / "n2_chamber_reversal_genes.tsv"

OUT_PG    = ROOT / "results" / "figure4_reversal_chronicity_per_gene.tsv"
OUT_SUMM  = ROOT / "results" / "figure4_reversal_chronicity_summary.tsv"
FIG_OUT   = ROOT / "figures" / "fig4_reversal_chronicity.pdf"

N_PERM = 1000
RNG_SEED = 0


def main() -> None:
    expr = pd.read_csv(EXPR_PATH, sep="\t", index_col=0)
    md = pd.read_csv(META_PATH, sep="\t")
    md_laa = md[md["tissue"] == "LAA"].copy()
    md_laa["duration_months"] = pd.to_numeric(md_laa["duration_months"],
                                               errors="coerce")
    af_laa = md_laa[(md_laa["rhythm"] == "AF")
                     & md_laa["duration_months"].notna()
                     & md_laa["gsm"].isin(expr.columns)]
    sr_laa = md_laa[(md_laa["rhythm"] == "SR")
                     & md_laa["gsm"].isin(expr.columns)]
    print(f"AF LAA samples with numeric duration: {len(af_laa)}")
    print(f"SR LAA visual anchors                : {len(sr_laa)}")

    af_gsms = af_laa["gsm"].tolist()
    durations = af_laa.set_index("gsm")["duration_months"]
    af_expr = expr[af_gsms]

    rev = pd.read_csv(REV_PATH, sep="\t")
    rev["expected_dir_in_AF"] = np.where(
        rev["healthy_direction"] == "RA-biased", +1, -1
    )
    print(f"Reversal genes loaded: {len(rev)}")

    # per-gene rho vs duration
    rows = []
    for _, r in rev.iterrows():
        g = r["gene"]
        if g not in af_expr.index:
            continue
        x = durations.values
        y = af_expr.loc[g].astype(float).values
        rho, p = spearmanr(x, y, nan_policy="omit")
        rows.append({
            "gene":              g,
            "healthy_direction": r["healthy_direction"],
            "expected_dir_in_AF": int(r["expected_dir_in_AF"]),
            "rho_vs_duration":   round(float(rho), 4),
            "p_spearman":        round(float(p), 4),
            "rho_aligned":       round(float(rho) * int(r["expected_dir_in_AF"]), 4),
            "n_AF":              len(af_gsms),
        })
    pg = pd.DataFrame(rows)
    pg = pg.sort_values("rho_aligned", ascending=False).reset_index(drop=True)
    OUT_PG.parent.mkdir(parents=True, exist_ok=True)
    pg.to_csv(OUT_PG, sep="\t", index=False)
    print(f"wrote {OUT_PG.relative_to(ROOT)}")

    # composite signature rho
    obs_signature = float(pg["rho_aligned"].mean())
    obs_frac_pos  = float((pg["rho_aligned"] > 0).mean())
    n_signature   = len(pg)

    # 1000-permutation null: pick random size-matched gene set, compute
    # mean of *unaligned* spearman rhos signed by expected direction;
    # since random genes have no expected direction, we sample directions
    # randomly from the same {+1, -1} multiset to be conservative.
    print(f"Running {N_PERM} size-matched random gene permutations ...")
    rng = np.random.default_rng(RNG_SEED)
    all_genes = af_expr.index.tolist()
    # pre-compute rho for all genes once for speed
    print("  computing all-gene rhos for null pool ...")
    all_rho = {}
    x = durations.values.astype(float)
    Y = af_expr.values.astype(float)  # gene x sample
    # row-wise spearman: rank x and each y row
    rx = pd.Series(x).rank().values
    rY = pd.DataFrame(Y).rank(axis=1).values
    n = len(x)
    rx_c = rx - rx.mean()
    rY_c = rY - rY.mean(axis=1, keepdims=True)
    num = (rY_c * rx_c).sum(axis=1)
    denom = np.sqrt((rY_c**2).sum(axis=1) * (rx_c**2).sum())
    rho_all = np.where(denom > 0, num / denom, 0.0)
    rho_pool = pd.Series(rho_all, index=af_expr.index)
    rho_pool = rho_pool.replace([np.inf, -np.inf], 0.0).fillna(0.0)

    expected_signs = pg["expected_dir_in_AF"].values
    null_signatures = np.zeros(N_PERM)
    pool_index = rho_pool.index.tolist()
    pool_size = len(pool_index)
    sig_size = n_signature
    rng_idx = rng.integers(0, pool_size, size=(N_PERM, sig_size))
    for k in range(N_PERM):
        # sample sig_size genes with replacement (large pool, near
        # without-replacement effectively)
        idx = rng_idx[k]
        random_rhos = rho_pool.iloc[idx].values
        # random direction labels matching observed direction multiset
        rng.shuffle(expected_signs)  # in-place shuffle for randomness
        null_signatures[k] = float((random_rhos * expected_signs).mean())
    emp_p = float((null_signatures >= obs_signature).sum() + 1) / (N_PERM + 1)

    summary = pd.DataFrame([{
        "n_signature_genes":     n_signature,
        "n_AF_LAA_samples":      len(af_gsms),
        "duration_min_months":   float(durations.min()),
        "duration_max_months":   float(durations.max()),
        "signature_rho_aligned": round(obs_signature, 4),
        "frac_genes_aligned_positive": round(obs_frac_pos, 3),
        "n_perm":                N_PERM,
        "empirical_p":           round(emp_p, 4),
        "null_mean":             round(float(null_signatures.mean()), 4),
        "null_sd":               round(float(null_signatures.std()), 4),
    }])
    summary.to_csv(OUT_SUMM, sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"wrote {OUT_SUMM.relative_to(ROOT)}")

    # ---- figure -----------------------------------------------------------
    fig = plt.figure(figsize=(COL_DOUBLE, 0.5 * COL_DOUBLE))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.50,
                          height_ratios=[1.0, 1.0])

    # panel a: schematic + composite signature ρ + empirical null histogram
    ax = fig.add_subplot(gs[0, 0])
    ax.hist(null_signatures, bins=40, color="#CCCCCC",
            edgecolor="white", linewidth=0.3)
    ax.axvline(obs_signature, color="#B2182B", lw=1.5, ymax=0.95,
               label=f"observed = {obs_signature:.3f}")
    ax.set_xlabel("composite signature ρ (sign-aligned)")
    ax.set_ylabel("permutations")
    ax.legend(loc="upper left", frameon=False, fontsize=5)
    ax.set_title(
        f"a  Matrix-invasion signature ρ vs random null\n"
        f"empirical p = {emp_p:.3f} ({N_PERM} permutations)",
        loc="left", pad=4, fontsize=6.4,
    )
    strip_spines(ax)
    panel_label(ax, "a", x=-0.20, y=1.06)

    # panel b: per-gene aligned ρ ranked
    ax = fig.add_subplot(gs[0, 1:])
    n_show = len(pg)
    y = np.arange(n_show)[::-1]
    colors = [RA_COLOR if d == "RA-biased" else LA_COLOR
              for d in pg["healthy_direction"]]
    ax.barh(y, pg["rho_aligned"], color=colors,
            edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="#444", lw=0.6)
    ax.axvline(obs_signature, color="#B2182B", ls="--", lw=0.8,
               label=f"signature mean = {obs_signature:.3f}")
    for yi, gene, ra in zip(y, pg["gene"], pg["rho_aligned"]):
        ax.text(-0.02 if ra >= 0 else 0.02, yi, gene,
                fontsize=4.6, ha="right" if ra >= 0 else "left",
                va="center", color="#1A1A1A")
    ax.set_yticks([]); ax.set_xlabel("Spearman ρ vs AF duration (sign-aligned)")
    ax.set_xlim(min(pg["rho_aligned"].min(), -0.6) - 0.1,
                max(pg["rho_aligned"].max(), 0.6) + 0.1)
    ax.legend(loc="lower right", frameon=False, fontsize=5)
    ax.set_title(
        f"b  Per-gene contributions to the matrix-invasion signature "
        f"({(pg['rho_aligned']>0).sum()}/{len(pg)} aligned positive; "
        f"red = RA→LA matrix gene, blue = LA→RA)",
        loc="left", pad=4, fontsize=6.4,
    )
    strip_spines(ax)
    panel_label(ax, "b", x=-0.04, y=1.06)

    # panel c (bottom row): example genes - top 3 RA->LA aligned positive
    top = pg.head(3)
    sr_gsms = sr_laa[sr_laa["gsm"].isin(expr.columns)]["gsm"].tolist()
    for k, (_, row) in enumerate(top.iterrows()):
        ax = fig.add_subplot(gs[1, k])
        g = row["gene"]
        af_y = expr.loc[g, af_gsms].astype(float).values
        sr_y = expr.loc[g, sr_gsms].astype(float).values
        ax.scatter([0]*len(sr_y), sr_y, color=LA_COLOR, s=12,
                   edgecolor="white", linewidth=0.3, alpha=0.6,
                   label="SR (anchor at 0)")
        ax.scatter(durations.values, af_y, color=RA_COLOR, s=18,
                   edgecolor="white", linewidth=0.3, label="AF")
        # linear fit on AF only
        mask = ~np.isnan(durations.values) & ~np.isnan(af_y)
        if mask.sum() >= 3:
            xs = durations.values[mask]; ys = af_y[mask]
            slope, intercept = np.polyfit(xs, ys, 1)
            xline = np.linspace(xs.min(), xs.max(), 50)
            ax.plot(xline, slope * xline + intercept, color="#444", lw=0.8)
        ax.set_xlabel("AF duration (months)")
        if k == 0:
            ax.set_ylabel("LA expression (RMA log2)")
        ax.set_title(
            f"c{['', '', ''][k]}  {g}: aligned ρ = {row['rho_aligned']:+.2f}",
            loc="left", pad=2, fontsize=6.2,
        )
        if k == 0:
            ax.legend(loc="upper left", frameon=False, fontsize=5)
        strip_spines(ax)
    panel_label(fig.axes[2], "c", x=-0.20, y=1.06)

    fig.suptitle(
        "Matrix-invasion chamber-reversal signature intensifies with AF "
        "duration in the LA appendage",
        y=1.00, fontsize=7.2, fontweight="bold",
    )

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
