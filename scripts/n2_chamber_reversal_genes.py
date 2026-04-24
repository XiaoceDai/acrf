#!/usr/bin/env python3
"""
N2 — Discover chamber-identity REVERSAL genes.

Conceptual contribution:
    A reversal gene is a gene whose chamber-identity in healthy donors
    (LA-biased or RA-biased) is *flipped* in atrial fibrillation (AF)
    donors of paired-chamber cohorts, while sinus-rhythm (SR) donors
    of the same cohorts preserve the healthy direction.

    Sample-level limma cannot identify reversal genes because it does
    not ask the chamber-direction question, and pathway-GSEA does not
    expose the per-gene LA/RA flip.  This object class is uniquely
    enabled by the donor-level chamber-aware metric of CRF.

Inputs (gene-symbol space; log-scale or auto-logged here):
    bulk/Project_P01/GSE112339/GSE112339_HTSeq_hg19_Ensembl.txt.gz
        (healthy reference; 8 paired LA / RA donors)
    Project_P1/GSE128188/expr_symbol.tsv  + sample_metadata.tsv
    Project_P1/GSE79768/expr_symbol.tsv   + sample_metadata.tsv
    Project_P1/GSE115574/expr_symbol.tsv  + sample_metadata.tsv

Per-cohort, per-donor:
    aCAI_g(d) = log2(score_LA(d, g)) - log2(score_RA(d, g))
For a paired cohort, AF_dir(g) = sign(mean(aCAI_g | AF donors)),
SR_dir(g) = sign(mean(aCAI_g | SR donors));
Healthy_dir(g) = sign(mean(aCAI_g | 8 healthy donors)).

Reversal call per gene per cohort:
    af_reversed       : Healthy_dir x AF_dir < 0
    sr_preserved      : Healthy_dir x SR_dir > 0   (excludes batch-only flips)
    af_specific_flip  : af_reversed AND sr_preserved
    above_magnitude   : |aCAI_healthy| >= MAG_HEALTHY AND |aCAI_AF| >= MAG_AF
Final reversal gene call:
    af_specific_flip AND above_magnitude in >= MIN_COHORTS / 3 cohorts.

Outputs:
    results/n2_chamber_reversal_long.tsv         (gene x cohort rows)
    results/n2_chamber_reversal_genes.tsv        (filtered final list)
    figures/figN2_chamber_reversal.pdf           (summary panel)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import (
    COL_DOUBLE, LA_COLOR, RA_COLOR, panel_label,
    set_nature_style, strip_spines,
)

set_nature_style()

OUT_LONG  = ROOT / "results" / "n2_chamber_reversal_long.tsv"
OUT_GENES = ROOT / "results" / "n2_chamber_reversal_genes.tsv"
FIG_OUT   = ROOT / "figures" / "figN2_chamber_reversal.pdf"

# --- thresholds ----------------------------------------------------------
# Calibrated against per-cohort empirical |aCAI| distributions
# (see results/n2_chamber_reversal_long.tsv summary).
# |aCAI_healthy| >= 0.3 log2 captures the top ~15 % of chamber-biased
# genes in the healthy reference (~1.23-fold chamber bias and stronger).
# AF reversal magnitudes are typically smaller than healthy bias
# magnitudes; |aCAI_AF| >= 0.15 log2 (~1.11-fold reversed bias) is set
# at roughly the per-cohort 80th percentile of |aCAI_AF|.
MAG_HEALTHY = 0.3
MAG_AF      = 0.15
MIN_COHORTS = 2       # require AF-specific reversal in >= 2 of 3 paired cohorts


# ---------------------------------------------------------------------------
# Loaders.
# ---------------------------------------------------------------------------
def load_healthy_GSE112339() -> tuple[pd.DataFrame, list[tuple[str, str]]]:
    """Return (gene x sample) log2 matrix and list of (LA_sample, RA_sample)
    paired donor tuples."""
    f = ROOT / "bulk" / "Project_P01" / "GSE112339" \
        / "GSE112339_HTSeq_hg19_Ensembl.txt.gz"
    raw = pd.read_csv(f, sep="\t", index_col=0)
    # collapse duplicated gene symbols by max-mean (rare)
    raw = raw.groupby(raw.index).max()
    # restrict to LA / RA samples only
    keep = [c for c in raw.columns if c.startswith("LA_") or c.startswith("RA_")]
    raw = raw[keep]
    # log2-CPM-like transform on counts
    cpm = (raw.div(raw.sum(axis=0), axis=1) * 1e6)
    log = np.log2(cpm + 1.0)
    # build LA<->RA pairs by donor suffix
    la = sorted([c for c in log.columns if c.startswith("LA_")])
    ra = sorted([c for c in log.columns if c.startswith("RA_")])
    pairs = []
    for a in la:
        donor = a[3:]
        b = f"RA_{donor}"
        if b in ra:
            pairs.append((a, b))
    return log, pairs


def load_paired_cohort(cohort: str) -> tuple[
        pd.DataFrame, list[tuple[str, str, str]]]:
    """Return (gene x sample) log2 matrix and list of
    (LA_gsm, RA_gsm, rhythm) per donor that has both chambers present."""
    expr = pd.read_csv(ROOT / "Project_P1" / cohort / "expr_symbol.tsv",
                       sep="\t", index_col=0)
    expr = expr.groupby(expr.index).max()
    md = pd.read_csv(ROOT / "Project_P1" / cohort / "sample_metadata.tsv",
                     sep="\t")
    md = md[md["gsm"].isin(expr.columns)]

    # Auto-log if not already log-scale (RNA-seq counts will have huge max)
    if expr.values.max() > 100:
        expr = np.log2(expr.clip(lower=0) + 1.0)

    # build paired LA/RA per donor with rhythm
    pairs = []
    for donor, sub in md.groupby("donor_id"):
        la_rows = sub[sub["chamber"] == "LA"]
        ra_rows = sub[sub["chamber"] == "RA"]
        if len(la_rows) >= 1 and len(ra_rows) >= 1:
            rhythm = la_rows["rhythm"].iloc[0]
            pairs.append((la_rows["gsm"].iloc[0],
                          ra_rows["gsm"].iloc[0], rhythm))
    return expr, pairs


def per_donor_aCAI(expr: pd.DataFrame,
                   pairs: list[tuple[str, ...]]) -> pd.DataFrame:
    """gene x donor matrix of (log_LA - log_RA) per pair."""
    cols = {}
    for p in pairs:
        la, ra = p[0], p[1]
        if la in expr.columns and ra in expr.columns:
            label = f"{la}__{ra}"
            cols[label] = expr[la] - expr[ra]
    return pd.DataFrame(cols, index=expr.index)


# ---------------------------------------------------------------------------
# Core analysis.
# ---------------------------------------------------------------------------
def main() -> None:
    print("Loading healthy GSE112339 ...")
    h_expr, h_pairs = load_healthy_GSE112339()
    h_aCAI = per_donor_aCAI(h_expr, h_pairs)
    h_mean = h_aCAI.mean(axis=1)
    h_dir  = np.sign(h_mean)
    print(f"  GSE112339 healthy: {h_expr.shape[1]} samples, "
          f"{len(h_pairs)} paired donors, "
          f"{h_aCAI.shape[0]} genes")

    paired = {}
    for cohort in ["GSE128188", "GSE79768", "GSE115574"]:
        print(f"Loading paired cohort {cohort} ...")
        expr, pairs = load_paired_cohort(cohort)
        aCAI = per_donor_aCAI(expr, pairs)
        af_pairs = [p for p in pairs if p[2] == "AF"]
        sr_pairs = [p for p in pairs if p[2] == "SR"]
        af_cols = [f"{p[0]}__{p[1]}" for p in af_pairs
                   if f"{p[0]}__{p[1]}" in aCAI.columns]
        sr_cols = [f"{p[0]}__{p[1]}" for p in sr_pairs
                   if f"{p[0]}__{p[1]}" in aCAI.columns]
        paired[cohort] = {
            "aCAI": aCAI,
            "af_cols": af_cols, "sr_cols": sr_cols,
            "n_AF": len(af_cols), "n_SR": len(sr_cols),
        }
        print(f"  {cohort}: {len(af_cols)} AF + {len(sr_cols)} SR paired donors")

    # gene set common to healthy + all 3 paired cohorts
    common = set(h_aCAI.index)
    for c in paired.values():
        common &= set(c["aCAI"].index)
    common = sorted(common)
    print(f"\nCommon gene symbols across all 4 cohorts: {len(common)}")

    h_mean = h_mean.loc[common]
    h_dir  = h_dir.loc[common]

    # build long table
    rows = []
    for cohort, blob in paired.items():
        a = blob["aCAI"].loc[common]
        af_mean = a[blob["af_cols"]].mean(axis=1) if blob["af_cols"] else pd.Series(np.nan, index=common)
        sr_mean = a[blob["sr_cols"]].mean(axis=1) if blob["sr_cols"] else pd.Series(np.nan, index=common)
        af_dir = np.sign(af_mean)
        sr_dir = np.sign(sr_mean)

        af_reversed     = (h_dir * af_dir < 0)
        sr_preserved    = (h_dir * sr_dir > 0)
        af_specific     = af_reversed & sr_preserved
        above_magnitude = (h_mean.abs() >= MAG_HEALTHY) & (af_mean.abs() >= MAG_AF)
        call            = af_specific & above_magnitude

        df = pd.DataFrame({
            "gene":           common,
            "cohort":         cohort,
            "n_AF_donors":    blob["n_AF"],
            "n_SR_donors":    blob["n_SR"],
            "aCAI_healthy":   h_mean.values,
            "healthy_dir":    h_dir.values.astype(int),
            "aCAI_AF":        af_mean.values,
            "AF_dir":         af_dir.values.astype(int),
            "aCAI_SR":        sr_mean.values,
            "SR_dir":         sr_dir.values.astype(int),
            "af_reversed":    af_reversed.values,
            "sr_preserved":   sr_preserved.values,
            "af_specific":    af_specific.values,
            "above_magnitude": above_magnitude.values,
            "call_in_cohort": call.values,
        })
        rows.append(df)

    long = pd.concat(rows, ignore_index=True)
    OUT_LONG.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(OUT_LONG, sep="\t", index=False)
    print(f"wrote {OUT_LONG.relative_to(ROOT)} ({len(long)} rows)")

    # cross-cohort consistency: count cohorts where call is True per gene
    n_cohorts_call = (long.groupby("gene")["call_in_cohort"].sum()
                          .astype(int).rename("n_cohorts_call"))

    # also total magnitude across cohorts (mean |aCAI_AF|)
    af_mag = (long.groupby("gene")["aCAI_AF"].apply(lambda x: x.abs().mean())
                  .rename("mean_abs_aCAI_AF"))
    h_mag  = h_mean.abs().reindex(n_cohorts_call.index).rename("abs_aCAI_healthy")
    h_dir_per_gene = h_dir.reindex(n_cohorts_call.index)
    healthy_label = h_dir_per_gene.map({1.0: "LA-biased",
                                         -1.0: "RA-biased",
                                         0.0: "ambiguous"})

    # the final reversal gene list
    keep_genes = [str(g) for g in
                  n_cohorts_call[n_cohorts_call >= MIN_COHORTS].index.tolist()]
    final = pd.DataFrame({
        "gene": pd.Series(keep_genes, dtype=object),
        "n_cohorts_call":  n_cohorts_call.loc[keep_genes].values,
        "healthy_direction": healthy_label.loc[keep_genes].values,
        "aCAI_healthy":  h_mean.loc[keep_genes].round(3).values,
        "mean_abs_aCAI_AF": af_mag.loc[keep_genes].round(3).values,
    })
    # add per-cohort AF-mean aCAI columns
    for cohort in paired.keys():
        sub = (long.query("cohort == @cohort")
                  .set_index("gene")
                  .reindex(keep_genes)
                  ["aCAI_AF"].round(3)
                  .rename(f"aCAI_AF_{cohort}")
                  .reset_index())
        sub["gene"] = sub["gene"].astype(str)
        final = final.merge(sub, on="gene", how="left")
    # sort by mean |AF rewiring magnitude|
    final = final.sort_values(
        ["n_cohorts_call", "mean_abs_aCAI_AF"],
        ascending=[False, False],
    ).reset_index(drop=True)

    OUT_GENES.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUT_GENES, sep="\t", index=False)
    print(f"wrote {OUT_GENES.relative_to(ROOT)} "
          f"({len(final)} reversal genes; "
          f"{(final['healthy_direction']=='LA-biased').sum()} LA→RA, "
          f"{(final['healthy_direction']=='RA-biased').sum()} RA→LA)")

    # ---- summary figure --------------------------------------------------
    fig = plt.figure(figsize=(COL_DOUBLE, 0.55 * COL_DOUBLE))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.55,
                          width_ratios=[1.05, 1.05, 1.0])

    # panel a: healthy vs AF aCAI per gene scatter (one panel per cohort)
    for k, cohort in enumerate(paired.keys()):
        ax = fig.add_subplot(gs[0, k])
        sub = long.query("cohort == @cohort").set_index("gene").loc[common]
        x = sub["aCAI_healthy"].values
        y = sub["aCAI_AF"].values
        # background grey
        ax.scatter(x, y, s=1.2, color="#CCCCCC", alpha=0.4,
                   edgecolor="none", zorder=1)
        # cross-cohort reversal genes highlighted
        is_kept = sub.index.isin(keep_genes)
        is_la_to_ra = is_kept & (sub["healthy_dir"] > 0)
        is_ra_to_la = is_kept & (sub["healthy_dir"] < 0)
        ax.scatter(x[is_la_to_ra], y[is_la_to_ra],
                   s=4, color=LA_COLOR, alpha=0.85,
                   edgecolor="white", linewidth=0.15, zorder=3,
                   label="LA→RA")
        ax.scatter(x[is_ra_to_la], y[is_ra_to_la],
                   s=4, color=RA_COLOR, alpha=0.85,
                   edgecolor="white", linewidth=0.15, zorder=3,
                   label="RA→LA")
        # axes guides
        lim = max(np.percentile(np.abs(x), 99.5),
                  np.percentile(np.abs(y), 99.5))
        ax.axhline(0, color="#888888", lw=0.4, zorder=0)
        ax.axvline(0, color="#888888", lw=0.4, zorder=0)
        ax.plot([-lim, lim], [-lim, lim], color="#888888",
                lw=0.4, ls="--", zorder=0)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xlabel("aCAI healthy (log2 LA − RA)")
        ax.set_ylabel("aCAI AF (log2 LA − RA)")
        n_af = paired[cohort]["n_AF"]; n_sr = paired[cohort]["n_SR"]
        ax.set_title(f"{cohort}  (AF n={n_af}, SR n={n_sr})", pad=2)
        ax.legend(loc="upper left", frameon=False, fontsize=5)
        strip_spines(ax)
        if k == 0:
            panel_label(ax, "a", x=-0.18, y=1.06)

    # panel b: top-N reversal genes ranked bar
    ax = fig.add_subplot(gs[1, :2])
    n_show = min(40, len(final))
    top = final.head(n_show).copy()
    y_pos = np.arange(n_show)[::-1]
    colors = [LA_COLOR if d == "LA-biased" else RA_COLOR
              for d in top["healthy_direction"]]
    ax.barh(y_pos, top["mean_abs_aCAI_AF"], color=colors,
            edgecolor="white", linewidth=0.2)
    for yi, gene, n_c in zip(y_pos, top["gene"], top["n_cohorts_call"]):
        ax.text(0.02, yi, f"{gene}", fontsize=5,
                ha="left", va="center", color="white",
                fontweight="bold")
        ax.text(top.iloc[len(y_pos) - yi - 1]["mean_abs_aCAI_AF"] + 0.05,
                yi, f"({n_c}/3)", fontsize=4.5,
                ha="left", va="center", color="#444444")
    ax.set_yticks([]); ax.set_xlabel("mean |aCAI AF| across paired cohorts")
    ax.set_title(
        f"b  Top {n_show} chamber-reversal genes (cross-cohort consistent in "
        f"≥ {MIN_COHORTS}/3 cohorts)",
        loc="left", pad=4,
    )
    strip_spines(ax)

    # panel c: counts per healthy_direction × n_cohorts_call
    ax = fig.add_subplot(gs[1, 2])
    counts = (final.groupby(["healthy_direction", "n_cohorts_call"]).size()
                   .unstack(fill_value=0).T)
    cohorts_order = sorted(counts.index.tolist())
    x = np.arange(len(cohorts_order))
    width = 0.35
    if "LA-biased" in counts.columns:
        ax.bar(x - width/2, counts["LA-biased"].reindex(cohorts_order, fill_value=0),
               width=width, color=LA_COLOR, label="LA→RA")
    if "RA-biased" in counts.columns:
        ax.bar(x + width/2, counts["RA-biased"].reindex(cohorts_order, fill_value=0),
               width=width, color=RA_COLOR, label="RA→LA")
    ax.set_xticks(x); ax.set_xticklabels(cohorts_order)
    ax.set_xlabel("cross-cohort consistency (n cohorts)")
    ax.set_ylabel("reversal genes (n)")
    ax.set_title("c  Direction × consistency", loc="left", pad=4)
    ax.legend(loc="upper right", frameon=False, fontsize=5)
    strip_spines(ax)

    fig.suptitle(
        "Chamber-identity reversal genes — a CRF-specific discovery class",
        y=1.00, fontsize=8,
    )

    FIG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT)
    plt.close(fig)
    print(f"wrote {FIG_OUT.relative_to(ROOT)}")

    # ---- console summary -------------------------------------------------
    print("\n=== N2 summary ===")
    print(f"genes in healthy GSE112339              : {len(h_aCAI)}")
    print(f"genes common across 4 cohorts            : {len(common)}")
    print(f"chamber-reversal call (n_cohorts >= {MIN_COHORTS}): {len(final)}")
    if len(final):
        print(f"  LA→RA reversal: {(final['healthy_direction']=='LA-biased').sum()}")
        print(f"  RA→LA reversal: {(final['healthy_direction']=='RA-biased').sum()}")
        print(f"  3/3 consistent : {(final['n_cohorts_call']==3).sum()}")
        print(f"\nTop 15 by mean |aCAI AF|:")
        print(final.head(15)[["gene", "healthy_direction",
                              "n_cohorts_call", "aCAI_healthy",
                              "mean_abs_aCAI_AF"]].to_string(index=False))


if __name__ == "__main__":
    main()
