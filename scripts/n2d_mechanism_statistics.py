#!/usr/bin/env python3
"""N2d  Three quantitative statistics supporting the dual-axis sub-mechanism"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
GENES_TSV = ROOT / "results" / "n2_chamber_reversal_genes.tsv"
LONG_TSV  = ROOT / "results" / "n2_chamber_reversal_long.tsv"
SCORES_PRIMARY = ROOT / "Project_P3" / "GSE255612" / "program_scores_primary.tsv"

OUT_TSV   = ROOT / "results" / "n2d_mechanism_statistics.tsv"


# Canonical SLRP family (small leucine-rich proteoglycans).
SLRP_FAMILY = {
    # Class I
    "DCN", "BGN", "ASPN", "ECM2", "ECMX",
    # Class II
    "FMOD", "LUM", "PRELP", "KERA", "OMD",
    # Class III
    "EPYC", "OGN", "OPTC",
    # Class IV / V (less canonical but commonly grouped)
    "CHAD", "NYX", "TSKU", "PODN", "PODNL1",
}


def test_directional_asymmetry(genes: pd.DataFrame) -> dict:
    n_total = len(genes)
    n_la2ra = int((genes["healthy_direction"] == "LA-biased").sum())
    n_ra2la = int((genes["healthy_direction"] == "RA-biased").sum())
    minor = min(n_la2ra, n_ra2la)
    bt = stats.binomtest(minor, n_total, p=0.5, alternative="two-sided")
    return {
        "test":          "directional_asymmetry",
        "metric":        "binomtest (RA→LA vs LA→RA, p_null = 0.5)",
        "n_total":       n_total,
        "n_RA_to_LA":    n_ra2la,
        "n_LA_to_RA":    n_la2ra,
        "ratio":         f"{n_ra2la}:{n_la2ra}",
        "statistic":     bt.statistic,
        "p_value":       float(bt.pvalue),
    }


def test_slrp_enrichment(genes: pd.DataFrame, long_df: pd.DataFrame) -> dict:
    universe = set(long_df["gene"].astype(str).unique())
    foreground = set(genes["gene"].astype(str).tolist())
    universe |= foreground

    fg_in_slrp = foreground & SLRP_FAMILY
    bg_in_slrp = (universe - foreground) & SLRP_FAMILY

    a = len(fg_in_slrp)
    b = len(foreground) - a
    c = len(bg_in_slrp)
    d = (len(universe) - len(foreground)) - c

    # 2x2:  [[fg_slrp, fg_other], [bg_slrp, bg_other]]
    table = [[a, b], [c, d]]
    odds, p = stats.fisher_exact(table, alternative="greater")
    return {
        "test":          "SLRP_enrichment",
        "metric":        "Fisher's exact (one-sided greater)",
        "n_total":       len(universe),
        "fg_size":       len(foreground),
        "fg_in_SLRP":    a,
        "fg_SLRP_genes": ";".join(sorted(fg_in_slrp)),
        "bg_in_SLRP":    c,
        "odds_ratio":    float(odds),
        "p_value":       float(p),
    }


def test_bmp_bidirectionality() -> dict:
    """Use the singscore_v2 long table (donor × ct × program) to extract
    F11 BMP Δ per donor in adipocyte vs lymph endo, then test."""
    df = pd.read_csv(SCORES_PRIMARY, sep="\t")
    bmp_id = "REACTOME__REACTOME_SIGNALING_BY_BMP"
    sub = df[df["program_id"] == bmp_id]
    cells = ["adipocyte", "lymphatic_endothelial"]
    sub = sub[sub["compartment"].isin(cells)]

    rec = {}
    for ct in cells:
        ss = sub[sub["compartment"] == ct]
        af = ss.loc[ss["rhythm"] == "AF",   "score"].to_numpy()
        ct_ = ss.loc[ss["rhythm"] == "Ctrl", "score"].to_numpy()
        rec[f"{ct}_AF_mean"]  = float(af.mean())
        rec[f"{ct}_Ctl_mean"] = float(ct_.mean())
        rec[f"{ct}_delta"]    = float(af.mean() - ct_.mean())
        rec[f"{ct}_n_AF"]     = int(len(af))
        rec[f"{ct}_n_Ctl"]    = int(len(ct_))

    # Welch t-test:  is Δ_adipocyte significantly different from Δ_lymph?
    # Combine the two compartments AF / Ctl values per donor.  Build per-
    # donor Δ vectors:  for each donor that has both AF / Ctl in the same
    # compartment, take the difference of mean singscore vs the cohort mean
    # (use single-sample residuals, then test paired across compartments).
    # Simpler:  compare per-cell Δ distributions via Welch t on the donor
    # AF singscores referenced to the per-cell Ctl mean.

    def per_donor_delta(ct):
        ss = sub[sub["compartment"] == ct]
        # use AF singscores minus the Ctl mean of the same cell type
        ct_mean = ss.loc[ss["rhythm"] == "Ctrl", "score"].mean()
        af = ss.loc[ss["rhythm"] == "AF",
                       ["donor_id", "score"]].copy()
        af["delta"] = af["score"] - ct_mean
        return af.set_index("donor_id")["delta"]

    da = per_donor_delta("adipocyte")
    dl = per_donor_delta("lymphatic_endothelial")
    common = sorted(set(da.index) & set(dl.index))
    da_c = da.loc[common].to_numpy()
    dl_c = dl.loc[common].to_numpy()

    # paired Wilcoxon (sign-test variant) on Δ_adipo − Δ_lymph
    if len(common) >= 5:
        diff = da_c - dl_c
        wilc = stats.wilcoxon(diff, alternative="two-sided")
        wilc_p = float(wilc.pvalue)
        n_pos = int((diff > 0).sum()); n_neg = int((diff < 0).sum())
        # binomtest on sign
        bt = stats.binomtest(n_pos, n_pos + n_neg, p=0.5,
                                alternative="two-sided")
        sign_p = float(bt.pvalue)
    else:
        wilc_p = float("nan"); sign_p = float("nan")
        n_pos = n_neg = 0

    rec.update({
        "test":           "BMP_bidirectionality",
        "metric":         "paired Wilcoxon  +  sign-test on  Δ_adipo − Δ_lymph",
        "n_paired_donors": len(common),
        "n_pos":           n_pos,
        "n_neg":           n_neg,
        "wilcoxon_p":      wilc_p,
        "sign_p":          sign_p,
    })
    return rec


def main() -> None:
    genes  = pd.read_csv(GENES_TSV, sep="\t")
    long   = pd.read_csv(LONG_TSV, sep="\t")

    rows = [
        test_directional_asymmetry(genes),
        test_slrp_enrichment(genes, long),
        test_bmp_bidirectionality(),
    ]
    out = pd.DataFrame(rows)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)}")
    for r in rows:
        print()
        for k, v in r.items():
            print(f"  {k:22s}  {v}")


if __name__ == "__main__":
    main()
