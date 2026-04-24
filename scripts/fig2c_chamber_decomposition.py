#!/usr/bin/env python3
"""
Figure 2C — Rewiring state map (conclusion-layer classification).

This panel converts continuous donor-level aCAI signals (Fig 2B) into
discrete, interpretable rewiring states for each family × cohort pair:
    Preserved / Attenuated / Reversed / Amplified

Inputs:
    results/fig2_primary_acai_long.tsv
    results/fig2_primary_axis_families.tsv
Outputs:
    figures/fig2c_chamber_decomposition.pdf
    results/fig2c_chamber_decomposition.tsv
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
DEFAULT_FAM_TSV = ROOT / "results" / "fig2_primary_axis_families.tsv"
DEFAULT_ATTR_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
DEFAULT_OUT_TSV = ROOT / "results" / "fig2c_chamber_decomposition.tsv"
DEFAULT_FIG_OUT = ROOT / "figures" / "fig2c_chamber_decomposition.pdf"

# NOTE: Fig 2C row order intentionally uses only Fig 1F-derived metadata
# (family_id within direction).  Sorting by Integrated Atlas PriorityScore
# would leak Fig 2D / 2E results into Fig 2C, breaking the figure's narrative
# order; do not reintroduce that sort key here.

LA_COLOR = "#C97A74"  # LA = red
RA_COLOR = "#7D99B8"  # RA = blue
GREY_MID = "#9C9C9C"

STATE_COLORS = {
    "Preserved": "#A9C9B5",
    "Attenuated": "#E9C27B",
    "Reversed": "#C6A2D7",
    "Amplified": "#E59595",
}
STATE_SHORT = {
    "Preserved": "Pre",
    "Attenuated": "Att",
    "Reversed": "Rev",
    "Amplified": "Amp",
}

# Classification thresholds for aligned AF aCAI versus reference SR aCAI.
EPS_REVERSED = 0.02
EPS_FLOOR = 0.02
RATIO_ATTENUATED = 0.67
RATIO_AMPLIFIED = 1.33

# Fig1F compartment palette (v3 baseline attribution).
COMPARTMENTS = [
    "atrialCM",
    "fibroblast",
    "endothelial",
    "lymphatic_endothelial",
    "macrophage",
    "lymphocyte",
    "mast_cell",
    "B_cell",
    "pericyte",
    "vascular_smooth_muscle",
    "adipocyte",
    "epicardial",
    "neuronal",
]
COMPARTMENT_PALETTE = [
    "#4E79A7",
    "#A0CBE8",
    "#F28E2B",
    "#FFBE7D",
    "#59A14F",
    "#8CD17D",
    "#B6992D",
    "#F1CE63",
    "#499894",
    "#86BCB6",
    "#E15759",
    "#FF9D9A",
    "#79706E",
]
COMPARTMENT_COLOR_MAP = dict(zip(COMPARTMENTS, COMPARTMENT_PALETTE))
COMPARTMENT_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "lymphatic_endothelial": "L-Endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "mast_cell": "Mast",
    "B_cell": "B",
    "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC",
    "adipocyte": "Adipo",
    "epicardial": "Epi",
    "neuronal": "Neuro",
}
CARRIER_MIXED_LABEL = "Mixed/Unresolved"
CARRIER_MIXED_COLOR = "#BDBDBD"
CARRIER_MIN_Z = 1.0
CARRIER_MIN_GAP = 0.25


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
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 8.5,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Draw rewiring state map.")
    p.add_argument("--acai-tsv", default=str(DEFAULT_ACAI_TSV), help="Input donor-level aCAI table.")
    p.add_argument("--fam-tsv", default=str(DEFAULT_FAM_TSV), help="Input pathway-family mapping table.")
    p.add_argument("--out-tsv", default=str(DEFAULT_OUT_TSV), help="Output state table TSV.")
    p.add_argument("--out-pdf", default=str(DEFAULT_FIG_OUT), help="Output PDF path.")
    return p.parse_args()


def classify_state(af_aligned: float, ref_aligned: float) -> str:
    if af_aligned <= -EPS_REVERSED:
        return "Reversed"
    ratio = af_aligned / max(ref_aligned, EPS_FLOOR)
    if ratio < RATIO_ATTENUATED:
        return "Attenuated"
    if ratio <= RATIO_AMPLIFIED:
        return "Preserved"
    return "Amplified"


def load_chamber_carriers(attr_tsv: Path) -> pd.DataFrame:
    """Load family-level chamber-resolved carriers from the HCA Fig1F output.

    Expects columns family_id, LA_carrier, RA_carrier, consensus_class.
    Returns a DataFrame with one row per family, missing families are added
    later as Mixed/Unresolved.
    """
    if not attr_tsv.exists():
        return pd.DataFrame(columns=["family_id", "LA_carrier", "RA_carrier", "consensus_class"])
    attr = pd.read_csv(attr_tsv, sep="\t")
    needed = {"family_id", "LA_carrier", "RA_carrier", "consensus_class"}
    missing = needed - set(attr.columns)
    if missing:
        raise ValueError(
            f"{attr_tsv.name} is missing required columns {missing}; "
            "expected the HCA chamber-resolved family carrier table."
        )
    out = attr[["family_id", "LA_carrier", "RA_carrier", "consensus_class"]].copy()
    out["family_id"] = out["family_id"].astype(int)
    return out


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
    fam = pd.read_csv(fam_tsv, sep="\t")[["pathway", "family_id"]].rename(columns={"pathway": "program_id"})
    long_df = long_df.merge(fam, on="program_id", how="left")
    long_df = long_df.dropna(subset=["family_id"]).copy()
    long_df["family_id"] = long_df["family_id"].astype(int)

    cohorts = ["GSE128188", "GSE79768", "GSE115574"]
    dir_sign = {"LA-biased": 1.0, "RA-biased": -1.0}

    # Representative pathway label per family (largest |mean aCAI|).
    rep = (
        long_df.groupby(["family_id", "program_id", "direction"], as_index=False)["aCAI"]
        .mean()
        .assign(abs_mean=lambda d: d["aCAI"].abs())
        .sort_values(["family_id", "abs_mean"], ascending=[True, False])
        .groupby("family_id", as_index=False)
        .first()[["family_id", "program_id", "direction"]]
        .rename(columns={"program_id": "rep_pathway"})
    )
    rep["short"] = rep["rep_pathway"].map(lambda s: sanitize_label(short_pathway_label(s, 24)))

    # Pooled SR aligned baseline per family (fallback reference for weak/negative cohort SR).
    pooled_sr = (
        long_df[long_df["rhythm"] == "SR"]
        .groupby("family_id", as_index=False)["aCAI"]
        .mean()
        .rename(columns={"aCAI": "sr_pooled"})
    )
    pooled_sr = pooled_sr.merge(rep[["family_id", "direction"]], on="family_id", how="left")
    pooled_sr["sr_pooled_aligned"] = pooled_sr["sr_pooled"] * pooled_sr["direction"].map(dir_sign)
    pooled_ref = dict(zip(pooled_sr["family_id"], pooled_sr["sr_pooled_aligned"]))

    rows = []
    for (fid, cohort), g in long_df.groupby(["family_id", "cohort"]):
        if cohort not in cohorts:
            continue
        donor_scores = g.groupby(["donor_id", "rhythm"], as_index=False)["aCAI"].mean()
        af = donor_scores[donor_scores["rhythm"] == "AF"]["aCAI"]
        sr = donor_scores[donor_scores["rhythm"] == "SR"]["aCAI"]
        if len(af) == 0 or len(sr) == 0:
            continue

        direction = g["direction"].iloc[0]
        sign = dir_sign[direction]
        sr_mean = float(sr.mean())
        af_mean = float(af.mean())
        sr_aligned = sign * sr_mean
        af_aligned = sign * af_mean

        ref_aligned = sr_aligned
        ref_source = "cohort_SR"
        if not np.isfinite(ref_aligned) or ref_aligned <= EPS_FLOOR:
            ref_aligned = pooled_ref.get(fid, np.nan)
            ref_source = "pooled_SR"
        if not np.isfinite(ref_aligned) or ref_aligned <= EPS_FLOOR:
            ref_aligned = EPS_FLOOR
            ref_source = "floor"

        ratio = af_aligned / ref_aligned
        state = classify_state(af_aligned, ref_aligned)
        rows.append({
            "family_id": fid,
            "cohort": cohort,
            "direction": direction,
            "n_SR": int((donor_scores["rhythm"] == "SR").sum()),
            "n_AF": int((donor_scores["rhythm"] == "AF").sum()),
            "mean_aCAI_SR": round(sr_mean, 4),
            "mean_aCAI_AF": round(af_mean, 4),
            "aligned_SR": round(sr_aligned, 4),
            "aligned_AF": round(af_aligned, 4),
            "reference_aligned_SR": round(float(ref_aligned), 4),
            "reference_source": ref_source,
            "af_over_ref": round(float(ratio), 4),
            "state": state,
        })

    state_df = pd.DataFrame(rows)
    state_df = state_df.merge(rep[["family_id", "short"]], on="family_id", how="left")
    state_df["cohort"] = pd.Categorical(state_df["cohort"], categories=cohorts, ordered=True)
    state_df["dir_rank"] = state_df["direction"].map({"LA-biased": 0, "RA-biased": 1}).fillna(2)
    state_df = state_df.sort_values(["dir_rank", "family_id", "cohort"]).reset_index(drop=True)
    state_df = state_df.drop(columns=["dir_rank"])

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    state_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"wrote {out_tsv.relative_to(ROOT)}  ({len(state_df)} rows)")

    # ---- Draw state map ----------------------------------------------------
    row_meta = (
        state_df[["family_id", "direction", "short"]]
        .drop_duplicates()
        .sort_values(["direction", "family_id"], ascending=[True, True])
        .reset_index(drop=True)
    )
    row_meta["dir_rank"] = row_meta["direction"].map({"LA-biased": 0, "RA-biased": 1})
    row_meta = (
        row_meta.sort_values(["dir_rank", "family_id"]).drop(columns=["dir_rank"]).reset_index(drop=True)
    )

    row_ids = row_meta["family_id"].tolist()
    n_rows = len(row_ids)
    n_cols = len(cohorts)
    row_idx = {fid: i for i, fid in enumerate(row_ids)}
    col_idx = {c: j for j, c in enumerate(cohorts)}

    state_mat = np.empty((n_rows, n_cols), dtype=object)
    ratio_mat = np.full((n_rows, n_cols), np.nan)
    src_mat = np.empty((n_rows, n_cols), dtype=object)
    for _, r in state_df.iterrows():
        i = row_idx[int(r["family_id"])]
        j = col_idx[str(r["cohort"])]
        state_mat[i, j] = r["state"]
        ratio_mat[i, j] = float(r["af_over_ref"])
        src_mat[i, j] = r["reference_source"]

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(190 * MM, 150 * MM))
    # Three-column layout: LEFT (LA/RA group strip) | CENTER (state heatmap)
    # | RIGHT (legend + cohort counts).  Carrier is intentionally NOT shown
    # in Fig 2C; LA/RA carrier display lives in Fig 1F (HCA) and the
    # Integrated Atlas only.
    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[1.2, 6.8, 2.0],
        left=0.07,
        right=0.985,
        top=0.90,
        bottom=0.10,
        wspace=0.05,
    )
    ax_dir = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    ax_leg = fig.add_subplot(gs[0, 2])

    # Direction strip (left): one solid block per direction group with white
    # rotated text; high contrast for print legibility.
    n_la_rows = int((row_meta["direction"] == "LA-biased").sum())
    n_ra_rows = n_rows - n_la_rows
    LA_DEEP = "#9C5550"  # deep desaturated red
    RA_DEEP = "#4E6B89"  # deep desaturated blue
    if n_la_rows > 0:
        ax_dir.add_patch(mpatches.Rectangle((0, 0), 1, n_la_rows, facecolor=LA_DEEP, edgecolor="none"))
        ax_dir.text(0.5, n_la_rows / 2.0, "LA-biased", ha="center", va="center",
                    color="white", fontsize=8.0, fontweight="bold", rotation=90)
    if n_ra_rows > 0:
        ax_dir.add_patch(mpatches.Rectangle((0, n_la_rows), 1, n_ra_rows, facecolor=RA_DEEP, edgecolor="none"))
        ax_dir.text(0.5, n_la_rows + n_ra_rows / 2.0, "RA-biased", ha="center", va="center",
                    color="white", fontsize=8.0, fontweight="bold", rotation=90)
    ax_dir.set_xlim(0, 1)
    ax_dir.set_ylim(n_rows, 0)
    ax_dir.axis("off")

    # Main discrete state grid.
    for i in range(n_rows):
        for j in range(n_cols):
            state = state_mat[i, j]
            fc = STATE_COLORS.get(state, "#EEEEEE")
            ax.add_patch(mpatches.Rectangle((j, i), 1, 1, facecolor=fc, edgecolor="white", linewidth=0.8))
            tag = STATE_SHORT.get(state, "")
            star = "*" if src_mat[i, j] != "cohort_SR" else ""
            ax.text(j + 0.5, i + 0.43, f"{tag}{star}", ha="center", va="center", fontsize=6.4, color="#1F1F1F", fontweight="bold")
            if np.isfinite(ratio_mat[i, j]):
                ax.text(j + 0.5, i + 0.72, f"{ratio_mat[i, j]:.2f}", ha="center", va="center", fontsize=5.3, color="#404040")

    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(["GSE128188", "GSE79768", "GSE115574"], fontsize=7.0, fontweight="bold")
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels([f"F{fid}  {row_meta.loc[i, 'short']}" for i, fid in enumerate(row_ids)], fontsize=6.3)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # LA / RA group boundary line (across main heatmap + LA/RA strip).
    if 0 < n_la_rows < n_rows:
        ax.axhline(n_la_rows, color="#1A1A1A", lw=1.6, zorder=10)
        ax_dir.axhline(n_la_rows, color="#1A1A1A", lw=1.6, zorder=10, clip_on=False)

    # Right legend block: state legend (with thresholds) + cohort counts
    # + star annotation.  No carrier info here -- carrier lives in Fig 1F /
    # Atlas only.
    ax_leg.axis("off")
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(1, 0)

    state_threshold = {
        "Preserved": "AF/ref  0.67 - 1.33",
        "Attenuated": "AF/ref  < 0.67",
        "Reversed": "aligned AF \u2264 -0.02",
        "Amplified": "AF/ref  > 1.33",
    }

    # ---- Upper: Rewiring states with thresholds ----
    ax_leg.text(0.0, 0.02, "Rewiring states", fontsize=7.8, fontweight="bold", ha="left", va="top")
    y = 0.075
    dy_state = 0.06
    for state in ["Preserved", "Attenuated", "Reversed", "Amplified"]:
        ax_leg.add_patch(mpatches.Rectangle((0.00, y), 0.10, 0.040, facecolor=STATE_COLORS[state], edgecolor="none"))
        ax_leg.text(0.14, y + 0.013, state, fontsize=6.8, ha="left", va="center", color="#1A1A1A", fontweight="bold")
        ax_leg.text(0.14, y + 0.030, state_threshold[state], fontsize=5.5, ha="left", va="center", color="#5A5A5A", style="italic")
        y += dy_state

    # ---- Divider ----
    ax_leg.add_patch(mpatches.Rectangle((0.0, 0.36), 0.95, 0.0035, facecolor="#D5D5D5", edgecolor="none"))

    # ---- Middle: State counts by cohort ----
    ax_leg.text(0.0, 0.40, "State counts by cohort", fontsize=7.4, fontweight="bold", ha="left", va="top")
    ctab = pd.crosstab(state_df["cohort"], state_df["state"]).reindex(cohorts).fillna(0).astype(int)
    cohort_short = {"GSE128188": "128", "GSE79768": "797", "GSE115574": "115"}
    y = 0.46
    for cohort in cohorts:
        row = ctab.loc[cohort]
        ax_leg.text(0.0, y, cohort, fontsize=6.4, fontweight="bold", color="#1A1A1A", ha="left", va="center")
        y += 0.030
        ax_leg.text(
            0.04, y,
            f"Pre {row.get('Preserved', 0)}   Att {row.get('Attenuated', 0)}   "
            f"Rev {row.get('Reversed', 0)}   Amp {row.get('Amplified', 0)}",
            fontsize=5.9, color="#333333", ha="left", va="center",
        )
        y += 0.044

    # ---- Divider ----
    ax_leg.add_patch(mpatches.Rectangle((0.0, 0.86), 0.95, 0.0025, facecolor="#E5E5E5", edgecolor="none"))

    # ---- Bottom: star annotation ----
    n_star = int((state_df["reference_source"] != "cohort_SR").sum())
    ax_leg.text(
        0.0, 0.90,
        f"*  uses pooled / floor SR reference\n   ({n_star} / {len(state_df)} cells)",
        fontsize=5.8, ha="left", va="top", color=GREY_MID,
    )

    fig.suptitle(
        "Fig 2C  Rewiring state map (family \u00d7 cohort)",
        x=0.07,
        y=0.965,
        ha="left",
        fontsize=8.8,
        fontweight="bold",
    )
    fig.text(
        0.07,
        0.025,
        "Each cell is one family in one cohort.  State call uses healthy-aligned aCAI ratio "
        "(AF / SR baseline; aligned by healthy direction so LA-biased = + and RA-biased = -).  "
        "Subtext = AF / ref ratio.  Cells marked * fall back to a pooled / floor SR reference "
        "when the cohort-specific SR baseline is too weak.",
        fontsize=5.9,
        color="#444444",
        ha="left",
        va="bottom",
    )

    fig_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_out, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {fig_out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
