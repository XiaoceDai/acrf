#!/usr/bin/env python3
"""Integrated AF rewiring atlas (decision matrix)."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_style import LA_COLOR, RA_COLOR, set_nature_style

ROOT = Path(__file__).resolve().parents[2]

ATTR_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
DRIVER_TSV = ROOT / "results" / "fig2d_family_forest_v3.tsv"
CONS_TSV = ROOT / "results" / "fig2e_replication_lollipop_v3.tsv"

OUT_TSV = ROOT / "results" / "fig_integrated_af_rewiring_atlas.tsv"
OUT_PDF = ROOT / "figures" / "fig_integrated_af_rewiring_atlas.pdf"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"

PRIORITY_COMPONENTS_DOC = (
    "PriorityScore=100*(0.40*E+0.30*R+0.20*C+0.10*H) with "
    "E=min(1,|pooled_g|/0.8), R={high:1,moderate:0.6,low:0.2}, "
    "C=state_concordance_n/3, H=1-min(1,I2). "
    "priority_score/rank are for display ordering and atlas PDF layout only; "
    "not used in Fig 5b/5c Path A supportive Boolean gates."
)

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

STATE_COLORS = {
    "Preserved": "#A9C9B5",
    "Attenuated": "#E9C27B",
    "Reversed": "#C6A2D7",
    "Amplified": "#E59595",
}
STATE_SHORT = {"Preserved": "Pre", "Attenuated": "Att", "Reversed": "Rev", "Amplified": "Amp"}

DRIVER_COLORS = {"LA-driven": LA_COLOR, "RA-driven": RA_COLOR, "Bilateral": "#9E9E9E"}
# Replication palette: purple-grey single hue ramp.  Decoupled from any other
# semantic axis in the atlas (state / driver / bias / carrier) to prevent
# accidental color collision.  Atlas-only; Fig 2E keeps its own REP_CLASS_COLORS.
REPL_COLORS = {"high": "#5D4E6D", "moderate": "#9B8AA5", "low": "#C5B9CC"}

# Bias / Driver use the same LA/RA brand colors as the rest of the figures,
# but distinguished by alpha: low-alpha = historical baseline (Bias),
# full-alpha = AF-time driver (Driver).
BIAS_ALPHA = 0.55
DRIVER_ALPHA = 1.0


CONSENSUS_COLORS = {
    "SHARED": "#1B5E20",
    "SPLIT": "#B07AA1",
    "LA_ONLY": LA_COLOR,
    "RA_ONLY": RA_COLOR,
    "MIXED": "#9C9C9C",
}
CONSENSUS_SHORT = {
    "SHARED": "S",
    "SPLIT": "X",
    "LA_ONLY": "L",
    "RA_ONLY": "R",
    "MIXED": "M",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build integrated AF rewiring atlas.")
    p.add_argument("--attr-tsv", default=str(ATTR_TSV), help="Fig1F v3 HCA family carrier TSV.")
    p.add_argument("--driver-tsv", default=str(DRIVER_TSV), help="Fig2D family decomposition TSV.")
    p.add_argument("--cons-tsv", default=str(CONS_TSV), help="Fig2E consolidation TSV.")
    p.add_argument("--out-tsv", default=str(OUT_TSV), help="Output integrated TSV.")
    p.add_argument("--out-pdf", default=str(OUT_PDF), help="Output integrated PDF.")
    return p.parse_args()


def resolve_path(path_str: str) -> Path:
    p = Path(path_str)
    return p if p.is_absolute() else ROOT / p


def load_chamber_carriers(attr_tsv: Path) -> pd.DataFrame:
    """Read the HCA chamber-resolved family carrier table emitted by
    scripts/v3/fig1f_v3_HCA_baseline_attribution.py."""
    df = pd.read_csv(attr_tsv, sep="\t")
    keep = ["family_id", "LA_carrier", "RA_carrier", "consensus_class", "LA_top_z", "RA_top_z"]
    have = [c for c in keep if c in df.columns]
    return df[have].copy()


def stars_from_score(score: float) -> str:
    if score >= 75.0:
        return "***"
    if score >= 50.0:
        return "**"
    return "*"


def class_from_score(score: float) -> str:
    if score >= 75.0:
        return "high"
    if score >= 50.0:
        return "medium"
    return "low"


def build_integrated_table(attr_tsv: Path, driver_tsv: Path, cons_tsv: Path) -> pd.DataFrame:
    carrier_df = load_chamber_carriers(attr_tsv)
    driver = pd.read_csv(driver_tsv, sep="\t")
    cons = pd.read_csv(cons_tsv, sep="\t")

    driver_pool = (
        driver[driver["cohort"] == "POOLED"][["family_id", "driver"]]
        .drop_duplicates(subset=["family_id"])
        .copy()
    )
    base = cons.copy()
    base = base.merge(driver_pool, on="family_id", how="left")
    base = base.merge(carrier_df, on="family_id", how="left")
    base["LA_carrier"] = base["LA_carrier"].fillna(CARRIER_MIXED_LABEL)
    base["RA_carrier"] = base["RA_carrier"].fillna(CARRIER_MIXED_LABEL)
    base["consensus_class"] = base["consensus_class"].fillna("MIXED")
    base["driver"] = base["driver"].fillna("Bilateral")

    # Priority score components.
    base["E_strength"] = np.minimum(1.0, np.abs(base["pooled_g"].astype(float)) / 0.8)
    base["R_repl"] = base["replication_class"].astype(str).str.lower().map(
        {"high": 1.0, "moderate": 0.6, "low": 0.2}
    ).fillna(0.2)
    base["C_concord"] = (base["state_concordance_n"].astype(float) / 3.0).clip(0.0, 1.0)
    base["I2_01"] = base["I2"].astype(float).clip(0.0, 1.0)
    base["H_hetero"] = 1.0 - np.minimum(1.0, base["I2_01"])
    base["PriorityScore"] = 100.0 * (
        0.40 * base["E_strength"] +
        0.30 * base["R_repl"] +
        0.20 * base["C_concord"] +
        0.10 * base["H_hetero"]
    )
    base["PriorityScore"] = base["PriorityScore"].round(1)
    base["PriorityStars"] = base["PriorityScore"].map(stars_from_score)
    base["PriorityClass"] = base["PriorityScore"].map(class_from_score)
    base["state_consensus"] = base["dominant_state"].astype(str)
    base["state_concordance"] = base["state_concordance_n"].astype(int).astype(str) + "/3"
    base["ReplicationSupport"] = base["replication_class"].astype(str).str.capitalize()

    base["Family"] = "F" + base["family_id"].astype(int).astype(str)
    base["LA_carrier_short"] = base["LA_carrier"].map(lambda x: "Mix" if str(x) == CARRIER_MIXED_LABEL else COMPARTMENT_SHORT.get(str(x), str(x)[:4]))
    base["RA_carrier_short"] = base["RA_carrier"].map(lambda x: "Mix" if str(x) == CARRIER_MIXED_LABEL else COMPARTMENT_SHORT.get(str(x), str(x)[:4]))
    base["HealthyBias"] = base["direction"].astype(str)
    base["ChamberDriver"] = base["driver"].astype(str)

    base["abs_pooled_g"] = np.abs(base["pooled_g"].astype(float))
    out = base.sort_values(
        ["PriorityScore", "abs_pooled_g", "family_id"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    if LABELS_TSV.is_file():
        labs = pd.read_csv(LABELS_TSV, sep="\t")[["family_id", "family_label"]].drop_duplicates(
            subset=["family_id"], keep="first"
        )
        out = out.merge(labs, on="family_id", how="left")
    else:
        out["family_label"] = ""
    out["priority_score"] = out["PriorityScore"]
    out["rank"] = np.arange(1, len(out) + 1, dtype=int)
    out["components"] = PRIORITY_COMPONENTS_DOC
    cols = [
        "family_id",
        "Family",
        "family_label",
        "priority_score",
        "rank",
        "components",
        "short",
        "LA_carrier",
        "LA_carrier_short",
        "RA_carrier",
        "RA_carrier_short",
        "consensus_class",
        "LA_top_z",
        "RA_top_z",
        "HealthyBias",
        "state_consensus",
        "state_concordance",
        "state_concordance_n",
        "ChamberDriver",
        "ReplicationSupport",
        "replication_class",
        "pooled_g",
        "abs_pooled_g",
        "pooled_ci_lo",
        "pooled_ci_hi",
        "I2_01",
        "PriorityScore",
        "PriorityStars",
        "PriorityClass",
        "E_strength",
        "R_repl",
        "C_concord",
        "H_hetero",
    ]
    return out[cols].copy()


def plot_decision_matrix(df: pd.DataFrame, out_pdf: Path) -> None:
    n = len(df)
    set_nature_style()
    fig_h = max(4.8, 0.33 * n + 1.8)
    fig = plt.figure(figsize=(11.2, fig_h))
    gs = fig.add_gridspec(1, 2, width_ratios=[7.4, 2.2], left=0.27, right=0.985, top=0.90, bottom=0.10, wspace=0.04)
    ax = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])
    ax_leg.axis("off")

    columns = [
        "LA baseline\ncarrier",
        "RA baseline\ncarrier",
        "Baseline carrier\nclass (LA vs RA)",
        "Healthy\nbias",
        "AF rewiring\nconsensus",
        "Chamber\ndriver",
        "Replication\nsupport",
        "Priority",
    ]
    col_widths = [0.9, 0.9, 0.55, 1.0, 1.0, 1.0, 1.0, 1.0]
    x_edges = np.concatenate([[0.0], np.cumsum(col_widths)])
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    x_total = float(x_edges[-1])

    # Row background bands.
    for i in range(n):
        bg = "#FBFBFB" if i % 2 == 0 else "#F3F3F3"
        ax.add_patch(mpatches.Rectangle((0, i), x_total, 1, facecolor=bg, edgecolor="none", zorder=0))

    def _draw_carrier(col_idx: int, i: int, carrier: str, label: str) -> None:
        if carrier == CARRIER_MIXED_LABEL:
            ax.add_patch(
                mpatches.Rectangle(
                    (x_edges[col_idx], i),
                    col_widths[col_idx],
                    1,
                    facecolor=CARRIER_MIXED_COLOR,
                    edgecolor="white",
                    linewidth=0.8,
                    hatch="///",
                )
            )
        else:
            ax.add_patch(
                mpatches.Rectangle(
                    (x_edges[col_idx], i),
                    col_widths[col_idx],
                    1,
                    facecolor=COMPARTMENT_COLOR_MAP.get(carrier, CARRIER_MIXED_COLOR),
                    edgecolor="white",
                    linewidth=0.8,
                )
            )
        ax.text(x_centers[col_idx], i + 0.53, label, ha="center", va="center", fontsize=5.8, color="#1F1F1F")

    # Draw cells.
    for i, r in df.iterrows():
        _draw_carrier(0, i, str(r["LA_carrier"]), str(r["LA_carrier_short"]))
        _draw_carrier(1, i, str(r["RA_carrier"]), str(r["RA_carrier_short"]))

        # Carrier consensus class
        cls = str(r["consensus_class"])
        ax.add_patch(
            mpatches.Rectangle(
                (x_edges[2], i), col_widths[2], 1,
                facecolor=CONSENSUS_COLORS.get(cls, "#9C9C9C"), edgecolor="white", linewidth=0.8,
            )
        )
        ax.text(
            x_centers[2], i + 0.53,
            CONSENSUS_SHORT.get(cls, "?"), ha="center", va="center",
            fontsize=6.0, color="white", fontweight="bold",
        )

        # Healthy bias
        hb = str(r["HealthyBias"])
        hb_col = LA_COLOR if hb == "LA-biased" else RA_COLOR
        hb_tag = "LA" if hb == "LA-biased" else "RA"
        ax.add_patch(
            mpatches.Rectangle((x_edges[3], i), col_widths[3], 1, facecolor=hb_col, edgecolor="white", linewidth=0.8, alpha=BIAS_ALPHA)
        )
        ax.text(x_centers[3], i + 0.53, hb_tag, ha="center", va="center", fontsize=6.2, color="#1F1F1F", fontweight="bold")

        # State consensus
        st = str(r["state_consensus"])
        ax.add_patch(
            mpatches.Rectangle((x_edges[4], i), col_widths[4], 1, facecolor=STATE_COLORS.get(st, "#E8E8E8"), edgecolor="white", linewidth=0.8)
        )
        ax.text(x_centers[4], i + 0.50, STATE_SHORT.get(st, st[:3]), ha="center", va="center", fontsize=6.0, color="#1F1F1F", fontweight="bold")
        ax.text(
            x_edges[5] - 0.07,
            i + 0.20,
            str(r["state_concordance"]),
            ha="right",
            va="center",
            fontsize=5.3,
            color="#2F2F2F",
            fontweight="bold",
        )

        # Chamber driver
        drv = str(r["ChamberDriver"])
        drv_col = DRIVER_COLORS.get(drv, "#B0B0B0")
        drv_tag = {"LA-driven": "LA", "RA-driven": "RA", "Bilateral": "Bi"}.get(drv, "Bi")
        ax.add_patch(
            mpatches.Rectangle((x_edges[5], i), col_widths[5], 1, facecolor=drv_col, edgecolor="white", linewidth=0.8, alpha=DRIVER_ALPHA)
        )
        ax.text(x_centers[5], i + 0.53, drv_tag, ha="center", va="center", fontsize=6.2, color="white", fontweight="bold")

        # Replication support
        repl = str(r["replication_class"]).lower()
        repl_col = REPL_COLORS.get(repl, "#B0B0B0")
        repl_tag = {"high": "H", "moderate": "M", "low": "L"}.get(repl, "L")
        ax.add_patch(
            mpatches.Rectangle((x_edges[6], i), col_widths[6], 1, facecolor=repl_col, edgecolor="white", linewidth=0.8)
        )
        # White text on dark purple-grey, dark text on light purple-grey for legibility.
        repl_text_color = "white" if repl == "high" else "#1F1F1F"
        ax.text(x_centers[6], i + 0.53, repl_tag, ha="center", va="center", fontsize=6.2, color=repl_text_color, fontweight="bold")

        # Priority bubble
        ax.add_patch(mpatches.Rectangle((x_edges[7], i), col_widths[7], 1, facecolor="#F6F2EA", edgecolor="white", linewidth=0.8))
        score = float(r["PriorityScore"])
        bubble_s = 35.0 + 280.0 * (score / 100.0) ** 1.35
        ax.scatter([x_centers[7] - 0.10], [i + 0.50], s=bubble_s, c="#C8A076", edgecolors="#6C5842", linewidths=0.6, zorder=3)
        ax.text(
            x_centers[7] - 0.10,
            i + 0.49,
            f"{score:.1f}",
            ha="center",
            va="center",
            fontsize=5.2,
            color="white",
            fontweight="bold",
            zorder=4,
        )
        ax.text(
            x_edges[8] - 0.12,
            i + 0.52,
            str(r["PriorityStars"]),
            ha="right",
            va="center",
            fontsize=5.5,
            color="#3E2F22",
            fontweight="bold",
            zorder=4,
        )

    # Grid separators.
    for x in x_edges:
        ax.plot([x, x], [0, n], color="#FFFFFF", lw=1.0, zorder=5)
    for y in range(n + 1):
        ax.plot([0, x_total], [y, y], color="#FFFFFF", lw=1.0, zorder=5)

    y_labels = []
    for fam, short, cls in zip(df["Family"], df["short"], df["consensus_class"]):
        mark = " \u2020" if str(cls) in {"MIXED", "SPLIT"} else ""
        y_labels.append(f"{fam}{mark}  {short}")

    ax.set_xlim(0, x_total)
    ax.set_ylim(n, 0)
    ax.set_xticks(x_centers)
    ax.set_xticklabels(columns, fontsize=6.6, fontweight="bold")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=5)
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_yticklabels(y_labels, fontsize=6.2)
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Legend block: vertically stacked sections; only show carrier classes
    # that actually occur in this panel.
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(1, 0)

    def _section_title(text: str, y_in: float) -> float:
        ax_leg.text(0.0, y_in, text, fontsize=6.6, fontweight="bold", color="#1A1A1A")
        return y_in + 0.038

    def _legend_row(text: str, color: str, y_in: float, *, hatch: str | None = None) -> float:
        rect = mpatches.Rectangle(
            (0.00, y_in - 0.018), 0.08, 0.03,
            facecolor=color, edgecolor="#8A8A8A" if hatch else "none",
            linewidth=0.4 if hatch else 0.0,
        )
        if hatch:
            rect.set_hatch(hatch)
        ax_leg.add_patch(rect)
        ax_leg.text(0.11, y_in - 0.003, text, fontsize=5.8, ha="left", va="center", color="#222222")
        return y_in + 0.034

    y = 0.04

    # 1. AF rewiring state
    y = _section_title("AF rewiring state", y)
    for st in ["Amplified", "Preserved", "Attenuated", "Reversed"]:
        y = _legend_row(f"{STATE_SHORT[st]}  =  {st}", STATE_COLORS[st], y)

    # 2. Carrier class — only those present in this panel.
    y += 0.014
    present_classes = list(dict.fromkeys(df["consensus_class"].astype(str).tolist()))
    cls_full_map = {
        "SHARED": "S = LA & RA same carrier",
        "LA_ONLY": "L = only LA passes",
        "RA_ONLY": "R = only RA passes",
        "SPLIT": "X = LA & RA disagree",
        "MIXED": "M = neither side passes",
    }
    cls_order = [c for c in ["SHARED", "LA_ONLY", "RA_ONLY", "SPLIT", "MIXED"] if c in present_classes]
    if cls_order:
        y = _section_title("Carrier class (LA vs RA agree?)", y)
        for cls in cls_order:
            y = _legend_row(cls_full_map[cls], CONSENSUS_COLORS[cls], y)

    # 3. Chamber driver
    y += 0.014
    y = _section_title("Chamber driver", y)
    for name, c in [("LA-driven", LA_COLOR), ("RA-driven", RA_COLOR), ("Bilateral", "#9E9E9E")]:
        y = _legend_row(name, c, y)

    # 4. Replication support
    y += 0.014
    y = _section_title("Replication support", y)
    for name, c in [("High", REPL_COLORS["high"]), ("Moderate", REPL_COLORS["moderate"]), ("Low", REPL_COLORS["low"])]:
        y = _legend_row(name, c, y)

    # 5. Priority score bubbles
    y += 0.018
    y = _section_title("Priority score", y)
    for s, lbl, thresh in [(85, "***", "\u2265 75"), (62, "**", "50\u201374"), (38, "*", "< 50")]:
        size = 35.0 + 280.0 * (s / 100.0) ** 1.35
        ax_leg.scatter([0.04], [y - 0.003], s=size, c="#C8A076", edgecolors="#6C5842", linewidths=0.6)
        ax_leg.text(0.11, y - 0.003, f"{lbl}    {thresh}", ha="left", va="center", fontsize=6.0, color="#222222")
        y += 0.052

    ax_leg.text(
        0.0,
        0.97,
        "PriorityScore = 100 \u00d7 (0.40\u00b7E + 0.30\u00b7R + 0.20\u00b7C + 0.10\u00b7H)\n"
        "E = min(1, |pooled_g| / 0.8);  R = {high:1, mod:0.6, low:0.2}\n"
        "C = state_concordance_n / 3;  H = 1 \u2212 min(1, I\u00b2)\n"
        "Carrier from HCA (paired LA + RA).  \u2020 = SPLIT / MIXED row.",
        fontsize=5.4,
        color="#555555",
        ha="left",
        va="bottom",
    )

    fig.suptitle(
        "Integrated AF rewiring atlas: family-level decision matrix",
        fontsize=9.0,
        fontweight="bold",
        x=0.27,
        y=0.965,
        ha="left",
    )
    fig.text(
        0.27,
        0.03,
        "Rows ordered by PriorityScore (high to low).  State concordance = dominant-state frequency across 3 cohorts.\n"
        "Baseline carrier columns show the healthy executor (Fig 1F HCA, paired LA + RA snRNA from healthy donors).  "
        "It may differ from the AF-time cell type whose activity changes most -- that is validated by Fig 3.",
        fontsize=6.1,
        color="#444444",
        ha="left",
        va="bottom",
    )

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    attr_tsv = resolve_path(args.attr_tsv)
    driver_tsv = resolve_path(args.driver_tsv)
    cons_tsv = resolve_path(args.cons_tsv)
    out_tsv = resolve_path(args.out_tsv)
    out_pdf = resolve_path(args.out_pdf)

    for p in [attr_tsv, driver_tsv, cons_tsv]:
        if not p.exists():
            raise FileNotFoundError(f"Missing input: {p}")

    out_df = build_integrated_table(attr_tsv, driver_tsv, cons_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"wrote {out_tsv.relative_to(ROOT)} ({len(out_df)} rows)")

    plot_decision_matrix(out_df, out_pdf)
    print(f"wrote {out_pdf.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
