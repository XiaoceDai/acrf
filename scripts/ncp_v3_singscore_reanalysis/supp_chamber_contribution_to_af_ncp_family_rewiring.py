#!/usr/bin/env python3
"""Supplementary figure: Chamber contribution to AF-associated NCP family rewiring."""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fig2b_family_donor_heatmap as fig2b


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"

OUT_SOURCE = RESULTS / "supp_chamber_contribution_to_af_ncp_family_rewiring_source.tsv"
OUT_QC = RESULTS / "supp_chamber_contribution_to_af_ncp_family_rewiring_qc.md"
OUT_BASE = FIGURES / "supp_chamber_contribution_to_af_ncp_family_rewiring"

COHORT_ORDER = ["GSE128188", "GSE79768", "GSE115574"]
EXPECTED_FAMILY_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16]

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"
BILATERAL_COLOR = "#A6A6A6"
GREY_LIGHT = "#E6E6E6"
GREY_MID = "#9C9C9C"
GREY_DARK = "#444444"
FG = "#202020"

RATIO_THRESH = 1.25
BUBBLE_SCALE = 0.55
SQRT2 = np.sqrt(2.0)
MM = 1.0 / 25.4


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify_chamber_source(delta_la: float, delta_ra: float, *, ratio_thresh: float = RATIO_THRESH) -> str:
    abs_la = abs(delta_la)
    abs_ra = abs(delta_ra)
    if abs_la >= abs_ra * ratio_thresh:
        return "LA-dominant"
    if abs_ra >= abs_la * ratio_thresh:
        return "RA-dominant"
    return "Bilateral"


def classify_quadrant(delta_la: float, delta_ra: float) -> str:
    if delta_la >= 0 and delta_ra >= 0:
        return "Q1 (+,+)"
    if delta_la < 0 <= delta_ra:
        return "Q2 (-,+)"
    if delta_la < 0 and delta_ra < 0:
        return "Q3 (-,-)"
    return "Q4 (+,-)"


def compute_chamber_contribution(
    family_acairel: pd.DataFrame,
    families: pd.DataFrame,
) -> pd.DataFrame:
    """Compute LA and RA chamber contributions per family and cohort."""
    family_acairel = family_acairel.loc[family_acairel["dataset"].isin(COHORT_ORDER)].copy()
    family_acairel["donor_id"] = family_acairel["donor_id"].astype(str)

    rows: list[dict[str, Any]] = []
    for family_id in EXPECTED_FAMILY_ORDER:
        family_group = family_acairel.loc[family_acairel["family_id"] == family_id]
        family_meta = families.loc[families["family_id"] == family_id].iloc[0]
        family_label = str(family_meta["family_label"])
        canonical_direction = str(family_meta["canonical_direction"])
        canonical_order = int(family_meta["canonical_family_display_order"])

        for cohort in COHORT_ORDER:
            cohort_group = family_group.loc[family_group["dataset"] == cohort]
            af = cohort_group.loc[cohort_group["rhythm"] == "AF"]
            sr = cohort_group.loc[cohort_group["rhythm"] == "SR"]

            delta_la = float(af["family_aligned_score_LA_rel"].mean() - sr["family_aligned_score_LA_rel"].mean())
            delta_ra = float(af["family_aligned_score_RA_rel"].mean() - sr["family_aligned_score_RA_rel"].mean())
            delta_acai = float(af["family_aCAIrel"].mean() - sr["family_aCAIrel"].mean())

            rows.append(
                {
                    "family_id": int(family_id),
                    "family_label": family_label,
                    "canonical_direction": canonical_direction,
                    "canonical_order": canonical_order,
                    "cohort": cohort,
                    "n_AF": len(af),
                    "n_SR": len(sr),
                    "delta_LA_rel": delta_la,
                    "delta_RA_rel": delta_ra,
                    "delta_aCAIrel": delta_acai,
                    "abs_delta_aCAIrel": abs(delta_acai),
                    "chamber_source": classify_chamber_source(delta_la, delta_ra),
                    "quadrant": classify_quadrant(delta_la, delta_ra),
                }
            )

    source = pd.DataFrame(rows)
    require(
        source["cohort"].isin(COHORT_ORDER).all(),
        "Unexpected cohort values in chamber contribution source",
    )
    require(
        source.groupby("family_id").size().eq(3).all(),
        "Each family must have three cohort contributions",
    )
    return source


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 8,
            "axes.titlesize": 8,
            "axes.labelsize": 7.2,
            "xtick.labelsize": 6.0,
            "ytick.labelsize": 6.0,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "lines.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def render(source: pd.DataFrame) -> plt.Figure:
    set_style()

    plot_df = source.copy()
    max_abs = float(
        np.nanmax(np.abs(np.r_[plot_df["delta_LA_rel"].values, plot_df["delta_RA_rel"].values]))
    )
    lim = max(0.08, max_abs * 1.22)
    label_dx = lim * 0.014
    label_dy = lim * 0.010

    fig = plt.figure(figsize=(188 * MM, 120 * MM))
    gs = fig.add_gridspec(
        1, 4, width_ratios=[1.0, 1.0, 1.0, 0.90], left=0.07, right=0.98, top=0.90, bottom=0.13, wspace=0.14
    )
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax_leg = fig.add_subplot(gs[0, 3])

    quadrant_bg = {
        "Q1": "#FAFAFA",
        "Q2": "#FFF4F0",
        "Q3": "#F7FAFD",
        "Q4": "#F2F8F4",
    }
    quadrant_labels = {
        "Q1": "LA up\nRA up",
        "Q2": "LA down\nRA up",
        "Q3": "LA down\nRA down",
        "Q4": "LA up\nRA down",
    }
    quadrant_pos = {"Q1": (0.03, 0.97), "Q2": (0.03, 0.03), "Q3": (0.60, 0.03), "Q4": (0.60, 0.97)}

    driver_colors = {
        "LA-dominant": LA_COLOR,
        "RA-dominant": RA_COLOR,
        "Bilateral": BILATERAL_COLOR,
    }
    direction_markers = {"LA": "o", "RA": "s"}

    for ax, cohort in zip(axes, COHORT_ORDER):
        sub = plot_df[plot_df["cohort"] == cohort].copy()

        ax.add_patch(mpatches.Rectangle((-lim, 0), lim, lim, facecolor=quadrant_bg["Q2"], edgecolor="none", zorder=0))
        ax.add_patch(mpatches.Rectangle((0, 0), lim, lim, facecolor=quadrant_bg["Q1"], edgecolor="none", zorder=0))
        ax.add_patch(mpatches.Rectangle((-lim, -lim), lim, lim, facecolor=quadrant_bg["Q3"], edgecolor="none", zorder=0))
        ax.add_patch(mpatches.Rectangle((0, -lim), lim, lim, facecolor=quadrant_bg["Q4"], edgecolor="none", zorder=0))

        ax.axhline(0, color=GREY_MID, lw=0.7, zorder=1)
        ax.axvline(0, color=GREY_MID, lw=0.7, zorder=1)
        ax.plot([-lim, lim], [-lim, lim], color=GREY_LIGHT, lw=0.6, ls="--", zorder=1)
        ax.plot([-lim, lim], [lim, -lim], color=GREY_LIGHT, lw=0.6, ls="--", zorder=1)

        for q, (xp, yp) in quadrant_pos.items():
            ax.text(
                xp,
                yp,
                quadrant_labels[q],
                transform=ax.transAxes,
                ha="left",
                va="top" if yp > 0.5 else "bottom",
                fontsize=5.0,
                color=GREY_DARK,
                linespacing=1.0,
            )

        for _, row in sub.iterrows():
            family_id = int(row["family_id"])
            x = float(row["delta_LA_rel"])
            y = float(row["delta_RA_rel"])
            size = (35 + 260 * (abs(float(row["delta_aCAIrel"])) / (max_abs + 1e-9))) * BUBBLE_SCALE
            face = driver_colors.get(row["chamber_source"], BILATERAL_COLOR)
            marker = direction_markers.get(row["canonical_direction"][:2], "o")
            ax.scatter(x, y, s=size, marker=marker, facecolor=face, edgecolor="black", linewidth=0.35, alpha=0.92, zorder=3)
            ax.text(x + label_dx, y + label_dy, f"F{family_id}", fontsize=5.3, color="#222222", ha="left", va="center", zorder=4)

        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{cohort}  (n = {len(sub)})", fontsize=7.2, fontweight="bold", pad=4)
        ax.tick_params(axis="both", labelsize=6.0, length=2.5)

    axes[0].set_ylabel(r"$\Delta\mathrm{aCAI}_{RA\ rel}$  (AF - SR)", fontsize=7.2, color=RA_COLOR)
    for ax in axes:
        ax.set_xlabel(r"$\Delta\mathrm{aCAI}_{LA\ rel}$  (AF - SR)", fontsize=7.2, color=LA_COLOR)

    # Legend / summary panel
    ax_leg.axis("off")
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(1.22, 0)
    ax_leg.text(0.0, 0.02, "Chamber source", fontsize=7.2, fontweight="bold", ha="left", va="top", color=FG)

    y = 0.09
    for source in ["LA-dominant", "RA-dominant", "Bilateral"]:
        ax_leg.scatter([0.04], [y - 0.005], s=40, marker="o", facecolor=driver_colors[source], edgecolor="#555555", linewidth=0.45)
        ax_leg.text(0.11, y - 0.005, source, fontsize=6.3, ha="left", va="center", color=GREY_DARK)
        y += 0.06

    ax_leg.plot([0.0, 0.88], [y - 0.01, y - 0.01], color=GREY_LIGHT, lw=0.6)
    y += 0.01
    ax_leg.text(0.0, y, r"Bubble size = $|\Delta aCAI_{rel}|$", fontsize=6.8, fontweight="bold", ha="left", va="top", color=FG)
    y += 0.05
    for frac in [0.25, 0.60, 1.00]:
        size = (35 + 260 * frac) * BUBBLE_SCALE
        val = frac * max_abs
        ax_leg.scatter([0.04], [y], s=size, marker="o", facecolor="#F2F2F2", edgecolor="#7F7F7F", linewidth=0.5)
        ax_leg.text(0.11, y, f"{val:.2f}", fontsize=6.0, ha="left", va="center", color=GREY_DARK)
        y += 0.05

    ax_leg.plot([0.0, 0.88], [y - 0.01, y - 0.01], color=GREY_LIGHT, lw=0.6)
    ax_leg.text(0.0, y, "Baseline direction marker", fontsize=6.8, fontweight="bold", ha="left", va="top", color=FG)
    y += 0.05
    ax_leg.scatter([0.04], [y], s=32, marker="o", facecolor="white", edgecolor="black", linewidth=0.5)
    ax_leg.text(0.11, y, "LA-biased family", fontsize=6.2, ha="left", va="center", color=GREY_DARK)
    y += 0.05
    ax_leg.scatter([0.04], [y], s=32, marker="s", facecolor="white", edgecolor="black", linewidth=0.5)
    ax_leg.text(0.11, y, "RA-biased family", fontsize=6.2, ha="left", va="center", color=GREY_DARK)

    y += 0.05
    ax_leg.plot([0.0, 0.88], [y - 0.01, y - 0.01], color=GREY_LIGHT, lw=0.6)
    ax_leg.text(0.0, y, "Diagonal guides", fontsize=6.8, fontweight="bold", ha="left", va="top", color=FG)
    y += 0.04
    ax_leg.text(0.0, y, "y = x  : symmetric co-shift", fontsize=5.8, ha="left", va="top", color=GREY_DARK)
    y += 0.035
    ax_leg.text(0.0, y, "y = -x : antagonistic shift", fontsize=5.8, ha="left", va="top", color=GREY_DARK)

    y += 0.05
    ax_leg.text(0.0, y, "Source counts by cohort", fontsize=6.8, fontweight="bold", ha="left", va="top", color=FG)
    y += 0.04
    ctab = pd.crosstab(plot_df["cohort"], plot_df["chamber_source"]).reindex(COHORT_ORDER).fillna(0).astype(int)
    for cohort in COHORT_ORDER:
        row = ctab.loc[cohort]
        txt = f"{cohort}: LA {row.get('LA-dominant', 0)} | RA {row.get('RA-dominant', 0)} | Bi {row.get('Bilateral', 0)}"
        ax_leg.text(0.0, y, txt, fontsize=5.8, ha="left", va="top", color=GREY_DARK)
        y += 0.055

    fig.suptitle(
        "Chamber contribution to AF-associated NCP family rewiring",
        x=0.07,
        y=0.965,
        ha="left",
        fontsize=8.6,
        fontweight="bold",
    )
    fig.text(
        0.07,
        0.03,
        "Each point is one family in one cohort. "
        "x = ΔaCAIrel_LA, y = ΔaCAIrel_RA (AF - SR). "
        "Point colour denotes the dominant chamber contribution "
        f"(|LA| ≥ {RATIO_THRESH:.2f}|RA|: LA-dominant; |RA| ≥ {RATIO_THRESH:.2f}|LA|: RA-dominant; else bilateral). "
        "Dashed diagonals are y = x (co-shift) and y = -x (antagonistic shift).",
        fontsize=5.8,
        color=GREY_DARK,
        ha="left",
        va="bottom",
    )

    return fig


def save_outputs(fig: plt.Figure) -> dict[str, Path]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    outputs = {
        "pdf": OUT_BASE.with_suffix(".pdf"),
        "svg": OUT_BASE.with_suffix(".svg"),
        "png": OUT_BASE.with_suffix(".png"),
    }
    for key, path in outputs.items():
        fig.savefig(path, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    return outputs


def write_qc(
    manifest: pd.DataFrame,
    paths: dict[str, Path],
    source: pd.DataFrame,
    outputs: dict[str, Path],
) -> None:
    input_lines = [
        f"- `{path.relative_to(ROOT)}`: `{sha256(path)}`"
        for path in sorted(paths.values(), key=lambda item: str(item))
    ]
    output_lines = [
        f"- `{path.relative_to(ROOT)}`: `{sha256(path)}`"
        for path in [OUT_SOURCE, *outputs.values()]
    ]

    ctab = pd.crosstab(source["cohort"], source["chamber_source"]).reindex(COHORT_ORDER).fillna(0).astype(int)
    count_lines = []
    for cohort in COHORT_ORDER:
        row = ctab.loc[cohort]
        count_lines.append(
            f"- {cohort}: LA-dominant={row.get('LA-dominant', 0)}, "
            f"RA-dominant={row.get('RA-dominant', 0)}, Bilateral={row.get('Bilateral', 0)}"
        )

    text = "\n".join(
        [
            "# Supplementary figure: Chamber contribution to AF-associated NCP family rewiring",
            "",
            "Status: **PASS**",
            "",
            "## Scope",
            "- This figure migrates the corrected chamber-decomposition analysis to the supplement.",
            "- It is not included in the main Figure 2 consistency classification.",
            "- It is not labelled as a mechanism-direction validation panel.",
            "",
            "## Inputs",
            "- Donor-level family_aCAIrel and aligned chamber-relative scores from the corrected Figure 2 manifest.",
            "- No legacy effect or consistency table is read as a plotting input.",
            "",
            "## Computed columns",
            "- delta_LA_rel = mean(family_aligned_score_LA_rel | AF) - mean(family_aligned_score_LA_rel | SR).",
            "- delta_RA_rel = mean(family_aligned_score_RA_rel | AF) - mean(family_aligned_score_RA_rel | SR).",
            "- chamber_source: LA-dominant, RA-dominant, or Bilateral (ratio threshold = 1.25).",
            "",
            "## Chamber-source counts by cohort",
            *count_lines,
            "",
            "## Input SHA256",
            *input_lines,
            "",
            "## Output SHA256",
            *output_lines,
            "",
        ]
    )
    OUT_QC.write_text(text, encoding="utf-8")


def main() -> None:
    manifest, paths = fig2b.load_manifest()
    donor_source, families, donor_stats = fig2b.build_source(paths)
    require(
        float(donor_stats["max_verification_error"]) <= 1e-12,
        "Donor-family validation failed upstream",
    )

    family_acairel = pd.read_csv(paths["family_acairel"], sep="\t")
    require(
        {"family_aligned_score_LA_rel", "family_aligned_score_RA_rel"}.issubset(family_acairel.columns),
        "family_acairel table missing chamber-relative score columns",
    )

    source = compute_chamber_contribution(family_acairel, families)
    source = source.sort_values(["canonical_order", "cohort"])

    RESULTS.mkdir(parents=True, exist_ok=True)
    source.to_csv(OUT_SOURCE, sep="\t", index=False, float_format="%.15g")

    fig = render(source)
    outputs = save_outputs(fig)
    write_qc(manifest, paths, source, outputs)

    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)} ({len(source)} rows)")
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
