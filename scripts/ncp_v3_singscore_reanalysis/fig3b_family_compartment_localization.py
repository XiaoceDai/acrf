#!/usr/bin/env python3
"""Figure 3B: LA family × compartment localization (locked corrected lineage)."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"

DELTA_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
CLS_TSV = ROOT / "results" / "fig3_v3_carrier_fidelity_class.tsv"
CANONICAL_TSV = ROOT / "results" / "ncp_v3_canonical_display_order.tsv"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"
LOCK_NOTE = ROOT / "results" / "fig3_publication_authority_lock.md"

OUT_BASE = FIGURES / "fig3b_family_compartment_localization"
OUT_PREVIEW = FIGURES / "fig3b_family_compartment_localization_preview.png"
OUT_SOURCE = RESULTS / "fig3b_family_compartment_localization_source.tsv"
OUT_SUMMARY = RESULTS / "fig3b_family_compartment_localization_summary.tsv"
OUT_SCALE_QC = RESULTS / "fig3b_family_compartment_localization_scale_qc.tsv"
OUT_QC = RESULTS / "fig3b_family_compartment_localization_qc.md"
APPROVED_SOURCE_SHA256 = "dcc7bf8e69f36790ee99efec8a0cd54ec6e77edf639492a6407366c6a5135ee5"
APPROVED_SCALE_M = 0.0481313711

EXPECTED_LA_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14]
EXPECTED_RA_ORDER = [9, 10, 11, 15, 16]
EXPECTED_FAM_ORDER = EXPECTED_LA_ORDER + EXPECTED_RA_ORDER

COMPARTMENTS = [
    "atrialCM",
    "fibroblast",
    "endothelial",
    "macrophage",
    "lymphocyte",
    "pericyte",
    "vascular_smooth_muscle",
    "adipocyte",
    "neuronal",
]
COMP_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "pericyte": "Peri",
    "vascular_smooth_muscle": "VSMC",
    "adipocyte": "Adipo",
    "neuronal": "Neuro",
}
CARRIER_TO_COMP = {
    "atrialCM": "atrialCM",
    "fibroblast": "fibroblast",
    "endothelial": "endothelial",
    "macrophage": "macrophage",
    "lymphocyte": "lymphocyte",
    "pericyte": "pericyte",
    "vascular_smooth_muscle": "vascular_smooth_muscle",
    "adipocyte": "adipocyte",
    "neuronal": "neuronal",
}

CLASS_ORDER = ["Fidelity", "Relay", "Extension", "Below_detection"]
CLASS_DISPLAY = {
    "Fidelity": "Fidelity",
    "Relay": "Relay",
    "Extension": "Extension",
    "Below_detection": "Below",
}
CLASS_COLORS = {
    "Fidelity": "#8FAF7E",
    "Relay": "#A898B8",
    "Extension": "#7FA8AA",
    "Below_detection": "#CECECE",
}

VIOLET = "#7A6B8F"
OCHRE = "#B8956A"
LA_COLOR = "#9C5550"
RA_COLOR = "#4E6B89"
FG_COLOR = "#1A1A1A"
MUTED_COLOR = "#595959"
FAINT_COLOR = "#7A7A7A"
WHITE = "#FFFFFF"
HATCH_BG = "#E8E8E8"

# Exact short publication labels shared with Figures 1D and 2.
# Must not be truncated or altered.
PUBLICATION_LABELS: dict[int, str] = {
    1: "Translation / RQC",
    2: "Mito Ca²⁺ transport",
    3: "Adipo–lipid axis",
    4: "Mito proteostasis / FA–BCAA",
    5: "tRNA aminoacylation",
    6: "Mito translation / import",
    7: "OXPHOS / ETC",
    8: "Complex I biogenesis",
    12: "rRNA processing",
    13: "Amino-acid metabolism",
    14: "Translation initiation",
    9: "HS-GAG / bone / GPCR",
    10: "GAG disease / TGF-β",
    11: "BMP signalling",
    15: "NTRK2 / TrkB",
    16: "Keratan sulfate",
}

MM = 1 / 25.4
WIDTH_MM = 178
MARGIN_MM = 3
TITLE_H_MM = 11
LEGEND_H_MM = 24
ROW_H_MM = 4.2
HEADER_H_MM = 8
MIN_FONT_PT = 6.5
SIG_MARKER_SIZE = 2.8
CARRIER_LINEWIDTH = 0.6


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 6.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": WHITE,
            "savefig.facecolor": WHITE,
        }
    )


def short_family_label(full_label: str, max_chars: int = 28) -> str:
    s = str(full_label).strip()
    if len(s) <= max_chars:
        return s
    return s[: max_chars - 1] + "…"


def load_tables() -> tuple[pd.DataFrame, pd.DataFrame, list[int], dict[int, str], dict[int, str], dict[int, str]]:
    require(LOCK_NOTE.is_file(), f"Missing authority lock: {LOCK_NOTE}")
    require("LOCKED" in LOCK_NOTE.read_text(encoding="utf-8"), "Authority lock note not LOCKED")

    delta = pd.read_csv(DELTA_TSV, sep="\t")
    cls = pd.read_csv(CLS_TSV, sep="\t")
    canonical = pd.read_csv(CANONICAL_TSV, sep="\t")
    labels = pd.read_csv(LABELS_TSV, sep="\t")
    hca = pd.read_csv(HCA_TSV, sep="\t")

    fam_order = (
        canonical[["family_display_order", "family_id", "canonical_direction"]]
        .drop_duplicates("family_id")
        .sort_values("family_display_order")["family_id"]
        .astype(int)
        .tolist()
    )
    require(fam_order == EXPECTED_FAM_ORDER, f"Canonical order drifted: {fam_order}")
    missing_labels = [fid for fid in fam_order if fid not in PUBLICATION_LABELS]
    require(not missing_labels, f"Missing publication labels for families: {missing_labels}")

    dir_lookup = dict(
        zip(
            canonical.drop_duplicates("family_id")["family_id"].astype(int),
            canonical.drop_duplicates("family_id")["canonical_direction"].astype(str),
        )
    )
    # Authoritative short labels are fixed by the publication list above.
    label_lookup = {fid: PUBLICATION_LABELS[fid] for fid in fam_order}
    carrier_lookup = dict(zip(hca["family_id"].astype(int), hca["LA_carrier"].astype(str)))
    return delta, cls, fam_order, dir_lookup, label_lookup, carrier_lookup


def assert_locked_results(delta: pd.DataFrame, cls: pd.DataFrame) -> None:
    require(len(delta) == 144, f"Expected 144 family×compartment rows, got {len(delta)}")
    require(set(delta["compartment"]) == set(COMPARTMENTS), "Compartment set drifted")

    counts = cls["fidelity_class"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int)
    require(int(counts["Fidelity"]) == 1, "Fidelity count != 1")
    require(int(counts["Relay"]) == 2, "Relay count != 2")
    require(int(counts["Extension"]) == 1, "Extension count != 1")
    require(int(counts["Below_detection"]) == 12, "Below-detection count != 12")

    f10 = cls.loc[cls["family_id"] == 10].iloc[0]
    require(str(f10["fidelity_class"]) == "Fidelity", "F10 not Fidelity")

    f11 = cls.loc[cls["family_id"] == 11].iloc[0]
    require(str(f11["fidelity_class"]) == "Relay", "F11 not Relay")
    require("adipocyte" in str(f11["sig_compartments"]), "F11 adipocyte not supported")

    f15 = cls.loc[cls["family_id"] == 15].iloc[0]
    require(str(f15["fidelity_class"]) == "Relay", "F15 not Relay")
    require(str(f15["sig_compartments"]) == "endothelial", f"F15 support not endothelial-only: {f15['sig_compartments']}")
    f15_lymph = delta.loc[(delta["family_id"] == 15) & (delta["compartment"] == "lymphocyte")].iloc[0]
    require(float(f15_lymph["p_perm"]) >= 0.05, "F15 lymphocyte incorrectly significant")

    f16 = cls.loc[cls["family_id"] == 16].iloc[0]
    require(str(f16["fidelity_class"]) == "Extension", "F16 not Extension")
    require(set(str(f16["sig_compartments"]).split(",")) == {"atrialCM", "fibroblast"}, "F16 supports drifted")
    f16_ad = delta.loc[(delta["family_id"] == 16) & (delta["compartment"] == "adipocyte")].iloc[0]
    require(float(f16_ad["p_perm"]) > 0.05, "F16 adipocyte should not be marked p<0.05")
    require(abs(float(f16_ad["p_perm"]) - 0.0542) < 0.002, "F16 adipocyte p drifted from locked 0.0542")

    f4 = cls.loc[cls["family_id"] == 4].iloc[0]
    require(str(f4["fidelity_class"]) == "Below_detection", "F4 not Below-detection")
    require(delta.loc[(delta["family_id"] == 4) & (delta["p_perm"] < 0.05)].empty, "F4 has family-level p<0.05")


def build_source(
    delta: pd.DataFrame,
    cls: pd.DataFrame,
    fam_order: list[int],
    dir_lookup: dict[int, str],
    label_lookup: dict[int, str],
    carrier_lookup: dict[int, str],
) -> tuple[pd.DataFrame, float]:
    class_lookup = dict(zip(cls["family_id"].astype(int), cls["fidelity_class"].astype(str)))
    rows = []
    for disp_i, fid in enumerate(fam_order, start=1):
        carrier = carrier_lookup.get(fid, "Mixed/Unresolved")
        for comp in COMPARTMENTS:
            row = delta.loc[(delta["family_id"] == fid) & (delta["compartment"] == comp)].iloc[0]
            p = float(row["p_perm"]) if pd.notna(row["p_perm"]) else float("nan")
            d = float(row["dScore"]) if pd.notna(row["dScore"]) else float("nan")
            included = bool(row["ct_included"])
            is_carrier = (
                carrier != "Mixed/Unresolved"
                and CARRIER_TO_COMP.get(carrier) == comp
            )
            rows.append(
                {
                    "family_display_order": disp_i,
                    "family_id": fid,
                    "family_label": label_lookup.get(fid, ""),
                    "canonical_direction": dir_lookup.get(fid, ""),
                    "compartment": comp,
                    "compartment_short": COMP_SHORT[comp],
                    "dScore": d,
                    "p_perm": p,
                    "ct_included": included,
                    "n_af": int(row["n_af"]),
                    "n_ctrl": int(row["n_ctrl"]),
                    "is_carrier_cell": is_carrier,
                    "is_sig_p_lt_0_05": bool(pd.notna(p) and p < 0.05),
                    "fidelity_class": class_lookup.get(fid, "Below_detection"),
                    "hca_LA_carrier": carrier,
                }
            )
    source = pd.DataFrame(rows)
    finite = source.loc[source["dScore"].notna(), "dScore"].astype(float)
    m_scale = float(finite.abs().max())
    require(m_scale > 0, "Scale M must be positive")
    return source, m_scale


def donor_counts(source: pd.DataFrame) -> dict[str, tuple[int, int]]:
    out: dict[str, tuple[int, int]] = {}
    for comp in COMPARTMENTS:
        sub = source.loc[source["compartment"] == comp].iloc[0]
        out[comp] = (int(sub["n_af"]), int(sub["n_ctrl"]))
    return out


def assert_approved_source(source: pd.DataFrame, m_scale: float) -> None:
    require(source.shape == (144, 15), "Figure 3B source dimensions changed")
    require(abs(m_scale - APPROVED_SCALE_M) <= 1e-9, "Approved global symmetric scale changed")
    require(source["family_id"].drop_duplicates().tolist() == EXPECTED_FAM_ORDER, "Family order changed")
    require(source["compartment"].drop_duplicates().tolist() == COMPARTMENTS, "Compartment order changed")
    require(int((~source["ct_included"]).sum()) == 0, "Approved eligibility calls changed")
    expected_support = {
        (10, "fibroblast"),
        (11, "adipocyte"),
        (15, "endothelial"),
        (16, "atrialCM"),
        (16, "fibroblast"),
    }
    observed_support = set(
        map(
            tuple,
            source.loc[source["is_sig_p_lt_0_05"], ["family_id", "compartment"]]
            .to_numpy(),
        )
    )
    require(observed_support == expected_support, f"Approved support calls changed: {observed_support}")
    f10 = source.loc[(source["family_id"] == 10) & (source["compartment"] == "fibroblast")].iloc[0]
    require(bool(f10["is_sig_p_lt_0_05"]) and bool(f10["is_carrier_cell"]), "F10 carrier/support overlap changed")
    require(source.loc[source["family_id"] == 10, "fidelity_class"].eq("Fidelity").all(), "F10 class changed")
    require(source.loc[source["family_id"] == 11, "fidelity_class"].eq("Relay").all(), "F11 class changed")
    require(source.loc[source["family_id"] == 15, "fidelity_class"].eq("Relay").all(), "F15 class changed")
    require(source.loc[source["family_id"] == 16, "fidelity_class"].eq("Extension").all(), "F16 class changed")
    require(not source.loc[(source["family_id"] == 15) & (source["compartment"] == "lymphocyte"), "is_sig_p_lt_0_05"].iloc[0], "F15 lymphocyte support changed")
    require(not source.loc[(source["family_id"] == 16) & (source["compartment"] == "adipocyte"), "is_sig_p_lt_0_05"].iloc[0], "F16 adipocyte support changed")
    require(not source.loc[source["family_id"] == 4, "is_sig_p_lt_0_05"].any(), "F4 support changed")
    for comp, group in source.groupby("compartment", sort=False):
        require(group["n_af"].nunique() == 1 and group["n_ctrl"].nunique() == 1, f"{comp} donor counts vary by family")


def render(source: pd.DataFrame, m_scale: float) -> plt.Figure:
    set_style()
    n_rows = len(EXPECTED_FAM_ORDER)
    n_cols = len(COMPARTMENTS)
    matrix_h = n_rows * ROW_H_MM
    height_mm = MARGIN_MM + TITLE_H_MM + HEADER_H_MM + matrix_h + LEGEND_H_MM + MARGIN_MM

    fam_order = source["family_id"].drop_duplicates().tolist()
    delta_pivot = source.pivot(index="family_id", columns="compartment", values="dScore")
    delta_pivot = delta_pivot.reindex(index=fam_order, columns=COMPARTMENTS)
    p_pivot = source.pivot(index="family_id", columns="compartment", values="p_perm")
    p_pivot = p_pivot.reindex(index=fam_order, columns=COMPARTMENTS)
    incl_pivot = source.pivot(index="family_id", columns="compartment", values="ct_included")
    incl_pivot = incl_pivot.reindex(index=fam_order, columns=COMPARTMENTS)
    carrier_pivot = source.pivot(index="family_id", columns="compartment", values="is_carrier_cell")
    carrier_pivot = carrier_pivot.reindex(index=fam_order, columns=COMPARTMENTS)
    class_lookup = dict(zip(source["family_id"], source["fidelity_class"]))
    dir_lookup = dict(zip(source["family_id"], source["canonical_direction"]))
    label_lookup = dict(zip(source["family_id"], source["family_label"]))
    n_la = sum(1 for fid in fam_order if dir_lookup[fid] == "LA-biased")
    n_ra = n_rows - n_la

    cmap = LinearSegmentedColormap.from_list(
        "fig3b_sym",
        [(0.0, VIOLET), (0.5, WHITE), (1.0, OCHRE)],
    )
    norm = TwoSlopeNorm(vmin=-m_scale, vcenter=0.0, vmax=m_scale)

    fig = plt.figure(figsize=(WIDTH_MM * MM, height_mm * MM))
    fig.patch.set_facecolor(WHITE)

    top_y = height_mm - MARGIN_MM
    # Title band
    ax_title = fig.add_axes([MARGIN_MM / WIDTH_MM, (top_y - TITLE_H_MM) / height_mm, 1 - 2 * MARGIN_MM / WIDTH_MM, TITLE_H_MM / height_mm])
    ax_title.set_xlim(0, WIDTH_MM - 2 * MARGIN_MM)
    ax_title.set_ylim(0, TITLE_H_MM)
    ax_title.axis("off")
    ax_title.text(0, TITLE_H_MM - 0.5, "B", fontsize=11, fontweight="bold", color=FG_COLOR, va="top")
    ax_title.text(7, TITLE_H_MM - 0.5, "LA compartment localization and carrier-relative classification",
                  fontsize=8.8, fontweight="bold", color=FG_COLOR, va="top")
    ax_title.text(7, TITLE_H_MM - 4.8,
                  "16 frozen NCP families × 9 LA compartments · family-level AF–control localization",
                  fontsize=6.5, color=MUTED_COLOR, va="top")
    ax_title.text(7, TITLE_H_MM - 7.8,
                  "1 Fidelity · 2 Relay · 1 Extension · 12 Below-detection",
                  fontsize=6.5, fontweight="bold", color=FG_COLOR, va="top")

    matrix_bottom = MARGIN_MM + LEGEND_H_MM
    matrix_top = matrix_bottom + matrix_h
    left = MARGIN_MM
    w_dir, w_fid, w_lbl = 3.5, 6.0, 36.0
    w_matrix = 96.0
    w_class = 16.0
    x_dir = left
    x_fid = x_dir + w_dir + 1.0
    x_lbl = x_fid + w_fid + 1.0
    x_mat = x_lbl + w_lbl + 1.5
    x_cls = x_mat + w_matrix + 0.6

    def mm_axes(x_mm: float, y_mm: float, w_mm: float, h_mm: float) -> plt.Axes:
        return fig.add_axes([x_mm / WIDTH_MM, y_mm / height_mm, w_mm / WIDTH_MM, h_mm / height_mm])

    # Direction strip: one continuous LA/RA color band
    ax_dir = mm_axes(x_dir, matrix_bottom, w_dir, matrix_h)
    ax_dir.set_xlim(0, 1)
    ax_dir.set_ylim(n_rows, 0)
    ax_dir.axis("off")
    for i, fid in enumerate(fam_order):
        col = LA_COLOR if dir_lookup[fid] == "LA-biased" else RA_COLOR
        ax_dir.add_patch(mpatches.Rectangle((0, i), 1, 1, facecolor=col, edgecolor=WHITE, linewidth=0.4))
    if n_la > 0:
        ax_dir.text(0.5, n_la / 2, "LA", ha="center", va="center", fontsize=6.5, fontweight="bold", color=WHITE)
    if n_ra > 0:
        ax_dir.text(0.5, n_la + n_ra / 2, "RA", ha="center", va="center", fontsize=6.5, fontweight="bold", color=WHITE)

    # Family ID (neutral)
    ax_fid = mm_axes(x_fid, matrix_bottom, w_fid, matrix_h)
    ax_fid.set_xlim(0, 1)
    ax_fid.set_ylim(n_rows, 0)
    ax_fid.axis("off")
    for i, fid in enumerate(fam_order):
        ax_fid.text(0.95, i + 0.5, f"F{fid}", ha="right", va="center", fontsize=6.5, color=FG_COLOR, fontweight="bold")

    # Labels (neutral)
    ax_lbl = mm_axes(x_lbl, matrix_bottom, w_lbl, matrix_h)
    ax_lbl.set_label("family_labels")
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(n_rows, 0)
    ax_lbl.axis("off")
    for i, fid in enumerate(fam_order):
        ax_lbl.text(0.98, i + 0.5, label_lookup[fid], ha="right", va="center", fontsize=6.5, color=FG_COLOR)

    # Main matrix
    ax_main = mm_axes(x_mat, matrix_bottom, w_matrix, matrix_h)
    ax_main.set_label("family_compartment_matrix")
    mat = delta_pivot.to_numpy(dtype=float)
    pmat = p_pivot.to_numpy(dtype=float)
    incl = incl_pivot.to_numpy(dtype=object)
    carriers = carrier_pivot.to_numpy(dtype=object)
    eligible = np.asarray(incl, dtype=bool)
    masked = np.ma.masked_where(~eligible | ~np.isfinite(mat), mat)
    require(
        np.ma.getmaskarray(masked)[~eligible].all(),
        "An ineligible cell remains visible as a numeric heatmap value",
    )
    im = ax_main.imshow(
        masked, cmap=cmap, norm=norm, aspect="auto",
        extent=[-0.5, n_cols - 0.5, n_rows - 0.5, -0.5],
        interpolation="nearest",
    )
    n_ineligible = int((incl == False).sum())
    hatch_drawn = 0
    for i in range(n_rows):
        for j in range(n_cols):
            included = incl[i, j]
            if included is False or (isinstance(included, (bool, np.bool_)) and not bool(included)):
                ax_main.add_patch(mpatches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    facecolor=HATCH_BG, edgecolor="#AAAAAA", linewidth=0.3, hatch="///", zorder=3,
                ))
                hatch_drawn += 1
                continue
            p = pmat[i, j]
            if np.isfinite(p) and p < 0.05:
                ax_main.plot(j, i, "o", markersize=SIG_MARKER_SIZE, color=FG_COLOR,
                             markeredgecolor=WHITE, markeredgewidth=0.3, zorder=5)
            if carriers[i, j] is True or carriers[i, j] == True:
                ax_main.add_patch(mpatches.Rectangle(
                    (j - 0.46, i - 0.46), 0.92, 0.92,
                    facecolor="none", edgecolor="#555555", linewidth=CARRIER_LINEWIDTH, zorder=4,
                ))
    if n_ineligible > 0:
        require(hatch_drawn == n_ineligible, f"Ineligible hatch mismatch: {hatch_drawn} drawn vs {n_ineligible} ineligible")
    else:
        require(hatch_drawn == 0, f"Unexpected hatch cells when n_ineligible=0: {hatch_drawn}")

    ax_main.set_xlim(-0.5, n_cols - 0.5)
    ax_main.set_ylim(n_rows - 0.5, -0.5)
    ax_main.set_xticks([])
    ax_main.set_yticks([])
    for s in ax_main.spines.values():
        s.set_visible(False)

    # Matrix header: two-level (compartment label, then AF/Ctrl n)
    dcounts = donor_counts(source)
    ax_header = mm_axes(x_mat, matrix_top, w_matrix, HEADER_H_MM)
    ax_header.set_label("compartment_header")
    ax_header.set_xlim(-0.5, n_cols - 0.5)
    ax_header.set_ylim(0, HEADER_H_MM)
    ax_header.axis("off")
    for j, comp in enumerate(COMPARTMENTS):
        ax_header.text(j, 5.8, COMP_SHORT[comp],
                       ha="center", va="baseline", fontsize=6.5, fontweight="bold", color=FG_COLOR)
        n_af, n_ctrl = dcounts[comp]
        ax_header.text(j, 1.6, f"{n_af}/{n_ctrl}",
                       ha="center", va="baseline", fontsize=6.5, color=FAINT_COLOR)

    # Class column heading
    ax_cls_head = mm_axes(x_cls, matrix_top, w_class, HEADER_H_MM)
    ax_cls_head.set_label("class_header")
    ax_cls_head.set_xlim(0, 1)
    ax_cls_head.set_ylim(0, HEADER_H_MM)
    ax_cls_head.axis("off")
    ax_cls_head.text(0.5, 5.8, "Carrier-relative",
                     ha="center", va="baseline", fontsize=6.5, fontweight="semibold",
                     fontfamily="Roboto Condensed", color=FG_COLOR)
    ax_cls_head.text(0.5, 1.6, "class",
                     ha="center", va="baseline", fontsize=6.5, color=FG_COLOR)

    # Class column
    ax_cls = mm_axes(x_cls, matrix_bottom, w_class, matrix_h)
    ax_cls.set_xlim(0, 1)
    ax_cls.set_ylim(n_rows, 0)
    ax_cls.axis("off")
    for i, fid in enumerate(fam_order):
        cname = class_lookup[fid]
        ax_cls.add_patch(mpatches.Rectangle(
            (0.05, i + 0.08), 0.9, 0.84,
            facecolor=CLASS_COLORS.get(cname, "#CCCCCC"),
            edgecolor=WHITE, linewidth=0.5,
        ))
        ax_cls.text(
            0.5, i + 0.5, CLASS_DISPLAY.get(cname, cname),
            ha="center", va="center", fontsize=6.5,
            fontweight="semibold" if cname == "Below_detection" else "bold",
            color=WHITE if cname != "Below_detection" else "#333333",
        )

    # Legend band
    ax_leg = mm_axes(MARGIN_MM, MARGIN_MM, WIDTH_MM - 2 * MARGIN_MM, LEGEND_H_MM)
    ax_leg.set_label("figure_footer")
    leg_w = WIDTH_MM - 2 * MARGIN_MM
    ax_leg.set_xlim(0, leg_w)
    ax_leg.set_ylim(0, LEGEND_H_MM)
    ax_leg.axis("off")
    footer_x = x_mat - MARGIN_MM
    footer_right = footer_x + w_matrix

    # Footer row 1: labels aligned exactly with the 9-compartment matrix
    ax_leg.text(
        footer_x,
        21.8,
        "AF–control family Δsingscore",
        fontsize=6.5,
        fontweight="semibold",
        color=FG_COLOR,
        ha="left",
        va="center",
    )
    ax_leg.text(
        footer_right,
        21.8,
        "Global symmetric scale; no row standardization",
        fontsize=6.5,
        color=FAINT_COLOR,
        ha="right",
        va="center",
    )

    # Footer row 2: colorbar only
    cax = fig.add_axes([
        x_mat / WIDTH_MM,
        (MARGIN_MM + 14.0) / height_mm,
        w_matrix / WIDTH_MM,
        3.0 / height_mm,
    ])
    cax.set_label("figure_colorbar")
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_ticks([-m_scale, 0.0, m_scale])
    cbar.ax.tick_params(labelsize=6.5, length=2, width=0.4, pad=2)
    cbar.outline.set_linewidth(0.4)

    # Footer row 3: legend symbols, separate from colorbar ticks
    sym_y = 7.0
    sym_x = footer_x
    # Carrier outline
    ax_leg.add_patch(mpatches.Rectangle(
        (sym_x + 0.2, sym_y - 0.2), 2.0, 2.0,
        facecolor="none", edgecolor="#555555", linewidth=CARRIER_LINEWIDTH,
    ))
    ax_leg.text(sym_x + 2.8, sym_y + 0.8, "healthy HCA LA carrier",
                fontsize=6.5, va="center", color=MUTED_COLOR)
    # p < 0.05 dot
    ax_leg.plot(sym_x + 32.0 + 1.1, sym_y + 0.8, "o", markersize=SIG_MARKER_SIZE, color=FG_COLOR,
                markeredgecolor=WHITE, markeredgewidth=0.3)
    ax_leg.text(sym_x + 34.5, sym_y + 0.8, "nominal family-level permutation p < 0.05",
                fontsize=6.5, va="center", color=MUTED_COLOR)
    # Hatch (conditional)
    if n_ineligible > 0:
        ax_leg.add_patch(mpatches.Rectangle(
            (sym_x + 82.0, sym_y - 0.2), 2.2, 2.2,
            facecolor=HATCH_BG, edgecolor="#AAAAAA", linewidth=0.4, hatch="///",
        ))
        ax_leg.text(sym_x + 84.8, sym_y + 0.8, "ineligible",
                    fontsize=6.5, va="center", color=MUTED_COLOR)

    # FDR note
    ax_leg.text(sym_x, 1.5,
                "Nominal support is a localization cue, not FDR-controlled discovery.",
                fontsize=6.5, color=FAINT_COLOR, va="bottom")

    return fig


def render_composite_panel(source: pd.DataFrame, m_scale: float) -> plt.Figure:
    """Render exact 178 × 86 mm panel B for the final Figure 3 composite."""
    set_style()
    width_mm, height_mm = 178.0, 86.0
    margin = 2.3
    title_h, header_h, row_h, footer_h = 9.0, 6.2, 3.2, 15.0
    n_rows, n_cols = len(EXPECTED_FAM_ORDER), len(COMPARTMENTS)
    matrix_h = n_rows * row_h
    require(
        abs(2 * margin + title_h + header_h + matrix_h + footer_h - height_mm) <= 1e-9,
        "Panel B vertical geometry changed",
    )
    fam_order = source["family_id"].drop_duplicates().astype(int).tolist()
    require(fam_order == EXPECTED_FAM_ORDER, "Panel B canonical family order changed")
    delta = source.pivot(index="family_id", columns="compartment", values="dScore").reindex(
        index=fam_order, columns=COMPARTMENTS
    )
    pvals = source.pivot(index="family_id", columns="compartment", values="p_perm").reindex(
        index=fam_order, columns=COMPARTMENTS
    )
    included = source.pivot(index="family_id", columns="compartment", values="ct_included").reindex(
        index=fam_order, columns=COMPARTMENTS
    )
    carriers = source.pivot(index="family_id", columns="compartment", values="is_carrier_cell").reindex(
        index=fam_order, columns=COMPARTMENTS
    )
    class_lookup = dict(zip(source["family_id"], source["fidelity_class"]))
    direction_lookup = dict(zip(source["family_id"], source["canonical_direction"]))
    label_lookup = dict(zip(source["family_id"], source["family_label"]))
    counts = pd.Series(class_lookup).value_counts()
    require(
        [int(counts.get(name, 0)) for name in CLASS_ORDER] == [1, 2, 1, 12],
        "Panel B class counts changed",
    )

    fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
    fig.patch.set_facecolor(WHITE)

    def mm_axes(x: float, bottom: float, width: float, height: float) -> plt.Axes:
        return fig.add_axes([x / width_mm, bottom / height_mm, width / width_mm, height / height_mm])

    title_bottom = height_mm - margin - title_h
    title_ax = mm_axes(margin, title_bottom, width_mm - 2 * margin, title_h)
    title_ax.axis("off")
    title_ax.set_xlim(0, width_mm - 2 * margin)
    title_ax.set_ylim(0, title_h)
    title_ax.text(0, 8.5, "B", fontsize=11, fontweight="bold", va="top", color=FG_COLOR)
    title_ax.text(
        7,
        8.5,
        "LA compartment localization and carrier-relative classification",
        fontsize=8.8,
        fontweight="bold",
        va="top",
        color=FG_COLOR,
    )
    title_ax.text(
        7,
        4.5,
        "16 frozen NCP families × 9 LA compartments · family-level AF–control localization",
        fontsize=6.5,
        va="top",
        color=MUTED_COLOR,
    )
    title_ax.text(
        7,
        1.6,
        "1 Fidelity · 2 Relay · 1 Extension · 12 Below-detection",
        fontsize=6.0,
        fontweight="bold",
        va="top",
        color=FG_COLOR,
    )

    matrix_bottom = margin + footer_h
    matrix_top = matrix_bottom + matrix_h
    x_dir, w_dir = margin, 3.5
    x_fid, w_fid = x_dir + w_dir + 1.0, 6.0
    x_lbl, w_lbl = x_fid + w_fid + 1.0, 39.0
    x_mat, w_matrix = x_lbl + w_lbl + 1.5, 103.0
    x_cls, w_class = x_mat + w_matrix + 0.7, 17.0
    require(x_cls + w_class <= width_mm - margin + 1e-9, "Panel B right edge clipped")

    ax_dir = mm_axes(x_dir, matrix_bottom, w_dir, matrix_h)
    ax_dir.set_xlim(0, 1)
    ax_dir.set_ylim(n_rows, 0)
    ax_dir.axis("off")
    n_la = sum(direction_lookup[fid] == "LA-biased" for fid in fam_order)
    for row_index, family_id in enumerate(fam_order):
        color = LA_COLOR if direction_lookup[family_id] == "LA-biased" else RA_COLOR
        ax_dir.add_patch(mpatches.Rectangle((0, row_index), 1, 1, facecolor=color, edgecolor=WHITE, linewidth=0.25))
    ax_dir.text(0.5, n_la / 2, "LA", fontsize=6.0, color=WHITE, fontweight="bold", ha="center", va="center")
    ax_dir.text(0.5, n_la + (n_rows - n_la) / 2, "RA", fontsize=6.0, color=WHITE, fontweight="bold", ha="center", va="center")

    ax_fid = mm_axes(x_fid, matrix_bottom, w_fid, matrix_h)
    ax_fid.set_xlim(0, 1)
    ax_fid.set_ylim(n_rows, 0)
    ax_fid.axis("off")
    ax_lbl = mm_axes(x_lbl, matrix_bottom, w_lbl, matrix_h)
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(n_rows, 0)
    ax_lbl.axis("off")
    for row_index, family_id in enumerate(fam_order):
        ax_fid.text(0.95, row_index + 0.5, f"F{family_id}", fontsize=6.0, fontweight="bold", ha="right", va="center")
        ax_lbl.text(0.98, row_index + 0.5, label_lookup[family_id], fontsize=6.0, ha="right", va="center")

    cmap = LinearSegmentedColormap.from_list("fig3b_composite", [(0, VIOLET), (0.5, WHITE), (1, OCHRE)])
    norm = TwoSlopeNorm(vmin=-m_scale, vcenter=0, vmax=m_scale)
    ax_main = mm_axes(x_mat, matrix_bottom, w_matrix, matrix_h)
    ax_main.set_gid("fig3b-composite-matrix")
    matrix = delta.to_numpy(float)
    p_matrix = pvals.to_numpy(float)
    eligible = included.to_numpy(bool)
    carrier_matrix = carriers.to_numpy(bool)
    image = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    image.set_array(matrix[np.isfinite(matrix)])
    for row_index in range(n_rows):
        for column_index in range(n_cols):
            if not eligible[row_index, column_index]:
                ax_main.add_patch(
                    mpatches.Rectangle(
                        (column_index - 0.5, row_index - 0.5),
                        1,
                        1,
                        facecolor=HATCH_BG,
                        edgecolor="#AAAAAA",
                        linewidth=0.25,
                        hatch="///",
                    )
                )
                continue
            if np.isfinite(matrix[row_index, column_index]):
                ax_main.add_patch(
                    mpatches.Rectangle(
                        (column_index - 0.5, row_index - 0.5),
                        1,
                        1,
                        facecolor=cmap(norm(matrix[row_index, column_index])),
                        edgecolor=WHITE,
                        linewidth=0.12,
                        zorder=1,
                    )
                )
            if np.isfinite(p_matrix[row_index, column_index]) and p_matrix[row_index, column_index] < 0.05:
                ax_main.plot(
                    column_index,
                    row_index,
                    "o",
                    markersize=2.5,
                    color=FG_COLOR,
                    markeredgecolor=WHITE,
                    markeredgewidth=0.25,
                    zorder=5,
                )
            if carrier_matrix[row_index, column_index]:
                ax_main.add_patch(
                    mpatches.Rectangle(
                        (column_index - 0.46, row_index - 0.46),
                        0.92,
                        0.92,
                        facecolor="none",
                        edgecolor="#555555",
                        linewidth=0.55,
                        zorder=4,
                    )
                )
    ax_main.set_xlim(-0.5, n_cols - 0.5)
    ax_main.set_ylim(n_rows - 0.5, -0.5)
    ax_main.set_xticks([])
    ax_main.set_yticks([])
    for spine in ax_main.spines.values():
        spine.set_visible(False)

    donor_count_lookup = donor_counts(source)
    ax_header = mm_axes(x_mat, matrix_top, w_matrix, header_h)
    ax_header.set_xlim(-0.5, n_cols - 0.5)
    ax_header.set_ylim(0, header_h)
    ax_header.axis("off")
    for column_index, compartment in enumerate(COMPARTMENTS):
        ax_header.text(column_index, 4.5, COMP_SHORT[compartment], fontsize=6.2, fontweight="bold", ha="center", va="baseline")
        n_af, n_ctrl = donor_count_lookup[compartment]
        ax_header.text(column_index, 1.2, f"{n_af}/{n_ctrl}", fontsize=6.0, color=FAINT_COLOR, ha="center", va="baseline")

    ax_cls_head = mm_axes(x_cls, matrix_top, w_class, header_h)
    ax_cls_head.set_xlim(0, 1)
    ax_cls_head.set_ylim(0, header_h)
    ax_cls_head.axis("off")
    ax_cls_head.text(0.5, 4.5, "Carrier-relative", fontsize=6.0, fontweight="semibold", ha="center", va="baseline")
    ax_cls_head.text(0.5, 1.2, "class", fontsize=6.0, ha="center", va="baseline")
    ax_cls = mm_axes(x_cls, matrix_bottom, w_class, matrix_h)
    ax_cls.set_xlim(0, 1)
    ax_cls.set_ylim(n_rows, 0)
    ax_cls.axis("off")
    for row_index, family_id in enumerate(fam_order):
        class_name = class_lookup[family_id]
        patch = mpatches.Rectangle(
            (0.04, row_index + 0.08),
            0.92,
            0.84,
            facecolor=CLASS_COLORS[class_name],
            edgecolor=WHITE,
            linewidth=0.3,
        )
        patch.set_gid(f"fig3b-class-{family_id}")
        ax_cls.add_patch(patch)
        ax_cls.text(
            0.5,
            row_index + 0.5,
            CLASS_DISPLAY[class_name],
            fontsize=6.0,
            fontweight="semibold" if class_name == "Below_detection" else "bold",
            color=WHITE if class_name != "Below_detection" else "#333333",
            ha="center",
            va="center",
        )

    footer = mm_axes(x_mat, margin, w_matrix, footer_h)
    footer.set_gid("fig3b-composite-footer")
    footer.set_xlim(0, w_matrix)
    footer.set_ylim(0, footer_h)
    footer.axis("off")
    footer.text(0, 13.7, "AF–control family Δsingscore", fontsize=6.2, fontweight="semibold", va="top")
    footer.text(w_matrix, 13.7, "Global symmetric scale; no row standardization", fontsize=6.0, color=FAINT_COLOR, ha="right", va="top")
    cax = mm_axes(x_mat, margin + 7.6, w_matrix, 2.2)
    cax.set_xlim(-m_scale, m_scale)
    cax.set_ylim(0, 1)
    cax.set_yticks([])
    color_edges = np.linspace(-m_scale, m_scale, 97)
    for start, end in zip(color_edges[:-1], color_edges[1:]):
        cax.add_patch(
            mpatches.Rectangle(
                (start, 0),
                end - start,
                1,
                facecolor=cmap(norm((start + end) / 2)),
                edgecolor="none",
            )
        )
    cax.set_xticks([-m_scale, 0, m_scale])
    cax.tick_params(axis="x", labelsize=6.0, length=1.8, width=0.35, pad=1.2)
    for spine in cax.spines.values():
        spine.set_linewidth(0.35)
    footer.add_patch(mpatches.Rectangle((0.2, 0.9), 1.8, 1.8, facecolor="none", edgecolor="#555555", linewidth=0.55))
    footer.text(2.6, 1.8, "healthy HCA LA carrier", fontsize=6.0, color=MUTED_COLOR, va="center")
    footer.plot(31.5, 1.8, "o", markersize=2.5, color=FG_COLOR, markeredgecolor=WHITE, markeredgewidth=0.25)
    footer.text(33.4, 1.8, "nominal family-level permutation p < 0.05", fontsize=6.0, color=MUTED_COLOR, va="center")
    return fig


def palette_accessibility_qc() -> dict[str, float]:
    def rgb(color: str) -> np.ndarray:
        return np.asarray(mpl.colors.to_rgb(color), dtype=float)

    luminance_weights = np.asarray([0.2126, 0.7152, 0.0722])
    heat = [rgb(VIOLET), rgb(WHITE), rgb(OCHRE)]
    classes = [rgb(CLASS_COLORS[name]) for name in CLASS_ORDER]
    heat_luminance = [float(np.dot(color, luminance_weights)) for color in heat]
    grayscale_sign_distance = abs(heat_luminance[0] - heat_luminance[2])
    dot_luminance = float(np.dot(rgb(FG_COLOR), luminance_weights))
    outline_luminance = float(np.dot(rgb("#555555"), luminance_weights))
    min_dot_contrast = min(abs(dot_luminance - value) for value in heat_luminance)
    min_outline_contrast = min(abs(outline_luminance - value) for value in heat_luminance)
    simulations = {
        "protanopia": np.asarray(
            [[0.56667, 0.43333, 0.0], [0.55833, 0.44167, 0.0], [0.0, 0.24167, 0.75833]]
        ),
        "deuteranopia": np.asarray(
            [[0.625, 0.375, 0.0], [0.700, 0.300, 0.0], [0.0, 0.300, 0.700]]
        ),
        "tritanopia": np.asarray(
            [[0.950, 0.050, 0.0], [0.0, 0.43333, 0.56667], [0.0, 0.475, 0.525]]
        ),
    }
    min_cvd_heat_sign = np.inf
    min_cvd_class = np.inf
    for matrix in simulations.values():
        simulated_heat = [matrix @ color for color in heat]
        simulated_classes = [matrix @ color for color in classes]
        min_cvd_heat_sign = min(
            min_cvd_heat_sign,
            float(np.linalg.norm(simulated_heat[0] - simulated_heat[2])),
        )
        min_cvd_class = min(
            min_cvd_class,
            min(
                float(np.linalg.norm(simulated_classes[left] - simulated_classes[right]))
                for left in range(len(simulated_classes))
                for right in range(left)
            ),
        )
    require(grayscale_sign_distance >= 0.12, "Heatmap signs collapse in grayscale")
    require(min_dot_contrast >= 0.30, "Nominal-support dot loses grayscale contrast")
    require(min_outline_contrast >= 0.08, "Carrier outline loses grayscale contrast")
    require(min_cvd_heat_sign >= 0.20, "Heatmap signs collapse under CVD simulation")
    require(min_cvd_class >= 0.08, "Class fills collapse under CVD simulation")
    return {
        "grayscale_heat_sign": grayscale_sign_distance,
        "grayscale_dot_contrast": min_dot_contrast,
        "grayscale_outline_contrast": min_outline_contrast,
        "cvd_heat_sign": float(min_cvd_heat_sign),
        "cvd_class": float(min_cvd_class),
    }


def validate_figure_layout(fig: plt.Figure, source: pd.DataFrame) -> dict[str, float | int]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [
        text
        for text in fig.findobj(mpl.text.Text)
        if text.get_visible() and text.get_text().strip()
    ]
    require(texts, "Figure 3B contains no visible text")
    min_font = min(float(text.get_fontsize()) for text in texts)
    require(min_font >= MIN_FONT_PT, f"Visible text below {MIN_FONT_PT:.1f} pt: {min_font:.2f}")
    canvas = fig.bbox
    tolerance = 1.0
    outside = []
    for text in texts:
        bbox = text.get_window_extent(renderer=renderer)
        if (
            bbox.x0 < canvas.x0 - tolerance
            or bbox.y0 < canvas.y0 - tolerance
            or bbox.x1 > canvas.x1 + tolerance
            or bbox.y1 > canvas.y1 + tolerance
        ):
            outside.append(text.get_text())
    require(not outside, f"Visible text outside canvas: {outside}")

    axes = {ax.get_label(): ax for ax in fig.axes}
    header = axes["compartment_header"]
    header_labels = [text for text in header.texts if text.get_text() in COMP_SHORT.values()]
    count_values = {f"{n_af}/{n_ctrl}" for n_af, n_ctrl in donor_counts(source).values()}
    header_counts = [text for text in header.texts if text.get_text() in count_values]
    require(len(header_labels) == 9 and len(header_counts) == 9, "Two-level header is incomplete")
    for label, count in zip(header_labels, header_counts):
        require(
            not label.get_window_extent(renderer).overlaps(count.get_window_extent(renderer)),
            f"Header label/count overlap: {label.get_text()}",
        )

    label_axis = axes["family_labels"]
    require(
        all(
            text.get_window_extent(renderer).x0 >= label_axis.bbox.x0 - tolerance
            and text.get_window_extent(renderer).x1 <= label_axis.bbox.x1 + tolerance
            for text in label_axis.texts
        ),
        "A family label exceeds its allotted column",
    )
    global_note = next(
        text for text in texts
        if text.get_text() == "Global symmetric scale; no row standardization"
    )
    global_bbox = global_note.get_window_extent(renderer)
    require(global_bbox.x1 <= canvas.x1 + tolerance, "Global scale note is clipped")

    colorbar_axis = axes["figure_colorbar"]
    tick_boxes = [
        tick.get_window_extent(renderer)
        for tick in colorbar_axis.get_xticklabels()
        if tick.get_visible() and tick.get_text().strip()
    ]
    footer = axes["figure_footer"]
    symbol_text = [
        text for text in footer.texts
        if text.get_text() in {
            "healthy HCA LA carrier",
            "nominal family-level permutation p < 0.05",
            "ineligible",
        }
    ]
    require(
        not any(
            text.get_window_extent(renderer).overlaps(tick_box)
            for text in symbol_text
            for tick_box in tick_boxes
        ),
        "Symbol legend overlaps colorbar tick labels",
    )

    n_ineligible = int((~source["ct_included"]).sum())
    hatch_patches = [
        patch
        for ax in fig.axes
        for patch in ax.patches
        if patch.get_hatch()
    ]
    expected_hatches = n_ineligible + (1 if n_ineligible > 0 else 0)
    require(len(hatch_patches) == expected_hatches, "Ineligible hatch/legend parity failed")
    require(
        ("ineligible" in [text.get_text() for text in symbol_text]) == (n_ineligible > 0),
        "Conditional hatch legend visibility changed",
    )
    return {
        "min_font_pt": min_font,
        "n_visible_text": len(texts),
        "n_bbox_violations": len(outside),
        "n_header_labels": len(header_labels),
        "n_header_counts": len(header_counts),
        "n_hatch_patches": len(hatch_patches),
    }


def write_outputs(source: pd.DataFrame, m_scale: float) -> dict[str, str]:
    counts = source.drop_duplicates("family_id")["fidelity_class"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int)
    n_sig = int(source["is_sig_p_lt_0_05"].sum())
    n_hatch = int((~source["ct_included"]).sum())
    n_outline = int(source["is_carrier_cell"].sum())
    finite = source.loc[source["dScore"].notna(), ["family_id", "compartment", "dScore"]]
    min_row = finite.loc[finite["dScore"].idxmin()]
    max_row = finite.loc[finite["dScore"].idxmax()]

    summary = pd.DataFrame(
        [
            {"metric": "n_families", "value": 16},
            {"metric": "n_compartments", "value": 9},
            {"metric": "order_provenance", "value": "ncp_v3_canonical_display_order family_display_order"},
            {"metric": "n_fidelity", "value": int(counts["Fidelity"])},
            {"metric": "n_relay", "value": int(counts["Relay"])},
            {"metric": "n_extension", "value": int(counts["Extension"])},
            {"metric": "n_below_detection", "value": int(counts["Below_detection"])},
            {"metric": "n_sig_cells_p_lt_0_05", "value": n_sig},
            {"metric": "scale_M_abs_max_dScore", "value": m_scale},
        ]
    )
    scale_qc = pd.DataFrame(
        [
            {"vmin": -m_scale, "vcenter": 0.0, "vmax": m_scale, "row_standardized": False},
            {"metric": "n_hatch", "value": n_hatch},
            {"metric": "n_carrier_outline", "value": n_outline},
            {"metric": "n_sig_dot", "value": n_sig},
            {"metric": "min_cell", "value": f"F{int(min_row.family_id)}×{min_row.compartment}={float(min_row.dScore):.6g}"},
            {"metric": "max_cell", "value": f"F{int(max_row.family_id)}×{max_row.compartment}={float(max_row.dScore):.6g}"},
        ]
    )

    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    source.to_csv(OUT_SOURCE, sep="\t", index=False, float_format="%.10g")
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    fig = render(source, m_scale)
    layout_qc = validate_figure_layout(fig, source)
    palette_qc = palette_accessibility_qc()
    extra_qc = pd.DataFrame(
        [
            {"metric": key, "value": value}
            for key, value in {**layout_qc, **palette_qc}.items()
        ]
    )
    scale_qc = pd.concat([scale_qc, extra_qc], ignore_index=True)
    scale_qc.to_csv(OUT_SCALE_QC, sep="\t", index=False)
    pdf = OUT_BASE.with_suffix(".pdf")
    svg = OUT_BASE.with_suffix(".svg")
    png = OUT_BASE.with_suffix(".png")
    fig.savefig(pdf, bbox_inches=None, pad_inches=0)
    fig.savefig(svg, bbox_inches=None, pad_inches=0)
    fig.savefig(png, dpi=600, bbox_inches=None, pad_inches=0)
    fig.savefig(OUT_PREVIEW, dpi=150, bbox_inches=None, pad_inches=0)
    plt.close(fig)

    output_hashes = {
        str(p.relative_to(ROOT)): sha256(p)
        for p in (pdf, svg, png, OUT_PREVIEW, OUT_SOURCE, OUT_SUMMARY, OUT_SCALE_QC)
    }
    OUT_QC.write_text(
        "\n".join(
            [
                "# Figure 3B family compartment localization QC",
                "",
                "Status: **PASS**",
                "",
                "## Assertions",
                "- canonical order: LA F1–F8,F12–F14 then RA F9,F10,F11,F15,F16",
                "- class counts: 1 Fidelity | 2 Relay | 1 Extension | 12 Below-detection",
                "- F10 Fidelity; F11 Relay adipocyte; F15 Relay endothelial only; F16 Extension fibroblast+atrialCM",
                "- F16 adipocyte p≈0.0542 not marked; F4 Below-detection; BCAA label retained only in family name",
                "- no PriorityScore order; no p<0.10 trends",
                f"- n_ineligible = {n_hatch}; hatch legend {'shown' if n_hatch > 0 else 'omitted'}",
                "",
                f"- scale M = {m_scale:.6g}",
                f"- n_sig p<0.05 = {n_sig}",
                f"- minimum visible font = {float(layout_qc['min_font_pt']):.2f} pt",
                f"- visible-text bbox violations = {int(layout_qc['n_bbox_violations'])}",
                f"- two-level header labels/counts = {int(layout_qc['n_header_labels'])}/{int(layout_qc['n_header_counts'])}",
                f"- grayscale heat-sign distance = {palette_qc['grayscale_heat_sign']:.3f}",
                f"- grayscale dot/outline contrast = {palette_qc['grayscale_dot_contrast']:.3f}/{palette_qc['grayscale_outline_contrast']:.3f}",
                f"- minimum simulated-CVD heat/class distance = {palette_qc['cvd_heat_sign']:.3f}/{palette_qc['cvd_class']:.3f}",
                "",
                "## Outputs",
                *[f"- `{k}` — `{v}`" for k, v in output_hashes.items()],
                "",
            ]
        ),
        encoding="utf-8",
    )
    return output_hashes


def validate_pdf(pdf: Path) -> None:
    pdftotext = shutil.which("pdftotext")
    pdfinfo = shutil.which("pdfinfo")
    require(pdftotext is not None and pdfinfo is not None, "PDF QC tools unavailable")
    info = subprocess.run(
        [pdfinfo, str(pdf)], check=True, capture_output=True, text=True
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read Figure 3B PDF dimensions")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    expected_height_mm = (
        MARGIN_MM + TITLE_H_MM + HEADER_H_MM
        + len(EXPECTED_FAM_ORDER) * ROW_H_MM + LEGEND_H_MM + MARGIN_MM
    )
    require(
        abs(width_mm - WIDTH_MM) < 0.05 and abs(height_mm - expected_height_mm) < 0.05,
        f"Figure 3B PDF dimensions changed: {width_mm:.2f} × {height_mm:.2f} mm",
    )
    visible = subprocess.run(
        [pdftotext, "-layout", str(pdf), "-"],
        check=True, capture_output=True, text=True,
    ).stdout
    visible_flat = re.sub(r"\s+", " ", visible)
    required = [
        "LA compartment localization and carrier-relative classification",
        "1 Fidelity",
        "2 Relay",
        "1 Extension",
        "12 Below-detection",
        "AF–control family Δsingscore",
        "Global symmetric scale; no row standardization",
        "healthy HCA LA carrier",
        "nominal family-level permutation p < 0.05",
        "Nominal support is a localization cue",
        "Carrier-relative",
        "class",
        "Below",
        "18/16",
        "16/15",
        "16/10",
        "15/13",
        "12/10",
    ]
    missing = [p for p in required if p not in visible_flat]
    require(not missing, f"Missing visible phrases: {missing}")
    forbidden = ["PriorityScore", "p < 0.10", "histogram", "NCP_v3"]
    hits = [t for t in forbidden if t in visible]
    require(not hits, f"Forbidden text visible: {hits}")
    svg = OUT_BASE.with_suffix(".svg")
    png = OUT_BASE.with_suffix(".png")
    require(svg.is_file() and svg.stat().st_size > 0, "SVG export missing")
    svg_text = svg.read_text(encoding="utf-8")
    require(
        "Global symmetric scale; no row standardization" in svg_text,
        "Global scale note missing or clipped in SVG",
    )
    image = plt.imread(png)
    expected_width_px = round(WIDTH_MM / 25.4 * 600)
    expected_height_px = round(expected_height_mm / 25.4 * 600)
    require(
        abs(image.shape[1] - expected_width_px) <= 1
        and abs(image.shape[0] - expected_height_px) <= 1,
        f"600-dpi PNG dimensions changed: {image.shape[1]} × {image.shape[0]}",
    )


def main() -> None:
    require(OUT_SOURCE.is_file(), "Approved Figure 3B source snapshot is missing")
    require(
        sha256(OUT_SOURCE) == APPROVED_SOURCE_SHA256,
        "Approved Figure 3B source changed before layout-only rendering",
    )
    delta, cls, fam_order, dir_lookup, label_lookup, carrier_lookup = load_tables()
    assert_locked_results(delta, cls)
    source, m_scale = build_source(delta, cls, fam_order, dir_lookup, label_lookup, carrier_lookup)
    assert_approved_source(source, m_scale)
    output_hashes = write_outputs(source, m_scale)
    require(
        sha256(OUT_SOURCE) == APPROVED_SOURCE_SHA256,
        "Layout polish modified the approved Figure 3B source table",
    )
    validate_pdf(OUT_BASE.with_suffix(".pdf"))
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_SOURCE.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_SCALE_QC.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")
    print(f"scale_M={m_scale:.6g}; outputs={len(output_hashes)}")


if __name__ == "__main__":
    main()
