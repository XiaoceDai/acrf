#!/usr/bin/env python3
"""
Assemble final Supplementary Figure S3 — Full LA single-nucleus family and program localization.

Phase 5B3A panel roles (authoritative for this assembly):
  A — 16-family × 9-compartment corrected localization map
  B — 30-program × 9-compartment corrected localization map
  C — corrected class / compartment detail (1/2/1/12; F15 endothelial only)

Phase 5B2 candidate PDFs (stability / stacked fullmatrix / VHS) do not match these
roles and are NOT used as A–C panels. Numerical values are read from the same
corrected frozen fig3_v3_* / fidelity-class lineage (no scientific recomputation).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "supplementary"))

from _pathway_labels import short_pathway_label  # noqa: E402
from _export import save_fig_bundle, set_pub_style  # noqa: E402

FAM_TSV = ROOT / "results" / "fig3_v3_family_celltype_delta.tsv"
PROG_TSV = ROOT / "results" / "fig3_v3_program_celltype_delta.tsv"
CLASS_TSV = ROOT / "results" / "fig3c_family_fidelity_classes.tsv"
ORDER_TSV = ROOT / "results" / "ncp_v3_canonical_display_order.tsv"
LABEL_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
COMPARTMENTS = [
    "atrialCM", "fibroblast", "endothelial", "macrophage", "lymphocyte",
    "pericyte", "vascular_smooth_muscle", "adipocyte", "neuronal",
]
COMP_SHORT = {
    "atrialCM": "aCM", "fibroblast": "Fibro", "endothelial": "Endo",
    "macrophage": "Macro", "lymphocyte": "Lymph", "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC", "adipocyte": "Adipo", "neuronal": "Neuro",
}
CLASS_COLORS = {
    "Fidelity": "#2E7D32",
    "Relay": "#F9A825",
    "Extension": "#6A1B9A",
    "Below_detection": "#BDBDBD",
}


def draw_heatmap(ax, mat, pmat, vmax, row_labels, col_labels):
    cmap = LinearSegmentedColormap.from_list(
        "AF_d_supp", [(0.0, RA_DEEP), (0.5, "#FFFFFF"), (1.0, LA_DEEP)]
    )
    cmap.set_bad("#E0E0E0")
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    masked = np.ma.masked_invalid(mat)
    n_rows, n_cols = mat.shape
    im = ax.imshow(
        masked, cmap=cmap, norm=norm, aspect="auto",
        extent=[-0.5, n_cols - 0.5, n_rows - 0.5, -0.5],
    )
    for i in range(n_rows):
        for j in range(n_cols):
            p = pmat[i, j]
            if not np.isfinite(p):
                continue
            if p < 0.05:
                ax.add_patch(
                    mpatches.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1,
                        facecolor="none", edgecolor="#1A1A1A", linewidth=0.9,
                    )
                )
                ax.text(j, i + 0.12, "*", ha="center", va="center",
                        fontsize=7.5, fontweight="bold", color="#1A1A1A")
            elif p < 0.10:
                ax.text(j, i + 0.12, ".", ha="center", va="center",
                        fontsize=8, fontweight="bold", color="#444444")
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([COMP_SHORT[c] for c in col_labels], rotation=35, ha="right", fontsize=6.2)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels, fontsize=5.6)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def assert_classes(cls: pd.DataFrame) -> None:
    vc = cls["fidelity_class"].value_counts()
    assert int(vc.get("Fidelity", 0)) == 1
    assert int(vc.get("Relay", 0)) == 2
    assert int(vc.get("Extension", 0)) == 1
    assert int(vc.get("Below_detection", 0)) == 12
    f10 = cls.loc[cls.family_id == 10].iloc[0]
    f11 = cls.loc[cls.family_id == 11].iloc[0]
    f15 = cls.loc[cls.family_id == 15].iloc[0]
    f16 = cls.loc[cls.family_id == 16].iloc[0]
    f4 = cls.loc[cls.family_id == 4].iloc[0]
    assert f10.fidelity_class == "Fidelity" and str(f10.snRNA_supported_compartment) == "fibroblast"
    assert f11.fidelity_class == "Relay" and str(f11.snRNA_supported_compartment) == "adipocyte"
    assert f15.fidelity_class == "Relay" and str(f15.snRNA_supported_compartment) == "endothelial"
    assert f16.fidelity_class == "Extension"
    comps = {c.strip() for c in str(f16.snRNA_supported_compartment).split(",")}
    assert comps == {"fibroblast", "atrialCM"}
    assert f4.fidelity_class == "Below_detection"


def main() -> dict:
    set_pub_style()
    fam = pd.read_csv(FAM_TSV, sep="\t")
    prog = pd.read_csv(PROG_TSV, sep="\t")
    cls = pd.read_csv(CLASS_TSV, sep="\t")
    canonical = pd.read_csv(ORDER_TSV, sep="\t")
    labels = pd.read_csv(LABEL_TSV, sep="\t")
    assert_classes(cls)

    fam_order = (
        canonical[["family_display_order", "family_id"]]
        .drop_duplicates("family_id")
        .sort_values("family_display_order")["family_id"]
        .astype(int)
        .tolist()
    )
    assert len(fam_order) == 16
    assert set(COMPARTMENTS) <= set(fam["compartment"].unique())

    fam_pivot = fam.pivot(index="family_id", columns="compartment", values="dScore").reindex(
        index=fam_order, columns=COMPARTMENTS
    )
    fam_pmat = fam.pivot(index="family_id", columns="compartment", values="p_perm").reindex(
        index=fam_order, columns=COMPARTMENTS
    )
    from _export import ascii_pub_label

    fam_short = {
        int(r.family_id): ascii_pub_label(short_pathway_label(str(r.family_label), 34))
        for r in labels.itertuples(index=False)
    }
    # Prefer short_label when present (publication-facing), still ASCII-safe
    if "short_label" in labels.columns:
        for r in labels.itertuples(index=False):
            sl = str(r.short_label).strip()
            if sl and sl.lower() != "nan":
                fam_short[int(r.family_id)] = ascii_pub_label(short_pathway_label(sl, 34))
    fam_labels = [f"F{fid}  {fam_short.get(fid, '')}" for fid in fam_order]

    prog = prog.copy()
    prog["family_id"] = prog["family_id"].astype(int)
    prog["program_order"] = prog["program_order"].astype(int)
    prog["fam_rank"] = prog["family_id"].map({fid: i for i, fid in enumerate(fam_order)}).fillna(99).astype(int)
    prog["program_short"] = prog["program_id"].map(
        lambda s: ascii_pub_label(short_pathway_label(str(s), 32))
    )
    prog_order = (
        prog.drop_duplicates("program_id")
        .sort_values(["fam_rank", "program_order"])["program_id"]
        .tolist()
    )
    assert len(prog_order) == 30
    prog_pivot = prog.pivot(index="program_id", columns="compartment", values="dScore").reindex(
        index=prog_order, columns=COMPARTMENTS
    )
    prog_pmat = prog.pivot(index="program_id", columns="compartment", values="p_perm").reindex(
        index=prog_order, columns=COMPARTMENTS
    )
    prog_lab = dict(zip(prog["program_id"], prog["program_short"]))
    prog_fam = dict(zip(prog["program_id"], prog["family_id"]))
    prog_labels = [f"F{prog_fam[p]}  {prog_lab[p]}" for p in prog_order]

    fam_arr = fam_pivot.to_numpy(dtype=float)
    prog_arr = prog_pivot.to_numpy(dtype=float)
    vmax_fam = max(0.04, float(np.nanmax(np.abs(fam_arr))) * 1.05)
    vmax_prog = max(0.04, float(np.nanmax(np.abs(prog_arr))) * 1.05)

    MM = 1 / 25.4
    fig = plt.figure(figsize=(190 * MM, 300 * MM))
    gs = fig.add_gridspec(
        3, 2, height_ratios=[2.4, 4.6, 1.7], width_ratios=[12, 0.55],
        left=0.28, right=0.97, top=0.93, bottom=0.05, hspace=0.22, wspace=0.05,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a_cb = fig.add_subplot(gs[0, 1])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_b_cb = fig.add_subplot(gs[1, 1])
    ax_c = fig.add_subplot(gs[2, :])

    im_a = draw_heatmap(ax_a, fam_arr, fam_pmat.to_numpy(dtype=float), vmax_fam, fam_labels, COMPARTMENTS)
    cba = fig.colorbar(im_a, cax=ax_a_cb)
    cba.set_label("Δsingscore (AF − Ctrl)", fontsize=6.0)
    cba.ax.tick_params(labelsize=5.2)
    ax_a.set_title("A  16 NCP families × 9 LA compartments", loc="left", fontsize=7.5, fontweight="bold", pad=3)

    im_b = draw_heatmap(ax_b, prog_arr, prog_pmat.to_numpy(dtype=float), vmax_prog, prog_labels, COMPARTMENTS)
    cbb = fig.colorbar(im_b, cax=ax_b_cb)
    cbb.set_label("Δsingscore (AF − Ctrl)", fontsize=6.0)
    cbb.ax.tick_params(labelsize=5.2)
    ax_b.set_title("B  30 NCP programs × 9 LA compartments", loc="left", fontsize=7.5, fontweight="bold", pad=3)

    # Panel C — class table
    ax_c.axis("off")
    ax_c.set_title(
        "C  Carrier-relative classes (nominal permutation; LA-only; not FDR-controlled)",
        loc="left", fontsize=7.5, fontweight="bold", pad=2,
    )
    cls_ord = cls.set_index("family_id").loc[fam_order].reset_index()
    rows = []
    for r in cls_ord.itertuples(index=False):
        rows.append([
            f"F{int(r.family_id)}",
            str(r.fidelity_class).replace("_", "-"),
            str(r.snRNA_supported_compartment) if pd.notna(r.snRNA_supported_compartment) and str(r.snRNA_supported_compartment) not in {"", "nan"} else "—",
        ])
    table = ax_c.table(
        cellText=rows,
        colLabels=["Family", "Class", "Supported compartment(s)"],
        cellLoc="left",
        colLoc="left",
        bbox=[0.0, 0.02, 1.0, 0.90],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(5.4)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#CCCCCC")
        cell.set_linewidth(0.3)
        if r == 0:
            cell.set_facecolor("#F0F0F0")
            cell.set_text_props(fontweight="bold")
        elif c == 1 and r > 0:
            cls_name = rows[r - 1][1].replace("-", "_")
            # map display back
            key = rows[r - 1][1]
            key = {"Below-detection": "Below_detection"}.get(key, key)
            cell.set_facecolor(CLASS_COLORS.get(key, "#FFFFFF"))

    fig.suptitle(
        "Full LA single-nucleus family and program localization",
        x=0.28, y=0.975, ha="left", fontsize=9.0, fontweight="bold",
    )
    fig.text(
        0.28, 0.955,
        "NCP · GSE255612 LA · exact-reference singscore lineage · "
        "Fidelity F10 (fibroblast); Relay F11 (adipocyte), F15 (endothelial only); "
        "Extension F16 (fibroblast + atrialCM); Below-detection 12 · "
        "F4 remains Below-detection (no fifth non-Below class) · "
        "nominal p (boxed * <0.05; dot <0.10); not FDR discovery · LA-only",
        fontsize=5.4, color="#333333", ha="left", va="top", wrap=True,
    )

    paths = save_fig_bundle(fig, "figS3_la_single_nucleus_full_localization", pad_inches=0.03)
    return {
        "paths": paths,
        "panel_structure": "A=16x9 family; B=30x9 program; C=class/compartment",
        "notes": (
            "Phase 5B2 PDF candidates (stability/fullmatrix/VHS) remapped: "
            "assembled A–C from corrected frozen fig3_v3_* and fidelity-class sources"
        ),
    }


if __name__ == "__main__":
    info = main()
    print("wrote", info["paths"]["pdf"])
