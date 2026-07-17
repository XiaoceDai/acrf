#!/usr/bin/env python3
"""
Final Supplementary Figure S8 — Full NCP family overlap with external AF gene sets.

All 64 tested pairs from frozen fig6c_family_af_overlap_source.tsv.
No scientific recomputation. No continuous FDR color scale.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "supplementary"))
from _export import save_fig_bundle, set_pub_style  # noqa: E402

SRC = ROOT / "results" / "fig6c_family_af_overlap_source.tsv"
ORDER = ROOT / "results" / "ncp_v3_canonical_display_order.tsv"

SET_ORDER = [
    ("yuan2025_cross_population", "Yuan 2025\nPrioritized"),
    ("roselli2018_loci", "Roselli 2018\nRisk genes"),
    ("selewa2023_regulatory", "Selewa 2023\nRegulatory"),
    ("choi_mendelian_rare", "Choi 2025\nRare / SV"),
]


def main() -> dict:
    set_pub_style()
    df = pd.read_csv(SRC, sep="\t")
    assert len(df) == 64
    assert int((df["n_overlap"] > 0).sum()) == 5
    assert int(df["fdr_significant"].sum()) == 0
    assert int(df["universe_n"].iloc[0]) == 14635

    genes = sorted({
        g
        for s in df.loc[df["n_overlap"] > 0, "overlap_genes"].astype(str)
        for g in s.split(";")
        if g and g != "nan"
    })
    assert genes == ["CAV1", "HEY2"]

    # Exact nonzero cells
    want = {
        (3, "yuan2025_cross_population", "CAV1"),
        (3, "roselli2018_loci", "CAV1"),
        (13, "yuan2025_cross_population", "CAV1"),
        (13, "roselli2018_loci", "CAV1"),
        (9, "selewa2023_regulatory", "HEY2"),
    }
    got = set()
    for r in df.itertuples(index=False):
        if int(r.n_overlap) > 0:
            got.add((int(r.family_id), str(r.external_set_id), str(r.overlap_genes)))
    assert got == want

    order = (
        pd.read_csv(ORDER, sep="\t")[["family_display_order", "family_id"]]
        .drop_duplicates("family_id")
        .sort_values("family_display_order")
    )
    fam_ids = order["family_id"].astype(int).tolist()
    assert len(fam_ids) == 16

    from _export import ascii_pub_label

    # short labels from source (ASCII-safe display only; source table untouched)
    lab = {
        int(k): ascii_pub_label(v)
        for k, v in (
            df[["family_id", "family_short_label"]]
            .drop_duplicates("family_id")
            .set_index("family_id")["family_short_label"]
            .to_dict()
            .items()
        )
    }
    # If short labels pass through pathway helper elsewhere, keep final ASCII form
    lab = {k: ascii_pub_label(v) for k, v in lab.items()}

    MM = 1 / 25.4
    fig = plt.figure(figsize=(170 * MM, 150 * MM))
    ax = fig.add_axes([0.28, 0.10, 0.68, 0.78])
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 16)
    ax.invert_yaxis()
    ax.set_xticks([0.5, 1.5, 2.5, 3.5])
    ax.set_xticklabels([s[1] for s in SET_ORDER], fontsize=6.0)
    ax.set_yticks([i + 0.5 for i in range(16)])
    ax.set_yticklabels([f"F{fid}  {lab.get(fid, '')}" for fid in fam_ids], fontsize=6.0)

    # grid
    for i in range(17):
        ax.axhline(i, color="#DDDDDD", lw=0.4)
    for j in range(5):
        ax.axvline(j, color="#DDDDDD", lw=0.4)

    set_ids = [s[0] for s in SET_ORDER]
    for i, fid in enumerate(fam_ids):
        for j, sid in enumerate(set_ids):
            row = df[(df.family_id == fid) & (df.external_set_id == sid)]
            assert len(row) == 1
            n = int(row["n_overlap"].iloc[0])
            if n > 0:
                gene = str(row["overlap_genes"].iloc[0])
                # discrete pill — no continuous significance scale
                ax.add_patch(
                    mpatches.FancyBboxPatch(
                        (j + 0.12, i + 0.18), 0.76, 0.64,
                        boxstyle="round,pad=0.02,rounding_size=0.08",
                        facecolor="#F3E5F5", edgecolor="#6A1B9A", linewidth=0.7,
                    )
                )
                ax.text(j + 0.5, i + 0.5, gene, ha="center", va="center",
                        fontsize=7.0, fontweight="bold", color="#4A148C")

    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("All 64 tested family × external-set pairs", loc="left", fontsize=7.0, pad=4)

    fig.suptitle(
        "Full NCP family overlap with external AF gene sets",
        x=0.05, y=0.96, ha="left", fontsize=9.0, fontweight="bold",
    )
    fig.text(
        0.05, 0.915,
        "16 families × 4 external sets · 5 nonzero cells · genes CAV1, HEY2 only · "
        "0/64 BH-FDR-significant enrichments · universe N = 14,635 · "
        "member overlap ≠ family enrichment · blanks = tested zero overlap · "
        "not a Path A gate",
        fontsize=5.4, color="#333333", ha="left",
    )
    fig.text(
        0.05, 0.03,
        "Exact nonzero cells: F3×Yuan CAV1; F3×Roselli CAV1; F13×Yuan CAV1; "
        "F13×Roselli CAV1; F9×Selewa HEY2. Choi column tested with zero overlap. "
        "Repeated CAV1 cells are not independent genetic findings.",
        fontsize=5.2, color="#444444", ha="left",
    )

    paths = save_fig_bundle(fig, "figS8_full_external_af_gene_overlap", pad_inches=0.02)
    return {"paths": paths, "panel_structure": "16x4 full matrix", "notes": "Frozen source; discrete gene pills only"}


if __name__ == "__main__":
    print(main()["paths"]["pdf"])
