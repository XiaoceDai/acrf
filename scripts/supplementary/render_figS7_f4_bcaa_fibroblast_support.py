#!/usr/bin/env python3
"""
Final Supplementary Figure S7 — F4 BCAA constituent-program fibroblast support.

Reads corrected exact-reference program×compartment delta table (frozen).
No scientific recomputation. Publication-clean boundaries.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "supplementary"))
from _export import save_fig_bundle, set_pub_style  # noqa: E402

PROG_TSV = ROOT / "results" / "fig3_v3_program_celltype_delta.tsv"
CLASS_TSV = ROOT / "results" / "fig3c_family_fidelity_classes.tsv"
HCA_TSV = ROOT / "results" / "fig1f_v3_HCA_family_carriers.tsv"

PROGRAM_ID = "REACTOME__REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM"
LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
COMPARTMENT_SHORT = {
    "atrialCM": "aCM",
    "fibroblast": "Fibro",
    "endothelial": "Endo",
    "macrophage": "Macro",
    "lymphocyte": "Lymph",
    "pericyte": "Peri",
    "vascular_smooth_muscle": "vSMC",
    "adipocyte": "Adipo",
    "neuronal": "Neuro",
}


def main() -> dict:
    set_pub_style()
    df = pd.read_csv(PROG_TSV, sep="\t")
    df = df[df["program_id"] == PROGRAM_ID].copy()
    df = df.sort_values("dScore", ascending=False).reset_index(drop=True)
    assert len(df) == 9

    cls = pd.read_csv(CLASS_TSV, sep="\t")
    f4 = cls.loc[cls.family_id == 4].iloc[0]
    assert f4.fidelity_class == "Below_detection"

    fibro = df.loc[df.compartment == "fibroblast"].iloc[0]
    p_fib = float(fibro["p_perm"])
    d_fib = float(fibro["dScore"])
    assert p_fib < 0.05

    hca = pd.read_csv(HCA_TSV, sep="\t")
    family_carrier = str(hca.loc[hca["family_id"] == 4, "LA_carrier"].iloc[0])

    MM = 1 / 25.4
    fig = plt.figure(figsize=(170 * MM, 115 * MM))
    ax = fig.add_axes([0.18, 0.34, 0.78, 0.48])
    for i, r in df.iterrows():
        v = float(r["dScore"])
        p = float(r["p_perm"])
        ct = str(r["compartment"])
        color = LA_DEEP if v >= 0 else RA_DEEP
        ax.barh(i, v, height=0.65, color=color, edgecolor="#1A1A1A", linewidth=0.4, alpha=0.92)
        if ct == family_carrier:
            ax.text(
                -0.024, i, "< family HCA carrier",
                fontsize=5.4, color="#888888", style="italic", ha="left", va="center",
            )
        if p < 0.05:
            ax.text(
                v + (0.001 if v >= 0 else -0.001), i,
                f"  *  p={p:.3f}",
                fontsize=6.2, fontweight="bold", color="#1A1A1A",
                ha="left" if v >= 0 else "right", va="center",
            )
        elif p < 0.10:
            ax.text(
                v + (0.001 if v >= 0 else -0.001), i,
                f"  .  p={p:.3f}",
                fontsize=5.6, color="#444444",
                ha="left" if v >= 0 else "right", va="center",
            )
    ax.axvline(0, color="#1A1A1A", lw=0.5)
    ax.set_yticks(list(range(len(df))))
    ax.set_yticklabels([COMPARTMENT_SHORT[c] for c in df["compartment"]], fontsize=7.0, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Δsingscore (AF − Ctrl)", fontsize=7.2)
    ax.set_xlim(-0.045, 0.045)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.suptitle(
        "F4 BCAA constituent-program fibroblast support",
        x=0.05, y=0.95, ha="left", fontsize=9.0, fontweight="bold",
    )
    fig.text(
        0.05, 0.89,
        "Program-level example only · parent F4 family remains Below-detection · "
        "not a fifth non-Below family · not FDR-controlled discovery · "
        "does not alter Figure 3 family classes or Figure 5 support counts",
        fontsize=5.5, color="#333333", ha="left",
    )
    fig.text(
        0.05, 0.20,
        f"Corrected exact-reference singscore: BCAA constituent program "
        f"(REACTOME branched-chain amino-acid catabolism) shows nominal fibroblast "
        f"support (Δ={d_fib:.4f}, p_perm={p_fib:.4f}; n_ctrl={int(fibro.n_ctrl)}, "
        f"n_af={int(fibro.n_af)}). Family-level F4 class = Below-detection. "
        f"HCA LA carrier prior for F4 = {COMPARTMENT_SHORT.get(family_carrier, family_carrier)}. "
        "Color encodes direction only (red AF>Ctrl; blue AF<Ctrl), not fidelity class.",
        fontsize=5.5, color="#333333", ha="left", wrap=True,
    )

    paths = save_fig_bundle(fig, "figS7_f4_bcaa_fibroblast_support", pad_inches=0.03)
    return {
        "paths": paths,
        "panel_structure": "single program×compartment bar panel",
        "p_fibroblast": p_fib,
        "d_fibroblast": d_fib,
        "notes": "Rerendered from corrected fig3_v3_program_celltype_delta.tsv with publication boundaries",
    }


if __name__ == "__main__":
    print(main()["paths"]["pdf"])
