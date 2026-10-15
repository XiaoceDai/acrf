"""Redraw Frontiers Figure 1 from the BMC reference figure and the frozen M table.

Panel A states the paired contrast. Panel B shows the healthy reference:
eight donor-level aligned aCAI values and the program-specific magnitude M.
Donor coordinates are read from the vector marks in the BMC Figure 1 and
checked against migration_M_program.tsv. Median M must be 0.0284.
"""

import io
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymupdf
from matplotlib.lines import Line2D

ROOT = Path("/home/bio/Project_P02")
BMC_FIG = ROOT / "submission/bmc_bioinformatics/final_pack/figures/Figure_1.pdf"
PROG_ZIP = ROOT / "submission/frontiers_bioinformatics/supplementary/Supplementary_program_definitions.zip"
OUT = ROOT / "submission/frontiers_bioinformatics/figures/Figure_1.pdf"

INK = "#1C1C1C"
MUTED = "#5C5C5C"
RULE = "#D0D0D0"
LA = "#C97A74"
RA = "#7D99B8"
DONOR = "#858585"

# Row order is the frozen canonical display order, top to bottom as in BMC Figure 1B.
LABELS = [
    "REACTOME__REACTOME_TRANSLATION",
    "REACTOME__REACTOME_RIBOSOME_ASSOCIATED_QUALITY_CONTROL",
    "REACTOME__REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT",
    "HALLMARK__HALLMARK_ADIPOGENESIS",
    "REACTOME__REACTOME_TRIGLYCERIDE_METABOLISM",
    "HALLMARK__HALLMARK_FATTY_ACID_METABOLISM",
    "REACTOME__REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION",
    "REACTOME__REACTOME_PROTEIN_LOCALIZATION",
    "REACTOME__REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
    "REACTOME__REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM",
    "REACTOME__REACTOME_TRNA_AMINOACYLATION",
    "REACTOME__REACTOME_MITOCHONDRIAL_TRANSLATION",
    "REACTOME__REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT",
    "REACTOME__REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
    "HALLMARK__HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "REACTOME__REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
    "REACTOME__REACTOME_COMPLEX_IV_ASSEMBLY",
    "REACTOME__REACTOME_COMPLEX_III_ASSEMBLY",
    "REACTOME__REACTOME_COMPLEX_I_BIOGENESIS",
    "REACTOME__REACTOME_RRNA_PROCESSING",
    "REACTOME__REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
    "REACTOME__REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
    "REACTOME__REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM",
    "REACTOME__REACTOME_RUNX2_REGULATES_BONE_DEVELOPMENT",
    "REACTOME__REACTOME_CLASS_B_2_SECRETIN_FAMILY_RECEPTORS",
    "REACTOME__REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM",
    "HALLMARK__HALLMARK_TGF_BETA_SIGNALING",
    "REACTOME__REACTOME_SIGNALING_BY_BMP",
    "REACTOME__REACTOME_SIGNALING_BY_NTRK2_TRKB",
    "REACTOME__REACTOME_KERATAN_SULFATE_BIOSYNTHESIS",
]
SHORT = {
    "REACTOME__REACTOME_TRANSLATION": "Translation, cytosolic",
    "REACTOME__REACTOME_RIBOSOME_ASSOCIATED_QUALITY_CONTROL": "Ribosome quality control",
    "REACTOME__REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT": "Mitochondrial Ca²⁺ transport",
    "HALLMARK__HALLMARK_ADIPOGENESIS": "Adipogenesis",
    "REACTOME__REACTOME_TRIGLYCERIDE_METABOLISM": "Triglyceride metabolism",
    "HALLMARK__HALLMARK_FATTY_ACID_METABOLISM": "Fatty-acid metabolism",
    "REACTOME__REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION": "Mitochondrial proteostasis",
    "REACTOME__REACTOME_PROTEIN_LOCALIZATION": "Protein localization",
    "REACTOME__REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION": "Mitochondrial β-oxidation",
    "REACTOME__REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM": "BCAA catabolism",
    "REACTOME__REACTOME_TRNA_AMINOACYLATION": "tRNA aminoacylation",
    "REACTOME__REACTOME_MITOCHONDRIAL_TRANSLATION": "Mitochondrial translation",
    "REACTOME__REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT": "Mitochondrial protein import",
    "REACTOME__REACTOME_RESPIRATORY_ELECTRON_TRANSPORT": "Electron transport",
    "HALLMARK__HALLMARK_OXIDATIVE_PHOSPHORYLATION": "Oxidative phosphorylation",
    "REACTOME__REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE": "TCA cycle",
    "REACTOME__REACTOME_COMPLEX_IV_ASSEMBLY": "Complex IV assembly",
    "REACTOME__REACTOME_COMPLEX_III_ASSEMBLY": "Complex III assembly",
    "REACTOME__REACTOME_COMPLEX_I_BIOGENESIS": "Complex I biogenesis",
    "REACTOME__REACTOME_RRNA_PROCESSING": "rRNA processing",
    "REACTOME__REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES": "Amino-acid metabolism",
    "REACTOME__REACTOME_EUKARYOTIC_TRANSLATION_INITIATION": "Translation initiation",
    "REACTOME__REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM": "HS-GAG metabolism",
    "REACTOME__REACTOME_RUNX2_REGULATES_BONE_DEVELOPMENT": "RUNX2, bone development",
    "REACTOME__REACTOME_CLASS_B_2_SECRETIN_FAMILY_RECEPTORS": "Secretin-family receptors",
    "REACTOME__REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM": "GAG metabolism, disease",
    "HALLMARK__HALLMARK_TGF_BETA_SIGNALING": "TGF-β signalling",
    "REACTOME__REACTOME_SIGNALING_BY_BMP": "BMP signalling",
    "REACTOME__REACTOME_SIGNALING_BY_NTRK2_TRKB": "NTRK2 / TrkB signalling",
    "REACTOME__REACTOME_KERATAN_SULFATE_BIOSYNTHESIS": "Keratan sulfate biosynthesis",
}


def _table() -> pd.DataFrame:
    with zipfile.ZipFile(PROG_ZIP) as zf:
        prog = pd.read_csv(io.BytesIO(zf.read("frozen_program_metadata.tsv")), sep="\t")
        mtab = pd.read_csv(io.BytesIO(zf.read("migration_M_program.tsv")), sep="\t")
    df = prog.merge(mtab[["program_id", "M"]], on="program_id", how="left", validate="one_to_one")
    df = df.set_index("program_id").loc[LABELS].reset_index()
    if abs(float(df["M"].median()) - 0.0284) > 0.0002:
        raise SystemExit(f"median M is {df['M'].median()}")
    if list(df["canonical_display_order"]) != list(range(1, 31)):
        raise SystemExit("canonical display order is not 1..30 in label order")
    return df


def _donor_marks(df: pd.DataFrame) -> list[dict]:
    """Read the 240 donor dots and 30 magnitude diamonds from BMC Figure 1B."""
    page = pymupdf.open(BMC_FIG)[0]
    label_y = []
    for block in page.get_text("dict")["blocks"]:
        if block.get("type") != 0:
            continue
        for line in block["lines"]:
            bb = line["bbox"]
            text = "".join(span["text"] for span in line["spans"]).strip()
            if abs(bb[0] - 279.1) < 2 and 100 < bb[1] < 275:
                label_y.append(((bb[1] + bb[3]) / 2, text))
    label_y.sort()
    if len(label_y) != 30:
        raise SystemExit(f"expected 30 program labels in BMC Figure 1, found {len(label_y)}")

    x0, x1 = 378.0, 407.81
    scale = (x1 - x0) / 0.1
    rows = []
    for (y, _pdf_label), rec in zip(label_y, df.itertuples(index=False)):
        donors = []
        diamonds = []
        for drawing in page.get_drawings():
            rect = drawing["rect"]
            fill = drawing.get("fill")
            if rect.is_empty or fill is None or rect.width > 12:
                continue
            cx = (rect.x0 + rect.x1) / 2
            cy = (rect.y0 + rect.y1) / 2
            if abs(cy - y) > 3.2 or not (370 < cx < 510):
                continue
            value = (cx - x0) / scale
            if 3.1 < rect.width < 3.36:
                donors.append(value)
            elif rect.width > 7.0:
                kind = "LA-biased" if fill[0] > 0.7 else "RA-biased"
                diamonds.append((value, kind))
        if len(donors) != 8 or len(diamonds) != 1:
            raise SystemExit(f"{rec.program_id}: donors={len(donors)} diamonds={len(diamonds)}")
        value, kind = diamonds[0]
        if abs(value - float(rec.M)) > 0.001:
            raise SystemExit(f"{rec.program_id}: drawn M {value} != table {rec.M}")
        if kind != rec.canonical_direction:
            raise SystemExit(f"{rec.program_id}: drawn direction {kind}")
        if min(donors) <= 0:
            raise SystemExit(f"{rec.program_id}: non-positive donor aCAI")
        rows.append(
            {
                "program_id": rec.program_id,
                "name": SHORT[rec.program_id],
                "family": int(rec.family_id),
                "direction": rec.canonical_direction,
                "M": float(rec.M),
                "donors": donors,
            }
        )
    if sum(row["direction"] == "LA-biased" for row in rows) != 22:
        raise SystemExit("expected 22 LA-biased programs")
    return rows


def _panel_a(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.0, 0.985, "A", fontsize=12, fontweight="bold", color=INK, ha="left", va="top")
    ax.text(0.09, 0.985, "Paired contrast", fontsize=9.5, fontweight="bold", color=INK, ha="left", va="top")
    ax.plot([0.0, 0.98], [0.945, 0.945], color=RULE, lw=0.6)

    ax.text(0.55, 0.905, "Non-failing reference donors", fontsize=8, fontweight="bold", color=INK, ha="center")
    ax.text(0.38, 0.855, "LA score", fontsize=7.2, color=LA, fontweight="bold", ha="center")
    ax.text(0.78, 0.855, "RA score", fontsize=7.2, color=RA, fontweight="bold", ha="center")

    def pair(y, label):
        ax.text(0.0, y, label, fontsize=7.4, color=INK, ha="left", va="center")
        ax.plot([0.42, 0.74], [y, y], color="#B5B5B5", lw=0.8, zorder=1)
        ax.scatter([0.38], [y], s=155, marker="s", color=LA, zorder=2)
        ax.scatter([0.78], [y], s=155, marker="s", color=RA, zorder=2)
        ax.text(0.38, y, "LA", fontsize=6.2, color="white", fontweight="bold", ha="center", va="center", zorder=3)
        ax.text(0.78, y, "RA", fontsize=6.2, color="white", fontweight="bold", ha="center", va="center", zorder=3)

    pair(0.785, "Donor 1")
    pair(0.710, "Donor 2")
    pair(0.635, "Donor 3")
    ax.text(0.55, 0.575, "⋮", fontsize=12, color=MUTED, ha="center", va="center")
    ax.text(0.0, 0.528, "8 donors, 16 tissue samples", fontsize=7.6, fontweight="bold", color=INK, ha="left", va="center")
    ax.text(0.0, 0.488, "Paired LA and RA", fontsize=7.6, fontweight="bold", color=INK, ha="left", va="center")
    ax.text(0.0, 0.448, "GSE112339", fontsize=7.2, color=MUTED, ha="left", va="center")

    ax.plot([0.0, 0.98], [0.435, 0.435], color=RULE, lw=0.6)
    equations = [
        r"$\Delta$  =  LA − RA",
        r"aCAI  =  $s$ × $\Delta$",
        r"$M$  =  mean |aCAI|",
        r"aCAIrel  =  aCAI / $M$",
    ]
    for i, eq in enumerate(equations):
        ax.text(0.0, 0.385 - i * 0.055, eq, fontsize=8.2, color=INK, ha="left", va="center")

    ax.text(
        0.0,
        0.115,
        "s = +1 for an LA-biased reference,\n"
        "and −1 for an RA-biased reference.\n"
        "Each contrast uses the donor's paired samples.\n"
        "Panel B shows the eight reference donors.",
        fontsize=7.0,
        color=MUTED,
        ha="left",
        va="center",
        linespacing=1.35,
    )


def _panel_b(ax, rows: list[dict]) -> None:
    n = len(rows)
    jit = np.linspace(-0.22, 0.22, 8)
    for i, row in enumerate(rows):
        y = n - 1 - i
        color = LA if row["direction"] == "LA-biased" else RA
        donors = sorted(row["donors"])
        ax.plot([min(donors), max(donors)], [y, y], color="#C8C8C8", lw=1.1, zorder=1, solid_capstyle="round")
        ax.scatter(donors, y + jit, s=9, c=DONOR, linewidths=0, zorder=2)
        ax.scatter([row["M"]], [y], s=28, marker="D", facecolors=color, edgecolors="white", linewidths=0.4, zorder=3)
        ax.text(-0.02, y, row["name"], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.7, color=INK, clip_on=False)
        show_family = i == 0 or row["family"] != rows[i - 1]["family"]
        if show_family:
            ax.text(-0.62, y, f"F{row['family']}", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.7, color=color, fontweight="bold", clip_on=False)

    ax.set_ylim(-0.8, n - 0.2)
    ax.set_xlim(0, 0.42)
    ax.set_yticks([])
    ax.set_xlabel("Aligned reference contrast, aCAI", fontsize=8)
    ax.tick_params(axis="x", labelsize=7.2, length=3, width=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(RULE)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])

    ax.text(-0.62, 1.055, "B", fontsize=12, fontweight="bold", color=INK, transform=ax.transAxes, ha="left", va="bottom")
    ax.text(-0.50, 1.055, "Non-failing reference, 30 programs", fontsize=9.5, fontweight="bold", color=INK, transform=ax.transAxes, ha="left", va="bottom")
    ax.text(
        -0.62,
        1.012,
        "22 LA-biased, 8 RA-biased. Median M = 0.0284.",
        fontsize=6.4,
        color=MUTED,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
    )
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DONOR, markersize=4.2, label="Donor"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=LA, markersize=5.2, label="M, LA-biased"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=RA, markersize=5.2, label="M, RA-biased"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=6.4, loc="upper right", borderpad=0.1, handletextpad=0.4)


def main() -> None:
    df = _table()
    rows = _donor_marks(df)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.6,
        "axes.edgecolor": INK,
        "text.color": INK,
        "xtick.color": INK,
    })
    fig = plt.figure(figsize=(7.40, 6.55))
    ax_a = fig.add_axes([0.02, 0.05, 0.27, 0.90])
    ax_b = fig.add_axes([0.60, 0.08, 0.36, 0.80])
    _panel_a(ax_a)
    _panel_b(ax_b, rows)
    fig.savefig(OUT)
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
