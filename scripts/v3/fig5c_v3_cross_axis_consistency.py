#!/usr/bin/env python3
"""Fig 5c v3 - Cross-axis convergence classification (Path A strict)."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_style import enforce_font_floor, FS_PANEL_TITLE, FS_TICK, FS_FOOTNOTE

IN_TSV = ROOT / "results" / "fig5b_v3_triangulation.tsv"
LABELS_TSV = ROOT / "results" / "ncp_v3_family_labels.tsv"
OUT_TSV = ROOT / "results" / "fig5c_v3_cross_axis.tsv"
OUT_PDF = ROOT / "figures" / "fig5c_cross_axis.pdf"

CLASS_ORDER = ["Triple-positive", "Double-positive", "Single-axis", "Silent"]
CLASS_COLORS = {
    "Triple-positive":  "#2E7B5E",
    "Double-positive":  "#7BB36C",
    "Single-axis":      "#C8C39A",
    "Silent":           "#BBBBBB",
}
SUPPORT_YELLOW = "#F5C518"
LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"


def manuscript_cross_axis_table(full: pd.DataFrame) -> pd.DataFrame:
    m = full.copy()
    if LABELS_TSV.is_file() and ("family_label" not in m.columns or "axis_group" not in m.columns):
        labs = pd.read_csv(LABELS_TSV, sep="\t")[["family_id", "family_label", "axis_group"]].drop_duplicates(
            subset=["family_id"], keep="first"
        )
        for col in ("family_label", "axis_group"):
            if col in m.columns:
                m = m.drop(columns=[col])
        m = m.merge(labs, on="family_id", how="left")
    if "family_label" not in m.columns:
        m["family_label"] = ""
    if "axis_group" not in m.columns:
        m["axis_group"] = ""
    out = pd.DataFrame({
        "family_id": m["family_id"].astype(int),
        "family_label": m["family_label"],
        "axis_group": m["axis_group"],
        "bulk_bit": m["L1_supportive"].astype(int),
        "snRNA_bit": m["L2_supportive"].astype(int),
        "duration_bit": m["L3_supportive"].astype(int),
        "cross_axis_class": m["cross_axis_class"],
    })
    out["support_pattern"] = (
        out["bulk_bit"].astype(str) + out["snRNA_bit"].astype(str) + out["duration_bit"].astype(str)
    )
    return out[
        ["family_id", "family_label", "axis_group", "bulk_bit", "snRNA_bit", "duration_bit",
         "support_pattern", "cross_axis_class"]
    ]


def assign_class(row: pd.Series) -> str:
    n = int(row["L1_supportive"]) + int(row["L2_supportive"]) + int(row["L3_supportive"])
    if n == 3:
        return "Triple-positive"
    if n == 2:
        return "Double-positive"
    if n == 1:
        return "Single-axis"
    return "Silent"


def class_subtype(row: pd.Series) -> str:
    """Detail string showing which axes are supportive."""
    parts = []
    for k, lbl in [("L1_supportive", "L1"),
                     ("L2_supportive", "L2"),
                     ("L3_supportive", "L3")]:
        if int(row[k]) == 1:
            parts.append(lbl)
    if not parts:
        return "none"
    return "+".join(parts)


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": enforce_font_floor(FS_TICK),
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def render_pdf(df: pd.DataFrame) -> None:
    set_style()

    n = len(df)
    counts = df["cross_axis_class"].value_counts().reindex(CLASS_ORDER, fill_value=0)
    n_triple = int(counts["Triple-positive"])
    n_double = int(counts["Double-positive"])
    n_single = int(counts["Single-axis"])
    n_silent = int(counts["Silent"])

    MM = 1.0 / 25.4
    fig_w_mm, fig_h_mm = 200, 185
    fig = plt.figure(figsize=(fig_w_mm * MM, fig_h_mm * MM))
    gs = fig.add_gridspec(
        3, 1,
        height_ratios=[1.6, 4.0, 4.4],
        left=0.06, right=0.96, top=0.87, bottom=0.04, hspace=0.55,
    )
    ax_bar = fig.add_subplot(gs[0])
    ax_grid = fig.add_subplot(gs[1])
    ax_list = fig.add_subplot(gs[2])

    ax_bar.set_xlim(0, n)
    ax_bar.set_ylim(0, 1.8)
    ax_bar.axis("off")

    # Section header well above bar
    ax_bar.text(0, 1.75, f"a   Cross-axis convergence  (n = {n} families, Path A strict)",
                  fontsize=enforce_font_floor(FS_PANEL_TITLE), fontweight="bold", color="#1A1A1A", va="top")

    bar_y0, bar_h = 0.18, 0.45
    label_levels = [1.20, 0.95]  # high vs low row, alternated for narrow segments

    x_cur = 0.0
    label_alt = 0
    for cls_name in CLASS_ORDER:
        c = int(counts[cls_name])
        if c == 0:
            continue
        ax_bar.add_patch(mpatches.Rectangle(
            (x_cur, bar_y0), c, bar_h,
            facecolor=CLASS_COLORS[cls_name], edgecolor="white", linewidth=1.2,
        ))
        # Inside the bar: count only (always fits)
        ax_bar.text(x_cur + c / 2.0, bar_y0 + bar_h / 2.0, f"{c}",
                      ha="center", va="center",
                      fontsize=enforce_font_floor(FS_PANEL_TITLE + 1), fontweight="bold",
                      color="white" if cls_name in ("Triple-positive", "Silent") else "#1A1A1A")
        # Class label above bar; stagger height for narrow segments to avoid collision
        if c <= 2:
            ly = label_levels[label_alt]
            label_alt = 1 - label_alt
            label_x = x_cur + c / 2.0
            ax_bar.text(label_x, ly, cls_name,
                          ha="center", va="bottom",
                          fontsize=enforce_font_floor(FS_TICK), fontweight="bold",
                          color=CLASS_COLORS[cls_name])
            ax_bar.plot([label_x, label_x],
                          [ly - 0.02, bar_y0 + bar_h + 0.02],
                          color=CLASS_COLORS[cls_name], lw=0.7)
        else:
            ax_bar.text(x_cur + c / 2.0, bar_y0 - 0.05, cls_name,
                          ha="center", va="top",
                          fontsize=enforce_font_floor(FS_TICK + 0.2), fontweight="bold",
                          color=CLASS_COLORS[cls_name])
        x_cur += c

    # "Path A" definition moved to global subtitle (see fig.text below)

    ax_grid.set_xlim(-0.2, 7.5)
    ax_grid.set_ylim(0, n + 1.0)
    ax_grid.invert_yaxis()
    ax_grid.axis("off")
    ax_grid.text(-0.2, -0.5, "b   Per-family evidence stack",
                   fontsize=enforce_font_floor(FS_PANEL_TITLE), fontweight="bold", color="#1A1A1A", va="bottom")

    # column headers
    col_headers = [("L1", 4.0), ("L2", 4.9), ("L3", 5.8)]
    ax_grid.text(0.0, -0.05, "Family", fontsize=enforce_font_floor(FS_TICK), fontweight="bold", color="#1A1A1A", va="bottom")
    ax_grid.text(2.45, -0.05, "Cross-axis class", fontsize=enforce_font_floor(FS_TICK), fontweight="bold", color="#1A1A1A",
                   ha="center", va="bottom")
    for col_label, col_x in col_headers:
        ax_grid.text(col_x, -0.05, col_label, ha="center", va="bottom",
                       fontsize=enforce_font_floor(FS_TICK), fontweight="bold", color="#1A1A1A")
    ax_grid.text(6.6, -0.05, "n.sup", ha="center", va="bottom",
                   fontsize=enforce_font_floor(FS_TICK), fontweight="bold", color="#1A1A1A")
    ax_grid.text(7.0, -0.05, "Subtype", ha="left", va="bottom",
                   fontsize=enforce_font_floor(FS_TICK), fontweight="bold", color="#1A1A1A")

    for i, row in df.iterrows():
        y = i + 0.5
        y0 = i + 0.5 - 0.40
        # Family ID + short
        bias_color = LA_DEEP if str(row["HealthyBias"]) == "LA-biased" else RA_DEEP
        ax_grid.text(0.0, y, f"F{int(row['family_id'])}",
                       ha="left", va="center", fontsize=enforce_font_floor(FS_TICK), fontweight="bold", color=bias_color)
        ax_grid.text(0.55, y, str(row["short"])[:22],
                       ha="left", va="center", fontsize=enforce_font_floor(FS_TICK - 0.5), color="#1A1A1A")

        # Class block (coloured pill)
        cls_name = str(row["cross_axis_class"])
        ax_grid.add_patch(mpatches.FancyBboxPatch(
            (1.6, y0), 1.7, 0.80,
            boxstyle="round,pad=0.02,rounding_size=0.10",
            facecolor=CLASS_COLORS[cls_name], edgecolor="none",
        ))
        ax_grid.text(2.45, y, cls_name,
                       ha="center", va="center",
                       fontsize=enforce_font_floor(FS_TICK - 0.5), fontweight="bold",
                       color="white" if cls_name in ("Triple-positive", "Silent") else "#1A1A1A")

        # L1 / L2 / L3 supportive bitmap
        for col_x, key in [(4.0, "L1_supportive"),
                              (4.9, "L2_supportive"),
                              (5.8, "L3_supportive")]:
            sup = int(row[key])
            face = SUPPORT_YELLOW if sup == 1 else "#F0F0F0"
            edge = "#1A1A1A" if sup == 1 else "#999999"
            ax_grid.add_patch(mpatches.Rectangle(
                (col_x - 0.30, y0), 0.60, 0.80,
                facecolor=face, edgecolor=edge, linewidth=0.5,
            ))
            if sup == 1:
                ax_grid.text(col_x, y, "+", ha="center", va="center",
                               fontsize=enforce_font_floor(FS_PANEL_TITLE + 1), color="#1A1A1A", fontweight="bold")

        # n.sup
        n_sup = int(row["n_supportive"])
        ax_grid.text(6.6, y, f"{n_sup}/3",
                       ha="center", va="center", fontsize=enforce_font_floor(FS_TICK), fontweight="bold",
                       color="#1A1A1A")
        # Subtype
        ax_grid.text(7.0, y, str(row["cross_axis_subtype"]),
                       ha="left", va="center", fontsize=enforce_font_floor(FS_TICK - 0.5),
                       color="#444444", family="monospace")

    # Light row separators
    for i in range(n + 1):
        ax_grid.plot([-0.2, 7.5], [i, i], color="#EEEEEE", lw=0.3, zorder=0)

    ax_list.set_xlim(0, 1)
    ax_list.set_ylim(0, 1)
    ax_list.invert_yaxis()
    ax_list.axis("off")
    ax_list.text(0.0, 0.02, "c   Members per cross-axis class  (with supportive-axis subtype)",
                   fontsize=enforce_font_floor(FS_PANEL_TITLE), fontweight="bold", color="#1A1A1A")

    y_cur = 0.10
    line_h = 0.065
    for cls_name in CLASS_ORDER:
        sub = df[df["cross_axis_class"] == cls_name].sort_values(
            ["n_supportive", "PriorityScore"], ascending=[False, False])
        n_cls = len(sub)
        # Class header
        ax_list.text(0.00, y_cur,
                       f"\u25A0  {cls_name}  ({n_cls})",
                       fontsize=enforce_font_floor(FS_TICK), fontweight="bold", color=CLASS_COLORS[cls_name])
        threshold_strip = {
            "Triple-positive":  "all three axes supportive",
            "Double-positive":  "two of three axes supportive",
            "Single-axis":      "exactly one axis supportive",
            "Silent":           "no axis supportive  (strict null;  absence of evidence \u2260 evidence of absence)",
        }[cls_name]
        ax_list.text(0.30, y_cur, threshold_strip,
                       fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5), color="#666666", fontstyle="italic")
        y_cur += line_h * 0.85
        if n_cls == 0:
            ax_list.text(0.04, y_cur, "(none)", fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5), color="#888888",
                           family="monospace")
            y_cur += line_h * 0.75
            continue
        # Members per row, 3 per line
        members = []
        for _, r in sub.iterrows():
            fid = int(r["family_id"])
            short = str(r["short"])[:20]
            sub_lbl = str(r["cross_axis_subtype"])
            members.append(f"F{fid:2d} {short:20s} [{sub_lbl}]")
        per_line = 2
        for k in range(0, len(members), per_line):
            chunk = members[k:k + per_line]
            txt = "    ".join(chunk)
            ax_list.text(0.04, y_cur, txt,
                           fontsize=enforce_font_floor(FS_TICK - 0.5), color="#1A1A1A", family="monospace")
            y_cur += line_h * 0.75
        y_cur += line_h * 0.45

    fig.text(0.06, 0.965,
              "Fig 5c  Cross-axis convergence classification  (Path A strict)",
              fontsize=enforce_font_floor(FS_PANEL_TITLE + 2), fontweight="bold", color="#1A1A1A")
    fig.text(0.06, 0.945,
              "Integrative take:  AF molecular endotypes are fragmentary at the family level; "
              f"{n_triple} families reach triple-axis convergence (F10, F15; Path A) "
              f"({n_double} double, {n_single} single, {n_silent} silent).",
              fontsize=enforce_font_floor(FS_TICK), color="#444444", style="italic")
    fig.text(0.06, 0.923,
              "Path A definition:  L1 supportive = replication_class = high (CI excludes 0);  "
              "L2 supportive = fidelity_class \u2260 Below_detection;  "
              "L3 supportive = chronicity_class \u2260 Stable.",
              fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5), color="#555555")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


def main():
    df = pd.read_csv(IN_TSV, sep="\t")
    df["cross_axis_class"] = df.apply(assign_class, axis=1)
    df["cross_axis_subtype"] = df.apply(class_subtype, axis=1)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    slim = manuscript_cross_axis_table(df)
    slim.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(slim)} rows)")
    counts = df["cross_axis_class"].value_counts().reindex(CLASS_ORDER, fill_value=0)
    print("Cross-axis class counts:")
    for cls, c in counts.items():
        print(f"  {cls}:  {c}")
    render_pdf(df)


if __name__ == "__main__":
    main()
