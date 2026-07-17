#!/usr/bin/env python3
"""Fig 4 d/e/f v3 - per-family duration scatter case studies (LAA AF only)."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_style import enforce_font_floor, FS_PANEL_TITLE, FS_TICK, FS_AXIS_LABEL, FS_FOOTNOTE

LONG_TSV = ROOT / "results" / "fig4_v3_singscore_per_sample.tsv"
META_TSV = ROOT / "Project_P4" / "GSE41177" / "sample_metadata.tsv"
CLS_TSV = ROOT / "results" / "fig4_v3_chronicity_class.tsv"
OUT_DIR = ROOT / "figures"

LA_DEEP = "#9C5550"
RA_DEEP = "#4E6B89"
SR_COLOR = "#7D99B8"
AF_COLOR = "#9C5550"

CASES = {
    "fig4d": {
        "family_id": 9,
        "label": "F9 HS-GAG metabolism",
        "class_name": "Strong-trend (4/4 sensitivity)",
        "class_color": "#C44569",
        "panel_letter": "4d",
        "out_pdf": OUT_DIR / "fig4_d_F9_HSGAG.pdf",
        "interpretation_paragraphs": [
            ("Data:",
             "F9 (Heparan / Glycosaminoglycan metabolism, RA-biased baseline) shows a robust negative duration "
             "association (\u03c1 = -0.42, 4/4 sensitivity concordance, all four codings same-sign)."),
            ("Biology:",
             "HS-GAG biosynthetic / catabolic activity declines with AF duration in LAA, compatible with "
             "progressive depletion / saturation of stromal proteoglycan turnover during chronic AF."),
        ],
    },
    "fig4e": {
        "family_id": 2,
        "label": "F2 Mito Ca\u00b2\u207a transport",
        "class_name": "Strong-trend (biological anchor)",
        "class_color": "#C44569",
        "panel_letter": "4e",
        "out_pdf": OUT_DIR / "fig4_e_F2_MitoCa.pdf",
        "interpretation_paragraphs": [
            ("Data:",
             "F2 (Mitochondrial Ca\u00b2\u207a transport, LA-biased baseline) shows the strongest single-family "
             "duration trend (\u03c1 = -0.56, p_perm = 0.059, marginal). Three of four sensitivity codings agree "
             "in direction (3/4); the LAPV coding reverses sign (\u03c1 = +0.20)."),
            ("Literature anchor:",
             "The negative direction matches the Ca\u00b2\u207a-handling remodelling of cardiomyocytes during "
             "sustained AF documented in Voigt et al. 2014 (Circulation), Heijman et al. 2014 (Circ Res) and "
             "Sun et al. 2018, supporting this trend as a genuine progression signal rather than noise."),
        ],
    },
    "fig4f": {
        "family_id": 15,
        "label": "F15 NTRK2 / TrkB signalling",
        "class_name": "Weak-trend (cross-Fig 3D link)",
        "class_color": "#F8B500",
        "panel_letter": "4f",
        "out_pdf": OUT_DIR / "fig4_f_F15_NTRK2.pdf",
        "interpretation_paragraphs": [
            ("Data:",
             "F15 (NTRK2 / TrkB signalling, RA-biased baseline) sits just below the Strong-trend threshold "
             "(\u03c1 = -0.39 in primary, 4/4 sensitivity concordance)."),
            ("Cross-figure link:",
             "F15 was classified as cellular relay in Fig 3D (significant downshift in lymphocyte + endothelial, "
             "carrier fibroblast unchanged). The convergence of single-cell relay (Fig 3D) and bulk duration "
             "decline (Fig 4f) strengthens F15's atlas prioritisation despite being individually marginal in "
             "either axis alone."),
        ],
    },
}


def _wrap_paragraph(text: str, max_chars: int = 140) -> list[str]:
    """Soft-wrap a paragraph at word boundaries to ~max_chars per line."""
    words = text.split()
    lines: list[str] = []
    cur = ""
    for w in words:
        if not cur:
            cur = w
        elif len(cur) + 1 + len(w) <= max_chars:
            cur = cur + " " + w
        else:
            lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return lines


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": enforce_font_floor(FS_TICK),
        "axes.linewidth": 0.6,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def render_case(scores_long, meta, cls_df, key, cfg):
    fid = cfg["family_id"]
    sub = scores_long[(scores_long["family_id"] == fid) & (scores_long["tissue"] == "LAA")].copy()
    if sub.empty:
        print(f"[SKIP] {key}: no LAA family score")
        return
    # family-level: average across programs in family
    fam_scores = sub.groupby("gsm", as_index=False)["score"].mean()
    merged = fam_scores.merge(meta[["gsm", "donor_id", "rhythm", "duration_months"]], on="gsm", how="left")
    merged["duration_months"] = pd.to_numeric(merged["duration_months"], errors="coerce")

    af = merged[(merged["rhythm"] == "AF") & merged["duration_months"].notna()].copy()
    sr = merged[merged["rhythm"] == "SR"].copy()
    if len(af) < 4:
        print(f"[SKIP] {key}: insufficient AF samples")
        return

    cls_row = cls_df.set_index("family_id").loc[fid]
    rho = float(cls_row["rho_primary"])
    p_perm = float(cls_row["p_perm_primary"])
    ci_lo = float(cls_row["ci_lo_primary"])
    ci_hi = float(cls_row["ci_hi_primary"])
    sens = {
        "primary": float(cls_row["rho_primary"]),
        "gt1m_eq1": float(cls_row["rho_gt1m_eq1"]),
        "ordinal": float(cls_row["rho_ordinal"]),
        "lapv": float(cls_row["rho_lapv"]),
    }
    sens_n = int(cls_row["sens_same_sign_n"])

    MM = 1.0 / 25.4
    fig = plt.figure(figsize=(170 * MM, 150 * MM))  # taller for larger fonts
    gs = fig.add_gridspec(
        1, 1, left=0.13, right=0.96, top=0.76, bottom=0.40,
    )
    ax = fig.add_subplot(gs[0, 0])

    # SR anchor (visual only; not in fit)
    sr_handle = None
    if not sr.empty:
        sr_handle = ax.scatter(sr["duration_months"], sr["score"], s=42, marker="s",
                                 facecolor="white", edgecolor=SR_COLOR, linewidth=1.2, zorder=2)
    # AF samples used in fit
    af_handle = ax.scatter(af["duration_months"], af["score"], s=58, marker="o",
                             facecolor=AF_COLOR, edgecolor="#1A1A1A", linewidth=0.6,
                             alpha=0.85, zorder=3)

    # Linear fit on AF only (visualisation; ρ is rank-based)
    if len(af) >= 3:
        x = af["duration_months"].to_numpy(float)
        y = af["score"].to_numpy(float)
        slope, intercept = np.polyfit(x, y, 1)
        x_line = np.linspace(x.min(), x.max(), 50)
        ax.plot(x_line, slope * x_line + intercept, color="#1A1A1A", lw=1.2,
                ls="--", alpha=0.7, zorder=2)

    ax.set_xlabel("AF duration (months, LAA)", fontsize=enforce_font_floor(FS_AXIS_LABEL))
    ax.set_ylabel(f"Family-level singscore  (F{fid})", fontsize=enforce_font_floor(FS_AXIS_LABEL))
    ax.tick_params(labelsize=enforce_font_floor(FS_TICK))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(0.4)
    ax.spines["bottom"].set_linewidth(0.4)

    legend_handles = [af_handle]
    legend_labels = [f"AF donors (n = {len(af)}, finite duration; fit + ρ)"]
    if sr_handle is not None:
        legend_handles.append(sr_handle)
        legend_labels.append(f"SR controls (blue open square, n = {len(sr)}; visual anchor only)")
    fig.legend(
        legend_handles, legend_labels,
        loc="upper left", bbox_to_anchor=(0.13, 0.86),
        fontsize=enforce_font_floor(FS_TICK - 0.5), frameon=False, ncol=2, handlelength=1.1, columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.text(
        0.13, 0.795,
        "AF donors with non-numeric duration (GEO '>1M') are not plotted; ρ, CI, p_perm, chronicity class: "
        f"n = {len(af)} AF LAA donors with finite duration.",
        fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5),
        color="#444444",
        ha="left",
        va="top",
    )

    ci_crosses_zero = (ci_lo <= 0.0 <= ci_hi)
    ci_line = f"95% CI: [{ci_lo:+.2f}, {ci_hi:+.2f}]"
    if ci_crosses_zero:
        ci_line = ci_line + "  -- includes 0"
    p_marginal = (0.05 < p_perm < 0.10)
    p_line = f"p_perm = {p_perm:.3f}"
    if p_marginal:
        p_line = p_line + "  \u2020  marginal"
    stat_lines = [
        f"\u03c1 = {rho:+.2f}",
        ci_line,
        p_line,
        f"n = {len(af)} AF donors",
    ]
    box_text = "\n".join(stat_lines)
    ax.text(
        0.97, 0.97, box_text,
        transform=ax.transAxes,
        fontsize=enforce_font_floor(FS_TICK), color="#1A1A1A", ha="right", va="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#888", linewidth=0.5),
    )
    note_lines = ["Bootstrap 5000 resamples;  permutation 10,000 shuffles."]
    if ci_crosses_zero:
        note_lines.append("CI crosses 0 at n = 12 -- reported transparently.")
    if p_marginal:
        note_lines.append("\u2020 0.05 < p_perm < 0.10:  treat as a directional trend, not a strong test.")
    ax.text(
        0.97, 0.69, "\n".join(note_lines),
        transform=ax.transAxes,
        fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5),
        color="#C44569" if p_marginal else "#666666",
        ha="right", va="top", style="italic",
    )

    coding_label = {
        "primary":  "primary  (LAA, n = 12)",
        "gt1m_eq1": "gt1m_eq1 (LAA, n = 12)",
        "ordinal":  "ordinal  (LAA, n = 12)",
        "lapv":     "LAPV     (region, n = 12)",
    }
    body_rows = []
    reverser_keys = []
    for k, v in sens.items():
        same_sign = (np.sign(v) == np.sign(rho) and rho != 0)
        marker = "[same]" if same_sign else "[opp]"
        if not same_sign and k != "primary":
            reverser_keys.append(k)
        body_rows.append((coding_label[k], f"\u03c1 = {v:+.2f}", marker,
                            k == "primary", same_sign))

    x_left = 0.36
    y_top = 0.30
    line_h = 0.045
    ax.text(x_left, y_top, "Sensitivity codings",
              transform=ax.transAxes,
              fontsize=enforce_font_floor(FS_TICK), fontweight="bold", color="#1A1A1A", ha="left", va="top")
    for i, (lab, rho_str, marker, is_primary, same_sign) in enumerate(body_rows):
        y = y_top - line_h * (i + 1)
        if is_primary:
            weight, color, style = "bold", "#1A1A1A", "normal"
        elif not same_sign:
            weight, color, style = "normal", "#888888", "italic"
        else:
            weight, color, style = "normal", "#333333", "normal"
        ax.text(x_left + 0.005, y, lab,
                  transform=ax.transAxes,
                  fontsize=enforce_font_floor(FS_TICK - 0.5), color=color, fontweight=weight, fontstyle=style,
                  ha="left", va="top", family="monospace")
        ax.text(x_left + 0.42, y, rho_str,
                  transform=ax.transAxes,
                  fontsize=enforce_font_floor(FS_TICK - 0.5), color=color, fontweight=weight, fontstyle=style,
                  ha="right", va="top", family="monospace")
        ax.text(x_left + 0.45, y, marker,
                  transform=ax.transAxes,
                  fontsize=enforce_font_floor(FS_TICK - 0.5), color=color, fontstyle=style,
                  ha="left", va="top")
    y_conc = y_top - line_h * (len(body_rows) + 1)
    conc_text = f"Same-sign concordance:  {sens_n} / 4"
    if reverser_keys:
        rev_str = ", ".join(k.upper() for k in reverser_keys)
        conc_text = conc_text + f"   ({rev_str} reverses sign)"
    ax.text(x_left, y_conc, conc_text,
              transform=ax.transAxes,
              fontsize=enforce_font_floor(FS_TICK - 0.5), color="#1A1A1A", fontweight="bold", ha="left", va="top")
    box_w, box_h = 0.58, line_h * (len(body_rows) + 2) + 0.02
    rect = mpatches.FancyBboxPatch(
        (x_left - 0.02, y_conc - 0.025), box_w, box_h,
        boxstyle="round,pad=0.005,rounding_size=0.01",
        facecolor="#F5F5F5", edgecolor="#888", linewidth=0.4,
        transform=ax.transAxes, zorder=0.5,
    )
    ax.add_patch(rect)

    fig.text(
        0.05, 0.92,
        f"Fig {cfg['panel_letter']}  {cfg['label']}  \u2014  {cfg['class_name']}",
        fontsize=enforce_font_floor(FS_PANEL_TITLE + 2), fontweight="bold", color="#1A1A1A",
    )
    fig.add_artist(
        mpatches.FancyBboxPatch(
            (0.74, 0.905), 0.22, 0.040,
            boxstyle="round,pad=0.01,rounding_size=0.012",
            facecolor=cfg["class_color"], edgecolor="none",
            transform=fig.transFigure,
        )
    )
    fig.text(0.85, 0.925, cfg["class_name"], fontsize=enforce_font_floor(FS_TICK), fontweight="bold",
             color="white", ha="center", va="center")

    y_cur = 0.34
    line_h = 0.022
    header_paragraph = (
        "GSE41177 LAA bulk Affymetrix expression.  "
        f"Panels show SR controls at month 0 (n = {len(sr)}, anchor only; excluded from correlation) plus "
        f"AF donors with finite duration (n = {len(af)}).  "
        "Additional AF LAA donors with censored duration strings ('>1M' in GEO metadata) are retained in the cohort "
        f"but omitted from scatter panels; reported Spearman ρ, bootstrap CI, permutation p, and chronicity class "
        f"use only those n = {len(af)} AF donors.  Linear trend line is visualisation only (rank-based ρ).  "
        "LAPV sensitivity row uses paired LAPV cohort (same donors), not LAA-pooled."
    )
    for j, ln in enumerate(_wrap_paragraph(header_paragraph, max_chars=130)):
        fig.text(0.05, y_cur - line_h * j, ln, fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5), color="#444444")
    y_cur -= line_h * (len(_wrap_paragraph(header_paragraph, max_chars=130)) + 0.6)
    for tag, body in cfg["interpretation_paragraphs"]:
        # Tag on its own line (avoids collision with wrapped body)
        fig.text(0.05, y_cur, tag, fontsize=enforce_font_floor(FS_FOOTNOTE), fontweight="bold", color="#1A1A1A")
        y_cur -= line_h * 0.95
        wrapped = _wrap_paragraph(body, max_chars=130)
        for j, ln in enumerate(wrapped):
            fig.text(0.07, y_cur - line_h * j, ln, fontsize=enforce_font_floor(FS_FOOTNOTE - 0.5), color="#444444")
        y_cur -= line_h * (len(wrapped) + 0.4)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(cfg["out_pdf"], bbox_inches=None, pad_inches=0.04)
    plt.close(fig)
    print("wrote " + str(cfg["out_pdf"].relative_to(ROOT)))


def main():
    set_style()
    scores_long = pd.read_csv(LONG_TSV, sep="\t")
    meta = pd.read_csv(META_TSV, sep="\t")
    cls_df = pd.read_csv(CLS_TSV, sep="\t")
    for key, cfg in CASES.items():
        render_case(scores_long, meta, cls_df, key, cfg)


if __name__ == "__main__":
    main()
