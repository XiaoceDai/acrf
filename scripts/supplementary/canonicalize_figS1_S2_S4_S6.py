#!/usr/bin/env python3
"""
Canonicalize S1, S2, S4, S6 into figures/supplementary/ with publication footers.

Uses existing renderers with output paths redirected so legacy assets and
authoritative source tables under results/ are not overwritten.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "supplementary"))

from _export import PREVIEW, SUPP, save_fig_bundle, set_pub_style  # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _savefig_to_bundle(fig: plt.Figure, stem: str, footer: str | None = None) -> dict:
    if footer:
        fig.text(
            0.5, 0.01, footer,
            ha="center", va="bottom", fontsize=5.2, color="#333333",
            wrap=True, transform=fig.transFigure,
        )
    return save_fig_bundle(fig, stem, pad_inches=0.04)


def _hooked_bundle(stem: str, title: str | None, footer: str, captured: dict):
    """Return a Figure.savefig replacement that writes a canonical bundle once."""
    orig = plt.Figure.savefig

    def hooked(self, fname, *args, **kwargs):  # type: ignore[no-untyped-def]
        # Prevent recursion: restore original before bundle export
        plt.Figure.savefig = orig  # type: ignore[assignment]
        if title:
            if getattr(self, "_suptitle", None) is None:
                self.suptitle(title, fontsize=8.5, fontweight="bold", y=0.995)
            else:
                try:
                    self._suptitle.set_text(title)
                except Exception:
                    self.suptitle(title, fontsize=8.5, fontweight="bold", y=0.995)
        if footer:
            self.text(
                0.5, 0.01, footer,
                ha="center", va="bottom", fontsize=5.0, color="#333333",
                transform=self.transFigure,
            )
        captured["paths"] = save_fig_bundle(self, stem, pad_inches=0.04)
        return None

    return hooked


def build_s1() -> dict:
    set_pub_style()
    mod = _load("figS1_cohort_overview", ROOT / "scripts" / "figS1_cohort_overview.py")
    tmp_sum = ROOT / "results" / "manuscript_audit" / "_phase5b3a_tmp_figS1_summary.tsv"
    mod.OUT = SUPP / "figS1_cohort_overview.pdf"
    mod.OUT_SUMMARY = tmp_sum
    captured: dict = {}
    orig = plt.Figure.savefig
    plt.Figure.savefig = _hooked_bundle(  # type: ignore[assignment]
        "figS1_cohort_overview",
        "Cohort overview and analysis roles of public datasets",
        "Counts use the analysis unit appropriate to each dataset and are not a pooled donor total. "
        "Publication-facing term: NCP.",
        captured,
    )
    try:
        mod.main()
    finally:
        plt.Figure.savefig = orig  # type: ignore[assignment]
    if tmp_sum.exists():
        tmp_sum.unlink()
    return {
        "paths": captured["paths"],
        "panel_structure": "a–d cohort overview",
        "notes": "rerendered via figS1 with redirected outputs + analysis-unit footer",
    }


def build_s2() -> dict:
    set_pub_style()
    mod = _load("figS2_pairing", ROOT / "scripts" / "figS2_pairing_attribute_heatmap.py")
    tmp_sum = ROOT / "results" / "manuscript_audit" / "_phase5b3a_tmp_figS2_summary.tsv"
    mod.OUT = SUPP / "figS2_pairing_sample_attributes.pdf"
    mod.SUMM = tmp_sum
    # Resolve relocated metadata without editing the legacy renderer file
    path_map = {
        "GSE128188": ROOT / "Project_P2" / "Project_P1" / "GSE128188" / "sample_metadata.tsv",
        "GSE79768": ROOT / "Project_P2" / "Project_P1" / "GSE79768" / "sample_metadata.tsv",
        "GSE115574": ROOT / "Project_P2" / "Project_P1" / "GSE115574" / "sample_metadata.tsv",
        "GSE41177": ROOT / "Project_P4" / "GSE41177" / "sample_metadata.tsv",
    }
    for cfg in mod.COHORTS:
        p = path_map.get(cfg["name"])
        if p is not None and p.exists():
            cfg["path"] = p
        elif not Path(cfg["path"]).exists():
            raise FileNotFoundError(f"S2 metadata missing for {cfg['name']}: {cfg['path']}")
    captured: dict = {}
    orig = plt.Figure.savefig
    plt.Figure.savefig = _hooked_bundle(  # type: ignore[assignment]
        "figS2_pairing_sample_attributes",
        "Donor pairing and sample-attribute structure",
        "GSE41177 SR samples document cohort structure only; they were not used as duration anchors "
        "in the final Figure 4 inference or final Figure 4D panels. Datasets are not a pooled donor denominator.",
        captured,
    )
    try:
        mod.main()
    finally:
        plt.Figure.savefig = orig  # type: ignore[assignment]
    if tmp_sum.exists():
        tmp_sum.unlink()
    return {
        "paths": captured["paths"],
        "panel_structure": "a–d pairing matrices",
        "notes": "SR duration-anchor boundary footer added",
    }


def build_s4() -> dict:
    """Replot S4 from frozen results/figS4_per_cohort_pca.tsv — no PCA recomputation."""
    import numpy as np
    import pandas as pd
    from matplotlib.lines import Line2D

    set_pub_style()
    df = pd.read_csv(ROOT / "results" / "figS4_per_cohort_pca.tsv", sep="\t")
    AF_COLOR, SR_COLOR = "#B2182B", "#2166AC"
    HEALTHY_LA, HEALTHY_RA = "#9C5550", "#4E6B89"

    MM = 1 / 25.4
    fig = plt.figure(figsize=(180 * MM, 120 * MM))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.45, left=0.08, right=0.98, top=0.88, bottom=0.12)
    panels = [
        ("a", "GSE112339"),
        ("b", "GSE128188"),
        ("c", "GSE79768"),
        ("d", "GSE115574"),
        ("e", "GSE41177"),
    ]
    coords = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
    for (lab, cohort), (r, c) in zip(panels, coords):
        ax = fig.add_subplot(gs[r, c])
        sub = df[df.cohort == cohort]
        for _, row in sub.iterrows():
            lab_s = str(row.label)
            if row.category == "healthy_chamber":
                col = HEALTHY_LA if lab_s == "LA" else HEALTHY_RA
                marker = "o"
            elif lab_s in {"AF", "SR"}:
                col = AF_COLOR if lab_s == "AF" else SR_COLOR
                marker = "o"
            else:
                is_af = lab_s.startswith("AF")
                is_ra = lab_s.endswith("RA") or lab_s.endswith("LAPV")
                col = AF_COLOR if is_af else SR_COLOR
                marker = "^" if is_ra else "o"
            ax.scatter(
                row.PC1, row.PC2, c=col, marker=marker, s=12,
                edgecolors="white", linewidths=0.3, zorder=3,
            )
        ax.set_title(f"{lab}  {cohort}", loc="left", fontsize=7.0, fontweight="bold")
        ax.set_xlabel("PC1", fontsize=6)
        ax.set_ylabel("PC2", fontsize=6)
        ax.tick_params(labelsize=5.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    leg = fig.add_subplot(gs[1, 2])
    leg.axis("off")
    handles = [
        Line2D([], [], marker="o", color="w", markerfacecolor=HEALTHY_LA, markersize=6, label="healthy LA"),
        Line2D([], [], marker="o", color="w", markerfacecolor=HEALTHY_RA, markersize=6, label="healthy RA"),
        Line2D([], [], marker="o", color="w", markerfacecolor=SR_COLOR, markersize=6, label="SR — LA/LAA"),
        Line2D([], [], marker="^", color="w", markerfacecolor=SR_COLOR, markersize=6, label="SR — RA/LAPV"),
        Line2D([], [], marker="o", color="w", markerfacecolor=AF_COLOR, markersize=6, label="AF — LA/LAA"),
        Line2D([], [], marker="^", color="w", markerfacecolor=AF_COLOR, markersize=6, label="AF — RA/LAPV"),
    ]
    leg.legend(handles=handles, loc="center", frameon=False, fontsize=6)
    fig.suptitle(
        "Per-cohort principal-component overview of expression datasets",
        x=0.08, y=0.96, ha="left", fontsize=9.0, fontweight="bold",
    )
    fig.text(
        0.5, 0.02,
        "QC and sample-context visualization only; not an independent AF-remodeling test. "
        "Coordinates from frozen results/figS4_per_cohort_pca.tsv (no PCA recomputation).",
        ha="center", va="bottom", fontsize=5.2, color="#333333",
    )
    paths = save_fig_bundle(fig, "figS4_cohort_pca_overview", pad_inches=0.03)
    return {
        "paths": paths,
        "panel_structure": "a–e per-cohort PCA",
        "notes": "Replotted from frozen figS4_per_cohort_pca.tsv; no PCA recomputation",
    }


def build_s6() -> dict:
    """Promote chamber-contribution asset with optional title confirmation via copy+export."""
    from _export import export_from_pdf

    src = ROOT / "figures" / "ncp_v3_singscore_reanalysis" / "supp_chamber_contribution_to_af_ncp_family_rewiring.pdf"
    # Prefer existing SVG/PNG if present, but force canonical names + 600 dpi PNG
    paths = export_from_pdf(src, "figS6_chamber_contribution_af_ncp_rewiring")
    # Also copy native SVG if available (higher fidelity than pdf→svg)
    native_svg = src.with_suffix(".svg")
    if native_svg.exists():
        paths["svg"].write_bytes(native_svg.read_bytes())
    return {
        "paths": paths,
        "panel_structure": "three-cohort chamber contribution bubbles",
        "notes": "Promoted QC-PASS asset; no scientific recomputation",
    }


if __name__ == "__main__":
    for fn in (build_s1, build_s2, build_s4, build_s6):
        info = fn()
        print(fn.__name__, info["paths"]["pdf"])
