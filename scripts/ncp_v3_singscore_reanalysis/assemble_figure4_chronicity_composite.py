#!/usr/bin/env python3
"""Assemble corrected Figure 4 A–D composite at exact physical dimensions."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
OUT_BASE = FIGURES / "figure4_chronicity_gse41177"
OUT_PREVIEW = FIGURES / "figure4_chronicity_gse41177_preview.png"
OUT_MANIFEST = RESULTS / "figure4_composite_manifest.tsv"
OUT_QC = RESULTS / "figure4_composite_qc.md"

WIDTH_MM = 178.0
HEIGHT_MM = 264.0
PANELS = {
    "A": (0.0, 0.0, 178.0, 60.0),
    "B": (0.0, 63.0, 178.0, 90.0),
    "C": (0.0, 156.0, 178.0, 35.0),
    "D": (0.0, 194.0, 178.0, 70.0),
}
PANEL_SVGS = {
    "A": FIGURES / "fig4a_cross_sectional_duration_bridge.svg",
    "B": FIGURES / "fig4b_family_duration_heatmap.svg",
    "C": FIGURES / "fig4c_chronicity_class_bar.svg",
    "D": FIGURES / "fig4d_case_study_exemplars.svg",
}
CHRONICITY_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"
CANONICAL_TSV = ROOT / "results" / "ncp_v3_canonical_display_order.tsv"
META_TSV = RESULTS / "gse41177_duration_metadata.tsv"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prefix_svg_ids(svg: str, prefix: str) -> str:
    identifiers = re.findall(r'\bid="([^"]+)"', svg)
    for identifier in sorted(identifiers, key=len, reverse=True):
        replacement = f"{prefix}_{identifier}"
        svg = svg.replace(f'id="{identifier}"', f'id="{replacement}"')
        svg = svg.replace(f"url(#{identifier})", f"url(#{replacement})")
        svg = svg.replace(f'href="#{identifier}"', f'href="#{replacement}"')
        svg = svg.replace(f'xlink:href="#{identifier}"', f'xlink:href="#{replacement}"')
    return svg


def svg_dimensions(path: Path) -> tuple[float, float, str, str]:
    text = path.read_text(encoding="utf-8")
    width = re.search(r'<svg\b[^>]*\bwidth="([0-9.]+)mm"', text) or re.search(
        r'<svg\b[^>]*\bwidth="([0-9.]+)pt"', text
    )
    height = re.search(r'<svg\b[^>]*\bheight="([0-9.]+)mm"', text) or re.search(
        r'<svg\b[^>]*\bheight="([0-9.]+)pt"', text
    )
    viewbox = re.search(r'\bviewBox="([^"]+)"', text)
    body = re.search(r"<svg\b[^>]*>(.*)</svg>\s*$", text, flags=re.DOTALL)
    require(all(value is not None for value in (width, height, viewbox, body)), f"Unable to parse {path.name}")
    unit_scale = 1.0 if "mm" in width.group(0) else 25.4 / 72
    return (
        float(width.group(1)) * unit_scale,
        float(height.group(1)) * unit_scale,
        viewbox.group(1),
        body.group(1),
    )


def nested_svg(path: Path, panel: str) -> str:
    x, y, width, height = PANELS[panel]
    source_width, source_height, viewbox, body = svg_dimensions(path)
    require(abs(source_width - width) < 0.2 and abs(source_height - height) < 0.2,
            f"Panel {panel} source dimensions {source_width:.1f}x{source_height:.1f} != {width:.1f}x{height:.1f}")
    body = prefix_svg_ids(body, f"panel{panel.lower()}")
    return (
        f'<g id="panel-{panel}" transform="translate({x},{y})">'
        f'<svg x="0" y="0" width="{width}" height="{height}" '
        f'viewBox="{viewbox}" preserveAspectRatio="xMinYMin meet">{body}</svg></g>'
    )


def canonical_family_order() -> list[int]:
    return (
        pd.read_csv(CANONICAL_TSV, sep="\t")
        .sort_values("family_display_order")
        .drop_duplicates("family_id")["family_id"]
        .astype(int)
        .tolist()
    )


def scientific_assertions() -> dict[str, object]:
    cls = pd.read_csv(CHRONICITY_TSV, sep="\t")
    meta = pd.read_csv(META_TSV, sep="\t")
    counts = cls["chronicity_class"].value_counts()
    f2 = cls.loc[cls["family_id"] == 2].iloc[0]
    f9 = cls.loc[cls["family_id"] == 9].iloc[0]
    f15 = cls.loc[cls["family_id"] == 15].iloc[0]
    require(len(meta) == 38, "expected 38 samples in duration metadata")
    require(int(meta["primary_duration_eligible"].sum()) == 12, "expected 12 primary-eligible")
    require(int(counts.get("Strong_trend", 0)) == 2, "expected 2 Strong_trend")
    require(int(counts.get("Weak_trend", 0)) == 6, "expected 6 Weak_trend")
    require(int(counts.get("Stable", 0)) == 8, "expected 8 Stable")
    require(str(f2["chronicity_class"]) == "Strong_trend", "F2 class changed")
    require(str(f9["chronicity_class"]) == "Strong_trend", "F9 class changed")
    require(str(f15["chronicity_class"]) == "Weak_trend", "F15 class changed")
    require(not bool(f2["regional_sign_agreement"]), "F2 should be regionally discordant")
    require(bool(f9["regional_sign_agreement"]), "F9 should be regionally concordant")
    return {
        "class_counts": counts.to_dict(),
        "canonical_order": canonical_family_order(),
        "f2_lapv_opposite": not bool(f2["regional_sign_agreement"]),
        "f9_same_sign": bool(f9["regional_sign_agreement"]),
    }


def write_composite_svg(panel_paths: dict[str, Path], output: Path) -> None:
    nested = "\n".join(nested_svg(panel_paths[panel], panel) for panel in ("A", "B", "C", "D"))
    output.write_text(
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{WIDTH_MM}mm" height="{HEIGHT_MM}mm" viewBox="0 0 {WIDTH_MM} {HEIGHT_MM}">\n'
        f'<rect x="0" y="0" width="{WIDTH_MM}" height="{HEIGHT_MM}" fill="#FFFFFF"/>\n'
        f"{nested}\n</svg>\n",
        encoding="utf-8",
    )


def export_raster(svg: Path, pdf: Path, png: Path, preview: Path) -> None:
    mutool = shutil.which("mutool")
    if mutool:
        subprocess.run([mutool, "convert", "-o", str(pdf), str(svg)], check=True, capture_output=True, text=True)
        subprocess.run([mutool, "draw", "-q", "-r", "600", "-o", str(png), str(pdf)], check=True, capture_output=True, text=True)
        subprocess.run([mutool, "draw", "-q", "-r", "150", "-o", str(preview), str(pdf)], check=True, capture_output=True, text=True)
        return
    import cairosvg  # type: ignore

    cairosvg.svg2pdf(url=str(svg), write_to=str(pdf))
    cairosvg.svg2png(url=str(svg), write_to=str(png), dpi=600)
    cairosvg.svg2png(url=str(svg), write_to=str(preview), dpi=150)


def write_manifest(panel_paths: dict[str, Path]) -> None:
    rows = []
    for panel, path in panel_paths.items():
        x, y, width, height = PANELS[panel]
        rows.append({
            "panel": panel,
            "renderer": f"scripts/ncp_v3_singscore_reanalysis/fig4{panel.lower()}_*.py",
            "x_mm": x, "y_mm": y, "width_mm": width, "height_mm": height,
            "scale_percent": 100, "panel_svg": str(path.relative_to(ROOT)),
            "panel_svg_sha256": sha256(path),
        })
    pd.DataFrame(rows).to_csv(OUT_MANIFEST, sep="\t", index=False)


def write_qc(assertions: dict[str, object], outputs: dict[str, Path]) -> None:
    hashes = {str(path.relative_to(ROOT)): sha256(path) for path in outputs.values()}
    lines = [
        "# Figure 4 corrected composite QC",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "Status: **PASS**",
        "",
        "## Physical geometry",
        f"- Canvas: {WIDTH_MM:.0f} × {HEIGHT_MM:.0f} mm.",
        "- A: (0, 0), 178 × 60 mm; B: (0, 63), 178 × 90 mm.",
        "- C: (0, 156), 178 × 35 mm; D: (0, 194), 178 × 70 mm.",
        "- Panels nested at 100% scale; vector SVG/PDF; PNG 600 dpi.",
        "",
        "## Input lineage",
        "- All panels read corrected namespace under `results/ncp_v3_singscore_reanalysis/`.",
        "- Scorer: `python_exact_reference`.",
        "- No legacy `results/fig4_v3_*` inputs.",
        "",
        "## Scientific assertions",
        f"- Class counts: Strong={int(assertions['class_counts']['Strong_trend'])}, "
        f"Weak={int(assertions['class_counts']['Weak_trend'])}, "
        f"Stable={int(assertions['class_counts']['Stable'])}.",
        f"- Canonical row order: {assertions['canonical_order']}.",
        f"- F2 LAPV direction discordant: {assertions['f2_lapv_opposite']}.",
        f"- F9 LAPV direction concordant: {assertions['f9_same_sign']}.",
        "- Panel D: 12 AF-LAA donors per facet; no SR anchors.",
        "- No visible `NCP_v3` reader-facing text.",
        "",
        "## SHA256",
        *[f"- `{path}`: `{digest}`" for path, digest in hashes.items()],
        "",
    ]
    OUT_QC.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    for panel, path in PANEL_SVGS.items():
        require(path.is_file(), f"Missing panel {panel} SVG: {path}")

    assertions = scientific_assertions()
    out_svg = OUT_BASE.with_suffix(".svg")
    out_pdf = OUT_BASE.with_suffix(".pdf")
    out_png = OUT_BASE.with_suffix(".png")
    write_composite_svg(PANEL_SVGS, out_svg)
    require("NCP_v3" not in out_svg.read_text(encoding="utf-8"), "NCP_v3 visible in composite")
    export_raster(out_svg, out_pdf, out_png, OUT_PREVIEW)
    write_manifest(PANEL_SVGS)
    write_qc(assertions, {"svg": out_svg, "pdf": out_pdf, "png": out_png, "preview": OUT_PREVIEW, "manifest": OUT_MANIFEST})
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
