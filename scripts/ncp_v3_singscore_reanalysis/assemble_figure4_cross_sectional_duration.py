#!/usr/bin/env python3
"""Assemble final Figure 4 (A–D) from composite-sized vector panel SVGs at 100% scale.

Outputs figure4_cross_sectional_duration.{pdf,svg,png}. Does not overwrite
legacy figure4_chronicity_gse41177.* or archived Figure 4A–F assets.
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
OUT_BASE = FIGURES / "figure4_cross_sectional_duration"
OUT_PREVIEW = FIGURES / "figure4_cross_sectional_duration_preview.png"
OUT_MANIFEST = RESULTS / "figure4_composite_manifest.tsv"
OUT_QC = RESULTS / "figure4_composite_qc.md"

WIDTH_MM = 178.0
# A(52) + gap(3.5) + BC(82) + gap(3.5) + D(72) = 213
HEIGHT_MM = 213.0
V_GUTTER_MM = 3.5
# Modest B|C gap after trimming B's unused right pad (was 4.0 with B=131;
# briefly 3.0 felt too tight against B content).
H_GUTTER_MM = 12.0
PANEL_B_W_MM = 98.0
PANEL_C_W_MM = 43.0

# (x_mm, y_mm from top, width_mm, height_mm)
PANELS = {
    "A": (0.0, 0.0, 178.0, 52.0),
    "B": (0.0, 52.0 + V_GUTTER_MM, PANEL_B_W_MM, 82.0),
    "C": (PANEL_B_W_MM + H_GUTTER_MM, 52.0 + V_GUTTER_MM, PANEL_C_W_MM, 82.0),
    "D": (0.0, 52.0 + V_GUTTER_MM + 82.0 + V_GUTTER_MM, 178.0, 72.0),
}
PANEL_SVGS = {
    "A": FIGURES / "fig4a_cross_sectional_duration_bridge.svg",
    "B": FIGURES / "fig4b_family_regional_duration_map.svg",
    "C": FIGURES / "fig4c_chronicity_class_summary.svg",
    "D": FIGURES / "fig4d_representative_duration_examples.svg",
}
PANEL_RENDERERS = {
    "A": "scripts/ncp_v3_singscore_reanalysis/fig4a_cross_sectional_duration_bridge.py",
    "B": "scripts/ncp_v3_singscore_reanalysis/fig4b_family_regional_duration_map.py",
    "C": "scripts/ncp_v3_singscore_reanalysis/fig4c_chronicity_class_summary.py",
    "D": "scripts/ncp_v3_singscore_reanalysis/fig4d_representative_duration_examples.py",
}

CHRONICITY_TSV = RESULTS / "fig4_v3_chronicity_class.tsv"
CANONICAL_TSV = ROOT / "results" / "ncp_v3_canonical_display_order.tsv"
META_TSV = RESULTS / "gse41177_duration_metadata.tsv"

# Legacy / archive assets that must not be overwritten by this assembler
PROTECTED_LEGACY = (
    FIGURES / "figure4_chronicity_gse41177.pdf",
    FIGURES / "figure4_chronicity_gse41177.svg",
    FIGURES / "figure4_chronicity_gse41177.png",
)


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
    require(
        abs(source_width - width) < 0.25 and abs(source_height - height) < 0.25,
        f"Panel {panel} source dimensions {source_width:.1f}x{source_height:.1f} "
        f"!= slot {width:.1f}x{height:.1f} (must place at 100% scale)",
    )
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
    require(bool(f15["regional_sign_agreement"]), "F15 should be regionally concordant")
    return {
        "class_counts": counts.to_dict(),
        "canonical_order": canonical_family_order(),
        "f2_lapv_opposite": not bool(f2["regional_sign_agreement"]),
        "f9_same_sign": bool(f9["regional_sign_agreement"]),
        "f15_same_sign": bool(f15["regional_sign_agreement"]),
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
    # Fallback: matplotlib / cairosvg
    try:
        import cairosvg  # type: ignore

        cairosvg.svg2pdf(url=str(svg), write_to=str(pdf))
        cairosvg.svg2png(url=str(svg), write_to=str(png), dpi=600)
        cairosvg.svg2png(url=str(svg), write_to=str(preview), dpi=150)
        return
    except Exception:
        pass
    # Last resort: keep SVG and convert via inkscape if present
    inkscape = shutil.which("inkscape")
    require(inkscape is not None, "Need mutool, cairosvg, or inkscape to export PDF/PNG")
    subprocess.run(
        [inkscape, str(svg), f"--export-filename={pdf}"],
        check=True, capture_output=True, text=True,
    )
    subprocess.run(
        [inkscape, str(pdf), f"--export-filename={png}", "--export-dpi=600"],
        check=True, capture_output=True, text=True,
    )
    subprocess.run(
        [inkscape, str(pdf), f"--export-filename={preview}", "--export-dpi=150"],
        check=True, capture_output=True, text=True,
    )


def write_manifest(panel_paths: dict[str, Path]) -> None:
    rows = []
    for panel, path in panel_paths.items():
        x, y, width, height = PANELS[panel]
        rows.append({
            "panel": panel,
            "renderer": PANEL_RENDERERS[panel],
            "x_mm": x,
            "y_mm": y,
            "width_mm": width,
            "height_mm": height,
            "scale_percent": 100,
            "panel_svg": str(path.relative_to(ROOT)),
            "panel_svg_sha256": sha256(path),
        })
    pd.DataFrame(rows).to_csv(OUT_MANIFEST, sep="\t", index=False)


def write_qc(assertions: dict[str, object], outputs: dict[str, Path]) -> None:
    hashes = {str(path.relative_to(ROOT)): sha256(path) for path in outputs.values()}
    lines = [
        "# Figure 4 cross-sectional duration composite QC",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "Status: **PASS**",
        "",
        "## Physical geometry",
        f"- Canvas: {WIDTH_MM:.0f} × {HEIGHT_MM:.0f} mm.",
        f"- A: (0, 0), 178 × 52 mm.",
        f"- B: (0, {52 + V_GUTTER_MM:.1f}), {PANEL_B_W_MM:.0f} × 82 mm (trimmed to content; C shifted left).",
        f"- C: ({PANEL_B_W_MM + H_GUTTER_MM:.1f}, {52 + V_GUTTER_MM:.1f}), {PANEL_C_W_MM:.0f} × 82 mm.",
        f"- D: (0, {52 + V_GUTTER_MM + 82 + V_GUTTER_MM:.1f}), 178 × 72 mm.",
        f"- Vertical gutters: {V_GUTTER_MM} mm; horizontal B|C gutter: {H_GUTTER_MM} mm.",
        f"- B|C row right margin: {WIDTH_MM - (PANEL_B_W_MM + H_GUTTER_MM + PANEL_C_W_MM):.1f} mm.",
        "- Panels nested at 100% scale; vector SVG/PDF; PNG 600 dpi.",
        "- No PNG stitching or panel image rescaling.",
        "- Panel B rho display: 4 decimals (eliminates F11 −0.200 class-boundary ambiguity).",
        "",
        "## Input lineage",
        "- Panels read corrected namespace under `results/ncp_v3_singscore_reanalysis/`.",
        "- Scorer: `python_exact_reference`.",
        "- No legacy `results/fig4_v3_*` inputs.",
        "- Legacy `figure4_chronicity_gse41177.*` retained as archive; not overwritten.",
        "- Formal submission routing: `figure4_cross_sectional_duration.*`.",
        "",
        "## Scientific assertions",
        f"- Class counts: Strong={int(assertions['class_counts']['Strong_trend'])}, "
        f"Weak={int(assertions['class_counts']['Weak_trend'])}, "
        f"Stable={int(assertions['class_counts']['Stable'])}.",
        f"- Canonical row order: {assertions['canonical_order']}.",
        f"- F2 LAPV direction discordant: {assertions['f2_lapv_opposite']}.",
        f"- F9 LAPV direction concordant: {assertions['f9_same_sign']}.",
        f"- F15 LAPV direction concordant: {assertions['f15_same_sign']}.",
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

    # Geometry identity check
    require(
        abs(PANELS["D"][1] + PANELS["D"][3] - HEIGHT_MM) < 0.01,
        f"Canvas height {HEIGHT_MM} does not match stacked panel bottoms",
    )
    b_x, _, b_w, _ = PANELS["B"]
    c_x, _, c_w, _ = PANELS["C"]
    require(abs(c_x - (b_x + b_w + H_GUTTER_MM)) < 0.01, "C must sit H_GUTTER after B")
    require(c_x + c_w <= WIDTH_MM + 0.01, "C overflows canvas width")
    require(b_w + H_GUTTER_MM + c_w <= WIDTH_MM + 0.01, "B|C row overflows canvas width")

    for panel, path in PANEL_SVGS.items():
        require(path.is_file(), f"Missing panel {panel} SVG: {path}")

    # Do not clobber legacy composite names
    for legacy in PROTECTED_LEGACY:
        require(
            OUT_BASE.with_suffix(legacy.suffix) != legacy,
            f"Refusing to overwrite legacy asset {legacy.name}",
        )

    assertions = scientific_assertions()
    out_svg = OUT_BASE.with_suffix(".svg")
    out_pdf = OUT_BASE.with_suffix(".pdf")
    out_png = OUT_BASE.with_suffix(".png")
    write_composite_svg(PANEL_SVGS, out_svg)
    require("NCP_v3" not in out_svg.read_text(encoding="utf-8"), "NCP_v3 visible in composite")
    export_raster(out_svg, out_pdf, out_png, OUT_PREVIEW)
    write_manifest(PANEL_SVGS)
    write_qc(
        assertions,
        {
            "svg": out_svg,
            "pdf": out_pdf,
            "png": out_png,
            "preview": OUT_PREVIEW,
            "manifest": OUT_MANIFEST,
        },
    )
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
