#!/usr/bin/env python3
"""Assemble Figure 5 A–D composite from composite-specific vector panel renders.

Panels are re-rendered at exact final physical sizes and nested at 100% scale.
Standalone PNG assets are never stitched or scaled.
Formal root / submission-pack Figure 5 assets are not replaced.
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import fig5a_cross_layer_triangulation_contract as fig5a
import fig5b_cross_layer_triangulation_matrix as fig5b
import fig5c_pathA_gate_architecture as fig5c
import fig5d_geomx_roi_spatial_evidence as fig5d

FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"

OUT_BASE = FIGURES / "figure5_cross_layer_triangulation"
OUT_PREVIEW = FIGURES / "figure5_cross_layer_triangulation_preview.png"
OUT_MANIFEST = RESULTS / "figure5_composite_manifest.tsv"
OUT_QC = RESULTS / "figure5_composite_qc.md"
OUT_CD_GEOM = RESULTS / "figure5cd_composite_geometry.tsv"

# Canvas: A / gutter / B / gutter / (C | gutter | D)
WIDTH_MM = 178.0
HEIGHT_MM = 239.0  # 48 + 3.5 + 96 + 3.5 + 88
GUTTER_V = 3.5
GUTTER_H = 4.0
PANEL_A_H = 48.0
PANEL_B_H = 96.0
PANEL_CD_H = 88.0

PANELS = {
    # (x, y_from_bottom, w, h)
    "A": (0.0, HEIGHT_MM - PANEL_A_H, 178.0, PANEL_A_H),
    "B": (0.0, HEIGHT_MM - PANEL_A_H - GUTTER_V - PANEL_B_H, 178.0, PANEL_B_H),
    "C": (0.0, 0.0, 76.0, PANEL_CD_H),
    "D": (76.0 + GUTTER_H, 0.0, 98.0, PANEL_CD_H),
}

# Convert to SVG top-left coordinates (y increases downward)
def panel_svg_xy(panel: str) -> tuple[float, float, float, float]:
    x, y_bottom, w, h = PANELS[panel]
    y_top = HEIGHT_MM - y_bottom - h
    return x, y_top, w, h


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prefix_svg_ids(svg: str, prefix: str) -> str:
    identifiers = re.findall(r'\bid="([^"]+)"', svg)
    for identifier in sorted(set(identifiers), key=len, reverse=True):
        replacement = f"{prefix}_{identifier}"
        svg = svg.replace(f'id="{identifier}"', f'id="{replacement}"')
        svg = svg.replace(f"url(#{identifier})", f"url(#{replacement})")
        svg = svg.replace(f'href="#{identifier}"', f'href="#{replacement}"')
        svg = svg.replace(f'xlink:href="#{identifier}"', f'xlink:href="#{replacement}"')
    return svg


def svg_body_and_viewbox(path: Path) -> tuple[str, str, float, float]:
    text = path.read_text(encoding="utf-8")
    width = re.search(r'<svg\b[^>]*\bwidth="([0-9.]+)(mm|pt)"', text)
    height = re.search(r'<svg\b[^>]*\bheight="([0-9.]+)(mm|pt)"', text)
    viewbox = re.search(r'\bviewBox="([^"]+)"', text)
    body = re.search(r"<svg\b[^>]*>(.*)</svg>\s*$", text, flags=re.DOTALL)
    require(all(v is not None for v in (width, height, viewbox, body)), f"Unable to parse {path.name}")
    assert width and height and viewbox and body
    scale = 1.0 if width.group(2) == "mm" else 25.4 / 72.0
    return (
        body.group(1),
        viewbox.group(1),
        float(width.group(1)) * scale,
        float(height.group(1)) * scale,
    )


def nested_svg(path: Path, panel: str) -> str:
    x, y, w, h = panel_svg_xy(panel)
    body, viewbox, src_w, src_h = svg_body_and_viewbox(path)
    require(
        abs(src_w - w) < 0.25 and abs(src_h - h) < 0.25,
        f"Panel {panel} source {src_w:.2f}×{src_h:.2f} != slot {w:.2f}×{h:.2f}",
    )
    body = prefix_svg_ids(body, f"panel{panel.lower()}")
    return (
        f'<g id="panel-{panel}" transform="translate({x},{y})">'
        f'<svg x="0" y="0" width="{w}" height="{h}" '
        f'viewBox="{viewbox}" preserveAspectRatio="xMinYMin meet">{body}</svg></g>'
    )


def write_composite_svg(panel_paths: dict[str, Path], output: Path) -> None:
    nested = "\n".join(nested_svg(panel_paths[p], p) for p in ("A", "B", "C", "D"))
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
    require(mutool is not None, "mutool is required for vector composite export")
    subprocess.run([mutool, "convert", "-o", str(pdf), str(svg)], check=True, capture_output=True, text=True)
    subprocess.run([mutool, "draw", "-q", "-r", "600", "-o", str(png), str(pdf)], check=True, capture_output=True, text=True)
    subprocess.run([mutool, "draw", "-q", "-r", "150", "-o", str(preview), str(pdf)], check=True, capture_output=True, text=True)


def render_panel_a(out_svg: Path) -> None:
    fig5a.apply_layout("composite")
    require(fig5a.STRIP_H_MM > 25.0, f"5A strip too small: {fig5a.STRIP_H_MM:.2f}")
    if not fig5a.OUT_CANONICAL.is_file():
        fig5a.write_manifest()
        fig5a.build_canonical_triangulation()
    fig, title_ax, col_axes, col_w, geomx_ax = fig5a.render()
    fig5a.validate_layout(fig, title_ax, col_axes, col_w, geomx_ax)
    fig.savefig(out_svg, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    fig5a.apply_layout("standalone")


def render_panel_b(out_svg: Path) -> None:
    fig5b.apply_layout("composite")
    _manifest, df = fig5b.load_and_assert()
    source = fig5b.build_source(df)
    fig, ax, rule_segments = fig5b.render(source)
    fig5b.validate_layout(fig, ax, rule_segments)
    fig.savefig(out_svg, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    fig5b.apply_layout("standalone")


def render_panel_c(out_svg: Path) -> None:
    fig5c.apply_layout("composite")
    combo = fig5c.load_and_assert()
    source, combo_source, detail = fig5c.build_sources(combo)
    fig, ax = fig5c.render(source, combo_source, detail)
    fig5c.validate_layout(fig, ax)
    fig.savefig(out_svg, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    fig5c.apply_layout("standalone")


def render_panel_d(out_svg: Path) -> None:
    fig5d.apply_layout("composite")
    df = fig5d.load_and_assert()
    W, meta, rows, markers = fig5d.build_sources(df)
    fig, ax = fig5d.render(W, meta, rows, markers)
    fig5d.validate_layout(fig, ax)
    fig.savefig(out_svg, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    fig5d.apply_layout("standalone")


def scientific_assertions() -> dict[str, object]:
    canonical = RESULTS / "figure5_triangulation_canonical.tsv"
    require(canonical.is_file(), f"Missing canonical: {canonical}")
    df = pd.read_csv(canonical, sep="\t")
    require(len(df) == 16, "expected 16 families")
    path_a = {int(x) for x in df.loc[df["Path_A"], "family_id"]}
    require(path_a == {10, 15}, f"Path A families changed: {path_a}")

    support = {int(r.family_id): int(r.support_count) for r in df.itertuples(index=False)}
    require(support[10] == 3 and support[15] == 3, "F10/F15 must be 3/3")
    require(support[9] == 2, "F9 must be 2/3")
    for f in (4, 5, 6):
        require(support[f] == 0, f"F{f} must be 0/3")
    counts = df["support_count"].value_counts().to_dict()
    require(int(counts.get(3, 0)) == 2, "3/3 count")
    require(int(counts.get(2, 0)) == 1, "2/3 count")
    require(int(counts.get(1, 0)) == 10, "1/3 count")
    require(int(counts.get(0, 0)) == 3, "0/3 count")

    summary_d = RESULTS / "fig5d_geomx_roi_spatial_evidence_summary.tsv"
    require(summary_d.is_file(), "5D summary missing")
    dsum = pd.read_csv(summary_d, sep="\t")
    require(int(dsum["n_rois"].iloc[0]) == 24, "5D n_rois")
    require(int(dsum["n_patients"].iloc[0]) == 4, "5D n_patients")
    f9 = dsum.loc[dsum["name"] == "F9"].iloc[0]
    require(int(f9["sr_n"]) == 12 and int(f9["af_n"]) == 12, "F9 SR/AF counts")

    return {
        "n_families": 16,
        "path_a": [f"F{i}" for i in sorted(path_a)],
        "support_counts": {int(k): int(v) for k, v in counts.items()},
        "f9_sr_af": (12, 12),
    }


def write_manifest(panel_paths: dict[str, Path]) -> None:
    rows = []
    for panel, path in panel_paths.items():
        x, y, w, h = panel_svg_xy(panel)
        rows.append(
            {
                "panel": panel,
                "renderer": f"scripts/ncp_v3_singscore_reanalysis/fig5{panel.lower()}_*.py",
                "layout_mode": "composite",
                "x_mm": x,
                "y_mm": y,
                "width_mm": w,
                "height_mm": h,
                "scale_percent": 100,
                "panel_svg": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                "panel_svg_sha256": sha256(path),
            }
        )
    pd.DataFrame(rows).to_csv(OUT_MANIFEST, sep="\t", index=False)


def write_qc(assertions: dict[str, object], outputs: dict[str, Path]) -> None:
    hashes = {str(p.relative_to(ROOT)): sha256(p) for p in outputs.values()}
    OUT_QC.write_text(
        "\n".join(
            [
                "# Figure 5 cross-layer triangulation composite QC",
                "",
                f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
                "",
                "Status: **PASS**",
                "",
                "## Physical geometry",
                f"- Canvas: {WIDTH_MM:.0f} × {HEIGHT_MM:.0f} mm.",
                "- A: 178 × 48 mm (full width).",
                "- B: 178 × 96 mm (full width; +4 mm for one-step type bump).",
                "- C: 76 × 88 mm; D: 98 × 88 mm; horizontal gutter 4 mm.",
                "- Vertical gutters: 3.5 mm.",
                "- Panels nested at 100% scale from composite-specific vector SVGs.",
                "- PNG exported at 600 dpi; no standalone PNG stitching.",
                "",
                "## C/D shared geometry contract",
                f"- `{OUT_CD_GEOM.relative_to(ROOT)}`",
                "- Shared panel-letter / title / subtitle / footer type system.",
                "- Shared top margin (2.5 mm) for letter baseline alignment.",
                "",
                "## Scientific assertions",
                f"- 16 canonical families; Path A = {assertions['path_a']}.",
                f"- Support-count distribution: {assertions['support_counts']}.",
                "- F10/F15 = 3/3 Path A; F9 = 2/3; F4/F5/F6 = 0/3.",
                "- Panel A: three Path A gates; HCA / dominant state / GeoMx are non-gates.",
                "- Panel D: 24 ROIs × 4 patients; F9 SR/AF = 12/12; GeoMx outside Path A.",
                "- No visible `NCP_v3`; no PriorityScore ordering.",
                "- Formal root / submission-pack Figure 5 assets were not replaced.",
                "",
                "## SHA256",
                *[f"- `{path}`: `{digest}`" for path, digest in sorted(hashes.items())],
                "",
            ]
        ),
        encoding="utf-8",
    )


def validate_composite(svg: Path, pdf: Path, png: Path) -> None:
    info = subprocess.run(["pdfinfo", str(pdf)], check=True, capture_output=True, text=True).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read composite PDF dimensions")
    w_mm = float(match.group(1)) * 25.4 / 72
    h_mm = float(match.group(2)) * 25.4 / 72
    require(abs(w_mm - WIDTH_MM) < 0.15, f"PDF width {w_mm:.2f} != {WIDTH_MM}")
    require(abs(h_mm - HEIGHT_MM) < 0.15, f"PDF height {h_mm:.2f} != {HEIGHT_MM}")

    svg_text = svg.read_text(encoding="utf-8")
    require("NCP_v3" not in svg_text, "NCP_v3 visible in composite")
    require("PriorityScore" not in svg_text, "PriorityScore visible in composite")
    require("Not a Path A gate" not in svg_text, "5D badge must remain omitted")
    for phrase in ("Cross-layer triangulation", "Rule-based Path A", "GeoMx ROI spatial context"):
        require(phrase in svg_text, f"missing panel title fragment: {phrase}")

    from PIL import Image

    img = Image.open(png)
    expect_w = round(WIDTH_MM / 25.4 * 600)
    expect_h = round(HEIGHT_MM / 25.4 * 600)
    require(abs(img.size[0] - expect_w) <= 2, f"PNG width px {img.size[0]} != {expect_w}")
    require(abs(img.size[1] - expect_h) <= 2, f"PNG height px {img.size[1]} != {expect_h}")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    require(OUT_CD_GEOM.is_file(), f"Missing C/D geometry contract: {OUT_CD_GEOM}")

    out_svg = OUT_BASE.with_suffix(".svg")
    out_pdf = OUT_BASE.with_suffix(".pdf")
    out_png = OUT_BASE.with_suffix(".png")

    with tempfile.TemporaryDirectory(prefix="fig5comp_") as tmp:
        tmpdir = Path(tmp)
        panel_paths = {
            "A": tmpdir / "fig5a_composite.svg",
            "B": tmpdir / "fig5b_composite.svg",
            "C": tmpdir / "fig5c_composite.svg",
            "D": tmpdir / "fig5d_composite.svg",
        }
        # Also keep durable composite-specific panel SVGs under figures/
        durable = {
            "A": FIGURES / "fig5a_cross_layer_triangulation_contract_composite.svg",
            "B": FIGURES / "fig5b_cross_layer_triangulation_matrix_composite.svg",
            "C": FIGURES / "fig5c_pathA_gate_architecture_composite.svg",
            "D": FIGURES / "fig5d_geomx_roi_spatial_evidence_composite.svg",
        }

        print("Rendering composite-specific panel A…")
        render_panel_a(panel_paths["A"])
        print("Rendering composite-specific panel B…")
        render_panel_b(panel_paths["B"])
        print("Rendering composite-specific panel C…")
        render_panel_c(panel_paths["C"])
        print("Rendering composite-specific panel D…")
        render_panel_d(panel_paths["D"])

        for key, src in panel_paths.items():
            shutil.copy2(src, durable[key])

        print("Assembling composite SVG…")
        write_composite_svg(panel_paths, out_svg)
        export_raster(out_svg, out_pdf, out_png, OUT_PREVIEW)

    assertions = scientific_assertions()
    validate_composite(out_svg, out_pdf, out_png)
    write_manifest(durable)
    write_qc(
        assertions,
        {
            "svg": out_svg,
            "pdf": out_pdf,
            "png": out_png,
            "preview": OUT_PREVIEW,
            "manifest": OUT_MANIFEST,
            "cd_geometry": OUT_CD_GEOM,
            **{f"panel_{k.lower()}_svg": v for k, v in durable.items()},
        },
    )
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_MANIFEST.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
