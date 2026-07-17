#!/usr/bin/env python3
"""Build and validate the exact-size final Figure 3 vector composite."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "v3"))

import fig3a_v3_schematic as fig3a
import fig3b_family_compartment_localization as fig3b
import fig3c_fidelity_relay_family_donor_examples as fig3c
import fig3d_relay_extension_family_donor_examples as fig3d

FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
OUT_BASE = FIGURES / "figure3_la_single_nucleus_localization"
OUT_PREVIEW = FIGURES / "figure3_la_single_nucleus_localization_preview.png"
OUT_MANIFEST = RESULTS / "figure3_composite_manifest.tsv"
OUT_QC = RESULTS / "figure3_composite_qc.md"
COMPOSITE_GEOMETRY = RESULTS / "figure3cd_case_study_composite_geometry.tsv"
AXIS_CONTRACT = RESULTS / "figure3cd_case_study_axis_contract.tsv"

WIDTH_MM = 178.0
HEIGHT_MM = 225.0
PANELS = {
    "A": (0.0, 0.0, 178.0, 48.0),
    "B": (0.0, 51.5, 178.0, 86.0),
    "C": (0.0, 141.0, 87.0, 84.0),
    "D": (91.0, 141.0, 87.0, 84.0),
}
PANEL_BASES = {
    "A": FIGURES / "fig3a_localization_bridge_composite",
    "B": FIGURES / "fig3b_family_compartment_localization_composite",
    "C": FIGURES / "fig3c_fidelity_relay_composite_specific",
    "D": FIGURES / "fig3d_relay_extension_composite_specific",
}
FOOTER = (
    "Donor×compartment family scores; nominal permutation p values are "
    "localization cues, not FDR-controlled discoveries."
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
    width = re.search(r'<svg\b[^>]*\bwidth="([0-9.]+)pt"', text)
    height = re.search(r'<svg\b[^>]*\bheight="([0-9.]+)pt"', text)
    viewbox = re.search(r'\bviewBox="([^"]+)"', text)
    body = re.search(r"<svg\b[^>]*>(.*)</svg>\s*$", text, flags=re.DOTALL)
    require(all(value is not None for value in (width, height, viewbox, body)), f"Unable to parse {path.name}")
    return (
        float(width.group(1)) * 25.4 / 72,
        float(height.group(1)) * 25.4 / 72,
        viewbox.group(1),
        body.group(1),
    )


def nested_svg(path: Path, panel: str) -> str:
    x, y, width, height = PANELS[panel]
    source_width, source_height, viewbox, body = svg_dimensions(path)
    require(abs(source_width - width) < 0.01 and abs(source_height - height) < 0.01, f"Panel {panel} source dimensions changed")
    values = [float(value) for value in viewbox.split()]
    require(len(values) == 4, f"Panel {panel} SVG viewBox is malformed")
    scale_x = width / values[2]
    scale_y = height / values[3]
    require(abs(scale_x - scale_y) < 1e-8, f"Panel {panel} would be non-uniformly scaled")
    body = prefix_svg_ids(body, f"panel{panel.lower()}")
    return (
        f'<g id="panel-{panel}" transform="translate({x},{y})">'
        f'<svg x="0" y="0" width="{width}" height="{height}" '
        f'viewBox="{viewbox}" preserveAspectRatio="xMinYMin meet">{body}</svg></g>'
    )


def text_layout_qc(fig: plt.Figure, panel: str) -> dict[str, float | int]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    texts = [
        text
        for text in fig.findobj(mpl.text.Text)
        if text.get_visible() and text.get_text().strip()
    ]
    min_font = min(float(text.get_fontsize()) for text in texts)
    violations = 0
    for text in texts:
        box = text.get_window_extent(renderer=renderer)
        if box.x0 < canvas.x0 - 1.0 or box.y0 < canvas.y0 - 1.0 or box.x1 > canvas.x1 + 1.0 or box.y1 > canvas.y1 + 1.0:
            violations += 1
    annotation_boxes = [
        text.get_window_extent(renderer=renderer)
        for text in texts
        if text.get_text().startswith("Ctrl/AF")
    ]
    annotation_overlaps = sum(
        first.overlaps(second)
        for index, first in enumerate(annotation_boxes)
        for second in annotation_boxes[index + 1 :]
    )
    require(min_font >= 6.0, f"Panel {panel} contains text below 6 pt")
    require(violations == 0, f"Panel {panel} contains {violations} clipped text elements")
    require(annotation_overlaps == 0, f"Panel {panel} contains overlapping numerical annotations")
    return {
        "min_font_pt": min_font,
        "bbox_violations": violations,
        "annotation_overlaps": annotation_overlaps,
    }


def save_panel(fig: plt.Figure, panel: str) -> tuple[Path, dict[str, float | int]]:
    base = PANEL_BASES[panel]
    svg, pdf, png = base.with_suffix(".svg"), base.with_suffix(".pdf"), base.with_suffix(".png")
    layout = text_layout_qc(fig, panel)
    fig.savefig(svg, bbox_inches=None, pad_inches=0)
    fig.savefig(pdf, bbox_inches=None, pad_inches=0)
    fig.savefig(png, dpi=600, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    width_mm, height_mm, _, _ = svg_dimensions(svg)
    expected = PANELS[panel]
    require(abs(width_mm - expected[2]) < 0.01 and abs(height_mm - expected[3]) < 0.01, f"Panel {panel} export dimensions changed")
    require("<image" not in svg.read_text(encoding="utf-8"), f"Panel {panel} unexpectedly contains raster content")
    return svg, layout


def load_and_render_panels() -> tuple[dict[str, Path], dict[str, dict[str, float | int]], dict[str, object]]:
    facts_a = fig3a.assert_authority()
    require(int(facts_a["n_families"]) == 16, "Panel A family count changed")
    require(int(facts_a["n_donors"]) == 34 and int(facts_a["n_compartments"]) == 9, "Panel A cohort dimensions changed")
    require(int(facts_a["n_perm"]) == 10000, "Panel A permutation count changed")

    source_b = pd.read_csv(fig3b.OUT_SOURCE, sep="\t")
    fig3b.assert_approved_source(source_b, fig3b.APPROVED_SCALE_M)
    require(source_b["family_id"].drop_duplicates().astype(int).tolist() == fig3b.EXPECTED_FAM_ORDER, "Panel B family order changed")
    require(source_b["compartment"].drop_duplicates().tolist() == fig3b.COMPARTMENTS, "Panel B compartment order changed")

    bridge_c, source3b_c, classes_c, hca_c = fig3c.load_family_inputs()
    summary_c = fig3c.build_summary(bridge_c, source3b_c, classes_c, hca_c)
    fig3c.load_sensitivity_qc(summary_c)
    geometry = pd.read_csv(COMPOSITE_GEOMETRY, sep="\t")
    geometry_c = geometry.loc[geometry["target_panel"].eq("Fig3C")].copy()
    require(set(geometry.loc[geometry["record_type"].eq("facet"), "width_mm"]) == {24.5}, "Composite facet width changed")

    donor_d, contract_d, metadata_d, _ = fig3d.load_inputs()
    source_d, summary_d = fig3d.build_outputs(donor_d, contract_d, metadata_d)
    axis_table = pd.read_csv(AXIS_CONTRACT, sep="\t")
    axis_rows = axis_table.loc[axis_table["contract_id"].eq("figure3cd_composite_v2")]
    require(len(axis_rows) == 1, "Composite C/D axis contract is missing")
    axis = axis_rows.iloc[0]
    require(float(axis["panel_width_mm"]) == 87.0 and float(axis["panel_height_mm"]) == 84.0, "Composite C/D dimensions changed")
    require(float(axis["annotation_font_pt"]) >= 6.0, "Composite annotation font is below 6 pt")

    figures = {
        "A": fig3a.render_composite_panel(),
        "B": fig3b.render_composite_panel(source_b, fig3b.APPROVED_SCALE_M),
        "C": fig3c.render(
            bridge_c,
            summary_c,
            geometry_c,
            composite=True,
            panel_width_mm=87.0,
            panel_height_mm=84.0,
            annotation_font_pt=float(axis["annotation_font_pt"]),
        ),
        "D": fig3d.render(source_d, summary_d, axis, composite=True),
    }
    paths: dict[str, Path] = {}
    layout: dict[str, dict[str, float | int]] = {}
    for panel, figure in figures.items():
        paths[panel], layout[panel] = save_panel(figure, panel)

    class_counts = source_b[["family_id", "fidelity_class"]].drop_duplicates()["fidelity_class"].value_counts()
    require([int(class_counts.get(name, 0)) for name in fig3b.CLASS_ORDER] == [1, 2, 1, 12], "Panel B class counts changed")
    scientific_qc(summary_c, summary_d)
    return paths, layout, {
        "facts_a": facts_a,
        "source_b": source_b,
        "summary_c": summary_c,
        "summary_d": summary_d,
        "axis": axis,
    }


def scientific_qc(summary_c: pd.DataFrame, summary_d: pd.DataFrame) -> None:
    def c_row(family: int, compartment: str) -> pd.Series:
        rows = summary_c.loc[summary_c["family_id"].eq(family) & summary_c["compartment"].eq(compartment)]
        require(len(rows) == 1, f"Missing F{family} {compartment} assertion row")
        return rows.iloc[0]

    require(str(c_row(10, "fibroblast")["fidelity_class"]) == "Fidelity", "F10 class changed")
    require(bool(c_row(10, "fibroblast")["is_carrier"]) and bool(c_row(10, "fibroblast")["is_supported"]), "F10 carrier/support changed")
    require(str(c_row(11, "endothelial")["fidelity_class"]) == "Relay", "F11 class changed")
    require(bool(c_row(11, "endothelial")["is_carrier"]) and not bool(c_row(11, "endothelial")["is_supported"]), "F11 Endothelial call changed")
    require(not bool(c_row(11, "adipocyte")["is_carrier"]) and bool(c_row(11, "adipocyte")["is_supported"]), "F11 Adipocyte call changed")

    def d_row(family: int, compartment: str) -> pd.Series:
        rows = summary_d.loc[summary_d["family_id"].eq(family) & summary_d["compartment"].eq(compartment)]
        require(len(rows) == 1, f"Missing F{family} {compartment} assertion row")
        return rows.iloc[0]

    require(str(d_row(15, "fibroblast")["fidelity_class"]) == "Relay", "F15 class changed")
    require(not bool(d_row(15, "fibroblast")["is_supported"]), "F15 Fibroblast support changed")
    require(bool(d_row(15, "endothelial")["is_supported"]), "F15 Endothelial support changed")
    require(not bool(d_row(15, "lymphocyte")["is_supported"]), "F15 Lymphocyte support changed")
    require(str(d_row(16, "fibroblast")["fidelity_class"]) == "Extension", "F16 class changed")
    require(bool(d_row(16, "fibroblast")["is_supported"]) and bool(d_row(16, "atrialCM")["is_supported"]), "F16 supported set changed")
    require(not bool(d_row(16, "adipocyte")["is_supported"]), "F16 Adipocyte support changed")
    require(set(summary_c["n_perm"].astype(int)) == {10000} and set(summary_d["n_perm"].astype(int)) == {10000}, "Primary permutation count changed")
    require(int(d_row(15, "endothelial")["n_perm_sensitivity"]) == 1000000, "F15 precision-sensitivity audit missing")


def write_composite_svg(panel_paths: dict[str, Path], output: Path) -> None:
    nested = "\n".join(nested_svg(panel_paths[panel], panel) for panel in ("A", "B", "C", "D"))
    footer = (
        '<text id="shared-cd-footer" x="89" y="222.2" text-anchor="middle" '
        'font-family="DejaVu Sans" font-size="2.1167" fill="#5D5D5D">'
        f"{FOOTER}</text>"
    )
    output.write_text(
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{WIDTH_MM}mm" height="{HEIGHT_MM}mm" viewBox="0 0 {WIDTH_MM} {HEIGHT_MM}">\n'
        f'<rect x="0" y="0" width="{WIDTH_MM}" height="{HEIGHT_MM}" fill="#FFFFFF"/>\n'
        f"{nested}\n{footer}\n</svg>\n",
        encoding="utf-8",
    )


def palette_qc() -> dict[str, float]:
    colors = [fig3b.CLASS_COLORS[name] for name in fig3b.CLASS_ORDER]
    rgb = np.asarray([mpl.colors.to_rgb(color) for color in colors])
    luminance = rgb @ np.asarray([0.2126, 0.7152, 0.0722])
    pairwise = [float(np.linalg.norm(rgb[i] - rgb[j])) for i in range(len(rgb)) for j in range(i + 1, len(rgb))]
    gray_pairwise = [float(abs(luminance[i] - luminance[j])) for i in range(len(rgb)) for j in range(i + 1, len(rgb))]
    b_qc = fig3b.palette_accessibility_qc()
    return {
        "class_rgb_min_distance": min(pairwise),
        "class_gray_min_distance": min(gray_pairwise),
        "heat_grayscale_sign_distance": float(b_qc["grayscale_heat_sign"]),
        "heat_min_cvd_distance": min(float(b_qc["cvd_heat_sign"]), float(b_qc["cvd_class"])),
    }


def export_and_validate(panel_paths: dict[str, Path], layout: dict[str, dict[str, float | int]], data: dict[str, object]) -> None:
    mutool = shutil.which("mutool")
    pdfinfo = shutil.which("pdfinfo")
    pdftotext = shutil.which("pdftotext")
    require(all(tool is not None for tool in (mutool, pdfinfo, pdftotext)), "Required PDF export/QC tools are unavailable")
    out_svg = OUT_BASE.with_suffix(".svg")
    out_pdf = OUT_BASE.with_suffix(".pdf")
    out_png = OUT_BASE.with_suffix(".png")
    write_composite_svg(panel_paths, out_svg)
    subprocess.run([mutool, "convert", "-o", str(out_pdf), str(out_svg)], check=True, capture_output=True, text=True)
    subprocess.run([mutool, "draw", "-q", "-r", "600", "-o", str(out_png), str(out_pdf)], check=True, capture_output=True, text=True)
    subprocess.run([mutool, "draw", "-q", "-r", "150", "-o", str(OUT_PREVIEW), str(out_pdf)], check=True, capture_output=True, text=True)

    info = subprocess.run([pdfinfo, str(out_pdf)], check=True, capture_output=True, text=True).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read final PDF dimensions")
    pdf_width = float(match.group(1)) * 25.4 / 72
    pdf_height = float(match.group(2)) * 25.4 / 72
    require(abs(pdf_width - WIDTH_MM) < 0.05 and abs(pdf_height - HEIGHT_MM) < 0.05, "Final PDF dimensions changed")
    image = plt.imread(out_png)
    require(
        abs(image.shape[1] - round(WIDTH_MM / 25.4 * 600)) <= 1
        and abs(image.shape[0] - round(HEIGHT_MM / 25.4 * 600)) <= 1,
        "Final 600-dpi PNG dimensions changed",
    )
    preview = plt.imread(OUT_PREVIEW)
    require(
        abs(preview.shape[1] - round(WIDTH_MM / 25.4 * 150)) <= 1
        and abs(preview.shape[0] - round(HEIGHT_MM / 25.4 * 150)) <= 1,
        "Final preview dimensions changed",
    )
    visible = subprocess.run([pdftotext, "-layout", str(out_pdf), "-"], check=True, capture_output=True, text=True).stdout
    normalized = " ".join(visible.split())
    require("NCP_v3" not in visible, "Visible NCP_v3 text remains")
    require("program-level" not in visible.lower() and "gene-level" not in visible.lower(), "Forbidden analysis-level result remains visible")
    require(normalized.count(FOOTER) == 1, "Shared C/D footer is missing or duplicated")
    require(normalized.count("Family score (centered rank units)") == 1, "C/D y-axis title is missing or duplicated")
    for phrase in (
        "LA single-nucleus localization and carrier-relative classification",
        "LA compartment localization and carrier-relative classification",
        "Fidelity and Relay examples",
        "Relay and Extension examples",
    ):
        require(phrase in normalized, f"Visible title is incomplete: {phrase}")
    svg_text = out_svg.read_text(encoding="utf-8")
    require("<image" not in svg_text, "Final SVG contains unexpected raster content")
    for badge_id in (
        "panelc_fidelity-badge",
        "panelc_relay-badge",
        "paneld_relay-badge",
        "paneld_extension-badge",
    ):
        require(svg_text.count(f'id="{badge_id}"') == 1, f"Badge {badge_id} missing or duplicated")
    for panel in ("A", "B", "C", "D"):
        require(svg_text.count(f'id="panel-{panel}"') == 1, f"Panel {panel} is missing or duplicated")
    require(svg_text.count('id="shared-cd-footer"') == 1, "SVG shared footer is missing or duplicated")
    write_manifest(panel_paths)
    write_qc(layout, data, out_svg, out_pdf, out_png)


def write_manifest(panel_paths: dict[str, Path]) -> None:
    sources = {
        "A": "results/v3/fig3a_v3_schematic_source.tsv",
        "B": str(fig3b.OUT_SOURCE.relative_to(ROOT)),
        "C": f"{fig3c.BRIDGE.relative_to(ROOT)};{fig3c.OUT_SUMMARY.relative_to(ROOT)}",
        "D": f"{fig3d.OUT_SOURCE.relative_to(ROOT)};{fig3d.OUT_SUMMARY.relative_to(ROOT)}",
    }
    renderers = {
        "A": "scripts/v3/fig3a_v3_schematic.py",
        "B": "scripts/ncp_v3_singscore_reanalysis/fig3b_family_compartment_localization.py",
        "C": "scripts/ncp_v3_singscore_reanalysis/fig3c_fidelity_relay_family_donor_examples.py",
        "D": "scripts/ncp_v3_singscore_reanalysis/fig3d_relay_extension_family_donor_examples.py",
    }
    rows = []
    for panel, path in panel_paths.items():
        x, y, width, height = PANELS[panel]
        rows.append(
            {
                "panel": panel,
                "source": sources[panel],
                "renderer": renderers[panel],
                "x_mm": x,
                "y_mm": y,
                "width_mm": width,
                "height_mm": height,
                "scale_percent": 100,
                "scale_x": 1.0,
                "scale_y": 1.0,
                "asset_status": "vector",
                "panel_svg": str(path.relative_to(ROOT)),
                "panel_svg_sha256": sha256(path),
            }
        )
    pd.DataFrame(rows).to_csv(OUT_MANIFEST, sep="\t", index=False)


def write_qc(
    layout: dict[str, dict[str, float | int]],
    data: dict[str, object],
    svg: Path,
    pdf: Path,
    png: Path,
) -> None:
    palette = palette_qc()
    source_b = data["source_b"]
    class_counts = source_b[["family_id", "fidelity_class"]].drop_duplicates()["fidelity_class"].value_counts()
    hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (svg, pdf, png, OUT_PREVIEW, OUT_MANIFEST, COMPOSITE_GEOMETRY, AXIS_CONTRACT)
    }
    lines = [
        "# Final Figure 3 composite QC",
        "",
        "Status: **PASS**",
        "",
        "## Phase 6B0-R display correction",
        "- scientific_data_changed = NO",
        "- scientific_rerun = NO",
        "- display_text_changed = YES",
        "- change_reason = Figure 3 score-computability threshold clarification",
        "- decision_code = PRODUCTION_3_ONLY_UPDATE_DISPLAY_NO_RERUN",
        "",
        "## Physical geometry and vector assembly",
        f"- Final canvas: {WIDTH_MM:.0f} × {HEIGHT_MM:.0f} mm.",
        "- A: (0, 0), 178 × 48 mm; B: (0, 51.5), 178 × 86 mm.",
        "- C: (0, 141), 87 × 84 mm; D: (91, 141), 87 × 84 mm.",
        "- All panels nested at 100% physical scale with equal x/y transform; no stretching.",
        "- PDF/SVG retain vector text, symbols, donor plots, and the B matrix; PNG is 600 dpi.",
        "- Exactly one C/D y-axis title, one shared footer, and one badge per case block.",
        "",
        "## Typography and containment",
        *[
            f"- Panel {panel}: minimum visible font {float(values['min_font_pt']):.1f} pt; "
            f"text bbox violations {int(values['bbox_violations'])}; "
            f"numerical-annotation overlaps {int(values['annotation_overlaps'])}."
            for panel, values in layout.items()
        ],
        "- Right-edge labels, class badges, module glyphs, and panel titles are contained.",
        "",
        "## Scientific assertions",
        "- A: 16-family axis; 34 donors; 9 compartments; 10,000 primary permutations.",
        "- B: canonical 16-family order and 9-compartment order retained; carrier outlines and nominal-support dots retained.",
        f"- B class counts: Fidelity {int(class_counts['Fidelity'])}; Relay {int(class_counts['Relay'])}; "
        f"Extension {int(class_counts['Extension'])}; Below-detection {int(class_counts['Below_detection'])}.",
        "- F10 Fidelity: Fibroblast is the sole carrier and supported compartment.",
        "- F11 Relay: Endothelial is carrier/not supported; Adipocyte is non-carrier/supported.",
        "- F15 Relay: Fibroblast FALSE, Endothelial TRUE, Lymphocyte FALSE support; precision-sensitivity audit retained.",
        "- F16 Extension: Fibroblast and aCM TRUE, Adipocyte FALSE support.",
        "- C/D share y limits, ticks, point/mean/zero-line marks, block heights, separator, badge dimensions, and 6-pt annotation style.",
        "- Renderers open frozen family-level sources only; no gene- or program-level result is shown.",
        "- No visible `NCP_v3` text.",
        "",
        "## Palette accessibility",
        "- Semantic class palette matches Figure 3B in A–D.",
        f"- Minimum class RGB distance: {palette['class_rgb_min_distance']:.4f}.",
        f"- Minimum class grayscale luminance distance: {palette['class_gray_min_distance']:.4f}.",
        f"- Heatmap grayscale sign distance: {palette['heat_grayscale_sign_distance']:.4f}.",
        f"- Heatmap minimum simulated CVD distance: {palette['heat_min_cvd_distance']:.4f}.",
        "",
        "## SHA256",
        *[f"- `{path}`: `{digest}`" for path, digest in hashes.items()],
        "",
    ]
    OUT_QC.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    panel_paths, layout, data = load_and_render_panels()
    export_and_validate(panel_paths, layout, data)
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_MANIFEST.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
