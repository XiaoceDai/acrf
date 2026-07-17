#!/usr/bin/env python3
"""Assemble vector composite-specific Figure 3C/3D donor-example panels."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt

import fig3c_fidelity_relay_family_donor_examples as fig3c
import fig3d_relay_extension_family_donor_examples as fig3d


ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "figures" / "ncp_v3_singscore_reanalysis"
RESULTS = ROOT / "results" / "ncp_v3_singscore_reanalysis"
OUT_BASE = FIGURES / "fig3cd_family_donor_examples_composite"
OUT_QC = RESULTS / "fig3cd_family_donor_examples_composite_qc.md"
WIDTH_MM = 176.0
HEIGHT_MM = 82.0


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
    ids = re.findall(r'\bid="([^"]+)"', svg)
    for identifier in sorted(ids, key=len, reverse=True):
        replacement = f"{prefix}_{identifier}"
        svg = svg.replace(f'id="{identifier}"', f'id="{replacement}"')
        svg = svg.replace(f"url(#{identifier})", f"url(#{replacement})")
        svg = svg.replace(f'href="#{identifier}"', f'href="#{replacement}"')
        svg = svg.replace(f'xlink:href="#{identifier}"', f'xlink:href="#{replacement}"')
    return svg


def nested_svg(path: Path, prefix: str, x_mm: float) -> str:
    text = prefix_svg_ids(path.read_text(encoding="utf-8"), prefix)
    viewbox = re.search(r'viewBox="([^"]+)"', text)
    body = re.search(r"<svg\b[^>]*>(.*)</svg>\s*$", text, flags=re.DOTALL)
    require(viewbox is not None and body is not None, f"Unable to parse {path.name}")
    return (
        f'<g transform="translate({x_mm},0)">'
        f'<svg x="0" y="0" width="88" height="82" '
        f'viewBox="{viewbox.group(1)}" preserveAspectRatio="none">'
        f"{body.group(1)}</svg></g>"
    )


def write_composite_svg(panel_c: Path, panel_d: Path, output: Path) -> None:
    footer = (
        '<text x="88" y="76.0" text-anchor="middle" font-family="DejaVu Sans" '
        'font-size="2.1167" fill="#5D5D5D">'
        "Donor×compartment family scores; nominal permutation p values are localization"
        '<tspan x="88" dy="2.75">cues, not FDR-controlled discoveries.</tspan>'
        "</text>"
    )
    svg = (
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" '
        'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{WIDTH_MM}mm" height="{HEIGHT_MM}mm" '
        f'viewBox="0 0 {WIDTH_MM} {HEIGHT_MM}">\n'
        '<rect x="0" y="0" width="176" height="82" fill="#FFFFFF"/>\n'
        f"{nested_svg(panel_c, 'panelc', 0.0)}\n"
        f"{nested_svg(panel_d, 'paneld', 88.0)}\n"
        f"{footer}\n"
        "</svg>\n"
    )
    output.write_text(svg, encoding="utf-8")


def render_panels(panel_c: Path, panel_d: Path) -> None:
    bridge, fig3b_source, classes, hca = fig3c.load_family_inputs()
    summary_c = fig3c.build_summary(bridge, fig3b_source, classes, hca)
    geometry = fig3c.load_geometry()
    fig_c = fig3c.render(bridge, summary_c, geometry, composite=True)
    fig3c.validate_figure(fig_c, summary_c)
    fig_c.savefig(panel_c, bbox_inches=None, pad_inches=0)
    plt.close(fig_c)

    donor, contract, metadata, axis = fig3d.load_inputs()
    source_d, summary_d = fig3d.build_outputs(donor, contract, metadata)
    fig_d = fig3d.render(source_d, summary_d, axis, composite=True)
    fig3d.validate_figure(fig_d)
    fig_d.savefig(panel_d, bbox_inches=None, pad_inches=0)
    plt.close(fig_d)


def validate_outputs(svg: Path, pdf: Path, png: Path) -> None:
    pdfinfo = shutil.which("pdfinfo")
    pdftotext = shutil.which("pdftotext")
    require(pdfinfo is not None and pdftotext is not None, "PDF QC tools unavailable")
    info = subprocess.run(
        [pdfinfo, str(pdf)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read composite PDF dimensions")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    require(
        abs(width_mm - WIDTH_MM) < 0.05 and abs(height_mm - HEIGHT_MM) < 0.05,
        "Composite PDF dimensions changed",
    )
    text = subprocess.run(
        [pdftotext, "-layout", str(pdf), "-"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    normalized = " ".join(text.split())
    footer = (
        "Donor×compartment family scores; nominal permutation p values are "
        "localization cues, not FDR-controlled discoveries."
    )
    require(normalized.count(footer) == 1, "Composite shared footer is missing or duplicated")
    require(
        normalized.count("Family score (centered rank units)") == 1,
        "Composite shared y-axis title is missing or duplicated",
    )
    require(
        "C Fidelity and Relay examples" in normalized
        and "D Relay and Extension examples" in normalized,
        "Composite panel titles are incomplete",
    )
    image = plt.imread(png)
    expected_width = round(WIDTH_MM / 25.4 * 600)
    expected_height = round(HEIGHT_MM / 25.4 * 600)
    require(
        abs(image.shape[1] - expected_width) <= 1
        and abs(image.shape[0] - expected_height) <= 1,
        "Composite 600-dpi PNG dimensions changed",
    )
    svg_text = svg.read_text(encoding="utf-8")
    require(svg_text.count("FDR-controlled discoveries.") == 1, "SVG footer duplicated")
    require(svg_text.count('id="paneld_extension-badge"') == 1, "Composite SVG Extension badge duplicated")


def main() -> None:
    mutool = shutil.which("mutool")
    require(mutool is not None, "mutool is required for vector composite export")
    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    out_svg = OUT_BASE.with_suffix(".svg")
    out_pdf = OUT_BASE.with_suffix(".pdf")
    out_png = OUT_BASE.with_suffix(".png")
    with tempfile.TemporaryDirectory(prefix="fig3cd_") as tmp:
        tmpdir = Path(tmp)
        panel_c = tmpdir / "panel_c.svg"
        panel_d = tmpdir / "panel_d.svg"
        render_panels(panel_c, panel_d)
        write_composite_svg(panel_c, panel_d, out_svg)
    subprocess.run(
        [mutool, "convert", "-o", str(out_pdf), str(out_svg)],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        [mutool, "draw", "-q", "-r", "600", "-o", str(out_png), str(out_pdf)],
        check=True,
        capture_output=True,
        text=True,
    )
    validate_outputs(out_svg, out_pdf, out_png)
    hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (out_pdf, out_svg, out_png)
    }
    OUT_QC.write_text(
        "\n".join(
            [
                "# Figure 3C/3D composite-specific QC",
                "",
                "Status: **PASS**",
                "",
                "- Figure 3C carries the sole shared y-axis title.",
                "- Figure 3D omits the repeated y-axis title.",
                "- Standalone panel footers are omitted.",
                "- Exactly one shared footer is rendered beneath both panels.",
                "- Composite dimensions: 176 × 82 mm.",
                "- PDF and SVG are vector exports; PNG is 600 dpi.",
                "- Extension badge duplicate/clipping assertions inherited from Figure 3D and rechecked in SVG.",
                "",
                "## SHA256",
                *[f"- `{path}`: `{digest}`" for path, digest in hashes.items()],
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
