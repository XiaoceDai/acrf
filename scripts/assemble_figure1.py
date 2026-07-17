#!/usr/bin/env python3
"""Assemble the frozen Figure 1A-D vector panels into one publication figure."""
from __future__ import annotations

import base64
import hashlib
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
OUT_BASE = FIGURES / "figure1_composite"
OUT_PREVIEW = FIGURES / "figure1_composite_preview.png"
OUT_MANIFEST = RESULTS / "figure1_composite_manifest.tsv"
OUT_QC = RESULTS / "figure1_composite_qc.md"

PAGE_WIDTH_MM = 183.0
PAGE_HEIGHT_MM = 254.5
OUTER_MARGIN_MM = 4.0
GUTTER_MM = 4.5
LEFT_WIDTH_MM = 68.2
RIGHT_WIDTH_MM = 102.3
TOP_HEIGHT_MM = 105.0
BOTTOM_HEIGHT_MM = 137.0
LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"

PANELS = [
    {
        "panel": "A",
        "pdf": FIGURES / "fig1a_donor_paired_estimand.pdf",
        "svg": FIGURES / "fig1a_donor_paired_estimand.svg",
        "x_mm": OUTER_MARGIN_MM,
        "y_mm": OUTER_MARGIN_MM,
        "width_mm": LEFT_WIDTH_MM,
        "height_mm": TOP_HEIGHT_MM,
    },
    {
        "panel": "B",
        "pdf": FIGURES
        / "ncp_v3_singscore_reanalysis"
        / "fig1b_healthy_reference_direction_magnitude.pdf",
        "svg": FIGURES
        / "ncp_v3_singscore_reanalysis"
        / "fig1b_healthy_reference_direction_magnitude.svg",
        "x_mm": OUTER_MARGIN_MM + LEFT_WIDTH_MM + GUTTER_MM,
        "y_mm": OUTER_MARGIN_MM,
        "width_mm": RIGHT_WIDTH_MM,
        "height_mm": TOP_HEIGHT_MM,
    },
    {
        "panel": "C",
        "pdf": FIGURES / "fig1c_selection_provenance.pdf",
        "svg": FIGURES / "fig1c_selection_provenance.svg",
        "x_mm": OUTER_MARGIN_MM,
        "y_mm": OUTER_MARGIN_MM + TOP_HEIGHT_MM + GUTTER_MM,
        "width_mm": LEFT_WIDTH_MM,
        "height_mm": BOTTOM_HEIGHT_MM,
    },
    {
        "panel": "D",
        "pdf": FIGURES / "fig1d_program_family_topology.pdf",
        "svg": FIGURES / "fig1d_program_family_topology.svg",
        "x_mm": OUTER_MARGIN_MM + LEFT_WIDTH_MM + GUTTER_MM,
        "y_mm": OUTER_MARGIN_MM + TOP_HEIGHT_MM + GUTTER_MM,
        "width_mm": RIGHT_WIDTH_MM,
        "height_mm": BOTTOM_HEIGHT_MM,
    },
]

FROZEN_INPUTS = [
    RESULTS / "innovation_ncp_v3_panel.tsv",
    RESULTS / "innovation_ncp_v3_axis_families.tsv",
    RESULTS / "ncp_v3_canonical_display_order.tsv",
    RESULTS
    / "ncp_v3_singscore_reanalysis"
    / "fig1b_canonical_display_order.tsv",
]
PANEL_RENDERERS = [
    ROOT / "scripts" / "fig1a_donor_paired_estimand.py",
    ROOT
    / "scripts"
    / "ncp_v3_singscore_reanalysis"
    / "fig1b_healthy_reference_direction_magnitude.py",
    ROOT / "scripts" / "fig1c_healthy_selection_frozen_deployment.py",
    ROOT / "scripts" / "fig1d_program_family_topology.py",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def svg_dimensions_mm(path: Path) -> tuple[float, float]:
    root = ET.parse(path).getroot()

    def parse_mm(value: str) -> float:
        match = re.fullmatch(r"\s*([0-9.]+)\s*(mm|pt)?\s*", value)
        require(match is not None, f"Unsupported SVG dimension in {path}: {value}")
        number = float(match.group(1))
        return number * 25.4 / 72.0 if match.group(2) == "pt" else number

    return parse_mm(root.attrib["width"]), parse_mm(root.attrib["height"])


def validate_geometry() -> None:
    require(
        abs(
            2 * OUTER_MARGIN_MM
            + LEFT_WIDTH_MM
            + GUTTER_MM
            + RIGHT_WIDTH_MM
            - PAGE_WIDTH_MM
        )
        < 1e-9,
        "Figure 1 horizontal geometry does not close",
    )
    require(
        abs(
            2 * OUTER_MARGIN_MM
            + TOP_HEIGHT_MM
            + GUTTER_MM
            + BOTTOM_HEIGHT_MM
            - PAGE_HEIGHT_MM
        )
        < 1e-9,
        "Figure 1 vertical geometry does not close",
    )
    require(
        abs(LEFT_WIDTH_MM / (LEFT_WIDTH_MM + RIGHT_WIDTH_MM) - 0.40) < 0.001,
        "Figure 1 columns are not in the requested 0.40:0.60 ratio",
    )
    require(BOTTOM_HEIGHT_MM > TOP_HEIGHT_MM, "Bottom row must be taller than top row")

    for panel in PANELS:
        for kind in ("pdf", "svg"):
            require(panel[kind].is_file(), f"Missing panel {panel['panel']} {kind}: {panel[kind]}")
        actual_width, actual_height = svg_dimensions_mm(panel["svg"])
        require(
            abs(actual_width - panel["width_mm"]) < 0.05
            and abs(actual_height - panel["height_mm"]) < 0.05,
            f"Panel {panel['panel']} SVG dimensions are {actual_width} x {actual_height} mm, "
            f"expected {panel['width_mm']} x {panel['height_mm']} mm",
        )


def validate_frozen_content() -> None:
    for path in FROZEN_INPUTS:
        require(path.is_file(), f"Missing frozen Figure 1 input: {path}")

    b_summary = pd.read_csv(
        RESULTS
        / "ncp_v3_singscore_reanalysis"
        / "fig1b_healthy_reference_direction_magnitude_summary.tsv",
        sep="\t",
    ).iloc[0]
    c_summary_frame = pd.read_csv(RESULTS / "fig1c_selection_provenance_summary.tsv", sep="\t")
    d_summary_frame = pd.read_csv(RESULTS / "fig1d_program_family_topology_summary.tsv", sep="\t")
    c_summary = dict(zip(c_summary_frame["metric"], c_summary_frame["value"]))
    d_summary = dict(zip(d_summary_frame["metric"], d_summary_frame["value"]))

    require(int(b_summary["n_programs"]) == 30, "Figure 1B program count changed")
    require(
        int(b_summary["n_la_biased"]) == 22 and int(b_summary["n_ra_biased"]) == 8,
        "Figure 1B direction split changed",
    )
    require(
        int(c_summary["candidate_pathways"]) == 80
        and int(c_summary["strict_eligible"]) == 38
        and int(c_summary["strict_selected"]) == 25
        and int(c_summary["fallback_selected"]) == 5
        and int(c_summary["final_frozen"]) == 30,
        "Figure 1C frozen selection counts changed",
    )
    require(
        int(d_summary["frozen_programs"]) == 30
        and int(d_summary["fixed_families"]) == 16,
        "Figure 1D frozen topology counts changed",
    )

    b_source = pd.read_csv(
        RESULTS
        / "ncp_v3_singscore_reanalysis"
        / "fig1b_healthy_reference_direction_magnitude_source.tsv",
        sep="\t",
    )
    d_source = pd.read_csv(RESULTS / "fig1d_program_family_topology_source.tsv", sep="\t")
    b_order = (
        b_source.sort_values("figure_display_order")
        .drop_duplicates("program_id")["program_id"]
        .tolist()
    )
    d_order = d_source.sort_values("canonical_display_order")["program_id"].tolist()
    require(b_order == d_order, "Figure 1B and 1D program order no longer matches")


def validate_presentation_sources() -> None:
    for path in PANEL_RENDERERS:
        text = path.read_text(encoding="utf-8")
        font_sizes = [
            float(match)
            for match in re.findall(r"fontsize\s*=\s*([0-9]+(?:\.[0-9]+)?)", text)
        ]
        require(font_sizes and min(font_sizes) >= 6.5, f"Text below 6.5 pt detected in {path}")
        require("DejaVu Sans" in text, f"Shared DejaVu Sans font is not declared in {path}")
    for path in PANEL_RENDERERS[1:]:
        text = path.read_text(encoding="utf-8")
        require(LA_COLOR in text and RA_COLOR in text, f"Canonical LA/RA colors missing in {path}")


def write_svg() -> None:
    image_nodes = []
    for panel in PANELS:
        encoded = base64.b64encode(panel["svg"].read_bytes()).decode("ascii")
        image_nodes.append(
            f'  <image id="panel-{panel["panel"]}" x="{panel["x_mm"]}mm" '
            f'y="{panel["y_mm"]}mm" width="{panel["width_mm"]}mm" '
            f'height="{panel["height_mm"]}mm" preserveAspectRatio="none" '
            f'href="data:image/svg+xml;base64,{encoded}"/>'
        )
    svg = (
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{PAGE_WIDTH_MM}mm" height="{PAGE_HEIGHT_MM}mm" '
        f'viewBox="0 0 {PAGE_WIDTH_MM} {PAGE_HEIGHT_MM}">\n'
        f'  <rect width="{PAGE_WIDTH_MM}" height="{PAGE_HEIGHT_MM}" fill="#FFFFFF"/>\n'
        + "\n".join(image_nodes)
        + "\n</svg>\n"
    )
    OUT_BASE.with_suffix(".svg").write_text(svg, encoding="utf-8")


def write_pdf() -> None:
    pdflatex = shutil.which("pdflatex")
    require(pdflatex is not None, "pdflatex is required for vector PDF assembly")
    with tempfile.TemporaryDirectory(prefix="figure1-composite-") as tmp:
        tmp_path = Path(tmp)
        tex_path = tmp_path / "figure1_composite.tex"
        puts = []
        for panel in PANELS:
            y_from_bottom = PAGE_HEIGHT_MM - panel["y_mm"] - panel["height_mm"]
            puts.append(
                rf"\put(\LenToUnit{{{panel['x_mm']}mm}},\LenToUnit{{{y_from_bottom}mm}})"
                rf"{{\includegraphics[width={panel['width_mm']}mm,height={panel['height_mm']}mm]"
                rf"{{\detokenize{{{panel['pdf']}}}}}}}"
            )
        tex = rf"""\documentclass{{article}}
\usepackage[paperwidth={PAGE_WIDTH_MM}mm,paperheight={PAGE_HEIGHT_MM}mm,margin=0mm]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{eso-pic}}
\pagestyle{{empty}}
\setlength{{\unitlength}}{{1mm}}
\begin{{document}}
\AddToShipoutPictureFG*{{\AtPageLowerLeft{{\begin{{picture}}(0,0)
{chr(10).join(puts)}
\end{{picture}}}}}}
\null
\end{{document}}
"""
        tex_path.write_text(tex, encoding="utf-8")
        completed = subprocess.run(
            [pdflatex, "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
            cwd=tmp_path,
            check=False,
            capture_output=True,
            text=True,
        )
        require(completed.returncode == 0, f"pdflatex failed:\n{completed.stdout}\n{completed.stderr}")
        shutil.copy2(tmp_path / "figure1_composite.pdf", OUT_BASE.with_suffix(".pdf"))


def rasterize() -> None:
    pdftocairo = shutil.which("pdftocairo")
    require(pdftocairo is not None, "pdftocairo is required for Figure 1 PNG exports")
    subprocess.run(
        [
            pdftocairo,
            "-png",
            "-singlefile",
            "-r",
            "600",
            str(OUT_BASE.with_suffix(".pdf")),
            str(OUT_BASE),
        ],
        check=True,
    )
    preview_stem = OUT_PREVIEW.with_suffix("")
    subprocess.run(
        [
            pdftocairo,
            "-png",
            "-singlefile",
            "-r",
            "150",
            str(OUT_BASE.with_suffix(".pdf")),
            str(preview_stem),
        ],
        check=True,
    )


def validate_outputs() -> None:
    pdfinfo = shutil.which("pdfinfo")
    pdftotext = shutil.which("pdftotext")
    require(pdfinfo is not None and pdftotext is not None, "PDF QC tools are required")
    info = subprocess.run(
        [pdfinfo, str(OUT_BASE.with_suffix(".pdf"))],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "Unable to read composite PDF dimensions")
    width_mm = float(match.group(1)) * 25.4 / 72.0
    height_mm = float(match.group(2)) * 25.4 / 72.0
    require(
        abs(width_mm - PAGE_WIDTH_MM) < 0.05 and abs(height_mm - PAGE_HEIGHT_MM) < 0.05,
        f"Composite PDF dimensions are {width_mm:.2f} x {height_mm:.2f} mm",
    )

    text = subprocess.run(
        [pdftotext, "-layout", str(OUT_BASE.with_suffix(".pdf")), "-"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    required_titles = [
        "Donor-paired chamber estimand",
        "Healthy-reference direction and magnitude",
        "Healthy-only selection",
        "Frozen NCP program-to-family topology",
    ]
    require(all(title in text for title in required_titles), "Composite PDF is missing a panel title")
    require(
        not any(token in text for token in ("NCP_v3", "results/", "scripts/")),
        "Composite PDF contains forbidden publication-facing text",
    )

    composite_svg = OUT_BASE.with_suffix(".svg").read_text(encoding="utf-8")
    require(
        composite_svg.count("data:image/svg+xml;base64,") == 4,
        "Composite SVG does not contain four self-contained vector panels",
    )
    require(
        not re.search(r'href="(?:file:|/)', composite_svg),
        "Composite SVG contains an external file reference",
    )


def write_manifest_and_qc() -> None:
    rows = []
    for panel in PANELS:
        rows.append(
            {
                "asset_role": f"panel_{panel['panel']}",
                "asset_path": str(panel["pdf"].relative_to(ROOT)),
                "sha256": sha256(panel["pdf"]),
                "x_mm": panel["x_mm"],
                "y_mm": panel["y_mm"],
                "width_mm": panel["width_mm"],
                "height_mm": panel["height_mm"],
            }
        )
    for path in FROZEN_INPUTS:
        rows.append(
            {
                "asset_role": "frozen_source",
                "asset_path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "x_mm": "",
                "y_mm": "",
                "width_mm": "",
                "height_mm": "",
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(OUT_MANIFEST, sep="\t", index=False)

    output_lines = "\n".join(
        f"- `{path.relative_to(ROOT)}` — `{sha256(path)}`"
        for path in [
            OUT_BASE.with_suffix(".pdf"),
            OUT_BASE.with_suffix(".svg"),
            OUT_BASE.with_suffix(".png"),
            OUT_PREVIEW,
        ]
    )
    input_lines = "\n".join(
        f"- {row.asset_role}: `{row.asset_path}` — `{row.sha256}`"
        for row in manifest.itertuples(index=False)
    )
    OUT_QC.write_text(
        f"""# Figure 1 composite QC

Status: **PASS**

## Assembly geometry
- final size: {PAGE_WIDTH_MM} x {PAGE_HEIGHT_MM} mm
- outer margins: {OUTER_MARGIN_MM} mm
- horizontal and vertical gutters: {GUTTER_MM} mm
- columns: {LEFT_WIDTH_MM} / {RIGHT_WIDTH_MM} mm (0.40 / 0.60)
- rows: {TOP_HEIGHT_MM} / {BOTTOM_HEIGHT_MM} mm
- reading order: A -> B -> C -> D
- vector assembly: panel PDFs into composite PDF; self-contained panel SVGs into composite SVG
- raster export: 600 dpi
- intended manuscript width: {PAGE_WIDTH_MM} mm
- font family: DejaVu Sans
- minimum publication-facing text size: 6.5 pt
- PDF panel titles and forbidden visible-path text: checked
- SVG external file references: none

## Frozen-content assertions
- programs: 30
- program direction composition: 22 LA / 8 RA
- selection lineage: 80 -> 38 -> 25 strict + 5 fallback -> 30
- fixed topology: 30 programs -> 16 families
- Figure 1B and 1D canonical program order: identical
- canonical direction colors: LA {LA_COLOR}; RA {RA_COLOR}
- panel source data or frozen identifiers modified by assembler: no

## Inputs
{input_lines}

## Outputs
{output_lines}
""",
        encoding="utf-8",
    )


def main() -> None:
    OUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    validate_geometry()
    validate_frozen_content()
    validate_presentation_sources()
    write_svg()
    write_pdf()
    rasterize()
    validate_outputs()
    write_manifest_and_qc()
    print(f"Wrote {OUT_BASE.relative_to(ROOT)}.[pdf|svg|png]")
    print(f"Wrote {OUT_PREVIEW.relative_to(ROOT)}")
    print(f"Wrote {OUT_MANIFEST.relative_to(ROOT)}")
    print(f"Wrote {OUT_QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
