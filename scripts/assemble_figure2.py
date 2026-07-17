#!/usr/bin/env python3
"""Assemble the corrected Figure 2A-D vector panels into one publication figure.

Layout (reading order A -> B -> C|D):

  A  full width (178 x 55 mm)
  B  full width (178 x 105 mm)
  C | D  side-by-side composite variants (88 x 134 mm each), shared row contract
"""
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
NCP_FIG = FIGURES / "ncp_v3_singscore_reanalysis"
NCP_RES = RESULTS / "ncp_v3_singscore_reanalysis"

OUT_BASE = FIGURES / "figure2_composite"
OUT_PREVIEW = FIGURES / "figure2_composite_preview.png"
OUT_MANIFEST = RESULTS / "figure2_composite_manifest.tsv"
OUT_QC = RESULTS / "figure2_composite_qc.md"

PAGE_WIDTH_MM = 183.0
GUTTER_MM = 4.5
PANEL_A_W_MM = 178.0
PANEL_A_H_MM = 55.0
PANEL_B_W_MM = 178.0
PANEL_B_H_MM = 105.0
PANEL_CD_W_MM = 88.0
PANEL_CD_H_MM = 134.0
CD_ROW_W_MM = PANEL_CD_W_MM + GUTTER_MM + PANEL_CD_W_MM  # 180.5

OUTER_MARGIN_X_AB_MM = (PAGE_WIDTH_MM - PANEL_A_W_MM) / 2.0  # 2.5
OUTER_MARGIN_X_CD_MM = (PAGE_WIDTH_MM - CD_ROW_W_MM) / 2.0  # 1.25
OUTER_MARGIN_Y_MM = 4.0

PAGE_HEIGHT_MM = (
    OUTER_MARGIN_Y_MM
    + PANEL_A_H_MM
    + GUTTER_MM
    + PANEL_B_H_MM
    + GUTTER_MM
    + PANEL_CD_H_MM
    + OUTER_MARGIN_Y_MM
)

Y_A_MM = OUTER_MARGIN_Y_MM
Y_B_MM = Y_A_MM + PANEL_A_H_MM + GUTTER_MM
Y_CD_MM = Y_B_MM + PANEL_B_H_MM + GUTTER_MM

LA_COLOR = "#C97A74"
RA_COLOR = "#7D99B8"

PANELS = [
    {
        "panel": "A",
        "pdf": NCP_FIG / "fig2a_frozen_ncp_bulk_deployment.pdf",
        "svg": NCP_FIG / "fig2a_frozen_ncp_bulk_deployment.svg",
        "x_mm": OUTER_MARGIN_X_AB_MM,
        "y_mm": Y_A_MM,
        "width_mm": PANEL_A_W_MM,
        "height_mm": PANEL_A_H_MM,
    },
    {
        "panel": "B",
        "pdf": FIGURES / "fig2b_donor_acai_heatmap.pdf",
        "svg": FIGURES / "fig2b_donor_acai_heatmap.svg",
        "x_mm": OUTER_MARGIN_X_AB_MM,
        "y_mm": Y_B_MM,
        "width_mm": PANEL_B_W_MM,
        "height_mm": PANEL_B_H_MM,
    },
    {
        "panel": "C",
        "pdf": NCP_FIG / "fig2c_family_cohort_state_map_composite.pdf",
        "svg": NCP_FIG / "fig2c_family_cohort_state_map_composite.svg",
        "x_mm": OUTER_MARGIN_X_CD_MM,
        "y_mm": Y_CD_MM,
        "width_mm": PANEL_CD_W_MM,
        "height_mm": PANEL_CD_H_MM,
    },
    {
        "panel": "D",
        "pdf": NCP_FIG / "fig2d_family_effects_consistency_composite.pdf",
        "svg": NCP_FIG / "fig2d_family_effects_consistency_composite.svg",
        "x_mm": OUTER_MARGIN_X_CD_MM + PANEL_CD_W_MM + GUTTER_MM,
        "y_mm": Y_CD_MM,
        "width_mm": PANEL_CD_W_MM,
        "height_mm": PANEL_CD_H_MM,
    },
]

COMPOSITE_GEOMETRY = NCP_RES / "figure2_family_row_geometry_composite.tsv"
PANEL_RENDERERS = [
    ROOT / "scripts" / "ncp_v3_singscore_reanalysis" / "fig2a_bulk_deployment_bridge.py",
    ROOT / "scripts" / "ncp_v3_singscore_reanalysis" / "fig2b_family_donor_heatmap.py",
    ROOT / "scripts" / "ncp_v3_singscore_reanalysis" / "fig2c_family_cohort_state_map.py",
    ROOT / "scripts" / "ncp_v3_singscore_reanalysis" / "fig2d_family_effects_consistency.py",
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
        abs(PAGE_HEIGHT_MM - (
            OUTER_MARGIN_Y_MM
            + PANEL_A_H_MM
            + GUTTER_MM
            + PANEL_B_H_MM
            + GUTTER_MM
            + PANEL_CD_H_MM
            + OUTER_MARGIN_Y_MM
        ))
        < 1e-9,
        "Figure 2 vertical geometry does not close",
    )
    require(
        abs(CD_ROW_W_MM - (PANEL_CD_W_MM + GUTTER_MM + PANEL_CD_W_MM)) < 1e-9,
        "Figure 2 C|D row width does not close",
    )
    require(PAGE_WIDTH_MM >= CD_ROW_W_MM, "Page narrower than C|D row")
    require(PAGE_WIDTH_MM >= PANEL_A_W_MM, "Page narrower than A/B panels")

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


def validate_cd_row_contract() -> None:
    require(COMPOSITE_GEOMETRY.is_file(), f"Missing composite geometry: {COMPOSITE_GEOMETRY}")
    geometry = pd.read_csv(COMPOSITE_GEOMETRY, sep="\t")
    require(len(geometry) == 16, "Composite geometry is not 16 rows")
    require(
        abs(float(geometry["row_height"].iloc[0]) - 5.60) < 1e-6,
        "Composite row height is not 5.60 mm",
    )
    require(
        abs(float(geometry["matrix_body_top"].iloc[0]) - 27.0) < 1e-6,
        "Composite TABLE_TOP is not 27.0 mm",
    )
    require(
        geometry["direction_block"].tolist() == ["LA"] * 11 + ["RA"] * 5,
        "Composite LA/RA blocks changed",
    )


def validate_presentation_sources() -> None:
    for path in PANEL_RENDERERS:
        text = path.read_text(encoding="utf-8")
        require("DejaVu Sans" in text, f"Shared DejaVu Sans font is not declared in {path}")
    for path in PANEL_RENDERERS[2:]:
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
    with tempfile.TemporaryDirectory(prefix="figure2-composite-") as tmp:
        tmp_path = Path(tmp)
        tex_path = tmp_path / "figure2_composite.tex"
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
        shutil.copy2(tmp_path / "figure2_composite.pdf", OUT_BASE.with_suffix(".pdf"))


def rasterize() -> None:
    pdftocairo = shutil.which("pdftocairo")
    require(pdftocairo is not None, "pdftocairo is required for Figure 2 PNG exports")
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
        "Cohort-specific family rewiring states",
        "Cross-cohort AF-SR effects and consistency",
        "Hedges' g (SMD)",
    ]
    # Panel A/B titles vary; require at least C/D collage markers and cohort IDs.
    require(all(title in text for title in required_titles), "Composite PDF is missing a panel title")
    require("GSE128188" in text and "GSE79768" in text and "GSE115574" in text, "Cohort IDs missing")
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
    rows.append(
        {
            "asset_role": "composite_geometry",
            "asset_path": str(COMPOSITE_GEOMETRY.relative_to(ROOT)),
            "sha256": sha256(COMPOSITE_GEOMETRY),
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
        f"""# Figure 2 composite QC

Status: **PASS**

## Assembly geometry
- final size: {PAGE_WIDTH_MM} x {PAGE_HEIGHT_MM} mm
- vertical margins: {OUTER_MARGIN_Y_MM} mm
- horizontal gutters: {GUTTER_MM} mm
- A/B width: {PANEL_A_W_MM} mm (centered with {OUTER_MARGIN_X_AB_MM} mm side margins)
- C|D: {PANEL_CD_W_MM} + {GUTTER_MM} + {PANEL_CD_W_MM} = {CD_ROW_W_MM} mm (centered with {OUTER_MARGIN_X_CD_MM} mm side margins)
- row heights: A {PANEL_A_H_MM} / B {PANEL_B_H_MM} / C|D {PANEL_CD_H_MM} mm
- reading order: A -> B -> C|D
- C and D use composite panel variants with shared row contract (row 5.60 mm; TABLE_TOP 27.0 mm; panel 134 mm)
- vector assembly: panel PDFs into composite PDF; self-contained panel SVGs into composite SVG
- raster export: 600 dpi
- intended manuscript width: {PAGE_WIDTH_MM} mm
- font family: DejaVu Sans
- PDF panel titles and forbidden visible-path text: checked
- SVG external file references: none

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
    validate_cd_row_contract()
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
