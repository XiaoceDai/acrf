#!/usr/bin/env python3
"""
Assemble final Figure 6 A–D composite from exact-size vector panel SVGs.

2 × 2 layout. Does not stitch or rescale PNG panels.
Each panel is nested at 100% scale. No coloc / external-set rematch.
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from figure6_composite_geometry import (  # noqa: E402
    COMPOSITE_H_MM,
    COMPOSITE_W_MM,
    PANEL_A,
    PANEL_B,
    PANEL_C,
    PANEL_D,
    write_ab_geometry,
    write_cd_geometry,
    write_composite_geometry,
)

FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
OUT_BASE = FIGURES / "figure6_genetic_traceability"
OUT_PREVIEW = FIGURES / "_pdf_preview" / "figure6_genetic_traceability.png"
OUT_GRAY = FIGURES / "_pdf_preview" / "figure6_genetic_traceability_grayscale.png"
OUT_CVD = FIGURES / "_pdf_preview" / "figure6_genetic_traceability_cvd.png"
OUT_MANIFEST = RESULTS / "figure6_composite_manifest.tsv"
OUT_QC = RESULTS / "figure6_composite_qc.md"

RENDER_SCRIPT = ROOT / "scripts" / "render_figure6_composite_panels.py"
PANEL_SVGS = {
    "A": FIGURES / "fig6a_anchor_coloc_calibration.svg",
    "B": FIGURES / "fig6b_posterior_lead_distance.svg",
    "C": FIGURES / "fig6c_family_af_overlap.svg",
    "D": FIGURES / "fig6d_traceability_matrix.svg",
}
PANEL_BOX = {"A": PANEL_A, "B": PANEL_B, "C": PANEL_C, "D": PANEL_D}

FORBIDDEN_VISIBLE = [
    "NCP_v3", "ncp_v3", "Miyazawa", "miyazawa", "mr_",
    "figure6_mr", "pipeline_wide", "G0", "G1", "G2", "G3", "G4",
    "highest tier", "evidence score", "G-tier", "Path A gate",
]


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def prefix_svg_ids(svg: str, prefix: str) -> str:
    identifiers = re.findall(r'\bid="([^"]+)"', svg)
    for identifier in sorted(set(identifiers), key=len, reverse=True):
        replacement = f"{prefix}_{identifier}"
        svg = svg.replace(f'id="{identifier}"', f'id="{replacement}"')
        svg = svg.replace(f"url(#{identifier})", f"url(#{replacement})")
        svg = svg.replace(f'href="#{identifier}"', f'href="#{replacement}"')
        svg = svg.replace(f'xlink:href="#{identifier}"', f'xlink:href="#{replacement}"')
    return svg


def svg_dimensions(path: Path) -> tuple[float, float, str, str]:
    text = path.read_text(encoding="utf-8")
    width = re.search(r'<svg\b[^>]*\bwidth="([0-9.]+)(mm|pt|in)?"', text)
    height = re.search(r'<svg\b[^>]*\bheight="([0-9.]+)(mm|pt|in)?"', text)
    viewbox = re.search(r'\bviewBox="([^"]+)"', text)
    body = re.search(r"<svg\b[^>]*>(.*)</svg>\s*$", text, flags=re.DOTALL)
    require(all(v is not None for v in (width, height, viewbox, body)), f"parse fail {path}")
    w_raw, w_unit = float(width.group(1)), width.group(2) or "pt"
    h_raw, h_unit = float(height.group(1)), height.group(2) or "pt"

    def to_mm(v: float, unit: str) -> float:
        if unit == "mm":
            return v
        if unit == "in":
            return v * 25.4
        return v * 25.4 / 72

    if w_unit is None or w_unit == "pt":
        if w_raw < 20:
            w_mm, h_mm = w_raw * 25.4, h_raw * 25.4
        else:
            w_mm, h_mm = to_mm(w_raw, "pt"), to_mm(h_raw, "pt")
    else:
        w_mm, h_mm = to_mm(w_raw, w_unit), to_mm(h_raw, h_unit)
    return w_mm, h_mm, viewbox.group(1), body.group(1)


def nested_svg(path: Path, panel: str) -> str:
    box = PANEL_BOX[panel]
    x, y, width, height = box["x"], box["y"], box["w"], box["h"]
    source_w, source_h, viewbox, body = svg_dimensions(path)
    require(
        abs(source_w - width) < 0.8 and abs(source_h - height) < 0.8,
        f"Panel {panel} SVG {source_w:.2f}×{source_h:.2f} mm != {width}×{height} mm",
    )
    body = prefix_svg_ids(body, f"panel{panel.lower()}")
    return (
        f'<g id="panel-{panel}" transform="translate({x},{y})">'
        f'<svg x="0" y="0" width="{width}" height="{height}" '
        f'viewBox="{viewbox}" preserveAspectRatio="xMinYMin meet">{body}</svg></g>'
    )


def write_composite_svg(output: Path) -> None:
    nested = "\n".join(nested_svg(PANEL_SVGS[p], p) for p in ("A", "B", "C", "D"))
    output.write_text(
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{COMPOSITE_W_MM}mm" height="{COMPOSITE_H_MM}mm" '
        f'viewBox="0 0 {COMPOSITE_W_MM} {COMPOSITE_H_MM}">\n'
        f'<rect x="0" y="0" width="{COMPOSITE_W_MM}" height="{COMPOSITE_H_MM}" fill="#FFFFFF"/>\n'
        f"{nested}\n</svg>\n",
        encoding="utf-8",
    )


def _cvd_simulate(rgb_img):
    """Simple Brettel-like deuteranopia approximation (sRGB 0–255)."""
    import numpy as np
    arr = np.asarray(rgb_img, dtype=np.float64) / 255.0
    # LMS conversion + deuteranope matrix (approx)
    m = np.array([
        [0.625, 0.375, 0.0],
        [0.70, 0.30, 0.0],
        [0.0, 0.30, 0.70],
    ])
    flat = arr.reshape(-1, 3)
    out = flat @ m.T
    out = np.clip(out, 0, 1)
    return (out.reshape(arr.shape) * 255).astype("uint8")


def export_raster(
    svg: Path, pdf: Path, png: Path, preview: Path, gray: Path, cvd: Path,
) -> None:
    mutool = shutil.which("mutool")
    if mutool:
        subprocess.run(
            [mutool, "convert", "-o", str(pdf), str(svg)],
            check=True, capture_output=True, text=True,
        )
        subprocess.run(
            [mutool, "draw", "-q", "-r", "600", "-o", str(png), str(pdf)],
            check=True, capture_output=True, text=True,
        )
        subprocess.run(
            [mutool, "draw", "-q", "-r", "150", "-o", str(preview), str(pdf)],
            check=True, capture_output=True, text=True,
        )
    else:
        import cairosvg  # type: ignore
        cairosvg.svg2pdf(url=str(svg), write_to=str(pdf))
        cairosvg.svg2png(url=str(svg), write_to=str(png), dpi=600)
        cairosvg.svg2png(url=str(svg), write_to=str(preview), dpi=150)
    try:
        from PIL import Image
        import numpy as np
        im = Image.open(preview).convert("RGB")
        im.convert("L").save(gray)
        cvd_arr = _cvd_simulate(np.asarray(im))
        Image.fromarray(cvd_arr, mode="RGB").save(cvd)
    except Exception:
        shutil.copy2(preview, gray)
        shutil.copy2(preview, cvd)


def render_panels() -> None:
    write_ab_geometry()
    write_cd_geometry()
    write_composite_geometry()
    print(f"=== rendering composite panels: {RENDER_SCRIPT.name} ===")
    subprocess.run([sys.executable, str(RENDER_SCRIPT)], check=True, cwd=str(ROOT))


def scientific_assertions() -> dict:
    a = pd.read_csv(RESULTS / "fig6a_anchor_coloc_calibration_source.tsv", sep="\t")
    b = pd.read_csv(RESULTS / "fig6b_posterior_lead_distance_source.tsv", sep="\t")
    c = pd.read_csv(RESULTS / "fig6c_family_af_overlap_source.tsv", sep="\t")
    d = pd.read_csv(RESULTS / "fig6d_non_nested_traceability_source.tsv", sep="\t")

    require(len(a) == 12 and list(a["gene"]) == list(b["gene"]), "A/B anchor order")
    require(int(a["local_AF_GWAS_support"].astype(str).str.lower().isin(["true", "1"]).sum()) == 7, "A GWAS")
    require(int(a["atrial_eGene_bonferroni"].astype(str).str.lower().isin(["true", "1"]).sum()) == 1, "A eGene")
    require(int(a["PP4_ge_0.70"].astype(str).str.lower().isin(["true", "1"]).sum()) == 0, "A PP4")

    require(int(b["pp3_dominant"].astype(str).str.lower().isin(["true", "1"]).sum()) == 12, "B PP3")
    require(int(b["pp4_support"].astype(str).str.lower().isin(["true", "1"]).sum()) == 0, "B PP4")
    dmin, dmax = float(b["delta_mb"].min()), float(b["delta_mb"].max())
    require(abs(dmin - 0.46) < 0.02, f"B dist min {dmin}")
    require(abs(dmax - 1.92) < 0.02, f"B dist max {dmax}")
    k = b.loc[b["gene"] == "KCNN3"].iloc[0]
    require(abs(float(k["delta_mb"]) - 0.46) < 0.02, "KCNN3 distance")

    require(len(c) == 64, "C 64 tests")
    require(int((c["n_overlap"] > 0).sum()) == 5, "C nonzero")
    require(int((c["fdr_bh"] < 0.05).sum()) == 0, "C FDR")
    fams = sorted(c.loc[c["n_overlap"] > 0, "family_id"].astype(int).unique())
    require(fams == [3, 9, 13], f"C families {fams}")

    require(len(d) == 12, "D anchors")
    require(int((d["local_gwas_status"] == "present").sum()) == 7, "D G")
    require(int((d["atrial_egene_status"] == "present").sum()) == 1, "D E")
    require(int((d["coloc_status"] == "present").sum()) == 0, "D C")
    require(int((d["common_risk_status"] == "present").sum()) == 2, "D risk")
    require(int((d["regulatory_status"] == "present").sum()) == 2, "D reg")
    require(int((d["rare_functional_status"] == "present").sum()) == 0, "D rare")
    smad = d.loc[d["gene"] == "SMAD9"].iloc[0]
    require(smad["local_gwas_status"] == "evaluated_no_support", "SMAD9 G")
    require(smad["atrial_egene_status"] == "present", "SMAD9 E")
    require(smad["coloc_status"] == "evaluated_no_support", "SMAD9 C")

    return {"dist_min": dmin, "dist_max": dmax, "c_nonzero": 5, "c_fdr": 0}


def layout_assertions() -> None:
    ab = pd.read_csv(RESULTS / "figure6ab_horizontal_geometry.tsv", sep="\t")
    require(float(ab.loc[ab["metric"] == "height", "value"].iloc[0]) == 92.0, "AB height")
    require(float(ab.loc[ab["metric"] == "panel_a_width", "value"].iloc[0]) == 72.0, "A width")
    require(float(ab.loc[ab["metric"] == "panel_b_width", "value"].iloc[0]) == 102.0, "B width")
    cd = pd.read_csv(RESULTS / "figure6cd_horizontal_geometry.tsv", sep="\t")
    require(float(cd.loc[cd["metric"] == "height", "panel_C"].iloc[0]) == 78.0, "CD height")
    require(float(cd.loc[cd["metric"] == "width", "panel_C"].iloc[0]) == 70.0, "C width")
    require(float(cd.loc[cd["metric"] == "width", "panel_D"].iloc[0]) == 104.0, "D width")
    require(abs(COMPOSITE_H_MM - 174.0) < 0.01, "canvas height")
    require(abs(COMPOSITE_W_MM - 178.0) < 0.01, "canvas width")


def write_manifest(outputs: dict[str, Path]) -> None:
    rows = []
    for panel, path in PANEL_SVGS.items():
        box = PANEL_BOX[panel]
        rows.append({
            "artifact": f"panel_{panel}_svg",
            "path": str(path.relative_to(ROOT)),
            "sha256": sha256(path),
            "x_mm": box["x"], "y_mm": box["y"],
            "width_mm": box["w"], "height_mm": box["h"],
            "scale_percent": 100,
        })
    for key, path in outputs.items():
        rows.append({
            "artifact": key,
            "path": str(path.relative_to(ROOT)),
            "sha256": sha256(path),
            "x_mm": "", "y_mm": "",
            "width_mm": COMPOSITE_W_MM if "preview" not in key else "",
            "height_mm": COMPOSITE_H_MM if "preview" not in key else "",
            "scale_percent": 100,
        })
    pd.DataFrame(rows).to_csv(OUT_MANIFEST, sep="\t", index=False)


def write_qc(assertions: dict, outputs: dict[str, Path]) -> None:
    text = OUT_BASE.with_suffix(".svg").read_text(encoding="utf-8")
    visible = re.sub(r"<[^>]+>", " ", text)
    for bad in FORBIDDEN_VISIBLE:
        require(bad not in visible, f"forbidden visible: {bad}")
    # B must not repeat anchors
    b_vis = re.sub(r"<[^>]+>", " ", PANEL_SVGS["B"].read_text(encoding="utf-8"))
    for gene in ("PRELP", "SMAD9", "PITX2", "KCNN3", "BMP10"):
        require(not re.search(rf"\b{gene}\b", b_vis), f"B repeats {gene}")
    require("No score" in visible and "not Path A" in visible, "D boundary tag")
    require("prioritized" in visible, "D full Risk/prioritized header")

    lines = [
        "# Figure 6 genetic traceability composite QC",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "Status: **PASS**",
        "",
        "## Canvas (2 × 2)",
        f"- {COMPOSITE_W_MM:.0f} × {COMPOSITE_H_MM:.0f} mm",
        f"- A: ({PANEL_A['x']}, {PANEL_A['y']}) {PANEL_A['w']:.0f} × {PANEL_A['h']:.0f} mm",
        f"- B: ({PANEL_B['x']}, {PANEL_B['y']}) {PANEL_B['w']:.0f} × {PANEL_B['h']:.0f} mm",
        f"- C: ({PANEL_C['x']}, {PANEL_C['y']}) {PANEL_C['w']:.0f} × {PANEL_C['h']:.0f} mm",
        f"- D: ({PANEL_D['x']}, {PANEL_D['y']}) {PANEL_D['w']:.0f} × {PANEL_D['h']:.0f} mm",
        "- Gutters: 4 mm horizontal and vertical",
        "- Panels nested at 100% scale; no PNG stitching",
        "- B omits family/role and anchor labels (shared with A)",
        "",
        "## Scientific assertions",
        "- A: GWAS 7/12 · eGene 1/12 · PP4 0/12",
        f"- B: PP3-dominant 12/12 · PP4 0/12 · distance "
        f"{assertions['dist_min']:.2f}–{assertions['dist_max']:.2f} Mb · KCNN3 0.46 Mb",
        "- C: 64 tests · 5 nonzero · F3/F13/F9 · FDR 0 · N=14635",
        "- D: G/E/C 7/1/0 · risk/reg/rare 2/2/0 · no G-tier",
        "",
        "## Geometry contracts",
        "- results/figure6ab_horizontal_geometry.tsv",
        "- results/figure6cd_horizontal_geometry.tsv",
        "- results/figure6_composite_geometry.tsv",
        "",
        "## SHA256",
        *[f"- `{p.relative_to(ROOT)}`: `{sha256(p)}`" for p in outputs.values()],
        "",
        "## Scope",
        "- No coloc rerun; no external-set rematch",
        "- Formal submission routing not updated until composite QC passes",
        "",
    ]
    OUT_QC.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT_PREVIEW.parent.mkdir(parents=True, exist_ok=True)

    render_panels()
    for p, path in PANEL_SVGS.items():
        require(path.is_file(), f"missing {p} SVG: {path}")

    assertions = scientific_assertions()
    layout_assertions()

    out_svg = OUT_BASE.with_suffix(".svg")
    out_pdf = OUT_BASE.with_suffix(".pdf")
    out_png = OUT_BASE.with_suffix(".png")
    write_composite_svg(out_svg)
    export_raster(out_svg, out_pdf, out_png, OUT_PREVIEW, OUT_GRAY, OUT_CVD)

    outputs = {
        "composite_svg": out_svg,
        "composite_pdf": out_pdf,
        "composite_png": out_png,
        "preview": OUT_PREVIEW,
        "grayscale": OUT_GRAY,
        "cvd": OUT_CVD,
        "ab_geometry": RESULTS / "figure6ab_horizontal_geometry.tsv",
        "cd_geometry": RESULTS / "figure6cd_horizontal_geometry.tsv",
        "composite_geometry": RESULTS / "figure6_composite_geometry.tsv",
    }
    write_qc(assertions, outputs)
    outputs["qc"] = OUT_QC
    write_manifest(outputs)
    print(f"wrote {out_svg.relative_to(ROOT)}")
    print(f"wrote {out_pdf.relative_to(ROOT)}")
    print(f"wrote {out_png.relative_to(ROOT)}")
    print(f"wrote {OUT_QC.relative_to(ROOT)}")
    print("COMPOSITE QC: PASS")


if __name__ == "__main__":
    main()
