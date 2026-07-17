#!/usr/bin/env python3
"""
Audit figure layout / readability for all figures/*.pdf.

Outputs
-------
results/figure_layout_audit.tsv
    one row per PDF with driver mapping, static risk flags, min observed fontsize
figures/_audit_png/
    PNG renderings of every PDF for visual spot-checking
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"
OUT = ROOT / "results"
AUDIT_PNG = FIG / "_audit_png"
OUT_TSV = OUT / "figure_layout_audit.tsv"

MAIN: list[tuple[str, Path]] = []
SUPP: list[tuple[str, Path]] = []

# import submission lists for driver mapping
import importlib.util
spec = importlib.util.spec_from_file_location("build_plos_submission_folder", str(ROOT / "scripts" / "build_plos_submission_folder.py"))
if spec and spec.loader:
    bsp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bsp)
    for label, p in getattr(bsp, "MAIN_FIGURES", []):
        MAIN.append((label, p))
    for label, p in getattr(bsp, "SUPPORTING_FIGURES", []):
        SUPP.append((label, p))

SUBMISSION_PDFS = {p.name for _, p in MAIN + SUPP}

# Read reproducibility table for driver mapping
REPRO = {}
repro_tsv = OUT / "reproducibility_table.tsv"
if repro_tsv.is_file():
    import pandas as pd
    rdf = pd.read_csv(repro_tsv, sep="\t", dtype=str)
    for _, row in rdf.iterrows():
        for part in str(row["output"]).split(" + "):
            if part.endswith(".pdf"):
                REPRO[Path(part).name] = str(row["script"])

# Fontsize regex for python scripts
FONT_RE = re.compile(r"fontsize\s*=\s*(\d+(?:\.\d+)?)")
FIG_TEXT_RE = re.compile(r"fig\.text\s*\(\s*0\.(\d+)")
SUPTITLE_RE = re.compile(r"fig\.suptitle\s*\(|")
GRIDS_RE = re.compile(r"top\s*=\s*(0\.\d+)")


def find_pdf_driver(pdf: Path) -> str | None:
    # 1) reproducibility table
    name = pdf.name
    if name in REPRO:
        return REPRO[name]
    # 2) grep figures/ references in scripts
    for script in sorted((ROOT / "scripts").rglob("*.py")):
        if script.name == "audit_figure_layout.py":
            continue
        text = script.read_text(errors="ignore")
        if f'"{name}"' in text or f"'/{name}'" in text or f'/{name}"' in text or f'/{name}\'' in text:
            if f"fig.savefig" in text or "OUT_PDF" in text or "FIG_OUT" in text:
                return str(script.relative_to(ROOT))
    return None


def render_png(pdf: Path) -> Path:
    AUDIT_PNG.mkdir(parents=True, exist_ok=True)
    out_png = AUDIT_PNG / f"{pdf.stem}.png"
    if not out_png.is_file():
        subprocess.run(
            ["pdftoppm", "-png", "-r", "150", "-singlefile", str(pdf), str(out_png.with_suffix(""))],
            check=True, capture_output=True,
        )
    return out_png


def parse_pdf_fonts(pdf: Path) -> dict[str, any]:
    """Use pdftotext -bbox to estimate character positions and detect text overlap."""
    info = {"n_chars": 0, "overlap_candidates": 0, "overlap_zones": []}
    try:
        proc = subprocess.run(
            ["pdftotext", "-bbox", "-bbox-layout", str(pdf), "-"],
            check=True, capture_output=True, text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return info
    xml = proc.stdout
    # xMin/yMin/xMax/yMax are in page-space points; gather all textlines
    lines = re.findall(r'<textline[^>]*>(.*?)</textline>', xml, flags=re.DOTALL)
    boxes: list[tuple[float, float, float, float, str]] = []
    for line in lines:
        for m in re.finditer(r'<word\s+xMin="([0-9.]+)"\s+yMin="([0-9.]+)"\s+xMax="([0-9.]+)"\s+yMax="([0-9.]+)">([^<]+)</word>', line):
            x0, y0, x1, y1, txt = float(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4)), m.group(5)
            if txt.strip():
                boxes.append((x0, y0, x1, y1, txt.strip()))
    info["n_chars"] = len(boxes)
    # Simple overlap: boxes that overlap in both x and y by >20% of smaller area and are not substrings
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            x0a, y0a, x1a, y1a, ta = boxes[i]
            x0b, y0b, x1b, y1b, tb = boxes[j]
            iw = max(0.0, min(x1a, x1b) - max(x0a, x0b))
            ih = max(0.0, min(y1a, y1b) - max(y0a, y0b))
            if iw <= 0 or ih <= 0:
                continue
            area_a = (x1a - x0a) * (y1a - y0a)
            area_b = (x1b - x0b) * (y1b - y0b)
            inter = iw * ih
            if area_a and area_b and inter / min(area_a, area_b) > 0.20:
                if not (ta in tb or tb in ta):
                    info["overlap_candidates"] += 1
                    info["overlap_zones"].append(f"{ta}|{tb}")
    return info


def static_risk(driver: str | None) -> dict[str, str]:
    """Return risk flags from the driver script."""
    risk = {
        "min_fontsize": "",
        "fig_text_count": "0",
        "suptitle": "no",
        "gridspec_top": "",
        "readability_risk": "",
        "header_overlap_risk": "",
    }
    if not driver:
        return risk
    p = ROOT / driver
    if not p.is_file():
        return risk
    text = p.read_text(errors="ignore")
    sizes = [float(m) for m in FONT_RE.findall(text)]
    risk["min_fontsize"] = f"{min(sizes):.2f}" if sizes else ""
    risk["fig_text_count"] = str(len(FIG_TEXT_RE.findall(text)))
    risk["suptitle"] = "yes" if SUPTITLE_RE.search(text) else "no"
    tops = [float(m) for m in GRIDS_RE.findall(text)]
    risk["gridspec_top"] = f"{min(tops):.3f}" if tops else ""
    # readability risk: sub-6.5 or sub-7 in main-text figures
    if sizes:
        if min(sizes) < 6.5:
            risk["readability_risk"] = "sub_6.5"
        elif min(sizes) < 7.0:
            risk["readability_risk"] = "sub_7"
    # header overlap risk: fig.text near top + gridspec top < 0.90
    near_top = any(float(y) >= 92 for y in FIG_TEXT_RE.findall(text))
    if near_top and tops and min(tops) < 0.90:
        risk["header_overlap_risk"] = "high"
    elif near_top:
        risk["header_overlap_risk"] = "medium"
    return risk


def main() -> None:
    pdfs = sorted(FIG.glob("*.pdf"))
    rows = []
    for pdf in pdfs:
        driver = find_pdf_driver(pdf)
        risk = static_risk(driver)
        label = ""
        if pdf.name in {p.name for _, p in MAIN}:
            label = "main"
        elif pdf.name in {p.name for _, p in SUPP}:
            label = "supporting"
        elif pdf.name.startswith("fig1") or pdf.name.startswith("fig2"):
            label = "legacy"
        else:
            label = "exploratory"
        render_png(pdf)
        overlap = parse_pdf_fonts(pdf)
        overlap_sev = ""
        if overlap["overlap_candidates"] == 0:
            overlap_sev = "none"
        elif overlap["overlap_candidates"] <= 3:
            overlap_sev = "low"
        elif overlap["overlap_candidates"] <= 8:
            overlap_sev = "medium"
        else:
            overlap_sev = "high"
        rows.append({
            "pdf": pdf.name,
            "kind": label,
            "driver_script": driver or "manual_or_legacy",
            "min_fontsize_pt": risk["min_fontsize"],
            "fig_text_count": risk["fig_text_count"],
            "suptitle": risk["suptitle"],
            "gridspec_top": risk["gridspec_top"],
            "readability_risk": risk["readability_risk"],
            "header_overlap_risk": risk["header_overlap_risk"],
            "overlap_severity": overlap_sev,
            "overlap_candidates": overlap["overlap_candidates"],
            "overlap_zones": "; ".join(overlap["overlap_zones"][:10]),
            "fix_action": "",
        })

    import pandas as pd
    df = pd.DataFrame(rows)
    df = df.sort_values(["kind", "pdf"])
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(df)} PDFs)")
    print(f"rendered PNGs to {AUDIT_PNG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
