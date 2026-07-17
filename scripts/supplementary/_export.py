#!/usr/bin/env python3
"""Shared export helpers for Phase 5B3A canonical SI figures."""
from __future__ import annotations

import hashlib
import re
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SUPP = ROOT / "figures" / "supplementary"
PREVIEW = ROOT / "results" / "manuscript_audit" / "phase5b3a_previews"

FORBIDDEN = [
    "NCP_v3",
    "Normal Chamber Prior version 3",
    "regulator-tier",
    "regulator tier",
    "tier inference",
    "G-tier",
    "G0",
    "G1",
    "G2",
    "G3",
    "G4",
    "Figure 7",
    "Fig. 7",
    "Figure 5E",
    "Fig 5E",
    "Figure 6E",
    "mechanism-direction validation",
    "Atlas PriorityScore",
    "genetic validation",
]


def sha256(path: Path | str) -> str:
    p = Path(path)
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else ""


def ascii_pub_label(text: str) -> str:
    """Map common unicode science glyphs to ASCII for Liberation Sans print safety."""
    return (
        str(text)
        .replace("\u207a", "+")  # ⁺
        .replace("\u207b", "-")  # ⁻
        .replace("\u00b2", "2")  # ²
        .replace("\u00b3", "3")  # ³
        .replace("\u2212", "-")  # −
        .replace("\u2013", "-")  # –
        .replace("\u2014", "-")  # —
    )


def set_pub_style(min_pt: float = 5.5) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
            "font.size": max(min_pt, 7.0),
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def save_fig_bundle(fig: plt.Figure, stem: str, *, pad_inches: float = 0.02) -> dict:
    """Write PDF/SVG/PNG (600 dpi) under figures/supplementary/ and a preview PNG."""
    SUPP.mkdir(parents=True, exist_ok=True)
    PREVIEW.mkdir(parents=True, exist_ok=True)
    pdf = SUPP / f"{stem}.pdf"
    svg = SUPP / f"{stem}.svg"
    png = SUPP / f"{stem}.png"
    prev = PREVIEW / f"{stem.split('_')[0]}_preview.png" if stem.startswith("figS") else PREVIEW / f"{stem}_preview.png"
    # preview names: figS1_preview.png etc.
    m = re.match(r"(figS\d+)_", stem)
    if m:
        prev = PREVIEW / f"{m.group(1)}_preview.png"

    fig.savefig(pdf, format="pdf", bbox_inches="tight", pad_inches=pad_inches)
    fig.savefig(svg, format="svg", bbox_inches="tight", pad_inches=pad_inches)
    fig.savefig(png, format="png", dpi=600, bbox_inches="tight", pad_inches=pad_inches)
    fig.savefig(prev, format="png", dpi=150, bbox_inches="tight", pad_inches=pad_inches)
    plt.close(fig)
    return {"pdf": pdf, "svg": svg, "png": png, "preview": prev}


def export_from_pdf(src_pdf: Path, stem: str) -> dict:
    """Copy/export an existing vector PDF to canonical names without mutating source."""
    SUPP.mkdir(parents=True, exist_ok=True)
    PREVIEW.mkdir(parents=True, exist_ok=True)
    pdf = SUPP / f"{stem}.pdf"
    svg = SUPP / f"{stem}.svg"
    png = SUPP / f"{stem}.png"
    m = re.match(r"(figS\d+)_", stem)
    prev = PREVIEW / (f"{m.group(1)}_preview.png" if m else f"{stem}_preview.png")

    pdf.write_bytes(src_pdf.read_bytes())
    # SVG via poppler (prefix form writes <prefix>.svg)
    svg_prefix = SUPP / stem
    subprocess.check_call(
        ["pdftocairo", "-svg", str(pdf), str(svg_prefix)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    # Some poppler builds write the prefix path with no extension
    bare = SUPP / stem
    if not svg.exists() and bare.exists() and bare.is_file():
        bare.rename(svg)
    if bare.exists() and bare.is_file() and svg.exists() and bare.resolve() != svg.resolve():
        bare.unlink()
    if not svg.exists():
        cands = [p for p in SUPP.glob(f"{stem}*.svg") if p.name != svg.name]
        if cands:
            cands[0].rename(svg)
    subprocess.check_call(
        ["pdftocairo", "-png", "-r", "600", "-singlefile", str(pdf), str(png.with_suffix(""))],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.check_call(
        ["pdftocairo", "-png", "-r", "150", "-singlefile", str(pdf), str(prev.with_suffix(""))],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return {"pdf": pdf, "svg": svg, "png": png, "preview": prev}


def pdf_text(path: Path) -> str:
    try:
        return subprocess.check_output(["pdftotext", "-layout", str(path), "-"], text=True, errors="ignore")
    except Exception:
        return ""


def scan_forbidden(text: str) -> list[str]:
    hits = []
    for tok in FORBIDDEN:
        if tok in ("G0", "G1", "G2", "G3", "G4"):
            if re.search(rf"\b{tok}\b", text):
                hits.append(tok)
        elif tok.lower() in text.lower() if tok.islower() or " " in tok else tok in text:
            # case-sensitive for NCP_v3 etc.; case-insensitive for phrases
            if tok in text or (any(c.islower() for c in tok) and tok.lower() in text.lower()):
                hits.append(tok)
            elif tok not in text and tok.lower() in text.lower() and (" " in tok or tok.islower()):
                hits.append(tok)
    # cleaner pass
    hits = []
    for tok in FORBIDDEN:
        if tok in ("G0", "G1", "G2", "G3", "G4"):
            if re.search(rf"\b{tok}\b", text):
                hits.append(tok)
        elif re.search(re.escape(tok), text, flags=re.IGNORECASE if tok[:3] not in {"NCP", "Fig"} else 0):
            # NCP_v3 exact; phrases case-insensitive
            if tok.startswith("NCP_") or tok.startswith("Figure") or tok.startswith("Fig"):
                if tok in text or (tok.startswith("Fig") and tok in text):
                    hits.append(tok)
                elif tok.lower() in text.lower() and "figure" in tok.lower():
                    hits.append(tok)
            else:
                if tok.lower() in text.lower():
                    hits.append(tok)
        elif tok.lower() in text.lower():
            hits.append(tok)
    return sorted(set(hits))


def page_geometry_mm(pdf: Path) -> tuple[float, float, int]:
    out = subprocess.check_output(["pdfinfo", str(pdf)], text=True, stderr=subprocess.DEVNULL)
    pages = 1
    w = h = float("nan")
    for line in out.splitlines():
        if line.startswith("Pages:"):
            pages = int(line.split(":", 1)[1].strip())
        if line.startswith("Page size:"):
            # e.g. Page size: 538.583 x 793.701 pts
            m = re.search(r"([\d.]+)\s*x\s*([\d.]+)\s*pts", line)
            if m:
                w = float(m.group(1)) * 25.4 / 72.0
                h = float(m.group(2)) * 25.4 / 72.0
    return w, h, pages


def png_size(png: Path) -> tuple[int, int]:
    import struct

    data = Path(png).read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"not a PNG: {png}")
    w, h = struct.unpack(">II", data[16:24])
    return int(w), int(h)
