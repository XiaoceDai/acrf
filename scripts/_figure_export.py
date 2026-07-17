"""Shared PDF + PNG export helpers for figure drivers."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt


def pdf_safe_text(text: str) -> str:
    """Replace glyphs that break lightweight PDF viewers (Cursor remote, etc.)."""
    repl = {
        "²⁺": "2+",
        "β": "beta",
        "–": "-",
        "−": "-",
        "≈": "~",
        "•": "-",
        "·": "/",
        "≥": ">=",
        "≤": "<=",
        "×": "x",
        "—": "-",
    }
    out = str(text)
    for src, dst in repl.items():
        out = out.replace(src, dst)
    return out


def save_figure_bundle(
    fig,
    pdf_path: Path,
    *,
    dpi: int = 300,
    png: bool = True,
    gs_rewrite: bool = True,
    **savefig_kw,
) -> tuple[Path, Path | None]:
    """Save publication PDF and a Cursor-friendly PNG preview."""
    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    defaults = {
        "bbox_inches": "tight",
        "pad_inches": 0.04,
        "facecolor": "white",
        "edgecolor": "none",
    }
    defaults.update(savefig_kw)
    fig.savefig(pdf_path, **defaults)

    png_path: Path | None = None
    if png:
        png_path = pdf_path.with_suffix(".png")
        fig.savefig(png_path, dpi=dpi, **defaults)
        audit_dir = pdf_path.parent / "_audit_png"
        if audit_dir.is_dir():
            shutil.copy2(png_path, audit_dir / png_path.name)

    plt.close(fig)

    if gs_rewrite and _gs_available():
        _rewrite_pdf_gs(pdf_path)

    return pdf_path, png_path


def _gs_available() -> bool:
    try:
        subprocess.run(["gs", "--version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def _rewrite_pdf_gs(pdf_path: Path) -> None:
    tmp = pdf_path.with_suffix(".compat.pdf")
    try:
        subprocess.run(
            [
                "gs", "-q", "-sDEVICE=pdfwrite",
                "-dCompatibilityLevel=1.4",
                "-dNOPAUSE", "-dBATCH",
                f"-sOutputFile={tmp}",
                str(pdf_path),
            ],
            check=True,
            capture_output=True,
        )
        tmp.replace(pdf_path)
    except (subprocess.CalledProcessError, OSError):
        if tmp.is_file():
            tmp.unlink(missing_ok=True)
