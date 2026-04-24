"""Shared docx styling helper for the manuscript_drafts/*.docx artifacts."""
from __future__ import annotations

from docx.document import Document as _DocType
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_LINE_SPACING, WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Cm, Pt, RGBColor


def apply_base_style(doc: _DocType) -> None:
    """Apply NC-style page margins, fonts, line spacing, and heading styles."""
    for section in doc.sections:
        section.top_margin    = Cm(2.5)
        section.bottom_margin = Cm(2.5)
        section.left_margin   = Cm(2.5)
        section.right_margin  = Cm(2.5)

    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(11)
    normal.font.color.rgb = RGBColor(0x10, 0x10, 0x10)
    pf = normal.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    pf.space_before = Pt(0)
    pf.space_after  = Pt(6)

    h1 = doc.styles["Heading 1"]
    h1.font.name = "Arial"
    h1.font.size = Pt(13)
    h1.font.bold = True
    h1.font.color.rgb = RGBColor(0x1A, 0x1A, 0x1A)

    h2 = doc.styles["Heading 2"]
    h2.font.name = "Arial"
    h2.font.size = Pt(11.5)
    h2.font.bold = True
    h2.font.color.rgb = RGBColor(0x1A, 0x1A, 0x1A)
    h2.paragraph_format.space_before = Pt(14)
    h2.paragraph_format.space_after  = Pt(4)

    h3 = doc.styles["Heading 3"]
    h3.font.name = "Arial"
    h3.font.size = Pt(11)
    h3.font.bold = True
    h3.font.color.rgb = RGBColor(0x1A, 0x1A, 0x1A)
    h3.paragraph_format.space_before = Pt(10)
    h3.paragraph_format.space_after  = Pt(2)


def add_body(doc: _DocType, text: str, *, italic: bool = False) -> None:
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    run = p.add_run(text)
    if italic:
        run.italic = True


def add_label_run(p, label: str, body: str) -> None:
    """Append `**label** body` style mixed-format inline text to paragraph p."""
    run_lab = p.add_run(label)
    run_lab.bold = True
    p.add_run(body)


def add_caption_paragraph(doc: _DocType, panels: list[tuple[str, str]],
                          *, trailing_note: str | None = None) -> None:
    """Compose a figure-caption paragraph from [(panel_letter, body_text), ...].

    Format: **a** body. **b** body. ...
    Optional trailing italic note appended after a single space.
    """
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    for i, (lab, body) in enumerate(panels):
        if i:
            p.add_run("  ")
        run_l = p.add_run(lab)
        run_l.bold = True
        p.add_run(" " + body)
    if trailing_note:
        p.add_run("  ")
        run_n = p.add_run(trailing_note)
        run_n.italic = True


def add_table(doc: _DocType, rows: list[list[str]],
              *, col_widths_cm: list[float] | None = None) -> None:
    """Create a simple grid-style table with first row as bold header."""
    n_cols = len(rows[0])
    table = doc.add_table(rows=len(rows), cols=n_cols)
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER

    if col_widths_cm:
        for col_idx, w in enumerate(col_widths_cm):
            for cell in table.columns[col_idx].cells:
                cell.width = Cm(w)

    for r_idx, row in enumerate(rows):
        for c_idx, val in enumerate(row):
            cell = table.cell(r_idx, c_idx)
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            p.paragraph_format.space_after = Pt(0)
            run = p.add_run(val)
            run.font.size = Pt(9)
            if r_idx == 0:
                run.bold = True

    # spacing after table
    doc.add_paragraph()
