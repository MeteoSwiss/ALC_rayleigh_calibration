#!/usr/bin/env python3
"""
Minimal Markdown -> PDF for REPORT.md (reportlab + PIL only).

Handles the subset used in the report: #/##/### headings, pipe tables, images
(paths resolved relative to the .md file), bullet lists, horizontal rules,
**bold**, *italic*, `code`, and [links](url). Not a general Markdown engine.

Usage:  python md_to_pdf.py <input.md> [output.pdf]
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, Image, HRFlowable)

# Register full-Unicode fonts (DejaVu ships with matplotlib) so math glyphs
# (×, ², ≲, σ, superscripts, em-dash) render instead of boxes.
import matplotlib
_fd = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
FONT, FONT_B, FONT_I, FONT_M = "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Courier"
try:
    pdfmetrics.registerFont(TTFont("DejaVu", _fd / "DejaVuSans.ttf"))
    pdfmetrics.registerFont(TTFont("DejaVu-Bold", _fd / "DejaVuSans-Bold.ttf"))
    pdfmetrics.registerFont(TTFont("DejaVu-Obl", _fd / "DejaVuSans-Oblique.ttf"))
    pdfmetrics.registerFont(TTFont("DejaVuMono", _fd / "DejaVuSansMono.ttf"))
    pdfmetrics.registerFontFamily("DejaVu", normal="DejaVu", bold="DejaVu-Bold", italic="DejaVu-Obl", boldItalic="DejaVu-Bold")
    FONT, FONT_B, FONT_I, FONT_M = "DejaVu", "DejaVu-Bold", "DejaVu-Obl", "DejaVuMono"
except Exception as e:
    print("font registration failed, falling back to Helvetica:", e)

PAGE = A4
LMARGIN = RMARGIN = 16 * mm
TMARGIN = BMARGIN = 15 * mm
CONTENT_W = PAGE[0] - LMARGIN - RMARGIN


def _inline(t: str) -> str:
    """Markdown inline -> reportlab mini-HTML."""
    t = t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    t = re.sub(r"`([^`]+)`", rf'<font face="{FONT_M}" size="8">\1</font>', t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", t)
    t = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<link href="\2" color="blue">\1</link>', t)
    return t


def build(md_path: Path, pdf_path: Path) -> None:
    base = md_path.parent
    lines = md_path.read_text(encoding="utf-8").splitlines()

    ss = getSampleStyleSheet()
    styles = {
        "h1": ParagraphStyle("h1", parent=ss["Heading1"], fontName=FONT_B, fontSize=17, spaceBefore=10, spaceAfter=8, textColor=colors.HexColor("#1a3c6e")),
        "h2": ParagraphStyle("h2", parent=ss["Heading2"], fontName=FONT_B, fontSize=13.5, spaceBefore=12, spaceAfter=6, textColor=colors.HexColor("#1a3c6e")),
        "h3": ParagraphStyle("h3", parent=ss["Heading3"], fontName=FONT_B, fontSize=11.5, spaceBefore=9, spaceAfter=4, textColor=colors.HexColor("#33527a")),
        "body": ParagraphStyle("body", parent=ss["BodyText"], fontName=FONT, fontSize=9.5, leading=13, spaceAfter=5),
        "bullet": ParagraphStyle("bullet", parent=ss["BodyText"], fontName=FONT, fontSize=9.5, leading=13, leftIndent=14, bulletIndent=4, spaceAfter=2),
        "cell": ParagraphStyle("cell", parent=ss["BodyText"], fontName=FONT, fontSize=8.2, leading=10),
        "cellh": ParagraphStyle("cellh", parent=ss["BodyText"], fontName=FONT_B, fontSize=8.2, leading=10, textColor=colors.white),
        "cap": ParagraphStyle("cap", parent=ss["BodyText"], fontName=FONT_I, fontSize=8, leading=10, textColor=colors.grey, spaceBefore=2, spaceAfter=10, alignment=1),
    }

    flow = []
    i = 0
    para_buf: list[str] = []

    def flush_para():
        if para_buf:
            flow.append(Paragraph(_inline(" ".join(para_buf)), styles["body"]))
            para_buf.clear()

    def add_image(alt, rel):
        p = (base / rel).resolve()
        if not p.exists():
            flow.append(Paragraph(f"[missing image: {rel}]", styles["body"])); return
        try:
            w, h = PILImage.open(p).size
        except Exception:
            flow.append(Paragraph(f"[unreadable image: {rel}]", styles["body"])); return
        dw = CONTENT_W
        dh = dw * h / w
        max_h = 150 * mm
        if dh > max_h:
            dh = max_h; dw = dh * w / h
        flow.append(Image(str(p), width=dw, height=dh))
        if alt:
            flow.append(Paragraph(_inline(alt), styles["cap"]))

    def add_table(rows):
        # rows: list of list[str]; first row header
        data = []
        for r_i, row in enumerate(rows):
            style = styles["cellh"] if r_i == 0 else styles["cell"]
            data.append([Paragraph(_inline(c), style) for c in row])
        ncol = max(len(r) for r in data)
        for r in data:
            while len(r) < ncol:
                r.append(Paragraph("", styles["cell"]))
        t = Table(data, colWidths=[CONTENT_W / ncol] * ncol, repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#33527a")),
            ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#b0b0b0")),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#eef2f7")]),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ]))
        flow.append(t)
        flow.append(Spacer(1, 6))

    while i < len(lines):
        ln = lines[i].rstrip()
        # image
        m = re.match(r"^!\[(.*?)\]\((.*?)\)\s*$", ln)
        if m:
            flush_para(); add_image(m.group(1), m.group(2)); i += 1; continue
        # table block
        if ln.startswith("|") and i + 1 < len(lines) and re.match(r"^\|[ :\-\|]+\|?\s*$", lines[i + 1]):
            flush_para()
            rows = []
            header = [c.strip() for c in ln.strip("|").split("|")]
            rows.append(header)
            i += 2  # skip separator
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            add_table(rows); continue
        # heading
        m = re.match(r"^(#{1,3})\s+(.*)$", ln)
        if m:
            flush_para()
            lvl = len(m.group(1))
            flow.append(Paragraph(_inline(m.group(2)), styles[f"h{lvl}"]))
            i += 1; continue
        # hr
        if re.match(r"^---+\s*$", ln):
            flush_para(); flow.append(HRFlowable(width="100%", thickness=0.6, color=colors.HexColor("#b0b0b0"), spaceBefore=6, spaceAfter=6)); i += 1; continue
        # bullet
        m = re.match(r"^\s*[-*]\s+(.*)$", ln)
        if m:
            flush_para()
            flow.append(Paragraph(_inline(m.group(1)), styles["bullet"], bulletText="•"))
            i += 1; continue
        # blank
        if not ln.strip():
            flush_para(); i += 1; continue
        # paragraph text
        para_buf.append(ln.strip())
        i += 1
    flush_para()

    doc = SimpleDocTemplate(str(pdf_path), pagesize=PAGE,
                            leftMargin=LMARGIN, rightMargin=RMARGIN,
                            topMargin=TMARGIN, bottomMargin=BMARGIN,
                            title="Rayleigh calibration report")
    doc.build(flow)
    print(f"wrote {pdf_path}")


if __name__ == "__main__":
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/REPORT.md")
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else src.with_suffix(".pdf")
    build(src, out)
