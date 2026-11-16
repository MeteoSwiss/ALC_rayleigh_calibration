"""Build the 20-slide calibration-status presentation with python-pptx.

Theme: "night-sky lidar" — deep navy with section accents (teal = status/Rayleigh,
green = Mini-MPL/532 nm, red = the sign error). Wavelength trio (1064/910/532 nm)
is the recurring motif. All numbers come from the project validation reports.
"""
import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from PIL import Image

# ----------------------------------------------------------------- palette
INK    = "0E1B33"   # darkest navy (dark backgrounds)
NAVY   = "13294B"   # primary dark
SKY    = "2E7DA1"   # teal-blue  (status / Rayleigh section)
GREEN  = "2FBF71"   # 532 nm / Mini-MPL / "best"
GREEND = "1E7A48"   # darker green for text on light
AMBER  = "E9A23B"
RED    = "D64550"
LIGHT  = "F5F7FB"   # light page background
CARD   = "FFFFFF"
CARDLN = "E2E8F0"
MUTED  = "5A6B82"
INKTXT = "1A2A44"
WHITE  = "FFFFFF"
ICE    = "AEBED2"   # light text on dark

HEAD = "Trebuchet MS"
BODY = "Calibri"
MONO = "Consolas"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIG  = os.path.join(HERE, "figs")
SAMP = os.path.join(ROOT, "tests", "sample_plots", "plots")

F_RANK   = os.path.join(FIG, "fig_method_ranking.png")
F_V2     = os.path.join(FIG, "fig_v2_vs_v11.png")
F_MPL    = os.path.join(FIG, "fig_minimpl_precision.png")
F_SIGN   = os.path.join(FIG, "fig_signfix_impact.png")
F_SCHEM  = os.path.join(FIG, "fig_rayleigh_schematic.png")
F_NEG    = os.path.join(FIG, "fig_signerror_negative.png")
P_CHM    = os.path.join(SAMP, "0-20000-0-06610", "2026", "20260225_0-20000-0-06610_rayleigh_diag_compact.png")
P_MPL    = os.path.join(SAMP, "0-20000-0-07014", "2026", "20260423_0-20000-0-07014_rayleigh_diag_compact.png")


def _crop_top(src, frac=0.035):
    """Trim the raw debug super-title strip off a diagnostic PNG."""
    im = Image.open(src)
    w, h = im.size
    out = os.path.join(FIG, "_crop_" + os.path.basename(src))
    im.crop((0, int(h * frac), w, h)).save(out)
    return out


P_CHM = _crop_top(P_CHM)
P_MPL = _crop_top(P_MPL)

prs = Presentation()
prs.slide_width  = Inches(13.333)
prs.slide_height = Inches(7.5)
SW, SH = 13.333, 7.5
BLANK = prs.slide_layouts[6]


# ------------------------------------------------------------- helpers
def hexc(s):
    return RGBColor.from_string(s)


def slide():
    return prs.slides.add_slide(BLANK)


def set_bg(s, color):
    s.background.fill.solid()
    s.background.fill.fore_color.rgb = hexc(color)


def no_line(shp):
    shp.line.fill.background()


def rect(s, x, y, w, h, fill, line=None, line_w=0.75, rounded=False, radius=0.06):
    shape = MSO_SHAPE.ROUNDED_RECTANGLE if rounded else MSO_SHAPE.RECTANGLE
    shp = s.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    shp.fill.solid()
    shp.fill.fore_color.rgb = hexc(fill)
    if line:
        shp.line.color.rgb = hexc(line)
        shp.line.width = Pt(line_w)
    else:
        no_line(shp)
    shp.shadow.inherit = False
    if rounded:
        try:
            shp.adjustments[0] = radius
        except Exception:
            pass
    return shp


def oval(s, x, y, d, fill, line=None, line_w=1.0):
    shp = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(x), Inches(y), Inches(d), Inches(d))
    shp.fill.solid()
    shp.fill.fore_color.rgb = hexc(fill)
    if line:
        shp.line.color.rgb = hexc(line)
        shp.line.width = Pt(line_w)
    else:
        no_line(shp)
    shp.shadow.inherit = False
    return shp


def _set_indent(p, marl_in, indent_in):
    pPr = p._p.get_or_add_pPr()
    pPr.set("marL", str(int(Inches(marl_in))))
    pPr.set("indent", str(int(Inches(indent_in))))


def paras(s, x, y, w, h, items, size=15, color=INKTXT, font=BODY,
          anchor=MSO_ANCHOR.TOP, align=PP_ALIGN.LEFT):
    """items: list of dicts. Each: {runs|text, bullet, bcolor, glyph, level,
       size, color, bold, italic, font, space_before, space_after, line_spacing, align}"""
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    for m in ("margin_left", "margin_right", "margin_top", "margin_bottom"):
        setattr(tf, m, 0)
    for i, it in enumerate(items):
        if isinstance(it, str):
            it = {"text": it}
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = it.get("align", align)
        if it.get("space_before") is not None:
            p.space_before = Pt(it["space_before"])
        p.space_after = Pt(it.get("space_after", 6))
        if it.get("line_spacing"):
            p.line_spacing = it["line_spacing"]
        level = it.get("level", 0)
        runs = it.get("runs")
        if runs is None:
            runs = [(it.get("text", ""), {})]
        bullet = it.get("bullet", False)
        if bullet:
            base = 0.30 * level
            _set_indent(p, base + 0.26, -0.26)
            g = p.add_run()
            g.text = it.get("glyph", "▪") + "  "
            g.font.size = Pt(it.get("size", size))
            g.font.bold = True
            g.font.name = font
            g.font.color.rgb = hexc(it.get("bcolor", SKY))
        elif level:
            _set_indent(p, 0.30 * level, 0)
        for t, opt in runs:
            r = p.add_run()
            r.text = t
            f = r.font
            f.size = Pt(opt.get("size", it.get("size", size)))
            f.bold = opt.get("bold", it.get("bold", False))
            f.italic = opt.get("italic", it.get("italic", False))
            f.name = opt.get("font", it.get("font", font))
            f.color.rgb = hexc(opt.get("color", it.get("color", color)))
    return tb


def text(s, x, y, w, h, content, size=15, color=INKTXT, bold=False, italic=False,
         font=BODY, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP):
    return paras(s, x, y, w, h,
                 [{"text": content, "size": size, "color": color, "bold": bold,
                   "italic": italic, "font": font, "align": align, "space_after": 0}],
                 anchor=anchor)


def img_fit(s, path, x, y, w, h):
    iw, ih = Image.open(path).size
    ar = iw / ih
    box_ar = w / h
    if ar > box_ar:
        dw, dh = w, w / ar
    else:
        dh, dw = h, h * ar
    px = x + (w - dw) / 2
    py = y + (h - dh) / 2
    return s.shapes.add_picture(path, Inches(px), Inches(py), Inches(dw), Inches(dh))


def img_frame(s, path, x, y, w, h, pad=0.0):
    """white rounded frame behind a contained image"""
    rect(s, x - pad, y - pad, w + 2 * pad, h + 2 * pad, CARD, line=CARDLN, rounded=True, radius=0.03)
    img_fit(s, path, x, y, w, h)


SECTION = {"status": SKY, "rayleigh": SKY, "minimpl": GREEN, "sign": RED, "close": GREEN}


def accent_bar(s, color, w=0.14):
    rect(s, 0, 0, w, SH, color)


def wl_dots(s, x, y, dark=True, d=0.16, gap=1.55, size=10.5):
    trio = [("1064 nm", "7FA8D0" if dark else NAVY),
            ("910 nm", SKY),
            ("532 nm", GREEN)]
    cx = x
    for lbl, col in trio:
        oval(s, cx, y, d, col)
        text(s, cx + d + 0.07, y - 0.05, 1.3, 0.3, lbl, size=size,
             color=ICE if dark else MUTED, bold=True, font=BODY, anchor=MSO_ANCHOR.MIDDLE)
        cx += gap


def footer(s, idx, dark=False):
    c = "55698A" if dark else "9AA9BE"
    text(s, 0.62, 7.08, 9, 0.3, "ALC calibration status  ·  E-PROFILE  ·  M. Hervo, MeteoSwiss  ·  June 2026",
         size=8.5, color=c)
    text(s, 12.1, 7.06, 0.7, 0.3, f"{idx:02d}", size=10, color=c, bold=True, align=PP_ALIGN.RIGHT)


def header(s, kicker, title, section, idx, title_size=29):
    set_bg(s, LIGHT)
    accent = SECTION[section]
    accent_bar(s, accent)
    text(s, 0.62, 0.40, 11.8, 0.3, kicker.upper(), size=12.5, color=accent, bold=True, font=HEAD)
    text(s, 0.60, 0.70, 12.2, 1.0, title, size=title_size, color=INK, bold=True, font=HEAD)
    footer(s, idx)
    return accent


def divider(s, idx, kicker, title, section, subtitle, big_no):
    set_bg(s, INK)
    accent = SECTION[section]
    accent_bar(s, accent, w=0.24)
    text(s, 0.9, 1.95, 3, 2.2, big_no, size=150, color=NAVY, bold=True, font=HEAD)
    text(s, 0.95, 2.78, 11, 0.4, kicker.upper(), size=15, color=accent, bold=True, font=HEAD)
    text(s, 0.90, 3.18, 11.6, 1.5, title, size=40, color=WHITE, bold=True, font=HEAD)
    paras(s, 0.95, 4.62, 10.8, 1.2,
          [{"text": subtitle, "size": 16, "color": ICE, "line_spacing": 1.15}])
    wl_dots(s, 0.95, 6.35, dark=True)
    footer(s, idx, dark=True)


def chip(s, x, y, w, h, label, fill, txtcolor=WHITE, size=12.5):
    rect(s, x, y, w, h, fill, rounded=True, radius=0.22)
    text(s, x, y, w, h, label, size=size, color=txtcolor, bold=True, font=HEAD,
         align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)


def stat_card(s, x, y, w, h, number, label, color, num_size=33):
    rect(s, x, y, w, h, CARD, line=CARDLN, rounded=True, radius=0.07)
    rect(s, x, y, 0.10, h, color)  # thin left accent (sharp, aligns with card edge)
    text(s, x + 0.22, y + 0.14, w - 0.34, 0.7, number, size=num_size, color=color, bold=True, font=HEAD)
    paras(s, x + 0.24, y + 0.86, w - 0.42, h - 0.9,
          [{"text": label, "size": 11.5, "color": MUTED, "line_spacing": 1.0}])


# ===================================================================== 1 TITLE
s = slide()
set_bg(s, INK)
# subtle motif: three vertical "beam" bars on the right, section colors
for bx, col, htop in [(11.0, "1E3A5F", 0.0), (11.55, SKY, 0.0), (12.1, GREEN, 0.0)]:
    rect(s, bx, 0, 0.14, SH, col if col != SKY else "1C4E6B")
rect(s, 11.0, 0, 1.24, SH, INK)  # mask -> remove (kept simple); overwrite beams subtly
# (re-add slim beams over the mask)
rect(s, 11.62, 1.4, 0.06, 4.7, "26527A")
rect(s, 11.95, 0.9, 0.06, 5.7, "2F6E5B")
rect(s, 12.28, 2.0, 0.06, 3.6, "2E7DA1")
accent_bar(s, GREEN, w=0.24)
text(s, 0.95, 1.55, 9, 0.4, "METEOSWISS  ·  E-PROFILE  ·  AUTOMATIC LIDARS & CEILOMETERS",
     size=13, color=SKY, bold=True, font=HEAD)
paras(s, 0.92, 2.05, 10.4, 2.4,
      [{"runs": [("ALC Calibration", {})], "size": 50, "color": WHITE, "bold": True, "font": HEAD,
        "space_after": 2, "line_spacing": 1.0},
       {"runs": [("Status & Rayleigh developments", {})], "size": 50, "color": GREEN, "bold": True,
        "font": HEAD, "line_spacing": 1.0}])
paras(s, 0.95, 4.35, 10.2, 1.4,
      [{"text": "Molecular calibration across the network · a focus on Mini-MPL (532 nm) ·",
        "size": 17, "color": ICE, "line_spacing": 1.2, "space_after": 2},
       {"runs": [("and the Klett ", {"color": ICE}),
                 ("sign-error correction", {"color": AMBER, "bold": True}),
                 ("  (E-PROF v1.0 → v1.1)", {"color": ICE})],
        "size": 17, "line_spacing": 1.2}])
wl_dots(s, 0.97, 5.95, dark=True)
text(s, 0.95, 6.85, 11, 0.3, "M. Hervo, MeteoSwiss  ·  June 2026  ·  eprofile-calibration v2.0.0",
     size=12, color="6E83A6", bold=True)

# ===================================================================== 2 AGENDA
s = slide()
header(s, "Overview", "What this deck covers", "status", 2)
cards = [
    ("1", "Calibration status", "Where the E-PROFILE ALC calibration stands in June 2026 — "
     "instruments, methods and validated results.", SKY),
    ("2", "Rayleigh & Mini-MPL", "Developments in molecular calibration and the methods that drive "
     "it — with a special focus on the 532 nm Mini-MPL.", GREEN),
    ("3", "The Klett sign error", "What the E-PROF v1.0 sign error was, and how the results changed "
     "once it was fixed in v1.1.", RED),
]
cw, gap = 3.78, 0.39
x0 = 0.62
y0 = 2.15
for i, (no, title, body, col) in enumerate(cards):
    x = x0 + i * (cw + gap)
    rect(s, x, y0, cw, 3.9, CARD, line=CARDLN, rounded=True, radius=0.05)
    rect(s, x, y0, cw, 0.16, col)
    oval(s, x + 0.35, y0 + 0.5, 0.95, col)
    text(s, x + 0.35, y0 + 0.5, 0.95, 0.95, no, size=34, color=WHITE, bold=True, font=HEAD,
         align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
    text(s, x + 0.35, y0 + 1.75, cw - 0.7, 0.6, title, size=19, color=INK, bold=True, font=HEAD)
    paras(s, x + 0.37, y0 + 2.35, cw - 0.72, 1.4,
          [{"text": body, "size": 13, "color": MUTED, "line_spacing": 1.18}])
text(s, 0.62, 6.35, 12, 0.4,
     "Three wavelengths, one calibration framework:  CHM15k 1064 nm  ·  CL61 / CL51 / CL31 910 nm  ·  Mini-MPL 532 nm",
     size=12.5, color=MUTED, italic=True)

# ===================================================================== 3 BIG PICTURE
s = slide()
header(s, "Calibration status", "From raw signal to calibrated backscatter", "status", 3)
flow = [("Raw signal", "range-corrected\nRCS / r²", NAVY),
        ("Lidar constant Cₗ", "Rayleigh  or\nliquid-cloud", SKY),
        ("Calibrated β_att", "Kalman-smoothed\ndaily series", GREEN),
        ("E-PROFILE L2", "operational\nL1→L2 hub", AMBER)]
fw, fgap = 2.72, 0.55
fx = 0.7
fy = 2.05
for i, (t, sub, col) in enumerate(flow):
    x = fx + i * (fw + fgap)
    rect(s, x, fy, fw, 1.5, CARD, line=CARDLN, rounded=True, radius=0.08)
    rect(s, x, fy, fw, 0.13, col)
    text(s, x + 0.18, fy + 0.28, fw - 0.36, 0.5, t, size=15.5, color=INK, bold=True, font=HEAD)
    paras(s, x + 0.2, fy + 0.74, fw - 0.4, 0.7,
          [{"text": sub.replace("\n", "  "), "size": 11.5, "color": MUTED, "line_spacing": 1.0}])
    if i < 3:
        text(s, x + fw + 0.04, fy + 0.35, 0.5, 0.8, "›", size=34, color=col, bold=True,
             font=HEAD, align=PP_ALIGN.CENTER)
# convention + hub defaults
rect(s, 0.7, 4.05, 6.05, 2.55, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 0.95, 4.25, 5.6, 0.4, "What the calibration delivers", size=15, color=SKY, bold=True, font=HEAD)
paras(s, 0.97, 4.72, 5.6, 1.9,
      [{"text": "Absolute attenuated backscatter β_att for every ALC in the network — the L2 product.",
        "bullet": True, "bcolor": SKY, "size": 13, "color": INKTXT, "space_after": 7},
       {"runs": [("Wiegner lidar constant convention:  ", {"color": INKTXT}),
                 ("Cₗ = RCS / β_att", {"color": INK, "bold": True, "font": MONO})],
        "bullet": True, "bcolor": SKY, "size": 13, "space_after": 7},
       {"text": "Per-night constants are fed to a Kalman filter for a smooth, drift-aware series.",
        "bullet": True, "bcolor": SKY, "size": 13, "color": INKTXT}])
rect(s, 7.0, 4.05, 5.6, 2.55, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 7.25, 4.25, 5.2, 0.4, "Operational L1→L2 hub — default Cₗ", size=15, color=AMBER, bold=True, font=HEAD)
paras(s, 7.27, 4.72, 5.2, 1.9,
      [{"runs": [("CHM15k 3·10¹¹   ·   CL61 1   ·   Mini-MPL 5·10⁵", {"color": INKTXT, "font": MONO, "size": 12.5})],
        "bullet": True, "bcolor": AMBER, "space_after": 7},
       {"runs": [("CL31 / CL51 1·10⁸ ", {"color": INKTXT, "font": MONO, "size": 12.5}),
                 ("— still uncalibrated", {"color": RED, "bold": True, "size": 12.5})],
        "bullet": True, "bcolor": AMBER, "space_after": 7},
       {"text": "Closing that gap (cloud calibration of CL31/CL51) is part of this work.",
        "bullet": True, "bcolor": AMBER, "size": 12.5, "color": INKTXT}])

# ===================================================================== 4 INSTRUMENTS & METHODS
s = slide()
header(s, "Calibration status", "Instruments & calibration methods", "status", 4)
# manual table on the left
tx, ty, tw = 0.7, 2.1, 7.5
cols = [tx + 0.15, tx + 2.55, tx + 3.7, tx + 5.0, tx + 6.15]
rect(s, tx, ty, tw, 0.55, NAVY)
hdrs = ["Instrument", "λ", "Rayleigh", "Cloud", "Notes"]
hx = [cols[0], cols[1], cols[2], cols[3], cols[4]]
for hxi, hd in zip(hx, hdrs):
    text(s, hxi, ty, 2.3, 0.55, hd, size=12.5, color=WHITE, bold=True, font=HEAD, anchor=MSO_ANCHOR.MIDDLE)
rows = [
    ("CHM15k", "1064 nm", "Y", "-", "molecular workhorse"),
    ("CL61", "910 nm", "Y", "Y", "new: full network"),
    ("CL51", "910 nm", "-", "Y", "cloud only"),
    ("CL31", "910 nm", "-", "Y", "cloud only"),
    ("Mini-MPL", "532 nm", "Y", "-", "Rayleigh only"),
]
rh = 0.74
for i, (nm, lam, ray, cld, note) in enumerate(rows):
    ry = ty + 0.55 + i * rh
    rect(s, tx, ry, tw, rh, CARD if i % 2 == 0 else "EEF2F8")
    hi = GREEN if nm in ("CHM15k", "CL61", "Mini-MPL") else MUTED
    text(s, cols[0], ry, 2.3, rh, nm, size=13.5, color=INK, bold=True, font=HEAD, anchor=MSO_ANCHOR.MIDDLE)
    text(s, cols[1], ry, 1.1, rh, lam, size=12.5, color=INKTXT, anchor=MSO_ANCHOR.MIDDLE)
    text(s, cols[2], ry, 1.2, rh, "●" if ray == "Y" else "—", size=14,
         color=GREEN if ray == "Y" else "C2CCDA", bold=True, anchor=MSO_ANCHOR.MIDDLE)
    text(s, cols[3], ry, 1.2, rh, "●" if cld == "Y" else "—", size=14,
         color=SKY if cld == "Y" else "C2CCDA", bold=True, anchor=MSO_ANCHOR.MIDDLE)
    text(s, cols[4], ry, 1.5, rh, note, size=11, color=MUTED, anchor=MSO_ANCHOR.MIDDLE)
# clean outline around the whole table (no fill, so rows/header stay visible)
outline = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(tx), Inches(ty), Inches(tw), Inches(0.55 + len(rows) * rh))
outline.fill.background(); outline.line.color.rgb = hexc(CARDLN); outline.line.width = Pt(1.0); outline.shadow.inherit = False

# method family cards on the right
mx = 8.55
rect(s, mx, 2.1, 4.05, 2.05, CARD, line=CARDLN, rounded=True, radius=0.06)
rect(s, mx, 2.1, 0.12, 2.05, SKY)
text(s, mx + 0.28, 2.28, 3.6, 0.4, "Rayleigh (molecular)", size=15, color=SKY, bold=True, font=HEAD)
paras(s, mx + 0.3, 2.78, 3.55, 1.3,
      [{"text": "Molecular scattering in clean air as a natural reference.", "size": 12, "color": INKTXT,
        "line_spacing": 1.12, "space_after": 5},
       {"text": "CHM15k · CL61 · Mini-MPL · (CHM8k)", "size": 11.5, "color": MUTED, "italic": True}])
rect(s, mx, 4.35, 4.05, 2.05, CARD, line=CARDLN, rounded=True, radius=0.06)
rect(s, mx, 4.35, 0.12, 2.05, GREEN)
text(s, mx + 0.28, 4.53, 3.6, 0.4, "Liquid-cloud (O'Connor)", size=15, color=GREEND, bold=True, font=HEAD)
paras(s, mx + 0.3, 5.03, 3.55, 1.3,
      [{"runs": [("Opaque liquid clouds, theoretical S = 18.8 sr. ", {"color": INKTXT})],
        "size": 12, "line_spacing": 1.12, "space_after": 5},
       {"text": "Python port bit-for-bit vs MATLAB (≤1.4·10⁻⁵).", "size": 11.5, "color": MUTED, "italic": True}])

# ===================================================================== 5 STATUS SNAPSHOT
s = slide()
header(s, "Calibration status", "Where the calibration stands — June 2026", "status", 5)
stats = [
    ("24", "instruments validated on L1 2026\n(10 CHM15k · 4 Mini-MPL · 10 CL61)", SKY),
    ("10", "CL61 newly validated network-wide\nmedian σ_SD 7.7% — most precise group", GREEN),
    ("6", "selectable molecular-window methods\n(E-PROF v0.25 → v2)", NAVY),
    ("v2.0", "eprofile-calibration package\nRayleigh + cloud + WV correction", AMBER),
]
sw = 2.95
for i, (no, lab, col) in enumerate(stats):
    stat_card(s, 0.62 + i * (sw + 0.16), 2.1, sw, 1.55, no, lab.replace("\n", " "), col)
# two columns
rect(s, 0.62, 3.95, 6.0, 2.7, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 0.88, 4.13, 5.5, 0.4, "✓  Working / validated", size=15, color=GREEND, bold=True, font=HEAD)
paras(s, 0.9, 4.62, 5.55, 2.0,
      [{"text": "Rayleigh calibration — 6 methods, full network", "bullet": True, "bcolor": GREEN, "size": 12.5, "space_after": 6},
       {"text": "Water-vapour correction (Python ↔ MATLAB parity)", "bullet": True, "bcolor": GREEN, "size": 12.5, "space_after": 6},
       {"text": "Liquid-cloud calibration (bit-for-bit vs MATLAB)", "bullet": True, "bcolor": GREEN, "size": 12.5, "space_after": 6},
       {"text": "CL61 validated + independent cloud cross-check", "bullet": True, "bcolor": GREEN, "size": 12.5, "space_after": 6},
       {"text": "Static monitoring dashboard", "bullet": True, "bcolor": GREEN, "size": 12.5}])
rect(s, 6.75, 3.95, 5.85, 2.7, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 7.0, 4.13, 5.4, 0.4, "⟳  In progress", size=15, color=AMBER, bold=True, font=HEAD)
paras(s, 7.02, 4.62, 5.45, 2.0,
      [{"text": "Upgrade operations from v1.0 → E-PROF v2 (sign-corrected)", "bullet": True, "bcolor": AMBER, "size": 12.5, "space_after": 6},
       {"text": "Calibrate CL31 / CL51 via cloud method (fill the gap)", "bullet": True, "bcolor": AMBER, "size": 12.5, "space_after": 6},
       {"text": "Cloudnet CL61 pipeline (ACTRIS cross-network)", "bullet": True, "bcolor": AMBER, "size": 12.5, "space_after": 6},
       {"text": "Paper figures (E-PROFILE ALC paper)", "bullet": True, "bcolor": AMBER, "size": 12.5}])

# ===================================================================== 6 DIVIDER RAYLEIGH
s = slide()
divider(s, 6, "Part 1", "Rayleigh calibration\ndevelopments", "rayleigh",
        "How molecular calibration works, the six method versions, and which one wins on precision.",
        "01")

# ===================================================================== 7 HOW RAYLEIGH WORKS
s = slide()
header(s, "Rayleigh developments", "Rayleigh (molecular) calibration in a nutshell", "rayleigh", 7)
paras(s, 0.7, 2.05, 5.5, 4.6,
      [{"text": "Select clear, aerosol-free nights (darkness-adaptive, cloud/fog screened).",
        "bullet": True, "bcolor": SKY, "size": 14, "space_after": 11},
       {"runs": [("Build a molecular reference profile β_mol from a model atmosphere ",
                  {"color": INKTXT}),
                 ("(US 1976 or CAMS, Bucholtz 1995).", {"color": MUTED})],
        "bullet": True, "bcolor": SKY, "size": 14, "space_after": 11},
       {"text": "Find a clean molecular window (2–6 km) where signal tracks β_mol.",
        "bullet": True, "bcolor": SKY, "size": 14, "space_after": 11},
       {"runs": [("Fit  signal ∝ β_mol  → slope gives the lidar constant ", {"color": INKTXT}),
                 ("Cₗ", {"color": INK, "bold": True})],
        "bullet": True, "bcolor": SKY, "size": 14, "space_after": 11},
       {"text": "Klett inversion below the window recovers the aerosol profile.",
        "bullet": True, "bcolor": SKY, "size": 14, "space_after": 11},
       {"runs": [("Sensitivity test ", {"color": INKTXT, "bold": True}),
                 ("(LR ±20 sr, altitude ±200 m, time subsets) → median Cₗ + 2σ.",
                  {"color": INKTXT})],
        "bullet": True, "bcolor": SKY, "size": 14}])
img_frame(s, F_SCHEM, 6.5, 2.0, 6.2, 3.6, pad=0.12)
text(s, 6.5, 5.78, 6.2, 0.5,
     "Molecular window (green) → linear fit → Cₗ; Klett + sensitivity give the uncertainty.",
     size=11, color=MUTED, italic=True, align=PP_ALIGN.CENTER)

# ===================================================================== 8 CHM15k EXAMPLE
s = slide()
header(s, "Rayleigh developments", "Worked example — CHM15k (1064 nm), Payerne", "rayleigh", 8)
img_frame(s, P_CHM, 0.7, 2.0, 8.4, 4.75, pad=0.1)
rect(s, 9.35, 2.0, 3.3, 4.75, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 9.58, 2.2, 2.9, 0.4, "What you see", size=15, color=SKY, bold=True, font=HEAD)
paras(s, 9.6, 2.72, 2.85, 4.0,
      [{"text": "Night-mean signal vs the molecular reference.", "bullet": True, "bcolor": SKY, "size": 12, "space_after": 8},
       {"text": "Detected molecular window highlighted.", "bullet": True, "bcolor": SKY, "size": 12, "space_after": 8},
       {"text": "Linear fit → Cₗ for the night.", "bullet": True, "bcolor": SKY, "size": 12, "space_after": 8},
       {"text": "Retrieved aerosol backscatter below.", "bullet": True, "bcolor": SKY, "size": 12, "space_after": 8},
       {"text": "Quality flag + 2σ uncertainty.", "bullet": True, "bcolor": SKY, "size": 12, "space_after": 8},
       {"text": "2026-02-25 · sample data shipped with the package.", "size": 11, "color": MUTED, "italic": True, "space_before": 6}])

# ===================================================================== 9 SIX METHODS
s = slide()
header(s, "Rayleigh developments", "Six molecular-window methods, by E-PROF version", "rayleigh", 9)
methods = [
    ("v0.25", "MATLAB Auto_Calib_25 — forced intercept, min ΣRMSE", NAVY, ""),
    ("v1.0", "Legacy window + the Klett sign error — current operational version", RED, "operational"),
    ("v1.1", "Legacy window, sign-corrected — but admits aerosol", AMBER, "sign fixed"),
    ("v1.2", "Above-aerosol + R²/shape/slope gates, max R²", SKY, ""),
    ("v2", "All gates + temporal aerosol rejection + cell flagging", GREEN, "recommended"),
    ("EARLINET", "Shape + SNR gate, lowest qualifying window", "6D7E96", ""),
    ("Bellini", "ALICENET 3–7 km, Breusch–Godfrey residual test", "9AA7B8", ""),
]
my = 2.1
rhh = 0.62
for i, (ver, desc, col, tag) in enumerate(methods):
    y = my + i * (rhh + 0.04)
    rect(s, 0.7, y, 11.9, rhh, CARD, line=CARDLN, rounded=True, radius=0.12)
    disp = ver if ver in ("EARLINET", "Bellini") else "E-PROF " + ver
    chip(s, 0.85, y + 0.1, 1.8, rhh - 0.2, disp, col, size=12)
    text(s, 2.88, y, 7.4, rhh, desc, size=13, color=INKTXT, anchor=MSO_ANCHOR.MIDDLE)
    if tag:
        tcol = {"operational": RED, "sign error": RED, "sign fixed": AMBER,
                "current default": SKY, "recommended": GREEND}.get(tag, MUTED)
        text(s, 10.2, y, 2.3, rhh, "◀ " + tag, size=11.5, color=tcol, bold=True, italic=True,
             anchor=MSO_ANCHOR.MIDDLE)

# ===================================================================== 10 PRECISION RANKING
s = slide()
header(s, "Rayleigh developments", "Which method is most precise?", "rayleigh", 10)
img_frame(s, F_RANK, 0.7, 1.95, 8.5, 4.45, pad=0.1)
rect(s, 9.45, 1.95, 3.2, 4.45, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 9.68, 2.13, 2.8, 0.4, "Takeaways", size=15, color=SKY, bold=True, font=HEAD)
paras(s, 9.7, 2.64, 2.8, 3.7,
      [{"runs": [("E-PROF v2 wins", {"bold": True, "color": GREEND}), ("  (σ_SD 9.5%).", {"color": INKTXT})],
        "bullet": True, "bcolor": GREEN, "size": 12.5, "space_after": 9},
       {"text": "EARLINET a close second at higher yield.", "bullet": True, "bcolor": SKY, "size": 12.5, "space_after": 9},
       {"runs": [("v1.1 / v0.25 look stable but admit aerosol ", {"color": INKTXT}),
                 ("(high within-night spread).", {"color": MUTED})],
        "bullet": True, "bcolor": AMBER, "size": 12.5, "space_after": 9},
       {"runs": [("Recommendation: ", {"bold": True, "color": INK}),
                 ("move operations from v1.0 → v2 (sign-corrected).", {"color": INKTXT})],
        "bullet": True, "bcolor": GREEN, "size": 12.5}])

# ===================================================================== 11 V2 vs V11
s = slide()
header(s, "Rayleigh developments", "E-PROF v2 vs v1.1 across the network", "rayleigh", 11)
img_frame(s, F_V2, 0.7, 1.95, 8.5, 4.5, pad=0.1)
rect(s, 9.45, 1.95, 3.2, 4.5, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 9.68, 2.13, 2.8, 0.4, "Takeaways", size=15, color=SKY, bold=True, font=HEAD)
paras(s, 9.7, 2.64, 2.8, 3.7,
      [{"text": "v2 lifts CHM15k yield by ~+10 pp.", "bullet": True, "bcolor": GREEN, "size": 12.5, "space_after": 9},
       {"runs": [("CL61 σ_SD 9.7 → 7.0 %", {"bold": True, "color": GREEND}), (".", {})],
        "bullet": True, "bcolor": GREEN, "size": 12.5, "space_after": 9},
       {"text": "Mini-MPL already precise (6.4–6.8%).", "bullet": True, "bcolor": SKY, "size": 12.5, "space_after": 9},
       {"runs": [("L1 ≡ L2", {"bold": True, "color": INK}),
                 (" once native L1 is binned to the L2 grid.", {"color": INKTXT})],
        "bullet": True, "bcolor": SKY, "size": 12.5}])

# ===================================================================== 12 DIVIDER MINI-MPL
s = slide()
divider(s, 12, "Part 2 · Special focus", "Mini-MPL\nat 532 nm", "minimpl",
        "A green-wavelength lidar in a 910/1064 nm network — why it is special, and how well it calibrates.",
        "02")

# ===================================================================== 13 MINI-MPL WHY SPECIAL
s = slide()
header(s, "Focus · Mini-MPL", "A 532 nm system in a 910 / 1064 nm network", "minimpl", 13)
# three stat callouts
m_stats = [("532 nm", "green Nd:YAG\nwavelength", GREEN),
           ("≈16×", "stronger molecular\nsignal vs 1064 nm (λ⁻⁴)", SKY),
           ("4", "network units: Mobile,\nAléria, Lille, Brest", NAVY)]
for i, (no, lab, col) in enumerate(m_stats):
    stat_card(s, 0.62 + i * 3.05, 2.1, 2.9, 1.5, no, lab.replace("\n", " "), col, num_size=30)
# feature rows
feats = [
    ("Easy clean-air detection", "Strong λ⁻⁴ molecular return makes the molecular window easy to find.", GREEN),
    ("No water-vapour correction", "532 nm sits outside the 910 nm WV absorption band — unlike the ceilometers.", SKY),
    ("More aerosol-sensitive", "Aerosol scattering is much stronger at 532 nm — which matters for the sign error.", AMBER),
    ("Rayleigh only", "No liquid-cloud method here — molecular calibration is the route for Mini-MPL.", NAVY),
]
fy = 3.95
for i, (t, d, col) in enumerate(feats):
    y = fy + i * 0.68
    oval(s, 0.72, y + 0.03, 0.32, col)
    text(s, 1.2, y - 0.06, 4.0, 0.5, t, size=13.5, color=INK, bold=True, font=HEAD, anchor=MSO_ANCHOR.MIDDLE)
    text(s, 5.2, y - 0.06, 7.4, 0.6, d, size=12, color=MUTED, anchor=MSO_ANCHOR.MIDDLE)

# ===================================================================== 14 MINI-MPL EXAMPLE
s = slide()
header(s, "Focus · Mini-MPL", "Worked example — Mini-MPL (532 nm), Lille", "minimpl", 14)
img_frame(s, P_MPL, 0.7, 2.0, 8.4, 4.75, pad=0.1)
rect(s, 9.35, 2.0, 3.3, 4.75, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 9.58, 2.2, 2.9, 0.4, "Same method, 532 nm", size=14.5, color=GREEND, bold=True, font=HEAD)
paras(s, 9.6, 2.72, 2.85, 4.0,
      [{"text": "Identical molecular-window pipeline as the ceilometers.", "bullet": True, "bcolor": GREEN, "size": 12, "space_after": 8},
       {"text": "Molecular signal is strong and clean at 532 nm.", "bullet": True, "bcolor": GREEN, "size": 12, "space_after": 8},
       {"text": "No WV correction needed at this wavelength.", "bullet": True, "bcolor": GREEN, "size": 12, "space_after": 8},
       {"text": "Lille (FR), night of 2026-04-23 — sample data.", "size": 11, "color": MUTED, "italic": True, "space_before": 6}])

# ===================================================================== 15 MINI-MPL RESULTS
s = slide()
header(s, "Focus · Mini-MPL", "Mini-MPL calibration results", "minimpl", 15)
img_frame(s, F_MPL, 0.7, 1.95, 8.5, 4.5, pad=0.1)
rect(s, 9.45, 1.95, 3.2, 4.5, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 9.68, 2.13, 2.8, 0.4, "Takeaways", size=15, color=GREEND, bold=True, font=HEAD)
paras(s, 9.7, 2.64, 2.8, 3.7,
      [{"runs": [("Median σ_SD 8.5%", {"bold": True, "color": GREEND}), (" — CHM15k-class.", {"color": INKTXT})],
        "bullet": True, "bcolor": GREEN, "size": 12.5, "space_after": 9},
       {"text": "v2 is clearly best for Mini-MPL (others 16–18%).", "bullet": True, "bcolor": GREEN, "size": 12.5, "space_after": 9},
       {"runs": [("Cross-wavelength tracker at SIRTA: ", {"color": INKTXT}),
                 ("r = 0.79", {"bold": True, "color": INK}),
                 (" vs 1064 nm.", {"color": INKTXT})],
        "bullet": True, "bcolor": SKY, "size": 12.5, "space_after": 9},
       {"runs": [("Brest", {"bold": True, "color": RED}),
                 (" is the one problem unit (background).", {"color": INKTXT})],
        "bullet": True, "bcolor": RED, "size": 12.5}])

# ===================================================================== 16 DIVIDER SIGN ERROR
s = slide()
divider(s, 16, "Part 3", "The Klett\nsign error", "sign",
        "What the E-PROF v1.0 sign error was — and how the results changed when v1.1 fixed it.",
        "03")

# ===================================================================== 17 WHAT WAS THE ERROR
s = slide()
header(s, "Klett sign error", "What was the sign error?", "sign", 17)
paras(s, 0.7, 2.05, 5.7, 4.7,
      [{"runs": [("A sign error in the ", {"color": INKTXT}),
                 ("Klett–Fernald inversion", {"bold": True, "color": INK}),
                 (" — the classic retrieval error.", {"color": INKTXT})],
        "bullet": True, "bcolor": RED, "size": 14, "space_after": 11},
       {"runs": [("The transmission term used ", {"color": INKTXT}),
                 ("exp(−2·qₜ)", {"font": MONO, "bold": True, "color": RED}),
                 (" instead of ", {"color": INKTXT}),
                 ("exp(+2·qₜ)", {"font": MONO, "bold": True, "color": GREEND}),
                 (".", {"color": INKTXT})],
        "bullet": True, "bcolor": RED, "size": 14, "space_after": 11},
       {"text": "It drives the retrieved aerosol backscatter negative — physically impossible.",
        "bullet": True, "bcolor": RED, "size": 14, "space_after": 11},
       {"runs": [("Documented by ", {"color": INKTXT}),
                 ("Speidel & Vogelmann (2023)", {"bold": True, "color": INK}),
                 (", Applied Optics 62(4), 861.", {"color": INKTXT})],
        "bullet": True, "bcolor": RED, "size": 14, "space_after": 11},
       {"runs": [("Fixed in commit a4e7140; ", {"color": INKTXT, "font": MONO}),
                 ("now corrected in 3 active codebases (Python + 2 MATLAB).",
                  {"color": INKTXT, "font": BODY})],
        "bullet": True, "bcolor": RED, "size": 13.5}])
# code cards
text(s, 6.7, 2.0, 6, 0.35, "In the Klett inversion (atmosphere.py)", size=12.5, color=MUTED, bold=True, font=HEAD)
rect(s, 6.7, 2.45, 5.95, 1.5, "2A1418", line=RED, line_w=1.0, rounded=True, radius=0.05)
text(s, 6.95, 2.62, 5.6, 0.4, "E-PROF v1.0  —  sign error", size=12.5, color=RED, bold=True, font=HEAD)
paras(s, 6.95, 3.12, 5.55, 0.9,
      [{"runs": [("numerator = beta_att * np.exp(", {"color": "E7C9CC", "font": MONO, "size": 12.5}),
                 ("-2", {"color": RED, "font": MONO, "bold": True, "size": 13}),
                 (" * qt)", {"color": "E7C9CC", "font": MONO, "size": 12.5})], "space_after": 2},
       {"runs": [("# aerosol backscatter → negative", {"color": "B07A80", "font": MONO, "size": 10.5})]}])
rect(s, 6.7, 4.2, 5.95, 1.5, "0F2A1C", line=GREEN, line_w=1.0, rounded=True, radius=0.05)
text(s, 6.95, 4.37, 5.6, 0.4, "E-PROF v1.1+  —  corrected", size=12.5, color=GREEN, bold=True, font=HEAD)
paras(s, 6.95, 4.87, 5.55, 0.9,
      [{"runs": [("numerator = beta_att * np.exp(", {"color": "C7E7D5", "font": MONO, "size": 12.5}),
                 ("+2", {"color": GREEN, "font": MONO, "bold": True, "size": 13}),
                 (" * qt)", {"color": "C7E7D5", "font": MONO, "size": 12.5})], "space_after": 2},
       {"runs": [("# physical, positive profile", {"color": "7FB295", "font": MONO, "size": 10.5})]}])
text(s, 6.7, 5.95, 5.95, 0.6,
     "(+ a truncated denominator integral and a total- vs aerosol-OD slope check were corrected too.)",
     size=10.5, color=MUTED, italic=True)

# ===================================================================== 18 GOES NEGATIVE
s = slide()
header(s, "Klett sign error", "The symptom: the retrieval goes negative", "sign", 18)
img_frame(s, F_NEG, 0.7, 2.0, 8.5, 4.6, pad=0.1)
rect(s, 9.45, 2.0, 3.2, 4.6, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 9.68, 2.2, 2.8, 0.4, "Why it matters", size=15, color=RED, bold=True, font=HEAD)
paras(s, 9.7, 2.72, 2.8, 3.8,
      [{"runs": [("The −2 sign", {"bold": True, "color": RED}),
                 (" pushes aerosol backscatter below zero aloft.", {"color": INKTXT})],
        "bullet": True, "bcolor": RED, "size": 12.5, "space_after": 10},
       {"runs": [("It corrupts the downstream ", {"color": INKTXT}),
                 ("β_att / AOD product", {"bold": True, "color": INK}),
                 (".", {"color": INKTXT})],
        "bullet": True, "bcolor": RED, "size": 12.5, "space_after": 10},
       {"runs": [("The +2 fix", {"bold": True, "color": GREEND}),
                 (" restores a physical, positive profile.", {"color": INKTXT})],
        "bullet": True, "bcolor": GREEN, "size": 12.5}])

# ===================================================================== 19 DIFFERENT RESULTS
s = slide()
header(s, "Klett sign error", "v1.0 → v1.1: what actually changed", "sign", 19)
img_frame(s, F_SIGN, 0.7, 1.95, 8.5, 4.55, pad=0.1)
rect(s, 9.45, 1.95, 3.2, 4.55, CARD, line=CARDLN, rounded=True, radius=0.05)
text(s, 9.68, 2.13, 2.85, 0.4, "The key result", size=15, color=RED, bold=True, font=HEAD)
paras(s, 9.7, 2.64, 2.85, 3.8,
      [{"runs": [("Cₗ stability barely moves ", {"color": INKTXT}),
                 ("(σ_SD 10.3 ≈ 11.6%).", {"color": MUTED})],
        "bullet": True, "bcolor": SKY, "size": 12.5, "space_after": 9},
       {"text": "1064 / 910 nm absolute Cₗ shift ~0.5% — negligible.",
        "bullet": True, "bcolor": SKY, "size": 12.5, "space_after": 9},
       {"runs": [("532 nm Mini-MPL: ", {"color": INKTXT}),
                 ("+8 to +18%", {"bold": True, "color": GREEND}),
                 (" — the fingerprint of the fix.", {"color": INKTXT})],
        "bullet": True, "bcolor": GREEN, "size": 12.5, "space_after": 9},
       {"text": "Effect scales with aerosol optical depth → largest at 532 nm.",
        "bullet": True, "bcolor": AMBER, "size": 12}])

# ===================================================================== 20 CONCLUSIONS
s = slide()
set_bg(s, INK)
accent_bar(s, GREEN, w=0.24)
text(s, 0.62, 0.45, 11.8, 0.3, "CONCLUSIONS", size=13, color=GREEN, bold=True, font=HEAD)
text(s, 0.60, 0.78, 12.2, 0.9, "Status, and where it goes next", size=30, color=WHITE, bold=True, font=HEAD)
# left: takeaways
PH = 3.78
rect(s, 0.62, 1.95, 6.05, PH, "152B4A", rounded=True, radius=0.04)
text(s, 0.9, 2.2, 5.5, 0.4, "Key takeaways", size=17, color=SKY, bold=True, font=HEAD)
paras(s, 0.92, 2.85, 5.55, PH - 0.9,
      [{"text": "Rayleigh calibration is mature across CHM15k, CL61 and Mini-MPL.", "bullet": True, "bcolor": SKY, "size": 14, "color": "E8EEF6", "space_after": 14, "line_spacing": 1.1},
       {"runs": [("CL61 is the most precise group ", {"color": "E8EEF6"}), ("(median σ_SD 7.7%).", {"color": ICE})],
        "bullet": True, "bcolor": SKY, "size": 14, "space_after": 14, "line_spacing": 1.1},
       {"runs": [("Mini-MPL calibrates to CHM15k class ", {"color": "E8EEF6"}), ("(σ_SD 8.5%).", {"color": ICE})],
        "bullet": True, "bcolor": GREEN, "size": 14, "space_after": 14, "line_spacing": 1.1},
       {"text": "E-PROF v2 is the most precise method (σ_SD 9.5%).", "bullet": True, "bcolor": GREEN, "size": 14, "color": "E8EEF6", "space_after": 14, "line_spacing": 1.1},
       {"runs": [("Sign error fixed (v1.0→v1.1): ", {"color": "E8EEF6"}),
                 ("corrupts β_att, not Cₗ stability; +8–18% at 532 nm.", {"color": ICE})],
        "bullet": True, "bcolor": RED, "size": 14, "line_spacing": 1.1}])
# right: recommendations
rect(s, 6.83, 1.95, 5.85, PH, "152B4A", rounded=True, radius=0.04)
text(s, 7.1, 2.2, 5.4, 0.4, "Next steps", size=17, color=GREEN, bold=True, font=HEAD)
paras(s, 7.12, 2.85, 5.4, PH - 0.9,
      [{"text": "Upgrade operations from v1.0 to E-PROF v2.", "bullet": True, "bcolor": GREEN, "size": 14, "color": "E8EEF6", "space_after": 14, "line_spacing": 1.1},
       {"text": "Deploy CL61 Rayleigh calibration operationally.", "bullet": True, "bcolor": GREEN, "size": 14, "color": "E8EEF6", "space_after": 14, "line_spacing": 1.1},
       {"text": "Calibrate CL31 / CL51 via the cloud method.", "bullet": True, "bcolor": GREEN, "size": 14, "color": "E8EEF6", "space_after": 14, "line_spacing": 1.1},
       {"text": "Feed per-night constants to the Kalman filter.", "bullet": True, "bcolor": GREEN, "size": 14, "color": "E8EEF6", "space_after": 14, "line_spacing": 1.1},
       {"text": "Finish the E-PROFILE ALC paper figures.", "bullet": True, "bcolor": GREEN, "size": 14, "color": "E8EEF6"}])
# closing banner
rect(s, 0.62, 6.02, 12.06, 0.82, GREEN, rounded=True, radius=0.14)
text(s, 0.62, 6.02, 12.06, 0.82,
     "Three wavelengths, one sign-corrected framework — E-PROF v2 + CL61 + Mini-MPL, ready to deploy.",
     size=15, color=INK, bold=True, font=HEAD, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
footer(s, 20, dark=True)

out = os.path.join(HERE, "ALC_calibration_status_2026.pptx")
prs.save(out)
print("saved", out, "with", len(prs.slides._sldIdLst), "slides")
