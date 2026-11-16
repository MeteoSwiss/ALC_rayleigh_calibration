"""Build the E-PROFILE ALC intercomparison-progression deck (Payerne) with python-pptx.

Reuses the "night-sky lidar" MeteoSwiss theme + helpers from build_deck.py. Each processing step
(Rayleigh -> cloud -> multiple-scattering -> water vapour -> wavelength -> overlap -> electronic
/dark -> SNR -> T/p/q source -> Cloudnet -> OmB -> dashboard) gets a slide with the challenge, what
it changes, the numbers, and the supporting intercomparison figure. Figures live in figs_intercomp/.
All numbers come from doc/reports/05_attbsc_validation.md and the companion reports 01-10.
"""
import os
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from PIL import Image

# ----------------------------------------------------------------- palette (night-sky lidar)
INK    = "0E1B33"
NAVY   = "13294B"
SKY    = "2E7DA1"   # CL61 / teal
GREEN  = "2FBF71"   # cloud / "good"
GREEND = "1E7A48"
AMBER  = "E9A23B"   # CL31 / 910 nm ceilometer
RED    = "C0392B"   # CHM15k / reference
LIGHT  = "F5F7FB"
CARD   = "FFFFFF"
CARDLN = "E2E8F0"
MUTED  = "5A6B82"
INKTXT = "1A2A44"
WHITE  = "FFFFFF"
ICE    = "AEBED2"

# instrument colours (match the validation figures)
CHM_C, CL61_C, CL31_C = RED, SKY, AMBER

HEAD = "Trebuchet MS"
BODY = "Calibri"
MONO = "Consolas"

HERE = os.path.dirname(os.path.abspath(__file__))
FIG  = os.path.join(HERE, "figs_intercomp")


def F(name):
    return os.path.join(FIG, name)


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
    rect(s, x - pad, y - pad, w + 2 * pad, h + 2 * pad, CARD, line=CARDLN, rounded=True, radius=0.03)
    img_fit(s, path, x, y, w, h)


SECTION = {"method": SKY, "rayleigh": SKY, "cloud": GREEN, "corr": AMBER, "hw": RED,
           "product": NAVY, "close": GREEN}


def accent_bar(s, color, w=0.14):
    rect(s, 0, 0, w, SH, color)


def instr_dots(s, x, y, dark=True, gap=2.2, d=0.16, size=10.5):
    trio = [("CHM15k · 1064 nm", CHM_C), ("CL61 · 910 nm", CL61_C), ("CL31 · 910 nm", CL31_C)]
    cx = x
    for lbl, col in trio:
        oval(s, cx, y, d, col)
        text(s, cx + d + 0.07, y - 0.05, 1.9, 0.3, lbl, size=size,
             color=ICE if dark else MUTED, bold=True, font=BODY, anchor=MSO_ANCHOR.MIDDLE)
        cx += gap


def footer(s, idx, dark=False):
    c = "55698A" if dark else "9AA9BE"
    text(s, 0.62, 7.08, 10, 0.3,
         "E-PROFILE ALC calibration  ·  Payerne intercomparison  ·  M. Hervo, MeteoSwiss  ·  July 2026",
         size=8.5, color=c)
    text(s, 12.1, 7.06, 0.7, 0.3, f"{idx:02d}", size=10, color=c, bold=True, align=PP_ALIGN.RIGHT)


def header(s, kicker, title, section, idx, title_size=27):
    set_bg(s, LIGHT)
    accent = SECTION[section]
    accent_bar(s, accent)
    text(s, 0.62, 0.36, 11.9, 0.3, kicker.upper(), size=12.5, color=accent, bold=True, font=HEAD)
    text(s, 0.60, 0.66, 12.4, 1.0, title, size=title_size, color=INK, bold=True, font=HEAD)
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
    instr_dots(s, 0.95, 6.35, dark=True)
    footer(s, idx, dark=True)


def chip(s, x, y, w, h, label, fill, txtcolor=WHITE, size=12.5):
    rect(s, x, y, w, h, fill, rounded=True, radius=0.22)
    text(s, x, y, w, h, label, size=size, color=txtcolor, bold=True, font=HEAD,
         align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)


def stat_card(s, x, y, w, h, number, label, color, num_size=30):
    rect(s, x, y, w, h, CARD, line=CARDLN, rounded=True, radius=0.07)
    rect(s, x, y, 0.10, h, color)
    text(s, x + 0.22, y + 0.12, w - 0.34, 0.7, number, size=num_size, color=color, bold=True, font=HEAD)
    paras(s, x + 0.24, y + 0.80, w - 0.42, h - 0.85,
          [{"text": label, "size": 11, "color": MUTED, "line_spacing": 1.0}])


def side_text(s, items, x=0.66, y=1.95, w=5.55, h=4.9, size=13.5):
    """Left column of challenge/what/numbers bullets on a content slide."""
    paras(s, x, y, w, h, items, size=size)


def fig_caption(s, path, x, y, w, h, cap, capcolor=MUTED, pad=0.12):
    img_frame(s, path, x, y, w, h, pad=pad)
    text(s, x, y + h + 0.14, w, 0.5, cap, size=10.5, color=capcolor, italic=True, align=PP_ALIGN.CENTER)


# ===================================================================== 1 TITLE
s = slide()
set_bg(s, INK)
rect(s, 11.62, 1.4, 0.06, 4.7, "26527A")
rect(s, 11.95, 0.9, 0.06, 5.7, "2F6E5B")
rect(s, 12.28, 2.0, 0.06, 3.6, "B07A2E")
accent_bar(s, GREEN, w=0.24)
text(s, 0.95, 1.45, 11, 0.4, "METEOSWISS  ·  E-PROFILE  ·  AUTOMATIC LIDARS & CEILOMETERS",
     size=13, color=SKY, bold=True, font=HEAD)
paras(s, 0.92, 1.95, 11.4, 2.4,
      [{"runs": [("Calibrating the ALC network", {})], "size": 46, "color": WHITE, "bold": True,
        "font": HEAD, "space_after": 2, "line_spacing": 1.0},
       {"runs": [("one correction at a time", {})], "size": 46, "color": GREEN, "bold": True,
        "font": HEAD, "line_spacing": 1.0}])
paras(s, 0.95, 4.25, 10.6, 1.4,
      [{"runs": [("A step-by-step tour of the ", {"color": ICE}),
                 ("Payerne intercomparison", {"color": WHITE, "bold": True}),
                 (" — three co-located ceilometers,", {"color": ICE})],
        "size": 17, "line_spacing": 1.2, "space_after": 2},
       {"runs": [("and what each calibration step changes:  ", {"color": ICE}),
                 ("from ±20 % scatter to ±3 % agreement.", {"color": AMBER, "bold": True})],
        "size": 17, "line_spacing": 1.2}])
instr_dots(s, 0.97, 5.9, dark=True)
text(s, 0.95, 6.85, 11, 0.3, "M. Hervo, MeteoSwiss  ·  July 2026  ·  eprofile-calibration v2.0  ·  L1 uniform validation",
     size=12, color="6E83A6", bold=True)

# ===================================================================== 2 THE TEST-BENCH
s = slide()
header(s, "The Payerne test-bench", "Three co-located instruments, one reference", "method", 2)
side_text(s, [
    {"runs": [("Payerne (06610): ", {"bold": True, "color": INK}),
              ("CHM15k 1064 nm (molecular reference), CL31 & CL61 at 910 nm — same air, same time.",
               {})], "bullet": True, "bcolor": SKY, "space_after": 10},
    {"text": "Each channel is calibrated INDEPENDENTLY from the native L1 signal (no L2 product, no "
             "vendor constant) — so the comparison tests the calibration itself.",
     "bullet": True, "bcolor": SKY, "space_after": 10},
    {"runs": [("Compared pairwise vs the CHM15k over ", {}),
              ("500–3000 m AGL", {"bold": True, "color": INK}),
              (", Mar–Jun 2026 (N ≈ 150–230 k gate-pairs).", {})],
     "bullet": True, "bcolor": SKY, "space_after": 10},
    {"runs": [("Headline metric: ", {}), ("median relative bias", {"bold": True, "color": INK}),
              (" + ", {}), ("log-space r", {"bold": True, "color": INK}),
              (" (robust across two decades of β).", {})],
     "bullet": True, "bcolor": SKY, "space_after": 10},
    {"runs": [("Panels: ", {"bold": True, "color": INK}),
              ("(a) median±IQR profiles · (b) scatter vs ref · (c) difference pdf · "
               "(d–g) curtains (grey = flagged, colour = kept).", {})],
     "bullet": True, "bcolor": SKY},
])
fig_caption(s, F("fig_payerne.png"), 6.45, 1.95, 6.35, 4.55,
            "The endpoint, every correction ON — CL61 +0.1 % (cloud) / −0.6 % (Rayleigh), log r 0.94.  "
            "We now rewind and add the corrections one by one.")

# ===================================================================== 3 THE LADDER
s = slide()
header(s, "What this deck covers", "The correction ladder", "method", 3)
ladder = [
    ("Rayleigh v2", "temporal aerosol rejection in the molecular fit", "σ_SD 11.6 → 9.5 %", SKY),
    ("Liquid-cloud cal.", "calibrate the CL31/CL51 that Rayleigh can't reach", "yield ×4–9", GREEN),
    ("Multiple scattering", "per-instrument PVC η tables (a_G = 5.5 µm)", "CL61 gap 12 % → <1 %", GREEN),
    ("Water vapour", "two-way 910 nm transmission (CAMS L137)", "+14–35 % on β", AMBER),
    ("Wavelength → 1064", "molecular-aware, not a single Ångström", "CL61 +21 % → +0.1 %", AMBER),
    ("Overlap (T model)", "temperature-dependent near-range (Hervo 2016)", "<0.1 % in-band", SKY),
    ("Overlap prototype", "reconstruct true O(z) to 1 km (noise / neighbours)", "ongoing", MUTED),
    ("Electronic / dark", "hood-measured detector background", "CL31 pedestal exposed", RED),
    ("SNR filter", "none vs SNR≥3 vs Cloudnet", "keep 'none'", RED),
    ("T/p/q source", "CAMS 0.4° / 1° / ERA5 / sonde", "~1 % flat; valleys dominate", NAVY),
]
my = 1.95
rhh = 0.44
for i, (step, desc, eff, col) in enumerate(ladder):
    y = my + i * (rhh + 0.055)
    rect(s, 0.66, y, 11.98, rhh, CARD, line=CARDLN, rounded=True, radius=0.16)
    oval(s, 0.80, y + 0.10, rhh - 0.20, col)
    text(s, 0.72, y, 0.5, rhh, str(i), size=11, color=WHITE, bold=True, font=HEAD,
         align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
    text(s, 1.35, y, 2.7, rhh, step, size=13, color=INK, bold=True, font=HEAD, anchor=MSO_ANCHOR.MIDDLE)
    text(s, 4.15, y, 5.6, rhh, desc, size=11.5, color=MUTED, anchor=MSO_ANCHOR.MIDDLE)
    chip(s, 9.95, y + 0.06, 2.55, rhh - 0.12, eff, col, size=10.5)
text(s, 0.66, 6.95, 12, 0.3,
     "Each row below gets its own slide: the challenge, what it changes, and the supporting figure.",
     size=11, color=MUTED, italic=True)

# ===================================================================== 4 DIVIDER — CALIBRATION
s = slide()
divider(s, 4, "Part 1", "Absolute calibration", "rayleigh",
        "Molecular (Rayleigh) and liquid-cloud calibration, and the multiple-scattering fix that "
        "reconciles them.", "01")

# ===================================================================== 5 STEP 0 RAYLEIGH ONLY
s = slide()
header(s, "Rung 0 · Before this work", "The starting point: the operational L2 product", "hw", 5)
side_text(s, [
    {"runs": [("This is the E-PROFILE ", {}),
              ("L2 attenuated-backscatter product exactly as distributed today", {"bold": True, "color": INK}),
              (" — no water-vapour, no wavelength, no overlap harmonisation.", {})],
     "bullet": True, "bcolor": RED, "space_after": 11},
    {"runs": [("Each instrument carries its operational constant: ", {}),
              ("CHM15k the Rayleigh v1.0 calibration; CL61 Vaisala's internal constant; "
               "CL31 the 1·10⁸ default — effectively uncalibrated.", {"bold": True, "color": INK})],
     "bullet": True, "bcolor": RED, "space_after": 11},
    {"runs": [("Result: ", {"bold": True, "color": INK}),
              ("CL31 reads ", {}), ("−67 %", {"bold": True, "color": RED}),
              (" (uncalibrated) and CL61 ", {}), ("+33 %", {"bold": True, "color": RED}),
              (" (Vaisala constant, no WV / wavelength). Neither 910 nm ceilometer is on the "
               "CHM15k scale.", {})],
     "bullet": True, "bcolor": RED, "space_after": 11},
    {"runs": [("→ ", {"bold": True}),
              ("Everything that follows is what our L1 calibration + harmonisation add to this baseline.",
               {"bold": True, "color": INK})],
     "bullet": True, "bcolor": AMBER},
])
fig_caption(s, F("ladder/00_before_L2.png"), 6.35, 1.95, 6.5, 4.5,
            "Rung 0 — the L2 product AS DISTRIBUTED, no harmonisation: CL31 −67 % (orange, too low) "
            "and CL61 +33 % (grey, too high) vs the CHM15k reference (red).")

# ===================================================================== 6 RAYLEIGH v2
s = slide()
header(s, "Step 1 · Better molecular selection", "Rayleigh v2: reject aerosol that looks molecular", "rayleigh", 6)
side_text(s, [
    {"runs": [("Challenge: ", {"bold": True, "color": RED}),
              ("a smooth aerosol layer can look linear in the night-mean and be mistaken for "
               "molecular → biased lidar constant.", {})],
     "bullet": True, "bcolor": SKY, "space_after": 12},
    {"runs": [("v2 uses the ", {}), ("full-night time series", {"bold": True, "color": INK}),
              (": molecular scattering is steady (only tracks T/p) while aerosol advects and "
               "fluctuates. Temporally-variable windows are rejected; contaminated time-height "
               "cells are flagged and the fit uses the cleaned mean.", {})],
     "bullet": True, "bcolor": SKY, "space_after": 12},
    {"text": "The gate relaxation (config C8) also roughly doubles the usable clear-night fraction.",
     "bullet": True, "bcolor": SKY},
], h=3.0)
cards = [("9.5 %", "σ_SD (drift-insensitive precision)\nvs 11.6 % for v1.1", GREEN),
         ("+34 pp", "CHM15k valid-night yield\n38 → 72 %", SKY),
         ("+30 pp", "CL61 valid-night yield\n45 → 75 %", SKY),
         ("6", "selectable molecular-window\nmethods (v0.25 → v2)", NAVY)]
for i, (no, lab, col) in enumerate(cards):
    cx = 6.55 + (i % 2) * 3.15
    cy = 1.95 + (i // 2) * 1.75
    stat_card(s, cx, cy, 2.95, 1.55, no, lab.replace("\n", " "), col)
text(s, 6.55, 5.55, 6.2, 1.3,
     "Note: the intercomparison bias is the same under v1 or v2 (CHM15k and CL61 agree either way) — "
     "molecular selection changes which nights pass and how tightly, i.e. precision and yield, not the "
     "mean bias. That is why this rung has no distinct intercomparison figure of its own.",
     size=11.5, color=MUTED, italic=True)

# ===================================================================== 7 + CLOUD
s = slide()
header(s, "Step 2 · Liquid-cloud calibration", "Cloud calibration fills the CL31 / CL51 gap", "cloud", 7)
side_text(s, [
    {"runs": [("O'Connor / Hopkin: ", {"bold": True, "color": INK}),
              ("an opaque liquid cloud fully attenuates the beam → the integrated β is fixed by a "
               "theoretical lidar ratio S = 18.8 sr.", {})],
     "bullet": True, "bcolor": GREEN, "space_after": 10},
    {"runs": [("Calibrates the ceilometers Rayleigh can't reach ", {}),
              ("(CL31/CL51, otherwise stuck at the 1·10⁸ default).", {"color": INK, "bold": True})],
     "bullet": True, "bcolor": GREEN, "space_after": 10},
    {"runs": [("Fixed ", {}), ("100–2400 m integration gate", {"bold": True, "color": INK}),
              (" (not the reported cloud base) → immune to Vaisala CBH bias.", {})],
     "bullet": True, "bcolor": GREEN, "space_after": 10},
    {"runs": [("Native ", {}), ("30 s cadence", {"bold": True, "color": INK}),
              (" → valid-day yield ×4–9 (CL31 16 → 53 %, CL51 11 → 45 %).", {})],
     "bullet": True, "bcolor": GREEN},
])
fig_caption(s, F("fig_cloud_calibration_example.png"), 6.4, 1.95, 6.45, 4.5,
            "One night, one ceilometer: cloud detection, the 0.1–2.4 km integral, WV correction, "
            "apparent S ≈ 16 sr, and the calibration-coefficient histogram.")

# ===================================================================== 8 MULTIPLE SCATTERING
s = slide()
header(s, "Step 3 · Multiple scattering", "Each ceilometer gets its own MS table", "cloud", 8)
side_text(s, [
    {"runs": [("η ≤ 1 in ", {}), ("B = 1/(2·η·S)", {"font": MONO, "bold": True, "color": INK}),
              (" absorbs the multiple-scattering enhancement. The legacy flat 8 µm ladder made "
               "CL31 ≈ CL51 and forced a spurious CL51 factor onto the CHM15k.", {})],
     "bullet": True, "bcolor": GREEN, "space_after": 11},
    {"runs": [("Now: reference-exact ", {}),
              ("Hogan-2006 PVC tables at the calibration-scene droplet size a_G = 5.5 µm", {"bold": True, "color": INK}),
              (", one per receiver FOV. CL61 no longer borrows CL51's.", {})],
     "bullet": True, "bcolor": GREEN, "space_after": 11},
    {"runs": [("η ≈ 0.95 (was 0.83) at low cloud base → the cloud-base tilt in C drops 3–4×.", {})],
     "bullet": True, "bcolor": GREEN, "space_after": 11},
    {"runs": [("Regenerated intercomparison — with the pre-2026-07 cloud cal the ", {}),
              ("CL61 cloud reads −18 % and splits from its own Rayleigh (−1 %)", {"bold": True, "color": RED}),
              (";  the recalibration is the next slide.", {})],
     "bullet": True, "bcolor": SKY},
])
fig_caption(s, F("ladder/01_add_cloud_legacyMS.png"), 6.5, 1.95, 6.35, 4.5,
            "Rung 1 (regenerated) — cloud channels added with the legacy 8 µm MS table: CL61 cloud "
            "−18 % (scatter, panel b) splits from CL61 Rayleigh −1 %.")

# ===================================================================== 9 CONVERGENCE
s = slide()
header(s, "Step 3 · Result", "CL61's two calibrations now agree", "cloud", 9)
side_text(s, [
    {"text": "The multiple-scattering fix removes the last free parameter that split the two methods.",
     "bullet": True, "bcolor": GREEN, "space_after": 12},
    {"runs": [("Payerne CL61: ", {"bold": True, "color": INK}),
              ("C_L(cloud) 1.207 vs C_L(Rayleigh) 1.216 — < 1 %.", {})],
     "bullet": True, "bcolor": GREEN, "space_after": 12},
    {"runs": [("The β rows: ", {}), ("cloud +0.1 % vs Rayleigh −0.6 %", {"bold": True, "color": GREEND}),
              (" — 0.7 points apart.", {})],
     "bullet": True, "bcolor": GREEN, "space_after": 12},
    {"runs": [("Two independent calibrations — ", {}),
              ("strong-signal cloud and weak-signal molecular — land on the same lidar constant", {"bold": True, "color": INK}),
              (". A self-consistency proof, not a coincidence.", {})],
     "bullet": True, "bcolor": SKY},
])
fig_caption(s, F("ladder/02_per_instrument_MS.png"), 6.5, 1.95, 6.35, 4.5,
            "Rung 2 (regenerated) — with the per-instrument MS tables the CL61 cloud and Rayleigh rows "
            "sit on top of each other and on the CHM15k.")

# ===================================================================== 10 DIVIDER — HARMONISATION
s = slide()
divider(s, 10, "Part 2", "Harmonisation\ncorrections", "corr",
        "Making a 910 nm ceilometer comparable to a 1064 nm reference: water vapour, wavelength, and "
        "overlap.", "02")

# ===================================================================== 11 WATER VAPOUR
s = slide()
header(s, "Step 4 · Water vapour", "The 910 nm water-vapour correction", "corr", 11)
side_text(s, [
    {"runs": [("CL31/CL51/CL61 emit in the ", {}), ("905–911 nm H₂O absorption band", {"bold": True, "color": INK}),
              ("; the 1064 nm CHM15k does not. Their β is attenuated by water vapour the reference "
               "never sees.", {})],
     "bullet": True, "bcolor": AMBER, "space_after": 11},
    {"runs": [("Two-way ", {}), ("T²_wv", {"font": MONO, "bold": True, "color": INK}),
              (" from CAMS L137 humidity × the measured laser spectrum; β ÷ T²_wv. ", {}),
              ("Mandatory", {"bold": True, "color": RED}), (" — a no-WV mode is rejected.", {})],
     "bullet": True, "bcolor": AMBER, "space_after": 11},
    {"runs": [("In-band β increase ", {}), ("+14–35 %", {"bold": True, "color": INK}),
              (" (Payerne CL61 +30 %, CL31 +23 %).", {})],
     "bullet": True, "bcolor": AMBER, "space_after": 11},
    {"runs": [("Built-in check: against a 910 nm reference the correction cancels (Uccle Δ +1.1 %) — "
               "it is a real spectral effect, not a fudge.", {})],
     "bullet": True, "bcolor": SKY},
])
fig_caption(s, F("ladder/03_no_water_vapour.png"), 6.5, 1.95, 6.35, 4.5,
            "Regenerated with WV OFF — every 910 nm channel collapses to ~−22 % (CL61 −22.5 %, CL31 "
            "−22.8 %). Turning WV ON lifts them back to agreement (+22.6 / +14.9 pts).")

# ===================================================================== 12 WAVELENGTH SIMPLE
s = slide()
header(s, "Step 5 · Wavelength conversion", "The naive way: a single Ångström exponent", "corr", 12)
side_text(s, [
    {"runs": [("To compare 910 vs 1064 nm you must convert. The simple approach scales the ", {}),
              ("whole signal", {"bold": True, "color": INK}), (" by (910/1064)^α.", {})],
     "bullet": True, "bcolor": AMBER, "space_after": 12},
    {"runs": [("Problem: the molecular part follows ", {}), ("λ⁻⁴", {"font": MONO, "bold": True, "color": RED}),
              (" (ratio 1.877), not λ⁻ᵅ. A flat α leaves it ~1.6× over-scaled.", {})],
     "bullet": True, "bcolor": AMBER, "space_after": 12},
    {"runs": [("Result: CL61 sits ", {}), ("+21 %", {"bold": True, "color": RED}),
              (" high at Payerne (+32 Lindenberg, +28 Camborne), ", {}),
              ("growing with altitude", {"bold": True, "color": INK}),
              (" and at night / in winter — the fingerprint of a molecular mis-scale, not a "
               "calibration error.", {})],
     "bullet": True, "bcolor": AMBER},
], h=3.2)
fig_caption(s, F("ladder/05_wavelength_angstrom.png"), 6.5, 1.95, 6.35, 4.5,
            "Rung 4 (regenerated) — single-Ångström conversion: the CL61 rides high and the difference "
            "histogram is displaced positive.")

# ===================================================================== 13 WAVELENGTH MOLECULAR-AWARE
s = slide()
header(s, "Step 5 · Wavelength conversion", "The right way: separate molecular & aerosol", "corr", 13)
side_text(s, [
    {"runs": [("Component-separated conversion: ", {"bold": True, "color": INK}),
              ("compute the analytic molecular β from CAMS p/T, remove it, apply Ångström (α = 1) to "
               "the aerosol residual ONLY, then re-add the molecular at 1064 nm.", {})],
     "bullet": True, "bcolor": GREEN, "space_after": 12},
    {"runs": [("CL61 → ", {}), ("+0.1 % (Payerne), −1.7 % (Lindenberg), −2.9 % (Camborne)", {"bold": True, "color": GREEND}),
              (": ±3 % agreement.", {})],
     "bullet": True, "bcolor": GREEN, "space_after": 12},
    {"runs": [("The altitude tilt collapses ", {}),
              ("(Lindenberg across-band spread 12.4 → 0.2 points)", {"bold": True, "color": INK}),
              (" — something no multiplicative calibration-scale fix can do.", {})],
     "bullet": True, "bcolor": GREEN, "space_after": 12},
    {"runs": [("This is the paper's central methodological result.", {"italic": True, "color": MUTED})]},
], h=3.4)
fig_caption(s, F("ladder/07_hood_dark.png"), 6.5, 1.95, 6.35, 4.5,
            "Fully-corrected (regenerated) — molecular-aware conversion: CL61 +0.1 %, difference "
            "histogram on zero. (Light-blue = the CL31 hood-dark twin, slide 18.)")

# ===================================================================== 14 OVERLAP T MODEL
s = slide()
header(s, "Step 6 · Near range", "Overlap: the operational temperature model", "corr", 14)
side_text(s, [
    {"runs": [("The CHM15k near range needs an overlap correction; the residual is ", {}),
              ("temperature-dependent", {"bold": True, "color": INK}),
              ("  (Hervo 2016): rcs_0 × (1 + (a(z)·T + b(z))/100).", {})],
     "bullet": True, "bcolor": SKY, "space_after": 11},
    {"runs": [("No-op above ~700 m: in the 500–3000 m band it moves the statistics ", {}),
              ("< 0.1 %", {"bold": True, "color": INK}),
              (" — it matters for the near-range PRODUCT, not the intercomparison numbers.", {})],
     "bullet": True, "bcolor": SKY, "space_after": 11},
    {"runs": [("At 25 °C: ", {}), ("+59 % @ 200 m, +14 % @ 300 m, +2 % @ 400 m.", {"bold": True, "color": INK})],
     "bullet": True, "bcolor": SKY, "space_after": 11},
    {"text": "Unit-specific in shape AND sign below ~400 m — Amsterdam's four co-located units make "
             "this visible.",
     "bullet": True, "bcolor": SKY},
])
fig_caption(s, F("fig_amsterdam_overlap_impact.png"), 6.5, 1.95, 6.35, 4.5,
            "Four co-located CHM15k: unit-specific overlap factors, corrected vs uncorrected medians, "
            "and the realised near-range impact.")

# ===================================================================== 15 OVERLAP PROTOTYPE
s = slide()
header(s, "Step 7 · Ongoing", "Correcting the overlap higher up", "corr", 15)
side_text(s, [
    {"runs": [("The factory reference overlap (TUB120011) completes ~750 m, but real overlaps reach ", {}),
              ("~1500 m", {"bold": True, "color": RED}), (" → a near-range tilt survives the whole chain.", {})],
     "bullet": True, "bcolor": AMBER, "space_after": 11},
    {"runs": [("Prototype: reconstruct the ", {}), ("true O(z) up to 1 km", {"bold": True, "color": INK}),
              (" from a clear-night molecular baseline fitted INSIDE a homogeneous residual layer "
               "(STRATFINDER picks the scene; night residual-layer median 764 m, p90 1214 m).", {})],
     "bullet": True, "bcolor": AMBER, "space_after": 11},
    {"runs": [("Two independent reconstructions — CL61-tilt and CHM-only — ", {}),
              ("coincide at 600 m (0.650 vs 0.649)", {"bold": True, "color": GREEND}), (".", {})],
     "bullet": True, "bcolor": AMBER},
])
fig_caption(s, F("12_three_overlap_functions.png"), 6.5, 1.95, 6.35, 4.5,
            "Three overlap functions: the (wrong) generic factory file, the onboard overlap from hood "
            "noise, and the true overlap reconstructed from the CL61 tilt.")

# ===================================================================== 16 OVERLAP WITHOUT A MODEL
s = slide()
header(s, "Step 7 · Ongoing", "Overlap without a factory model — noise & neighbours", "corr", 16)
side_text(s, [
    {"runs": [("At multi-instrument sites the CHM15k overlap need not come from a factory .cfg: ", {}),
              ("the co-located CL61 + a terminal-hood dark reconstruct the true O(z) directly", {"bold": True, "color": INK}),
              ("  (the Payerne natural-experiment triangle).", {})],
     "bullet": True, "bcolor": SKY, "space_after": 12},
    {"runs": [("For single-instrument sites, the overlap is recovered from ", {}),
              ("clear-sky NOISE", {"bold": True, "color": INK}),
              (" (electronic noise ∝ 1/O) — aerosol-immune and 2× tighter than regression "
               "(RMSE 0.016 vs 0.030 vs the hood truth).", {})],
     "bullet": True, "bcolor": SKY, "space_after": 12},
    {"runs": [("Runs network-wide (132 CHM15k) with a day-to-day variability QC flag.", {})],
     "bullet": True, "bcolor": SKY},
], h=3.4)
fig_caption(s, F("12_network_spectral_vs_m2.png"), 6.5, 1.95, 6.35, 4.5,
            "Network spectral-noise overlap vs the regression method across 132 CHM15k: the median is "
            "invariant, the per-station aerosol bias removed.")

# ===================================================================== 17 DIVIDER — HARDWARE & SAMPLING
s = slide()
divider(s, 17, "Part 3", "Hardware, sampling\n& model input", "hw",
        "Electronic distortions, the SNR filter, and the temperature/humidity profile source — the "
        "things that bite below the calibration.", "03")

# ===================================================================== 18 ELECTRONIC / DARK
s = slide()
header(s, "Step 8 · Hardware", "Electronic distortions & the dark measurement", "hw", 18)
side_text(s, [
    {"runs": [("Covered-telescope (hood) sessions at Payerne isolate the ", {}),
              ("detector background P = rcs_0/z²", {"bold": True, "color": INK}), (" (the dark).", {})],
     "bullet": True, "bcolor": RED, "space_after": 10},
    {"runs": [("Per SAMPLE it's symmetric noise (median ~0); the ", {}),
              ("per-GATE median exposes an instrument-specific distortion", {"bold": True, "color": INK}),
              (" — invisible in a single profile.", {})],
     "bullet": True, "bcolor": RED, "space_after": 10},
    {"runs": [("CL31: two under-damped resonances (142 kHz + 30 kHz, R² 0.98) — a large near-range "
               "pedestal. CL61: AC-coupling undershoot. CHM15k: a small relaxation.", {})],
     "bullet": True, "bcolor": RED, "space_after": 10},
    {"runs": [("The CL31's native '+21 %' was partly its own dark pedestal inflating β. ", {}),
              ("The cloud method is offset-immune", {"bold": True, "color": GREEND}),
              (" (strong signal ≫ dark).", {})],
     "bullet": True, "bcolor": GREEN},
])
fig_caption(s, F("fig_dark_binscatter_payerne.png"), 6.4, 2.0, 6.5, 4.2,
            "NEW — per-sample β/r² binscatter (top) and per-gate median (bottom) for the three "
            "instruments. CL31's electronic distortion is ~20–100× the others'.")

# ===================================================================== 19 SNR FILTER
s = slide()
header(s, "Step 9 · Sampling", "SNR filtering: keep 'none'", "hw", 19)
side_text(s, [
    {"runs": [("Three modes tested: ", {}), ("none / SNR ≥ 3 / Cloudnet per-profile", {"bold": True, "color": INK}), (".", {})],
     "bullet": True, "bcolor": RED, "space_after": 11},
    {"runs": [("The Cloudnet per-profile filter conditions on positive noise excursions before "
               "averaging → ", {}), ("destructive for a bias comparison", {"bold": True, "color": RED}),
              (":", {})],
     "bullet": True, "bcolor": RED, "space_after": 6},
    {"runs": [("good CL61  +0.1 % → −10.5 %  (log r 0.94 → 0.76);  noisy CL31  −7.9 % → +85.6 %", {"bold": True, "color": INK})],
     "level": 1, "space_after": 11, "size": 12.5},
    {"runs": [("It keeps only ~3.7 % of the CL31 samples (199 k → 7.5 k).", {})],
     "bullet": True, "bcolor": RED, "space_after": 11},
    {"runs": [("'none' preserves the ", {}),
              ("symmetric-noise-averages-out property", {"bold": True, "color": GREEND}),
              (" an unbiased instrument-vs-instrument comparison needs.", {})],
     "bullet": True, "bcolor": GREEN},
])
fig_caption(s, F("fig_payerne_noisefilter_comparison.png"), 6.5, 2.05, 6.35, 4.1,
            "Median relative bias per channel under the three modes — the Cloudnet filter swings the "
            "CL31 to +86 % and drags the clean CL61 negative.")

# ===================================================================== 20 T/P/Q SOURCE
s = slide()
header(s, "Step 10 · Model input", "T/p/q source: resolution vs orography", "product", 20)
side_text(s, [
    {"runs": [("The WV and molecular terms need a T/p/q profile. Tested: ", {}),
              ("CAMS 0.4°, CAMS 1°, ERA5 (19-level), radiosonde.", {"bold": True, "color": INK})],
     "bullet": True, "bcolor": NAVY, "space_after": 11},
    {"runs": [("Flat / lowland sites: ", {}), ("ΔC ≈ 1 % or better", {"bold": True, "color": GREEND}),
              (" — CAMS 0.4° is adequate (Payerne, Uccle, Köln…).", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 11},
    {"runs": [("The leading error is ", {}), ("orography, not resolution", {"bold": True, "color": RED}),
              (": valley sites where the model grid-cell altitude is wrong — Aosta's CAMS cell sits "
               "+1627 m above the station → several %.", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 11},
    {"runs": [("Temporal 1 h vs 3 h is negligible (RMS < 0.4 %).", {})],
     "bullet": True, "bcolor": SKY},
])
fig_caption(s, F("fig_wv_10stations_impact.png"), 6.5, 1.95, 6.35, 4.5,
            "Ten stations: grid-cell orography mismatch (left) and cloud-calibration error (right) — "
            "flat sites cluster < 1 %, valley sites spread wide.")

# ===================================================================== 21 DIVIDER — PRODUCTS
s = slide()
divider(s, 21, "Part 4", "Products &\nmonitoring", "product",
        "Beyond the calibration: target classification, an independent model cross-check, and the "
        "operational dashboard.", "04")

# ===================================================================== 22 CLOUDNET
s = slide()
header(s, "Cross-network", "Cloudnet target classification on E-PROFILE L1", "product", 22)
side_text(s, [
    {"runs": [("CloudnetPy-style ", {}), ("aerosol / cloud / ice / drizzle classification", {"bold": True, "color": INK}),
              (" ported to E-PROFILE L1 (CAMS temperature + attenuated backscatter).", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 11},
    {"runs": [("Per-instrument rcs_0 → β factors: ", {}),
              ("CHM15k ×3·10⁻¹², CL31 ×1·10⁻⁸, CL61 ×1.0.", {"font": MONO, "bold": True, "color": INK})],
     "bullet": True, "bcolor": NAVY, "space_after": 11},
    {"runs": [("Payerne 06610: ", {}), ("A = CHM15k, B = CL31, C = CL61.", {"bold": True, "color": INK})],
     "bullet": True, "bcolor": NAVY, "space_after": 11},
    {"text": "Feeds QC, adaptive detection thresholds and multi-instrument target agreement; runs as "
             "a standalone batch over ~433 streams.",
     "bullet": True, "bcolor": NAVY},
])
fig_caption(s, F("fig_cloudnet_classification_cl61.png"), 7.35, 1.62, 4.95, 5.15,
            "CL61 Payerne 2026-03-10: screened β, depolarization, target classification "
            "(aerosol / ice / supercooled / drizzle / liquid), and the cloud-aerosol β histogram.")

# ===================================================================== 23 OmB
s = slide()
header(s, "Independent check", "OmB: calibrated β vs the CAMS aerosol forecast", "product", 23)
side_text(s, [
    {"runs": [("Observation-minus-Background compares each ", {}),
              ("calibrated β to the CAMS aerosol forecast", {"bold": True, "color": INK}),
              (" (per wavelength; native at 1064/532, Ångström at 910; WV-corrected).", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 12},
    {"text": "A daily per-station bias / RMS product — detects window contamination, drift and "
             "calibration-scale errors independently of the co-located reference.",
     "bullet": True, "bcolor": NAVY, "space_after": 12},
    {"runs": [("The Kalman lidar constant removes the large CL31 negative bias: ", {}),
              ("−0.272 → −0.047 Mm⁻¹sr⁻¹.", {"bold": True, "color": GREEND})],
     "bullet": True, "bcolor": GREEN},
])
fig_caption(s, F("omb_CL31.png"), 6.5, 1.95, 6.35, 4.5,
            "Payerne CL31 OmB vs CAMS (910 nm, Ångström-interpolated + WV-corrected): the calibrated "
            "product tracks the model background through the boundary layer.")

# ===================================================================== 24 DASHBOARD
s = slide()
header(s, "Monitoring", "The operational monitoring dashboard", "product", 24)
side_text(s, [
    {"runs": [("Static ", {}), ("SQLite → Plotly → Jinja2", {"font": MONO, "bold": True, "color": INK}),
              (", rebuilt daily for ~436 stations on CSCS and published to the European Weather Cloud.", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 12},
    {"runs": [("Per-station: ", {"bold": True, "color": INK}),
              ("C_L time series (Rayleigh + cloud + Kalman), OmB, ICAO sensitivity — plus network maps.", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 12},
    {"runs": [("Live: ", {}),
              ("alc-calib.ch-meteoswiss-emermet.f.ewcloud.host", {"font": MONO, "bold": True, "color": SKY})],
     "bullet": True, "bcolor": NAVY},
], h=3.2)
fig_caption(s, F("dashboard_main.png"), 6.35, 2.0, 6.5, 4.35,
            "Live network summary: 433 instruments, 217 k calibrations — active instruments by type "
            "and two maps of median C_L vs the theoretical / operational value.")

# ===================================================================== 25 DASHBOARD STATION
s = slide()
header(s, "Monitoring", "… and one page per instrument", "product", 25)
side_text(s, [
    {"runs": [("Every stream has its own page: ", {}),
              ("data-availability & status calendar, the Rayleigh-vs-cloud C_L scatter, and the "
               "Kalman-smoothed series", {"bold": True, "color": INK}), (".", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 12},
    {"runs": [("Downloads: the calibration NetCDF and a link to the CEDA L2 data.", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 12},
    {"runs": [("Shown: ", {}), ("Payerne CL61", {"bold": True, "color": SKY}),
              (" — the two calibrations now sit together at C_L ≈ 1.25 (the convergence of slide 9, "
               "live).", {})],
     "bullet": True, "bcolor": NAVY},
], h=3.4)
fig_caption(s, F("dashboard_station_payerne.png"), 6.35, 2.0, 6.5, 4.35,
            "Payerne CL61 station page: availability calendar (green = pass), then the "
            "Rayleigh (blue) vs liquid-cloud (green) lidar-constant series.")

# ===================================================================== 26 EWC EMAIL NOTIFICATIONS
s = slide()
header(s, "Monitoring · in discussion", "Partner email notifications on the EWC VM", "product", 26)
side_text(s, [
    {"runs": [("Goal: replace the manual partner emails ", {}),
              ("(per-instrument issue alerts + a monthly dashboard digest)", {"color": MUTED}),
              (" with a self-serve subscribe / notify service.", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 10},
    {"runs": [("The alc-calib VM is a full VM behind nginx → it can host a tiny FastAPI + SQLite "
               "backend for ", {}),
              ("/subscribe /confirm /manage /unsubscribe", {"font": MONO, "color": INK}),
              (" — no login, double opt-in + capability-token links.", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 10},
    {"runs": [("Sending: outbound ", {}), ("SMTP 587 is BLOCKED", {"bold": True, "color": RED}),
              ("; 465 and 443 are open → use an ", {}),
              ("HTTPS transactional API (Brevo / Mailjet, EU/GDPR)", {"bold": True, "color": INK}),
              (", not raw SMTP.", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 10},
    {"runs": [("Subscription scope = ", {}),
              ("all | countries[] | instruments[]", {"font": MONO, "color": INK}),
              (" resolved from the census; the alert trigger reuses the dashboard watchlist "
               "(drift / failure-streak / low success).", {})],
     "bullet": True, "bcolor": NAVY, "space_after": 10},
    {"runs": [("One external dependency: ", {}),
              ("SPF / DKIM on the From-domain", {"bold": True, "color": AMBER}),
              (" (MeteoSwiss IT, once) for deliverability.", {})],
     "bullet": True, "bcolor": AMBER},
], x=0.66, y=1.9, w=12.0, h=4.9, size=14)
# a compact flow strip along the bottom
flow = [("Subscribe", SKY), ("Double opt-in", SKY), ("Watchlist trigger", AMBER),
        ("HTTPS email API", GREEN), ("Alert / monthly digest", GREEND)]
fx = 0.66
for i, (lbl, col) in enumerate(flow):
    chip(s, fx, 6.05, 2.15, 0.5, lbl, col, size=12)
    if i < len(flow) - 1:
        text(s, fx + 2.15, 6.05, 0.36, 0.5, "›", size=22, color=MUTED, bold=True,
             font=HEAD, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
    fx += 2.51

# ===================================================================== 27 LESSONS
s = slide()
header(s, "Beyond Payerne", "What co-location taught us", "hw", 27)
side_text(s, [
    {"runs": [("The chain validates across the network: ", {}),
              ("CHM15k units mutually ±5 %; CL61 vs CHM ±3 %", {"bold": True, "color": GREEND}),
              (" at Payerne / Lindenberg / Camborne; Magurele's two units agree to 1.7 % vs EARLINET.", {})],
     "bullet": True, "bcolor": SKY, "space_after": 12},
    {"runs": [("A bad unit hides in plain sight: ", {}),
              ("Amsterdam CHM15k-B reads +17 %", {"bold": True, "color": RED}),
              (" (daytime near-range) yet self-reports clean — only co-location caught it.", {})],
     "bullet": True, "bcolor": RED, "space_after": 12},
    {"runs": [("Practical QC: a ", {}), ("day − night bias split", {"bold": True, "color": INK}),
              (" is a cheap, calibration-independent detector (B: +14 pts vs ≤5 for healthy units).", {})],
     "bullet": True, "bcolor": AMBER},
], h=3.4)
fig_caption(s, F("fig_amsterdam_b_investigation.png"), 6.5, 1.95, 6.35, 4.5,
            "Amsterdam unit B, six probes: a smooth diurnal +32 % bias peaking at 15 UT, near-range "
            "weighted, stationary over the quarter — a signal-handling defect, not calibration.")

# ===================================================================== 26 SUMMARY
s = slide()
set_bg(s, INK)
accent_bar(s, GREEN, w=0.24)
text(s, 0.62, 0.45, 11.8, 0.3, "SUMMARY", size=13, color=GREEN, bold=True, font=HEAD)
text(s, 0.60, 0.78, 12.2, 0.9, "From ±20 % scatter to ±3 % agreement", size=28, color=WHITE, bold=True, font=HEAD)
PH = 3.95
rect(s, 0.62, 1.9, 6.05, PH, "152B4A", rounded=True, radius=0.04)
text(s, 0.9, 2.13, 5.5, 0.4, "Key takeaways", size=17, color=SKY, bold=True, font=HEAD)
paras(s, 0.92, 2.72, 5.55, PH - 0.9,
      [{"runs": [("CL61 calibrates two independent ways to the same constant ", {"color": "E8EEF6"}),
                 ("(< 1 %).", {"color": ICE})], "bullet": True, "bcolor": GREEN, "size": 13.5, "space_after": 12, "line_spacing": 1.1},
       {"runs": [("Cloud calibration + per-instrument MS finally calibrate the CL31 / CL51.", {"color": "E8EEF6"})],
        "bullet": True, "bcolor": GREEN, "size": 13.5, "space_after": 12, "line_spacing": 1.1},
       {"runs": [("Molecular-aware wavelength conversion is the central fix ", {"color": "E8EEF6"}),
                 ("(+21 % → +0.1 %).", {"color": ICE})], "bullet": True, "bcolor": AMBER, "size": 13.5, "space_after": 12, "line_spacing": 1.1},
       {"runs": [("Water vapour is mandatory at 910 nm; keep the noise filter 'none'.", {"color": "E8EEF6"})],
        "bullet": True, "bcolor": AMBER, "size": 13.5, "space_after": 12, "line_spacing": 1.1},
       {"runs": [("The CL31 is honest only below ~0.7 km — report a detection-limit profile.", {"color": "E8EEF6"})],
        "bullet": True, "bcolor": RED, "size": 13.5, "line_spacing": 1.1}])
rect(s, 6.83, 1.9, 5.85, PH, "152B4A", rounded=True, radius=0.04)
text(s, 7.1, 2.13, 5.4, 0.4, "Next steps", size=17, color=GREEN, bold=True, font=HEAD)
paras(s, 7.12, 2.72, 5.4, PH - 0.9,
      [{"text": "Roll v2 + cloud + WV into operations (from E-PROF v1.0).", "bullet": True, "bcolor": GREEN, "size": 13.5, "color": "E8EEF6", "space_after": 13, "line_spacing": 1.1},
       {"text": "Finish the overlap-from-noise prototype network-wide.", "bullet": True, "bcolor": GREEN, "size": 13.5, "color": "E8EEF6", "space_after": 13, "line_spacing": 1.1},
       {"text": "Fold Cloudnet classification into the daily runner.", "bullet": True, "bcolor": GREEN, "size": 13.5, "color": "E8EEF6", "space_after": 13, "line_spacing": 1.1},
       {"text": "Extend the hood / dark campaign beyond Payerne.", "bullet": True, "bcolor": GREEN, "size": 13.5, "color": "E8EEF6", "space_after": 13, "line_spacing": 1.1},
       {"text": "Publish the E-PROFILE ALC calibration paper.", "bullet": True, "bcolor": GREEN, "size": 13.5, "color": "E8EEF6"}])
rect(s, 0.62, 6.05, 12.06, 0.82, GREEN, rounded=True, radius=0.14)
text(s, 0.62, 6.05, 12.06, 0.82,
     "One L1 framework, three wavelengths, every correction explicit — CL61 to ±3 %, CL31/CL51 "
     "calibrated, ready to deploy.",
     size=14.5, color=INK, bold=True, font=HEAD, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
footer(s, 28, dark=True)

out = os.path.join(HERE, "ALC_Payerne_intercomparison_2026.pptx")
prs.save(out)
print("saved", out, "with", len(prs.slides._sldIdLst), "slides")
