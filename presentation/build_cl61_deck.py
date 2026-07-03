# -*- coding: utf-8 -*-
"""Presentation FR : offset electronique du CL61 et calibration Rayleigh.
Diagnostic hood -> modele physique (reponse a couplage AC) -> dependance temperature ->
correction / implications reseau. Figures = doc/reports/figs_paper_report/ (versionnees).
Rebuild: python presentation/build_cl61_deck.py
"""
from pptx import Presentation
from pptx.util import Inches as I, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

CHAR = RGBColor(0x36, 0x45, 0x4F); BLUE = RGBColor(0x1F, 0x77, 0xB4)
RED = RGBColor(0xD6, 0x27, 0x28); GREY = RGBColor(0x6E, 0x7B, 0x84)
DARK = RGBColor(0x22, 0x2E, 0x3A); WHITE = RGBColor(0xFF, 0xFF, 0xFF)
ORANGE = RGBColor(0xE1, 0x7A, 0x1D)
F = "C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_paper_report/"
prs = Presentation(); prs.slide_width = I(13.333); prs.slide_height = I(7.5)
BL = prs.slide_layouts[6]

def slide(dark=False):
    s = prs.slides.add_slide(BL)
    if dark:
        s.background.fill.solid(); s.background.fill.fore_color.rgb = DARK
    return s

def txt(s, x, y, w, h, lines, size=14, color=CHAR, bold=False, font="Calibri", align=PP_ALIGN.LEFT):
    tb = s.shapes.add_textbox(I(x), I(y), I(w), I(h)); tf = tb.text_frame
    tf.word_wrap = True
    for i, (t, kw) in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = t; p.alignment = kw.get("align", align)
        p.font.size = Pt(kw.get("size", size)); p.font.bold = kw.get("bold", bold)
        p.font.color.rgb = kw.get("color", color); p.font.name = kw.get("font", font)
        p.space_after = Pt(kw.get("after", 4))
    return tb

def title(s, t, sub=None, color=CHAR):
    txt(s, 0.55, 0.28, 12.2, 0.9, [(t, dict(size=30, bold=True, font="Cambria"))], color=color)
    if sub:
        txt(s, 0.55, 0.95, 12.2, 0.5, [(sub, dict(size=13, color=GREY))])

def pic(s, path, x, y, w):
    s.shapes.add_picture(F + path, I(x), I(y), width=I(w))

def table(s, x, y, w, rows, col_w=None, size=11, hdr=BLUE):
    nr, nc = len(rows), len(rows[0])
    shp = s.shapes.add_table(nr, nc, I(x), I(y), I(w), I(0.32 * nr)).table
    for j in range(nc):
        if col_w: shp.columns[j].width = I(col_w[j])
    for i, r in enumerate(rows):
        for j, v in enumerate(r):
            c = shp.cell(i, j); c.text = str(v)
            p = c.text_frame.paragraphs[0]; p.font.size = Pt(size)
            p.font.bold = (i == 0); p.font.color.rgb = WHITE if i == 0 else CHAR
            c.fill.solid(); c.fill.fore_color.rgb = hdr if i == 0 else (RGBColor(0xF2, 0xF5, 0xF7) if i % 2 else WHITE)

# 1 — titre (sombre)
s = slide(dark=True)
txt(s, 0.8, 2.3, 11.7, 1.6, [("Offset electronique du CL61", dict(size=44, bold=True, font="Cambria", color=WHITE)),
    ("et calibration Rayleigh : diagnostic, modele physique, correction, implications reseau", dict(size=21, color=RGBColor(0xCA,0xDC,0xFC)))])
txt(s, 0.8, 5.6, 11.7, 1.2, [("E-PROFILE / MeteoSwiss — Payerne, Aoste, Camborne, Lindenberg, Uccle", dict(size=14, color=RGBColor(0x9F,0xB4,0xC2))),
    ("Analyse 2026-07 · mesures hood mai–juin 2026 · rapport : cl61_rayleigh_investigation.md", dict(size=12, color=RGBColor(0x9F,0xB4,0xC2)))])

# 2 — probleme
s = slide(); title(s, "Le probleme : deux calibrations, deux constantes", "Payerne CL61 (ident C), mars–juin 2026")
txt(s, 0.55, 1.6, 4.4, 3.2, [("−12 %", dict(size=64, bold=True, color=RED)),
    ("ecart entre C_L(Rayleigh) = 1.25\net C_L(nuage) = 1.43", dict(size=15)),
    ("La methode nuage est ancree : accord −0.6 % avec le CHM15k co-localise (calibre independamment).", dict(size=13, color=GREY))])
table(s, 5.3, 1.7, 7.4, [
    ["Hypothese testee", "Verdict"],
    ["Correction vapeur d'eau absente (cal. / validation)", "Non — appliquee, essentielle (T²wv≈0.78 → +28 % sur C_L)"],
    ["Longueur d'onde laser (910.74 vs 910.55 nm)", "Non — +0.5 % ; enveloppe spectrale ≤ 2.4 %"],
    ["Contamination aerosol de la fenetre moleculaire", "Non — mauvais signe (gonflerait C_L)"],
    ["Rapport lidar (52 sr, transmission Klett)", "Non — 2σ ≈ 2–6 % (scan ±20 sr)"],
    ["CAMS 1° vs 0.4° / vs radiosondage", "Non — ≤0.3 % / −3.6 % (secondaire)"],
    ["Offset electronique (signal residuel negatif)", "CONFIRME — cause principale"]], col_w=[4.5, 2.9], size=11)

# 3 — mesures hood
s = slide(); title(s, "Mesure directe : telescope couvert (hood)", "3 fenetres operationnelles retrouvees dans l'archive L1 — aucune atmosphere : la moyenne EST l'offset")
pic(s, "fig_cl61_dark_windows.png", 0.4, 1.55, 8.6)
table(s, 9.25, 1.8, 3.7, [["z", "b_dark [Mm⁻¹sr⁻¹]"], ["3 km", "−0.008"], ["5 km", "−0.017"], ["8 km", "−0.027"]], col_w=[1.2, 2.5], size=12)
txt(s, 9.25, 3.6, 3.7, 3.2, [("Signature :", dict(bold=True, size=13)),
    ("• negatif, croissant avec l'altitude, remontant vers zero (~11 km)", dict(size=12)),
    ("• en signal non corrige de la distance β/z² : lobe positif proche puis undershoot negatif — reponse d'un filtre passe-haut (couplage AC), pas un offset constant amplifie par z²", dict(size=12)),
    ("• histogrammes quasi-gaussiens, centre decale : un decalage de distribution, pas des valeurs aberrantes", dict(size=12))])

# 4 — impact moleculaire
s = slide(); title(s, "Impact sur la calibration moleculaire", "L'offset vaut 3–13 % du signal moleculaire que la methode Rayleigh ajuste")
pic(s, "fig_cl61_dark_vs_molecular.png", 0.4, 1.6, 8.3)
table(s, 8.95, 1.8, 4.0, [["fenetre a", "biais C_L"], ["2 km", "−3.1 %"], ["3 km", "−7.8 %"], ["4 km", "−10.6 %"], ["5 km", "−13.0 %"], ["6 km", "−44.7 %"]], col_w=[1.7, 2.3], size=12)
txt(s, 8.95, 4.1, 4.0, 2.6, [("L'ecart observe (−12 %) correspond a une fenetre effective a 4.5–5 km — exactement ou l'algorithme place ses fenetres (2.6–6 km).", dict(size=13)),
    ("L'offset est tres inferieur au bruit par profil → invisible profil par profil, il n'emerge que dans la moyenne nocturne.", dict(size=12, color=GREY))])

# 5 — preuve profils
s = slide(); title(s, "Preuve causale (1/2) : profils avant / apres correction", "Nuits claires de calibration — la correction aligne le profil sur le modele « constante nuage »")
pic(s, "fig_cl61_profiles_beforeafter.png", 0.75, 1.6, 11.8)
txt(s, 0.75, 6.75, 11.8, 0.6, [("Corrige (bleu) ↔ modele C_L(nuage)·βmol·T² (rouge) dans la fenetre 3–5 km ; l'original (gris) suit le modele Rayleigh biaise bas : l'offset est le coin entre les deux methodes.", dict(size=12, color=GREY))])

# 6 — recalibration
s = slide(); title(s, "Preuve causale (2/2) : recalibration complete", "Soustraction de b_dark(z) puis recalibration eprof_v2 (WV active) nuit par nuit")
pic(s, "fig_cl61_offset_correction.png", 0.4, 1.6, 8.3)
txt(s, 8.9, 1.9, 4.0, 4.8, [("−26.5 % → −8.3 %", dict(size=30, bold=True, color=BLUE)),
    ("ecart median a la constante nuage, nuits communes (rampe lineaire)", dict(size=13)),
    ("+4 / +22 / +23 % par nuit appariee", dict(size=14, bold=True)),
    ("9 nuits eligibles au lieu de 7 : l'offset faisait aussi echouer le critere de fenetre (|intercept|).", dict(size=12)),
    ("Le modele physique (couplage AC, diapo suivante) fait aussi bien : −7.2 %, dispersion nuit-a-nuit 0.1 %.", dict(size=12, color=GREY))])

# 7 — modele physique (NOUVEAU)
s = slide(); title(s, "Modele physique : reponse d'une chaine a couplage AC", "Lobe positif rapide − undershoot lent remontant vers zero = reponse impulsionnelle d'un filtre passe-haut")
pic(s, "fig_cl61_offset_physical_model.png", 0.35, 1.55, 8.5)
txt(s, 8.95, 1.55, 4.05, 2.0, [
    ("P(r) = A_p e^(−r/L_p) − A_u e^(−r/L_u) + b_∞", dict(size=13, bold=True, color=BLUE)),
    ("Le couplage AC bloque le DC (fond + courant d'obscurite, ×1000 le signal). Sa reponse a une impulsion : l'impulsion MOINS un undershoot d'aire egale (integrale nulle).", dict(size=11.5))])
table(s, 8.95, 3.7, 4.05, [["terme", "τ = 2L/c", "origine"],
    ["lobe +", "4.7 µs", "transitoire proche"],
    ["undershoot −", "30.5 µs", "recovery couplage AC (RC)"]], col_w=[1.35, 1.0, 1.7], size=10)
txt(s, 8.95, 4.95, 4.05, 2.3, [
    ("• 3 regimes : lobe + (0–1.7 km), undershoot − (min 2.8 km), remontee lente puis plat (>11 km)", dict(size=11)),
    ("• la rampe lineaire = limite petit-r du terme lent ; τ_u≈30 µs vs fenetre ~100 µs → aspect quasi-lineaire", dict(size=11)),
    ("• le SIGNE tranche : undershoot AC = negatif remontant (nos hood) ; l'afterpulse donnerait une traine POSITIVE → exclu", dict(size=11)),
    ("• RMSE ÷ 2 vs la rampe ; recalibration −7.2 %", dict(size=11, bold=True))])

# 8 — dependance temperature (NOUVEAU)
s = slide(); title(s, "Dependance en temperature : offset surtout STABLE", "3 fenetres hood, temperature interne 22–44 °C — la constante RC ne derive pas, l'amplitude si (modestement)")
pic(s, "fig_cl61_offset_vs_temperature.png", 0.35, 1.55, 8.6)
txt(s, 9.05, 1.7, 3.95, 5.4, [
    ("Constante RC τ_u stable (~37 µs)", dict(size=13, bold=True, color=BLUE)),
    ("→ la forme de la remontee ne depend pas de la temperature.", dict(size=11.5, after=8)),
    ("Amplitude A_u : +quand l'electronique refroidit", dict(size=13, bold=True, color=BLUE)),
    ("~−1.5 %/K — coherent avec le coefficient de temperature NEGATIF du gain APD (plus froid → plus de gain → undershoot plus profond).", dict(size=11.5, after=8)),
    ("In-session ×1.13, poole ×1.32", dict(size=13, bold=True, color=BLUE)),
    ("dans la seule session de 25.5 h (seule la temperature varie) : ×1.13 ; le ×1.32 poole inclut une derive ENTRE sessions.", dict(size=11.5, after=8)),
    ("Consequence : un seul b(z) mesure retire l'essentiel du biais ; re-caracterisation periodique > indexation temperature instantanee.", dict(size=12, bold=True))])

# 9 — reseau
s = slide(); title(s, "Coherence reseau : le signe de l'offset explique tout", "Proxy mensuel 9–13 km (attendu ≈ +0.03) · firmware 1.2.7 constant partout — aucun evenement de configuration")
pic(s, "fig_cl61_network_coherence.png", 0.4, 1.55, 7.6)
table(s, 8.2, 1.7, 4.75, [
    ["Station", "Offset", "ratio ray/nuage"],
    ["Payerne", "negatif stable (−0.015)", "0.87 ✓"],
    ["Camborne", "negatif (−0.03…−0.05)", "0.945 ✓"],
    ["Aoste", "change de signe (hiver − / ete +)", "0.89 + cycle ✓"],
    ["Uccle", "±0.02 autour de zero", "serie Rayleigh inutilisable ✓"],
    ["Lindenberg", "POSITIF croissant (→ +0.10)", "1.06 (ray HAUT) ✓"]], col_w=[1.3, 2.05, 1.4], size=10)
txt(s, 8.2, 4.15, 4.75, 2.6, [("Aoste : la convergence Rayleigh–nuage a partir d'avril 2026 = relaxation de l'offset hivernal (dependance a la temperature interne) + montee du levier vapeur d'eau (corr +0.93). Pas de changement materiel/firmware.", dict(size=12))])

# 10 — universalite multi-instrument (NOUVEAU)
s = slide(); title(s, "Le phenomene est universel : les 3 instruments co-localises",
                   "Payerne CHM15k / CL31 / CL61 sous capot (4 sessions) — offset en % du signal moleculaire propre a chaque instrument")
pic(s, "fig_hood_fractional_bias.png", 0.3, 1.6, 8.7)
txt(s, 9.15, 1.65, 3.95, 5.4, [
    ("CHM15k (la reference !)", dict(size=13, bold=True, color=RED)),
    ("−15…−25 % dans la fenetre 3-5 km — PLUS GRAND que le CL61. Verifie : offset retire → C_L +11.5 % (vrai biais de calibration).", dict(size=11)),
    ("Mais l'ancre tient", dict(size=13, bold=True, color=BLUE)),
    ("la constante NUAGE du CHM15k vient du signal fort proche (0.1-2.4 km), immunisee → aucune constante Rayleigh de ceilometre n'est une reference propre.", dict(size=11)),
    ("CL31 : non calibrable en Rayleigh", dict(size=13, bold=True, color=RGBColor(0x94, 0x67, 0xbd))),
    ("offset ≳ signal moleculaire a 3-5 km → explique le statut operationnel CL31/CL51.", dict(size=11)),
    ("Mecanisme different (CHM15k photon-counting = sur-soustraction fond/overlap ; CL61 = undershoot AC), mais effet net partage : offset negatif en signal faible.", dict(size=10.5, color=GREY))])

# 11 — litterature
s = slide(); title(s, "Ce que dit la litterature", "Le phenomene est documente — le modele physique et la quantification C_L sont notre apport")
for x, t, b in ((0.55, "Kotthaus et al. 2016 (AMT) — CL31", "APD couplee en AC avec « zero variable » : baseline negative dependante de la distance (« cosmetic shift »), non nulle jusqu'a ~5.5 km, dependante de la temperature ; correction hood OBLIGATOIRE avant calibration. Le cas d'ecole du couplage AC."),
                (4.85, "Le et al. 2026 (AMT, en revision) — CL61", "« Residual background components may still remain » apres la correction interne ; soustraction hood P(r,T) indexee en temperature, repetee tous les ~6 mois ; une unite biaisee des 5 km — dans la fenetre Rayleigh."),
                (9.15, "Hopkin 2019 · Looschelders 2025 · electronique nucleaire", "« Doivent etre corriges pour les signaux faibles… negligeables pour les nuages » → la divergence exacte de nos 2 methodes. Le baseline shift des chaines AC a fort taux est un classique (Knoll) : remede = baseline restorer / compensation pole-zero.")):
    txt(s, x, 1.7, 3.9, 4.8, [(t, dict(size=14, bold=True, color=BLUE)), (b, dict(size=12))])
txt(s, 0.55, 6.5, 12.2, 0.7, [("Nouveaute : (i) quantifier le biais de C_L induit par l'offset dans un ajustement moleculaire 2–6 km ; (ii) un MODELE physique de l'offset (reponse a couplage AC, τ_p / τ_u) — extrapolable en distance, indexable en temperature.", dict(size=13, bold=True))])

# 11 — solutions
s = slide(); title(s, "Solutions", "Materielles (constructeur) et operationnelles (reseau)")
txt(s, 0.55, 1.6, 5.9, 5.4, [("Materiel / firmware (Vaisala)", dict(size=16, bold=True, color=BLUE)),
    ("• Restauration de ligne de base / compensation pole-zero (baseline restorer) dans l'etage analogique — annule l'undershoot A LA SOURCE (remede standard)", dict(size=13)),
    ("• Allonger τ (RC ≫ fenetre) ou couplage DC + ADC grande dynamique → undershoot negligeable/constant, trivial a soustraire", dict(size=13)),
    ("• Gating / masquage du flash de champ proche (moins de charge du condensateur)", dict(size=13)),
    ("• Echantillonnage pre-trigger de la ligne de base a chaque tir", dict(size=13)),
    ("• Stabilisation thermique (TEC) → (A_u, τ_u) constants, corrigeables une fois pour toutes", dict(size=13)),
    ("• Publier le signal BRUT (ou la baseline estimee) pour permettre la correction utilisateur", dict(size=13))])
txt(s, 6.85, 1.6, 6.0, 5.4, [("Operationnel (E-PROFILE, des maintenant)", dict(size=16, bold=True, color=BLUE)),
    ("• Constante NUAGE comme reference CL61 (deja le cas) ; Rayleigh natif = biais −10…−15 %", dict(size=13)),
    ("• Soustraire le modele physique b(z)=[A_p e^−r/Lp − A_u e^−r/Lu]·z² (extrapolable, indexable en T)", dict(size=13)),
    ("• Tests hood tous les ~3 mois par unite (la derive ENTRE sessions domine le terme temperature)", dict(size=13)),
    ("• Methode « residu nuits claires » : residu du profil vs modele moleculaire complet, 3–6 km (validee : −0.021 ciel clair vs −0.015 hood)", dict(size=13)),
    ("• Suivi mensuel du proxy 9–13 km sur le dashboard (derive d'offset par unite)", dict(size=13)),
    ("• Escalade Vaisala avec le dossier complet", dict(size=13))])

# 12 — conclusions (sombre)
s = slide(dark=True)
txt(s, 0.8, 0.6, 11.7, 0.9, [("Conclusions", dict(size=34, bold=True, font="Cambria", color=WHITE))])
txt(s, 0.8, 1.6, 11.7, 5.5, [
    ("1.  L'ecart Rayleigh–nuage du CL61 (−12 %) est cause aux ⅔ par un offset electronique negatif, mesure directement au telescope couvert (−0.008…−0.027 Mm⁻¹sr⁻¹).", dict(size=15, color=WHITE, after=9)),
    ("2.  La preuve est causale : offset soustrait → ecart −26.5 % → −8.3 % (rampe) / −7.2 % (modele physique), dispersion nuit-a-nuit 0.1 %, et plus de nuits calibrables.", dict(size=15, color=WHITE, after=9)),
    ("3.  Le phenomene est general : coherent sur tout le reseau CL61 (Lindenberg, Aoste, Uccle) ET universel entre types d'instruments — le CHM15k (ancre) porte lui-meme −15…−25 % (biais verifie +11.5 %), mais sa constante NUAGE reste immunisee.", dict(size=15, color=WHITE, after=9)),
    ("4.  Mecanisme : reponse d'une chaine a couplage AC (undershoot passe-haut), τ_p≈4.7 µs / τ_u≈30.5 µs ; l'afterpulse est exclu par le signe. Dependance temperature modeste, dans l'amplitude (gain APD) ; RC stable.", dict(size=15, color=WHITE, after=9)),
    ("5.  A faire : re-caracterisation hood ~trimestrielle, correction par modele physique, methode residu en routine, dossier Vaisala (baseline restorer).", dict(size=15, color=RGBColor(0xCA,0xDC,0xFC), after=9))])

import os
out = r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/presentation/CL61_offset_electronique_FR.pptx"
os.makedirs(os.path.dirname(out), exist_ok=True)
prs.save(out)
print("saved", out, "-", len(prs.slides._sldIdLst), "slides")
