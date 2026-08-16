# -*- coding: utf-8 -*-
"""v3 control matrix: which instrument offers which calibration METHOD, which VARIANTS exist for
each method, and where each variant's per-night constants come from.

This is the file to edit when a new run lands.  `source_for()` returns either a source spec the
builder can read, or a STRING explaining why that (instrument, method, variant) does not exist --
the page greys the option out and shows the string as its tooltip.  Nothing is ever silently
carried over from another variant.

Method availability is a property of the instrument type:
  CHM15k   Rayleigh only  -- it saturates in liquid cloud, so a cloud calibration is meaningless
                             (doc/reports, "CHM15k: no cloud calibration")
  CL31     cloud only     -- no Rayleigh night survives its molecular-window SNR at these stations
  CL61     both           -- two independent retrievals of the same constant
"""
from __future__ import annotations
from pathlib import Path

DATA = Path(r"C:/DATA/Projects/202606_E-PROFILE_calibration")
NETWORK_V22 = DATA / "calout_v22_04"          # network eprof_v2.2 run, 0.4 deg CAMS
DARK_EST = DATA / "dark_test"                 # same runner + ESTIMATED clear-night dark
DARK_MEAS = DATA / "diag_v22_dark"            # Payerne: MEASURED covered-telescope dark b(z)
NOWV = DATA / "diag_v22_nowv"                 # same runner, water-vapour correction switched off

# --------------------------------------------------------------------------- CL61 WV spectrum
# The water-vapour correction integrates the LASER LINE over the 910 nm absorption band, so its
# amplitude is set by (lambda0, FWHM).  The operational model uses 910.74 nm / FWHM 1.0 nm.  The
# authors of the Le & O'Connor preprint report, from a Vaisala personal communication (response to
# referee RC1), the manufacturer's true spectrum: lambda0 = 910.55 nm, sigma = 0.08 nm, i.e.
# FWHM = 2*sqrt(2*ln2)*sigma = 0.188 nm.  That line is narrow enough to sit BETWEEN the strong
# absorption features, so the operational model over-corrects by roughly a factor 9 -- corroborated
# three independent ways (the authors' own HITRAN line-by-line gives the CL51 ~10x more water-vapour
# attenuation than the CL61; our own LUT gives 9.3x; our empirical cloud-scene measurement says the
# raw CL61 signal carries only ~8 % of the modelled absorption).
#
# The 910.55 / sigma 0.08 run is therefore the PHYSICALLY CORRECT hypothesis and is flagged as the
# page's reference; the others (910.74 operational, 910.55/FWHM 1.0, 910.55/FWHM 0.1, no WV at all)
# stay as the comparison ladder that brackets it.
LAM_RUNS = {                                  # variant suffix -> runner output tree
    "l55s008": DATA / "diag_v22_l55s008",     # 910.55 nm, sigma 0.08 -> FWHM 0.188  (constructeur)
    "l55w10":  DATA / "diag_v22_l55w10",      # 910.55 nm, FWHM 1.0   (line position only)
    "l55w01":  DATA / "diag_v22_l55w01",      # 910.55 nm, FWHM 0.1   (narrower than the truth)
}
REFERENCE_SUFFIX = "l55s008"
# Only the CL61 streams were re-run: the spectrum in question is the CL61's.
LAM_TYPES = ("CL61",)

# Ordered variant lists per method (the page renders them in this order).
VARIANTS = {
    "rayleigh": ["v2.0", "v2.2", "v2.2dark", "v2.2sansWV",
                 "v2.2_l55s008", "v2.2_l55w10", "v2.2_l55w01"],
    "cloud":    ["cloudWV", "cloudNoWV",
                 "cloud_l55s008", "cloud_l55w10", "cloud_l55w01"],
}
VARIANT_LABEL = {
    "v2.0":         "v2.0 (opérationnel)",
    "v2.2":         "v2.2 — WV λ910,74 / FWHM 1,0 (modèle actuel)",
    "v2.2dark":     "v2.2 + fond soustrait (dark)",
    "v2.2sansWV":   "v2.2 — sans correction WV",
    "v2.2_l55s008": "v2.2 — WV λ910,55 σ0,08 (spectre constructeur) ★",
    "v2.2_l55w10":  "v2.2 — WV λ910,55 / FWHM 1,0",
    "v2.2_l55w01":  "v2.2 — WV λ910,55 / FWHM 0,1",
    "cloudWV":      "nuage — WV λ910,74 / FWHM 1,0 (modèle actuel)",
    "cloudNoWV":    "nuage — sans correction WV",
    "cloud_l55s008": "nuage — WV λ910,55 σ0,08 (spectre constructeur) ★",
    "cloud_l55w10":  "nuage — WV λ910,55 / FWHM 1,0",
    "cloud_l55w01":  "nuage — WV λ910,55 / FWHM 0,1",
}
# A 1064 nm instrument has no water-vapour term at all, so naming its variants after a laser line
# would be nonsense — the page swaps in these labels for the CHM15k rows.
# Per-type label overrides.  The generic labels above carry the CL61 spectrum (910.74/1.0)
# because the lambda LADDER is a CL61 story; every other type must display ITS OWN water-vapour
# model -- CL31 909.7/6.0, CL51 910.0/3.4 (Qmini campaign / Wiegner) -- and the 1064 nm CHM15k
# no spectrum at all (operator rule: never name a 1064 nm variant after a laser line).
VARIANT_LABEL_BY_TYPE = {
    "CHM15k": {"v2.2": "v2.2"},
    "CL31": {"v2.2":    "v2.2 — WV λ909,7 / FWHM 6,0 (modèle actuel)",
             "cloudWV": "nuage — WV λ909,7 / FWHM 6,0 (modèle actuel)"},
    "CL51": {"v2.2":    "v2.2 — WV λ910,0 / FWHM 3,4 (modèle actuel)",
             "cloudWV": "nuage — WV λ910,0 / FWHM 3,4 (modèle actuel)"},
}
VARIANT_SHORT_BY_TYPE = {
    "CHM15k": {"v2.2": "v2.2"},
    "CL31": {"v2.2": "v2.2 (λ909,7)", "cloudWV": "nuage (λ909,7)"},
    "CL51": {"v2.2": "v2.2 (λ910,0)", "cloudWV": "nuage (λ910,0)"},
}
# retro-compat (anciens noms utilises par le rendu)
VARIANT_LABEL_1064 = VARIANT_LABEL_BY_TYPE["CHM15k"]
VARIANT_SHORT_1064 = VARIANT_SHORT_BY_TYPE["CHM15k"]


def variants_for(itype, method):
    """The variants that are PHYSICALLY MEANINGFUL for this instrument type and method.

    A dropdown must only offer entries that mean something for the instrument in front of it;
    "greyed out with a tooltip" is reserved for the different situation of a meaningful variant
    whose run has not been produced yet.  Water vapour does not touch a 1064 nm unit at all, and
    the laser-line ladder is a statement about the CL61's emission spectrum specifically.
    """
    full = VARIANTS.get(method, [])
    if itype == "CHM15k":                       # 1064 nm: no WV term exists in its constant
        return [v for v in full if v not in ("v2.2sansWV",) and v.split("_")[-1] not in LAM_RUNS]
    if itype not in LAM_TYPES:                  # 910 nm but not the re-run type (CL31, CL51)
        return [v for v in full if v.split("_")[-1] not in LAM_RUNS]
    return list(full)


# Profile-side measured-dark checkbox: the constants-side twin to switch to when it is ticked, so
# a dark-corrected profile is never divided by a dark-free constant (mixing two signal definitions
# measurably WORSENS the comparison).  Anything not listed has no dark-corrected run.
DARK_TWIN = {"v2.0": "v2.2dark", "v2.2": "v2.2dark", "v2.2dark": "v2.2dark"}
DARK_NO_TWIN_WHY = ("aucun run d'étalonnage dark-corrigé pour cette variante : seul le PROFIL est "
                    "corrigé ici, la constante reste dark-free. État hybride — à lire comme un "
                    "diagnostic, pas comme une comparaison cohérente. (Pour la méthode nuage "
                    "l'effet est de toute façon négligeable : l'intégration 100–2400 m se fait là "
                    "où dark/signal ≈ 1e-3.)")
DARK_NO_MEAS_WHY = ("pas de mesure capot pour cette unité (campagne télescope couvert à Payerne "
                    "uniquement)")
DARK_HINT = "cohérence : profils et constantes corrigés ensemble"
# Short forms for legends and chart titles (the long labels do not fit a plotly legend).
VARIANT_SHORT = {
    "v2.0": "v2.0", "v2.2": "v2.2 (λ910,74)", "v2.2dark": "v2.2+dark",
    "v2.2sansWV": "v2.2 sans WV", "v2.2_l55s008": "v2.2 λ910,55 σ0,08 ★",
    "v2.2_l55w10": "v2.2 λ910,55 F1,0", "v2.2_l55w01": "v2.2 λ910,55 F0,1",
    "cloudWV": "nuage (λ910,74)", "cloudNoWV": "nuage sans WV",
    "cloud_l55s008": "nuage λ910,55 σ0,08 ★", "cloud_l55w10": "nuage λ910,55 F1,0",
    "cloud_l55w01": "nuage λ910,55 F0,1",
}
# Variants that carry the manufacturer spectrum — the page badges them as the reference hypothesis.
REFERENCE_VARIANTS = ("v2.2_l55s008", "cloud_l55s008")
REFERENCE_NOTE = ("Spectre d'émission constructeur (Vaisala, comm. pers. citée dans la réponse au "
                  "relecteur RC1 du preprint Le & O'Connor) : lambda0 = 910,55 nm, sigma = 0,08 nm, "
                  "soit FWHM = 0,188 nm. Le modèle opérationnel (910,74 / FWHM 1,0) sur-corrige "
                  "l'absorption d'environ un facteur 9.")
# Per-PROFILE cloud-calibration dumps (rayleigh_availability/cloud_profile_dump.py, one output
# folder per water-vapour configuration).  They feed the Hopkin-style CBH heatmap, which is the
# per-profile view of the same constant the page divides by -- so the panel follows the operator's
# CLOUD variant choice instead of being frozen.
CLOUD_DUMP_TAG = {
    "cloudWV": "nominal", "cloudNoWV": "nowv", "cloud_l55s008": "l55s008",
    "cloud_l55w10": "l55w10", "cloud_l55w01": "l55w01",
}


def dump_npz(variant, wmo, ident):
    tag = CLOUD_DUMP_TAG.get(variant)
    if tag is None:
        return None
    return DATA / f"cloud_profile_dump_{tag}" / f"{wmo}_{ident}_profiles.npz"


METHOD_LABEL = {"rayleigh": "Rayleigh", "cloud": "Nuages"}
METHOD_BY_TYPE = {"CHM15k": ["rayleigh"], "CL31": ["cloud"], "CL61": ["rayleigh", "cloud"]}
METHOD_WHY = {
    ("CHM15k", "cloud"): "le CHM15k sature dans les nuages liquides — étalonnage nuage impossible",
    ("CL31", "rayleigh"): "aucune nuit Rayleigh exploitable pour le CL31 (SNR insuffisant dans la "
                          "fenêtre moléculaire) — étalonnage nuage uniquement",
}

# Instruments at 1064 nm have no water-vapour correction to remove.
_NO_WV_TO_REMOVE = "instrument à 1064 nm : aucune correction vapeur d'eau n'entre dans sa " \
                   "constante — les variantes « sans WV » et « spectre laser » sont sans objet"
_NO_NC = "pas de NetCDF de calibration v2.0 publié pour cette station sur ce poste"


def _csv(root, wmo, ident):
    return root / f"{wmo}_{ident}" / f"{wmo}_{ident}_cal.csv"


def source_for(site_key, ident, itype, method, variant, wmo):
    """-> ("nc"|"json"|"csv", arg) or a French string explaining the absence."""
    is1064 = itype == "CHM15k"

    # --- the CL61 water-vapour-spectrum ladder --------------------------------------------------
    suffix = variant.split("_")[-1] if "_" in variant else None
    if suffix in LAM_RUNS:
        if is1064:
            return _NO_WV_TO_REMOVE
        if itype not in LAM_TYPES:
            return (f"le rejeu du spectre laser ne couvre que les {', '.join(LAM_TYPES)} — "
                    f"la raie du {itype} est différente et n'a pas été rejouée")
        root = LAM_RUNS[suffix]
        if not _csv(root, wmo, ident).exists():
            return (f"rejeu « {suffix} » pas encore écrit pour ce flux "
                    f"({root.name}/{wmo}_{ident})")
        return ("csv", _csv(root, wmo, ident))

    if method == "cloud":
        if site_key == "payerne":
            # "avec WV" = the published v2.0/v2.2 cloud constants (the cloud retrieval is identical
            # between those two: v2.2 only moves the Rayleigh gates), "sans WV" = the diag_v22_nowv
            # rerun.  The two therefore come from different runs -- stated on the page.
            return ("nc", None) if variant == "cloudWV" else ("csv", _csv(NOWV, wmo, ident))
        if variant == "cloudWV":
            return ("csv", _csv(NETWORK_V22, wmo, ident))
        return ("csv", _csv(NOWV, wmo, ident))

    # --- Rayleigh -------------------------------------------------------------------------------
    if variant == "v2.0":
        if site_key != "payerne":
            return _NO_NC
        return ("nc", None)
    if variant == "v2.2":
        if site_key == "payerne":
            return ("json", "A" if ident == "A" else "Cr")
        return ("csv", _csv(NETWORK_V22, wmo, ident))
    if variant == "v2.2dark":
        if site_key == "payerne":
            return ("csv", _csv(DARK_MEAS, wmo, ident))
        sub = {"amsterdam": "amst_dark", "lindenberg": "lind_dark"}.get(site_key)
        if sub is None:
            return "aucun run « dark » pour cette station"
        return ("csv", _csv(DARK_EST / sub, wmo, ident))
    if variant == "v2.2sansWV":
        if is1064:
            return _NO_WV_TO_REMOVE
        return ("csv", _csv(NOWV, wmo, ident))
    return f"variante « {variant} » inconnue"


_SPECTRUM_WARN = """<b>★ Spectre d'émission du CL61 — hypothèse de référence.</b> Les auteurs
    du preprint Le &amp; O'Connor rapportent, d'une communication personnelle Vaisala (réponse au
    relecteur RC1), la raie réelle du CL61 : <b>λ₀ = 910,55 nm, σ = 0,08 nm</b>, soit
    <b>FWHM = 0,188 nm</b>. Le modèle opérationnel (910,74 nm / FWHM 1,0) <b>sur-corrige
    l'absorption d'un facteur ≈ 9</b> — corroboré trois fois indépendamment : le HITRAN
    raie-par-raie des auteurs donne au CL51 ≈ 10× l'atténuation WV du CL61, notre propre LUT donne
    9,3×, et notre mesure empirique sur scènes nuageuses dit que le signal brut du CL61 ne porte
    que ~8 % de l'absorption modélisée. Les variantes marquées <b>★</b> (constantes) et le mode
    <b>« λ910,55 σ0,08 »</b> de la correction WV de comparaison portent ce spectre : les choisir des
    <b>deux côtés</b> est l'état physiquement cohérent. Les autres λ ne sont là que pour encadrer."""

_PAY_WARN = [
    """<b>CL31 : bloc optique remplacé le 2026-07-07 ~13:00.</b> Sa constante passe de ~3.1e7
    (jan–mai) / 4.1e7 (juin) à 9.1e7 (juil–août). La fenêtre couvre la bascule <b>volontairement</b> :
    choisir <b>juin 2026</b> pour l'état d'avant, <b>juillet</b>/<b>août</b> pour l'état d'après. Le
    lissage est le filtre opérationnel (<code>monitoring.kalman</code>), qui absorbe la marche en
    quelques jours ; le panneau L2, lui, porte le défaut fixe 1e8 sur toute la période.""",
    """<b>La référence CHM15k est peu contrainte sur cette fenêtre</b> : une seule nuit Rayleigh
    acceptée entre le 8 juillet et le 13 août, après deux mois de trou. Le filtre bouge donc à peine
    et juillet/août tournent sur la valeur de fin bloquée — une part du résidu CL61 est dans la
    référence, pas dans le CL61.""",
    """<b>Nuages « avec WV » vs « sans WV » :</b> les constantes « avec WV » sont celles des
    NetCDF v2.0/v2.2 publiés (le calcul nuage est identique entre ces deux versions), les
    constantes « sans WV » viennent du rerun <code>diag_v22_nowv</code>. Le test <b>cohérent</b> =
    variante « sans WV » avec la case « <b>WV comparaison</b> » <b>décochée</b> ; la vue historique
    = « avec WV » + case <b>cochée</b>. Les états croisés restent des diagnostics.""",
    """<b>v2.2 + dark (Payerne)</b> : le fond électronique <b>mesuré</b> (télescope couvert) est
    soustrait du profil L1 <i>et</i> la constante vient du run dark-corrigé — les deux côtés ou
    aucun. Les stations d'Amsterdam et de Lindenberg utilisent un dark <b>estimé</b> depuis les
    nuits claires, qui ne touche que la constante (variante de démonstration).""",
    _SPECTRUM_WARN,
]

_AMS_WARN = [
    """<b>Artefact instrumental connu :</b> la fonction de recouvrement appliquée à l'unité B est
    ~25 % trop basse entre 500 et 1000 m (elle revient à 1 vers 1600 m) : l'excès de B en champ
    proche est une erreur de recouvrement, pas d'étalonnage.""",
    """Les quatre unités sont à 1064 nm : les corrections vapeur d'eau et longueur d'onde sont des
    <b>no-ops exacts</b> et les commandes correspondantes sont désactivées pour cette station.""",
]

_LIN_WARN = [
    """<b>Fenêtre et constantes bornées :</b> le L1 journalier local du CL61 s'arrête le
    30 juin 2026, la fenêtre appariée est donc juin uniquement. Son enregistrement v2.2 (run réseau)
    s'arrête le 15 juin ; la série de Kalman est maintenue à sa dernière valeur ensuite.""",
    _SPECTRUM_WARN,
]

SITE_V3 = {
    "payerne": dict(
        instruments=[
            dict(ident="A", itype="CHM15k", label="CHM15k (A)", color="#1f77b4"),
            dict(ident="B", itype="CL31", label="CL31 (B)", color="#ff7f0e"),
            dict(ident="C", itype="CL61", label="CL61 (C)", color="#2ca02c"),
        ],
        default={"A": ("rayleigh", "v2.2"), "B": ("cloud", "cloudWV"),
                 "C": ("cloud", "cloudWV")},
        iref=0,
        # Each state = (constants variant, profile-side WV mode).  A COHERENT state uses the same
        # laser-line hypothesis on both sides; a crossed one mixes two definitions of the signal
        # and is drawn pale.  A flat residual-vs-PWV slope is the signature of the right spectrum.
        pwv=dict(ident="C", ref="A", method="rayleigh",
                 states=[["v2.2", "nom", "λ910,74 FWHM 1,0 des deux côtés — modèle actuel",
                          "#d62728", "circle"],
                         ["v2.2_l55s008", "ctor", "λ910,55 σ0,08 des deux côtés — spectre "
                          "constructeur ★", "#2f9e44", "diamond"],
                         ["v2.2sansWV", "none", "sans WV des deux côtés", "#1f77b4", "square"],
                         ["v2.2", "none", "cal. avec WV × comparaison sans WV — croisé",
                          "#f1a8a8", "circle-open"],
                         ["v2.2sansWV", "nom", "cal. sans WV × comparaison avec WV — croisé",
                          "#a8d8b0", "square-open"]]),
        title="Payerne — trois céilomètres co-localisés : L1 + étalonnage v2 face au produit L2",
        subtitle="""Rétrodiffusion atténuée sur les mêmes heures, la même grille d'altitude et les
    mêmes corrections des deux côtés. <span class="muted">Gauche = L1 <code>rcs_0</code> divisé par
    la constante de Kalman choisie ci-dessus, instrument par instrument. Droite = L2
    <code>attenuated_backscatter_0</code> tel que distribué.</span>""",
        warnings=_PAY_WARN,
    ),
    "amsterdam": dict(
        instruments=[
            dict(ident="A", itype="CHM15k", label="CHM15k (A)", color="#1f77b4"),
            dict(ident="B", itype="CHM15k", label="CHM15k (B)", color="#ff7f0e"),
            dict(ident="C", itype="CHM15k", label="CHM15k (C)", color="#2ca02c"),
            dict(ident="D", itype="CHM15k", label="CHM15k (D)", color="#d62728"),
        ],
        default={i: ("rayleigh", "v2.2") for i in "ABCD"},
        iref=0,
        title="Amsterdam — quatre CHM15k co-localisés : L1 + étalonnage v2.2",
        subtitle="""Rétrodiffusion atténuée sur les mêmes heures et la même grille d'altitude pour
    les quatre CHM15k de Schiphol. <span class="muted">Aucune archive L2 locale : pas de panneau
    « produit distribué ».</span>""",
        warnings=_AMS_WARN,
    ),
    "lindenberg": dict(
        instruments=[
            dict(ident="0", itype="CHM15k", label="CHM15k (0)", color="#1f77b4"),
            dict(ident="C", itype="CL61", label="CL61 (C)", color="#2ca02c"),
        ],
        default={"0": ("rayleigh", "v2.2"), "C": ("cloud", "cloudWV")},
        iref=0,
        title="Lindenberg — CHM15k face au CL61 : L1 + étalonnage v2.2",
        subtitle="""Rétrodiffusion atténuée sur les mêmes heures et la même grille d'altitude pour
    le CHM15k (1064 nm) et le CL61 (910 nm) co-localisés. <span class="muted">Aucune archive L2
    locale.</span>""",
        warnings=_LIN_WARN,
    ),
}
