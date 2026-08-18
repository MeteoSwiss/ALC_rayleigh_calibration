# -*- coding: utf-8 -*-
"""Site table for the L1-vs-L2 inter-comparison dashboard.

One SiteSpec dict per station; build_l1_l2_dashboard / l1_l2_calib / render_l1_l2_dashboard read
everything station-specific from here, so adding a station (or a calibration variant to one) is an
entry in this file, not an edit of the pipeline.

Fields:
  wmo/lat/lon/alt      station identity (lat/lon/alt from validation.paper.calib_benchmark)
  channels             ident/itype/calib/label/color; `chan` disambiguates two calibrations of the
                       same physical unit (the CL61 Rayleigh-vs-cloud trick, Payerne C/Cr)
  iref                 index into channels of the reference everything is compared to
  sources              ("L1","L2") or ("L1",) — the render stage hides every L2 element when absent
  start/end            archive read window (what goes into the stream cache)
  win_start/win_end    comparison window (clamped to what all streams actually cover)
  variants             ORDERED calibration-variant names; variants[0] is the base/canonical one
  calib_dirs           variant -> folder of the smoothed <wmo>_<ident>_<calib>_L1.csv series
  calib_builder        "payerne" = the bespoke l1_l2_calib.main() (NetCDFs + study outputs);
                       "network" = l1_l2_calib.main_site() from the network-runner CSVs (run_dirs)
  run_dirs             ("network" builder only) variant -> runner output root holding
                       <wmo>_<ident>/<wmo>_<ident>_cal.csv — a future "v2.2dark" run for these
                       sites is ONE more entry in variants + run_dirs + calib_dirs
  dark_variants        variants that subtract the measured dark baseline from the L1 profiles
  dark_npz             the measured b(z) npz those variants read (Payerne only has one)
  collapse_combos      True = every channel is 1064 nm, WV/wavelength are no-ops -> compute one
                       combo per variant and hide those segmented controls in the page
  title/subtitle/warnings   HTML fragments for the page header (warnings = list of .warn boxes)
"""
from __future__ import annotations
from pathlib import Path

DATA_ROOT = Path(r"C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard")
# The network-wide eprof_v2.2 run (0.4 deg CAMS): per-stream <wmo>_<ident>/<wmo>_<ident>_cal.csv.
NETWORK_V22 = Path(r"C:/DATA/Projects/202606_E-PROFILE_calibration/calout_v22_04")
# Same runner, with the ESTIMATED (clear-night) dark subtracted inside the calibration.
# Demonstration variant only — see the per-site note.
DARK_EST = Path(r"C:/DATA/Projects/202606_E-PROFILE_calibration/dark_test")
# Payerne rerun WITHOUT the water-vapour correction inside the calibration (eprof_v2.2 settings
# otherwise identical; 20250101-20260813). Motivated by the measured CL61 emission line
# (910.74 ± 0.10 nm, FWHM < 1.5 nm) sitting between the strong WV absorption features: the
# semi-analytic replay shows the CL61 RAW cloud signal carries only ~8 % of the modelled WV path
# absorption (CL51 45 %, CL31 80 %), so the pipeline's WV correction is suspected ~12x too strong
# for the CL61 — this variant lets the page pair no-WV constants with the WV-free profile view.
NOWV_RUN = Path(r"C:/DATA/Projects/202606_E-PROFILE_calibration/diag_v22_nowv")

# Note appended to the sites that carry the demonstration variant.
_DARK_EST_NOTE = """<b>v2.2 dark-est</b> = soustraction du dark ESTIMÉ des nuits claires
    (θ ≈ 0 à cette station — variante de démonstration, non recommandée opérationnellement)."""

_PAYERNE_SUBTITLE = """Attenuated backscatter from the same hours, the same altitude grid and the
    same corrections on both sides. <span class="muted">Left = L1 <code>rcs_0</code> divided by the
    daily Kalman C<sub>L</sub> from the latest v2.0 calibration NetCDFs. Right = L2
    <code>attenuated_backscatter_0</code> exactly as distributed (CHM15k v1.0 Rayleigh, CL31 the
    1e8 default, CL61 Vaisala's internal constant). Read from the daily
    <code>E-PROFILE_L1_2026</code> / <code>E-PROFILE_L2_2026</code> archives, falling back to the
    5-minute granules on days with no concatenated file.</span>"""

_PAYERNE_WARNINGS = [
    """<b>CL31 optical block replaced 2026-07-07 ~13:00</b> (after that day's dark
    measurement). Its per-night constant steps from ~3.1e7 (Jan–May) / 4.1e7 (Jun) to 9.1e7 (Jul–Aug).
    <b>This window spans the swap on purpose</b> — set the period to <b>Jun 2026</b> for the state
    before it and <b>Jul</b>/<b>Aug</b> for after. The smoothing is the operational dashboard's own
    filter (<code>monitoring.kalman</code>: median-normalised, 4 % day-to-day drift, 15 %/yr over
    gaps, rolling-IQR outlier rejection), which absorbs the step within days — 4.11e7 on 1 Jul,
    8.46e7 by 20 Jul, 8.00e7 on 6 Aug, identical to the station page — so the L1 panel stays valid
    either side of it. The L2 panel does not: it carries the fixed 1e8 default throughout, which
    happens to sit near the post-swap truth and far from the pre-swap one.""",
    """<b>The CHM15k reference is under-constrained over this window.</b> Exactly
    <b>one</b> Rayleigh night is accepted between 8 Jul and 13 Aug (22 Jul, 8.28e11), after a
    two-month gap from 22 May (6.91e11), so the random-walk filter barely moves: it holds
    <b>6.84e11</b> — and the filter's rolling-IQR test rejects that night outright, so its series
    ends 22 May and July/August run on the clamped end value. Scale the reference to the 22 Jul night
    and CL61 moves from −14 % to about +3 %, so most of the residual sits in the reference rather
    than in CL61, whose own cloud calibration has 82 nights and tracks well. Constants come only from
    the v2.0 NetCDFs, all available years (43 CHM15k / 264 CL31 / 82 CL61 nights); the raw
    operational <code>_cal.csv</code> is deliberately not used, as its flag=0.5 nights include two
    wild CHM15k outliers (20–21 Jun, 3.4–4.1e12) that the published NetCDF drops.""",
    """<b>v2.2 sans WV(cal)</b> : constantes calculées <b>SANS</b> la correction de vapeur d'eau
    (raie CL61 mesurée à 910.74 nm, FWHM &lt; 1.5 nm — l'absorption WV est atténuée par
    conception spectrale, le signal brut n'en porte que ~8 % du modèle). Le test <b>cohérent</b> =
    cette variante avec la case «&nbsp;WV comparaison&nbsp;» <b>décochée</b> ; la vue historique =
    <b>v2.2</b> avec «&nbsp;WV comparaison&nbsp;» <b>cochée</b>. Les états croisés (sans-WV(cal) ×
    comparaison ON, ou v2.2 × comparaison OFF) sont des diagnostics seulement — ils mélangent deux
    définitions du signal.""",
]

SITES = {
    "payerne": dict(
        key="payerne", name="Payerne", wmo="0-20000-0-06610",
        lat=46.8137, lon=6.9425, alt=491.0,
        channels=[
            dict(ident="A", itype="CHM15k", calib="rayleigh", label="CHM15k (A)", color="#1f77b4"),
            dict(ident="B", itype="CL31",   calib="cloud",    label="CL31 (B)",   color="#ff7f0e"),
            dict(ident="C", itype="CL61",   calib="cloud",    label="CL61 (C)",   color="#2ca02c"),
            # The CL61 carries a Rayleigh calibration as well as the cloud one; it is a separate
            # retrieval of the same constant, and the only CL61 channel a molecular-gate change moves.
            dict(ident="C", itype="CL61",   calib="rayleigh", label="CL61 (C, Rayleigh)",
                 color="#17becf", chan="Cr"),
        ],
        iref=0,
        sources=("L1", "L2"),
        start="20260601", end="20260814",
        win_start="20260601", win_end="20260813",
        variants=["v2.0", "v2.2", "v2.2dark", "v2.2sansWV"],
        variant_labels={"v2.0": "v2.0 (operational)", "v2.2": "v2.2 (noise-aware gates)",
                        "v2.2dark": "v2.2 + dark (electronic baseline subtracted)",
                        "v2.2sansWV": "v2.2 sans WV(cal)"},
        variant_short={"v2.2dark": "v2.2+dark", "v2.2sansWV": "v2.2 sansWV"},
        calib_dirs={"v2.0": DATA_ROOT / "calib", "v2.2": DATA_ROOT / "calib_v22",
                    "v2.2dark": DATA_ROOT / "calib_v22dark",
                    "v2.2sansWV": DATA_ROOT / "calib_v22nowv"},
        calib_builder="payerne",
        dark_variants=("v2.2dark",),
        dark_npz=str(DATA_ROOT.parent / "rayleigh_availability" / "dark_profiles_payerne.npz"),
        collapse_combos=False,
        # Page-size containment: the single-Ångström wavelength mode is dropped from the COMBOS
        # (4 variants x 2 WV x 2 modes = 16 combos instead of 18). The other sites keep the default.
        wl_modes=["none", "molecular"],
        # The discriminating test: per-hour CL61-Rayleigh residual vs the CHM15k reference, against
        # the CAMS PWV of that hour, for the 4 coherent/crossed (variant x WV-comparison) states.
        pwv=dict(chan="Cr", cals=["v2.2", "v2.2sansWV"], wl="molecular"),
        title="Payerne — three co-located ceilometers: L1 + v2 calibration vs the L2 product",
        subtitle=_PAYERNE_SUBTITLE,
        warnings=_PAYERNE_WARNINGS,
        profiles_note="""Each curve stops where that instrument runs out of signal, not at a fixed
    height; the grid runs to 15 km (CHM15k 15.3, CL61 15.7, CL31 7.7 km native). Default view
    0–4 km — zoom out for the rest.""",
    ),

    "amsterdam": dict(
        key="amsterdam", name="Amsterdam", wmo="0-20000-0-06240",
        lat=52.317, lon=4.8037, alt=6.0,
        # Four co-located CHM15k units, all Rayleigh-calibrated. Same type everywhere -> the
        # default palette instead of the per-type colours (validation.paper.figures rule).
        channels=[
            dict(ident="A", itype="CHM15k", calib="rayleigh", label="CHM15k (A)", color="#1f77b4"),
            dict(ident="B", itype="CHM15k", calib="rayleigh", label="CHM15k (B)", color="#ff7f0e"),
            dict(ident="C", itype="CHM15k", calib="rayleigh", label="CHM15k (C)", color="#2ca02c"),
            dict(ident="D", itype="CHM15k", calib="rayleigh", label="CHM15k (D)", color="#d62728"),
        ],
        iref=0,
        sources=("L1",),                    # no local L2 archive for this station
        start="20260601", end="20260814",
        win_start="20260601", win_end="20260813",
        variants=["v2.2", "v2.2dark_est"],
        variant_labels={"v2.2": "v2.2 (network run)", "v2.2dark_est": "v2.2 dark-est"},
        variant_short={"v2.2dark_est": "v2.2 dark-est"},
        calib_dirs={"v2.2": DATA_ROOT / "calib_amsterdam" / "v2.2",
                    "v2.2dark_est": DATA_ROOT / "calib_amsterdam" / "v2.2dark_est"},
        calib_builder="network",
        run_dirs={"v2.2": NETWORK_V22, "v2.2dark_est": DARK_EST / "amst_dark"},
        collapse_combos=True,               # all 1064 nm: WV + wavelength conversions are no-ops
        title="Amsterdam — four co-located CHM15k: L1 + v2.2 calibration",
        subtitle="""Attenuated backscatter from the same hours and the same altitude grid for the
    four CHM15k units at Schiphol. <span class="muted">L1 <code>rcs_0</code> divided by the daily
    Kalman C<sub>L</sub> smoothed from the network eprof_v2.2 run
    (<code>calout_v22_04</code>). All four units are 1064 nm, so the water-vapour and wavelength
    corrections are exact no-ops and are not offered. No L2 archive is on disk for this station,
    so there is no distributed-product panel.</span>""",
        warnings=[
            """<b>Known instrument artefact:</b> unit B's applied static overlap function is ~25 %
    too low at 500–1000 m (back to 1 by ~1600 m), so a near-range excess of B over the other three
    units is the overlap error, not a calibration difference. The local L1 archive covers
    1 Jun – 12 Jul 2026; the page window clamps to what all four units share. """ + _DARK_EST_NOTE,
        ],
    ),

    "lindenberg": dict(
        key="lindenberg", name="Lindenberg", wmo="0-20000-0-10393",
        lat=52.21, lon=14.12, alt=123.0,
        channels=[
            dict(ident="0", itype="CHM15k", calib="rayleigh", label="CHM15k (0)", color="#1f77b4"),
            dict(ident="C", itype="CL61",   calib="cloud",    label="CL61 (C)",   color="#2ca02c"),
            # Same CL61, Rayleigh retrieval of the same constant (the Payerne C/Cr trick).
            dict(ident="C", itype="CL61",   calib="rayleigh", label="CL61 (C, Rayleigh)",
                 color="#17becf", chan="Cr"),
        ],
        iref=0,
        sources=("L1",),
        start="20260601", end="20260814",
        win_start="20260601", win_end="20260813",
        variants=["v2.2", "v2.2dark_est"],
        variant_labels={"v2.2": "v2.2 (network run)", "v2.2dark_est": "v2.2 dark-est"},
        variant_short={"v2.2dark_est": "v2.2 dark-est"},
        calib_dirs={"v2.2": DATA_ROOT / "calib_lindenberg" / "v2.2",
                    "v2.2dark_est": DATA_ROOT / "calib_lindenberg" / "v2.2dark_est"},
        calib_builder="network",
        run_dirs={"v2.2": NETWORK_V22, "v2.2dark_est": DARK_EST / "lind_dark"},
        collapse_combos=False,              # the CL61 is 910 nm: WV + wavelength combos stay live
        title="Lindenberg — CHM15k vs CL61: L1 + v2.2 calibration",
        subtitle="""Attenuated backscatter from the same hours and the same altitude grid for the
    co-located CHM15k (1064 nm, Rayleigh-calibrated) and CL61 (910.55 nm, cloud- and
    Rayleigh-calibrated). <span class="muted">L1 <code>rcs_0</code> divided by the daily Kalman
    C<sub>L</sub> smoothed from the network eprof_v2.2 run (<code>calout_v22_04</code>). The CL61
    channels are corrected to 1064 nm (water vapour + wavelength conversion, per-day CAMS). No L2
    archive is on disk for this station, so there is no distributed-product panel.</span>""",
        warnings=[
            """<b>Window and constants clamp:</b> the CL61's local daily L1 ends 30 Jun 2026, so the
    paired window is June only. Its v2.2 network-run record ends 15 Jun 2026 — the Kalman series is
    held at its last value for 16–30 Jun (clamped interpolation), on both the cloud and the Rayleigh
    channel. The CHM15k record runs through the whole window. """ + _DARK_EST_NOTE,
        ],
    ),
}


def get_site(key):
    if key not in SITES:
        raise KeyError(f"unknown site '{key}' — one of {', '.join(SITES)}")
    return SITES[key]
