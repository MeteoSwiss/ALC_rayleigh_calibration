# -*- coding: utf-8 -*-
"""Shared configuration of the reference-site agreement study (ref23).

Question: are co-located, independently calibrated instruments consistent with each other within
the per-night v2.3 calibration uncertainty (doc/v2.3/uncertainty.md), what is their bias, and is
the v2.3 uncertainty better calibrated than the v2.2 ensemble uncertainty against this independent
truth?  First of all at Payerne, against the CL61 (cloud calibration, as in the paper).

This module only holds the configuration (sites, streams, pairs, paths, display grid) and small
helpers shared by the harness (ref23_run.py) and the analysis (ref23_analysis.py).
"""
from __future__ import annotations

import json
from datetime import date, datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = Path("C:/DATA/Projects/202610_v23/reference_agreement")
RUN_DIR = OUT / "runs"
# Rayleigh re-run with the REVIEWED v2.3 budget (window-mean noise, vertical cap 15 %, molecular
# climatology at the window height ASL; merged 2026-10-21). RUN_DIR holds the pre-review budget.
RUN_DIR_REVIEWED = OUT / "runs_rayleigh_reviewed"
FIG_DIR = OUT
L1_ROOT = Path("D:/E-PROFILE_L1_2026")
# 0.4 deg CAMS only (monthly 2025-01..2026-06, daily 2026-06-01..): never a 1-degree folder
# (D:/CAMS, D:/CAMS_run_v20: the grid-point orography is ~900 m off at Payerne).
CAMS = "D:/CAMS_Monthly_04;D:/CAMS_daily"
CAMS_MONTHLY_04 = "D:/CAMS_Monthly_04"
CAMS_DAILY = "D:/CAMS_daily"
CENSUS = REPO / "validation" / "scope_l1_2026_census.json"
WV_LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
# The published per-night constants (E-PROFILE dashboard NetCDFs, fetched 2026-10-09 by
# validation/paper/fetch_dashboard_calib.py): Rayleigh calibration_version 220 (eprof_v2.2,
# ensemble uncertainty, 2 sigma) and cloud calibration_version 1000 (O'Connor).
PUBLISHED_NC = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/"
                    "paper_v08/calib_dashboard/nc")
# The operational CL61 water-vapour line (manufacturer spectrum, operational since 2026-08).
WV_SPECTRUM = '{"CL61": [910.55, 0.188]}'

START = date(2025, 1, 1)
EPOCH = datetime(1970, 1, 1)
SUCCESS = (1.0, 0.5)

# ---------------------------------------------------------------------------------- sites
# streams: ident -> (instrument type, calibration methods run for it)
SITES = {
    "payerne": dict(wmo="0-20000-0-06610", name="Payerne",
                    streams={"A": ("CHM15k", ("ray",)), "B": ("CL31", ("cloud",)),
                             "C": ("CL61", ("ray", "cloud"))}),
    "lindenberg": dict(wmo="0-20000-0-10393", name="Lindenberg",
                       streams={"0": ("CHM15k", ("ray",)), "C": ("CL61", ("ray", "cloud"))}),
    "aosta": dict(wmo="0-380-5-1", name="Aosta",
                  streams={"0": ("CHM15k", ("ray",)), "B": ("CL61", ("ray", "cloud"))}),
    "schiphol": dict(wmo="0-20000-0-06240", name="Amsterdam Schiphol",
                     streams={i: ("CHM15k", ("ray",)) for i in "ABCD"}),
    "budapest": dict(wmo="0-20000-0-12843", name="Budapest",
                     streams={i: ("CHM15k", ("ray",)) for i in "AB"}),
    "magurele": dict(wmo="0-20008-0-INO", name="Magurele",
                     streams={i: ("CHM15k", ("ray",)) for i in "AB"}),
}

# Pairs X / R: ln rho = ln(beta_X / beta_R).  R is the reference (the CL61 at the mixed sites).
# A channel is (ident, method); "ray" = Rayleigh, "cloud" = liquid-cloud (O'Connor).
PAIRS = {
    "payerne": [(("A", "ray"), ("C", "cloud")), (("A", "ray"), ("C", "ray")),
                (("C", "ray"), ("C", "cloud")), (("B", "cloud"), ("C", "cloud")),
                (("A", "ray"), ("B", "cloud"))],
    "lindenberg": [(("0", "ray"), ("C", "cloud")), (("0", "ray"), ("C", "ray")),
                   (("C", "ray"), ("C", "cloud"))],
    "aosta": [(("0", "ray"), ("B", "cloud")), (("0", "ray"), ("B", "ray")),
              (("B", "ray"), ("B", "cloud"))],
    "schiphol": [((a, "ray"), (b, "ray")) for a, b in
                 (("B", "A"), ("C", "A"), ("D", "A"), ("C", "B"), ("D", "B"), ("D", "C"))],
    "budapest": [(("B", "ray"), ("A", "ray"))],
    "magurele": [(("B", "ray"), ("A", "ray"))],
}

# Instrument pairs whose L1 ratio ingredients the harness prepares (the cloud/Rayleigh pair of one
# instrument needs no L1: rho = C_ray / C_cloud for the same photons).
def l1_pairs(site_key):
    out = []
    for (x, _mx), (r, _mr) in PAIRS[site_key]:
        if x != r and (x, r) not in out:
            out.append((x, r))
    return out


# ---------------------------------------------------------------------------------- grid
# Display grid of the v3 inter-comparison dashboard (inter-comparison_dashboard/build_v3.py, Z_AGL):
# 30 m up to 3 km (the 500-3000 m statistics band is bit-identical to the v2 uniform grid),
# 90 m up to 8 km.  Heights are m above ground (the instruments point vertically: range = AGL).
Z_AGL = np.concatenate([np.arange(0.0, 3000.0 + 30.0, 30.0),
                        np.arange(3090.0, 8010.0 + 90.0, 90.0)])
BANDS = {"b500_3000": (500.0, 3000.0),      # v3 dashboard statistics band
         "b1500_3000": (1500.0, 3000.0)}    # above full overlap of every unit
WINDOWS_S = (300, 3600)                      # SNR evaluation windows (v3 default 5 min; 60 min)
SNR_MIN = 3.0
NF_FMIN = 0.5                                # an hour is admitted when >= 50 % of its data detect
MIN_COV_S = 1800.0                           # hourly median needs >= 30 min of night samples
CLOUD_EXPAND_S = 900.0                       # cloud/fog exclusion expanded by +-15 min (v3 screen)
TARGET_NM = 1064.0                           # common wavelength of the comparison (v3)
ALPHAS = (0.0, 0.5, 1.0, 1.5, 2.0)           # Angstrom exponents of the lambda-conversion test
ALPHA_NOMINAL = 1.0


# ---------------------------------------------------------------------------------- helpers
def census():
    return {(s["wmo"], s["ident"]): s for s in json.load(open(CENSUS, encoding="utf-8"))}


def l1_file(wmo, ident, d):
    return L1_ROOT / wmo / f"{d:%Y}" / f"{d:%m}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"


def site_nights(site_key, d0=START, d1=None):
    """Nights D (date) for which at least one stream of the site has its day-D L1 file."""
    s = SITES[site_key]
    d1 = d1 or date.today()
    out = []
    d = d0
    while d <= d1:
        if any(l1_file(s["wmo"], i, d).exists() for i in s["streams"]):
            out.append(d)
        d += timedelta(days=1)
    return out


def ymd(d):
    return d.strftime("%Y%m%d")


def season_of(month):
    return {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
            6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}[int(month)]


def night_mask(time_days, lat, lon, date_str, sza_threshold=100.0):
    """The operational night of ``date_str`` (calibration.io.data_loader.filter_time_range with
    use_sza_night): solar day D-1 after 12:00 solar or solar day D before 12:00 solar, and solar
    zenith angle > sza_threshold, with the declination of day D.  Returns a boolean mask."""
    t = np.asarray(time_days, float)
    utc_days = np.floor(t)
    utc_hours = (t - utc_days) * 24.0
    solar_hours = utc_hours + lon / 15.0
    solar_days = utc_days.copy()
    fwd = solar_hours >= 24.0
    bwd = solar_hours < 0.0
    solar_hours[fwd] -= 24.0
    solar_days[fwd] += 1.0
    solar_hours[bwd] += 24.0
    solar_days[bwd] -= 1.0
    cur = datetime.strptime(date_str, "%Y%m%d")
    d_cur = float((cur - EPOCH).days)
    in_night = (((solar_days == d_cur - 1.0) & (solar_hours >= 12.0))
                | ((solar_days == d_cur) & (solar_hours < 12.0)))
    doy = cur.timetuple().tm_yday
    decl = np.radians(-23.44) * np.cos(2.0 * np.pi * (doy + 10) / 365.25)
    la = np.radians(lat)
    ha = np.radians((solar_hours - 12.0) * 15.0)
    cos_sza = np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)
    sza = np.degrees(np.arccos(np.clip(cos_sza, -1.0, 1.0)))
    return in_night & (sza > sza_threshold)
