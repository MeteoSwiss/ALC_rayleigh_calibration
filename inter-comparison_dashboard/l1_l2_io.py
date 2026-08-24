# -*- coding: utf-8 -*-
"""Readers + day-granular corrections for the Payerne L1-vs-L2 dashboard.

Reads the DAILY archives on D: rather than the monthly L2 mirror:

    D:/E-PROFILE_L1_2026/<wmo>/<YYYY>/<MM>/L1_<wmo>_<ident><YYYYMMDD>.nc
    D:/E-PROFILE_L2_2026/<wmo>/<YYYY>/<MM>/L2_<wmo>_<ident><YYYYMMDD>.nc

and falls back, per day, to the 5-minute granules that live in a day sub-directory:

    .../<MM>/<DD>/L2_<wmo>_<ident><YYYYMMDDHHMM>.nc      (one profile per file)

The fallback is what makes CL61 usable: its daily L2 concatenation is missing before ~2026-06-11,
and single days elsewhere are missing too. A granule carries exactly the same variables as the
daily file, so the two paths produce the same dict.

The corrections are re-implemented at DAY granularity (intercompare applies them per month, using
one monthly CAMS file). June-August 2026 have no monthly CAMS at all -- only the daily
CAMS_Beta_<YYYYMMDD>.nc cache -- so a monthly resolver would silently drop the water-vapour
correction and fall back to the US-standard molecular profile over exactly the period of interest.
Per-day CAMS is also simply more accurate than a month-mean.
"""
from __future__ import annotations
import glob
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from validation.paper import intercompare as IC
from calibration.io.cams import find_cams_file
from calibration.water_vapor_correction.water_vapor import (
    two_way_wv_transmission, cams_point_too_far)

L1_ROOT = Path(os.environ.get("ALC_DASH_L1_ROOT", "D:/E-PROFILE_L1_2026"))
L2_ROOT = Path(os.environ.get("ALC_DASH_L2_ROOT", "D:/E-PROFILE_L2_2026"))
# CAMS: the daily operational cache first (it is the only thing covering Jun-Aug 2026), then the
# monthly archives intercompare uses (0.4 deg primary, 1 deg fallback).
CAMS_DAILY = [Path(p) for p in os.environ.get(
    "ALC_DASH_CAMS_DAILY", "D:/CAMS_daily").split(";") if p]
# NOTE deliberately NOT listing D:/CAMS here: that folder holds 1-DEGREE monthlies (2018-...), and
# find_cams_file prefers a monthly over a daily within each folder -- for any window before the
# daily cache starts, the 1-degree grid (grid-point orography 894 m too high at Payerne, PWV -26 %)
# would silently win. Missing days fall back to find_cams_month = ALC_VAL_CAMS_04 (0.4 deg).


# --------------------------------------------------------------------------- file discovery
def day_files(root, level, wmo, ident, day):
    """Files covering ONE day: the concatenated daily file if it exists, else that day's 5-minute
    granules. Returns [] when the day is absent from the archive entirely."""
    mdir = root / wmo / f"{day.year}" / f"{day.month:02d}"
    daily = mdir / f"{level}_{wmo}_{ident}{day:%Y%m%d}.nc"
    if daily.exists():
        return [str(daily)]
    gl = sorted(glob.glob(str(mdir / f"{day.day:02d}" / f"{level}_{wmo}_{ident}{day:%Y%m%d}*.nc")))
    # guard against a prefix collision (…_A2026081400.nc must not match ident 'A2' etc.)
    pat = re.compile(rf"{re.escape(level)}_{re.escape(wmo)}_{re.escape(ident)}\d{{12}}\.nc$")
    return [g for g in gl if pat.search(g.replace("\\", "/"))]


def _daterange(start, end):
    d = datetime.strptime(start, "%Y%m%d")
    d1 = datetime.strptime(end, "%Y%m%d")
    while d <= d1:
        yield d
        d += timedelta(days=1)


# --------------------------------------------------------------------------- L2
def _l2_files(files):
    """Concatenate one day's L2 file(s) into raw (native-resolution) arrays."""
    times, betas, calc, qf, cbh, vv = [], [], [], [], [], []
    meta = None
    for f in files:
        try:
            with Dataset(f) as nc:
                tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
                tt = IC._decode_time(np.asarray(nc.variables["time"][:], "f8"), tu)
                a = np.asarray(nc.variables["altitude"][:], "f8")
                ab = IC._clean(nc.variables["attenuated_backscatter_0"][:])
                if ab.shape != (tt.size, a.size) and ab.shape == (a.size, tt.size):
                    ab = ab.T
                if ab.shape != (tt.size, a.size):
                    continue
                cc = IC._clean(nc.variables["calibration_constant_0"][:]).ravel() \
                    if "calibration_constant_0" in nc.variables else np.full(tt.size, np.nan)
                if cc.size != tt.size:
                    cc = np.full(tt.size, np.nan if cc.size == 0 else float(cc.ravel()[0]))
                q = IC._read2d(nc, "quality_flag", tt.size, a.size)
                cb = IC._read2d(nc, "cloud_base_height", tt.size, a.size)
                v = IC._read1d(nc, "vertical_visibility", tt.size)
                if meta is None:
                    meta = dict(alt=a,
                                salt=float(np.ravel(nc.variables["station_altitude"][:])[0]),
                                lat=float(np.ravel(nc.variables["station_latitude"][:])[0]),
                                lon=float(np.ravel(nc.variables["station_longitude"][:])[0]),
                                wl=IC._scalar(nc, "l0_wavelength", np.nan),
                                itype=IC._str(nc, "instrument_type"))
                elif a.size != meta["alt"].size:
                    continue                       # different range grid -> skip this granule
        except Exception:
            continue
        times.append(tt); betas.append(ab); calc.append(cc); qf.append(q); cbh.append(cb); vv.append(v)
    if meta is None or not times:
        return None
    return dict(meta=meta, time=np.concatenate(times), beta=np.concatenate(betas, axis=0),
                calc=np.concatenate(calc), qf=np.concatenate(qf, axis=0),
                cbh=np.concatenate(cbh, axis=0), vv=np.concatenate(vv))


def read_l2(wmo, ident, start, end, workers=12):
    """Native-resolution L2 over [start,end] from the daily archive (+5-minute fallback).
    Same dict shape as intercompare.read_l2, plus `sources` = per-day provenance counts."""
    days = list(_daterange(start, end))
    jobs = [(d, day_files(L2_ROOT, "L2", wmo, ident, d)) for d in days]
    src = {"daily": 0, "granules": 0, "missing": 0}
    for d, fs in jobs:
        src["missing" if not fs else ("daily" if len(fs) == 1 else "granules")] += 1
    parts = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(lambda j: _l2_files(j[1]) if j[1] else None, jobs):
            if r is not None:
                parts.append(r)
    if not parts:
        return None
    meta = parts[0]["meta"]
    nA = meta["alt"].size
    parts = [p for p in parts if p["meta"]["alt"].size == nA]
    time = np.concatenate([p["time"] for p in parts])
    o = np.argsort(time)
    out = dict(time=time[o], alt=meta["alt"],
               beta=np.concatenate([p["beta"] for p in parts], axis=0)[o],
               calc=np.concatenate([p["calc"] for p in parts])[o],
               qf=np.concatenate([p["qf"] for p in parts], axis=0)[o],
               cbh=np.concatenate([p["cbh"] for p in parts], axis=0)[o],
               vv=np.concatenate([p["vv"] for p in parts])[o],
               station_alt=meta["salt"], lat=meta["lat"], lon=meta["lon"],
               wavelength=meta["wl"], itype=meta["itype"], wmo=wmo, ident=ident, sources=src)
    return out


# --------------------------------------------------------------------------- L1
                     # --------------------------------------------------- dark baseline (measured)
_DARK_NPZ = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability"
                 "/dark_profiles_payerne.npz")
_DARK_CACHE = {}


def dark_profile(wmo, ident):
    """Measured electronic baseline b(z), in `rcs_0 / z^2` units, or None.

    From the Payerne dark campaign (telescopes covered, May-July 2026): the detector reads a
    NON-ZERO, NON-FLAT baseline. It matters for an inter-comparison because it is not the same
    shape on the two instruments -- measured here: CHM15k -17.0 % of the molecular night signal
    with a +1.2 %/km slope (essentially a scale error), CL61 -16.0 % with a -5.7 %/km slope
    (-3.0 % at 2 km down to -38.4 % at 6 km). The LEVELS nearly cancel in a CHM15k/CL61 ratio;
    the ~7 %/km SLOPE difference does not, and it opens with altitude -- which is exactly the
    divergence the profile panels show.

    Only Payerne (0-20000-0-06610) has been measured; other stations return None and are left
    untouched. This is a MEASURED profile, not the fitted intercept that `subtract_background`
    removes -- that one is atmosphere in disguise (see doc/reports/altitude_forward_model_study.md).
    """
    key = (str(wmo), str(ident))
    if key in _DARK_CACHE:
        return _DARK_CACHE[key]
    out = None
    if str(wmo) == "0-20000-0-06610" and _DARK_NPZ.exists():
        try:
            with np.load(_DARK_NPZ) as z:
                if f"{ident}_b" in z:
                    out = (np.asarray(z[f"{ident}_range"], "f8"), np.asarray(z[f"{ident}_b"], "f8"))
        except Exception:
            out = None
    _DARK_CACHE[key] = out
    return out


def _subtract_dark(rcs, rng, dark):
    """rcs_0 - b(z)*z^2, with b resampled onto this file's range grid."""
    if dark is None:
        return rcs
    rd, bd = dark
    b = np.interp(rng, rd, bd, left=np.nan, right=np.nan)
    ok = np.isfinite(b)
    if not ok.any():
        return rcs
    corr = np.zeros_like(rng)
    corr[ok] = b[ok] * rng[ok] ** 2
    return rcs - corr[None, :]


def _l1_files(files, dark=None):
    """One day's L1 file(s) -> hourly-median retimed arrays (mirrors intercompare._l1_day, but over
    a LIST so a day of 5-minute granules retimes as one block)."""
    ts, rcss, cbhs, vvs, temps = [], [], [], [], []
    rng = None
    for f in files:
        try:
            with Dataset(f) as nc:
                tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
                t = IC._decode_time(np.asarray(nc.variables["time"][:], "f8"), tu)
                r = np.asarray(nc.variables["range"][:], "f8")
                rcs = IC._clean(nc.variables["rcs_0"][:])
                if rcs.shape != (t.size, r.size):
                    rcs = rcs.T if rcs.shape == (r.size, t.size) else None
                if rcs is None or t.size == 0:
                    continue
                # Before any retiming or filtering: the baseline is a property of the raw counts.
                rcs = _subtract_dark(rcs, r, dark)
                if rng is None:
                    rng = r
                elif r.size != rng.size:
                    continue
                cbh = IC._read2d(nc, "cloud_base_height", t.size, r.size)
                vv = IC._read1d(nc, "vertical_visibility", t.size)
                temp = IC._read1d(nc, "temp_int", t.size) - 273.15      # K -> degC
        except Exception:
            continue
        ts.append(t); rcss.append(rcs); cbhs.append(cbh); vvs.append(vv); temps.append(temp)
    if rng is None or not ts:
        return None
    t = np.concatenate(ts)
    o = np.argsort(t)
    rcs = np.concatenate(rcss, axis=0)[o]
    cbh = np.concatenate(cbhs, axis=0)[o]
    vv = np.concatenate(vvs)[o]
    temp = np.concatenate(temps)[o]
    t = t[o]
    mode = IC.noise_filter_mode()
    snr_idx = ()
    if mode == "cloudnet":
        rcs = IC.apply_cloudnet_filter(rcs, rng)
    elif mode == "snr3":
        snr_idx = (0,)
    grid, (rcsh, cbhh, vvh, temph) = IC.retime_hourly(
        t, [rcs, cbh, vv, temp], min_cov_s=IC.MIN_AVG_S, snr_idx=snr_idx)
    return grid, rcsh, cbhh, vvh, temph, rng


def _scan_meta(jobs):
    """First readable file's station metadata + range grid (shared by the hourly and native
    readers so both see the same grid)."""
    for _, fs in jobs:
        for f in fs:
            try:
                with Dataset(f) as nc:
                    return dict(salt=float(np.ravel(nc.variables["station_altitude"][:])[0]),
                                lat=float(np.ravel(nc.variables["station_latitude"][:])[0]),
                                lon=float(np.ravel(nc.variables["station_longitude"][:])[0]),
                                wl=IC._scalar(nc, "l0_wavelength", np.nan),
                                rng=np.asarray(nc.variables["range"][:], "f8"))
            except Exception:
                continue
    return None


def _l1_files_native(files):
    """One day's L1 file(s) at NATIVE time resolution — no retiming, no dark subtraction, no
    noise filter.  The nf_v3 window-SNR statistics need the raw per-sample stream; everything
    else (decoding, fill handling, screening INPUTS cbh/vv) is identical to _l1_files.

    Each file is retried a few times: under the threaded reader, HDF5 sporadically fails a
    clean open (HDF5-DIAG "unable to open file"), and a silently skipped day makes the read
    NON-DETERMINISTIC — the build and its verification (check_v3) would then see different
    hours and the mask comparison fails on data vintage, not wiring."""
    ts, rcss, cbhs, vvs = [], [], [], []
    rng = None
    for f in files:
        got = None
        for attempt in range(3):
            try:
                with Dataset(f) as nc:
                    tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
                    t = IC._decode_time(np.asarray(nc.variables["time"][:], "f8"), tu)
                    r = np.asarray(nc.variables["range"][:], "f8")
                    rcs = IC._clean(nc.variables["rcs_0"][:])
                    if rcs.shape != (t.size, r.size):
                        rcs = rcs.T if rcs.shape == (r.size, t.size) else None
                    if rcs is None or t.size == 0:
                        got = ()
                        break
                    cbh = IC._read2d(nc, "cloud_base_height", t.size, r.size)
                    vv = IC._read1d(nc, "vertical_visibility", t.size)
                    got = (t, r, rcs, cbh, vv)
                    break
            except Exception:
                time.sleep(0.3 * (attempt + 1))
        if not got:
            continue
        t, r, rcs, cbh, vv = got
        if rng is None:
            rng = r
        elif r.size != rng.size:
            continue
        ts.append(t)
        rcss.append(np.asarray(rcs, "f4"))          # f4: 60+ days of CL61 must fit in memory
        cbhs.append(np.asarray(cbh, "f4"))
        vvs.append(np.asarray(vv, "f4"))
    if rng is None or not ts:
        return None
    t = np.concatenate(ts)
    o = np.argsort(t)
    return (t[o], np.concatenate(rcss, axis=0)[o], np.concatenate(cbhs, axis=0)[o],
            np.concatenate(vvs)[o], rng)


def read_l1_native(wmo, ident, start, end, workers=8):
    """Native-resolution twin of read_l1, for the noise-filter build (nf_v3): same daily files,
    same decoding, but the raw sample stream instead of hourly medians."""
    days = list(_daterange(start, end))
    jobs = [(d, day_files(L1_ROOT, "L1", wmo, ident, d)) for d in days]
    meta = _scan_meta(jobs)
    if meta is None:
        return None
    nR = meta["rng"].size
    parts = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(lambda j: _l1_files_native(j[1]) if j[1] else None, jobs):
            if res is not None and res[4].size == nR:
                parts.append(res[:4])
    if not parts:
        return None
    time = np.concatenate([p[0] for p in parts])
    o = np.argsort(time)
    return dict(time=time[o], alt=meta["rng"] + meta["salt"],
                beta=np.concatenate([p[1] for p in parts], axis=0)[o],
                cbh=np.concatenate([p[2] for p in parts], axis=0)[o],
                vv=np.concatenate([p[3] for p in parts])[o],
                station_alt=meta["salt"], lat=meta["lat"], lon=meta["lon"],
                wavelength=meta["wl"], wmo=wmo, ident=ident)


def read_l1(wmo, ident, start, end, workers=12):
    """Hourly-median L1 rcs_0 over [start,end] from the daily archive (+5-minute fallback).
    Same dict shape as intercompare.read_l1, plus `sources`."""
    days = list(_daterange(start, end))
    jobs = [(d, day_files(L1_ROOT, "L1", wmo, ident, d)) for d in days]
    src = {"daily": 0, "granules": 0, "missing": 0}
    for d, fs in jobs:
        src["missing" if not fs else ("daily" if len(fs) == 1 else "granules")] += 1
    meta = _scan_meta(jobs)
    if meta is None:
        return None
    nR = meta["rng"].size
    parts = []
    # ALC_DASH_DARK=1 forces the baseline subtraction at CACHE time -- default OFF: the dashboard's
    # "v2.2dark" variant applies it per-combo instead (exact, since the baseline is time-constant
    # and the cache holds hourly MEDIANS), which lets one cache serve corrected and uncorrected
    # views. Correcting only the profile while dividing by an uncorrected constant mixes two signal
    # definitions and WORSENS the comparison (measured 2026-08-15) -- the variant corrects both.
    dark = dark_profile(wmo, ident) if os.environ.get("ALC_DASH_DARK", "0") == "1" else None
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(lambda j: _l1_files(j[1], dark) if j[1] else None, jobs):
            if res is not None and res[5].size == nR:
                parts.append(res[:5])
    if not parts:
        return None
    time = np.concatenate([p[0] for p in parts])
    o = np.argsort(time)
    return dict(time=time[o], alt=meta["rng"] + meta["salt"],
                beta=np.concatenate([p[1] for p in parts], axis=0)[o], calc=None,
                qf=np.zeros((time.size, nR), "f8"),
                cbh=np.concatenate([p[2] for p in parts], axis=0)[o],
                vv=np.concatenate([p[3] for p in parts])[o],
                temp_int=np.concatenate([p[4] for p in parts])[o],
                station_alt=meta["salt"], lat=meta["lat"], lon=meta["lon"],
                wavelength=meta["wl"], itype="", wmo=wmo, ident=ident, sources=src)


# --------------------------------------------------------------------------- day-granular CAMS
_DAY_CAMS = {}


def cams_for_day(day):
    """CAMS file covering ONE day: the daily operational cache first, then the monthly archives.
    Returns None if neither has it (the caller then excludes / falls back, as before)."""
    ds = day.strftime("%Y%m%d") if hasattr(day, "strftime") else str(day)
    if ds in _DAY_CAMS:
        return _DAY_CAMS[ds]
    f = None
    for folder in CAMS_DAILY:
        f = find_cams_file(folder, ds)
        if f is not None:
            break
    if f is None:
        f = IC.find_cams_month(int(ds[:4]), int(ds[4:6]))
    _DAY_CAMS[ds] = f
    return f


def _day_groups(t):
    """[(YYYYMMDD, boolean row mask)] for a time axis."""
    days = pd.to_datetime(t).normalize()
    uniq = pd.unique(days)
    return [(pd.Timestamp(u), np.asarray(days == u)) for u in uniq]


def _win(day):
    """CAMS extraction window around a day (the night before included, as the monthly path does)."""
    d = np.datetime64(pd.Timestamp(day).date())
    return d - np.timedelta64(1, "D"), d + np.timedelta64(2, "D")


def apply_wv(beta, d, lam0, fwhm):
    """Two-way water-vapour transmission correction, resolved per DAY. Mirrors
    intercompare.apply_wv (including the 'CAMS grid point too far' guard, which excludes rather
    than reports uncorrected) but never falls back to a month-mean profile."""
    z_asl = np.asarray(d["alt"], "f8")
    out = beta.copy()
    info = {"days": 0, "days_excluded": 0, "median_t2": []}
    for day, sel in _day_groups(d["time"]):
        cams = cams_for_day(day)
        if cams is None or cams_point_too_far(cams, d["lat"], d["lon"]):
            out[sel, :] = np.nan; info["days_excluded"] += 1; continue
        t0, t1 = _win(day)
        prof = IC._cams_wv_profile_cached(str(cams), d["lat"], d["lon"], t0, t1)
        if prof is None:
            out[sel, :] = np.nan; info["days_excluded"] += 1; continue
        h_wv, n_wv = prof
        t2 = np.asarray(two_way_wv_transmission(z_asl, d["station_alt"], h_wv, n_wv,
                                                IC.WV_LUT, lam0, fwhm), "f8")
        if t2.size == z_asl.size and np.any(np.isfinite(t2)):
            out[sel, :] = beta[sel, :] / t2[None, :]
            info["days"] += 1; info["median_t2"].append(float(np.nanmedian(t2)))
        else:
            out[sel, :] = np.nan; info["days_excluded"] += 1
    return out, info


def wavelength_correct_molecular(beta, d, lam, target, alpha):
    """Component-separated lam->target conversion with the molecular term from CAMS T/p, resolved
    per DAY. A day with no CAMS falls back to the US-standard molecular profile, as before."""
    if not np.isfinite(lam) or abs(lam - target) < 1.0:
        return beta
    alt_asl = np.asarray(d["alt"], "f8")
    z_agl = alt_asl - d["station_alt"]
    f = (lam / target) ** alpha
    out = beta.copy()
    us_l = us_t = None
    for day, sel in _day_groups(d["time"]):
        cams = cams_for_day(day)
        bml = bmt = None
        if cams is not None:
            t0, t1 = _win(day)
            bml = IC._mol_att_cams(str(cams), d["lat"], d["lon"], t0, t1, alt_asl,
                                   d["station_alt"], lam)
            bmt = IC._mol_att_cams(str(cams), d["lat"], d["lon"], t0, t1, alt_asl,
                                   d["station_alt"], target)
        if bml is None or bmt is None:
            if us_l is None:
                us_l = IC._molecular_beta(z_agl, d["station_alt"], lam)
                us_t = IC._molecular_beta(z_agl, d["station_alt"], target)
            bml, bmt = us_l, us_t
        out[sel, :] = bmt[None, :] + (beta[sel, :] - bml[None, :]) * f
    return out
