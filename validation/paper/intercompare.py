"""
intercompare.py — Python port of the MATLAB paper_val_process.m multi-channel inter-comparison.
For each channel: read L2_monthly attenuated backscatter, apply the (Kalman) calibration produced by
calib_benchmark.py, the water-vapour correction (910 nm), and the wavelength normalisation; build a
quality-screened stream; retime to an hourly median grid; average onto a common altitude grid; and
compute bias / RMSE / relative-bias / Pearson-r vs a reference channel over 500-3000 m AGL.

Faithful to _PORT_SPEC.md sections 1.4-1.11. Reuses the operational calibration package for the WV and
molecular pieces. Returns an R dict mirroring the MATLAB R struct (channels, altGrid, time_sync, beta,
beta_disp, cbh_disp, station, stats).
"""
from __future__ import annotations
import csv
import os
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from calibration.io.cams import ensure_cams_file, find_cams_file
from calibration.cloud.calibration import INSTRUMENT_CAL_DEFAULT
from calibration.water_vapor_correction.water_vapor import (
    cams_water_vapor_profile, two_way_wv_transmission, laser_spectrum_for, in_water_vapor_band,
    cams_point_too_far, cams_temperature_pressure_profile)
from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties
from calibration.config import InstrumentType

REPO = Path(__file__).resolve().parents[2]
L2_MONTHLY = Path("A:/E-PROFILE_L2_monthly")
CALIB = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/calib")
# CAMS for the WV + molecular (Rayleigh) profiles: 0.4 deg monthly primary, 1 deg monthly fallback for
# months the 0.4 deg archive has not covered yet (both L137 model levels). Env-overridable.
CAMS = "D:/CAMS"
CAMS_04 = Path(os.environ.get("ALC_VAL_CAMS_04", "D:/CAMS_Monthly_04"))
CAMS_1DEG = Path(os.environ.get("ALC_VAL_CAMS_1DEG", "D:/CAMS"))
# Daily 0.4 deg CAMS (CAMS_Beta_<YYYYMMDD>.nc, the operational prefetch) — the only source for recent
# months whose monthly file does not exist yet (e.g. Jun-Aug 2026). Used per day as a last fallback.
CAMS_DAILY = Path(os.environ.get("ALC_VAL_CAMS_DAILY", "D:/CAMS_daily"))
_G0, _RD, _LAPSE = 9.80665, 287.058, 0.0065   # std gravity, dry-air gas const, std lapse rate
# --- L1 data path (validate by applying the CSCS calibration to the raw range-corrected signal) ---
# rcs_0 is read from the native daily L1 archive and turned into attenuated backscatter with
# beta_att [Mm^-1 sr^-1] = rcs_0 / C_L * 1e6, where C_L is the CSCS Kalman lidar constant from the
# calout <key>_kalman.csv. (No L2 product, no raw vendor files.) Both env-overridable.
L1_ROOT = Path(os.environ.get("ALC_VAL_L1_ROOT", "D:/E-PROFILE_L1_2026"))

# --- Noise filtering mode for the L1 validation ---------------------------------------------------
# The detection gate is applied on the NATIVE (sub-hourly) profiles inside _l1_day, BEFORE the hourly
# median — the only stage where the raw profiles still exist. Three modes, selectable at run time:
#   "none"     (default) — no gate. Symmetric noise averages out in the hourly median; keeping the
#              UNGATED signal avoids conditioning on positive noise (which biases the retained median
#              HIGH at altitude, instrument-dependently). Detection limits are a separate sensitivity
#              profile, not a censoring gate.
#   "snr3"     — per-gate window SNR>=3 (|median|/(robust_sigma/sqrt n)) inside retime_hourly; the
#              operational sensitivity gate. (Legacy switch ALC_VAL_L1_SNR=1 still selects this.)
#   "cloudnet" — the Cloudnet / ceilopyter per-profile filter (ceilopyter.noise.remove_noise): mask
#              each sample whose non-range-corrected signal rcs_0/r^2 falls below CLOUDNET_SNR_LIMIT x
#              the per-profile far-range std, then median the survivors.
# The mode is read LIVE from ALC_VAL_NOISE_FILTER (or set_noise_filter()) so a runner can switch it
# before dispatching work — including across the site ProcessPool, whose Windows-spawn children
# inherit os.environ. Cloudnet knobs: ALC_VAL_CLOUDNET_SNR (default 5), ALC_VAL_CLOUDNET_FLOOR (0).
_VALID_FILTERS = ("none", "snr3", "cloudnet")


def noise_filter_mode():
    """Active noise-filter mode ('none'|'snr3'|'cloudnet'), read live from the environment."""
    m = os.environ.get("ALC_VAL_NOISE_FILTER", "").strip().lower()
    if not m:
        m = "snr3" if os.environ.get("ALC_VAL_L1_SNR", "0") != "0" else "none"
    return m if m in _VALID_FILTERS else "none"


def set_noise_filter(mode):
    """Select the noise-filter mode programmatically (sets ALC_VAL_NOISE_FILTER for this process and
    any ProcessPool children spawned afterwards). Returns the normalised mode."""
    mode = str(mode).strip().lower()
    if mode not in _VALID_FILTERS:
        raise ValueError(f"noise filter must be one of {_VALID_FILTERS}, got {mode!r}")
    os.environ["ALC_VAL_NOISE_FILTER"] = mode
    return mode
CALOUT = Path(os.environ.get(
    "ALC_VAL_CALOUT", "C:/DATA/Projects/202606_E-PROFILE_calibration/E_PROFILE_calout_2025_2026"))
WV_LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
STD_ATM = REPO / "calibration" / "data" / "standard_atmosphere_US_1976_50km.csv"
# CL61: the manufacturer emission line (910.55 nm, sigma 0.08 -> FWHM 0.188 nm) is the OPERATIONAL
# spectrum since 2026-08 (E-PROFILE per-night NetCDFs reproduce it night-for-night); the historical
# Qmini 910.74/1.0 over-corrects WV ~x9 and is retired. The beta WV correction must use the SAME
# spectrum as the calibration constant it is paired with.
WV_PARAMS = {"CL31": (909.7, 6.0), "CL51": (910.0, 3.4), "CL61": (910.55, 0.188),
             "CHM15k": (1064.47, 0.5), "Mini-MPL": (532.0, 0.1)}


# --------------------------------------------------------------------------- L2 reading
def read_l2(wmo, ident, start, end):
    """Read + concatenate L2_monthly files. beta = attenuated_backscatter_0 [Mm^-1 sr^-1]
    (stored 1e-6 1/(m sr) == 1 Mm^-1 sr^-1). Returns a dict of arrays (time x range)."""
    d0 = datetime.strptime(start, "%Y%m%d"); d1 = datetime.strptime(end, "%Y%m%d")
    months = pd.period_range(d0, d1, freq="M")
    times, betas, calc, qf, cbh, vv = [], [], [], [], [], []
    alt = lat = lon = salt = None; wl = np.nan; itype = ""
    for p in months:
        f = L2_MONTHLY / wmo / f"{p.year}" / f"L2_{wmo}_{ident}{p.year}{p.month:02d}.nc"
        if not f.exists():
            continue
        with Dataset(f) as nc:
            t = np.asarray(nc.variables["time"][:], "f8")  # days since 1970 (MATLAB datenum-like? assume CF)
            tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
            tt = _decode_time(t, tu)
            ab = _clean(nc.variables["attenuated_backscatter_0"][:])
            a = np.asarray(nc.variables["altitude"][:], "f8")
            if ab.shape[0] == a.size and ab.shape[0] != tt.size:
                ab = ab.T  # -> (time, range)
            elif ab.shape != (tt.size, a.size) and ab.shape == (a.size, tt.size):
                ab = ab.T
            cc = _clean(nc.variables["calibration_constant_0"][:]).ravel()
            q = _read2d(nc, "quality_flag", tt.size, a.size)
            cb = _read2d(nc, "cloud_base_height", tt.size, a.size)
            v = _read1d(nc, "vertical_visibility", tt.size)
            if alt is None:
                alt = a; salt = float(np.ravel(nc.variables["station_altitude"][:])[0])
                lat = float(np.ravel(nc.variables["station_latitude"][:])[0])
                lon = float(np.ravel(nc.variables["station_longitude"][:])[0])
                wl = _scalar(nc, "l0_wavelength", np.nan)
                itype = _str(nc, "instrument_type")
        times.append(tt); betas.append(ab); calc.append(cc); qf.append(q); cbh.append(cb); vv.append(v)
    if not times:
        return None
    time = np.concatenate(times)
    beta = np.concatenate(betas, axis=0)
    out = dict(time=time, alt=alt, beta=beta, calc=np.concatenate(calc),
               qf=np.concatenate(qf, axis=0), cbh=np.concatenate(cbh, axis=0),
               vv=np.concatenate(vv), station_alt=salt, lat=lat, lon=lon, wavelength=wl, itype=itype,
               wmo=wmo, ident=ident)
    o = np.argsort(time)
    for k in ("time", "calc", "vv"):
        out[k] = out[k][o]
    for k in ("beta", "qf", "cbh"):
        out[k] = out[k][o]
    return out


def _decode_time(t, units):
    # CF "days since 1970-01-01" (E-PROFILE L2). Return datetime64[ns].
    base = np.datetime64("1970-01-01")
    if "second" in units:
        return base + (t * 1e9).astype("timedelta64[ns]")
    return base + (t * 86400e9).astype("timedelta64[ns]")


def _clean(v):
    """netCDF var -> float array with _FillValue (auto-mask or |x|>1e30) and non-finite -> NaN."""
    if np.ma.isMaskedArray(v):
        x = np.ma.filled(v.astype("f8"), np.nan)
    else:
        x = np.asarray(v, "f8")
    x = np.where(np.abs(x) > 1e30, np.nan, x)
    return x


def _read2d(nc, name, nt, na):
    if name not in nc.variables:
        return np.full((nt, na), np.nan)
    x = _clean(nc.variables[name][:])
    if x.ndim == 1:
        return x.reshape(nt, -1) if x.size == nt else np.full((nt, na), np.nan)
    if x.shape[0] != nt and x.shape[1] == nt:
        x = x.T
    return x


def _read1d(nc, name, nt):
    if name not in nc.variables:
        return np.full(nt, np.nan)
    x = _clean(nc.variables[name][:]).ravel()
    return x if x.size == nt else np.full(nt, np.nan)


def _scalar(nc, name, default):
    if name in nc.variables:
        v = np.ravel(np.asarray(nc.variables[name][:], "f8"))
        if v.size and np.isfinite(v[0]):
            return float(v[0])
    return default


def _str(nc, name):
    if name in nc.variables:
        try:
            return "".join(np.asarray(nc.variables[name][:]).astype(str).ravel()).strip()
        except Exception:
            return ""
    return ""


# --------------------------------------------------------------------------- calibration
def load_calib_series(key, level=None):
    """Read the Kalman calibration series for a channel. level in {'L1','L2',None}: prefer
    <key>_<level>.csv, then the L2 eprof_v2 series, then the legacy un-suffixed <key>.csv (so an
    L1 request for a channel with no L1 record falls back to L2, NOT to the old eprof_v1.2 file)."""
    cands = []
    if level:
        cands.append(CALIB / f"{key}_{level}.csv")
    cands.append(CALIB / f"{key}_L2.csv")
    cands.append(CALIB / f"{key}.csv")
    f = next((c for c in cands if c.exists()), None)
    if f is None:
        return None
    rows = list(csv.DictReader(open(f, encoding="utf-8")))
    dd = np.array([np.datetime64(r["time"][:10]) for r in rows])
    ck = np.array([float(r["C_kalman"]) if r["C_kalman"] not in ("", "nan") else np.nan for r in rows])
    m = np.isfinite(ck)
    return (dd[m], ck[m]) if m.any() else None


def interp_calib(dd, val, t):
    """Linear interp of daily calibration onto profile times, clamped at the ends."""
    tf = t.astype("datetime64[ns]").astype("f8")
    df = dd.astype("datetime64[ns]").astype("f8")
    if val.size == 1:
        return np.full(t.size, val[0])
    v = np.interp(tf, df, val)            # np.interp already clamps to end values
    return v


# --------------------------------------------------------------------------- L1 reading (rcs_0)
def _l1_day(fp):
    """Read ONE daily L1 file and median-retime rcs_0 / cloud_base_height / vertical_visibility /
    internal temperature to its hourly grid (keeps memory bounded over 2 years). temp_int is
    converted K->degC for the CHM15k overlap correction. Returns (grid, rcs_h, cbh_h, vv_h, temp_h,
    range) or None."""
    try:
        with Dataset(fp) as nc:
            tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
            t = _decode_time(np.asarray(nc.variables["time"][:], "f8"), tu)
            rng = np.asarray(nc.variables["range"][:], "f8")
            rcs = _clean(nc.variables["rcs_0"][:])
            if rcs.shape != (t.size, rng.size):
                rcs = rcs.T if rcs.shape == (rng.size, t.size) else None
            if rcs is None or t.size == 0:
                return None
            cbh = _read2d(nc, "cloud_base_height", t.size, rng.size)
            vv = _read1d(nc, "vertical_visibility", t.size)
            temp = _read1d(nc, "temp_int", t.size) - 273.15   # internal temperature K -> degC (NaN if absent)
    except Exception:
        return None
    # Native-profile noise filtering (see noise_filter_mode): the Cloudnet filter masks noisy samples
    # in place (then the hourly median below runs over the survivors); the SNR3 gate is applied per
    # hourly window inside retime_hourly. "none" leaves the ungated signal.
    mode = noise_filter_mode()
    snr_idx = ()
    if mode == "cloudnet":
        rcs = apply_cloudnet_filter(rcs, rng)
    elif mode == "snr3":
        snr_idx = (0,)
    grid, (rcsh, cbhh, vvh, temph) = retime_hourly(t, [rcs, cbh, vv, temp], min_cov_s=MIN_AVG_S,
                                                   snr_idx=snr_idx)
    return grid, rcsh, cbhh, vvh, temph, rng


def read_l1(wmo, ident, start, end, workers=16):
    """Read native L1 rcs_0 over [start,end] (daily files), hourly-median retimed. Returns the SAME
    dict shape as read_l2 with beta := rcs_0 (UNCALIBRATED) and alt = range + station_altitude (L1
    range is AGL). No calibration_constant_0 / quality_flag in L1 (qf -> zeros)."""
    d0 = datetime.strptime(start, "%Y%m%d"); d1 = datetime.strptime(end, "%Y%m%d")
    files, d = [], d0
    while d <= d1:
        f = L1_ROOT / wmo / f"{d.year}" / f"{d.month:02d}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"
        if f.exists():
            files.append(str(f))
        d += timedelta(days=1)
    if not files:
        return None
    meta = None
    for f in files:
        try:
            with Dataset(f) as nc:
                meta = dict(salt=float(np.ravel(nc.variables["station_altitude"][:])[0]),
                            lat=float(np.ravel(nc.variables["station_latitude"][:])[0]),
                            lon=float(np.ravel(nc.variables["station_longitude"][:])[0]),
                            wl=_scalar(nc, "l0_wavelength", np.nan),
                            rng=np.asarray(nc.variables["range"][:], "f8"))
            break
        except Exception:
            continue
    if meta is None:
        return None
    nR = meta["rng"].size
    parts = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(_l1_day, files):
            if res is None or res[5].size != nR:   # skip unreadable / different range grid
                continue
            parts.append(res[:5])
    if not parts:
        return None
    time = np.concatenate([p[0] for p in parts])
    rcs = np.concatenate([p[1] for p in parts], axis=0)
    cbh = np.concatenate([p[2] for p in parts], axis=0)
    vv = np.concatenate([p[3] for p in parts])
    temp = np.concatenate([p[4] for p in parts])
    o = np.argsort(time)
    return dict(time=time[o], alt=meta["rng"] + meta["salt"], beta=rcs[o], calc=None,
                qf=np.zeros(rcs.shape, dtype="f8")[o], cbh=cbh[o], vv=vv[o], temp_int=temp[o],
                station_alt=meta["salt"], lat=meta["lat"], lon=meta["lon"],
                wavelength=meta["wl"], itype="", wmo=wmo, ident=ident)


def load_calout_kalman(wmo, ident, method):
    """C_L Kalman series (Wiegner lidar constant) straight from the CSCS calout <key>_kalman.csv,
    method-filtered (rayleigh|cloud) -> (dates, C_L). For L1: beta_att = rcs_0 / C_L. Both methods
    store C_L in the calout, so this is a true lidar constant (NOT the report's O'Connor C)."""
    f = CALOUT / f"{wmo}_{ident}" / f"{wmo}_{ident}_kalman.csv"
    if not f.is_file():
        return None
    dd, ck = [], []
    for r in csv.DictReader(open(f, encoding="utf-8")):
        if r.get("method") != method:
            continue
        v = r.get("kalman")
        if v in ("", "nan", None):
            continue
        try:
            val = float(v)
        except ValueError:
            continue
        s = str(r.get("date", ""))
        if val > 0 and val == val and len(s) >= 8:
            dd.append(np.datetime64(f"{s[:4]}-{s[4:6]}-{s[6:8]}")); ck.append(val)
    if not dd:
        return None
    o = np.argsort(dd)
    return np.array(dd)[o], np.array(ck)[o]


# --------------------------------------------------------------------------- WV correction
def cams_periods(t):
    """Split the profile times t into CAMS periods: (sel, cams_file_or_None, tstart, tend, label).
    A month with a monthly file (0.4 deg, else 1 deg) is one period averaged over the month; a month
    without one is split into DAYS served by the daily 0.4 deg files (CAMS_DAILY, the file of the day
    itself, else of the next day — a daily file holds the night ending on its date). A period with no
    file at all is returned with cams_file=None so the caller can NaN-mask it."""
    tt = pd.to_datetime(t)
    months = pd.PeriodIndex(tt, freq="M")
    for p in months.unique():
        selm = np.asarray(months == p)
        cams = find_cams_month(p.year, p.month)
        if cams is not None:
            tstart = np.datetime64(f"{p.year}-{p.month:02d}-01") - np.timedelta64(1, "D")
            tend = np.datetime64(f"{p.year}-{p.month:02d}-01") + np.timedelta64(40, "D")
            yield selm, cams, tstart, tend, str(p)
            continue
        days = pd.PeriodIndex(tt, freq="D")
        for d in days[selm].unique():
            sel = np.asarray(days == d)
            d64 = np.datetime64(f"{d.year}-{d.month:02d}-{d.day:02d}")
            f = None
            for cand in (d64, d64 + np.timedelta64(1, "D")):
                f = find_cams_file(CAMS_DAILY, str(cand).replace("-", ""))
                if f is not None:
                    break
            yield sel, f, d64, d64 + np.timedelta64(1, "D"), str(d)


def apply_wv(beta, l2, lam0, fwhm):
    """Divide beta by the two-way WV transmission, per CAMS period (monthly file, else daily files;
    no CAMS -> NaN)."""
    z_asl = l2["alt"]
    z_agl = z_asl - l2["station_alt"]
    out = beta.copy()
    info = {"months": [], "months_excluded": [], "median_t2": []}
    for sel, cams, tstart, tend, p in cams_periods(l2["time"]):
        if cams is None:
            out[sel, :] = np.nan
            info["months_excluded"].append(str(p)); continue
        # same guard as the calibration (flag -10 there): NEVER correct with a far-away
        # domain-edge grid point — exclude the data instead of reporting it uncorrected.
        if cams_point_too_far(cams, l2["lat"], l2["lon"]):
            out[sel, :] = np.nan
            info["months_excluded"].append(f"{p} (CAMS too far)"); continue
        prof = _cams_wv_profile_cached(cams, l2["lat"], l2["lon"], tstart, tend)
        if prof is None:
            out[sel, :] = np.nan; info["months_excluded"].append(str(p)); continue
        h_wv, n_wv = prof
        t2 = two_way_wv_transmission(l2["station_alt"] + z_agl, l2["station_alt"], h_wv, n_wv,
                                     WV_LUT, lam0, fwhm)
        t2 = np.asarray(t2, "f8")
        if t2.size == z_agl.size and np.any(np.isfinite(t2)):
            out[sel, :] = beta[sel, :] / t2[None, :]
            info["months"].append(str(p)); info["median_t2"].append(float(np.nanmedian(t2)))
        else:
            out[sel, :] = np.nan; info["months_excluded"].append(str(p))
    return out, info


# --------------------------------------------------------------------------- wavelength
def _molecular_beta(z_agl, station_alt, wavelength_nm):
    """Molecular ATTENUATED backscatter beta_mol*T^2_mol [Mm^-1 sr^-1] on z_agl from the US
    standard atmosphere, evaluated at the station's ABSOLUTE altitude (z_agl + station_alt),
    with the two-way molecular transmission integrated from the instrument upward — the same
    convention as _mol_att_cams. (Fixed 2026-07: the previous version ignored station_alt and
    interpolated the sea-level profile at z_AGL, overestimating beta_mol by ~exp(alt/8400)
    (Payerne +6 %) and including the below-station layer in T^2. The error largely cancels in
    the 910/1064 M-term (<1 % of signal) but not for the 532 nm molaer path.)
    The two-way molecular transmission matters: the measured signal the
    molaer model subtracts from is attenuated, and at 532 nm T^2_mol is already ~0.85-0.90 by
    2-3 km — subtracting the UNattenuated beta_mol there biases the extracted aerosol low by
    ~beta_mol*(1-T^2), which after the 532->1064 recombination produced a spurious ~-40 % on the
    Mini-MPL (the native 532 nm EARLINET comparison shows the instrument itself is within a few %)."""
    grid = np.arange(0.0, 20001.0, 30.0)                       # ASL grid (US-Std file reaches 50 km)
    atm = load_standard_atmosphere(STD_ATM, grid)
    z_asl = np.asarray(z_agl, "f8") + float(station_alt)
    T = np.interp(z_asl, grid, atm.temperature, left=np.nan, right=np.nan)
    P = np.interp(z_asl, grid, atm.pressure, left=np.nan, right=np.nan)
    mol = calculate_molecular_properties(T, P, np.asarray(z_agl, "f8"), wavelength_nm * 1e-9)
    return mol.beta_mol * mol.transmission * 1e6               # (m^-1 sr^-1 -> Mm^-1 sr^-1) x T^2


def wavelength_correct(beta, l2, lam, target, alpha, model):
    if not np.isfinite(lam) or abs(lam - target) < 1.0:
        return beta
    if model == "molaer":
        z_agl = l2["alt"] - l2["station_alt"]
        bml = _molecular_beta(z_agl, l2["station_alt"], lam)
        bmt = _molecular_beta(z_agl, l2["station_alt"], target)
        wl_corr = (lam / target) ** (-alpha)
        beta_aer = beta - bml[None, :]
        return bmt[None, :] + beta_aer / wl_corr
    wl_corr = (lam / target) ** (-alpha)
    return beta / wl_corr


# ---------------------------------------------------------- CAMS molecular (optimal wavelength conv.)
def find_cams_month(year, month):
    """Monthly CAMS file for (year,month): 0.4 deg primary (CAMS_04), 1 deg fallback (CAMS_1DEG)."""
    ds = f"{year}{month:02d}01"
    f = ensure_cams_file(CAMS_04, ds, auto_download=False)
    if f is None:
        f = ensure_cams_file(CAMS_1DEG, ds, auto_download=False)
    return f


# The CAMS water-vapour and T/p profiles depend ONLY on (file, lat, lon, month window), never on the
# instrument signal or the noise-filter mode — but the pipeline extracts them once per channel per
# WV-on/off twin (and again per noise-filter mode when a runner sweeps modes). Memoise them so a
# multi-mode / multi-channel sweep pays each monthly extraction once. Safe: the returned arrays are
# read-only downstream (interpolated, never mutated in place). Single-process (run_site); a ProcessPool
# worker keeps its own cache.
_CAMS_WV_CACHE = {}
_CAMS_TP_CACHE = {}


def _cams_wv_profile_cached(cams, lat, lon, tstart, tend):
    key = (str(cams), float(lat), float(lon), str(tstart), str(tend))
    if key not in _CAMS_WV_CACHE:
        _CAMS_WV_CACHE[key] = cams_water_vapor_profile(cams, lat, lon, tstart, tend)
    return _CAMS_WV_CACHE[key]


def _cams_tp_profile_cached(cams_file, lat, lon, tstart, tend):
    key = (str(cams_file), float(lat), float(lon), str(tstart), str(tend))
    if key not in _CAMS_TP_CACHE:
        _CAMS_TP_CACHE[key] = cams_temperature_pressure_profile(cams_file, lat, lon, tstart, tend)
    return _CAMS_TP_CACHE[key]


def _mol_att_cams(cams_file, lat, lon, tstart, tend, alt_asl, station_alt, wavelength_nm):
    """Molecular attenuated backscatter beta_mol*T2_mol [Mm^-1 sr^-1] on the instrument grid from the
    CAMS T/p profile (nearest grid point, night-mean). BELOW the lowest CAMS level (elevated / valley
    stations where the CAMS grid-cell surface sits above the real station -- Aosta ~ +1627 m) T and p
    are extrapolated HYDROSTATICALLY with the standard lapse rate, not held constant: constant-fill
    underestimates the air density (hence beta_mol) by up to ~20 % over Aosta's below-grid gates.
    Returns None if CAMS has no usable profile there."""
    prof = _cams_tp_profile_cached(cams_file, lat, lon, tstart, tend)
    if prof is None:
        return None
    h_cams, t_cams, p_cams = np.asarray(prof[0], "f8"), np.asarray(prof[1], "f8"), np.asarray(prof[2], "f8")
    if h_cams.size < 2 or not np.all(np.isfinite([h_cams[0], t_cams[0], p_cams[0]])):
        return None
    T = np.interp(alt_asl, h_cams, t_cams, left=np.nan, right=np.nan)
    P = np.interp(alt_asl, h_cams, p_cams, left=np.nan, right=np.nan)
    below = alt_asl < h_cams[0]
    if below.any():
        dz = h_cams[0] - alt_asl[below]                              # metres below the lowest level (>0)
        Tb = t_cams[0] + _LAPSE * dz                                 # standard lapse -> warmer downward
        T[below] = Tb
        P[below] = p_cams[0] * (Tb / t_cams[0]) ** (_G0 / (_RD * _LAPSE))   # barometric w/ lapse rate
    if not np.isfinite(T).any():
        return None
    mol = calculate_molecular_properties(T, P, alt_asl - station_alt, wavelength_nm * 1e-9)
    return mol.beta_mol * mol.transmission * 1e6


def wavelength_correct_molecular(beta, l2, lam, target, alpha):
    """Optimal 910->target conversion: WV-corrected beta in, analytic molecular (Rayleigh) from CAMS
    T/p removed, Angstrom (lam/target)^alpha applied to the AEROSOL RESIDUAL only, molecular re-added at
    the target line:  beta_target = beta_mol_att|target + [beta - beta_mol_att|lam] * (lam/target)^alpha.
    Per CAMS period (monthly 0.4 deg / 1 deg, else daily 0.4 deg); a period with no CAMS falls back to
    the US-std molecular."""
    if not np.isfinite(lam) or abs(lam - target) < 1.0:
        return beta
    alt_asl = np.asarray(l2["alt"], "f8")
    z_agl = alt_asl - l2["station_alt"]
    f = (lam / target) ** alpha
    out = beta.copy()
    for sel, cf, tstart, tend, _p in cams_periods(l2["time"]):
        bml = bmt = None
        if cf is not None:
            bml = _mol_att_cams(cf, l2["lat"], l2["lon"], tstart, tend, alt_asl, l2["station_alt"], lam)
            bmt = _mol_att_cams(cf, l2["lat"], l2["lon"], tstart, tend, alt_asl, l2["station_alt"], target)
        if bml is None or bmt is None:              # no CAMS this month -> US-std molecular
            bml = _molecular_beta(z_agl, l2["station_alt"], lam)
            bmt = _molecular_beta(z_agl, l2["station_alt"], target)
        out[sel, :] = bmt[None, :] + (beta[sel, :] - bml[None, :]) * f
    return out


# --------------------------------------------------------------------------- screening
def screen(beta, l2):
    b = beta.copy()
    b[l2["qf"] > 0] = np.nan
    cbh = l2["cbh"]
    has_cloud = np.any(np.isfinite(cbh) & (cbh > 0) & (cbh < 20000), axis=1)
    vv = l2["vv"].copy(); vv[vv < 0] = np.nan   # negative = "no fog" sentinel (L2 uses -1, CL61 L1 uses -99)
    has_excl = has_cloud | np.isfinite(vv)
    t = pd.to_datetime(l2["time"])
    if t.size > 1:
        dt = np.median(np.diff(t.values).astype("timedelta64[s]").astype(float))
        if dt > 0:
            win = int(2 * round(15 * 60 / dt) + 1)
            expanded = pd.Series(has_excl.astype(float)).rolling(win, center=True, min_periods=1).max().values > 0
        else:
            expanded = has_excl
    else:
        expanded = has_excl
    b[expanded, :] = np.nan
    return b


# --------------------------------------------------------------------------- gridding + stats
MIN_AVG_S = 1800.0    # minimum temporal coverage per averaging window [s] (>= 30 min requirement)
SNR_MIN = 3.0         # per-gate detection threshold over the averaging window (SNR3, cf.
                      # calibration/sensitivity: beta_min = SNR * sigma(tau))
_SNR_MIN_SAMPLES = 5  # need at least this many samples to estimate the noise at all


def snr_mask(X, snr_min=SNR_MIN):
    """Per-gate SNR of the window median: |median| / (robust sigma / sqrt(n_finite)).
    X is (time x range) over ONE averaging window. Returns a boolean keep-mask per gate
    (True = detected at >= snr_min). Windows too short to estimate noise keep everything.
    Uses the same robust scale (1.4826*MAD) as the operational sensitivity product."""
    from calibration.sensitivity.noise import robust_std
    X = np.asarray(X, "f8")
    if X.ndim != 2 or X.shape[0] < _SNR_MIN_SAMPLES:
        return np.ones(X.shape[-1], bool)
    with np.errstate(all="ignore"):
        med = np.nanmedian(X, axis=0)
        sig = robust_std(X, axis=0)
        nfin = np.isfinite(X).sum(axis=0)
        snr = med / (sig / np.sqrt(np.maximum(nfin, 1)))
    # keep where the estimate is undefined (sig==0 with finite med) or passes the threshold
    keep = ~np.isfinite(snr) | (snr >= snr_min)
    keep[~np.isfinite(med)] = True   # NaN gates stay NaN anyway; do not turn them into "removed"
    return keep


# --- Cloudnet / ceilopyter per-profile noise filter (the "cloudnet" mode) --------------------------
CLOUDNET_TOP_FRACTION = 0.10   # fraction of far-range gates used to estimate the per-profile noise


def cloudnet_noise_mask(P, snr_limit=None, noise_floor=None):
    """Boolean is_noise mask (time x range) — the Cloudnet filter of ceilopyter.noise.remove_noise.
    P is the NON-range-corrected signal rcs_0/r^2. Per profile the noise is the std over the top
    CLOUDNET_TOP_FRACTION (far-range) gates; a sample is noise where P < snr_limit*max(noise, floor).

    Reimplemented rather than imported to (a) fix ceilopyter's empty top slice when there are no
    trailing all-zero gates — E-PROFILE rcs_0 has none, so the upstream function would silently no-op
    — and (b) tolerate NaN (nanstd; NaN samples stay unflagged, i.e. kept as NaN). With the default
    floor 0 the far-range std alone sets a scale-free per-profile SNR, so the same code serves the
    three instruments whose rcs_0 magnitudes span ~1e12 (CHM ~9e5, CL31 ~2e2, CL61 ~7e-7)."""
    snr_limit = (float(os.environ.get("ALC_VAL_CLOUDNET_SNR", "5")) if snr_limit is None
                 else snr_limit)
    noise_floor = (float(os.environ.get("ALC_VAL_CLOUDNET_FLOOR", "0")) if noise_floor is None
                   else noise_floor)
    P = np.asarray(P, "f8")
    if P.ndim != 2 or P.shape[1] == 0:
        return np.zeros(P.shape, bool)
    # Match ceilopyter: drop trailing all-zero range gates from the noise reference. (Its own slice
    # beta[:, -n_top-n_zeros:-n_zeros] becomes EMPTY when there are none — the E-PROFILE rcs_0 case —
    # so it no-ops; here the reference then falls back to the last n_top gates.)
    zero_ranges = np.all(P == 0, axis=0)
    n_zeros = int(np.argmax(~zero_ranges[::-1]))
    n_top = max(1, int(round(P.shape[1] * CLOUDNET_TOP_FRACTION)))
    end = P.shape[1] - n_zeros
    with np.errstate(all="ignore"):
        noise = np.nanstd(P[:, max(0, end - n_top):end], axis=1)   # per-profile far-range noise
    noise = np.where(np.isfinite(noise), noise, 0.0)
    noise = np.maximum(noise, noise_floor)
    with np.errstate(all="ignore"):
        snr = P / noise[:, None]                              # noise==0 -> +/-inf -> keep(+)/flag(-)
    return snr < snr_limit                                    # NaN -> False -> kept as NaN (faithful)


def apply_cloudnet_filter(rcs, rng):
    """Return rcs_0 with the Cloudnet-flagged noise samples set to NaN. The filter runs on the
    non-range-corrected signal rcs_0/r^2 (r in km); gates with r<=0 (the CL61 first gate at 0 m) are
    excluded from the estimate and never flagged."""
    r2 = (np.asarray(rng, "f8") * 1e-3) ** 2
    ok = r2 > 0
    P = np.full(rcs.shape, np.nan)
    P[:, ok] = rcs[:, ok] / r2[ok][None, :]
    out = rcs.copy()
    out[cloudnet_noise_mask(P)] = np.nan
    return out


def retime_hourly(time, arrays, min_cov_s=None, snr_idx=()):
    """Median-aggregate each (time x range) or (time,) array onto a regular 60-min grid.
    min_cov_s: bins whose samples span less than this coverage [s] are left NaN (enforces the
    >= 30-min temporal averaging). snr_idx: indices of 2D arrays whose gates with window
    SNR < SNR_MIN are removed (NaN) — the science stream, not the display stream."""
    t = pd.to_datetime(time)
    idx = pd.DatetimeIndex(t)
    grid = pd.date_range(idx.min().floor("h"), idx.max().ceil("h"), freq="60min")
    binid = np.clip(np.searchsorted(grid.values, idx.values, side="right") - 1, 0, len(grid) - 1)
    dt = float(np.median(np.diff(idx.values).astype("timedelta64[s]").astype(float))) if t.size > 1 else np.nan
    out = []
    for j, A in enumerate(arrays):
        if A.ndim == 1:
            G = np.full(len(grid), np.nan)
            for b in range(len(grid)):
                sel = binid == b
                if sel.any():
                    if min_cov_s and np.isfinite(dt) and sel.sum() * dt < min_cov_s:
                        continue
                    G[b] = np.nanmedian(A[sel])
        else:
            G = np.full((len(grid), A.shape[1]), np.nan)
            for b in range(len(grid)):
                sel = binid == b
                if sel.any():
                    if min_cov_s and np.isfinite(dt) and sel.sum() * dt < min_cov_s:
                        continue
                    with np.errstate(all="ignore"):
                        G[b] = np.nanmedian(A[sel], axis=0)
                    if j in snr_idx:
                        G[b, ~snr_mask(A[sel])] = np.nan
        out.append(G)
    return grid.values, out


def build_common_grid(alts):
    z0, z1, dz = -np.inf, np.inf, 0.0
    for a in alts:
        d = np.median(np.diff(a))
        dz = max(dz, d); z0 = max(z0, np.min(a)); z1 = min(z1, np.max(a))
    dz = max(round(dz), 1)
    return np.arange(z0, z1 + dz, dz)


def regrid(beta, alt_src, altGrid):
    half = np.median(np.diff(altGrid)) / 2
    G = np.full((beta.shape[0], altGrid.size), np.nan)
    for iz, z in enumerate(altGrid):
        sel = np.abs(alt_src - z) < half
        if sel.any():
            with np.errstate(all="ignore"):
                G[:, iz] = np.nanmean(beta[:, sel], axis=1)
    return G


def _stats(cur, ref, zmask):
    """Pairwise metrics over the altitude band. Linear metrics (bias, relbias, r) are kept for
    continuity; the log-space r and the median relative bias are robust to the log-distributed
    beta (molecular floor + rare large aerosol/cloud values dominate the linear moments)."""
    a = cur[:, zmask]; b = ref[:, zmask]
    m = np.isfinite(a) & np.isfinite(b)
    a = a[m]; b = b[m]; n = a.size
    if n <= 2:
        return dict(n=int(n), bias=np.nan, medbias=np.nan, rmse=np.nan, std=np.nan,
                    relbias_pct=np.nan, medrelbias_pct=np.nan, r=np.nan, r_log=np.nan, n_log=0)
    d = a - b
    mb = float(np.mean(b))
    pos = (a > 0) & (b > 0)
    n_log = int(pos.sum())
    r_log = float(np.corrcoef(np.log10(a[pos]), np.log10(b[pos]))[0, 1]) if n_log > 2 else np.nan
    bpos = b > 0
    medrel = float(100 * np.median(d[bpos] / b[bpos])) if bpos.any() else np.nan
    return dict(n=int(n), bias=float(np.mean(d)), medbias=float(np.median(d)),
                rmse=float(np.sqrt(np.mean(d**2))), std=float(np.std(d, ddof=1)),
                relbias_pct=float(100 * np.mean(d) / mb) if mb != 0 else np.nan,
                medrelbias_pct=medrel,
                r=float(np.corrcoef(a, b)[0, 1]),
                r_log=r_log, n_log=n_log)


# --------------------------------------------------------------------------- main process
def process(cfg):
    """cfg: dict(wmo, start, end, referenceChannel(0-based), channels[list], lambda_target, alpha, zMin, zMax).
    Each channel: dict(wmo, ident, calib, label, key, wavelengthModel?)."""
    target = cfg.get("lambda_target", 1064.0); alpha = cfg.get("alpha", 1.0)
    zmin = cfg.get("zMin", 500.0); zmax = cfg.get("zMax", 3000.0)
    chans = []
    level_l1 = cfg.get("dataLevel") == "L1"
    for ch in cfg["channels"]:
        l2 = (read_l1(ch["wmo"], ch["ident"], cfg["start"], cfg["end"]) if level_l1
              else read_l2(ch["wmo"], ch["ident"], cfg["start"], cfg["end"]))
        if l2 is None:
            chans.append(None); continue
        beta = l2["beta"].copy()
        # optional per-channel hood-offset removal: subtract the measured b_dark(z) from L1 rcs_0
        # (physical model P^bgi x z^2) BEFORE the calibration scaling. Used by the CL61/CHM15k/CL31
        # offset-corr channels. bdark_scale converts the stored b_phys units to the rcs_0 units
        # (CL61: 1e-6 Mm->SI; CHM15k & CL31: 1.0, already rcs_0 units).
        if ch.get("bdark") is not None:
            _bd = np.load(ch["bdark"])
            _ragl = l2["alt"] - l2.get("station_alt", 0.0)
            _sc = ch.get("bdark_scale", 1.0)
            beta = beta - _sc * np.interp(_ragl, _bd["rng"], _bd["b_phys"], left=0.0, right=0.0)[None, :]
        if level_l1:
            # APPLY the calibration to the L1 raw signal: beta_att [Mm^-1 sr^-1] = rcs_0 / C_L * 1e6.
            # Native channels use the OPERATIONAL calout Kalman C_L (same as the dashboard/balfrin);
            # offset-corr Rayleigh channels use their recalibrated series via calib_key (same scale,
            # verified ratio 1.000). Both methods store C_L (a true lidar constant).
            if ch.get("calib_key"):
                cal = load_calib_series(ch["calib_key"], "L1")
            else:
                cal = load_calout_kalman(ch["wmo"], ch["ident"], ch["calib"]) if ch["calib"] != "none" else None
            if ch["calib"] != "none" and cal is None:
                print(f"    [skip] {ch['label']}: no calibration ({ch.get('calib_key') or ch['wmo']+'_'+ch['ident']+' '+ch['calib']})")
                chans.append(None); continue
            if cal is not None:
                ck = interp_calib(cal[0], cal[1], l2["time"])
                beta = beta / ck[:, None] * 1e6
                med_corr = float(np.nanmedian(1e6 / ck))
            else:
                med_corr = np.nan
        else:
            # L2 path: re-scale the L2 product. The calibration CSVs hold the absolute lidar
            # constant C_L for BOTH methods (single physical constant everywhere).
            #   Rayleigh: corr = calibration_constant_0 / C_L  (undo the provider constant)
            #   Cloud:    corr = INSTRUMENT_CAL_DEFAULT / C_L  (the O'Connor multiplier; the
            #             calout C_L is expressed vs the default applied in the L1 calibration)
            cal = load_calib_series(ch["key"], cfg.get("calibLevel"))
            if ch["calib"] != "none" and cal is None:
                print(f"    [skip] {ch['label']}: no calibration series ({ch['key']})")
                chans.append(None); continue
            if cal is not None and ch["calib"] != "none":
                ck = interp_calib(cal[0], cal[1], l2["time"])
                if ch["calib"] == "rayleigh" and ch.get("itype") != "CL61":
                    corr = l2["calc"] / ck
                else:
                    # cloud (all types) and CL61 (both methods): the calout C_L is defined against
                    # the L1 rcs_0 with the applied default; the CL61 L2 calibration_constant_0 is
                    # Vaisala's internal factor in a DIFFERENT unit system and must not be used
                    # (validated against the convention-free L1 path: CL61-Rayleigh reads +13 %
                    # there, while calc/C_L produced a spurious -43 %).
                    corr = INSTRUMENT_CAL_DEFAULT.get(ch.get("itype", ""), 1.0) / ck
                beta = beta * corr[:, None]
                med_corr = float(np.nanmedian(corr))
            else:
                med_corr = np.nan
        # WV (910 nm)
        wl = l2["wavelength"] if np.isfinite(l2["wavelength"]) else WV_PARAMS.get(ch.get("itype", ""), (np.nan,))[0]
        if np.isfinite(wl) and in_water_vapor_band(wl):
            lam0, fwhm = WV_PARAMS.get(ch.get("itype", ""), (910.0, 3.4))
            beta, wvinfo = apply_wv(beta, l2, lam0, fwhm)
            if wvinfo["months_excluded"]:
                print(f"    [wv] {ch['label']}: no CAMS for {','.join(wvinfo['months_excluded'])}"
                      " -> month(s) NaN-masked", flush=True)
        # wavelength normalisation
        beta = wavelength_correct(beta, l2, wl if np.isfinite(wl) else target, target, alpha,
                                  ch.get("wavelengthModel", "angstrom"))
        # streams
        beta_disp = beta.copy(); beta_disp[l2["qf"] > 0] = np.nan
        beta_scr = screen(beta, l2)
        cbh_low = np.nanmin(np.where((l2["cbh"] > 0) & (l2["cbh"] < 20000), l2["cbh"], np.nan), axis=1) \
            if l2["cbh"].ndim == 2 else l2["cbh"]
        chans.append(dict(ch=ch, l2=l2, beta_scr=beta_scr, beta_disp=beta_disp, cbh=cbh_low, med_corr=med_corr))

    valid = [c for c in chans if c is not None]
    if not valid:
        return None
    # temporal sync: per-channel hourly median, then union grid. The science stream (scr) gets
    # the >=30-min coverage requirement and the per-gate SNR>=3 removal; the display stream keeps
    # everything visible.
    gridded = []
    for c in valid:
        g, arrs = retime_hourly(c["l2"]["time"], [c["beta_scr"], c["beta_disp"], c["cbh"]],
                                min_cov_s=MIN_AVG_S, snr_idx=(0,))
        gridded.append(dict(c=c, grid=g, scr=arrs[0], disp=arrs[1], cbh=arrs[2]))
    union = np.unique(np.concatenate([g["grid"] for g in gridded]))
    # reindex each to union time grid
    for g in gridded:
        idx = {t: i for i, t in enumerate(g["grid"])}
        pos = np.array([idx.get(t, -1) for t in union])
        def take(A):
            if A.ndim == 1:
                out = np.full(union.size, np.nan)
            else:
                out = np.full((union.size, A.shape[1]), np.nan)
            ok = pos >= 0
            out[ok] = A[pos[ok]]
            return out
        g["scrU"] = take(g["scr"]); g["dispU"] = take(g["disp"]); g["cbhU"] = take(g["cbh"])
    # common altitude grid
    altGrid = build_common_grid([g["c"]["l2"]["alt"] for g in gridded])
    for g in gridded:
        g["betaC"] = regrid(g["scrU"], g["c"]["l2"]["alt"], altGrid)
        g["dispC"] = regrid(g["dispU"], g["c"]["l2"]["alt"], altGrid)
    # stats vs reference
    station_alt = valid[0]["l2"]["station_alt"]
    zmask = (altGrid >= zmin + station_alt) & (altGrid <= zmax + station_alt)
    ref_ch = cfg["channels"][cfg.get("referenceChannel", 0)]
    iref = next((i for i, g in enumerate(gridded) if g["c"]["ch"] is ref_ch), 0)
    ref = gridded[iref]["betaC"]
    R = dict(altGrid=altGrid, time_sync=union, station=dict(altitude=station_alt,
             lat=valid[0]["l2"]["lat"], lon=valid[0]["l2"]["lon"]), channels=[], beta=[], beta_disp=[], cbh=[], stats=[])
    for g in gridded:
        R["channels"].append(dict(label=g["c"]["ch"]["label"], calib=g["c"]["ch"]["calib"],
                                  itype=g["c"]["ch"].get("itype", ""),
                                  wavelength=g["c"]["l2"]["wavelength"], med_corr=g["c"]["med_corr"]))
        R["beta"].append(g["betaC"]); R["beta_disp"].append(g["dispC"]); R["cbh"].append(g["cbhU"])
        R["stats"].append(_stats(g["betaC"], ref, zmask))
    return R
