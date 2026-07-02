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

from calibration.io.cams import ensure_cams_file
from calibration.cloud.calibration import INSTRUMENT_CAL_DEFAULT
from calibration.water_vapor_correction.water_vapor import (
    cams_water_vapor_profile, two_way_wv_transmission, laser_spectrum_for, in_water_vapor_band,
    cams_point_too_far)
from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties
from calibration.config import InstrumentType

REPO = Path(__file__).resolve().parents[2]
L2_MONTHLY = Path("A:/E-PROFILE_L2_monthly")
CALIB = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/calib")
CAMS = "D:/CAMS"
# --- L1 data path (validate by applying the CSCS calibration to the raw range-corrected signal) ---
# rcs_0 is read from the native daily L1 archive and turned into attenuated backscatter with
# beta_att [Mm^-1 sr^-1] = rcs_0 / C_L * 1e6, where C_L is the CSCS Kalman lidar constant from the
# calout <key>_kalman.csv. (No L2 product, no raw vendor files.) Both env-overridable.
L1_ROOT = Path(os.environ.get("ALC_VAL_L1_ROOT", "D:/E-PROFILE_L1_2026"))
CALOUT = Path(os.environ.get(
    "ALC_VAL_CALOUT", "C:/DATA/Projects/202606_E-PROFILE_calibration/E_PROFILE_calout_2025_2026"))
WV_LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
STD_ATM = REPO / "calibration" / "data" / "standard_atmosphere_US_1976_50km.csv"
WV_PARAMS = {"CL31": (909.7, 6.0), "CL51": (910.0, 3.4), "CL61": (910.74, 1.0),
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
    """Read ONE daily L1 file and median-retime rcs_0 / cloud_base_height / vertical_visibility to its
    hourly grid (keeps memory bounded over 2 years). Returns (grid, rcs_h, cbh_h, vv_h, range) or None."""
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
    except Exception:
        return None
    grid, (rcsh, cbhh, vvh) = retime_hourly(t, [rcs, cbh, vv], min_cov_s=MIN_AVG_S, snr_idx=(0,))
    return grid, rcsh, cbhh, vvh, rng


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
            if res is None or res[4].size != nR:   # skip unreadable / different range grid
                continue
            parts.append(res[:4])
    if not parts:
        return None
    time = np.concatenate([p[0] for p in parts])
    rcs = np.concatenate([p[1] for p in parts], axis=0)
    cbh = np.concatenate([p[2] for p in parts], axis=0)
    vv = np.concatenate([p[3] for p in parts])
    o = np.argsort(time)
    return dict(time=time[o], alt=meta["rng"] + meta["salt"], beta=rcs[o], calc=None,
                qf=np.zeros(rcs.shape, dtype="f8")[o], cbh=cbh[o], vv=vv[o],
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
def apply_wv(beta, l2, lam0, fwhm):
    """Divide beta by the two-way WV transmission, per month (exact-month CAMS; missing month -> NaN)."""
    t = l2["time"]
    z_asl = l2["alt"]
    z_agl = z_asl - l2["station_alt"]
    months = pd.PeriodIndex(pd.to_datetime(t), freq="M")
    out = beta.copy()
    info = {"months": [], "months_excluded": [], "median_t2": []}
    for p in months.unique():
        sel = np.asarray(months == p)
        ds0 = f"{p.year}{p.month:02d}01"
        cams = ensure_cams_file(Path(CAMS), ds0, auto_download=False)
        if cams is None:
            out[sel, :] = np.nan
            info["months_excluded"].append(str(p)); continue
        # same guard as the calibration (flag -10 there): NEVER correct with a far-away
        # domain-edge grid point — exclude the data instead of reporting it uncorrected.
        if cams_point_too_far(cams, l2["lat"], l2["lon"]):
            out[sel, :] = np.nan
            info["months_excluded"].append(f"{p} (CAMS too far)"); continue
        tstart = np.datetime64(f"{p.year}-{p.month:02d}-01") - np.timedelta64(1, "D")
        tend = (np.datetime64(f"{p.year}-{p.month:02d}-01") + np.timedelta64(40, "D"))
        prof = cams_water_vapor_profile(cams, l2["lat"], l2["lon"], tstart, tend)
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
    standard atmosphere. The two-way molecular transmission matters: the measured signal the
    molaer model subtracts from is attenuated, and at 532 nm T^2_mol is already ~0.85-0.90 by
    2-3 km — subtracting the UNattenuated beta_mol there biases the extracted aerosol low by
    ~beta_mol*(1-T^2), which after the 532->1064 recombination produced a spurious ~-40 % on the
    Mini-MPL (the native 532 nm EARLINET comparison shows the instrument itself is within a few %)."""
    grid = np.arange(0, 15001, 30.0)
    atm = load_standard_atmosphere(STD_ATM, grid)
    mol = calculate_molecular_properties(atm.temperature, atm.pressure, grid, wavelength_nm * 1e-9)
    bmol_att = mol.beta_mol * mol.transmission * 1e6   # (m^-1 sr^-1 -> Mm^-1 sr^-1) x two-way T^2
    return np.interp(z_agl, grid, bmol_att, left=np.nan, right=np.nan)


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
        if level_l1:
            # APPLY the CSCS calibration to the raw signal: beta_att [Mm^-1 sr^-1] = rcs_0 / C_L * 1e6,
            # C_L = the calout Kalman lidar constant (same physical constant for Rayleigh and cloud).
            cal = load_calout_kalman(ch["wmo"], ch["ident"], ch["calib"]) if ch["calib"] != "none" else None
            if ch["calib"] != "none" and cal is None:
                print(f"    [skip] {ch['label']}: no calout Kalman ({ch['wmo']}_{ch['ident']} {ch['calib']})")
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
