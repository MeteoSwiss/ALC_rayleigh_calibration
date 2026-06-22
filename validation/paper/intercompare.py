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
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from calibration.io.cams import ensure_cams_file
from calibration.water_vapor_correction.water_vapor import (
    cams_water_vapor_profile, two_way_wv_transmission, laser_spectrum_for, in_water_vapor_band)
from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties
from calibration.config import InstrumentType

REPO = Path(__file__).resolve().parents[2]
L2_MONTHLY = Path("A:/E-PROFILE_L2_monthly")
CALIB = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/calib")
CAMS = "D:/CAMS"
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
    <key>_<level>.csv, fall back to <key>.csv (the un-suffixed/legacy series)."""
    f = None
    for cand in ([CALIB / f"{key}_{level}.csv"] if level else []) + [CALIB / f"{key}.csv"]:
        if cand.exists():
            f = cand
            break
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
_STD = None
def _molecular_beta(z_agl, station_alt, wavelength_nm):
    """Molecular attenuated-ish backscatter [Mm^-1 sr^-1] on z_agl from the US standard atmosphere."""
    global _STD
    if _STD is None:
        _STD = load_standard_atmosphere(STD_ATM, np.arange(0, 15001, 30.0))
    grid = np.arange(0, 15001, 30.0)
    atm = load_standard_atmosphere(STD_ATM, grid)
    mol = calculate_molecular_properties(atm.temperature, atm.pressure, grid, wavelength_nm * 1e-9)
    bmol_grid = mol.beta_mol * 1e6  # m^-1 sr^-1 -> Mm^-1 sr^-1
    return np.interp(z_agl, grid, bmol_grid, left=np.nan, right=np.nan)


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
    vv = l2["vv"].copy(); vv[vv == -1] = np.nan
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
def retime_hourly(time, arrays):
    """Median-aggregate each (time x range) or (time,) array onto a regular 60-min grid."""
    t = pd.to_datetime(time)
    idx = pd.DatetimeIndex(t)
    grid = pd.date_range(idx.min().floor("h"), idx.max().ceil("h"), freq="60min")
    binid = np.clip(np.searchsorted(grid.values, idx.values, side="right") - 1, 0, len(grid) - 1)
    out = []
    for A in arrays:
        if A.ndim == 1:
            G = np.full(len(grid), np.nan)
            for b in range(len(grid)):
                sel = binid == b
                if sel.any():
                    G[b] = np.nanmedian(A[sel])
        else:
            G = np.full((len(grid), A.shape[1]), np.nan)
            for b in range(len(grid)):
                sel = binid == b
                if sel.any():
                    with np.errstate(all="ignore"):
                        G[b] = np.nanmedian(A[sel], axis=0)
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
    a = cur[:, zmask]; b = ref[:, zmask]
    m = np.isfinite(a) & np.isfinite(b)
    a = a[m]; b = b[m]; n = a.size
    if n <= 2:
        return dict(n=int(n), bias=np.nan, medbias=np.nan, rmse=np.nan, std=np.nan, relbias_pct=np.nan, r=np.nan)
    d = a - b
    return dict(n=int(n), bias=float(np.mean(d)), medbias=float(np.median(d)),
                rmse=float(np.sqrt(np.mean(d**2))), std=float(np.std(d, ddof=1)),
                relbias_pct=float(100 * np.mean(d) / np.mean(b)),
                r=float(np.corrcoef(a, b)[0, 1]))


# --------------------------------------------------------------------------- main process
def process(cfg):
    """cfg: dict(wmo, start, end, referenceChannel(0-based), channels[list], lambda_target, alpha, zMin, zMax).
    Each channel: dict(wmo, ident, calib, label, key, wavelengthModel?)."""
    target = cfg.get("lambda_target", 1064.0); alpha = cfg.get("alpha", 1.0)
    zmin = cfg.get("zMin", 500.0); zmax = cfg.get("zMax", 3000.0)
    chans = []
    for ch in cfg["channels"]:
        l2 = read_l2(ch["wmo"], ch["ident"], cfg["start"], cfg["end"])
        if l2 is None:
            chans.append(None); continue
        beta = l2["beta"].copy()
        # calibration (from the L1- or L2-derived Kalman series, per cfg['calibLevel']).
        # Only RAYLEIGH swaps L1<->L2: its lidar constant is on the same scale for both sources
        # (binned L1 == L2 to ~1%). The cloud O'Connor coefficient is on the INPUT scale (raw
        # rcs_0 V*m^2 for L1 vs the attbsc_0 Mm^-1 sr^-1 for L2), so it is NOT transferable onto
        # the L2 beta - cloud always uses the L2 coefficient here.
        clevel = cfg.get("calibLevel") if ch["calib"] == "rayleigh" else "L2"
        cal = load_calib_series(ch["key"], clevel)
        if ch["calib"] != "none" and cal is None:
            print(f"    [skip] {ch['label']}: no calibration series ({ch['key']} / {clevel})")
            chans.append(None); continue
        if cal is not None and ch["calib"] != "none":
            ck = interp_calib(cal[0], cal[1], l2["time"])
            # Rayleigh: correction = calibration_constant_0 / C_L (divide by the stored lidar constant).
            # Cloud: the O'Connor coefficient C calibrates beta_true = C * beta_obs. The cloud reader
            # integrates the L2 attbsc_0 on its stored scale (units "1e-6*1/(m*sr)" = Mm^-1 sr^-1, i.e.
            # beta_factor 1), so the returned C is 1e6x too small relative to the physical 1/(m*sr)
            # integral; multiplying by 1e6 restores the physical multiplier on attbsc_0 [Mm^-1 sr^-1].
            # Uniform across instruments (CL31/CL51/CL61 L2 share the units string); calibration_constant_0
            # does NOT enter (the cloud method uses attbsc_0, not the reconstructed RCS).
            corr = (l2["calc"] / ck) if ch["calib"] == "rayleigh" else (ck * 1e6)
            beta = beta * corr[:, None]
            med_corr = float(np.nanmedian(corr))
        else:
            med_corr = np.nan
        # WV (910 nm)
        wl = l2["wavelength"] if np.isfinite(l2["wavelength"]) else WV_PARAMS.get(ch.get("itype", ""), (np.nan,))[0]
        if np.isfinite(wl) and in_water_vapor_band(wl):
            lam0, fwhm = WV_PARAMS.get(ch.get("itype", ""), (910.0, 3.4))
            beta, _ = apply_wv(beta, l2, lam0, fwhm)
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
    # temporal sync: per-channel hourly median, then union grid
    gridded = []
    for c in valid:
        g, arrs = retime_hourly(c["l2"]["time"], [c["beta_scr"], c["beta_disp"], c["cbh"]])
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
                                  wavelength=g["c"]["l2"]["wavelength"], med_corr=g["c"]["med_corr"]))
        R["beta"].append(g["betaC"]); R["beta_disp"].append(g["dispC"]); R["cbh"].append(g["cbhU"])
        R["stats"].append(_stats(g["betaC"], ref, zmask))
    return R
