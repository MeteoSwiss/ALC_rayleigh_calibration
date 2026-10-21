# -*- coding: utf-8 -*-
"""Is the residual CHM15k annual gain cycle driven by sunlight / day length? (report 19, q16 agent b)

After the kappa x T_int normalisation the CHM15k Rayleigh constant C_L keeps an annual cycle of ~+-9 %,
height-independent, common to co-located CHM15k, with a peak near the solstice (DOY ~150-177) while the
calibration-night T_int peaks at DOY ~206. Hypothesis tested here: the residual follows sunlight (day length,
twilight proximity of the calibration window, insolation of the preceding day) or the daytime solar heating of
the instrument (with a memory of hours to weeks), not the night temperature.

Per calibration night of every CHM15k unit (main optical module, the q15f nightly table) we compute:
  * geometry: day length, hours from sunset to the window start / centre, hours from the window end to sunrise,
    mean and minimum solar zenith angle in the window, night length;
  * insolation of the preceding 24 h: top-of-atmosphere, clear-sky surface (Haurwitz) and a cloud-aware proxy
    (clear-sky x Kasten-Czeplak factor with the ceilometer's own clear-profile fraction); exponentially
    weighted (EW) means with time constants 1-21 d; measured global radiation for Bern (SwissMetNet OGD);
  * daytime instrument heating: max T_int of the preceding daylight, day max - night T_int, T_int - T_ext,
    hours with temp_lom out of regulation (>= +1.5 K above its night set point) and with the detector off its
    set point (> 0.5 K), the same as EW means over 6 h - 21 d; EW means of T_int (thermal memory).
Tests:
  1. phase: annual harmonic of each driver vs the residual r = ln C - kappa (T_int - median), with the intra-night
     kappa (q2a b_T_ctrl) and the multi-day kappa (q1 b2);
  2. synoptic: 31-day anomalies, ln C on T_int + driver per unit (partial correlation), pooled with a bootstrap
     over units; transfer of the synoptic gain to the annual scale; latitude contrast (day-length amplitude grows
     x2.5 from 40 to 60 N) and twilight-proximity contrast;
  3. joint model ln C = a + kappa T_int + gamma D per unit era (kappa fixed and free), its effect on the annual
     amplitude, and an out-of-sample test (fit 2025, test 2026) per unit and with one network gamma; co-located
     pairs (differential drivers vs the A/B ratio);
  4. q15e open point: synoptic kappa at 7-8 km (~0) vs in the window (+2 %/K) - decomposition of the T_int
     anomaly into its ambient (T_ext) and internal parts and control for sunshine / aerosol.

Run:  .venv/Scripts/python.exe validation/chm15k_stability/q16b_insolation.py [--rebuild]
Outputs: C:/DATA/Projects/202610_CHM15k_stability/q16/q16b_*
"""
from __future__ import annotations

import argparse
import csv
import itertools
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.signal import lfilter
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
from common import BLOCKS_PER_DAY, OUT, PAIR_SITES, load_cache, load_constants  # noqa: E402
from q1_temperature import ANOM_WINDOW, running_median  # noqa: E402
from validation.paper.kalman_improve_alc import flag_outliers  # noqa: E402

warnings.filterwarnings("ignore")
Q15 = OUT / "q15"
Q16 = OUT / "q16"
Q16.mkdir(exist_ok=True)
W = 2 * np.pi / 365.25
PAYERNE = "0-20000-0-06610_A"
N_WORKERS = 4
SPLIT = 20454.0                      # 2026-01-01 (days since 1970): fit 2025, test 2026
SUNSET_COSZ = np.cos(np.deg2rad(90.833))
LOM_EXC = 1.5                        # K above the night set point = out of regulation (whole-degree reporting)
DET_EXC = 0.5                        # K off the detector set point
TAUS_LOM = (0.25, 1.0, 3.0, 7.0, 21.0)
TAUS_T = (1.0, 3.0, 7.0, 21.0)
TAUS_GHI = (3.0, 7.0, 21.0)
SMN = {"0-20008-0-BRN_A": ("ber", ["historical_2020-2029", "recent", "now"]),
       "0-756-2-ERISWIL_A": ("nap", ["recent", "now"])}
SMN_DIR = Path("D:/OGD_SMN")

# Drivers, their family and unit. Order = table order.
DRIVERS = {
    "daylen_h": ("geometry", "h"),
    "night_len_h": ("geometry", "h"),
    "sunset_to_t0_h": ("geometry (twilight)", "h"),
    "sunset_to_tmid_h": ("geometry (twilight)", "h"),
    "t1_to_sunrise_h": ("geometry (twilight)", "h"),
    "sza_win_mean": ("geometry (twilight)", "deg"),
    "sza_win_min": ("geometry (twilight)", "deg"),
    "toa_prev": ("insolation", "MJ/m2/d"),
    "csky_prev": ("insolation", "MJ/m2/d"),
    "ghi_proxy_prev": ("insolation", "MJ/m2/d"),
    "cf_day": ("insolation (sunshine)", "frac"),
    "ghi_proxy_ew3d": ("insolation", "W/m2"),
    "ghi_proxy_ew7d": ("insolation", "W/m2"),
    "ghi_proxy_ew21d": ("insolation", "W/m2"),
    "smn_gre_prev": ("insolation (measured)", "MJ/m2/d"),
    "Tint_daymax": ("instrument heating", "degC"),
    "dT_day_night": ("instrument heating", "K"),
    "Tint_Text_night": ("instrument heating", "K"),
    "Tint_Text_daymax": ("instrument heating", "K"),
    "lom_hours_prev": ("instrument heating (LOM)", "h"),
    "lom_maxexc_prev": ("instrument heating (LOM)", "K"),
    "det_hours_prev": ("instrument heating (detector)", "h"),
    "lom_ew0.25d": ("LOM memory", "h/d"),
    "lom_ew1d": ("LOM memory", "h/d"),
    "lom_ew3d": ("LOM memory", "h/d"),
    "lom_ew7d": ("LOM memory", "h/d"),
    "lom_ew21d": ("LOM memory", "h/d"),
    "Tint_ew1d": ("thermal memory", "degC"),
    "Tint_ew3d": ("thermal memory", "degC"),
    "Tint_ew7d": ("thermal memory", "degC"),
    "Tint_ew21d": ("thermal memory", "degC"),
}
SYNOPTIC_NONE = ("daylen_h", "night_len_h", "toa_prev", "csky_prev")   # no synoptic variance


# ---------------------------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------------------------
def harmonic(t, y):
    """Fit y = m + a cos(w doy) + b sin(w doy). Returns (phasor a+ib, residuals, fitted)."""
    ok = np.isfinite(t) & np.isfinite(y)
    res = np.full(len(y), np.nan)
    if ok.sum() < 30:
        return np.nan + 0j, res, res.copy()
    d = np.mod(t[ok], 365.25)
    X = np.column_stack([np.ones(ok.sum()), np.cos(W * d), np.sin(W * d)])
    c, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    res[ok] = y[ok] - X @ c
    da = np.mod(t, 365.25)
    return c[1] + 1j * c[2], res, c[0] + c[1] * np.cos(W * da) + c[2] * np.sin(W * da)


def peak_doy(p):
    return float(np.mod(np.angle(p) / W, 365.25)) if np.isfinite(p) else np.nan


def wrap_doy(d):
    return (d + 182.625) % 365.25 - 182.625


def anomaly(t, x):
    return x - running_median(t, x, ANOM_WINDOW)


def ols(X, y):
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ c
    dof = max(len(y) - X.shape[1], 1)
    try:
        se = np.sqrt(np.diag(res @ res / dof * np.linalg.inv(X.T @ X)))
    except np.linalg.LinAlgError:
        se = np.full(X.shape[1], np.nan)
    return c, se, res


def circ_median_doy(d):
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    if d.size == 0:
        return np.nan
    m = np.angle(np.mean(np.exp(1j * W * d))) / W
    return float(np.mod(m + np.median(wrap_doy(d - m)), 365.25))


# ---------------------------------------------------------------------------------------------
# Solar geometry (NOAA low-precision almanac, ~0.01 deg; enough for SZA and sunset times)
# ---------------------------------------------------------------------------------------------
def solar(t, lat, lon):
    """cos(SZA), declination [rad] and the Sun-Earth distance factor (1/r^2) for UTC days since 1970."""
    n = np.asarray(t, float) + 2440587.5 - 2451545.0
    L = np.mod(280.460 + 0.9856474 * n, 360.0)
    g = np.deg2rad(np.mod(357.528 + 0.9856003 * n, 360.0))
    lam = np.deg2rad(L + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    eps = np.deg2rad(23.439 - 4e-7 * n)
    dec = np.arcsin(np.sin(eps) * np.sin(lam))
    ra = np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam))
    gmst = np.deg2rad(np.mod(280.46061837 + 360.98564736629 * n, 360.0))
    H = gmst + np.deg2rad(lon) - ra
    phi = np.deg2rad(lat)
    cosz = np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.cos(H)
    rfac = 1.0 / (1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)) ** 2
    return cosz, dec, rfac


def day_length(t, lat):
    """Hours with the Sun's upper limb above the horizon (refraction 0.833 deg) on the day of t."""
    _, dec, _ = solar(t, lat, 0.0)
    phi = np.deg2rad(lat)
    x = (np.sin(np.deg2rad(-0.833)) - np.sin(phi) * np.sin(dec)) / (np.cos(phi) * np.cos(dec))
    return 2.0 * np.rad2deg(np.arccos(np.clip(x, -1, 1))) / 15.0


def ghi_clear(cosz):
    """Haurwitz clear-sky global horizontal irradiance [W/m2]."""
    c = np.clip(cosz, 1e-3, None)
    return np.where(cosz > 0.01, 1098.0 * c * np.exp(-0.057 / c), 0.0)


def ew_mean(x, valid, tau_blocks):
    """Causal exponentially weighted mean of x over the valid samples (gaps excluded, not zero)."""
    a = np.exp(-1.0 / tau_blocks)
    xv = np.where(valid, x, 0.0)
    num = lfilter([1 - a], [1, -a], xv)
    den = lfilter([1 - a], [1, -a], valid.astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0.05, num / den, np.nan)


# ---------------------------------------------------------------------------------------------
# Measured global radiation (SwissMetNet OGD 10-min gre000z0)
# ---------------------------------------------------------------------------------------------
def read_smn_gre(abbr, parts):
    t_all, g_all = [], []
    for p in parts:
        f = SMN_DIR / f"ogd-smn_{abbr}_t_{p}.csv"
        if not f.exists():
            continue
        with open(f, encoding="latin-1") as fh:
            head = fh.readline().strip().split(";")
            ig, it = head.index("gre000z0"), head.index("reference_timestamp")
            for line in fh:
                s = line.rstrip("\n").split(";")
                ts = s[it]
                if not (ts[6:10] in ("2025", "2026")):
                    continue
                if s[ig] == "":
                    continue
                from datetime import datetime
                d = datetime.strptime(ts, "%d.%m.%Y %H:%M")
                t_all.append((d - datetime(1970, 1, 1)).total_seconds() / 86400.0)
                g_all.append(float(s[ig]))
    if not t_all:
        return None
    t, g = np.array(t_all), np.array(g_all)
    t, idx = np.unique(t, return_index=True)
    return t, g[idx]


def smn_prev24(tsmn, gsmn, t0):
    """Mean measured global radiation over the 24 h before t0 [MJ/m2/d]; needs >= 90 % coverage."""
    out = np.full(len(t0), np.nan)
    for k, x in enumerate(t0):
        i0, i1 = np.searchsorted(tsmn, x - 1.0), np.searchsorted(tsmn, x)
        if i1 - i0 >= 0.9 * 144:
            out[k] = np.mean(gsmn[i0:i1]) * 86400 / 1e6
    return out


# ---------------------------------------------------------------------------------------------
# Stage 1: per-night drivers
# ---------------------------------------------------------------------------------------------
def unit_drivers(args):
    meta, rows = args
    s = meta["stream"]
    wmo, ident = s.rsplit("_", 1)
    lat, lon = float(meta["lat"]), float(meta["lon"])
    cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
    if cons is None or cache is None:
        return None
    m = (cons["module"] == meta["module"]) & (cons["c"] > 0) & np.isfinite(cons["c"])
    cons = {k: v[m] for k, v in cons.items()}
    cons = {k: v[~flag_outliers(cons["c"])] for k, v in cons.items()}
    tm = 0.5 * (cons["t0"] + cons["t1"])
    order = np.argsort(tm)
    tm, t0a, t1a = tm[order], cons["t0"][order], cons["t1"][order]
    j = np.searchsorted(tm, rows["t"] - 1e-7)
    j = np.clip(j, 0, len(tm) - 1)
    okm = np.abs(tm[j] - rows["t"]) < 1e-5
    t0 = np.where(okm, t0a[j], np.nan)
    t1 = np.where(okm, t1a[j], np.nan)
    n = len(rows["t"])
    D = {k: np.full(n, np.nan) for k in DRIVERS}

    # --- geometry -------------------------------------------------------------------------------
    tl = t0 + lon / 360.0                                                     # local solar time
    D["daylen_h"] = day_length(np.floor(tl - 0.5) + 0.5 - lon / 360.0, lat)  # preceding local noon
    step = 2.0 / 1440.0
    back = np.arange(-720, 1) * step            # 24 h before t0
    fwd = np.arange(0, 721) * step              # 24 h after t1
    for k in range(n):
        if not okm[k]:
            continue
        cz_b, _, _ = solar(t0[k] + back, lat, lon)
        up = np.where(cz_b > SUNSET_COSZ)[0]
        sunset = t0[k] + back[up[-1]] if up.size else np.nan
        cz_f, _, _ = solar(t1[k] + fwd, lat, lon)
        upf = np.where(cz_f > SUNSET_COSZ)[0]
        sunrise = t1[k] + fwd[upf[0]] if upf.size else np.nan
        D["sunset_to_t0_h"][k] = (t0[k] - sunset) * 24
        D["sunset_to_tmid_h"][k] = (rows["t"][k] - sunset) * 24
        D["t1_to_sunrise_h"][k] = (sunrise - t1[k]) * 24
        D["night_len_h"][k] = (sunrise - sunset) * 24
        tw = np.linspace(t0[k], t1[k], 25)
        czw, _, _ = solar(tw, lat, lon)
        sza = np.rad2deg(np.arccos(np.clip(czw, -1, 1)))
        D["sza_win_mean"][k] = sza.mean()
        D["sza_win_min"][k] = sza.min()

    # --- regular 10-min grid of the cache -------------------------------------------------------
    day0, day1 = int(cache["day"].min()), int(cache["day"].max())
    nd = day1 - day0 + 1
    names = list(cache["hk_vars"])
    grid = {v: np.full((nd, BLOCKS_PER_DAY), np.nan) for v in ("Tint", "Text", "Tlom", "Tdet", "cf")}
    di = cache["day"].astype(int) - day0
    for v, hv in (("Tint", "temp_int"), ("Text", "temp_ext"), ("Tlom", "temperature_optical_module"),
                  ("Tdet", "temperature_detector")):
        grid[v][di] = cache["hk"][:, :, names.index(hv)]
    na = cache["n_all"].astype(float)
    nc = cache["n_clear"].astype(float)
    grid["cf"][di] = np.where(na > 0, nc / np.maximum(na, 1), np.nan)
    grid = {k: v.ravel() for k, v in grid.items()}
    tb = day0 + (np.arange(nd * BLOCKS_PER_DAY) + 0.5) / BLOCKS_PER_DAY
    cz, _, rf = solar(tb, lat, lon)
    sun = cz > 0.0
    night = cz < np.cos(np.deg2rad(100.0))
    gcs = ghi_clear(cz)
    toa = 1361.0 * rf * np.clip(cz, 0, None)
    # Kasten-Czeplak cloud factor with the cloud fraction 1 - clear fraction (ceilometer's own sky view)
    kfac = 1.0 - 0.75 * np.clip(1.0 - grid["cf"], 0, 1) ** 3.4
    lom_sp = np.nanmedian(grid["Tlom"][night])
    det_sp = np.nanmedian(grid["Tdet"][night])
    lom_x = grid["Tlom"] - lom_sp
    lom_out = lom_x >= LOM_EXC
    det_out = np.abs(grid["Tdet"] - det_sp) > DET_EXC
    v_lom = np.isfinite(lom_x)
    v_det = np.isfinite(grid["Tdet"])
    v_T = np.isfinite(grid["Tint"])
    # proxy irradiance per block: clear-sky x cloud factor; gaps keep the clear-sky value scaled later
    v_cf = np.isfinite(kfac) & sun
    ghi_p = np.where(v_cf, gcs * np.where(np.isfinite(kfac), kfac, 0.0), 0.0)

    ew = {}
    for tau in TAUS_LOM:
        ew[f"lom_ew{tau:g}d"] = ew_mean(lom_out.astype(float), v_lom, tau * BLOCKS_PER_DAY) * 24.0
    for tau in TAUS_T:
        ew[f"Tint_ew{tau:g}d"] = ew_mean(np.nan_to_num(grid["Tint"]), v_T, tau * BLOCKS_PER_DAY)
    # EW of the proxy irradiance: weight gaps out through the clear-sky-weighted cloud factor
    for tau in TAUS_GHI:
        kk = ew_mean(np.where(np.isfinite(kfac), kfac, 0.0) * gcs, np.isfinite(kfac) | ~sun, tau * BLOCKS_PER_DAY)
        ew[f"ghi_proxy_ew{tau:g}d"] = kk

    for k in range(n):
        if not okm[k]:
            continue
        i1 = int(np.floor((t0[k] - day0) * BLOCKS_PER_DAY))      # first block at/after t0 (exclusive end)
        i0 = i1 - BLOCKS_PER_DAY
        if i0 < 0 or i1 > len(tb):
            continue
        sl = slice(i0, i1)
        s_sun = sun[sl]
        D["toa_prev"][k] = toa[sl].sum() * 600 / 1e6
        D["csky_prev"][k] = gcs[sl].sum() * 600 / 1e6
        w = gcs[sl] * np.isfinite(kfac[sl])
        if w.sum() > 0.7 * gcs[sl].sum() and gcs[sl].sum() > 0:
            kf = np.nansum(np.where(np.isfinite(kfac[sl]), kfac[sl], 0) * gcs[sl]) / w.sum()
            D["ghi_proxy_prev"][k] = D["csky_prev"][k] * kf
            D["cf_day"][k] = np.nansum(np.where(np.isfinite(grid["cf"][sl]), grid["cf"][sl], 0) * gcs[sl]) / w.sum()
        ti = grid["Tint"][sl][s_sun]
        te = grid["Text"][sl][s_sun]
        if np.isfinite(ti).sum() >= 0.5 * max(s_sun.sum(), 1) and s_sun.sum() > 0:
            D["Tint_daymax"][k] = np.nanmax(ti)
            D["dT_day_night"][k] = np.nanmax(ti) - rows["T_int"][k]
            dd = ti - te
            if np.isfinite(dd).any():
                D["Tint_Text_daymax"][k] = np.nanmax(dd)
        if v_lom[sl].sum() >= 0.7 * BLOCKS_PER_DAY:
            D["lom_hours_prev"][k] = lom_out[sl][v_lom[sl]].sum() / 6.0 * BLOCKS_PER_DAY / v_lom[sl].sum()
            D["lom_maxexc_prev"][k] = np.nanmax(lom_x[sl])
        if v_det[sl].sum() >= 0.7 * BLOCKS_PER_DAY:
            D["det_hours_prev"][k] = det_out[sl][v_det[sl]].sum() / 6.0 * BLOCKS_PER_DAY / v_det[sl].sum()
        for name, arr in ew.items():
            D[name][k] = arr[i1 - 1]
    D["Tint_Text_night"] = rows["T_int"] - rows["T_ext"]
    if s in SMN:
        abbr, parts = SMN[s]
        r = read_smn_gre(abbr, parts)
        if r is not None:
            D["smn_gre_prev"] = smn_prev24(r[0], r[1], np.where(okm, t0, np.nan))
    D["t0"], D["t1"] = t0, t1
    info = dict(lom_sp=float(lom_sp), det_sp=float(det_sp),
                lom_out_frac_day=float(np.nanmean(lom_out[sun & v_lom])) if (sun & v_lom).any() else np.nan)
    return s, D, info


def build_nights(units, rebuild=False):
    f = Q16 / "q16b_nights.npz"
    if f.exists() and not rebuild:
        z = np.load(f, allow_pickle=True)
        return z["drivers"].item(), z["info"].item()
    with Pool(N_WORKERS) as pool:
        res = [r for r in pool.map(unit_drivers, units, chunksize=2) if r is not None]
    drivers = {s: D for s, D, _ in res}
    info = {s: i for s, _, i in res}
    np.savez_compressed(f, drivers=np.array(drivers, dtype=object), info=np.array(info, dtype=object))
    by = {m["stream"]: (m, r) for m, r in units}
    with open(Q16 / "q16b_nights.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        cols = ["t", "lnC", "T_int", "T_ext"]
        dcols = ["t0", "t1"] + list(DRIVERS)
        wr.writerow(["stream", "module"] + cols + dcols)
        for s, D in drivers.items():
            m, r = by[s]
            for i in range(len(r["t"])):
                wr.writerow([s, m["module"]] + [f"{r[c][i]:.6g}" for c in cols] + [f"{D[c][i]:.5g}" for c in dcols])
    print(f"drivers: {len(drivers)} units -> {f}", flush=True)
    return drivers, info


# ---------------------------------------------------------------------------------------------
# kappa per unit (same sources as q15f)
# ---------------------------------------------------------------------------------------------
def load_kappas():
    k_in = {r["stream"]: float(r["b_T_ctrl"]) for r in csv.DictReader(open(OUT / "q2a_intranight.csv", encoding="utf-8"))}
    k_m2 = {}
    for r in csv.DictReader(open(OUT / "q1_slopes.csv", encoding="utf-8")):
        if r["arm"] == "A" and r["tvar"] == "temp_int" and r["b2"] not in ("", "nan"):
            k_m2[(r["stream"], r["era"])] = float(r["b2"])
    noon = {r["stream"]: float(r["noon_dT"]) for r in csv.DictReader(open(OUT / "q3_daytime_excess.csv", encoding="utf-8"))
            if r["noon_dT"] not in ("", "nan")}
    return k_in, k_m2, noon


def residuals(r, kappa):
    T = r["T_int"]
    return r["lnC"] - kappa / 100.0 * (T - np.nanmedian(T))


# ---------------------------------------------------------------------------------------------
# Tests 1-3 per unit
# ---------------------------------------------------------------------------------------------
def analyse_unit(meta, r, D, kin, km2, noon_dT):
    s = meta["stream"]
    t = r["t"]
    T = r["T_int"]
    out_units = dict(stream=s, site=meta["site"], lat=meta["lat"], lon=meta["lon"], module=meta["module"],
                     n=int(np.isfinite(r["lnC"]).sum()), kappa_in=kin, kappa_m2=km2)
    pC, _, _ = harmonic(t, r["lnC"])
    pT, _, _ = harmonic(t, T)
    out_units.update(amp_C=100 * abs(pC), peak_C=peak_doy(pC), amp_T=abs(pT), peak_T=peak_doy(pT))
    res = {}
    for lab, k in (("in", kin), ("m2", km2)):
        if np.isfinite(k):
            res[lab] = residuals(r, k)
            p, _, _ = harmonic(t, res[lab])
            out_units[f"amp_r{lab}"] = 100 * abs(p)
            out_units[f"peak_r{lab}"] = peak_doy(p)
            out_units[f"p_r{lab}"] = p
    la = anomaly(t, r["lnC"])
    ta = anomaly(t, T)
    drows = []
    for name in DRIVERS:
        x = D[name]
        if np.isfinite(x).sum() < 60 or np.nanstd(x) < 1e-6:
            continue
        pD, _, _ = harmonic(t, x)
        row = dict(stream=s, driver=name, n=int(np.isfinite(x).sum()), lat=meta["lat"],
                   amp_D=abs(pD), peak_D=peak_doy(pD), mean_D=float(np.nanmean(x)))
        # synoptic: ln C anomaly on T anomaly + driver anomaly
        xa = anomaly(t, x) if name not in SYNOPTIC_NONE else np.full(len(t), np.nan)
        ok = np.isfinite(la) & np.isfinite(ta) & np.isfinite(xa)
        if ok.sum() >= 50 and np.std(xa[ok]) > 1e-6:
            sd = np.std(xa[ok])
            X = np.column_stack([np.ones(ok.sum()), ta[ok], xa[ok] / sd])
            c, se, _ = ols(X, la[ok])
            # partial correlation of ln C and the driver, T_int controlled
            e1 = ols(np.column_stack([np.ones(ok.sum()), ta[ok]]), la[ok])[2]
            e2 = ols(np.column_stack([np.ones(ok.sum()), ta[ok]]), xa[ok])[2]
            row.update(g_syn_sd=100 * c[2], se_syn_sd=100 * se[2], g_syn=100 * c[2] / sd, kappa_syn=100 * c[1],
                       pcorr=float(np.corrcoef(e1, e2)[0, 1]), sd_anom=sd, n_syn=int(ok.sum()))
            # transfer: synoptic gain x driver annual harmonic, projected on the residual's phasor
            for lab in res:
                pr = out_units[f"p_r{lab}"]
                pred = c[2] / sd * pD
                row[f"transfer_{lab}"] = 100 * np.real(pred * np.conj(pr)) / abs(pr)
        for lab, y in res.items():
            pr = out_units[f"p_r{lab}"]
            dphi = wrap_doy(peak_doy(pD) - peak_doy(pr))
            row[f"dphase_{lab}"] = dphi
            row[f"cos_{lab}"] = float(np.cos(W * dphi))
            # joint model, kappa fixed: residual = a + gamma D
            ok = np.isfinite(y) & np.isfinite(x)
            if ok.sum() < 60:
                continue
            c, se, e = ols(np.column_stack([np.ones(ok.sum()), x[ok]]), y[ok])
            p_after, _, _ = harmonic(t[ok], e)
            p_before, _, _ = harmonic(t[ok], y[ok])
            row[f"gamma_{lab}"] = 100 * c[1]
            row[f"gamma_se_{lab}"] = 100 * se[1]
            row[f"amp_before_{lab}"] = 100 * abs(p_before)
            row[f"amp_after_{lab}"] = 100 * abs(p_after)
            row[f"scat_before_{lab}"] = 100 * np.std(y[ok])
            row[f"scat_after_{lab}"] = 100 * np.std(e)
            # out-of-sample: fit 2025, test 2026
            row.update({f"{k}_{lab}": v for k, v in oos(t[ok], y[ok], x[ok]).items()})
        # joint model, kappa free (raw): does the driver bring the raw kappa back to the anomaly kappa?
        ok = np.isfinite(r["lnC"]) & np.isfinite(T) & np.isfinite(x)
        if ok.sum() >= 60:
            c0, _, _ = ols(np.column_stack([np.ones(ok.sum()), T[ok]]), r["lnC"][ok])
            c1, se1, e1 = ols(np.column_stack([np.ones(ok.sum()), T[ok], x[ok]]), r["lnC"][ok])
            p1, _, _ = harmonic(t[ok], e1)
            row.update(kappa_raw=100 * c0[1], kappa_joint=100 * c1[1], kappa_joint_se=100 * se1[1],
                       gamma_free=100 * c1[2], amp_after_free=100 * abs(p1))
            if np.isfinite(noon_dT):
                row["day_excess_kappa_m2"] = km2 * noon_dT if np.isfinite(km2) else np.nan
                row["day_excess_joint"] = 100 * c1[1] * noon_dT
        drows.append(row)
    return out_units, drows


def monthly_rms(t, e):
    mon = np.floor(t / 30.4375).astype(int)
    means = [np.mean(e[mon == m]) for m in np.unique(mon) if (mon == m).sum() >= 5]
    return float(np.sqrt(np.mean(np.square(np.array(means) - np.mean(means))))) if len(means) >= 4 else np.nan


def oos(t, y, x):
    """Fit on 2025, test on 2026. Metrics on the test set after removing its mean (offset/ageing)."""
    fit, test = t < SPLIT, t >= SPLIT
    if fit.sum() < 60 or test.sum() < 40:
        return {}
    d = np.mod(t, 365.25)
    H = np.column_stack([np.cos(W * d), np.sin(W * d)])
    out = {}
    models = {"k": np.ones((len(t), 1)),
              "kD": np.column_stack([np.ones(len(t)), x]),
              "kH": np.column_stack([np.ones(len(t)), H])}
    for name, X in models.items():
        c, *_ = np.linalg.lstsq(X[fit], y[fit], rcond=None)
        e = y[test] - X[test] @ c
        e = e - e.mean()
        out[f"oos_mon_{name}"] = 100 * monthly_rms(t[test], e)
        out[f"oos_sd_{name}"] = 100 * np.std(e)
    return out


# ---------------------------------------------------------------------------------------------
# Network-level summaries
# ---------------------------------------------------------------------------------------------
def pooled_gamma(units, drivers, k_m2, name, lab="m2", k_in=None):
    """One network gamma: residual (kappa fixed) on the driver with unit fixed effects, fit on 2025,
    applied to 2026 of every unit. Also the bootstrap SE over units of the synoptic pooled gain."""
    fit_y, fit_x, test = [], [], []
    for m, r in units:
        s = m["stream"]
        if s == PAYERNE or s not in drivers:
            continue
        k = k_m2.get((s, m["module"]), np.nan) if lab == "m2" else k_in.get(s, np.nan)
        if not np.isfinite(k):
            continue
        y = residuals(r, k)
        x = drivers[s][name]
        t = r["t"]
        ok = np.isfinite(y) & np.isfinite(x)
        f, te = ok & (t < SPLIT), ok & (t >= SPLIT)
        if f.sum() < 60 or te.sum() < 40:
            continue
        fit_y.append(y[f] - y[f].mean())
        fit_x.append(x[f] - x[f].mean())
        test.append((t[te], y[te], x[te], y[f].mean(), x[f].mean()))
    if not fit_y:
        return {}
    Y, X = np.concatenate(fit_y), np.concatenate(fit_x)
    g = float(X @ Y / (X @ X))
    red, red_h = [], []
    for tt, yy, xx, my, mx in test:
        e0 = yy - yy.mean()
        e1 = yy - g * xx
        e1 = e1 - e1.mean()
        a, b = monthly_rms(tt, e0), monthly_rms(tt, e1)
        if np.isfinite(a) and np.isfinite(b):
            red.append(1 - b / a)
    return dict(gamma_pooled=100 * g, oos_pooled_monthly_reduction=float(np.median(red)) if red else np.nan,
                oos_pooled_n=len(red))


def pooled_synoptic(units, drivers, name, n_boot=300, seed=1):
    """Pooled synoptic gain (ln C anomaly on T and driver anomalies) with a bootstrap over units."""
    blocks = []
    for m, r in units:
        s = m["stream"]
        if s == PAYERNE or s not in drivers or name in SYNOPTIC_NONE:
            continue
        t = r["t"]
        la, ta, xa = anomaly(t, r["lnC"]), anomaly(t, r["T_int"]), anomaly(t, drivers[s][name])
        ok = np.isfinite(la) & np.isfinite(ta) & np.isfinite(xa)
        if ok.sum() >= 50:
            blocks.append((la[ok], ta[ok], xa[ok]))
    if len(blocks) < 10:
        return {}

    def fit(bl):
        y = np.concatenate([b[0] for b in bl])
        X = np.column_stack([np.ones(len(y)), np.concatenate([b[1] for b in bl]), np.concatenate([b[2] for b in bl])])
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        return c
    c = fit(blocks)
    rng = np.random.default_rng(seed)
    boots = np.array([fit([blocks[i] for i in rng.integers(0, len(blocks), len(blocks))]) for _ in range(n_boot)])
    sd_x = np.std(np.concatenate([b[2] for b in blocks]))
    return dict(g_pool=100 * c[2], g_pool_se=100 * boots[:, 2].std(), g_pool_sd=100 * c[2] * sd_x,
                kappa_pool=100 * c[1], n_units_pool=len(blocks))


def pairs_test(units, drivers, k_m2):
    """Co-located twins: A/B ratio (kappa-normalised) synoptic anomalies on driver differences."""
    by = {m["stream"]: (m, r) for m, r in units}
    rows = []
    for wmo, idents in PAIR_SITES.items():
        members = [f"{wmo}_{i}" for i in idents if f"{wmo}_{i}" in by and f"{wmo}_{i}" in drivers]
        for a, b in itertools.combinations(members, 2):
            (ma, ra), (mb, rb) = by[a], by[b]
            ka, kb = k_m2.get((a, ma["module"]), np.nan), k_m2.get((b, mb["module"]), np.nan)
            if not (np.isfinite(ka) and np.isfinite(kb)):
                continue
            common, ia, ib = np.intersect1d(ra["night"], rb["night"], return_indices=True)
            if len(common) < 60:
                continue
            t = 0.5 * (ra["t"][ia] + rb["t"][ib])
            o = np.argsort(t)
            t = t[o]
            y = (residuals(ra, ka)[ia] - residuals(rb, kb)[ib])[o]
            ya = anomaly(t, y)
            pr, _, _ = harmonic(t, y)
            for name in ("Tint_daymax", "dT_day_night", "Tint_Text_night", "lom_hours_prev", "lom_ew7d",
                         "det_hours_prev", "Tint_ew7d", "sunset_to_t0_h", "sza_win_mean"):
                dx = (drivers[a][name][ia] - drivers[b][name][ib])[o]
                row = dict(pair=f"{a}/{b.rsplit('_', 1)[1]}", driver=name, n=int(len(t)),
                           amp_ratio=100 * abs(pr), peak_ratio=peak_doy(pr), sd_diff=float(np.nanstd(dx)))
                xa = anomaly(t, dx)
                ok = np.isfinite(ya) & np.isfinite(xa)
                if ok.sum() >= 40 and np.std(xa[ok]) > 1e-6:
                    c, se, _ = ols(np.column_stack([np.ones(ok.sum()), xa[ok]]), ya[ok])
                    row.update(g_syn=100 * c[1], se_syn=100 * se[1], r=float(np.corrcoef(xa[ok], ya[ok])[0, 1]))
                ok = np.isfinite(y) & np.isfinite(dx)
                if ok.sum() >= 40 and np.std(dx[ok]) > 1e-6:
                    pd_, _, _ = harmonic(t[ok], dx[ok])
                    c, _, e = ols(np.column_stack([np.ones(ok.sum()), dx[ok]]), y[ok])
                    pe, _, _ = harmonic(t[ok], e)
                    row.update(amp_diff_driver=abs(pd_), g_full=100 * c[1], amp_ratio_after=100 * abs(pe))
                rows.append(row)
    return rows


# ---------------------------------------------------------------------------------------------
# Test 4: the q15e open point (synoptic kappa at 7-8 km vs in the window)
# ---------------------------------------------------------------------------------------------
BANDS = [(0.2, 0.5), (0.5, 0.8), (0.8, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.0),
         (3.0, 4.0), (4.0, 5.0), (5.0, 6.0), (6.0, 7.0), (7.0, 8.0)]


T4_DRIVERS = ("ghi_proxy_prev", "cf_day", "sunset_to_t0_h", "sza_win_mean", "dT_day_night", "Tint_daymax",
              "Tint_Text_night", "lom_ew7d", "Tint_ew7d", "Tint_ew21d")


def test4(units, drivers):
    ser = np.load(Q15 / "q15e_series.npz")
    keys = set(k.split("__")[0] for k in ser.files)
    specs = {
        "T_int": [],
        "T_int + T_ext": ["T_ext"],
        "T_int + insolation": ["ghi_proxy_prev", "cf_day"],
        "T_int + day heating": ["dT_day_night", "Tint_Text_night"],
        "T_int + aerosol": ["aer"],
        "T_int + all": ["T_ext", "ghi_proxy_prev", "cf_day", "dT_day_night", "aer"],
    }
    pooled = {sp: {b: [] for b in ["win", "lnc"] + list(range(len(BANDS)))} for sp in specs}
    split_rows = {b: [] for b in ["win", "lnc"] + list(range(len(BANDS)))}
    unit_rows = []
    for m, r in units:
        s = m["stream"]
        if s == PAYERNE or s not in keys or s not in drivers:
            continue
        lnc_e = ser[f"{s}__lnc"]
        idx = {round(float(v), 9): i for i, v in enumerate(lnc_e)}
        j = np.array([idx.get(round(float(v), 9), -1) for v in r["lnC"]])
        if (j >= 0).sum() < 60:
            continue
        t = r["t"]
        lnx = np.where(j[:, None] >= 0, ser[f"{s}__lnx"][np.clip(j, 0, None)], np.nan)
        lnx_win = np.where(j >= 0, ser[f"{s}__lnx_win"][np.clip(j, 0, None)], np.nan)
        cov = dict(T_int=r["T_int"], T_ext=r["T_ext"], aer=r["aer"], **{k: drivers[s][k] for k in
                                                                        ("ghi_proxy_prev", "cf_day", "dT_day_night",
                                                                         "Tint_Text_night")})
        A = {k: anomaly(t, v) for k, v in cov.items()}
        Y = {"win": anomaly(t, lnx_win), "lnc": anomaly(t, r["lnC"])}
        for jb in range(len(BANDS)):
            Y[jb] = anomaly(t, lnx[:, jb])
        for sp, extra in specs.items():
            cols = ["T_int"] + extra
            for b, y in Y.items():
                ok = np.isfinite(y) & np.all([np.isfinite(A[c]) for c in cols], axis=0)
                if ok.sum() >= 50:
                    pooled[sp][b].append((y[ok], np.column_stack([A[c][ok] for c in cols])))
        # T_int anomaly split into its ambient part (explained by T_ext) and the internal remainder
        ok0 = np.isfinite(A["T_int"]) & np.isfinite(A["T_ext"])
        if ok0.sum() >= 50:
            c, *_ = np.linalg.lstsq(np.column_stack([np.ones(ok0.sum()), A["T_ext"][ok0]]), A["T_int"][ok0], rcond=None)
            amb = np.full(len(t), np.nan)
            amb[ok0] = c[0] + c[1] * A["T_ext"][ok0]
            intern = A["T_int"] - amb
            unit_rows.append(dict(stream=s, tint_on_text=c[1], r2_tint_text=float(np.corrcoef(A["T_int"][ok0], A["T_ext"][ok0])[0, 1] ** 2),
                                  sd_tint_anom=float(np.nanstd(A["T_int"][ok0])), sd_internal=float(np.nanstd(intern[ok0]))))
            for b, y in Y.items():
                ok = np.isfinite(y) & np.isfinite(amb) & np.isfinite(intern)
                if ok.sum() >= 50:
                    split_rows[b].append((y[ok], np.column_stack([amb[ok], intern[ok]])))
    rng = np.random.default_rng(3)
    # sunlight / heating drivers on the clear-sky molecular-normalised signal (aerosol-free 7-8 km vs window)
    drv_blocks = {}
    for m, r in units:
        s = m["stream"]
        if s == PAYERNE or s not in keys or s not in drivers:
            continue
        lnc_e = ser[f"{s}__lnc"]
        idx = {round(float(v), 9): i for i, v in enumerate(lnc_e)}
        j = np.array([idx.get(round(float(v), 9), -1) for v in r["lnC"]])
        if (j >= 0).sum() < 60:
            continue
        t = r["t"]
        x78 = np.where(j >= 0, ser[f"{s}__lnx"][np.clip(j, 0, None), len(BANDS) - 1], np.nan)
        xw = np.where(j >= 0, ser[f"{s}__lnx_win"][np.clip(j, 0, None)], np.nan)
        ta = anomaly(t, r["T_int"])
        for b, yv in (("7.0-8.0", x78), ("win", xw), ("lnc", r["lnC"])):
            ya = anomaly(t, yv)
            for dn in T4_DRIVERS:
                xa = anomaly(t, drivers[s][dn])
                ok = np.isfinite(ya) & np.isfinite(ta) & np.isfinite(xa)
                if ok.sum() >= 50 and np.std(xa[ok]) > 1e-6:
                    drv_blocks.setdefault((b, dn), []).append((ya[ok], np.column_stack([ta[ok], xa[ok]])))

    def pfit(blocks, n_boot=200):
        if len(blocks) < 10:
            return np.nan, np.nan, np.nan, np.nan, 0
        def f(bl):
            y = np.concatenate([b[0] for b in bl])
            X = np.vstack([b[1] for b in bl])
            X = np.column_stack([np.ones(len(y)), X])
            c, *_ = np.linalg.lstsq(X, y, rcond=None)
            return c[1:3] if X.shape[1] > 2 else np.r_[c[1], np.nan]
        c = f(blocks)
        bs = np.array([f([blocks[i] for i in rng.integers(0, len(blocks), len(blocks))]) for _ in range(n_boot)])
        return 100 * c[0], 100 * bs[:, 0].std(), 100 * c[1], 100 * bs[:, 1].std(), len(blocks)
    rows = []
    for sp in specs:
        for b in ["win", "lnc"] + list(range(len(BANDS))):
            k, se, _, _, nu = pfit(pooled[sp][b])
            band = b if isinstance(b, str) else f"{BANDS[b][0]:.1f}-{BANDS[b][1]:.1f}"
            rows.append(dict(spec=sp, band=band, kappa=k, se=se, n_units=nu))
    for (b, dn), bl in drv_blocks.items():
        kt, st, kd, sd_, nu = pfit(bl)
        sdx = float(np.std(np.concatenate([x[1][:, 1] for x in bl])))
        rows.append(dict(spec=f"driver: {dn}", band=b, kappa=kt, se=st, gamma=kd, gamma_se=sd_,
                         gamma_per_sd=kd * sdx, gamma_per_sd_se=sd_ * sdx, n_units=nu))
    for b in ["win", "lnc"] + list(range(len(BANDS))):
        ka, sa, ki, si, nu = pfit(split_rows[b])
        band = b if isinstance(b, str) else f"{BANDS[b][0]:.1f}-{BANDS[b][1]:.1f}"
        rows.append(dict(spec="split: ambient part", band=band, kappa=ka, se=sa, n_units=nu))
        rows.append(dict(spec="split: internal part", band=band, kappa=ki, se=si, n_units=nu))
    return rows, unit_rows


# ---------------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------------
def composites(units, drivers, k_m2, names, lab="m2"):
    """Network composite by month of year: unit-demeaned residual and drivers (driver scaled by the pooled
    full-sample gamma so it is in % of C_L)."""
    mons = np.arange(12)
    acc_r = [[] for _ in mons]
    acc = {n: [[] for _ in mons] for n in names}
    ys, xs = {n: [] for n in names}, {n: [] for n in names}
    for m, r in units:
        s = m["stream"]
        if s == PAYERNE or s not in drivers:
            continue
        k = k_m2.get((s, m["module"]), np.nan)
        if not np.isfinite(k):
            continue
        y = residuals(r, k)
        mo = ((np.mod(r["t"], 365.25) / 30.4375).astype(int)) % 12
        yc = y - np.nanmean(y)
        for i in mons:
            sel = (mo == i) & np.isfinite(yc)
            if sel.sum() >= 3:
                acc_r[i].append(np.mean(yc[sel]))
        for n in names:
            x = drivers[s][n]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 60:
                continue
            xc = x - np.nanmean(x[ok])
            ys[n].append(y[ok] - y[ok].mean())
            xs[n].append(xc[ok])
            for i in mons:
                sel = (mo == i) & ok
                if sel.sum() >= 3:
                    acc[n][i].append(np.mean(xc[sel]))
    comp_r = np.array([100 * np.median(a) if a else np.nan for a in acc_r])
    out = {}
    for n in names:
        if not xs[n]:
            continue
        X, Y = np.concatenate(xs[n]), np.concatenate(ys[n])
        g = X @ Y / (X @ X)
        out[n] = (100 * g * np.array([np.median(a) if a else np.nan for a in acc[n]]), 100 * g)
    return comp_r, out


def make_figures(units, drivers, k_m2, unit_rows, drow, t4rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mon_mid = (np.arange(12) + 0.5) * 30.4375
    # Fig 1: composites
    geo = ["daylen_h", "csky_prev", "ghi_proxy_ew21d", "sza_win_mean", "sunset_to_t0_h"]
    heat = ["Tint_daymax", "dT_day_night", "lom_ew21d", "Tint_ew21d", "Tint_Text_night"]
    comp_r, out = composites(units, drivers, k_m2, geo + heat)
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.2), sharey=True)
    for ax, names, title in ((axs[0], geo, "sunlight / geometry drivers"), (axs[1], heat, "instrument heating drivers")):
        ax.plot(mon_mid, comp_r, "k-o", lw=2.5, label="residual ln C - kappa_m2 T_int")
        for n in names:
            if n in out:
                ax.plot(mon_mid, out[n][0], "-s", ms=4, label=f"{n} (gamma fit, {out[n][1]:+.3g} %/unit)")
        ax.axvline(172, color="0.6", ls=":")
        ax.axhline(0, color="0.8")
        ax.set_xlabel("day of year")
        ax.set_title(title)
        ax.legend(fontsize=7.5, loc="lower center")
        ax.set_xlim(0, 365)
    axs[0].set_ylabel("network composite [% of C_L] (unit means removed)")
    fig.suptitle("q16b: residual annual cycle of the CHM15k constant vs candidate sunlight / heating drivers")
    fig.tight_layout()
    fig.savefig(Q16 / "q16b_composites.png", dpi=130)
    plt.close(fig)

    # Fig 2: phases (driver peak DOY median + IQR) vs residual
    fig, ax = plt.subplots(figsize=(14, 6.5))
    ylab = []
    yy = 0
    for lab, key in (("residual (kappa_in)", "peak_rin"), ("residual (kappa_m2)", "peak_rm2"), ("ln C (raw)", "peak_C"),
                     ("night T_int", "peak_T")):
        v = np.array([u.get(key, np.nan) for u in unit_rows if u["stream"] != PAYERNE], float)
        med = circ_median_doy(v)
        q = np.nanpercentile(wrap_doy(v - med), [25, 75]) + med
        ax.errorbar(med, yy, xerr=[[med - q[0]], [q[1] - med]], fmt="ks", ms=8)
        ylab.append(lab)
        yy += 1
    for name in DRIVERS:
        v = np.array([d["peak_D"] for d in drow if d["driver"] == name and d["stream"] != PAYERNE], float)
        if v.size < 10:
            continue
        a = np.array([d["amp_D"] for d in drow if d["driver"] == name and d["stream"] != PAYERNE], float)
        # drivers whose cycle is tiny or anti-seasonal: still plot the phase of the maximum
        med = circ_median_doy(v)
        q = np.nanpercentile(wrap_doy(v - med), [25, 75]) + med
        ax.errorbar(med, yy, xerr=[[med - q[0]], [q[1] - med]], fmt="o", color="C0" if "geometry" in DRIVERS[name][0]
                    or "insolation" in DRIVERS[name][0] else "C3", ms=6)
        ylab.append(f"{name} (amp {np.nanmedian(a):.3g} {DRIVERS[name][1]})")
        yy += 1
    ax.set_yticks(range(len(ylab)))
    ax.set_yticklabels(ylab, fontsize=7.5)
    ax.axvspan(148, 177, color="0.85", zorder=0)
    ax.set_xlabel("DOY of the annual maximum (network median, bar = IQR)")
    ax.set_xlim(0, 365)
    ax.invert_yaxis()
    ax.set_title("Phase test: annual maximum of each driver vs the residual (grey band = residual peak range)")
    fig.tight_layout()
    fig.savefig(Q16 / "q16b_phases.png", dpi=130)
    plt.close(fig)

    # Fig 3: latitude contrast
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    u = [x for x in unit_rows if x["stream"] != PAYERNE and np.isfinite(x.get("amp_rm2", np.nan))]
    lat = np.array([x["lat"] for x in u])
    p_r = np.array([x["p_rm2"] for x in u])
    inph = 100 * np.real(p_r * np.exp(-1j * W * 172.0))
    dl = {d["stream"]: d for d in drow if d["driver"] == "daylen_h"}
    dla = np.array([dl[x["stream"]]["amp_D"] if x["stream"] in dl else np.nan for x in u])
    axs[0].scatter(lat, inph, s=14)
    axs[0].set_xlabel("latitude [deg N]")
    axs[0].set_ylabel("residual amplitude in phase with DOY 172 [%]")
    axs[0].set_title("residual (kappa_m2) vs latitude")
    axs[1].scatter(dla, inph, s=14)
    ok = np.isfinite(dla) & np.isfinite(inph)
    g = np.sum(dla[ok] * inph[ok]) / np.sum(dla[ok] ** 2)
    xx = np.linspace(0, np.nanmax(dla), 10)
    axs[1].plot(xx, g * xx, "r-", label=f"proportional (day-length hypothesis): {g:.2f} %/h")
    axs[1].set_xlabel("day-length annual amplitude [h]")
    axs[1].set_title("latitude contrast: residual vs day-length amplitude")
    axs[1].legend(fontsize=8)
    sz = {d["stream"]: d for d in drow if d["driver"] == "sza_win_mean"}
    sza = np.array([sz[x["stream"]]["amp_D"] if x["stream"] in sz else np.nan for x in u])
    axs[2].scatter(sza, inph, s=14)
    axs[2].set_xlabel("annual amplitude of the window-mean SZA [deg] (twilight proximity)")
    axs[2].set_title("twilight-proximity contrast")
    for ax in axs:
        ax.axhline(0, color="0.8")
    fig.tight_layout()
    fig.savefig(Q16 / "q16b_latitude.png", dpi=130)
    plt.close(fig)

    # Fig 4: q15e open point - synoptic kappa per band (altitude on Y)
    fig, ax = plt.subplots(figsize=(13, 5.5))
    bands = [f"{a:.1f}-{b:.1f}" for a, b in BANDS]
    zmid = np.array([0.5 * (a + b) for a, b in BANDS])
    for sp, st in (("T_int", "k-o"), ("T_int + T_ext", "C1-s"), ("T_int + insolation", "C0-^"),
                   ("T_int + aerosol", "C2-v"), ("T_int + all", "C3-D"), ("split: ambient part", "C4--"),
                   ("split: internal part", "C5--")):
        k = np.array([next((r["kappa"] for r in t4rows if r["spec"] == sp and r["band"] == b), np.nan) for b in bands])
        se = np.array([next((r["se"] for r in t4rows if r["spec"] == sp and r["band"] == b), np.nan) for b in bands])
        ax.errorbar(k, zmid, xerr=se, fmt=st, ms=5, capsize=2, label=sp)
        kw = next((r["kappa"] for r in t4rows if r["spec"] == sp and r["band"] == "win"), np.nan)
        ax.plot([kw], [3.9], st[:2] + "*", ms=14)
    ax.axvline(0, color="0.7")
    ax.set_xlabel("synoptic kappa of the clear-sky molecular-normalised signal [%/K] (pooled, bootstrap SE)")
    ax.set_ylabel("range band centre [km] (stars at 3.9 km = calibration window)")
    ax.set_title("q15e open point: synoptic T_int sensitivity by altitude, with controls")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_xlim(-6, 12)
    fig.tight_layout()
    fig.savefig(Q16 / "q16b_q15e_kappa_z.png", dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------
def write_csv(rows, path):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys and not k.startswith("p_")]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()
                        if not k.startswith("p_")})


def med_iqr(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan, np.nan
    return tuple(np.percentile(v, [50, 25, 75]))


def solstice_curvature(units=None, drivers=None, n_boot=1000):
    """Short-baseline test independent of the harmonic: June (DOY 152-181) minus the mean of May and July
    (DOY 121-151, 182-212), per unit, network median with a bootstrap over units. A sunlight driver peaking at
    the solstice must make June stand out above its neighbours; T_int (peak ~206) does so through curvature."""
    if units is None:
        units = list(np.load(Q15 / "q15f_nightly.npz", allow_pickle=True)["units"])
    if drivers is None:
        drivers = np.load(Q16 / "q16b_nights.npz", allow_pickle=True)["drivers"].item()
    k_in, k_m2, _ = load_kappas()
    rng = np.random.default_rng(0)
    rows = []
    series = {"ln C (raw)": lambda m, r: r["lnC"],
              "residual kappa_in": lambda m, r: residuals(r, k_in.get(m["stream"], np.nan)),
              "residual kappa_m2": lambda m, r: residuals(r, k_m2.get((m["stream"], m["module"]), np.nan)),
              "T_int [K]": lambda m, r: r["T_int"] / 100.0}
    for name in ("daylen_h", "csky_prev", "ghi_proxy_prev", "sunset_to_t0_h", "sza_win_mean", "Tint_daymax",
                 "dT_day_night"):
        series[f"{name} [{DRIVERS[name][1]}]"] = (lambda nm: (lambda m, r: drivers[m["stream"]][nm] / 100.0))(name)
    for lab, fn in series.items():
        d = []
        for m, r in units:
            if m["stream"] == PAYERNE or m["stream"] not in drivers:
                continue
            y = fn(m, r)
            doy = np.mod(r["t"], 365.25)
            ok = np.isfinite(y)
            a = y[(doy >= 152) & (doy < 182) & ok]
            b = y[(((doy >= 121) & (doy < 152)) | ((doy >= 182) & (doy < 213))) & ok]
            if len(a) >= 5 and len(b) >= 10:
                d.append(100 * (a.mean() - b.mean()))
        d = np.array(d)
        if d.size < 10:
            continue
        bs = [np.median(rng.choice(d, d.size)) for _ in range(n_boot)]
        rows.append(dict(series=lab, june_minus_may_july=float(np.median(d)), se=float(np.std(bs)), n_units=int(d.size),
                         frac_negative=float(np.mean(d < 0))))
        print(f"curvature {lab:28s} {np.median(d):+.2f} +- {np.std(bs):.2f} (n={d.size})")
    write_csv(rows, Q16 / "q16b_solstice_curvature.csv")
    return rows


def mirror_test(units=None, drivers=None, n_boot=1000):
    """Solstice-mirror (hysteresis) test. Every purely astronomical driver (day length, TOA / clear-sky
    insolation, sunset time, SZA in the window, twilight proximity) is symmetric about the solstice (DOY 172):
    DOY 172-x and 172+x have the same value. Any sunlight-only explanation therefore requires the residual to be
    the same on the spring and on the mirrored autumn date. Compare spring (DOY 60-141) with its mirror
    (DOY 203-284) per unit, and plot the composite against the distance from the solstice for both branches."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if units is None:
        units = list(np.load(Q15 / "q15f_nightly.npz", allow_pickle=True)["units"])
    if drivers is None:
        drivers = np.load(Q16 / "q16b_nights.npz", allow_pickle=True)["drivers"].item()
    k_in, k_m2, _ = load_kappas()
    rng = np.random.default_rng(0)
    series = {"ln C (raw)": (lambda m, r: r["lnC"], 100.0),
              "residual kappa_in": (lambda m, r: residuals(r, k_in.get(m["stream"], np.nan)), 100.0),
              "residual kappa_m2": (lambda m, r: residuals(r, k_m2.get((m["stream"], m["module"]), np.nan)), 100.0),
              "T_int [K]": (lambda m, r: r["T_int"], 1.0)}
    for name in ("daylen_h", "csky_prev", "ghi_proxy_prev", "sunset_to_t0_h", "sza_win_mean", "Tint_daymax",
                 "dT_day_night", "lom_ew7d", "cf_day"):
        series[f"{name} [{DRIVERS[name][1]}]"] = ((lambda nm: (lambda m, r: drivers[m["stream"]][nm]))(name), 1.0)
    # clear-sky molecular-normalised signal at 7-8 km (aerosol-free) and in the window, from the q15e series
    ser = np.load(Q15 / "q15e_series.npz")
    skeys = set(k.split("__")[0] for k in ser.files)

    def x_series(col):
        def fn(m, r):
            s_ = m["stream"]
            if s_ not in skeys:
                return np.full(len(r["t"]), np.nan)
            idx = {round(float(v), 9): i for i, v in enumerate(ser[f"{s_}__lnc"])}
            j = np.array([idx.get(round(float(v), 9), -1) for v in r["lnC"]])
            src = ser[f"{s_}__lnx"][:, -1] if col == "x78" else ser[f"{s_}__lnx_win"]
            return np.where(j >= 0, src[np.clip(j, 0, None)], np.nan)
        return fn
    series["ln X 7-8 km (clear-sky)"] = (x_series("x78"), 100.0)
    series["ln X window (clear-sky)"] = (x_series("win"), 100.0)
    edges = np.arange(0, 181, 15.0)                      # distance from the solstice [d]
    rows, curves = [], {}
    for lab, (fn, scale) in series.items():
        diffs = []
        comp = {"spring": [[] for _ in edges[:-1]], "autumn": [[] for _ in edges[:-1]]}
        for m, r in units:
            if m["stream"] == PAYERNE or m["stream"] not in drivers:
                continue
            y = fn(m, r)
            ok = np.isfinite(y)
            if ok.sum() < 60:
                continue
            y = y - np.nanmean(y)
            doy = np.mod(r["t"], 365.25)
            a = y[(doy >= 60) & (doy < 141) & ok]
            b = y[(doy >= 203) & (doy < 284) & ok]
            if len(a) >= 8 and len(b) >= 8:
                diffs.append(scale * (a.mean() - b.mean()))
            dist = np.abs(wrap_doy(doy - 172.0))
            branch = np.where(wrap_doy(doy - 172.0) < 0, "spring", "autumn")
            for br in ("spring", "autumn"):
                for i in range(len(edges) - 1):
                    sel = ok & (branch == br) & (dist >= edges[i]) & (dist < edges[i + 1])
                    if sel.sum() >= 3:
                        comp[br][i].append(scale * y[sel].mean())
        d = np.array(diffs)
        if d.size < 10:
            continue
        bs = [np.median(rng.choice(d, d.size)) for _ in range(n_boot)]
        rows.append(dict(series=lab, spring_minus_autumn=float(np.median(d)), se=float(np.std(bs)),
                         n_units=int(d.size), frac_positive=float(np.mean(d > 0))))
        curves[lab] = {br: np.array([np.median(v) if len(v) >= 10 else np.nan for v in comp[br]]) for br in comp}
        print(f"mirror {lab:28s} spring-autumn {np.median(d):+.2f} +- {np.std(bs):.2f} (n={d.size})")
    write_csv(rows, Q16 / "q16b_mirror.csv")
    mid = 0.5 * (edges[:-1] + edges[1:])
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    for ax, labs, yl in ((axs[0], ["ln C (raw)", "residual kappa_m2", "ln X 7-8 km (clear-sky)"], "% (unit mean removed)"),
                         (axs[1], ["T_int [K]", "Tint_daymax [degC]", "dT_day_night [K]"], "K (unit mean removed)"),
                         (axs[2], ["daylen_h [h]", "sunset_to_t0_h [h]"], "h (unit mean removed)")):
        for k, lab in enumerate(labs):
            if lab not in curves:
                continue
            ax.plot(mid, curves[lab]["spring"], "-o", color=f"C{k}", label=f"{lab} - spring (before solstice)")
            ax.plot(mid, curves[lab]["autumn"], "--s", color=f"C{k}", label=f"{lab} - autumn (after solstice)")
        ax.set_xlabel("distance from the summer solstice [days]")
        ax.set_ylabel(yl)
        ax.axhline(0, color="0.8")
        ax.legend(fontsize=7)
    axs[0].set_title("C_L and residuals: spring vs autumn branch")
    axs[1].set_title("instrument temperatures")
    axs[2].set_title("astronomical drivers (symmetric by construction)")
    fig.suptitle("Solstice-mirror test: a sunlight-only driver needs the solid and dashed curves to coincide")
    fig.tight_layout()
    fig.savefig(Q16 / "q16b_mirror.png", dpi=130)
    plt.close(fig)
    return rows


def interannual_test(units=None, drivers=None, n_boot=500):
    """Same calendar month, different year (2026 - 2025), paired within each unit. Astronomy (day length,
    insolation at the top of the atmosphere, sunset times, twilight) is identical in both years, so any
    difference of the residual is weather/instrument driven. The 31-day anomaly tests are blind at this time
    scale; this test sees monthly-mean weather (e.g. the record-sunny spring 2025) and multi-week memories."""
    from datetime import datetime, timedelta
    if units is None:
        units = list(np.load(Q15 / "q15f_nightly.npz", allow_pickle=True)["units"])
    if drivers is None:
        drivers = np.load(Q16 / "q16b_nights.npz", allow_pickle=True)["drivers"].item()
    k_in, k_m2, _ = load_kappas()
    dnames = ("T_int", "ghi_proxy_prev", "ghi_proxy_ew7d", "ghi_proxy_ew21d", "cf_day", "Tint_daymax", "dT_day_night",
              "Tint_Text_night", "lom_ew7d", "Tint_ew21d", "aer", "zb_km", "daylen_h")
    recs = []                       # (unit idx, month, d_raw, d_rin, d_rm2, d_driver...)
    for ui, (m, r) in enumerate(units):
        s = m["stream"]
        if s == PAYERNE or s not in drivers:
            continue
        ys = {"raw": r["lnC"], "in": residuals(r, k_in.get(s, np.nan)),
              "m2": residuals(r, k_m2.get((s, m["module"]), np.nan))}
        xs = {d: (r[d] if d in r else drivers[s][d]) for d in dnames}
        dts = [datetime(1970, 1, 1) + timedelta(days=float(x)) for x in r["t"]]
        yr = np.array([d.year for d in dts])
        mo = np.array([d.month for d in dts])
        for mm in range(1, 13):
            a, b = (yr == 2025) & (mo == mm), (yr == 2026) & (mo == mm)
            if a.sum() < 4 or b.sum() < 4:
                continue
            rec = dict(unit=ui, month=mm)
            for k, v in list(ys.items()) + list(xs.items()):
                va, vb = v[a][np.isfinite(v[a])], v[b][np.isfinite(v[b])]
                rec[k] = (vb.mean() - va.mean()) if (va.size >= 3 and vb.size >= 3) else np.nan
            recs.append(rec)
    rng = np.random.default_rng(5)
    uidx = np.array([q["unit"] for q in recs])
    uu = np.unique(uidx)
    rows = []
    # per-month network medians (paired differences)
    for mm in sorted(set(q["month"] for q in recs)):
        sel = [q for q in recs if q["month"] == mm]
        if len(sel) < 20:
            continue
        row = dict(test="month median 2026-2025", month=mm, n_units=len(sel))
        for k in ("raw", "in", "m2") + dnames:
            v = np.array([q[k] for q in sel], float)
            row[k] = float(np.nanmedian(v)) * (100 if k in ("raw", "in", "m2") else 1)
        rows.append(row)
        print("interannual", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in row.items()})
    # pooled regressions on unit-month pairs, bootstrap over units; intercept absorbs ageing / offsets

    def fit(yk, xk, idx):
        y = np.array([recs[i][yk] for i in idx], float)
        x = np.array([recs[i][xk] for i in idx], float)
        mo = np.array([recs[i]["month"] for i in idx])
        ok = np.isfinite(y) & np.isfinite(x)
        if ok.sum() < 30 or np.std(x[ok]) < 1e-9:
            return np.nan, np.nan, np.nan
        c, *_ = np.linalg.lstsq(np.column_stack([np.ones(ok.sum()), x[ok]]), y[ok], rcond=None)
        r_ = np.corrcoef(x[ok], y[ok])[0, 1]
        # within-month version (month fixed effects): only between-unit weather differences remain
        X = np.column_stack([x[ok]] + [(mo[ok] == k).astype(float) for k in np.unique(mo[ok])])
        cw, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
        return c[1], r_, cw[0]
    by_unit = {u: [i for i, q in enumerate(recs) if q["unit"] == u] for u in uu}
    for yk in ("raw", "in", "m2"):
        for xk in dnames:
            if xk == "daylen_h":
                continue
            allidx = list(range(len(recs)))
            g, r_, gw = fit(yk, xk, allidx)
            bs = []
            for _ in range(n_boot):
                pick = rng.choice(uu, len(uu))
                bs.append(fit(yk, xk, [i for u in pick for i in by_unit[u]])[0])
            sdx = float(np.nanstd([q[xk] for q in recs]))
            rows.append(dict(test="pooled unit-month", y=yk, driver=xk, gamma=100 * g, gamma_se=100 * np.nanstd(bs),
                             r=r_, gamma_within_month=100 * gw, sd_dx=sdx, n=len(recs)))
            print(f"interannual pooled d{yk} on d{xk}: gamma {100 * g:+.3g} +- {100 * np.nanstd(bs):.2g} "
                  f"(r {r_:+.2f}; within-month {100 * gw:+.3g}); sd dx {sdx:.3g}")
    write_csv(rows, Q16 / "q16b_interannual.csv")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true")
    a = ap.parse_args()
    units = list(np.load(Q15 / "q15f_nightly.npz", allow_pickle=True)["units"])
    drivers, info = build_nights(units, a.rebuild)
    k_in, k_m2, noon = load_kappas()
    unit_rows, drow = [], []
    for m, r in units:
        s = m["stream"]
        if s not in drivers:
            continue
        u, d = analyse_unit(m, r, drivers[s], k_in.get(s, np.nan), k_m2.get((s, m["module"]), np.nan),
                            noon.get(s, np.nan))
        u.update(info[s])
        unit_rows.append(u)
        drow += d
    write_csv(unit_rows, Q16 / "q16b_units.csv")
    write_csv(drow, Q16 / "q16b_unit_drivers.csv")

    # network summary per driver
    summ = []
    net = [x for x in unit_rows if x["stream"] != PAYERNE]
    for lab in ("in", "m2"):
        a_ = med_iqr([x.get(f"amp_r{lab}") for x in net])
        print(f"residual {lab}: amp {a_[0]:.2f} [{a_[1]:.2f},{a_[2]:.2f}] %, peak {circ_median_doy([x.get(f'peak_r{lab}') for x in net]):.0f}")
    for name in DRIVERS:
        dd = [d for d in drow if d["driver"] == name and d["stream"] != PAYERNE]
        if len(dd) < 5:
            continue
        row = dict(driver=name, family=DRIVERS[name][0], unit=DRIVERS[name][1], n_units=len(dd))
        row["amp_D"], row["amp_D_q25"], row["amp_D_q75"] = med_iqr([d["amp_D"] for d in dd])
        row["peak_D"] = circ_median_doy([d["peak_D"] for d in dd])
        for lab in ("in", "m2"):
            row[f"dphase_{lab}"] = med_iqr([d.get(f"dphase_{lab}") for d in dd])[0]
            row[f"cos_{lab}"] = med_iqr([d.get(f"cos_{lab}") for d in dd])[0]
            row[f"gamma_{lab}"] = med_iqr([d.get(f"gamma_{lab}") for d in dd])[0]
            row[f"amp_before_{lab}"] = med_iqr([d.get(f"amp_before_{lab}") for d in dd])[0]
            row[f"amp_after_{lab}"] = med_iqr([d.get(f"amp_after_{lab}") for d in dd])[0]
            row[f"scat_after_minus_before_{lab}"] = med_iqr([d.get(f"scat_after_{lab}", np.nan) - d.get(f"scat_before_{lab}", np.nan) for d in dd])[0]
            row[f"transfer_{lab}"] = med_iqr([d.get(f"transfer_{lab}") for d in dd])[0]
            k0 = np.array([d.get(f"oos_mon_k_{lab}", np.nan) for d in dd])
            k1 = np.array([d.get(f"oos_mon_kD_{lab}", np.nan) for d in dd])
            kh = np.array([d.get(f"oos_mon_kH_{lab}", np.nan) for d in dd])
            ok = np.isfinite(k0) & np.isfinite(k1) & np.isfinite(kh)
            row[f"oos_n_{lab}"] = int(ok.sum())
            row[f"oos_mon_k_{lab}"] = float(np.median(k0[ok])) if ok.any() else np.nan
            row[f"oos_mon_kD_{lab}"] = float(np.median(k1[ok])) if ok.any() else np.nan
            row[f"oos_mon_kH_{lab}"] = float(np.median(kh[ok])) if ok.any() else np.nan
            row[f"oos_red_{lab}"] = float(np.median(1 - k1[ok] / k0[ok])) if ok.any() else np.nan
            row[f"oos_red_H_{lab}"] = float(np.median(1 - kh[ok] / k0[ok])) if ok.any() else np.nan
            row[f"oos_frac_better_{lab}"] = float(np.mean(k1[ok] < k0[ok])) if ok.any() else np.nan
            # night-level variance reduction on the 2026 test set (less noisy than monthly means):
            # share of the calendar-predictable variance (harmonic fitted on 2025) the driver captures
            s0 = np.array([d.get(f"oos_sd_k_{lab}", np.nan) for d in dd]) ** 2
            s1 = np.array([d.get(f"oos_sd_kD_{lab}", np.nan) for d in dd]) ** 2
            sh = np.array([d.get(f"oos_sd_kH_{lab}", np.nan) for d in dd]) ** 2
            ok = np.isfinite(s0) & np.isfinite(s1) & np.isfinite(sh)
            row[f"oos_dvar_D_{lab}"] = float(np.sum(s0[ok] - s1[ok]) / ok.sum()) if ok.any() else np.nan
            row[f"oos_dvar_H_{lab}"] = float(np.sum(s0[ok] - sh[ok]) / ok.sum()) if ok.any() else np.nan
            row[f"oos_capture_{lab}"] = row[f"oos_dvar_D_{lab}"] / row[f"oos_dvar_H_{lab}"] if ok.any() else np.nan
            # annual gain the driver would need to carry the in-phase residual (sign included)
            need = []
            for d in dd:
                u_ = next((x for x in net if x["stream"] == d["stream"]), None)
                if u_ is None or not np.isfinite(u_.get(f"p_r{lab}", np.nan)):
                    continue
                pr_ = u_[f"p_r{lab}"]
                need.append(100 * abs(pr_) * np.cos(W * wrap_doy(d["peak_D"] - peak_doy(pr_))) / d["amp_D"])
            row[f"gamma_needed_{lab}"] = float(np.median(need)) if need else np.nan
        row["g_syn_sd"], row["g_syn_sd_q25"], row["g_syn_sd_q75"] = med_iqr([d.get("g_syn_sd") for d in dd])
        row["pcorr"] = med_iqr([d.get("pcorr") for d in dd])[0]
        tt = np.array([d.get("g_syn_sd", np.nan) / d.get("se_syn_sd", np.nan) for d in dd])
        row["frac_t_gt2"] = float(np.nanmean(np.abs(tt[np.isfinite(tt)]) > 2)) if np.isfinite(tt).any() else np.nan
        row["kappa_raw"] = med_iqr([d.get("kappa_raw") for d in dd])[0]
        row["kappa_joint"] = med_iqr([d.get("kappa_joint") for d in dd])[0]
        row["amp_after_free"] = med_iqr([d.get("amp_after_free") for d in dd])[0]
        row["day_excess_kappa_m2"] = med_iqr([d.get("day_excess_kappa_m2") for d in dd])[0]
        row["day_excess_joint"] = med_iqr([d.get("day_excess_joint") for d in dd])[0]
        # latitude contrast: does the residual's in-phase amplitude scale with the driver amplitude?
        pr = {x["stream"]: x.get("p_rm2", np.nan) for x in net}
        aD = np.array([d["amp_D"] for d in dd])
        ip = np.array([100 * np.real(pr.get(d["stream"], np.nan) * np.exp(-1j * W * d["peak_D"])) for d in dd])
        ok = np.isfinite(aD) & np.isfinite(ip)
        if ok.sum() > 10:
            row["rho_amp"] = float(spearmanr(aD[ok], ip[ok])[0])
            c, se, _ = ols(np.column_stack([np.ones(ok.sum()), aD[ok]]), ip[ok])
            row["xsec_slope"], row["xsec_slope_se"], row["xsec_icpt"] = c[1], se[1], c[0]
        row.update(pooled_synoptic(units, drivers, name))
        row.update(pooled_gamma(units, drivers, k_m2, name, "m2"))
        summ.append(row)
        print(f"{name:18s} peak {row['peak_D']:5.0f} cos {row['cos_m2']:+.2f} amp {row['amp_before_m2']:.1f}->{row['amp_after_m2']:.1f} "
              f"oos {row['oos_red_m2']:+.2f} (H {row['oos_red_H_m2']:+.2f}) gsyn {row['g_syn_sd']:+.2f} pc {row['pcorr']:+.3f} "
              f"pool {row.get('g_pool', np.nan):+.3g}+-{row.get('g_pool_se', np.nan):.2g} kj {row['kappa_joint']:+.2f} "
              f"rho {row.get('rho_amp', np.nan):+.2f} need {row['gamma_needed_m2']:+.3g} capture {row['oos_capture_m2']:+.2f} "
              f"(dvarH {row['oos_dvar_H_m2']:.1f})", flush=True)
    write_csv(summ, Q16 / "q16b_drivers.csv")
    # latitude contrast in bins: residual in-phase amplitude vs the day-length / TOA amplitude
    dl = {d["stream"]: d["amp_D"] for d in drow if d["driver"] == "daylen_h"}
    toa = {d["stream"]: d["amp_D"] for d in drow if d["driver"] == "toa_prev"}
    lrows = []
    for lo, hi in ((0, 45), (45, 50), (50, 53), (53, 57), (57, 90)):
        sel = [x for x in net if lo <= x["lat"] < hi and np.isfinite(x.get("p_rm2", np.nan))]
        if not sel:
            continue
        for lab in ("in", "m2"):
            ip = [100 * np.real(x[f"p_r{lab}"] * np.exp(-1j * W * 172.0)) for x in sel if np.isfinite(x.get(f"p_r{lab}", np.nan))]
            lrows.append(dict(lat_bin=f"{lo}-{hi}", variant=lab, n=len(ip), resid_inphase172=float(np.median(ip)),
                              resid_q25=float(np.percentile(ip, 25)), resid_q75=float(np.percentile(ip, 75)),
                              daylen_amp_h=float(np.median([dl[x["stream"]] for x in sel if x["stream"] in dl])),
                              toa_amp=float(np.median([toa[x["stream"]] for x in sel if x["stream"] in toa])),
                              amp_C=float(np.median([x["amp_C"] for x in sel])),
                              T_amp=float(np.median([x["amp_T"] for x in sel]))))
            print("lat", lrows[-1])
    write_csv(lrows, Q16 / "q16b_latitude_bins.csv")

    prow = pairs_test(units, drivers, k_m2)
    write_csv(prow, Q16 / "q16b_pairs.csv")
    t4, t4u = test4(units, drivers)
    write_csv(t4, Q16 / "q16b_q15e_kappa.csv")
    write_csv(t4u, Q16 / "q16b_q15e_units.csv")
    for r in t4:
        if r["band"] in ("win", "lnc", "7.0-8.0", "2.5-3.0", "4.0-5.0"):
            print(f"t4 {r['spec']:28s} {r['band']:8s} kappa {r['kappa']:+.2f} +- {r['se']:.2f} "
                  f"gamma/sd {r.get('gamma_per_sd', np.nan):+.2f} +- {r.get('gamma_per_sd_se', np.nan):.2f} ({r['n_units']})")
    make_figures(units, drivers, k_m2, unit_rows, drow, t4)
    solstice_curvature(units, drivers)
    mirror_test(units, drivers)
    interannual_test(units, drivers)
    # Bern measured radiation vs the ceilometer proxy
    units_by = {m["stream"]: r for m, r in units}
    for s in SMN:
        if s in drivers:
            x, y = drivers[s]["smn_gre_prev"], drivers[s]["ghi_proxy_prev"]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() > 20:
                tt = units_by[s]["t"]
                xa, ya = anomaly(tt, x), anomaly(tt, y)
                oka = np.isfinite(xa) & np.isfinite(ya)
                print(f"{s}: measured vs proxy insolation r = {np.corrcoef(x[ok], y[ok])[0, 1]:.3f} "
                      f"(31-d anomalies r = {np.corrcoef(xa[oka], ya[oka])[0, 1]:.3f}), n = {ok.sum()}, "
                      f"mean measured {np.mean(x[ok]):.2f} proxy {np.mean(y[ok]):.2f} MJ/m2/d")
    # Payerne separately
    pay = [d for d in drow if d["stream"] == PAYERNE]
    for d in pay:
        if d["driver"] in ("daylen_h", "csky_prev", "lom_ew7d", "Tint_daymax"):
            print("Payerne", d["driver"], {k: round(v, 3) for k, v in d.items() if isinstance(v, float) and k in
                                           ("cos_m2", "amp_before_m2", "amp_after_m2", "g_syn_sd", "pcorr")})


if __name__ == "__main__":
    main()
