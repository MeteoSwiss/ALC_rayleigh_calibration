# -*- coding: utf-8 -*-
"""q16c - molecular detectability and noise in the CHM15k Rayleigh calibration.

Operator question: "Can it be that the instrument sees the molecular signal only up to a certain
altitude, and above that it is only noise?"  If so, the Rayleigh fit (and the 7-8 km band used in
q15e) could be fitting noise or background residuals, and a noise-driven bias could follow the
season (short summer nights = fewer profiles, more noise).

Stages (run: .venv/Scripts/python.exe validation/chm15k_stability/q16c_detectability.py <stage>)
  l1       native L1 on ~26 units x ~40 calibration nights (operational window t0..t1): per-gate noise
           measured by temporal first differences, binned to the operational 300 s x 30 m grid, the
           night-mean molecular SNR profile, the SNR=1/3/5 heights, in-window SNR, and the profile
           shape X(z) = S/(C p_mol) to 15 km (does it follow the molecular slope or flatten?).
  sim      Monte-Carlo of the OPERATIONAL estimator (calibration.rayleigh: screen + mean, eprof_v2.2
           window search with temporal stack and measured sigma, 125-member Klett/LR/shift/time
           ensemble, median, QC -2/-3/-6/-7/-8/-9) on synthetic molecular profiles with the noise of
           every sampled night (and a controlled noise-multiplier scan) -> bias of C vs SNR.
  net      network-wide prediction (all CHM15k units, q15f nightly table): noise transfer
           sd -> night-mean sigma(z) calibrated on the L1 stage -> SNR heights, in-window SNR,
           predicted noise-only u/C and the predicted noise bias per night; seasonal statistics.
  scatter  does the night-to-night scatter of ln C scale with the predicted noise?
  summary  figures + q16c_summary.md

Nothing in calibration/ is modified: the operational functions are imported and called as-is.
Outputs: C:/DATA/Projects/202610_CHM15k_stability/q16/q16c_*
"""
from __future__ import annotations

import csv
import json
import sys
import warnings
from datetime import datetime, timedelta
from multiprocessing import Pool
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (CENSUS, L1_ROOT, OUT, chm15k_streams, load_constants,  # noqa: E402
                    stream_key)

warnings.filterwarnings("ignore")
Q16 = OUT / "q16"
Q15 = OUT / "q15"
PAYERNE = "0-20000-0-06610_A"
N_WORKERS = 4
WAVELENGTH = 1064e-9
W_YEAR = 2 * np.pi / 365.25
N_NIGHTS = 40
N_UNITS = 26

# Operational settings (options.json + config.py defaults + calibration.py constants)
OP = dict(range_start_m=2000.0, range_end_m=6000.0, increment_bins=8,
          half_length_options_m=tuple(range(250, 2000, 240)), lr=52.0, threshold_quality=15.0,
          aerosol_scattering_threshold=2.0, nmad=4.0, z_low_cloud=6000.0, method="eprof_v2.2",
          block_s=300.0, gates_per_bin=3)
LR_DELTAS = (-20, -10, 0, +10, +20)
ALT_SHIFTS_M = (-200, -100, 0, +100, +200)
N_TIME_SAMPLES = 4
TIME_SUBSET_FRACTION = 0.7

SEASONS = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}
SNR_LEVELS = (1.0, 3.0, 5.0)
X_EDGES = np.arange(1000.0, 15001.0, 250.0)
X_MID = 0.5 * (X_EDGES[:-1] + X_EDGES[1:])


# ---------------------------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------------------------
def season_of(t_days):
    m = (datetime(1970, 1, 1) + timedelta(days=float(t_days))).month
    for k, v in SEASONS.items():
        if m in v:
            return k
    return "?"


def doy_of(t_days):
    return (datetime(1970, 1, 1) + timedelta(days=float(t_days))).timetuple().tm_yday


def station_meta():
    out = {}
    for s in json.load(open(CENSUS, encoding="utf-8")):
        out[stream_key(s["wmo"], s["ident"])] = dict(alt=float(s.get("alt", 0.0) or 0.0),
                                                     lat=s["lat"], lon=s["lon"], site=s["site"])
    return out


def molecular(alt, rng):
    """beta_mol, p_mol (= beta_mol T_mol^2 / r^2) on the range grid, as the operational code does
    (US-1976 standard atmosphere on altitude + range, calculate_molecular_properties)."""
    from calibration.rayleigh.atmosphere import calculate_molecular_properties, load_standard_atmosphere
    atm = load_standard_atmosphere(None, alt + rng)
    mp = calculate_molecular_properties(atm.temperature, atm.pressure, rng, WAVELENGTH)
    return mp.beta_mol, mp.p_mol


def annual_fit(t, y):
    """Offset + annual harmonic (+ trend when > 400 d). Returns amplitude, peak DOY."""
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = np.asarray(t)[ok], np.asarray(y)[ok]
    if len(y) < 20 or t.max() - t.min() < 250:
        return np.nan, np.nan
    ph = W_YEAR * (t % 365.25)
    cols = [np.ones_like(t), np.sin(ph), np.cos(ph)]
    if t.max() - t.min() > 400:
        cols.append((t - t.mean()) / 365.25)
    c, *_ = np.linalg.lstsq(np.column_stack(cols), y, rcond=None)
    return float(np.hypot(c[1], c[2])), float((np.arctan2(c[1], c[2]) / W_YEAR) % 365.25 + 1.0)


def first_cross(z, snr, level, zmin=1000.0):
    """Lowest altitude above zmin where the (smoothed) SNR falls below level (stays below 150 m)."""
    ok = z >= zmin
    zz, ss = z[ok], snr[ok]
    below = ss < level
    for i in range(len(zz) - 5):
        if below[i] and below[i:i + 5].all():
            return float(zz[i])
    return np.nan


def run_med(x, k=9):
    """Running median over k gates (NaN-aware) - used to smooth per-gate sigma estimates."""
    from scipy.ndimage import median_filter
    y = np.where(np.isfinite(x), x, np.nanmedian(x))
    return median_filter(y, size=k, mode="nearest")


# ---------------------------------------------------------------------------------------------
# Emulation of the operational Rayleigh estimator (calibration/rayleigh/calibration.py steps 6-10)
# ---------------------------------------------------------------------------------------------
def emulate_operational(stack_rcs, rng, beta_mol, p_mol, sigma_bin, seed=42, method=None,
                        fixed_window=None):
    """Run the operational estimator on a (n_blocks x n_range) rcs stack on the 30 m fit grid.

    sigma_bin = 1-sigma noise of ONE binned block profile on the signal (rcs/r^2) scale; the
    operational code passes sigma/sqrt(n_profiles) of the night mean to the noise-aware gates.
    Returns a dict with the final constant (median of the perturbation ensemble), the nominal
    member, simpler variants (window ratio-mean), the QC flag and the chosen window.
    """
    from calibration.rayleigh.calibration import _compute_cl_for_perturbation
    from calibration.rayleigh.rayleigh_fit import find_optimal_molecular_window
    from calibration.rayleigh.atmosphere import klett_inversion, MOLECULAR_LIDAR_RATIO

    method = method or OP["method"]
    out = dict(flag=0, C=np.nan, C_nom=np.nan, C_ratio=np.nan, C_slope=np.nan, zb=np.nan, zt=np.nan,
               n_use=0)
    r2 = rng ** 2
    # (6a) profile outlier screen on the 2-6 km band metric
    keep = np.arange(stack_rcs.shape[0])
    if stack_rcs.shape[0] >= 10:
        band = (rng >= OP["range_start_m"]) & (rng <= OP["range_end_m"])
        metric = np.nanmedian(stack_rcs[:, band] / r2[band], axis=1)
        fin = np.isfinite(metric)
        if fin.sum() >= 10:
            cen = np.median(metric[fin])
            mad = np.median(np.abs(metric[fin] - cen))
            thr = OP["nmad"] * 1.4826 * mad
            k = fin & (np.abs(metric - cen) <= max(thr, 1e-30))
            if k.sum() >= max(3, int(0.5 * stack_rcs.shape[0])):
                keep = np.where(k)[0]
    rcs_use = stack_rcs[keep]
    out["n_use"] = rcs_use.shape[0]
    rcs_mean = np.nanmean(rcs_use, axis=0)
    signal = rcs_mean / r2
    signal_stack = rcs_use / r2[None, :]
    sigma_signal = None if sigma_bin is None else sigma_bin / np.sqrt(max(rcs_use.shape[0], 1))

    if fixed_window is not None:
        w = (rng >= fixed_window[0]) & (rng <= fixed_window[1])
        out["C_fixed"] = float(np.nanmean(signal[w] / p_mol[w]))

    fit = find_optimal_molecular_window(
        signal=signal, p_mol=p_mol, range_alc=rng,
        half_length_options_m=OP["half_length_options_m"], range_start_m=OP["range_start_m"],
        range_end_m=OP["range_end_m"], increment_bins=OP["increment_bins"], method=method,
        signal_stack=signal_stack, sigma_signal=sigma_signal)
    if not np.isfinite(fit.slope) or not np.isfinite(fit.relative_error):
        out["flag"] = -2
        return out
    if not fit.is_valid:
        out["flag"] = -7 if fit.slope <= 0 else -8
        return out
    out["zb"], out["zt"] = fit.range_start_m, fit.range_end_m
    w = (rng >= fit.range_start_m) & (rng <= fit.range_end_m)
    out["C_ratio"] = float(np.nanmean(signal[w] / p_mol[w]))
    out["r2"] = float(fit.r_squared)
    if fit.relative_error > OP["threshold_quality"]:
        out["flag"] = -2
        return out

    # (8-10) perturbation ensemble: 5 time samples x 5 LR x 5 shifts, median
    rs = np.random.default_rng(seed)
    n_p = rcs_use.shape[0]
    n_keep = min(n_p, max(3, int(n_p * TIME_SUBSET_FRACTION)))
    samples = [rcs_mean]
    for _ in range(N_TIME_SAMPLES):
        idx = np.sort(rs.choice(n_p, size=n_keep, replace=False))
        samples.append(np.nanmean(rcs_use[idx], axis=0))
    cl = []
    for rm in samples:
        for d_lr in LR_DELTAS:
            for sh in ALT_SHIFTS_M:
                pr = _compute_cl_for_perturbation(
                    rcs_mean=rm, range_alc=rng, beta_mol=beta_mol, fit_result=fit,
                    lidar_ratio_aerosol=OP["lr"] + d_lr, altitude_shift_m=sh,
                    subtract_background=False, consider_points_lower_than_molecular=True)
                if pr is not None:
                    cl.append(pr.lidar_constant)
                    if rm is rcs_mean and d_lr == 0 and sh == 0:
                        out["C_nom"] = float(pr.lidar_constant)
    if not cl:
        out["flag"] = -5
        return out
    cl = np.array(cl)
    c_med = float(np.median(cl))
    q25, q75 = np.percentile(cl, [25, 75])
    mad = float(np.median(np.abs(cl - c_med)))
    unc = 2 * max((q75 - q25) / 1.349, mad / 0.6745)
    out["C"], out["u"] = c_med, unc

    # QC -3 (slope method), -6 (uncertainty), -9 (lower-signal layer)
    sig_n = signal / fit.slope
    batt = sig_n * r2
    ii = np.where(w & np.isfinite(rcs_mean))[0]
    ref_idx = int((ii[0] + ii[-1]) / 2)
    ref_val = np.nanmean(batt[w] / beta_mol[w])
    try:
        _, _, ext_aer = klett_inversion(beta_att=batt, beta_mol=beta_mol, range_alc=rng,
                                        reference_index=ref_idx, lidar_ratio_aerosol=OP["lr"],
                                        reference_value=ref_val, i_start=0, i_end=ii[-1],
                                        sign_error_v10=False)
        od = np.trapezoid(ext_aer[:ii[0]], rng[:ii[0]])
        cl_slope = fit.slope * np.exp(2 * od)
    except Exception:
        cl_slope = np.nan
    out["C_slope"] = float(cl_slope)
    err = abs((cl_slope - c_med) / c_med * 100) if np.isfinite(cl_slope) else np.inf
    if err > OP["threshold_quality"]:
        out["flag"] = -3
    elif unc > c_med:
        out["flag"] = -6
    else:
        d = fit.search_diagnostics
        if d is not None:
            slp, r2g = np.asarray(d.slopes, float), np.asarray(d.r_squared, float)
            clean = np.isfinite(slp) & (slp > 0) & (r2g >= 0.5)
            if clean.sum() >= 3:
                ci = int(np.argmin(np.abs(np.asarray(d.range_bin_m) - fit.center_range_m)))
                hj = int(np.argmin(np.abs(np.asarray(d.half_length_m) - fit.half_length_m)))
                c_ref = float(np.nanpercentile(np.where(clean, slp, np.nan), 10))
                if c_ref > 0 and slp[ci, hj] / c_ref > OP["aerosol_scattering_threshold"]:
                    out["flag"] = -9
    return out


# ---------------------------------------------------------------------------------------------
# Stage l1: native L1, sampled calibration nights
# ---------------------------------------------------------------------------------------------
def l1_path(wmo, ident, day):
    d = datetime(1970, 1, 1) + timedelta(days=int(day))
    return L1_ROOT / wmo / f"{d:%Y}" / f"{d:%m}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"


def read_night(wmo, ident, t0, t1):
    import netCDF4
    parts = []
    meta = None
    for day in range(int(np.floor(t0)), int(np.floor(t1)) + 1):
        f = l1_path(wmo, ident, day)
        if not f.exists():
            continue
        try:
            with netCDF4.Dataset(f) as ds:
                t = np.ma.filled(ds["time"][:].astype("f8"), np.nan)
                i0, i1 = np.searchsorted(t, t0), np.searchsorted(t, t1, side="right")
                if i1 <= i0:
                    continue
                if meta is None:
                    meta = dict(rng=np.ma.filled(ds["range"][:].astype("f8"), np.nan),
                                alt=float(np.ma.filled(ds["station_altitude"][...], np.nan)),
                                tilt=float(np.ma.filled(ds["tilt_angle"][...], 0.0)))
                parts.append(dict(
                    t=t[i0:i1],
                    rcs=np.ma.filled(ds["rcs_0"][i0:i1, :].astype("f8"), np.nan),
                    cbh=np.ma.filled(ds["cloud_base_height"][i0:i1, 0].astype("f8"), np.nan),
                    sd=np.ma.filled(ds["stddev"][i0:i1].astype("f8"), np.nan),
                    bg=np.ma.filled(ds["bckgrd_rcs_0"][i0:i1].astype("f8"), np.nan)))
        except (OSError, RuntimeError, KeyError, IndexError):
            continue
    if not parts or meta is None:
        return None
    d = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    d.update(meta)
    return d


def process_night(wmo, ident, alt_census, t0, t1, c_op, zb, zt):
    """All per-night SNR / detectability diagnostics on the operational grid."""
    from calibration.sensitivity.noise import first_difference_sigma
    from calibration.rayleigh.molecular_methods import flag_contaminated_cells
    d = read_night(wmo, ident, t0, t1)
    if d is None or d["rcs"].shape[0] < 100:
        return None
    alt = d["alt"] if np.isfinite(d["alt"]) else alt_census
    rng_n = d["rng"]
    t_s = (d["t"] - t0) * 86400.0
    dt = float(np.median(np.diff(t_s)))
    # low-cloud screen as filter_cloudy_profiles: cloud base <= z_low_cloud, dilated by 15 min
    low = np.isfinite(d["cbh"]) & (d["cbh"] > 0) & (d["cbh"] <= OP["z_low_cloud"])
    from scipy.ndimage import maximum_filter1d
    win = max(1, int(round(15 * 60 / dt)))
    contam = maximum_filter1d(low.astype(np.uint8), size=2 * win - 1, mode="constant") > 0
    keep = ~contam & np.all(np.isfinite(d["rcs"][:, :50]), axis=1)
    n_nat = int(keep.sum())
    if n_nat < 100:
        return None
    sig_nat = d["rcs"] / rng_n[None, :] ** 2
    sigma_nat, _ = first_difference_sigma(sig_nat, t_s, dt, sel=keep)
    # bin to the operational 300 s x 30 m grid
    g = OP["gates_per_bin"]
    nz = len(rng_n) // g
    rng = rng_n[:nz * g].reshape(nz, g).mean(axis=1)
    blk = np.floor(t_s / OP["block_s"]).astype(int)
    rows, tb = [], []
    for b in np.unique(blk[keep]):
        s = keep & (blk == b)
        if s.sum() < 0.5 * OP["block_s"] / dt:
            continue
        rows.append(np.nanmean(d["rcs"][s, :nz * g], axis=0).reshape(nz, g).mean(axis=1))
        tb.append(b * OP["block_s"])
    if len(rows) < 10:
        return None
    stack = np.array(rows)
    tb = np.array(tb, float)
    beta_mol, p_mol = molecular(alt, rng)
    r2 = rng ** 2
    # operational profile screen (same as emulate_operational step 6a)
    band = (rng >= OP["range_start_m"]) & (rng <= OP["range_end_m"])
    metric = np.nanmedian(stack[:, band] / r2[band], axis=1)
    cen = np.median(metric)
    mad = np.median(np.abs(metric - cen))
    ks = np.abs(metric - cen) <= max(OP["nmad"] * 1.4826 * mad, 1e-30)
    if ks.sum() < max(3, int(0.5 * len(metric))):
        ks[:] = True
    stack, tb = stack[ks], tb[ks]
    nb = stack.shape[0]
    S_stack = stack / r2[None, :]
    flag = flag_contaminated_cells(S_stack, p_mol, rng, nmad=4.0, min_excess=0.25)
    S_clean = np.nanmean(np.where(flag, np.nan, S_stack), axis=0)
    n_eff = np.maximum((~flag).sum(axis=0), 1)
    S_mean = np.nanmean(S_stack, axis=0)
    # noise: empirical on the binned stack, and propagated from the native gates (white-noise test)
    sigma_bin, _ = first_difference_sigma(S_stack, tb, OP["block_s"])
    sn = sigma_nat[:nz * g].reshape(nz, g)
    n_per_blk = OP["block_s"] / dt
    sigma_prop = np.sqrt(np.nanmean(sn ** 2, axis=1) / g) / np.sqrt(n_per_blk)
    sb_s = run_med(sigma_bin)
    sigma_mean = sb_s / np.sqrt(n_eff)
    # window (operational, AGL range) and the night constant in L1 units
    zb_r, zt_r = zb - alt, zt - alt
    w = (rng >= zb_r) & (rng <= zt_r)
    if w.sum() < 5:
        return None
    C = float(np.nanmean(S_mean[w] / p_mol[w]))
    if not C > 0:
        return None
    snr = C * p_mol / sigma_mean                       # molecular SNR of the night mean, 30 m gate
    res = dict(stream=stream_key(wmo, ident), t=0.5 * (t0 + t1), doy=doy_of(0.5 * (t0 + t1)),
               season=season_of(0.5 * (t0 + t1)), alt=alt, tilt=d["tilt"], n_nat=n_nat, nb=nb,
               dur_h=(t1 - t0) * 24, c_op=c_op, c_l1=C, c_ratio=C / c_op, zb=zb_r, zt=zt_r,
               sd=float(np.nanmedian(d["sd"][keep])), bg=float(np.nanmedian(d["bg"][keep])))
    for lev in SNR_LEVELS:
        res[f"z_snr{int(lev)}"] = first_cross(rng, snr, lev)
        res[f"z250_snr{int(lev)}"] = first_cross(rng, snr * np.sqrt(250.0 / 29.97), lev)
    iw = np.where(w)[0]
    res["snr_zb"] = float(np.nanmedian(snr[iw[:5]]))
    res["snr_zt"] = float(np.nanmedian(snr[iw[-5:]]))
    res["frac_win_snr3"] = float(np.mean(snr[w] < 3))
    res["frac_win_snr1"] = float(np.mean(snr[w] < 1))
    res["u_noise"] = float(np.sqrt(np.nansum((sigma_mean[w] / (C * p_mol[w])) ** 2)) / w.sum())
    # noise at fixed bands, white-noise check and transfer from the HK stddev
    for a, b in ((3000, 4000), (5000, 6000), (7000, 8000), (12000, 15000)):
        s = (rng >= a) & (rng < b)
        res[f"snr_{a // 1000}_{b // 1000}"] = float(np.nanmedian(snr[s]))
        res[f"sbin_{a // 1000}_{b // 1000}"] = float(np.nanmedian(sigma_bin[s]))
        res[f"wn_{a // 1000}_{b // 1000}"] = float(np.nanmedian(sigma_bin[s] / sigma_prop[s]))
    # rho(z): night-mean 30 m sigma * sqrt(3 n_nat) / sd  (network transfer function)
    rho = sigma_mean * np.sqrt(g * n_nat) / res["sd"]
    # profile shape X(z) = S/(C p_mol) on 250 m bins to 15 km (clean mean and plain mean)
    gb = np.digitize(rng, X_EDGES) - 1
    X = np.full(len(X_MID), np.nan)
    Xm = np.full(len(X_MID), np.nan)
    sX = np.full(len(X_MID), np.nan)
    for k in range(len(X_MID)):
        s = gb == k
        if s.sum():
            X[k] = np.nanmean(S_clean[s] / (C * p_mol[s]))
            Xm[k] = np.nanmean(S_mean[s] / (C * p_mol[s]))
            sX[k] = np.sqrt(np.nansum((sigma_mean[s] / (C * p_mol[s])) ** 2)) / s.sum()
    # molecular-slope tests: weighted S = A p_mol + B  (i) in the window, (ii) in 7-8 km,
    # (iii) from the window bottom to 15 km (offset floor / residual baseline)
    for lab, sel in (("win", w), ("78", (rng >= 7000) & (rng < 8000)),
                     ("long", (rng >= zb_r) & (rng <= 15000))):
        x, y, s = p_mol[sel], S_mean[sel], sigma_mean[sel]
        ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(s) & (s > 0)
        if ok.sum() < 8:
            continue
        Wt = 1 / s[ok] ** 2
        Xd = np.column_stack([x[ok], np.ones(ok.sum())])
        cov = np.linalg.inv(Xd.T @ (Xd * Wt[:, None]))
        coef = cov @ (Xd.T @ (Wt * y[ok]))
        resid = (y[ok] - Xd @ coef) / s[ok]
        chi2 = float(np.sum(resid ** 2) / max(ok.sum() - 2, 1))
        A, B = coef
        sA, sB = np.sqrt(np.diag(cov) * max(chi2, 1.0))
        pm = np.nanmean(x[ok])
        res[f"slopeA_{lab}"] = A / C            # 1 = molecular slope matches the window constant
        res[f"slopeA_{lab}_se"] = sA / C
        res[f"offB_{lab}"] = B / (C * pm)       # offset as a fraction of the molecular signal in band
        res[f"offB_{lab}_se"] = sB / (C * pm)
        res[f"chi2_{lab}"] = chi2
        if lab == "long":
            res["offB_long_at_zt"] = B / (C * p_mol[iw[-1]])
            res["offB_long_at_75"] = B / (C * np.interp(7500.0, rng, p_mol))
            res["offB_long_win_bias"] = float(np.nanmean(B / (C * p_mol[w])))
    arrays = dict(snr=snr.astype("f4"), sig_rel=(sb_s / C).astype("f4"), rho=rho.astype("f4"),
                  X=X.astype("f4"), Xm=Xm.astype("f4"), sX=sX.astype("f4"))
    return res, arrays, rng


def pick_units():
    """Geographically spread CHM15k sample (one unit per 4-degree cell) + Payerne (reported apart)."""
    from common import load_cache
    streams = chm15k_streams()
    cells, out = {}, []
    rs = np.random.default_rng(7)
    order = rs.permutation(len(streams))
    for i in order:
        w, ident, s, la, lo = streams[i]
        key = stream_key(w, ident)
        if key == PAYERNE:
            continue
        c = (int(la // 4), int(lo // 4))
        if c in cells or load_cache(w, ident) is None:
            continue
        cons = load_constants(w, ident)
        if cons is None:
            continue
        mod = max(set(cons["module"]), key=list(cons["module"]).count)
        m = (cons["module"] == mod) & (cons["c"] > 0)
        if m.sum() < 120:
            continue
        cells[c] = key
        out.append((w, ident))
    out = sorted(out)[:N_UNITS]
    out.append(tuple(PAYERNE.rsplit("_", 1)))
    return out


def l1_unit(args):
    wmo, ident = args
    meta = station_meta().get(stream_key(wmo, ident), dict(alt=0.0))
    cons = load_constants(wmo, ident)
    mod = max(set(cons["module"]), key=list(cons["module"]).count)
    m = np.where((cons["module"] == mod) & (cons["c"] > 0) & np.isfinite(cons["t0"]))[0]
    months = np.array([(datetime(1970, 1, 1) + timedelta(days=float(x))).month for x in cons["t"][m]])
    pick = []
    per_month = int(np.ceil(N_NIGHTS / 12))
    for mo in range(1, 13):
        mm = m[months == mo]
        if len(mm):
            pick += list(mm[np.linspace(0, len(mm) - 1, min(per_month, len(mm))).astype(int)])
    pick = sorted(set(pick))
    if len(pick) > N_NIGHTS:
        pick = [pick[i] for i in np.linspace(0, len(pick) - 1, N_NIGHTS).astype(int)]
    rows, arrs, rng = [], [], None
    for k in pick:
        try:
            r = process_night(wmo, ident, meta["alt"], cons["t0"][k], cons["t1"][k], cons["c"][k],
                              cons["zb"][k], cons["zt"][k])
        except Exception as e:  # keep going: one broken file must not kill the unit
            print(f"  {wmo}_{ident} night {k}: {e}", flush=True)
            r = None
        if r is None:
            continue
        rows.append(r[0])
        arrs.append(r[1])
        rng = r[2]
    print(f"{wmo}_{ident}: {len(rows)} nights", flush=True)
    return rows, arrs, rng


def stage_l1():
    units = pick_units()
    print(f"l1 stage: {len(units)} units", flush=True)
    with Pool(N_WORKERS) as pool:
        res = pool.map(l1_unit, units, chunksize=1)
    rows = [r for rr in res for r in rr[0]]
    arrs = [a for rr in res for a in rr[1]]
    rng = next(rr[2] for rr in res if rr[2] is not None)
    nz = min(len(a["snr"]) for a in arrs)
    keys = list(rows[0].keys())
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(Q16 / "q16c_l1_nights.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    np.savez_compressed(Q16 / "q16c_l1_profiles.npz", rng=rng[:nz], x_mid=X_MID,
                        **{k: np.array([a[k][:nz] if a[k].size > len(X_MID) else a[k] for a in arrs])
                           for k in arrs[0]})
    print(f"l1 stage: {len(rows)} nights -> q16c_l1_nights.csv", flush=True)


def load_l1():
    import pandas as pd
    df = pd.read_csv(Q16 / "q16c_l1_nights.csv")
    z = np.load(Q16 / "q16c_l1_profiles.npz")
    return df, {k: z[k] for k in z.files}


# ---------------------------------------------------------------------------------------------
# Stage sim: Monte-Carlo of the operational estimator
# ---------------------------------------------------------------------------------------------
SIM_ZMAX = 9500.0


def synth_profile(alt, rng, aer=8e-7):
    """Noise-free signal (C = 1) with a boundary-layer aerosol (AOD ~0.07 at 1064 nm, LR 52;
    aer = 0 gives a purely molecular atmosphere)."""
    beta_mol, p_mol = molecular(alt, rng)
    beta_aer = aer * np.clip(1.0 - (rng - 1200.0) / 600.0, 0.0, 1.0)
    alpha_aer = 52.0 * beta_aer
    tau = np.concatenate([[0.0], np.cumsum(0.5 * (alpha_aer[1:] + alpha_aer[:-1]) * np.diff(rng))])
    sig = (beta_mol + beta_aer) * np.exp(-2 * tau) * (p_mol / beta_mol)   # p_mol/beta_mol = T_mol^2/r^2
    return beta_mol, p_mol, sig


def sim_night(args):
    """R realisations of one sampled night: its own noise profile (sigma_bin/C) and block count."""
    idx, alt, nb, sig_rel, rng_full, R, mult, seed = args[:8]
    aer = args[8] if len(args) > 8 else 8e-7
    sel = rng_full <= SIM_ZMAX
    rng = rng_full[sel]
    sig_rel = np.asarray(sig_rel, float)[sel] * mult
    beta_mol, p_mol, s0 = synth_profile(alt, rng, aer)
    rs = np.random.default_rng(seed)
    out = []
    for r in range(R):
        stack_sig = s0[None, :] + sig_rel[None, :] * rs.standard_normal((int(nb), len(rng)))
        res = emulate_operational(stack_sig * rng[None, :] ** 2, rng, beta_mol, p_mol, sig_rel,
                                  seed=42 + r, fixed_window=(3000.0, 5500.0))
        res.update(idx=idx, rep=r, mult=mult)
        out.append(res)
    return out


def stage_sim():
    df, P = load_l1()
    rng = P["rng"]
    jobs = []
    # (a) every sampled night with its own noise, R = 12
    for i, r in df.iterrows():
        if r["stream"] == PAYERNE:
            continue
        jobs.append((int(i), float(r["alt"]), int(r["nb"]), P["sig_rel"][i], rng, 12, 1.0, 1000 + int(i)))
    # (b) controlled scan around the network-median night, R = 120 per multiplier
    ok = df["stream"] != PAYERNE
    med_sig = np.nanmedian(P["sig_rel"][ok.values], axis=0)
    nb_med = int(np.median(df.loc[ok, "nb"]))
    for mult in (0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        for chunk in range(6):
            jobs.append((-1, 450.0, nb_med, med_sig, rng, 20, mult, 50000 + int(mult * 100) + chunk * 7))
    # (c) season presets: network-median noise profile and median block count per season
    for sea in SEASONS:
        s = ok & (df["season"] == sea)
        ms = np.nanmedian(P["sig_rel"][s.values], axis=0)
        for chunk in range(6):
            jobs.append((-10 - list(SEASONS).index(sea), 450.0, int(np.median(df.loc[s, "nb"])), ms, rng,
                         20, 1.0, 70000 + chunk * 11 + list(SEASONS).index(sea)))
    print(f"sim stage: {len(jobs)} jobs", flush=True)
    with Pool(N_WORKERS) as pool:
        res = pool.map(sim_night, jobs, chunksize=2)
    rows = [x for rr in res for x in rr]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(Q16 / "q16c_sim.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    # noise-free reference of the synthetic estimator (intrinsic method bias, no noise)
    beta_mol, p_mol, s0 = synth_profile(450.0, rng[rng <= SIM_ZMAX])
    rr = rng[rng <= SIM_ZMAX]
    tiny = 1e-6 * s0
    ref = emulate_operational(np.repeat((s0 * rr ** 2)[None, :], nb_med, axis=0)
                              + (tiny * rr ** 2)[None, :] * np.random.default_rng(1).standard_normal((nb_med, len(rr))),
                              rr, beta_mol, p_mol, tiny, fixed_window=(3000.0, 5500.0))
    json.dump({k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in ref.items()},
              open(Q16 / "q16c_sim_noisefree.json", "w"), indent=1)
    print(f"sim stage: {len(rows)} realisations; noise-free C = {ref['C']:.5f}", flush=True)




def stage_sim0():
    """Pure-noise control: same scan and season presets on a purely molecular atmosphere."""
    df, P = load_l1()
    rng = P["rng"]
    ok = df["stream"] != PAYERNE
    med_sig = np.nanmedian(P["sig_rel"][ok.values], axis=0)
    nb_med = int(np.median(df.loc[ok, "nb"]))
    jobs = []
    for mult in (0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        for chunk in range(6):
            jobs.append((-1, 450.0, nb_med, med_sig, rng, 20, mult, 90000 + int(mult * 100) + chunk * 7, 0.0))
    for sea in SEASONS:
        s = ok & (df["season"] == sea)
        ms = np.nanmedian(P["sig_rel"][s.values], axis=0)
        for chunk in range(6):
            jobs.append((-10 - list(SEASONS).index(sea), 450.0, int(np.median(df.loc[s, "nb"])), ms, rng,
                         20, 1.0, 95000 + chunk * 11 + list(SEASONS).index(sea), 0.0))
    with Pool(N_WORKERS) as pool:
        res = pool.map(sim_night, jobs, chunksize=2)
    rows = [x for rr in res for x in rr]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(Q16 / "q16c_sim0.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"sim0 stage: {len(rows)} realisations", flush=True)


# ---------------------------------------------------------------------------------------------
# Stage cf: strictly cloud-free profiles (no cloud layer at ANY height) -> profile shape to 15 km
# ---------------------------------------------------------------------------------------------
def pick_nights(cons):
    mod = max(set(cons["module"]), key=list(cons["module"]).count)
    m = np.where((cons["module"] == mod) & (cons["c"] > 0) & np.isfinite(cons["t0"]))[0]
    months = np.array([(datetime(1970, 1, 1) + timedelta(days=float(x))).month for x in cons["t"][m]])
    pick = []
    per_month = int(np.ceil(N_NIGHTS / 12))
    for mo in range(1, 13):
        mm = m[months == mo]
        if len(mm):
            pick += list(mm[np.linspace(0, len(mm) - 1, min(per_month, len(mm))).astype(int)])
    pick = sorted(set(pick))
    if len(pick) > N_NIGHTS:
        pick = [pick[i] for i in np.linspace(0, len(pick) - 1, N_NIGHTS).astype(int)]
    return pick


def cf_night(wmo, ident, alt_census, t0, t1, zb, zt):
    import netCDF4
    prof, cnt, n_cf, n_all, n_hi = None, None, 0, 0, 0
    sd_list, alt, rng_n = [], alt_census, None
    for day in range(int(np.floor(t0)), int(np.floor(t1)) + 1):
        f = l1_path(wmo, ident, day)
        if not f.exists():
            continue
        try:
            with netCDF4.Dataset(f) as ds:
                t = np.ma.filled(ds["time"][:].astype("f8"), np.nan)
                i0, i1 = np.searchsorted(t, t0), np.searchsorted(t, t1, side="right")
                if i1 <= i0:
                    continue
                rng_n = np.ma.filled(ds["range"][:].astype("f8"), np.nan)
                a = float(np.ma.filled(ds["station_altitude"][...], np.nan))
                alt = a if np.isfinite(a) else alt_census
                rcs = np.ma.filled(ds["rcs_0"][i0:i1, :].astype("f8"), np.nan)
                cbh = np.ma.filled(ds["cloud_base_height"][i0:i1, :].astype("f8"), np.nan)
        except (OSError, RuntimeError, KeyError, IndexError):
            continue
        anyc = np.any(np.isfinite(cbh) & (cbh > 0), axis=1)
        low = np.isfinite(cbh[:, 0]) & (cbh[:, 0] > 0) & (cbh[:, 0] <= OP["z_low_cloud"])
        from scipy.ndimage import maximum_filter1d
        lowd = maximum_filter1d(low.astype(np.uint8), size=149, mode="constant") > 0
        anyd = maximum_filter1d(anyc.astype(np.uint8), size=149, mode="constant") > 0
        op = ~lowd
        cf = ~anyd
        n_all += int(op.sum())
        n_hi += int((op & anyc).sum())
        n_cf += int(cf.sum())
        if cf.sum():
            x = rcs[cf]
            fin = np.isfinite(x)
            s_ = np.where(fin, x, 0).sum(axis=0)
            c_ = fin.sum(axis=0)
            if prof is None or len(prof) != len(s_):
                prof, cnt = s_, c_
            else:
                prof, cnt = prof + s_, cnt + c_
    if prof is None or n_cf < 300:
        return dict(n_cf=n_cf, n_op=n_all, hi_frac=n_hi / max(n_all, 1)), None
    S = prof / np.maximum(cnt, 1) / rng_n ** 2
    beta_mol, p_mol = molecular(alt, rng_n)
    w = (rng_n >= zb - alt) & (rng_n <= zt - alt)
    C = np.nanmean(S[w] / p_mol[w])
    gb = np.digitize(rng_n, X_EDGES) - 1
    X = np.array([np.nanmean(S[gb == k] / (C * p_mol[gb == k])) if np.any(gb == k) else np.nan
                  for k in range(len(X_MID))])
    return dict(n_cf=n_cf, n_op=n_all, hi_frac=n_hi / max(n_all, 1), c_cf=C), X.astype("f4")


def cf_unit(args):
    wmo, ident = args
    meta = station_meta().get(stream_key(wmo, ident), dict(alt=0.0))
    cons = load_constants(wmo, ident)
    rows, X = [], []
    for k in pick_nights(cons):
        try:
            r, x = cf_night(wmo, ident, meta["alt"], cons["t0"][k], cons["t1"][k], cons["zb"][k], cons["zt"][k])
        except Exception as e:
            print(f"  cf {wmo}_{ident} {k}: {e}", flush=True)
            continue
        r.update(stream=stream_key(wmo, ident), t=0.5 * (cons["t0"][k] + cons["t1"][k]),
                 season=season_of(0.5 * (cons["t0"][k] + cons["t1"][k])))
        rows.append(r)
        X.append(x if x is not None else np.full(len(X_MID), np.nan, "f4"))
    print(f"cf {wmo}_{ident}: {sum(np.isfinite(x).any() for x in X)}/{len(X)} cloud-free nights", flush=True)
    return rows, X


def stage_cf():
    units = pick_units()
    with Pool(N_WORKERS) as pool:
        res = pool.map(cf_unit, units, chunksize=1)
    rows = [r for rr in res for r in rr[0]]
    X = np.array([x for rr in res for x in rr[1]])
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(Q16 / "q16c_cf_nights.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    np.savez_compressed(Q16 / "q16c_cf_profiles.npz", x_mid=X_MID, X=X)
    print(f"cf stage: {len(rows)} nights", flush=True)


# ---------------------------------------------------------------------------------------------
# Post-processing helpers
# ---------------------------------------------------------------------------------------------
def nm(x):
    return float(np.nanmedian(np.asarray(x, float)))


def iqr(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (np.nan, np.nan)
    return tuple(np.percentile(x, [25, 75]))


def transfer_function(df, P):
    """rho(z) = sigma_mean(z) sqrt(3 n_nat) / sd: the HK `stddev` -> night-mean 30 m noise transfer.
    Returned as a network shape (normalised at 7-8 km) and a level per sampled unit."""
    rng = P["rng"]
    ok = (df["stream"] != PAYERNE).values
    b56 = (rng >= 5000) & (rng < 6000)
    lvl = np.nanmedian(P["rho"][:, b56], axis=1)
    shape = np.nanmedian(P["rho"][ok] / lvl[ok, None], axis=0)
    # the stored profiles stop at the shortest L1 range (~10.2 km): above, the noise is
    # background-limited and flat, so hold the last value to 15 km
    ext = np.arange(rng[-1] + 29.97, 15000.0, 29.97)
    shape = np.r_[shape, np.full(ext.size, np.nanmedian(shape[-30:]))]
    unit = df.assign(lvl=lvl).groupby("stream")["lvl"].median().to_dict()
    net = nm([v for k, v in unit.items() if k != PAYERNE])
    return shape, unit, net, lvl


def sim_tables():
    """Per-night and scan results of the Monte-Carlo, as bias relative to the noise-free estimator."""
    import pandas as pd
    s = pd.read_csv(Q16 / "q16c_sim.csv")
    ref = json.load(open(Q16 / "q16c_sim_noisefree.json"))
    s["ok"] = s["flag"] == 0
    s["bC"] = 100 * (s["C"] / ref["C"] - 1)
    s["bnom"] = 100 * (s["C_nom"] / ref["C_nom"] - 1)
    s["bratio"] = 100 * (s["C_ratio"] / ref["C_ratio"] - 1)
    s["bfixed"] = 100 * (s["C_fixed"] / ref["C_fixed"] - 1)
    return s, ref


def stage_net():
    """Network-wide (all CHM15k, q15f nightly table) predicted SNR heights, in-window SNR, u/C, noise bias."""
    import pandas as pd
    df, P = load_l1()
    shape, unit_lvl, net_lvl, _ = transfer_function(df, P)
    rng = np.r_[P["rng"], np.arange(P["rng"][-1] + 29.97, 15000.0, 29.97)]
    cr = nm(df.loc[df["stream"] != PAYERNE, "c_ratio"])
    s, ref = sim_tables()
    bias_fit = sim_bias_fit(s, df)
    nt = pd.read_csv(Q15 / "q15f_nightly.csv")
    meta = station_meta()
    out = []
    for key, g in nt.groupby("stream"):
        g = g.sort_values("t")
        alt = meta.get(key, dict(alt=0.0))["alt"]
        _, p_mol = molecular(alt, rng)
        lvl = unit_lvl.get(key, net_lvl)
        C1 = np.exp(g["lnC"].values) * cr
        sd = np.exp(g["ln_sd"].values)
        n = np.exp(g["ln_nclear"].values)
        sigma = lvl * shape[None, :] * sd[:, None] / np.sqrt(3 * n[:, None])
        snr = C1[:, None] * p_mol[None, :] / sigma
        zb, zt = g["zb_km"].values * 1000, g["zt_km"].values * 1000
        rec = dict(stream=np.repeat(key, len(g)), t=g["t"].values, lnC=g["lnC"].values,
                   rel_u=g["rel_u"].values, zb=zb, zt=zt, unit_lvl=np.repeat(key in unit_lvl, len(g)))
        above = rng[None, :] >= 1000
        for lev in SNR_LEVELS:
            below = (snr < lev) & above
            first = np.argmax(below, axis=1)
            rec[f"z_snr{int(lev)}"] = np.where(below.any(axis=1), rng[first], np.nan)
        win = (rng[None, :] >= zb[:, None]) & (rng[None, :] <= zt[:, None])
        with np.errstate(invalid="ignore", divide="ignore"):
            rec["snr_zt"] = np.array([np.interp(z, rng, row) for z, row in zip(zt, snr)])
            rec["snr_78"] = np.nanmedian(np.where((rng >= 7000) & (rng < 8000), snr, np.nan), axis=1)
            rec["frac_win_snr3"] = np.sum(win & (snr < 3), axis=1) / np.maximum(win.sum(axis=1), 1)
            rec["u_pred"] = np.sqrt(np.sum(np.where(win, 1 / snr ** 2, 0), axis=1)) / np.maximum(win.sum(axis=1), 1)
        rec["bias_pred"] = bias_fit(100 * rec["u_pred"])
        out.append(pd.DataFrame(rec))
    net = pd.concat(out, ignore_index=True)
    net["season"] = [season_of(t) for t in net["t"]]
    net["month"] = [(datetime(1970, 1, 1) + timedelta(days=float(t))).month for t in net["t"]]
    net.to_csv(Q16 / "q16c_net_nights.csv", index=False)
    # per-unit annual harmonics of the predicted noise quantities
    rows = []
    for key, g in net.groupby("stream"):
        r = dict(stream=key, n=len(g), payerne=key == PAYERNE)
        for col in ("z_snr3", "snr_zt", "u_pred", "bias_pred", "zt", "frac_win_snr3"):
            y = g[col].values
            if col == "u_pred":
                y = 100 * y
            amp, pk = annual_fit(g["t"].values, y)
            r[f"{col}_amp"], r[f"{col}_peak"], r[f"{col}_med"] = amp, pk, nm(y)
        a_c, p_c = annual_fit(g["t"].values, 100 * g["lnC"].values)
        r["lnC_amp"], r["lnC_peak"] = a_c, p_c
        # in-phase projection of the predicted bias onto the unit's own C_L cycle
        r["bias_inphase"] = (r["bias_pred_amp"] * np.cos(W_YEAR * (r["bias_pred_peak"] - p_c))
                             if np.isfinite(p_c) else np.nan)
        rows.append(r)
    pd.DataFrame(rows).to_csv(Q16 / "q16c_net_units.csv", index=False)
    print(f"net stage: {len(net)} nights, {len(rows)} units; transfer level net {net_lvl:.3g}, "
          f"c_ratio {cr:.3f}", flush=True)


def sim_bias_fit(s, df):
    """Noise-induced bias of the operational C (%) vs the noise-only u/C of its window (%), from the
    controlled scan with boundary-layer aerosol (the realistic, larger case): bias(mult) minus the
    bias at the lowest noise (aerosol geometry, not noise), linearly interpolated in u."""
    sc = s[(s["idx"] == -1) & s["ok"]]
    ok = (df["stream"] != PAYERNE).values
    u0 = 100 * nm(df.loc[ok, "u_noise"])
    g = sc.groupby("mult")["bC"].mean()
    u = np.r_[0.0, u0 * g.index.values]
    b = np.r_[0.0, g.values - g.values[0]]
    return lambda x: np.interp(np.asarray(x, float), u, b)


def stage_scatter():
    """Night-to-night scatter of ln C vs the predicted noise (task 4)."""
    import pandas as pd
    from q1_temperature import robust_mask, running_median
    net = pd.read_csv(Q16 / "q16c_net_nights.csv")
    s, _ = sim_tables()
    df, _ = load_l1()
    # noise -> scatter of the OPERATIONAL constant: std over repeats / window-mean u (per sim night)
    pn = s[(s["idx"] >= 0) & s["ok"]].groupby("idx")["bC"].std()
    u_n = 100 * df.loc[pn.index, "u_noise"].values
    k_sim = nm(pn.values / u_n)
    rows, pooled = [], []
    for key, g in net.groupby("stream"):
        g = g.sort_values("t")
        y = 100 * g["lnC"].values
        a = y - running_median(g["t"].values, y, 31.0)
        m = robust_mask(a) & np.isfinite(g["u_pred"].values)
        if m.sum() < 40:
            continue
        a, u = a[m], 100 * g["u_pred"].values[m] * k_sim
        sea = g["season"].values[m]
        rel_u = 100 * g["rel_u"].values[m]
        r = dict(stream=key, n=int(m.sum()), sd_anom=1.4826 * nm(np.abs(a - nm(a))), u_noise=nm(u),
                 rel_u_op=nm(rel_u))
        for sn in SEASONS:
            ss = sea == sn
            if ss.sum() >= 10:
                r[f"sd_{sn}"] = 1.4826 * nm(np.abs(a[ss] - nm(a[ss])))
                r[f"u_{sn}"] = nm(u[ss])
        # within-unit: |anomaly| vs predicted noise (anomalies normalised by the unit's scatter)
        pooled.append(pd.DataFrame(dict(stream=key, a=a - nm(a), u=u, urel=u / nm(u), sea=sea,
                                        rel_u=rel_u)))
        rows.append(r)
    U = pd.DataFrame(rows)
    U.to_csv(Q16 / "q16c_scatter_units.csv", index=False)
    Pd = pd.concat(pooled, ignore_index=True)
    Pd.to_csv(Q16 / "q16c_scatter_nights.csv", index=False)
    print(f"scatter stage: k_sim {k_sim:.2f}, {len(U)} units", flush=True)
    json.dump(dict(k_sim=k_sim), open(Q16 / "q16c_scatter_meta.json", "w"))

def harmonic_pooled(doy, y):
    """Pooled annual harmonic (offset + sin + cos on DOY): amplitude, peak DOY, amplitude s.e."""
    doy = np.asarray(doy, float)
    y = np.asarray(y, float)
    ok = np.isfinite(doy) & np.isfinite(y)
    ph = W_YEAR * (doy[ok] - 1.0)
    X = np.column_stack([np.ones(ok.sum()), np.sin(ph), np.cos(ph)])
    c, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    res = y[ok] - X @ c
    cov = res @ res / max(ok.sum() - 3, 1) * np.linalg.inv(X.T @ X)
    amp = float(np.hypot(c[1], c[2]))
    se = float(np.sqrt(0.5 * (cov[1, 1] + cov[2, 2])))
    return amp, float((np.arctan2(c[1], c[2]) / W_YEAR) % 365.25 + 1.0), se


def stage_summary():
    """Tables (csv), figures (png) and the pickled numbers used to write q16c_summary.md."""
    import pickle
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    df, P = load_l1()
    rng = P["rng"]
    pay = (df["stream"] == PAYERNE).values
    d = df[~pay]
    s, ref = sim_tables()
    s0 = pd.read_csv(Q16 / "q16c_sim0.csv")
    s0["ok"] = s0["flag"] == 0
    for c in ("C", "C_nom", "C_ratio", "C_fixed"):
        s0["b" + c] = 100 * (s0[c] - 1)
    net = pd.read_csv(Q16 / "q16c_net_nights.csv")
    nu = pd.read_csv(Q16 / "q16c_net_units.csv")
    cf = pd.read_csv(Q16 / "q16c_cf_nights.csv")
    cfX = np.load(Q16 / "q16c_cf_profiles.npz")["X"]
    U = pd.read_csv(Q16 / "q16c_scatter_units.csv")
    Nn = pd.read_csv(Q16 / "q16c_scatter_nights.csv")
    k_sim = json.load(open(Q16 / "q16c_scatter_meta.json"))["k_sim"]
    u0 = 100 * nm(d["u_noise"])
    seas = list(SEASONS)
    colors = {"DJF": "#1f77b4", "MAM": "#2ca02c", "JJA": "#d62728", "SON": "#ff7f0e"}
    mm = np.arange(1, 13)

    # ---------- task 1: SNR heights (L1 sample) ----------
    cols1 = ["z_snr1", "z_snr3", "z_snr5", "z250_snr1", "z250_snr3", "z250_snr5", "zb", "zt",
             "snr_zb", "snr_zt", "snr_7_8", "u_noise", "nb", "dur_h", "c_ratio"]
    t1 = d.groupby("season")[cols1].median().reindex(seas)
    t1.loc["all"] = d[cols1].median()
    t1.loc["Payerne A"] = df.loc[pay, cols1].median()
    t1["frac_nights_win_lt3"] = [*(d.groupby("season")["frac_win_snr3"].apply(lambda x: (x > 0).mean())
                                   .reindex(seas)), (d["frac_win_snr3"] > 0).mean(),
                                 (df.loc[pay, "frac_win_snr3"] > 0).mean()]
    t1["mean_frac_win_gates_lt3"] = [*(d.groupby("season")["frac_win_snr3"].mean().reindex(seas)),
                                     d["frac_win_snr3"].mean(), df.loc[pay, "frac_win_snr3"].mean()]
    t1.to_csv(Q16 / "q16c_table_snr_heights.csv")
    # network (cache HK, all nights) - sampled units with their own transfer, and all units
    sampled = net[net["unit_lvl"] & (net["stream"] != PAYERNE)]
    allnet = net[net["stream"] != PAYERNE]
    coln = ["z_snr1", "z_snr3", "z_snr5", "zb", "zt", "snr_zt", "snr_78", "u_pred"]
    tn = allnet.groupby("season")[coln].median().reindex(seas)
    tn["frac_nights_win_lt3"] = allnet.groupby("season")["frac_win_snr3"].apply(lambda x: (x > 0).mean()).reindex(seas)
    tn["mean_frac_win_gates_lt3"] = allnet.groupby("season")["frac_win_snr3"].mean().reindex(seas)
    tn.loc["all"] = [*allnet[coln].median(), (allnet["frac_win_snr3"] > 0).mean(), allnet["frac_win_snr3"].mean()]
    tn.to_csv(Q16 / "q16c_table_network.csv")
    nuu = nu[~nu["payerne"]]

    # ---------- task 2: profile shape on strictly cloud-free nights ----------
    xm = X_MID
    okcf = np.isfinite(cfX).any(axis=1) & (cf["stream"] != PAYERNE).values
    for lab, (lo, hi) in (("X78", (7000, 8000)), ("X910", (9000, 10000)), ("X1112", (11000, 12000)),
                          ("X1315", (13000, 15000)), ("X45", (4000, 5000))):
        b = (xm >= lo) & (xm < hi)
        cf[lab] = np.nanmean(cfX[:, b], axis=1)
    cf["doy"] = [doy_of(t) for t in cf["t"]]
    c2 = cf[okcf]
    colx = ["X45", "X78", "X910", "X1112", "X1315", "hi_frac"]
    t2 = c2.groupby("season")[colx].median().reindex(seas)
    t2.loc["all"] = c2[colx].median()
    t2["n"] = [*c2.groupby("season").size().reindex(seas), len(c2)]
    t2.to_csv(Q16 / "q16c_table_profile_shape.csv")
    a78 = harmonic_pooled(c2["doy"].values, 100 * np.log(np.clip(c2["X78"].values, 0.2, None)))
    t2b = d.groupby("season")[["slopeA_win", "offB_win", "chi2_win", "slopeA_78", "offB_78",
                               "offB_long_win_bias"]].median().reindex(seas)
    t2b.to_csv(Q16 / "q16c_table_slope_tests.csv")

    # ---------- task 3: simulation ----------
    def scan(tab, cols):
        okk = tab[tab["ok"]]
        g = okk.groupby("mult")[cols].mean()
        e = okk.groupby("mult")[cols].std() / np.sqrt(okk.groupby("mult").size().values[:, None])
        g["acc"] = tab.groupby("mult")["ok"].mean()
        g["zt"] = okk.groupby("mult")["zt"].mean()
        g["zb"] = okk.groupby("mult")["zb"].mean()
        g["sdC"] = okk.groupby("mult")[cols[0]].std()
        g["u"] = u0 * g.index
        return g, e
    g_a, e_a = scan(s[s["idx"] == -1], ["bC", "bnom", "bratio", "bfixed"])
    g_0, e_0 = scan(s0[s0["idx"] == -1], ["bC", "bC_nom", "bC_ratio", "bC_fixed"])
    g_a.join(e_a, rsuffix="_se").to_csv(Q16 / "q16c_sim_scan_aerosol.csv")
    g_0.join(e_0, rsuffix="_se").to_csv(Q16 / "q16c_sim_scan_purenoise.csv")
    pn = s[s["idx"] >= 0]
    a = pn[pn["ok"]].groupby("idx")[["bC", "bnom", "bratio", "bfixed"]].mean()
    a["acc"] = pn.groupby("idx")["ok"].mean()
    a = a.join(df[["stream", "season", "doy", "u_noise", "nb"]])
    a["month"] = [(datetime(2025, 1, 1) + timedelta(days=int(x) - 1)).month for x in a["doy"]]
    a.to_csv(Q16 / "q16c_sim_nights.csv")
    sim_sea = a.groupby("season")[["bC", "bnom", "bratio", "acc", "u_noise", "nb"]].mean().reindex(seas)
    sim_sea_se = a.groupby("season")[["bC", "bnom"]].sem().reindex(seas)
    hb = harmonic_pooled(a["doy"].values, a["bC"].values)
    hn = harmonic_pooled(a["doy"].values, a["bnom"].values)
    sea0 = s0[s0["idx"] <= -10].copy()
    sea0["season"] = [seas[-10 - int(i)] for i in sea0["idx"]]
    sea0g = sea0[sea0["ok"]].groupby("season")["bC"].agg(["mean", "sem"]).reindex(seas)

    # ---------- task 4: scatter ----------
    Uu = U[U["stream"] != PAYERNE]
    r_units = float(np.corrcoef(Uu["sd_anom"], Uu["u_noise"])[0, 1])
    Xd = np.column_stack([np.ones(len(Uu)), Uu["u_noise"] ** 2])
    cu, *_ = np.linalg.lstsq(Xd, Uu["sd_anom"] ** 2, rcond=None)
    share = float(np.median(Uu["u_noise"] ** 2 / Uu["sd_anom"] ** 2))
    Nd = Nn[Nn["stream"] != PAYERNE].copy()
    Nd["urel_s"] = Nd["u"] / Nd.groupby(["stream", "sea"])["u"].transform("median")
    Nd["a_s"] = Nd["a"] - Nd.groupby(["stream", "sea"])["a"].transform("median")
    Nd["q"] = pd.qcut(Nd["urel_s"], 8, duplicates="drop")
    oc = Nd.groupby("q", observed=True).agg(sd=("a_s", lambda x: 1.4826 * np.median(np.abs(x))),
                                            u=("u", "median"), urel=("urel_s", "median"), n=("a_s", "size"))
    oc.to_csv(Q16 / "q16c_scatter_octiles.csv")
    Xo = np.column_stack([np.ones(len(oc)), oc["u"] ** 2])
    co, *_ = np.linalg.lstsq(Xo, oc["sd"] ** 2, rcond=None)
    sea_sc = pd.DataFrame({sn: [Uu[f"sd_{sn}"].median(), Uu[f"u_{sn}"].median()] for sn in seas},
                          index=["sd_anom", "u_noise"]).T

    # ---------- figure 1: detectability ----------
    fig, ax = plt.subplots(1, 3, figsize=(17, 6.2))
    zk = rng / 1000
    for sn in seas:
        m = (df["season"] == sn).values & ~pay
        prof = np.nanmedian(P["snr"][m], axis=0)
        ax[0].semilogx(prof, zk, color=colors[sn], lw=1.8, label=f"{sn} (30 m gate)")
        ax[0].semilogx(prof * np.sqrt(250 / 29.97), zk, color=colors[sn], lw=1, ls="--")
    ax[0].semilogx(np.nanmedian(P["snr"][pay], axis=0), zk, color="k", lw=1.2, ls=":", label="Payerne A (30 m)")
    for lev in SNR_LEVELS:
        ax[0].axvline(lev, color="grey", lw=0.8, ls=":")
    ax[0].axhspan(nm(d["zb"]) / 1000, nm(d["zt"]) / 1000, color="0.85", zorder=0, label="median op. window")
    ax[0].axhspan(7, 8, color="#c9b3e6", alpha=0.5, zorder=0, label="7-8 km band (q15e)")
    ax[0].set_xlim(0.1, 300)
    ax[0].set_ylim(0, 10.2)
    ax[0].set_xlabel("molecular SNR of the night mean (C p_mol / sigma)")
    ax[0].set_ylabel("range above instrument [km]")
    ax[0].set_title("(a) night-mean molecular SNR (solid 30 m, dashed 250 m)")
    ax[0].legend(fontsize=7, loc="upper right")
    mon = sampled.groupby("month")
    for col, lab, cc in (("z_snr1", "SNR = 1", "0.55"), ("z_snr3", "SNR = 3", "k"), ("z_snr5", "SNR = 5", "0.3")):
        ax[1].plot(mm, mon[col].median().reindex(mm) / 1000, color=cc, marker="o", label=f"height of {lab} (30 m gate)")
    q1 = mon["z_snr3"].quantile(0.25).reindex(mm) / 1000
    q3 = mon["z_snr3"].quantile(0.75).reindex(mm) / 1000
    ax[1].fill_between(mm, q1, q3, color="k", alpha=0.12)
    ax[1].plot(mm, mon["zt"].median().reindex(mm) / 1000, color="tab:blue", marker="s", label="window top zt")
    ax[1].plot(mm, mon["zb"].median().reindex(mm) / 1000, color="tab:blue", marker="v", ls="--",
               label="window bottom zb")
    ax[1].axhspan(7, 8, color="#c9b3e6", alpha=0.4, zorder=0)
    ax[1].set_xticks(mm)
    ax[1].set_xlabel("month")
    ax[1].set_ylabel("range above instrument [km]")
    ax[1].set_ylim(0, 11)
    ax[1].set_title("(b) SNR heights vs window (26 units, all nights, cache HK)")
    ax[1].legend(fontsize=7, loc="lower left")
    for sn in seas:
        m = okcf & (cf["season"] == sn).values
        ax[2].plot(np.nanmedian(cfX[m], axis=0), xm / 1000, color=colors[sn], lw=1.6, label=f"{sn} (n={m.sum()})")
    q1, q3 = np.nanpercentile(cfX[okcf], [25, 75], axis=0)
    ax[2].fill_betweenx(xm / 1000, q1, q3, color="0.6", alpha=0.25, label="IQR, all nights")
    ax[2].axvline(1, color="k", lw=0.8)
    ax[2].axhspan(7, 8, color="#c9b3e6", alpha=0.4, zorder=0)
    ax[2].set_xlim(-0.5, 2.5)
    ax[2].set_ylim(1, 15)
    ax[2].set_xlabel("X = S / (C_window p_mol)   (1 = molecular slope)")
    ax[2].set_ylabel("range above instrument [km]")
    ax[2].set_title("(c) profile shape, strictly cloud-free nights")
    ax[2].legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(Q16 / "q16c_detectability.png", dpi=130)
    plt.close(fig)

    # ---------- figure 2: simulation ----------
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.6))
    names = ("operational C (ensemble median)", "nominal member", "ratio in chosen window", "fixed 3-5.5 km window")
    for (g, e, cols, lab, ls) in ((g_0, e_0, ["bC", "bC_nom", "bC_ratio", "bC_fixed"], "pure molecular", "-"),
                                  (g_a, e_a, ["bC", "bnom", "bratio", "bfixed"], "with BL aerosol", "--")):
        for c, nmz, cc in zip(cols, names, ("k", "tab:blue", "tab:green", "tab:purple")):
            ax[0].errorbar(g["u"], g[c], yerr=e[c], color=cc, ls=ls, marker="o", ms=3, capsize=2,
                           label=f"{nmz}, {lab}")
    for sn in ("DJF", "JJA"):
        ax[0].axvline(100 * nm(d.loc[d["season"] == sn, "u_noise"]), color=colors[sn], lw=1, ls=":")
    ax[0].axvline(100 * nm(df.loc[pay, "u_noise"]), color="0.4", lw=1, ls="-.")
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("noise-only u/C of the window mean [%] (dotted DJF/JJA, dash-dot Payerne)")
    ax[0].set_ylabel("bias of C [%] (vs noise-free)")
    ax[0].set_title("(a) noise bias of the operational estimator")
    ax[0].legend(fontsize=6, ncol=2)
    ax[1].plot(g_0["u"], 100 * g_0["acc"], "k-o", ms=3, label="accepted [%], pure molecular")
    ax[1].plot(g_a["u"], 100 * g_a["acc"], "k--o", ms=3, label="accepted [%], with aerosol")
    ax[1].set_xscale("log")
    ax[1].set_ylabel("accepted nights [%]")
    ax[1].set_xlabel("noise-only u/C [%]")
    axb = ax[1].twinx()
    axb.plot(g_0["u"], g_0["zt"] / 1000, "b-s", ms=3, label="window top (pure)")
    axb.plot(g_0["u"], g_0["zb"] / 1000, "b-v", ms=3, label="window bottom (pure)")
    axb.set_ylabel("chosen window [km]", color="b")
    ax[1].set_title("(b) acceptance and chosen window vs noise")
    ax[1].legend(fontsize=7, loc="lower left")
    axb.legend(fontsize=7, loc="center left")
    mg = a.groupby("month")[["bC", "bnom"]].agg(["mean", "sem"])
    ax[2].errorbar(mm, mg[("bC", "mean")].reindex(mm), yerr=mg[("bC", "sem")].reindex(mm), color="k", marker="o",
                   capsize=2, label="operational C")
    ax[2].errorbar(mm, mg[("bnom", "mean")].reindex(mm), yerr=mg[("bnom", "sem")].reindex(mm), color="tab:blue",
                   marker="s", capsize=2, label="nominal member")
    ax[2].axhline(0, color="k", lw=0.6)
    ax[2].set_xticks(mm)
    ax[2].set_xlabel("month")
    ax[2].set_ylabel("simulated bias [%] (each night's own noise)")
    ax[2].set_ylim(-3, 3)
    ax[2].set_title("(c) noise bias by month (observed C_L cycle: +-16 %)")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(Q16 / "q16c_sim_bias.png", dpi=130)
    plt.close(fig)

    # ---------- figure 3: scatter ----------
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.4))
    ax[0].scatter(Uu["u_noise"], Uu["sd_anom"], s=12, color="k", label="CHM15k units")
    pr = U[U["stream"] == PAYERNE]
    ax[0].scatter(pr["u_noise"], pr["sd_anom"], s=30, color="r", label="Payerne A")
    uu = np.linspace(0, max(Uu["u_noise"].max(), float(pr["u_noise"].max())) * 1.05, 50)
    base = np.sqrt(max(cu[0], 0))
    ax[0].plot(uu, np.sqrt(base ** 2 + uu ** 2), "b-", label="noise added in quadrature")
    ax[0].plot(uu, uu, "k:", lw=0.8, label="1:1 (noise only)")
    ax[0].set_xlabel("predicted noise scatter of C [%] (unit median night)")
    ax[0].set_ylabel("observed robust scatter of ln C anomalies [%]")
    ax[0].set_title(f"(a) across units (r = {r_units:.2f})")
    ax[0].legend(fontsize=7)
    ax[1].plot(oc["u"], oc["sd"], "ko-", label="observed (unit x season controlled)")
    b0 = np.sqrt(max(oc["sd"].iloc[0] ** 2 - oc["u"].iloc[0] ** 2, 0))
    ax[1].plot(oc["u"], np.sqrt(b0 ** 2 + oc["u"] ** 2), "b-", label="expected from noise")
    ax[1].set_xlabel("predicted noise scatter of C [%] (octiles of relative noise)")
    ax[1].set_ylabel("robust scatter of ln C anomalies [%]")
    ax[1].set_title("(b) within unit and season: noisier nights")
    ax[1].legend(fontsize=7)
    x = np.arange(4)
    ax[2].bar(x - 0.2, sea_sc["sd_anom"], 0.4, color="0.5", label="observed scatter")
    ax[2].bar(x + 0.2, sea_sc["u_noise"], 0.4, color="tab:blue", label="predicted noise")
    ax[2].set_xticks(x)
    ax[2].set_xticklabels(seas)
    ax[2].set_ylabel("[%] (network median of units)")
    ax[2].set_title("(c) night-to-night scatter vs noise by season")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(Q16 / "q16c_scatter.png", dpi=130)
    plt.close(fig)

    res = dict(t1=t1, tn=tn, nuu=nuu, t2=t2, a78=a78, t2b=t2b, g_a=g_a, e_a=e_a, g_0=g_0, e_0=e_0,
               sim_sea=sim_sea, sim_sea_se=sim_sea_se, hb=hb, hn=hn, sea0g=sea0g, r_units=r_units, cu=cu,
               share=share, oc=oc, co=co, sea_sc=sea_sc, k_sim=k_sim, u0=u0, Uu=Uu, U=U)
    pickle.dump(res, open(Q16 / "q16c_summary_inputs.pkl", "wb"))
    print("summary: figures and tables written", flush=True)



if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("l1", "all"):
        stage_l1()
    if stage in ("cf", "all"):
        stage_cf()
    if stage in ("sim", "all"):
        stage_sim()
    if stage in ("sim0", "all"):
        stage_sim0()
    if stage in ("net", "all"):
        stage_net()
    if stage in ("scatter", "all"):
        stage_scatter()
    if stage in ("summary", "all"):
        stage_summary()
