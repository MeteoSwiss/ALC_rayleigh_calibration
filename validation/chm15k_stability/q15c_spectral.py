# -*- coding: utf-8 -*-
"""q15c - spectral / wavelength factors of the CHM15k calibration-constant variability (1064 nm).

Question (q15 brief): the CHM15k Rayleigh constant C_L has an annual cycle of +-16 % (peak ~DOY 176) but
kappa x T_int(night) removes only ~35 % of it. Which spectral factors could explain more?

Factors (one part each):
  1. filter passband drifting against the laser line, driven by the internal temperature: linear kappa x
     T_int vs (a) other driver temperatures (temp_ext, temperature_optical_module, temperature_detector),
     (b) exponential moving averages of T_int / temp_ext with lags of 3 h ... 8 d (a filter that integrates
     the daytime heat), (c) a T_int + temp_ext mix, (d) a Gaussian-passband (non-linear) response. For each
     driver the coefficient is estimated on 31-day anomalies only (never on the raw seasonal series), and the
     annual cycle it predicts is subtracted from C_L before refitting the paper's Kalman annual sinusoid.
  2. laser-line drift when the laser-optical module (LOM) leaves regulation (temp_lom excursions).
  3. water-vapour absorption at 1064 nm (HITRAN 2020 lines via HAPI + a bounded continuum) on the two-way
     path to the Rayleigh calibration window, with CAMS 0.4 deg L137 humidity per calibration night.
  4. molecular cross-section / refractive-index / rotational-Raman effects of a wavelength shift (analytic).

Run (each part caches its CSV; `summary` writes the markdown + figures):
  python validation/chm15k_stability/q15c_spectral.py hk        # factors 1 + 2 (10-min HK caches)
  python validation/chm15k_stability/q15c_spectral.py wv        # factor 3 (needs hapi on PYTHONPATH)
  python validation/chm15k_stability/q15c_spectral.py summary   # factor 4 + table + figures
Outputs: C:/DATA/Projects/202610_CHM15k_stability/q15/q15c_*
"""
from __future__ import annotations

import csv
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
from common import BLOCKS_PER_DAY, CENSUS, OUT, chm15k_streams, load_cache, load_constants, solar_hour, stream_key  # noqa: E402
from q1_temperature import ANOM_WINDOW, running_median  # noqa: E402
from q13_kalman_network import annual  # noqa: E402

warnings.filterwarnings("ignore")
Q15 = OUT / "q15"
HITRAN_DIR = Q15 / "q15c_hitran"
EPOCH = datetime(1970, 1, 1)
LN2 = np.log(2.0)
W_DOY = 2 * np.pi / 365.25
PAYERNE = "0-20000-0-06610_A"
N_WORKERS = 4

# Spectroscopy (operator working assumption + literature, see brief)
LAMBDA_VAC = 1064.49          # nm, Nd:YAG line in vacuum (CIPM I2-stabilised standard 1064.4906 nm)
N_AIR = 1.000274              # refractive index of air at 1064 nm (Ciddor, 15 degC, 101.325 kPa)
FILTER_FWHM = 1.0             # nm, Lufft datasheet
ALPHA_W_FILTER = 0.020        # 1/K, q11 inversion network median (alpha/W)
DELTA0_W = 0.19               # laser red of the filter centre, in units of W (q11 median)
LASER_TUNING = {"cavity mode (microchip)": 0.013, "gain centre": 0.0042}   # nm/K (gain 3...5.5 pm/K)

TAUS_H = [0.0, 3.0, 12.0, 24.0, 48.0, 96.0, 192.0]


# ---------------------------------------------------------------------------------------------------
# generic helpers
# ---------------------------------------------------------------------------------------------------
def harmonic(t_days, y):
    """Annual harmonic of y(t) by least squares: (amplitude, DOY of maximum). y in any unit."""
    ok = np.isfinite(y) & np.isfinite(t_days)
    if ok.sum() < 30:
        return np.nan, np.nan
    doy = np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in t_days[ok]], float)
    X = np.column_stack([np.ones(ok.sum()), np.sin(W_DOY * doy), np.cos(W_DOY * doy)])
    c, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    amp = float(np.hypot(c[1], c[2]))
    peak = (np.arctan2(c[1], c[2]) / W_DOY) % 365.25
    return amp, float(peak)


def anomaly(t, x):
    return x - running_median(t, x, ANOM_WINDOW)


def ols_slope(x, y):
    """Slope, standard error and r^2 of y on x (4-MAD screen on y)."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 40 or np.std(x[ok]) == 0:
        return np.nan, np.nan, np.nan
    med = np.median(y[ok])
    ok &= np.abs(y - med) <= 4 * 1.4826 * np.median(np.abs(y[ok] - med)) + 1e-12
    X = np.column_stack([np.ones(ok.sum()), x[ok]])
    c, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    res = y[ok] - X @ c
    s2 = res @ res / max(ok.sum() - 2, 1)
    se = np.sqrt(s2 * np.linalg.pinv(X.T @ X)[1, 1])
    r2 = np.corrcoef(x[ok], y[ok])[0, 1] ** 2
    return float(c[1]), float(se), float(r2)


def ema_regular(x, tau_h, dt_h=1.0 / 6.0):
    """First-order low-pass (thermal lag tau) of a regular 10-min series; gaps are held at the last value."""
    from scipy.signal import lfilter
    y = np.asarray(x, float).copy()
    ok = np.isfinite(y)
    if not ok.any():
        return y
    idx = np.where(ok, np.arange(len(y)), 0)
    np.maximum.accumulate(idx, out=idx)
    y = y[idx]
    y[: np.argmax(ok)] = y[np.argmax(ok)]
    if tau_h <= 0:
        return y
    a = 1.0 - np.exp(-dt_h / tau_h)
    out, _ = lfilter([a], [1.0, -(1.0 - a)], y, zi=[(1.0 - a) * y[0]])
    return out


# ---------------------------------------------------------------------------------------------------
# factors 1 + 2: housekeeping-driven spectral drift
# ---------------------------------------------------------------------------------------------------
def load_kappas():
    k_in = {r["stream"]: float(r["b_T_ctrl"]) for r in csv.DictReader(open(OUT / "q2a_intranight.csv", encoding="utf-8"))}
    k_m2 = {}
    for r in csv.DictReader(open(OUT / "q1_slopes.csv", encoding="utf-8")):
        if r["arm"] == "A" and r["tvar"] == "temp_int" and r["b2"] not in ("", "nan"):
            k_m2[(r["stream"], r["era"])] = float(r["b2"])
    q11 = {}
    for r in csv.DictReader(open(OUT / "q11_filter_inversion.csv", encoding="utf-8")):
        if r["alpha_over_W"] not in ("", "nan"):
            q11[(r["stream"], r["era"])] = (float(r["alpha_over_W"]), float(r["delta0_over_W"]), float(r["T_mean"]),
                                            r["c_significant"] == "True")
    tlom_night = {r["stream"]: float(r["b_Tlom"]) for r in csv.DictReader(open(OUT / "q7_tlom_night.csv", encoding="utf-8"))}
    return k_in, k_m2, q11, tlom_night


def regular_grid(cache):
    """Regular 10-min grid of the cache (missing days -> NaN). Returns t (n,), hk (n, nvar)."""
    day = np.round(cache["day"]).astype(int)
    d0, d1 = day.min(), day.max()
    nd = d1 - d0 + 1
    hk = np.full((nd, BLOCKS_PER_DAY, cache["hk"].shape[-1]), np.nan, np.float32)
    hk[day - d0] = cache["hk"]
    t = ((np.arange(d0, d1 + 1)[:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY)).ravel()
    return t, hk.reshape(-1, hk.shape[-1]).astype(float)


def window_median(t, x, t0, t1, min_frac=0.5):
    i0 = np.searchsorted(t, t0)
    i1 = np.searchsorted(t, t1, side="right")
    out = np.full(len(t0), np.nan)
    for k in range(len(t0)):
        seg = x[i0[k]:i1[k]]
        if seg.size and np.isfinite(seg).mean() >= min_frac:
            out[k] = np.nanmedian(seg)
    return out


def noon_minus_night(t, x, raw_ok, t0, t1, lon, months=None):
    """Median over days of [x at local solar 12-16 h] - [x in the calibration window of the same date's night].

    The night is the calibration window that ENDS on the morning of that date (the constant applied to the day)."""
    sh = solar_hour(t, lon)
    day_of = np.floor(t + lon / 360.0).astype(int)
    noon = (sh >= 12) & (sh < 16) & raw_ok & np.isfinite(x)
    d0 = day_of.min()
    idx = day_of[noon] - d0
    n = np.bincount(idx, minlength=day_of.max() - d0 + 1)
    s = np.bincount(idx, weights=x[noon], minlength=day_of.max() - d0 + 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        dmean = np.where(n >= 12, s / np.maximum(n, 1), np.nan)
    night = window_median(t, np.where(raw_ok, x, np.nan), t0, t1)
    diffs = []
    for e, v in zip(t1, night):
        k = int(np.floor(e + lon / 360.0)) - d0
        if not np.isfinite(v) or k < 0 or k >= len(dmean) or not np.isfinite(dmean[k]):
            continue
        if months is not None and (EPOCH + timedelta(days=float(k + d0))).month not in months:
            continue
        diffs.append(dmean[k] - v)
    return float(np.median(diffs)) if len(diffs) >= 10 else np.nan


def unit_hk(args):
    """Factors 1 and 2 for one CHM15k stream (main optical module)."""
    wmo, ident, site, lat, lon, k_in, k_m2_all, q11_all, b_tlom = args
    key = stream_key(wmo, ident)
    cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
    if cons is None or cache is None:
        return None
    mod = max(set(cons["module"]), key=list(cons["module"]).count)
    m = (cons["module"] == mod) & (cons["c"] > 0) & np.isfinite(cons["t0"]) & np.isfinite(cons["t1"])
    if m.sum() < 100:
        return None
    td, t0, t1, c = cons["t"][m], cons["t0"][m], cons["t1"][m], cons["c"][m]
    names = list(cache["hk_vars"])
    tg, hk = regular_grid(cache)
    raw = {v: hk[:, names.index(v)] for v in ("temp_int", "temp_ext", "temperature_optical_module",
                                              "temperature_detector")}
    # night drivers
    drivers = {}
    for tau in TAUS_H:
        drivers[f"Tint_ema{int(tau)}h"] = ema_regular(raw["temp_int"], tau)
    for tau in (0.0, 24.0, 96.0):
        drivers[f"Text_ema{int(tau)}h"] = ema_regular(raw["temp_ext"], tau)
    drivers["Tlom"] = raw["temperature_optical_module"]
    drivers["Tdet"] = raw["temperature_detector"]
    raw_ok = np.isfinite(raw["temp_int"])
    ext_ok = np.isfinite(raw["temp_ext"])
    night = {}
    for name, x in drivers.items():
        # lagged series are gap-held; only nights with the raw sensor present are used
        okm = raw_ok if name.startswith("Tint") else (ext_ok if name.startswith("Text") else np.isfinite(x))
        night[name] = window_median(tg, np.where(okm, x, np.nan), t0, t1)
    # previous-afternoon heat: max of T_int in the 18 h preceding the calibration window
    i0 = np.searchsorted(tg, t0 - 0.75)
    i1 = np.searchsorted(tg, t0)
    night["Tint_prevmax"] = np.array([np.nanmax(raw["temp_int"][a:b]) if b > a and np.isfinite(raw["temp_int"][a:b]).any()
                                      else np.nan for a, b in zip(i0, i1)])

    times = [EPOCH + timedelta(days=float(x)) for x in td]
    logc = np.log(c)
    la = anomaly(td, logc)
    a_op, p_op, s_op = annual(times, c)
    amp_lnC, peak_lnC = harmonic(td, 100 * (logc - np.nanmedian(logc)))
    row = dict(stream=key, site=site, module=mod, n=int(m.sum()), amp_op=a_op, peak_op=p_op, scat_op=s_op,
               amp_lnC=amp_lnC, peak_lnC=peak_lnC, kappa_in=k_in, kappa_m2=k_m2_all.get((key, mod), np.nan))
    sd_la = np.nanstd(la[np.abs(la - np.nanmedian(la)) < 4 * 1.4826 * np.nanmedian(np.abs(la - np.nanmedian(la)))])
    row["sd_anom_lnC"] = 100 * sd_la

    def evaluate(tag, pred_pct, kappa=np.nan, r2=np.nan, xa=None, dnoon=np.nan, dnoon_jja=np.nan):
        """pred_pct: predicted ln C effect (%) per calibration night."""
        ok = np.isfinite(pred_pct)
        if ok.sum() < 100:
            return
        cn = np.where(ok, c * np.exp(-pred_pct / 100.0), np.nan)
        tt = [x for x, o in zip(times, ok) if o]
        a_after, _, s_after = annual(tt, cn[ok])
        a_ref, _, _ = annual(tt, c[ok])
        ap, pp = harmonic(td, pred_pct)
        row[f"{tag}_kappa"] = kappa
        row[f"{tag}_r2"] = r2
        row[f"{tag}_amp_pred"] = ap
        row[f"{tag}_peak_pred"] = pp
        row[f"{tag}_amp_after"] = a_after
        row[f"{tag}_removed"] = 1 - a_after / a_ref if np.isfinite(a_ref) and a_ref > 0 else np.nan
        row[f"{tag}_scat_after"] = s_after
        pa = anomaly(td, pred_pct)
        row[f"{tag}_n2n_sd"] = float(np.nanstd(pa))
        row[f"{tag}_diurnal"] = kappa * dnoon if np.isfinite(kappa) else np.nan
        row[f"{tag}_diurnal_jja"] = kappa * dnoon_jja if np.isfinite(kappa) else np.nan

    # (0) reference: intra-night kappa x night T_int (should reproduce q13 amp_in) and multi-day kappa
    T = night["Tint_ema0h"]
    tref = np.nanmedian(T)
    dn_T = noon_minus_night(tg, raw["temp_int"], raw_ok, t0, t1, lon)
    dn_Tj = noon_minus_night(tg, raw["temp_int"], raw_ok, t0, t1, lon, months=(6, 7, 8))
    evaluate("lin_in", k_in * (T - tref), k_in, np.nan, dnoon=dn_T, dnoon_jja=dn_Tj)
    if np.isfinite(row["kappa_m2"]):
        evaluate("lin_m2", row["kappa_m2"] * (T - tref), row["kappa_m2"], dnoon=dn_T, dnoon_jja=dn_Tj)
    aT, pT = harmonic(td, T)
    row["Tint_amp"], row["Tint_peak"] = aT, pT
    row["kappa_needed_Tint"] = amp_lnC / aT if aT > 0 else np.nan

    # (1a,b) every driver with its own anomaly coefficient
    for name, x in night.items():
        xa = anomaly(td, x)
        k, se, r2 = ols_slope(xa, 100 * la)
        if not np.isfinite(k):
            continue
        ax_, px_ = harmonic(td, x)
        row[f"{name}_amp"], row[f"{name}_peak"] = ax_, px_
        row[f"{name}_kneeded"] = amp_lnC / ax_ if ax_ > 0 else np.nan
        if name in drivers:
            dn = noon_minus_night(tg, drivers[name], raw_ok, t0, t1, lon)
            dnj = noon_minus_night(tg, drivers[name], raw_ok, t0, t1, lon, months=(6, 7, 8))
        else:
            dn = dnj = np.nan
        evaluate(name, k * (x - np.nanmedian(x)), k, r2, dnoon=dn, dnoon_jja=dnj)
        row[f"{name}_se"] = se

    # (1c) T_int + temp_ext mix (joint anomaly regression)
    xi, xe = anomaly(td, night["Tint_ema0h"]), anomaly(td, night["Text_ema0h"])
    ok = np.isfinite(xi) & np.isfinite(xe) & np.isfinite(la)
    if ok.sum() > 60:
        y = 100 * la
        med = np.median(y[ok])
        ok &= np.abs(y - med) <= 4 * 1.4826 * np.median(np.abs(y[ok] - med))
        X = np.column_stack([np.ones(ok.sum()), xi[ok], xe[ok]])
        cf, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
        r2 = 1 - np.var(y[ok] - X @ cf) / np.var(y[ok])
        pred = cf[1] * (night["Tint_ema0h"] - np.nanmedian(night["Tint_ema0h"])) + \
            cf[2] * (night["Text_ema0h"] - np.nanmedian(night["Text_ema0h"]))
        row["mix_k_int"], row["mix_k_ext"] = cf[1], cf[2]
        evaluate("mix", pred, np.nan, r2)

    # (1c') fast + slow: intra-night kappa x T_int (fixed, measured within nights) plus a slow lagged
    # component fitted on the residual night-to-night anomalies (a filter with two thermal paths)
    for tau in (24, 96, 192):
        xs = night[f"Tint_ema{tau}h"]
        resid = la - k_in / 100 * anomaly(td, T)
        ks, _, r2s = ols_slope(anomaly(td, xs), 100 * resid)
        if np.isfinite(ks):
            evaluate(f"fastslow{tau}h", k_in * (T - np.nanmedian(T)) + ks * (xs - np.nanmedian(xs)), ks, r2s,
                     dnoon=dn_T * k_in / ks if ks else np.nan, dnoon_jja=dn_Tj * k_in / ks if ks else np.nan)

    # (1d) Gaussian passband: same local slope as the linear model, plus the passband curvature
    def gauss_pred(aw, d0w, Tm, x):
        # filter centre drifts red by alpha per K; Delta = laser - filter centre shrinks when Delta0 > 0
        return 100 * (-4 * LN2) * ((d0w - aw * (x - Tm)) ** 2 - d0w ** 2)
    q = q11_all.get((key, mod))
    if q is not None:
        aw, d0w, Tm, sig = q
        evaluate("gauss_q11", gauss_pred(aw, d0w, Tm, T))
        evaluate("lin_q11", 100 * 8 * LN2 * aw * d0w * (T - Tm))     # its tangent at T_mean
        row["gauss_q11_sig"] = sig
    km = row["kappa_m2"] if np.isfinite(row["kappa_m2"]) else k_in
    d0w = km / 100 / (8 * LN2 * ALPHA_W_FILTER)
    evaluate("gauss_gen", gauss_pred(ALPHA_W_FILTER, d0w, tref, T))
    # Gaussian on the best lagged driver (4 d EMA of T_int)
    xl = night["Tint_ema96h"]
    kl = row.get("Tint_ema96h_kappa", np.nan)
    if np.isfinite(kl):
        d0l = kl / 100 / (8 * LN2 * ALPHA_W_FILTER)
        evaluate("gauss_ema96", gauss_pred(ALPHA_W_FILTER, d0l, np.nanmedian(xl), xl))

    # ---- factor 2: LOM regulation ----
    lom = raw["temperature_optical_module"]
    nl = night["Tlom"]
    if np.isfinite(nl).sum() > 50:
        vals, cnt = np.unique(np.round(nl[np.isfinite(nl)] * 20) / 20, return_counts=True)
        sp = float(vals[np.argmax(cnt)])
        exc_n = nl - sp
        row["lom_setpoint"] = sp
        row["lom_night_frac_out"] = float(np.mean(np.abs(exc_n[np.isfinite(exc_n)]) > 0.5))
        row["lom_night_amp"], row["lom_night_peak"] = harmonic(td, exc_n)
        row["lom_night_anom_sd"] = float(np.nanstd(anomaly(td, exc_n)))
        sh = solar_hour(tg, lon)
        aft = (sh >= 12) & (sh < 18) & np.isfinite(lom)
        mon = np.array([(EPOCH + timedelta(days=float(d))).month for d in np.floor(tg[aft])])
        e_aft = lom[aft] - sp
        row["lom_day_frac_out"] = float(np.mean(e_aft > 0.5))
        jja = np.isin(mon, (6, 7, 8))
        row["lom_day_frac_out_jja"] = float(np.mean(e_aft[jja] > 0.5)) if jja.any() else np.nan
        row["lom_day_exc_mean"] = float(np.mean(np.clip(e_aft, 0, None)))
        row["lom_day_exc_jja"] = float(np.mean(np.clip(e_aft[jja], 0, None))) if jja.any() else np.nan
        row["lom_day_exc_p95_jja"] = float(np.percentile(e_aft[jja], 95)) if jja.any() else np.nan
        row["b_Tlom_q7"] = b_tlom
    return row


def run_hk():
    k_in, k_m2, q11, tlom = load_kappas()
    jobs = [(w, i, s, la, lo, k_in[stream_key(w, i)], k_m2, q11, tlom.get(stream_key(w, i), np.nan))
            for w, i, s, la, lo in chm15k_streams() if stream_key(w, i) in k_in]
    rows = []
    with ProcessPoolExecutor(N_WORKERS) as ex:
        for r in ex.map(unit_hk, jobs, chunksize=2):
            if r is not None:
                rows.append(r)
                print(f"  {r['stream']:28s} amp_op {r['amp_op']:5.1f}  lin_in removed {r.get('lin_in_removed', np.nan):+.2f}"
                      f"  ema96 removed {r.get('Tint_ema96h_removed', np.nan):+.2f}", flush=True)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("stream", "site", "module"), k))
    with open(Q15 / "q15c_hk_units.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"-> {Q15 / 'q15c_hk_units.csv'} ({len(rows)} units)")


# ---------------------------------------------------------------------------------------------------
# factor 3: water-vapour absorption at 1064 nm
# ---------------------------------------------------------------------------------------------------
P_GRID = np.array([0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.05])     # atm
T_GRID = np.array([220.0, 240.0, 255.0, 270.0, 285.0, 300.0, 315.0])   # K
LAM_SCAN = np.round(np.arange(1064.29, 1064.7101, 0.01), 3)            # vacuum nm (laser positions scanned)


def build_cross_sections():
    """sigma(lambda, p, T) [cm2/molecule] of H2O lines (HITRAN, Voigt, 25 cm-1 wing cut-off) on a LUT,
    for the laser positions LAM_SCAN and the 1-nm filter-weighted band. Cached as npz."""
    f = Q15 / "q15c_h2o_sigma_lut.npz"
    if f.exists():
        with np.load(f) as z:
            return {k: z[k] for k in z.files}
    import hapi
    hapi.db_begin(str(HITRAN_DIR))
    nu_scan = 1e7 / LAM_SCAN
    sig = np.zeros((len(LAM_SCAN), len(P_GRID), len(T_GRID)))
    sig_band = np.zeros((len(P_GRID), len(T_GRID)))
    spec_ref = None
    for ip, p in enumerate(P_GRID):
        for it, T in enumerate(T_GRID):
            nu, k = hapi.absorptionCoefficient_Voigt(
                SourceTables="H2O_all", WavenumberRange=[9360, 9430], WavenumberStep=0.002,
                Environment={"p": float(p), "T": float(T)}, Diluent={"air": 0.99, "self": 0.01},
                OmegaWing=25, HITRAN_units=True)
            nu, k = np.asarray(nu), np.asarray(k)
            sig[:, ip, it] = np.interp(nu_scan, nu, k)
            lam = 1e7 / nu
            # filter: Gaussian FWHM 1 nm centred DELTA0_W x W blue of the laser line
            lc = LAMBDA_VAC - DELTA0_W * FILTER_FWHM
            w = np.exp(-4 * LN2 * ((lam - lc) / FILTER_FWHM) ** 2)
            sig_band[ip, it] = np.sum(w * k) / np.sum(w)
            if abs(p - 0.9) < 1e-6 and T == 285.0:
                spec_ref = np.column_stack([lam, k])
    out = dict(sig=sig, sig_band=sig_band, lam_scan=LAM_SCAN, p=P_GRID, T=T_GRID, spec_ref=spec_ref)
    np.savez(f, **out)
    return out


def sigma_interp(lut2d, p_atm, T):
    from scipy.interpolate import RegularGridInterpolator
    f = RegularGridInterpolator((P_GRID, T_GRID), lut2d, bounds_error=False, fill_value=None)
    pts = np.column_stack([np.clip(p_atm, P_GRID[0], P_GRID[-1]), np.clip(T, T_GRID[0], T_GRID[-1])])
    return f(pts)


# continuum (MT_CKD-type, cm2 molec-1 atm-1); literature bounds for the 1.06 um window
CONT = {"central": (2e-26, 2e-27), "upper": (1e-25, 1e-26)}   # (self Cs, foreign Cf)

WV_SAMPLE = None   # set in run_wv


def wv_months():
    out = []
    for y in (2025, 2026):
        for mth in range(1, 13):
            f = Path(f"D:/CAMS_Monthly_04/CAMS_Beta_{y}{mth:02d}.nc")
            if f.exists():
                out.append((y, mth, f))
    return out


def unit_wv(args):
    """Per calibration night: PWV, column below the window, two-way H2O transmission bias at the laser
    line(s), for the filter-weighted band and for the continuum bounds."""
    wmo, ident, site, lat, lon, alt, lut = args
    from calibration.water_vapor_correction.water_vapor import _hydrostatic_z_p, _read_cams_columns, _wv_number_density
    cons = load_constants(wmo, ident)
    if cons is None:
        return None
    key = stream_key(wmo, ident)
    i_line = int(np.argmin(np.abs(lut["lam_scan"] - LAMBDA_VAC)))
    rows = []
    for y, mth, f in wv_months():
        d0 = (datetime(y, mth, 1) - EPOCH).days
        d1 = (datetime(y + (mth == 12), mth % 12 + 1, 1) - EPOCH).days
        sel = np.where((cons["t"] >= d0) & (cons["t"] < d1) & (cons["c"] > 0))[0]
        if sel.size == 0:
            continue
        try:
            tdt, level, T, q, zs, lnsp = _read_cams_columns(str(f), float(lat), float(lon))
        except Exception as e:     # noqa: BLE001
            print(f"  {key} {f.name}: {e}", flush=True)
            continue
        z, P = _hydrostatic_z_p(level, T, q, zs, lnsp)       # (lev, nt)
        nwv = _wv_number_density(q, P, T)                     # m-3
        e_atm = nwv * 1.380649e-23 * T / 101325.0
        tnum = (tdt - np.datetime64("1970-01-01")) / np.timedelta64(1, "D")
        for k in sel:
            it = np.where((tnum >= cons["t0"][k] - 0.0625) & (tnum <= cons["t1"][k] + 0.0625))[0]
            if it.size == 0:
                it = np.array([int(np.argmin(np.abs(tnum - cons["t"][k])))])
            res = {}
            for j in it:
                o = np.argsort(z[:, j])
                zz, pp, tt, nn, ee = z[o, j], P[o, j] / 101325.0, T[o, j], nwv[o, j], e_atm[o, j]
                zg = np.arange(alt, min(zz.max(), 12000.0), 20.0)
                n_g = np.interp(zg, zz, nn)          # below the lowest model level: constant fill
                p_g = np.interp(zg, zz, pp)
                t_g = np.interp(zg, zz, tt)
                e_g = np.interp(zg, zz, ee)
                dz_cm = 20.0 * 100.0
                n_cm = n_g * 1e-6                    # molec cm-3
                win = (zg >= cons["zb"][k]) & (zg <= cons["zt"][k])
                if win.sum() < 2:
                    continue

                def two_way(sig_profile):
                    tau = np.cumsum(sig_profile * n_cm) * dz_cm
                    return float(np.log(np.mean(np.exp(-2 * tau[win]))))      # ln <T^2> over the window

                acc = {"pwv": np.sum(n_g) * 20.0 / 6.022e23 * 0.018015,      # kg m-2 = mm
                       "col_below_zb": np.sum(n_g[zg < cons["zb"][k]]) * 20.0 / 6.022e23 * 0.018015}
                for il, lam in enumerate(lut["lam_scan"]):
                    if il % 2 and il != i_line:
                        continue
                    acc[f"b_{lam:.2f}"] = two_way(sigma_interp(lut["sig"][il], p_g, t_g))
                acc["b_line"] = acc[f"b_{lut['lam_scan'][i_line]:.2f}"]
                acc["b_band"] = two_way(sigma_interp(lut["sig_band"], p_g, t_g))
                for tag, (cs, cf) in CONT.items():
                    acc[f"b_cont_{tag}"] = two_way(cs * e_g + cf * (p_g - e_g))
                for kk, v in acc.items():
                    res.setdefault(kk, []).append(v)
            if not res:
                continue
            r = dict(stream=key, site=site, t=cons["t"][k], c=cons["c"][k], module=cons["module"][k])
            r.update({kk: float(np.mean(v)) for kk, v in res.items()})
            rows.append(r)
    return rows


def geographic_sample(n=30):
    """Spread the WV sample over the network (one unit per ~grid cell of 4 deg), Payerne always included."""
    al = {(s["wmo"], s["ident"]): float(s["alt"]) for s in json.load(open(CENSUS, encoding="utf-8"))}
    streams = [(w, i, s, la, lo, al[(w, i)]) for w, i, s, la, lo in chm15k_streams()
               if (OUT / "hk" / f"{w}_{i}.npz").exists() and 27 < la < 80 and -30 < lo < 45]
    cells = {}
    for st in streams:
        cells.setdefault((int(st[3] // 4), int(st[4] // 4)), []).append(st)
    pick = [v[0] for v in cells.values()]
    if not any(p[0] == "0-20000-0-06610" and p[1] == "A" for p in pick):
        pick += [st for st in streams if st[0] == "0-20000-0-06610" and st[1] == "A"]
    return pick[:n] if len(pick) > n else pick


def run_wv():
    lut = build_cross_sections()
    sample = geographic_sample(40)
    print(f"WV sample: {len(sample)} units")
    rows = []
    with ProcessPoolExecutor(N_WORKERS) as ex:
        for rr in ex.map(unit_wv, [(*s, lut) for s in sample]):
            if rr:
                rows.extend(rr)
                print(f"  {rr[0]['stream']}: {len(rr)} nights, median PWV {np.median([r['pwv'] for r in rr]):.1f} mm,"
                      f" line bias {100 * np.median([r['b_line'] for r in rr]):+.3f} %", flush=True)
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(Q15 / "q15c_wv_nights.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"-> {Q15 / 'q15c_wv_nights.csv'} ({len(rows)} nights)")


# ---------------------------------------------------------------------------------------------------
# summary: factor 4 (analytic) + table + figures
# ---------------------------------------------------------------------------------------------------
def med_iqr(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan, np.nan, 0
    return float(np.median(x)), float(np.percentile(x, 25)), float(np.percentile(x, 75)), int(x.size)


def fmt(x, nd=1, sign=True):
    m, lo, hi, n = med_iqr(x)
    s = "+" if sign else ""
    return f"{m:{s}.{nd}f} [{lo:{s}.{nd}f}, {hi:{s}.{nd}f}]"


def factor4():
    """Molecular backscatter / refractive index / rotational Raman for a wavelength shift (analytic)."""
    # beta_mol ~ lambda^-(4+eps); eps ~ 0.08 at 1064 nm (dispersion of n-1 and of the King factor, Bucholtz 1995)
    expo = 4.08
    out = {}
    for lab, dl in (("air-vs-vacuum nominal (0.29 nm)", 0.29), ("filter detuning 0.19 W (0.19 nm)", 0.19),
                    ("laser cavity tuning, 1 K LOM", 0.013), ("laser gain centre, 1 K LOM", 0.0042)):
        out[lab] = -expo * dl / 1064.0 * 100.0
    # refractive index of air: d(n-1)/dlambda ~ -1.6e-6 /um at 1 um -> path / range error negligible
    out["refractive index (range scale, 0.3 nm)"] = (1.6e-6 * 0.0003) * 100
    return out


def rr_leakage():
    """Rotational-Raman leakage through a 1-nm Gaussian filter detuned 0.19 W blue of the laser line,
    as a fraction of the total (Cabannes + RR) molecular backscatter, and its temperature dependence."""
    h, cl, kB = 6.62607e-34, 2.99792458e10, 1.380649e-23
    tot = {}
    for T in (230.0, 260.0, 290.0):
        frac_pass, frac_rr = 0.0, 0.0
        for gas, B, gI, x, gam2 in (("N2", 1.98957, (6, 3), 0.7808, 0.51), ("O2", 1.43768, (0, 1), 0.2095, 1.26)):
            # anisotropy weights (relative); Stokes / anti-Stokes S-branch lines
            Js = np.arange(0, 60)
            for J in Js:
                g = gI[J % 2]
                if g == 0:
                    continue
                E = B * J * (J + 1) * h * cl
                pop = g * (2 * J + 1) * np.exp(-E / (kB * T))
                for sgn, dJ in ((+1, 2), (-1, -2)):
                    Jf = J + dJ
                    if Jf < 0:
                        continue
                    if dJ == 2:
                        pt = 3 * (J + 1) * (J + 2) / (2 * (2 * J + 1) * (2 * J + 3))
                        dnu = -B * (4 * J + 6)
                    else:
                        pt = 3 * J * (J - 1) / (2 * (2 * J + 1) * (2 * J - 1))
                        dnu = B * (4 * J - 2)
                    wgt = x * gam2 * pop * pt
                    dlam = -dnu * (LAMBDA_VAC ** 2) * 1e-7          # nm (Stokes = red = +)
                    off = dlam + DELTA0_W * FILTER_FWHM               # relative to filter centre
                    frac_rr += wgt
                    frac_pass += wgt * np.exp(-4 * LN2 * (off / FILTER_FWHM) ** 2)
        tot[T] = frac_pass / frac_rr
    # RR share of the total molecular backscatter at 180 deg ~ 3.6 % (Behrendt & Nakamura 2002, Cabannes vs total)
    rr_share = 0.036
    return {T: 100 * rr_share * v for T, v in tot.items()}, rr_share


def run_summary():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    H = list(csv.DictReader(open(Q15 / "q15c_hk_units.csv", encoding="utf-8")))

    def col(name, rows=None, excl_pay=True):
        rows = H if rows is None else rows
        return np.array([float(r[name]) if r.get(name, "") not in ("", "nan") else np.nan for r in rows
                         if not (excl_pay and r["stream"] == PAYERNE)], float)
    pay = [r for r in H if r["stream"] == PAYERNE]
    lines = []
    P = lambda *a: lines.append(" ".join(str(x) for x in a))  # noqa: E731

    # ---------------- factor 1 ----------------
    drivers = ["lin_in", "lin_m2"] + [f"Tint_ema{int(t)}h" for t in TAUS_H] + \
        ["Tint_prevmax", "Text_ema0h", "Text_ema24h", "Text_ema96h", "Tlom", "Tdet", "mix", "fastslow24h", "fastslow96h",
         "fastslow192h", "gauss_q11", "lin_q11",
         "gauss_gen", "gauss_ema96"]
    tab1 = []
    amp_op = col("amp_op")
    for d in drivers:
        if f"{d}_removed" not in H[0]:
            continue
        rem = col(f"{d}_removed")
        tab1.append(dict(driver=d, n=int(np.isfinite(rem).sum()), kappa=med_iqr(col(f"{d}_kappa"))[0],
                         r2=med_iqr(col(f"{d}_r2"))[0], amp_pred=med_iqr(col(f"{d}_amp_pred"))[0],
                         peak_pred=med_iqr(col(f"{d}_peak_pred"))[0], amp_after=med_iqr(col(f"{d}_amp_after"))[0],
                         removed=med_iqr(rem)[0], removed_p25=med_iqr(rem)[1], removed_p75=med_iqr(rem)[2],
                         n2n_sd=med_iqr(col(f"{d}_n2n_sd"))[0], diurnal=med_iqr(col(f"{d}_diurnal"))[0],
                         diurnal_jja=med_iqr(col(f"{d}_diurnal_jja"))[0],
                         driver_amp=med_iqr(col(f"{d}_amp"))[0] if f"{d}_amp" in H[0] else np.nan,
                         driver_peak=med_iqr(col(f"{d}_peak"))[0] if f"{d}_peak" in H[0] else np.nan,
                         kneeded=med_iqr(col(f"{d}_kneeded"))[0] if f"{d}_kneeded" in H[0] else np.nan))
    with open(Q15 / "q15c_factor1_drivers.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tab1[0].keys()))
        w.writeheader()
        w.writerows(tab1)
    T1 = {r["driver"]: r for r in tab1}

    # ---------------- factor 2 ----------------
    k_in = col("kappa_in")
    frac_tune = {lab: a / ALPHA_W_FILTER / FILTER_FWHM for lab, a in LASER_TUNING.items()}   # alpha_L / alpha_f
    s_cav = -k_in * frac_tune["cavity mode (microchip)"]
    s_gain = -k_in * frac_tune["gain centre"]
    lom_amp = col("lom_night_amp")
    lom_sd = col("lom_night_anom_sd")
    lom_dexc = col("lom_day_exc_mean")
    lom_dexc_j = col("lom_day_exc_jja")
    b7 = col("b_Tlom_q7")

    # ---------------- factor 3 ----------------
    WVR = list(csv.DictReader(open(Q15 / "q15c_wv_nights.csv", encoding="utf-8")))
    wv_units = {}
    for r in WVR:
        wv_units.setdefault(r["stream"], []).append(r)
    lam_cols = sorted([k for k in WVR[0] if k.startswith("b_10")], key=lambda s: float(s[2:]))
    wv_rows = []
    for st, rr in wv_units.items():
        t = np.array([float(r["t"]) for r in rr])
        o = np.argsort(t)
        t = t[o]
        g = {k: np.array([float(rr[i][k]) for i in o]) for k in ["pwv", "col_below_zb", "b_line", "b_band",
                                                                   "b_cont_central", "b_cont_upper"] + lam_cols}
        if len(t) < 60:
            continue
        row = dict(stream=st, site=rr[0]["site"], n=len(t), pwv_med=np.median(g["pwv"]),
                   col_below_zb_med=np.median(g["col_below_zb"]))
        row["pwv_amp"], row["pwv_peak"] = harmonic(t, g["pwv"])
        for k in ["b_line", "b_band", "b_cont_central", "b_cont_upper"] + lam_cols:
            y = 100 * g[k]
            row[f"{k}_mean"] = float(np.mean(y))
            row[f"{k}_amp"], row[f"{k}_peak"] = harmonic(t, y)
            row[f"{k}_n2n"] = float(np.nanstd(anomaly(t, y)))
        tot = 100 * (g["b_line"] + g["b_cont_central"])
        row["total_mean"] = float(np.mean(tot))
        row["total_amp"], row["total_peak"] = harmonic(t, tot)
        row["total_n2n"] = float(np.nanstd(anomaly(t, tot)))
        tot_u = 100 * (g["b_line"] + g["b_cont_upper"])
        row["totalU_amp"], _ = harmonic(t, tot_u)
        row["totalU_n2n"] = float(np.nanstd(anomaly(t, tot_u)))
        # share of the night-to-night C_L variance it could explain (correlate with the ln C anomaly)
        lc = np.log(np.array([float(rr[i]["c"]) for i in o]))
        la = anomaly(t, lc)
        ok = np.isfinite(la) & np.isfinite(tot)
        row["corr_lnC_anom"] = float(np.corrcoef(anomaly(t, tot)[ok], la[ok])[0, 1]) if ok.sum() > 30 else np.nan
        wv_rows.append(row)
    with open(Q15 / "q15c_wv_units.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(wv_rows[0].keys()))
        w.writeheader()
        w.writerows(wv_rows)
    WU = [r for r in wv_rows if r["stream"] != PAYERNE]
    wcol = lambda k: np.array([r[k] for r in WU], float)  # noqa: E731
    lut = build_cross_sections()

    # ---------------- factor 4 ----------------
    f4 = factor4()
    rr, rr_share = rr_leakage()

    # ---------------- table ----------------
    lin = T1["lin_in"]
    m2 = T1.get("lin_m2", lin)
    best_ema = max([d for d in T1 if d.startswith("Tint_ema")], key=lambda d: T1[d]["removed"])
    be = T1[best_ema]
    peak_c = med_iqr(col("peak_op"))[0]
    rows_md = []
    rows_md.append(f"| 1a. Filter passband vs laser line, linear κ × night T_int (intra-night κ) | + (C_L rises with T_int) | "
                   f"κ = {fmt(k_in, 2)} %/K | {lin['n2n_sd']:.1f} (r² of κ·ΔT on ln C anomaly: see 1b) | "
                   f"+{lin['diurnal']:.1f} (JJA +{lin['diurnal_jja']:.1f}) | {lin['amp_pred']:.1f} (removes "
                   f"{100 * lin['removed']:.0f} % of {med_iqr(amp_op)[0]:.1f}) | {lin['peak_pred']:.0f} (C_L {peak_c:.0f}) | "
                   f"measured, {lin['n']} units | 1 |")
    rows_md.append(f"| 1b. Same, κ re-estimated on night-to-night anomalies of the best lagged driver ({best_ema}) | + | "
                   f"κ = {be['kappa']:+.2f} %/K (r² {be['r2']:.2f}) | {be['n2n_sd']:.1f} | +{be['diurnal']:.1f} | "
                   f"{be['amp_pred']:.1f} (removes {100 * be['removed']:.0f} %, i.e. "
                   f"{100 * (be['removed'] - T1['Tint_ema0h']['removed']):+.0f} pts vs τ = 0 on the same κ method) | "
                   f"{be['peak_pred']:.0f} | measured, {be['n']} units | 2 |")
    mx = T1.get("mix")
    if mx:
        rows_md.append(f"| 1c. T_int + temp_ext mix | + / {np.nanmedian(col('mix_k_ext')):+.2f} %/K ext | "
                       f"κ_int {np.nanmedian(col('mix_k_int')):+.2f}, κ_ext {np.nanmedian(col('mix_k_ext')):+.2f} %/K | "
                       f"{mx['n2n_sd']:.1f} | – | {mx['amp_pred']:.1f} (removes {100 * mx['removed']:.0f} %) | "
                       f"{mx['peak_pred']:.0f} | measured, {mx['n']} units | 3 |")
    gq = T1.get("gauss_q11")
    gg = T1.get("gauss_gen")
    rows_md.append(f"| 1d. Gaussian (non-linear) passband instead of linear | concave: less gain at the warm end | "
                   f"½d²lnC/dT² = −4ln2(α/W)² = {-4 * LN2 * ALPHA_W_FILTER ** 2 * 100:.2f} %/K² | ≈ same as 1a | "
                   f"≈ 1a | generic: {gg['amp_pred']:.1f} (removes {100 * gg['removed']:.0f} % vs "
                   f"{100 * m2['removed']:.0f} % linear multi-day κ); q11 units: {gq['removed'] * 100 if gq else np.nan:.0f} % vs "
                   f"{100 * T1['lin_q11']['removed'] if 'lin_q11' in T1 else np.nan:.0f} % | {gg['peak_pred']:.0f} | "
                   f"model on measured T, {gg['n']} units | 4 |")
    rows_md.append(f"| 2. Laser line drift when the LOM leaves regulation | − for Δ0>0 (laser moves further red of the "
                   f"filter); empirical sign unit-specific | model −{frac_tune['cavity mode (microchip)']:.2f}·κ = "
                   f"{np.nanmedian(s_cav):+.2f} %/K (cavity) / {np.nanmedian(s_gain):+.2f} %/K (gain); "
                   f"q7 night {np.nanmedian(b7):+.2f}, day −2.0…+4.2 %/K | "
                   f"{np.nanmedian(np.abs(s_cav) * lom_sd):.2f} | {np.nanmedian(s_cav * lom_dexc):+.2f} (JJA "
                   f"{np.nanmedian(s_cav * lom_dexc_j):+.2f}; hot units to ±{np.nanpercentile(np.abs(s_cav * lom_dexc_j), 90):.1f}) | "
                   f"{np.nanmedian(np.abs(s_cav) * lom_amp):.2f} | {med_iqr(col('lom_night_peak'))[0]:.0f} (opposite sign) | "
                   f"measured excursions, model slope | 5 |")
    rows_md.append(f"| 3. H2O absorption at 1064 nm, two-way to the calibration window (not corrected for 1064 nm) | "
                   f"− (more PWV → lower C_L) | line {np.median(wcol('b_line_mean')) / np.median(wcol('pwv_med')):+.4f} %/mm, "
                   f"+continuum ≈ {np.median(wcol('total_mean')) / np.median(wcol('pwv_med')):+.4f} %/mm (upper "
                   f"{(np.median(wcol('b_line_mean')) + np.median(wcol('b_cont_upper_mean'))) / np.median(wcol('pwv_med')):+.4f}) | "
                   f"{np.median(wcol('total_n2n')):.3f} (upper {np.median(wcol('totalU_n2n')):.3f}) | "
                   f"≈0 (BL PWV day≈night) | {np.median(wcol('total_amp')):.3f} (upper {np.median(wcol('totalU_amp')):.3f}) | "
                   f"C_L effect max at {np.median(wcol('total_peak')):.0f}, i.e. C_L MIN at the PWV peak (DOY {np.median(wcol('pwv_peak')):.0f}): anti-phase | HITRAN+CAMS, {len(WU)} units | 6 |")
    rows_md.append(f"| 4. Molecular λ⁻⁴ / refractive index / rotational Raman for a wavelength shift | − (red shift → "
                   f"smaller β_mol) | {f4['laser cavity tuning, 1 K LOM']:+.4f} %/K LOM (cavity), "
                   f"{f4['filter detuning 0.19 W (0.19 nm)']:+.3f} % per 0.19 nm; RR leakage "
                   f"{rr[290.0]:.3f} % of β_mol | <0.001 | <0.01 | <0.005 | – | analytic | 7 |")

    # ---------------- figures ----------------
    fig, axs = plt.subplots(1, 3, figsize=(20, 5.6))
    ax = axs[0]
    taus = [t for t in TAUS_H if f"Tint_ema{int(t)}h" in T1]
    rem = [100 * T1[f"Tint_ema{int(t)}h"]["removed"] for t in taus]
    kk = [T1[f"Tint_ema{int(t)}h"]["kappa"] for t in taus]
    xt = [max(t, 1.0) for t in taus]
    ax.semilogx(xt, rem, "o-", color="#d62728", lw=2, label="T_int EMA: annual amplitude removed")
    for d, c_, lab in (("Text_ema0h", "#1f77b4", "temp_ext"), ("Text_ema24h", "#1f77b4", None),
                       ("Text_ema96h", "#1f77b4", None)):
        pass
    te = [(max(float(d.split("ema")[1][:-1]), 1.0), 100 * T1[d]["removed"]) for d in ("Text_ema0h", "Text_ema24h", "Text_ema96h")
          if d in T1]
    ax.semilogx([a for a, b in te], [b for a, b in te], "s--", color="#1f77b4", lw=1.5, label="temp_ext EMA")
    ax.axhline(100 * lin["removed"], color="k", ls=":", lw=1.2, label=f"intra-night κ × T_int ({100 * lin['removed']:.0f} %)")
    ax.axhline(100 * m2["removed"], color="grey", ls="--", lw=1.2, label=f"multi-day κ × T_int ({100 * m2['removed']:.0f} %)")
    ax.set_xlabel("thermal lag τ of the driver [h] (τ=0 plotted at 1 h)")
    ax.set_ylabel("share of the C_L annual amplitude removed [%] (network median)")
    ax2 = ax.twinx()
    ax2.semilogx(xt, kk, "^:", color="#2ca02c", lw=1)
    ax2.set_ylabel("anomaly κ of the lagged T_int [%/K]", color="#2ca02c")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("(a) factor 1: lagged / alternative driver temperatures", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[1]
    pk_c = col("peak_op")
    for d, colr, lab in (("Tint_ema0h", "#d62728", "T_int night"), (best_ema, "#ff7f0e", best_ema),
                         ("Text_ema0h", "#1f77b4", "temp_ext night"), ("Tint_prevmax", "#9467bd", "T_int prev. afternoon max")):
        if f"{d}_peak" in H[0]:
            ax.hist(col(f"{d}_peak"), bins=np.arange(100, 280, 6), histtype="step", lw=2, color=colr, label=lab)
    ax.hist(pk_c, bins=np.arange(100, 280, 6), histtype="stepfilled", alpha=0.35, color="k", label="C_L (Kalman fit)")
    ax.set_xlabel("day of year of the annual maximum")
    ax.set_ylabel("number of units")
    ax.legend(fontsize=8)
    ax.set_title("(b) phase: C_L peaks ~3-4 weeks before the temperatures", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[2]
    x = col("lin_in_amp_pred")
    y = amp_op
    ax.scatter(x, y, s=18, color="#7f7f7f", label="κ_in × T_int (predicted)")
    ax.scatter(col(f"{best_ema}_amp_pred"), y, s=18, color="#ff7f0e", label=f"{best_ema} (anomaly κ)")
    ax.scatter(col("gauss_gen_amp_pred"), y, s=12, marker="x", color="#2ca02c", label="Gaussian passband (multi-day κ)")
    lim = [0, np.nanpercentile(y, 98)]
    ax.plot(lim, lim, "k--", lw=0.8)
    ax.set_xlabel("predicted annual amplitude of ln C_L [%]")
    ax.set_ylabel("observed annual amplitude of C_L [%] (Kalman sinusoid)")
    ax.set_title("(c) predicted vs observed annual amplitude per unit", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q15 / "q15c_factor1_drivers.png", dpi=120)
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(20, 5.6))
    ax = axs[0]
    sp = lut["spec_ref"]
    o = np.argsort(sp[:, 0])
    ax.semilogy(sp[o, 0], sp[o, 1], color="#1f77b4", lw=0.8, label="H2O σ (HITRAN, 900 hPa, 285 K)")
    ax.axvline(LAMBDA_VAC, color="#d62728", lw=1.5, label=f"Nd:YAG {LAMBDA_VAC} nm (vac) = 1064.20 nm air")
    ax.axvspan(LAMBDA_VAC - 0.05, LAMBDA_VAC + 0.15, color="#d62728", alpha=0.12, label="plausible laser line range")
    lc = LAMBDA_VAC - DELTA0_W * FILTER_FWHM
    lamf = np.linspace(1062.5, 1066.5, 400)
    ax3 = ax.twinx()
    ax3.plot(lamf, np.exp(-4 * LN2 * ((lamf - lc) / FILTER_FWHM) ** 2), color="k", ls="--", lw=1)
    ax3.set_ylabel("filter transmission (1 nm Gaussian, Δ0 = 0.19 W)")
    ax.set_xlim(1062.5, 1066.5)
    ax.set_xlabel("vacuum wavelength [nm]")
    ax.set_ylabel("cross-section [cm² molecule⁻¹]")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("(a) H2O lines around the 1064 nm line", fontsize=10)
    ax = axs[1]
    lams = [float(k[2:]) for k in lam_cols]
    medb = [np.median(wcol(f"{k}_mean")) for k in lam_cols]
    ampb = [np.median(wcol(f"{k}_amp")) for k in lam_cols]
    ax.semilogy(lams, np.abs(medb), "o-", color="#1f77b4", label="|mean two-way bias| [%]")
    ax.semilogy(lams, ampb, "s-", color="#d62728", label="annual amplitude [%]")
    ax.axvline(LAMBDA_VAC, color="k", ls=":")
    ax.set_xlabel("assumed laser line (vacuum) [nm]")
    ax.set_ylabel("C_L bias from H2O lines [%] (network median)")
    ax.legend(fontsize=8)
    ax.set_title("(b) sensitivity to the exact line position", fontsize=10)
    ax.grid(alpha=0.3, which="both")
    ax = axs[2]
    pwv = np.array([float(r["pwv"]) for r in WVR if r["stream"] != PAYERNE])
    bl = np.array([100 * (float(r["b_line"]) + float(r["b_cont_central"])) for r in WVR if r["stream"] != PAYERNE])
    bu = np.array([100 * (float(r["b_line"]) + float(r["b_cont_upper"])) for r in WVR if r["stream"] != PAYERNE])
    ax.scatter(pwv, bl, s=3, alpha=0.3, color="#1f77b4", label="line + central continuum")
    ax.scatter(pwv, bu, s=3, alpha=0.3, color="#ff7f0e", label="line + upper continuum")
    ax.set_xlabel("CAMS precipitable water vapour [mm]")
    ax.set_ylabel("ln C_L bias [%]")
    ax.legend(fontsize=8, markerscale=4)
    ax.set_title(f"(c) per calibration night ({len(pwv)} nights, {len(WU)} units)", fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q15 / "q15c_factor3_watervapour.png", dpi=120)
    plt.close(fig)

    # ---------------- markdown ----------------
    md = []
    md.append("# q15c – spectral / wavelength factors of the CHM15k C_L variability (1064 nm)\n")
    md.append("| factor | sign | sensitivity | night-to-night [%] | diurnal [%] | annual amplitude [%] | phase DOY | evidence | rank |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    md += rows_md
    md.append("")
    md.append(f"Network medians over {len(H) - len(pay)} CHM15k units (main optical module, 2025–2026; Payerne "
              f"0-20000-0-06610_A excluded). Observed: C_L annual amplitude {fmt(amp_op)} %, peak DOY "
              f"{fmt(pk_c, 0, False)}; night-to-night sd of ln C_L about its 31-day median {fmt(col('sd_anom_lnC'))} %.")
    md.append("")
    md.append("## Factor 1 – driver comparison (network medians)\n")
    md.append("Each driver's coefficient is fitted on 31-day anomalies (night-to-night only); its predicted cycle is then "
              "removed from C_L and the paper Kalman sinusoid refitted (`annual` of q13).\n")
    md.append("| driver | κ anomaly [%/unit] | r² (n2n) | driver annual amp | driver peak DOY | κ needed for 100 % | "
              "predicted amp [%] | peak DOY | amp after [%] | removed [%] (IQR) | n2n sd [%] | diurnal [%] (JJA) |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in tab1:
        md.append(f"| {r['driver']} | {r['kappa']:+.2f} | {r['r2']:.3f} | {r['driver_amp']:.2f} | {r['driver_peak']:.0f} | "
                  f"{r['kneeded']:.2f} | {r['amp_pred']:.1f} | {r['peak_pred']:.0f} | {r['amp_after']:.1f} | "
                  f"{100 * r['removed']:.0f} ({100 * r['removed_p25']:.0f}…{100 * r['removed_p75']:.0f}) | {r['n2n_sd']:.2f} | "
                  f"{r['diurnal']:+.1f} ({r['diurnal_jja']:+.1f}) |")
    md.append("")
    dphi = col("peak_lnC") - col("Tint_peak")
    cap = 1 - np.abs(np.sin(W_DOY * dphi))
    kin_phase = col("amp_lnC") * np.cos(W_DOY * dphi) / col("Tint_amp")
    dphi_e = col("peak_lnC") - col(f"{best_ema}_peak")
    md.append(f"**Phase limit.** ln C_L peaks at DOY {fmt(col('peak_lnC'), 0, False)}, the calibration-night T_int at "
              f"{fmt(col('Tint_peak'), 0, False)} (lag {fmt(-dphi, 0, False)} d; {best_ema}: {fmt(-dphi_e, 0, False)} d). "
              f"Any driver proportional to a temperature that peaks in late July can remove at most 1 − |sin Δφ| = "
              f"{fmt(100 * cap, 0, False)} % of the C_L annual amplitude, whatever its coefficient. The in-phase part alone "
              f"would need κ = {fmt(kin_phase, 2, False)} %/K on T_int, vs {fmt(k_in, 2, False)} (intra-night) and "
              f"{fmt(col('kappa_m2'), 2, False)} (multi-day) measured.\n")
    md.append(f"**Diurnal definition.** Diurnal = κ × [T at 12–16 h solar − T in the calibration window of the night ending "
              f"that morning]: ΔT = {fmt(col('lin_in_diurnal') / k_in, 1, False)} K (JJA {fmt(col('lin_in_diurnal_jja') / k_in, 1, False)} K). "
              f"Report 19 §4.4 (noon anomaly of the mean diurnal cycle, 3.7 K) gives +4.8 % (JJA +7.5 %). A lagged driver "
              f"(τ ≥ 12 h) predicts NO daytime bias, contradicting the measured day−night pair ratios (report 19 §4.2–4.4, r = 0.98 "
              f"with the τ = 0 prediction) and the intra-night lag scan (q11, best τ ≈ 0).\n")
    md.append("![Factor 1](q15c_factor1_drivers.png)\n")
    md.append("## Factor 2 – LOM regulation (laser line)\n")
    md.append(f"* Night LOM set-point dither: fraction of calibration nights > 0.5 K off set-point {fmt(col('lom_night_frac_out'), 3, False)}; "
              f"annual amplitude of the night excursion {fmt(lom_amp, 2, False)} K, anomaly sd {fmt(lom_sd, 2, False)} K.")
    md.append(f"* Afternoon (12–18 h solar) out-of-regulation fraction {fmt(col('lom_day_frac_out'), 3, False)} (JJA "
              f"{fmt(col('lom_day_frac_out_jja'), 3, False)}); mean positive excursion {fmt(lom_dexc, 2, False)} K (JJA "
              f"{fmt(lom_dexc_j, 2, False)} K, JJA p95 {fmt(col('lom_day_exc_p95_jja'), 1, False)} K).")
    md.append(f"* Slope model: with Δ0/W = κ/(8 ln2 α_f/W), a laser tuning α_L gives d lnC/dT_lom = −κ α_L/α_f: "
              f"{np.nanmedian(s_cav):+.2f} %/K (13 pm/K cavity) or {np.nanmedian(s_gain):+.2f} %/K (4.2 pm/K gain centre). "
              f"Measured (q7): night partial slope {fmt(b7, 2)} %/K; daytime Aosta +4.2, Budapest +1.4, Măgurele +0.4, "
              f"Payerne −2.0 %/K (mode hops / unit-specific detuning).")
    md.append("")
    md.append("## Factor 3 – water vapour at 1064 nm\n")
    md.append(f"* HITRAN 2020 H2O (6 isotopologues, 9330–9460 cm⁻¹, Voigt, 25 cm⁻¹ cut-off): σ at 1064.49 nm vac = "
              f"{lut['sig'][int(np.argmin(np.abs(lut['lam_scan'] - LAMBDA_VAC))), 5, 4]:.2e} cm² (900 hPa, 285 K), dominated by "
              f"the weak line at 1064.494 nm (S = 4.3e-27); the nearest strong line (S = 1.1e-25) is at 1064.736 nm vac "
              f"(+0.25 nm), the strongest within the 1-nm filter at 1064.68 nm vac. Filter-weighted band σ ≈ "
              f"{lut['sig_band'][5, 4]:.2e} cm² (relevant for the solar background only, not for the elastic signal).")
    md.append(f"* Continuum (not in HITRAN lines): MT_CKD-type self/foreign coefficients bounded at Cs = 2e-26…1e-25, "
              f"Cf = 2e-27…1e-26 cm² molec⁻¹ atm⁻¹ (literature, no measurement at 1.06 µm available here).")
    md.append(f"* CAMS 0.4° L137 per calibration night, {len(WU)} geographically spread units, {len(pwv)} nights: PWV median "
              f"{fmt(wcol('pwv_med'), 1, False)} mm, annual amplitude {fmt(wcol('pwv_amp'), 1, False)} mm, peak DOY "
              f"{fmt(wcol('pwv_peak'), 0, False)}.")
    md.append(f"* Two-way bias of ln C_L (window-averaged): line {fmt(wcol('b_line_mean'), 3)} %, band {fmt(wcol('b_band_mean'), 3)} %, "
              f"continuum central {fmt(wcol('b_cont_central_mean'), 3)} % / upper {fmt(wcol('b_cont_upper_mean'), 3)} %. "
              f"Annual amplitude (line+central) {fmt(wcol('total_amp'), 3, False)} %, night-to-night sd {fmt(wcol('total_n2n'), 3, False)} %, "
              f"correlation with the ln C_L anomaly {fmt(wcol('corr_lnC_anom'), 2)}.")
    worst = max(lam_cols, key=lambda k: np.median(wcol(f"{k}_amp")))
    md.append(f"* Line position scan 1064.29–1064.71 nm vac: worst case at {worst[2:]} nm, annual amplitude "
              f"{np.median(wcol(f'{worst}_amp')):.2f} % (mean bias {np.median(wcol(f'{worst}_mean')):+.2f} %); "
              f"within the plausible laser range 1064.44–1064.64 nm the amplitude stays ≤ "
              f"{max(np.median(wcol(f'{k}_amp')) for k in lam_cols if 1064.44 <= float(k[2:]) <= 1064.64):.3f} %.")
    md.append("")
    md.append("![Factor 3](q15c_factor3_watervapour.png)\n")
    md.append("## Factor 4 – molecular λ⁻⁴, refractive index, rotational Raman\n")
    for k, v in f4.items():
        md.append(f"* {k}: {v:+.4f} % of β_mol")
    md.append(f"* Rotational-Raman leakage through the 1-nm filter (Gaussian, detuned 0.19 W blue), RR share of molecular "
              f"backscatter {100 * rr_share:.1f} %: transmitted RR = " + ", ".join(f"{v:.3f} % at {T:.0f} K" for T, v in rr.items())
              + " of the total β_mol. If the molecular model is the total (Cabannes + RR) Bucholtz cross-section, it over-predicts the received "
              f"molecular signal by ≈ {100 * rr_share - rr[290.0]:.1f} % – a static offset, temperature dependence < 0.01 %.")
    md.append("")
    md.append("## Payerne (excluded from medians)\n")
    for r in pay:
        md.append(f"* {r['stream']}: amp_op {float(r['amp_op']):.1f} %, lin_in removed {100 * float(r.get('lin_in_removed') or 'nan'):.0f} %, "
                  f"{best_ema} removed {100 * float(r.get(f'{best_ema}_removed') or 'nan'):.0f} %, LOM afternoon out-of-regulation "
                  f"(JJA) {float(r.get('lom_day_frac_out_jja') or 'nan'):.2f}, JJA mean excursion {float(r.get('lom_day_exc_jja') or 'nan'):.2f} K.")
    md.append("")
    md.append("## Example units\n")
    md.append("| unit | site | amp_op | peak | κ_in | lin_in removed | " + best_ema + " removed | Text removed | gauss_gen removed |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    ex = ["0-20000-0-06240_A", "0-20000-0-06240_B", "0-20000-0-12843_A", "0-20008-0-INO_A", "0-20000-0-16054_A",
          "0-20000-0-10393_A", "0-20000-0-03808_A"]
    for r in H:
        if r["stream"] in ex:
            g = lambda k: float(r.get(k) or "nan")  # noqa: E731
            md.append(f"| {r['stream']} | {r['site']} | {g('amp_op'):.1f} | {g('peak_op'):.0f} | {g('kappa_in'):+.2f} | "
                      f"{100 * g('lin_in_removed'):.0f} % | {100 * g(f'{best_ema}_removed'):.0f} % | "
                      f"{100 * g('Text_ema0h_removed'):.0f} % | {100 * g('gauss_gen_removed'):.0f} % |")
    md.append("")
    md += ["## Notes", ""] + [f"* {x}" for x in lines]
    (Q15 / "q15c_summary.md").write_text("\n".join(md), encoding="utf-8")
    sys.stdout.reconfigure(encoding="utf-8")
    print("\n".join(md[:14]))
    print(f"-> {Q15 / 'q15c_summary.md'}")


if __name__ == "__main__":
    part = sys.argv[1] if len(sys.argv) > 1 else "all"
    if part in ("hk", "all"):
        run_hk()
    if part in ("wv", "all"):
        run_wv()
    if part in ("summary", "all"):
        run_summary()
