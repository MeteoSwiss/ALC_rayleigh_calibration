# -*- coding: utf-8 -*-
"""q15a: can the CHM15k LASER TRANSMITTER explain the annual cycle of the Rayleigh constant?

Part of the q15 attribution study (brief: C:/DATA/Projects/202610_CHM15k_stability/q15/BRIEF.md).
The CHM15k laser is a diode-pumped, passively Q-switched Nd:YAG microchip laser in the heater-only
laser-optical module (LOM). The firmware divides the summed raw signal by the number of pulses of the
record, so the normalised signal (and C_L) is proportional to the energy PER PULSE. For a passively
Q-switched laser the average power is E x PRF; at constant absorbed pump power, E ~ 1/PRF. The pulse
count per record (`laser_pulses`) is therefore the only transmitter monitor in the L1 files.

Factors quantified per unit (main optical-module era of the 2025-2026 nightly constants):
  1. pulse energy vs temperature: d ln(pulses)/dT_int on 31-day anomalies -> implied energy coefficient
     (E ~ 1/PRF), plus the LOM regulation loss (temperature_optical_module above set-point);
  2. laser ageing: linear trend of ln C_L per era (%/yr and %/1000 h of laser_life_time) and how much
     apparent annual amplitude a trend produces in the Kalman annual-sinusoid fit over ~18 months
     (unit's own sampling); 11-year archive (13 stations, serlom eras) for multi-year ageing rates and an
     18-month sliding-window aliasing test on real data;
  3. pump-diode temperature: bounded from the PRF response to T_int and to LOM regulation loss;
  4. pulses-per-record normalisation: seasonal cycle of laser_pulses and its predicted C_L effect;
  5. laser servicing events: laser_life_time discontinuities (16-bit counter WRAPS at 65536 h, which is
     not an event), module swaps, laser-off gaps; C_L steps at events and their aliasing into the annual fit.

Run:  .venv/Scripts/python.exe validation/chm15k_stability/q15a_laser.py
Outputs: C:/DATA/Projects/202610_CHM15k_stability/q15/q15a_*.csv, q15a_*.png
"""
from __future__ import annotations

import csv
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from validation.paper.kalman_improve_alc import fit_model, flag_outliers  # noqa: E402
from common import (BLOCKS_PER_DAY, OUT, chm15k_streams, load_cache, load_constants,  # noqa: E402
                    stream_key)
from q1_temperature import ANOM_WINDOW, night_hk, ols, robust_mask, running_median  # noqa: E402

warnings.filterwarnings("ignore")
Q15 = OUT / "q15"
EPOCH = datetime(1970, 1, 1)
PAYERNE = "0-20000-0-06610_A"
WRAP = 65536.0                  # laser_life_time is a 16-bit hour counter
MIN_NIGHTS = 100
N_WORKERS = 4
KAPPA_LOM_NIGHT = 0.30          # %/K, q7 network partial slope on temp_lom at night (report 19 sec. 7.3)
W = 2 * np.pi / 365.25


# ---------------------------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------------------------
def harmonic(t, y, trend=False):
    """OLS y = a + s sin + c cos (+ b t_yr). Returns amp (in y units), peak DOY, trend/yr, se trend, resid."""
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if len(y) < 30:
        return np.nan, np.nan, np.nan, np.nan, None
    doy = W * (t % 365.25)
    cols = [np.ones_like(t), np.sin(doy), np.cos(doy)]
    if trend:
        cols.append((t - t.mean()) / 365.25)
    coef, se, res = ols(np.column_stack(cols), y)
    amp = np.hypot(coef[1], coef[2])
    # peak of s sin(x) + c cos(x) at x = atan2(s, c); DOY = x / W
    peak = (np.arctan2(coef[1], coef[2]) / W) % 365.25
    return amp, peak, (coef[3] if trend else np.nan), (se[3] if trend else np.nan), res


def kalman_amp(t_days, values):
    """Annual amplitude (% of offset) and peak DOY of the paper Kalman model, as in q13 (amp_op)."""
    v = np.asarray(values, float)
    out = flag_outliers(v)
    times = [EPOCH + timedelta(days=float(x)) for x, o in zip(t_days, out) if not o]
    y = v[~out]
    if len(y) < 60:
        return np.nan, np.nan
    _, params, model = fit_model(times, y)
    if model != "annual cycle":
        return np.nan, np.nan
    amp, phase, off = params
    peak = ((np.pi / 2 if amp > 0 else 3 * np.pi / 2) - phase) / (2 * np.pi / 365) % 365
    return 100 * abs(amp) / off, peak


def kalman_amp_noscreen(t_days, values):
    """Same Kalman sinusoid without the IQR screen (for synthetic, outlier-free series)."""
    times = [EPOCH + timedelta(days=float(x)) for x in t_days]
    v = np.asarray(values, float)
    if not np.isfinite(v).all():
        return np.nan, np.nan
    _, params, model = fit_model(times, v)
    amp, phase, off = params
    peak = ((np.pi / 2 if amp > 0 else 3 * np.pi / 2) - phase) / (2 * np.pi / 365) % 365
    return 100 * abs(amp) / off, peak


def unwrap_lifetime(x):
    """Map the 16-bit laser hour counter to [0, 65536) (some files carry it shifted by -65536)."""
    return np.mod(x, WRAP)


# ---------------------------------------------------------------------------------------------
# Laser life-time events from the 10-min cache
# ---------------------------------------------------------------------------------------------
def lifetime_events(cache, i_lt):
    """Discontinuities of laser_life_time, using the first valid 10-min block of each cached day.

    The laser runs 24/7 (median counting rate 1.000 h/h on every unit), so the expected increment
    between two samples is the elapsed time. Classified events:
      (wrap)     counter passes 65535 -> 0: bookkeeping, not an event (not listed);
      glitch     the counter leaves its trajectory and returns to it within 3 samples (e.g. a day of 0):
                 a read-out artefact, not a laser change;
      encoding   persistent offset that is a multiple of 4096 h (bit flips 2^12..2^15): decoding artefact;
      new_laser  counter restarts near 0 (< 2000 h) and keeps counting from there -> laser replaced;
      off        increment short of elapsed time by > 48 h -> laser (or instrument) off;
      jump       any other persistent discontinuity > 48 h -> counter of another (used) laser.
    """
    lt = unwrap_lifetime(cache["hk"][:, :, i_lt])
    tb = cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY
    first = np.argmax(np.isfinite(lt), axis=1)
    has = np.isfinite(lt).any(axis=1)
    t = tb[np.arange(len(lt)), first][has]
    v = lt[np.arange(len(lt)), first][has]
    mod = cache["module"][has]

    def on_track(j, i0):
        return abs(((v[j] - v[i0]) % WRAP) - 24.0 * (t[j] - t[i0])) <= 48

    ev = []
    ref, i = 0, 1
    while i < len(v):
        exp_inc = 24.0 * (t[i] - t[ref])
        inc = (v[i] - v[ref]) % WRAP
        dev = inc - exp_inc
        swap = mod[i] != mod[ref]
        if abs(dev) <= 48 and not swap:
            ref, i = i, i + 1
            continue
        rec = dict(t=float(t[i]), lt_before=float(v[ref]), lt_after=float(v[i]), gap_h=float(exp_inc),
                   dev_h=float(dev), swap=bool(swap), mod_before=str(mod[ref]), mod_after=str(mod[i]))
        back = [j for j in range(i + 1, min(i + 4, len(v))) if on_track(j, ref)]
        persist = i + 1 < len(v) and on_track(i + 1, i)
        if not swap and (back or not persist):
            # the counter returns to its old trajectory, or the new value itself is an isolated outlier
            ev.append(dict(rec, kind="glitch"))
            i = back[0] if back else i + 1
            continue
        r4096 = abs(dev) % 4096
        if v[i] < 2000 and abs(dev) > 48:
            kind = "new_laser"
        elif -exp_inc - 48 <= dev < -48:
            kind = "off"
        elif abs(dev) >= 4000 and min(r4096, 4096 - r4096) <= 60:
            kind = "encoding"
        elif abs(dev) > 48:
            kind = "jump"
        else:
            kind = "swap_only"
        ev.append(dict(rec, kind=kind))
        ref, i = i, i + 1
    duty = np.nan
    if len(v) > 10:
        inc = np.diff(v) % WRAP
        el = 24.0 * np.diff(t)
        okd = np.abs(inc - el) < 48
        duty = float(inc[okd].sum() / el[okd].sum()) if okd.any() else np.nan
    return ev, duty


# ---------------------------------------------------------------------------------------------
# Diurnal cycle of pulses / LOM temperature from the 10-min cache
# ---------------------------------------------------------------------------------------------
def diurnal(cache, lon, iv, main_mod):
    """Noon (solar 11-15 h) minus same-date night (solar 0-4 h) per day; medians all / JJA.

    Also the daytime PRF response to LOM regulation loss: OLS of (ln P - night ln P) on
    (T_lom - set-point)+ with (T_int - night T_int) as control, over daytime blocks.
    """
    nd = len(cache["day"])
    tb = cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY
    tloc = tb + lon / 360.0
    ldate = np.floor(tloc).astype(np.int64)
    sh = (tloc - ldate) * 24.0
    use = np.isin(cache["module"], [main_mod, "TUB" + main_mod, main_mod.replace("TUB", "")])
    lp = np.log(cache["hk"][:, :, iv["laser_pulses"]])
    tl = cache["hk"][:, :, iv["temperature_optical_module"]]
    ti = cache["hk"][:, :, iv["temp_int"]]
    setp = np.nanmedian(tl[use][(sh[use] < 4)]) if use.any() else np.nan
    if use.sum() < 30:
        return {}
    flat = {k: a[use].ravel() for k, a in (("lp", lp), ("tl", tl), ("ti", ti), ("sh", sh), ("ld", ldate))}
    nightm = flat["sh"] < 4
    noonm = (flat["sh"] >= 11) & (flat["sh"] < 15)
    dayblk = (flat["sh"] >= 8) & (flat["sh"] < 18)
    # per local date medians
    dates = np.unique(flat["ld"])
    di = np.searchsorted(dates, flat["ld"]).astype(np.int64)
    nmed = {}
    for name in ("lp", "tl", "ti"):
        for lab, msk in (("night", nightm), ("noon", noonm)):
            arr = np.full(len(dates), np.nan)
            vals = flat[name]
            # group medians via sorting
            sel = msk & np.isfinite(vals)
            order = np.argsort(di[sel], kind="stable")
            g = di[sel][order]
            v = vals[sel][order]
            bounds = np.flatnonzero(np.diff(g)) + 1
            for gi, seg in zip(np.split(g, bounds), np.split(v, bounds)):
                if seg.size >= 6:
                    arr[gi[0]] = np.median(seg)
            nmed[(name, lab)] = arr
    month = np.array([(EPOCH + timedelta(days=int(d))).month for d in dates])
    jja = np.isin(month, (6, 7, 8))
    out = {}
    for name, scale in (("lp", 100.0), ("tl", 1.0), ("ti", 1.0)):
        d = nmed[(name, "noon")] - nmed[(name, "night")]
        out[f"d_{name}_noon"] = float(np.nanmedian(d) * scale) if np.isfinite(d).sum() > 20 else np.nan
        out[f"d_{name}_noon_jja"] = float(np.nanmedian(d[jja]) * scale) if np.isfinite(d[jja]).sum() > 10 else np.nan
    # LOM above set-point by day: fraction of daytime blocks > set-point + 0.5 K
    exc = flat["tl"] - setp
    out["lom_setpoint"] = float(setp)
    out["lom_day_exc_frac"] = float(np.nanmean(exc[dayblk] > 0.5)) if dayblk.any() else np.nan
    out["lom_noon_exc_jja"] = float(np.nanmedian(np.clip(nmed[("tl", "noon")][jja] - setp, 0, None))) \
        if np.isfinite(nmed[("tl", "noon")][jja]).sum() > 10 else np.nan
    # PRF response to LOM regulation loss (daytime blocks, referenced to the same date's night)
    sel = dayblk & np.isfinite(flat["lp"]) & np.isfinite(flat["tl"]) & np.isfinite(flat["ti"])
    dlp = flat["lp"][sel] - nmed[("lp", "night")][di[sel]]
    dti = flat["ti"][sel] - nmed[("ti", "night")][di[sel]]
    ex = np.clip(flat["tl"][sel] - setp, 0, None)
    ok = np.isfinite(dlp) & np.isfinite(dti)
    out["dlnP_dTlom_exc"] = out["dlnP_dTlom_exc_se"] = out["dlnP_dTint_day"] = np.nan
    if ok.sum() > 500 and np.mean(ex[ok] > 0.5) > 0.02:
        c, se, _ = ols(np.column_stack([np.ones(ok.sum()), ex[ok], dti[ok]]), dlp[ok])
        out["dlnP_dTlom_exc"], out["dlnP_dTlom_exc_se"], out["dlnP_dTint_day"] = 100 * c[1], 100 * se[1], 100 * c[2]
    return out


# ---------------------------------------------------------------------------------------------
# Per-stream analysis (2025-2026 network)
# ---------------------------------------------------------------------------------------------
def analyse_stream(args):
    wmo, ident, site, lat, lon = args
    key = stream_key(wmo, ident)
    cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
    if cons is None or cache is None:
        return None
    names = list(cache["hk_vars"])
    iv = {n: names.index(n) for n in names}
    nh = night_hk(cache, cons["t0"], cons["t1"])
    T = nh[:, iv["temp_int"]]
    Tl = nh[:, iv["temperature_optical_module"]]
    P = nh[:, iv["laser_pulses"]]
    LT = unwrap_lifetime(nh[:, iv["laser_life_time"]])
    mods = cons["module"]
    main = max(set(mods), key=list(mods).count)
    m = (mods == main) & (cons["c"] > 0) & np.isfinite(T) & np.isfinite(P) & (P > 0)
    if m.sum() < MIN_NIGHTS:
        return None
    t, c, T, Tl, P, LT = cons["t"][m], cons["c"][m], T[m], Tl[m], P[m], LT[m]
    lc = np.log(c)
    keep = robust_mask(lc)                 # per-module screen (one module here)
    lp = np.log(P)
    t, c, lc, T, Tl, lp, LT = (a[keep] for a in (t, c, lc, T, Tl, lp, LT))
    if len(t) < MIN_NIGHTS:
        return None
    # pulses-per-record configuration steps (averaging-time change, factor ~2-4) are not laser physics:
    # blank the nights off the dominant configuration for the pulse analyses only
    poff = np.abs(lp - np.median(lp)) >= 0.1
    lp = np.where(poff, np.nan, lp)
    # laser age: the 16-bit counter wraps every 65536 h (7.5 yr). Wraps inside the record are unwrapped
    # (lower bound); wraps before 2025 are invisible, so an upper bound adds as many wraps as the module's
    # production year allows (laser assumed as old as the module, continuous operation, +1 yr slack).
    LTu = LT + WRAP * np.cumsum(np.r_[0, np.diff(LT) < -WRAP / 2])
    year = 2000 + int(str(main).replace("TUB", "")[:2])
    mod_age_h = (np.median(t) - (datetime(year, 7, 1) - EPOCH).days) * 24.0
    k_wrap = max(0, int(np.floor((mod_age_h + 8766 - np.median(LTu)) / WRAP)))
    r = dict(stream=key, site=site, module=main, n=int(len(t)), pulses_med=float(np.exp(np.nanmedian(lp))),
             frac_pulse_offconfig=float(poff.mean()), module_year=year,
             laser_age_kh=float(np.median(LTu) / 1000), laser_age_max_kh=float((np.median(LTu) + k_wrap * WRAP) / 1000),
             wrap_in_record=bool(np.any(np.diff(LT) < -WRAP / 2)))
    LT = LTu

    # ---- short-term (31-day anomalies) ----
    la_c = lc - running_median(t, lc, ANOM_WINDOW)
    la_p = lp - running_median(t, lp, ANOM_WINDOW)
    ta = T - running_median(t, T, ANOM_WINDOW)
    tla = Tl - running_median(t, Tl, ANOM_WINDOW)
    ok = np.isfinite(la_c) & np.isfinite(la_p) & np.isfinite(ta)
    r["sd_lnC_anom"] = 100 * float(np.std(la_c[ok]))
    r["sd_lnP_anom"] = 100 * float(np.std(la_p[ok]))
    r["sd_lnP_anom_mad"] = 100 * 1.4826 * float(np.median(np.abs(la_p[ok] - np.median(la_p[ok]))))
    r["share_var_P"] = float(np.var(la_p[ok]) / np.var(la_c[ok]))    # predicted share if E ~ 1/PRF
    cf, se, _ = ols(np.column_stack([np.ones(ok.sum()), ta[ok]]), la_p[ok])
    r["dlnP_dT"], r["dlnP_dT_se"] = 100 * cf[1], 100 * se[1]          # %/K
    r["kappa_E_implied"] = -r["dlnP_dT"]                                # %/K if E ~ 1/PRF
    cf, se, _ = ols(np.column_stack([np.ones(ok.sum()), ta[ok], la_p[ok]]), la_c[ok])
    r["g_lnC_lnP"], r["g_lnC_lnP_se"], r["kappa_ctrl"] = cf[2], se[2], 100 * cf[1]
    r["r_lnC_lnP_anom"] = float(np.corrcoef(la_c[ok], la_p[ok])[0, 1])
    okl = ok & np.isfinite(tla)
    r["sd_Tlom_anom"] = float(np.std(tla[okl])) if okl.sum() > 30 else np.nan
    r["sd_lnC_from_Tlom"] = KAPPA_LOM_NIGHT * r["sd_Tlom_anom"]          # % (night slope q7)

    # ---- annual: pulses (predicted C effect = -ln P) and LOM ----
    a, pk, tr, trse, _ = harmonic(t, -lp, trend=True)
    r["amp_pred_P"], r["peak_pred_P"], r["trend_lnP_pct_yr"] = 100 * a, pk, -100 * tr
    a, pk, _, _, _ = harmonic(t, Tl)
    r["Tlom_amp_K"], r["Tlom_peak"] = a, pk
    r["amp_pred_Tlom"] = KAPPA_LOM_NIGHT * a
    setp = np.nanmedian(Tl)
    r["lom_night_exc_frac"] = float(np.nanmean(Tl > setp + 0.5))

    # ---- ageing (trend of ln C in the era, harmonic controlled) ----
    a0, pk0, _, _, _ = harmonic(t, lc)
    a1, pk1, tr1, se1, _ = harmonic(t, lc, trend=True)
    r["amp_harm"], r["peak_harm"] = 100 * a0, pk0
    r["amp_harm_trend"], r["peak_harm_trend"] = 100 * a1, pk1
    r["trend_pct_yr"], r["trend_se"] = 100 * tr1, 100 * se1
    r["span_d"] = float(t[-1] - t[0])
    # collinearity of the trend regressor with the harmonic over this sampling
    doy = W * (t % 365.25)
    H = np.column_stack([np.ones_like(t), np.sin(doy), np.cos(doy)])
    cs, _, res = ols(H, t - t.mean())
    r["R2_trend_on_harm"] = float(1 - np.var(res) / np.var(t - t.mean()))
    # Kalman model (paper, amp_op definition) on raw vs detrended
    r["amp_kal"], r["peak_kal"] = kalman_amp(t, c)
    tr1 = tr1 if np.isfinite(tr1) else 0.0
    r["amp_kal_detr"], r["peak_kal_detr"] = kalman_amp(t, c * np.exp(-tr1 * (t - t.mean()) / 365.25))
    # synthetic: a pure trend of the unit's own rate and of +-5 %/yr through the Kalman sinusoid
    r["alias_own"], r["alias_own_peak"] = kalman_amp_noscreen(t, np.exp(tr1 * (t - t.mean()) / 365.25))
    r["alias_5pct"], r["alias_5pct_peak"] = kalman_amp_noscreen(t, np.exp(0.05 * (t - t.mean()) / 365.25))

    # ---- laser events (whole cache, all modules) and their C_L steps ----
    ev, duty = lifetime_events(cache, iv["laser_life_time"])
    r["duty"] = duty
    lc_all = np.log(np.where(cons["c"] > 0, cons["c"], np.nan))
    for e in ev:
        before = (cons["t"] < e["t"]) & (cons["t"] >= e["t"] - 45)
        after = (cons["t"] >= e["t"]) & (cons["t"] < e["t"] + 45)
        e["dlnC_pct"] = 100 * float(np.nanmedian(lc_all[after]) - np.nanmedian(lc_all[before])) \
            if np.isfinite(lc_all[before]).sum() >= 5 and np.isfinite(lc_all[after]).sum() >= 5 else np.nan
        e["in_main_era"] = bool(t[0] < e["t"] < t[-1] and e["mod_before"].endswith(main) and e["mod_after"].endswith(main))
        e["stream"], e["site"] = key, site
    # null distribution of the same 45-day step statistic at ordinary dates of the main era
    null = []
    for td in np.arange(t[0] + 45, t[-1] - 45, 10.0):
        b = (t < td) & (t >= td - 45)
        a = (t >= td) & (t < td + 45)
        if b.sum() >= 5 and a.sum() >= 5:
            null.append(100 * (np.median(lc[a]) - np.median(lc[b])))
    r["null_step_absmed"] = float(np.median(np.abs(null))) if null else np.nan
    for e in ev:
        e["null_step_absmed"] = r["null_step_absmed"]
    inera = [e for e in ev if e["in_main_era"] and e["kind"] in ("new_laser", "jump", "off")]
    r["n_events_main_era"] = len(inera)
    r["n_glitch"] = sum(e["kind"] in ("glitch", "encoding") for e in ev)
    r["n_modules"] = len(set(mods))
    # effect of in-era laser changes on the annual fit: harmonic with vs without per-segment offsets
    r["amp_harm_seg"] = np.nan
    r["kal_alias_steps"] = np.nan
    big = [e for e in inera if e["kind"] in ("new_laser", "jump")]
    if big:
        seg = np.searchsorted([e["t"] for e in big], t)
        doy = W * (t % 365.25)
        cols = [np.sin(doy), np.cos(doy)] + [(seg == k).astype(float) for k in range(len(big) + 1)]
        X = np.column_stack(cols)
        keepc = X[:, 2:].sum(axis=0) > 0
        X = np.column_stack([X[:, :2], X[:, 2:][:, keepc]])
        cf, _, _ = ols(X, lc)
        r["amp_harm_seg"] = 100 * float(np.hypot(cf[0], cf[1]))
        # apparent Kalman amplitude produced by the step pattern alone (segment offsets as signal)
        offs = X[:, 2:] @ cf[2:]
        r["kal_alias_steps"], _ = kalman_amp_noscreen(t, np.exp(offs - offs.mean()))

    # ---- diurnal ----
    r.update(diurnal(cache, lon, iv, main))
    r["diurnal_pred_P"] = -r.get("d_lp_noon", np.nan)
    r["diurnal_pred_P_jja"] = -r.get("d_lp_noon_jja", np.nan)
    series = dict(t=t, lc=lc, lp=lp, T=T, Tl=Tl, LT=LT)
    return r, ev, series


# ---------------------------------------------------------------------------------------------
# 11-year archive: multi-year ageing and an 18-month aliasing test on real data
# ---------------------------------------------------------------------------------------------
def archive_arm():
    rows, wins = [], []
    for f in sorted((OUT / "archive_11yr").glob("chm_*_1d.csv")):
        stn = f.stem.split("_")[1]
        t, c, T, ser, P = [], [], [], [], []
        for row in csv.DictReader(open(f, encoding="utf-8")):
            try:
                d = datetime.strptime(row["Time"], "%d-%b-%Y")
            except ValueError:
                continue
            fv = lambda k: float(row[k]) if row[k] not in ("", "NaN") else np.nan  # noqa: E731
            t.append((d - EPOCH).days + 0.5)
            c.append(fv("Lidar_constant"))
            T.append(fv("temp_int") - 273.15)
            ser.append(row["serlom"])
            P.append(fv("laser_pulses"))
        t, c, T, ser, P = map(np.array, (t, c, T, ser, P))
        good = (T > -23) & (c > 0) & (c < 1e14)
        for mod in np.unique(ser[good]):
            m = (ser == mod) & good
            if m.sum() < 200:
                continue
            tm, lc, Pm = t[m], np.log(c[m]), P[m]
            k = robust_mask(lc)
            tm, lc, Pm, cm = tm[k], lc[k], Pm[k], c[m][k]
            span = tm[-1] - tm[0]
            if span < 2 * 365:
                continue
            a, pk, tr, se, _ = harmonic(tm, lc, trend=True)
            lpm = np.log(Pm)
            okp = np.isfinite(lpm) & (np.abs(lpm - np.nanmedian(lpm)) < 0.1)
            ap, pkp, trp, _, _ = harmonic(tm[okp], -lpm[okp], trend=True) if okp.sum() > 100 else (np.nan,) * 5
            kal_full, kal_peak = kalman_amp(tm, cm)
            rows.append(dict(station=stn, module=mod, n=int(len(tm)), span_yr=span / 365.25,
                             trend_pct_yr=100 * tr, trend_se=100 * se, amp_harm=100 * a, peak_harm=pk,
                             amp_kal_full=kal_full, peak_kal_full=kal_peak,
                             amp_pred_P=100 * ap, peak_pred_P=pkp, trend_lnP_pct_yr=-100 * trp))
            # sliding 18-month windows (start every 91 d): Kalman sinusoid raw vs detrended (window trend
            # is NOT used: the multi-year module trend is the 'truth' we remove)
            for s in np.arange(tm[0], tm[-1] - 548, 91):
                w = (tm >= s) & (tm < s + 548)
                if w.sum() < 150:
                    continue
                ar, _ = kalman_amp(tm[w], cm[w])
                ad, _ = kalman_amp(tm[w], cm[w] * np.exp(-tr * (tm[w] - tm[w].mean()) / 365.25))
                wins.append(dict(station=stn, module=mod, start=float(s), n=int(w.sum()), amp_raw=ar,
                                 amp_detr=ad, amp_full=kal_full, trend_pct_yr=100 * tr))
    return rows, wins


# ---------------------------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------------------------
def write(rows, name):
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(Q15 / name, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def q(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return "n/a"
    return f"{np.median(x):+.2f} [{np.percentile(x, 25):+.2f}, {np.percentile(x, 75):+.2f}] (n {x.size})"


def main():
    Q15.mkdir(exist_ok=True)
    streams = chm15k_streams()
    rows, events, series = [], [], {}
    with ProcessPoolExecutor(N_WORKERS) as ex:
        for res in ex.map(analyse_stream, streams, chunksize=4):
            if res is None:
                continue
            r, ev, s = res
            rows.append(r)
            events.extend(ev)
            series[r["stream"]] = s
            print(f"  {r['stream']:28s} {r['site'][:16]:16s} n {r['n']:4d} dlnP/dT {r['dlnP_dT']:+.3f} %/K  "
                  f"ampP {r['amp_pred_P']:.2f} %  trend {r['trend_pct_yr']:+.1f} %/yr  ev {r['n_events_main_era']}",
                  flush=True)
    write(rows, "q15a_units.csv")
    write(events, "q15a_laser_events.csv")
    arows, wins = archive_arm()
    write(arows, "q15a_archive_ageing.csv")
    write(wins, "q15a_archive_windows.csv")
    np.savez_compressed(Q15 / "q15a_series.npz", **{f"{k}__{v}": a for k, s in series.items() for v, a in s.items()})

    net = [r for r in rows if r["stream"] != PAYERNE]
    A = {k: np.array([r.get(k, np.nan) for r in net], float) for k in net[0] if k not in ("stream", "site", "module")}
    print(f"\nunits {len(net)} (Payerne excluded)")
    for k in ("sd_lnC_anom", "sd_lnP_anom", "sd_lnP_anom_mad", "share_var_P", "dlnP_dT", "kappa_E_implied",
              "g_lnC_lnP", "r_lnC_lnP_anom", "kappa_ctrl", "sd_Tlom_anom", "sd_lnC_from_Tlom",
              "amp_pred_P", "peak_pred_P", "trend_lnP_pct_yr", "Tlom_amp_K", "amp_pred_Tlom", "lom_night_exc_frac",
              "amp_harm", "peak_harm", "amp_harm_trend", "trend_pct_yr", "trend_se", "R2_trend_on_harm",
              "amp_kal", "amp_kal_detr", "alias_own", "alias_5pct", "alias_5pct_peak", "duty",
              "n_events_main_era", "amp_harm_seg", "kal_alias_steps", "d_lp_noon", "d_lp_noon_jja",
              "d_tl_noon", "d_tl_noon_jja", "d_ti_noon_jja", "lom_day_exc_frac", "lom_noon_exc_jja",
              "dlnP_dTlom_exc", "dlnP_dTint_day", "diurnal_pred_P", "diurnal_pred_P_jja", "laser_age_kh"):
        print(f"  {k:22s} {q(A[k])}")
    sig = np.abs(A["trend_pct_yr"]) > 2 * A["trend_se"]
    print(f"  trend significant (2 se) in {sig.sum()} / {np.isfinite(A['trend_pct_yr']).sum()} units; "
          f"positive {np.sum(sig & (A['trend_pct_yr'] > 0))}, negative {np.sum(sig & (A['trend_pct_yr'] < 0))}")
    from scipy.stats import spearmanr
    ok = np.isfinite(A["laser_age_kh"]) & np.isfinite(A["amp_kal"])
    print(f"  Spearman amp_kal vs laser age: {spearmanr(A['laser_age_kh'][ok], A['amp_kal'][ok])}")
    ok = np.isfinite(A["trend_pct_yr"]) & np.isfinite(A["laser_age_kh"])
    print(f"  Spearman trend vs laser age: {spearmanr(A['laser_age_kh'][ok], A['trend_pct_yr'][ok])}")
    ok = np.isfinite(A["amp_pred_P"]) & np.isfinite(A["amp_kal"])
    print(f"  Spearman amp_kal vs amp_pred_P: {spearmanr(A['amp_pred_P'][ok], A['amp_kal'][ok])}")
    ev_main = [e for e in events if e["in_main_era"] and e["kind"] != "swap_only"]
    kinds = {}
    for e in events:
        kinds[e["kind"]] = kinds.get(e["kind"], 0) + 1
    print(f"  events all {kinds}; in main eras {len(ev_main)}")
    for kind in ("new_laser", "jump", "off", "glitch", "encoding"):
        st = [e["dlnC_pct"] for e in events if e["kind"] == kind and not e["swap"]]
        print(f"    C step at {kind:5s} (same module): {q(st)}")
    st = [e["dlnC_pct"] for e in events if e["swap"]]
    print(f"    C step at module swap: {q(st)}")
    if PAYERNE in [r["stream"] for r in rows]:
        p = [r for r in rows if r["stream"] == PAYERNE][0]
        print("  PAYERNE:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in p.items()})
    print("\narchive modules")
    for r in arows:
        print(f"  {r['station']} {r['module']} span {r['span_yr']:.1f} yr trend {r['trend_pct_yr']:+.2f}+/-{r['trend_se']:.2f}"
              f" %/yr  amp_harm {r['amp_harm']:.1f} peak {r['peak_harm']:.0f} amp_kal {r['amp_kal_full']:.1f}  "
              f"ampP {r['amp_pred_P']:.2f} peakP {r['peak_pred_P']:.0f} trendP {r['trend_lnP_pct_yr']:+.2f}")
    if wins:
        Wd = {k: np.array([w[k] for w in wins], float) for k in ("amp_raw", "amp_detr", "amp_full")}
        print(f"  18-month windows {len(wins)}: amp raw {q(Wd['amp_raw'])}, detrended {q(Wd['amp_detr'])}, "
              f"raw-detr {q(Wd['amp_raw'] - Wd['amp_detr'])}, raw-full {q(Wd['amp_raw'] - Wd['amp_full'])}")
    make_figures(rows, events, arows, wins, series)


# ---------------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------------
def make_figures(rows, events, arows, wins, series):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    net = [r for r in rows if r["stream"] != PAYERNE]
    A = {k: np.array([r.get(k, np.nan) for r in net], float) for k in net[0] if k not in ("stream", "site", "module")}
    fig, axs = plt.subplots(2, 3, figsize=(18, 9))
    ax = axs[0, 0]
    ax.hist(A["dlnP_dT"], bins=40, color="#1f77b4")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("d ln(pulses)/dT$_{int}$, 31-d anomalies [%/K]")
    ax.set_ylabel("units")
    ax.set_title(f"PRF vs temperature: median {np.nanmedian(A['dlnP_dT']):+.3f} %/K\n"
                 f"(energy if E~1/PRF: {np.nanmedian(A['kappa_E_implied']):+.3f} %/K vs κ ≈ +1.5 %/K)", fontsize=10)
    ax = axs[0, 1]
    ax.scatter(A["amp_pred_P"], A["amp_kal"], s=18, c=A["peak_pred_P"], cmap="twilight", vmin=0, vmax=365)
    ax.plot([0, 30], [0, 30], "k--", lw=0.8)
    ax.set_xlim(0, max(2, np.nanpercentile(A["amp_pred_P"], 99) * 1.2))
    ax.set_xlabel("annual amplitude of predicted effect −ln(pulses) [%]")
    ax.set_ylabel("C$_L$ annual amplitude (Kalman, amp_op) [%]")
    ax.set_title("pulse-count seasonal cycle is ~40× too small", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[0, 2]
    d = A["d_lp_noon_jja"]
    ax.hist(-d[np.isfinite(d)], bins=40, color="#ff7f0e", label="JJA")
    ax.hist(-A["d_lp_noon"][np.isfinite(A["d_lp_noon"])], bins=40, color="#1f77b4", alpha=0.6, label="all")
    ax.set_xlabel("predicted diurnal C$_L$ effect −Δln(pulses) noon−night [%]")
    ax.set_ylabel("units")
    ax.legend()
    ax.set_title("diurnal: PRF change at noon (E~1/PRF)", fontsize=10)
    ax = axs[1, 0]
    ax.errorbar(np.arange(len(arows)), [r["trend_pct_yr"] for r in arows], yerr=[2 * r["trend_se"] for r in arows],
                fmt="o", color="#2ca02c")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(np.arange(len(arows)))
    ax.set_xticklabels([f"{r['station']}\n{r['module']}\n{r['span_yr']:.0f}y" for r in arows], fontsize=7)
    ax.set_ylabel("ln C$_L$ trend, harmonic controlled [%/yr]")
    ax.set_title("ageing per optical module, 11-yr archive (±2 se)", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[1, 1]
    ax.scatter(np.abs(A["trend_pct_yr"]), A["alias_own"], s=18, color="#9467bd", label="unit's own trend")
    ax.scatter(np.full_like(A["alias_5pct"], 5.0), A["alias_5pct"], s=10, color="#8c564b", alpha=0.4,
               label="±5 %/yr synthetic")
    ax.set_xlabel("|linear trend| [%/yr]")
    ax.set_ylabel("apparent Kalman annual amplitude from trend alone [%]")
    ax.set_title("aliasing of a trend into the 18-month annual fit", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axs[1, 2]
    if wins:
        wr = np.array([w["amp_raw"] for w in wins])
        wd = np.array([w["amp_detr"] for w in wins])
        ax.scatter(wd, wr, s=12, color="#17becf")
        lim = [0, np.nanpercentile(np.r_[wr, wd], 99) * 1.05]
        ax.plot(lim, lim, "k--", lw=0.8)
        ax.set_xlabel("18-month Kalman amplitude, module trend removed [%]")
        ax.set_ylabel("18-month Kalman amplitude, raw [%]")
        ax.set_title(f"archive sliding windows (n {len(wins)}): raw − detrended median "
                     f"{np.nanmedian(wr - wd):+.1f} %", fontsize=10)
        ax.grid(alpha=0.3)
    fig.suptitle("q15a – CHM15k laser transmitter: pulse rate, ageing, aliasing", fontsize=12)
    fig.tight_layout()
    fig.savefig(Q15 / "q15a_laser_overview.png", dpi=120)
    plt.close(fig)

    # example units: ln C anomaly-free series, -ln P (x10 for visibility) and laser lifetime
    ex = [k for k in ("0-20000-0-12843_A", "0-20000-0-06240_A", "0-20008-0-INO_A", "0-100-20000-0000001_A")
          if k in series]
    fig, axs = plt.subplots(1, len(ex), figsize=(5 * len(ex), 4.5), squeeze=False)
    for ax, k in zip(axs[0], ex):
        s = series[k]
        x = [EPOCH + timedelta(days=float(v)) for v in s["t"]]
        ax.plot(x, 100 * (s["lc"] - np.median(s["lc"])), ".", ms=3, color="#1f77b4", label="ln C$_L$ [%]")
        ax.plot(x, 100 * 10 * (-(s["lp"] - np.median(s["lp"]))), ".", ms=3, color="#d62728",
                label="−ln pulses ×10 [%]")
        ax2 = ax.twinx()
        ax2.plot(x, s["LT"] / 1000, "-", color="0.5", lw=1)
        ax2.set_ylabel("laser_life_time [kh]", color="0.4")
        ax.set_title(k, fontsize=9)
        ax.set_ylabel("deviation from median [%]")
        ax.legend(fontsize=7, loc="upper left")
        ax.tick_params(axis="x", labelrotation=30, labelsize=7)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q15 / "q15a_examples.png", dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Post-hoc network statistics (reads the CSVs written by main)
# ---------------------------------------------------------------------------------------------
def post():
    import pandas as pd
    from scipy.stats import rankdata, spearmanr
    d = pd.read_csv(Q15 / "q15a_units.csv")
    pay = d[d.stream == PAYERNE]
    d = d[d.stream != PAYERNE].copy()
    gen = pd.read_csv(OUT / "q14_generation.csv")[["stream", "year"]]
    d = d.merge(gen, on="stream", how="left")
    lines = []
    out = lines.append
    # E ~ 1/PRF test: inverse-variance weighted network coefficient of ln C anomalies on ln P anomalies
    ok = np.isfinite(d.g_lnC_lnP) & (d.g_lnC_lnP_se > 0)
    w = 1 / d.g_lnC_lnP_se[ok] ** 2
    gm = float(np.sum(w * d.g_lnC_lnP[ok]) / np.sum(w))
    out(f"weighted network d lnC / d lnP (31-d anomalies, T_int controlled): {gm:+.2f} +/- {1 / np.sqrt(w.sum()):.2f}"
        f"  (E~1/PRF predicts -1, constant E predicts 0; n {ok.sum()})")
    ok = np.isfinite(d.dlnP_dT) & (d.dlnP_dT_se > 0)
    w = 1 / d.dlnP_dT_se[ok] ** 2
    out(f"weighted network d lnP/dT_int: {np.sum(w * d.dlnP_dT[ok]) / w.sum():+.3f} %/K; "
        f"negative in {np.mean(d.dlnP_dT[ok] < 0) * 100:.0f} % of units; "
        f"Spearman with kappa_ctrl {spearmanr(d.dlnP_dT[ok], d.kappa_ctrl[ok])[0]:+.2f}")

    def partial_spearman(x, y, z):
        m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        rx, ry, rz = (rankdata(v[m]) for v in (x, y, z))
        ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
        ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
        return float(np.corrcoef(ex, ey)[0, 1]), int(m.sum())
    x, y, z = d.laser_age_kh.values, d.amp_kal.values, d.year.values.astype(float)
    out(f"Spearman amp_kal vs laser age {spearmanr(x, y, nan_policy='omit')[0]:+.2f}; vs module year "
        f"{spearmanr(z, y, nan_policy='omit')[0]:+.2f}; partial (age | year) {partial_spearman(x, y, z)}")
    out(f"Spearman ln C trend vs ln P trend {spearmanr(d.trend_pct_yr, d.trend_lnP_pct_yr, nan_policy='omit')[0]:+.2f};"
        f" vs laser age {spearmanr(d.trend_pct_yr, d.laser_age_kh, nan_policy='omit')[0]:+.2f};"
        f" vs module year {spearmanr(d.trend_pct_yr, d.year, nan_policy='omit')[0]:+.2f}")
    out(f"ln C trend 2025-26: median {d.trend_pct_yr.median():+.2f} %/yr = {d.trend_pct_yr.median() * 1000 / 8766:+.2f}"
        f" %/1000 h (duty {d.duty.median():.3f}); weighted mean "
        f"{np.average(d.trend_pct_yr, weights=1 / d.trend_se ** 2):+.2f} %/yr")
    out(f"Kalman amplitude raw - detrended (own trend): {q(d.amp_kal - d.amp_kal_detr)}")
    out(f"laser counter wrapped (65536 h) inside 2025-26 in {int(d.wrap_in_record.sum())} units; laser age lower "
        f"bound {q(d.laser_age_kh)} kh, upper bound {q(d.laser_age_max_kh)} kh")
    for col in ("laser_age_kh", "laser_age_max_kh"):
        out(f"  {col}: Spearman vs amp_kal {spearmanr(d[col], d.amp_kal, nan_policy='omit')[0]:+.2f} "
            f"(partial | module year {partial_spearman(d[col].values, d.amp_kal.values, z)[0]:+.2f}); vs trend "
            f"{spearmanr(d[col], d.trend_pct_yr, nan_policy='omit')[0]:+.2f}; vs kappa_ctrl "
            f"{spearmanr(d[col], d.kappa_ctrl, nan_policy='omit')[0]:+.2f}")
    ev_all = pd.read_csv(Q15 / "q15a_laser_events.csv")
    ev_all = ev_all[ev_all.stream != PAYERNE]
    out(f"lifetime events (all eras): {ev_all.groupby('kind').size().to_dict()}")
    for kind in ("new_laser", "jump", "off", "glitch", "encoding"):
        e = ev_all[(ev_all.kind == kind) & ~ev_all.swap & np.isfinite(ev_all.dlnC_pct)]
        if len(e):
            out(f"  |C step| at {kind:9s}: median {e.dlnC_pct.abs().median():.1f} % vs null (same units) "
                f"{e.null_step_absmed.median():.1f} %  (n {len(e)})")
    e = ev_all[ev_all.swap & np.isfinite(ev_all.dlnC_pct)]
    out(f"  |C step| at module swaps: median {e.dlnC_pct.abs().median():.1f} %, signed {q(e.dlnC_pct)}")
    # what a module/laser swap step would do to the annual fit if the series were NOT split at the swap
    tt = np.arange(20089.0, 20647.0, 2.0)            # every 2nd night, 2025-01-01 .. 2026-07-12
    step = float(e.dlnC_pct.median()) / 100
    al = []
    for ts in np.arange(tt[0] + 30, tt[-1] - 30, 15.0):
        a_, pk_ = kalman_amp_noscreen(tt, np.exp(np.where(tt >= ts, step, 0.0)))
        al.append((ts, a_, pk_))
    al = np.array(al)
    out(f"  synthetic: a {100 * step:+.1f} % swap step NOT split, anywhere in 2025-01..2026-07 -> Kalman annual amp "
        f"median {np.median(al[:, 1]):.1f} %, max {al[:, 1].max():.1f} % (vs C_L amp ~16 %)")
    ev = d[np.isfinite(d.amp_harm_seg)]
    out(f"units with in-era laser events (reset/jump/off>14 d): {len(ev)}; harmonic amp without vs with segment "
        f"offsets: {q(ev.amp_harm - ev.amp_harm_seg)}; Kalman amp of the step pattern alone {q(ev.kal_alias_steps)}")
    for _, r in ev.sort_values("kal_alias_steps", ascending=False).head(6).iterrows():
        out(f"   {r.stream:26s} {r.site[:16]:16s} amp {r.amp_harm:5.1f} -> seg {r.amp_harm_seg:5.1f}  "
            f"step-alias {r.kal_alias_steps:4.1f} %")
    hot = d[d.lom_noon_exc_jja > 0]
    out(f"units whose LOM leaves regulation at summer noon (median JJA noon excess > 0): {len(hot)}: "
        + ", ".join(f"{r.site[:12]} {r.lom_noon_exc_jja:.1f} K" for _, r in hot.iterrows()))
    out(f"daytime LOM excursion fraction (blocks > set-point + 0.5 K, 08-18 solar): {q(d.lom_day_exc_frac * 100)} %")
    out(f"d lnP / d T_lom excess by day: {q(d.dlnP_dTlom_exc)} %/K")
    for k in ("0-20000-0-12843_A", "0-20000-0-06240_A", "0-20008-0-INO_A", "0-100-20000-0000001_A",
              "0-20000-0-03391_A"):
        if (d.stream == k).any():
            r = d[d.stream == k].iloc[0]
            out(f"  example {k} {r.site[:14]}: dlnP/dT {r.dlnP_dT:+.3f} %/K, ampP {r.amp_pred_P:.2f} % (DOY "
                f"{r.peak_pred_P:.0f}), amp_kal {r.amp_kal:.1f} %, trend {r.trend_pct_yr:+.1f}+/-{r.trend_se:.1f} %/yr, "
                f"alias {r.alias_own:.2f} %, diurnal -dlnP JJA {r.diurnal_pred_P_jja:+.2f} %, laser {r.laser_age_kh:.0f} kh")
    if len(pay):
        r = pay.iloc[0]
        out(f"  PAYERNE (excluded): dlnP/dT {r.dlnP_dT:+.3f} %/K, ampP {r.amp_pred_P:.2f} %, trend {r.trend_pct_yr:+.1f}"
            f" %/yr, diurnal -dlnP JJA {r.diurnal_pred_P_jja:+.2f} %, LOM JJA noon excess {r.lom_noon_exc_jja:.1f} K, "
            f"LOM day excursion {100 * r.lom_day_exc_frac:.0f} %")
    wn = pd.read_csv(Q15 / "q15a_archive_windows.csv")
    wn["diff"] = wn.amp_raw - wn.amp_detr
    g = wn.groupby(["station", "module"]).agg(trend=("trend_pct_yr", "first"), dmed=("diff", "median"),
                                              dmax=("diff", lambda v: np.nanmax(np.abs(v))), n=("diff", "size"))
    out("archive 18-month windows, |raw - detrended| Kalman amplitude per module:")
    for (st, mo), r in g.iterrows():
        out(f"   {st} {mo}: trend {r.trend:+6.2f} %/yr  median diff {r.dmed:+.2f}  max |diff| {r.dmax:.2f}  (n {r.n})")
    a = pd.read_csv(Q15 / "q15a_archive_ageing.csv")
    out(f"archive module trends: {q(a.trend_pct_yr)} %/yr; |trend| {q(a.trend_pct_yr.abs())}; ln P trend "
        f"{q(a.trend_lnP_pct_yr)}; -lnP annual amp {q(a.amp_pred_P)} peak {q(a.peak_pred_P)}")
    txt = "\n".join(lines)
    print(txt)
    (Q15 / "q15a_post.txt").write_text(txt, encoding="utf-8")


if __name__ == "__main__":
    if "--post" not in sys.argv:
        main()
    post()
