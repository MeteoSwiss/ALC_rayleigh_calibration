# -*- coding: utf-8 -*-
"""q15b: receiver gain chain and internal reference of the CHM15k vs the calibration-constant variability.

Question (q15 brief): the CHM15k Rayleigh constant C_L has an annual cycle of +-16 % (peak ~DOY 176) and a
night-to-night scatter of ~13 %, of which the temperature normalisation (kappa x T_int) removes only ~35 %.
This script quantifies the receiver-side factors:

  1. the internal reference pulse p_calc (`calibration_pulse`). The firmware divides the signal by it:
       beta_raw = ((P/lp) - b) / (cs * o(r) * p_calc) * r^2
     so C_L (calibrated on the divided signal) is proportional to G_signal / p_calc. If p_calc follows the
     real receiver gain, the division removes the gain change and ln C_L is independent of ln p_calc
     (slope 0). If p_calc moves for reasons the backscatter path does not share (reference-source ageing
     or temperature), the variation is imported into C_L with slope -1.
     Tests: (a) nightly anomalies (31-day running median removed) of ln C_L on ln p_calc with the T_int
     anomaly held fixed; (b) the same with co-located same-wavelength pairs (atmosphere cancels), night
     and day, with temperature and solar background held fixed; (c) annual harmonic and trend of the night
     p_calc, and the C_L cycle it would imply under slope -1 (bound) and under the measured slope.
  2. detector temperature regulation (`temperature_detector`): excursions from the set-point, combined with
     literature sensitivities of a Si APD at 1064 nm.
  3. window transmission (`window_transmission`): anomaly regression, annual harmonic, cleaning steps.
  4. `status_detector` (informative or not).

Outputs: C:/DATA/Projects/202610_CHM15k_stability/q15/q15b_*.csv, q15b_*.png, q15b_summary.md (written by hand
from the printed numbers).

Run:  .venv/Scripts/python.exe validation/chm15k_stability/q15b_receiver.py
"""
from __future__ import annotations

import csv
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (BLOCKS_PER_DAY, OUT, chm15k_streams, load_cache, load_constants,  # noqa: E402
                    solar_hour, stream_key)
from q1_temperature import ANOM_WINDOW, night_hk, robust_mask, running_median  # noqa: E402

warnings.filterwarnings("ignore")

Q15 = OUT / "q15"
PAYERNE = "0-20000-0-06610_A"          # abnormal unit: reported separately, excluded from medians
MIN_NIGHTS = 100                        # nights of the main optical module
MIN_SPAN = 330                          # days, for an annual harmonic
YEAR = 365.25
C_PEAK_NETWORK = 176.0                  # DOY of the network C_L maximum (q13)
STEP_PTS = 5.0                          # %-points, window-transmission step (cleaning / sudden soiling)
STEP_WIN = 3                            # days on each side for the step detection
STEP_C_WIN = 15.0                       # days on each side for the C_L step


# ---------------------------------------------------------------------------------------------
# Small statistics helpers
# ---------------------------------------------------------------------------------------------
def ols(X, y):
    """OLS coefficients, standard errors and R^2."""
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    dof = max(len(y) - X.shape[1], 1)
    s2 = res @ res / dof
    try:
        se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    except np.linalg.LinAlgError:
        se = np.full(X.shape[1], np.nan)
    r2 = 1 - (res @ res) / np.sum((y - y.mean()) ** 2)
    return coef, se, r2


def harmonic(t, y, trend=False):
    """Annual harmonic of y (natural units of y). Returns amplitude, DOY of maximum, trend per year.

    y = a + A sin(w d) + B cos(w d) [+ b yrs]  ->  amplitude sqrt(A^2+B^2), maximum at atan2(A, B)/w.
    The Kalman model of the brief has no trend, so the comparable amplitude is the trend-free one.
    """
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if len(t) < 60 or (t.max() - t.min()) < MIN_SPAN:
        return np.nan, np.nan, np.nan
    w = 2 * np.pi / YEAR
    cols = [np.ones_like(t), np.sin(w * t), np.cos(w * t)]
    if trend:
        cols.append((t - t.mean()) / YEAR)
    coef, _, _ = ols(np.column_stack(cols), y)
    # t is days since 1970: day-of-year phase of the epoch is 0 (1 Jan) up to leap-year drift (<1 d).
    amp = float(np.hypot(coef[1], coef[2]))
    peak = float((np.arctan2(coef[1], coef[2]) / w) % YEAR) + 1.0
    return amp, peak, float(coef[3]) if trend else np.nan


def project(amp, peak, ref_peak):
    """Signed component of a harmonic (amp, peak) along a reference harmonic peaking at ref_peak."""
    return amp * np.cos(2 * np.pi * (peak - ref_peak) / YEAR)


def anomaly(t, x):
    return x - running_median(t, x, ANOM_WINDOW)


# ---------------------------------------------------------------------------------------------
# Per-unit analysis
# ---------------------------------------------------------------------------------------------
def block_times(cache):
    return (cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY)


def diurnal_hk(cache, lon, names):
    """Noon (11-13 solar h) minus night (22-04 solar h) per calendar day, median over days. All / JJA."""
    hv = list(cache["hk_vars"])
    sh = solar_hour(block_times(cache), lon)
    noon = (sh >= 11) & (sh < 13)
    night = (sh >= 22) | (sh < 4)
    from datetime import datetime, timedelta
    mon = np.array([(datetime(1970, 1, 1) + timedelta(days=float(d))).month for d in cache["day"]])
    jja = np.isin(mon, (6, 7, 8))
    out = {}
    for name, log in names:
        x = cache["hk"][..., hv.index(name)].astype(float)
        if log:
            x = np.where(x > 0, np.log(x), np.nan)
        with np.errstate(invalid="ignore"):
            dn = np.nanmedian(np.where(noon, x, np.nan), axis=1) - np.nanmedian(np.where(night, x, np.nan), axis=1)
        scale = 100.0 if log else 1.0
        out[f"{name}_noon"] = scale * np.nanmedian(dn)
        out[f"{name}_noon_JJA"] = scale * np.nanmedian(dn[jja]) if np.isfinite(dn[jja]).sum() > 10 else np.nan
    return out


def detector_regulation(cache, lon):
    """Detector temperature vs its set-point (mode over the record). Fractions of 10-min blocks."""
    hv = list(cache["hk_vars"])
    td = cache["hk"][..., hv.index("temperature_detector")].astype(float)
    ti = cache["hk"][..., hv.index("temp_int")].astype(float)
    ok = np.isfinite(td)
    if ok.sum() < 1000:
        return dict(det_setpoint=np.nan)
    vals, cnt = np.unique(np.round(td[ok], 2), return_counts=True)
    sp = float(vals[np.argmax(cnt)])
    dev = td - sp
    sh = solar_hour(block_times(cache), lon)
    night = (sh >= 22) | (sh < 4)
    from datetime import datetime, timedelta
    mon = np.array([(datetime(1970, 1, 1) + timedelta(days=float(d))).month for d in cache["day"]])
    jja_day = np.isin(mon, (6, 7, 8))[:, None] & (sh >= 10) & (sh < 18)
    exc = np.abs(dev) >= 0.5
    sd = cache["hk"][..., hv.index("status_detector")].astype(float)
    return dict(det_setpoint=sp, det_n_values=int(len(vals)),
                det_exc_frac=float(np.mean(exc[ok])),
                det_exc_frac_night=float(np.mean(exc[ok & night])) if (ok & night).any() else np.nan,
                det_exc_frac_jja_day=float(np.mean(exc[ok & jja_day])) if (ok & jja_day).any() else np.nan,
                det_max_dev=float(np.nanmax(np.abs(dev))),
                det_tint_p99=float(np.nanpercentile(ti, 99)),
                status_det_frac_not100=float(np.mean(sd[np.isfinite(sd)] != 100)) if np.isfinite(sd).any() else np.nan)


def window_steps(cache, lon, t_c, lnc_corr):
    """Sudden window-transmission steps (night daily medians) and the C_L step across them."""
    hv = list(cache["hk_vars"])
    w = cache["hk"][..., hv.index("window_transmission")].astype(float)
    sh = solar_hour(block_times(cache), lon)
    night = (sh >= 22) | (sh < 4)
    wd = np.nanmedian(np.where(night, w, np.nan), axis=1)
    day = cache["day"]
    events = []
    n = len(day)
    for i in range(STEP_WIN, n - STEP_WIN):
        # contiguous days only (no gap in the cache around the candidate step)
        if day[i + STEP_WIN - 1] - day[i - STEP_WIN] > 2 * STEP_WIN + 1:
            continue
        before = wd[i - STEP_WIN:i]
        after = wd[i:i + STEP_WIN]
        if np.isfinite(before).sum() < 2 or np.isfinite(after).sum() < 2:
            continue
        events.append((i, np.nanmedian(after) - np.nanmedian(before)))
    out = []
    if not events:
        return out
    idx = np.array([e[0] for e in events])
    # placebo: the same C_L step statistic on every 7th day with no window change (|dw| < 1 pt)
    for k in range(0, len(events), 7):
        if abs(events[k][1]) < 1.0:
            d0 = day[events[k][0]]
            b = (t_c >= d0 - STEP_C_WIN) & (t_c < d0) & np.isfinite(lnc_corr)
            a = (t_c >= d0 + 0.5) & (t_c < d0 + STEP_C_WIN) & np.isfinite(lnc_corr)
            if b.sum() >= 3 and a.sum() >= 3:
                out.append(dict(day=float(d0), w_before=np.nan, w_after=np.nan, dlnw=0.0,
                                dlnC=100 * (np.median(lnc_corr[a]) - np.median(lnc_corr[b])), nb=int(b.sum()),
                                na=int(a.sum()), placebo=1))
    dw = np.array([e[1] for e in events])
    # keep local extrema of |dw| above threshold, at least 10 days apart
    cand = np.where(np.abs(dw) >= STEP_PTS)[0]
    taken = []
    for k in cand[np.argsort(-np.abs(dw[cand]))]:
        if all(abs(day[idx[k]] - day[idx[j]]) > 10 for j in taken):
            taken.append(k)
    for k in taken:
        d0 = day[idx[k]]
        wb = np.nanmedian(wd[max(idx[k] - STEP_WIN, 0):idx[k]])
        wa = np.nanmedian(wd[idx[k]:idx[k] + STEP_WIN])
        b = (t_c >= d0 - STEP_C_WIN) & (t_c < d0) & np.isfinite(lnc_corr)
        a = (t_c >= d0 + 0.5) & (t_c < d0 + STEP_C_WIN) & np.isfinite(lnc_corr)
        if b.sum() >= 3 and a.sum() >= 3 and wb > 0 and wa > 0:
            out.append(dict(day=float(d0), w_before=float(wb), w_after=float(wa), dlnw=100 * np.log(wa / wb),
                            dlnC=100 * (np.median(lnc_corr[a]) - np.median(lnc_corr[b])), nb=int(b.sum()),
                            na=int(a.sum()), placebo=0))
    return out


def analyse_unit(args):
    wmo, ident, site, lat, lon = args
    key = stream_key(wmo, ident)
    cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
    if cons is None or cache is None:
        return None
    hv = list(cache["hk_vars"])
    mod = max(set(cons["module"]), key=list(cons["module"]).count)
    hk = night_hk(cache, cons["t0"], cons["t1"])
    T = hk[:, hv.index("temp_int")]
    pc = hk[:, hv.index("calibration_pulse")]
    w = hk[:, hv.index("window_transmission")]
    m = (cons["module"] == mod) & np.isfinite(T) & (cons["c"] > 0) & (pc > 0) & (w > 0)
    if m.sum() < MIN_NIGHTS:
        return None
    t, c, T, pc, w = cons["t"][m], cons["c"][m], T[m], pc[m], w[m]
    lnc = np.log(c)
    keep = robust_mask(lnc)                         # one screen per optical module (main module only here)
    t, lnc, T, lnp, lnw = t[keep], lnc[keep], T[keep], np.log(pc[keep]), np.log(w[keep])
    o = np.argsort(t)
    t, lnc, T, lnp, lnw = t[o], lnc[o], T[o], lnp[o], lnw[o]
    la, ta, pa, wa = anomaly(t, lnc), anomaly(t, T), anomaly(t, lnp), anomaly(t, lnw)
    ok = np.isfinite(la) & np.isfinite(ta) & np.isfinite(pa) & np.isfinite(wa)
    row = dict(stream=key, site=site, module=mod, lat=lat, lon=lon, n=int(ok.sum()),
               span_d=float(t[-1] - t[0]))
    one = np.ones(ok.sum())
    # --- 1a: p_calc vs T_int (anomalies) and ln C vs ln p_calc
    cp, sp, _ = ols(np.column_stack([one, ta[ok]]), pa[ok])
    row.update(kappa_pcal=100 * cp[1], se_kappa_pcal=100 * sp[1])
    c0, s0, r2_t = ols(np.column_stack([one, ta[ok]]), la[ok])
    c1, s1, _ = ols(np.column_stack([one, pa[ok]]), la[ok])
    c2, s2, r2_tp = ols(np.column_stack([one, ta[ok], pa[ok]]), la[ok])
    c3, s3, r2_tpw = ols(np.column_stack([one, ta[ok], pa[ok], wa[ok]]), la[ok])
    c4, s4, r2_tw = ols(np.column_stack([one, ta[ok], wa[ok]]), la[ok])
    row.update(kappa_T=100 * c0[1], se_kappa_T=100 * s0[1],
               g_pcal_simple=c1[1], se_g_pcal_simple=s1[1],
               g_pcal=c2[2], se_g_pcal=s2[2], kappa_T_ctrl=100 * c2[1],
               g_win=c3[3], se_g_win=s3[3], g_pcal_full=c3[2],
               sd_la=100 * np.std(la[ok]), sd_pa=100 * np.std(pa[ok]), sd_wa=100 * np.std(wa[ok]),
               sd_ta=float(np.std(ta[ok])),
               dr2_pcal=r2_tp - r2_t, dr2_win=r2_tw - r2_t, r2_T=r2_t,
               corr_la_pa=float(np.corrcoef(la[ok], pa[ok])[0, 1]),
               corr_ta_pa=float(np.corrcoef(ta[ok], pa[ok])[0, 1]),
               corr_la_wa=float(np.corrcoef(la[ok], wa[ok])[0, 1]))
    # night-to-night contributions in % std: artefact bound (|slope| 1) and measured slope
    row.update(n2n_pcal_bound=100 * np.std(pa[ok]), n2n_pcal_meas=100 * abs(c2[2]) * np.std(pa[ok]),
               n2n_win_meas=100 * abs(c3[3]) * np.std(wa[ok]), n2n_win_2way=200 * np.std(wa[ok]))
    # --- 1c / 3: annual harmonics (trend-free, as the Kalman model) and trends
    a_c, p_c, _ = harmonic(t, lnc)
    _, _, tr_c = harmonic(t, lnc, trend=True)
    a_p, p_p, _ = harmonic(t, lnp)
    _, _, tr_p = harmonic(t, lnp, trend=True)
    a_w, p_w, _ = harmonic(t, lnw)
    _, _, tr_w = harmonic(t, lnw, trend=True)
    a_T, p_T, _ = harmonic(t, T)
    row.update(amp_C=100 * a_c, peak_C=p_c, trend_C=100 * tr_c,
               amp_pcal=100 * a_p, peak_pcal=p_p, trend_pcal=100 * tr_p,
               amp_win=100 * a_w, peak_win=p_w, trend_win=100 * tr_w, amp_T=a_T, peak_T=p_T)
    # implied C_L harmonic. Bound (slope -1): amplitude a_p, maximum at the p_calc minimum.
    if np.isfinite(a_p):
        pk_bound = (p_p + YEAR / 2) % YEAR
        row.update(pcal_C_amp_bound=100 * a_p, pcal_C_peak_bound=pk_bound,
                   pcal_C_proj_bound=project(100 * a_p, pk_bound, C_PEAK_NETWORK),
                   pcal_C_proj_bound_own=project(100 * a_p, pk_bound, p_c) if np.isfinite(p_c) else np.nan)
        g = c2[2]
        pk_meas = p_p if g > 0 else pk_bound
        row.update(pcal_C_amp_meas=100 * abs(g) * a_p, pcal_C_peak_meas=pk_meas,
                   pcal_C_proj_meas=project(100 * abs(g) * a_p, pk_meas, C_PEAK_NETWORK))
    if np.isfinite(a_w):
        g = c3[3]
        pk_meas = p_w if g > 0 else (p_w + YEAR / 2) % YEAR
        row.update(win_C_amp_2way=200 * a_w, win_C_peak_2way=p_w,
                   win_C_proj_2way=project(200 * a_w, p_w, C_PEAK_NETWORK),
                   win_C_amp_meas=100 * abs(g) * a_w, win_C_peak_meas=pk_meas,
                   win_C_proj_meas=project(100 * abs(g) * a_w, pk_meas, C_PEAK_NETWORK))
    # median night levels
    row.update(pcal_night_med=float(np.exp(np.median(lnp))), win_night_med=float(np.exp(np.median(lnw))))
    # --- diurnal HK, detector regulation
    row.update(diurnal_hk(cache, lon, (("calibration_pulse", True), ("window_transmission", True),
                                       ("temperature_detector", False), ("bckgrd_rcs_0", True))))
    row.update(detector_regulation(cache, lon))
    # --- window steps, with C_L temperature-normalised by the unit's own anomaly kappa
    lnc_corr = lnc - c0[1] * (T - np.median(T))
    steps = window_steps(cache, lon, t, lnc_corr)
    for s in steps:
        s.update(stream=key, site=site)
    # series for the figure (Schiphol / Budapest / examples)
    series = dict(t=t, lnc=lnc, lnp=lnp, lnw=lnw, T=T, la=la, pa=pa, ta=ta, wa=wa)
    return row, steps, series


# ---------------------------------------------------------------------------------------------
# Pairs: atmosphere cancels; night and day regression of the ratio on the p_calc ratio
# ---------------------------------------------------------------------------------------------
def pair_tests():
    from q2b_pairs import is_night, relative          # reuse the night-referenced anomaly definition
    z = np.load(OUT / "q2b_frames.npz", allow_pickle=False)
    labels = sorted({k.split("|")[0] for k in z.files})
    rows = []
    band = "1.5-2.5 km"
    for lab in labels:
        fr = {k.split("|")[1]: z[k] for k in z.files if k.split("|")[0] == lab}
        rel = fr["rel_" + band]
        dlnp = fr["rel_lnp"]
        dT = fr["temp_int_A"] - fr["temp_int_B"]
        fr["dT"] = dT
        dT_rel = relative(fr, "dT")
        with np.errstate(invalid="ignore", divide="ignore"):
            fr["lnw"] = np.log(fr["window_transmission_A"] / fr["window_transmission_B"])
            bgA = fr["bckgrd_rcs_0_A"] / np.nanmedian(fr["bckgrd_rcs_0_A"])
            bgB = fr["bckgrd_rcs_0_B"] / np.nanmedian(fr["bckgrd_rcs_0_B"])
        dlnw = relative(fr, "lnw")
        dbg = bgA - bgB
        night = is_night(fr["sh"])
        day = (fr["sh"] >= 10) & (fr["sh"] < 16)
        for name, sel in (("night", night), ("day", day)):
            for ctrl in ("none", "T", "T+bg"):
                cols = [dlnp]
                if ctrl in ("T", "T+bg"):
                    cols.append(dT_rel)
                if ctrl == "T+bg" and name == "day":
                    cols.append(dbg)
                ok = sel & np.isfinite(rel) & np.all([np.isfinite(x) for x in cols], axis=0)
                if ok.sum() < 50:
                    continue
                X = np.column_stack([np.ones(ok.sum())] + [x[ok] for x in cols])
                cf, se, r2 = ols(X, rel[ok])
                rows.append(dict(pair=lab, hours=name, control=ctrl, n=int(ok.sum()), g_pcal=cf[1],
                                 se_g_pcal=se[1], sd_dlnp=100 * np.std(dlnp[ok]), sd_rel=100 * np.std(rel[ok]),
                                 kappa_dT=100 * cf[2] if len(cf) > 2 else np.nan))
            okw = sel & np.isfinite(rel) & np.isfinite(dlnw) & np.isfinite(dT_rel) & np.isfinite(dlnp)
            if okw.sum() >= 50 and np.std(dlnw[okw]) > 0:
                X = np.column_stack([np.ones(okw.sum()), dlnw[okw], dT_rel[okw], dlnp[okw]])
                cf, se, _ = ols(X, rel[okw])
                rows.append(dict(pair=lab, hours=name, control="window|T,pcal", n=int(okw.sum()), g_pcal=cf[1],
                                 se_g_pcal=se[1], sd_dlnp=100 * np.std(dlnw[okw]), sd_rel=100 * np.std(rel[okw]),
                                 kappa_dT=100 * cf[2]))
    return rows


def pair_unit_decomposition():
    """Each unit's own coefficient on its own p_calc within a pair (T_int anomalies of both held fixed).

    rel_AB = a_A dlnp_A - a_B dlnp_B + k_A dT_A - k_B dT_B; a = 0 when p_calc follows the real gain,
    a = -1 when the p_calc variation is not shared by the backscatter path.
    """
    from q2b_pairs import is_night, relative
    z = np.load(OUT / "q2b_frames.npz", allow_pickle=False)
    rows = []
    for lab in sorted({k.split("|")[0] for k in z.files}):
        fr = {k.split("|")[1]: z[k] for k in z.files if k.split("|")[0] == lab}
        with np.errstate(invalid="ignore", divide="ignore"):
            fr["lpA"] = np.log(fr["calibration_pulse_A"])
            fr["lpB"] = np.log(fr["calibration_pulse_B"])
        fr["TA"] = fr["temp_int_A"].astype(float)
        fr["TB"] = fr["temp_int_B"].astype(float)
        v = {k: relative(fr, k) for k in ("lpA", "lpB", "TA", "TB")}
        rel = fr["rel_1.5-2.5 km"]
        for name, sel in (("night", is_night(fr["sh"])), ("day", (fr["sh"] >= 10) & (fr["sh"] < 16))):
            X = [v["lpA"], v["lpB"], v["TA"], v["TB"]]
            ok = sel & np.isfinite(rel) & np.all([np.isfinite(x) for x in X], axis=0)
            if ok.sum() < 50:
                continue
            cf, se, _ = ols(np.column_stack([np.ones(ok.sum())] + [x[ok] for x in X]), rel[ok])
            rows.append(dict(pair=lab, hours=name, n=int(ok.sum()), a_A=cf[1], se_a_A=se[1], a_B=-cf[2],
                             se_a_B=se[2], kappa_A=100 * cf[3], kappa_B=-100 * cf[4],
                             sd_lpA=100 * np.std(v["lpA"][ok]), sd_lpB=100 * np.std(v["lpB"][ok])))
    return rows


# ---------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------
def q(x, p):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return np.percentile(x, p) if x.size else np.nan


def summarise(rows, key):
    x = np.array([r.get(key, np.nan) for r in rows], float)
    return f"{q(x, 50):+.3g} [IQR {q(x, 25):+.3g} .. {q(x, 75):+.3g}] (n {np.isfinite(x).sum()})"


def figures(rows, series, steps, pair_rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    net = [r for r in rows if r["stream"] != PAYERNE]
    A = {k: np.array([r.get(k, np.nan) for r in net], float)
         for k in ("g_pcal", "se_g_pcal", "kappa_pcal", "amp_pcal", "peak_pcal", "amp_C", "peak_C", "g_win",
                   "amp_win", "peak_win", "trend_pcal", "calibration_pulse_noon", "calibration_pulse_noon_JJA")}
    fig, axs = plt.subplots(2, 3, figsize=(18, 9))
    ax = axs[0, 0]
    g = A["g_pcal"]
    ax.hist(g[np.isfinite(g)], bins=np.linspace(-2, 2, 41), color="#1f77b4", alpha=0.8)
    ax.axvline(0, color="g", lw=1.5, label="0: p_calc = real gain (division correct)")
    ax.axvline(-1, color="r", lw=1.5, label="-1: p_calc artefact imported")
    ax.axvline(np.nanmedian(g), color="k", ls="--", label=f"median {np.nanmedian(g):+.2f}")
    ax.set_xlabel("slope d ln C_L / d ln p_calc (night anomalies, T_int held fixed)")
    ax.set_ylabel("units")
    ax.legend(fontsize=8)
    ax.set_title("1a  nightly C_L anomalies vs p_calc anomalies", fontsize=10)
    ax = axs[0, 1]
    pr = [r for r in pair_rows if (r["hours"] == "night" and r["control"] == "T") or (r["hours"] == "day" and r["control"] == "T+bg")]
    labs = [f"{r['pair'].replace('AMSTERDAM_AP_SCHIP', 'Schiphol')} {r['hours']}" for r in pr]
    yy = np.arange(len(pr))
    ax.errorbar([r["g_pcal"] for r in pr], yy, xerr=[2 * r["se_g_pcal"] for r in pr], fmt="o",
                color="#444")
    ax.axvline(0, color="g")
    ax.axvline(-1, color="r")
    ax.set_yticks(yy)
    ax.set_yticklabels(labs, fontsize=7)
    ax.set_xlabel("slope d ln(R_A/R_B) / d ln(p_A/p_B)  (+-2 se)")
    ax.set_title("1b  co-located pairs (atmosphere cancels)", fontsize=10)
    ax = axs[0, 2]
    okp = np.isfinite(A["amp_pcal"])
    sc = ax.scatter(A["peak_pcal"][okp], A["amp_pcal"][okp], c=A["trend_pcal"][okp], cmap="coolwarm",
                    vmin=-5, vmax=5, s=22)
    ax.axvline(C_PEAK_NETWORK, color="k", ls=":", label="C_L peak (DOY 176)")
    ax.axvline((C_PEAK_NETWORK + YEAR / 2) % YEAR, color="r", ls=":", label="p_calc minimum needed to\nexplain C_L (slope -1)")
    ax.set_xlabel("DOY of night p_calc maximum")
    ax.set_ylabel("annual amplitude of night ln p_calc [%]")
    plt.colorbar(sc, ax=ax, label="p_calc trend [%/yr]")
    ax.legend(fontsize=7)
    ax.set_title("1c  seasonal cycle of the night reference pulse", fontsize=10)
    ax = axs[1, 0]
    tp = np.array([r.get("trend_pcal", np.nan) for r in net])
    tc = np.array([r.get("trend_C", np.nan) for r in net])
    ax.plot(tp, tc, "o", ms=4, alpha=0.6)
    for r in net:
        if abs(r.get("trend_pcal", 0) or 0) > 10:
            ax.annotate(r["site"][:10].title(), (r["trend_pcal"], r["trend_C"]), fontsize=7)
    xx = np.array([-16, 16])
    ax.plot(xx, -xx + np.nanmedian(tc), "r--", lw=1, label="slope -1 (drift imported)")
    ax.axhline(np.nanmedian(tc), color="g", ls="--", lw=1, label="slope 0")
    ax.set_xlabel("trend of night ln p_calc [%/yr]")
    ax.set_ylabel("trend of ln C_L [%/yr]")
    ax.legend(fontsize=8)
    ax.set_title("1d  multi-month drift: C_L vs reference pulse", fontsize=10)
    ax = axs[1, 1]
    st = [s for s in steps if s["stream"] != PAYERNE and not s.get("placebo")]
    if st:
        x = np.array([s["dlnw"] for s in st])
        y = np.array([s["dlnC"] for s in st])
        ax.plot(x, y, "o", ms=4, alpha=0.6)
        xx = np.linspace(-40, 40, 3)
        ax.plot(xx, 2 * xx, "r--", lw=1, label="two-way physical (slope 2)")
        ax.plot(xx, 0 * xx, "g--", lw=1, label="no effect")
        ax.set_xlim(-40, 40)
        ax.set_ylim(-60, 60)
        ax.legend(fontsize=8)
    ax.set_xlabel("window-transmission step Δln w [%]")
    ax.set_ylabel("C_L step (T-normalised, ±15 d) [%]")
    ax.set_title("3  sudden window-transmission steps", fontsize=10)
    ax = axs[1, 2]
    d = A["calibration_pulse_noon"]
    dj = A["calibration_pulse_noon_JJA"]
    ax.hist(d[np.isfinite(d)], bins=30, alpha=0.6, label=f"all, median {np.nanmedian(d):.1f} %")
    ax.hist(dj[np.isfinite(dj)], bins=30, alpha=0.6, label=f"JJA, median {np.nanmedian(dj):.1f} %")
    ax.set_xlabel("noon minus night ln p_calc [%]")
    ax.set_ylabel("units")
    ax.legend(fontsize=8)
    ax.set_title("daytime drop of the reference pulse", fontsize=10)
    for a in axs.ravel():
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q15 / "q15b_receiver.png", dpi=120)
    plt.close(fig)


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main():
    Q15.mkdir(exist_ok=True)
    streams = chm15k_streams()
    rows, steps, series = [], [], {}
    with ProcessPoolExecutor(max_workers=4) as ex:
        for res in ex.map(analyse_unit, streams, chunksize=4):
            if res is None:
                continue
            r, s, ser = res
            rows.append(r)
            steps.extend(s)
            if r["stream"] in ("0-20000-0-06240_A", "0-20000-0-12843_A", PAYERNE):
                series[r["stream"]] = ser
    rows.sort(key=lambda r: r["stream"])
    write_csv(Q15 / "q15b_units.csv", rows)
    write_csv(Q15 / "q15b_window_steps.csv", steps)
    pair_rows = pair_tests()
    write_csv(Q15 / "q15b_pairs.csv", pair_rows)
    dec = pair_unit_decomposition()
    write_csv(Q15 / "q15b_pairs_units.csv", dec)
    for r in dec:
        print(f"  unit-decomp {r['pair']:24s} {r['hours']:5s} a_A {r['a_A']:+.2f}+-{r['se_a_A']:.2f} "
              f"a_B {r['a_B']:+.2f}+-{r['se_a_B']:.2f} (sd lnp {r['sd_lpA']:.1f}/{r['sd_lpB']:.1f} %)")
    net = [r for r in rows if r["stream"] != PAYERNE]
    print(f"units analysed {len(rows)} (network {len(net)} without Payerne)")
    for k in ("kappa_T", "kappa_pcal", "g_pcal_simple", "g_pcal", "kappa_T_ctrl", "g_win", "g_pcal_full",
              "corr_la_pa", "corr_ta_pa", "corr_la_wa", "dr2_pcal", "dr2_win",
              "sd_la", "sd_pa", "sd_wa", "n2n_pcal_bound", "n2n_pcal_meas", "n2n_win_meas", "n2n_win_2way",
              "amp_C", "peak_C", "trend_C", "amp_pcal", "peak_pcal", "trend_pcal", "pcal_C_proj_bound",
              "pcal_C_proj_bound_own", "pcal_C_amp_meas", "pcal_C_proj_meas",
              "amp_win", "peak_win", "trend_win", "win_C_amp_2way", "win_C_proj_2way", "win_C_amp_meas",
              "win_C_proj_meas", "pcal_night_med", "win_night_med",
              "calibration_pulse_noon", "calibration_pulse_noon_JJA", "window_transmission_noon",
              "window_transmission_noon_JJA", "temperature_detector_noon", "temperature_detector_noon_JJA",
              "det_n_values", "det_exc_frac", "det_exc_frac_night", "det_exc_frac_jja_day", "det_max_dev",
              "status_det_frac_not100"):
        print(f"  {k:28s} {summarise(net, k)}")
    # significance count of the p_calc slope
    g = np.array([r["g_pcal"] for r in net])
    se = np.array([r["se_g_pcal"] for r in net])
    print(f"  g_pcal: |g/se|>2 in {np.mean(np.abs(g / se) > 2) * 100:.0f} %; consistent with -1 (|g+1|<2se) "
          f"in {np.mean(np.abs(g + 1) < 2 * se) * 100:.0f} %; with 0 in {np.mean(np.abs(g) < 2 * se) * 100:.0f} %")
    w_ = 1 / se ** 2
    okw = np.isfinite(g) & np.isfinite(w_)
    print(f"  g_pcal inverse-variance weighted mean {np.sum(g[okw] * w_[okw]) / np.sum(w_[okw]):+.3f} "
          f"+- {1 / np.sqrt(np.sum(w_[okw])):.3f}")
    gw = np.array([r["g_win"] for r in net])
    sw = np.array([r["se_g_win"] for r in net])
    ww = 1 / sw ** 2
    okw = np.isfinite(gw) & np.isfinite(ww)
    print(f"  g_win inverse-variance weighted mean {np.sum(gw[okw] * ww[okw]) / np.sum(ww[okw]):+.3f} "
          f"+- {1 / np.sqrt(np.sum(ww[okw])):.3f}; |g/se|>2 in {np.mean(np.abs(gw / sw) > 2) * 100:.0f} %")
    sp = [r for r in rows if r["stream"] == PAYERNE]
    if sp:
        r = sp[0]
        print("  PAYERNE: " + ", ".join(f"{k} {r.get(k, np.nan):.3g}" for k in (
            "kappa_pcal", "g_pcal", "se_g_pcal", "g_win", "amp_pcal", "peak_pcal", "trend_pcal", "amp_win",
            "peak_win", "calibration_pulse_noon", "det_exc_frac", "det_exc_frac_jja_day")))
    # long-term: cross-unit trend of C_L vs trend of night p_calc (expected slope -1 if drifts are imported)
    from scipy.stats import spearmanr, theilslopes
    tp = np.array([r.get("trend_pcal", np.nan) for r in net])
    tc = np.array([r.get("trend_C", np.nan) for r in net])
    ok = np.isfinite(tp) & np.isfinite(tc)
    ts = theilslopes(tc[ok], tp[ok])
    print(f"  trend_C vs trend_pcal across units: Theil-Sen {ts.slope:+.2f} [{ts.low_slope:+.2f}, {ts.high_slope:+.2f}],"
          f" Spearman {spearmanr(tp[ok], tc[ok])[0]:+.2f} (p {spearmanr(tp[ok], tc[ok])[1]:.2f}); units with |trend_pcal|>10 %/yr:")
    for r in net:
        if abs(r.get("trend_pcal", 0) or 0) > 10:
            print(f"     {r['stream']:24s} {r['site'][:14]:14s} trend p_calc {r['trend_pcal']:+.1f} C_L {r['trend_C']:+.1f} %/yr")
    # circular (amplitude-weighted) mean phases
    for k_a, k_p in (("amp_pcal", "peak_pcal"), ("amp_win", "peak_win"), ("amp_C", "peak_C")):
        a = np.array([r.get(k_a, np.nan) for r in net])
        ph = np.array([r.get(k_p, np.nan) for r in net]) * 2 * np.pi / YEAR
        ok = np.isfinite(a) & np.isfinite(ph)
        zc = np.sum(a[ok] * np.exp(1j * ph[ok]))
        coh = abs(zc) / np.sum(a[ok])
        print(f"  {k_p}: amplitude-weighted circular mean DOY {(np.angle(zc) / 2 / np.pi * YEAR) % YEAR:.0f}, "
              f"phase coherence {coh:.2f}, coherent network amplitude {abs(zc) / ok.sum():.2f} %")
    st = [s for s in steps if s["stream"] != PAYERNE and not s.get("placebo")]
    pl = np.array([s["dlnC"] for s in steps if s["stream"] != PAYERNE and s.get("placebo")])
    if pl.size:
        print(f"  window-step placebo (no window change): n {pl.size}, median dlnC {np.median(pl):+.1f} %, "
              f"IQR {np.percentile(pl, 25):+.1f} .. {np.percentile(pl, 75):+.1f}")
    if st:
        xs = np.array([s["dlnw"] for s in st])
        ys = np.array([s["dlnC"] for s in st])
        ts = theilslopes(ys, xs)
        print(f"  window steps Theil-Sen slope {ts.slope:+.2f} [{ts.low_slope:+.2f}, {ts.high_slope:+.2f}]")
        x = np.array([s["dlnw"] for s in st])
        y = np.array([s["dlnC"] for s in st])
        for name, sel in (("up (cleaning)", x > 0), ("down (soiling)", x < 0), ("all", np.ones_like(x, bool))):
            if sel.sum() >= 5:
                cf, se_, _ = ols(np.column_stack([np.ones(sel.sum()), x[sel]]), y[sel])
                print(f"  window steps {name}: n {sel.sum()} units {len({s['stream'] for s, k in zip(st, sel) if k})}"
                      f" median dlnw {np.median(x[sel]):+.1f} %, median dlnC {np.median(y[sel]):+.1f} %, "
                      f"slope {cf[1]:+.2f} +- {se_[1]:.2f}")
    for r in pair_rows:
        print(f"  pair {r['pair']:24s} {r['hours']:5s} ctrl {r['control']:13s} n {r['n']:5d} slope {r['g_pcal']:+.3f}"
              f" +- {r['se_g_pcal']:.3f}  sd_x {r['sd_dlnp']:.1f} %  kappa_dT {r['kappa_dT']:+.2f}")
    # named examples
    for key in ("0-20000-0-06240_A", "0-20000-0-06240_B", "0-20000-0-12843_A", "0-20000-0-12843_B",
                "0-20008-0-INO_A", "0-20000-0-03808_A", "0-20000-0-10393_A", "0-100-20000-000627_A"):
        for r in rows:
            if r["stream"] == key:
                print(f"  {key:24s} {r['site'][:16]:16s} g_pcal {r['g_pcal']:+.2f}+-{r['se_g_pcal']:.2f} "
                      f"kappa_pcal {r['kappa_pcal']:+.2f} amp_pcal {r['amp_pcal']:.1f}@{r['peak_pcal']:.0f} "
                      f"trend {r['trend_pcal']:+.1f}%/yr amp_C {r['amp_C']:.1f}@{r['peak_C']:.0f} "
                      f"g_win {r['g_win']:+.2f} amp_win {r['amp_win']:.1f}@{r['peak_win']:.0f} "
                      f"noon pcal {r['calibration_pulse_noon']:.1f}")
    figures(rows, series, steps, pair_rows)
    print(f"-> {Q15}")


if __name__ == "__main__":
    main()
