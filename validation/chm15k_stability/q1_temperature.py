# -*- coding: utf-8 -*-
"""Q1 of report 19: does the CHM15k Rayleigh constant depend on instrument temperature?

The internal temperature follows the season (74 % annual cycle at Payerne), and so do the clear-night
aerosol load, night length and molecular SNR. A naive regression log C ~ T therefore mixes season
with temperature. Three models are fitted per unit and optical-module era:

  M0  naive        log C = a + b T
  M1  season       log C = a + b T + annual harmonic (sin/cos DOY) + linear trend (laser ageing)
  M2  anomaly      (log C - 31-d running median)  vs  (T - 31-d running median)

M2 isolates synoptic temperature anomalies (warm or cold spells relative to the local seasonal
state) and is the cleanest estimate of a real temperature effect on the constant.

Two arms:
  A  network 2025-2026: operational nightly constants (dashboard NetCDF) matched to the 10-min
     housekeeping cache (hk_cache.py) over each night's calibration window, ~150 CHM15k;
  B  the 11-year daily archive A:/CHM15k (13 stations, 2013-2024) exported by export_11yr_archive.m,
     split at optical-module swaps (serlom).

Run:  python validation/chm15k_stability/q1_temperature.py
"""
from __future__ import annotations

import csv
import sys
from datetime import datetime, timedelta

EPOCH_DT = datetime(1970, 1, 1)
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (BLOCKS_PER_DAY, FIGS, OUT, chm15k_streams, load_cache, load_constants,  # noqa: E402
                    stream_key)

MIN_NIGHTS = 40          # per unit-era, after outlier screening
MIN_T_RANGE = 5.0        # K, p5-p95 of the regressor; below this a slope is not identifiable
ANOM_WINDOW = 31.0       # days, running-median window of the M2 anomaly model
ARCHIVE = OUT / "archive_11yr"
TVARS = ("temp_int", "temp_ext", "temperature_detector", "temperature_optical_module")
EXTRA = ("calibration_pulse", "window_transmission", "status_laser")


# ---------------------------------------------------------------------------------------------
# Regression helpers
# ---------------------------------------------------------------------------------------------
def robust_mask(y, nmad=4.0):
    """Drop points more than nmad MADs from the median (one pass, on log C)."""
    med = np.nanmedian(y)
    mad = 1.4826 * np.nanmedian(np.abs(y - med))
    if not np.isfinite(mad) or mad == 0:
        return np.isfinite(y)
    return np.isfinite(y) & (np.abs(y - med) <= nmad * mad)


def ols(X, y):
    """OLS with classical standard errors. Returns coef, se, residuals."""
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    dof = max(len(y) - X.shape[1], 1)
    s2 = res @ res / dof
    try:
        cov = s2 * np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        se = np.full(X.shape[1], np.nan)
    return coef, se, res


def running_median(t, x, window):
    """Centred running median of x over +-window/2 days (irregular sampling)."""
    out = np.full_like(x, np.nan, dtype=float)
    half = window / 2.0
    lo = np.searchsorted(t, t - half)
    hi = np.searchsorted(t, t + half, side="right")
    for i in range(len(t)):
        seg = x[lo[i]:hi[i]]
        seg = seg[np.isfinite(seg)]
        if seg.size >= 5:
            out[i] = np.median(seg)
    return out


def fit_models(t, logc, T, extra=None):
    """Fit M0/M1/M2 (and M1+extra regressors) of log C on T. t in days. Slopes in %/K."""
    ok = np.isfinite(t) & np.isfinite(logc) & np.isfinite(T)
    t, logc, T = t[ok], logc[ok], T[ok]
    extra = None if extra is None else {k: v[ok] for k, v in extra.items()}
    keep = robust_mask(logc)
    t, logc, T = t[keep], logc[keep], T[keep]
    extra = None if extra is None else {k: v[keep] for k, v in extra.items()}
    n = len(t)
    if n < MIN_NIGHTS:
        return None
    trange = np.percentile(T, 95) - np.percentile(T, 5)
    if trange < MIN_T_RANGE:
        return None
    order = np.argsort(t)
    t, logc, T = t[order], logc[order], T[order]
    extra = None if extra is None else {k: v[order] for k, v in extra.items()}
    one = np.ones(n)
    doy = 2 * np.pi * (t % 365.25) / 365.25
    yrs = (t - t.mean()) / 365.25
    out = dict(n=n, t_range=trange, t_mean=float(T.mean()),
               span_d=float(t[-1] - t[0]), std_logc=100 * float(np.std(logc)))

    c, se, _ = ols(np.column_stack([one, T]), logc)
    out.update(b0=100 * c[1], se0=100 * se[1], r0=float(np.corrcoef(T, logc)[0, 1]))

    # Season + ageing control. The harmonic needs at least ~1 year to be identifiable;
    # for shorter records only the trend is kept (season then partly loads on T - flagged).
    full_year = out["span_d"] >= 330
    cols = [one, T, yrs] + ([np.sin(doy), np.cos(doy)] if full_year else [])
    c, se, _ = ols(np.column_stack(cols), logc)
    out.update(b1=100 * c[1], se1=100 * se[1], season_ctrl=bool(full_year))

    # Partial correlation of T with the seasonal harmonic: how collinear the two are.
    if full_year:
        cs, _, _ = ols(np.column_stack([one, np.sin(doy), np.cos(doy)]), T)
        out["t_seasonal_frac"] = float(1 - np.var(T - np.column_stack([one, np.sin(doy), np.cos(doy)]) @ cs)
                                       / np.var(T))
    else:
        out["t_seasonal_frac"] = np.nan

    la = logc - running_median(t, logc, ANOM_WINDOW)
    ta = T - running_median(t, T, ANOM_WINDOW)
    ok2 = np.isfinite(la) & np.isfinite(ta)
    if ok2.sum() >= MIN_NIGHTS and np.std(ta[ok2]) > 0.3:
        c, se, _ = ols(np.column_stack([np.ones(ok2.sum()), ta[ok2]]), la[ok2])
        out.update(b2=100 * c[1], se2=100 * se[1], n2=int(ok2.sum()),
                   ta_std=float(np.std(ta[ok2])), r2=float(np.corrcoef(ta[ok2], la[ok2])[0, 1]))
    else:
        out.update(b2=np.nan, se2=np.nan, n2=int(ok2.sum()), ta_std=np.nan, r2=np.nan)

    # M2 + atmospheric confound controls: warm spells are often hazier, and aerosol leaking into the
    # "molecular" window inflates C. Anomalies of the window bottom height and of the low-level
    # aerosol contrast ln(S_1-2.5km / S_window) are added; if b2 survives them it is not haze.
    if extra and "zb" in extra and ok2.sum() >= MIN_NIGHTS:
        za = extra["zb"] - running_median(t, extra["zb"], ANOM_WINDOW)
        aa = extra["aer"] - running_median(t, extra["aer"], ANOM_WINDOW)
        ok3 = ok2 & np.isfinite(za) & np.isfinite(aa)
        if ok3.sum() >= MIN_NIGHTS:
            X = np.column_stack([np.ones(ok3.sum()), ta[ok3], za[ok3] / 1000.0, aa[ok3]])
            c, se, _ = ols(X, la[ok3])
            out.update(b2x=100 * c[1], se2x=100 * se[1], n2x=int(ok3.sum()),
                       g_zb=100 * c[2], g_aer=100 * c[3],
                       r_ta_aer=float(np.corrcoef(ta[ok3], aa[ok3])[0, 1]))
        extra = {k: v for k, v in extra.items() if k not in ("zb", "aer")}

    # M1 + instrumental regressors (standardised) -> partial T slope with HK held fixed.
    if extra:
        xs = []
        for k, v in extra.items():
            if np.isfinite(v).all() and np.std(v) > 0:
                xs.append((v - v.mean()) / np.std(v))
        if xs:
            c, se, _ = ols(np.column_stack(cols + xs), logc)
            out.update(b1x=100 * c[1], se1x=100 * se[1])
    out["_series"] = (t, logc, T, la, ta)
    return out


# ---------------------------------------------------------------------------------------------
# Arm A: network 2025-2026
# ---------------------------------------------------------------------------------------------
def night_hk(cache, t0, t1):
    """Median of every HK variable over the calibration window [t0, t1] (days since 1970)."""
    tb = (cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY).ravel()
    hk = cache["hk"].reshape(-1, cache["hk"].shape[-1])
    order = np.argsort(tb)
    tb, hk = tb[order], hk[order]
    out = np.full((len(t0), hk.shape[1]), np.nan)
    i0 = np.searchsorted(tb, t0)
    i1 = np.searchsorted(tb, t1, side="right")
    for k in range(len(t0)):
        if i1[k] > i0[k]:
            out[k] = np.nanmedian(hk[i0[k]:i1[k]], axis=0)
    return out


def night_aerosol_contrast(cache, cons, alt):
    """ln(clear-sky signal at 1.0-2.5 km / signal in the night's molecular window), per night."""
    nb = cache["rcs"].shape[-1]
    tb = (cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY).ravel()
    order = np.argsort(tb)
    tb = tb[order]
    rcs = cache["rcs"].reshape(-1, nb)[order]
    clear = (cache["n_clear"].ravel()[order] >= 0.9 * np.maximum(cache["n_all"].ravel()[order], 1))
    rmid = 0.5 * (cache["range_edges"][:-1] + cache["range_edges"][1:])
    low = (rmid >= 1000) & (rmid <= 2500)
    out = np.full(len(cons["t"]), np.nan)
    i0 = np.searchsorted(tb, cons["t0"])
    i1 = np.searchsorted(tb, cons["t1"], side="right")
    for k in range(len(out)):
        win = (rmid >= cons["zb"][k] - alt) & (rmid <= cons["zt"][k] - alt)
        if i1[k] <= i0[k] or win.sum() < 1:
            continue
        seg = rcs[i0[k]:i1[k]][clear[i0[k]:i1[k]]]
        if seg.shape[0] < 6:
            continue
        with np.errstate(invalid="ignore", divide="ignore"):
            s_low = np.nanmean(seg[:, low])
            s_win = np.nanmean(seg[:, win])
            if s_low > 0 and s_win > 0:
                out[k] = np.log(s_low / s_win)
    return out


def arm_a():
    import json
    from common import CENSUS
    alts = {(s["wmo"], s["ident"]): float(s["alt"]) for s in json.load(open(CENSUS, encoding="utf-8"))}
    rows, series = [], {}
    for wmo, ident, site, lat, lon in chm15k_streams():
        cons = load_constants(wmo, ident)
        cache = load_cache(wmo, ident)
        if cons is None or cache is None:
            continue
        aer = night_aerosol_contrast(cache, cons, alts[(wmo, ident)])
        hk = night_hk(cache, cons["t0"], cons["t1"])
        names = list(cache["hk_vars"])
        col = {v: hk[:, names.index(v)] for v in names}
        # Regulation-loss excursion: optical module / detector above the unit's own set-point,
        # i.e. above the 20th percentile of the whole record (the regulated floor).
        allhk = cache["hk"].reshape(-1, len(names))
        for v in ("temperature_optical_module", "temperature_detector"):
            floor = np.nanpercentile(allhk[:, names.index(v)], 20)
            col[v + "_exc"] = np.clip(col[v] - floor, 0, None)
        logc = np.log(np.where(cons["c"] > 0, cons["c"], np.nan))
        for mod in np.unique(cons["module"]):
            m = cons["module"] == mod
            for tv in TVARS:
                ex = None
                if tv == "temp_int":
                    ex = {k: col[k][m] for k in EXTRA}
                    ex.update(zb=cons["zb"][m], aer=aer[m])
                res = fit_models(cons["t"][m], logc[m], col[tv][m], extra=ex)
                if res is None:
                    continue
                key = (stream_key(wmo, ident), mod, tv)
                series[key] = res.pop("_series")
                rows.append(dict(arm="A", stream=stream_key(wmo, ident), site=site, lat=lat, lon=lon,
                                 era=mod, tvar=tv, **res,
                                 lom_exc_frac=float(np.mean(col["temperature_optical_module_exc"][m] > 1.0)),
                                 det_exc_frac=float(np.mean(col["temperature_detector_exc"][m] > 1.0))))
        print(f"  A {stream_key(wmo, ident):28s} {site[:20]:20s} nights {len(cons['t']):4d}", flush=True)
    return rows, series


# ---------------------------------------------------------------------------------------------
# Payerne dark arm: operational vs dark + background corrected constants (the only hooded CHM15k)
# ---------------------------------------------------------------------------------------------
PAY_CORR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/diag_v22_darkbg_v2/0-20000-0-06610_A/"
                "0-20000-0-06610_A_cal.csv")


def arm_payerne_dark():
    """Same three models on the Payerne CHM15k, operational constant vs measured-dark-corrected one.

    The hood dark has a slow thermal component (D_slow ~21 rcs/K at 60-500 m); an additive baseline
    that moves with temperature moves the Rayleigh constant. If the T slope changes when the measured
    dark is removed, that part of it is the dark, not the gain.
    """
    wmo, ident = "0-20000-0-06610", "A"
    cons = load_constants(wmo, ident)
    cache = load_cache(wmo, ident)
    if cons is None or cache is None or not PAY_CORR.exists():
        return []
    corr = {}
    for r in csv.DictReader(open(PAY_CORR, encoding="utf-8")):
        if r["method"] == "rayleigh" and r["flag"] in ("1", "0.5"):
            corr[r["date"][:8]] = float(r["cal_value"])
    # Same date key as rayleigh_availability/chm15k_constant_vs_temperature.py: the date of the
    # calibration's central time.
    days = [(EPOCH_DT + timedelta(days=float(x))).strftime("%Y%m%d") for x in cons["t"]]
    c_corr = np.array([corr.get(d, np.nan) for d in days])
    hk = night_hk(cache, cons["t0"], cons["t1"])
    T = hk[:, list(cache["hk_vars"]).index("temp_int")]
    rows = []
    for label, c in (("operational", cons["c"]), ("dark+bg corrected", c_corr)):
        both = np.isfinite(cons["c"]) & np.isfinite(c_corr)       # identical night sample for both
        logc = np.log(np.where((c > 0) & both, c, np.nan))
        for mod in np.unique(cons["module"]):
            m = cons["module"] == mod
            res = fit_models(cons["t"][m], logc[m], T[m])
            if res is None:
                continue
            res.pop("_series")
            rows.append(dict(arm="PAYDARK", stream=f"{wmo}_{ident}", site=f"PAYERNE {label}", lat=np.nan,
                             lon=np.nan, era=mod, tvar="temp_int", **res, lom_exc_frac=np.nan, det_exc_frac=np.nan))
            print(f"  PAY {label:18s} {mod}: n={res['n']} b0 {res['b0']:+.2f}+/-{res['se0']:.2f}  "
                  f"b1 {res['b1']:+.2f}+/-{res['se1']:.2f}  b2 {res['b2']:+.2f}+/-{res['se2']:.2f} %/K", flush=True)
    return rows


# ---------------------------------------------------------------------------------------------
# Arm B: 11-year archive
# ---------------------------------------------------------------------------------------------
def arm_b():
    rows, series = [], {}
    for f in sorted(ARCHIVE.glob("chm_*_1d.csv")):
        stn = f.stem.split("_")[1]
        t, c, T, ser, pc, so = [], [], [], [], [], []
        for r in csv.DictReader(open(f, encoding="utf-8")):
            try:
                d = datetime.strptime(r["Time"], "%d-%b-%Y")
            except ValueError:
                continue
            t.append((d - datetime(1970, 1, 1)).days + 0.5)
            c.append(float(r["Lidar_constant"]) if r["Lidar_constant"] not in ("", "NaN") else np.nan)
            T.append(float(r["temp_int"]) - 273.15 if r["temp_int"] not in ("", "NaN") else np.nan)
            ser.append(r["serlom"])
            pc.append(float(r["p_calc"]) if r["p_calc"] not in ("", "NaN") else np.nan)
            so.append(float(r["state_optics"]) if r["state_optics"] not in ("", "NaN") else np.nan)
        t, c, T, ser, pc, so = map(np.array, (t, c, T, ser, pc, so))
        good = (T > -23) & (c > 0) & (c < 1e14)        # same physical screens as the MATLAB script
        logc = np.log(np.where(good, c, np.nan))
        for mod in np.unique(ser[good]):
            m = (ser == mod) & good
            res = fit_models(t[m], logc[m], T[m], extra={"calibration_pulse": pc[m], "window_transmission": so[m]})
            if res is None:
                continue
            series[(stn, mod, "temp_int")] = res.pop("_series")
            rows.append(dict(arm="B", stream=stn, site=stn.upper(), lat=np.nan, lon=np.nan, era=mod,
                             tvar="temp_int", **res, lom_exc_frac=np.nan, det_exc_frac=np.nan))
        print(f"  B {stn}: {good.sum()} days", flush=True)
    return rows, series


# ---------------------------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------------------------
def write_csv(rows, path):
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("arm", "stream", "site", "era", "tvar"), k))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def summarise(rows):
    print("\n=== Slope of log C_L on temperature [%/K], per unit-era: median [IQR], n significant (|b|>3se) ===")
    for arm in ("A", "B"):
        for tv in TVARS:
            sel = [r for r in rows if r["arm"] == arm and r["tvar"] == tv]
            if not sel:
                continue
            line = f"  arm {arm} {tv:27s} N={len(sel):3d}"
            for m in ("0", "1", "2", "2x", "1x"):
                b = np.array([r.get(f"b{m}", np.nan) for r in sel], float)
                se = np.array([r.get(f"se{m}", np.nan) for r in sel], float)
                ok = np.isfinite(b) & np.isfinite(se)
                if not ok.any():
                    continue
                q = np.percentile(b[ok], [25, 50, 75])
                nsig_p = int(np.sum(b[ok] > 3 * se[ok]))
                nsig_n = int(np.sum(b[ok] < -3 * se[ok]))
                # inverse-variance pooled estimate + its standard error
                w = 1 / se[ok] ** 2
                pooled = np.sum(w * b[ok]) / np.sum(w)
                line += (f" | M{m} {q[1]:+.2f} [{q[0]:+.2f},{q[2]:+.2f}] pooled {pooled:+.2f}"
                         f"+/-{1 / np.sqrt(np.sum(w)):.2f} sig+{nsig_p}/-{nsig_n}")
            print(line)


def figures(rows, series_a, series_b):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIGS.mkdir(parents=True, exist_ok=True)

    # Fig 1: per-unit slopes M0 / M1 / M2 (temp_int), arm A then arm B, sorted by M0.
    fig, axs = plt.subplots(2, 1, figsize=(20, 9), gridspec_kw=dict(height_ratios=[1, 0.6]))
    for ax, arm in zip(axs, ("A", "B")):
        sel = [r for r in rows if r["arm"] == arm and r["tvar"] == "temp_int"]
        sel.sort(key=lambda r: r["b0"])
        x = np.arange(len(sel))
        for m, col, lab, dx in (("0", "#d62728", "M0 naive", -0.25), ("1", "#1f77b4", "M1 season+trend", 0.0),
                                ("2", "#2ca02c", "M2 31-d anomaly", 0.25)):
            b = np.array([r[f"b{m}"] for r in sel], float)
            se = np.array([r[f"se{m}"] for r in sel], float)
            ax.errorbar(x + dx, b, yerr=se, fmt="o", ms=3, lw=0.8, color=col, label=lab, alpha=0.85)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{r['site'][:14]} {r['stream'][-1] if arm == 'A' else ''} {r['era'][-6:]}"
                            for r in sel], rotation=90, fontsize=5 if arm == "A" else 8)
        ax.set_ylabel("d ln C$_L$ / dT$_{int}$ [%/K]")
        ax.set_ylim(-8, 8)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, ncol=3, loc="upper left")
        ax.set_title("Arm A: network 2025-2026, operational nightly constants" if arm == "A" else
                     "Arm B: 11-year archive 2013-2024 (A:/CHM15k), split at optical-module swaps", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / "q1_slopes_per_unit.png", dpi=130)
    plt.close(fig)

    # Fig 2: distribution of slopes per model, temp_int + other temperatures (arm A).
    # Detector and optical-module temperatures are regulated (never a 5 K range) -> no slope; the
    # third panel shows the 11-year archive instead.
    fig, axs = plt.subplots(1, 3, figsize=(20, 5), sharey=False)
    bins = np.linspace(-8, 8, 41)
    for ax, (arm, tv, title) in zip(axs, (("A", "temp_int", "internal temperature, network 2025-2026"),
                                          ("A", "temp_ext", "external temperature, network 2025-2026"),
                                          ("B", "temp_int", "internal temperature, 11-yr archive 2013-2024"))):
        for m, col, lab in (("0", "#d62728", "M0 naive"), ("1", "#1f77b4", "M1 season"),
                            ("2", "#2ca02c", "M2 anomaly")):
            b = np.array([r[f"b{m}"] for r in rows if r["arm"] == arm and r["tvar"] == tv], float)
            b = b[np.isfinite(b)]
            if b.size:
                ax.hist(np.clip(b, -8, 8), bins=bins, histtype="step", lw=1.6, color=col,
                        label=f"{lab}: median {np.median(b):+.2f}")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("slope [%/K]")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("number of unit-eras")
    fig.suptitle("CHM15k Rayleigh constant vs instrument temperatures - network 2025-2026", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "q1_slope_histograms.png", dpi=130)
    plt.close(fig)

    # Fig 3: the season confound at a glance - examples, raw vs anomaly scatter.
    picks = []
    for want in ("0-20000-0-06610_A", "0-20000-0-10393_0", "0-20000-0-06240_A", "0-20008-0-INO_B",
                 "0-20000-0-12843_A"):
        for k in series_a:
            if k[0] == want and k[2] == "temp_int":
                picks.append(k)
                break
    if picks:
        fig, axs = plt.subplots(2, len(picks), figsize=(4.2 * len(picks), 7.5), squeeze=False)
        for j, k in enumerate(picks):
            t, logc, T, la, ta = series_a[k]
            doy = np.array([(EPOCH_DT + timedelta(days=float(x))).timetuple().tm_yday for x in t])
            ax = axs[0, j]
            sc = ax.scatter(T, 100 * (logc - np.median(logc)), c=doy, cmap="twilight", s=10, vmin=1, vmax=366)
            ax.set_xlabel("T$_{int}$ [\u00b0C]")
            ax.set_ylabel("ln C$_L$ - median [%]")
            ax.set_title(f"{k[0]}  {k[1]}\nraw (colour = day of year)", fontsize=9)
            ax.grid(alpha=0.3)
            ax = axs[1, j]
            ok = np.isfinite(la) & np.isfinite(ta)
            ax.scatter(ta[ok], 100 * la[ok], s=10, color="#2ca02c")
            if ok.sum() > 5:
                p = np.polyfit(ta[ok], 100 * la[ok], 1)
                xx = np.linspace(ta[ok].min(), ta[ok].max(), 10)
                ax.plot(xx, np.polyval(p, xx), "k--", lw=1, label=f"{p[0]:+.2f} %/K")
                ax.legend(fontsize=8)
            ax.set_xlabel("T$_{int}$ anomaly vs 31-d median [K]")
            ax.set_ylabel("ln C$_L$ anomaly [%]")
            ax.set_title("anomaly model M2", fontsize=9)
            ax.grid(alpha=0.3)
        plt.colorbar(sc, ax=axs[0, -1], label="day of year")
        fig.tight_layout()
        fig.savefig(FIGS / "q1_examples_raw_vs_anomaly.png", dpi=130)
        plt.close(fig)

    # Fig 4: 11-yr Payerne series (arm B) - constant and temperature, to show the seasonal co-variation.
    pk = [k for k in series_b if k[0] == "pay"]
    if pk:
        fig, ax = plt.subplots(figsize=(20, 4.8))
        ax2 = ax.twinx()
        for k in pk:
            t, logc, T, la, ta = series_b[k]
            d = [(EPOCH_DT + timedelta(days=float(x))) for x in t]
            ax.plot(d, np.exp(logc) / 1e11, ".", ms=2.5, label=f"module {k[1]}")
            ax2.plot(d, T, "-", color="#ff7f0e", lw=0.5, alpha=0.6)
        ax.set_ylabel("C$_L$ / 1e11 (archive)")
        ax2.set_ylabel("T$_{int}$ [\u00b0C]", color="#ff7f0e")
        ax.legend(fontsize=8, loc="upper left", ncol=4)
        ax.grid(alpha=0.3)
        ax.set_title("Payerne CHM15k 2013-2024: daily Rayleigh constant (dots, per module) and internal temperature",
                     fontsize=10)
        fig.tight_layout()
        fig.savefig(FIGS / "q1_payerne_11yr.png", dpi=130)
        plt.close(fig)


def main():
    print("Arm A (network 2025-2026)")
    rows_a, series_a = arm_a()
    print("Arm B (11-year archive)")
    rows_b, series_b = arm_b()
    print("Payerne dark arm")
    rows_p = arm_payerne_dark()
    rows = rows_a + rows_b + rows_p
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(rows, OUT / "q1_slopes.csv")
    summarise(rows)
    figures(rows, series_a, series_b)
    print(f"-> {OUT / 'q1_slopes.csv'}  figs in {FIGS}")


if __name__ == "__main__":
    main()
