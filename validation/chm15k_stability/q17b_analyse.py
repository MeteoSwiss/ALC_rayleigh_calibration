# -*- coding: utf-8 -*-
"""q17b (analysis): the wavelength-scaling test CHM15k (1064 nm) vs CL61 (910.55 nm).

Model. For a night n and a height band z, with ln X = ln(rcs / (beta_mol T_mol^2 [T_wv^2])):
    ln X_CHM = g_CHM(t) + a(z, t)                       a = ln R_1064 (aerosol in the 'molecular' band)
    ln X_CL61 = g_CL61(t) + r * a(z, t),   r = (910.55/1064)^(4-AE) = 0.535 (AE 0) .. 0.62 (AE 1) .. 0.68 (AE 1.5)
An aerosol cycle gives CL61/CHM15k = r ~ 0.6 for the seasonal (and synoptic) variations; an instrument
gain of the CHM15k gives ~0. The fraction of the CHM15k cycle that is aerosol is f = s / r, where s is
the CL61 response per unit CHM15k response (in-phase harmonic ratio, monthly regression slope).

Per site and band (common clear nights, same blocks):
  * annual harmonic (offset + cos + sin [+ trend if span > 540 d]) of ln X_CHM (raw and kappa-normalised
    with the unit's intra-night kappa), ln X_CL61 (WV-corrected, + variants), and the ratio; CL61 in-phase
    amplitude projected on the CHM15k phase;
  * seasonal slope: monthly medians, OLS CL61 on CHM15k, month bootstrap CI;
  * synoptic slope: anomalies vs a centred 31-day running median, OLS + errors-in-variables correction,
    correlation; event nights (top-decile CHM15k positive anomalies): mean CL61 / mean CHM15k anomaly;
  * variants: CL61 without WV correction, WV humidity fill 'q', CL61 dark template, AERONET below-band
    aerosol transmission, CL61 gain monitored by its own cloud (liquid-cloud, aerosol-immune) constants.
Also: CL61 Rayleigh vs cloud constants (task 4) and the CHM15k full-record cycle (CHM15k-only nights).

Outputs: q17/q17b_bands.csv, q17b_attribution.csv, q17b_constants.csv, q17b_*.png
"""
from __future__ import annotations

import csv
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from common import NC_DIR  # noqa: E402
from q17b_extract import Q17, SITES  # noqa: E402
from q17b_nights import BAND_NAMES, BANDS  # noqa: E402

warnings.filterwarnings("ignore")
W = 2 * np.pi / 365.25
EPOCH = datetime(1970, 1, 1)
R_AE = {"AE0": (910.55 / 1064.0) ** 4 / (910.55 / 1064.0) ** 0,   # placeholder, replaced below
        }
LAM = 910.55 / 1064.0


# CL61 artefact periods excluded from the comparison. Lindenberg Jan-Apr 2025: CL61 ln X is +50..+134 %
# above its own median in EVERY band (2-10 km) while the co-located CHM15k is in its winter low and no CL61
# cloud constant exists - an instrument/processing state of the CL61, not air (q16g's 'spring artefact').
CL61_EXCLUDE_BEFORE = {"Lindenberg": 20209.0}      # 2025-05-01 (days since 1970)


def cl61_valid(site, t):
    return np.asarray(t, float) >= CL61_EXCLUDE_BEFORE.get(site, -np.inf)


def r_of_ae(ae):
    """CL61 / CHM15k ratio of the aerosol excess R-1 (small-excess limit): (lambda_CL/lambda_CHM)^(4-AE)."""
    return LAM ** (4.0 - ae)


R_LO, R_MID, R_HI = r_of_ae(0.0), r_of_ae(1.0), r_of_ae(1.5)
rng = np.random.default_rng(17)


def flag_iqr(v, k=1.5):
    v = np.asarray(v, float)
    q25, q75 = np.nanpercentile(v, 25), np.nanpercentile(v, 75)
    c = k * (q75 - q25)
    return (v < q25 - c) | (v > q75 + c)


def harmonic(t, ys, trend):
    """Joint least squares of several series on the same nights. Returns list of complex annual
    phasors (cos + i sin coefficient: amplitude |p|, peak DOY = angle/W)."""
    cols = [np.ones_like(t), np.cos(W * np.mod(t, 365.25)), np.sin(W * np.mod(t, 365.25))]
    if trend:
        cols.append((t - t.mean()) / 365.25)
    M = np.column_stack(cols)
    out = []
    for y in ys:
        c, *_ = np.linalg.lstsq(M, y, rcond=None)
        out.append(c[1] + 1j * c[2])
    return out


def peak(p):
    return float(np.mod(np.angle(p) / W, 365.25)) if np.isfinite(p) else np.nan


def running_median(t, y, half=15.5, min_n=5):
    out = np.full(len(t), np.nan)
    for k in range(len(t)):
        m = np.abs(t - t[k]) <= half
        if m.sum() >= min_n:
            out[k] = np.median(y[m])
    return out


def month_key(t):
    return np.array([(EPOCH + timedelta(days=float(x))).strftime("%Y%m") for x in t])


def monthly(t, cols, min_n=4):
    mk = month_key(t)
    rows = []
    for m in np.unique(mk):
        s = mk == m
        if s.sum() >= min_n:
            rows.append([np.median(t[s])] + [np.median(c[s]) for c in cols] + [s.sum()])
    return np.array(rows) if rows else np.zeros((0, len(cols) + 2))


def ols_slope(x, y):
    x = x - x.mean()
    y = y - y.mean()
    vx = np.sum(x * x)
    return float(np.sum(x * y) / vx) if vx > 0 else np.nan


def load_kappa_in():
    from q16b_insolation import load_kappas
    k_in, _, _ = load_kappas()
    return k_in


def read_constants(stream, method):
    from q16g_cl61 import read_nc
    return read_nc(NC_DIR / f"ALC_calibration_{stream}20*.nc", method)


def cloud_gain(site, t):
    """CL61 gain from its own liquid-cloud constants: 45-day running median of ln C_cloud at time t."""
    wmo, _, icl, _, _ = SITES[site]
    d = read_constants(f"{wmo}_{icl}", 1)
    if d is None or len(d["t"]) < 15:
        return np.full(len(t), np.nan)
    tc, lc = d["t"], np.log(d["c"])
    bad = flag_iqr(lc)
    tc, lc = tc[~bad], lc[~bad]
    out = np.full(len(t), np.nan)
    for k, x in enumerate(t):
        m = np.abs(tc - x) <= 22.5
        if m.sum() >= 4:
            out[k] = np.median(lc[m])
    return out


# ---------------------------------------------------------------------------------------------
def analyse_series(t, x, y, sx, sy, trend):
    """All the scaling statistics of CL61 (y) vs CHM15k (x) on common nights."""
    res = {}
    px, py, pr = harmonic(t, [x, y, x - y], trend)
    res.update(amp_x=100 * abs(px), peak_x=peak(px), amp_y=100 * abs(py), peak_y=peak(py),
               amp_ratio=100 * abs(pr), peak_ratio=peak(pr),
               amp_y_inphase=100 * float(np.real(py * np.conj(px)) / abs(px)),
               s_harm=float(np.real(py * np.conj(px)) / abs(px) ** 2))
    # monthly regression
    mo = monthly(t, [x, y])
    res["n_months"] = len(mo)
    if len(mo) >= 5:
        res["s_month"] = ols_slope(mo[:, 1], mo[:, 2])
        bs = []
        for _ in range(2000):
            i = rng.integers(0, len(mo), len(mo))
            if np.ptp(mo[i, 1]) > 0:
                bs.append(ols_slope(mo[i, 1], mo[i, 2]))
        res["s_month_lo"], res["s_month_hi"] = np.percentile(bs, [2.5, 97.5])
        res["range_month_x"] = 100 * float(np.ptp(mo[:, 1]))
        res["r_month"] = float(np.corrcoef(mo[:, 1], mo[:, 2])[0, 1])
        # harmonic-ratio bootstrap by month (resample calendar months with their nights)
        mk = month_key(t)
        um = np.unique(mk)
        bsh = []
        for _ in range(500):
            pick = rng.choice(um, len(um))
            idx = np.concatenate([np.where(mk == m)[0] for m in pick])
            if np.ptp(t[idx]) < 200:
                continue
            a, b = harmonic(t[idx], [x[idx], y[idx]], trend and np.ptp(t[idx]) > 540)
            bsh.append(float(np.real(b * np.conj(a)) / abs(a) ** 2))
        if bsh:
            res["s_harm_lo"], res["s_harm_hi"] = np.percentile(bsh, [2.5, 97.5])
    # synoptic anomalies
    ax = x - running_median(t, x)
    ay = y - running_median(t, y)
    ok = np.isfinite(ax) & np.isfinite(ay)
    if ok.sum() >= 30:
        a_, b_ = ax[ok], ay[ok]
        res["n_syn"] = int(ok.sum())
        res["sd_syn_x"] = 100 * float(np.std(a_))
        res["sd_syn_y"] = 100 * float(np.std(b_))
        res["r_syn"] = float(np.corrcoef(a_, b_)[0, 1])
        res["s_syn"] = ols_slope(a_, b_)
        noise = np.nanmean(sx[ok] ** 2) if sx is not None else 0.0
        vt = np.var(a_) - noise
        res["noise_frac_x"] = float(noise / np.var(a_))
        res["s_syn_eiv"] = float(np.cov(a_, b_)[0, 1] / vt) if vt > 0.2 * np.var(a_) else np.nan
        thr = np.percentile(a_, 90)
        ev = a_ >= thr
        res["ev_x"] = 100 * float(a_[ev].mean())
        res["ev_y"] = 100 * float(b_[ev].mean())
        res["s_event"] = float(b_[ev].mean() / a_[ev].mean())
        bse = []
        for _ in range(1000):
            i = rng.integers(0, len(a_), len(a_))
            aa, bb = a_[i], b_[i]
            e = aa >= np.percentile(aa, 90)
            bse.append(bb[e].mean() / aa[e].mean())
        res["s_event_lo"], res["s_event_hi"] = np.percentile(bse, [2.5, 97.5])
        bss = []
        for _ in range(1000):
            i = rng.integers(0, len(a_), len(a_))
            bss.append(ols_slope(a_[i], b_[i]))
        res["s_syn_lo"], res["s_syn_hi"] = np.percentile(bss, [2.5, 97.5])
    return res


def site_band(df, site, bn, ycol, kappa, extra=None):
    d = df[df.site == site]
    t = d["t"].to_numpy(float)
    x = d[f"chm_{bn}"].to_numpy(float)
    if site == "Payerne" and f"chmdark_{bn}" in d:
        pass
    y = d[f"{ycol}_{bn}"].to_numpy(float)
    sx = d[f"sig_chm_{bn}"].to_numpy(float)
    sy = d[f"sig_cl_{bn}"].to_numpy(float)
    T = d["T_int"].to_numpy(float)
    if extra == "aeronet":
        x = x + 2 * d["tau1064"].to_numpy(float)
        y = y + 2 * d["tau910"].to_numpy(float)
    if extra == "cloudgain":
        y = y - cloud_gain(site, t)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(t) & cl61_valid(site, t)
    t, x, y, sx, sy, T = t[ok], x[ok], y[ok], sx[ok], sy[ok], T[ok]
    bad = flag_iqr(x) | flag_iqr(y)
    t, x, y, sx, sy, T = t[~bad], x[~bad], y[~bad], sx[~bad], sy[~bad], T[~bad]
    if len(t) < 25 or np.ptp(t) < 150:
        return None
    trend = np.ptp(t) > 540
    res = dict(site=site, band=bn, ycol=ycol + ("" if extra is None else "+" + extra), n=len(t),
               span_d=float(np.ptp(t)), t_first=(EPOCH + timedelta(days=float(t.min()))).strftime("%Y-%m-%d"),
               t_last=(EPOCH + timedelta(days=float(t.max()))).strftime("%Y-%m-%d"),
               sig_x_med=100 * float(np.nanmedian(sx)), sig_y_med=100 * float(np.nanmedian(sy)))
    res.update(analyse_series(t, x, y, sx, sy, trend))
    if np.isfinite(kappa):
        xk = x - kappa / 100.0 * (T - np.nanmedian(T))
        rk = analyse_series(t, xk, y, sx, sy, trend)
        res.update({f"k_{k}": v for k, v in rk.items() if k.startswith(("amp_", "peak_", "s_harm", "s_month"))})
    return res


def constants_check():
    """Task 4: CL61 Rayleigh vs cloud constants, and the CHM15k vs CL61 cloud constant, per site."""
    from q16g_cl61 import harm, match_days, screened
    from common import load_constants
    rows = []
    for site, (wmo, ich, icl, _, _) in SITES.items():
        ray = read_constants(f"{wmo}_{icl}", 0)
        cld = read_constants(f"{wmo}_{icl}", 1)
        if ray is None or cld is None:
            continue
        tr, lr = screened(0.5 * (ray["t0"] + ray["t1"]), ray["c"])
        tc, lc = screened(cld["t"], cld["c"])
        row = dict(site=site, n_ray=len(tr), n_cloud=len(tc))
        for tag, (t, y) in (("cl61_ray", (tr, lr)), ("cl61_cloud", (tc, lc))):
            a, p, n = harm(t, y, trend=False)
            row[f"amp_{tag}"], row[f"peak_{tag}"] = a, p
        ia, ib = match_days(tr, tc, 3.0)
        if len(ia) >= 20:
            a, p, n = harm(tr[ia], lr[ia] - lc[ib], trend=False)
            row.update(amp_ray_over_cloud=a, peak_ray_over_cloud=p, n_ray_over_cloud=n,
                       med_ray_over_cloud=float(np.exp(np.median(lr[ia] - lc[ib]))))
            # the in-phase (solstice) component of Rayleigh/cloud
            t = tr[ia]
            yy = lr[ia] - lc[ib]
            ph = harmonic(t, [yy], False)[0]
            row["inphase176_ray_over_cloud"] = 100 * float(np.real(ph * np.exp(-1j * W * 176.0)))
        ch = load_constants(wmo, ich)
        if ch is not None:
            tch, lch = screened(0.5 * (ch["t0"] + ch["t1"]), ch["c"])
            ia, ib = match_days(tch, tc, 3.0)
            if len(ia) >= 20:
                a, p, n = harm(tch[ia], lch[ia] - lc[ib], trend=False)
                row.update(amp_chm_over_cl61cloud=a, peak_chm_over_cl61cloud=p, n_chm_over_cl61cloud=n)
                a, p, n = harm(tch[ia], lch[ia], trend=False)
                row.update(amp_chm_same_nights=a, peak_chm_same_nights=p)
        rows.append(row)
        print("constants", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    return rows


def chm_full_record(dfa, k_in):
    rows = []
    for site, (wmo, ich, _, _, _) in SITES.items():
        d = dfa[dfa.site == site]
        kap = k_in.get(f"{wmo}_{ich}", np.nan)
        for bn in BAND_NAMES:
            t = d["t"].to_numpy(float)
            x = d[f"chm_{bn}"].to_numpy(float)
            T = d["T_int"].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(T)
            t, x, T = t[ok], x[ok], T[ok]
            bad = flag_iqr(x)
            t, x, T = t[~bad], x[~bad], T[~bad]
            if len(t) < 50:
                continue
            xk = x - kap / 100 * (T - np.median(T)) if np.isfinite(kap) else x * np.nan
            p, pk, pT = harmonic(t, [x, xk, T], np.ptp(t) > 540)
            rows.append(dict(site=site, band=bn, n=len(t), kappa_in=kap, amp=100 * abs(p), peak=peak(p),
                             amp_k=100 * abs(pk), peak_k=peak(pk), amp_T=abs(pT), peak_T=peak(pT)))
    return rows


def write_csv(rows, path):
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
    df = pd.read_csv(Q17 / "q17b_nights.csv")
    dfa = pd.read_csv(Q17 / "q17b_nights_chm.csv")
    k_in = load_kappa_in()
    rows = []
    for site, (wmo, ich, _, _, _) in SITES.items():
        kap = k_in.get(f"{wmo}_{ich}", np.nan)
        for bn in BAND_NAMES:
            for ycol, extra in (("cl_wv", None), ("cl_nowv", None), ("cl_wvq", None), ("cl_wvdark", None),
                                ("cl_wv", "aeronet"), ("cl_wv", "cloudgain")):
                r = site_band(df, site, bn, ycol, kap, extra)
                if r is not None:
                    r["kappa_in"] = kap
                    rows.append(r)
                    if ycol == "cl_wv" and extra is None:
                        print(site, bn, {k: round(r[k], 2) for k in ("n", "amp_x", "peak_x", "amp_y", "peak_y",
                                                                       "s_harm", "s_month", "s_syn", "r_syn",
                                                                       "s_event") if k in r}, flush=True)
    write_csv(rows, Q17 / "q17b_bands.csv")
    write_csv(chm_full_record(dfa, k_in), Q17 / "q17b_chm_fullrecord.csv")
    write_csv(constants_check(), Q17 / "q17b_constants.csv")


if __name__ == "__main__":
    main()
