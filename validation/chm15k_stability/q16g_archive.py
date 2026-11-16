# -*- coding: utf-8 -*-
"""q16g (task 1a): is the solstice-phased CHM15k gain cycle present EVERY year of the 11-year archive?

The 11-year archive (A:/CHM15k/<stn>, 13 stations, 2013-2024, exported to CSV by export_11yr_archive.m)
is an INDEPENDENT processing chain:
  * input = the raw Lufft NetCDF files (beta_raw), not E-PROFILE L1 / raw2l1;
  * calibration = the MeteoSwiss MATLAB Rayleigh code (auto_calib v0.19-0.252, CeiloCalib_*.nc), not the
    Python repository; molecular profile from NCEP T/p (get_TP_NCEP), not US-1976;
  * firmware 0.720 ... 1.130 over the period.
If the solstice-phased cycle is present every year with a stable amplitude/phase there, the 2025-2026
processing (raw2l1, hermes, Python calibration, US-1976) cannot be its cause.

Per station x optical module (serlom era) x calendar year (>= 60 nights spanning >= 300 d, >= 8 months):
  ln C = a + b (t - t_mid) + annual harmonic  -> amplitude, peak DOY (trend absorbs ageing within the year)
for raw ln C, kappa-normalised ln C (kappa = per-era M2, 31-day anomalies, on daily T_int), T_int, p_calc,
base (dark). Plus, per era over all years: the solstice-mirror test (spring DOY 60-141 minus autumn
DOY 203-284 of the de-trended series) and the network-coherent month anomalies (to compare with the
2025-2026 same-month jumps of q16b).

Outputs: q16g_archive_years.csv, q16g_archive_eras.csv, q16g_archive_months.csv, q16g_archive.png
"""
from __future__ import annotations

import csv
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from q1_temperature import ARCHIVE, robust_mask, running_median  # noqa: E402

Q16 = Path("C:/DATA/Projects/202610_CHM15k_stability/q16")
W = 2 * np.pi / 365.25
LAT = {"ber": 46.95, "cbw": 51.97, "dav": 46.81, "eri": 47.08, "gal": 45.44, "ino": 44.35, "joy": 50.91,
       "kse": 46.55, "lei": 51.35, "lin": 52.21, "mhd": 53.33, "pay": 46.81, "sir": 48.71}


def load_station(f):
    t, c, T, ser, pc, base = [], [], [], [], [], []
    for r in csv.DictReader(open(f, encoding="utf-8")):
        try:
            d = datetime.strptime(r["Time"], "%d-%b-%Y")
        except ValueError:
            continue

        def num(k):
            v = r.get(k, "")
            return float(v) if v not in ("", "NaN") else np.nan
        t.append((d - datetime(1970, 1, 1)).days + 0.5)
        c.append(num("Lidar_constant"))
        T.append(num("temp_int") - 273.15)
        ser.append(r["serlom"])
        pc.append(num("p_calc"))
        base.append(num("base"))
    return {k: np.array(v) for k, v in dict(t=t, c=c, T=T, ser=ser, pc=pc, base=base).items()}


def doy_of(t):
    return np.array([datetime.utcfromtimestamp(x * 86400).timetuple().tm_yday for x in t])


def harmonic_trend(t, y, trend=True):
    """Offset (+ linear trend) + annual harmonic. Returns amp (same units as y), peak DOY, coef, residual."""
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    ph = W * (doy_of(t) - 1.0)
    cols = [np.ones_like(t), np.sin(ph), np.cos(ph)]
    if trend:
        cols.append((t - t.mean()) / 365.25)
    X = np.column_stack(cols)
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    amp = float(np.hypot(c[1], c[2]))
    peak = float((np.arctan2(c[1], c[2]) / W) % 365.25 + 1.0)
    return amp, peak, c, y - X @ c


def kappa_m2(t, lnc, T):
    """Multi-day kappa on 31-day running-median anomalies (report 19 M2), in 1/K."""
    la = lnc - running_median(t, lnc, 31.0)
    ta = T - running_median(t, T, 31.0)
    ok = np.isfinite(la) & np.isfinite(ta)
    if ok.sum() < 60:
        return np.nan
    X = np.column_stack([np.ones(ok.sum()), ta[ok]])
    c, *_ = np.linalg.lstsq(X, la[ok], rcond=None)
    return float(c[1])


def circ_mean_doy(peaks, weights=None):
    a = W * (np.asarray(peaks) - 1.0)
    w = np.ones_like(a) if weights is None else np.asarray(weights)
    z = np.sum(w * np.exp(1j * a)) / np.sum(w)
    return float((np.angle(z) / W) % 365.25 + 1.0), float(np.abs(z))


def main():
    years_rows, era_rows, series = [], [], []
    for f in sorted(ARCHIVE.glob("chm_*_1d.csv")):
        stn = f.stem.split("_")[1]
        d = load_station(f)
        good = (d["T"] > -23) & (d["c"] > 0) & (d["c"] < 1e14)
        for mod in np.unique(d["ser"][good]):
            m = (d["ser"] == mod) & good
            t, lnc, T, pc, base = d["t"][m], np.log(d["c"][m]), d["T"][m], d["pc"][m], d["base"][m]
            keep = robust_mask(lnc)            # per module, never per station
            t, lnc, T, pc, base = t[keep], lnc[keep], T[keep], pc[keep], base[keep]
            o = np.argsort(t)
            t, lnc, T, pc, base = t[o], lnc[o], T[o], pc[o], base[o]
            if len(t) < 60 or t[-1] - t[0] < 300:
                continue
            k = kappa_m2(t, lnc, T)
            lnk = lnc - (k if np.isfinite(k) else 0.0) * T
            # era-wide fit: harmonic + trend
            amp, pk, _, res = harmonic_trend(t, lnc)
            ampk, pkk, _, resk = harmonic_trend(t, lnk)
            ampT, pkT, _, _ = harmonic_trend(t, T)
            # solstice-mirror test on the de-trended series (harmonic removed? no: trend only)
            _, _, ctr, _ = harmonic_trend(t, lnc)
            detr = lnc - ctr[0] - ctr[3] * (t - t.mean()) / 365.25
            detrk = lnk - harmonic_trend(t, lnk)[2][0] - harmonic_trend(t, lnk)[2][3] * (t - t.mean()) / 365.25
            dy = doy_of(t)
            spr = (dy >= 60) & (dy <= 141)
            aut = (dy >= 203) & (dy <= 284)
            Td = T
            mirror = 100 * (np.mean(detr[spr]) - np.mean(detr[aut])) if spr.sum() > 10 and aut.sum() > 10 else np.nan
            mirrork = 100 * (np.mean(detrk[spr]) - np.mean(detrk[aut])) if spr.sum() > 10 and aut.sum() > 10 else np.nan
            mirrorT = (np.mean(Td[spr]) - np.mean(Td[aut])) if spr.sum() > 10 and aut.sum() > 10 else np.nan
            yrs = np.array([datetime.utcfromtimestamp(x * 86400).year for x in t])
            era_rows.append(dict(station=stn, lat=LAT.get(stn, np.nan), module=mod, n=len(t),
                                 start=f"{datetime.utcfromtimestamp(t[0] * 86400):%Y-%m-%d}",
                                 end=f"{datetime.utcfromtimestamp(t[-1] * 86400):%Y-%m-%d}",
                                 n_years=len(np.unique(yrs)), kappa_m2=100 * k,
                                 amp=100 * amp, peak=pk, amp_k=100 * ampk, peak_k=pkk, amp_T=ampT, peak_T=pkT,
                                 mirror_spring_minus_autumn=mirror, mirror_k=mirrork, mirror_T=mirrorT,
                                 scatter=100 * float(np.std(res))))
            series.append(dict(station=stn, module=mod, t=t, lnc=lnc, lnk=lnk, T=T, res=res))
            for y in np.unique(yrs):
                s = yrs == y
                months = len(np.unique([datetime.utcfromtimestamp(x * 86400).month for x in t[s]]))
                if s.sum() < 60 or t[s].max() - t[s].min() < 300 or months < 8:
                    continue
                r = dict(station=stn, lat=LAT.get(stn, np.nan), module=mod, year=int(y), n=int(s.sum()),
                         months=months, kappa_m2=100 * k)
                for tag, arr, scale in (("lnc", lnc, 100), ("lnk", lnk, 100), ("T", T, 1),
                                        ("lnpcal", np.log(np.where(pc > 0, pc, np.nan)), 100),
                                        ("lnbase", np.log(np.where(base > 0, base, np.nan)), 100)):
                    if np.isfinite(arr[s]).sum() < 60:
                        continue
                    a, p, _, rr = harmonic_trend(t[s], arr[s])
                    r[f"amp_{tag}"] = scale * a
                    r[f"peak_{tag}"] = p
                    if tag == "lnc":
                        r["scatter_lnc"] = 100 * float(np.std(rr))
                        # in-phase amplitude with the solstice (DOY 172) - signed
                        r["inphase172_lnc"] = 100 * a * np.cos(W * (p - 172.0))
                years_rows.append(r)
        print(f"{stn}: {sum(1 for r in years_rows if r['station'] == stn)} station-years", flush=True)

    # ---------------------------------------------------------------- network-coherent month anomalies
    # residual of each era (harmonic + trend removed) averaged per calendar month; then, per (year, month),
    # the median across stations and the fraction of stations sharing its sign
    mrows = []
    acc = {}
    for s in series:
        ym = np.array([(datetime.utcfromtimestamp(x * 86400).year, datetime.utcfromtimestamp(x * 86400).month)
                       for x in s["t"]])
        for key in set(map(tuple, ym)):
            sel = (ym[:, 0] == key[0]) & (ym[:, 1] == key[1])
            if sel.sum() >= 5:
                acc.setdefault(key, []).append(100 * float(np.median(s["res"][sel])))
    for key in sorted(acc):
        v = np.array(acc[key])
        if len(v) >= 4:
            med = float(np.median(v))
            mrows.append(dict(year=key[0], month=key[1], n_stations=len(v), median_resid=med,
                              frac_same_sign=float(np.mean(np.sign(v) == np.sign(med))),
                              iqr=float(np.subtract(*np.percentile(v, [75, 25])))))
    write_csv(years_rows, Q16 / "q16g_archive_years.csv")
    write_csv(era_rows, Q16 / "q16g_archive_eras.csv")
    write_csv(mrows, Q16 / "q16g_archive_months.csv")
    summarise(years_rows, era_rows, mrows)
    figure(years_rows, era_rows, series, mrows)


def summarise(years_rows, era_rows, mrows):
    yr = [r for r in years_rows if "amp_lnc" in r]
    yr_np = [r for r in yr if r["station"] != "pay" or r["module"] != "140016"]
    amp = np.array([r["amp_lnc"] for r in yr_np])
    pk = np.array([r["peak_lnc"] for r in yr_np])
    ampk = np.array([r.get("amp_lnk", np.nan) for r in yr_np])
    pkk = np.array([r.get("peak_lnk", np.nan) for r in yr_np])
    pT = np.array([r.get("peak_T", np.nan) for r in yr_np])
    ip = np.array([r["inphase172_lnc"] for r in yr_np])
    lines = [f"station-years (excl. Payerne TUB140016): {len(yr_np)}",
             f"raw ln C: amp median {np.median(amp):.1f} % IQR [{np.percentile(amp, 25):.1f}, {np.percentile(amp, 75):.1f}],"
             f" circular mean peak {circ_mean_doy(pk)[0]:.0f} (R={circ_mean_doy(pk)[1]:.2f}),"
             f" peaks within DOY 120-230: {np.mean((pk > 120) & (pk < 230)):.0%}, in-phase(172) > 0: {np.mean(ip > 0):.0%}",
             f"kappa-normalised: amp median {np.nanmedian(ampk):.1f} %, peak {circ_mean_doy(pkk[np.isfinite(pkk)])[0]:.0f}"
             f" (R={circ_mean_doy(pkk[np.isfinite(pkk)])[1]:.2f})",
             f"T_int (daily) peak {circ_mean_doy(pT[np.isfinite(pT)])[0]:.0f}"]
    for y in sorted({r["year"] for r in yr_np}):
        rr = [r for r in yr_np if r["year"] == y]
        a = np.array([r["amp_lnc"] for r in rr])
        p = np.array([r["peak_lnc"] for r in rr])
        lines.append(f"  {y}: n={len(rr)} amp {np.median(a):.1f} % peak {circ_mean_doy(p)[0]:.0f} (R={circ_mean_doy(p)[1]:.2f})")
    era = [r for r in era_rows if r["n_years"] >= 2]
    mir = np.array([r["mirror_spring_minus_autumn"] for r in era if not (r["station"] == "pay" and r["module"] == "140016")])
    mT = np.array([r["mirror_T"] for r in era if not (r["station"] == "pay" and r["module"] == "140016")])
    lines.append(f"mirror test (eras >= 2 yr, n={len(mir)}): ln C spring-autumn {np.nanmedian(mir):+.1f} % "
                 f"(mean {np.nanmean(mir):+.1f} +- {np.nanstd(mir) / np.sqrt(np.isfinite(mir).sum()):.1f}); "
                 f"T_int {np.nanmedian(mT):+.1f} K")
    mm = np.array([r["median_resid"] for r in mrows])
    big = [r for r in mrows if abs(r["median_resid"]) >= 5 and r["frac_same_sign"] >= 0.75]
    lines.append(f"network-coherent month anomalies: sd of the station-median {np.std(mm):.1f} %; "
                 f"|med| >= 5 % with >= 75 % same sign: {len(big)} of {len(mrows)} months")
    for r in big:
        lines.append(f"   {r['year']}-{r['month']:02d}: {r['median_resid']:+.1f} % ({r['n_stations']} stations, "
                     f"{r['frac_same_sign']:.0%} same sign)")
    txt = "\n".join(lines)
    print(txt)
    (Q16 / "q16g_archive_summary.txt").write_text(txt, encoding="utf-8")


def figure(years_rows, era_rows, series, mrows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    yr = [r for r in years_rows if "amp_lnc" in r and not (r["station"] == "pay" and r["module"] == "140016")]
    stations = sorted({r["station"] for r in yr})
    cm = plt.get_cmap("tab20")
    for i, s in enumerate(stations):
        rr = [r for r in yr if r["station"] == s]
        ax[0].scatter([r["year"] + 0.04 * (i - 6) for r in rr], [r["peak_lnc"] for r in rr],
                      s=[6 * r["amp_lnc"] for r in rr], color=cm(i), label=s.upper(), alpha=0.8)
    ax[0].axhline(172, color="k", ls="--", lw=0.8)
    ax[0].set_ylim(0, 366)
    ax[0].set_xlabel("year")
    ax[0].set_ylabel("peak DOY of the raw ln C_L harmonic (marker size ~ amplitude)")
    ax[0].set_title("(a) per station-year: phase of the annual cycle")
    ax[0].legend(fontsize=7, ncol=2, loc="lower left")
    for i, s in enumerate(stations):
        rr = [r for r in yr if r["station"] == s]
        ax[1].plot([r["year"] for r in rr], [r["amp_lnc"] for r in rr], "o-", color=cm(i), ms=4, lw=0.8)
    ax[1].set_xlabel("year")
    ax[1].set_ylabel("annual amplitude of ln C_L [%]")
    ax[1].set_title("(b) amplitude per station-year (raw, MATLAB chain on raw Lufft)")
    # composite: per era residual + harmonic vs DOY, binned 15 d, median across eras
    bins = np.arange(1, 367, 15)
    comp = []
    for s in series:
        if s["station"] == "pay" and s["module"] == "140016":
            continue
        _, _, c, _ = harmonic_trend(s["t"], s["lnc"])
        detr = s["lnc"] - c[0] - c[3] * (s["t"] - s["t"].mean()) / 365.25
        dy = doy_of(s["t"])
        comp.append([np.nanmedian(detr[(dy >= b) & (dy < b + 15)]) if ((dy >= b) & (dy < b + 15)).sum() > 3 else np.nan
                     for b in bins])
        _, _, cT, _ = harmonic_trend(s["t"], s["T"])
    comp = 100 * np.array(comp)
    ax[2].plot(bins + 7, comp.T, color="0.8", lw=0.6)
    ax[2].plot(bins + 7, np.nanmedian(comp, axis=0), "k-o", ms=4, label="median over module eras")
    ax[2].axvline(172, color="r", ls="--", lw=0.8, label="solstice")
    ax[2].set_xlabel("day of year")
    ax[2].set_ylabel("de-trended ln C_L [%]")
    ax[2].set_title("(c) composite seasonal course, all years (2013-2024)")
    ax[2].legend(fontsize=8)
    fig.suptitle("q16g - 11-year archive (raw Lufft beta_raw + MATLAB/NCEP Rayleigh, independent of raw2l1 and the Python chain)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(Q16 / "q16g_archive.png", dpi=130)


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


if __name__ == "__main__":
    main()
