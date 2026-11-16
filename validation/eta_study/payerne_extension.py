# -*- coding: utf-8 -*-
"""Payerne CL61 vs CHM15k, June-October 2026: the inter-comparison dashboard arithmetic replayed
on an L1-only hourly cache that extends past the dashboard window (2026-06-11..08-14), with the
calibration series of the dark-corrected reruns (hood dark v3, 10 sessions) and the reruns
of this study for the period after 2026-08-13.

Why: the "+7.8 %" of the handoff (CL61 Rayleigh+dark vs CHM15k v2.2+dark, page default) is the
median over a window in which the gap drifts from +19 % (June) to -3 % (August); the autumn months
(lower molecular windows, less dark leverage, cleaner CHM15k nights) decide whether the summer level
or the August level is the representative one.

Inputs
  cache         C:/DATA/Projects/202611_eta_study/ext/_streams_payerne_ext.npz
                (build_l1_l2_dashboard.build_stream_cache, SOURCES=("L1",), 2026-06-01..2026-10-31)
  dark profile  rayleigh_availability/dark_profiles_payerne_v3.npz (the profile the v3 reruns subtract)
  CL61 Rayleigh+dark   diag_v22_dark_v3/0-20000-0-06610_C         (eprof_v2.2, 910.55/0.188, dark v3)
  CL61 Rayleigh, cloud diag_v22_l55s008 (to 2026-08-13) + ext/c_nodark (from 2026-08-14), no dark
  CHM15k v2.2+dark     diag_v22_dark_v3/0-20000-0-06610_A
  CHM15k v2.2          calib_raw/payerne_A_v2.2.json (to 2026-08-14) + ext/a_v22 (from 2026-08-14)
Statistics: exactly the page's medrel (median of the per-pixel relative difference over the
500-3000 m band), per month, night (sun < -10 deg) / day (sun > 0), and per altitude band; the
CL61 profile side uses the manufacturer line (the page's "ctor" state) and the molecular-aware
910 -> 1064 nm conversion with alpha = 1.

Outputs (C:/DATA/Projects/202611_eta_study/ext/):
  monthly_gap.csv     medrel per (cur, ref, month, subset)
  raw_ratio_monthly.csv   uncalibrated signal ratios CHM15k/CL61 and CL31/CL61 per month
  nights_C.csv, nights_A.csv   per-night constants with the molecular window of each series
  summary.txt

Run:  python validation/eta_study/payerne_extension.py [--dark v1|v3]   (v3 default; v1 = the dashboard's 6-session profile)
"""
from __future__ import annotations

import csv
import json
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "inter-comparison_dashboard"))

from validation.paper import intercompare as IC            # noqa: E402
from monitoring.kalman import kalman_best_estimate          # noqa: E402
import build_l1_l2_dashboard as BD                          # noqa: E402
import build_v3 as B3                                       # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
EXT = Path("C:/DATA/Projects/202611_eta_study/ext")
CACHE = EXT / "_streams_payerne_ext.npz"
DARK_VERSION = "v1" if "--dark" in sys.argv and sys.argv[sys.argv.index("--dark") + 1] == "v1" else "v3"
# v3 = 10 hood sessions May-Sep 2026, detected windows, P-smoothed (build_dark_v3.py, the paper reruns);
# v1 = the 6-session logbook-window profile the inter-comparison dashboard was built with (2026-08).
DARK_NPZ = DATA / "rayleigh_availability" / ("dark_profiles_payerne_v3.npz" if DARK_VERSION == "v3" else "dark_profiles_payerne.npz")
RUN_C_DARK = DATA / ("diag_v22_dark_v3" if DARK_VERSION == "v3" else "diag_v22_l55s008dark")
RUN_A_DARK = DATA / ("diag_v22_dark_v3" if DARK_VERSION == "v3" else "diag_v22_dark")
SUFFIX = "" if DARK_VERSION == "v3" else "_darkv1"
WMO = "0-20000-0-06610"
TARGET = 1064.0
MONTHS = ["2026-06", "2026-07", "2026-08", "2026-09", "2026-10"]
SUCCESS = ("1", "1.0", "0.5")


# ----------------------------------------------------------------------------- calibration series
def csv_nights(path, method):
    out = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["method"] != method or r["flag"] not in SUCCESS:
                continue
            try:
                v = float(r["cal_value"])
            except ValueError:
                continue
            if v > 0:
                out[r["date"]] = (v, float(r["uncertainty"] or "nan"),
                                  float(r["bottom_height"] or "nan"), float(r["top_height"] or "nan"))
    return out


def json_nights(path):
    rec = json.loads(Path(path).read_text())
    out = {}
    for day, v in rec.items():
        if v[0] in (1.0, 0.5) and v[1] and v[1] > 0:
            out[day] = (float(v[1]), float(v[2]) if v[2] else np.nan,
                        float(v[3]) if v[3] else np.nan, float(v[4]) if v[4] else np.nan)
    return out


def merged(*parts):
    """Later parts override earlier ones on common dates."""
    out = {}
    for p in parts:
        out.update(p)
    return dict(sorted(out.items()))


def kalman_on(hours, nights):
    dates = [datetime.strptime(d, "%Y%m%d") for d in nights]
    C = np.array([v[0] for v in nights.values()])
    U = np.array([v[1] for v in nights.values()])
    kt, ks, _ = kalman_best_estimate(dates, C, uncertainties=U)
    return IC.interp_calib(kt, ks, hours)


def series():
    c_dark = csv_nights(RUN_C_DARK / f"{WMO}_C" / f"{WMO}_C_cal.csv", "rayleigh")
    c_nodk = merged(csv_nights(DATA / "diag_v22_l55s008" / f"{WMO}_C" / f"{WMO}_C_cal.csv", "rayleigh"),
                    csv_nights(EXT / "c_nodark" / f"{WMO}_C" / f"{WMO}_C_cal.csv", "rayleigh"))
    c_cld = merged(csv_nights(DATA / "diag_v22_l55s008" / f"{WMO}_C" / f"{WMO}_C_cal.csv", "cloud"),
                   csv_nights(EXT / "c_nodark" / f"{WMO}_C" / f"{WMO}_C_cal.csv", "cloud"))
    a_dark = csv_nights(RUN_A_DARK / f"{WMO}_A" / f"{WMO}_A_cal.csv", "rayleigh")
    a_v22 = merged(json_nights(DATA / "inter-comparison_dashboard" / "calib_raw" / "payerne_A_v2.2.json"),
                   csv_nights(EXT / "a_v22" / f"{WMO}_A" / f"{WMO}_A_cal.csv", "rayleigh"))
    return dict(C_raydark=c_dark, C_raynodark=c_nodk, C_cloud=c_cld, A_v22dark=a_dark, A_v22=a_v22)


# ----------------------------------------------------------------------------- sun
def sun_elev(t, lat_deg, lon_deg):
    lat = np.deg2rad(lat_deg)
    t = t.astype("datetime64[m]")
    doy = (t.astype("datetime64[D]") - t.astype("datetime64[Y]")).astype(int) + 1
    hr = (t - t.astype("datetime64[D]")).astype("timedelta64[m]").astype(float) / 60.0 + 0.5
    g = 2 * np.pi / 365.0 * (doy - 1 + (hr - 12) / 24)
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g)
                    - 0.040849 * np.sin(2 * g))
    ha = np.deg2rad((hr * 60 + eqt + 4 * lon_deg) / 4 - 180)
    return np.rad2deg(np.arcsin(np.sin(lat) * np.sin(decl) + np.cos(lat) * np.cos(decl) * np.cos(ha)))


# ----------------------------------------------------------------------------- main
def main():
    BD.set_site("payerne")
    BD.SOURCES = ("L1",)
    s = {ch["ident"]: ch["itype"] for ch in BD.CHANNELS}
    station_alt = float(BD.STA["alt"]) if isinstance(BD.STA, dict) and "alt" in BD.STA else 491.0
    raw = BD._read_cache(CACHE)
    z_agl = B3.Z_AGL
    alt_asl = z_agl + station_alt
    band = np.where(B3.BAND)[0]
    blocks, hours_of = {}, {}
    for ident in ("A", "B", "C"):
        d = dict(raw[f"L1_{ident}"])
        d["itype"] = s[ident]
        g, S, _, _ = B3.stream_block(d, alt_asl)
        blocks[ident] = S
        hours_of[ident] = g
        raw[f"L1_{ident}"] = d
    union = np.unique(np.concatenate(list(hours_of.values())))
    al = {}
    for ident, S in blocks.items():
        A = np.full((union.size, z_agl.size), np.nan)
        A[np.searchsorted(union, hours_of[ident])] = S
        al[ident] = A
    have = np.any(np.isfinite(al["A"][:, band]), axis=1) & np.any(np.isfinite(al["C"][:, band]), axis=1)
    hours = union[have]
    print(f"paired hours A&C: {hours.size} over {str(hours[0])[:10]}..{str(hours[-1])[:10]}", flush=True)
    days = [str(x).replace("-", "") for x in np.unique(hours.astype("datetime64[D]"))]
    day_index = {d: i for i, d in enumerate(days)}
    hday = np.array([day_index[str(x).replace("-", "")] for x in hours.astype("datetime64[D]")])
    hmon = np.array([str(x) for x in hours.astype("datetime64[M]")])
    lat, lon = float(raw["L1_C"]["lat"]), float(raw["L1_C"]["lon"])
    elev = sun_elev(hours, lat, lon)
    night, day = elev < -10, elev > 0

    # per-day correction ingredients (WV at both lines, molecular at 910.55 and 1064)
    corr = {ident: B3.corr_arrays(raw[f"L1_{ident}"], station_alt, alt_asl, days) for ident in ("B", "C")}
    print("corr C: lam", corr["C"]["lam"], "wv days with CAMS:",
          int(np.isfinite(corr["C"]["wv2"]).any(axis=1).sum()), "/", len(days), flush=True)

    # dark profile v3 on the display grid (rcs_0 units)
    dk = np.load(DARK_NPZ)
    dark = {}
    for ident in ("A", "B", "C"):
        r, b = np.asarray(dk[f"{ident}_range"], "f8"), np.asarray(dk[f"{ident}_b_rcs"], "f8")
        ok = np.isfinite(r) & np.isfinite(b)
        dark[ident] = np.nan_to_num(np.interp(z_agl, r[ok], b[ok], left=np.nan, right=np.nan))

    ser = series()
    cal = {k: kalman_on(hours, v) for k, v in ser.items()}

    def channel(ident, calkey=None, dark_on=False, wvkey="wv2", wl="molecular", alpha=1.0):
        v = al[ident][have].copy()
        if dark_on:
            v = v - dark[ident][None, :]
        if calkey is not None:
            v = v * 1e6 / cal[calkey][:, None]
        c = corr.get(ident)
        if c is None:
            return v
        if wvkey and c.get(wvkey) is not None:
            v = v / c[wvkey][hday]
        if wl == "molecular" and c.get("bmt") is not None:
            f = (c["lam"] / TARGET) ** alpha
            v = c["bmt"][hday] + (v - c["bml"][hday]) * f
        return v

    def medrel(A, Bm, rows, gates=band):
        a, b = A[np.ix_(rows, gates)], Bm[np.ix_(rows, gates)]
        m = np.isfinite(a) & np.isfinite(b) & (b > 0)
        return (float(np.median((a[m] - b[m]) / b[m] * 100)) if m.sum() > 20 else np.nan), int(m.sum())

    curs = {
        "C Rayleigh+dark (v3)": channel("C", "C_raydark", dark_on=True),
        "C Rayleigh no dark": channel("C", "C_raynodark"),
        "C cloud (oper eta)": channel("C", "C_cloud", dark_on=True),
        "C cloud (proposed eta, x0.90)": channel("C", "C_cloud", dark_on=True) / 0.90,
    }
    refs = {
        "A v2.2+dark (v3)": channel("A", "A_v22dark", dark_on=True),
        "A v2.2": channel("A", "A_v22"),
    }
    bands = {"0.5-1 km": band[(z_agl[band] >= 500) & (z_agl[band] < 1000)],
             "1-2 km": band[(z_agl[band] >= 1000) & (z_agl[band] < 2000)],
             "2-3 km": band[(z_agl[band] >= 2000) & (z_agl[band] <= 3000)]}
    rows_out = []
    lines = []
    for rk, R in refs.items():
        for ck, Cc in curs.items():
            for m in MONTHS:
                r = np.flatnonzero(hmon == m)
                for sub, rr in (("all", r), ("night", r[night[r]]), ("day", r[day[r]])):
                    v, n = medrel(Cc, R, rr)
                    rows_out.append(dict(cur=ck, ref=rk, month=m, subset=sub, band="0.5-3 km", medrel=v, n=n, hours=rr.size))
                for bk, g in bands.items():
                    v, n = medrel(Cc, R, r, g)
                    rows_out.append(dict(cur=ck, ref=rk, month=m, subset="all", band=bk, medrel=v, n=n, hours=r.size))
            for sub, rr in (("all", np.arange(hours.size)),):
                v, n = medrel(Cc, R, rr)
                rows_out.append(dict(cur=ck, ref=rk, month="Jun-Oct", subset=sub, band="0.5-3 km", medrel=v, n=n, hours=rr.size))
    with open(EXT / ("monthly_gap" + SUFFIX + ".csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        w.writerows(rows_out)
    for rk in refs:
        lines.append(f"\n== medrel % vs {rk} (band 0.5-3 km): month: all | night | day ==")
        for ck in curs:
            parts = []
            for m in MONTHS:
                vals = {r_["subset"]: r_ for r_ in rows_out if r_["cur"] == ck and r_["ref"] == rk and r_["month"] == m and r_["band"] == "0.5-3 km"}
                parts.append(f"{m[5:]}: {vals['all']['medrel']:6.1f} | {vals['night']['medrel']:6.1f} | {vals['day']['medrel']:6.1f} ({vals['all']['hours']} h)")
            lines.append(f"  {ck:30s} " + "   ".join(parts))
        lines.append("  per band (all hours): " + " ".join(
            f"{m[5:]}: " + "/".join(f"{[r_['medrel'] for r_ in rows_out if r_['cur'] == 'C Rayleigh+dark (v3)' and r_['ref'] == rk and r_['month'] == m and r_['band'] == bk][0]:5.1f}" for bk in bands)
            for m in MONTHS) + "   [Rayleigh+dark, 0.5-1/1-2/2-3 km]")

    # uncalibrated signal ratios
    rawA = al["A"][have] - dark["A"][None, :]
    rawC = al["C"][have] - dark["C"][None, :]
    rawB = al["B"][have] - dark["B"][None, :]

    def q(X, Y, rr, g):
        a, b = X[np.ix_(rr, g)], Y[np.ix_(rr, g)]
        m = np.isfinite(a) & np.isfinite(b) & (b > 0) & (a > 0)
        return float(np.median(a[m] / b[m])) if m.sum() > 20 else np.nan

    g12 = bands["1-2 km"]
    g23 = bands["2-3 km"]
    rr_out = []
    lines.append("\n== uncalibrated signal ratios (hood dark v3 subtracted): A/C x1e-11, B/C x1e-8 ==")
    for m in MONTHS:
        r = np.flatnonzero(hmon == m)
        rn, rd = r[night[r]], r[day[r]]
        rec = dict(month=m, AC_12_night=q(rawA, rawC, rn, g12) * 1e-11, AC_12_day=q(rawA, rawC, rd, g12) * 1e-11,
                   AC_23_night=q(rawA, rawC, rn, g23) * 1e-11, AC_23_day=q(rawA, rawC, rd, g23) * 1e-11,
                   BC_12_night=q(rawB, rawC, rn, g12) * 1e-8, BC_12_day=q(rawB, rawC, rd, g12) * 1e-8,
                   hours=r.size, night_hours=rn.size)
        rr_out.append(rec)
        lines.append(f"  {m}: A/C 1-2 km night {rec['AC_12_night']:.3f} day {rec['AC_12_day']:.3f} | 2-3 km night {rec['AC_23_night']:.3f} "
                     f"day {rec['AC_23_day']:.3f} | B/C 1-2 km night {rec['BC_12_night']:.3f} day {rec['BC_12_day']:.3f}  ({r.size} h, {rn.size} night)")
    with open(EXT / ("raw_ratio_monthly" + SUFFIX + ".csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rr_out[0].keys()))
        w.writeheader()
        w.writerows(rr_out)

    # per-night constants and their monthly medians
    for key, fn in (("C", f"nights_C{SUFFIX}.csv"), ("A", f"nights_A{SUFFIX}.csv")):
        keys = [k for k in ser if k.startswith(key)]
        alld = sorted(set().union(*[set(ser[k]) for k in keys]))
        with open(EXT / fn, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["date"] + [f"{k}_{c}" for k in keys for c in ("C", "unc", "bot", "top")])
            for d in alld:
                row = [d]
                for k in keys:
                    v = ser[k].get(d)
                    row += list(v) if v else [np.nan] * 4
                w.writerow(row)
        lines.append(f"\n== monthly medians of the per-night constants ({key}) ==")
        for m in ["2026-02", "2026-03", "2026-04", "2026-05"] + MONTHS:
            mm = m.replace("-", "")
            parts = []
            for k in keys:
                vals = [v[0] for d, v in ser[k].items() if d.startswith(mm)]
                zc = [0.5 * (v[2] + v[3]) for d, v in ser[k].items() if d.startswith(mm) and np.isfinite(v[2])]
                if vals:
                    sc = 1e-11 if key == "A" else 1.0
                    parts.append(f"{k} {np.median(vals) * sc:6.3f} (n{len(vals):2d}, zc {np.median(zc) if zc else np.nan:4.0f} m)")
            lines.append(f"  {m}: " + " | ".join(parts))
    # Kalman monthly medians
    lines.append("\n== Kalman series interpolated on the paired hours, monthly medians ==")
    for k in cal:
        sc = 1e-11 if k.startswith("A") else 1.0
        lines.append(f"  {k:12s} " + "  ".join(f"{m[5:]}: {np.median(cal[k][hmon == m]) * sc:6.3f}" for m in MONTHS if (hmon == m).any()))
    txt = "\n".join(lines)
    print(txt)
    (EXT / ("summary" + SUFFIX + ".txt")).write_text(txt, encoding="utf-8")
    np.savez(EXT / ("ext_hours" + SUFFIX + ".npz"), hours=hours.astype("datetime64[h]").astype("i8"), elev=elev,
             **{f"cal_{k}": v for k, v in cal.items()})


if __name__ == "__main__":
    main()
