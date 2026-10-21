# -*- coding: utf-8 -*-
"""q16g (task 2): is ABSOLUTE HUMIDITY the driver of the residual CHM15k annual gain cycle?

Hypothesis: water adsorbed by a porous (evaporated) interference-filter coating, or humidity inside the
housing, shifts the filter pass band / transmission. Absolute humidity peaks in July-August, closer to the
solstice than the night T_int (~DOY 206)?

Stage 1 (extraction, 4 workers over CAMS files): CAMS 0.4 deg L137 (D:/CAMS_Monthly_04 2025-01..2026-06,
then D:/CAMS_daily), 3-hourly, at the model level closest to the station altitude, nearest grid cell:
T, q, p -> vapour pressure e, absolute humidity AH = e/(R_v T) [g m-3], relative humidity RH (Magnus, water).
Per station -> q16g_cams_surface.npz.

Stage 2 (per calibration night, all units of q15f_nightly): AH and e over the calibration window (+-1.5 h),
AH 24 h before the window end, exponentially weighted AH memories (3/7/21 d), ambient RH, and the
INSIDE-housing RH proxy RH_in = e / e_s(T_int) (housing air exchanges water with the outside but is warmer).

Tests (as q16b): phase vs the residual; solstice-mirror (spring DOY 60-141 minus autumn 203-284);
synoptic partial gain (ln C anomaly on T_int + driver anomalies, 31-day running-median anomalies) and its
transfer to the annual harmonic; the same for the aerosol-free clear-sky X(7-8 km) (q15e series) and the
window X; latitude / climate test (unit AH amplitude vs residual amplitude); interannual same-month
2026-2025 differences.

Outputs: q16g_humidity_nights.csv, q16g_humidity_units.csv, q16g_humidity_drivers.csv,
q16g_humidity_mirror.csv, q16g_humidity_interannual.csv, q16g_humidity.png
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
from common import CENSUS, load_constants  # noqa: E402
from q16b_insolation import (PAYERNE, Q15, Q16, W, anomaly, circ_median_doy, harmonic, load_kappas, ols,  # noqa: E402
                             peak_doy, residuals, wrap_doy)

N_WORKERS = 4
N_LEV = 30                    # lowest 30 model levels (~2 km above the model ground)
RV = 461.5
EPS = 0.622
EPOCH = datetime(1970, 1, 1)
CAMS_NPZ = Q16 / "q16g_cams_surface.npz"
DRIVERS = ("AH_night", "AH_prev24", "AH_ew3d", "AH_ew7d", "AH_ew21d", "e_night", "RH_night", "RH_in_night",
           "T_cams_night")
UNITS = {"AH_night": "g/m3", "AH_prev24": "g/m3", "AH_ew3d": "g/m3", "AH_ew7d": "g/m3", "AH_ew21d": "g/m3",
         "e_night": "hPa", "RH_night": "%", "RH_in_night": "%", "T_cams_night": "degC"}


# ---------------------------------------------------------------------------------------------
# Stage 1: CAMS extraction
# ---------------------------------------------------------------------------------------------
def cams_files():
    out = []
    for y in (2025, 2026):
        for m in range(1, 13):
            f = Path(f"D:/CAMS_Monthly_04/CAMS_Beta_{y}{m:02d}.nc")
            if f.exists():
                out.append(("monthly", f))
    for f in sorted(Path("D:/CAMS_daily").glob("CAMS_Beta_2026[0-9][0-9][0-9][0-9].nc")):
        out.append(("daily", f))
    return out


def stations():
    out = []
    for s in json.load(open(CENSUS, encoding="utf-8")):
        if s["type"].upper().startswith("CHM"):
            out.append((f"{s['wmo']}_{s['ident']}", float(s["lat"]), float(s["lon"]), float(s.get("alt", 0.0) or 0.0)))
    return out


def e_sat(Tc):
    """Saturation vapour pressure over water [Pa] (Magnus, Alduchov & Eskridge 1996)."""
    return 610.94 * np.exp(17.625 * Tc / (Tc + 243.04))


def extract_file(args):
    import netCDF4
    from calibration.water_vapor_correction.water_vapor import _A137, _B137
    kind, f, sts = args
    try:
        with netCDF4.Dataset(f) as d:
            lat = np.asarray(d["latitude"][:], dtype="f8")
            lon = np.asarray(d["longitude"][:], dtype="f8")
            lev = np.asarray(d["level"][:], dtype=int)[-N_LEV:]
            tt = np.asarray(d["time"][:], dtype="f8") / 24.0 + (datetime(1900, 1, 1) - EPOCH).days
            iy = np.array([int(np.argmin(np.abs(lat - s[1]))) for s in sts])
            ix = np.array([int(np.argmin(np.abs(lon - s[2]))) for s in sts])
            ok = (np.abs(lat[iy] - np.array([s[1] for s in sts])) < 0.3) & (np.abs(lon[ix] - np.array([s[2] for s in sts])) < 0.3)
            y0, y1, x0, x1 = iy.min(), iy.max() + 1, ix.min(), ix.max() + 1
            q = np.asarray(d["q"][:, -N_LEV:, y0:y1, x0:x1], dtype="f4")
            T = np.asarray(d["t"][:, -N_LEV:, y0:y1, x0:x1], dtype="f4")
            alt = np.asarray(d["altitude"][:, -N_LEV:, y0:y1, x0:x1], dtype="f4")
            lnsp = np.asarray(d["lnsp"][:, 0, y0:y1, x0:x1], dtype="f8")   # surface field stored at the first level only
    except Exception as e:  # noqa: BLE001
        print(f"  {f.name}: {e}", flush=True)
        return None
    out = {}
    for k, s in enumerate(sts):
        if not ok[k]:
            continue
        j, i = iy[k] - y0, ix[k] - x0
        a = alt[:, :, j, i]                                  # (nt, nlev)
        il = np.argmin(np.abs(a - s[3]), axis=1)             # level closest to the station altitude
        nt = len(tt)
        qq = q[np.arange(nt), il, j, i].astype("f8")
        TT = T[np.arange(nt), il, j, i].astype("f8")
        ps = np.exp(lnsp[:, j, i])
        L = lev[il]
        p = 0.5 * (_A137[L - 1] + _B137[L - 1] * ps + _A137[L] + _B137[L] * ps)
        e = qq * p / (EPS + (1 - EPS) * qq)
        ah = 1000.0 * e / (RV * TT)
        rh = 100.0 * e / e_sat(TT - 273.15)
        out[s[0]] = np.column_stack([tt, ah, e / 100.0, rh, TT - 273.15, a[np.arange(nt), il] - s[3]])
    print(f"  {f.name}: {len(out)} stations", flush=True)
    return kind, out


def stage_extract():
    sts = stations()
    files = cams_files()
    acc = {}
    with ProcessPoolExecutor(N_WORKERS) as ex:
        for res in ex.map(extract_file, [(k, f, sts) for k, f in files]):
            if res is None:
                continue
            kind, out = res
            for s, arr in out.items():
                acc.setdefault(s, {"monthly": [], "daily": []})[kind].append(arr)
    final = {}
    for s, dd in acc.items():
        mo = np.concatenate(dd["monthly"]) if dd["monthly"] else np.zeros((0, 6))
        da = np.concatenate(dd["daily"]) if dd["daily"] else np.zeros((0, 6))
        tmax = mo[:, 0].max() if len(mo) else -np.inf
        da = da[da[:, 0] > tmax + 1e-6]
        arr = np.concatenate([mo, da])
        arr = arr[np.argsort(arr[:, 0])]
        _, u = np.unique(np.round(arr[:, 0] * 8), return_index=True)
        final[s] = arr[u]
    np.savez_compressed(CAMS_NPZ, **final)
    print(f"-> {CAMS_NPZ}: {len(final)} stations", flush=True)


# ---------------------------------------------------------------------------------------------
# Stage 2: per-night drivers
# ---------------------------------------------------------------------------------------------
def ew_mean(t, x, tau):
    """Causal exponentially weighted mean on the 3-hourly series (gaps ignored)."""
    out = np.full_like(x, np.nan)
    acc, wsum, tprev = 0.0, 0.0, None
    for i in range(len(t)):
        if tprev is not None:
            dec = np.exp(-(t[i] - tprev) / tau)
            acc *= dec
            wsum *= dec
        if np.isfinite(x[i]):
            acc += x[i]
            wsum += 1.0
        out[i] = acc / wsum if wsum > 0 else np.nan
        tprev = t[i]
    return out


def night_drivers(m, r, cams):
    s = m["stream"]
    if s not in cams:
        return None
    c = cams[s]
    tc, ah, e, rh, Tc = c[:, 0], c[:, 1], c[:, 2], c[:, 3], c[:, 4]
    cons = load_constants(*s.rsplit("_", 1))
    # window start/end of every night of the q15f table (matched on the calibration time)
    idx = {round(float(v), 4): i for i, v in enumerate(cons["t"])}
    j = np.array([idx.get(round(float(v), 4), -1) for v in r["t"]])
    t0 = np.where(j >= 0, cons["t0"][np.clip(j, 0, None)], r["t"] - 0.1)
    t1 = np.where(j >= 0, cons["t1"][np.clip(j, 0, None)], r["t"] + 0.1)
    ews = {tau: ew_mean(tc, ah, tau) for tau in (3.0, 7.0, 21.0)}
    n = len(r["t"])
    D = {k: np.full(n, np.nan) for k in DRIVERS}
    for i in range(n):
        w = (tc >= t0[i] - 1.5 / 24) & (tc <= t1[i] + 1.5 / 24)
        if not w.any():
            w = np.abs(tc - r["t"][i]) <= 1.5 / 24
        if w.any():
            D["AH_night"][i] = ah[w].mean()
            D["e_night"][i] = e[w].mean()
            D["RH_night"][i] = rh[w].mean()
            D["T_cams_night"][i] = Tc[w].mean()
            Tin = r["T_int"][i]
            if np.isfinite(Tin):
                D["RH_in_night"][i] = 100.0 * (100.0 * e[w].mean()) / (610.94 * np.exp(17.625 * Tin / (Tin + 243.04)))
        p = (tc > t1[i] - 1.0) & (tc <= t1[i])
        if p.sum() >= 4:
            D["AH_prev24"][i] = ah[p].mean()
        k = np.searchsorted(tc, t0[i]) - 1
        if 0 <= k < len(tc) and t0[i] - tc[k] < 0.5:
            for tau, v in ews.items():
                D[f"AH_ew{int(tau)}d"][i] = v[k]
    return D


def x_series(ser, skeys, s, r, col):
    if s not in skeys:
        return np.full(len(r["t"]), np.nan)
    idx = {round(float(v), 9): i for i, v in enumerate(ser[f"{s}__lnc"])}
    j = np.array([idx.get(round(float(v), 9), -1) for v in r["lnC"]])
    src = ser[f"{s}__lnx"][:, -1] if col == "x78" else ser[f"{s}__lnx_win"]
    return np.where(j >= 0, src[np.clip(j, 0, None)], np.nan)


def stage_analyse():
    units = list(np.load(Q15 / "q15f_nightly.npz", allow_pickle=True)["units"])
    cams = dict(np.load(CAMS_NPZ))
    k_in, k_m2, _ = load_kappas()
    ser = np.load(Q15 / "q15e_series.npz")
    skeys = set(k.split("__")[0] for k in ser.files)
    unit_rows, drows, nights = [], [], []
    store = {}
    for m, r in units:
        s = m["stream"]
        D = night_drivers(m, r, cams)
        if D is None:
            continue
        t = r["t"]
        Y = {"lnC": r["lnC"], "res_in": residuals(r, k_in.get(s, np.nan)),
             "res_m2": residuals(r, k_m2.get((s, m["module"]), np.nan)),
             "x78": x_series(ser, skeys, s, r, "x78"), "xwin": x_series(ser, skeys, s, r, "win")}
        store[s] = (m, r, D, Y)
        for i in range(len(t)):
            nights.append(dict(stream=s, module=m["module"], t=round(float(t[i]), 4), lnC=float(r["lnC"][i]),
                               T_int=float(r["T_int"][i]), **{k: round(float(D[k][i]), 4) for k in DRIVERS}))
        u = dict(stream=s, site=m["site"], lat=m["lat"], lon=m["lon"], alt=m["alt"],
                 kappa_in=k_in.get(s, np.nan), dz_cams=float(np.nanmedian(cams[s][:, 5])))
        ph = {}
        for k, y in Y.items():
            p, _, _ = harmonic(t, y)
            ph[k] = p
            u[f"amp_{k}"] = 100 * abs(p)
            u[f"peak_{k}"] = peak_doy(p)
        pT, _, _ = harmonic(t, r["T_int"])
        u["amp_Tint"], u["peak_Tint"] = abs(pT), peak_doy(pT)
        ta = anomaly(t, r["T_int"])
        for name in DRIVERS:
            x = D[name]
            if np.isfinite(x).sum() < 60:
                continue
            pD, _, _ = harmonic(t, x)
            u[f"amp_{name}"], u[f"peak_{name}"] = abs(pD), peak_doy(pD)
            row = dict(stream=s, driver=name, lat=m["lat"], amp_D=abs(pD), peak_D=peak_doy(pD),
                       mean_D=float(np.nanmean(x)))
            xa = anomaly(t, x)
            for yk in ("lnC", "x78", "xwin"):
                ya = anomaly(t, Y[yk])
                ok = np.isfinite(ya) & np.isfinite(ta) & np.isfinite(xa)
                if ok.sum() < 50 or np.std(xa[ok]) < 1e-9:
                    continue
                sd = np.std(xa[ok])
                c, se, _ = ols(np.column_stack([np.ones(ok.sum()), ta[ok], xa[ok] / sd]), ya[ok])
                e1 = ols(np.column_stack([np.ones(ok.sum()), ta[ok]]), ya[ok])[2]
                e2 = ols(np.column_stack([np.ones(ok.sum()), ta[ok]]), xa[ok])[2]
                row[f"g_sd_{yk}"] = 100 * c[2]
                row[f"g_{yk}"] = 100 * c[2] / sd
                row[f"pcorr_{yk}"] = float(np.corrcoef(e1, e2)[0, 1])
                row["sd_anom"] = sd
                # transfer of the synoptic gain to the annual harmonic, projected on the residual phasor
                for rk in ("res_in", "res_m2"):
                    pr = ph[rk]
                    if np.isfinite(pr):
                        row[f"transfer_{yk}_{rk}"] = 100 * np.real(c[2] / sd * pD * np.conj(pr)) / abs(pr)
            for rk in ("res_in", "res_m2", "x78"):
                pr = ph[rk]
                if np.isfinite(pr) and np.isfinite(pD):
                    dphi = wrap_doy(peak_doy(pD) - peak_doy(pr))
                    row[f"dphase_{rk}"] = dphi
                    # gain needed for the driver's in-phase harmonic to carry the residual
                    row[f"needed_{rk}"] = 100 * abs(pr) / abs(pD)
            drows.append(row)
        unit_rows.append(u)
    write_csv(nights, Q16 / "q16g_humidity_nights.csv")
    write_csv(unit_rows, Q16 / "q16g_humidity_units.csv")
    write_csv(drows, Q16 / "q16g_humidity_drivers.csv")
    summary(unit_rows, drows, store)


def pooled(store, name, yk, n_boot=300, seed=2):
    blocks = []
    for s, (m, r, D, Y) in store.items():
        if s == PAYERNE:
            continue
        t = r["t"]
        ya, ta, xa = anomaly(t, Y[yk]), anomaly(t, r["T_int"]), anomaly(t, D[name])
        ok = np.isfinite(ya) & np.isfinite(ta) & np.isfinite(xa)
        if ok.sum() >= 50:
            blocks.append((ya[ok], ta[ok], xa[ok]))
    if len(blocks) < 10:
        return np.nan, np.nan, np.nan, 0

    def fit(bl):
        y = np.concatenate([b[0] for b in bl])
        X = np.column_stack([np.ones(len(y)), np.concatenate([b[1] for b in bl]), np.concatenate([b[2] for b in bl])])
        return np.linalg.lstsq(X, y, rcond=None)[0]
    c = fit(blocks)
    rng = np.random.default_rng(seed)
    bs = np.array([fit([blocks[i] for i in rng.integers(0, len(blocks), len(blocks))])[2] for _ in range(n_boot)])
    sdx = np.std(np.concatenate([b[2] for b in blocks]))
    return 100 * c[2], 100 * bs.std(), 100 * c[2] * sdx, len(blocks)


def mirror(store, n_boot=1000):
    rng = np.random.default_rng(0)
    rows = []
    fns = {"ln C (raw)": lambda m, r, D, Y: 100 * Y["lnC"], "residual kappa_in": lambda m, r, D, Y: 100 * Y["res_in"],
           "residual kappa_m2": lambda m, r, D, Y: 100 * Y["res_m2"], "ln X 7-8 km": lambda m, r, D, Y: 100 * Y["x78"],
           "T_int [K]": lambda m, r, D, Y: r["T_int"]}
    for name in DRIVERS:
        fns[f"{name} [{UNITS[name]}]"] = (lambda nm: (lambda m, r, D, Y: D[nm]))(name)
    curves = {}
    edges = np.arange(0, 181, 15.0)
    for lab, fn in fns.items():
        diffs = []
        comp = {"spring": [[] for _ in edges[:-1]], "autumn": [[] for _ in edges[:-1]]}
        for s, (m, r, D, Y) in store.items():
            if s == PAYERNE:
                continue
            y = fn(m, r, D, Y)
            ok = np.isfinite(y)
            if ok.sum() < 60:
                continue
            y = y - np.nanmean(y)
            doy = np.mod(r["t"], 365.25)
            a = y[(doy >= 60) & (doy < 141) & ok]
            b = y[(doy >= 203) & (doy < 284) & ok]
            if len(a) >= 8 and len(b) >= 8:
                diffs.append(a.mean() - b.mean())
            dist = np.abs(wrap_doy(doy - 172.0))
            br = np.where(wrap_doy(doy - 172.0) < 0, "spring", "autumn")
            for b_ in ("spring", "autumn"):
                for i in range(len(edges) - 1):
                    sel = ok & (br == b_) & (dist >= edges[i]) & (dist < edges[i + 1])
                    if sel.sum() >= 3:
                        comp[b_][i].append(y[sel].mean())
        d = np.array(diffs)
        if d.size < 10:
            continue
        bs = [np.median(rng.choice(d, d.size)) for _ in range(n_boot)]
        rows.append(dict(series=lab, spring_minus_autumn=float(np.median(d)), se=float(np.std(bs)), n_units=int(d.size)))
        curves[lab] = {b_: np.array([np.median(v) if len(v) >= 10 else np.nan for v in comp[b_]]) for b_ in comp}
    write_csv(rows, Q16 / "q16g_humidity_mirror.csv")
    return rows, curves, 0.5 * (edges[:-1] + edges[1:])


def interannual(store, n_boot=500):
    from datetime import timedelta
    recs = []
    for ui, (s, (m, r, D, Y)) in enumerate(store.items()):
        if s == PAYERNE:
            continue
        dts = [EPOCH + timedelta(days=float(x)) for x in r["t"]]
        yr = np.array([d.year for d in dts])
        mo = np.array([d.month for d in dts])
        for mm in range(1, 13):
            a, b = (yr == 2025) & (mo == mm), (yr == 2026) & (mo == mm)
            if a.sum() < 4 or b.sum() < 4:
                continue
            rec = dict(unit=ui, month=mm)
            for k, v in list(Y.items()) + [(k, D[k]) for k in ("AH_night", "AH_ew21d", "RH_in_night", "RH_night")] + [("T_int", r["T_int"])]:
                va, vb = v[a][np.isfinite(v[a])], v[b][np.isfinite(v[b])]
                rec[k] = (vb.mean() - va.mean()) if (va.size >= 3 and vb.size >= 3) else np.nan
            recs.append(rec)
    rows = []
    for mm in range(1, 13):
        sel = [q for q in recs if q["month"] == mm]
        if len(sel) < 20:
            continue
        row = dict(test="month median 2026-2025", month=mm, n_units=len(sel))
        for k in ("lnC", "res_in", "x78", "AH_night", "AH_ew21d", "RH_in_night", "RH_night", "T_int"):
            v = np.array([q[k] for q in sel], float)
            row[k] = float(np.nanmedian(v)) * (100 if k in ("lnC", "res_in", "x78") else 1)
        rows.append(row)
    rng = np.random.default_rng(5)
    uu = np.unique([q["unit"] for q in recs])
    by_unit = {u: [i for i, q in enumerate(recs) if q["unit"] == u] for u in uu}

    def fit(yk, xk, idx):
        y = np.array([recs[i][yk] for i in idx], float)
        x = np.array([recs[i][xk] for i in idx], float)
        ok = np.isfinite(y) & np.isfinite(x)
        if ok.sum() < 30:
            return np.nan, np.nan
        c = np.linalg.lstsq(np.column_stack([np.ones(ok.sum()), x[ok]]), y[ok], rcond=None)[0]
        return c[1], np.corrcoef(x[ok], y[ok])[0, 1]
    for yk in ("res_in", "x78"):
        for xk in ("AH_night", "AH_ew21d", "RH_in_night"):
            g, r_ = fit(yk, xk, range(len(recs)))
            bs = [fit(yk, xk, [i for u in rng.choice(uu, len(uu)) for i in by_unit[u]])[0] for _ in range(n_boot)]
            rows.append(dict(test="pooled unit-month", y=yk, driver=xk, gamma=100 * g, gamma_se=100 * np.nanstd(bs), r=r_,
                             n=len(recs)))
    # month-level: network median dy vs network median dAH across months
    mm = [r_ for r_ in rows if r_["test"].startswith("month")]
    if len(mm) >= 5:
        for yk in ("res_in", "x78"):
            a = np.array([q[yk] for q in mm])
            b = np.array([q["AH_night"] for q in mm])
            rows.append(dict(test="months: network-median d" + yk + " vs dAH_night", r=float(np.corrcoef(a, b)[0, 1]),
                             gamma=float(np.polyfit(b, a, 1)[0]), n=len(mm)))
    write_csv(rows, Q16 / "q16g_humidity_interannual.csv")
    return rows


def med_iqr(v):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
    if v.size == 0:
        return np.nan, np.nan, np.nan
    return float(np.median(v)), float(np.percentile(v, 25)), float(np.percentile(v, 75))


def summary(unit_rows, drows, store):
    net = [u for u in unit_rows if u["stream"] != PAYERNE]
    lines = [f"units {len(net)} (Payerne A excluded)"]
    for k in ("lnC", "res_in", "res_m2", "x78"):
        a = med_iqr([u.get(f"amp_{k}") for u in net])
        lines.append(f"{k}: amp {a[0]:.1f} [{a[1]:.1f}, {a[2]:.1f}] %, peak {circ_median_doy([u.get(f'peak_{k}') for u in net]):.0f}")
    lines.append(f"T_int night: amp {med_iqr([u.get('amp_Tint') for u in net])[0]:.1f} K, peak "
                 f"{circ_median_doy([u.get('peak_Tint') for u in net]):.0f}")
    table = []
    for name in DRIVERS:
        dd = [d for d in drows if d["driver"] == name and d["stream"] != PAYERNE]
        if len(dd) < 10:
            continue
        row = dict(driver=name, unit=UNITS[name], n_units=len(dd))
        row["amp_D"] = med_iqr([d["amp_D"] for d in dd])[0]
        row["peak_D"] = circ_median_doy([d["peak_D"] for d in dd])
        row["mean_D"] = med_iqr([d["mean_D"] for d in dd])[0]
        for rk in ("res_in", "res_m2", "x78"):
            row[f"dphase_{rk}"] = med_iqr([d.get(f"dphase_{rk}") for d in dd])[0]
            row[f"needed_{rk}"] = med_iqr([d.get(f"needed_{rk}") for d in dd])[0]
        for yk in ("lnC", "x78", "xwin"):
            row[f"pcorr_{yk}"] = med_iqr([d.get(f"pcorr_{yk}") for d in dd])[0]
            g, se, gsd, nu = pooled(store, name, yk)
            row[f"g_pool_{yk}"], row[f"g_pool_se_{yk}"], row[f"g_pool_sd_{yk}"] = g, se, gsd
            row[f"transfer_{yk}_res_in"] = med_iqr([d.get(f"transfer_{yk}_res_in") for d in dd])[0]
        # latitude / climate: across units, residual amplitude vs the driver's amplitude
        from scipy.stats import spearmanr
        ua = {u["stream"]: u for u in net}
        a = [(ua[d["stream"]].get("amp_res_in"), d["amp_D"]) for d in dd if d["stream"] in ua]
        a = np.array([x for x in a if x[0] is not None and np.isfinite(x[0])], float)
        if len(a) > 10:
            row["rho_amp_res_vs_ampD"] = float(spearmanr(a[:, 0], a[:, 1]).correlation)
        a = np.array([(ua[d["stream"]].get("amp_x78"), d["amp_D"]) for d in dd if d["stream"] in ua], float)
        a = a[np.isfinite(a).all(axis=1)]
        if len(a) > 10:
            row["rho_amp_x78_vs_ampD"] = float(spearmanr(a[:, 0], a[:, 1]).correlation)
            row["ampD_ratio_hi_lo"] = float(np.percentile(a[:, 1], 75) / np.percentile(a[:, 1], 25))
            row["amp_x78_ratio_hi_lo_by_D"] = float(np.median(a[a[:, 1] >= np.median(a[:, 1]), 0]) /
                                                    np.median(a[a[:, 1] < np.median(a[:, 1]), 0]))
        table.append(row)
        lines.append(" | ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}" for k, v in row.items()))
    write_csv(table, Q16 / "q16g_humidity_table.csv")
    mrows, curves, mid = mirror(store)
    for r_ in mrows:
        lines.append(f"mirror {r_['series']:28s} {r_['spring_minus_autumn']:+.2f} +- {r_['se']:.2f} (n={r_['n_units']})")
    irows = interannual(store)
    for r_ in irows:
        lines.append("interannual " + ", ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}" for k, v in r_.items()))
    txt = "\n".join(lines)
    print(txt)
    (Q16 / "q16g_humidity_summary.txt").write_text(txt, encoding="utf-8")
    figure(store, unit_rows, curves, mid)


def figure(store, unit_rows, curves, mid):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
    # (a) composite annual course (unit mean removed, 15-day bins), normalised to the harmonic amplitude
    bins = np.arange(0, 366, 15)
    comp = {k: [] for k in ("res_in", "x78", "AH_night", "T_int", "RH_in_night")}
    for s, (m, r, D, Y) in store.items():
        if s == PAYERNE:
            continue
        doy = np.mod(r["t"], 365.25)
        for k in comp:
            y = Y[k] if k in Y else (r["T_int"] if k == "T_int" else D[k])
            if np.isfinite(y).sum() < 60:
                continue
            p, _, _ = harmonic(r["t"], y)
            if not np.isfinite(p) or abs(p) == 0:
                continue
            y = (y - np.nanmean(y)) / abs(p)
            comp[k].append([np.nanmean(y[(doy >= b) & (doy < b + 15)]) if ((doy >= b) & (doy < b + 15)).sum() >= 2 else np.nan
                            for b in bins])
    labs = {"res_in": "residual (kappa_in)", "x78": "clear-sky X(7-8 km)", "AH_night": "absolute humidity (night)",
            "T_int": "T_int (night)", "RH_in_night": "RH inside proxy e/e_s(T_int)"}
    for k, c in zip(comp, ("C0", "C2", "C3", "C1", "C4")):
        ax[0].plot(bins + 7.5, np.nanmedian(np.array(comp[k]), axis=0), "-o", ms=3, color=c, label=labs[k])
    ax[0].axvline(172, color="k", ls="--", lw=0.8)
    ax[0].set_xlabel("day of year")
    ax[0].set_ylabel("unit mean removed / own annual amplitude")
    ax[0].set_title("(a) normalised composite seasonal courses")
    ax[0].legend(fontsize=8)
    # (b) mirror composite: AH vs residual vs x78
    for lab, c in (("residual kappa_in", "C0"), ("ln X 7-8 km", "C2")):
        if lab in curves:
            ax[1].plot(mid, curves[lab]["spring"], "-o", color=c, ms=3, label=f"{lab}: spring")
            ax[1].plot(mid, curves[lab]["autumn"], "--s", color=c, ms=3, label=f"{lab}: autumn")
    a2 = ax[1].twinx()
    lab = "AH_night [g/m3]"
    if lab in curves:
        a2.plot(mid, curves[lab]["spring"], "-o", color="C3", ms=3, label="AH: spring")
        a2.plot(mid, curves[lab]["autumn"], "--s", color="C3", ms=3, label="AH: autumn")
    a2.set_ylabel("absolute humidity anomaly [g m-3]", color="C3")
    ax[1].set_xlabel("days from the solstice")
    ax[1].set_ylabel("% (unit mean removed)")
    ax[1].set_title("(b) solstice-mirror: spring (solid) vs autumn (dashed)")
    ax[1].legend(fontsize=7, loc="lower left")
    a2.legend(fontsize=7, loc="upper right")
    # (c) per-unit AH amplitude vs X(7-8 km) amplitude
    net = [u for u in unit_rows if u["stream"] != PAYERNE]
    x = np.array([u.get("amp_AH_night", np.nan) for u in net], float)
    y = np.array([u.get("amp_x78", np.nan) for u in net], float)
    lat = np.array([u["lat"] for u in net], float)
    sc = ax[2].scatter(x, y, c=lat, cmap="viridis", s=14)
    plt.colorbar(sc, ax=ax[2], label="latitude [deg N]")
    ax[2].set_xlabel("annual amplitude of night absolute humidity [g m-3]")
    ax[2].set_ylabel("annual amplitude of clear-sky X(7-8 km) [%]")
    ax[2].set_title("(c) climate contrast: humidity amplitude vs gain amplitude")
    fig.suptitle("q16g - absolute humidity (CAMS 0.4 deg, station level) as a driver of the CHM15k gain cycle", fontsize=10)
    fig.tight_layout()
    fig.savefig(Q16 / "q16g_humidity.png", dpi=130)


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
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", action="store_true")
    a = ap.parse_args()
    import warnings
    warnings.filterwarnings("ignore")
    if a.extract or not CAMS_NPZ.exists():
        stage_extract()
    stage_analyse()
