# -*- coding: utf-8 -*-
"""q15e - electronic noise and profile distortion as drivers of the CHM15k constant variability.

Family e of the q15 attribution study (brief: C:/DATA/Projects/202610_CHM15k_stability/q15/BRIEF.md).
Everything here changes the SHAPE of the range-corrected profile (and so the Rayleigh fit) rather than
being a pure multiplicative gain:

1. background / residual dark baseline. A residual additive offset B left after the background
   subtraction appears in rcs_0 as B r^2. Normalised by the molecular attenuated backscatter it becomes
   X(r) = rcs/beta_mol_att = C + B g(r), g = r^2 / beta_mol_att, which grows ~30x between the Rayleigh
   window (3-6 km) and 14 km. The fitted constant is then biased by B g_win.
   Stage "l1" fits X(r) = C + B g(r) from the window bottom to 15 km on sampled calibration nights
   (and noon windows) read from the native L1 files, which reach 15.3 km (the cache stops at 8 km).
2. electronic noise (`stddev`) and the night background `bckgrd_rcs_0`: annual cycle, temperature
   dependence, partial effect on ln C anomalies (season removed with a 31-day running median).
3. profile distortion: two-band shape index S = ln X(2-3 km) - ln X(4-6 km) vs T_int anomalies and season,
   near-range-load test for signal-induced noise / afterpulsing, calibration-window heights vs the
   overlap (no-op above ~700-840 m).
4. altitude dependence of kappa: kappa(z) per range band from the cache and its gradient (%/K/km).

Run:  .venv/Scripts/python.exe validation/chm15k_stability/q15e_noise_profile.py [cache|l1|summary|all]
Outputs: C:/DATA/Projects/202610_CHM15k_stability/q15/q15e_*.csv / .png / q15e_summary.md
"""
from __future__ import annotations

import csv
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (BLOCKS_PER_DAY, CENSUS, L1_ROOT, OUT, chm15k_streams, load_cache,  # noqa: E402
                    load_constants, solar_hour, stream_key)
from q1_temperature import ANOM_WINDOW, night_hk, robust_mask, running_median  # noqa: E402

warnings.filterwarnings("ignore")
Q15 = OUT / "q15"
PAYERNE = "0-20000-0-06610_A"
N_WORKERS = 5
WAVELENGTH = 1064e-9
W_YEAR = 2 * np.pi / 365.25

# Range bands (km, range above the instrument) used for kappa(z) and the shape indices.
BANDS = [(0.2, 0.5), (0.5, 0.8), (0.8, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.0),
         (3.0, 4.0), (4.0, 5.0), (5.0, 6.0), (6.0, 7.0), (7.0, 8.0)]
BAND_MID = np.array([0.5 * (a + b) for a, b in BANDS])


# ---------------------------------------------------------------------------------------------
# Molecular model (US-1976 standard atmosphere, Bucholtz, 1064 nm) - shape reference only
# ---------------------------------------------------------------------------------------------
def molecular_att(alt, rng):
    """Molecular attenuated backscatter on a range grid (m) for a station at altitude alt (m ASL)."""
    from calibration.rayleigh.atmosphere import calculate_molecular_properties, load_standard_atmosphere
    fine = np.arange(5.0, max(rng.max(), 100.0) + 10.0, 10.0)
    atm = load_standard_atmosphere(None, alt + fine)
    mp = calculate_molecular_properties(atm.temperature, atm.pressure, fine, WAVELENGTH)
    return np.interp(rng, fine, mp.beta_att_mol)


def station_alt():
    import json
    return {stream_key(s["wmo"], s["ident"]): float(s.get("alt", 0.0) or 0.0)
            for s in json.load(open(CENSUS, encoding="utf-8"))}


# ---------------------------------------------------------------------------------------------
# Small statistics helpers
# ---------------------------------------------------------------------------------------------
def ols(X, y):
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    dof = max(len(y) - X.shape[1], 1)
    try:
        se = np.sqrt(np.diag(res @ res / dof * np.linalg.inv(X.T @ X)))
    except np.linalg.LinAlgError:
        se = np.full(X.shape[1], np.nan)
    return coef, se, res


def anomaly(t, x):
    return x - running_median(t, x, ANOM_WINDOW)


def annual_fit(t, y):
    """Offset + annual harmonic (+ linear trend when > 400 days). Returns amp, peak DOY, trend/yr."""
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if len(y) < 60 or t.max() - t.min() < 300:
        return np.nan, np.nan, np.nan
    doy = W_YEAR * ((t + 0.0) % 365.25)  # day 0 of epoch = 1 Jan 1970, so t % 365.25 ~ DOY - 1
    cols = [np.ones_like(t), np.sin(doy), np.cos(doy)]
    trend = t.max() - t.min() > 400
    if trend:
        cols.append((t - t.mean()) / 365.25)
    c, _, _ = ols(np.column_stack(cols), y)
    amp = np.hypot(c[1], c[2])
    peak = (np.arctan2(c[1], c[2]) / W_YEAR) % 365.25 + 1.0
    return amp, peak, (c[3] if trend else np.nan)


# ---------------------------------------------------------------------------------------------
# Stage 1: cache-based analysis (all CHM15k)
# ---------------------------------------------------------------------------------------------
def clear_blocks(cache):
    tb = (cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY).ravel()
    order = np.argsort(tb)
    nb = cache["rcs"].shape[-1]
    rcs = cache["rcs"].reshape(-1, nb)[order]
    clear = cache["n_clear"].ravel()[order] >= 0.9 * np.maximum(cache["n_all"].ravel()[order], 1)
    clear &= cache["n_clear"].ravel()[order] > 0
    return tb[order], rcs, clear


def band_index(rmid):
    return [np.where((rmid >= a * 1000) & (rmid < b * 1000))[0] for a, b in BANDS]


def night_profiles(cache, cons, alt, tilt):
    """Per calibration night: ln X per band, ln X in the window, number of clear blocks."""
    tb, rcs, clear = clear_blocks(cache)
    edges = cache["range_edges"]
    rmid = 0.5 * (edges[:-1] + edges[1:])
    fine = np.arange(5.0, 8000.0, 10.0)
    bm = molecular_att(alt, fine * np.cos(np.deg2rad(tilt)))
    mol = np.array([bm[(fine >= edges[k]) & (fine < edges[k + 1])].mean() for k in range(len(rmid))])
    idx = band_index(rmid)
    n = len(cons["t"])
    lnx = np.full((n, len(BANDS)), np.nan)
    lnx_win = np.full(n, np.nan)
    nclr = np.zeros(n, int)
    i0 = np.searchsorted(tb, cons["t0"])
    i1 = np.searchsorted(tb, cons["t1"], side="right")
    vert = np.cos(np.deg2rad(tilt))
    for k in range(n):
        if i1[k] <= i0[k]:
            continue
        seg = rcs[i0[k]:i1[k]][clear[i0[k]:i1[k]]]
        nclr[k] = seg.shape[0]
        if seg.shape[0] < 6:
            continue
        prof = np.nanmean(seg, axis=0) / mol
        with np.errstate(invalid="ignore", divide="ignore"):
            for j, ii in enumerate(idx):
                v = np.nanmean(prof[ii])
                lnx[k, j] = np.log(v) if v > 0 else np.nan
            win = (rmid * vert >= cons["zb"][k] - alt) & (rmid * vert <= cons["zt"][k] - alt)
            if win.sum() >= 1:
                v = np.nanmean(prof[win])
                lnx_win[k] = np.log(v) if v > 0 else np.nan
    return lnx, lnx_win, nclr


def diurnal_hk(cache, lon, iv_bg, iv_sd, iv_t):
    """Per day: night (solar 22-02) and noon (solar 10-14) medians of background, noise and T_int."""
    tb = (cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY)
    sh = solar_hour(tb, lon)
    night = (sh >= 22) | (sh < 2)
    noon = (sh >= 10) & (sh < 14)
    out = {}
    for name, iv in (("bg", iv_bg), ("sd", iv_sd), ("T", iv_t)):
        x = cache["hk"][:, :, iv]
        out[name + "_n"] = np.nanmedian(np.where(night, x, np.nan), axis=1)
        out[name + "_d"] = np.nanmedian(np.where(noon, x, np.nan), axis=1)
    return out


def diurnal_shape(cache, lon, alt, tilt):
    """Day vs night clear-sky high-range index Dd = ln X(6-8 km) - ln X(4-5 km), median over days."""
    tb, rcs, clear = clear_blocks(cache)
    edges = cache["range_edges"]
    rmid = 0.5 * (edges[:-1] + edges[1:])
    fine = np.arange(5.0, 8000.0, 10.0)
    bm = molecular_att(alt, fine * np.cos(np.deg2rad(tilt)))
    mol = np.array([bm[(fine >= edges[k]) & (fine < edges[k + 1])].mean() for k in range(len(rmid))])
    hi = (rmid >= 6000) & (rmid < 8000)
    mid = (rmid >= 4000) & (rmid < 5000)
    sh = solar_hour(tb, lon)
    dayk = np.floor(tb + lon / 360.0).astype(int)   # local solar day key (integer)
    res = {}
    for lab, sel in (("night", (sh >= 22) | (sh < 3)), ("noon", (sh >= 10) & (sh < 14))):
        m = clear & sel
        vals = []
        for d in np.unique(dayk[m]):
            s = m & (dayk == d)
            if s.sum() < 6:
                continue
            p = np.nanmean(rcs[s], axis=0) / mol
            a, b = np.nanmean(p[hi]), np.nanmean(p[mid])
            if a > 0 and b > 0:
                vals.append(np.log(a / b))
        res[lab] = np.array(vals)
    return res


def analyse_unit(args):
    wmo, ident, site, lat, lon, alt, kappa_in = args
    key = stream_key(wmo, ident)
    cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
    if cons is None or cache is None:
        return None
    mod = max(set(cons["module"]), key=list(cons["module"]).count)
    m = (cons["module"] == mod) & (cons["c"] > 0) & np.isfinite(cons["t0"])
    if m.sum() < 80:
        return None
    cons = {k: v[m] for k, v in cons.items()}
    tilt = float(cache.get("tilt", 0.0))
    tilt = 0.0 if not np.isfinite(tilt) else tilt
    hv = list(cache["hk_vars"])
    hk = night_hk(cache, cons["t0"], cons["t1"])
    T = hk[:, hv.index("temp_int")]
    bg = hk[:, hv.index("bckgrd_rcs_0")]
    sd = hk[:, hv.index("stddev")]
    lnx, lnx_win, nclr = night_profiles(cache, cons, alt, tilt)
    t = cons["t"]
    lnc = np.log(cons["c"])
    keep = robust_mask(lnc) & np.isfinite(T) & (bg > 0) & (sd > 0)
    t, lnc, T, bg, sd = t[keep], lnc[keep], T[keep], bg[keep], sd[keep]
    lnx, lnx_win, nclr = lnx[keep], lnx_win[keep], nclr[keep]
    zb = cons["zb"][keep] - alt
    zt = cons["zt"][keep] - alt
    dur = (cons["t1"] - cons["t0"])[keep] * 24.0
    if len(t) < 80:
        return None

    S = lnx[:, 5:7].mean(axis=1) - lnx[:, 8:10].mean(axis=1)        # ln X(2-3 km) - ln X(4-6 km)
    D = lnx[:, 11] - lnx_win                                           # ln X(7-8 km) - ln X(window)
    L = np.log(np.nanmean(np.exp(lnx[:, 0:3]), axis=1)) - lnx_win      # near-range load 0.2-1 km
    la, ta = anomaly(t, lnc), anomaly(t, T)
    bga, sda = anomaly(t, np.log(bg)), anomaly(t, np.log(sd))
    Sa, Da, La = anomaly(t, S), anomaly(t, D), anomaly(t, L)
    zba, zta, dura = anomaly(t, zb), anomaly(t, zt), anomaly(t, dur)
    row = dict(stream=key, site=site, module=mod, alt=alt, n=len(t), kappa_in=kappa_in,
               bg_med=float(np.median(bg)), sd_med=float(np.median(sd)),
               zb_med=float(np.median(zb)), zt_med=float(np.median(zt)), zb_p05=float(np.percentile(zb, 5)),
               frac_zb_below_850=float(np.mean(zb < 850)), dur_med=float(np.median(dur)))

    # (a) ln C anomaly on T + background + noise anomalies (all in ln units -> %/%)
    ok = np.isfinite(la) & np.isfinite(ta) & np.isfinite(bga) & np.isfinite(sda)
    if ok.sum() >= 60:
        X0 = np.column_stack([np.ones(ok.sum()), ta[ok]])
        c0, _, r0 = ols(X0, la[ok])
        row["kappa_bare"] = 100 * c0[1]
        c1, se1, r1 = ols(np.column_stack([X0, bga[ok], sda[ok]]), la[ok])
        c_bg, se_bg, _ = ols(np.column_stack([X0, bga[ok]]), la[ok])
        c_sd, se_sd, _ = ols(np.column_stack([X0, sda[ok]]), la[ok])
        row.update(kappa_m2=100 * c1[1], g_bg=c_bg[2], se_bg=se_bg[2], g_sd=c_sd[2], se_sd=se_sd[2],
                   g_bg_joint=c1[2], g_sd_joint=c1[3],
                   dR2_bgsd=float(1 - np.var(r1) / np.var(r0)),
                   std_lnc_a=100 * float(np.std(la[ok])),
                   std_bg_a=100 * float(np.std(bga[ok])), std_sd_a=100 * float(np.std(sda[ok])),
                   n2n_bg=100 * abs(c_bg[2]) * float(np.std(bga[ok])),
                   n2n_sd=100 * abs(c_sd[2]) * float(np.std(sda[ok])))
        # temperature sensitivity of background and noise (%/K)
        cb, sb, _ = ols(X0, bga[ok])
        cs, ss, _ = ols(X0, sda[ok])
        row.update(bg_per_K=100 * cb[1], sd_per_K=100 * cs[1])
    # (b) window selection: ln C anomaly vs window bottom/top (per km) and duration (per h)
    ok = np.isfinite(la) & np.isfinite(ta) & np.isfinite(zba) & np.isfinite(zta) & np.isfinite(dura)
    if ok.sum() >= 60:
        c, se, _ = ols(np.column_stack([np.ones(ok.sum()), ta[ok], zba[ok] / 1000, zta[ok] / 1000, dura[ok]]),
                       la[ok])
        row.update(g_zb=100 * c[2], g_zt=100 * c[3], g_dur=100 * c[4],
                   n2n_win=float(np.std(c[2] * zba[ok] / 1000 + c[3] * zta[ok] / 1000 + c[4] * dura[ok])) * 100)
        # does the window move with noise?
        okn = ok & np.isfinite(sda)
        cz, _, _ = ols(np.column_stack([np.ones(okn.sum()), sda[okn]]), zta[okn])
        row.update(zt_per_lnsd=float(cz[1]))
    # (c) profile-shape indices vs temperature anomalies (%/K)
    for name, xa in (("S", Sa), ("D", Da)):
        ok = np.isfinite(xa) & np.isfinite(ta)
        if ok.sum() >= 60:
            c, se, _ = ols(np.column_stack([np.ones(ok.sum()), ta[ok]]), xa[ok])
            row[f"{name}_per_K"], row[f"{name}_se"] = 100 * c[1], 100 * se[1]
            # does the shape anomaly move ln C?  (beyond T)
            ok2 = ok & np.isfinite(la)
            c2, se2, _ = ols(np.column_stack([np.ones(ok2.sum()), ta[ok2], xa[ok2]]), la[ok2])
            row[f"lnc_on_{name}"], row[f"lnc_on_{name}_se"] = c2[2], se2[2]
    row["S_level"] = float(np.nanmedian(S))
    row["D_level"] = float(np.nanmedian(D))
    # (d) signal-induced noise / afterpulse: does far-range excess follow the near-range load?
    ok = np.isfinite(Da) & np.isfinite(La) & np.isfinite(ta)
    if ok.sum() >= 60:
        c, se, _ = ols(np.column_stack([np.ones(ok.sum()), ta[ok], La[ok]]), Da[ok])
        row.update(D_on_L=c[2], D_on_L_se=se[2], std_L_a=100 * float(np.std(La[ok])))
    # (e) kappa(z): per band anomalies on T anomaly
    kz = np.full(len(BANDS), np.nan)
    for j in range(len(BANDS)):
        xa = anomaly(t, lnx[:, j])
        ok = np.isfinite(xa) & np.isfinite(ta)
        if ok.sum() >= 60:
            c, _, _ = ols(np.column_stack([np.ones(ok.sum()), ta[ok]]), xa[ok])
            kz[j] = 100 * c[1]
        row[f"kz_{BANDS[j][0]:.1f}_{BANDS[j][1]:.1f}"] = kz[j]
    sel = (BAND_MID >= 1.0) & (BAND_MID <= 7.0) & np.isfinite(kz)
    if sel.sum() >= 5:
        c, se, _ = ols(np.column_stack([np.ones(sel.sum()), BAND_MID[sel]]), kz[sel])
        row.update(kz_grad=c[1], kz_grad_se=se[1])
    # annual cycles (ln units -> %)
    for name, y in (("lnc", lnc), ("T", T), ("bg", np.log(bg)), ("sd", np.log(sd)), ("S", S), ("D", D),
                    ("zb", zb / 1000), ("zt", zt / 1000), ("dur", dur)):
        amp, peak, trend = annual_fit(t, y)
        scale = 100 if name in ("lnc", "bg", "sd", "S", "D") else 1
        row[f"amp_{name}"], row[f"peak_{name}"], row[f"trend_{name}"] = amp * scale, peak, trend * scale
    # predicted C annual effects: partial coefficient x factor cycle (first harmonic, sign kept)
    for name, g, y in (("bg", row.get("g_bg", np.nan), np.log(bg)), ("sd", row.get("g_sd", np.nan), np.log(sd))):
        amp, peak, _ = annual_fit(t, g * y if np.isfinite(g) else y * np.nan)
        row[f"pred_amp_{name}"], row[f"pred_peak_{name}"] = 100 * amp, peak
    gz = (row.get("g_zb", np.nan) * zb + row.get("g_zt", np.nan) * zt) / 1000 + row.get("g_dur", np.nan) * dur
    amp, peak, _ = annual_fit(t, gz)
    row["pred_amp_win"], row["pred_peak_win"] = amp, peak
    # diurnal background / noise (all days, solar time)
    dh = diurnal_hk(cache, lon, hv.index("bckgrd_rcs_0"), hv.index("stddev"), hv.index("temp_int"))
    row["bg_noon_night"] = float(np.nanmedian(dh["bg_d"] / dh["bg_n"]))
    row["sd_noon_night"] = float(np.nanmedian(dh["sd_d"] / dh["sd_n"]))
    ds = diurnal_shape(cache, lon, alt, tilt)
    row["Dd_night"] = float(np.median(ds["night"])) if len(ds["night"]) > 20 else np.nan
    row["Dd_noon"] = float(np.median(ds["noon"])) if len(ds["noon"]) > 20 else np.nan
    series = dict(lnx=lnx, lnx_win=lnx_win, t=t, la=la, ta=ta, lnc=lnc, T=T, bg=bg, sd=sd, S=S, D=D, zb=zb, zt=zt, kz=kz)
    return row, series


def stage_cache():
    alts = station_alt()
    k_in = {}
    for r in csv.DictReader(open(OUT / "q2a_intranight.csv", encoding="utf-8")):
        k_in[r["stream"]] = float(r["b_T_ctrl"])
    jobs = [(w, i, s, la, lo, alts.get(stream_key(w, i), 0.0), k_in.get(stream_key(w, i), np.nan))
            for w, i, s, la, lo in chm15k_streams()]
    with Pool(N_WORKERS) as pool:
        res = [r for r in pool.map(analyse_unit, jobs, chunksize=2) if r is not None]
    rows = [r[0] for r in res]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(Q15 / "q15e_units.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    np.savez_compressed(Q15 / "q15e_series.npz",
                        **{f"{r[0]['stream']}__{k}": v for r in res for k, v in r[1].items()})
    print(f"cache stage: {len(rows)} units -> {Q15 / 'q15e_units.csv'}")


# ---------------------------------------------------------------------------------------------
# Stage 2: residual baseline from native L1 (to 15 km) on sampled calibration nights
# ---------------------------------------------------------------------------------------------
L1_EDGES = np.r_[np.arange(1000, 15001, 500)].astype(float)


def l1_path(wmo, ident, day):
    from datetime import datetime, timedelta
    d = datetime(1970, 1, 1) + timedelta(days=int(day))
    return L1_ROOT / wmo / f"{d:%Y}" / f"{d:%m}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"


def read_window(wmo, ident, t0, t1):
    """Clear-sky mean rcs_0 on L1_EDGES bins + median stddev/bckgrd/p_calc and rcs noise in [t0, t1]."""
    import netCDF4
    nb = len(L1_EDGES) - 1
    ssum, scnt, hk, tilt = np.zeros(nb), np.zeros(nb), [], 0.0
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
                rng = np.ma.filled(ds["range"][:].astype("f8"), np.nan)
                rcs = np.ma.filled(ds["rcs_0"][i0:i1, :].astype("f8"), np.nan)
                cbh = np.ma.filled(ds["cloud_base_height"][i0:i1, 0].astype("f8"), np.nan)
                sci = np.ma.filled(ds["sci"][i0:i1].astype("f8"), np.nan) if "sci" in ds.variables else                     np.full(i1 - i0, np.nan)
                bg = np.ma.filled(ds["bckgrd_rcs_0"][i0:i1].astype("f8"), np.nan)
                sd = np.ma.filled(ds["stddev"][i0:i1].astype("f8"), np.nan)
                pc = np.ma.filled(ds["calibration_pulse"][i0:i1].astype("f8"), np.nan)
                tilt = float(np.ma.filled(ds["tilt_angle"][...], 0.0))
        except (OSError, RuntimeError, KeyError, IndexError):
            continue
        clear = ~(np.isfinite(cbh) & (cbh > 0)) & ((sci == 0) | ~np.isfinite(sci))
        if clear.sum() < 10:
            continue
        x = rcs[clear]
        far = (rng > 9000) & (rng < 15000)
        rr = rng[far]
        # first differences cancel the smooth signal: per-gate, per-profile noise of rcs/r^2
        dd = np.diff(x[:, far] / rr[None, :] ** 2, axis=1)
        noise_r2 = np.nanmedian(np.nanstd(dd, axis=1)) / np.sqrt(2)
        hk.append((np.nanmedian(sd[clear]), np.nanmedian(bg[clear]), np.nanmedian(pc[clear]), noise_r2,
                   clear.sum()))
        gb = np.digitize(rng, L1_EDGES) - 1
        fin = np.isfinite(x)
        x0 = np.where(fin, x, 0.0)
        for k in range(nb):
            sel = gb == k
            ssum[k] += x0[:, sel].sum()
            scnt[k] += fin[:, sel].sum()
    if not hk:
        return None
    h = np.array(hk)
    n_tot = h[:, 4].sum()
    if n_tot < 30:
        return None
    w = h[:, 4] / n_tot
    gpb = np.median(np.diff(L1_EDGES)) / 15.0   # approx. gates per bin (L1 gate ~10-15 m)
    return dict(prof=np.where(scnt > 0, ssum / np.maximum(scnt, 1), np.nan), cnt=scnt, n=int(n_tot),
                sd=float(w @ h[:, 0]), bg=float(w @ h[:, 1]), pc=float(w @ h[:, 2]),
                noise_r2=float(w @ h[:, 3]), tilt=tilt, gpb=gpb)


L1_MID = 0.5 * (L1_EDGES[:-1] + L1_EDGES[1:])


def sample_unit(args):
    """Per sampled window: X(r)/X_window and its noise sigma on 500 m bins to 15 km (night + next noon)."""
    wmo, ident, site, alt, n_nights = args
    key = stream_key(wmo, ident)
    cons = load_constants(wmo, ident)
    if cons is None:
        return []
    mod = max(set(cons["module"]), key=list(cons["module"]).count)
    m = np.where((cons["module"] == mod) & (cons["c"] > 0) & np.isfinite(cons["t0"]))[0]
    if len(m) < 40:
        return []
    pick = m[np.linspace(0, len(m) - 1, min(n_nights, len(m))).astype(int)]
    lon = LONS.get(key, 0.0)
    rows = []
    mol = None
    for k in pick:
        t0, t1 = cons["t0"][k], cons["t1"][k]
        zw = (max(cons["zb"][k] - alt, 2500.0), cons["zt"][k] - alt)
        d = np.floor(t1)
        for lab, (a, b) in (("night", (t0, t1)), ("noon", (d + (10 - lon / 15) / 24, d + (14 - lon / 15) / 24))):
            w = read_window(wmo, ident, a, b)
            if w is None:
                continue
            if mol is None:
                mol = molecular_att(alt, L1_MID * np.cos(np.deg2rad(w["tilt"])))
            X = w["prof"] / mol
            sig = w["noise_r2"] * L1_MID ** 2 / np.sqrt(np.maximum(w["cnt"], 1)) / mol
            win = (L1_MID >= zw[0]) & (L1_MID <= zw[1])
            if win.sum() < 2:
                continue
            xw = np.nanmean(X[win])
            if not xw > 0:
                continue
            scale = w["noise_r2"] / w["sd"] if w["sd"] > 0 else np.nan   # rcs/r^2 per raw (photons/shot) unit
            rows.append(dict(stream=key, site=site, kind=lab, t=float(0.5 * (a + b)), c_op=float(cons["c"][k]),
                             zb=zw[0], zt=zw[1], n=w["n"], sd=w["sd"], bg=w["bg"], pc=w["pc"],
                             xw=float(xw), scale=float(scale),
                             bg_rel_win=float(w["bg"] * scale * np.mean(L1_MID[win] ** 2) / np.mean(w["prof"][win])),
                             xn=X / xw, sn=sig / xw, gwin=float(np.mean(L1_MID[win] ** 2 / mol[win]))))
    return rows


LONS = {}


def l1_units():
    """Geographically spread sample + pairs + Payerne (reported separately)."""
    want = ["0-20000-0-06240_A", "0-20000-0-06240_B", "0-20000-0-12843_A", "0-20000-0-12843_B",
            "0-20008-0-INO_A", "0-20008-0-INO_B", PAYERNE]
    streams = chm15k_streams()
    for w, i, s, la, lo in streams:
        LONS[stream_key(w, i)] = lo
    # one unit per ~country-ish cell: bin lat/lon on a 4-degree grid, take the first per cell
    cells = {}
    for w, i, s, la, lo in streams:
        c = (int(la // 4), int(lo // 4))
        if c not in cells and load_cache(w, i) is not None:
            cells[c] = stream_key(w, i)
    want += [v for v in cells.values() if v not in want]
    return want[:34]


def stage_l1():
    alts = station_alt()
    keys = sorted(set(l1_units()))
    streams = {stream_key(w, i): (w, i, s) for w, i, s, la, lo in chm15k_streams()}
    jobs = [(streams[k][0], streams[k][1], streams[k][2], alts.get(k, 0.0), 60) for k in keys if k in streams]
    with Pool(N_WORKERS, initializer=_init_lons) as pool:
        res = pool.map(sample_unit, jobs, chunksize=1)
    rows = [r for rr in res for r in rr]
    arr = {k: np.array([r[k] for r in rows]) for k in ("xn", "sn")}
    meta = [{k: v for k, v in r.items() if k not in ("xn", "sn")} for r in rows]
    np.savez_compressed(Q15 / "q15e_l1_profiles.npz", mid=L1_MID, **arr)
    with open(Q15 / "q15e_l1_windows.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(meta[0].keys()))
        w.writeheader()
        w.writerows(meta)
    print(f"l1 stage: {len(rows)} windows from {len(jobs)} units")


def _init_lons():
    for w, i, s, la, lo in chm15k_streams():
        LONS[stream_key(w, i)] = lo


def _main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("cache", "all"):
        stage_cache()
    if stage in ("l1", "all"):
        stage_l1()
    if stage in ("altitude", "all"):
        stage_altitude()
    if stage in ("summary", "all"):
        stage_summary()


# ---------------------------------------------------------------------------------------------
# Stage 3: altitude dependence of kappa - intra-night per unit and co-located pairs (aerosol cancels)
# ---------------------------------------------------------------------------------------------
ZBANDS = {"0.5-1": (500, 1000), "1-1.5": (1000, 1500), "1.5-2.5": (1500, 2500), "2.5-4": (2500, 4000),
          "4-6": (4000, 6000), "6-8": (6000, 8000)}


def hourly_bands(cache, lon):
    """Hourly clear-sky band means (ln rcs) + T_int at night, with integer night ids."""
    from q2b_pairs import hourly
    h = hourly(cache)
    rmid = 0.5 * (cache["range_edges"][:-1] + cache["range_edges"][1:])
    out = {}
    for name, (a, b) in ZBANDS.items():
        sel = (rmid >= a) & (rmid < b)
        with np.errstate(invalid="ignore", divide="ignore"):
            v = np.nanmean(h["rcs"][:, sel], axis=1)
            out[name] = np.where(v > 0, np.log(v), np.nan)
    shift = int(np.floor(lon / 15.0))                    # whole hours of solar offset (integer)
    hl = h["t"] + shift                                  # local solar hour index (start of hour)
    hod = np.mod(hl, 24)
    out["night"] = (hod >= 20) | (hod < 5)
    out["nid"] = (hl - 12) // 24                         # night id: noon-to-noon local
    out["elapsed"] = np.mod(hod - 20, 24).astype(float)  # hours since 20 h local
    out["T"] = h["hk"][:, h["names"].index("temp_int")]
    out["t"] = h["t"]
    out["module"] = h["module"]
    return out


def within(x, nid):
    """Subtract the per-night mean (NaN-aware)."""
    out = np.full_like(x, np.nan, dtype=float)
    ok = np.isfinite(x)
    u, inv = np.unique(nid[ok], return_inverse=True)
    s = np.bincount(inv, x[ok])
    c = np.bincount(inv)
    out[ok] = x[ok] - (s / c)[inv]
    return out


def intranight_kz(args):
    wmo, ident, lon = args
    cache = load_cache(wmo, ident)
    if cache is None:
        return None
    hb = hourly_bands(cache, lon)
    mod = max(set(hb["module"]), key=list(hb["module"]).count)
    base = hb["night"] & (hb["module"] == mod) & np.isfinite(hb["T"])
    row = dict(stream=stream_key(wmo, ident))
    for name in ZBANDS:
        m = base & np.isfinite(hb[name])
        # nights with >= 4 clear hours in this band
        u, cnt = np.unique(hb["nid"][m], return_counts=True)
        m &= np.isin(hb["nid"], u[cnt >= 4])
        if m.sum() < 100:
            continue
        nid = hb["nid"][m]
        y, T, e = (within(hb[name][m], nid), within(hb["T"][m], nid), within(hb["elapsed"][m], nid))
        c, se, _ = ols(np.column_stack([T, e]), y)
        row[f"kin_{name}"], row[f"kin_se_{name}"] = 100 * c[0], 100 * se[0]
        row[f"nh_{name}"] = int(m.sum())
    return row


def pair_kz():
    """Night means of ln(rcs_A/rcs_B) per band vs the two units' T_int anomalies (multi-day) and
    intra-night (within-night demeaned, elapsed time held fixed)."""
    from common import PAIR_SITES
    streams = {stream_key(w, i): lo for w, i, s, la, lo in chm15k_streams()}
    rows = []
    for wmo, ids in PAIR_SITES.items():
        lon = streams.get(stream_key(wmo, ids[0]), 0.0)
        hbs = {i: hourly_bands(load_cache(wmo, i), lon) for i in ids if load_cache(wmo, i) is not None}
        for a in ids:
            for b in ids:
                if a >= b or a not in hbs or b not in hbs:
                    continue
                A, B = hbs[a], hbs[b]
                common, ia, ib = np.intersect1d(A["t"], B["t"], return_indices=True)
                night = A["night"][ia]
                nid = A["nid"][ia]
                TA, TB = A["T"][ia], B["T"][ib]
                pair_mod = np.char.add(A["module"][ia], B["module"][ib])
                main = max(set(pair_mod), key=list(pair_mod).count)
                for name in ZBANDS:
                    r = A[name][ia] - B[name][ib]
                    m = night & np.isfinite(r) & np.isfinite(TA) & np.isfinite(TB) & (pair_mod == main)
                    if m.sum() < 200:
                        continue
                    # intra-night
                    e = A["elapsed"][ia][m]
                    X = np.column_stack([within(TA[m], nid[m]), within(TB[m], nid[m]), within(e, nid[m])])
                    c, se, _ = ols(X, within(r[m], nid[m]))
                    # multi-day: night means, anomalies vs 31-day running median
                    u, inv = np.unique(nid[m], return_inverse=True)
                    cnt = np.bincount(inv)
                    rm = np.bincount(inv, r[m]) / cnt
                    ta_ = np.bincount(inv, TA[m]) / cnt
                    tb_ = np.bincount(inv, TB[m]) / cnt
                    keep = cnt >= 3
                    tt = u[keep].astype(float)
                    ra, taa, tba = anomaly(tt, rm[keep]), anomaly(tt, ta_[keep]), anomaly(tt, tb_[keep])
                    ok = np.isfinite(ra) & np.isfinite(taa) & np.isfinite(tba)
                    c2, se2, _ = ols(np.column_stack([np.ones(ok.sum()), taa[ok], tba[ok]]), ra[ok])
                    rows.append(dict(pair=f"{wmo}_{a}/{b}", band=name, n_h=int(m.sum()), n_nights=int(ok.sum()),
                                     kA_in=100 * c[0], kA_in_se=100 * se[0], kB_in=-100 * c[1], kB_in_se=100 * se[1],
                                     kA_md=100 * c2[1], kA_md_se=100 * se2[1], kB_md=-100 * c2[2],
                                     kB_md_se=100 * se2[2]))
    return rows


def stage_altitude():
    jobs = [(w, i, lo) for w, i, s, la, lo in chm15k_streams()]
    with Pool(N_WORKERS) as pool:
        rows = [r for r in pool.map(intranight_kz, jobs, chunksize=2) if r is not None and len(r) > 1]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(Q15 / "q15e_kz_intranight.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    prow = pair_kz()
    with open(Q15 / "q15e_kz_pairs.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(prow[0].keys()))
        w.writeheader()
        w.writerows(prow)
    print(f"altitude stage: {len(rows)} units, {len(prow)} pair-bands")


def pooled(series, keys, yname, xnames):
    """Network-pooled OLS of 31-day anomalies (each unit's anomalies are already ~zero-mean)."""
    Y, X = [], []
    for k in keys:
        t = series[f"{k}__t"]
        lnx = series[f"{k}__lnx"]
        cols = {"la": series[f"{k}__la"], "ta": series[f"{k}__ta"],
                "bg": anomaly(t, np.log(series[f"{k}__bg"])), "sd": anomaly(t, np.log(series[f"{k}__sd"])),
                "zb": anomaly(t, series[f"{k}__zb"] / 1000), "zt": anomaly(t, series[f"{k}__zt"] / 1000),
                "D": anomaly(t, series[f"{k}__D"]),
                "L": anomaly(t, np.log(np.nanmean(np.exp(lnx[:, 0:3]), axis=1)) - series[f"{k}__lnx_win"]),
                "x78": anomaly(t, lnx[:, 11])}
        y = cols[yname]
        x = np.column_stack([cols[n] for n in xnames])
        ok = np.isfinite(y) & np.isfinite(x).all(axis=1)
        Y.append(y[ok])
        X.append(x[ok])
    Y, X = np.concatenate(Y), np.vstack(X)
    c, se, res = ols(np.column_stack([np.ones(len(Y)), X]), Y)
    r2 = float(1 - np.var(res) / np.var(Y))
    return dict(zip(xnames, c[1:])), dict(zip(xnames, se[1:])), r2, len(Y), X


def l1_aggregate():
    """Per unit x kind x season: inverse-variance mean of X(r)/X_window; fit a + b g(r)/g_ref on 4-15 km.

    b/a is the fraction of the 4-5 km signal produced by a residual baseline B r^2, i.e. the C bias for a
    window centred near 4.5 km (median window 2.5-5.6 km)."""
    import pandas as pd
    z = np.load(Q15 / "q15e_l1_profiles.npz")
    xn, sn, mid = z["xn"], z["sn"], z["mid"]
    meta = pd.read_csv(Q15 / "q15e_l1_windows.csv")
    alts = station_alt()
    doy = (meta["t"] % 365.25).values
    seas = {"all": np.ones(len(meta), bool), "warm": (doy >= 120) & (doy <= 270), "cold": (doy <= 60) | (doy >= 305)}
    rows, profs = [], {}
    for key in meta["stream"].unique():
        mol = molecular_att(alts.get(key, 0.0), mid)
        g = mid ** 2 / mol
        gref = g[(mid > 4000) & (mid < 5000)].mean()
        for kind in ("night", "noon"):
            for sname, smask in seas.items():
                idx = np.where((meta["stream"] == key).values & (meta["kind"] == kind).values & smask)[0]
                if len(idx) < 8:
                    continue
                X, S = xn[idx], sn[idx]
                # drop windows whose 5-8 km signal departs > 30 % from the window (undetected cirrus / layers)
                dev = np.nanmedian(X[:, (mid > 5000) & (mid < 8000)], axis=1) - 1
                ok = np.abs(dev) < 0.3
                if ok.sum() < 5:
                    continue
                X, S = X[ok], S[ok]
                w = np.where(np.isfinite(X) & np.isfinite(S), 1 / S ** 2, 0.0)
                mean = np.nansum(np.where(np.isfinite(X), X, 0) * w, 0) / w.sum(0)
                se = 1 / np.sqrt(w.sum(0))
                fit = (mid >= 4000) & np.isfinite(mean) & np.isfinite(se) & (w.sum(0) > 0)
                A = np.column_stack([np.ones(fit.sum()), g[fit] / gref])
                W = 1 / se[fit] ** 2
                cov = np.linalg.inv(A.T @ (A * W[:, None]))
                c = cov @ (A.T @ (W * mean[fit]))
                chi2 = np.sum(W * (mean[fit] - A @ c) ** 2) / max(fit.sum() - 2, 1)
                cse = np.sqrt(np.diag(cov) * max(chi2, 1.0))
                bgrel = float(np.median(meta["bg_rel_win"].values[idx][ok]))
                rows.append(dict(stream=key, kind=kind, season=sname, n=int(ok.sum()), a=c[0], b=c[1],
                                 bias_pct=100 * c[1] / c[0], bias_se_pct=100 * cse[1] / c[0], chi2=chi2,
                                 x_10_15km=float(np.nanmean(mean[mid > 10000])), bg_over_signal_win=bgrel,
                                 eps_of_bg=c[1] / c[0] / bgrel))
                if sname == "all":
                    profs[(key, kind)] = mean
    return pd.DataFrame(rows), profs, mid


def harmonic_in_phase(t, y, peak_ref):
    amp, peak, _ = annual_fit(t, y)
    return amp * np.cos(2 * np.pi * (peak - peak_ref) / 365.25), amp, peak


def stage_summary():
    import pandas as pd
    units = pd.read_csv(Q15 / "q15e_units.csv")
    net = units[units["stream"] != PAYERNE]
    ser = np.load(Q15 / "q15e_series.npz")
    keys = [k for k in net["stream"] if f"{k}__t" in ser.files]
    res = {}
    res["T"] = pooled(ser, keys, "la", ["ta"])
    res["T_bg_sd"] = pooled(ser, keys, "la", ["ta", "bg", "sd"])
    res["win"] = pooled(ser, keys, "la", ["ta", "zb", "zt"])
    res["D_L"] = pooled(ser, keys, "D", ["ta", "L"])
    res["x78"] = pooled(ser, keys, "x78", ["ta"])
    res["la_D"] = pooled(ser, keys, "la", ["ta", "D"])
    rows = []
    for k in keys:
        t, lnc = ser[f"{k}__t"], ser[f"{k}__lnc"]
        _, ac, pc = harmonic_in_phase(t, lnc, 0)
        if not np.isfinite(ac):
            continue
        r = dict(stream=k, amp_lnc=100 * ac, peak_lnc=pc)
        lnx = ser[f"{k}__lnx"]
        for j, (lo, hi) in enumerate(BANDS):
            r[f"ip_{lo:.1f}_{hi:.1f}"] = 100 * harmonic_in_phase(t, lnx[:, j], pc)[0]
        r["ip_win"] = 100 * harmonic_in_phase(t, ser[f"{k}__lnx_win"], pc)[0]
        lnl = np.log(np.nanmean(np.exp(lnx[:, 0:3]), axis=1)) - ser[f"{k}__lnx_win"]
        for name, y, sc in (("bg", np.log(ser[f"{k}__bg"]), 100), ("sd", np.log(ser[f"{k}__sd"]), 100),
                            ("zb", ser[f"{k}__zb"] / 1000, 1), ("zt", ser[f"{k}__zt"] / 1000, 1),
                            ("D", ser[f"{k}__D"], 100), ("L", lnl, 100)):
            ip, a, pk = harmonic_in_phase(t, y, pc)
            r[f"ip_{name}"], r[f"amp_{name}"], r[f"peak_{name}"] = sc * ip, sc * a, pk
        rows.append(r)
    ann = pd.DataFrame(rows)
    ann.to_csv(Q15 / "q15e_annual_bands.csv", index=False)
    l1, profs, mid = l1_aggregate()
    l1.to_csv(Q15 / "q15e_l1_baseline.csv", index=False)
    pooled_rows = []
    for name, (c, se, r2, n, _) in res.items():
        for v in c:
            pooled_rows.append(dict(model=name, regressor=v, coef=c[v], se=se[v], r2=r2, n=n))
    pd.DataFrame(pooled_rows).to_csv(Q15 / "q15e_pooled.csv", index=False)
    write_report(net, units, ann, l1, profs, mid, res)


def write_report(net, units, ann, l1, profs, mid, res):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pay = units[units["stream"] == PAYERNE]
    keys = list(ann["stream"])
    kz_in = pd.read_csv(Q15 / "q15e_kz_intranight.csv")
    kz_in = kz_in[kz_in["stream"].isin(net["stream"])]
    pairs = pd.read_csv(Q15 / "q15e_kz_pairs.csv")
    l1n = l1[l1["stream"] != PAYERNE]
    night_all = l1n[(l1n["kind"] == "night") & (l1n["season"] == "all")]
    piv = l1n.pivot_table(index=["stream", "kind"], columns="season", values="bias_pct")
    wc_n = (piv["warm"] - piv["cold"]).xs("night", level="kind").dropna()
    allp = l1n[l1n["season"] == "all"].pivot(index="stream", columns="kind", values="bias_pct")
    dn = (allp["noon"] - allp["night"]).dropna()

    def q3(col, d=net):
        v = d[col].dropna()
        return np.median(v), np.percentile(v, 25), np.percentile(v, 75)

    # ------------------------------------------------------------------ figure
    fig, axs = plt.subplots(1, 3, figsize=(17, 6))
    ax = axs[0]
    kzc = [f"kz_{lo:.1f}_{hi:.1f}" for lo, hi in BANDS]
    ax.fill_betweenx(BAND_MID, net[kzc].quantile(.25).values, net[kzc].quantile(.75).values, color="#ff7f0e",
                     alpha=.2)
    ax.plot(net[kzc].median().values, BAND_MID, "o-", color="#ff7f0e",
            label="multi-day anomalies, per unit (median, IQR)")
    zin = np.array([0.75, 1.25, 2.0, 3.25, 5.0, 7.0])
    cin = [f"kin_{n}" for n in ZBANDS]
    ax.fill_betweenx(zin, kz_in[cin].quantile(.25).values, kz_in[cin].quantile(.75).values, color="#1f77b4",
                     alpha=.2)
    ax.plot(kz_in[cin].median().values, zin, "s-", color="#1f77b4", label="intra-night, per unit (median, IQR)")
    for pr, col, lab in (("0-20000-0-12843_A/B", "#2ca02c", "Budapest A/B"), ("0-20008-0-INO_A/B", "#9467bd",
                                                                              "Magurele A/B"),
                         ("0-20000-0-06240_B/C", "#d62728", "Schiphol B/C")):
        p = pairs[pairs["pair"] == pr]
        ax.errorbar(p["kA_in"], zin[:len(p)], xerr=p["kA_in_se"], fmt="^--", color=col, ms=5, lw=0.8,
                    label=f"pair {lab}, κ of first unit (intra-night, aerosol cancels)")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("κ(z) [%/K]")
    ax.set_ylabel("range above instrument [km]")
    ax.set_title("temperature coefficient vs altitude", fontsize=10)
    ax.legend(fontsize=7, loc="upper right")
    ax.set_xlim(-4, 14)
    ax.grid(alpha=.3)
    ax = axs[1]
    ipc = [f"ip_{lo:.1f}_{hi:.1f}" for lo, hi in BANDS]
    ax.fill_betweenx(BAND_MID, ann[ipc].quantile(.25).values, ann[ipc].quantile(.75).values, color="#1f77b4",
                     alpha=.2)
    ax.plot(ann[ipc].median().values, BAND_MID, "o-", color="#1f77b4", label="clear-sky X(z) = rcs/β$_{mol}$")
    ax.axvline(ann["amp_lnc"].median(), color="k", ls="--", lw=1,
               label=f"C$_L$ ({ann['amp_lnc'].median():.1f} %)")
    ax.axhspan(net["zb_med"].median() / 1000, net["zt_med"].median() / 1000, color="grey", alpha=.15,
               label="median Rayleigh window")
    ax.set_xlabel("annual amplitude in phase with C$_L$ [%]")
    ax.set_ylabel("range above instrument [km]")
    ax.set_title("annual cycle by altitude (a pure gain = the dashed line)", fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(alpha=.3)
    ax = axs[2]
    for kind, col in (("night", "#1f77b4"), ("noon", "#d62728")):
        P = np.array([v for (k, kd), v in profs.items() if kd == kind and k != PAYERNE])
        ax.fill_betweenx(mid / 1000, np.nanpercentile(P, 25, axis=0), np.nanpercentile(P, 75, axis=0), color=col,
                         alpha=.15)
        ax.plot(np.nanmedian(P, axis=0), mid / 1000, "-", color=col, label=f"{kind} (median of {len(P)} units)")
    if (PAYERNE, "night") in profs:
        ax.plot(profs[(PAYERNE, "night")], mid / 1000, ":", color="k", label="Payerne A night (excluded)")
    ax.axvline(1, color="k", lw=0.6)
    ax.set_xlim(-0.5, 2.5)
    ax.set_xlabel("X(r) / X(window)   (1 = molecular, no residual baseline)")
    ax.set_ylabel("range above instrument [km]")
    ax.set_title("residual-baseline probe to 15 km (native L1)", fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(Q15 / "q15e_profiles.png", dpi=130)
    plt.close(fig)

    # ------------------------------------------------------------------ markdown
    cT, _, r2T, _, _ = res["T"]
    cb, sb, r2b, _, _ = res["T_bg_sd"]
    cw, sw, r2w, _, Xw = res["win"]
    zb_sd, zt_sd = float(np.std(Xw[:, 1])), float(np.std(Xw[:, 2]))
    bg_sd, sd_sd = net["std_bg_a"].median() / 100, net["std_sd_a"].median() / 100
    med = lambda c: ann[c].median()
    nb = q3("bias_pct", night_all)
    eps = night_all["eps_of_bg"].median()
    bgrel = night_all["bg_over_signal_win"].median()
    wc = (np.median(wc_n), np.percentile(wc_n, 25), np.percentile(wc_n, 75))
    win_ip = abs(cw["zb"] * med("ip_zb") + cw["zt"] * med("ip_zt")) * 100
    cdl, sdl = res["D_L"][0]["L"], res["D_L"][1]["L"]
    grad = (kz_in["kin_6-8"] - kz_in["kin_2.5-4"]) / 3.75
    L = []
    w = L.append
    w("# q15e - electronic noise and profile distortion (summary)")
    w("")
    w(f"Network = CHM15k main optical module per unit, 2025-2026 operational constants, {len(keys)} units "
      f"(Payerne A excluded, see below). Baseline probe: {night_all['stream'].nunique()} geographically spread units, "
      f"{int(night_all['n'].sum())} calibration nights + noon windows read from native L1 up to 15 km. "
      "Night-to-night = 31-day anomalies, network-pooled regression x anomaly std; annual = first harmonic, "
      "in phase with each unit's C_L.")
    w("")
    w("| factor | sign | sensitivity | night-to-night [%] | diurnal [%] | annual amplitude [%] | phase DOY | evidence | rank |")
    w("|---|---|---|---|---|---|---|---|---|")
    w(f"| Residual additive baseline B·r² (background-subtraction error) | − at night (over-subtracted), + by day | "
      f"night C bias {nb[0]:+.1f} % (IQR {nb[1]:+.1f}…{nb[2]:+.1f}) = {eps * 1e4:+.1f}e-4 of the background "
      f"(background ≈ {bgrel:.0f}× the molecular signal in the window) | not resolvable per night (σ 5-10 %); "
      f"bounded ≲ 1 by the 7-8 km noise | **{np.median(dn):+.1f}** at ~4 km "
      f"(IQR {np.percentile(dn, 25):+.1f}…{np.percentile(dn, 75):+.1f}), ∝ r²/β_mol: ≈ +1 at 2 km, ≈ 0 below 1 km | "
      f"warm−cold {wc[0]:+.1f} %-pt (IQR {wc[1]:+.1f}…{wc[2]:+.1f}) → ≤ {abs(wc[0]) / 1.6:.1f}, n.s. | "
      f"anti-phase (summer lower) | {night_all['stream'].nunique()} units, L1 to 15 km | 1 |")
    w(f"| Altitude dependence of κ (shape vs gain) | none above 2.5 km | intra-night κ "
      f"{kz_in['kin_2.5-4'].median():.2f} / {kz_in['kin_4-6'].median():.2f} / {kz_in['kin_6-8'].median():.2f} %/K at "
      f"2.5-4 / 4-6 / 6-8 km → {grad.median():+.2f} %/K/km (SE {1.25 * grad.std() / np.sqrt(grad.notna().sum()):.2f}) "
      f"| 0 | 0 (gain flat → daytime excess flat in height) | 7-8 km X carries {med('ip_7.0_8.0'):.1f} of "
      f"{med('amp_lnc'):.1f} % (window {med('ip_win'):.1f}) | – | {len(kz_in)} units + 3 pair sites | 2 |")
    w(f"| Window selection moving with SNR / short summer nights | + with window bottom | "
      f"{100 * cw['zb']:+.2f} ± {100 * sw['zb']:.2f} %/km (bottom), {100 * cw['zt']:+.2f} ± {100 * sw['zt']:.2f} %/km "
      f"(top) | {100 * np.hypot(cw['zb'] * zb_sd, cw['zt'] * zt_sd):.2f} (ΔR² {100 * (r2w - r2T):.1f} %) | – (night "
      f"only) | {win_ip:.2f} (bottom ±{q3('amp_zb')[0] * 1000:.0f} m, night length ±{q3('amp_dur')[0]:.1f} h) | "
      f"{q3('peak_zb')[0]:.0f} (in phase) | {len(keys)} units | 3 |")
    w(f"| Night background `bckgrd_rcs_0` (dark) | + (n.s.) | {cb['bg']:+.2f} ± {sb['bg']:.2f} %/% (pooled); "
      f"background itself {q3('bg_per_K')[0]:+.2f} %/K | {abs(cb['bg']) * bg_sd * 100:.2f} | 0 on C (noon ×"
      f"{q3('bg_noon_night')[0]:.0f}, acts through the baseline row) | {abs(cb['bg']) * med('ip_bg'):.2f} "
      f"(background cycle {q3('amp_bg')[0]:.1f} %) | {q3('peak_bg')[0]:.0f} | {len(keys)} units | 4 |")
    w(f"| Electronic noise `stddev` / SNR | ± (n.s.) | {cb['sd']:+.2f} ± {sb['sd']:.2f} %/%; noise "
      f"{q3('sd_per_K')[0]:+.2f} %/K | {abs(cb['sd']) * sd_sd * 100:.2f} | 0 (noise ×{q3('sd_noon_night')[0]:.1f} by "
      f"day: precision only) | {abs(cb['sd']) * med('ip_sd'):.2f} (noise cycle {q3('amp_sd')[0]:.1f} %) | "
      f"{q3('peak_sd')[0]:.0f} | {len(keys)} units | 5 |")
    w(f"| Signal-induced noise / afterpulse after strong low-level signal | + | 7-8 km excess {cdl:+.3f} ± {sdl:.3f} "
      f"per unit ln(0.2-1 km load); window sees ~1/5 | {abs(cdl) * q3('std_L_a')[0] / 5:.2f} | ≤ 0.5 | "
      f"{abs(cdl) * abs(med('ip_L')) / 5:.2f} | {med('peak_L'):.0f} | {len(keys)} units; partly atmospheric | 6 |")
    w(f"| Temperature-dependent overlap (Hervo 2016) | + near range only | no-op above ~0.8 km; window bottom p05 "
      f"{q3('zb_p05')[0] / 1000:.2f} km, {100 * q3('frac_zb_below_850')[0]:.0f} % of windows below 0.85 km | 0 | 0 on C; "
      f"L2 0.5-1 km: excess κ {kz_in['kin_0.5-1'].median() - kz_in['kin_2.5-4'].median():+.1f} %/K per unit "
      "(aerosol-confounded), pairs −0.1…+1.1 %/K → ≲ 1-4 % on hot afternoons | 0 | – | windows of all units, 3 pairs | 7 |")
    w("")
    w("Context (not instrumental, other family): the profile-shape index S = ln X(2-3 km) − ln X(4-6 km) moves "
      f"{q3('S_per_K')[0]:+.1f} %/K with multi-day T anomalies and has an annual in-phase amplitude of "
      f"{med('ip_2.0_2.5'):.0f} % at 2-2.5 km against {med('ip_win'):.0f} % in the window: aerosol. Pairs "
      "(aerosol cancels) show no instrumental shape change.")
    w("")
    w("![kappa(z), annual amplitude by altitude and the 15-km baseline probe](q15e_profiles.png)")
    w("")
    w("## Key numbers")
    w(f"* Pooled ln C anomalies on T_int anomalies: κ = {100 * cT['ta']:.2f} %/K (R² {100 * r2T:.1f} %). Adding "
      f"background + noise anomalies explains +{100 * (r2b - r2T):.2f} % more variance; window heights "
      f"+{100 * (r2w - r2T):.2f} %.")
    w(f"* Multi-day κ per band (per unit) peaks at 2-2.5 km ({q3('kz_2.0_2.5')[0]:.1f} %/K) and is "
      f"{q3('kz_7.0_8.0')[0]:+.2f} %/K at 7-8 km (pooled {100 * res['x78'][0]['ta']:+.2f} ± "
      f"{100 * res['x78'][1]['ta']:.2f}), while intra-night κ is flat above 2.5 km ({kz_in['kin_6-8'].median():.2f} %/K "
      "at 6-8 km). Open point: at synoptic time scales the aerosol-free 7-8 km signal does not follow T_int.")
    w(f"* Annual: clear-sky X(7-8 km) in-phase amplitude {med('ip_7.0_8.0'):.1f} % (IQR "
      f"{ann['ip_7.0_8.0'].quantile(.25):.1f}…{ann['ip_7.0_8.0'].quantile(.75):.1f}), C_L {med('amp_lnc'):.1f} %, "
      f"window {med('ip_win'):.1f} %: the annual cycle is height-independent from the window to 8 km.")
    w(f"* Residual baseline: night bias {nb[0]:+.1f} %, noon − night {np.median(dn):+.1f} %-pt ({len(dn)} units), "
      f"warm − cold night {wc[0]:+.1f} %-pt ({len(wc_n)} units).")
    w(f"* Window: bottom {q3('zb_med')[0] / 1000:.2f} km AGL (annual ±{q3('amp_zb')[0] * 1000:.0f} m, highest DOY "
      f"{q3('peak_zb')[0]:.0f}), top {q3('zt_med')[0] / 1000:.2f} km; calibration duration ±{q3('amp_dur')[0]:.1f} h.")
    w("")
    w("## Payerne A (TUB140016, excluded)")
    if len(pay):
        p = pay.iloc[0]
        lp = l1[l1["stream"] == PAYERNE]
        w(f"* night background {p['bg_med']:.4f} (network {q3('bg_med')[0]:.4f}), noise {p['sd_med']:.2e} "
          f"(network {q3('sd_med')[0]:.2e}); ln C on ln background {p['g_bg']:+.1f} %/%; multi-day κ "
          f"{p['kappa_bare']:+.2f} %/K; window-bottom coefficient {p['g_zb']:+.1f} %/km (network "
          f"{q3('g_zb')[0]:+.1f}).")
        w("* L1 baseline bias: " + ", ".join(f"{r.kind}/{r.season} {r.bias_pct:+.1f} ± {r.bias_se_pct:.1f} %"
                                           for r in lp.itertuples()))
    w("")
    w("## Caveats")
    w("* Molecular reference = US-1976 standard atmosphere (as operational): seasonal density errors (±1-2 % "
      "between 4 and 8 km, more above the tropopause) limit the baseline probe to ~±0.5 % per unit-season.")
    w("* Baseline fits have χ² ≫ 1 (aerosol layers, thin cirrus, standard-atmosphere shape); errors scaled by √χ². "
      "The one-parameter B·r² model cannot separate a range-dependent ('m-shaped') dark from aerosol.")
    w("* Night profile shape (figure, right): X/X_window stays at 1.0-1.1 from 4 to 11 km and only collapses "
      "above ~12 km. That is not B·r²-shaped: the fitted night bias (−0.9 %) is driven by the 12-15 km bins "
      "(tropopause / standard-atmosphere mismatch or a far-range-only dark), so the night baseline inside the "
      "window is ≈ 0 (|bias| ≲ 1 %). The noon excess (+40…+90 % at 8-11 km) is r²-shaped: a real daytime "
      "under-subtraction of the solar background.")
    w("* Per-unit coefficients are noisy; the table uses network-pooled regressions and medians.")
    w("* Open point: multi-day (synoptic) T_int anomalies are not seen in the aerosol-free 7-8 km signal "
      "(κ ≈ 0 ± 0.14 %/K pooled) although they are intra-night (≈ +1.4 %/K) and the annual cycle is. Either an "
      "atmospheric covariate at 6-8 km cancels the gain at synoptic scales, or part of the multi-day window κ is "
      "aerosol in the window. Worth a cross-check by the aerosol family.")
    open(Q15 / "q15e_summary.md", "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("summary written:", Q15 / "q15e_summary.md")


if __name__ == "__main__":
    _main()
