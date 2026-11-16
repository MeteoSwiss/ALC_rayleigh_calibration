# -*- coding: utf-8 -*-
"""q15d: atmosphere and Rayleigh-method factors behind the CHM15k constant's annual cycle (report 19, q15, agent d).

The Rayleigh constant is the slope of the clear-sky signal against the MODELLED molecular signal in a window
zb..zt (linregress with a free intercept, 1/r^2 weighting, molecular_source "standard" = US-1976 by ASL
altitude, not scaled to the station), multiplied by the inverse two-way aerosol transmission below the window
from a Klett inversion with LR = 52 sr. Anything that makes the real atmosphere differ from that model biases
C_L without the instrument changing. Per night this script predicts, on the night's own window:

  F1  aerosol inside the window:   slope of (beta_m + beta_a) T_a,win^2 vs beta_m   (CAMS beta_a at 1064 nm)
  F2  aerosol extinction below:    the Klett transmission factor the operational chain applied, re-derived
                                   from the cached clear-sky rcs (tau_K), vs AERONET AOD_1064 and CAMS
  F3  molecular model:             slope of beta_m,CAMS T_m,CAMS^2 vs beta_m,US76 T_m,US76^2
  F4  sampling / selection:        window centre and night length (short-term regressions x annual swing),
                                   CAMS aerosol of accepted vs all nights
plus the lidar's own profile shape inside the window (ln R gradient), which tests the CAMS aerosol cycle.

Each predicted ln-bias series is reduced like the observed ln C: annual harmonic (+ linear trend) and
anomalies against a 31-day running median. Co-located pairs test whether the residual cycle is common.

Inputs: q15d_cams_extract.py output, common.load_constants/load_cache, AERONET daily files (q5).
Run:    python validation/chm15k_stability/q15d_atmosphere.py
"""
from __future__ import annotations

import csv
import json
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (CENSUS, OUT, PAIR_SITES, RANGE_EDGES, RANGE_MID, chm15k_streams, load_cache,  # noqa: E402
                    load_constants, stream_key)
from q1_temperature import ANOM_WINDOW, night_hk, robust_mask, running_median  # noqa: E402
from q5_aeronet import read_aeronet_daily  # noqa: E402
from calibration.rayleigh.atmosphere import calculate_molecular_properties, load_standard_atmosphere  # noqa: E402

warnings.filterwarnings("ignore")
Q15 = OUT / "q15"
CAMS_DIR = Q15 / "q15d_cams"
SCRATCH = Path("C:/Users/hervo/AppData/Local/Temp/claude/C--Users-hervo-OneDrive-Documents-ALC-rayleigh-calibration/"
               "e6120125-e8c9-4e0b-8352-26777f3067de/scratchpad")
AERONET_DIR = SCRATCH / "aeronet_raw"
PAYERNE = "0-20000-0-06610_A"
EPOCH = datetime(1970, 1, 1)
WL = 1064e-9
LR_KLETT = 52.0          # options.json LRaer (operational Klett)
LR_CAMS = 50.0           # extinction-to-backscatter ratio used to turn CAMS beta_a into extinction
S_MOL = 8 * np.pi / 3
DZ = 30.0                # fine grid for the CAMS-based window simulation [m]
W_YEAR = 2 * np.pi / 365.25


# ---------------------------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------------------------
def beta_mol(T, p):
    """Bucholtz molecular backscatter at 1064 nm, identical to the calibration chain."""
    return calculate_molecular_properties(T, p, np.zeros_like(T), WL).beta_mol


def two_way(ext, dz, start=0):
    """Two-way transmission from index `start` upward on a uniform grid (1 below start)."""
    tau = np.zeros_like(ext)
    tau[start + 1:] = np.cumsum(0.5 * (ext[start:-1] + ext[start + 1:])) * dz
    return np.exp(-2 * tau)


def fit_slope(x, y):
    """linregress slope (free intercept), as the molecular-window search does."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5:
        return np.nan
    x, y = x[ok], y[ok]
    xm = x - x.mean()
    return float(np.sum(xm * (y - y.mean())) / np.sum(xm * xm))


def harmonic(t, y):
    """Annual harmonic + linear trend by least squares. Returns (amp, peak DOY, trend per yr, fitted)."""
    ok = np.isfinite(t) & np.isfinite(y)
    if ok.sum() < 40 or (t[ok].max() - t[ok].min()) < 300:
        return np.nan, np.nan, np.nan, None
    doy = np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in t[ok]], float)
    tt = (t[ok] - t[ok].mean()) / 365.25
    X = np.column_stack([np.ones(ok.sum()), np.cos(W_YEAR * doy), np.sin(W_YEAR * doy), tt])
    coef, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    a, b = coef[1], coef[2]
    amp = float(np.hypot(a, b))
    peak = float((np.arctan2(b, a) / W_YEAR) % 365.25)
    return amp, peak, float(coef[3]), (ok, X @ coef)


def anomaly(t, y):
    return y - running_median(t, y, ANOM_WINDOW)


def doy_of(t):
    return np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in t], float)


# ---------------------------------------------------------------------------------------------
# CAMS columns (from q15d_cams_extract.py)
# ---------------------------------------------------------------------------------------------
def load_cams():
    files = sorted(CAMS_DIR.glob("cams_*.npz"))
    acc = {}
    keys = None
    for f in files:
        with np.load(f) as z:
            keys = list(z["keys"])
            for v in ("time", "aer1064", "T", "alt", "p"):
                acc.setdefault(v, []).append(z[v] if v == "time" else z[v])
    time = np.concatenate(acc["time"])
    order = np.argsort(time, kind="stable")
    time = time[order]
    keep = np.r_[True, np.diff(time) > 1e-6]          # daily files overlap: keep first occurrence
    sel = order[keep]
    out = {"keys": keys, "time": time[keep]}
    for v in ("aer1064", "T", "alt", "p"):
        out[v] = np.concatenate(acc[v], axis=1)[:, sel, :]
    return out


def cams_aerosol_from_attenuated(alt, att, T, p):
    """Aerosol backscatter from the CAMS 'backscatter from ground' field.

    The field is the ATTENUATED TOTAL backscatter seen from the model surface (molecular + aerosol,
    two-way attenuated): at 20-25 km it equals beta_mol * T_mol^2 within a few % (checked on 4 sites).
    Forward retrieval from the surface with LR_CAMS for the aerosol extinction; negative residuals are
    kept (no positive bias from clipping) but do not attenuate.
    """
    bm = beta_mol(T, p)
    ba = np.zeros_like(att)
    tau = 0.0
    for i in range(len(alt)):
        if i > 0:
            dz = alt[i] - alt[i - 1]
            tau += dz * 0.5 * (S_MOL * (bm[i] + bm[i - 1]) + LR_CAMS * (max(ba[i - 1], 0.0) + max(ba[i - 1], 0.0)))
        ba[i] = att[i] / np.exp(-2 * tau) - bm[i]
    return ba


def cams_night(cams, ip, t0, t1):
    """Night-mean CAMS profile (ascending altitude) at stream index ip, or None.

    Returns (alt, beta_aer, T, p); beta_aer is the de-attenuated aerosol backscatter at 1064 nm.
    """
    m = (cams["time"] >= t0 - 0.07) & (cams["time"] <= t1 + 0.07)
    if not m.any():
        return None
    sl = lambda v: np.nanmean(cams[v][ip, m, :], axis=0)[::-1].astype(float)   # noqa: E731
    alt, att, T, p = sl("alt"), sl("aer1064"), sl("T"), sl("p")
    if not np.all(np.isfinite(alt)) or not np.all(np.isfinite(att)):
        return None
    return alt, cams_aerosol_from_attenuated(alt, att, T, p), T, p


# ---------------------------------------------------------------------------------------------
# Per-night predicted biases from CAMS
# ---------------------------------------------------------------------------------------------
class StdAtm:
    def __init__(self):
        self.z = np.arange(0, 12001, DZ)
        a = load_standard_atmosphere(None, self.z)
        self.bm = beta_mol(a.temperature, a.pressure)


STD = None


def night_cams_terms(prof, alt_stn, zb, zt):
    """Predicted ln-biases of C_L for one night (CAMS profile, window zb..zt ASL)."""
    alt, aer, T, p = prof
    z = np.arange(alt_stn + DZ, max(zt, 6000) + DZ, DZ)
    r = z - alt_stn
    lnp = np.interp(z, alt, np.log(p))
    Tz = np.interp(z, alt, T)
    ba = np.interp(z, alt, aer, left=aer[0])
    bm_c = beta_mol(Tz, np.exp(lnp))
    bm_s = np.interp(z, STD.z, STD.bm)
    w = (z >= zb) & (z <= zt)
    if w.sum() < 5:
        return None
    iw = int(np.argmax(w))
    tm_s = two_way(S_MOL * bm_s, DZ)
    tm_c = two_way(S_MOL * bm_c, DZ)
    ta_win = two_way(LR_CAMS * np.clip(ba, 0, None), DZ, start=iw)          # aerosol attenuation INSIDE the window only
    x = bm_s * tm_s / r ** 2
    # F1 aerosol in window (US76 molecular, as operational)
    b_aer = fit_slope(x[w], ((bm_s + ba) * tm_s * ta_win / r ** 2)[w]) - 1.0
    ratio_aer = float(np.mean(ba[w] / bm_s[w]))
    # F3 molecular: real (CAMS) molecular signal vs the US-1976 model
    b_mol = fit_slope(x[w], (bm_c * tm_c / r ** 2)[w]) - 1.0
    ratio_mol = float(np.mean(bm_c[w] / bm_s[w])) - 1.0
    # in-window profile shape the aerosol would imprint (d ln R / dz per km), CAMS
    lnR = np.log((bm_s + ba) * ta_win / bm_s)[w]
    g_cams = np.polyfit(z[w] / 1000.0, lnR, 1)[0]
    # F2 aerosol optical depth below the window (CAMS, LR_CAMS) and total column to 11 km
    tau_below = LR_CAMS * np.trapezoid(np.clip(ba[:iw], 0, None), z[:iw]) if iw > 1 else 0.0
    return dict(b_aer=np.log1p(b_aer), r_aer=ratio_aer, b_mol=np.log1p(b_mol), r_mol=ratio_mol,
                g_cams=g_cams, tau_below_cams=tau_below)


def cams_column_aod(prof, alt_stn):
    alt, aer, T, p = prof
    ok = alt >= alt_stn
    return LR_CAMS * np.trapezoid(np.clip(aer[ok], 0, None), alt[ok])


# ---------------------------------------------------------------------------------------------
# Lidar-side terms from the cached clear-sky rcs
# ---------------------------------------------------------------------------------------------
def night_rcs(cache, t0, t1):
    """Night median of the 10-min clear-sky rcs profiles inside [t0, t1] (38 range bins)."""
    tb = (cache["day"][:, None] + (np.arange(144)[None, :] + 0.5) / 144).ravel()
    rcs = cache["rcs"].reshape(-1, cache["rcs"].shape[-1])
    ncl = cache["n_clear"].ravel()
    order = np.argsort(tb)
    tb, rcs, ncl = tb[order], rcs[order], ncl[order]
    i0, i1 = np.searchsorted(tb, t0), np.searchsorted(tb, t1, side="right")
    out = np.full((len(t0), rcs.shape[1]), np.nan)
    nb = np.zeros(len(t0), int)
    for k in range(len(t0)):
        seg, n = rcs[i0[k]:i1[k]], ncl[i0[k]:i1[k]]
        seg = seg[n > 0]
        if len(seg) >= 3:
            out[k] = np.nanmedian(seg, axis=0)
            nb[k] = len(seg)
    return out, nb


def lidar_terms(rcs, c, alt_stn, tilt, zb, zt):
    """Klett-like forward aerosol retrieval below the window (LR 52) and the in-window ln R gradient."""
    z = alt_stn + RANGE_MID * np.cos(np.deg2rad(tilt))
    dzb = np.diff(RANGE_EDGES) * np.cos(np.deg2rad(tilt))
    bm = np.interp(z, STD.z, STD.bm)
    batt = rcs / c
    w = (z >= zb) & (z <= zt) & np.isfinite(batt)
    if w.sum() < 3:
        return None
    # forward integration from the ground: tau accumulated with the extinction of the bins below
    tau_a = 0.0
    tau_m = 0.0
    ba = np.zeros_like(z)
    for i in range(len(z)):
        if z[i] >= zb:
            break
        t2 = np.exp(-2 * (tau_a + tau_m))
        ba[i] = max(batt[i] / t2 - bm[i], 0.0) if np.isfinite(batt[i]) else 0.0
        tau_a += LR_KLETT * ba[i] * dzb[i]
        tau_m += S_MOL * bm[i] * dzb[i]
    # in-window shape: ln of the normalised ratio, d/dz per km
    tm2 = np.exp(-2 * np.cumsum(S_MOL * bm * dzb))
    R = batt * np.exp(2 * tau_a) / (bm * tm2)
    g = np.polyfit(z[w] / 1000.0, np.log(np.clip(R[w], 1e-3, None)), 1)[0] if w.sum() >= 3 else np.nan
    return dict(tau_klett=tau_a, g_obs=g, R_win=float(np.nanmedian(R[w]))), R, z


# ---------------------------------------------------------------------------------------------
# AERONET
# ---------------------------------------------------------------------------------------------
def aeronet_by_stream():
    f = SCRATCH / "aeronet_match.json"
    if not f.exists():
        return {}
    out = {}
    for key, site, ae_site, dist in json.load(open(f, encoding="utf-8")):
        p = AERONET_DIR / f"AOD15_daily_{ae_site}.txt"
        if p.exists():
            d = read_aeronet_daily(p)
            if len(d) >= 50:
                out[key] = d
    return out


def night_aod(aer, t0, t1):
    vals = [aer[x] for x in {(EPOCH + timedelta(days=float(t0))).date(),
                             (EPOCH + timedelta(days=float(t1))).date()} if x in aer]
    return np.mean([v[0] for v in vals]) if vals else np.nan


# ---------------------------------------------------------------------------------------------
# Main per-unit loop
# ---------------------------------------------------------------------------------------------
def kappa_table():
    k = {}
    for r in csv.DictReader(open(OUT / "q1_slopes.csv", encoding="utf-8")):
        if r["arm"] == "A" and r["tvar"] == "temp_int" and r["b2"] not in ("", "nan"):
            k[(r["stream"], r["era"])] = float(r["b2"]) / 100
    return k


def process_unit(wmo, ident, site, alt_stn, ip, cams, kap, aeronet):
    key = stream_key(wmo, ident)
    cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
    if cons is None or cache is None:
        return None
    mod = max(set(cons["module"]), key=list(cons["module"]).count)
    m = (cons["module"] == mod) & (cons["c"] > 0) & np.isfinite(cons["zb"]) & (cons["t"] <= cache["day"][-1] + 1)
    if m.sum() < 100:
        return None
    lc = np.log(cons["c"][m])
    keep = robust_mask(lc)
    idx = np.where(m)[0][keep]
    t, t0, t1 = cons["t"][idx], cons["t0"][idx], cons["t1"][idx]
    c, zb, zt, u = cons["c"][idx], cons["zb"][idx], cons["zt"][idx], cons["u"][idx]
    lc = np.log(c)
    hk = night_hk(cache, t0, t1)
    T = hk[:, list(cache["hk_vars"]).index("temp_int")]
    k = kap.get((key, mod), np.nan)
    lc_T = lc - k * (T - np.nanmedian(T)) if np.isfinite(k) else np.full_like(lc, np.nan)
    rcs, nblk = night_rcs(cache, t0, t1)
    tilt = float(cache["tilt"]) if "tilt" in cache else 0.0
    n = len(t)
    Rprof = np.full((n, len(RANGE_MID)), np.nan)
    zbins = alt_stn + RANGE_MID * np.cos(np.deg2rad(tilt))
    rec = {v: np.full(n, np.nan) for v in ("b_aer", "r_aer", "b_mol", "r_mol", "g_cams", "tau_below_cams",
                                             "aod_cams", "tau_klett", "g_obs", "R_win", "aod_aeronet")}
    for i in range(n):
        if ip is not None:
            prof = cams_night(cams, ip, t0[i], t1[i])
            if prof is not None:
                r = night_cams_terms(prof, alt_stn, zb[i], zt[i])
                if r is not None:
                    for kk, v in r.items():
                        rec[kk][i] = v
                rec["aod_cams"][i] = cams_column_aod(prof, alt_stn)
        if np.all(np.isfinite(rcs[i][:20])) or np.isfinite(rcs[i]).sum() > 20:
            r = lidar_terms(rcs[i], c[i], alt_stn, tilt, zb[i], zt[i])
            if r is not None:
                r, Rprof[i], zbins = r
                for kk, v in r.items():
                    rec[kk][i] = v
        if key in aeronet:
            rec["aod_aeronet"][i] = night_aod(aeronet[key], t0[i], t1[i])
    # CAMS aerosol of ALL nights of the record (selection test): same fixed window = the unit's median window
    sel_all = None
    if ip is not None:
        zbm, ztm = float(np.median(zb)), float(np.median(zt))
        days = np.arange(np.floor(t.min()), np.floor(t.max()) + 1)
        rr = []
        for d in days:
            prof = cams_night(cams, ip, d - 0.125, d + 0.125)     # 21-03 UTC
            if prof is None:
                rr.append(np.nan)
                continue
            r = night_cams_terms(prof, alt_stn, zbm, ztm)
            rr.append(np.nan if r is None else r["b_aer"])
        sel_all = (days, np.array(rr), zbm, ztm)
    return dict(key=key, site=site, module=mod, t=t, t0=t0, t1=t1, c=c, lc=lc, lc_T=lc_T, kappa=k, T=T, zb=zb,
                zt=zt, u=u, nblk=nblk, rec=rec, sel_all=sel_all, Rprof=Rprof, zbins=zbins, alt=alt_stn)


def summarise_unit(U):
    """Per-unit metrics: observed and predicted annual harmonics, short-term regressions."""
    t, lc = U["t"], U["lc"]
    R = U["rec"]
    row = dict(stream=U["key"], site=U["site"], module=U["module"], n=len(t), kappa=100 * U["kappa"])
    a, p, tr, _ = harmonic(t, lc)
    row.update(amp_obs=100 * a, peak_obs=p, trend_obs=100 * tr)
    a, p, tr, _ = harmonic(t, U["lc_T"])
    row.update(amp_obsT=100 * a, peak_obsT=p)
    la = anomaly(t, lc)
    laT = anomaly(t, U["lc_T"])
    row["scat_obs"] = 100 * np.nanstd(la)
    zc = 0.5 * (U["zb"] + U["zt"]) / 1000.0
    dur = (U["t1"] - U["t0"]) * 24
    series = dict(aer=R["b_aer"], mol=R["b_mol"], kl=2 * R["tau_klett"], zc=zc, dur=dur, g_obs=R["g_obs"],
                  g_cams=R["g_cams"], tauc=R["tau_below_cams"], aodc=R["aod_cams"], aoda=R["aod_aeronet"],
                  urel=U["u"] / U["c"], zb=U["zb"] / 1000.0, nblk=U["nblk"].astype(float),
                  rmol=R["r_mol"], raer=R["r_aer"])
    for name, s in series.items():
        a, p, tr, _ = harmonic(t, s)
        scale = 100 if name in ("aer", "mol", "kl", "rmol", "raer") else 1
        row[f"amp_{name}"] = scale * a
        row[f"peak_{name}"] = p
        row[f"trend_{name}"] = scale * tr
        row[f"mean_{name}"] = scale * np.nanmean(s)
        sa = anomaly(t, s)
        row[f"std_{name}"] = scale * np.nanstd(sa)
        ok = np.isfinite(sa) & np.isfinite(laT)
        if ok.sum() >= 30 and np.nanstd(sa[ok]) > 0:
            b = np.polyfit(sa[ok], laT[ok], 1)[0]
            rr = np.corrcoef(sa[ok], laT[ok])[0, 1]
            row[f"slope_{name}"] = b                     # d lnC / d(series) on 31-day anomalies (kappa-normalised)
            row[f"r_{name}"] = rr
        else:
            row[f"slope_{name}"] = np.nan
            row[f"r_{name}"] = np.nan
    # signed projection of a predicted cycle onto the observed one: how much of the kappa-normalised
    # annual amplitude disappears when the predicted bias is removed
    for name, s in (("aer", R["b_aer"]), ("mol", R["b_mol"]), ("both", R["b_aer"] + R["b_mol"])):
        a0, _, _, _ = harmonic(t, U["lc_T"])
        a1, _, _, _ = harmonic(t, U["lc_T"] - s)
        row[f"amp_obsT_minus_{name}"] = 100 * a1
    # selection: CAMS aerosol bias of accepted vs all nights (same fixed window), JJA and DJF
    if U["sel_all"] is not None:
        days, ball, zbm, ztm = U["sel_all"]
        acc = np.isin(days, np.floor(t))
        dd = doy_of(days)
        jja = (dd >= 152) & (dd < 244)
        djf = (dd >= 335) | (dd < 60)
        for lab, sm in (("jja", jja), ("djf", djf)):
            row[f"aer_all_{lab}"] = 100 * np.nanmedian(ball[sm])
            row[f"aer_acc_{lab}"] = 100 * np.nanmedian(ball[sm & acc])
            row[f"frac_acc_{lab}"] = np.mean(acc[sm]) if sm.any() else np.nan
        a, p, _, _ = harmonic(days.astype(float), ball)
        row["amp_aer_allnights"] = 100 * a
        row["peak_aer_allnights"] = p
    # seasonal scatter: short summer nights noisier?
    dd = doy_of(t)
    for lab, sm in (("jja", (dd >= 152) & (dd < 244)), ("djf", (dd >= 335) | (dd < 60))):
        row[f"scat_{lab}"] = 100 * np.nanstd(la[sm]) if sm.sum() >= 10 else np.nan
        row[f"nblk_{lab}"] = np.nanmedian(U["nblk"][sm]) if sm.sum() else np.nan
    return row


def pair_table(units):
    """Co-located same-wavelength pairs: annual amplitude of same-night ln(C_A/C_B), raw and kappa-normalised."""
    rows = []
    for wmo, ids in PAIR_SITES.items():
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                A, B = units.get(f"{wmo}_{ids[i]}"), units.get(f"{wmo}_{ids[j]}")
                if A is None or B is None:
                    continue
                da, db = np.floor(A["t"]), np.floor(B["t"])
                common, ia, ib = np.intersect1d(da, db, return_indices=True)
                if len(common) < 40:
                    continue
                t = A["t"][ia]
                r = A["lc"][ia] - B["lc"][ib]
                rT = A["lc_T"][ia] - B["lc_T"][ib]
                aa, pa, _, _ = harmonic(t, A["lc"][ia])
                ab, pb, _, _ = harmonic(t, B["lc"][ib])
                aaT, _, _, _ = harmonic(t, A["lc_T"][ia])
                abT, _, _, _ = harmonic(t, B["lc_T"][ib])
                ar, pr, _, _ = harmonic(t, r)
                arT, prT, _, _ = harmonic(t, rT)
                aer_d = A["rec"]["b_aer"][ia] - B["rec"]["b_aer"][ib]
                xa, xb = anomaly(t, A["lc_T"][ia]), anomaly(t, B["lc_T"][ib])
                okc = np.isfinite(xa) & np.isfinite(xb)
                corr = float(np.corrcoef(xa[okc], xb[okc])[0, 1]) if okc.sum() > 10 else np.nan
                rows.append(dict(pair=f"{A['site']} {ids[i]}/{ids[j]}", n=len(common),
                                 amp_A=100 * aa, peak_A=pa, amp_B=100 * ab, peak_B=pb,
                                 amp_ratio=100 * ar, peak_ratio=pr,
                                 amp_A_T=100 * aaT, amp_B_T=100 * abT, amp_ratio_T=100 * arT,
                                 std_ratio_anom=100 * np.nanstd(anomaly(t, r)),
                                 std_A_anom=100 * np.nanstd(anomaly(t, A["lc"][ia])),
                                 corr_anom_T=corr,
                                 cams_aer_diff_std=100 * np.nanstd(aer_d)))
    return rows


def main():
    global STD
    STD = StdAtm()
    cams = load_cams()
    cams_idx = {k: i for i, k in enumerate(cams["keys"])}
    alt_of = {stream_key(s["wmo"], s["ident"]): float(s["alt"]) for s in json.load(open(CENSUS, encoding="utf-8"))}
    kap = kappa_table()
    aeronet = aeronet_by_stream()
    print(f"CAMS steps {len(cams['time'])}, AERONET streams {len(aeronet)}", flush=True)
    rows, units, nights = [], {}, []
    for wmo, ident, site, lat, lon in chm15k_streams():
        key = stream_key(wmo, ident)
        U = process_unit(wmo, ident, site, alt_of.get(key, 0.0), cams_idx.get(key), cams, kap, aeronet)
        if U is None:
            continue
        units[key] = U
        row = summarise_unit(U)
        row["alt"] = alt_of.get(key, np.nan)
        row["lat"], row["lon"] = lat, lon
        rows.append(row)
        for i in range(len(U["t"])):
            nights.append(dict(stream=key, t=U["t"][i], lnc=U["lc"][i], lnc_T=U["lc_T"][i], zb=U["zb"][i],
                               zt=U["zt"][i], T_int=U["T"][i], nblk=U["nblk"][i],
                               **{k: U["rec"][k][i] for k in U["rec"]}))
        print(f"{key:28s} {site[:16]:16s} n {row['n']:4d} amp obs {row['amp_obs']:5.1f}% (T-norm {row['amp_obsT']:5.1f})"
              f" | aer {row['amp_aer']:4.1f}% pk {row['peak_aer']:5.0f} mean {row['mean_aer']:4.1f}"
              f" | mol {row['amp_mol']:4.1f}% pk {row['peak_mol']:5.0f} | Klett 2tau {row['amp_kl']:4.1f}%"
              f" | zc {row['amp_zc']:.2f} km slope {row['slope_zc']:+.3f}/km", flush=True)
    with open(Q15 / "q15d_units.csv", "w", newline="", encoding="utf-8") as fh:
        keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 999)
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    with open(Q15 / "q15d_nights.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(nights[0].keys()))
        w.writeheader()
        w.writerows(nights)
    prs = pair_table(units)
    if prs:
        with open(Q15 / "q15d_pairs.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(prs[0].keys()))
            w.writeheader()
            w.writerows(prs)
    np.save(Q15 / "q15d_units.npy", np.array([dict(key=k, t=U["t"], lc=U["lc"], lc_T=U["lc_T"], zb=U["zb"],
                                                    zt=U["zt"], rec=U["rec"], Rprof=U["Rprof"], zbins=U["zbins"], alt=U["alt"],
                                                    T=U["T"], kappa=U["kappa"], site=U["site"])
                                               for k, U in units.items()],
                                             dtype=object), allow_pickle=True)
    print(f"-> {Q15 / 'q15d_units.csv'} ({len(rows)} units), q15d_nights.csv, q15d_pairs.csv")


if __name__ == "__main__":
    main()
