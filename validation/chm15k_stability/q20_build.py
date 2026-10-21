# -*- coding: utf-8 -*-
"""q20 stage 1 - per-night input table for the correction budget (report 19, q20).

For the 20 study units and Payerne A, every published nightly Rayleigh constant with a night start in
2025-01-01 .. 2026-07-31 (common.load_constants, calibration_method 0) is joined with:
  * housekeeping medians over the calibration window (HK cache, hk_cache.py): temp_int, temp_ext,
    temperature_optical_module, calibration_pulse, window_transmission, laser_life_time, and the number of
    clear 10-min profiles;
  * CAMS AOD by species at 00 UTC of the night (q17a_cams_points.npz, same convention as q17a);
  * AERONET L1.5 daily AOD at 1064 nm (mean of the evening and morning days, as q5/q15f) where matched;
  * the molecular-model term b_mol = ln slope(real CAMS molecular signal vs US-1976) in the night's own
    window (q15d F3, same helpers), from the CAMS 0.4 deg L137 T/p extract of q15d.
The tier label (q19 record run) is joined at analysis time (q20_corrections.py), because that run may still
be in progress when this table is built.

Run:    .venv/Scripts/python.exe validation/chm15k_stability/q20_build.py
Output: C:/DATA/Projects/202610_CHM15k_stability/q20/q20_nights.csv
"""
from __future__ import annotations

import json
import sys
import warnings
from datetime import timedelta
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from q20_lib import ALL_UNITS, EPOCH, Q20, T_END, T_START  # noqa: E402
from common import BLOCKS_PER_DAY, CENSUS, OUT, load_cache, load_constants  # noqa: E402
import q15d_atmosphere as Q15D  # noqa: E402  (beta_mol, two_way, fit_slope, StdAtm: the q15d F3 code path)
from q5_aeronet import read_aeronet_daily  # noqa: E402

warnings.filterwarnings("ignore")
CAMS_SPECIES = OUT / "q17" / "q17a_cams_points.npz"
CAMS_TP_DIR = OUT / "q15" / "q15d_cams"
AERONET_DIR = Q15D.AERONET_DIR
HK_KEEP = {"T_int": "temp_int", "T_ext": "temp_ext", "T_lom": "temperature_optical_module",
           "pcal": "calibration_pulse", "wt": "window_transmission", "life_h": "laser_life_time",
           "pulses": "laser_pulses"}


def census():
    return {f"{s['wmo']}_{s['ident']}": s for s in json.load(open(CENSUS, encoding="utf-8"))}


def window_hk(cache, t0, t1):
    """Median of every HK variable and the sum of clear profiles over each night's window."""
    tb = (cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY).ravel()
    nv = cache["hk"].shape[-1]
    hk = cache["hk"].reshape(-1, nv)
    nc = cache["n_clear"].ravel().astype(float)
    order = np.argsort(tb)
    tb, hk, nc = tb[order], hk[order], nc[order]
    med = np.full((len(t0), nv), np.nan)
    sc = np.full(len(t0), np.nan)
    i0 = np.searchsorted(tb, t0)
    i1 = np.searchsorted(tb, t1, side="right")
    for k in range(len(t0)):
        if i1[k] > i0[k]:
            med[k] = np.nanmedian(hk[i0[k]:i1[k]], axis=0)
            sc[k] = np.nansum(nc[i0[k]:i1[k]])
    return med, sc, float(cache["day"].max()) + 1.0


def load_cams_tp(keys):
    """T, p, alt profiles [station, time, level] for the requested stream keys, all CAMS files."""
    acc = {v: [] for v in ("time", "T", "alt", "p")}
    idx = None
    for f in sorted(CAMS_TP_DIR.glob("cams_*.npz")):
        with np.load(f, allow_pickle=True) as z:
            if idx is None:
                kk = list(z["keys"])
                idx = [kk.index(k) if k in kk else -1 for k in keys]
            sel = np.array([max(i, 0) for i in idx])
            acc["time"].append(z["time"])
            for v in ("T", "alt", "p"):
                acc[v].append(z[v][sel])
    time = np.concatenate(acc["time"])
    order = np.argsort(time, kind="stable")
    time = time[order]
    keep = np.r_[True, np.diff(time) > 1e-6]               # daily files overlap: keep the first occurrence
    sel = order[keep]
    out = {"time": time[keep]}
    for v in ("T", "alt", "p"):
        out[v] = np.concatenate(acc[v], axis=1)[:, sel, :]
    for j, i in enumerate(idx):
        if i < 0:
            for v in ("T", "alt", "p"):
                out[v][j] = np.nan
    return out


def b_mol_night(cams, j, t0, t1, alt_stn, zb, zt):
    """q15d F3: ln slope of the CAMS molecular signal against the US-1976 model in the window [zb, zt] ASL."""
    m = (cams["time"] >= t0 - 0.07) & (cams["time"] <= t1 + 0.07)
    if not m.any() or not (np.isfinite(zb) and np.isfinite(zt)):
        return np.nan
    alt = np.nanmean(cams["alt"][j, m, :], axis=0)[::-1].astype(float)
    T = np.nanmean(cams["T"][j, m, :], axis=0)[::-1].astype(float)
    p = np.nanmean(cams["p"][j, m, :], axis=0)[::-1].astype(float)
    if not (np.all(np.isfinite(alt)) and np.all(np.isfinite(T)) and np.all(np.isfinite(p))):
        return np.nan
    dz = Q15D.DZ
    z = np.arange(alt_stn + dz, max(zt, 6000) + dz, dz)
    r = z - alt_stn
    bm_c = Q15D.beta_mol(np.interp(z, alt, T), np.exp(np.interp(z, alt, np.log(p))))
    bm_s = np.interp(z, Q15D.STD.z, Q15D.STD.bm)
    w = (z >= zb) & (z <= zt)
    if w.sum() < 5:
        return np.nan
    tm_s = Q15D.two_way(Q15D.S_MOL * bm_s, dz)
    tm_c = Q15D.two_way(Q15D.S_MOL * bm_c, dz)
    x = bm_s * tm_s / r ** 2
    b = Q15D.fit_slope(x[w], (bm_c * tm_c / r ** 2)[w]) - 1.0
    return float(np.log1p(b)) if np.isfinite(b) else np.nan


def aeronet_tables():
    match = json.load(open(AERONET_DIR.parent / "aeronet_match.json", encoding="utf-8"))
    out = {}
    for key, _site, ae_site, dist in match:
        f = AERONET_DIR / f"AOD15_daily_{ae_site}.txt"
        if f.exists():
            out[key] = (ae_site, float(dist), read_aeronet_daily(f))
    return out


def main():
    Q20.mkdir(parents=True, exist_ok=True)
    cen = census()
    keys = [k for k, *_ in ALL_UNITS]
    print("loading CAMS T/p ...", flush=True)
    cams = load_cams_tp(keys)
    Q15D.STD = Q15D.StdAtm()
    zs = np.load(CAMS_SPECIES, allow_pickle=True)
    sp_keys = list(zs["keys"])
    sp_names = list(zs["sl_names"])
    sp_t = np.round(zs["sl_t"]).astype(int)
    sp_pos = {t: i for i, t in enumerate(sp_t)}
    aer = aeronet_tables()
    rows = []
    for j, (key, label, country, setting) in enumerate(ALL_UNITS):
        wmo, ident = key.rsplit("_", 1)
        s = cen[key]
        alt_stn = float(s["alt"])
        cons = load_constants(wmo, ident)
        sel = (cons["t0"] >= T_START) & (cons["t0"] < T_END) & np.isfinite(cons["c"]) & (cons["c"] > 0)
        cons = {k: v[sel] for k, v in cons.items()}
        n = len(cons["t"])
        cache = load_cache(wmo, ident)
        med, nclear, cache_end = window_hk(cache, cons["t0"], cons["t1"])
        names = list(cache["hk_vars"])
        night = np.floor(cons["t0"] - 0.5).astype(int)
        d = dict(stream=[key] * n, label=[label] * n, country=[country] * n, setting=[setting] * n,
                 t=cons["t"], t0=cons["t0"], t1=cons["t1"], night=night,
                 date=[(EPOCH + timedelta(days=float(np.floor(x + 0.5)))).strftime("%Y%m%d") for x in cons["t0"]],
                 c=cons["c"], u=cons["u"], zb_asl=cons["zb"], zt_asl=cons["zt"], module=cons["module"],
                 alt=np.full(n, alt_stn), lat=np.full(n, s["lat"]), lon=np.full(n, s["lon"]),
                 in_cache=(cons["t1"] <= cache_end))
        for col, var in HK_KEEP.items():
            d[col] = med[:, names.index(var)] if var in names else np.full(n, np.nan)
        d["n_clear"] = nclear
        # CAMS AOD550 by species at 00 UTC of the night (q17a convention: night key + 1)
        jj = sp_keys.index(key) if key in sp_keys else -1
        ii = np.array([sp_pos.get(int(x) + 1, -1) for x in night])
        for vi, nm in enumerate(sp_names):
            v = np.full(n, np.nan)
            if jj >= 0:
                ok = ii >= 0
                v[ok] = zs["sl_data"][ii[ok], jj, vi]
            d[nm] = v
        # AERONET: mean of the days of the window start and end
        a1064 = np.full(n, np.nan)
        if key in aer:
            daily = aer[key][2]
            for k in range(n):
                days = {(EPOCH + timedelta(days=float(cons["t0"][k]))).date(),
                        (EPOCH + timedelta(days=float(cons["t1"][k]))).date()}
                vals = [daily[x][0] for x in days if x in daily]
                if vals:
                    a1064[k] = float(np.mean(vals))
        d["aod1064_aeronet"] = a1064
        d["aeronet_site"] = [aer[key][0] if key in aer else ""] * n
        # molecular-model term (q15d F3)
        d["b_mol"] = np.array([b_mol_night(cams, j, cons["t0"][k], cons["t1"][k], alt_stn,
                                           cons["zb"][k], cons["zt"][k]) for k in range(n)])
        df = pd.DataFrame(d)
        rows.append(df)
        print(f"{key:32s} {label:17s} nights {n:4d}  HK {np.isfinite(df.T_int).sum():4d}  CAMS sp "
              f"{np.isfinite(df.duaod550).sum():4d}  AERONET {np.isfinite(df.aod1064_aeronet).sum():4d}  "
              f"b_mol {np.isfinite(df.b_mol).sum():4d} (median {np.nanmedian(df.b_mol) * 100:+.2f} %)  "
              f"modules {sorted(set(cons['module']))}", flush=True)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(Q20 / "q20_nights.csv", index=False)
    print(f"-> {Q20 / 'q20_nights.csv'} ({len(out)} nights)")


if __name__ == "__main__":
    main()
