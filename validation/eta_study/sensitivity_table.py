"""SENSITIVITIES (not fits) of the cloud C_L to the inputs of the proposed eta, per Vaisala type.

Around the physics-based central values (PVC, a_G 3.3 um from Cloudnet, alpha = network type median
alpha_k from the ceilometer tail slope, type optics from validation/multiple_scattering_eta.OPTICS,
operational WV line, S = 18.8 sr) each input is perturbed ONE at a time and the change of
  - the de-seasoned CBH slope of ln C_L (per profile, unit x month fixed effects, day-clustered CI), and
  - the C_L level at 1 km (profiles with CBH 0.9-1.1 km)
is reported. This documents the uncertainty budget; it does NOT choose inputs (a flat C_L is a check,
not a fitting target). Perturbations:
  a_G +-0.5 um; alpha +-5 /km; rho_t x sqrt(2) and / sqrt(2) (half-angle vs full-angle / 1/e conventions);
  rho_l x sqrt(2) and / sqrt(2); CL31 WV line centre +-1 nm (FWHM 6) and FWHM 3 / 9 nm (centre 909.7):
  the stored T2 belongs to the operational line and humidity is not stored, so a line change enters as
  T2 -> T2**R with R = ln T2_new / ln T2_oper at 1 km from the HITRAN LUT on a standard humidity profile
  (PWV 1.5 cm, scale height 2 km) -- an approximation; S +-0.5 sr changes the level only (C_L ~ S_theo).
Per-profile C_L of a variant from the stored windows as in cbh_slope_windowed.py.

Run:  python validation/eta_study/sensitivity_table.py --stats --src DIR --out DIR --tag NAME
      python validation/eta_study/sensitivity_table.py --pool --stats-files F1;F2 --out DIR
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

A_G = 3.3
ALPHA = {"CL31": 27.8, "CL51": 30.2, "CL61": 31.0}           # network type medians (unit medians)
S_THEO = 18.8
CL31_LINE = (909.7, 6.0)                                      # operational LASER_SPECTRUM["CL31"]
LEVEL_BAND = (0.9, 1.1)                                       # km, "level at 1 km"


def wv_strength(l0, fwhm, ref=CL31_LINE, z_km=1.0):
    """R = ln T2(line) / ln T2(ref) at z_km on a standard exponential humidity profile (HITRAN LUT)."""
    from calibration.water_vapor_correction.water_vapor import load_abs_cross_section, DEFAULT_ABS_CROSS_SECTION
    wl, hgt, cs = load_abs_cross_section(DEFAULT_ABS_CROSS_SECTION)
    z = np.arange(0.0, z_km * 1000 + 1, 10.0)
    n = 15.0 / (2000.0 * 2.99e-26) * np.exp(-z / 2000.0) * 1e-6                # molecules cm^-3
    tau = (np.array([np.interp(z, hgt, cs[i]) for i in range(wl.size)]) * n[None, :] * 1000.0).sum(1)

    def lnT2(lc, fw):
        g = np.exp(-0.5 * ((wl - lc) / (fw / 2.3548)) ** 2)
        return np.log((g * np.exp(-2 * tau)).sum() / g.sum())
    return float(lnT2(l0, fwhm) / lnT2(*ref))


def variants(itype):
    from validation.multiple_scattering_eta import OPTICS
    rt0, rl0 = OPTICS[itype]
    c = dict(rho_t=rt0, rho_l=rl0, alpha=ALPHA[itype], a_g=A_G, wv_r=1.0)
    v = [dict(c, name="operational"), dict(c, name="central")]
    v += [dict(c, name="a_G -0.5 um", a_g=A_G - 0.5), dict(c, name="a_G +0.5 um", a_g=A_G + 0.5),
          dict(c, name="alpha -5 /km", alpha=ALPHA[itype] - 5), dict(c, name="alpha +5 /km", alpha=ALPHA[itype] + 5),
          dict(c, name="rho_t / sqrt2", rho_t=rt0 / np.sqrt(2)), dict(c, name="rho_t x sqrt2", rho_t=rt0 * np.sqrt(2)),
          dict(c, name="rho_l / sqrt2", rho_l=rl0 / np.sqrt(2)), dict(c, name="rho_l x sqrt2", rho_l=rl0 * np.sqrt(2))]
    if itype == "CL31":
        for name, (l0, fw) in (("WV centre -1 nm", (908.7, 6.0)), ("WV centre +1 nm", (910.7, 6.0)),
                               ("WV FWHM 3 nm", (909.7, 3.0)), ("WV FWHM 9 nm", (909.7, 9.0))):
            v.append(dict(c, name=name, wv_r=wv_strength(l0, fw)))
    return v


def table(job):
    itype, v = job
    import calibration.cloud._filters as CF
    from validation.multiple_scattering_eta import WAVELENGTH, eta_for_cloud_base
    oper = {"CL31": CF._ETA_CL31, "CL51": CF._ETA_CL51, "CL61": CF._ETA_CL61}[itype]
    if v["name"] == "operational":
        return oper
    th = WAVELENGTH[itype] / (np.pi * v["a_g"] * 1e-6)
    return np.column_stack([oper[:, 0], [eta_for_cloud_base(z * 1000, v["rho_t"], v["rho_l"], th, alpha=v["alpha"] * 1e-3)
                                         for z in oper[:, 0]]])


def unit_stats(job):
    """Per-day sufficient statistics of y_v = ln C_L,v (unit x month demeaned) and level sums at 1 km."""
    key, src, tabs, wvr = job
    from calibration.cloud._filters import _interp1_linear_extrap
    z = np.load(Path(src) / f"{key}.npz", allow_pickle=True)
    B = z["beta_before"].astype("float64")
    T2 = z["t2_wv"].astype("float64")
    r_abs = (z["peak_range_m"].astype("float64")[:, None] + z["range_rel_m"].astype("float64")[None, :]) / 1000.0
    good = np.isfinite(B) & np.isfinite(T2) & (T2 > 0)
    Bz, T2s = np.where(good, B, 0.0), np.where(good, T2, 1.0)
    c = z["c_profile"].astype("float64")
    cbh = z["cbh_m"].astype("float64") / 1000.0
    I = np.array([(Bz / T2s ** r * _interp1_linear_extrap(t[:, 0], t[:, 1], r_abs)).sum(1) for t, r in zip(tabs, wvr)])
    ok = np.isfinite(c) & (c > 0) & np.isfinite(cbh) & np.all(np.isfinite(I) & (I > 0), axis=0)
    y = -np.log(c[ok])[None, :] + np.log(I[:, ok] / I[0, ok])          # variant 0 = operational
    x = cbh[ok]
    date = z["date"].astype(int)[ok]
    month = date // 100
    lev = (x >= LEVEL_BAND[0]) & (x <= LEVEL_BAND[1])
    dlev = (y[:, lev] - y[1:2, lev]).sum(1)                           # sum of ln(C_L,v / C_L,central) near 1 km
    xt, yt = x.copy(), y.copy()
    for m in np.unique(month):
        s = month == m
        xt[s] -= x[s].mean()
        yt[:, s] -= y[:, s].mean(axis=1, keepdims=True)
    days = np.unique(date)
    idx = np.searchsorted(days, date)
    sxx = np.bincount(idx, weights=xt * xt, minlength=days.size)
    sxy = np.array([np.bincount(idx, weights=xt * yt[k], minlength=days.size) for k in range(y.shape[0])])
    return key, sxx, sxy, dlev, int(lev.sum())


def run_stats(a):
    src = Path(a.src)
    keys = sorted(p.stem for p in src.glob("*.npz"))
    itype = {k: json.loads(str(np.load(src / f"{k}.npz", allow_pickle=True)["meta"]))["itype"] for k in keys}
    V = {it: variants(it) for it in ALPHA}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        T = {it: list(ex.map(table, [(it, v) for v in V[it]])) for it in ALPHA}
        jobs = [(k, str(src), T[itype[k]], [v["wv_r"] for v in V[itype[k]]]) for k in keys if itype[k] in V]
        res = list(ex.map(unit_stats, jobs, chunksize=1))
    arr = {}
    for key, sxx, sxy, dlev, nlev in res:
        arr[f"{key}|itype"] = np.array(itype[key])
        arr[f"{key}|sxx"], arr[f"{key}|sxy"] = sxx, sxy
        arr[f"{key}|dlev"], arr[f"{key}|nlev"] = dlev, np.array(nlev)
    arr["variants"] = np.array(json.dumps({it: [{k: (float(x) if isinstance(x, (float, np.floating)) else x)
                                                 for k, x in d.items()} for d in V[it]] for it in V}))
    np.savez_compressed(out / f"sensitivity_stats_{a.tag}.npz", **arr)
    print(f"{len(res)} units -> {out}", flush=True)


def run_pool(a):
    import pandas as pd
    data, V = {}, None
    for f in a.stats_files.split(";"):
        z = np.load(f, allow_pickle=True)
        V = json.loads(str(z["variants"]))
        for k in {n.split("|")[0] for n in z.files if "|" in n}:
            data[k] = {n: z[f"{k}|{n}"] for n in ("itype", "sxx", "sxy", "dlev", "nlev")}
    rows = []
    for it in ("CL31", "CL51", "CL61"):
        ks = [k for k, d in data.items() if str(d["itype"]) == it]
        sxx = np.concatenate([data[k]["sxx"] for k in ks])
        sxy = np.concatenate([data[k]["sxy"] for k in ks], axis=1)
        dlev = np.sum([data[k]["dlev"] for k in ks], axis=0)
        nlev = int(np.sum([data[k]["nlev"] for k in ks]))
        b = sxy.sum(1) / sxx.sum()
        db = b - b[1]                                                    # change vs central
        g = (sxy - sxy[1:2]) - db[:, None] * sxx[None, :]
        G = sxx.size
        se = np.sqrt(G / (G - 1) * (g ** 2).sum(1)) / sxx.sum()
        for v, bb, d, s, dl in zip(V[it], b, db, se, dlev):
            rows.append(dict(itype=it, n_units=len(ks), input=v["name"], wv_R=v["wv_r"],
                             slope_pct_km=100 * bb, d_slope_pct_km=100 * d, d_slope_ci95=196 * s,
                             d_level_1km_pct=100 * (np.exp(dl / max(nlev, 1)) - 1), n_profiles_1km=nlev))
        for sgn in (-1, 1):
            rows.append(dict(itype=it, n_units=len(ks), input=f"S {'-' if sgn < 0 else '+'}0.5 sr", wv_R=np.nan,
                             slope_pct_km=np.nan, d_slope_pct_km=0.0, d_slope_ci95=0.0,
                             d_level_1km_pct=100 * sgn * 0.5 / S_THEO, n_profiles_1km=nlev))
    t = pd.DataFrame(rows)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / "sensitivity_table.csv", index=False, float_format="%.3f")
    print(t.round(2).to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", action="store_true")
    ap.add_argument("--pool", action="store_true")
    ap.add_argument("--src", default="C:/DATA/Projects/202611_eta_study/incloud")
    ap.add_argument("--out", default="C:/DATA/Projects/202611_eta_study/sensitivity")
    ap.add_argument("--tag", default="pc")
    ap.add_argument("--stats-files", default="")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    if a.stats:
        run_stats(a)
    if a.pool:
        run_pool(a)


if __name__ == "__main__":
    main()
