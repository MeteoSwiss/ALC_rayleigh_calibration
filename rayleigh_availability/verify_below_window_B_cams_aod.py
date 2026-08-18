# -*- coding: utf-8 -*-
"""ADVERSARIAL CHECK B -- is the aerosol the mechanism needs actually present on the real nights?

Fast netCDF4 rewrite of verif_B (xarray's .sel/.transpose pulled whole 344 MB fields; direct
index slicing reads one column in 0.6 s).

The claim needs AOD_1064 ~ 0.38 BELOW the fit window, with a TRUE aerosol lidar ratio S = 90 sr,
to explain the -22.5 % Payerne offset -- and needs S < 52 sr at the sites whose offset is positive.
It judged that implausible against a GENERIC AOD_1064 = 0.05. The right benchmark is what CAMS
actually carries over each station on the very nights the study calls v2.0-KEPT and
v2.2-RECOVERED, since the mechanism only has to explain the difference between those populations.

Measured here, per night, from D:/CAMS_run_v20 (L137 model levels):
  * AOD_1064 = INT aerext1064 dz   [-]   over the whole column, below 1500 m AGL, and 2-6 km AGL;
  * the local aerosol lidar ratio S_1064 = aerext1064 / aerbackscatgnd1064  [sr], median over the
    first 1 km AGL where the two-way attenuation of the "from the ground" backscatter is still
    within a few per cent of 1 (higher up ext/bsc overestimates S by exp(+2*tau)).

Altitudes: L137 hydrostatic integration (calibration.water_vapor_correction._hydrostatic_z_p) --
in these files ``z``/``lnsp`` are SURFACE fields, so reading ``z`` as a per-level geopotential
returns a constant (-3341 m at Payerne) and silently zeroes every integral.
Night sampling: the 00 and 03 UTC CAMS steps of the calibration date.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

REPO = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
sys.path.insert(0, str(REPO))
from calibration.water_vapor_correction.water_vapor import _hydrostatic_z_p  # noqa: E402

CAMS_DIR = Path("D:/CAMS_run_v20")
RES = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
OUT = Path(__file__).parent / "verif_B.json"

STATIONS = {
    "PAYERNE_CHM15k_A":    dict(lat=46.81369, lon=6.942547, alt=490.0, offset_pct=-22.5,
                                base="base_eprof_v2_PAYERNE_CHM15k_A.json",
                                cand="cand_N2.5_PAYERNE_CHM15k_A.json"),
    "LINDENBERG_CHM15k_0": dict(lat=52.21, lon=14.12, alt=123.0, offset_pct=+23.2,
                                base="base_eprof_v2_LINDENBERG_CHM15k_0.json",
                                cand="cand_N2.5_LINDENBERG_CHM15k_0.json"),
    "PALAISEAU_CHM15k_B":  dict(lat=48.71810, lon=2.20740, alt=156.0, offset_pct=+7.4,
                                base="base_eprof_v2_PALAISEAU_CHM15k_B.json",
                                cand="cand_N2.5_PALAISEAU_CHM15k_B.json"),
}

_files: dict = {}
_cols: dict = {}


def _columns(date: str, lat: float, lon: float):
    """(times 'YYYYMMDDHH', alt_asl [m] (n_lev,n_t), ext [m^-1], bsc [m^-1 sr^-1]) or None."""
    path = None
    for name in (f"CAMS_Beta_{date[:6]}.nc", f"CAMS_Beta_{date}.nc"):
        if (CAMS_DIR / name).exists():
            path = CAMS_DIR / name
            break
    if path is None:
        return None
    key = (str(path), lat, lon)
    if key in _cols:
        return _cols[key]

    if str(path) not in _files:
        _files[str(path)] = nc.Dataset(path)
    ds = _files[str(path)]
    if "aerext1064" not in ds.variables:
        _cols[key] = None
        return None

    lats = np.asarray(ds.variables["latitude"][:], float)
    lons = np.asarray(ds.variables["longitude"][:], float)
    j = int(np.argmin(np.abs(lats - lat)))
    i = int(np.argmin(np.abs(lons - lon)))
    tt = nc.num2date(ds.variables["time"][:], ds.variables["time"].units,
                     only_use_cftime_datetimes=False)
    stamps = np.array([t.strftime("%Y%m%d%H") for t in tt])

    level = np.asarray(ds.variables["level"][:], int)
    order = np.argsort(level)                                   # ascending: 1 = top
    ext = np.asarray(ds.variables["aerext1064"][:, :, j, i], float).T[order]      # (n_lev,n_t)
    if "aerbackscatgnd1064" in ds.variables:
        bsc = np.asarray(ds.variables["aerbackscatgnd1064"][:, :, j, i], float).T[order]
    else:
        bsc = np.full_like(ext, np.nan)
    T = np.asarray(ds.variables["t"][:, :, j, i], float).T[order]
    q = np.asarray(ds.variables["q"][:, :, j, i], float).T[order]
    z_raw = np.asarray(ds.variables["z"][:, :, j, i], float).T[order]
    lnsp_raw = np.asarray(ds.variables["lnsp"][:, :, j, i], float).T[order]

    # z / lnsp are SURFACE fields written into one level slot: take the first finite per time.
    def _surface(a):
        out = np.full(a.shape[1], np.nan)
        for k in range(a.shape[1]):
            f = np.where(np.isfinite(a[:, k]))[0]
            if f.size:
                out[k] = a[f[0], k]
        return out

    z_surf, lnsp = _surface(z_raw), _surface(lnsp_raw)
    good = np.isfinite(z_surf) & np.isfinite(lnsp)
    alt = np.full_like(T, np.nan)
    if good.any():
        gi = np.where(good)[0]
        alt[:, gi], _ = _hydrostatic_z_p(level[order], T[:, gi], q[:, gi],
                                         z_surf[gi], lnsp[gi])
    out = (stamps, alt, ext, bsc)
    _cols[key] = out
    return out


def night(date, lat, lon):
    got = _columns(date, lat, lon)
    if got is None:
        return None
    stamps, alt, ext, bsc = got
    idx = np.where((stamps == date + "00") | (stamps == date + "03"))[0]
    if idx.size == 0:
        return None
    a = np.nanmean(alt[:, idx], axis=1)
    e = np.nanmean(ext[:, idx], axis=1)
    b = np.nanmean(bsc[:, idx], axis=1)
    m = np.isfinite(a) & np.isfinite(e)
    if m.sum() < 10:
        return None
    a, e, b = a[m], e[m], b[m]
    o = np.argsort(a)
    return a[o], e[o], b[o]


def aod(a, e, lo, hi):
    m = (a >= lo) & (a <= hi)
    return float(np.trapezoid(e[m], a[m])) if m.sum() >= 2 else np.nan


def s_local(a, e, b, lo, hi):
    m = (a >= lo) & (a <= hi) & np.isfinite(b) & (b > 0) & (e > 0)
    return float(np.median(e[m] / b[m])) if m.sum() >= 2 else np.nan


def main():
    out = {}
    for label, cfg in STATIONS.items():
        pb, pc = RES / "baselines" / cfg["base"], RES / "candidates" / cfg["cand"]
        base = json.loads(pb.read_text(encoding="utf-8"))
        cand = json.loads(pc.read_text(encoding="utf-8"))
        ok = lambda r: r is not None and r[0] in (1.0, 0.5)   # noqa: E731

        groups = {"kept": [], "recovered": []}
        for date, rec in base.items():
            if ok(rec):
                groups["kept"].append(date)
            elif ok(cand.get(date)):
                groups["recovered"].append(date)

        rows = {}
        for g, dates in groups.items():
            rows[g] = []
            for date in dates:
                p = night(date, cfg["lat"], cfg["lon"])
                if p is None:
                    continue
                a, e, b = p
                rows[g].append(dict(
                    date=date,
                    aod_col=aod(a, e, cfg["alt"], 30000.0),
                    aod_below1500=aod(a, e, cfg["alt"], cfg["alt"] + 1500.0),
                    aod_2_6km=aod(a, e, cfg["alt"] + 2000.0, cfg["alt"] + 6000.0),
                    S_local_0_1km=s_local(a, e, b, cfg["alt"], cfg["alt"] + 1000.0)))

        def st(g, k):
            v = np.array([r[k] for r in rows[g] if np.isfinite(r[k])])
            if v.size == 0:
                return dict(n=0)
            return dict(n=int(v.size), median=float(np.median(v)),
                        p90=float(np.percentile(v, 90)), p99=float(np.percentile(v, 99)),
                        max=float(v.max()))

        summ = {g: {k: st(g, k) for k in
                    ("aod_col", "aod_below1500", "aod_2_6km", "S_local_0_1km")}
                for g in ("kept", "recovered")}
        out[label] = dict(observed_offset_pct=cfg["offset_pct"], summary=summ, rows=rows,
                          n_kept=len(groups["kept"]), n_recovered=len(groups["recovered"]))

        print(f"\n=== {label}   observed recovered-minus-kept C_L offset "
              f"{cfg['offset_pct']:+.1f} %", flush=True)
        for k in ("aod_col", "aod_below1500", "aod_2_6km", "S_local_0_1km"):
            a_, b_ = summ["kept"][k], summ["recovered"][k]
            if a_.get("n") and b_.get("n"):
                unit = "sr" if k.startswith("S") else "  "
                print(f"    {k:14s} {unit}  KEPT      n={a_['n']:3d} med {a_['median']:8.4f} "
                      f"p90 {a_['p90']:8.4f} p99 {a_['p99']:8.4f} max {a_['max']:8.4f}",
                      flush=True)
                print(f"    {'':14s} {unit}  RECOVERED n={b_['n']:3d} med {b_['median']:8.4f} "
                      f"p90 {b_['p90']:8.4f} p99 {b_['p99']:8.4f} max {b_['max']:8.4f}",
                      flush=True)

    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
