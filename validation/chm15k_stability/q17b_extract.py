# -*- coding: utf-8 -*-
"""q17b (extraction): data for the wavelength-scaling test CHM15k (1064 nm) vs CL61 (910.55 nm).

Question (q17b): is the residual solstice-phased annual cycle of the CHM15k Rayleigh constant an
instrument gain, or aerosol inside the 'molecular' reference? At co-located sites the CL61 sees the
same air; an aerosol layer gives an excess R-1 at 910 nm of ~0.55-0.68x the 1064 nm one, a CHM15k
gain gives ~0 at the CL61. This script only EXTRACTS; the test itself is in q17b_analyse.py.

Stages
  cams    CAMS 0.4 deg L137 columns (t, q, model-level altitude, lnsp -> full-level pressure) at the
          4 sites, levels 77..137 (surface .. ~11 km), every time step of D:/CAMS_Monthly_04
          (2025-01 .. 2026-06) and of the main-domain daily files D:/CAMS_daily (2026-06 ..).
          Never the 1 deg folders (wrong orography).
  blocks  native L1 (D:/E-PROFILE_L1_2026; Lindenberg CL61 fallback D:/E-PROFILE_L1_CL61_LIND) of the
          CHM15k and the CL61 of each site -> per UTC day, 5-min blocks x 100 m vertical bins
          (0-12 km above the instrument, z = range * cos(tilt)): block-mean rcs_0, profile counts,
          a cloud flag (any cloud_base_height > 0 in the block; CL61 also vertical visibility) and
          block-mean housekeeping (temp_int; CL61 laser temperature / energy / window transmission).

Outputs: C:/DATA/Projects/202610_CHM15k_stability/q17/q17b_cams_<site>.npz,
         q17/q17b_blocks/<site>_<instr>.npz
Run:     .venv/Scripts/python.exe validation/chm15k_stability/q17b_extract.py cams|blocks [--workers 5]
"""
from __future__ import annotations

import argparse
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import date, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from common import OUT  # noqa: E402
from calibration.water_vapor_correction.water_vapor import _A137, _B137  # noqa: E402

warnings.filterwarnings("ignore")
Q17 = OUT / "q17"
BLK = Q17 / "q17b_blocks"
L1_MAIN = Path("D:/E-PROFILE_L1_2026")
L1_LIND = Path("D:/E-PROFILE_L1_CL61_LIND")

# site -> (wmo, CHM15k ident, CL61 ident, lat, lon)
SITES = {
    "Aosta": ("0-380-5-1", "0", "B", 45.7422, 7.357),
    "Lindenberg": ("0-20000-0-10393", "0", "C", 52.21, 14.12),
    "Camborne": ("0-20000-0-03808", "A", "C", 50.2183, -5.3275),
    "Payerne": ("0-20000-0-06610", "A", "C", 46.81369, 6.942547),
}
FIRST_DAY = date(2025, 1, 1)
LAST_DAY = date(2026, 10, 15)
NBLK = 288                                   # 5-min blocks per UTC day
Z_EDGES = np.arange(0.0, 12000.1, 100.0)      # vertical bins above the instrument [m]
NZ = len(Z_EDGES) - 1
LEV_MIN = 77                                  # model level 77 ~ 11 km; 137 = surface


# ---------------------------------------------------------------------------------------------
# CAMS
# ---------------------------------------------------------------------------------------------
def _cams_files():
    monthly = sorted(Path("D:/CAMS_Monthly_04").glob("CAMS_Beta_2*.nc"))
    daily = sorted(f for f in Path("D:/CAMS_daily").glob("CAMS_Beta_2026*.nc") if len(f.stem) == 18)
    return monthly, daily


def _cams_one(path):
    """All time steps of one CAMS file at the 4 sites: t, q, altitude, p on levels >= LEV_MIN."""
    import netCDF4
    with netCDF4.Dataset(path) as ds:
        lat = np.asarray(ds["latitude"][:], float)
        lon = np.asarray(ds["longitude"][:], float)
        level = np.asarray(ds["level"][:], float).astype(int)
        keep = np.where(level >= LEV_MIN)[0]
        lv = level[keep]
        tv = ds["time"]
        tt = netCDF4.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)
        tdays = np.array([(np.datetime64(x) - np.datetime64("1970-01-01")) / np.timedelta64(1, "D") for x in tt])
        out = {}
        for site, (_, _, _, la, lo) in SITES.items():
            i = int(np.argmin(np.abs(lat - la)))
            j = int(np.argmin(np.abs(lon - lo)))
            T = np.ma.filled(ds["t"][:, keep[0]:keep[-1] + 1, i, j], np.nan).astype("f4")
            q = np.ma.filled(ds["q"][:, keep[0]:keep[-1] + 1, i, j], np.nan).astype("f4")
            alt = np.ma.filled(ds["altitude"][:, keep[0]:keep[-1] + 1, i, j], np.nan).astype("f4")
            lnsp_all = np.ma.filled(ds["lnsp"][:, :, i, j], np.nan)
            slot = int(np.where(np.all(np.isfinite(lnsp_all), axis=0))[0][0])
            ps = np.exp(lnsp_all[:, slot])
            ph_lo = _A137[lv][None, :] + _B137[lv][None, :] * ps[:, None]
            ph_hi = _A137[lv - 1][None, :] + _B137[lv - 1][None, :] * ps[:, None]
            p = (0.5 * (ph_lo + ph_hi)).astype("f4")
            out[site] = dict(time=tdays, T=T, q=q, alt=alt, p=p, glat=lat[i], glon=lon[j])
    return str(path), out


def stage_cams(workers):
    monthly, daily = _cams_files()
    files = [str(f) for f in monthly + daily]
    acc = {s: [] for s in SITES}
    with ProcessPoolExecutor(workers) as ex:
        for path, res in ex.map(_cams_one, files):
            daily_file = "CAMS_daily" in path
            for s, d in res.items():
                d["daily"] = daily_file
                acc[s].append(d)
            print("cams", Path(path).name, flush=True)
    for s, lst in acc.items():
        # monthly first; a daily step is used only for times the monthly files do not cover
        t_m = np.concatenate([d["time"] for d in lst if not d["daily"]])
        tmax_m = t_m.max()
        parts = {k: [] for k in ("time", "T", "q", "alt", "p")}
        seen = set()
        for d in sorted(lst, key=lambda x: (x["daily"], x["time"][0])):
            for k in range(len(d["time"])):
                tk = round(float(d["time"][k]) * 8) / 8
                if tk in seen or (d["daily"] and tk <= tmax_m):
                    continue
                seen.add(tk)
                for v in parts:
                    parts[v].append(d[v][k] if v != "time" else tk)
        o = np.argsort(parts["time"])
        arr = {v: np.asarray(parts[v])[o] for v in parts}
        np.savez_compressed(Q17 / f"q17b_cams_{s}.npz", glat=lst[0]["glat"], glon=lst[0]["glon"], **arr)
        print(s, "CAMS steps", len(arr["time"]), flush=True)


# ---------------------------------------------------------------------------------------------
# L1 blocks
# ---------------------------------------------------------------------------------------------
def l1_file(wmo, ident, d):
    for root in (L1_MAIN, L1_LIND):
        f = root / wmo / f"{d:%Y}" / f"{d:%m}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"
        if f.exists():
            return f
    return None


HK = {"CHM15k": ("temp_int",), "CL61": ("temp_int", "temperature_laser", "laser_energy", "window_transmission")}


def _read_day(f, instr):
    import netCDF4
    with netCDF4.Dataset(f) as ds:
        t = np.ma.filled(ds["time"][:].astype("f8"), np.nan)
        rng = np.ma.filled(ds["range"][:].astype("f8"), np.nan)
        if "tilt_angle" in ds.variables:
            tl = np.ma.filled(ds["tilt_angle"][...].astype("f8"), np.nan)
            tilt = float(np.nanmedian(tl)) if np.size(tl) else 0.0
        else:
            tilt = 0.0
        tilt = 0.0 if not np.isfinite(tilt) else tilt
        z = rng * np.cos(np.deg2rad(tilt))
        ng = int(np.searchsorted(z, Z_EDGES[-1]))
        rcs = np.ma.filled(ds["rcs_0"][:, :ng].astype("f4"), np.nan)
        cbh = np.ma.filled(ds["cloud_base_height"][:].astype("f4"), np.nan)
        cloud = np.any(np.isfinite(cbh) & (cbh > 0), axis=1)
        if instr == "CL61" and "vertical_visibility" in ds.variables:
            vv = np.ma.filled(ds["vertical_visibility"][:].astype("f4"), np.nan)
            cloud |= np.isfinite(vv) & (vv > 0)
        hk = {}
        for v in HK[instr]:
            hk[v] = np.ma.filled(ds[v][:].astype("f8"), np.nan) if v in ds.variables else np.full(len(t), np.nan)
        alt = float(np.ma.filled(ds["station_altitude"][...], np.nan))
    return t, z[:ng], rcs, cloud, hk, tilt, alt


def _blocks_task(args):
    site, instr, ident, wmo, days = args
    out = []
    for d in days:
        f = l1_file(wmo, ident, d)
        if f is None:
            continue
        try:
            t, z, rcs, cloud, hk, tilt, alt = _read_day(f, instr)
        except Exception as e:  # noqa: BLE001 - a corrupt file must not kill the run
            print("  read fail", f.name, e, flush=True)
            continue
        day0 = (d - date(1970, 1, 1)).days
        ib = np.floor((t - day0) * NBLK).astype(int)
        ok = (ib >= 0) & (ib < NBLK) & np.isfinite(t)
        zb = np.digitize(z, Z_EDGES) - 1
        # gate -> bin sums (rcs NaN -> excluded)
        fin = np.isfinite(rcs)
        r0 = np.where(fin, rcs, 0.0).astype("f8")
        S = np.zeros((NBLK, NZ))
        N = np.zeros((NBLK, NZ))
        G = np.zeros((NBLK, NZ))
        # sum over gates in each bin first (per profile), then over profiles in each block
        gb_s = np.zeros((rcs.shape[0], NZ))
        gb_n = np.zeros((rcs.shape[0], NZ))
        for k in range(NZ):
            g = zb == k
            if g.any():
                gb_s[:, k] = r0[:, g].sum(axis=1)
                gb_n[:, k] = fin[:, g].sum(axis=1)
        np.add.at(S, ib[ok], gb_s[ok])
        np.add.at(N, ib[ok], gb_n[ok])
        prof_n = np.bincount(ib[ok], minlength=NBLK).astype("f4")
        cl = np.bincount(ib[ok], weights=cloud[ok].astype(float), minlength=NBLK)
        hkb = {}
        for v, x in hk.items():
            xx = np.where(np.isfinite(x[ok]), x[ok], 0.0)
            nn = np.bincount(ib[ok], weights=np.isfinite(x[ok]).astype(float), minlength=NBLK)
            hkb[v] = (np.bincount(ib[ok], weights=xx, minlength=NBLK) / np.maximum(nn, 1)).astype("f4")
            hkb[v][nn == 0] = np.nan
        G[:] = 0
        mean = (S / np.where(N > 0, N, np.nan)).astype("f4")
        out.append(dict(day=day0, rcs=mean, nprof=prof_n, ncloud=cl.astype("f4"), tilt=tilt, alt=alt,
                        hk=hkb))
    return site, instr, out


def stage_blocks(workers, only=None):
    BLK.mkdir(parents=True, exist_ok=True)
    tasks = []
    for site, (wmo, ich, icl, _, _) in SITES.items():
        if only and site not in only:
            continue
        for instr, ident in (("CHM15k", ich), ("CL61", icl)):
            d = FIRST_DAY
            month_days = {}
            while d <= LAST_DAY:
                month_days.setdefault((d.year, d.month), []).append(d)
                d += timedelta(days=1)
            for days in month_days.values():
                tasks.append((site, instr, ident, wmo, days))
    acc = {}
    with ProcessPoolExecutor(workers) as ex:
        for site, instr, out in ex.map(_blocks_task, tasks):
            acc.setdefault((site, instr), []).extend(out)
            print(f"blocks {site} {instr}: +{len(out)} days", flush=True)
    for (site, instr), lst in acc.items():
        if not lst:
            continue
        lst.sort(key=lambda x: x["day"])
        hkv = HK[instr]
        np.savez_compressed(
            BLK / f"{site}_{instr}.npz", day=np.array([x["day"] for x in lst]),
            rcs=np.stack([x["rcs"] for x in lst]), nprof=np.stack([x["nprof"] for x in lst]),
            ncloud=np.stack([x["ncloud"] for x in lst]), tilt=np.array([x["tilt"] for x in lst]),
            alt=np.array([x["alt"] for x in lst]), z_edges=Z_EDGES, hk_vars=np.array(hkv),
            hk=np.stack([np.stack([x["hk"][v] for v in hkv], axis=-1) for x in lst]))
        print("saved", site, instr, len(lst), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=("cams", "blocks"))
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--sites", nargs="*")
    a = ap.parse_args()
    Q17.mkdir(parents=True, exist_ok=True)
    if a.stage == "cams":
        stage_cams(a.workers)
    else:
        stage_blocks(a.workers, a.sites)
