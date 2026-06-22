"""
earlinet.py — read EARLINET L2 1064 nm backscatter, convert to attenuated backscatter (assumed lidar
ratio 50 sr, molecular Rayleigh + two-way transmission), and compare to the colocated CHM15k ceilometer
(operational Python Rayleigh + Kalman calibration). Port of read_earlinet_att_backscatter.m + the
paper_val_earlinet matching (CHM profiles within +/-30 min of each EARLINET time, median, interpolated to
the EARLINET grid; statistics over 500-5000 m AGL). Faithful to _PORT_SPEC.md section 3.
"""
from __future__ import annotations
import glob
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties
from validation.paper import intercompare as IC

EARLINET = Path("A:/EARLINET")
RANGE_REF = np.arange(0, 15001, 15.0)            # m AGL
LR_1064 = 50.0                                   # assumed lidar ratio [sr]
# code -> CHM channel (wmo, ident) + overlap-min [m AGL]
SITES = {
    "sir": dict(wmo="0-250-1001-07151", ident="B", overlap=2000.0),
    "lei": dict(wmo="0-20000-0-10471", ident="0", overlap=800.0),
    "cbw": dict(wmo="0-20000-0-06348", ident="A", overlap=1000.0),
    "ino": dict(wmo="0-20008-0-INO", ident="B", overlap=1100.0),
    "ari": dict(wmo="0-20000-0-10471", ident="0", overlap=800.0),
}


def read_earlinet(code, start, end, station_alt, overlap):
    folder = EARLINET / code
    files = glob.glob(str(folder / f"*{code}*b1064*.nc"))
    # dedup: key = start_end timestamps, keep highest version
    best = {}
    for f in files:
        parts = Path(f).stem.split("_")
        if len(parts) < 9:
            continue
        key = parts[5] + "_" + parts[6]
        try:
            ver = int(parts[7].lstrip("v"))
        except ValueError:
            ver = 0
        if key not in best or ver > best[key][0]:
            best[key] = (ver, f, parts[5])
    # molecular Rayleigh on the uniform AGL grid (US standard atmosphere)
    grid = RANGE_REF
    atm = load_standard_atmosphere(IC.STD_ATM, grid + station_alt)
    mol = calculate_molecular_properties(atm.temperature, atm.pressure, grid, 1064e-9)
    alpha_mol = mol.alpha_mol               # m^-1 (for optical depth in metres)
    d0 = datetime.strptime(start, "%Y%m%d"); d1 = datetime.strptime(end, "%Y%m%d")
    times, atts = [], []
    for key, (_, f, tstr) in sorted(best.items()):
        try:
            ftime = datetime.strptime(tstr, "%Y%m%d%H%M")
        except ValueError:
            continue
        if ftime < d0 or ftime > d1:
            continue
        try:
            with Dataset(f) as nc:
                alt_asl = np.asarray(nc.variables["altitude"][:], "f8").ravel()
                tr = np.asarray(nc.variables["time"][:], "f8").ravel()
                bsc = IC._clean(nc.variables["backscatter"][:]).ravel()    # m^-1 sr^-1
        except Exception:
            continue
        rng = alt_asl - station_alt
        good = np.isfinite(rng) & np.isfinite(bsc)
        if good.sum() < 5:
            continue
        b = np.interp(grid, rng[good], bsc[good], left=np.nan, right=np.nan)  # m^-1 sr^-1
        # overlap: fill below overlap-min with the first valid value above it
        ov = grid < overlap
        if ov.any():
            fill_idx = np.where(~ov & np.isfinite(b))[0]
            if fill_idx.size:
                b[ov] = b[fill_idx[0]]
        ext = np.nan_to_num(b * LR_1064)                      # m^-1
        od = np.concatenate([[0], np.cumsum((ext[1:] + alpha_mol[1:] + ext[:-1] + alpha_mol[:-1]) / 2 * np.diff(grid))])
        trans = np.exp(-od)
        att = (b + mol.beta_mol) * trans * trans * 1e6        # Mm^-1 sr^-1
        times.append(np.datetime64(ftime)); atts.append(att)
    if not times:
        return None
    return dict(time=np.array(times), grid=grid, att=np.array(atts), station_alt=station_alt)


def compare(code, start, end):
    site = SITES[code]
    l2 = IC.read_l2(site["wmo"], site["ident"], start, end)
    if l2 is None:
        return None
    # calibrate CHM (rayleigh)
    cal = IC.load_calib_series(f"{site['wmo']}_{site['ident']}_rayleigh")
    if cal is None:
        return dict(error="no CHM calibration series")
    ck = IC.interp_calib(cal[0], cal[1], l2["time"])
    beta = l2["beta"] * (l2["calc"] / ck)[:, None]            # Mm^-1 sr^-1, on l2['alt'] (ASL)
    ea = read_earlinet(code, start, end, l2["station_alt"], site["overlap"])
    if ea is None:
        return dict(error="no EARLINET profiles in window")
    # match: for each EARLINET time, median CHM within +/-30 min, interp to EARLINET grid (AGL)
    chm_t = pd.to_datetime(l2["time"])
    z_chm_agl = l2["alt"] - l2["station_alt"]
    pairs_e, pairs_c = [], []
    for te, ae in zip(ea["time"], ea["att"]):
        lo = pd.Timestamp(te) - pd.Timedelta(minutes=30); hi = pd.Timestamp(te) + pd.Timedelta(minutes=30)
        sel = (chm_t >= lo) & (chm_t <= hi)
        if sel.sum() == 0:
            continue
        with np.errstate(all="ignore"):
            cprof = np.nanmedian(beta[np.asarray(sel)], axis=0)
        ci = np.interp(ea["grid"], z_chm_agl, cprof, left=np.nan, right=np.nan)
        pairs_e.append(ae); pairs_c.append(ci)
    if not pairs_e:
        return dict(error="no temporal matches")
    E = np.array(pairs_e); C = np.array(pairs_c)
    zmask = (ea["grid"] >= 500) & (ea["grid"] <= 5000)
    s = IC._stats(C, E, zmask)
    s["matched"] = len(pairs_e)
    return s
