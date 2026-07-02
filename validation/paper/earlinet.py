"""
earlinet.py — read EARLINET L2 1064 nm backscatter, convert to attenuated backscatter (per-file
assumed lidar ratio, molecular Rayleigh + two-way transmission), and compare to the colocated CHM15k
ceilometer (operational Python Rayleigh + Kalman calibration). Started as a port of
read_earlinet_att_backscatter.m + the paper_val_earlinet matching (CHM profiles within +/-30 min of
each EARLINET time, median, interpolated to the EARLINET grid; statistics over 500-5000 m AGL), then
hardened beyond the MATLAB for the paper:
  - EARLINET gates below the instrument overlap are excluded (NaN), not constant-filled;
  - the transmission integral extends the lowest trusted extinction to the ground instead of
    assuming an aerosol-free boundary layer;
  - the CHM stream is screened exactly like the station intercomparison (quality flag, clouds,
    fog, +/-15 min expansion) via intercompare.screen();
  - profile times are the averaging-window midpoints (time_bounds), dedup prefers higher QC level.
"""
from __future__ import annotations
import glob
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties
from validation.paper import intercompare as IC

EARLINET = Path("A:/EARLINET")
RANGE_REF = np.arange(0, 15001, 15.0)            # m AGL
LR_1064 = 50.0                                   # fallback lidar ratio [sr] when the file has none
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
    # dedup: key = start_end timestamps; keep the highest (qc level, version)
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
        try:
            qc = int(parts[8].lstrip("qc"))
        except ValueError:
            qc = 0
        rank = (qc, ver)
        if key not in best or rank > best[key][0]:
            best[key] = (rank, f, parts[5])
    # molecular Rayleigh on the uniform AGL grid (US standard atmosphere)
    grid = RANGE_REF
    atm = load_standard_atmosphere(IC.STD_ATM, grid + station_alt)
    mol = calculate_molecular_properties(atm.temperature, atm.pressure, grid, 1064e-9)
    alpha_mol = mol.alpha_mol               # m^-1 (for optical depth in metres)
    d0 = datetime.strptime(start, "%Y%m%d")
    d1 = datetime.strptime(end, "%Y%m%d") + timedelta(days=1)   # include the end day fully
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
                bsc = IC._clean(nc.variables["backscatter"][:]).ravel()    # m^-1 sr^-1
                # profile time = midpoint of the averaging window (time_bounds, s since 1970);
                # fall back to the filename start time if absent or implausible
                tmid = None
                if "time_bounds" in nc.variables:
                    tb = np.asarray(nc.variables["time_bounds"][:], "f8").ravel()
                    if tb.size >= 2 and np.all(np.isfinite(tb[:2])):
                        tmid = np.datetime64("1970-01-01") + np.timedelta64(int(round(tb[:2].mean())), "s")
                        if abs((tmid - np.datetime64(ftime)) / np.timedelta64(1, "h")) > 24:
                            tmid = None
                # lidar ratio the SCC retrieval assumed for this scene (elastic products);
                # fall back to the fixed literature 50 sr
                lr = LR_1064
                if "assumed_particle_lidar_ratio" in nc.variables:
                    v = np.ravel(np.asarray(nc.variables["assumed_particle_lidar_ratio"][:], "f8"))
                    if v.size and np.isfinite(v[0]) and v[0] > 0:
                        lr = float(v[0])
        except Exception:
            continue
        rng = alt_asl - station_alt
        good = np.isfinite(rng) & np.isfinite(bsc)
        if good.sum() < 5:
            continue
        b = np.interp(grid, rng[good], bsc[good], left=np.nan, right=np.nan)  # m^-1 sr^-1
        # Gates below the instrument overlap are unreliable: exclude them (NaN) rather than
        # constant-fill them like the MATLAB did. Nothing is fabricated in the backscatter —
        # NaN gates are simply left out of the comparison.
        b[grid < overlap] = np.nan
        ext = b * lr                                  # m^-1 aerosol extinction, NaN where no data
        fin = np.where(np.isfinite(ext))[0]
        if fin.size == 0:
            continue
        # Two-way transmission: the beam IS attenuated by the boundary-layer aerosol even where
        # the lidar has no trusted retrieval, so for the OD integral only we extend the lowest
        # trusted gate's extinction (clamped >= 0) down to the ground, fill interior gaps
        # linearly, and assume zero above the top gate (free troposphere).
        ext_od = np.interp(grid, grid[fin], ext[fin], left=max(ext[fin[0]], 0.0), right=0.0)
        od = np.concatenate([[0], np.cumsum((ext_od[1:] + alpha_mol[1:] + ext_od[:-1] + alpha_mol[:-1]) / 2 * np.diff(grid))])
        trans = np.exp(-od)
        att = (b + mol.beta_mol) * trans * trans * 1e6        # Mm^-1 sr^-1
        times.append(tmid if tmid is not None else np.datetime64(ftime)); atts.append(att)
    if not times:
        return None
    return dict(time=np.array(times), grid=grid, att=np.array(atts), station_alt=station_alt)


def compare(code, start, end, return_profiles=False):
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
    # paper-consistent screening (same policy as the station intercomparison): quality_flag,
    # clouds (any CBH 0-20 km), fog/vertical visibility, +/-15 min temporal expansion.
    # EARLINET profiles are already cloud-screened by the SCC.
    beta = IC.screen(beta, l2)
    ea = read_earlinet(code, start, end, l2["station_alt"], site["overlap"])
    if ea is None:
        return dict(error="no EARLINET profiles in window")
    # match: for each EARLINET time, median CHM within +/-30 min, interp to EARLINET grid (AGL)
    chm_t = pd.to_datetime(l2["time"])
    z_chm_agl = l2["alt"] - l2["station_alt"]
    pairs_e, pairs_c, pairs_t = [], [], []
    for te, ae in zip(ea["time"], ea["att"]):
        lo = pd.Timestamp(te) - pd.Timedelta(minutes=30); hi = pd.Timestamp(te) + pd.Timedelta(minutes=30)
        sel = (chm_t >= lo) & (chm_t <= hi)
        if sel.sum() == 0:
            continue
        with np.errstate(all="ignore"):
            cprof = np.nanmedian(beta[np.asarray(sel)], axis=0)
        if not np.isfinite(cprof).any():
            continue    # every CHM profile in the window was screened out (clouds/fog/qf)
        ci = np.interp(ea["grid"], z_chm_agl, cprof, left=np.nan, right=np.nan)
        pairs_e.append(ae); pairs_c.append(ci); pairs_t.append(te)
    if not pairs_e:
        return dict(error="no temporal matches")
    E = np.array(pairs_e); C = np.array(pairs_c)
    zmask = (ea["grid"] >= 500) & (ea["grid"] <= 5000)
    s = IC._stats(C, E, zmask)
    s["matched"] = len(pairs_e)
    if return_profiles:
        s["betaE"] = E; s["betaC"] = C; s["grid"] = ea["grid"]; s["times"] = np.array(pairs_t)
    return s


def load_matlab_earlinet(code):
    """relbias/r/N from the MATLAB R_earlinet_<code>.mat (different period; reference only)."""
    import scipy.io as sio
    f = Path("C:/Users/hervo/OneDrive/Documents/MATLAB/ALC/figs_paper_validation") / f"R_earlinet_{code}.mat"
    if not f.is_file():
        return {}
    try:
        m = sio.loadmat(str(f), squeeze_me=True, struct_as_record=False)
        st = m["Re"].stats
        return dict(relbias=float(getattr(st, "relbias_pct", np.nan)), r=float(getattr(st, "r", np.nan)),
                    n=int(getattr(st, "n", 0)))
    except Exception:
        return {}


def run_all(start="20250101", end="20260630"):
    """Compare every EARLINET site to its colocated CHM15k; return rows for the report."""
    rows = []
    for code in ("sir", "lei", "cbw", "ino", "ari"):
        try:
            s = compare(code, start, end)
        except Exception as e:
            s = {"error": repr(e)}
        mm = load_matlab_earlinet(code)
        rows.append((code, s, mm))
        if s and "error" not in s:
            print("  %s: PY relbias=%+.1f%% (med %+.1f%%) r=%.2f (log %.2f) matched=%d | MAT relbias=%+.1f%% r=%.2f n=%d"
                  % (code, s["relbias_pct"], s.get("medrelbias_pct", float("nan")), s["r"],
                     s.get("r_log", float("nan")), s.get("matched", 0),
                     mm.get("relbias", float("nan")), mm.get("r", float("nan")), mm.get("n", 0)), flush=True)
        else:
            print("  %s: %s" % (code, s), flush=True)
    return rows


if __name__ == "__main__":
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    import warnings
    warnings.filterwarnings("ignore")
    run_all()
