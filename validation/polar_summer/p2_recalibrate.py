# -*- coding: utf-8 -*-
"""p2 - polar-summer evaluation: recalibrate every night of the affected and control streams.

What
  The operational Rayleigh code (calibration.rayleigh.calibrate_rayleigh, options.json + molecular
  method eprof_v2.2 = the published calibration_version 220) is run on native L1
  (D:/E-PROFILE_L1_2026) for every night of
    * the POLAR streams (CHM15k / CL61 north of 56.6 N, where no profile reaches SZA 100 around
      the solstice), 2025-01..2026-08, and
    * CONTROL streams (52-55 N short summer nights, mid- and low-latitude references, two CL61),
      April-September,
  with these variants per night, all from ONE in-memory copy of the night:
    op          operational (SZA > 100), polar fallback off            -> validation, reference
    polar       v2.3 low-sun fallback, candidate defaults (delta 3 deg, >= 4 h, SZA >= 78), QC
                recorded but not applied (polar_max_residual_pct = inf)
    pol_h2/h6   fallback window = 2 h / 6 h centred on solar midnight (delta 0)
    pol_force   fallback forced on nights with < 3 h of darkness (trigger study)
    t98..t90    operational selection with sza_night_threshold lowered to 98, 96, 94, 92, 90
    b80_85 .. b98_100, b100
                the SAME night restricted to SZA bands [80,85) [85,90) [90,94) [94,98) [98,100)
                and SZA > 100 (min_time_range 1 h), on nights that have >= 1 h of darkness ->
                bias versus solar elevation
  Per run: flag, C_L, uncertainty, window, the profiles actually used (count, SZA, T_int, solar
  background bckgrd_rcs_0 for the CHM15k), the fit, the far-range residual (implied window bias,
  calibration.rayleigh.polar.far_range_residual), the molecular-normalised signal X in height bands
  and the far-range photon noise relative to the window signal.

I/O discipline
  Each L1 day file is read ONCE: a per-stream sliding cache assembles the night [D-1, D] in memory
  with exactly load_l1_data's concatenation (checked bit-identical by `selftest`), and every variant
  runs on that copy (preloaded_data). At most 8 worker processes (shared SATA disk).

The package is not modified at run time except three names of calibration.rayleigh.calibration that
are monkeypatched in-process (write -> no-op; window fit -> capture; time filter -> SZA bands),
exactly as validation/chm15k_stability/q16a_recalibration.py does.

Sub-commands:  selftest | run [--workers 8] [--units ...] [--limit N] | time --unit W:I --date D
Outputs: C:/DATA/Projects/202610_v23/polar/runs/<wmo>_<ident>/<first-night>.csv
"""
from __future__ import annotations

import argparse
import copy
import csv
import dataclasses
import json
import os
import re
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta, date as date_type
from pathlib import Path

import numpy as np

# Operational CL61 laser spectrum (910.55 nm, FWHM 0.188 nm): what the published v2.2 files use.
os.environ.setdefault("ALC_WV_SPECTRUM", '{"CL61": [910.55, 0.188]}')

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

OUT = Path("C:/DATA/Projects/202610_v23/polar")
RUN_DIR = OUT / "runs"
L1_ROOT = Path("D:/E-PROFILE_L1_2026")
CAMS = "D:/CAMS_Monthly_04;D:/CAMS_daily"
CENSUS = REPO / "validation/scope_l1_2026_census.json"
EPOCH = datetime(1970, 1, 1)

POLAR_UNITS = [("0-20000-0-01062", "A"), ("0-20000-0-01028", "A"), ("0-20000-0-01001", "A"),
               ("0-20000-0-01010", "A"), ("0-208-0-06011", "A"), ("0-20000-0-01311", "A"),
               ("0-20000-0-03005", "A"), ("0-20000-0-01492", "A"), ("0-20000-0-03018", "A"),
               ("0-20000-0-03023", "A"), ("0-20008-0-BIR", "A")]
# Controls: short summer nights at 52-55 N, mid/low-latitude references, spread over countries,
# plus two CL61 (Birkenes is the only CL61 north of 56.6 N).
CONTROL_UNITS = [("0-20000-0-03162", "A"), ("0-20000-0-10020", "0"), ("0-20000-0-10035", "0"),
                 ("0-20000-0-10393", "0"), ("0-20000-0-06240", "A"), ("0-250-1001-07151", "B"),
                 ("0-20000-0-12843", "A"), ("0-20008-0-UGR", "A"), ("0-20000-0-10393", "C"),
                 ("0-20000-0-06447", "B")]
PERIOD_POLAR = ("20250101", "20260831")
CONTROL_MONTHS = (4, 5, 6, 7, 8, 9)

THRESHOLDS = (98.0, 96.0, 94.0, 92.0, 90.0)
BANDS = {"b80_85": (80.0, 85.0), "b85_90": (85.0, 90.0), "b90_94": (90.0, 94.0),
         "b94_98": (94.0, 98.0), "b98_100": (98.0, 100.0), "b100": (100.0, None)}
POLAR_DEFAULT = dict(polar_night_fallback=True, polar_trigger_dark_hours=0.0, polar_sza_delta=3.0,
                     polar_min_hours=4.0, polar_min_sza=78.0, polar_max_residual_pct=np.inf)
RUN_SPECS = {
    "op": {},
    "polar": dict(POLAR_DEFAULT),
    "pol_h2": dict(POLAR_DEFAULT, polar_sza_delta=0.0, polar_min_hours=2.0),
    "pol_h6": dict(POLAR_DEFAULT, polar_sza_delta=0.0, polar_min_hours=6.0),
    "pol_force": dict(POLAR_DEFAULT, polar_trigger_dark_hours=24.0),
    # the recommended configuration = the CalibrationOptions defaults with the option switched on
    "final": dict(polar_night_fallback=True),
    **{f"t{int(t)}": dict(sza_night_threshold=t) for t in THRESHOLDS},
    **{k: dict(sza_night_threshold=lo, min_time_range=1.0) for k, (lo, hi) in BANDS.items()},
}
SUCCESS = (1.0, 0.5, 0.25)
XBANDS = ((2000, 3000), (3000, 4000), (4000, 5000), (5000, 6000), (6000, 8000), (8000, 10000),
          (10000, 12000), (12000, 14000))
FIELDS = (["wmo", "ident", "type", "date", "run", "sza_max_night", "dark_h", "module", "flag", "c",
           "u", "zb_agl", "zt_agl", "zc_agl", "n_used", "t_mean", "t_first", "t_last", "sza_mean",
           "sza_min", "sza_max", "tint_mean", "bg_mean", "slope", "r2", "rel_err", "err_pct",
           "resid_pct", "resid_q30", "resid_tmed", "x_far", "noise_far_rel", "ls_floor",
           "ls_window_h", "ls_n"]
          + [f"x_{a // 1000}_{b // 1000}" for a, b in XBANDS] + ["cpu_s", "message"])

# --------------------------------------------------------------------------------------------
# Monkeypatching of calibration.rayleigh.calibration (this process only)
# --------------------------------------------------------------------------------------------
STATE = {"band": None, "capture": {}}
_PATCHED = {"done": False}


def profile_sza(data, date_str):
    """SZA of every profile, exactly as filter_time_range computes it (no equation of time)."""
    n = len(data.time_datetime)
    utc_h = np.empty(n)
    for i, dt in enumerate(data.time_datetime):
        utc_h[i] = dt.hour + dt.minute / 60.0 + dt.second / 3600.0
    sh = utc_h + data.longitude / 15.0
    sh = np.where(sh >= 24, sh - 24, np.where(sh < 0, sh + 24, sh))
    current = datetime.strptime(date_str, "%Y%m%d").date()
    doy = current.timetuple().tm_yday
    decl = np.radians(-23.44) * np.cos(2.0 * np.pi * (doy + 10) / 365.25)
    lat = np.radians(data.latitude)
    cz = np.sin(lat) * np.sin(decl) + np.cos(lat) * np.cos(decl) * np.cos(np.radians((sh - 12.0) * 15.0))
    return np.degrees(np.arccos(np.clip(cz, -1.0, 1.0)))


def _subset(data, mask):
    n = len(data.time)
    kw = {}
    for f in dataclasses.fields(data):
        v = getattr(data, f.name)
        if isinstance(v, np.ndarray) and v.ndim >= 1 and v.shape[0] == n and f.name not in (
                "range_alc", "altitude_grid"):
            kw[f.name] = v[mask]
        elif isinstance(v, list) and len(v) == n:
            kw[f.name] = [x for x, k in zip(v, mask) if k]
    return dataclasses.replace(data, **kw)


def _install_patches():
    if _PATCHED["done"]:
        return
    import calibration.rayleigh.calibration as RC
    orig_find = RC.find_optimal_molecular_window
    orig_filter = RC.filter_time_range
    orig_avg = RC.average_ceilometer_data

    def find(**kw):
        res = orig_find(**kw)
        STATE["capture"]["fit"] = (float(res.slope), float(res.r_squared), float(res.relative_error),
                                   float(res.center_range_m))
        return res

    def filt(data, date_str, options):
        out = orig_filter(data, date_str, options)
        band = STATE["band"]
        if band is None or len(out.time) == 0:
            return out
        lo, hi = band
        if hi is None:
            return out
        return _subset(out, profile_sza(out, date_str) < hi)

    def avg(data, *a, **k):
        STATE["capture"]["times"] = np.asarray(data.time, float).copy()
        return orig_avg(data, *a, **k)

    RC.write_calibration_result = lambda *a, **k: None
    RC.find_optimal_molecular_window = find
    RC.filter_time_range = filt
    RC.average_ceilometer_data = avg
    _PATCHED["done"] = True


# --------------------------------------------------------------------------------------------
# Read-once L1 reader: per-file raw items, assembled exactly like load_l1_data
# --------------------------------------------------------------------------------------------
def l1_path(wmo, ident, ds):
    return L1_ROOT / wmo / ds[:4] / ds[4:6] / f"L1_{wmo}_{ident}{ds}.nc"


def read_l1_file(path, itype):
    """One L1 file -> the items load_l1_data's loop body collects (same netCDF4 calls), plus the
    side variables used for statistics only (temp_int in degC, bckgrd_rcs_0)."""
    from netCDF4 import Dataset
    from calibration.io.data_loader import _extract_instrument_metadata
    from calibration.config import InstrumentType
    d = {}
    with Dataset(path, "r") as data:
        d["time"] = data.variables["time"][:]
        d["calendar"] = data.variables["time"].calendar
        d["time_units"] = data.variables["time"].units
        d["range"] = data.variables["range"][:]
        d["range_shape"] = data.variables["range"].shape
        d["altitude"] = float(data.variables["station_altitude"][:])
        d["latitude"] = float(data.variables["station_latitude"][:])
        d["longitude"] = float(data.variables["station_longitude"][:])
        d["rcs"] = data.variables["rcs_0"][:]
        if "linear_depol_ratio" in data.variables:
            d["depol"] = np.ma.filled(np.ma.masked_invalid(
                np.asarray(data.variables["linear_depol_ratio"][:], dtype="f8")), np.nan)
        else:
            d["depol"] = None
        d["cbh"] = data.variables["cloud_base_height"][:]
        n_t = len(d["time"])

        def _hk(*names):
            for nm in names:
                if nm in data.variables:
                    return data.variables[nm][:]
            return np.full(n_t, np.nan)

        is_mini_mpl = itype == InstrumentType.MINI_MPL
        d["temp_om"] = _hk("temperature_optical_module", "temperature_laser")
        if is_mini_mpl:
            for k in ("window", "status_laser", "status_detector", "laser_life", "cal_pulse",
                      "laser_energy"):
                d[k] = np.full(n_t, np.nan)
        else:
            d["window"] = _hk("window_transmission")
            d["status_laser"] = _hk("status_laser", "state_laser")
            d["status_detector"] = _hk("status_detector", "state_detector")
            d["laser_life"] = _hk("laser_life_time")
            d["cal_pulse"] = _hk("calibration_pulse")
            d["laser_energy"] = _hk("laser_energy", "laser_pulse_energy")
        d["vert_vis"] = _hk("vertical_visibility")
        d["meta"] = _extract_instrument_metadata(data, itype)
        # side variables (statistics only)
        tint = _plain(_hk("temp_int"))
        d["side_tint"] = np.where(tint > 150.0, tint - 273.15, tint)      # K -> degC
        d["side_bg"] = _plain(_hk("bckgrd_rcs_0"))
        d["side_time"] = _plain(d["time"])
    return d


def _plain(x):
    """float64 ndarray with masked / fill values as NaN."""
    return np.asarray(np.ma.filled(np.ma.asarray(x, dtype="f8"), np.nan), float)


def assemble(items, itype):
    """Concatenate per-file items exactly like calibration.io.data_loader.load_l1_data."""
    from netCDF4 import num2date
    from calibration.io.data_loader import CeilometerData
    if not items:
        return None
    first = items[0]
    for it in items[1:]:
        if it["range_shape"] != first["range_shape"]:
            raise ValueError("Range shape changed between files")
    time_ = np.concatenate([it["time"] for it in items])
    rcs = np.concatenate([it["rcs"] for it in items], axis=0)
    cbh = np.concatenate([it["cbh"] for it in items], axis=0)
    no_cloud = itype.no_cloud_value
    cbh = np.ma.filled(np.ma.masked_invalid(cbh.astype("f8")), np.nan)
    cbh[~np.isfinite(cbh) | (np.abs(cbh) > 1e30) | (cbh < 0)] = no_cloud
    cat = {k: np.concatenate([it[k] for it in items]) for k in
           ("temp_om", "window", "status_laser", "status_detector", "laser_life", "cal_pulse",
            "vert_vis", "laser_energy")}
    depol_list = [it["depol"] for it in items]
    depol = (np.concatenate(depol_list, axis=0)
             if depol_list and all(x is not None for x in depol_list) else None)
    try:
        time_datetime = list(num2date(time_, first["time_units"], first["calendar"]))
    except OverflowError:
        return None
    hours_since_start = (time_ - time_.min()) * 24
    range_alc = first["range"]
    om_id, serial, firmware = first["meta"]
    return CeilometerData(
        time=time_, time_datetime=time_datetime, hours_since_start=hours_since_start,
        range_alc=range_alc, altitude=first["altitude"], altitude_grid=range_alc + first["altitude"],
        latitude=first["latitude"], longitude=first["longitude"], rcs=rcs, cbh=cbh,
        temperature_optical_module=cat["temp_om"], window_transmission=cat["window"],
        status_laser=cat["status_laser"], status_detector=cat["status_detector"],
        laser_life_time=cat["laser_life"], calibration_pulse=cat["cal_pulse"],
        vertical_visibility=cat["vert_vis"], laser_energy=cat["laser_energy"], depol=depol,
        optical_module_id=om_id, instrument_serial_number=serial,
        instrument_firmware_version=firmware, calendar=first["calendar"],
        time_units=first["time_units"])


class NightReader:
    """Sliding read-once cache: night D = files [D-1, D]; each file is read once per stream."""

    def __init__(self, wmo, ident, itype):
        self.wmo, self.ident, self.itype = wmo, ident, itype
        self.cache = {}
        self.n_read = 0

    def _item(self, ds):
        if ds not in self.cache:
            p = l1_path(self.wmo, self.ident, ds)
            self.cache[ds] = read_l1_file(p, self.itype) if p.exists() else None
            if self.cache[ds] is not None:
                self.n_read += 1
            for k in sorted(self.cache)[:-2]:      # keep only the last two days in memory
                del self.cache[k]
        return self.cache[ds]

    def night(self, ds):
        prev = (datetime.strptime(ds, "%Y%m%d") - timedelta(days=1)).strftime("%Y%m%d")
        items = [x for x in (self._item(prev), self._item(ds)) if x is not None]
        if not items:
            return None, None
        side = dict(time=np.concatenate([x["side_time"] for x in items]),
                    tint=np.concatenate([x["side_tint"] for x in items]),
                    bg=np.concatenate([x["side_bg"] for x in items]))
        return assemble(items, self.itype), side


# --------------------------------------------------------------------------------------------
# Per-night calibration
# --------------------------------------------------------------------------------------------
_CENSUS = {}


def census_entry(wmo, ident):
    if not _CENSUS:
        for s in json.load(open(CENSUS, encoding="utf-8")):
            _CENSUS[(s["wmo"], s["ident"])] = s
    return _CENSUS[(wmo, ident)]


def base_options():
    from calibration import CalibrationOptions, DataLevel
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.folder_output = OUT / "_nowrite"            # never written (write is patched out)
    o.plot_all = False
    o.plot_main = False
    o.molecular_method = "eprof_v2.2"             # calibration_version 220 of the published files
    o.cams_folder = CAMS
    return o


def make_info(s):
    from calibration import InstrumentInfo
    from calibration.config import InstrumentType
    return InstrumentInfo(site_name=s.get("site", s["wmo"]), wmo_id=s["wmo"], identifier=s["ident"],
                          instrument_type=InstrumentType(s["type"]), latitude=s["lat"],
                          longitude=s["lon"], altitude=s.get("alt", 0.0))


def runs_for_night(geo_sza_max, dark_h):
    """Which variants make sense on this night (cheap geometric pre-selection).

    'polar' runs on every night: where a real night exists it must equal 'op' exactly (real-data
    check of the bit-identity); the other fallback windows only differ from 'op' without darkness.
    """
    runs = ["op", "polar"]
    if dark_h <= 0.0:
        runs += ["pol_h2", "pol_h6"]
    if dark_h < 3.0:
        runs.append("pol_force")
    if dark_h < 6.0:
        runs += [f"t{int(t)}" for t in THRESHOLDS if geo_sza_max > t]
    if dark_h >= 1.0:
        runs += [k for k, (lo, hi) in BANDS.items() if geo_sza_max > lo]
    return runs


def residual_variants(fi, slope, zc, band=(8000.0, 12000.0)):
    """Two cirrus-robust alternatives to the night-mean median residual (implied window bias, %):
    q30 = 30th percentile over the band gates de-biased by the measured per-gate photon noise
          (E[q30] of N(delta, s) = delta - 0.524 s), unaffected by a layer covering < 70 % of gates;
    tmed = median over the binned profiles of each profile's median band residual, unaffected by a
          layer present in < 50 % of the profiles."""
    out = {}
    try:
        sig = np.asarray(fi["signal"], float)
        pm = np.asarray(fi["p_mol"], float)
        rng = np.asarray(fi["range_alc"], float)
        okp = np.isfinite(pm) & (pm > 0)
        pmc = float(np.interp(zc, rng[okp], pm[okp]))
        m = (rng >= band[0]) & (rng <= band[1]) & np.isfinite(sig) & okp
        sg = fi.get("sigma_signal")
        if sg is not None and m.sum() >= 10:
            sg = np.asarray(sg, float)
            mm = m & np.isfinite(sg)
            if mm.sum() >= 10:
                d = np.percentile(sig[mm] - slope * pm[mm], 30) + 0.5244 * np.median(sg[mm])
                out["resid_q30"] = 100.0 * d / (slope * pmc)
        st = fi.get("signal_stack")
        if st is not None and np.ndim(st) == 2 and np.shape(st)[0] >= 3:
            st = np.asarray(st, float)
            mb = (rng >= band[0]) & (rng <= band[1]) & okp
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                per = np.nanmedian(st[:, mb] - slope * pm[mb][None, :], axis=1)
            per = per[np.isfinite(per)]
            if per.size >= 3:
                out["resid_tmed"] = 100.0 * float(np.median(per)) / (slope * pmc)
    except Exception:  # noqa: BLE001 - diagnostics only
        pass
    return out


def _side_mean(side, times, key):
    if side is None or times is None or not len(times):
        return np.nan
    idx = np.searchsorted(side["time"], times)
    idx = np.clip(idx, 0, side["time"].size - 1)
    ok = np.abs(side["time"][idx] - times) < 1e-7
    v = side[key][idx[ok]]
    v = v[np.isfinite(v)]
    return float(np.mean(v)) if v.size else np.nan


def calibrate_night(s, date_str, native, side, runs=None):
    _install_patches()
    from calibration import calibrate_rayleigh
    from calibration.rayleigh import polar
    info = make_info(s)
    o_base = base_options()
    geo = polar.night_geometry(date_str, info.latitude, 100.0)
    if runs is None:
        runs = runs_for_night(geo.sza_max, geo.dark_hours)
    rows = []
    module = getattr(native, "optical_module_id", None)
    for run in runs:
        spec = RUN_SPECS[run]
        o = copy.copy(o_base)
        for k, v in spec.items():
            setattr(o, k, v)
        STATE["band"] = BANDS.get(run)
        STATE["capture"] = {}
        fi = {}
        t0 = time.process_time()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = calibrate_rayleigh(date_str, info, o, preloaded_data=native, fit_inputs_out=fi)
            flag, c, u = float(r.flag), float(r.lidar_constant), float(r.uncertainty)
            zb = (r.calibration_bottom_height - info.altitude
                  if r.calibration_bottom_height is not None else np.nan)
            zt = (r.calibration_top_height - info.altitude
                  if r.calibration_top_height is not None else np.nan)
            msg = r.message
        except Exception as exc:  # noqa: BLE001 - one bad variant must not lose the night
            flag, c, u, zb, zt = -99.0, np.nan, np.nan, np.nan, np.nan
            msg = f"{type(exc).__name__}: {exc}"
        cpu = time.process_time() - t0
        ok = flag in SUCCESS
        row = dict(wmo=s["wmo"], ident=s["ident"], type=s["type"], date=date_str, run=run,
                   sza_max_night=round(geo.sza_max, 3), dark_h=round(geo.dark_hours, 3),
                   module=module, flag=flag, c=c if ok else np.nan, u=u if ok else np.nan,
                   zb_agl=zb, zt_agl=zt, cpu_s=round(cpu, 3), message=str(msg)[:160])
        m = re.search(r"method err: ([0-9.]+)%", str(msg))
        row["err_pct"] = float(m.group(1)) if m else np.nan
        times = STATE["capture"].get("times")
        if times is not None and len(times):
            sub = _subset_times(native, times)
            sza = profile_sza(sub, date_str) if sub is not None else np.full(len(times), np.nan)
            row.update(n_used=len(times), t_mean=float(np.mean(times)), t_first=float(times.min()),
                       t_last=float(times.max()), sza_mean=float(np.nanmean(sza)),
                       sza_min=float(np.nanmin(sza)), sza_max=float(np.nanmax(sza)),
                       tint_mean=_side_mean(side, times, "tint"), bg_mean=_side_mean(side, times, "bg"))
        fit = STATE["capture"].get("fit")
        if fit is not None:
            row.update(slope=fit[0], r2=fit[1], rel_err=fit[2], zc_agl=fit[3])
        if fi and fit is not None and np.isfinite(fit[0]) and fit[0] > 0:
            sig = np.asarray(fi["signal"], float)
            pm = np.asarray(fi["p_mol"], float)
            rng = np.asarray(fi["range_alc"], float)
            bias, xf = polar.far_range_residual(sig, pm, rng, fit[0], fit[3])
            row.update(resid_pct=bias, x_far=xf)
            row.update(residual_variants(fi, fit[0], fit[3]))
            for a, b in XBANDS:
                mm = (rng >= a) & (rng < b) & np.isfinite(sig) & (pm > 0)
                if mm.sum() >= 5:
                    row[f"x_{a // 1000}_{b // 1000}"] = float(
                        np.mean(sig[mm]) / (fit[0] * np.mean(pm[mm])))
            sn = fi.get("sigma_native")
            if sn and np.isfinite(fit[3]):
                rn = np.asarray(sn["range_alc"], float)
                sg = np.asarray(sn["sigma"], float)
                mm = (rn >= 10000) & (rn <= 14000) & np.isfinite(sg)
                okp = np.isfinite(pm) & (pm > 0)
                pmc = float(np.interp(fit[3], rng[okp], pm[okp]))
                if mm.sum() >= 5 and pmc > 0:
                    row["noise_far_rel"] = float(np.median(sg[mm]) / (fit[0] * pmc))
        ls = fi.get("low_sun") if fi else None
        if ls is not None:
            # resid_pct / resid_tmed stay the harness's own estimators (night-mean median / time
            # median) for every run, so that dark and low-sun nights are compared like for like
            row.update(ls_floor=ls.sza_floor, ls_window_h=ls.window_hours, ls_n=ls.n_profiles)
        rows.append(row)
    STATE["band"] = None
    return rows


def _subset_times(native, times):
    t = np.asarray(native.time, float)
    idx = np.clip(np.searchsorted(t, times), 0, t.size - 1)
    ok = np.abs(t[idx] - times) < 1e-7
    mask = np.zeros(t.size, bool)
    mask[idx[ok]] = True
    if not mask.any():
        return None
    return _subset(native, mask)


# --------------------------------------------------------------------------------------------
# Tasks
# --------------------------------------------------------------------------------------------
def _write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})
    os.replace(tmp, path)


def work_chunk(wmo, ident, dates, out_path):
    """Worker: one stream, consecutive nights, each L1 file read once -> one CSV (resumable)."""
    if out_path.exists():
        return str(out_path), 0, 0.0, 0
    warnings.filterwarnings("ignore")
    import logging
    logging.disable(logging.CRITICAL)
    from calibration.config import InstrumentType
    s = census_entry(wmo, ident)
    reader = NightReader(wmo, ident, InstrumentType(s["type"]))
    rows = []
    t0 = time.time()
    for d in dates:
        try:
            native, side = reader.night(d)
            if native is None:
                continue
            rows.extend(calibrate_night(s, d, native, side))
        except Exception as exc:  # noqa: BLE001
            rows.append(dict(wmo=wmo, ident=ident, type=s["type"], date=d, run="op", flag=-98,
                             message=f"{type(exc).__name__}: {exc}"[:160]))
    _write_rows(out_path, rows)
    return str(out_path), len(rows), time.time() - t0, reader.n_read


def build_tasks(units=None, block_days=61):
    tasks = []
    if units is None:
        units = [(u, "polar") for u in POLAR_UNITS] + [(u, "control") for u in CONTROL_UNITS]
    for (wmo, ident), kind in units:
        d0 = datetime.strptime(PERIOD_POLAR[0], "%Y%m%d")
        d1 = datetime.strptime(PERIOD_POLAR[1], "%Y%m%d")
        dates = []
        for k in range((d1 - d0).days + 1):
            d = d0 + timedelta(days=k)
            if kind == "control" and d.month not in CONTROL_MONTHS:
                continue
            ds = d.strftime("%Y%m%d")
            if l1_path(wmo, ident, ds).exists():
                dates.append(ds)
        # consecutive blocks (one extra file read per block for its first night's D-1)
        for i in range(0, len(dates), block_days):
            blk = dates[i:i + block_days]
            tasks.append((wmo, ident, blk, RUN_DIR / f"{wmo}_{ident}" / f"{blk[0]}.csv"))
    return tasks


FINAL_DIR = OUT / "runs_final"


def work_final(wmo, ident, dates, out_path):
    """Worker: the recommended configuration ('final') on a list of nights of one stream."""
    warnings.filterwarnings("ignore")
    import logging
    logging.disable(logging.CRITICAL)
    from calibration.config import InstrumentType
    s = census_entry(wmo, ident)
    reader = NightReader(wmo, ident, InstrumentType(s["type"]))
    rows = []
    t0 = time.time()
    for d in dates:
        native, side = reader.night(d)
        if native is not None:
            rows.extend(calibrate_night(s, d, native, side, runs=["final"]))
    _write_rows(out_path, rows)
    return str(out_path), len(rows), time.time() - t0, reader.n_read


def rerun_final(workers):
    """'final' differs from the evaluated 'pol_h6' variant only where the SZA floor shortens the 6-h
    window (Hopen / Bjornoya around the solstice): pol_h6 rejected those nights (window < 6 h),
    'final' accepts a window down to min_time_range. Re-run exactly those nights (read once)."""
    import pandas as pd
    frames = [pd.read_csv(f, dtype={"date": str, "ident": str}) for f in sorted(RUN_DIR.glob("*/*.csv"))]
    df = pd.concat(frames, ignore_index=True)
    rej = df[(df.run == "pol_h6") & (df.flag == -14)]
    tasks = []
    for (wmo, ident), g in rej.groupby(["wmo", "ident"]):
        dates = sorted(g.date.unique())
        tasks.append((wmo, ident, dates, FINAL_DIR / f"{wmo}_{ident}.csv"))
    print(f"final re-run: {sum(len(t[2]) for t in tasks)} nights on {len(tasks)} streams", flush=True)
    with ProcessPoolExecutor(max_workers=min(workers, 8)) as ex:
        for f in as_completed([ex.submit(work_final, *t) for t in tasks]):
            p, n, dt, nr = f.result()
            print(f"{Path(p).name}: {n} rows, {nr} files, {dt:.0f} s", flush=True)


def pool_run(tasks, workers):
    t0 = time.time()
    done = 0
    n_files = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(work_chunk, *t) for t in tasks]
        for f in as_completed(futs):
            try:
                p, n, dt, nr = f.result()
            except Exception:  # noqa: BLE001
                traceback.print_exc()
                continue
            done += 1
            n_files += nr
            print(f"[{done}/{len(tasks)}] {Path(p).parent.name}/{Path(p).name}: {n} rows, "
                  f"{nr} files, {dt:.0f} s (elapsed {time.time() - t0:.0f} s)", flush=True)
    print(f"total L1 files read: {n_files}", flush=True)


# --------------------------------------------------------------------------------------------
# selftest: the read-once assembly is identical to load_l1_data
# --------------------------------------------------------------------------------------------
def selftest():
    from calibration.config import InstrumentType
    from calibration.io.data_loader import load_l1_data
    cases = [("0-20000-0-01010", "A", "20250616"), ("0-20008-0-BIR", "A", "20260612"),
             ("0-20000-0-03005", "A", "20250101")]
    for wmo, ident, ds in cases:
        s = census_entry(wmo, ident)
        it = InstrumentType(s["type"])
        reader = NightReader(wmo, ident, it)
        prev = (datetime.strptime(ds, "%Y%m%d") - timedelta(days=2)).strftime("%Y%m%d")
        reader.night(prev)                                  # warm the sliding cache like a run
        a, _ = reader.night((datetime.strptime(prev, "%Y%m%d") + timedelta(days=1)).strftime("%Y%m%d"))
        a, _ = reader.night(ds)
        prev1 = (datetime.strptime(ds, "%Y%m%d") - timedelta(days=1)).strftime("%Y%m%d")
        files = [p for p in (l1_path(wmo, ident, prev1), l1_path(wmo, ident, ds)) if p.exists()]
        b = load_l1_data(files, it)
        if a is None or b is None:
            print(f"{wmo}_{ident} {ds}: read-once {'None' if a is None else 'data'}, "
                  f"load_l1_data {'None' if b is None else 'data'} ({len(files)} files)")
            continue
        bad = []
        for f in dataclasses.fields(b):
            va, vb = getattr(a, f.name), getattr(b, f.name)
            if isinstance(vb, np.ndarray):
                same = (va.shape == vb.shape and va.dtype == vb.dtype
                        and np.array_equal(np.asarray(va), np.asarray(vb), equal_nan=True)
                        and type(va) is type(vb))
            else:
                same = va == vb
            if not same:
                bad.append(f.name)
        print(f"{wmo}_{ident} {ds}: files read {reader.n_read}; differing fields: {bad or 'none'}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    fr = sub.add_parser("final")
    fr.add_argument("--workers", type=int, default=2)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=8)
    r.add_argument("--units", default="", help="comma list wmo:ident (default: all)")
    r.add_argument("--limit", type=int, default=0)
    t = sub.add_parser("time")
    t.add_argument("--unit", default="0-20000-0-01010:A")
    t.add_argument("--date", default="20250616")
    args = ap.parse_args()
    if args.cmd == "selftest":
        selftest()
    elif args.cmd == "final":
        rerun_final(args.workers)
    elif args.cmd == "run":
        units = None
        if args.units:
            sel = {tuple(u.split(":")) for u in args.units.split(",")}
            units = [(u, "polar") for u in POLAR_UNITS if u in sel] + \
                    [(u, "control") for u in CONTROL_UNITS if u in sel]
        tasks = build_tasks(units)
        if args.limit:
            tasks = tasks[:args.limit]
        print(f"{len(tasks)} tasks, {sum(len(t[2]) for t in tasks)} unit-nights, "
              f"{args.workers} workers", flush=True)
        pool_run(tasks, min(args.workers, 8))
    elif args.cmd == "time":
        from calibration.config import InstrumentType
        wmo, ident = args.unit.split(":")
        s = census_entry(wmo, ident)
        reader = NightReader(wmo, ident, InstrumentType(s["type"]))
        t0 = time.time()
        native, side = reader.night(args.date)
        t1 = time.time()
        rows = calibrate_night(s, args.date, native, side)
        for rw in rows:
            print({k: (round(rw[k], 4) if isinstance(rw.get(k), float) else rw.get(k))
                   for k in ("run", "flag", "c", "n_used", "sza_min", "sza_max", "tint_mean",
                             "bg_mean", "resid_pct", "x_far", "noise_far_rel", "err_pct", "cpu_s")})
        print(f"read {t1 - t0:.1f} s, calibrate {time.time() - t1:.1f} s")


if __name__ == "__main__":
    main()
