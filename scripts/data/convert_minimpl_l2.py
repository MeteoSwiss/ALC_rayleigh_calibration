#!/usr/bin/env python3
"""
Reconstruct daily L1 files for the Toulouse miniMPL from A: L2-monthly files.

The E-PROFILE L2 product stores attenuated_backscatter_0 = rcs_0 / calibration_constant_0.
To recover the original range-corrected signal we MULTIPLY BACK by the original
(operational) calibration constant:

        rcs_0[t, z] = attenuated_backscatter_0[z, t] * calibration_constant_0[t]

This restores exactly the L1 range-corrected signal (in the L2's native units, whose
calibration_constant scale ~1e5 matches the report's miniMPL C_L). The reconstructed
daily files are written in the layout the stock loader expects so the unmodified
pipeline (in either worktree) can read them.

Mapping:
  altitude (ASL)            -> range = altitude - station_altitude
  cloud_base_height fill    -> -999.9  (Mini-MPL no-cloud sentinel)
  temperature_optical_module-> NaN     (required by load_l1_data; not used here)

Usage:  python convert_minimpl_l2.py --start 202201 --end 202512
"""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

WMO = "0-20000-0-07617"
IDENT = "A"
SRC = Path("A:/E-PROFILE_L2_monthly") / WMO
OUT_ROOT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/_minimpl_L1")
NO_CLOUD = -999.9

# Unit reconciliation: attenuated_backscatter_0 is stored in "1E-6/(m*sr)" and
# calibration_constant_0 in "1E6*...". Their literal product (atten*C) lands 1e6 above
# the stored rcs scale the operational pipeline / report use (C_L ~1e5, matching
# calibration_constant_0). Multiplying by 1e-6 puts the reconstructed rcs — and hence the
# re-derived lidar constant — on the report's scale. (AOD / with-vs-without are scale-free.)
RCS_UNIT_FACTOR = 1e-6


def _month_files(start: str, end: str) -> list[Path]:
    files = []
    for f in sorted(SRC.glob(f"*/L2_{WMO}_{IDENT}*.nc")):
        ym = "".join(ch for ch in f.stem.split("_")[-1] if ch.isdigit())[-6:]  # YYYYMM
        if start <= ym <= end:
            files.append(f)
    return files


def _write_day(out_path: Path, time_d, range_m, rcs, cbh, lat, lon, alt, units, cal):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds = Dataset(out_path, "w", format="NETCDF4")
    ds.instrument_type = "Mini-MPL"
    ds.wigos_station_id = WMO
    ds.history = "Reconstructed L1 (rcs = attenuated_backscatter_0 * calibration_constant_0) from L2-monthly"
    ds.createDimension("time", len(time_d))
    ds.createDimension("range", len(range_m))
    ds.createDimension("layer", cbh.shape[1])

    vt = ds.createVariable("time", "f8", ("time",))
    vt.units = units
    vt.calendar = cal
    vt[:] = time_d

    vr = ds.createVariable("range", "f4", ("range",))
    vr.units = "m"
    vr[:] = range_m

    vrcs = ds.createVariable("rcs_0", "f8", ("time", "range"), zlib=True)
    vrcs.long_name = "range corrected signal (reconstructed = atten_bsc * calib_const)"
    vrcs[:] = rcs

    vc = ds.createVariable("cloud_base_height", "f8", ("time", "layer"),
                           fill_value=NO_CLOUD, zlib=True)
    vc[:] = cbh

    vtom = ds.createVariable("temperature_optical_module", "f8", ("time",), fill_value=np.nan)
    vtom[:] = np.full(len(time_d), np.nan)

    for name, val in [("station_latitude", lat), ("station_longitude", lon),
                      ("station_altitude", alt)]:
        v = ds.createVariable(name, "f8")
        v[:] = val
    ds.close()


def convert(start: str, end: str) -> None:
    files = _month_files(start, end)
    print(f"{len(files)} monthly files in [{start},{end}]")
    n_days = 0
    for mf in files:
        d = Dataset(mf)
        tvar = d.variables["time"]
        units, cal = tvar.units, getattr(tvar, "calendar", "standard")
        tnum = np.array(tvar[:])
        alt = np.array(d.variables["altitude"][:])
        station_alt = float(d.variables["station_altitude"][:])
        lat = float(d.variables["station_latitude"][:])
        lon = float(d.variables["station_longitude"][:])
        range_m = alt - station_alt

        atten = d.variables["attenuated_backscatter_0"][:]      # (altitude, time) masked
        cc = np.array(d.variables["calibration_constant_0"][:])  # (time,)
        cbh = d.variables["cloud_base_height"][:]                # (time, layer) masked
        d.close()

        # rcs[t, z] = atten[z, t] * C[t]  (x unit factor -> report/operational scale)
        atten_f = np.ma.filled(atten.astype("f8"), np.nan)       # (alt, time)
        rcs = (atten_f.T) * cc[:, None] * RCS_UNIT_FACTOR        # (time, range)

        # cbh: masked / fill / nan -> NO_CLOUD sentinel
        cbh_f = np.ma.filled(cbh.astype("f8"), np.nan)
        cbh_f[~np.isfinite(cbh_f) | (cbh_f > 1e30)] = NO_CLOUD

        dts = num2date(tnum, units, cal)
        ymd = np.array([f"{x.year:04d}{x.month:02d}{x.day:02d}" for x in dts])
        for day in np.unique(ymd):
            m = ymd == day
            out = OUT_ROOT / WMO / day[:4] / day[4:6] / f"L1_{WMO}_{IDENT}{day}.nc"
            _write_day(out, tnum[m], range_m, rcs[m], cbh_f[m], lat, lon, station_alt, units, cal)
            n_days += 1
        print(f"  {mf.name}: {int(np.unique(ymd).size)} days", flush=True)
    print(f"done: {n_days} daily files -> {OUT_ROOT/WMO}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="202201", help="YYYYMM")
    ap.add_argument("--end", default="202512", help="YYYYMM")
    args = ap.parse_args()
    convert(args.start, args.end)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
