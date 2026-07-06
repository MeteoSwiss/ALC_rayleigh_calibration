"""Prefetch ERA5 humidity profiles for every station into a single local cache, for the
ERA5 water-vapour correction (wv_source='era5'). Pulls q/t/z on pressure levels from the DestinE
Earth Data Hub (ARCO/Zarr, lazy point access) for ALL census stations at once, subsampled to a
fixed hourly step, and writes a combined NetCDF the calibration reads offline on the CSCS compute
nodes (which have no internet).

Cache format (read by calibration.cloud.calibration._era5_levels_all_times):
  dims   : station, time, level
  coords : lat(station), lon(station), key(station), time(datetime64), level(pressure hPa)
  vars   : q [kg/kg], t [K], z [geopotential m^2/s^2]

Run on a machine WITH internet (home PC or the CSCS login node), then copy the cache to /scratch.
Needs the EDH token in ~/.destinationearth and: pip install zarr aiohttp xarray netCDF4.

ERA5 is HOURLY on the Hub; the default keeps every hour (--freq 1). A coarser step can be requested
for a quick test, but the production cache is hourly (best temporal fidelity for the WV correction).

Example:
  python scripts/prefetch_era5_edh.py --start 20250101 --end 20260630 --freq 1 \
      --out /scratch/mch/mhrvo/ERA5_WV/era5_wv_2025_2026.nc
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import xarray as xr

EDH_URL = "https://edh:{tok}@api.earthdatahub.destine.eu/era5/reanalysis-era5-pressure-levels-v0.zarr"
REPO = Path(__file__).resolve().parents[1]
DEFAULT_CENSUS = REPO / "validation" / "scope_l1_2026_census.json"


def load_stations(census_path: Path, vaisala_only: bool):
    """Unique (key, lat, lon) station locations from the census. Water vapour only affects the
    910 nm Vaisala instruments, so by default we cache just their locations (deduped by lat/lon)."""
    rows = json.loads(Path(census_path).read_text())
    seen, out = set(), []
    for r in rows:
        if vaisala_only and r.get("type") not in ("CL31", "CL51", "CL61"):
            continue
        lat, lon = round(float(r["lat"]), 3), round(float(r["lon"]), 3)
        if (lat, lon) in seen:
            continue
        seen.add((lat, lon))
        out.append((f"{r['wmo']}_{r['ident']}", lat, lon))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", default=str(DEFAULT_CENSUS))
    ap.add_argument("--start", default="20250101")
    ap.add_argument("--end", default="20260630")
    ap.add_argument("--freq", type=int, default=1, help="hourly step to keep (1 = hourly, ERA5 native; 3 = 3-hourly for a quick test)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0, help="cap n stations (test)")
    ap.add_argument("--all-types", action="store_true", help="cache all station types, not just Vaisala")
    args = ap.parse_args()

    tok = Path.home().joinpath(".destinationearth").read_text().strip()
    stations = load_stations(Path(args.census), vaisala_only=not args.all_types)
    if args.limit:
        stations = stations[: args.limit]
    keys = [s[0] for s in stations]
    lats = np.array([s[1] for s in stations], float)
    lons = np.array([s[2] for s in stations], float)
    print(f"{len(stations)} stations; period {args.start}..{args.end}; {args.freq}-hourly", flush=True)

    t0 = time.time()
    ds = xr.open_dataset(EDH_URL.format(tok=tok), chunks={}, engine="zarr")
    tcoord = "valid_time" if "valid_time" in ds.coords else "time"
    lcoord = "isobaricInhPa" if "isobaricInhPa" in ds.coords else "pressure_level"
    start = f"{args.start[:4]}-{args.start[4:6]}-{args.start[6:8]}"
    end = f"{args.end[:4]}-{args.end[4:6]}-{args.end[6:8]}"
    sel_lat = xr.DataArray(lats, dims="station")
    sel_lon = xr.DataArray(lons, dims="station")
    sub = (ds[["q", "t", "z"]]
           .sel({tcoord: slice(start, end)})
           .sel(latitude=sel_lat, longitude=sel_lon, method="nearest"))
    hrs = sub[tcoord].dt.hour
    sub = sub.sel({tcoord: hrs.isin(list(range(0, 24, args.freq)))})
    print(f"  lazy shape: {dict(sub.sizes)}; opening took {time.time()-t0:.1f}s. Loading ...", flush=True)
    t1 = time.time()
    sub = sub.load()
    print(f"  loaded in {time.time()-t1:.0f}s", flush=True)

    # station axis is dim 'station'; ensure (station, time, level) order
    q = sub["q"].transpose("station", tcoord, lcoord).values.astype("float32")
    t = sub["t"].transpose("station", tcoord, lcoord).values.astype("float32")
    z = sub["z"].transpose("station", tcoord, lcoord).values.astype("float32")
    out = xr.Dataset(
        {"q": (("station", "time", "level"), q),
         "t": (("station", "time", "level"), t),
         "z": (("station", "time", "level"), z)},
        coords={"lat": ("station", lats), "lon": ("station", lons),
                "key": ("station", keys),
                "time": ("time", sub[tcoord].values),
                "level": ("level", np.asarray(sub[lcoord].values, float))})
    out.attrs["source"] = "ERA5 pressure levels via DestinE Earth Data Hub"
    out.attrs["built"] = f"{args.start}..{args.end} {args.freq}-hourly"
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    enc = {v: {"zlib": True, "complevel": 4} for v in ("q", "t", "z")}
    out.to_netcdf(args.out, encoding=enc)
    mb = Path(args.out).stat().st_size / 1e6
    print(f"  saved {args.out} ({mb:.0f} MB, {len(sub[tcoord])} times, {len(out.level)} levels) "
          f"total {time.time()-t0:.0f}s", flush=True)
    print("ERA5_PREFETCH_DONE")


if __name__ == "__main__":
    main()
