"""Fetch ERA5 pressure-level q/t/z (0.25 deg, hourly, independent finer reference) for the box
covering the 10 WV-study stations, winter and summer, 00/12 UT. Cached to OUT/era5_{season}.nc."""
import sys
from pathlib import Path

import earthkit.data as ekd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from validation.paper.wv_resolution_lib import OUT  # noqa: E402

ekd.settings.set("cache-policy", "user")   # persistent cache -> no Windows temp-lock at exit
AREA = [56, -1, 43, 14]                     # N, W, S, E  (covers all 10 stations)
LEVELS = ["200", "225", "250", "300", "350", "400", "450", "500", "550", "600", "650", "700",
          "750", "775", "800", "825", "850", "875", "900", "925", "950", "975", "1000"]
DAYS = [f"{d:02d}" for d in range(10, 21)]  # 11 days per season
SEASONS = {"winter": "202501", "summer": "202507"}


def fetch(ym, out):
    ds = ekd.from_source("cds", "reanalysis-era5-pressure-levels", {
        "product_type": "reanalysis",
        "variable": ["specific_humidity", "temperature", "geopotential"],
        "pressure_level": LEVELS,
        "year": ym[:4], "month": ym[4:6], "day": DAYS, "time": ["00:00", "12:00"],
        "area": AREA, "grid": [0.25, 0.25], "data_format": "netcdf",
    })
    x = ds.to_xarray()
    x.to_netcdf(out)
    print(f"  saved {out}  dims={dict(x.sizes)}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for season, ym in SEASONS.items():
        out = OUT / f"era5_{season}.nc"
        if out.is_file() and out.stat().st_size > 1000:
            print(f"  {out} exists, skip"); continue
        print(f"fetching ERA5 {season} {ym} ...", flush=True)
        fetch(ym, out)
    print("ERA5_FETCH_DONE")


if __name__ == "__main__":
    main()
