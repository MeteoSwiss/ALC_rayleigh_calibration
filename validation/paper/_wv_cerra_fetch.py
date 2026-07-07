"""Fetch CERRA (5.5 km) for the valley/vertical-resolution WV test, winter overlap day 2025-01-15 12 UT.
Reduced to keep the CDS/CERRA queue+download tractable:
  - model-levels (lower 40 of 106, terrain-following, dense in the boundary layer): q, t
  - single-levels: surface_pressure, orography (anchor the model-level height integration)
  - pressure-levels (~15 p-levels): q, t, z  (for the model-vs-pressure comparison)
Full European native grid (no area subsetting). CERRA queues for minutes -> run in background.
"""
import sys
from pathlib import Path

import earthkit.data as ekd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from validation.paper.wv_resolution_lib import OUT  # noqa: E402

ekd.settings.set("cache-policy", "user")
Y, M, D, T = ["2025"], ["01"], ["15"], ["12:00"]
ML = [str(i) for i in range(67, 107)]          # lower 40 model levels (106 = surface)
PL = ["1000", "975", "950", "925", "900", "875", "850", "825", "800", "750", "700", "600", "500", "400", "300"]


def get(dataset, req, out):
    if out.is_file() and out.stat().st_size > 10000:
        print(f"  {out.name} exists, skip", flush=True); return
    print(f"fetching {dataset} -> {out.name} ...", flush=True)
    ekd.from_source("cds", dataset, req).save(str(out))
    print(f"  saved {out.name} ({out.stat().st_size/1e6:.0f} MB)", flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    get("reanalysis-cerra-model-levels",
        dict(variable=["specific_humidity", "temperature"], model_level=ML,
             data_type=["reanalysis"], year=Y, month=M, day=D, time=T, data_format="grib"),
        OUT / "cerra_ml_20250115.grib")
    get("reanalysis-cerra-pressure-levels",
        dict(variable=["specific_humidity", "temperature", "geopotential"], pressure_level=PL,
             data_type=["reanalysis"], product_type="analysis", year=Y, month=M, day=D, time=T, data_format="grib"),
        OUT / "cerra_pl_20250115.grib")
    get("reanalysis-cerra-single-levels",
        dict(variable=["surface_pressure", "orography"], data_type=["reanalysis"], product_type="analysis",
             level_type="surface_or_atmosphere", year=Y, month=M, day=D, time=T, data_format="grib"),
        OUT / "cerra_sfc_20250115.grib")
    print("CERRA_FETCH_DONE")


if __name__ == "__main__":
    main()
