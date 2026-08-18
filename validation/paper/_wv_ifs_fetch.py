"""Fetch operational IFS (9 km, ~0.1 deg) humidity profiles over a Europe box covering the 10 WV
stations, via ECMWF Polytope (address polytope.ecmwf.int, EmailKey auth from ~/.polytopeapirc).
IFS (~9 km) is finer than ERA5 (~28 km) -> the extra rung between ERA5 and CERRA in the ladder.
Winter (Jan) + summer (Jul) 2025, valid 00 and 12 UT (00Z run, steps 0 and 12). Small GRIB per season.
"""
import sys
from pathlib import Path

import earthkit.data as ekd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from validation.paper.wv_resolution_lib import OUT  # noqa: E402

AREA = "56/-1/43/14"           # N/W/S/E covering all 10 stations
GRID = "0.1/0.1"               # ~9 km (native IFS)
PL = ["1000", "950", "925", "900", "850", "800", "750", "700", "600", "500"]
# one representative day per season; per-(date,step) requests avoid MARS strict cross-product errors
SEASONS = {"winter": ["20250114", "20250115", "20250116"], "summer": ["20250714", "20250715", "20250716"]}
STEPS = ["0", "12"]            # valid 00 and 12 UT from the 00Z run


def fetch_one(date, step):
    req = {"class": "od", "stream": "oper", "expver": "0001", "type": "fc", "levtype": "pl",
           "levelist": PL, "param": ["133", "130", "129"], "date": date, "time": "0000",
           "step": step, "area": AREA, "grid": GRID}
    return ekd.from_source("polytope", "ecmwf-mars", req, address="polytope.ecmwf.int", stream=False)


def fetch(dates, out):
    if out.is_file() and out.stat().st_size > 10000:
        print(f"  {out.name} exists, skip", flush=True); return
    import earthkit.data as ekd2
    parts = []
    for d in dates:
        for s in STEPS:
            try:
                parts.append(fetch_one(d, s))
                print(f"  got {d} step {s}", flush=True)
            except Exception as e:
                print(f"  MISS {d} step {s}: {str(e)[:90]}", flush=True)
    if not parts:
        print(f"  no IFS for {out.name}"); return
    merged = ekd2.from_source("multi", parts)
    merged.save(str(out))
    print(f"  saved {out.name} ({out.stat().st_size/1e6:.1f} MB)", flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for season, dates in SEASONS.items():
        fetch(dates, OUT / f"ifs_{season}.grib")
    print("IFS_FETCH_DONE")


if __name__ == "__main__":
    main()
