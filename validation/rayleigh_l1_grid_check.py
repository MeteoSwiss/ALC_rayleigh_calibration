"""
rayleigh_l1_grid_check.py — is the L1<->L2 Rayleigh divergence a method/grid interaction?
For a few CHM15k streams (1064 nm -> no CAMS needed) over several days, compare the success flag of:
  L1_v11_native : v1.1 on native L1 (1024 bins x 15 s)
  L1_v2_native  : v2   on native L1
  L1_v2_binned  : v2   on L1 pre-binned to the L2 grid (30 m x 300 s)
  L2_v2         : v2   on L2 (30 m x 300 s as delivered)
Hypothesis: v1.1 robust on native L1; v2 rejects native L1 but works once binned to the L2 grid,
and then agrees with L2. Prints per-combo valid fraction + median C_L.
"""
from __future__ import annotations
import json, logging, os, sys, warnings, tempfile
from datetime import datetime, timedelta
from pathlib import Path
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
from netCDF4 import Dataset
from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo
from calibration.config import InstrumentType, DataLevel

REPO = Path(__file__).resolve().parents[1]
SCOPE = json.loads((REPO / "validation" / "scope_network_2026.json").read_text())
L1ROOT = Path("D:/E-PROFILE_L1_2026"); L2ROOT = Path("D:/E-PROFILE_L2_2026")
NSTREAM = int(os.environ.get("RC_NTYPE", "5"))
DAYSTEP = int(os.environ.get("RC_DAYSTEP", "11"))
COMBOS = [
    ("L1_v11_native", "L1", L1ROOT, DataLevel.L1, "eprof_v1.1", None, None),
    ("L1_v2_native",  "L1", L1ROOT, DataLevel.L1, "eprof_v2",  None, None),
    ("L1_v2_binned",  "L1", L1ROOT, DataLevel.L1, "eprof_v2",  30.0, 300.0),
    ("L2_v2",         "L2", L2ROOT, DataLevel.L2_DAILY, "eprof_v2", None, None),
]


def info_for(wmo, ident, f):
    with Dataset(f) as ds:
        return InstrumentInfo(site_name=wmo, wmo_id=wmo, identifier=ident, instrument_type=InstrumentType.CHM15k,
                              latitude=float(ds.variables["station_latitude"][...]),
                              longitude=float(ds.variables["station_longitude"][...]),
                              altitude=float(ds.variables["station_altitude"][...]))


def run_one(level, root, dl, wmo, ident, ds, method, ar, at):
    f = root / wmo / "2026" / ds[4:6] / f"{level}_{wmo}_{ident}{ds}.nc"
    if not f.exists():
        return None
    o = CalibrationOptions.from_json(str(REPO / "options.json"))
    o.data_level = dl; o.folder_root = root; o.molecular_method = method
    o.average_range_m = ar; o.average_time_s = at
    o.abs_cs_lookup_table = Path(""); o.apply_wv_correction = False
    o.folder_output = Path(tempfile.mkdtemp()); o.plot_main = False; o.plot_all = False
    try:
        r = calibrate_rayleigh(ds, info_for(wmo, ident, f), o)
        return (int(r.flag), float(r.lidar_constant))
    except Exception:
        return (-99, float("nan"))


def main():
    warnings.filterwarnings("ignore"); logging.getLogger().setLevel(logging.CRITICAL)
    chm = sorted([m for m in SCOPE if m.get("group") == "CHM15k"], key=lambda m: -m.get("n_days", 0))[:NSTREAM]
    print(f"{len(chm)} CHM15k streams, daystep={DAYSTEP}")
    agg = {c[0]: {"ok": 0, "tot": 0, "cl": []} for c in COMBOS}
    for m in chm:
        wmo, ident = m["wmo"], m["ident"]
        d0 = datetime.strptime(m["first"], "%Y%m%d"); d1 = datetime.strptime(m["last"], "%Y%m%d")
        d = d0
        while d <= d1:
            ds = d.strftime("%Y%m%d"); d += timedelta(days=DAYSTEP)
            for name, level, root, dl, method, ar, at in COMBOS:
                res = run_one(level, root, dl, wmo, ident, ds, method, ar, at)
                if res is None:
                    continue
                flag, cl = res
                agg[name]["tot"] += 1
                if flag == 1 and cl > 0:
                    agg[name]["ok"] += 1; agg[name]["cl"].append(cl)
        print(f"  done {wmo}", flush=True)
    print(f"\n{'combo':16s} {'valid/total':>12s} {'valid%':>7s} {'median C_L':>12s}")
    for name in (c[0] for c in COMBOS):
        a = agg[name]; pct = 100.0 * a["ok"] / a["tot"] if a["tot"] else float("nan")
        med = np.median(a["cl"]) if a["cl"] else float("nan")
        print(f"{name:16s} {a['ok']:5d}/{a['tot']:<6d} {pct:6.0f}% {med:12.3e}")
    print("RAYLEIGH_GRID_CHECK_DONE")


if __name__ == "__main__":
    main()
