#!/usr/bin/env python3
"""
Scan A:/E-PROFILE_L2_monthly and build a manifest of all CHM15k / CL61 / Mini-MPL
*instruments* (the types that support Rayleigh calibration). A single WMO folder may
host several co-located instruments, distinguished by the **identifier** letter in the
filename (e.g. Payerne 0-20000-0-06610: A=CL31, B=CHM15k, C=CL61). Each (WMO, identifier)
pair is therefore one entry, matching the MATLAB per-instrument files rayleigh_<WMO>_<id>.mat.

Saves stations_l2_manifest.json.
"""
import glob
import json
import os
from collections import defaultdict
from netCDF4 import Dataset

BASE = "A:/E-PROFILE_L2_monthly"
WANTED = {"CHM15k", "CL61", "Mini-MPL"}


def _identifier(path):
    """Filename token is '<identifier><YYYYMM>', e.g. 'A202501' or 'j201612'."""
    tok = os.path.basename(path).split("_")[-1].replace(".nc", "")
    return tok[:-6] if len(tok) > 6 else tok[:1]


def main():
    # Group every L2 file by (WMO, identifier) = one physical instrument.
    inst = defaultdict(list)
    for wmo in sorted(os.listdir(BASE)):
        for f in glob.glob(f"{BASE}/{wmo}/*/L2_*.nc"):
            inst[(wmo, _identifier(f))].append(f)

    out = []
    total_months = 0
    for (wmo, ident), files in inst.items():
        files = sorted(files)
        try:
            d = Dataset(files[0])
            itype = str(getattr(d, "instrument_type", "?"))
            if itype not in WANTED:
                d.close()
                continue
            lat = float(d.variables["station_latitude"][:]) if "station_latitude" in d.variables else 0.0
            lon = float(d.variables["station_longitude"][:]) if "station_longitude" in d.variables else 0.0
            alt = float(d.variables["station_altitude"][:]) if "station_altitude" in d.variables else 0.0
            ok = "attenuated_backscatter_0" in d.variables and "calibration_constant_0" in d.variables
            d.close()
        except Exception:
            continue
        if not ok:
            continue
        out.append(dict(wmo=wmo, identifier=ident, itype=itype,
                        lat=lat, lon=lon, alt=alt, n_months=len(files)))
        total_months += len(files)

    out.sort(key=lambda s: (s["itype"], s["wmo"], s["identifier"]))
    json.dump(out, open("C:/DATA/Projects/202606_E-PROFILE_calibration/stations_l2_manifest.json", "w"), indent=1)
    from collections import Counter
    print("instruments:", len(out), "| by type:", dict(Counter(s["itype"] for s in out)))
    print("total month-files:", total_months, "(~%d nights)" % (total_months * 30))
    print("saved stations_l2_manifest.json")


if __name__ == "__main__":
    main()
