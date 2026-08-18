"""
scope_cloud_2026.py — enumerate every CL31 / CL51 / CL61 instrument stream in the 2026 L1 archive
(per WMO x identifier), with type + lat/lon/alt + first/last date + n_days. Writes
validation/scope_cloud_2026.json — the manifest for the cloud-calibration optimization sweep.
"""
from __future__ import annotations
import collections
import glob
import json
import os
import re
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
from netCDF4 import Dataset

ROOT = Path("D:/E-PROFILE_L1_2026")
TARGET = {"CL31", "CL51", "CL61"}
REPO = Path(__file__).resolve().parents[1]
DATE_RE = re.compile(r"(\d{8})\.nc$")


def scalar(d, *names, default=0.0):
    for n in names:
        if n in d.variables:
            a = np.asarray(d.variables[n][:]).ravel()
            if a.size and np.isfinite(a[0]):
                return float(a[0])
    return default


def main():
    manifest = []
    for wmo in sorted(os.listdir(ROOT)):
        files = glob.glob(str(ROOT / wmo / "2026" / "*" / f"L1_{wmo}_*.nc"))
        if not files:
            continue
        by_id = collections.defaultdict(list)
        for f in files:
            rest = os.path.basename(f)[len(f"L1_{wmo}_"):-3]
            by_id[rest[:-8]].append(f)
        for ident, fs in by_id.items():
            fs = sorted(fs)
            try:
                with Dataset(fs[0]) as d:
                    it = str(getattr(d, "instrument_type", "?"))
                    if it not in TARGET:
                        continue
                    lat = scalar(d, "station_latitude", "latitude")
                    lon = scalar(d, "station_longitude", "longitude")
                    alt = scalar(d, "station_altitude", "altitude")
                    site = str(getattr(d, "site_location", getattr(d, "wigos_station_id", wmo))).split(",")[0]
            except Exception:
                continue
            dates = sorted(m.group(1) for m in (DATE_RE.search(os.path.basename(f)) for f in fs) if m)
            if not dates:
                continue
            manifest.append(dict(wmo=wmo, ident=ident, type=it, group=it, site=site,
                                 lat=lat, lon=lon, alt=alt,
                                 first=dates[0], last=dates[-1], n_days=len(dates), label=f"{wmo}_{ident}"))
    (REPO / "validation" / "scope_cloud_2026.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    c = collections.Counter(m["group"] for m in manifest)
    print(f"cloud manifest: {len(manifest)} streams -> {dict(c)}")
    # top 10 by n_days per type = the phase-1 subset
    for t in ("CL31", "CL51", "CL61"):
        top = sorted([m for m in manifest if m["group"] == t], key=lambda m: -m["n_days"])[:10]
        print(f"  {t}: phase-1 picks (top {len(top)} by n_days): "
              + ", ".join(f"{m['site'][:10]}({m['n_days']})" for m in top[:10]))


if __name__ == "__main__":
    main()
