# -*- coding: utf-8 -*-
"""Instrument-type changes INSIDE a stream (same WIGOS id + ident) over the release window.

The census carries one type per stream, taken from the newest L1 file (scripts/refresh_census.py), and
the runner calibrates the whole record with it. A unit swap under the same ident (Libus and Tusimice:
CL31 -> CL51 in July 2026, found in the v2.3 release comparison) is therefore calibrated with the wrong
type on one side of the swap. This scan reads the `instrument_type` attribute (and the range extent) of
the first file of every month of every census stream in the local L1 archive and lists the streams whose
type changes, with the month of the change.

Output: C:/DATA/Projects/202611_v23_release/type_changes.csv (one row per stream-month) and a summary on
stdout.  Run:  python validation/v23_release/scan_type_changes.py
"""
from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOTS = ["D:/E-PROFILE_L1_2025", "D:/E-PROFILE_L1_2026"]
OUT = Path("C:/DATA/Projects/202611_v23_release/type_changes.csv")


def main():
    c = json.loads((REPO / "validation" / "scope_l1_2026_census.json").read_text(encoding="utf-8"))
    rows = c if isinstance(c, list) else c.get("streams", [])
    rec = []
    for r in rows:
        wmo, ident = r["wmo"], r["ident"]
        months = {}
        for root in ROOTS:
            for f in glob.glob(f"{root}/{wmo}/*/*/L1_{wmo}_{ident}2*.nc"):
                ym = os.path.basename(f)[-11:-5]
                if ym not in months or os.path.basename(f) < os.path.basename(months[ym]):
                    months[ym] = f
        for ym, f in sorted(months.items()):
            try:
                ds = netCDF4.Dataset(f)
                t = str(getattr(ds, "instrument_type", "?"))
                rng = ds["range"][:] if "range" in ds.variables else np.array([np.nan])
                rmax = float(np.nanmax(rng))
                ngates = int(rng.size)
                fw = str(getattr(ds, "instrument_firmware_version", ""))
                ds.close()
            except Exception as exc:  # noqa: BLE001
                t, rmax, ngates, fw = f"ERR {type(exc).__name__}", np.nan, 0, ""
            rec.append(dict(key=f"{wmo}_{ident}", site=r.get("site"), census_type=r["type"], month=ym, itype=t,
                            max_range_m=rmax, n_gates=ngates, firmware=fw, file=os.path.basename(f)))
    df = pd.DataFrame(rec)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"{len(df)} stream-months, {df.key.nunique()} streams")
    chg = []
    for k, g in df.groupby("key"):
        g = g.sort_values("month")
        types = [t for t in g.itype if not t.startswith("ERR")]
        if len(set(types)) > 1:
            first = g.itype.iloc[0]
            flip = g[g.itype != first].month.iloc[0]
            chg.append((k, g.site.iloc[0], g.census_type.iloc[0], " -> ".join(dict.fromkeys(types)), flip))
    print(f"streams whose instrument_type changes inside the record: {len(chg)}")
    for x in chg:
        print("  ", x)
    gates = []
    for k, g in df.groupby("key"):
        ng = set(int(x) for x in g.n_gates if x > 0)
        if len(ng) > 1:
            gates.append((k, g.site.iloc[0], g.census_type.iloc[0], sorted(ng)))
    print(f"streams whose number of range gates changes inside the record: {len(gates)}")
    for x in gates:
        print("  ", x)


if __name__ == "__main__":
    main()
