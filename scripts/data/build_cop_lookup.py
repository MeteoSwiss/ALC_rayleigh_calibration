#!/usr/bin/env python3
"""
Build a per-instrument lookup of the operational calibration constant (median of
calibration_constant_0) from the A:/E-PROFILE_L2_monthly files. Used to convert
re-derived lidar constants into dimensionless calibration coefficients
(re-derived C_L / operational C_op, ideal = 1.0) for the MATLAB-vs-Python comparison.

A single WMO folder may host several co-located instruments, distinguished by the
identifier letter (e.g. Payerne A=CL31, B=CHM15k, C=CL61). The lookup is therefore
keyed per instrument '<WMO>_<identifier>' so each instrument's C_op is its own — a
per-WMO median would otherwise mix instruments. This key matches the Python output
folders and the MATLAB rayleigh_<WMO>_<id>.mat files.

Samples up to N_SAMPLE monthly files per instrument (evenly spaced) — C_op is
near-constant per instrument, so this is adequate for a ranked-overview figure.

Output: cop_lookup.json  { "<WMO>_<id>": {"cop_median": .., "unit_factor": .., "n_files": ..} }
"""
import glob
import json
import os
import re
from collections import defaultdict
import numpy as np
from netCDF4 import Dataset

BASE = "A:/E-PROFILE_L2_monthly"
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/cop_lookup.json"
N_SAMPLE = 8


def unit_factor(units: str) -> float:
    m = re.match(r"\s*([0-9.]+[eE][-+]?[0-9]+)", units or "")
    return float(m.group(1)) if m else 1.0


def _identifier(path):
    tok = os.path.basename(path).split("_")[-1].replace(".nc", "")
    return tok[:-6] if len(tok) > 6 else tok[:1]


def main():
    # Group files by (WMO, identifier) = one instrument.
    inst = defaultdict(list)
    for wmo in sorted(os.listdir(BASE)):
        for f in glob.glob(f"{BASE}/{wmo}/*/L2_*.nc"):
            inst[f"{wmo}_{_identifier(f)}"].append(f)

    lookup = {}
    keys = sorted(inst)
    for k, key in enumerate(keys):
        files = sorted(inst[key])
        idx = np.unique(np.linspace(0, len(files) - 1, min(N_SAMPLE, len(files))).astype(int))
        cops, uf = [], 1.0
        for i in idx:
            try:
                d = Dataset(files[i])
                cc = np.asarray(d.variables["calibration_constant_0"][:], dtype="f8")
                cc = cc[np.isfinite(cc) & (cc != 0)]
                if cc.size:
                    cops.append(float(np.median(cc)))
                uf = unit_factor(getattr(d.variables["attenuated_backscatter_0"], "units", ""))
                d.close()
            except Exception:
                continue
        if cops:
            lookup[key] = dict(cop_median=float(np.median(cops)),
                               unit_factor=uf, n_files=len(files))
        if (k + 1) % 50 == 0:
            print(f"  {k+1}/{len(keys)} instruments scanned", flush=True)
    json.dump(lookup, open(OUT, "w"), indent=1)
    print(f"saved {OUT} with {len(lookup)} instruments")


if __name__ == "__main__":
    main()
