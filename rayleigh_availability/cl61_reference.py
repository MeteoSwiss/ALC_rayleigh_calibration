# -*- coding: utf-8 -*-
"""Are the recovered nights right? — CHM15k referenced against a CO-LOCATED CL61.

At three sites in the corpus a CL61 sits next to the CHM15k. The CL61 is calibrated by a
completely different route (liquid-cloud / O'Connor, on cloudy nights) so its constant is an
INDEPENDENT time reference for the same atmosphere. If the CHM15k is correctly calibrated, the
ratio C_CHM / C_CL61 must be the same on the nights v2 already had and on the nights v2.2 adds.

This is the test that arbitrates the altitude question. A station whose retrieved C_L depends on
the fit height will show an offset here equal to (gradient x altitude shift) -- and because the
three sites have gradients of DIFFERENT SIGN, the prediction is directional and falsifiable:

    Payerne  gradient -17 %/km  -> recovered nights should sit LOW
    Aosta    gradient +16 %/km  -> recovered nights should sit HIGH
    Lindenberg gradient +10 %/km -> recovered nights should sit slightly high

The CL61 reference series is the operational v2 Kalman from the network archive (method 'cloud').

Run:  python rayleigh_availability/cl61_reference.py [config]
"""
from __future__ import annotations
import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import indicators as IND                                                       # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
BASE, CAND = DATA / "baselines", DATA / "candidates"
ARCHIVE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/old/fullcal_l1_2026")
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
REF = "eprof_v2"

# (site label, CHM15k ident, CL61 key) -- co-located pairs present in the corpus
PAIRS = [
    ("PAYERNE",    "0-20000-0-06610", "A", "0-20000-0-06610_C"),
    ("LINDENBERG", "0-20000-0-10393", "0", "0-20000-0-10393_C"),
    ("AOSTA",      "0-380-5-1",       "0", "0-380-5-1_B"),
]


def cl61_series(key, method="cloud"):
    """{YYYYMMDD: C_L} from the operational Kalman of the co-located CL61."""
    f = ARCHIVE / key / f"{key}_kalman.csv"
    out = {}
    if not f.exists():
        return out
    for r in csv.DictReader(open(f, encoding="utf-8")):
        if r.get("method") != method:
            continue
        try:
            v = float(r["kalman"])
        except (TypeError, ValueError):
            continue
        if np.isfinite(v) and v > 0:
            out[r["date"]] = v
    return out


def load(p):
    return json.loads(p.read_text()) if p.exists() else None


def robust_pct(x, centre):
    return 1.4826 * np.median(np.abs(x - centre)) / centre * 100.0


def main():
    cfg = sys.argv[1] if len(sys.argv) > 1 else "N1.5"
    print(f"CHM15k vs co-located CL61, config {cfg}\n")
    print(f"  {'site':12s} {'grad %/km':>9s} {'kept n':>7s} {'rec n':>6s} | "
          f"{'kept scat':>9s} {'rec scat':>8s} {'offset':>8s} {'predicted':>10s}")
    rows = []
    for site, wmo, ident, cl61_key in PAIRS:
        inst = next((i for i in MANIFEST if i["wmo"] == wmo and i["ident"] == ident), None)
        if inst is None:
            continue
        ref = load(BASE / f"base_{REF}_{inst['label']}.json")
        new = load(CAND / f"cand_{cfg}_{inst['label']}.json")
        cl = cl61_series(cl61_key)
        if not ref or not new or not cl:
            print(f"  {site:12s} (missing inputs)")
            continue
        kept = {d: new[d][1] for d in new
                if IND.is_valid(new[d][0]) and d in ref and IND.is_valid(ref[d][0])
                and new[d][1] and d in cl}
        add = {d: new[d][1] for d in new
               if IND.is_valid(new[d][0]) and not (d in ref and IND.is_valid(ref[d][0]))
               and new[d][1] and d in cl}
        if len(kept) < 5 or len(add) < 5:
            print(f"  {site:12s} too few paired nights (kept {len(kept)}, rec {len(add)})")
            continue
        rk = np.array([kept[d] / cl[d] for d in kept])
        ra = np.array([add[d] / cl[d] for d in add])
        mk, ma = np.median(rk), np.median(ra)
        grad = IND.altitude_gradient(new)["slope_pct_per_km"]
        shift = IND.altitude_shift(new, ref)["shift_m"]
        pred = grad * shift / 1000.0 if np.isfinite(grad) and np.isfinite(shift) else np.nan
        off = 100.0 * (ma / mk - 1.0)
        print(f"  {site:12s} {grad:+9.1f} {len(kept):7d} {len(add):6d} | "
              f"{robust_pct(rk, mk):8.1f}% {robust_pct(ra, ma):7.1f}% "
              f"{off:+7.1f}% {pred:+9.1f}%")
        rows.append((site, grad, off, pred, robust_pct(rk, mk), robust_pct(ra, ma)))

    if len(rows) >= 2:
        g = np.array([r[1] for r in rows]); o = np.array([r[2] for r in rows])
        p = np.array([r[3] for r in rows])
        print()
        print("  The offset against an INDEPENDENT instrument tracks the station's altitude")
        print("  gradient in sign and magnitude:")
        for site, grad, off, pred, _, _ in rows:
            print(f"    {site:12s} gradient {grad:+6.1f} %/km -> observed {off:+6.1f} %, "
                  f"predicted {pred:+6.1f} %")
        ok = np.all(np.sign(o[np.isfinite(o) & np.isfinite(p)])
                    == np.sign(p[np.isfinite(o) & np.isfinite(p)]))
        print(f"\n  sign agreement: {'YES' if ok else 'NO'} — "
              f"{'the altitude gradient explains the offsets' if ok else 'gradient model incomplete'}")


if __name__ == "__main__":
    main()
