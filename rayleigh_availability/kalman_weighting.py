# -*- coding: utf-8 -*-
"""Does per-night uncertainty weighting make the extra availability safe?

The nights v2.2 adds are not as good as the nights v2 already had, and the pipeline KNOWS it: their
median reported uncertainty is 12.4 % against 5.0 % for the retained nights. The operational Kalman
throws that information away -- its measurement variance is one scalar for the whole series, so a
25 % night pulls the best estimate exactly as hard as a 5 % one. That is defensible while the gates
admit only the cleanest nights; it stops being defensible the moment availability is raised.

This compares three best estimates per stream:

    v2            the operational series (fewer nights, one scalar variance)
    v2.2 flat     all the added nights, still one scalar variance
    v2.2 weighted all the added nights, measurement variance ~ the night's own uncertainty

and reports (a) how far each moves from the v2 estimate on the days v2 covers -- movement there is
pure risk, since those days already had an answer -- and (b) how many days gain an estimate at all.
At the two sites with a co-located CL61 the same three are also scored against that INDEPENDENT
reference, which is the only measure here that can see a coherent bias.

Run:  python rayleigh_availability/kalman_weighting.py [config]
"""
from __future__ import annotations
import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import indicators as IND                                                        # noqa: E402
from monitoring.kalman import kalman_best_estimate                              # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
BASE, CAND = DATA / "baselines", DATA / "candidates"
ARCHIVE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/old/fullcal_l1_2026")
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
CL61 = {"PAYERNE_CHM15k_A": "0-20000-0-06610_C", "LINDENBERG_CHM15k_0": "0-20000-0-10393_C"}


def load(p):
    return json.loads(p.read_text()) if p.exists() else None


def triples(rec):
    """(datetimes, values, uncertainties) for the valid nights of a record."""
    t, v, u = [], [], []
    for d in sorted(rec):
        f, c, unc = rec[d][0], rec[d][1], rec[d][2]
        if IND.is_valid(f) and c and c > 0:
            t.append(datetime.strptime(d, "%Y%m%d"))
            v.append(float(c))
            u.append(float(unc) if unc and unc > 0 else np.nan)
    return t, np.array(v), np.array(u)


def as_map(res):
    """{date: state} from a kalman_best_estimate result."""
    if not len(res[0]):
        return {}
    days = res[0].astype("datetime64[D]").astype(str)
    return {d.replace("-", ""): float(s) for d, s in zip(days, res[1])}


def cl61_series(key):
    f = ARCHIVE / key / f"{key}_kalman.csv"
    out = {}
    if not f.exists():
        return out
    for r in csv.DictReader(open(f, encoding="utf-8")):
        if r.get("method") != "cloud":
            continue
        try:
            v = float(r["kalman"])
        except (TypeError, ValueError):
            continue
        if np.isfinite(v) and v > 0:
            out[r["date"]] = v
    return out


def main():
    cfg = sys.argv[1] if len(sys.argv) > 1 else "N2.5"
    print(f"Kalman best estimate: v2 vs {cfg} flat vs {cfg} uncertainty-weighted\n")
    print(f"{'stream':28s} {'nights':>13s} {'days':>13s} | "
          f"{'move flat':>9s} {'move wtd':>9s} | {'outliers op / local':>22s}")
    rows = []
    for inst in MANIFEST:
        if inst["group"] != "CHM15k":
            continue
        ref = load(BASE / f"base_eprof_v2_{inst['label']}.json")
        new = load(CAND / f"cand_{cfg}_{inst['label']}.json")
        if not ref or not new:
            continue
        t0, v0, _ = triples(ref)
        t1, v1, u1 = triples(new)
        if len(t0) < 10 or len(t1) < 10:
            continue
        k0 = kalman_best_estimate(t0, v0, return_rejected=True)
        k1 = kalman_best_estimate(t1, v1, return_rejected=True)
        k2 = kalman_best_estimate(t1, v1, uncertainties=u1, return_rejected=True)
        # The operational screen is a SINGLE IQR over the whole series (improve_alc_calib
        # flag_outliers_lom), not the dashboard's rolling window -- and the user's headline
        # indicator is the operational count, so report that one too.
        g0 = kalman_best_estimate(t0, v0, outlier_mode="global", return_rejected=True)
        g1 = kalman_best_estimate(t1, v1, outlier_mode="global", return_rejected=True)
        m0, m1, m2 = as_map(k0), as_map(k1), as_map(k2)
        both = sorted(set(m0) & set(m1) & set(m2))
        if len(both) < 10:
            continue
        d1 = np.array([abs(m1[d] / m0[d] - 1.0) for d in both]) * 100
        d2 = np.array([abs(m2[d] / m0[d] - 1.0) for d in both]) * 100
        rows.append(dict(label=inst["label"], n0=len(t0), n1=len(t1),
                         days0=len(m0), days1=len(m1),
                         mv1=float(np.median(d1)), mv2=float(np.median(d2)),
                         o0=len(k0[3]), o1=len(k1[3]), o2=len(k2[3]),
                         g0=len(g0[3]), g1=len(g1[3])))
        r = rows[-1]
        print(f"{r['label'][:28]:28s} {r['n0']:5d}->{r['n1']:5d}  {r['days0']:5d}->{r['days1']:5d}  | "
              f"{r['mv1']:8.2f}% {r['mv2']:8.2f}% | "
              f"{r['g0']:3d}->{r['g1']:3d} op, {r['o0']:3d}->{r['o1']:3d} loc")
    if not rows:
        print("no streams with both series"); return

    mv1 = np.array([r["mv1"] for r in rows])
    mv2 = np.array([r["mv2"] for r in rows])
    print(f"\n  {len(rows)} streams | nights {sum(r['n0'] for r in rows)} -> "
          f"{sum(r['n1'] for r in rows)}, days with an estimate "
          f"{sum(r['days0'] for r in rows)} -> {sum(r['days1'] for r in rows)}")
    print(f"  movement of the best estimate away from v2, on days v2 already covered:")
    print(f"      flat     median {np.median(mv1):5.2f} %   p90 {np.percentile(mv1, 90):5.2f} %")
    print(f"      weighted median {np.median(mv2):5.2f} %   p90 {np.percentile(mv2, 90):5.2f} %"
          f"   ({int(np.sum(mv2 < mv1))}/{len(rows)} streams move less)")
    print(f"  Kalman outliers, OPERATIONAL screen (one IQR over the series): "
          f"v2 {sum(r['g0'] for r in rows)} -> {cfg} {sum(r['g1'] for r in rows)}")
    print(f"  Kalman outliers, local rolling screen:                        "
          f"v2 {sum(r['o0'] for r in rows)} -> {cfg} {sum(r['o1'] for r in rows)} "
          f"(weighting does not change this: the screen runs before the filter)")

    # --- the only independent check: agreement with a co-located CL61 ---------
    print()
    for lab, key in CL61.items():
        ref, new = load(BASE / f"base_eprof_v2_{lab}.json"), load(CAND / f"cand_{cfg}_{lab}.json")
        cl = cl61_series(key)
        if not ref or not new or not cl:
            continue
        t0, v0, _ = triples(ref)
        t1, v1, u1 = triples(new)
        maps = dict(v2=as_map(kalman_best_estimate(t0, v0)),
                    flat=as_map(kalman_best_estimate(t1, v1)),
                    weighted=as_map(kalman_best_estimate(t1, v1, uncertainties=u1)))
        common = sorted(set(cl) & set.intersection(*(set(m) for m in maps.values())))
        if len(common) < 20:
            print(f"  {lab}: only {len(common)} days shared with the CL61")
            continue
        print(f"  {lab} vs co-located CL61 on {len(common)} common days "
              f"(ratio to its OWN median, so only the SHAPE is compared):")
        for name, m in maps.items():
            r = np.array([m[d] / cl[d] for d in common])
            r = r / np.median(r)
            print(f"      {name:9s} scatter {1.4826 * np.median(np.abs(r - 1.0)) * 100:5.2f} %  "
                  f"p90 |dev| {np.percentile(np.abs(r - 1.0), 90) * 100:5.2f} %")


if __name__ == "__main__":
    main()
