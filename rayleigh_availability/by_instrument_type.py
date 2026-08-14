# -*- coding: utf-8 -*-
"""What the noise-aware gates do to each INSTRUMENT TYPE, not just the CHM15k.

The study was scoped to the CHM15k because that is where v2 over-rejects, but the corpus also
carries CL61 and Mini-MPL streams and the gates are type-agnostic — so the change has to be checked
where it was not aimed. A Mini-MPL is a photon-counting instrument with a very different noise
character from a CHM15k, which is exactly the case a NOISE-relative gate could get wrong.

Note the classification mask is not applicable to the Mini-MPL at all (ceiloclass has no reader for
it), so for that group only the gate change is in play.

Run:  python rayleigh_availability/by_instrument_type.py [config]
"""
from __future__ import annotations
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
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())


def load(p):
    return json.loads(p.read_text()) if p.exists() else None


def series(rec):
    t, v, u = [], [], []
    for d in sorted(rec):
        f, c, unc = rec[d][0], rec[d][1], rec[d][2]
        if IND.is_valid(f) and c and c > 0:
            t.append(datetime.strptime(d, "%Y%m%d"))
            v.append(float(c))
            u.append(float(unc) if unc and unc > 0 else np.nan)
    return t, np.array(v), np.array(u)


def main():
    cfg = sys.argv[1] if len(sys.argv) > 1 else "N2.5"
    print(f"Impact of eprof_v2.2 / {cfg} by instrument type\n")
    print(f"{'type':10s} {'streams':>7s} {'clear nights':>12s} | {'v1.1':>7s} {'v2':>7s} "
          f"{'v2.2':>7s} | {'gained':>6s} {'moved':>6s} {'lost':>5s} | "
          f"{'outliers op':>12s} {'rel unc %':>16s}")
    for group in ("CHM15k", "CL61", "Mini-MPL"):
        insts = [i for i in MANIFEST if i["group"] == group]
        nclear = nv1 = nv2 = nvn = 0
        gained = moved = lost = 0
        o_ref = o_new = 0
        uk, ur = [], []
        used = 0
        for i in insts:
            b1 = load(BASE / f"base_eprof_v1.1_{i['label']}.json")
            b2 = load(BASE / f"base_eprof_v2_{i['label']}.json")
            cn = load(CAND / f"cand_{cfg}_{i['label']}.json")
            if not b2 or not cn:
                continue
            used += 1
            common = set(b2) & set(cn)
            for d in common:
                if not IND.is_clear(b2[d][0]):
                    continue
                nclear += 1
                a_ok, b_ok = IND.is_valid(b2[d][0]), IND.is_valid(cn[d][0])
                nv2 += a_ok
                nvn += b_ok
                if b1 and d in b1 and IND.is_valid(b1[d][0]):
                    nv1 += 1
                if b_ok and not a_ok:
                    gained += 1
                    if cn[d][2] and cn[d][1]:
                        ur.append(100 * cn[d][2] / cn[d][1])
                elif a_ok and not b_ok:
                    lost += 1
                elif a_ok and b_ok and b2[d][1] and cn[d][1]:
                    if abs(cn[d][1] / b2[d][1] - 1.0) > 0.01:
                        moved += 1
                    if cn[d][2]:
                        uk.append(100 * cn[d][2] / cn[d][1])
            for rec, acc in ((b2, "ref"), (cn, "new")):
                t, v, _ = series(rec)
                if len(t) >= 10:
                    k = kalman_best_estimate(t, v, outlier_mode="global", return_rejected=True)
                    if acc == "ref":
                        o_ref += len(k[3])
                    else:
                        o_new += len(k[3])
        if not nclear:
            print(f"{group:10s} (no outputs)")
            continue
        print(f"{group:10s} {used:7d} {nclear:12d} | {100*nv1/nclear:6.1f}% {100*nv2/nclear:6.1f}% "
              f"{100*nvn/nclear:6.1f}% | {gained:6d} {moved:6d} {lost:5d} | "
              f"{o_ref:5d}->{o_new:5d} | kept {np.median(uk) if uk else np.nan:4.1f} "
              f"new {np.median(ur) if ur else np.nan:4.1f}")
    print("\n  availability is per CLEAR night (the denominator the gates act on); 'moved' counts "
          "nights valid under BOTH where |dC_L/C_L| > 1 %; outliers use the OPERATIONAL screen.")
    print("  The classification mask does not apply to Mini-MPL (ceiloclass has no reader for it).")


if __name__ == "__main__":
    main()
