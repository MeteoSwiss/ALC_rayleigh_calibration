# -*- coding: utf-8 -*-
"""Is there really no aerosol? — screen the CHM15k with the CO-LOCATED CL61's classification.

The pre-fit mask built from the CHM15k's OWN classification flags almost nothing in 2-4 km on the
nights v2.2 recovers, which is why it does not bring their fit windows back down. That result has
one weak spot: a single-channel ceilometer classifier is threshold-based, so it could simply be
blind to a layer that is nonetheless strong enough to bias a Rayleigh fit.

The CL61 next to it is not blind in the same way -- it measures depolarization, so its Cloudnet
classification separates layers a CHM15k cannot. Two questions, in order:

  1. does the CL61 classification see contamination in 2-4 km on the CHM15k's recovered nights,
     where the CHM15k's own classification sees none?
  2. if the CL61's mask is applied to the CHM15k fit, do the low windows become eligible again?

If the answer is no twice, "residual aerosol in the low windows" is not why those nights need a
high window, and the altitude offset of Phase 2 has an instrumental cause instead.

Run:  python rayleigh_availability/cl61_crossmask.py [config] [--apply]
"""
from __future__ import annotations
import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import indicators as IND                                                        # noqa: E402
import run_classification as RC                                                 # noqa: E402
from calibration.io.classification import (                                     # noqa: E402
    classification_files_for_night, read_classification_curtain)
from calibration.rayleigh.calibration import PREFIT_CONTAM_CODES                # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
BASE, CAND, CLASSIF = DATA / "baselines", DATA / "candidates", DATA / "classification"
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
BANDS = ((2000, 4000), (4000, 6000), (6000, 8000))
SITES = [("PAYERNE", "PAYERNE_CHM15k_A", "PAYERNE_CL61_C"),
         ("LINDENBERG", "LINDENBERG_CHM15k_0", "LINDENBERG_CL61_C")]


def load(p):
    return json.loads(p.read_text()) if p.exists() else None


def curtain(label, wmo, ds):
    f = classification_files_for_night(CLASSIF, label, wmo, ds)
    return read_classification_curtain(f) if len(f) == 2 else None


def night_bands(cur, ds, half_h=4.0):
    """Contaminated fraction per altitude band over the dark hours around solar midnight."""
    t, r, c = cur
    d0 = np.datetime64(f"{ds[:4]}-{ds[4:6]}-{ds[6:]}", "s")
    night = (t >= d0 - np.timedelta64(int(half_h * 3600), "s")) & \
            (t <= d0 + np.timedelta64(int(half_h * 3600), "s"))
    if night.sum() < 50:
        return None
    m = np.isin(c[night], PREFIT_CONTAM_CODES)
    out = []
    for lo, hi in BANDS:
        b = (r >= lo) & (r < hi)
        out.append(float(m[:, b].mean()) if b.any() else np.nan)
    return out


def groups(chm_label, cfg):
    ref = load(BASE / f"base_eprof_v2_{chm_label}.json") or {}
    new = load(CAND / f"cand_{cfg}_{chm_label}.json") or {}
    kept, rec, dead = [], [], []
    for d, v in new.items():
        if IND.is_valid(v[0]):
            (kept if IND.is_valid(ref.get(d, [0])[0]) else rec).append(d)
        elif v[0] == -2:
            dead.append(d)
    return ref, new, dict(kept=sorted(kept), recovered=sorted(rec), still_rejected=sorted(dead))


def main():
    cfg = next((a for a in sys.argv[1:] if not a.startswith("--")), "N2.5")
    do_apply = "--apply" in sys.argv
    print(f"CHM15k screened by its own vs the co-located CL61 classification (config {cfg})\n")

    for site, chm_label, cl61_label in SITES:
        inst = next((i for i in MANIFEST if i["label"] == chm_label), None)
        cl61 = next((e for e in RC.extra_streams() if e["label"] == cl61_label), None)
        if inst is None or cl61 is None or not (CLASSIF / cl61_label).is_dir():
            print(f"{site}: no CL61 classification yet\n")
            continue
        ref, new, grp = groups(chm_label, cfg)
        print(f"== {site} ==")
        print(f"  {'group':16s} {'n':>4s} | {'own 2-4km':>10s} {'CL61 2-4km':>11s} | "
              f"{'own 4-6km':>10s} {'CL61 4-6km':>11s}")
        rows = {}
        for name, dates in grp.items():
            own, ref61 = [], []
            for ds in dates:
                a = curtain(chm_label, inst["wmo"], ds)
                b = curtain(cl61_label, cl61["wmo"], ds)
                if a is None or b is None:
                    continue
                ba, bb = night_bands(a, ds), night_bands(b, ds)
                if ba is None or bb is None:
                    continue
                own.append(ba)
                ref61.append(bb)
            if not own:
                continue
            o, c = np.array(own), np.array(ref61)
            rows[name] = (o, c, dates)
            print(f"  {name:16s} {len(o):4d} | {np.median(o[:, 0]) * 100:9.1f}% "
                  f"{np.median(c[:, 0]) * 100:10.1f}% | {np.median(o[:, 1]) * 100:9.1f}% "
                  f"{np.median(c[:, 1]) * 100:10.1f}%")
        if "recovered" in rows and "kept" in rows:
            r_own, r_61, _ = rows["recovered"]
            k_own, k_61, _ = rows["kept"]
            print(f"  -> in 2-4 km the CL61 sees {np.median(r_61[:, 0]) * 100:.1f}% on recovered "
                  f"nights vs {np.median(k_61[:, 0]) * 100:.1f}% on kept nights "
                  f"(p90 recovered {np.percentile(r_61[:, 0], 90) * 100:.0f}%)")
        if do_apply and "recovered" in rows:
            _apply(site, inst, cl61_label, cl61, rows["recovered"][2], cfg)
        print()


def _apply(site, inst, cl61_label, cl61, dates, cfg):
    """Re-fit the recovered nights with the CL61's mask and report whether the window comes down."""
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    from calibration import (calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel)
    from calibration.config import InstrumentType
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = Path("D:/E-PROFILE_L1_2026")
    o.data_level = DataLevel.L1
    o.molecular_method = "eprof_v2.2"
    o.molecular_params = dict(max_chi2red=2.5)
    o.plot_main = o.plot_all = False
    o.folder_output = DATA / "tmp"
    info = InstrumentInfo(site_name=inst["label"], wmo_id=inst["wmo"], identifier=inst["ident"],
                          instrument_type=InstrumentType(inst["type"]),
                          latitude=inst["lat"], longitude=inst["lon"], altitude=inst["alt"])
    dh, dc, n = [], [], 0
    for ds in dates:
        cur = curtain(cl61_label, cl61["wmo"], ds)
        if cur is None:
            continue
        a = calibrate_rayleigh(ds, info, o)
        b = calibrate_rayleigh(ds, info, o, classification=cur)
        if not (IND.is_valid(a.flag) and IND.is_valid(b.flag)):
            continue
        dh.append(b.calibration_bottom_height - a.calibration_bottom_height)
        dc.append(b.lidar_constant / a.lidar_constant - 1.0)
        n += 1
    if n:
        print(f"  applying the CL61 mask to {n} recovered nights: window bottom "
              f"{np.median(dh):+.0f} m, C_L {np.median(dc) * 100:+.1f}% "
              f"({int(np.sum(np.array(dh) < -100))} nights moved down >100 m)")


if __name__ == "__main__":
    main()
