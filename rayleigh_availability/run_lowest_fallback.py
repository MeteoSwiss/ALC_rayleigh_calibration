# -*- coding: utf-8 -*-
"""Re-run ONLY the recovered nights with fallback_objective='lowest'.

The registered eprof_v2.2 inherits _select_optimal's default fallback_objective='score': a
noise-tier recovery picks the best-scoring window anywhere in the search range, and on hazy nights
the score (R2-led) prefers the cleaner-looking HIGH windows -- maximising the window lift that,
multiplied by the station's own C-vs-height slope, produces the recovered-night offset. The
'lowest' mode (fit as low as the gates permit; written for exactly this, comment at
molecular_methods.py:658) was never swept.

Eligibility is identical under both objectives, so availability cannot change; only WHERE the
recovered night fits -- and therefore its constant -- can. Downstream QC (-3 method agreement, -9
layer backstop) sees the new window, so a few nights may legitimately change flag; that is part of
the answer.

Out: <DATA>/candidates/cand_N2.5low_<label>.json  (recovered dates only)
"""
from __future__ import annotations
import json
import logging
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel  # noqa: E402
from calibration.config import InstrumentType  # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
L1_ROOT = Path("D:/E-PROFILE_L1_2026")
MAN = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
VALID = (1.0, 0.5)
# every CL61 stream plus the four CHM15k streams that carry the offset analysis
WANT = [i for i in MAN if i["group"] == "CL61"
        or any(s in i["label"] for s in ("PAYERNE", "LINDENBERG", "GOTTFRIEDING", "GUADIANA"))]


def recovered_dates(lab):
    b = json.loads((DATA / "baselines" / f"base_eprof_v2_{lab}.json").read_text())
    c = json.loads((DATA / "candidates" / f"cand_N2.5_{lab}.json").read_text())
    bv = {d for d, v in b.items() if v[0] in VALID and v[1]}
    return sorted(d for d, v in c.items() if v[0] in VALID and v[1] and d not in bv)


def _f(v):
    try:
        f = float(v)
        return f if f == f else None
    except (TypeError, ValueError):
        return None


def run_stream(inst):
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    lab = inst["label"]
    out_f = DATA / "candidates" / f"cand_N2.5low_{lab}.json"
    rec = json.loads(out_f.read_text()) if out_f.exists() else {}
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.molecular_method = "eprof_v2.2"
    o.molecular_params = dict(max_chi2red=2.5, fallback_objective="lowest")
    o.plot_main = o.plot_all = False
    o.folder_output = DATA / "tmp_lowest" / lab          # per-stream dir: no write collisions
    info = InstrumentInfo(site_name=lab, wmo_id=inst["wmo"], identifier=inst["ident"],
                          instrument_type=InstrumentType(inst["type"]),
                          latitude=inst["lat"], longitude=inst["lon"], altitude=inst["alt"])
    for ds in recovered_dates(lab):
        if ds in rec:
            continue
        try:
            r = calibrate_rayleigh(ds, info, o)
        except Exception as exc:
            rec[ds] = [-99, None, None, None, None, f"EXC {type(exc).__name__}"[:80]]
            continue
        if r is None:
            continue
        rec[ds] = [_f(r.flag), _f(r.lidar_constant), _f(r.uncertainty),
                   _f(r.calibration_bottom_height), _f(r.calibration_top_height),
                   str(r.message or "")[:100]]
    out_f.write_text(json.dumps(rec), encoding="utf-8")
    n_ok = sum(1 for v in rec.values() if v[0] in VALID)
    return lab, len(rec), n_ok


def main():
    jobs = [i for i in WANT if recovered_dates(i["label"])]
    print(f"{len(jobs)} streams, recovered nights only, fallback_objective='lowest'", flush=True)
    with ProcessPoolExecutor(max_workers=int(os.environ.get("RA_WORKERS", "10"))) as ex:
        for fut in as_completed([ex.submit(run_stream, i) for i in jobs]):
            lab, n, n_ok = fut.result()
            print(f"  {lab}: {n} nights, {n_ok} valid", flush=True)
    print("LOWEST_DONE", flush=True)


if __name__ == "__main__":
    main()
