"""
run_eprof_v10.py — generate the E-PROF v1.0 baseline (the pre-a4e7140 sign-error
calibration) for the long-run method comparison.

v1.0 = the same 'main' molecular window as E-PROF v1.1, but with the historical bug
re-enabled (options.sign_error_v10=True -> Klett sign error + total-OD reference). It is
run on the SAME fit-nights already stored in each results_<label>.json (the 10 CHM15k +
4 Mini-MPL suitable instruments; the 910 nm instruments are excluded), and injected back
into those JSONs under the method key 'eprof_v10' so the comparison/precision tooling
picks it up alongside the other versions.
"""
from __future__ import annotations
import os, sys, json, logging, warnings
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
from rayleigh_calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from longrun_methods import build_instruments, OUT, ROOT

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.ERROR)

TMP = OUT / "_v10_tmp"
TMP.mkdir(parents=True, exist_ok=True)


def opt():
    o = CalibrationOptions.from_json(Path("options.json"))
    o.folder_root = ROOT
    o.data_level = DataLevel.L2_MONTHLY
    o.molecular_source = "standard"
    o.molecular_method = "main"      # v1.0 uses the legacy 'main' window (= v1.1's window)
    o.sign_error_v10 = True          # ... but with the pre-a4e7140 sign-error math
    o.plot_main = False
    o.plot_all = False
    o.folder_output = TMP
    return o


def to_list(r):
    """Map a CalibrationResult to the 10-field list used in results_*.json."""
    ok = r.flag in (1, 1.0, 0.5)
    cl = float(r.lidar_constant)
    err = float(r.uncertainty)
    rel = abs(100.0 * err / cl) if cl > 0 else 999.0
    bot = float(getattr(r, "bottom_height", 0.0) or 0.0)
    top = float(getattr(r, "top_height", 0.0) or 0.0)
    nan = float("nan")
    return [bool(ok), cl, err, rel, nan, nan, nan, bot, top, (bot + top) / 2.0]


FAIL = [False, -1.0, -1.0, 999.0, float("nan"), float("nan"), float("nan"), 0.0, 0.0, 0.0]


def main():
    insts = build_instruments()
    o = opt()
    total = 0
    for inst in insts:
        jf = OUT / f"results_{inst['label']}.json"
        if not jf.exists():
            print(f"  {inst['label']}: no results JSON, skipping", flush=True)
            continue
        data = json.load(open(jf))
        info = InstrumentInfo(site_name=inst["label"], wmo_id=inst["wmo"], identifier=inst["ident"],
                              instrument_type=inst["itype"], latitude=inst["lat"],
                              longitude=inst["lon"], altitude=inst["alt"])
        ncal = 0
        for ds in sorted(data):
            try:
                r = calibrate_rayleigh(ds, info, o)
                row = to_list(r)
            except Exception:
                row = list(FAIL)
            data[ds]["eprof_v10"] = row
            if row[0]:
                ncal += 1
            total += 1
        jf.write_text(json.dumps(data), encoding="utf-8")
        print(f"  {inst['label']}: {len(data)} fit-nights, {ncal} v1.0-calibrated -> eprof_v10 injected",
              flush=True)
    print(f"EPROF_V10_DONE: {total} night-calibrations across {len(insts)} instruments")


if __name__ == "__main__":
    main()
