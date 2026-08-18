"""
run_rayleigh_diagnostics.py — per-night diagnostic over ALL CHM15k + Mini-MPL L1 streams, to (a)
build the constant time series and (b) attribute WHY a station is problematic. For each night:

  * night fate (from calibrate_rayleigh flag): success(1) / fit_failed(-2,-6,-7,-8 = clear but no
    valid window) / cloudy(-1) / nodata(0) / cams(-4) / nan(-5)
  * the v2 lidar constant + the window's scattering ratio (aerosol-in-FT proxy), when fit-reaching
  * instrument health from a light read of the L1 file (L2 has no HK):
      - signal strength  = median raw P in 0.2-1.0 km   (low -> weak laser / overlap)
      - background noise = std of raw P in 13-15 km      (high -> electronic background)
      - SNR proxy        = signal strength / background noise
      - laser_life_time, status_laser, window_transmission (CHM15k HK; NaN for Mini-MPL)

Aggregating per station then lets us classify the dominant cause: insufficient clear sky
(low clear fraction), FT aerosol (high scattering ratio), low laser (low signal / high laser age),
electronic background (high far-range noise / low SNR).

Output: figs_paper_validation/rayleigh_diag/diag_<label>.json = {ds: [flag, cl, scat, sig, bg, laser, wtrans]}
"""
from __future__ import annotations
import json
import logging
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np
from netCDF4 import Dataset

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType
from calibration.rayleigh.molecular_methods import select_molecular_window, compute_window_grid

REPO = Path(__file__).resolve().parents[1]
MANIFEST = [m for m in json.loads((REPO / "validation" / "scope_network_2026.json").read_text())
            if m["group"] in ("CHM15k", "Mini-MPL")]
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/rayleigh_diag")
OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path("D:/E-PROFILE_L1_2026")
HALF = tuple(range(250, 2000, 240))
QC = 15.0


def base_options():
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = ROOT
    o.data_level = DataLevel.L1
    o.molecular_method = "eprof_v2"
    o.apply_wv_correction = True
    o.plot_main = False
    o.plot_all = False
    o.folder_output = OUT
    return o


def hk_from_file(wmo, ident, ds):
    """Light read of the current-day L1 file: signal strength, far background, laser HK."""
    mo = ds[4:6]
    fp = ROOT / wmo / "2026" / mo / f"L1_{wmo}_{ident}{ds}.nc"
    if not fp.exists():
        return [np.nan] * 4
    try:
        with Dataset(fp, "r") as d:
            rng = np.asarray(d.variables["range"][:], float)
            near = (rng >= 200) & (rng <= 1000)
            far = (rng >= 13000) & (rng <= 15000)
            # raw P = rcs_0 / r^2  (rcs is range-corrected); guard r=0
            r2 = np.where(rng > 0, rng ** 2, np.nan)
            sig = bg = np.nan
            if near.any():
                rc = np.ma.filled(np.ma.masked_invalid(np.asarray(d.variables["rcs_0"][:, near], float)), np.nan)
                sig = float(np.nanmedian(rc / r2[near]))
            if far.any():
                rc = np.ma.filled(np.ma.masked_invalid(np.asarray(d.variables["rcs_0"][:, far], float)), np.nan)
                bg = float(np.nanstd(rc / r2[far]))
            def med(*names):
                for n in names:
                    if n in d.variables:
                        v = np.ma.filled(np.ma.masked_invalid(np.asarray(d.variables[n][:], float)), np.nan)
                        if np.any(np.isfinite(v)):
                            return float(np.nanmedian(v))
                return np.nan
            laser = med("laser_life_time")
            wtrans = med("window_transmission")
        return [sig, bg, laser, wtrans]
    except Exception:
        return [np.nan] * 4


def run_one(args):
    inst = args
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.CRITICAL)
    o = base_options()
    info = InstrumentInfo(site_name=inst["label"], wmo_id=inst["wmo"], identifier=inst["ident"],
                          instrument_type=InstrumentType(inst["type"]),
                          latitude=inst["lat"], longitude=inst["lon"], altitude=inst["alt"])
    d0 = datetime.strptime(inst["first"], "%Y%m%d")
    d1 = datetime.strptime(inst["last"], "%Y%m%d")
    per_night = {}
    d = d0
    while d <= d1:
        ds = d.strftime("%Y%m%d"); d += timedelta(days=1)
        fin = {}
        try:
            res = calibrate_rayleigh(ds, info, o, fit_inputs_out=fin)
        except Exception:
            continue
        flag = float(res.flag)
        cl = scat = np.nan
        if fin:
            try:
                grid = compute_window_grid(fin["signal"], fin["p_mol"], fin["range_alc"], HALF,
                                           range_start_m=2000, range_end_m=6000, increment_bins=8,
                                           signal_stack=fin["signal_stack"])
                scat = float(np.nanmedian(grid.scattering_ratio))
                w = select_molecular_window("eprof_v2", fin["signal"], fin["p_mol"], fin["range_alc"],
                                            HALF, grid=grid)
                if w.ok and np.isfinite(w.rel_error) and w.rel_error <= QC and w.cl > 0:
                    cl = float(w.cl)
            except Exception:
                pass
        sig, bg, laser, wtrans = hk_from_file(inst["wmo"], inst["ident"], ds)
        per_night[ds] = [flag, cl, scat, sig, bg, laser, wtrans]
    (OUT / f"diag_{inst['label']}.json").write_text(json.dumps(per_night), encoding="utf-8")
    n_data = sum(1 for v in per_night.values() if v[0] != 0)
    n_clear = sum(1 for v in per_night.values() if v[0] in (1, -2, -6, -7, -8))
    n_ok = sum(1 for v in per_night.values() if v[0] == 1)
    return inst["label"], inst["group"], len(per_night), n_data, n_clear, n_ok


def main():
    jobs = list(MANIFEST)
    workers = int(os.environ.get("RC_WORKERS", "14"))
    print(f"rayleigh diagnostics: {len(jobs)} streams (CHM15k+Mini-MPL), L1, {workers} workers")
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_one, j): j for j in jobs}
        for fut in as_completed(futs):
            try:
                label, group, n, nd, nc, nok = fut.result()
                done += 1
                if done % 20 == 0:
                    print(f"  [{done}/{len(jobs)}] {label:22s} ({group:8s}) nights={n} data={nd} clear={nc} ok={nok}")
            except Exception as e:
                print(f"  FAILED: {e}")
    print("RAYLEIGH_DIAG_DONE")


if __name__ == "__main__":
    main()
