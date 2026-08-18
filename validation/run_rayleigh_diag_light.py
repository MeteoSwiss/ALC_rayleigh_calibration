"""
run_rayleigh_diag_light.py — lightweight per-station instrument-health + FT-aerosol diagnostics for
all CHM15k + Mini-MPL streams. The calibration metrics (constant time series, valid%, sigma_SD,
outlier%, clear-sky fraction) already exist in the network run (net_L1_*.json); this only adds the
cause metrics, by SAMPLING ~12 nights per station (disk-light) instead of re-running every night:

  * signal strength = median near-range (0.2-1 km) raw P            (low -> weak laser / overlap)
  * background noise = std of far-range (13-15 km) raw P            (high -> electronic background)
  * laser_life_time, window_transmission                            (CHM15k HK; NaN for Mini-MPL)
  * FT scattering ratio = median window scattering on sampled fit-nights (high -> FT aerosol)

Output: figs_paper_validation/rayleigh_diag/diaglite_<label>.json
        = {sig, bg, laser, wtrans, scat, n_hk, n_scat}
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
from calibration.rayleigh.molecular_methods import compute_window_grid

REPO = Path(__file__).resolve().parents[1]
MANIFEST = [m for m in json.loads((REPO / "validation" / "scope_network_2026.json").read_text())
            if m["group"] in ("CHM15k", "Mini-MPL")]
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/rayleigh_diag")
NET = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/network_v2_v11")
OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path("D:/E-PROFILE_L1_2026")
HALF = tuple(range(250, 2000, 240))
N_HK = 12          # nights sampled for HK/background
N_SCAT = 10        # fit-nights sampled for scattering


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
    fp = ROOT / wmo / "2026" / ds[4:6] / f"L1_{wmo}_{ident}{ds}.nc"
    if not fp.exists():
        return None
    try:
        with Dataset(fp, "r") as d:
            rng = np.asarray(d.variables["range"][:], float)
            r2 = np.where(rng > 0, rng ** 2, np.nan)
            near = (rng >= 200) & (rng <= 1000)
            far = (rng >= 13000) & (rng <= 15000)
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
            return [sig, bg, med("laser_life_time"), med("window_transmission")]
    except Exception:
        return None


def sample(seq, k):
    seq = list(seq)
    if len(seq) <= k:
        return seq
    idx = np.linspace(0, len(seq) - 1, k).round().astype(int)
    return [seq[i] for i in sorted(set(idx))]


def run_one(inst):
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.CRITICAL)
    label = inst["label"]
    # all candidate dates (for HK sampling) and the fit-nights (for scattering sampling)
    d0 = datetime.strptime(inst["first"], "%Y%m%d"); d1 = datetime.strptime(inst["last"], "%Y%m%d")
    all_days = []
    d = d0
    while d <= d1:
        all_days.append(d.strftime("%Y%m%d")); d += timedelta(days=1)
    netf = NET / f"net_L1_{label}.json"
    fit_nights = sorted(json.loads(netf.read_text()).keys()) if netf.exists() else []

    sig, bg, laser, wt = [], [], [], []
    for ds in sample(all_days, N_HK):
        hk = hk_from_file(inst["wmo"], inst["ident"], ds)
        if hk:
            sig.append(hk[0]); bg.append(hk[1]); laser.append(hk[2]); wt.append(hk[3])

    scat = []
    if fit_nights:
        o = base_options()
        info = InstrumentInfo(site_name=label, wmo_id=inst["wmo"], identifier=inst["ident"],
                              instrument_type=InstrumentType(inst["type"]),
                              latitude=inst["lat"], longitude=inst["lon"], altitude=inst["alt"])
        for ds in sample(fit_nights, N_SCAT):
            fin = {}
            try:
                calibrate_rayleigh(ds, info, o, fit_inputs_out=fin)
                if fin:
                    g = compute_window_grid(fin["signal"], fin["p_mol"], fin["range_alc"], HALF,
                                            range_start_m=2000, range_end_m=6000, increment_bins=8,
                                            signal_stack=fin["signal_stack"])
                    scat.append(float(np.nanmedian(g.scattering_ratio)))
            except Exception:
                continue

    def m(x):
        x = [v for v in x if np.isfinite(v)]
        return float(np.median(x)) if x else float("nan")
    rec = dict(sig=m(sig), bg=m(bg), laser=m(laser), wtrans=m(wt), scat=m(scat),
               n_hk=len(sig), n_scat=len(scat))
    (OUT / f"diaglite_{label}.json").write_text(json.dumps(rec), encoding="utf-8")
    return label, inst["group"], rec["n_hk"], rec["n_scat"]


def main():
    workers = int(os.environ.get("RC_WORKERS", "6"))
    print(f"rayleigh light diag: {len(MANIFEST)} streams, sample {N_HK} HK + {N_SCAT} scat nights, {workers} workers")
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_one, m): m for m in MANIFEST}
        for fut in as_completed(futs):
            try:
                label, grp, nhk, nsc = fut.result(); done += 1
                if done % 25 == 0:
                    print(f"  [{done}/{len(MANIFEST)}] {label} hk={nhk} scat={nsc}")
            except Exception as e:
                print(f"  FAILED: {e}")
    print("DIAGLITE_DONE")


if __name__ == "__main__":
    main()
