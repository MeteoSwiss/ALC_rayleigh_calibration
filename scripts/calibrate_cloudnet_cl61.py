"""
calibrate_cloudnet_cl61.py — calibrate the downloaded Cloudnet CL61 raw (Lindenberg +
Hyytiala) over the available period with the 7 Rayleigh molecular-window methods (via the
RAW reader, WV-corrected) and — once the Python cloud port is verified — the liquid-cloud
method. Handles daily and 5-min files (RAW reader concatenates per day-folder).

Rayleigh part is ready now. Cloud part imports cloud_calibration if present.
"""
from __future__ import annotations
import os, sys, glob, json, logging, warnings
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
from rayleigh_calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from rayleigh_calibration.config import InstrumentType
from compare_molecular_methods import METHODS, run_methods, calibrates

ROOT = Path("R:/CL61/RAW_cloudnet_dl")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloudnet_cl61")
OUT.mkdir(parents=True, exist_ok=True)
SAMPLE = [3, 9, 15, 21, 27]
SITES = [
    dict(label="Lindenberg_CL61", wmo="lindenberg", lat=52.21, lon=14.12, alt=123.0),
    dict(label="Hyytiala_CL61",   wmo="hyy",        lat=61.845, lon=24.287, alt=179.0),  # existing folder is 'hyy' (1-min files)
]


def avail_days(wmo):
    """Days with data under <ROOT>/<wmo>/ — either a daily .nc or a non-empty day folder."""
    root = ROOT / wmo
    days = set()
    for p in root.glob("????????.nc"):
        if p.stat().st_size > 0:
            days.add(p.stem)
    for p in root.glob("2*"):
        if p.is_dir() and any(p.glob("*.nc")):
            days.add(p.name)
    return sorted(days)


def base_options():
    o = CalibrationOptions.from_json(Path("options.json"))
    o.folder_root = ROOT
    o.data_level = DataLevel.RAW
    o.molecular_source = "standard"
    o.apply_wv_correction = True          # CL61 = 910 nm -> WV required (CAMS 2025-2026 present)
    o.plot_main = False
    o.plot_all = False
    o.folder_output = OUT
    return o


def main():
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    all_results = {}
    for s in SITES:
        info = InstrumentInfo(site_name=s["label"], wmo_id=s["wmo"], identifier="",
                              instrument_type=InstrumentType.CL61,
                              latitude=s["lat"], longitude=s["lon"], altitude=s["alt"])
        days = set(avail_days(s["wmo"]))
        out_json = OUT / f"rayleigh_{s['label']}.json"
        if not days:
            print(f"{s['label']}: no day-folders under {ROOT/s['wmo']} - skipping (no redownload)")
            continue
        if out_json.exists() and os.environ.get("FORCE_RECAL", "") != "1":
            print(f"{s['label']}: {out_json.name} exists - skipping (set FORCE_RECAL=1 to recompute)")
            continue
        # calibrate sampled nights for which BOTH the night-date and its prev-day have data
        from datetime import date, timedelta
        nights = sorted(
            d for d in days
            if int(d[6:8]) in SAMPLE
            and (date(int(d[:4]), int(d[4:6]), int(d[6:8])) - timedelta(days=1)).strftime("%Y%m%d") in days
        )
        per_night = {}
        n_cal = 0
        for ds in nights:
            o = base_options()
            fin = {}
            try:
                calibrate_rayleigh(ds, info, o, fit_inputs_out=fin)
            except Exception:
                continue
            if not fin:
                continue
            res = run_methods(fin["signal"], fin["p_mol"], fin["range_alc"], fin["signal_stack"])
            per_night[ds] = res
            if any(calibrates(w) for w in res.values()):
                n_cal += 1
        all_results[s["label"]] = per_night
        dump = {ds: {m: [bool(r[m].ok), float(r[m].cl), float(r[m].cl_err), float(r[m].rel_error)]
                     for m in METHODS} for ds, r in per_night.items()}
        out_json.write_text(json.dumps(dump), encoding="utf-8")
        print(f"{s['label']}: {len(days)} day-folders, {len(nights)} complete sampled nights, "
              f"{n_cal} with >=1 calibration -> rayleigh_{s['label']}.json")

    # cloud calibration — run on ALL available days (cloud cal wants maximum data coverage)
    try:
        from rayleigh_calibration.cloud_calibration import liquid_cloud_calibration, CloudCalConfig
        import json as _json
        opt = _json.loads(Path("options.json").read_text())
        for s in SITES:
            out_cloud = OUT / f"cloud_{s['label']}.json"
            if out_cloud.exists() and os.environ.get("FORCE_RECAL", "") != "1":
                print(f"{s['label']} cloud: {out_cloud.name} exists - skipping")
                continue
            all_days = avail_days(s["wmo"])
            cloud_results = {}
            n_cloud = 0
            for ymd in all_days:
                # prefer daily file, else grab first file in the day folder
                daily_nc = ROOT / s["wmo"] / f"{ymd}.nc"
                if not daily_nc.is_file():
                    folder = ROOT / s["wmo"] / ymd
                    fs = sorted(folder.glob("*.nc")) if folder.is_dir() else []
                    if not fs:
                        continue
                    daily_nc = fs[0]   # cloud cal on first file if no daily concat
                cfg = CloudCalConfig(
                    nc_file=str(daily_nc),
                    instrument="CL61",
                    apply_wv_correction=True,
                    cams_folder=opt.get("cams_folder", "D:/CAMS/"),
                    abs_cs_lookup_table=opt.get("abs_cs_lookup_table", ""),
                    station_latitude=s["lat"],
                    station_longitude=s["lon"],
                    aerosol_lidar_ratio=50.0,
                )
                try:
                    res = liquid_cloud_calibration(cfg)
                    coef = float(res.cal_mean) if res.n_profiles > 0 else float("nan")
                    std  = float(res.cal_std)  if res.n_profiles > 0 else float("nan")
                    ok = bool(res.n_profiles > 0 and np.isfinite(coef) and coef > 0)
                    cloud_results[ymd] = [ok, coef, std]
                    if ok:
                        n_cloud += 1
                except Exception:
                    cloud_results[ymd] = [False, float("nan"), float("nan")]
            out_cloud.write_text(_json.dumps(cloud_results), encoding="utf-8")
            print(f"{s['label']} cloud: {len(all_days)} days, {n_cloud} with valid C -> cloud_{s['label']}.json")
    except ImportError:
        print("cloud_calibration module not importable - cloud cal skipped.")
    print("CLOUDNET_CL61_DONE")


if __name__ == "__main__":
    main()
