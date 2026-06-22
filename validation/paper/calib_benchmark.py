"""
calib_benchmark.py — per-channel calibration + Kalman for the paper-validation benchmark, using the
operational Python routines only (no MATLAB). Rayleigh channels (CHM15k/CL61/Mini-MPL) are calibrated
per night from L2_monthly via calibrate_rayleigh; cloud channels (CL31/CL51/CL61) per day from L2_daily
via liquid_cloud_calibration. Each raw nightly/daily series is Kalman-smoothed to a daily grid by the
SAME bridge the MATLAB used (run_kalman_from_matlab.py): for Rayleigh the lidar constant is normalised by
its median before the Kalman (O(1) coefficient) and rescaled after, exactly as load_rayleigh_python_kalman.

Outputs (figs_paper_validation/paper_python/calib/<key>.csv):
    columns time, C_daily, C_daily_std, C_kalman, C_kalman_std   (C_kalman in calibration units)

Usage:  python -m validation.paper.calib_benchmark [station ...]   (default: all benchmark stations)
"""
from __future__ import annotations
import csv
import os
import subprocess
import sys
import tempfile
import warnings
from datetime import datetime, timedelta
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType
from calibration.cloud import liquid_cloud_calibration, CloudCalConfig

REPO = Path(__file__).resolve().parents[2]
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/calib")
OUT.mkdir(parents=True, exist_ok=True)
L2_MONTHLY = Path("A:/E-PROFILE_L2_monthly")
L1_ROOT = Path("D:/E-PROFILE_L1_2026")     # native L1 archive (2015-2026), <wmo>/YYYY/MM/L1_<wmo>_<id><date>.nc
CAMS = "D:/CAMS"
WV_LUT = str(REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc")
KALMAN_PY = "C:/Users/hervo/OneDrive/Documents/Python/improve_alc_calib/run_kalman_from_matlab.py"
RAYLEIGH_METHOD = "eprof_v2"                 # operational molecular-window method (C8 defaults)
CLOUD_AVG_RANGE_M = 30.0                     # bin cloud input to the L2 grid (native L1 over-rejects finer)
# Daily-L2 archives by year (for the per-day cloud calibration).
L2_DAILY = {2026: Path("D:/E-PROFILE_L2_2026"), 2025: Path("D:/E-PROFILE_L2_2025")}
L2_DAILY_FALLBACK = Path("D:/E-PROFILE_L2_2021-2025")
IT = {"CHM15k": InstrumentType.CHM15k, "CL31": InstrumentType.CL31, "CL51": InstrumentType.CL51,
      "CL61": InstrumentType.CL61, "Mini-MPL": InstrumentType.MINI_MPL}

# --- Benchmark channels (mirrors MATLAB make_validation_figures.m) ---------------------------
def _ch(wmo, ident, itype, calib, label, lat, lon, alt, raw=None):
    # raw: folder of flat native files (e.g. Vaisala CL61 "A:/CL61_PAY"); when set, BOTH the
    # Rayleigh (RAW reader, binned to the L2 grid) and cloud (liquid_cloud_calibration_from_data
    # on the day's concatenated profiles) calibrations read these instead of the L2 product.
    return dict(wmo=wmo, ident=ident, itype=itype, calib=calib, label=label, lat=lat, lon=lon, alt=alt, raw=raw)

PAY = (46.8137, 6.9425, 491.0)
CL61_PAY_RAW = "A:/CL61_PAY"   # native Vaisala CL61 (the L2 product has no daily files + bad lat/lon)
BENCHMARK = {
    "payerne": dict(start="20260301", end="20260531", channels=[
        _ch("0-20000-0-06610", "A", "CHM15k", "rayleigh", "CHM15k (Rayleigh)", *PAY),
        _ch("0-20000-0-06610", "B", "CL31", "cloud", "CL31 (cloud)", *PAY),
        _ch("0-20000-0-06610", "C", "CL61", "cloud", "CL61 (cloud)", *PAY, raw=CL61_PAY_RAW),
        # Cloud is calibrated from the native raw (no daily L2 + bad L2 lat/lon); the Rayleigh
        # molecular fit fails on the native CL61 signal (no eligible molecular window) but succeeds
        # on the E-PROFILE L2 product, so the Rayleigh coefficient keeps the L2_monthly source.
        _ch("0-20000-0-06610", "C", "CL61", "rayleigh", "CL61 (Rayleigh)", *PAY),
    ]),
    "amsterdam": dict(start="20260301", end="20260531", channels=[
        _ch("0-20000-0-06240", "A", "CHM15k", "rayleigh", "CHM15k A", 52.317, 4.8037, 6.0),
        _ch("0-20000-0-06240", "B", "CHM15k", "rayleigh", "CHM15k B", 52.317, 4.8037, 6.0),
        _ch("0-20000-0-06240", "C", "CHM15k", "rayleigh", "CHM15k C", 52.317, 4.8037, 6.0),
        _ch("0-20000-0-06240", "D", "CHM15k", "rayleigh", "CHM15k D", 52.317, 4.8037, 6.0),
    ]),
    "uccle": dict(start="20260301", end="20260531", channels=[
        _ch("0-20000-0-06447", "A", "CL51", "cloud", "CL51 (cloud)", 50.8, 4.35, 100.0),
        _ch("0-20000-0-06447", "B", "CL61", "cloud", "CL61 (cloud)", 50.8, 4.35, 100.0),
        _ch("0-20000-0-06447", "B", "CL61", "rayleigh", "CL61 (Rayleigh)", 50.8, 4.35, 100.0),
    ]),
    "sirta": dict(start="20250301", end="20260228", channels=[
        _ch("0-250-1001-07151", "B", "CHM15k", "rayleigh", "CHM15k (Rayleigh)", 48.71, 2.21, 156.0),
        _ch("0-250-1001-07151", "A", "CL31", "cloud", "CL31 (cloud)", 48.71, 2.21, 156.0),
        _ch("0-20000-0-07145", "A", "Mini-MPL", "rayleigh", "Mini-MPL (Rayleigh)", 48.71, 2.21, 156.0),
    ]),
    # EARLINET CHM15k references (for earlinet.py): Leipzig (lei/ari), Cabauw (cbw), Magurele (ino).
    "earlinet": dict(start="20250101", end="20260630", channels=[
        _ch("0-20000-0-10471", "0", "CHM15k", "rayleigh", "Leipzig CHM (lei/ari)", 51.35, 12.43, 125.0),
        _ch("0-20000-0-06348", "A", "CHM15k", "rayleigh", "Cabauw CHM (cbw)", 51.97, 4.93, 1.0),
        _ch("0-20008-0-INO", "B", "CHM15k", "rayleigh", "Magurele CHM (ino)", 44.35, 26.03, 93.0),
    ]),
}


def key_of(ch):
    return f"{ch['wmo']}_{ch['ident']}_{ch['calib']}"


def _daterange(start, end):
    d0 = datetime.strptime(start, "%Y%m%d"); d1 = datetime.strptime(end, "%Y%m%d")
    d = d0
    while d <= d1:
        yield d
        d += timedelta(days=1)


def _daily_file(ch, d, level="L2"):
    """Locate the daily L1 or L2 file for a cloud channel on date d."""
    ds = d.strftime("%Y%m%d")
    if level == "L1":
        f = L1_ROOT / ch["wmo"] / f"{d.year}" / f"{d.month:02d}" / f"L1_{ch['wmo']}_{ch['ident']}{ds}.nc"
        return f if f.exists() else None
    for root in (L2_DAILY.get(d.year), L2_DAILY_FALLBACK):
        if root is None:
            continue
        f = root / ch["wmo"] / f"{d.year}" / f"{d.month:02d}" / f"L2_{ch['wmo']}_{ch['ident']}{ds}.nc"
        if f.exists():
            return f
    return None


def raw_rayleigh(ch, start, end, level="L2"):
    """Per-night lidar constant (calibrate_rayleigh, eprof_v2). Source: native RAW files when
    ch['raw'] is set, else L1 (D:/E-PROFILE_L1_2026, binned to the L2 grid) or L2_monthly.
    Returns (dates, C, Cstd)."""
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.molecular_method = RAYLEIGH_METHOD
    if ch.get("raw"):
        rp = Path(ch["raw"])
        info = InstrumentInfo(site_name=rp.name, wmo_id=rp.name, identifier="",
                              instrument_type=IT[ch["itype"]], latitude=ch["lat"], longitude=ch["lon"], altitude=ch["alt"])
        o.folder_root = rp.parent; o.data_level = DataLevel.RAW
    else:
        info = InstrumentInfo(site_name=ch["wmo"], wmo_id=ch["wmo"], identifier=ch["ident"],
                              instrument_type=IT[ch["itype"]], latitude=ch["lat"], longitude=ch["lon"], altitude=ch["alt"])
        if level == "L1":
            o.folder_root = L1_ROOT; o.data_level = DataLevel.L1       # l1_bin_to_l2_grid -> 30 m x 300 s
        else:
            o.folder_root = L2_MONTHLY; o.data_level = DataLevel.L2_MONTHLY
    o.cams_folder = Path(CAMS); o.abs_cs_lookup_table = Path("")
    o.apply_wv_correction = (ch["itype"] in ("CL31", "CL51", "CL61"))
    o.folder_output = Path(tempfile.mkdtemp()); o.plot_main = o.plot_all = False
    dates, C, Cstd = [], [], []
    for d in _daterange(start, end):
        try:
            r = calibrate_rayleigh(d.strftime("%Y%m%d"), info, o)
        except Exception:
            continue
        if r.flag in (1, 0.5) and np.isfinite(r.lidar_constant) and r.lidar_constant > 0:
            dates.append(d); C.append(float(r.lidar_constant)); Cstd.append(float(r.uncertainty))
    return dates, np.array(C), np.array(Cstd)


def _cloud_one(f, ch):
    """liquid_cloud_calibration on one L1/L2 file (binned to the L2 grid); (cal_median, cal_std) or None."""
    try:
        res = liquid_cloud_calibration(CloudCalConfig(
            nc_file=str(f), instrument=ch["itype"], apply_wv_correction=True,
            cams_folder=CAMS, abs_cs_lookup_table=WV_LUT, average_range_m=CLOUD_AVG_RANGE_M,
            station_latitude=ch["lat"], station_longitude=ch["lon"], aerosol_lidar_ratio=50.0))
    except Exception:
        return None
    if res.n_profiles > 0 and np.isfinite(res.cal_median) and res.cal_median > 0:
        return float(res.cal_median), float(res.cal_std)
    return None


def _raw_cl61_ceilodata(raw_root, date, ch):
    """Concatenate a day of flat native CL61 files into a cloud CeiloData on the E-PROFILE
    Mm^-1 sr^-1 scale (so the returned C matches the L2-based cloud channels)."""
    from calibration.io.data_loader import load_raw_data
    from calibration.cloud.calibration import CeiloData
    files = sorted(Path(raw_root).glob(f"*{date.strftime('%Y%m%d')}*.nc"))
    if not files:
        return None
    d = load_raw_data(files, IT[ch["itype"]])
    if d is None:
        return None
    te = np.asarray(d.time, float)            # epoch seconds
    if te.size == 0:
        return None
    rng = np.asarray(d.range_alc, float)
    beta = np.asarray(d.rcs, float)
    if beta.shape[0] == te.size:
        beta = beta.T                          # -> (range, time)
    beta = beta * 1e6                          # 1/(m sr) -> Mm^-1 sr^-1 (E-PROFILE L2 stored scale)
    cbh = np.asarray(d.cbh, float) if getattr(d, "cbh", None) is not None else np.full(te.size, np.nan)
    if cbh.ndim > 1:
        cbh = cbh[:, 0] if cbh.shape[0] == te.size else cbh[0, :]
    t64 = np.datetime64("1970-01-01T00:00:00", "ns") + (te * 1e9).astype("timedelta64[ns]")
    tnum = te / 86400.0 + 719529.0             # MATLAB datenum (datenum(1970,1,1)=719529)
    return CeiloData(time=t64, time_num=tnum, station_altitude=ch["alt"], station_latitude=ch["lat"],
                     station_longitude=ch["lon"], range=rng, range_resol=float(rng[1] - rng[0]),
                     beta=np.ascontiguousarray(beta), cbh=cbh, quality_flag=None,
                     window_transmission=None, laser_energy=None)


def raw_cloud_from_raw(ch, raw_root, start, end):
    """Per-day cloud coefficient from native CL61 files (liquid_cloud_calibration_from_data),
    binned to the 30 m grid where the O'Connor method succeeds (the native 4.8 m over-rejects)."""
    from calibration.cloud.calibration import liquid_cloud_calibration_from_data
    dates, C, Cstd = [], [], []
    for d in _daterange(start, end):
        data = _raw_cl61_ceilodata(raw_root, d, ch)
        if data is None:
            continue
        cfg = CloudCalConfig(nc_file="", instrument=ch["itype"], apply_wv_correction=True,
                             cams_folder=CAMS, abs_cs_lookup_table=WV_LUT,
                             station_latitude=ch["lat"], station_longitude=ch["lon"],
                             aerosol_lidar_ratio=50.0, average_range_m=30.0, average_time_s=300.0)
        try:
            res = liquid_cloud_calibration_from_data(data, cfg)
        except Exception:
            continue
        if res.n_profiles > 0 and np.isfinite(res.cal_median) and res.cal_median > 0:
            dates.append(d); C.append(float(res.cal_median)); Cstd.append(float(res.cal_std))
    return dates, np.array(C), np.array(Cstd)


def raw_cloud(ch, start, end, level="L2"):
    """Per-day cloud coefficient from L1 or L2 daily files (liquid_cloud_calibration). (dates, C, Cstd).

    Dispatches to the native-file path when ch['raw'] is set (e.g. Payerne CL61); else falls back
    to the monthly L2 file (one value per month) when a channel has NO daily archive.
    NOTE: the L1 cloud coefficient is on the raw rcs_0 scale (V*m^2, factor 1), the L2 one on the
    attenuated-backscatter scale (Mm^-1 sr^-1) - the two are reconciled as a lidar constant downstream.
    """
    if ch.get("raw"):
        return raw_cloud_from_raw(ch, ch["raw"], start, end)
    dates, C, Cstd = [], [], []
    for d in _daterange(start, end):
        f = _daily_file(ch, d, level)
        if f is None:
            continue
        r = _cloud_one(f, ch)
        if r is not None:
            dates.append(d); C.append(r[0]); Cstd.append(r[1])
    if dates or level == "L1":
        return dates, np.array(C), np.array(Cstd)
    # monthly fallback (L2 only): one cloud value per L2_monthly file, dated at mid-month
    seen = set()
    for d in _daterange(start, end):
        ym = (d.year, d.month)
        if ym in seen:
            continue
        seen.add(ym)
        mf = L2_MONTHLY / ch["wmo"] / f"{d.year}" / f"L2_{ch['wmo']}_{ch['ident']}{d.year}{d.month:02d}.nc"
        if not mf.exists():
            continue
        r = _cloud_one(mf, ch)
        if r is not None:
            dates.append(datetime(d.year, d.month, 15)); C.append(r[0]); Cstd.append(r[1])
    return dates, np.array(C), np.array(Cstd)


def kalman(dates, C, Cstd, normalise):
    """Run the operational Kalman bridge. For Rayleigh, normalise by median (O(1)) then rescale."""
    if len(C) < 3:
        return None
    med = float(np.median(C)) if normalise else 1.0
    cin = C / med
    sin = Cstd / med
    with tempfile.TemporaryDirectory() as td:
        fin = Path(td) / "in.csv"; fout = Path(td) / "out.csv"
        with open(fin, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["time", "C", "C_std"])
            for dt, c, s in zip(dates, cin, sin):
                w.writerow([dt.strftime("%Y-%m-%d"), c, s if np.isfinite(s) else ""])
        try:
            subprocess.run([sys.executable, KALMAN_PY, str(fin), str(fout)], check=True,
                           capture_output=True, text=True, timeout=600)
        except Exception as exc:
            print(f"    kalman failed: {exc}"); return None
        rows = list(csv.DictReader(open(fout, encoding="utf-8")))
    out = []
    for r in rows:
        def g(k):
            v = r.get(k, "")
            return float(v) if v not in ("", "nan", "NaN") else float("nan")
        out.append((r["time"][:10], g("C_daily") * med, g("C_daily_std") * med,
                    g("C_kalman") * med, g("C_kalman_std") * med))
    # Fallback (load_rayleigh_python_kalman step 9): if the Kalman output is essentially all NaN
    # (sparse/gappy record the filter can't handle), hold the constant median calibration.
    if out and sum(1 for o in out if np.isfinite(o[3])) < 2:
        fb = float(np.median(C)); fbs = float(np.std(C))
        out = [(o[0], o[1], o[2], fb, fbs) for o in out]
    return out


def run_channel(ch, start, end, level="L2"):
    k = key_of(ch); outp = OUT / f"{k}_{level}.csv"
    # native-raw channels (Payerne CL61) ignore the level (always the same raw source) -> compute once (L2 tag)
    if ch.get("raw") and level == "L1":
        print(f"  {k} [L1]: native-raw channel, see {k}_L2.csv"); return
    if outp.exists() and os.environ.get("FORCE_RECAL", "") != "1":
        print(f"  {k} [{level}]: exists, skip"); return
    if ch["calib"] == "rayleigh":
        dates, C, Cstd = raw_rayleigh(ch, start, end, level); normalise = True
    else:
        dates, C, Cstd = raw_cloud(ch, start, end, level); normalise = False
    print(f"  {k} [{level}] ({ch['label']}): {len(C)} raw {'nights' if ch['calib']=='rayleigh' else 'days'}"
          f"{' median C=%.3e' % np.median(C) if len(C) else ''}", flush=True)
    res = kalman(dates, C, Cstd, normalise)
    if res is None:
        print(f"    -> no Kalman series"); return
    with open(outp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["time", "C_daily", "C_daily_std", "C_kalman", "C_kalman_std"])
        w.writerows(res)
    nk = sum(1 for r in res if np.isfinite(r[3]))
    print(f"    -> {len(res)} daily rows ({nk} finite Kalman) -> {outp.name}", flush=True)


def main():
    warnings.filterwarnings("ignore")
    args = sys.argv[1:]
    levels = [a for a in args if a in ("L1", "L2")] or ["L2", "L1"]
    stations = [a for a in args if a not in ("L1", "L2")] or list(BENCHMARK)
    for level in levels:
        for st in stations:
            if st not in BENCHMARK:
                print(f"unknown station {st}"); continue
            cfg = BENCHMARK[st]
            print(f"== {st} [{level}] ({cfg['start']}-{cfg['end']}) ==", flush=True)
            for ch in cfg["channels"]:
                run_channel(ch, cfg["start"], cfg["end"], level)
    print("CALIB_BENCHMARK_DONE", flush=True)


if __name__ == "__main__":
    main()
