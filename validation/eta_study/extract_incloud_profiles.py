"""Step 2 input: in-cloud attenuated-backscatter windows around the peak of every ACCEPTED
liquid-cloud calibration profile, before and after the WV and eta (multiple-scattering) corrections.

The per-profile dump (rayleigh_availability/cloud_profile_dump.py) keeps scalars only, so this
re-runs today's operational cloud calibration -- the configuration of
scripts/run_network_calibration.py::_do_cloud as replicated in cloud_profile_dump.process_day:
shared read-once loader (`load_instrument_day` -> `slice_to_date` -> 30 s / 10 m working grid),
WV correction mandatory (CAMS model levels, 0.4 deg), PVC eta tables, aerosol transmission
correction -- and harvests, for each accepted profile i (finite S_consistent):

  working grid (what the calibration integrates; 10 m for CL31/CL51, 9.6 m = 2 x 4.8 m for CL61):
    beta_before = attenuated beta as fed to the calibration (rcs / applied constant)
    t2_wv       = two-way WV transmission used by the calibration
    eta         = eta(range) multiplier of the PVC table
    beta_after  = beta_before / t2_wv * eta   (exactly the array the cloud filters and S use)
  on gates peak-60 m .. peak+300 m, peak = argmax of beta_after above cal_minheight (the same
  peak the cloud filters use).
  native L1 gates (only when the working grid coarsened the RANGE, i.e. the CL61): the same
  window on the 4.8 m gates, time-averaged over the native profiles of the 30-s block, with
  t2_wv interpolated in range and eta evaluated on the native ranges.

Per profile scalars: time, date, cbh (instrument), peak range, S_apparent, S_consistent and the
per-profile coefficient C (after the aerosol-transmission correction, as operational).

Units: Payerne CL61 + CL31 (flagged 'distorted'), Lindenberg CL61 (>= 2025-05-01), Palaiseau CL31,
every other CL61 with cloud days, and a seeded, country-stratified random sample of CL31 / CL51
units inside the Europe CAMS box. Days: the unit's known cloud-calibration days
(rayleigh_availability/cbh_native/cbh_scenes.csv) in 2025-01 .. 2026-06 (--all-days: every L1 day).

Jobs = one unit-month; each writes <OUT>/_parts/<key>/<YYYYMM>.npz (resumable: an existing part is
skipped), then the parts of a unit are merged into <OUT>/<key>.npz (compressed).

Run:  python validation/eta_study/extract_incloud_profiles.py --list            (show the units)
      python validation/eta_study/extract_incloud_profiles.py --test <key> <YYYYMM>
      python validation/eta_study/extract_incloud_profiles.py --workers 24 --go
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
# Operational physics, pinned explicitly (also inherited by the spawned workers): the CL61 line
# 910.55 nm / FWHM 0.188 nm, and none of the R&D switches.
os.environ["ALC_WV_SPECTRUM"] = '{"CL61": [910.55, 0.188]}'
for _v in ("ALC_ETA_TABLES", "ALC_WV_DISABLE", "ALC_DARK_PROFILE"):
    os.environ.pop(_v, None)

import argparse
import csv
import json
import subprocess
import sys
import threading
import time as time_mod
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

# Paths: home-PC defaults; the CSCS array (validation/eta_study/cscs/) overrides them by env.
L1_ROOT = Path(os.environ.get("ETA_L1_ROOT", "D:/E-PROFILE_L1_2026"))
# 0.4 deg only -- never A: or the 1-deg folders (D:/CAMS, D:/CAMS_run_v20, CSCS CAMS_Monthly)
CAMS_FOLDERS = os.environ.get("ETA_CAMS", "D:/CAMS_Monthly_04;D:/CAMS_daily")
_FORBIDDEN_CAMS = {"CAMS", "CAMS_run_v20", "CAMS_Monthly"}       # 1-deg sets (orography bias)
if any(Path(f.strip()).name in _FORBIDDEN_CAMS for f in CAMS_FOLDERS.split(";") if f.strip()):
    raise SystemExit(f"refusing a 1-deg CAMS folder in ETA_CAMS={CAMS_FOLDERS!r}")
WV_LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
CENSUS = REPO / "validation" / "scope_l1_2026_census.json"
SCENES_CSV = REPO / "rayleigh_availability" / "cbh_native" / "cbh_scenes.csv"
OUT = Path(os.environ.get("ETA_OUT", "C:/DATA/Projects/202611_eta_study/incloud"))
PERIOD = ("20250101", "20260630")
BELOW_M, ABOVE_M = 60.0, 300.0
EUROPE_BOX = dict(lat=(27.0, 80.0), lon=(-30.0, 45.0))   # default CAMS area (ALC_CAMS_AREA)
N_CL31, N_CL51 = 20, 10
MIN_CLOUD_DAYS = 30
SEED = 20261102

FIXED = {  # key -> (site label, distorted flag, earliest date)
    "0-20000-0-06610_C": ("Payerne", False, None),
    "0-20000-0-06610_B": ("Payerne", True, None),
    "0-20000-0-10393_C": ("Lindenberg", False, "20250501"),   # ~2x hardware step at 2025-05
    "0-250-1001-07151_A": ("Palaiseau", False, None),
}


# =============================================================================
#  Unit selection
# =============================================================================
def _country(wmo: str) -> str:
    """Country proxy for stratified sampling: WMO block number for 0-20000-0-NNNNN ids, the
    ISO-3166 numeric issuer otherwise."""
    parts = wmo.split("-")
    if parts[1] == "20000" and len(parts) > 3:
        return "wmo" + parts[3][:2]
    return parts[1]


def cloud_days():
    """(key -> (type, sorted cloud days within PERIOD)) from the per-day calibration table."""
    out = defaultdict(lambda: [None, set()])
    with open(SCENES_CSV, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if PERIOD[0] <= r["date"] <= PERIOD[1]:
                out[r["key"]][0] = r["type"]
                out[r["key"]][1].add(r["date"])
    return {k: (t, sorted(d)) for k, (t, d) in out.items()}


def select_units(network=False):
    """The 43-unit study sample (default) or, with network=True, every CL31/CL51/CL61 unit that
    has cloud-calibration days in PERIOD and lies inside the Europe CAMS box (the CSCS array)."""
    census = {f"{s['wmo']}_{s['ident']}": s for s in json.load(open(CENSUS, encoding="utf-8"))}
    cd = cloud_days()
    units = {}

    def add(key, role, label=None, distorted=False, earliest=None):
        s = census[key]
        itype, days = cd[key]
        if earliest:
            days = [d for d in days if d >= earliest]
        units[key] = dict(key=key, itype=itype, site=label or s["site"], lat=float(s["lat"]),
                          lon=float(s["lon"]), alt=float(s["alt"]), distorted=distorted, role=role,
                          n_cloud_days=len(days), days=days)

    for key, (label, distorted, earliest) in FIXED.items():
        add(key, "fixed", label, distorted, earliest)
    rng = np.random.default_rng(SEED)

    def inside(key):
        s = census.get(key)
        return (s is not None and EUROPE_BOX["lat"][0] <= s["lat"] <= EUROPE_BOX["lat"][1]
                and EUROPE_BOX["lon"][0] <= s["lon"] <= EUROPE_BOX["lon"][1])

    if network:
        for key, (itype, days) in sorted(cd.items()):
            if itype in ("CL31", "CL51", "CL61") and key not in units and key in census and inside(key):
                add(key, "network")
        return units
    for key, (itype, days) in sorted(cd.items()):
        if itype == "CL61" and key not in units and inside(key):
            add(key, "all_other_cl61")
    for itype, n in (("CL31", N_CL31), ("CL51", N_CL51)):
        pool = [k for k, (t, d) in sorted(cd.items())
                if t == itype and k not in units and inside(k) and len(d) >= MIN_CLOUD_DAYS]
        by_c = defaultdict(list)
        for k in pool:
            by_c[_country(k.rsplit("_", 1)[0])].append(k)
        for c in by_c:
            rng.shuffle(by_c[c])
        countries = sorted(by_c)
        rng.shuffle(countries)
        picked, sites = [], {k.rsplit("_", 1)[0] for k in units}
        while len(picked) < n and any(by_c.values()):      # round-robin over countries
            for c in countries:
                while by_c[c] and len(picked) < n:
                    k = by_c[c].pop()
                    if k.rsplit("_", 1)[0] not in sites:   # one unit per site (no co-located twins)
                        picked.append(k)
                        sites.add(k.rsplit("_", 1)[0])
                        break
        for k in picked:
            add(k, f"random_{itype}")
    return units


# =============================================================================
#  Worker: one unit-month
# =============================================================================
def _window(offs_n, n_gates, j):
    idx = j + offs_n
    ok = (idx >= 0) & (idx < n_gates)
    return idx, ok


def process_day(u, ds):
    """Run the operational cloud calibration on one unit-day; return the per-profile harvest."""
    from calibration.cloud import CloudCalConfig
    from calibration.cloud.calibration import (
        liquid_cloud_calibration_from_data, set_defaults, build_cloud_input_from_day)
    from calibration.cloud._filters import apply_multiple_scattering_correction
    from calibration.io.instrument_day import load_instrument_day

    wmo, ident = u["key"].rsplit("_", 1)
    fp = L1_ROOT / wmo / ds[:4] / ds[4:6] / f"L1_{wmo}_{ident}{ds}.nc"
    if not fp.exists():
        return None, "no L1 file"
    itype = u["itype"]
    cfg = set_defaults(CloudCalConfig(
        nc_file=str(fp), instrument=itype, apply_wv_correction=True,
        apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
        cams_folder=CAMS_FOLDERS, abs_cs_lookup_table=str(WV_LUT), cams_folder_fallback="",
        wv_source="cams", era5_cache="",
        station_latitude=u["lat"], station_longitude=u["lon"],
        average_time_s=0.0, average_range_m=0.0))
    d = datetime.strptime(ds, "%Y%m%d").date()
    idd = load_instrument_day([str(fp)], itype, CAMS_FOLDERS, cams_folder_fallback="",
                              read_cams=False, build_working=False)
    if idd is None:
        return None, "L1 unreadable"
    day = idd.slice_to_date(d)
    data = build_cloud_input_from_day(day, cfg)
    if data.beta is None or not np.any(np.isfinite(data.beta)):
        return None, "no usable signal"
    beta0 = np.array(data.beta, dtype="float64", copy=True)            # (range, time), before
    rng = np.asarray(data.range, dtype="float64")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = liquid_cloud_calibration_from_data(data, cfg)
    s_cons = np.asarray(res.S_consistent, dtype="float64") if res.S_consistent is not None else None
    if s_cons is None or int(getattr(res, "n_profiles", 0)) == 0:
        return None, "no valid cloud profile"
    vidx = np.where(np.isfinite(s_cons))[0]
    t2 = np.asarray(res.trans2_wv, dtype="float64")
    eta = apply_multiple_scattering_correction(np.ones((rng.size, 1)), rng, cfg)[:, 0]
    after = beta0[:, vidx] / t2[:, vidx] * eta[:, None]
    lower_gate = int(np.argmax(rng >= cfg.cal_minheight))
    roi = after.copy()
    roi[:lower_gate + 1] = 0.0
    roi = np.where(np.isfinite(roi), roi, -np.inf)
    jp = np.argmax(roi, axis=0)

    dr = float(rng[1] - rng[0])
    offs = np.arange(-int(np.ceil(BELOW_M / dr - 1e-6)), int(np.ceil(ABOVE_M / dr - 1e-6)) + 1)
    n, K = vidx.size, offs.size
    w = {k: np.full((n, K), np.nan, dtype="float32") for k in ("beta_before", "beta_after", "t2_wv", "eta")}
    for m, (i, j) in enumerate(zip(vidx, jp)):
        idx, ok = _window(offs, rng.size, j)
        w["beta_before"][m, ok] = beta0[idx[ok], i]
        w["beta_after"][m, ok] = after[idx[ok], m]
        w["t2_wv"][m, ok] = t2[idx[ok], i]
        w["eta"][m, ok] = eta[idx[ok]]

    out = dict(time_ns=np.asarray(res.time, dtype="datetime64[ns]").astype("int64"),
               date=np.full(n, int(ds), dtype=np.int32),
               cbh_m=np.asarray(res.cbh, dtype="float32"),
               peak_range_m=rng[jp].astype("float32"),
               s_apparent=np.asarray(res.S_apparent, dtype="float64")[vidx].astype("float32"),
               s_consistent=s_cons[vidx].astype("float32"),
               c_profile=np.asarray(res.all_coefficients, dtype="float64")[vidx].astype("float32"),
               **w)
    out["_dr"] = dr
    out["_offs"] = offs

    # --- native L1 gates when the working grid coarsened the range (CL61: 2 x 4.8 m) ---
    tf, rf = _factors(day)
    if rf > 1:
        nat = day.native
        nrng = np.asarray(nat.range_alc, dtype="float64")
        nrcs = np.ma.filled(np.ma.masked_invalid(np.asarray(nat.rcs, dtype="float64")), np.nan)  # (time, range)
        wrcs = np.ma.filled(np.ma.masked_invalid(np.asarray(day.working.rcs, dtype="float64")), np.nan)
        with np.errstate(all="ignore"):
            ratio = beta0 / wrcs.T
        scale = float(np.nanmedian(ratio[np.isfinite(ratio) & (ratio != 0)]))
        ndr = float(nrng[1] - nrng[0])
        noffs = np.arange(-int(np.ceil(BELOW_M / ndr - 1e-6)), int(np.ceil(ABOVE_M / ndr - 1e-6)) + 1)
        neta = apply_multiple_scattering_correction(np.ones((nrng.size, 1)), nrng, cfg)[:, 0]
        nb = np.full((n, noffs.size), np.nan, dtype="float32")
        na = np.full((n, noffs.size), np.nan, dtype="float32")
        npk = np.full(n, np.nan, dtype="float32")
        for m, (i, j) in enumerate(zip(vidx, jp)):
            prof = np.nanmean(nrcs[i * tf:(i + 1) * tf], axis=0) * scale
            t2n = np.interp(nrng, rng, t2[:, i])
            prof_a = prof / t2n * neta
            lo, hi = max(0, (j - 1) * rf), min(nrng.size, (j + 2) * rf)   # +-1 working gate
            seg = prof_a[lo:hi]
            if not np.isfinite(seg).any():
                continue
            jn = lo + int(np.nanargmax(seg))
            idx, ok = _window(noffs, nrng.size, jn)
            nb[m, ok] = prof[idx[ok]]
            na[m, ok] = prof_a[idx[ok]]
            npk[m] = nrng[jn]
        out.update(native_beta_before=nb, native_beta_after=na, native_peak_range_m=npk)
        out["_ndr"] = ndr
        out["_noffs"] = noffs
    out["_factors"] = (tf, rf)
    return out, f"OK {n}"


def _factors(day):
    from calibration.io.instrument_day import _applied_factors
    return _applied_factors(day.native, day.target_range_m, day.target_time_s)


def process_unit_month(job):
    """job = (unit dict, YYYYMM, day list) -> writes the part file; returns a log line."""
    u, ym, days = job
    part = OUT / "_parts" / u["key"] / f"{ym}.npz"
    if part.exists():
        return u["key"], ym, "exists", 0, 0.0
    t0 = time_mod.time()
    harvest, msgs = [], defaultdict(int)
    for ds in days:
        try:
            out, msg = process_day(u, ds)
        except Exception as exc:  # noqa: BLE001 - one bad day must not kill the month
            out, msg = None, f"{type(exc).__name__}: {str(exc)[:80]}"
        msgs[msg.split()[0] if msg.startswith("OK") else msg] += 1
        if out is not None:
            harvest.append(out)
    part.parent.mkdir(parents=True, exist_ok=True)
    if harvest:
        keys = [k for k in harvest[0] if not k.startswith("_")]
        arrays = {k: np.concatenate([h[k] for h in harvest]) for k in keys}
        arrays["range_rel_m"] = (harvest[0]["_offs"] * harvest[0]["_dr"]).astype("float32")
        arrays["dr_m"] = np.float32(harvest[0]["_dr"])
        if "_ndr" in harvest[0]:
            arrays["native_range_rel_m"] = (harvest[0]["_noffs"] * harvest[0]["_ndr"]).astype("float32")
            arrays["native_dr_m"] = np.float32(harvest[0]["_ndr"])
        arrays["coarsen_factors"] = np.array(harvest[0]["_factors"], dtype=np.int16)
        np.savez_compressed(part, **arrays)
        n = int(arrays["date"].size)
    else:
        np.savez_compressed(part, empty=np.array(True))
        n = 0
    return u["key"], ym, dict(msgs), n, time_mod.time() - t0


def merge_unit(u):
    parts = sorted((OUT / "_parts" / u["key"]).glob("*.npz"))
    zs = [np.load(p) for p in parts]
    zs = [z for z in zs if "empty" not in z.files]
    if not zs:
        return 0
    per_profile = [k for k in zs[0].files if zs[0][k].ndim >= 1 and k not in
                   ("range_rel_m", "native_range_rel_m", "coarsen_factors")]
    arrays = {k: np.concatenate([z[k] for z in zs]) for k in per_profile}
    for k in ("range_rel_m", "native_range_rel_m", "dr_m", "native_dr_m"):
        if k in zs[0].files:
            arrays[k] = zs[0][k]
    arrays["coarsen_factors"] = np.array([z["coarsen_factors"] for z in zs])
    meta = {k: v for k, v in u.items() if k != "days"}
    meta.update(period=PERIOD, cams=CAMS_FOLDERS, wv_spectrum=os.environ["ALC_WV_SPECTRUM"],
                window=f"peak-{BELOW_M:.0f} m .. peak+{ABOVE_M:.0f} m", git=_git_head(),
                created=datetime.now().isoformat(timespec="seconds"),
                note="beta_after = beta_before / t2_wv * eta (working grid); native_* only when the "
                     "working grid coarsened the range (CL61)")
    arrays["meta"] = np.array(json.dumps(meta))
    np.savez_compressed(OUT / f"{u['key']}.npz", **arrays)
    return int(arrays["date"].size)


def _git_head():
    try:
        return subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, timeout=10).stdout.strip() + " (+ working tree)"
    except Exception:  # noqa: BLE001
        return "unknown"


def jobs_for(units, all_days=False):
    jobs = []
    for u in units.values():
        days = u["days"]
        if all_days:
            wmo, ident = u["key"].rsplit("_", 1)
            days = sorted(p.name[-11:-3] for p in L1_ROOT.joinpath(wmo).glob(f"*/*/L1_{wmo}_{ident}*.nc")
                          if PERIOD[0] <= p.name[-11:-3] <= PERIOD[1])
            earliest = FIXED.get(u["key"], (None, None, None))[2]
            if earliest:
                days = [d for d in days if d >= earliest]
        by_m = defaultdict(list)
        for d in days:
            by_m[d[:6]].append(d)
        jobs += [(u, ym, dd) for ym, dd in sorted(by_m.items())]
    return jobs


def check_spectrum():
    from calibration.water_vapor_correction.water_vapor import laser_spectrum_for
    got = tuple(laser_spectrum_for("CL61", 910.0))
    if not np.allclose(got, (910.55, 0.188)):
        raise RuntimeError(f"CL61 spectrum seen by the code is {got}, expected (910.55, 0.188)")
    return got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--test", nargs=2, metavar=("KEY", "YYYYMM"))
    ap.add_argument("--all-days", action="store_true")
    ap.add_argument("--go", action="store_true", help="start the full extraction")
    ap.add_argument("--network", action="store_true",
                    help="every CL31/CL51/CL61 unit with cloud days inside the Europe CAMS box")
    ap.add_argument("--keylist", help="with --network --list: write the unit keys (one per line) here")
    ap.add_argument("--unit", help="CSCS array task: all months of ONE unit, local pool, then merge")
    ap.add_argument("--months", help="with --unit: restrict to these YYYYMM (comma list; dry runs)")
    a = ap.parse_args()
    units = select_units(network=a.network or bool(a.unit))
    print(f"CL61 spectrum seen by the code: {check_spectrum()}", flush=True)
    if a.list:
        for u in units.values():
            print(f"{u['key']:24s} {u['itype']:5s} {u['role']:15s} {u['site'][:22]:22s} "
                  f"{u['lat']:6.2f} {u['lon']:7.2f} cloud_days={u['n_cloud_days']}"
                  + ("  [distorted unit]" if u["distorted"] else ""), flush=True)
        js = jobs_for(units, a.all_days)
        print(f"{len(units)} units, {len(js)} unit-months, {sum(len(j[2]) for j in js)} unit-days", flush=True)
        if a.keylist:
            # largest first, so the array's first (throttled) tasks are the long ones
            keys = sorted(units, key=lambda k: -units[k]["n_cloud_days"] * (2 if units[k]["itype"] == "CL61" else 1))
            # LF endings whatever the OS: the keylist is read by bash on Linux (sed), where a CRLF
            # file yields keys ending in '\r'
            with open(a.keylist, "w", encoding="utf-8", newline="\n") as f:
                f.write("\n".join(keys) + "\n")
            print(f"keylist -> {a.keylist}", flush=True)
        return
    OUT.mkdir(parents=True, exist_ok=True)
    if a.unit:
        u = units[a.unit]
        (OUT / "units").mkdir(parents=True, exist_ok=True)
        with open(OUT / "units" / f"{a.unit}.json", "w", encoding="utf-8") as f:
            json.dump({kk: vv for kk, vv in u.items() if kk != "days"}, f, indent=1)
        js = jobs_for({a.unit: u}, a.all_days)
        if a.months:
            keep = set(a.months.split(","))
            js = [j for j in js if j[1] in keep]
        t0 = time_mod.time()
        print(f"{a.unit}: {len(js)} unit-months, {sum(len(j[2]) for j in js)} days -> {a.workers} processes",
              flush=True)
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for key, ym, msgs, n, dt in ex.map(process_unit_month, js, chunksize=1):
                print(f"  {key} {ym}: {n} profiles ({dt:.0f} s) {msgs}", flush=True)
        print(f"  merged {a.unit}: {merge_unit(u)} profiles in {time_mod.time() - t0:.0f} s", flush=True)
        return
    if a.test:
        key, ym = a.test
        js = [j for j in jobs_for({key: units[key]}, a.all_days) if j[1] == ym]
        for j in js:
            print(process_unit_month(j), flush=True)
        return
    if not a.go:
        raise SystemExit("full run needs --go (run --test first to size it)")
    js = jobs_for(units, a.all_days)
    # largest first (by L1 file count = days in the job), CL61 days ~2x the cost of CL31/CL51
    js.sort(key=lambda j: -len(j[2]) * (2 if j[0]["itype"] == "CL61" else 1))
    with open(OUT / "units.json", "w", encoding="utf-8") as f:
        json.dump({k: {kk: vv for kk, vv in u.items() if kk != "days"} for k, u in units.items()}, f, indent=1)
    print(f"{len(units)} units, {len(js)} unit-months -> {a.workers} processes", flush=True)
    t0 = time_mod.time()
    remaining = defaultdict(int)
    for u, _, _ in js:
        remaining[u["key"]] += 1
    state = dict(done=0, profiles=0)
    stop = threading.Event()

    def progress():
        """Every minute: jobs done / total, profiles so far, GB written -> progress.log."""
        while not stop.wait(60.0):
            gb = sum(p.stat().st_size for p in OUT.rglob("*.npz")) / 1e9
            line = (f"{datetime.now():%Y-%m-%d %H:%M:%S}  jobs {state['done']}/{len(js)}  "
                    f"profiles {state['profiles']}  written {gb:.3f} GB  elapsed {time_mod.time() - t0:.0f} s")
            with open(OUT / "progress.log", "a", encoding="utf-8") as fl:
                fl.write(line + "\n")

    th = threading.Thread(target=progress, daemon=True)
    th.start()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(process_unit_month, j) for j in js]
        for f in as_completed(futs):
            key, ym, msgs, n, dt = f.result()
            state["done"] += 1
            state["profiles"] += n
            remaining[key] -= 1
            print(f"[{state['done']}/{len(js)} {time_mod.time() - t0:7.0f} s] {key} {ym}: {n} profiles "
                  f"({dt:.0f} s) {msgs}", flush=True)
            if remaining[key] == 0:
                print(f"  merged {key}: {merge_unit(units[key])} profiles", flush=True)
    stop.set()
    with open(OUT / "progress.log", "a", encoding="utf-8") as fl:
        fl.write(f"{datetime.now():%Y-%m-%d %H:%M:%S}  DONE jobs {state['done']}/{len(js)}  "
                 f"profiles {state['profiles']}  total {time_mod.time() - t0:.0f} s\n")
    print(f"total {time_mod.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
