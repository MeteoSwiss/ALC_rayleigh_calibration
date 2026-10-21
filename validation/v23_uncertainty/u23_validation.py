# -*- coding: utf-8 -*-
"""v2.3 uncertainty budget -- validation on real data (doc/v2.3/uncertainty.md).

Harness
  The q16a harness (validation/chm15k_stability/q16a_recalibration.py), which reproduces the
  published nightly constants bit-identically: calibrate_rayleigh unmodified, options.json +
  molecular_method eprof_v2.2, the NetCDF write patched to a no-op. Here it runs with
  ``uncertainty_model = "v2.3"`` and a RECORDING configuration that is wider than the defaults
  (sub-windows from 1 km, 1.5 km under the window, up to a per-gate SNR of 1, no QC limit), so that
  every design choice of the vertical term can be re-evaluated offline from the stored sub-window
  constants without recalibrating. C_L itself is the operational one (same ensemble median).

Sub-commands
  run      recalibrate the unit list, 2025-01..2026-07 (every 2nd night; every night for the twins)
  verify   rerun a subset with the shipped defaults and check the offline recomputation == code
  payerne  the Payerne 2026-09-07 night (window 6.06-6.54 km ASL)
  analyse  pull distributions, twins, per-site ranking, QC limit, figures, tables
  kalman   effect of weighting the Kalman best estimate by the new uncertainty
  cost     runtime of the budget (shipped defaults)

Outputs only under C:/DATA/Projects/202610_v23/uncertainty/.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
STAB = REPO / "validation" / "chm15k_stability"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(STAB))

# The operational v2.2 configuration (ops/config.sh, v2.2 cutover) runs the CL61 water-vapour
# correction with the manufacturer line, not the code default 910.74/1.0 nm (which over-corrects
# WV x9 and creates a spurious height gradient). Set it here so the recalibration reproduces the
# published CL61 constants; inherited by the worker processes. Only the CL61 reads it.
os.environ.setdefault("ALC_WV_SPECTRUM", '{"CL61": [910.55, 0.188]}')

OUT = Path("C:/DATA/Projects/202610_v23/uncertainty")
RUN_DIR = OUT / "runs"
L1_ROOT = Path("D:/E-PROFILE_L1_2026")
CAMS = "D:/CAMS_Monthly_04;D:/CAMS_daily"
CENSUS = REPO / "validation" / "scope_l1_2026_census.json"
EPOCH = datetime(1970, 1, 1)
SUCCESS = (1.0, 0.5)

# Recording configuration: wider than the shipped defaults so the choices can be tuned offline.
RECORD_PARAMS = dict(extend_below_m=1500.0, min_altitude_m=1000.0, extend_above_m=1e9,
                     snr_min=1.0, max_rel_uncertainty=1e9)

# Units. CHM15k (30): the 20 q16a units (Schiphol A-D, Budapest A/B, Magurele A/B twins,
# Lindenberg, Granada, Evora, Potenza, Sofia, Aosta, Hohenpeissenberg, Palaiseau, Kleine Scheidegg,
# Wattisham, Schleswig, Oslo) + Payerne A (abnormal receiver), Schleiz, Messina + geographic spread
# (Aldergrove, Zagreb, Torshavn, Vasarosnameny, Gottfrieding, Madrid, Montsec).
# CL61 (3): Lindenberg, Payerne, Aosta. Mini-MPL (3): Lille, Brest, Aleria.
UNITS = [
    ("0-20000-0-06240", "A"), ("0-20000-0-06240", "B"), ("0-20000-0-06240", "C"),
    ("0-20000-0-06240", "D"), ("0-20000-0-12843", "A"), ("0-20000-0-12843", "B"),
    ("0-20008-0-INO", "A"), ("0-20008-0-INO", "B"),
    ("0-20000-0-10393", "0"), ("0-20008-0-UGR", "A"), ("0-620-3704-5480", "B"),
    ("0-20008-0-POT", "A"), ("0-100-20000-0000001", "A"), ("0-380-5-1", "0"),
    ("0-20000-0-10962", "0"), ("0-250-1001-07151", "B"), ("0-20000-0-06736", "A"),
    ("0-20000-0-03590", "A"), ("0-20000-0-10035", "0"), ("0-20000-0-01492", "A"),
    ("0-20000-0-06610", "A"), ("0-20000-0-10564", "0"), ("0-20000-0-00203", "A"),
    ("0-20000-0-03917", "A"), ("0-20000-0-14240", "A"), ("0-208-0-06011", "A"),
    ("0-348-1-72618", "A"), ("0-20000-0-10872", "0"), ("0-724-60000-0000000CIEMAT001", "A"),
    ("0-20008-0-MSA", "A"),
    # CL61
    ("0-20000-0-10393", "C"), ("0-20000-0-06610", "C"), ("0-380-5-1", "B"),
    # Mini-MPL
    ("0-20000-0-07014", "A"), ("0-20000-0-07110", "A"), ("0-20000-0-07774", "A"),
]
# Same-site twins (same wavelength, same night): the atmosphere and the common-mode gain cancel.
TWINS = {"0-20000-0-06240": ("A", "B", "C", "D"), "0-20000-0-12843": ("A", "B"),
         "0-20008-0-INO": ("A", "B")}
# Night sampling, chosen for the shared L1 disk: blocks of 7 CONSECUTIVE nights (a night reads the
# files of D-1 and D, so night D+1 re-uses file D from the OS cache: 8 files for 7 nights), one block
# every 2 weeks for the twins/key units (50 % of the nights) and every 3 weeks for the others (33 %).
# Every unit uses the same calendar blocks, so the twins are sampled on the same nights. The 31-day
# running median of the pull test is taken from the PUBLISHED (all-night) series, which this harness
# reproduces bit-identically, so the sparser sample does not degrade it.
KEY_UNITS = ({("0-20000-0-06240", i) for i in "ABCD"}
             | {("0-20000-0-12843", "A"), ("0-20000-0-12843", "B"), ("0-20008-0-INO", "A"),
                ("0-20008-0-INO", "B"), ("0-20000-0-06610", "A"), ("0-20000-0-10564", "0"),
                ("0-20000-0-00203", "A"), ("0-20000-0-10393", "0"), ("0-20000-0-07014", "A"),
                ("0-20000-0-07110", "A"), ("0-20000-0-07774", "A")})

FIELDS = ["wmo", "ident", "type", "date", "flag", "c", "u_ens", "u_rec", "sig_ens", "sig_vert_rec",
          "sig_mol", "sig_noise", "g_rec", "g_sigma_rec", "chi2_rec", "resid_rec", "mol_bias",
          "n_sub", "z_lo", "z_hi", "zb_agl", "zt_agl", "n_native", "t_mean", "tint_mean",
          "snr_2_3", "snr_3_4", "snr_4_5", "snr_5_6", "snr_6_8", "cpu_s", "cpu_budget_s",
          "sub", "message"]


def census():
    return {(s["wmo"], s["ident"]): s for s in json.load(open(CENSUS, encoding="utf-8"))}


def _info(s):
    from calibration import InstrumentInfo
    from calibration.config import InstrumentType
    itype = {"CHM15k": InstrumentType.CHM15k, "CL61": InstrumentType.CL61,
             "Mini-MPL": InstrumentType.MINI_MPL}[s["type"]]
    return InstrumentInfo(site_name=s.get("site", s["wmo"]), wmo_id=s["wmo"], identifier=s["ident"],
                          instrument_type=itype, latitude=s["lat"], longitude=s["lon"],
                          altitude=s.get("alt", 0.0))


def _options(params=None, model="v2.3"):
    from calibration import CalibrationOptions, DataLevel
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.folder_output = OUT / "nowrite"          # never written (the write is patched out)
    o.cams_folder = CAMS
    o.plot_all = False
    o.plot_main = False
    o.molecular_method = "eprof_v2.2"          # calibration_version 220 of the published files
    o.uncertainty_model = model
    o.uncertainty_params = dict(RECORD_PARAMS if params is None else params)
    return o


_TIMER = {"t": 0.0}


def _install():
    """q16a patches (no NetCDF write, profile-time capture) + a timer around the budget."""
    import q16a_recalibration as Q
    Q._install_patches()
    import calibration.rayleigh.uncertainty as U
    if not getattr(U, "_timed", False):
        orig = U.budget_for_night

        def timed(**kw):
            # perf_counter: process_time has a 15.6 ms granularity on Windows, coarser than the call
            t0 = time.perf_counter()
            try:
                return orig(**kw)
            finally:
                _TIMER["t"] += time.perf_counter() - t0
        U.budget_for_night = timed
        U._timed = True
    return Q


def _enc_sub(b, c):
    if b is None or b.sub_cl.size == 0:
        return ""
    return "|".join([",".join(f"{z:.0f}" for z in b.sub_centre_m),
                     ",".join(f"{x / c:.5f}" for x in b.sub_cl),
                     ",".join(f"{s:.5f}" for s in b.sub_sigma_rel)])


def calibrate_night(s, date_str, hk, params=None, model="v2.3"):
    """One night -> one row dict (or None when there is no L1)."""
    Q = _install()
    from calibration import calibrate_rayleigh
    from calibration.io.data_loader import build_file_paths, load_data
    info = _info(s)
    o = _options(params, model)
    files = [f for f in build_file_paths(date_str, info, o) if f.exists()]
    if not files:
        return None
    t0 = time.process_time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        native = load_data(files, info.instrument_type, o.data_level)
    if native is None:
        return None
    Q.STATE["band"] = None
    Q.STATE["sza"] = None
    Q.STATE["capture"] = {}
    _TIMER["t"] = 0.0
    fi = {}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = calibrate_rayleigh(date_str, info, o, preloaded_data=native, fit_inputs_out=fi)
        flag, c, u = float(r.flag), float(r.lidar_constant), float(r.uncertainty)
        zb = (r.calibration_bottom_height - info.altitude
              if r.calibration_bottom_height is not None else np.nan)
        zt = (r.calibration_top_height - info.altitude
              if r.calibration_top_height is not None else np.nan)
        msg = r.message
    except Exception as exc:  # noqa: BLE001
        flag, c, u, zb, zt = -99.0, np.nan, np.nan, np.nan, np.nan
        msg = f"{type(exc).__name__}: {exc}"
    cpu = time.process_time() - t0
    b = fi.get("uncertainty_budget")
    row = dict(wmo=s["wmo"], ident=s["ident"], type=s["type"], date=date_str, flag=flag,
               c=c if flag in SUCCESS else np.nan, cpu_s=round(cpu, 3),
               cpu_budget_s=round(_TIMER["t"], 4), message=str(msg)[:140])
    if b is not None:
        # The budget is computed before the QC gates: keep it also for rejected nights
        cc = b.lidar_constant
        row.update(c_all=cc, u_ens=b.u_ensemble, u_rec=b.u_total, sig_ens=b.sigma_ensemble_rel,
                   sig_vert_rec=b.sigma_vertical_rel, sig_mol=b.sigma_molecular_rel,
                   sig_noise=b.sigma_noise_rel, g_rec=b.gradient_per_km,
                   g_sigma_rec=b.gradient_sigma_per_km, chi2_rec=b.chi2_red,
                   resid_rec=b.residual_rel, mol_bias=b.molecular_bias_rel, n_sub=b.n_subwindows,
                   z_lo=b.z_lo_m, z_hi=b.z_hi_m, sub=_enc_sub(b, cc))
        if not np.isfinite(zb):
            zb, zt = np.nan, np.nan
    row.update(zb_agl=zb, zt_agl=zt)
    times = Q.STATE["capture"].get("times")
    if times is not None and len(times):
        row.update(n_native=len(times), t_mean=float(np.mean(times)))
        if hk is not None:
            tm, _bg = hk.mean(times)
            row["tint_mean"] = tm
    if fi:
        row.update(Q._snr_bands(fi))
    return row


def _write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS + ["c_all"], extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS + ["c_all"]})
    os.replace(tmp, path)


def work_chunk(wmo, ident, dates, out_path, params=None):
    if out_path.exists():
        return str(out_path), 0, 0.0
    warnings.filterwarnings("ignore")
    import logging
    logging.disable(logging.CRITICAL)
    import q16a_recalibration as Q
    s = census()[(wmo, ident)]
    hk = Q.HK(wmo, ident)
    hk = hk if hk.ok else None
    rows = []
    t0 = time.time()
    for d in dates:
        try:
            r = calibrate_night(s, d, hk, params)
            if r is not None:
                rows.append(r)
        except Exception as exc:  # noqa: BLE001
            rows.append(dict(wmo=wmo, ident=ident, type=s["type"], date=d, flag=-98,
                             message=f"{type(exc).__name__}: {exc}"[:140]))
    _write_rows(out_path, rows)
    return str(out_path), len(rows), time.time() - t0


def _pool(tasks, workers):
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(work_chunk, *t) for t in tasks]
        for f in as_completed(futs):
            try:
                p, n, dt = f.result()
            except Exception:  # noqa: BLE001
                traceback.print_exc()
                continue
            done += 1
            print(f"[{done}/{len(tasks)}] {Path(p).parent.name}/{Path(p).name}: {n} rows in "
                  f"{dt:.0f} s (elapsed {time.time() - t0:.0f} s)", flush=True)


def unit_dates(wmo, ident, period_weeks, d0=datetime(2025, 1, 1), d1=datetime(2026, 7, 31)):
    """Nights D (YYYYMMDD) in 7-night blocks, one block every ``period_weeks`` weeks."""
    out = []
    for k in range((d1 - d0).days + 1):
        if (k // 7) % int(period_weeks) != 0:
            continue
        d = d0 + timedelta(days=k)
        f = L1_ROOT / wmo / f"{d:%Y}" / f"{d:%m}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"
        if f.exists():
            out.append(f"{d:%Y%m%d}")
    return out


def cmd_run(args):
    units = UNITS if not args.units else [tuple(u.split(":")) for u in args.units.split(",")]
    tasks = []
    for wmo, ident in units:
        period = 2 if (wmo, ident) in KEY_UNITS else 3
        dates = unit_dates(wmo, ident, period)
        # one task per 7-night block (keeps a block's shared files inside one worker)
        blocks = {}
        for d in dates:
            k = (datetime.strptime(d, "%Y%m%d") - datetime(2025, 1, 1)).days // 7
            blocks.setdefault(k, []).append(d)
        for k, ds in sorted(blocks.items()):
            tasks.append((wmo, ident, ds, RUN_DIR / f"{wmo}_{ident}" / f"{ds[0]}.csv"))
    if args.limit:
        tasks = tasks[:args.limit]
    print(f"{len(tasks)} tasks, {sum(len(t[2]) for t in tasks)} unit-nights", flush=True)
    _pool(tasks, args.workers)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=8)
    r.add_argument("--units", default="")
    r.add_argument("--limit", type=int, default=0)
    sub.add_parser("analyse")
    v = sub.add_parser("verify")
    v.add_argument("--workers", type=int, default=8)
    sub.add_parser("payerne")
    sub.add_parser("kalman")
    c = sub.add_parser("cost")
    c.add_argument("--n", type=int, default=60)
    args = ap.parse_args()
    if args.cmd == "run":
        cmd_run(args)
    else:
        import u23_analysis as A
        getattr(A, f"cmd_{args.cmd}")(args)


if __name__ == "__main__":
    main()
