# -*- coding: utf-8 -*-
"""q16e check 3: re-run the operational Rayleigh chain on sampled CHM15k nights and decompose C_L.

The calibration package is NOT modified: `_compute_cl_for_perturbation` is wrapped (monkeypatched in
this process only) so the nominal Klett solution of each night can be inspected. Outputs go ONLY to
C:/DATA/Projects/202610_CHM15k_stability/q16/.

Per night we record
  * the re-run constant vs the published one (does the code reproduce calibration_version 220?);
  * the slope of the molecular fit (= C_L x T_aer^2 of the window) and C_L / slope = exp(2 tau_aer);
  * tau_aer split by range band (0-100 m, 100-300 m, 300-1000 m, 1 km-window bottom, in window);
  * the constant obtained when the aerosol extinction below z0 (100 m, 300 m) is replaced by its value
    at z0 -- the option `z_start_ext` (100 m) exists in CalibrationOptions but is never used, the
    Klett integral starts at the first range gate;
  * the night length, the number of profiles, the window bottom/top.
Annual harmonics of each quantity are then fitted per unit (DOY), so an in-code SEASONAL term (Klett
transmission, near-range garbage) can be compared with the ~+-9 % unexplained cycle.
"""
from __future__ import annotations

import csv
import logging
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import L1_ROOT, OUT, chm15k_streams, load_constants  # noqa: E402

Q16 = OUT / "q16"
EPOCH = datetime(1970, 1, 1)
N_UNITS = 10
NIGHTS_PER_MONTH = 2
PAYERNE = "0-20000-0-06610_A"


def _trap_cum(y, x):
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


def one_night(args):
    warnings.filterwarnings("ignore")
    logging.disable(logging.CRITICAL)
    import calibration.rayleigh.calibration as cal
    from calibration import CalibrationOptions, InstrumentInfo
    from calibration.config import DataLevel, InstrumentType
    from calibration.rayleigh.atmosphere import MOLECULAR_LIDAR_RATIO

    s, ds, c_pub = args
    rec = []
    orig = cal._compute_cl_for_perturbation

    def wrapped(**kw):
        kw2 = dict(kw)
        kw2["return_diagnostics"] = True
        pr = orig(**kw2)
        rec.append((kw, pr))
        return pr

    cal._compute_cl_for_perturbation = wrapped
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.folder_output = Q16 / "q16e_recal"
    o.plot_main = False
    o.plot_all = False
    info = InstrumentInfo(site_name=s["site"], wmo_id=s["wmo"], identifier=s["ident"],
                          instrument_type=InstrumentType.CHM15k, latitude=s["lat"], longitude=s["lon"],
                          altitude=s.get("alt", 0.0))
    try:
        r = cal.calibrate_rayleigh(ds, info, o)
    except Exception as exc:  # noqa: BLE001
        return dict(stream=f"{s['wmo']}_{s['ident']}", date=ds, err=str(exc)[:80])
    row = dict(stream=f"{s['wmo']}_{s['ident']}", date=ds, lat=s["lat"], c_pub=c_pub,
               c_rerun=r.lidar_constant, flag=r.flag)
    if not rec or rec[-1][1] is None or r.lidar_constant <= 0:
        return row
    kw, pr = rec[-1]                         # the nominal (LR, shift 0, all profiles) call
    rng = np.asarray(kw["range_alc"], float)
    bm = np.asarray(kw["beta_mol"], float)
    fr = kw["fit_result"]
    rcs = np.asarray(kw["rcs_mean"], float)
    ext_aer = np.asarray(pr.ext_tot, float) - bm * MOLECULAR_LIDAR_RATIO
    win = (rng >= fr.range_start_m) & (rng <= fr.range_end_m) & np.isfinite(rcs)
    i0 = int(np.where(win)[0][0])
    i1 = int(np.where(win)[0][-1])
    tau = _trap_cum(np.nan_to_num(ext_aer), rng)

    def band(lo, hi):
        sel = (rng >= lo) & (rng < hi) & (np.arange(rng.size) < i0)
        if sel.sum() < 2:
            return 0.0
        return float(np.trapezoid(np.nan_to_num(ext_aer[sel]), rng[sel]))

    def cl_with(ext_a):
        ext_t = ext_a + bm * MOLECULAR_LIDAR_RATIO
        res = cal.calculate_lidar_constant(rcs_mean=rcs, beta_tot=pr.beta_tot, ext_tot=ext_t, range_alc=rng,
                                           molecular_mask=win, fit_result=fr, subtract_background=False)
        return res.lidar_constant

    c_nom = cl_with(ext_aer)
    alt = {}
    for z0 in (100.0, 300.0):
        e = ext_aer.copy()
        j = int(np.searchsorted(rng, z0))
        e[:j] = e[j]
        alt[z0] = cl_with(e)
    c_noaer = cl_with(np.zeros_like(ext_aer))
    t_prof = [kw2 for kw2, _ in rec]
    row.update(
        c_nom=c_nom, slope=fr.slope, zb=fr.range_start_m, zt=fr.range_end_m,
        tau_win_bottom=float(tau[i0]), tau_win_mid=float(tau[(i0 + i1) // 2]),
        tau_0_100=band(0, 100), tau_100_300=band(100, 300), tau_300_1000=band(300, 1000),
        tau_1000_win=band(1000, 1e6),
        r_z100=alt[100.0] / c_nom, r_z300=alt[300.0] / c_nom, r_noaer=c_noaer / c_nom,
        first_gate=float(rng[0]), n_calls=len(t_prof),
        rcs_0_100_over_win=float(np.nanmean(rcs[rng < 100] / bm[rng < 100]) / np.nanmean(rcs[win] / bm[win])),
    )
    return row


def pick_jobs():
    import json
    census = {f"{s['wmo']}_{s['ident']}": s for s in json.load(open(REPO / "validation/scope_l1_2026_census.json",
                                                                  encoding="utf-8"))}
    cands = []
    for wmo, ident, site, lat, lon in chm15k_streams():
        key = f"{wmo}_{ident}"
        if key == PAYERNE:
            continue
        c = load_constants(wmo, ident)
        if c is None:
            continue
        ok = (c["c"] > 0) & (c["t"] >= (datetime(2025, 1, 1) - EPOCH).days) & (c["t"] < (datetime(2026, 1, 1) - EPOCH).days)
        mon = np.array([(EPOCH + timedelta(days=float(x))).month for x in c["t"][ok]])
        if ok.sum() >= 150 and len(np.unique(mon)) == 12:
            cands.append((lon, key, c, ok))
    cands.sort()
    idx = np.linspace(0, len(cands) - 1, N_UNITS).round().astype(int)
    rng = np.random.default_rng(16)
    jobs = []
    for k in idx:
        lon, key, c, ok = cands[k]
        s = census[key]
        t, t0, cc = c["t"][ok], c["t0"][ok], c["c"][ok]
        mon = np.array([(EPOCH + timedelta(days=float(x))).month for x in t])
        for m in range(1, 13):
            sel = np.where(mon == m)[0]
            for j in rng.choice(sel, size=min(NIGHTS_PER_MONTH, sel.size), replace=False):
                d = EPOCH + timedelta(days=float(np.floor(t0[j] + 0.5)))
                jobs.append((s, d.strftime("%Y%m%d"), float(cc[j])))
    return jobs


def annual_fit(doy, y):
    w = 2 * np.pi / 365.25
    A = np.c_[np.ones_like(doy), np.cos(w * doy), np.sin(w * doy)]
    ok = np.isfinite(y)
    if ok.sum() < 8:
        return np.nan, np.nan
    p, *_ = np.linalg.lstsq(A[ok], y[ok], rcond=None)
    amp = np.hypot(p[1], p[2])
    peak = (np.arctan2(p[2], p[1]) / w) % 365.25
    return amp, peak


def main():
    (Q16 / "q16e_recal").mkdir(parents=True, exist_ok=True)
    jobs = pick_jobs()
    print(f"{len(jobs)} nights from {len({j[0]['wmo'] + j[0]['ident'] for j in jobs})} units")
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = list(ex.map(one_night, jobs))
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(Q16 / "q16e_klett_rerun.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    ok = [r for r in rows if "c_nom" in r]
    print(f"re-run succeeded on {len(ok)} / {len(rows)} nights")
    rat = np.array([r["c_rerun"] / r["c_pub"] for r in ok])
    print(f"re-run / published: median {np.median(rat):.4f}, within 2 %: {np.mean(np.abs(rat - 1) < 0.02) * 100:.0f} %")
    doy = np.array([datetime.strptime(r["date"], "%Y%m%d").timetuple().tm_yday for r in ok], float)
    units = sorted({r["stream"] for r in ok})
    print("quantity: network median of per-unit annual amplitude (peak DOY)")
    for q, f in (("ln C/slope = 2 tau_aer (in-code aerosol transmission)", lambda r: np.log(r["c_nom"] / r["slope"])),
                 ("2 tau_aer 0-100 m", lambda r: 2 * r["tau_0_100"]),
                 ("2 tau_aer 0-300 m", lambda r: 2 * (r["tau_0_100"] + r["tau_100_300"])),
                 ("ln r_z100", lambda r: np.log(r["r_z100"])),
                 ("ln r_z300", lambda r: np.log(r["r_z300"])),
                 ("ln r_noaer", lambda r: np.log(r["r_noaer"])),
                 ("ln slope", lambda r: np.log(r["slope"])),
                 ("ln C_pub", lambda r: np.log(r["c_pub"]))):
        amps, peaks, meds = [], [], []
        for u in units:
            s = [i for i, r in enumerate(ok) if r["stream"] == u]
            y = np.array([f(ok[i]) for i in s], float)
            a, p = annual_fit(doy[s], y)
            amps.append(a)
            peaks.append(p)
            meds.append(np.nanmedian(y))
        print(f"  {q:55s} median {100 * np.nanmedian(meds):+6.2f} %  amp {100 * np.nanmedian(amps):5.2f} %  "
              f"peak DOY {np.nanmedian(peaks):5.0f}")


if __name__ == "__main__":
    main()
