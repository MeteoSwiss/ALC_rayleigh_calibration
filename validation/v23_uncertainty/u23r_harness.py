# -*- coding: utf-8 -*-
"""Independent review of the v2.3 uncertainty budget -- recalibration harness (u23r_*).

Why a second harness
--------------------
The u23 validation (u23_validation.py) recorded, per night, the eprof_v2.2 constant, the budget
ingredients and the C_L of 500 m sub-windows over a wide range. The review needs more, from the SAME
in-memory night (each night's L1 is read ONCE, every variant is computed in memory):

* both window methods: ``eprof_v2.2`` (the published constants) and ``eprof_v2.3`` (the release
  method, repaired tier 2), each with ``uncertainty_model = "v2.3"`` and the wide RECORDING range;
* the tier of the selected window (tier 1 = strict gates, tier 2 = noise-aware fallback), from the
  window-selection message;
* the EXACT shipped-default budget, obtained by calling the budget again on the captured inputs with
  ``uncertainty_params = {}`` (no offline emulation of the sub-window range);
* a measured-noise variant of the photon-noise term: the noise of the night mean measured on the
  binned stack (molecular_methods._binned_mean_noise, the eprof_v2.3 estimator) instead of the
  propagated sigma_signal / sqrt(k);
* split-night constants -- odd vs even 300-s profiles and first vs second half of the night -- in the
  operational window and in every recorded sub-window: (C_odd - C_even) / 2 is an independent
  realisation of the noise of C_L (fast noise only), (C_h1 - C_h2) / 2 adds slow changes;
* the three dimensions of the v2.2 ensemble separately (window shift, lidar ratio, time subsets with
  the code's own RNG seed), to quantify overlaps between the ensemble and the new terms.

Units
-----
``orig``   the 36 units of the u23 validation, same nights (blocks of 7, every 2nd/3rd week);
``fresh``  units NOT used to design or tune the budget: 27 CHM15k (14 countries), 8 CL61, 2 Mini-MPL.

Outputs only under C:/DATA/Projects/202610_v23/uncertainty_review/runs/<set>/<method>/<unit>/.

Usage
-----
    python validation/v23_uncertainty/u23r_harness.py run --set orig --workers 12
    python validation/v23_uncertainty/u23r_harness.py run --set fresh --workers 12
"""
from __future__ import annotations

import argparse
import copy
import csv
import os
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "chm15k_stability"))
sys.path.insert(0, str(HERE))

import u23_validation as V  # noqa: E402  (sets ALC_WV_SPECTRUM = the operational CL61 line)

OUT = Path("C:/DATA/Projects/202610_v23/uncertainty_review")
RUNS = OUT / "runs"
METHODS = ("eprof_v2.2", "eprof_v2.3")
SUCCESS = (1.0, 0.5)

# Fresh units: never used to design or tune the budget. Spread over countries and latitudes
# (12-69 N), coasts, plains and two mountain sites; not chosen by data volume.
FRESH_CHM15K = [
    ("0-100-20000-000627", "A"),          # Plovdiv (BG)
    ("0-20000-0-01010", "A"),             # Andoya (NO), 69 N
    ("0-20000-0-01311", "A"),             # Bergen Flesland (NO)
    ("0-20000-0-03005", "A"),             # Lerwick (UK)
    ("0-20000-0-03162", "A"),             # Eskdalemuir (UK)
    ("0-20000-0-03808", "A"),             # Camborne (UK), co-located CL61 C
    ("0-20000-0-06348", "A"),             # Cabauw (NL)
    ("0-20000-0-06784", "A"),             # Davos (CH), 1594 m
    ("0-20008-0-BRN", "A"),               # Bern (CH)
    ("0-20000-0-10140", "0"),             # Hamburg (DE)
    ("0-20000-0-10471", "0"),             # Leipzig (DE)
    ("0-20000-0-10739", "0"),             # Stuttgart (DE)
    ("0-20000-0-10895", "0"),             # Fuerstenzell (DE)
    ("0-20008-0-MEL", "A"),               # Melpitz (DE)
    ("0-20000-0-12982", "B"),             # Szeged (HU)
    ("0-348-1-64310", "A"),               # Debrecen (HU)
    ("0-20000-0-14280", "A"),             # Osijek (HR)
    ("0-20000-0-14328", "A"),             # Ogulin (HR)
    ("0-380-6-1", "A"),                   # Rome (IT)
    ("0-20008-0-CGR", "A"),               # Capo Granitola (IT), 37.6 N
    ("0-20008-0-CMN", "A"),               # Monte Cimone (IT), 2165 m
    ("0-20008-0-UVA", "A"),               # Valladolid (ES)
    ("0-20008-0-UEX", "A"),               # Badajoz (ES)
    ("0-724-60000-000000000UGRA001", "A"),  # Guadiana-UGR (ES), 37.9 N
    ("0-250-1001-75105006", "A"),         # QualAir Paris (FR)
    ("0-20000-0-73009", "A"),             # Western / London ON (CA)
    ("0-20000-0-78990", "A"),             # Bonaire (BQ), 12 N
]
FRESH_CL61 = [
    ("0-20000-0-03808", "C"), ("0-20000-0-06418", "B"), ("0-20000-0-06447", "B"),
    ("0-20008-0-BIR", "A"), ("0-203-10-LNG", "A"), ("0-20000-0-11538", "B"),
    ("0-756-4-EERLCL61", "A"), ("0-20008-0-LAU", "A"),
]
FRESH_MPL = [("0-20000-0-07145", "A"), ("0-20000-0-07617", "A")]

FIELDS = (list(V.FIELDS) + ["c_all", "method", "tier", "sel_msg", "n_prof", "k_time",
                            "native_per_bin", "sig_noise_prop", "sig_noise_meas", "noise_ratio_36",
                            "def_u", "def_sig_vert", "def_g", "def_gs", "def_res", "def_chi2",
                            "def_nsub", "def_zlo", "def_zhi", "def_sig_noise", "def_sig_mol",
                            "def_sig_aer", "def_exceeds", "c_odd", "c_even", "c_h1", "c_h2",
                            "sub_w", "sub_odd", "sub_even", "sub_h1", "sub_h2",
                            "ens_shift", "ens_lr", "ens_time", "sig_noise_wm", "sub_sig_wm",
                            "sub_lo", "sub_hi", "rev_u", "rev_sig_vert", "rev_g", "rev_gs",
                            "rev_res", "rev_nsub", "rev_sig_noise", "rev_sig_mol", "rev_exceeds"])
# The budget as first shipped (stored as def_*), pinned explicitly so the harness stays reproducible
# after the review changed DEFAULT_PARAMS (at the time of the review runs DEFAULT_PARAMS were these).
PRE_REVIEW = dict(molecular_height_reference="agl", noise_estimator="per_gate", vertical_cap_rel=None)
# The review's revision (exact code values, stored as rev_*): the molecular climatology at the window
# centre ASL (its fit convention) and the window-mean noise estimator. The saturation of the vertical
# term (vertical_cap_rel = 0.15) is a deterministic transform of rev_sig_vert and is applied offline
# (u23r_saturation.saturate), so the uncapped value stays available; rev_u / rev_exceeds are uncapped.
REVISED = dict(molecular_height_reference="asl", noise_estimator="window_mean", vertical_cap_rel=None)

# --------------------------------------------------------------------------------------------
# Capture of the budget inputs, the fit result and the window-selection message
# --------------------------------------------------------------------------------------------
CAP: dict = {}
_P = {"done": False, "orig_budget": None}


def install():
    """q16a patches (no NetCDF write, profile-time capture) + three capture wrappers."""
    Q = V._install()
    if _P["done"]:
        return Q
    import calibration.rayleigh.calibration as RC
    import calibration.rayleigh.molecular_methods as MM
    import calibration.rayleigh.uncertainty as U

    orig_find = RC.find_optimal_molecular_window

    def find(**kw):
        res = orig_find(**kw)
        CAP["fit"] = res
        return res

    orig_sel = MM.select_molecular_window

    def sel(*a, **k):
        mw = orig_sel(*a, **k)
        CAP["sel_msg"] = str(getattr(mw, "message", ""))
        return mw

    # Tier of the selected window. eprof_v2.2 gives a tier-2 window the same message as a tier-1
    # one, so the tier is decided directly: the strict tier alone (use_noise=False disables the
    # noise-aware fallback; the strict mask never reads the noise) either finds a window or not.
    orig_opt = MM._select_optimal

    def opt(g, **kw):
        mw = orig_opt(g, **kw)
        try:
            strict = orig_opt(g, **dict(kw, use_noise=False))
            CAP["tier"] = 1.0 if strict.ok else (2.0 if mw.ok else np.nan)
        except Exception:                  # noqa: BLE001 - a diagnostic must never lose a night
            pass
        return mw

    orig_budget = U.budget_for_night

    def bud(**kw):
        CAP["budget_kw"] = kw
        return orig_budget(**kw)

    RC.find_optimal_molecular_window = find
    MM.select_molecular_window = sel
    MM._select_optimal = opt
    U.budget_for_night = bud
    _P["orig_budget"] = orig_budget
    _P["done"] = True
    return Q


def options(method):
    o = V._options(None, "v2.3")                       # recording range, no QC limit
    o.molecular_method = method
    return o


def _enc(x):
    return ",".join("nan" if not np.isfinite(v) else f"{v:.5f}" for v in np.asarray(x, float))


# --------------------------------------------------------------------------------------------
# Extras computed in memory on the captured inputs
# --------------------------------------------------------------------------------------------
def extras(row, fi, kw, b):
    """Shipped-default budget, noise variants, split-night constants and ensemble parts."""
    import calibration.rayleigh.uncertainty as U
    from calibration.rayleigh.calibration import (_compute_cl_for_perturbation, LR_DELTAS,
                                                  ALT_SHIFTS_M, N_TIME_SAMPLES,
                                                  TIME_SUBSET_FRACTION)
    from calibration.rayleigh.molecular_methods import _binned_mean_noise

    o_rec = kw["options"]
    rng = np.asarray(kw["range_alc"], float)
    rcs_mean = np.asarray(kw["rcs_mean"], float)
    fit = kw["fit_result"]
    mp = kw["mol_props"]
    cl = float(kw["cl"])
    zb, zt = float(fit.range_start_m), float(fit.range_end_m)

    # (1) the shipped (pre-review) budget, on the identical inputs (exact, no offline emulation)
    o_def = copy.copy(o_rec)
    o_def.uncertainty_params = dict(PRE_REVIEW)
    bd = _P["orig_budget"](**dict(kw, options=o_def))
    row.update(def_u=bd.u_total, def_sig_vert=bd.sigma_vertical_rel, def_g=bd.gradient_per_km,
               def_gs=bd.gradient_sigma_per_km, def_res=bd.residual_rel, def_chi2=bd.chi2_red,
               def_nsub=bd.n_subwindows, def_zlo=bd.z_lo_m, def_zhi=bd.z_hi_m,
               def_sig_noise=bd.sigma_noise_rel, def_sig_mol=bd.sigma_molecular_rel,
               def_sig_aer=bd.sigma_aerosol_rel, def_exceeds=int(bool(bd.exceeds_limit)))
    # the candidate revision, same inputs (exact code values)
    o_rev = copy.copy(o_rec)
    o_rev.uncertainty_params = dict(REVISED)
    br = _P["orig_budget"](**dict(kw, options=o_rev))
    row.update(rev_u=br.u_total, rev_sig_vert=br.sigma_vertical_rel, rev_g=br.gradient_per_km,
               rev_gs=br.gradient_sigma_per_km, rev_res=br.residual_rel, rev_nsub=br.n_subwindows,
               rev_sig_noise=br.sigma_noise_rel, rev_sig_mol=br.sigma_molecular_rel,
               rev_exceeds=int(bool(br.exceeds_limit)))

    # (2) photon noise: propagated (budget) vs measured on the binned stack
    k = U.time_binning_factor(kw.get("sigma_native"), kw.get("hours"))
    row.update(k_time=k, native_per_bin=fi.get("native_per_bin"))
    signal = rcs_mean / rng ** 2
    in_win = (rng >= zb) & (rng <= zt) & np.isfinite(rcs_mean)
    sg = kw.get("sigma_signal")
    sg_c = None if sg is None else np.asarray(sg, float) / np.sqrt(k)
    gn = U.gate_relative_noise(signal, sg_c, mp.p_mol, in_win)
    if gn is None:
        gn = U._fallback_gate_noise(signal, mp.p_mol, rng, zb, zt)
    row["sig_noise_prop"] = U.window_noise_rel(gn, in_win)
    stack = fi.get("signal_stack")
    if stack is not None and np.ndim(stack) == 2:
        meas = _binned_mean_noise(np.asarray(stack, float))
        if meas is not None and meas.size == rng.size:
            e = U.gate_relative_noise(signal, meas, mp.p_mol, in_win)
            row["sig_noise_meas"] = U.window_noise_rel(e, in_win)
            if sg_c is not None:
                band = (rng >= 3000) & (rng <= 6000) & np.isfinite(meas) & (meas > 0) & np.isfinite(sg_c)
                if band.sum() >= 5:
                    row["noise_ratio_36"] = float(np.median(sg_c[band] / meas[band]))

    # (3) split-night constants in the window and in the recorded sub-windows
    if stack is not None and np.ndim(stack) == 2 and np.shape(stack)[0] >= 4:
        rcs_stack = np.asarray(stack, float) * rng[None, :] ** 2
        n = rcs_stack.shape[0]
        idx = np.arange(n)
        params = U.resolve_params(o_rec)
        cen_rec = np.asarray(b.sub_centre_m, float) if b is not None else np.array([])
        w = 500.0
        if cen_rec.size >= 2 and np.median(np.diff(cen_rec)) < 400.0:
            w = 250.0
        elif cen_rec.size < 3:
            w = 250.0
        row["sub_w"] = w
        valid = np.isfinite(rcs_mean)
        # window-mean noise measured on the binned stack (range-correlated noise kept), for the
        # operational window and every recorded sub-window (same tiling, boundaries recorded)
        med = getattr(o_rec, "time_aggregation", "mean") == "median"
        row["sig_noise_wm"] = U.window_mean_noise_rel(stack, mp.p_mol, in_win, med)
        if cen_rec.size:
            bounds = []

            def probe(a, bb):
                bounds.append((a, bb))
                return 1.0
            U.subwindow_constants(probe, rng, gn, valid, zb, zt, dict(params, subwindow_m=w))
            lut = {int(round(0.5 * (a + bb))): (a, bb) for a, bb in bounds}
            lo, hi, sw = [], [], []
            for c in cen_rec:
                a, bb = lut.get(int(round(c)), (np.nan, np.nan))
                lo.append(a)
                hi.append(bb)
                sw.append(U.window_mean_noise_rel(stack, mp.p_mol, (rng >= a) & (rng <= bb), med)
                          if np.isfinite(a) else np.nan)
            row["sub_lo"] = _enc(lo)
            row["sub_hi"] = _enc(hi)
            row["sub_sig_wm"] = _enc(sw)
        for name, ii in (("odd", idx[0::2]), ("even", idx[1::2]), ("h1", idx[: n // 2]),
                         ("h2", idx[n // 2:])):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rm = np.nanmean(rcs_stack[ii], axis=0)
            cfun = U.operational_cl_in_window(rm, rng, mp.beta_mol, fit, o_rec)
            c0 = cfun(zb, zt)
            row[f"c_{name}"] = (c0 / cl) if c0 is not None else np.nan
            if cen_rec.size:
                try:
                    cen, cls, _s, _a, _b = U.subwindow_constants(
                        cfun, rng, gn, valid, zb, zt, dict(params, subwindow_m=w))
                    lut = {int(round(c)): v for c, v in zip(cen, cls)}
                    vals = [lut.get(int(round(c)), np.nan) / cl for c in cen_rec]
                    row[f"sub_{name}"] = _enc(vals)
                except Exception:          # noqa: BLE001 - a diagnostic must never lose a night
                    pass

        # (4) the v2.2 ensemble dimensions, one at a time (nominal values for the others)
        def one(rm, lr, shift):
            pr = _compute_cl_for_perturbation(
                rcs_mean=rm, range_alc=rng, beta_mol=mp.beta_mol, fit_result=fit,
                lidar_ratio_aerosol=lr, altitude_shift_m=shift,
                subtract_background=o_rec.subtract_background,
                consider_points_lower_than_molecular=o_rec.consider_points_lower_than_molecular,
                sign_error_v10=o_rec.sign_error_v10)
            return np.nan if pr is None else pr.lidar_constant / cl
        lr0 = o_rec.lidar_ratio_aerosol
        row["ens_shift"] = _enc([one(rcs_mean, lr0, s) for s in ALT_SHIFTS_M])
        row["ens_lr"] = _enc([one(rcs_mean, lr0 + d, 0) if lr0 + d > 0 else np.nan
                              for d in LR_DELTAS])
        gen = np.random.default_rng(42)                 # the code's own seed and draw order
        n_keep = min(n, max(3, int(n * TIME_SUBSET_FRACTION)))
        tv = [one(rcs_mean, lr0, 0)]
        for _ in range(N_TIME_SAMPLES):
            sel = np.sort(gen.choice(n, size=n_keep, replace=False))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tv.append(one(np.nanmean(rcs_stack[sel], axis=0), lr0, 0))
        row["ens_time"] = _enc(tv)


# --------------------------------------------------------------------------------------------
# One night, both methods, one L1 read
# --------------------------------------------------------------------------------------------
def night_rows(s, date_str, hk, Q):
    from calibration import calibrate_rayleigh
    from calibration.io.data_loader import build_file_paths, load_data
    info = V._info(s)
    files = [f for f in build_file_paths(date_str, info, options(METHODS[0])) if f.exists()]
    if not files:
        return []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        native = load_data(files, info.instrument_type, options(METHODS[0]).data_level)
    if native is None:
        return []
    rows = []
    for method in METHODS:
        o = options(method)
        CAP.clear()
        Q.STATE["band"] = None
        Q.STATE["sza"] = None
        Q.STATE["capture"] = {}
        fi = {}
        t0 = time.process_time()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = calibrate_rayleigh(date_str, info, o, preloaded_data=copy.deepcopy(native),
                                       fit_inputs_out=fi)
            flag, c = float(r.flag), float(r.lidar_constant)
            zb = (r.calibration_bottom_height - info.altitude
                  if r.calibration_bottom_height is not None else np.nan)
            zt = (r.calibration_top_height - info.altitude
                  if r.calibration_top_height is not None else np.nan)
            msg = r.message
        except Exception as exc:  # noqa: BLE001
            flag, c, zb, zt = -99.0, np.nan, np.nan, np.nan
            msg = f"{type(exc).__name__}: {exc}"
        row = dict(wmo=s["wmo"], ident=s["ident"], type=s["type"], date=date_str, flag=flag,
                   c=c if flag in SUCCESS else np.nan, cpu_s=round(time.process_time() - t0, 3),
                   message=str(msg)[:140], method=method, tier=CAP.get("tier", np.nan),
                   sel_msg=str(CAP.get("sel_msg", ""))[:100])
        b = fi.get("uncertainty_budget")
        if b is not None:
            row.update(c_all=b.lidar_constant, u_ens=b.u_ensemble, u_rec=b.u_total,
                       sig_ens=b.sigma_ensemble_rel, sig_vert_rec=b.sigma_vertical_rel,
                       sig_mol=b.sigma_molecular_rel, sig_noise=b.sigma_noise_rel,
                       g_rec=b.gradient_per_km, g_sigma_rec=b.gradient_sigma_per_km,
                       chi2_rec=b.chi2_red, resid_rec=b.residual_rel,
                       mol_bias=b.molecular_bias_rel, n_sub=b.n_subwindows, z_lo=b.z_lo_m,
                       z_hi=b.z_hi_m, sub=V._enc_sub(b, b.lidar_constant))
            # the window of the budget (also for nights a later QC gate rejected)
            fr = CAP.get("budget_kw", {}).get("fit_result")
            if fr is not None and not np.isfinite(zb):
                zb, zt = float(fr.range_start_m), float(fr.range_end_m)
        row.update(zb_agl=zb, zt_agl=zt)
        if fi.get("signal_stack") is not None:
            row["n_prof"] = int(np.shape(fi["signal_stack"])[0])
        times = Q.STATE["capture"].get("times")
        if times is not None and len(times):
            row.update(n_native=len(times), t_mean=float(np.mean(times)))
            if hk is not None:
                row["tint_mean"] = hk.mean(times)[0]
        if fi:
            row.update(Q._snr_bands(fi))
        kw = CAP.get("budget_kw")
        if kw is not None and b is not None:
            try:
                extras(row, fi, kw, b)
            except Exception as exc:  # noqa: BLE001
                row["sel_msg"] = (row["sel_msg"] + f" | extras failed: {type(exc).__name__}")[:140]
        rows.append(row)
    return rows


def _write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})
    os.replace(tmp, path)


def work_chunk(set_name, wmo, ident, dates):
    """One 7-night block of one unit -> one CSV per method."""
    outs = {m: RUNS / set_name / m / f"{wmo}_{ident}" / f"{dates[0]}.csv" for m in METHODS}
    if all(p.exists() for p in outs.values()):
        return f"{wmo}_{ident}/{dates[0]}", 0, 0.0
    warnings.filterwarnings("ignore")
    import logging
    logging.disable(logging.CRITICAL)
    Q = install()
    s = V.census()[(wmo, ident)]
    hk = Q.HK(wmo, ident)
    hk = hk if hk.ok else None
    acc = {m: [] for m in METHODS}
    t0 = time.time()
    for d in dates:
        try:
            for r in night_rows(s, d, hk, Q):
                acc[r["method"]].append(r)
        except Exception as exc:  # noqa: BLE001
            for m in METHODS:
                acc[m].append(dict(wmo=wmo, ident=ident, type=s["type"], date=d, flag=-98,
                                   method=m, message=f"{type(exc).__name__}: {exc}"[:140]))
    for m, p in outs.items():
        _write_rows(p, acc[m])
    return f"{wmo}_{ident}/{dates[0]}", sum(len(v) for v in acc.values()), time.time() - t0


# --------------------------------------------------------------------------------------------
# Night sampling (same calendar blocks as the u23 validation)
# --------------------------------------------------------------------------------------------
def unit_dates(wmo, ident, period_weeks, d0, d1):
    out = []
    for k in range((d1 - d0).days + 1):
        if (k // 7) % int(period_weeks) != 0:
            continue
        d = d0 + timedelta(days=k)
        f = V.L1_ROOT / wmo / f"{d:%Y}" / f"{d:%m}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"
        if f.exists():
            out.append(f"{d:%Y%m%d}")
    return out


def tasks_for(set_name):
    d0 = datetime(2025, 1, 1)
    plan = []
    if set_name in ("orig", "origx"):          # origx = same nights, harness with every field
        for wmo, ident in V.UNITS:
            period = 2 if (wmo, ident) in V.KEY_UNITS else 3
            plan.append((wmo, ident, period, d0, datetime(2026, 7, 31)))
    elif set_name == "fresh":
        for wmo, ident in FRESH_CHM15K:
            plan.append((wmo, ident, 2, d0, datetime(2026, 8, 15)))
        for wmo, ident in FRESH_CL61 + FRESH_MPL:
            # Lauder (NZ) is outside the European CAMS box; its own box exists from 2026-06-24 only
            # (a 910 nm night without CAMS is never calibrated), so earlier nights would only cost
            # L1 reads.
            first = datetime(2026, 6, 24) if wmo == "0-20008-0-LAU" else d0
            plan.append((wmo, ident, 1, first, datetime(2026, 8, 15)))
    else:
        raise ValueError(set_name)
    tasks = []
    for wmo, ident, period, a, b in plan:
        blocks = {}
        for d in unit_dates(wmo, ident, period, a, b):
            k = (datetime.strptime(d, "%Y%m%d") - d0).days // 7
            blocks.setdefault(k, []).append(d)
        for _k, ds in sorted(blocks.items()):
            tasks.append((set_name, wmo, ident, ds))
    return tasks


def cmd_run(args):
    tasks = tasks_for(args.set)
    if args.units:
        keep = set(args.units.split(","))
        tasks = [t for t in tasks if f"{t[1]}_{t[2]}" in keep]
    if args.limit:
        tasks = tasks[:args.limit]
    # interleave units so concurrent workers read different stations' files
    tasks.sort(key=lambda t: (t[3][0], t[1], t[2]))
    print(f"{len(tasks)} tasks, {sum(len(t[3]) for t in tasks)} unit-nights", flush=True)
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(work_chunk, *t) for t in tasks]
        for f in as_completed(futs):
            try:
                name, n, dt = f.result()
            except Exception:  # noqa: BLE001
                traceback.print_exc()
                continue
            done += 1
            if done % 20 == 0 or done == len(tasks):
                print(f"[{done}/{len(tasks)}] {name}: {n} rows in {dt:.0f} s "
                      f"(elapsed {time.time() - t0:.0f} s)", flush=True)


def cmd_check(args):
    """Blocks to recompute: fewer rows than expected nights (an L1 file missing while the block
    ran), or written inside a time window for a station whose L1 folder was incomplete then
    (2026-10-21 20:15:20-20:16:17: the Schiphol folder was moved aside and back). --delete removes
    them (both methods) so that a new `run` recomputes exactly those blocks."""
    t_from = datetime.fromisoformat(args.t_from).timestamp() if args.t_from else None
    t_to = datetime.fromisoformat(args.t_to).timestamp() if args.t_to else None
    bad = []
    for set_name, wmo, ident, dates in tasks_for(args.set):
        for m in METHODS:
            p = RUNS / args.set / m / f"{wmo}_{ident}" / f"{dates[0]}.csv"
            if not p.exists():
                continue
            with open(p, encoding="utf-8") as f:
                n_rows = sum(1 for _ in f) - 1
            in_window = (t_from is not None and wmo.startswith(args.prefix or "\0")
                         and t_from <= p.stat().st_mtime <= t_to)
            if n_rows < len(dates) or in_window:
                bad.append((wmo, ident, dates[0], n_rows, len(dates), in_window))
                break
    for b in bad:
        print(f"{b[0]}_{b[1]} {b[2]}: {b[3]}/{b[4]} rows{' (written in the window)' if b[5] else ''}")
    print(f"{len(bad)} blocks to recompute")
    if args.delete:
        for wmo, ident, d0, *_ in bad:
            for m in METHODS:
                p = RUNS / args.set / m / f"{wmo}_{ident}" / f"{d0}.csv"
                if p.exists():
                    p.unlink()
        print("deleted")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--set", choices=("orig", "origx", "fresh"), required=True)
    r.add_argument("--workers", type=int, default=12)
    r.add_argument("--units", default="")
    r.add_argument("--limit", type=int, default=0)
    c = sub.add_parser("check")
    c.add_argument("--set", choices=("orig", "origx", "fresh"), required=True)
    c.add_argument("--prefix", default="", help="station prefix of the time-window test")
    c.add_argument("--t_from", default="", help="local time, ISO format")
    c.add_argument("--t_to", default="")
    c.add_argument("--delete", action="store_true")
    args = ap.parse_args()
    if args.cmd == "run":
        cmd_run(args)
    elif args.cmd == "check":
        cmd_check(args)


if __name__ == "__main__":
    main()
