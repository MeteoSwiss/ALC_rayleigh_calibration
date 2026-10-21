# -*- coding: utf-8 -*-
"""ref23 harness: one read of every co-located instrument per night, everything computed in memory.

For every night D (2025-01-01 -> latest L1) of the six reference sites (ref23_common.SITES) each
instrument's L1 files [D-1, D] are read ONCE (calibration.io.instrument_day.load_instrument_day, the
runner's read-once path) and, in memory:

  1. Rayleigh (CHM15k, CL61), twice, with the operational code (calibrate_rayleigh unmodified,
     options.json, the NetCDF write patched out as in the q16a / u23 harnesses):
       v22  molecular_method eprof_v2.2, uncertainty_model "v2.3" with the relative QC limit
            switched off (max_rel_uncertainty 1e9): the constant, window and flags of the published
            v2.2 calibration (calibration_version 220) PLUS the v2.3 budget on that window. Its v2.2
            ensemble uncertainty is the budget's u_ensemble (= the published
            lidar_constant_uncertainty); the v2.2 'u > C' gate (-6) is re-applied offline.
       v23  molecular_method eprof_v2.3, uncertainty_model "v2.3", shipped defaults (the release
            candidate: repaired tier 2, full budget, relative QC limit 60 %).
     Every budget component (1 sigma, relative) is stored, so the totals can be recomputed if the
     budget is revised. The tier (strict = 1 / noise-aware fallback = 2) is captured from the window
     search message.
  2. Liquid-cloud (O'Connor) calibration of UTC day D (CL61, CL31), exactly as
     scripts/run_network_calibration._do_cloud does it (shared read sliced to day D), with the
     per-profile coefficients kept so the reported uncertainty can be examined.
  3. The inter-comparison ingredients of night D (v3 dashboard conventions, restricted to the
     calibration night = SZA > 100 deg): hourly medians of rcs_0 on the v3 display grid after a
     SITE cloud/fog screen (union of every co-located instrument's cloud base / Vaisala vertical
     visibility, expanded +-15 min), the SNR >= 3 admission masks (5-min and 60-min windows,
     nf_v3 semantics, per instrument and per pair = intersection), and the per-night CAMS
     corrections of the 910 nm units (two-way water-vapour transmission with the unit's own laser
     line, molecular attenuated backscatter at its wavelength and at 1064 nm).

The CL61 runs with the operational water-vapour line ALC_WV_SPECTRUM = {"CL61": [910.55, 0.188]}
(dark-free, as operational).  CAMS 0.4 deg only ("D:/CAMS_Monthly_04;D:/CAMS_daily").

Output: one pickle per site-month under C:/DATA/Projects/202610_v23/reference_agreement/runs/
(resumable: an existing pickle is skipped).

Usage
  python validation/v23_reference/ref23_run.py run --workers 12 [--sites payerne,lindenberg]
  python validation/v23_reference/ref23_run.py test --site payerne --date 20260715
"""
from __future__ import annotations

import os

# Environment BEFORE any numpy / calibration import (spawned workers re-import this module).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
os.environ["ALC_WV_SPECTRUM"] = '{"CL61": [910.55, 0.188]}'      # operational CL61 line
os.environ["ALC_VAL_CAMS_04"] = "D:/CAMS_Monthly_04"               # intercompare monthly 0.4 deg
os.environ["ALC_VAL_CAMS_1DEG"] = "D:/CAMS_Monthly_04"             # never the 1-degree D:/CAMS
os.environ["ALC_VAL_CAMS_DAILY"] = "D:/CAMS_daily"
os.environ["ALC_VAL_NOISE_FILTER"] = "none"

import argparse  # noqa: E402
import logging  # noqa: E402
import math  # noqa: E402
import pickle  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from datetime import datetime, timedelta  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for p in (REPO, REPO / "validation" / "chm15k_stability", REPO / "inter-comparison_dashboard", HERE):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import ref23_common as RC  # noqa: E402

# Recording configuration of the two Rayleigh variants (see the module docstring).
VARIANTS = {
    "v22": dict(method="eprof_v2.2", params={"max_rel_uncertainty": 1e9}),
    "v23": dict(method="eprof_v2.3", params={}),
}
_TIER = {"msg": "", "ok": False, "strict_ok": None}
_PATCHED = {"done": False}


# ------------------------------------------------------------------------------------- set-up
def _install():
    """q16a patches (no NetCDF write, profile-time capture) + capture of the window-search verdict.

    The tier cannot be read from the selection message: under eprof_v2.2 a fallback (tier-2) pick
    and a strict (tier-1) pick both end with "max composite quality". The selector is therefore
    wrapped and re-run with the noise-aware tier disabled (use_noise=False): if that strict-only
    selection finds a window the night is tier 1 (and it is the same window, since tier 1 never
    reads the noise), otherwise the calibrated window came from the fallback tier (tier 2)."""
    import q16a_recalibration as Q
    Q._install_patches()
    if not _PATCHED["done"]:
        from calibration.rayleigh import molecular_methods as MM
        orig = MM._select_optimal

        def select(g, **kw):
            mw = orig(g, **kw)
            _TIER["msg"] = str(getattr(mw, "message", ""))
            _TIER["ok"] = bool(getattr(mw, "ok", False))
            try:
                strict = orig(g, **dict(kw, use_noise=False))
                _TIER["strict_ok"] = bool(strict.ok)
                if strict.ok and mw.ok and not (strict.start_m == mw.start_m
                                                and strict.end_m == mw.end_m):
                    _TIER["strict_ok"] = None          # never expected: flag it
            except Exception:  # noqa: BLE001 - a diagnostic must never fail the night
                _TIER["strict_ok"] = None
            return mw
        MM._select_optimal = select
        for k, v in list(MM._SELECTORS.items()):
            if v is orig:
                MM._SELECTORS[k] = select
        _PATCHED["done"] = True
    return Q


def _info(s):
    from calibration import InstrumentInfo
    from calibration.config import InstrumentType
    itype = {"CHM15k": InstrumentType.CHM15k, "CL61": InstrumentType.CL61,
             "CL31": InstrumentType.CL31}[s["type"]]
    return InstrumentInfo(site_name=s.get("site", s["wmo"]), wmo_id=s["wmo"], identifier=s["ident"],
                          instrument_type=itype, latitude=s["lat"], longitude=s["lon"],
                          altitude=s.get("alt", 0.0))


def _options(variant):
    from calibration import CalibrationOptions, DataLevel
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = RC.L1_ROOT
    o.data_level = DataLevel.L1
    o.folder_output = RC.OUT / "nowrite"          # never written (the write is patched out)
    o.cams_folder = RC.CAMS
    o.plot_all = False
    o.plot_main = False
    o.molecular_method = VARIANTS[variant]["method"]
    o.uncertainty_model = "v2.3"
    o.uncertainty_params = dict(VARIANTS[variant]["params"])
    return o


# ------------------------------------------------------------------------------------- Rayleigh
def rayleigh_night(native, date_str, s, variant):
    """One Rayleigh calibration of an in-memory night -> row dict (all budget components)."""
    Q = _install()
    from calibration import calibrate_rayleigh
    info = _info(s)
    o = _options(variant)
    Q.STATE["band"] = None
    Q.STATE["sza"] = None
    Q.STATE["capture"] = {}
    _TIER["msg"], _TIER["ok"], _TIER["strict_ok"] = "", False, None
    fi = {}
    t0 = time.perf_counter()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = calibrate_rayleigh(date_str, info, o, preloaded_data=native, fit_inputs_out=fi)
        flag, c, msg = float(r.flag), float(r.lidar_constant), str(r.message)
        zb = (r.calibration_bottom_height - info.altitude
              if r.calibration_bottom_height is not None else np.nan)
        zt = (r.calibration_top_height - info.altitude
              if r.calibration_top_height is not None else np.nan)
    except Exception as exc:  # noqa: BLE001 - one bad night must not lose the task
        flag, c, zb, zt = -99.0, np.nan, np.nan, np.nan
        msg = f"{type(exc).__name__}: {exc}"
    wmsg = _TIER["msg"]
    # 1 = strict tier, 2 = noise-aware fallback tier, 0 = no window, -1 = undetermined
    if not _TIER["ok"]:
        tier = 0
    elif _TIER["strict_ok"] is None:
        tier = -1
    else:
        tier = 1 if _TIER["strict_ok"] else 2
    row = dict(flag=flag, c=(c if flag in RC.SUCCESS else np.nan), zb=zb, zt=zt, tier=tier,
               window_msg=wmsg[:160], msg=msg[:200], cpu=round(time.perf_counter() - t0, 3))
    b = fi.get("uncertainty_budget")
    if b is not None:
        # computed before the QC gates: kept for rejected nights too (c_all)
        row.update(c_all=float(b.lidar_constant), u_ens=float(b.u_ensemble),
                   u_total=float(b.u_total), sig_ens=float(b.sigma_ensemble_rel),
                   sig_vert=float(b.sigma_vertical_rel), sig_mol=float(b.sigma_molecular_rel),
                   sig_aer=float(b.sigma_aerosol_rel), sig_noise=float(b.sigma_noise_rel),
                   g=float(b.gradient_per_km), g_sigma=float(b.gradient_sigma_per_km),
                   resid=float(b.residual_rel), chi2=float(b.chi2_red),
                   mol_bias=float(b.molecular_bias_rel), n_sub=int(b.n_subwindows),
                   z_lo=float(b.z_lo_m), z_hi=float(b.z_hi_m), exceeds=bool(b.exceeds_limit),
                   budget_model=str(b.model),
                   sub_z=np.asarray(b.sub_centre_m, "f4"),
                   sub_c=np.asarray(b.sub_cl / b.lidar_constant if b.lidar_constant > 0
                                    else b.sub_cl, "f4"),
                   sub_s=np.asarray(b.sub_sigma_rel, "f4"))
    times = Q.STATE["capture"].get("times")
    if times is not None and len(times):
        row.update(n_prof=int(len(times)), t_first=float(np.min(times)),
                   t_last=float(np.max(times)), t_mean=float(np.mean(times)))
    return row


# ------------------------------------------------------------------------------------- cloud
def cloud_day(idd, d, s, fp):
    """Liquid-cloud calibration of UTC day d, as scripts/run_network_calibration._do_cloud."""
    from calibration.cloud import CloudCalConfig
    from calibration.cloud.calibration import (liquid_cloud_calibration_from_data, set_defaults,
                                               build_cloud_input_from_day)
    from calibration.flags import cloud_flag, dominant_cloud_reject_flag
    t0 = time.perf_counter()
    row = dict(flag=np.nan, c=np.nan, u=np.nan, n=0)
    try:
        cfg = set_defaults(CloudCalConfig(
            nc_file=str(fp), instrument=s["type"], apply_wv_correction=True,
            apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
            cams_folder=RC.CAMS, abs_cs_lookup_table=str(RC.WV_LUT), cams_folder_fallback="",
            wv_source="cams", era5_cache="",
            station_latitude=s["lat"], station_longitude=s["lon"],
            average_time_s=0.0, average_range_m=0.0))
        data = build_cloud_input_from_day(idd.slice_to_date(d), cfg)
        beta = getattr(data, "beta", None)
        if beta is None or not np.any(np.isfinite(np.asarray(beta, float))):
            row.update(flag=0.0, msg="No data (no usable signal)")
            return row
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = liquid_cloud_calibration_from_data(data, cfg)
        n = int(getattr(res, "n_profiles", 0))
        coef = float(res.cal_median)
        std = float(res.cal_std)
        flag = cloud_flag(n, coef, std)
        if flag == -1:
            flag, _reason, _rej = dominant_cloud_reject_flag(
                getattr(res, "filter_stats", None), getattr(res, "cloud_stats", None),
                getattr(res, "consistency_stats", None))
            if flag <= -20:
                cb = np.asarray(getattr(data, "cbh", []), dtype=float)
                if not np.any(np.isfinite(cb) & (cb > 0)):
                    flag = -1.0
        cl = float(res.lidar_constant)
        ok = float(flag) in RC.SUCCESS and math.isfinite(cl) and cl > 0
        row.update(flag=float(flag), c=(cl if ok else np.nan),
                   u=(cl * std / coef if ok and math.isfinite(coef) and coef != 0 else np.nan),
                   n=n, coef=coef, std=std, cc_applied=float(getattr(res, "calibration_constant_applied",
                                                                     np.nan)),
                   wl=float(cfg.wavelength), fwhm=float(cfg.laser_fwhm))
        vc = np.asarray(getattr(res, "valid_coefficients", []), "f8")
        if vc.size:
            tv = np.asarray(getattr(res, "time", []))
            try:
                tv = (np.asarray(tv, "datetime64[ns]").astype("int64") / 86400e9).astype("f8")
            except (TypeError, ValueError):
                tv = np.asarray(tv, "f8")
            row.update(vc=vc.astype("f4"), vt=tv,
                       vcbh=np.asarray(getattr(res, "cbh", []), "f4"))
    except Exception as exc:  # noqa: BLE001
        msg = str(exc)
        low = msg.lower()
        flag = (-10.0 if ("cams data too far" in low or "outside cams domain" in low)
                else -4.0 if ("no cams file" in low or "cams_beta" in low) else -99.0)
        row.update(flag=flag, msg=f"{type(exc).__name__}: {msg[:160]}")
    row["cpu"] = round(time.perf_counter() - t0, 3)
    return row


# ------------------------------------------------------------------------------------- ingredients
def _filled(a, dtype="f4"):
    return np.ma.filled(np.ma.masked_invalid(np.ma.asarray(a).astype(dtype)), np.nan)


def _exclusion_flags(native, itype):
    """Per-sample cloud/fog flag (v3 screen): any cloud base in (0, 20 km), or -- Vaisala only,
    as the operational screen (InstrumentType.reports_vv_obscuration) -- a vertical visibility."""
    cbh = np.asarray(native.cbh, float)
    if cbh.ndim == 1:
        cbh = cbh[:, None]
    has_cloud = np.any(np.isfinite(cbh) & (cbh > 0) & (cbh < 20000), axis=1)
    if itype in ("CL31", "CL61") and native.vertical_visibility is not None:
        vv = np.asarray(native.vertical_visibility, float).reshape(-1)
        if vv.size == has_cloud.size:
            has_cloud |= np.isfinite(vv) & (vv >= 0)
    return has_cloud


def _hourly(t_days, X, excl, dt_s):
    """Hourly medians of the screened night samples (intercompare.retime_hourly semantics: an hour
    needs >= 30 min of samples; the median runs over the non-excluded ones)."""
    hid = np.floor(np.asarray(t_days) * 24.0 + 1e-9).astype("i8")
    hours = np.unique(hid)
    H = np.full((hours.size, X.shape[1]), np.nan, "f4")
    nsamp = np.zeros(hours.size, "i4")
    nclear = np.zeros(hours.size, "i4")
    for k, h in enumerate(hours):
        sel = hid == h
        nsamp[k] = int(sel.sum())
        if not np.isfinite(dt_s) or sel.sum() * dt_s < RC.MIN_COV_S:
            continue
        cl = sel & ~excl
        nclear[k] = int(cl.sum())
        if cl.any():
            with np.errstate(all="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore")
                H[k] = np.nanmedian(X[cl], axis=0)
    return hours, H, nsamp, nclear


def ingredients(loaded, d, site_key, meta):
    """Ratio ingredients of night d for every instrument of the site (see the module docstring)."""
    import nf_v3 as NF
    from validation.paper import intercompare as IC
    from calibration.io.cams import find_cams_file
    from calibration.water_vapor_correction.water_vapor import (
        two_way_wv_transmission, cams_point_too_far, laser_spectrum_for)
    ds = RC.ymd(d)
    z = RC.Z_AGL
    lat, lon, alt = meta["lat"], meta["lon"], meta["alt"]
    # --- site cloud/fog screen: union of every instrument's flags, expanded +-15 min ---------------
    ft = []
    for ident, (itype, idd) in loaded.items():
        n = idd.native
        fl = _exclusion_flags(n, itype)
        ft.append(np.asarray(n.time, float)[fl])
    flag_t = np.sort(np.concatenate(ft)) if ft else np.zeros(0)

    def excluded(t):
        if flag_t.size == 0:
            return np.zeros(t.size, bool)
        j = np.searchsorted(flag_t, t)
        lo = flag_t[np.clip(j - 1, 0, flag_t.size - 1)]
        hi = flag_t[np.clip(j, 0, flag_t.size - 1)]
        dmin = np.minimum(np.abs(t - lo), np.abs(hi - t)) * 86400.0
        return dmin <= RC.CLOUD_EXPAND_S

    inst, win = {}, {}
    for ident, (itype, idd) in loaded.items():
        n = idd.native
        t = np.asarray(n.time, float)
        nm = RC.night_mask(t, lat, lon, ds)
        rng = np.asarray(n.range_alc, float)
        g = rng <= z[-1] + 60.0
        if nm.sum() < 5:
            continue
        tn = t[nm]
        X = _filled(np.asarray(n.rcs)[nm][:, g] if not np.ma.isMaskedArray(n.rcs)
                    else n.rcs[nm][:, g])
        ex = excluded(tn)
        o = np.argsort(tn, kind="stable")
        tn, X, ex = tn[o], X[o], ex[o]
        dt_s = float(np.median(np.diff(tn)) * 86400.0) if tn.size > 1 else np.nan
        hours, Hn, nsamp, nclear = _hourly(tn, X, ex, dt_s)
        Hg = IC.regrid(Hn.astype("f8"), rng[g], z).astype("f4")
        Xs = X.copy()
        Xs[ex] = np.nan
        tsec = np.round(tn * 86400.0).astype("i8")
        w = {}
        for W in RC.WINDOWS_S:
            wids, med, sig, nfin, small = NF.stream_windows(tsec, Xs, W)
            keep = NF.window_keep(med, sig, nfin, small, None)
            rk = IC.regrid(keep.astype("f8") * nfin, rng[g], z)
            rn = IC.regrid(nfin, rng[g], z)
            rn0 = np.nan_to_num(rn)
            with np.errstate(all="ignore"):
                keepb = (rn0 > 0) & (np.nan_to_num(rk) / np.where(rn0 > 0, rn0, 1.0) >= RC.NF_FMIN)
            w[W] = (wids.astype("i8"), keepb, rn0.astype("f4"))
        inst[ident] = dict(itype=itype, hours=hours, H=Hg, nsamp=nsamp, nclear=nclear, dt=dt_s,
                           wl=(float(idd.wavelength_nm_file) if idd.wavelength_nm_file is not None
                               else float(idd.wavelength_nm)),
                           excl_frac=float(ex.mean()) if ex.size else np.nan)
        win[ident] = w
    if not inst:
        return None
    hours = np.unique(np.concatenate([v["hours"] for v in inst.values()]))
    nh, nz = hours.size, z.size
    for ident, v in inst.items():
        A = np.full((nh, nz), np.nan, "f4")
        A[np.searchsorted(hours, v["hours"])] = v["H"]
        v["H"] = A
        cnt = np.zeros(nh, "i4")
        cnt[np.searchsorted(hours, v["hours"])] = v["nclear"]
        v["nclear"] = cnt
        del v["hours"], v["nsamp"]

    def hour_admit(wids, W, num_w, den_w):
        num = np.zeros((nh, nz))
        den = np.zeros((nh, nz))
        hpos = {int(h): i for i, h in enumerate(hours)}
        hh = wids * W // 3600
        for j in range(wids.size):
            u = hpos.get(int(hh[j]))
            if u is None:
                continue
            num[u] += num_w[j]
            den[u] += den_w[j]
        with np.errstate(all="ignore"):
            return (den > 0) & (num / np.where(den > 0, den, 1.0) >= RC.NF_FMIN)

    # per-instrument admission (the literal "par instrument" filter)
    for ident, v in inst.items():
        for W in RC.WINDOWS_S:
            wids, keepb, rn = win[ident][W]
            v[f"adm{W}"] = hour_admit(wids, W, keepb * rn, rn)
    # per-pair admission: the v3 "common" (intersection) rule on the windows of the two units
    pairs = {}
    for x, r in RC.l1_pairs(site_key):
        if x not in inst or r not in inst:
            continue
        pp = {}
        for W in RC.WINDOWS_S:
            wx, kx, nx = win[x][W]
            wr, kr, nr = win[r][W]
            allw = np.union1d(wx, wr)
            px = np.zeros((allw.size, nz), bool)
            pr = np.zeros((allw.size, nz), bool)
            wt = np.zeros((allw.size, nz))
            ix = np.searchsorted(allw, wx)
            ir = np.searchsorted(allw, wr)
            px[ix] = (nx > 0) & kx
            pr[ir] = (nr > 0) & kr
            wt[ix] += nx
            wt[ir] += nr
            pp[f"adm{W}"] = hour_admit(allw, W, (px & pr) * wt, wt)
        pairs[f"{x}|{r}"] = pp
    # --- per-night CAMS corrections of the 910 nm units ---------------------------------------------
    corr = {}
    cams = find_cams_file(RC.CAMS, ds)
    t0 = np.datetime64(d - timedelta(days=1)) + np.timedelta64(12, "h")
    t1 = np.datetime64(d) + np.timedelta64(12, "h")
    for ident, v in inst.items():
        lam = v["wl"]
        if not (np.isfinite(lam) and abs(lam - RC.TARGET_NM) >= 1.0):
            continue
        c = dict(cams=(str(cams) if cams is not None else ""), t2=None, bml=None, bmt=None)
        alt_asl = z + alt
        if cams is not None and not cams_point_too_far(cams, lat, lon):
            try:
                prof = IC._cams_wv_profile_cached(str(cams), lat, lon, t0, t1)
                if prof is not None:
                    lam0, fwhm = laser_spectrum_for(v["itype"], 910.0)
                    t2 = np.asarray(two_way_wv_transmission(alt_asl, alt, prof[0], prof[1],
                                                            IC.WV_LUT, lam0, fwhm), "f8")
                    if t2.size == z.size and np.any(np.isfinite(t2)):
                        c.update(t2=t2.astype("f4"), lam0=lam0, fwhm=fwhm)
                bml = IC._mol_att_cams(str(cams), lat, lon, t0, t1, alt_asl, alt, lam)
                bmt = IC._mol_att_cams(str(cams), lat, lon, t0, t1, alt_asl, alt, RC.TARGET_NM)
                if bml is not None and bmt is not None:
                    c.update(bml=np.asarray(bml, "f4"), bmt=np.asarray(bmt, "f4"), mol_src="cams")
            except Exception as exc:  # noqa: BLE001
                c["err"] = f"{type(exc).__name__}: {exc}"[:160]
        if c.get("bml") is None:
            c.update(bml=np.asarray(IC._molecular_beta(z, alt, lam), "f4"),
                     bmt=np.asarray(IC._molecular_beta(z, alt, RC.TARGET_NM), "f4"),
                     mol_src="us1976")
        corr[ident] = c
    return dict(hours=hours, inst=inst, pairs=pairs, corr=corr)


# ------------------------------------------------------------------------------------- one night
def process_night(site_key, d, cen):
    """Read every instrument of the site once for night d and compute everything."""
    from calibration.io.instrument_day import load_instrument_day
    site = RC.SITES[site_key]
    wmo = site["wmo"]
    ds = RC.ymd(d)
    rec = dict(site=site_key, date=ds, ray={}, cloud={}, ing=None, t_read={}, err={})
    loaded = {}
    meta = None
    for ident, (itype, _methods) in site["streams"].items():
        fD = RC.l1_file(wmo, ident, d)
        if not fD.exists():
            continue
        files = [str(f) for f in (RC.l1_file(wmo, ident, d - timedelta(days=1)), fD) if f.exists()]
        t0 = time.perf_counter()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                idd = load_instrument_day(files, itype, RC.CAMS, read_cams=False,
                                          build_working=False)
        except Exception as exc:  # noqa: BLE001
            rec["err"][ident] = f"load: {type(exc).__name__}: {exc}"[:200]
            idd = None
        rec["t_read"][ident] = round(time.perf_counter() - t0, 3)
        if idd is None:
            continue
        loaded[ident] = (itype, idd)
        if meta is None:
            s = cen[(wmo, ident)]
            meta = dict(lat=float(s["lat"]), lon=float(s["lon"]), alt=float(idd.native.altitude))
    if not loaded:
        return rec
    # ingredients first: they only read the in-memory arrays
    t0 = time.perf_counter()
    try:
        rec["ing"] = ingredients(loaded, d, site_key, meta)
    except Exception as exc:  # noqa: BLE001
        rec["err"]["ing"] = f"{type(exc).__name__}: {exc}"[:300] + " | " + \
            traceback.format_exc()[-600:]
    rec["t_ing"] = round(time.perf_counter() - t0, 3)
    for ident, (itype, idd) in loaded.items():
        s = cen[(wmo, ident)]
        methods = site["streams"][ident][1]
        if "ray" in methods:
            rec["ray"][ident] = {v: rayleigh_night(idd.native, ds, s, v) for v in VARIANTS}
        if "cloud" in methods:
            rec["cloud"][ident] = cloud_day(idd, d, s, RC.l1_file(wmo, ident, d))
    return rec


def run_task(site_key, month, dates, out_path):
    if out_path.exists():
        return str(out_path), 0, 0.0
    warnings.filterwarnings("ignore")
    logging.disable(logging.CRITICAL)
    cen = RC.census()
    t0 = time.time()
    recs = []
    for ds in dates:
        d = datetime.strptime(ds, "%Y%m%d").date()
        try:
            recs.append(process_night(site_key, d, cen))
        except Exception as exc:  # noqa: BLE001
            recs.append(dict(site=site_key, date=ds, fatal=f"{type(exc).__name__}: {exc}",
                             tb=traceback.format_exc()[-1500:]))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tmp")
    with open(tmp, "wb") as f:
        pickle.dump(recs, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, out_path)
    return str(out_path), len(recs), time.time() - t0


# ------------------------------------------------------------------------------------- CLI
def cmd_run(args):
    sites = args.sites.split(",") if args.sites else list(RC.SITES)
    tasks = []
    for sk in sites:
        nights = RC.site_nights(sk)
        months = {}
        for d in nights:
            months.setdefault(d.strftime("%Y%m"), []).append(RC.ymd(d))
        for m, ds in sorted(months.items()):
            tasks.append((sk, m, ds, RC.RUN_DIR / sk / f"{m}.pkl"))
    todo = [t for t in tasks if not t[3].exists()]
    # longest first (the CL61 sites), so the pool does not end on one slow task
    weight = {"payerne": 3, "lindenberg": 3, "aosta": 3, "schiphol": 2, "budapest": 1, "magurele": 1}
    todo.sort(key=lambda t: -weight.get(t[0], 1) * len(t[2]))
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(tasks)} site-month tasks, {len(todo)} to run "
          f"({sum(len(t[2]) for t in todo)} site-nights), {args.workers} workers", flush=True)
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_task, *t): t for t in todo}
        for f in as_completed(futs):
            t = futs[f]
            try:
                p, n, dt = f.result()
            except Exception:  # noqa: BLE001
                print(f"FAILED {t[0]} {t[1]}", flush=True)
                traceback.print_exc()
                continue
            done += 1
            print(f"[{done}/{len(todo)}] {t[0]} {t[1]}: {n} nights in {dt:.0f} s "
                  f"(elapsed {time.time() - t0:.0f} s)", flush=True)


def cmd_test(args):
    logging.disable(logging.CRITICAL)
    warnings.filterwarnings("ignore")
    d = datetime.strptime(args.date, "%Y%m%d").date()
    t0 = time.time()
    rec = process_night(args.site, d, RC.census())
    print(f"wall {time.time() - t0:.1f} s, read {rec['t_read']}, ingredients {rec.get('t_ing')} s, "
          f"err {rec['err']}")
    for ident, vv in rec["ray"].items():
        for v, r in vv.items():
            print(ident, v, {k: (round(x, 5) if isinstance(x, float) else x) for k, x in r.items()
                             if not isinstance(x, np.ndarray)})
    for ident, r in rec["cloud"].items():
        print(ident, "cloud", {k: (round(x, 5) if isinstance(x, float) else x) for k, x in r.items()
                               if not isinstance(x, np.ndarray)})
    ing = rec["ing"]
    if ing:
        print("hours", ing["hours"].size, [str(np.datetime64(int(h), "h")) for h in ing["hours"]])
        for ident, v in ing["inst"].items():
            print(ident, v["itype"], "wl", v["wl"], "H finite", int(np.isfinite(v["H"]).sum()),
                  "adm300", int(v["adm300"].sum()), "excl", round(v["excl_frac"], 3))
        for k, v in ing["pairs"].items():
            print("pair", k, {w: int(a.sum()) for w, a in v.items()})
        for k, v in ing["corr"].items():
            print("corr", k, v.get("cams"), v.get("mol_src"),
                  None if v.get("t2") is None else float(np.nanmin(v["t2"])))
    if args.save:
        with open(args.save, "wb") as f:
            pickle.dump(rec, f)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=12)
    r.add_argument("--sites", default="")
    r.add_argument("--limit", type=int, default=0)
    t = sub.add_parser("test")
    t.add_argument("--site", default="payerne")
    t.add_argument("--date", default="20260715")
    t.add_argument("--save", default="")
    args = ap.parse_args()
    {"run": cmd_run, "test": cmd_test}[args.cmd](args)


if __name__ == "__main__":
    main()
