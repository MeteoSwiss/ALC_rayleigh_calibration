# -*- coding: utf-8 -*-
"""Inter-comparison dashboard v3 — compute stage (ALL sites into ONE payload).

WHY THIS EXISTS (and how it differs from build_l1_l2_dashboard.py)
-----------------------------------------------------------------
v2 precomputed EVERY (calibration variant x WV x wavelength) combination x every contiguous month
range, server-side.  That is why the Payerne page was 39 MB and why per-instrument freedom was
impossible: giving each of 4 channels its own method and variant is a combinatorial product of
hundreds of states, and each state carried a full copy of every profile, histogram and statistic.

v3 ships the INGREDIENTS instead, and recombines them in the browser.  The key algebraic facts
that make this exact rather than approximate:

  * the calibration divides:            beta = rcs_0 / C_L(t) * 1e6
  * the water-vapour correction is a per-(day, gate) DIVISION:   beta / T2(day, z)
  * the wavelength conversions are AFFINE in beta, with coefficients that depend only on
    (day, gate):   molecular  ->  beta_mol_1064(day,z) + (beta - beta_mol_lam(day,z)) * f
                   Angstrom   ->  beta * f                       f = (lam/1064)^alpha  (scalar)
  * the measured dark baseline is a per-gate CONSTANT subtracted from rcs_0.

So the whole chain is  out(h,z) = A(day(h),z) * (raw(h,z) - b(z)) / C(t_h) + B(day(h),z), and the
browser can evaluate it per hour for ANY per-instrument choice of method/variant/corrections.  The
payload therefore holds, per site:

  * ONE raw block per (source, instrument): the hourly-median, cloud-screened, regridded rcs_0
    (L1) / attenuated_backscatter_0 (L2) on the display grid, on the paired-hour axis.  No
    calibration, no corrections applied.  (log-quantised uint16 -> base64)
  * per instrument and per DAY: T2, beta_mol(lam), beta_mol(1064) on the same display grid, plus
    the scalar Angstrom factor.  (float32 -> base64)
  * per (instrument, method, variant): C_L interpolated onto the paired-hour axis, plus the
    per-night points and the Kalman line for the constants charts, plus an availability flag and
    the REASON when a variant does not exist for that channel.

Month-range filtering, the reference-channel choice, medians, IQR, histograms, the statistics
table and the noise-floor mask are all per-hour arithmetic and are done client-side.

Correction ordering note: v2 applied the corrections on the instrument's NATIVE range grid and then
regridded; v3 regrids first and applies the (day, gate) factors on the display grid.  Median and
regrid are both linear/affine-commuting operators and the factors are constant within a day, so the
only difference is the second-order term of interpolating a product instead of multiplying two
interpolants -- measured at < 0.02 % on the statistics band (see check_v3_regression.py).

Run:  python inter-comparison_dashboard/build_v3.py [site ...]
Out:  C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard/v3_<site>.json
"""
from __future__ import annotations
import base64
import json
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
sys.path.insert(0, str(Path(__file__).resolve().parent))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
os.environ.setdefault("ALC_VAL_CAMS_04", "A:/CAMS_Monthly_04")

import numpy as np

from validation.paper import intercompare as IC
from validation.paper import run_paper_validation as RP
from monitoring.kalman import kalman_best_estimate
import build_l1_l2_dashboard as BD          # cache reader, dark_for, pwv, l2_applied
import l1_l2_calib as CAL
import l1_l2_io as IO
import sites
import variants_v3 as V3

TARGET, ALPHA = 1064.0, 1.0
ZMIN, ZMAX = RP.ZMIN, RP.ZMAX               # 500-3000 m AGL statistics band
OUT = sites.DATA_ROOT

# ---------------------------------------------------------------------------- display grid
# 30 m up to 3 km so the 500-3000 m statistics band is BIT-IDENTICAL to the v2 uniform grid (the
# regrid half-width is the MEDIAN step, which stays 30 m), then progressively coarser -- above the
# band the panels are a display, and the noise-floor mask truncates most curves well below 8 km.
Z_AGL = np.concatenate([
    np.arange(0.0, 3000.0 + 30.0, 30.0),        # 101 gates, 30 m   (band = gates 17..100)
    np.arange(3090.0, 8010.0 + 90.0, 90.0),     #  55 gates, 90 m
    np.arange(8370.0, 15090.0, 360.0),          #  19 gates, 360 m
])
BAND = (Z_AGL >= ZMIN) & (Z_AGL <= ZMAX)
# Curtain payload: every 2nd display gate up to 8 km and every 2nd hour (the v2 page kept every
# hour, which alone was 10 MB of the 39 MB).
PC_ZSEL = np.where(Z_AGL <= 8000.0)[0][::2]
PC_HSTEP = 2


# ---------------------------------------------------------------------------- encoders
def b64(arr, dtype):
    return base64.b64encode(np.ascontiguousarray(arr, dtype=dtype).tobytes()).decode("ascii")


def quantize(A):
    """(n_h x n_z) float -> log-quantised uint16, per gate.

    Relative precision is what the ratio statistics care about, and one gate spans decades between
    a clear molecular night and a hazy hour, so a LINEAR quantiser would be coarse exactly where
    the signal is small.  Log spacing gives a uniform ~0.007 % relative step over the whole range.
    Codes: 0 = NaN, 32768 = exact zero, 1..32767 = positive magnitude, 32769..65535 = negative.
    """
    n_h, n_z = A.shape
    Q = np.zeros((n_h, n_z), "u2")
    lo = np.ones(n_z, "f8")
    hi = np.ones(n_z, "f8")
    for g in range(n_z):
        v = A[:, g]
        fin = np.isfinite(v)
        mag = np.abs(v[fin])
        pos = mag[mag > 0]
        if pos.size:
            h = float(pos.max())
            l = max(float(pos.min()), h * 1e-9)
            if not (h > l):
                l = h * (1 - 1e-9)
            lo[g], hi[g] = l, h
            with np.errstate(all="ignore"):
                mm = np.clip(np.nan_to_num(np.abs(v), nan=l), l, h)
                q = np.round(np.log(mm / l) / np.log(h / l) * 32766.0).astype("i8") + 1
            q = np.clip(q, 1, 32767)
            col = np.zeros(n_h, "u2")
            nzm = fin & (np.abs(v) > 0)
            col[nzm] = (q[nzm] + np.where(v[nzm] < 0, 32768, 0)).astype("u2")
            col[fin & (np.abs(v) == 0)] = 32768
            Q[:, g] = col
        else:
            Q[fin, g] = 32768                          # all-zero gate
    return dict(b=b64(Q, "<u2"), lo=[float(x) for x in lo], hi=[float(x) for x in hi],
                nh=int(n_h), nz=int(n_z))


def f32block(A):
    """(n x m) float -> base64 float32 (NaN survives the round trip)."""
    A = np.asarray(A, "f4")
    return dict(b=b64(A, "<f4"), n=int(A.shape[0]), m=int(A.shape[1]))


# ---------------------------------------------------------------------------- panels
def stream_block(d, ch_alt_grid):
    """One cached stream -> (hours, screened block, qf-only block) on the display grid.

    Mirrors build_l1_l2_dashboard.build_panel with the calibration and the corrections REMOVED --
    what is left is exactly the raw signal the browser will transform.
    """
    beta = np.asarray(d["beta"], "f8")
    scr = IC.screen(beta, d)
    disp = beta.copy()
    disp[d["qf"] > 0] = np.nan
    g, (S, Dq) = IC.retime_hourly(d["time"], [scr, disp], min_cov_s=IC.MIN_AVG_S, snr_idx=())
    alt = np.asarray(d["alt"], "f8")
    return g, IC.regrid(S, alt, ch_alt_grid), IC.regrid(Dq, alt, ch_alt_grid)


def corr_arrays(d, station_alt, alt_asl, days):
    """Per-day correction ingredients for ONE stream, evaluated on the DISPLAY grid.

    Returns dict(lam, f, wv=<n_days x n_z or None>, bml, bmt) where
      wv  = the two-way water-vapour transmission T2 (NaN for a day with no usable CAMS, which is
            how the v2 pipeline excludes those days),
      bml/bmt = molecular attenuated backscatter at the instrument's and the target wavelength
            (CAMS T/p, US-standard fallback), i.e. the affine offset of the "advanced" conversion.
    """
    lam = float(d["wavelength"]) if np.isfinite(d.get("wavelength", np.nan)) else \
        IC.WV_PARAMS.get(d["itype"], (np.nan,))[0]
    out = dict(lam=(None if not np.isfinite(lam) else float(lam)),
               f=float((lam / TARGET) ** ALPHA) if np.isfinite(lam) else 1.0,
               wv=None, wv2=None, bml=None, bmt=None)
    lat, lon = float(d["lat"]), float(d["lon"])
    z_agl = alt_asl - station_alt
    nz = alt_asl.size

    if np.isfinite(lam) and IC.in_water_vapor_band(lam):
        from calibration.water_vapor_correction.water_vapor import (
            two_way_wv_transmission, cams_point_too_far)
        # Two laser-line hypotheses: the pipeline's operational (lambda0, FWHM), and -- for the
        # CL61 only -- the manufacturer's measured line (910.55 nm, sigma 0.08 -> FWHM 0.188),
        # which is ~9x less absorbed.  Instruments the finding does not concern keep one array and
        # the page's "constructeur" option falls back to it.
        specs = {"nom": IC.WV_PARAMS.get(d["itype"], (910.0, 3.4))}
        if d["itype"] in V3.LAM_TYPES:
            specs["ctor"] = (910.55, 0.188)
        W = {k: np.full((len(days), nz), np.nan) for k in specs}
        for i, ds in enumerate(days):
            cams = IO.cams_for_day(ds)
            if cams is None or cams_point_too_far(cams, lat, lon):
                continue
            t0, t1 = IO._win(ds)
            prof = IC._cams_wv_profile_cached(str(cams), lat, lon, t0, t1)
            if prof is None:
                continue
            for k, (lam0, fwhm) in specs.items():
                t2 = np.asarray(two_way_wv_transmission(alt_asl, station_alt, prof[0], prof[1],
                                                        IC.WV_LUT, lam0, fwhm), "f8")
                if t2.size == nz and np.any(np.isfinite(t2)):
                    W[k][i] = t2
        out["wv"] = W["nom"]
        out["wv2"] = W.get("ctor")

    if np.isfinite(lam) and abs(lam - TARGET) >= 1.0:
        Bl = np.full((len(days), nz), np.nan)
        Bt = np.full((len(days), nz), np.nan)
        us_l = us_t = None
        for i, ds in enumerate(days):
            cams = IO.cams_for_day(ds)
            bml = bmt = None
            if cams is not None:
                t0, t1 = IO._win(ds)
                bml = IC._mol_att_cams(str(cams), lat, lon, t0, t1, alt_asl, station_alt, lam)
                bmt = IC._mol_att_cams(str(cams), lat, lon, t0, t1, alt_asl, station_alt, TARGET)
            if bml is None or bmt is None:
                if us_l is None:
                    us_l = IC._molecular_beta(z_agl, station_alt, lam)
                    us_t = IC._molecular_beta(z_agl, station_alt, TARGET)
                bml, bmt = us_l, us_t
            Bl[i], Bt[i] = bml, bmt
        out["bml"], out["bmt"] = Bl, Bt
    return out


# ---------------------------------------------------------------------------- calibration series
def _kalman_record(dates, C, U, label_method):
    kt, ks, kstd = kalman_best_estimate(dates, C, uncertainties=U)
    if not len(kt):
        return None
    return dict(
        nights=int(len(C)),
        points=[dict(d=d.strftime("%Y-%m-%d"), v=float(f"{c:.6g}"),
                     s=(float(f"{s:.4g}") if np.isfinite(s) else None))
                for d, c, s in zip(dates, C, U)],
        method=label_method,
        kal=dict(d=[str(t)[:10] for t in kt], v=[float(f"{v:.6g}") for v in ks],
                 s=[float(f"{s:.4g}") if np.isfinite(s) else 0.0 for s in kstd]),
        _kt=kt, _ks=ks)


def calib_for(site_key, ident, itype, method, variant, wmo):
    """(record | None, reason) for one (instrument, method, variant) triple.

    Never carries a series over from another variant: a channel a run does not cover is reported as
    UNAVAILABLE with the reason, so the page can grey it out instead of showing someone else's
    numbers under this variant's name.
    """
    spec = V3.source_for(site_key, ident, itype, method, variant, wmo)
    if isinstance(spec, str):
        return None, spec
    kind, arg = spec
    if kind == "nc":                                  # v2.0 published NetCDFs (Payerne only)
        try:
            dates, C, U, M, created, src = CAL.all_nights(ident)
        except FileNotFoundError:
            return None, f"aucun NetCDF ALC_calibration_{wmo}_{ident}<année>.nc sur ce poste"
        keep = np.isin(M, np.asarray([0 if method == "rayleigh" else 1]))
        if keep.sum() < 5:
            return None, (f"seulement {int(keep.sum())} nuit(s) « {method} » dans les NetCDF v2.0")
        dates = [d for d, k in zip(dates, keep) if k]
        C, U = C[keep], U[keep]
    elif kind == "json":                              # per-night json from the availability corpus
        nights = CAL.v22_nights(arg)
        if nights is None:
            return None, f"pas de sortie eprof_v2.2 par nuit ({arg}) pour ce canal"
        dates, C, U = nights
    else:                                             # kind == "csv": a network-runner tree
        nights = CAL.run_csv_nights(arg, method)
        if nights is None:
            return None, (f"aucune nuit « {method} » exploitable dans "
                          f"{Path(arg).parent.parent.name}/{Path(arg).parent.name}")
        dates, C, U = nights
    rec = _kalman_record(dates, C, U, "Rayleigh" if method == "rayleigh" else "Nuages liquides")
    if rec is None:
        return None, f"trop peu de nuits ({len(C)}) pour le filtre de Kalman"
    return rec, None


# ---------------------------------------------------------------------------- one site
def build_site(site_key):
    s = sites.get_site(site_key)
    v3 = V3.SITE_V3[site_key]
    BD.set_site(site_key)
    npz = str(BD.CACHE_NPZ)
    if not Path(npz).exists():
        raise SystemExit(f"stream cache {npz} missing — run build_l1_l2_dashboard.py {site_key} first")
    t0, t1 = BD.common_window(npz)
    cache = BD.load_cache(npz, t0, t1)
    sources = tuple(s["sources"])
    idents = list(dict.fromkeys(i["ident"] for i in v3["instruments"]))
    station_alt = float(s["alt"])
    alt_asl = Z_AGL + station_alt
    print(f"== {s['name']} ==  window {np.datetime64(t0,'D')} .. {np.datetime64(t1,'D')}  "
          f"streams {len(sources)}x{len(idents)}", flush=True)

    # --- raw panels ---------------------------------------------------------------------------
    raw, hours_of = {}, {}
    for src in sources:
        for ident in idents:
            g, S, Dq = stream_block(cache[(src, ident)], alt_asl)
            raw[(src, ident)] = (S, Dq)
            hours_of[(src, ident)] = g
    union = np.unique(np.concatenate(list(hours_of.values())))
    aligned, aligned_disp = {}, {}
    for k, (S, Dq) in raw.items():
        pos = np.searchsorted(union, hours_of[k])
        A = np.full((union.size, Z_AGL.size), np.nan)
        A[pos] = S
        aligned[k] = A
        if k[0] == "L1":
            Ad = np.full((union.size, Z_AGL.size), np.nan)
            Ad[pos] = Dq
            aligned_disp[k] = Ad

    # Paired hours: every instrument, in every source, has at least one finite gate in the band.
    # Identical to the v2 `have` mask (the corrections never turn a finite value non-finite except
    # on a whole CAMS-less day, which the page reports through the same NaN path).
    have = np.logical_and.reduce([np.any(np.isfinite(aligned[k][:, BAND]), axis=1)
                                  for k in aligned])
    hours = union[have]
    print(f"   paired hours: {hours.size}", flush=True)

    # Day axis over the UNPAIRED union: the curtains use the full axis, and both must index the
    # same per-day correction arrays.
    days = [str(x).replace("-", "") for x in np.unique(union.astype("datetime64[D]"))]
    day_index = {d: i for i, d in enumerate(days)}
    hour_day = [day_index[str(x).replace("-", "")] for x in hours.astype("datetime64[D]")]
    union_day = np.array([day_index[str(x).replace("-", "")]
                          for x in union.astype("datetime64[D]")])
    months = [str(m) for m in np.unique(hours.astype("datetime64[M]"))]
    month_index = {m: i for i, m in enumerate(months)}
    hour_month = [month_index[str(m)] for m in hours.astype("datetime64[M]")]

    # --- correction ingredients (per instrument; L1/L2 share them when the wavelength agrees) ---
    corr, corr_of = {}, {}
    for ident in idents:
        for src in sources:
            d = cache[(src, ident)]
            lam = float(d["wavelength"]) if np.isfinite(d.get("wavelength", np.nan)) else np.nan
            same = next((k for k in corr if k.split("_")[0] == ident
                         and ((corr[k]["lam"] is None and not np.isfinite(lam))
                              or (corr[k]["lam"] is not None and np.isfinite(lam)
                                  and abs(corr[k]["lam"] - lam) < 1e-6))), None)
            if same is None:
                same = ident if ident not in corr else f"{ident}_{src}"
                corr[same] = corr_arrays(d, station_alt, alt_asl, days)
                print(f"   corr {same}: lam={corr[same]['lam']} "
                      f"wv={'yes' if corr[same]['wv'] is not None else 'no'} "
                      f"wl={'yes' if corr[same]['bml'] is not None else 'no'}", flush=True)
            corr_of[(src, ident)] = same

    # --- calibration series -------------------------------------------------------------------
    calib = {}
    for inst in v3["instruments"]:
        for method in V3.METHOD_BY_TYPE.get(inst["itype"], ["rayleigh"]):
            for variant in V3.variants_for(inst["itype"], method):
                rec, why = calib_for(site_key, inst["ident"], inst["itype"], method, variant,
                                     s["wmo"])
                ck = f'{inst["ident"]}|{method}|{variant}'
                if rec is None:
                    calib[ck] = dict(ok=False, why=why)
                    continue
                c = IC.interp_calib(rec.pop("_kt"), rec.pop("_ks"), hours)
                rec["c"] = [float(f"{v:.7g}") for v in c]
                rec["ok"] = True
                calib[ck] = rec
                print(f"   {ck:34s} {rec['nights']:4d} nights  C_L[0]={c[0]:.4g}", flush=True)

    # Warning prose: interpolate the numbers from the series just computed (see prose_tokens).
    tok = prose_tokens(site_key, calib, t0, t1)
    warns = [w.format_map(tok) for w in v3.get("warnings", [])]

    # --- dark baselines (measured; Payerne only) -----------------------------------------------
    dark = {}
    for ident in idents:
        b = BD.dark_for(ident, Z_AGL) if s.get("dark_npz") else None
        dark[ident] = None if b is None else [float(f"{x:.6g}") for x in b]

    # --- payload blocks -------------------------------------------------------------------------
    streams = {}
    for (src, ident), A in aligned.items():
        streams[f"{src}|{ident}"] = quantize(A[have])

    corr_out = {}
    for k, c in corr.items():
        corr_out[k] = dict(lam=c["lam"], f=c["f"],
                           wv=(None if c["wv"] is None else f32block(c["wv"])),
                           wv2=(None if c.get("wv2") is None else f32block(c["wv2"])),
                           bml=(None if c["bml"] is None else f32block(c["bml"])),
                           bmt=(None if c["bmt"] is None else f32block(c["bmt"])))

    payload = dict(
        key=site_key, name=s["name"], wmo=s["wmo"], lat=s["lat"], lon=s["lon"], alt=station_alt,
        start=str(np.datetime64(t0, "D")), end=str(np.datetime64(t1, "D")),
        sources=list(sources), zmin=ZMIN, zmax=ZMAX, target=TARGET, alpha=ALPHA,
        z=[float(x) for x in Z_AGL], band=[int(np.where(BAND)[0][0]), int(np.where(BAND)[0][-1])],
        hours=[str(t)[:13] for t in hours.astype("datetime64[h]")],
        hour_day=hour_day, hour_month=hour_month, days=days, months=months,
        instruments=v3["instruments"], default=v3["default"],
        default_dark=v3.get("default_dark", {}), iref=v3.get("iref", 0),
        title=v3["title"], subtitle=v3["subtitle"], warnings=warns,
        streams=streams, corr=corr_out, corr_of={f"{a}|{b}": v for (a, b), v in corr_of.items()},
        dark=dark, dark_kind=V3.dark_kind(site_key), calib=calib,
        l2_applied=BD.l2_applied_constants(npz, t0, t1),
    )
    payload["hopkin"] = hopkin_payload(v3, s["wmo"], calib)
    payload["pcolor"] = pcolor_payload(v3, aligned, aligned_disp, union, union_day, corr, corr_of,
                                       calib, dark)
    if v3.get("pwv"):
        payload["pwv"] = pwv_payload(v3, aligned, have, hours, hour_day, corr, corr_of, calib,
                                     dark)
    out = OUT / f"v3_{site_key}.json"
    out.write_text(json.dumps(payload), encoding="utf-8")
    print(f"-> {out.name}  ({out.stat().st_size/1e6:.2f} MB)", flush=True)
    return payload


# ---------------------------------------------------------------------------- transform (server)
def transform(rawblock, cval, darkb, c, wv, wl, hdayidx):
    """The EXACT arithmetic the browser performs, in numpy — used for the curtains, the PWV panel
    and the non-regression check, so the two implementations are provably the same formula."""
    v = np.asarray(rawblock, "f8").copy()
    if darkb is not None:
        v = v - np.asarray(darkb, "f8")[None, :]
    if cval is not None:
        v = v * 1e6 / np.asarray(cval, "f8")[:, None]
    wv = {True: "nom", False: "none"}.get(wv, wv)
    if wv != "none" and c.get("wv") is not None:
        W = c["wv2"] if (wv == "ctor" and c.get("wv2") is not None) else c["wv"]
        v = v / np.asarray(W, "f8")[hdayidx]
    need_wl = c.get("lam") is not None and abs(c["lam"] - TARGET) >= 1.0
    if need_wl and wl == "molecular" and c.get("bmt") is not None:
        v = np.asarray(c["bmt"], "f8")[hdayidx] + (v - np.asarray(c["bml"], "f8")[hdayidx]) * c["f"]
    elif need_wl and wl == "angstrom":
        v = v * c["f"]
    return v


# ---------------------------------------------------------------------------- curtains
def pcolor_payload(v3, aligned, aligned_disp, union, union_day, corr, corr_of, calib, dark):
    """Static time-height curtains, one per instrument, in the page's DEFAULT state.

    Built from the UNPAIRED L1 axis (the strict pairing drops every hour any instrument is
    cloud-screened, which would riddle the curtain with holes), every 2nd hour and every 2nd
    display gate up to 8 km, log10 quantised to one byte (0.014 decade per step — invisible on a
    Viridis scale, and one quarter of the v2 payload).
    """
    h = union.astype("datetime64[h]")
    full = np.arange(h.min(), h.max() + np.timedelta64(1, "h"), np.timedelta64(1, "h"))
    pos = np.searchsorted(full, h)
    axis = full[::PC_HSTEP]
    out = []
    for inst in v3["instruments"]:
        ident = inst["ident"]
        method, variant = v3["default"][ident]
        rec = calib.get(f"{ident}|{method}|{variant}")
        if rec and rec.get("ok"):        # constants re-interpolated onto the UNPAIRED axis
            kd = np.array([np.datetime64(x) for x in rec["kal"]["d"]])
            cval = IC.interp_calib(kd, np.asarray(rec["kal"]["v"], "f8"), union)
        else:
            cval = np.ones(union.size)
        c = corr[corr_of[("L1", ident)]]
        # the boot state's profile-side dark, so the curtain matches what the page shows on load
        db = dark.get(ident) if v3.get("default_dark", {}).get(ident) else None
        wvmode = "ctor" if variant in getattr(V3, "REFERENCE_VARIANTS", ()) else True
        S = transform(aligned[("L1", ident)], cval, db, c, wvmode, "molecular", union_day)
        D = transform(aligned_disp[("L1", ident)], cval, db, c, wvmode, "molecular", union_day)
        med = np.nanmedian(S, axis=0)
        n = np.sum(np.isfinite(S), axis=0)
        nprof = int(np.any(np.isfinite(S), axis=1).sum())
        kp = BD._keep_mask(med, n, Z_AGL, max(nprof, 1))
        S[:, ~kp] = np.nan
        D[:, ~kp] = np.nan
        D[np.isfinite(S)] = np.nan
        blocks = {}
        for name, M in (("scr", S), ("disp", D)):
            G = np.full((full.size, PC_ZSEL.size), np.nan)
            G[pos] = M[:, PC_ZSEL]
            with np.errstate(all="ignore"):
                L = np.log10(G[::PC_HSTEP])
            q = np.full(L.shape, 255, "u1")
            ok = np.isfinite(L)
            q[ok] = np.clip(np.round((np.clip(L[ok], -2.0, 1.5) + 2.0) / 3.5 * 254.0), 0,
                            254).astype("u1")
            blocks[name] = b64(q.T, "u1")            # transposed: plotly wants z[y][x]
        out.append(dict(ident=ident, **blocks))
    return dict(hours=[str(t) for t in axis], alt=[float(Z_AGL[i]) for i in PC_ZSEL],
                nx=int(axis.size), ny=int(PC_ZSEL.size), lo=-2.0, hi=1.5, ch=out)


# ---------------------------------------------------------------------------- PWV panel
def pwv_payload(v3, aligned, have, hours, hour_day, corr, corr_of, calib, dark):
    """Per-hour CL61-Rayleigh residual vs the reference, against the day's CAMS PWV, for the four
    (constants variant x WV-comparison) states — the WV design-mitigation test.  Unchanged in
    substance from v2; only the arithmetic route is new."""
    spec = v3["pwv"]
    ident, ref_ident = spec["ident"], spec["ref"]
    hd = np.asarray(hour_day)
    band = np.where(BAND)[0]
    ref_m, ref_v = v3["default"][ref_ident]
    ref_rec = calib.get(f"{ref_ident}|{ref_m}|{ref_v}")
    if not (ref_rec and ref_rec.get("ok")):
        return None
    Aref = aligned[("L1", ref_ident)][have]
    Acur = aligned[("L1", ident)][have]
    cref = np.asarray(ref_rec["c"], "f8")
    # COHERENCE: when the site default reference is a dark-corrected run, its measured b(z) must
    # be subtracted from the reference profile too — dark-run constants over an un-darked profile
    # is exactly the hybrid half-state the page's DARK_TWIN machinery forbids (found by the
    # 2026-08-16 adversarial review: the hybrid inflated every displayed slope by +0.4..+0.8 %/mm).
    ref_db = dark.get(ref_ident) if ref_v in getattr(V3, "DARK_RUNS", ()) else None
    parts, fits = {}, {}
    for st in spec["states"]:
        variant, wv = st[0], st[1]
        rec = calib.get(f"{ident}|{spec.get('method','rayleigh')}|{variant}")
        if not (rec and rec.get("ok")):
            continue
        ccur = np.asarray(rec["c"], "f8")
        R = transform(Aref, cref, ref_db, corr[corr_of[("L1", ref_ident)]], wv, "molecular", hd)
        C = transform(Acur, ccur, None, corr[corr_of[("L1", ident)]], wv, "molecular", hd)
        a, b = C[:, band], R[:, band]
        m = np.isfinite(a) & np.isfinite(b) & (b > 0)
        rel = np.where(m, (a - b) / np.where(m, b, 1.0) * 100.0, np.nan)
        parts[f"{variant}|{wv}"] = np.nanmedian(rel, axis=1)
    x = np.array([BD.pwv_mm_for_day(str(h)[:10].replace("-", "")) for h in hours])
    resid = {}
    for st, y in parts.items():
        resid[st] = [None if not np.isfinite(v) else round(float(v), 3) for v in y]
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() >= 3:
            # Theil-Sen, the SAME statistic the page displays (an OLS here used to log the
            # opposite state ranking to the published panel — outlier leverage)
            xm, ym = x[m], y[m]
            dx = xm[None, :] - xm[:, None]
            iu = np.triu_indices(xm.size, 1)
            pair = (ym[None, :] - ym[:, None])[iu][dx[iu] != 0] / dx[iu][dx[iu] != 0]
            slope = float(np.median(pair))
            icept = float(np.median(ym - slope * xm))
            fits[st] = dict(slope=round(slope, 3), intercept=round(icept, 3), n=int(m.sum()))
            print(f"   pwv {st:<18s} n={int(m.sum()):4d}  TS slope {slope:+.3f} %/mm", flush=True)
    return dict(pwv_mm=[None if not np.isfinite(v) else round(float(v), 3) for v in x],
                resid=resid, fits=fits, ident=ident, ref=ref_ident,
                band=[float(ZMIN), float(ZMAX)], states=spec["states"])


# ---------------------------------------------------------------------------- prose tokens
class _SafeTok(dict):
    """format_map dict that leaves unknown {tokens} verbatim, so a warning without tokens (or a
    token another site does not define) can never crash the build."""
    def __missing__(self, key):
        return "{" + key + "}"


def _fmt_c(v):
    """4.12e7 -> '4.1e7' (the compact scientific style the prose uses)."""
    from math import floor, log10
    if not np.isfinite(v) or v <= 0:
        return "—"
    e = int(floor(log10(v)))
    return f"{v / 10**e:.1f}e{e}"


def prose_tokens(site_key, calib, t0, t1):
    """Numbers the warning prose interpolates, computed from the SAME series the page uses.

    Any count or level that describes the current window lives here, not hard-coded in
    variants_v3 — a re-run with new data re-derives the prose automatically (the stale
    'une seule nuit Rayleigh' claim of 2026-08 was exactly this class of rot).
    """
    tok = {}

    def nights(key, d0=None, d1=None):
        rec = calib.get(key)
        if not (rec and rec.get("ok")):
            return "—"
        ds = [p["d"] for p in rec["points"]]
        if d0:
            ds = [d for d in ds if d0 <= d <= d1]
        return len(ds)

    w0, w1 = str(np.datetime64(t0, "D")), str(np.datetime64(t1, "D"))
    if site_key == "payerne":
        for v, name in (("v2.0", "V20"), ("v2.2", "V22"), ("v2.2dark", "DARK")):
            tok[f"A_{name}_WIN"] = nights(f"A|rayleigh|{v}", w0, w1)
            tok[f"A_{name}_JA"] = nights(f"A|rayleigh|{v}", "2026-07-08", "2026-08-13")
        rec = calib.get("B|cloud|cloudWV")
        if rec and rec.get("ok"):
            pts = [(p["d"], p["v"]) for p in rec["points"]]
            for name, d0, d1 in (("SPRING", "2026-01-01", "2026-05-31"),
                                 ("JUNE", "2026-06-01", "2026-06-30"),
                                 ("SUMMER", "2026-07-08", "2026-08-31")):
                v = [x for d, x in pts if d0 <= d <= d1]
                tok[f"B_LVL_{name}"] = _fmt_c(float(np.median(v))) if v else "—"
    if site_key == "lindenberg":
        for key, name in (("C|cloud|cloudWV", "LIN_V22_END"),
                          ("C|cloud|cloud_l55s008", "LIN_RERUN_END")):
            rec = calib.get(key)
            tok[name] = rec["points"][-1]["d"] if rec and rec.get("ok") else "—"
    return _SafeTok(tok)


# ---------------------------------------------------------------------------- Hopkin heatmap
# Per-PROFILE view of the cloud constant against cloud-base height (Hopkin et al. 2019, AMT, Fig. 6
# conventions, as in rayleigh_availability/cloud_profile_dump.py): x = the profile's 1/C as a % of
# the unit median, y = CBH.  A tilted cloud means the retrieved constant depends on how far the
# beam travelled through the atmosphere -- exactly what an over-strong water-vapour correction
# produces, so the panel is the per-profile twin of the WV-spectrum question.
HK_X = np.arange(60.0, 140.0 + 1e-9, 2.0)
HK_Y = np.arange(0.25, 2.5 + 1e-9, 0.1)
HK_BAND = 0.25


def _cluster_ols(x, y, cluster):
    """OLS slope with CR1 cluster-robust SE, clusters = calendar days (profiles inside a day share
    the weather and one calibration state).  Same estimator as cloud_profile_dump.cluster_ols."""
    n = x.size
    if n < 5:
        return float("nan"), float("nan"), 0
    X = np.column_stack([np.ones(n), x])
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    e = y - X @ beta
    meat = np.zeros((2, 2))
    groups = np.unique(cluster)
    for g in groups:
        m = cluster == g
        u = X[m].T @ e[m]
        meat += np.outer(u, u)
    c = (groups.size / max(groups.size - 1, 1)) * ((n - 1) / max(n - 2, 1))
    V = c * XtX_inv @ meat @ XtX_inv
    return float(beta[1]), float(np.sqrt(V[1, 1])), int(groups.size)


def hopkin_payload(v3, wmo, calib):
    """{ "<ident>|<cloud variant>": {hist, bands, slope, ...} } for every dump on disk.

    A variant whose CONSTANTS are unavailable for this channel is skipped even when a dump folder
    exists: the CL31 was not part of the laser-line rerun, so its "l55*" dumps are byte-identical
    copies of the nominal one and showing them under a lambda label would be a lie.
    """
    out = {}
    for inst in v3["instruments"]:
        if "cloud" not in V3.METHOD_BY_TYPE.get(inst["itype"], []):
            continue
        for variant in V3.variants_for(inst["itype"], "cloud"):
            if not (calib.get(f'{inst["ident"]}|cloud|{variant}') or {}).get("ok"):
                continue
            f = V3.dump_npz(variant, wmo, inst["ident"])
            if f is None or not Path(f).exists():
                continue
            with np.load(f, allow_pickle=True) as z:
                cbh = np.asarray(z["cbh_m"], "f8") / 1000.0
                cinv = 1.0 / np.asarray(z["c_oconnor"], "f8")
                day = np.asarray(z["date"])
            ok = np.isfinite(cbh) & np.isfinite(cinv) & (cinv > 0)
            cbh, cinv, day = cbh[ok], cinv[ok], day[ok]
            if cbh.size < 50:
                continue
            pct = 100.0 * cinv / np.median(cinv)
            lnc = 100.0 * np.log(cinv)
            slope, se, ndays = _cluster_ols(cbh, lnc, day)
            inwin = (pct >= HK_X[0]) & (pct <= HK_X[-1])
            slope_w, se_w, _ = _cluster_ols(cbh[inwin], lnc[inwin], day[inwin])
            H, _, _ = np.histogram2d(pct, cbh, bins=[HK_X, HK_Y])
            bands = []
            edges = np.arange(HK_Y[0], HK_Y[-1] + 1e-9, HK_BAND)
            for lo, hi in zip(edges[:-1], edges[1:]):
                v = pct[inwin & (cbh >= lo) & (cbh < hi)]
                if v.size >= 5:
                    bands.append([round(float((lo + hi) / 2), 3), round(float(v.mean()), 2),
                                  round(float(v.std()), 2), int(v.size)])
            out[f'{inst["ident"]}|{variant}'] = dict(
                hist=[[int(x) for x in H[:, j]] for j in range(H.shape[1])],   # z[y][x]
                bands=bands, n=int(cbh.size), ndays=ndays,
                slope=round(slope, 2), se=round(se, 2),
                slope_win=round(slope_w, 2), se_win=round(se_w, 2),
                frac_out=round(float(1.0 - inwin.mean()), 4))
            print(f"   hopkin {inst['ident']}|{variant:14s} n={cbh.size:6d} "
                  f"slope {slope:+6.2f}±{se:.2f} %/km", flush=True)
    return dict(x=[float(v) for v in HK_X], y=[float(v) for v in HK_Y], panels=out)


def main():
    keys = [a for a in sys.argv[1:] if not a.startswith("--")] or list(V3.SITE_V3)
    OUT.mkdir(parents=True, exist_ok=True)
    for k in keys:
        build_site(k)


if __name__ == "__main__":
    main()
