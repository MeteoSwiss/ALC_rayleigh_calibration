# -*- coding: utf-8 -*-
"""MECHANISM SCAN 3 -- can an OVERLAP error make the lidar constant altitude-dependent?

CONTEXT
-------
The observational study shows the retrieved Rayleigh constant depends on where the molecular
window sits: v2.2-recovered nights (which fit 1-2 km higher than v2.0-kept ones) give -22.5 % at
Payerne CHM15k, -17.5 % at Payerne CL61, +23.2 % at Lindenberg, +31.7 / +32.1 % at Aosta,
+7.4 % at SIRTA. "Bad overlap" is the first suspicion of any ceilometer operator, and this repo
has MEASURED overlap defects on exactly these instruments
(``doc/reports/08_overlap_nearrange_offset.md``): the generic reference table TUB120011 completes
at ~750 m while real CHM15k overlaps complete near 1300-1500 m, and the Payerne TUB140016
module's applied overlap is ~14 % off its true overlap at 400-500 m.

This module answers three questions with the SHIPPED retrieval:

 (a) does a residual overlap error at 0.5-1.5 km reach a constant fitted ABOVE 2 km at all?
 (b) does it produce a SLOPE dC_L/dz, or only an OFFSET?
 (c) what would the CHM15k-vs-CL61 overlap difference alone do to a co-located pair?

THE FORWARD EQUATION (identical to forward_model.py, one extra factor)
---------------------------------------------------------------------
Truth, single scattering, range z in m:

    rcs_true(z) = C_true * O_true(z) * beta_tot(z) * exp(-2 * INT_0^z ext_tot dz')      [1]
    beta_tot = beta_mol + beta_aer                          [m^-1 sr^-1]
    ext_tot  = beta_mol * 8*pi/3 + beta_aer * S_aer         [m^-1],  S_aer = 52 sr

The CHM15k firmware (and the L1 producer generally) divides by its OWN overlap table
``O_applied``, so what the calibration pipeline reads as ``rcs_0`` is

    rcs_L1(z) = rcs_true(z) / O_applied(z)
              = C_true * O_err(z) * beta_tot(z) * exp(-2*OD(z)),  O_err = O_true / O_applied  [2]

``O_err`` is the RESIDUAL overlap error, dimensionless, exactly 1 where the correction is right.
It is the only quantity that matters here, and it is passed to
``forward_model.forward_signal(..., overlap=O_err)``. Everything else is the closed forward model
whose molecular closure residual is 0.0046 % (``forward_model.run_closure_tests``).

WHY AN OVERLAP ERROR CAN REACH A HIGH WINDOW AT ALL
---------------------------------------------------
The window fit only sees gates above 2 km, where O_err = 1 for every ALC here, so the fitted
slope is untouched. But ``_compute_cl_for_perturbation`` runs the Klett inversion from
``i_start = 0`` (the GROUND), and ``calculate_lidar_constant`` integrates the optical depth from
the ground too:

    C_L(z) = rcs_L1(z) / beta_tot_klett(z) * exp(2 * INT_0^z ext_tot_klett dz')          [3]

A signal deficit at 500 m makes Klett return a beta_tot -- and hence an ext_aer -- that is too
low there, so the optical depth accumulated below the overlap-completion height z_full is wrong
by Delta_OD, and every C_L above z_full is multiplied by exp(2*Delta_OD).

That is the propagation path the task asks about, and it is a MULTIPLICATIVE CONSTANT. The Klett
backward solution below any height is fixed by the (correct) solution AT that height, so once the
reference window sits above z_full the whole retrieved profile below z_full -- and therefore
Delta_OD -- is identical no matter where the window is. The prediction is a pure offset; the
simulation below measures how pure.

HOW EVERY NUMBER IS ATTRIBUTED
------------------------------
Every scenario is run TWICE in the same atmosphere: once with ``O_err = 1`` (baseline) and once
with the residual error. Only the DIFFERENCE is attributed to overlap. This matters, because a
boundary-layer aerosol that reaches into the lowest fit windows produces a large slope of its own
(-4 %/km with the profile used here) which has nothing to do with overlap -- see
``aerosol_reaching_the_window`` in the JSON, and mechanism scan on aerosol.

OUTPUT
------
Figures in ``doc/reports/figs_altitude_audit/`` and every number in
``rayleigh_availability/overlap_scan_results.json``. Run::

    python rayleigh_availability/overlap_scan.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from calibration.rayleigh.atmosphere import MOLECULAR_LIDAR_RATIO      # noqa: E402
from rayleigh_availability.forward_model import (                      # noqa: E402
    PAYERNE_CHM15K,
    default_grid,
    default_options,
    forward_signal,
    make_atmosphere,
    retrieve,
)

FIGDIR = REPO / "doc" / "reports" / "figs_altitude_audit"
TUB_NPZ = REPO / "doc" / "reports" / "figs_cl61_chm_tilt" / "tub140016_overlap_noise.npz"
NETWORK_NPZ = REPO / "validation" / "overlap_clearsky" / "network_overlaps.npz"

# Ladder: forced molecular windows of +-490 m centred every 500 m from 2.5 to 6.5 km AGL.
CENTRES_M = np.arange(2500.0, 6501.0, 500.0)
HALF_WIDTH_M = 490.0

# Boundary-layer aerosol of the reference "realistic" night. It is deliberately confined BELOW
# the fit band (raised-cosine taper to exactly zero at 1900 m, first window gate 2010 m) so that
# the ladder baseline closes on C_true and the only thing left is the overlap error. The aerosol
# still fills the whole overlap-error region, which is what amplifies the defect.
BL_BETA_SURFACE = 1.2e-6      # m^-1 sr^-1 at the ground
BL_SCALE_HEIGHT = 700.0       # m
BL_TAPER_START = 1500.0       # m, where the raised-cosine cut begins
BL_TAPER_END = 1900.0         # m, where beta_aer is exactly 0


# ---------------------------------------------------------------------------
# Aerosol confined below the fit band
# ---------------------------------------------------------------------------
def capped_exponential_aerosol(z, beta_surface, scale_height_m, lidar_ratio_sr,
                               taper_start_m, taper_end_m):
    """``beta_aer = beta_surface * exp(-z/H)`` cut to exactly zero with a raised cosine.

    A hard truncation would put a step in beta_tot; the raised cosine over
    ``[taper_start_m, taper_end_m]`` leaves the profile C1-smooth and, critically, EXACTLY zero
    above ``taper_end_m`` so the molecular fit band is uncontaminated.
    Units: beta in m^-1 sr^-1, ext = beta * S_aer in m^-1.
    """
    z = np.asarray(z, float)
    beta = float(beta_surface) * np.exp(-z / float(scale_height_m))
    w = np.ones_like(z)
    ramp = (z > taper_start_m) & (z < taper_end_m)
    w[ramp] = 0.5 * (1.0 + np.cos(np.pi * (z[ramp] - taper_start_m)
                                  / (taper_end_m - taper_start_m)))
    w[z >= taper_end_m] = 0.0
    beta = beta * w
    return beta, beta * float(lidar_ratio_sr)


# ---------------------------------------------------------------------------
# Residual-overlap shapes
# ---------------------------------------------------------------------------
def parametric_overlap_error(z, depth: float, z_full_m: float) -> NDArray[np.float64]:
    """Smooth residual overlap error ``O_err(z) = 1 - depth * w(z)``.

    ``w`` is a raised cosine: 1 at the ground, 0 at ``z_full_m``, exactly 0 above -- the shape a
    real geometric overlap RATIO has (monotone, smooth, complete at a finite height). A POSITIVE
    ``depth`` means the applied overlap over-corrects so the L1 signal is too LOW near the ground
    (the documented network case: the generic 750 m table used where the true overlap only
    completes at 1300-1500 m). A negative ``depth`` is the opposite sign.
    """
    z = np.asarray(z, float)
    w = np.where(z < z_full_m, 0.5 * (1.0 + np.cos(np.pi * np.clip(z / z_full_m, 0, 1))), 0.0)
    return 1.0 - float(depth) * w


def measured_overlap_errors(z) -> dict:
    """The two MEASURED Payerne CHM15k residual-overlap errors, interpolated onto ``z``.

    Source ``doc/reports/figs_cl61_chm_tilt/tub140016_overlap_noise.npz`` (hood-noise retrieval,
    N = 8540 hood profiles, native 1024-gate CHM15k grid), holding three curves:

      ``O_applied``            what the TUB140016 firmware divides into the raw signal
                               (0.63 at 500 m, complete only ~1.5-2 km);
      ``O_ref_file``           the GENERIC TUB120011 table (complete ~750 m), the network default
                               where no per-module .cfg exists;
      ``O_true_reconstructed`` the module's current true overlap, reconstructed against the
                               co-located CL61 (0.55 at 500 m).

    Two residual errors follow:

      ``tilt``    = O_true_reconstructed / O_applied -- the REAL documented Payerne defect (the
                    2014 module has aged since its onboard calibration): -14 % at 400-500 m.
      ``generic`` = O_applied / O_ref_file           -- the network-wide scenario in which the
                    750 m generic scaffold is used on a unit completing at 1500 m: -25 % at
                    500 m, -19 % at 750 m, back to 1 by ~1500 m.
    """
    data = np.load(TUB_NPZ, allow_pickle=True)
    r = np.asarray(data["rng"], float)
    o_app = np.asarray(data["O_applied"], float)
    o_ref = np.asarray(data["O_ref_file"], float)
    o_true = np.asarray(data["O_true_reconstructed"], float)

    def _ratio(num, den):
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(den > 0.05, num / den, np.nan)
        # Both curves are 1 above ~2 km (they agree to < 3e-3 there); force exactly 1 so the fit
        # band carries no residual-interpolation noise and the attribution stays clean.
        ratio = np.where(r > 2000.0, 1.0, ratio)
        good = np.isfinite(ratio)
        return np.interp(z, r[good], ratio[good], left=float(ratio[good][0]), right=1.0)

    return dict(tilt=_ratio(o_true, o_app), generic=_ratio(o_app, o_ref))


def network_overlap_stats(z) -> dict:
    """Median / percentile completion height of the 132 measured network CHM15k overlaps.

    Source ``validation/overlap_clearsky/network_overlaps.npz`` (M2 clear-sky-noise retrieval,
    Apr-Jun 2026, 10 m grid 15-2495 m). Completion = lowest range above which O stays > 0.99.
    """
    data = np.load(NETWORK_NPZ, allow_pickle=True)
    r = np.asarray(data["rng"], float)
    curves = np.asarray(data["curves"], float)
    z99 = []
    for c in curves:
        idx = np.where(np.isfinite(c) & (c > 0.99))[0]
        z99.append(r[idx[0]] if idx.size else np.nan)
    z99 = np.asarray(z99, float)
    return dict(rng=r, curves=curves, z99=z99,
                z99_median=float(np.nanmedian(z99)),
                z99_p5=float(np.nanpercentile(z99, 5)),
                z99_p95=float(np.nanpercentile(z99, 95)),
                z99_max=float(np.nanmax(z99)), n=int(curves.shape[0]))


# ---------------------------------------------------------------------------
# One scenario -> one C_L ladder
# ---------------------------------------------------------------------------
def ladder(z, atmosphere, C_true, o_err=None, beta_aer=None, ext_aer=None,
           options=None, centres=CENTRES_M, half_width=HALF_WIDTH_M):
    """Forced-window C_L ladder. Returns ``(centres_m, C_L)`` in the units of ``C_true``."""
    rcs = forward_signal(z, C_true, atmosphere["beta_mol"], beta_aer=beta_aer, ext_aer=ext_aer,
                         overlap=o_err)
    values = []
    for centre in centres:
        res = retrieve(z, rcs, atmosphere, options=options,
                       window=(centre - half_width, centre + half_width))
        values.append(res["C_L"] if res["ok"] else np.nan)
    return np.asarray(centres, float), np.asarray(values, float)


def stats(centres, values, C_true):
    """(mean bias %, LS slope %/km over the ladder, peak-to-peak spread in percentage points)."""
    rel = (np.asarray(values, float) / float(C_true) - 1.0) * 100.0
    good = np.isfinite(rel)
    if good.sum() < 3:
        return np.nan, np.nan, np.nan
    slope = float(np.polyfit(np.asarray(centres, float)[good] / 1000.0, rel[good], 1)[0])
    return float(np.mean(rel[good])), slope, float(np.ptp(rel[good]))


def attribute(centres, cl_base, cl_err, C_true):
    """Overlap-attributed effect: the per-window ratio C_L(with error) / C_L(baseline).

    Returned as
      ``offset_pct``      mean of (ratio - 1) * 100 over the ladder  [%]
      ``slope_pct_km``    least-squares slope of (ratio - 1)*100 vs window centre  [%/km]
      ``window_spread_pp`` peak-to-peak of (ratio-1)*100 over the ladder  [percentage points]
    ``window_spread_pp`` is THE discriminator: 0 means a pure offset.
    """
    ratio = (np.asarray(cl_err, float) / np.asarray(cl_base, float) - 1.0) * 100.0
    good = np.isfinite(ratio)
    slope = float(np.polyfit(np.asarray(centres, float)[good] / 1000.0, ratio[good], 1)[0])
    return dict(offset_pct=float(np.mean(ratio[good])),
                slope_pct_km=slope,
                window_spread_pp=float(np.ptp(ratio[good])),
                per_window_pct=ratio.tolist())


# ---------------------------------------------------------------------------
# Mechanistic cross-check: predicted offset = exp(2 * Delta_OD) below z_full
# ---------------------------------------------------------------------------
def predicted_offset_from_delta_od(z, atmosphere, C_true, o_err, beta_aer, ext_aer, options,
                                   centre_m=4500.0, half_width_m=HALF_WIDTH_M):
    """Integrate the Klett extinction error below the window and predict the C_L offset.

    Equation [3] says C_L is multiplied by exp(2*Delta_OD) where
    ``Delta_OD = INT_0^z (ext_tot_klett[with error] - ext_tot_klett[baseline]) dz'``.
    We take Delta_OD at the window bottom; if the mechanism is understood, exp(2*Delta_OD)
    must equal the measured offset ratio. Returns (delta_od, predicted_offset_pct).
    """
    win = (centre_m - half_width_m, centre_m + half_width_m)
    rcs_b = forward_signal(z, C_true, atmosphere["beta_mol"], beta_aer=beta_aer,
                           ext_aer=ext_aer, overlap=None)
    rcs_e = forward_signal(z, C_true, atmosphere["beta_mol"], beta_aer=beta_aer,
                           ext_aer=ext_aer, overlap=o_err)
    rb = retrieve(z, rcs_b, atmosphere, options=options, window=win)
    re = retrieve(z, rcs_e, atmosphere, options=options, window=win)
    i_bot = int(np.searchsorted(z, win[0]))
    d_ext = np.asarray(re["ext_tot"], float) - np.asarray(rb["ext_tot"], float)
    d_od = float(np.trapezoid(d_ext[:i_bot], z[:i_bot]))
    return d_od, float((np.exp(2.0 * d_od) - 1.0) * 100.0)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = {}

    z = default_grid()
    C_true = PAYERNE_CHM15K["C_true"]
    atm = make_atmosphere(z, station_alt_m=PAYERNE_CHM15K["station_alt_m"],
                          wavelength_nm=PAYERNE_CHM15K["wavelength_nm"])
    opts = default_options()
    S_aer = float(opts.lidar_ratio_aerosol)

    beta_bl, ext_bl = capped_exponential_aerosol(z, BL_BETA_SURFACE, BL_SCALE_HEIGHT, S_aer,
                                                 BL_TAPER_START, BL_TAPER_END)
    aod_bl = float(np.trapezoid(ext_bl, z))
    i500 = int(np.argmin(np.abs(z - 500.0)))
    sr500 = float(1.0 + beta_bl[i500] / atm["beta_mol"][i500])
    print(f"grid: {z.size} gates, dz = {np.diff(z)[0]:.2f} m, z[0] = {z[0]:.2f} m; "
          f"C_true = {C_true:.3e}; S_aer = {S_aer:.0f} sr; ladder 2.5-6.5 km, +-490 m windows")
    print(f"BL aerosol: beta_s = {BL_BETA_SURFACE:.2e} m^-1 sr^-1, H = {BL_SCALE_HEIGHT:.0f} m, "
          f"tapered to 0 at {BL_TAPER_END:.0f} m; scattering ratio at 500 m = {sr500:.2f}; "
          f"AOD(1064 nm) = {aod_bl:.4f}")
    out["setup"] = dict(n_gates=int(z.size), dz_m=float(np.diff(z)[0]), z0_m=float(z[0]),
                        C_true=C_true, S_aer_sr=S_aer, wavelength_nm=1064.0,
                        ladder_centres_m=CENTRES_M.tolist(), half_width_m=HALF_WIDTH_M,
                        bl_beta_surface=BL_BETA_SURFACE, bl_scale_height_m=BL_SCALE_HEIGHT,
                        bl_taper_end_m=BL_TAPER_END,
                        bl_scattering_ratio_500m=sr500, bl_aod_1064=aod_bl)

    atmos_cases = dict(molecular=(None, None), bl_aerosol=(beta_bl, ext_bl))

    # ---- 0. CLOSURE ------------------------------------------------------------------------
    baselines = {}
    print("\n--- CLOSURE (O_err = 1) ------------------------------------------------------")
    for label, (ba, ea) in atmos_cases.items():
        c, v = ladder(z, atm, C_true, None, ba, ea, opts)
        baselines[label] = v
        bias, slope, spread = stats(c, v, C_true)
        print(f"  {label:11s}: mean bias {bias:+.4f} %, slope {slope:+.5f} %/km, "
              f"ladder spread {spread:.4f} pp")
        out.setdefault("closure", {})[label] = dict(mean_bias_pct=bias,
                                                    slope_pct_per_km=slope,
                                                    ladder_spread_pp=spread, cl=v.tolist())

    # A DOCUMENTED CONFOUND, kept for the record: the same aerosol NOT capped below the band.
    beta_in, ext_in = capped_exponential_aerosol(z, BL_BETA_SURFACE, BL_SCALE_HEIGHT, S_aer,
                                                 2200.0, 3000.0)
    c, v = ladder(z, atm, C_true, None, beta_in, ext_in, opts)
    bias, slope, spread = stats(c, v, C_true)
    print(f"  [confound] same aerosol reaching 3.0 km INTO the fit band: bias {bias:+.3f} %, "
          f"slope {slope:+.3f} %/km  <- aerosol, NOT overlap")
    out["aerosol_reaching_the_window"] = dict(taper_start_m=2200.0, taper_end_m=3000.0,
                                              mean_bias_pct=bias, slope_pct_per_km=slope,
                                              ladder_spread_pp=spread, cl=v.tolist(),
                                              note="belongs to the aerosol mechanism scan; "
                                                   "listed here only so it is not mistaken for "
                                                   "an overlap effect")

    # ---- 1. MEASURED Payerne residual overlaps ----------------------------------------------
    measured = measured_overlap_errors(z)
    net = network_overlap_stats(z)
    print(f"\n132 network CHM15k overlaps, completion height (O > 0.99): "
          f"median {net['z99_median']:.0f} m, p5 {net['z99_p5']:.0f} m, "
          f"p95 {net['z99_p95']:.0f} m, max {net['z99_max']:.0f} m")
    out["network_overlap_completion_m"] = {k: net[k] for k in
                                           ("z99_median", "z99_p5", "z99_p95", "z99_max", "n")}

    print("\n--- MEASURED overlap defects (effect attributed vs same-atmosphere baseline) ---")
    print("   scenario            atmosphere   offset [%]   slope [%/km]   window spread [pp]")
    measured_res = {}
    for name, o_err in measured.items():
        row = dict(O_err_at={f"{h:.0f}m": float(np.interp(h, z, o_err))
                             for h in (300, 400, 500, 750, 1000, 1500, 2000, 2500)})
        for atm_label, (ba, ea) in atmos_cases.items():
            c, v = ladder(z, atm, C_true, o_err, ba, ea, opts)
            att = attribute(c, baselines[atm_label], v, C_true)
            att["absolute_bias_pct"] = stats(c, v, C_true)[0]
            att["cl"] = v.tolist()
            row[atm_label] = att
            print(f"   {name:18s}  {atm_label:11s}  {att['offset_pct']:+9.3f}  "
                  f"{att['slope_pct_km']:+12.5f}   {att['window_spread_pp']:14.5f}")
        measured_res[name] = row
    out["measured"] = measured_res

    # mechanistic cross-check on the biggest measured case
    d_od, pred = predicted_offset_from_delta_od(z, atm, C_true, measured["generic"],
                                                beta_bl, ext_bl, opts)
    got = measured_res["generic"]["bl_aerosol"]["offset_pct"]
    print(f"\n   mechanism check (generic defect, BL aerosol, window 4.01-4.99 km):")
    print(f"     Delta_OD below the window = {d_od:.5e}  ->  exp(2*Delta_OD) - 1 = {pred:+.3f} %")
    print(f"     measured offset           = {got:+.3f} %   (agreement "
          f"{abs(pred-got):.3f} pp)")
    out["mechanism_check"] = dict(delta_od=d_od, predicted_offset_pct=pred,
                                  measured_offset_pct=got, residual_pp=abs(pred - got))

    # ---- 2. Parametric sweep: depth x completion height ---------------------------------------
    depths = [0.05, 0.10, 0.25]
    z_fulls = [500.0, 750.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0, 6000.0]
    print("\n--- PARAMETRIC SWEEP, BL-aerosol truth: offset [%] / slope [%/km] ---------------")
    print("  z_full [m] |" + "".join(f"    depth {d*100:.0f} %       |" for d in depths))
    sweep = {}
    for zf in z_fulls:
        cells = []
        for d in depths:
            o_err = parametric_overlap_error(z, d, zf)
            c, v = ladder(z, atm, C_true, o_err, beta_bl, ext_bl, opts)
            att = attribute(c, baselines["bl_aerosol"], v, C_true)
            att.update(z_full_m=zf, depth=d, cl=v.tolist())
            sweep[f"zfull{zf:.0f}_depth{d*100:.0f}"] = att
            cells.append(f" {att['offset_pct']:+8.3f} / {att['slope_pct_km']:+7.4f} |")
        print(f"  {zf:9.0f} |" + "".join(cells))
    out["parametric_sweep_bl"] = sweep

    sweep_mol = {}
    print("\n--- same sweep, PURE MOLECULAR truth (no aerosol to amplify the defect) --------")
    print("  z_full [m] |" + "".join(f"    depth {d*100:.0f} %       |" for d in depths))
    for zf in z_fulls:
        cells = []
        for d in depths:
            o_err = parametric_overlap_error(z, d, zf)
            c, v = ladder(z, atm, C_true, o_err, None, None, opts)
            att = attribute(c, baselines["molecular"], v, C_true)
            att.update(z_full_m=zf, depth=d)
            sweep_mol[f"zfull{zf:.0f}_depth{d*100:.0f}"] = att
            cells.append(f" {att['offset_pct']:+8.3f} / {att['slope_pct_km']:+7.4f} |")
        print(f"  {zf:9.0f} |" + "".join(cells))
    out["parametric_sweep_molecular"] = sweep_mol

    # ---- 3. Opposite sign --------------------------------------------------------------------
    print("\n--- OPPOSITE SIGN (applied overlap UNDER-corrects, depth -25 %) -----------------")
    neg = {}
    for zf in (1000.0, 1500.0, 2500.0, 4000.0):
        o_err = parametric_overlap_error(z, -0.25, zf)
        c, v = ladder(z, atm, C_true, o_err, beta_bl, ext_bl, opts)
        att = attribute(c, baselines["bl_aerosol"], v, C_true)
        att.update(z_full_m=zf, depth=-0.25)
        neg[f"zfull{zf:.0f}"] = att
        print(f"  z_full {zf:6.0f} m: offset {att['offset_pct']:+8.3f} %, "
              f"slope {att['slope_pct_km']:+.4f} %/km, spread {att['window_spread_pp']:.4f} pp")
    out["opposite_sign"] = neg

    # ---- 4. How big would the residual have to be to make -15 %/km? -------------------------
    # Push a residual overlap that is STILL evolving through the whole fit band. This is
    # physically impossible for a ceilometer (geometric overlap is complete well below 2 km on
    # every unit measured here) but it bounds the mechanism.
    print("\n--- BOUNDING: a residual still evolving through the fit band (unphysical) -------")
    bound = {}
    for zf, d in ((8000.0, 0.25), (10000.0, 0.50), (10000.0, 0.90)):
        o_err = parametric_overlap_error(z, d, zf)
        c, v = ladder(z, atm, C_true, o_err, beta_bl, ext_bl, opts)
        att = attribute(c, baselines["bl_aerosol"], v, C_true)
        att.update(z_full_m=zf, depth=d,
                   O_err_2km=float(np.interp(2000.0, z, o_err)),
                   O_err_6km=float(np.interp(6000.0, z, o_err)))
        bound[f"zfull{zf:.0f}_depth{d*100:.0f}"] = att
        print(f"  z_full {zf:6.0f} m, depth {d*100:3.0f} %  (O_err {att['O_err_2km']:.3f} at "
              f"2 km -> {att['O_err_6km']:.3f} at 6 km): slope {att['slope_pct_km']:+.3f} %/km")
    out["bounding_unphysical"] = bound

    # ---- 5. Co-located CHM15k vs CL61 (part c) ----------------------------------------------
    # CL61: same L2 grid (300 s x 30 m), 910.55 nm, C_true = 1 (the L1 rcs_0 of a CL61 already is
    # an attenuated backscatter). Its firmware overlap completes by ~350 m and the hood-noise
    # retrieval reproduces it to +-1 %, so its residual error is at most a few % over a very short
    # path. The CHM15k carries the measured defects.
    atm_cl61 = make_atmosphere(z, station_alt_m=PAYERNE_CHM15K["station_alt_m"],
                               wavelength_nm=910.55)
    beta_bl61, ext_bl61 = capped_exponential_aerosol(z, BL_BETA_SURFACE * 1.15, BL_SCALE_HEIGHT,
                                                     S_aer, BL_TAPER_START, BL_TAPER_END)
    c, base_cl61 = ladder(z, atm_cl61, 1.0, None, beta_bl61, ext_bl61, opts)
    bias61, slope61, _ = stats(c, base_cl61, 1.0)
    print(f"\n--- PART (c) CO-LOCATED PAIR ---------------------------------------------------")
    print(f"  CL61 baseline closure (O_err = 1, 910.55 nm): bias {bias61:+.4f} %, "
          f"slope {slope61:+.5f} %/km")
    pair = dict(_cl61_baseline=dict(mean_bias_pct=bias61, slope_pct_per_km=slope61,
                                    cl=base_cl61.tolist()),
                _chm_baseline=dict(cl=baselines["bl_aerosol"].tolist()))
    scen = {
        "CHM15k_generic_defect": (atm, C_true, measured["generic"], beta_bl, ext_bl,
                                  baselines["bl_aerosol"]),
        "CHM15k_tilt_defect": (atm, C_true, measured["tilt"], beta_bl, ext_bl,
                               baselines["bl_aerosol"]),
        "CL61_1pct_to_350m": (atm_cl61, 1.0, parametric_overlap_error(z, 0.01, 350.0),
                              beta_bl61, ext_bl61, base_cl61),
        "CL61_5pct_to_350m": (atm_cl61, 1.0, parametric_overlap_error(z, 0.05, 350.0),
                              beta_bl61, ext_bl61, base_cl61),
    }
    for name, (aa, ct, oe, ba, ea, base) in scen.items():
        c, v = ladder(z, aa, ct, oe, ba, ea, opts)
        att = attribute(c, base, v, ct)
        att["cl"] = v.tolist()
        att["C_true"] = ct
        pair[name] = att
        print(f"  {name:24s}: offset {att['offset_pct']:+8.3f} %, "
              f"slope {att['slope_pct_km']:+.5f} %/km, spread {att['window_spread_pp']:.5f} pp")
    d_off = pair["CHM15k_generic_defect"]["offset_pct"] - pair["CL61_5pct_to_350m"]["offset_pct"]
    d_slp = (pair["CHM15k_generic_defect"]["slope_pct_km"]
             - pair["CL61_5pct_to_350m"]["slope_pct_km"])
    pair["difference_chm_minus_cl61"] = dict(offset_pp=d_off, slope_pct_km=d_slp)
    print(f"  CHM15k(generic) - CL61(5 %): offset difference {d_off:+.3f} pp, "
          f"slope difference {d_slp:+.5f} %/km")
    out["colocated_pair"] = pair

    # ---- 6. Figures ---------------------------------------------------------------------------
    _fig_overlap_shapes(plt, z, measured, net, out)
    _fig_ladders(plt, z, out, measured_res, sweep, bound, baselines, C_true)
    _fig_pair(plt, out, pair, C_true)

    json_path = REPO / "rayleigh_availability" / "overlap_scan_results.json"
    json_path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"\nJSON -> {json_path}")
    return out


# ---------------------------------------------------------------------------
# Figures (altitude always on Y, landscape, dpi 150)
# ---------------------------------------------------------------------------
def _fig_overlap_shapes(plt, z, measured, net, out):
    fig, axes = plt.subplots(1, 3, figsize=(16.6, 6.3))

    ax = axes[0]
    data = np.load(TUB_NPZ, allow_pickle=True)
    r = np.asarray(data["rng"], float)
    for key, color, label in (("O_applied", "C0",
                               "$O_{applied}$ - TUB140016 firmware (hood noise)"),
                              ("O_ref_file", "C1", "$O_{ref}$ - generic TUB120011 table"),
                              ("O_true_reconstructed", "C3",
                               "$O_{true}$ - reconstructed vs co-located CL61")):
        c = np.asarray(data[key], float)
        good = np.isfinite(c) & (r < 3000)
        ax.plot(c[good], r[good] / 1000.0, color=color, lw=1.9, label=label)
    ax.set_xlim(0, 1.15)
    ax.set_ylim(0, 2.5)
    ax.set_xlabel("overlap  $O(z)$   [-]")
    ax.set_ylabel("range above the instrument  [km AGL]")
    ax.set_title("MEASURED Payerne CHM15k overlaps\n(doc/reports/08, hood-noise retrieval)",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=7.6)

    ax = axes[1]
    for name, color, label in (("generic", "C1",
                                "generic 750 m table on a 1500 m overlap\n($-$25 % at 500 m)"),
                               ("tilt", "C3",
                                "Payerne module tilt\n($-$14 % at 450 m)")):
        ax.plot(measured[name], z / 1000.0, color=color, lw=2.0, label=label)
    for depth, ls in ((0.05, ":"), (0.10, "--"), (0.25, "-")):
        ax.plot(parametric_overlap_error(z, depth, 1500.0), z / 1000.0, color="0.45", lw=1.1,
                ls=ls, label=f"parametric depth {depth*100:.0f} %, $z_{{full}}$ = 1.5 km")
    ax.axhspan(2.01, 6.99, color="C2", alpha=0.10)
    ax.text(0.735, 4.3, "molecular fit band\nevery window lives here", fontsize=8.5, color="C2")
    ax.set_xlim(0.7, 1.05)
    ax.set_ylim(0, 7)
    ax.axvline(1.0, color="k", lw=0.9)
    ax.set_xlabel("residual overlap error  $O_{err}=O_{true}/O_{applied}$   [-]")
    ax.set_ylabel("range above the instrument  [km AGL]")
    ax.set_title("What the pipeline actually sees\n$O_{err}$ is EXACTLY 1 above 2 km",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=7.2)

    ax = axes[2]
    r2, curves = net["rng"], net["curves"]
    for c in curves[::3]:
        ax.plot(c, r2 / 1000.0, color="0.78", lw=0.4, alpha=0.7)
    ax.plot(np.nanmedian(curves, axis=0), r2 / 1000.0, color="C0", lw=2.3,
            label=f"network median (n = {net['n']})")
    ax.plot(np.nanpercentile(curves, 5, axis=0), r2 / 1000.0, color="C3", lw=1.4, ls="--",
            label="5th percentile")
    ax.plot(np.nanpercentile(curves, 95, axis=0), r2 / 1000.0, color="C2", lw=1.4, ls="--",
            label="95th percentile")
    ax.axhline(2.0, color="C2", lw=1.6)
    ax.text(0.05, 2.06, "bottom of the fit band", fontsize=8, color="C2")
    ax.set_xlim(0, 1.1)
    ax.set_ylim(0, 2.5)
    ax.set_xlabel("applied overlap  $O(z)$   [-]")
    ax.set_ylabel("range above the instrument  [km AGL]")
    ax.set_title(f"{net['n']} network CHM15k applied overlaps\n"
                 f"completion ($O>0.99$): median {net['z99_median']:.0f} m,\n"
                 f"p95 {net['z99_p95']:.0f} m, worst unit {net['z99_max']:.0f} m", fontsize=10.5)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)

    fig.suptitle("Overlap functions in play - measured on the real instruments, not assumed",
                 fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = FIGDIR / "overlap_01_functions.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    out.setdefault("figures", []).append(str(path))
    print(f"figure -> {path}")


def _fig_ladders(plt, z, out, measured_res, sweep, bound, baselines, C_true):
    fig, axes = plt.subplots(1, 3, figsize=(16.6, 6.3))
    y = CENTRES_M / 1000.0

    ax = axes[0]
    ax.axvline(0.0, color="k", lw=1.0)
    for name, color, label in (("generic", "C1", "generic-table defect ($-$25 % at 500 m)"),
                               ("tilt", "C3", "Payerne module tilt ($-$14 % at 450 m)")):
        row = measured_res[name]["bl_aerosol"]
        ax.plot(row["per_window_pct"], y, "o-", color=color, ms=6, lw=1.8,
                label=f"{label}\n  offset {row['offset_pct']:+.2f} %, "
                      f"slope {row['slope_pct_km']:+.5f} %/km")
        rowm = measured_res[name]["molecular"]
        ax.plot(rowm["per_window_pct"], y, "s--", color=color, ms=4, lw=1.1, alpha=0.55,
                label=f"  pure-molecular truth: {rowm['offset_pct']:+.2f} %")
    ax.set_xlabel("overlap-attributed change of $C_L$   [%]\n"
                  "(same night with vs without the residual overlap error)")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title("MEASURED overlap defects\na perfectly vertical line = a pure OFFSET",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left", fontsize=7.2)

    ax = axes[1]
    ax.axvline(0.0, color="k", lw=1.0)
    cmap = plt.get_cmap("viridis")
    zfs = sorted({v["z_full_m"] for v in sweep.values()})
    for i, zf in enumerate(zfs):
        row = sweep[f"zfull{zf:.0f}_depth25"]
        ax.plot(row["per_window_pct"], y, "o-", ms=4, lw=1.5,
                color=cmap(i / max(1, len(zfs) - 1)),
                label=f"$z_{{full}}$ = {zf/1000:.1f} km  ({row['slope_pct_km']:+.3f} %/km)")
    ax.set_xlabel("overlap-attributed change of $C_L$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title("Residual depth 25 %, completion height swept\n"
                 "(BL-aerosol truth, aerosol capped at 1.9 km)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left", fontsize=7)

    ax = axes[2]
    for depth, marker in ((0.05, "^"), (0.10, "s"), (0.25, "o")):
        xs = [sweep[f"zfull{zf:.0f}_depth{depth*100:.0f}"]["slope_pct_km"] for zf in zfs]
        ax.plot(xs, [zf / 1000.0 for zf in zfs], marker + "-", ms=5, lw=1.5,
                label=f"residual depth {depth*100:.0f} %")
    for key, row in bound.items():
        ax.plot([row["slope_pct_km"]], [row["z_full_m"] / 1000.0], "*", ms=13, color="0.3")
    ax.annotate("unphysical bound:\n$O_{err}$ still evolving\nthrough the whole band",
                xy=(bound[list(bound)[-1]]["slope_pct_km"],
                    bound[list(bound)[-1]]["z_full_m"] / 1000.0),
                xytext=(0.30, 0.62), textcoords="axes fraction", fontsize=7.5, color="0.3",
                arrowprops=dict(arrowstyle="->", color="0.3", lw=0.9))
    ax.axvline(0.0, color="k", lw=1.0)
    ax.axvline(-15.0, color="C3", lw=2.2, ls="--", label="observed Payerne CHM15k $-$15 %/km")
    ax.axvline(+15.0, color="C0", lw=2.2, ls=":", label="observed Aosta $\\approx$ +15 %/km")
    ax.axhspan(0.0, 2.0, color="C2", alpha=0.08)
    ax.text(0.02, 0.03, "every measured ALC overlap\ncompletes inside this band",
            transform=ax.transAxes, fontsize=7.5, color="C2")
    ax.set_xlabel("overlap-attributed ladder slope  $dC_L/dz$   [% per km]")
    ax.set_ylabel("overlap completion height $z_{full}$  [km AGL]")
    ax.set_title("THE DISCRIMINATOR\nslope produced vs slope observed", fontsize=11)
    ax.set_xscale("symlog", linthresh=0.001)
    ax.set_ylim(0, 10.6)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=7.4)

    fig.suptitle("Overlap error -> $C_L$ ladder: a pure OFFSET, no slope   "
                 "(forced windows $\\pm$490 m, centres 2.5-6.5 km AGL)", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = FIGDIR / "overlap_02_ladders.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    out.setdefault("figures", []).append(str(path))
    print(f"figure -> {path}")


def _fig_pair(plt, out, pair, C_true):
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 6.3),
                             gridspec_kw=dict(width_ratios=[1.0, 1.05]))
    y = CENTRES_M / 1000.0

    ax = axes[0]
    ax.axvline(0.0, color="k", lw=1.0)
    style = {"CHM15k_generic_defect": ("C1", "o-", "CHM15k 1064 nm, generic-table defect"),
             "CHM15k_tilt_defect": ("C3", "o-", "CHM15k 1064 nm, measured module tilt"),
             "CL61_5pct_to_350m": ("C2", "s-", "CL61 910 nm, 5 % residual to 350 m"),
             "CL61_1pct_to_350m": ("C4", "s--", "CL61 910 nm, 1 % residual to 350 m")}
    for name, (color, fmt, label) in style.items():
        row = pair[name]
        ax.plot(row["per_window_pct"], y, fmt, color=color, ms=6, lw=1.7,
                label=f"{label}\n  {row['offset_pct']:+.2f} %, "
                      f"{row['slope_pct_km']:+.5f} %/km")
    ax.set_xlabel("overlap-attributed change of $C_L$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title("Co-located pair: CHM15k (overlap to 1.5 km)\n"
                 "vs CL61 (overlap to 0.35 km)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left", fontsize=7.4)

    ax = axes[1]
    ax.axis("off")
    d = pair["difference_chm_minus_cl61"]
    text = (
        "PART (c)   CHM15k minus CL61, overlap difference ALONE\n"
        "=========================================================\n"
        f"  offset difference : {d['offset_pp']:+.3f} percentage points\n"
        f"  slope  difference : {d['slope_pct_km']:+.5f} %/km\n"
        "\n"
        "OBSERVED co-located differences (recovered vs kept nights)\n"
        "=========================================================\n"
        "  Payerne  CHM15k A  -22.5 %  |  CL61 C  -17.5 %  ->  5.0 pp apart\n"
        "  Aosta    CHM15k 0  +31.7 %  |  CL61 B  +32.1 %  ->  0.4 pp apart\n"
        "\n"
        "WHY OVERLAP CANNOT CONTRIBUTE AT ALL\n"
        "=========================================================\n"
        "  The observed quantity is a RATIO within ONE instrument\n"
        "  (recovered nights / kept nights).  A residual overlap error\n"
        "  is a multiplicative constant on that instrument's C_L, so it\n"
        "  cancels EXACTLY in the ratio - whatever its size.\n"
        "\n"
        "  And the two co-located units, whose overlaps differ by ~1.2 km\n"
        "  in completion height and by 25 percentage points at 500 m,\n"
        "  move TOGETHER: same sign, 0.4-5.0 pp apart.  If overlap drove\n"
        "  the altitude dependence they would diverge by tens of percent."
    )
    ax.text(0.0, 0.99, text, transform=ax.transAxes, fontsize=9.6, va="top",
            family="monospace",
            bbox=dict(boxstyle="round,pad=0.6", fc="#f6f6f6", ec="0.6"))

    fig.suptitle("Part (c) - would the CHM15k / CL61 overlap difference explain the "
                 "co-located spread?", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path = FIGDIR / "overlap_03_colocated_pair.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    out.setdefault("figures", []).append(str(path))
    print(f"figure -> {path}")


if __name__ == "__main__":
    main()
