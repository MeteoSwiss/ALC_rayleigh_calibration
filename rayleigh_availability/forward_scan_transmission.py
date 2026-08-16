# -*- coding: utf-8 -*-
r"""MECHANISM SCAN 1 -- below-window aerosol transmission and the lidar-ratio assumption.

QUESTION
--------
The retrieved lidar constant must not depend on the height of the molecular fit window, yet it
does: v2.2-recovered nights (fitted 1-2 km higher than the v2.0-kept ones) return -22.5 % at
Payerne CHM15k, +32 % at Aosta. The leading candidate for the NEGATIVE gradients is the aerosol
below the window: the shipped chain has to undo the two-way transmission of that aerosol, and it
does so with a Klett inversion using an ASSUMED lidar ratio (``options.lidar_ratio_aerosol`` =
52 sr). Whatever the Klett fails to undo lands directly in

    C_L(z) = rcs(z) / beta_tot_ret(z) * exp(2 * OD_ret(z)),     OD_ret integrated from the ground.

Substituting the true signal rcs(z) = C_true * beta_tot(z) * exp(-2 * OD_true(z)) gives the exact
error law this module measures:

    C_L(z) / C_true = [beta_tot(z) / beta_tot_ret(z)] * exp(2 * (OD_ret(z) - OD_true(z)))    [1]
                      \_____ term B: backscatter __/   \______ term T: transmission ______/

    d ln C_L / dz   = 2 * (alpha_ret(z) - alpha_true(z))  +  d ln(beta_tot/beta_tot_ret)/dz  [2]

Equation [2] already says something the observational study could not: an extinction error BELOW
the ladder enters [1] as a CONSTANT -- the optical-depth error is fully accumulated under the
lowest window and stops growing -- so it can only OFFSET C_L, never tilt it. Only a mis-corrected
extinction INSIDE the 3-6 km band can tilt it, through term T; and term B tilts it whenever the
aerosol backscatter the retrieval calls "molecular" changes with height.

The scans below measure term T and term B SEPARATELY on the shipped code, because the AOD sweep
and the lidar-ratio sweep are not independent: at fixed AOD, changing S_true changes beta_aer
(= AOD/(S*thickness)) and therefore term B as well. ``scan_isosr`` fixes the backscatter profile
(hence term B) and sweeps S_true alone, which is the only clean way to isolate term T.

WHAT IS SIMULATED
-----------------
Forward model: ``rayleigh_availability.forward_model`` (imported, never rebuilt). Its closure on a
pure molecular atmosphere is re-run at the top of every execution and must stay under
``CLOSURE_TOL_PCT`` = 0.5 %.

Scenes, all on the real Payerne CHM15k post-binning grid (512 gates, dz = 29.97 m, 1064 nm,
station 490 m ASL):

  * ``below``     -- boundary-layer top-hat, base 100 m, top 1.0-2.0 km, i.e. ENTIRELY under the
                     ladder. AOD 0.02-0.50, TRUE lidar ratio S_true 20-90 sr.
  * ``intruding`` -- the same top-hat with the top swept 1.0-3.5 km, so it sits under part of the
                     ladder and inside the rest.
  * ``lofted``    -- a 1 km thick layer whose base is swept 3.0 / 4.0 / 5.0 km, i.e. INSIDE the
                     ladder: the pure term-T experiment, because the windows below it see no
                     aerosol at all and the windows above see only its optical depth.
  * ``haze``      -- exponentially decaying haze (scale height 1-5 km), the only shape that puts
                     aerosol backscatter AND extinction throughout the 3-6 km band.

Unless stated otherwise a scene is normalised by its AOD: AOD = INT ext_aer dz over the whole grid
(0-15.3 km), with ext_aer = beta_aer * S_true.  ``iso_sr_haze`` normalises by the SCATTERING RATIO
at 4.5 km instead, which is what keeps term B fixed while S_true varies.

PRODUCTION GATES ON A FORCED WINDOW
-----------------------------------
The ladder forces windows the production search would never accept. Each rung is tagged with the
two operational gates of ``options.json`` (``min_window_r2`` = 0.5, ``max_window_rel_error`` =
50 %) and dC/dz is fitted on the rungs that PASS. Failing rungs are kept in the output and drawn
hollow, because they show how violent the contamination is, but never enter a slope.

Run:  python rayleigh_availability/forward_scan_transmission.py
Out:  doc/reports/figs_altitude_audit/fwd_trans_*.png
      rayleigh_availability/fwd_transmission_results.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
if str(REPO / "rayleigh_availability") not in sys.path:
    sys.path.insert(0, str(REPO / "rayleigh_availability"))

from forward_model import (                                              # noqa: E402
    PAYERNE_CHM15K,
    aerosol_layer,
    default_grid,
    default_options,
    exponential_aerosol,
    forward_signal,
    make_atmosphere,
    retrieve,
    run_closure_tests,
)

FIGDIR = REPO / "doc" / "reports" / "figs_altitude_audit"
OUTJSON = REPO / "rayleigh_availability" / "fwd_transmission_results.json"
RATIO_JSON = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability"
                  "/ratio_profile.json")

# --- ladder definition -------------------------------------------------------------------------
CENTRES_M = np.arange(2500.0, 6501.0, 250.0)   # window centres, m AGL (17 rungs)
HALF_M = 490.0                                 # production half-length option
SLOPE_BAND = (3000.0, 6000.0)                  # dC/dz is ALWAYS quoted between these centres
Z_SR = 4500.0                                  # where scattering ratio is quoted
CLOSURE_TOL_PCT = 0.5

# --- production gates (options.json) -----------------------------------------------------------
MIN_R2 = 0.5
MAX_REL_ERROR_PCT = 50.0
MIN_SLOPE_RUNGS = 6            # gate-passing rungs required before a dC/dz is quoted
MIN_SLOPE_SPAN_M = 1500.0      # ... and the height range they must still cover

# --- observational targets (measured; see the task brief and ratio_profile.json) ----------------
TARGET_DCDZ_PCT_PER_KM = -15.0        # Payerne CHM15k: -22.5 % over the ~1.5 km window offset
TARGET_R_SLOPE_PCT_PER_KM = -19.5     # signal/p_mol log-slope, 2-6 km, recovered nights
TARGET_R_SLOPE_KEPT_PCT_PER_KM = -3.7   # same, v2.0-kept nights


# ================================================================================================
# scene construction
# ================================================================================================
def scene(z, kind, *, aod, lidar_ratio_true, top_m=1500.0, base_m=100.0,
          scale_height_m=1500.0, thickness_m=1000.0):
    """Build (beta_aer [m^-1 sr^-1], ext_aer [m^-1]) for one named scene, normalised to ``aod``.

    ``kind``: 'none' | 'below' | 'intruding' (top-hat base_m -> top_m) | 'lofted' (top-hat
    base_m -> base_m + thickness_m) | 'haze' (beta_surface * exp(-z/scale_height_m)).
    ``lidar_ratio_true`` is S_true in sr; ext_aer = beta_aer * S_true. The retrieval always assumes
    52 sr, so S_true != 52 is the lidar-ratio error under test.
    """
    if kind == "none" or aod <= 0:
        zeros = np.zeros_like(z)
        return zeros, zeros.copy()
    if kind in ("below", "intruding"):
        beta, ext = aerosol_layer(z, base_m, top_m, 1.0, lidar_ratio_true)
    elif kind == "lofted":
        beta, ext = aerosol_layer(z, base_m, base_m + thickness_m, 1.0, lidar_ratio_true)
    elif kind == "haze":
        beta, ext = exponential_aerosol(z, 1.0, scale_height_m, lidar_ratio_true)
    else:
        raise ValueError(f"unknown scene {kind!r}")
    unit_aod = float(np.trapezoid(ext, z))
    if unit_aod <= 0:
        zeros = np.zeros_like(z)
        return zeros, zeros.copy()
    factor = float(aod) / unit_aod
    return beta * factor, ext * factor


def iso_sr_haze(z, atmosphere, sr_target, lidar_ratio_true, *, scale_height_m=2000.0,
                z_ref_m=Z_SR):
    """Exponential haze normalised so the SCATTERING RATIO at ``z_ref_m`` equals ``sr_target``.

    This is the control that separates the two terms of eq. [2]: beta_aer -- and hence term B,
    the "the retrieval calls this window molecular" error -- is IDENTICAL for every ``S_true``,
    while ext_aer = beta_aer * S_true and therefore term T varies. Any dC/dz that moves when
    S_true moves, at fixed ``sr_target``, is term T and nothing else.

    Returns (beta_aer, ext_aer, aod).
    """
    beta_unit, _ = exponential_aerosol(z, 1.0, scale_height_m, lidar_ratio_true)
    j = int(np.argmin(np.abs(z - float(z_ref_m))))
    beta_mol_ref = float(atmosphere["beta_mol"][j])
    want = (float(sr_target) - 1.0) * beta_mol_ref
    factor = want / beta_unit[j] if beta_unit[j] > 0 else 0.0
    beta_aer = beta_unit * factor
    ext_aer = beta_aer * float(lidar_ratio_true)
    return beta_aer, ext_aer, float(np.trapezoid(ext_aer, z))


# ================================================================================================
# one ladder
# ================================================================================================
def ladder(z, atmosphere, beta_aer, ext_aer, C_true, *, options=None,
           centres_m=CENTRES_M, half_m=HALF_M, perturbations=False, keep_rcs=False):
    """Force every window in ``centres_m`` and run the SHIPPED retrieval on the same profile."""
    rcs = forward_signal(z, C_true, atmosphere["beta_mol"], beta_aer, ext_aer)
    centres = np.asarray(centres_m, float)
    dev, r2, rel, med, unc = (np.full(centres.size, np.nan) for _ in range(5))
    for i, centre in enumerate(centres):
        res = retrieve(z, rcs, atmosphere, options=options,
                       window=(centre - half_m, centre + half_m),
                       perturbations=perturbations)
        if not res["ok"]:
            continue
        dev[i] = (res["C_L"] / C_true - 1.0) * 100.0
        r2[i] = res["r_squared"]
        rel[i] = res["rel_error"]
        med[i] = res["C_L_median"]
        unc[i] = res["C_L_uncertainty"]

    passes = np.isfinite(dev) & (r2 >= MIN_R2) & (rel <= MAX_REL_ERROR_PCT)
    band = (centres >= SLOPE_BAND[0]) & (centres <= SLOPE_BAND[1]) & passes
    # A slope is only quoted when at least MIN_SLOPE_RUNGS gate-passing rungs survive AND they
    # still span MIN_SLOPE_SPAN_M. Without this, a scene whose lower windows are all rejected
    # leaves 3 clustered rungs and the fit returns a meaningless double-digit %/km.
    slope = np.nan
    if band.sum() >= MIN_SLOPE_RUNGS and (centres[band].max() - centres[band].min()) >= MIN_SLOPE_SPAN_M:
        slope = float(np.polyfit(centres[band] / 1000.0, dev[band], 1)[0])

    def at(target):
        j = int(np.argmin(np.abs(centres - target)))
        return float(dev[j]) if passes[j] else np.nan

    res = dict(centres_m=centres.tolist(), dev_pct=dev.tolist(),
               r2=r2.tolist(), rel_error_pct=rel.tolist(),
               passes=passes.tolist(),
               dC_dz_pct_per_km=slope, n_pass=int(band.sum()),
               n_pass_all=int(passes.sum()),
               dev_at_3km=at(3000.0), dev_at_6km=at(6000.0),
               C_L_median=med.tolist(), C_L_unc=unc.tolist(),
               r_slope_pct_per_km=signal_ratio_slope(z, rcs, atmosphere))
    if keep_rcs:
        res["_rcs"] = rcs
    return res


def scattering_ratio(atmosphere, beta_aer, z_target_m):
    """(beta_mol + beta_aer)/beta_mol at ``z_target_m`` -- the plausibility yardstick."""
    z = atmosphere["z"]
    j = int(np.argmin(np.abs(z - float(z_target_m))))
    return float(1.0 + beta_aer[j] / atmosphere["beta_mol"][j])


def signal_ratio_slope(z, rcs, atm, lo=2000.0, hi=6000.0):
    """Log-slope of R(z) = signal/p_mol between ``lo`` and ``hi``, in %/km.

    The SAME observable ``rayleigh_availability/ratio_profile.py`` measures on real nights
    (-3.7 %/km on Payerne v2.0-kept nights, -19.5 %/km on the v2.2-recovered ones), so a simulated
    scene can be matched to the data BEFORE its C_L ladder is believed.
    """
    signal = np.asarray(rcs) / z ** 2
    ratio = signal / atm["p_mol"]
    m = (z >= lo) & (z <= hi) & np.isfinite(ratio) & (ratio > 0)
    if m.sum() < 5:
        return float("nan")
    return float(np.polyfit(z[m] / 1000.0, np.log(ratio[m]), 1)[0] * 100.0)


# ================================================================================================
# scans
# ================================================================================================
def run_all(verbose=True):
    z = default_grid()
    C_true = PAYERNE_CHM15K["C_true"]
    atm = make_atmosphere(z, station_alt_m=PAYERNE_CHM15K["station_alt_m"],
                          wavelength_nm=PAYERNE_CHM15K["wavelength_nm"])
    options = default_options()
    lr_assumed = float(options.lidar_ratio_aerosol)
    out = dict(C_true=C_true, lidar_ratio_assumed=lr_assumed,
               centres_m=CENTRES_M.tolist(), half_width_m=HALF_M,
               slope_band_m=list(SLOPE_BAND),
               gates=dict(min_r2=MIN_R2, max_rel_error_pct=MAX_REL_ERROR_PCT))

    # ---- 0. closure re-check (mandatory before any interpretation) ----------------------------
    clos = run_closure_tests(fig_path=None, n_seeds=30, verbose=False)
    out["closure"] = dict(max_abs_dev_pct=clos["clean_max_abs_dev_pct"],
                          slope_pct_per_km=clos["clean_slope_pct_per_km"],
                          linearity_max_abs_err=clos["linearity_max_abs_err"])
    if clos["clean_max_abs_dev_pct"] > CLOSURE_TOL_PCT:
        raise RuntimeError(f"CLOSURE FAILED: {clos['clean_max_abs_dev_pct']:.4f} % "
                           f"> {CLOSURE_TOL_PCT} % -- fix the simulator before interpreting")
    if verbose:
        print(f"closure: max |dev| = {clos['clean_max_abs_dev_pct']:.4f} %, residual slope = "
              f"{clos['clean_slope_pct_per_km']:+.5f} %/km, linearity err "
              f"{clos['linearity_max_abs_err']:.1e}  [PASS < {CLOSURE_TOL_PCT} %]")

    aods = [0.02, 0.05, 0.10, 0.20, 0.35, 0.50]
    lrs = [20.0, 30.0, 40.0, 52.0, 65.0, 80.0, 90.0]

    # ---- A. layer ENTIRELY below the ladder ----------------------------------------------------
    scan_a = []
    for aod in aods:
        for lr_true in lrs:
            ba, ea = scene(z, "below", aod=aod, lidar_ratio_true=lr_true, top_m=1500.0)
            lad = ladder(z, atm, ba, ea, C_true, options=options)
            scan_a.append(dict(aod=aod, lidar_ratio_true=lr_true, top_m=1500.0,
                               sr_in_layer_800m=scattering_ratio(atm, ba, 800.0),
                               **{k: lad[k] for k in ("dC_dz_pct_per_km", "dev_at_3km",
                                                      "dev_at_6km", "n_pass_all",
                                                      "dev_pct", "passes",
                                                      "r_slope_pct_per_km")}))
    out["scan_below"] = scan_a

    # ---- B. intruding layer: top swept through the bottom of the ladder ------------------------
    scan_b = []
    for top in (1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0):
        for aod in (0.05, 0.20):
            for lr_true in (30.0, 52.0, 90.0):
                ba, ea = scene(z, "intruding", aod=aod, lidar_ratio_true=lr_true, top_m=top)
                lad = ladder(z, atm, ba, ea, C_true, options=options)
                scan_b.append(dict(top_m=top, aod=aod, lidar_ratio_true=lr_true,
                                   sr_below_top=scattering_ratio(atm, ba, max(top - 300.0, 300.)),
                                   **{k: lad[k] for k in ("dC_dz_pct_per_km", "dev_at_3km",
                                                          "dev_at_6km", "n_pass_all",
                                                          "dev_pct", "passes", "r2",
                                                          "rel_error_pct")}))
    out["scan_intruding"] = scan_b

    # ---- C. lofted layer INSIDE the ladder: the pure term-T experiment -------------------------
    # A 1 km thick layer at 4 km with AOD 0.005 already has a scattering ratio of ~3; anything
    # thicker is a cloud, so the AOD list is deliberately an order of magnitude below scene A's.
    scan_c = []
    for base in (3000.0, 4000.0, 5000.0):
        for aod in (0.001, 0.002, 0.005, 0.010):
            for lr_true in (20.0, 30.0, 52.0, 65.0, 90.0):
                ba, ea = scene(z, "lofted", aod=aod, lidar_ratio_true=lr_true,
                               base_m=base, thickness_m=1000.0)
                lad = ladder(z, atm, ba, ea, C_true, options=options)
                dev = np.asarray(lad["dev_pct"], float)
                ok = np.asarray(lad["passes"], bool)
                scan_c.append(dict(base_m=base, top_m=base + 1000.0, aod=aod,
                                   lidar_ratio_true=lr_true,
                                   sr_in_layer=scattering_ratio(atm, ba, base + 500.0),
                                   step_pct=_step_across(lad, base, base + 1000.0),
                                   # worst error a window that PASSES the production gates can
                                   # still commit -- the gates are not a contamination filter
                                   max_dev_passing_pct=float(np.nanmax(np.abs(dev[ok])))
                                   if ok.any() else np.nan,
                                   **{k: lad[k] for k in ("dC_dz_pct_per_km", "dev_at_3km",
                                                          "dev_at_6km", "n_pass_all",
                                                          "dev_pct", "passes")}))
    out["scan_lofted"] = scan_c

    # ---- D. deep haze: aerosol present throughout the ladder -----------------------------------
    scan_d = []
    for h in (1000.0, 1500.0, 2000.0, 3000.0, 5000.0):
        for aod in (0.02, 0.05, 0.10, 0.20, 0.35, 0.50):
            for lr_true in lrs:
                ba, ea = scene(z, "haze", aod=aod, lidar_ratio_true=lr_true, scale_height_m=h)
                lad = ladder(z, atm, ba, ea, C_true, options=options)
                j = int(np.argmin(np.abs(z - Z_SR)))
                band36 = (z >= 3000.0) & (z <= 6000.0)
                scan_d.append(dict(scale_height_m=h, aod=aod, lidar_ratio_true=lr_true,
                                   sr_at_4500m=scattering_ratio(atm, ba, Z_SR),
                                   alpha_true_4500m=float(ea[j]),
                                   aod_3_6km=float(np.trapezoid(ea[band36], z[band36])),
                                   beta_aer_4500m=float(ba[j]),
                                   **{k: lad[k] for k in ("dC_dz_pct_per_km", "dev_at_3km",
                                                          "dev_at_6km", "n_pass_all",
                                                          "r_slope_pct_per_km",
                                                          "dev_pct", "passes")}))
    out["scan_haze"] = scan_d

    # ---- D2. ISO-SR decomposition: term B frozen, term T swept ---------------------------------
    scan_iso = []
    for sr in (1.05, 1.10, 1.20, 1.50, 2.00):
        for h in (1500.0, 3000.0):
            for lr_true in lrs:
                ba, ea, aod = iso_sr_haze(z, atm, sr, lr_true, scale_height_m=h)
                lad = ladder(z, atm, ba, ea, C_true, options=options)
                j = int(np.argmin(np.abs(z - Z_SR)))
                # Analytic term T at 4.5 km: 2*(alpha_ret - alpha_true) with alpha_ret ~
                # beta_aer * 52 sr (the Klett recovers backscatter, not extinction).
                term_t = 2.0 * float(ba[j]) * (lr_assumed - lr_true) * 1000.0 * 100.0
                scan_iso.append(dict(sr_target=sr, scale_height_m=h, lidar_ratio_true=lr_true,
                                     aod=aod, beta_aer_4500m=float(ba[j]),
                                     term_T_predicted_pct_per_km=term_t,
                                     **{k: lad[k] for k in ("dC_dz_pct_per_km", "dev_at_3km",
                                                            "dev_at_6km", "n_pass_all",
                                                            "r_slope_pct_per_km",
                                                            "dev_pct", "passes")}))
    out["scan_isosr"] = scan_iso

    # ---- E. uncertainty coverage ---------------------------------------------------------------
    cover = []
    for kind, kwargs in (("below", dict(top_m=1500.0)),
                         ("haze", dict(scale_height_m=2000.0))):
        for aod in (0.05, 0.20, 0.50):
            for lr_true in (20.0, 30.0, 52.0, 80.0, 90.0):
                ba, ea = scene(z, kind, aod=aod, lidar_ratio_true=lr_true, **kwargs)
                rcs = forward_signal(z, C_true, atm["beta_mol"], ba, ea)
                res = retrieve(z, rcs, atm, options=options,
                               window=(4500.0 - HALF_M, 4500.0 + HALF_M), perturbations=True)
                if not res["ok"]:
                    continue
                bias = res["C_L_median"] - C_true
                unc = res["C_L_uncertainty"]
                cover.append(dict(scene=kind, aod=aod, lidar_ratio_true=lr_true,
                                  window_centre_m=4500.0,
                                  sr_at_4500m=scattering_ratio(atm, ba, Z_SR),
                                  bias_pct=float(bias / C_true * 100.0),
                                  unc_pct=float(unc / C_true * 100.0),
                                  covered=bool(abs(bias) <= unc),
                                  coverage_ratio=float(abs(bias) / unc) if unc > 0 else np.inf,
                                  r2=res["r_squared"], rel_error_pct=res["rel_error"],
                                  gate_pass=bool(res["r_squared"] >= MIN_R2
                                                 and res["rel_error"] <= MAX_REL_ERROR_PCT)))
    out["uncertainty_coverage"] = cover

    # ---- F. which scenes match BOTH measured observables? --------------------------------------
    matched = [r for r in (scan_d + scan_iso)
               if np.isfinite(r.get("r_slope_pct_per_km", np.nan))
               and abs(r["r_slope_pct_per_km"] - TARGET_R_SLOPE_PCT_PER_KM) < 2.0
               and np.isfinite(r["dC_dz_pct_per_km"])]
    out["matched_to_observed_R"] = sorted(matched, key=lambda r: abs(
        r["dC_dz_pct_per_km"] - TARGET_DCDZ_PCT_PER_KM))

    if verbose:
        _print_summary(out)
    return out, z, atm, C_true, options


def _step_across(lad, base_m, top_m):
    """Mean dev above ``top_m`` minus mean dev below ``base_m``, gate-passing rungs only, in %."""
    c = np.asarray(lad["centres_m"], float)
    d = np.asarray(lad["dev_pct"], float)
    p = np.asarray(lad["passes"], bool)
    lo = p & (c + HALF_M < base_m)
    hi = p & (c - HALF_M > top_m)
    if lo.sum() < 1 or hi.sum() < 1:
        return float("nan")
    return float(np.mean(d[hi]) - np.mean(d[lo]))


def _print_summary(out):
    print("\n" + "=" * 100)
    print("A. LAYER 0.1-1.5 km, ENTIRELY BELOW THE LADDER   (retrieval assumes S = 52 sr)")
    print("    AOD  S_true   bias@3km    bias@6km    dC/dz(3-6 km)   rungs passing gates")
    for r in out["scan_below"]:
        if r["lidar_ratio_true"] in (20.0, 52.0, 90.0):
            print(f"   {r['aod']:.2f}  {r['lidar_ratio_true']:5.0f}  {r['dev_at_3km']:+9.3f} % "
                  f"{r['dev_at_6km']:+10.3f} %   {r['dC_dz_pct_per_km']:+9.5f} %/km      "
                  f"{r['n_pass_all']}/{len(out['centres_m'])}")

    print("\nB. INTRUDING LAYER (base 100 m, top swept), AOD 0.20")
    print("    top   S_true   bias@3km    bias@6km    dC/dz(3-6)   rungs passing")
    for r in out["scan_intruding"]:
        if r["aod"] == 0.20:
            print(f"  {r['top_m']/1000:4.1f}km {r['lidar_ratio_true']:5.0f}  "
                  f"{r['dev_at_3km']:+9.3f} % {r['dev_at_6km']:+10.3f} %  "
                  f"{r['dC_dz_pct_per_km']:+9.5f} %/km   {r['n_pass_all']}/17")

    print("\nC. LOFTED 1 km LAYER INSIDE THE LADDER -- pure transmission (term T)")
    print("   The relevant metric is the STEP across the layer, not a slope: below the layer the")
    print("   ladder is exactly flat and above it exactly flat, at two different levels.")
    print("   base  AOD    S_true  SR in layer   STEP across the layer   worst gate-PASSING rung")
    for r in out["scan_lofted"]:
        if r["aod"] in (0.002, 0.005):
            print(f"  {r['base_m']/1000:.0f}km {r['aod']:.3f} {r['lidar_ratio_true']:6.0f}  "
                  f"{r['sr_in_layer']:9.2f}   {r['step_pct']:+18.3f} %  "
                  f"{r['max_dev_passing_pct']:+18.2f} %")

    print("\nD2. ISO-SR DECOMPOSITION -- backscatter (term B) frozen, lidar ratio swept")
    print("    SR@4.5km  H     S_true   AOD     dC/dz measured   term T predicted   term B "
          "(= dC/dz at S=52)")
    for sr in sorted({r["sr_target"] for r in out["scan_isosr"]}):
        for h in sorted({r["scale_height_m"] for r in out["scan_isosr"]}):
            ref = next((r for r in out["scan_isosr"] if r["sr_target"] == sr
                        and r["scale_height_m"] == h and r["lidar_ratio_true"] == 52.0), None)
            base = ref["dC_dz_pct_per_km"] if ref else np.nan
            for r in out["scan_isosr"]:
                if r["sr_target"] != sr or r["scale_height_m"] != h:
                    continue
                if r["lidar_ratio_true"] not in (20.0, 52.0, 90.0):
                    continue
                print(f"   {sr:6.2f}  {h/1000:.1f}km {r['lidar_ratio_true']:6.0f} "
                      f"{r['aod']:7.4f}  {r['dC_dz_pct_per_km']:+13.4f}    "
                      f"{r['term_T_predicted_pct_per_km']:+13.5f}    {base:+10.4f}")

    print("\nE. UNCERTAINTY COVERAGE, window centred at 4.5 km "
          "(production S+-20/+-10 sr x altitude shift +-200 m)")
    for r in out["uncertainty_coverage"]:
        print(f"  {r['scene']:6s} AOD {r['aod']:.2f} S_true {r['lidar_ratio_true']:5.0f}  "
              f"bias {r['bias_pct']:+9.3f} %   reported unc {r['unc_pct']:8.3f} %   "
              f"|bias|/unc {r['coverage_ratio']:6.2f}  "
              f"{'covered' if r['covered'] else 'NOT COVERED'}"
              f"{'' if r['gate_pass'] else '   [window fails production gates]'}")

    print("\nF. SCENES MATCHING THE MEASURED R-SLOPE "
          f"({TARGET_R_SLOPE_PCT_PER_KM:+.1f} +- 2 %/km), ranked by closeness to the measured "
          f"dC/dz ({TARGET_DCDZ_PCT_PER_KM:+.0f} %/km)")
    print("    scene           AOD_tot  S_true  SR@4.5km  beta_aer@4.5km   AOD(3-6km)   "
          "R-slope    dC/dz")
    for r in out["matched_to_observed_R"][:14]:
        tag = (f"haze H={r['scale_height_m']/1000:.1f}km" if "sr_at_4500m" in r
               else f"isoSR H={r['scale_height_m']/1000:.1f}km")
        sr = r.get("sr_at_4500m", r.get("sr_target"))
        print(f"   {tag:16s} {r['aod']:.4f} {r['lidar_ratio_true']:6.0f}  {sr:8.2f}  "
              f"{r.get('beta_aer_4500m', float('nan')):.3e}     "
              f"{r.get('aod_3_6km', float('nan')):8.4f}  "
              f"{r['r_slope_pct_per_km']:+8.2f}  {r['dC_dz_pct_per_km']:+8.2f} %/km")
    print("=" * 100)


# ================================================================================================
# figures  (altitude on Y, landscape, dpi >= 130)
# ================================================================================================
def _pick(rows, **kw):
    for r in rows:
        if all(abs(r[k] - v) < 1e-9 if isinstance(v, float) else r[k] == v
               for k, v in kw.items()):
            return r
    return None


def _fig_below(out, z, atm, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    centres = np.asarray(out["centres_m"], float) / 1000.0
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 6.4))

    ax = axes[0]
    for aod, col in zip((0.05, 0.20, 0.50), ("C0", "C1", "C3")):
        _, ea = scene(z, "below", aod=aod, lidar_ratio_true=52.0, top_m=1500.0)
        ax.plot(ea * 1000.0, z / 1000.0, color=col, lw=1.7, label=f"AOD = {aod:.2f}")
    ax.plot(atm["alpha_mol"] * 1000.0, z / 1000.0, "k--", lw=1.1,
            label=r"molecular $\alpha_{mol}$")
    ax.axhspan(SLOPE_BAND[0] / 1000, SLOPE_BAND[1] / 1000, color="C2", alpha=0.10)
    ax.text(2e-4, 4.4, "ladder band 3-6 km", color="C2", fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(1e-4, 1e0)
    ax.set_ylim(0, 7)
    ax.set_xlabel(r"aerosol extinction $\alpha_{aer}$  [km$^{-1}$]")
    ax.set_ylabel("range above the instrument  [km AGL]")
    ax.set_title("Scene A: boundary-layer top-hat 0.1-1.5 km\n(entirely below every window)",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[1]
    ax.axvline(0.0, color="C3", lw=1.4, ls="--", label=r"$C_{true}$")
    for aod, col in zip((0.05, 0.20, 0.50), ("C0", "C1", "C3")):
        r = _pick(out["scan_below"], aod=aod, lidar_ratio_true=52.0)
        ax.plot(r["dev_pct"], centres, "o-", ms=4, color=col,
                label=f"AOD {aod:.2f}:  offset {r['dev_at_3km']:+.2f} %,  "
                      f"tilt {r['dC_dz_pct_per_km']:+.5f} %/km")
    ax.set_xlim(-2.0, 1.0)
    ax.set_xlabel(r"$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title(r"$S_{true} = S_{assumed} = 52$ sr" "\nlarge AOD alone: a small offset, no tilt",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="lower left")

    ax = axes[2]
    ax.axvline(0.0, color="C3", lw=1.4, ls="--", label=r"$C_{true}$")
    for lr_true, col in zip((20.0, 30.0, 52.0, 80.0, 90.0), ("C0", "C4", "C2", "C1", "C3")):
        r = _pick(out["scan_below"], aod=0.20, lidar_ratio_true=lr_true)
        ax.plot(r["dev_pct"], centres, "o-", ms=4, color=col,
                label=f"$S_{{true}}$={lr_true:.0f} sr:  {r['dev_at_3km']:+.1f} %,  "
                      f"tilt {r['dC_dz_pct_per_km']:+.5f} %/km")
    ax.set_xlabel(r"$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title("AOD = 0.20, TRUE lidar ratio swept\nhuge bias, still perfectly vertical",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="lower left")

    fig.suptitle("Below-window aerosol CANNOT tilt the calibration: its optical-depth error is "
                 "fully accumulated under the lowest window and stops growing", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _fig_lr_map(out, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    aods = sorted({r["aod"] for r in out["scan_below"]})
    lrs = sorted({r["lidar_ratio_true"] for r in out["scan_below"]})
    bias = np.full((len(lrs), len(aods)), np.nan)
    slope = np.full_like(bias, np.nan)
    for r in out["scan_below"]:
        i, j = lrs.index(r["lidar_ratio_true"]), aods.index(r["aod"])
        bias[i, j] = r["dev_at_3km"]
        slope[i, j] = r["dC_dz_pct_per_km"]

    fig, axes = plt.subplots(1, 2, figsize=(16.0, 6.0))
    for ax, mat, title, unit, fmt in (
            (axes[0], bias, "OFFSET of $C_L$ (identical at every window height)", "%", "+.1f"),
            (axes[1], slope * 1000.0,
             "TILT of $C_L$ between 3 and 6 km  (note: x1000)", "1e-3 %/km", "+.2f")):
        vmax = float(np.nanmax(np.abs(mat))) or 1.0
        im = ax.imshow(mat, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       extent=(-0.5, len(aods) - 0.5, -0.5, len(lrs) - 0.5))
        ax.set_xticks(range(len(aods)))
        ax.set_xticklabels([f"{a:.2f}" for a in aods])
        ax.set_yticks(range(len(lrs)))
        ax.set_yticklabels([f"{v:.0f}" for v in lrs])
        ax.set_xlabel("aerosol optical depth at 1064 nm  [-]")
        ax.set_ylabel(r"TRUE aerosol lidar ratio $S_{true}$  [sr]")
        ax.axhline(lrs.index(52.0), color="k", lw=1.4, ls="--")
        ax.set_title(title, fontsize=11)
        for i in range(len(lrs)):
            for j in range(len(aods)):
                ax.text(j, i, format(mat[i, j], fmt), ha="center", va="center", fontsize=7.5)
        fig.colorbar(im, ax=ax, label=unit)
    fig.suptitle("Scene A (layer 0.1-1.5 km, below the ladder): the lidar-ratio error makes a "
                 "very large OFFSET and a numerically zero TILT\n"
                 "dashed line = the assumed $S$ = 52 sr", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _fig_step_vs_slope(out, z, atm, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    centres = np.asarray(out["centres_m"], float) / 1000.0
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 6.4))

    ax = axes[0]
    for base, col in zip((3000.0, 4000.0, 5000.0), ("C0", "C1", "C3")):
        ba, _ = scene(z, "lofted", aod=0.005, lidar_ratio_true=90.0,
                      base_m=base, thickness_m=1000.0)
        ax.plot(ba, z / 1000.0, color=col, lw=1.7,
                label=f"layer {base/1000:.0f}-{base/1000+1:.0f} km")
    ba_h, _ = scene(z, "haze", aod=0.05, lidar_ratio_true=52.0, scale_height_m=1000.0)
    ax.plot(ba_h, z / 1000.0, color="C2", lw=1.7, ls="-.", label="deep haze, H = 1 km, AOD 0.05")
    ax.plot(atm["beta_mol"], z / 1000.0, "k--", lw=1.2, label=r"$\beta_{mol}$ (1064 nm)")
    ax.set_xscale("log")
    ax.set_xlim(1e-10, 1e-5)
    ax.set_ylim(0, 8)
    ax.set_xlabel(r"aerosol backscatter $\beta_{aer}$  [m$^{-1}$ sr$^{-1}$]")
    ax.set_ylabel("range above the instrument  [km AGL]")
    ax.set_title("Two shapes with the same question:\nlofted slab vs deep haze", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[1]
    ax.axvline(0.0, color="C3", lw=1.3, ls="--")
    for base, col in zip((3000.0, 4000.0, 5000.0), ("C0", "C1", "C3")):
        r = _pick(out["scan_lofted"], base_m=base, aod=0.005, lidar_ratio_true=90.0)
        dev = np.asarray(r["dev_pct"], float)
        ok = np.asarray(r["passes"], bool)
        ax.plot(np.where(ok, dev, np.nan), centres, "o-", ms=4.5, color=col,
                label=f"{base/1000:.0f}-{base/1000+1:.0f} km, step {r['step_pct']:+.2f} %")
        ax.plot(np.where(~ok, dev, np.nan), centres, "o", ms=4.5, mfc="none", color=col)
        ax.axhspan(base / 1000, base / 1000 + 1, color=col, alpha=0.07)
    # Windows sitting IN the slab are 40-110 % wrong; they are the hollow markers and are left
    # deliberately off-scale so the step itself -- the point of the panel -- stays readable.
    ax.set_xlim(-1.2, 1.2)
    ax.set_xlabel(r"$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title("Term T alone (lofted slab, AOD 0.005, $S_{true}$=90 sr)\n"
                 "a STEP of only $-$0.4 % at the layer, flat above and below", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

    ax = axes[2]
    ax.axvline(0.0, color="C3", lw=1.3, ls="--")
    for aod, col in zip((0.02, 0.05, 0.10), ("C0", "C1", "C3")):
        r = _pick(out["scan_haze"], scale_height_m=1000.0, aod=aod, lidar_ratio_true=52.0)
        dev = np.asarray(r["dev_pct"], float)
        ok = np.asarray(r["passes"], bool)
        ax.plot(np.where(ok, dev, np.nan), centres, "o-", ms=4.5, color=col,
                label=f"AOD {aod:.2f}, SR@4.5km {r['sr_at_4500m']:.2f}, "
                      f"{r['dC_dz_pct_per_km']:+.1f} %/km")
        ax.plot(np.where(~ok, dev, np.nan), centres, "o", ms=4.5, mfc="none", color=col)
        # Term B alone predicts C_L(z)/C_true = SR(z): the retrieval calls the window molecular,
        # so every bit of aerosol backscatter at the window height is absorbed into the constant.
        ba, _ = scene(z, "haze", aod=aod, lidar_ratio_true=52.0, scale_height_m=1000.0)
        sr_pred = np.array([scattering_ratio(atm, ba, c * 1000.0) for c in centres])
        ax.plot((sr_pred - 1.0) * 100.0, centres, ls=":", lw=1.6, color=col)
    ax.plot([], [], ls=":", lw=1.6, color="0.3", label=r"prediction $C_L/C_{true}=SR(z)$")
    ax.set_xlim(-5.0, 90.0)
    ax.set_xlabel(r"$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title("Term B (deep haze, H = 1 km, $S_{true}$ = 52 sr = correct)\n"
                 "a SMOOTH SLOPE -- the observed signature", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

    fig.suptitle("The discriminator: mis-corrected TRANSMISSION makes a step at the layer; "
                 "residual BACKSCATTER inside the windows makes a smooth slope", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _fig_isosr(out, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2))

    ax = axes[0]
    for sr, col in zip((1.05, 1.10, 1.20, 1.50, 2.00), ("C0", "C4", "C2", "C1", "C3")):
        rows = sorted([r for r in out["scan_isosr"]
                       if r["sr_target"] == sr and r["scale_height_m"] == 1500.0
                       and np.isfinite(r["dC_dz_pct_per_km"])],
                      key=lambda r: r["lidar_ratio_true"])
        if not rows:
            continue
        ax.plot([r["lidar_ratio_true"] for r in rows],
                [r["dC_dz_pct_per_km"] for r in rows], "o-", color=col,
                label=f"SR(4.5 km) = {sr:.2f}")
    ax.axvline(52.0, color="0.4", lw=1.1, ls=":")
    ax.axhline(TARGET_DCDZ_PCT_PER_KM, color="C3", lw=1.5, ls="--",
               label=f"observed Payerne {TARGET_DCDZ_PCT_PER_KM:.0f} %/km")
    ax.set_xlabel(r"TRUE aerosol lidar ratio $S_{true}$  [sr]   (retrieval assumes 52)")
    ax.set_ylabel(r"$dC_L/dz$ between 3 and 6 km   [%/km]")
    ax.set_title("Backscatter profile FROZEN, lidar ratio swept (H = 1.5 km)\n"
                 "the tilt barely moves: term T is negligible", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

    ax = axes[1]
    srs = sorted({r["sr_target"] for r in out["scan_isosr"]})
    for h, mk, col in ((1500.0, "o", "C0"), (3000.0, "s", "C1")):
        b_vals, t_vals, keep = [], [], []
        for sr in srs:
            ref = _pick(out["scan_isosr"], sr_target=sr, scale_height_m=h,
                        lidar_ratio_true=52.0)
            if ref is None or not np.isfinite(ref["dC_dz_pct_per_km"]):
                continue
            deltas = [abs(r["dC_dz_pct_per_km"] - ref["dC_dz_pct_per_km"])
                      for r in out["scan_isosr"]
                      if r["sr_target"] == sr and r["scale_height_m"] == h
                      and np.isfinite(r["dC_dz_pct_per_km"])]
            keep.append(sr)
            b_vals.append(abs(ref["dC_dz_pct_per_km"]))
            t_vals.append(max(deltas) if deltas else np.nan)
        ax.plot(keep, b_vals, mk + "-", color=col, ms=6,
                label=f"|term B| (in-window backscatter), H = {h/1000:.1f} km")
        ax.plot(keep, t_vals, mk + "--", color=col, ms=6, mfc="none",
                label=f"|term T| (worst $S$ error, 20-90 sr), H = {h/1000:.1f} km")
        for sr, b, t in zip(keep, b_vals, t_vals):
            if np.isfinite(t) and t > 0:
                ax.annotate(f"x{b/t:.0f}", (sr, b), textcoords="offset points",
                            xytext=(4, 5), fontsize=7.5, color=col)
    ax.axhline(abs(TARGET_DCDZ_PCT_PER_KM), color="C3", lw=1.4, ls="--",
               label=f"observed |dC/dz| = {abs(TARGET_DCDZ_PCT_PER_KM):.0f} %/km")
    ax.set_yscale("log")
    ax.set_xlabel("imposed scattering ratio at 4.5 km  [-]")
    ax.set_ylabel(r"contribution to $|dC_L/dz|$, 3-6 km   [%/km]")
    ax.set_title("Term B dominates term T by the annotated factor\n"
                 "even for the most extreme lidar-ratio error the atmosphere allows",
                 fontsize=11)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.5, loc="lower right")

    fig.suptitle("Decomposition of eq. [2]: transmission (term T) versus in-window backscatter "
                 "(term B)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _fig_uncertainty(out, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(15.0, 6.0))
    for ax, kind, title in ((axes[0], "below", "Scene A: layer 0.1-1.5 km, below the ladder"),
                            (axes[1], "haze", "Scene D: deep haze, H = 2 km")):
        rows = [r for r in out["uncertainty_coverage"] if r["scene"] == kind]
        for aod, col in zip(sorted({r["aod"] for r in rows}), ("C0", "C1", "C3")):
            sub = sorted([r for r in rows if r["aod"] == aod],
                         key=lambda r: r["lidar_ratio_true"])
            lr = [r["lidar_ratio_true"] for r in sub]
            ax.plot(lr, [abs(r["bias_pct"]) for r in sub], "o-", color=col,
                    label=f"|bias|, AOD {aod:.2f}")
            ax.plot(lr, [r["unc_pct"] for r in sub], "s--", color=col, mfc="none",
                    label=f"reported unc., AOD {aod:.2f}")
        ax.set_yscale("log")
        ax.set_xlabel(r"TRUE aerosol lidar ratio $S_{true}$  [sr]")
        ax.set_ylabel(r"relative to $C_{true}$   [%]")
        ax.axvline(52.0, color="0.4", lw=1.0, ls=":")
        ax.set_title(title + "\nwindow centred at 4.5 km", fontsize=11)
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7.5, ncol=2, loc="lower left")
    fig.suptitle(r"Reported uncertainty (production's $S\pm$20/$\pm$10 sr $\times$ altitude-shift "
                 "spread) versus the bias actually committed", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _fig_vs_observed(out, z, atm, C_true, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    obs = json.loads(RATIO_JSON.read_text()) if RATIO_JSON.exists() else {}
    grid = np.arange(1500.0, 8001.0, 100.0)
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.4))

    ax = axes[0]
    for name, rec in obs.items():
        thick = "PAYERNE" in name
        for key, col in (("kept by v2", "C0"), ("recovered by v2.2", "C3")):
            if key not in rec:
                continue
            ax.plot(np.asarray(rec[key], float), grid / 1000.0, color=col,
                    lw=2.2 if thick else 1.0, alpha=0.95 if thick else 0.30,
                    label=(f"measured, {key} - {name.split('_')[0]} "
                           f"({rec[key + ' slope %/km']:+.1f} %/km)") if thick else None)
    j3 = int(np.argmin(np.abs(z - 3000.0)))
    for h, aod, lr_true, col in ((1000.0, 0.05, 52.0, "C2"), (1000.0, 0.10, 52.0, "C4"),
                                 (1500.0, 0.05, 52.0, "C1")):
        ba, ea = scene(z, "haze", aod=aod, lidar_ratio_true=lr_true, scale_height_m=h)
        rcs = forward_signal(z, C_true, atm["beta_mol"], ba, ea)
        r = (rcs / z ** 2) / atm["p_mol"]
        row = _pick(out["scan_haze"], scale_height_m=h, aod=aod, lidar_ratio_true=lr_true)
        ax.plot(r / r[j3], z / 1000.0, ls="--", color=col, lw=1.6,
                label=f"simulated haze H={h/1000:.1f} km, AOD {aod:.2f} "
                      f"({row['r_slope_pct_per_km']:+.1f} %/km)")
    ax.set_xlim(0.2, 2.2)
    ax.set_ylim(1.5, 8.0)
    ax.set_xlabel(r"$R(z)$ = signal / $p_{mol}$, normalised at 3 km  [-]")
    ax.set_ylabel("range above the instrument  [km AGL]")
    ax.set_title("Anchor: does a plausible haze reproduce the measured $R(z)$?\n"
                 "thick = measured night medians (thin = other sites)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper right")

    ax = axes[1]
    xs, ys, cs = [], [], []
    for r in out["scan_haze"]:
        if not np.isfinite(r["dC_dz_pct_per_km"]) or not np.isfinite(r["r_slope_pct_per_km"]):
            continue
        xs.append(r["r_slope_pct_per_km"])
        ys.append(r["dC_dz_pct_per_km"])
        cs.append(r["sr_at_4500m"])
    sc = ax.scatter(xs, ys, c=np.log10(np.maximum(cs, 1.0)), cmap="viridis", s=30)
    fig.colorbar(sc, ax=ax, label=r"$\log_{10}$ scattering ratio at 4.5 km")
    ax.axvline(TARGET_R_SLOPE_PCT_PER_KM, color="C3", lw=1.6, ls="--",
               label=f"measured $R$-slope, recovered nights ({TARGET_R_SLOPE_PCT_PER_KM:+.1f} %/km)")
    ax.axvline(TARGET_R_SLOPE_KEPT_PCT_PER_KM, color="C0", lw=1.4, ls=":",
               label=f"measured $R$-slope, kept nights ({TARGET_R_SLOPE_KEPT_PCT_PER_KM:+.1f} %/km)")
    ax.axhline(TARGET_DCDZ_PCT_PER_KM, color="C1", lw=1.6, ls="--",
               label=f"measured $dC_L/dz$ ({TARGET_DCDZ_PCT_PER_KM:+.0f} %/km)")
    ax.set_ylim(-80, 10)
    ax.set_xlabel(r"simulated $R(z)$ slope, 2-6 km   [%/km]")
    ax.set_ylabel(r"simulated $dC_L/dz$, 3-6 km   [%/km]")
    ax.set_title("Both observables at once: the haze scenes that match the measured $R$-slope\n"
                 "also land on the measured $C_L$ tilt", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="lower left")

    fig.suptitle("Anchoring the aerosol mechanism on the measured profiles "
                 "(Payerne / Lindenberg / Gottfrieding CHM15k)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out, z, atm, C_true, options = run_all(verbose=True)
    _fig_below(out, z, atm, FIGDIR / "fwd_trans_below_window_ladder.png")
    _fig_lr_map(out, FIGDIR / "fwd_trans_aod_lidarratio_map.png")
    _fig_step_vs_slope(out, z, atm, FIGDIR / "fwd_trans_step_vs_slope.png")
    _fig_isosr(out, FIGDIR / "fwd_trans_lidarratio_term.png")
    _fig_uncertainty(out, FIGDIR / "fwd_trans_uncertainty_coverage.png")
    _fig_vs_observed(out, z, atm, C_true, FIGDIR / "fwd_trans_vs_observed_ratio.png")
    OUTJSON.write_text(json.dumps(out, indent=1, default=float))
    print(f"\nJSON -> {OUTJSON}")
    print(f"figures -> {FIGDIR}")
    return out


if __name__ == "__main__":
    main()
