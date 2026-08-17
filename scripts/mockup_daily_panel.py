#!/usr/bin/env python3
"""Mockup: ONE daily calibration panel that serves every calibration method.

The per-night diagnostics are currently three unrelated matplotlib figures (Rayleigh success,
Rayleigh rejection, cloud) with different layouts, different axes and no navigation. This mockup
replaces all three with a single client-side panel:

  * top row -- the time-height curtain (2/3 width) with its masks, bands and cloud base, and a
    profile panel (1/3) whose SELECTOR switches between Raw RCS, calibrated signal and attenuated
    backscatter, with a linear/log x toggle, the molecular reference and the calibration area;
  * below   -- the method's diagnostic panels and its messages/flags;
  * left    -- a month calendar, day arrows and a method switch, so a night is one click away.

The point of a single panel is that the *same* code draws both methods. That is only defensible
because the physics lines up: the Rayleigh figure plots ``beta_att = signal_normalized * range**2``
(plotting.py:899), i.e. beta_att = RCS / C_L and signal_normalized = beta_att / range**2, and the
cloud path carries physical ``beta`` = beta_att directly (calibration/cloud/calibration.py
``build_cloud_input``). So the three views are one definition, computed from whichever quantity the
method happens to hold.

Everything is captured from REAL calibrations: the production plotting functions are wrapped, so a
panel can only show arrays the PNG was drawn from. Days where no figure is drawn at all (some
Rayleigh rejections draw nothing) still appear, carrying their flag and message -- the image
product cannot do that.

  python scripts/mockup_daily_panel.py --key 0-20000-0-06610_C \
      --start 20260701 --end 20260709 --l1-root D:/E-PROFILE_L1_2026 --out mock_daily.html
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# The curtain encoder is shared with the feasibility mockup rather than copied: block-averaging in
# linear space and the uint8 quantisation are subtle enough that two divergent copies would rot.
from scripts.mockup_dynamic_diag import (  # noqa: E402
    _curtain, _f32, _jsonable, _runs,
)
import scripts.mockup_dynamic_diag as DYN  # noqa: E402

DYN.HIRES = False          # block averaging measured decisively better than a native-resolution
                           # opt-in (42.5 % of subsampled cells were falsely empty), so the second
                           # payload is gone: there is one resolution and it is the averaged one.

CAPTURED: list = []
SCALE = 1e6                # Mm^-1, exactly the PNG's own axis scaling


# ================================================================== capture of the real renders
_CAPTURE_INSTALLED = False


def _install_capture() -> None:
    """Wrap the three production plotting functions and keep their arguments.

    IDEMPOTENT, and that matters: build_payload() calls this every time, so a batch run that builds
    900 payloads used to wrap the already-wrapped functions 900 times over. Every plot call then
    descended 900 nested wrappers and appended 900 entries to CAPTURED -- the results stayed correct
    (the reader takes the last entry) but the cost grew linearly through the run.
    """
    global _CAPTURE_INSTALLED
    if _CAPTURE_INSTALLED:
        return
    _CAPTURE_INSTALLED = True
    import calibration.plotting as P
    import calibration.rayleigh.calibration as RC

    def wrap(real, kind, is_cloud=False):
        def w(*a, **kw):
            fig = real(*a, **kw)
            CAPTURED.append({"kind": kind, "args": a, "kwargs": kw})
            return fig
        return w

    # The window selector is wrapped too, so the "why this window" card quotes the ACTUAL grid and
    # the ACTUAL parameters the run used rather than a description of the algorithm. eprof_v2/v2.2
    # call _select_optimal directly as a module global (molecular_methods.py:865), so rebinding the
    # attribute is enough -- the _SELECTORS dict is not on that path.
    import calibration.rayleigh.molecular_methods as MM
    real_sel = MM._select_optimal

    def sel(g, **kw):
        mw = real_sel(g, **kw)
        CAPTURED.append({"kind": "window", "grid": g, "mw": mw, "params": kw})
        return mw
    MM._select_optimal = sel

    P.plot_rayleigh_diagnostics_compact = wrap(P.plot_rayleigh_diagnostics_compact, "rayleigh_ok")
    RC.plot_rayleigh_diagnostics_compact = P.plot_rayleigh_diagnostics_compact
    P.plot_rayleigh_diagnostics_failure = wrap(P.plot_rayleigh_diagnostics_failure, "rayleigh_fail")
    RC.plot_rayleigh_diagnostics_failure = P.plot_rayleigh_diagnostics_failure
    P.plot_cloud_diagnostics_compact = wrap(P.plot_cloud_diagnostics_compact, "cloud", True)


# ============================================================================= profile builders
def _views(rcs, sig, beta, mol_sig, mol_beta) -> dict:
    """The three profile views, in the PNG's own units.

    ``rcs`` is instrument-native and has no molecular counterpart (nothing to compare an
    uncalibrated signal to) -- that absence is meaningful and is preserved rather than papered over
    with a rescaled curve.
    """
    out = {}
    if rcs is not None and np.isfinite(rcs).any():
        out["rcs"] = {"x": _f32(rcs), "label": "Raw RCS (a.u.)", "logx": False, "mol": None}
    if sig is not None and np.isfinite(sig).any():
        out["sig"] = {"x": _f32(sig * SCALE), "label": "Signal (Mm⁻¹)", "logx": True,
                      "mol": _f32(mol_sig * SCALE) if mol_sig is not None else None}
    if beta is not None and np.isfinite(beta).any():
        # LOG by default, and this is a deliberate departure from the PNG, whose beta_att panel is
        # linear (plotting.py:900-908) while its cloud counterpart is log. beta_att spans orders of
        # magnitude between the molecular tail and a cloud, so linear collapses everything below the
        # peak into the axis; the operator can still flip it. The cloud figure had it right.
        out["beta"] = {"x": _f32(beta * SCALE), "label": "β_att (Mm⁻¹ sr⁻¹)", "logx": True,
                       "mol": _f32(mol_beta * SCALE) if mol_beta is not None else None}
    return out


def _profile_rayleigh(kw: dict) -> dict:
    """Rayleigh: all three views exist as first-class arrays on the fit grid."""
    rng = np.asarray(kw.get("range_alc", []), float)
    alt = float(kw.get("altitude", 0.0) or 0.0)
    z_km = rng * 1e-3                       # AGL, matching the curtain and the cloud method
    sig = np.asarray(kw.get("signal_normalized", []), float)
    p_mol = np.asarray(kw.get("p_mol", []), float)
    b_mol = np.asarray(kw.get("beta_att_mol", []), float)
    beta = sig * rng ** 2 if sig.size == rng.size else None   # derived in the plot, never passed
    # fit_altitude_start/end are ASL (the PNG draws them on an ASL axis); subtract the station
    # altitude so the band lands on the same gates as the curtain's fitted layer.
    z0 = (float(kw.get("fit_altitude_start", 0.0) or 0.0) - alt) * 1e-3
    z1 = (float(kw.get("fit_altitude_end", 0.0) or 0.0) - alt) * 1e-3
    out = {
        "y": _f32(z_km), "y_label": "Range (km AGL)",
        "y_max": float(min(z1 + 3.0, float(z_km.max()) if z_km.size else z1 + 3.0)),
        "views": _views(np.asarray(kw.get("rcs_mean", []), float), sig, beta, p_mol, b_mol),
        "cal_band": [z0, z1], "cal_band_label": "Rayleigh fit window",
    }
    # These curves are ALREADY unflagged-only: rcs_mean is the aggregate of rcs_use =
    # data.rcs[keep_idx] (calibration.py:810/826), and used_profile_indices IS that same keep_idx
    # (calibration.py:1288). Nothing to filter here -- but say so, because "mean profile" on its own
    # does not tell the operator whether the screened profiles are in it.
    upi, rcs_m = kw.get("used_profile_indices"), kw.get("rcs")
    if upi is not None and rcs_m is not None:
        out["mean_of"] = "mean of %d/%d unflagged profiles" % (int(np.size(upi)),
                                                               int(np.asarray(rcs_m).shape[0]))
    return out


def _profile_rayleigh_fail(kw: dict) -> dict:
    """Rejection: only the night-mean signal exists, and the molecular curve is SCALED to it.

    The failure figure fits nothing, so there is no calibration constant and no absolute
    normalisation. Its molecular curve is matched to the observed signal over a fixed range band
    purely so the two shapes can be compared by eye (plotting.py:1053-1162) -- shipping the scaled
    product rather than p_mol keeps that honest, because the scale factor is not a calibration.
    """
    rng = np.asarray(kw.get("range_alc", []), float)
    z_km = rng * 1e-3                       # AGL everywhere on this page
    m = kw.get("rcs")
    if m is None:
        return {}
    m = np.asarray(m, float)
    # ONLY UNFLAGGED PROFILES. The failure PNG averages the entire night, including the very cloudy
    # profiles that caused the rejection -- which is most of why its "mean profile" looks nothing
    # like molecular and reads as a worse night than it was. The same low-cloud test the success
    # path uses is available in these kwargs (cbh + z_low_cloud), so apply it here too.
    n_all = int(m.shape[0])
    keep = np.ones(n_all, bool)
    note = "mean of all %d profiles (no screen available)" % n_all
    cbh, z_cut = kw.get("cbh"), kw.get("z_low_cloud")
    if cbh is not None and z_cut is not None and np.isfinite(float(z_cut)):
        a = np.asarray(cbh, float)
        a = a[:, 0] if a.ndim > 1 else a
        ncv = float(kw.get("no_cloud_value", -9.0))
        a = np.where((a == ncv) | (a <= 0), np.nan, a)
        cand = ~(np.isfinite(a) & (a < float(z_cut)))
        if cand.any():                      # an all-cloud night has nothing to screen down to
            keep = cand
            note = "mean of %d/%d unflagged profiles" % (int(keep.sum()), n_all)
        else:
            note = "mean of all %d profiles (every one is low-cloud flagged)" % n_all
    with np.errstate(all="ignore"):
        rcs_mean = np.nanmean(m[keep], axis=0)
        sig = rcs_mean / (rng ** 2)
    finite = np.isfinite(sig)
    if not finite.any():
        return {}
    # NOT "calibrated signal": nothing was calibrated on this night. The selector shows this
    # view's own name so the label cannot promise a constant that does not exist.
    views = {"sig": {"x": _f32(np.where(finite, sig, np.nan)), "label": "Signal / range² (a.u.)",
                     "logx": True, "mol": None, "name": "Signal (uncalibrated)"},
             "rcs": {"x": _f32(rcs_mean), "label": "Raw RCS (a.u.)", "logx": False, "mol": None}}
    pm = kw.get("p_mol")
    if pm is not None:
        pm = np.asarray(pm, float)
        r0 = float(kw.get("range_start_m", 2000.0))
        r1 = float(kw.get("range_end_m", 6000.0))
        band = (rng >= r0) & (rng <= r1) & finite & np.isfinite(pm)
        if band.any() and np.nanmedian(pm[band]) > 0:
            scale = float(np.nanmedian(sig[band]) / np.nanmedian(pm[band]))
            views["sig"]["mol"] = _f32(pm * scale)
            views["sig"]["mol_name"] = "molecular (shape-matched, not calibrated)"
    win = kw.get("molecular_window")
    band_km = None
    if win is not None and np.all(np.isfinite(np.asarray(win, float))):
        w = np.asarray(win, float)          # molecular_window is in RANGE metres already
        band_km = [float(w[0]) * 1e-3, float(w[1]) * 1e-3]
    return {"y": _f32(z_km), "y_label": "Range (km AGL)",
            "y_max": float(z_km.max()) if np.isfinite(z_km).any() else 1.0,
            "views": views, "cal_band": band_km, "mean_of": note,
            "cal_band_label": "molecular window (attempted)"}


def _profile_cloud(data, res, sel: np.ndarray, cal_lo: float, cal_hi: float) -> dict:
    """Cloud: everything derives from the physical beta the method works on.

    beta IS the attenuated backscatter, so the other two views are exact algebra rather than a
    second retrieval: signal = beta / range**2, and Raw RCS = beta * C_applied (the constant the
    loader divided out). On the L1 path that product is the raw RCS up to a fixed unit factor,
    which is why the axis says a.u. -- the SHAPE, which is what this panel is read for, is exact.
    """
    rng = np.asarray(data.range, float)
    beta_2d = np.asarray(data.beta, float)                 # (n_range, n_time)
    use = sel if sel.any() else np.ones(beta_2d.shape[1], bool)
    with np.errstate(all="ignore"):
        prof = np.nanmean(beta_2d[:, use], axis=1)
        sig = np.where(rng > 0, prof / np.where(rng > 0, rng ** 2, np.nan), np.nan)
    cc = getattr(data, "calibration_constant_applied", None)
    rcs = prof * float(cc) if (cc and np.isfinite(cc)) else None
    y_max = float(min(cal_hi * 1e-3 + 1.0, float(rng.max()) * 1e-3))
    cbh = np.asarray(getattr(data, "cbh", []), float)
    med_cbh = (float(np.nanmedian(cbh[sel])) if (sel.any() and cbh.size == beta_2d.shape[1]
                                                 and np.isfinite(cbh[sel]).any()) else float("nan"))
    views = _views(rcs, sig, prof, None, None)   # beta already defaults to log, as this method's
                                                 # own figure does (plotting.py:176)
    out = {"y": _f32(rng * 1e-3), "y_label": "Range (km AGL)", "y_max": y_max,
           "views": views,
           "cal_band": [cal_lo * 1e-3, cal_hi * 1e-3],
           "cal_band_label": "O'Connor integration gate (100–2400 m)",
           "n_used": int(sel.sum()),
           "mean_of": ("mean of %d selected profiles" % int(sel.sum())) if sel.any()
                      else "night mean (no profile was selected)"}
    if np.isfinite(med_cbh):
        # Only the median cloud base, drawn as a line and labelled as the summary statistic it is.
        # The PNG also shades a "where B accumulates (CBH -> +300 m)" band there; that band is a
        # PER-PROFILE truth and this curve is a mean over profiles with DIFFERENT cloud bases, so
        # the accumulation layer is already smeared across the average. Shading one 300 m slab at
        # the median asserts a coherence the averaging destroyed, so it is not drawn here.
        out["med_cbh_km"] = med_cbh * 1e-3
        out["med_cbh_label"] = "median cloud base of the selected profiles"
    return out


# ============================================================================= curtain builders
def _iso(tdt, step: int) -> list | None:
    if tdt is None or not len(tdt):
        return None
    return [str(t)[:19] for t in np.asarray(tdt)[::step]]


def _curtain_rayleigh(kw: dict, failed: bool) -> dict:
    """Rayleigh curtain on RANGE AGL -- the single vertical coordinate used across the whole page.

    The PNG mixes the two: its curtain is range AGL while its molecular panels are altitude ASL, so
    the same feature sits at two different heights in one figure (491 m apart at Payerne). AGL is the
    side to converge on -- it is the instrument's own coordinate, it is what the cloud method already
    uses end to end, and it is what the fit window is searched over (best_range_m). The profile's
    fit band is therefore converted OUT of ASL, in _profile_rayleigh.
    """
    rcs = kw.get("rcs")
    if rcs is None:
        return {}
    m = np.asarray(rcs, float)
    rng = np.asarray(kw.get("range_alc", []), float)
    rr = np.asarray(kw.get("rcs_range_alc", rng), float)
    if failed:
        rr = rng
    c = _curtain(m, rr)
    st = c["st_t"]
    c["time"] = _iso(kw.get("time_datetime"), st)
    if c["time"] is None:
        c["hours"] = _jsonable(np.asarray(kw.get("hours_since_start", []), float)[::st])
    c["y_label"] = "Range (km AGL)"
    c["y_max_km"] = float(rr.max()) * 1e-3 if rr.size else None
    # The curtain holds RAW RCS, and the view selector drives BOTH panels, so the page has to be
    # able to turn it into the other two quantities. Both are exact algebra on the stored log10
    # field -- beta_att = RCS/C_L is a constant offset, the calibrated signal divides by range**2 on
    # top of that -- so one curtain serves all three views and nothing extra is shipped. Without a
    # constant (a rejected night) only the raw view is defined, and log_c stays null.
    c["base"] = "rcs"
    cl = kw.get("cl_median")
    c["log_c"] = (float(np.log10(float(cl))) if (cl and np.isfinite(float(cl)) and float(cl) > 0)
                  else None)

    n_t = int(m.shape[0])
    cbh0 = None
    key = "cbh" if failed else "cloud_base_height"
    if kw.get(key) is not None:
        a = np.asarray(kw[key], float)
        a = a[:, 0] if a.ndim > 1 else a
        ncv = float(kw.get("no_cloud_value", -9.0))
        cbh0 = np.where((a == ncv) | (a <= 0), np.nan, a)      # mask FIRST, then stride
        if np.isfinite(cbh0).any():
            c["cbh"] = [{"y": _jsonable(cbh0[::st] * 1e-3), "name": "cloud base",
                         "color": "#ffffff", "size": 5}]

    c["bands"] = []
    br, bh = kw.get("best_range_m"), kw.get("best_half_m")
    if not failed:
        flagged = np.zeros(n_t, bool)
        if cbh0 is not None:
            z_cut = kw.get("z_low_cloud")
            if z_cut is None:
                z_cut = float(br) - float(bh) if (br is not None and bh is not None) else 0.0
            flagged = np.isfinite(cbh0) & (cbh0 < float(z_cut))
        used = np.zeros(n_t, bool)
        upi = kw.get("used_profile_indices")
        if upi is not None and np.size(upi) > 0:
            used[np.asarray(upi, int)] = True
        nu = ~used
        # Both masks are the SAME screen (calibration.py:797-809: profiles more than ~4 robust
        # sigma from the night's median range-normalised signal). They are split by whether the
        # profile also carries a cloud base below the fit window, so the names say which of the two
        # it is rather than "screened / not used", which described neither.
        c["bands"] = [
            {"runs": _runs(nu[::st] & flagged[::st]), "color": "rgba(214,39,40,0.42)",
             "name": "excluded — low cloud below the window"},
            {"runs": _runs(nu[::st] & ~flagged[::st]), "color": "rgba(90,90,90,0.40)",
             "name": "excluded — signal outlier, no low cloud"},
        ]
        if cbh0 is not None:
            hc = used & np.isfinite(cbh0) & ~flagged
            y_hi = float(rr.max()) * 1e-3
            # EXACTLY the PNG's expression (plotting.py:1030): outside the mask y_lo == y_hi, a
            # ZERO-height band. Emitting null there instead looks equivalent but is not --
            # fill:'tonexty' bridges across nulls and welds distant points into spurious wedges.
            c["high_cloud_ylo"] = [round(float(np.clip((v - 500.0) * 1e-3, 0.0, y_hi)), 4)
                                   if (h and np.isfinite(v)) else round(y_hi, 4)
                                   for h, v in zip(hc[::st], cbh0[::st])]
            c["has_high_cloud"] = bool(hc.any())
            c["y_hi_km"] = y_hi
        if br is not None and bh is not None:
            c["layer_km"] = [(float(br) - float(bh)) * 1e-3, (float(br) + float(bh)) * 1e-3]
            c["layer_name"] = "molecular layer (fitted)"
    else:
        win = kw.get("molecular_window")
        if win is not None and np.all(np.isfinite(np.asarray(win, float))):
            w = np.asarray(win, float)
            c["layer_km"] = [float(w[0]) * 1e-3, float(w[1]) * 1e-3]
            c["layer_name"] = "molecular window (attempted)"
    return c


def _curtain_cloud(data, sel: np.ndarray, cal_lo: float, cal_hi: float, y_max: float) -> dict:
    rng = np.asarray(data.range, float)
    beta = np.asarray(data.beta, float)                    # (n_range, n_time)
    # This panel's y axis is fixed at ~3.4 km while the mesh is full-resolution, so every gate
    # above is drawn and then clipped. Cropping there is free fidelity, and it happens AFTER the
    # colour percentiles are taken so the colours still match the PNG.
    c = _curtain(beta.T, rng, crop_km=y_max)
    st = c["st_t"]
    c["time"] = _iso(np.asarray(data.time), st)
    if c["time"] is None:
        c["hours"] = _jsonable(np.arange(beta.shape[1], dtype=float)[::st])
    c["y_label"] = "Range (km AGL)"
    c["y_max_km"] = y_max
    # Mirror image of the Rayleigh case: the cloud path's curtain is already physical beta_att, so
    # the raw view multiplies by the applied constant instead of dividing.
    c["base"] = "beta"
    cc = getattr(data, "calibration_constant_applied", None)
    c["log_c"] = (float(np.log10(float(cc))) if (cc and np.isfinite(float(cc)) and float(cc) > 0)
                  else None)
    c["bands"] = [{"runs": _runs(sel[::st]), "color": "rgba(44,160,44,0.35)",
                   "name": "used for calibration"}]
    c["layer_km"] = [cal_lo * 1e-3, cal_hi * 1e-3]
    c["layer_name"] = "integration gate"
    cbh = np.asarray(getattr(data, "cbh", []), float)
    if cbh.size == beta.shape[1] and np.isfinite(cbh).any():
        ok = cbh > 0
        c["cbh"] = [
            {"y": _jsonable(np.where(ok & ~sel, cbh, np.nan)[::st] * 1e-3), "name": "cloud base",
             "color": "#ffffff", "size": 5},
            {"y": _jsonable(np.where(ok & sel, cbh, np.nan)[::st] * 1e-3),
             "name": "cloud base (used)", "color": "#ff2d95", "size": 8},
        ]
    return c


# ================================================================================ per-day payload
def _share_y(cur: dict, prof: dict) -> None:
    """Put the curtain and the profile on ONE vertical axis.

    They must agree in QUANTITY (range AGL, for every method and every panel -- the builders take
    care of that) and in LIMITS, or a feature seen in the curtain cannot be read across to the
    profile at the same height. The shared limit is the CURTAIN's full extent, not the
    profile's tighter fit window: cropping the curtain would hide precisely the high cloud the
    screening reacts to, whereas an over-tall profile costs nothing but white space -- and the page
    zoom-links the two, so one drag narrows both.
    """
    y = cur.get("y_max_km")
    if y and prof:
        prof["y_max"] = float(y)


def _why_window(cap: dict) -> dict:
    """Reconstruct WHY the selector chose this window, from the grid it actually scored.

    `_select_optimal` maximises a composite quality (molecular_methods.py:679-688)

        Q = R² − w_ratio·|SR−1| − w_resid·resid − w_snr·σ_ratio − w_tvar·CV_t − w_rel·relerr
            + w_npts·n̂

    over the windows that passed the gates. Recomputing each term at the winning cell is what turns
    "the algorithm picked 3.8-7.1 km" into "it picked it because the temporal CV cost it 0.06 and
    everything else was nearly free" -- and the weights come from the captured call, not from this
    file's idea of the defaults, so a re-tuned deployment cannot make this card lie.
    """
    g, mw, prm = cap["grid"], cap["mw"], cap["params"]
    ci = int(np.argmin(np.abs(np.asarray(g.center_m, float) - float(mw.center_m))))
    hj = int(np.argmin(np.abs(np.asarray(g.half_m, float) - float(mw.half_m))))
    at = lambda a: (float(np.asarray(a, float)[ci, hj])                          # noqa: E731
                    if np.ndim(a) == 2 else float("nan"))
    w = {k: float(prm.get(k, d)) for k, d in
         (("w_ratio", 0.25), ("w_resid", 0.20), ("w_snr", 0.10),
          ("w_tvar", 0.35), ("w_rel", 0.20), ("w_npts", 0.10))}

    have_tvar = bool(np.any(np.isfinite(g.temporal_cv)))
    npts = np.asarray(g.n_pts, float)
    n_norm = at(npts) / (float(np.nanmax(npts)) if np.nanmax(npts) > 0 else 1.0)
    cv = at(g.temporal_cv) if have_tvar else 0.0
    cv = 0.0 if not np.isfinite(cv) else cv
    rel = at(g.rel_error)
    rel = 0.0 if not np.isfinite(rel) else rel
    sr, resid, rstd, r2 = at(g.scattering_ratio), at(g.residual_pct), at(g.ratio_std), at(g.r2)

    terms = [
        {"name": "R² of the molecular fit", "sym": "R²", "value": f"{r2:.4f}",
         "contrib": r2, "why": "how straight signal-vs-molecular is in this window"},
        {"name": "aerosol load", "sym": "−w·|SR−1|", "value": f"{sr:.3f}",
         "contrib": -w["w_ratio"] * abs(sr - 1.0),
         "why": "scattering ratio; 1.00 is aerosol-free"},
        {"name": "Rayleigh-shape residual", "sym": "−w·resid", "value": f"{resid:.1f} %",
         "contrib": -w["w_resid"] * resid / 100.0,
         "why": "relative RMSE of the forced-through-zero fit"},
        {"name": "in-window scatter", "sym": "−w·σ_ratio", "value": f"{rstd:.3f}",
         "contrib": -w["w_snr"] * rstd, "why": "SNR proxy: spread of signal/molecular"},
        {"name": "temporal variability", "sym": "−w·CV_t", "value":
            ("—" if not have_tvar else f"{cv:.3f}"), "contrib": -w["w_tvar"] * cv,
         "why": "molecular is steady through the night; aerosol advects"},
        {"name": "curvature", "sym": "−w·relerr", "value": f"{rel:.1f} %",
         "contrib": -w["w_rel"] * rel / 100.0,
         "why": "|slope − median ratio|, a spatial aerosol signature"},
        {"name": "window length", "sym": "+w·n̂", "value": f"{at(npts):.0f} gates",
         "contrib": w["w_npts"] * n_norm, "why": "longer windows are rewarded, mildly"},
    ]
    total = sum(t["contrib"] for t in terms)
    # The driver is the biggest PENALTY, not the biggest term: R² is always the largest number and
    # naming it every time would say nothing about this night.
    pen = [t for t in terms if t["contrib"] < 0]
    driver = min(pen, key=lambda t: t["contrib"])["name"] if pen else None

    elig = np.asarray(getattr(mw, "eligible", np.zeros_like(g.r2, bool)), bool)
    gates = [("window start", f"≥ {prm.get('min_window_start_m', 2000.0):.0f} m",
              f"{float(mw.start_m):.0f} m"),
             ("R²", f"≥ {prm.get('min_r2', 0.5):.2f}", f"{r2:.3f}"),
             ("Rayleigh residual", f"≤ {prm.get('max_residual_pct', 12.0):.0f} %",
              f"{resid:.1f} %"),
             ("scattering ratio", f"≤ {prm.get('max_scattering_ratio', 1.1):.2f}", f"{sr:.3f}"),
             ("in-window scatter", f"≤ {prm.get('max_ratio_std', 0.30):.2f}", f"{rstd:.3f}"),
             ("temporal CV", f"≤ {prm.get('max_temporal_cv', 0.5):.2f}",
              "—" if not have_tvar else f"{cv:.3f}"),
             ("curvature", f"≤ {prm.get('max_rel_error', 15.0):.0f} %", f"{rel:.1f} %")]
    return {
        "message": str(getattr(mw, "message", "") or ""),
        "window_km": [float(mw.start_m) * 1e-3, float(mw.end_m) * 1e-3],
        "centre_km": float(mw.center_m) * 1e-3, "half_km": float(mw.half_m) * 1e-3,
        "terms": [{**t, "contrib": round(t["contrib"], 4)} for t in terms],
        "total": round(total, 4), "driver": driver,
        "n_eligible": int(elig.sum()), "n_total": int(elig.size),
        "gates": [{"name": a, "limit": b, "value": c} for a, b, c in gates],
        "weights": w,
    }


def _pack_rayleigh_ok(kw: dict) -> dict:
    cur, prof = _curtain_rayleigh(kw, False), _profile_rayleigh(kw)
    _share_y(cur, prof)
    out = {"kind": "rayleigh_ok", "curtain": cur, "profile": prof, "diag": {}}
    grids = {}
    for name in ("slopes", "intercepts", "r_squared"):
        v = kw.get(name)
        if v is None:
            continue
        a = np.abs(np.asarray(v, float)) if name == "intercepts" else np.asarray(v, float)
        grids[name] = {"b64": _f32(a), "shape": list(a.shape)}
    out["diag"]["grids"] = grids
    for name, key in (("range_bin_m", "range_km"), ("half_length_m", "half_km")):
        if kw.get(name) is not None:
            out["diag"][key] = _f32(np.asarray(kw[name], float) * 1e-3)
    br, bh = kw.get("best_range_m"), kw.get("best_half_m")
    if br is not None and bh is not None:
        out["diag"]["best_km"] = [float(br) * 1e-3, float(bh) * 1e-3]
    cl_m, cl_med = kw.get("cl_matrix"), kw.get("cl_median")
    if cl_m is not None and cl_med:
        cm = np.asarray(cl_m, float)
        rel = (cm - float(cl_med)) / float(cl_med) * 100.0
        out["diag"]["sens"] = {
            "b64": _f32(rel), "shape": list(rel.shape),
            "vabs": float(max(np.nanmax(np.abs(rel)) if np.isfinite(rel).any() else 1.0, 1.0)),
            "lr": _jsonable(kw.get("lr_values", [])), "shift": _jsonable(kw.get("alt_shifts", []))}
        out["diag"]["spread"] = {"vals": _jsonable(cm[np.isfinite(cm)].ravel()),
                                 "median": float(cl_med),
                                 "unc": float(kw.get("cl_uncertainty", 0.0) or 0.0)}
    return out


def _pack_rayleigh_fail(kw: dict) -> dict:
    cur, prof = _curtain_rayleigh(kw, True), _profile_rayleigh_fail(kw)
    _share_y(cur, prof)
    return {"kind": "rayleigh_fail", "reason": str(kw.get("reason", "")),
            "curtain": cur, "profile": prof, "diag": {}}


def _pack_cloud(data, res) -> dict:
    n_time = int(np.asarray(data.time).size)
    g = lambda n, d=None: getattr(res, n, d)                                    # noqa: E731
    S_app = np.asarray(g("S_apparent"), float) if g("S_apparent") is not None else np.full(n_time, np.nan)
    S_con = np.asarray(g("S_consistent"), float) if g("S_consistent") is not None else np.full(n_time, np.nan)
    coeffs = np.asarray(g("all_coefficients"), float) if g("all_coefficients") is not None else np.full(n_time, np.nan)
    # A REJECTED scene runs the same code path with an empty result, so `sel` is all-False and every
    # green element legitimately disappears -- that absence IS the rejection signal.
    sel = np.isfinite(S_con) if S_con.size == n_time else np.zeros(n_time, bool)
    cfg = g("config")
    cal_lo = float(getattr(cfg, "cal_minheight", 100.0)) if cfg else 100.0
    cal_hi = float(getattr(cfg, "cal_maxheight", 2400.0)) if cfg else 2400.0
    rng = np.asarray(data.range, float)
    y_max = float(min(cal_hi * 1e-3 + 1.0, float(rng.max()) * 1e-3))
    cal_med = g("cal_median", float("nan"))
    s_theo = (float(np.nanmedian(S_con[sel])) / cal_med) if (sel.any() and cal_med) else 18.8

    cur = _curtain_cloud(data, sel, cal_lo, cal_hi, y_max)
    st = cur["st_t"]
    upper = s_theo * 3.0
    if sel.any() and np.any(np.isfinite(S_con[sel])):
        upper = max(upper, 1.3 * float(np.nanmax(S_con[sel])))
    cal_std = g("cal_std", float("nan"))
    diag = {"x": cur.get("time") or cur.get("hours"),
            "S_app": _jsonable(S_app[::st]), "S_con": _jsonable(np.where(sel, S_con, np.nan)[::st]),
            "s_theo": float(s_theo), "s_ymax": float(upper),
            "coeffs": _jsonable(coeffs[np.isfinite(coeffs)]),
            "cal_median": None if not np.isfinite(cal_med) else float(cal_med),
            "cal_std": None if not np.isfinite(cal_std) else float(cal_std)}
    # A prose description of WHY, in place of the stats dump -- the same treatment the Rayleigh
    # window gets. The three filters are a funnel, so the useful question is which stage removed
    # the profiles, not what each counter reads.
    def _stage(label, what, d):
        det = {str(k): int(v) for k, v in (d or {}).items()
               if isinstance(v, (int, float)) and v} if isinstance(d, dict) else {}
        return {"name": label, "what": what, "removed": int(sum(det.values())), "detail": det}

    stages = [_stage("instrument health", "window transmission, laser energy, quality flag",
                     g("filter_stats")),
              _stage("cloud scene", "peak sharpness, ±300 m around the peak, aerosol ratio, "
                                    "cloud-base range", g("cloud_stats")),
              _stage("temporal consistency", "N consecutive profiles within ±X % of their mean",
                     g("consistency_stats"))]
    remaining, n_used = n_time, int(sel.sum())
    for s in stages:
        s["before"] = remaining
        remaining = max(remaining - s["removed"], 0)
        s["after"] = remaining
    worst = max(stages, key=lambda s: s["removed"]) if stages else None
    s_med = float(np.nanmedian(S_con[sel])) if sel.any() else float("nan")
    diag["why"] = {
        "kind": "cloud",
        "n_total": n_time, "n_used": n_used,
        "stages": stages,
        "driver": (worst["name"] if worst and worst["removed"] else None),
        "gate_km": [cal_lo * 1e-3, cal_hi * 1e-3],
        "s_theo": float(s_theo),
        "s_med": None if not np.isfinite(s_med) else s_med,
        "C": None if not np.isfinite(cal_med) else float(cal_med),
        "C_std": None if not np.isfinite(cal_std) else float(cal_std),
        "C_L": (float(g("lidar_constant")) if np.isfinite(g("lidar_constant", float("nan")))
                else None),
        "applied": (float(getattr(data, "calibration_constant_applied", float("nan")))
                    if getattr(data, "calibration_constant_applied", None) else None),
        "wv": bool(getattr(data, "trans2_wv", None) is not None),
    }
    prof = _profile_cloud(data, res, sel, cal_lo, cal_hi)
    _share_y(cur, prof)          # already agree for cloud (both AGL, both capped at cal_hi + 1 km)
    return {"kind": "cloud", "curtain": cur, "profile": prof, "diag": diag}


# ============================================================================================ run
def _flag_label(flag) -> str | None:
    try:
        from calibration.flags import FLAG_MEANINGS
        return FLAG_MEANINGS.get(float(flag))
    except Exception:                                                           # noqa: BLE001
        return None


def _run_one(rnc, stream: dict, d: datetime, method: str) -> dict:
    """One (day, method): calibrate in-process and pack whatever was drawn."""
    CAPTURED.clear()
    try:
        rows = rnc._do_cloud(stream, d, d) if method == "cloud" else rnc._do_rayleigh(stream, d, d)
    except Exception as exc:                                                    # noqa: BLE001
        return {"kind": "error", "message": f"{type(exc).__name__}: {exc}"}
    ds = d.strftime("%Y%m%d")
    meta = {}
    for r in (rows or []):
        if str(r.get("date")) == ds:
            # A rejection writes the SENTINEL cal_value = -1 (run_network_calibration.py:352/406/470),
            # not a null. So "was this night calibrated?" is decided HERE, on the constant, and is
            # never inferred downstream from which figure happened to be drawn -- the two disagree
            # (flag 0.5 draws the success figure and does produce a constant; a cloud rejection
            # draws the full cloud figure and does not).
            try:
                cv = float(r.get("cal_value"))
            except (TypeError, ValueError):
                cv = float("nan")
            meta = {"flag": r.get("flag"), "message": r.get("message"),
                    "constant": cv if (np.isfinite(cv) and cv > 0) else None,
                    "uncertainty": r.get("uncertainty")}
            meta["flag_label"] = _flag_label(meta.get("flag"))
            break
    plots = [c for c in CAPTURED if c["kind"] != "window"]
    wins = [c for c in CAPTURED if c["kind"] == "window"]
    if not plots:
        # Some Rayleigh rejections (-4, -10) draw no figure at all. The operator still wants the
        # verdict, so the day exists in the panel with its flag and message and no plots.
        out = {"kind": "none"}
        out.update(meta)
        return out
    c = plots[-1]
    if c["kind"] == "cloud":
        a = c["args"]
        out = _pack_cloud(a[0], a[1])
    elif c["kind"] == "rayleigh_fail":
        out = _pack_rayleigh_fail(c["kwargs"])
    else:
        out = _pack_rayleigh_ok(c["kwargs"])
        if wins and getattr(wins[-1]["mw"], "ok", False):
            try:
                out["why"] = _why_window(wins[-1])
            except Exception as exc:                                            # noqa: BLE001
                print(f"    (why-window unavailable: {type(exc).__name__}: {exc})", flush=True)
    out.update(meta)
    return out


def build_payload(key: str, day_list: list, methods: list) -> tuple[dict, float]:
    """Calibrate every (day, method) in-process and pack the panel payload.

    Separate from main() so scripts/mockup_station_page.py can build the same payload for the panel
    it embeds. Expects the ALC_* environment to be set by the caller.
    """
    _install_capture()
    import importlib
    rnc = importlib.import_module("scripts.run_network_calibration")
    census = json.loads((REPO / "validation/scope_l1_2026_census.json").read_text(encoding="utf-8"))
    rows = census if isinstance(census, list) else census.get("streams", [])
    stream = next(r for r in rows if f"{r['wmo']}_{r['ident']}" == key)

    days, dates = {}, []
    t_start = time.perf_counter()
    for d in day_list:
        ds = d.strftime("%Y%m%d")
        entry = {}
        for m in methods:
            t0 = time.perf_counter()
            p = _run_one(rnc, stream, d, m)
            entry[m] = p
            print(f"  {ds} {m:<8} {p['kind']:<14} flag={p.get('flag')}  "
                  f"{time.perf_counter() - t0:5.1f}s", flush=True)
        if any(v.get("kind") != "error" for v in entry.values()):
            days[ds] = entry
            dates.append(ds)
    payload = {"station": {"key": key, "wmo": stream["wmo"], "ident": stream["ident"],
                           "type": stream.get("type", ""),
                           "name": stream.get("name") or stream.get("station_name") or ""},
               "methods": methods, "dates": dates, "days": days}
    return payload, t_start


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--key", required=True)
    ap.add_argument("--start")
    ap.add_argument("--end")
    # Rayleigh and cloud are close to mutually exclusive by construction (one needs a clear night,
    # the other a liquid cloud), so no short contiguous window shows every outcome class. An explicit
    # day list lets the mockup cover them all without carrying a fortnight of payload.
    ap.add_argument("--dates", default="", help="comma-separated YYYYMMDD, instead of --start/--end")
    ap.add_argument("--methods", default="rayleigh,cloud")
    ap.add_argument("--l1-root", required=True)
    ap.add_argument("--cams", default="")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=Path("./_daily_panel_work"))
    args = ap.parse_args()

    os.environ["ALC_L1_ROOT"] = args.l1_root
    os.environ["ALC_FULLCAL_DIR"] = str(args.work.resolve())
    os.environ["PLOTS"] = "1"
    if args.cams:
        os.environ["ALC_CAMS_DIR"] = args.cams
    args.work.mkdir(parents=True, exist_ok=True)

    if args.dates:
        day_list = [datetime.strptime(s.strip(), "%Y%m%d")
                    for s in args.dates.split(",") if s.strip()]
    elif args.start and args.end:
        d, day_list = datetime.strptime(args.start, "%Y%m%d"), []
        d1 = datetime.strptime(args.end, "%Y%m%d")
        while d <= d1:
            day_list.append(d)
            d += timedelta(days=1)
    else:
        ap.error("give either --dates or both --start and --end")
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    payload, t_start = build_payload(args.key, day_list, methods)
    if not payload["dates"]:
        print("nothing to show", file=sys.stderr)
        sys.exit(2)
    dates = payload["dates"]
    raw = json.dumps(payload, separators=(",", ":"))
    gz_kb = len(gzip.compress(raw.encode())) / 1024

    import plotly.offline as pyo
    html = (HTML.replace("__META__", json.dumps(
        {"n_days": len(dates), "gz_kb": round(gz_kb, 1),
         # The unit production would fetch is ONE day and ONE method -- that is what to compare
         # against a PNG, not the whole embedded bundle this self-contained mockup carries.
         "per_unit_kb": round(gz_kb / max(len(dates) * len(methods), 1), 1),
         "curated": bool(args.dates),
         "elapsed_s": round(time.perf_counter() - t_start, 1)}))
        .replace("__PAYLOAD__", raw)
        .replace("__PLOTLY__", pyo.get_plotlyjs()))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html, encoding="utf-8")
    print(f"\n{len(dates)} days x {len(methods)} methods   payload {len(raw)/1024:.0f} KB "
          f"({gz_kb:.0f} KB gzip, {gz_kb/max(len(dates),1):.0f} KB per day)")
    print(f"mockup: {args.out}")


# The panel is published as three fragments -- CSS, body, JS -- rather than one page string, so
# scripts/mockup_station_page.py can embed the SAME panel inside the full station page instead of
# reimplementing it. Two copies of this much Plotly wiring would diverge within a week.
PANEL_CSS = r"""
 :root { --line:#dbe3ea; --ink:#1a2530; --dim:#66707a; --ok:#2ca02c; --bad:#b00020; }
 * { box-sizing:border-box; }
 body { font:14px/1.5 -apple-system,"Segoe UI",Roboto,Arial,sans-serif; margin:0; color:var(--ink);
        background:#f7f9fb; }
 .hdr { background:#fff3cd; border-bottom:1px solid #ffe08a; padding:8px 18px; font-size:13px; }
 .top { display:flex; align-items:baseline; gap:14px; flex-wrap:wrap; padding:12px 18px 0; }
 .top h1 { font-size:19px; margin:0; }
 .sub { color:var(--dim); font-size:13px; }
 .layout { display:grid; grid-template-columns:214px 1fr; gap:16px; padding:12px 18px 26px;
           align-items:start; }
 .rail { position:sticky; top:10px; }
 .card { background:#fff; border:1px solid var(--line); border-radius:10px; padding:9px; }
 .card + .card { margin-top:12px; }
 .toolbar { display:flex; align-items:center; gap:10px; flex-wrap:wrap; margin-bottom:12px; }
 button { font:inherit; font-size:13px; padding:4px 11px; border:1px solid #c3ceda;
          background:#fff; border-radius:7px; cursor:pointer; color:var(--ink); }
 button:hover:not(:disabled) { background:#eef4fb; border-color:#2a5a82; }
 button:disabled { opacity:.4; cursor:default; }
 .seg { display:inline-flex; border:1px solid #c3ceda; border-radius:7px; overflow:hidden; }
 .seg button { border:0; border-radius:0; }
 .seg button + button { border-left:1px solid #c3ceda; }
 .seg button.on { background:#2a5a82; color:#fff; }
 .datelbl { font-weight:600; font-size:15px; min-width:118px; text-align:center; }
 .chip { display:inline-block; padding:2px 9px; border-radius:11px; font-size:12px;
         font-weight:600; border:1px solid transparent; }
 .chip.ok { background:#e8f6ea; color:#1a6b28; border-color:#b6e0bd; }
 .chip.bad { background:#fdecef; color:var(--bad); border-color:#f5c2c7; }
 .chip.na { background:#eef1f4; color:var(--dim); border-color:#dbe3ea; }
 .toprow { display:grid; grid-template-columns:2fr 1fr; gap:12px; align-items:start; }
 .diagrow { display:grid; grid-template-columns:repeat(3,1fr); gap:12px; margin-top:12px; }
 h2 { font-size:13px; margin:2px 0 8px; text-transform:uppercase; letter-spacing:.04em;
      color:var(--dim); }
 .cal { display:grid; grid-template-columns:repeat(7,1fr); gap:3px; }
 .cal .dow { font-size:10px; color:var(--dim); text-align:center; }
 .cal .d { aspect-ratio:1; border-radius:5px; border:1px solid var(--line); font-size:11px;
           display:flex; align-items:center; justify-content:center; cursor:default;
           background:#f2f4f7; color:#b6bec6; }
 .cal .d.has { cursor:pointer; color:#fff; font-weight:600; border-color:transparent; }
 .cal .d.has:hover { outline:2px solid #2a5a82; }
 .cal .d.cur { outline:2.5px solid #1a2530; }
 .mon { font-weight:600; font-size:12px; margin:2px 0 6px; }
 .legend { font-size:11px; color:var(--dim); margin-top:9px; line-height:1.7; }
 .sw { display:inline-block; width:10px; height:10px; border-radius:3px; margin-right:5px;
       vertical-align:-1px; }
 .msg { background:#fff; border:1px solid var(--line); border-radius:10px; padding:11px 13px;
        margin-top:12px; }
 .msg.bad { border-left:5px solid var(--bad); background:#fdecef; }
 .msg.ok { border-left:5px solid var(--ok); }
 .msg b { font-size:14px; }
 .msg .det { color:var(--dim); font-size:12.5px; margin-top:4px; }
 pre.summary { font:12px/1.45 ui-monospace,Consolas,monospace; margin:0; white-space:pre-wrap; }
 /* min-height, not padding: both cards' control rows must occupy exactly the same vertical space,
    whether they hold buttons or a line of text, or the two plots below them stop lining up. */
 .ctrls { display:flex; gap:8px; align-items:center; flex-wrap:wrap; margin-bottom:6px;
          min-height:30px; }
 .ctrls label { font-size:11.5px; color:var(--dim); }
 .note { font-size:11.5px; color:var(--dim); }
 .notebar { margin:4px 2px 2px; border-top:1px solid #eef2f6; padding-top:6px; }
 .hd { font-size:12px; text-transform:uppercase; letter-spacing:.04em; color:var(--dim);
       margin-right:4px; }
 .axlbl { display:inline-flex; align-items:center; gap:5px; font-size:11px; color:var(--dim); }
 .empty { text-align:center; color:var(--dim); padding:60px 10px; }
 .whycard { margin-top:12px; }
 .wtop { margin:0 0 8px; font-size:13px; }
 .formula { font:12.5px/1.5 ui-monospace,Consolas,monospace; background:#f4f7fa;
            border:1px solid #e3e9ef; border-radius:7px; padding:8px 11px; margin:0 0 9px;
            overflow-x:auto; white-space:pre; }
 .wtab { border-collapse:collapse; width:100%; font-size:12.5px; }
 .wtab th { text-align:left; font-size:11px; text-transform:uppercase; letter-spacing:.03em;
            color:var(--dim); border-bottom:1px solid var(--line); padding:3px 8px; }
 .wtab td { padding:3px 8px; border-bottom:1px solid #eef2f6; }
 .wtab td:first-child { font-family:ui-monospace,Consolas,monospace; color:#2a5a82; }
 .wtab .num { text-align:right; font-variant-numeric:tabular-nums; }
 .wtab .wy { color:var(--dim); font-size:11.5px; }
 .wtab tr.drive { background:#fff8e6; }
 .wtab tr.tot td { border-bottom:0; border-top:1px solid var(--line); }
 .gates { margin-top:9px; display:flex; flex-wrap:wrap; gap:6px; }
 .gate { font-size:11px; background:#eef4fb; border:1px solid #d6e4f2; border-radius:11px;
         padding:2px 9px; color:#0b3d61; }
 .gate .gv { color:#2a5a82; font-variant-numeric:tabular-nums; }
 @media (max-width:1100px) { .toprow { grid-template-columns:1fr; }
   .layout { grid-template-columns:1fr; } .rail { position:static; }
   .diagrow { grid-template-columns:1fr 1fr; } }
"""

#: Chrome for the STANDALONE mockup only. The station page brings its own header, which is why the
#: panel JS treats #meta / #stitle / #ssub as optional.
PANEL_HDR = r"""<div class="hdr"><b>MOCKUP — unified daily calibration panel.</b> Every plot is drawn
in your browser from data captured from <b>real</b> calibrations; no PNG exists.
<span id="meta"></span></div>

<div class="top"><h1 id="stitle"></h1><span class="sub" id="ssub"></span></div>"""

PANEL_BODY = r"""<div class="layout">
 <div class="rail">
  <div class="card">
   <h2>Calendar</h2>
   <div id="cal"></div>
   <div class="legend" id="legend"></div>
  </div>
  <div class="card">
   <h2>Flags this day</h2>
   <div id="flags"></div>
  </div>
 </div>

 <div>
  <div class="toolbar">
   <button id="prev">← previous</button>
   <span class="datelbl" id="datelbl"></span>
   <button id="next">next →</button>
   <span class="seg" id="mswitch"></span>
   <span id="cst" class="sub"></span>
  </div>

  <!-- ONE card, ONE Plotly figure, two subplots on a SHARED y axis. Two separate figures could
       only be aligned by hand-matching margins, control rows and card heights, and any later edit
       silently broke it; here the alignment is structural and the y zoom is shared for free. -->
  <div class="card">
   <div class="ctrls">
    <b class="hd">Time–height &amp; profile</b>
    <span class="seg" id="vswitch"></span>
    <!-- Scoped to the profile ON PURPOSE: the curtain is quantised in log space and bands on a
         linear ramp, so it stays log10. Saying so on the control is cheaper than explaining it. -->
    <label class="axlbl">profile x
      <span class="seg" id="xswitch">
        <button data-x="lin">linear</button><button data-x="log">log</button></span></label>
   </div>
   <div id="panel"></div>
   <!-- Provenance goes BELOW the figure: it describes what was drawn, so above the plot it just
        pushed the panel down and was read before there was anything to read it about. -->
   <div class="notebar"><span class="note" id="curtitle"></span></div>
  </div>

  <div id="msg"></div>
  <div id="why"></div>
  <div class="diagrow" id="diag"></div>
 </div>
</div>

"""

PANEL_JS = r"""
const D = JSON.parse(document.getElementById('payload').textContent);
const META = __META__;
// #meta / #stitle / #ssub belong to the standalone mockup's chrome; the station page supplies its
// own header and does not define them, so every write here is optional by construction.
const setText = (id, t) => { const e = document.getElementById(id); if (e) e.textContent = t; };
setText('meta',
  ` ${META.n_days} days × ${D.methods.length} methods, all embedded = ${META.gz_kb} KB gzipped. ` +
  `The unit production would fetch is one day and one method: ~${META.per_unit_kb} KB, against ` +
  `~2300 KB for the PNG it replaces.` + (META.curated
    ? ' Days were picked to cover every outcome (both methods, each alone, and rejected):' +
      ' Rayleigh needs a clear night and the cloud method a liquid cloud, so no short run of' +
      ' consecutive days contains them all.' : ''));
const S = D.station;
setText('stitle', (S.name || S.wmo) + ' · ' + S.type);
setText('ssub', `${S.wmo} · unit ${S.ident} · ${D.dates.length} days loaded`);

const dec = (b64, type) => {
  const s = atob(b64), u = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i);
  return type === 'u8' ? u : new Float32Array(u.buffer);
};
const f32 = b => Array.from(dec(b, 'f32'));
const LAY = { margin:{l:62,r:10,t:16,b:44}, template:'plotly_white', font:{size:11} };
const CFG = { responsive:true, displaylogo:false };
// Matching axis VALUES is not enough for the curtain and the profile to line up on screen: the
// plotting box is what is left after the margins, so different top/bottom margins put the same
// height at different pixels. These two numbers are shared by both panels for that reason -- change
// them together or the alignment silently breaks.
const PANEL_H = 470, PANEL_T = 62, PANEL_B = 70;

// uint8 curtain -> 2-D physical log10 values (255 = missing -> null, so gaps stay gaps)
function grid(c) {
  const q = dec(c.b64, 'u8'), [ny, nx] = c.shape, zz = [];
  for (let j = 0; j < ny; j++) {
    const row = new Array(nx);
    for (let i = 0; i < nx; i++) { const v = q[j*nx+i];
      row[i] = v === 255 ? null : c.lo + (v/254)*(c.hi - c.lo); }
    zz.push(row);
  }
  return zz;
}
const xOf = c => c.time || c.hours;
const isDate = c => !!c.time;

// ---------------------------------------------------------------- state
let curDate = D.dates[0];
let curMethod = null;
let curView = 'beta';        // resolved on every render by buildViewSwitch
let userView = null;         // set ONLY by an explicit click on the selector
// null = "follow the view's own default" (RCS is linear, signal is log, exactly as the PNG draws
// them). An explicit click on the linear/log switch pins it until the view changes again.
let curLogX = null;

const VIEW_NAMES = { rcs:'Raw RCS', sig:'Calibrated signal', beta:'Attenuated backscatter' };
const VIEW_ORDER = ['rcs', 'sig', 'beta'];

// One reading rule for every window-search panel: GREEN IS THE GOOD END. R² rises towards 1 (good)
// so it runs white->green; |intercept| is best at 0 so it runs green->white. The slope is the odd
// one out -- it is good at ZERO and bad in BOTH directions -- so it gets a diverging scale pinned
// symmetrically about zero rather than a sequential one. Without that, a slope grid coloured
// sequentially reads as "more is better", which is the opposite of the truth.
const GREEN_UP   = [[0,'#ffffff'],[0.5,'#a8ddb5'],[1,'#0b6b2e']];   // 0 = bad  -> 1 = good
const GREEN_DOWN = [[0,'#0b6b2e'],[0.5,'#a8ddb5'],[1,'#ffffff']];   // 0 = good -> high = bad
const GRID_SPEC = {
  slopes:     { title:'Slope',       diverging:true },
  intercepts: { title:'|Intercept|', scale:GREEN_DOWN, zmin:0, exp:true },
  r_squared:  { title:'R²',          scale:GREEN_UP,   zmin:0, zmax:1 },
};

const dayOf = (d, m) => (D.days[d] || {})[m] || null;
// A six-month page cannot embed every day's payload (~340 KB each), so the station dashboard ships
// a tiny per-day INDEX (flag, constant, has_fig) and fetches the full payload only for the day
// being looked at. Everything that just needs to know HOW a day turned out -- the calendar, the
// flag chips, which methods ran -- reads the index, so it is complete from the first paint;
// only the panel itself waits for a fetch.
const summaryOf = (d, m) => dayOf(d, m) || ((D.index || {})[d] || {})[m] || null;
// Calibrated iff a constant was produced. Not "which figure was drawn" -- those disagree: flag 0.5
// draws the success figure AND yields a constant, a cloud rejection draws the full cloud figure and
// yields none. See the sentinel note in _run_one.
const isOK = p => !!(p && p.constant != null);
const hasFig = p => !!(p && p.curtain && p.curtain.b64);
// THREE greens: which method calibrated is the thing an operator scans the calendar for, and a
// single green hid it. Deepest = both methods; the two single-method greens are separated by
// lightness, not hue, so the "both is more" reading survives. The pale one takes dark text --
// white on it is not legible at 11 px.
// The METHOD colours, not a green ramp: blue is Rayleigh and green is liquid-cloud everywhere else
// on the dashboard (monitoring/config.py METHOD_COLORS), so the calendar has to agree or the same
// colour would mean two things on one page. Cyan is the additive mix, for the nights both produced.
const CAL_COL = {
  both:     { bg:'#17a2b8', fg:'#ffffff', lbl:'calibrated — Rayleigh <b>and</b> cloud' },
  rayleigh: { bg:'#1f77b4', fg:'#ffffff', lbl:'calibrated — Rayleigh only' },
  cloud:    { bg:'#2ca02c', fg:'#ffffff', lbl:'calibrated — cloud only' },
  fig:      { bg:'#d98c00', fg:'#ffffff', lbl:'rejected — diagnostics shown' },
  none:     { bg:'#b00020', fg:'#ffffff', lbl:'rejected — no figure drawn' },
};
function dayColour(d) {
  const present = D.methods.map(m => summaryOf(d, m)).filter(Boolean);
  if (!present.length) return null;
  const okm = D.methods.filter(m => isOK(summaryOf(d, m)));
  if (!okm.length) return present.some(hasFig) ? CAL_COL.fig : CAL_COL.none;
  // "both" means every method CONFIGURED for this station, so a Rayleigh-only instrument (a CHM15k
  // never gets a cloud calibration -- it saturates in liquid cloud) still reads as a full success
  // rather than as a permanent half-result.
  if (okm.length >= D.methods.length) return CAL_COL.both;
  return CAL_COL[okm[0]] || CAL_COL.both;
}
function buildLegend() {
  const rows = D.methods.length > 1
    ? [CAL_COL.both, CAL_COL.rayleigh, CAL_COL.cloud]
    : [{ bg:CAL_COL.both.bg, lbl:'calibrated' }];
  rows.push(CAL_COL.fig, CAL_COL.none, { bg:'#f2f4f7', lbl:'no data' });
  document.getElementById('legend').innerHTML = rows.map(r =>
    `<span class="sw" style="background:${r.bg};border:1px solid #dbe3ea"></span>${r.lbl}`
  ).join('<br>');
}

// ---------------------------------------------------------------- calendar
function buildCal() {
  const host = document.getElementById('cal');
  host.innerHTML = '';
  const months = [...new Set(D.dates.map(d => d.slice(0,6)))];
  months.forEach(mo => {
    const y = +mo.slice(0,4), m = +mo.slice(4,6);
    const lbl = document.createElement('div');
    lbl.className = 'mon';
    lbl.textContent = new Date(y, m-1, 1).toLocaleString('en', {month:'long', year:'numeric'});
    host.appendChild(lbl);
    const g = document.createElement('div');
    g.className = 'cal';
    ['M','T','W','T','F','S','S'].forEach(t => {
      const e = document.createElement('div'); e.className = 'dow'; e.textContent = t; g.appendChild(e); });
    const first = new Date(y, m-1, 1), lead = (first.getDay() + 6) % 7;   // Monday-first
    for (let i = 0; i < lead; i++) g.appendChild(document.createElement('div'));
    const ndays = new Date(y, m, 0).getDate();
    for (let dd = 1; dd <= ndays; dd++) {
      const ds = `${y}${String(m).padStart(2,'0')}${String(dd).padStart(2,'0')}`;
      const cell = document.createElement('div');
      cell.className = 'd'; cell.textContent = dd;
      const col = dayColour(ds);
      if (col) {
        cell.classList.add('has');
        cell.style.background = col.bg;
        cell.style.color = col.fg;
        cell.title = ds + ' — ' + D.methods.map(mm => {
          const p = dayOf(ds, mm);
          return p ? `${mm}: flag ${p.flag} (${isOK(p) ? 'calibrated' : 'rejected'})`
                   : `${mm}: not run`;
        }).join(', ');
        cell.addEventListener('click', () => { curDate = ds; render(); });
      }
      cell.dataset.ds = ds;
      g.appendChild(cell);
    }
    host.appendChild(g);
  });
}
function markCal() {
  document.querySelectorAll('.cal .d').forEach(e =>
    e.classList.toggle('cur', e.dataset.ds === curDate));
}

// ---------------------------------------------------------------- switches
function buildMethodSwitch() {
  const host = document.getElementById('mswitch');
  host.innerHTML = '';
  const avail = D.methods.filter(m => summaryOf(curDate, m));
  if (!avail.includes(curMethod)) curMethod = avail[0] || D.methods[0];
  // The switch only appears when the day genuinely has two calibrations -- a lone disabled
  // button would imply a second product exists and is broken, which is not what happened.
  if (avail.length < 2) { host.style.display = 'none'; return; }
  host.style.display = '';
  avail.forEach(m => {
    const b = document.createElement('button');
    b.textContent = m === 'cloud' ? 'Cloud (O\u2019Connor)' : 'Rayleigh';
    b.className = m === curMethod ? 'on' : '';
    b.addEventListener('click', () => { curMethod = m; render(); });
    host.appendChild(b);
  });
}
function buildViewSwitch(views) {
  const host = document.getElementById('vswitch');
  host.innerHTML = '';
  const have = VIEW_ORDER.filter(v => views && views[v]);
  // Attenuated backscatter is THE default whenever the night can form it: it is the physical
  // product and the one view both methods share. It has to be re-preferred rather than remembered,
  // because a rejected night has no constant and so cannot form it -- without this the selector
  // would stay stuck on the fallback view after the operator moved back to a calibrated night.
  curView = (userView && have.includes(userView)) ? userView
          : (have.includes('beta') ? 'beta' : (have[have.length - 1] || 'beta'));
  have.forEach(v => {
    const b = document.createElement('button');
    b.textContent = views[v].name || VIEW_NAMES[v];
    b.className = v === curView ? 'on' : '';
    b.addEventListener('click', () => { userView = v; curLogX = null; render(); });
    host.appendChild(b);
  });
  const eff = curLogX === null ? !!(views[curView] && views[curView].logx) : curLogX;
  document.querySelectorAll('#xswitch button').forEach(b =>
    b.className = ((b.dataset.x === 'log') === eff) ? 'on' : '');
}
document.querySelectorAll('#xswitch button').forEach(b =>
  b.addEventListener('click', () => { curLogX = b.dataset.x === 'log'; render(); }));

// ------------------------------------------------- curtain + profile, ONE figure, SHARED y axis
// The view selector drives BOTH panels. That is possible without shipping three curtains because
// the three quantities are exact algebra on the stored log10 field:
//     beta_att = RCS / C_L          -> a constant offset in log space
//     signal   = beta_att / range^2 -> that offset, plus a per-GATE offset
// so the curtain is transformed in place. Quantisation happened in the base field's log space and
// both transforms are affine there, so nothing is lost relative to the stored curtain.
// The curtain is ALWAYS shown on log10. That is not a preference: it is quantised to 254 levels
// spaced uniformly in LOG space, so a linear ramp puts ~37 of those levels across the top decade
// and crushes everything below into the first few colours -- it bands visibly. The profile has no
// such constraint (it ships as float32), which is why the linear/log switch applies to it alone.
// A genuinely linear curtain would need a second, linearly-quantised payload.
function curtainZ(c, view) {
  const zz = grid(c), rk = f32(c.range_km);
  const lc = (c.log_c === null || c.log_c === undefined) ? null : +c.log_c;
  let off = 0, byGate = false;
  if (view === 'rcs')  off = (c.base === 'rcs') ? 0 : (lc === null ? 0 : lc);
  if (view === 'beta') off = (c.base === 'rcs') ? (lc === null ? 0 : -lc) : 0;
  if (view === 'sig') { off = (c.base === 'rcs') ? (lc === null ? 0 : -lc) : 0; byGate = true; }
  let lo = c.lo + off, hi = c.hi + off;    // keep the PNG's own 5-95 percentile limits
  const sample = [];
  for (let j = 0; j < zz.length; j++) {
    // rows are range gates; clamp so the first gate cannot give log10(0) = -Infinity
    const extra = byGate ? -2 * Math.log10(Math.max(rk[j] * 1000, 1)) : 0;
    const row = zz[j];
    for (let i = 0; i < row.length; i++) {
      if (row[i] === null) continue;
      row[i] += off + extra;
      if (byGate && ((i + j) % 7 === 0)) sample.push(row[i]);
    }
  }
  if (byGate && sample.length > 20) {
    // The per-gate shift is not a constant, so the stored percentiles no longer apply. Take fresh
    // 5/95 percentiles off a subsample rather than min/max -- a single hot gate would otherwise
    // stretch the scale and posterise the whole panel, which is the same failure in another guise.
    sample.sort((a, b) => a - b);
    lo = sample[Math.floor(0.05 * (sample.length - 1))];
    hi = sample[Math.floor(0.95 * (sample.length - 1))];
  }
  if (!(hi > lo)) { lo = c.lo + off; hi = c.hi + off; }
  return { z: zz, lo: lo, hi: hi };
}

const Z_LABEL = { rcs:'RCS', sig:'signal', beta:'β_att' };

function drawPanel(p) {
  const host = document.getElementById('panel');
  const c = p && p.curtain, pr = p && p.profile;
  buildViewSwitch(pr && pr.views);
  if (!c || !c.b64) {
    host.innerHTML = '<p class="empty">no time–height data for this night</p>'; return; }
  const views = (pr && pr.views) || {};
  const v = views[curView];
  const useLog = curLogX === null ? !!(v && v.logx) : curLogX;
  // A view the PROFILE cannot form (no constant on a rejected night) must not be formed for the
  // curtain either, or the two halves of one card would be showing different quantities.
  const cv = v ? curView : (c.base === 'rcs' ? 'rcs' : 'beta');
  const X = xOf(c), rk = f32(c.range_km);
  const Z = curtainZ(c, cv);

  const tr = [{ type:'heatmap', z:Z.z, x:X, y:rk, colorscale:'Viridis', zmin:Z.lo, zmax:Z.hi,
      xaxis:'x', yaxis:'y',
      colorbar:{ title:{ text:'log₁₀(' + Z_LABEL[cv] + ')' }, thickness:12, x:0.615, len:0.92 },
      hovertemplate:'%{x}<br>%{y:.2f} km<br>log₁₀ %{z:.2f}<extra></extra>' }];
  const shapes = [];
  (c.bands || []).forEach(b => {
    const runs = b.runs || [];
    if (!runs.length) return;
    // Drawn as a TRACE, not a layout shape. Shapes cannot be legend items, so these masks used to
    // be a shape plus a dummy null-point trace carrying the name -- and clicking that legend entry
    // toggled the dummy, leaving the mask on screen. One filled trace per mask makes the legend
    // click do what it looks like it does. Null separators split it into one polygon per run.
    const bx = [], by = [], y0 = 0, y1 = c.y_max_km;
    runs.forEach(r => {
      const a = X[r[0]], z = X[Math.min(r[1] + 1, X.length - 1)];
      bx.push(a, z, z, a, a, null);
      by.push(y0, y0, y1, y1, y0, null);
    });
    tr.push({ x:bx, y:by, mode:'lines', fill:'toself', fillcolor:b.color, line:{width:0},
        xaxis:'x', yaxis:'y', name:b.name, hoverinfo:'skip' });
  });
  if (c.high_cloud_ylo && c.has_high_cloud) {
    // Two gap-free traces: ceiling at y_hi, floor at y_lo (== y_hi outside the mask, so the band
    // collapses to zero height there instead of being bridged across a null). 'hvh' = step mid.
    tr.push({ x:X, y:X.map(() => c.y_hi_km), mode:'lines', line:{width:0, shape:'hvh'},
        xaxis:'x', yaxis:'y', hoverinfo:'skip', showlegend:false });
    tr.push({ x:X, y:c.high_cloud_ylo, mode:'lines', fill:'tonexty', xaxis:'x', yaxis:'y',
        fillcolor:'rgba(25,211,243,0.22)', line:{width:0, shape:'hvh'},
        name:'high cloud (masked above fit)', hoverinfo:'skip' });
  }
  (c.cbh || []).forEach(s => tr.push({ x:X, y:s.y, mode:'markers', name:s.name,
      xaxis:'x', yaxis:'y',
      marker:{ size:s.size, color:s.color, line:{color:'#000', width:0.7} },
      hovertemplate:'%{x}<br>cloud base %{y:.2f} km<extra></extra>' }));
  if (c.layer_km) {
    // Drawn across BOTH subplots: the calibration layer is the feature the two panels are meant to
    // be compared through, so it must be one continuous band, not two that nearly line up.
    ['x domain','x2 domain'].forEach((xr, i) => shapes.push({ type:'rect', xref:xr,
        yref: i ? 'y2' : 'y', x0:0, x1:1, y0:c.layer_km[0], y1:c.layer_km[1],
        fillcolor:'rgba(255,193,7,0.18)', line:{color:'#ffc107', width:1.4}, layer:'above' }));
    tr.push({ x:[null], y:[null], mode:'lines', xaxis:'x', yaxis:'y',
        line:{color:'#ffc107', width:2}, name:c.layer_name || 'calibration layer' });
  }

  // ---- profile, on the second subplot
  let xa2 = { title:{ text:'(no profile)' }, domain:[0.72, 1.0], anchor:'y2' };
  if (v) {
    const y = f32(pr.y);
    tr.push({ x:f32(v.x), y:y, mode:'lines', name:pr.mean_of || 'observed', xaxis:'x2', yaxis:'y2',
        line:{color:'#1f77b4', width:1.1}, connectgaps:false });
    if (v.mol) tr.push({ x:f32(v.mol), y:y, mode:'lines', xaxis:'x2', yaxis:'y2',
        name:v.mol_name || 'molecular', line:{color:'#d62728', width:1.1, dash:'dash'} });
    if (pr.med_cbh_km !== undefined) {
      shapes.push({ type:'line', xref:'x2 domain', yref:'y2', x0:0, x1:1,
          y0:pr.med_cbh_km, y1:pr.med_cbh_km, line:{color:'#d62728', width:1.1, dash:'dot'} });
      tr.push({ x:[null], y:[null], mode:'lines', xaxis:'x2', yaxis:'y2',
          line:{color:'#d62728', width:1.4, dash:'dot'},
          name:pr.med_cbh_label || 'median cloud base' });
    }
    if (pr.cal_band && !c.layer_km) shapes.push({ type:'rect', xref:'x2 domain', yref:'y2',
        x0:0, x1:1, y0:pr.cal_band[0], y1:pr.cal_band[1], fillcolor:'rgba(255,193,7,0.20)',
        line:{color:'#ffc107', width:1.2}, layer:'below' });
    xa2 = { title:{ text:v.label }, domain:[0.72, 1.0], anchor:'y2' };
    if (useLog) { xa2.type = 'log'; xa2.exponentformat = 'e'; }
  }

  Plotly.newPlot(host, tr, Object.assign({}, LAY, {
    height:PANEL_H, shapes:shapes, margin:{l:66, r:10, t:PANEL_T, b:PANEL_B},
    xaxis:{ title:{ text: isDate(c) ? 'Time (UTC)' : 'Hours since start' }, domain:[0, 0.60],
            anchor:'y', automargin:true, type: isDate(c) ? 'date' : 'linear' },
    yaxis:{ title:{ text:c.y_label }, range:[0, c.y_max_km], anchor:'x' },
    xaxis2: xa2,
    // matches:'y' is what makes the alignment structural -- one y zoom moves both panels, and no
    // relayout handler is needed to keep them in step.
    yaxis2:{ matches:'y', anchor:'x2', showticklabels:false },
    legend:{ orientation:'h', y:1.10, x:0, xanchor:'left', font:{size:9.5} } }), CFG);
  const el = document.getElementById('curtitle');
  if (el) el.textContent =
    c.shape[1] + ' × ' + c.shape[0] + ' block-averaged from the native grid (' +
    c.st_t + '×' + c.st_r + ' cells per block) · curtain colour is log₁₀ (it is quantised in log ' +
    'space) · both panels share one vertical axis' +
    (pr && pr.mean_of ? ' · ' + pr.mean_of : '') +
    ((c.bands || []).some(b => (b.runs || []).length)
      ? ' · click a mask in the legend to hide it' : '');
}

// ------------------------------------------------------- why THIS window (Rayleigh only)
// The window grids show WHERE the search looked; they do not say why one cell won. This card
// re-states the objective the selector maximises, fills in the winning cell's numbers, and names
// the term that actually cost the most -- which is the question an operator asks when a night
// fits somewhere unexpected.
function drawWhy(p) {
  const host = document.getElementById('why');
  if (!host) return;
  const w = (p && p.why) || (p && p.diag && p.diag.why);
  if (!w) { host.innerHTML = ''; return; }
  if (w.kind === 'cloud') { drawWhyCloud(host, w); return; }
  const sgn = v => (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(3);
  const rows = w.terms.map(t => {
    const drive = t.name === w.driver;
    return `<tr class="${drive ? 'drive' : ''}">
      <td>${t.sym}</td><td>${t.name}${drive ? ' <b>← the deciding penalty</b>' : ''}</td>
      <td class="num">${t.value}</td><td class="num">${sgn(t.contrib)}</td>
      <td class="wy">${t.why}</td></tr>`;
  }).join('');
  const gates = w.gates.map(g =>
    `<span class="gate"><b>${g.name}</b> ${g.limit} <span class="gv">→ ${g.value}</span></span>`
  ).join('');
  host.innerHTML = `<div class="card whycard">
    <h2>Why this window?</h2>
    <p class="wtop">The fit ran over <b>${w.window_km[0].toFixed(2)}–${w.window_km[1].toFixed(2)} km</b>
      (centre ${w.centre_km.toFixed(2)} km, half-length ${w.half_km.toFixed(2)} km).
      It was chosen as the <b>${w.message || 'best-scoring window'}</b> among
      <b>${w.n_eligible}</b> of ${w.n_total} candidate windows that passed every gate.</p>
    <pre class="formula">Q = R² − ${w.weights.w_ratio}·|SR−1| − ${w.weights.w_resid}·resid − ${w.weights.w_snr}·σ_ratio − ${w.weights.w_tvar}·CV_t − ${w.weights.w_rel}·relerr + ${w.weights.w_npts}·n̂</pre>
    <table class="wtab"><tr><th>term</th><th>what it measures</th><th class="num">value</th>
      <th class="num">contributes</th><th>why it is in the score</th></tr>${rows}
      <tr class="tot"><td></td><td><b>Q for this window</b></td><td></td>
        <td class="num"><b>${w.total.toFixed(3)}</b></td><td></td></tr></table>
    <div class="gates">${gates}</div></div>`;
}

// The cloud equivalent: the O'Connor chain stated as a chain, and the three filters shown as the
// FUNNEL they are. The old panel printed each stage's counter dict, which answers "what did each
// filter count" when the operator's question is "which stage cost me the night".
function drawWhyCloud(host, w) {
  const n = v => (v === null || v === undefined || !isFinite(v)) ? '—' : (+v).toPrecision(4);
  const rows = w.stages.map(s => {
    const drive = s.name === w.driver;
    const det = Object.keys(s.detail || {}).length
      ? Object.entries(s.detail).map(([k, v]) => `${k}: ${v}`).join(' · ') : '—';
    return `<tr class="${drive ? 'drive' : ''}">
      <td>${s.name}${drive ? ' <b>← removed the most</b>' : ''}</td>
      <td class="wy">${s.what}</td>
      <td class="num">−${s.removed}</td><td class="num">${s.after}</td>
      <td class="wy">${det}</td></tr>`;
  }).join('');
  const verdict = w.n_used
    ? `<b>${w.n_used}</b> of ${w.n_total} profiles survived all three filters and set the constant.`
    : `<b>No profile survived</b> all three filters, so the night produced no constant — the ` +
      `absence of any green in the panel above IS the rejection.`;
  host.innerHTML = `<div class="card whycard">
    <h2>Why this cloud calibration?</h2>
    <p class="wtop">A fully-attenuating liquid cloud returns a KNOWN integrated backscatter, so the
      integral over the fixed ${w.gate_km[0].toFixed(2)}–${w.gate_km[1].toFixed(2)} km gate measures
      the calibration rather than the cloud. ${verdict}</p>
    <pre class="formula">β corrected for ${w.wv ? 'water vapour (two-way) and ' : ''}multiple scattering η(range)
S_apparent = 1 / (2 · ∫ β dz)          over the ${w.gate_km[0].toFixed(2)}–${w.gate_km[1].toFixed(2)} km gate
C          = S_consistent / ${w.s_theo.toFixed(2)} sr      (theoretical S for liquid water)
C_L        = C_applied / C             (Wiegner, comparable with the Rayleigh constant)</pre>
    <table class="wtab"><tr><th>filter stage</th><th>what it tests</th><th class="num">removed</th>
      <th class="num">left</th><th>breakdown</th></tr>${rows}</table>
    <div class="gates">
      <span class="gate"><b>median S</b> <span class="gv">${n(w.s_med)} sr</span></span>
      <span class="gate"><b>theoretical S</b> <span class="gv">${w.s_theo.toFixed(2)} sr</span></span>
      <span class="gate"><b>C = S/S_theo</b> <span class="gv">${n(w.C)}${
        w.C_std !== null ? ' ± ' + n(w.C_std) : ''}</span></span>
      <span class="gate"><b>C applied</b> <span class="gv">${n(w.applied)}</span></span>
      <span class="gate"><b>C_L = C_applied/C</b> <span class="gv">${n(w.C_L)}</span></span>
    </div></div>`;
}

// ---------------------------------------------------------------- diagnostics
function card(title, id) {
  return `<div class="card"><h2>${title}</h2><div id="${id}"></div></div>`;
}
function drawDiag(p) {
  const host = document.getElementById('diag');
  host.innerHTML = '';
  if (!p) return;
  if (p.kind === 'rayleigh_ok') {
    const d = p.diag;
    host.innerHTML = card('Window: slope','g_slopes') + card('Window: |intercept|','g_intercepts')
      + card('Window: R²','g_r_squared') + card('Sensitivity grid','d_sens')
      + card('Lidar constant C_L spread','d_spread');
    const rkm = d.range_km ? f32(d.range_km) : null, hkm = d.half_km ? f32(d.half_km) : null;
    ['slopes','intercepts','r_squared'].forEach(key => {
      const G = (d.grids||{})[key]; if (!G) return;
      const spec = GRID_SPEC[key];
      const a = dec(G.b64,'f32'), [n1,n2] = G.shape, zz = [];
      for (let j = 0; j < n2; j++) { const row = new Array(n1);
        for (let i = 0; i < n1; i++) row[i] = a[i*n2+j]; zz.push(row); }
      const t = { type:'heatmap', z:zz, x:rkm, y:hkm,
        colorbar:{ title:{text:spec.title}, thickness:11,
                   // Plotly's default SI prefixes render 2.5e-13 as "250f" (f = femto), which no
                   // one reads as a number. These grids are tiny-valued, so force exponents.
                   exponentformat: spec.exp ? 'e' : undefined },
        hovertemplate:'centre %{x:.2f} km<br>half %{y:.2f} km<br>%{z:.3g}<extra></extra>' };
      if (spec.diverging) {
        // Symmetric limits, so the neutral colour lands exactly on zero: for the slope, zero IS
        // the good answer (signal parallel to molecular) and it must be findable at a glance.
        let M = 0;
        for (const r of zz) for (const v of r) if (isFinite(v)) M = Math.max(M, Math.abs(v));
        Object.assign(t, { colorscale:'RdBu', reversescale:true, zmin:-M, zmax:M, zmid:0 });
      } else {
        Object.assign(t, { colorscale:spec.scale, zmin:spec.zmin, zmax:spec.zmax });
      }
      const tr = [t];
      if (d.best_km) tr.push({ x:[d.best_km[0]], y:[d.best_km[1]], mode:'markers',
        marker:{symbol:'x', size:11, color:'#d62728', line:{width:2}}, name:'chosen' });
      Plotly.newPlot('g_'+key, tr, Object.assign({}, LAY, { height:270, showlegend:false,
        xaxis:{title:{text:'Centre range (km)'}}, yaxis:{title:{text:'Half-length (km)'}} }), CFG);
    });
    if (d.sens) {
      const a = dec(d.sens.b64,'f32'), [n1,n2] = d.sens.shape, zz = [];
      for (let i = 0; i < n1; i++) { const row = new Array(n2);
        for (let j = 0; j < n2; j++) row[j] = a[i*n2+j]; zz.push(row); }
      Plotly.newPlot('d_sens', [{ type:'heatmap', z:zz,
        x:d.sens.shift.map(v => (v>0?'+':'')+v.toFixed(0)), y:d.sens.lr.map(v => v.toFixed(0)),
        colorscale:'RdBu', reversescale:true, zmin:-d.sens.vabs, zmax:d.sens.vabs,
        colorbar:{title:{text:'Dev (%)'}, thickness:11},
        hovertemplate:'LR %{y} sr<br>shift %{x} m<br>%{z:+.2f} %<extra></extra>' }],
        Object.assign({}, LAY, { height:270, xaxis:{title:{text:'Alt shift (m)'}, type:'category'},
          yaxis:{title:{text:'LR (sr)'}, type:'category'} }), CFG);
    }
    if (d.spread) {
      const s = d.spread, n = s.vals.length;
      Plotly.newPlot('d_spread', [{ x:s.vals.map((_,i)=>1+(n>1?-0.08+0.16*i/(n-1):0)), y:s.vals,
        mode:'markers', marker:{size:7, opacity:0.55, color:'#4c78a8'},
        hovertemplate:'C_L %{y:.4g}<extra></extra>' }],
        Object.assign({}, LAY, { height:270, showlegend:false,
          shapes:[{ type:'rect', xref:'paper', x0:0, x1:1, y0:s.median-s.unc, y1:s.median+s.unc,
                    fillcolor:'rgba(214,39,40,0.12)', line:{width:0}, layer:'below' },
                  { type:'line', xref:'paper', x0:0, x1:1, y0:s.median, y1:s.median,
                    line:{color:'#d62728', width:1.4} }],
          xaxis:{visible:false, range:[0.7,1.3]},
          yaxis:{title:{text:'C_L'}, exponentformat:'e'} }), CFG);
    }
  } else if (p.kind === 'cloud') {
    const d = p.diag;
    // No 'Summary' card: its contents are now the "Why this cloud calibration?" description above.
    host.innerHTML = card('Apparent vs consistent lidar ratio','d_S')
      + card('Coefficient distribution','d_hist');
    const dt = d.x && typeof d.x[0] === 'string';
    Plotly.newPlot('d_S', [
      { x:d.x, y:d.S_app, mode:'markers', name:'apparent S',
        marker:{size:4, color:'#999', opacity:0.5} },
      { x:d.x, y:d.S_con, mode:'markers', name:'consistent S (used)',
        marker:{size:6, color:'#2ca02c'} }],
      Object.assign({}, LAY, { height:270,
        shapes:[{ type:'line', xref:'paper', x0:0, x1:1, y0:d.s_theo, y1:d.s_theo,
                  line:{color:'#d62728', width:1.2, dash:'dash'} }],
        xaxis:{title:{text: dt ? 'Time (UTC)' : 'Hours since start'}, type: dt?'date':'linear'},
        yaxis:{title:{text:'S (sr)'}, range:[0, d.s_ymax]},
        legend:{orientation:'h', y:1.13, font:{size:9}} }), CFG);
    const hs = [];
    if (d.cal_median !== null) {
      if (d.cal_std !== null) hs.push({ type:'rect', yref:'paper', y0:0, y1:1,
          x0:d.cal_median-d.cal_std, x1:d.cal_median+d.cal_std,
          fillcolor:'rgba(214,39,40,0.12)', line:{width:0}, layer:'below' });
      hs.push({ type:'line', yref:'paper', y0:0, y1:1, x0:d.cal_median, x1:d.cal_median,
          line:{color:'#d62728', width:1.4} });
    }
    Plotly.newPlot('d_hist', [{ type:'histogram', x:d.coeffs, marker:{color:'#2ca02c'},
        opacity:0.75, nbinsx:Math.min(30, Math.max(5, Math.floor(d.coeffs.length/2))) }],
      Object.assign({}, LAY, { height:270, shapes:hs, showlegend:false,
        xaxis:{title:{text:'coefficient C (per profile)'}}, yaxis:{title:{text:'count'}} }), CFG);
  } else if (p.kind === 'rayleigh_fail') {
    host.innerHTML = '<div class="card" style="grid-column:1/-1"><p class="empty">' +
      'The window search never ran — the night was rejected before a fit was attempted, so there ' +
      'is no slope / intercept / R² grid and no sensitivity to show. The curtain and the ' +
      'night-mean profile above are the whole diagnostic.</p></div>';
  }
}

// ---------------------------------------------------------------- messages + flags
function drawMsg(p) {
  const host = document.getElementById('msg');
  if (!p) { host.innerHTML = ''; return; }
  const ok = isOK(p);
  const bits = [];
  if (p.kind === 'none')
    bits.push('<b>No diagnostic figure is produced for this rejection.</b>');
  else if (p.kind === 'rayleigh_fail')
    bits.push('<b>NOT CALIBRATED — ' + (p.reason || 'rejected') + '</b>');
  else if (ok) bits.push('<b>Calibrated</b>');
  else bits.push('<b>No calibration constant was produced</b>');
  const det = [];
  if (p.flag !== undefined && p.flag !== null)
    det.push('flag <b>' + p.flag + '</b>' + (p.flag_label ? ' — ' + p.flag_label : ''));
  if (p.message) det.push('message: “' + p.message + '”');
  if (p.constant) det.push('C_L = ' + (+p.constant).toPrecision(5) +
    (p.uncertainty ? ' ± ' + (+p.uncertainty).toPrecision(3) : ''));
  if (p.profile && p.profile.mean_of) det.push(p.profile.mean_of);
  host.className = 'msg ' + (ok ? 'ok' : 'bad');
  host.innerHTML = bits.join('') + (det.length ? '<div class="det">' + det.join(' · ') + '</div>' : '');
}
function drawFlags() {
  const host = document.getElementById('flags');
  host.innerHTML = D.methods.map(m => {
    const p = summaryOf(curDate, m);
    const lbl = m === 'cloud' ? 'Cloud' : 'Rayleigh';
    if (!p) return `<div style="margin:4px 0"><span class="chip na">${lbl}: not run</span></div>`;
    const cls = isOK(p) ? 'ok' : 'bad';
    return `<div style="margin:4px 0"><span class="chip ${cls}">${lbl}: flag ${p.flag}</span>` +
      (p.flag_label ? `<div style="font-size:11px;color:#66707a;margin:2px 0 0 2px">${p.flag_label}</div>` : '') +
      `</div>`;
  }).join('');
}

// ---------------------------------------------------------------- render
function render() {
  const i = D.dates.indexOf(curDate);
  document.getElementById('prev').disabled = i <= 0;
  document.getElementById('next').disabled = i < 0 || i >= D.dates.length - 1;
  document.getElementById('datelbl').textContent =
    curDate.slice(0,4) + '-' + curDate.slice(4,6) + '-' + curDate.slice(6,8);
  buildMethodSwitch();
  const p = dayOf(curDate, curMethod);
  // On the lazy station dashboard the index knows this day exists before its payload is fetched.
  // Ask the host page for it; the plots below simply draw empty until it lands.
  if (!p && summaryOf(curDate, curMethod) && window.__onMissingDay) {
    window.__onMissingDay(curDate, curMethod);
  }
  document.getElementById('cst').textContent =
    p && p.constant ? 'C_L = ' + (+p.constant).toPrecision(5) : '';
  drawPanel(p);                    // curtain + profile, one figure, shared y
  drawWhy(p);
  drawDiag(p);
  drawMsg(p);
  drawFlags();
  markCal();
}
document.getElementById('prev').addEventListener('click', () => {
  const i = D.dates.indexOf(curDate); if (i > 0) { curDate = D.dates[i-1]; render(); } });
document.getElementById('next').addEventListener('click', () => {
  const i = D.dates.indexOf(curDate); if (i < D.dates.length-1) { curDate = D.dates[i+1]; render(); } });
document.addEventListener('keydown', e => {
  if (e.target.tagName === 'INPUT') return;
  if (e.key === 'ArrowLeft')  document.getElementById('prev').click();
  if (e.key === 'ArrowRight') document.getElementById('next').click();
  if (e.key === 'm') { const a = D.methods.filter(m => summaryOf(curDate, m));
    if (a.length > 1) { curMethod = a[(a.indexOf(curMethod)+1) % a.length]; render(); } }
});
buildLegend();
buildCal();
render();
"""

HTML = ("""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Daily calibration panel — mockup</title>
<script>__PLOTLY__</script>
<style>""" + PANEL_CSS + """</style></head><body>
""" + PANEL_HDR + "\n" + PANEL_BODY + """
<script id="payload" type="application/json">__PAYLOAD__</script>
<script>""" + PANEL_JS + """</script></body></html>""")


if __name__ == "__main__":
    main()
