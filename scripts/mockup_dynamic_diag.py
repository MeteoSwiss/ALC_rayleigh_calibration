#!/usr/bin/env python3
"""Feasibility test: can the per-night diagnostics be CLIENT-SIDE plots instead of server PNGs?

The per-night diagnostic image is the expensive product in this pipeline and the reason the release
recompute runs with PLOTS=0. This script answers whether the same information can be shipped as
DATA and drawn in the browser, for all THREE diagnostics the pipeline emits:

  rayleigh_ok    the 9-panel success dashboard (3 molecular profiles, 3 window grids, the annotated
                 RCS curtain, the sensitivity grid and the C_L spread)
  rayleigh_fail  the 2-panel rejection figure, whose ONLY marker of failure is a dark-red suptitle
                 "NOT CALIBRATED: <reason>" -- the numeric flag is never drawn on the PNG at all
  cloud          the 5-panel liquid-cloud (O'Connor) dashboard

It measures rather than assumes: it times one real night with images off and on, wraps the
production plotting functions so the payload is built from EXACTLY the arrays the PNG was drawn
from, and renders those arrays with Plotly into a standalone page.

Nothing is deployed: this writes a mockup and prints a comparison.

  python scripts/mockup_dynamic_diag.py --key 0-208-0-06011_A --date 20260405 \
      --l1-root D:/E-PROFILE_L1_2026 --out mock.html
  python scripts/mockup_dynamic_diag.py --key 0-20000-0-06610_C --date 20260701 \
      --method cloud --l1-root D:/E-PROFILE_L1_2026 --out mock_cloud.html
"""
from __future__ import annotations

import argparse
import base64
import gzip
import json
import os
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

CAPTURED: dict = {}
#: emit the opt-in native-resolution curtain alongside the default one (--no-hires to skip)
HIRES = True


def _b64(a: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")


def _f32(a) -> str:
    return _b64(np.asarray(a, dtype="float32"))


def _jsonable(a) -> list:
    """float list with NaN -> None (JSON has no NaN)."""
    return [None if not np.isfinite(v) else round(float(v), 6) for v in np.asarray(a, dtype=float)]


def _quant_u8(m: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """2-D log-ish field -> uint8, 255 = missing. Same trick the v3 dashboard uses for its
    time-height curtains: viridis has 256 entries anyway, so this is visually lossless."""
    q = np.full(m.shape, 255, dtype=np.uint8)
    ok = np.isfinite(m)
    if ok.any():
        v = np.clip((m[ok] - lo) / max(hi - lo, 1e-9), 0, 1)
        q[ok] = np.clip(np.round(v * 254.0), 0, 254).astype(np.uint8)
    return q


#: Curtain grid, chosen to match what the PNG itself resolves. plot_rayleigh_diagnostics_compact
#: strides the mesh to <=800x800 (plotting.py:1095-1096), so at this grid the client-side version
#: shows EXACTLY what the image shows and loses nothing relative to today's product. Measured on a
#: real CL31 night (5759x1024): 265 KB gzipped vs a 2542 KB PNG -- still 10x cheaper. Going finer is
#: pointless (the PNG cannot show it) and full native resolution is 4186 KB gzipped, i.e. LARGER
#: than the image it replaces, which is why "ship everything and let the browser zoom" is not an
#: option. True zoom-in beyond the PNG would need an opt-in second fetch of a full-resolution
#: payload for the one night being interrogated.
PNG_GRID_T = 800
PNG_GRID_R = 800


def _block_mean(a: np.ndarray, bt: int, br: int) -> np.ndarray:
    """Block-average a (n_t, n_r) array by (bt, br), ignoring NaN.

    Averaging is done by the CALLER in LINEAR space, before the log10: the mean of logs is a
    geometric mean, which is biased low on noisy data and would darken the curtain. A block that is
    entirely NaN stays NaN, so genuine dropouts remain gaps instead of being filled in.
    The tail rows/columns that do not fill a whole block are averaged over what they have.
    """
    if bt == 1 and br == 1:
        return a
    n_t = int(np.ceil(a.shape[0] / bt)) * bt
    n_r = int(np.ceil(a.shape[1] / br)) * br
    pad = np.full((n_t, n_r), np.nan, dtype=float)
    pad[:a.shape[0], :a.shape[1]] = a
    blocks = pad.reshape(n_t // bt, bt, n_r // br, br)
    with warnings.catch_warnings():                 # an all-NaN block is expected, not an error
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(blocks, axis=(1, 3))


def _curtain(mat, rng_m, st_target_t=PNG_GRID_T, st_target_r=PNG_GRID_R,
             vmin=None, vmax=None, crop_km=None, reduce="mean") -> dict:
    """(n_time, n_range) -> quantised, strided curtain block.

    Colour limits are the 5th/95th percentile of the FINITE log10 values computed on the FULL field,
    exactly as the PNG does, and they travel with the block: recomputing them from a cropped or
    strided field would shift every colour.

    *crop_km* drops range gates the panel never displays (the cloud diagnostic fixes its y-axis at
    ~3.4 km while its pcolormesh is full-resolution, so most gates are off-screen). Cropping happens
    AFTER the percentiles are taken, so the colours still match the PNG.
    """
    m = np.asarray(mat, dtype=float)
    rng = np.asarray(rng_m, dtype=float)
    lin = np.where(m > 0, m, np.nan)                # non-positive -> gap, as the PNG does
    with np.errstate(all="ignore"):
        fin = np.log10(lin[np.isfinite(lin)])
    lo = vmin if vmin is not None else (float(np.percentile(fin, 5)) if fin.size else 0.0)
    hi = vmax if vmax is not None else (float(np.percentile(fin, 95)) if fin.size else 6.0)
    if crop_km is not None:
        keep = rng <= float(crop_km) * 1e3
        if keep.any():
            lin, rng = lin[:, keep], rng[keep]
    st_t = max(1, int(np.ceil(lin.shape[0] / st_target_t)))
    st_r = max(1, int(np.ceil(lin.shape[1] / st_target_r)))
    if reduce == "mean":
        red = _block_mean(lin, st_t, st_r)          # average IN LINEAR SPACE, then log
        r_km = _block_mean(rng.reshape(1, -1), 1, st_r).ravel()
    else:
        red, r_km = lin[::st_t, ::st_r], rng[::st_r]
    with np.errstate(all="ignore"):
        sub = np.log10(red)
    return {"b64": _b64(_quant_u8(sub.T, lo, hi)), "shape": list(sub.T.shape), "lo": lo, "hi": hi,
            "range_km": _f32(r_km * 1e-3), "st_t": st_t, "st_r": st_r, "reduce": reduce}


def _runs(mask: np.ndarray) -> list:
    """Contiguous True runs as [start, end_inclusive]. Merging them is what makes the bands read
    as regions rather than as a comb of per-profile slivers (mirrors the PNG's _hatched_runs)."""
    idx = np.where(mask)[0]
    if idx.size == 0:
        return []
    return [[int(r[0]), int(r[-1])] for r in np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)]


# =============================================================================== payload builders
def _payload_rayleigh_ok(kw: dict) -> dict:
    """The 9-panel success dashboard. Written against calibration/plotting.py:812-1052; an earlier
    version of this silently dropped the retrieved beta_att, the cloud base, the excluded regions
    and two whole panels, so every element below is deliberate."""
    out = {"kind": "rayleigh_ok"}
    rng = np.asarray(kw.get("range_alc", []), dtype=float)
    alt = float(kw.get("altitude", 0.0) or 0.0)
    z_km = (rng + alt) * 1e-3
    SCALE = 1e6

    sig_n = np.asarray(kw.get("signal_normalized", []), dtype=float)
    prof = {"rcs_mean": np.asarray(kw.get("rcs_mean", []), dtype=float),
            "signal_normalized": sig_n * SCALE,
            "p_mol": np.asarray(kw.get("p_mol", []), dtype=float) * SCALE,
            "beta_att_mol": np.asarray(kw.get("beta_att_mol", []), dtype=float) * SCALE,
            # DERIVED in the plot, never passed in -- this is the "retrieved beta_att" curve
            "beta_att": (sig_n * rng ** 2) * SCALE if sig_n.size == rng.size else np.array([])}
    out["prof"] = {k: _f32(v) for k, v in prof.items() if v.size}
    out["z_km"] = _f32(z_km)
    z0 = float(kw.get("fit_altitude_start", 0.0) or 0.0) * 1e-3
    z1 = float(kw.get("fit_altitude_end", 0.0) or 0.0) * 1e-3
    out["fit_alt_km"] = [z0, z1]
    out["z_max_km"] = float(min(z1 + 3.0, float(z_km.max()) if z_km.size else z1 + 3.0))

    grids = {}
    for name in ("slopes", "intercepts", "r_squared"):
        v = kw.get(name)
        if v is None:
            continue
        a = np.abs(np.asarray(v, float)) if name == "intercepts" else np.asarray(v, float)
        grids[name] = {"b64": _f32(a), "shape": list(a.shape)}
    out["grids"] = grids
    for name, key in (("range_bin_m", "range_km"), ("half_length_m", "half_km")):
        if kw.get(name) is not None:
            out[key] = _f32(np.asarray(kw[name], float) * 1e-3)
    br, bh = kw.get("best_range_m"), kw.get("best_half_m")
    if br is not None and bh is not None:
        out["best_km"] = [float(br) * 1e-3, float(bh) * 1e-3]

    cl_m, cl_med = kw.get("cl_matrix"), kw.get("cl_median")
    if cl_m is not None and cl_med:
        cm = np.asarray(cl_m, float)
        rel = (cm - float(cl_med)) / float(cl_med) * 100.0
        out["sens"] = {"b64": _f32(rel), "shape": list(rel.shape),
                       "vabs": float(max(np.nanmax(np.abs(rel)) if np.isfinite(rel).any() else 1.0, 1.0)),
                       "lr": _jsonable(kw.get("lr_values", [])),
                       "shift": _jsonable(kw.get("alt_shifts", []))}
        out["cl_spread"] = {"vals": _jsonable(cm[np.isfinite(cm)].ravel()),
                            "median": float(cl_med),
                            "unc": float(kw.get("cl_uncertainty", 0.0) or 0.0)}

    rcs = kw.get("rcs")
    if rcs is not None:
        m = np.asarray(rcs, float)
        rr = np.asarray(kw.get("rcs_range_alc", rng), float)
        cur = _curtain(m, rr)
        st_t = cur["st_t"]
        tdt = kw.get("time_datetime")
        if tdt is not None and len(tdt):
            cur["time"] = [str(t)[:19] for t in np.asarray(tdt)[::st_t]]
        else:
            hrs = np.asarray(kw.get("hours_since_start", []), float)
            cur["hours"] = _jsonable(hrs[::st_t]) if hrs.size else []

        n_t = int(m.shape[0])
        cbh0 = None
        if kw.get("cloud_base_height") is not None:
            c = np.asarray(kw["cloud_base_height"], float)
            c = c[:, 0] if c.ndim > 1 else c
            ncv = float(kw.get("no_cloud_value", -9.0))
            cbh0 = np.where((c == ncv) | (c <= 0), np.nan, c)      # mask FIRST, then stride
            cur["cbh_km"] = _jsonable(cbh0[::st_t] * 1e-3)
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
        not_used = ~used
        cur["excl_lowcloud"] = _runs(not_used[::st_t] & flagged[::st_t])
        cur["excl_screened"] = _runs(not_used[::st_t] & ~flagged[::st_t])
        if cbh0 is not None:
            hc = used & np.isfinite(cbh0) & ~flagged
            y_hi = float(rr.max()) * 1e-3
            # EXACTLY the PNG's expression (plotting.py:1030):
            #   y_lo = np.where(high_cloud, clip((cbh - 500 m), 0, y_hi), y_hi)
            # Outside the mask y_lo == y_hi, i.e. a ZERO-HEIGHT band. Emitting null there instead
            # looks equivalent but is not: Plotly's fill:"tonexty" bridges across nulls and joins
            # distant points into spurious wedges (the operator saw triangles).
            cur["high_cloud_ylo"] = [round(float(np.clip((c - 500.0) * 1e-3, 0.0, y_hi)), 4)
                                     if (h and np.isfinite(c)) else round(y_hi, 4)
                                     for h, c in zip(hc[::st_t], cbh0[::st_t])]
            cur["has_high_cloud"] = bool(hc.any())
            cur["y_hi_km"] = y_hi
        if br is not None and bh is not None:
            cur["mol_layer_km"] = [(float(br) - float(bh)) * 1e-3, (float(br) + float(bh)) * 1e-3]
        out["curtain"] = cur
        if HIRES:
            # The opt-in payload: NATIVE resolution, no striding. This is the one that only makes
            # sense fetched on demand for the single night being interrogated -- it is larger than
            # the PNG it replaces, so it can never be the default.
            hi = _curtain(m, rr, st_target_t=m.shape[0], st_target_r=m.shape[1],
                          vmin=cur["lo"], vmax=cur["hi"])
            hi["time"] = ([str(t)[:19] for t in np.asarray(tdt)] if (tdt is not None and len(tdt))
                          else None)
            if hi["time"] is None:
                hrs = np.asarray(kw.get("hours_since_start", []), float)
                hi["hours"] = _jsonable(hrs) if hrs.size else []
            out["curtain_hi"] = hi
    out["title"] = str(kw.get("title", ""))
    return out


def _payload_rayleigh_fail(kw: dict) -> dict:
    """The 2-panel rejection figure (calibration/plotting.py:1053-1162).

    The PNG communicates the rejection with NOTHING but a dark-red suptitle
    "<title>  —  NOT CALIBRATED: <reason>". The numeric flag never appears on it, so the payload
    carries the flag and its label explicitly -- the operator asked to still see the flag result.
    """
    out = {"kind": "rayleigh_fail", "reason": str(kw.get("reason", "")),
           "title": str(kw.get("title", "")), "banner_color": "#b00020"}
    rng = np.asarray(kw.get("range_alc", []), float)
    alt = float(kw.get("altitude", 0.0) or 0.0)
    z_km = (rng + alt) * 1e-3

    win = kw.get("molecular_window")
    has_window = win is not None and np.all(np.isfinite(np.asarray(win, float)))
    if has_window:
        w = np.asarray(win, float)
        out["window_range_km"] = [float(w[0]) * 1e-3, float(w[1]) * 1e-3]          # ax_r: AGL
        out["window_alt_km"] = [float(w[0] + alt) * 1e-3, float(w[1] + alt) * 1e-3]  # ax_p: ASL

    rcs = kw.get("rcs")
    if rcs is not None:
        m = np.asarray(rcs, float)
        cur = _curtain(m, rng)
        st_t = cur["st_t"]
        tdt = kw.get("time_datetime")
        if tdt is not None and len(tdt):
            cur["time"] = [str(t)[:19] for t in np.asarray(tdt)[::st_t]]
        else:
            hrs = np.asarray(kw.get("hours_since_start", []), float)
            cur["hours"] = _jsonable(hrs[::st_t]) if hrs.size else []
        cbh = kw.get("cbh")
        if cbh is not None:
            c = np.asarray(cbh, float)
            c = c[:, 0] if c.ndim > 1 else c
            ncv = float(kw.get("no_cloud_value", -9.0))
            c = np.where((c == ncv) | (c <= 0), np.nan, c)
            if np.isfinite(c).any():
                cur["cbh_km"] = _jsonable(c[::st_t] * 1e-3)
        cur["y_max_km"] = float(rng.max()) * 1e-3 if rng.size else None
        out["curtain"] = cur
        if HIRES:
            hi = _curtain(m, rng, st_target_t=m.shape[0], st_target_r=m.shape[1],
                          vmin=cur["lo"], vmax=cur["hi"])
            if tdt is not None and len(tdt):
                hi["time"] = [str(t)[:19] for t in np.asarray(tdt)]
            else:
                hrs = np.asarray(kw.get("hours_since_start", []), float)
                hi["hours"] = _jsonable(hrs) if hrs.size else []
            out["curtain_hi"] = hi
        # the profile panel is computed on the FULL matrix, not the strided one
        with np.errstate(all="ignore"):
            rcs_mean = np.nanmean(m, axis=0)
            sig = rcs_mean / (rng ** 2)
        finite = np.isfinite(sig)
        out["has_profile"] = bool(finite.any())
        if finite.any():
            out["sig"] = _f32(np.where(finite, sig, np.nan))
            out["z_km"] = _f32(z_km)
            out["z_max_p_km"] = float(z_km.max()) if np.isfinite(z_km).any() else 1.0
            pm = kw.get("p_mol")
            if pm is not None:
                pm = np.asarray(pm, float)
                r0 = float(kw.get("range_start_m", 2000.0))
                r1 = float(kw.get("range_end_m", 6000.0))
                band = (rng >= r0) & (rng <= r1) & finite & np.isfinite(pm)
                if band.any() and np.nanmedian(pm[band]) > 0:
                    scale = float(np.nanmedian(sig[band]) / np.nanmedian(pm[band]))
                    out["pm_scaled"] = _f32(pm * scale)   # ship the PRODUCT, not p_mol + the band
    return out


def _payload_cloud(data, res, title: str) -> dict:
    """The 5-panel liquid-cloud dashboard (calibration/plotting.py:78-236).

    A REJECTED scene (flags -20..-26) goes through the same code path with an empty result, so
    `sel` is all-False and every green element disappears -- that absence IS the rejection signal,
    and the summary text legitimately prints "nan". Both are reproduced rather than tidied away.
    """
    out = {"kind": "cloud", "title": str(title or "")}
    t = np.asarray(data.time)
    rng = np.asarray(data.range, float)
    beta = np.asarray(data.beta, float)          # (n_range, n_time), already WV-corrected in place
    n_time = int(t.size)
    try:
        hrs = (t.astype("datetime64[s]").astype("float64")
               - t[0].astype("datetime64[s]").astype("float64")) / 3600.0
    except Exception:                            # noqa: BLE001 - mirrors the plot's own fallback
        hrs = np.arange(n_time, dtype=float)

    g = lambda n, d=None: getattr(res, n, d)     # noqa: E731
    S_app = np.asarray(g("S_apparent"), float) if g("S_apparent") is not None else np.full(n_time, np.nan)
    S_con = np.asarray(g("S_consistent"), float) if g("S_consistent") is not None else np.full(n_time, np.nan)
    coeffs = np.asarray(g("all_coefficients"), float) if g("all_coefficients") is not None else np.full(n_time, np.nan)
    sel = np.isfinite(S_con) if S_con.size == n_time else np.zeros(n_time, bool)
    cfg = g("config")
    cal_lo = float(getattr(cfg, "cal_minheight", 100.0)) if cfg else 100.0
    cal_hi = float(getattr(cfg, "cal_maxheight", 2400.0)) if cfg else 2400.0
    valid_c = coeffs[np.isfinite(coeffs)]
    cal_med = g("cal_median", float("nan"))
    s_theo = (float(np.nanmedian(S_con[sel])) / cal_med) if (sel.any() and cal_med) else 18.8

    # This panel's y axis is fixed at min(cal_hi/1000 + 1, top gate) ~ 3.4 km, so every gate above
    # that is drawn by matplotlib and then clipped away. Cropping there is free fidelity.
    y_max_km = float(min(cal_hi * 1e-3 + 1.0, float(rng.max()) * 1e-3))
    cur = _curtain(beta.T, rng, crop_km=y_max_km)   # _curtain wants (n_time, n_range)
    st_t = cur["st_t"]
    cur["hours"] = _jsonable(hrs[::st_t])
    out["curtain"] = cur
    out["y_max_km"] = y_max_km
    if HIRES:
        hi = _curtain(beta.T, rng, st_target_t=beta.shape[1], st_target_r=beta.shape[0],
                      vmin=cur["lo"], vmax=cur["hi"], crop_km=y_max_km)
        hi["hours"] = _jsonable(hrs)
        out["curtain_hi"] = hi
    out["gate_km"] = [cal_lo * 1e-3, cal_hi * 1e-3]

    cbh = np.asarray(getattr(data, "cbh", []), float)
    if cbh.size == n_time and np.isfinite(cbh).any():
        ok = cbh > 0
        out["cbh_unsel_km"] = _jsonable(np.where(ok & ~sel, cbh, np.nan)[::st_t] * 1e-3)
        out["cbh_sel_km"] = _jsonable(np.where(ok & sel, cbh, np.nan)[::st_t] * 1e-3)
    out["sel_runs"] = _runs(sel[::st_t])         # green bands over the profiles actually used
    out["n_sel"] = int(sel.sum())

    if sel.any():
        with np.errstate(all="ignore"):
            out["prof"] = _f32(np.nanmean(beta[:, sel], axis=1))
        out["rng_km"] = _f32(rng * 1e-3)
        med_cbh = float(np.nanmedian(cbh[sel])) if cbh.size == n_time else float("nan")
        if np.isfinite(med_cbh):
            out["med_cbh_km"] = med_cbh * 1e-3
            out["b_layer_km"] = [med_cbh * 1e-3, (med_cbh + 300.0) * 1e-3]

    out["hrs"] = _jsonable(hrs[::st_t])
    out["S_app"] = _jsonable(S_app[::st_t])
    out["S_con"] = _jsonable(np.where(sel, S_con, np.nan)[::st_t])
    out["s_theo"] = float(s_theo)
    upper = s_theo * 3.0
    if sel.any() and np.any(np.isfinite(S_con[sel])):
        upper = max(upper, 1.3 * float(np.nanmax(S_con[sel])))
    out["s_ymax"] = float(upper)
    out["coeffs"] = _jsonable(valid_c)
    out["cal_median"] = None if not np.isfinite(cal_med) else float(cal_med)
    cal_std = g("cal_std", float("nan"))
    out["cal_std"] = None if not np.isfinite(cal_std) else float(cal_std)

    # the Summary panel, reproduced verbatim (including the literal "nan" a rejected night prints)
    lines = [f"C_L (Wiegner)   = {g('lidar_constant', float('nan')):.4g}",
             f"coefficient C   = {cal_med:.4g}",
             f"std(C)          = {cal_std:.3g}",
             (f"rel. unc        = {100 * cal_std / cal_med:.1f} %" if cal_med else "rel. unc = —"),
             f"n profiles used = {g('n_profiles', 0)}",
             f"theoretical S   = {s_theo:.2f} sr"]
    for name, d in (("instrument filter", g("filter_stats")), ("cloud filter", g("cloud_stats")),
                    ("consistency", g("consistency_stats"))):
        if isinstance(d, dict) and d:
            lines.append(f"{name} rejections:")
            for k, v in d.items():
                lines.append(f"   {k}: {v}")
    out["summary"] = lines
    return out


# ==================================================================================== capture
def _install_capture():
    """Wrap the three production plotting functions: time the real render AND keep the arguments,
    so the payload is provably built from the same arrays the PNG was drawn from."""
    import calibration.plotting as P
    import calibration.rayleigh.calibration as RC

    def wrap(real, kind, is_cloud=False):
        def w(*a, **kw):
            t0 = time.perf_counter()
            fig = real(*a, **kw)
            CAPTURED.setdefault("renders", []).append(time.perf_counter() - t0)
            CAPTURED["kind"] = kind
            CAPTURED["args"], CAPTURED["kwargs"] = a, kw
            sp = kw.get("save_path") or (a[3] if is_cloud and len(a) > 3 else None)
            if sp and Path(sp).exists():
                CAPTURED["png"] = Path(sp)
            return fig
        return w

    P.plot_rayleigh_diagnostics_compact = wrap(P.plot_rayleigh_diagnostics_compact, "rayleigh_ok")
    RC.plot_rayleigh_diagnostics_compact = P.plot_rayleigh_diagnostics_compact
    P.plot_rayleigh_diagnostics_failure = wrap(P.plot_rayleigh_diagnostics_failure, "rayleigh_fail")
    RC.plot_rayleigh_diagnostics_failure = P.plot_rayleigh_diagnostics_failure
    P.plot_cloud_diagnostics_compact = wrap(P.plot_cloud_diagnostics_compact, "cloud", True)


def _time_runner(key: str, date: str, plots: str, method: str) -> float:
    env = dict(os.environ, PLOTS=plots)
    t0 = time.perf_counter()
    subprocess.run([sys.executable, str(REPO / "scripts/run_network_calibration.py"),
                    "--stream", key, "--start", date, "--end", date, "--methods", method,
                    "--force", "--per-type", "0", "--ignore-coverage", "--workers", "1"],
                   env=env, cwd=str(REPO), capture_output=True, text=True, timeout=1800)
    return time.perf_counter() - t0


BODY = {
    "rayleigh_ok": """
<h3>Molecular panels <span class="win">(gold band = Rayleigh fit window)</span></h3>
<div class="grid"><div class="card"><div id="p1"></div></div>
<div class="card"><div id="p2"></div></div><div class="card"><div id="p3"></div></div></div>
<h3>Window search <span class="win">(red &times; = chosen window)</span></h3>
<div class="grid"><div class="card"><div id="g_slopes"></div></div>
<div class="card"><div id="g_intercepts"></div></div>
<div class="card"><div id="g_r_squared"></div></div></div>
<h3>Range-corrected signal &mdash; full night
  <button class="resbtn" id="res-btn">full resolution</button>
  <span class="res" id="res"></span></h3>
<div class="card"><div id="curtain"></div></div>
<h3>Sensitivity</h3>
<div class="grid2"><div class="card"><div id="sens"></div></div>
<div class="card"><div id="spread"></div></div></div>""",
    "rayleigh_fail": """
<div class="banner" id="banner"></div>
<h3>Range-corrected signal &mdash; full night
  <button class="resbtn" id="res-btn">full resolution</button>
  <span class="res" id="res"></span></h3>
<div class="grid23"><div class="card"><div id="curtain"></div></div>
<div class="card"><div id="failprof"></div></div></div>""",
    "cloud": """
<h3>Attenuated backscatter
  <button class="resbtn" id="res-btn">full resolution</button>
  <span class="res" id="res"></span></h3>
<div class="grid23"><div class="card"><div id="cl_curtain"></div></div>
<div class="card"><div id="cl_prof"></div></div></div>
<div class="grid"><div class="card"><div id="cl_S"></div></div>
<div class="card"><div id="cl_hist"></div></div>
<div class="card"><pre id="cl_summary" class="summary"></pre></div></div>""",
}

HTML = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Dynamic diagnostic — feasibility mockup</title>
<script>{plotly}</script>
<style>
 body {{ font: 14px/1.5 -apple-system,"Segoe UI",Roboto,Arial,sans-serif; margin: 22px 28px;
        color:#1a2530; }}
 .hdr {{ background:#fff3cd; border:1px solid #ffe08a; border-radius:8px; padding:10px 14px;
         margin-bottom:16px; }}
 .banner {{ background:#fdecef; border:1px solid #f5c2c7; border-left:5px solid #b00020;
            color:#b00020; font-weight:600; border-radius:8px; padding:10px 14px; margin:12px 0; }}
 .banner .flag {{ font-weight:400; color:#7a2030; display:block; margin-top:3px; font-size:13px; }}
 .grid {{ display:grid; grid-template-columns:1fr 1fr 1fr; gap:14px; }}
 .grid2 {{ display:grid; grid-template-columns:1fr 1fr; gap:14px; }}
 .grid23 {{ display:grid; grid-template-columns:2fr 1fr; gap:14px; }}
 .card {{ border:1px solid #dbe3ea; border-radius:10px; padding:8px; }}
 .summary {{ font:12px/1.45 ui-monospace,"Cascadia Code",Consolas,monospace; margin:6px 8px;
             white-space:pre; }}
 table {{ border-collapse:collapse; margin:14px 0; }}
 td,th {{ border:1px solid #dbe3ea; padding:5px 10px; text-align:right; }}
 th:first-child, td:first-child {{ text-align:left; }}
 .win {{ font-weight:600; color:#b7791f; }}
 .resbtn {{ font:inherit; font-size:12px; font-weight:600; margin-left:14px; padding:3px 10px;
            border:1px solid #2a5a82; border-radius:6px; background:#eef4fb; color:#0b3d61;
            cursor:pointer; }}
 .resbtn:hover {{ background:#dceaf7; }}
 .res {{ font-size:12px; font-weight:400; color:#66707a; margin-left:10px; }}
 .res.on {{ color:#0b7285; font-weight:600; }}
</style></head><body>
<div class="hdr"><b>FEASIBILITY MOCKUP</b> ({kind}) — every plot below is drawn <b>in your
browser</b> from a {size_kb:.0f} KB payload ({gz_kb:.0f} KB gzipped) captured from the real
calibration of <b>{key} {date}</b>. No PNG is involved. The equivalent server-rendered image is
{png_kb:.0f} KB and took <b>{render_s:.2f} s</b> of matplotlib time.</div>
<table>
<tr><th>quantity</th><th>server PNG</th><th>client-side payload</th></tr>
<tr><td>bytes on the wire</td><td>{png_kb:.0f} KB</td><td>{gz_kb:.0f} KB (gzip)</td></tr>
<tr><td>CPU to produce</td><td>{render_s:.2f} s</td><td>{dump_ms:.0f} ms</td></tr>
<tr><td>night with images vs without</td><td colspan="2">{ab}</td></tr>
</table>
{body}
<script id="payload" type="application/json">{payload}</script>
<script>
const D = JSON.parse(document.getElementById('payload').textContent);
const dec = (b64, type) => {{
  const s = atob(b64), u = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i);
  return type === 'u8' ? u : new Float32Array(u.buffer);
}};
const LAY = {{ margin:{{l:64,r:12,t:32,b:46}}, height:320, template:'plotly_white',
               font:{{size:11}} }};
const CFG = {{responsive:true}};
// uint8 curtain -> 2-D array of physical log10 values (255 = missing -> null)
function grid(c) {{
  const q = dec(c.b64, 'u8'), [ny, nx] = c.shape, zz = [];
  for (let j = 0; j < ny; j++) {{
    const row = new Array(nx);
    for (let i = 0; i < nx; i++) {{
      const v = q[j * nx + i];
      row[i] = v === 255 ? null : c.lo + (v / 254) * (c.hi - c.lo);
    }}
    zz.push(row);
  }}
  return zz;
}}
const xOf = c => c.time || c.hours;
const xIsDate = c => !!c.time;

// --------------------------------------------------------------- adaptive resolution
// The default curtain matches the PNG's own grid (<=800x800), so nothing is lost relative to
// today's image. The native-resolution block is an OPT-IN second payload: it is larger than the
// PNG it replaces, so it can never be the default -- but once the operator zooms in, only the
// visible columns are decoded, which is cheap. In production this second block is a separate
// fetch; in this mockup it is embedded so the file stays self-contained.
function attachHiRes(divId, cLo, cHi, label) {{
  const gd = document.getElementById(divId);
  if (!cHi) return;
  const xLo = xOf(cLo), xHi = xOf(cHi);
  const asNum = v => (typeof v === 'string' ? Date.parse(v) : +v);
  const numHi = xHi.map(asNum);
  const full = [asNum(xLo[0]), asNum(xLo[xLo.length - 1])];
  const rkHi = Array.from(dec(cHi.range_km, 'f32'));
  const qHi = dec(cHi.b64, 'u8');
  const [nyHi, nxHi] = cHi.shape;
  let mode = 'png';

  function sliceHi(x0, x1) {{               // decode ONLY the visible columns
    let i0 = 0, i1 = nxHi - 1;
    while (i0 < nxHi - 1 && numHi[i0] < x0) i0++;
    while (i1 > 0 && numHi[i1] > x1) i1--;
    i0 = Math.max(0, i0 - 1); i1 = Math.min(nxHi - 1, i1 + 1);
    const w = i1 - i0 + 1, zz = [];
    for (let j = 0; j < nyHi; j++) {{
      const row = new Array(w);
      for (let i = 0; i < w; i++) {{
        const v = qHi[j * nxHi + (i0 + i)];
        row[i] = v === 255 ? null : cHi.lo + (v / 254) * (cHi.hi - cHi.lo);
      }}
      zz.push(row);
    }}
    return {{ z: zz, x: xHi.slice(i0, i1 + 1), y: rkHi, n: w }};
  }}
  function setPNG() {{
    if (mode === 'png') return;
    Plotly.restyle(gd, {{ z:[grid(cLo)], x:[xLo], y:[Array.from(dec(cLo.range_km,'f32'))] }}, [0]);
    mode = 'png'; badge();
  }}
  function setHi(x0, x1) {{
    const s = sliceHi(x0, x1);
    Plotly.restyle(gd, {{ z:[s.z], x:[s.x], y:[s.y] }}, [0]);
    mode = 'hi'; badge(s.n);
  }}
  function badge(n) {{
    const el = document.getElementById(label);
    if (!el) return;
    el.textContent = mode === 'png'
      ? `resolution: PNG grid (${{cLo.shape[1]}} x ${{cLo.shape[0]}}) — zoom in for native`
      : `resolution: NATIVE (${{n}} x ${{nyHi}} in view, of ${{nxHi}} x ${{nyHi}})`;
    el.className = mode === 'png' ? 'res' : 'res on';
  }}
  badge();

  gd.on('plotly_relayout', function (e) {{
    if (e['xaxis.autorange'] || e['autosize']) {{ setPNG(); return; }}
    const a = e['xaxis.range[0]'], b = e['xaxis.range[1]'];
    if (a === undefined || b === undefined) return;
    const x0 = asNum(a), x1 = asNum(b);
    // switch once the view is a small enough slice that native columns actually add detail
    if ((x1 - x0) / (full[1] - full[0]) < 0.5) setHi(x0, x1); else setPNG();
  }});
  const btn = document.getElementById(label + '-btn');
  if (btn) btn.addEventListener('click', function () {{
    if (mode === 'hi') {{ setPNG(); Plotly.relayout(gd, {{'xaxis.autorange': true}}); }}
    else setHi(full[0], full[1]);
    btn.textContent = mode === 'hi' ? 'back to PNG grid' : 'full resolution';
  }});
}}

// ------------------------------------------------------------------ RAYLEIGH SUCCESS
if (D.kind === 'rayleigh_ok') {{
  const z = Array.from(dec(D.z_km, 'f32'));
  const band = () => (D.fit_alt_km ? [{{ type:'rect', xref:'paper', x0:0, x1:1,
      y0:D.fit_alt_km[0], y1:D.fit_alt_km[1], fillcolor:'rgba(255,193,7,0.18)',
      line:{{width:0}}, layer:'below' }}] : []);
  const profile = (div, series, title, xlab, logx) => {{
    const tr = series.filter(s => D.prof[s.k]).map(s => ({{
      x: Array.from(dec(D.prof[s.k], 'f32')), y: z, mode:'lines', name: s.n,
      line:{{ color:s.c, width:s.w || 1.1, dash:s.d }} }}));
    const xa = {{ title:{{text:xlab}} }};
    if (logx) {{ xa.type = 'log'; xa.exponentformat = 'e'; }}
    Plotly.newPlot(div, tr, Object.assign({{}}, LAY, {{ title:{{text:title,font:{{size:12}}}},
      shapes:band(), xaxis:xa,
      yaxis:{{title:{{text:'Altitude (km ASL)'}}, range:[0, D.z_max_km]}},
      legend:{{orientation:'h', y:1.13, font:{{size:10}}}} }}), CFG);
  }};
  profile('p1', [{{k:'rcs_mean', n:'RCS', c:'#1f77b4', w:0.9}}], 'Molecular: RCS', 'Raw RCS', false);
  profile('p2', [{{k:'signal_normalized', n:'normalised', c:'#1f77b4', w:0.9}},
                 {{k:'p_mol', n:'molecular', c:'#d62728', d:'dash'}}],
          'Molecular: Signal vs Theory', 'Signal (Mm⁻¹)', true);
  profile('p3', [{{k:'beta_att', n:'retrieved β_att', c:'#1f77b4', w:0.9}},
                 {{k:'beta_att_mol', n:'molecular β_att', c:'#d62728', d:'dash'}}],
          'Molecular: Attenuated Backscatter', 'β_att (Mm⁻¹ sr⁻¹)', false);

  const rkm = D.range_km ? Array.from(dec(D.range_km, 'f32')) : null;
  const hkm = D.half_km ? Array.from(dec(D.half_km, 'f32')) : null;
  [['slopes','Window: Slope','RdBu','Slope'],['intercepts','Window: Intercept','Magma','|Intercept|'],
   ['r_squared','Window: R²','YlGnBu','R²']].forEach(g => {{
    const G = D.grids[g[0]]; if (!G) return;
    const a = dec(G.b64,'f32'), [n1,n2] = G.shape, zz = [];
    for (let j = 0; j < n2; j++) {{ const row = new Array(n1);
      for (let i = 0; i < n1; i++) row[i] = a[i*n2+j]; zz.push(row); }}
    const tr = [{{ type:'heatmap', z:zz, x:rkm, y:hkm, colorscale:g[2],
        zmin: g[0]==='r_squared'?0:undefined, zmax: g[0]==='r_squared'?1:undefined,
        colorbar:{{title:{{text:g[3]}}, thickness:12}},
        hovertemplate:'centre %{{x:.2f}} km<br>half %{{y:.2f}} km<br>%{{z:.3g}}<extra></extra>' }}];
    if (D.best_km) tr.push({{ x:[D.best_km[0]], y:[D.best_km[1]], mode:'markers',
        marker:{{symbol:'x', size:11, color:'#d62728', line:{{width:2}}}}, name:'chosen window' }});
    Plotly.newPlot('g_'+g[0], tr, Object.assign({{}}, LAY, {{ title:{{text:g[1],font:{{size:12}}}},
      showlegend:false, xaxis:{{title:{{text:'Centre range (km)'}}}},
      yaxis:{{title:{{text:'Half-length (km)'}}}} }}), CFG);
  }});

  const c = D.curtain, X = xOf(c), rk = Array.from(dec(c.range_km,'f32'));
  const traces = [{{ type:'heatmap', z:grid(c), x:X, y:rk, colorscale:'Viridis',
      zmin:c.lo, zmax:c.hi, colorbar:{{title:{{text:'log₁₀(RCS)'}}, thickness:13}},
      hovertemplate:'%{{x}}<br>%{{y:.2f}} km<br>log₁₀RCS %{{z:.2f}}<extra></extra>' }}];
  const shapes = [];
  const bandRuns = (runs, color, name) => {{
    (runs||[]).forEach(r => shapes.push({{ type:'rect', yref:'paper', y0:0, y1:1,
        x0:X[r[0]], x1:X[Math.min(r[1]+1, X.length-1)], fillcolor:color,
        line:{{width:0}}, layer:'above' }}));
    if ((runs||[]).length) traces.push({{ x:[null], y:[null], mode:'markers', name:name,
        marker:{{size:9, color:color, symbol:'square'}} }});
  }};
  bandRuns(c.excl_lowcloud, 'rgba(214,39,40,0.28)', 'excluded (low cloud)');
  bandRuns(c.excl_screened, 'rgba(70,70,70,0.26)', 'screened / not used');
  if (c.high_cloud_ylo && c.y_hi_km && c.has_high_cloud) {{
    // Two traces, both gap-free: the ceiling at y_hi and the floor at y_lo (== y_hi outside the
    // mask, so the band collapses to nothing there rather than being bridged across a null).
    // shape 'hvh' reproduces matplotlib's step="mid".
    traces.push({{ x:X, y:X.map(() => c.y_hi_km), mode:'lines', line:{{width:0, shape:'hvh'}},
        hoverinfo:'skip', showlegend:false }});
    traces.push({{ x:X, y:c.high_cloud_ylo, mode:'lines', fill:'tonexty',
        fillcolor:'rgba(25,211,243,0.22)', line:{{width:0, shape:'hvh'}},
        name:'high cloud (masked above fit)', hoverinfo:'skip' }});
  }}
  if (c.cbh_km) traces.push({{ x:X, y:c.cbh_km, mode:'markers', name:'cloud base',
      marker:{{size:3.5, color:'white', line:{{color:'black', width:0.4}}}},
      hovertemplate:'%{{x}}<br>cloud base %{{y:.2f}} km<extra></extra>' }});
  if (c.mol_layer_km) {{
    shapes.push({{ type:'rect', xref:'paper', x0:0, x1:1, y0:c.mol_layer_km[0], y1:c.mol_layer_km[1],
        fillcolor:'rgba(255,193,7,0.20)', line:{{color:'#ffc107', width:1.6}}, layer:'above' }});
    traces.push({{ x:[null], y:[null], mode:'lines', line:{{color:'#ffc107',width:2}},
        name:'molecular layer' }});
  }}
  // No descriptive title here: the legend already names every overlay, and the two collided.
  Plotly.newPlot('curtain', traces, Object.assign({{}}, LAY, {{ height:450, shapes:shapes,
    margin:{{l:74, r:12, t:44, b:48}},
    xaxis:{{title:{{text: xIsDate(c)?'Time (UTC)':'Hours since start'}}, type: xIsDate(c)?'date':'linear'}},
    yaxis:{{title:{{text:'Range (km)'}}}},
    legend:{{orientation:'h', y:1.12, x:0, xanchor:'left', font:{{size:9.5}}}} }}), CFG);
  attachHiRes('curtain', c, D.curtain_hi, 'res');

  if (D.sens) {{
    const S = D.sens, a = dec(S.b64,'f32'), [n1,n2] = S.shape, zz = [];
    for (let i = 0; i < n1; i++) {{ const row = new Array(n2);
      for (let j = 0; j < n2; j++) row[j] = a[i*n2+j]; zz.push(row); }}
    Plotly.newPlot('sens', [{{ type:'heatmap', z:zz,
      x:S.shift.map(v => (v>0?'+':'')+v.toFixed(0)), y:S.lr.map(v => v.toFixed(0)),
      colorscale:'RdBu', reversescale:true, zmin:-S.vabs, zmax:S.vabs,
      colorbar:{{title:{{text:'Deviation (%)'}}, thickness:12}},
      hovertemplate:'LR %{{y}} sr<br>shift %{{x}} m<br>%{{z:+.2f}} %<extra></extra>' }}],
      Object.assign({{}}, LAY, {{ title:{{text:'Sensitivity Grid',font:{{size:12}}}},
        xaxis:{{title:{{text:'Alt shift (m)'}}, type:'category'}},
        yaxis:{{title:{{text:'LR (sr)'}}, type:'category'}} }}), CFG);
  }}
  if (D.cl_spread) {{
    const S = D.cl_spread, n = S.vals.length;
    Plotly.newPlot('spread', [{{ x:S.vals.map((_,i)=>1+(n>1?-0.08+0.16*i/(n-1):0)), y:S.vals,
      mode:'markers', name:'LR×shift combos',
      marker:{{size:7, opacity:0.55, color:'#4c78a8'}},
      hovertemplate:'C_L %{{y:.4g}}<extra></extra>' }}],
      Object.assign({{}}, LAY, {{ title:{{text:'Lidar constant C_L spread',font:{{size:12}}}},
        shapes:[{{ type:'rect', xref:'paper', x0:0, x1:1, y0:S.median-S.unc, y1:S.median+S.unc,
                   fillcolor:'rgba(214,39,40,0.12)', line:{{width:0}}, layer:'below' }},
                {{ type:'line', xref:'paper', x0:0, x1:1, y0:S.median, y1:S.median,
                   line:{{color:'#d62728', width:1.4}} }}],
        xaxis:{{visible:false, range:[0.7,1.3]}}, showlegend:false,
        yaxis:{{title:{{text:'Lidar constant C_L'}}, exponentformat:'e'}} }}), CFG);
  }}
}}

// ------------------------------------------------------------------ RAYLEIGH FAILURE
if (D.kind === 'rayleigh_fail') {{
  // the PNG says this ONLY through a dark-red suptitle; here we also surface the flag code, which
  // the image never shows at all
  document.getElementById('banner').innerHTML =
    'NOT CALIBRATED: ' + D.reason +
    (D.flag !== undefined && D.flag !== null
      ? '<span class="flag">flag ' + D.flag + (D.flag_label ? ' — ' + D.flag_label : '') +
        (D.message ? '  ·  message: "' + D.message + '"' : '') + '</span>' : '');
  const c = D.curtain, X = xOf(c), rk = Array.from(dec(c.range_km,'f32'));
  const tr = [{{ type:'heatmap', z:grid(c), x:X, y:rk, colorscale:'Viridis', zmin:c.lo, zmax:c.hi,
      colorbar:{{title:{{text:'log₁₀(RCS)'}}, thickness:13}},
      hovertemplate:'%{{x}}<br>%{{y:.2f}} km<br>log₁₀RCS %{{z:.2f}}<extra></extra>' }}];
  const shapes = [];
  if (c.cbh_km) tr.push({{ x:X, y:c.cbh_km, mode:'markers', name:'cloud base',
      marker:{{size:3.5, color:'white', line:{{color:'black', width:0.4}}}} }});
  if (D.window_range_km) {{
    shapes.push({{ type:'rect', xref:'paper', x0:0, x1:1, y0:D.window_range_km[0],
        y1:D.window_range_km[1], fillcolor:'rgba(255,193,7,0.20)',
        line:{{color:'#ffc107', width:1.4}}, layer:'above' }});
    tr.push({{ x:[null], y:[null], mode:'lines', line:{{color:'#ffc107',width:2}},
        name:'molecular window (attempted)' }});
  }}
  Plotly.newPlot('curtain', tr, Object.assign({{}}, LAY, {{ height:430, shapes:shapes,
    title:{{text:'Range-corrected signal — full night', font:{{size:12}}}},
    xaxis:{{title:{{text: xIsDate(c)?'Time (UTC)':'Hours since start'}}, type: xIsDate(c)?'date':'linear'}},
    yaxis:{{title:{{text:'Range (km)'}}, range:[0, c.y_max_km]}},
    legend:{{orientation:'h', y:1.10, font:{{size:10}}}} }}), CFG);
  attachHiRes('curtain', c, D.curtain_hi, 'res');

  if (D.has_profile) {{
    const z = Array.from(dec(D.z_km,'f32'));
    const t2 = [{{ x:Array.from(dec(D.sig,'f32')), y:z, mode:'lines', name:'signal / range²',
        line:{{color:'#1f77b4', width:0.9}}, connectgaps:true }}];
    if (D.pm_scaled) t2.push({{ x:Array.from(dec(D.pm_scaled,'f32')), y:z, mode:'lines',
        name:'molecular (scaled)', line:{{color:'#d62728', width:1.0, dash:'dash'}} }});
    Plotly.newPlot('failprof', t2, Object.assign({{}}, LAY, {{ height:430,
      title:{{text:'Night-mean profile vs molecular', font:{{size:12}}}},
      shapes: D.window_alt_km ? [{{ type:'rect', xref:'paper', x0:0, x1:1, y0:D.window_alt_km[0],
          y1:D.window_alt_km[1], fillcolor:'rgba(255,193,7,0.15)', line:{{width:0}},
          layer:'below' }}] : [],
      xaxis:{{title:{{text:'Signal (a.u.)'}}, type:'log', exponentformat:'e'}},
      yaxis:{{title:{{text:'Altitude (km ASL)'}}, range:[0, D.z_max_p_km]}},
      legend:{{orientation:'h', y:1.10, font:{{size:10}}}} }}), CFG);
  }} else {{
    document.getElementById('failprof').innerHTML =
      '<p style="text-align:center;padding:150px 0;color:#66707a">no usable profiles</p>';
  }}
}}

// ------------------------------------------------------------------ CLOUD
if (D.kind === 'cloud') {{
  const c = D.curtain, X = D.hrs, rk = Array.from(dec(c.range_km,'f32'));
  const tr = [{{ type:'heatmap', z:grid(c), x:X, y:rk, colorscale:'Viridis', zmin:c.lo, zmax:c.hi,
      colorbar:{{title:{{text:'log₁₀(β)'}}, thickness:13}},
      hovertemplate:'%{{x:.2f}} h<br>%{{y:.2f}} km<br>log₁₀β %{{z:.2f}}<extra></extra>' }}];
  const shapes = [{{ type:'rect', xref:'paper', x0:0, x1:1, y0:D.gate_km[0], y1:D.gate_km[1],
      fillcolor:'rgba(255,255,255,0.10)', line:{{width:0}}, layer:'above' }}];
  (D.sel_runs||[]).forEach(r => shapes.push({{ type:'rect', yref:'paper', y0:0, y1:1,
      x0:X[r[0]], x1:X[Math.min(r[1]+1, X.length-1)], fillcolor:'rgba(44,160,44,0.18)',
      line:{{width:0}}, layer:'above' }}));
  if (D.sel_runs && D.sel_runs.length) tr.push({{ x:[null], y:[null], mode:'markers',
      name:'used for calibration', marker:{{size:9, color:'rgba(44,160,44,0.5)', symbol:'square'}} }});
  if (D.cbh_unsel_km) tr.push({{ x:X, y:D.cbh_unsel_km, mode:'markers', name:'cloud base',
      marker:{{size:3.5, color:'white', line:{{color:'black', width:0.4}}}} }});
  if (D.cbh_sel_km) tr.push({{ x:X, y:D.cbh_sel_km, mode:'markers', name:'cloud base (used)',
      marker:{{size:6, color:'#2ca02c', line:{{color:'black', width:0.5}}}} }});
  Plotly.newPlot('cl_curtain', tr, Object.assign({{}}, LAY, {{ height:420, shapes:shapes,
    title:{{text:'Attenuated backscatter — cloud base (dots), profiles used (green)', font:{{size:12}}}},
    xaxis:{{title:{{text:'Hours since start'}}}},
    yaxis:{{title:{{text:'Range (km AGL)'}}, range:[0, D.y_max_km]}},
    legend:{{orientation:'h', y:1.10, font:{{size:9.5}}}} }}), CFG);
  attachHiRes('cl_curtain', c, D.curtain_hi, 'res');

  const pt = [];
  if (D.prof) pt.push({{ x:Array.from(dec(D.prof,'f32')), y:Array.from(dec(D.rng_km,'f32')),
      mode:'lines', name:'mean of ' + D.n_sel + ' selected', line:{{color:'#2ca02c', width:1.1}} }});
  const pshapes = [{{ type:'rect', xref:'paper', x0:0, x1:1, y0:D.gate_km[0], y1:D.gate_km[1],
      fillcolor:'rgba(255,193,7,0.14)', line:{{width:0}}, layer:'below' }}];
  if (D.b_layer_km) pshapes.push({{ type:'rect', xref:'paper', x0:0, x1:1, y0:D.b_layer_km[0],
      y1:D.b_layer_km[1], fillcolor:'rgba(44,160,44,0.22)', line:{{width:0}}, layer:'below' }});
  if (D.med_cbh_km !== undefined) pshapes.push({{ type:'line', xref:'paper', x0:0, x1:1,
      y0:D.med_cbh_km, y1:D.med_cbh_km, line:{{color:'#d62728', width:1.2, dash:'dash'}} }});
  Plotly.newPlot('cl_prof', pt, Object.assign({{}}, LAY, {{ height:420, shapes:pshapes,
    title:{{text:'Representative profile — gold = 100–2400 m gate, green = B accumulates', font:{{size:11}}}},
    xaxis:{{title:{{text:'β (m⁻¹ sr⁻¹)'}}, type:'log', exponentformat:'e'}},
    yaxis:{{title:{{text:'Range (km AGL)'}}, range:[0, D.y_max_km]}},
    legend:{{orientation:'h', y:1.10, font:{{size:10}}}} }}), CFG);

  Plotly.newPlot('cl_S', [
    {{ x:D.hrs, y:D.S_app, mode:'markers', name:'apparent S',
       marker:{{size:4, color:'#999', opacity:0.5}} }},
    {{ x:D.hrs, y:D.S_con, mode:'markers', name:'consistent S (used)',
       marker:{{size:6, color:'#2ca02c'}} }}],
    Object.assign({{}}, LAY, {{ title:{{text:'Apparent vs consistent lidar ratio', font:{{size:12}}}},
      shapes:[{{ type:'line', xref:'paper', x0:0, x1:1, y0:D.s_theo, y1:D.s_theo,
                 line:{{color:'#d62728', width:1.2, dash:'dash'}} }}],
      xaxis:{{title:{{text:'Hours since start'}}}},
      yaxis:{{title:{{text:'Lidar ratio S (sr)'}}, range:[0, D.s_ymax]}},
      legend:{{orientation:'h', y:1.13, font:{{size:9.5}}}} }}), CFG);

  const hs = [];
  if (D.cal_median !== null) {{
    if (D.cal_std !== null) hs.push({{ type:'rect', yref:'paper', y0:0, y1:1,
        x0:D.cal_median-D.cal_std, x1:D.cal_median+D.cal_std,
        fillcolor:'rgba(214,39,40,0.12)', line:{{width:0}}, layer:'below' }});
    hs.push({{ type:'line', yref:'paper', y0:0, y1:1, x0:D.cal_median, x1:D.cal_median,
        line:{{color:'#d62728', width:1.4}} }});
  }}
  Plotly.newPlot('cl_hist', [{{ type:'histogram', x:D.coeffs, marker:{{color:'#2ca02c'}},
      opacity:0.75, nbinsx:Math.min(30, Math.max(5, Math.floor(D.coeffs.length/2))) }}],
    Object.assign({{}}, LAY, {{ title:{{text:'Calibration coefficient distribution', font:{{size:12}}}},
      shapes:hs, xaxis:{{title:{{text:"O'Connor coefficient C (per profile)"}}}},
      yaxis:{{title:{{text:'count'}}}}, showlegend:false }}), CFG);

  document.getElementById('cl_summary').textContent = (D.summary || []).join('\\n');
}}
</script></body></html>"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--key", required=True)
    ap.add_argument("--date", required=True)
    ap.add_argument("--method", default="rayleigh", choices=["rayleigh", "cloud"])
    ap.add_argument("--l1-root", required=True)
    ap.add_argument("--cams", default="")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=Path("./_dyn_diag_work"))
    ap.add_argument("--skip-ab", action="store_true")
    ap.add_argument("--no-hires", action="store_true",
                    help="omit the opt-in native-resolution curtain")
    args = ap.parse_args()

    global HIRES
    HIRES = not args.no_hires
    os.environ["ALC_L1_ROOT"] = args.l1_root
    os.environ["ALC_FULLCAL_DIR"] = str(args.work.resolve())
    if args.cams:
        os.environ["ALC_CAMS_DIR"] = args.cams
    args.work.mkdir(parents=True, exist_ok=True)

    ab = "not measured (--skip-ab)"
    if not args.skip_ab:
        print("timing the same night with images OFF then ON ...", flush=True)
        t_off = _time_runner(args.key, args.date, "0", args.method)
        t_on = _time_runner(args.key, args.date, "1", args.method)
        ab = (f"{t_off:.1f} s without images vs {t_on:.1f} s with -> images cost "
              f"{t_on - t_off:+.1f} s ({100 * (t_on - t_off) / max(t_off, 1e-9):+.0f} %)")
        print("  " + ab, flush=True)

    os.environ["PLOTS"] = "1"
    _install_capture()
    import importlib
    from datetime import datetime
    rnc = importlib.import_module("scripts.run_network_calibration")
    census = json.loads((REPO / "validation/scope_l1_2026_census.json").read_text(encoding="utf-8"))
    rows = census if isinstance(census, list) else census.get("streams", [])
    s = next(r for r in rows if f"{r['wmo']}_{r['ident']}" == args.key)
    d = datetime.strptime(args.date, "%Y%m%d")
    print(f"calibrating {args.key} {args.date} ({args.method}) in-process ...", flush=True)
    res_rows = (rnc._do_cloud(s, d, d) if args.method == "cloud" else rnc._do_rayleigh(s, d, d))

    if "kind" not in CAPTURED:
        print("\nNo diagnostic was drawn. A cloud PNG is only drawn for a success or a cloud "
              "rejection (flag <= -20); some Rayleigh rejections (-4, -10) draw nothing at all.",
              file=sys.stderr)
        sys.exit(2)

    t0 = time.perf_counter()
    kind = CAPTURED["kind"]
    if kind == "cloud":
        a = CAPTURED["args"]
        payload = _payload_cloud(a[0], a[1], (a[2] if len(a) > 2 else CAPTURED["kwargs"].get("title", "")))
    elif kind == "rayleigh_fail":
        payload = _payload_rayleigh_fail(CAPTURED["kwargs"])
        # the flag never appears on the PNG; attach it from the calibration result
        for r in (res_rows or []):
            if str(r.get("date")) == args.date:
                payload["flag"] = r.get("flag")
                payload["message"] = r.get("message")
                try:
                    from calibration.flags import FLAG_MEANINGS
                    payload["flag_label"] = FLAG_MEANINGS.get(float(r.get("flag")))
                except Exception:  # noqa: BLE001
                    pass
                break
    else:
        payload = _payload_rayleigh_ok(CAPTURED["kwargs"])
    raw = json.dumps(payload, separators=(",", ":"))
    dump_ms = (time.perf_counter() - t0) * 1e3
    gz_kb = len(gzip.compress(raw.encode())) / 1024
    png = CAPTURED.get("png")
    png_kb = (png.stat().st_size / 1024) if png and png.exists() else float("nan")
    render_s = max(CAPTURED.get("renders", [float("nan")]))

    import plotly.offline as pyo
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(HTML.format(plotly=pyo.get_plotlyjs(), payload=raw, key=args.key,
                                    date=args.date, kind=kind, body=BODY[kind],
                                    size_kb=len(raw) / 1024, gz_kb=gz_kb, png_kb=png_kb,
                                    render_s=render_s, dump_ms=dump_ms, ab=ab), encoding="utf-8")
    print(f"\n{'':-<64}\nkind          : {kind}")
    print(f"PNG           : {png_kb:8.0f} KB   rendered in {render_s:6.2f} s")
    print(f"payload       : {len(raw)/1024:8.0f} KB   built in    {dump_ms:6.0f} ms")
    print(f"payload gzip  : {gz_kb:8.0f} KB   ({png_kb / max(gz_kb, 1e-9):.1f}x smaller)")
    print(f"mockup        : {args.out}")


if __name__ == "__main__":
    main()
