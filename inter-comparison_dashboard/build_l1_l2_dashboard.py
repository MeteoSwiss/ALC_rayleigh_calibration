# -*- coding: utf-8 -*-
"""L1-vs-L2 inter-comparison dashboard — compute stage (per SITE, see sites.py).

Builds every combination of the three optional corrections for every available source and dumps one
JSON the static HTML page then toggles through client-side. Nothing is recomputed in the browser.

  LEFT  panel  "L1 + calibration"     : L1 rcs_0 / C_L(t) * 1e6, C_L = the daily Kalman series of
                                        the site's calibration variants (l1_l2_calib.py)
  RIGHT panel  "L2 as distributed"    : attenuated_backscatter_0 straight out of the L2 product —
                                        only for sites whose spec lists an "L2" source (Payerne)

Optional corrections, applied IDENTICALLY to both panels so the panels stay comparable:
  * water vapour        - divide by the two-way WV transmission at the laser line (910 nm units only)
  * wavelength, simple  - single Angstrom exponent (alpha=1) to the 1064 nm reference
  * wavelength, advanced- component-separated molecular/aerosol conversion using the CAMS T/p profile

All panels are reduced to the SAME hourly grid and the SAME altitude grid, then restricted to the
hours where every instrument has data in every source, so every curve and every difference on the
page comes from one strictly paired sample. The canonical combo (first variant, WV on, molecular
wavelength) additionally emits a coarsened time-height "pcolor" block for the per-channel curtains.

Run:  python inter-comparison_dashboard/build_l1_l2_dashboard.py [payerne|amsterdam|lindenberg] [--reread]
Out:  C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard/data_<site>.json
"""
from __future__ import annotations
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")                # validation.paper, monitoring
sys.path.insert(0, str(Path(__file__).resolve().parent))   # sibling modules (folder name is not a valid package name)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")
# The 0.4 deg CAMS monthlies (L137 model levels) — finer than the 1 deg archive, which matters at
# Payerne where the 1 deg grid-cell surface sits ~900 m above the station.
os.environ.setdefault("ALC_VAL_CAMS_04", "A:/CAMS_Monthly_04")

import numpy as np

from validation.paper import intercompare as IC
from validation.paper import run_paper_validation as RP
from validation.paper.calib_benchmark import key_of
import l1_l2_calib as CAL
import l1_l2_io as IO
import sites

TARGET, ALPHA = 1064.0, 1.0
ZMIN, ZMAX = RP.ZMIN, RP.ZMAX               # 500-3000 m AGL statistics band
# Display grid. The instruments reach 15.3 km (CHM15k), 15.7 km (CL61) and 7.7 km (CL31), so the
# grid runs to 15 km and each curve is truncated at its own noise floor instead of at a fixed height.
ZTOP_PLOT = 15000.0                         # display range [m AGL]
DZ = 30.0                                   # display altitude step [m] (the L2 native grid)
# Pcolor payload coarsening: only up to 8 km, every 2nd display gate, one canonical combo.
PCOLOR_ZTOP = 8000.0
PCOLOR_ZSTEP = 2

OUT = sites.DATA_ROOT                       # data + rendered page live outside the repo

WL_MODES = ["none", "angstrom", "molecular"]

# ---------------------------------------------------------------------------- per-site state
# set_site() fills these from sites.py; the pool initializer re-applies it in every worker (the
# workers re-import this module, which defaults to Payerne).
SITE = None
WMO = None
STA = None
CHANNELS = None
IREF = 0
SOURCES = ("L1", "L2")
START = END = WIN_START = WIN_END = None
CALIB_VARIANTS = {}
COMBOS = []
DARK_NPZ = None
CACHE_NPZ = None
_DARKP = {}


def set_site(site_key):
    """Point every module-level knob at ONE site from sites.py."""
    global SITE, WMO, STA, CHANNELS, IREF, SOURCES, START, END, WIN_START, WIN_END
    global CALIB_VARIANTS, COMBOS, DARK_NPZ, CACHE_NPZ
    s = sites.get_site(site_key)
    SITE = s
    WMO = s["wmo"]
    STA = (s["lat"], s["lon"], s["alt"])
    CHANNELS = s["channels"]
    IREF = s["iref"]
    SOURCES = tuple(s["sources"])
    START, END = s["start"], s["end"]
    WIN_START, WIN_END = s["win_start"], s["win_end"]
    # Calibration variant is a third dimension: it changes the CONSTANT the L1 panel divides by, so
    # it cannot be applied after the hourly medians are taken (median(beta/c) != median(beta)/
    # median(c)). The L2 panel is independent of it and is recomputed identically for each variant
    # -- wasteful, but it keeps one combo key for the whole payload and the pool absorbs it.
    CALIB_VARIANTS = {v: Path(s["calib_dirs"][v]) for v in s["variants"]}
    if s.get("collapse_combos"):
        # Every channel at 1064 nm: WV and both wavelength conversions are exact no-ops, so one
        # combo per variant carries the whole page (the render stage hides those controls).
        COMBOS = [(v, True, "molecular") for v in s["variants"]]
    else:
        # A site may restrict the wavelength modes (sites.py wl_modes) to contain the page size --
        # Payerne drops the single-Angstrom mode: 4 variants x 2 WV x 2 modes = 16 combos.
        wl_list = [wl for wl in s.get("wl_modes", WL_MODES) if wl in WL_MODES]
        COMBOS = [(v, wv, wl) for v in s["variants"] for wv in (False, True) for wl in wl_list]
    DARK_NPZ = Path(os.environ.get("ALC_DASH_DARK_NPZ", s["dark_npz"])) \
        if s.get("dark_npz") else None
    CACHE_NPZ = OUT / f"_streams_{site_key}.npz"
    _DARKP.clear()


# Measured dark baseline b(z), rcs_0 units, from the covered-telescope campaign. A "dark" variant
# (sites.py dark_variants) subtracts it from the L1 profiles AND uses constants from the
# dark-corrected calibration run -- both sides or neither: a corrected profile divided by an
# uncorrected constant (or vice versa) mixes two signal definitions and WORSENS the comparison
# (measured 2026-08-15). Subtracting after the hourly median is exact because the baseline is
# constant in time.
def dark_for(ident, rng):
    """b(z) resampled on this channel's range grid (rcs_0 units), or None."""
    key = str(ident)
    if key not in _DARKP:
        prof = None
        if DARK_NPZ is not None and DARK_NPZ.exists():
            try:
                with np.load(DARK_NPZ) as z:
                    if f"{ident}_b_rcs" in z:
                        prof = (np.asarray(z[f"{ident}_range"], "f8"),
                                np.asarray(z[f"{ident}_b_rcs"], "f8"))
            except Exception:
                prof = None
        _DARKP[key] = prof
    prof = _DARKP[key]
    if prof is None:
        return None
    rd, bd = prof
    ok = np.isfinite(rd) & np.isfinite(bd)
    b = np.interp(np.asarray(rng, "f8"), rd[ok], bd[ok], left=np.nan, right=np.nan)
    return np.nan_to_num(b, nan=0.0)


set_site("payerne")                          # default; main()/_init() re-point per run


def combo_key(cal, wv, wl):
    return f"{cal}_wv{int(wv)}_{wl}"


def range_key(i, j):
    """Key for the contiguous month range months[i..j] (inclusive)."""
    return f"r{i}_{j}"


# --------------------------------------------------------------------------- corrections
def apply_corrections(beta, d, itype, wv, wl):
    """The three optional corrections, in the pipeline's order (WV first, then wavelength).
    Both are no-ops for a 1064 nm unit; WV is a no-op outside the 910 nm absorption band.
    WV and the advanced wavelength conversion go through the DAY-granular CAMS path (l1_l2_io) --
    Jun-Aug 2026 has only the daily CAMS cache, which a monthly resolver would miss entirely."""
    lam = float(d["wavelength"]) if np.isfinite(d["wavelength"]) else \
        IC.WV_PARAMS.get(itype, (np.nan,))[0]
    if wv and np.isfinite(lam) and IC.in_water_vapor_band(lam):
        lam0, fwhm = IC.WV_PARAMS.get(itype, (910.0, 3.4))
        beta, _ = IO.apply_wv(beta, d, lam0, fwhm)
    lam_use = lam if np.isfinite(lam) else TARGET
    if abs(lam_use - TARGET) >= 1.0:
        if wl == "molecular":
            beta = IO.wavelength_correct_molecular(beta, d, lam_use, TARGET, ALPHA)
        elif wl == "angstrom":
            beta = IC.wavelength_correct(beta, d, lam_use, TARGET, ALPHA, "angstrom")
    return beta


def process_channel(d, ch, source, wv, wl, cal="v2.0"):
    """One channel, one source, one correction combo -> the dict grid_and_stats() consumes."""
    beta = d["beta"].copy()
    med_corr = 1.0
    if source == "L1" and cal in SITE.get("dark_variants", ()):
        b = dark_for(ch["ident"], np.asarray(d["alt"], "f8") - float(d["station_alt"]))
        if b is not None:
            beta = beta - b[None, :]
    if source == "L1":
        saved = IC.CALIB
        try:
            IC.CALIB = CALIB_VARIANTS[cal]
            ser = IC.load_calib_series(key_of(dict(wmo=WMO, ident=ch["ident"], calib=ch["calib"])), "L1")
        finally:
            IC.CALIB = saved
        if ser is None:
            raise RuntimeError(f"no {cal} calibration series for {ch['label']}")
        ck = IC.interp_calib(ser[0], ser[1], d["time"])
        beta = beta / ck[:, None] * 1e6
        med_corr = float(np.nanmedian(1e6 / ck))
    beta = apply_corrections(beta, d, ch["itype"], wv, wl)
    disp = beta.copy()
    disp[d["qf"] > 0] = np.nan
    return dict(ch=dict(ch, wmo=WMO, lat=STA[0], lon=STA[1], alt=STA[2]), l1=d,
                beta_scr=IC.screen(beta, d), beta_disp=disp, beta_scr_noovl=None,
                med_corr=med_corr)


# --------------------------------------------------------------------------- one panel
def _grid(items, want_disp=False):
    """Hourly-median grid + common altitude grid for one panel.

    This is run_paper_validation.grid_and_stats trimmed to what the dashboard needs, with ONE
    deliberate difference: `snr_idx=()` instead of its hard-coded `(0,)`, i.e. the per-gate SNR>=3
    detection gate is NOT applied.

    That gate would make the two panels incomparable. It is evaluated inside each hourly window and
    needs >=5 samples: L2 arrives at native 5-minute resolution (~12 samples/hour, gate active)
    while L1 has already been hourly-median retimed by the reader (1 sample/hour, gate inert). An
    SNR gate keeps only gates whose median is significantly non-zero, so where it IS active it
    conditions on positive noise and biases the retained median high and rising with altitude —
    which is precisely the CL31 L1-vs-L2 floor discrepancy it produced. Leaving it off for both is
    also the project's validated choice for a bias comparison (noise filter 'none').

    want_disp additionally grids the qf-masked-only stream (clouds kept) — only the canonical
    pcolor combo pays for it."""
    gridded = []
    for it in items:
        l1 = it["l1"]
        cbh = l1["cbh"]
        cbh_low = (np.nanmin(np.where((cbh > 0) & (cbh < 20000), cbh, np.nan), axis=1)
                   if cbh.ndim == 2 else cbh)
        arrays = [it["beta_scr"], cbh_low] + ([it["beta_disp"]] if want_disp else [])
        g, arrs = IC.retime_hourly(l1["time"], arrays, min_cov_s=IC.MIN_AVG_S, snr_idx=())
        gridded.append(dict(l1=l1, grid=g, scr=arrs[0],
                            disp=(arrs[2] if want_disp else None)))

    union = np.unique(np.concatenate([g["grid"] for g in gridded]))
    for g in gridded:
        idx = {t: i for i, t in enumerate(g["grid"])}
        pos = np.array([idx.get(t, -1) for t in union])
        g["scrU"] = RP._reindex(g["scr"], pos, union.size)
        g["dispU"] = RP._reindex(g["disp"], pos, union.size) if want_disp else None
    # Deliberately NOT intercompare.build_common_grid: it takes the INTERSECTION of the channels'
    # vertical extents (z1 = min of the maxima), so CL31's 7.7 km ceiling would truncate CHM15k
    # (15.3 km) and CL61 (15.7 km) as well. Each channel is regridded from its OWN native altitude
    # straight onto the fixed display grid instead, so every instrument keeps its full range and the
    # short one simply carries NaN above its own top.
    return union, [(g["scrU"], g["dispU"], g["l1"]["alt"]) for g in gridded]


def build_panel(cache, source, wv, wl, cal="v2.0", want_disp=False):
    """One source + one correction combo, regridded onto the shared display grid.
    Returns (time_sync, beta[nch] on the display grid, alt_asl, disp[nch] or None)."""
    items = [process_channel(cache[(source, ch["ident"])], ch, source, wv, wl, cal)
             for ch in CHANNELS]
    tsync, streams = _grid(items, want_disp)
    alt = np.arange(0.0, ZTOP_PLOT + DZ, DZ) + STA[2]          # ASL display grid
    B = [IC.regrid(b, a, alt) for b, _, a in streams]
    Bd = [IC.regrid(d_, a, alt) for _, d_, a in streams] if want_disp else None
    return tsync, B, alt, Bd


# --------------------------------------------------------------------------- statistics helpers
def _median_iqr(B):
    """Column-wise median / 25th / 75th percentile of a (time x alt) block, NaN-safe."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return (np.nanmedian(B, axis=0), np.nanpercentile(B, 25, axis=0),
                np.nanpercentile(B, 75, axis=0), np.sum(np.isfinite(B), axis=0))


def _keep_mask(med, n, z_agl, nprof, minfrac=0.20, rise=1.5):
    """Where the median profile is still a measurement rather than a noise/offset floor.

    Four cuts, all applied above ZMIN and all truncating everything above the first offending gate:

    * coverage < 20 % of the sample (the rule validation.paper.figures._truncate_noise_floor uses);
    * a non-positive median -- the channel has run out of signal entirely;
    * a median that has turned back UP by more than `rise` x its running minimum AND never comes
      back down. An elevated aerosol layer also lifts the median, but the profile drops again above
      it; an r^2-amplified electronic offset only ever grows, so requiring the rise to persist to
      the top separates the two;
    KNOWN LIMIT: this catches a floor that RISES, not one that goes flat. CL31 sits flat around
    0.33-0.40 from ~1.5 km up and only crosses zero near 4.9 km, so its curve is drawn well past the
    ~1.5-2 km where it stops measuring backscatter. Both panels now behave identically (the earlier
    L1-vs-L2 asymmetry was the SNR gate, see _grid), so the filter is at least self-consistent, but a
    principled ceiling needs the detection-limit profile from calibration/sensitivity rather than a
    shape heuristic -- tried and rejected here: a "stalled decay" test never fires (the floor drifts
    down through zero, so new running minima keep appearing) and an "at least molecular decay" test
    cuts at 510 m (any fit window low enough to start in real signal already contains the floor).

    The mask is reported alongside the untruncated curves, never baked into them, so the page can
    switch the filter off and show the raw median.
    """
    good = np.isfinite(med) & (n >= max(10, minfrac * nprof))
    keep = np.ones(med.size, bool)

    bad = np.where((~good | (med <= 0)) & (z_agl > ZMIN))[0]
    cut = bad[0] if bad.size else med.size

    idx = np.array([i for i in np.where(z_agl > ZMIN)[0] if i < cut and good[i] and med[i] > 0])
    if idx.size:
        v = med[idx]
        # sustained rise: scan down from the top for the lowest gate from which every gate stays over
        sustained = np.logical_and.accumulate((v > rise * np.minimum.accumulate(v))[::-1])[::-1]
        w = np.where(sustained)[0]
        if w.size:
            cut = min(cut, int(idx[w[0]]))
    keep[cut:] = False
    return keep & good


# Relative-difference histogram: fixed edges so every combo/range/source shares one axis, and the
# same quantity the stats table reports the median of.
HIST_EDGES = np.arange(-100.0, 200.0 + 2.5, 2.5)


def _rel_diff(cur, ref):
    """Per-sample relative difference (cur-ref)/ref in %, on intercompare._stats' convention:
    both finite and ref > 0, but cur is NOT required positive. Keeping the noise-floor negatives is
    what makes the histogram median equal the `medrelbias_pct` in the statistics table."""
    m = np.isfinite(cur) & np.isfinite(ref) & (ref > 0)
    return np.where(m, (cur - ref) / np.where(m, ref, 1.0) * 100.0, np.nan)


def _rel_hist(cur, ref, band):
    """Distribution of the per-sample relative difference (cur vs ref) in %, over the stats band."""
    v = _rel_diff(cur[:, band], ref[:, band]).ravel()
    v = v[np.isfinite(v)]
    if not v.size:
        return [0] * (HIST_EDGES.size - 1), None
    h, _ = np.histogram(v, HIST_EDGES)
    return [int(x) for x in h], round(float(np.median(v)), 3)


def _rel_profile(cur, ref):
    """Median per-hour relative difference profile (cur vs ref) in %, NaN-safe.
    Ratios, not a ratio of medians, so a few bright hours cannot dominate. Same sample convention
    as _rel_diff / intercompare._stats."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmedian(_rel_diff(cur, ref), axis=0)


def _clean(a):
    """numpy -> JSON, NaN as null."""
    return [None if not np.isfinite(x) else round(float(x), 6) for x in np.asarray(a, float)]


# --------------------------------------------------------------------------- worker
_CACHE = None


def _init(site_key, npz_path, t0, t1):
    global _CACHE
    set_site(site_key)
    _CACHE = load_cache(npz_path, t0, t1)


def _stats_for(B, sel, band, z_agl):
    """Every profile / difference / statistic for ONE subset of the paired hours.

    Curves are emitted UNTRUNCATED, with a per-gate `keep` mask alongside: the page applies the
    noise-floor filter client-side so it can be switched off and the raw median inspected."""
    out = {}
    nh = int(sel.sum())
    Bs = {s: [b[sel] for b in B[s]] for s in SOURCES}
    keep = {}
    for s in SOURCES:
        prof, stats, keep[s] = [], [], []
        for k in range(len(CHANNELS)):
            med, q1, q3, n = _median_iqr(Bs[s][k])
            kp = _keep_mask(med, n, z_agl, nh)
            keep[s].append(kp)
            prof.append(dict(med=_clean(med), q1=_clean(q1), q3=_clean(q3), n=_clean(n),
                             keep=[int(x) for x in kp],
                             ztop=(float(z_agl[kp].max()) if kp.any() else None)))
            stats.append({kk: (None if not np.isfinite(vv) else round(float(vv), 4))
                          for kk, vv in IC._stats(Bs[s][k], Bs[s][IREF], band).items()})
        out[s] = dict(profiles=prof, stats=stats)
    # a difference is only meaningful where BOTH curves in it are still measurements
    for s in SOURCES:
        out[s]["diff_vs_ref"] = [_clean(_rel_profile(Bs[s][k], Bs[s][IREF]))
                                 for k in range(len(CHANNELS))]
        out[s]["diff_keep"] = [[int(x) for x in (keep[s][k] & keep[s][IREF])]
                               for k in range(len(CHANNELS))]
        hists = [_rel_hist(Bs[s][k], Bs[s][IREF], band) for k in range(len(CHANNELS))]
        out[s]["hist"] = [h for h, _ in hists]
        out[s]["hist_med"] = [m for _, m in hists]
    if "L2" in SOURCES:
        out["l1_vs_l2"] = [_clean(_rel_profile(Bs["L1"][k], Bs["L2"][k]))
                           for k in range(len(CHANNELS))]
        out["l1_vs_l2_keep"] = [[int(x) for x in (keep["L1"][k] & keep["L2"][k])]
                                for k in range(len(CHANNELS))]
    out["n_hours"] = nh
    return out


# --------------------------------------------------------------------------- pcolor payload
def _pcolor_z(M):
    """One coarse (time x alt) block -> JSON rows PER ALTITUDE (plotly heatmap z[y][x]):
    log10 of the value, 2 decimals, null where empty/non-positive."""
    with np.errstate(all="ignore"):
        L = np.log10(M)
    L[~np.isfinite(L)] = np.nan
    return [[None if not np.isfinite(v) else round(float(v), 2) for v in L[:, iz]]
            for iz in range(L.shape[1])]


def _pcolor_payload(hours, z_agl, scr, disp):
    """Time-height block for the per-channel curtains — the CANONICAL combo only, coarsened so the
    page stays small: every PCOLOR_ZSTEP-th display gate up to PCOLOR_ZTOP, hourly columns on a
    CONTINUOUS axis (missing hours are null columns, so plotly draws uniform cells), log10 values
    rounded to 2 decimals, and everything above the channel's own keep-mask ceiling nulled.
    `disp` (the qf-masked-only stream, clouds kept) is emitted ONLY where the screened stream is
    empty — it is drawn as a grey backdrop underneath, so anywhere else it would be invisible."""
    zsel = np.where(z_agl <= PCOLOR_ZTOP)[0][::PCOLOR_ZSTEP]
    h = hours.astype("datetime64[h]")
    axis = np.arange(h.min(), h.max() + np.timedelta64(1, "h"), np.timedelta64(1, "h"))
    pos = np.searchsorted(axis, h)
    ch_out = []
    for k in range(len(CHANNELS)):
        S = np.full((axis.size, zsel.size), np.nan)
        S[pos] = scr[k][:, zsel]
        Dq = np.full((axis.size, zsel.size), np.nan)
        Dq[pos] = disp[k][:, zsel]
        # keep-mask ceiling from the full-window screened stream (same rule as the profile curves).
        # nprof = the hours THIS channel actually has data, not the axis length: the unpaired axis
        # includes every hour any instrument reports, and 20 % of that would wipe out a channel
        # whose clear-sky availability is low (Lindenberg CL61: ~140 clear hours on a 721 h axis).
        med, _, _, n = _median_iqr(scr[k])
        nprof = int(np.any(np.isfinite(scr[k]), axis=1).sum())
        kp = _keep_mask(med, n, z_agl, max(nprof, 1))
        S[:, ~kp[zsel]] = np.nan
        Dq[:, ~kp[zsel]] = np.nan
        Dq[np.isfinite(S)] = np.nan
        ch_out.append(dict(scr=_pcolor_z(S), disp=_pcolor_z(Dq)))
    return dict(hours=[str(t) for t in axis],
                alt=[float(z) for z in z_agl[zsel]], ch=ch_out)


def _run(spec):
    cal, wv, wl = spec
    # canonical combo = the page's default view; it alone carries the pcolor block
    canonical = (cal == SITE["variants"][0] and wv and wl == "molecular")
    out = {}
    panels = {}
    for source in SOURCES:
        panels[source] = build_panel(_CACHE, source, wv, wl, cal,
                                     want_disp=(canonical and source == "L1"))
    # --- pair the panels on the hours every source shares -------------------------------------
    common_t = panels[SOURCES[0]][0]
    for s in SOURCES[1:]:
        common_t = np.intersect1d(common_t, panels[s][0])
    idx = {s: np.searchsorted(panels[s][0], common_t) for s in SOURCES}
    alt = panels["L1"][2]
    band = (alt - STA[2] >= ZMIN) & (alt - STA[2] <= ZMAX)
    B = {s: [b[idx[s]] for b in panels[s][1]] for s in SOURCES}
    # hours where every instrument has data in EVERY source -> one strictly paired sample
    have = np.logical_and.reduce([np.any(np.isfinite(B[s][k][:, band]), axis=1)
                                  for s in SOURCES for k in range(len(CHANNELS))])
    for s in SOURCES:
        B[s] = [b[have] for b in B[s]]
    hours = common_t[have]
    z_agl = alt - STA[2]

    # Pcolor from the UNPAIRED L1 panel: the strict pairing above drops every hour ANY instrument
    # is cloud-screened, which would empty the grey backdrop and riddle the curtain with holes --
    # as an overview the full per-channel hours are the informative view.
    if canonical and panels["L1"][0].size:
        out["pcolor"] = _pcolor_payload(panels["L1"][0], z_agl, panels["L1"][1], panels["L1"][3])

    # Per-hour residual of the PWV test channel vs the reference (sites.py `pwv`): the median
    # relative difference over the stats band, one scalar per paired hour, for the four
    # coherent/crossed (constants-variant x WV-comparison) states at the canonical wavelength mode.
    # main() joins these to the per-day CAMS PWV -- the discriminating scatter of the page.
    pwv_spec = SITE.get("pwv")
    if pwv_spec and cal in pwv_spec["cals"] and wl == pwv_spec.get("wl", "molecular"):
        k_c = next((k for k, c in enumerate(CHANNELS)
                    if c.get("chan", c["ident"]) == pwv_spec["chan"]), None)
        if k_c is not None and k_c != IREF and hours.size:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                resid = np.nanmedian(
                    _rel_diff(B["L1"][k_c][:, band], B["L1"][IREF][:, band]), axis=1)
            out["pwv_resid"] = dict(
                hours=[str(t) for t in hours.astype("datetime64[h]")],
                resid=[None if not np.isfinite(v) else round(float(v), 3) for v in resid])

    # Every contiguous run of whole months in the window, so the page can offer a from/to selector.
    # The corrections + gridding above are independent of the subset, so a range costs only medians.
    months = [str(m) for m in np.unique(hours.astype("datetime64[M]"))]
    hm = hours.astype("datetime64[M]").astype(str)
    for i in range(len(months)):
        for j in range(i, len(months)):
            sel = np.isin(hm, months[i:j + 1])
            out[range_key(i, j)] = _stats_for(B, sel, band, z_agl)
    out["months"] = months
    return combo_key(cal, wv, wl), out


# --------------------------------------------------------------------------- cache
_SCALAR = ("station_alt", "lat", "lon", "wavelength", "itype", "wmo", "ident")


def _streams():
    """(source, ident) pairs the cache holds — sources from the site spec, idents de-duplicated
    (the CL61 cloud/Rayleigh channels share one physical stream)."""
    return [(s, i) for s in SOURCES for i in dict.fromkeys(ch["ident"] for ch in CHANNELS)]


def build_stream_cache(path=None):
    """Read all the site's streams ONCE from the daily archives (with the 5-minute fallback) and
    cache them, so re-running the correction combos never touches the archive again."""
    path = CACHE_NPZ if path is None else path
    flat = {}
    for source, ident in _streams():
        rd = IO.read_l1 if source == "L1" else IO.read_l2
        d = rd(WMO, ident, START, END)
        if d is None:
            raise RuntimeError(f"no {source} data for {ident} over {START}..{END}")
        print("  %-3s %s: %-14s %s..%s  days %s" % (
            source, ident, str(d["beta"].shape), str(np.datetime64(d["time"].min(), "D")),
            str(np.datetime64(d["time"].max(), "D")), d["sources"]), flush=True)
        for k, v in d.items():
            if k == "sources" or v is None:
                continue
            flat[f"{source}_{ident}###{k}"] = np.asarray(v)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **flat)
    print(f"cached -> {path.name} ({path.stat().st_size/1e6:.0f} MB)", flush=True)


def _read_cache(path):
    d = np.load(path, allow_pickle=True)
    cache = {}
    for key in d.files:
        stream, k = key.split("###")
        cache.setdefault(stream, {})[k] = d[key]
    for st in cache.values():
        for k, v in list(st.items()):
            if k in _SCALAR and getattr(v, "ndim", None) == 0:
                st[k] = v.item()
    return cache


def load_cache(npz_path, t0, t1):
    """The site's streams, clipped to [t0, t1] so every panel covers the same period."""
    raw = _read_cache(npz_path)
    cache = {}
    for ch in CHANNELS:
        for source in SOURCES:
            d = dict(raw[f"{source}_{ch['ident']}"])
            t = np.asarray(d["time"])
            m = (t >= t0) & (t <= t1)
            for f in ("time", "beta", "qf", "cbh", "vv", "temp_int", "calc"):
                if f in d and d[f] is not None and np.ndim(d[f]) >= 1 and len(d[f]) == m.size:
                    d[f] = np.asarray(d[f])[m]
            d["itype"] = ch["itype"]
            cache[(source, ch["ident"])] = d
    return cache


def common_window(npz_path):
    """The period every source covers for all instruments, clamped to [WIN_START, WIN_END]."""
    raw = _read_cache(npz_path)
    t0 = max(np.asarray(raw[f"{s}_{i}"]["time"]).min() for s, i in _streams())
    t1 = min(np.asarray(raw[f"{s}_{i}"]["time"]).max() for s, i in _streams())
    t0 = max(t0, np.datetime64(datetime.strptime(WIN_START, "%Y%m%d")))
    t1 = min(t1, np.datetime64(datetime.strptime(WIN_END, "%Y%m%d")) + np.timedelta64(1, "D"))
    return t0, t1


# --------------------------------------------------------------------------- L2 applied constant
def l2_applied_constants(npz_path, t0, t1):
    """Daily median of calibration_constant_0 from the L2 product — the constant the RIGHT panel
    is effectively calibrated with, for the time-series charts. Empty for the L1-only sites."""
    if "L2" not in SOURCES:
        return {}
    raw = _read_cache(npz_path)
    out = {}
    for ch in CHANNELS:
        d = raw["L2_" + ch["ident"]]
        c = np.asarray(d.get("calc"), float) if d.get("calc") is not None else None
        t = np.asarray(d["time"])
        if c is None or c.size != t.size:
            continue
        m = np.isfinite(c) & (c > 0) & (t >= t0) & (t <= t1)
        if not m.any():
            continue
        day = t[m].astype("datetime64[D]")
        uniq = np.unique(day)
        vals = [float(np.nanmedian(c[m][day == u])) for u in uniq]
        out[ch["ident"]] = dict(date=[str(u) for u in uniq], value=[round(v, 6) for v in vals])
    return out


# --------------------------------------------------------------------------- PWV (CAMS)
_PWV_DAY = {}
# kg per H2O molecule; a column integral of the number density then lands in kg m^-2, and
# 1 kg m^-2 of liquid water is exactly 1 mm of precipitable water.
_M_H2O = 18.015e-3 / 6.02214076e23


def pwv_mm_for_day(ds):
    """Precipitable water vapour [mm] above the station for one day ('YYYYMMDD'), from the SAME
    CAMS file and the SAME day window the WV correction itself uses (l1_l2_io.cams_for_day +
    intercompare's memoised profile), so the scatter's x-axis is exactly the water the correction
    saw. NaN when no CAMS covers the day (those hours are dropped from the fit)."""
    if ds in _PWV_DAY:
        return _PWV_DAY[ds]
    val = np.nan
    cams = IO.cams_for_day(ds)
    if cams is not None:
        t0, t1 = IO._win(ds)
        prof = IC._cams_wv_profile_cached(str(cams), STA[0], STA[1], t0, t1)
        if prof is not None:
            h = np.asarray(prof[0], "f8")
            n = np.asarray(prof[1], "f8")
            ok = np.isfinite(h) & np.isfinite(n)
            h, n = h[ok], n[ok]
            if h.size >= 2:
                # integrate the column ABOVE the station (below-surface levels are constant fill)
                if h[0] < STA[2] < h[-1]:
                    n0 = np.interp(STA[2], h, n)
                    m = h > STA[2]
                    h = np.concatenate([[STA[2]], h[m]])
                    n = np.concatenate([[n0], n[m]])
                val = float(np.trapz(n, h) * _M_H2O)
    _PWV_DAY[ds] = val
    return val


def assemble_pwv(parts, spec):
    """Join the per-state hourly residual series (from _run) to the per-day CAMS PWV.

    parts: {state_key ('cal_wv0'/'cal_wv1'): {hours, resid}}. States share the hour axis only where
    their pairing agrees (WV-off keeps days the WV exclusion drops), so the union axis carries None
    where a state lacks the hour. The per-state linear fit (slope in % per mm) is computed here
    once, not client-side."""
    hours = sorted({h for p in parts.values() for h in p["hours"]})
    x = np.array([pwv_mm_for_day(h[:10].replace("-", "")) for h in hours])
    resid_by_state, fits = {}, {}
    for st, p in sorted(parts.items()):
        idx = {h: i for i, h in enumerate(p["hours"])}
        r = [p["resid"][idx[h]] if h in idx else None for h in hours]
        resid_by_state[st] = r
        y = np.array([np.nan if v is None else float(v) for v in r])
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() >= 3:
            slope, icept = np.polyfit(x[m], y[m], 1)
            fits[st] = dict(slope=round(float(slope), 3), intercept=round(float(icept), 3),
                            n=int(m.sum()))
            print(f"   pwv {st:<16s} n={int(m.sum()):4d}  slope {slope:+.3f} %/mm  "
                  f"intercept {icept:+.2f} %", flush=True)
    return dict(hours=hours,
                pwv_mm=[None if not np.isfinite(v) else round(float(v), 3) for v in x],
                resid_by_state=resid_by_state, fits=fits,
                chan=spec["chan"], ref=CHANNELS[IREF]["label"],
                band=[float(ZMIN), float(ZMAX)])


def main(site_key=None):
    if site_key is None:
        args = [a for a in sys.argv[1:] if not a.startswith("--")]
        site_key = args[0] if args else "payerne"
    set_site(site_key)
    OUT.mkdir(parents=True, exist_ok=True)
    # Adopt the pre-generalisation Payerne cache under its per-site name (one rename, no re-read).
    legacy = OUT / "_streams.npz"
    if site_key == "payerne" and not CACHE_NPZ.exists() and legacy.exists():
        os.replace(legacy, CACHE_NPZ)
        print(f"adopted legacy cache {legacy.name} -> {CACHE_NPZ.name}", flush=True)
    if not CACHE_NPZ.exists() or "--reread" in sys.argv[1:]:
        print(f"reading {'+'.join(SOURCES)} daily archives {START}..{END} ...", flush=True)
        build_stream_cache()
    npz = str(CACHE_NPZ)
    t0, t1 = common_window(npz)
    print(f"common window: {t0} .. {t1}", flush=True)

    # {variant: {chan: record}} — Payerne keeps its bespoke builder (NetCDFs + study outputs),
    # the other sites smooth the network-runner CSVs (sites.py run_dirs).
    if SITE.get("calib_builder") == "network":
        cal_series = CAL.main_site(SITE)
    else:
        cal_series = CAL.main()
    l2c = l2_applied_constants(npz, t0, t1)

    workers = min(len(COMBOS), max(1, (os.cpu_count() or 8) - 2), 6)
    print(f"computing {len(COMBOS)} correction combos x {len(SOURCES)} sources over "
          f"{workers} workers ...", flush=True)
    combos, months, pcolor, pwv_parts = {}, [], None, {}
    others = [k for k in range(len(CHANNELS)) if k != IREF]
    with ProcessPoolExecutor(max_workers=workers, initializer=_init,
                             initargs=(site_key, npz, t0, t1)) as ex:
        for key, res in ex.map(_run, COMBOS):
            if "pcolor" in res:
                pcolor = res.pop("pcolor")
            if "pwv_resid" in res:
                # state key = variant + WV toggle; the wavelength mode is fixed by the pwv spec
                pwv_parts[key.rsplit("_", 1)[0]] = res.pop("pwv_resid")
            combos[key] = res
            months = res["months"]
            full = res[range_key(0, len(months) - 1)]           # the whole window
            st = full["L1"]["stats"]
            print("   %-16s N=%4d h   L1 vs %s: %s" % (
                key, full["n_hours"], CHANNELS[IREF]["label"].split()[0],
                "  ".join("%s %s" % (CHANNELS[k]["label"].split()[0],
                                     ("%+.1f%%" % st[k]["medrelbias_pct"])
                                     if st[k]["medrelbias_pct"] is not None else "—")
                          for k in others)), flush=True)

    pwv = assemble_pwv(pwv_parts, SITE["pwv"]) if pwv_parts else None

    alt_agl = (np.arange(0.0, ZTOP_PLOT + DZ, DZ)).tolist()
    payload = dict(
        meta=dict(site=site_key, wmo=WMO, station=SITE["name"], lat=STA[0], lon=STA[1], alt=STA[2],
                  start=str(np.datetime64(t0, "D")), end=str(np.datetime64(t1, "D")),
                  zmin=ZMIN, zmax=ZMAX, target=TARGET, alpha=ALPHA,
                  cams=os.environ.get("ALC_VAL_CAMS_04", ""),
                  sources=list(SOURCES), variants=list(SITE["variants"]),
                  variant_labels=SITE.get("variant_labels", {}),
                  variant_short=SITE.get("variant_short", {}),
                  collapse=bool(SITE.get("collapse_combos", False)),
                  wl_modes=[wl for wl in SITE.get("wl_modes", WL_MODES) if wl in WL_MODES],
                  title=SITE["title"], subtitle=SITE["subtitle"],
                  warnings=SITE.get("warnings", []),
                  profiles_note=SITE.get("profiles_note")),
        alt_agl=alt_agl,
        channels=[dict(label=c["label"], itype=c["itype"], ident=c["ident"], color=c["color"],
                       calib=c["calib"], chan=c.get("chan", c["ident"])) for c in CHANNELS],
        iref=IREF, months=months, combos=combos, calib_series=cal_series, l2_applied=l2c,
        hist_edges=[round(float(x), 3) for x in HIST_EDGES],
        pcolor=pcolor, pwv=pwv,
    )
    out_json = OUT / f"data_{site_key}.json"
    out_json.write_text(json.dumps(payload), encoding="utf-8")
    print(f"-> {out_json} ({out_json.stat().st_size/1e6:.1f} MB)")
    if site_key == "payerne":
        # backward-compat copy under the historic name
        (OUT / "data.json").write_text(json.dumps(payload), encoding="utf-8")


if __name__ == "__main__":
    main()
