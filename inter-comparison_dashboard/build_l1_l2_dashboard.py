# -*- coding: utf-8 -*-
"""Payerne L1-vs-L2 dashboard — compute stage.

Builds every combination of the three optional corrections for BOTH sources and dumps one JSON the
static HTML page then toggles through client-side. Nothing is recomputed in the browser.

  LEFT  panel  "L1 + v2 calibration"  : L1 rcs_0 / C_L(t) * 1e6, C_L = the daily Kalman series from
                                        the latest v2.0 calibration NetCDFs (presentation/l1_l2_calib.py)
  RIGHT panel  "L2 as distributed"    : attenuated_backscatter_0 straight out of the L2 product
                                        (CHM15k = v1.0 Rayleigh, CL31 = the 1e8 default, CL61 = Vaisala)

Optional corrections, applied IDENTICALLY to both panels so the panels stay comparable:
  * water vapour        - divide by the two-way WV transmission at the laser line (910 nm units only)
  * wavelength, simple  - single Angstrom exponent (alpha=1) to the 1064 nm reference
  * wavelength, advanced- component-separated molecular/aerosol conversion using the CAMS T/p profile

Both panels are reduced to the SAME hourly grid and the SAME altitude grid, then restricted to the
hours where all three instruments have data in both panels, so every curve and every difference on
the page comes from one strictly paired sample.

Run:  python inter-comparison_dashboard/build_l1_l2_dashboard.py
Out:  C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard/data.json
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

WMO = "0-20000-0-06610"
STA = (46.8137, 6.9425, 491.0)
# Archive read window (what goes in the stream cache). The daily L2 archive carries CHM15k/CL31 back
# to January, but the CL61 daily concatenation only starts ~2026-06-11 (before that there is neither
# a daily file nor a 5-minute granule), and a paired three-instrument comparison can only use hours
# all three have.
START, END = "20260601", "20260814"
# Comparison window. Spans the CL31 optical-block replacement (2026-07-07 ~13:00) deliberately, so
# its impact is visible: pick Jun 2026 in the period selector for the pre-swap state and Jul/Aug for
# after. The operational Kalman follows the step (see l1_l2_calib), so the L1 panel stays valid
# across it. Changing this does NOT need a cache re-read while it sits inside [START, END]; the
# effective start is clamped to where CL61's daily L2 begins (2026-06-11).
WIN_START, WIN_END = "20260601", "20260813"
TARGET, ALPHA = 1064.0, 1.0
ZMIN, ZMAX = RP.ZMIN, RP.ZMAX               # 500-3000 m AGL statistics band
# Display grid. The instruments reach 15.3 km (CHM15k), 15.7 km (CL61) and 7.7 km (CL31), so the
# grid runs to 15 km and each curve is truncated at its own noise floor instead of at a fixed height.
ZTOP_PLOT = 15000.0                         # display range [m AGL]
DZ = 30.0                                   # display altitude step [m] (the L2 native grid)

OUT = Path(r"C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard")            # data + rendered page live outside the repo
CALIB_V2 = CAL.OUT                          # the Kalman CSVs written from the v2 NetCDFs

# One channel per physical instrument. calib= picks which key_of() series the L1 panel loads.
CHANNELS = [
    dict(ident="A", itype="CHM15k", calib="rayleigh", label="CHM15k (A)", color="#1f77b4"),
    dict(ident="B", itype="CL31",   calib="cloud",    label="CL31 (B)",   color="#ff7f0e"),
    dict(ident="C", itype="CL61",   calib="cloud",    label="CL61 (C)",   color="#2ca02c"),
    # The CL61 carries a Rayleigh calibration as well as the cloud one; it is a separate retrieval
    # of the same constant, and the only CL61 channel a molecular-gate change can move.
    dict(ident="C", itype="CL61",   calib="rayleigh", label="CL61 (C, Rayleigh)", color="#17becf",
         chan="Cr"),
]
IREF = 0                                    # CHM15k is the 1064 nm reference everywhere

WL_MODES = ["none", "angstrom", "molecular"]
# Calibration variant is a third dimension: it changes the CONSTANT the L1 panel divides by, so it
# cannot be applied after the hourly medians are taken (median(beta/c) != median(beta)/median(c)).
# The L2 panel is independent of it and is recomputed identically for each variant -- wasteful, but
# it keeps one combo key for the whole payload and the pool absorbs it.
CALIB_VARIANTS = {"v2.0": CAL.OUT, "v2.2": CAL.OUT_V22}
COMBOS = [(cal, wv, wl) for cal in CALIB_VARIANTS for wv in (False, True) for wl in WL_MODES]


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
def _grid(items):
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
    also the project's validated choice for a bias comparison (noise filter 'none')."""
    gridded = []
    for it in items:
        l1 = it["l1"]
        cbh = l1["cbh"]
        cbh_low = (np.nanmin(np.where((cbh > 0) & (cbh < 20000), cbh, np.nan), axis=1)
                   if cbh.ndim == 2 else cbh)
        g, arrs = IC.retime_hourly(l1["time"], [it["beta_scr"], cbh_low],
                                   min_cov_s=IC.MIN_AVG_S, snr_idx=())
        gridded.append(dict(l1=l1, grid=g, scr=arrs[0]))

    union = np.unique(np.concatenate([g["grid"] for g in gridded]))
    for g in gridded:
        idx = {t: i for i, t in enumerate(g["grid"])}
        pos = np.array([idx.get(t, -1) for t in union])
        g["scrU"] = RP._reindex(g["scr"], pos, union.size)
    # Deliberately NOT intercompare.build_common_grid: it takes the INTERSECTION of the channels'
    # vertical extents (z1 = min of the maxima), so CL31's 7.7 km ceiling would truncate CHM15k
    # (15.3 km) and CL61 (15.7 km) as well. Each channel is regridded from its OWN native altitude
    # straight onto the fixed display grid instead, so every instrument keeps its full range and the
    # short one simply carries NaN above its own top.
    return union, [(g["scrU"], g["l1"]["alt"]) for g in gridded]


def build_panel(cache, source, wv, wl, cal="v2.0"):
    """One source + one correction combo, regridded onto the shared display grid.
    Returns (time_sync, beta[nch] on the display grid, alt_asl)."""
    items = [process_channel(cache[(source, ch["ident"])], ch, source, wv, wl, cal)
             for ch in CHANNELS]
    tsync, streams = _grid(items)
    alt = np.arange(0.0, ZTOP_PLOT + DZ, DZ) + STA[2]          # ASL display grid
    return tsync, [IC.regrid(b, a, alt) for b, a in streams], alt


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


def _init(npz_path, t0, t1):
    global _CACHE
    _CACHE = load_cache(npz_path, t0, t1)


def _stats_for(B, sel, band, z_agl):
    """Every profile / difference / statistic for ONE subset of the paired hours.

    Curves are emitted UNTRUNCATED, with a per-gate `keep` mask alongside: the page applies the
    noise-floor filter client-side so it can be switched off and the raw median inspected."""
    out = {}
    nh = int(sel.sum())
    Bs = {s: [b[sel] for b in B[s]] for s in ("L1", "L2")}
    keep = {}
    for s in ("L1", "L2"):
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
    for s in ("L1", "L2"):
        out[s]["diff_vs_ref"] = [_clean(_rel_profile(Bs[s][k], Bs[s][IREF]))
                                 for k in range(len(CHANNELS))]
        out[s]["diff_keep"] = [[int(x) for x in (keep[s][k] & keep[s][IREF])]
                               for k in range(len(CHANNELS))]
        hists = [_rel_hist(Bs[s][k], Bs[s][IREF], band) for k in range(len(CHANNELS))]
        out[s]["hist"] = [h for h, _ in hists]
        out[s]["hist_med"] = [m for _, m in hists]
    out["l1_vs_l2"] = [_clean(_rel_profile(Bs["L1"][k], Bs["L2"][k])) for k in range(len(CHANNELS))]
    out["l1_vs_l2_keep"] = [[int(x) for x in (keep["L1"][k] & keep["L2"][k])]
                            for k in range(len(CHANNELS))]
    out["n_hours"] = nh
    return out


def _run(spec):
    cal, wv, wl = spec
    out = {}
    panels = {}
    for source in ("L1", "L2"):
        panels[source] = build_panel(_CACHE, source, wv, wl, cal)
    # --- pair the two panels on the hours they share ------------------------------------------
    tL1, tL2 = panels["L1"][0], panels["L2"][0]
    common_t = np.intersect1d(tL1, tL2)
    iL1 = np.searchsorted(tL1, common_t)
    iL2 = np.searchsorted(tL2, common_t)
    alt = panels["L1"][2]
    band = (alt - STA[2] >= ZMIN) & (alt - STA[2] <= ZMAX)
    B = {"L1": [b[iL1] for b in panels["L1"][1]], "L2": [b[iL2] for b in panels["L2"][1]]}
    # hours where every instrument has data in BOTH panels -> one strictly paired sample
    have = np.logical_and.reduce([np.any(np.isfinite(B[s][k][:, band]), axis=1)
                                  for s in ("L1", "L2") for k in range(len(CHANNELS))])
    for s in ("L1", "L2"):
        B[s] = [b[have] for b in B[s]]
    hours = common_t[have]
    z_agl = alt - STA[2]

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
CACHE_NPZ = OUT / "_streams.npz"
_SCALAR = ("station_alt", "lat", "lon", "wavelength", "itype", "wmo", "ident")
_STREAMS = [(s, i) for s in ("L1", "L2")
            for i in dict.fromkeys(ch["ident"] for ch in CHANNELS)]


def build_stream_cache(path=CACHE_NPZ):
    """Read all six streams ONCE from the daily archives (with the 5-minute fallback) and cache
    them, so re-running the correction combos never touches the archive again."""
    flat = {}
    for source, ident in _STREAMS:
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
    """L1 + L2 streams for A/B/C, clipped to [t0, t1] so both panels cover the same period."""
    raw = _read_cache(npz_path)
    cache = {}
    for ch in CHANNELS:
        for source in ("L1", "L2"):
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
    """The period both sources cover for all three instruments, clamped to [WIN_START, WIN_END]."""
    raw = _read_cache(npz_path)
    t0 = max(np.asarray(raw[f"{s}_{i}"]["time"]).min() for s, i in _STREAMS)
    t1 = min(np.asarray(raw[f"{s}_{i}"]["time"]).max() for s, i in _STREAMS)
    t0 = max(t0, np.datetime64(datetime.strptime(WIN_START, "%Y%m%d")))
    t1 = min(t1, np.datetime64(datetime.strptime(WIN_END, "%Y%m%d")) + np.timedelta64(1, "D"))
    return t0, t1


# --------------------------------------------------------------------------- L2 applied constant
def l2_applied_constants(npz_path, t0, t1):
    """Daily median of calibration_constant_0 from the L2 product — the constant the RIGHT panel
    is effectively calibrated with, for the time-series charts."""
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


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if not CACHE_NPZ.exists() or "--reread" in sys.argv[1:]:
        print(f"reading L1 + L2 daily archives {START}..{END} ...", flush=True)
        build_stream_cache()
    npz = str(CACHE_NPZ)
    t0, t1 = common_window(npz)
    print(f"common window: {t0} .. {t1}", flush=True)

    cal_all = CAL.main()                                 # {'v2.0': {...}, 'v2.2': {...}}
    calib, calib22 = cal_all['v2.0'], cal_all['v2.2']
    l2c = l2_applied_constants(npz, t0, t1)

    workers = min(len(COMBOS), max(1, (os.cpu_count() or 8) - 2), 6)
    print(f"computing {len(COMBOS)} correction combos x 2 sources over {workers} workers ...", flush=True)
    combos, months = {}, []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init, initargs=(npz, t0, t1)) as ex:
        for key, res in ex.map(_run, COMBOS):
            combos[key] = res
            months = res["months"]
            full = res[range_key(0, len(months) - 1)]           # the whole window
            st = full["L1"]["stats"]
            print("   %-16s N=%4d h   L1 vs CHM15k: %s" % (
                key, full["n_hours"],
                "  ".join("%s %+.1f%%" % (CHANNELS[k]["label"].split()[0], st[k]["medrelbias_pct"])
                          for k in (1, 2))), flush=True)

    alt_agl = (np.arange(0.0, ZTOP_PLOT + DZ, DZ)).tolist()
    payload = dict(
        meta=dict(wmo=WMO, station="Payerne", lat=STA[0], lon=STA[1], alt=STA[2],
                  start=str(np.datetime64(t0, "D")), end=str(np.datetime64(t1, "D")),
                  zmin=ZMIN, zmax=ZMAX, target=TARGET, alpha=ALPHA,
                  cams=os.environ.get("ALC_VAL_CAMS_04", "")),
        alt_agl=alt_agl,
        channels=[dict(label=c["label"], itype=c["itype"], ident=c["ident"], color=c["color"],
                       calib=c["calib"], chan=c.get("chan", c["ident"])) for c in CHANNELS],
        iref=IREF, months=months, combos=combos, calib=calib, calib22=calib22, l2_applied=l2c,
        hist_edges=[round(float(x), 3) for x in HIST_EDGES],
    )
    (OUT / "data.json").write_text(json.dumps(payload), encoding="utf-8")
    print(f"-> {OUT / 'data.json'} ({(OUT / 'data.json').stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
