# -*- coding: utf-8 -*-
"""Noise-filter (SNR) ingredients for the v3 dashboard — build stage.

The page's noise filter removes data whose signal is not DETECTED at SNR >= 3 over a chosen
evaluation window (5 min / 30 min / 60 min / 3 h).  The payload is hourly, so every sub-hourly
decision — and every decision that combines instruments — is precomputed here on the NATIVE
streams and shipped as per-hour bit masks; the browser only looks bits up.  Four inclusion rules
(the mode dropdown), because SNR screening conditions the sample on signal strength and the
noisy CL31 would otherwise average only its high-aerosol scenes:

  inst    per instrument — each channel screened by its own window SNR (the literal filter;
          deliberately exhibits the selection bias, quantified on-page by the sampling-bias line)
  common  intersection — a (window x gate) pixel survives only if EVERY instrument of the site
          detects it: identical atmospheric sample for all, coverage set by the noisiest unit
  scene   the mask is built from the site's REFEREE (the most sensitive co-located instrument,
          CHM15k per operator decision 2026-08) against one of three thresholds, and applied
          identically to every instrument — the selection depends on the atmosphere, not on each
          instrument's own noise realisation
  p2      "average first, filter after" — nothing is removed before averaging; the browser masks
          the displayed period profile where the AGGREGATE SNR < 3, using the monthly noise sums
          shipped here (sn = sum n, ss = sum n*sigma^2 per month x gate)

Window semantics per (wall-clock aligned) window, replicating intercompare.snr_mask exactly —
signed median over the window, robust sigma = calibration.sensitivity.noise.robust_std, windows
with < 5 samples keep everything — with ONE documented extension: where a MEASURED covered-
telescope dark b(z) exists, the SNR is evaluated on the dark-subtracted signal (median - b), so
an r^2-amplified electronic offset cannot pass as detection.  The user-facing dark checkbox does
NOT change the masks: detection is a property of the photons, not of the display state.

Sub-hourly windows are aggregated to the hour as a sample-weighted fraction of data lying in
detected windows; the hour is admitted when that fraction >= NF_FMIN (0.5).  The exact
"median of surviving samples" alternative would need one full value block per (mode x window)
and is deliberately not shipped; check_v3.py quantifies the approximation.

Scene thresholds (operator decision 2026-08-24): three ABSOLUTE attenuated-backscatter levels,
0.1 / 0.25 / 0.5 Mm-1 sr-1 at the target wavelength, measured AND detected by the referee —
a window passes when the referee's calibrated window median satisfies
beta_att >= max(threshold, 3 sigma_ref/sqrt(n)), i.e. the scene is at least that bright and
the referee actually detects it.  The referee's constants are FROZEN at the site's default
variant so the mask cannot drift with the user's variant clicks.
"""
from __future__ import annotations
import base64
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from validation.paper import intercompare as IC
from calibration.sensitivity.noise import robust_std
import l1_l2_io as IO

NF_WINS = (300, 1800, 3600, 10800)          # 5 min / 30 min / 60 min / 3 h
NF_FMIN = 0.5                               # sub-hourly admission: >= half the hour's data
SNR_MIN = IC.SNR_MIN                        # 3.0 — same gate as the operational product
MIN_SAMPLES = IC._SNR_MIN_SAMPLES           # < 5 samples: too short to estimate noise, keep all
THR_KEYS = ("b010", "b025", "b050")         # bit order in the scene block (thr*4 + win)
THR_BETA = {"b010": 0.10, "b025": 0.25, "b050": 0.50}   # Mm-1 sr-1 at the target wavelength


# ---------------------------------------------------------------------------- encoders
def _bits(A, dtype):
    A = np.ascontiguousarray(A, dtype)
    return dict(b=base64.b64encode(A.tobytes()).decode("ascii"),
                nh=int(A.shape[0]), nz=int(A.shape[1]))


# ---------------------------------------------------------------------------- screening
def nf_screen(beta, cbh, vv, time):
    """intercompare.screen minus the qf step (E-PROFILE L1 carries no qf — the hourly reader
    fills zeros), applied at NATIVE resolution: the same cloud/fog exclusion, expanded by
    +-15 min, that produced the hourly blocks."""
    b = beta.copy()
    has_cloud = np.any(np.isfinite(cbh) & (cbh > 0) & (cbh < 20000), axis=1)
    v = np.asarray(vv, "f8").copy()
    v[v < 0] = np.nan                       # negative = "no fog" sentinel
    has_excl = has_cloud | np.isfinite(v)
    t = pd.to_datetime(time)
    if t.size > 1:
        dt = np.median(np.diff(t.values).astype("timedelta64[s]").astype(float))
        if dt > 0:
            win = int(2 * round(15 * 60 / dt) + 1)
            expanded = pd.Series(has_excl.astype(float)).rolling(
                win, center=True, min_periods=1).max().values > 0
        else:
            expanded = has_excl
    else:
        expanded = has_excl
    b[expanded, :] = np.nan
    return b


# ---------------------------------------------------------------------------- window statistics
def stream_windows(tsec, X, W):
    """Wall-clock W-second windows of one native stream (time-sorted).

    Returns (wids, med, sig, nfin, small): window ids (epoch//W), per-gate signed median,
    robust sigma and finite count (n_w x n_z), and the < MIN_SAMPLES flag per window.
    Vectorised by bucketing windows of equal row count (nanmedian over a 3D stack), chunked to
    bound memory.
    """
    wid = tsec // W
    change = np.flatnonzero(np.diff(wid)) + 1
    starts = np.concatenate([[0], change])
    ends = np.concatenate([change, [wid.size]])
    wids = wid[starts]
    counts = ends - starts
    nw, nz = wids.size, X.shape[1]
    med = np.full((nw, nz), np.nan)
    sig = np.full((nw, nz), np.nan)
    nfin = np.zeros((nw, nz), "f8")
    small = counts < MIN_SAMPLES
    for c in np.unique(counts):
        idx = np.flatnonzero(counts == c)
        chunk = max(1, int(2e8 / (int(c) * nz * 8)))
        for i0 in range(0, idx.size, chunk):
            ii = idx[i0:i0 + chunk]
            S = np.stack([X[starts[j]:ends[j]] for j in ii])
            with np.errstate(all="ignore"):
                med[ii] = np.nanmedian(S, axis=1)
                sig[ii] = robust_std(S, axis=1)
                nfin[ii] = np.isfinite(S).sum(axis=1)
    return wids, med, sig, nfin, small


def window_keep(med, sig, nfin, small, dark_nat):
    """Per-(window, native gate) detection mask — intercompare.snr_mask semantics, evaluated on
    the dark-subtracted signal where a measured b(z) exists."""
    m = med - dark_nat[None, :] if dark_nat is not None else med
    with np.errstate(all="ignore"):
        snr = m / (sig / np.sqrt(np.maximum(nfin, 1)))
    keep = ~np.isfinite(snr) | (snr >= SNR_MIN)
    keep[~np.isfinite(m)] = True            # NaN gates stay NaN anyway; not "removed"
    keep[small] = True
    return keep


# ---------------------------------------------------------------------------- hour aggregation
def _hour_accumulate(union_h, wids, W, num_w, den_w):
    """Sum window-level (num, den) onto the union hour axis.  W <= 3600: each window belongs to
    one hour; W = 10800: each hour inherits its containing 3 h window."""
    n_u, nz = union_h.size, num_w.shape[1]
    num = np.zeros((n_u, nz))
    den = np.zeros((n_u, nz))
    if W <= 3600:
        hpos = {int(h): i for i, h in enumerate(union_h)}
        hh = wids * W // 3600
        for j in range(wids.size):
            u = hpos.get(int(hh[j]))
            if u is None:
                continue
            num[u] += np.nan_to_num(num_w[j])
            den[u] += np.nan_to_num(den_w[j])
    else:
        wpos = {int(w): j for j, w in enumerate(wids)}
        for u in range(n_u):
            j = wpos.get(int(union_h[u] * 3600 // W))
            if j is not None:
                num[u] = np.nan_to_num(num_w[j])
                den[u] = np.nan_to_num(den_w[j])
    return num, den


def _admit(num, den):
    with np.errstate(all="ignore"):
        return (den > 0) & (num / np.where(den > 0, den, 1.0) >= NF_FMIN)


# ---------------------------------------------------------------------------- one stream (worker)
def _stream_stats(site_key, wmo, ident, d0, d1, t0, t1, z_agl, station_alt, use_dark):
    """One instrument's window statistics on the display grid, all four windows.

    Module-level so a ProcessPoolExecutor can run the site's instruments concurrently — the
    window medians are single-core numpy, and the per-window results (display grid, ~175
    gates) are small enough to send back over IPC.  Returns (per_w, printed lines) or None
    when the native stream cannot be read."""
    import build_l1_l2_dashboard as BD
    BD.set_site(site_key)
    lines = []
    d = IO.read_l1_native(wmo, ident, d0, d1)
    if d is None:
        return None
    tarr = np.asarray(d["time"])
    m = (tarr >= t0) & (tarr <= t1)              # same clip as load_cache
    scr = nf_screen(np.asarray(d["beta"], "f4")[m], np.asarray(d["cbh"])[m],
                    np.asarray(d["vv"])[m], tarr[m])
    alt_nat = np.asarray(d["alt"], "f8")
    tsec = tarr[m].astype("datetime64[s]").astype("i8")
    order = np.argsort(tsec, kind="stable")
    tsec, scr = tsec[order], scr[order]
    dt_med = float(np.median(np.diff(tsec))) if tsec.size > 1 else np.nan
    dark_nat = BD.dark_for(ident, alt_nat - station_alt) if use_dark else None
    alt_asl = z_agl + station_alt
    lines.append(f"   nf {ident}: {tsec.size} native profiles, dt={dt_med:.0f}s, "
                 f"{alt_nat.size} gates, dark={'yes' if dark_nat is not None else 'no'}")
    per_w = {}
    for W in NF_WINS:
        wids, med, sig, nfin, small = stream_windows(tsec, scr, W)
        keep = window_keep(med, sig, nfin, small, dark_nat)
        rk = IC.regrid(keep.astype("f8") * nfin, alt_nat, alt_asl)
        rn = IC.regrid(nfin, alt_nat, alt_asl)
        entry = dict(wids=wids, rk=rk, rn=rn,
                     med=IC.regrid(med, alt_nat, alt_asl),
                     sig=IC.regrid(sig, alt_nat, alt_asl))
        with np.errstate(all="ignore"):
            entry["keepb"] = (np.nan_to_num(rn) > 0) & \
                (np.nan_to_num(rk) / np.where(np.nan_to_num(rn) > 0,
                                              np.nan_to_num(rn), 1.0) >= NF_FMIN)
        per_w[W] = entry
        del med, sig, nfin, small, keep
    lines.append(f"   nf {ident}: windows " +
                 " ".join(f"{W}s:{per_w[W]['wids'].size}" for W in NF_WINS))
    return per_w, lines


# ---------------------------------------------------------------------------- main entry
def compute_nf(v3, site_key, wmo, t0, t1, union, have, calib,
               dark_disp, z_agl, station_alt, hour_month, months, use_dark):
    """-> the payload's "nf" dict (see module docstring).

    Reads the NATIVE L1 streams itself (l1_l2_io.read_l1_native): the _streams_* cache is
    hourly by design, and an hourly stream has exactly one sample per window — no noise
    statistics.  The per-instrument statistics run in a process pool (one core each).
    Raises when a stream the site displays cannot be read natively."""
    idents = [i["ident"] for i in v3["instruments"]]
    ref_id = v3.get("nf_ref")
    union_h = union.astype("datetime64[h]").astype("i8")
    n_u, nz = union_h.size, z_agl.size
    nW = len(NF_WINS)

    def cseries(ident):
        """Frozen constants C(t) of the site-default variant, as (kal dates, kal values)."""
        method, variant = v3["default"][ident]
        rec = calib.get(f"{ident}|{method}|{variant}")
        if not (rec and rec.get("ok")):
            return None, f"{method}|{variant}"
        kd = np.array([np.datetime64(x) for x in rec["kal"]["d"]])
        return (kd, np.asarray(rec["kal"]["v"], "f8")), f"{method}|{variant}"

    # --- per-stream window statistics on the display grid, all four windows (parallel) --------
    d0 = str(np.datetime64(t0, "D")).replace("-", "")
    d1 = str(np.datetime64(t1, "D")).replace("-", "")
    stats = {}
    with ProcessPoolExecutor(max_workers=min(len(idents), 6)) as ex:
        futs = {ident: ex.submit(_stream_stats, site_key, wmo, ident, d0, d1,
                                 np.datetime64(t0), np.datetime64(t1), np.asarray(z_agl, "f8"),
                                 float(station_alt), bool(use_dark))
                for ident in idents}
        for ident in idents:
            res = futs[ident].result()
            if res is None:
                raise SystemExit(f"nf: no native L1 for {wmo}_{ident} {d0}..{d1} — "
                                 f"cannot build the noise filter")
            stats[ident], lines = res
            for ln in lines:
                print(ln, flush=True)

    # --- per-instrument admission bits --------------------------------------------------------
    inst_bits = {}
    for ident in idents:
        bits = np.zeros((n_u, nz), "u1")
        for w, W in enumerate(NF_WINS):
            e = stats[ident][W]
            num, den = _hour_accumulate(union_h, e["wids"], W,
                                        e["rk"], e["rn"])
            bits |= (_admit(num, den).astype("u1") << w)
        inst_bits[ident] = bits

    # --- common (intersection) bits -----------------------------------------------------------
    common = np.zeros((n_u, nz), "u1")
    for w, W in enumerate(NF_WINS):
        all_wids = sorted(set().union(*[set(stats[i][W]["wids"].tolist()) for i in idents]))
        nw = len(all_wids)
        wpos = {wd: j for j, wd in enumerate(all_wids)}
        passing = np.ones((nw, nz), bool)
        weight = np.zeros((nw, nz))
        for ident in idents:
            e = stats[ident][W]
            pres = np.zeros((nw, nz), bool)
            kb = np.zeros((nw, nz), bool)
            rows = np.array([wpos[wd] for wd in e["wids"].tolist()])
            rn = np.nan_to_num(e["rn"])
            pres[rows] = rn > 0
            kb[rows] = e["keepb"]
            weight[rows] += rn
            passing &= pres & kb
        num, den = _hour_accumulate(union_h, np.asarray(all_wids, "i8"), W,
                                    passing * weight, weight)
        common |= (_admit(num, den).astype("u1") << w)

    # --- scene bits: the referee's calibrated beta_att vs three ABSOLUTE thresholds -----------
    # A window passes threshold T when the referee both MEASURES beta_att >= T and DETECTS it
    # (>= 3 sigma/sqrt(n) of its own window noise) — "beta_att detected by the reference".
    scene = None
    scene_why = None
    consts = {}
    if ref_id is None:
        scene_why = "nf_ref non déclaré pour ce site"
    else:
        ref_c, ref_lbl = cseries(ref_id)
        consts[ref_id] = ref_lbl
        if ref_c is None:
            scene_why = (f"constantes indisponibles pour la référence {ref_id} ({ref_lbl}) — "
                         f"masque scène impossible")
    if scene_why is None:
        scene = np.zeros((n_u, nz), "u2")
        for w, W in enumerate(NF_WINS):
            e = stats[ref_id][W]
            wids = e["wids"]
            tc = (wids * W + W // 2).astype("datetime64[s]")
            cref = IC.interp_calib(ref_c[0], ref_c[1], tc)
            rn_ref = np.nan_to_num(e["rn"])
            dref = np.asarray(dark_disp[ref_id], "f8") if dark_disp.get(ref_id) else \
                np.zeros(nz)
            with np.errstate(all="ignore"):
                beta_cal = (np.nan_to_num(e["med"]) - dref[None, :]) * 1e6 / cref[:, None]
                det = SNR_MIN * np.nan_to_num(e["sig"]) * 1e6 / cref[:, None] / \
                    np.sqrt(np.maximum(rn_ref, 1))
            for t, key in enumerate(THR_KEYS):
                thrv = np.maximum(THR_BETA[key], det)
                passing = (rn_ref > 0) & (beta_cal >= thrv)
                num, den = _hour_accumulate(union_h, wids, W, passing * rn_ref, rn_ref)
                scene |= (_admit(num, den).astype("u2") << (t * 4 + w))

    # --- p2 monthly noise sums (hourly sigma & n on the paired axis) --------------------------
    p2 = {}
    n_m = len(months)
    hm = np.asarray(hour_month)
    for ident in idents:
        e = stats[ident][3600]
        sigH = np.zeros((n_u, nz))
        nH = np.zeros((n_u, nz))
        hpos = {int(h): i for i, h in enumerate(union_h)}
        hh = e["wids"]                          # W = 3600: wid IS the epoch hour
        for j in range(hh.size):
            u = hpos.get(int(hh[j]))
            if u is None:
                continue
            sigH[u] = np.nan_to_num(e["sig"][j])
            nH[u] = np.nan_to_num(e["rn"][j])
        sigP, nP = sigH[have], nH[have]
        sn = np.zeros((n_m, nz), "f8")
        ss = np.zeros((n_m, nz), "f8")
        for m in range(n_m):
            sel = hm == m
            sn[m] = nP[sel].sum(axis=0)
            ss[m] = (nP[sel] * sigP[sel] ** 2).sum(axis=0)
        p2[ident] = dict(sn=_f32(sn), ss=_f32(ss))

    out = dict(
        wins=list(NF_WINS), fmin=NF_FMIN, snr=float(SNR_MIN),
        ref=ref_id, noisiest=v3.get("nf_noisiest"), consts=consts,
        inst={i: _bits(inst_bits[i][have], "u1") for i in idents},
        common=_bits(common[have], "u1"),
        scene=(None if scene is None else _bits(scene[have], "u2")),
        scene_why=scene_why,
        p2=p2,
    )
    # sanity print: retained fraction per instrument at each window, in the statistics band
    band = (z_agl >= 500.0) & (z_agl <= 3000.0)
    for ident in idents:
        b = inst_bits[ident][have][:, band]
        fr = [(b >> w & 1).mean() * 100 for w in range(nW)]
        print(f"   nf {ident}: retenu bande 0.5-3 km  " +
              "  ".join(f"{W}s {f:.0f}%" for W, f in zip(NF_WINS, fr)), flush=True)
    cb = common[have][:, band]
    print("   nf commun: " + "  ".join(f"{W}s {(cb >> w & 1).mean() * 100:.0f}%"
                                       for w, W in enumerate(NF_WINS)), flush=True)
    if scene is not None:
        sb = scene[have][:, band]
        for t, key in enumerate(THR_KEYS):
            print(f"   nf scène {key}: " +
                  "  ".join(f"{W}s {(sb >> (t * 4 + w) & 1).mean() * 100:.0f}%"
                            for w, W in enumerate(NF_WINS)), flush=True)
    return out


def _f32(A):
    A = np.asarray(A, "f4")
    return dict(b=base64.b64encode(np.ascontiguousarray(A, "<f4").tobytes()).decode("ascii"),
                n=int(A.shape[0]), m=int(A.shape[1]))
