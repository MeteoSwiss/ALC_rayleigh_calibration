# -*- coding: utf-8 -*-
"""q19 - build the PROPOSED (not applied) package change from the frozen gate rule.

Reads calibration/rayleigh/{molecular_methods,calibration}.py READ-ONLY, writes patched copies under
C:/DATA/Projects/202610_CHM15k_stability/q19/patch_work/b/ and the unified diff
q19/q19_proposed.patch (paths relative to the repository root). The thresholds come from
q19_gate_variants.json (variant given as argument, default fix_vc2 = the final rule).
Usage: python q19_make_patch.py [variant]
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
Q19 = Path("C:/DATA/Projects/202610_CHM15k_stability/q19")
W = Q19 / "patch_work"
A = W / "a" / "calibration" / "rayleigh"
B = W / "b" / "calibration" / "rayleigh"
A.mkdir(parents=True, exist_ok=True)
B.mkdir(parents=True, exist_ok=True)
for f in ("molecular_methods.py", "calibration.py"):
    shutil.copy(REPO / "calibration" / "rayleigh" / f, A / f)
spec = json.load(open(Q19 / "q19_gate_variants.json"))
name = sys.argv[1] if len(sys.argv) > 1 else "fix_vc2"
G = spec[name]["gate"]

mm = (A / "molecular_methods.py").read_text(encoding="utf-8")


def rep(src, old, new, count=1):
    assert src.count(old) == count, (old[:80], src.count(old))
    return src.replace(old, new)


# --- 1. WindowGrid carries the noise of `signal` ------------------------------------------------
mm = rep(mm, '''    range_alc: Optional[NDArray[np.float64]] = None  # range (m)
''', '''    range_alc: Optional[NDArray[np.float64]] = None  # range (m)
    sigma: Optional[NDArray[np.float64]] = None      # 1-sigma noise of `signal` per gate (None = unknown);
                                                     # read by the tier-2 vertical-consistency test
''')
mm = rep(mm, '''        range_alc=np.asarray(range_alc, float),
    )
''', '''        range_alc=np.asarray(range_alc, float),
        sigma=None if sig_n is None else np.asarray(sig_n, float),
    )
''')

# --- 2. defaults -----------------------------------------------------------------------------------
mm = rep(mm, '''                       use_noise=True, max_chi2red=2.5, ref_chi2max=3.0, ref_pct=10.0,
                       max_ratio_std_excess=0.40, max_temporal_cv_excess=0.8),''',
         f'''                       use_noise=True, max_chi2red=2.5, ref_chi2max=3.0, ref_pct=10.0,
                       max_ratio_std_excess=0.40, max_temporal_cv_excess=0.8,
                       # Tier-2 vertical-consistency test (q19, 2026-10): a rescued window must not be
                       # brighter in R = signal/p_mol than any {G.get("block_m", 750):g}-m block of the column
                       # above it or of the 1 km below it by more than vc_eps_tol + vc_k sigma, nor darker than
                       # a block of that 1 km below by more than vc_tol_lo + vc_k sigma (aerosol tail running
                       # into the window). See _vertical_consistency.
                       use_vc=True, vc_k={G["k"]:g}, vc_eps_tol={G["eps_tol"]:g}, vc_tol_lo={G["tol_lo"]:g},
                       vc_s_max={("np.inf" if G["s_max"] == float("inf") else repr(G["s_max"]))}, vc_block_m={G.get("block_m", 750.0):g}),''')

# --- 3. the test itself (before _select_optimal) ---------------------------------------------------
vc_code = '''def _vertical_consistency(g: WindowGrid, i: int, j: int, vc_k=2.0, vc_eps_tol=0.10, vc_tol_lo=0.10,
                          vc_s_max=np.inf, vc_block_m=750.0, vc_snr_block=3.0, vc_z_cap_m=9000.0,
                          vc_z_floor_m=1500.0, vc_below_m=1000.0, vc_min_gates=10, **_):
    """Is window (i, j) consistent with the rest of the column? -> (ok, verdict).

    Physics: aerosol only ADDS backscatter, and at ceilometer wavelengths the two-way extinction of a
    residual layer inside the 2-9 km band is negligible (R-1 = 0.3 at 4 km -> alpha ~ 1e-3 /km at 1064 nm),
    so on any night the cleanest part of the column has the LOWEST R = signal/p_mol. A window whose R
    exceeds that of some block of the column above it (or within 1 km below it) by more than its noise
    allows contains aerosol (or the profile is distorted in height), and its constant is biased.
    Verdicts: 'accept'; 'aerosol' (some block, above the window or within 1 km below it, darker than the
    window by > vc_eps_tol + vc_k sigma); 'tail' (a block within 1 km below brighter than the window by
    > vc_tol_lo + vc_k sigma: the window sits in an aerosol tail still decaying into it -- the windows that
    the strict tier rejects are typically just above such a tail, +3..5 % on average); 'unverifiable'
    (no block at all; or, if vc_s_max is finite, no block >= 1 km above it measured to vc_s_max).
    The usable column ends where the EXPECTED 300-m block SNR (H0: R = window level) drops below
    vc_snr_block, <= vc_z_cap_m. A block ABOVE that is brighter (cirrus, an elevated layer) never rejects.
    sigma is the measured noise of the night mean (WindowGrid.sigma); with the propagated sigma of
    eprof_v2.2 before 2026-10 (x4.2-4.4 too large) the test loses its power (q19).
    """
    if g.sigma is None or g.signal is None or g.p_mol is None or g.range_alc is None:
        return True, "no-noise"                       # nothing measured -> behave as before
    z = np.asarray(g.range_alc, float)
    sig = np.asarray(g.signal, float)
    p = np.asarray(g.p_mol, float)
    sg = np.asarray(g.sigma, float)
    ok = np.isfinite(sig) & np.isfinite(p) & (p > 0) & np.isfinite(sg) & (sg > 0)
    dz = float(np.median(np.diff(z))) if z.size > 1 else 1.0
    cb = int(round(g.center_m[i] / dz))
    hb = int(round(g.half_m[j] / dz))
    s, e = max(cb - hb, 0), min(cb + hb, z.size)
    win = np.zeros(z.size, bool)
    win[s:e] = True
    win &= ok
    if win.sum() < 3:
        return False, "unverifiable"
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(ok, sig / p, np.nan)
        w = np.where(ok, (p / np.where(ok, sg, 1.0)) ** 2, 0.0)      # 1/sigma_R^2
    sw = float(np.sum(w[win]))
    if not sw > 0:
        return False, "unverifiable"
    Rw = float(np.sum(w[win] * R[win]) / sw)
    sRw = 1.0 / np.sqrt(sw)
    if not Rw > 0:
        return False, "unverifiable"
    zs, ze = float(z[s]), float(z[e - 1])
    # top of the usable column (11-gate running median of sigma keeps one noisy gate from cutting it)
    h = 5
    sgp = np.concatenate([np.full(h, np.nan), np.where(ok, sg, np.nan), np.full(h, np.nan)])
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sg_s = np.nanmedian(np.lib.stride_tricks.sliding_window_view(sgp, 2 * h + 1), axis=1)
        snr_blk = Rw * p / sg_s * np.sqrt(max(1, int(round(300.0 / dz))))
    above = np.where(z > ze)[0]
    bad = ~(np.isfinite(snr_blk[above]) & (snr_blk[above] >= vc_snr_block) & (z[above] <= vc_z_cap_m))
    n_up = int(np.argmax(bad)) if bad.any() else above.size
    z_hi = float(z[above[n_up - 1]]) if n_up > 0 else ze
    nb = max(1, int(round(vc_block_m / dz)))
    up = np.where(ok & (z > ze) & (z <= z_hi))[0]
    blocks = [up[k:k + nb] for k in range(0, up.size, nb)]
    if len(blocks) > 1 and blocks[-1].size < vc_min_gates:
        blocks[-2] = np.concatenate([blocks[-2], blocks[-1]])
        blocks = blocks[:-1]
    is_up = [True] * len(blocks)
    lo = np.where(ok & (z >= max(vc_z_floor_m, zs - vc_below_m)) & (z < zs))[0][::-1]
    bl = [lo[k:k + nb] for k in range(0, lo.size, nb)]
    blocks += bl
    is_up += [False] * len(bl)
    exc, sgm, zc, up_flag = [], [], [], []
    for blk, u in zip(blocks, is_up):
        if blk.size < vc_min_gates:
            continue
        swb = float(np.sum(w[blk]))
        if not swb > 0:
            continue
        Rb = float(np.sum(w[blk] * R[blk]) / swb)
        if not Rb > 0:
            continue
        r = Rw / Rb
        exc.append(r - 1.0)
        sgm.append(r * np.hypot(sRw / Rw, (1.0 / np.sqrt(swb)) / Rb))
        zc.append(float(np.mean(z[blk])))
        up_flag.append(u)
    exc, sgm, zc, up_flag = (np.asarray(exc), np.asarray(sgm), np.asarray(zc),
                             np.asarray(up_flag, bool))
    if exc.size == 0:
        return False, "unverifiable"                   # neither above nor below: nothing to compare
    # A window reaching the SNR limit (no block above) is NOT rejected for that alone: its flatness is
    # tested by the tier-2 chi2red with the measured noise, and it is judged on the 1 km below it (the
    # first design rejected it: ~30 % of clean noisy nights lost in the injection self-test, q19).
    far = up_flag & (zc >= ze + 1000.0)
    if np.isfinite(vc_s_max) and not np.any(sgm[far] <= vc_s_max):
        return False, "unverifiable"
    if np.max(exc - vc_k * sgm) > vc_eps_tol:
        return False, "aerosol"
    below = ~up_flag
    if below.any() and np.isfinite(vc_tol_lo) and np.min(exc[below] + vc_k * sgm[below]) < -vc_tol_lo:
        return False, "tail"
    return True, "accept"


'''
mm = rep(mm, '''def _select_optimal(g: WindowGrid,''', vc_code + '''def _select_optimal(g: WindowGrid,''')
mm = rep(mm, '''                    max_ratio_std_excess=0.30, max_temporal_cv_excess=0.5,
                    w_chi2=0.0, **_) -> MethodWindow:''', '''                    max_ratio_std_excess=0.30, max_temporal_cv_excess=0.5,
                    w_chi2=0.0, use_vc=False, **vc) -> MethodWindow:''')

# --- 4. tier-2: the vertical-consistency test joins the eligibility, applied best-first ------------
mm = rep(mm, '''    if not np.any(elig):
        return _fail(g, "optimal", elig, "no eligible window")

    if fallback and fallback_objective == "lowest":''', '''    if not np.any(elig):
        return _fail(g, "optimal", elig, "no eligible window")

    if fallback and use_vc and g.sigma is not None:
        # Tier-2 vertical-consistency test. The fallback admits windows on nights the strict gates
        # rejected, and those nights are aerosol-loaded (q16e: +9.4 % vs the unit's tier-1 nights; q19).
        # Candidates are taken best-first by the SAME objective; the first one consistent with the column
        # is kept, so the test can move a rescue upward out of an aerosol tail instead of losing the night.
        if fallback_objective == "lowest":
            order_key = np.where(elig, g.start_m, np.inf)
            pick = lambda cand: np.unravel_index(int(np.argmin(np.where(cand, order_key, np.inf))),
                                                 cand.shape)
        else:
            n_norm = g.n_pts / (np.nanmax(g.n_pts) if np.nanmax(g.n_pts) > 0 else 1.0)
            tvar_pen = np.where(np.isfinite(g.temporal_cv), g.temporal_cv, 0.0) if have_tvar else 0.0
            rel_pen = np.where(np.isfinite(g.rel_error), g.rel_error, 0.0) / 100.0
            quality = (g.r2 - w_ratio * np.abs(g.scattering_ratio - 1.0)
                       - w_resid * (g.residual_pct / 100.0) - w_snr * g.ratio_std
                       - w_tvar * tvar_pen - w_rel * rel_pen + w_npts * n_norm)
            pick = lambda cand: _argmax2d(np.where(cand, quality, -np.inf))
        cand = elig.copy()
        first = None
        while np.any(cand):
            i, j = (int(x) for x in pick(cand))
            ok_vc, verdict = _vertical_consistency(g, i, j, **vc)
            first = first or verdict
            if ok_vc:
                return _pack(g, i, j, "optimal", elig,
                             f"noise-aware fallback + vertical consistency (first candidate: {first})")
            cand[i, j] = False
        return _fail(g, "optimal", elig,
                     f"tier-2 windows failed the vertical-consistency test (first candidate: {first})")

    if fallback and fallback_objective == "lowest":''')

mm = rep(mm, '''from dataclasses import dataclass
from typing import Optional, Dict, Any
''', '''import warnings
from dataclasses import dataclass
from typing import Optional, Dict, Any
''')
mm = rep(mm, '''# Keyed by E-PROF version; the private selectors keep their algorithm names.''', '''def _binned_mean_noise(stack, min_rows=8, smooth_gates=11):
    """1-sigma noise of the NIGHT-MEAN signal per gate, measured on the binned (masked) stack it is the mean of.

    calibration._sigma_on_fit_grid propagated the native first-difference sigma with 1/sqrt(n_binned)
    only, so sigma_signal was sqrt(300 s / dt_native) too large (x4.2-4.4 measured on 15-s CHM15k, x2.9 on
    CL61 whose native noise is not white) and every noise-relative tier-2 gate (chi2red, both excess
    floors, the de-biased reference) was ~18x too permissive (q19, 2026-10). Measured here directly, no
    white-noise assumption: robust first differences of consecutive binned profiles / sqrt(2) = one
    binned profile, / sqrt(n_valid) = the night mean, 11-gate running median. Consecutive 5-min
    differences also carry some atmospheric change: an upper bound on the noise. None if too short."""
    try:
        st = np.asarray(stack, float)
        if st.ndim != 2 or st.shape[0] < min_rows:
            return None
        d = np.diff(st, axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            med = np.nanmedian(d, axis=0)
            sb = 1.4826 * np.nanmedian(np.abs(d - med[None, :]), axis=0) / np.sqrt(2.0)
            nd = np.sum(np.isfinite(d), axis=0)
            nv = np.sum(np.isfinite(st), axis=0)
            sb = np.where(nd >= min_rows - 1, sb, np.nan)
            h = smooth_gates // 2
            sp = np.concatenate([np.full(h, np.nan), sb, np.full(h, np.nan)])
            sb = np.nanmedian(np.lib.stride_tricks.sliding_window_view(sp, smooth_gates), axis=1)
        out = sb / np.sqrt(np.maximum(nv, 1))
        return out if np.any(np.isfinite(out) & (out > 0)) else None
    except Exception:                      # a diagnostic must never cost a calibration
        return None


# Keyed by E-PROF version; the private selectors keep their algorithm names.''')
mm = rep(mm, '''        masked = np.where(flag, np.nan, np.asarray(signal_stack, float))
        with np.errstate(invalid="ignore"):
            clean_signal = np.nanmean(masked, axis=0)
''', '''        masked = np.where(flag, np.nan, np.asarray(signal_stack, float))
        with np.errstate(invalid="ignore"):
            clean_signal = np.nanmean(masked, axis=0)
        if sigma_signal is not None:
            # the noise of THIS mean, measured on the stack it averages (falls back to the propagated one)
            meas = _binned_mean_noise(masked)
            if meas is not None and np.size(meas) == np.size(sigma_signal):
                sigma_signal = np.where(np.isfinite(meas), meas, np.asarray(sigma_signal, float))
''')
(B / "molecular_methods.py").write_text(mm, encoding="utf-8")

# ------------------------------------------------------------------------------------------------
# calibration.py: measure the night-mean noise on the binned stack the fit uses
# ------------------------------------------------------------------------------------------------
cal = (A / "calibration.py").read_text(encoding="utf-8")
cal = rep(cal, '''    sigma_native = _native_noise(data)
    data = average_ceilometer_data(data, average_time_s=avg_time_s, average_range_m=avg_range_m)''',
          '''    sigma_native = _native_noise(data)
    n_native = len(data.time)
    data = average_ceilometer_data(data, average_time_s=avg_time_s, average_range_m=avg_range_m)
    # native profiles averaged into each fit profile (white-noise fallback of the noise propagation)
    native_per_bin = max(1.0, n_native / max(len(data.time), 1))''')
cal = rep(cal, '''    sigma_signal = _sigma_on_fit_grid(sigma_native, data.range_alc, rcs_use.shape[0], options)''',
          '''    # The night mean averages rcs_use.shape[0] binned profiles of native_per_bin native profiles each:
    # passing only the binned count (as before 2026-10) left sigma sqrt(300 s / dt_native) too large
    # (x4.2-4.4 measured on 15-s CHM15k) and every noise-relative tier-2 gate ~18x too permissive (q19).
    # This white-noise propagation is now only the FALLBACK: molecular_methods.select_molecular_window
    # measures the noise on the binned (masked) stack itself (_binned_mean_noise). Tier 1 never reads it.
    sigma_signal = _sigma_on_fit_grid(sigma_native, data.range_alc,
                                      int(round(rcs_use.shape[0] * native_per_bin)), options)''')
(B / "calibration.py").write_text(cal, encoding="utf-8")
print("patched copy written for variant", name, G)


def _diff(fn):
    r = subprocess.run(["git", "diff", "--no-index", "--src-prefix=x/", "--dst-prefix=y/",
                        (A / fn).relative_to(W).as_posix(), (B / fn).relative_to(W).as_posix()],
                       cwd=str(W), capture_output=True, text=True)
    return r.stdout.replace("x/a/calibration", "a/calibration").replace("y/b/calibration", "b/calibration")


with open(Q19 / "q19_proposed.patch", "w", encoding="utf-8", newline="") as fh:
    fh.write(_diff("calibration.py") + _diff("molecular_methods.py"))
print("written", Q19 / "q19_proposed.patch")
