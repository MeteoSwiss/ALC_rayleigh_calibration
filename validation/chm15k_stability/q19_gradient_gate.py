# -*- coding: utf-8 -*-
"""q19 - proof of concept: a vertical-consistency (gradient) gate for the eprof_v2.2 fallback tier.

Problem
  eprof_v2.2 = the strict eprof_v2 gates (tier 1) + a noise-relative fallback (tier 2) that is consulted
  only when tier 1 leaves no eligible molecular window. Its aerosol test is self-referential (scattering
  ratio vs the 10th percentile of the SAME night's windows) and chi2red <= 2.5 tolerates aerosol when the
  noise is high. q16e: tier-2 nights are 13 % of the published CHM15k series (33 % in June) and sit +9.4 %
  (median) above the unit's tier-1 nights within +-45 d.

Idea (proposition (c))
  On a molecular night R(z) = signal/p_mol is FLAT above the aerosol: at 1064 nm the two-way extinction of
  a residual layer is negligible (R-1 = 0.3 at 4 km -> alpha ~ 1e-3 /km), so the only way R can change
  with height in and above a valid window is (i) aerosol BACKSCATTER decaying with height (R falls),
  (ii) an elevated layer above (R rises then falls) or (iii) an instrumental distortion (additive baseline
  -> R ~ C + b/p_mol, non-linearity). The gate fits ln R(z) = a + g z over an EXTENDED band around the
  candidate window with the per-gate noise that v2.2 already propagates (sigma_signal), and accepts a
  tier-2 window only when the profile is vertically consistent within noise + a physical tolerance.

Implementation (the calibration package and options.json are NEVER edited)
  Inside this process only:
    * molecular_methods.compute_window_grid -> wrapper that attaches sigma_signal to the WindowGrid;
    * molecular_methods._select_optimal (+ the _SELECTORS entry) -> wrapper that
        - reproduces the operational selection unchanged (tier 1 is bit-identical by construction),
        - records the tier, the profile, the statistics of the selected window and of every relaxed
          (tier-2) eligible window (record mode),
        - in gate mode, on tier-2 nights only, adds the vertical-consistency test to the tier-2
          eligibility and re-selects with the operational objective (max composite quality);
    * calibration.rayleigh.calibration.write_calibration_result -> no-op.
  Everything else (cloud screen, L2-grid binning, native noise, cell flagging, Klett LR 52, the 125-member
  ensemble, all QC) is the operational code.

Sub-commands
  record   : every night of the selected units, operational v2.2 + statistics (no gating)
  inject   : injection self-test on clean tier-1 nights (noise / aerosol / baseline / uniform aerosol)
  gate     : re-run the tier-2 nights with the gate active (exact constants of re-selected windows)
  payerne  : the Payerne CHM15k 2026-09-07 case (operational classification screen reproduced)
  analyse  : tables, figures, summary inputs (q19_analyse.py)

Outputs only under C:/DATA/Projects/202610_CHM15k_stability/q19/.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import OUT, L1_ROOT  # noqa: E402

Q19 = OUT / "q19"
REC_DIR = Q19 / "record"
GATE_DIR = Q19 / "gate"
INJ_DIR = Q19 / "inject"
CENSUS = REPO / "validation/scope_l1_2026_census.json"
EPOCH = datetime(1970, 1, 1)
# CL61 WV correction: 0.4 deg monthly L137 (to 2026-06) then the daily files (never the 1 deg folders).
CAMS_FOLDERS = "D:/CAMS_Monthly_04;D:/CAMS_daily"
# The CL61 manufacturer spectrum is operational since 2026-08 (CLAUDE.md, 2026-08-31 verification).
os.environ.setdefault("ALC_WV_SPECTRUM", '{"CL61": [910.55, 0.188]}')

# --------------------------------------------------------------------------------------------
# Units. CHM15k first: the q16a set (incl. the co-located twins Schiphol A-D, Budapest A/B,
# Magurele A/B, Lindenberg) + the q16e units, spread over 13 countries. Payerne A (abnormal receiver)
# is run but always reported separately. A few CL61 (Rayleigh is operational for them) and two CL31
# (Rayleigh "indicative only") show what the gate does on other instruments.
# --------------------------------------------------------------------------------------------
CHM_UNITS = [
    ("0-20000-0-06240", "A"), ("0-20000-0-06240", "B"), ("0-20000-0-06240", "C"),
    ("0-20000-0-06240", "D"), ("0-20000-0-10393", "0"), ("0-20008-0-UGR", "A"),
    ("0-620-3704-5480", "B"), ("0-20008-0-POT", "A"), ("0-100-20000-0000001", "A"),
    ("0-20008-0-INO", "A"), ("0-20008-0-INO", "B"), ("0-380-5-1", "0"),
    ("0-20000-0-12843", "A"), ("0-20000-0-12843", "B"), ("0-20000-0-10962", "0"),
    ("0-250-1001-07151", "B"), ("0-20000-0-06736", "A"), ("0-20000-0-03590", "A"),
    ("0-20000-0-10035", "0"), ("0-20000-0-01492", "A"),
    ("0-20000-0-12944", "A"), ("0-20000-0-12992", "A"), ("0-20008-0-UEX", "A"),
    ("0-20008-0-UVA", "A"), ("0-348-1-57202", "A"), ("0-724-60000-0000000CIEMAT001", "A"),
]
PAYERNE = ("0-20000-0-06610", "A")
CL61_UNITS = [("0-20000-0-10393", "C"), ("0-20000-0-06610", "C"), ("0-380-5-1", "B"),
              ("0-20000-0-03808", "C"), ("0-20000-0-06447", "B"), ("0-20008-0-BIR", "A")]
CL31_UNITS = [("0-20000-0-06610", "B")]
TWINS = {"0-20000-0-06240": ("A", "B", "C", "D"), "0-20000-0-12843": ("A", "B"),
         "0-20008-0-INO": ("A", "B")}

D_START = datetime(2025, 1, 1)
D_END = datetime(2026, 10, 10)
SUCCESS = (1.0, 0.5)

# --------------------------------------------------------------------------------------------
# The vertical-consistency statistic
# --------------------------------------------------------------------------------------------
# Band definition. The extended band runs from max(Z_FLOOR, window bottom - BELOW) to z_hi, the height
# where the EXPECTED signal-to-noise ratio of a 300 m block (10 fit gates) drops to SNR_BLOCK, computed
# under H0 (R constant = the window level) from the measured per-gate noise -- i.e. it adapts to every
# instrument and night. Z_CAP keeps the band below the region where undetected cirrus and the far-range
# baseline collapse (q16c: X(9-10 km) = 1.03, X(11-12 km) = 0.82) start to matter.
BAND = dict(z_floor_m=1500.0, below_m=1000.0, block_m=300.0, snr_block=3.0, z_cap_m=9000.0)

STAT_KEYS = ("zs", "ze", "z_lo", "z_hi", "n_ext", "n_up", "n_above", "n_below", "Rw", "sRw", "snr_win",
             "g_ext", "sg_ext", "chi_ext", "c_ext", "sc_ext",
             "g_up", "sg_up", "chi_up", "g_lo", "sg_lo",
             "eps_up", "seps_up", "eps_lo", "seps_lo", "chi_flat_up", "chi_flat_ext")


def _running_median(x, n=11):
    """NaN-aware running median (odd window, edges padded with NaN)."""
    from numpy.lib.stride_tricks import sliding_window_view
    x = np.asarray(x, float)
    h = n // 2
    xp = np.concatenate([np.full(h, np.nan), x, np.full(h, np.nan)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmedian(sliding_window_view(xp, n), axis=1)


def _wls(y, X, w):
    """Weighted least squares. Returns (coef, cov, chi2, dof) or None."""
    A = X.T @ (w[:, None] * X)
    try:
        cov = np.linalg.inv(A)
    except np.linalg.LinAlgError:
        return None
    coef = cov @ (X.T @ (w * y))
    res = y - X @ coef
    chi2 = float(np.sum(w * res * res))
    return coef, cov, chi2, int(y.size - X.shape[1])


def _grad(sig, p, z_km, w, zc_km, quad=False):
    """d ln R / dz at zc [% per km] from signal = p*(c0 + c1 (z-zc) [+ c2 (z-zc)^2]) (linear in c).

    The linearised model gives the local log-derivative at zc (g = c1/c0) and keeps the fit linear, so
    the noise propagates analytically. sigma_g is inflated by the Birge ratio sqrt(max(1, chi2red)) so
    unmodelled structure or correlated noise never makes the test over-confident.
    Returns (g, sigma_g, chi2red, curvature, sigma_curvature); curvature in % per km^2."""
    dzk = z_km - zc_km
    cols = [p, p * dzk] + ([p * dzk * dzk] if quad else [])
    X = np.column_stack(cols)
    r = _wls(sig, X, w)
    if r is None:
        return (np.nan,) * 5
    c, cov, chi2, dof = r
    if not (c[0] > 0) or dof < 3:
        return (np.nan,) * 5
    chi2red = chi2 / dof
    birge = np.sqrt(max(1.0, chi2red))
    g = c[1] / c[0]
    var_g = (cov[1, 1] - 2 * g * cov[0, 1] + g * g * cov[0, 0]) / c[0] ** 2
    cur, scur = np.nan, np.nan
    if quad:
        k = c[2] / c[0]
        var_k = (cov[2, 2] - 2 * k * cov[0, 2] + k * k * cov[0, 0]) / c[0] ** 2
        cur, scur = 100 * k, 100 * np.sqrt(max(var_k, 0.0)) * birge
    return 100 * g, 100 * np.sqrt(max(var_g, 0.0)) * birge, chi2red, cur, scur


def _wmean(R, w, m):
    sw = float(np.sum(w[m]))
    if sw <= 0:
        return np.nan, np.nan
    return float(np.sum(w[m] * R[m]) / sw), float(1.0 / np.sqrt(sw))


def prep_profile(signal, p_mol, rng, sigma):
    """Per-night quantities shared by every candidate window (computed once)."""
    if sigma is None:
        return None
    sig = np.asarray(signal, float)
    p = np.asarray(p_mol, float)
    z = np.asarray(rng, float)
    sg = np.asarray(sigma, float)
    if not (sig.size == p.size == z.size == sg.size) or z.size < 20:
        return None
    ok = np.isfinite(sig) & np.isfinite(p) & (p > 0) & np.isfinite(sg) & (sg > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(ok, sig / p, np.nan)
        sR = np.where(ok, sg / p, np.nan)
        wR = np.where(ok, 1.0 / sR ** 2, 0.0)                         # weights of R
        wS = np.where(ok, 1.0 / np.where(ok, sg, 1.0) ** 2, 0.0)      # weights of the signal
    return dict(sig=sig, p=p, z=z, zk=z / 1e3, sg=sg, ok=ok, R=R, sR=sR, wR=wR, wS=wS,
                sg_s=_running_median(np.where(ok, sg, np.nan), 11),
                dz=float(np.median(np.diff(z))))


def window_stats(P, s_idx, e_idx, band=None):
    """Vertical-consistency statistics of one candidate window [s_idx, e_idx) of the fit grid.

    P = prep_profile(...) of the time-cleaned night-mean signal (rcs/r^2), the molecular power, range
    (m AGL) and the 1-sigma noise of the signal per fit gate, exactly as the window search sees them.
    Bands: window [zs, ze]; below [z_lo, zs) with z_lo = max(1.5 km, zs - 1 km); above (ze, z_hi] with z_hi
    where the expected 300-m block SNR under H0 (R = window level) falls below 3 (cap 9 km).
    Gradients in % per km, excesses as fractions (window / band - 1).
    """
    b = dict(BAND, **(band or {}))
    out = dict.fromkeys(STAT_KEYS, np.nan)
    if P is None:
        return out
    z, ok = P["z"], P["ok"]
    s_idx, e_idx = max(int(s_idx), 0), min(int(e_idx), z.size)
    if e_idx - s_idx < 3:
        return out
    zs, ze = float(z[s_idx]), float(z[e_idx - 1])
    out.update(zs=zs, ze=ze)
    win = np.zeros(z.size, bool)
    win[s_idx:e_idx] = True
    win &= ok
    if win.sum() < 3:
        return out
    R, wR = P["R"], P["wR"]
    Rw, sRw = _wmean(R, wR, win)
    if not (np.isfinite(Rw) and Rw > 0):
        return out
    out.update(Rw=Rw, sRw=sRw, snr_win=Rw / sRw)
    nblk = max(1, int(round(b["block_m"] / P["dz"])))
    with np.errstate(divide="ignore", invalid="ignore"):
        snr_blk = Rw * P["p"] / P["sg_s"] * np.sqrt(nblk)
    above = np.where(z > ze)[0]
    bad = ~(np.isfinite(snr_blk[above]) & (snr_blk[above] >= b["snr_block"]) & (z[above] <= b["z_cap_m"]))
    n_ok_above = int(np.argmax(bad)) if bad.any() else above.size
    z_hi = float(z[above[n_ok_above - 1]]) if n_ok_above > 0 else ze
    z_lo = max(b["z_floor_m"], zs - b["below_m"])
    out.update(z_lo=z_lo, z_hi=z_hi)
    sig, p, zk, wS = P["sig"], P["p"], P["zk"], P["wS"]
    zc = 0.5 * (zs + ze) / 1e3
    ext = ok & (z >= z_lo) & (z <= z_hi)
    up = ok & (z >= zs) & (z <= z_hi)
    lo_w = ok & (z >= z_lo) & (z <= ze)
    ab = ok & (z > ze) & (z <= z_hi)
    bl = ok & (z >= z_lo) & (z < zs)
    out.update(n_ext=int(ext.sum()), n_up=int(up.sum()), n_above=int(ab.sum()), n_below=int(bl.sum()))
    if ext.sum() >= 10:
        _, _, _, cur, scur = _grad(sig[ext], p[ext], zk[ext], wS[ext], zc, quad=True)
        g1, sg1, chi1, _, _ = _grad(sig[ext], p[ext], zk[ext], wS[ext], zc)
        out.update(g_ext=g1, sg_ext=sg1, chi_ext=chi1, c_ext=cur, sc_ext=scur)
        res = (R[ext] - Rw) * np.sqrt(wR[ext])
        out["chi_flat_ext"] = float(np.sum(res * res) / max(ext.sum() - 1, 1))
    if up.sum() >= 10:
        g, sgg, chi, _, _ = _grad(sig[up], p[up], zk[up], wS[up], zc)
        res = (R[up] - Rw) * np.sqrt(wR[up])
        out.update(g_up=g, sg_up=sgg, chi_up=chi,
                   chi_flat_up=float(np.sum(res * res) / max(up.sum() - 1, 1)))
    if lo_w.sum() >= 10 and bl.sum() >= 5:
        g, sgg, _, _, _ = _grad(sig[lo_w], p[lo_w], zk[lo_w], wS[lo_w], zc)
        out.update(g_lo=g, sg_lo=sgg)
    if ab.sum() >= 5:
        Ra, sRa = _wmean(R, wR, ab)
        if np.isfinite(Ra) and Ra > 0:
            e = Rw / Ra
            out.update(eps_up=e - 1.0, seps_up=e * np.hypot(sRw / Rw, sRa / Ra))
    if bl.sum() >= 5:
        Rl, sRl = _wmean(R, wR, bl)
        if np.isfinite(Rl) and Rl > 0:
            e = Rw / Rl
            out.update(eps_lo=e - 1.0, seps_lo=e * np.hypot(sRw / Rw, sRl / Rl))
    return out


def vertical_stats(signal, p_mol, rng, sigma, s_idx, e_idx, band=None):
    """Convenience: prep + window_stats for a single window."""
    return window_stats(prep_profile(signal, p_mol, rng, sigma), s_idx, e_idx, band)


# Cleanest-block test. Aerosol only ADDS backscatter and, at 1064 nm, the two-way extinction of a residual
# layer inside the band is negligible (R-1 = 0.3 at 4 km -> alpha ~ 1e-3 /km), so on any night the
# CLEANEST part of the column has the lowest R = signal/p_mol, and a window whose R exceeds that of some
# block of the column above it by more than the noise contains aerosol (or the profile is distorted).
# Blocks of BLOCK_M are formed from the window top up to z_hi (SNR-limited, <= 9 km); optionally also
# below the window down to max(1.5 km, zs - 1 km) (a boundary layer there only RAISES R, so it cannot
# trigger the test; a too-low R there would mean an overlap/near-range defect).
BLOCK_M = 750.0


def block_test(P, s_idx, e_idx, k=2.0, block_m=BLOCK_M, include_below=False, band=None, min_gates=10):
    """Window excess over every block of the column above it (and optionally below it).

    Returns dict: n_blk, lcb (max over blocks of excess - k sigma), lcb_exc / lcb_z / lcb_s (excess,
    height and sigma of that block), ucb (min over blocks of excess + k sigma), best_s (smallest sigma of
    an excess among blocks >= 1 km above the window top), eps_top / seps_top (excess over the top 1.5 km).
    Excess = R_window / R_block - 1 (fraction); positive = window brighter than the block.
    """
    out = dict(n_blk=0, n_up_blk=0, lcb=np.nan, lcb_exc=np.nan, lcb_z=np.nan, lcb_s=np.nan, ucb=np.nan,
               ucb_exc=np.nan, ucb_z=np.nan, ucb_lo=np.nan, ucb_lo_exc=np.nan, best_s=np.nan,
               eps_top=np.nan, seps_top=np.nan)
    st = window_stats(P, s_idx, e_idx, band)
    if P is None or not np.isfinite(st["Rw"]):
        return out, st
    z, ok, R, wR = P["z"], P["ok"], P["R"], P["wR"]
    Rw, sRw, ze, zs, z_hi, z_lo = st["Rw"], st["sRw"], st["ze"], st["zs"], st["z_hi"], st["z_lo"]
    nb = max(1, int(round(block_m / P["dz"])))
    groups = []
    up = np.where(ok & (z > ze) & (z <= z_hi))[0]
    gu = [up[i:i + nb] for i in range(0, up.size, nb)]
    if len(gu) > 1 and gu[-1].size < min_gates:
        gu[-2] = np.concatenate([gu[-2], gu[-1]])
        gu = gu[:-1]
    groups += [(g, True) for g in gu if g.size >= min_gates]
    if include_below:
        lo = np.where(ok & (z >= z_lo) & (z < zs))[0][::-1]
        gl = [lo[i:i + nb] for i in range(0, lo.size, nb)]
        groups += [(g, False) for g in gl if g.size >= min_gates]
    if not groups:
        return out, st
    exc, sig, zc, isup = [], [], [], []
    for g, u in groups:
        w = wR[g]
        sw = float(np.sum(w))
        if sw <= 0:
            continue
        Rb = float(np.sum(w * R[g]) / sw)
        sRb = float(1.0 / np.sqrt(sw))
        if not Rb > 0:
            continue
        r = Rw / Rb
        exc.append(r - 1.0)
        sig.append(r * np.hypot(sRw / Rw, sRb / Rb))
        zc.append(float(np.mean(z[g])))
        isup.append(u)
    if not exc:
        return out, st
    exc, sig, zc, isup = np.array(exc), np.array(sig), np.array(zc), np.array(isup, bool)
    lcb = exc - k * sig
    ucb = exc + k * sig
    i = int(np.argmax(lcb))
    out.update(n_blk=int(exc.size), n_up_blk=int(isup.sum()), lcb=float(lcb[i]), lcb_exc=float(exc[i]),
               lcb_z=float(zc[i]), lcb_s=float(sig[i]))
    if isup.any():
        j = int(np.argmin(np.where(isup, ucb, np.inf)))
        far = isup & (zc >= ze + 1000.0)
        out.update(ucb=float(ucb[j]), ucb_exc=float(exc[j]), ucb_z=float(zc[j]),
                   best_s=float(np.min(sig[far])) if far.any() else np.nan)
    if (~isup).any():
        j = int(np.argmin(np.where(~isup, ucb, np.inf)))
        out.update(ucb_lo=float(ucb[j]), ucb_lo_exc=float(exc[j]))
    top = ok & (z > max(ze, z_hi - 1500.0)) & (z <= z_hi)
    if top.sum() >= min_gates:
        Rt, sRt = _wmean(R, wR, top)
        if Rt > 0:
            e = Rw / Rt
            out.update(eps_top=e - 1.0, seps_top=e * np.hypot(sRw / Rw, sRt / Rt))
    return out, st


# --------------------------------------------------------------------------------------------
# Gate decision (thresholds frozen after the record-run design step; see q19_summary.md)
# --------------------------------------------------------------------------------------------
GATE = dict(
    stat="up",           # gradient band: "up" = window bottom .. z_hi (default); "ext" = from zb - 1 km
    g_ref=0.0,           # reference gradient [%/km] (0 = molecular)
    g_tol=5.0,           # physical tolerance [%/km]
    k=2.0,               # noise allowance in units of the (Birge-inflated) sigma_g
    s_max=np.inf,        # power requirement on sigma_g [%/km] (inf = none)
    min_above_m=0.0,     # the band must reach at least this far above the window top
    reject_rising=True,  # reject significantly RISING R(z) too (elevated layer above / baseline)
)


def gate_decision(st, gate=None):
    """Classify one candidate window from its statistics (see q19_lib.rule_verdict, same logic).

    Returns (accept: bool, verdict: str) with verdict in 'accept', 'reject_aerosol' (R falls with
    height), 'reject_rising' (R rises: elevated layer above or instrumental baseline),
    'reject_unverifiable' (no band above the window, or sigma_g above the power limit)."""
    G = dict(GATE, **(gate or {}))
    key = G["stat"]
    if key == "clean":
        return _clean_decision(st, G)
    g, sg = st.get(f"g_{key}", np.nan), st.get(f"sg_{key}", np.nan)
    above = st.get("z_hi", np.nan) - st.get("ze", np.nan)
    if not (np.isfinite(g) and np.isfinite(sg)):
        return False, "reject_unverifiable"
    if G["min_above_m"] > 0 and not (np.isfinite(above) and above >= G["min_above_m"]):
        return False, "reject_unverifiable"
    if sg > G["s_max"]:
        return False, "reject_unverifiable"
    d = g - G["g_ref"]
    lim = G["g_tol"] + G["k"] * sg
    if d < -lim:
        return False, "reject_aerosol"
    if d > lim:
        return (not G["reject_rising"]), "reject_rising" if G["reject_rising"] else "accept_rising"
    return True, "accept"


def _clean_decision(st, G):
    """Vertical-consistency verdict from block_test() output merged into st (same logic as
    q19_lib.vc_verdict): unverifiable (no block above the window / power), aerosol (some block darker
    than the window by > eps_tol + k sigma), tail (a block within 1 km below brighter by > tol_lo +
    k sigma), [rising: a block above brighter, only if reject_rising], accept."""
    if not np.isfinite(st.get("lcb", np.nan)):
        return False, "reject_unverifiable"          # no block at all
    if G.get("require_above", True) and not (st.get("n_up_blk", 0) or 0) > 0:
        # first design: a window reaching the SNR limit (no block above it) was rejected. The injection
        # self-test showed this rejects ~30 % of clean noisy nights (wide v2.2 windows reach the limit), while
        # the in-window flatness is already tested by tier 2's chi2red once the noise is right -> the final
        # rule (require_above=False) judges such a window on the blocks within 1 km below it.
        return False, "reject_unverifiable"
    if np.isfinite(G.get("s_max", np.inf)) and not (st.get("best_s", np.nan) <= G["s_max"]):
        return False, "reject_unverifiable"
    if st["lcb"] > G["eps_tol"]:
        return False, "reject_aerosol"
    tol_lo = G.get("tol_lo", np.inf)
    if np.isfinite(tol_lo) and np.isfinite(st.get("ucb_lo", np.nan)) and st["ucb_lo"] < -tol_lo:
        return False, "reject_tail"
    if G.get("reject_rising", False) and np.isfinite(st.get("ucb", np.nan)) and st["ucb"] < -G["eps_tol"]:
        return False, "reject_rising"
    return True, "accept"


def candidate_stats(P, s, e, gate=None):
    """Statistics of one candidate window for the gate in use (window_stats [+ block_test])."""
    G = dict(GATE, **(gate or {}))
    if G.get("stat") == "clean":
        bt, st = block_test(P, s, e, k=G["k"], block_m=G.get("block_m", BLOCK_M),
                            include_below=G.get("include_below", False))
        return dict(st, **bt)
    return window_stats(P, s, e)


# --------------------------------------------------------------------------------------------
# Monkeypatching (this process only)
# --------------------------------------------------------------------------------------------
STATE = {"mode": "record", "gate": None, "capture": {}, "cells": True, "sigma_fix": False,
         "native": {}}


def empirical_sigma(stack, sigma_op=None, native=None, bin_s=300.0, min_rows=8):
    """1-sigma noise of the NIGHT-MEAN signal per fit gate, measured on the binned stack the fit uses.

    Why: eprof_v2.2 propagates the native first-difference sigma with /sqrt(n_binned_profiles) only
    (calibration._sigma_on_fit_grid) and never divides out the time averaging INSIDE each 300-s bin, so
    its sigma_signal is sqrt(300 s / dt_native) too large (x4.5 for a 15-s CHM15k; measured x4.2-4.4) and
    every chi2red is ~18x too small. Here the noise is measured directly: robust first differences of
    consecutive binned profiles / sqrt(2) -> one binned profile; / sqrt(n_valid) -> the night mean; then
    an 11-gate running median (the per-gate estimate from ~50-130 differences scatters by ~10 %).
    Consecutive-bin differences also carry 5-min atmospheric variability, so this is an UPPER bound on
    the noise (conservative for a rejection test). Fallback for short stacks: the operational sigma with
    the missing within-bin factor divided out (white-noise assumption).
    """
    if stack is not None and np.ndim(stack) == 2 and np.shape(stack)[0] >= min_rows:
        st = np.asarray(stack, float)
        d = np.diff(st, axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            med = np.nanmedian(d, axis=0)
            sb = 1.4826 * np.nanmedian(np.abs(d - med[None, :]), axis=0) / np.sqrt(2.0)
            nd = np.sum(np.isfinite(d), axis=0)
            nv = np.sum(np.isfinite(st), axis=0)
        sb = np.where(nd >= min_rows - 1, sb, np.nan)
        out = _running_median(sb, 11) / np.sqrt(np.maximum(nv, 1))
        if np.any(np.isfinite(out) & (out > 0)):
            return out, "empirical"
    if sigma_op is not None and native and np.isfinite(native.get("dt_s", np.nan)):
        per_bin = max(1.0, bin_s / float(native["dt_s"]))
        return np.asarray(sigma_op, float) / np.sqrt(per_bin), "native_corrected"
    return None, "none"
_PATCHED = {"done": False}


def _masks(g, kw):
    """Exact replicas of the tier-1 (strict) and tier-2 (relaxed) eligibility of _select_optimal."""
    have_tvar = bool(np.any(np.isfinite(g.temporal_cv)))
    base = (np.isfinite(g.r2) & (g.start_m >= kw["min_window_start_m"]) & (g.slope > 0)
            & (np.abs(g.intercept) < g.slope)
            & (np.isfinite(g.rel_error) & (g.rel_error <= kw["max_rel_error"])))

    def _tvar(mask, excess):
        if not have_tvar:
            return mask
        if excess and np.any(np.isfinite(g.temporal_cv_excess)):
            return mask & (~np.isfinite(g.temporal_cv_excess)
                           | (g.temporal_cv_excess <= kw["max_temporal_cv_excess"]))
        return mask & (~np.isfinite(g.temporal_cv) | (g.temporal_cv <= kw["max_temporal_cv"]))

    with np.errstate(invalid="ignore"):
        strict = (base & (g.r2 >= kw["min_r2"])
                  & (np.isfinite(g.residual_pct) & (g.residual_pct <= kw["max_residual_pct"]))
                  & (g.ratio_std <= kw["max_ratio_std"])
                  & (np.isfinite(g.scattering_ratio) & (g.scattering_ratio <= kw["max_scattering_ratio"])))
        strict = _tvar(strict, excess=False)
        sr = np.where(np.isfinite(g.scattering_ratio_dbz), g.scattering_ratio_dbz, g.scattering_ratio)
        relaxed = (base & (np.isfinite(g.chi2red) & (g.chi2red <= kw["max_chi2red"]))
                   & (~np.isfinite(g.ratio_std_excess) | (g.ratio_std_excess <= kw["max_ratio_std_excess"]))
                   & (np.isfinite(sr) & (sr <= kw["max_scattering_ratio"])))
        relaxed = _tvar(relaxed, excess=True)
    return strict, relaxed


def _quality(g, kw):
    """Exact replica of the composite quality score of _select_optimal (objective 'score')."""
    have_tvar = bool(np.any(np.isfinite(g.temporal_cv)))
    n_norm = g.n_pts / (np.nanmax(g.n_pts) if np.nanmax(g.n_pts) > 0 else 1.0)
    tvar_pen = np.where(np.isfinite(g.temporal_cv), g.temporal_cv, 0.0) if have_tvar else 0.0
    rel_pen = np.where(np.isfinite(g.rel_error), g.rel_error, 0.0) / 100.0
    return (g.r2 - kw["w_ratio"] * np.abs(g.scattering_ratio - 1.0)
            - kw["w_resid"] * (g.residual_pct / 100.0) - kw["w_snr"] * g.ratio_std
            - kw["w_tvar"] * tvar_pen - kw["w_rel"] * rel_pen + kw["w_npts"] * n_norm)


def _cell_idx(g, i, j):
    """Fit-grid index span [s, e) of grid cell (i, j), as compute_window_grid builds it."""
    rng = g.range_alc
    dz = float(np.abs(rng[1] - rng[0])) if len(rng) > 1 else 1.0
    cb = int(round(g.center_m[i] / dz))
    hb = int(round(g.half_m[j] / dz))
    return cb - hb, min(cb + hb, len(rng))


def _install_patches():
    if _PATCHED["done"]:
        return
    import calibration.rayleigh.calibration as RC
    from calibration.rayleigh import molecular_methods as MM

    orig_grid = MM.compute_window_grid
    orig_sel = MM._select_optimal

    orig_sfg = RC._sigma_on_fit_grid

    def sfg_w(sigma_native, range_fit, n_profiles, options):
        # stash the native sampling step (fallback of empirical_sigma); return value unchanged
        STATE["native"] = dict(dt_s=(float(sigma_native["dt_s"]) if sigma_native else np.nan),
                               n_bin=int(n_profiles))
        return orig_sfg(sigma_native, range_fit, n_profiles, options)

    def grid_w(*a, **k):
        sig_op = k.get("sigma_signal")
        emp, how = empirical_sigma(k.get("signal_stack"), sig_op, STATE.get("native"))
        if STATE.get("sigma_fix") and emp is not None and sig_op is not None:
            # variant "fix": the EXISTING tier-2 gates (chi2red, excess floors, de-biased reference) see
            # the correctly scaled noise. Tier 1 never reads sigma, so it stays bit-identical.
            k = dict(k, sigma_signal=np.where(np.isfinite(emp), emp, np.asarray(sig_op, float)))
        g = orig_grid(*a, **k)
        g._q19_sigma = emp                 # the gate statistic ALWAYS uses the measured noise
        g._q19_sigma_op = sig_op
        g._q19_sigma_how = how
        return g

    def sel_w(g, **kw):
        if not kw.get("use_noise", False):
            return orig_sel(g, **kw)                       # plain eprof_v2: untouched
        assert kw.get("objective", "score") == "score" and kw.get("fallback_objective", "score") == "score"
        assert not kw.get("use_bg", False)
        if STATE.get("force_t2"):
            # injection self-test: behave as if the strict tier had found nothing, so the noise-aware
            # tier (and the gate) must handle the night (r2 >= inf is never true -> strict mask empty)
            kw = dict(kw, min_r2=np.inf)
        cap = STATE["capture"]
        mw = orig_sel(g, **kw)                             # the operational v2.2 selection
        strict, relaxed = _masks(g, kw)
        noisy_ok = bool(np.any(np.isfinite(g.chi2red)))
        tier = 1 if strict.any() else (2 if mw.ok else 0)
        sigma = getattr(g, "_q19_sigma", None)
        P = prep_profile(g.signal, g.p_mol, g.range_alc, sigma)
        cap["tier"] = tier
        cap["noisy_ok"] = noisy_ok
        cap["sigma_how"] = getattr(g, "_q19_sigma_how", "")
        sop = getattr(g, "_q19_sigma_op", None)
        if sop is not None and sigma is not None:
            zz = np.asarray(g.range_alc, float)
            m = (zz >= 2000) & (zz <= 5000) & np.isfinite(sigma) & (sigma > 0) & np.isfinite(sop)
            cap["sigma_ratio"] = float(np.median(np.asarray(sop)[m] / sigma[m])) if m.sum() >= 5 else np.nan
        cap["sigma_op"] = None if sop is None else np.asarray(sop, float)
        cap["n_relaxed"] = int(relaxed.sum())
        cap["profile"] = (np.asarray(g.signal, float), np.asarray(g.p_mol, float),
                          np.asarray(g.range_alc, float),
                          None if sigma is None else np.asarray(sigma, float))
        if tier == 2:
            # sanity: the replicated relaxed mask must be the operational tier-2 eligibility
            cap["mask_ok"] = bool(mw.eligible is not None and np.array_equal(mw.eligible, relaxed))
        if mw.ok:
            i = int(np.argmin(np.abs(g.center_m - mw.center_m)))
            j = int(np.argmin(np.abs(g.half_m - mw.half_m)))
            s, e = _cell_idx(g, i, j)
            cap["sel"] = dict(i=i, j=j, s=s, e=e, ratio_med=float(g.ratio_med[i, j]),
                              **window_stats(P, s, e))
        # what the fallback WOULD select had tier 1 failed (on tier-1 nights; same relaxed gates)
        if tier == 1 and noisy_ok and relaxed.any():
            mw2 = orig_sel(g, **dict(kw, min_r2=np.inf))
            if mw2.ok:
                i = int(np.argmin(np.abs(g.center_m - mw2.center_m)))
                j = int(np.argmin(np.abs(g.half_m - mw2.half_m)))
                s, e = _cell_idx(g, i, j)
                cap["t2sel"] = dict(i=i, j=j, s=s, e=e, ratio_med=float(g.ratio_med[i, j]),
                                    **window_stats(P, s, e))
        # every relaxed-eligible cell, with its quality and statistics (offline re-selection)
        if STATE.get("cells", True) and tier == 2 and relaxed.any() and noisy_ok:
            q = _quality(g, kw)
            rows = []
            for i, j in zip(*np.where(relaxed)):
                s, e = _cell_idx(g, i, j)
                st = window_stats(P, s, e)
                rows.append([i, j, s, e, float(g.ratio_med[i, j]), float(q[i, j])]
                            + [st[k] for k in CELL_KEYS])
            cap["cells"] = np.asarray(rows, float)
        # ---- gate mode: tier 2 only -------------------------------------------------------------------
        if STATE["mode"] == "gate" and tier == 2:
            q = _quality(g, kw)
            cand = relaxed.copy()
            tested = []
            while np.any(cand):
                i, j = MM._argmax2d(np.where(cand, q, -np.inf))
                s, e = _cell_idx(g, i, j)
                st = candidate_stats(P, s, e, STATE["gate"])
                acc, verdict = gate_decision(st, STATE["gate"])
                tested.append(verdict)
                if acc:
                    cap["gate"] = dict(n_tested=len(tested), verdict="accept",
                                       first_verdict=tested[0], moved=len(tested) > 1)
                    cap["gsel"] = dict(i=i, j=j, s=s, e=e, ratio_med=float(g.ratio_med[i, j]), **st)
                    return MM._pack(g, i, j, "optimal", relaxed,
                                    "noise-aware fallback + vertical-consistency gate")
                cand[i, j] = False
            cap["gate"] = dict(n_tested=len(tested), verdict=tested[0] if tested else "none",
                               first_verdict=tested[0] if tested else "none", moved=False)
            return MM._fail(g, "optimal", relaxed, "tier-2 windows failed the vertical-consistency gate")
        return mw

    RC._sigma_on_fit_grid = sfg_w
    MM.compute_window_grid = grid_w
    MM._select_optimal = sel_w
    MM._SELECTORS["eprof_v2.2"] = sel_w          # the short path (< 5 profiles) uses the dict
    RC.write_calibration_result = lambda *a, **k: None
    _PATCHED["done"] = True


CELL_KEYS = ("zs", "ze", "z_lo", "z_hi", "g_ext", "sg_ext", "c_ext", "sc_ext", "g_up", "sg_up",
             "g_lo", "sg_lo", "eps_up", "seps_up", "eps_lo", "seps_lo", "chi_flat_up", "chi_flat_ext")
CELL_COLS = ("i", "j", "s", "e", "ratio_med", "quality") + CELL_KEYS

# --------------------------------------------------------------------------------------------
# Harness (same options as the published series: options.json + eprof_v2.2 = calibration_version 220)
# --------------------------------------------------------------------------------------------
_CENSUS = {}


def census_entry(wmo, ident):
    if not _CENSUS:
        for s in json.load(open(CENSUS, encoding="utf-8")):
            _CENSUS[(s["wmo"], s["ident"])] = s
    return _CENSUS[(wmo, ident)]


def make_options():
    from calibration import CalibrationOptions, DataLevel
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.folder_output = Q19 / "nowrite"          # never written (write is patched out)
    o.plot_all = False
    o.plot_main = False
    o.molecular_method = "eprof_v2.2"
    o.cams_folder = CAMS_FOLDERS
    return o


def make_info(s):
    from calibration import InstrumentInfo
    from calibration.config import InstrumentType
    itype = {"CHM15k": InstrumentType.CHM15k, "CL61": InstrumentType.CL61,
             "CL31": InstrumentType.CL31, "CL51": InstrumentType.CL51}[s["type"]]
    return InstrumentInfo(site_name=s.get("site", s["wmo"]), wmo_id=s["wmo"], identifier=s["ident"],
                          instrument_type=itype, latitude=s["lat"], longitude=s["lon"],
                          altitude=s.get("alt", 0.0))


def load_native(s, date_str):
    from calibration.io.data_loader import build_file_paths, load_data
    info = make_info(s)
    o = make_options()
    files = [f for f in build_file_paths(date_str, info, o) if f.exists()]
    if not files:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_data(files, info.instrument_type, o.data_level)


def run_one(s, date_str, native, mode="record", gate=None, classification=None, cells=True,
            sigma_fix=False, force_t2=False):
    """One calibrate_rayleigh call; returns (result dict, capture dict).

    sigma_fix=False: operational eprof_v2.2 (its own, over-estimated sigma_signal feeds the existing
    tier-2 gates); True: the same code with the measured night-mean noise (variant "fix")."""
    _install_patches()
    from calibration import calibrate_rayleigh
    info = make_info(s)
    o = make_options()
    STATE["mode"] = mode
    STATE["gate"] = gate
    STATE["cells"] = cells
    STATE["sigma_fix"] = bool(sigma_fix)
    STATE["force_t2"] = bool(force_t2)
    STATE["capture"] = {}
    STATE["native"] = {}
    fi = {}
    t0 = time.process_time()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = calibrate_rayleigh(date_str, info, o, preloaded_data=native, fit_inputs_out=fi,
                                   classification=classification)
        flag, c, u, msg = float(r.flag), float(r.lidar_constant), float(r.uncertainty), r.message
        zb = (r.calibration_bottom_height - info.altitude) if r.calibration_bottom_height is not None else np.nan
        zt = (r.calibration_top_height - info.altitude) if r.calibration_top_height is not None else np.nan
    except Exception as exc:  # noqa: BLE001 - one bad night must not kill the unit
        flag, c, u, zb, zt = -99.0, np.nan, np.nan, np.nan, np.nan
        msg = f"{type(exc).__name__}: {exc}"
    cap = STATE["capture"]
    res = dict(flag=flag, c=c if flag in SUCCESS else np.nan, u=u if flag in SUCCESS else np.nan,
               zb=zb if flag in SUCCESS else np.nan, zt=zt if flag in SUCCESS else np.nan,
               message=str(msg)[:100], cpu=time.process_time() - t0,
               n_profiles=fi.get("n_profiles", np.nan))
    return res, cap


# --------------------------------------------------------------------------------------------
# record
# --------------------------------------------------------------------------------------------
REC_FIELDS = (["wmo", "ident", "type", "date", "variant", "flag", "c", "u", "zb", "zt", "tier",
               "n_relaxed", "n_profiles", "mask_ok", "sigma_how", "sigma_ratio", "message", "cpu"]
              + [f"sel_{k}" for k in ("s", "e", "ratio_med") + STAT_KEYS]
              + [f"t2_{k}" for k in ("s", "e", "ratio_med") + STAT_KEYS])
Z_KEEP = 12000.0      # profile arrays stored up to 12 km


def _row(s, date_str, res, cap, variant="op"):
    row = dict(wmo=s["wmo"], ident=s["ident"], type=s["type"], date=date_str, variant=variant, **res)
    row.update(tier=cap.get("tier", np.nan), n_relaxed=cap.get("n_relaxed", np.nan),
               mask_ok=cap.get("mask_ok", ""), sigma_how=cap.get("sigma_how", ""),
               sigma_ratio=cap.get("sigma_ratio", np.nan))
    for pre, key in (("sel", "sel"), ("t2", "t2sel")):
        d = cap.get(key)
        if d:
            for k in ("s", "e", "ratio_med") + STAT_KEYS:
                row[f"{pre}_{k}"] = d.get(k, np.nan)
    return row


def _write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})
    os.replace(tmp, path)


def record_chunk(wmo, ident, dates, out_csv):
    """Worker: one unit, a list of nights -> CSV + npz (profiles of every night, relaxed cells)."""
    if out_csv.exists():
        return str(out_csv), 0, 0.0
    warnings.filterwarnings("ignore")
    import logging
    logging.disable(logging.CRITICAL)
    s = census_entry(wmo, ident)
    rows, prof, cells = [], {}, {}
    t0 = time.time()
    for d in dates:
        try:
            native = load_native(s, d)
        except Exception as exc:  # noqa: BLE001
            rows.append(dict(wmo=wmo, ident=ident, type=s["type"], date=d, variant="op", flag=-98,
                             message=f"{type(exc).__name__}: {exc}"[:100]))
            continue
        if native is None:
            continue
        for variant, fix in (("op", False), ("fix", True)):
            try:
                res, cap = run_one(s, d, native, "record", sigma_fix=fix)
            except Exception as exc:  # noqa: BLE001
                rows.append(dict(wmo=wmo, ident=ident, type=s["type"], date=d, variant=variant,
                                 flag=-98, message=f"{type(exc).__name__}: {exc}"[:100]))
                continue
            rows.append(_row(s, d, res, cap, variant))
            if variant == "op" and "profile" in cap and cap["profile"][3] is not None:
                sig, p, z, sg = cap["profile"]
                sop = cap.get("sigma_op")
                sop = sop if (sop is not None and np.size(sop) == z.size) else np.full(z.size, np.nan)
                m = z <= Z_KEEP
                prof[d] = np.vstack([z[m], sig[m], p[m], sg[m], sop[m]]).astype("f8")
            if isinstance(cap.get("cells"), np.ndarray):
                cells[f"{variant}_{d}"] = cap["cells"].astype("f4")
            if res["flag"] in (0.0, -1.0):
                break                     # no data / not clear: identical in both variants
    _write_csv(out_csv, rows, REC_FIELDS)
    np.savez_compressed(out_csv.with_suffix(".npz"),
                        **{f"p_{k}": v for k, v in prof.items()}, **{f"c_{k}": v for k, v in cells.items()})
    return str(out_csv), len(rows), time.time() - t0


def unit_dates(wmo, ident, d0=D_START, d1=D_END):
    out = []
    for k in range((d1 - d0).days + 1):
        d = d0 + timedelta(days=k)
        f = L1_ROOT / wmo / f"{d:%Y}" / f"{d:%m}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"
        if f.exists():
            out.append(f"{d:%Y%m%d}")
    return out


def _pool(fn, tasks, workers):
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(fn, *t) for t in tasks]
        for f in as_completed(futs):
            try:
                p, n, dt = f.result()
            except Exception:  # noqa: BLE001
                traceback.print_exc()
                continue
            done += 1
            print(f"[{done}/{len(tasks)}] {Path(p).parent.name}/{Path(p).name}: {n} rows in {dt:.0f} s "
                  f"(elapsed {time.time() - t0:.0f} s)", flush=True)


def all_units(which="all"):
    u = []
    if which in ("all", "chm"):
        u += CHM_UNITS + [PAYERNE]
    if which in ("all", "cl61"):
        u += CL61_UNITS
    if which in ("all", "cl31"):
        u += CL31_UNITS
    return u


def cmd_record(args):
    units = all_units(args.which) if not args.units else [tuple(x.split(":")) for x in args.units.split(",")]
    tasks = []
    for wmo, ident in units:
        dates = unit_dates(wmo, ident)
        if args.stride > 1:
            dates = dates[::args.stride]
        for m in sorted({d[:6] for d in dates}):
            ds = [d for d in dates if d[:6] == m]
            tasks.append((wmo, ident, ds, REC_DIR / f"{wmo}_{ident}" / f"{m}.csv"))
    # longest units first is irrelevant (month chunks are similar); interleave units for I/O spread
    tasks.sort(key=lambda t: (t[3].name, t[0], t[1]))
    if args.limit:
        tasks = tasks[:args.limit]
    print(f"{len(tasks)} tasks, {sum(len(t[2]) for t in tasks)} unit-nights", flush=True)
    _pool(record_chunk, tasks, args.workers)


# --------------------------------------------------------------------------------------------
# gate: re-run the tier-2 nights with the gate active
# --------------------------------------------------------------------------------------------
GATE_FIELDS = (["wmo", "ident", "type", "date", "variant", "flag", "c", "u", "zb", "zt", "tier",
                "verdict", "first_verdict", "n_tested", "moved", "message", "cpu"]
               + [f"g_{k}" for k in ("s", "e", "ratio_med") + STAT_KEYS])


def gate_chunk(wmo, ident, dates, variants, out_csv):
    """Worker: re-run nights with each gate variant (name, {"base": "op"|"fix", "gate": {...}})."""
    if out_csv.exists():
        return str(out_csv), 0, 0.0
    warnings.filterwarnings("ignore")
    import logging
    logging.disable(logging.CRITICAL)
    s = census_entry(wmo, ident)
    rows = []
    t0 = time.time()
    for d in dates:
        native = load_native(s, d)
        if native is None:
            continue
        for name, spec in variants:
            res, cap = run_one(s, d, native, "gate", gate=spec["gate"], cells=False,
                               sigma_fix=(spec.get("base", "op") == "fix"))
            row = dict(wmo=wmo, ident=ident, type=s["type"], date=d, variant=name, **res,
                       tier=cap.get("tier", np.nan))
            gd = cap.get("gate") or {}
            row.update(verdict=gd.get("verdict", ""), first_verdict=gd.get("first_verdict", ""),
                       n_tested=gd.get("n_tested", ""), moved=gd.get("moved", ""))
            gs = cap.get("gsel") or {}
            for k in ("s", "e", "ratio_med") + STAT_KEYS:
                row[f"g_{k}"] = gs.get(k, np.nan)
            rows.append(row)
    _write_csv(out_csv, rows, GATE_FIELDS)
    return str(out_csv), len(rows), time.time() - t0


def load_gate_variants():
    """Gate variants frozen by the design step in q19_gate_variants.json:
    {name: {"base": "op"|"fix", "gate": {...GATE keys...}}}."""
    f = Q19 / "q19_gate_variants.json"
    if f.exists():
        return [(k, v) for k, v in json.load(open(f, encoding="utf-8")).items()]
    return [("main", {"base": "op", "gate": dict(GATE)})]


def cmd_gate(args):
    rec = load_record()
    t2 = rec[rec.tier == 2]
    variants = load_gate_variants()
    print("variants:", [v[0] for v in variants])
    tasks = []
    for (wmo, ident), g in t2.groupby(["wmo", "ident"]):
        dates = sorted(g.date.unique())
        for k in range(0, len(dates), 25):
            tasks.append((wmo, ident, dates[k:k + 25], variants,
                          GATE_DIR / f"{wmo}_{ident}" / f"chunk{k // 25:03d}.csv"))
    print(f"{len(tasks)} tasks, {t2.date.size} tier-2 unit-night-variants", flush=True)
    _pool(gate_chunk, tasks, args.workers)


def load_record():
    import pandas as pd
    fr = []
    for f in sorted(REC_DIR.glob("*/*.csv")):
        try:
            fr.append(pd.read_csv(f, dtype={"date": str, "ident": str, "wmo": str}))
        except Exception:  # noqa: BLE001
            continue
    df = pd.concat(fr, ignore_index=True)
    df["key"] = df["wmo"] + "_" + df["ident"]
    return df


# --------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("record")
    r.add_argument("--workers", type=int, default=12)
    r.add_argument("--which", default="all")
    r.add_argument("--units", default="")
    r.add_argument("--stride", type=int, default=1)
    r.add_argument("--limit", type=int, default=0)
    g = sub.add_parser("gate")
    g.add_argument("--workers", type=int, default=12)
    t = sub.add_parser("one")
    t.add_argument("unit")
    t.add_argument("date")
    t.add_argument("--mode", default="record")
    t.add_argument("--fix", action="store_true")
    args = ap.parse_args()
    if args.cmd == "record":
        cmd_record(args)
    elif args.cmd == "gate":
        cmd_gate(args)
    elif args.cmd == "one":
        wmo, ident = args.unit.split(":")
        s = census_entry(wmo, ident)
        native = load_native(s, args.date)
        res, cap = run_one(s, args.date, native, args.mode, sigma_fix=args.fix)
        print(res)
        for k, v in cap.items():
            if k in ("profile", "cells", "sigma_op"):
                continue
            if isinstance(v, dict):
                print(k, {kk: (round(float(vv), 3) if isinstance(vv, (float, np.floating)) else vv)
                          for kk, vv in v.items()})
            else:
                print(k, v)


if __name__ == "__main__":
    main()
