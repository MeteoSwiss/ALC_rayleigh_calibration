"""Step 2 forward model: PVC eta and the apparent in-cloud attenuated backscatter around cloud base,
as each Vaisala type would record it, tabulated over (alpha, a_G, CBH) and an adiabatic variant.

Physics (validation/multiple_scattering_eta.py, Hogan 2006 PVC, line-for-line port of
multiscatter 1.2.11 small_angle.c, OPTICS/WAVELENGTH of that module):
  beta_app(r) = (alpha(r) / S) exp(-2 tau(r)) E(r),   S = 18.8 sr,  E = PVC enhancement
  eta         = [2 sum_strips e^-2tau_prev (1 - e^-2 sigma dr)/2 E]^-1 (eta_for_cloud_base, dr = 1 m)
computed on a 1 m grid from the instrument to max(tau = 6, base + 420 m), stored from 60 m below to
300 m above the cloud base (z_rel = -60..300 m).

Instrument response (ASSUMPTIONS, nothing in the repo describes the receiver):
  * pulse: Gaussian in time, FWHM 110 ns (CL31/CL51; brief value) and 160 ns (CL61: CL61 User Guide
    M212475EN-E spec table as quoted by Le et al. 2026 Table 1, see report 03 section 2.1); range
    FWHM = c FWHM_t / 2 (16.5 m / 24.0 m). Every type is ALSO tabulated at 0, 60, 110, 160 ns so
    the step-2 fit can treat the pulse as free. A real diode pulse is asymmetric -- not modelled.
  * receiver: ideal -- no detector/amplifier bandwidth beyond the range-gate integration (each
    sample = mean of the convolved profile over its gate). Flagged assumption.
  * gates: native 10 m (CL31/CL51), 4.8 m (CL61); working grid = the 10 m block average of the
    calibration loader (CL31/CL51: 1 native gate = 10 m; CL61: 2 native gates = 9.6 m). The gate
    phase w.r.t. cloud base is unknown in real data -> 5 phases (gate centres at (k + p) * dg,
    p = 0, 0.2, .., 0.8) are stored, to be averaged or peak-aligned in step 2.

Grid: alpha 3..60 /km (16, log), a_G {2,2.5,3,3.5,4,5,6,8} um, CBH 0.2..3 km (15).
Adiabatic variant: N {50,100,200,400} cm-3, f_ad {0.6,1.0}, same CBH: LWC(h) = f_ad x moist-adiabatic
condensate from a saturated base on the US-standard atmosphere at the CBH (sea-level site), lognormal
sigma = 0.35 droplets -> r_e(h) = r_v exp(sigma^2), a_G(h) = r_e exp(-1.5 sigma^2),
alpha(h) = 2 pi N r_e^2 exp(-3 sigma^2) (same relations as step 1d); the forward-lobe width follows
a_G(h) gate by gate (pvc_enhancement_profile below = the reference code with theta per scattering
gate; verified identical to pvc_enhancement for a constant theta).

Run:  python validation/eta_study/forward_lut.py [--workers 6]
Output: C:/DATA/Projects/202611_eta_study/forward_lut.npz (+ forward_lut_checks.txt)
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import sys
import time as time_mod
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from validation.multiple_scattering_eta import (  # noqa: E402
    OPTICS, WAVELENGTH, pvc_enhancement, eta_for_cloud_base)

OUT = Path("C:/DATA/Projects/202611_eta_study")
TYPES = ("CL31", "CL51", "CL61")
S_LIQ = 18.8
C_LIGHT = 299_792_458.0
DR = 1.0
# Window relative to cloud base. Env overrides (inherited by the spawned workers) build the
# extended table the step-2 fit needs (peak-60..peak+300 m with the model peak ~10-16 m above base):
# ETA_LUT_ZLO=-100 ETA_LUT_ZHI=380 ETA_LUT_NAME=forward_lut_ext.npz
Z_LO = float(os.environ.get("ETA_LUT_ZLO", -60.0))
Z_HI = float(os.environ.get("ETA_LUT_ZHI", 300.0))
LUT_NAME = os.environ.get("ETA_LUT_NAME", "forward_lut.npz")
PAD = 120.0                                    # computed beyond +300 m for the convolution edge
TAU_MAX = 6.0
ALPHA_KM = np.logspace(np.log10(3.0), np.log10(60.0), 16)
AG_UM = np.array([2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0])
CBH_M = np.linspace(200.0, 3000.0, 15)
PULSE_NS = np.array([0.0, 60.0, 110.0, 160.0])
PULSE_NOMINAL_NS = {"CL31": 110.0, "CL51": 110.0, "CL61": 160.0}
NATIVE_GATE_M = {"CL31": 10.0, "CL51": 10.0, "CL61": 4.8}
WORK_BLOCK = {"CL31": 1, "CL51": 1, "CL61": 2}
PHASES = np.arange(5) / 5.0
N_CM3 = np.array([50.0, 100.0, 200.0, 400.0])
F_AD = np.array([0.6, 1.0])
SIGMA_LN = 0.35
Z_REL = np.arange(Z_LO, Z_HI + DR / 2, DR)


# =============================================================================
#  PVC with a per-gate forward-lobe width (adiabatic droplets grow with height)
# =============================================================================
def pvc_enhancement_profile(r, sigma, rho_t, rho_l, theta):
    """pvc_enhancement (validation/multiple_scattering_eta.py) with theta given PER GATE: the lobe
    width of the photons scattered at gate i is theta[i]. Statement order identical to the
    reference port; reduces to it exactly when theta is constant."""
    n = r.size
    dr = r[1] - r[0]
    theta = np.broadcast_to(np.asarray(theta, dtype=float), (n,))
    mu2_t = rho_l ** 2
    fov_factor = 1.0 / (1.0 - np.exp(-(rho_t ** 2) / mu2_t))
    E1 = np.zeros(n); E1s2 = np.zeros(n); E1t2 = np.zeros(n); E1cv = np.zeros(n)
    Em = np.zeros(n); Ems2 = np.zeros(n); Emt2 = np.zeros(n); Emcv = np.zeros(n)
    D = np.zeros(n)
    fs_od = np.cumsum(0.5 * (sigma + np.concatenate(([0.0], sigma[:-1]))) * dr)
    for i in np.nonzero(sigma > 0)[0]:
        th2 = theta[i] ** 2
        ri = r[i]
        w_all = sigma[i] * dr
        w = 0.5 * w_all
        D[i] += 0.5 * sigma[i] * np.exp(-(rho_t ** 2) / mu2_t)
        E1[i] += w
        E1s2[i] += w * mu2_t * ri ** 2
        E1t2[i] += w * (mu2_t + th2)
        E1cv[i] += w * mu2_t * ri
        if E1[i] > 0.0:
            Ems2[i] += w * ((E1s2[i] + Ems2[i]) + 2.0 * (E1cv[i] + Emcv[i]) * 0.0)
            Emt2[i] += w * (E1t2[i] + Emt2[i] + th2 * (E1[i] + Em[i]))
            Emcv[i] += w * (E1cv[i] + Emcv[i])
            Em[i] += w * (E1[i] + Em[i])
        e1, e1s2, e1t2, e1cv = E1[i], E1s2[i], E1t2[i], E1cv[i]
        em, ems2, emt2, emcv = Em[i], Ems2[i], Emt2[i], Emcv[i]
        j = slice(i + 1, n)
        rj = r[j]
        dj = rj - ri
        D[j] += sigma[i] * np.exp(-(rho_t * rj) ** 2 / (th2 * dj ** 2 + mu2_t * rj ** 2))
        E1[j] += w_all
        E1s2[j] += w_all * (mu2_t * rj ** 2 + th2 * dj ** 2)
        E1t2[j] += w_all * (mu2_t + th2)
        E1cv[j] += w_all * (mu2_t * rj + th2 * dj)
        if e1 > 0.0:
            Ems2[j] += w_all * ((e1s2 + ems2) + (e1t2 + emt2 + th2 * (e1 + em)) * dj ** 2
                                + 2.0 * (e1cv + emcv) * dj)
            Emt2[j] += w_all * (e1t2 + emt2 + th2 * (e1 + em))
            Emcv[j] += w_all * (e1cv + emcv + (th2 * (e1 + em) + e1t2 + emt2) * dj)
            Em[j] += w_all * (e1 + em)
    double_ratio = fov_factor * (fs_od - D * dr)
    with np.errstate(divide="ignore", invalid="ignore"):
        fm = np.where(Em > 0, 1.0 - np.exp(-(rho_t * r) ** 2 * Em / np.where(Ems2 > 0, Ems2, 1.0)), 0.0)
    return 1.0 + double_ratio + fov_factor * Em * fm


# =============================================================================
#  Forward model for one cloud
# =============================================================================
def forward_profile(itype, cbh, sigma_of_h, theta_of_h, depth):
    """Return (eta, E(z_rel), beta_app(z_rel) on -60..300+PAD m) for a cloud with base `cbh` and
    extinction / lobe width given as functions of height above base h >= 0."""
    rho_t, rho_l = OPTICS[itype]
    r = np.arange(DR, cbh + depth + DR / 2, DR)
    h = r - cbh
    inside = h >= 0
    sigma = np.where(inside, sigma_of_h(np.maximum(h, 0.0)), 0.0)
    theta = np.where(inside, theta_of_h(np.maximum(h, 0.0)), 1.0)
    if np.ndim(theta_of_h(np.array([0.0, 1.0]))) == 0 or np.ptp(theta[inside]) == 0:
        E = pvc_enhancement(r, sigma, rho_t, rho_l, float(theta[inside][0]))
    else:
        E = pvc_enhancement_profile(r, sigma, rho_t, rho_l, theta)
    tau_prev = np.concatenate(([0.0], np.cumsum(sigma[:-1]) * DR))
    strip = np.exp(-2 * tau_prev) * (1.0 - np.exp(-2 * sigma * DR)) / 2.0
    eta = 1.0 / (2.0 * np.sum(strip * E))
    tau_mid = tau_prev + 0.5 * sigma * DR
    beta = sigma / S_LIQ * np.exp(-2.0 * tau_mid) * E
    zr = np.arange(Z_LO, Z_HI + PAD + DR / 2, DR)
    Ez = np.interp(cbh + zr, r, E, left=1.0, right=np.nan)
    bz = np.interp(cbh + zr, r, beta, left=0.0, right=np.nan)
    return eta, Ez, bz


def convolve_pulse(bz, pulse_ns):
    if pulse_ns <= 0:
        return bz.copy()
    sig = C_LIGHT * pulse_ns * 1e-9 / 2.0 / 2.354820045
    k = np.arange(-int(5 * sig), int(5 * sig) + 1) * DR
    w = np.exp(-0.5 * (k / sig) ** 2)
    w /= w.sum()
    return np.convolve(bz, w, mode="same")


def gate_sample(bz, dg):
    """Gate means of the 1 m profile (z = Z_LO..Z_HI+PAD) for gates of width dg whose centres sit
    at (k + p) dg above the cloud base, p in PHASES; gates kept when centred in [Z_LO, Z_HI]."""
    zr = np.arange(Z_LO, Z_HI + PAD + DR / 2, DR)
    edges = np.concatenate((zr - DR / 2, [zr[-1] + DR / 2]))
    cum = np.concatenate(([0.0], np.cumsum(np.nan_to_num(bz) * DR)))
    k = np.arange(int(np.floor(Z_LO / dg)) - 1, int(np.ceil(Z_HI / dg)) + 2)
    out = []
    for p in PHASES:
        c = (k + p) * dg
        G = np.interp(c + dg / 2, edges, cum) - np.interp(c - dg / 2, edges, cum)
        out.append(G / dg)
    centres = np.array([(k + p) * dg for p in PHASES])
    keep = (k * dg >= Z_LO - dg) & (k * dg <= Z_HI)
    return np.array(out)[:, keep], centres[:, keep]


def post(itype, bz):
    """Pulse convolution then native and working gate sampling for one profile."""
    n_ix = (Z_REL - Z_LO).astype(int)
    conv, nat, wrk = [], [], []
    for pns in PULSE_NS:
        cz = convolve_pulse(bz, pns)
        conv.append(cz[n_ix])
        g, _ = gate_sample(cz, NATIVE_GATE_M[itype])
        nat.append(g)
        g, _ = gate_sample(cz, NATIVE_GATE_M[itype] * WORK_BLOCK[itype])
        wrk.append(g)
    return np.array(conv), np.array(nat), np.array(wrk)


def theta_of(itype, a_um):
    return WAVELENGTH[itype] / (np.pi * a_um * 1e-6)


# =============================================================================
#  Adiabatic cloud microphysics
# =============================================================================
def adiabatic_lwc(cbh, depth, f_ad):
    """LWC(h) [g m-3] on h = 0..depth (1 m) for a saturated parcel lifted from the CBH on the
    US-standard atmosphere (sea-level site), x f_ad."""
    g, cp, Rd, Rv, L = 9.80665, 1004.0, 287.04, 461.5, 2.501e6
    eps = Rd / Rv
    T = 288.15 - 0.0065 * cbh
    p = 101325.0 * (T / 288.15) ** 5.2559

    def rs(T, p):
        es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
        return eps * es / (p - es)

    n = int(round(depth / DR)) + 1
    lwc = np.zeros(n)
    r0 = rs(T, p)
    for i in range(1, n):
        r_s = rs(T, p)
        gam = g * (1 + L * r_s / (Rd * T)) / (cp + L ** 2 * r_s * eps / (Rd * T ** 2))
        rho = p / (Rd * T * (1 + 0.61 * r_s))
        T -= gam * DR
        p -= rho * g * DR
        rho = p / (Rd * T * (1 + 0.61 * rs(T, p)))
        lwc[i] = max(r0 - rs(T, p), 0.0) * rho * 1e3
    return f_ad * lwc


def micro_from_lwc(lwc_gm3, n_cm3):
    """(r_e [m], a_G [m], alpha [m-1]) from LWC and N for lognormal sigma = SIGMA_LN droplets."""
    rv = (3.0 * lwc_gm3 / (4.0 * np.pi * 1.0e6 * n_cm3 * 1e6)) ** (1.0 / 3.0)
    re = rv * np.exp(SIGMA_LN ** 2)
    ag = re * np.exp(-1.5 * SIGMA_LN ** 2)
    alpha = 2.0 * np.pi * n_cm3 * 1e6 * re ** 2 * np.exp(-3.0 * SIGMA_LN ** 2)
    return re, ag, alpha


# =============================================================================
#  Jobs
# =============================================================================
def job_grid(args):
    """One (type, alpha) slice of the homogeneous grid."""
    itype, ia = args
    al = ALPHA_KM[ia] * 1e-3
    depth = max(TAU_MAX / al, Z_HI + PAD) + 2.0
    n_ag, n_c = AG_UM.size, CBH_M.size
    eta = np.zeros((n_ag, n_c))
    E = np.zeros((n_ag, n_c, Z_REL.size), np.float32)
    b1 = np.zeros((n_ag, n_c, Z_REL.size), np.float32)
    conv, nat, wrk = [], [], []
    for ig, ag in enumerate(AG_UM):
        th = theta_of(itype, ag)
        for ic, cbh in enumerate(CBH_M):
            e, Ez, bz = forward_profile(itype, cbh, lambda h: np.full_like(h, al),
                                        lambda h: th, depth)
            n_ix = (Z_REL - Z_LO).astype(int)
            eta[ig, ic] = e
            E[ig, ic] = Ez[n_ix]
            b1[ig, ic] = bz[n_ix]
            c_, n_, w_ = post(itype, bz)
            conv.append(c_); nat.append(n_); wrk.append(w_)
    sh = (n_ag, n_c)
    return (itype, ia, eta, E, b1,
            np.array(conv, np.float32).reshape(sh + conv[0].shape),
            np.array(nat, np.float32).reshape(sh + nat[0].shape),
            np.array(wrk, np.float32).reshape(sh + wrk[0].shape))


def job_adiabatic(args):
    itype, inn, ifa = args
    n_cm3, f_ad = N_CM3[inn], F_AD[ifa]
    res = []
    for cbh in CBH_M:
        depth_max = 3000.0
        lwc = adiabatic_lwc(cbh, depth_max, f_ad)
        re, ag, al = micro_from_lwc(lwc, n_cm3)
        ag = np.maximum(ag, 1.0e-6)                 # lobe width floor at the very base
        hh = np.arange(lwc.size) * DR
        tau = np.cumsum(al) * DR
        depth = min(depth_max, max(float(hh[np.searchsorted(tau, TAU_MAX)]) if tau[-1] > TAU_MAX
                                   else depth_max, Z_HI + PAD) + 2.0)
        th_h = theta_of(itype, ag * 1e6)
        e, Ez, bz = forward_profile(itype, cbh, lambda h: np.interp(h, hh, al),
                                    lambda h: np.interp(h, hh, th_h), depth)
        n_ix = (Z_REL - Z_LO).astype(int)
        c_, n_, w_ = post(itype, bz)
        zz = np.clip(Z_REL, 0, None)
        res.append(dict(eta=e, E=Ez[n_ix], beta=bz[n_ix], conv=c_, nat=n_, wrk=w_,
                        alpha=np.where(Z_REL >= 0, np.interp(zz, hh, al), 0.0),
                        aG=np.where(Z_REL >= 0, np.interp(zz, hh, ag), np.nan),
                        lwc=np.where(Z_REL >= 0, np.interp(zz, hh, lwc), 0.0)))
    return itype, inn, ifa, res


def checks():
    """Sanity: (1) the per-gate-theta PVC equals the reference for a constant theta;
    (2) eta at dr = 1 m vs the operational PVC tables (a_G = 5.5 um, alpha = 10 /km, dr = 2 m)."""
    from calibration.cloud._filters import _ETA_CL31, _ETA_CL51, _ETA_CL61
    lines = []
    rho_t, rho_l = OPTICS["CL31"]
    r = np.arange(1.0, 1400.0, 1.0)
    sig = np.where(r >= 1000, 20e-3, 0.0)
    th = theta_of("CL31", 4.0)
    d = np.max(np.abs(pvc_enhancement_profile(r, sig, rho_t, rho_l, np.full(r.size, th))
                      - pvc_enhancement(r, sig, rho_t, rho_l, th)))
    lines.append(f"per-gate-theta PVC vs reference (constant theta): max |dE| = {d:.3e}")
    for itype, tab in (("CL31", _ETA_CL31), ("CL51", _ETA_CL51), ("CL61", _ETA_CL61)):
        for cbh_km, ref in tab[[1, 3, 6, 9]]:
            e1, _, _ = forward_profile(itype, cbh_km * 1000, lambda h: np.full_like(h, 10e-3),
                                       lambda h: theta_of(itype, 5.5), TAU_MAX / 10e-3 + 2)
            e2 = eta_for_cloud_base(cbh_km * 1000, *OPTICS[itype], theta_of(itype, 5.5), alpha=10e-3)
            lines.append(f"{itype} CBH {cbh_km:.3f} km: table {ref:.5f}  dr=2 m {e2:.5f}  "
                         f"dr=1 m (this LUT) {e1:.5f}  ({100 * (e1 / ref - 1):+.2f} % vs table)")
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time_mod.time()
    chk = checks()
    print("\n".join(chk), flush=True)
    (OUT / LUT_NAME.replace(".npz", "_checks.txt")).write_text("\n".join(chk) + "\n", encoding="utf-8")

    nt, na, ng, nc, nz, npul = len(TYPES), ALPHA_KM.size, AG_UM.size, CBH_M.size, Z_REL.size, PULSE_NS.size
    eta = np.zeros((nt, na, ng, nc))
    E = np.zeros((nt, na, ng, nc, nz), np.float32)
    b1 = np.zeros((nt, na, ng, nc, nz), np.float32)
    conv = np.zeros((nt, na, ng, nc, npul, nz), np.float32)
    nat, wrk = {}, {}
    jobs = [(t, ia) for t in TYPES for ia in range(na)]
    jobs.sort(key=lambda j: j[1])                   # low alpha (deep clouds, slowest) first
    ajobs = [(t, i, j) for t in TYPES for i in range(N_CM3.size) for j in range(F_AD.size)]
    adi = {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs_a = [ex.submit(job_adiabatic, j) for j in ajobs]
        for k, (itype, ia, e, Ez, bz, cz, nz_, wz) in enumerate(ex.map(job_grid, jobs), 1):
            it = TYPES.index(itype)
            eta[it, ia], E[it, ia], b1[it, ia], conv[it, ia] = e, Ez, bz, cz
            nat.setdefault(itype, np.zeros((na,) + nz_.shape, np.float32))[ia] = nz_
            wrk.setdefault(itype, np.zeros((na,) + wz.shape, np.float32))[ia] = wz
            print(f"  grid {k}/{len(jobs)} ({itype}, alpha {ALPHA_KM[ia]:.1f}/km) "
                  f"{time_mod.time() - t0:.0f} s", flush=True)
        for f in futs_a:
            itype, inn, ifa, res = f.result()
            adi[(itype, inn, ifa)] = res
        print(f"  adiabatic {len(ajobs)} jobs done {time_mod.time() - t0:.0f} s", flush=True)

    out = dict(types=np.array(TYPES), alpha_km=ALPHA_KM, aG_um=AG_UM, cbh_m=CBH_M, z_rel_m=Z_REL,
               pulse_ns=PULSE_NS, phases=PHASES, s_liq=np.float64(S_LIQ),
               pulse_nominal_ns=np.array([PULSE_NOMINAL_NS[t] for t in TYPES]),
               eta=eta, E=E, beta_1m=b1, beta_conv=conv,
               adi_N_cm3=N_CM3, adi_f_ad=F_AD, adi_sigma_ln=np.float64(SIGMA_LN))
    for t in TYPES:
        dg = NATIVE_GATE_M[t]
        _, cn = gate_sample(np.zeros(int((Z_HI + PAD - Z_LO) / DR) + 1), dg)
        _, cw = gate_sample(np.zeros(int((Z_HI + PAD - Z_LO) / DR) + 1), dg * WORK_BLOCK[t])
        out[f"{t}_native_gate_m"] = np.float64(dg)
        out[f"{t}_native_centres_m"] = cn            # (phase, gate)
        out[f"{t}_native_beta"] = nat[t]             # (alpha, aG, cbh, pulse, phase, gate)
        out[f"{t}_work_gate_m"] = np.float64(dg * WORK_BLOCK[t])
        out[f"{t}_work_centres_m"] = cw
        out[f"{t}_work_beta"] = wrk[t]
        for key in ("eta", "E", "beta", "conv", "alpha", "aG", "lwc"):
            arr = np.array([[[adi[(t, i, j)][ic][key] for ic in range(nc)]
                             for j in range(F_AD.size)] for i in range(N_CM3.size)])
            out[f"{t}_adi_{key}"] = arr.astype(np.float32) if arr.ndim > 3 else arr
        out[f"{t}_adi_native_beta"] = np.array([[[adi[(t, i, j)][ic]["nat"] for ic in range(nc)]
                                                 for j in range(F_AD.size)] for i in range(N_CM3.size)],
                                               np.float32)
        out[f"{t}_adi_work_beta"] = np.array([[[adi[(t, i, j)][ic]["wrk"] for ic in range(nc)]
                                               for j in range(F_AD.size)] for i in range(N_CM3.size)],
                                             np.float32)
    out["README"] = np.array(__doc__)
    np.savez_compressed(OUT / LUT_NAME, **out)
    print(f"{LUT_NAME} written ({(OUT / LUT_NAME).stat().st_size / 1e6:.1f} MB) "
          f"in {time_mod.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
