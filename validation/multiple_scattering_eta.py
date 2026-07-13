"""Multiple-scattering correction factor eta(range) for the liquid-cloud calibration.

Recomputes the O'Connor/Hopkin multiple-scattering factor eta -- defined through the
saturated integrated attenuated backscatter of a fully-attenuating liquid cloud,
B = 1 / (2 eta S) -- for every instrument type of the E-PROFILE cloud calibration,
using the Photon Variance-Covariance (PVC) method of Hogan (2006, Appl. Opt. 45,
5984), the same fast multiple-scattering model Hopkin et al. (2019, AMT 12, 4131)
used to derive the operational CL31/CL51 corrections.

Model summary (Hogan 2006, quasi-small-angle / equivalent-medium):
- Half of the droplet extinction is diffracted into a narrow forward lobe,
  approximated as a Gaussian of 1/e angular half-width Theta = lambda / (pi a_G)
  (a_G = equivalent-area droplet radius; ~30-60 mrad for 5-10 um droplets, i.e. much
  WIDER than any ceilometer FOV, which is why the FOV footprint matters).
- The two-way problem is folded into a one-way transport in an "equivalent medium"
  where the forward-lobe scattering rate is the full extinction sigma (2 x sigma/2)
  and the return journey is in vacuum (reciprocity theorem).
- Three photon populations per gate (unscattered u, singly forward-scattered s,
  multiply forward-scattered m), each summarized by relative energy P (vs the
  unscattered beam), lateral second moment s2 = <x^2+y^2>, angular second moment
  th2, and covariance sth (Hogan Eqs. 14-26).
- A telescope of half-angle FOV rho_t accepts photons backscattered at range R with
  lateral offset s <= rho_t R; for a Gaussian population the accepted fraction is
  f = 1 - exp(-(rho_t R)^2 / s2)  (Hogan Eqs. 6 and 11).
- The multiple-scattering enhancement of the apparent backscatter at R is then
  E(R) = [f_u + P_s f_s + P_m f_m] / f_u, and for a deep cloud starting at z_c
  eta(z_c) = [2 * integral( sigma exp(-2 tau) E dz )]^-1.

Validation: with a single droplet size in the 8-20 um DIAMETER band used by
Le & O'Connor (2026) for the CL61, the model must reproduce the operational
CL31/CL51 eta tables (calibration/cloud/calibration.py, from Hopkin et al. 2019)
-- that both validates this implementation and pins the effective droplet size,
which is then reused for the instruments that never had a proper table:
CHM15k, Mini-MPL, MPL (currently mis-corrected, see report).

Instrument optics (half-angles, 1/e):
- CL31   rho_t = 0.83 mrad (Wiegner et al. 2014 Table 1; Vaisala datasheet),
         divergence +-0.4 x +-0.7 mrad -> rho_l = 0.57 mrad, 910 nm
- CL51   rho_t = 0.56 mrad, divergence +-0.15 x +-0.25 -> rho_l = 0.21 mrad, 910 nm
- CL61   rho_t = 0.56 mrad (receiver "field-of-view divergence +-0.56 mrad", Vaisala
         CL61 User Guide M212475EN, receiving specifications -- IDENTICAL to the CL51),
         divergence +-0.2 x +-0.35 mrad (Le & O'Connor 2026 Table 1) -> 0.28 mrad
- CHM15k rho_t = 0.23 mrad (Wiegner et al. 2014 Table 1), divergence 0.15 mrad,
         1064 nm
- Mini-MPL FOV 220 urad (operations manual) -> rho_t = 0.11 mrad, shared-telescope
         divergence ~0.055 mrad, 532 nm
- MPL    FOV ~100 urad (Campbell et al. 2002) -> rho_t = 0.05 mrad, div 0.025, 532 nm

Outputs (landscape figures + the python tables to paste into calibration/cloud):
  C:/DATA/Projects/202606_E-PROFILE_calibration/figs_multiple_scattering/  (or --out)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# The operational tables to validate against (Hopkin et al. 2019 values in the code)
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from calibration.cloud.calibration import _ETA_CL31, _ETA_CL51  # noqa: E402


# --------------------------------------------------------------------------- PVC core
def pvc_enhancement(r, sigma, rho_t, rho_l, theta):
    """Multiple-scattering enhancement E(r) of apparent backscatter (Hogan 2006).

    Line-for-line port of the reference implementation in Hogan's multiscatter
    1.2.11 ``src/small_angle.c`` (``-algorithms original``), i.e. including the
    two refinements the paper only sketches:

    - **Eloranta's exact double scattering**: the once-scattered received energy
      is an exact per-origin integral, beta2/beta1 = fov_factor * (fs_od(R) -
      sum_i sigma_i exp(-rho_t^2 / (theta^2 dr_ij^2/R^2 + rho_l^2)) dr), not a
      pooled-Gaussian acceptance;
    - **within-gate multiple scattering**: forward scattering at gate i also
      feeds gate i itself with trapezoidal half weight (the C code's j-loop
      starts at j = i with trapez_int_factor = 0.5).

    Triple-and-higher orders keep the Gaussian moment parameterisation, fed by
    (once + multi) jointly, with the reference's exact statement order at j = i.

    r     : (n,) range grid [m], regular spacing
    sigma : (n,) extinction coefficient [m^-1] (cloud only; molecular negligible)
    rho_t : telescope half-angle field of view [rad] (top-hat)
    rho_l : laser 1/e divergence half-angle [rad]
    theta : forward-lobe 1/e half-width [rad] (lambda / (pi a_G))

    Returns E(r) >= 1 with E = 1 where no multiple scattering contributes.
    """
    n = r.size
    dr = r[1] - r[0]
    th2 = theta ** 2
    mu2_t = rho_l ** 2                     # transmitter Gaussian variance [rad^2]
    fov_factor = 1.0 / (1.0 - np.exp(-(rho_t ** 2) / mu2_t))

    # E-weighted moments of the once- and multiply-forward-scattered photons
    # (lateral second moment, pointing-angle second moment, covariance) --
    # identical bookkeeping to ms_E_once/ms_Emu2_once/ms_Pmu2_once/ms_EPcov_once.
    E1 = np.zeros(n); E1s2 = np.zeros(n); E1t2 = np.zeros(n); E1cv = np.zeros(n)
    Em = np.zeros(n); Ems2 = np.zeros(n); Emt2 = np.zeros(n); Emcv = np.zeros(n)
    D = np.zeros(n)                        # Eloranta accumulator (escaped part)

    # forward-scattering optical depth at each gate, trapezoidal including half
    # of the gate itself (the C code's fs_optical_depth default path)
    fs_od = np.cumsum(0.5 * (sigma + np.concatenate(([0.0], sigma[:-1]))) * dr)

    for i in np.nonzero(sigma > 0)[0]:
        ri = r[i]
        w_all = sigma[i] * dr

        # ---- j = i: within-gate contribution, trapezoid weight 0.5 ----------
        w = 0.5 * w_all
        D[i] += 0.5 * sigma[i] * np.exp(-(rho_t ** 2) / mu2_t)   # dr_ij = 0
        E1[i] += w
        E1s2[i] += w * mu2_t * ri ** 2
        E1t2[i] += w * (mu2_t + th2)
        E1cv[i] += w * mu2_t * ri
        # multi feed at j = i, using the JUST-updated once moments and the C
        # statement order (EPcov reads the already-updated Pmu2_multi at j = i)
        if E1[i] > 0.0:
            Ems2[i] += w * ((E1s2[i] + Ems2[i])
                            + 2.0 * (E1cv[i] + Emcv[i]) * 0.0)
            Emt2[i] += w * (E1t2[i] + Emt2[i] + th2 * (E1[i] + Em[i]))
            Emcv[i] += w * (E1cv[i] + Emcv[i])
            Em[i] += w * (E1[i] + Em[i])

        # frozen post-self-update gate-i values for the j > i transfers
        e1, e1s2, e1t2, e1cv = E1[i], E1s2[i], E1t2[i], E1cv[i]
        em, ems2, emt2, emcv = Em[i], Ems2[i], Emt2[i], Emcv[i]

        # ---- j > i: vectorised, trapezoid weight 1.0 -------------------------
        j = slice(i + 1, n)
        rj = r[j]
        dj = rj - ri
        D[j] += sigma[i] * np.exp(-(rho_t * rj) ** 2 / (th2 * dj ** 2 + mu2_t * rj ** 2))
        E1[j] += w_all
        E1s2[j] += w_all * (mu2_t * rj ** 2 + th2 * dj ** 2)
        E1t2[j] += w_all * (mu2_t + th2)
        E1cv[j] += w_all * (mu2_t * rj + th2 * dj)
        if e1 > 0.0:
            Ems2[j] += w_all * ((e1s2 + ems2)
                                + (e1t2 + emt2 + th2 * (e1 + em)) * dj ** 2
                                + 2.0 * (e1cv + emcv) * dj)
            Emt2[j] += w_all * (e1t2 + emt2 + th2 * (e1 + em))
            Emcv[j] += w_all * (e1cv + emcv
                                + (th2 * (e1 + em) + e1t2 + emt2) * dj)
            Em[j] += w_all * (e1 + em)

    # ---- received components (relative to single scattering) ----------------
    # exact Eloranta double scattering: fov_factor * (fs_od - D*dr)
    double_ratio = fov_factor * (fs_od - D * dr)
    # triple+ via the pooled-Gaussian top-hat integral (INT_TOPHAT in the C):
    # ratio = fov_factor * Em * (1 - exp(-(rho_t R)^2 * Em / Ems2))
    with np.errstate(divide="ignore", invalid="ignore"):
        fm = np.where(Em > 0,
                      1.0 - np.exp(-(rho_t * r) ** 2
                                   * Em / np.where(Ems2 > 0, Ems2, 1.0)),
                      0.0)
    return 1.0 + double_ratio + fov_factor * Em * fm


def eta_for_cloud_base(cbh, rho_t, rho_l, theta, alpha=30e-3, dr=2.0, tau_max=6.0):
    """eta for a deep homogeneous liquid cloud with base at range *cbh* [m].

    alpha : in-cloud extinction [m^-1] (default 30 km^-1, typical stratocumulus).
    The integral saturates by tau_max; eta = 1 / (2 * integral(sigma e^-2tau E dr)).
    """
    depth = tau_max / alpha
    r = np.arange(dr, cbh + depth + dr, dr)
    sigma = np.where(r >= cbh, alpha, 0.0)
    E = pvc_enhancement(r, sigma, rho_t, rho_l, theta)
    # Exact per-strip single-scatter integral (Hogan Eq. 9): each gate contributes
    # e^-2tau_prev (1 - e^-2 sigma dr)/2, so pure single scattering integrates to
    # exactly (1 - e^-2 tau_max)/2 regardless of dr (a plain left-edge Riemann sum
    # biases eta a few % HIGH, even above 1 at short range).
    tau_prev = np.concatenate(([0.0], np.cumsum(sigma[:-1]) * dr))
    strip = np.exp(-2 * tau_prev) * (1.0 - np.exp(-2 * sigma * dr)) / 2.0
    B_rel = np.sum(strip * E)                # = 1/2 for single scattering (E = 1)
    return 1.0 / (2.0 * B_rel)


def eloranta_double_check(rho_t, rho_l, theta, cbh=1000.0, alpha=30e-3):
    """Analytic double-scattering ratio (Hogan Eq. 11 = Eloranta 1998) vs the
    implementation's E - 1 on a VERY thin slab (tau = 0.06: triple+ orders are
    < 3 % of double, so E - 1 is essentially the double term) -- sanity check.
    (The authoritative validation is the digit-level match against Hogan's
    multiscatter 1.2.11 reference binary; see the report.)"""
    dr = 0.25
    depth = 2.0
    r = np.arange(dr, cbh + depth + dr, dr)
    sigma = np.where(r >= cbh, alpha, 0.0)
    R = r[-1]
    # analytic Eq. (11), midpoint integral
    fu = 1 - np.exp(-rho_t ** 2 / rho_l ** 2)
    s2 = rho_l ** 2 * R ** 2 + theta ** 2 * (R - r) ** 2
    fs = 1 - np.exp(-(rho_t * R) ** 2 / s2)
    beta2_over_beta1 = np.sum(sigma * fs) * dr / fu
    E = pvc_enhancement(r, sigma, rho_t, rho_l, theta)
    return beta2_over_beta1, E[-1] - 1.0


# --------------------------------------------------------------------------- config
WAVELENGTH = {"CL31": 910e-9, "CL51": 910e-9, "CL61": 910.55e-9,
              "CHM15k": 1064e-9, "Mini-MPL": 532e-9, "MPL": 532e-9}

OPTICS = {  # rho_t, rho_l [rad] (half-angles, 1/e)
    "CL31":     (0.83e-3, 0.57e-3),
    "CL51":     (0.56e-3, 0.21e-3),
    "CL61":     (0.56e-3, 0.28e-3),   # FOV +-0.56 mrad (CL61 User Guide) = same as CL51
    "CHM15k":   (0.23e-3, 0.15e-3),
    "Mini-MPL": (0.11e-3, 0.055e-3),
    "MPL":      (0.05e-3, 0.025e-3),
}

# range grid of the operational tables (km -> m), extended to 4 km
CBH_GRID_KM = np.array([0.25, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875,
                        2.125, 2.375, 2.75, 3.25, 3.75])


def theta_lobe(instr, a_um):
    return WAVELENGTH[instr] / (np.pi * a_um * 1e-6)


def eta_table(instr, a_um, alpha=30e-3, rho_t=None):
    rt, rl = OPTICS[instr]
    if rho_t is not None:
        rt = rho_t
    th = theta_lobe(instr, a_um)
    return np.array([eta_for_cloud_base(c * 1000.0, rt, rl, th, alpha=alpha)
                     for c in CBH_GRID_KM])


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path,
                    default=Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_multiple_scattering"))
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ---- 0. implementation cross-check vs the analytic double-scattering formula
    b2, pvc = eloranta_double_check(*OPTICS["CL31"], theta_lobe("CL31", 6.0))
    print(f"[check] thin-cloud double scattering: analytic={b2:.4f}  PVC(E_s)={pvc:.4f} "
          f"(agree within {100 * abs(b2 - pvc) / b2:.1f}%)")

    # ---- 1. validation: joint (droplet size, extinction) fit against the
    # operational CL31/CL51 tables. alpha matters because the MS lateral spread
    # grows with the penetration depth ~tau/alpha, so it shapes eta(range).
    tables = {"CL31": _ETA_CL31, "CL51": _ETA_CL51}
    a_scan = np.arange(3.0, 15.5, 1.0)
    al_scan = np.array([10e-3, 20e-3, 30e-3, 45e-3, 60e-3])
    best = (None, None, np.inf)
    for a in a_scan:
        for al in al_scan:
            tot = 0.0
            for k, tab in tables.items():
                model = np.interp(tab[:, 0], CBH_GRID_KM, eta_table(k, a, alpha=al))
                tot += np.sqrt(np.mean((model - tab[:, 1]) ** 2))
            if tot < best[2]:
                best = (a, al, tot)
    a_best, al_best, _ = best
    print(f"[fit] best a_G = {a_best:.1f} um (diameter {2 * a_best:.0f} um), "
          f"alpha = {al_best * 1e3:.0f} /km")
    for k, tab in tables.items():
        model = np.interp(tab[:, 0], CBH_GRID_KM, eta_table(k, a_best, alpha=al_best))
        print(f"      {k}: rms(model - operational table) = "
              f"{np.sqrt(np.mean((model - tab[:, 1]) ** 2)) * 100:.2f} % "
              f"(max dev {np.max(np.abs(model - tab[:, 1])) * 100:.2f} %)")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    for ax, (k, tab) in zip(axes, tables.items()):
        for a, ls in ((4.0, ":"), (a_best, "-"), (10.0, "--")):
            ax.plot(eta_table(k, a, alpha=al_best), CBH_GRID_KM, ls, color="C0",
                    label=f"PVC model, a_G={a:.0f} um")
        ax.plot(tab[:, 1], tab[:, 0], "ks", ms=6, label="operational table (Hopkin)")
        ax.set_title(k); ax.set_xlabel("eta"); ax.grid(alpha=0.3)
        ax.set_ylabel("cloud-base range [km]")
        ax.legend(fontsize=8)
    ax = axes[2]  # alpha sensitivity on CL51
    for al, c in ((10e-3, "C1"), (30e-3, "C0"), (60e-3, "C2")):
        ax.plot(eta_table("CL51", a_best, alpha=al), CBH_GRID_KM, "-", color=c,
                label=f"CL51, alpha={al * 1e3:.0f} /km")
    ax.plot(_ETA_CL51[:, 1], _ETA_CL51[:, 0], "ks", ms=6, label="operational table")
    ax.set_title("extinction sensitivity"); ax.set_xlabel("eta"); ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.suptitle("PVC multiple-scattering model vs operational CL31/CL51 eta tables "
                 "(validation + droplet-size fit)", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out / "fig_eta_validation_cl31_cl51.png", dpi=160)
    plt.close(fig)

    # ---- 2. eta(range) for every instrument at the fitted droplet size (+ band)
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.6))
    ax = axes[0]
    results = {}
    for k, c in (("CL31", "C0"), ("CL51", "C1"), ("CL61", "C3"),
                 ("CHM15k", "C2"), ("Mini-MPL", "C4"), ("MPL", "C5")):
        lo = eta_table(k, 10.0, alpha=al_best)
        hi = eta_table(k, 4.0, alpha=al_best)
        mid = eta_table(k, a_best, alpha=al_best)
        results[k] = mid
        ax.fill_betweenx(CBH_GRID_KM, lo, hi, color=c, alpha=0.15)
        ax.plot(mid, CBH_GRID_KM, "-", color=c, lw=2, label=k)
    ax.set_xlabel("eta (droplet radius 4-10 um band, line = fitted "
                  f"{a_best:.0f} um)")
    ax.set_ylabel("cloud-base range [km]"); ax.grid(alpha=0.3); ax.legend(fontsize=9)
    ax.set_title("Multiple-scattering factor per instrument (PVC, Hogan 2006)")

    ax = axes[1]  # CL61 vs CL51: identical receiver FOV (0.56 mrad), only the
    # divergence differs (0.28 vs 0.21 mrad) -> eta differs by < 0.01 everywhere,
    # which justifies applying the operational CL51 table to the CL61.
    ax.plot(eta_table("CL61", a_best, alpha=al_best), CBH_GRID_KM, "-", color="C3",
            label="CL61 (FOV 0.56, div 0.28 mrad)")
    ax.plot(eta_table("CL51", a_best, alpha=al_best), CBH_GRID_KM, "--", color="C1",
            label="CL51 optics (FOV 0.56, div 0.21 mrad)")
    ax.plot(_ETA_CL51[:, 1], _ETA_CL51[:, 0], "ks", ms=6,
            label="operational CL51 table (applied to CL61)")
    ax.set_xlabel("eta"); ax.grid(alpha=0.3); ax.legend(fontsize=9)
    ax.set_title("CL61 = CL51 receiver FOV (0.56 mrad) -> CL51 table justified")
    fig.suptitle("Recomputed multiple-scattering corrections "
                 "(910 / 1064 / 532 nm, stratocumulus alpha=30 /km)", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(out / "fig_eta_all_instruments.png", dpi=160)
    plt.close(fig)

    # ---- 3. python tables ready to paste into calibration/cloud/calibration.py
    with open(out / "eta_tables.py", "w", encoding="utf-8") as f:
        f.write(f"# PVC (Hogan 2006) eta tables, a_G={a_best:.1f} um, "
                f"alpha={al_best * 1e3:.0f}/km\n")
        for k, v in results.items():
            rows = ", ".join(f"[{c:.3f}, {e:.5f}]" for c, e in zip(CBH_GRID_KM, v))
            f.write(f"_ETA_{k.upper().replace('-', '')} = np.array([{rows}])\n")
    print(f"[out] figures + eta_tables.py -> {out}")

    for k, v in results.items():
        print(f"  {k:9s} eta(0.25 km)={v[0]:.3f}  eta(1.125 km)={v[4]:.3f}  "
              f"eta(2.375 km)={v[9]:.3f}")


if __name__ == "__main__":
    main()
