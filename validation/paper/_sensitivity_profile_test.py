"""
_sensitivity_profile_test.py — PROTOTYPE (Payerne; standalone). Instead of censoring gates with an SNR
gate, show a **detection-limit profile per instrument** so the reader sees up to what altitude each
ceilometer is usable — the CL31 question. Built like the operational sensitivity product
(`calibration.sensitivity.detection`):

    beta_att_min(r) = SNR * sigma(r, tau),   sigma(r, tau) = sigma0(r) * sqrt(dt / tau)

with sigma0(r) the per-gate NATIVE-sampling noise of the calibrated beta (robust lag-1 estimate on the
native profiles, so the slow atmospheric signal cancels), SNR = 3, tau = the validation averaging
(30 min). The **maximum usable altitude** is the highest range at which the detection limit still sits
below the clear-air molecular signal beta_mol(r) (from CAMS) — above it the instrument cannot even see
the Rayleigh floor, so its beta there is noise, not a measurement.

Usage:  python -m validation.paper._sensitivity_profile_test
Output: fig_sensitivity_profile_payerne.png + printed max-usable-altitude table.
"""
from __future__ import annotations
import os, sys, glob, warnings
from datetime import datetime, timedelta
from pathlib import Path
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from validation.paper import intercompare as IC
from validation.paper.calib_benchmark import key_of
from calibration.sensitivity.detection import scale_sigma_to_tau

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
L1_ROOT = Path("D:/E-PROFILE_L1_2026")
SNR, TAU = 3.0, 1800.0                       # detection SNR and averaging time (30 min, = validation)
ZMAX = 6000.0
# Payerne channels: (wmo, ident, itype, calib, laser_nm, colour). rcs_0/C_L*1e6 = calibrated beta_att at
# the instrument's OWN wavelength -> the molecular floor must be evaluated at that wavelength too
# (β_mol at 910 nm is ~1.877× the 1064 nm value).
CH = [("0-20000-0-06610", "A", "CHM15k", "rayleigh", 1064.47, "#d62728"),
      ("0-20000-0-06610", "B", "CL31", "cloud", 909.7, "#ff7f0e"),
      ("0-20000-0-06610", "C", "CL61", "cloud", 910.74, "#1f77b4")]
LAT, LON, SALT = 46.8137, 6.9425, 491.0


def _robust_std(x, axis=0):
    med = np.nanmedian(x, axis=axis)
    return 1.4826 * np.nanmedian(np.abs(x - np.expand_dims(med, axis)), axis=axis)


def native_noise(wmo, ident, days):
    """sigma0(r) [rcs units] from the robust lag-1 profile-to-profile scatter of the native rcs_0
    (slow signal cancels in the difference), median over sample days; + median rcs_0(r), dt, range."""
    sig, sigcols, rcsmed, rng, dt = [], [], [], None, None
    for d in days:
        f = L1_ROOT / wmo / f"{d.year}" / f"{d.month:02d}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"
        if not f.exists():
            continue
        try:
            with Dataset(f) as nc:
                t = np.asarray(nc.variables["time"][:], "f8")
                r = np.asarray(nc.variables["range"][:], "f8")
                rcs = np.asarray(nc.variables["rcs_0"][:], "f8")
                if rcs.shape != (t.size, r.size):
                    rcs = rcs.T
        except Exception:
            continue
        rcs[~np.isfinite(rcs)] = np.nan
        if rng is None:
            rng = r
            dt = float(np.median(np.diff(t)) * 86400.0)      # native cadence [s]
        if r.size != rng.size:
            continue
        drcs = np.diff(rcs, axis=0)                          # lag-1: atmospheric signal ~cancels
        sigcols.append(_robust_std(drcs, axis=0) / np.sqrt(2.0))
        rcsmed.append(np.nanmedian(rcs, axis=0))
    if not sigcols:
        return None
    return (np.nanmedian(np.vstack(sigcols), axis=0), np.nanmedian(np.vstack(rcsmed), axis=0), rng, dt)


def main():
    warnings.filterwarnings("ignore")
    days = [datetime(2026, 6, 1) + timedelta(days=k) for k in range(0, 30, 3)]   # 10 June sample days
    cams = IC.find_cams_month(2026, 6)
    fig, ax = plt.subplots(figsize=(10.5, 8.4))
    print(f"Payerne detection-limit prototype — SNR={SNR:.0f}, tau={TAU:.0f}s, {len(days)} sample days\n")
    print(f"{'instrument':>10} | native dt | max usable altitude (β_att_min < β_mol at own λ)")
    ann = []
    for wmo, ident, itype, calib, laser_nm, col in CH:
        res = native_noise(wmo, ident, days)
        if res is None:
            print(f"  {itype}: no L1"); continue
        sigma0_rcs, rcsmed, rng, dt = res
        alt_asl = rng + SALT
        # calibrate the noise + signal to beta_att (Mm^-1 sr^-1) with the daily-median C_L
        ch = dict(wmo=wmo, ident=ident, itype=itype, calib=calib)
        cal = IC.load_calib_series(key_of(ch), "L1")
        C_L = float(np.nanmedian(cal[1])) if cal is not None else 1.0
        sigma0_b = sigma0_rcs / C_L * 1e6                              # native-sampling noise on beta
        beta_min = SNR * scale_sigma_to_tau(sigma0_b, dt, TAU)        # detection limit at tau [Mm^-1 sr^-1]
        beta_sig = np.abs(rcsmed) / C_L * 1e6                         # median measured beta_att
        # clear-air molecular floor at the INSTRUMENT'S OWN wavelength (the weakest signal it must see)
        bmol = IC._mol_att_cams(cams, LAT, LON, np.datetime64("2026-05-31"), np.datetime64("2026-07-10"),
                                alt_asl, SALT, laser_nm) if cams is not None else np.full(rng.size, np.nan)
        z = rng / 1000.0
        m = (rng <= ZMAX) & np.isfinite(beta_min)
        # max usable altitude = highest gate (in-overlap, >200 m) where beta_min < the molecular floor
        usable = m & (rng > 200) & np.isfinite(bmol) & (beta_min < bmol)
        zmax_use = float(np.nanmax(rng[usable])) if usable.any() else np.nan
        ax.plot(beta_min[m], z[m], "-", color=col, lw=2.4, label=f"{itype} ({laser_nm:.0f} nm): β_att_min")
        # median measured signal only where the instrument is usable (above -> noise, not shown)
        mv = m & np.isfinite(bmol) & (rng <= (zmax_use if np.isfinite(zmax_use) else 0))
        ax.plot(beta_sig[mv], z[mv], ":", color=col, lw=1.3, alpha=0.75)
        if np.isfinite(zmax_use):
            ax.axhline(zmax_use / 1000.0, color=col, ls="--", lw=1.1, alpha=0.65)
            ann.append((zmax_use / 1000.0, col, f"{itype}: {zmax_use/1000:.2f} km"))
        print(f"  {itype:>8} | {dt:5.0f} s   | {zmax_use/1000:.2f} km" if np.isfinite(zmax_use)
              else f"  {itype:>8} | {dt:5.0f} s   | (below floor everywhere)")
    # molecular floors: the target each instrument must beat — 1064 (CHM) and 910 (CL31/CL61)
    zf = np.arange(0, ZMAX, 30.0)
    for lam, cc, ls, lab in ((1064.47, "k", "-", "β_mol 1064 nm (CHM floor)"),
                             (910.0, "#555", "--", "β_mol 910 nm (CL31/CL61 floor)")):
        bf = IC._mol_att_cams(cams, LAT, LON, np.datetime64("2026-05-31"), np.datetime64("2026-07-10"),
                              zf + SALT, SALT, lam) if cams is not None else None
        if bf is not None:
            ax.plot(bf, zf / 1000.0, ls, color=cc, lw=1.7, label=lab)
    for zk, col, txt in ann:                                          # usable-altitude callouts
        ax.text(3.4, zk, txt, color=col, fontsize=9.5, fontweight="bold", va="center", ha="right")
    ax.set_xscale("log"); ax.set_xlim(1e-3, 5)
    ax.set_xlabel(r"attenuated backscatter  [Mm$^{-1}$sr$^{-1}$]")
    ax.set_ylabel("altitude AGL  [km]"); ax.set_ylim(0, ZMAX / 1000.0)
    ax.set_title("Detection-limit profile — up to what altitude is each ceilometer usable? (Payerne, June 2026)",
                 fontsize=12, fontweight="bold")
    ax.grid(alpha=0.3, which="both")
    ax.text(0.015, 0.02, "solid = β_att_min (noise floor, SNR 3 @ 30 min)\ndotted = median measured β "
            "(shown only where usable)\nblack/grey = clear-air molecular floor\ndashed line = max usable "
            "altitude (β_att_min = β_mol)", transform=ax.transAxes, fontsize=8.5, va="bottom",
            bbox=dict(boxstyle="round", fc="#f5f5f5", ec="#999"))
    ax.legend(fontsize=8.8, loc="upper right")
    fig.tight_layout()
    outp = OUT / "fig_sensitivity_profile_payerne.png"
    fig.savefig(outp, dpi=170); plt.close(fig)
    print("\n->", outp)


if __name__ == "__main__":
    main()
