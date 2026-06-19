#!/usr/bin/env python3
"""
Reproduce the MeteoFrance "Klett method implementation report" examples on the
Toulouse miniMPL (reconstructed L1 from L2-monthly).

Implements the THREE code versions the report compares, as standalone Klett variants:
  baseline           : original E-PROFILE code  (exp(-2..), truncated denominator)
  sign               : report sign correction   (exp(+2..), truncated denominator)
  sign+integration   : full fix (= repo HEAD)    (exp(+2..), full-window denominator)

Outputs (C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/minimpl/report/):
  ext_profiles_<date>.png     extinction vs altitude for several lidar ratios, 3 versions
  aod_table.csv               Klett AOD baseline vs corrected vs AERONET, the 4 report nights

Run with the withfix worktree importable (it provides the preprocessing helpers; only
the Klett step differs between versions and is implemented here).
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORKTREE = Path("D:/rc_worktrees/withfix")
sys.path.insert(0, str(WORKTREE))

from calibration import CalibrationOptions, InstrumentInfo
from calibration.config import InstrumentType
from calibration.io.data_loader import (
    build_file_paths, load_l1_data, filter_time_range, filter_cloudy_profiles)
from calibration.rayleigh.atmosphere import (
    load_standard_atmosphere, calculate_molecular_properties, MOLECULAR_LIDAR_RATIO)
from calibration.rayleigh.rayleigh_fit import find_optimal_molecular_window

sys.path.insert(0, str(Path(__file__).resolve().parent))
from get_aeronet import nightly_aod_532

MINIMPL_ROOT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/_minimpl_L1")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/minimpl/report")
STD_ATM = WORKTREE / "standard_atmosphere_US_1976_50km.csv"
REPORT_NIGHTS = ["20221010", "20230305", "20230929", "20250114"]
REPORT_PHOTO = {"20221010": 0.153, "20230305": 0.423, "20230929": 0.091, "20250114": 0.091}
LRS = [30, 40, 50, 60, 70]


def _info():
    return InstrumentInfo(site_name="MINI_MPL Toulouse", wmo_id="0-20000-0-07617",
                          identifier="A", instrument_type=InstrumentType.MINI_MPL,
                          latitude=43.578, longitude=1.374, altitude=154.64)


def _options():
    o = CalibrationOptions.from_json(WORKTREE / "options.json")
    o.folder_root = MINIMPL_ROOT
    return o


# --------------------------------------------------------------------------- #
# Three Klett variants. Each returns signed beta_aer (NO zeroing) so the report's
# "before zero-forcing" extinction can be shown.
# --------------------------------------------------------------------------- #
def _klett(beta_att, beta_mol, range_alc, ref_idx, lr, ref_val, i_start, i_end, version):
    n = len(beta_att)
    dz = abs(range_alc[1] - range_alc[0])
    lr_diff = lr - MOLECULAR_LIDAR_RATIO
    cum = np.zeros(n + 1)
    np.cumsum(beta_mol * lr_diff * dz, out=cum[1:])
    R = np.arange(i_start, i_end)
    qt = cum[ref_idx] - cum[R]
    sign = -2.0 if version == "baseline" else 2.0
    numerator = beta_att[R] * np.exp(sign * qt)

    if version == "sign+integration":
        # full-window forward integral (repo HEAD)
        T = 2.0 * (cum[ref_idx] - cum[R])
        weighted = beta_att[R] * np.exp(T) * lr * dz
        cum_den = np.zeros(n + 1)
        np.cumsum(weighted, out=cum_den[i_start + 1:i_end + 1])
        denom_int = cum_den[ref_idx] - cum_den[R]
        denominator = ref_val + 2 * denom_int
    else:
        # truncated denominator (baseline / sign-only), with sign per version
        T = sign * (cum[ref_idx] - cum[i_start:ref_idx])
        weighted = beta_att[i_start:ref_idx] * np.exp(T) * lr * dz
        rev = np.cumsum(weighted[::-1])[::-1]
        denom_sum = np.zeros(len(R))
        off = R - i_start
        valid = off < len(rev)
        denom_sum[valid] = rev[off[valid]]
        denominator = ref_val + 2 * denom_sum

    beta_aer = np.full(n, np.nan)
    beta_aer[R] = numerator / denominator - beta_mol[R]
    return beta_aer


def preprocess(date_str):
    """Replicate the pipeline up to the Rayleigh window; return night arrays."""
    info, options = _info(), _options()
    files = [f for f in build_file_paths(date_str, info, options) if f.exists()]
    if not files:
        return None
    data = load_l1_data(files, info.instrument_type)
    if data is None:
        return None
    data = filter_time_range(data, date_str, options)
    if len(data.time) == 0:
        return None
    data, clear, _ = filter_cloudy_profiles(data, options, info.instrument_type.no_cloud_value)
    if not clear:
        return None
    atm = load_standard_atmosphere(STD_ATM, np.asarray(data.altitude_grid))
    mol = calculate_molecular_properties(atm.temperature, atm.pressure, data.range_alc,
                                         info.instrument_type.wavelength_nm * 1e-9)
    rcs_mean = np.nanmean(data.rcs, axis=0)
    if np.all(np.isnan(rcs_mean)):
        return None
    signal = rcs_mean / data.range_alc ** 2
    fit = find_optimal_molecular_window(
        signal=signal, p_mol=mol.p_mol, range_alc=data.range_alc,
        half_length_options_m=options.half_length_options_m,
        range_start_m=options.range_start_m, range_end_m=options.range_end_m,
        increment_bins=options.fit_range_increment_bins)
    return dict(range=data.range_alc, altitude_grid=np.asarray(data.altitude_grid),
                rcs_mean=rcs_mean, beta_mol=mol.beta_mol, fit=fit)


def night_beta_att(pre, lr):
    """beta_att, reference setup for the nominal molecular window."""
    rng, rcs, beta_mol, fit = pre["range"], pre["rcs_mean"], pre["beta_mol"], pre["fit"]
    beta_att = (rcs / rng ** 2) / fit.slope * rng ** 2
    mask = (rng >= fit.range_start_m) & (rng <= fit.range_end_m) & ~np.isnan(rcs)
    idx = np.where(mask)[0]
    ref_idx = int((idx[0] + idx[-1]) / 2)
    ref_val = np.nanmean(beta_att[mask] / beta_mol[mask])
    return beta_att, ref_idx, ref_val, idx[-1]


def fig_ext_profiles(date_str):
    pre = preprocess(date_str)
    if pre is None:
        print(f"  {date_str}: night not usable for profile figure"); return
    rng, beta_mol, alt = pre["range"], pre["beta_mol"], pre["altitude_grid"] / 1000.0
    versions = ["baseline", "sign", "sign+integration"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
    for ax, ver in zip(axes, versions):
        for lr in LRS:
            beta_att, ref_idx, ref_val, i_end = night_beta_att(pre, lr)
            beta_aer = _klett(beta_att, beta_mol, rng, ref_idx, lr, ref_val, 0, i_end, ver)
            ext = beta_aer * lr  # signed, before zeroing
            ax.plot(ext * 1e6, alt, lw=1.0, label=f"S={lr} sr")
        ax.axvline(0, color="k", lw=0.6)
        ax.set_title(ver)
        ax.set_xlabel("aerosol extinction  [Mm$^{-1}$]")
        ax.grid(alpha=0.3)
        ax.set_ylim(0, min(8, alt.max()))
    axes[0].set_ylabel("altitude ASL [km]")
    axes[0].legend(fontsize=8)
    fig.suptitle(f"Toulouse miniMPL {date_str} — aerosol extinction (before zeroing), "
                 f"3 code versions (cf. report Fig. 8)")
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / f"ext_profiles_{date_str}.png"
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {p}")


def aod_table():
    rows = []
    for d in REPORT_NIGHTS:
        pre = preprocess(d)
        rec = {"date": d, "AOD_baseline": None, "AOD_corrected": None,
               "AOD_aeronet_mine": None, "AOD_photometer_report": REPORT_PHOTO[d]}
        if pre is not None:
            rng, beta_mol = pre["range"], pre["beta_mol"]
            for ver, key in [("baseline", "AOD_baseline"), ("sign+integration", "AOD_corrected")]:
                beta_att, ref_idx, ref_val, i_end = night_beta_att(pre, 50)
                beta_aer = _klett(beta_att, beta_mol, rng, ref_idx, 50, ref_val, 0, i_end, ver)
                ext = np.maximum(0, beta_aer) * 50  # zeroed, as report does for AOD
                rec[key] = float(np.trapz(np.nan_to_num(ext), rng))
        try:
            rec["AOD_aeronet_mine"] = nightly_aod_532("Toulouse_MF", d)
        except Exception as exc:
            print(f"   AERONET fetch failed for {d}: {exc}")
        rows.append(rec)
        print(f"  {d}: baseline={rec['AOD_baseline']}, corrected={rec['AOD_corrected']}, "
              f"AERONET(mine)={rec['AOD_aeronet_mine']}, report_photo={rec['AOD_photometer_report']}")
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "aod_table.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"  saved {OUT/'aod_table.csv'}")


if __name__ == "__main__":
    print("== extinction profiles (centerpiece: 2025-01-14) ==")
    for d in REPORT_NIGHTS:
        fig_ext_profiles(d)
    print("== AOD table ==")
    aod_table()
