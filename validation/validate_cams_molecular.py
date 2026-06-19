"""validate_cams_molecular.py

Validate the new CAMS-molecular option (molecular_source='cams') against:
  (a) the US Standard 1976 atmosphere (the previous default, identical to the
      MATLAB reference which always uses the standard atmosphere); and
  (b) the MATLAB molecular formula calcMolecularProperties.m, by saving the CAMS
      T/p + Python beta_mol to a .mat for an identical-input cross-check.

For the Payerne CL61 (910.74 nm) it builds the molecular reference both ways on a
30 m grid, computes beta_mol with calculate_molecular_properties, and reports the
window-mean CAMS/std ratio (the leading-order effect on the lidar constant).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.io import savemat

# The package is pip-installed and ships the US Standard Atmosphere as package data,
# so no sys.path / cwd juggling is needed here.
from calibration.rayleigh.atmosphere import (
    load_standard_atmosphere,
    load_cams_atmosphere,
    calculate_molecular_properties,
)

# Payerne CL61 (0-20000-0-06610_C)
LAT, LON, ALT = 46.813, 6.944, 491.0
WAVELENGTH_M = 910.74e-9
STD_ATM = None  # use the US Standard Atmosphere shipped as package data
CAMS_DIR = Path("D:/CAMS")

# 30 m E-PROFILE grid; molecular fit window 2-6 km AGL
range_alc = np.arange(0.0, 12000.0, 30.0)
altitude_grid = ALT + range_alc
win = (range_alc >= 2000) & (range_alc <= 6000)

atm_std = load_standard_atmosphere(STD_ATM, altitude_grid)
mol_std = calculate_molecular_properties(
    atm_std.temperature, atm_std.pressure, range_alc, WAVELENGTH_M
)

print("Payerne CL61 molecular reference: CAMS T/p vs US Standard 1976")
print(f"{'month':7s} {'bmol(CAMS)/std':>15s} {'impliedCL(CAMS)/std':>20s} "
      f"{'T_std-T_cams@4km':>17s} {'P_cams/P_std@4km':>17s}")

i4 = int(np.argmin(np.abs(range_alc - 4000)))
cl_ratios = []
saved = False
for ym, day in [("202602", 15), ("202603", 15), ("202604", 15), ("202605", 15)]:
    cams_file = CAMS_DIR / f"CAMS_Beta_{ym}.nc"
    if not cams_file.exists():
        print(f"{ym:7s}  CAMS file missing")
        continue
    night = np.datetime64(f"{ym[:4]}-{ym[4:6]}-{day:02d}")
    t_start = night - np.timedelta64(4, "h")
    t_end = night + np.timedelta64(4, "h")
    atm_cams = load_cams_atmosphere(cams_file, LAT, LON, t_start, t_end, altitude_grid)
    if atm_cams is None:
        print(f"{ym:7s}  CAMS profile unavailable")
        continue
    mol_cams = calculate_molecular_properties(
        atm_cams.temperature, atm_cams.pressure, range_alc, WAVELENGTH_M
    )
    ratio = float(np.nanmean(mol_cams.beta_mol[win] / mol_std.beta_mol[win]))
    # The Rayleigh fit scales signal/range^2 onto p_mol (~ beta_mol), so the lidar
    # constant scales ~ 1/ratio to first order.
    cl_ratio = 1.0 / ratio
    cl_ratios.append(cl_ratio)
    dT = float(atm_std.temperature[i4] - atm_cams.temperature[i4])
    pr = float(atm_cams.pressure[i4] / atm_std.pressure[i4])
    print(f"{ym:7s} {ratio:15.4f} {cl_ratio:20.4f} {dT:+17.1f} {pr:17.4f}")
    if not saved:
        savemat("C:/DATA/Projects/202606_E-PROFILE_calibration/validate_cams_molecular.mat", {
            "range_alc": range_alc,
            "altitude_grid": altitude_grid,
            "wavelength_m": WAVELENGTH_M,
            "T_cams": atm_cams.temperature,
            "P_cams": atm_cams.pressure,
            "T_std": atm_std.temperature,
            "P_std": atm_std.pressure,
            "beta_mol_cams_py": mol_cams.beta_mol,
            "beta_mol_std_py": mol_std.beta_mol,
            "month": ym,
        })
        saved = True

if cl_ratios:
    m = float(np.mean(cl_ratios))
    print(f"\nmean implied CL(CAMS)/std = {m:.4f}  ({(m - 1) * 100:+.2f} %)")
    print("Saved validate_cams_molecular.mat for the MATLAB cross-check.")
print("VALIDATE_DONE")
