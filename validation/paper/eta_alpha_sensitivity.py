"""eta_alpha_sensitivity.py - Appendix 4: sensitivity of the multiple-scattering factor eta.

Paper text: "for the CL31 at 1 km, eta = 0.84, 0.76 and 0.71 for alpha = 10, 20 and 30 km-1, whereas a change
of the radius by +-0.5 um changes eta by about 1 %", and Table A4-1 (eta at 0.5 / 1 / 2 km, droplet radius
5.5 um, alpha = 10 km-1; the operational tables of calibration/cloud/_filters.py).

Output on 2026-10-10:  alpha 10: 0.894 0.835 0.759 | alpha 20: 0.836 0.759 0.675 | alpha 30: 0.794 0.71 0.63
                       radius 5.0 / 5.5 / 6.0 um at 1 km: 0.845 / 0.835 / 0.825
Run:  python -m validation.paper.eta_alpha_sensitivity
History: inline checks of 2026-10-10, gathered here.
"""
import numpy as np

from validation.multiple_scattering_eta import eta_table, CBH_GRID_KM
from calibration.cloud import _filters as F

HEIGHTS = (0.5, 1.0, 2.0)


def at(z, e):
    return [round(float(np.interp(x, z, e)), 3) for x in HEIGHTS]


def main():
    z = np.asarray(CBH_GRID_KM, float)
    for a in (10e-3, 20e-3, 30e-3):
        print(f"alpha {a * 1e3:2.0f} /km  eta CL31 (0.5, 1, 2 km) =", at(z, eta_table("CL31", 5.5, alpha=a)))
    for r in (5.0, 5.5, 6.0):
        print(f"radius {r} um, alpha 10 /km: eta CL31 at 1 km = {at(z, eta_table('CL31', r, alpha=10e-3))[1]}")
    print("Table A4-1 (operational tables, 0.5 / 1 / 2 km):")
    for name in ("CL31", "CL51", "CL61", "CHM15K", "MINIMPL"):
        t = getattr(F, f"_ETA_{name}", None)
        if t is None:
            print(f"  {name}: no operational table")
            continue
        t = np.asarray(t, float)
        print(f"  {name}:", at(t[:, 0], t[:, 1]))


if __name__ == "__main__":
    main()
