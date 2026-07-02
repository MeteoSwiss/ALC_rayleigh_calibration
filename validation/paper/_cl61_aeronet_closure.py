"""AERONET forward-Klett closure for the Payerne CL61: which lidar constant is right?

For every AERONET L2 observation (2026-03..06), take the CL61 signal +/-30 min, calibrate it with
  C_L(Rayleigh) = 1.251   and   C_L(cloud) = 1.425,
run an iterative FORWARD inversion at 910 nm (beta_tot = beta_att/T2, T2 from aerosol LR_aer +
molecular + CAMS water vapour), integrate the aerosol extinction to 6 km -> AOD_910, and compare
with the AERONET AOD interpolated to 910 nm (Angstrom 440-870). The constant whose AOD closes on
AERONET is the physically correct one. LR_aer: same-day AERONET .lid 1020 nm inversion when valid,
else 50 sr. (Method after extinction_ceilometer_v2_VPROFILE_AERONET.m / extinction_pay_v2.m,
adapted to 910 nm + WV and to a forward scheme driven by the calibrated attenuated backscatter.)
"""
import sys, csv, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties
from calibration.water_vapor_correction.water_vapor import (
    cams_water_vapor_profile, two_way_wv_transmission, DEFAULT_ABS_CROSS_SECTION)
from calibration.io.cams import ensure_cams_file

WMO, IDENT, ALT, LAT, LON = "0-20000-0-06610", "C", 491.0, 46.813, 6.943
L1 = Path("D:/E-PROFILE_L1_2026")
STD = Path("calibration/data/standard_atmosphere_US_1976_50km.csv")
LEV20 = Path(r"C:/Users/hervo/Downloads/20260101_20261231_Payerne/20260101_20261231_Payerne.lev15")
LID = Path(r"C:/Users/hervo/Downloads/aeronet_2026/Payerne_2026_LID15.txt")
CL_RAY, CL_CLD = 1.251, 1.425
LR_MOL = 8 * np.pi / 3
D0, D1 = datetime(2026, 3, 1), datetime(2026, 6, 20)
ZTOP = 6000.0


def parse_aeronet(path, header_marker, date_key="Date(dd:mm:yyyy)", time_key="Time(hh:mm:ss)"):
    """Generic AERONET text parser -> list of dicts keyed by header names. The header line is the
    first line CONTAINING header_marker (works for both the download portal and web-API formats)."""
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    hi = next(i for i, l in enumerate(lines) if header_marker in l and "Date(dd:mm:yyyy)" in l)
    hdr = [h.strip() for h in lines[hi].split(",")]
    for l in lines[hi + 1:]:
        p = l.split(",")
        if len(p) < len(hdr) - 2:
            continue
        d = dict(zip(hdr, p))
        try:
            d["_t"] = datetime.strptime(d[date_key] + " " + d[time_key], "%d:%m:%Y %H:%M:%S")
        except Exception:
            continue
        rows.append(d)
    return rows


def fnum(x):
    try:
        v = float(x)
        return np.nan if v <= -999 else v
    except Exception:
        return np.nan


# ---- AERONET AOD -> 910 nm ----
aod_rows = parse_aeronet(LEV20, "Date(")
aods = []
for r in aod_rows:
    if not (D0 <= r["_t"] <= D1):
        continue
    a1020 = fnum(r.get("AOD_1020nm")); ang = fnum(r.get("440-870_Angstrom_Exponent"))
    if np.isfinite(a1020) and np.isfinite(ang):
        aods.append((r["_t"], a1020 * (1020.0 / 910.0) ** ang, ang))
print(f"AERONET obs in window: {len(aods)}")

# ---- AERONET lidar ratio by day (1020 nm L2 inversions; sparse) ----
lid_rows = parse_aeronet(LID, "Site,")
lr_by_day = {}
for r in lid_rows:
    lr = fnum(r.get("Lidar_Ratio[1020nm]"))
    if np.isfinite(lr):        # L1.5 almucantar inversions (2026 has no L2 yet)
        lr_by_day.setdefault(r["_t"].date(), []).append(lr)
print(f"days with valid AERONET LR(1020): {len(lr_by_day)}")

# ---- static profiles (molecular at 910 nm on the CL61 grid; WV per month) ----
_grid_cache = {}
def profiles_for(rng, month_ds):
    key = (rng.size, month_ds[:6])
    if key in _grid_cache:
        return _grid_cache[key]
    atm = load_standard_atmosphere(STD, rng + ALT)
    mol = calculate_molecular_properties(atm.temperature, atm.pressure, rng, 910.0e-9)
    cams = ensure_cams_file(Path("D:/CAMS"), month_ds, auto_download=False)
    h, nwv = cams_water_vapor_profile(cams, LAT, LON,
                                      np.datetime64(f"{month_ds[:4]}-{month_ds[4:6]}-01"),
                                      np.datetime64(f"{month_ds[:4]}-{month_ds[4:6]}-28"))
    t2wv = two_way_wv_transmission(rng + ALT, ALT, h, nwv, DEFAULT_ABS_CROSS_SECTION, 910.74, 1.0)
    out = (mol.beta_mol * 1e6, mol.alpha_mol * 1e6, np.asarray(t2wv, "f8"))   # Mm^-1 (sr^-1)
    _grid_cache[key] = out
    return out


def forward_inversion(beta_att, rng, bmol, amol, t2wv, lr_aer, n_iter=8):
    """beta_att [Mm^-1 sr^-1] calibrated attenuated backscatter (incl. WV attenuation).
    Iterative forward: beta_tot = beta_att/(T2_scatt*T2_wv); alpha_aer = LR*(beta_tot-beta_mol).
    Returns (alpha_aer [Mm^-1], aod910) or None."""
    dz = np.median(np.diff(rng)) * 1e-6                 # Mm
    top = rng <= ZTOP
    alpha_aer = np.zeros_like(beta_att)
    for _ in range(n_iter):
        od = np.cumsum((alpha_aer + amol)) * dz         # one-way, to gate
        t2s = np.exp(-2 * (od - (alpha_aer + amol) * dz / 2))
        beta_tot = beta_att / (t2s * t2wv)
        ba = beta_tot - bmol
        ba = np.clip(ba, -0.05, None)                   # allow slight negative noise
        new = lr_aer * ba
        new[~np.isfinite(new)] = 0.0
        if np.nanmax(np.abs(new - alpha_aer)[top]) < 0.01:
            alpha_aer = new; break
        alpha_aer = new
    # extend the lowest reliable gate (~120 m, after overlap) down to the ground
    i0 = np.searchsorted(rng, 120.0)
    alpha_aer[:i0] = alpha_aer[i0]
    aod = float(np.nansum(alpha_aer[top]) * dz)
    if not np.isfinite(aod) or aod < 0 or aod > 2:
        return None
    return alpha_aer, aod


# ---- main loop over AERONET observations ----
res = []
day_cache = {}
for t_ae, aod910_ae, ang in aods:
    ds = t_ae.strftime("%Y%m%d")
    if ds not in day_cache:
        f = L1 / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
        if not f.exists():
            day_cache[ds] = None
        else:
            with Dataset(f) as nc:
                tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
                tv = np.asarray(nc.variables["time"][:], "f8")
                base = datetime(1970, 1, 1)
                tt = np.array([base + timedelta(seconds=x) for x in
                               (tv * 86400.0 if "day" in tu else tv)])
                rng = np.asarray(nc.variables["range"][:], "f8")
                rcs = np.asarray(nc.variables["rcs_0"][:], "f8")
                if rcs.shape[0] != tt.size:
                    rcs = rcs.T
            day_cache[ds] = (tt, rng, rcs)
    if day_cache[ds] is None:
        continue
    tt, rng, rcs = day_cache[ds]
    sel = np.abs((tt - t_ae) / timedelta(minutes=1)) <= 30
    if sel.sum() < 40:
        continue
    prof = np.nanmean(rcs[sel], axis=0) * 1e6            # vendor beta_att, Mm^-1 sr^-1 (C~1)
    bmol, amol, t2wv = profiles_for(rng, ds[:6] + "01")
    # cloud guard: any huge beta below 6 km -> skip
    if np.nanmax(prof[rng <= ZTOP]) > 30:
        continue
    lr = float(np.median(lr_by_day.get(t_ae.date(), [50.0])))
    lr = min(max(lr, 20.0), 100.0)
    row = dict(t=t_ae, aeronet=aod910_ae, lr=lr)
    # JOINT (C_L, baseline) closure: the dark probe found a NEGATIVE vendor baseline b in
    # beta_att; correct it (beta - b) BEFORE applying the constant. The (C_L, b) pair that
    # closes the AOD, cross-checked against the independent dark-probe b, is the answer.
    okall = True
    for tag, cl, b in (("ray", CL_RAY, 0.0), ("cld", CL_CLD, 0.0),
                       ("ray_b", CL_RAY, -0.03), ("cld_b1", CL_CLD, -0.02),
                       ("cld_b2", CL_CLD, -0.03), ("cld_b3", CL_CLD, -0.04)):
        r = forward_inversion((prof - b * cl) / cl, rng, bmol, amol, t2wv, lr)
        if r is None:
            okall = False; break
        row[tag] = r[1]
    if okall:
        res.append(row)

print(f"closure samples: {len(res)}")
if res:
    ae = np.array([r["aeronet"] for r in res])
    ray = np.array([r["ray"] for r in res]); cld = np.array([r["cld"] for r in res])
    lrs = np.array([r["lr"] for r in res])
    print(f"AERONET AOD910: median={np.median(ae):.3f}  (LR used: median={np.median(lrs):.0f} sr, "
          f"{np.sum(lrs != 50)} from AERONET inversions)")
    for tag, key in (("C_L=1.251 (Rayleigh), b=0    ", "ray"), ("C_L=1.425 (cloud),   b=0    ", "cld"),
                     ("C_L=1.251 (Rayleigh), b=-0.03", "ray_b"), ("C_L=1.425 (cloud),   b=-0.02", "cld_b1"),
                     ("C_L=1.425 (cloud),   b=-0.03", "cld_b2"), ("C_L=1.425 (cloud),   b=-0.04", "cld_b3")):
        x = np.array([r[key] for r in res])
        m = np.isfinite(x) & np.isfinite(ae) & (ae > 0.01)
        ratio = x[m] / ae[m]
        print(f"{tag}: AOD median={np.median(x[m]):.3f}  ratio lidar/AERONET: "
              f"median={np.median(ratio):.2f}  p25-p75={np.percentile(ratio,25):.2f}-{np.percentile(ratio,75):.2f}  n={m.sum()}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    for a, x, cl, nm, col in ((ax[0], ray, CL_RAY, "C$_L$(Rayleigh)=1.251", "#1f77b4"),
                              (ax[1], cld, CL_CLD, "C$_L$(cloud)=1.425", "#404040")):
        m = np.isfinite(x) & np.isfinite(ae)
        a.plot(ae[m], x[m], "o", ms=5, color=col, alpha=0.65)
        lim = (0, max(0.4, np.percentile(np.r_[ae[m], x[m]], 99) * 1.1))
        a.plot(lim, lim, "k--", lw=1)
        rat = np.median(x[m][ae[m] > 0.01] / ae[m][ae[m] > 0.01])
        a.set_xlim(*lim); a.set_ylim(*lim); a.grid(alpha=0.3)
        a.set_xlabel("AERONET AOD @910 nm"); a.set_ylabel("CL61 forward-Klett AOD @910 nm")
        a.set_title(f"{nm}:  median AOD ratio = {rat:.2f}  (n={m.sum()})")
    fig.suptitle("Payerne CL61 AERONET closure — which lidar constant closes the AOD?",
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_aeronet_closure.png")
    fig.savefig(out, dpi=150)
    print("saved", out)
print("AERONET_CLOSURE_DONE")
