"""AERONET forward-Klett closure for the Payerne CHM15k (Jan-Mar 2025, the AERONET L2 overlap).

Purpose: anchor the validation chain absolutely. The CHM15k is the station reference the
CL61-cloud agrees with (-0.6 %); if its Rayleigh-calibrated AOD closes on AERONET, then
C_L(cloud) is right and the CL61-Rayleigh -12 % is confirmed against an absolute standard.
1064 nm -> no WV. C_L: daily Kalman (calib CSV) interpolated. LR: AERONET .lid when valid,
else 50 sr (sensitivity: forward AOD scales ~linearly with LR in thin conditions)."""
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
from validation.paper._cl61_aeronet_closure import parse_aeronet, fnum, forward_inversion  # reuse

WMO, IDENT, ALT, LAT, LON = "0-20000-0-06610", "A", 491.0, 46.813, 6.943
L1 = Path("D:/E-PROFILE_L1_2026")
STD = Path("calibration/data/standard_atmosphere_US_1976_50km.csv")
LEV20 = Path(r"C:/Users/hervo/Downloads/20260101_20261231_Payerne/20260101_20261231_Payerne.lev15")
LID = Path(r"C:/Users/hervo/Downloads/aeronet_2026/Payerne_2026_LID15.txt")
CAL = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/calib/0-20000-0-06610_A_rayleigh_L1.csv")
D0, D1 = datetime(2026, 3, 1), datetime(2026, 6, 20)

# daily Kalman C_L
rows = list(csv.DictReader(open(CAL, encoding="utf-8")))
kd = np.array([np.datetime64(r["time"][:10]) for r in rows if r["C_kalman"] not in ("", "nan")])
kv = np.array([float(r["C_kalman"]) for r in rows if r["C_kalman"] not in ("", "nan")])
print(f"CHM15k Kalman C_L: {kv.size} days, 2026 window median = "
      f"{np.median(kv[(kd >= np.datetime64('2026-03-01')) & (kd <= np.datetime64('2026-06-20'))]):.3e}")

def cl_for(day):
    return float(np.interp(np.datetime64(day).astype("datetime64[D]").astype(float),
                           kd.astype("datetime64[D]").astype(float), kv))

aod_rows = parse_aeronet(LEV20, "Date(")
aods = []
for r in aod_rows:
    if not (D0 <= r["_t"] <= D1):
        continue
    a1020 = fnum(r.get("AOD_1020nm")); ang = fnum(r.get("440-870_Angstrom_Exponent"))
    if np.isfinite(a1020) and np.isfinite(ang):
        aods.append((r["_t"], a1020 * (1020.0 / 1064.0) ** ang))
print(f"AERONET obs in window: {len(aods)}")

lid_rows = parse_aeronet(LID, "Site,")
lr_by_day = {}
for r in lid_rows:
    lr = fnum(r.get("Lidar_Ratio[1020nm]"))
    if np.isfinite(lr):
        lr_by_day.setdefault(r["_t"].date(), []).append(lr)
print(f"days with valid AERONET LR(1020): {len(lr_by_day)}")

_cache = {}
def mol_for(rng):
    if rng.size not in _cache:
        atm = load_standard_atmosphere(STD, rng + ALT)
        mol = calculate_molecular_properties(atm.temperature, atm.pressure, rng, 1064.0e-9)
        _cache[rng.size] = (mol.beta_mol * 1e6, mol.alpha_mol * 1e6, np.ones(rng.size))
    return _cache[rng.size]

res, day_cache = [], {}
for t_ae, aod_ae in aods:
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
                tt = np.array([base + timedelta(seconds=x) for x in (tv * 86400.0 if "day" in tu else tv)])
                rng = np.asarray(nc.variables["range"][:], "f8")
                rcs = np.asarray(nc.variables["rcs_0"][:], "f8")
                if rcs.shape[0] != tt.size:
                    rcs = rcs.T
            day_cache[ds] = (tt, rng, rcs)
    if day_cache[ds] is None:
        continue
    tt, rng, rcs = day_cache[ds]
    sel = np.abs((tt - t_ae) / timedelta(minutes=1)) <= 30
    if sel.sum() < 20:
        continue
    cl = cl_for(t_ae.date())
    prof = np.nanmean(rcs[sel], axis=0) / cl * 1e6     # beta_att Mm^-1 sr^-1
    bmol, amol, t2wv = mol_for(rng)
    if np.nanmax(prof[rng <= 6000]) > 30:
        continue
    lr = float(np.median(lr_by_day.get(t_ae.date(), [50.0])))
    lr = min(max(lr, 20.0), 100.0)
    r = forward_inversion(prof, rng, bmol, amol, t2wv, lr)
    if r is not None:
        res.append(dict(t=t_ae, aeronet=aod_ae, lidar=r[1], lr=lr))

print(f"closure samples: {len(res)}")
if res:
    ae = np.array([r["aeronet"] for r in res]); li = np.array([r["lidar"] for r in res])
    m = np.isfinite(ae) & np.isfinite(li) & (ae > 0.01)
    ratio = li[m] / ae[m]
    print(f"AERONET AOD1064 median={np.median(ae[m]):.3f}   CHM15k forward-Klett AOD median={np.median(li[m]):.3f}")
    print(f"ratio lidar/AERONET: median={np.median(ratio):.2f}  p25-p75={np.percentile(ratio,25):.2f}-{np.percentile(ratio,75):.2f}  n={m.sum()}")
    print("interpretation: ratio ~1 -> the CHM15k Rayleigh C_L (the network anchor) is absolutely validated;")
    print("scaling: an alternative constant k*C_L would scale the AOD by ~1/k in thin conditions.")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.plot(ae[m], li[m], "o", ms=5, color="#d62728", alpha=0.65)
    lim = (0, max(0.25, float(np.percentile(np.r_[ae[m], li[m]], 99)) * 1.15))
    ax.plot(lim, lim, "k--", lw=1)
    ax.set_xlim(*lim); ax.set_ylim(*lim); ax.grid(alpha=0.3)
    ax.set_xlabel("AERONET AOD @1064 nm"); ax.set_ylabel("CHM15k forward-Klett AOD @1064 nm")
    ax.set_title("Payerne CHM15k AERONET closure, 2026 (same window as the CL61 closure) — median ratio %.2f (n=%d)"
                 % (np.median(ratio), m.sum()))
    fig.tight_layout()
    out = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_chm15k_aeronet_closure.png")
    fig.savefig(out, dpi=150)
    print("saved", out)
print("CHM_CLOSURE_DONE")
