"""Temporal variability of the CL61 offset WITHIN the 25.5-h hood window
(2026-05-26 11:45 -> 05-27 13:15): the direct night-time dark measurement.
Hourly band means of beta_att; Payerne late May: sun below horizon ~19:15-03:50 UT."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from validation.paper._cl61_dark_windows import load

t1, t2 = datetime(2026, 5, 26, 11, 45), datetime(2026, 5, 27, 13, 15)
X, rng = load(t1, t2)
# rebuild profile times (load() concatenates; reload times the same way)
from netCDF4 import Dataset
from pathlib import Path
tt = []
for ds in ("20260526", "20260527"):
    f = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610") / ds[:4] / ds[4:6] / f"L1_0-20000-0-06610_C{ds}.nc"
    with Dataset(f) as nc:
        tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
        tv = np.asarray(nc.variables["time"][:], "f8")
        from datetime import timedelta
        t = [datetime(1970, 1, 1) + timedelta(seconds=x) for x in (tv * 86400.0 if "day" in tu else tv)]
        tt += [x for x in t if t1 <= x <= t2]
tt = np.array(tt)
assert len(tt) == X.shape[0], (len(tt), X.shape)

BANDS = [(3000, 6000), (8000, 12000)]
hours = np.array([t.replace(minute=0, second=0, microsecond=0) for t in tt])
uh = sorted(set(hours))
series = {b: [] for b in BANDS}
for h in uh:
    sel = hours == h
    for b in BANDS:
        zb = (rng >= b[0]) & (rng <= b[1])
        series[b].append(float(np.nanmean(X[sel][:, zb]) * 1e6))
night = np.array([(h.hour >= 20 or h.hour < 4) for h in uh])
print("band        day mean   night mean   night/day")
for b in BANDS:
    v = np.array(series[b])
    d, n = np.nanmedian(v[~night]), np.nanmedian(v[night])
    print(f"{b[0]/1000:.0f}-{b[1]/1000:.0f} km   {d:+8.4f}   {n:+8.4f}     {n/d:+5.2f}" if d != 0 else "")
fig, ax = plt.subplots(figsize=(13, 5))
for b, c in zip(BANDS, ("#1f77b4", "#d62728")):
    ax.plot(uh, series[b], "o-", color=c, lw=1.4, label=f"{b[0]/1000:.0f}-{b[1]/1000:.0f} km")
ax.axvspan(datetime(2026, 5, 26, 19, 15), datetime(2026, 5, 27, 3, 50), color="0.85", label="night")
ax.axhline(0, color="k", lw=0.8)
ax.set_ylabel(r"hourly-mean $\beta_{att}$ (hood on) [Mm$^{-1}$ sr$^{-1}$]")
ax.grid(alpha=0.3); ax.legend(fontsize=9)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %Hh"))
ax.set_title("CL61 offset during the 25.5-h hood test (2026-05-26/27) - direct day/night comparison")
fig.tight_layout()
fig.savefig("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_hood_night.png", dpi=150)
print("saved fig_cl61_hood_night.png")
