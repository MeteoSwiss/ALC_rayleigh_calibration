# -*- coding: utf-8 -*-
"""Solar-leak test using the 9 Jun 2026 hood session: two black bin-bags were added over the
hood 10:34-11:39 UTC to improve solar-background hermeticity. If part of the daytime hood signal
is solar leak (not electronic offset), the FAR-range non-range-corrected level should step DOWN
while the bags are on. We plot the far-range band level vs time for all three instruments, bag
window shaded. Sun keeps rising through the morning, so a downward step *against* that trend is
unambiguous solar leak."""
import sys, warnings, math
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
LAT, LON = 46.8137, 6.9425
INST = {"CHM15k": ("A", "#d62728", (8000, 14000)),
        "CL31":   ("B", "#9467bd", (5000, 7000)),
        "CL61":   ("C", "#1f77b4", (8000, 14000))}
BAG1, BAG2 = datetime(2026, 6, 9, 10, 34), datetime(2026, 6, 9, 11, 39)
T0, T1 = datetime(2026, 6, 9, 9, 0), datetime(2026, 6, 9, 12, 15)

def solar_elev(dt):
    # NOAA low-precision solar elevation (deg)
    d = (dt - datetime(2026, 1, 1)).total_seconds() / 86400.0
    g = math.radians((357.529 + 0.98560028 * (d + 0.5)) % 360)
    L = math.radians((280.459 + 0.98564736 * (d + 0.5)) % 360)
    lam = L + math.radians(1.915 * math.sin(g) + 0.020 * math.sin(2 * g))
    eps = math.radians(23.439)
    dec = math.asin(math.sin(eps) * math.sin(lam))
    hutc = dt.hour + dt.minute / 60 + dt.second / 3600
    # equation of time (min) ~ small; use mean sun
    ha = math.radians(15 * (hutc - 12) + LON)
    la = math.radians(LAT)
    el = math.asin(math.sin(la) * math.sin(dec) + math.cos(la) * math.cos(dec) * math.cos(ha))
    return math.degrees(el)

fig, ax = plt.subplots(1, 3, figsize=(17, 5.6), sharex=True)
print("instrument  far band      median(no-bag) median(bag)  step   verdict")
for k, (inst, (ident, col, (lo, hi))) in enumerate(INST.items()):
    f = L1 / "2026" / "06" / f"L1_0-20000-0-06610_{ident}20260609.nc"
    with Dataset(f) as nc:
        tv = np.asarray(nc.variables["time"][:], "f8")
        tt = np.array([datetime(1970, 1, 1) + timedelta(days=x) for x in tv])
        r = np.asarray(nc.variables["range"][:], "f8")
        x = np.asarray(nc.variables["rcs_0"][:], "f8")
        if x.shape[0] != tt.size:
            x = x.T
    zkm = r / 1000.0
    band = (r >= lo) & (r <= hi)
    far = np.nanmean(x[:, band] / zkm[None, band] ** 2, axis=1)     # per-profile far-range level
    sel = (tt >= T0) & (tt <= T1)
    tt, far = tt[sel], far[sel]
    # 5-min median series
    tb, mb = [], []
    b0 = T0
    while b0 < T1:
        m = (tt >= b0) & (tt < b0 + timedelta(minutes=5))
        if m.sum() > 3:
            tb.append(b0 + timedelta(minutes=2.5)); mb.append(np.nanmedian(far[m]))
        b0 += timedelta(minutes=5)
    tb = np.array(tb); mb = np.array(mb)
    hrs = np.array([t.hour + t.minute / 60 for t in tb])
    a = ax[k]
    a.axvspan(10 + 34 / 60, 11 + 39 / 60, color="0.5", alpha=0.18, label="bags on")
    a.plot(hrs, mb, "-o", ms=3, color=col)
    a.axhline(0, color="0.7", lw=0.7)
    a.set_title(f"{inst}  ({lo/1000:.0f}-{hi/1000:.0f} km)")
    a.set_xlabel("time UTC [h]"); a.grid(alpha=0.3)
    if k == 0:
        a.set_ylabel("far-range mean  rcs_0/z²  (native units)")
    a.legend(fontsize=8, loc="upper left")
    inbag = np.array([BAG1 <= t <= BAG2 for t in tb])
    nobag = ~inbag & (hrs > 9.3)
    mnb, mib = np.nanmedian(mb[nobag]), np.nanmedian(mb[inbag])
    rel = 100 * (mib - mnb) / abs(mnb) if mnb else float("nan")
    verdict = "solar leak (down)" if mib < mnb else "no drop"
    print(f"{inst:9s}  {lo/1000:.0f}-{hi/1000:.0f} km    {mnb:+12.4g} {mib:+12.4g}  {rel:+6.0f}%  {verdict}")
print(f"\nsolar elevation: {solar_elev(datetime(2026,6,9,9,30)):.0f}° (09:30) -> "
      f"{solar_elev(datetime(2026,6,9,11,0)):.0f}° (11:00) -> {solar_elev(datetime(2026,6,9,12,0)):.0f}° (12:00) "
      f"(rising: a downward step during bags is against-trend = solar)")
fig.suptitle("9 Jun 2026 hood session — far-range level vs time; bin-bags added 10:34-11:39 "
             "(solar-hermeticity test). A downward step while bagged = solar leak.", fontweight="bold", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_hood_trashbag_test.png"
fig.savefig(OUT, dpi=150); print("saved", OUT)
