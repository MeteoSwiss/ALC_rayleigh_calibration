"""All-CL61 coherence + the Aosta seasonal-cycle question.

For every CL61 with both calibration methods (Payerne, Lindenberg, Aosta, Camborne):
  - monthly medians of C_L(rayleigh) and C_L(cloud) from the calout kalman series;
  - the monthly two-way WV transmission T2_wv at a 3-5 km window from the 1 deg CAMS;
  - correlation of the monthly method ratio with T2_wv -> is the Rayleigh seasonality the
    seasonally-varying WV (over-)correction? (C_L_ray ~ 1/T2_wv; the cloud method integrates
    0.1-2.4 km where the WV lever is ~4x smaller.)
Figure: per station, monthly C_L (blue=Rayleigh, dark grey=cloud) + T2_wv; plus ratio-vs-T2 scatter.
"""
import sys, csv, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from calibration.water_vapor_correction.water_vapor import (
    cams_water_vapor_profile, two_way_wv_transmission, DEFAULT_ABS_CROSS_SECTION)
from calibration.io.cams import ensure_cams_file

CALOUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E_PROFILE_calout_2025_2026")
STATIONS = {
    "Payerne":    dict(key="0-20000-0-06610_C", lat=46.81, lon=6.94, alt=491.0),
    "Lindenberg": dict(key="0-20000-0-10393_C", lat=52.21, lon=14.12, alt=98.0),
    "Aosta":      dict(key="0-380-5-1_B",       lat=45.74, lon=7.36, alt=570.0),
    "Camborne":   dict(key="0-20000-0-03808_C", lat=50.22, lon=-5.33, alt=88.0),
}
RAY_COL, CLD_COL = "#1f77b4", "#404040"


def kalman_series(key, method):
    f = CALOUT / key / f"{key}_kalman.csv"
    dd, vv = [], []
    for r in csv.DictReader(open(f, encoding="utf-8")):
        if r.get("method") != method:
            continue
        try:
            v = float(r["kalman"]); s = str(r["date"])
        except (TypeError, ValueError):
            continue
        if v > 0 and len(s) >= 8:
            dd.append(np.datetime64(f"{s[:4]}-{s[4:6]}-{s[6:8]}")); vv.append(v)
    return np.array(dd), np.array(vv)


def monthly_median(dd, vv):
    out = {}
    mm = dd.astype("datetime64[M]")
    for m in np.unique(mm):
        out[str(m)] = float(np.median(vv[mm == m]))
    return out


_t2cache = {}
def t2_month(lat, lon, alt, month):
    key = (round(lat, 2), month)
    if key in _t2cache:
        return _t2cache[key]
    ds = month.replace("-", "") + "01"
    cams = ensure_cams_file(Path("D:/CAMS"), ds, auto_download=False)
    val = np.nan
    if cams is not None:
        m0 = np.datetime64(month + "-01")
        prof = cams_water_vapor_profile(cams, lat, lon, m0, m0 + np.timedelta64(27, "D"))
        if prof is not None:
            h, n = prof
            zg = np.arange(0, 6001, 50.0) + alt
            t2 = two_way_wv_transmission(zg, alt, h, n, DEFAULT_ABS_CROSS_SECTION, 910.74, 1.0)
            zm = (zg - alt >= 3000) & (zg - alt <= 5000)
            val = float(np.nanmean(np.asarray(t2, "f8")[zm]))
    _t2cache[key] = val
    return val


fig, axes = plt.subplots(2, len(STATIONS), figsize=(5.6 * len(STATIONS), 9))
print("station     n_mon  ratio(med)  corr(ratio,1/T2)  corr(CLray,1/T2)  corr(CLcld,1/T2)")
summary = {}
for i, (name, st) in enumerate(STATIONS.items()):
    dr, vr = kalman_series(st["key"], "rayleigh")
    dc, vc = kalman_series(st["key"], "cloud")
    if dr.size < 10 or dc.size < 10:
        print(f"{name:10s}  insufficient series (ray {dr.size}, cloud {dc.size})")
        axes[0][i].set_title(f"{name}: no dual series"); continue
    mr, mc = monthly_median(dr, vr), monthly_median(dc, vc)
    months = sorted(set(mr) & set(mc))
    t2 = np.array([t2_month(st["lat"], st["lon"], st["alt"], m) for m in months])
    r = np.array([mr[m] for m in months]); c = np.array([mc[m] for m in months])
    ratio = r / c
    ok = np.isfinite(t2)
    inv = 1.0 / t2[ok]
    def cor(a, b):
        return float(np.corrcoef(a, b)[0, 1]) if a.size > 3 else np.nan
    cr = cor(ratio[ok], inv); cra = cor(r[ok], inv); ccl = cor(c[ok], inv)
    print(f"{name:10s}  {len(months):4d}   {np.median(ratio):.3f}      {cr:+.2f}            {cra:+.2f}            {ccl:+.2f}")
    summary[name] = dict(months=len(months), ratio=float(np.median(ratio)), corr_ratio=cr,
                         corr_ray=cra, corr_cloud=ccl)
    x = np.array([np.datetime64(m) for m in months])
    ax = axes[0][i]
    ax.plot(x, r, "o-", color=RAY_COL, lw=1.6, label="C$_L$ Rayleigh")
    ax.plot(x, c, "s-", color=CLD_COL, lw=1.6, label="C$_L$ cloud")
    ax2 = ax.twinx()
    ax2.plot(x, t2, "^--", color="#2ca02c", lw=1.2, alpha=0.8, label="T$^2_{wv}$ (3-5 km)")
    ax2.set_ylabel("T$^2_{wv}$", color="#2ca02c"); ax2.tick_params(axis="y", labelcolor="#2ca02c")
    ax.set_title(name); ax.grid(alpha=0.3); ax.set_ylabel("monthly median C$_L$")
    if i == 0:
        ax.legend(fontsize=8, loc="upper left")
    for lab in ax.get_xticklabels():
        lab.set_rotation(30); lab.set_fontsize(7)
    axb = axes[1][i]
    axb.plot(1.0 / t2[ok], ratio[ok], "o", color="#d62728")
    axb.set_xlabel("1 / T$^2_{wv}$ (3-5 km)"); axb.set_ylabel("C$_L$ ray / cloud")
    axb.grid(alpha=0.3)
    axb.set_title(f"ratio vs WV lever: corr={cr:+.2f}")
fig.suptitle("CL61 network coherence: monthly lidar constants (blue=Rayleigh, dark grey=cloud) vs the "
             "seasonal water-vapour lever", fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_network_coherence.png")
fig.savefig(out, dpi=150)
print("saved", out)
import json
Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_coherence.json").write_text(json.dumps(summary, indent=1))
print("COHERENCE_DONE")
