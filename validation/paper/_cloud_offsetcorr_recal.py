# -*- coding: utf-8 -*-
"""Rerun the LIQUID-CLOUD calibration with the electronic-offset (bias) correction subtracted from
L1 rcs_0, for the Payerne CL31 (two-resonance model, cl31_b_dark.npz) and the Uccle CL51 (fixed 40 m
digitizer ripple, cl51_b_dark.npz). L1-ONLY, operational pipeline: the same liquid_cloud_calibration
(O'Connor) the dashboard/balfrin runs, reading D:/E-PROFILE_L1_2026 daily files. The offset is
injected by wrapping read_ceilometer_data: after it returns physical beta = rcs_0/C (C reported as
data.calibration_constant_applied), subtract b_phys/C (b_phys is in rcs_0 units) so the corrected
beta = (rcs_0 - b_phys)/C -- unit-exact.

Question answered: does removing the electronic offset move the cloud-derived calibration constant?
Expectation: little -- the O'Connor method integrates the strong near-range liquid-cloud return, so
it is largely offset-immune (unlike the weak-signal Rayleigh fit). This QUANTIFIES that immunity."""
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
import calibration.cloud.calibration as CC
from validation.paper.calib_benchmark import raw_cloud, kalman, OUT

D = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/"
START, END = "20260301", "20260630"      # 4 months -> a robust paired cloud-day sample
FIG = D + "fig_cloud_offsetcorr_recal.png"
FIG_REPORT = "C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_paper_report/fig_cloud_offsetcorr_recal.png"

# instrument, key, station lat/lon/alt, offset file
JOBS = [
    dict(name="Payerne CL31", itype="CL31", wmo="0-20000-0-06610", ident="B",
         lat=46.8137, lon=6.9425, alt=491.0, npz="cl31_b_dark.npz", key="0-20000-0-06610_B_cloud"),
    dict(name="Uccle CL51", itype="CL51", wmo="0-20000-0-06447", ident="A",
         lat=50.8, lon=4.35, alt=100.0, npz="cl51_b_dark.npz", key="0-20000-0-06447_A_cloud"),
]

# ---- offset injection: wrap build_cloud_input (the CeiloData builder that
# liquid_cloud_calibration now feeds; read_ceilometer_data was retired to tests/_ceilo_reader).
STATE = {"rng": None, "b": None}
_orig_adapt = CC.build_cloud_input
def _patched_adapt(cd, config, rcs_units, cal_const_applied=None):
    data = _orig_adapt(cd, config, rcs_units, cal_const_applied)
    if STATE["b"] is not None and data.beta.size and data.range.size:
        Cc = data.calibration_constant_applied
        if Cc and np.isfinite(Cc) and Cc != 0:
            off = np.interp(data.range, STATE["rng"], STATE["b"], left=0.0, right=0.0)  # rcs_0 units
            data.beta = data.beta - (off / Cc)[:, None]     # beta is (range, time)
    return data
CC.build_cloud_input = _patched_adapt


def _ch(j):
    return dict(wmo=j["wmo"], ident=j["ident"], itype=j["itype"], calib="cloud",
                label=j["name"], lat=j["lat"], lon=j["lon"], alt=j["alt"], raw=None)


results = []
for j in JOBS:
    ch = _ch(j)
    b = np.load(D + j["npz"])
    print(f"== {j['name']} ({j['key']}) ==", flush=True)
    STATE["rng"], STATE["b"] = None, None
    dn, Cn, Sn = raw_cloud(ch, START, END, level="L1")
    STATE["rng"], STATE["b"] = b["rng"].astype(float), b["b_phys"].astype(float)
    dc, Cc, Sc = raw_cloud(ch, START, END, level="L1")
    STATE["b"] = None
    # pair by date
    mn = {d: c for d, c in zip(dn, Cn)}; mc = {d: c for d, c in zip(dc, Cc)}
    days = sorted(set(mn) & set(mc))
    cn = np.array([mn[d] for d in days]); cc = np.array([mc[d] for d in days])
    dpair = np.array([(cc[i] / cn[i] - 1) * 100 for i in range(len(days))])
    med_shift = float(np.median(dpair)) if len(days) else np.nan
    p95 = float(np.nanpercentile(np.abs(dpair), 95)) if len(days) else np.nan
    print(f"  native: {len(Cn)} cloud days (median C={np.median(Cn):.3f})" if len(Cn) else "  native: 0")
    print(f"  corrected: {len(Cc)} cloud days (median C={np.median(Cc):.3f})" if len(Cc) else "  corrected: 0")
    print(f"  PAIRED {len(days)} days: median shift {med_shift:+.2f}%  |95pct| {p95:.2f}%  "
          f"max|shift| {np.nanmax(np.abs(dpair)) if len(days) else np.nan:.2f}%", flush=True)
    # Kalman both -> operational-style series
    kn = kalman(dn, Cn, Sn, normalise=False); kc = kalman(dc, Cc, Sc, normalise=False)
    for tag, res in (("native", kn), ("offsetcorr", kc)):
        if res:
            with open(OUT / f"{j['key']}_{tag}_L1.csv", "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f); w.writerow(["time", "C_daily", "C_daily_std", "C_kalman", "C_kalman_std"])
                w.writerows(res)
    results.append(dict(job=j, days=days, cn=cn, cc=cc, dpair=dpair, kn=kn, kc=kc,
                        med_shift=med_shift, p95=p95))

# ---------- figure: 2 rows (CL31, CL51), left=daily C native vs corrected, right=shift histogram
fig, ax = plt.subplots(len(results), 2, figsize=(15, 4.6 * len(results)), squeeze=False)
for i, R in enumerate(results):
    j = R["job"]; days = R["days"]
    al = ax[i][0]
    al.plot(days, R["cn"], "o", ms=4, color="0.55", label="native cloud $C$", alpha=0.8)
    al.plot(days, R["cc"], ".", ms=6, color="#1f77b4", label="offset-corrected $C$")
    if R["kn"]:
        import datetime as _dt
        tk = [_dt.datetime.strptime(r[0], "%Y-%m-%d") for r in R["kn"]]
        al.plot(tk, [r[3] for r in R["kn"]], "-", color="0.4", lw=1.4, label="Kalman (native)")
        al.plot(tk, [r[3] for r in R["kc"]], "-", color="#d62728", lw=1.4, label="Kalman (offset-corr)")
    al.set_ylabel("cloud calibration $C$")
    al.set_title(f"{j['name']} — cloud $C$: native vs offset-corrected (L1, operational O'Connor)")
    al.legend(fontsize=8, ncol=2); al.grid(alpha=0.3)
    for lb in al.get_xticklabels():
        lb.set_rotation(20); lb.set_ha("right")
    ar = ax[i][1]
    ar.hist(R["dpair"], bins=np.linspace(-3, 3, 31), color="#1f77b4", edgecolor="white")
    ar.axvline(0, color="0.5", lw=1)
    ar.axvline(R["med_shift"], color="#d62728", lw=2, ls="--",
               label=f"median {R['med_shift']:+.2f}%\n|95pct| {R['p95']:.2f}%")
    ar.set_xlabel("per-day shift  $C_{corr}/C_{native}-1$  [%]"); ar.set_ylabel("cloud days")
    ar.set_title(f"{j['name']} — offset correction barely moves cloud $C$ (offset-immune)")
    ar.legend(fontsize=9); ar.grid(alpha=0.3)
fig.suptitle("Liquid-cloud calibration with the electronic-offset (bias) correction — the strong-signal "
             "O'Connor method is offset-immune (unlike the weak-signal Rayleigh fit)",
             fontweight="bold", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(FIG, dpi=150); Path(FIG_REPORT).parent.mkdir(parents=True, exist_ok=True); fig.savefig(FIG_REPORT, dpi=150)
print("saved", FIG); print("saved", FIG_REPORT); print("CLOUD_OFFSETCORR_RECAL_DONE")
