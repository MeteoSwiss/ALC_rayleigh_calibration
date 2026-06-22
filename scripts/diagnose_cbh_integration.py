"""Diagnose whether the (high-biased) Vaisala CBH harms the O'Connor cloud-calibration integral.

Vaisala ceilometers report a cloud base where the *integrated* signal crosses a threshold
("a pilot can no longer see the ground"), i.e. some optical depth INTO the cloud -- so the
reported CBH sits ABOVE the first cloud particles. The worry: is "the beginning of the cloud"
left out of the calibration integral?

Key code fact this script tests empirically (calibration/cloud/calibration.py):
  * calculate_lidar_ratio (:1069) integrates beta over the FIXED gate window
    [cal_minheight, cal_maxheight] -- it does NOT integrate from CBH. So the lower cloud should
    be captured regardless of the reported CBH.
  * apply_cloud_filters keys the peak/aerosol filters off the beta PEAK (max_idx, :1127), and
    uses the reported CBH only as a coarse height gate (Filter 4, :1170).

For every profile SELECTED for calibration we measure:
  reported CBH ; beta-peak height ; cumulative-B accumulation heights (10/50/90% of the integral)
  ; the FRACTION of the integral B that lies BELOW the reported CBH ; a beta-based "optical onset".
If a large fraction of B sits below the reported CBH yet the integral still captures it, the
calibration is robust to the CBH definition (the worry is only about the DIAGNOSTIC annotation).

Run:  python scripts/diagnose_cbh_integration.py [OUTDIR]
Outputs: <OUTDIR>/cbh_integration_profile.png, cbh_integration_aggregate.png, cbh_integration_stats.json
"""
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import run_all_l1_2026 as R  # noqa: E402
from calibration.cloud.calibration import (  # noqa: E402
    CloudCalConfig, liquid_cloud_calibration_from_data, read_ceilometer_data, set_defaults)

OUTDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "C:/DATA/Projects/202606_E-PROFILE_calibration/cbh_diag")
OUTDIR.mkdir(parents=True, exist_ok=True)
census = json.loads(R.CENSUS.read_text(encoding="utf-8"))

KEY_DETAIL = "0-20000-0-06479"     # KLEINE_BROGEL, Vaisala CL51 -- the worked example
DATE_DETAIL = datetime(2026, 3, 1)


def _cfg(s, fp):
    return set_defaults(CloudCalConfig(
        nc_file=str(fp), instrument=s["type"], apply_wv_correction=True,
        apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
        cams_folder=str(R.CAMS), abs_cs_lookup_table=str(R.WV_LUT),
        station_latitude=s["lat"], station_longitude=s["lon"],
        average_time_s=300.0, average_range_m=10.0))


def _run(s, d):
    fp = R._l1_file(s["wmo"], s["ident"], d)
    if not fp.exists():
        return None
    cfg = _cfg(s, fp)
    data, status = read_ceilometer_data(cfg.nc_file, cfg)
    beta = getattr(data, "beta", None) if data is not None else None
    if status != 0 or beta is None or not np.any(np.isfinite(np.asarray(beta, dtype=float))):
        return None
    res = liquid_cloud_calibration_from_data(data, cfg)
    return data, res, cfg


def _profile_metrics(col, rng, cbh, cal_lo, cal_hi):
    """Per-profile geometry of the integral vs the reported CBH."""
    gate = (rng >= cal_lo) & (rng <= cal_hi)
    r = rng[gate]
    b = np.array(col[gate], dtype=float)
    b[~np.isfinite(b)] = 0.0
    if r.size < 5 or not np.any(b > 0):
        return None
    pk = int(np.argmax(b))
    peak_h, peak_b = float(r[pk]), float(b[pk])
    cumB = np.concatenate([[0.0], np.cumsum((b[1:] + b[:-1]) / 2.0 * np.diff(r))])
    totB = float(cumB[-1])
    if totB <= 0:
        return None

    def h_at(frac):
        return float(r[min(int(np.searchsorted(cumB, frac * totB)), r.size - 1)])

    frac_below_cbh = np.nan
    if cbh is not None and np.isfinite(cbh):
        j = int(np.clip(np.searchsorted(r, cbh), 0, r.size - 1))
        frac_below_cbh = float(cumB[j] / totB)
    onset_h = np.nan  # beta-based optical onset: first gate exceeding 10% of the peak
    above = np.where(b >= 0.1 * peak_b)[0]
    if above.size:
        onset_h = float(r[above[0]])
    return dict(r=r, b=b, cumB=cumB, totB=totB, peak_h=peak_h, peak_b=peak_b,
                h10=h_at(0.1), h50=h_at(0.5), h90=h_at(0.9),
                frac_below_cbh=frac_below_cbh, onset_h=onset_h,
                cbh=(float(cbh) if cbh is not None and np.isfinite(cbh) else np.nan))


def _selected(data, res):
    sc = np.asarray(getattr(res, "S_consistent", None), dtype=float)
    n_time = data.beta.shape[1]
    if sc is None or sc.size != n_time:
        return np.zeros(n_time, dtype=bool)
    return np.isfinite(sc)


def collect(s, d, sink):
    out = _run(s, d)
    if not out:
        return None
    data, res, cfg = out
    sel = _selected(data, res)
    rng = np.asarray(data.range, dtype=float)
    cbh = np.asarray(data.cbh, dtype=float) if getattr(data, "cbh", None) is not None else np.full(sel.size, np.nan)
    recs = []
    for i in np.where(sel)[0]:
        m = _profile_metrics(data.beta[:, i], rng, cbh[i], cfg.cal_minheight, cfg.cal_maxheight)
        if m:
            recs.append(m)
            sink.append(m)
    return data, res, cfg, recs


# ---------------------------------------------------------------- aggregate scan
agg = []
s_detail = next(x for x in census if R._key(x).startswith(KEY_DETAIL))
detail = collect(s_detail, DATE_DETAIL, agg)
# extend the sample across the period (same Vaisala station) until we have enough profiles
k = 0
while len(agg) < 200 and k < 92:
    d = datetime(2026, 3, 1) + timedelta(days=k)
    k += 1
    if d == DATE_DETAIL:
        continue
    try:
        collect(s_detail, d, agg)
    except Exception as exc:  # noqa: BLE001
        print(f"skip {d:%Y%m%d}: {type(exc).__name__}", flush=True)

print(f"collected {len(agg)} selected calibration profiles", flush=True)

# ---------------------------------------------------------------- stats
def _arr(key):
    return np.array([m[key] for m in agg if np.isfinite(m.get(key, np.nan))], dtype=float)

cbh = _arr("cbh"); peak = np.array([m["peak_h"] for m in agg], dtype=float)
offs = np.array([m["cbh"] - m["peak_h"] for m in agg if np.isfinite(m["cbh"])], dtype=float)
frac = _arr("frac_below_cbh")
onset_off = np.array([m["cbh"] - m["onset_h"] for m in agg
                      if np.isfinite(m["cbh"]) and np.isfinite(m["onset_h"])], dtype=float)
stats = {
    "n_profiles": len(agg),
    "cbh_minus_peak_m": {"median": float(np.median(offs)) if offs.size else None,
                          "p10": float(np.percentile(offs, 10)) if offs.size else None,
                          "p90": float(np.percentile(offs, 90)) if offs.size else None},
    "cbh_minus_onset_m": {"median": float(np.median(onset_off)) if onset_off.size else None,
                          "p90": float(np.percentile(onset_off, 90)) if onset_off.size else None},
    "frac_of_B_below_CBH": {"median": float(np.median(frac)) if frac.size else None,
                            "p10": float(np.percentile(frac, 10)) if frac.size else None,
                            "p90": float(np.percentile(frac, 90)) if frac.size else None},
}
(OUTDIR / "cbh_integration_stats.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
print(json.dumps(stats, indent=2), flush=True)

# ---------------------------------------------------------------- figures
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# (1) representative profile on the detail day
if detail and detail[3]:
    data, res, cfg, recs = detail
    m = sorted(recs, key=lambda x: -x["totB"])[len(recs) // 2]  # a median-strength profile
    r_km = m["r"] * 1e-3
    fig, (axb, axc) = plt.subplots(1, 2, figsize=(20, 8))
    # left: beta profile + markers
    axb.plot(m["b"], r_km, color="#2ca02c", lw=1.1, label=r"$\beta$ (selected profile)")
    axb.axhspan(cfg.cal_minheight * 1e-3, cfg.cal_maxheight * 1e-3, color="gold", alpha=0.08,
                label="integration window (gate range)")
    if np.isfinite(m["cbh"]):
        axb.axhline(m["cbh"] * 1e-3, color="#d62728", ls="--", lw=1.4, label=f"reported CBH = {m['cbh']:.0f} m")
    axb.axhline(m["peak_h"] * 1e-3, color="orange", ls="-", lw=1.4, label=f"β-peak = {m['peak_h']:.0f} m")
    if np.isfinite(m["onset_h"]):
        axb.axhline(m["onset_h"] * 1e-3, color="#1f77b4", ls=":", lw=1.6, label=f"β optical onset = {m['onset_h']:.0f} m")
    axb.set_xscale("log"); axb.set_xlabel(r"$\beta$ (m$^{-1}$ sr$^{-1}$)"); axb.set_ylabel("Range (km AGL)")
    axb.set_ylim(0, min((m["peak_h"] + 1500) * 1e-3, r_km.max()))
    axb.set_title(f"{s_detail.get('site','')} — {DATE_DETAIL:%Y%m%d} — selected profile")
    axb.legend(fontsize=9, loc="upper right"); axb.grid(True, alpha=0.25)
    # right: cumulative B vs height with accumulation markers
    axc.plot(m["cumB"] / m["totB"] * 100.0, r_km, color="k", lw=1.4, label="cumulative B (% of total)")
    for frac_lvl, h, c in [(10, m["h10"], "#9ecae1"), (50, m["h50"], "#3182bd"), (90, m["h90"], "#08519c")]:
        axc.axhline(h * 1e-3, color=c, ls="-", lw=1.1, label=f"{frac_lvl}% of B by {h:.0f} m")
    if np.isfinite(m["cbh"]):
        axc.axhline(m["cbh"] * 1e-3, color="#d62728", ls="--", lw=1.4,
                    label=f"reported CBH ({m['frac_below_cbh']*100:.0f}% of B below it)")
    axc.set_xlabel("cumulative B (% of total integral)"); axc.set_ylabel("Range (km AGL)")
    axc.set_ylim(0, min((m["peak_h"] + 1500) * 1e-3, r_km.max()))
    axc.set_title("Where the calibration integral B accumulates")
    axc.legend(fontsize=9, loc="lower right"); axc.grid(True, alpha=0.25)
    fig.tight_layout(); fig.savefig(OUTDIR / "cbh_integration_profile.png", dpi=150); plt.close(fig)
    print("wrote cbh_integration_profile.png", flush=True)

# (2) aggregate distributions
fig, axes = plt.subplots(1, 3, figsize=(21, 6))
if offs.size:
    axes[0].hist(offs, bins=30, color="#4c78a8", alpha=0.8)
    axes[0].axvline(float(np.median(offs)), color="r", lw=1.4, label=f"median {np.median(offs):.0f} m")
    axes[0].set_xlabel("reported CBH − β-peak (m)"); axes[0].set_ylabel("count")
    axes[0].set_title("CBH vs β-peak offset"); axes[0].legend(); axes[0].grid(True, alpha=0.25)
if frac.size:
    axes[1].hist(frac * 100.0, bins=30, color="#54a24b", alpha=0.8)
    axes[1].axvline(float(np.median(frac)) * 100.0, color="r", lw=1.4, label=f"median {np.median(frac)*100:.0f}%")
    axes[1].set_xlabel("fraction of integral B below reported CBH (%)"); axes[1].set_ylabel("count")
    axes[1].set_title("How much of B lies below the reported CBH"); axes[1].legend(); axes[1].grid(True, alpha=0.25)
if onset_off.size:
    axes[2].hist(onset_off, bins=30, color="#e45756", alpha=0.8)
    axes[2].axvline(float(np.median(onset_off)), color="k", lw=1.4, label=f"median {np.median(onset_off):.0f} m")
    axes[2].set_xlabel("reported CBH − β optical onset (m)"); axes[2].set_ylabel("count")
    axes[2].set_title("CBH bias above the optical onset"); axes[2].legend(); axes[2].grid(True, alpha=0.25)
fig.suptitle(f"CBH-vs-integral diagnostic — {s_detail.get('site','')} (Vaisala {s_detail.get('type','')}), "
             f"{len(agg)} selected profiles", fontsize=13)
fig.tight_layout(); fig.savefig(OUTDIR / "cbh_integration_aggregate.png", dpi=150); plt.close(fig)
print("wrote cbh_integration_aggregate.png", flush=True)
print("DONE", flush=True)
