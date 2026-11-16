# -*- coding: utf-8 -*-
"""Q2c of report 19: diurnal cycle of the CHM15k housekeeping, network-wide.

Which instrument quantities move between day and night? For every cached stream, the median
diurnal composite (local solar hour) of each housekeeping variable is computed relative to its
night-time value (22-04 solar). Of special interest:
  * calibration_pulse - the internal reference (photons/shot); if it drops by day while the signal
    is normalised by it, the apparent constant changes with time of day;
  * bckgrd_rcs_0 - solar background, the driver of photon-counting dead time;
  * temperature_optical_module / temperature_detector - regulated, but saturate on hot days
    (Rosu et al. 2026);
  * window_transmission, status_laser.

Also: per stream, the correlation of the hourly calibration-pulse anomaly with the background and
with the internal temperature, to tell a background artefact from a thermal one.

Run:  python validation/chm15k_stability/q2c_hk_diurnal.py
"""
from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FIGS, OUT, chm15k_streams, load_cache, solar_hour, stream_key  # noqa: E402

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice")

VARS = ("calibration_pulse", "bckgrd_rcs_0", "temp_int", "temperature_optical_module", "temperature_detector",
        "window_transmission", "status_laser")
RELATIVE = ("calibration_pulse", "bckgrd_rcs_0", "window_transmission", "status_laser")  # ratio to night
SEASONS = {"DJF": (12, 1, 2), "JJA": (6, 7, 8)}


def composites(cache, lon):
    names = list(cache["hk_vars"])
    nd = len(cache["day"])
    t = (cache["day"][:, None] + (np.arange(144)[None, :] + 0.5) / 144).ravel()
    sh = solar_hour(t, lon)
    hb = np.floor(sh).astype(int)
    months = ((cache["day"][:, None] * np.ones((1, 144))).ravel() / 30.4375 % 12).astype(int) + 1  # approx
    hk = cache["hk"].reshape(nd * 144, -1)
    night = (sh >= 22) | (sh < 4)
    out = {}
    for v in VARS:
        x = hk[:, names.index(v)].astype(float)
        # Reference = the same day's night median (keeps slow drifts and swaps out of the composite).
        day_idx = np.repeat(np.arange(nd), 144)
        ref = np.full(nd, np.nan)
        for d in range(nd):
            seg = x[d * 144:(d + 1) * 144][night[d * 144:(d + 1) * 144]]
            seg = seg[np.isfinite(seg)]
            if seg.size > 6:
                ref[d] = np.median(seg)
        r = ref[day_idx]
        a = 100 * (x / r - 1) if v in RELATIVE else x - r
        out[v] = dict(all=np.array([np.nanmedian(a[hb == h]) for h in range(24)]))
        for s, mm in SEASONS.items():
            sel = np.isin(months, mm)
            out[v][s] = np.array([np.nanmedian(a[(hb == h) & sel]) for h in range(24)])
        out[v]["_anom"] = a
    # Correlation of the calibration-pulse anomaly with background and T_int anomaly (daytime hours).
    day = (sh >= 9) & (sh < 17)
    cp = out["calibration_pulse"]["_anom"]
    stats = {}
    for v in ("bckgrd_rcs_0", "temp_int", "temperature_optical_module", "temperature_detector"):
        z = out[v]["_anom"]
        ok = day & np.isfinite(cp) & np.isfinite(z)
        stats["r_cp_" + v] = float(np.corrcoef(cp[ok], z[ok])[0, 1]) if ok.sum() > 100 and np.std(z[ok]) > 0 else np.nan
    for v in VARS:
        out[v].pop("_anom")
    return out, stats


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows, comps = [], {}
    for wmo, ident, site, lat, lon in chm15k_streams():
        c = load_cache(wmo, ident)
        if c is None or len(c["day"]) < 60:
            continue
        comp, st = composites(c, lon)
        comps[stream_key(wmo, ident)] = (site, lat, comp)
        row = dict(stream=stream_key(wmo, ident), site=site, lat=lat, lon=lon, **st)
        for v in VARS:
            row[f"{v}_noon"] = float(np.nanmean(comp[v]["all"][11:14]))
            row[f"{v}_noon_JJA"] = float(np.nanmean(comp[v]["JJA"][11:14]))
        rows.append(row)
        print(f"  {row['stream']:28s} {site[:18]:18s} cal_pulse noon {row['calibration_pulse_noon']:+6.1f} % "
              f"(JJA {row['calibration_pulse_noon_JJA']:+6.1f})  bg x{1 + row['bckgrd_rcs_0_noon'] / 100:6.1f}  "
              f"dT_int {row['temp_int_noon']:+5.1f} K  dT_lom {row['temperature_optical_module_noon']:+5.2f} K  "
              f"dT_det {row['temperature_detector_noon']:+5.2f} K  r(cp,bg) {st['r_cp_bckgrd_rcs_0']:+.2f} "
              f"r(cp,Tint) {st['r_cp_temp_int']:+.2f}", flush=True)
    with open(OUT / "q2c_hk_diurnal.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    FIGS.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(2, 4, figsize=(21, 8.5))
    units = {"calibration_pulse": "% vs night", "bckgrd_rcs_0": "% vs night", "temp_int": "K vs night",
             "temperature_optical_module": "K vs night", "temperature_detector": "K vs night",
             "window_transmission": "% vs night", "status_laser": "% vs night"}
    for ax, v in zip(axs.ravel(), VARS):
        allc = []
        for key, (site, lat, comp) in comps.items():
            y = comp[v]["JJA"]
            ax.plot(np.arange(24) + 0.5, y, "-", lw=0.5, alpha=0.35, color="#ff7f0e")
            allc.append(comp[v]["all"])
        allc = np.array(allc)
        ax.plot(np.arange(24) + 0.5, np.nanmedian(allc, axis=0), "k-", lw=2.2, label="network median (all seasons)")
        ax.plot([], [], "-", color="#ff7f0e", lw=1, label="individual units, JJA")
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_xlim(0, 24)
        ax.set_xticks(range(0, 25, 3))
        ax.set_xlabel("local solar hour")
        ax.set_ylabel(units[v])
        ax.set_title(v, fontsize=10)
        lo, hi = np.nanpercentile(allc, [1, 99])
        if v == "bckgrd_rcs_0":
            ax.set_yscale("symlog", linthresh=100)
        else:
            pad = 0.5 * (hi - lo) + 0.5
            ax.set_ylim(lo - pad, hi + pad)
        ax.grid(alpha=0.3)
    axs.ravel()[0].legend(fontsize=8)
    ax = axs.ravel()[-1]
    cp = np.array([r["calibration_pulse_noon_JJA"] for r in rows], float)
    bg = np.array([r["bckgrd_rcs_0_noon_JJA"] for r in rows], float)
    ax.scatter(np.log10(1 + bg / 100), cp, s=14)
    ax.set_xlabel("log10(noon / night background), JJA")
    ax.set_ylabel("calibration pulse noon vs night [%], JJA")
    ax.set_title("internal reference vs solar background", fontsize=10)
    ax.grid(alpha=0.3)
    fig.suptitle("CHM15k housekeeping: diurnal cycle relative to the same day's night (network, 2025-2026)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "q2c_hk_diurnal.png", dpi=130)
    print(f"-> {FIGS / 'q2c_hk_diurnal.png'}")


if __name__ == "__main__":
    main()
