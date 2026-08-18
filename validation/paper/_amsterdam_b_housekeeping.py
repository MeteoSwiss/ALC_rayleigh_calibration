"""
_amsterdam_b_housekeeping.py — is anything in unit B's HOUSEKEEPING different from A/C/D, and
does it explain the daytime near-range excess (§8.2)?

The CHM15k L1 carries a full housekeeping suite. This script compares it four-ways over the
validation window and correlates unit B's hourly bias against its own housekeeping:

  read (per unit A-D): bckgrd_rcs_0, calibration_pulse, window_transmission, status_laser,
  status_detector, temp_int / temperature_detector / temperature_optical_module / temp_ext,
  laser_pulses, average_time, laser_life_time, stddev, sci, error_ext  -> hourly medians;

  (a-d) diurnal composites of the solar-sensitive channels (background, calibration pulse,
        detector temperature, raw noise std) — a unit whose DETECTION CHAIN reacts differently
        to solar load shows up here;
  (e,f) daily medians of window transmission and laser/detector quality — contamination/aging;
  (g)   service-code statistics (fraction of records with error bits set, per unit);
  (h)   |Spearman| ranking of unit-B hourly in-band bias (vs A, exact paper matrices) against
        unit-B housekeeping, day hours only (night hours carry no excess to explain — and using
        day-only also breaks the trivial "everything is diurnal" confounding).

Medians per unit + the correlation table go to discrepancy_analysis.json ("amsterdam_b_hk").
Usage: python -m validation.paper._amsterdam_b_housekeeping
"""
from __future__ import annotations
import json
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from netCDF4 import Dataset

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
RESULTS = OUT / "discrepancy_analysis.json"
WMO = "0-20000-0-06240"
BAND = (500.0, 3000.0)

# (variable, label, K->degC?)  — everything 1-D (time,) in the CHM15k L1
HKVARS = [
    ("bckgrd_rcs_0", "background [photons/shot]", False),
    ("calibration_pulse", "calibration pulse [photons/shot]", False),
    ("window_transmission", "window transmission [%]", False),
    ("status_laser", "laser quality [%]", False),
    ("status_detector", "detector quality [%]", False),
    ("temp_int", "internal T [degC]", True),
    ("temperature_detector", "detector T [degC]", True),
    ("temperature_optical_module", "optical module T [degC]", True),
    ("temp_ext", "external T [degC]", True),
    ("laser_pulses", "laser pulses / record", False),
    ("average_time", "average time [ms]", False),
    ("laser_life_time", "laser lifetime [h]", False),
    ("stddev", "raw noise std [photons/shot]", False),
    ("sci", "sky condition index", False),
    ("error_ext", "service code", False),
]


def read_hk(ident, start, end):
    """Hourly-median housekeeping for one unit. Returns DataFrame indexed by hour."""
    from validation.paper import intercompare as IC
    d0 = datetime.strptime(start, "%Y%m%d"); d1 = datetime.strptime(end, "%Y%m%d")
    rows = {v: [] for v, _, _ in HKVARS}
    times = []
    d = d0
    while d <= d1:
        f = IC.L1_ROOT / WMO / f"{d.year}" / f"{d.month:02d}" / f"L1_{WMO}_{ident}{d:%Y%m%d}.nc"
        d += timedelta(days=1)
        if not f.exists():
            continue
        try:
            with Dataset(f) as nc:
                tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
                t = IC._decode_time(np.asarray(nc.variables["time"][:], "f8"), tu)
                times.append(t)
                for v, _, conv in HKVARS:
                    if v in nc.variables:
                        x = IC._clean(nc.variables[v][:]).ravel()
                        x = x if x.size == t.size else np.full(t.size, np.nan)
                    else:
                        x = np.full(t.size, np.nan)
                    rows[v].append(x - 273.15 if conv else x)
        except Exception:
            for v, _, _ in HKVARS:
                if rows[v] and len(rows[v]) > len(times) - 1:
                    rows[v].pop()
            continue
    t = pd.DatetimeIndex(np.concatenate(times))
    df = pd.DataFrame({v: np.concatenate(rows[v]) for v, _, _ in HKVARS}, index=t).sort_index()
    # error codes: keep raw for bit statistics, hourly max otherwise loses rare bits
    err = df.pop("error_ext")
    hourly = df.resample("1h").median()
    hourly["error_any"] = (err.fillna(0) != 0).resample("1h").mean()  # fraction of records w/ error
    return hourly, err


def main():
    warnings.filterwarnings("ignore")
    from validation.paper import run_paper_validation as RPV
    from validation.paper import intercompare as IC
    from validation.paper import figures as FIG
    from validation.paper.calib_benchmark import BENCHMARK
    from calibration.sensitivity.noise import solar_elevation
    from scipy.stats import spearmanr

    st = BENCHMARK["amsterdam"]
    idents = ["A", "B", "C", "D"]

    # ---------------- housekeeping, all four units --------------------------------------------
    HK, ERR = {}, {}
    for ident in idents:
        print(f"housekeeping {ident} ...", flush=True)
        HK[ident], ERR[ident] = read_hk(ident, st["start"], st["end"])

    # ---------------- hourly bias B vs A on the exact paper matrices --------------------------
    print("paper pipeline A + B (bias series) ...", flush=True)
    items = []
    for ch in st["channels"]:
        if ch["ident"] not in ("A", "B"):
            continue
        l1 = IC.read_l1(ch["wmo"], ch["ident"], st["start"], st["end"])
        on = RPV.channel_beta(l1, ch, 1064.0, apply_wv=True)
        items.append(dict(ch=ch, l1=l1, **on))
    R = RPV.grid_and_stats(items, iref=0)
    z = np.asarray(R["altGrid"]) - R["station"]["altitude"]
    band = (z >= BAND[0]) & (z <= BAND[1])
    A, B = R["beta"][0], R["beta"][1]
    tt = pd.DatetimeIndex(np.asarray(R["time_sync"]).astype("datetime64[s]"))
    with np.errstate(all="ignore"):
        rel = (B[:, band] - A[:, band]) / A[:, band]
        rel[~(np.isfinite(rel) & (A[:, band] > 0))] = np.nan
        bias = 100 * np.nanmedian(rel, axis=1)                    # hourly med relbias B vs A [%]
    elev = solar_elevation(tt.values.astype("datetime64[s]"), R["station"]["lat"], R["station"]["lon"])
    isday = elev > 5.0

    cols = {"A": "#1f77b4", "B": "#d62728", "C": "#2ca02c", "D": "#9467bd"}
    res = {"medians": {}, "spearman_day": {}, "error_any_pct": {}}

    fig, axes = plt.subplots(2, 4, figsize=(20, 9.5))

    # (a-d) diurnal composites of the solar-sensitive channels
    diurnal_vars = [("bckgrd_rcs_0", "(a) Solar background", "background [photons/shot]"),
                    ("calibration_pulse", "(b) Calibration pulse", "photons/shot"),
                    ("temperature_detector", "(c) Detector temperature", "degC"),
                    ("stddev", "(d) Raw noise std", "photons/shot")]
    for j, (v, ttl, ylab) in enumerate(diurnal_vars):
        ax = axes[0][j]
        for ident in idents:
            h = HK[ident]
            comp = h[v].groupby(h.index.hour).median()
            ax.plot(comp.index, comp.values, "o-", color=cols[ident], ms=3,
                    label=f"unit {ident}", lw=1.8 if ident == "B" else 1.1)
        ax.set_xlabel("hour [UTC]"); ax.set_ylabel(ylab)
        ax.grid(alpha=0.3); ax.set_title(ttl, fontsize=10)
        if j == 0:
            ax.legend(fontsize=8)

    # (e) window transmission + (f) laser/detector quality — daily medians
    ax = axes[1][0]
    for ident in idents:
        d = HK[ident]["window_transmission"].resample("1D").median()
        ax.plot(d.index, d.values, "-", color=cols[ident], lw=1.8 if ident == "B" else 1.1)
    ax.set_ylabel("window transmission [%]"); ax.grid(alpha=0.3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d%b")); ax.set_title("(e) Window transmission", fontsize=10)

    ax = axes[1][1]
    for ident in idents:
        dl = HK[ident]["status_laser"].resample("1D").median()
        dd = HK[ident]["status_detector"].resample("1D").median()
        ax.plot(dl.index, dl.values, "-", color=cols[ident], lw=1.8 if ident == "B" else 1.1)
        ax.plot(dd.index, dd.values, "--", color=cols[ident], lw=1.4 if ident == "B" else 0.9)
    ax.set_ylabel("quality [%]  (solid laser / dashed detector)"); ax.grid(alpha=0.3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d%b")); ax.set_title("(f) Laser / detector quality", fontsize=10)

    # (g) service-code statistics
    ax = axes[1][2]
    frac = []
    for ident in idents:
        e = ERR[ident].fillna(0).astype("int64")
        pct = 100.0 * float((e != 0).mean())
        frac.append(pct)
        res["error_any_pct"][ident] = pct
    ax.bar(range(4), frac, color=[cols[i] for i in idents])
    ax.set_xticks(range(4)); ax.set_xticklabels([f"unit {i}" for i in idents])
    ax.set_ylabel("records with any service bit set [%]")
    ax.grid(alpha=0.3, axis="y"); ax.set_title("(g) Service codes (error_ext != 0)", fontsize=10)

    # (h) |Spearman| of unit-B DAYTIME hourly bias vs unit-B housekeeping
    ax = axes[1][3]
    hkB = HK["B"].reindex(tt.floor("h"))
    rows = []
    for v, lab, _ in HKVARS:
        if v == "error_ext":
            continue
        x = hkB[v].values.astype(float)
        m = isday & np.isfinite(bias) & np.isfinite(x)
        if m.sum() < 50 or np.nanstd(x[m]) == 0:
            continue
        rho = float(spearmanr(x[m], bias[m]).statistic)
        rows.append((lab, rho, int(m.sum())))
        res["spearman_day"][v] = dict(rho=rho, n=int(m.sum()))
    rows.sort(key=lambda r: -abs(r[1]))
    top = rows[:8]
    ax.barh(range(len(top)), [r[1] for r in top],
            color=["#d62728" if abs(r[1]) > 0.4 else "0.6" for r in top])
    ax.set_yticks(range(len(top))); ax.set_yticklabels([r[0] for r in top], fontsize=8)
    ax.invert_yaxis(); ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Spearman rho (day hours, bias B vs A)")
    ax.grid(alpha=0.3, axis="x"); ax.set_title("(h) What tracks the daytime bias?", fontsize=10)

    # medians table (stdout + JSON)
    print("\n=== housekeeping medians (validation window) ===")
    hdr = "%-34s" + "%12s" * 4
    print(hdr % (("variable",) + tuple(f"unit {i}" for i in idents)))
    for v, lab, _ in HKVARS:
        if v == "error_ext":
            continue
        meds = [float(np.nanmedian(HK[i][v])) for i in idents]
        res["medians"][v] = dict(zip(idents, meds))
        print(("%-34s" + "%12.3f" * 4) % ((lab,) + tuple(meds)))
    print("error_ext != 0 [%]:", {i: round(res["error_any_pct"][i], 2) for i in idents})
    print("\ntop daytime-bias correlates:", [(r[0], round(r[1], 2)) for r in rows[:5]], flush=True)

    fig.suptitle("Amsterdam CHM15k A-D housekeeping: what is different about unit B?",
                 fontweight="bold", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "fig_amsterdam_b_housekeeping.png", dpi=160)
    plt.close(fig)
    print("-> fig_amsterdam_b_housekeeping.png", flush=True)

    J = {}
    if RESULTS.is_file():
        try:
            J = json.loads(RESULTS.read_text())
        except Exception:
            J = {}
    J["amsterdam_b_hk"] = res
    RESULTS.write_text(json.dumps(J, indent=1))
    print("AMSTERDAM_B_HK_DONE", flush=True)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    main()
