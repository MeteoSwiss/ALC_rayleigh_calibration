# -*- coding: utf-8 -*-
"""Detect the exact covered-telescope (dark) intervals of the Payerne campaign in the L1 data.

Why: the logbook times (Confluence "ALC dark noise measurement campaign") are approximate, and a
single atmospheric profile inside a dark window biases the dark profile b(z) by orders of magnitude
(the atmospheric near-range signal is 10^2-10^4 times the dark level). So the hood-on / hood-off
edges are taken from the data, not from the logbook.

Method, per instrument and session:
  * read the L1 rcs_0 from 90 min before the logged start to 90 min after the logged end;
  * indicator = median rcs_0 over 300-1500 m per profile (where the atmosphere is always strong,
    while under the hood only the small electronic baseline remains);
  * reference levels: S = median |indicator| outside the logged window (atmosphere), D = median
    |indicator| over the central 60 % of the logged window (dark); a profile is "covered" when the
    15-profile running median of |indicator| is below sqrt(D*S);
  * the dark window = the longest contiguous covered run that overlaps the logged window, trimmed by
    MARGIN_S on both sides (hood placement / removal, detector settling).
The result is written as JSON (windows per ident) and a review figure (one panel per session).

Run:  python rayleigh_availability/dark_windows_detect.py
"""
from __future__ import annotations
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4

REPO = Path(__file__).resolve().parents[1]
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
OUT_JSON = DATA / "dark_windows_detected.json"
OUT_FIG = REPO / "doc" / "reports" / "figs_remote_dark" / "dark_windows_detected.png"

PAD = timedelta(minutes=90)      # searched / reference margin around the logged window
# covered if the smoothed |indicator| is below the GEOMETRIC MEAN of the dark level D (core of the
# logged window) and of the atmospheric level S (outside it): robust to both the large per-profile
# noise of the CHM15k dark and the small (~10x) dark/atmosphere contrast of the CL31.
MARGIN_S = 120.0                 # trimmed inside the detected edges
SMOOTH_N = 15                    # running-median length on the indicator [profiles]
MERGE_GAP = timedelta(minutes=5) # covered runs separated by a shorter excursion are one session (CL31 noise)

# Logbook (Confluence page "ALC dark noise measurement campaign", version 31, UTC). End=None: not
# logged (removed by someone else) -> the end is searched in the data up to the end of the day.
LOGBOOK = {
    "A": [("2026-05-12 09:24", "2026-05-12 14:53"), ("2026-05-26 11:42", "2026-05-27 13:21"),
          ("2026-06-09 09:16", "2026-06-09 11:55"), ("2026-06-23 12:32", "2026-06-23 15:09"),
          ("2026-07-07 08:47", "2026-07-07 11:32"), ("2026-07-21 07:58", "2026-07-21 12:10"),
          ("2026-08-18 07:53", "2026-08-18 10:53"), ("2026-09-01 12:43", "2026-09-01 14:48"),
          ("2026-09-15 07:47", "2026-09-15 12:00"), ("2026-09-29 07:40", "2026-09-29 11:34")],
    "B": [("2026-05-12 09:33", "2026-05-12 14:55"), ("2026-05-26 11:35", "2026-05-27 13:16"),
          ("2026-06-09 09:16", "2026-06-09 11:57"), ("2026-06-23 10:12", "2026-06-23 12:13"),
          ("2026-07-07 08:46", "2026-07-07 11:19"), ("2026-07-09 09:45", "2026-07-10 06:45"),
          ("2026-07-21 08:06", "2026-07-21 12:18"), ("2026-08-18 07:50", "2026-08-18 13:25"),
          ("2026-09-01 10:25", "2026-09-01 12:31"), ("2026-09-15 07:38", "2026-09-15 10:32"),
          ("2026-09-29 07:46", "2026-09-29 09:00")],   # end not logged; the data show 09:39
    "C": [("2026-05-12 09:17", "2026-05-12 14:58"), ("2026-05-26 11:37", "2026-05-27 13:26"),
          ("2026-06-09 09:10", "2026-06-09 11:51"), ("2026-06-23 10:00", "2026-06-23 12:25"),
          ("2026-07-07 11:35", "2026-07-07 14:04"), ("2026-07-21 07:54", "2026-07-21 12:06"),
          ("2026-08-18 10:57", "2026-08-18 13:38"), ("2026-09-01 10:29", "2026-09-01 12:37"),
          ("2026-09-15 12:52", "2026-09-15 15:50"), ("2026-09-29 07:35", "2026-09-29 11:25")],
}
NAMES = {"A": "CHM15k", "B": "CL31", "C": "CL61"}
# Sessions whose logged END is only an upper bound (hood removed by someone else, time not logged):
# the dark reference level D is then taken from the first 15 min after the hood went on (the core of
# the bounded window may already be atmospheric), and the off edge comes from the data.
APPROX_END = {("B", "2026-09-29 07:46")}


def _p(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M")


def read_span(ident, t0, t1):
    """(times, indicator) for all L1 profiles in [t0, t1]."""
    T, X = [], []
    day = t0.date()
    while day <= t1.date():
        f = L1 / f"{day:%Y}" / f"{day:%m}" / f"L1_0-20000-0-06610_{ident}{day:%Y%m%d}.nc"
        if f.exists():
            with netCDF4.Dataset(f) as ds:
                t = np.ma.filled(ds["time"][:].astype("f8"), np.nan)
                rng = np.ma.filled(ds["range"][:].astype("f8"), np.nan)
                rcs = np.ma.filled(ds["rcs_0"][:].astype("f8"), np.nan)
            if rcs.shape[0] != t.size:
                rcs = rcs.T
            band = (rng >= 300) & (rng <= 1500)
            tt = np.array([datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t])
            m = (tt >= t0) & (tt <= t1)
            T.append(tt[m])
            X.append(np.nanmedian(rcs[m][:, band], axis=1))
        day += timedelta(days=1)
    if not T:
        return np.array([]), np.array([])
    return np.concatenate(T), np.concatenate(X)


def _running_median(x, n):
    h = n // 2
    pad = np.pad(x, h, mode="edge")
    return np.array([np.nanmedian(pad[i:i + n]) for i in range(x.size)])


def detect(ident, w0, w1, approx_end=False):
    lo = _p(w0)
    hi = _p(w1) if w1 else datetime.combine(lo.date(), datetime.max.time())
    t, x = read_span(ident, lo - PAD, hi + PAD)
    if t.size == 0:
        return dict(log=(w0, w1), status="no data")
    outside = (t < lo) | (t > hi)
    if approx_end:
        outside = t < lo                     # after the bound the hood is already off: atmosphere
    S = np.nanmedian(np.abs(x[outside])) if outside.any() else np.nan
    span = (hi - lo)
    core = (t >= lo + 0.2 * span) & (t <= hi - 0.2 * span)
    if approx_end:
        core = (t >= lo + timedelta(minutes=2)) & (t <= lo + timedelta(minutes=15))
    D = np.nanmedian(np.abs(x[core])) if core.any() else np.nan
    thr = np.sqrt(max(D, 1e-30) * S)
    xs = _running_median(np.abs(x), SMOOTH_N)
    covered = xs < thr
    # contiguous covered runs; keep the longest one overlapping the logged window
    runs, i = [], 0
    while i < t.size:
        if covered[i]:
            j = i
            while j + 1 < t.size and covered[j + 1]:
                j += 1
            runs.append((i, j))
            i = j + 1
        else:
            i += 1
    merged = []
    for r in runs:
        if merged and (t[r[0]] - t[merged[-1][1]]) <= MERGE_GAP:
            merged[-1] = (merged[-1][0], r[1])
        else:
            merged.append(r)
    runs = [r for r in merged if t[r[1]] >= lo and t[r[0]] <= hi]
    res = dict(log=(w0, w1), S=float(S), thr=float(thr), t=t, x=x, covered=covered)
    if not runs:
        res["status"] = "no covered run found"
        return res
    a, b = max(runs, key=lambda r: (t[r[1]] - t[r[0]]).total_seconds())
    on, off = t[a], t[b]
    use0, use1 = on + timedelta(seconds=MARGIN_S), off - timedelta(seconds=MARGIN_S)
    inside = (t >= use0) & (t <= use1)
    # n_uncovered_inside counts profiles of a merged gap (short uncovered excursions) that the window
    # still contains; only the pre-swap CL31 12 May session has any, and it is not pooled.
    res.update(status="ok", on=on, off=off, use=(use0, use1), n=int(inside.sum()),
               n_uncovered_inside=int((inside & ~covered).sum()),
               dark_level=float(np.nanmedian(x[inside])) if inside.any() else np.nan)
    return res


def main():
    out, results = {}, []
    for ident in ("A", "B", "C"):
        out[ident] = []
        print(f"== {NAMES[ident]} ({ident})")
        for w0, w1 in LOGBOOK[ident]:
            r = detect(ident, w0, w1, approx_end=(ident, w0) in APPROX_END)
            results.append((ident, r))
            if r["status"] != "ok":
                print(f"   {w0} -> {r['status']}")
                continue
            d_on = (r["on"] - _p(w0)).total_seconds() / 60
            d_off = ((r["off"] - _p(w1)).total_seconds() / 60) if w1 else float("nan")
            print(f"   log {w0} -> {w1 or 'n/a':16s} | data on {r['on']:%m-%d %H:%M:%S} off "
                  f"{r['off']:%m-%d %H:%M:%S} (dt on {d_on:+5.1f} min, off {d_off:+5.1f} min) | "
                  f"used {r['n']} prof, uncovered inside={r['n_uncovered_inside']}, "
                  f"dark/atm = {r['dark_level'] / r['S']:+.1e}")
            out[ident].append([f"{r['use'][0]:%Y-%m-%d %H:%M:%S}", f"{r['use'][1]:%Y-%m-%d %H:%M:%S}"])
    OUT_JSON.write_text(json.dumps(out, indent=1))
    print(f"-> {OUT_JSON}")
    figure(results)


def figure(results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ok = [(i, r) for i, r in results if "t" in r]
    ncol = 6
    nrow = int(np.ceil(len(ok) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(24, 3.0 * nrow), squeeze=False)
    for ax, (ident, r) in zip(axs.ravel(), ok):
        ax.semilogy(r["t"], np.abs(r["x"]), ".", ms=1.5, color="0.4")
        lo = _p(r["log"][0])
        hi = _p(r["log"][1]) if r["log"][1] else r["t"][-1]
        ax.axvspan(lo, hi, color="orange", alpha=0.18, lw=0, label="logbook")
        if r.get("status") == "ok":
            ax.axvspan(r["use"][0], r["use"][1], color="green", alpha=0.25, lw=0, label="used (data)")
        ax.axhline(r["thr"], color="r", lw=0.8, ls="--")
        ax.set_title(f"{NAMES[ident]} {r['log'][0][:10]}", fontsize=9)
        ax.tick_params(labelsize=7)
        for lab in ax.get_xticklabels():
            lab.set_rotation(30)
    for ax in axs.ravel()[len(ok):]:
        ax.axis("off")
    axs[0, 0].legend(fontsize=7)
    fig.suptitle("Payerne dark sessions: |median rcs_0| 300-1500 m per profile; orange = logbook, "
                 "green = covered interval detected in the data (trimmed), red = detection threshold",
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=110)
    print(f"-> {OUT_FIG}")


if __name__ == "__main__":
    main()
