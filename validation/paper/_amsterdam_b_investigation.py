"""
_amsterdam_b_investigation.py — why does Amsterdam CHM15k unit B read +15 % against units A/C/D?

The station splits (discrepancy_analysis B1/B3) show the B excess is day-weighted (+23 % day vs
+9 % night) and near-range-weighted (+33 % at 0.5-1 km fading to +7 % at 2-3 km). This script digs
into the mechanism on the exact paper matrices (run_paper_validation.run_site) and produces one
landscape figure with six probes:

  (a) diurnal cycle       — median relative bias vs unit A per hour of day, B against C/D controls;
  (b) bias profile        — median (X/A - 1) per altitude gate, day vs night, B vs C;
  (c) stability           — daily in-band median relative bias time series (drift / step check);
  (d) solar elevation     — in-band bias binned by solar elevation (the "solar load" axis);
  (e) internal temperature— in-band bias binned by unit B's own internal temperature [degC];
      (d) vs (e) disentangles direct solar/stray-light effects from thermal ones: the elevation
      and temperature axes correlate, so the sharper of the two dependences points at the driver;
  (f) additive test       — median absolute difference (B - A) per gate, day vs night: a signal
      GAIN error scales with beta (relative bias flat, absolute difference decaying like beta);
      an ADDITIVE daytime offset (background/afterpulse mis-subtraction) gives an absolute
      difference that does NOT follow the signal shape.

Numbers are appended to discrepancy_analysis.json under "amsterdam_b".
Usage: python -m validation.paper._amsterdam_b_investigation
"""
from __future__ import annotations
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
RESULTS = OUT / "discrepancy_analysis.json"
BAND = (500.0, 3000.0)


def _medrel(X, A, rows=None, gates=None):
    """Median relative bias 100*med((X-A)/A) over selected rows/gates (A>0, both finite)."""
    if rows is not None:
        X, A = X[rows], A[rows]
    if gates is not None:
        X, A = X[:, gates], A[:, gates]
    m = np.isfinite(X) & np.isfinite(A) & (A > 0)
    if m.sum() < 50:
        return np.nan
    return float(100 * np.median((X[m] - A[m]) / A[m]))


def main():
    warnings.filterwarnings("ignore")
    from validation.paper import run_paper_validation as RPV
    from validation.paper import intercompare as IC
    from validation.paper import figures as FIG
    from validation.paper.calib_benchmark import BENCHMARK
    from calibration.sensitivity.noise import solar_elevation

    print("processing amsterdam through the paper pipeline ...", flush=True)
    R, cfg, _rows = RPV.run_site("amsterdam")
    z = np.asarray(R["altGrid"]) - R["station"]["altitude"]
    band = (z >= BAND[0]) & (z <= BAND[1])
    tt = np.asarray(R["time_sync"]).astype("datetime64[s]")
    hours = tt.astype("datetime64[h]").astype("int64") % 24
    elev = solar_elevation(tt, R["station"]["lat"], R["station"]["lon"])
    isday = elev > 5.0
    lab = [c["label"] for c in R["channels"]]
    iA, iB = lab.index("CHM15k A"), lab.index("CHM15k B")
    iC, iD = lab.index("CHM15k C"), lab.index("CHM15k D")
    cols = FIG.channel_colors(R["channels"])
    A = R["beta"][iA]

    # unit B internal temperature on the union hourly grid (re-read L1: run_site does not keep it)
    st = BENCHMARK["amsterdam"]
    l1b = IC.read_l1("0-20000-0-06240", "B", st["start"], st["end"])
    th = {np.datetime64(t, "h"): v for t, v in zip(np.asarray(l1b["time"]), l1b["temp_int"])}
    tempB = np.array([th.get(np.datetime64(t, "h"), np.nan) for t in tt])

    res = {}
    fig, axes = plt.subplots(2, 3, figsize=(19, 9.5))

    # (a) diurnal cycle of the in-band bias, B vs C/D controls -----------------------------------
    ax = axes[0][0]
    diurnal = {}
    for k in (iB, iC, iD):
        v = [_medrel(R["beta"][k], A, rows=(hours == h), gates=band) for h in range(24)]
        diurnal[lab[k]] = v
        ax.plot(range(24), v, "o-", color=cols[k], label=lab[k])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("hour of day [UTC]"); ax.set_ylabel("med relbias vs A [%]")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title("(a) Diurnal cycle of the bias (0.5–3 km)", fontsize=10)
    res["diurnal"] = diurnal

    # (b) bias profile day vs night, B vs C ------------------------------------------------------
    ax = axes[0][1]
    zplot = (z >= 0) & (z <= 4000)
    for k, ls_day, ls_ngt in ((iB, "-", "--"), (iC, "-", "--")):
        for sel, ls, tag in ((isday, ls_day, "day"), (~isday, ls_ngt, "night")):
            prof = np.array([_medrel(R["beta"][k][:, [i]], A[:, [i]], rows=sel)
                             for i in np.where(zplot)[0]])
            lw = 1.8 if k == iB else 1.0
            ax.plot(prof, z[zplot], ls, color=cols[k], lw=lw, label=f"{lab[k]} {tag}")
    ax.axvline(0, color="k", lw=0.8)
    for zb in BAND:
        ax.axhline(zb, ls=":", color="k", lw=0.8)
    ax.set_xlim(-25, 60); ax.set_ylim(0, 4000)
    ax.set_xlabel("med relbias vs A [%]"); ax.set_ylabel("altitude a.g.l. [m]")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title("(b) Bias profile: day (solid) vs night (dashed)", fontsize=10)

    # (c) daily stability -------------------------------------------------------------------------
    ax = axes[0][2]
    days = tt.astype("datetime64[D]")
    udays = np.unique(days)
    for k in (iB, iC, iD):
        v = [_medrel(R["beta"][k], A, rows=(days == d), gates=band) for d in udays]
        ax.plot(udays, v, ".-", color=cols[k], lw=0.8, ms=3, label=lab[k])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("daily med relbias vs A [%]"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.set_title("(c) Daily bias: stable pattern, no drift/step", fontsize=10)

    # (d) bias vs solar elevation -----------------------------------------------------------------
    ax = axes[1][0]
    edges = np.array([-90, -5, 5, 15, 25, 35, 45, 60])
    centers = 0.5 * (edges[:-1] + edges[1:])
    binned = {}
    for k in (iB, iC):
        v = [_medrel(R["beta"][k], A, rows=(elev >= e0) & (elev < e1), gates=band)
             for e0, e1 in zip(edges[:-1], edges[1:])]
        binned[lab[k]] = v
        ax.plot(centers, v, "o-", color=cols[k], label=lab[k])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("solar elevation [deg]"); ax.set_ylabel("med relbias vs A [%]")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title("(d) Bias vs solar elevation (0.5–3 km)", fontsize=10)
    res["vs_elevation"] = dict(centers=list(map(float, centers)), **binned)

    # (e) bias vs unit B internal temperature -----------------------------------------------------
    ax = axes[1][1]
    tq = np.nanpercentile(tempB, [0, 20, 40, 60, 80, 100])
    tc, tv = [], []
    for t0, t1 in zip(tq[:-1], tq[1:]):
        sel = (tempB >= t0) & (tempB < t1 + 1e-9)
        tc.append(0.5 * (t0 + t1))
        tv.append(_medrel(R["beta"][iB], A, rows=sel, gates=band))
    ax.plot(tc, tv, "o-", color=cols[iB], label="CHM15k B")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("unit B internal temperature [degC]"); ax.set_ylabel("med relbias vs A [%]")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title("(e) Bias vs internal temperature (quintiles)", fontsize=10)
    res["vs_temperature"] = dict(temp_C=list(map(float, tc)), medrel=list(map(float, tv)))

    # (f) additive-vs-gain test: median absolute difference per gate, day vs night ---------------
    ax = axes[1][2]
    for sel, ls, tag in ((isday, "-", "day"), (~isday, "--", "night")):
        dif = np.array([np.nanmedian((R["beta"][iB][:, i] - A[:, i])[sel]) for i in np.where(zplot)[0]])
        ax.plot(dif, z[zplot], ls, color=cols[iB], lw=1.8, label=f"B - A {tag}")
        medA = np.array([np.nanmedian(A[sel][:, i]) for i in np.where(zplot)[0]])
        ax.plot(medA * 0.155, z[zplot], ls, color="0.55", lw=1.0,
                label=f"0.155 x median A {tag} (pure-gain shape)")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_ylim(0, 4000)
    ax.set_xlabel(r"median $\beta$ difference [Mm$^{-1}$ sr$^{-1}$]"); ax.set_ylabel("altitude a.g.l. [m]")
    ax.grid(alpha=0.3); ax.legend(fontsize=7)
    ax.set_title("(f) Additive-vs-gain: |difference| profile against the pure-gain shape", fontsize=10)

    # headline numbers
    res["summary"] = dict(
        all=_medrel(R["beta"][iB], A, gates=band),
        day=_medrel(R["beta"][iB], A, rows=isday, gates=band),
        night=_medrel(R["beta"][iB], A, rows=~isday, gates=band),
        peak_hour=int(np.nanargmax(np.array(diurnal["CHM15k B"], dtype=float))),
        peak_hour_bias=float(np.nanmax(np.array(diurnal["CHM15k B"], dtype=float))),
    )
    print("summary:", res["summary"], flush=True)

    fig.suptitle(f"Amsterdam CHM15k unit B ({res['summary']['all']:+.0f} % vs unit A): a daytime, "
                 "near-range, unit-specific signal effect — not calibration, not overlap-model residual",
                 fontweight="bold", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "fig_amsterdam_b_investigation.png", dpi=160)
    plt.close(fig)
    print(f"-> fig_amsterdam_b_investigation.png", flush=True)

    # merge into the discrepancy JSON
    J = {}
    if RESULTS.is_file():
        try:
            J = json.loads(RESULTS.read_text())
        except Exception:
            J = {}
    J["amsterdam_b"] = res
    RESULTS.write_text(json.dumps(J, indent=1))
    print("AMSTERDAM_B_DONE", flush=True)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    main()
