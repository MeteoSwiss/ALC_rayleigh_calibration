"""fig_v09_kalman_kleine_scheidegg.py - paper Figure 7 (manuscript v0.9 onwards).

Kalman best estimate of the lidar constant C_L of the CHM15k at Kleine Scheidegg (WIGOS
0-20000-0-06736, instrument A), 1 January 2025 - 7 October 2026: nightly Rayleigh constants
(eprof_v2.2, as published by the operational dashboard), nights flagged as outliers by the global IQR
test, the annual-cycle prediction model and the improve_alc Kalman best estimate with its 1-sigma band.

Input : paper_v08/calib_dashboard/0-20000-0-06736_A_rayleigh.csv (dashboard nightly constants fetched
        2026-10-09 by rayleigh_availability/background_bias_network.py through
        validation/paper/fetch_dashboard_calib.write_stream; only the C_daily column is used).
Filter: validation/paper/kalman_improve_alc.py (operational improve_alc_calib Kalman, type CHM15k).
Output: paper_v09/fig_v09_kalman_kleine_scheidegg.png (or the path given as first argument). It is
        embedded unchanged in the manuscript.

Prints the statistics quoted in Section 5.1.3: number of nights and outliers, model, annual amplitude
A, offset C0, peak-to-peak 2A/C0 and the measurement scatter sigma_m relative to the median.
On 2026-10-10: 280 nights, 9 outliers, annual cycle, 2A/C0 = 36 %, sigma_m/median = 11.0 %.

Usage (repository root, project venv):
    python -m validation.paper.fig_v09_kalman_kleine_scheidegg [<output.png>]

History: the figure was first made on 2026-10-10 by an inline snippet that wrote to
paper_v08/fig_v09_kalman_kleine_scheidegg.png (then copied to paper_v09/). This script is that code,
transcribed unchanged apart from layout, the docstring and the output-path argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from validation.paper.kalman_improve_alc import fitfunc_sin, kalman_improve_alc  # noqa: E402

PAPER = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation")
CALIB = PAPER / "paper_v08" / "calib_dashboard" / "0-20000-0-06736_A_rayleigh.csv"
DEFAULT_OUT = PAPER / "paper_v09" / "fig_v09_kalman_kleine_scheidegg.png"


def main(argv):
    out_png = Path(argv[0]) if argv else DEFAULT_OUT
    d = pd.read_csv(CALIB, parse_dates=["time"])
    valid = d.C_daily.notna()
    t = np.array(list(d.time[valid]))
    v = d.C_daily[valid].values
    r = kalman_improve_alc(list(d.time[valid]), v, "CHM15k")
    outlier_days = set(r["outlier_dates"])
    out = np.array([pd.Timestamp(x).date() in outlier_days for x in t])
    g = pd.to_datetime(r["grid"])
    st = r["state"]
    sd = r["std"]

    if r["model"] == "annual cycle":
        amp, _phase, c0 = r["params"]
    else:
        amp, c0 = np.nan, r["params"][0]
    print(f"Kleine Scheidegg: {int(valid.sum())} nights, {int(out.sum())} outliers, model={r['model']}, "
          f"A={amp:.3e}, C0={c0:.3e}, 2A/C0={2 * amp / c0 * 100:.0f} %, "
          f"sigma_m/median={np.sqrt(r['var_meas']) / np.median(v) * 100:.1f} %")

    fig, ax = plt.subplots(figsize=(16, 4.8), dpi=200)
    ax.fill_between(g, (st - sd) / 1e11, (st + sd) / 1e11, color="tab:green", alpha=0.25, lw=0,
                    label="best estimate ±1σ")
    ax.plot(g, st / 1e11, color="tab:green", lw=1.6, label="Kalman best estimate")
    if r["model"] == "annual cycle":
        doy = np.array([x.timetuple().tm_yday for x in g])
        ax.plot(g, fitfunc_sin(doy, *r["params"]) / 1e11, ":", color="k", lw=1.2,
                label="annual-cycle prediction model")
    ax.plot(t[~out], v[~out] / 1e11, "o", color="tab:red", ms=3.5,
            label=f"nightly Rayleigh constant, version 2.2 ({int((~out).sum())} nights)")
    ax.plot(t[out], v[out] / 1e11, "x", color="0.45", ms=6, mew=1.2,
            label=f"flagged as outlier, IQR criterion ({int(out.sum())} nights)")
    ax.set_ylabel("Lidar constant $C_L$  (10$^{11}$ a.u.)")
    ax.set_xlim(pd.Timestamp("2024-12-01"), pd.Timestamp("2026-10-15"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", ncol=3, fontsize=9, framealpha=0.9)
    ax.set_title("CHM15k Kleine Scheidegg (WIGOS 0-20000-0-06736, instrument A), "
                 "1 January 2025 - 7 October 2026", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png)
    print("figure saved:", out_png)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
