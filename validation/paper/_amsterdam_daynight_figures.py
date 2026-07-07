"""
_amsterdam_daynight_figures.py — the Amsterdam four-CHM15k station figure (report Figure 3)
rendered for DAY-only and NIGHT-only hours, to make the unit-B daytime excess visible in the
station-figure format itself (profiles, scatter, histogram, curtains).

The site is processed ONCE through the paper pipeline (run_paper_validation.run_site); the day
and night variants NaN-mask the synchronized hourly matrices on solar elevation (day: elev > 5
deg, night: the complement — the SAME convention as the section-8 splits, so the numbers match
Figure 15). Masking (rather than subsetting) keeps the curtain time axis honest: excluded hours
are blank. Per-channel statistics are recomputed on each masked set and a comparison table
(all / day / night) is printed and stored in discrepancy_analysis.json ("amsterdam_daynight").

Usage: python -m validation.paper._amsterdam_daynight_figures
"""
from __future__ import annotations
import json
import sys
import warnings
from pathlib import Path

import numpy as np

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
RESULTS = OUT / "discrepancy_analysis.json"


def mask_R(R, keep, zmin, zmax, iref):
    """NaN-mask all time rows where keep is False; recompute the per-channel stats."""
    from validation.paper import intercompare as IC
    R2 = dict(R)
    kc = keep[:, None]
    R2["beta"] = [np.where(kc, b, np.nan) for b in R["beta"]]
    R2["beta_disp"] = [np.where(kc, b, np.nan) for b in R["beta_disp"]]
    R2["beta_noovl"] = [np.where(kc, b, np.nan) if b is not None else None
                        for b in R.get("beta_noovl", [None] * len(R["beta"]))]
    R2["cbh"] = [np.where(keep, c, np.nan) for c in R["cbh"]]
    z = np.asarray(R["altGrid"])
    zmask = (z >= zmin + R["station"]["altitude"]) & (z <= zmax + R["station"]["altitude"])
    ref = R2["beta"][iref]
    R2["stats"] = [IC._stats(b, ref, zmask) for b in R2["beta"]]
    return R2


def main():
    warnings.filterwarnings("ignore")
    from validation.paper import run_paper_validation as RPV
    from validation.paper import figures as FIG
    from calibration.sensitivity.noise import solar_elevation

    print("processing amsterdam through the paper pipeline ...", flush=True)
    R, cfg, _rows = RPV.run_site("amsterdam")
    tt = np.asarray(R["time_sync"]).astype("datetime64[s]")
    elev = solar_elevation(tt, R["station"]["lat"], R["station"]["lon"])
    isday = elev > 5.0
    iref = cfg["referenceChannel"]
    zmin, zmax = cfg["zMin"], cfg["zMax"]
    base = "Amsterdam (0-20000-0-06240)  —  2026-03-01 to 2026-05-31"

    res = {"convention": "day: solar elevation > 5 deg; night: complement", "channels": {}}
    variants = [("all", None, None),
                ("day", isday, OUT / "fig_amsterdam_day.png"),
                ("night", ~isday, OUT / "fig_amsterdam_night.png")]
    Rsets = {}
    for name, keep, png in variants:
        Rv = R if keep is None else mask_R(R, keep, zmin, zmax, iref)
        Rsets[name] = Rv
        if png is not None:
            tag = "DAY only (solar elev > 5°)" if name == "day" else "NIGHT only (solar elev ≤ 5°)"
            FIG.fig_multi_alc(Rv, cfg, png, f"{base}  —  {tag}")
            print(f"-> {png.name}", flush=True)

    # comparison table
    print("\n%-12s %-7s %12s %8s %12s %8s %10s" %
          ("channel", "hours", "med relbias", "log r", "relbias", "r", "N"))
    for k, ch in enumerate(R["channels"]):
        if k == iref:
            continue
        ent = {}
        for name in ("all", "day", "night"):
            s = Rsets[name]["stats"][k]
            ent[name] = dict(medrel=s["medrelbias_pct"], rlog=s["r_log"],
                             relbias=s["relbias_pct"], r=s["r"], n=s["n"])
            print("%-12s %-7s %+11.1f%% %8.2f %+11.1f%% %8.2f %10d" %
                  (ch["label"], name, s["medrelbias_pct"], s["r_log"],
                   s["relbias_pct"], s["r"], s["n"]))
        res["channels"][ch["label"]] = ent

    J = {}
    if RESULTS.is_file():
        try:
            J = json.loads(RESULTS.read_text())
        except Exception:
            J = {}
    J["amsterdam_daynight"] = res
    RESULTS.write_text(json.dumps(J, indent=1))
    print("AMSTERDAM_DAYNIGHT_DONE", flush=True)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    main()
