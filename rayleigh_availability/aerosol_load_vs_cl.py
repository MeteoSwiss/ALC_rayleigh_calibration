# -*- coding: utf-8 -*-
"""Does the aerosol load explain the recovered nights' constant? — per-night test.

The ratio-profile diagnostic shows the nights v2.2 recovers are not clean nights that noise got
rejected: on all sites tested, signal/molecular falls ~19 %/km right through the fit band, while on
the nights v2 already had it is flat. A decline of that size is what two-way aerosol extinction
looks like (19 %/km => alpha ~ 0.1 /km), so those nights carry a real aerosol column BELOW the
window the fit ends up using.

If the retrieved constant is biased by that column, the bias must scale with it. This measures, per
night: the slope of R(z) over the 2-6 km band (the aerosol proxy, independent of the window search)
against the night's constant relative to the stream's own level. A clear negative correlation means
the transmission below the window is under-corrected and the bias is PREDICTABLE -- which makes it
correctable, or at least reportable per night. No correlation means the aerosol column is not what
sets the constant and the offset has another cause.

Run:  python rayleigh_availability/aerosol_load_vs_cl.py [stream-substring ...]
Out:  rayleigh_availability/figs/aerosol_load_vs_cl.png
"""
from __future__ import annotations
import json
import logging
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import indicators as IND                                                        # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
BASE, CAND = DATA / "baselines", DATA / "candidates"
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
FIG = REPO / "rayleigh_availability" / "figs"
CFG = "N2.5"
BAND = (2000.0, 6000.0)
MAX_NIGHTS = 200


def _slope(args):
    """(date, R-slope in %/km over BAND) for one night, from the fit inputs alone."""
    label, wmo, ident, typ, lat, lon, alt, ds = args
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
    from calibration.config import InstrumentType
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = Path("D:/E-PROFILE_L1_2026")
    o.data_level = DataLevel.L1
    o.molecular_method = "eprof_v2.2"
    o.molecular_params = dict(max_chi2red=2.5)
    o.plot_main = o.plot_all = False
    o.folder_output = DATA / "tmp"
    info = InstrumentInfo(site_name=label, wmo_id=wmo, identifier=ident,
                          instrument_type=InstrumentType(typ),
                          latitude=lat, longitude=lon, altitude=alt)
    fit = {}
    try:
        calibrate_rayleigh(ds, info, o, fit_inputs_out=fit)
    except Exception:
        return ds, np.nan
    sig, pmol, rng = fit.get("signal"), fit.get("p_mol"), fit.get("range_alc")
    if sig is None or pmol is None or rng is None:
        return ds, np.nan
    r = np.asarray(rng, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.asarray(sig, float) / np.asarray(pmol, float)
    ok = np.isfinite(ratio) & (r >= BAND[0]) & (r <= BAND[1]) & (ratio > 0)
    if ok.sum() < 20:
        return ds, np.nan
    # slope in %/km of the level itself, so it is comparable across instruments
    a = np.polyfit(r[ok] / 1000.0, ratio[ok], 1)
    return ds, float(a[0] / np.median(ratio[ok]) * 100.0)


def main():
    pats = [a.upper() for a in sys.argv[1:]] or ["PAYERNE", "LINDENBERG", "GOTTFRIEDING"]
    insts = [i for i in MANIFEST
             if i["group"] == "CHM15k" and any(p in i["label"].upper() for p in pats)]
    workers = int(os.environ.get("RA_WORKERS", str(max(1, (os.cpu_count() or 8) - 2))))
    fig, axes = plt.subplots(1, len(insts), figsize=(4.3 * len(insts), 4.4), sharey=True)
    axes = np.atleast_1d(axes)
    summary = {}
    for ax, inst in zip(axes, insts):
        ref = json.loads((BASE / f"base_eprof_v2_{inst['label']}.json").read_text())
        new = json.loads((CAND / f"cand_{CFG}_{inst['label']}.json").read_text())
        val = {d: v for d, v in new.items() if IND.is_valid(v[0]) and v[1]}
        dates = sorted(val)
        step = max(1, len(dates) // MAX_NIGHTS)
        dates = dates[::step]
        args = [(inst["label"], inst["wmo"], inst["ident"], inst["type"],
                 inst["lat"], inst["lon"], inst["alt"], d) for d in dates]
        sl = {}
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for fut in as_completed([ex.submit(_slope, a) for a in args]):
                d, s = fut.result()
                if np.isfinite(s):
                    sl[d] = s
        if not sl:
            continue
        lvl = np.median([val[d][1] for d in sl])
        x = np.array([sl[d] for d in sorted(sl)])
        y = np.array([val[d][1] / lvl for d in sorted(sl)])
        isrec = np.array([not IND.is_valid(ref.get(d, [0])[0]) for d in sorted(sl)])
        ax.plot(x[~isrec], y[~isrec], "o", ms=4, color="#777", alpha=.75,
                label=f"kept by v2 (n={int((~isrec).sum())})")
        ax.plot(x[isrec], y[isrec], "o", ms=4, color="#1f77b4", alpha=.75,
                label=f"recovered (n={int(isrec.sum())})")
        ok = np.isfinite(x) & np.isfinite(y)
        r_all = np.corrcoef(x[ok], y[ok])[0, 1] if ok.sum() > 5 else np.nan
        if ok.sum() > 5:
            p = np.polyfit(x[ok], y[ok], 1)
            xs = np.linspace(np.min(x[ok]), np.max(x[ok]), 20)
            ax.plot(xs, np.polyval(p, xs), "-", color="#d62728", lw=1.6)
        summary[inst["label"]] = dict(r=float(r_all), n=int(ok.sum()),
                                      med_slope_kept=float(np.median(x[~isrec])) if (~isrec).any() else np.nan,
                                      med_slope_rec=float(np.median(x[isrec])) if isrec.any() else np.nan)
        ax.axhline(1.0, color="k", lw=.6, ls=":")
        ax.set_title(f"{inst['site'][:18]}   r = {r_all:+.2f}", fontsize=10)
        ax.set_xlabel("slope of signal/molecular, 2-6 km (%/km)\n(more negative = more aerosol)")
        ax.grid(alpha=.3)
        ax.legend(fontsize=8)
    axes[0].set_ylabel(r"$C_L$ / stream median")
    fig.suptitle("Aerosol load below the fit window vs the retrieved constant",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "aerosol_load_vs_cl.png", dpi=140)
    plt.close(fig)
    (DATA / "aerosol_load_vs_cl.json").write_text(json.dumps(summary, indent=1))
    for k, v in summary.items():
        print(f"{k:30s} r={v['r']:+.2f} n={v['n']:4d} | median slope kept "
              f"{v['med_slope_kept']:+6.1f} %/km, recovered {v['med_slope_rec']:+6.1f} %/km")
    print(f"-> {FIG / 'aerosol_load_vs_cl.png'}")


if __name__ == "__main__":
    main()
