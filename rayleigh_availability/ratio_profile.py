# -*- coding: utf-8 -*-
"""Where does dC_L/dz come from? — the shape of signal/p_mol, before any window is chosen.

A lidar constant must not depend on the altitude it is fitted at, yet across this network the
retrieved C_L does: median |dC_L/dz| ~7 %/km, and that is what makes the nights v2.2 recovers
(fitted ~1.5 km higher) sit off the ones v2 already had. Phase 3's classification mask does not
explain it -- on recovered nights the classification flags almost nothing in the band whose windows
are ineligible.

This diagnostic drops the window machinery entirely and looks at the quantity every method fits:

    R(z) = signal(z) / p_mol(z)     which is C_L wherever the atmosphere is purely molecular

R(z) is normalised per night at 3 km, then combined across nights. The shape says which explanation
survives:

  * flat                      -> no gradient; the fit height cannot matter
  * elevated below ~3-4 km,
    flat above                -> residual AEROSOL: a real layer adds backscatter low down
  * a smooth monotone trend
    across the whole 2-8 km   -> INSTRUMENTAL: a range-dependent error in the signal itself
                                 (background/afterpulse residual, overlap or r^2 normalisation),
                                 which no atmospheric screening can remove

Kept and recovered nights are drawn separately: if the two agree, the gradient is a property of the
INSTRUMENT, not of the nights v2.2 adds -- i.e. v2's own nights carry the same bias, they just never
show it because they are all fitted at the same height.

Run:  python rayleigh_availability/ratio_profile.py [stream-substring ...]
Out:  rayleigh_availability/figs/ratio_profile.png
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
GRID = np.arange(1500.0, 8001.0, 100.0)      # AGL, the band every method searches (2-6 km) + margin
Z_REF = 3000.0
MAX_NIGHTS = 60                              # per stream per group: enough for a robust median


def _night_ratio(args):
    """R(z) = signal/p_mol on GRID for one night, normalised at Z_REF. None if unusable."""
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
        return None
    sig, pmol, rng = (fit.get("signal"), fit.get("p_mol"), fit.get("range_alc"))
    if sig is None or pmol is None or rng is None:
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.asarray(sig, float) / np.asarray(pmol, float)
    ok = np.isfinite(ratio) & (np.asarray(rng, float) > 0)
    if ok.sum() < 50:
        return None
    r = np.interp(GRID, np.asarray(rng, float)[ok], ratio[ok], left=np.nan, right=np.nan)
    ref = np.interp(Z_REF, np.asarray(rng, float)[ok], ratio[ok])
    if not np.isfinite(ref) or ref <= 0:
        return None
    return r / ref


def groups_for(inst):
    ref = json.loads((BASE / f"base_eprof_v2_{inst['label']}.json").read_text()) \
        if (BASE / f"base_eprof_v2_{inst['label']}.json").exists() else {}
    new = json.loads((CAND / f"cand_{CFG}_{inst['label']}.json").read_text()) \
        if (CAND / f"cand_{CFG}_{inst['label']}.json").exists() else {}
    kept = [d for d, v in new.items() if IND.is_valid(v[0]) and IND.is_valid(ref.get(d, [0])[0])]
    rec = [d for d, v in new.items() if IND.is_valid(v[0]) and not IND.is_valid(ref.get(d, [0])[0])]
    step = max(1, len(kept) // MAX_NIGHTS)
    step_r = max(1, len(rec) // MAX_NIGHTS)
    return sorted(kept)[::step][:MAX_NIGHTS], sorted(rec)[::step_r][:MAX_NIGHTS]


def main():
    pats = [a.upper() for a in sys.argv[1:]] or ["PAYERNE", "LINDENBERG", "GOTTFRIEDING",
                                                 "AMSTERDAM_AP_SCHIP_CHM15k_A"]
    insts = [i for i in MANIFEST
             if i["group"] == "CHM15k" and any(p in i["label"].upper() for p in pats)]
    workers = int(os.environ.get("RA_WORKERS", str(max(1, (os.cpu_count() or 8) - 2))))
    fig, axes = plt.subplots(1, len(insts), figsize=(4.0 * len(insts), 5.0), sharey=True)
    axes = np.atleast_1d(axes)
    out = {}
    for ax, inst in zip(axes, insts):
        kept, rec = groups_for(inst)
        jobs = {"kept by v2": kept, "recovered by v2.2": rec}
        for (name, dates), col in zip(jobs.items(), ("#777777", "#1f77b4")):
            if not dates:
                continue
            args = [(inst["label"], inst["wmo"], inst["ident"], inst["type"],
                     inst["lat"], inst["lon"], inst["alt"], d) for d in dates]
            prof = []
            with ProcessPoolExecutor(max_workers=workers) as ex:
                for fut in as_completed([ex.submit(_night_ratio, a) for a in args]):
                    r = fut.result()
                    if r is not None:
                        prof.append(r)
            if not prof:
                continue
            p = np.array(prof)
            med = np.nanmedian(p, axis=0)
            lo = np.nanpercentile(p, 25, axis=0)
            hi = np.nanpercentile(p, 75, axis=0)
            ax.fill_betweenx(GRID, lo, hi, color=col, alpha=0.18)
            ax.plot(med, GRID, "-", color=col, lw=2.0, label=f"{name} (n={len(p)})")
            out.setdefault(inst["label"], {})[name] = med.tolist()
            # slope of R over 2-6 km, in %/km -- directly comparable with dC_L/dz
            band = (GRID >= 2000) & (GRID <= 6000) & np.isfinite(med)
            if band.sum() > 5:
                s = np.polyfit(GRID[band] / 1000.0, med[band], 1)[0] * 100.0
                out[inst["label"]][name + " slope %/km"] = float(s)
        g = IND.altitude_gradient(json.loads(
            (CAND / f"cand_{CFG}_{inst['label']}.json").read_text())
        )["slope_pct_per_km"] if (CAND / f"cand_{CFG}_{inst['label']}.json").exists() else np.nan
        ax.axvline(1.0, color="k", lw=0.8, ls=":")
        ax.axhline(Z_REF, color="k", lw=0.6, ls=":")
        ax.set_title(f"{inst['site'][:18]}\ndC$_L$/dz = {g:+.1f} %/km", fontsize=10)
        ax.set_xlabel(r"$R(z)/R(3\,\mathrm{km})$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")
        ax.set_xlim(0.6, 1.6)
    axes[0].set_ylabel("range AGL (m)")
    fig.suptitle("Signal / molecular, normalised at 3 km — a flat line means the fit height "
                 "cannot matter", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "ratio_profile.png", dpi=140)
    plt.close(fig)
    (DATA / "ratio_profile.json").write_text(json.dumps(out, indent=1))
    for lab, d in out.items():
        print(lab, {k: round(v, 2) for k, v in d.items() if k.endswith("%/km")})
    print(f"-> {FIG / 'ratio_profile.png'}")


if __name__ == "__main__":
    main()
