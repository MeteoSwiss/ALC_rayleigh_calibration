# -*- coding: utf-8 -*-
"""What differs between CHM15k optical-module generations? (report 19, part 14)

The optical-module serial TUBYYxxxx carries the production year 20YY. For every stream (main module of
2025-2026) this collects, from the 10-min HK cache (night medians), the calibration NetCDF and one L1
file: the temperature coefficient kappa (intra-night and multi-day), the Rayleigh constant level, the
internal reference pulse p_calc, the pulse count per record, the night background and noise, the
regulated detector and optical-module temperatures, the internal temperature, the laser age, the
daytime drop of p_calc, the curvature, and the firmware version. Each is related to the production
year (Spearman) and summarised for <= 2015, 2016-2019 and >= 2020 modules.

Run:  python validation/chm15k_stability/q14_generation.py   (after q1, q2a, q2c, q10)
"""
from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FIGS, L1_ROOT, OUT, chm15k_streams, load_cache, load_constants, stream_key  # noqa: E402

warnings.filterwarnings("ignore")
NIGHT_BLOCKS = np.r_[0:24, 132:144]          # 00-04 and 22-24 UTC


def first_l1_attr(wmo, ident, attr):
    import netCDF4
    for f in sorted((L1_ROOT / wmo).glob(f"*/*/L1_{wmo}_{ident}*.nc"))[-1:]:
        with netCDF4.Dataset(f) as ds:
            return str(ds.getncattr(attr)) if attr in ds.ncattrs() else ""
    return ""


def read(name, key="stream"):
    return {r[key]: r for r in csv.DictReader(open(OUT / name, encoding="utf-8"))}


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    qa, qc = read("q2a_intranight.csv"), read("q2c_hk_diurnal.csv")
    m2, curv = {}, {}
    for r in csv.DictReader(open(OUT / "q1_slopes.csv", encoding="utf-8")):
        if r["arm"] == "A" and r["tvar"] == "temp_int":
            m2[(r["stream"], r["era"])] = float(r["b2"]) if r["b2"] not in ("", "nan") else np.nan
    for r in csv.DictReader(open(OUT / "q10_curvature.csv", encoding="utf-8")):
        if r["arm"] == "A":
            curv[(r["stream"], r["era"])] = float(r["c_quad"])
    rows = []
    for wmo, ident, site, lat, lon in chm15k_streams():
        key = stream_key(wmo, ident)
        cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
        if cons is None or cache is None:
            continue
        mod = max(set(cons["module"]), key=list(cons["module"]).count)
        # serials appear as "TUB140016" (L1) or "140016" (calibration NetCDF): YY = first two digits
        digits = mod[3:] if mod.upper().startswith("TUB") else mod
        if len(digits) < 6 or not digits[:2].isdigit():
            continue
        year = 2000 + int(digits[:2])
        names = list(cache["hk_vars"])
        day_mod = np.char.replace(cache["module"].astype(str), "TUB", "") == digits
        hk = cache["hk"][day_mod][:, NIGHT_BLOCKS, :]

        def med(v):
            return float(np.nanmedian(hk[..., names.index(v)]))
        cm = cons["module"] == mod
        row = dict(stream=key, site=site, module=mod, year=year,
                   kappa_in=float(qa[key]["b_T_ctrl"]) if key in qa else np.nan,
                   kappa_m2=m2.get((key, mod), np.nan), curvature=curv.get((key, mod), np.nan),
                   C_L=float(np.nanmedian(cons["c"][cm])) / 1e11,
                   p_calc=med("calibration_pulse"), pulses=med("laser_pulses"), background=med("bckgrd_rcs_0"),
                   noise=med("stddev"), T_det=med("temperature_detector"), T_lom=med("temperature_optical_module"),
                   T_int=med("temp_int"), laser_age_kh=med("laser_life_time") / 1000,
                   window=med("window_transmission"),
                   pcal_noon_drop=float(qc[key]["calibration_pulse_noon"]) if key in qc else np.nan,
                   firmware=first_l1_attr(wmo, ident, "instrument_firmware_version"))
        rows.append(row)
    with open(OUT / "q14_generation.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    yr = np.array([r["year"] for r in rows])
    groups = (("<=2015", yr <= 2015), ("2016-19", (yr >= 2016) & (yr <= 2019)), (">=2020", yr >= 2020))
    feats = ["kappa_in", "kappa_m2", "curvature", "C_L", "p_calc", "pulses", "background", "noise", "T_det", "T_lom",
             "T_int", "laser_age_kh", "window", "pcal_noon_drop"]
    print(f"{len(rows)} streams; modules per group: " + ", ".join(f"{g} {m.sum()}" for g, m in groups))
    print(f"{'feature':16s} " + " ".join(f"{g:>10s}" for g, _ in groups) + "   Spearman(year)   p")
    for f in feats:
        x = np.array([r[f] for r in rows], float)
        ok = np.isfinite(x)
        rho, p = spearmanr(yr[ok], x[ok])
        print(f"{f:16s} " + " ".join(f"{np.nanmedian(x[m]):10.4g}" for _, m in groups) + f"   {rho:+.2f}          {p:.1e}")
    fw = {}
    for r in rows:
        fw.setdefault(r["firmware"], []).append(r["year"])
    print("firmware: " + "; ".join(f"{k or '?'}: n={len(v)}, years {min(v)}-{max(v)}" for k, v in sorted(fw.items())))

    # figure: the features that matter, vs production year
    show = [("kappa_in", "κ intra-night [%/K]"), ("C_L", "Rayleigh constant C$_L$ [10$^{11}$]"),
            ("p_calc", "reference pulse p$_{calc}$ (night)"), ("pulses", "laser pulses per record (night)"),
            ("background", "night background"), ("pcal_noon_drop", "p$_{calc}$ drop at noon [%]")]
    fig, axs = plt.subplots(2, 3, figsize=(19, 8))
    for ax, (f, lab) in zip(axs.ravel(), show):
        x = np.array([r[f] for r in rows], float)
        ax.plot(yr + np.random.default_rng(0).uniform(-0.2, 0.2, len(yr)), x, ".", color="#1f77b4", alpha=0.6)
        ys = sorted(set(yr))
        ax.plot(ys, [np.nanmedian(x[yr == y]) if np.sum(yr == y) >= 2 else np.nan for y in ys], "o-", color="#d62728",
                lw=2, label="median per production year")
        for xv in (2015.5, 2019.5):
            ax.axvline(xv, color="grey", ls=":", lw=1)
        ok = np.isfinite(x)
        rho, p = spearmanr(yr[ok], x[ok])
        ax.set_title(f"{lab}   (Spearman ρ = {rho:+.2f})", fontsize=10)
        ax.set_xlabel("optical-module production year (TUB serial)")
        if f in ("background",):
            ax.set_yscale("log")
        lo, hi = np.nanpercentile(x, [2, 98])
        if f not in ("background",):
            ax.set_ylim(lo - 0.1 * (hi - lo), hi + 0.1 * (hi - lo))
        ax.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=8)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "q14_generation.png", dpi=130)
    print(f"-> {FIGS / 'q14_generation.png'}")


if __name__ == "__main__":
    main()
