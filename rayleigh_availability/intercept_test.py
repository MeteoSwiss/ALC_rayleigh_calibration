# -*- coding: utf-8 -*-
"""Free-intercept fix test on the Payerne CHM15k recovered nights.

WHY. The window fit estimates signal = a*p_mol + b with a FREE intercept b -- an additive signal
residual (electronic baseline / imperfect dark subtraction) that the fit separates cleanly from
the molecular slope. But the recorded constant then ignores b: with subtract_background=False
(the default; options.json does not set it), calculate_lidar_constant divides the RAW rcs_mean by
beta_tot, so inside the window cl_profile(z) = C + b/(a) / p_mol(z) * C ... i.e. a relative bias
b/(a*p_mol(z)) that GROWS with window height (p_mol falls steeply). Recovered (noise-tier) nights
fit ~1.5 km higher than kept nights, which is exactly where that bias explodes -- the candidate
mechanism for the v2.2-lower-than-v2 offset and the -15 %/km window-height gradient at Payerne.

THE TEST. subtract_background=True is the existing switch that subtracts fit_result.intercept in
both the Klett preparation and the C_L profile. For every recovered and kept Payerne night, run
the calibration BOTH ways (same options otherwise -> identical window, identical fit) and compare:

  C_raw      the shipped value (must reproduce the recorded candidate constant)
  C_fix      with the intercept subtracted
  rel_pred   median over the window of b/(a*p_mol)  -- the predicted relative bias of C_raw
  rel_obs    C_raw/C_fix - 1                        -- the observed one

If rel_obs tracks rel_pred night by night, the mechanism is PROVEN.
If, in addition, C_fix on recovered nights closes on the kept-night level and the C-vs-window-height
gradient collapses, the fix is VALIDATED.
Kept nights double as the no-regression control: their windows are low (p_mol large), so
C_fix ~ C_raw there -- the fix must not move the population v2 already calibrates.

Run:  python rayleigh_availability/intercept_test.py [--sample N]
Out:  <DATA>/intercept_test_payerne.json  +  rayleigh_availability/figs/intercept_fix_payerne.png
"""
from __future__ import annotations
import json
import logging
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel  # noqa: E402
from calibration.config import InstrumentType  # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
L1_ROOT = Path("D:/E-PROFILE_L1_2026")
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
FIG = REPO / "rayleigh_availability" / "figs"
OUT_JSON = DATA / "intercept_test_payerne.json"
VALID = (1.0, 0.5)


def _inst():
    for i in MANIFEST:
        if "PAYERNE" in i["site"].upper() and i["group"] == "CHM15k":
            return i
    raise SystemExit("Payerne CHM15k not in manifest")


def _options(fix: bool) -> CalibrationOptions:
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.molecular_method = "eprof_v2.2"
    o.molecular_params = dict(max_chi2red=2.5)          # the registered N2.5 candidate
    o.plot_main = o.plot_all = False
    o.folder_output = DATA / "tmp_intercept"
    o.subtract_background = bool(fix)
    return o


def run_night(args):
    ds, pop, c_rec = args
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    inst = _inst()
    info = InstrumentInfo(site_name=inst["label"], wmo_id=inst["wmo"], identifier=inst["ident"],
                          instrument_type=InstrumentType(inst["type"]),
                          latitude=inst["lat"], longitude=inst["lon"], altitude=inst["alt"])
    row = dict(date=ds, pop=pop, c_rec=c_rec)
    fit = {}
    try:
        r_raw = calibrate_rayleigh(ds, info, _options(False), fit_inputs_out=fit)
        r_fix = calibrate_rayleigh(ds, info, _options(True))
    except Exception as exc:
        row["error"] = f"{type(exc).__name__}: {exc}"[:120]
        return row
    row["flag_raw"] = float(r_raw.flag)
    row["flag_fix"] = float(r_fix.flag)
    if r_raw.flag in VALID:
        row["c_raw"] = float(r_raw.lidar_constant)
        row["u_raw"] = float(r_raw.uncertainty)
        row["bot_asl"] = float(r_raw.calibration_bottom_height)
        row["top_asl"] = float(r_raw.calibration_top_height)
    if r_fix.flag in VALID:
        row["c_fix"] = float(r_fix.lidar_constant)
        row["u_fix"] = float(r_fix.uncertainty)
    # In-window free/forced fits on the identical inputs the pipeline used.
    if fit and r_raw.flag in VALID:
        rng = np.asarray(fit["range_alc"], float)
        sig = np.asarray(fit["signal"], float)
        pm = np.asarray(fit["p_mol"], float)
        alt = float(fit["altitude"])
        lo, hi = row["bot_asl"] - alt, row["top_asl"] - alt
        m = (rng >= lo) & (rng <= hi) & np.isfinite(sig) & np.isfinite(pm)
        if m.sum() >= 5:
            x, y = pm[m], sig[m]
            a, b = np.polyfit(x, y, 1)
            row["slope_free"] = float(a)
            row["intercept"] = float(b)
            row["slope_forced"] = float(np.sum(x * y) / np.sum(x * x))
            with np.errstate(divide="ignore", invalid="ignore"):
                row["rel_pred"] = float(np.median(b / (a * x)))
            row["win_mid_agl"] = float(0.5 * (lo + hi))
    return row


def main():
    inst = _inst()
    lab = inst["label"]
    base = json.loads((DATA / "baselines" / f"base_eprof_v2_{lab}.json").read_text())
    cand = json.loads((DATA / "candidates" / f"cand_N2.5_{lab}.json").read_text())

    def ok(rec, d):
        return d in rec and rec[d][0] in VALID and rec[d][1]

    recovered = sorted(d for d in cand if ok(cand, d) and not ok(base, d))
    kept = sorted(d for d in cand if ok(cand, d) and ok(base, d))
    ns = int(sys.argv[sys.argv.index("--sample") + 1]) if "--sample" in sys.argv else 0
    if ns:
        recovered = recovered[:: max(1, len(recovered) // ns)][:ns]
        kept = kept[:: max(1, len(kept) // ns)][:ns]
    jobs = [(d, "recovered", float(cand[d][1])) for d in recovered] + \
           [(d, "kept", float(cand[d][1])) for d in kept]
    print(f"{lab}: {len(recovered)} recovered + {len(kept)} kept nights, 2 runs each", flush=True)

    rows = []
    with ProcessPoolExecutor(max_workers=int(os.environ.get("RA_WORKERS", "10"))) as ex:
        futs = [ex.submit(run_night, j) for j in jobs]
        for k, fut in enumerate(as_completed(futs), 1):
            rows.append(fut.result())
            if k % 20 == 0:
                print(f"  {k}/{len(jobs)}", flush=True)
    rows.sort(key=lambda r: r["date"])
    OUT_JSON.write_text(json.dumps(rows, indent=1), encoding="utf-8")

    # ---- summary ---------------------------------------------------------
    def med(vals):
        v = [x for x in vals if x is not None and np.isfinite(x)]
        return float(np.median(v)) if v else np.nan

    good = [r for r in rows if "c_raw" in r]
    rep = [abs(r["c_raw"] / r["c_rec"] - 1) for r in good if r.get("c_rec")]
    print(f"\nreproduction |C_raw/C_rec-1|: median {med(rep):.2e}, max {max(rep) if rep else np.nan:.2e}")

    for pop in ("kept", "recovered"):
        g = [r for r in good if r["pop"] == pop]
        both = [r for r in g if "c_fix" in r]
        print(f"\n{pop} ({len(g)} nights, {len(both)} valid both ways)")
        print(f"  median C_raw {med([r['c_raw'] for r in g]):.3e}   "
              f"median C_fix {med([r.get('c_fix') for r in g]):.3e}")
        print(f"  median rel_obs (C_raw/C_fix-1) {med([r['c_raw']/r['c_fix']-1 for r in both])*100:+.1f}%   "
              f"median rel_pred (b/(a*p_mol))    {med([r.get('rel_pred') for r in g])*100:+.1f}%")
        fl = [r for r in g if r.get("flag_fix") not in VALID]
        if fl:
            print(f"  fix-run lost {len(fl)} nights: " +
                  ", ".join(f"{r['date']}(flag {r.get('flag_fix')})" for r in fl[:8]))

    # mechanism closure: rel_obs vs rel_pred night by night
    pairs = [(r["c_raw"] / r["c_fix"] - 1, r["rel_pred"]) for r in good
             if "c_fix" in r and "rel_pred" in r]
    if len(pairs) >= 8:
        o, p = np.array(pairs).T
        cc = np.corrcoef(o, p)[0, 1]
        print(f"\nmechanism closure over {len(pairs)} nights: corr(rel_obs, rel_pred) = {cc:.3f}, "
              f"median |obs-pred| = {np.median(np.abs(o-p))*100:.1f} pts")

    # gradient vs window height, before and after (recovered+kept pooled)
    for tag, key in (("raw", "c_raw"), ("fix", "c_fix")):
        pts = [(r["win_mid_agl"], r[key]) for r in good if key in r and "win_mid_agl" in r]
        if len(pts) >= 8:
            z, c = np.array(pts).T
            sl = np.polyfit(z / 1000.0, np.log(c), 1)[0] * 100
            print(f"dC/dz ({tag}, all nights): {sl:+.1f} %/km over {len(pts)} nights")

    # offset recovered-vs-kept, before and after
    for tag, key in (("raw", "c_raw"), ("fix", "c_fix")):
        mk = med([r.get(key) for r in good if r["pop"] == "kept"])
        mr = med([r.get(key) for r in good if r["pop"] == "recovered"])
        if np.isfinite(mk) and np.isfinite(mr):
            print(f"recovered vs kept ({tag}): {(mr/mk-1)*100:+.1f}%")

    _figure(rows)
    print(f"\nrows -> {OUT_JSON}")


def _figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    good = [r for r in rows if "c_raw" in r and "c_fix" in r]
    if not good:
        return
    col = {"kept": "#777777", "recovered": "#1f77b4"}
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15.5, 4.6))

    # (1) time series raw vs fix
    for pop in ("kept", "recovered"):
        g = sorted((r for r in good if r["pop"] == pop), key=lambda r: r["date"])
        t = [datetime.strptime(r["date"], "%Y%m%d") for r in g]
        ax1.plot(t, [r["c_raw"] for r in g], "x", ms=5, color=col[pop], alpha=0.55,
                 label=f"{pop} — raw (shipped)")
        ax1.plot(t, [r["c_fix"] for r in g], "o", ms=4, color=col[pop],
                 label=f"{pop} — intercept-subtracted")
    kmed = np.median([r["c_fix"] for r in good if r["pop"] == "kept"])
    ax1.axhline(kmed, color="#d62728", lw=1.2, ls=":", label="kept median (fix)")
    ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax1.set_ylabel(r"$C_L$")
    ax1.set_title("Payerne CHM15k — raw vs intercept-subtracted")
    ax1.legend(fontsize=7.5)
    ax1.grid(alpha=0.3)

    # (2) C vs window mid-height: the gradient must collapse
    for key, mk, lab in (("c_raw", "x", "raw"), ("c_fix", "o", "fix")):
        pts = [(r["win_mid_agl"], r[key], r["pop"]) for r in good if "win_mid_agl" in r]
        z = np.array([p[0] for p in pts]) / 1000.0
        c = np.array([p[1] for p in pts])
        cl = [col[p[2]] for p in pts]
        ax2.scatter(c, z, marker=mk, s=22, c=cl, alpha=0.7,
                    label=f"{lab}: {np.polyfit(z, np.log(c), 1)[0]*100:+.1f} %/km")
    ax2.set_xlabel(r"$C_L$")
    ax2.set_ylabel("window mid-height AGL [km]")
    ax2.set_title("Window-height dependence")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    # (3) mechanism closure: predicted vs observed relative bias
    o = [(r["c_raw"] / r["c_fix"] - 1) * 100 for r in good if "rel_pred" in r]
    p = [r["rel_pred"] * 100 for r in good if "rel_pred" in r]
    cl = [col[r["pop"]] for r in good if "rel_pred" in r]
    ax3.scatter(p, o, s=22, c=cl, alpha=0.7)
    lim = [min(p + o + [-5]), max(p + o + [5])]
    ax3.plot(lim, lim, "-", color="#d62728", lw=1.0, label="1:1")
    ax3.set_xlabel("predicted bias  median $b/(a\\,p_{mol})$  [%]")
    ax3.set_ylabel("observed bias  $C_{raw}/C_{fix}-1$  [%]")
    ax3.set_title("Mechanism closure")
    ax3.legend(fontsize=8)
    ax3.grid(alpha=0.3)

    fig.suptitle("Additive signal residual: ignoring the fit intercept biases the pointwise $C_L$",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / "intercept_fix_payerne.png"
    fig.savefig(out, dpi=140)
    print(f"figure -> {out}")


if __name__ == "__main__":
    main()
