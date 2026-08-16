#!/usr/bin/env python3
"""Per-PROFILE cloud-calibration dump + Hopkin-style heatmaps for 10 selected units.

The existing CBH-heterogeneity study (doc/reports/cbh_hopkin_study.md) works on per-DAY
calibration medians (cbh_scenes.csv). Hopkin et al. 2019 (AMT, Fig. 6) instead plots every
individual calibrated PROFILE. This driver re-runs the operational liquid-cloud calibration
(same config as scripts/run_network_calibration.py::_do_cloud, replicated verbatim -- WV
correction mandatory at 910 nm, aerosol transmission correction on, native PVC eta tables)
for 10 hand-picked units on their known cloud-calibration days only, harvests the
per-profile arrays that CloudCalResults already carries (all_coefficients / cbh / time /
S_apparent), and draws true per-profile Hopkin heatmaps comparable to the per-day gallery
(C:/DATA/Projects/202606_E-PROFILE_calibration/cbh_hopkin_gallery/).

Outputs
-------
- <DUMP_DIR>/<key>_profiles.npz     one per stream (time_ns, date, cbh_m, c_oconnor, s_apparent)
- gallery <key>_perprofile.png      one Hopkin heatmap per unit (per-day gallery conventions)
- doc/reports/figs_cbh_heterogeneity/perprofile_hopkin_10units.png   2x5 landscape composite
- <DUMP_DIR>/perprofile_stats.csv   per-unit slope table (day-clustered SE) + leverage check

Usage:  python rayleigh_availability/cloud_profile_dump.py [--workers 10] [--plots-only]
Resumable: a stream whose npz already exists is not re-computed (delete the npz to redo).
"""
from __future__ import annotations

# Single-threaded BLAS BEFORE numpy (the main process AND every spawned worker re-import
# this module on Windows, so the env is set everywhere before numpy loads).
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import csv
import json
import sys
import time as time_mod
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

# --- Paths (local Windows dev machine) ---------------------------------------
L1_ROOT = Path("D:/E-PROFILE_L1_2026")            # holds BOTH 2025 and 2026 (year from date)
# 0.4-deg monthlies first, then the daily archive (Jun-Aug 2026 only exists as dailies).
# NEVER D:/CAMS or D:/CAMS_run_v20 -- those are 1-degree (bad orography -> PWV bias).
CAMS_FOLDERS = "A:/CAMS_Monthly_04;D:/CAMS_daily"
WV_LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
CENSUS = REPO / "validation" / "scope_l1_2026_census.json"
SCENES_CSV = REPO / "rayleigh_availability" / "cbh_native" / "cbh_scenes.csv"
SLOPES_CSV = REPO / "rayleigh_availability" / "cbh_native" / "per_flux_slopes.csv"
# ALC_DUMP_TAG suffixe le dossier de sortie : un dossier par configuration nuage (WV nominale,
# sans WV, 910,55/1,0, 910,55/0,1) pour que le dashboard puisse afficher la heatmap Hopkin
# correspondant a la configuration choisie par l'operateur.
_TAG = os.environ.get("ALC_DUMP_TAG", "").strip()
DUMP_DIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/cloud_profile_dump"
                + (f"_{_TAG}" if _TAG else ""))
GALLERY_DIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/cbh_hopkin_gallery")
REPORT_FIG_DIR = REPO / "doc" / "reports" / "figs_cbh_heterogeneity"

# The 10 units (wmo, ident) -- see the study report for the selection rationale.
UNITS = [
    ("0-20000-0-06610", "C"),   # CL61  Payerne  (site du dashboard)
    ("0-20000-0-03809", "A"),   # CL31  network extreme negative slope
    ("0-20000-0-07627", "A"),   # CL31  network extreme positive slope
    ("0-20000-0-02055", "A"),   # CL31  flat control (NAIMAKKA)
    ("0-20000-0-06610", "B"),   # CL31  Payerne (hood reference site)
    ("0-20000-0-11487", "A"),   # CL51  3-flag suspect
    ("0-20000-0-11679", "A"),   # CL51  typical (largest n within +-1 %/km of type mean +4.5)
    ("0-196-0-CYPN", "A"),      # CL51  strong positive (Cyprus)
    ("0-20000-0-10393", "C"),   # CL61  Lindenberg (largest CL61 record)
    ("0-20000-0-03808", "C"),   # CL61  Camborne
    ("0-380-5-1", "B"),         # CL61  Aosta (alpine)
]

# ALC_DUMP_UNITS restreint la liste (ex. aux seules unites des sites du dashboard) :
# "<wmo>_<ident>,<wmo>_<ident>". Vide = les 10 unites de l'etude.
_ONLY = [s.strip() for s in os.environ.get("ALC_DUMP_UNITS", "").split(",") if s.strip()]
if _ONLY:
    _extra = [tuple(k.rsplit("_", 1)) for k in _ONLY
              if tuple(k.rsplit("_", 1)) not in UNITS]
    UNITS = [u for u in UNITS if f"{u[0]}_{u[1]}" in _ONLY] + _extra

# Plot domain -- SAME conventions as the per-day gallery (make_hopkin_heatmaps.py):
# x = per-profile constant (1/C_oconnor) as % of the unit median, y = CBH km (Y axis ALWAYS).
X_LO, X_HI = 60.0, 140.0
Y_LO, Y_HI = 0.25, 2.5
BAND = 0.25                       # mean+-SD annotation bands, km


def _l1_file(wmo: str, ident: str, ds: str) -> Path:
    """L1 path for date string YYYYMMDD (year taken from the date itself)."""
    return L1_ROOT / wmo / ds[:4] / ds[4:6] / f"L1_{wmo}_{ident}{ds}.nc"


# =============================================================================
#  Worker: one instrument-day -> per-profile arrays
# =============================================================================
def process_day(task):
    """Run the operational cloud calibration on one instrument-day and harvest the
    per-profile rows. Returns (key, ds, payload_dict_or_None, message)."""
    wmo, ident, itype, lat, lon, ds = task
    key = f"{wmo}_{ident}"
    fp = _l1_file(wmo, ident, ds)
    if not fp.exists():
        return key, ds, None, "no L1 file"
    try:
        # Imports inside the worker keep the (heavy) calibration import off the parent's
        # critical path; the config below mirrors run_network_calibration._do_cloud verbatim.
        from calibration.cloud import CloudCalConfig
        from calibration.cloud.calibration import (
            liquid_cloud_calibration_from_data, set_defaults, build_cloud_input_from_day)
        from calibration.io.instrument_day import load_instrument_day

        # ALC_WV_DISABLE=1 -> sans correction WV ; ALC_WV_SPECTRUM -> autre (lambda0, FWHM).
        # set_defaults lit la table partagee via laser_spectrum_for, donc la surcharge de
        # spectre s'applique ici aussi : un meme script produit toutes les variantes nuage.
        cfg = set_defaults(CloudCalConfig(
            nc_file=str(fp), instrument=itype,
            apply_wv_correction=(os.environ.get("ALC_WV_DISABLE") != "1"),
            apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
            cams_folder=CAMS_FOLDERS, abs_cs_lookup_table=str(WV_LUT),
            cams_folder_fallback="",
            wv_source="cams", era5_cache="",
            station_latitude=lat, station_longitude=lon,
            # No separate averaging: the shared read is already the coarse 30 s/10 m grid.
            average_time_s=0.0, average_range_m=0.0,
        ))
        d = datetime.strptime(ds, "%Y%m%d").date()
        idd = load_instrument_day([str(fp)], itype, CAMS_FOLDERS,
                                  cams_folder_fallback="",
                                  read_cams=False, build_working=False)
        if idd is None:
            return key, ds, None, "L1 unreadable"
        data = build_cloud_input_from_day(idd.slice_to_date(d), cfg)
        beta = getattr(data, "beta", None)
        if beta is None or not np.any(np.isfinite(np.asarray(beta, dtype=float))):
            return key, ds, None, "no usable signal"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = liquid_cloud_calibration_from_data(data, cfg)
        if int(getattr(res, "n_profiles", 0)) == 0 or np.asarray(res.time).size == 0:
            return key, ds, None, "no valid cloud profile"

        # Per-profile harvest. res.time / res.cbh are already valid-only slices; the full-
        # length arrays (all_coefficients, S_apparent, S_consistent) share the day's time
        # axis, and valid_idx = finite(S_consistent) is exactly the mask the core applied
        # (coefficients = S_consistent / S_theoretical) -- reconstruct it to align them.
        s_cons = np.asarray(res.S_consistent, dtype=float)
        vidx = np.isfinite(s_cons)
        c_full = np.asarray(res.all_coefficients, dtype=float)
        s_app_full = (np.asarray(res.S_apparent, dtype=float)
                      if res.S_apparent is not None else np.full_like(c_full, np.nan))
        if vidx.sum() != np.asarray(res.time).size:      # defensive: should never happen
            return key, ds, None, "valid-index mismatch"
        c = c_full[vidx]
        s_app = s_app_full[vidx]
        t_ns = np.asarray(res.time, dtype="datetime64[ns]").astype("int64")
        cbh = np.asarray(res.cbh, dtype=float)
        # keep only finite, physical pairs (C > 0 so ln(1/C) is defined downstream)
        keep = np.isfinite(c) & (c > 0) & np.isfinite(cbh)
        if not np.any(keep):
            return key, ds, None, "no finite pair"
        payload = dict(time_ns=t_ns[keep], cbh_m=cbh[keep],
                       c_oconnor=c[keep], s_apparent=s_app[keep],
                       date=np.full(int(keep.sum()), int(ds), dtype=np.int32))
        return key, ds, payload, f"OK ({int(keep.sum())} profiles)"
    except Exception as exc:  # noqa: BLE001 - one bad day must not kill the stream
        return key, ds, None, f"{type(exc).__name__}: {exc}"


# =============================================================================
#  Dump phase
# =============================================================================
def run_dump(workers: int):
    census = json.load(open(CENSUS, encoding="utf-8"))
    byk = {(s["wmo"], s["ident"]): s for s in census}

    # Day list per stream = its known cloud-calibration days (huge speedup vs full archive).
    scenes = {}
    with open(SCENES_CSV, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            scenes.setdefault(r["key"], set()).add(r["date"])

    DUMP_DIR.mkdir(parents=True, exist_ok=True)
    tasks, meta = [], {}
    for wmo, ident in UNITS:
        key = f"{wmo}_{ident}"
        s = byk.get((wmo, ident))
        if s is None:
            print(f"{key}: NOT IN CENSUS - skipped", flush=True)
            continue
        meta[key] = s
        if (DUMP_DIR / f"{key}_profiles.npz").exists():
            print(f"{key}: npz exists - skipped (resumable)", flush=True)
            continue
        days = sorted(scenes.get(key, ()))
        if not days:
            print(f"{key}: no scene days in {SCENES_CSV.name}", flush=True)
            continue
        tasks += [(wmo, ident, s["type"], s["lat"], s["lon"], ds) for ds in days]

    if tasks:
        print(f"dump: {len(tasks)} instrument-days on {workers} workers", flush=True)
        acc = {}         # key -> list of payloads
        fails = {}       # key -> {message: count}
        pending = {}     # key -> number of days still in flight (save as soon as it hits 0)
        for t in tasks:
            k = f"{t[0]}_{t[1]}"
            pending[k] = pending.get(k, 0) + 1

        def _save_stream(key):
            payloads = acc.get(key, [])
            if not payloads:
                print(f"{key}: ALL {sum(fails.get(key, {}).values())} days failed: "
                      f"{fails.get(key, {})}", flush=True)
                return
            arrays = {f: np.concatenate([p[f] for p in payloads])
                      for f in ("time_ns", "cbh_m", "c_oconnor", "s_apparent", "date")}
            order = np.argsort(arrays["time_ns"])
            arrays = {f: a[order] for f, a in arrays.items()}
            s = meta[key]
            np.savez_compressed(
                DUMP_DIR / f"{key}_profiles.npz", **arrays,
                meta=json.dumps(dict(key=key, type=s["type"], site=s["site"],
                                     lat=s["lat"], lon=s["lon"],
                                     n_days_ok=len(payloads),
                                     fails=fails.get(key, {}))))
            n = arrays["cbh_m"].size
            note = "  (< 5000 profiles - thin!)" if n < 5000 else ""
            print(f"{key}: {n} profiles over {len(payloads)} days -> npz{note}", flush=True)

        t0, done = time_mod.time(), 0
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = [pool.submit(process_day, t) for t in tasks]
            for fut in as_completed(futs):
                key, ds, payload, msg = fut.result()
                done += 1
                if payload is not None:
                    acc.setdefault(key, []).append(payload)
                else:
                    fails.setdefault(key, {})
                    fails[key][msg] = fails[key].get(msg, 0) + 1
                pending[key] -= 1
                if pending[key] == 0:      # stream complete -> persist NOW (resumable)
                    _save_stream(key)
                if done % 100 == 0 or done == len(tasks):
                    dt = time_mod.time() - t0
                    eta = dt / done * (len(tasks) - done)
                    print(f"  {done}/{len(tasks)}  ({dt:.0f}s, ETA {eta:.0f}s)", flush=True)
    return meta


# =============================================================================
#  Statistics: day-clustered slope, intra-day leverage
# =============================================================================
def cluster_ols(x, y, cluster):
    """OLS y = a + b x with CR1 cluster-robust SE on b (clusters = calendar days:
    profiles within a day share weather + one calibration state, so they are correlated
    and naive per-profile SEs would be absurdly small)."""
    n = x.size
    X = np.column_stack([np.ones(n), x])
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    e = y - X @ beta
    meat = np.zeros((2, 2))
    groups = np.unique(cluster)
    for g in groups:
        m = cluster == g
        u = X[m].T @ e[m]
        meat += np.outer(u, u)
    G, k = groups.size, 2
    c = (G / max(G - 1, 1)) * ((n - 1) / max(n - k, 1))   # CR1 small-sample factor
    V = c * XtX_inv @ meat @ XtX_inv
    return float(beta[1]), float(np.sqrt(V[1, 1])), G


def load_stream(key):
    z = np.load(DUMP_DIR / f"{key}_profiles.npz")
    meta = json.loads(str(z["meta"]))
    cinv = 1.0 / z["c_oconnor"]                    # per-profile constant (Wiegner sense)
    pct = 100.0 * cinv / np.median(cinv)           # % of the unit median
    return dict(meta=meta, cbh_km=z["cbh_m"] / 1000.0, pct=pct,
                lnc = 100.0 * np.log(cinv),        # 100*ln(1/C): slope in %/km
                day=z["date"], n=z["cbh_m"].size)


def compute_stats(day_slopes):
    rows = []
    for wmo, ident in UNITS:
        key = f"{wmo}_{ident}"
        fp = DUMP_DIR / f"{key}_profiles.npz"
        if not fp.exists():
            continue
        d = load_stream(key)
        slope, se, ndays = cluster_ols(d["cbh_km"], d["lnc"], d["day"])
        # Robustness: same fit restricted to the displayed [X_LO, X_HI] % window. The full-sample
        # OLS is fragile to the extreme-C tail (individual profiles reach ~30-350 % of the median
        # -- exactly the tail the per-DAY median statistic clips), so the windowed slope is the
        # one comparable to the heatmap's band means.
        inwin = (d["pct"] >= X_LO) & (d["pct"] <= X_HI)
        slope_w, se_w, _ = cluster_ols(d["cbh_km"][inwin], d["lnc"][inwin], d["day"][inwin])
        ds_row = day_slopes.get(key, {})
        dslope = float(ds_row.get("slope", np.nan))
        dse = float(ds_row.get("se", np.nan))
        comb = np.sqrt(se ** 2 + dse ** 2)
        agree = bool(abs(slope - dslope) <= 2.0 * comb) if np.isfinite(comb) else None
        # Intra-day leverage: per-day p90-p10 of profile CBH; median across days.
        spreads, day_med = [], []
        for g in np.unique(d["day"]):
            v = d["cbh_km"][d["day"] == g]
            day_med.append(np.median(v))
            if v.size >= 5:
                spreads.append(np.percentile(v, 90) - np.percentile(v, 10))
        intra = float(np.median(spreads)) if spreads else np.nan
        inter = float(np.percentile(day_med, 90) - np.percentile(day_med, 10))
        rows.append(dict(key=key, type=d["meta"]["type"], site=d["meta"]["site"],
                         n_profiles=d["n"], n_days=ndays,
                         slope_profile=slope, se_profile=se,
                         slope_profile_win=slope_w, se_profile_win=se_w,
                         frac_outside_win=float(1.0 - inwin.mean()),
                         slope_day=dslope, se_day=dse, agree_2se=agree,
                         intraday_p90p10_km=intra, interday_p90p10_km=inter))
    return rows


# =============================================================================
#  Plots (per-day gallery conventions; CBH ALWAYS on Y)
# =============================================================================
def make_plots(stats_rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    x_bins = np.arange(X_LO, X_HI + 1e-9, 2.0)     # 2 % bins
    y_bins = np.arange(Y_LO, Y_HI + 1e-9, 0.1)     # 100 m bins
    cmap = plt.get_cmap("jet").copy()
    cmap.set_bad("white")
    by_key = {r["key"]: r for r in stats_rows}
    REPORT_FIG_DIR.mkdir(parents=True, exist_ok=True)

    def band_stats(pct, cbh_km):
        sel = (pct >= X_LO) & (pct <= X_HI)
        edges = np.arange(Y_LO, Y_HI + 1e-9, BAND)
        out = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            v = pct[sel & (cbh_km >= lo) & (cbh_km < hi)]
            if v.size >= 5:
                out.append(((lo + hi) / 2, v.mean(), v.std(), v.size))
        return np.array(out) if out else np.empty((0, 4))

    def draw_panel(ax, d, r, annot_fs=7.0):
        h, xe, ye, im = ax.hist2d(d["pct"], d["cbh_km"], bins=[x_bins, y_bins],
                                  cmap=cmap, norm=LogNorm(vmin=1))
        im.set_array(np.ma.masked_less(im.get_array(), 1))
        bs = band_stats(d["pct"], d["cbh_km"])
        if len(bs):
            ax.errorbar(bs[:, 1], bs[:, 0], xerr=bs[:, 2], fmt="o", ms=3.5,
                        color="k", mfc="w", mew=0.9, elinewidth=1.0,
                        capsize=2.0, zorder=5)
            for yc, m, s, _n in bs:
                ax.text(X_HI - 1.0, yc, f"{m:.0f}$\\pm${s:.0f}",
                        ha="right", va="center", fontsize=annot_fs, zorder=6,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white",
                                  ec="none", alpha=0.75))
        ax.set_xlim(X_LO, X_HI)
        ax.set_ylim(Y_LO, Y_HI)
        ax.axvline(100.0, color="0.35", lw=0.8, ls="--", zorder=4)
        ttl = (f"{r['key']}  ({r['type']}, {r['site']})\n"
               f"per-profile {r['slope_profile']:+.1f}$\\pm${r['se_profile']:.1f} %/km "
               f"(n={r['n_profiles']})   per-day {r['slope_day']:+.1f} %/km")
        ax.set_title(ttl, fontsize=8)
        return im

    # --- one PNG per unit -> the per-day gallery, suffix _perprofile ---------
    for wmo, ident in UNITS:
        key = f"{wmo}_{ident}"
        if key not in by_key:
            continue
        d, r = load_stream(key), by_key[key]
        fig, ax = plt.subplots(figsize=(9.5, 6.0))
        im = draw_panel(ax, d, r, annot_fs=8.0)
        ax.set_xlabel("per-profile constant (1/C) / unit median  [%]", fontsize=9)
        ax.set_ylabel("Cloud-base height  [km]", fontsize=9)
        ax.tick_params(labelsize=8)
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label("Number of profiles", fontsize=9)
        cb.ax.tick_params(labelsize=8)
        fig.tight_layout()
        fig.savefig(GALLERY_DIR / f"{key}_perprofile.png", dpi=140)
        plt.close(fig)
    print(f"gallery: per-profile PNGs -> {GALLERY_DIR}", flush=True)

    # --- 2x5 landscape composite -> report figs ------------------------------
    fig, axes = plt.subplots(2, 5, figsize=(24.0, 9.6), sharex=True, sharey=True)
    im = None
    for ax, (wmo, ident) in zip(axes.ravel(), UNITS):
        key = f"{wmo}_{ident}"
        if key not in by_key:
            ax.set_visible(False)
            continue
        d, r = load_stream(key), by_key[key]
        im = draw_panel(ax, d, r, annot_fs=6.5)
    for i, ax in enumerate(axes.ravel()):
        if i // 5 == 1:
            ax.set_xlabel("per-profile constant (1/C) / unit median  [%]", fontsize=9)
        if i % 5 == 0:
            ax.set_ylabel("Cloud-base height  [km]", fontsize=9)
        ax.tick_params(labelsize=8)
    if im is not None:
        cb = fig.colorbar(im, ax=axes, pad=0.012, fraction=0.025)
        cb.set_label("Number of profiles", fontsize=9)
    fig.suptitle("Per-PROFILE Hopkin heatmaps -- 10 selected units "
                 "(every calibrated profile; slope = day-clustered OLS of 100$\\cdot$ln(1/C) vs CBH)",
                 fontsize=13, y=0.995)
    out = REPORT_FIG_DIR / "perprofile_hopkin_10units.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"composite -> {out}", flush=True)


# =============================================================================
#  Main
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--plots-only", action="store_true",
                    help="skip the dump (use existing npz), only stats + figures")
    args = ap.parse_args()

    if not args.plots_only:
        run_dump(args.workers)

    day_slopes = {}
    with open(SLOPES_CSV, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            day_slopes[r["key"]] = r

    rows = compute_stats(day_slopes)
    out_csv = DUMP_DIR / "perprofile_stats.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"stats -> {out_csv}", flush=True)
    for r in rows:
        print(f"  {r['key']:<22} {r['type']:<5} n={r['n_profiles']:>7} "
              f"prof {r['slope_profile']:+6.2f}+-{r['se_profile']:.2f}  "
              f"win {r['slope_profile_win']:+6.2f}+-{r['se_profile_win']:.2f} "
              f"(out {100 * r['frac_outside_win']:.1f}%)  "
              f"day {r['slope_day']:+6.2f}+-{r['se_day']:.2f}  "
              f"agree2SE={r['agree_2se']}  intra_p90p10={r['intraday_p90p10_km']:.2f} km  "
              f"inter={r['interday_p90p10_km']:.2f} km", flush=True)

    make_plots(rows)


if __name__ == "__main__":
    main()
