# -*- coding: utf-8 -*-
"""Why each calibration method works or fails, per instrument type (illustration for Section 6.1).

One representative unit per type, night-time L1 data, 1 February - 29 June 2026:
  top row    (Rayleigh)  clear-sky profiles (no cloud detected below 10 km): range-corrected signal
             divided by the molecular attenuated backscatter (US standard atmosphere), normalised to
             its median over 3-5 km; median and IQR of hourly medians. A usable Rayleigh calibration
             needs a flat ratio with small spread in the free troposphere.
  bottom row (cloud)     opaque low clouds (first base 500-2000 m, no second layer): signal in the
             cloud-relative frame (height above cloud base), normalised by each profile's peak; median
             and IQR. The liquid-cloud method needs a sharp, unsaturated peak that drops by >= 20x within
             300 m above it (the screening test of Section 5.1.2).

Run:  python validation/paper/fig_method_illustration.py
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from remote_dark.common import l1_file, solar_elevation_deg   # noqa: E402

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v08")
CACHE = OUT / "method_illustration_cache.npz"
UNITS = [("CL31", "0-20000-0-06610_B"), ("CL51", "0-20000-0-11035_A"), ("CL61", "0-20000-0-10393_C"),
         ("CHM15k", "0-20000-0-10393_0"), ("Mini-MPL", "0-20000-0-07774_A")]
TCOL = {"CL31": "#1b9e77", "CL51": "#66a61e", "CL61": "#7570b3", "CHM15k": "#d95f02", "Mini-MPL": "#e7298a"}
D0, D1 = datetime(2026, 2, 1), datetime(2026, 6, 29)
ZCLR = np.arange(250.0, 10001.0, 100.0)            # bin centres, clear-sky row
ZREL = np.arange(-300.0, 1201.0, 15.0)             # height above cloud base, cloud row


def unit_meta(key):
    c = json.loads((REPO / "validation/scope_l1_2026_census.json").read_text())
    for x in c:
        if f"{x['wmo']}_{x['ident']}" == key:
            return x
    raise KeyError(key)


def day_job(args):
    wmo, ident, lat, lon, day = args
    f = l1_file(wmo, ident, day)
    if not f.exists():
        return None
    try:
        with netCDF4.Dataset(f) as ds:
            t = np.ma.filled(ds.variables["time"][:].astype("f8"), np.nan)
            cbh = np.ma.filled(ds.variables["cloud_base_height"][:].astype("f8"), np.nan)
            if cbh.ndim == 1:
                cbh = cbh[:, None]
            rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
            rcs = np.ma.filled(ds.variables["rcs_0"][:].astype("f8"), np.nan)
    except Exception:                                                      # noqa: BLE001
        return None
    times = np.array([datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t])
    night = solar_elevation_deg(lat, lon, times) < -6.0
    anycloud = np.any((cbh > 0) & (cbh < 10000), axis=1)
    out = {}
    # clear sky: hourly medians of the cloud-free night profiles, binned to 100 m
    clr = night & ~anycloud
    hours = np.array([x.hour for x in times])
    prof = []
    for h in np.unique(hours[clr]):
        m = clr & (hours == h)
        if m.sum() < 10:
            continue
        p = np.nanmedian(rcs[m], axis=0)
        prof.append([np.nanmean(p[(rng >= z - 50) & (rng < z + 50)]) for z in ZCLR])
    out["clear"] = np.array(prof, "f4") if prof else None
    # opaque low cloud, cloud-relative frame, peak-normalised
    second = (cbh[:, 1] > 0) & (cbh[:, 1] < 6000) if cbh.shape[1] > 1 else np.zeros(t.size, bool)
    cl = night & (cbh[:, 0] >= 500) & (cbh[:, 0] <= 2000) & ~second
    cprof = []
    for i in np.where(cl)[0][::3]:                                    # every 3rd profile is plenty
        # align on the cloud PEAK (searched -300..+600 m around the reported base): the Mini-MPL L1
        # (5-min averages, 30 m gates) reports a base well below its smeared peak
        win = (rng >= cbh[i, 0] - 300) & (rng <= cbh[i, 0] + 600)
        if not win.any() or not np.isfinite(rcs[i, win]).any():
            continue
        j = np.where(win)[0][np.nanargmax(rcs[i, win])]
        pk = rcs[i, j]
        if not (np.isfinite(pk) and pk > 0):
            continue
        cprof.append(np.interp(ZREL, rng - rng[j], rcs[i] / pk, left=np.nan, right=np.nan))
    out["cloud"] = np.array(cprof, "f4") if cprof else None
    return out


LIND = ("0-20000-0-10393", "C", "0")         # Lindenberg: CL61 (ident C) and CHM15k (ident 0)
SAT_CACHE = OUT / "method_illustration_saturation.npz"


def _minute_medians(t, x):
    m = np.floor(t * 1440.0)
    keys, inv = np.unique(m, return_inverse=True)
    out = np.full((keys.size,) + x.shape[1:], np.nan)
    for k in range(keys.size):
        with np.errstate(all="ignore"):
            out[k] = np.nanmedian(x[inv == k], axis=0)
    return keys, out


def lind_day(day):
    """Time-matched (1-min) integrated cloud backscatter of the Lindenberg CL61 and CHM15k in opaque
    low clouds (opacity judged on the CL61)."""
    got = {}
    for ident in (LIND[1], LIND[2]):
        f = l1_file(LIND[0], ident, day)
        if not f.exists():
            return None
        try:
            with netCDF4.Dataset(f) as ds:
                t = np.ma.filled(ds.variables["time"][:].astype("f8"), np.nan)
                rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
                rcs = np.ma.filled(ds.variables["rcs_0"][:].astype("f4"), np.nan)
                cbh = np.ma.filled(ds.variables["cloud_base_height"][:].astype("f8"), np.nan)
        except Exception:                                                  # noqa: BLE001
            return None
        cbh = cbh if cbh.ndim == 2 else cbh[:, None]
        sel = (cbh[:, 0] >= 300) & (cbh[:, 0] <= 2000)
        if sel.sum() < 30:
            return None
        km, prof = _minute_medians(t[sel], rcs[sel])
        _, cb = _minute_medians(t[sel], cbh[sel, :2])
        got[ident] = (km, prof, cb, rng)
    k61, p61, c61, r61 = got[LIND[1]]
    kch, pch, _, rch = got[LIND[2]]
    common, i61, ich = np.intersect1d(k61, kch, return_indices=True)
    B61, Bch, CB = [], [], []
    for a, b in zip(i61, ich):
        c = c61[a, 0]
        if not np.isfinite(c) or (np.isfinite(c61[a, 1]) and 0 < c61[a, 1] < 6000):
            continue
        y = p61[a]
        pk = np.nanmax(y[(r61 >= c - 150) & (r61 <= c + 150)])
        above = np.nanmedian(y[(r61 >= c + 300) & (r61 <= c + 400)])
        if not (np.isfinite(pk) and pk > 0 and above < pk / 20.0):               # opaque for the CL61
            continue
        m61 = (r61 >= c - 100) & (r61 <= c + 400)
        mch = (rch >= c - 100) & (rch <= c + 400)
        B61.append(np.nansum(y[m61]) * np.median(np.diff(r61)))
        Bch.append(np.nansum(pch[b, mch]) * np.median(np.diff(rch)))
        CB.append(c)
    return np.array(B61), np.array(Bch), np.array(CB)


def saturation():
    if SAT_CACHE.exists():
        z = np.load(SAT_CACHE)
        return z["B61"], z["Bch"], z["CB"]
    days = [D0 + timedelta(days=k) for k in range((D1 - D0).days + 1)]
    A, B, C = [], [], []
    with ThreadPoolExecutor(16) as pool:
        for r in pool.map(lind_day, days):
            if r is not None and r[0].size:
                A.append(r[0])
                B.append(r[1])
                C.append(r[2])
    B61, Bch, CB = np.concatenate(A), np.concatenate(B), np.concatenate(C)
    np.savez_compressed(SAT_CACHE, B61=B61, Bch=Bch, CB=CB)
    print(f"Lindenberg co-located opaque-cloud minutes: {B61.size}")
    return B61, Bch, CB


def median_constant(key, method):
    import csv
    f = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/method_test_2026") / key / f"{key}_cal.csv"
    v = [float(r["cal_value"]) for r in csv.DictReader(open(f, encoding="utf-8"))
         if r["method"] == method and r["flag"] in ("1", "1.0", "0.5") and float(r["cal_value"]) > 0]
    return float(np.median(v))


def collect():
    data = {}
    with ThreadPoolExecutor(16) as pool:   # threads: a process pool broke on Windows
        for typ, key in UNITS:
            u = unit_meta(key)
            days = [D0 + timedelta(days=k) for k in range((D1 - D0).days + 1)]
            jobs = [(u["wmo"], u["ident"], u["lat"], u["lon"], d) for d in days]
            C, K = [], []
            for r in pool.map(day_job, jobs, chunksize=4):
                if r is None:
                    continue
                if r["clear"] is not None:
                    C.append(r["clear"])
                if r["cloud"] is not None:
                    K.append(r["cloud"])
            data[f"{typ}_clear"] = np.vstack(C) if C else np.empty((0, ZCLR.size))
            data[f"{typ}_cloud"] = np.vstack(K) if K else np.empty((0, ZREL.size))
            print(f"{typ} {key}: {data[f'{typ}_clear'].shape[0]} clear night hours, "
                  f"{data[f'{typ}_cloud'].shape[0]} opaque-cloud profiles")
    np.savez_compressed(CACHE, **data)
    return data


def main():
    data = dict(np.load(CACHE)) if CACHE.exists() else collect()
    from validation.paper import intercompare as IC
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(21, 14))
    gs = fig.add_gridspec(3, 5, height_ratios=[1, 1, 0.85], hspace=0.45)
    axs = np.array([[fig.add_subplot(gs[r, c]) for c in range(5)] for r in range(2)])
    verdict_r = {"CL31": "Rayleigh: NOT usable\n(noise > molecular signal)",
                 "CL51": "Rayleigh: biased\n(signal artefacts, ratio drifts)",
                 "CL61": "Rayleigh: usable\n(dark signal permitting)",
                 "CHM15k": "Rayleigh: usable", "Mini-MPL": "Rayleigh: usable"}
    verdict_c = {"CL31": "cloud: usable", "CL51": "cloud: usable", "CL61": "cloud: usable",
                 "CHM15k": "cloud: NOT usable\n(saturated: amplitude low, see bottom)",
                 "Mini-MPL": "cloud: NOT usable\n(no sharp drop above peak)"}
    for j, (typ, key) in enumerate(UNITS):
        u = unit_meta(key)
        col = TCOL[typ]
        lam = 1064.0 if typ == "CHM15k" else (532.0 if typ == "Mini-MPL" else 910.0)
        mol = IC._molecular_beta(ZCLR, u["alt"], lam)
        ax = axs[0, j]
        C = data[f"{typ}_clear"]
        if C.shape[0]:
            R = C / mol[None, :]
            ref = np.nanmedian(R[:, (ZCLR >= 3000) & (ZCLR <= 5000)], axis=1)
            R = R / ref[:, None]
            with np.errstate(all="ignore"):
                med, q1, q3 = (np.nanmedian(R, axis=0), np.nanpercentile(R, 25, axis=0),
                               np.nanpercentile(R, 75, axis=0))
            ax.fill_betweenx(ZCLR / 1e3, q1, q3, color=col, alpha=0.25, lw=0)
            ax.plot(med, ZCLR / 1e3, color=col, lw=2)
        ax.axvline(1, color="k", lw=0.9)
        ax.axhspan(3, 5, color="0.9", zorder=0)
        ax.set_xlim(-1, 3)
        ax.set_ylim(0.25, 10)
        ax.set_title(f"{typ} ({u['site'].title()})\n{verdict_r[typ]}", fontsize=10.5, color=col)
        ax.set_xlabel("signal / molecular (norm. 3-5 km)")
        ax.grid(alpha=0.3)
        ax = axs[1, j]
        K = data[f"{typ}_cloud"]
        if K.shape[0]:
            with np.errstate(all="ignore"):
                med, q1, q3 = (np.nanmedian(K, axis=0), np.nanpercentile(K, 25, axis=0),
                               np.nanpercentile(K, 75, axis=0))
            ax.fill_betweenx(ZREL, q1, q3, color=col, alpha=0.25, lw=0)
            ax.plot(med, ZREL, color=col, lw=2)
        ax.axvline(1 / 20, color="k", ls="--", lw=0.9)
        ax.axvline(0, color="k", lw=0.6)
        ax.axhspan(300, 400, color="0.85", zorder=0)
        ax.set_xscale("symlog", linthresh=1e-2)
        ax.set_xlim(-0.05, 2)
        ax.set_ylim(-300, 1200)
        ax.set_title(f"{verdict_c[typ]}  (n = {K.shape[0]})", fontsize=9.5, color=col)
        ax.set_xlabel("signal / peak signal (symlog)")
        ax.grid(alpha=0.3)
    axs[0, 0].set_ylabel("Clear night: altitude a.g.l. [km]")
    axs[1, 0].set_ylabel("Opaque cloud: height above the cloud peak [m]")
    # third row: CHM15k saturation, co-located with a CL61 at Lindenberg
    B61, Bch, CB = saturation()
    c61 = median_constant("0-20000-0-10393_C", "cloud")
    cch = median_constant("0-20000-0-10393_0", "rayleigh")
    b61, bch = B61 / c61, Bch / cch
    ok = (b61 > 0) & (bch > 0)
    ax = fig.add_subplot(gs[2, 0:3])
    hb = ax.hexbin(b61[ok], bch[ok] / b61[ok], xscale="log", yscale="log", gridsize=60, mincnt=3, cmap="Oranges",
                   bins="log", extent=(np.log10(3e-3), np.log10(0.2), np.log10(0.1), np.log10(3)))
    edges = np.logspace(np.log10(3e-3), np.log10(0.2), 16)
    mid = np.sqrt(edges[1:] * edges[:-1])
    med = [np.median((bch / b61)[ok & (b61 >= a) & (b61 < b)]) if np.sum(ok & (b61 >= a) & (b61 < b)) > 20 else np.nan
           for a, b in zip(edges[:-1], edges[1:])]
    ax.plot(mid, med, "o-", color="#d95f02", lw=2, label="median per brightness bin")
    ax.axhline(1, color="k", lw=0.9)
    ax.axvline(1 / (2 * 18.8), color="0.4", ls="--", lw=1, label="opaque liquid cloud, 1/(2 S), S = 18.8 sr")
    ax.set_xlabel(r"CL61 integrated cloud backscatter B [sr$^{-1}$] (cloud-calibrated)")
    ax.set_ylabel(r"B$_{CHM15k}$ / B$_{CL61}$  (CHM15k Rayleigh-calibrated)")
    ax.set_title(f"CHM15k saturation (Lindenberg, co-located CL61): {ok.sum()} time-matched 1-min opaque "
                 "low clouds, B integrated over base-100 .. base+400 m", fontsize=10.5, color=TCOL["CHM15k"],
                 loc="left")
    ax.legend(fontsize=9, loc="lower left")
    ax.grid(alpha=0.3, which="both")
    fig.colorbar(hb, ax=ax, label="count", pad=0.01)
    print(f"saturation: n={ok.sum()}, ratio medians per bin: {np.round(med, 2)}")
    ax = fig.add_subplot(gs[2, 3:5])
    ax.axis("off")
    ax.text(0, 0.95, "Reading the figure", fontsize=11, fontweight="bold", va="top")
    ax.text(0, 0.8, "Top: a usable Rayleigh calibration needs signal/molecular flat, with small\n"
            "spread, around the fit range. CL31: noise everywhere above ~2 km.\n"
            "CL51: the ratio drifts with altitude (instrument artefacts).\n\n"
            "Middle: the cloud method needs a sharp peak that drops >= 20x within\n"
            "300-400 m. The Mini-MPL signal stays high far above the peak: its L1\n"
            "(5-min averages, 30 m gates) smears the cloud, and detector recovery\n"
            "may add to it - the two cannot be separated with these data.\n\n"
            "Bottom: in the same clouds a linear CHM15k would give B_CHM/B_CL61 ~ 1\n"
            "at every brightness; a ratio below 1 that falls as the cloud gets brighter\n"
            "is the signature of photon-counting saturation.", fontsize=9.5, va="top", family="monospace")
    fig.suptitle("Why each calibration method works or fails, by instrument type (night-time L1, 1 Feb - 29 Jun 2026)."
                 "  Top: clear-sky signal relative to the molecular profile (median, IQR of hourly profiles; grey: "
                 "3-5 km normalisation).\nMiddle: opaque low clouds, signal relative to its peak (median, IQR); the "
                 "dashed line is the 1/20 drop required 300-400 m above the base (grey band) by the cloud screening. "
                 "Bottom: CHM15k vs co-located CL61 in the same clouds.",
                 fontsize=11, y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "fig_v08_method_illustration.png", dpi=160, bbox_inches="tight")
    print(f"-> {OUT / 'fig_v08_method_illustration.png'}")


if __name__ == "__main__":
    main()
