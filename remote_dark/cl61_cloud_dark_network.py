# -*- coding: utf-8 -*-
"""Dark signal of EVERY network CL61 from opaque low clouds (the "natural hood", fig16), and the check
whether the units with a low Rayleigh/cloud ratio (Uccle 0.79, Payerne 0.83) are the ones with a large
dark.

Per profile (night, sun < -6 deg): first cloud base 300-1500 m, no second layer below 6 km, and an
OPAQUE cloud (signal 300-400 m above the base < 1/20 of the peak). Above an opaque cloud no atmospheric
signal returns, so the recorded signal is the instrument's own additive baseline. Only gates
>= CBH + BUFFER_M are used, to stay clear of the signal-induced undershoot just above the cloud
(fig18: decay ~550 m). Per 100 m gate bin the dark is the median of P = rcs_0 / r^2 over all selected
profiles; uncertainty from a bootstrap over days.

The estimate is NOT rescaled. At Payerne it is 1.67x the hood dark (2.5-7 km), i.e. the strong cloud
return still adds a negative offset even 2 km above the cloud. An extrapolation to zero cloud brightness
(per-gate regression on the cloud charge) was tested and removed (2026-10-09): it moved Payerne only to
1.46x the hood, but was unstable elsewhere (Camborne +2 -> +26 %, Lindenberg -3 -> +18 %, Uccle -20 ->
-5 %) and predicted the observed Rayleigh/cloud ratios worse (Spearman 0.61 vs 0.68, rms 0.54 vs 0.21).

A positive offset is NOT automatically contamination: instrument backgrounds can be positive (Kotthaus
et al., 2016; Wiegner et al., 2019), but above a cloud that is not fully opaque residual atmospheric
signal is positive too, so positive values are ambiguous while negative ones cannot be atmospheric.

The dark is expressed relative to the molecular attenuated backscatter at 910.55 nm (US standard
atmosphere) with the unit's median cloud constant (method_test_2026 run, else the dashboard series);
its mean over 3-6 km (where the Rayleigh fit sits) is compared with the observed C_Rayleigh/C_cloud of
the method test (paper Section 6.1).

Run:  python -m remote_dark.cl61_cloud_dark_network            (from the repo root)
      python -m remote_dark.cl61_cloud_dark_network --replot   (figure from summary.json)
"""
from __future__ import annotations

import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from remote_dark.common import l1_file, solar_elevation_deg   # noqa: E402

D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
OUT = D / "remote_dark" / "cl61_cloud_dark"
FIG = REPO / "doc/reports/figs_remote_dark/cl61_cloud_dark_network.png"
HOOD = D / "rayleigh_availability/dark_profiles_payerne_v3.npz"
METHOD_RUN = D / "method_test_2026"
D0, D1 = datetime(2025, 10, 1), datetime(2026, 10, 8)
CBH_MIN, CBH_MAX = 300.0, 1500.0
BUFFER_M = 2000.0
EDGES = np.arange(2000.0, 12001.0, 100.0)
ZC = 0.5 * (EDGES[1:] + EDGES[:-1])
MIN_PROF = 50
NBOOT = 100
WORKERS = 28
BAND = (ZC >= 3000) & (ZC <= 6000)        # where the Rayleigh fit sits


def census_cl61():
    c = json.loads((REPO / "validation/scope_l1_2026_census.json").read_text())
    return [x for x in c if x["type"] == "CL61"]


def day_job(args):
    """One unit-day -> (CBH, P binned (n x nb)) for the opaque-cloud night profiles, or None."""
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
            times = np.array([datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t])
            el = solar_elevation_deg(lat, lon, times)
            second = (cbh[:, 1] > 0) & (cbh[:, 1] < 6000.0) if cbh.shape[1] > 1 else np.zeros(len(t), bool)
            sel = (el < -6.0) & (cbh[:, 0] >= CBH_MIN) & (cbh[:, 0] <= CBH_MAX) & ~second
            if sel.sum() < 10:
                return None
            rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
            # read the whole variable, then select: boolean indexing on the NetCDF variable itself
            # reads row by row and is prohibitively slow on 10-s CL61 files (Lindenberg)
            rcs = np.ma.filled(ds.variables["rcs_0"][:], np.nan)[sel].astype("f8")
    except Exception:                                                       # noqa: BLE001
        return None
    cb = cbh[sel, 0]
    keep = np.zeros(cb.size, bool)
    for i, c in enumerate(cb):
        m_pk = (rng >= c - 150) & (rng <= c + 150)
        m_ab = (rng >= c + 300) & (rng <= c + 400)
        pk = np.nanmax(rcs[i, m_pk]) if m_pk.any() else np.nan
        keep[i] = np.isfinite(pk) and pk > 0 and np.nanmedian(rcs[i, m_ab]) < pk / 20.0   # opaque
    if keep.sum() < 10:
        return None
    rcs, cb = rcs[keep], cb[keep]
    with np.errstate(all="ignore"):
        P = rcs / np.where(rng > 0, rng, np.nan)[None, :] ** 2
    idx = np.digitize(rng, EDGES) - 1
    Pb = np.full((P.shape[0], ZC.size), np.nan, "f4")
    for b in range(ZC.size):
        m = idx == b
        if m.any():
            with np.errstate(all="ignore"):
                Pb[:, b] = np.nanmean(P[:, m], axis=1)
    Pb[ZC[None, :] < (cb[:, None] + BUFFER_M)] = np.nan          # only well above the cloud
    return cb.astype("f4"), Pb


def collect(unit, pool):
    days, d = [], D0
    while d <= D1:
        days.append(d)
        d += timedelta(days=1)
    jobs = [(unit["wmo"], unit["ident"], unit["lat"], unit["lon"], d) for d in days]
    Cs, Ps, day_id = [], [], []
    for k, r in enumerate(pool.map(day_job, jobs, chunksize=4)):
        if r is None:
            continue
        Cs.append(r[0])
        Ps.append(r[1])
        day_id.append(np.full(r[0].size, k))
    if not Cs:
        return None
    return np.concatenate(Cs), np.vstack(Ps), np.concatenate(day_id)


def dark_median(P):
    """Per gate bin: median of P over the selected profiles (NaN where fewer than MIN_PROF)."""
    n = np.isfinite(P).sum(axis=0)
    with np.errstate(all="ignore"):
        med = np.nanmedian(P, axis=0)
    med[n < MIN_PROF] = np.nan
    return med


def molecular(z, alt, lam=910.55):
    from validation.paper import intercompare as IC
    return IC._molecular_beta(z, alt, lam)


def cloud_constant(key):
    """Median cloud C_L of the unit: method-test run, else the dashboard series."""
    f = METHOD_RUN / key / f"{key}_cal.csv"
    vals = []
    if f.exists():
        for r in csv.DictReader(open(f, encoding="utf-8")):
            if r["method"] == "cloud" and r["flag"] in ("1", "1.0", "0.5") and float(r["cal_value"]) > 0:
                vals.append(float(r["cal_value"]))
    if len(vals) >= 5:
        return float(np.median(vals)), len(vals)
    g = D / "figs_paper_validation/paper_v08/calib_dashboard" / f"{key}_cloud.csv"
    if g.exists():
        v = [float(r["C_daily"]) for r in csv.DictReader(open(g)) if r["C_daily"] not in ("", "nan")]
        if v:
            return float(np.median(v)), len(v)
    return np.nan, 0


def observed_ratio(key):
    s = json.loads((D / "figs_paper_validation/paper_v08/method_comparison_summary.json").read_text())
    r = s.get(key)
    return (r["r_med"], r["r_q1"], r["r_q3"]) if r else (np.nan, np.nan, np.nan)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    units = census_cl61()
    res = {}
    rng_ = np.random.default_rng(3)
    with ProcessPoolExecutor(WORKERS) as pool:
        for u in units:
            key = f"{u['wmo']}_{u['ident']}"
            cache = OUT / f"{key}.npz"
            if cache.exists():
                z = np.load(cache)
                CB, P, day = z["CB"], z["P"], z["day"]
            else:
                got = collect(u, pool)
                if got is None:
                    print(f"{key} {u['site']}: no opaque-cloud night profiles")
                    continue
                CB, P, day = got
                np.savez_compressed(cache, CB=CB, P=P, day=day)
            dark = dark_median(P)
            uniq = np.unique(day)
            groups = [np.where(day == k)[0] for k in uniq]
            boots = np.array([dark_median(P[np.concatenate([groups[i] for i in rng_.integers(0, uniq.size, uniq.size)])])
                              for _ in range(NBOOT)])
            C, nC = cloud_constant(key)
            mol = molecular(ZC, u["alt"])
            to_frac = (ZC ** 2 / C * 1e6 / mol) if np.isfinite(C) else np.full(ZC.size, np.nan)
            frac = dark * to_frac
            frac_sd = np.nanstd(boots, axis=0) * to_frac
            obs = observed_ratio(key)
            # CL61 rcs_0 is manufacturer-calibrated attenuated backscatter (m-1 sr-1): the dark in absolute
            # units needs neither the cloud constant nor a molecular model
            to_beta = ZC ** 2 * 1e6                                           # P -> Mm-1 sr-1
            res[key] = dict(beta=(dark * to_beta).tolist(),
                            dark_beta_3_6km=float(np.nanmean((dark * to_beta)[BAND])),
                            dark_beta_3_6km_sd=float(np.nanstd([np.nanmean((b * to_beta)[BAND]) for b in boots])))
            res[key].update(site=u["site"], n_prof=int(P.shape[0]), n_days=int(uniq.size), C_cloud=C, n_C=nC,
                            dark_3_6km=float(np.nanmean(frac[BAND])),
                            dark_3_6km_sd=float(np.nanstd([np.nanmean((b * to_frac)[BAND]) for b in boots])),
                            obs_ratio=obs, frac=frac.tolist(), frac_sd=frac_sd.tolist(), P=dark.tolist())
            r = res[key]
            print(f"{key:22s} {u['site'][:16]:16s} {r['n_prof']:6d} prof / {r['n_days']:3d} d | dark/mol 3-6 km "
                  f"{100 * r['dark_3_6km']:+6.1f} +/- {100 * r['dark_3_6km_sd']:4.1f} % | "
                  f"observed C_Ray/C_cloud {obs[0]:.2f}", flush=True)
    # Payerne validation against the hood (P view)
    h = np.load(HOOD)
    hz, hP = h["C_range"], h["C_P"]
    hb = np.array([np.nanmean(hP[(hz >= e0) & (hz < e1)]) for e0, e1 in zip(EDGES[:-1], EDGES[1:])])
    pk = "0-20000-0-06610_C"
    if pk in res:
        v = np.array(res[pk]["P"])
        ok = (ZC >= 2500) & (ZC <= 7000) & np.isfinite(v) & np.isfinite(hb)
        res["payerne_vs_hood"] = float(np.sum(v[ok] * hb[ok]) / np.sum(hb[ok] ** 2))
        res["payerne_hood_P"] = hb.tolist()
        print(f"Payerne CL61 low-cloud dark / hood dark (2.5-7 km): {res['payerne_vs_hood']:.2f}")
    res["z"] = ZC.tolist()
    (OUT / "summary.json").write_text(json.dumps(res, indent=1, default=float))
    figure(res)


def figure(res):
    """(a) dark of every CL61 from opaque low clouds, in absolute attenuated backscatter (the CL61 rcs_0 is
    manufacturer-calibrated, m-1 sr-1: no cloud constant, no molecular model); (b) its mean over 3-6 km vs
    the observed Rayleigh-to-cloud ratio of the method test."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    z = np.array(res["z"]) / 1e3
    keys = sorted([k for k in res if k.startswith("0-") and isinstance(res[k], dict)
                   and np.isfinite(res[k]["dark_beta_3_6km"])], key=lambda k: res[k]["site"])
    cmap = plt.get_cmap("tab20")
    col = {k: cmap((2 * i) % 20 + (i // 10)) for i, k in enumerate(keys)}

    def sm(x):
        """500 m running mean, NaN-aware (no edge fall-off to zero); display only."""
        x = np.asarray(x, float)
        ok = np.isfinite(x)
        num = np.convolve(np.where(ok, x, 0.0), np.ones(5), mode="same")
        den = np.convolve(ok.astype(float), np.ones(5), mode="same")
        with np.errstate(all="ignore"):
            return np.where(den > 0, num / den, np.nan)

    fig, axs = plt.subplots(1, 2, figsize=(17, 7), gridspec_kw=dict(width_ratios=[1, 1]))
    ax = axs[0]
    for k in keys:
        lw = 2.8 if ("06447" in k or "06610" in k) else 1.3
        ax.plot(sm(res[k]["beta"]), z, "-", color=col[k], lw=lw, label=res[k]["site"].title()[:16])
    ax.axvline(0, color="k", lw=0.7)
    ax.axhspan(3, 6, color="0.9", zorder=0)
    ax.set_xlim(-0.08, 0.10)
    ax.set_ylim(2, 10)
    ax.set_xlabel(r"signal above opaque low clouds, attenuated backscatter  [Mm$^{-1}$ sr$^{-1}$]")
    ax.set_ylabel("Altitude a.g.l. [km]")
    ax.set_title("(a) Dark of every CL61 from opaque low clouds (night, gates >= 2 km above cloud base;\n"
                 "median, 500 m running mean); grey: Rayleigh fit range", fontsize=10)
    ax.legend(fontsize=8, loc="lower right", ncol=1)
    ax.grid(alpha=0.3)
    ax = axs[1]
    xs, ys = [], []
    for k in keys:
        r = res[k]
        o = r["obs_ratio"]
        if not np.isfinite(o[0]):
            continue
        x = r["dark_beta_3_6km"]
        ax.errorbar(x, o[0], xerr=r["dark_beta_3_6km_sd"], yerr=[[o[0] - o[1]], [o[2] - o[0]]], fmt="o", ms=9,
                    capsize=3, color=col[k])
        ax.annotate(r["site"].title()[:12], (x, o[0]), fontsize=8.5, xytext=(5, 4), textcoords="offset points")
        xs.append(x)
        ys.append(o[0])
    rho, p = spearmanr(xs, ys)
    ax.axhline(1, color="0.5", lw=0.7)
    ax.axvline(0, color="0.5", lw=0.7)
    ax.set_xlabel(r"signal above opaque low clouds, mean over 3-6 km  [Mm$^{-1}$ sr$^{-1}$]")
    ax.set_ylabel(r"observed C$_{L,Rayleigh}$ / C$_{L,cloud}$ (Feb-Jun 2026, median and IQR)")
    ax.set_title(f"(b) Dark vs observed Rayleigh-to-cloud ratio of each CL61\n(Spearman rho = {rho:.2f}, "
                 f"p = {p:.2f}, n = {len(xs)})", fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=200)
    print(f"-> {FIG}  (Spearman {rho:.2f}, p {p:.3f})")


def replot():
    figure(json.loads((OUT / "summary.json").read_text()))


if __name__ == "__main__":
    replot() if "--replot" in sys.argv else main()
