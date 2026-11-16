# -*- coding: utf-8 -*-
"""Mechanism of the CHM15k temperature coefficient (report 19, part 4).

Candidate drivers and the fingerprint each one leaves:
  * internal reference (calibration pulse): the firmware divides beta_raw by p_calc, so a reference
    that drifts with temperature makes the signal move with slope -d ln p_calc / dT;
  * laser pulse energy: status_laser moves with T, effect identical at all heights;
  * detector: regulated temperature, gain changes are seen by p_calc and cancel;
  * optics / alignment (thermal expansion): effect strongest in the near range, decays with height;
  * additive electronic baseline (dark): relative effect grows where the signal is weak (far range,
    clean air) and shrinks under strong aerosol signal; night background moves with T;
  * window: window_transmission moves with T.

For each unit, on the clear calibration nights, the same intra-night regression as q2a
(ln x = a + b dT_int + c hours) is applied to
  * every housekeeping variable,
  * the clear-sky signal in every range bin (kappa(z)),
  * the signal UN-normalised by p_calc (rcs_0 x p_calc) in the Rayleigh window,
  * the window signal split by signal level (aerosol load) terciles.

Run:  python validation/chm15k_stability/q4_mechanism.py
"""
from __future__ import annotations

import csv
import json
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import CENSUS, FIGS, OUT, RANGE_MID, chm15k_streams, load_cache, load_constants, stream_key  # noqa: E402
from q2a_intranight import HOUR_BLOCKS, MIN_CLEAR_FRAC, flatten  # noqa: E402

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice")
warnings.filterwarnings("ignore", message="invalid value encountered")
warnings.filterwarnings("ignore", message="divide by zero")

HK_LOG = ("calibration_pulse", "bckgrd_rcs_0", "stddev")              # fitted as ln -> %/K
HK_LIN = ("status_laser", "window_transmission", "status_detector",     # fitted linearly -> unit/K
          "temperature_detector", "temperature_optical_module", "laser_pulses")
MIN_NIGHTS = 25


def joint_fit(y, dT, hrs):
    ok = np.isfinite(y) & np.isfinite(dT) & np.isfinite(hrs)
    if ok.sum() < 40 or np.std(dT[ok]) < 0.2:
        return np.nan, np.nan
    y, dT, hrs = y[ok], dT[ok], hrs[ok]
    med = np.median(y)
    rsd = 1.4826 * np.median(np.abs(y - med))
    keep = np.abs(y - med) <= 4 * max(rsd, 1e-6)
    X = np.column_stack([np.ones(keep.sum()), dT[keep], hrs[keep]])
    coef, *_ = np.linalg.lstsq(X, y[keep], rcond=None)
    res = y[keep] - X @ coef
    cov = (res @ res / max(keep.sum() - 3, 1)) * np.linalg.inv(X.T @ X)
    return coef[1], np.sqrt(cov[1, 1])


def collect(wmo, ident, alt):
    cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
    if cons is None or cache is None:
        return None
    names = list(cache["hk_vars"])
    fl = flatten(cache)
    cos_t = np.cos(np.deg2rad(float(cache["tilt"]) if np.isfinite(cache["tilt"]) else 0.0))
    acc = {k: [] for k in ("dT", "hrs", "win", "win_raw", "level", "prof") + HK_LOG + HK_LIN}
    nn = 0
    for k in range(len(cons["t"])):
        i0 = np.searchsorted(fl["t"], cons["t0"][k])
        i1 = np.searchsorted(fl["t"], cons["t1"][k], side="right")
        nh = (i1 - i0) // HOUR_BLOCKS
        if nh < 4:
            continue
        rmin, rmax = (cons["zb"][k] - alt) / cos_t, (cons["zt"][k] - alt) / cos_t
        win = (RANGE_MID >= rmin) & (RANGE_MID <= rmax)
        low = (RANGE_MID >= 300) & (RANGE_MID <= 1500)
        if win.sum() < 1:
            continue
        rcs = fl["rcs"][i0:i0 + nh * HOUR_BLOCKS].reshape(nh, HOUR_BLOCKS, -1)
        clear = (fl["n_clear"][i0:i0 + nh * HOUR_BLOCKS] >= MIN_CLEAR_FRAC *
                 np.maximum(fl["n_all"][i0:i0 + nh * HOUR_BLOCKS], 1)).reshape(nh, HOUR_BLOCKS)
        rcs = np.where(clear[..., None], rcs, np.nan)
        prof = np.nanmean(rcs, axis=1)                                  # (nh, nbin) hourly profiles
        prof[clear.sum(axis=1) < 4] = np.nan
        hk = np.nanmedian(fl["hk"][i0:i0 + nh * HOUR_BLOCKS].reshape(nh, HOUR_BLOCKS, -1), axis=1)
        T = hk[:, names.index("temp_int")]
        pc = hk[:, names.index("calibration_pulse")]
        s_win = np.nanmean(prof[:, win], axis=1)
        if np.isfinite(s_win).sum() < 4 or np.nanmedian(s_win) <= 0:
            continue
        nn += 1
        tt = fl["t"][i0:i0 + nh * HOUR_BLOCKS].reshape(nh, HOUR_BLOCKS).mean(axis=1)
        acc["dT"].append(T - np.nanmedian(T))
        acc["hrs"].append((tt - tt[0]) * 24)
        acc["win"].append(np.log(s_win / np.nanmedian(s_win)))
        raw = s_win * pc
        acc["win_raw"].append(np.log(raw / np.nanmedian(raw)))
        # aerosol load of the night: low-level signal relative to the window, one value per night
        acc["level"].append(np.full(nh, np.log(np.nanmedian(np.nanmean(prof[:, low], axis=1)) / np.nanmedian(s_win))))
        with np.errstate(divide="ignore", invalid="ignore"):
            acc["prof"].append(np.log(prof / np.nanmedian(prof, axis=0)))
        for v in HK_LOG:
            x = hk[:, names.index(v)]
            acc[v].append(np.log(x / np.nanmedian(x)))
        for v in HK_LIN:
            x = hk[:, names.index(v)]
            acc[v].append(x - np.nanmedian(x))
    if nn < MIN_NIGHTS:
        return None
    return {k: (np.concatenate(v, axis=0) if v else None) for k, v in acc.items()} | {"n_nights": nn}


def pair_background_test():
    """Receiver vs transmitter: in daytime, does the solar-background ratio of two co-located units
    follow their temperature difference like the signal ratio does?

    The background enters through the same receiver (telescope, filter band, detector, electronics)
    but owes nothing to the laser. A receiver gain change moves both ratios alike; a laser-energy or
    laser-line/filter-match change moves only the signal ratio. Anomalies are taken against a
    +-15-day median at the same solar hour, which removes the sun-geometry cycle.
    """
    import itertools
    from common import PAIR_SITES
    from q2b_pairs import MAIN_BAND, hourly, pair_frame
    lon = {(w, i): lo for w, i, s, la, lo in chm15k_streams()}

    def anom_same_hour(fr, x):
        out = np.full_like(x, np.nan)
        hb = np.floor(fr["sh"]).astype(int)
        for h in range(24):
            idx = np.where((hb == h) & np.isfinite(x))[0]
            t = fr["t"][idx]
            for k, i in enumerate(idx):
                m = np.abs(t - t[k]) <= 15
                if m.sum() >= 5:
                    out[i] = x[i] - np.median(x[idx[m]])
        return out

    def fit(x, y):
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 50 or np.std(x[ok]) == 0:
            return np.nan, np.nan, int(ok.sum())
        p, c = np.polyfit(x[ok], y[ok], 1, cov=True)
        return 100 * p[0], 100 * np.sqrt(c[0, 0]), int(ok.sum())

    rows = []
    for wmo, units in PAIR_SITES.items():
        hrs = {u: hourly(load_cache(wmo, u)) for u in units}
        for a, b in itertools.combinations(sorted(hrs), 2):
            fr = pair_frame(hrs[a], hrs[b], lon[(wmo, a)])
            sel = (fr["sh"] >= 9) & (fr["sh"] < 17) & np.isfinite(fr[MAIN_BAND])
            with np.errstate(all="ignore"):
                lnbg = np.where(sel, np.log(fr["bckgrd_rcs_0_A"] / fr["bckgrd_rcs_0_B"]), np.nan)
                lnsig = np.where(sel, fr[MAIN_BAND], np.nan)
                dT = np.where(sel, fr["temp_int_A"] - fr["temp_int_B"], np.nan)
            dTa = anom_same_hour(fr, dT)
            bs, ss, n = fit(dTa, anom_same_hour(fr, lnsig))
            bb, sb, _ = fit(dTa, anom_same_hour(fr, lnbg))
            rows.append(dict(pair=f"{wmo}_{a}/{b}", n=n, dT_sd=float(np.nanstd(dTa)),
                             b_signal=bs, se_signal=ss, b_background=bb, se_background=sb))
            print(f"  {wmo}_{a}/{b}: daytime signal ratio {bs:+.2f}+/-{ss:.2f} %/K, "
                  f"background ratio {bb:+.2f}+/-{sb:.2f} %/K (n {n})", flush=True)
    with open(OUT / "q4_pairs_background.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return rows


def main():
    print("Receiver vs transmitter (pair background test)")
    pair_background_test()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    alts = {(s["wmo"], s["ident"]): float(s["alt"]) for s in json.load(open(CENSUS, encoding="utf-8"))}
    rows, profs = [], []
    for wmo, ident, site, lat, lon in chm15k_streams():
        d = collect(wmo, ident, alts[(wmo, ident)])
        if d is None:
            continue
        dT, hrs = d["dT"], d["hrs"]
        row = dict(stream=stream_key(wmo, ident), site=site, n_nights=d["n_nights"])
        row["b_C"], row["se_C"] = (100 * x for x in joint_fit(d["win"], dT, hrs))
        row["b_C_unnorm"], row["se_C_unnorm"] = (100 * x for x in joint_fit(d["win_raw"], dT, hrs))
        for v in HK_LOG:
            row[f"b_{v}"], row[f"se_{v}"] = (100 * x for x in joint_fit(d[v], dT, hrs))
        for v in HK_LIN:
            row[f"b_{v}"], row[f"se_{v}"] = joint_fit(d[v], dT, hrs)
        # additive vs multiplicative: slope in the clean and in the aerosol-loaded third of nights
        lev = d["level"]
        q1, q2 = np.nanpercentile(lev, [33, 67])
        for name, m in (("clean", lev <= q1), ("loaded", lev >= q2)):
            row[f"b_C_{name}"], row[f"se_C_{name}"] = (100 * x for x in joint_fit(d["win"][m], dT[m], hrs[m]))
        prof = np.array([100 * joint_fit(d["prof"][:, j], dT, hrs)[0] for j in range(len(RANGE_MID))])
        profs.append(prof)
        rows.append(row)
        print(f"  {row['stream']:28s} {site[:16]:16s} n {d['n_nights']:4d}  C {row['b_C']:+5.2f}  "
              f"C_unnorm {row['b_C_unnorm']:+5.2f}  p_cal {row['b_calibration_pulse']:+5.2f}  "
              f"bg {row['b_bckgrd_rcs_0']:+6.2f} %/K  laser {row['b_status_laser']:+.3f} %/K  "
              f"window {row['b_window_transmission']:+.3f} %/K  T_det {row['b_temperature_detector']:+.3f}  "
              f"T_lom {row['b_temperature_optical_module']:+.3f} K/K  clean {row['b_C_clean']:+5.2f} "
              f"loaded {row['b_C_loaded']:+5.2f}", flush=True)
    with open(OUT / "q4_mechanism.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    def med(k):
        x = np.array([r[k] for r in rows], float)
        return np.nanmedian(x), np.nanpercentile(x, 25), np.nanpercentile(x, 75)
    print("\nNetwork medians [IQR] of intra-night slopes (elapsed time held fixed):")
    for k, u in (("b_C", "%/K"), ("b_C_unnorm", "%/K"), ("b_calibration_pulse", "%/K"), ("b_bckgrd_rcs_0", "%/K"),
                 ("b_stddev", "%/K"), ("b_status_laser", "%-pt/K"), ("b_window_transmission", "%-pt/K"),
                 ("b_status_detector", "%-pt/K"), ("b_temperature_detector", "K/K"),
                 ("b_temperature_optical_module", "K/K"), ("b_C_clean", "%/K"), ("b_C_loaded", "%/K")):
        m, lo, hi = med(k)
        print(f"  {k:32s} {m:+7.3f} [{lo:+.3f}, {hi:+.3f}] {u}")
    bc = np.array([r["b_C"] for r in rows], float)
    bp = np.array([r["b_calibration_pulse"] for r in rows], float)
    ok = np.isfinite(bc) & np.isfinite(bp)
    print(f"  per-unit corr(b_C, b_p_cal) = {np.corrcoef(bc[ok], bp[ok])[0, 1]:+.2f}")

    # figure
    P = np.array(profs)
    fig, axs = plt.subplots(1, 3, figsize=(20, 6.2))
    ax = axs[0]
    for p in P:
        ax.plot(p, RANGE_MID / 1000, "-", color="#1f77b4", lw=0.5, alpha=0.25)
    ax.plot(np.nanmedian(P, axis=0), RANGE_MID / 1000, "k-", lw=2.4, label="network median")
    ax.plot(np.nanpercentile(P, 25, axis=0), RANGE_MID / 1000, "k--", lw=1, label="IQR")
    ax.plot(np.nanpercentile(P, 75, axis=0), RANGE_MID / 1000, "k--", lw=1)
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlim(-6, 10)
    ax.set_ylim(0, 8)
    ax.set_xlabel("intra-night d ln(signal) / dT$_{int}$ [%/K]")
    ax.set_ylabel("range [km]")
    ax.set_title("(a) height dependence of the temperature coefficient", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axs[1]
    keys = ["b_C", "b_C_unnorm", "b_calibration_pulse", "b_bckgrd_rcs_0", "b_C_clean", "b_C_loaded"]
    labels = ["C (window)", "C un-normalised\n(rcs x p_cal)", "ln p_cal", "ln night\nbackground", "C, clean\nnights",
              "C, aerosol-\nloaded nights"]
    data = [np.array([r[k] for r in rows], float) for k in keys]
    data = [d_[np.isfinite(d_)] for d_ in data]
    ax.boxplot(data, showfliers=False)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=8)
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_ylabel("%/K")
    ax.set_title("(b) intra-night slopes per unit (network)", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[2]
    keys = ["b_status_laser", "b_window_transmission", "b_status_detector", "b_temperature_detector",
            "b_temperature_optical_module"]
    labels = ["laser energy\n[%-pt/K]", "window\ntransm. [%-pt/K]", "detector\nstatus [%-pt/K]", "T detector\n[K/K]",
              "T optical\nmodule [K/K]"]
    data = [np.array([r[k] for r in rows], float) for k in keys]
    data = [d_[np.isfinite(d_)] for d_ in data]
    ax.boxplot(data, showfliers=False)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=8)
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_title("(c) housekeeping response to T$_{int}$ (intra-night)", fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "q4_mechanism.png", dpi=130)
    print(f"-> {FIGS / 'q4_mechanism.png'}")


if __name__ == "__main__":
    main()
