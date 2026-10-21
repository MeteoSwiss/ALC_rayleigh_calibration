# -*- coding: utf-8 -*-
"""Q2a of report 19: is the CHM15k Rayleigh constant stable WITHIN a night?

For every operational calibration night, the clear-sky signal in that night's own molecular window
(calibration_bottom/top_height, frozen for the whole night) is followed hour by hour from the
10-min cache. Over one night the molecular backscatter in the window changes by well under 1 %
(T/p at 3-6 km), so the hourly signal ratio r(h) = S(h) / median_night(S) is a direct measure of the
relative lidar constant C(h) / C_night.

Per unit we report
  * the intra-night slope d ln r / dT_int (night cooling of the instrument, a few K),
  * the drift with time since the start of the night,
  * the excess intra-night scatter above the photon-noise expectation (estimated from the
    10-min block-to-block scatter), i.e. the part of the hourly variability that is real.

Run:  python validation/chm15k_stability/q2a_intranight.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import warnings

import numpy as np

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (BLOCKS_PER_DAY, CENSUS, FIGS, OUT, RANGE_MID, chm15k_streams, load_cache,  # noqa: E402
                    load_constants, stream_key)

HOUR_BLOCKS = 6                 # 6 x 10 min
MIN_CLEAR_FRAC = 0.9            # a block enters only if >= 90 % of its profiles are clear
MIN_HOURS = 4                   # nights with fewer usable hours are skipped
MIN_NIGHTS = 15


def station_alt():
    return {(s["wmo"], s["ident"]): float(s["alt"]) for s in json.load(open(CENSUS, encoding="utf-8"))}


def flatten(cache):
    """Absolute block times (days since 1970) and the per-block arrays, time-sorted."""
    tb = (cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY).ravel()
    order = np.argsort(tb)
    nb = cache["rcs"].shape[-1]
    return dict(t=tb[order],
                rcs=cache["rcs"].reshape(-1, nb)[order],
                hk=cache["hk"].reshape(-1, cache["hk"].shape[-1])[order],
                n_all=cache["n_all"].ravel()[order].astype(float),
                n_clear=cache["n_clear"].ravel()[order].astype(float))


def night_series(fl, names, t0, t1, rmin, rmax):
    """Hourly window signal, T_int and noise estimate over one calibration night."""
    i0 = np.searchsorted(fl["t"], t0)
    i1 = np.searchsorted(fl["t"], t1, side="right")
    if i1 - i0 < MIN_HOURS * HOUR_BLOCKS:
        return None
    sel = (RANGE_MID >= rmin) & (RANGE_MID <= rmax)
    if sel.sum() < 2:
        return None
    clear = fl["n_clear"][i0:i1] >= MIN_CLEAR_FRAC * np.maximum(fl["n_all"][i0:i1], 1)
    with np.errstate(invalid="ignore"):
        s = np.nanmean(fl["rcs"][i0:i1][:, sel], axis=1)
    s[~clear] = np.nan
    tint = fl["hk"][i0:i1, names.index("temp_int")]
    tt = fl["t"][i0:i1]
    nh = (i1 - i0) // HOUR_BLOCKS
    hs, ht, hT, hn = [], [], [], []
    for h in range(nh):
        seg = slice(h * HOUR_BLOCKS, (h + 1) * HOUR_BLOCKS)
        v = s[seg]
        v = v[np.isfinite(v)]
        if v.size < 4:
            continue
        hs.append(np.mean(v))
        # photon-noise estimate of the hourly mean from the 10-min scatter
        hn.append(np.std(v, ddof=1) / np.sqrt(v.size) / np.mean(v) if np.mean(v) > 0 else np.nan)
        ht.append(np.mean(tt[seg]))
        hT.append(np.nanmedian(tint[seg]))
    if len(hs) < MIN_HOURS:
        return None
    hs = np.array(hs)
    if np.median(hs) <= 0:
        return None
    lr = np.log(np.clip(hs / np.median(hs), 1e-3, None))
    return dict(lr=lr, t=np.array(ht), T=np.array(hT), noise=np.array(hn),
                hours=(np.array(ht) - ht[0]) * 24)


def analyse_stream(wmo, ident, alt):
    cons = load_constants(wmo, ident)
    cache = load_cache(wmo, ident)
    if cons is None or cache is None:
        return None, None
    names = list(cache["hk_vars"])
    fl = flatten(cache)
    cos_tilt = np.cos(np.deg2rad(float(cache["tilt"]) if np.isfinite(cache["tilt"]) else 0.0))
    nights = []
    for k in range(len(cons["t"])):
        rmin = (cons["zb"][k] - alt) / cos_tilt
        rmax = (cons["zt"][k] - alt) / cos_tilt
        if not (np.isfinite(rmin) and np.isfinite(rmax)):
            continue
        ns = night_series(fl, names, cons["t0"][k], cons["t1"][k], rmin, rmax)
        if ns is not None:
            ns["module"] = cons["module"][k]
            nights.append(ns)
    if len(nights) < MIN_NIGHTS:
        return None, nights
    lr = np.concatenate([n["lr"] for n in nights])
    dT = np.concatenate([n["T"] - np.nanmedian(n["T"]) for n in nights])
    hrs = np.concatenate([n["hours"] for n in nights])
    noise = np.concatenate([n["noise"] for n in nights])
    # Robust screen: thin cloud / aerosol layers that slip past the cloud-base screen show up as
    # +-50 % hours; clip at 4 robust SDs so they neither drive the slope nor the scatter.
    fin = np.isfinite(lr) & np.isfinite(dT)
    rsd = 1.4826 * np.median(np.abs(lr[fin] - np.median(lr[fin])))
    ok = fin & (np.abs(lr - np.median(lr[fin])) <= 4 * max(rsd, 1e-3))
    out = dict(stream=stream_key(wmo, ident), n_nights=len(nights), n_hours=int(ok.sum()),
               dT_std=float(np.std(dT[ok])))
    if np.std(dT[ok]) > 0.2:
        p, cov = np.polyfit(dT[ok], lr[ok], 1, cov=True)
        out.update(b_T=100 * p[0], se_T=100 * np.sqrt(cov[0, 0]))
    else:
        out.update(b_T=np.nan, se_T=np.nan)
    p, cov = np.polyfit(hrs[ok], lr[ok], 1, cov=True)
    out.update(b_hour=100 * p[0], se_hour=100 * np.sqrt(cov[0, 0]))
    # Joint fit: the instrument cools steadily through the night, so T and elapsed time are
    # collinear, and a slow atmospheric evolution in the window (aerosol settling) would masquerade as
    # a temperature effect. Nights cool at different rates (wind, cloud before/after), which makes
    # the partial T slope identifiable with elapsed time held fixed.
    X = np.column_stack([np.ones(ok.sum()), dT[ok], hrs[ok]])
    coef, *_ = np.linalg.lstsq(X, lr[ok], rcond=None)
    res = lr[ok] - X @ coef
    cov = (res @ res / max(ok.sum() - 3, 1)) * np.linalg.inv(X.T @ X)
    out.update(b_T_ctrl=100 * coef[1], se_T_ctrl=100 * np.sqrt(cov[1, 1]),
               b_hour_ctrl=100 * coef[2], r_T_hour=float(np.corrcoef(dT[ok], hrs[ok])[0, 1]))
    # Excess scatter: observed hourly scatter vs median photon-noise expectation.
    obs = 1.4826 * np.median(np.abs(lr[ok] - np.median(lr[ok])))
    exp_ = np.sqrt(np.nanmedian(noise[ok] ** 2))
    out.update(sd_obs=100 * obs, sd_noise=100 * exp_, frac_clipped=float(1 - ok.sum() / max(fin.sum(), 1)),
               sd_excess=100 * np.sqrt(max(obs ** 2 - exp_ ** 2, 0.0)))
    return out, nights


def figures(rows, examples):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIGS.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(1, 3, figsize=(20, 5.2))
    b = np.array([r["b_T"] for r in rows], float)
    se = np.array([r["se_T"] for r in rows], float)
    ok = np.isfinite(b)
    axs[0].hist(np.clip(b[ok], -15, 15), bins=np.linspace(-15, 15, 41), color="#1f77b4")
    w = 1 / se[ok] ** 2
    pooled = np.sum(w * b[ok]) / np.sum(w)
    axs[0].axvline(0, color="k", lw=0.8)
    axs[0].set_title(f"intra-night d ln C / dT$_{{int}}$ per unit\npooled {pooled:+.2f} ± "
                     f"{1 / np.sqrt(np.sum(w)):.2f} %/K, N={ok.sum()}", fontsize=10)
    axs[0].set_xlabel("%/K")
    hb = np.array([r["b_hour"] for r in rows], float)
    axs[1].hist(np.clip(hb, -5, 5), bins=np.linspace(-5, 5, 41), color="#2ca02c")
    axs[1].axvline(0, color="k", lw=0.8)
    axs[1].set_title(f"drift through the night per unit\nmedian {np.nanmedian(hb):+.2f} %/h", fontsize=10)
    axs[1].set_xlabel("%/h")
    so = np.array([r["sd_obs"] for r in rows], float)
    sn = np.array([r["sd_noise"] for r in rows], float)
    axs[2].scatter(sn, so, s=14)
    mx = np.nanmax(so) * 1.05
    axs[2].plot([0, mx], [0, mx], "k--", lw=1, label="pure noise")
    axs[2].set_xlabel("expected hourly photon noise [%]")
    axs[2].set_ylabel("observed hourly scatter of C(h)/C$_{night}$ [%]")
    axs[2].set_title(f"intra-night scatter vs noise (median excess {np.nanmedian([r['sd_excess'] for r in rows]):.1f} %)",
                     fontsize=10)
    axs[2].legend(fontsize=8)
    for ax in axs:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "q2a_intranight_summary.png", dpi=130)
    plt.close(fig)

    if examples:
        fig, axs = plt.subplots(1, len(examples), figsize=(4.4 * len(examples), 4.6), squeeze=False)
        for ax, (key, nights) in zip(axs[0], examples):
            for n in nights:
                ax.plot(n["T"] - np.nanmedian(n["T"]), 100 * n["lr"], ".", ms=3, alpha=0.4, color="#1f77b4")
            dT = np.concatenate([n["T"] - np.nanmedian(n["T"]) for n in nights])
            lr = np.concatenate([n["lr"] for n in nights])
            ok = np.isfinite(dT) & np.isfinite(lr) & (np.abs(lr) < 0.7)
            if ok.sum() > 10:
                edges = np.percentile(dT[ok], np.linspace(0, 100, 9))
                mid, med = [], []
                for a, b2 in zip(edges[:-1], edges[1:]):
                    m = ok & (dT >= a) & (dT <= b2)
                    if m.sum() > 5:
                        mid.append(np.median(dT[m]))
                        med.append(100 * np.median(lr[m]))
                ax.plot(mid, med, "o-", color="#d62728", lw=2, label="binned median")
            ax.axhline(0, color="k", lw=0.8)
            ax.set_ylim(-30, 30)
            ax.set_xlabel("T$_{int}$ - night median [K]")
            ax.set_ylabel("C(h)/C$_{night}$ - 1 [%]")
            ax.set_title(f"{key} ({len(nights)} nights)", fontsize=9)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIGS / "q2a_intranight_examples.png", dpi=130)
        plt.close(fig)


def main():
    alts = station_alt()
    rows, examples = [], []
    want = {"0-20000-0-06610_A", "0-20000-0-10393_0", "0-20000-0-06240_A", "0-20008-0-INO_B", "0-20000-0-12843_A"}
    for wmo, ident, site, lat, lon in chm15k_streams():
        res, nights = analyse_stream(wmo, ident, alts[(wmo, ident)])
        if res is None:
            continue
        res["site"] = site
        rows.append(res)
        if stream_key(wmo, ident) in want:
            examples.append((f"{site} {ident}", nights))
        print(f"  {res['stream']:28s} {site[:18]:18s} nights {res['n_nights']:4d}  dT sd {res['dT_std']:.2f} K  "
              f"b_T {res['b_T']:+6.2f}+/-{res['se_T']:.2f} %/K  drift {res['b_hour']:+.2f} %/h  "
              f"| joint b_T {res['b_T_ctrl']:+6.2f}+/-{res['se_T_ctrl']:.2f} drift {res['b_hour_ctrl']:+.2f} "
              f"r(T,h) {res['r_T_hour']:+.2f} |  "
              f"scatter {res['sd_obs']:.1f} % (noise {res['sd_noise']:.1f}, excess {res['sd_excess']:.1f})", flush=True)
    with open(OUT / "q2a_intranight.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    b = np.array([r["b_T"] for r in rows], float)
    se = np.array([r["se_T"] for r in rows], float)
    ok = np.isfinite(b) & np.isfinite(se)
    wts = 1 / se[ok] ** 2
    print(f"\nIntra-night T slope: median {np.median(b[ok]):+.2f} %/K, pooled {np.sum(wts * b[ok]) / np.sum(wts):+.2f}"
          f" +/- {1 / np.sqrt(np.sum(wts)):.2f}, sig +{np.sum(b[ok] > 3 * se[ok])} / -{np.sum(b[ok] < -3 * se[ok])}"
          f" of {ok.sum()}")
    bc = np.array([r["b_T_ctrl"] for r in rows], float)
    sc = np.array([r["se_T_ctrl"] for r in rows], float)
    okc = np.isfinite(bc) & np.isfinite(sc)
    wc = 1 / sc[okc] ** 2
    print(f"Intra-night T slope, elapsed time held fixed: median {np.median(bc[okc]):+.2f} %/K, pooled "
          f"{np.sum(wc * bc[okc]) / np.sum(wc):+.2f} +/- {1 / np.sqrt(np.sum(wc)):.2f}, sig +{np.sum(bc[okc] > 3 * sc[okc])}"
          f" / -{np.sum(bc[okc] < -3 * sc[okc])} of {okc.sum()}")
    print(f"Drift: median {np.nanmedian([r['b_hour'] for r in rows]):+.2f} %/h; "
          f"excess scatter median {np.nanmedian([r['sd_excess'] for r in rows]):.1f} %")
    figures(rows, examples)


if __name__ == "__main__":
    main()
