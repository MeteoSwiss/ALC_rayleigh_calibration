# -*- coding: utf-8 -*-
"""Does the temperature-corrected CHM15k constant agree better with CL61 and EARLINET? (report 19, part 9)

The paper-v0.9 validation (validation/paper/run_paper_v09.py) is re-run UNCHANGED, writing its caches to
this study's folder (the paper outputs are not touched), with three extra CL61 sites (Lindenberg,
Camborne, Aosta; CL61 cloud constant = reference, as at Payerne). The temperature correction is a pure
per-hour multiplicative factor on the CHM15k channel, so it is applied afterwards to the cached hourly
profiles:

    beta_corr(t) = beta(t) * exp(-kappa_u/100 * (T_int(t) - T_ref(t)))     (beta = backscatter)

kappa_u  = the unit's intra-night coefficient (q2a, elapsed time held fixed; the deployable estimate);
T_int(t) = hourly median of the 10-min HK cache;
T_ref(t) = +-15-day running median of T_int over the calibration windows of the nightly constants
           (the temperature the daily Kalman constant represents).
Statistics = run_paper_v09.pair_stats over the paper bands (stations 500-3000 m, EARLINET 500-5000 m),
all hours, night (22-02 UTC) and day (10-14 UTC), as in the paper figures.

Run:  python validation/chm15k_stability/q9_validation_tcorr.py            (process + score)
      python validation/chm15k_stability/q9_validation_tcorr.py --score    (score cached runs only)
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import BLOCKS_PER_DAY, FIGS, OUT, load_cache, load_constants  # noqa: E402
from q1_temperature import night_hk  # noqa: E402
import validation.paper.run_paper_v09 as V  # noqa: E402

warnings.filterwarnings("ignore")
VOUT = OUT / "validation_v09"
L1_END = "20260712"            # local L1 archive end for the non-Payerne sites
EXTRA_STATIONS = {
    "lindenberg": dict(name="Lindenberg", wmo="0-20000-0-10393", start="20250101", end="20260630", ref=0, channels=[
        V._ch("0-20000-0-10393", "C", "CL61", "cloud", "CL61 (cloud)"),
        V._ch("0-20000-0-10393", "0", "CHM15k", "rayleigh", "CHM15k (Rayleigh)")]),
    "camborne": dict(name="Camborne", wmo="0-20000-0-03808", start="20250910", end=L1_END, ref=0, channels=[
        V._ch("0-20000-0-03808", "C", "CL61", "cloud", "CL61 (cloud)"),
        V._ch("0-20000-0-03808", "A", "CHM15k", "rayleigh", "CHM15k (Rayleigh)")]),
    "aosta": dict(name="Aosta", wmo="0-380-5-1", start="20251031", end=L1_END, ref=0, channels=[
        V._ch("0-380-5-1", "B", "CL61", "cloud", "CL61 (cloud)"),
        V._ch("0-380-5-1", "0", "CHM15k", "rayleigh", "CHM15k (Rayleigh)")]),
}
PAYERNE = dict(V.STATIONS["payerne"])
PAYERNE["channels"] = PAYERNE["channels"][:2]            # CL61 cloud (ref) + CHM15k; CL31 not needed here
STATIONS = {"payerne": PAYERNE, **EXTRA_STATIONS}
EARLINET_CHM = {"leipzig": ("0-20000-0-10471", "0"), "magurele": ("0-20008-0-INO", "B"),
                "magurele_a": ("0-20008-0-INO", "A")}


def kappas():
    return {r["stream"]: float(r["b_T_ctrl"]) for r in csv.DictReader(open(OUT / "q2a_intranight.csv", encoding="utf-8"))}


def temperature_factor(wmo, ident, times, kappa_pct):
    """exp(-kappa (T_int(t) - T_ref(t))) at the given datetime64 times (hourly profiles)."""
    cache, cons = load_cache(wmo, ident), load_constants(wmo, ident)
    names = list(cache["hk_vars"])
    tb = (cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY).ravel()
    T = cache["hk"][..., names.index("temp_int")].ravel()
    order = np.argsort(tb)
    tb, T = tb[order], T[order]
    td = (np.asarray(times).astype("datetime64[s]").astype("int64") / 86400.0)
    # hourly median T_int around each profile time (+-30 min)
    Tt = np.full(td.size, np.nan)
    lo = np.searchsorted(tb, td - 1 / 48)
    hi = np.searchsorted(tb, td + 1 / 48)
    for i in range(td.size):
        seg = T[lo[i]:hi[i]]
        seg = seg[np.isfinite(seg)]
        if seg.size:
            Tt[i] = np.median(seg)
    # reference: T_int over the calibration windows, +-15-day running median
    Tn = night_hk(cache, cons["t0"], cons["t1"])[:, names.index("temp_int")]
    tc = cons["t"]
    ok = np.isfinite(Tn)
    tc, Tn = tc[ok], Tn[ok]
    Tref = np.full(td.size, np.nan)
    for i in range(td.size):
        m = np.abs(tc - td[i]) <= 15
        if m.sum() >= 3:
            Tref[i] = np.median(Tn[m])
    f = np.exp(-kappa_pct / 100.0 * (Tt - Tref))
    return np.where(np.isfinite(f), f, 1.0), Tt - Tref


def process():
    VOUT.mkdir(parents=True, exist_ok=True)
    V.OUT = VOUT                                   # caches go to this study's folder, paper outputs untouched
    V.ALL_STATIONS.update(STATIONS)
    for key in STATIONS:
        if not (VOUT / f"cache_{key}.npz").exists():
            print(f"processing {key} ...", flush=True)
            V.process_station(key)
    for key in EARLINET_CHM:
        if not (VOUT / f"cache_{key}.npz").exists():
            print(f"processing EARLINET {key} ...", flush=True)
            V.process_earlinet(key)


def hours_utc(times):
    return (np.asarray(times).astype("datetime64[h]").astype("int64") % 24)


def stats_sets(x, ref, z, band, times):
    h = hours_utc(times)
    out = {}
    for name, sel in (("all", np.ones(h.size, bool)), ("night", (h >= 22) | (h < 2)), ("day", (h >= 10) & (h < 14))):
        s = V.pair_stats(x[sel], ref[sel], z, band, [])
        out[name] = dict(med=s["med"], iqr=s["q75"] - s["q25"], r_log=s["r_log"], n=s["n_profiles"])
    return out


def score():
    B = kappas()
    rows = []
    for key, st in STATIONS.items():
        c = np.load(VOUT / f"cache_{key}.npz")
        z, t, beta, labels = c["z"], c["time"], c["beta"], list(c["labels"])
        ich = [i for i, l in enumerate(labels) if "CHM15k" in l][0]
        ch = st["channels"][ich]
        f, dT = temperature_factor(ch["wmo"], ch["ident"], t, B[f"{ch['wmo']}_{ch['ident']}"])
        ref = beta[st["ref"]].astype("f8")
        x0 = beta[ich].astype("f8")
        x1 = x0 * f[:, None]
        for var, x in (("operational", x0), ("T-corrected", x1)):
            s = stats_sets(x, ref, z, V.ZBAND_STATION, t)
            rows.append(dict(comparison=f"{st['name']}: CHM15k vs CL61 cloud", variant=var,
                             kappa=B[f"{ch['wmo']}_{ch['ident']}"], **{f"{m}_{k}": v[m] for k, v in s.items() for m in v}))
    for key, (wmo, ident) in EARLINET_CHM.items():
        c = np.load(VOUT / f"cache_{key}.npz")
        z, t, E, C = c["z"], c["time"], c["E"].astype("f8"), c["C"].astype("f8")
        f, dT = temperature_factor(wmo, ident, t, B[f"{wmo}_{ident}"])
        for var, x in (("operational", C), ("T-corrected", C * f[:, None])):
            s = stats_sets(x, E, z, V.ZBAND_EARLINET, t)
            rows.append(dict(comparison=f"{V.EARLINET[key]['name']} {ident}: CHM15k vs EARLINET {V.EARLINET[key]['lidar']}",
                             variant=var, kappa=B[f"{wmo}_{ident}"],
                             **{f"{m}_{k}": v[m] for k, v in s.items() for m in v}))
    for r in rows:
        r["day_minus_night"] = r["med_day"] - r["med_night"]
        print(f"  {r['comparison']:52s} {r['variant']:12s} kappa {r['kappa']:+.2f}  bias all {r['med_all']:+6.1f} % "
              f"(IQR {r['iqr_all']:5.1f}, r_log {r['r_log_all']:.3f}, n {r['n_all']})  night {r['med_night']:+6.1f}  "
              f"day {r['med_day']:+6.1f}  day-night {r['day_minus_night']:+6.1f}")
    with open(OUT / "q9_validation_tcorr.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    figure(rows)


def figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    comps = list(dict.fromkeys(r["comparison"] for r in rows))
    fig, axs = plt.subplots(1, 3, figsize=(21, 6))
    x = np.arange(len(comps))
    for ax, key, title in ((axs[0], "day_minus_night", "(a) day - night bias [%-pt] (0 = no diurnal bias)"),
                           (axs[1], "med_all", "(b) overall median bias vs reference [%]"),
                           (axs[2], "iqr_all", "(c) spread: IQR of the relative difference [%-pt]")):
        for var, col, dx in (("operational", "#7f7f7f", -0.18), ("T-corrected", "#d62728", 0.18)):
            v = [next(r[key] for r in rows if r["comparison"] == c and r["variant"] == var) for c in comps]
            ax.bar(x + dx, v, width=0.36, color=col, label=var)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x)
        short = []
        for c in comps:
            site, rest = c.split(":", 1)
            site = site.replace("Bucharest-Măgurele", "Măgurele").replace("Leipzig 0", "Leipzig")
            ref_name = "CL61" if "CL61" in rest else rest.strip().split()[-1]   # CL61 or the EARLINET lidar
            short.append(f"{site}\nvs {ref_name}")
        ax.set_xticklabels(short, fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=8)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "q9_validation_tcorr.png", dpi=130)
    print(f"-> {FIGS / 'q9_validation_tcorr.png'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--score", action="store_true")
    a = ap.parse_args()
    if not a.score:
        process()
    score()


if __name__ == "__main__":
    main()
