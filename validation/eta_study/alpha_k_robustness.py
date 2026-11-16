"""Robustness of the tail-slope extinction alpha_k (step 2, part B).

(i)   Fit-window sensitivity: k_app on +20..+80, +20..+150 and +40..+150 m above the peak (each
      inverted with the model k_app computed on the SAME window), plus the median SNR at +80 and
      +150 m (noise = std of the peak-normalised profile over +200..+300 m, where the cloud signal
      has died out for alpha >~ 10 /km).
(ii)  Causal-receiver test: per type, mean peak-normalised profiles in alpha_k quintiles ->
      10-90 % rise distance below the peak and tail k (+20..+150 m). A rise that stays put while
      the tail changes = atmosphere; both moving together = response/jitter. Then a fit of
      Gaussian pulse (nominal) x one-sided exponential receiver (tau_r free) x PVC cloud
      (a_G 3.3 um) to the type-mean and per-quintile mean profiles: does tau_r separate from alpha?
(iii) alpha_k by month per unit, and daily alpha_k vs Cloudnet LWP (rebuilt from der_scaled and
      N_scaled over the liquid column, lognormal sigma 0.35) and N_scaled (lowest 100 m, status 1)
      at Payerne / Palaiseau / Lindenberg, +-15 min around each profile.

Inputs: step2/profiles/*_working.parquet (fit_alpha_shape.py), incloud/<key>.npz, forward_lut_ext.npz.
Run:    python validation/eta_study/alpha_k_robustness.py [--workers 8] [--sample 60000]
Output: C:/DATA/Projects/202611_eta_study/step2b/
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_alpha_shape as F  # noqa: E402

OUT = F.BASE / "step2b"
WINDOWS = ((20.0, 80.0), (20.0, 150.0), (40.0, 150.0))
NOISE_WIN = (200.0, 300.0)
TAU_R_NS = np.array([0.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0])
C_LIGHT = 299_792_458.0


def kapp_window(P, zg, lo, hi):
    """Tail slope (1/km) of the peak-normalised profiles on gates lo..hi m above the peak."""
    k0 = int(np.argmin(np.abs(zg)))
    y = P / P[..., k0:k0 + 1]
    sel = (zg >= lo - 1e-6) & (zg <= hi + 1e-6)
    ys, zs = y[..., sel], zg[sel]
    ok = ys > F.KAPP_MIN_FRAC
    lny = np.where(ok, np.log(np.where(ok, ys, 1.0)), 0.0)
    n = ok.sum(-1)
    sz, szz = (ok * zs).sum(-1), (ok * zs ** 2).sum(-1)
    sy, szy = lny.sum(-1), (lny * zs).sum(-1)
    with np.errstate(all="ignore"):
        slope = (n * szy - sz * sy) / (n * szz - sz ** 2)
    return np.where(n >= 3, -slope * 1e3, np.nan)


def invert(kobs, cbh, kmod):
    """alpha from k via the model table kmod (cbhN, alphaF)."""
    node = np.clip(np.round((cbh - F.CBH_NODES[0]) / 50.0).astype(int), 0, F.CBH_NODES.size - 1)
    out = np.full(kobs.size, np.nan)
    for c in np.unique(node):
        m = node == c
        kk = kmod[c]
        ok = np.isfinite(kk)
        o = np.argsort(kk[ok])
        out[m] = np.interp(kobs[m], kk[ok][o], F.ALPHA_FINE[ok][o], left=np.nan, right=np.nan)
    return out


def load_unit(key, n_max, rng):
    """(P unit-sum windows, cbh, alpha_k, date, itype, dg, distorted) for up to n_max valid profiles."""
    import pandas as pd
    df = pd.read_parquet(F.OUT / "profiles" / f"{key}_working.parquet",
                         columns=["itype", "distorted", "valid", "cbh_m", "alpha_k", "date", "time_ns", "k_app"])
    z = np.load(F.SRC / f"{key}.npz", allow_pickle=True)
    dg = float(z["dr_m"])
    md = F.models(df.itype.iloc[0], "working", dg)
    kd = np.round(z["range_rel_m"] / dg).astype(int)
    cols = np.array([int(np.nonzero(kd == k)[0][0]) for k in md["ks"]])
    iv = np.nonzero(df.valid.to_numpy())[0]
    if iv.size > n_max:
        iv = np.sort(rng.choice(iv, n_max, replace=False))
    W = z["beta_before"][iv][:, cols].astype("float64")
    P = W / W.sum(1, keepdims=True)
    cb = df.cbh_m.to_numpy()[iv]
    return dict(P=P, cbh=np.where(np.isfinite(cb), cb, 1000.0), alpha_k=df.alpha_k.to_numpy()[iv],
                date=df.date.to_numpy()[iv], time_ns=df.time_ns.to_numpy()[iv], itype=df.itype.iloc[0],
                dg=dg, distorted=bool(df.distorted.iloc[0]), zg=md["zg"], md=md)


# ----------------------------------------------------------------------------- (ii) receiver model
def receiver_models(itype, dg, cbh, tau_ns):
    """Unit-sum peak-aligned model shapes (alphaF, K) for Gaussian pulse (nominal) x causal
    exponential receiver (tau) x PVC (a_G 3.3), at the CBH node closest to `cbh`."""
    L = F.lut()
    z = L["z_rel_m"]
    ti = list(L["types"]).index(itype)
    ag = L["aG_um"]
    j = np.searchsorted(ag, 3.3)
    w = (3.3 - ag[j - 1]) / (ag[j] - ag[j - 1])
    ic = int(np.argmin(np.abs(L["cbh_m"] - cbh)))
    b = L["beta_1m"][ti][:, j - 1, ic] * (1 - w) + L["beta_1m"][ti][:, j, ic] * w      # (alpha16, z)
    # extend the stored 1 m profiles (base-100..base+380 m) to +700 m with their own exponential
    # tail (log-slope of the last 50 m): the causal receiver shifts the peak upward, so the
    # peak+300 m window needs gates beyond +380 m (negligible signal except at the lowest alpha)
    z_ext = np.arange(z[-1] + 1.0, 700.0 + 0.5, 1.0)
    with np.errstate(all="ignore"):
        slope = (np.log(np.maximum(b[:, -1], 1e-300)) - np.log(np.maximum(b[:, -51], 1e-300))) / 50.0
    slope = np.minimum(np.nan_to_num(slope, nan=-1.0), 0.0)
    b = np.concatenate([b, b[:, -1:] * np.exp(slope[:, None] * (z_ext[None, :] - z[-1]))], axis=1)
    z = np.concatenate([z, z_ext])
    sig = C_LIGHT * F.NOMINAL[itype] * 1e-9 / 2 / 2.354820045
    kk = np.arange(-int(5 * sig), int(5 * sig) + 1)
    g = np.exp(-0.5 * (kk / sig) ** 2)
    g /= g.sum()
    b = np.apply_along_axis(lambda x: np.convolve(x, g, mode="same"), -1, b)
    if tau_ns > 0:
        Lr = C_LIGHT * tau_ns * 1e-9 / 2.0                                       # m
        kr = np.arange(0, int(8 * Lr) + 1)
        e = np.exp(-kr / Lr)
        e /= e.sum()
        # causal: the receiver delays the signal -> it appears at LARGER range (smear upward)
        b = np.apply_along_axis(lambda x: np.convolve(x, e, mode="full")[:x.size], -1, b)
    ks = F.gate_offsets(dg)
    S = F.peak_aligned(F.gate_sample_conv(b, z, dg), ks)                         # (alpha16, K)
    S = F._interp_axis(S, np.log(L["alpha_km"]), np.log(F.ALPHA_FINE), axis=0)
    return S / S.sum(-1, keepdims=True)


def rise_distance(mean_prof, zg):
    """10-90 % rise distance below the peak of a peak-normalised mean profile (m)."""
    k0 = int(np.argmin(np.abs(zg)))
    y = mean_prof / mean_prof[k0]
    zb, yb = zg[:k0 + 1], y[:k0 + 1]                                    # below the peak, ascending z

    def cross(level):
        idx = np.where(yb >= level)[0]
        if idx.size == 0 or idx[0] == 0:
            return np.nan
        i = idx[0]
        return zb[i - 1] + (zb[i] - zb[i - 1]) * (level - yb[i - 1]) / (yb[i] - yb[i - 1])
    return cross(0.9) - cross(0.1)


def receiver_fit(itype, dg, cbh, Pmean):
    """Grid fit over (alpha, tau_r): best alpha and chi2 per tau_r."""
    rows = []
    for tau in TAU_R_NS:
        S = receiver_models(itype, dg, cbh, tau)
        sse = ((S - Pmean[None, :]) ** 2).sum(-1)
        b = int(np.argmin(sse))
        rows.append(dict(tau_r_ns=tau, alpha_best=float(F.ALPHA_FINE[b]), chi2=float(sse[b]),
                         alpha_edge=bool(b in (0, F.ALPHA_FINE.size - 1))))
    return rows


# ----------------------------------------------------------------------------- (iii) Cloudnet LWP / N
def cloudnet_lwp_day(job):
    """Per ceilometer profile: median over +-15 min of the Cloudnet liquid-column LWP (g m-2, rebuilt)
    and of N_scaled in the lowest 100 m of the liquid layer next to the CBH (status 1)."""
    from netCDF4 import Dataset
    cn_site, ds, rows = job
    fp = F.CLOUDNET / f"{ds}_{cn_site}_der.nc"
    if not fp.exists():
        return []
    day0 = np.datetime64(f"{ds[:4]}-{ds[4:6]}-{ds[6:]}T00:00:00").astype("datetime64[s]").astype(np.int64)
    ts = np.array([r[1] for r in rows], float)
    try:
        with Dataset(fp) as nc:
            cnt = np.asarray(nc.variables["time"][:], float) * 3600.0 + day0
            h = np.asarray(nc.variables["height"][:], float)
            salt = float(np.asarray(nc.variables["altitude"][:]).ravel()[0])
            i0 = int(max(np.searchsorted(cnt, ts.min() - F.CN_WINDOW_S) - 1, 0))
            i1 = int(min(np.searchsorted(cnt, ts.max() + F.CN_WINDOW_S) + 1, cnt.size))
            sl = (slice(i0, i1), slice(None))
            st = np.ma.filled(np.ma.asarray(nc.variables["der_retrieval_status"][sl], dtype=int), 0)
            rs = np.ma.filled(np.ma.asarray(nc.variables["der_scaled"][sl], dtype=float), np.nan)
            nn = np.ma.filled(np.ma.asarray(nc.variables["N_scaled"][sl], dtype=float), np.nan)
    except Exception:  # noqa: BLE001
        return []
    cnt = cnt[i0:i1]
    dz = float(np.median(np.diff(h)))
    liq = (st >= 1) & (st <= 3)
    rs = np.where((rs > 0) & (rs < 1e-3), rs, np.nan)
    nn = np.where(nn > 0, nn, np.nan)
    rv = rs * np.exp(-F.SIGMA_LN ** 2)
    lwc = 4.0 / 3.0 * np.pi * 1.0e6 * nn * rv ** 3                         # g m-3
    lwp = np.nansum(np.where(liq, lwc, 0.0), axis=1) * dz                  # g m-2 (liquid column)
    lwp_ok = np.isfinite(np.where(liq, lwc, np.nan)).any(axis=1)
    out = []
    for (row, t, cbh) in rows:
        lo, hi = np.searchsorted(cnt, t - F.CN_WINDOW_S), np.searchsorted(cnt, t + F.CN_WINDOW_S, side="right")
        cloud_msl = cbh + salt
        L_, N_ = [], []
        for r in range(lo, hi):
            if lwp_ok[r]:
                L_.append(lwp[r])
            cand = np.where(liq[r] & (np.abs(h - cloud_msl) <= F.BASE_SEARCH_M))[0]
            if cand.size:
                jb = int(cand.min())
                while jb > 0 and liq[r, jb - 1]:
                    jb -= 1
                g = np.arange(jb, h.size)
                g = g[((h[g] - (h[jb] - dz / 2)) < 100.0)]
                g = g[liq[r, g] & (st[r, g] == 1)]
                v = nn[r, g]
                v = v[np.isfinite(v)]
                if v.size:
                    N_.append(np.median(v))
        if L_ or N_:
            out.append((row, float(np.median(L_)) if L_ else np.nan,
                        float(np.median(N_)) * 1e-6 if N_ else np.nan))
    return out


def main():
    import pandas as pd
    from scipy.stats import spearmanr
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--sample", type=int, default=60000, help="max profiles per unit for (i)/(ii)")
    ap.add_argument("--skip-i", action="store_true", help="reuse window_sensitivity_by_unit.csv")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    keys = sorted(p.stem.replace("_working", "") for p in (F.OUT / "profiles").glob("*_working.parquet"))
    data = {k: load_unit(k, a.sample, rng) for k in keys}
    log = open(OUT / "summary.txt", "w", encoding="utf-8")

    def say(s=""):
        print(s, flush=True)
        log.write(s + "\n")

    # ---------------- (i) window sensitivity + SNR
    say("== (i) alpha_k vs fit window; SNR at +80 / +150 m (noise = std over +200..+300 m) ==")
    rows = []
    for k, d in ({} if a.skip_i else data).items():
        md, zg = d["md"], d["zg"]
        pi = list(F.PULSES).index(F.NOMINAL[d["itype"]])
        k0 = int(np.argmin(np.abs(zg)))
        y = d["P"] / d["P"][:, k0:k0 + 1]
        nz = (zg >= NOISE_WIN[0]) & (zg <= NOISE_WIN[1])
        sig = np.std(y[:, nz], axis=1)
        r = dict(key=k, itype=d["itype"], dg=d["dg"], distorted=d["distorted"], n=len(d["cbh"]),
                 snr80=float(np.nanmedian(np.interp(80.0, zg, np.nanmedian(y, 0)) / sig)),
                 snr150=float(np.nanmedian(np.interp(150.0, zg, np.nanmedian(y, 0)) / sig)))
        snr_p80 = y[:, np.argmin(np.abs(zg - 80))] / sig
        snr_p150 = y[:, np.argmin(np.abs(zg - 150))] / sig
        r["snr80_profile_med"] = float(np.nanmedian(snr_p80))
        r["snr150_profile_med"] = float(np.nanmedian(snr_p150))
        r["frac_snr150_below3"] = float(np.mean(snr_p150 < 3))
        for lo, hi in WINDOWS:
            kobs = kapp_window(d["P"], zg, lo, hi)
            kmod = kapp_window(md["Mh"][:, pi, 1], zg, lo, hi)                     # (cbhN, alphaF)
            al = invert(kobs, d["cbh"], kmod)
            r[f"alpha_k_{int(lo)}_{int(hi)}"] = float(np.nanmedian(al))
            r[f"frac_out_{int(lo)}_{int(hi)}"] = float(np.mean(~np.isfinite(al)))
        rows.append(r)
    if a.skip_i:
        wdf = pd.read_csv(OUT / "window_sensitivity_by_unit.csv")
        say("(reused window_sensitivity_by_unit.csv)")
    else:
        wdf = pd.DataFrame(rows)
        wdf.to_csv(OUT / "window_sensitivity_by_unit.csv", index=False, float_format="%.4g")
    agg = {c: "median" for c in wdf.columns if c.startswith(("alpha_k_", "snr", "frac_"))}
    t = wdf.groupby(["itype", "dg"]).agg({**agg, "n": "sum", "key": "count"}).rename(columns={"key": "n_units"})
    say(t.round(3).to_string())
    say("(medians of the per-unit medians; SNR medians are per-profile)")

    # ---------------- (ii) quintiles: rise distance and tail k; receiver fit
    say("\n== (ii) alpha_k quintiles: mean peak-normalised profile -> rise 10-90 % (m) and tail k (1/km) ==")
    qrows, frows = [], []
    for it in ("CL31", "CL51", "CL61"):
        ds = [d for d in data.values() if d["itype"] == it and not d["distorted"]
              and abs(d["dg"] - (9.6 if it == "CL61" else 10.0)) < 0.05]
        if not ds:
            continue
        P = np.vstack([d["P"] for d in ds])
        ak = np.concatenate([d["alpha_k"] for d in ds])
        cb = np.concatenate([d["cbh"] for d in ds])
        zg = ds[0]["zg"]
        dg = ds[0]["dg"]
        k0 = int(np.argmin(np.abs(zg)))
        ok = np.isfinite(ak)
        edges = np.nanpercentile(ak[ok], [0, 20, 40, 60, 80, 100])
        cbh_type = float(np.nanmedian(cb))
        # whole-type mean
        Pm = np.nanmean(P[ok], 0)
        Pm /= Pm.sum()
        for r in receiver_fit(it, dg, cbh_type, Pm):
            frows.append(dict(itype=it, quintile="all", **r))
        for qi in range(5):
            m = ok & (ak >= edges[qi]) & (ak <= edges[qi + 1])
            Pq = np.nanmean(P[m], 0)
            yq = Pq / Pq[k0]
            ktail = float(kapp_window(Pq[None, :], zg, 20.0, 150.0)[0])
            qrows.append(dict(itype=it, quintile=qi + 1, alpha_k_lo=edges[qi], alpha_k_hi=edges[qi + 1],
                              alpha_k_med=float(np.nanmedian(ak[m])), n=int(m.sum()),
                              cbh_med=float(np.nanmedian(cb[m])), rise_10_90_m=rise_distance(yq, zg),
                              tail_k=ktail, y_at_minus_gate=float(yq[k0 - 1]),
                              y_at_minus_2gate=float(yq[k0 - 2])))
            Pqn = Pq / Pq.sum()
            for r in receiver_fit(it, dg, float(np.nanmedian(cb[m])), Pqn):
                frows.append(dict(itype=it, quintile=qi + 1, **r))
    qdf = pd.DataFrame(qrows)
    qdf.to_csv(OUT / "quintile_rise_tail.csv", index=False, float_format="%.4g")
    say(qdf.round(3).to_string(index=False))
    fdf = pd.DataFrame(frows)
    fdf.to_csv(OUT / "receiver_fit.csv", index=False, float_format="%.4g")
    say("\n== (ii) Gaussian pulse (nominal) x causal exponential receiver (tau_r) x PVC: best alpha and chi2 per tau_r ==")
    say(fdf.pivot_table(index=["itype", "quintile"], columns="tau_r_ns", values="alpha_best").round(1).to_string())
    say(fdf.pivot_table(index=["itype", "quintile"], columns="tau_r_ns", values="chi2").applymap(lambda v: f"{v:.2e}").to_string())
    best = fdf.loc[fdf.groupby(["itype", "quintile"]).chi2.idxmin()]
    say("best (tau_r, alpha) per type/quintile:\n" + best[["itype", "quintile", "tau_r_ns", "alpha_best", "chi2"]].to_string(index=False))

    # ---------------- (iii) monthly alpha_k per unit; daily alpha_k vs Cloudnet LWP and N
    say("\n== (iii) alpha_k by month (per unit, all valid profiles) ==")
    mon = []
    for k in keys:
        df = pd.read_parquet(F.OUT / "profiles" / f"{k}_working.parquet",
                             columns=["itype", "site", "valid", "alpha_k", "date", "cbh_m"])
        df = df[df.valid]
        df["month"] = df.date // 100
        g = df.groupby("month").alpha_k.agg(["median", "size"]).reset_index()
        g.insert(0, "key", k)
        g.insert(1, "itype", df.itype.iloc[0])
        g.insert(2, "site", df.site.iloc[0])
        mon.append(g)
    mdf = pd.concat(mon, ignore_index=True)
    mdf.to_csv(OUT / "alpha_k_by_month.csv", index=False, float_format="%.4g")
    piv = mdf[mdf["size"] >= 200].pivot_table(index="month", columns="itype", values="median", aggfunc="median")
    say("type median of the unit monthly medians (months with >= 200 profiles):\n" + piv.round(1).to_string())
    # seasonal amplitude per unit
    amp = mdf[mdf["size"] >= 200].groupby("key")["median"].agg(lambda x: (x.max() - x.min()) / x.median() if len(x) >= 6 else np.nan)
    say(f"per-unit monthly range / median: median {np.nanmedian(amp):.2f}, P10 {np.nanpercentile(amp.dropna(), 10):.2f}, "
        f"P90 {np.nanpercentile(amp.dropna(), 90):.2f}")

    say("\n== (iii) daily alpha_k vs Cloudnet LWP (rebuilt) and N_scaled (lowest 100 m), +-15 min ==")
    jobs, meta = defaultdict(list), {}
    for k in F.CN_SITES:
        p = F.OUT / "profiles" / f"{k}_working.parquet"
        if not p.exists():
            continue
        df = pd.read_parquet(p, columns=["valid", "alpha_k", "date", "time_ns", "cbh_m"])
        df = df[df.valid & np.isfinite(df.cbh_m)]
        meta[k] = df
        for row, t, cb in zip(df.index, df.time_ns, df.cbh_m):
            tt = int(t // 1_000_000_000)
            ds = str(np.datetime64(tt, "s").astype("datetime64[D]")).replace("-", "")
            jobs[(k, F.CN_SITES[k], ds)].append((row, tt, float(cb)))
    keys_j = sorted(jobs)
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        res = list(ex.map(cloudnet_lwp_day, [(s, d, jobs[(k, s, d)]) for (k, s, d) in keys_j], chunksize=2))
    crow = []
    for (k, s, d), r in zip(keys_j, res):
        df = meta[k]
        for row, lwp, nsc in r:
            crow.append((k, row, lwp, nsc))
    cdf = pd.DataFrame(crow, columns=["key", "row", "lwp_gm2", "N_cm3"])
    rows = []
    for k, df in meta.items():
        c = cdf[cdf.key == k].set_index("row")
        j = df.join(c[["lwp_gm2", "N_cm3"]], how="inner")
        dd = j.groupby("date")[["alpha_k", "lwp_gm2", "N_cm3"]].median().dropna(subset=["alpha_k"])
        r = dict(key=k, n_profiles=len(j), n_days=len(dd),
                 rs_alpha_lwp_daily=spearmanr(dd.alpha_k, dd.lwp_gm2, nan_policy="omit").statistic,
                 rs_alpha_N_daily=spearmanr(dd.alpha_k, dd.N_cm3, nan_policy="omit").statistic,
                 rs_alpha_lwp_profile=spearmanr(j.alpha_k, j.lwp_gm2, nan_policy="omit").statistic,
                 rs_alpha_N_profile=spearmanr(j.alpha_k, j.N_cm3, nan_policy="omit").statistic,
                 lwp_med=float(np.nanmedian(j.lwp_gm2)), N_med=float(np.nanmedian(j.N_cm3)),
                 alpha_k_med=float(np.nanmedian(j.alpha_k)))
        rows.append(r)
        dd.to_csv(OUT / f"daily_alpha_k_lwp_N_{k}.csv", float_format="%.4g")
    ldf = pd.DataFrame(rows)
    ldf.to_csv(OUT / "alpha_k_vs_lwp_N.csv", index=False, float_format="%.4g")
    say(ldf.round(3).to_string(index=False))
    figure(qdf, fdf, wdf, mdf, ldf, say)
    log.close()


def figure(qdf, fdf, wdf, mdf, ldf, say):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = {"CL31": "C0", "CL51": "C2", "CL61": "C3"}
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.4))
    ax = axes[0]
    for it, g in qdf.groupby("itype"):
        ax.plot(g.alpha_k_med, g.rise_10_90_m, "o-", color=cols[it], label=f"{it} rise 10-90 %")
        ax2 = ax.twinx() if it == "CL31" else ax2
        ax2.plot(g.alpha_k_med, g.tail_k, "s--", color=cols[it], alpha=0.7, label=f"{it} tail k")
    ax.set_xlabel("alpha_k quintile median [km$^{-1}$]")
    ax.set_ylabel("rise distance 10-90 % below the peak [m] (solid)")
    ax2.set_ylabel("tail k on +20..+150 m [km$^{-1}$] (dashed)")
    ax.legend(loc="upper left", fontsize=8)
    ax2.legend(loc="lower right", fontsize=8)
    ax.set_title("(ii) rise vs tail per alpha_k quintile")
    ax.grid(alpha=0.3)
    ax = axes[1]
    for (it, qq), g in fdf[fdf.quintile.astype(str) == "all"].groupby(["itype", "quintile"]):
        ax.plot(g.tau_r_ns, g.chi2 / g.chi2.min(), "o-", color=cols[it], label=f"{it}: chi2 / min")
        for _, r in g.iterrows():
            ax.annotate(f"{r.alpha_best:.0f}", (r.tau_r_ns, r.chi2 / g.chi2.min()), fontsize=7, color=cols[it])
    ax.set_yscale("log")
    ax.set_xlabel("receiver time constant tau_r [ns]")
    ax.set_ylabel("chi2 / min chi2 (labels: best alpha)")
    ax.set_title("(ii) pulse x receiver x PVC fit to the type-mean profile")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax = axes[2]
    piv = mdf[mdf["size"] >= 200].pivot_table(index="month", columns="itype", values="median", aggfunc="median")
    for it in piv.columns:
        ax.plot([str(m) for m in piv.index], piv[it], "o-", color=cols.get(it, "k"), label=it)
    ax.set_xticks(range(0, len(piv.index), 2))
    ax.set_xticklabels([str(m) for m in piv.index][::2], rotation=45, fontsize=8)
    ax.set_ylabel("alpha_k [km$^{-1}$] (median of unit monthly medians)")
    ax.set_title("(iii) alpha_k by month")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.suptitle("alpha_k robustness: rise vs tail, receiver-response fit, seasonality", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "fig_alpha_k_robustness.png", dpi=130)
    plt.close(fig)
    say("figure written")


if __name__ == "__main__":
    main()
