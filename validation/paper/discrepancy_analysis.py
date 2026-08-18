"""
discrepancy_analysis.py — quantitative probes behind the "where do the differences come from"
section of the paper-validation report (all on the uniform L1 pipeline of run_paper_validation):

  A1. Payerne CL61 Rayleigh vs cloud calibration series: how thin is the Rayleigh series, what
      constant ratio does the Kalman hold, and how it drifts.
  A2. Uccle CL51 (and CL31/CL61) cloud-calibration oscillations: period, amplitude, seasonal cycle.
  B1. Day/night + seasonal splits of the intercomparison statistics per station/channel
      (post-hoc masks on the hourly synchronized matrices; no reprocessing).
  B2. Mini-MPL 532->1064 Angstrom-exponent sensitivity (analytic on the median profile).
  B3. Altitude-band dependence of the bias (500-1000 / 1000-2000 / 2000-3000 m): separates
      near-range/overlap-type mechanisms (bias largest low) from conversion/noise mechanisms
      (bias growing where the aerosol fraction shrinks).
  B4. Monthly evolution of the median relative bias: a flat line is a scale-like (calibration)
      offset; a seasonal cycle points at the 910 nm WV residual / laser-wavelength drift; a drift
      or step points at an instrument change.

Each station is processed ONCE through run_paper_validation.run_site (module-level cache) and the
B-probes are post-hoc masks/bands on the synchronized science matrices — so every number here is
exactly the paper pipeline, only split.

Outputs figures to OUT (landscape) + a JSON with the numbers used by the report.
Usage: python -m validation.paper.discrepancy_analysis [csv|full]
  csv  = A1+A2 only (fast, no station reprocessing)
  full = everything (runs the 7 benchmark stations through the L1 pipeline once)
"""
from __future__ import annotations
import csv as _csv
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
CALIB = OUT / "calib"
RESULTS = OUT / "discrepancy_analysis.json"


def _read_calib(key):
    """CSV -> dict of arrays (time, C_daily, C_daily_std, C_kalman, C_kalman_std)."""
    f = CALIB / f"{key}.csv"
    if not f.is_file():
        return None
    rows = list(_csv.DictReader(open(f, encoding="utf-8")))
    def col(name):
        return np.array([float(r[name]) if r.get(name) not in ("", "nan", None) else np.nan for r in rows])
    t = np.array([np.datetime64(r["time"][:10]) for r in rows])
    return dict(time=t, C_daily=col("C_daily"), C_kalman=col("C_kalman"),
                C_kalman_std=col("C_kalman_std") if "C_kalman_std" in rows[0] else np.full(t.size, np.nan))


# ---------------------------------------------------------------------------
# A1. Payerne CL61: Rayleigh vs cloud lidar-constant series (both TRUE C_L from the calout)
# ---------------------------------------------------------------------------
RAY_COL = "#1f77b4"    # general guideline: blue = Rayleigh
CLD_COL = "#404040"    # dark grey = cloud


def payerne_cl61(res):
    """Both series read as the absolute Wiegner lidar constant C_L (the calibration CSVs now hold
    C_L for BOTH methods). The two methods measure the same physical constant; their ratio is the
    method discrepancy that the intercomparison sees between the two CL61 entries."""
    ray = _read_calib("0-20000-0-06610_C_rayleigh_L1")
    cld = _read_calib("0-20000-0-06610_C_cloud_L1")
    if ray is None or cld is None:
        print("payerne_cl61: series missing"); return
    # align on common dates
    _, ir, ic = np.intersect1d(ray["time"], cld["time"], return_indices=True)
    kr, kc = ray["C_kalman"][ir], cld["C_kalman"][ic]
    ratio = kr / kc
    nray = int(np.isfinite(ray["C_daily"]).sum()); ncld = int(np.isfinite(cld["C_daily"]).sum())
    res["payerne_cl61"] = dict(
        n_raw_rayleigh_nights=nray, n_raw_cloud_days=ncld,
        span_days=int(ray["time"].size),
        median_CL_rayleigh=float(np.nanmedian(kr)), median_CL_cloud=float(np.nanmedian(kc)),
        median_ratio=float(np.nanmedian(ratio)),
        ratio_p10=float(np.nanpercentile(ratio, 10)), ratio_p90=float(np.nanpercentile(ratio, 90)),
        # beta_ray/beta_cloud = C_L(cloud)/C_L(rayleigh): the method discrepancy in beta space
        expected_relbias_pct=float(100 * (np.nanmedian(kc) / np.nanmedian(kr) - 1)),
    )
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    for d, c, nm in ((ray, RAY_COL, "Rayleigh"), (cld, CLD_COL, "cloud")):
        ax[0].plot(d["time"], d["C_daily"], "x", color=c, ms=4, alpha=0.6)
        ax[0].plot(d["time"], d["C_kalman"], "-", color=c, lw=1.8,
                   label=f"{nm} (n_raw={np.isfinite(d['C_daily']).sum()})")
    ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[0].set_ylabel("lidar constant C$_L$")
    ax[0].set_title("(a) Payerne CL61: the SAME physical C$_L$ from the two methods")
    ax[1].plot(ray["time"][ir], ratio, "-", color="C3", lw=1.5)
    ax[1].axhline(np.nanmedian(ratio), color="k", ls="--", lw=1,
                  label=f"median {np.nanmedian(ratio):.2f}")
    ax[1].set_ylabel("C$_L$(Rayleigh) / C$_L$(cloud)"); ax[1].grid(alpha=0.3); ax[1].legend()
    ax[1].set_title("(b) Method ratio -> the CL61 Rayleigh-vs-cloud offset in the validation")
    for a in ax:
        a.xaxis.set_major_formatter(mdates.DateFormatter("%b %y"))
    fig.suptitle("Payerne CL61: Rayleigh vs cloud lidar constant (blue = Rayleigh, dark grey = cloud)",
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "fig_report_cl61ray_payerne.png", dpi=160); plt.close(fig)
    print("payerne_cl61:", res["payerne_cl61"])


# ---------------------------------------------------------------------------
# A2. CL51/CL31 cloud-calibration oscillations (Uccle CL51 ref; Payerne CL31 for comparison)
# ---------------------------------------------------------------------------
def _spectrum(t, y):
    """Lomb-Scargle-ish via FFT on the daily-interpolated series; returns period [d], power."""
    m = np.isfinite(y)
    if m.sum() < 60:
        return None
    td = (t - t[0]).astype("timedelta64[D]").astype(float)
    ti = np.arange(td[m].min(), td[m].max() + 1)
    yi = np.interp(ti, td[m], y[m])
    yi = yi - np.nanmean(yi)
    win = np.hanning(yi.size)
    P = np.abs(np.fft.rfft(yi * win)) ** 2
    fr = np.fft.rfftfreq(yi.size, d=1.0)
    per = np.full(fr.size, np.inf); per[1:] = 1.0 / fr[1:]
    return per[1:], P[1:]


def cl51_oscillation(res):
    # per-type colours: CL51 purple, CL31 orange, CL61-cloud dark grey (general colour rules)
    series = [("0-20000-0-06447_A_cloud_L1", "Uccle CL51 (cloud)", "#9467bd"),
              ("0-20000-0-06610_B_cloud_L1", "Payerne CL31 (cloud)", "#ff7f0e"),
              ("0-20000-0-06610_C_cloud_L1", "Payerne CL61 (cloud)", CLD_COL)]
    fig, axes = plt.subplots(len(series), 3, figsize=(19, 3.4 * len(series)))
    out = {}
    for i, (key, label, lcol) in enumerate(series):
        d = _read_calib(key)
        if d is None:
            continue
        t, cd, ck = d["time"], d["C_daily"], d["C_kalman"]
        med = np.nanmedian(ck)
        rel = 100 * (cd / med - 1)          # daily C_L departure from the series median [%]
        relk = 100 * (ck / med - 1)
        ax = axes[i]
        ax[0].plot(t, rel, ".", ms=2.5, color="0.6")
        ax[0].plot(t, relk, "-", color=lcol, lw=1.6)
        ax[0].set_ylabel("C$_L$ / median - 1 [%]"); ax[0].grid(alpha=0.3)
        ax[0].set_title(f"(a{i+1}) {label}: daily (dots) + Kalman")
        ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%b %y"))
        sp = _spectrum(t, cd)
        dom_period = np.nan
        if sp:
            per, P = sp
            sel = (per > 3) & (per < 400)
            ax[1].semilogx(per[sel], P[sel] / P[sel].max(), color="C1")
            dom_period = float(per[sel][np.argmax(P[sel])])
            ax[1].axvline(dom_period, color="k", ls="--", lw=1, label=f"dominant ~{dom_period:.0f} d")
            ax[1].legend(fontsize=8)
        ax[1].set_xlabel("period [days]"); ax[1].set_ylabel("norm. power"); ax[1].grid(alpha=0.3)
        ax[1].set_title(f"(b{i+1}) spectrum of the daily series")
        # seasonal cycle: monthly median departure
        months = t.astype("datetime64[M]").astype(int) % 12 + 1
        mm = [np.nanmedian(rel[(months == m) & np.isfinite(rel)]) if ((months == m) & np.isfinite(rel)).any()
              else np.nan for m in range(1, 13)]
        ax[2].bar(range(1, 13), mm, color="C2")
        ax[2].set_xticks(range(1, 13)); ax[2].set_xlabel("month"); ax[2].set_ylabel("median dep. [%]")
        ax[2].grid(alpha=0.3, axis="y"); ax[2].set_title(f"(c{i+1}) seasonal cycle")
        amp = float(np.nanmax(mm) - np.nanmin(mm)) if np.isfinite(mm).any() else np.nan
        out[key] = dict(label=label, n_raw=int(np.isfinite(cd).sum()),
                        daily_scatter_pct=float(np.nanstd(rel)),
                        kalman_p2p_pct=float(np.nanmax(relk) - np.nanmin(relk)),
                        dominant_period_days=dom_period, seasonal_amplitude_pct=amp)
        print(key, out[key])
    res["oscillations"] = out
    fig.suptitle("Cloud-calibration constant stability: time series, spectra, seasonal cycles",
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT / "fig_report_cl51_oscillation.png", dpi=160); plt.close(fig)


# ---------------------------------------------------------------------------
# Shared station cache: each station runs ONCE through the uniform L1 pipeline
# (run_paper_validation.run_site), all B-probes reuse the synchronized matrices.
# ---------------------------------------------------------------------------
STATIONS = ["payerne", "amsterdam", "uccle", "sirta", "lindenberg", "aosta", "camborne"]
_SITE_CACHE = {}


def _get_site(name):
    """(R, cfg) from the paper pipeline, computed once per process."""
    if name not in _SITE_CACHE:
        import warnings
        warnings.filterwarnings("ignore")
        from validation.paper import run_paper_validation as RPV
        out = RPV.run_site(name)
        _SITE_CACHE[name] = (out[0], out[1]) if out is not None else (None, None)
    return _SITE_CACHE[name]


def _channel_colors(R):
    from validation.paper import figures as FIG
    return FIG.channel_colors(R["channels"])


# ---------------------------------------------------------------------------
# B1. day/night + seasonal splits of the synchronized matrices
# ---------------------------------------------------------------------------
def daynight_seasonal(res):
    from validation.paper import intercompare as IC
    from calibration.sensitivity.noise import solar_elevation
    splits = {}
    fig, axes = plt.subplots(2, len(STATIONS), figsize=(5.0 * len(STATIONS), 8.5))
    for col, name in enumerate(STATIONS):
        R, cfg = _get_site(name)
        if R is None:
            continue
        z = np.asarray(R["altGrid"]) - R["station"]["altitude"]
        zmask = (z >= cfg["zMin"]) & (z <= cfg["zMax"])
        tt = np.asarray(R["time_sync"]).astype("datetime64[s]")
        elev = solar_elevation(tt, R["station"]["lat"], R["station"]["lon"])
        isday = elev > 5.0
        month = (tt.astype("datetime64[M]").astype(int) % 12) + 1
        seas = {"DJF": np.isin(month, (12, 1, 2)), "MAM": np.isin(month, (3, 4, 5)),
                "JJA": np.isin(month, (6, 7, 8)), "SON": np.isin(month, (9, 10, 11))}
        iref = cfg["referenceChannel"]; ref = R["beta"][iref]
        cols = _channel_colors(R)
        st = {}
        for k, ch in enumerate(R["channels"]):
            if k == iref:
                continue
            cur = R["beta"][k]
            def _sub(mask):
                return IC._stats(cur[mask], ref[mask], zmask)
            entry = dict(all=IC._stats(cur, ref, zmask), day=_sub(isday), night=_sub(~isday),
                         **{sn: _sub(sm) for sn, sm in seas.items() if sm.any()})
            st[ch["label"]] = entry
        splits[name] = st
        # panel: med relbias (top) + log r (bottom) per channel for day/night/season
        cats = ["all", "day", "night", "DJF", "MAM", "JJA", "SON"]
        labels = list(st.keys())
        ccols = [cols[k] for k, ch in enumerate(R["channels"]) if k != iref]
        w = 0.8 / max(len(labels), 1)
        for row, met, ylab in ((0, "medrelbias_pct", "med relbias [%]"), (1, "r_log", "log r")):
            ax = axes[row][col]
            for j, lb in enumerate(labels):
                v = [st[lb].get(cat, {}).get(met, np.nan) for cat in cats]
                ax.bar(np.arange(len(cats)) + j * w, v, w, color=ccols[j],
                       label=lb if row == 0 else None)
            ax.set_xticks(np.arange(len(cats)) + 0.4 - w / 2); ax.set_xticklabels(cats, fontsize=8)
            ax.grid(alpha=0.3, axis="y"); ax.set_ylabel(ylab)
            if row == 0:
                ax.set_title(name); ax.legend(fontsize=7)
    res["daynight_seasonal"] = _to_jsonable(splits)
    fig.suptitle("When do the results differ? day/night + seasonal splits of the intercomparison metrics",
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "fig_report_daynight_seasonal.png", dpi=160); plt.close(fig)
    print("daynight_seasonal: done", list(splits.keys()))


# ---------------------------------------------------------------------------
# B3. altitude-band dependence of the bias — WHERE in the profile does it sit?
# A near-range mechanism (overlap, incomplete offset removal) is largest in the lowest band and
# fades aloft; a conversion/noise mechanism (532->1064 molaer, 910->1064 in clean air, SNR floor)
# GROWS with altitude as the aerosol fraction shrinks; a pure calibration-scale offset is flat.
# ---------------------------------------------------------------------------
BANDS = ((500.0, 1000.0), (1000.0, 2000.0), (2000.0, 3000.0))


def altitude_bands(res):
    from validation.paper import intercompare as IC
    out = {}
    fig, axes = plt.subplots(1, len(STATIONS), figsize=(4.6 * len(STATIONS), 5.5))
    for col, name in enumerate(STATIONS):
        R, cfg = _get_site(name)
        if R is None:
            continue
        z = np.asarray(R["altGrid"]) - R["station"]["altitude"]
        iref = cfg["referenceChannel"]; ref = R["beta"][iref]
        cols = _channel_colors(R)
        st = {}
        ax = axes[col]
        xb = np.arange(len(BANDS))
        labels = [ch["label"] for k, ch in enumerate(R["channels"]) if k != iref]
        w = 0.8 / max(len(labels), 1)
        j = 0
        for k, ch in enumerate(R["channels"]):
            if k == iref:
                continue
            vals = []
            for z0, z1 in BANDS:
                zm = (z >= z0) & (z <= z1)
                s = IC._stats(R["beta"][k], ref, zm)
                vals.append(s["medrelbias_pct"])
            st[ch["label"]] = {f"{int(z0)}-{int(z1)}m": v for (z0, z1), v in zip(BANDS, vals)}
            ax.bar(xb + j * w, vals, w, color=cols[k], label=ch["label"])
            j += 1
        out[name] = st
        ax.set_xticks(xb + 0.4 - w / 2)
        ax.set_xticklabels([f"{int(z0/100)/10}-{int(z1/100)/10} km" for z0, z1 in BANDS], fontsize=8)
        ax.axhline(0, color="k", lw=0.8)
        ax.grid(alpha=0.3, axis="y"); ax.set_ylabel("med relbias [%]"); ax.set_title(name)
        ax.legend(fontsize=7)
    res["altitude_bands"] = _to_jsonable(out)
    fig.suptitle("Where in the profile does the bias sit? median relative bias per altitude band"
                 " (near-range mechanisms fade aloft; conversion/noise mechanisms grow)", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "fig_report_altitude_bands.png", dpi=160); plt.close(fig)
    print("altitude_bands: done", list(out.keys()))


# ---------------------------------------------------------------------------
# B4. monthly evolution of the median relative bias — WHEN does it move?
# Flat = scale-like (calibration constant); seasonal cycle = 910 nm WV residual / laser-wavelength
# drift (cf. A2 oscillations); drift/step = instrument or processing change.
# ---------------------------------------------------------------------------
def monthly_bias(res):
    from validation.paper import intercompare as IC
    out = {}
    fig, axes = plt.subplots(1, len(STATIONS), figsize=(4.6 * len(STATIONS), 5.0))
    for col, name in enumerate(STATIONS):
        R, cfg = _get_site(name)
        if R is None:
            continue
        z = np.asarray(R["altGrid"]) - R["station"]["altitude"]
        zmask = (z >= cfg["zMin"]) & (z <= cfg["zMax"])
        tt = np.asarray(R["time_sync"]).astype("datetime64[M]")
        months = np.unique(tt)
        iref = cfg["referenceChannel"]; ref = R["beta"][iref]
        cols = _channel_colors(R)
        ax = axes[col]
        st = {}
        for k, ch in enumerate(R["channels"]):
            if k == iref:
                continue
            mv = []
            for m in months:
                sel = tt == m
                s = IC._stats(R["beta"][k][sel], ref[sel], zmask)
                mv.append(s["medrelbias_pct"] if s["n"] > 200 else np.nan)   # skip thin months
            st[ch["label"]] = {str(m): v for m, v in zip(months, mv)}
            ax.plot(months.astype("datetime64[D]"), mv, "o-", color=cols[k], lw=1.5, ms=4,
                    label=ch["label"])
        out[name] = st
        ax.axhline(0, color="k", lw=0.8)
        ax.grid(alpha=0.3); ax.set_ylabel("med relbias [%]"); ax.set_title(name)
        ax.legend(fontsize=7)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b%y"))
        for lab in ax.get_xticklabels():
            lab.set_rotation(45); lab.set_fontsize(7)
    res["monthly_bias"] = _to_jsonable(out)
    fig.suptitle("When does the bias move? monthly median relative bias vs the site reference"
                 " (flat = scale-like; seasonal = WV/laser-drift residual; step = instrument change)",
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "fig_report_monthly_bias.png", dpi=160); plt.close(fig)
    print("monthly_bias: done", list(out.keys()))


# ---------------------------------------------------------------------------
# B2. Mini-MPL Angstrom-exponent sensitivity (analytic on the median profiles)
# ---------------------------------------------------------------------------
def minimpl_alpha(res):
    """The molaer conversion scales the AEROSOL part by (532/1064)^(-alpha) = 2^(-alpha) ... i.e.
    beta_aer(1064) = beta_aer(532) / 2^alpha. Changing alpha rescales beta_aer by 2^(1-alpha)
    relative to the alpha=1 reference, so the relbias moves by ~f_aer * (2^(1-alpha) - 1) where
    f_aer is the aerosol fraction of the (already converted) Mini-MPL signal in the stats band."""
    from validation.paper import intercompare as IC
    R, cfg = _get_site("sirta")
    if R is None:
        print("minimpl_alpha: no sirta"); return
    z = np.asarray(R["altGrid"]) - R["station"]["altitude"]
    zmask = (z >= cfg["zMin"]) & (z <= cfg["zMax"])
    k = next(i for i, c in enumerate(R["channels"]) if "MPL" in c["label"])
    iref = cfg["referenceChannel"]
    B, REF = R["beta"][k], R["beta"][iref]
    bmol1064 = IC._molecular_beta(z, R["station"]["altitude"], 1064.0)
    med = np.nanmedian(B[:, zmask], axis=0)
    baer = med - bmol1064[zmask]
    f_aer = float(np.nanmedian(np.clip(baer, 0, None) / med))
    base = IC._stats(B, REF, zmask)
    alphas = np.arange(0.3, 1.8001, 0.1)
    rel = []
    for a in alphas:
        scale = 2.0 ** (1.0 - a)             # beta_aer multiplier relative to alpha=1
        Bs = bmol1064[None, :] + (B - bmol1064[None, :]) * scale
        s = IC._stats(Bs, REF, zmask)
        rel.append(s["medrelbias_pct"])
    res["minimpl_alpha"] = dict(f_aer_med=f_aer, alphas=list(map(float, alphas)),
                                medrelbias_pct=list(map(float, rel)),
                                base_medrelbias_pct=float(base["medrelbias_pct"]),
                                alpha_for_zero=float(np.interp(0, rel[::-1], alphas[::-1]))
                                if (min(rel) < 0 < max(rel)) else np.nan)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(alphas, rel, "o-", color="C0")
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.axvline(1.0, color="C3", lw=1.2, ls=":", label="used $\\alpha$=1")
    ax.set_xlabel("Angstrom exponent $\\alpha$ (aerosol 532$\\to$1064 nm)")
    ax.set_ylabel("Mini-MPL median relative bias vs CHM15k [%]")
    ax.grid(alpha=0.3); ax.legend()
    ax.set_title("Palaiseau Mini-MPL: sensitivity of the validation bias to the assumed Angstrom exponent"
                 f"  (median aerosol fraction in band: {res['minimpl_alpha']['f_aer_med']:.2f})")
    fig.tight_layout()
    fig.savefig(OUT / "fig_report_minimpl_alpha.png", dpi=160); plt.close(fig)
    print("minimpl_alpha:", {k2: v for k2, v in res["minimpl_alpha"].items() if k2 != "alphas" and k2 != "medrelbias_pct"})


def _to_jsonable(x):
    if isinstance(x, dict):
        return {k: _to_jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_to_jsonable(v) for v in x]
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return v if np.isfinite(v) else None
    if isinstance(x, (np.integer, int)):
        return int(x)
    return x


def main(mode="csv"):
    res = {}
    if RESULTS.is_file():
        try:
            res = json.loads(RESULTS.read_text())
        except Exception:
            res = {}
    payerne_cl61(res)
    cl51_oscillation(res)
    if mode == "full":
        daynight_seasonal(res)
        altitude_bands(res)
        monthly_bias(res)
        minimpl_alpha(res)
    RESULTS.write_text(json.dumps(_to_jsonable(res), indent=1))
    print("wrote", RESULTS)
    print("DISCREPANCY_ANALYSIS_DONE", flush=True)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    main(sys.argv[1] if len(sys.argv) > 1 else "csv")
