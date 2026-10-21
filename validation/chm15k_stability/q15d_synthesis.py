# -*- coding: utf-8 -*-
"""q15d synthesis: network statistics and figures for the atmosphere / Rayleigh-method factors.

Reads the per-unit and per-night output of q15d_atmosphere.py and writes
  q15d_factors.csv           one row per factor (network median, IQR, pooled sensitivities)
  q15d_monthly.csv           network monthly composites (observed and predicted ln-biases)
  q15d_composites.png        monthly composites (landscape)
  q15d_profiles.png          seasonal normalised-ratio profiles R(z) inside/around the window (altitude on Y)
  q15d_amplitudes.png        in-phase annual amplitude per factor vs the observed cycle
Run:  python validation/chm15k_stability/q15d_synthesis.py
"""
from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import OUT  # noqa: E402
from q15d_atmosphere import EPOCH, PAYERNE, Q15, anomaly, doy_of, harmonic  # noqa: E402

warnings.filterwarnings("ignore")
PEAK_C = 176.0


def med_iqr(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan, np.nan, 0
    return float(np.median(x)), float(np.percentile(x, 25)), float(np.percentile(x, 75)), int(x.size)


def inphase(amp, peak):
    """Component of an annual harmonic in phase with the C_L cycle (peak DOY 176)."""
    return np.asarray(amp, float) * np.cos(2 * np.pi * (np.asarray(peak, float) - PEAK_C) / 365.25)


def circ_median_peak(peaks):
    p = np.asarray(peaks, float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return np.nan
    ang = 2 * np.pi * p / 365.25
    return float((np.arctan2(np.median(np.sin(ang)), np.median(np.cos(ang))) * 365.25 / (2 * np.pi)) % 365.25)


def read_csv(path):
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return rows


def col(rows, k):
    return np.array([r.get(k, np.nan) if isinstance(r.get(k, np.nan), float) else np.nan for r in rows], float)


def pooled_slope(units, series_fn, min_n=30):
    """Pooled OLS slope of kappa-normalised ln C anomalies on a series' anomalies (31-day running median)."""
    xs, ys = [], []
    for U in units:
        s = series_fn(U)
        if s is None:
            continue
        xa, ya = anomaly(U["t"], s), anomaly(U["t"], U["lc_T"])
        ok = np.isfinite(xa) & np.isfinite(ya)
        if ok.sum() >= min_n:
            xs.append(xa[ok])
            ys.append(ya[ok])
    if not xs:
        return np.nan, np.nan, np.nan, 0
    x, y = np.concatenate(xs), np.concatenate(ys)
    keep = np.abs(y - np.median(y)) < 5 * 1.4826 * np.median(np.abs(y - np.median(y)))
    x, y = x[keep], y[keep]
    X = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    se = np.sqrt(res @ res / (len(y) - 2) / np.sum((x - x.mean()) ** 2))
    r = np.corrcoef(x, y)[0, 1]
    return float(coef[1]), float(se), float(r), int(len(y))


def monthly(units, series_fn):
    """Network median (and IQR) of the per-unit monthly mean of a series minus the unit's annual mean."""
    mat = []
    for U in units:
        s = series_fn(U)
        if s is None:
            continue
        mon = np.array([(EPOCH.toordinal() + int(x)) for x in U["t"]])
        from datetime import date
        mo = np.array([date.fromordinal(int(m)).month for m in mon])
        mm = np.array([np.nanmean(s[mo == k]) if np.sum((mo == k) & np.isfinite(s)) >= 3 else np.nan
                       for k in range(1, 13)])
        if np.isfinite(mm).sum() >= 9:
            mat.append(mm - np.nanmean(mm))
    mat = np.array(mat)
    return (np.nanmedian(mat, 0), np.nanpercentile(mat, 25, 0), np.nanpercentile(mat, 75, 0), len(mat))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows_all = read_csv(Q15 / "q15d_units.csv")
    rows = [r for r in rows_all if r["stream"] != PAYERNE]
    pay = [r for r in rows_all if r["stream"] == PAYERNE]
    U_all = list(np.load(Q15 / "q15d_units.npy", allow_pickle=True))
    units = [U for U in U_all if U["key"] != PAYERNE and np.isfinite(U["kappa"])]
    q13 = {r["stream"]: r for r in read_csv(OUT / "q13_kalman_network.csv")}
    print(f"units {len(rows)} (+ Payerne separately), with kappa {len(units)}")

    out = []

    def add(factor, quantity, vals, unit, note=""):
        m, lo, hi, n = med_iqr(vals)
        out.append(dict(factor=factor, quantity=quantity, median=m, q25=lo, q75=hi, n=n, unit=unit, note=note))
        print(f"  {factor:14s} {quantity:44s} {m:9.3f} [{lo:8.3f} .. {hi:8.3f}] n {n:3d} {unit} {note}")

    # --- observed ------------------------------------------------------------------------------
    add("observed", "annual amplitude ln C (LSQ harmonic)", col(rows, "amp_obs"), "%")
    add("observed", "peak DOY", col(rows, "peak_obs"), "DOY", f"circular median {circ_median_peak(col(rows, 'peak_obs')):.0f}")
    add("observed", "annual amplitude after kappa (M2)", col(rows, "amp_obsT"), "%")
    add("observed", "peak DOY after kappa", col(rows, "peak_obsT"), "DOY",
        f"circular median {circ_median_peak(col(rows, 'peak_obsT')):.0f}")
    add("observed", "in-phase residual after kappa", inphase(col(rows, "amp_obsT"), col(rows, "peak_obsT")), "%")
    add("observed", "night-to-night scatter (31-d anomaly)", col(rows, "scat_obs"), "%")
    a_q13 = np.array([q13.get(r["stream"], {}).get("amp_op", np.nan) for r in rows], float)
    add("observed", "q13 Kalman amp_op (same units)", a_q13, "%")

    # --- F1 aerosol in the window --------------------------------------------------------------
    add("F1 aer-window", "CAMS mean bias ln(slope) (face value)", col(rows, "mean_aer"), "%")
    add("F1 aer-window", "CAMS mean R-1 in window (face value)", col(rows, "mean_raer"), "%")
    add("F1 aer-window", "CAMS annual amplitude (face value)", col(rows, "amp_aer"), "%")
    add("F1 aer-window", "CAMS peak DOY", col(rows, "peak_aer"), "DOY",
        f"circular median {circ_median_peak(col(rows, 'peak_aer')):.0f}")
    add("F1 aer-window", "CAMS night-to-night std (face value)", col(rows, "std_aer"), "%")
    add("F1 aer-window", "obs residual amp after removing CAMS", col(rows, "amp_obsT_minus_aer"), "%")
    b_aer, se_aer, r_aer, n_aer = pooled_slope(units, lambda U: U["rec"]["b_aer"])
    add("F1 aer-window", "pooled transfer d lnC / d CAMS-bias", [b_aer], "-", f"se {se_aer:.3f} r {r_aer:.3f} n {n_aer}")
    add("F1 aer-window", "per-unit transfer slope", col(rows, "slope_aer"), "-")
    sc = b_aer * col(rows, "amp_aer")
    add("F1 aer-window", "CAMS x pooled transfer: amplitude", sc, "%")
    add("F1 aer-window", "CAMS x pooled transfer: in-phase", inphase(sc, col(rows, "peak_aer")), "%")
    add("F1 aer-window", "CAMS x transfer: night-to-night std", b_aer * col(rows, "std_aer"), "%")
    add("F1 aer-window", "variance share r^2 (per unit)", col(rows, "r_aer") ** 2 * np.sign(col(rows, "r_aer")), "-")

    # lidar-side shape of the normalised ratio inside the window
    add("F1 shape", "in-window d lnR/dz mean (lidar)", col(rows, "mean_g_obs"), "1/km")
    add("F1 shape", "in-window d lnR/dz annual amplitude", col(rows, "amp_g_obs"), "1/km")
    add("F1 shape", "in-window d lnR/dz peak DOY", col(rows, "peak_g_obs"), "DOY",
        f"circular median {circ_median_peak(col(rows, 'peak_g_obs')):.0f}")
    add("F1 shape", "CAMS-predicted d lnR/dz mean", col(rows, "mean_g_cams"), "1/km")
    b_g, se_g, r_g, n_g = pooled_slope(units, lambda U: U["rec"]["g_obs"])
    add("F1 shape", "pooled d lnC / d(g_obs)", [b_g], "km", f"se {se_g:.3f} r {r_g:.3f} n {n_g}")

    # AERONET route (net aerosol: in-window excess and below-window transmission together)
    b_ae, se_ae, r_ae, n_ae = pooled_slope(units, lambda U: U["rec"]["aod_aeronet"]
                                           if np.isfinite(U["rec"]["aod_aeronet"]).sum() >= 30 else None)
    add("F1+F2 AERONET", "pooled d lnC / d AOD1064 (x0.1 = %/0.1)", [b_ae * 10], "%/0.1AOD",
        f"se {se_ae * 10:.2f} r {r_ae:.3f} n {n_ae}")
    add("F1+F2 AERONET", "AOD1064 annual amplitude (accepted nights)", col(rows, "amp_aoda"), "AOD")
    add("F1+F2 AERONET", "AOD1064 peak DOY", col(rows, "peak_aoda"), "DOY",
        f"circular median {circ_median_peak(col(rows, 'peak_aoda')):.0f}")
    pa = b_ae * col(rows, "amp_aoda") * 100
    add("F1+F2 AERONET", "predicted amplitude (pooled sens x AOD amp)", pa, "%")
    add("F1+F2 AERONET", "predicted in-phase", inphase(pa, col(rows, "peak_aoda")), "%")
    add("F1+F2 AERONET", "night-to-night std (sens x AOD anomaly std)", b_ae * col(rows, "std_aoda") * 100, "%")
    add("F1+F2 AERONET", "AOD1064 mean (accepted nights)", col(rows, "mean_aoda"), "AOD")

    # --- F2 extinction below the window --------------------------------------------------------
    add("F2 below", "applied Klett correction 2 tau_K mean", col(rows, "mean_kl"), "%")
    add("F2 below", "applied Klett correction amplitude", col(rows, "amp_kl"), "%")
    add("F2 below", "applied Klett correction peak DOY", col(rows, "peak_kl"), "DOY",
        f"circular median {circ_median_peak(col(rows, 'peak_kl')):.0f}")
    add("F2 below", "CAMS tau below window mean", col(rows, "mean_tauc"), "AOD")
    # residual (uncorrected) transmission: -2 (f_below * AOD_AERONET - tau_K), per night, units with AERONET
    res_amp, res_pk, res_mean, ratio_k = [], [], [], []
    for U in units:
        R = U["rec"]
        f = np.clip(R["tau_below_cams"] / R["aod_cams"], 0, 1)
        tb = f * R["aod_aeronet"]
        resid = -2 * (tb - R["tau_klett"])
        if np.isfinite(resid).sum() >= 40:
            a, p, _, _ = harmonic(U["t"], resid)
            res_amp.append(100 * a)
            res_pk.append(p)
            res_mean.append(100 * np.nanmean(resid))
            ratio_k.append(np.nanmedian(R["tau_klett"] / tb))
    add("F2 below", "tau_K / (f_below x AOD_AERONET)", ratio_k, "-", "share of below-window AOD the Klett recovers")
    add("F2 below", "residual -2(tau_true - tau_K) mean", res_mean, "%")
    add("F2 below", "residual amplitude", res_amp, "%")
    add("F2 below", "residual in-phase", inphase(res_amp, res_pk), "%")
    add("F2 below", "residual peak DOY", res_pk, "DOY", f"circular median {circ_median_peak(res_pk):.0f}")

    # --- F3 molecular model ----------------------------------------------------------------------
    add("F3 molecular", "CAMS/US76 bias mean", col(rows, "mean_mol"), "%")
    add("F3 molecular", "annual amplitude", col(rows, "amp_mol"), "%")
    add("F3 molecular", "peak DOY", col(rows, "peak_mol"), "DOY",
        f"circular median {circ_median_peak(col(rows, 'peak_mol')):.0f}")
    add("F3 molecular", "in-phase", inphase(col(rows, "amp_mol"), col(rows, "peak_mol")), "%")
    add("F3 molecular", "night-to-night std", col(rows, "std_mol"), "%")
    b_m, se_m, r_m, n_m = pooled_slope(units, lambda U: U["rec"]["b_mol"])
    add("F3 molecular", "pooled d lnC / d bias (expect 1)", [b_m], "-", f"se {se_m:.3f} r {r_m:.3f} n {n_m}")
    add("F3 molecular", "obs residual amp after removing mol", col(rows, "amp_obsT_minus_mol"), "%")
    add("F3 molecular", "trend", col(rows, "trend_mol"), "%/yr")

    # --- F4 sampling / selection ---------------------------------------------------------------
    add("F4 window", "window centre mean", col(rows, "mean_zc"), "km")
    add("F4 window", "window centre annual amplitude", col(rows, "amp_zc"), "km")
    add("F4 window", "window centre peak DOY", col(rows, "peak_zc"), "DOY",
        f"circular median {circ_median_peak(col(rows, 'peak_zc')):.0f}")
    b_z, se_z, r_z, n_z = pooled_slope(units, lambda U: 0.5 * (U["zb"] + U["zt"]) / 1000.0)
    add("F4 window", "pooled d lnC / d z_centre", [b_z * 100], "%/km", f"se {se_z * 100:.2f} r {r_z:.3f} n {n_z}")
    add("F4 window", "per-unit d lnC / d z_centre", col(rows, "slope_zc") * 100, "%/km")
    pz = col(rows, "slope_zc") * col(rows, "amp_zc") * 100
    add("F4 window", "per-unit slope x amplitude", pz, "%")
    add("F4 window", "per-unit slope x amplitude in-phase", inphase(pz, col(rows, "peak_zc")), "%")
    add("F4 window", "pooled slope x amplitude", b_z * col(rows, "amp_zc") * 100, "%")
    add("F4 window", "night-to-night (|slope| x zc anomaly std)", np.abs(col(rows, "slope_zc")) * col(rows, "std_zc") * 100,
        "%")
    add("F4 window", "variance share r^2 (per unit)", col(rows, "r_zc") ** 2, "-")
    add("F4 night", "night length mean", col(rows, "mean_dur"), "h")
    add("F4 night", "night length annual amplitude", col(rows, "amp_dur"), "h")
    add("F4 night", "per-unit d lnC / d(night length)", col(rows, "slope_dur") * 100, "%/h")
    pdur = col(rows, "slope_dur") * col(rows, "amp_dur") * 100
    add("F4 night", "per-unit slope x amplitude in-phase", inphase(pdur, (col(rows, "peak_dur") + 0)), "%")
    add("F4 night", "clear 10-min blocks JJA", col(rows, "nblk_jja"), "-")
    add("F4 night", "clear 10-min blocks DJF", col(rows, "nblk_djf"), "-")
    add("F4 night", "scatter JJA", col(rows, "scat_jja"), "%")
    add("F4 night", "scatter DJF", col(rows, "scat_djf"), "%")
    add("F4 select", "accepted fraction of nights JJA", col(rows, "frac_acc_jja"), "-")
    add("F4 select", "accepted fraction of nights DJF", col(rows, "frac_acc_djf"), "-")
    dsel = (col(rows, "aer_acc_jja") - col(rows, "aer_acc_djf")) - (col(rows, "aer_all_jja") - col(rows, "aer_all_djf"))
    add("F4 select", "CAMS JJA-DJF: accepted minus all nights", dsel, "% (CAMS)")
    add("F4 select", "same x pooled transfer", b_aer * dsel, "%")
    add("F4 select", "CAMS all-night amplitude (fixed window)", col(rows, "amp_aer_allnights"), "% (CAMS)")
    add("F4 select", "CAMS accepted-night amplitude", col(rows, "amp_aer"), "% (CAMS)")
    add("F4 noise", "relative uncertainty u/C amplitude", col(rows, "amp_urel") * 100, "%")
    add("F4 noise", "per-unit d lnC / d(u/C)", col(rows, "slope_urel"), "-")

    # cross-unit: does the residual cycle scale with the atmosphere?
    from scipy.stats import spearmanr
    res_ip = inphase(col(rows, "amp_obsT"), col(rows, "peak_obsT"))
    for name, x in (("CAMS in-phase aerosol", inphase(col(rows, "amp_aer"), col(rows, "peak_aer"))),
                    ("AERONET AOD amplitude", col(rows, "amp_aoda")),
                    ("window-centre amplitude", col(rows, "amp_zc")),
                    ("latitude", col(rows, "lat")),
                    ("molecular in-phase", inphase(col(rows, "amp_mol"), col(rows, "peak_mol")))):
        ok = np.isfinite(x) & np.isfinite(res_ip)
        if ok.sum() >= 10:
            rho, p = spearmanr(x[ok], res_ip[ok])
            add("cross-unit", f"Spearman residual in-phase vs {name}", [rho], "-", f"p {p:.3g} n {ok.sum()}")

    # Payerne separately
    for r in pay:
        print(f"  Payerne A: amp obs {r['amp_obs']:.1f} % pk {r['peak_obs']:.0f}, mol {r['amp_mol']:.1f} % pk "
              f"{r['peak_mol']:.0f}, CAMS aer {r['amp_aer']:.1f} %, zc amp {r['amp_zc']:.2f} km slope "
              f"{r['slope_zc'] * 100:+.1f} %/km, Klett {r['amp_kl']:.1f} %")
        out.append(dict(factor="Payerne A (excluded)", quantity="amp_obs/amp_mol/amp_aer/amp_zc/slope_zc",
                        median=r["amp_obs"], q25=r["amp_mol"], q75=r["amp_aer"], n=1, unit="%",
                        note=f"zc amp {r['amp_zc']:.2f} km, slope {r['slope_zc'] * 100:+.1f} %/km"))

    with open(Q15 / "q15d_factors.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    # --- pairs -------------------------------------------------------------------------------------
    if (Q15 / "q15d_pairs.csv").exists():
        print("pairs:")
        for r in read_csv(Q15 / "q15d_pairs.csv"):
            print(f"  {r['pair']:30s} n {r['n']:4.0f} amp A {r['amp_A']:5.1f} B {r['amp_B']:5.1f} ratio {r['amp_ratio']:5.1f} | "
                  f"kappa-norm A {r['amp_A_T']:5.1f} B {r['amp_B_T']:5.1f} ratio {r['amp_ratio_T']:5.1f} | "
                  f"anom corr {r['corr_anom_T']:+.2f}  std ratio-anom {r['std_ratio_anom']:.1f} vs A {r['std_A_anom']:.1f}")

    # --- monthly composites ----------------------------------------------------------------------
    comp = {
        "observed ln C": lambda U: 100 * U["lc"],
        "observed ln C, kappa-normalised": lambda U: 100 * U["lc_T"],
        "F3 molecular (CAMS vs US76)": lambda U: 100 * U["rec"]["b_mol"],
        "F1 CAMS in-window aerosol x transfer": lambda U: 100 * b_aer * U["rec"]["b_aer"],
        "F2 applied Klett 2 tau_K": lambda U: 100 * 2 * U["rec"]["tau_klett"],
        "F4 window centre x pooled slope": lambda U: 100 * b_z * 0.5 * (U["zb"] + U["zt"]) / 1000.0,
        "F1+F2 AERONET AOD x pooled sens": lambda U: (100 * b_ae * U["rec"]["aod_aeronet"]
                                                     if np.isfinite(U["rec"]["aod_aeronet"]).sum() >= 30 else None),
    }
    shape = {"lidar d lnR/dz in window [1/km]": lambda U: U["rec"]["g_obs"],
             "window centre [km]": lambda U: 0.5 * (U["zb"] + U["zt"]) / 1000.0,
             "CAMS in-window R-1 (face value) [-]": lambda U: U["rec"]["r_aer"]}
    mrows = []
    res = {}
    for name, fn in {**comp, **shape}.items():
        m, lo, hi, n = monthly(units, fn)
        res[name] = (m, lo, hi, n)
        for k in range(12):
            mrows.append(dict(series=name, month=k + 1, median=m[k], q25=lo[k], q75=hi[k], n_units=n))
    with open(Q15 / "q15d_monthly.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(mrows[0].keys()))
        w.writeheader()
        w.writerows(mrows)

    months = np.arange(1, 13)
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.6))
    ax = axs[0]
    for name, color in (("observed ln C", "k"), ("observed ln C, kappa-normalised", "0.5")):
        m, lo, hi, n = res[name]
        ax.plot(months, m, color=color, lw=2.2, label=f"{name} (n={n})")
        ax.fill_between(months, lo, hi, color=color, alpha=0.12)
    for name, color in (("F3 molecular (CAMS vs US76)", "#1f77b4"), ("F1 CAMS in-window aerosol x transfer", "#d62728"),
                        ("F2 applied Klett 2 tau_K", "#2ca02c"), ("F4 window centre x pooled slope", "#9467bd"),
                        ("F1+F2 AERONET AOD x pooled sens", "#ff7f0e")):
        m, lo, hi, n = res[name]
        ax.plot(months, m, color=color, lw=1.6, label=f"{name} (n={n})")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(months)
    ax.set_xlabel("month")
    ax.set_ylabel("monthly deviation from annual mean [% of C_L]")
    ax.set_title("Network monthly composite: observed cycle vs predicted atmospheric/method terms", fontsize=10)
    ax.legend(fontsize=7, loc="lower center")
    ax.grid(alpha=0.3)
    ax = axs[1]
    for name, color in (("lidar d lnR/dz in window [1/km]", "#8c564b"),):
        m, lo, hi, n = res[name]
        ax.plot(months, m, color=color, lw=2, label=f"{name} (n={n})")
        ax.fill_between(months, lo, hi, color=color, alpha=0.15)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(months)
    ax.set_xlabel("month")
    ax.set_ylabel("deviation from annual mean [1/km]")
    ax.set_title("Lidar's own in-window profile slope (aerosol would make it more negative)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axs[2]
    m, lo, hi, n = res["window centre [km]"]
    ax.plot(months, m, color="#9467bd", lw=2, label=f"window centre (n={n})")
    ax.fill_between(months, lo, hi, color="#9467bd", alpha=0.15)
    ax.set_xticks(months)
    ax.set_xlabel("month")
    ax.set_ylabel("window-centre height, deviation from mean [km]")
    ax2 = ax.twinx()
    m, lo, hi, n = res["CAMS in-window R-1 (face value) [-]"]
    ax2.plot(months, m, color="#d62728", lw=1.6, ls="--", label="CAMS R-1 in window (face value)")
    ax2.set_ylabel("CAMS R-1 deviation [-]", color="#d62728")
    ax.set_title("Window height and CAMS in-window aerosol, seasonal", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax2.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q15 / "q15d_composites.png", dpi=130)
    plt.close(fig)

    # --- seasonal normalised-ratio profiles (altitude on Y) --------------------------------------
    seasons = {"DJF": (335, 60), "MAM": (60, 152), "JJA": (152, 244), "SON": (244, 335)}
    zrel = np.arange(-2.0, 3.01, 0.25)            # km relative to the window bottom
    prof = {s: [] for s in seasons}
    for U in units:
        Rp, zb_ = U["Rprof"], U["zbins"]
        d = doy_of(U["t"])
        for s, (a, b) in seasons.items():
            sm = ((d >= a) | (d < b)) if a > b else ((d >= a) & (d < b))
            rows_s = []
            for i in np.where(sm)[0]:
                if not np.isfinite(Rp[i]).any():
                    continue
                rr = np.interp(zrel, (zb_ - U["zb"][i]) / 1000.0, np.log(np.clip(Rp[i], 1e-3, None)),
                               left=np.nan, right=np.nan)
                rows_s.append(rr)
            if len(rows_s) >= 5:
                prof[s].append(np.nanmedian(rows_s, 0))
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.6))
    ax = axs[0]
    for s, color in zip(seasons, ("#1f77b4", "#2ca02c", "#d62728", "#ff7f0e")):
        if prof[s]:
            m = 100 * (np.exp(np.nanmedian(prof[s], 0)) - 1)
            ax.plot(m, zrel, color=color, lw=2, label=f"{s} (n={len(prof[s])} units)")
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlim(-20, 60)
    ax.set_ylabel("height above the window bottom z_b [km]")
    ax.set_xlabel("normalised ratio R - 1 = rcs / (C_L beta_mol T^2) - 1 [%]  (network median)")
    ax.set_title("Clear-sky signal vs US76 molecular, per night normalised by its own C_L", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axs[1]
    if prof["JJA"] and prof["DJF"]:
        dj = np.nanmedian(prof["JJA"], 0) - np.nanmedian(prof["DJF"], 0)
        ax.plot(100 * dj, zrel, color="k", lw=2, label="JJA - DJF")
    if prof["MAM"] and prof["SON"]:
        ax.plot(100 * (np.nanmedian(prof["MAM"], 0) - np.nanmedian(prof["SON"], 0)), zrel, color="0.5", lw=1.5,
                label="MAM - SON")
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlim(-20, 40)
    ax.set_ylabel("height above the window bottom z_b [km]")
    ax.set_xlabel("seasonal difference of ln R [%]")
    ax.set_title("Summer minus winter shape: inside the window (z >= 0) it stays flat", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q15 / "q15d_profiles.png", dpi=130)
    plt.close(fig)

    # --- in-phase amplitude summary -------------------------------------------------------------
    labels, vals = [], []
    items = [("observed (raw)", inphase(col(rows, "amp_obs"), col(rows, "peak_obs"))),
             ("observed after kappa", res_ip),
             ("F1 CAMS x transfer", inphase(sc, col(rows, "peak_aer"))),
             ("F1+F2 AERONET x sens", inphase(pa, col(rows, "peak_aoda"))),
             ("F2 Klett residual", inphase(res_amp, res_pk)),
             ("F3 molecular", inphase(col(rows, "amp_mol"), col(rows, "peak_mol"))),
             ("F4 window height", inphase(pz, col(rows, "peak_zc")))]
    fig, ax = plt.subplots(figsize=(13, 5))
    for k, (lab, v) in enumerate(items):
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        if v.size:
            ax.boxplot(v, positions=[k], widths=0.6, showfliers=False)
            ax.scatter(np.full(v.size, k) + np.random.default_rng(0).uniform(-0.2, 0.2, v.size), v, s=6, alpha=0.4)
        labels.append(f"{lab}\n(n={v.size})")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(items)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("annual amplitude in phase with C_L (peak DOY 176) [%]")
    ax.set_title("Atmosphere / Rayleigh-method factors: in-phase annual amplitude per unit", fontsize=10)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(Q15 / "q15d_amplitudes.png", dpi=130)
    plt.close(fig)
    print("-> q15d_factors.csv, q15d_monthly.csv, q15d_composites.png, q15d_profiles.png, q15d_amplitudes.png")


if __name__ == "__main__":
    main()
