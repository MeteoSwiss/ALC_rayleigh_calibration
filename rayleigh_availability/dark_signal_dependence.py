# -*- coding: utf-8 -*-
"""Is the hood-measured dark the baseline of NORMAL operation, or does the detector add a signal-
dependent baseline (recovery after a strong near-range return) that the hood does not see?

Data: the hourly, cloud-screened, unfiltered matched profiles of the paper's Appendix 6 (Payerne,
4 March - 8 October 2026, beta_att at 1064 nm, operational constants; cache_payerne_cl61dark.npz).
Clear NIGHT hours only (sun < -6 deg).

Per instrument (CL61 cloud-calibrated: constant independent of the dark; CHM15k Rayleigh):
  * Q = near-range load of the hour = integral over 0.2-2 km of (beta - beta_mol) [sr-1], the
    return the detector has just handled before the free-troposphere gates;
  * hours split into Q terciles; per tercile the median profile over 3-8 km is fitted as
        beta(z) = T * beta_mol(z) + theta * d(z)
    with d = the hood dark in beta units (exactly what the "+ dark" configuration subtracts) and T
    absorbing the calibration error and the two-way aerosol transmission (which also grows with Q).
    theta ~ 1 in every tercile  -> the hood dark is the operational baseline;
    theta growing with Q        -> an extra signal-induced (recovery) baseline the hood misses.
  * the same fit on the dark-corrected configuration (expect theta ~ 0);
  * far-field offset: median of (beta - beta_mol) / r^2 over 10-14 km, per tercile (a constant
    offset in the non-range-corrected signal, e.g. background-subtraction error).
Uncertainties: bootstrap over NIGHTS (hours of a night resampled together).

Caveat: elevated free-tropospheric aerosol (more frequent in the high-Q, summer hours) also makes
beta/beta_mol decrease with altitude and biases theta upward; the 4-8 km fit is given as a check.

Run:  python rayleigh_availability/dark_signal_dependence.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation/paper"))
from validation.paper import intercompare as IC   # noqa: E402

D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
CACHE = D / "figs_paper_validation/paper_v08/cache_payerne_cl61dark.npz"
OUT = D / "rayleigh_availability"
FIG = REPO / "doc/reports/figs_remote_dark/dark_signal_dependence.png"
LAT, LON, ALT = 46.8137, 6.9425, 491.0
NBOOT = 300
PLOT_ONLY = False
rng_ = np.random.default_rng(1)


def sun_elevation(t):
    """Solar elevation [deg] (NOAA approximation) for datetime64 times (UTC)."""
    t = np.asarray(t, "datetime64[s]")
    doy = (t.astype("datetime64[D]") - t.astype("datetime64[Y]")).astype(int) + 1
    hr = (t - t.astype("datetime64[D]")).astype(float) / 3600.0
    g = 2 * np.pi / 365 * (doy - 1 + (hr - 12) / 24)
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g)
                    - 0.040849 * np.sin(2 * g))
    ha = np.radians((hr * 60 + eqt + 4 * LON) / 4 - 180)
    la = np.radians(LAT)
    return np.degrees(np.arcsin(np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)))


def molecular_per_hour(times, z):
    """CAMS molecular attenuated backscatter at 1064 nm for each hour (CAMS period of the hour)."""
    M = np.full((times.size, z.size), np.nan)
    for sel, cf, t0, t1, _p in IC.cams_periods(np.asarray(times).astype("datetime64[ns]")):
        b = IC._mol_att_cams(cf, LAT, LON, t0, t1, z + ALT, ALT, 1064.47) if cf is not None else None
        M[sel] = b if b is not None else IC._molecular_beta(z, ALT, 1064.47)
    return M


def fit_theta(B, M, d, z, zlo, zhi):
    """Fit median profile: median(B) = T * median(M) + theta * d over [zlo, zhi]. Returns (T, theta)."""
    m = (z >= zlo) & (z <= zhi)
    with np.errstate(all="ignore"):
        b = np.nanmedian(B[:, m], axis=0)
        mol = np.nanmedian(M[:, m], axis=0)
    ok = np.isfinite(b) & np.isfinite(mol) & np.isfinite(d[m])
    A = np.vstack([mol[ok], d[m][ok]]).T
    # weight by 1/mol so every altitude counts in relative terms (the ratio form of the model)
    w = 1.0 / mol[ok]
    coef, *_ = np.linalg.lstsq(A * w[:, None], b[ok] * w, rcond=None)
    return coef


def far_offset(B, M, z):
    m = (z >= 10000) & (z <= 14000)
    with np.errstate(all="ignore"):
        return float(np.nanmedian(np.nanmedian((B[:, m] - M[:, m]) / (z[m] / 1e3) ** 2, axis=0)))


def boot(fn, nights, *arrays):
    """fn(*arrays) on the full set, plus a bootstrap over nights -> (value, std)."""
    val = np.asarray(fn(*arrays))
    un = np.unique(nights)
    idx = {n: np.where(nights == n)[0] for n in un}
    res = []
    for _ in range(NBOOT):
        pick = np.concatenate([idx[n] for n in rng_.choice(un, un.size)])
        res.append(fn(*[a[pick] for a in arrays]))
    return val, np.std(np.asarray(res), axis=0)


def main():
    c = np.load(CACHE)
    z, beta, labels = c["z"], c["beta"].astype("f8"), [str(x) for x in c["labels"]]
    times = c["time"].astype("datetime64[s]")
    night = sun_elevation(times) < -6
    clear = np.all(np.isfinite(beta[:, :, (z >= 200) & (z <= 8000)]), axis=(0, 2))
    sel = night & clear
    # the night an hour belongs to (date of the evening): shift by 12 h
    nights = (times[sel] - np.timedelta64(12, "h")).astype("datetime64[D]")
    print(f"clear night hours: {sel.sum()} in {np.unique(nights).size} nights")
    M = molecular_per_hour(times[sel], z)
    out, rows = {}, []
    plot = {}
    for name, i_none, i_dark in (("CL61 cloud", 2, 3), ("CHM15k", 0, 1)):
        Bn, Bd = beta[i_none][sel], beta[i_dark][sel]
        d = np.nanmedian(Bn - Bd, axis=0)                      # hood dark in beta units (negative)
        nr = (z >= 200) & (z <= 2000)
        Q = np.nansum((Bn[:, nr] - M[:, nr]), axis=1) * np.median(np.diff(z)) * 1e-6   # sr-1
        edges = np.nanpercentile(Q, [0, 33.3, 66.7, 100])
        terc = np.clip(np.searchsorted(edges, Q, side="right") - 1, 0, 2)
        out[name] = {}
        plot[name] = dict(d=d, terc=[])
        for k in range(3):
            s = terc == k
            r = dict(n_hours=int(s.sum()), n_nights=int(np.unique(nights[s]).size),
                     Q_median=float(np.median(Q[s])))
            for tag, B in (("none", Bn), ("dark", Bd)):
                for zlo in (3000, 4000):
                    v, e = boot(lambda b, m: fit_theta(b, m, d, z, zlo, 8000), nights[s], B[s], M[s])
                    r[f"{tag}_{zlo // 1000}-8km"] = dict(T=float(v[0]), T_sd=float(e[0]),
                                                         theta=float(v[1]), theta_sd=float(e[1]))
                v, e = boot(lambda b, m: far_offset(b, m, z), nights[s], B[s], M[s])
                r[f"{tag}_far_offset"] = (float(v), float(e))
            # three-term check: hood dark + a constant offset in the NON-range-corrected signal (c r^2,
            # c = far-field offset left after the hood dark is subtracted, measured at 10-14 km)
            def fit_c(bn, bd, m):
                cr = far_offset(bd, m, z)
                return fit_theta(bn - cr * (z[None, :] / 1e3) ** 2, m, d, z, 3000, 8000)
            v, e = boot(fit_c, nights[s], Bn[s], Bd[s], M[s])
            r["none_3-8km_farcorr"] = dict(T=float(v[0]), T_sd=float(e[0]), theta=float(v[1]), theta_sd=float(e[1]))
            out[name][f"tercile{k + 1}"] = r
            with np.errstate(all="ignore"):
                plot[name]["terc"].append((np.nanmedian(Bn[s], axis=0) / np.nanmedian(M[s], axis=0),
                                           np.nanmedian(Bd[s], axis=0) / np.nanmedian(M[s], axis=0),
                                           float(np.median(Q[s]))))
            a, b = r["none_3-8km"], r["dark_3-8km"]
            print(f"{name:10s} Q-tercile {k + 1}: Q={r['Q_median']:.2e} sr-1  {r['n_hours']:4d} h / "
                  f"{r['n_nights']:3d} nights | no corr: theta {a['theta']:+.2f}+/-{a['theta_sd']:.2f} "
                  f"T {a['T']:.3f} | dark-corr: theta {b['theta']:+.2f}+/-{b['theta_sd']:.2f} | "
                  f"4-8 km theta {r['none_4-8km']['theta']:+.2f}+/-{r['none_4-8km']['theta_sd']:.2f} | "
                  f"far offset {r['none_far_offset'][0]:+.2e}+/-{r['none_far_offset'][1]:.1e} | "
                  f"theta with far offset removed {r['none_3-8km_farcorr']['theta']:+.2f}+/-"
                  f"{r['none_3-8km_farcorr']['theta_sd']:.2f}")
    (OUT / "dark_signal_dependence.json").write_text(json.dumps(out, indent=1))
    figure(z, plot, out)


def figure(z, plot, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 4, figsize=(20, 6.5), gridspec_kw=dict(width_ratios=[1, 1, 1, 1]))
    cols = ["#9ecae1", "#3182bd", "#08306b"]
    for j, name in enumerate(("CL61 cloud", "CHM15k")):
        ax = axs[2 * j]
        p = plot[name]
        sm = lambda x: np.convolve(np.nan_to_num(x, nan=1.0), np.ones(31) / 31, mode="same")   # 450 m, display only
        for k, (rn, rd, q) in enumerate(p["terc"]):
            rn, rd = sm(rn), sm(rd)
            ax.plot(100 * (rn - 1), z / 1e3, color=cols[k], lw=1.6, label=f"no correction, Q tercile {k + 1}")
            ax.plot(100 * (rd - 1), z / 1e3, "--", color=cols[k], lw=1.2, label=f"dark subtracted, tercile {k + 1}")
        ax.axvline(0, color="k", lw=0.7)
        ax.axhspan(3, 8, color="0.92", zorder=0)
        ax.set_xlim(-100, 60)
        ax.set_ylim(2, 10)
        ax.set_xlabel(r"median $\beta_{att}$ / molecular $-$ 1  [%]")
        ax.set_ylabel("Altitude a.g.l. [km]")
        ax.set_title(f"{name}: clear nights by near-range load Q (450 m running mean)", fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="lower left")
        ax = axs[2 * j + 1]
        q = [out[name][f"tercile{k + 1}"]["Q_median"] for k in range(3)]
        for tag, mk, lab in (("none_3-8km", "o-", "no correction, fit 3-8 km"),
                             ("none_4-8km", "s--", "no correction, fit 4-8 km"),
                             ("dark_3-8km", "^-", "dark subtracted, fit 3-8 km"),
                             ("none_3-8km_farcorr", "D-", "no correction, far-field offset c r$^2$ removed")):
            th = [out[name][f"tercile{k + 1}"][tag]["theta"] for k in range(3)]
            sd = [out[name][f"tercile{k + 1}"][tag]["theta_sd"] for k in range(3)]
            ax.errorbar(np.array(q) * 1e3, th, yerr=sd, fmt=mk, capsize=3, label=lab)
        ax.axhline(1, color="k", lw=0.8)
        ax.axhline(0, color="0.5", lw=0.8, ls=":")
        ax.set_xlabel(r"near-range load Q, $\int_{0.2}^{2\,km}(\beta-\beta_{mol})dz$  [10$^{-3}$ sr$^{-1}$]")
        ax.set_ylabel(r"$\theta$ = amplitude of the hood dark in the sky profile")
        ax.set_title(f"{name}: hood-dark amplitude vs near-range load\n(1 = hood dark fully present, 0 = absent)",
                     fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.5)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=120)
    print(f"-> {FIG}")


if __name__ == "__main__":
    main()
