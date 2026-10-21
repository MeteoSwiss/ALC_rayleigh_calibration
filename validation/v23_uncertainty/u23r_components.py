# -*- coding: utf-8 -*-
"""Review attack: do the budget components double-count, and is each one right on its own?

(a) ensemble vs the new terms: the v2.2 ensemble is split into its three dimensions (window shift,
    lidar ratio, time subsets; each with the other two at nominal). The shift dimension sees the
    local gradient (doc: robust sigma ~0.15 g), the time subsets see the photon noise (doc: ~0.65 of
    it) -- both are counted again by sigma_vert and sigma_noise.
(b) the photon-noise term (claim 7, the sqrt(k) correction): the split-night constants give an
    independent realisation of the noise of C_L, (C_odd - C_even)/2 (odd vs even 300-s profiles; the
    atmosphere is common to both halves). Its scatter divided by sigma_noise must be ~1. Also the
    measured-noise variant (noise of the night mean measured on the binned stack).
(c) |g| contains the noise of the gradient: (g_odd - g_even)/2 over the shipped sub-window range is a
    realisation of that noise; compared with |g| and with the formal g_sigma.
(d) molecular term vs the time-series truth: the 31-day median removes the slowly varying part of
    the US-1976 bias b, which the budget carries in full (sqrt(b^2 + 1.35 %^2)).

Usage: python u23r_components.py [orig|fresh] [eprof_v2.2|eprof_v2.3]
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u23r_common as C  # noqa: E402
import u23_analysis as A  # noqa: E402
from u23r_truth import prepare  # noqa: E402

warnings.filterwarnings("ignore")


def rsd5(v):
    """Robust sigma of a 5-member dimension, as the ensemble computes it (max IQR/1.349, MAD/0.6745)."""
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 3:
        return np.nan
    q25, q75 = np.percentile(v, [25, 75])
    mad = np.median(np.abs(v - np.median(v)))
    return max((q75 - q25) / 1.349, mad / 0.6745)


def shipped_cfg():
    return A.shipped_cfg()


def split_gradients(d):
    """g of the odd / even halves over the SHIPPED sub-window range (recorded centres)."""
    from calibration.rayleigh.uncertainty import fit_log_gradient
    cfg = shipped_cfg()
    go, ge = np.full(len(d), np.nan), np.full(len(d), np.nan)
    for i, (z, r, s, so, se, zb, zt) in enumerate(zip(d._z, d._r, d._s, d._sub_odd, d._sub_even,
                                                      d.zb_agl, d.zt_agl)):
        if z.size == 0 or so.size != z.size or se.size != z.size:
            continue
        zz, _, ss = A.select_subwindows(z, r, s, float(zb), float(zt), cfg)
        keep = np.isin(z, zz)
        for arr, out in ((so, go), (se, ge)):
            f = fit_log_gradient(z[keep], arr[keep], np.sqrt(2.0) * s[keep], 0.5 * (zb + zt),
                                 cfg["floor"])
            out[i] = f["g"]
    return go, ge


def main():
    import pandas as pd
    set_name = sys.argv[1] if len(sys.argv) > 1 else "orig"
    method = sys.argv[2] if len(sys.argv) > 2 else "eprof_v2.2"
    df = prepare(set_name, method)["df"]
    d = df[df.ok & np.isfinite(df.sig23)].copy()
    lines = []

    def say(s):
        print(s)
        lines.append(s)

    # (a) ensemble dimensions --------------------------------------------------------------------
    d["s_shift"] = [rsd5(v) for v in d._ens_shift]
    d["s_lr"] = [rsd5(v) for v in d._ens_lr]
    d["s_time"] = [rsd5(v) for v in d._ens_time]
    for typ, g in d.groupby("type"):
        ok = np.isfinite(g.s_shift) & np.isfinite(g.def_g)
        r_shift_g = np.median(g.s_shift[ok] / np.abs(g.def_g[ok]).clip(1e-4))
        say(f"[{typ}] n {len(g)}: median sigma_ens {100 * g.sig_ens.median():.2f} %; dimensions "
            f"(robust sigma, others nominal): shift {100 * g.s_shift.median():.2f} %, LR "
            f"{100 * g.s_lr.median():.2f} %, time {100 * g.s_time.median():.2f} %")
        say(f"   shift sigma / |g| (per km): median {r_shift_g:.3f} (doc: ~0.15); corr(shift sigma,"
            f" |g|) = {np.corrcoef(g.s_shift[ok], np.abs(g.def_g[ok]))[0, 1]:+.2f}")
        okn = np.isfinite(g.s_time) & (g.def_sig_noise > 0)
        say(f"   time-subset sigma / sigma_noise: median {np.median(g.s_time[okn] / g.def_sig_noise[okn]):.2f}"
            f" (doc: ~0.65 if the subsets saw only the noise)")
        # double-counted variance as a fraction of the total budget variance
        dc_shift = (g.s_shift ** 2).clip(upper=g.sig_ens ** 2)
        dc_time = np.minimum(g.s_time ** 2, g.def_sig_noise ** 2)
        frac = (dc_shift + dc_time) / g.sig23 ** 2
        say(f"   overlap (shift part + time part counted twice) / total variance: median "
            f"{100 * np.nanmedian(frac):.1f} %, p90 {100 * np.nanpercentile(frac, 90):.1f} %; total "
            f"sigma without the overlap: median ratio "
            f"{np.nanmedian(np.sqrt(np.clip(g.sig23 ** 2 - dc_shift - dc_time, 0, None)) / g.sig23):.3f}")

    # (b) photon noise vs split-night constants -----------------------------------------------
    d["n_split"] = 0.5 * (d.c_odd - d.c_even)          # relative (constants stored / C_L)
    d["h_split"] = 0.5 * (d.c_h1 - d.c_h2)
    rows = []
    for typ, g in d.groupby("type"):
        for kcls, gg in g.groupby(pd.cut(g.k_time, [0, 5, 12, 18, 30, 1e9],
                                         labels=["k<=5", "5<k<=12", "12<k<=18", "18<k<=30", "k>30"]),
                                  observed=True):
            ok = np.isfinite(gg.n_split) & (gg.def_sig_noise > 0)
            if ok.sum() < 10:
                continue
            zp = gg.n_split[ok] / gg.def_sig_noise[ok]
            zm = gg.n_split[ok] / gg.sig_noise_meas[ok]
            zh = gg.h_split[ok] / gg.def_sig_noise[ok]
            rows.append(dict(type=typ, k_class=kcls, n=int(ok.sum()), k_med=gg.k_time.median(),
                             sig_noise_med=100 * gg.def_sig_noise.median(),
                             sig_noise_meas_med=100 * gg.sig_noise_meas.median(),
                             split_rstd=100 * C.rstd(gg.n_split[ok]),
                             pull_noise=C.rstd(zp), pull_noise_meas=C.rstd(zm),
                             pull_halves=C.rstd(zh),
                             ratio_prop_meas=float(np.nanmedian(gg.noise_ratio_36))))
    nt = pd.DataFrame(rows)
    nt.to_csv(C.OUT / f"u23r_noise_split_{set_name}_{method}.csv", index=False)
    pd.set_option("display.width", 250)
    say(nt.round(3).to_string(index=False))

    # (c) gradient noise -------------------------------------------------------------------------
    go, ge = split_gradients(d)
    d["g_split"] = 0.5 * (go - ge)
    for typ, g in d.groupby("type"):
        ok = np.isfinite(g.g_split) & np.isfinite(g.def_g)
        if ok.sum() < 10:
            continue
        sg_emp = C.rstd(g.g_split[ok])
        say(f"[{typ}] gradient noise from odd/even halves: robust sigma {100 * sg_emp:.2f} %/km "
            f"(formal g_sigma median {100 * g.def_gs[ok].median():.2f} %/km); median |g| "
            f"{100 * np.median(np.abs(g.def_g[ok])):.2f} %/km; noise share of E[g^2]: "
            f"{100 * sg_emp ** 2 / np.mean(g.def_g[ok] ** 2):.1f} % (mean), median per night "
            f"{100 * np.median((g.g_split[ok] ** 2) / np.maximum(g.def_g[ok] ** 2, 1e-12)):.1f} %")

    # (d) molecular term vs the time-series truth ------------------------------------------------
    m = d.type.eq("CHM15k") & d.key.ne(C.PAYERNE) & np.isfinite(d.anom)
    g = d[m]
    s_nob = C.total(g, g.def_sig_vert, mol_col="def_sig_mol") ** 2 - g.def_sig_mol ** 2 + 0.0135 ** 2
    p_ship = C.pull_stats(g.anom / g.sig23)
    p_nob = C.pull_stats(g.anom / np.sqrt(s_nob))
    say(f"[CHM15k] molecular term: median sigma_mol {100 * g.def_sig_mol.median():.2f} % (|b| median "
        f"{100 * np.abs(g.mol_bias).median():.2f} %); pull with sqrt(b^2+1.35^2) {p_ship['rstd']:.3f} "
        f"vs with 1.35 % only {p_nob['rstd']:.3f}")
    # is the climatological bias visible in the anomalies? (it should not be: the 31-day median
    # removes its slow part)
    ok = np.isfinite(g.mol_bias)
    say(f"   corr(anomaly, expected bias b) = {np.corrcoef(g.anom[ok], g.mol_bias[ok])[0, 1]:+.3f} "
        f"(n {int(ok.sum())})")
    (C.OUT / f"u23r_components_{set_name}_{method}.txt").write_text("\n".join(lines),
                                                                      encoding="utf-8")
    d.drop(columns=[c for c in d.columns if c.startswith("_")]).to_csv(
        C.OUT / f"u23r_components_nights_{set_name}_{method}.csv", index=False)


if __name__ == "__main__":
    main()
