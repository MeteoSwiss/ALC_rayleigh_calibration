# -*- coding: utf-8 -*-
"""Review attack: is the budget calibrated as a FUNCTION of the reported sigma (reliability)?

A pooled pull of ~1 can hide over-coverage of quiet nights and under-coverage of noisy ones. The
nights are binned by the reported sigma (and by season, tier, noise level, window height, unit) and
the observed anomaly scatter of each bin is compared with the reported sigma of the bin. The truth's
own error (standard error of the running median) is shown separately.

Unit level: are the units exchangeable? Under "the budget is equally well calibrated everywhere" the
per-unit pull spreads differ only by sampling; a permutation test (pulls shuffled across units,
same n per unit) gives the expected dispersion of the per-unit robust pull std.

Usage: python u23r_reliability.py [orig|fresh] [eprof_v2.2|eprof_v2.3]
Writes u23r_reliability_<set>_<method>.csv / _units.csv and the figure u23r_reliability_<...>.png.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u23r_common as C  # noqa: E402
from u23r_truth import prepare  # noqa: E402

warnings.filterwarnings("ignore")


def bin_table(d, by, labels=None):
    """Per bin: n, median reported sigma (v2.2 / v2.3), observed robust sigma and RMS of the anomaly,
    robust pull std and fraction outside 2 sigma (v2.3), median truth error."""
    import pandas as pd
    rows = []
    for key, g in d.groupby(by, observed=True):
        if len(g) < 15:
            continue
        a = g.anom.values
        z23 = a / g.sig23.values
        z22 = a / g.sig_ens.values
        # bootstrap 68 % interval of the observed robust sigma
        rs = np.random.default_rng(0)
        boot = [C.rstd(rs.choice(a, a.size)) for _ in range(300)]
        rows.append(dict(bin=labels.get(key, key) if labels else key, n=len(g),
                         rep22=100 * np.median(g.sig_ens), rep23=100 * np.median(g.sig23),
                         obs_rstd=100 * C.rstd(a), obs_lo=100 * np.percentile(boot, 16),
                         obs_hi=100 * np.percentile(boot, 84),
                         obs_rms=100 * float(np.sqrt(np.mean(a ** 2))),
                         truth_se=100 * np.nanmedian(g.anom_se),
                         pull22=C.pull_stats(z22)["rstd"], pull23=C.pull_stats(z23)["rstd"],
                         out23=C.pull_stats(z23)["out2"],
                         out23_hi=float(np.mean(z23 > 2)), out23_lo=float(np.mean(z23 < -2))))
    return pd.DataFrame(rows)


def unit_dispersion(d, n_perm=2000, seed=1):
    """Observed spread of the per-unit robust pull std vs its permutation distribution."""
    import pandas as pd
    z = (d.anom / d.sig23).values
    keys = d.key.values
    units = [k for k in np.unique(keys) if np.sum(keys == k) >= 15]
    m = np.isin(keys, units)
    z, keys = z[m], keys[m]
    sizes = {k: int(np.sum(keys == k)) for k in units}
    obs = {k: C.pull_stats(z[keys == k])["rstd"] for k in units}
    obs_sd = float(np.std(np.log([obs[k] for k in units])))
    rs = np.random.default_rng(seed)
    sims = []
    per_unit = {k: [] for k in units}
    for _ in range(n_perm):
        p = rs.permutation(z)
        i0 = 0
        vals = []
        for k in units:
            v = C.pull_stats(p[i0:i0 + sizes[k]])["rstd"]
            per_unit[k].append(v)
            vals.append(v)
            i0 += sizes[k]
        sims.append(np.std(np.log(vals)))
    sims = np.array(sims)
    rows = []
    for k in units:
        sv = np.array(per_unit[k])
        rows.append(dict(key=k, n=sizes[k], pull23=obs[k], perm_p05=np.percentile(sv, 5),
                         perm_p95=np.percentile(sv, 95),
                         outside=bool(obs[k] < np.percentile(sv, 2.5) or obs[k] > np.percentile(sv, 97.5))))
    t = pd.DataFrame(rows).sort_values("pull23")
    return t, obs_sd, sims


def main():
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    set_name = sys.argv[1] if len(sys.argv) > 1 else "orig"
    method = sys.argv[2] if len(sys.argv) > 2 else "eprof_v2.2"
    df = prepare(set_name, method)["df"]
    d = df[df.ok & np.isfinite(df.anom) & np.isfinite(df.sig23) & (df.type == "CHM15k")
           & (df.key != C.PAYERNE)].copy()
    tag = f"{set_name}_{method}"
    out = []
    # (1) by reported sigma (deciles of the v2.3 sigma)
    d["q23"] = pd.qcut(d.sig23, 10, labels=False)
    t1 = bin_table(d, "q23")
    t1.insert(0, "by", "reported sigma decile")
    out.append(t1)
    # (2) season
    seas = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA",
            8: "JJA", 9: "SON", 10: "SON", 11: "SON"}
    d["season"] = d.month.map(seas)
    t2 = bin_table(d, "season")
    t2.insert(0, "by", "season")
    out.append(t2)
    # (3) tier
    t3 = bin_table(d, "tier", {1.0: "tier 1", 2.0: "tier 2"})
    t3.insert(0, "by", "tier")
    out.append(t3)
    # (4) photon-noise terciles, (5) window-centre height, (6) dominant component
    d["noise_t"] = pd.qcut(d.def_sig_noise, 3, labels=["low noise", "mid noise", "high noise"])
    t4 = bin_table(d, "noise_t")
    t4.insert(0, "by", "photon noise tercile")
    out.append(t4)
    d["zc_bin"] = pd.cut(d.zc_agl, [0, 3000, 4000, 5000, 20000],
                         labels=["< 3 km", "3-4 km", "4-5 km", "> 5 km"])
    t5 = bin_table(d, "zc_bin")
    t5.insert(0, "by", "window centre AGL")
    out.append(t5)
    comp = np.vstack([d.sig_ens ** 2, d.def_sig_vert ** 2, d.def_sig_mol ** 2 + d.def_sig_aer ** 2,
                      d.def_sig_noise ** 2]).T
    names = np.array(["ensemble", "vertical", "molecular+aerosol", "noise"])
    d["dominant"] = names[np.nanargmax(np.nan_to_num(comp), axis=1)]
    t6 = bin_table(d, "dominant")
    t6.insert(0, "by", "dominant component")
    out.append(t6)
    tab = pd.concat(out, ignore_index=True)
    tab.to_csv(C.OUT / f"u23r_reliability_{tag}.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(tab.round(3).to_string(index=False))
    # units
    ut, obs_sd, sims = unit_dispersion(d)
    ut.to_csv(C.OUT / f"u23r_reliability_units_{tag}.csv", index=False)
    print(ut.round(3).to_string(index=False))
    p = float(np.mean(sims >= obs_sd))
    print(f"per-unit pull spread: sd(log pull) observed {obs_sd:.3f}, permutation median "
          f"{np.median(sims):.3f} [95 % {np.percentile(sims, 95):.3f}], p = {p:.4f}; units outside "
          f"their 95 % permutation interval: {int(ut.outside.sum())} of {len(ut)}")

    # figure: reliability (sigma deciles), categories, units
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.6))
    a = ax[0]
    for col, lab, colr, mk in (("rep22", "v2.2 ensemble", "0.5", "o"),
                               ("rep23", "v2.3 budget", "tab:blue", "s")):
        dd = d.copy()
        dd["q"] = pd.qcut(dd["sig_ens" if col == "rep22" else "sig23"], 10, labels=False)
        tt = bin_table(dd, "q")
        a.errorbar(tt[col], tt.obs_rstd, yerr=[tt.obs_rstd - tt.obs_lo, tt.obs_hi - tt.obs_rstd],
                   fmt=mk + "-", color=colr, label=f"{lab}: observed robust sigma per decile")
        a.plot(tt[col], np.sqrt(tt[col] ** 2 + tt.truth_se ** 2), mk + ":", color=colr, mfc="none",
               label=f"{lab}: reported + truth error (quadrature)")
    lim = [0.5, 40]
    a.plot(lim, lim, "k--", lw=1, label="1:1")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(lim)
    a.set_ylim(2, 40)
    a.set_xlabel("reported 1-sigma [%] (decile median)")
    a.set_ylabel("observed anomaly scatter [%] (robust sigma, 68 % bootstrap)")
    a.set_title(f"Reliability, CHM15k ({set_name} units, {method})", fontsize=10)
    a.legend(fontsize=7, loc="upper left")
    a.grid(alpha=0.3, which="both")
    a = ax[1]
    cats = tab[tab.by != "reported sigma decile"].reset_index(drop=True)
    y = np.arange(len(cats))
    a.barh(y, cats.pull23, color="tab:blue", alpha=0.6, label="v2.3 robust pull std")
    a.plot(cats.pull22, y, "o", color="0.4", label="v2.2 robust pull std")
    a.axvline(1.0, color="k", lw=1)
    a.set_yticks(y)
    a.set_yticklabels([f"{r.by}: {r.bin} (n {r.n})" for r in cats.itertuples()], fontsize=7)
    a.set_xscale("log")
    a.set_xlabel("robust pull std (time-series truth)")
    a.set_title("Pull by category", fontsize=10)
    a.legend(fontsize=7, loc="lower right")
    a.grid(alpha=0.3, axis="x", which="both")
    a = ax[2]
    a.fill_between(np.arange(len(ut)), ut.perm_p05, ut.perm_p95, color="0.85",
                   label="90 % permutation interval (units exchangeable)")
    a.plot(np.arange(len(ut)), ut.pull23, "o", color="tab:blue", label="observed per-unit pull")
    a.axhline(1.0, color="k", lw=1)
    a.set_xticks(np.arange(len(ut)))
    a.set_xticklabels([k.split("_")[0][-5:] + "_" + k.rsplit("_", 1)[1] for k in ut.key],
                      rotation=80, fontsize=6)
    a.set_ylabel("robust pull std of the unit (v2.3)")
    a.set_title(f"Per unit: spread p = {p:.3f} (permutation)", fontsize=10)
    a.legend(fontsize=7, loc="upper left")
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(C.OUT / f"u23r_reliability_{tag}.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
