# -*- coding: utf-8 -*-
"""q19 - the eprof_v2.2 noise propagation: sigma_signal is sqrt(300 s / dt_native) too large.

calibration._sigma_on_fit_grid divides the native first-difference sigma by sqrt(native gates per fit
gate) and by sqrt(n_binned_profiles), but never by sqrt(native profiles per 300-s bin): the time
averaging inside each bin is lost. The night-mean noise is measured here directly on the binned stack
(robust first differences of consecutive 5-min profiles, q19_gradient_gate.empirical_sigma) and compared
with sigma_signal on every recorded night; the reduced chi-square of the forced-0 Rayleigh fit in the
selected window is recomputed with both noise estimates.
Writes q19_sigma_check.csv and q19_sigma_check.png (landscape).
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import q19_lib as L  # noqa: E402

Q = L.Q
warnings.filterwarnings("ignore")


def chi2_window(P_arr, s, e):
    """Reduced chi-square of the forced-0 fit signal = a p_mol in [s, e) with both noise estimates."""
    z, sig, p, sge, sgo = P_arr
    s, e = int(s), int(e)
    x, y = p[s:e], sig[s:e]
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(sge[s:e]) & np.isfinite(sgo[s:e])
    if ok.sum() < 5:
        return np.nan, np.nan
    x, y, se, so = x[ok], y[ok], sge[s:e][ok], sgo[s:e][ok]
    a = np.sum(x * y) / np.sum(x * x)
    r = y - a * x
    dof = max(ok.sum() - 1, 1)
    return float(np.sum((r / se) ** 2) / dof), float(np.sum((r / so) ** 2) / dof)


def main():
    df = L.load()
    op = df[df.variant == "op"].copy()
    rows = []
    for key, g in op.groupby("key"):
        for m, gm in g.groupby(g.date.str[:6]):
            f = Q.REC_DIR / key / f"{m}.npz"
            if not f.exists():
                continue
            with np.load(f) as zf:
                for r in gm.itertuples():
                    k = f"p_{r.date}"
                    if k not in zf.files or not np.isfinite(r.sel_s):
                        continue
                    ce, co = chi2_window(zf[k], r.sel_s, r.sel_e)
                    rows.append(dict(key=key, type=r.type, date=r.date, tier=r.tier, month=r.month,
                                     sigma_ratio=r.sigma_ratio, chi2_emp=ce, chi2_op=co, zs=r.sel_zs))
    S = pd.DataFrame(rows)
    S.to_csv(Q.Q19 / "q19_sigma_check.csv", index=False)
    lines = []
    for t, g in S.groupby("type"):
        lines.append(f"{t}: nights {len(g)}, sigma_op/sigma_measured median {g.sigma_ratio.median():.2f} "
                     f"[{g.sigma_ratio.quantile(.1):.2f}, {g.sigma_ratio.quantile(.9):.2f}] (10-90 %); "
                     f"window chi2red (forced-0 fit) with the measured noise: tier-1 "
                     f"{g[g.tier == 1].chi2_emp.median():.2f}, tier-2 {g[g.tier == 2].chi2_emp.median():.2f}; "
                     f"with the operational sigma: tier-1 {g[g.tier == 1].chi2_op.median():.3f}, "
                     f"tier-2 {g[g.tier == 2].chi2_op.median():.3f}; tier-2 windows with chi2red > 2.5 when "
                     f"measured: {100 * np.mean(g[g.tier == 2].chi2_emp > 2.5):.0f} %")
    # per-unit sigma ratio (native sampling differs between units)
    pu = S.groupby(["type", "key"]).sigma_ratio.median().reset_index()
    pu.to_csv(Q.Q19 / "q19_sigma_units.csv", index=False)
    open(Q.Q19 / "q19_sigma_check.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    a = ax[0]
    for t, col in (("CHM15k", "k"), ("CL61", "tab:blue"), ("CL31", "tab:orange")):
        x = S[S.type == t].sigma_ratio.dropna()
        if len(x):
            a.hist(x, bins=np.linspace(1, 6, 51), histtype="step", lw=1.5, color=col, density=True,
                   label=f"{t} (median {x.median():.2f})")
    for v, lab in ((np.sqrt(300 / 15), "sqrt(300/15)"), (np.sqrt(300 / 12), "sqrt(300/12)")):
        a.axvline(v, color="0.5", ls="--", lw=0.8)
        a.text(v, a.get_ylim()[1] * 0.95, lab, rotation=90, va="top", ha="right", fontsize=8)
    a.set_xlabel("operational sigma_signal / measured night-mean noise (2-5 km)")
    a.set_ylabel("density")
    a.set_title("eprof_v2.2 noise propagation is too large by ~sqrt(300 s / dt_native)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
    a = ax[1]
    chm = S[S.type == "CHM15k"]
    bins = np.logspace(-3, 3, 61)
    for tier, col in ((1, "tab:green"), (2, "tab:red")):
        g = chm[chm.tier == tier]
        a.hist(g.chi2_op.dropna(), bins=bins, histtype="step", color=col, ls="--", lw=1.2,
               label=f"tier {tier}, operational sigma")
        a.hist(g.chi2_emp.dropna(), bins=bins, histtype="step", color=col, lw=1.6,
               label=f"tier {tier}, measured noise")
    a.axvline(2.5, color="k", lw=1)
    a.text(2.6, a.get_ylim()[1] * 0.9, "max_chi2red 2.5", fontsize=8)
    a.set_xscale("log")
    a.set_xlabel("reduced chi-square of the window's forced-0 Rayleigh fit")
    a.set_ylabel("nights")
    a.set_title("CHM15k: the tier-2 chi2 gate is inert with the operational sigma")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)
    a = ax[2]
    # example profile: noise estimates vs height
    t2h = chm[(chm.tier == 2) & (chm.sigma_ratio > 4.0)]
    ex = t2h.iloc[len(t2h) // 2]
    P = L.profile_of(ex.key, ex.date)
    if P is not None:
        z, sig, p, sge, sgo = P
        m = (z > 300) & (z < 10000)
        a.plot(sig[m] / sge[m], z[m] / 1e3, color="k", lw=0.8, label="signal / measured noise")
        a.plot(sig[m] / sgo[m], z[m] / 1e3, color="tab:red", lw=0.8, label="signal / operational sigma")
        a.axvline(1, color="0.5", lw=0.8)
        a.set_xscale("log")
        a.set_ylabel("range above the instrument [km]")
        a.set_xlabel("per-gate SNR of the night mean (30 m)")
        a.set_title(f"example: {ex.key} {ex.date}")
        a.legend(fontsize=8)
        a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(Q.Q19 / "q19_sigma_check.png", dpi=120)


if __name__ == "__main__":
    main()
