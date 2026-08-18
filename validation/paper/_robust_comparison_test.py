"""
_robust_comparison_test.py — PROTOTYPE (Payerne only, standalone; does NOT touch the pipeline).

Compares three per-altitude-gate statistics of the test-vs-CHM15k agreement on the paired hourly
attenuated backscatter, to expose the SNR-selection bias and the coherent electronic offset:

  (1) SNR-gated median relative bias          — the CURRENT method. At altitude the SNR>=3 gate keeps
      only gates whose noisy beta exceeds 3*sigma -> conditions on positive noise -> biased HIGH,
      instrument-dependently.
  (2) ungated (NO SNR) median relative bias    — removes that selection bias (zero-mean noise no longer
      censored), but over a multi-month record it still sees the COHERENT electronic offset b(z) (the
      CL31 digitizer ripple / baseline), which does not average out and dominates where the true signal
      is weak.
  (3) errors-in-variables (Deming) regression WITH intercept, ungated, per gate:
          beta_test(z,t) = rho(z) * beta_ref(z,t) + alpha(z)
      slope     rho(z)   = calibration ratio c_test/c_ref  — robust to BOTH the SNR selection AND the
                           additive offset (the offset falls entirely into the intercept);
      intercept alpha(z) = the net electronic offset b_test - rho*b_ref (~ b_test, the CHM having ~no
                           digitizer ripple).
      Block (per-day) bootstrap CI on rho(z). Noise-variance ratio lambda(z) from a lag-1 estimator.

Usage:  python -m validation.paper._robust_comparison_test
Output: figs_paper_validation/paper_python/fig_robust_comparison_payerne.png  + printed summary.
"""
from __future__ import annotations
import os
import sys
import warnings
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from validation.paper import intercompare as IC
from validation.paper import run_paper_validation as RPV

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
ZMIN, ZMAX = 500.0, 3500.0            # extend to 3500 m — the noisy top the user flagged
N_BOOT = 200
TYPE_COLORS = {"CL31": "#ff7f0e", "CL51": "#9467bd", "CL61": "#1f77b4"}


# --------------------------------------------------------------------------- data on a common grid
def _read_gated_ungated(ch, start, end):
    """Read L1 twice: with the operational per-gate SNR>=3 gate (in _l1_day, on the native profiles),
    and WITHOUT it (IC.L1_SNR_GATE=False). Same hourly time grid; only the gate NaNs differ."""
    IC.L1_SNR_GATE = True
    l1g = IC.read_l1(ch["wmo"], ch["ident"], start, end)
    IC.L1_SNR_GATE = False
    l1u = IC.read_l1(ch["wmo"], ch["ident"], start, end)
    IC.L1_SNR_GATE = True
    return l1g, l1u


def build_payerne():
    st = RPV.BENCHMARK["payerne"]; sc = RPV.SITE["payerne"]; target = sc["target"]
    l1cache, items = {}, []
    for ch in st["channels"]:
        key = (ch["wmo"], ch["ident"])
        if key not in l1cache:
            l1cache[key] = _read_gated_ungated(ch, st["start"], st["end"])
        l1g, l1u = l1cache[key]
        if l1g is None or l1u is None:
            continue
        # channel_beta = full pipeline (overlap, calib, WV, CAMS-molecular wavelength) — beta_scr is the
        # cloud-screened hourly beta on l1["time"]. Gated vs ungated differ ONLY through the L1 SNR gate.
        rg = RPV.channel_beta(l1g, ch, target)
        ru = RPV.channel_beta(l1u, ch, target)
        if rg is None or ru is None:
            continue
        items.append(dict(ch=ch, time=l1u["time"], alt=l1u["alt"], salt=l1u["station_alt"],
                          betag=rg["beta_scr"], betau=ru["beta_scr"]))
    if not items:
        raise SystemExit("no channels")
    union = np.unique(np.concatenate([it["time"] for it in items]))
    for it in items:
        idx = {t: i for i, t in enumerate(it["time"])}
        pos = np.array([idx.get(t, -1) for t in union])
        ok = pos >= 0
        for tag, src in (("gatU", "betag"), ("ungU", "betau")):
            M = np.full((union.size, it[src].shape[1]), np.nan)
            M[ok] = it[src][pos[ok]]
            it[tag] = M
    altGrid = IC.build_common_grid([it["alt"] for it in items])
    for it in items:
        it["gatC"] = IC.regrid(it["gatU"], it["alt"], altGrid)
        it["ungC"] = IC.regrid(it["ungU"], it["alt"], altGrid)
    salt = items[0]["salt"]
    return items, altGrid, salt, union, sc["ref"]


# --------------------------------------------------------------------------- estimators
def _noise_col(Y):
    """Per-gate noise sigma(z) from the lag-1 hourly difference (robust): sigma = 1.4826*MAD(dY)/sqrt(2)."""
    d = np.diff(Y, axis=0)
    med = np.nanmedian(d, axis=0)
    mad = np.nanmedian(np.abs(d - med[None, :]), axis=0)
    return 1.4826 * mad / np.sqrt(2.0)


def deming_cols(X, Y, lam):
    """Per-gate Deming regression WITH intercept: Y = rho*X + alpha, errors in both (variance ratio
    lam = var(Y_noise)/var(X_noise)). X,Y are (n_time, n_gate); returns rho[gate], alpha[gate], n[gate]."""
    M = np.isfinite(X) & np.isfinite(Y)
    Xm = np.where(M, X, np.nan); Ym = np.where(M, Y, np.nan)
    n = M.sum(axis=0)
    with np.errstate(all="ignore"):
        xbar = np.nanmean(Xm, axis=0); ybar = np.nanmean(Ym, axis=0)
        dx = Xm - xbar[None, :]; dy = Ym - ybar[None, :]
        sxx = np.nansum(dx * dx, axis=0); syy = np.nansum(dy * dy, axis=0)
        sxy = np.nansum(dx * dy, axis=0)
        lam = np.broadcast_to(lam, sxx.shape)
        disc = (syy - lam * sxx) ** 2 + 4.0 * lam * sxy ** 2
        rho = (syy - lam * sxx + np.sqrt(disc)) / (2.0 * sxy)
        alpha = ybar - rho * xbar
    bad = (n < 10) | ~np.isfinite(sxy) | (np.abs(sxy) < 1e-30)
    rho[bad] = np.nan; alpha[bad] = np.nan
    return rho, alpha, n


def med_relbias_cols(X, Y):
    """Per-gate median relative bias 100*median((Y-X)/X) over cells with X>0 (the pipeline's statistic)."""
    out = np.full(X.shape[1], np.nan)
    for j in range(X.shape[1]):
        x = X[:, j]; y = Y[:, j]
        m = np.isfinite(x) & np.isfinite(y) & (x > 0)
        if m.sum() >= 10:
            out[j] = 100.0 * np.median((y[m] - x[m]) / x[m])
    return out


# --------------------------------------------------------------------------- main
def main():
    warnings.filterwarnings("ignore")
    items, altGrid, salt, union, iref = build_payerne()
    z_agl = altGrid - salt
    band = (z_agl >= ZMIN) & (z_agl <= ZMAX)
    zkm = z_agl[band] / 1000.0
    ref = items[iref]
    Xu = ref["ungC"][:, band]; Xg = ref["gatC"][:, band]
    days = union.astype("datetime64[D]")                    # per-hour day label for the block bootstrap
    uday = np.unique(days)

    tests = [it for k, it in enumerate(items) if k != iref and it["ch"]["itype"] in ("CL31", "CL51", "CL61")]
    # keep one entry per itype (the cloud channel) to avoid duplicate CL61 rows
    seen = {}; picked = []
    for it in tests:
        t = it["ch"]["itype"]
        if it["ch"]["calib"] == "cloud" or t not in seen:
            seen[t] = it;
    picked = list(seen.values())

    fig, axes = plt.subplots(1, len(picked) + 1, figsize=(6.2 * (len(picked) + 1), 6.4), squeeze=False)
    print(f"Payerne robust-comparison prototype — {len(union)} hours, {ZMIN:.0f}-{ZMAX:.0f} m, "
          f"{N_BOOT} block-bootstrap resamples\n")
    summary = {}
    for ci, it in enumerate(picked):
        t = it["ch"]["itype"]; col = TYPE_COLORS.get(t, "#333")
        Yu = it["ungC"][:, band]; Yg = it["gatC"][:, band]
        # (1) SNR-gated median, (2) ungated median
        med_gat = med_relbias_cols(Xg, Yg)
        med_ung = med_relbias_cols(Xu, Yu)
        # (3) Deming with intercept (ungated), noise ratio from lag-1
        lam = (_noise_col(Yu) / np.where(_noise_col(Xu) > 0, _noise_col(Xu), np.nan)) ** 2
        lam = np.where(np.isfinite(lam) & (lam > 0), lam, 1.0)
        rho, alpha, npg = deming_cols(Xu, Yu, lam)
        rel_dem = 100.0 * (rho - 1.0)
        # block bootstrap on rho (resample days)
        boot = np.full((N_BOOT, band.sum()), np.nan)
        day_idx = {d: np.where(days == d)[0] for d in uday}
        for b in range(N_BOOT):
            pick = np.random.default_rng(1234 + b).choice(uday, size=uday.size, replace=True)
            rows = np.concatenate([day_idx[d] for d in pick])
            r_b, _, _ = deming_cols(Xu[rows], Yu[rows], lam)
            boot[b] = 100.0 * (r_b - 1.0)
        lo = np.nanpercentile(boot, 2.5, axis=0); hi = np.nanpercentile(boot, 97.5, axis=0)

        ax = axes[0][ci]
        ax.plot(med_gat, zkm, "-", color="#d62728", lw=1.6, label="(1) median, SNR-gated (current)")
        ax.plot(med_ung, zkm, "-", color="#7f7f7f", lw=1.6, label="(2) median, no SNR filter")
        ax.plot(rel_dem, zkm, "-", color=col, lw=2.2, label="(3) Deming w/ intercept: ρ−1")
        ax.fill_betweenx(zkm, lo, hi, color=col, alpha=0.20, label="(3) 95% block-bootstrap")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("relative bias vs CHM15k  [%]"); ax.set_ylabel("altitude AGL  [km]")
        ax.set_title(f"Payerne {t} ({it['ch']['calib']}) vs CHM15k", fontsize=11, fontweight="bold")
        ax.set_xlim(-40, 60); ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper right")

        # offset panel accumulation
        axes[0][-1].plot(alpha, zkm, "-", color=col, lw=2.0, label=f"{t}: α(z) offset")
        # print band-mean summary
        b12 = (z_agl[band] >= 1000) & (z_agl[band] <= 2000)
        b23 = (z_agl[band] >= 2000) & (z_agl[band] <= 3000)
        b35 = (z_agl[band] >= 2500) & (z_agl[band] <= 3500)
        def bm(a, m): return float(np.nanmedian(a[m]))
        print(f"{t} ({it['ch']['calib']}) vs CHM15k — median rel bias [%] by band:")
        for lab, m in (("1-2 km", b12), ("2-3 km", b23), ("2.5-3.5 km", b35)):
            print(f"   {lab:>10}:  (1)SNR-gated {bm(med_gat,m):+6.1f}   (2)no-filter {bm(med_ung,m):+6.1f}"
                  f"   (3)Deming ρ−1 {bm(rel_dem,m):+6.1f}   | α offset {bm(alpha,m):+.2e} Mm⁻¹sr⁻¹")
        summary[t] = dict(med_gat=med_gat, med_ung=med_ung, rel_dem=rel_dem, alpha=alpha)
        print()

    axoff = axes[0][-1]
    axoff.axvline(0, color="k", lw=0.8)
    axoff.set_xlabel(r"intercept α(z)  [Mm$^{-1}$sr$^{-1}$]  (electronic offset)")
    axoff.set_ylabel("altitude AGL  [km]")
    axoff.set_title("Electronic offset from the Deming intercept", fontsize=11, fontweight="bold")
    axoff.grid(alpha=0.3); axoff.legend(fontsize=9, loc="best")
    fig.suptitle("Offset-robust, censoring-free instrument comparison — prototype (Payerne)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    outp = OUT / "fig_robust_comparison_payerne.png"
    fig.savefig(outp, dpi=170); plt.close(fig)
    print("->", outp)


if __name__ == "__main__":
    main()
