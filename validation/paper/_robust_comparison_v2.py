"""
_robust_comparison_v2.py — PROTOTYPE (Payerne, 1 month = June 2026; standalone).

The robust answer to "compare instruments objectively" that a naive per-gate line fit failed to give.
The calibration ratio rho = c_test/c_ref is ONE number (a single lidar constant), not a per-gate
quantity, so we estimate it POOLED across all (gate, time) cells with a PER-GATE OFFSET b(z) removed
and an inverse-noise-variance WEIGHT so the strong-signal cells dominate and the noise floor cannot
fabricate a slope:

    beta_test(z,t) = rho * beta_ref(z,t) + b(z) + noise           (fixed-effects, errors-in-variables)

  1. per-gate noise sigma(z) from the lag-1 hourly difference; weight w(z) = 1/(sig_test^2+sig_ref^2).
  2. de-mean each gate over time (removes b(z)); pool the residuals; weighted Deming slope -> ONE rho.
  3. b(z) = <beta_test>_z - rho*<beta_ref>_z  (the electronic-offset profile, well determined by averaging).
  4. block (per-day) bootstrap CI on rho.
Reported next to: the per-gate Deming rho(z) (a DIAGNOSTIC — honestly wide where SNR~1), and the
reference SNR(z) profile marking the altitude above which the comparison carries no calibration info.

Usage:  python -m validation.paper._robust_comparison_v2
Output: fig_robust_comparison_v2_payerne.png + printed summary.
"""
from __future__ import annotations
import os, sys, warnings
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
from validation.paper._robust_comparison_test import _read_gated_ungated, deming_cols, _noise_col, med_relbias_cols

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
START, END = "20260601", "20260630"        # 1 month for speed
ZMIN, ZMAX = 500.0, 3500.0
N_BOOT = 200
TYPE_COLORS = {"CL31": "#ff7f0e", "CL51": "#9467bd", "CL61": "#1f77b4"}


def build():
    sc = RPV.SITE["payerne"]; target = sc["target"]
    items, cache = [], {}
    for ch in RPV.BENCHMARK["payerne"]["channels"]:
        key = (ch["wmo"], ch["ident"])
        if key not in cache:
            cache[key] = _read_gated_ungated(ch, START, END)
        l1g, l1u = cache[key]
        if l1g is None or l1u is None:
            continue
        rg = RPV.channel_beta(l1g, ch, target); ru = RPV.channel_beta(l1u, ch, target)
        if rg is None or ru is None:
            continue
        items.append(dict(ch=ch, time=l1u["time"], alt=l1u["alt"], salt=l1u["station_alt"],
                          betag=rg["beta_scr"], betau=ru["beta_scr"]))
    union = np.unique(np.concatenate([it["time"] for it in items]))
    for it in items:
        idx = {t: i for i, t in enumerate(it["time"])}
        pos = np.array([idx.get(t, -1) for t in union]); ok = pos >= 0
        for tag, src in (("gatU", "betag"), ("ungU", "betau")):
            M = np.full((union.size, it[src].shape[1]), np.nan); M[ok] = it[src][pos[ok]]; it[tag] = M
    altGrid = IC.build_common_grid([it["alt"] for it in items])
    for it in items:
        it["gatC"] = IC.regrid(it["gatU"], it["alt"], altGrid)
        it["ungC"] = IC.regrid(it["ungU"], it["alt"], altGrid)
    return items, altGrid, items[0]["salt"], union, sc["ref"]


def _demean_pool(X, Y, w):
    """Per-gate de-mean over time (removes the gate offset b(z)), return pooled weighted residual
    sums Sxx,Syy,Sxy and the per-gate means. w is a per-gate (n_gate,) weight."""
    M = np.isfinite(X) & np.isfinite(Y)
    Xm = np.where(M, X, np.nan); Ym = np.where(M, Y, np.nan)
    xbar = np.nanmean(Xm, 0); ybar = np.nanmean(Ym, 0)
    dx = np.where(M, Xm - xbar[None, :], 0.0); dy = np.where(M, Ym - ybar[None, :], 0.0)
    W = np.where(M, w[None, :], 0.0)
    Sxx = np.nansum(W * dx * dx); Syy = np.nansum(W * dy * dy); Sxy = np.nansum(W * dx * dy)
    return Sxx, Syy, Sxy, xbar, ybar


def pooled_rho(X, Y, w, lam):
    Sxx, Syy, Sxy, xbar, ybar = _demean_pool(X, Y, w)
    if not np.isfinite(Sxy) or abs(Sxy) < 1e-30:
        return np.nan, xbar, ybar
    rho = (Syy - lam * Sxx + np.sqrt((Syy - lam * Sxx) ** 2 + 4 * lam * Sxy ** 2)) / (2 * Sxy)
    return rho, xbar, ybar


def main():
    warnings.filterwarnings("ignore")
    items, altGrid, salt, union, iref = build()
    z_agl = altGrid - salt
    band = (z_agl >= ZMIN) & (z_agl <= ZMAX)
    zkm = z_agl[band] / 1000.0
    ref = items[iref]
    Xu = ref["ungC"][:, band]; Xg = ref["gatC"][:, band]
    sx = _noise_col(Xu); snr_ref = np.nanmean(np.where(Xu > 0, Xu, np.nan), 0) / np.where(sx > 0, sx, np.nan)
    days = union.astype("datetime64[D]"); uday = np.unique(days)
    day_idx = {d: np.where(days == d)[0] for d in uday}
    tests = {}
    for k, it in enumerate(items):
        if k == iref or it["ch"]["itype"] not in ("CL31", "CL51", "CL61"):
            continue
        t = it["ch"]["itype"]
        if it["ch"]["calib"] == "cloud" or t not in tests:
            tests[t] = it
    picked = list(tests.values())

    print(f"Payerne robust v2 — {START}..{END}, {len(union)} h, {ZMIN:.0f}-{ZMAX:.0f} m, {N_BOOT} bootstraps\n")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.6))
    for it in picked:
        t = it["ch"]["itype"]; col = TYPE_COLORS.get(t, "#333")
        Yu = it["ungC"][:, band]
        sy = _noise_col(Yu)
        w = 1.0 / np.where((sx ** 2 + sy ** 2) > 0, sx ** 2 + sy ** 2, np.nan)       # inverse total noise var
        w = np.where(np.isfinite(w), w, 0.0)
        lam = float(np.nansum(w * sy ** 2) / max(np.nansum(w * sx ** 2), 1e-30))
        # (A) pooled robust rho + offset b(z)
        rho, xbar, ybar = pooled_rho(Xu, Yu, w, lam)
        bz = ybar - rho * xbar
        # bootstrap rho over days
        rb = []
        for b in range(N_BOOT):
            pick = np.random.default_rng(99 + b).choice(uday, uday.size, replace=True)
            rows = np.concatenate([day_idx[d] for d in pick])
            r, _, _ = pooled_rho(Xu[rows], Yu[rows], w, lam)
            rb.append(r)
        rlo, rhi = np.nanpercentile(rb, [2.5, 97.5])
        # (B) per-gate diagnostic rho(z) (equal weight, the "not robust" one) for contrast
        rg, ag, _ = deming_cols(Xu, Yu, (sy / np.where(sx > 0, sx, np.nan)) ** 2)
        # informative altitude = where reference SNR >= 1
        zinfo = zkm[snr_ref >= 1.0]
        z_top = float(zinfo.max()) if zinfo.size else np.nan

        print(f"{t} ({it['ch']['calib']}) vs CHM15k:")
        print(f"   POOLED robust calibration ratio  rho = {rho:.3f}  ->  bias {100*(rho-1):+.1f} %"
              f"   [95% CI {100*(rlo-1):+.1f} .. {100*(rhi-1):+.1f}]")
        print(f"   electronic offset b(z):  {np.nanmedian(bz[z_agl[band]<=1500]):+.2e} (<1.5km)  "
              f"{np.nanmedian(bz[(z_agl[band]>=2000)&(z_agl[band]<=2500)]):+.2e} (2-2.5km)  "
              f"{np.nanmedian(bz[z_agl[band]>=3000]):+.2e} (>3km) Mm⁻¹sr⁻¹")
        print(f"   reference SNR>=1 up to ~{z_top:.1f} km  -> calibration info only below that\n")

        # per-gate diagnostic + pooled line
        axes[0].plot(100 * (rg - 1), zkm, "-", color=col, lw=1.0, alpha=0.5,
                     label=f"{t}: per-gate ρ−1 (diagnostic, bruité)")
        axes[0].axvline(100 * (rho - 1), color=col, lw=2.2, ls="-",
                        label=f"{t}: ρ−1 mutualisé = {100*(rho-1):+.1f}%")
        axes[0].axvspan(100 * (rlo - 1), 100 * (rhi - 1), color=col, alpha=0.12)
        axes[1].plot(bz, zkm, "-", color=col, lw=2.0, label=f"{t}: offset b(z)")

    axes[0].axvline(0, color="k", lw=0.8); axes[0].set_xlim(-40, 40)
    axes[0].set_xlabel("bias vs CHM15k  [%]"); axes[0].set_ylabel("altitude AGL  [km]")
    axes[0].set_title("Calibration ratio: mutualisé robuste (barre) vs per-gate (bruité)", fontsize=10.5, fontweight="bold")
    axes[0].grid(alpha=0.3); axes[0].legend(fontsize=8, loc="upper right")
    axes[1].axvline(0, color="k", lw=0.8)
    axes[1].set_xlabel(r"offset b(z)  [Mm$^{-1}$sr$^{-1}$]"); axes[1].set_ylabel("altitude AGL  [km]")
    axes[1].set_title("Profil d'offset électronique (intercept par gate)", fontsize=10.5, fontweight="bold")
    axes[1].grid(alpha=0.3); axes[1].legend(fontsize=9, loc="best")
    fig.suptitle("Comparaison objective robuste — ρ mutualisé pondéré-SNR + offset b(z) (Payerne, juin 2026)",
                 fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    outp = OUT / "fig_robust_comparison_v2_payerne.png"
    fig.savefig(outp, dpi=170); plt.close(fig)
    print("->", outp)


if __name__ == "__main__":
    main()
