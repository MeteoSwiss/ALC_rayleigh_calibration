"""Slide 16 of the E-PROFILE Annual Meeting 2026 deck: Payerne with the measured hood dark removed.

Companion of slide 15 (presentation/slide15_payerne_wide.py): same one-row layout and fonts, built from
the paper's Appendix-6 cache (validation/paper/run_paper_v09.py, APPENDIX_DARK "payerne_cl61dark"):
  * CHM15k + dark : measured hood dark v2 subtracted from rcs_0, constant from the matching dark rerun
                    (diag_v22_dark_v2, improve_alc Kalman -> calib_improve_alc/<PAY>_A_rayleighdark.csv);
  * CL61 cloud + dark : same dark subtraction, operational cloud constant (cloud calibration is
                    dark-immune, so it is reused).
Reference = CL61 cloud (+ dark), as on slide 15; the "no dark" pair is drawn for contrast on the same
matched hours. The CL31 is not in this cache (the paper's dark processing covers the CHM15k and CL61).
Panel (d) shows the subtracted dark divided by the constant of the corrected channel, in % of the CAMS
molecular attenuated backscatter at the instrument's wavelength (x-axis scaled on 300-4000 m).

Run from the repo root:  python presentation/slide16_payerne_dark_wide.py
Output: presentation/figs_annual2026/fig_payerne_dark_wide.png (+ stats printed as JSON)
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from validation.paper import run_paper_v09 as P  # noqa: E402
from validation.paper import figures as FIG  # noqa: E402

OUT_PNG = REPO / "presentation" / "figs_annual2026" / "fig_payerne_dark_wide.png"
KEY = "payerne_cl61dark"


def main():
    plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "axes.labelsize": 15,
                         "xtick.labelsize": 14, "ytick.labelsize": 14})
    st = P.APPENDIX_DARK[KEY]
    path = P.OUT / f"cache_{KEY}.npz"
    if not path.exists():
        sys.exit(f"MISSING CACHE {path}")
    c = np.load(path)
    z, beta = c["z"], c["beta"].astype("f8")
    labels = list(c["labels"])
    n = 1 + len(P.VORDER)
    assert P.VORDER[-1] == "dark", P.VORDER
    i_chm, i_cld = 0, n                       # no correction
    i_chm_d, i_cld_d = n - 1, 2 * n - 1       # + measured dark
    print("labels:", labels, "->", labels[i_chm], "|", labels[i_chm_d], "|", labels[i_cld], "|", labels[i_cld_d])
    assert "dark" in labels[i_chm_d] and "dark" in labels[i_cld_d] and "cloud" in labels[i_cld]

    use = [i_chm, i_chm_d, i_cld, i_cld_d]
    zmax, band = 4000, P.ZBAND_STATION
    minn = max(30, int(0.02 * beta.shape[1]))
    common = np.all(np.isfinite(beta[use]), axis=0)
    zb = (z >= band[0]) & (z <= band[1])
    n_common_h = int(np.any(common[:, zb], axis=1).sum())

    def m(k):
        return np.where(common, beta[k], np.nan)

    c_chm, c_cl61 = FIG.TYPE_COLORS["CHM15k"], FIG.CL61_COLORS["cloud"]
    layers = [(100, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 5000)]
    s_dark = P.pair_stats(m(i_chm_d), m(i_cld_d), z, band, layers)
    s_none = P.pair_stats(m(i_chm), m(i_cld), z, band, layers)

    fig, axs = plt.subplots(1, 4, figsize=(27, 8.2),
                            gridspec_kw=dict(width_ratios=[1, 1, 1.12, 1], wspace=0.32),
                            constrained_layout=False)
    fig.subplots_adjust(left=0.04, right=0.985, top=0.86, bottom=0.13)

    # (a) median profiles, dark removed (solid) and without correction (dashed)
    prof = [(m(i_cld_d), c_cl61, "CL61 cloud, dark removed (ref.)"),
            (m(i_chm_d), c_chm, "CHM15k, dark removed")]
    P.panel_profiles(axs[0], z, prof, band, zmax, minn,
                     f"(a) Median profiles (± IQR)\n{n_common_h} common hours")
    for k, col, lab in ((i_cld, c_cl61, "CL61 cloud, no correction"), (i_chm, c_chm, "CHM15k, no correction")):
        med, _, _ = P._median_iqr(m(k), minn)
        axs[0].plot(med, z, "--", color=col, lw=1.5, alpha=0.75, label=lab)
    axs[0].legend(loc="upper right")

    # (b) CHM15k relative difference to the CL61, without and with the dark removed on both sides
    rd = [(P.rel_diff(m(i_chm), m(i_cld)), "#f4a582", "no correction", s_none["med"]),
          (P.rel_diff(m(i_chm_d), m(i_cld_d)), c_chm, "dark removed", s_dark["med"])]
    P.panel_reldiff(axs[1], z, rd, band, zmax, minn, "(b) CHM15k difference to the CL61 cloud")
    axs[1].set_ylabel("")

    # (c) bin scatter, dark removed
    P.panel_binscatter(axs[2], m(i_cld_d)[:, zb].ravel(), m(i_chm_d)[:, zb].ravel(), c_chm,
                       "CL61 cloud, dark removed  " + P.BLAB, "CHM15k, dark removed  " + P.BLAB,
                       f"(c) CHM15k vs CL61, dark removed, {int(band[0])}–{int(band[1])} m\n"
                       f"median {s_dark['med']:+.1f} %, log r = {s_dark['r_log']:.2f}")

    # (d) the subtracted darks in % of the molecular signal
    ax = axs[3]
    t0 = np.datetime64(f"{st['start'][:4]}-{st['start'][4:6]}-{st['start'][6:]}")
    t1 = np.datetime64(f"{st['end'][:4]}-{st['end'][4:6]}-{st['end'][6:]}")

    def kmed(fn):
        d, v = P._kalman(P.CALIB_DASH / fn)
        s = (d >= t0) & (d <= t1)
        return float(np.median(v[s]))

    cst = {"C": kmed(f"{P.PAY}_C_cloud.csv")}
    try:
        cst["A"] = kmed(f"{P.PAY}_A_rayleighdark.csv")
        a_src = "rayleighdark"
    except Exception as e:  # fall back to the operational constant (paper appendix choice)
        print("A_rayleighdark not usable, operational constant used:", e)
        cst["A"], a_src = kmed(f"{P.PAY}_A_rayleigh.csv"), "rayleigh"
    src = np.load(P.VARIANTS["dark"]["npz_src"])
    # dark / C_L expressed in % of the CAMS molecular attenuated backscatter at the instrument's own
    # wavelength (CL61 910.55 nm, CHM15k 1064 nm), the metric used on the "Dark signal" slide.
    molp, mol, dpct = {}, {}, {}
    for lam, ident, nm in ((910.55, "C", "910 nm"), (1064.47, "A", "1064 nm")):
        mm, _, _ = P.molecular_iqr(c["time"], z, lam)
        molp[ident] = mm
        mol[nm] = [round(float(np.interp(h, z, mm)), 4) for h in (500, 1000, 2000, 3000)]
    lo, hi = -5.0, 5.0
    for ident, col, nm in (("C", c_cl61, "CL61 (910 nm)"), ("A", c_chm, "CHM15k (1064 nm)")):
        dz = np.load(P._vnpz("dark", ident))
        r = np.asarray(dz["rng"], "f8")
        bd = np.asarray(dz["b_phys"], "f8") / cst[ident] * 1e6
        se = np.nan_to_num(np.asarray(src[f"{ident}_sem"], "f8")) / cst[ident] * 1e6
        mi = np.interp(r, z, molp[ident], left=np.nan, right=np.nan)
        pct, pse = 100.0 * bd / mi, 100.0 * se / mi
        ax.fill_betweenx(r, pct - pse, pct + pse, color=col, alpha=0.25, lw=0)
        ax.plot(pct, r, color=col, lw=2.0, label=f"{nm} (± s.e.)")
        s = (r >= 300) & (r <= zmax) & np.isfinite(pct)
        lo, hi = min(lo, np.nanmin(pct[s] - pse[s])), max(hi, np.nanmax(pct[s] + pse[s]))
        dpct[nm] = {f"{a}-{b}": round(float(np.nanmedian(pct[(r >= a) & (r <= b)])), 1)
                    for a, b in ((300, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 4000))}
    pad = 0.12 * (hi - lo)
    ax.axvline(0, color="k", lw=0.9)
    P._band(ax, band)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(0, zmax)
    ax.set_xlabel(r"Dark / C$_L$  [% of molecular $\beta_{att}$]")
    ax.set_ylabel("Altitude a.g.l. [m]")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")
    ax.set_title("(d) Measured hood dark\nin % of the molecular signal")
    mol["dark_pct_of_molecular"] = dpct

    for a in axs:
        leg = a.get_legend()
        if leg is not None:
            plt.setp(leg.get_texts(), fontsize=13)
    fig.suptitle(f"Payerne: {P._period(c['time'])}, measured hood dark removed (CHM15k and CL61)",
                 fontsize=18, fontweight="bold", y=0.985)
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=160)

    def slim(s):
        return {k: (round(v, 2) if isinstance(v, float) else v) for k, v in s.items() if k != "layers"} | \
            {"layers": {k: round(v, 1) for k, v in s["layers"].items()}}

    print("wrote", OUT_PNG)
    print("STATS", json.dumps(dict(period=P._period(c["time"]), n_common_hours=n_common_h,
                                   dark_removed=slim(s_dark), no_correction=slim(s_none),
                                   C_cloud_CL61=cst["C"], C_CHM15k=cst["A"], C_CHM15k_source=a_src,
                                   molecular=mol, dark_npz_keys=sorted(src.files))))


if __name__ == "__main__":
    main()
