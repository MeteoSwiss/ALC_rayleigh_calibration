"""Slide 16 of the E-PROFILE Annual Meeting 2026 deck, CL31 variant: "where is the CL31?".

Same processing as the paper Appendix-6 cache behind slide16_payerne_dark_wide.py
(validation/paper/run_paper_v09.py, APPENDIX_DARK "payerne_cl61dark"), restricted to the CL31 POST-swap era
(optical block replaced 2026-07-07 ~13 UTC, gain x2.255; the dashboard Kalman needs until 2026-07-18 to
re-converge, which is also slide 15's exclusion) and with the CL31 added:
  * CL61 cloud + dark  : reference (as on slide 16); the CL61 cloud without dark is kept for the
                         "both sides or neither" convention of slide 16;
  * CHM15k             : without dark (operational Rayleigh constant) and with the measured hood dark
                         (constant from the matching dark rerun, <PAY>_A_rayleighdark.csv);
  * CL31 cloud         : without dark and with the measured hood dark B_b_rcs of
                         dark_profiles_payerne_v2.npz. That profile pools the 6 POST-swap hood sessions only
                         (dark_profiles.pool, 2026-07-09 .. 2026-09-29), in rcs_0 units on the CL31 L1 grid.
                         The cloud constant is dark-immune (fixed 100-2400 m integration), so it is reused.
All channels are compared on the hours where every one of them is finite. The paper code is imported,
never edited; this script writes its own cache (OUT below), never the paper caches.

Run from the repo root:  python presentation/slide16_cl31_dark.py [--plot-only]
Output: presentation/figs_annual2026/fig_payerne_dark_wide_cl31.png (+ stats printed as JSON)
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
from matplotlib.lines import Line2D  # noqa: E402

from validation.paper import run_paper_v09 as P  # noqa: E402
from validation.paper import run_paper_validation as RPV  # noqa: E402
from validation.paper import intercompare as IC  # noqa: E402
from validation.paper import figures as FIG  # noqa: E402

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/slide16_cl31")
CACHE = OUT / "cache_payerne_cl31dark.npz"
CL31_DARK = OUT / "baseline_dark_B_rcs.npz"
OUT_PNG = REPO / "presentation" / "figs_annual2026" / "fig_payerne_dark_wide_cl31.png"
START, END = "20260718", "20261008"
PAY = P.PAY

# Channel order of the cache (indices used everywhere below).
I_CL61_D, I_CL61, I_CHM, I_CHM_D, I_CL31, I_CL31_D = range(6)


def _channels():
    return [
        dict(P._ch(PAY, "C", "CL61", "cloud", "CL61 cloud + dark"), bdark=str(P._vnpz("dark", "C"))),
        P._ch(PAY, "C", "CL61", "cloud", "CL61 cloud"),
        P._ch(PAY, "A", "CHM15k", "rayleigh", "CHM15k"),
        dict(P._ch(PAY, "A", "CHM15k", "rayleighdark", "CHM15k + dark"), bdark=str(P._vnpz("dark", "A"))),
        P._ch(PAY, "B", "CL31", "cloud", "CL31 cloud"),
        dict(P._ch(PAY, "B", "CL31", "cloud", "CL31 cloud + dark"), bdark=str(CL31_DARK)),
    ]


def prepare_cl31_dark():
    """CL31 post-swap hood dark in the {rng, b_phys} format channel_beta subtracts from rcs_0.
    Refuses (exits) when the profile cannot be applied cleanly."""
    src = np.load(P.DARK_NPZ_SRC)
    rng = np.asarray(src["B_range"], "f8")
    b = np.asarray(src["B_b_rcs"], "f8")
    nwin = int(src["B_nwin"])
    if nwin < 2:
        sys.exit(f"STOP: CL31 dark pooled on {nwin} window(s) only")
    if np.isnan(b).any():
        sys.exit(f"STOP: CL31 dark has {int(np.isnan(b).sum())} NaN gates")
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(CL31_DARK, rng=rng, b_phys=b)
    return dict(nwin=nwin, nprof=int(src["B_nprof"]), rng_min=float(rng[0]), rng_max=float(rng[-1]),
                n_gates=int(rng.size))


def process():
    """Mirror of run_paper_v09.process_station, writing to OUT. Also checks the CL31 L1 range grid."""
    P._configure()
    items, cache = [], {}
    for ch in _channels():
        k = (ch["wmo"], ch["ident"])
        if k not in cache:
            cache[k] = IC.read_l1(ch["wmo"], ch["ident"], START, END)
        l1 = cache[k]
        if l1 is None:
            sys.exit(f"STOP: no L1 for {ch['label']}")
        if ch["ident"] == "B" and ch.get("bdark"):
            bd = np.load(ch["bdark"])
            zagl = l1["alt"] - l1["station_alt"]
            off = np.abs(np.interp(zagl, bd["rng"], bd["rng"]) - zagl)
            inside = (zagl >= bd["rng"][0]) & (zagl <= bd["rng"][-1])
            print(f"  CL31 L1 grid: {zagl.size} gates {zagl[0]:.1f}..{zagl[-1]:.1f} m; dark grid "
                  f"{bd['rng'].size} gates {bd['rng'][0]:.0f}..{bd['rng'][-1]:.0f} m; "
                  f"{int(inside.sum())} L1 gates inside, max offset {off[inside].max():.2f} m", flush=True)
        res = RPV.channel_beta(l1, ch, P.TARGET_NM, apply_wv=True)
        if res is None:
            sys.exit(f"STOP: no calibration series for {ch['label']}")
        items.append(dict(ch=ch, l1=l1, **res))
        print(f"  {ch['label']}: {l1['time'].size} h, median 1e6/C_L = {res['med_corr']:.4g}", flush=True)
    R = RPV.grid_and_stats(items, 0)
    zagl = R["altGrid"] - R["station"]["altitude"]
    keep = (zagl >= 0) & (zagl <= P.ZTOP_CACHE)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE, z=zagl[keep], time=R["time_sync"].astype("datetime64[s]"),
                        beta=np.stack([b[:, keep] for b in R["beta"]]).astype("f4"),
                        labels=np.array([c["label"] for c in _channels()]))


def _kmed(fn, t0, t1):
    d, v = P._kalman(P.CALIB_DASH / fn)
    s = (d >= t0) & (d <= t1)
    return float(np.median(v[s]))


def figure(dark_info):
    plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "axes.labelsize": 15,
                         "xtick.labelsize": 14, "ytick.labelsize": 14})
    c = np.load(CACHE)
    z, beta = c["z"], c["beta"].astype("f8")
    labels = list(c["labels"])
    print("labels:", labels)
    zmax, band = 4000, P.ZBAND_STATION
    minn = max(30, int(0.02 * beta.shape[1]))
    common = np.all(np.isfinite(beta), axis=0)
    zb = (z >= band[0]) & (z <= band[1])
    n_common_h = int(np.any(common[:, zb], axis=1).sum())

    def m(k):
        return np.where(common, beta[k], np.nan)

    layers = [(100, 500), (500, 1000), (1000, 2000), (2000, 3000)]
    stats = {}
    for name, k in (("CL31_no_dark", I_CL31), ("CL31_dark", I_CL31_D),
                    ("CHM15k_no_dark", I_CHM), ("CHM15k_dark", I_CHM_D), ("CL61_no_dark", I_CL61)):
        stats[name + "_vs_CL61dark"] = P.pair_stats(m(k), m(I_CL61_D), z, band, layers)
    # slide-16 convention ("both sides or neither"): uncorrected channel vs uncorrected CL61
    for name, k in (("CL31_no_dark", I_CL31), ("CHM15k_no_dark", I_CHM)):
        stats[name + "_vs_CL61nodark"] = P.pair_stats(m(k), m(I_CL61), z, band, layers)

    c_chm, c_cl61, c_cl31 = FIG.TYPE_COLORS["CHM15k"], FIG.CL61_COLORS["cloud"], FIG.TYPE_COLORS["CL31"]
    period = P._period(c["time"])
    fig, axs = plt.subplots(1, 4, figsize=(27, 8.2),
                            gridspec_kw=dict(width_ratios=[1, 1, 1.12, 1], wspace=0.32),
                            constrained_layout=False)
    fig.subplots_adjust(left=0.04, right=0.985, top=0.80, bottom=0.13)

    # (a) median profiles: dark removed solid, no dark dashed
    prof = [(m(I_CL61_D), c_cl61, "CL61"), (m(I_CHM_D), c_chm, "CHM15k"), (m(I_CL31_D), c_cl31, "CL31")]
    P.panel_profiles(axs[0], z, prof, band, zmax, minn,
                     f"(a) Median profiles (± IQR)\n{n_common_h} common hours")
    for k, col in ((I_CHM, c_chm), (I_CL31, c_cl31)):
        med, _, _ = P._median_iqr(m(k), minn)
        axs[0].plot(med, z, "--", color=col, lw=1.6, alpha=0.8)
    axs[0].get_legend().remove()

    # (b) relative difference to the CL61 cloud + dark
    def st(name):
        return stats[name + "_vs_CL61dark"]["med"]

    rd = [(P.rel_diff(m(I_CL31), m(I_CL61_D)), "#fdd0a2", "CL31, no dark", st("CL31_no_dark")),
          (P.rel_diff(m(I_CL31_D), m(I_CL61_D)), c_cl31, "CL31, dark removed", st("CL31_dark")),
          (P.rel_diff(m(I_CHM), m(I_CL61_D)), "#f4a582", "CHM15k, no dark", st("CHM15k_no_dark")),
          (P.rel_diff(m(I_CHM_D), m(I_CL61_D)), c_chm, "CHM15k, dark removed", st("CHM15k_dark"))]
    P.panel_reldiff(axs[1], z, rd, band, zmax, minn, "(b) Difference to the CL61 cloud, dark removed")
    axs[1].set_ylabel("")
    axs[1].legend(loc="upper left", framealpha=0.92)

    # (c) bin scatter of the CL31 (dark removed) vs the CL61; CHM15k binned median overlaid
    s31 = stats["CL31_dark_vs_CL61dark"]
    x31 = m(I_CL61_D)[:, zb].ravel()
    P.panel_binscatter(axs[2], x31, m(I_CL31_D)[:, zb].ravel(), c_cl31,
                       "CL61 cloud, dark removed  " + P.BLAB, "Instrument, dark removed  " + P.BLAB,
                       f"(c) CL31 vs CL61, dark removed, {int(band[0])}–{int(band[1])} m\n"
                       f"median {s31['med']:+.1f} %, log r = {s31['r_log']:.2f}")
    y = m(I_CHM_D)[:, zb].ravel()
    ok = np.isfinite(x31) & np.isfinite(y) & (x31 > 0) & (y > 0)
    edges = np.logspace(-2.3, 1.0, 23)
    idx = np.digitize(x31[ok], edges)
    xc, ym = [], []
    for k in range(1, edges.size):
        s = idx == k
        if s.sum() >= 50:
            xc.append(np.sqrt(edges[k - 1] * edges[k]))
            ym.append(np.median(y[ok][s]))
    axs[2].plot(xc, ym, "s-", color=c_chm, mec="k", mew=0.6, ms=5, lw=1.4, label="CHM15k binned median")
    for ln in axs[2].containers:
        ln.set_label("CL31 binned median (IQR)")
    axs[2].legend(loc="upper left")

    # (d) the subtracted darks in % of the molecular signal (same logic as slide 16 panel d)
    ax = axs[3]
    t0 = np.datetime64(f"{START[:4]}-{START[4:6]}-{START[6:]}")
    t1 = np.datetime64(f"{END[:4]}-{END[4:6]}-{END[6:]}")
    cst = {"C": _kmed(f"{PAY}_C_cloud.csv", t0, t1), "A": _kmed(f"{PAY}_A_rayleighdark.csv", t0, t1),
           "B": _kmed(f"{PAY}_B_cloud.csv", t0, t1)}
    src = np.load(P.DARK_NPZ_SRC)
    molp, dpct = {}, {}
    for lam, ident in ((910.55, "C"), (1064.47, "A"), (910.0, "B")):
        molp[ident] = P.molecular_iqr(c["time"], z, lam)[0]
    lo, hi = -5.0, 5.0
    dark_files = {"C": P._vnpz("dark", "C"), "A": P._vnpz("dark", "A"), "B": CL31_DARK}
    for ident, col, nm in (("C", c_cl61, "CL61 (910 nm)"), ("A", c_chm, "CHM15k (1064 nm)"),
                           ("B", c_cl31, "CL31 (910 nm)")):
        dz = np.load(dark_files[ident])
        r = np.asarray(dz["rng"], "f8")
        bd = np.asarray(dz["b_phys"], "f8") / cst[ident] * 1e6
        se = np.nan_to_num(np.asarray(src[f"{ident}_sem"], "f8")) / cst[ident] * 1e6
        mi = np.interp(r, z, molp[ident], left=np.nan, right=np.nan)
        pct, pse = 100.0 * bd / mi, 100.0 * se / mi
        ax.fill_betweenx(r, pct - pse, pct + pse, color=col, alpha=0.25, lw=0)
        ax.plot(pct, r, color=col, lw=2.0, label=nm)
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
    ax.set_title("(d) Measured hood dark (± s.e.)\nin % of the molecular signal")

    for a in axs:
        leg = a.get_legend()
        if leg is not None:
            plt.setp(leg.get_texts(), fontsize=14)
    # one shared key for colours and line styles, under the title
    key = [Line2D([], [], color=c_cl61, lw=2.5, label="CL61 cloud, dark removed (reference)"),
           Line2D([], [], color=c_chm, lw=2.5, label="CHM15k Rayleigh"),
           Line2D([], [], color=c_cl31, lw=2.5, label="CL31 cloud (post-swap dark, 6 hood sessions)"),
           Line2D([], [], color="0.3", lw=2, ls="-", label="dark removed"),
           Line2D([], [], color="0.3", lw=2, ls="--", label="no dark")]
    fig.legend(handles=key, loc="upper center", ncol=5, fontsize=15, frameon=False,
               bbox_to_anchor=(0.5, 0.925))
    fig.suptitle(f"Payerne, CL31 post-swap era {period}: measured hood dark removed (CL61, CHM15k, CL31)",
                 fontsize=18, fontweight="bold", y=0.985)
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=160)

    def slim(s):
        return {"band_500-3000": round(s["med"], 1), "n_hours": s["n_profiles"], "r_log": round(s["r_log"], 2)} | \
            {k: round(v, 1) for k, v in s["layers"].items()}

    print("wrote", OUT_PNG)
    print("STATS", json.dumps(dict(period=period, n_common_hours=n_common_h, cl31_dark=dark_info,
                                   constants_median=cst,
                                   rel_diff={k: slim(v) for k, v in stats.items()},
                                   dark_pct_of_molecular=dpct), indent=1))


def main():
    info = prepare_cl31_dark()
    print("CL31 dark:", info, flush=True)
    if "--plot-only" not in sys.argv or not CACHE.exists():
        process()
    figure(info)


if __name__ == "__main__":
    main()
