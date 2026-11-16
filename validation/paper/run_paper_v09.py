"""
run_paper_v09.py — Section 6 (Validation) figures and statistics for paper version 0.9.

Version 0.9 = version 0.8 with ONE change: the daily best-estimate constants are the dashboard nightly
values (v2.2 Rayleigh / O'Connor cloud, fetched 2026-10-09) smoothed with the OPERATIONAL
improve_alc_calib Kalman filter (validation/paper/kalman_improve_alc.py; calib_improve_alc/ built from
paper_v08/calib_dashboard/) instead of the dashboard display filter. Everything else is identical.

What changes with respect to run_paper_validation.py (v0.7 figures):
  * Calibration = ONLY the constants published on the operational dashboard
    (https://alc-calib.ch-meteoswiss-emermet.f.ewcloud.host/): the daily Kalman C_L of each stream,
    "v2.2" for Rayleigh, "O'Connor" for the liquid-cloud method, fetched by fetch_dashboard_calib.py.
    No local reruns, no measured-dark variants (dark subtraction is not operational).
  * CL61 water-vapour correction with the operational manufacturer spectrum (910.55 nm, FWHM 0.188 nm),
    i.e. the spectrum the dashboard constants were computed with.
  * One figure family for every comparison:
        (a) median of the matched profiles (+/- IQR),
        (b) bin scatter: 2-D histogram of the compared vs reference beta_att with the binned median,
        (c) median of the relative difference 100*(X - ref)/ref per altitude (+/- IQR).
    EARLINET sites (Leipzig, Magurele) and Palaiseau (CHM15k vs CL31) use one row (a, b, c); the
    multi-instrument sites (Payerne with the CL61 as reference, Amsterdam with four CHM15k) show (a)
    and the relative-difference panel on top and one bin-scatter panel per compared instrument below.

Processing per channel is unchanged (run_paper_validation.channel_beta): L1 rcs_0 -> CHM15k overlap
temperature correction (Hervo et al., 2016) -> beta_att = rcs_0 / C_L -> 910 nm WV correction ->
910 -> 1064 nm conversion (CAMS molecular + Angstrom 1 on the aerosol part) -> cloud/fog screening ->
hourly median -> common altitude grid. Statistics over 500-3000 m a.g.l. (stations) and 500-5000 m
a.g.l. (EARLINET), headline metric = median relative difference.

The expensive part (L1 reading) is cached per site in OUT/cache_<site>.npz; re-plot with --plot-only.

Usage:  python -m validation.paper.run_paper_v09 [site ...] [--plot-only]
        (paper run of 2026-10-10: ALC_VAL_SITE_WORKERS=8 python -m validation.paper.run_paper_v09,
         log in OUT/run_v09.log; needs calib_improve_alc/ from build_calib_improve_alc.py and the Payerne
         dark reruns diag_v22_dark_v2 from rayleigh_availability/run_dark_v2.ps1)
        sites: payerne amsterdam palaiseau leipzig magurele magurele_a   (default: all)
"""
from __future__ import annotations
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
# Operational CL61 emission spectrum for every calibration routine imported below (single source of
# truth water_vapor.laser_spectrum_for); intercompare.WV_PARAMS carries the same values.
os.environ.setdefault("ALC_WV_SPECTRUM", '{"CL61": [910.55, 0.188]}')

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from validation.paper import intercompare as IC
from validation.paper import earlinet as EA
from validation.paper import run_paper_validation as RPV
from validation.paper import figures as FIG
from validation.paper.kalman_improve_alc import kalman_improve_alc
import pandas as pd

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v09")
CALIB_DASH = OUT / "calib_improve_alc"   # dashboard nightly values + improve_alc Kalman (SOURCE.txt)
TARGET_NM = 1064.0                 # every channel is expressed at 1064 nm (CHM15k / EARLINET line)
ZBAND_STATION = (500.0, 3000.0)    # statistics band, stations [m a.g.l.]
ZBAND_EARLINET = (500.0, 5000.0)   # statistics band, EARLINET [m a.g.l.]
ZTOP_CACHE = 8000.0                # cached altitude range [m a.g.l.]


def _ch(wmo, ident, itype, calib, label, exclude=()):
    return dict(wmo=wmo, ident=ident, itype=itype, calib=calib, label=label, exclude=list(exclude))


PAY, AMS, SIR = "0-20000-0-06610", "0-20000-0-06240", "0-250-1001-07151"
# Station comparisons. ref = index of the reference channel. Periods = the overlap of the dashboard
# calibration record with the local L1 archive (Payerne CL61 from its first Rayleigh Kalman day).
STATIONS = {
    "payerne": dict(
        name="Payerne", wmo=PAY, start="20260304", end="20261008", ref=0, channels=[
            _ch(PAY, "C", "CL61", "cloud", "CL61 (cloud)"),
            _ch(PAY, "A", "CHM15k", "rayleigh", "CHM15k (Rayleigh)"),
            # CL31 optical block replaced 2026-07-07 ~13 UTC: the dashboard Kalman keeps the
            # pre-swap constant until the first post-swap cloud calibration (2026-07-15/17), so
            # these days carry a ~2x-wrong constant and are excluded (calibration transition).
            _ch(PAY, "B", "CL31", "cloud", "CL31 (cloud)",
                exclude=[("2026-07-07T00:00", "2026-07-18T00:00")]),
            # The CL61 Rayleigh constant is NOT shown here: at Payerne it is contaminated by the
            # CL61 dark noise (electronic baseline) in the molecular fit window. It is analysed,
            # with and without the measured dark correction, in the appendix (APPENDIX_DARK).
        ]),
    "amsterdam": dict(
        name="Amsterdam Schiphol", wmo=AMS, start="20250111", end="20261008", ref=3, channels=[
            _ch(AMS, "A", "CHM15k", "rayleigh", "CHM15k A"),
            _ch(AMS, "B", "CHM15k", "rayleigh", "CHM15k B"),
            _ch(AMS, "C", "CHM15k", "rayleigh", "CHM15k C"),
            _ch(AMS, "D", "CHM15k", "rayleigh", "CHM15k D"),
        ]),
    "palaiseau": dict(
        name="Palaiseau (SIRTA)", wmo=SIR, start="20250103", end="20261008", ref=0, channels=[
            _ch(SIR, "B", "CHM15k", "rayleigh", "CHM15k (Rayleigh)"),
            _ch(SIR, "A", "CL31", "cloud", "CL31 (cloud)"),
        ]),
}
# Appendix 6: Payerne CL61, cloud vs Rayleigh calibration, without and with corrections of the signal
# baseline, the CHM15k being the independent reference. Three corrections (VARIANTS), each subtracted
# from rcs_0 (b(z) in rcs_0 units) AND, for Rayleigh, paired with constants from a calibration rerun
# that applied the same correction (identical to the operational v2.2 otherwise, CL61 with the
# operational 910.55/0.188 spectrum):
#   bg     : operational far-field background bias only -- the clear-night median of [measured -
#            modelled molecular] P at 12-15 km, added back as c*r^2 (no hood needed);
#   dark   : measured hood dark, v2 profile (10 sessions, windows verified in the L1 data);
#   darkbg : both.
# The cloud constant is insensitive to both (integration 100-2400 m, cloud signal >> b), so it is reused.
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
RA = DATA / "rayleigh_availability"
VARIANTS = {
    "bg": dict(label="background", npz_src=RA / "background_only_payerne.npz", run=DATA / "diag_v22_bgonly"),
    "dark": dict(label="dark", npz_src=RA / "dark_profiles_payerne_v2.npz", run=DATA / "diag_v22_dark_v2"),
    "darkbg": dict(label="dark + background", npz_src=RA / "dark_profiles_payerne_v2_bg.npz",
                   run=DATA / "diag_v22_darkbg_v2"),
}
# Variants used in the paper. The far-field background bias was found negligible for the paper
# (operator decision 2026-10-09); "bg"/"darkbg" stay available for diagnostics.
VORDER = ("dark",)
DARK_NPZ_SRC = VARIANTS["dark"]["npz_src"]


def _vnpz(v, ident):
    return OUT / f"baseline_{v}_{ident}_rcs.npz"


def _vrun(v, ident, kind):
    return VARIANTS[v]["run"] / f"{PAY}_{ident}" / f"{PAY}_{ident}_{kind}.csv"


def _appendix_channels():
    """Index layout (n = 1 + len(VORDER)): CHM15k 0..n-1, CL61 cloud n..2n-1, CL61 Rayleigh 2n..3n-1,
    each block ordered (no correction, *VORDER)."""
    ch = [_ch(PAY, "A", "CHM15k", "rayleigh", "CHM15k")]
    ch += [dict(_ch(PAY, "A", "CHM15k", "rayleigh" + v, f"CHM15k + {VARIANTS[v]['label']}"),
                bdark=str(_vnpz(v, "A"))) for v in VORDER]
    ch += [_ch(PAY, "C", "CL61", "cloud", "CL61 cloud")]
    ch += [dict(_ch(PAY, "C", "CL61", "cloud", f"CL61 cloud + {VARIANTS[v]['label']}"),
                bdark=str(_vnpz(v, "C"))) for v in VORDER]
    ch += [_ch(PAY, "C", "CL61", "rayleigh", "CL61 Rayleigh")]
    ch += [dict(_ch(PAY, "C", "CL61", "rayleigh" + v, f"CL61 Rayleigh + {VARIANTS[v]['label']}"),
                bdark=str(_vnpz(v, "C"))) for v in VORDER]
    return ch


APPENDIX_DARK = {
    "payerne_cl61dark": dict(
        name="Payerne CL61 baseline appendix", wmo=PAY, start="20260304", end="20261008", ref=0,
        ztop=15500.0, zplot=15000.0,
        channels=_appendix_channels()),
}
# "Both sides or neither": each CL61 configuration is compared with the CHM15k processed the same way.
_NV = 1 + len(VORDER)
DARK_PAIRS = [(_NV + i, i) for i in range(_NV)] + [(2 * _NV + i, i) for i in range(_NV)]
ALL_STATIONS = {**STATIONS, **APPENDIX_DARK}
# EARLINET comparisons: earlinet.SITES code + display names.
EARLINET = {
    "leipzig": dict(code="ari", name="Leipzig", lidar="PollyXT", instr="CHM15k (Rayleigh)",
                    start="20250101", end="20261008"),
    "magurele": dict(code="ino", name="Bucharest-M\u0103gurele", lidar="RALI",
                     instr="CHM15k B (Rayleigh)", start="20250101", end="20261008"),
    "magurele_a": dict(code="ino_a", name="Bucharest-M\u0103gurele", lidar="RALI",
                       instr="CHM15k A (Rayleigh)", start="20250101", end="20261008"),
}


def _configure():
    """Point the shared library at the dashboard constants (re-run in every worker process, since a
    spawned child re-imports intercompare with its defaults)."""
    IC.CALIB = CALIB_DASH
    IC.WV_PARAMS["CL61"] = (910.55, 0.188)
    warnings.filterwarnings("ignore")


def prepare_dark_inputs():
    """Appendix 6 inputs per correction variant: (1) b(z) in the {rng, b_phys} npz format that
    run_paper_validation.channel_beta subtracts from rcs_0; (2) the Rayleigh Kalman series of the
    matching rerun as <stream>_rayleigh<variant>.csv next to the dashboard series (NOT dashboard values,
    recorded in SOURCE_rayleigh_variants.txt)."""
    import csv
    src = []
    for v in VORDER:
        z = np.load(VARIANTS[v]["npz_src"])
        for ident in ("A", "C"):
            b = np.nan_to_num(np.asarray(z[f"{ident}_b_rcs"], "f8"), nan=0.0)
            np.savez(_vnpz(v, ident), rng=np.asarray(z[f"{ident}_range"], "f8"), b_phys=b)
            rows = ["time,C_daily,C_daily_std,C_kalman,C_kalman_std"]
            # v0.9: the operational improve_alc Kalman on the nightly constants of the local rerun
            dd, vv = _calib_points(_vrun(v, ident, "cal"), "rayleigh")
            kr = kalman_improve_alc([pd.Timestamp(x).to_pydatetime() for x in dd], vv,
                                    "CL61" if ident == "C" else "CHM15k")
            for gd, gs in zip(kr["grid"], kr["state"]):
                rows.append(f"{gd:%Y-%m-%d},,,{gs:.6e},")
            (CALIB_DASH / f"{PAY}_{ident}_rayleigh{v}.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
            src.append(f"{PAY}_{ident}_rayleigh{v}.csv = improve_alc Kalman of the local rerun {_vrun(v, ident, 'cal')}")
    (CALIB_DASH / "SOURCE_rayleigh_variants.txt").write_text(
        "\n".join(src) + "\n(appendix only; not operational/dashboard products)\n", encoding="utf-8")


# --------------------------------------------------------------------------- processing
def _apply_exclusions(res, l1, ch):
    for a, b in ch.get("exclude", []):
        m = (l1["time"] >= np.datetime64(a)) & (l1["time"] < np.datetime64(b))
        res["beta_scr"][m] = np.nan
        res["beta_disp"][m] = np.nan


def process_station(key):
    _configure()
    st = ALL_STATIONS[key]
    items, cache = [], {}
    for ch in st["channels"]:
        k = (ch["wmo"], ch["ident"])
        if k not in cache:
            cache[k] = IC.read_l1(ch["wmo"], ch["ident"], st["start"], st["end"])
        l1 = cache[k]
        if l1 is None:
            raise RuntimeError(f"{key}: no L1 for {ch['label']}")
        res = RPV.channel_beta(l1, ch, TARGET_NM, apply_wv=True)
        if res is None:
            raise RuntimeError(f"{key}: no dashboard calibration for {ch['label']}")
        _apply_exclusions(res, l1, ch)
        items.append(dict(ch=ch, l1=l1, **res))
        print(f"  [{key}] {ch['label']}: {l1['time'].size} h, median 1e6/C_L = {res['med_corr']:.4g}",
              flush=True)
    R = RPV.grid_and_stats(items, st["ref"])
    salt = R["station"]["altitude"]
    zagl = R["altGrid"] - salt
    keep = (zagl >= 0) & (zagl <= st.get("ztop", ZTOP_CACHE))
    np.savez_compressed(
        OUT / f"cache_{key}.npz", z=zagl[keep], time=R["time_sync"].astype("datetime64[s]"),
        beta=np.stack([b[:, keep] for b in R["beta"]]).astype("f4"),
        labels=np.array([c["label"] for c in st["channels"]]),
        itypes=np.array([c["itype"] for c in st["channels"]]),
        calibs=np.array([c["calib"] for c in st["channels"]]), ref=st["ref"])
    return key


def process_earlinet(key):
    _configure()
    e = EARLINET[key]
    s = EA.compare(e["code"], e["start"], e["end"], return_profiles=True)
    if not s or "error" in s:
        raise RuntimeError(f"{key}: {s}")
    np.savez_compressed(OUT / f"cache_{key}.npz", z=s["grid"], E=s["betaE"].astype("f4"),
                        C=s["betaC"].astype("f4"), time=np.array(s["times"]).astype("datetime64[s]"))
    return key


# --------------------------------------------------------------------------- statistics
def rel_diff(x, ref):
    """100*(x-ref)/ref where both finite and ref>0, else NaN."""
    with np.errstate(all="ignore"):
        d = 100.0 * (x - ref) / ref
    d[~(np.isfinite(x) & np.isfinite(ref) & (ref > 0))] = np.nan
    return d


def pair_stats(x, ref, z, band, layers):
    """Median relative difference (+IQR), log-r and N over the band, and per altitude layer."""
    zb = (z >= band[0]) & (z <= band[1])
    d = rel_diff(x[:, zb], ref[:, zb])
    m = np.isfinite(d)
    a, b = x[:, zb][m], ref[:, zb][m]
    pos = (a > 0) & (b > 0)
    out = dict(n=int(m.sum()), n_profiles=int(np.any(m, axis=1).sum()),
               med=float(np.nanmedian(d)) if m.any() else np.nan,
               q25=float(np.nanpercentile(d, 25)) if m.any() else np.nan,
               q75=float(np.nanpercentile(d, 75)) if m.any() else np.nan,
               r_log=float(np.corrcoef(np.log10(a[pos]), np.log10(b[pos]))[0, 1]) if pos.sum() > 2 else np.nan,
               ratio_of_medians=float(np.nanmedian(a) / np.nanmedian(b)) if m.any() else np.nan,
               layers={})
    for lo, hi in layers:
        zl = (z >= lo) & (z < hi)
        dl = rel_diff(x[:, zl], ref[:, zl])
        out["layers"][f"{int(lo)}-{int(hi)}"] = (float(np.nanmedian(dl))
                                                 if np.isfinite(dl).any() else np.nan)
    return out


# --------------------------------------------------------------------------- figures
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec

plt.rcParams.update({"font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10})
BLAB = r"$\beta_{att}$ at 1064 nm [Mm$^{-1}$ sr$^{-1}$]"
LOGBINS = np.logspace(-2.3, 1.0, 67)


def _median_iqr(B, minn):
    with np.errstate(all="ignore"):
        med = np.nanmedian(B, axis=0)
        q1 = np.nanpercentile(B, 25, axis=0)
        q3 = np.nanpercentile(B, 75, axis=0)
    n = np.isfinite(B).sum(axis=0)
    bad = n < minn
    med[bad] = q1[bad] = q3[bad] = np.nan
    return med, q1, q3


def _band(ax, band):
    ax.axhspan(band[0], band[1], color="0.5", alpha=0.08, lw=0, zorder=0)


def panel_profiles(ax, z, series, band, zmax, minn, title):
    """series: list of (B[time x z], colour, label). B already restricted to the matched sample."""
    xmax = 0.0
    for B, col, lab in series:
        med, q1, q3 = _median_iqr(B, minn)
        ax.fill_betweenx(z, q1, q3, color=col, alpha=0.13, lw=0)
        ax.plot(med, z, color=col, lw=1.8, label=lab)
        sel = (z > 150) & (z < zmax)
        if np.isfinite(med[sel]).any():
            xmax = max(xmax, np.nanpercentile(med[sel], 99))   # scale on medians, not a noisy IQR
    _band(ax, band)
    ax.set_xlim(0, max(0.2, 1.6 * xmax))
    ax.set_ylim(0, zmax)
    ax.set_xlabel(BLAB)
    ax.set_ylabel("Altitude a.g.l. [m]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper right")
    ax.set_title(title)


def panel_reldiff(ax, z, series, band, zmax, minn, title, xlim=60):
    """series: list of (D[time x z] relative difference %, colour, label, band median)."""
    for D, col, lab, med_band in series:
        med, q1, q3 = _median_iqr(D, minn)
        ax.fill_betweenx(z, q1, q3, color=col, alpha=0.12, lw=0)
        ax.plot(med, z, color=col, lw=1.8, label=f"{lab}  (median {med_band:+.1f} %)")
    ax.axvline(0, color="k", lw=0.9)
    for v in (-10, 10):
        ax.axvline(v, color="0.4", lw=0.7, ls=":")
    _band(ax, band)
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(0, zmax)
    ax.set_xlabel("Relative difference to the reference [%]")
    ax.set_ylabel("Altitude a.g.l. [m]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper right")
    ax.set_title(title)


def panel_binscatter(ax, x, y, col, xlab, ylab, title, cbar=True):
    """2-D histogram (log-log) of y vs x with the binned median of y (+ IQR bars) and the 1:1 line."""
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x, y = x[m], y[m]
    h = ax.hist2d(x, y, bins=[LOGBINS, LOGBINS], cmap="Blues", norm=LogNorm(vmin=1), rasterized=True)
    if cbar:
        cb = plt.colorbar(h[3], ax=ax, pad=0.02)
        cb.set_label("counts")
    edges = np.logspace(-2.3, 1.0, 23)
    idx = np.digitize(x, edges)
    xc, ym, y1, y3 = [], [], [], []
    for k in range(1, edges.size):
        s = idx == k
        if s.sum() >= 50:
            xc.append(np.sqrt(edges[k - 1] * edges[k]))
            ym.append(np.median(y[s]))
            y1.append(np.percentile(y[s], 25))
            y3.append(np.percentile(y[s], 75))
    xc, ym, y1, y3 = map(np.array, (xc, ym, y1, y3))
    if xc.size:
        ax.errorbar(xc, ym, yerr=[ym - y1, y3 - ym], color=col, mec="k", mew=0.6, marker="o", ms=4.5,
                    lw=1.4, capsize=2, label="binned median (IQR)")
    ax.plot([1e-3, 1e2], [1e-3, 1e2], "k--", lw=1, label="1:1")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(LOGBINS[0], LOGBINS[-1])
    ax.set_ylim(LOGBINS[0], LOGBINS[-1])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.grid(alpha=0.3, which="major")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title(title)


def _period(t):
    # hourly bins are labelled by their start; drop the empty midnight bin that closes the last day
    t = np.asarray(t).astype("datetime64[s]")
    first = t.min().astype("datetime64[D]")
    last = (t.max() - np.timedelta64(1, "s")).astype("datetime64[D]")
    return f"{first} to {last}"


def fig_station(key):
    st = STATIONS[key]
    c = np.load(OUT / f"cache_{key}.npz")
    z, beta, ref = c["z"], c["beta"].astype("f8"), int(c["ref"])
    labels, chans = list(c["labels"]), st["channels"]
    cols = FIG.channel_colors(chans)
    comps = [k for k in range(len(chans)) if k != ref]
    minn = max(30, int(0.02 * beta.shape[1]))
    # (a) matched sample: gate-wise intersection of all channels
    common = np.all(np.isfinite(beta), axis=0)
    n_common_h = int(np.any(common[:, (z >= ZBAND_STATION[0]) & (z <= ZBAND_STATION[1])], axis=1).sum())
    prof = [(np.where(common, beta[k], np.nan), cols[k], labels[k] + (" (reference)" if k == ref else ""))
            for k in range(len(chans))]
    layers = [(100, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 5000)]
    # night (22-02 UTC, dark at every site/season) vs midday (10-14 UTC): separates calibration effects
    # (identical day and night) from background/solar artefacts
    hour = c["time"].astype("datetime64[h]").astype(int) % 24
    sel = {"night": (hour >= 22) | (hour <= 2), "day": (hour >= 10) & (hour <= 14)}
    stats, rd = {}, []
    for k in comps:
        s = pair_stats(beta[k], beta[ref], z, ZBAND_STATION, layers)
        for nm, m in sel.items():
            s[nm] = pair_stats(beta[k][m], beta[ref][m], z, ZBAND_STATION, layers)
        stats[labels[k]] = s
        rd.append((rel_diff(beta[k], beta[ref]), cols[k], labels[k], s["med"]))
    single = len(comps) == 1
    zmax = 4000
    if single:
        fig, axs = plt.subplots(1, 3, figsize=(17, 5.6), gridspec_kw=dict(width_ratios=[1, 1.15, 1]))
        ax_a, ax_s, ax_c = axs
        scat_axes = [ax_s]
    else:
        fig = plt.figure(figsize=(18, 10.5))
        gs = GridSpec(2, 6, figure=fig, height_ratios=[1.05, 1], hspace=0.32, wspace=0.55)
        ax_a = fig.add_subplot(gs[0, 0:3])
        ax_c = fig.add_subplot(gs[0, 3:6])
        scat_axes = [fig.add_subplot(gs[1, 2 * i:2 * i + 2]) for i in range(len(comps))]
    panel_profiles(ax_a, z, prof, ZBAND_STATION, zmax, minn,
                   f"(a) Median of matched profiles (\u00b1 IQR), {n_common_h} common hours")
    letters = "bcdefg"
    zb = (z >= ZBAND_STATION[0]) & (z <= ZBAND_STATION[1])
    for i, (k, ax) in enumerate(zip(comps, scat_axes)):
        s = stats[labels[k]]
        lt = "b" if single else letters[1 + i]
        panel_binscatter(ax, beta[ref][:, zb].ravel(), beta[k][:, zb].ravel(), cols[k],
                         f"{labels[ref]}  " + BLAB, f"{labels[k]}  " + BLAB,
                         f"({lt}) {labels[k]} vs {labels[ref]}, {int(ZBAND_STATION[0])}\u2013"
                         f"{int(ZBAND_STATION[1])} m\nmedian rel. diff. {s['med']:+.1f} %, "
                         f"log r = {s['r_log']:.2f}, N = {s['n']}", cbar=True)
    panel_reldiff(ax_c, z, rd, ZBAND_STATION, zmax, minn,
                  f"({'c' if single else 'b'}) Median relative difference to {labels[ref]} (\u00b1 IQR)")
    fig.suptitle(f"{st['name']} ({st['wmo']}): {_period(c['time'])}, operational dashboard calibration",
                 fontsize=12.5, fontweight="bold", y=0.995)
    fig.savefig(OUT / f"fig_v09_{key}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return dict(site=st["name"], reference=labels[ref], period=_period(c["time"]),
                n_common_hours=n_common_h, pairs=stats)


def fig_earlinet(key):
    e = EARLINET[key]
    c = np.load(OUT / f"cache_{key}.npz")
    z, E, C = c["z"], c["E"].astype("f8"), c["C"].astype("f8")
    layers = [(500, 1000), (1000, 2000), (2000, 3000), (3000, 5000)]
    s = pair_stats(C, E, z, ZBAND_EARLINET, layers)
    s["matched_profiles"] = int(E.shape[0])
    minn = max(20, int(0.05 * E.shape[0]))
    lowest_e = float(z[np.argmax(np.isfinite(E).sum(axis=0) >= minn)])
    zmax = 6000
    colC = FIG.TYPE_COLORS["CHM15k"]
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.6), gridspec_kw=dict(width_ratios=[1, 1.15, 1]))
    panel_profiles(axs[0], z, [(E, "k", f"EARLINET {e['lidar']} ({e['name']})"), (C, colC, e["instr"])],
                   ZBAND_EARLINET, zmax, minn, f"(a) Median of matched profiles (\u00b1 IQR), "
                   f"{E.shape[0]} matched profiles")
    zb = (z >= ZBAND_EARLINET[0]) & (z <= ZBAND_EARLINET[1])
    panel_binscatter(axs[1], E[:, zb].ravel(), C[:, zb].ravel(), colC, f"EARLINET {e['lidar']}  " + BLAB,
                     f"{e['instr']}  " + BLAB,
                     f"(b) {int(ZBAND_EARLINET[0])}\u2013{int(ZBAND_EARLINET[1])} m: median rel. diff. "
                     f"{s['med']:+.1f} %, log r = {s['r_log']:.2f}\nN = {s['n']} matched range bins")
    panel_reldiff(axs[2], z, [(rel_diff(C, E), colC, e["instr"], s["med"])], ZBAND_EARLINET, zmax, minn,
                  "(c) Median relative difference to EARLINET (\u00b1 IQR)")
    fig.suptitle(f"{e['name']}: EARLINET {e['lidar']} vs {e['instr']}, {_period(c['time'])}, "
                 f"operational dashboard calibration", fontsize=12.5, fontweight="bold", y=1.0)
    fig.savefig(OUT / f"fig_v09_{key}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return dict(site=e["name"], lidar=e["lidar"], instrument=e["instr"], period=_period(c["time"]),
                lowest_earlinet_m=lowest_e, stats=s)


# --------------------------------------------------------------------------- main
def _calib_points(path, method=None):
    """(dates, C_L) of successful nightly constants from a calib CSV (C_daily) or a *_cal.csv."""
    import csv
    d, v = [], []
    for r in csv.DictReader(open(path, encoding="utf-8")):
        if "cal_value" in r:
            if r["method"] != method or r["flag"] not in ("1", "0.5"):
                continue
            s, x = r["date"], r["cal_value"]
            s = f"{s[:4]}-{s[4:6]}-{s[6:8]}"
        else:
            s, x = r["time"], r["C_daily"]
        if x not in ("", "nan") and float(x) > 0:
            d.append(np.datetime64(s))
            v.append(float(x))
    return np.array(d), np.array(v)


def _kalman(path):
    """(dates, C_kalman) from a calib CSV."""
    import csv
    rows = [r for r in csv.DictReader(open(path, encoding="utf-8")) if r["C_kalman"] not in ("", "nan")]
    return (np.array([np.datetime64(r["time"]) for r in rows]),
            np.array([float(r["C_kalman"]) for r in rows]))


def molecular_iqr(times, z_agl, lam_nm, lat=46.8137, lon=6.9425, alt=491.0):
    """Median and IQR over the period of the molecular attenuated backscatter [Mm-1 sr-1] at lam_nm,
    from the CAMS T/p profiles of each CAMS period (monthly, else daily), weighted by the hours."""
    l2 = dict(time=np.asarray(times).astype("datetime64[ns]"), lat=lat, lon=lon)
    profs, w = [], []
    for sel, cf, tstart, tend, _p in IC.cams_periods(l2["time"]):
        if cf is None:
            continue
        b = IC._mol_att_cams(cf, lat, lon, tstart, tend, z_agl + alt, alt, lam_nm)
        if b is not None:
            profs.append(b)
            w.append(int(sel.sum()))
    P = np.repeat(np.array(profs), w, axis=0)
    return np.nanmedian(P, axis=0), np.nanpercentile(P, 25, axis=0), np.nanpercentile(P, 75, axis=0)


def fig_appendix_dark(key="payerne_cl61dark"):
    """Appendix 6 figure: (a) CL61 constants (cloud; Rayleigh without / with the correction), (b)
    Rayleigh-to-cloud ratio per configuration, (c) the measured darks in beta units with the molecular
    signal (median and IQR) for scale, (d) median matched profiles with the molecular reference, (e)
    median relative difference to the CHM15k processed the same way."""
    st = APPENDIX_DARK[key]
    c = np.load(OUT / f"cache_{key}.npz")
    z, beta = c["z"], c["beta"].astype("f8")
    labels = list(c["labels"])
    zplot = st.get("zplot", 6000.0)
    t0 = np.datetime64(f"{st['start'][:4]}-{st['start'][4:6]}-{st['start'][6:]}")
    t1 = np.datetime64(f"{st['end'][:4]}-{st['end'][4:6]}-{st['end'][6:]}")
    n = 1 + len(VORDER)                                   # configurations per instrument/method
    fin = VORDER[-1]                                      # the fully corrected configuration
    i_chm, i_cld, i_ray = 0, n, 2 * n                     # "none" indices; + k for VORDER[k-1]
    vcol = {"none": "#1f77b4", "bg": "#9467bd", "dark": "#17becf", "darkbg": "#2ca02c"}
    vname = {"none": "no correction (operational)", "bg": "+ background", "dark": "+ dark",
             "darkbg": "+ dark + background"}
    out = dict(period=f"{t0} to {t1}", constants={}, ratio={}, chm15k_over_operational={})
    kc_d, kc_v = _kalman(CALIB_DASH / f"{PAY}_C_cloud.csv")
    fig = plt.figure(figsize=(18, 10.5))
    gs = GridSpec(2, 6, figure=fig, height_ratios=[0.85, 1.1], hspace=0.33, wspace=0.6)
    ax_t = fig.add_subplot(gs[0, 0:4])
    ax_r = fig.add_subplot(gs[0, 4:6])
    d, v = _calib_points(CALIB_DASH / f"{PAY}_C_cloud.csv")
    s = (d >= t0) & (d <= t1)
    ks = (kc_d >= t0) & (kc_d <= t1)
    ax_t.plot(d[s], v[s], "s", color="#404040", ms=4, alpha=0.75, mec="none")
    ax_t.plot(kc_d[ks], kc_v[ks], "-", color="#404040", lw=2, label=f"CL61 cloud (O'Connor) ({s.sum()} days)")
    out["constants"]["cloud"] = dict(n=int(s.sum()), median_kalman=float(np.median(kc_v[ks])))
    bins = np.arange(0.5, 1.4, 0.025)
    for vk in ("none",) + VORDER:
        if vk == "none":
            pts, kal = (CALIB_DASH / f"{PAY}_C_rayleigh.csv", None), CALIB_DASH / f"{PAY}_C_rayleigh.csv"
        else:
            pts, kal = (_vrun(vk, "C", "cal"), "rayleigh"), CALIB_DASH / f"{PAY}_C_rayleigh{vk}.csv"
        d, v = _calib_points(*pts)
        s = (d >= t0) & (d <= t1)
        d, v = d[s], v[s]
        kd, kv = _kalman(kal)
        ks = (kd >= t0) & (kd <= t1)
        if vk in ("none", fin):
            ax_t.plot(d, v, "o", color=vcol[vk], ms=4, alpha=0.75, mec="none")
            ax_t.plot(kd[ks], kv[ks], "-", color=vcol[vk], lw=2,
                      label=f"CL61 Rayleigh {vname[vk]} ({s.sum()} nights)")
        cl = np.interp(d.astype("datetime64[D]").astype("f8"), kc_d.astype("datetime64[D]").astype("f8"), kc_v)
        rat = v / cl
        out["constants"][vk] = dict(n_nights=int(s.sum()), median_kalman=float(np.median(kv[ks])))
        out["ratio"][vk] = dict(median=float(np.median(rat)), q25=float(np.percentile(rat, 25)),
                                q75=float(np.percentile(rat, 75)), n=int(rat.size))
        ax_r.hist(rat, bins=bins, histtype="step", lw=1.8, color=vcol[vk],
                  label=f"{vname[vk]}: {np.median(rat):.2f} ({rat.size} n.)")
    ka_d, ka_v = _kalman(CALIB_DASH / f"{PAY}_A_rayleigh.csv")
    for vk in VORDER:
        kd, kv = _kalman(CALIB_DASH / f"{PAY}_A_rayleigh{vk}.csv")
        m = (kd >= t0) & (kd <= t1)
        base = np.interp(kd[m].astype("f8"), ka_d.astype("datetime64[D]").astype("f8"), ka_v)
        out["chm15k_over_operational"][vk] = float(np.median(kv[m] / base))
    ax_t.set_ylabel("CL61 lidar constant C$_L$")
    ax_t.set_title("(a) CL61 calibration constants: nightly values (markers) and Kalman estimate (lines)")
    ax_t.grid(alpha=0.3)
    ax_t.legend(fontsize=8.5, loc="lower left")
    ax_r.axvline(1, color="k", lw=0.9)
    ax_r.set_xlabel("C$_{L,Rayleigh}$ (night) / C$_{L,cloud}$ (Kalman)")
    ax_r.set_ylabel("Number of nights")
    ax_r.set_title("(b) CL61 Rayleigh-to-cloud constant ratio (median)")
    ax_r.legend(fontsize=8, loc="upper left")
    ax_r.grid(alpha=0.3)
    # (c) measured darks in beta units (linear axis) with the molecular signal (median + IQR) for scale
    ax_d = fig.add_subplot(gs[1, 0:2])
    cl61 = out["constants"]["cloud"]["median_kalman"]
    chm = float(np.median(ka_v[(ka_d >= t0) & (ka_d <= t1)]))
    src = np.load(VARIANTS["dark"]["npz_src"])
    for lam, col, nm in ((910.55, "#7f3c8d", "910.55 nm"), (1064.47, "#d62728", "1064 nm")):
        mm, q1, q3 = molecular_iqr(c["time"], z, lam)
        ax_d.fill_betweenx(z, q1, q3, color=col, alpha=0.3, lw=0)
        ax_d.plot(mm, z, ":", color=col, lw=1.4, label=f"molecular {nm} (median, IQR)")
        out[f"molecular_median_{int(lam)}"] = [float(np.interp(h, z, mm)) for h in (3000, 5000, 8000)]
    for ident, cc, col, nm in (("C", cl61, "#7f3c8d", "CL61"), ("A", chm, "#d62728", "CHM15k")):
        dz = np.load(_vnpz("dark", ident))
        bd = dz["b_phys"] / cc * 1e6
        se = np.nan_to_num(np.asarray(src[f"{ident}_sem"], "f8")) / cc * 1e6
        ax_d.fill_betweenx(dz["rng"], bd - se, bd + se, color=col, alpha=0.25, lw=0)
        ax_d.plot(bd, dz["rng"], color=col, lw=1.6, label=f"{nm}: measured dark / C$_L$ (± s.e.)")
    ax_d.axvline(0, color="k", lw=0.8)
    ax_d.set_xlim(*((-0.03, 0.03) if zplot <= 6000 else (-0.06, 0.12)))
    ax_d.set_ylim(0, zplot)
    ax_d.set_xlabel(r"[Mm$^{-1}$ sr$^{-1}$]")
    ax_d.set_ylabel("Altitude a.g.l. [m]")
    ax_d.set_title("(c) Measured dark signals and molecular signal (linear scale)")
    ax_d.grid(alpha=0.3)
    ax_d.legend(fontsize=7.5, loc="upper right")
    # (d, e)
    ax_p = fig.add_subplot(gs[1, 2:4])
    ax_c = fig.add_subplot(gs[1, 4:6])
    minn = max(30, int(0.02 * beta.shape[1]))
    common = np.all(np.isfinite(beta), axis=0)
    show = {i_chm: ("#d62728", "-"), i_chm + n - 1: ("#f4a582", "--"), i_cld: ("#404040", "-"),
            i_cld + n - 1: ("#9a9a9a", "--"), i_ray: ("#1f77b4", "-"), i_ray + n - 1: ("#17becf", "--")}
    xmax = 0.0
    for k, (col, ls) in show.items():
        med, q1, q3 = _median_iqr(np.where(common, beta[k], np.nan), minn)
        ax_p.plot(med, z, ls, color=col, lw=1.7, label=labels[k])
        xmax = max(xmax, np.nanpercentile(med[(z > 150) & (z < zplot)], 99))
    _band(ax_p, ZBAND_STATION)
    if zplot > 6000:
        mm, q1, q3 = molecular_iqr(c["time"], z, 1064.47)
        ax_p.plot(mm, z, ":", color="k", lw=1.6, label=r"molecular $\beta_{att}$ (1064 nm, CAMS)")
        ax_p.set_xscale("log")
        ax_p.set_xlim(3e-3, 1.5)
    else:
        ax_p.set_xlim(0, 1.4 * xmax)
    ax_p.set_ylim(0, zplot)
    ax_p.set_xlabel(BLAB)
    ax_p.set_ylabel("Altitude a.g.l. [m]")
    ax_p.grid(alpha=0.3)
    ax_p.legend(fontsize=7.5, loc="upper right")
    ax_p.set_title(f"(d) Median of matched profiles (no correction / {vname[fin][2:]})")
    layers = [(100, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 5000), (5000, 8000)]
    rd, out["pairs"] = [], {}
    pcol = {i_cld: "#404040", i_cld + n - 1: "#9a9a9a", i_ray: "#1f77b4", i_ray + n - 1: "#17becf"}
    for k, r in DARK_PAIRS:
        s = pair_stats(beta[k], beta[r], z, ZBAND_STATION, layers)
        out["pairs"][f"{labels[k]} vs {labels[r]}"] = s
        if k in pcol:
            short = ("CL61 cloud" if k < i_ray else "CL61 Rayleigh") + \
                (", no correction" if k in (i_cld, i_ray) else f", {vname[fin][2:]}")
            rd.append((rel_diff(beta[k], beta[r]), pcol[k], short, s["med"]))
    panel_reldiff(ax_c, z, rd, ZBAND_STATION, zplot, minn,
                  "(e) Median relative difference to the CHM15k\n(CHM15k with the same correction)")
    ax_c.legend(fontsize=7.8, loc="upper left")
    fig.suptitle(f"Payerne CL61: cloud vs Rayleigh calibration with and without the measured dark signal, "
                 f"{_period(c['time'])}", fontsize=12.5, fontweight="bold", y=0.995)
    fig.savefig(OUT / f"fig_v09_{key}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    out["n_common_hours"] = int(np.any(common[:, (z >= 500) & (z <= 3000)], axis=1).sum())
    return out


def _work(key):
    warnings.filterwarnings("ignore")
    return process_station(key) if key in ALL_STATIONS else process_earlinet(key)


def _figure(key):
    if key in APPENDIX_DARK:
        return fig_appendix_dark(key)
    return fig_station(key) if key in STATIONS else fig_earlinet(key)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    plot_only = "--plot-only" in sys.argv
    keys = args or list(STATIONS) + list(EARLINET) + list(APPENDIX_DARK)
    OUT.mkdir(parents=True, exist_ok=True)
    _configure()
    if any(k in APPENDIX_DARK for k in keys):
        prepare_dark_inputs()
    if not plot_only:
        workers = max(1, min(len(keys), int(os.environ.get("ALC_VAL_SITE_WORKERS", "6"))))
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for k in ex.map(_work, keys):
                print(f"== processed {k}", flush=True)
    summary_f = OUT / "summary_v09.json"
    summary = json.loads(summary_f.read_text()) if summary_f.exists() else {}
    for k in keys:
        summary[k] = _figure(k)
        print(f"-> fig_v09_{k}.png", flush=True)
    summary_f.write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))
    print("PAPER_V09_DONE", flush=True)


if __name__ == "__main__":
    main()
