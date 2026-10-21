# -*- coding: utf-8 -*-
"""Payerne intercomparison — the correction LADDER, one intercomparison figure per configuration.

Reads the three co-located Payerne L1 streams ONCE, CACHES them to a fast .npz (so every later run
skips the slow archive read), then renders each rung of the correction ladder IN PARALLEL (one
worker process per configuration).  Everything downstream of the L1 read is a cheap post-calibration
correction re-applied to the cached data:

  * calibration source   -> IC.CALIB (calib/ = new per-instrument MS,  calib.bak = legacy MS)
  * overlap (CHM15k)     -> OV.correct_rcs on/off
  * electronic dark      -> subtract b_phys(z) (CL31 hood-dark twin)
  * water vapour (910)   -> IC.apply_wv on/off
  * wavelength -> 1064   -> 'molecular' (component-separated) | 'angstrom' (single exponent) | 'none'

Figure layout (per request): a 2x2 grid with NO curtains — the median+/-IQR profile on the LEFT
(full height, x fixed 0..1.5 Mm-1 sr-1), the scatter-vs-reference TOP-RIGHT, and the difference
histogram BOTTOM-RIGHT (x fixed -0.5..0.5).  CL31 is present from rung 0.

Run:  python presentation/make_ladder_figures.py            (uses the npz cache if present)
      python presentation/make_ladder_figures.py --reread   (force a fresh L1 read + rewrite cache)
Output: presentation/figs_intercomp/ladder/NN_<label>.png  + ladder_summary.csv
"""
import os
import sys
import csv
import warnings
from concurrent.futures import ProcessPoolExecutor

warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
# cap per-worker BLAS threads so the process pool does not oversubscribe the CPU
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "3")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from validation.paper import intercompare as IC
from validation.paper import overlap as OV
from validation.paper import run_paper_validation as RP
from validation.paper import figures as FIG
from validation.paper.calib_benchmark import key_of

WMO = "0-20000-0-06610"
STA = (46.8137, 6.9425, 491.0)
START, END = "20260301", "20260630"
TARGET, ALPHA = 1064.0, 1.0
ZMIN, ZMAX = RP.ZMIN, RP.ZMAX
XMAX_PROF = 1.5        # profile x-axis [Mm-1 sr-1], fixed on every figure
HX = 0.5               # difference-histogram x-axis half-width, fixed on every figure

ROOT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python"
CALIB_NEW = Path(ROOT) / "calib"
CALIB_BAK = Path(ROOT) / "calib.bak_20260707_pre_newcloud"
HOOD_NPZ  = str(Path(ROOT) / "cl31_b_dark.npz")

HERE = Path(r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/presentation")
OUT = HERE / "figs_intercomp" / "ladder"
OUT.mkdir(parents=True, exist_ok=True)
L1_NPZ = HERE / "figs_intercomp" / "_l1_cache_payerne.npz"


# --------------------------------------------------------------------- L1 npz cache
_SCALAR = ("station_alt", "lat", "lon", "wavelength", "itype", "wmo", "ident")


def save_l1_cache(cache, path):
    flat = {}
    for ident, l1 in cache.items():
        for k, v in l1.items():
            if v is None:
                continue
            flat[f"{ident}###{k}"] = np.asarray(v)
    np.savez(path, **flat)          # uncompressed = fastest reload


def load_l1_cache(path):
    d = np.load(path, allow_pickle=True)
    cache = {}
    for key in d.files:
        ident, k = key.split("###")
        cache.setdefault(ident, {})[k] = d[key]
    for l1 in cache.values():
        for k, v in list(l1.items()):
            if k in _SCALAR and getattr(v, "ndim", None) == 0:
                l1[k] = v.item()
        l1.setdefault("calc", None)
    return cache


def get_l1(reread=False):
    if L1_NPZ.exists() and not reread:
        cache = load_l1_cache(L1_NPZ)
        if any(k.startswith("L2_") for k in cache):
            print(f"loaded cached L1+L2 from {L1_NPZ.name}", flush=True)
            return cache
        print("cache lacks L2 -> re-reading", flush=True)
    print("reading Payerne L1 + operational L2 (A=CHM15k, B=CL31, C=CL61) ...", flush=True)
    cache = {}
    for ident in ("A", "B", "C"):
        cache[ident] = IC.read_l1(WMO, ident, START, END)
        cache["L2_" + ident] = IC.read_l2(WMO, ident, START, END)
        print(f"   {ident}: L1 {cache[ident]['beta'].shape}  L2 {cache['L2_'+ident]['beta'].shape}",
              flush=True)
    save_l1_cache(cache, L1_NPZ)
    print(f"cached -> {L1_NPZ.name} ({L1_NPZ.stat().st_size/1e6:.0f} MB)", flush=True)
    return cache


# --------------------------------------------------------------------- channel defs
def _ch(ident, itype, calib, label, **extra):
    d = dict(wmo=WMO, ident=ident, itype=itype, calib=calib, label=label,
             lat=STA[0], lon=STA[1], alt=STA[2])
    d.update(extra)
    return d


def cl(ch, calib_dir):
    d = dict(ch); d["calib_dir"] = str(calib_dir); return d


CHM   = _ch("A", "CHM15k", "rayleigh", "CHM15k (Rayleigh)")
CL31C = _ch("B", "CL31",  "cloud",    "CL31 (cloud)")
CL61C = _ch("C", "CL61",  "cloud",    "CL61 (cloud)")
CL61R = _ch("C", "CL61",  "rayleigh", "CL61 (Rayleigh)")
CL31H = _ch("B", "CL31",  "cloud",    "CL31 (cloud, hood-dark corr)", bdark=HOOD_NPZ)

# Operational L2 product channels (rung 0 = the network BEFORE this work): the E-PROFILE L2
# attenuated_backscatter_0 as distributed, with its operational calibration_constant_0
# (CHM15k Rayleigh v1.0; CL31 the 1e8 default = uncalibrated; CL61 Vaisala's internal constant),
# and NO water-vapour / wavelength / overlap harmonisation.
CHM_L2  = _ch("A", "CHM15k", "rayleigh", "CHM15k (L2)", source="l2")
CL31_L2 = _ch("B", "CL31",  "cloud",    "CL31 (L2)",    source="l2")
CL61_L2 = _ch("C", "CL61",  "cloud",    "CL61 (L2)",    source="l2")


# --------------------------------------------------------------------- one channel, one config
def process(l1, ch, cfg, apply_overlap=True):
    beta = l1["beta"].copy()
    if apply_overlap and cfg["overlap"] and ch["itype"] in ("CHM15k", "CHM8k"):
        model = OV.load_overlap_model(ch["wmo"], ch["ident"])
        if model is not None:
            beta = OV.correct_rcs(beta, l1["alt"] - l1["station_alt"], l1["temp_int"], model)
    if cfg["dark"] and ch.get("bdark"):
        bd = np.load(ch["bdark"])
        beta = beta - np.interp(l1["alt"] - l1["station_alt"], bd["rng"], bd["b_phys"],
                                left=0.0, right=0.0)[None, :]
    saved = IC.CALIB
    try:
        IC.CALIB = Path(ch.get("calib_dir", CALIB_NEW))
        cal = IC.load_calib_series(key_of(ch), "L1")
    finally:
        IC.CALIB = saved
    if cal is None:
        return None
    ck = IC.interp_calib(cal[0], cal[1], l1["time"])
    beta = beta / ck[:, None] * 1e6
    med_corr = float(np.nanmedian(1e6 / ck))
    wl = l1["wavelength"] if np.isfinite(l1["wavelength"]) else IC.WV_PARAMS.get(ch["itype"], (np.nan,))[0]
    is_wv = bool(np.isfinite(wl) and IC.in_water_vapor_band(wl))
    if cfg["wv"] and is_wv:
        lam0, fwhm = IC.WV_PARAMS.get(ch["itype"], (910.0, 3.4))
        beta, _ = IC.apply_wv(beta, l1, lam0, fwhm)
    wl_use = wl if np.isfinite(wl) else TARGET
    if ch["itype"] in ("CL31", "CL51", "CL61") and abs(wl_use - TARGET) >= 1.0:
        if cfg["wl"] == "molecular":
            beta = IC.wavelength_correct_molecular(beta, l1, wl_use, TARGET, ALPHA)
        elif cfg["wl"] == "angstrom":
            beta = IC.wavelength_correct(beta, l1, wl_use, TARGET, ALPHA, "angstrom")
    beta_disp = beta.copy()
    beta_disp[l1["qf"] > 0] = np.nan
    beta_scr = IC.screen(beta, l1)
    return dict(beta_scr=beta_scr, beta_disp=beta_disp, med_corr=med_corr)


def process_l2(l2, ch):
    """Rung 0: the operational L2 product AS DISTRIBUTED — no calibration, no WV, no wavelength,
    no overlap. Just the provider attenuated_backscatter_0, screened for cloud/quality."""
    beta = l2["beta"].copy()
    beta_disp = beta.copy()
    beta_disp[l2["qf"] > 0] = np.nan
    beta_scr = IC.screen(beta, l2)
    return dict(beta_scr=beta_scr, beta_disp=beta_disp, med_corr=1.0)


def build_R(l1_cache, channels, cfg):
    IC.set_noise_filter("none")
    items = []
    for ch in channels:
        if ch.get("source") == "l2":
            l1 = l1_cache["L2_" + ch["ident"]]
            r = process_l2(l1, ch); noovl = None
        else:
            l1 = l1_cache[ch["ident"]]
            r = process(l1, ch, cfg)
            noovl = None
            if r is not None and ch["itype"] in ("CHM15k", "CHM8k") \
                    and OV.load_overlap_model(ch["wmo"], ch["ident"]):
                n = process(l1, ch, cfg, apply_overlap=False)
                noovl = n["beta_scr"] if n else None
        if r is None:
            print(f"    [skip] {ch['label']}: no data/series"); return None
        items.append(dict(ch=ch, l1=l1, beta_scr_noovl=noovl, **r))
    return RP.grid_and_stats(items, 0)


# --------------------------------------------------------------------- the 2x2 figure (no curtains)
def fig_ladder(R, out_png, title, zmax_plot=6000):
    z = np.asarray(R["altGrid"]) - R["station"]["altitude"]
    zmask = (z >= 0) & (z <= zmax_plot); zc = z[zmask]
    zmin, zmax = ZMIN, ZMAX
    band = (z >= zmin) & (z <= zmax)
    nch = len(R["channels"]); cols = FIG.channel_colors(R["channels"])

    fig = plt.figure(figsize=(15, 8.2))
    gs = GridSpec(2, 2, width_ratios=[1.12, 1.0], hspace=0.34, wspace=0.22,
                  left=0.065, right=0.975, top=0.90, bottom=0.09)

    # (a) profile — LEFT, full height, fixed x 0..XMAX_PROF
    axp = fig.add_subplot(gs[:, 0])
    have = [np.any(np.isfinite(R["beta"][k][:, zmask]), axis=1) for k in range(nch)]
    common = np.logical_and.reduce(have); ncom = int(common.sum()); synced = ncom > 0
    for k in range(nch):
        B = R["beta"][k][:, zmask]
        if synced:
            B = B[common]
        med, q1, q3, nz = FIG._median_iqr(B)
        keep = FIG._truncate_noise_floor(med, nz, zc, zmin, B.shape[0])
        c = cols[k]
        for q in (q1, q3):
            m = keep & np.isfinite(q) & (q > 0)
            axp.plot(q[m], zc[m], ":", color=c, lw=1.0)
        axp.plot(med[keep], zc[keep], "-", color=c, lw=1.9, label=R["channels"][k]["label"])
        Bn = R.get("beta_noovl", [None] * nch)[k]
        if Bn is not None:
            Bn = Bn[:, zmask]
            if synced:
                Bn = Bn[common]
            medn, _, _, nzn = FIG._median_iqr(Bn)
            keepn = FIG._truncate_noise_floor(medn, nzn, zc, zmin, Bn.shape[0])
            axp.plot(medn[keepn], zc[keepn], "--", color=c, lw=1.3,
                     label=R["channels"][k]["label"] + " (no overlap corr)")
    axp.axhline(zmin, ls="--", color="k", lw=0.8); axp.axhline(zmax, ls="--", color="k", lw=0.8)
    axp.set_xlim(0, XMAX_PROF); axp.set_ylim(0, zmax_plot)
    axp.grid(alpha=0.3)
    axp.set_xlabel(r"$\beta_{att}$ [Mm$^{-1}$ sr$^{-1}$]"); axp.set_ylabel("Altitude a.g.l. [m]")
    axp.set_title("(a) Median (solid) $\\pm$ IQR (dotted)   N=%d common hours" % ncom, fontsize=11)
    axp.legend(loc="upper right", fontsize=9)

    # (b) scatter vs reference — TOP-RIGHT
    axs = fig.add_subplot(gs[0, 1])
    ref = R["beta"][0][:, band].ravel()
    allv = []
    for k in range(1, nch):
        cur = R["beta"][k][:, band].ravel()
        m = np.isfinite(cur) & np.isfinite(ref) & (cur > 0) & (ref > 0)
        a, b = ref[m], cur[m]
        if a.size > 6000:
            sel = np.random.default_rng(k).permutation(a.size)[:6000]; a, b = a[sel], b[sel]
        st = R["stats"][k]
        axs.scatter(a, b, 4, color=cols[k], alpha=0.25, edgecolors="none",
                    label="%s (r=%.2f, %+.0f%%)" % (R["channels"][k]["label"], st["r"], st["relbias_pct"]))
        allv += [a, b]
    if allv:
        allv = np.concatenate(allv)
        lim = (max(np.nanmin(allv), 1e-2), np.nanpercentile(allv, 99.8))
        axs.plot(lim, lim, "k--", lw=1.0); axs.set_xlim(*lim); axs.set_ylim(*lim)
    axs.set_xscale("log"); axs.set_yscale("log"); axs.grid(alpha=0.3)
    axs.set_xlabel("%s [Mm$^{-1}$ sr$^{-1}$]" % R["channels"][0]["label"])
    axs.set_ylabel("channel [Mm$^{-1}$ sr$^{-1}$]")
    axs.set_title("(b) Scatter vs %s (%.0f-%.0f m)" % (R["channels"][0]["label"], zmin, zmax), fontsize=11)
    axs.legend(loc="lower right", fontsize=8)

    # (c) difference histogram — BOTTOM-RIGHT, fixed x -HX..HX
    axh = fig.add_subplot(gs[1, 1])
    refb = R["beta"][0][:, band]
    edges = np.linspace(-HX, HX, 61)
    for k in range(1, nch):
        d = (R["beta"][k][:, band] - refb).ravel(); d = d[np.isfinite(d)]
        axh.hist(d, edges, density=True, histtype="step", color=cols[k], lw=1.7,
                 label="%s (med %+.2f)" % (R["channels"][k]["label"], np.median(d)))
    axh.axvline(0, ls="--", color="k", lw=0.8); axh.set_xlim(-HX, HX)
    axh.grid(alpha=0.3); axh.set_xlabel(r"$\beta_{att}$ difference [Mm$^{-1}$ sr$^{-1}$]")
    axh.set_ylabel("pdf")
    axh.set_title("(c) Difference vs %s" % R["channels"][0]["label"], fontsize=11)
    axh.legend(loc="upper right", fontsize=8)

    fig.suptitle(title, fontweight="bold", fontsize=12)
    fig.savefig(out_png, dpi=160); plt.close(fig)


# --------------------------------------------------------------------- ladder (distinct rungs only)
def _cfg(**kw):
    c = dict(overlap=True, dark=False, wv=False, wl="none"); c.update(kw); return c


LADDER = [
    ("00_before_L2", "operational L2 product, no corrections (network status BEFORE this work)",
     [CHM_L2, CL31_L2, CL61_L2], _cfg()),

    ("01_add_cloud_legacyMS", "+ CL61 liquid-cloud calibration (legacy 8 um MS table)",
     [CHM, cl(CL31C, CALIB_NEW), cl(CL61C, CALIB_BAK), cl(CL61R, CALIB_NEW)], _cfg(wl="molecular", wv=True)),

    ("02_per_instrument_MS", "+ per-instrument multiple-scattering (PVC 5.5 um)",
     [CHM, cl(CL31C, CALIB_NEW), cl(CL61C, CALIB_NEW), cl(CL61R, CALIB_NEW)], _cfg(wl="molecular", wv=True)),

    ("03_no_water_vapour", "water-vapour correction OFF (counterfactual)",
     [CHM, cl(CL31C, CALIB_NEW), cl(CL61C, CALIB_NEW), cl(CL61R, CALIB_NEW)], _cfg(wl="molecular", wv=False)),

    ("05_wavelength_angstrom", "wavelength -> 1064 by a single Angstrom exponent (naive)",
     [CHM, cl(CL31C, CALIB_NEW), cl(CL61C, CALIB_NEW), cl(CL61R, CALIB_NEW)], _cfg(wl="angstrom", wv=True)),

    ("07_hood_dark", "molecular-aware wavelength + hood-dark on CL31 (fully corrected)",
     [CHM, cl(CL31C, CALIB_NEW), cl(CL31H, CALIB_NEW), cl(CL61C, CALIB_NEW), cl(CL61R, CALIB_NEW)],
     _cfg(wl="molecular", wv=True, dark=True)),
]

_WORKER_L1 = None


def _init_worker(npz_path):
    global _WORKER_L1
    matplotlib.use("Agg")
    _WORKER_L1 = load_l1_cache(npz_path)


def _render_one(spec):
    tag, title, channels, cfg = spec["tag"], spec["title"], spec["channels"], spec["cfg"]
    R = build_R(_WORKER_L1, channels, cfg)
    if R is None:
        return tag, []
    fig_ladder(R, OUT / f"{tag}.png",
               f"Payerne (0-20000-0-06610)  —  {title}  —  2026-03-01 to 2026-06-30")
    rows = []
    for k, ch in enumerate(channels):
        st = R["stats"][k]
        rows.append(dict(config=tag, channel=ch["label"], med_relbias_pct=round(st["medrelbias_pct"], 2),
                         log_r=round(st["r_log"], 3), relbias_pct=round(st["relbias_pct"], 2),
                         r=round(st["r"], 3), n=st["n"]))
    return tag, rows


def main():
    reread = "--reread" in sys.argv[1:]
    get_l1(reread=reread)                                   # ensure the npz cache exists
    specs = [dict(tag=t, title=ti, channels=ch, cfg=c) for (t, ti, ch, c) in LADDER]
    workers = min(len(specs), max(1, (os.cpu_count() or 8) - 2), 6)
    print(f"rendering {len(specs)} configs over {workers} worker processes ...", flush=True)
    all_rows = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(str(L1_NPZ),)) as ex:
        for tag, rows in ex.map(_render_one, specs):
            print(f"   -> {tag}.png", flush=True)
            for r in rows:
                print("        %-30s med %+6.1f%%  logr %.2f  N=%d"
                      % (r["channel"], r["med_relbias_pct"], r["log_r"], r["n"]), flush=True)
            all_rows += rows
    with open(OUT / "ladder_summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["config", "channel", "med_relbias_pct", "log_r",
                                          "relbias_pct", "r", "n"])
        w.writeheader(); w.writerows(all_rows)
    print("LADDER_DONE ->", OUT, flush=True)


if __name__ == "__main__":
    main()
