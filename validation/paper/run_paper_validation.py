"""
run_paper_validation.py — uniform, all-Level-1 attenuated-backscatter validation for the paper.

ONE methodology is applied identically to every instrument at every site (no L2 shortcuts, no MATLAB
reference). For each channel the calibrated attenuated backscatter is built from the native Level-1
range-corrected signal in these steps:

  1. Read L1 rcs_0 from the daily archive (intercompare.read_l1), which also carries the internal
     temperature temp_int [degC] and the cloud/fog fields used for screening.
  2. Overlap correction (CHM15k only): rcs_0 *= 1 + (a(z)*T + b(z))/100, the temperature-dependent
     Hervo-2016 model (overlap.correct_rcs). A no-op above ~700 m (full overlap).
  2b. Electronic-offset correction (flagged units only, shown as an additional "offset-corr" channel
     next to the native one): rcs_0 -= b_phys(z), the fixed digitizer-ripple pattern extracted from
     clear nights (doc/reports/network_offset_coefficients.md). The cloud calibration constant is
     offset-immune, so the same C_L applies to both variants.
  3. Calibration: beta_att [Mm^-1 sr^-1] = rcs_0 / C_L * 1e6, with C_L the daily Kalman lidar
     constant (calib_benchmark.py, L1-derived) interpolated to the profile times. The SAME formula
     serves the Rayleigh and the liquid-cloud constants — both are true lidar constants on the rcs_0
     scale, so no L2 provider constant and no CL61 unit exception are involved.
  4. Water-vapour correction (910 nm: CL31/CL51/CL61): divide beta by the two-way WV transmission
     from monthly CAMS at the instrument laser line (intercompare.apply_wv). Its per-instrument
     impact is quantified by re-running the channel with the correction OFF.
  5. Wavelength normalisation to the site reference wavelength (Angstrom alpha=1; Mini-MPL 532->1064
     via the molaer model) — intercompare.wavelength_correct.
  6. Screening (quality flag, any cloud base 0-20 km, fog/vertical visibility, +/-15 min expansion)
     into a science stream; a display stream keeps clouds visible — intercompare.screen.
  7. 60-min median grid with the >=30-min coverage and per-gate SNR>=3 detection requirements
     (intercompare.retime_hourly), then a common altitude grid (build_common_grid / regrid).
  8. Statistics vs the site reference channel over 500-3000 m AGL (intercompare._stats): the robust
     pair median-relative-bias / log-r is the headline; linear relbias / r are kept for continuity.

The EARLINET (ceilometer vs research-lidar) comparison is run through the same L1 methodology by
earlinet.compare(). Outputs: per-site figures, the calibration time series, a water-vapour-impact
figure, the EARLINET figures, and summary_stats.csv. No markdown report, no .mat files.

Usage:  python -m validation.paper.run_paper_validation [site ...]     (default: all sites)
"""
from __future__ import annotations
import csv
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path

# Cap the per-worker BLAS/thread count so the site ProcessPool (below) does not oversubscribe the CPU
# (each worker still reads L1 with its own thread pool). Set before numpy imports its BLAS backend.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from validation.paper import intercompare as IC
from validation.paper import overlap as OV
from validation.paper import figures as FIG
from validation.paper import earlinet as EA
from validation.paper.calib_benchmark import BENCHMARK, key_of

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
CALIB = IC.CALIB                                   # where calib_benchmark.py wrote <key>_L1.csv

# Comparison band and aerosol Angstrom exponent (shared by every site).
ZMIN, ZMAX, ALPHA = 500.0, 3000.0, 1.0

# Per site: reference channel index (into BENCHMARK[name]['channels']) and the reference wavelength
# all channels are normalised to. 1064 nm (CHM15k reference) everywhere except Uccle (CL51, 910 nm).
SITE = {
    "payerne":    dict(ref=0, target=1064.0, wmo="0-20000-0-06610"),
    "amsterdam":  dict(ref=0, target=1064.0, wmo="0-20000-0-06240"),
    "uccle":      dict(ref=0, target=910.0,  wmo="0-20000-0-06447"),
    "sirta":      dict(ref=0, target=1064.0, wmo="0-250-1001-07151"),
    "lindenberg": dict(ref=0, target=1064.0, wmo="0-20000-0-10393"),
    "aosta":      dict(ref=0, target=1064.0, wmo="0-380-5-1"),
    "camborne":   dict(ref=0, target=1064.0, wmo="0-20000-0-03808"),
}
SITE_NAME = {"payerne": "Payerne", "amsterdam": "Amsterdam", "uccle": "Uccle", "sirta": "Palaiseau",
             "lindenberg": "Lindenberg", "aosta": "Aosta", "camborne": "Camborne"}

# Electronic-offset corrections (doc/reports/network_offset_coefficients.md): units flagged
# correctable get an ADDITIONAL "offset-corr" channel — the same instrument with the clear-night
# digitizer-ripple pattern b_phys(z) subtracted from rcs_0 — so native and corrected are both in the
# results. The Uccle CL51 carries the network's strongest offset (40 m ripple, 9.4 % of the
# mid-range signal, split-half repro 1.00). Same cloud C_L (the O'Connor method is offset-immune).
NETWORK_OFFSET = OUT / "network_offset"
EXTRA_CHANNELS = {
    "uccle": [dict(wmo="0-20000-0-06447", ident="A", itype="CL51", calib="cloud",
                   label="CL51 (cloud, offset-corr)",
                   bdark=str(NETWORK_OFFSET / "0-20000-0-06447_A_CL51.npz"))],
    # Payerne CL31 offset-corr: subtract the MEASURED terminal-hood dark, not the clear-night ripple.
    # cl31_b_dark.npz stores b_phys = P_dark(z)*z^2 with P_dark = median(rcs_0/z^2) over the 4 pooled
    # hood sessions (incl. the ~25 h 26-27 May), fitted to a two-resonance physical model (_cl31_offset_
    # physical_model.py). Subtracting b_phys from rcs_0 BEFORE calibration == removing the dark in the
    # raw P = beta/r^2 space, then re-forming beta and calibrating (the working CHM15k recipe).
    "payerne": [dict(wmo="0-20000-0-06610", ident="B", itype="CL31", calib="cloud",
                     label="CL31 (cloud, hood-dark corr)",
                     bdark=str(OUT / "cl31_b_dark.npz"))],
}
EARLINET_LABEL = {"sir": "Palaiseau", "ino": "Magurele", "ino_a": "Magurele", "ari": "Leipzig",
                  "lei": "Leipzig", "cbw": "Cabauw", "sir_532": "Palaiseau 532 nm"}


# --------------------------------------------------------------------------- per-channel pipeline
def channel_beta(l1, ch, target, apply_wv=True, apply_overlap=True):
    """Steps 2-6 for one channel on an already-read L1 dict. Returns a dict with the screened
    (science) and display beta streams, the median calibration multiplier, and whether the channel
    carries a water-vapour correction. Returns None if the channel has no calibration series.
    apply_wv / apply_overlap switch steps 4 / 2 off for the impact-quantification twins."""
    beta = l1["beta"].copy()                                   # L1 rcs_0

    # (2) overlap correction — CHM15k only, using the internal temperature carried by read_l1.
    if apply_overlap and ch["itype"] in ("CHM15k", "CHM8k"):
        model = OV.load_overlap_model(ch["wmo"], ch["ident"])
        if model is not None:
            beta = OV.correct_rcs(beta, l1["alt"] - l1["station_alt"], l1["temp_int"], model)
        else:
            print(f"    [overlap] {ch['label']}: no model for {ch['wmo']}/{ch['ident']} -> skipped")

    # (2b) electronic-offset correction — "offset-corr" channels only: subtract the fixed
    # digitizer-ripple pattern b_phys(z) (rcs_0 units, network_offset_coefficients.md) from rcs_0.
    if ch.get("bdark"):
        bd = np.load(ch["bdark"])
        beta = beta - np.interp(l1["alt"] - l1["station_alt"], bd["rng"], bd["b_phys"],
                                left=0.0, right=0.0)[None, :]

    # (3) calibration — beta_att = rcs_0 / C_L * 1e6 with the L1-derived Kalman constant.
    cal = IC.load_calib_series(key_of(ch), "L1")
    if cal is None:
        return None
    ck = IC.interp_calib(cal[0], cal[1], l1["time"])
    beta = beta / ck[:, None] * 1e6
    med_corr = float(np.nanmedian(1e6 / ck))

    # (4) water vapour — 910 nm family only.
    wl = l1["wavelength"] if np.isfinite(l1["wavelength"]) else IC.WV_PARAMS.get(ch["itype"], (np.nan,))[0]
    is_wv = bool(np.isfinite(wl) and IC.in_water_vapor_band(wl))
    if apply_wv and is_wv:
        lam0, fwhm = IC.WV_PARAMS.get(ch["itype"], (910.0, 3.4))
        beta, info = IC.apply_wv(beta, l1, lam0, fwhm)
        if info["months_excluded"]:
            print(f"    [wv] {ch['label']}: no usable CAMS for {','.join(info['months_excluded'])}"
                  " -> month(s) NaN-masked")

    # (5) wavelength normalisation to the site reference. The 910 nm Vaisalas (CL31/CL51/CL61) use the
    # OPTIMAL component-separated conversion: the analytic molecular (Rayleigh) part is computed from
    # CAMS T/p (0.4 deg, 1 deg fallback; hydrostatic below the lowest CAMS level) and the Angstrom law is
    # applied to the aerosol residual only — physically correct where the molecular fraction is large.
    # The 532 nm Mini-MPL keeps its molaer (US-std) path; a same-wavelength pair is a no-op either way.
    wl_use = wl if np.isfinite(wl) else target
    if ch["itype"] in ("CL31", "CL51", "CL61") and abs(wl_use - target) >= 1.0:
        beta = IC.wavelength_correct_molecular(beta, l1, wl_use, target, ALPHA)
    else:
        model = "molaer" if ch["itype"] in ("Mini-MPL", "MPL") else "angstrom"
        beta = IC.wavelength_correct(beta, l1, wl_use, target, ALPHA, model)

    # (6) science (screened) + display (clouds visible) streams.
    beta_disp = beta.copy()
    beta_disp[l1["qf"] > 0] = np.nan
    beta_scr = IC.screen(beta, l1)
    return dict(beta_scr=beta_scr, beta_disp=beta_disp, med_corr=med_corr, is_wv=is_wv)


# --------------------------------------------------------------------------- gridding + stats
def _reindex(A, pos, n):
    """Place rows of A onto the union time grid (pos<0 -> NaN row)."""
    out = np.full(n if A.ndim == 1 else (n, A.shape[1]), np.nan)
    ok = pos >= 0
    out[ok] = A[pos[ok]]
    return out


def grid_and_stats(items, iref):
    """Grid a list of processed channels (steps 7-8) and return the R dict fig_multi_alc expects.
    Each item carries its L1 dict (for time/alt/cbh) and the screened + display beta streams.
    Mirrors intercompare.process but L1-only."""
    gridded = []
    for it in items:
        l1 = it["l1"]
        cbh = l1["cbh"]
        cbh_low = (np.nanmin(np.where((cbh > 0) & (cbh < 20000), cbh, np.nan), axis=1)
                   if cbh.ndim == 2 else cbh)
        arrays = [it["beta_scr"], it["beta_disp"], cbh_low]
        snr_idx = (0,)
        if it.get("beta_scr_noovl") is not None:      # overlap-uncorrected twin, same gates
            arrays.append(it["beta_scr_noovl"])
            snr_idx = (0, 3)
        g, arrs = IC.retime_hourly(l1["time"], arrays, min_cov_s=IC.MIN_AVG_S, snr_idx=snr_idx)
        gridded.append(dict(it=it, l1=l1, grid=g, scr=arrs[0], disp=arrs[1], cbh=arrs[2],
                            scrn=(arrs[3] if len(arrs) > 3 else None)))

    union = np.unique(np.concatenate([g["grid"] for g in gridded]))
    for g in gridded:
        idx = {t: i for i, t in enumerate(g["grid"])}
        pos = np.array([idx.get(t, -1) for t in union])
        g["scrU"] = _reindex(g["scr"], pos, union.size)
        g["dispU"] = _reindex(g["disp"], pos, union.size)
        g["cbhU"] = _reindex(g["cbh"], pos, union.size)
        g["scrnU"] = _reindex(g["scrn"], pos, union.size) if g["scrn"] is not None else None

    altGrid = IC.build_common_grid([g["l1"]["alt"] for g in gridded])
    for g in gridded:
        g["betaC"] = IC.regrid(g["scrU"], g["l1"]["alt"], altGrid)
        g["dispC"] = IC.regrid(g["dispU"], g["l1"]["alt"], altGrid)
        g["noovlC"] = IC.regrid(g["scrnU"], g["l1"]["alt"], altGrid) if g["scrnU"] is not None else None

    station_alt = gridded[0]["l1"]["station_alt"]
    zmask = (altGrid >= ZMIN + station_alt) & (altGrid <= ZMAX + station_alt)
    ref = gridded[iref]["betaC"]
    R = dict(altGrid=altGrid, time_sync=union, station=dict(
        altitude=station_alt, lat=gridded[0]["l1"]["lat"], lon=gridded[0]["l1"]["lon"]),
        channels=[], beta=[], beta_disp=[], beta_noovl=[], cbh=[], stats=[])
    for g in gridded:
        ch = g["it"]["ch"]
        R["channels"].append(dict(label=ch["label"], calib=ch["calib"], itype=ch["itype"],
                                  wavelength=g["l1"]["wavelength"], med_corr=g["it"]["med_corr"]))
        R["beta"].append(g["betaC"])
        R["beta_disp"].append(g["dispC"])
        R["beta_noovl"].append(g["noovlC"])
        R["cbh"].append(g["cbhU"])
        R["stats"].append(IC._stats(g["betaC"], ref, zmask))
    return R


# --------------------------------------------------------------------------- one site
def run_site(name):
    """Read L1 once per channel, run the pipeline (WV on, and WV off for 910 nm channels), grid and
    compute the statistics. Returns (R, cfg, rows) or None. R + cfg feed fig_multi_alc; rows feed the
    summary. The water-vapour impact is the in-band change in beta between the WV-on and WV-off runs
    (median beta_on/beta_off ~ median 1/T^2), plus the change in median relative bias vs the reference."""
    st = BENCHMARK[name]
    sc = SITE[name]
    channels = list(st["channels"]) + EXTRA_CHANNELS.get(name, [])
    items_on, items_off = [], []
    l1_cache = {}                      # (wmo, ident) -> L1 dict; offset-corr twins and the CL61's
    for ch in channels:                # two methods share one physical instrument -> read L1 once
        key = (ch["wmo"], ch["ident"])
        if key not in l1_cache:
            l1_cache[key] = IC.read_l1(ch["wmo"], ch["ident"], st["start"], st["end"])
        l1 = l1_cache[key]
        if l1 is None:
            print(f"    [skip] {ch['label']}: no L1 in window")
            continue
        on = channel_beta(l1, ch, sc["target"], apply_wv=True)
        if on is None:
            print(f"    [skip] {ch['label']}: no calibration series ({key_of(ch)}_L1)")
            continue
        # overlap-uncorrected twin (CHM15k with a model only): same pipeline, step 2 off — shown
        # as the dashed median in the station figure so the correction's effect is visible.
        noovl = None
        if ch["itype"] in ("CHM15k", "CHM8k") and OV.load_overlap_model(ch["wmo"], ch["ident"]):
            noovl = channel_beta(l1, ch, sc["target"], apply_wv=True, apply_overlap=False)
        items_on.append(dict(ch=ch, l1=l1, **on,
                             beta_scr_noovl=(noovl["beta_scr"] if noovl else None)))
        # WV-off twin for the impact quantification (non-WV channels reuse the same result so the
        # channel list, the reference and the grids stay aligned between the two passes).
        off = channel_beta(l1, ch, sc["target"], apply_wv=False) if on["is_wv"] else on
        items_off.append(dict(ch=ch, l1=l1, **off))
    if not items_on:
        return None

    cfg = dict(referenceChannel=sc["ref"], zMin=ZMIN, zMax=ZMAX, start=st["start"], end=st["end"],
               channels=[it["ch"] for it in items_on])
    R = grid_and_stats(items_on, sc["ref"])
    R_off = grid_and_stats(items_off, sc["ref"])

    band = ((R["altGrid"] - R["station"]["altitude"]) >= ZMIN) & \
           ((R["altGrid"] - R["station"]["altitude"]) <= ZMAX)
    rows = []
    for k, it in enumerate(items_on):
        s = R["stats"][k]
        wv_pct = wv_dmedrel = np.nan
        if it["is_wv"]:
            on = R["beta"][k][:, band]
            off = R_off["beta"][k][:, band]
            m = np.isfinite(on) & np.isfinite(off) & (off > 0)
            if m.any():
                wv_pct = 100.0 * (float(np.nanmedian(on[m] / off[m])) - 1.0)   # in-band beta increase
            wv_dmedrel = s["medrelbias_pct"] - R_off["stats"][k]["medrelbias_pct"]
        rows.append(dict(kind="station", site=name, label=it["ch"]["label"], itype=it["ch"]["itype"],
                         calib=it["ch"]["calib"], ref=(k == sc["ref"]),
                         medrel=s["medrelbias_pct"], rlog=s["r_log"], relbias=s["relbias_pct"],
                         r=s["r"], n=s["n"], is_wv=it["is_wv"], wv_pct=wv_pct, wv_dmedrel=wv_dmedrel))
    return R, cfg, rows


# --------------------------------------------------------------------------- EARLINET
def run_earlinet(all_rows):
    """Ceilometer vs EARLINET research lidar (CHM/Mini-MPL side now read from L1 via earlinet.compare).
    Appends one summary row per site to all_rows and writes the 2x2 figures."""
    for code in ("sir", "ino", "ino_a", "ari", "lei", "cbw", "sir_532"):
        try:
            s = EA.compare(code, "20250101", "20260630", return_profiles=True)
        except Exception as exc:
            s = {"error": repr(exc)}
        if not s or "error" in s:
            print(f"  EARLINET {code}: {s}")
            continue
        site = EA.SITES[code]
        FIG.fig_earlinet(code, EARLINET_LABEL[code], s["betaE"], s["betaC"], s["grid"], s["times"],
                         s, {}, OUT / f"fig_earlinet_{code}.png", betaC_raw=s.get("betaC_raw"),
                         instr=site.get("instr", "CHM15k (Rayleigh)"), itype=site.get("itype", "CHM15k"))
        print("  EARLINET %-8s med %+6.1f%%  logr %.2f  relbias %+6.1f%%  r %.2f  matched=%d -> fig_earlinet_%s.png"
              % (code, s.get("medrelbias_pct", np.nan), s.get("r_log", np.nan),
                 s["relbias_pct"], s.get("r", np.nan), s["matched"], code), flush=True)
        all_rows.append(dict(kind="earlinet", site=code, label=EARLINET_LABEL[code],
                             itype=site.get("itype", "CHM15k"), calib="rayleigh", ref=False,
                             medrel=s.get("medrelbias_pct", np.nan), rlog=s.get("r_log", np.nan),
                             relbias=s["relbias_pct"], r=s["r"], n=s["matched"],
                             is_wv=False, wv_pct=np.nan, wv_dmedrel=np.nan))


# --------------------------------------------------------------------------- figures / summary
def calib_channel_list():
    """One entry per instrument for the combined C_L time-series grid (both CL61 methods share a
    panel). Reads the calib CSVs written by calib_benchmark.py."""
    order, seen = [], {}
    for name, st in BENCHMARK.items():
        for c in st["channels"]:
            stream = f"{c['wmo']}_{c['ident']}"
            k = key_of(c)
            if stream in seen:
                if k not in [s["key"] for s in seen[stream]["series"]]:
                    seen[stream]["series"].append(dict(key=k, calib=c["calib"], itype=c["itype"]))
                continue
            site = SITE_NAME.get(name) or c["label"].split(" ")[0]
            base = c["label"].split(" (")[0]
            title = f"{site} {base}" if SITE_NAME.get(name) else base
            entry = dict(title=title, unit="C$_L$", series=[dict(key=k, calib=c["calib"], itype=c["itype"])])
            seen[stream] = entry
            order.append(entry)
    return order


def make_wv_impact_figure(wv_rows, out_png):
    """Landscape horizontal-bar figure: per 910 nm channel, the in-band % by which the water-vapour
    correction raises beta_att (= median 1/T^2 - 1), annotated with the change in median relative
    bias vs the reference."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = [f"{SITE_NAME.get(r['site'], r['site'])} — {r['label']}" for r in wv_rows]
    vals = [r["wv_pct"] for r in wv_rows]
    colors = [FIG.TYPE_COLORS.get(r["itype"], "#1f77b4") for r in wv_rows]
    y = np.arange(len(wv_rows))
    fig, ax = plt.subplots(figsize=(12, 0.6 * len(wv_rows) + 2))
    ax.barh(y, vals, color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel(r"Water-vapour impact on $\beta_{att}$ in 500–3000 m  [%]  (median $1/T^2_{wv}-1$)")
    ax.set_title("Impact of the 910 nm water-vapour correction, per instrument", fontsize=12, fontweight="bold")
    for yi, r in zip(y, wv_rows):
        note = "" if not np.isfinite(r["wv_dmedrel"]) else f"  Δmed rel bias {r['wv_dmedrel']:+.1f}%"
        ax.text(r["wv_pct"] + (0.4 if r["wv_pct"] >= 0 else -0.4), yi,
                f"+{r['wv_pct']:.1f}%{note}", va="center",
                ha="left" if r["wv_pct"] >= 0 else "right", fontsize=8)
    vmax = max([v for v in vals if np.isfinite(v)] + [1.0])
    ax.set_xlim(right=vmax * 1.45)          # headroom so the annotations stay inside the axes
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"   -> {out_png.name}", flush=True)


def write_summary(rows, out_csv):
    cols = ["kind", "site", "label", "itype", "calib", "ref", "medrel", "rlog", "relbias", "r", "n",
            "is_wv", "wv_pct", "wv_dmedrel"]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    print(f"   -> {out_csv.name} ({len(rows)} rows)", flush=True)


def _d(yyyymmdd):
    return datetime.strptime(yyyymmdd, "%Y%m%d").strftime("%Y-%m-%d")


# --------------------------------------------------------------------------- main
def _process_site(name):
    """ProcessPool worker: process one site end-to-end (read L1, grid, stats, write its station figure)
    and return (name, rows). The big R arrays stay in the worker; only the small summary rows return."""
    warnings.filterwarnings("ignore")
    res = run_site(name)
    if res is None:
        return name, None
    R, cfg, rows = res
    title = "%s (%s)  —  %s to %s" % (SITE_NAME[name], SITE[name]["wmo"],
                                      _d(cfg["start"]), _d(cfg["end"]))
    FIG.fig_multi_alc(R, cfg, OUT / f"fig_{name}.png", title)
    return name, rows


def main():
    warnings.filterwarnings("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    # Args: site names run just those stations; the pseudo-site "earlinet" runs just the EARLINET
    # comparison; no args = full run (all stations + EARLINET + calibration-series figure).
    req = sys.argv[1:]
    full = not req
    sites = list(SITE) if full else [a for a in req if a in SITE]
    do_earlinet = full or ("earlinet" in req)
    all_rows, wv_rows = [], []

    # Sites are independent -> process them in parallel. ALC_VAL_SITE_WORKERS caps the pool (default 6).
    workers = max(1, min(len(sites), int(os.environ.get("ALC_VAL_SITE_WORKERS", "6"))))
    site_rows = {}
    if workers > 1 and len(sites) > 1:
        print(f"[parallel] {len(sites)} sites over {workers} workers", flush=True)
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for name, rows in ex.map(_process_site, sites):
                site_rows[name] = rows
    else:
        for name in sites:
            n, rows = _process_site(name)
            site_rows[n] = rows
    for name in sites:
        rows = site_rows.get(name)
        if not rows:
            print(f"== {name} ==   no data", flush=True)
            continue
        print(f"== {name} ==", flush=True)
        for r in rows:
            wv = (" | WV +%.1f%% (Δmed rel %+.1f%%)" % (r["wv_pct"], r["wv_dmedrel"])) if r["is_wv"] else ""
            tag = " (ref)" if r["ref"] else ""
            print("   %-28s med %+6.1f%%  logr %.2f  relbias %+6.1f%%  r %.2f  N=%7d%s%s"
                  % (r["label"] + tag, r["medrel"], r["rlog"], r["relbias"], r["r"], r["n"], wv, ""),
                  flush=True)
            all_rows.append(r)
            if r["is_wv"]:
                wv_rows.append(r)
        print(f"   -> fig_{name}.png", flush=True)

    if wv_rows:
        make_wv_impact_figure(wv_rows, OUT / "fig_wv_impact.png")

    # Global products: only on a full run — a targeted `... payerne` or `... earlinet` stays fast.
    if full:
        FIG.fig_calib_timeseries(calib_channel_list(), CALIB, OUT / "fig_calib_timeseries.png")
        print("   -> fig_calib_timeseries.png", flush=True)
    if do_earlinet:
        print("== EARLINET ==", flush=True)
        run_earlinet(all_rows)

    suffix = "" if full else "_" + "_".join(req)
    write_summary(all_rows, OUT / f"summary_stats{suffix}.csv")
    print("PAPER_VALIDATION_DONE", flush=True)


if __name__ == "__main__":
    main()
